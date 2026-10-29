"""Per-replicate count records of the M7-13b novel-ORF gate (docs/specs/M7-13-calibration-gate.md
rules 1, 5-8, 10, § Attacks a, b, d, GC; M7-11b rule 11; lead decision 1 ``GATE_R``). Test
plumbing owned by the test-author: each record simulates ``novel_orf_power`` and runs GIPPA through
the ``novel_orfs`` step, then reduces the written ``tables/novel_orfs.tsv`` and the side table
``novel_class_fdr`` with the pure helpers of ``tests/calibration/_m7_gate.py`` (implementer), so
no check reruns GIPPA and the harness logs the counts (shards merge from them).

The chain (``run_gippa``) is ``gippa run --executor local`` with the local runner's targets cut to
``novel_orfs`` (every layer is precomputed by gippasim; no tool runs; plots and report are not
built). A test ``$GIPPA_DB_DIR`` is used (never the real one). The side table is read from
``work/novel_class_fdr`` (else the published ``novel/novel_class_fdr.parquet``); it exists only
once M7-11b (#256) lands, so every record raises ``FileNotFoundError`` before that.

Seeds (rule 1): replicate ``i`` simulates with ``int(replicate_rng(GATE_BASE, i).integers(0,
2**31 - 1))``; every arm of replicate i (database, de novo, decoys stripped, G4-k3, GC, RBS
mutant, clustered null) re-simulates that same seed, so the arms are paired. Test-side draws use
``replicate_rng(K3_BASE | GC_BASE | CL_BASE, i)``.

Test-side transforms (no simulator change):

* ``inject_false_orfs`` (G4-k3, rule 8): the M7-12 non-coding six-frame pool
  (``gippasim.proteomics.sixframe_noncoding_pool`` rebuilt from the Bakta-like JSON and the novel
  truth); a *segment* is a stop-free stretch of one frame closed by a stop codon holding >= 5
  distinct pool products. Per class of ``SIX_CLASSES`` 20 segments are drawn; each gets 3-5 of
  its products as wrong de novo reads, half of them (``k // 2``) carrying one reversed internal
  stretch of 3-6 residues (the §8.5 attack 5 error type; spec "half of them as reversed
  stretches of the segment"). ``nterm_extension`` gets none: its pool products are upstream of a
  gene, and a false N-terminal candidate is an extension of that gene, not a segment.
* ``gc_poor_peptides`` (attack GC): 40 intergenic segments drawn from the GC-poorest quartile of
  segments with >= 2 products, 1-3 exact products each.
* ``inject_null_entries`` (the clustered-null arm of attacks a and d, M7-13b lead decision 7;
  statistician probe, Slurm 1874863): ``CL_ENTRIES`` = 400 null ``antisense`` makedb entries
  appended to the database ``peptides.tsv``. An entry is a segment with >= 5 products whose
  makedb ID (``makedb_id``) the simulation did not hit (no target or ``rev_`` row) and that
  shares no key with, and does not overlap, any planted ORF (either strand). Each entry goes to
  the target or the decoy side with p = 1/2 and gets k = 3-5 distinct products as PSM rows
  (decoy: the product reversed, protein ``rev_<id>``), each with a ``(q_value, pep, score)``
  triple resampled from the simulation's ``rev_GIPPAORF_`` rows (the null component); target
  intensity resampled from the ``GIPPAORF_`` target rows, ``run1``-``run3``, charge 2-3.

Appended reads copy the de novo TSV columns: fresh ``spectrum_id`` ``inj.<tag>.<n>`` in ``run1``,
charge 2, score and intensity resampled from the file's own rows. Every record is appended as one
JSON line to ``<log_dir()>/novel_orfs_gate_replicates.jsonl``.
"""

from __future__ import annotations

import dataclasses
import json
import math
import re
import shutil
import tempfile
import time
from collections.abc import Mapping, Sequence
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

GATE_BASE, K3_BASE, GC_BASE = "novel_orfs.gate", "novel_orfs.g4k3", "novel_orfs.gc"
CL_BASE = "novel_orfs.cluster"
RECORDS_FILE = "novel_orfs_gate_replicates.jsonl"
REGIMES = ("picked_decoy", "denovo_subgroup")
#: the four gated classes (restated so test modules collect before ``_m7_gate`` has them)
NOVEL_CLASSES = ("intergenic", "cds_outframe", "antisense", "nterm_extension")
SIX_CLASSES = ("intergenic", "cds_outframe", "antisense")
K3_PER_CLASS, K3_PEPTIDES, K3_MIN_PRODUCTS = 20, (3, 5), 5
GC_SEGMENTS, GC_PEPTIDES = 40, (1, 3)
CL_CLASS, CL_ENTRIES, CL_PEPTIDES, CL_MIN_PRODUCTS, CL_DECOY_P = "antisense", 400, (3, 5), 5, 0.5
RBS_MUTANT_QUANTILE = 0.50
FDR_KEYS = ("n_called", "n_false", "n_targets")


def _C() -> Any:
    import _criteria

    return _criteria


def _G() -> Any:
    import _m7_gate

    return _m7_gate


def sim_seed(i: int) -> int:
    return int(_C().replicate_rng(GATE_BASE, i).integers(0, 2**31 - 1))


def levels() -> tuple[float, ...]:
    return tuple(_C().load_tolerances().levels)


def lvl(q: float) -> str:
    return f"{q:g}"


def g4_key(regime: str, cls: str, q: float) -> str:
    return f"{regime}.{cls}.{lvl(q)}"


def _log(record: Mapping[str, Any]) -> None:
    d = _C().log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / RECORDS_FILE).open("a") as fh:
        fh.write(json.dumps(record, default=float) + "\n")


def _tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


# --------------------------------------------------------------------------- chain


@dataclasses.dataclass
class NovelRun:
    out: Path
    novel_orfs: pd.DataFrame  # tables/novel_orfs.tsv as written (strings)
    side: pd.DataFrame  # novel_class_fdr (M7-11b rule 11)

    @property
    def work(self) -> Path:
        return self.out / "work"


def simulate(i: int, size: str, root: Path) -> Any:
    from gippasim import scenarios as sc
    from gippasim.simulate import simulate as run_sim

    n = sc.NOVEL_PER_CLASS_GATE if size == "1M" else None
    return run_sim(sc.novel_orf_power(sim_seed(i), size, n_per_class=n), sim_seed(i), root / "sim")


def run_gippa(res: Any, peptides: Path, out: Path, seed: int, sets: Sequence[str] = ()) -> NovelRun:
    """``gippa run`` (local executor) through ``novel_orfs`` only; RuntimeError on a non-zero exit
    (never an AssertionError, so a broken chain cannot pass as a strict xfail)."""
    from click.testing import CliRunner

    from gippa.cli import cli
    from gippa.runner import local
    from helpers.novel_steps_fixture import db_env

    orig = local.run_local

    def only_novel(cfg: Any, wd: Path, targets: Any, *a: Any, **k: Any) -> None:
        orig(cfg, wd, ["novel_orfs"], *a, **k)

    args = ["run", "--genome", str(res.genome_fasta), "--peptides", str(peptides),
            "--annotation", str(res.bakta_json), "--config", str(res.config_yaml),
            "--outdir", str(out), "--modules", "core", "--seed", str(seed),
            "--executor", "local"]  # fmt: skip
    for s in sets:
        args += ["--set", s]
    with db_env(out.parent), pytest.MonkeyPatch.context() as mp:
        mp.setattr(local, "run_local", only_novel)
        r = CliRunner().invoke(cli, args)
    if r.exception is not None and not isinstance(r.exception, SystemExit):
        raise RuntimeError(f"gippa run raised {r.exception!r}") from r.exception
    if r.exit_code != 0:
        raise RuntimeError(f"gippa run exit {r.exit_code}: {r.output[-3000:]}")
    side_path = out / "work" / "novel_class_fdr"
    if not side_path.exists():
        side_path = out / "novel" / "novel_class_fdr.parquet"
    if not side_path.exists():
        raise FileNotFoundError(f"no novel_class_fdr side table under {out} (M7-11b, #256)")
    return NovelRun(out, _tsv(out / "tables" / "novel_orfs.tsv"), pd.read_parquet(side_path))


def _lengths(res: Any) -> dict[str, int]:
    return {k: len(v) for k, v in _genome(res).items()}


def _genome(res: Any) -> dict[str, str]:
    from gippa.io.fasta import read_genome_fasta

    return {r.id: r.sequence.upper() for r in read_genome_fasta(Path(res.genome_fasta))}


def _truth(res: Any) -> pd.DataFrame:
    from gippasim.truth import read_novel_orfs

    return read_novel_orfs(res.outdir)


# --------------------------------------------------------------------------- reductions


def _pool_counts(by_class: Mapping[str, Mapping[str, int]]) -> dict[str, int]:
    """Sum over the four classes, plus ``n_truth_<class>`` (G3 guard keys, rule 7)."""
    G = _G()
    out = {k: sum(int(by_class[c][k]) for c in G.NOVEL_CLASSES) for k in G.NOVEL_COUNT_KEYS}
    out.update({f"n_truth_{c}": int(by_class[c]["n_truth"]) for c in G.NOVEL_CLASSES})
    return out


def _fdr_counts(m: Mapping[str, int], n_targets: int) -> dict[str, int]:
    return {"n_called": int(m["n_pred"]), "n_false": int(m["n_pred"] - m["tp_p"]),
            "n_targets": int(n_targets)}  # fmt: skip


def g4_counts(run: NovelRun, truth: pd.DataFrame, regime: str, lengths: Mapping[str, int],
              side: pd.DataFrame | None = None) -> dict[str, dict[str, int]]:  # fmt: skip
    """``{g4_key(regime, cls, q): {n_called, n_false, n_targets}}`` at every level."""
    G = _G()
    s = run.side if side is None else side
    targets = G.list_targets(s, regime)
    out: dict[str, dict[str, int]] = {}
    for q in levels():
        m = G.match_novel(G.g4_predictions(run.novel_orfs, s, q, regime=regime, lengths=lengths),
                          truth, lengths)  # fmt: skip
        for c in G.NOVEL_CLASSES:
            out[g4_key(regime, c, q)] = _fdr_counts(m[c], targets[c])
    return out


def g3_counts(run: NovelRun, truth: pd.DataFrame, lengths: Mapping[str, int]) -> dict[str, Any]:
    G = _G()
    return dict(G.match_novel(G.g3_predictions(run.novel_orfs, lengths), truth, lengths))


def false_candidates(run: NovelRun, truth: pd.DataFrame, lengths: Mapping[str, int]) -> dict:
    """Rows per class (called or not) at no truth key of their class (recorded, rule 10)."""
    G = _G()
    every = run.novel_orfs.assign(called="true")
    m = G.match_novel(G.g3_predictions(every, lengths), truth, lengths)
    return {c: int(m[c]["n_pred"] - m[c]["tp_p"]) for c in G.NOVEL_CLASSES}


def _diag(run: NovelRun, rows: pd.DataFrame, truth: pd.DataFrame, lengths: Mapping[str, int],
          genes: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:  # fmt: skip
    G = _G()
    return {"reasons": G.reason_counts(run.side), "fallback_calls": G.fallback_calls(run.side),
            SHARED_NO_PAIR: shared_no_pair_counts(run.novel_orfs, run.side, rows, truth,
                                                  lengths),
            NONATG_NTERM: nonatg_nterm_counts(run.novel_orfs, run.side, rows,
                                              genes)}  # fmt: skip


# --------------------------------------------------------------------------- shared peptides

#: lead decision 11: ``no_makedb_pair`` ORFs whose every peptide row names another planted ORF's
#: makedb entry (chance recurrence of a short peptide at another ATG; recorded, not a bug)
SHARED_NO_PAIR = "shared_peptide_no_pair"
NO_PAIR = "no_makedb_pair"
_MAKEDB_RE = re.compile(
    r"^GIPPAORF_(intergenic|cds_outframe|antisense|nterm_extension)_(.+)_(\d+)_(\d+)_([pm])$"
)


def _il(p: str) -> str:
    """I/L/J folded (M7-13 rule 3 / example 3: mass-equivalent residues compare equal)."""
    return str(p).strip().upper().replace("I", "L").replace("J", "L")


#: injection reach (review B2 on #271): the injected keys that reached GIPPA's table, summed over
#: the replicates run, must be at least this share of the injected ones. Fast tier (Slurm 1877255
#: / 1878529, R = 10-20): G4-k3 600 of 600 (1.0 in every replicate), GC 399 of 400 (minimum
#: 39 of 40 = 0.975 in one replicate); 0.90 leaves room below the lowest replicate seen.
REACH_MIN_SHARE = 0.90


def check_reach(rows: int, injected: int, what: str) -> None:
    """RuntimeError (never an AssertionError, so a strict xfail cannot absorb it) when the
    injected false ORFs did not reach GIPPA: ``rows`` is 0 or below ``REACH_MIN_SHARE`` of
    ``injected`` (``injected`` < 1 is itself a broken transform)."""
    if injected < 1 or rows < 1 or rows < REACH_MIN_SHARE * injected:
        raise RuntimeError(
            f"{what}: injection did not reach GIPPA ({rows} rows at {injected} injected keys,"
            f" need >= {REACH_MIN_SHARE:.2f} of them)"
        )


def parse_makedb_id(pid: str) -> tuple[str, str, int, int, str] | None:
    """``(cls, contig, a, b, strand)`` of ``GIPPAORF_<cls>_<contig>_<a>_<b>_<p|m>`` (1-based closed
    plus-strand span without the stop codon), else None (gene loci, ``rev_`` decoys)."""
    m = _MAKEDB_RE.match(str(pid).strip())
    if m is None:
        return None
    cls, contig, a, b, s = m.groups()
    return cls, contig, int(a), int(b), "+" if s == "p" else "-"


def planted_in_entry(entry: tuple[str, str, int, int, str],
                     truth: pd.DataFrame) -> set[tuple[str, str, int]]:  # fmt: skip
    """``(contig, strand, stop)`` of the truth ORFs whose start codon (``start`` on ``+``,
    ``end - 2`` on ``-``) lies in frame inside the entry on its strand (origin-crossing entries
    (``a > b``) are not unwrapped: no ORF)."""
    _, contig, a, b, strand = entry
    out: set[tuple[str, str, int]] = set()
    if a > b:  # origin-crossing entry: not unwrapped, so it never exempts (fails closed)
        return out
    for c, st, s, e, stop in zip(truth["contig"], truth["strand"], truth["start"], truth["end"],
                                 truth["stop"], strict=True):  # fmt: skip
        if str(c) != contig or str(st) != strand:
            continue
        p = int(s) if strand == "+" else int(e) - 2  # lowest base of the start codon
        if a <= p and p + 2 <= b and ((p - a) if strand == "+" else (b - p - 2)) % 3 == 0:
            out.add((contig, strand, int(stop)))
    return out


def shared_peptide_no_pair(peptides: Sequence[str], rows: pd.DataFrame,
                           own_key: tuple[str, str, int], truth: pd.DataFrame) -> bool:  # fmt: skip
    """Lead decision 11: at least one input row carries one of ``peptides`` (I/L-insensitive),
    and every such row names at least one ``GIPPAORF_`` entry, each holding a planted ORF other
    than ``own_key``. Gene loci, decoys and non-planted entries veto."""
    want = {_il(p) for p in peptides if str(p).strip()}
    if not want:
        return False
    sub = rows[rows["peptide"].map(_il).isin(want)]
    if sub.empty:
        return False
    own = (str(own_key[0]), str(own_key[1]), int(own_key[2]))
    for ids in sub["protein_ids"]:
        names = [x.strip() for x in re.split(r"[;,]", str(ids)) if x.strip()]
        if not names:
            return False
        for x in names:
            e = parse_makedb_id(x)
            if e is None or not (planted_in_entry(e, truth) - {own}):
                return False
    return True


def shared_no_pair_counts(novel_orfs: pd.DataFrame, side: pd.DataFrame, rows: pd.DataFrame,
                          truth: pd.DataFrame,
                          lengths: Mapping[str, int]) -> dict[str, int]:  # fmt: skip
    """Side rows (four classes, by side ``class``) with ``reason == no_makedb_pair`` that pass
    ``shared_peptide_no_pair`` on their ``novel_orfs.tsv`` row."""
    G = _G()
    out = dict.fromkeys(NOVEL_CLASSES, 0)
    by_id = {str(r["orf_id"]): r for r in novel_orfs.to_dict("records")}
    for oid, cls, why in zip(side["orf_id"], side["class"], side["reason"], strict=True):
        if cls not in out or why != NO_PAIR or str(oid) not in by_id:
            continue
        o = by_id[str(oid)]
        peps = [p for p in re.split(r"[;,]", str(o.get("peptides", ""))) if p.strip()]
        c, st = str(o["contig"]), str(o["strand"])
        key = (c, st, int(G.novel_stop(st, int(o["end"]), lengths.get(c))))
        out[str(cls)] += int(shared_peptide_no_pair(peps, rows, key, truth))
    return out


#: lead decision 11 amendment (M11-16f, statistician-approved on #481): ``no_makedb_pair`` ORFs
#: whose every peptide row names a non-ATG annotated gene by locus tag, the peptide being that
#: gene's initiator-M N-terminal tryptic peptide, at a locus not overlapping the gene (six-frame
#: reads the GTG / TTG start as V / L, so the row misses its gene and lands on a chance ATG copy;
#: recorded separately from ``SHARED_NO_PAIR``, not a bug)
NONATG_NTERM = "nonatg_nterm_recurrence"
_CODON_RE = re.compile(r"^[ACGT]{3}$")


def gene_table(bakta_json: Path | str) -> dict[str, dict[str, Any]]:
    """The Bakta-like JSON's features by ``locus`` (features without one are dropped)."""
    data = json.loads(Path(bakta_json).read_text(encoding="utf-8"))
    return {str(f["locus"]): dict(f) for f in data.get("features", []) if f.get("locus")}


def nterm_tryptic(peptide: str, aa: str) -> bool:
    """``peptide`` is an N-terminal Trypsin/P product of ``aa`` (I/L-folded): a non-empty prefix
    that ends after K/R or is the whole protein (missed cleavages allowed, as gippasim digests)."""
    p, a = _il(peptide), _il(aa)
    return bool(p) and a.startswith(p) and (len(p) == len(a) or p[-1] in "KR")


def _gene_overlaps(g: Mapping[str, Any], contig: str, start: int, end: int) -> bool:
    """1-based closed spans on one contig, either strand; origin-crossing spans (``start >
    stop``) and missing coordinates count as overlapping (fails closed)."""
    try:
        gs, ge = int(g["start"]), int(g["stop"])
    except (KeyError, TypeError, ValueError):
        return True
    if str(g.get("sequence", "")) != str(contig):
        return False
    if gs > ge or start > end:
        return True
    return gs <= end and start <= ge


def nonatg_nterm_recurrence(peptides: Sequence[str], rows: pd.DataFrame,
                            orf: tuple[str, int, int],
                            genes: Mapping[str, Mapping[str, Any]]) -> bool:  # fmt: skip
    """The amendment's narrow match: at least one input row carries one of ``peptides``
    (I/L-folded), and every name of every such row is a locus tag in ``genes`` whose
    ``start_type`` is a codon other than ATG, whose ``aa`` begins with M and has the row's
    peptide as its N-terminal tryptic peptide, and which does not overlap ``orf`` = ``(contig,
    lo, hi)`` (1-based closed, ``orf_span``). Anything else vetoes."""
    want = {_il(p) for p in peptides if str(p).strip()}
    if not want:
        return False
    sub = rows[rows["peptide"].map(_il).isin(want)]
    if sub.empty:
        return False
    contig, start, end = str(orf[0]), int(orf[1]), int(orf[2])
    for pep, ids in zip(sub["peptide"], sub["protein_ids"], strict=True):
        names = [x.strip() for x in re.split(r"[;,]", str(ids)) if x.strip()]
        if not names:
            return False
        for x in names:
            g = genes.get(x)
            if g is None:
                return False
            st = str(g.get("start_type") or "").strip().upper()
            aa = str(g.get("aa") or "").strip().upper()
            if not _CODON_RE.match(st) or st == "ATG" or not aa.startswith("M"):
                return False
            if not nterm_tryptic(str(pep), aa) or _gene_overlaps(g, contig, start, end):
                return False
    return True


def orf_span(contig: str, strand: str, start: object, end: object) -> tuple[str, int, int]:
    """``(contig, lo, hi)`` 1-based closed of a ``novel_orfs.tsv`` row from ``start_longest``
    (the widest span, so the overlap test is conservative) and ``end`` (on ``-`` the start is
    the upper base); an origin-crossing ORF gives ``lo > hi`` (counts as overlapping)."""
    s, e = int(str(start)), int(str(end))
    return (contig, s, e) if strand == "+" else (contig, e, s)


def nonatg_nterm_counts(novel_orfs: pd.DataFrame, side: pd.DataFrame, rows: pd.DataFrame,
                        genes: Mapping[str, Mapping[str, Any]]) -> dict[str, int]:  # fmt: skip
    """Side rows (four classes, by side ``class``) with ``reason == no_makedb_pair`` that pass
    ``nonatg_nterm_recurrence`` on their ``novel_orfs.tsv`` row (span: ``orf_span``)."""
    out = dict.fromkeys(NOVEL_CLASSES, 0)
    by_id = {str(r["orf_id"]): r for r in novel_orfs.to_dict("records")}
    for oid, cls, why in zip(side["orf_id"], side["class"], side["reason"], strict=True):
        if cls not in out or why != NO_PAIR or str(oid) not in by_id:
            continue
        o = by_id[str(oid)]
        peps = [p for p in re.split(r"[;,]", str(o.get("peptides", ""))) if p.strip()]
        orf = orf_span(str(o["contig"]), str(o["strand"]), o["start_longest"], o["end"])
        out[str(cls)] += int(nonatg_nterm_recurrence(peps, rows, orf, genes))
    return out


def unestimated_bugs(reasons: Mapping[str, int], shared: int,
                     bug_reasons: frozenset[str]) -> dict[str, int]:  # fmt: skip
    """The ``bug_reasons`` counts left once ``no_makedb_pair`` loses the ``shared`` ORFs."""
    n = int(reasons.get(NO_PAIR, 0))
    if not 0 <= int(shared) <= n:
        raise ValueError(f"{SHARED_NO_PAIR} {shared} outside 0..{NO_PAIR} {n}")
    out = {k: int(v) for k, v in reasons.items() if k in bug_reasons and int(v) > 0}
    if NO_PAIR in out:
        out[NO_PAIR] -= int(shared)
        if out[NO_PAIR] == 0:
            del out[NO_PAIR]
    return out


# --------------------------------------------------------------------------- attack d (pooled)


def pooled_side(run: NovelRun, res: Any) -> tuple[pd.DataFrame, int]:
    """Attack d's mutant: the side table with picked ``q`` from ONE list over the three six-frame
    classes (``intergenic``, ``cds_outframe``, ``antisense``; ``nterm_extension`` keeps its
    locus-tag list, it cannot pair by stop). The ORFs are rebuilt from the run's own work files
    (M7-06 ``candidate_rows`` + ``reconstruct``); ``(side, n side picked rows not rebuilt)``."""
    from gippa.config import resolve_config
    from gippa.io.fasta import read_genome_fasta
    from gippa.io.work_tables import _read_tsv_table
    from gippa.novel import class_call, fdr
    from gippa.novel.candidates import candidate_rows, reconstruct
    from gippa.steps.mapping import read_gene_model

    w = run.work
    cfg = resolve_config(Path(res.config_yaml))
    feats = w / "features"
    gm = read_gene_model(w / "genes", w / "replicons", feats if feats.exists() else None)
    genome = {r.id: r.sequence.upper() for r in read_genome_fasta(w / "genome")}
    peptides = pd.read_parquet(w / "peptides")
    hits = _read_tsv_table("peptide_hits", w / "mapping" / "peptide_hits.tsv")
    orfs = reconstruct(candidate_rows(hits, peptides, None), gm, genome, cfg)
    prefixes = tuple(cfg.input.decoy_prefixes)
    idm = class_call.makedb_id_map(peptides, gm, genome, prefixes)
    idm.loc[idm["class"].isin(SIX_CLASSES), "class"] = "intergenic"
    pooled = [dataclasses.replace(o, cls="intergenic") if o.cls in SIX_CLASSES else o
              for o in orfs]  # fmt: skip
    lengths = {k: r.length for k, r in gm.replicons.items()}
    res_p = fdr.picked_class_fdr(pooled, peptides, idm, decoy_prefixes=prefixes,
                                 replicon_lengths=lengths)  # fmt: skip
    side = run.side.copy()
    picked = side["regime"] == "picked_decoy"
    side.loc[picked, "q"] = [res_p.q.get(o, math.nan) for o in side.loc[picked, "orf_id"]]
    missing = int(sum(o not in res_p.q for o in side.loc[picked, "orf_id"]))
    return side, missing


# --------------------------------------------------------------------------- injection


def _segments(res: Any, classes: Sequence[str]) -> dict[str, list[dict[str, Any]]]:
    """Stop-closed segments of the non-coding six-frame pool with their distinct products."""
    from gippasim.proteomics import ProteomicsConfig, sixframe_noncoding_pool

    js = json.loads(Path(res.bakta_json).read_text())
    cds = [f for f in js["features"] if f.get("type") in ("cds", "sorf")]
    annotated = [(str(f["sequence"]), int(f["start"]) - 1, int(f["stop"]), str(f["strand"]))
                 for f in cds]  # fmt: skip
    tr = _truth(res)
    coding = annotated + [
        (str(c), int(s) - 1, int(e), str(st))
        for c, s, e, st in zip(tr["contig"], tr["start"], tr["end"], tr["strand"], strict=True)
    ]
    starts = [(c, s0 if st == "+" else e0, st) for c, s0, e0, st in annotated]
    pc = ProteomicsConfig()
    pool = sixframe_noncoding_pool(_genome(res), coding, annotated, pc.digest_len,
                                   pc.missed_cleavages, starts)  # fmt: skip
    out: dict[str, list[dict[str, Any]]] = {}
    for cls in classes:
        segs: dict[tuple[int, int], dict[str, Any]] = {}
        for f, a0, b0 in pool.idx[cls].tolist():
            cid, strand, off, aa, n = pool.frames[f]
            lo = aa.rfind(b"*", 0, a0) + 1
            hi = aa.find(b"*", b0)
            if hi < 0:
                continue
            seg = segs.setdefault((f, lo), {"products": set(), "frame": f, "lo": lo, "hi": hi})
            seg["products"].add(aa[a0:b0].decode())
        rows = []
        for (f, lo), seg in sorted(segs.items()):
            cid, strand, off, aa, n = pool.frames[f]
            t = off + 3 * seg["hi"]  # strand-local first base of the stop codon
            stop0 = t if strand == "+" else n - (t + 3)
            rows.append(
                {
                    "contig": cid,
                    "strand": strand,
                    "stop": stop0 + 1,
                    "products": sorted(seg["products"]),
                    "nt": (off + 3 * lo, t, n),
                }
            )
        out[cls] = rows
    return out


def _gc(seg: Mapping[str, Any], genome: Mapping[str, str]) -> float:
    nt0, nt1, n = seg["nt"]
    s0, e0 = (nt0, nt1) if seg["strand"] == "+" else (n - nt1, n - nt0)
    s = genome[seg["contig"]][s0:e0]
    return (s.count("G") + s.count("C")) / max(len(s), 1)


def _reverse_stretch(p: str, rng: np.random.Generator) -> str:
    if len(p) < 5:
        return p
    ln = int(rng.integers(3, min(6, len(p) - 2) + 1))
    a = int(rng.integers(1, len(p) - ln))
    return p[:a] + p[a : a + ln][::-1] + p[a + ln :]


def _append(path: Path, peptides: Sequence[tuple[str, str]], rng: np.random.Generator,
            dst: Path) -> Path:  # fmt: skip
    """``path`` plus one de novo read per ``(tag, peptide)``; written to ``dst``."""
    dn = _tsv(path)
    pick = rng.integers(0, len(dn), size=len(peptides))
    add = pd.DataFrame({c: ["NA"] * len(peptides) for c in dn.columns})
    add["peptide"] = [p for _, p in peptides]
    add["score"] = dn["score"].to_numpy()[pick]
    add["intensity"] = dn["intensity"].to_numpy()[pick]
    add["is_decoy"], add["run"], add["charge"], add["source"] = "false", "run1", "2", "denovo"
    add["spectrum_id"] = [f"inj.{t}.{k}" for k, (t, _) in enumerate(peptides)]
    pd.concat([dn, add], ignore_index=True).to_csv(dst, sep="\t", index=False)
    return dst


def inject_false_orfs(res: Any, rng: np.random.Generator) -> tuple[Path, dict[str, Any]]:
    """The G4-k3 transform (module doc); returns the de novo file and the injected keys."""
    segs = _segments(res, SIX_CLASSES)
    peps: list[tuple[str, str]] = []
    keys: dict[str, list[list[Any]]] = {}
    for cls in SIX_CLASSES:
        ok = [s for s in segs[cls] if len(s["products"]) >= K3_MIN_PRODUCTS]
        take = rng.choice(len(ok), size=min(K3_PER_CLASS, len(ok)), replace=False) if ok else []
        keys[cls] = []
        for j in sorted(int(x) for x in take):
            s = ok[j]
            k = int(rng.integers(K3_PEPTIDES[0], K3_PEPTIDES[1] + 1))
            chosen = [s["products"][int(x)] for x in rng.choice(len(s["products"]), size=k,
                                                                  replace=False)]  # fmt: skip
            for m, p in enumerate(chosen):
                peps.append((cls, _reverse_stretch(p, rng) if m < k // 2 else p))
            keys[cls].append([s["contig"], s["strand"], s["stop"]])
    src = Path(res.outdir) / "peptides_denovo.tsv"
    dst = _append(src, peps, rng, src.with_name("peptides_denovo.k3.tsv"))
    return dst, {"keys": keys, "n_peptides": len(peps)}


def gc_poor_peptides(res: Any, rng: np.random.Generator) -> tuple[Path, dict[str, Any]]:
    """The GC attack's transform (module doc)."""
    genome = _genome(res)
    segs = [s for s in _segments(res, ("intergenic",))["intergenic"] if len(s["products"]) >= 2]
    gc = np.array([_gc(s, genome) for s in segs])
    poor = [s for s, g in zip(segs, gc, strict=True) if g <= np.quantile(gc, 0.25)] if segs else []
    take = rng.choice(len(poor), size=min(GC_SEGMENTS, len(poor)), replace=False) if poor else []
    peps: list[tuple[str, str]] = []
    keys = []
    for j in sorted(int(x) for x in take):
        s = poor[j]
        k = min(int(rng.integers(GC_PEPTIDES[0], GC_PEPTIDES[1] + 1)), len(s["products"]))
        for x in rng.choice(len(s["products"]), size=k, replace=False):
            peps.append(("gc", s["products"][int(x)]))
        keys.append([s["contig"], s["strand"], s["stop"]])
    src = Path(res.outdir) / "peptides_denovo.tsv"
    dst = _append(src, peps, rng, src.with_name("peptides_denovo.gc.tsv"))
    return dst, {"keys": {"intergenic": keys}, "n_peptides": len(peps),
                 "mean_gc": float(gc.mean()) if len(gc) else None,
                 "poor_gc_max": float(np.quantile(gc, 0.25)) if len(gc) else None}  # fmt: skip


def makedb_id(cls: str, seg: Mapping[str, Any]) -> str:
    """gippasim's makedb ID of a segment (stop codon excluded, 1-based closed, plus-strand
    coordinates): ``GIPPAORF_<cls>_<contig>_<a>_<b>_<p|m>`` with ``a = nt0 + 1, b = t`` on ``+``
    and ``a = n - t + 1, b = n - nt0`` on ``-`` (``seg["nt"] = (nt0, t, n)`` strand-local)."""
    nt0, t, n = seg["nt"]
    if seg["strand"] == "+":
        return f"GIPPAORF_{cls}_{seg['contig']}_{nt0 + 1}_{t}_p"
    return f"GIPPAORF_{cls}_{seg['contig']}_{n - t + 1}_{n - nt0}_m"


def segment_span(seg: Mapping[str, Any]) -> tuple[int, int]:
    """0-based half-open plus-strand span of a segment with its stop codon."""
    nt0, t, n = seg["nt"]
    return (nt0, t + 3) if seg["strand"] == "+" else (n - t - 3, n - nt0)


def _planted(res: Any) -> tuple[set[tuple[str, str, int]], dict[str, list[tuple[int, int]]]]:
    tr = _truth(res)
    keys = {(str(c), str(s), int(p)) for c, s, p in zip(tr["contig"], tr["strand"], tr["stop"],
                                                       strict=True)}  # fmt: skip
    spans: dict[str, list[tuple[int, int]]] = {}
    for c, a, b in zip(tr["contig"], tr["start"], tr["end"], strict=True):
        lo, hi = sorted((int(a), int(b)))
        spans.setdefault(str(c), []).append((lo - 1, hi))
    return keys, spans


def null_segments(res: Any, cls: str = CL_CLASS) -> list[dict[str, Any]]:
    """Segments of ``cls`` eligible as clustered-null entries (module doc), with ``id``."""
    pep = _tsv(Path(res.outdir) / "peptides.tsv")
    seen = {p.removeprefix("rev_") for p in pep["protein_ids"] if "GIPPAORF_" in p}
    keys, spans = _planted(res)
    out = []
    for s in _segments(res, (cls,))[cls]:
        oid = makedb_id(cls, s)
        if len(s["products"]) < CL_MIN_PRODUCTS or oid in seen:
            continue
        if (s["contig"], s["strand"], s["stop"]) in keys:
            continue
        lo, hi = segment_span(s)
        if any(a < hi and lo < b for a, b in spans.get(s["contig"], ())):
            continue
        out.append({**s, "id": oid})
    return out


def inject_null_entries(res: Any, rng: np.random.Generator,
                        m: int = CL_ENTRIES) -> tuple[Path, dict[str, Any]]:  # fmt: skip
    """The clustered-null transform (module doc); the database file and the injected entries."""
    pep = _tsv(Path(res.outdir) / "peptides.tsv")
    ids = pep["protein_ids"]
    null = pep[ids.str.startswith("rev_GIPPAORF_")]
    tgt_int = pep.loc[ids.str.startswith("GIPPAORF_"), "intensity"].to_numpy()
    if null.empty or not len(tgt_int):  # not an AssertionError: a broken transform
        raise RuntimeError("clustered null: no GIPPAORF_ target or rev_GIPPAORF_ rows to resample")
    segs = null_segments(res)
    take = rng.choice(len(segs), size=min(m, len(segs)), replace=False)
    rows: list[dict[str, Any]] = []
    entries: list[dict[str, Any]] = []
    for j in sorted(int(x) for x in take):
        s = segs[j]
        decoy = bool(rng.random() < CL_DECOY_P)
        k = int(rng.integers(CL_PEPTIDES[0], CL_PEPTIDES[1] + 1))
        pick = rng.choice(len(s["products"]), size=k, replace=False)
        prods = [s["products"][int(x)] for x in pick]
        for p in prods:
            r = null.iloc[int(rng.integers(0, len(null)))]
            rows.append({"peptide": p[::-1] if decoy else p,
                         "protein_ids": ("rev_" if decoy else "") + s["id"],
                         "q_value": r["q_value"], "pep": r["pep"], "score": r["score"],
                         "is_decoy": "true" if decoy else "false",
                         "intensity": "NA" if decoy else tgt_int[int(rng.integers(0,
                                                                                  len(tgt_int)))],
                         "run": f"run{int(rng.integers(1, 4))}",
                         "charge": str(int(rng.integers(2, 4))),
                         "spectrum_id": f"inj.cl.{len(rows)}", "source": "database",
                         "search_db": "same"})  # fmt: skip
        entries.append({"id": s["id"], "decoy": decoy, "key": [s["contig"], s["strand"], s["stop"]],
                        "span": list(segment_span(s)), "contig": s["contig"],
                        "products": prods})  # fmt: skip
    add = pd.DataFrame(rows).reindex(columns=pep.columns, fill_value="NA")
    dst = Path(res.outdir) / "peptides.cluster.tsv"
    pd.concat([pep, add], ignore_index=True).to_csv(dst, sep="\t", index=False)
    targets = [e["key"] for e in entries if not e["decoy"]]
    return dst, {"keys": {CL_CLASS: targets}, "entries": entries, "n_eligible": len(segs),
                 "n_target_entries": len(targets), "n_rows": len(rows)}  # fmt: skip


def entered(run: NovelRun, keys: Mapping[str, Sequence[Sequence[Any]]],
            lengths: Mapping[str, int]) -> dict[str, int]:  # fmt: skip
    """Rows of the run at an injected key (any class): total, with a finite q, called."""
    G = _G()
    want = {(str(c), str(s), int(p)) for v in keys.values() for c, s, p in v}
    t = run.novel_orfs
    at = [(c, s, G.novel_stop(s, int(e), lengths.get(c))) in want
          for c, s, e in zip(t["contig"], t["strand"], t["end"], strict=True)]  # fmt: skip
    ids = set(t.loc[at, "orf_id"])
    s = run.side[run.side["orf_id"].isin(ids)]
    return {"rows": len(ids), "finite_q": int(np.isfinite(s["q"].astype(float)).sum()),
            "called": int(s["called"].astype(bool).sum())}  # fmt: skip


# --------------------------------------------------------------------------- records


def _root(tag: str) -> Path:
    return Path(tempfile.mkdtemp(prefix=f"m7gate-novel-{tag}-"))


@cache
def gate_record(i: int, size: str = "1M") -> dict[str, Any]:
    """G3, G4 (both regimes), attacks a and d, the recorded arms of replicate ``i``."""
    G = _G()
    root, t0 = _root("gate"), time.perf_counter()
    try:
        res = simulate(i, size, root)
        seed, lengths, truth = sim_seed(i), _lengths(res), _truth(res)
        db = run_gippa(res, Path(res.outdir) / "peptides.tsv", root / "db" / "OUT", seed)
        dn = run_gippa(res, Path(res.outdir) / "peptides_denovo.tsv", root / "dn" / "OUT", seed)
        pep = _tsv(Path(res.outdir) / "peptides.tsv")
        dn_rows = _tsv(Path(res.outdir) / "peptides_denovo.tsv")
        stripped = root / "peptides.nodecoy.tsv"
        pep[pep["is_decoy"].str.lower() != "true"].to_csv(stripped, sep="\t", index=False)
        nd = run_gippa(res, stripped, root / "nd" / "OUT", seed)
        genes = gene_table(res.bakta_json)
        g3 = g3_counts(db, truth, lengths)
        g4 = {
            **g4_counts(db, truth, "picked_decoy", lengths),
            **g4_counts(dn, truth, "denovo_subgroup", lengths),
        }
        a = G.match_novel(G.g4_predictions(db.novel_orfs, db.side, None, regime="picked_decoy",
                                           lengths=lengths), truth, lengths)  # fmt: skip
        targets = G.list_targets(db.side, "picked_decoy")
        pside, missing = pooled_side(db, res)
        d = g4_counts(db, truth, "picked_decoy", lengths, side=pside)
        report = json.loads((Path(res.outdir) / "truth" / "novel_orfs.json").read_text())
        rec = {
            "stream": GATE_BASE, "replicate": i, "seed": seed, "size": size,
            "g3": _pool_counts(g3), "g3_by_class": g3,
            "g3_denovo": _pool_counts(g3_counts(dn, truth, lengths)),
            "g3_unestimated": _pool_counts(g3_counts(nd, truth, lengths)),
            "g4": g4,
            "attack_a": {c: _fdr_counts(a[c], targets[c]) for c in G.NOVEL_CLASSES},
            "attack_d": d, "attack_d_unrebuilt": missing,
            "diag": {"picked_decoy": _diag(db, pep, truth, lengths, genes),
                     "denovo_subgroup": _diag(dn, dn_rows, truth, lengths, genes)},
            "realised": report["realised"],
            "false_candidates": {"picked_decoy": false_candidates(db, truth, lengths),
                                 "denovo_subgroup": false_candidates(dn, truth, lengths)},
            "seconds": time.perf_counter() - t0,
        }  # fmt: skip
    finally:
        shutil.rmtree(root, ignore_errors=True)
    _log(rec)
    return rec


@cache
def k3_record(i: int, size: str = "1M") -> dict[str, Any]:
    """G4-k3: the de novo arm of replicate ``i`` with ``inject_false_orfs`` (rule 8)."""
    root, t0 = _root("k3"), time.perf_counter()
    try:
        res = simulate(i, size, root)
        seed, lengths, truth = sim_seed(i), _lengths(res), _truth(res)
        path, info = inject_false_orfs(res, _C().replicate_rng(K3_BASE, i))
        dn = run_gippa(res, path, root / "k3" / "OUT", seed)
        rec = {"stream": K3_BASE, "replicate": i, "seed": seed, "size": size,
               "g4k3": g4_counts(dn, truth, "denovo_subgroup", lengths),
               "injected": {c: len(v) for c, v in info["keys"].items()},
               "n_peptides": info["n_peptides"], "entered": entered(dn, info["keys"], lengths),
               "diag": _diag(dn, _tsv(path), truth, lengths, gene_table(res.bakta_json)),
               "seconds": time.perf_counter() - t0}  # fmt: skip
    finally:
        shutil.rmtree(root, ignore_errors=True)
    _log(rec)
    return rec


@cache
def gc_record(i: int, size: str = "1M") -> dict[str, Any]:
    """Attack GC: the de novo arm of replicate ``i`` with ``gc_poor_peptides``."""
    root, t0 = _root("gc"), time.perf_counter()
    try:
        res = simulate(i, size, root)
        seed, lengths, truth = sim_seed(i), _lengths(res), _truth(res)
        path, info = gc_poor_peptides(res, _C().replicate_rng(GC_BASE, i))
        dn = run_gippa(res, path, root / "gc" / "OUT", seed)
        g4 = g4_counts(dn, truth, "denovo_subgroup", lengths)
        rec = {"stream": GC_BASE, "replicate": i, "seed": seed, "size": size,
               "gc": {k: v for k, v in g4.items() if ".intergenic." in k},
               "injected": len(info["keys"]["intergenic"]), "n_peptides": info["n_peptides"],
               "mean_gc": info["mean_gc"], "poor_gc_max": info["poor_gc_max"],
               "entered": entered(dn, info["keys"], lengths),
               "seconds": time.perf_counter() - t0}  # fmt: skip
    finally:
        shutil.rmtree(root, ignore_errors=True)
    _log(rec)
    return rec


@cache
def rbs_record(i: int, size: str = "1M") -> dict[str, Any]:
    """Attack b's mutant arm: the database arm of replicate ``i`` at ``rbs_min_quantile`` 0.50."""
    root, t0 = _root("rbs"), time.perf_counter()
    try:
        res = simulate(i, size, root)
        seed, lengths, truth = sim_seed(i), _lengths(res), _truth(res)
        db = run_gippa(res, Path(res.outdir) / "peptides.tsv", root / "rbs" / "OUT", seed,
                      sets=[f"novel_orfs.rbs_min_quantile={RBS_MUTANT_QUANTILE}"])  # fmt: skip
        rec = {"stream": f"{GATE_BASE}.rbs", "replicate": i, "seed": seed, "size": size,
               "g3_by_class": g3_counts(db, truth, lengths),
               "seconds": time.perf_counter() - t0}  # fmt: skip
    finally:
        shutil.rmtree(root, ignore_errors=True)
    _log(rec)
    return rec


@cache
def cluster_record(i: int, size: str = "1M") -> dict[str, Any]:
    """The clustered-null arm of replicate ``i`` (lead decision 7): the database arm on
    ``inject_null_entries``; attacks a and d and the arm's filtered G4 (recorded, not gated)."""
    G = _G()
    root, t0 = _root("cl"), time.perf_counter()
    try:
        res = simulate(i, size, root)
        seed, lengths, truth = sim_seed(i), _lengths(res), _truth(res)
        path, info = inject_null_entries(res, _C().replicate_rng(CL_BASE, i))
        db = run_gippa(res, path, root / "cl" / "OUT", seed)
        a = G.match_novel(G.g4_predictions(db.novel_orfs, db.side, None, regime="picked_decoy",
                                           lengths=lengths), truth, lengths)  # fmt: skip
        targets = G.list_targets(db.side, "picked_decoy")
        pside, missing = pooled_side(db, res)
        rec = {"stream": CL_BASE, "replicate": i, "seed": seed, "size": size,
               "cls": CL_CLASS, "m": CL_ENTRIES, "k": list(CL_PEPTIDES),
               "g4": g4_counts(db, truth, "picked_decoy", lengths),
               "attack_a": {c: _fdr_counts(a[c], targets[c]) for c in G.NOVEL_CLASSES},
               "attack_d": g4_counts(db, truth, "picked_decoy", lengths, side=pside),
               "attack_d_unrebuilt": missing,
               "n_entries": len(info["entries"]), "n_eligible": info["n_eligible"],
               "n_target_entries": info["n_target_entries"], "n_rows": info["n_rows"],
               "entered": entered(db, info["keys"], lengths),
               "diag": _diag(db, _tsv(path), truth, lengths, gene_table(res.bakta_json)),
               "seconds": time.perf_counter() - t0}  # fmt: skip
    finally:
        shutil.rmtree(root, ignore_errors=True)
    _log(rec)
    return rec
