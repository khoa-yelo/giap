"""Per-replicate count records of the M7-13a variants gate (docs/specs/M7-13-calibration-gate.md
rules 1, 3, 4, 9, 10; lead decisions 1-2). Test plumbing owned by the test-author: each record
runs one giapsim scenario through GIAP (``helpers.m7_gate_chain.run_chain``) and reduces the
written ``variants.tsv`` with the pure helpers of ``tests/calibration/_m7_gate.py`` (implementer),
so no check reruns GIAP and the harness logs the counts (shards merge from them).

Seeds (rule 1): replicate ``i`` of a stream base simulates with
``int(replicate_rng(base, i).integers(0, 2**31 - 1))``. Bases: ``variants.g1`` (G1),
``variants.g2`` (G2), ``variants.a6`` (A6a, A6b, A6c and their mutants share it: A6c reruns
GIAP on the same simulation with the systematic transform keyed by
``replicate_rng("variants.a6c", i)``), ``variants.ref_arm`` (reference-proteome arm).

Test-side transforms (no simulator change):

* ``systematic_swaps`` (A6c, rule 9): ``n_sites`` protein sites ``(locus_tag, p)`` are drawn among
  the adjacent pairs ``p, p+1`` (residues differing after I/L folding) covered by at least two
  distinct true peptides that have de novo reads; every de novo read whose true peptide covers a
  drawn site is rewritten as its true peptide with residues ``p, p+1`` exchanged (the same swap
  for every read of that site; the first drawn site in sorted order when a read covers several).
* ``reference_proteome_only`` (reference-proteome arm, lead decision 2 / Q8): reference
  identifications whose peptide (I/L folded) is not a substring of a reference-strain protein
  (``reference_bakta_like.json`` CDS/sORF ``aa``) are dropped; a database search against the
  reference proteome can only return such peptides. Substring, not an enzyme digest: giapsim
  peptides carry missed cleavages, so a strict digest would drop true reference peptides too.

Every record is appended as one JSON line to ``<log_dir()>/variants_gate_replicates.jsonl``.
"""

from __future__ import annotations

import json
import tempfile
import time
from collections import defaultdict
from collections.abc import Mapping
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from helpers.m7_gate_chain import GateRun, cleanup, run_chain

G1_BASE, G2_BASE, A6_BASE, A6C_BASE, REF_BASE = (
    "variants.g1", "variants.g2", "variants.a6", "variants.a6c", "variants.ref_arm",
)  # fmt: skip
A6C_SITES = 30
RECORDS_FILE = "variants_gate_replicates.jsonl"
G1_COUNT_KEYS = ("n_pred", "tp_p", "n_truth", "tp_r", "masked")


def _C() -> Any:
    import _criteria

    return _criteria


def _G() -> Any:
    import _m7_gate

    return _m7_gate


def sim_seed(base: str, i: int) -> int:
    return int(_C().replicate_rng(base, i).integers(0, 2**31 - 1))


def _root(tag: str) -> Path:
    return Path(tempfile.mkdtemp(prefix=f"m7gate-{tag}-"))


def _log(record: Mapping[str, Any]) -> None:
    d = _C().log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / RECORDS_FILE).open("a") as fh:
        fh.write(json.dumps(record, default=float) + "\n")


def _tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _fold(s: str) -> str:
    return s.replace("I", "L").replace("J", "L")


def strip_flag(variants: pd.DataFrame, flag: str) -> pd.DataFrame:
    """``variants`` with ``flag`` removed from every row (mutants)."""
    out = variants.copy()
    out["flags"] = [";".join(f for f in str(v).split(";") if f and f != flag)
                    for v in out["flags"]]  # fmt: skip
    return out


def with_confidence(variants: pd.DataFrame, confidence: str) -> pd.DataFrame:
    """``variants`` with every row's confidence set to ``confidence`` (mutants)."""
    out = variants.copy()
    out["confidence"] = confidence
    return out


# --------------------------------------------------------------------------- transforms


def _fasta(path: Path) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name is not None:
            out[name].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def systematic_swaps(
    path: Path, res: Any, rng: np.random.Generator, n_sites: int = A6C_SITES
) -> tuple[Path, dict[str, int]]:
    """The A6c transform (module doc); returns the rewritten de novo file and its counts."""
    out = Path(res.outdir)
    prot = _fasta(out / "proteins.faa")
    tp = _tsv(out / "truth" / "peptides.tsv")
    tp = tp[(tp["origin"] == "gene") & (tp["is_wrong"].str.lower() != "true")]
    locus = dict(zip(tp["peptide"], tp["locus_tag"], strict=True))
    dn = _tsv(path)
    tr = _tsv(out / "truth" / "denovo.tsv")
    true_of = {(r, s): t for r, s, t in zip(tr["run"], tr["spectrum_id"], tr["true_peptide"],
                                            strict=True)}  # fmt: skip
    is_dn = dn["source"] == "denovo"
    trues = [true_of.get((r, s)) for r, s in zip(dn["run"], dn["spectrum_id"], strict=True)]
    loc: dict[str, tuple[str, int]] = {}
    for t in {t for t, d in zip(trues, is_dn, strict=True) if t and d}:
        lt = locus.get(t)
        if lt is None or lt not in prot:
            continue
        off = prot[lt].find(t)
        if off >= 0:
            loc[t] = (lt, off)
    cover: dict[tuple[str, int], set[str]] = defaultdict(set)
    for t, (lt, off) in loc.items():
        seq = prot[lt]
        for p in range(off, off + len(t) - 1):
            if _fold(seq[p]) != _fold(seq[p + 1]):
                cover[(lt, p)].add(t)
    cands = sorted(k for k, v in cover.items() if len(v) >= 2)
    k = min(n_sites, len(cands))
    sites = sorted(cands[j] for j in rng.choice(len(cands), size=k, replace=False)) if k else []
    peptides = dn["peptide"].tolist()
    n_rewritten = 0
    for j, (t, d) in enumerate(zip(trues, is_dn, strict=True)):
        if not d or t not in loc:
            continue
        lt, off = loc[t]
        for slt, p in sites:
            if slt == lt and off <= p and p + 1 < off + len(t):
                q = p - off
                peptides[j] = t[:q] + t[q + 1] + t[q] + t[q + 2 :]
                n_rewritten += 1
                break
    dn["peptide"] = peptides
    dst = path.with_name("peptides_denovo.a6c.tsv")
    dn.to_csv(dst, sep="\t", index=False)
    return dst, {"n_sites": len(sites), "n_candidates": len(cands), "n_rewritten": n_rewritten}


def reference_proteome_only(path: Path, res: Any) -> tuple[Path, int]:
    """The reference-proteome transform (module doc); returns the file and the dropped count."""
    js = json.loads(Path(res.bakta_json).read_text())  # the reference strain's (run_chain)
    proteome = "|".join(_fold(str(f["aa"])) for f in js["features"]
                        if f.get("type") in ("cds", "sorf") and f.get("aa"))  # fmt: skip
    ref = _tsv(path)
    keep = [_fold(p) in proteome for p in ref["peptide"]]
    dst = path.with_name("reference_peptides.refproteome.tsv")
    ref[keep].to_csv(dst, sep="\t", index=False)
    return dst, int(len(keep) - sum(keep))


# --------------------------------------------------------------------------- records


def _truth(run: GateRun) -> tuple[pd.DataFrame, pd.DataFrame]:
    from giapsim.truth import read_alt_sites, read_indel_sites

    return read_alt_sites(run.outdir), read_indel_sites(run.outdir)


def _high_other(variants: pd.DataFrame) -> dict[str, int]:
    hi = variants[(variants["confidence"] == "high") & (variants["alt_type"] != "substitution")]
    return {str(k): int(v) for k, v in hi["alt_type"].value_counts().sort_index().items()}


@cache
def g1_record(i: int, size: str) -> dict[str, Any]:
    """G1 counts of replicate ``i`` (``variant_f1`` at ``size``, de novo + reference peptides,
    ``genome_relation = related``): the gated set, the recorded non-``mass_equivalent`` subset,
    the mutant set (every substitution row, any confidence), the high non-substitution rows."""
    from giapsim import scenarios as sc

    G = _G()
    seed, root, t0 = sim_seed(G1_BASE, i), _root("g1"), time.perf_counter()
    run = run_chain(sc.variant_f1(seed, size), seed, root, overrides=["genome_relation=related"])
    alt, ind = _truth(run)
    pred = G.g1_predictions(run.variants)
    mutant = G.g1_predictions(with_confidence(run.variants, "high"))
    rec = {
        "stream": G1_BASE, "replicate": i, "seed": seed, "size": size,
        "g1": G.match_variants(pred, alt, ind),
        "g1_non_mass_equivalent": G.match_variants(pred[~pred["mass_equivalent"]], alt, ind),
        "mutant_any_confidence": G.match_variants(mutant, alt, ind),
        "high_other": _high_other(run.variants),
        "seconds": time.perf_counter() - t0,
    }  # fmt: skip
    rec["masked"] = rec["g1"]["masked"]
    cleanup(root)
    _log(rec)
    return rec


@cache
def g2_record(i: int) -> dict[str, Any]:
    """G2 counts of replicate ``i``: ``same_strain_insilico`` database peptides (``search_db =
    same``), no reference peptides, ``genome_relation = same``."""
    from giapsim import scenarios as sc

    G = _G()
    seed, root, t0 = sim_seed(G2_BASE, i), _root("g2"), time.perf_counter()
    run = run_chain(sc.same_strain_insilico(seed), seed, root, peptides="peptides.tsv",
                    reference=None, overrides=["genome_relation=same"])  # fmt: skip
    v = run.variants
    rec = {"stream": G2_BASE, "replicate": i, "seed": seed,
           "violations": int(G.g2_violations(v)),
           "high": int((v["confidence"] == "high").sum()),
           "start_codon_high": int(G.start_codon_high(v)),
           "seconds": time.perf_counter() - t0}  # fmt: skip
    cleanup(root)
    _log(rec)
    return rec


def _a6_chain(i: int, tag: str, rewrite: Mapping[str, Any] | None = None) -> tuple[GateRun, Path]:
    from giapsim import scenarios as sc

    seed, root = sim_seed(A6_BASE, i), _root(tag)
    run = run_chain(sc.isobaric_attack(seed), seed, root, overrides=["genome_relation=same"],
                    rewrite=rewrite)  # fmt: skip
    return run, root


@cache
def a6_record(i: int) -> dict[str, Any]:
    """A6a/A6b counts of replicate ``i`` (``isobaric_attack``, de novo + reference peptides,
    ``genome_relation = same``) and the counts of their mutants: A6a with ``mass_equivalent``
    stripped, A6b with every swap row taken as ``high``, G2's check on these error-laden reads."""
    G = _G()
    t0 = time.perf_counter()
    run, root = _a6_chain(i, "a6")
    v = run.variants
    truth = _tsv(run.outdir / "truth" / "denovo.tsv")
    sp = run.site_peptides()
    rec = {
        "stream": A6_BASE, "replicate": i, "seed": sim_seed(A6_BASE, i),
        "a6a": G.a6a_unflagged(v, sp, truth),
        "a6a_mutant": G.a6a_unflagged(strip_flag(v, "mass_equivalent"), sp, truth),
        "a6b": G.a6_high_swaps(v),
        "a6b_mutant": G.a6_high_swaps(with_confidence(v, "high")),
        "g2_mutant_violations": int(G.g2_violations(v)),
        "seconds": time.perf_counter() - t0,
    }  # fmt: skip
    cleanup(root)
    _log(rec)
    return rec


@cache
def a6c_record(i: int) -> dict[str, Any]:
    """A6c counts of replicate ``i``: the A6 simulation with the systematic swap transform."""
    G = _G()
    t0 = time.perf_counter()
    info: dict[str, int] = {}

    def rewrite(path: Path, res: Any) -> Path:
        dst, counts = systematic_swaps(path, res, _C().replicate_rng(A6C_BASE, i))
        info.update(counts)
        return dst

    run, root = _a6_chain(i, "a6c", {"peptides_denovo.tsv": rewrite})
    rec = {"stream": A6C_BASE, "replicate": i, "seed": sim_seed(A6_BASE, i),
           "a6c": G.a6_high_swaps(run.variants), "transform": info,
           "seconds": time.perf_counter() - t0}  # fmt: skip
    cleanup(root)
    _log(rec)
    return rec


@cache
def ref_arm_record(i: int) -> dict[str, Any]:
    """Reference-proteome arm of replicate ``i``: ``variant_f1`` 100k, de novo reads, reference
    identifications restricted to the reference proteome; ``{class: [rate, n]}`` as written."""
    from giapsim import scenarios as sc
    from helpers.class_error_chain import measured_rates

    seed, root, t0 = sim_seed(REF_BASE, i), _root("ref"), time.perf_counter()
    dropped: list[int] = []

    def rewrite(path: Path, res: Any) -> Path:
        dst, n = reference_proteome_only(path, res)
        dropped.append(n)
        return dst

    run = run_chain(sc.variant_f1(seed, "100k"), seed, root,
                    overrides=["genome_relation=related"],
                    rewrite={"reference_peptides.tsv": rewrite})  # fmt: skip
    rates = measured_rates(run.variants)
    rec = {"stream": REF_BASE, "replicate": i, "seed": seed, "n_dropped": dropped[0],
           "rates": {f"{c[0]}|{c[1]}": [r, n] for c, (r, n) in sorted(rates.items())},
           "seconds": time.perf_counter() - t0}  # fmt: skip
    cleanup(root)
    _log(rec)
    return rec
