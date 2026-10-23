"""Oracles for MVP-04 part 2: ``decoys = "scored_null"`` (docs/specs/MVP-04-ecoli-evidence-run.md,
pipeline steps 5 and 6, examples 2-4, Q2, lead decisions 1 and 8; M7-12 decisions 21, 23, 24).

Nothing here imports giapsim. The views read ``peptides.tsv``, ``truth/peptides.tsv``,
``truth/proteomics.json`` and ``provenance.json`` of one ``simulate_on_genome`` output directory.

Vocabulary:

* a **null peptide** is a truth row with ``is_wrong = true`` and ``origin != decoy``. Its **kind**
  is ``sixframe`` when ``origin = sixframe_noncoding`` and ``protein`` otherwise (the shuffled
  source attributed to a Bakta locus, step 5);
* a **null PSM** is one (null peptide, run) pair. Rule 6: it is written exactly once, either as
  the target row (the peptide, its entry ID) or as the decoy row of its entry (the reversed
  peptide, ``rev_<entry>``), with the q / PEP / score of the same null draw;
* the **entry** of a null peptide is the Bakta locus tag (protein kind) or the makedb six-frame
  ``GIAPORF_{class}_{contig}_{start1}_{end1}_{p|m}`` ID (six-frame kind).

``gene_frames`` builds the C9-shaped frames that ``giap.genetable.detection.pair_decoys`` and
``choose_tier`` take; ``decoy_wins`` is the §7.3.2 step-3 picked competition (key: q asc, PEP asc,
score desc; a tie goes to the decoy) between each gene's best paired decoy and its best target that
GIAP's mapping would accept (a correct peptide of that gene, or a target peptide that occurs in the
gene's recorded protein with I = L), both at q <= ``q_value_max``.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from scipy import stats

ORF_DECOY = re.compile(
    r"^rev_GIAPORF_(intergenic|cds_outframe|antisense|nterm_extension)_(.+)_(\d+)_(\d+)_([pm])$"
)
ORF_ID = re.compile(
    r"^GIAPORF_(intergenic|cds_outframe|antisense|nterm_extension)_(.+)_(\d+)_(\d+)_([pm])$"
)
Q_VALUE_MAX = 0.01  # input.q_value_max default (spec "The run")
DECOY_PREFIXES = ("rev_", "REV_", "DECOY_", "decoy_", "XXX_", "##")  # input.decoy_prefixes default
#: Q2 / lead decision 1: below about this many gene-level null PSMs the split falls back to 50/50.
FALLBACK_PSMS = 50
MIN_ORF_AA = 7  # the shortest simulated peptide (as tests/helpers/m7_round3_oracle.py)


def il(s: str) -> str:
    return s.replace("I", "L")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ids(cell: str) -> list[str]:
    return [x for x in cell.split(";") if x]


def fnum(x: str) -> float:
    return float("nan") if x in ("", "NA", "nan") else float(x)


def finite(x: str) -> bool:
    return math.isfinite(fnum(x))


# ---------------------------------------------------------------- views of one output


@dataclass
class Null:
    peptide: str
    kind: str  # protein | sixframe
    origin: str
    targets: dict[str, dict[str, str]] = field(default_factory=dict)  # run -> target row
    decoys: dict[str, dict[str, str]] = field(default_factory=dict)  # run -> decoy row

    def entries(self) -> set[str]:
        """Entry IDs this null peptide is written under (``rev_`` stripped from decoys)."""
        out = {x for r in self.targets.values() for x in ids(r["protein_ids"])}
        for r in self.decoys.values():
            out |= {x[len("rev_") :] if x.startswith("rev_") else x for x in ids(r["protein_ids"])}
        return out


@dataclass
class View:
    out: Path
    peptides: list[dict[str, str]]
    truth: list[dict[str, str]]
    meta: dict[str, Any]
    provenance: dict[str, Any]
    runs: list[str]
    nulls: dict[str, Null]
    correct: set[str]  # correct (true) peptides
    correct_loci: dict[str, set[str]]  # correct peptide -> its truth loci

    @property
    def decoy_rows(self) -> list[dict[str, str]]:
        return [r for r in self.peptides if r["is_decoy"] == "true"]

    @property
    def target_rows(self) -> list[dict[str, str]]:
        return [r for r in self.peptides if r["is_decoy"] == "false"]


def view(out: Path) -> View:
    peps = rows(out / "peptides.tsv")
    truth = rows(out / "truth" / "peptides.tsv")
    meta = json.loads((out / "truth" / "proteomics.json").read_text())
    prov = json.loads((out / "provenance.json").read_text())
    runs = sorted({r["run"] for r in peps})
    nulls: dict[str, Null] = {}
    correct: set[str] = set()
    loci: dict[str, set[str]] = {}
    for r in truth:
        if r["origin"] == "decoy":
            continue
        if r["is_wrong"] == "true":
            kind = "sixframe" if r["origin"] == "sixframe_noncoding" else "protein"
            nulls.setdefault(r["peptide"], Null(r["peptide"], kind, r["origin"]))
        else:
            correct.add(r["peptide"])
            if r["locus_tag"] not in ("", "NA"):
                loci.setdefault(r["peptide"], set()).add(r["locus_tag"])
    by_rev = {n.peptide[::-1]: n for n in nulls.values()}
    for r in peps:
        if r["is_decoy"] == "false" and r["peptide"] in nulls:
            nulls[r["peptide"]].targets[r["run"]] = r
        elif r["is_decoy"] == "true" and r["peptide"] in by_rev:
            by_rev[r["peptide"]].decoys[r["run"]] = r
    return View(out, peps, truth, meta, prov, runs, nulls, correct, loci)


def null_psms(v: View) -> list[tuple[Null, str, dict[str, str], bool]]:
    """Every written null PSM: (null peptide, run, its row, is_decoy)."""
    out = []
    for n in v.nulls.values():
        for run, r in n.targets.items():
            out.append((n, run, r, False))
        for run, r in n.decoys.items():
            out.append((n, run, r, True))
    return out


def accounting_violations(v: View) -> list[str]:
    """Rule 6 bookkeeping: each null PSM is written once (target xor decoy, per run); every decoy
    row is the reversal of a null peptide in the same run (no ``match_wrong`` extras); a null
    peptide keeps one entry across its rows."""
    bad = []
    counts: dict[tuple[str, str, str], int] = {}
    for r in v.peptides:
        k = (r["peptide"], r["run"], r["is_decoy"])
        counts[k] = counts.get(k, 0) + 1
    for n in v.nulls.values():
        for run in v.runs:
            t = counts.get((n.peptide, run, "false"), 0)
            d = counts.get((n.peptide[::-1], run, "true"), 0)
            if t + d != 1:
                bad.append(f"{n.peptide} {run}: {t} target + {d} decoy rows (want exactly 1)")
        if len(n.entries()) != 1:
            bad.append(f"{n.peptide}: written under several entries {sorted(n.entries())}")
    rev = {n.peptide[::-1] for n in v.nulls.values()}
    for r in v.decoy_rows:
        if r["peptide"] not in rev:
            bad.append(f"decoy {r['peptide']} {r['protein_ids']} {r['run']}: not a null PSM")
    return bad


# ---------------------------------------------------------------- annotation


def annotation(path: Path) -> dict[str, dict[str, Any]]:
    """Recorded CDS / sORF features of a (gzipped) Bakta JSON by locus tag."""
    raw = Path(path).read_bytes()
    data = json.loads(gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw)
    return {str(f["locus"]): f for f in data["features"] if f.get("type") in ("cds", "sorf")}


def expressible(f: Mapping[str, Any]) -> bool:
    return not (f.get("truncated") or f.get("pseudo") or f.get("pseudogene"))


# ---------------------------------------------------------------- mixture (M3-19 formulas)


def mixture_pep(s: float, pi: float, mix: tuple[float, ...]) -> float:
    mu1, sd1, mu0, sd0 = mix
    l0 = math.log(pi) + float(stats.norm.logpdf(s, mu0, sd0))
    l1 = math.log1p(-pi) + float(stats.norm.logpdf(s, mu1, sd1))
    return math.exp(l0 - max(l0, l1) - math.log1p(math.exp(-abs(l0 - l1))))


def mixture_q(s: float, pi: float, mix: tuple[float, ...]) -> float:
    """The mixture FDR at threshold ``s``; for an equal-SD mixture with mu1 > mu0 it is decreasing
    in ``s``, so the monotone q equals it."""
    mu1, sd1, mu0, sd0 = mix
    l0 = math.log(pi) + float(stats.norm.logsf(s, mu0, sd0))
    l1 = math.log1p(-pi) + float(stats.norm.logsf(s, mu1, sd1))
    return math.exp(l0 - max(l0, l1) - math.log1p(math.exp(-abs(l0 - l1))))


# ---------------------------------------------------------------- GIAP detection frames


def gene_frames(v: View, genes: Iterable[str]) -> tuple[Any, Any, Any, Any]:
    """(C2-like genes, C9 peptides, empty C11, empty C12) for ``pair_decoys`` / ``choose_tier``."""
    import pandas as pd

    g = pd.DataFrame({"locus_tag": sorted(set(genes))})
    p = pd.DataFrame(
        {
            "peptide": [r["peptide"] for r in v.peptides],
            "protein_ids": [r["protein_ids"] for r in v.peptides],
            "q_value": [fnum(r["q_value"]) for r in v.peptides],
            "pep": [fnum(r["pep"]) for r in v.peptides],
            "score": [fnum(r["score"]) for r in v.peptides],
            "is_decoy": [r["is_decoy"] == "true" for r in v.peptides],
        }
    )
    hits = pd.DataFrame({"protein_id": [], "gene": [], "rank": [], "peptide": [],
                         "uniqueness": []})  # fmt: skip
    return g, p, hits, hits.copy()


@dataclass(frozen=True)
class Tier:
    tier: str
    warnings: tuple[str, ...]
    n_pairing: int
    paired_share: float
    pairing: Any  # the pair_decoys frame


def giap_tier(v: View, genes: Iterable[str]) -> Tier:
    """GIAP's own ``pair_decoys`` + ``choose_tier`` on the simulated peptides (defaults of
    ``input.decoy_prefixes``, ``input.q_value_max`` and ``MIN_DECOY_PAIRING``)."""
    from giap.genetable.detection import choose_tier, pair_decoys, protein_gene_map

    g, p, ph, prh = gene_frames(v, genes)
    pmap = protein_gene_map(g, p, ph, prh)
    pairing = pair_decoys(p, pmap, list(DECOY_PREFIXES), q_value_max=Q_VALUE_MAX)
    tier, warns = choose_tier(p, ph, prh, pairing)
    share = float(pairing["gene"].notna().mean()) if len(pairing) else float("nan")
    return Tier(tier, tuple(warns), len(pairing), share, pairing)


def _key(r: Mapping[str, str]) -> tuple[float, float, float]:
    return (fnum(r["q_value"]), fnum(r["pep"]), -fnum(r["score"]))


def decoy_wins(v: View, proteins: Mapping[str, str], eligible: set[str]) -> dict[str, bool]:
    """``gene -> decoy won`` for every eligible gene with a paired decoy at q <= ``Q_VALUE_MAX``
    (a ``rev_<gene>`` row). The gene's target is the best target row at q <= ``Q_VALUE_MAX``
    that names the gene and that GIAP would map to it (a correct peptide of the gene, or a
    peptide inside its recorded protein, I = L). No target, or a decoy key <= the target key
    (tie -> decoy), is a decoy win."""
    best_d: dict[str, tuple[float, float, float]] = {}
    for r in v.decoy_rows:
        if not finite(r["q_value"]) or fnum(r["q_value"]) > Q_VALUE_MAX:
            continue
        for x in ids(r["protein_ids"]):
            g = x[len("rev_") :] if x.startswith("rev_") else x
            if g in eligible:
                k = _key(r)
                best_d[g] = min(best_d.get(g, k), k)
    hay = {g: il(aa) for g, aa in proteins.items()}
    best_t: dict[str, tuple[float, float, float]] = {}
    for r in v.target_rows:
        if not finite(r["q_value"]) or fnum(r["q_value"]) > Q_VALUE_MAX:
            continue
        pep = r["peptide"]
        for g in ids(r["protein_ids"]):
            if g not in best_d:
                continue
            ok = g in v.correct_loci.get(pep, set()) or il(pep) in hay.get(g, "")
            if ok:
                k = _key(r)
                best_t[g] = min(best_t.get(g, k), k)
    return {g: (g not in best_t) or (d <= best_t[g]) for g, d in best_d.items()}


def target_won_genes(v: View, eligible: set[str]) -> set[str]:
    """Eligible genes whose best accepted target (a correct peptide unique to the gene in the
    truth, at q <= ``Q_VALUE_MAX``) beats the best paired decoy, or that have no paired decoy:
    the T of the picked gene FDR (acceptance 4; 1 / q_floor)."""
    best_t: dict[str, tuple[float, float, float]] = {}
    for r in v.target_rows:
        loci = v.correct_loci.get(r["peptide"], set())
        if len(loci) != 1 or not finite(r["q_value"]) or fnum(r["q_value"]) > Q_VALUE_MAX:
            continue
        g = next(iter(loci))
        if g in eligible:
            k = _key(r)
            best_t[g] = min(best_t.get(g, k), k)
    best_d: dict[str, tuple[float, float, float]] = {}
    for r in v.decoy_rows:
        if not finite(r["q_value"]) or fnum(r["q_value"]) > Q_VALUE_MAX:
            continue
        for x in ids(r["protein_ids"]):
            g = x[len("rev_") :]
            if g in best_t:
                k = _key(r)
                best_d[g] = min(best_d.get(g, k), k)
    return {g for g, t in best_t.items() if g not in best_d or t < best_d[g]}


# ---------------------------------------------------------------- search space (Q2)


def tryptic(protein: str, lo: int = 7, hi: int = 30, missed: int = 1) -> set[str]:
    """Distinct Trypsin/P products (cut after every K/R) of ``protein`` with ``lo``-``hi`` aa."""
    p = protein.rstrip("*")
    cuts = [0, *(i + 1 for i, c in enumerate(p[:-1]) if c in "KR"), len(p)]
    out = set()
    for i in range(len(cuts) - 1):
        for j in range(i + 1, min(i + 2 + missed, len(cuts))):
            if lo <= cuts[j] - cuts[i] <= hi:
                out.add(p[cuts[i] : cuts[j]])
    return out


def makedb_search_space(annotation_json: Path, genome: Path, workdir: Path) -> dict[str, Any]:
    """``giap makedb --six-frame-orfs`` on (annotation, genome): the entries by ``db_id`` and the
    distinct tryptic products of the protein entries and of the six-frame (``GIAPORF_``) entries
    not also produced by a protein entry; ``protein_share`` = protein / (protein + six-frame)."""
    from giap.makedb.build import MakedbOptions, build_strain_db

    workdir.mkdir(parents=True, exist_ok=True)
    raw = Path(annotation_json).read_bytes()
    plain = workdir / "annotation.json"
    plain.write_bytes(gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw)
    res = build_strain_db(MakedbOptions(annotation=plain, genome=genome, outdir=workdir / "db",
                                        six_frame_orfs=True, min_orf_aa=MIN_ORF_AA,
                                        contaminants=None))  # fmt: skip
    seqs: dict[str, list[str]] = {}
    cur = ""
    for line in Path(res.fasta).read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif line.strip():
            seqs[cur].append(line.strip())
    entries = {k: "".join(v) for k, v in seqs.items()}
    prot, six = set(), set()
    for k, s in entries.items():
        (six if k.startswith("GIAPORF_") else prot).update(tryptic(s))
    six -= prot
    return {"entries": entries, "n_protein": len(prot), "n_sixframe": len(six),
            "protein_share": len(prot) / (len(prot) + len(six))}  # fmt: skip


def write_sample_fasta(out: Path, dest: Path) -> Path:
    """``truth/sample_genome.fna`` with the ``[topology=circular]`` tag (the sample replicon is
    circular like K-12), for makedb on the sample (spec step 5)."""
    seqs: dict[str, list[str]] = {}
    cur = ""
    for line in (out / "truth" / "sample_genome.fna").read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif line.strip():
            seqs[cur].append(line.strip())
    dest.write_text("".join(f">{c} sample [topology=circular]\n" + "".join(s) + "\n"
                            for c, s in seqs.items()))  # fmt: skip
    return dest
