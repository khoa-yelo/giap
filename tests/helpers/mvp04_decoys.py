"""Oracles for MVP-04 part 2: ``decoys = "scored_null"`` (docs/specs/MVP-04-ecoli-evidence-run.md,
pipeline steps 5 and 6, examples 2-4, Q2, lead decisions 1 and 8; M7-12 decisions 21, 23, 24).

Nothing here imports gippasim. The views read ``peptides.tsv``, ``truth/peptides.tsv``,
``truth/proteomics.json`` and ``provenance.json`` of one ``simulate_on_genome`` output directory.

Vocabulary:

* a **null peptide** is a truth row with ``is_wrong = true`` and ``origin != decoy``. Its **kind**
  is ``sixframe`` when ``origin = sixframe_noncoding`` and ``protein`` otherwise (the shuffled
  source attributed to a Bakta locus, step 5);
* a **null PSM** is one (null peptide, run) pair. Rule 6: it is written exactly once, either as
  the target row (the peptide, its entry ID) or as the decoy row of its entry (the reversed
  peptide, ``rev_<entry>``), with the q / PEP / score of the same null draw;
* lead decision 9: a protein-entry null is a real Trypsin/P product of a recorded locus with no
  true peptide, truth ``origin = PROTEIN_NULL_ORIGIN`` (``undetected_product``), its truth row(s)
  placed like correct rows (``locus_tag``, ``contig``, 1-based closed ``start``/``end``,
  ``strand`` on the sample); so it maps to the genome and can give a gene false evidence;
* the **entry** of a null peptide is the Bakta locus tag (protein kind) or the makedb six-frame
  ``GIPPAORF_{class}_{contig}_{start1}_{end1}_{p|m}`` ID (six-frame kind).

Lead decision 10 (statistician review of PR #247, B1): the protein / six-frame null split is a
fixed 50/50; the search-space share is still recorded in ``provenance.json``.

``gene_frames`` builds the C9-shaped frames that ``gippa.genetable.detection.pair_decoys`` and
``choose_tier`` take; ``decoy_wins`` is the §7.3.2 step-3 picked competition (key: q asc, PEP asc,
score desc; a tie goes to the decoy) between each gene's best paired decoy and its best target that
GIPPA's mapping would accept (a correct peptide of that gene, or a target peptide that occurs in the
gene's recorded protein with I = L), both at q <= ``q_value_max``.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from scipy import stats

ORF_DECOY = re.compile(
    r"^rev_GIPPAORF_(intergenic|cds_outframe|antisense|nterm_extension)_(.+)_(\d+)_(\d+)_([pm])$"
)
ORF_ID = re.compile(
    r"^GIPPAORF_(intergenic|cds_outframe|antisense|nterm_extension)_(.+)_(\d+)_(\d+)_([pm])$"
)
Q_VALUE_MAX = 0.01  # input.q_value_max default (spec "The run")
DECOY_PREFIXES = ("rev_", "REV_", "DECOY_", "decoy_", "XXX_", "##")  # input.decoy_prefixes default
PROTEIN_NULL_ORIGIN = "undetected_product"  # lead decision 9 (MVP-04 part 2)
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
    truth: list[dict[str, str]] = field(default_factory=list)  # its truth/peptides.tsv rows

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
            nulls.setdefault(r["peptide"], Null(r["peptide"], kind, r["origin"])).truth.append(r)
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


# ---------------------------------------------------------------- GIPPA detection frames


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


def gippa_tier(v: View, genes: Iterable[str]) -> Tier:
    """GIPPA's own ``pair_decoys`` + ``choose_tier`` on the simulated peptides (defaults of
    ``input.decoy_prefixes``, ``input.q_value_max`` and ``MIN_DECOY_PAIRING``)."""
    from gippa.genetable.detection import choose_tier, pair_decoys, protein_gene_map

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
    that names the gene and that GIPPA would map to it (a correct peptide of the gene, or a
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
    """``gippa makedb --six-frame-orfs`` on (annotation, genome): the entries by ``db_id`` and the
    distinct tryptic products of the protein entries and of the six-frame (``GIPPAORF_``) entries
    not also produced by a protein entry; ``protein_share`` = protein / (protein + six-frame)."""
    from gippa.makedb.build import MakedbOptions, build_strain_db

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
        (six if k.startswith("GIPPAORF_") else prot).update(tryptic(s))
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


# ---------------------------------------------------------------- decision 9: protein nulls

_BASES = "TCAG"
_AA11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE11 = {a + b + c: _AA11[16 * i + 4 * j + k] for i, a in enumerate(_BASES)
          for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}  # fmt: skip
_COMP = str.maketrans("ACGT", "TGCA")
EXT_BP = 3000  # a placement may lie this far outside the recorded span (N-terminal extension)


def _seg(seq: str, s1: int, e1: int, strand: str) -> str:
    nt = seq[s1 - 1 : e1]
    return nt.translate(_COMP)[::-1] if strand == "-" else nt


def _tr(nt: str) -> str:
    return "".join(CODE11.get(nt[i : i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def read_seq(path: Path) -> str:
    """The one sequence of a plain FASTA (``truth/sample_genome.fna``)."""
    return "".join(x.strip() for x in path.read_text().splitlines()[1:]).upper()


def locus_proteins(out: Path, recorded: Mapping[str, Mapping[str, Any]], sample: str
                   ) -> dict[str, list[str]]:  # fmt: skip
    """``locus -> sample proteins``: the recorded span translated on the sample (initiator read
    as M), plus the N-terminally extended protein when ``truth/novel_orfs.tsv`` has an
    ``nterm_extension`` row sharing the locus's stop. Expressible, non-origin-crossing loci."""
    ext: dict[tuple[str, int], str] = {}
    for r in rows(out / "truth" / "novel_orfs.tsv"):
        if r["class"] == "nterm_extension":
            key = (r["strand"], int(r["end"]) if r["strand"] == "+" else int(r["start"]))
            ext[key] = r["aa"]
    res: dict[str, list[str]] = {}
    for tag, f in recorded.items():
        s1, e1, st = int(f["start"]), int(f["stop"]), str(f["strand"])
        if not expressible(f) or e1 < s1:
            continue
        aa = ("M" + _tr(_seg(sample, s1, e1, st))[1:]).rstrip("*")
        res[tag] = [aa]
        e = ext.get((st, e1 if st == "+" else s1))
        if e:
            res[tag].append(e.rstrip("*"))
    return res


def true_loci(v: View) -> set[str]:
    return {g for loci in v.correct_loci.values() for g in loci}


def candidate_loci(v: View, proteins: Mapping[str, list[str]]) -> dict[str, int]:
    """Decision 9's draw pool: ``locus -> number of its undetected products`` (distinct
    Trypsin/P products of its sample proteins, I = L, that are no correct peptide) for every
    expressible recorded locus with no true peptide and at least one such product."""
    seen = {il(p) for p in v.correct}
    done = true_loci(v)
    out = {}
    for tag, prots in proteins.items():
        if tag in done:
            continue
        k = len({il(x) for aa in prots for x in tryptic(aa)} - seen)
        if k:
            out[tag] = k
    return out


def protein_null_violations(
    v: View, recorded: Mapping[str, Mapping[str, Any]], proteins: Mapping[str, list[str]],
    sample: str, genome: str,
) -> dict[str, list[str]]:  # fmt: skip
    """Decision 9 per protein-entry null, by rule: ``origin``; ``product`` (a Trypsin/P product
    of a sample protein of its attributed locus, an expressible recorded locus); ``undetected``
    (the locus has no true peptide and the peptide is no correct peptide, I = L); ``maps`` (a
    truth row places it on the attributed locus: the sample segment translates to it, initiator
    M allowed, and the genome GIPPA gets differs in at most 2 residues, SNVs only)."""
    bad: dict[str, list[str]] = {k: [] for k in ("origin", "product", "undetected", "maps")}
    done = true_loci(v)
    seen = {il(p) for p in v.correct}
    for n in v.nulls.values():
        if n.kind != "protein":
            continue
        if n.origin != PROTEIN_NULL_ORIGIN:
            bad["origin"].append(f"{n.peptide}: origin {n.origin}")
        ent = sorted(n.entries())
        tag = ent[0] if len(ent) == 1 else ""
        f = recorded.get(tag)
        if f is None or not expressible(f):
            bad["product"].append(f"{n.peptide}: entry {ent} is no expressible recorded locus")
            continue
        if not any(n.peptide in tryptic(aa) for aa in proteins.get(tag, [])):
            bad["product"].append(f"{n.peptide}: no Trypsin/P product of {tag} on the sample")
        if tag in done or il(n.peptide) in seen:
            bad["undetected"].append(f"{n.peptide}: {tag} has a true peptide / peptide detected")
        ok = False
        for r in n.truth:
            if r["locus_tag"] != tag or r["contig"] in ("", "NA"):
                continue
            s1, e1, st = int(r["start"]), int(r["end"]), r["strand"]
            lo, hi = int(f["start"]) - EXT_BP, int(f["stop"]) + EXT_BP  # N-term extensions
            if st != f["strand"] or s1 < lo or e1 > hi:
                continue
            on_s, on_g = _tr(_seg(sample, s1, e1, st)), _tr(_seg(genome, s1, e1, st))
            if len(on_s) != len(n.peptide):
                continue
            same = on_s == n.peptide or ("M" + on_s[1:] == n.peptide)
            diff = sum(a != b for a, b in zip(on_g[1:], n.peptide[1:], strict=True))
            if same and diff <= 2:
                ok = True
        if not ok:
            bad["maps"].append(f"{n.peptide}: no truth placement on {tag} that maps")
    return bad


NTERM = "GIPPAORF_nterm_extension_"


def _accepted(r: Mapping[str, str]) -> bool:
    return finite(r["q_value"]) and fnum(r["q_value"]) <= Q_VALUE_MAX


def nterm_gene(entry: str, genes: Sequence[Any]) -> str | None:
    """The eligible gene a ``GIPPAORF_nterm_extension_*`` makedb entry overlaps on its strand, if
    exactly one (``genes``: objects with ``locus``, ``eligible``, ``strand`` and 0-based half-open
    ``start`` / ``end``, as ``mvp04_data.Gene``)."""
    *_, s1, e1, st = entry.split("_")
    strand = "+" if st == "p" else "-"
    lo, hi = int(s1) - 1, int(e1)
    hit = [
        g.locus for g in genes if g.eligible and g.strand == strand and g.start < hi and lo < g.end
    ]
    return hit[0] if len(hit) == 1 else None


@dataclass(frozen=True)
class Evidence:
    """One accepted (q <= ``Q_VALUE_MAX``) row giving gene-level evidence."""

    gene: str
    is_decoy: bool
    key: tuple[float, float, float]  # M4-24 shared key: q asc, PEP asc, score desc
    peptide: str
    channel: str  # correct | protein_null | nterm_null


def _identical(loci: Sequence[str], proteins: Mapping[str, list[str]]) -> bool:
    return len({tuple(proteins.get(g, [])) for g in loci}) == 1


def gene_evidence(
    v: View, genes: Sequence[Any], proteins: Mapping[str, list[str]]
) -> list[Evidence]:
    """§7.3.2 gene-level evidence per eligible gene, mirroring ``gippa.genetable.detection`` without
    running the mapping (statistician review of PR #247, B3(b); MVP-04 lead decision 11):

    * **correct** peptides give target evidence to their truth loci: one locus (``unique``), or
      several whose sample proteins are identical (``multi_identical``); ``shared`` gives none;
    * **protein-entry nulls** (decision 9) give target evidence to every eligible gene whose
      sample protein holds the peptide (I = L), when that is one gene or several with identical
      proteins (``multi_identical``); their ``rev_<locus>`` decoys give decoy evidence;
    * **six-frame nulls on a ``GIPPAORF_nterm_extension_*`` entry** place in the gene's
      ``nterm_extension`` locus class, which assigns: target evidence for the one eligible gene
      the entry overlaps on its strand; their ``rev_GIPPAORF_nterm_extension_*`` decoys pair with
      that gene only when the entry was voted to it, i.e. the peptide has an accepted target row
      (``protein_gene_map`` step 3, M7-12 decision 28);
    * other six-frame nulls and their decoys reach no gene (M7-12 decision 23).
    """
    eligible = {g.locus for g in genes if g.eligible}
    prot_nulls = {n.peptide for n in v.nulls.values() if n.kind == "protein"}
    nterm: dict[str, str] = {}
    voted: set[str] = set()
    for n in v.nulls.values():
        ent = sorted(n.entries())
        if n.kind == "sixframe" and len(ent) == 1 and ent[0].startswith(NTERM):
            g = nterm_gene(ent[0], genes)
            if g is not None:
                nterm[n.peptide] = g
                if any(_accepted(r) for r in n.targets.values()):
                    voted.add(n.peptide)
    hay = {g: [il(a) for a in aa] for g, aa in proteins.items() if g in eligible}
    null_loci: dict[str, list[str]] = {}
    out: list[Evidence] = []
    for r in v.peptides:
        if not _accepted(r):
            continue
        k = _key(r)
        pep = r["peptide"]
        if r["is_decoy"] == "true":
            back = pep[::-1]
            if r["protein_ids"].startswith("rev_" + NTERM):
                gs, ch = ([nterm[back]] if back in voted else []), "nterm_null"
            else:
                gs = [x[len("rev_") :] for x in ids(r["protein_ids"])
                      if x.startswith("rev_") and not x.startswith("rev_GIPPAORF_")]  # fmt: skip
                ch = "protein_null"
            out += [Evidence(g, True, k, pep, ch) for g in gs if g in eligible]
            continue
        if pep in nterm:
            gs, ch = [nterm[pep]], "nterm_null"
        elif pep in prot_nulls:
            if pep not in null_loci:
                p = il(pep)
                locs = sorted(g for g, aa in hay.items() if any(p in a for a in aa))
                null_loci[pep] = locs if locs and _identical(locs, proteins) else []
            gs, ch = null_loci[pep], "protein_null"
        elif pep in v.nulls:
            continue
        else:
            loci = sorted(v.correct_loci.get(pep, set()))
            ok = len(loci) == 1 or (len(loci) > 1 and _identical(loci, proteins))
            gs, ch = (loci if ok else []), "correct"
        out += [Evidence(g, False, k, pep, ch) for g in gs if g in eligible]
    return out


@dataclass(frozen=True)
class Competition:
    gene: str
    has_true: bool
    target_won: bool


def gene_competition(
    v: View, genes: Sequence[Any], proteins: Mapping[str, list[str]]
) -> list[Competition]:
    """§7.3.2 picked competition per eligible gene with evidence (``gene_evidence``): the best
    target key against the best decoy key; a tie goes to the decoy. ``has_true``: the gene has a
    true peptide (a truth row with ``is_wrong = false``, MVP-04 lead decision 11)."""
    best: dict[tuple[str, bool], tuple[float, float, float]] = {}
    for e in gene_evidence(v, genes, proteins):
        best[(e.gene, e.is_decoy)] = min(best.get((e.gene, e.is_decoy), e.key), e.key)
    done = true_loci(v)
    out = []
    for g in sorted({g for g, _ in best}):
        t, d = best.get((g, False)), best.get((g, True))
        won = t is not None and (d is None or t < d)
        out.append(Competition(g, g in done, won))
    return out
