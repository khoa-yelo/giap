"""Oracle for M7-12 round 3: lead decisions 21 and 22 of
docs/specs/M7-12-sim-variant-novel-scenarios.md (the statistician's round-2 block on PR #225,
comment 6358601933).

Nothing here imports gippasim. The decoy IDs are checked against what GIPPA itself would build:
``gippa.makedb.build.build_strain_db`` with ``six_frame_orfs = True`` and ``min_orf_aa = 7`` (the
shortest simulated peptide) on the simulator's ``bakta_like.json`` + the genome GIPPA receives.
That is the ``gippa makedb`` database M7-11's picked class FDR pairs decoys against (M7-11 rule 2:
"decoy rows whose stripped ID is ``<decoy_prefix>GIPPAORF_…`` of that entry").

* ``makedb_entries``: every class-tagged makedb entry (six-frame ORFs and ``nterm_extension``
  entries) by ``db_id``, with its 1-based closed unwrapped interval and its protein.
* ``entry_of``: the makedb entry a six-frame product lies in (same contig and strand, inside the
  interval, in its frame).
* ``wrong_class4``: decision 22's four wrong-peptide territories. ``nterm_extension`` when the
  product is in frame with an annotated CDS/sORF, lies entirely upstream of that gene's **true**
  start (``truth/features.tsv`` ``true_start``, else the annotated start) and no in-frame stop
  separates it from that start; otherwise the round-2 rule ``m7_round2_oracle.wrong_class``.
* ``mixture_pep`` / ``mixture_fdr``: the M3-19 two-component normal mixture, written from the
  formulas (for an equal-SD mixture with ``mu1 > mu0`` the FDR is decreasing in the score, so the
  monotone q equals it).
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

from scipy import stats

from helpers import m7_round2_oracle as R
from helpers import m7_sim_oracle as O

CLASSES4 = ("intergenic", "cds_outframe", "antisense", "nterm_extension")
DECOY_ID = re.compile(
    r"^rev_GIPPAORF_(intergenic|cds_outframe|antisense|nterm_extension)_(.+)_(\d+)_(\d+)_([pm])$"
)
SIXFRAME_DECOY_PREFIX = "rev_GIPPAORF_"
MIN_ORF_AA = 7


@dataclass(frozen=True)
class Entry:
    db_id: str
    cls: str
    contig: str
    s1: int  # 1-based closed, unwrapped (e1 may exceed the contig length)
    e1: int
    strand: str
    aa: str
    locus_tag: str | None


def ids(cell: str) -> list[str]:
    return [x for x in cell.split(";") if x]


def makedb_entries(out: Path, workdir: Path) -> dict[str, Entry]:
    """Class-tagged entries of ``gippa makedb --six-frame-orfs`` on the simulator outputs."""
    from gippa.makedb.build import MakedbOptions, build_strain_db

    ref = out / "reference.fasta"
    genome = ref if ref.exists() else out / "genome.fasta"
    ann = out / ("reference_bakta_like.json" if ref.exists() else "bakta_like.json")
    res = build_strain_db(MakedbOptions(annotation=ann, genome=genome, outdir=workdir,
                                        six_frame_orfs=True, min_orf_aa=MIN_ORF_AA,
                                        contaminants=None))  # fmt: skip
    seqs = O.read_fasta(res.fasta)
    entries = {}
    for r in O.rows(res.id_map):
        if not r["db_id"].startswith("GIPPAORF_"):
            continue
        tag = r["locus_tag"] if r["locus_tag"] not in ("", "NA") else None
        entries[r["db_id"]] = Entry(r["db_id"], r["class"], r["contig"], int(r["start"]),
                                    int(r["end"]), r["strand"], seqs[r["db_id"]], tag)  # fmt: skip
    return entries


class EntryIndex:
    """Entries by (contig, strand) for ``entry_of``."""

    def __init__(self, entries: dict[str, Entry], lengths: dict[str, int]) -> None:
        self.by: dict[tuple[str, str], list[Entry]] = {}
        self.lengths = lengths
        for e in entries.values():
            self.by.setdefault((e.contig, e.strand), []).append(e)

    def entry_of(self, contig: str, s1: int, e1: int, strand: str) -> list[Entry]:
        n = self.lengths[contig]
        out = []
        for e in self.by.get((contig, strand), []):
            for d in (0, n):
                a, b = s1 + d, e1 + d
                if e.s1 <= a and b <= e.e1:
                    frame = (a - e.s1) % 3 == 0 if strand == "+" else (e.e1 - b) % 3 == 0
                    if frame:
                        out.append(e)
                        break
        return out


def lengths(out: Path) -> dict[str, int]:
    return {c: len(s) for c, s in O.receiving_genome(out).items()}


# ---------------------------------------------------------------- decision 22 territories


@dataclass(frozen=True)
class Gene:
    contig: str
    strand: str
    true0: int  # + : 0-based first base of the true start codon; - : 0-based end (exclusive)


def true_starts(out: Path) -> list[Gene]:
    """Annotated CDS/sORFs with their true start (truth/features.tsv), origin-crossing skipped."""
    feats = {r["feature_id"]: r for r in O.rows(out / "truth" / "features.tsv")}
    genome = O.receiving_genome(out)
    genes = []
    for tag, g in O.bakta_features(out / "bakta_like.json").items():
        s1, e1, st = int(g["start"]), int(g["stop"]), g["strand"]
        if e1 < s1 or e1 > len(genome[g["contig"]]):
            continue
        ts = feats.get(tag, {}).get("true_start", "")
        if st == "+":
            genes.append(Gene(g["contig"], st, (int(ts) if ts else s1) - 1))
        else:
            genes.append(Gene(g["contig"], st, int(ts) if ts else e1))
    return genes


def nterm_upstream(span: R.Span, genes: list[Gene], genome: dict[str, str]) -> bool:
    """In frame with a gene, wholly upstream of its true start, no in-frame stop in between."""
    seq = genome[span.contig]
    s0, e0 = span.s1 - 1, span.e1
    for g in genes:
        if g.contig != span.contig or g.strand != span.strand:
            continue
        if g.strand == "+":
            if e0 > g.true0 or (g.true0 - s0) % 3:
                continue
            gap = range(e0, g.true0, 3)
        else:
            if s0 < g.true0 or (s0 - g.true0) % 3:
                continue
            gap = range(g.true0, s0, 3)
        if all(O.strand_slice(seq, c, c + 3, g.strand) not in O.STOPS for c in gap):
            return True
    return False


def wrong_class4(span: R.Span, cds: list[R.Span], genes: list[Gene],
                 genome: dict[str, str]) -> str:  # fmt: skip
    if nterm_upstream(span, genes, genome):
        return "nterm_extension"
    return R.wrong_class(span, cds)


# ---------------------------------------------------------------- mixture


def mixture_pep(s: float, pi: float, mix: tuple[float, float, float, float]) -> float:
    mu1, sd1, mu0, sd0 = mix
    l0 = math.log(pi) + float(stats.norm.logpdf(s, mu0, sd0))
    l1 = math.log1p(-pi) + float(stats.norm.logpdf(s, mu1, sd1))
    return math.exp(l0 - max(l0, l1) - math.log1p(math.exp(-abs(l0 - l1))))


def mixture_fdr(s: float, pi: float, mix: tuple[float, float, float, float]) -> float:
    mu1, sd1, mu0, sd0 = mix
    l0 = math.log(pi) + float(stats.norm.logsf(s, mu0, sd0))
    l1 = math.log1p(-pi) + float(stats.norm.logsf(s, mu1, sd1))
    return math.exp(l0 - max(l0, l1) - math.log1p(math.exp(-abs(l0 - l1))))


# ---------------------------------------------------------------- peptides.tsv views


def wrong_peptides(out: Path) -> set[str]:
    return {r["peptide"] for r in O.rows(out / "truth" / "peptides.tsv")
            if O.truthy(r["is_wrong"]) and r["origin"] == "sixframe_noncoding"}  # fmt: skip


def wrong_loci(out: Path) -> dict[str, R.Span]:
    return {r["peptide"]: R.Span(r["contig"], int(r["start"]), int(r["end"]), r["strand"])
            for r in O.rows(out / "truth" / "peptides.tsv")
            if O.truthy(r["is_wrong"]) and r["origin"] == "sixframe_noncoding"}  # fmt: skip


def split_rows(out: Path) -> tuple[list[dict[str, str]], list[dict[str, str]],
                                   list[dict[str, str]]]:  # fmt: skip
    """(wrong six-frame target rows, six-frame decoy rows, protein decoy rows) of peptides.tsv."""
    wrong = wrong_peptides(out)
    t: list[dict[str, str]] = []
    d6: list[dict[str, str]] = []
    dp: list[dict[str, str]] = []
    for r in O.rows(out / "peptides.tsv"):
        if r["is_decoy"] == "true":
            (d6 if r["protein_ids"].startswith(SIXFRAME_DECOY_PREFIX) else dp).append(r)
        elif r["peptide"] in wrong:
            t.append(r)
    return t, d6, dp
