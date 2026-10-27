"""Independent oracle for M7-12 round 2 (lead decisions 11-16 of
docs/specs/M7-12-sim-variant-novel-scenarios.md, after the statistician's block on PR #225).

Like tests/helpers/m7_sim_oracle.py this imports nothing from gippasim: every expectation is
recomputed from the written FASTA/JSON/TSV files.

* ``sample_frame_alt_sites`` recomputes ``truth/alt_sites.tsv`` the way decision 11 says: the SNV
  codon is taken in the **sample** gene's frame (``bakta_like.json`` on ``genome.fasta``) and lifted
  to the reference by the local offset ``pos_reference - pos_sample``; ``ref_aa`` is the reference
  translation of the lifted codon, ``alt_aa`` the sample translation; synonymous rows and rows with
  a stop on either side are dropped. SNV codons near an indel, at a gene's first codon or in a gene
  crossing the origin are *ambiguous*: the oracle reports them so the tests leave them out on both
  sides.
* ``indels_in_genes`` lists every ins/del with a gene, with ``frameshift`` and ``junction``.
* The six-frame helpers check decision 12's wrong peptides on the genome.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from helpers import m7_sim_oracle as O

AltKey = tuple[str, int, str, str, str, str]  # contig, codon_start, strand, gene, ref_aa, alt_aa
SiteKey = tuple[str, int, str]  # contig, codon_start (1-based, reference), gene


@dataclass
class FrameTruth:
    """Expected alt_sites rows, ambiguous reference codons and expected ``frameshift_upstream``."""

    want: set[AltKey] = field(default_factory=set)
    ambiguous: set[SiteKey] = field(default_factory=set)
    #: per expected site: (any frameshifting indel upstream in the gene, net shift mod 3)
    shift: dict[SiteKey, tuple[bool, int]] = field(default_factory=dict)
    #: per expected site: the strain row of (one of) its SNV(s)
    snv: dict[SiteKey, dict[str, str]] = field(default_factory=dict)


def indels_in_genes(outdir: Path) -> list[dict[str, Any]]:
    """Every ins/del row of ``strain_differences.tsv`` with a gene, plus ``frameshift`` (length
    change not a multiple of 3) and ``junction`` (decision 11's junction: the 1-based reference
    coordinate of the first reference base after the VCF anchor, ``pos_reference + 1``)."""
    out = []
    for r in O.rows(outdir / "truth" / "strain_differences.tsv"):
        if r["type"] in ("ins", "del") and r["gene"]:
            d = len(r["alt"]) - len(r["ref"])
            out.append({**r, "shift": d, "frameshift": d % 3 != 0,
                        "junction": int(r["pos_reference"]) + 1})  # fmt: skip
    return out


def _gene_span(g: Mapping[str, Any], n: int) -> tuple[int, int] | None:
    s0, e0 = int(g["start"]) - 1, int(g["stop"])
    if e0 <= s0 or e0 > n:
        return None  # crosses the origin
    return s0, e0


def sample_frame_alt_sites(outdir: Path, near: int = 3) -> FrameTruth:
    ref = O.read_fasta(outdir / "reference.fasta")
    smp = O.read_fasta(outdir / "genome.fasta")
    genes = O.bakta_features(outdir / "bakta_like.json")  # sample annotation
    indels = indels_in_genes(outdir)
    all_indels = [r for r in O.rows(outdir / "truth" / "strain_differences.tsv")
                  if r["type"] in ("ins", "del")]  # fmt: skip
    ft = FrameTruth()
    for r in O.rows(outdir / "truth" / "strain_differences.tsv"):
        if r["type"] != "snv" or not r["gene"]:
            continue
        g, c = genes[r["gene"]], r["contig"]
        ps, pr = int(r["pos_sample"]) - 1, int(r["pos_reference"]) - 1
        span = _gene_span(g, len(smp[c]))
        if span is None:
            k0 = O.ref_codon_start(g, pr + 1, len(ref[c]))  # best guess; excluded anyway
            ft.ambiguous |= {(c, k0 + d, r["gene"]) for d in range(-3, 4)}
            continue
        s0, e0 = span
        cs = s0 + 3 * ((ps - s0) // 3) if g["strand"] == "+" else e0 - 3 * ((e0 - 1 - ps) // 3 + 1)
        rc = cs + (pr - ps)
        key: SiteKey = (c, rc + 1, r["gene"])
        first = s0 if g["strand"] == "+" else e0 - 3
        close = any(
            ix["contig"] == c
            and abs(int(ix["pos_sample"]) - 1 - ps) <= max(len(ix["ref"]), len(ix["alt"])) + near
            for ix in all_indels
        )
        if cs == first or close:
            ft.ambiguous.add(key)
            continue
        cod_r = O.strand_slice(ref[c], rc, rc + 3, g["strand"])
        cod_s = O.strand_slice(smp[c], cs, cs + 3, g["strand"])
        ra, aa = O.CODE11.get(cod_r, "?"), O.CODE11.get(cod_s, "?")
        if ra == aa or "*" in (ra, aa):
            continue
        ft.want.add((c, rc + 1, g["strand"], r["gene"], ra, aa))
        up = [ix for ix in indels if ix["gene"] == r["gene"] and ix["contig"] == c and (
              int(ix["pos_sample"]) - 1 < ps if g["strand"] == "+"
              else int(ix["pos_sample"]) - 1 > ps)]  # fmt: skip
        fs = [ix for ix in up if ix["frameshift"]]
        ft.shift[key] = (bool(fs), sum(int(ix["shift"]) for ix in up) % 3)
        ft.snv[key] = r
    return ft


# ---------------------------------------------------------------- six-frame (decision 12)


def translate_nt(nt: str) -> str:
    return "".join(O.CODE11.get(nt[i : i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def codon_aa(seq: str, c0: int, strand: str) -> str | None:
    """Code-11 residue of the codon whose lowest 0-based base is ``c0`` (None off a contig end)."""
    if c0 < 0 or c0 + 3 > len(seq):
        return None
    return O.CODE11.get(O.strand_slice(seq, c0, c0 + 3, strand), "X")


@dataclass(frozen=True)
class Span:
    contig: str
    s1: int  # 1-based closed
    e1: int
    strand: str


def cds_spans(outdir: Path) -> list[Span]:
    """Annotated CDS/sORFs of the sample annotation (origin-crossing ``end`` unwrapped)."""
    seqs = O.read_fasta(outdir / "genome.fasta")
    out = []
    for g in O.bakta_features(outdir / "bakta_like.json").values():
        s1, e1 = int(g["start"]), int(g["stop"])
        if e1 < s1:
            e1 += len(seqs[g["contig"]])
        out.append(Span(g["contig"], s1, e1, g["strand"]))
    return out


def coding_spans(outdir: Path) -> list[Span]:
    """``cds_spans`` plus every row of ``truth/novel_orfs.tsv``: the reading frames a wrong
    *non-coding* peptide must avoid."""
    out = cds_spans(outdir)
    for r in O.rows(outdir / "truth" / "novel_orfs.tsv"):
        out.append(Span(r["contig"], int(r["start"]), int(r["end"]), r["strand"]))
    return out


def in_frame_with(p: Span, g: Span) -> bool:
    if p.contig != g.contig or p.strand != g.strand or p.e1 < g.s1 or p.s1 > g.e1:
        return False
    return (p.s1 - g.s1) % 3 == 0 if p.strand == "+" else (g.e1 - p.e1) % 3 == 0


def wrong_class(p: Span, cds: list[Span]) -> str:
    """intergenic (no CDS overlap), cds_outframe (same-strand CDS overlap, other frame) or
    antisense (opposite-strand CDS overlap only)."""
    hits = [g for g in cds if g.contig == p.contig and not (p.e1 < g.s1 or p.s1 > g.e1)]
    if not hits:
        return "intergenic"
    if any(g.strand == p.strand for g in hits):
        return "cds_outframe"
    return "antisense"
