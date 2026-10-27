"""Hand-built genomes and alignments for M7-18 (docs/specs/M7-18-variants-ambiguous-residues.md).

* ``c1``: linear, 66 nt. G1 ``+`` [0, 36) translates to ``AKCREFGHILWD``; codons A@0 K@3 C@6
  R@9 E@12 F@15 G@18 H@21 I@24 L@27 W@30 D@33.
* ``c2`` (the spec's ``c1'``): circular, 200 nt. G2 ``-`` [100, 136) translates to the same
  protein; residue ``i`` sits at codon key ``136 - 3 * (i + 1)`` (residue 2 at 127).
* ``c3``: circular, 200 nt. G3 ``+`` [190, 226) crosses the origin; residue 3 sits at the codon
  [199, 202), which reads bases 199, 0 and 1.

``overrides`` rewrite the gene-strand bases of single residues (for example the C codon to
``TGA``), so a hand-built op can sit over a stop codon.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from gippa.model import Replicon
from gippa.model.coords import Interval
from gippa.model.genes import Gene, GeneModel

CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
    "I": "ATT", "K": "AAA", "L": "CTT", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA",
    "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}  # fmt: skip
STOP = {"TAA": "*", "TAG": "*", "TGA": "*"}
P = "AKCREFGHILWD"
NONSTANDARD = "BJOUXZ"
IL = frozenset("IL")
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


@dataclass(frozen=True)
class Layout:
    """One gene on one contig: where residue ``i`` sits and what the genome reads there."""

    contig: str
    length: int
    topology: str
    start: int  # unwrapped
    strand: str
    protein: str
    overrides: tuple[tuple[int, str], ...] = ()

    def codon(self, i: int) -> int:
        """Lowest (unwrapped) coordinate of residue ``i``'s codon."""
        if self.strand == "+":
            return self.start + 3 * i
        return self.start + 3 * len(self.protein) - 3 * (i + 1)

    def key(self, i: int) -> int:
        c = self.codon(i)
        return c % self.length if self.topology == "circular" else c

    def gene_nt(self) -> str:
        ov = dict(self.overrides)
        return "".join(ov.get(i, CODON[a]) for i, a in enumerate(self.protein))

    def ref(self, i: int) -> str:
        nt = dict(self.overrides).get(i)
        if nt is None:
            return self.protein[i]
        return STOP.get(nt) or next(a for a, c in CODON.items() if c == nt)

    def sequence(self) -> str:
        plus = self.gene_nt() if self.strand == "+" else revcomp(self.gene_nt())
        n, s = self.length, self.start
        seq = ["C"] * n
        for j, b in enumerate(plus):
            seq[(s + j) % n] = b
        return "".join(seq)

    def replicon(self) -> Replicon:
        return Replicon(self.contig, self.length, self.topology)  # type: ignore[arg-type]

    def gene(self, tag: str = "G1") -> Gene:
        iv = Interval(self.contig, self.start, self.start + 3 * len(self.protein), self.strand)  # type: ignore[arg-type]
        return Gene(tag, "cds", iv, self.protein, None, "ATG")

    def gm(self) -> GeneModel:
        return GeneModel({self.contig: self.replicon()}, (self.gene(),), ())

    def genome(self) -> dict[str, str]:
        return {self.contig: self.sequence()}

    def frame(self) -> int:
        if self.strand == "+":
            return 1 + (self.start % self.length) % 3
        end = self.start + 3 * len(self.protein)
        return -(1 + ((self.length - end) % self.length) % 3)


def linear(protein: str = P, overrides: Mapping[int, str] | None = None) -> Layout:
    return Layout("c1", 66, "linear", 0, "+", protein, tuple(sorted((overrides or {}).items())))


def minus(protein: str = P, overrides: Mapping[int, str] | None = None) -> Layout:
    return Layout("c2", 200, "circular", 100, "-", protein,
                  tuple(sorted((overrides or {}).items())))  # fmt: skip


def origin(protein: str = P, overrides: Mapping[int, str] | None = None) -> Layout:
    return Layout("c3", 200, "circular", 190, "+", protein,
                  tuple(sorted((overrides or {}).items())))  # fmt: skip


def kind_of(ref: str, alt: str) -> str:
    if ref == alt:
        return "match"
    if ref in IL and alt in IL:
        return "il"
    return "substitution"


def hit(
    al: Any, lay: Layout, peptide: str, lo: int = 0, shift: int = 0, peptide_name: str | None = None
) -> Any:
    """A ``HitAlignment`` of ``peptide`` over residues ``[lo, lo + len(peptide))`` of ``lay``,
    one op per residue (match / il / substitution). ``shift`` adds to every codon (``-length``
    gives the plain twin of an unwrapped circular hit)."""
    ops = []
    for k, a in enumerate(peptide):
        i = lo + k
        ref = lay.ref(i)
        ops.append(al.ResidueOp(kind_of(ref, a), (k, k + 1), (lay.codon(i) + shift,), ref, a))
    name = peptide if peptide_name is None else peptide_name
    return al.HitAlignment(name, lay.contig, lay.strand, lay.frame(), tuple(ops), None)


def with_x(protein: str, i: int, letter: str = "X") -> str:
    return protein[:i] + letter + protein[i + 1 :]
