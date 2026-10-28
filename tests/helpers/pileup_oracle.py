"""Builders and an independent oracle for the M7-02 pileup tests (docs/specs/M7-02-*.md).

Written by the test author from the spec text; the oracle does not import ``gippa.variants.pileup``.

* :data:`G1P` / :data:`G2P` / :data:`G3P` are the proteins written into the circular 1,000-nt
  contig ``c1`` of the examples: G1 ``+`` [100, 400) (codon 112 = ``N``, aa_pos 5; an annotated
  ``GTG`` start), G2 ``-`` [500, 800) (codon 785 = ``K``, aa_pos 5) and G3 ``+`` [950, 1100)
  across the origin.
* :func:`gene_hit` builds a ``HitAlignment`` (from M7-01's dataclasses) over residues ``[a, b)``
  of a gene, reading the genome residue except where ``reads`` says otherwise.
* :func:`expected` is the rule 1-7 oracle: per codon key ``(contig, strand, codon mod L)`` it
  returns ``n_ref``, ``n_il``, the ALT counts and the call, counting distinct I/L-folded peptides.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from helpers.align_oracle import locus_nt, place, translate
from helpers.fs_genomes import background, codons_for

L = 1000
IL = frozenset("IL")
DIFF = ("substitution", "swap", "block", "near_block")


def _protein(n: int, overrides: Mapping[int, str], seed: int) -> str:
    filler = "AGSTEKRDQPFYHW"
    out = [filler[(seed + 5 * i) % len(filler)] for i in range(n)]
    for i, a in overrides.items():
        out[i] = a
    return "".join(out)


# G1: index 0 V (GTG start), 1 K, 2 T, 3 A, 4 N (codon 112), 6-7 GG (codons 118, 121),
#     10 L (codon 130), 11 V (codon 133), 13 W (codon 139)
G1P = _protein(
    100, {0: "V", 1: "K", 2: "T", 3: "A", 4: "N", 5: "E", 6: "G", 7: "G", 8: "S", 9: "R",
          10: "L", 11: "V", 12: "E", 13: "W", 14: "K"}, seed=1,
)  # fmt: skip
# G2 (minus): index 4 K (codon 785), 5-6 GG (codons 782, 779)
G2P = _protein(100, {4: "K", 5: "G", 6: "G"}, seed=3)
G3P = _protein(50, {}, seed=5)


def nt(protein: str) -> str:
    return "".join(codons_for(protein))


def build_genome() -> str:
    seq = background(L, seed=11)
    seq = place(seq, 100, "GTG" + nt(G1P[1:]), "+")
    seq = place(seq, 500, nt(G2P), "-")
    seq = place(seq, 950, nt(G3P), "+")
    return seq


SEQ = build_genome()


def codon_nt(seq: str, codon: int, strand: str) -> str:
    """The codon at lowest coordinate ``codon`` (unwrapped allowed), 5'->3' on ``strand``."""
    return locus_nt(seq, codon, codon + 3, strand)


def codon_aa(seq: str, codon: int, strand: str) -> str:
    return translate(codon_nt(seq, codon, strand))


@dataclass(frozen=True)
class GeneSpec:
    tag: str
    start: int
    end: int
    strand: str
    protein: str

    def codon(self, i: int) -> int:
        return self.start + 3 * i if self.strand == "+" else self.end - 3 * (i + 1)


G1 = GeneSpec("G1", 100, 400, "+", G1P)
G2 = GeneSpec("G2", 500, 800, "-", G2P)
G3 = GeneSpec("G3", 950, 1100, "+", G3P)


def frame_of(codons: Sequence[int], strand: str, length: int = L) -> int:
    if strand == "+":
        return 1 + (codons[0] % length) % 3
    end = codons[0] + 3
    return -(1 + ((length - end) % length) % 3)


def kind_of(ref: str, alt: str) -> str:
    if ref == alt:
        return "match"
    if ref in IL and alt in IL:
        return "il"
    return "substitution"


def make_hit(
    align_mod: Any,
    contig: str,
    strand: str,
    ops: Sequence[tuple[str, tuple[int, ...], str, str]],
    frameshift: Any = None,
    length: int = L,
) -> Any:
    """A ``HitAlignment`` from ``(kind, codons, ref, alt)`` ops in peptide order."""
    out, k = [], 0
    for kind, codons, ref, alt in ops:
        out.append(align_mod.ResidueOp(kind, (k, k + len(alt)), tuple(codons), ref, alt))
        k += len(alt)
    peptide = "".join(o.alt for o in out)
    return align_mod.HitAlignment(
        peptide, contig, strand, frame_of(ops[0][1], strand, length), tuple(out), frameshift
    )


def gene_ops(
    g: GeneSpec, a: int, b: int, reads: Mapping[int, str] | None = None, seq: str = SEQ
) -> list[tuple[str, tuple[int, ...], str, str]]:
    reads = reads or {}
    ops = []
    for i in range(a, b):
        c = g.codon(i)
        ref = codon_aa(seq, c, g.strand)
        alt = reads.get(i, ref)
        ops.append((kind_of(ref, alt), (c,), ref, alt))
    return ops


def gene_hit(
    align_mod: Any, g: GeneSpec, a: int, b: int, reads: Mapping[int, str] | None = None
) -> Any:
    return make_hit(align_mod, "c1", g.strand, gene_ops(g, a, b, reads))


# --------------------------------------------------------------------------- oracle


def fold(s: str) -> str:
    return s.replace("I", "L")


@dataclass(frozen=True)
class Sup:
    source: str = "denovo"
    search_db: str = "same"
    q: float = 0.001
    uniqueness: str = "unique"


def _eligible_ops(
    aligns: Iterable[Any], support: Mapping[str, Sup], relation: str, alpha: float
) -> list[tuple[Any, list[Any]]]:
    out = []
    for a in aligns:
        s = support.get(a.peptide)
        if s is None or s.q > alpha or s.uniqueness == "shared":
            continue
        no_alt = relation == "same" and s.source == "database" and s.search_db == "other"
        ops = [o for o in a.ops if not (no_alt and o.kind in DIFF)]
        out.append((a, ops))
    return out


def _category(op: Any) -> tuple[str, str] | None:
    if op.kind == "match":
        return ("REF", "")
    if op.kind == "il":
        return ("IL", "")
    if op.kind == "substitution":
        if op.ref in IL and op.alt in IL:
            return ("IL", "")
        return ("ALT", "J" if op.alt in IL else op.alt)
    return ("ALT", f"{op.ref}>{op.alt}")


def expected(
    aligns: Iterable[Any],
    support: Mapping[str, Sup],
    relation: str = "related",
    alpha: float = 0.01,
    length: int = L,
) -> dict[tuple[str, str, int], dict[str, Any]]:
    """Rules 1-7 per codon key; a folded peptide reading both REF and IL at a codon is IL."""
    cats: dict[tuple[str, str, int], dict[str, set[tuple[str, str]]]] = {}
    for a, ops in _eligible_ops(aligns, support, relation, alpha):
        for op in ops:
            key = (a.contig, a.strand, op.codons[0] % length)
            cat = _category(op)
            assert cat is not None
            cats.setdefault(key, {}).setdefault(fold(a.peptide), set()).add(cat)
    out: dict[tuple[str, str, int], dict[str, Any]] = {}
    for key, per_pep in cats.items():
        n_ref = n_il = 0
        alts: dict[str, int] = {}
        for c in per_pep.values():
            if c == {("REF", "")}:
                n_ref += 1
            elif c <= {("REF", ""), ("IL", "")}:
                n_il += 1
            else:
                assert len(c) == 1, f"oracle: ambiguous peptide categories {c}"
                ((_, k),) = c
                alts[k] = alts.get(k, 0) + 1
        if n_ref and alts:
            call = "mixed"
        elif alts:
            call = "ALT"
        elif n_ref:
            call = "REF"
        else:
            call = "no_call"
        out[key] = {"n_ref": n_ref, "n_il": n_il, "alts": alts, "call": call}
    return out


def alt_counts_text(alts: Mapping[str, int]) -> str:
    items = sorted(alts.items(), key=lambda kv: (-kv[1], kv[0]))
    return ";".join(f"{k}:{n}" for k, n in items)
