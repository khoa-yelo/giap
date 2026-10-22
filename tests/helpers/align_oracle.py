"""Independent oracle for the M7-01 residue aligner (docs/specs/M7-01-residue-alignment.md).

Written by the test author from the spec text (GIAP_SPEC.md §7.2.4, §7.5.2 step 1 and the M7-01
rules); it does not import ``giap``. It works on letters:

* :func:`locus_aa` translates ``[start, end)`` (unwrapped; read across the origin through
  ``seq + seq``) with table 11, reverse-complementing on ``-``; no I -> L folding, ``X`` for a
  codon with a non-ACGT base.
* :func:`best_ops` enumerates every global alignment of the peptide to the locus translation
  with the §7.2.4 cost-1 operations (substitution, adjacent swap, ``BLOCKS``, ``NEAR_BLOCKS``),
  I folded into L for the cost only, every cell inside the ``±2`` band around diagonal 0, and
  picks the M7-01 rule-4 optimum: lowest cost, then the lexicographically smallest diff-type
  multiset under the ``DIFF_TYPES`` precedence (substitution < swap < block < near_block), then
  the **leftmost** placement: the lexicographically smallest tuple of the differences' peptide
  start positions, in peptide order (test-author reading of lead decision 2; the first
  difference as far left as possible, then the next one, ...).
* :func:`expected_ops` turns that into ``(kind, pep, codons, ref, alt)`` tuples with the rule-8
  codon coordinates.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any

from giapsim.codons import CODE_11

PRECEDENCE: tuple[str, ...] = ("substitution", "swap", "block", "near_block")
RANK = {name: i for i, name in enumerate(PRECEDENCE)}
BAND = 2
_BLOCK_PAIRS = (("N", "GG"), ("Q", "AG"), ("GT", "AS"), ("AV", "GL"), ("AN", "GQ"))
_NEAR_PAIRS = (("W", "AD"), ("W", "EG"), ("R", "GV"))
_RC = str.maketrans("ACGTN", "TGCAN")
UNUSABLE = frozenset("*X")

#: (kind, (pep_i, pep_j), codons, ref, alt) -- the comparable form of a ``ResidueOp``
OpTuple = tuple[str, tuple[int, int], tuple[int, ...], str, str]
#: (kind, pep_i, pep_j, ref_i, ref_j) -- one step of an alignment path
Step = tuple[str, int, int, int, int]


def rc(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def fold(s: str) -> str:
    return s.replace("I", "L")


def _comp(s: str) -> str:
    return "".join(sorted(s))


BLOCK_SET = frozenset(frozenset({_comp(a), _comp(b)}) for a, b in _BLOCK_PAIRS)
NEAR_SET = frozenset(frozenset({_comp(a), _comp(b)}) for a, b in _NEAR_PAIRS)


def translate(nt: str) -> str:
    """Table 11 codon by codon, ``X`` for a non-ACGT codon; a trailing partial codon dropped."""
    out = []
    for i in range(0, len(nt) - 2, 3):
        out.append(CODE_11.get(nt[i : i + 3], "X"))
    return "".join(out)


def locus_nt(seq: str, start: int, end: int, strand: str) -> str:
    nt = (seq + seq)[start:end]
    return rc(nt) if strand == "-" else nt


def locus_aa(seq: str, start: int, end: int, strand: str) -> str:
    return translate(locus_nt(seq, start, end, strand))


def codon_at(i: int, start: int, end: int, strand: str) -> int:
    """Rule 8: lowest coordinate (unwrapped) of the codon of locus residue ``i``."""
    return start + 3 * i if strand == "+" else end - 3 * (i + 1)


def place(seq: str, start: int, nt: str, strand: str) -> str:
    """``seq`` with ``nt`` (5'->3' on ``strand``) written at ``[start, start+len)``, wrapping."""
    body = rc(nt) if strand == "-" else nt
    out = list(seq)
    for k, ch in enumerate(body):
        out[(start + k) % len(seq)] = ch
    return "".join(out)


def _exchange(pseg: str, rseg: str, pairs: frozenset[frozenset[str]]) -> bool:
    key = frozenset({_comp(fold(pseg)), _comp(fold(rseg))})
    return len(key) == 2 and key in pairs


def _steps(p: str, r: str, i: int, j: int) -> list[tuple[int, int, str]]:
    out: list[tuple[int, int, str]] = []
    if i < len(p) and j < len(r) and r[j] not in UNUSABLE:
        out.append((1, 1, "match" if fold(p[i]) == fold(r[j]) else "substitution"))
    if i + 2 <= len(p) and j + 2 <= len(r) and not set(r[j : j + 2]) & UNUSABLE:
        a, b = fold(p[i]), fold(p[i + 1])
        x, y = fold(r[j]), fold(r[j + 1])
        if a != b and a == y and b == x:
            out.append((2, 2, "swap"))
    for pa, rb in ((1, 2), (2, 1), (2, 2)):
        if i + pa <= len(p) and j + rb <= len(r) and not set(r[j : j + rb]) & UNUSABLE:
            ps, rs = p[i : i + pa], r[j : j + rb]
            if _exchange(ps, rs, BLOCK_SET):
                out.append((pa, rb, "block"))
            if (pa, rb) != (2, 2) and _exchange(ps, rs, NEAR_SET):
                out.append((pa, rb, "near_block"))
    return out


def all_alignments(p: str, r: str, kmax: int = 2) -> list[list[Step]]:
    """Every global banded alignment of peptide ``p`` to locus residues ``r`` at cost <= kmax."""
    found: list[list[Step]] = []
    path: list[Step] = []

    def rec(i: int, j: int, cost: int) -> None:
        if cost > kmax or abs(j - i) > BAND:
            return
        if i == len(p) and j == len(r):
            found.append(list(path))
            return
        for a, b, kind in _steps(p, r, i, j):
            path.append((kind, i, i + a, j, j + b))
            rec(i + a, j + b, cost + (kind != "match"))
            path.pop()

    rec(0, 0, 0)
    return found


def _key(path: Sequence[Step]) -> tuple[int, tuple[int, ...], tuple[int, ...]]:
    diffs = [s for s in path if s[0] != "match"]
    return (len(diffs), tuple(sorted(RANK[s[0]] for s in diffs)), tuple(s[1] for s in diffs))


def best_path(p: str, r: str, kmax: int = 2) -> list[Step] | None:
    """The rule-4 optimum (see the module docstring), ``None`` when the cost exceeds ``kmax``."""
    paths = all_alignments(p, r, kmax)
    if not paths:
        return None
    return min(paths, key=_key)


def min_cost(p: str, r: str, kmax: int = 2) -> int | None:
    path = best_path(p, r, kmax)
    return None if path is None else _key(path)[0]


def diff_types(path: Iterable[Step]) -> tuple[str, ...]:
    return tuple(sorted((s[0] for s in path if s[0] != "match"), key=RANK.__getitem__))


def to_ops(
    path: Sequence[Step],
    p: str,
    r: str,
    codons: Sequence[int],
    pep_offset: int = 0,
) -> list[OpTuple]:
    """Op tuples of ``path``; matched letters differing only as I/L are labelled ``il``.

    ``codons[j]`` is the codon coordinate of locus residue ``j``; ``pep_offset`` shifts the
    peptide indices (suffix of a frameshift split).
    """
    out: list[OpTuple] = []
    for kind, i, i2, j, j2 in path:
        ref, alt = r[j:j2], p[i:i2]
        if kind == "match" and ref != alt:
            kind = "il"
        out.append((kind, (i + pep_offset, i2 + pep_offset), tuple(codons[j:j2]), ref, alt))
    return out


def expected_ops(
    peptide: str, seq: str, start: int, end: int, strand: str, kmax: int = 2
) -> list[OpTuple]:
    """Rule-4 ops of a non-frameshift hit ``[start, end)`` on ``strand`` (unwrapped)."""
    r = locus_aa(seq, start, end, strand)
    path = best_path(peptide, r, kmax)
    assert path is not None, f"oracle: {peptide} does not align to {r} within {kmax}"
    codons = [codon_at(j, start, end, strand) for j in range(len(r))]
    return to_ops(path, peptide, r, codons)


def op_tuples(alignment: Any) -> list[OpTuple]:
    """The comparable form of ``HitAlignment.ops``."""
    return [
        (
            op.kind,
            (int(op.pep[0]), int(op.pep[1])),
            tuple(int(c) for c in op.codons),
            op.ref,
            op.alt,
        )
        for op in alignment.ops
    ]


def cost_of(ops: Iterable[OpTuple]) -> int:
    return sum(1 for o in ops if o[0] not in ("match", "il"))
