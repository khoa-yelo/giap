"""Independent oracle for the giapsim de novo error tests (issue M3-20; test-author code).

Not part of ``giapsim``: it decides, from the true and the observed peptide alone, whether an
observed de novo read is explained by a claimed error type (GIAP_SPEC.md §11.3 "De novo errors";
§7.2.4 block list), computes monoisotopic peptide masses from its own residue table (Unimod /
GIAP_SPEC.md Appendix B.5) and bins scores into deciles.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np

#: Monoisotopic residue masses (Da), Unimod / GIAP_SPEC.md Appendix B.5 (I and L share a mass).
RESIDUE_MASS: dict[str, float] = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919, "L": 113.08406, "I": 113.08406, "N": 114.04293,
    "D": 115.02694, "Q": 128.05858, "K": 128.09496, "E": 129.04259, "M": 131.04049,
    "H": 137.05891, "F": 147.06841, "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}  # fmt: skip
AA20 = "ACDEFGHIKLMNPQRSTVWY"
#: Unimod monoisotopic deltas used by the modification misreads.
OXIDATION, DEAMIDATED, METHYL = 15.994915, 0.984016, 14.015650

#: GIAP_SPEC.md §7.2.4, verbatim: mass-equivalent block exchanges between residue compositions
#: (both orders of each pair) and the near-equal exchanges.
SPEC_BLOCKS: frozenset[frozenset[str]] = frozenset(
    frozenset({a, b})
    for a, b in (("N", "GG"), ("Q", "AG"), ("GT", "AS"), ("AV", "GL"), ("AN", "GQ"))
)
SPEC_NEAR_BLOCKS: frozenset[frozenset[str]] = frozenset(
    frozenset({a, b}) for a, b in (("W", "AD"), ("W", "EG"), ("R", "GV"))
)
#: Modification misreads (true residue, observed residue): A+Oxidation = S, N+Deamidated = D,
#: G+Methyl = A.
MOD_READS: frozenset[tuple[str, str]] = frozenset({("A", "S"), ("N", "D"), ("G", "A")})

ERROR_TYPES = (
    "substitution", "swap", "block", "near_block", "reversed_stretch", "modification", "random",
)  # fmt: skip


def il(s: str) -> str:
    return s.replace("I", "L")


def comp(s: str) -> str:
    """Sorted composition with I folded into L."""
    return "".join(sorted(il(s)))


def norm_pairs(pairs: Iterable[Sequence[str]]) -> frozenset[frozenset[str]]:
    """Unordered pairs of sorted compositions, e.g. ``("GT", "AS")`` -> ``{"GT", "AS"}``."""
    return frozenset(frozenset({comp(a), comp(b)}) for a, b in pairs)


_BLOCKS = norm_pairs(tuple(p) for p in SPEC_BLOCKS)
_NEAR = norm_pairs(tuple(p) for p in SPEC_NEAR_BLOCKS)


def mass(peptide: str) -> float:
    """Monoisotopic residue mass sum (no water; differences are what the tests use)."""
    return sum(RESIDUE_MASS[c] for c in peptide)


def _diff_positions(a: str, b: str) -> list[int]:
    return [i for i, (x, y) in enumerate(zip(il(a), il(b), strict=True)) if x != y]


def _exchange(true: str, obs: str, table: frozenset[frozenset[str]]) -> bool:
    t, o = il(true), il(obs)
    for i in range(min(len(t), len(o)) + 1):
        if t[:i] != o[:i]:
            break
        for la in (1, 2):
            for lb in (1, 2):
                if la == lb == 1 or len(t) - la != len(o) - lb:
                    continue
                if i + la > len(t) or i + lb > len(o) or t[i + la :] != o[i + lb :]:
                    continue
                if frozenset({comp(t[i : i + la]), comp(o[i : i + lb])}) in table:
                    return True
    return False


def explains(true: str, obs: str, error_type: str) -> bool:
    """Whether ``obs`` is ``true`` after one error of ``error_type`` (I = L throughout)."""
    if error_type == "none":
        return il(true) == il(obs)
    if il(true) == il(obs) or not set(obs) <= set(AA20):
        return False
    if error_type == "random":
        return 7 <= len(obs) <= 20
    if error_type in ("block", "near_block"):
        return _exchange(true, obs, _BLOCKS if error_type == "block" else _NEAR)
    if len(true) != len(obs):
        return False
    d = _diff_positions(true, obs)
    if error_type == "substitution":
        return len(d) == 1
    if error_type == "modification":
        return len(d) == 1 and (true[d[0]], obs[d[0]]) in MOD_READS
    if error_type == "swap":
        return (
            len(d) == 2
            and d[1] == d[0] + 1
            and il(obs[d[0]]) == il(true[d[1]])
            and il(obs[d[1]]) == il(true[d[0]])
        )
    if error_type == "reversed_stretch":
        t, o = il(true), il(obs)
        for k in range(3, 7):
            for i in range(len(t) - k + 1):
                if o == t[:i] + t[i : i + k][::-1] + t[i + k :]:
                    return True
        return False
    raise ValueError(f"unknown error type {error_type!r}")


def decile_bins(scores: np.ndarray, n_bins: int = 10) -> np.ndarray:
    """Pooled quantile bins 0..n_bins-1 (edges at the k/n quantiles, right-closed)."""
    s = np.asarray(scores, dtype=float)
    edges = np.quantile(s, np.arange(1, n_bins) / n_bins)
    out: np.ndarray = np.searchsorted(edges, s, side="left")
    return out
