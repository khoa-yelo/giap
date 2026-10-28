"""MAYU-style gene-level FDR, closed form (GIPPA_SPEC.md §7.3.2, §8.3, Appendix B.3; ADR-0009).

Model: a false peptide lands on gene g with probability proportional to a_g (observable peptides,
floored at 1). Genes are split into equal-count rank bins of a_g; for bin b,
w_b = mean(a_g in b) / Σ_h a_h and p_b = 1 - (1 - w_b)^F(t). The expected number of false genes
solves the fixed point FP_b = (N_b - n_b + FP_b)·p_b, i.e. FP_b = (N_b - n_b)·p_b / (1 - p_b).
FDR(t) = Σ_b FP_b / Σ_b n_b (clipped to 1), and q_g = min over t ≥ best_q_g of FDR(t).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pandas as pd
from numpy.typing import ArrayLike

from gippa.stats.tda import monotone_q

__all__ = ["MayuResult", "mayu_gene_fdr"]

_CURVE_COLUMNS = ["t", "n_peptides", "F", "n_genes", "fp", "fdr", "fdr_2pep"]


#: libm element-wise: numpy's SIMD ``exp`` / ``power`` differ in the last bits between AVX2 and
#: AVX-512 nodes, which the cancellation in 1 - e^-λ(1 + λ) magnifies to ~1e-10 relative.
#: glibc's exp / pow agree across FMA-capable nodes of one glibc build (both kinds here), not
#: unconditionally
_exp = np.frompyfunc(math.exp, 1, 1)
_pow = np.frompyfunc(math.pow, 2, 1)


@dataclass(frozen=True, slots=True)
class MayuResult:
    """Per-gene q (NaN without an accepted peptide), the FDR curve and each gene's size bin."""

    gene_q: npt.NDArray[np.float64]
    curve: pd.DataFrame
    bin_of_gene: npt.NDArray[np.int64]


def _vec(x: ArrayLike, name: str) -> npt.NDArray[np.float64]:
    a = np.array(x, dtype=np.float64, ndmin=1)
    if a.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional, got shape {a.shape}")
    return a


def _rank_bins(a: npt.NDArray[np.float64], n_bins: int) -> npt.NDArray[np.int64]:
    """Equal-count bins by rank of a (ties by input order); min(n_bins, len(a)) bins."""
    n = a.size
    bins = np.zeros(n, dtype=np.int64)
    if n == 0:
        return bins
    k = min(n_bins, n)
    order = np.argsort(a, kind="stable")
    bins[order] = (np.arange(n, dtype=np.int64) * k) // n
    return bins


def _fp(missing: npt.NDArray[np.float64], p: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """FP_b = (N_b - n_b)·p_b / (1 - p_b), with FP_b = N_b - n_b when p_b = 1."""
    one = p >= 1.0
    denom = np.where(one, 1.0, 1.0 - p)
    return np.where(one, missing, missing * p / denom)


def mayu_gene_fdr(
    best_q: ArrayLike,
    peptide_q: ArrayLike,
    peptide_pep: ArrayLike | None,
    sizes: ArrayLike,
    n_bins: int,
) -> MayuResult:
    """MAYU closed-form gene FDR (§7.3.2) with the ≥ 2-peptide variant as ``fdr_2pep``.

    Args:
        best_q: best accepted peptide q per eligible gene (NaN = no accepted peptide).
        peptide_q: q of every distinct accepted peptide (defines P(t)).
        peptide_pep: PEP per accepted peptide, or None to use F(t) = t·|P(t)|.
        sizes: ``n_observable_peptides`` per eligible gene (floored at 1).
        n_bins: number of equal-count size bins (fewer when there are fewer genes).

    Raises:
        ValueError: on length mismatches, ``n_bins < 1`` or partly-NaN PEPs.
    """
    bq = _vec(best_q, "best_q")
    pq = _vec(peptide_q, "peptide_q")
    a = np.maximum(_vec(sizes, "sizes"), 1.0)
    if a.size != bq.size:
        raise ValueError(f"sizes has length {a.size}, best_q has length {bq.size}")
    if n_bins < 1:
        raise ValueError(f"n_bins must be >= 1, got {n_bins}")
    pep: npt.NDArray[np.float64] | None = None
    if peptide_pep is not None:
        pep = _vec(peptide_pep, "peptide_pep")
        if pep.size != pq.size:
            raise ValueError(f"peptide_pep has length {pep.size}, peptide_q has length {pq.size}")
        nan = np.isnan(pep)
        if nan.any() and not nan.all():
            raise ValueError("peptide_pep is partly NaN; use PEPs for all peptides or none")
        if nan.all() and pep.size > 0:
            pep = None

    bins = _rank_bins(a, n_bins)
    k = int(bins.max()) + 1 if bins.size else 0
    big_n = np.bincount(bins, minlength=k).astype(np.float64)
    total = a.sum() if a.size else 1.0
    w = np.bincount(bins, weights=a, minlength=k) / np.maximum(big_n, 1.0) / total

    finite = np.isfinite(bq)
    thresholds = np.unique(bq[finite])
    nt = thresholds.size
    pq_sorted = np.sort(pq[~np.isnan(pq)])
    n_pep = np.searchsorted(pq_sorted, thresholds, side="right").astype(np.int64)
    if pep is not None:
        order = np.argsort(pq, kind="stable")
        cum = np.r_[0.0, np.cumsum(pep[order])]
        f_t = cum[np.searchsorted(pq[order], thresholds, side="right")]
    else:
        f_t = thresholds * n_pep

    # n_b(t): genes in bin b with best_q <= t, via the index of each gene's threshold
    t_idx = np.searchsorted(thresholds, bq[finite])
    counts = np.zeros((nt, k), dtype=np.float64)
    np.add.at(counts, (t_idx, bins[finite]), 1.0)
    n_bt = np.cumsum(counts, axis=0)
    missing = big_n[None, :] - n_bt

    big_f = f_t[:, None]
    # (1 - w)^F with 0^0 = 1, so w = 1 and F = 0 gives p = 0 (no log1p(-1)·0 = NaN)
    p = np.clip(1.0 - _pow(1.0 - w[None, :], big_f).astype(np.float64), 0.0, 1.0)
    lam = big_f * w[None, :]
    p2 = np.clip(1.0 - _exp(-lam).astype(np.float64) * (1.0 + lam), 0.0, 1.0)
    p2 = np.minimum(p2, p)  # P(≥2 hits) ≤ P(≥1 hit); guards rounding

    fp = _fp(missing, p).sum(axis=1)
    fp2 = _fp(missing, p2).sum(axis=1)
    n_genes = n_bt.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        fdr = np.where(
            n_genes > 0, np.minimum(fp / np.where(n_genes > 0, n_genes, 1.0), 1.0), np.nan
        )
        fdr2 = np.where(
            n_genes > 0, np.minimum(fp2 / np.where(n_genes > 0, n_genes, 1.0), 1.0), np.nan
        )

    curve = pd.DataFrame(
        {
            "t": thresholds.astype(np.float64),
            "n_peptides": n_pep,
            "F": f_t.astype(np.float64),
            "n_genes": n_genes.astype(np.int64),
            "fp": fp.astype(np.float64),
            "fdr": fdr.astype(np.float64),
            "fdr_2pep": fdr2.astype(np.float64),
        },
        columns=_CURVE_COLUMNS,
    )

    gene_q = np.full(bq.size, np.nan, dtype=np.float64)
    if nt:
        q_t = monotone_q(fdr)
        gene_q[finite] = q_t[t_idx]
    return MayuResult(gene_q=gene_q, curve=curve, bin_of_gene=bins)
