"""Benjamini-Hochberg q-values (GIPPA_SPEC.md §8.2, §8.3, Appendix B.1)."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike

__all__ = ["bh"]


def bh(p: ArrayLike) -> npt.NDArray[np.float64]:
    """Benjamini-Hochberg step-up q-values, monotone in p and returned in input order.

    NaN entries stay NaN and do not count toward the number of tests m. Ties get equal q.

    Raises:
        ValueError: if ``p`` is not one-dimensional or any non-NaN value lies outside [0, 1].
    """
    pa = np.asarray(p, dtype=np.float64)
    if pa.ndim != 1:
        raise ValueError(f"p must be one-dimensional, got shape {pa.shape}")
    ok = ~np.isnan(pa)
    vals = pa[ok]
    if np.any((vals < 0.0) | (vals > 1.0)):
        raise ValueError("p values must lie in [0, 1]")
    q = np.full(pa.shape, np.nan, dtype=np.float64)
    m = vals.size
    if m == 0:
        return q
    order = np.argsort(vals, kind="stable")
    ecdf = np.arange(1, m + 1, dtype=np.float64) / m
    raw = vals[order] / ecdf
    q_sorted = np.minimum(np.minimum.accumulate(raw[::-1])[::-1], 1.0)
    q_ok = np.empty(m, dtype=np.float64)
    q_ok[order] = q_sorted
    q[ok] = q_ok
    return q
