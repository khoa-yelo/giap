"""Left-censored Blom normal scores (GIPPA_SPEC.md §7.3.5(d) step 1, §8.3; issue M5-01).

Censored entries (undetected genes) are tied at the lowest ranks: all n0 of them get the average
rank (n0 + 1) / 2, and every uncensored entry gets n0 + its average rank among the uncensored
values. The score is the Blom normal score Φ⁻¹((r - 3/8) / (n + 1/4)).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike
from scipy.special import ndtri
from scipy.stats import rankdata

__all__ = ["blom_scores"]


def blom_scores(x: ArrayLike, censored_mask: ArrayLike) -> npt.NDArray[np.float64]:
    """Blom normal scores with left-censoring and average ranks for ties.

    Parameters
    ----------
    x : array_like, shape (n,)
        Values to rank (log2 abundance). Values of censored entries are ignored (may be NaN).
    censored_mask : array_like of bool, shape (n,)
        True where the entry is left-censored (undetected).

    Returns
    -------
    numpy.ndarray of float64, shape (n,)
        Φ⁻¹((r - 3/8) / (n + 1/4)) in input order; empty for n = 0.

    Raises
    ------
    ValueError
        If the inputs are not 1-D of equal length, or an uncensored value is NaN or ±inf.
    """
    xa = np.asarray(x, dtype=np.float64)
    mask = np.asarray(censored_mask, dtype=bool)
    if xa.ndim != 1 or mask.ndim != 1 or xa.shape != mask.shape:
        raise ValueError("x and censored_mask must be 1-D arrays of equal length")
    n = xa.size
    if n == 0:
        return np.empty(0, dtype=np.float64)
    detected = ~mask
    if not np.all(np.isfinite(xa[detected])):
        raise ValueError("uncensored values must be finite")
    n0 = int(mask.sum())
    ranks = np.empty(n, dtype=np.float64)
    ranks[mask] = (n0 + 1) / 2.0
    if detected.any():
        ranks[detected] = n0 + rankdata(xa[detected], method="average")
    return np.asarray(ndtri((ranks - 0.375) / (n + 0.25)), dtype=np.float64)
