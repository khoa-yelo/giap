"""Target-decoy q-values and picked FDR (GIAP_SPEC.md §7.2.5, §7.3.2, §8.3, Appendix B.6).

Conventions: higher score = better; NaN and -inf = absent. Thresholds are the distinct observed
scores and counts include every score >= t, so a decoy tied with a target counts against it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike

__all__ = ["PickedResult", "monotone_q", "picked_fdr", "picked_fdr_rescaled", "tda_qvalues"]


def _scores(x: ArrayLike, name: str) -> npt.NDArray[np.float64]:
    """1-D float copy with -inf mapped to NaN (absent)."""
    a = np.array(x, dtype=np.float64, ndmin=1)
    if a.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional, got shape {a.shape}")
    a[np.isneginf(a)] = np.nan
    return a


def monotone_q(fdr_at_sorted: ArrayLike) -> npt.NDArray[np.float64]:
    """Monotone q from FDR at thresholds sorted strictest to loosest.

    q_k = min over j >= k of FDR_j (running minimum from the loose end), clipped to [0, 1].
    """
    f = np.clip(np.asarray(fdr_at_sorted, dtype=np.float64), 0.0, 1.0)
    if f.size == 0:
        return f.reshape(0)
    return np.minimum.accumulate(f[::-1])[::-1]


def tda_qvalues(
    target_scores: ArrayLike,
    decoy_scores: ArrayLike,
    pi: float = 1.0,
    plus_one: bool = True,
) -> npt.NDArray[np.float64]:
    """Separate-search target-decoy q-values, one per target in input order.

    FDR(t) = pi * (D(t) + 1) / T(t) (or pi * D(t) / T(t) when ``plus_one`` is False), where D(t)
    and T(t) count decoys and targets with score >= t; q_i = min over t <= s_i of FDR(t), clipped
    to 1. Absent targets (NaN, -inf) get NaN; absent decoys are ignored.

    Raises:
        ValueError: if ``pi`` is not in (0, 1] or an input is not one-dimensional.
    """
    if not (0.0 < pi <= 1.0):
        raise ValueError(f"pi must lie in (0, 1], got {pi}")
    t = _scores(target_scores, "target_scores")
    d = _scores(decoy_scores, "decoy_scores")
    q = np.full(t.shape, np.nan, dtype=np.float64)
    ok = ~np.isnan(t)
    s = t[ok]
    if s.size == 0:
        return q
    d_sorted = np.sort(d[~np.isnan(d)])
    order = np.argsort(-s, kind="stable")  # best first = strictest threshold first
    s_desc = s[order]
    s_asc = s_desc[::-1]
    n_t = np.float64(s.size - np.searchsorted(s_asc, s_desc, side="left"))
    n_d = d_sorted.size - np.searchsorted(d_sorted, s_desc, side="left")
    fdr = pi * (n_d + (1.0 if plus_one else 0.0)) / n_t
    q_ok = np.empty(s.size, dtype=np.float64)
    q_ok[order] = monotone_q(fdr)
    q[ok] = q_ok
    return q


@dataclass(frozen=True, slots=True)
class PickedResult:
    """Picked target-decoy competition per gene.

    Attributes:
        decoy_wins: per gene, True when the decoy wins (ties and decoy-only genes included).
        q: per gene q-value; NaN where the decoy wins or neither score exists.
        n_target_wins: number of genes won by the target.
        n_decoy_wins: number of genes won by the decoy.
    """

    decoy_wins: npt.NDArray[np.bool_]
    q: npt.NDArray[np.float64]
    n_target_wins: int
    n_decoy_wins: int


def picked_fdr(target_best: ArrayLike, decoy_best: ArrayLike) -> PickedResult:
    """Picked FDR (Savitski 2015): per gene keep max(target, decoy), a tie counts as a decoy.

    Gene q-values are ``tda_qvalues`` over the picked list (FDR = (D + 1)/T, made monotone).

    Raises:
        ValueError: if the inputs differ in length or are not one-dimensional.
    """
    t = _scores(target_best, "target_best")
    d = _scores(decoy_best, "decoy_best")
    if t.shape != d.shape:
        raise ValueError(f"target_best and decoy_best differ in length: {t.size} != {d.size}")
    has_t, has_d = ~np.isnan(t), ~np.isnan(d)
    decoy_wins = has_d & (~has_t | (d >= np.where(has_t, t, -np.inf)))
    target_wins = has_t & ~decoy_wins
    q = np.full(t.shape, np.nan, dtype=np.float64)
    q[target_wins] = tda_qvalues(t[target_wins], d[decoy_wins])
    return PickedResult(
        decoy_wins=decoy_wins,
        q=q,
        n_target_wins=int(target_wins.sum()),
        n_decoy_wins=int(decoy_wins.sum()),
    )


def picked_fdr_rescaled(
    target_best: ArrayLike, decoy_best: ArrayLike, paired_share: float
) -> PickedResult:
    """``picked_fdr`` with ADR-0009's rescaling: each decoy win counts 1/f, f = ``paired_share``.

    FDR(t) = (D(t)/f + 1)/T(t) over the picked list, made monotone; f = 1 gives ``picked_fdr``.

    Raises:
        ValueError: if ``paired_share`` is not in (0, 1], or as ``picked_fdr``.
    """
    if not (0.0 < paired_share <= 1.0):
        raise ValueError(f"paired_share must lie in (0, 1], got {paired_share}")
    res = picked_fdr(target_best, decoy_best)
    win = np.flatnonzero(~np.isnan(res.q))
    if not win.size:
        return res
    t = _scores(target_best, "target_best")[win]
    d = np.sort(_scores(decoy_best, "decoy_best")[res.decoy_wins])
    order = np.argsort(-t, kind="stable")
    t_desc = t[order]
    n_t = win.size - np.searchsorted(np.sort(t), t_desc, side="left")
    n_d = d.size - np.searchsorted(d, t_desc, side="left")
    q_win = np.empty(win.size)
    q_win[order] = monotone_q((n_d / paired_share + 1.0) / n_t)
    q = res.q.copy()
    q[win] = q_win
    return PickedResult(res.decoy_wins, q, res.n_target_wins, res.n_decoy_wins)
