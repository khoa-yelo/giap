"""Moran's I, join counts, gene-label and TU-block permutations, lag tests and block length.

GIAP_SPEC.md §7.3.5(c), §7.3.5(d) step 5, §7.3.6, §8.3, Appendix B.4; docs/specs/M5-02 (lead
decisions 4-13), ADR-0020 (residual join count) and ADR-0021 (block length). All functions are
pure; randomness comes only from the ``np.random.Generator`` passed in by the caller (obtained via
``giap.rng.for_step``).

Lags. Circular: pairs (i, (i + k) mod n), i = 0 ... n - 1, valid for 1 <= k <= (n - 1) // 2.
Linear: pairs (i, i + k), i + k < n, valid for 1 <= k <= n - 1. A statistic at an invalid lag does
not exist (NaN).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike

from giap.stats.permutation import perm_p

__all__ = [
    "LagTest",
    "block_len_prime",
    "choose_block_len",
    "join_count",
    "join_count_matrix",
    "label_permutations",
    "lag_pairs",
    "lag_test",
    "max_lag",
    "morans_i",
    "morans_i_matrix",
    "tu_block_permutations",
    "tu_blocks",
]

F64 = npt.NDArray[np.float64]
I64 = npt.NDArray[np.int64]
B8 = npt.NDArray[np.bool_]

_VAR_RTOL = 1e-24  # Σ z² <= 1e-24 · max(1, max x²) means "no variance"
_CHUNK_CELLS = 1 << 22  # rows per chunk = this // n (results do not depend on it)
_BLOCK_LEN_CAP = 20  # ADR-0021: l' <= 20


# ----------------------------------------------------------------------------- lags


def max_lag(n: int, circular: bool) -> int:
    """Largest valid lag on ``n`` genes: circular ``(n - 1) // 2``, linear ``n - 1``; 0 if n < 2."""
    if n < 2:
        return 0
    return (n - 1) // 2 if circular else n - 1


def _valid_lag(n: int, lag: int, circular: bool) -> bool:
    return 1 <= lag <= max_lag(n, circular)


def lag_pairs(n: int, lag: int, circular: bool) -> tuple[I64, I64]:
    """Index pairs ``(i, j)`` at lag ``lag``: j = i + lag (mod n when circular).

    Returns:
        Two int64 arrays of equal length; both empty if ``lag`` is not in [1, ``max_lag``].
    """
    if not _valid_lag(n, lag, circular):
        empty = np.zeros(0, dtype=np.int64)
        return empty, empty.copy()
    if circular:
        i = np.arange(n, dtype=np.int64)
        return i, (i + lag) % n
    i = np.arange(n - lag, dtype=np.int64)
    return i, i + lag


def _shifted(a: npt.NDArray[np.float64], k: int, circular: bool) -> tuple[F64, F64]:
    """Row-wise left/right operands of the lag-k pair products (``a`` padded when circular)."""
    n = a.shape[1] if not circular else a.shape[1] // 2
    if circular:
        return a[:, :n], a[:, k : k + n]
    return a[:, : n - k], a[:, k:]


def _row_chunks(n_rows: int, n: int) -> list[slice]:
    step = max(1, _CHUNK_CELLS // max(n, 1))
    return [slice(s, min(s + step, n_rows)) for s in range(0, n_rows, step)]


def _as_matrix(X: ArrayLike, name: str) -> F64:
    a = np.asarray(X, dtype=np.float64)
    if a.ndim != 2:
        raise ValueError(f"{name} must be 2-D (B, n), got shape {a.shape}")
    return a


# ----------------------------------------------------------------------------- Moran's I


def _moran_chunk(X: F64, lags: Sequence[int], circular: bool) -> F64:
    rows, n = X.shape
    out = np.full((rows, len(lags)), np.nan)
    valid = [c for c, k in enumerate(lags) if _valid_lag(n, int(k), circular)]
    if not valid:
        return out
    has = ~np.isnan(X)
    cnt = has.sum(axis=1)
    xs = np.where(has, X, 0.0)
    mean = xs.sum(axis=1) / np.maximum(cnt, 1)
    z = np.where(has, X - mean[:, None], 0.0)
    ss = np.einsum("bi,bi->b", z, z)
    mx2 = np.max(xs * xs, axis=1) if n else np.zeros(rows)
    ok = (cnt >= 2) & (ss > _VAR_RTOL * np.maximum(1.0, mx2))
    denom = np.where(ok, ss, 1.0) / np.maximum(cnt, 1)
    complete = bool(has.all())
    zz = np.concatenate([z, z], axis=1) if circular else z
    hh = has.astype(np.float64)
    hh = np.concatenate([hh, hh], axis=1) if circular else hh
    for c in valid:
        k = int(lags[c])
        a, b = _shifted(zz, k, circular)
        num = np.einsum("bi,bi->b", a, b)
        if complete:
            w = np.full(rows, float(a.shape[1]))
        else:
            ha, hb = _shifted(hh, k, circular)
            w = np.einsum("bi,bi->b", ha, hb)
        good = ok & (w > 0)
        out[good, c] = num[good] / w[good] / denom[good]
    return out


def morans_i_matrix(X: ArrayLike, lags: Sequence[int], circular: bool) -> F64:
    """Moran's I of every row of ``X`` at every lag.

    Parameters
    ----------
    X
        (B, n) values; NaN marks a gene without a value (e.g. ``x[orders]`` for permuted series).
    lags
        Lags to evaluate; invalid lags give NaN columns.
    circular
        Whether the gene index is circular.

    Returns
    -------
    numpy.ndarray
        (B, len(lags)) float64; row b equals ``morans_i(X[b], k, circular)`` for each k.
        Rows are processed in chunks; results do not depend on the chunking.
    """
    a = _as_matrix(X, "X")
    rows, n = a.shape
    out = np.full((rows, len(lags)), np.nan)
    for sl in _row_chunks(rows, n):
        out[sl] = _moran_chunk(a[sl], lags, circular)
    return out


def morans_i(x: ArrayLike, lag: int, circular: bool) -> float:
    """Moran's I at lag ``lag``: I_k = ((1/W_k) Σ z_i z_j) / ((1/n) Σ z_i²), z = x - x̄.

    x̄ and n run over entries with values (non-NaN); W_k counts lag-k pairs where both have values.

    Returns
    -------
    float
        NaN if the lag is invalid, n < 2, W_k = 0, or Σ z² <= 1e-24 · max(1, max x²).
    """
    row = np.asarray(x, dtype=np.float64).reshape(1, -1)
    return float(morans_i_matrix(row, [lag], circular)[0, 0])


# ----------------------------------------------------------------------------- join counts


def join_count_matrix(
    D: ArrayLike, lags: Sequence[int], circular: bool, P: ArrayLike | None = None
) -> F64:
    """Join count BB_k (or the residual statistic) of every row of ``D`` at every lag.

    Parameters
    ----------
    D
        (B, n) boolean detection over eligible genes.
    lags
        Lags to evaluate; invalid lags give NaN columns.
    circular
        Whether the gene index is circular.
    P
        Optional (B, n) detection probabilities, permuted together with ``D`` by the caller.
        Given P, the statistic is Σ over valid lag-k pairs of (d_i - p_i)(d_j - p_j)
        (ADR-0020 decision 3); otherwise Σ d_i d_j.

    Returns
    -------
    numpy.ndarray
        (B, len(lags)) float64; row b equals ``join_count(D[b], k, circular, p=P[b])``.

    Raises
    ------
    ValueError
        If ``D`` is not 2-D or ``P`` has a different shape.
    """
    d = np.asarray(D, dtype=bool)
    if d.ndim != 2:
        raise ValueError(f"D must be 2-D (B, n), got shape {d.shape}")
    r = d.astype(np.float64)
    if P is not None:
        p = _as_matrix(P, "P")
        if p.shape != d.shape:
            raise ValueError(f"P shape {p.shape} differs from D shape {d.shape}")
        r = r - p
    rows, n = r.shape
    out = np.full((rows, len(lags)), np.nan)
    valid = [c for c, k in enumerate(lags) if _valid_lag(n, int(k), circular)]
    for sl in _row_chunks(rows, n):
        rr = np.concatenate([r[sl], r[sl]], axis=1) if circular else r[sl]
        for c in valid:
            a, b = _shifted(rr, int(lags[c]), circular)
            out[sl, c] = np.einsum("bi,bi->b", a, b)
    return out


def join_count(d: ArrayLike, lag: int, circular: bool, p: ArrayLike | None = None) -> float:
    """Join count BB_k = Σ over valid lag-k pairs of d_i d_j (residual Σ (d - p)(d - p) given p).

    Returns
    -------
    float
        The count (0 is a value); NaN if the lag is invalid.
    """
    row = np.asarray(d, dtype=bool).reshape(1, -1)
    prow = None if p is None else np.asarray(p, dtype=np.float64).reshape(1, -1)
    return float(join_count_matrix(row, [lag], circular, P=prow)[0, 0])


# ----------------------------------------------------------------------------- permutations


def label_permutations(
    n: int, B: int, rng: np.random.Generator, movable: ArrayLike | None = None
) -> I64:
    """Gene-label null: ``B`` uniformly random orders of 0 ... n - 1.

    Parameters
    ----------
    n, B
        Series length and number of permutations.
    rng
        Generator supplying all randomness.
    movable
        Optional boolean mask of length n: only movable positions are permuted among
        themselves; the others stay fixed (identity).

    Returns
    -------
    numpy.ndarray
        (B, n) int64 orders; apply as ``x[orders]``.

    Raises
    ------
    ValueError
        If ``n`` or ``B`` is negative or ``movable`` has the wrong length.
    """
    if n < 0 or B < 0:
        raise ValueError(f"n and B must be >= 0, got n={n}, B={B}")
    out = np.tile(np.arange(n, dtype=np.int64), (B, 1))
    if movable is None:
        mov = np.arange(n, dtype=np.int64)
    else:
        m = np.asarray(movable, dtype=bool).ravel()
        if m.size != n:
            raise ValueError(f"movable has length {m.size}, expected {n}")
        mov = np.flatnonzero(m).astype(np.int64)
    if mov.size < 2 or B == 0:
        return out
    out[:, mov] = rng.permuted(np.tile(mov, (B, 1)), axis=1)
    return out


def tu_blocks(tu_labels: Sequence[str], circular: bool) -> I64:
    """Block id per position: maximal runs of equal consecutive TU labels, numbered from 0.

    On a circular replicon a run at the end with the label of the run at the start joins it (an
    origin-crossing TU). A TU interrupted by another TU's gene forms several blocks.
    """
    n = len(tu_labels)
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    lab = list(tu_labels)
    change = np.array([lab[i] != lab[i - 1] for i in range(1, n)], dtype=np.int64)
    ids = np.concatenate([[0], np.cumsum(change)]).astype(np.int64)
    if circular and ids[-1] > 0 and lab[0] == lab[-1]:
        ids[ids == ids[-1]] = 0
    return ids


def _block_segments(blocks: I64) -> tuple[I64, I64, I64]:
    """Positions grouped by block in emission order, with per-block offsets and lengths."""
    n = blocks.size
    _, inv = np.unique(blocks, return_inverse=True)
    inv = inv.astype(np.int64).ravel()
    seq = np.argsort(inv, kind="stable").astype(np.int64)
    lens = np.bincount(inv).astype(np.int64)
    offs = np.concatenate([[0], np.cumsum(lens)[:-1]]).astype(np.int64)
    same = inv[seq[1:]] == inv[seq[:-1]]
    gaps = np.flatnonzero(same & (np.diff(seq) != 1))
    if gaps.size:
        b0 = int(inv[0])
        lo, hi = int(offs[b0]), int(offs[b0] + lens[b0])
        g = int(gaps[0]) + 1
        if gaps.size > 1 or not lo < g < hi or int(seq[hi - 1]) != n - 1:
            raise ValueError("every block must be a run of positions (or one origin-crossing run)")
        seq[lo:hi] = np.concatenate([seq[g:hi], seq[lo:g]])
    return seq, offs, lens


def tu_block_permutations(blocks: ArrayLike, B: int, rng: np.random.Generator) -> I64:
    """TU-block null: ``B`` orders placing whole blocks in a uniformly random order.

    Each block keeps its internal order (orientation kept); a wrapped block is emitted as its
    end-part followed by its start-part.

    Parameters
    ----------
    blocks
        Block id per position (``tu_blocks``).
    B
        Number of permutations.
    rng
        Generator supplying all randomness.

    Returns
    -------
    numpy.ndarray
        (B, n) int64 orders.

    Raises
    ------
    ValueError
        If ``B`` is negative or a block is neither a run nor one origin-crossing run.
    """
    if B < 0:
        raise ValueError(f"B must be >= 0, got {B}")
    blk = np.asarray(blocks, dtype=np.int64).ravel()
    n = blk.size
    if n == 0:
        return np.zeros((B, 0), dtype=np.int64)
    seq, offs, lens = _block_segments(blk)
    nb = lens.size
    if nb == 1 or B == 0:
        return np.tile(seq, (B, 1))
    perm = rng.permuted(np.tile(np.arange(nb, dtype=np.int64), (B, 1)), axis=1).ravel()
    seg_len = lens[perm]
    out_start = np.cumsum(seg_len) - seg_len
    src = np.repeat(offs[perm] - out_start, seg_len) + np.arange(B * n, dtype=np.int64)
    return seq[src].reshape(B, n)


# ----------------------------------------------------------------------------- lag test


@dataclass(frozen=True, slots=True)
class LagTest:
    """Observed statistic against its permutation null at one lag.

    Attributes
    ----------
    value
        Observed statistic.
    null_mean, null_lo, null_hi
        Mean and linear 2.5 / 97.5 % quantiles of the valid (non-NaN) null.
    p_value
        Two-sided centred permutation p (+1 rule); NaN if it does not exist.
    n_valid
        Number of non-NaN null values.
    """

    value: float
    null_mean: float
    null_lo: float
    null_hi: float
    p_value: float
    n_valid: int


def lag_test(obs: float, null: ArrayLike) -> LagTest:
    """Test one lag statistic against its null.

    ``p_value`` = ``perm_p(obs, null, "two-sided")``. Does not exist: NaN ``obs`` or no valid null
    value (every field NaN except ``value`` and ``n_valid``); a constant null (no variation, e.g.
    no detected gene) gives ``p_value`` NaN whatever ``obs`` is (lead decision 6).
    """
    o = float(obs)
    z = np.asarray(null, dtype=np.float64).ravel()
    valid = z[~np.isnan(z)]
    n_valid = int(valid.size)
    nan = math.nan
    if math.isnan(o) or n_valid == 0:
        return LagTest(o, nan, nan, nan, nan, n_valid)
    if valid.min() == valid.max():
        c = float(valid[0])
        return LagTest(o, c, c, c, nan, n_valid)
    lo, hi = np.quantile(valid, [0.025, 0.975], method="linear")
    return LagTest(
        o, float(valid.mean()), float(lo), float(hi), perm_p(o, valid, "two-sided"), n_valid
    )


# ----------------------------------------------------------------------------- block length


def choose_block_len(
    z: ArrayLike,
    circular: bool,
    B: int,
    rng: np.random.Generator,
    alpha: float = 0.05,
    max_len: int = 20,
) -> int:
    """Block length l-hat: the first lag whose Moran's I is not significant (§7.3.5(d) step 5).

    For k = 1 ... min(``max_len``, ``max_lag(n)``), ``lag_test`` of Moran's I of ``z`` under ``B``
    gene-label permutations (one permutation matrix shared by every lag, so the draws depend on
    B and n only); l-hat = the first k with raw two-sided p > ``alpha`` or NaN. If every lag is
    significant, l-hat = min(``max_len``, ``max_lag(n)`` + 1, n). n < 3 gives 1.

    Raises
    ------
    ValueError
        If ``max_len`` < 1.
    """
    if max_len < 1:
        raise ValueError(f"max_len must be >= 1, got {max_len}")
    x = np.asarray(z, dtype=np.float64).ravel()
    n = x.size
    if n < 3:
        return 1
    orders = label_permutations(n, B, rng)
    top = min(max_len, max_lag(n, circular))
    lags = list(range(1, top + 1))
    obs = morans_i_matrix(x.reshape(1, -1), lags, circular)[0]
    null = morans_i_matrix(x[orders], lags, circular)
    for c, k in enumerate(lags):
        p = lag_test(float(obs[c]), null[:, c]).p_value
        if math.isnan(p) or p > alpha:
            return k
    return int(min(max_len, max_lag(n, circular) + 1, n))


def block_len_prime(l_hat: int, T: int) -> int:
    """Bootstrap block length l' = min(20, 2 l-hat + 1, T) (ADR-0021), shared by M5-06 and M5-08.

    Raises
    ------
    ValueError
        If ``l_hat`` < 1 or ``T`` < 1.
    """
    if l_hat < 1 or T < 1:
        raise ValueError(f"l_hat and T must be >= 1, got l_hat={l_hat}, T={T}")
    return int(min(_BLOCK_LEN_CAP, 2 * int(l_hat) + 1, int(T)))
