"""Autocorrelation rows per replicon: Moran's I and join counts under both nulls (M5-05).

GIAP_SPEC.md §7.3.5(c), §8.1, §8.4, §9.2, §9.3; docs/specs/M5-05-autocorrelation.md (lead
decisions 14-35), ADR-0020 with Amendment 1 (the Pearson residual join count) and ADR-0021;
docs/methods/autocorrelation.md "Statistics and nulls".

- **Values** (ADR-0020 decision 1, lead decision 22). Over the n genes of a replicon that have a
  value: y = centred log2 abundance, y_w = y - LOESS(y; span ``WIDE_SPAN``) and
  x = y_w - H y_w with H = ``loess_matrix`` at span ``SPAN`` over the ordinals 0 ... n - 1
  (``period`` = n on a circular replicon). The two spans are module constants, not config keys
  (decision 25): the calibration is conditional on them.
- **Nulls** (ADR-0020 decision 2). One permutation matrix per null and replicon, shared by all
  lags and both statistics. Every replicate permutes y_w and is detrended again with the same H,
  so observed and null are the same function of the data. The gene-label null keeps the
  missingness pattern (values move among the genes with values; detection and ``p_detect`` move
  among all eligible genes); the TU-block null moves values, gaps, d and p together.
- **Join count** (ADR-0020 Amendment 1, decision 17). ``value`` is the raw BB_k; the test
  statistic is the Pearson residual sum, so ``null_mean`` / ``null_lo`` / ``null_hi`` / ``effect``
  are on that scale.
- **Flags.** ``ci_unestimated`` on every row (ADR-0020 decision 4); ``not_estimable`` whenever the
  value or the p is NA; ``few_values`` on the ``morans_i`` rows of a replicon whose ``SPAN`` in
  points is below ``MIN_SPAN_POINTS`` (decisions 21-23: the detrend is degenerate there, and the
  join count never touches it).
- **Rows** are lexicographic on (contig in C1 order, statistic, null_model, lag) (decision 24);
  ``q_value`` is BH within (contig, statistic, null_model) (decision 16).

Pure functions; randomness only from the generator passed in (``giap.rng.for_step``
("context_stats.autocorrelation", seed)), drawn per replicon in C1 order as: the label matrix for
the values, the label matrix for detection, then the TU-block matrix.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd

from giap.config.model import GiapConfig
from giap.context.regions import replicon_index
from giap.stats import autocorr as ac
from giap.stats.multitest import bh
from giap.stats.smooth import loess_detrend, loess_matrix, loess_span_points

__all__ = [
    "AUTOCORR_COLUMNS",
    "METHOD",
    "MIN_SPAN_POINTS",
    "NULLS",
    "SPAN",
    "STATISTICS",
    "WIDE_SPAN",
    "autocorrelation_rows",
    "moran_values",
    "replicon_autocorrelation",
]

F64 = npt.NDArray[np.float64]
I64 = npt.NDArray[np.int64]

#: §9.2 test columns and the §9.3 autocorrelation columns, in table order
AUTOCORR_COLUMNS: tuple[str, ...] = (
    "contig", "statistic", "lag", "null_model", "value", "null_mean", "null_lo", "null_hi",
    "stat", "effect", "ci_low", "ci_high", "p_value", "q_value", "method", "n_perm", "flags",
)  # fmt: skip
#: row order within a replicon (lead decision 24), not the §9.3 enum order
STATISTICS: tuple[str, ...] = ("join_count_bb", "morans_i")
NULLS: tuple[str, ...] = ("gene_label", "tu_block")
METHOD: dict[str, str] = {
    "morans_i": "moran_i_permutation",
    "join_count_bb": "join_count_permutation",
}
#: ADR-0020 decision 1: the wide pre-trend removed before the detrend (a constant, decision 25)
WIDE_SPAN = 0.30
#: the detrend whose residuals are the values (a constant, decision 25)
SPAN = 0.10
#: below this span in points the detrend is degenerate (M5-01 decision 4, lead decisions 21-23)
MIN_SPAN_POINTS = 4
NOT_ESTIMABLE = "not_estimable"
FEW_VALUES = "few_values"
CI_UNESTIMATED = "ci_unestimated"
_FLOATS = ("value", "null_mean", "null_lo", "null_hi", "effect", "ci_low", "ci_high", "p_value",
           "q_value")  # fmt: skip
_INTS = ("lag", "n_perm")
_NA = math.nan


def _few_values(n_values: int) -> bool:
    """Whether the ``SPAN`` detrend of ``n_values`` points is degenerate (fit = y, x = 0)."""
    return int(loess_span_points(n_values, SPAN)) < MIN_SPAN_POINTS


def _wide_resid(y: F64, circular: bool) -> F64:
    """y_w: the centred values minus their ``WIDE_SPAN`` LOESS trend over the ordinals."""
    n = y.size
    period = float(n) if circular else None
    fit = loess_detrend(
        y - y.mean(), np.arange(n, dtype=float), WIDE_SPAN, circular=circular, period=period
    )
    return np.asarray(fit.resid, dtype=np.float64)


def _smoother(n: int, circular: bool, span_frac: float) -> F64:
    """The one ``span_frac`` smoother H of decision 22 (observed series and every replicate)."""
    return np.asarray(
        loess_matrix(
            np.arange(n, dtype=float),
            span_frac,
            circular=circular,
            period=float(n) if circular else None,
        ),
        dtype=np.float64,
    )


def moran_values(ix: pd.DataFrame, circular: bool, span_frac: float = SPAN) -> F64:
    """Detrended log2-abundance values of a replicon index (ADR-0020 decision 1, decision 22).

    Parameters
    ----------
    ix : pandas.DataFrame
        Replicon index of ``giap.context.regions.replicon_index``; column ``a`` holds the log2
        abundance of the detected genes and NaN elsewhere.
    circular : bool
        Whether the replicon is circular (circular padding and ``period`` = n in the detrend).
    span_frac : float, default ``SPAN``
        Span of the detrend whose residuals are the values; the ``WIDE_SPAN`` pre-trend is fixed.

    Returns
    -------
    numpy.ndarray of float64, shape (len(ix),)
        x at the genes with a value and NaN at the others. Over the n genes with a value:
        y = centred ``a``, y_w = y - LOESS(y; ``WIDE_SPAN``) and x = y_w - H y_w with H the
        ``span_frac`` smoother over the ordinals 0 ... n - 1. The residuals are identically 0
        (never NaN) when the span is below ``MIN_SPAN_POINTS`` points (decision 26): consumers
        test the span, never ``x == 0``.
    """
    a = _abundance(ix)
    has = np.isfinite(a)
    out = np.full(a.size, np.nan)
    n = int(has.sum())
    if n:
        y_w = _wide_resid(a[has], circular)
        out[has] = y_w - _smoother(n, circular, span_frac) @ y_w
    return out


def _abundance(ix: pd.DataFrame) -> F64:
    if not len(ix):
        return np.zeros(0, dtype=np.float64)
    return np.asarray(pd.to_numeric(ix["a"], errors="coerce").to_numpy(float), dtype=np.float64)


def _detrended_replicates(values: F64, h: F64) -> F64:
    """Rows of ``values`` (B, n) minus their H fit: the per-replicate detrend of ADR-0020."""
    return np.asarray(values - values @ h.T, dtype=np.float64)


def _label_replicates(y_w: F64, has: npt.NDArray[np.bool_], orders: I64, h: F64) -> F64:
    """(B, N) values of the gene-label null: y_w permuted among the genes with values."""
    out = np.full((orders.shape[0], has.size), np.nan)
    if y_w.size:
        out[:, has] = _detrended_replicates(y_w[orders], h)
    return out


def _tu_replicates(y_w_full: F64, orders: I64, h: F64) -> F64:
    """(B, N) values of the TU-block null: whole TUs moved, values and gaps together.

    Every replicate holds the same n values (a permutation of the positions), so the one
    smoother ``h`` detrends the compacted rows and the result is scattered back.
    """
    out = np.full(orders.shape, np.nan)
    if not orders.size:
        return out
    moved = y_w_full[orders]
    keep = np.isfinite(moved)
    n = int(keep[0].sum())
    if not n:
        return out
    out[keep] = _detrended_replicates(moved[keep].reshape(-1, n), h).ravel()
    return out


def _na_row(contig: str, statistic: str, null_model: str, lag: int,
            flags: Sequence[str]) -> dict[str, Any]:  # fmt: skip
    """A row with no statistic: every number NA, ``n_perm`` 0 and the given flags."""
    return {"contig": contig, "statistic": statistic, "lag": int(lag), "null_model": null_model,
            "value": _NA, "null_mean": _NA, "null_lo": _NA, "null_hi": _NA, "stat": statistic,
            "effect": _NA, "ci_low": _NA, "ci_high": _NA, "p_value": _NA, "q_value": _NA,
            "method": METHOD[statistic], "n_perm": 0, "flags": tuple(sorted(flags))}  # fmt: skip


def _row(contig: str, statistic: str, null_model: str, lag: int, value: float, observed: float,
         test: ac.LagTest, flags: Sequence[str]) -> dict[str, Any]:  # fmt: skip
    """One estimable row: ``value`` as reported, ``effect`` on the test statistic's scale."""
    return {"contig": contig, "statistic": statistic, "lag": int(lag), "null_model": null_model,
            "value": value, "null_mean": test.null_mean, "null_lo": test.null_lo,
            "null_hi": test.null_hi, "stat": statistic, "effect": observed - test.null_mean,
            "ci_low": _NA, "ci_high": _NA, "p_value": test.p_value, "q_value": _NA,
            "method": METHOD[statistic], "n_perm": int(test.n_valid),
            "flags": tuple(sorted(flags))}  # fmt: skip


def replicon_autocorrelation(
    ix: pd.DataFrame,
    contig: str,
    circular: bool,
    lags: Sequence[int],
    B: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Rows of one replicon: both statistics under both nulls at every lag (``q_value`` NA).

    Parameters
    ----------
    ix : pandas.DataFrame
        Replicon index of ``giap.context.regions.replicon_index`` (eligible genes in genome
        order, columns ``d``, ``p``, ``a``, ``tu_id``).
    contig : str
        Replicon name written to the ``contig`` column.
    circular : bool
        Whether the replicon is circular (lag pairs wrap and the detrend pads circularly).
    lags : sequence of int
        Lags to report; a lag above ``giap.stats.autocorr.max_lag`` has NA statistics and the
        flag ``not_estimable``, so every replicon has the same row grid.
    B : int
        Permutations per null (``stats.permutations``).
    rng : numpy.random.Generator
        Source of all randomness: the label matrix for the values, the label matrix for
        detection, then the TU-block matrix, in that order.

    Returns
    -------
    pandas.DataFrame
        ``AUTOCORR_COLUMNS`` with ``q_value`` NA, ordered by (statistic, null_model, lag).
    """
    lags = [int(k) for k in lags]
    n_all = len(ix)
    a = _abundance(ix)
    has = np.isfinite(a)
    n = int(has.sum())
    d = np.asarray(ix["d"].to_numpy(), dtype=bool) if n_all else np.zeros(0, dtype=bool)
    p = (
        np.asarray(pd.to_numeric(ix["p"], errors="coerce").to_numpy(float), dtype=np.float64)
        if n_all
        else np.zeros(0, dtype=np.float64)
    )
    tu = [str(t) for t in ix["tu_id"]] if n_all else []

    y_w_full = np.full(n_all, np.nan)
    h = np.zeros((0, 0))
    x = np.full(n_all, np.nan)
    if n:
        h = _smoother(n, circular, SPAN)
        y_w_full[has] = _wide_resid(a[has], circular)
        x[has] = y_w_full[has] - h @ y_w_full[has]

    order_x = ac.label_permutations(n, B, rng)
    order_d = ac.label_permutations(n_all, B, rng)
    order_tu = ac.tu_block_permutations(ac.tu_blocks(tu, circular), B, rng)

    moran_obs = ac.morans_i_matrix(x.reshape(1, -1), lags, circular)[0]
    x_label = _label_replicates(y_w_full[has], has, order_x, h)
    x_tu = _tu_replicates(y_w_full, order_tu, h)
    moran_null = {
        "gene_label": ac.morans_i_matrix(x_label, lags, circular),
        "tu_block": ac.morans_i_matrix(x_tu, lags, circular),
    }
    raw = ac.join_count_matrix(d.reshape(1, -1), lags, circular)[0]
    jc_obs = ac.join_count_matrix(d.reshape(1, -1), lags, circular, P=p.reshape(1, -1))[0]
    jc_null = {
        "gene_label": ac.join_count_matrix(d[order_d], lags, circular, P=p[order_d]),
        "tu_block": ac.join_count_matrix(d[order_tu], lags, circular, P=p[order_tu]),
    }

    degenerate = _few_values(n)  # decisions 21-23: the morans_i rows of this replicon are NA
    few = (CI_UNESTIMATED, FEW_VALUES, NOT_ESTIMABLE)
    rows: list[dict[str, Any]] = []
    for statistic in STATISTICS:
        moran = statistic == "morans_i"
        for null_model in NULLS:
            null = moran_null[null_model] if moran else jc_null[null_model]
            for col, lag in enumerate(lags):
                if moran and degenerate:
                    rows.append(_na_row(contig, statistic, null_model, lag, few))
                    continue
                value = float(moran_obs[col]) if moran else float(raw[col])
                observed = value if moran else float(jc_obs[col])
                test = ac.lag_test(observed, null[:, col])
                flags = [CI_UNESTIMATED]
                if math.isnan(value) or math.isnan(test.p_value):
                    flags.append(NOT_ESTIMABLE)
                rows.append(_row(contig, statistic, null_model, lag, value, observed, test, flags))
    return _frame(rows)


def _frame(rows: Sequence[dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(list(rows), columns=list(AUTOCORR_COLUMNS))
    df = df.astype({**dict.fromkeys(_FLOATS, "float64"), **dict.fromkeys(_INTS, "int64")})
    for c in ("contig", "statistic", "null_model", "stat", "method", "flags"):
        df[c] = df[c].astype(object)
    return df


def autocorrelation_rows(
    genes: pd.DataFrame, replicons: pd.DataFrame, cfg: GiapConfig, rng: np.random.Generator
) -> pd.DataFrame:
    """The ``autocorrelation`` table: every replicon in C1 order with its BH ``q_value``.

    Parameters
    ----------
    genes : pandas.DataFrame
        §9.3 ``genes`` columns (``contig``, ``start``, ``locus_tag``, ``partial``, ``tu_id``,
        ``detected``, ``p_detect``, ``abundance_log2``).
    replicons : pandas.DataFrame
        C1 replicons in table order; ``topology`` = ``circular`` makes a replicon circular.
    cfg : giap.config.model.GiapConfig
        ``context.moran_lags`` (an inclusive range) and ``stats.permutations``.
    rng : numpy.random.Generator
        Source of all randomness, consumed replicon by replicon in C1 order.

    Returns
    -------
    pandas.DataFrame
        ``AUTOCORR_COLUMNS``, sorted by (contig in C1 order, statistic, null_model, lag), with
        ``q_value`` = BH over the lags of each (contig, statistic, null_model) family.
    """
    lo, hi = (int(v) for v in cfg.context.moran_lags)
    lags = list(range(lo, hi + 1))
    B = int(cfg.stats.permutations)
    parts = [
        replicon_autocorrelation(
            replicon_index(genes, str(contig)),
            str(contig),
            str(topology) == "circular",
            lags,
            B,
            rng,
        )
        for contig, topology in zip(replicons["contig"], replicons["topology"], strict=True)
    ]
    if not parts:
        return _frame([])
    out = pd.concat(parts, ignore_index=True)
    q = np.full(len(out), np.nan)
    for _, sub in out.groupby(["contig", "statistic", "null_model"], sort=False):
        q[sub.index.to_numpy()] = bh(sub["p_value"].to_numpy(float))
    out["q_value"] = q
    return out
