"""Inputs and a naive oracle for the M5-05 autocorrelation table.

Inputs reuse the hand-built gene tables of ``region_tests_data`` (``Contig`` / ``gene_table``):
the §9.3 ``genes.tsv`` columns ``giap.context.regions.replicon_index`` reads (``contig``,
``start``, ``locus_tag``, ``partial``, ``tu_id``, ``detected``, ``p_detect``,
``abundance_log2``) plus the C1 ``replicons`` frame (``topology`` gives ``circular``).

The oracle recomputes one replicon's rows with the merged engines
(``giap.stats.autocorr``, ``giap.stats.smooth``, ``giap.stats.multitest``) in the order the spec
fixes, so a disagreement is in ``giap.context.autocorrelation``, not in the statistics:

* **values** (docs/specs/M5-05 lead decision 22 + ADR-0020 decision 1): over the n genes with a
  value, y = centred log2 abundance, y_w = y - LOESS(y; span 0.30) and x = y_w - H y_w with
  H = ``loess_matrix`` at span 0.10 over the ordinals 0 ... n - 1 (``period`` = n when circular);
  x sits at the genes with values, NaN elsewhere, and the lag statistics run over the full
  replicon index of N eligible genes;
* **nulls** (ADR-0020 decision 2): every replicate permutes y_w and is detrended again with the
  same H; the gene-label null keeps the missingness pattern (values move among genes with
  values, d and p move among all eligible genes), the TU-block null moves values, gaps, d and p
  together;
* **draws** (spec "RNG"): replicons in C1 order, per replicon the label matrix for x, the label
  matrix for d, then the TU-block matrix, all of ``cfg.stats.permutations`` rows;
* **join count** (ADR-0020 Amendment 1): ``value`` is the raw BB_k, the test statistic is the
  Pearson residual sum, so ``null_mean`` / ``null_lo`` / ``null_hi`` / ``effect`` are on the
  Pearson scale;
* **q_value**: BH within (contig, statistic, null_model) (lead decision 16).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd
import region_tests_data as rd
from region_tests_data import Contig, flag_set, gene_table, is_na  # noqa: F401  (re-exported)

#: ADR-0020 decision 1: the wide pre-trend removed before the 10 % detrend
WIDE_SPAN = 0.30
#: ``context.autocorr.loess_span``: the detrend whose residuals are the values
SPAN = 0.10
#: the detrend is degenerate below this span in points (M5-01 decision 4; k < 4 -> x = 0)
MIN_SPAN_POINTS = 4
STATISTICS = ("join_count_bb", "morans_i")
NULLS = ("gene_label", "tu_block")
METHOD = {"morans_i": "moran_i_permutation", "join_count_bb": "join_count_permutation"}
#: §9.2 test columns plus the §9.3 autocorrelation columns (spec "Public interface")
COLUMNS = (
    "contig",
    "statistic",
    "lag",
    "null_model",
    "value",
    "null_mean",
    "null_lo",
    "null_hi",
    "stat",
    "effect",
    "ci_low",
    "ci_high",
    "p_value",
    "q_value",
    "method",
    "n_perm",
    "flags",
)
NOT_ESTIMABLE = "not_estimable"
FEW_VALUES = "few_values"
CI_UNESTIMATED = "ci_unestimated"


def config(permutations: int = 19, lags: Sequence[int] = (1, 20)) -> Any:
    """A ``GiapConfig`` with ``stats.permutations`` and ``context.moran_lags``."""
    from giap.config.model import GiapConfig

    return GiapConfig.model_validate(
        {"stats": {"permutations": int(permutations)}, "context": {"moran_lags": list(lags)}}
    )


def lags_of(cfg: Any) -> list[int]:
    """``context.moran_lags`` as the inclusive range it denotes."""
    lo, hi = (int(v) for v in cfg.context.moran_lags)
    return list(range(lo, hi + 1))


def rows(genes: pd.DataFrame, replicons: pd.DataFrame, cfg: Any, seed: int = 7) -> pd.DataFrame:
    """``autocorrelation_rows`` with the spec's stream ``context_stats.autocorrelation``."""
    from giap.context.autocorrelation import autocorrelation_rows
    from giap.rng import for_step

    out = autocorrelation_rows(
        genes, replicons, cfg, for_step("context_stats.autocorrelation", seed)
    )
    return out  # type: ignore[no-any-return]


def oracle_rows(
    genes: pd.DataFrame, replicons: pd.DataFrame, cfg: Any, seed: int = 7
) -> pd.DataFrame:
    """The whole table, recomputed from the merged engines (see the module docstring)."""
    from giap.rng import for_step
    from giap.stats.multitest import bh

    rng = for_step("context_stats.autocorrelation", seed)
    lags = lags_of(cfg)
    b = int(cfg.stats.permutations)
    parts = [
        oracle_replicon(genes, str(c), str(t) == "circular", lags, b, rng)
        for c, t in zip(replicons["contig"], replicons["topology"], strict=True)
    ]
    out = pd.concat(parts, ignore_index=True)
    q = np.full(len(out), np.nan)
    for _, sub in out.groupby(["contig", "statistic", "null_model"], sort=False):
        q[sub.index.to_numpy()] = bh(sub["p_value"].to_numpy(float))
    out["q_value"] = q
    return out[list(COLUMNS)]


def span_points(n: int) -> int:
    """The 10 % span in points over ``n`` values (``giap.stats.smooth.loess_span_points``)."""
    from giap.stats.smooth import loess_span_points

    return int(loess_span_points(n, SPAN))


def few_values(n_values: int) -> bool:
    """Lead decisions 21/22: a 10 % span below 4 points makes the detrended values degenerate."""
    return span_points(n_values) < MIN_SPAN_POINTS


def oracle_values(y: Sequence[float] | np.ndarray, circular: bool) -> np.ndarray:
    """x over the n genes with values: the wide pre-trend and the 10 % detrend of decision 22."""
    from giap.stats.smooth import loess_detrend

    yy = np.asarray(y, dtype=float)
    n = yy.size
    if n == 0:
        return yy.copy()
    pos = np.arange(n, dtype=float)
    period = float(n) if circular else None
    centred = yy - yy.mean()
    wide = loess_detrend(centred, pos, WIDE_SPAN, circular=circular, period=period)
    h = smoother(n, circular)
    return np.asarray(wide.resid - h @ wide.resid, dtype=float)


def smoother(n: int, circular: bool) -> np.ndarray:
    """The one 10 % ``loess_matrix`` of decision 22 (observed series and every replicate)."""
    from giap.stats.smooth import loess_matrix

    return np.asarray(
        loess_matrix(
            np.arange(n, dtype=float),
            SPAN,
            circular=circular,
            period=float(n) if circular else None,
        )
    )


def _pre_detrended(y: np.ndarray, circular: bool) -> np.ndarray:
    from giap.stats.smooth import loess_detrend

    n = y.size
    period = float(n) if circular else None
    centred = y - y.mean()
    return np.asarray(
        loess_detrend(
            centred, np.arange(n, dtype=float), WIDE_SPAN, circular=circular, period=period
        ).resid
    )


def oracle_replicon(
    genes: pd.DataFrame,
    contig: str,
    circular: bool,
    lags: Sequence[int],
    b: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Rows of one replicon, drawing the three permutation matrices in the spec's order."""
    from giap.context.regions import replicon_index
    from giap.stats import autocorr as ac

    ix = replicon_index(genes, contig)
    n_all = len(ix)
    a = ix["a"].to_numpy(float) if n_all else np.zeros(0)
    has = np.isfinite(a)
    n = int(has.sum())
    d = ix["d"].to_numpy(bool) if n_all else np.zeros(0, dtype=bool)
    p = ix["p"].to_numpy(float) if n_all else np.zeros(0)
    tu = [str(t) for t in ix["tu_id"]]

    y_w = np.full(n_all, np.nan)
    x = np.full(n_all, np.nan)
    h = np.zeros((0, 0))
    if n:
        h = smoother(n, circular)
        y_w[has] = _pre_detrended(a[has], circular)
        x[has] = y_w[has] - h @ y_w[has]

    order_x = ac.label_permutations(n, b, rng)
    order_d = ac.label_permutations(n_all, b, rng)
    order_tu = ac.tu_block_permutations(ac.tu_blocks(tu, circular), b, rng)

    x_label = np.full((b, n_all), np.nan)
    if n:
        yl = y_w[has][order_x]
        x_label[:, has] = yl - yl @ h.T
    x_tu = np.full((b, n_all), np.nan)
    for row in range(b):
        moved = y_w[order_tu[row]]
        m = np.isfinite(moved)
        x_tu[row, m] = moved[m] - h @ moved[m]

    moran_obs = ac.morans_i_matrix(x.reshape(1, -1), lags, circular)[0]
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
    degenerate = few_values(n)

    out: list[dict[str, Any]] = []
    for statistic in STATISTICS:
        for null_model in NULLS:
            for col, lag in enumerate(lags):
                if statistic == "morans_i":
                    value = float(moran_obs[col])
                    stat_value = value
                    null = moran_null[null_model][:, col]
                else:
                    value = float(raw[col])
                    stat_value = float(jc_obs[col])
                    null = jc_null[null_model][:, col]
                test = ac.lag_test(stat_value, null)
                na = degenerate and statistic == "morans_i"
                flags = [CI_UNESTIMATED]
                if na:
                    flags += [FEW_VALUES, NOT_ESTIMABLE]
                elif math.isnan(value) or math.isnan(test.p_value):
                    flags.append(NOT_ESTIMABLE)
                out.append(
                    {
                        "contig": contig,
                        "statistic": statistic,
                        "lag": int(lag),
                        "null_model": null_model,
                        "value": math.nan if na else value,
                        "null_mean": math.nan if na else test.null_mean,
                        "null_lo": math.nan if na else test.null_lo,
                        "null_hi": math.nan if na else test.null_hi,
                        "stat": statistic,
                        "effect": math.nan if na else stat_value - test.null_mean,
                        "ci_low": math.nan,
                        "ci_high": math.nan,
                        "p_value": math.nan if na else test.p_value,
                        "q_value": math.nan,
                        "method": METHOD[statistic],
                        "n_perm": 0 if na else int(test.n_valid),
                        "flags": tuple(sorted(flags)),
                    }
                )
    return pd.DataFrame(out, columns=list(COLUMNS))


# --------------------------------------------------------------------------- comparisons


def key_frame(df: pd.DataFrame) -> pd.DataFrame:
    """``df`` indexed by (contig, statistic, null_model, lag) for row-wise comparison."""
    out = df.copy()
    out["lag"] = out["lag"].astype(int)
    return out.set_index(["contig", "statistic", "null_model", "lag"]).sort_index()


def same_value(a: Any, b: Any, tol: float = 1e-9) -> bool:
    """Equal cells, NA == NA; floats within ``tol``."""
    if is_na(a) or is_na(b):
        return bool(is_na(a) and is_na(b))
    if isinstance(a, (int, float, np.floating, np.integer)) and not isinstance(a, bool):
        return bool(math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=tol))
    return bool(a == b)


def differences(got: pd.DataFrame, want: pd.DataFrame, tol: float = 1e-9) -> list[str]:
    """Human-readable cell differences between two autocorrelation frames."""
    g, w = key_frame(got), key_frame(want)
    if list(g.index) != list(w.index):
        return [f"row keys differ: {sorted(set(w.index) ^ set(g.index))[:6]}"]
    bad: list[str] = []
    for col in COLUMNS:
        if col in ("contig", "statistic", "null_model", "lag"):
            continue
        for key, x, y in zip(g.index, g[col], w[col], strict=True):
            ok = flag_set(x) == flag_set(y) if col == "flags" else same_value(x, y, tol)
            if not ok:
                bad.append(f"{key} {col}: got {x!r}, want {y!r}")
    return bad[:12]


def pick(df: pd.DataFrame, contig: str, statistic: str, null_model: str, lag: int) -> pd.Series:
    """The one row of (contig, statistic, null_model, lag)."""
    m = (
        (df["contig"] == contig)
        & (df["statistic"] == statistic)
        & (df["null_model"] == null_model)
        & (df["lag"].astype(int) == int(lag))
    )
    hit = df[m]
    if len(hit) != 1:
        raise AssertionError(f"{len(hit)} rows for {contig}/{statistic}/{null_model}/lag {lag}")
    return hit.iloc[0]


def block(df: pd.DataFrame, contig: str, statistic: str, null_model: str) -> pd.DataFrame:
    """One BH family: the lags of (contig, statistic, null_model), in row order."""
    m = (df["contig"] == contig) & (df["statistic"] == statistic) & (df["null_model"] == null_model)
    return df[m]


# --------------------------------------------------------------------------- generators


def series_contig(
    name: str,
    values: Sequence[float],
    *,
    circular: bool = True,
    detected: Sequence[int] | None = None,
    p: float | Sequence[float] = 0.5,
    tu_size: int = 1,
) -> Contig:
    """A replicon whose detected genes carry ``values`` as ``abundance_log2``.

    ``detected`` defaults to all genes; TUs are runs of ``tu_size`` consecutive genes.
    """
    n = len(values)
    det = [1] * n if detected is None else [int(x) for x in detected]
    probs = [float(p)] * n if isinstance(p, (int, float)) else [float(x) for x in p]
    tus = [f"TU_{name}_{i // tu_size}" for i in range(n)]
    return Contig(name=name, d=det, p=probs, a=list(values), circular=circular, tu=tus)


def ramp(n: int, rng: np.random.Generator | None = None, scale: float = 1.0) -> np.ndarray:
    """Smooth trend + noise: a series whose detrended values are not degenerate."""
    g = np.random.default_rng(0) if rng is None else rng
    t = np.arange(n, dtype=float)
    return 10.0 + 2.0 * np.cos(2 * np.pi * t / n) + scale * g.standard_normal(n)


def genome(contigs: Sequence[Contig]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(genes, replicons) for the given contigs, in C1 order."""
    return rd.gene_table(list(contigs))
