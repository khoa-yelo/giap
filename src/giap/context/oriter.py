"""Origin-to-terminus gradient of one genome: oriC, replichores and descriptive slopes (M5-08).

GIAP_SPEC.md §7.3.5(e) and §9.3; docs/specs/M5-08-oriter.md (lead decisions 31-35 and the
decisions of 2026-10-21); ADR-0021 decision 1 **with Amendment 1** (the bootstrap block length
l' = min(20, 2 l-hat + 1, T); the validity domain is a dependence regime, and the
``uncalibrated`` flag is dropped -- no row here ever carries it); docs/methods/oriter.md.

- **Replicons.** Only circular chromosomes are analysed (lead decision 33): plasmids, linear
  replicons and ``unknown``-kind contigs get no row and NA gene columns.
- **oriC.** The C15 ``features`` rows with ``type = oric`` on the contig; the midpoint of each is
  (start + end) / 2 mod L (``end`` may be unwrapped, §6.5). One call gives its midpoint; several
  give their circular mean, but only when every pair lies within ``ORIC_SPREAD_FRAC`` * L
  circular distance -- farther apart the replicon is ``no_oric`` (lead decision 31). Both
  branches log ``context_stats.oriter.multiple_oric`` (decisions 4-5).
- **Geometry.** d = (m - o) mod L; d <= L/2 is the ``right`` replichore at position d / (L/2),
  otherwise ``left`` at (L - d) / (L/2). Positions lie in [0, 1]: 0 at the oriC, 1 at the
  terminus.
- **Points.** One point per TU with a value: y = mean ``abundance_log2`` of its detected eligible
  genes (decision 34), x = the replichore position of the TU midpoint.
- **Fits.** OLS y = a + b x over the pooled, right and left TUs, in that row order (decision 3).
  ``value`` = ``effect`` = b (log2 per replichore), ``p_value`` / ``q_value`` / ``lag`` and the
  null columns NA, flag ``descriptive``: gene placement is evolved, so there is no null here.
- **CI.** 95 % percentile interval of the slopes refitted on ``stats.bootstrap`` circular
  block-bootstrap index replicates of the fit's own TU sequence (pooled in genome order, an arm
  from oriC to terminus). The block length is clamped **per fit** (decision 1): an integer
  ``context.islands.block_len`` as min(config, T_s), else ``block_len_prime`` of the one l-hat
  estimated from the pooled fit's residuals and the fit's own T_s. Fewer than B/2 valid
  replicates drops the CI and adds ``not_estimable``.

Pure functions; randomness only from the generator passed in
(``giap.rng.for_step("context_stats.oriter", seed)``), consumed replicon by replicon in C1 order
as: the l-hat permutations (``auto`` block length and an estimable pooled fit only), then the
bootstrap replicates of the pooled, right and left fits, skipping a fit that does not exist
(decision 2).
"""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd

from giap.config.model import GiapConfig
from giap.context.autocorrelation import AUTOCORR_COLUMNS
from giap.context.regions import replicon_index
from giap.stats.autocorr import block_len_prime, choose_block_len
from giap.stats.bootstrap import circular_block_bootstrap, percentile_ci

__all__ = [
    "ARMS",
    "GENE_COLUMNS",
    "METHOD",
    "NULL_MODEL",
    "ORIC_SPREAD_FRAC",
    "STATISTIC",
    "gene_replichores",
    "oric_position",
    "oriter_rows",
    "replichore",
]

F64 = npt.NDArray[np.float64]
#: under ``giap.steps`` so the events handler records the ``multiple_oric`` event (MVP-03c)
LOGGER = logging.getLogger("giap.steps.context_stats")

#: §9.3 ``statistic`` of every row of this analysis
STATISTIC = "oriter_slope"
#: §9.3 ``method`` of every row
METHOD = "ols_replichore"
#: §9.3 ``null_model``: the CI comes from a circular block bootstrap, not a null
NULL_MODEL = "bootstrap"
#: fit order = row order = draw order (lead decisions 2-3)
ARMS: tuple[str, ...] = ("pooled", "right", "left")
#: the two replichore labels, ``right`` = increasing coordinate from the oriC (decision 32)
REPLICHORES: tuple[str, str] = ("right", "left")
#: columns of :func:`gene_replichores` (the §9.3 ``genes.tsv`` ori-ter columns)
GENE_COLUMNS: tuple[str, ...] = ("locus_tag", "replichore", "replichore_position")
#: C15 feature type of an origin of replication call
ORIC_TYPE = "oric"
#: lead decision 31: several oriC calls are averaged only inside this share of the replicon
ORIC_SPREAD_FRAC = 0.05
#: §6.7 event listing the oriC calls when a replicon has more than one (decisions 4-5)
MULTIPLE_ORIC_EVENT = "context_stats.oriter.multiple_oric"
#: a fit needs at least this many TUs with a value
MIN_TUS = 3
#: cap of the l-hat search (ADR-0021 decision 1)
MAX_BLOCK_LEN = 20
#: level of the percentile CI
CI_LEVEL = 0.95
DESCRIPTIVE = "descriptive"
NO_ORIC = "no_oric"
NOT_ESTIMABLE = "not_estimable"
_FLOATS = ("value", "null_mean", "null_lo", "null_hi", "effect", "ci_low", "ci_high",
           "p_value", "q_value")  # fmt: skip
_NA = math.nan


# --------------------------------------------------------------------------- geometry


def _midpoint(start: float, end: float, length: float) -> float:
    """Midpoint of an (possibly unwrapped) interval, reduced modulo ``length``."""
    return float(((float(start) + float(end)) / 2.0) % length)


def _circ_dist(a: F64, b: F64, length: float) -> F64:
    """Circular distance between two coordinate arrays on a circle of ``length``."""
    d = np.abs(a - b) % length
    return np.asarray(np.minimum(d, length - d), dtype=np.float64)


def replichore(mid: npt.ArrayLike, oric: float, length: int) -> tuple[npt.NDArray[np.str_], F64]:
    """Replichore label and position of points on a circular replicon (spec "Replichore").

    Parameters
    ----------
    mid : array_like of float
        Point coordinates (gene or TU midpoints), 0-based; they are reduced modulo ``length``.
    oric : float
        oriC midpoint on the same circle.
    length : int
        Replicon length L.

    Returns
    -------
    labels : numpy.ndarray of str
        ``right`` where the clockwise distance d = (m - oriC) mod L is at most L / 2,
        ``left`` otherwise.
    positions : numpy.ndarray of float64
        d / (L / 2) on the right arm and (L - d) / (L / 2) on the left arm: 0 at the oriC and
        1 at the terminus, always in [0, 1].
    """
    m = np.asarray(mid, dtype=np.float64).ravel()
    size = float(length)
    half = size / 2.0
    d = (m - float(oric)) % size
    right = d <= half
    pos = np.where(right, d / half, (size - d) / half)
    labels = np.where(right, REPLICHORES[0], REPLICHORES[1])
    return np.asarray(labels, dtype=np.str_), np.asarray(pos, dtype=np.float64)


def _oric_midpoints(features: pd.DataFrame | None, contig: str, length: int) -> F64:
    """Midpoints of the ``oric`` calls of ``contig`` in C15 table order."""
    if features is None or not len(features):
        return np.zeros(0, dtype=np.float64)
    sel = features[
        (features["type"].astype(str) == ORIC_TYPE) & (features["contig"].astype(str) == contig)
    ]
    return np.asarray(
        [
            _midpoint(float(s), float(e), float(length))
            for s, e in zip(sel["start"], sel["end"], strict=True)
        ],
        dtype=np.float64,
    )


def _circular_mean(midpoints: F64, length: float) -> float:
    """Circular mean of ``midpoints`` on a circle of ``length`` (angles 2 pi m / L)."""
    ang = 2.0 * math.pi * midpoints / length
    mean = math.atan2(float(np.sin(ang).mean()), float(np.cos(ang).mean()))
    return float((mean / (2.0 * math.pi) * length) % length)


def oric_position(features: pd.DataFrame | None, contig: str, length: int) -> float | None:
    """oriC midpoint used for ``contig``, or None when the analysis has no origin.

    Parameters
    ----------
    features : pandas.DataFrame or None
        C15 ``features`` table (columns ``type``, ``contig``, ``start``, ``end``); None when the
        run has no ``features.parquet``, which behaves like no oriC call.
    contig : str
        Replicon name.
    length : int
        Replicon length L (the calls may be unwrapped past it, §6.5).

    Returns
    -------
    float or None
        The single call's midpoint, or the circular mean of several calls when every pair lies
        within ``ORIC_SPREAD_FRAC`` * L of each other (lead decision 31). None when there is no
        call or when the calls are farther apart than that. More than one call always logs
        ``context_stats.oriter.multiple_oric`` with their midpoints, on both branches.
    """
    mids = _oric_midpoints(features, contig, length)
    if not mids.size:
        return None
    if mids.size == 1:
        return float(mids[0])
    spread = float(_circ_dist(mids[:, None], mids[None, :], float(length)).max())
    accepted = spread <= ORIC_SPREAD_FRAC * float(length)
    mean = _circular_mean(mids, float(length)) if accepted else None
    LOGGER.info(
        "%s: %d oriC calls, spread %.0f bp (%s)",
        contig,
        mids.size,
        spread,
        "averaged" if accepted else "refused",
        extra={
            "event": MULTIPLE_ORIC_EVENT,
            "data": {
                "contig": contig,
                "midpoints": [float(m) for m in mids],
                "spread": spread,
                "accepted": accepted,
                "oric": mean,
            },
        },
    )
    return mean


# --------------------------------------------------------------------------- gene columns


def _analysed(replicons: pd.DataFrame) -> list[tuple[str, int]]:
    """(contig, length) of the replicons the analysis covers: circular chromosomes, C1 order."""
    if not len(replicons):
        return []
    sel = (replicons["kind"].astype(str) == "chromosome") & (
        replicons["topology"].astype(str) == "circular"
    )
    rows = replicons.loc[sel, ["contig", "length"]]
    return [(str(c), int(n)) for c, n in zip(rows["contig"], rows["length"], strict=True)]


def gene_replichores(
    genes: pd.DataFrame, replicons: pd.DataFrame, features: pd.DataFrame | None
) -> pd.DataFrame:
    """The §9.3 ``genes.tsv`` ori-ter columns, one row per gene in gene-table order.

    Parameters
    ----------
    genes : pandas.DataFrame
        Gene table with ``locus_tag``, ``contig``, ``start`` and ``end`` (0-based half-open;
        ``end`` may be unwrapped).
    replicons : pandas.DataFrame
        C1 replicons (``contig``, ``length``, ``topology``, ``kind``).
    features : pandas.DataFrame or None
        C15 features; the ``oric`` rows give the origin of each replicon.

    Returns
    -------
    pandas.DataFrame
        ``locus_tag``, ``replichore`` (``right`` / ``left``, None off the analysis) and
        ``replichore_position`` (NaN off the analysis). Both are NA on plasmids, linear and
        ``unknown``-kind contigs (lead decision 33) and on replicons without an oriC call.
    """
    origins = {
        contig: (oric_position(features, contig, length), length)
        for contig, length in _analysed(replicons)
    }
    arms: list[str | None] = []
    positions: list[float] = []
    for contig, start, end in zip(
        genes["contig"].astype(str), genes["start"], genes["end"], strict=True
    ):
        oric, length = origins.get(contig, (None, 0))
        if oric is None:
            arms.append(None)
            positions.append(_NA)
            continue
        label, pos = replichore([_midpoint(float(start), float(end), float(length))], oric, length)
        arms.append(str(label[0]))
        positions.append(float(pos[0]))
    return pd.DataFrame(
        {
            "locus_tag": genes["locus_tag"].astype(str).to_numpy(),
            "replichore": arms,
            "replichore_position": positions,
        }
    )


# --------------------------------------------------------------------------- TU points


def _tu_values(genes: pd.DataFrame, contig: str) -> dict[str, float]:
    """y per TU: the mean ``abundance_log2`` of its detected eligible genes with a value."""
    ix = replicon_index(genes, contig)
    out: dict[str, float] = {}
    if not len(ix):
        return out
    a = pd.to_numeric(ix["a"], errors="coerce").to_numpy(float)
    for tu, sub in pd.DataFrame({"tu": ix["tu_id"].astype(str), "a": a}).groupby("tu", sort=False):
        vals = sub["a"].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        if vals.size:
            out[str(tu)] = float(vals.mean())
    return out


def _tu_points(
    genes: pd.DataFrame, tus: pd.DataFrame, contig: str, length: int, oric: float
) -> tuple[F64, F64, npt.NDArray[np.str_]]:
    """(x, y, arm) of the TUs of one replicon with a value, in genome order."""
    values = _tu_values(genes, contig)
    sel = tus[tus["contig"].astype(str) == contig].sort_values(["start", "tu_id"], kind="mergesort")
    mids: list[float] = []
    ys: list[float] = []
    for tu, start, end in zip(sel["tu_id"].astype(str), sel["start"], sel["end"], strict=True):
        if tu in values:
            mids.append(_midpoint(float(start), float(end), float(length)))
            ys.append(values[tu])
    arm, x = replichore(mids, oric, length)
    return x, np.asarray(ys, dtype=np.float64), arm


def _arm_order(x: F64, arm: npt.NDArray[np.str_], which: str) -> npt.NDArray[np.int64]:
    """Row positions of one fit's sequence: pooled in genome order, an arm oriC -> terminus."""
    if which == "pooled":
        return np.arange(x.size, dtype=np.int64)
    sel = np.flatnonzero(np.asarray(arm) == which)
    return np.asarray(sel[np.argsort(x[sel], kind="mergesort")], dtype=np.int64)


# --------------------------------------------------------------------------- the fit


def _ols_slope(x: F64, y: F64) -> float:
    """OLS slope of ``y`` on ``x``; NaN below ``MIN_TUS`` points or without variance in x."""
    if x.size < MIN_TUS:
        return _NA
    sxx = float(((x - x.mean()) ** 2).sum())
    if sxx <= 0.0:
        return _NA
    return float(((x - x.mean()) * (y - y.mean())).sum() / sxx)


def _ols_resid(x: F64, y: F64) -> F64:
    """Residuals of the OLS fit of ``y`` on ``x`` (empty when the slope does not exist)."""
    b = _ols_slope(x, y)
    if math.isnan(b):
        return np.zeros(0, dtype=np.float64)
    a = float(y.mean() - b * x.mean())
    return np.asarray(y - (a + b * x), dtype=np.float64)


def _block_len(config: int | str, l_hat: int, n_tus: int) -> int:
    """Block length of one fit: an integer config clamped to T, else l' (ADR-0021 decision 1)."""
    if isinstance(config, int):
        return int(min(int(config), n_tus))
    return int(block_len_prime(l_hat, n_tus))


def _slope_ci(
    x: F64, y: F64, block_len: int, B: int, rng: np.random.Generator
) -> tuple[float, float, int]:
    """(ci_low, ci_high, valid replicates) of the circular block bootstrap of the slope."""
    draws = circular_block_bootstrap(np.arange(x.size), block_len, B, rng).astype(int)
    slopes = np.asarray([_ols_slope(x[d], y[d]) for d in draws], dtype=np.float64)
    valid = slopes[np.isfinite(slopes)]
    if valid.size < B / 2:
        return _NA, _NA, int(valid.size)
    lo, hi = percentile_ci(valid, CI_LEVEL)
    return lo, hi, int(valid.size)


# --------------------------------------------------------------------------- rows


def _row(contig: str, arm: str, flags: Sequence[str], **values: float | int) -> dict[str, Any]:
    """One ori-ter row: every column NA unless ``values`` overrides it (``n_perm`` 0)."""
    row: dict[str, Any] = {
        "contig": contig, "statistic": STATISTIC, "lag": _NA, "null_model": NULL_MODEL,
        "value": _NA, "null_mean": _NA, "null_lo": _NA, "null_hi": _NA,
        "stat": f"{STATISTIC}_{arm}", "effect": _NA, "ci_low": _NA, "ci_high": _NA,
        "p_value": _NA, "q_value": _NA, "method": METHOD, "n_perm": 0,
        "flags": tuple(sorted(flags)),
    }  # fmt: skip
    row.update(values)
    return row


def _replicon_rows(
    genes: pd.DataFrame,
    tus: pd.DataFrame,
    contig: str,
    length: int,
    oric: float | None,
    cfg: GiapConfig,
    rng: np.random.Generator,
) -> list[dict[str, Any]]:
    """Rows of one circular chromosome, in the ``ARMS`` order (lead decisions 1-3)."""
    if oric is None:
        return [_row(contig, "pooled", (DESCRIPTIVE, NO_ORIC))]

    x_all, y_all, arm_all = _tu_points(genes, tus, contig, length, oric)
    B = int(cfg.stats.bootstrap)
    config = cfg.context.islands.block_len
    pooled = _arm_order(x_all, arm_all, "pooled")
    pooled_slope = _ols_slope(x_all[pooled], y_all[pooled])
    l_hat = 0
    if not math.isnan(pooled_slope) and not isinstance(config, int):
        l_hat = int(
            choose_block_len(
                _ols_resid(x_all[pooled], y_all[pooled]),
                True,
                int(cfg.stats.permutations),
                rng,
                alpha=0.05,
                max_len=MAX_BLOCK_LEN,
            )
        )

    rows: list[dict[str, Any]] = []
    for arm in ARMS:
        idx = _arm_order(x_all, arm_all, arm)
        x, y = x_all[idx], y_all[idx]
        slope = _ols_slope(x, y)
        if math.isnan(slope):  # no bootstrap draw for a fit that does not exist (decision 2)
            rows.append(_row(contig, arm, (DESCRIPTIVE, NOT_ESTIMABLE)))
            continue
        lo, hi, n_valid = _slope_ci(x, y, _block_len(config, l_hat, int(x.size)), B, rng)
        flags = [DESCRIPTIVE] if n_valid >= B / 2 else [DESCRIPTIVE, NOT_ESTIMABLE]
        rows.append(
            _row(
                contig,
                arm,
                flags,
                value=slope,
                effect=slope,
                ci_low=lo,
                ci_high=hi,
                n_perm=n_valid,
            )
        )
    return rows


def _frame(rows: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    """``AUTOCORR_COLUMNS`` frame with the table's dtypes.

    ``lag`` is the nullable ``Int64`` of the §9.3 schema (always NA here, but an integer column
    with missing values, so concatenating these rows with the M5-05 ones in
    :func:`giap.context.tables.autocorrelation_table` keeps the lags integers).
    """
    df = pd.DataFrame(list(rows), columns=list(AUTOCORR_COLUMNS))
    df = df.astype({**dict.fromkeys(_FLOATS, "float64"), "n_perm": "int64", "lag": "Int64"})
    for col in ("contig", "statistic", "null_model", "stat", "method", "flags"):
        df[col] = df[col].astype(object)
    return df


def oriter_rows(
    genes: pd.DataFrame,
    tus: pd.DataFrame,
    replicons: pd.DataFrame,
    features: pd.DataFrame | None,
    cfg: GiapConfig,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """The ori-ter rows of the ``autocorrelation`` table (§9.3, §7.3.5(e)).

    Parameters
    ----------
    genes : pandas.DataFrame
        §9.3 gene table (``locus_tag``, ``contig``, ``start``, ``end``, ``partial``, ``tu_id``,
        ``detected``, ``p_detect``, ``abundance_log2``); eligibility follows
        ``giap.context.regions.replicon_index``.
    tus : pandas.DataFrame
        C5 transcription units (``tu_id``, ``contig``, ``start``, ``end``), 0-based half-open.
    replicons : pandas.DataFrame
        C1 replicons in table order; only circular chromosomes are analysed.
    features : pandas.DataFrame or None
        C15 features; its ``oric`` rows give each replicon's origin. None (no
        ``features.parquet``) gives a ``no_oric`` row per circular chromosome, not an error.
    cfg : giap.config.model.GiapConfig
        ``stats.bootstrap``, ``stats.permutations`` and ``context.islands.block_len``.
    rng : numpy.random.Generator
        Source of all randomness (``giap.rng.for_step("context_stats.oriter", seed)``), consumed
        replicon by replicon in C1 order: the l-hat permutations, then the bootstrap replicates
        of the pooled, right and left fits, skipping fits that do not exist.

    Returns
    -------
    pandas.DataFrame
        ``giap.context.autocorrelation.AUTOCORR_COLUMNS``: three rows (pooled, right, left) per
        circular chromosome with an oriC call and one ``no_oric`` pooled row per circular
        chromosome without one. ``value`` = ``effect`` = the OLS slope in log2 per replichore,
        ``ci_low`` / ``ci_high`` the 95 % block-bootstrap percentile interval, ``n_perm`` the
        valid replicates, ``p_value`` / ``q_value`` / ``lag`` and the null columns NA, and the
        flag ``descriptive`` on every row (plus ``no_oric`` or ``not_estimable``).
    """
    rows: list[dict[str, Any]] = []
    for contig, length in _analysed(replicons):
        oric = oric_position(features, contig, length)
        rows.extend(_replicon_rows(genes, tus, contig, length, oric, cfg, rng))
    return _frame(rows)
