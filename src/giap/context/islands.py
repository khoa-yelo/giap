"""Expression islands along a replicon: TU scores, detrending, scan and block-bootstrap FWER.

GIAP_SPEC.md §7.3.5(d) steps 1-10, §9.3 ``tables/islands.tsv``, §6.5 IDs; issue M5-06
(docs/specs/M5-06-islands.md, lead decisions 21-27 and 40-51), ADR-0021 with Amendment 1;
docs/methods/islands.md.

- **Scores** (step 1). Within a replicon the eligible genes (not partial, ``p_detect`` present)
  are ranked by log2 abundance with the undetected ones left-censored at the lowest rank and
  converted to Blom normal scores; a detected gene without an abundance is not scored
  (decision 22).
- **TU series** (steps 2-3). z_t is the mean score of the scored genes of TU t, TUs follow the
  genome order of their first scored gene, and the series is detrended by LOESS over TU
  ordinals (span ``context.islands.detrend_span_frac``, circular padding on circular replicons).
- **Null and family** (steps 4-6). Each scanned replicon is resampled by a circular block
  bootstrap with l' = min(20, 2 l-hat + 1, T) (ADR-0021 decision 1), and M_b is the maximum of
  |S_w| over the windows of **every** scanned replicon: the family is genome-wide (decisions 21
  and 46), so one FWER covers the run. A replicon whose LOESS span gives fewer than
  ``MIN_SPAN_POINTS`` points, or that has fewer than 2 ``min_tus`` scored TUs, is not scanned:
  it carries the run flag ``few_values`` and enters neither maximum (decisions 40-41).
- **Islands** (steps 7-9). Significant same-direction windows are merged, and each island is
  reported by its **best window**, never by the union of its cluster (decision 42); the effect
  CI resamples the window's TUs (conditional on the selected window, ADR-0021).
- **Validity domain.** ADR-0021 Amendment 1: the FWER is claimed for the dependence regime the
  calibration covers, not for a replicon size, and no row is flagged ``uncalibrated``
  (decision 48). ``IslandScan`` carries l-hat, l' and the lag-1 Moran's I of the detrended
  scores per replicon so the §9.4 run report can show where the data sit.
- **The low direction is a structural zero at normal detection rates** (decision 50). Every
  undetected gene shares the censored Blom score c = Φ⁻¹(q / 2), so a fully silenced window of w
  TUs cannot reach below -√w·|c| — the ceiling of :func:`low_direction_ceiling`, about 3.4 at
  q = 0.45 and w = 20, against a genome-wide critical value of 4.1-4.5. Left censoring therefore
  caps low islands: this module runs a two-sided procedure, but its low half only resolves
  islands when the censored share is small (below about a quarter undetected).

Pure functions; randomness comes only from the generator passed in
(``giap.rng.for_step("context_stats.islands", seed)``), drawn as: the block-length permutations
of every scanned replicon in C1 order, then the bootstrap draws in C1 order, then the effect CIs
in island order.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.special import ndtri

from giap.config.model import GiapConfig
from giap.stats.autocorr import block_len_prime, choose_block_len, morans_i
from giap.stats.bootstrap import circular_block_bootstrap, tu_block_bootstrap
from giap.stats.scan import Cluster, ScanResult, fwer_p, merge_windows, scan_max, scan_null_max
from giap.stats.scores import blom_scores
from giap.stats.smooth import loess_detrend, loess_span_points

__all__ = [
    "CI_FEW_UNITS",
    "FEW_VALUES",
    "GENE_MAP_COLUMNS",
    "ISLAND_COLUMNS",
    "METHOD",
    "MIN_SPAN_POINTS",
    "NULL_MODEL",
    "STAT",
    "IslandScan",
    "TuSeries",
    "block_length",
    "gene_island_map",
    "gene_scores",
    "islands_frame",
    "low_direction_ceiling",
    "scan_replicons",
    "scan_with_nulls",
    "tu_series",
]

F64 = npt.NDArray[np.float64]
I64 = npt.NDArray[np.int64]

#: §9.3 ``islands.tsv`` columns, the §9.2 test columns and ``flags`` (decision 44), in table order
ISLAND_COLUMNS: tuple[str, ...] = (
    "island_id", "direction", "contig", "start", "end", "tu_first", "tu_last", "n_tus",
    "n_genes", "genes", "mean_abundance_log2", "detected_frac", "scan_s", "block_len",
    "stat", "effect", "ci_low", "ci_high", "p_value", "q_value", "method", "null_model",
    "n_perm", "flags",
)  # fmt: skip
#: columns of the gene map of step 9
GENE_MAP_COLUMNS: tuple[str, ...] = ("locus_tag", "island_id")
METHOD = "scan_max_stat"
NULL_MODEL = "circular_block_bootstrap"
#: §9.2 test column ``stat``: the *name* of the row's statistic (M5-09 lead decision 6); its
#: value is in ``scan_s``, which every row already carries.
STAT = "scan_max"
#: decisions 40-41: below this LOESS span in points the detrend is degenerate, so no scan
MIN_SPAN_POINTS = 4
#: run flag of a replicon that is not scanned (decision 40)
FEW_VALUES = "few_values"
#: M4-16 flag of a CI drawn from fewer than ``CI_MIN_UNITS`` resampling units
CI_FEW_UNITS = "ci_few_units"
CI_MIN_UNITS = 10
#: §7.3.5(d) step 5 keeps the literal 0.05 for the block-length rule (decision 24)
BLOCK_LEN_ALPHA = 0.05
#: ADR-0021: l-hat is searched over lags 1 ... 20
BLOCK_LEN_MAX_LAG = 20
_FLOATS = ("mean_abundance_log2", "detected_frac", "scan_s", "effect", "ci_low", "ci_high",
           "p_value", "q_value")  # fmt: skip
_INTS = ("start", "end", "n_tus", "n_genes", "block_len", "n_perm")
_OBJECTS = ("island_id", "direction", "contig", "tu_first", "tu_last", "genes", "stat",
            "method", "null_model", "flags")  # fmt: skip
_NA = math.nan


@dataclass(frozen=True, slots=True)
class TuSeries:
    """The scored TUs of one replicon in genome order (§7.3.5(d) steps 2-3, decision 43).

    Attributes
    ----------
    contig : str
        Replicon name.
    circular : bool
        Whether the replicon is circular (windows and the LOESS padding wrap).
    length : int
        Replicon length in bp (C1), used to unwrap the coordinates of a wrapping island.
    tu_ids : tuple of str
        Scored TU IDs in genome order.
    genes : tuple of tuple of str
        Scored genes of each TU, in genome order.
    start, end : numpy.ndarray of int64
        0-based half-open span of each TU (C5).
    score : numpy.ndarray of float64
        z_t, the mean Blom score of the TU's scored genes.
    detrended : numpy.ndarray of float64
        z_t minus its LOESS trend along TU order.
    """

    contig: str
    circular: bool
    length: int
    tu_ids: tuple[str, ...]
    genes: tuple[tuple[str, ...], ...]
    start: I64
    end: I64
    score: F64
    detrended: F64


@dataclass(frozen=True, slots=True)
class IslandScan:
    """The scan of every replicon against one genome-wide null maximum (steps 4-6).

    Attributes
    ----------
    series : tuple of TuSeries
        The replicons in C1 order, scanned or not.
    scanned : tuple of bool
        Whether each replicon was scanned (decisions 40-41).
    l_hat : tuple of int
        Per replicon, the Moran's I block length l-hat; 0 when the replicon is not scanned or
        ``context.islands.block_len`` fixes the block length (§9.4 run report, decision 48).
    block_len : tuple of int
        Per replicon, the bootstrap block length l'; 0 when the replicon is not scanned.
    morans_i_lag1 : tuple of float
        Per replicon, the lag-1 Moran's I of the detrended TU scores; NaN when not scanned.
    windows : tuple of giap.stats.scan.ScanResult
        Per replicon, the observed window table (empty when not scanned).
    p_value : tuple of numpy.ndarray of float64
        Per replicon, the FWER p of each window against ``null_max``.
    clusters : tuple of tuple of giap.stats.scan.Cluster
        Per replicon, the merged significant windows.
    null_max : numpy.ndarray of float64
        M_b, the maximum of |S_w| over the windows of every scanned replicon; empty when no
        replicon is scanned.
    n_perm : int
        Number of valid (non-NaN) M_b.
    flags : tuple of tuple of str
        Per replicon run flags (``few_values``).
    """

    series: tuple[TuSeries, ...]
    scanned: tuple[bool, ...]
    l_hat: tuple[int, ...]
    block_len: tuple[int, ...]
    morans_i_lag1: tuple[float, ...]
    windows: tuple[ScanResult, ...]
    p_value: tuple[F64, ...]
    clusters: tuple[tuple[Cluster, ...], ...]
    null_max: F64
    n_perm: int
    flags: tuple[tuple[str, ...], ...]


# ----------------------------------------------------------------------------- scores


def _scored(ix: pd.DataFrame) -> tuple[list[str], list[str], F64]:
    """Scored genes of a replicon index: their locus tags, TU IDs and Blom scores (step 1)."""
    if not len(ix):
        return [], [], np.zeros(0, dtype=np.float64)
    detected = np.asarray(ix["d"].to_numpy(), dtype=bool)
    a = np.asarray(pd.to_numeric(ix["a"], errors="coerce").to_numpy(float), dtype=np.float64)
    keep = ~detected | np.isfinite(a)  # decision 22: detected without a value is not scored
    det = detected[keep]
    z = blom_scores(np.where(det, a[keep], 0.0), ~det)
    tags = [str(t) for t in ix["locus_tag"].to_numpy()[keep]]
    tu_ids = [str(t) for t in ix["tu_id"].to_numpy()[keep]]
    return tags, tu_ids, z


def gene_scores(ix: pd.DataFrame) -> pd.Series:
    """Left-censored Blom scores of the scored genes of one replicon (§7.3.5(d) step 1).

    Parameters
    ----------
    ix : pandas.DataFrame
        Replicon index of ``giap.context.regions.replicon_index`` (eligible genes in genome
        order, columns ``locus_tag``, ``tu_id``, ``d``, ``a``).

    Returns
    -------
    pandas.Series
        Blom score per scored gene, indexed by ``locus_tag`` in genome order. Undetected genes
        are left-censored (tied at the lowest rank); a detected gene without a finite abundance
        is not scored and is absent from the result (decision 22).
    """
    tags, _, z = _scored(ix)
    return pd.Series(z, index=pd.Index(tags, name="locus_tag"), dtype="float64")


def _tu_spans(tus: pd.DataFrame, contig: str, tu_ids: Sequence[str]) -> tuple[I64, I64]:
    """0-based half-open (start, end) of each TU of ``contig`` from the C5 ``tus`` frame."""
    sub = tus[tus["contig"].astype(str) == contig]
    spans = {
        str(t): (int(s), int(e))
        for t, s, e in zip(sub["tu_id"], sub["start"], sub["end"], strict=True)
    }
    missing = [t for t in tu_ids if t not in spans]
    if missing:
        raise KeyError(f"TUs {missing[:3]} of contig {contig!r} are absent from the tus frame")
    start = np.array([spans[t][0] for t in tu_ids], dtype=np.int64)
    end = np.array([spans[t][1] for t in tu_ids], dtype=np.int64)
    return start, end


def tu_series(
    ix: pd.DataFrame,
    tus: pd.DataFrame,
    *,
    contig: str,
    circular: bool,
    length: int,
    span_frac: float,
) -> TuSeries:
    """TU scores of one replicon and their LOESS residuals (§7.3.5(d) steps 2-3).

    Parameters
    ----------
    ix : pandas.DataFrame
        Replicon index of ``giap.context.regions.replicon_index``.
    tus : pandas.DataFrame
        C5 ``tus`` frame (``tu_id``, ``contig``, ``start``, ``end``), 0-based half-open.
    contig : str
        Replicon name.
    circular : bool
        Whether the replicon is circular (the detrend pads circularly with period T).
    length : int
        Replicon length in bp (C1).
    span_frac : float
        LOESS span as a fraction of the scored TUs (``context.islands.detrend_span_frac``).

    Returns
    -------
    TuSeries
        TUs in genome order by C5 (``start``, ``tu_id``); a TU with no scored gene is dropped.

    Raises
    ------
    KeyError
        If a scored TU is absent from ``tus``.
    """
    tags, tu_of, z = _scored(ix)
    order: list[str] = []
    members: dict[str, list[str]] = {}
    values: dict[str, list[float]] = {}
    for tag, tu, score in zip(tags, tu_of, z.tolist(), strict=True):
        if tu not in members:
            order.append(tu)
            members[tu] = []
            values[tu] = []
        members[tu].append(tag)
        values[tu].append(float(score))
    # Genome order of the TUs by C5 (start, tu_id) (spec step 2): a TU crossing the origin
    # (start near length, end > length, §6.5) sorts last, so the window wrap rule holds.
    start, end = _tu_spans(tus, contig, order)
    rank = sorted(range(len(order)), key=lambda i: (int(start[i]), order[i]))
    order = [order[i] for i in rank]
    start, end = start[rank], end[rank]
    mean = np.array([float(np.mean(values[t])) for t in order], dtype=np.float64)
    n = mean.size
    if n:
        fit = loess_detrend(
            mean,
            np.arange(n, dtype=float),
            span_frac,
            circular=circular,
            period=float(n) if circular else None,
        )
        resid = np.asarray(fit.resid, dtype=np.float64)
    else:
        resid = np.zeros(0, dtype=np.float64)
    return TuSeries(
        contig=contig,
        circular=circular,
        length=int(length),
        tu_ids=tuple(order),
        genes=tuple(tuple(members[t]) for t in order),
        start=start,
        end=end,
        score=mean,
        detrended=resid,
    )


# ----------------------------------------------------------------------------- block length


def _block_lengths(ts: TuSeries, cfg: GiapConfig, rng: np.random.Generator) -> tuple[int, int]:
    """(l-hat, l') of one series; l-hat is 0 when the config fixes the block length."""
    t = int(ts.detrended.size)
    if t < 1:
        return 0, 0
    fixed = cfg.context.islands.block_len
    if isinstance(fixed, int):  # decision 25: an integer overrides l', bounded to T
        if fixed < 1:
            raise ValueError(f"context.islands.block_len must be >= 1, got {fixed}")
        return 0, min(int(fixed), t)
    l_hat = choose_block_len(
        ts.detrended,
        ts.circular,
        int(cfg.stats.permutations),
        rng,
        alpha=BLOCK_LEN_ALPHA,
        max_len=BLOCK_LEN_MAX_LAG,
    )
    return int(l_hat), block_len_prime(int(l_hat), t)


def block_length(ts: TuSeries, cfg: GiapConfig, rng: np.random.Generator) -> int:
    """Bootstrap block length l' of one replicon (§7.3.5(d) step 5, ADR-0021 decision 1).

    Parameters
    ----------
    ts : TuSeries
        The replicon's detrended TU scores.
    cfg : giap.config.model.GiapConfig
        ``context.islands.block_len`` (an integer overrides l', bounded to T) and
        ``stats.permutations`` (the gene-label permutations of the l-hat test).
    rng : numpy.random.Generator
        Source of the l-hat permutations; untouched when the config fixes the block length.

    Returns
    -------
    int
        l' = min(20, 2 l-hat + 1, T) with l-hat from
        ``giap.stats.autocorr.choose_block_len`` at the literal alpha 0.05 (decision 24);
        0 for a series with no TU.

    Raises
    ------
    ValueError
        If ``context.islands.block_len`` is an integer below 1.
    """
    return _block_lengths(ts, cfg, rng)[1]


def low_direction_ceiling(q: float, w: int) -> float:
    """Most negative S a fully undetected window of ``w`` TUs can reach (decision 50).

    Parameters
    ----------
    q : float
        Share of the replicon's scored genes that are undetected (left-censored).
    w : int
        Window length in TUs.

    Returns
    -------
    float
        √w·|Φ⁻¹(q / 2)|, the magnitude of the ceiling, before the LOESS attenuation; NaN when q
        is outside (0, 1). Every undetected gene shares the censored Blom score Φ⁻¹(q / 2), so
        no low window can exceed this, which at normal detection rates is below the genome-wide
        critical value: the low direction is then a structural zero, not merely low-powered.
    """
    if not 0.0 < float(q) < 1.0:
        return math.nan
    return float(math.sqrt(w) * abs(float(ndtri(float(q) / 2.0))))


# ----------------------------------------------------------------------------- the scan


def _is_scanned(ts: TuSeries, cfg: GiapConfig) -> bool:
    """Decisions 40-41: T >= 2 ``min_tus`` and a LOESS span of at least ``MIN_SPAN_POINTS``."""
    t = len(ts.tu_ids)
    if t < 2 * int(cfg.context.islands.min_tus):
        return False
    return loess_span_points(t, float(cfg.context.islands.detrend_span_frac)) >= MIN_SPAN_POINTS


def scan_replicons(
    series: Sequence[TuSeries], cfg: GiapConfig, rng: np.random.Generator
) -> IslandScan:
    """The scan (``scan_with_nulls`` without the per-replicon nulls); see there."""
    return scan_with_nulls(series, cfg, rng)[0]


def scan_with_nulls(
    series: Sequence[TuSeries], cfg: GiapConfig, rng: np.random.Generator
) -> tuple[IslandScan, tuple[F64, ...]]:
    """Scan every replicon against one genome-wide null maximum (§7.3.5(d) steps 4-6).

    Parameters
    ----------
    series : sequence of TuSeries
        The replicons in C1 order.
    cfg : giap.config.model.GiapConfig
        ``context.islands`` (``min_tus``, ``max_tus``, ``detrend_span_frac``, ``block_len``),
        ``stats.bootstrap``, ``stats.permutations`` and ``stats.alpha.context``.
    rng : numpy.random.Generator
        Source of all randomness: the block-length permutations of every scanned replicon in C1
        order, then the bootstrap draws in C1 order.

    Returns
    -------
    (IslandScan, per-replicon null maxima)
        Second item (M9-17 V18): per replicon its own max |S_w| per replicate (empty when not
        scanned); their ``fmax`` is ``null_max``. One ``null_max`` for the whole genome
        (decisions 21 and 46): M_b is the maximum of
        |S_w| over the windows of every scanned replicon, so a replicon that is not scanned
        contributes nothing and draws nothing. Windows are significant when ``p <= alpha``, so a
        NaN p never counts (decision 42).
    """
    isl = cfg.context.islands
    wmin, wmax = int(isl.min_tus), int(isl.max_tus)
    alpha = float(cfg.stats.alpha.context)
    n_boot = int(cfg.stats.bootstrap)
    scanned = tuple(_is_scanned(ts, cfg) for ts in series)

    l_hat: list[int] = []
    block_len: list[int] = []
    for ts, ok in zip(series, scanned, strict=True):
        hat, ell = _block_lengths(ts, cfg, rng) if ok else (0, 0)
        l_hat.append(hat)
        block_len.append(ell)

    null = np.full(n_boot, np.nan) if any(scanned) else np.zeros(0, dtype=np.float64)
    per: list[F64] = []
    for ts, ok, ell in zip(series, scanned, block_len, strict=True):
        if not ok:
            per.append(np.zeros(0, dtype=np.float64))
            continue
        reps = circular_block_bootstrap(ts.detrended, ell, n_boot, rng)
        per.append(np.asarray(scan_null_max(reps, wmin, wmax, ts.circular), dtype=np.float64))
        np.fmax(null, per[-1], out=null)

    windows: list[ScanResult] = []
    p_values: list[F64] = []
    clusters: list[tuple[Cluster, ...]] = []
    for ts, ok in zip(series, scanned, strict=True):
        if not ok:
            windows.append(scan_max(np.zeros(0), wmin, wmax, ts.circular))
            p_values.append(np.zeros(0, dtype=np.float64))
            clusters.append(())
            continue
        res = scan_max(ts.detrended, wmin, wmax, ts.circular)
        p = fwer_p(np.abs(res.s), null)  # decision 42: fwer_p raises on a negative |S|
        with np.errstate(invalid="ignore"):
            significant = p <= alpha  # a NaN p is never significant
        windows.append(res)
        p_values.append(p)
        clusters.append(tuple(merge_windows(res, significant, ts.detrended.size, ts.circular)))

    return IslandScan(
        series=tuple(series),
        scanned=scanned,
        l_hat=tuple(l_hat),
        block_len=tuple(block_len),
        morans_i_lag1=tuple(
            morans_i(ts.detrended, 1, ts.circular) if ok else math.nan
            for ts, ok in zip(series, scanned, strict=True)
        ),
        windows=tuple(windows),
        p_value=tuple(p_values),
        clusters=tuple(clusters),
        null_max=null,
        n_perm=int(np.count_nonzero(~np.isnan(null))),
        flags=tuple((FEW_VALUES,) if not ok else () for ok in scanned),
    ), tuple(per)


# ----------------------------------------------------------------------------- the table


def _eligible(genes: pd.DataFrame, contig: str) -> pd.DataFrame:
    """Eligible genes of ``contig`` in genome order (not partial, ``p_detect`` present)."""
    g = genes[genes["contig"].astype(str) == contig]
    p = pd.to_numeric(g["p_detect"], errors="coerce")
    g = g[~g["partial"].astype("boolean").fillna(False).to_numpy(bool) & p.notna().to_numpy()]
    return g.sort_values(["start", "locus_tag"], kind="mergesort").reset_index(drop=True)


def _window_positions(ts: TuSeries, start: int, length: int) -> list[int]:
    """TU positions covered by a window, wrapping on a circular replicon."""
    n = len(ts.tu_ids)
    return [(start + k) % n if ts.circular else start + k for k in range(length)]


def _island_row(
    ts: TuSeries, res: ScanResult, best: int, direction: str, block_len: int, p: float,
    n_perm: int, eligible: pd.DataFrame,
) -> dict[str, Any]:  # fmt: skip
    """One island from its best window (§7.3.5(d) step 7, decision 42)."""
    start_idx, width = int(res.start[best]), int(res.length[best])
    pos = _window_positions(ts, start_idx, width)
    tu_ids = [ts.tu_ids[j] for j in pos]
    eg = eligible[eligible["tu_id"].astype(str).isin(set(tu_ids)).to_numpy()]
    det = eg["detected"].astype("boolean").fillna(False).to_numpy(bool)
    a = pd.to_numeric(eg["abundance_log2"], errors="coerce").to_numpy(float)
    have = det & np.isfinite(a)
    wraps = ts.circular and start_idx + width > len(ts.tu_ids)
    return {
        "direction": direction,
        "contig": ts.contig,
        "start": int(ts.start[pos[0]]),
        "end": int(ts.end[pos[-1]]) + (ts.length if wraps else 0),
        "tu_first": tu_ids[0],
        "tu_last": tu_ids[-1],
        "n_tus": len(pos),
        "n_genes": len(eg),
        "genes": [str(t) for t in eg["locus_tag"]],
        "mean_abundance_log2": float(a[have].mean()) if have.any() else _NA,
        "detected_frac": float(det.mean()) if len(eg) else _NA,
        "scan_s": float(res.s[best]),
        "block_len": int(block_len),
        "stat": STAT,
        "effect": float(np.mean(ts.detrended[pos])),
        "p_value": float(p),
        "q_value": _NA,
        "method": METHOD,
        "null_model": NULL_MODEL,
        "n_perm": int(n_perm),
        "_positions": pos,
    }


def _frame(rows: Sequence[dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(list(rows), columns=list(ISLAND_COLUMNS))
    df = df.astype({**dict.fromkeys(_FLOATS, "float64"), **dict.fromkeys(_INTS, "int64")})
    for c in _OBJECTS:
        df[c] = df[c].astype(object)
    return df


def islands_frame(
    scan: IslandScan, genes: pd.DataFrame, cfg: GiapConfig, rng: np.random.Generator
) -> pd.DataFrame:
    """The ``islands`` table: one row per merged cluster, from its best window (steps 7-8).

    Parameters
    ----------
    scan : IslandScan
        Result of :func:`scan_replicons`.
    genes : pandas.DataFrame
        §9.3 ``genes`` columns (``contig``, ``start``, ``locus_tag``, ``partial``, ``tu_id``,
        ``detected``, ``p_detect``, ``abundance_log2``); the island's gene list holds every
        *eligible* gene of the window's TUs, not only the scored ones.
    cfg : giap.config.model.GiapConfig
        ``stats.bootstrap`` (the effect CI replicates).
    rng : numpy.random.Generator
        Source of the effect-CI draws, consumed in island order.

    Returns
    -------
    pandas.DataFrame
        ``ISLAND_COLUMNS``, ordered by contig (C1 order), then ``start``, then direction (high
        first), with IDs ``ISL_<n>`` in that order. ``stat`` is the §9.2 statistic *name*
        :data:`STAT` and its value is ``scan_s``; ``q_value`` is NA (§7.3.5(d) step 7) and
        ``flags`` carries ``ci_few_units`` when the best window has fewer than
        ``CI_MIN_UNITS`` TUs; ``uncalibrated`` is never emitted (decision 48).
    """
    rows: list[dict[str, Any]] = []
    for r, (ts, res, p, cl) in enumerate(
        zip(scan.series, scan.windows, scan.p_value, scan.clusters, strict=True)
    ):
        if not cl:
            continue
        eligible = _eligible(genes, ts.contig)
        for cluster in cl:
            row = _island_row(
                ts, res, cluster.best, cluster.direction, scan.block_len[r],
                float(p[cluster.best]), scan.n_perm, eligible,
            )  # fmt: skip
            row["_replicon"] = r
            rows.append(row)
    rows.sort(key=lambda r: (r["_replicon"], r["start"], r["direction"] != "high"))

    n_boot = int(cfg.stats.bootstrap)
    for n, row in enumerate(rows, start=1):
        ts = scan.series[row.pop("_replicon")]
        pos = row.pop("_positions")
        row["island_id"] = f"ISL_{n}"
        units = pd.DataFrame(
            {"tu_id": [ts.tu_ids[j] for j in pos], "z": [float(ts.detrended[j]) for j in pos]}
        )
        ci = tu_block_bootstrap(units, "tu_id", _mean_z, n_boot, rng)
        row["ci_low"], row["ci_high"] = ci.ci_low, ci.ci_high
        row["flags"] = (CI_FEW_UNITS,) if int(row["n_tus"]) < CI_MIN_UNITS else ()
    return _frame(rows)


def _mean_z(df: pd.DataFrame) -> float:
    """Effect statistic of the CI: the mean detrended score of the resampled TUs."""
    return float(df["z"].mean())


def _island_number(island_id: str) -> int:
    """``n`` of an ``ISL_<n>`` ID; a non-conforming ID sorts last."""
    tail = str(island_id).rsplit("_", 1)[-1]
    return int(tail) if tail.isdigit() else np.iinfo(np.int64).max


def gene_island_map(islands: pd.DataFrame) -> pd.DataFrame:
    """``locus_tag`` -> ``island_id`` for every gene of every island (§7.3.5(d) step 9).

    Parameters
    ----------
    islands : pandas.DataFrame
        The :func:`islands_frame` table (columns ``island_id``, ``p_value``, ``scan_s``,
        ``genes``).

    Returns
    -------
    pandas.DataFrame
        ``GENE_MAP_COLUMNS``, one row per gene in order of first appearance in ``islands``. A
        gene in two islands (only a high and a low island can overlap) takes the one with the
        smaller p, then the larger |S|, then the smaller ISL number.
    """
    best: dict[str, tuple[tuple[float, float, int], str]] = {}
    for row in islands[["island_id", "p_value", "scan_s", "genes"]].itertuples(index=False):
        island_id = str(row.island_id)
        p = float(row.p_value)
        key = (math.inf if math.isnan(p) else p, -abs(float(row.scan_s)),
               _island_number(island_id))  # fmt: skip
        for tag in (str(g) for g in row.genes):
            if tag not in best or key < best[tag][0]:
                best[tag] = (key, island_id)
    out = pd.DataFrame(
        {"locus_tag": list(best), "island_id": [v[1] for v in best.values()]},
        columns=list(GENE_MAP_COLUMNS),
    )
    return out.astype({"locus_tag": object, "island_id": object})
