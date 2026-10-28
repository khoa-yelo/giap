"""Region tests: O - E and Hodges-Lehmann shift under circular gene-index shifts (issue M4-17).

GIPPA_SPEC.md §7.3.5(b), §8.1, Appendix B.4; docs/specs/M4-17-region-tests.md. The index of a
replicon is its eligible genes (not partial, with a ``p_detect``) in genome order, 0 ... N - 1; a
region's index set I is moved to (I + k) mod N for every valid shift k (``circular_shifts``,
``choose_shifts``), both statistics are recomputed at each shift and the two-sided centred
``perm_p`` is reported with a TU-block bootstrap CI. q-values: BH within (region_type, stat).
ADR-0018: the O - E test statistic is studentised, sum(d - p) / sqrt(sum p (1 - p)) (the reported
effect stays O - E), and HL runs on log2 abundance residualised per replicon on the detectability
covariates (``residualise``; effect, CI and p on the residual scale). ADR-0018 amendment 1: the HL
p-value comes from T = HL / sqrt(V) at every shift; amendment 2: V is label-free, built from the
replicon's lag autocovariance of the residuals (``lag_autocovariance``, ``hl_variance``,
``studentised_hl_null``); V <= 0 gives an NA p flagged ``hl_variance_nonpositive``.
Pure functions: randomness only from the generator passed in (``for_step("context_stats", seed)``).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import astuple, dataclass, fields
from typing import Literal

import numpy as np
import numpy.typing as npt
import pandas as pd

from gippa.config.model import GippaConfig
from gippa.stats.bootstrap import tu_block_bootstrap
from gippa.stats.multitest import bh
from gippa.stats.permutation import choose_shifts, circular_shifts, perm_p
from gippa.stats.robust import hodges_lehmann

__all__ = ["COVARIATES", "LAG_MAX", "REGION_TYPES", "RegionNull", "RegionTest", "hl_null",
           "hl_variance", "lag_autocovariance", "o_minus_e_null", "region_tests",
           "replicon_index", "residualise", "studentised_hl_null", "studentised_null",
           "test_all_regions", "test_region"]  # fmt: skip

#: region types tested here (§7.3.5(b)); ``plasmid`` regions belong to M4-18
REGION_TYPES = ("prophage", "bgc", "crispr", "amr")
METHOD = {"o_minus_e": "circular_shift_o_minus_e", "hl_shift": "circular_shift_hodges_lehmann"}
FEW_SHIFTS_K = 19  # K < 19: the smallest attainable p, 1/(1 + K), exceeds 0.05
SMALL_REGION = 3
#: detectability covariates of genes.tsv (§9.3) residualised out of HL (ADR-0018, decision 10)
COVARIATES = ("length_aa", "n_observable_peptides", "gravy")
_ENUM = 64  # HL null: enumerate a shift's candidate differences once at most this many remain
_CHUNK = 4_000_000  # HL null: shifts per chunk such that K * m * m stays below this
#: decision 12: residual variance below (V_RTOL * data magnitude)^2 is fit rounding: V counts as 0
V_RTOL = 1e-7
#: ADR-0018 amendment 2: largest gene-index lag of the residual autocovariance
LAG_MAX = 30

F64 = npt.NDArray[np.float64]
I64 = npt.NDArray[np.int64]
B8 = npt.NDArray[np.bool_]


@dataclass(frozen=True, slots=True)
class RegionTest:
    """One (region, contig, statistic) row of ``regions_stats.tsv`` before BH."""

    region_id: str
    region_type: str
    contig: str
    stat: Literal["o_minus_e", "hl_shift"]
    observed: float
    expected: float | None
    effect: float
    ci_low: float
    ci_high: float
    p_value: float
    n_perm: int
    n_genes: int
    n_eligible: int
    flags: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class RegionNull:
    """The null a region test compared its statistic with (M9-17 V18): the studentised value
    at each used shift (NaN where undefined) and the observed studentised statistic."""

    region_id: str
    contig: str
    stat: Literal["o_minus_e", "hl_shift"]
    shifts: I64
    null: F64
    observed: float


def replicon_index(genes: pd.DataFrame, contig: str) -> pd.DataFrame:
    """Eligible genes of ``contig`` in genome order: ``idx`` (0 ... N - 1), ``locus_tag``,
    ``tu_id``, ``d`` (detected, 0/1), ``p`` (p_detect), ``a`` (log2 abundance of detected
    genes, NaN otherwise) and the ``COVARIATES`` (NaN when a column is absent). Eligible = not
    partial and ``p_detect`` present."""
    g = genes[genes["contig"].astype(str) == contig]
    p = pd.to_numeric(g["p_detect"], errors="coerce").astype(float)
    g = g[~g["partial"].astype("boolean").fillna(False).to_numpy(bool) & p.notna().to_numpy()]
    g = g.sort_values(["start", "locus_tag"], kind="mergesort")
    d = g["detected"].astype("boolean").fillna(False).to_numpy(bool)
    a = pd.to_numeric(g["abundance_log2"], errors="coerce").to_numpy(float)
    nan = np.full(len(g), np.nan)
    cov = {c: pd.to_numeric(g[c], errors="coerce").to_numpy(float) if c in g else nan
           for c in COVARIATES}  # fmt: skip
    return pd.DataFrame(
        {
            "idx": np.arange(len(g), dtype=np.int64),
            "locus_tag": g["locus_tag"].astype(str).to_numpy(),
            "tu_id": g["tu_id"].to_numpy(),
            "d": d.astype(np.int64),
            "p": pd.to_numeric(g["p_detect"]).to_numpy(float),
            "a": np.where(d, a, np.nan),
            **cov,
        }
    )


def residualise(a: npt.ArrayLike, det: npt.ArrayLike, cov: pd.DataFrame) -> F64:
    """ADR-0018 decision 1 / M4-17 decision 10: residuals of the OLS (with intercept) of ``a`` on
    log(length_aa), log1p(n_observable_peptides) and gravy (columns of ``cov``), fitted on the
    genes detected with a finite ``a``; NaN elsewhere. A non-finite transformed covariate takes
    that column's median over the fitted genes (0 if none is finite); rank-deficient designs use
    the minimum-norm least-squares solution, so constant covariates only centre ``a``."""
    av = np.asarray(a, dtype=np.float64).ravel()
    ok = np.asarray(det, dtype=bool).ravel() & np.isfinite(av)
    out = np.full(av.size, np.nan)
    if not ok.any():
        return out
    with np.errstate(divide="ignore", invalid="ignore"):
        cols = (
            np.log(cov["length_aa"].to_numpy(float)),
            np.log1p(cov["n_observable_peptides"].to_numpy(float)),
            cov["gravy"].to_numpy(float),
        )
    design = [np.ones(int(ok.sum()))]
    for c in cols:
        x = c[ok].copy()
        bad = ~np.isfinite(x)
        if bad.any():
            x[bad] = float(np.median(x[~bad])) if (~bad).any() else 0.0
        design.append(x)
    X = np.column_stack(design)
    beta = np.linalg.lstsq(X, av[ok], rcond=None)[0]
    out[ok] = av[ok] - X @ beta
    return out


def _sets(idx: npt.ArrayLike, shifts: npt.ArrayLike, n: int) -> I64:
    i = np.asarray(idx, dtype=np.int64).ravel()
    k = np.asarray(shifts, dtype=np.int64).ravel()
    return np.asarray((i[None, :] + k[:, None]) % n, dtype=np.int64)


def o_minus_e_null(r: npt.ArrayLike, idx: npt.ArrayLike, shifts: npt.ArrayLike) -> F64:
    """S_k = sum of ``r = d - p`` over the shifted set (I + k) mod N, for every shift k."""
    rv = np.asarray(r, dtype=np.float64).ravel()
    return np.asarray(rv[_sets(idx, shifts, rv.size)].sum(axis=1), dtype=np.float64)


def _student(s: F64, w: F64) -> F64:
    """s / sqrt(w); w = 0 gives 0 for s = 0 and +-inf otherwise (an extreme value)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(w > 0, s / np.sqrt(w), np.where(s == 0, 0.0, np.sign(s) * np.inf))
    return np.asarray(t, dtype=np.float64)


def studentised_null(
    r: npt.ArrayLike, v: npt.ArrayLike, idx: npt.ArrayLike, shifts: npt.ArrayLike
) -> F64:
    """ADR-0018 decision 2: T_k = sum r / sqrt(sum v) over (I + k) mod N, with ``r = d - p`` and
    ``v = p (1 - p)``, for every shift k."""
    rv = np.asarray(r, dtype=np.float64).ravel()
    vv = np.asarray(v, dtype=np.float64).ravel()
    sets = _sets(idx, shifts, rv.size)
    return _student(rv[sets].sum(axis=1), vv[sets].sum(axis=1))


def hl_null(a: npt.ArrayLike, det: npt.ArrayLike, idx: npt.ArrayLike, shifts: npt.ArrayLike) -> F64:
    """Hodges-Lehmann shift of the genes in (I + k) mod N against the replicon's other genes, both
    restricted to detected genes with a finite ``a``, for every shift k (NaN if a side is empty).

    Exact up to ~1e-14 relative: per shift, the k-th pairwise difference is selected by bisection
    on its value, counting pairs by binary search on the sorted values (inside-inside pairs
    subtracted), then by enumerating the few remaining candidates.
    """
    av = np.asarray(a, dtype=np.float64).ravel()
    ok = np.asarray(det, dtype=bool).ravel() & np.isfinite(av)
    sets = _sets(idx, shifts, av.size)
    out = np.full(sets.shape[0], np.nan)
    order = np.flatnonzero(ok)[np.argsort(av[ok], kind="stable")]
    zs = av[order]
    pos = np.full(av.size, -1, dtype=np.int64)
    pos[order] = np.arange(zs.size)
    step = max(1, _CHUNK // max(1, sets.shape[1] ** 2))
    for c in range(0, sets.shape[0], step):
        out[c : c + step] = _hl_rows(zs, pos[sets[c : c + step]])
    return out


def _hl_rows(zs: F64, rp: I64) -> F64:
    """HL per row of rank positions ``rp`` (-1 = not observed) against the rest of ``zs``."""
    valid = rp >= 0
    m = valid.sum(axis=1)
    total = m * (zs.size - m)
    out = np.full(rp.shape[0], np.nan)
    live = np.flatnonzero(total > 0)
    if live.size == 0:
        return out
    x = np.where(valid[live], zs[np.maximum(rp[live], 0)], np.nan)
    v = valid[live]
    k1, k2 = (total[live] - 1) // 2, total[live] // 2
    lo = np.nanmin(x, axis=1) - zs[-1] - 1.0
    hi = np.nanmax(x, axis=1) - zs[0] + 1.0
    s1 = _kth(zs, x, v, rp[live], k1, lo, hi, total[live])
    s2 = _kth(zs, x, v, rp[live], k2, lo, hi, total[live])
    out[live] = np.where(k1 == k2, s1, (s1 + s2) / 2.0)
    return out


def _count(zs: F64, x: F64, v: B8, t: F64) -> I64:
    """Per row, #{(i inside, j outside): zs[j] >= x_i - t}, i.e. pairs with x_i - z_j <= t."""
    q = np.where(v, x - t[:, None], np.inf)
    c_all = (zs.size - np.searchsorted(zs, q, side="left")).sum(axis=1)
    xin = np.where(v, x, -np.inf)
    c_in = (xin[:, None, :] >= q[:, :, None]).sum(axis=(1, 2))
    return np.asarray(c_all - c_in, dtype=np.int64)


def _kth(
    zs: F64, x: F64, v: B8, rp: I64, k: I64, lo0: F64, hi0: F64, total: I64
) -> F64:  # fmt: skip
    """k-th smallest (0-based) inside-minus-outside difference per row."""
    lo, hi = lo0.copy(), hi0.copy()
    c_lo, c_hi = np.zeros_like(k), total.copy()
    res = np.full(k.size, np.nan)
    tol = 1e-14 * max(1.0, float(np.abs(zs).max()))
    act = np.arange(k.size)
    enum: list[I64] = []
    while act.size:
        small = c_hi[act] - c_lo[act] <= _ENUM
        tiny = ~small & (hi[act] - lo[act] <= tol)
        res[act[tiny]] = hi[act[tiny]]
        enum.append(act[small])
        act = act[~small & ~tiny]
        if not act.size:
            break
        mid = lo[act] + 0.5 * (hi[act] - lo[act])
        c = _count(zs, x[act], v[act], mid)
        up = c <= k[act]
        lo[act[up]], c_lo[act[up]] = mid[up], c[up]
        hi[act[~up]], c_hi[act[~up]] = mid[~up], c[~up]
    e = np.concatenate(enum)
    if e.size == 0:
        return res
    ve = v[e]
    j0 = np.searchsorted(zs, np.where(ve, x[e] - hi[e][:, None], np.inf), side="left")
    j1 = np.searchsorted(zs, np.where(ve, x[e] - lo[e][:, None], np.inf), side="left")
    lens = (j1 - j0).ravel()
    row = np.repeat(np.repeat(np.arange(e.size), x.shape[1]), lens)
    xi = np.repeat(np.where(ve, x[e], 0.0).ravel(), lens)
    j = np.repeat(j0.ravel() - (np.cumsum(lens) - lens), lens) + np.arange(int(lens.sum()))
    inside = (np.arange(e.size)[:, None] * zs.size + rp[e])[ve]
    keep = ~np.isin(row * zs.size + j, inside)
    row, d = row[keep], xi[keep] - zs[j[keep]]
    srt = np.lexsort((d, row))
    first = np.searchsorted(row[srt], np.arange(e.size), side="left")
    res[e] = d[srt][first + (k[e] - c_lo[e])]
    return res


def lag_autocovariance(
    a: npt.ArrayLike, ok: npt.ArrayLike, max_lag: int = LAG_MAX, *, circular: bool = True
) -> F64:
    """ADR-0018 amendment 2: gamma(0 ... ``max_lag``), the autocovariance of ``a`` over the genes
    with ``ok`` and a finite value at gene-index lag h: x = a - mean on those genes, 0 elsewhere;
    gamma(h) = sum_i x_i x_(i+h) / #{i: both ok} (0 without a pair), i + h taken mod N on a
    ``circular`` replicon and only inside the sequence on a linear one. Truncated at the first
    h >= 1 with gamma(h) <= 0 (gamma = 0 from there on). No ok gene: all 0."""
    av = np.asarray(a, dtype=np.float64).ravel()
    okv = np.asarray(ok, dtype=bool).ravel() & np.isfinite(av)
    g = np.zeros(max_lag + 1)
    if not okv.any():
        return g
    x = np.where(okv, av - av[okv].mean(), 0.0)
    w = okv.astype(np.float64)
    n = av.size
    for h in range(max_lag + 1):
        if circular:
            num, den = float(x @ np.roll(x, -h)), float(w @ np.roll(w, -h))
        elif h < n:
            num, den = float(x[: n - h] @ x[h:]), float(w[: n - h] @ w[h:])
        else:
            num, den = 0.0, 0.0
        g[h] = num / den if den > 0 else 0.0
    bad = np.flatnonzero(g[1:] <= 0)
    if bad.size:
        g[1 + int(bad[0]) :] = 0.0
    return g


def _dist(i: I64, j: I64, n: int, circular: bool) -> I64:
    d = np.abs(i - j)
    return np.minimum(d, n - d) if circular else d


def hl_variance(
    ok: npt.ArrayLike, idx: npt.ArrayLike, shifts: npt.ArrayLike, gamma: npt.ArrayLike
) -> F64:
    """ADR-0018 amendment 2: V_k = Var(mean_in - mean_out) = S_in,in / m^2 + S_out,out / m'^2
    - 2 S_in,out / (m m') for every shift k, S_A,B the sum of gamma(d(i, j)) over ``ok`` genes
    i in A and j in B, d the circular gene-index distance (gamma = 0 beyond its length), A the set
    (I + k) mod N; NaN if m or m' is 0. No TU label enters V. ``test_region`` uses the linear
    distance |i - j| on linear replicons (decision 16)."""
    return _lag_variance(ok, idx, shifts, gamma, True)


def _lag_variance(
    ok: npt.ArrayLike, idx: npt.ArrayLike, shifts: npt.ArrayLike, gamma: npt.ArrayLike,
    circular: bool,
) -> F64:  # fmt: skip
    """``hl_variance`` on a circular or linear replicon."""
    okv = np.asarray(ok, dtype=bool).ravel()
    g = np.asarray(gamma, dtype=np.float64).ravel()
    n, h_max = okv.size, g.size - 1
    w = okv.astype(np.float64)
    # c_i = sum_j ok_j gamma(d(i, j)), each offset once (offsets with distance <= h_max)
    c = np.zeros(n)
    pos = np.arange(n, dtype=np.int64)
    lim = min(h_max, n - 1)
    offs = sorted({o % n for o in range(-lim, lim + 1)}) if circular else range(-lim, lim + 1)
    for o in offs:
        dist = min(o, n - o) if circular else abs(o)
        if dist > h_max or g[dist] == 0.0:
            continue
        j = pos + o
        if circular:
            c += g[dist] * w[j % n]
        else:
            inb = (j >= 0) & (j < n)
            c[inb] += g[dist] * w[j[inb]]
    s_all, m_all = float(w @ c), float(w.sum())
    sets = _sets(idx, shifts, n)
    out = np.full(sets.shape[0], np.nan)
    step = max(1, _CHUNK // max(1, sets.shape[1] ** 2))
    for b in range(0, sets.shape[0], step):
        st = sets[b : b + step]
        om = okv[st]
        m = om.sum(axis=1).astype(np.float64)
        s_in_all = np.where(om, c[st], 0.0).sum(axis=1)
        d = _dist(st[:, :, None], st[:, None, :], n, circular)
        gm = np.where(d <= h_max, g[np.minimum(d, h_max)], 0.0)
        s_ii = (gm * (om[:, :, None] & om[:, None, :])).sum(axis=(1, 2))
        s_io, s_oo = s_in_all - s_ii, s_all - 2.0 * s_in_all + s_ii
        mo = m_all - m
        live = (m > 0) & (mo > 0)
        mm, mmo = np.where(live, m, 1.0), np.where(live, mo, 1.0)
        v = s_ii / mm**2 + s_oo / mmo**2 - 2.0 * s_io / (mm * mmo)
        out[b : b + step] = np.where(live, v, np.nan)
    return out


def _hl_t(
    a: F64, ok: B8, idx: npt.ArrayLike, shifts: npt.ArrayLike, scale: float | None, circular: bool
) -> tuple[F64, bool]:
    """(T per shift, whether some finite V <= 0 occurred with gamma(0) above the guard)."""
    hl = hl_null(a, ok, idx, shifts)
    if not ok.any():
        return hl, False
    g = lag_autocovariance(a, ok, circular=circular)
    mag = float(np.abs(a[ok]).max()) if scale is None else scale
    v = _lag_variance(ok, idx, shifts, g, circular)
    na = np.isnan(hl) | np.isnan(v)
    if g[0] <= (V_RTOL * mag) ** 2:  # decision 12 guard: rounding noise, V counts as 0, T = 0
        return np.where(na, np.nan, 0.0), False
    nonpos = bool(np.any(v[~na] <= 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(v > 0, hl / np.sqrt(v), np.where(hl == 0, 0.0, np.sign(hl) * np.inf))
    return np.asarray(np.where(na, np.nan, t), dtype=np.float64), nonpos


def studentised_hl_null(
    a: npt.ArrayLike,
    ok: npt.ArrayLike,
    idx: npt.ArrayLike,
    shifts: npt.ArrayLike,
    *,
    scale: float | None = None,
    circular: bool = True,
) -> F64:
    """ADR-0018 amendment 2: T_k = HL_k / sqrt(V_k) for every shift k (k = 0 is the observed T),
    gamma estimated once from (``a``, ``ok``) (``lag_autocovariance``), V from ``hl_variance``;
    NaN where a side is empty. Decision 12: V counts as 0 (T = 0) when gamma(0) <=
    (V_RTOL * ``scale``)^2, ``scale`` the magnitude of the data before residualisation (default
    max |a| over ``ok``). V <= 0 otherwise gives T = 0 for HL = 0 and +-inf else; ``test_region``
    reports such a region's p as NA (decision 16)."""
    av = np.asarray(a, dtype=np.float64).ravel()
    okv = np.asarray(ok, dtype=bool).ravel() & np.isfinite(av)
    return _hl_t(av, okv, idx, shifts, scale, circular)[0]


def _flags(n_eligible: int, k: int, extra: tuple[str, ...]) -> tuple[str, ...]:
    out = ("small_region",) if n_eligible < SMALL_REGION else ()
    return (*out, *(("few_shifts",) if k < FEW_SHIFTS_K else ()), *extra)


def test_region(
    genes: pd.DataFrame,
    region_genes: Sequence[str],
    region_type: str,
    contig: str,
    cfg: GippaConfig,
    rng: np.random.Generator,
    *,
    region_id: str = "",
    circular: bool = True,
    flags: tuple[str, ...] = (),
) -> list[RegionTest]:
    """O - E and HL rows of one region on one replicon (``circular`` topology); ``flags`` are
    appended to each row's own flags (e.g. ``split_contigs``)."""
    return _test_region(genes, region_genes, region_type, contig, cfg, rng, region_id=region_id,
                        circular=circular, flags=flags)[0]  # fmt: skip


def _test_region(
    genes: pd.DataFrame,
    region_genes: Sequence[str],
    region_type: str,
    contig: str,
    cfg: GippaConfig,
    rng: np.random.Generator,
    *,
    region_id: str,
    circular: bool,
    flags: tuple[str, ...],
) -> tuple[list[RegionTest], list[RegionNull]]:
    """``test_region`` and the nulls of its rows with at least one null value."""
    on_contig = genes["contig"].astype(str) == contig
    if genes.loc[on_contig, "tu_id"].isna().any():  # decision 13: C6 makes tu_id non-null
        raise ValueError(f"missing tu_id on replicon {contig!r} (C6: every gene has a tu_id)")
    ix = replicon_index(genes, contig)
    members = set(map(str, region_genes))
    n_genes = int(genes.loc[on_contig, "locus_tag"].astype(str).isin(members).sum())
    inside = ix["locus_tag"].isin(members).to_numpy(bool)
    idx = np.flatnonzero(inside)
    n, m = len(ix), int(idx.size)
    base = (region_id, region_type, contig)
    na = math.nan
    if m == 0:
        f = _flags(0, 0, flags)
        return [
            RegionTest(*base, "o_minus_e", na, na, na, na, na, na, 0, n_genes, 0, f),
            RegionTest(*base, "hl_shift", na, None, na, na, na, na, 0, n_genes, 0, f),
        ], []
    valid = circular_shifts(n, idx, circular)
    shifts, _ = choose_shifts(
        valid, cfg.context.exhaustive_shift_max_genes, n, cfg.stats.permutations, rng
    )
    B = cfg.stats.bootstrap
    d, p, a = ix["d"].to_numpy(float), ix["p"].to_numpy(float), ix["a"].to_numpy(float)
    r, v = d - p, p * (1.0 - p)
    s_obs = float(r[idx].sum())
    t_obs = float(_student(np.array([s_obs]), np.array([v[idx].sum()]))[0])
    null_oe = studentised_null(r, v, idx, shifts)
    p_oe = perm_p(t_obs, null_oe) if shifts.size else na
    nulls: list[RegionNull] = []
    if shifts.size:
        nulls.append(RegionNull(region_id, contig, "o_minus_e", shifts, null_oe, t_obs))
    # bootstrap frames hold one float64 block (TU codes as floats): row takes stay cheap
    tu = pd.factorize(ix["tu_id"])[0].astype(np.float64)
    df = pd.DataFrame({"tu": tu[idx], "r": r[idx]})
    ci = tu_block_bootstrap(df, "tu", lambda t: float(t["r"].to_numpy().sum()), B, rng)
    oe = RegionTest(*base, "o_minus_e", float(d[idx].sum()), float(p[idx].sum()), s_obs,
                    ci.ci_low, ci.ci_high, p_oe, int(shifts.size), n_genes, m,
                    _flags(m, int(shifts.size), flags))  # fmt: skip

    raw = a[(d > 0) & np.isfinite(a)]
    a = residualise(a, d > 0, ix)
    ok = (d > 0) & np.isfinite(a)
    hl = hodges_lehmann(a[inside & ok], a[~inside & ok])
    k_hl, p_hl, lo, hi = 0, na, na, na
    hflag: tuple[str, ...] = ()
    if not math.isnan(hl):
        mag = float(np.abs(raw).max())
        t_all, nonpos = _hl_t(a, ok, idx, np.concatenate([[0], shifts]), mag, circular)
        null = t_all[1:]
        k_hl = int(np.count_nonzero(~np.isnan(null)))
        p_hl = perm_p(float(t_all[0]), null) if k_hl and not nonpos else na
        hflag = ("hl_variance_nonpositive",) if nonpos else ()
        if k_hl and not nonpos:  # no null to plot without a p (M9-17 decision 8)
            nulls.append(RegionNull(region_id, contig, "hl_shift", shifts, null, float(t_all[0])))
        hd = pd.DataFrame({"tu": tu[ok], "a": a[ok], "in": inside[ok].astype(np.float64)})
        hci = tu_block_bootstrap(hd, "tu", _hl_of, B, rng, strata_col="in")
        lo, hi = hci.ci_low, hci.ci_high
    hrow = RegionTest(*base, "hl_shift", hl, None, hl, lo, hi, p_hl, k_hl, n_genes, m,
                      _flags(m, k_hl, (*hflag, *flags)))  # fmt: skip
    return [oe, hrow], nulls


def _hl_of(t: pd.DataFrame) -> float:
    a, w = t["a"].to_numpy(float), t["in"].to_numpy(float) > 0.5
    return hodges_lehmann(a[w], a[~w])


def test_all_regions(
    genes: pd.DataFrame,
    regions: pd.DataFrame,
    region_genes: pd.DataFrame,
    replicons: pd.DataFrame,
    cfg: GippaConfig,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Region tests of every prophage, BGC, CRISPR and AMR region (sorted ``region_id``; a region
    with members on several replicons gets one test per replicon, flagged ``split_contigs``),
    with ``q_value`` (BH within (region_type, stat)), ``method`` and ``null_model``."""
    return region_tests(genes, regions, region_genes, replicons, cfg, rng)[0]


def region_tests(
    genes: pd.DataFrame,
    regions: pd.DataFrame,
    region_genes: pd.DataFrame,
    replicons: pd.DataFrame,
    cfg: GippaConfig,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, list[RegionNull]]:
    """``test_all_regions`` and every tested row's null (``RegionNull``, test order)."""
    types = [t for t in REGION_TYPES if t in cfg.context.region_types]
    reg = regions[regions["region_type"].isin(types)]
    rtype = dict(zip(reg["region_id"].astype(str), reg["region_type"].astype(str), strict=True))
    rcontig = dict(zip(reg["region_id"].astype(str), reg["contig"].astype(str), strict=True))
    contig_of = dict(zip(genes["locus_tag"].astype(str), genes["contig"].astype(str), strict=True))
    topo = dict(zip(replicons["contig"].astype(str), replicons["topology"].astype(str),
                    strict=True))  # fmt: skip
    rank = {c: i for i, c in enumerate(replicons["contig"].astype(str))}
    rg = region_genes[region_genes["region_id"].astype(str).isin(rtype)]
    by_region: dict[str, list[str]] = {rid: [] for rid in rtype}
    for rid, tag in zip(rg["region_id"].astype(str), rg["locus_tag"].astype(str), strict=True):
        by_region[rid].append(tag)
    rows: list[RegionTest] = []
    nulls: list[RegionNull] = []
    for rid in sorted(by_region):
        tags = sorted(set(by_region[rid]))
        contigs = sorted({contig_of[t] for t in tags if t in contig_of},
                         key=lambda c: (rank.get(c, len(rank)), c)) or [rcontig[rid]]  # fmt: skip
        extra = ("split_contigs",) if len(contigs) > 1 else ()
        for c in contigs:
            got, null = _test_region(genes, tags, rtype[rid], c, cfg, rng, region_id=rid,
                                     circular=topo.get(c) == "circular", flags=extra)  # fmt: skip
            rows += got
            nulls += null
    names = [f.name for f in fields(RegionTest)]
    out = pd.DataFrame([astuple(r) for r in rows], columns=names)
    floats = ("observed", "expected", "effect", "ci_low", "ci_high", "p_value")
    out = out.astype({**dict.fromkeys(floats, "float64"),
                      **dict.fromkeys(("n_perm", "n_genes", "n_eligible"), "int64")})  # fmt: skip
    out["q_value"] = np.nan
    for _, grp in out.groupby(["region_type", "stat"], sort=False):
        out.loc[grp.index, "q_value"] = bh(grp["p_value"].to_numpy(float))
    out["method"] = out["stat"].map(METHOD)
    out["null_model"] = "circular_shift"
    return out, nulls
