"""Region tests: O - E and Hodges-Lehmann shift under circular gene-index shifts (issue M4-17).

GIAP_SPEC.md §7.3.5(b), §8.1, Appendix B.4; docs/specs/M4-17-region-tests.md. The index of a
replicon is its eligible genes (not partial, with a ``p_detect``) in genome order, 0 ... N - 1; a
region's index set I is moved to (I + k) mod N for every valid shift k (``circular_shifts``,
``choose_shifts``), both statistics are recomputed at each shift and the two-sided centred
``perm_p`` is reported with a TU-block bootstrap CI. q-values: BH within (region_type, stat).
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

from giap.config.model import GiapConfig
from giap.stats.bootstrap import tu_block_bootstrap
from giap.stats.multitest import bh
from giap.stats.permutation import choose_shifts, circular_shifts, perm_p
from giap.stats.robust import hodges_lehmann

__all__ = ["REGION_TYPES", "RegionTest", "hl_null", "o_minus_e_null", "replicon_index",
           "test_all_regions", "test_region"]  # fmt: skip

#: region types tested here (§7.3.5(b)); ``plasmid`` regions belong to M4-18
REGION_TYPES = ("prophage", "bgc", "crispr", "amr")
METHOD = {"o_minus_e": "circular_shift_o_minus_e", "hl_shift": "circular_shift_hodges_lehmann"}
FEW_SHIFTS_K = 19  # K < 19: the smallest attainable p, 1/(1 + K), exceeds 0.05
SMALL_REGION = 3
_ENUM = 64  # HL null: enumerate a shift's candidate differences once at most this many remain
_CHUNK = 4_000_000  # HL null: shifts per chunk such that K * m * m stays below this

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


def replicon_index(genes: pd.DataFrame, contig: str) -> pd.DataFrame:
    """Eligible genes of ``contig`` in genome order: ``idx`` (0 ... N - 1), ``locus_tag``,
    ``tu_id``, ``d`` (detected, 0/1), ``p`` (p_detect) and ``a`` (log2 abundance of detected
    genes, NaN otherwise). Eligible = not partial and ``p_detect`` present."""
    g = genes[genes["contig"].astype(str) == contig]
    p = pd.to_numeric(g["p_detect"], errors="coerce").astype(float)
    g = g[~g["partial"].astype("boolean").fillna(False).to_numpy(bool) & p.notna().to_numpy()]
    g = g.sort_values(["start", "locus_tag"], kind="mergesort")
    d = g["detected"].astype("boolean").fillna(False).to_numpy(bool)
    a = pd.to_numeric(g["abundance_log2"], errors="coerce").to_numpy(float)
    return pd.DataFrame(
        {
            "idx": np.arange(len(g), dtype=np.int64),
            "locus_tag": g["locus_tag"].astype(str).to_numpy(),
            "tu_id": g["tu_id"].to_numpy(),
            "d": d.astype(np.int64),
            "p": pd.to_numeric(g["p_detect"]).to_numpy(float),
            "a": np.where(d, a, np.nan),
        }
    )


def _sets(idx: npt.ArrayLike, shifts: npt.ArrayLike, n: int) -> I64:
    i = np.asarray(idx, dtype=np.int64).ravel()
    k = np.asarray(shifts, dtype=np.int64).ravel()
    return np.asarray((i[None, :] + k[:, None]) % n, dtype=np.int64)


def o_minus_e_null(r: npt.ArrayLike, idx: npt.ArrayLike, shifts: npt.ArrayLike) -> F64:
    """S_k = sum of ``r = d - p`` over the shifted set (I + k) mod N, for every shift k."""
    rv = np.asarray(r, dtype=np.float64).ravel()
    return np.asarray(rv[_sets(idx, shifts, rv.size)].sum(axis=1), dtype=np.float64)


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


def _flags(n_eligible: int, k: int, extra: tuple[str, ...]) -> tuple[str, ...]:
    out = ("small_region",) if n_eligible < SMALL_REGION else ()
    return (*out, *(("few_shifts",) if k < FEW_SHIFTS_K else ()), *extra)


def test_region(
    genes: pd.DataFrame,
    region_genes: Sequence[str],
    region_type: str,
    contig: str,
    cfg: GiapConfig,
    rng: np.random.Generator,
    *,
    region_id: str = "",
    circular: bool = True,
    flags: tuple[str, ...] = (),
) -> list[RegionTest]:
    """O - E and HL rows of one region on one replicon (``circular`` topology); ``flags`` are
    appended to each row's own flags (e.g. ``split_contigs``)."""
    ix = replicon_index(genes, contig)
    members = set(map(str, region_genes))
    on_contig = genes["contig"].astype(str) == contig
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
        ]
    valid = circular_shifts(n, idx, circular)
    shifts, _ = choose_shifts(
        valid, cfg.context.exhaustive_shift_max_genes, n, cfg.stats.permutations, rng
    )
    B = cfg.stats.bootstrap
    d, p, a = ix["d"].to_numpy(float), ix["p"].to_numpy(float), ix["a"].to_numpy(float)
    r = d - p
    s_obs = float(r[idx].sum())
    p_oe = perm_p(s_obs, o_minus_e_null(r, idx, shifts)) if shifts.size else na
    # bootstrap frames hold one float64 block (TU codes as floats): row takes stay cheap
    tu = pd.factorize(ix["tu_id"])[0].astype(np.float64)
    df = pd.DataFrame({"tu": tu[idx], "r": r[idx]})
    ci = tu_block_bootstrap(df, "tu", lambda t: float(t["r"].to_numpy().sum()), B, rng)
    oe = RegionTest(*base, "o_minus_e", float(d[idx].sum()), float(p[idx].sum()), s_obs,
                    ci.ci_low, ci.ci_high, p_oe, int(shifts.size), n_genes, m,
                    _flags(m, int(shifts.size), flags))  # fmt: skip

    ok = (d > 0) & np.isfinite(a)
    hl = hodges_lehmann(a[inside & ok], a[~inside & ok])
    k_hl, p_hl, lo, hi = 0, na, na, na
    if not math.isnan(hl):
        null = hl_null(a, ok, idx, shifts)
        k_hl = int(np.count_nonzero(~np.isnan(null)))
        p_hl = perm_p(hl, null) if k_hl else na
        hd = pd.DataFrame({"tu": tu[ok], "a": a[ok], "in": inside[ok].astype(np.float64)})
        hci = tu_block_bootstrap(hd, "tu", _hl_of, B, rng, strata_col="in")
        lo, hi = hci.ci_low, hci.ci_high
    hrow = RegionTest(*base, "hl_shift", hl, None, hl, lo, hi, p_hl, k_hl, n_genes, m,
                      _flags(m, k_hl, flags))  # fmt: skip
    return [oe, hrow]


def _hl_of(t: pd.DataFrame) -> float:
    a, w = t["a"].to_numpy(float), t["in"].to_numpy(float) > 0.5
    return hodges_lehmann(a[w], a[~w])


def test_all_regions(
    genes: pd.DataFrame,
    regions: pd.DataFrame,
    region_genes: pd.DataFrame,
    replicons: pd.DataFrame,
    cfg: GiapConfig,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Region tests of every prophage, BGC, CRISPR and AMR region (sorted ``region_id``; a region
    with members on several replicons gets one test per replicon, flagged ``split_contigs``),
    with ``q_value`` (BH within (region_type, stat)), ``method`` and ``null_model``."""
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
    for rid in sorted(by_region):
        tags = sorted(set(by_region[rid]))
        contigs = sorted({contig_of[t] for t in tags if t in contig_of},
                         key=lambda c: (rank.get(c, len(rank)), c)) or [rcontig[rid]]  # fmt: skip
        extra = ("split_contigs",) if len(contigs) > 1 else ()
        for c in contigs:
            rows += test_region(genes, tags, rtype[rid], c, cfg, rng, region_id=rid,
                                circular=topo.get(c) == "circular", flags=extra)  # fmt: skip
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
    return out
