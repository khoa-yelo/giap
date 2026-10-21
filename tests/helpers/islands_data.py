"""Inputs, oracles and the step driver for M5-06 islands (docs/specs/M5-06-islands.md).

Only the last section ("driving the step") calls ``giap.context.islands``: the oracles are
written from the spec text (§7.3.5(d) steps 1-10, lead decisions 21-27 and 40-42, ADR-0021) on
top of the *merged* engines
``giap.stats.scores.blom_scores``, ``giap.stats.smooth.loess_detrend``,
``giap.stats.autocorr.choose_block_len`` / ``block_len_prime``, ``giap.stats.scan`` and
``giap.stats.bootstrap``, each tested by its own issue.

Inputs. Hand-built gene tables reuse the M4-17 builder (``region_tests_data.gene_table``,
§9.3 ``genes.tsv`` columns); ``tu_frame`` derives the C5 ``tus`` frame from them (0-based
half-open, like every internal table). ``sim_contig`` draws a larger replicon: TU-level AR(1)
log2 abundance, gene noise, an optional ori-ter-like gradient and detection independent of
abundance (the island null).

Interface assumed of ``giap.context.islands`` (spec "Public interface" plus what the reported
coordinates need; see the test module docstrings):

``TuSeries(contig, circular, length, tu_ids, genes, start, end, score, detrended)``
    ``length`` = replicon length (C1) and per-TU ``start`` / ``end`` (C5, 0-based half-open) are
    needed for the unwrapped island coordinates of step 7; everything else is the spec's.
``gene_scores(ix)``, ``tu_series(ix, tus, contig=, circular=, length=, span_frac=)``,
``block_length(ts, cfg, rng)``, ``scan_replicons(series, cfg, rng) -> IslandScan``,
``islands_frame(scan, genes, cfg, rng)``, ``gene_island_map(islands)``.
``IslandScan(series, scanned, block_len, windows, p_value, clusters, null_max, n_perm, flags)``.
"""

from __future__ import annotations

import math
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import region_tests_data as rd
from region_tests_data import Contig

#: §9.3 ``islands.tsv`` columns, the §9.2 test columns and ``flags``
ISLAND_COLUMNS = (
    "island_id", "direction", "contig", "start", "end", "tu_first", "tu_last", "n_tus",
    "n_genes", "genes", "mean_abundance_log2", "detected_frac", "scan_s", "block_len",
    "stat", "effect", "ci_low", "ci_high", "p_value", "q_value", "method", "null_model",
    "n_perm", "flags",
)  # fmt: skip
GENE_MAP_COLUMNS = ("locus_tag", "island_id")
METHOD = "scan_max_stat"
NULL_MODEL = "circular_block_bootstrap"
STEP = "context_stats.islands"  # §6.6 random stream of the island part
FEW_VALUES = "few_values"  # decision 40: replicon not scanned (span gives k < 4)
FEW_UNITS = "ci_few_units"  # M4-16: CI from < 10 resampling units
#: ADR-0021 Amendment 1 / decision 48: ``uncalibrated`` is dropped; the run report carries
#: l-hat, l' and the lag-1 Moran's I of the detrended TU scores per replicon instead.
RUN_REPORT_FIELDS = ("l_hat", "block_len", "morans_i_lag1")
MIN_SPAN_POINTS = 4  # decision 40/41: loess_span_points(T, frac) < 4 -> not scanned


def config(
    *,
    bootstrap: int = 199,
    permutations: int = 99,
    span_frac: float = 0.10,
    block_len: int | str = "auto",
    min_tus: int = 2,
    max_tus: int = 20,
    alpha: float = 0.05,
    seed: int = 7,
    abundance_missing: str | None = None,
) -> Any:
    """A ``GiapConfig`` with the §5.2 knobs the island step reads."""
    from giap.config.model import GiapConfig

    raw: dict[str, Any] = {
        "seed": seed,
        "stats": {
            "bootstrap": bootstrap,
            "permutations": permutations,
            "alpha": {"context": alpha},
        },
        "context": {
            "islands": {
                "min_tus": min_tus,
                "max_tus": max_tus,
                "detrend_span_frac": span_frac,
                "block_len": block_len,
            }
        },
    }
    if abundance_missing is not None:
        raw["abundance"] = {"missing": abundance_missing}
    return GiapConfig.model_validate(raw)


# --------------------------------------------------------------------------- input frames


def tu_frame(genes: pd.DataFrame) -> pd.DataFrame:
    """C5 ``tus`` from a gene table: one row per (contig, tu_id), 0-based half-open span of its
    genes, ``genes`` in genome order. These spans never cross the origin; :func:`cross_origin`
    rewrites a TU with genes at both ends of a circular replicon into the C5 crossing form."""
    rows = []
    for (contig, tu), g in genes.groupby(["contig", "tu_id"], sort=False):
        rows.append(
            (
                str(tu),
                str(contig),
                int(g["start"].min()),
                int(g["end"].max()),
                str(g["strand"].iloc[0]),
                [str(t) for t in g.sort_values("start")["locus_tag"]],
            )
        )
    out = pd.DataFrame(rows, columns=["tu_id", "contig", "start", "end", "strand", "genes"])
    out = out.sort_values(["contig", "start"], kind="mergesort").reset_index(drop=True)
    return out.astype({"start": "int64", "end": "int64"})


def cross_origin(
    genes: pd.DataFrame, tus: pd.DataFrame, replicons: pd.DataFrame, contig: str, tu_id: str
) -> pd.DataFrame:
    """``tus`` with TU ``tu_id`` of circular ``contig`` in the C5 origin-crossing form (C5
    ``build_tus`` links the last and first genes; spec 6.5 unwrapped): ``start`` is that of
    its genes past mid-replicon, ``end`` the end of its genes before it plus the length."""
    length = int(replicons.loc[replicons["contig"] == contig, "length"].iloc[0])
    g = genes[(genes["contig"].astype(str) == contig) & (genes["tu_id"].astype(str) == tu_id)]
    tail = g["start"].astype(int) >= length // 2
    assert tail.any() and (~tail).any(), "the TU must hold genes at both ends of the replicon"
    out = tus.copy()
    m = (out["contig"].astype(str) == contig) & (out["tu_id"].astype(str) == tu_id)
    out.loc[m, "start"] = int(g.loc[tail, "start"].min())
    out.loc[m, "end"] = int(g.loc[~tail, "end"].max()) + length
    out.loc[m, "genes"] = pd.Series(
        [[str(t) for t in pd.concat([g[tail], g[~tail]])["locus_tag"]]], index=out.index[m]
    )
    out = out.sort_values(["contig", "start"], kind="mergesort").reset_index(drop=True)
    return out.astype({"start": "int64", "end": "int64"})


def frames(contigs: Sequence[Contig]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """(genes, tus, replicons) for ``contigs`` (M4-17 ``gene_table`` + :func:`tu_frame`)."""
    genes, replicons = rd.gene_table(contigs)
    return genes, tu_frame(genes), replicons


def tu_labels(sizes: Sequence[int], prefix: str) -> list[str]:
    """TU label per gene for consecutive TUs of the given sizes."""
    return [f"{prefix}_TU_{t + 1}" for t, n in enumerate(sizes) for _ in range(n)]


def geometric_tu_sizes(rng: np.random.Generator, n_genes: int, p: float) -> list[int]:
    """TU sizes ``1 + Geometric(p) - 1`` summing to ``n_genes`` (mean 1/p)."""
    sizes: list[int] = []
    while sum(sizes) < n_genes:
        sizes.append(int(rng.geometric(p)))
    sizes[-1] -= sum(sizes) - n_genes
    if sizes[-1] == 0:
        sizes.pop()
    return sizes


def ar1(rng: np.random.Generator, n: int, rho: float, sd: float = 1.0) -> np.ndarray:
    """Stationary AR(1) series of length ``n``."""
    e = rng.standard_normal(n)
    x = np.empty(n)
    if n:
        x[0] = sd * e[0] / math.sqrt(1.0 - rho * rho)
        for t in range(1, n):
            x[t] = rho * x[t - 1] + sd * e[t]
    return x


def oriter_gradient(n_tus: int, amplitude: float, shape: str = "cosine") -> np.ndarray:
    """Ori-ter trend over TU order: +amplitude/2 at the origin, -amplitude/2 at L/2.

    ``cosine`` is the smooth version (what LOESS detrending is meant to remove), ``triangular``
    the literal replichore-distance version with kinks at the origin and the terminus.
    """
    u = np.arange(n_tus, dtype=float) / max(n_tus, 1)
    if shape == "triangular":
        return amplitude * (0.5 - 2.0 * np.minimum(u, 1.0 - u))
    if shape == "cosine":
        return 0.5 * amplitude * np.cos(2.0 * np.pi * u)
    raise ValueError(f"unknown gradient shape {shape!r}")


def sim_contig(
    rng: np.random.Generator,
    name: str,
    n_genes: int,
    *,
    rho: float = 0.0,
    tu_p: float = 1.0,
    detected_frac: float = 0.55,
    circular: bool = True,
    gradient: float = 0.0,
    gradient_shape: str = "cosine",
    island: tuple[int, int, float] | None = None,
) -> Contig:
    """A replicon under the island null: TU-level AR(1) log2 abundance (SD 1), gene noise 0.5,
    an optional ori-ter gradient and detection ``Bernoulli(detected_frac)`` independent of
    abundance. ``island = (first_tu, n_tus, effect)`` plants a shift of ``effect`` SD."""
    sizes = [1] * n_genes if tu_p >= 1.0 else geometric_tu_sizes(rng, n_genes, tu_p)
    tus = tu_labels(sizes, name)
    x = ar1(rng, len(sizes), rho) + oriter_gradient(len(sizes), gradient, gradient_shape)
    if island is not None:
        first, n_isl, eff = island
        x[first : first + n_isl] += eff
    a = 20.0 + np.repeat(x, sizes) + 0.5 * rng.standard_normal(n_genes)
    d = rng.random(n_genes) < detected_frac
    return Contig(
        name=name,
        d=[int(v) for v in d],
        p=[detected_frac] * n_genes,
        a=[float(v) for v in a],
        circular=circular,
        tu=tus,
    )


def scenario_contigs(rng: Any, pair: Any) -> list[Contig]:
    """A giapsim ``scenarios`` pair ``(LayoutConfig, ExpressionConfig)`` drawn as ``Contig``s.

    The layout and the expression truth are drawn from ``rng``; ``detected_true``,
    ``p_detect_true`` and ``abundance_log2_true`` become the gene table's evidence columns, as
    ``region_tests_data.tables_from`` does with ``p_source = "true"``.
    """
    import giapsim.expression as ex
    import giapsim.layout as la

    lcfg, ecfg = pair
    lay = la.simulate_layout(lcfg, rng)
    e = ex.simulate_expression(lay, ecfg, rng)
    g = lay.genes[["locus_tag", "contig", "start", "partial"]].copy()
    eg = e.genes.set_index("locus_tag").loc[g["locus_tag"]]
    tu_of = {gg: t for t, gs in zip(lay.tus["tu_id"], lay.tus["genes"], strict=True) for gg in gs}
    g["tu_id"] = g["locus_tag"].map(tu_of)
    topo = dict(zip(lay.replicons["contig"], lay.replicons["topology"], strict=True))
    out = []
    for name, gc in g.groupby("contig", sort=False):
        gc = gc.sort_values("start", kind="mergesort")
        tags = gc["locus_tag"]
        out.append(
            Contig(
                name=str(name),
                d=[int(v) for v in eg.loc[tags, "detected_true"].astype(int)],
                p=[float(v) for v in eg.loc[tags, "p_detect_true"]],
                a=[float(v) for v in eg.loc[tags, "abundance_log2_true"]],
                circular=str(topo[name]) == "circular",
                partial=[bool(v) for v in gc["partial"]],
                tu=[str(t) for t in gc["tu_id"]],
            )
        )
    return out


# --------------------------------------------------------------------------- score oracles


def censored_share(genes: pd.DataFrame, contig: str) -> float:
    """q = share of the scored genes of ``contig`` that are undetected (left-censored)."""
    t = scored_table(genes, contig)
    return float((~t["det"].to_numpy(bool)).mean()) if len(t) else math.nan


def ceiling_from_q(q: float, w: int) -> float:
    """Decision 50: the most negative S a window of ``w`` TUs can reach when every one of its
    genes is undetected, ``sqrt(w) * |Phi^-1(q / 2)|`` (before the LOESS attenuation)."""
    from scipy.special import ndtri

    if not 0.0 < q < 1.0:
        return math.nan
    return float(math.sqrt(w) * abs(ndtri(q / 2.0)))


def max_window(ts: Any, cfg: Any) -> int:
    """Longest window the scan uses on ``ts``: ``min(max_tus, T - 1)`` circular, ``min(.., T)``."""
    n = len(ts.tu_ids)
    return int(min(cfg.context.islands.max_tus, n - 1 if ts.circular else n))


def scored_table(genes: pd.DataFrame, contig: str) -> pd.DataFrame:
    """Scored genes of ``contig`` in genome order (spec step 1).

    Eligible = not partial and ``p_detect`` present. A detected gene with NaN abundance is not
    scored (decision 22); an undetected eligible gene is scored and censored.
    """
    g = genes[genes["contig"].astype(str) == contig].copy()
    p = pd.to_numeric(g["p_detect"], errors="coerce")
    elig = ~g["partial"].astype("boolean").fillna(False).to_numpy(bool) & p.notna().to_numpy()
    g = g[elig].sort_values(["start", "locus_tag"], kind="mergesort")
    det = g["detected"].astype("boolean").fillna(False).to_numpy(bool)
    a = pd.to_numeric(g["abundance_log2"], errors="coerce").to_numpy(float)
    keep = ~det | np.isfinite(a)
    out = g[keep].copy()
    out["det"] = det[keep]
    out["a"] = np.where(det[keep], a[keep], np.nan)
    return out.reset_index(drop=True)


def oracle_gene_scores(genes: pd.DataFrame, contig: str) -> pd.Series:
    """Blom scores of the scored genes of ``contig``, indexed by ``locus_tag`` (spec step 1)."""
    from giap.stats.scores import blom_scores

    t = scored_table(genes, contig)
    z = blom_scores(np.where(t["det"].to_numpy(bool), t["a"].to_numpy(float), 0.0),
                    ~t["det"].to_numpy(bool))  # fmt: skip
    return pd.Series(z, index=pd.Index(t["locus_tag"].astype(str), name="locus_tag"))


@dataclass(frozen=True)
class OracleSeries:
    """What the oracle knows about one replicon (mirrors ``TuSeries``)."""

    contig: str
    circular: bool
    length: int
    tu_ids: tuple[str, ...]
    genes: tuple[tuple[str, ...], ...]
    start: np.ndarray
    end: np.ndarray
    score: np.ndarray
    detrended: np.ndarray
    scanned: bool
    span_points: int


def oracle_series(
    genes: pd.DataFrame, tus: pd.DataFrame, replicons: pd.DataFrame, cfg: Any
) -> list[OracleSeries]:
    """One :class:`OracleSeries` per replicon, replicons in C1 order (spec steps 2-3)."""
    from giap.stats.smooth import loess_detrend, loess_span_points

    frac = cfg.context.islands.detrend_span_frac
    tu_idx = tus.set_index(["contig", "tu_id"])
    out = []
    for _, rep in replicons.iterrows():
        contig = str(rep["contig"])
        circular = str(rep["topology"]) == "circular"
        t = scored_table(genes, contig)
        members: dict[str, list[str]] = {}
        for tag, tu in zip(t["locus_tag"].astype(str), t["tu_id"].astype(str), strict=True):
            members.setdefault(tu, []).append(tag)
        # Genome order of the TUs: C5 (start, tu_id) (spec step 2); an origin-crossing TU is last.
        order = sorted(members, key=lambda tu: (int(tu_idx.loc[(contig, tu), "start"]), tu))
        scores = oracle_gene_scores(genes, contig)
        z = np.array([float(scores[members[tu]].mean()) for tu in order])
        n = z.size
        k = loess_span_points(n, frac) if n else 0
        scanned = n >= 2 * cfg.context.islands.min_tus and k >= MIN_SPAN_POINTS
        if n:
            fit = loess_detrend(z, np.arange(n, dtype=float), frac, circular=circular,
                                period=float(n))  # fmt: skip
            resid = fit.resid
        else:
            resid = np.zeros(0)
        out.append(
            OracleSeries(
                contig=contig,
                circular=circular,
                length=int(rep["length"]),
                tu_ids=tuple(order),
                genes=tuple(tuple(members[tu]) for tu in order),
                start=np.array([int(tu_idx.loc[(contig, tu), "start"]) for tu in order], dtype=int),
                end=np.array([int(tu_idx.loc[(contig, tu), "end"]) for tu in order], dtype=int),
                score=z,
                detrended=resid,
                scanned=scanned,
                span_points=k,
            )
        )
    return out


# --------------------------------------------------------------------------- scan oracle


@dataclass(frozen=True)
class OracleScan:
    """Block lengths, the genome-wide null maximum and the per-replicon windows."""

    series: tuple[OracleSeries, ...]
    l_hat: tuple[int, ...]
    block_len: tuple[int, ...]
    morans_i_lag1: tuple[float, ...]
    null_max: np.ndarray
    windows: tuple[Any, ...]
    replicon_null: tuple[np.ndarray, ...]
    p_value: tuple[np.ndarray, ...]
    clusters: tuple[tuple[Any, ...], ...]
    n_perm: int
    rng: Any


def oracle_scan(series: Sequence[OracleSeries], cfg: Any, rng: Any) -> OracleScan:
    """Spec steps 4-6 with the draw order of step 10 (block lengths for every scanned replicon in
    C1 order, then the bootstrap draws in C1 order); ``rng`` is advanced in place."""
    from giap.stats.autocorr import block_len_prime, choose_block_len, morans_i
    from giap.stats.bootstrap import circular_block_bootstrap
    from giap.stats.scan import fwer_p, merge_windows, scan_max, scan_null_max

    isl = cfg.context.islands
    wmin, wmax, alpha = isl.min_tus, isl.max_tus, cfg.stats.alpha.context
    lens: list[int] = []
    hats: list[int] = []
    for s in series:
        if not s.scanned:
            lens.append(0)
            hats.append(0)
            continue
        t = s.detrended.size
        if isinstance(isl.block_len, int):
            hats.append(0)  # no l-hat is estimated when the config fixes the block length
            lens.append(int(min(isl.block_len, t)))
        else:
            l_hat = choose_block_len(s.detrended, s.circular, cfg.stats.permutations, rng,
                                     alpha=0.05, max_len=20)  # fmt: skip
            hats.append(int(l_hat))
            lens.append(block_len_prime(l_hat, t))
    # M_b accumulates with np.fmax over the scanned replicons (NaN only when every term is NaN)
    null = np.full(cfg.stats.bootstrap, np.nan)
    per_replicon: list[np.ndarray] = []
    for s, ell in zip(series, lens, strict=True):
        if not s.scanned:
            per_replicon.append(np.zeros(0))
            continue
        Z = circular_block_bootstrap(s.detrended, ell, cfg.stats.bootstrap, rng)
        m = scan_null_max(Z, wmin, wmax, s.circular)
        per_replicon.append(m)
        null = np.fmax(null, m)
    if not any(x.scanned for x in series):
        null = np.zeros(0)
    windows, pvals, clusters = [], [], []
    for s in series:
        if not s.scanned:
            windows.append(scan_max(np.zeros(0), wmin, wmax, s.circular))
            pvals.append(np.zeros(0))
            clusters.append(())
            continue
        res = scan_max(s.detrended, wmin, wmax, s.circular)
        p = fwer_p(np.abs(res.s), null)
        with np.errstate(invalid="ignore"):
            sig = p <= alpha  # decision 42: NaN is never significant
        windows.append(res)
        pvals.append(p)
        clusters.append(tuple(merge_windows(res, sig, s.detrended.size, s.circular)))
    return OracleScan(
        series=tuple(series),
        l_hat=tuple(hats),
        block_len=tuple(lens),
        morans_i_lag1=tuple(
            morans_i(x.detrended, 1, x.circular) if x.scanned else math.nan for x in series
        ),
        null_max=null,
        windows=tuple(windows),
        replicon_null=tuple(per_replicon),
        p_value=tuple(pvals),
        clusters=tuple(clusters),
        n_perm=int(np.count_nonzero(~np.isnan(null))),
        rng=rng,
    )


def eligible_genes(genes: pd.DataFrame, contig: str, tu_ids: Sequence[str]) -> pd.DataFrame:
    """Eligible genes of the given TUs of ``contig``, in genome order (spec step 7)."""
    g = genes[genes["contig"].astype(str) == contig].copy()
    p = pd.to_numeric(g["p_detect"], errors="coerce")
    elig = ~g["partial"].astype("boolean").fillna(False).to_numpy(bool) & p.notna().to_numpy()
    g = g[elig & g["tu_id"].astype(str).isin(set(tu_ids)).to_numpy()]
    return g.sort_values(["start", "locus_tag"], kind="mergesort").reset_index(drop=True)


def window_tus(s: OracleSeries, start: int, length: int) -> list[int]:
    """TU positions of a window (wrapping on circular replicons)."""
    n = s.detrended.size
    return [(start + k) % n if s.circular else start + k for k in range(length)]


def oracle_frame(genes: pd.DataFrame, scan: OracleScan, cfg: Any, rng: Any) -> pd.DataFrame:
    """The expected ``islands.tsv`` rows (spec steps 7-8), drawing the effect CIs from ``rng``
    in island order."""
    from giap.stats.bootstrap import tu_block_bootstrap

    rows = []
    for ri, (s, res, p, cl) in enumerate(
        zip(scan.series, scan.windows, scan.p_value, scan.clusters, strict=True)
    ):
        for c in cl:
            i = c.best
            pos = window_tus(s, int(res.start[i]), int(res.length[i]))
            tu_ids = [s.tu_ids[j] for j in pos]
            eg = eligible_genes(genes, s.contig, tu_ids)
            det = eg["detected"].astype("boolean").fillna(False).to_numpy(bool)
            a = pd.to_numeric(eg["abundance_log2"], errors="coerce").to_numpy(float)
            have = det & np.isfinite(a)
            wraps = s.circular and int(res.start[i]) + int(res.length[i]) > s.detrended.size
            rows.append(
                {
                    "replicon": ri,
                    "direction": c.direction,
                    "contig": s.contig,
                    "start": int(s.start[pos[0]]),
                    "end": int(s.end[pos[-1]]) + (s.length if wraps else 0),
                    "tu_first": tu_ids[0],
                    "tu_last": tu_ids[-1],
                    "n_tus": len(pos),
                    "n_genes": len(eg),
                    "genes": [str(t) for t in eg["locus_tag"]],
                    "mean_abundance_log2": float(a[have].mean()) if have.any() else math.nan,
                    "detected_frac": (float(det.mean()) if len(eg) else math.nan),
                    "scan_s": float(res.s[i]),
                    "block_len": scan.block_len[ri],
                    "effect": float(np.mean([s.detrended[j] for j in pos])),
                    "p_value": float(p[i]),
                    "n_perm": scan.n_perm,
                    "_pos": pos,
                }
            )
    rows.sort(key=lambda r: (r["replicon"], r["start"], r["direction"] != "high"))
    for n, r in enumerate(rows, start=1):
        r["island_id"] = f"ISL_{n}"
        s = scan.series[r["replicon"]]
        df = pd.DataFrame({"tu_id": [s.tu_ids[j] for j in r["_pos"]],
                           "z": [s.detrended[j] for j in r["_pos"]]})  # fmt: skip
        ci = tu_block_bootstrap(df, "tu_id", lambda d: float(d["z"].mean()),
                                cfg.stats.bootstrap, rng)  # fmt: skip
        r["ci_low"], r["ci_high"] = ci.ci_low, ci.ci_high
        r["flags"] = (FEW_UNITS,) if r["n_tus"] < 10 else ()  # decision 48: no ``uncalibrated``
    return pd.DataFrame(rows)


def oracle_islands(
    genes: pd.DataFrame, tus: pd.DataFrame, replicons: pd.DataFrame, cfg: Any, seed: int
) -> tuple[pd.DataFrame, OracleScan]:
    """The whole oracle pipeline on one random stream (§6.6 ``context_stats.islands``)."""
    from giap.rng import for_step

    rng = for_step(STEP, seed)
    series = oracle_series(genes, tus, replicons, cfg)
    scan = oracle_scan(series, cfg, rng)
    return oracle_frame(genes, scan, cfg, rng), scan


# --------------------------------------------------------------------------- comparisons


def flag_set(cell: Any) -> set[str]:
    """Flags of a cell: tuple / list of names or a ``;``-separated string."""
    return rd.flag_set(cell)


def is_na(x: Any) -> bool:
    return rd.is_na(x)


def same(a: Any, b: Any, tol: float = 1e-9) -> bool:
    """Equality that treats NA as NA and compares floats with a tolerance."""
    if is_na(a) or is_na(b):
        return is_na(a) and is_na(b)
    if isinstance(a, (int, float, np.floating)) and isinstance(b, (int, float, np.floating)):
        return bool(math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=tol))
    if isinstance(a, (list, tuple, np.ndarray)) or isinstance(b, (list, tuple, np.ndarray)):
        return [str(x) for x in a] == [str(x) for x in b]
    return bool(str(a) == str(b))


def assert_matches_oracle(
    got: pd.DataFrame, want: pd.DataFrame, columns: Sequence[str] | None = None
) -> None:
    """Every expected island row, column by column (oracle column names = table names)."""
    cols = list(columns) if columns is not None else [
        "island_id", "direction", "contig", "start", "end", "tu_first", "tu_last", "n_tus",
        "n_genes", "genes", "mean_abundance_log2", "detected_frac", "scan_s", "block_len",
        "effect", "ci_low", "ci_high", "p_value", "n_perm",
    ]  # fmt: skip
    assert len(got) == len(want), f"{len(got)} islands, expected {len(want)}"
    for i in range(len(want)):
        for col in cols:
            a, b = got.iloc[i][col], want.iloc[i][col]
            assert same(a, b), f"row {i} column {col}: {a!r} != {b!r}"
        assert flag_set(got.iloc[i]["flags"]) == set(want.iloc[i]["flags"])


def rows_of(islands: pd.DataFrame, contig: str) -> pd.DataFrame:
    return islands[islands["contig"].astype(str) == contig].reset_index(drop=True)


def counts(genes: pd.DataFrame) -> Mapping[str, int]:
    """Scored TUs per contig (for the ``few_values`` expectations)."""
    return {
        str(c): int(scored_table(genes, str(c))["tu_id"].nunique())
        for c in pd.unique(genes["contig"].astype(str))
    }


# --------------------------------------------------------------------------- driving the step


def build_series(
    genes: pd.DataFrame, tus: pd.DataFrame, replicons: pd.DataFrame, cfg: Any
) -> list[Any]:
    """``TuSeries`` per replicon in C1 order (no randomness)."""
    from giap.context.islands import tu_series
    from giap.context.regions import replicon_index

    return [
        tu_series(
            replicon_index(genes, str(r["contig"])),
            tus,
            contig=str(r["contig"]),
            circular=str(r["topology"]) == "circular",
            length=int(r["length"]),
            span_frac=cfg.context.islands.detrend_span_frac,
        )
        for _, r in replicons.iterrows()
    ]


def step(
    genes: pd.DataFrame, series: Sequence[Any], cfg: Any, rng: Any
) -> tuple[pd.DataFrame, Any]:
    """``scan_replicons`` then ``islands_frame`` on one random stream."""
    from giap.context.islands import islands_frame, scan_replicons

    scan = scan_replicons(list(series), cfg, rng)
    return islands_frame(scan, genes, cfg, rng), scan


def min_p(scan: Any) -> float:
    """Smallest window p over every scanned replicon (NaN never counts, decision 42)."""
    valid = [float(np.nanmin(p)) for p in scan.p_value if p.size and not np.isnan(p).all()]
    return min(valid) if valid else float("inf")


def any_island(
    rng: Any,
    contigs: Sequence[Contig],
    levels: Sequence[float],
    block_len: int | str = "auto",
    bootstrap: int = 199,
    permutations: int = 199,
) -> np.ndarray:
    """Indicator of "at least one island" at every level, from one calibration replicate."""
    cfg = config(bootstrap=bootstrap, permutations=permutations, alpha=max(levels),
                 block_len=block_len)  # fmt: skip
    genes, tus, replicons = frames(contigs)
    out, scan = step(genes, build_series(genes, tus, replicons, cfg), cfg, rng)
    p = min_p(scan)
    assert (len(out) > 0) == (p <= max(levels)), (len(out), p)
    return np.array([float(p <= a) for a in levels])
