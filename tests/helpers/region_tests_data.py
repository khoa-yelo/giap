"""Inputs and oracles for the M4-17 region tests (docs/specs/M4-17-region-tests.md).

Two kinds of inputs, both in the column vocabulary the real producers write (nothing invented):

* hand-built gene tables (``gene_table``) with the §9.3 ``genes.tsv`` columns the region tests read:
  identity (C2: ``locus_tag``, ``feature_type``, ``contig``, ``start``, ``end``, ``strand``,
  ``partial``, ``crosses_origin``), context (C6: ``tu_id``, ``region_ids``), evidence (``detected``,
  ``detection_status``, ``p_detect``), abundance (``abundance_log2``) and the M4-10 protein
  covariates ``length_aa``, ``n_observable_peptides``, ``gravy`` (ADR-0018: HL residualises on
  them);
  C1 ``replicons``, C7 ``regions`` and C8 ``region_genes`` frames with the contract columns;
* simulated gene tables (``sim_tables`` / ``tables_from``): giapsim M4-07 layout (C1, C2, C7, C8;
  C6 via ``layers_frame``, or the whole contract written and re-read with M4-06
  ``load_m4_inputs`` in ``loaded_tables``) + M4-08 expression truth: ``detected = detected_true``,
  ``abundance_log2 = abundance_log2_true`` on detected genes (NaN otherwise: M4-14 fills only
  detected genes), ``p_detect`` = ``p_detect_true`` or the M4-13 cross-fitted model on M4-10
  covariates (``p_source="m4_13"``); covariate columns from M4-10 ``covariate_table``.

Eligibility follows M4-15/M4-16: partial genes have ``detection_status = partial_excluded``,
``detected = False`` and ``p_detect = NaN``. Nothing here imports ``giap.context``.

Oracles: ``naive_o_minus_e_null`` / ``naive_hl_null`` (one Python loop per shift), ``naive_perm_p``
and ``expected_shifts`` (shift sets from the spec text, independent of ``giap.stats``); ADR-0018:
``naive_studentised_null`` (O - E test statistic sum(d - p) / sqrt(sum p (1 - p)) per shift) and
``naive_residualise`` (log2 abundance residualised on the M4-13 transforms of the detectability
covariates among the replicon's detected genes). ADR-0018 amendment 2 (lead decision 15: V
label-free, superseding amendment 1's per-TU V): ``naive_lag_autocov`` (lag autocovariance
gamma(h) per replicon), ``naive_lag_hl_variance`` (per-shift V) and ``naive_studentised_hl_null``
(T = HL / sqrt(V)).
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

REGION_TYPES = ("prophage", "bgc", "crispr", "amr")
STATS = ("o_minus_e", "hl_shift")
METHOD = {"o_minus_e": "circular_shift_o_minus_e", "hl_shift": "circular_shift_hodges_lehmann"}
#: ``RegionTest`` fields (spec "Public interface"), then what ``test_all_regions`` adds
RESULT_FIELDS = (
    "region_id", "region_type", "contig", "stat", "observed", "expected", "effect", "ci_low",
    "ci_high", "p_value", "n_perm", "n_genes", "n_eligible", "flags",
)  # fmt: skip
TABLE_COLUMNS = (*RESULT_FIELDS, "q_value", "method", "null_model")
GENE_COLUMNS = (
    "locus_tag", "feature_type", "contig", "start", "end", "strand", "partial", "crosses_origin",
    "tu_id", "region_ids", "detected", "detection_status", "p_detect", "abundance_log2",
    "length_aa", "n_observable_peptides", "gravy",
)  # fmt: skip
#: ADR-0018 decision 1: the detectability covariates HL residualises on (§9.3 genes.tsv names)
COVARIATE_COLUMNS = ("length_aa", "n_observable_peptides", "gravy")
#: hand-built genes span 300 nt (99 codons + stop): constant covariates unless a Contig sets them
DEFAULT_COVARIATES = (99, 5, -0.2)
REGION_COLUMNS = ("region_id", "region_type", "contig", "start", "end", "source_tool", "score")
REGION_GENE_COLUMNS = ("region_id", "locus_tag", "overlap_frac", "subset", "core_biosynthetic")
FEW_SHIFTS_K = 19  # lead decision 5 / Q3: K < 19 -> the smallest p 1/(1+K) exceeds 0.05


def config(
    bootstrap: int = 50, permutations: int = 1000, exhaustive: int = 20000, **kw: Any
) -> Any:
    """A ``GiapConfig`` with ``stats.bootstrap``, ``stats.permutations``,
    ``context.exhaustive_shift_max_genes`` and optional ``abundance_missing`` set."""
    from giap.config.model import GiapConfig

    raw: dict[str, Any] = {
        "stats": {"bootstrap": bootstrap, "permutations": permutations},
        "context": {"exhaustive_shift_max_genes": exhaustive},
    }
    if "abundance_missing" in kw:
        raw["abundance"] = {"missing": kw.pop("abundance_missing")}
    if kw:
        raise TypeError(f"unknown config keys {sorted(kw)}")
    return GiapConfig.model_validate(raw)


def flag_set(cell: Any) -> set[str]:
    """Flags of one result cell: a tuple/list of names or a ``;``-separated string."""
    if cell is None:
        return set()
    if isinstance(cell, str):
        return {f for f in cell.split(";") if f}
    if isinstance(cell, float) and math.isnan(cell):
        return set()
    return {str(f) for f in cell}


def is_na(x: Any) -> bool:
    return x is None or x is pd.NA or (isinstance(x, float) and math.isnan(x))


# --------------------------------------------------------------------------- hand-built tables


@dataclass
class Contig:
    """One replicon for ``gene_table``: per gene detection d, p_detect p, log2 abundance a.

    ``partial`` marks partial genes (linear contig ends); ``tu`` gives TU labels (default: one TU
    per gene). ``a`` defaults to 10 + gene index on detected genes, NaN elsewhere. ``length_aa``,
    ``n_obs`` (``n_observable_peptides``) and ``gravy`` default to ``DEFAULT_COVARIATES`` (constant:
    residualising on them only centres the abundances, which leaves every HL difference unchanged).
    """

    name: str
    d: Sequence[int]
    p: Sequence[float]
    a: Sequence[float] | None = None
    circular: bool = True
    partial: Sequence[bool] | None = None
    tu: Sequence[str] | None = None
    kind: str | None = None
    length_aa: Sequence[int] | None = None
    n_obs: Sequence[int] | None = None
    gravy: Sequence[float] | None = None


def gene_table(contigs: Sequence[Contig]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(genes, replicons) for the given contigs; locus tags ``<contig>_g<i:05d>`` in genome order.

    Gene ``i`` spans [400 i + 100, 400 i + 400); replicon length = 400 n + 200.
    """
    rows: list[dict[str, Any]] = []
    reps = []
    for c in contigs:
        n = len(c.d)
        if len(c.p) != n or (c.a is not None and len(c.a) != n):
            raise ValueError(f"contig {c.name}: d, p, a lengths differ")
        part = [False] * n if c.partial is None else [bool(x) for x in c.partial]
        tus = [f"TU_{c.name}_{i}" for i in range(n)] if c.tu is None else list(c.tu)
        L0, N0, G0 = DEFAULT_COVARIATES
        la_ = [L0] * n if c.length_aa is None else list(c.length_aa)
        no_ = [N0] * n if c.n_obs is None else list(c.n_obs)
        gr_ = [G0] * n if c.gravy is None else list(c.gravy)
        if not len(la_) == len(no_) == len(gr_) == n:
            raise ValueError(f"contig {c.name}: covariate lengths differ from d")
        for i in range(n):
            det = bool(c.d[i]) and not part[i]
            if c.a is None:
                a = 10.0 + i if det else math.nan
            else:
                a = float(c.a[i]) if det else math.nan
            rows.append(
                {
                    "locus_tag": f"{c.name}_g{i:05d}",
                    "feature_type": "cds",
                    "contig": c.name,
                    "start": 400 * i + 100,
                    "end": 400 * i + 400,
                    "strand": "+",
                    "partial": part[i],
                    "crosses_origin": False,
                    "tu_id": tus[i],
                    "region_ids": [],
                    "detected": det,
                    "detection_status": (
                        "partial_excluded" if part[i] else ("detected" if det else "not_detected")
                    ),
                    "p_detect": math.nan if part[i] else float(c.p[i]),
                    "abundance_log2": a,
                    "length_aa": int(la_[i]),
                    "n_observable_peptides": int(no_[i]),
                    "gravy": float(gr_[i]),
                }
            )
        if c.circular:
            kind = c.kind or ("chromosome" if not reps else "plasmid")
        else:
            kind = c.kind or "unknown"
        reps.append((c.name, 400 * n + 200, "circular" if c.circular else "linear", kind))
    genes = pd.DataFrame(rows, columns=list(GENE_COLUMNS))
    genes = genes.astype({"start": "int64", "end": "int64", "partial": bool, "detected": bool,
                          "length_aa": "int64", "n_observable_peptides": "int64",
                          "gravy": "float64"})  # fmt: skip
    replicons = pd.DataFrame(reps, columns=["contig", "length", "topology", "kind"])
    return genes, replicons.astype({"length": "int64"})


def region_frames(
    genes: pd.DataFrame, members: Mapping[str, tuple[str, Sequence[str]]]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """C7 regions and C8 region_genes for ``{region_id: (region_type, locus_tags)}``; also returns
    ``genes`` with ``region_ids`` filled (sorted, as M2-20 writes them). C7 ``contig`` is the first
    member's contig; AMR rows have NA coordinates, the others the members' span."""
    g = genes.set_index("locus_tag")
    reg_rows, rg_rows = [], []
    of: dict[str, list[str]] = {}
    for rid, (typ, tags) in members.items():
        tags = list(tags)
        contig = str(g.loc[tags[0], "contig"]) if tags else "chr1"
        if typ == "amr" or not tags:
            s, e = None, None
        else:
            s, e = int(g.loc[tags, "start"].min()), int(g.loc[tags, "end"].max())
        reg_rows.append((rid, typ, contig, s, e, "test", None))
        for t in tags:
            rg_rows.append((rid, t, 1.0, None, False))
            of.setdefault(t, []).append(rid)
    regions = pd.DataFrame(reg_rows, columns=list(REGION_COLUMNS)).astype(
        {"start": "Int64", "end": "Int64", "score": "float64"}
    )
    rg = pd.DataFrame(rg_rows, columns=list(REGION_GENE_COLUMNS)).astype(
        {"overlap_frac": "float64", "core_biosynthetic": bool}
    )
    out = genes.copy()
    out["region_ids"] = [sorted(of.get(t, [])) for t in out["locus_tag"]]
    return out, regions, rg


def tags_of(genes: pd.DataFrame, contig: str, idx: Iterable[int]) -> list[str]:
    """Locus tags of the genes at positions ``idx`` (all genes, partial included) of ``contig``."""
    sub = genes[genes["contig"] == contig].reset_index(drop=True)
    return [str(sub.loc[i, "locus_tag"]) for i in idx]


def hot_contig(name: str, n_hot: int, n: int = 100, width: int = 3) -> Contig:
    """Circular contig of ``n`` genes whose region I = {0 .. width-1} has d - p = +0.5 per gene and
    ``n_hot`` disjoint copies of that pattern far from I and from each other. Exactly ``n_hot`` of
    the n - 1 shifts tie the observed O - E, all others are
    smaller, so the exhaustive two-sided p is (1 + n_hot) / n (spec example 8).

    Every other gene has d = 0, p = 0.01 (not 0: the ADR-0018 studentised statistic divides by
    sum p (1 - p) of the shifted set, which must stay positive); the p-values hold for both the
    plain sum and the studentised statistic."""
    d = [0] * n
    p = [0.01] * n
    starts = [0] + [20 * (j + 1) for j in range(n_hot)]
    if starts[-1] + width + width > n:
        raise ValueError("too many hot copies for n")
    for s in starts:
        for i in range(s, s + width):
            d[i], p[i] = 1, 0.5
    return Contig(name, d, p)


# --------------------------------------------------------------------------- oracles


def expected_shifts(n: int, idx: Sequence[int], circular: bool) -> list[int]:
    """Spec §7.3.5(b): k = 1..N-1; on linear contigs only k keeping (I + k) mod N inside [0, N)
    without wrapping (the shifted set is a translation of I)."""
    out = []
    lo, hi = min(idx), max(idx)
    for k in range(1, n):
        if circular or hi + k < n or lo + k >= n:
            out.append(k)
    return out


def naive_o_minus_e_null(r: np.ndarray, idx: Sequence[int], shifts: Sequence[int]) -> np.ndarray:
    n = len(r)
    return np.array([sum(float(r[(i + k) % n]) for i in idx) for k in shifts], dtype=float)


def naive_studentised_null(
    r: np.ndarray, v: np.ndarray, idx: Sequence[int], shifts: Sequence[int]
) -> np.ndarray:
    """ADR-0018 decision 2: T_k = sum r / sqrt(sum v) over (I + k) mod N, r = d - p,
    v = p (1 - p)."""
    n = len(r)
    out = []
    for k in shifts:
        s = sum(float(r[(i + k) % n]) for i in idx)
        w = sum(float(v[(i + k) % n]) for i in idx)
        out.append(s / math.sqrt(w))
    return np.array(out, dtype=float)


def studentised(r: np.ndarray, v: np.ndarray, idx: Sequence[int]) -> float:
    return float(sum(float(r[i]) for i in idx) / math.sqrt(sum(float(v[i]) for i in idx)))


def naive_residualise(
    a: np.ndarray,
    det: np.ndarray,
    length_aa: Sequence[float],
    n_obs: Sequence[float],
    gravy: Sequence[float],
) -> np.ndarray:
    """ADR-0018 decision 1 (test-author reading): among the replicon's detected genes with a finite
    log2 abundance, ordinary least squares of ``a`` on an intercept and the M4-13
    ``transform_covariates`` of the detectability covariates, log(length_aa),
    log1p(n_observable_peptides) and gravy (a non-finite transform -> the median of that column over
    the fitted genes); returns the residuals there and NaN elsewhere. Rank-deficient designs (e.g.
    constant covariates) use the minimum-norm solution, so constant covariates only centre ``a``."""
    av = np.asarray(a, dtype=float)
    ok = np.asarray(det, dtype=bool) & np.isfinite(av)
    out = np.full(av.size, np.nan)
    if not ok.any():
        return out
    with np.errstate(divide="ignore", invalid="ignore"):
        cols = [
            np.log(np.asarray(length_aa, dtype=float)),
            np.log1p(np.asarray(n_obs, dtype=float)),
            np.asarray(gravy, dtype=float),
        ]
    X = [np.ones(int(ok.sum()))]
    for c in cols:
        x = c[ok].copy()
        bad = ~np.isfinite(x)
        if bad.any():
            x[bad] = np.median(x[~bad]) if (~bad).any() else 0.0
        X.append(x)
    design = np.column_stack(X)
    beta, *_ = np.linalg.lstsq(design, av[ok], rcond=None)
    out[ok] = av[ok] - design @ beta
    return out


def naive_hl(x: Sequence[float], y: Sequence[float]) -> float:
    """Median of all pairwise differences x_i - y_j (NaN if either side is empty)."""
    if len(x) == 0 or len(y) == 0:
        return math.nan
    return float(np.median([xi - yj for xi in x for yj in y]))


def naive_hl_null(
    a: np.ndarray, det: np.ndarray, idx: Sequence[int], shifts: Sequence[int]
) -> np.ndarray:
    n = len(a)
    out = []
    for k in shifts:
        inside = {(i + k) % n for i in idx}
        ok = [j for j in range(n) if bool(det[j]) and np.isfinite(a[j])]
        x = [float(a[j]) for j in ok if j in inside]
        y = [float(a[j]) for j in ok if j not in inside]
        out.append(naive_hl(x, y))
    return np.array(out, dtype=float)


LAG_MAX = 30  # ADR-0018 amendment 2: gamma(h) for h <= 30 (statistician's v3, lagv.py H)


def _finite(v: Any) -> bool:
    return not (isinstance(v, float) and not math.isfinite(v))


def naive_lag_autocov(a: Sequence[Any], ok: Any, max_lag: int = LAG_MAX) -> list[Any]:
    """ADR-0018 amendment 2 (docs/methods/region-tests.md, review round 3, item 4): the replicon's
    autocovariance of the residualised detected abundances at circular gene-index lag
    h = 0 ... max_lag. x_i = a_i - mean(a over ok) on ``ok`` genes (flag set and value finite), 0
    elsewhere; gamma(h) = sum_i x_i x_{(i+h) mod N} / #{i : ok_i and ok_{(i+h) mod N}} (0 if that
    count is 0); truncated at its first non-positive value: gamma(h') = 0 for every h' >= the first
    h >= 1 with gamma(h) <= 0. Exact for Fraction inputs. No ok gene: every gamma 0."""
    n = len(a)
    okl = [bool(ok[i]) and _finite(a[i]) for i in range(n)]
    vals = [a[i] for i in range(n) if okl[i]]
    if not vals:
        return [0] * (max_lag + 1)
    mu = sum(vals) / len(vals)
    x = [a[i] - mu if okl[i] else 0 for i in range(n)]
    g: list[Any] = []
    for h in range(max_lag + 1):
        num = sum(x[i] * x[(i + h) % n] for i in range(n))
        den = sum(1 for i in range(n) if okl[i] and okl[(i + h) % n])
        g.append(num / den if den else 0)
    for h in range(1, max_lag + 1):
        if g[h] <= 0:
            g[h:] = [0] * (max_lag + 1 - h)
            break
    return g


def naive_lag_hl_variance(
    ok: Any, idx: Sequence[int], shifts: Sequence[int], gamma: Sequence[Any]
) -> list[Any]:
    """V_k = Var(mean_in - mean_out) = S_in,in / m^2 + S_out,out / m'^2 - 2 S_in,out / (m m'),
    S_A,B = sum over ok i in A and ok j in B of gamma(d(i, j)), d the circular gene-index distance
    min(|i - j|, N - |i - j|) (gamma = 0 beyond its length), for the set (I + k) mod N; m / m' the
    ok genes inside / outside; NaN if m or m' is 0. Exact for Fraction gamma."""
    n = len(ok)

    def gam(i: int, j: int) -> Any:
        d = abs(i - j)
        d = min(d, n - d)
        return gamma[d] if d < len(gamma) else 0

    out: list[Any] = []
    for k in shifts:
        inside = {(i + k) % n for i in idx}
        a_in = [j for j in range(n) if bool(ok[j]) and j in inside]
        a_out = [j for j in range(n) if bool(ok[j]) and j not in inside]
        m, mo = len(a_in), len(a_out)
        if m == 0 or mo == 0:
            out.append(math.nan)
            continue
        s_ii = sum(gam(i, j) for i in a_in for j in a_in)
        s_oo = sum(gam(i, j) for i in a_out for j in a_out)
        s_io = sum(gam(i, j) for i in a_in for j in a_out)
        out.append(s_ii / m**2 + s_oo / mo**2 - 2 * s_io / (m * mo))
    return out


def _lag_v(a: np.ndarray, ok: np.ndarray, idx: Sequence[int], shifts: Sequence[int]) -> list[Any]:
    okv = np.asarray(ok, dtype=bool) & np.isfinite(np.asarray(a, dtype=float))
    g = naive_lag_autocov([float(x) for x in a], okv)
    return naive_lag_hl_variance(okv, idx, shifts, g)


def naive_studentised_hl_null(
    a: np.ndarray, ok: np.ndarray, idx: Sequence[int], shifts: Sequence[int]
) -> np.ndarray:
    """ADR-0018 amendment 2: T_k = HL_k / sqrt(V_k) for every shift k (k = 0 gives the observed
    T); gamma estimated once from (a, ok) (``naive_lag_autocov``), V label-free
    (``naive_lag_hl_variance``); NaN where a side is empty; V = 0 with HL = 0 gives T = 0 (V <= 0
    with HL != 0 gives +-inf; comparisons with the implementation skip such draws,
    ``naive_v_positive``)."""
    okv = np.asarray(ok, dtype=bool) & np.isfinite(np.asarray(a, dtype=float))
    hl = naive_hl_null(a, okv, idx, shifts)
    v = _lag_v(a, okv, idx, shifts)
    out = []
    for h, w in zip(hl, v, strict=True):
        if math.isnan(h) or math.isnan(w):
            out.append(math.nan)
        elif w > 0:
            out.append(h / math.sqrt(w))
        else:
            out.append(0.0 if h == 0 else math.copysign(math.inf, h))
    return np.array(out, dtype=float)


def naive_v_positive(
    a: np.ndarray, ok: np.ndarray, idx: Sequence[int], shifts: Sequence[int]
) -> bool:
    """Every finite V_k > 0 (relative to gamma(0)). The truncated gamma need not be positive
    definite; V <= 0 is not pinned by the spec, so oracle comparisons skip such inputs."""
    okv = np.asarray(ok, dtype=bool) & np.isfinite(np.asarray(a, dtype=float))
    g0 = naive_lag_autocov([float(x) for x in a], okv, 0)[0]
    v = _lag_v(a, okv, idx, shifts)
    return all(math.isnan(w) or w > 1e-9 * max(g0, 1e-300) for w in v)


def naive_perm_p(obs: float, null: np.ndarray) -> float:
    """Two-sided p centred on the pooled mean (M4-03), NaN nulls dropped, ties within 1e-9."""
    z = np.asarray(null, dtype=float)
    z = z[~np.isnan(z)]
    if z.size == 0 or math.isnan(obs):
        return math.nan
    c = (z.sum() + obs) / (z.size + 1)
    d = abs(obs - c)
    tol = 1e-9 * max(1.0, d, abs(c))
    return (1 + int(np.count_nonzero(np.abs(z - c) >= d - tol))) / (1 + z.size)


# --------------------------------------------------------------------------- simulated tables


@dataclass
class SimTables:
    genes: pd.DataFrame
    regions: pd.DataFrame
    region_genes: pd.DataFrame
    replicons: pd.DataFrame
    layout: Any
    expression: Any


#: Calibration layout (spec acceptance "Calibration"): default regions + a 60-gene plasmid with its
#: own scattered AMR set (6 genes; size not fixed by the spec), origin crossing on (M4-07 lead
#: decision 2).
PLASMID_GENES = 60
PLASMID_AMR = 6


def calibration_layout_config(n_genes: int, **kw: Any) -> Any:
    """M4-07 ``LayoutConfig``: ``n_genes`` chromosome genes, default regions, plasmid ``p1`` of 60
    genes with its own AMR set (``RegionSpec.contig``, added by M4-17), origin crossing."""
    import giapsim.layout as la

    regions = (*la.DEFAULT_REGIONS, la.RegionSpec("amr", PLASMID_AMR, scattered=True, contig="p1"))
    base: dict[str, Any] = {
        "n_genes_chromosome": n_genes,
        "plasmids": (PLASMID_GENES,),
        "regions": regions,
        "origin_crossing": True,
    }
    base.update(kw)
    return la.LayoutConfig(**base)


def simulate(
    rng: np.random.Generator, layout_cfg: Any, expr_kw: Mapping[str, Any] | None = None
) -> tuple[Any, Any]:
    """(layout, expression) from one generator; detection ``covariates_only`` unless overridden."""
    import giapsim.expression as ex
    import giapsim.layout as la

    lay = la.simulate_layout(layout_cfg, rng)
    ecfg = ex.ExpressionConfig(**{"detection": "covariates_only", **dict(expr_kw or {})})
    return lay, ex.simulate_expression(lay, ecfg, rng)


def sim_tables(
    rng: np.random.Generator,
    layout_cfg: Any,
    expr_kw: Mapping[str, Any] | None = None,
    p_source: str = "true",
) -> SimTables:
    """Gene table, C1/C7/C8 from one giapsim layout + expression draw (see ``tables_from``)."""
    lay, expr = simulate(rng, layout_cfg, expr_kw)
    return tables_from(lay, expr, p_source, rng)


def tables_from(
    lay: Any,
    expr: Any,
    p_source: str | np.ndarray = "true",
    rng: np.random.Generator | None = None,
    genes_c2: pd.DataFrame | None = None,
    layers: pd.DataFrame | None = None,
) -> SimTables:
    """Gene table for a layout and expression truth.

    ``p_source``: ``"true"`` (``p_detect_true``), ``"m4_13"`` (cross-fitted M4-13 model on the
    M4-10 covariates of every eligible gene, seeded from ``rng``), ``"constant"`` (the eligible
    detected share; the naive model of the attack-2 teeth test) or an explicit array in layout gene
    order. ``genes_c2`` / ``layers`` default to the layout's C2 frame and ``layers_frame``.
    """
    from giap.config.model import GiapConfig
    from giap.genetable.covariates import covariate_table
    from giapsim.contract_tables import layers_frame

    e = expr.genes.set_index("locus_tag")
    g = lay.genes if genes_c2 is None else genes_c2
    lf = (layers_frame(lay) if layers is None else layers).set_index("locus_tag")
    tags = g["locus_tag"].astype(str)
    part = g["partial"].astype(bool).to_numpy()
    det = e.loc[tags, "detected_true"].astype(bool).to_numpy() & ~part
    ab = np.where(det, e.loc[tags, "abundance_log2_true"].to_numpy(float), np.nan)
    if isinstance(p_source, np.ndarray):
        p = np.asarray(p_source, dtype=float)
    elif p_source == "true":
        p = e.loc[tags, "p_detect_true"].to_numpy(float)
    elif p_source == "m4_13":
        if rng is None:
            raise ValueError("p_source='m4_13' needs rng")
        p = m4_13_p_detect(lay, det, part, rng)
    elif p_source == "constant":
        p = np.full(len(tags), float(det[~part].mean()))
    else:
        raise ValueError(f"unknown p_source {p_source!r}")
    p = np.where(part, np.nan, p)
    cov = covariate_table(g, GiapConfig())  # M4-10, the genes.tsv protein columns (§9.3)
    genes = pd.DataFrame(
        {
            "locus_tag": tags.to_numpy(),
            "feature_type": g["feature_type"].to_numpy(),
            "contig": g["contig"].to_numpy(),
            "start": g["start"].to_numpy(np.int64),
            "end": g["end"].to_numpy(np.int64),
            "strand": g["strand"].to_numpy(),
            "partial": part,
            "crosses_origin": g["crosses_origin"].astype(bool).to_numpy(),
            "tu_id": lf.loc[tags, "tu_id"].to_numpy(),
            "region_ids": [list(x) for x in lf.loc[tags, "region_ids"]],
            "detected": det,
            "detection_status": np.where(
                part, "partial_excluded", np.where(det, "detected", "not_detected")
            ),
            "p_detect": p,
            "abundance_log2": ab,
            "length_aa": cov["length_aa"].to_numpy(np.int64),
            "n_observable_peptides": cov["n_observable_peptides"].to_numpy(np.int64),
            "gravy": cov["gravy"].to_numpy(float),
        },
        columns=list(GENE_COLUMNS),
    )
    return SimTables(genes, lay.regions, lay.region_genes, lay.replicons, lay, expr)


def loaded_tables(lay: Any, expr: Any, work: Path) -> tuple[SimTables, Any]:
    """Write the layout's C1-C8 with giapsim ``write_contract_tables``, read them back with M4-06
    ``load_m4_inputs`` and build the gene table from the loaded C2/C6 (``p = p_detect_true``).
    Returns the tables (loaded C1/C7/C8) and the load issues."""
    from giap.io.work_tables import default_paths, load_m4_inputs
    from giapsim.contract_tables import write_contract_tables

    write_contract_tables(lay, work)
    inputs, issues = load_m4_inputs(default_paths(work))
    t = tables_from(lay, expr, "true", genes_c2=inputs.genes, layers=inputs.layers)
    return SimTables(t.genes, inputs.regions, inputs.region_genes, inputs.replicons, lay, expr), (
        issues
    )


def m4_13_p_detect(
    lay: Any, det: np.ndarray, part: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """Cross-fitted M4-13 ``p_detect`` (PR #83 interface) on M4-10 ``covariate_table``."""
    from giap.config.model import DetectabilityConfig, GiapConfig
    from giap.genetable.covariates import covariate_table
    from giap.genetable.detectability import fit_detectability  # M4-13 (PR #83)
    from giap.rng import for_step

    tags = lay.genes["locus_tag"].astype(str)
    cov = covariate_table(lay.genes, GiapConfig()).set_index("locus_tag", drop=False)
    d = pd.array(det, dtype="boolean")
    d[part] = pd.NA
    detected = pd.Series(d, index=tags.to_numpy(), name="detected")
    tu = pd.Series(lay.truth["tu_id"].to_numpy(), index=lay.truth["locus_tag"].astype(str))
    seed = int(rng.integers(0, 2**31 - 1))
    res = fit_detectability(
        detected, cov, tu, DetectabilityConfig(), for_step("detectability", seed)
    )
    return np.asarray(res.p_detect.reindex(tags.to_numpy()).to_numpy(float))


def run_region_tests(tables: SimTables, cfg: Any, rng: np.random.Generator) -> pd.DataFrame:
    """``giap.context.regions.test_all_regions`` on simulated tables (imported lazily)."""
    from giap.context import regions as R

    return R.test_all_regions(
        tables.genes, tables.regions, tables.region_genes, tables.replicons, cfg, rng
    )


def p_or_one(p: object) -> float:
    """A region test that cannot be computed (p NA) does not reject: count it as p = 1."""
    x = float(p)  # type: ignore[arg-type]
    return 1.0 if math.isnan(x) else x


def pick(res: pd.DataFrame, region_id: str, stat: str, contig: str | None = None) -> pd.Series:
    """The one result row of (region_id, stat[, contig])."""
    m = (res["region_id"] == region_id) & (res["stat"] == stat)
    if contig is not None:
        m &= res["contig"] == contig
    rows = res[m]
    if len(rows) != 1:
        raise AssertionError(f"{len(rows)} rows for {region_id}/{stat}/{contig}")
    return rows.iloc[0]
