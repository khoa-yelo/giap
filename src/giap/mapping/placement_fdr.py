"""Placement FDR (G3): strata, pi_wrong, q-values and ``mapping_fdr`` rows (issue M3-10).

GIAP_SPEC.md §7.2.5 steps 1-6 as amended by ADR-0010 (docs/methods/placement-fdr.md, "pi_wrong in
the stratum FDR"): within stratum m and decoy type, one list ordered by the ranking score E,
FDR_m(e) = (W_m(e) + c) / T_m(e), W = sum of the decoy hits' own pi_wrong with E <= e, c = max
pi_wrong over the searched peptides of that decoy type (1 for m >= 1). With constant pi this is
pi (D(e) + 1) / T(e). q-values use the M4-01 monotonisation (``giap.stats.tda.monotone_q``).

ADR-0019 supersedes §7.2.5 step 4 for database input in stratum 0 (docs/methods/placement-fdr.md,
"Cross-strain input and run collapse"): the pi_wrong of a database peptide is the estimate of its
length class, min(1, U95(k, n) / (1 - r)) over the stratum-0 database population, and c is the
largest class estimate. PEP and q-values are not read for it; the mean PEP only feeds the
``input_prefiltered`` run warning.

Pure functions; no file I/O.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.stats import beta, binom

from giap.config import GiapConfig
from giap.mapping.chance import ResidueFreqs, expected_hits
from giap.stats.tda import monotone_q

__all__ = [
    "DECOY_TYPE_OF",
    "PREFILTERED_FLAG",
    "STRATUM0_CLASS_COLUMNS",
    "PeptideCosts",
    "PlacementTables",
    "best_costs",
    "clopper_pearson_upper",
    "input_prefiltered",
    "is_low_complexity",
    "length_class_labels",
    "placement",
    "placement_tables",
    "stratum0_pi_classes",
    "stratum_qvalues",
]

STRATA: tuple[int, ...] = (0, 1, 2)
DECOY_TYPE_OF: dict[str, str] = {"database": "reverse", "denovo": "shuffle"}
FLAG_SEP = ";"
LOW_COMPLEXITY_MAX_DISTINCT = 3
CP_LEVEL = 0.95  # one-sided Clopper-Pearson level of the class estimate (ADR-0019 decision 1)
PREFILTERED_QUANTILE = 0.001  # binomial quantile of the run warning (ADR-0019 decision 4)
PREFILTERED_FLAG = "input_prefiltered"
STRATUM0_CLASS_COLUMNS: tuple[str, ...] = (
    "length_class",
    "n",
    "k",
    "r",
    "upper",
    "pi_wrong",
)

PER_PEPTIDE_COLUMNS: tuple[str, ...] = (
    "peptide",
    "source",
    "decoy_type",
    "length",
    "m_target",
    "e_target",
    "placement_q",
    "accepted",
    "flags",
)
MAPPING_FDR_COLUMNS: tuple[str, ...] = (
    "section",
    "decoy_type",
    "n_diff",
    "length_class",
    "targets",
    "decoys",
    "pi_wrong",
    "fdr_stratum",
    "n_accepted",
)
_NUMERIC = ("length", "m_target", "m_decoy", "e_target", "e_decoy", "pep", "q_value", "pi_denovo")


@dataclass(frozen=True, slots=True)
class PeptideCosts:
    """One mapped query peptide at its best target and decoy cost (§7.2.5 step 1)."""

    peptide: str
    length: int
    source: Literal["database", "denovo"]
    search_db: Literal["same", "other", "auto"]
    m_target: int | None
    m_decoy: int | None
    e_target: float | None
    e_decoy: float | None
    pep: float | None
    q_value: float | None
    pi_denovo: float | None
    is_contaminant: bool | None = None  # ADR-0019 / M3-25 decision 7: out of the pi-hat population


_COST_COLUMNS: tuple[str, ...] = tuple(PeptideCosts.__dataclass_fields__)


@dataclass(frozen=True, slots=True)
class PlacementTables:
    """Everything ``placement`` produces: the two §9.3 tables, the ADR-0019 stratum-0 class table
    and the run flags the calling step logs (``input_prefiltered``)."""

    per_peptide: pd.DataFrame
    mapping_fdr: pd.DataFrame
    stratum0_classes: pd.DataFrame
    run_flags: tuple[str, ...]


def _vec(x: npt.ArrayLike, name: str) -> npt.NDArray[np.float64]:
    a = np.asarray(x, dtype=np.float64)
    if a.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return a


def stratum_qvalues(
    e_target: npt.ArrayLike,
    e_decoy: npt.ArrayLike,
    pi_decoy: npt.ArrayLike,
    pi_plus: float,
) -> npt.NDArray[np.float64]:
    """Placement q-values of one stratum list, one per target in input order.

    T(e) = #targets with E <= e, W(e) = sum of ``pi_decoy`` over decoys with E <= e (a decoy tied
    with a target counts), FDR(e) = (W(e) + ``pi_plus``) / T(e) at every target E;
    q_i = min over e >= E_i of FDR(e), clipped to [0, 1]. NaN targets get NaN; NaN decoys are
    ignored.

    Raises:
        ValueError: on non-1-D input, mismatched decoy/weight sizes, weights outside [0, 1] or
            ``pi_plus`` outside [0, 1].
    """
    et = _vec(e_target, "e_target")
    ed = _vec(e_decoy, "e_decoy")
    w = _vec(pi_decoy, "pi_decoy")
    if ed.shape != w.shape:
        raise ValueError("e_decoy and pi_decoy must have the same length")
    if not 0.0 <= pi_plus <= 1.0:
        raise ValueError(f"pi_plus must lie in [0, 1], got {pi_plus}")
    if w.size and (np.any(np.isnan(w)) or w.min() < 0.0 or w.max() > 1.0):
        raise ValueError("pi_decoy must lie in [0, 1]")
    q = np.full(et.shape, np.nan, dtype=np.float64)
    ok = ~np.isnan(et)
    s = et[ok]
    if s.size == 0:
        return q
    keep = ~np.isnan(ed)
    d_order = np.argsort(ed[keep], kind="stable")
    d_sorted = ed[keep][d_order]
    cum_w = np.concatenate(([0.0], np.cumsum(w[keep][d_order])))
    order = np.argsort(s, kind="stable")
    s_sorted = s[order]
    n_t = np.searchsorted(s_sorted, s_sorted, side="right").astype(np.float64)
    big_w = cum_w[np.searchsorted(d_sorted, s_sorted, side="right")]
    q_sorted = monotone_q((big_w + float(pi_plus)) / n_t)
    out = np.empty(s.shape, dtype=np.float64)
    out[order] = q_sorted
    q[ok] = out
    return q


def clopper_pearson_upper(k: int, n: int, level: float = CP_LEVEL) -> float:
    """One-sided Clopper-Pearson upper bound of the binomial share k/n at ``level``.

    The smallest p with P(X <= k; n, p) = 1 - ``level``, i.e. the ``level`` quantile of
    Beta(k + 1, n - k); 1 when k = n.

    Parameters
    ----------
    k, n:
        Successes and trials, 0 <= k <= n, n >= 1.
    level:
        Confidence level, in (0, 1).

    Returns
    -------
    float
        The upper bound, in [k/n, 1].

    Raises
    ------
    ValueError
        On n < 1, k outside [0, n] or ``level`` outside (0, 1).
    """
    if n < 1:
        raise ValueError(f"n must be at least 1, got {n}")
    if not 0 <= k <= n:
        raise ValueError(f"k must lie in [0, {n}], got {k}")
    if not 0.0 < level < 1.0:
        raise ValueError(f"level must lie in (0, 1), got {level}")
    if k >= n:
        return 1.0
    return float(beta.ppf(level, k + 1, n - k))


def stratum0_pi_classes(costs: pd.DataFrame, classes: Sequence[Sequence[int]]) -> pd.DataFrame:
    """ADR-0019 decision 1: the stratum-0 database pi_wrong per length class.

    ``costs`` holds the stratum-0 database population P (the caller drops de novo rows,
    ``same_db`` exact peptides and contaminants) at its best target and decoy cost. Per class l of
    ``classes`` (a peptide belongs to the first class holding its length): n_l peptides, k_l of
    them without an exact target hit (``m_target`` != 0 or missing), r_l the share with an exact
    decoy hit, U_l the one-sided 95 % Clopper-Pearson upper bound of k_l/n_l and
    pi_l = 1 if r_l = 1 else min(1, U_l / (1 - r_l)).

    Parameters
    ----------
    costs:
        Population rows with ``length``, ``m_target`` and ``m_decoy`` (the ``best_costs`` schema;
        a missing cost means no hit at that search).
    classes:
        ``mapping.length_classes`` as (lo, hi) pairs, inclusive.

    Returns
    -------
    pandas.DataFrame
        One row per class in the order of ``classes``, columns ``STRATUM0_CLASS_COLUMNS``. An
        empty class has n = k = 0 and NaN ``r``, ``upper`` and ``pi_wrong``.

    Raises
    ------
    ValueError
        When ``classes`` is empty, a class is not an inclusive ``[lo, hi]`` pair, or ``costs``
        lacks ``length``, ``m_target`` or ``m_decoy``.
    """
    if not len(classes):
        raise ValueError("stratum0_pi_classes needs at least one length class")
    for item in classes:
        if len(item) != 2 or item[0] > item[1]:
            raise ValueError(f"length class {list(item)} is not an inclusive [lo, hi] range")
    missing = [c for c in ("length", "m_target", "m_decoy") if c not in costs.columns]
    if missing:
        raise ValueError(f"costs lacks the column(s) {missing} of the best_costs schema")
    bounds = [[int(lo), int(hi)] for lo, hi in classes]
    labels = length_class_labels(bounds)
    length = pd.to_numeric(costs["length"], errors="coerce").to_numpy(dtype=np.float64)
    exact_t = pd.to_numeric(costs["m_target"], errors="coerce").to_numpy(dtype=np.float64) == 0.0
    exact_d = pd.to_numeric(costs["m_decoy"], errors="coerce").to_numpy(dtype=np.float64) == 0.0
    cls = _class_index(length, bounds)
    rows: list[dict[str, object]] = []
    for i, label in enumerate(labels):
        sel = cls == i
        n = int(sel.sum())
        if n == 0:
            rows.append({"length_class": label, "n": 0, "k": 0, "r": math.nan,
                         "upper": math.nan, "pi_wrong": math.nan})  # fmt: skip
            continue
        k = n - int(exact_t[sel].sum())
        r = float(int(exact_d[sel].sum()) / n)
        upper = clopper_pearson_upper(k, n)
        pi = 1.0 if r >= 1.0 else min(1.0, upper / (1.0 - r))
        rows.append({"length_class": label, "n": n, "k": k, "r": r,
                     "upper": upper, "pi_wrong": pi})  # fmt: skip
    out = pd.DataFrame(rows, columns=list(STRATUM0_CLASS_COLUMNS))
    out["length_class"] = out["length_class"].astype(str)
    for col in ("n", "k"):
        out[col] = out[col].astype("int64")
    for col in ("r", "upper", "pi_wrong"):
        out[col] = out[col].astype("float64")
    return out


def input_prefiltered(n: int, k: int, mean_pep: float, r: float) -> bool:
    """ADR-0019 decision 4: has the input already been filtered to genome matches?

    An upstream pipeline that keeps only genome-matching peptides leaves no unmapped peptides to
    estimate pi_wrong from, and no genome-based rule can be calibrated on it. The warning fires
    when the number k of stratum-0 database peptides without an exact target hit falls below the
    0.001 quantile of Binomial(n, mean_pep (1 - r)), the count such an input would show if its
    PEP were honest.

    Parameters
    ----------
    n, k:
        Size of the stratum-0 database population and its number of peptides without an exact
        target hit, pooled over the length classes.
    mean_pep:
        Row-level mean PEP of that population; NaN when no peptide carries a PEP.
    r:
        Share of the population with an exact decoy hit.

    Returns
    -------
    bool
        True when the warning is due; False for an empty population or an input without PEP.
    """
    if n <= 0 or not math.isfinite(mean_pep):
        return False
    p = min(1.0, max(0.0, float(mean_pep) * (1.0 - r)))
    return bool(k < float(binom.ppf(PREFILTERED_QUANTILE, n, p)))


def is_low_complexity(peptide: str) -> bool:
    """``low_complexity``: at most 3 distinct residues (I and L counted as one)."""
    return len(set(peptide.upper().replace("I", "L"))) <= LOW_COMPLEXITY_MAX_DISTINCT


def length_class_labels(classes: list[list[int]]) -> list[str]:
    """TSV labels of ``mapping.length_classes``: ``lo`` when lo = hi, else ``lo-hi``."""
    return [str(lo) if lo == hi else f"{lo}-{hi}" for lo, hi in classes]


def best_costs(hits: pd.DataFrame, peptides: pd.DataFrame, freqs: ResidueFreqs) -> pd.DataFrame:
    """One row per mapped (non-decoy) input peptide at its best target and decoy cost.

    ``hits`` has the ``giap.mapping.hits`` schema; ``peptides`` carries ``peptide``, ``q_value``,
    ``pep``, ``is_decoy``, ``source``, ``search_db``, ``pi_denovo`` and optionally
    ``is_contaminant`` (missing = False; M3-25 decision 10: contaminants leave the ADR-0019
    pi-hat population). E at the best cost comes from ``expected_hits`` (M3-09). Input decoy rows
    are dropped. Rows sorted by peptide.
    """
    peps = peptides
    if "is_decoy" in peps.columns:
        peps = peps[~peps["is_decoy"].fillna(False).astype(bool)]
    peps = peps.drop_duplicates("peptide").sort_values("peptide", kind="stable")
    best: dict[str, dict[str, int]] = {"target": {}, "decoy": {}}
    if len(hits):
        g = hits.groupby(["search", "peptide"], sort=True)["n_diff"].min()
        for (search, pep), m in g.items():
            if search in best:
                best[str(search)][str(pep)] = int(m)
    rows: list[dict[str, object]] = []
    for r in peps.to_dict("records"):
        pep = str(r["peptide"])
        m_t = best["target"].get(pep)
        m_d = best["decoy"].get(pep)
        rows.append(
            {
                "peptide": pep,
                "length": len(pep),
                "source": str(r.get("source", "database")),
                "search_db": str(r.get("search_db", "other")),
                "m_target": np.nan if m_t is None else float(m_t),
                "m_decoy": np.nan if m_d is None else float(m_d),
                "e_target": np.nan if m_t is None else expected_hits(pep, m_t, freqs),
                "e_decoy": np.nan if m_d is None else expected_hits(pep, m_d, freqs),
                "pep": _num(r.get("pep")),
                "q_value": _num(r.get("q_value")),
                "pi_denovo": _num(r.get("pi_denovo")),
                "is_contaminant": _flag(r.get("is_contaminant")),
            }
        )
    out = pd.DataFrame(rows, columns=list(_COST_COLUMNS))
    out["length"] = out["length"].astype("int64")
    out["is_contaminant"] = out["is_contaminant"].astype("bool")
    for col in ("m_target", "m_decoy", "e_target", "e_decoy", "pep", "q_value", "pi_denovo"):
        out[col] = out[col].astype("float64")
    return out.reset_index(drop=True)


def _num(x: object) -> float:
    if x is None:
        return math.nan
    try:
        v = float(x)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return math.nan
    return v


def _flag(x: object) -> bool:
    """An optional boolean cell: missing (None, NaN, NA) counts as False."""
    return False if x is None or pd.isna(x) else bool(x)


def _flags(values: pd.Series) -> npt.NDArray[np.bool_]:
    """An optional boolean column as a mask; missing values count as False."""
    return np.array([_flag(v) for v in values], dtype=bool)


def _normalise_costs(costs: pd.DataFrame) -> pd.DataFrame:
    df = costs.copy()
    for col in _COST_COLUMNS:
        if col not in df.columns:
            df[col] = np.nan
    for col in _NUMERIC:
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("float64")
    df["peptide"] = df["peptide"].astype(str)
    df["source"] = df["source"].astype(str)
    df["search_db"] = df["search_db"].astype(str)
    df = df.sort_values("peptide", kind="stable").reset_index(drop=True)
    return df


def _class_index(
    lengths: npt.NDArray[np.float64], classes: list[list[int]]
) -> npt.NDArray[np.int64]:
    idx = np.full(lengths.shape, -1, dtype=np.int64)
    for k, (lo, hi) in enumerate(classes):
        idx[(lengths >= lo) & (lengths <= hi) & (idx < 0)] = k
    return idx


def _row(
    dt: str, m: int, label: str, t: int, d: int, w: float, c: float, n_acc: int
) -> dict[str, object]:
    return {
        "section": "placement",
        "decoy_type": dt,
        "n_diff": m,
        "length_class": label,
        "targets": t,
        "decoys": d,
        "pi_wrong": (w + c) / (d + 1),
        "fdr_stratum": (w + c) / t if t else math.nan,
        "n_accepted": n_acc,
    }


def _stratum0_pooled(
    df: pd.DataFrame, pop: npt.NDArray[np.bool_], table: pd.DataFrame
) -> tuple[int, int, float, float]:
    """(n, k, mean PEP, r) of the stratum-0 database population, pooled over the length classes
    (ADR-0019 decision 4 as read by M3-25 lead decision 4)."""
    n = int(table["n"].sum())
    k = int(table["k"].sum())
    weights = table["n"].to_numpy(dtype=np.float64)
    r = float(np.nansum(table["r"].to_numpy(dtype=np.float64) * weights) / n) if n else 0.0
    pep = df.loc[pop, "pep"].to_numpy(dtype=np.float64)
    pep = pep[~np.isnan(pep)]
    mean_pep = float(math.fsum(pep.tolist()) / pep.size) if pep.size else math.nan
    return n, k, mean_pep, r


def placement(costs: pd.DataFrame, cfg: GiapConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Placement q-values per peptide and the ``mapping_fdr`` rows (``section = placement``).

    Returns ``(per_peptide, mapping_fdr)``, the first two tables of ``placement_tables``.
    """
    out = placement_tables(costs, cfg)
    return out.per_peptide, out.mapping_fdr


def placement_tables(costs: pd.DataFrame, cfg: GiapConfig) -> PlacementTables:
    """Placement q-values, ``mapping_fdr`` rows, the stratum-0 class table and the run flags.

    ``per_peptide`` has one row per peptide with a target hit (columns ``PER_PEPTIDE_COLUMNS``,
    sorted by peptide); ``mapping_fdr`` one row per decoy type x m x length class plus ``all``,
    sorted by decoy type, m and class; ``stratum0_classes`` the ADR-0019 estimate per length
    class over the stratum-0 database population P (every database row except ``same_db`` exact
    peptides and contaminants). Each database peptide of stratum 0 carries the pi_wrong of its
    length class (1 outside every class) and c is the largest class estimate over the classes
    with n > 0; de novo rows keep ``pi_denovo`` and m >= 1 keeps pi_wrong = 1.

    Parameters
    ----------
    costs:
        One row per mapped query peptide (``best_costs``); an absent column is treated as
        missing throughout.
    cfg:
        Configuration; ``mapping.length_classes`` and ``stats.alpha.placement`` are read.

    Returns
    -------
    PlacementTables
        The three tables and the run flags (``input_prefiltered`` when ADR-0019 decision 4's
        binomial test says the input was pre-filtered to genome matches).
    """
    alpha = float(cfg.stats.alpha.placement)
    classes = [list(c) for c in cfg.mapping.length_classes]
    labels = length_class_labels(classes)
    df = _normalise_costs(costs)
    n = len(df)
    m_t = df["m_target"].to_numpy()
    m_d = df["m_decoy"].to_numpy()
    e_t = df["e_target"].to_numpy()
    e_d = df["e_decoy"].to_numpy()
    src = df["source"].to_numpy()
    cls = _class_index(df["length"].to_numpy(), classes)
    is_db = src == "database"
    # Same-strain shortcut: database rows only; de novo rows never (M3-10 lead decision 1).
    same_exact = is_db & (df["search_db"].to_numpy() == "same") & (m_t == 0)
    # ADR-0019 decision 1: π̂ per length class over the stratum-0 database population P, which
    # drops the same_db exact peptides (M3-10 lead decision 2) and the contaminants (M3-25
    # decision 7). PEP and q-values are not read for it (decisions 2 and 3).
    pop = is_db & ~same_exact & ~_flags(df["is_contaminant"])
    table = stratum0_pi_classes(df[pop], classes)
    pi_cls = np.nan_to_num(table["pi_wrong"].to_numpy(dtype=np.float64), nan=1.0)
    pi_db = np.where(cls >= 0, pi_cls[np.clip(cls, 0, None)], 1.0)
    present = table.loc[table["n"] > 0, "pi_wrong"]
    c_db = float(present.max()) if len(present) else 1.0
    flags = (PREFILTERED_FLAG,) if input_prefiltered(*_stratum0_pooled(df, pop, table)) else ()
    pi_dn = df["pi_denovo"].to_numpy()
    pi_dn = np.where(np.isnan(pi_dn), 1.0, np.clip(pi_dn, 0.0, 1.0))
    pi0 = np.where(is_db, pi_db, pi_dn)

    q = np.full(n, np.nan)
    q[same_exact] = 0.0
    rows: list[dict[str, object]] = []
    n_cls = len(classes)
    for source, dt in sorted(DECOY_TYPE_OF.items(), key=lambda kv: kv[1]):
        of_src = src == source
        for m in STRATA:
            searched = of_src & ~same_exact if m == 0 else of_src
            pi = pi0 if m == 0 else np.ones(n)
            if m != 0:
                c = 1.0
            elif source == "database":
                c = c_db  # the largest class π̂; a peptide outside every class does not enter c
            else:
                c = float(pi[searched].max()) if searched.any() else 1.0
            is_t = searched & (m_t == m)
            is_d = searched & (m_d == m)
            if is_t.any():
                q[is_t] = stratum_qvalues(e_t[is_t], e_d[is_d], pi[is_d], c)
            acc = is_t & (q <= alpha)
            t_k = np.bincount(cls[is_t & (cls >= 0)], minlength=n_cls)
            d_k = np.bincount(cls[is_d & (cls >= 0)], minlength=n_cls)
            a_k = np.bincount(cls[acc & (cls >= 0)], minlength=n_cls)
            for k, label in enumerate(labels):
                sel = is_d & (cls == k)
                w = math.fsum(pi[sel].tolist())
                rows.append(_row(dt, m, label, int(t_k[k]), int(d_k[k]), w, c, int(a_k[k])))
            w_all = math.fsum(pi[is_d].tolist())
            rows.append(
                _row(dt, m, "all", int(is_t.sum()), int(is_d.sum()), w_all, c, int(acc.sum()))
            )

    fdr = pd.DataFrame(rows, columns=list(MAPPING_FDR_COLUMNS))
    for col in ("n_diff", "targets", "decoys", "n_accepted"):
        fdr[col] = fdr[col].astype("int64")
    for col in ("pi_wrong", "fdr_stratum"):
        fdr[col] = fdr[col].astype("float64")

    placed = ~np.isnan(m_t)
    pep_flags = []
    for pep, same in zip(df["peptide"].to_numpy()[placed], same_exact[placed], strict=True):
        f = []
        if same:
            f.append("same_db")
        if is_low_complexity(str(pep)):
            f.append("low_complexity")
        pep_flags.append(FLAG_SEP.join(f))
    q_placed = q[placed]
    per = pd.DataFrame(
        {
            "peptide": df["peptide"].to_numpy()[placed],
            "source": src[placed],
            "decoy_type": [DECOY_TYPE_OF.get(s, "") for s in src[placed]],
            "length": df["length"].to_numpy()[placed].astype("int64"),
            "m_target": m_t[placed].astype("int64"),
            "e_target": e_t[placed],
            "placement_q": q_placed,
            "accepted": ~np.isnan(q_placed) & (q_placed <= alpha),
            "flags": pep_flags,
        },
        columns=list(PER_PEPTIDE_COLUMNS),
    )
    return PlacementTables(per, fdr, table, flags)
