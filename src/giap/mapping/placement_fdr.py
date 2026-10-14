"""Placement FDR (G3): strata, pi_wrong, q-values and ``mapping_fdr`` rows (issue M3-10).

GIAP_SPEC.md §7.2.5 steps 1-6 as amended by ADR-0010 (docs/methods/placement-fdr.md, "pi_wrong in
the stratum FDR"): within stratum m and decoy type, one list ordered by the ranking score E,
FDR_m(e) = (W_m(e) + c) / T_m(e), W = sum of the decoy hits' own pi_wrong with E <= e, c = max
pi_wrong over the searched peptides of that decoy type (1 for m >= 1). With constant pi this is
pi (D(e) + 1) / T(e). q-values use the M4-01 monotonisation (``giap.stats.tda.monotone_q``).

Pure functions; no file I/O.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
import numpy.typing as npt
import pandas as pd

from giap.config import GiapConfig
from giap.mapping.chance import ResidueFreqs, expected_hits
from giap.stats.tda import monotone_q

__all__ = [
    "DECOY_TYPE_OF",
    "PeptideCosts",
    "best_costs",
    "is_low_complexity",
    "length_class_labels",
    "pi_wrong_database",
    "placement",
    "stratum_qvalues",
]

STRATA: tuple[int, ...] = (0, 1, 2)
DECOY_TYPE_OF: dict[str, str] = {"database": "reverse", "denovo": "shuffle"}
FLAG_SEP = ";"
LOW_COMPLEXITY_MAX_DISTINCT = 3

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


_COST_COLUMNS: tuple[str, ...] = tuple(PeptideCosts.__dataclass_fields__)


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


def pi_wrong_database(peptides: pd.DataFrame, q_value_max: float) -> float:
    """pi_wrong for database input, stratum 0 (§7.2.5 step 4).

    Mean PEP over the peptides with a PEP; otherwise the maximum q-value <= ``q_value_max``;
    otherwise 1.
    """
    if len(peptides) == 0:
        return 1.0
    pep = pd.to_numeric(peptides.get("pep"), errors="coerce").to_numpy(dtype=np.float64)
    pep = pep[~np.isnan(pep)]
    if pep.size:
        return float(min(1.0, max(0.0, math.fsum(pep.tolist()) / pep.size)))
    qv = pd.to_numeric(peptides.get("q_value"), errors="coerce").to_numpy(dtype=np.float64)
    qv = qv[~np.isnan(qv) & (qv <= q_value_max)]
    if qv.size:
        return float(min(1.0, max(0.0, qv.max())))
    return 1.0


def is_low_complexity(peptide: str) -> bool:
    """``low_complexity``: at most 3 distinct residues (I and L counted as one)."""
    return len(set(peptide.upper().replace("I", "L"))) <= LOW_COMPLEXITY_MAX_DISTINCT


def length_class_labels(classes: list[list[int]]) -> list[str]:
    """TSV labels of ``mapping.length_classes``: ``lo`` when lo = hi, else ``lo-hi``."""
    return [str(lo) if lo == hi else f"{lo}-{hi}" for lo, hi in classes]


def best_costs(hits: pd.DataFrame, peptides: pd.DataFrame, freqs: ResidueFreqs) -> pd.DataFrame:
    """One row per mapped (non-decoy) input peptide at its best target and decoy cost.

    ``hits`` has the ``giap.mapping.hits`` schema; ``peptides`` carries ``peptide``, ``q_value``,
    ``pep``, ``is_decoy``, ``source``, ``search_db`` and ``pi_denovo``. E at the best cost comes
    from ``expected_hits`` (M3-09). Input decoy rows are dropped. Rows sorted by peptide.
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
            }
        )
    out = pd.DataFrame(rows, columns=list(_COST_COLUMNS))
    out["length"] = out["length"].astype("int64")
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


def placement(costs: pd.DataFrame, cfg: GiapConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Placement q-values per peptide and the ``mapping_fdr`` rows (``section = placement``).

    Returns ``(per_peptide, mapping_fdr)``. ``per_peptide`` has one row per peptide with a target
    hit (columns ``PER_PEPTIDE_COLUMNS``, sorted by peptide); ``mapping_fdr`` one row per decoy
    type x m x length class plus ``all``, sorted by decoy type, m and class.
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
    same_exact = (df["search_db"].to_numpy() == "same") & (m_t == 0)
    cls = _class_index(df["length"].to_numpy(), classes)
    is_db = src == "database"
    pi_db = pi_wrong_database(df[is_db], float(cfg.input.q_value_max))
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
            c = float(pi[searched].max()) if m == 0 and searched.any() else 1.0
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
    flags = []
    for pep, same in zip(df["peptide"].to_numpy()[placed], same_exact[placed], strict=True):
        f = []
        if same:
            f.append("same_db")
        if is_low_complexity(str(pep)):
            f.append("low_complexity")
        flags.append(FLAG_SEP.join(f))
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
            "flags": flags,
        },
        columns=list(PER_PEPTIDE_COLUMNS),
    )
    return per, fdr
