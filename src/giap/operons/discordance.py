"""Discordant OP pairs: the within-operon t test and its BH family (GIAP_SPEC.md §7.4.2, §7.4.3,
§8.1; docs/specs/M6-05-discordance.md, pass 1).

The abundance difference Delta of the usable OP pairs (both genes detected with a finite,
non-imputed ``abundance_log2``) is fitted by a Student-t (M6-01); every eligible both-detected OP
pair is then tested two-sided against that fit and the t p-values form one BH family.

Pass 1 (lead decisions 2026-10-21): fewer than ``MIN_OP_PAIRS`` usable pairs or a degenerate fit
leaves every OP pair unestimated with ``too_few_op_pairs``; the one-detected branch is unestimated
with ``no_detection_limit_model`` and ``P(miss)`` never enters the BH family (M6-index decision 1);
``ci_low``/``ci_high`` stay NA with ``unestimated``. Pure and deterministic: no I/O, no RNG.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Protocol

import numpy as np
import numpy.typing as npt
import pandas as pd

from giap.config import GiapConfig
from giap.genetable.abundance import MOLAR_TYPES
from giap.operons.pairs import TEST_COLUMNS
from giap.stats.multitest import bh
from giap.stats.tdist import StudentT, fit_student_t, student_t_two_sided

__all__ = ["MIN_OP_PAIRS", "Discordance", "PMissModel", "discordance", "fit_rows"]

F64 = npt.NDArray[np.float64]
log = logging.getLogger(__name__)

#: fewest usable OP pairs for the t fit (lead decision 1)
MIN_OP_PAIRS: int = 20
ROW_COLUMNS: tuple[str, ...] = ("pair_id", "test", *TEST_COLUMNS, "flags")

TOO_FEW = "too_few_op_pairs"
UNESTIMATED = "unestimated"
#: a tested row has p and q but no CI for a single pair (ADR-0020 decision 4 precedent)
CI_UNESTIMATED = "ci_unestimated"
NO_DL_MODEL = "no_detection_limit_model"
RELATIVE = "relative_scale"


class PMissModel(Protocol):
    """M6-02's detection-limit model: ``P(miss)`` of gene ``u`` given its partner's abundance
    ``a_v`` and the oriented ratio sample (pass 2, M6-05H; unused in pass 1)."""

    def __call__(self, u: str, a_v: float, ratios: F64) -> float: ...


@dataclass(frozen=True, slots=True)
class Discordance:
    """Result of :func:`discordance`.

    Attributes:
        fit: the within-operon t fit (NaN parameters, ``converged = False`` when not made).
        ratios: the Delta values the fit used, sorted ascending.
        rows: one row per OP pair: ``pair_id``, ``test``, the §9.2 test columns and ``flags``.
    """

    fit: StudentT
    ratios: F64
    rows: pd.DataFrame


def _bool(s: pd.Series) -> npt.NDArray[np.bool_]:
    return np.asarray(s.fillna(False).astype(bool).to_numpy(), dtype=np.bool_)


def _usable(pairs: pd.DataFrame) -> npt.NDArray[np.bool_]:
    """Both detected, both abundances finite, neither imputed."""
    a = pairs["abundance_log2_a"].astype(float).to_numpy()
    b = pairs["abundance_log2_b"].astype(float).to_numpy()
    imputed = np.zeros(len(pairs), dtype=bool)
    for col in ("abundance_imputed_a", "abundance_imputed_b"):
        if col in pairs.columns:
            imputed |= _bool(pairs[col])
    ok = _bool(pairs["detected_a"]) & _bool(pairs["detected_b"]) & np.isfinite(a) & np.isfinite(b)
    return np.asarray(ok & ~imputed, dtype=np.bool_)


def fit_rows(pairs: pd.DataFrame) -> F64:
    """Delta of the OP pairs usable for the fit, sorted ascending (so the fit is order-free)."""
    keep = (pairs["class"] == "OP").to_numpy() & _usable(pairs)
    d = pairs["delta_log2"].astype(float).to_numpy()[keep]
    return np.sort(d[np.isfinite(d)])


def _is_relative(genes: pd.DataFrame) -> bool:
    """The run's abundance scale is relative (§7.3.4)."""
    if "abundance_scale" in genes.columns:
        scales = set(genes["abundance_scale"].dropna().astype(str))
        return "relative" in scales
    if "abundance_type" in genes.columns:
        # ``unknown`` counts as relative (§7.3.4)
        types = set(genes["abundance_type"].dropna().astype(str)) - {"none"}
        return bool(types - MOLAR_TYPES)
    return False


def _per_gene(genes: pd.DataFrame, col: str, tags: pd.Series) -> npt.NDArray[np.bool_]:
    if col not in genes.columns:
        return np.zeros(len(tags), dtype=bool)
    m = genes.drop_duplicates("locus_tag").set_index("locus_tag")[col]
    return _bool(tags.map(m))


def _flags(v: object) -> list[str]:
    if isinstance(v, str):
        return [f for f in v.replace(",", ";").split(";") if f]
    if isinstance(v, list | tuple | np.ndarray):
        return [str(f) for f in v]
    return []


def _fit(delta: F64) -> StudentT:
    if delta.size < MIN_OP_PAIRS:
        return StudentT(math.nan, math.nan, math.nan, int(delta.size), False)
    fit = fit_student_t(delta)
    if not fit.converged:
        med = float(np.median(delta))
        mad = 1.4826 * float(np.median(np.abs(delta - med)))
        log.warning(
            "within-operon t fit did not converge on n=%d OP pairs (start loc=%.4g, scale=%.4g, "
            "df=8); every OP pair is unestimated (too_few_op_pairs)",
            delta.size,
            med,
            mad,
        )
    return fit


def discordance(
    pairs: pd.DataFrame, genes: pd.DataFrame, cfg: GiapConfig, p_miss: PMissModel | None
) -> Discordance:
    """Test every OP pair for discordant abundance (§7.4.3, pass 1).

    Args:
        pairs: the M6-04 pair frame (``giap.operons.pairs.pair_frame``).
        genes: the gene table; reads ``locus_tag``, ``partial``, ``abundance_imputed``,
            ``abundance_type`` and ``abundance_scale``.
        cfg: run configuration (``abundance.missing``).
        p_miss: M6-02's detection-limit model; unused in pass 1 (one-detected rows are
            ``unestimated (no_detection_limit_model)``).

    Returns:
        The fit, the Delta sample it used and one row per OP pair; ``q_value`` is BH over the
        rows with a p-value. The input frames are not modified.
    """
    del p_miss, cfg  # pass 1: no detection-limit model; NaN abundances are untested in every mode
    op = pairs.loc[pairs["class"] == "OP"].reset_index(drop=True)
    work = op.copy()
    for side in ("a", "b"):
        tags = work[f"gene_{side}"]
        work[f"abundance_imputed_{side}"] = _per_gene(genes, "abundance_imputed", tags)
        work[f"_partial_{side}"] = _per_gene(genes, "partial", tags)
    eligible = ~(work["_partial_a"].to_numpy() | work["_partial_b"].to_numpy())

    n = len(work)
    rows = pd.DataFrame({"pair_id": work["pair_id"].to_numpy(), "test": ["none"] * n})
    for col in TEST_COLUMNS:
        rows[col] = None if col in ("method", "null_model") else math.nan
    rows["n_perm"] = pd.array([pd.NA] * n, dtype="Int64")
    flags = [_flags(f) for f in (work["flags"] if "flags" in work.columns else [None] * n)]

    delta = fit_rows(work.loc[eligible])
    fit = _fit(delta)

    det_a, det_b = _bool(work["detected_a"]), _bool(work["detected_b"])
    a = work["abundance_log2_a"].astype(float).to_numpy()
    b = work["abundance_log2_b"].astype(float).to_numpy()
    d = work["delta_log2"].astype(float).to_numpy()
    both = eligible & det_a & det_b & np.isfinite(a) & np.isfinite(b) & np.isfinite(d)
    one = eligible & (det_a ^ det_b)

    if _is_relative(genes):
        extra = {i: [RELATIVE, UNESTIMATED] for i in range(n) if eligible[i]}
    elif not fit.converged:
        extra = {i: [TOO_FEW, UNESTIMATED] for i in range(n) if eligible[i]}
    else:
        idx = np.flatnonzero(both)
        p = student_t_two_sided(d[idx], fit)
        rows.loc[idx, "test"] = "both_detected_t"
        rows.loc[idx, "stat"] = d[idx]
        rows.loc[idx, "effect"] = d[idx] - fit.loc
        rows.loc[idx, "p_value"] = p
        rows.loc[idx, "q_value"] = bh(p)
        rows.loc[idx, "method"] = "student_t"
        rows.loc[idx, "null_model"] = "within_operon_t"
        extra = {int(i): [CI_UNESTIMATED] for i in idx}  # no CI for a single pair (decision 2)
        extra.update({int(i): [NO_DL_MODEL, UNESTIMATED] for i in np.flatnonzero(one)})
    for i, add in extra.items():
        flags[i] = flags[i] + [f for f in add if f not in flags[i]]
    rows["flags"] = flags
    return Discordance(fit=fit, ratios=delta, rows=rows[list(ROW_COLUMNS)])
