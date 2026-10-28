"""Discordant OP pairs: the within-operon t test and its BH family (GIPPA_SPEC.md §7.4.2, §7.4.3,
§8.1; docs/specs/M6-05-discordance.md, passes 1 and 2 / M6-05H).

The abundance difference Delta of the usable OP pairs (both genes detected with a finite,
non-imputed ``abundance_log2``) is fitted by a Student-t (M6-01); every eligible both-detected OP
pair is then tested two-sided against that fit and the t p-values form one BH family.

Lead decisions (2026-10-21): fewer than ``MIN_OP_PAIRS`` usable pairs or a degenerate fit leaves
every OP pair unestimated with ``too_few_op_pairs``; ``ci_low``/``ci_high`` stay NA. A pair with
exactly one gene seen (M6-02 lead decision 9: ``n_seen >= k_min`` observable peptides; without
``n_seen``, ``detected``) and no t test is scored by M6-02's ``P(miss)`` of the unseen gene u
(``one_detected_pmiss``, M6-05H) when a ``p_miss`` callable is given, else (or when it returns NaN)
it is unestimated with ``no_detection_limit_model`` (or the caller's ``no_model_flag``). A row
whose seen partner v has no abundance is untested with no flag. ``P(miss)`` is a probability, not
a p-value: it gets no q and never enters the BH family (M6-index decision 1, M6-02 decision 4).
M6-13 (ADR-0055 option (a)): by default the t scale is log-linear in the pair-mean log2 abundance
m (``sigma_model = "loglinear"``, :func:`gippa.stats.tdist.fit_student_t_scaled`); ``"constant"``
is the M6-05 null. P(miss) still receives the raw Delta of the fit pairs (M6-13 Q1).
Pure and deterministic: no I/O, no RNG.
"""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping
from dataclasses import dataclass, replace
from typing import Literal, Protocol

import numpy as np
import numpy.typing as npt
import pandas as pd

from gippa.config import GippaConfig
from gippa.genetable.abundance import MOLAR_TYPES
from gippa.operons.pairs import TEST_COLUMNS
from gippa.stats.multitest import bh
from gippa.stats.tdist import (
    ScaledStudentT,
    StudentT,
    fit_student_t,
    fit_student_t_scaled,
    scaled_t_two_sided,
)

__all__ = ["MIN_OP_PAIRS", "SIGMA_MODEL", "Discordance", "PMissFn", "discordance", "fit_pairs",
           "fit_rows"]  # fmt: skip

F64 = npt.NDArray[np.float64]
log = logging.getLogger(__name__)

#: fewest usable OP pairs for the t fit (lead decision 1)
MIN_OP_PAIRS: int = 20
#: the default scale model (M6-13 Q2: no config key yet)
SIGMA_MODEL: Literal["loglinear", "constant"] = "loglinear"
#: scaled-fit fallbacks whose step event is a WARNING (M6-13 step 5)
FALLBACK_WARN = frozenset({"not_converged", "collapsed", "worse_than_constant"})
ROW_COLUMNS: tuple[str, ...] = ("pair_id", "test", *TEST_COLUMNS, "flags")

TOO_FEW = "too_few_op_pairs"
UNESTIMATED = "unestimated"
#: a tested row has p and q but no CI for a single pair (ADR-0020 decision 4 precedent)
CI_UNESTIMATED = "ci_unestimated"
NO_DL_MODEL = "no_detection_limit_model"
RELATIVE = "relative_scale"
PMISS_TEST = "one_detected_pmiss"


class PMissFn(Protocol):
    """M6-02's detection-limit model as a callable: ``P(miss)`` of the unseen gene ``u`` given
    its seen partner ``v``'s abundance ``a_v`` and the oriented ratio sample r (§7.4.3); ``v``
    selects the out-of-fold model (M6-02 lead decision 9). NaN when it cannot be scored."""

    def __call__(self, u: str, a_v: float, ratios: F64, v: str, /) -> float: ...


@dataclass(frozen=True, slots=True)
class Discordance:
    """Result of :func:`discordance`.

    Attributes:
        fit: the constant within-operon t fit (NaN parameters, ``converged = False`` when not made).
        ratios: the Delta values the fit used, sorted ascending.
        rows: one row per OP pair: ``pair_id``, ``test``, the §9.2 test columns and ``flags``.
        scaled: the null of every both-detected p (M6-13); ``model = "constant"`` is ``fit``.
    """

    fit: StudentT
    ratios: F64
    rows: pd.DataFrame
    scaled: ScaledStudentT


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


def fit_pairs(pairs: pd.DataFrame) -> tuple[F64, F64]:
    """(Delta, m) of the OP pairs usable for the fit, m the pair-mean log2 abundance, sorted
    lexicographically on (Delta, m)."""
    keep = (pairs["class"] == "OP").to_numpy() & _usable(pairs)
    d = pairs["delta_log2"].astype(float).to_numpy()[keep]
    a = pairs["abundance_log2_a"].astype(float).to_numpy()[keep]
    m = (a + pairs["abundance_log2_b"].astype(float).to_numpy()[keep]) / 2
    ok = np.isfinite(d)
    order = np.lexsort((m[ok], d[ok]))
    return d[ok][order], m[ok][order]


def _scaled(
    x: F64, m: F64, sigma_model: Literal["loglinear", "constant"], relative: bool
) -> ScaledStudentT:
    """The M6-13 null; ``constant`` never runs the four-parameter fit (``reason = disabled``), a
    relative-scale run computes no p and has no null (``none``, ``relative_scale``)."""
    if sigma_model not in ("loglinear", "constant"):
        raise ValueError(f"sigma_model must be 'loglinear' or 'constant', got {sigma_model!r}")
    if relative:  # min_n above the pair count: no fit is made
        return replace(fit_student_t_scaled(x, m, min_n=x.size + 1), reason=RELATIVE)
    if sigma_model == "constant":
        s = fit_student_t_scaled(x, m, min_sd=math.inf, min_n=MIN_OP_PAIRS)
        return s if s.model == "none" else replace(s, reason="disabled")
    return fit_student_t_scaled(x, m, min_n=MIN_OP_PAIRS)


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


def _one_detected(
    rows: pd.DataFrame,
    work: pd.DataFrame,
    idx: npt.NDArray[np.intp],
    v_is_a: npt.NDArray[np.bool_],
    delta: F64,
    p_miss: PMissFn | None,
    no_model_flag: str,
) -> dict[int, list[str]]:
    """Score the one-detected rows ``idx`` in place with ``P(miss)`` (spec item 4) and return the
    flags of the rows left unestimated. The unseen gene u is ``gene_b`` (downstream, r = Δ, when
    ``v_is_a``) or ``gene_a`` (upstream, r = -Δ); ``a_v`` is the seen partner v's
    ``abundance_log2``. A row whose ``a_v`` is missing is untested with no flag, with or without a
    model (item 5); no q is written (decision 1)."""
    a = work["abundance_log2_a"].astype(float).to_numpy()
    b = work["abundance_log2_b"].astype(float).to_numpy()
    upstream = -delta[::-1]  # still ascending
    out: dict[int, list[str]] = {}
    for i in (int(j) for j in idx):
        down = bool(v_is_a[i])  # v = gene_a: u = gene_b is downstream
        a_v = float(a[i] if down else b[i])
        if not math.isfinite(a_v):
            continue
        if p_miss is None:
            out[i] = [no_model_flag, UNESTIMATED]
            continue
        u = str(work.at[i, "gene_b" if down else "gene_a"])
        v = str(work.at[i, "gene_a" if down else "gene_b"])
        p = float(p_miss(u, a_v, delta.copy() if down else upstream.copy(), v))
        if not math.isfinite(p):
            out[i] = [NO_DL_MODEL, UNESTIMATED]
            continue
        rows.loc[i, ["test", "stat", "p_value", "method", "null_model"]] = [
            PMISS_TEST, a_v, p, "detection_limit", "p_detect_model"]  # fmt: skip
        out[i] = [CI_UNESTIMATED]  # no effect or CI for P(miss) (§8; ADR-0056 decision 4)
    return out


def _seen(tags: pd.Series, n_seen: Mapping[str, int], k_min: int) -> npt.NDArray[np.bool_]:
    return np.array([int(n_seen.get(str(t), 0)) >= k_min for t in tags], dtype=bool)


def discordance(
    pairs: pd.DataFrame,
    genes: pd.DataFrame,
    cfg: GippaConfig,
    p_miss: PMissFn | None,
    *,
    n_seen: Mapping[str, int] | None = None,
    k_min: int = 1,
    no_model_flag: str = NO_DL_MODEL,
    sigma_model: Literal["loglinear", "constant"] = SIGMA_MODEL,
) -> Discordance:
    """Test every OP pair for discordant abundance (§7.4.3).

    Args:
        pairs: the M6-04 pair frame (``gippa.operons.pairs.pair_frame``).
        genes: the gene table; reads ``locus_tag``, ``partial``, ``abundance_imputed``,
            ``abundance_type`` and ``abundance_scale``.
        cfg: run configuration (``abundance.missing``).
        p_miss: M6-02's detection-limit model as ``p_miss(u, a_v, ratios, v)``, or ``None`` when
            there is none (one-detected rows are then ``unestimated`` with ``no_model_flag``).
        n_seen: locus tag -> distinct observable peptides seen in C11 (absent = 0); a gene is
            seen when ``n_seen >= k_min`` and a one-detected row has exactly one gene seen
            (M6-02 lead decision 9). ``None``: the event is ``detected_a xor detected_b``.
        k_min: the ``k_min`` of the model in force.
        no_model_flag: the reason flag of one-detected rows when ``p_miss`` is ``None``
            (``no_detection_limit_model``; the step passes ``pmiss_covariates_missing`` when the
            model in force needs covariates it cannot join).
        sigma_model: ``loglinear`` (M6-13, the default) or ``constant`` (the M6-05 null, bit for
            bit).

    Returns:
        The fit, the Delta sample it used and one row per OP pair; ``q_value`` is BH over the
        ``both_detected_t`` rows only (``one_detected_pmiss`` rows have no q). The input frames
        are not modified.
    """
    del cfg  # NaN abundances are untested in every ``abundance.missing`` mode
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
    relative = _is_relative(genes)
    scaled = _scaled(*fit_pairs(work.loc[eligible]), sigma_model, relative)

    det_a, det_b = _bool(work["detected_a"]), _bool(work["detected_b"])
    a = work["abundance_log2_a"].astype(float).to_numpy()
    b = work["abundance_log2_b"].astype(float).to_numpy()
    d = work["delta_log2"].astype(float).to_numpy()
    both = eligible & det_a & det_b & np.isfinite(a) & np.isfinite(b) & np.isfinite(d)
    if n_seen is None:
        seen_a, seen_b = det_a, det_b
    else:
        seen_a = _seen(work["gene_a"], n_seen, k_min)
        seen_b = _seen(work["gene_b"], n_seen, k_min)
    one = eligible & (seen_a ^ seen_b) & ~both  # the t test takes precedence

    if relative:
        extra = {i: [RELATIVE, UNESTIMATED] for i in range(n) if eligible[i]}
    elif not fit.converged:
        extra = {i: [TOO_FEW, UNESTIMATED] for i in range(n) if eligible[i]}
    else:
        idx = np.flatnonzero(both)
        p = scaled_t_two_sided(d[idx], (a[idx] + b[idx]) / 2, scaled)  # constant: = M6-05 p
        rows.loc[idx, "test"] = "both_detected_t"
        rows.loc[idx, "stat"] = d[idx]
        rows.loc[idx, "effect"] = d[idx] - scaled.loc
        rows.loc[idx, "p_value"] = p
        rows.loc[idx, "q_value"] = bh(p)  # the family: t p-values only (M6-index decision 1)
        rows.loc[idx, "method"] = "student_t"
        rows.loc[idx, "null_model"] = "within_operon_t"
        extra = {int(i): [CI_UNESTIMATED] for i in idx}  # no CI for a single pair (decision 2)
        extra.update(
            _one_detected(rows, work, np.flatnonzero(one), seen_a, delta, p_miss, no_model_flag)
        )
    for i, add in extra.items():
        flags[i] = flags[i] + [f for f in add if f not in flags[i]]
    rows["flags"] = flags
    return Discordance(fit=fit, ratios=delta, rows=rows[list(ROW_COLUMNS)], scaled=scaled)
