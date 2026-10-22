"""Detection enrichment: conditional Poisson-binomial ORA counted in TUs (spec §7.7.3; issue
M8-04; ADR-0029, ADR-0035).

Pure. The tested units are TUs. A TU's raw detection probability is
q_raw = 1 - prod(1 - p_detect) over its eligible (non-partial) genes, with p_detect clipped to
[1e-12, 1 - 1e-12]. Given M4-13's TU folds, q is a cross-fitted logistic recalibration of TU
detection on logit q_raw and log(eligible genes), so a TU's own detection never enters its fitted
value (decision 12). A TU is in a set when any of its eligible genes is. The p-value of a set C is
the law of its detected-TU count given the total N (``conditional_poisson_binomial`` on
``tilt_pair(q, N, comp=)``, tilted once per run). GO namespaces are tested ``over`` with elim on
TU sets; BH runs per (collection, direction). ``p_classic`` is the gene-level Wallenius tail
(M8-01, for goseq parity); it is never elim'd and never enters BH.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.special import expit

from giap.function.collections import COLLECTIONS, GeneSets
from giap.function.elim import elim
from giap.function.gsea import _closure, _depth, _natural_key
from giap.stats.multitest import bh
from giap.stats.poisson_binomial import CpbTails, conditional_poisson_binomial, tilt_pair
from giap.stats.wallenius import wallenius_ora

__all__ = ["COLUMNS", "TuUniverse", "detection_enrichment", "eligible_genes", "tu_members",
           "tu_universe"]  # fmt: skip

F64 = npt.NDArray[np.float64]
I64 = npt.NDArray[np.int64]
B8 = npt.NDArray[np.bool_]

CLIP = 1e-12  # p_detect is clipped to [CLIP, 1 - CLIP] (decision 12)
MIN_CLASS = 10  # detected and undetected training TUs each fold fit needs, else no recalibration
ELIM_THRESHOLD = 0.01  # strict: an elim p below it removes the term's TUs from its ancestors
GO_NAMESPACES = frozenset({"go_bp", "go_mf", "go_cc"})
DIRECTIONS = ("over", "under")
METHOD = "cpb_tu_ora"
NULL_MODEL = "conditional_poisson_binomial_tu"
#: Flags every row carries while a decision-12 attack arm fails after recalibration (ADR-0034
#: pattern). On since decision 13: arm 4b (c = 2.5, large, under) and 4d-4f fail at R = 1000;
#: lifted only on a pass of 4b-4f at R = 1000 in both blocks.
CALIBRATION_FLAGS: tuple[str, ...] = ("ora_uncalibrated",)
MAX_NEWTON = 50  # Newton iterations of the recalibration fit
NEWTON_TOL = 1e-10  # converged once max |step| is below it
MAX_HALVINGS = 60  # step halvings per iteration before the fit is declared stalled
COLUMNS = (
    "collection",
    "term",
    "name",
    "direction",
    "n_universe",
    "n_in_set",
    "n_detected",
    "k_detected_in_set",
    "n_in_set_genes",
    "k_detected_in_set_genes",
    "odds_weight",
    "p_classic",
    "p_noelim",
    "parent_group",
    "coverage",
    "stat",
    "effect",
    "ci_low",
    "ci_high",
    "p_value",
    "q_value",
    "method",
    "null_model",
    "n_perm",
    "flags",
)
_FLOAT = ("odds_weight", "p_classic", "p_noelim", "coverage", "stat", "effect", "ci_low",
          "ci_high", "p_value", "q_value")  # fmt: skip
_INT = ("n_universe", "n_in_set", "n_detected", "k_detected_in_set", "n_in_set_genes",
        "k_detected_in_set_genes", "n_perm")  # fmt: skip
_ORDER = {c: i for i, c in enumerate(COLLECTIONS)}


@dataclass(frozen=True, slots=True)
class TuUniverse:
    """TUs with >= 1 eligible gene, in natural TU-ID order."""

    tu_ids: tuple[str, ...]
    q: F64  # detection probability the null uses (recalibrated, or q_raw)
    comp: F64  # 1 - q, computed directly
    q_raw: F64  # 1 - prod(1 - p_detect) over the TU's eligible genes, p_detect clipped
    log_comp: F64  # sum of log1p(-p_detect), clipped
    n_genes: I64  # eligible genes per TU
    detected: B8  # any eligible gene detected
    clipped: B8  # a gene's p_detect was moved by the clip
    gene_tu: Mapping[str, int]  # eligible locus_tag -> TU index
    flags: tuple[str, ...]  # tu_recal_skipped, p_detect_clipped


def _missing_tu(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v)) or not str(v)


def eligible_genes(genes: pd.DataFrame) -> pd.DataFrame:
    """Non-partial genes sorted by (TU, locus tag) in natural order; fail closed (spec §2).

    Raises:
        ValueError: naming the eligible genes with a NaN ``p_detect``, a ``p_detect`` outside
            [0, 1] (decision 13; exact 0 and 1 are clipped later) or no ``tu_id``.
    """
    keep = np.ones(len(genes), dtype=bool)
    if "partial" in genes:
        keep = ~genes["partial"].astype("boolean").fillna(False).to_numpy(bool)
    el = genes[keep]
    tags = el["locus_tag"].astype(str)
    p = pd.to_numeric(el["p_detect"], errors="coerce").to_numpy(float)
    finite = np.isfinite(p)
    no_p = sorted(tags[~finite])
    out_p = sorted(tags[finite & ((p < 0.0) | (p > 1.0))])
    no_tu = sorted(t for t, v in zip(tags, el["tu_id"], strict=True) if _missing_tu(v))
    checks = (("NaN p_detect", no_p), ("p_detect outside [0, 1]", out_p), ("no tu_id", no_tu))
    problems = [f"{what}: {', '.join(bad)}" for what, bad in checks if bad]
    if problems:
        raise ValueError(f"eligible genes with {'; '.join(problems)}")
    key = [(_natural_key(str(t)), _natural_key(g)) for t, g in zip(el["tu_id"], tags, strict=True)]
    order = sorted(range(len(el)), key=key.__getitem__)
    return el.iloc[order].reset_index(drop=True)


def _loglik(x: F64, y: F64, beta: F64) -> float:
    eta = x @ beta
    return float(np.sum(y * eta - np.logaddexp(0.0, eta)))


def _logistic_mle(x: F64, y: F64) -> F64 | None:
    """Unpenalised logistic MLE by Newton from the identity on logit q, each step halved until
    the log-likelihood does not decrease (decision 13); converged once the full Newton step is
    below ``NEWTON_TOL``. None on a singular Hessian, a non-finite step, a line search that
    finds no increase in ``MAX_HALVINGS`` halvings or no convergence in ``MAX_NEWTON`` steps."""
    beta = np.zeros(x.shape[1])
    beta[1] = 1.0
    ll = _loglik(x, y, beta)
    for _ in range(MAX_NEWTON):
        mu = expit(x @ beta)
        try:
            step = np.linalg.solve(x.T @ (x * (mu * (1.0 - mu))[:, None]), x.T @ (y - mu))
        except np.linalg.LinAlgError:
            return None
        if not np.all(np.isfinite(step)):
            return None
        if np.max(np.abs(step)) < NEWTON_TOL:  # the full step, so halving cannot fake it
            return beta + step
        for _ in range(MAX_HALVINGS):
            new = beta + step
            ll_new = _loglik(x, y, new)
            if ll_new >= ll:  # near the optimum a step can tie on float noise
                break
            step = 0.5 * step
        else:
            return None
        beta, ll = new, ll_new
    return None


def _recalibrate(logit_raw: F64, n_genes: I64, detected: B8, fold: I64) -> tuple[F64, F64] | None:
    """Cross-fitted (q, comp) of d ~ logistic(b0 + b1 logit q_raw + b2 log n) per fold; None when
    a fold's training TUs have fewer than ``MIN_CLASS`` of either class or the fit fails. The size
    term is dropped when the training TUs all have the same size."""
    y = detected.astype(np.float64)
    log_n = np.log(n_genes.astype(np.float64))
    eta = np.empty_like(logit_raw)
    for j in np.unique(fold):
        train, test = fold != j, fold == j
        n_pos = int(y[train].sum())
        if min(n_pos, int(train.sum()) - n_pos) < MIN_CLASS:
            return None
        cols = [np.ones_like(logit_raw), logit_raw]
        if np.ptp(log_n[train]) > 0:
            cols.append(log_n)
        x = np.column_stack(cols)
        beta = _logistic_mle(x[train], y[train])
        if beta is None:
            return None
        eta[test] = x[test] @ beta
    return expit(eta), expit(-eta)


def tu_universe(genes: pd.DataFrame, folds: Mapping[str, int] | None) -> TuUniverse:
    """The TU units, raw and recalibrated q (decision 12); ``folds`` is M4-13's ``fold_of_tu``.

    ``folds`` has no default (decision 13): None (or empty) is the caller's explicit choice to
    skip the fit. Skipped, q is ``q_raw`` and ``flags`` holds ``tu_recal_skipped``, as when a
    fold is short of a class or the fit fails.

    Raises:
        ValueError: as ``eligible_genes``; naming the TUs missing from a non-empty ``folds``.
    """
    el = eligible_genes(genes)
    tags = el["locus_tag"].astype(str).tolist()
    tus = el["tu_id"].astype(str).tolist()
    ids = tuple(dict.fromkeys(tus))  # already in natural order
    index = {t: i for i, t in enumerate(ids)}
    tu_of = np.array([index[t] for t in tus], dtype=np.int64)
    n_tu = len(ids)
    p = pd.to_numeric(el["p_detect"], errors="coerce").to_numpy(float)
    p_clip = np.clip(p, CLIP, 1.0 - CLIP)
    det = el["detected"].astype("boolean").fillna(False).to_numpy(bool)
    log_comp = np.bincount(tu_of, weights=np.log1p(-p_clip), minlength=n_tu)
    n_genes = np.bincount(tu_of, minlength=n_tu).astype(np.int64)
    detected = np.bincount(tu_of, weights=det.astype(np.float64), minlength=n_tu) > 0
    clipped = np.bincount(tu_of, weights=(p_clip != p).astype(np.float64), minlength=n_tu) > 0
    q_raw = -np.expm1(log_comp)

    fitted = None
    if folds:
        missing = [t for t in ids if t not in folds]
        if missing:
            raise ValueError(f"TUs without a fold: {', '.join(missing)}")
        fold = np.array([int(folds[t]) for t in ids], dtype=np.int64)
        fitted = _recalibrate(np.log(q_raw) - log_comp, n_genes, detected, fold)
    flags = []
    if fitted is None:
        q, comp = q_raw, np.exp(log_comp)
        flags.append("tu_recal_skipped")
    else:
        q, comp = fitted
    if clipped.any():
        flags.append("p_detect_clipped")
    return TuUniverse(
        tu_ids=ids,
        q=q,
        comp=comp,
        q_raw=q_raw,
        log_comp=log_comp,
        n_genes=n_genes,
        detected=detected,
        clipped=clipped,
        gene_tu=dict(zip(tags, tu_of.tolist(), strict=True)),
        flags=tuple(flags),
    )


def tu_members(u: TuUniverse, genes: frozenset[str]) -> frozenset[int]:
    """Indices of the TUs with any of ``genes`` (eligible genes outside ``u`` are ignored)."""
    return frozenset(u.gene_tu[g] for g in genes if g in u.gene_tu)


def _unit(p: float) -> float:
    return min(max(p, 0.0), 1.0)


class _Null:
    """The conditional law of a TU set's detected count given N, on the run's tilted (q, comp)."""

    def __init__(self, u: TuUniverse) -> None:
        self.detected = u.detected
        self.n = int(u.detected.sum())
        # an undetected certain TU or a detected impossible one: the observation has probability 0
        certain_miss = np.any((u.comp == 0.0) & ~u.detected)
        self.degenerate = bool(certain_miss or np.any((u.q == 0.0) & u.detected))
        self.q: F64 = u.q
        self.comp: F64 = u.comp
        if not self.degenerate:
            self.q, self.comp = tilt_pair(u.q, self.n, comp=u.comp)
        self._memo: dict[frozenset[int], CpbTails | None] = {}

    def k(self, tus: frozenset[int]) -> int:
        return int(self.detected[sorted(tus)].sum())

    def tails(self, tus: frozenset[int]) -> CpbTails | None:
        """None when not a test: an empty set, the whole universe or a degenerate run. Both tails
        are clamped to [0, 1] here too, so BH never sees 1 + ulp (decision 13)."""
        if tus not in self._memo:
            law = None
            if tus and len(tus) < self.q.size and not self.degenerate:
                m = np.zeros(self.q.size, dtype=bool)
                m[sorted(tus)] = True
                law = conditional_poisson_binomial(
                    self.q[m], self.q[~m], self.k(tus), self.n,
                    comp_in=self.comp[m], comp_out=self.comp[~m],
                )  # fmt: skip
                law = replace(law, p_over=_unit(law.p_over), p_under=_unit(law.p_under))
            self._memo[tus] = law
        return self._memo[tus]

    def p_over(self, tus: frozenset[int]) -> float:
        law = self.tails(tus)
        return math.nan if law is None else law.p_over


def _elim_sets(
    tus: Mapping[str, frozenset[int]],
    parents: Mapping[str, frozenset[str]],
    p_elim: Mapping[str, float],
) -> dict[str, frozenset[int]]:
    """The TU set each term was elim-tested on: its TUs minus those of its significant
    descendants' own tested sets (M8-01's walker, replayed from its p-values)."""
    children: dict[str, set[str]] = {}
    for t, ps in parents.items():
        for a in ps:
            children.setdefault(a, set()).add(t)
    below = _closure(children)
    memo: dict[str, frozenset[int]] = {}

    def left(t: str) -> frozenset[int]:
        if t not in memo:
            drop = [left(d) for d in sorted(below(t)) if p_elim[d] < ELIM_THRESHOLD]
            memo[t] = tus[t].difference(*drop)
        return memo[t]

    return {t: left(t) for t in tus}


def _classic(det: B8, p: F64, in_cat: B8) -> tuple[float, dict[str, float], list[str]]:
    """Gene-level Wallenius (odds, {direction: p}, flags); NA with ``odds_out_of_range`` when the
    M8-01 primitive fails closed."""
    try:
        w = wallenius_ora(det, p, in_cat)
    except ValueError:
        return math.nan, dict.fromkeys(DIRECTIONS, math.nan), ["odds_out_of_range"]
    return w.odds, {"over": w.p_over, "under": w.p_under}, []


def _collection_rows(
    coll: GeneSets,
    u: TuUniverse,
    null: _Null,
    size: tuple[int, int],
    el: pd.DataFrame,
) -> list[dict[str, Any]]:
    tus = {t: tu_members(u, s) for t, s in coll.sets.items()}
    tested = sorted(t for t, s in tus.items() if size[0] <= len(s) <= size[1])
    go = coll.collection in GO_NAMESPACES
    tested_sets = {t: tus[t] for t in tested}
    elim_set: dict[str, frozenset[int]] = {}
    if go and tested:
        above = _closure(coll.parents)
        keep = set(tested)
        par = {t: frozenset(above(t) & keep) for t in tested}
        walked = elim(tested_sets, par, null.p_over, ELIM_THRESHOLD)
        p_elim = {t: pe for t, (_, pe) in walked.items()}
        elim_set = _elim_sets(tested_sets, par, p_elim)

    gene_pos = {g: i for i, g in enumerate(el["locus_tag"].astype(str))}
    gene_det = el["detected"].astype("boolean").fillna(False).to_numpy(bool)
    gene_p = pd.to_numeric(el["p_detect"], errors="coerce").to_numpy(float)
    rows: list[dict[str, Any]] = []
    for t in tested:
        full = tus[t]
        k = null.k(full)
        in_cat = np.zeros(gene_p.size, dtype=bool)
        in_cat[[gene_pos[g] for g in coll.sets[t] if g in gene_pos]] = True
        odds, p_classic, base_flags = _classic(gene_det, gene_p, in_cat)
        full_law = null.tails(full)
        for direction in DIRECTIONS:
            flags = ["ci_unestimated", *base_flags]
            p_noelim = math.nan
            if full_law is not None:
                p_noelim = full_law.p_over if direction == "over" else full_law.p_under
            elim_row = go and direction == "over"
            tested_tus = elim_set[t] if elim_row else full
            law = null.tails(tested_tus)
            if elim_row and not tested_tus:
                flags.append("elim_empty")
            p_value = math.nan
            if law is not None:
                p_value = law.p_over if direction == "over" else law.p_under
            effect = math.nan if law is None else null.k(tested_tus) - law.mean
            if math.isnan(p_value):
                flags.append("unestimated")
            if null.degenerate:
                flags.append("null_degenerate")
            if u.clipped[sorted(full)].any():
                flags.append("p_detect_clipped")
            if "tu_recal_skipped" in u.flags:
                flags.append("tu_recal_skipped")
            rows.append(
                {
                    "collection": coll.collection,
                    "term": t,
                    "name": coll.names.get(t, "NA"),
                    "direction": direction,
                    "n_universe": len(u.tu_ids),
                    "n_in_set": len(full),
                    "n_detected": null.n,
                    "k_detected_in_set": k,
                    "n_in_set_genes": int(in_cat.sum()),
                    "k_detected_in_set_genes": int((in_cat & gene_det).sum()),
                    "odds_weight": odds,
                    "p_classic": p_classic[direction],
                    "p_noelim": p_noelim,
                    "parent_group": None,
                    "coverage": float(coll.coverage),
                    "stat": float(k),
                    "effect": effect,
                    "ci_low": math.nan,
                    "ci_high": math.nan,
                    "p_value": p_value,
                    "q_value": math.nan,
                    "method": METHOD,
                    "null_model": NULL_MODEL,
                    "n_perm": 0,
                    "flags": [*flags, *CALIBRATION_FLAGS],
                }
            )
    return rows


def _adjust(rows: list[dict[str, Any]], coll: GeneSets, alpha: float) -> None:
    """BH per direction on ``p_value`` (NaN stays out), then ``parent_group`` for GO rows."""
    for direction in DIRECTIONS:
        fam = [r for r in rows if r["direction"] == direction]
        q = bh(np.array([r["p_value"] for r in fam], dtype=np.float64))
        for r, qv in zip(fam, q, strict=True):
            r["q_value"] = float(qv)
        if coll.collection not in GO_NAMESPACES:
            continue
        sig = {r["term"] for r in fam if r["q_value"] <= alpha}
        above, depth = _closure(coll.parents), _depth(coll.parents)
        for r in fam:
            if r["term"] in sig:
                anc = above(r["term"]) & sig
                r["parent_group"] = min(anc, key=lambda a: (depth(a), a)) if anc else r["term"]


def detection_enrichment(
    genes: pd.DataFrame,
    colls: Sequence[GeneSets],
    size: tuple[int, int],
    alpha: float,
    *,
    folds: Mapping[str, int] | None,
) -> pd.DataFrame:
    """The enrichment_detection rows of every collection (§9.3, spec M8-04 §4), sorted by
    collection, term and direction; run flags in ``attrs["run_flags"]``. Terms whose full TU set
    lies outside ``size`` write no row and do not enter elim. ``folds`` is required (decision
    13); None skips the recalibration and every row gets ``tu_recal_skipped``. A degenerate run
    (the observation has probability 0 under q) emits one ``UserWarning``.

    Raises:
        ValueError: as ``tu_universe``.
    """
    u = tu_universe(genes, folds)
    null = _Null(u)
    if null.degenerate:
        warnings.warn(
            "detection ORA: the observed detection has probability 0 under the TU null "
            "(an undetected TU with q = 1); every row is NA with null_degenerate",
            UserWarning,
            stacklevel=2,
        )
    el = eligible_genes(genes)
    rows: list[dict[str, Any]] = []
    for coll in colls:
        coll_rows = _collection_rows(coll, u, null, size, el)
        _adjust(coll_rows, coll, alpha)
        rows += coll_rows
    rows.sort(key=lambda r: (_ORDER[r["collection"]], r["term"], r["direction"]))
    df = pd.DataFrame(rows, columns=list(COLUMNS))
    for c in _FLOAT:
        df[c] = df[c].astype("float64")
    for c in _INT:
        df[c] = df[c].astype("int64")
    df["parent_group"] = df["parent_group"].astype(object)
    df.attrs["run_flags"] = [*u.flags, *(["null_degenerate"] if null.degenerate else [])]
    return df
