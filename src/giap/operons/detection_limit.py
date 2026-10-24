"""Detection-limit model: P(detect | abundance) and P(miss) for the one-detected pair test
(GIAP_SPEC.md §7.4.3 "One detected, one not", §7.3.3, §8.7; docs/specs/M6-02; method
docs/methods/detection-limit.md).

Both training sets are selected on the response, so both models are fitted by an unpenalised
conditional (truncated) likelihood (lead decision 5; Huggins 1989):

- peptide model: ``logit p_j = b0 + b_a a + b_L L_j + b_G GRAVY_j`` over the observable peptides
  of the training genes (finite ``abundance_log2``, >= ``k_min`` observable peptides seen in
  C11), conditioned on ``X_g >= k_min``; cross-fitted by TU folds;
- replicate model: per-run ``logit q_g = c0 + c_a a_g + covariates`` over genes seen in
  1..R-1 of R runs, a doubly truncated binomial; not identified at R = 2 (``replicate_two_runs``).

Separation or non-convergence gives ``usable = False`` (no divergent iterate is reported).
``p_miss`` is a probability, never a p-value (lead decision 4). Pure: no I/O.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import numpy.typing as npt
import pandas as pd
from numpy.typing import ArrayLike
from scipy.special import expit, log_expit
from scipy.stats import rankdata

from giap.config import GiapConfig
from giap.genetable.covariates import digest
from giap.genetable.detectability import COLUMNS, transform_covariates, tu_folds
from giap.stats.poisson_binomial import poisson_binomial_sf

__all__ = [
    "DetectionLimitFit",
    "PMissModel",
    "fit_detection_limit",
    "fit_peptide_model",
    "fit_replicate_model",
    "model_for_tu",
    "observable_peptides",
    "p_detect",
    "p_miss",
    "p_miss_for_genes",
]

Array = npt.NDArray[np.float64]
Kind = Literal["peptide", "replicate"]
MIN_ROWS = 50
MAX_ITER = 100
COEF_BOUND = 20.0  # |standardised coefficient| beyond this: (quasi-)separated, no MLE
SECTION = "detection_limit"
PEPTIDE_FEATURES = ("abundance_log2", "length", "gravy")
# Kyte & Doolittle (1982) hydropathy, the ProtParam GRAVY scale M4-10 uses.
_KD: Mapping[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}  # fmt: skip


@dataclass(frozen=True, slots=True)
class PMissModel:
    """A fitted detection-limit model. ``fold_coef[j]`` was trained on the TUs whose fold
    (``fold_of_tu``) is not ``j``; ``reason`` says why a model is not usable."""

    kind: Kind
    coef: Mapping[str, float]
    k_min: int
    n_fit: int
    auc: float
    usable: bool
    fold_of_tu: Mapping[str, int] = field(default_factory=dict)
    fold_coef: tuple[Mapping[str, float], ...] = ()
    reason: str | None = None


@dataclass(frozen=True, slots=True)
class DetectionLimitFit:
    model: PMissModel  # the model in force
    peptide: PMissModel | None
    replicate: PMissModel | None
    summary: pd.DataFrame  # detection_summary rows, section "detection_limit"


# ------------------------------------------------------------------------------ scoring


def _gravy(pep: str) -> float:
    vals = [_KD[c] for c in pep if c in _KD]
    return sum(vals) / len(vals) if vals else 0.0


def observable_peptides(aa: str, cfg: GiapConfig) -> list[tuple[int, float]]:
    """(length, GRAVY) of the distinct 0-missed-cleavage peptides in ``observable_len``."""
    return [(len(p), _gravy(p)) for p in _observable_seqs(aa, cfg)]


def _observable_seqs(aa: str, cfg: GiapConfig) -> list[str]:
    lo, hi = cfg.abundance.observable_len
    return list(dict.fromkeys(p for p in digest(aa, cfg.abundance.enzyme) if lo <= len(p) <= hi))


def _eta(coef: Mapping[str, float], x: Mapping[str, float]) -> float:
    return coef["intercept"] + sum(v * x[k] for k, v in coef.items() if k != "intercept")


def _cov_value(name: str, value: float) -> float:
    """M4-13 ``transform_covariates`` for one value (log1p / log; NaN stays NaN)."""
    if not math.isfinite(value):
        return math.nan
    if name == "n_observable_peptides":
        return math.log1p(value) if value > -1 else math.nan
    if name == "length_aa":
        return math.log(value) if value > 0 else math.nan
    return value


def p_detect(
    model: PMissModel,
    a: float,
    peps: Sequence[tuple[int, float]],
    *,
    cov: Mapping[str, float] | None = None,
) -> float:
    """P(gene called | abundance ``a``): Poisson-binomial P(>= k_min of ``peps`` seen) for the
    peptide model, the per-run probability for the replicate model; NaN if not usable."""
    if not model.usable or not math.isfinite(a):
        return math.nan
    if model.kind == "replicate":
        x = {"abundance_log2": a}
        for k in model.coef:
            if k not in ("intercept", "abundance_log2"):
                x[k] = _cov_value(k, float((cov or {}).get(k, math.nan)))
        return float(expit(_eta(model.coef, x)))
    if not peps:
        return math.nan
    p = [expit(_eta(model.coef, {"abundance_log2": a, "length": n, "gravy": g})) for n, g in peps]
    return min(1.0, max(0.0, poisson_binomial_sf(p, model.k_min)))  # pmf-sum rounding


def p_miss(
    model: PMissModel,
    a_v: float,
    ratios: ArrayLike,
    peps: Sequence[tuple[int, float]],
    *,
    cov: Mapping[str, float] | None = None,
) -> float:
    """mean over the oriented log-ratios r of 1 - p_detect(a_v + r); NaN for empty ``ratios``."""
    r = np.asarray(ratios, dtype=float).ravel()
    if r.size == 0:
        return math.nan
    return float(np.mean([1.0 - p_detect(model, a_v + x, peps, cov=cov) for x in r]))


def model_for_tu(model: PMissModel, tu_id: str) -> PMissModel:
    """The out-of-fold model of ``tu_id``; ``model`` itself for a TU that trained nothing."""
    j = model.fold_of_tu.get(tu_id)
    if j is None or not model.fold_coef:
        return model
    return dataclasses.replace(model, coef=model.fold_coef[j])


def p_miss_for_genes(
    model: PMissModel, genes: pd.DataFrame, cfg: GiapConfig
) -> Callable[[str, float, ArrayLike], float]:
    """``p_miss(u, a_v, ratios)`` for the genes of ``genes`` (``locus_tag``, ``tu_id``, ``aa``):
    the out-of-fold model of u's TU (:func:`model_for_tu`) on u's :func:`observable_peptides`
    (the M6-05 ``PMissModel`` callable, M6-05H). NaN for an unknown gene or one without a
    protein sequence; ``observable_peptides`` is computed once per gene."""
    g = genes.drop_duplicates("locus_tag").set_index("locus_tag")
    tu = g["tu_id"] if "tu_id" in g.columns else pd.Series(None, index=g.index, dtype=object)
    aa = g["aa"] if "aa" in g.columns else pd.Series(None, index=g.index, dtype=object)
    peps: dict[str, list[tuple[int, float]]] = {}

    def fn(u: str, a_v: float, ratios: ArrayLike) -> float:
        seq = aa.get(u)
        if not isinstance(seq, str):
            return math.nan
        if u not in peps:
            peps[u] = observable_peptides(seq, cfg)
        t = tu.get(u)
        m = model_for_tu(model, str(t)) if isinstance(t, str) else model
        return p_miss(m, a_v, ratios, peps[u])

    return fn


# ------------------------------------------------------------------------------ likelihoods

# A likelihood maps the linear predictor of every row to (log-lik, d log-lik / d eta).
LogLik = Callable[[Array], tuple[float, Array]]


def _states(p: Array, k: int) -> Array:
    """Forward counts capped at k (state k = ">= k"): F[:, j] = law of peptides < j."""
    g, n = p.shape
    f = np.zeros((g, n + 1, k + 1))
    f[:, 0, 0] = 1.0
    for j in range(n):
        q, prev = p[:, j : j + 1], f[:, j]
        new = prev * (1.0 - q)
        new[:, 1:] += prev[:, :-1] * q
        new[:, k] += prev[:, k] * q[:, 0]
        f[:, j + 1] = new
    return f


def _peptide_loglik(o: Array, mask: npt.NDArray[np.bool_], k: int) -> LogLik:
    """Sum_g [sum_j (o eta - log(1 + e^eta)) - log P(X_g >= k)] on padded (gene, peptide)."""

    def ll(eta: Array) -> tuple[float, Array]:
        eta = np.where(mask, eta, -np.inf)
        p = expit(eta)
        fwd = _states(p, k)
        bwd = _states(p[:, ::-1], k)[:, ::-1]  # bwd[:, j] = law of peptides >= j
        s = fwd[:, -1, k]
        cap = np.minimum(np.add.outer(np.arange(k + 1), np.arange(k + 1)), k) >= k - 1
        loo = np.einsum("gja,gjb,ab->gj", fwd[:, :-1], bwd[:, 1:], cap.astype(float))
        with np.errstate(divide="ignore", invalid="ignore"):
            cond = np.where(mask, p * loo / s[:, None], 0.0)  # P(O_j = 1 | X >= k)
            val = np.sum(np.where(mask, o * eta + log_expit(-eta), 0.0)) - np.sum(np.log(s))
        return float(val), o - cond

    return ll


def _replicate_loglik(s: Array, r: int) -> LogLik:
    """Sum_g [S log q + (R - S) log(1 - q) - log(1 - (1 - q)^R - q^R)]."""

    def ll(eta: Array) -> tuple[float, Array]:
        lq, l1q = log_expit(eta), log_expit(-eta)
        d = -np.expm1(r * lq) - np.exp(r * l1q)  # 1 - q^R - (1 - q)^R
        q = np.exp(lq)
        dd = r * q * (1.0 - q) * (np.exp((r - 1) * l1q) - np.exp((r - 1) * lq))
        with np.errstate(divide="ignore", invalid="ignore"):
            val = float(np.sum(s * lq + (r - s) * l1q - np.log(d)))
            return val, s - r * q - dd / d

    return ll


def _newton(x: Array, ll: LogLik) -> Array | None:
    """Unpenalised MLE of the standardised design ``x`` (rows x params, padded rows allowed):
    damped Newton with a central-difference Hessian of the analytic score. None when the
    iterate diverges (|coef| > COEF_BOUND), is non-finite or does not converge."""
    nb = x.shape[-1]

    def grad(b: Array) -> tuple[float, Array]:
        val, g = ll(x @ b)
        return val, np.tensordot(np.where(np.isfinite(g), g, 0.0), x, axes=g.ndim)

    beta = np.zeros(nb)
    val, g = grad(beta)
    for _ in range(MAX_ITER):
        h = np.empty((nb, nb))
        for i in range(nb):
            e = np.zeros(nb)
            e[i] = 1e-5
            h[i] = (grad(beta + e)[1] - grad(beta - e)[1]) / 2e-5
        neg_h = -(h + h.T) / 2
        try:
            step = np.linalg.solve(neg_h, g)
        except np.linalg.LinAlgError:
            return None
        full = step
        for _ in range(60):
            cand = beta + step
            v_c, g_c = grad(cand)
            if np.isfinite(v_c) and v_c >= val - 1e-12 * max(1.0, abs(val)):
                break
            step = step / 2
        else:
            return None
        beta, val, g = cand, v_c, g_c
        if not np.all(np.isfinite(beta)) or np.max(np.abs(beta)) > COEF_BOUND:
            return None
        if np.max(np.abs(step)) < 1e-9:
            # a step halved to nothing is not convergence: the Newton step itself must vanish at
            # a local maximum (-H positive definite), else the point is not an MLE
            return beta if np.max(np.abs(full)) < 1e-6 and _pos_def(neg_h) else None
    return None


def _pos_def(m: Array) -> bool:
    try:
        np.linalg.cholesky(m)
    except np.linalg.LinAlgError:
        return False
    return True


def _standardised_fit(
    z: Array, valid: npt.NDArray[np.bool_], ll: LogLik, names: Sequence[str]
) -> dict[str, float] | None:
    """Centre and scale the raw features ``z`` (..., f) over the ``valid`` cells, fit, and
    return the coefficients on the raw scale (keys ``intercept`` + ``names``)."""
    rows = z[valid]
    mu = rows.mean(axis=0)
    sd = rows.std(axis=0)
    sd = np.where(sd > 0, sd, 1.0)
    xs = np.where(valid[..., None], (z - mu) / sd, 0.0)
    x = np.concatenate([np.ones((*z.shape[:-1], 1)), xs], axis=-1)
    gamma = _newton(x, ll)
    if gamma is None:
        return None
    b = gamma[1:] / sd
    return {"intercept": float(gamma[0] - b @ mu)} | {
        n: float(v) for n, v in zip(names, b, strict=True)
    }


def _auc(y: Array, score: Array) -> float:
    pos = y > 0
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 == 0 or n0 == 0:
        return math.nan
    rk = rankdata(score)
    return float((rk[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def _cross_fit(
    tus: npt.NDArray[np.str_],
    fit: Callable[[npt.NDArray[np.bool_]], dict[str, float] | None],
    cfg: GiapConfig,
    rng: np.random.Generator,
) -> tuple[dict[str, float] | None, dict[str, int], list[dict[str, float]] | None, str | None]:
    """Full fit plus one fit per TU fold (trained on the other folds)."""
    k = int(cfg.detectability.folds)
    folds = tu_folds(list(tus), k, rng)
    if len(folds) < k:
        return None, folds, None, "too_few_tus"
    full = fit(np.ones(len(tus), dtype=bool))
    f = np.array([folds[t] for t in tus])
    per = [fit(f != j) for j in range(k)]
    if full is None or any(c is None for c in per):
        return None, folds, None, "separated_or_not_converged"
    return full, folds, [c for c in per if c is not None], None


def _unusable(
    kind: Kind, k_min: int, n_fit: int, reason: str, folds: Mapping[str, int] | None = None
) -> PMissModel:
    return PMissModel(kind, {}, k_min, n_fit, math.nan, False, dict(folds or {}), (), reason)


def _numeric(s: pd.Series) -> Array:
    return np.asarray(pd.to_numeric(s, errors="coerce"), dtype=float)


# ------------------------------------------------------------------------------ peptide model


def _k_min(genes: pd.DataFrame, cfg: GiapConfig) -> int:
    two = "detection_tier" in genes.columns and bool(
        (genes["detection_tier"].astype(object) == "two_peptide").any()
    )
    return int(cfg.detection.min_peptides_sequences_only) if two else 1


def _assigned(hits: pd.DataFrame | None) -> set[tuple[str, str]]:
    """(gene, peptide) pairs of C11 (any route, any uniqueness; I and L equivalent)."""
    if hits is None or len(hits) == 0:
        return set()
    h = hits[hits["gene"].notna()]
    return {
        (str(g), str(p).replace("I", "L")) for g, p in zip(h["gene"], h["peptide"], strict=True)
    }


def fit_peptide_model(
    genes: pd.DataFrame,
    peptide_hits: pd.DataFrame | None,
    cfg: GiapConfig,
    rng: np.random.Generator,
) -> PMissModel:
    """Conditional-likelihood peptide model, cross-fitted by TU folds (spec 2; decision 5)."""
    k = _k_min(genes, cfg)
    seen = _assigned(peptide_hits)
    # training = finite abundance and >= k_min observable peptides seen, whatever the gene-level
    # ``detected`` call says: the likelihood conditions on X_obs >= k_min, and selecting on the
    # FDR call (a strict subset in the mayu/decoy_picked tiers) biases beta_a (review of #157)
    det = genes[np.isfinite(_numeric(genes["abundance_log2"]))]
    det = det.sort_values("locus_tag", kind="stable")
    per_gene: list[tuple[str, float, list[tuple[int, float]], list[bool]]] = []
    for tag, tu, aa, a in zip(det["locus_tag"], det["tu_id"], det["aa"], det["abundance_log2"],
                              strict=True):  # fmt: skip
        seqs = _observable_seqs(str(aa), cfg) if isinstance(aa, str) else []
        o = [(str(tag), s.replace("I", "L")) in seen for s in seqs]
        if sum(o) >= k:
            per_gene.append((str(tu), float(a), [(len(s), _gravy(s)) for s in seqs], o))
    n_fit = sum(len(o) for *_, o in per_gene)
    if n_fit < MIN_ROWS:
        return _unusable("peptide", k, n_fit, "too_few_rows")
    obs = np.concatenate([np.asarray(o, float) for *_, o in per_gene])
    if obs.min() == obs.max():
        return _unusable("peptide", k, n_fit, "one_class")
    width = max(len(o) for *_, o in per_gene)
    mask = np.zeros((len(per_gene), width), dtype=bool)
    z = np.zeros((len(per_gene), width, 3))
    o_pad = np.zeros((len(per_gene), width))
    for i, (_, a, peps, o) in enumerate(per_gene):
        n = len(o)
        mask[i, :n] = True
        z[i, :n] = [(a, ln, g) for ln, g in peps]
        o_pad[i, :n] = o
    tus = np.array([t for t, *_ in per_gene])

    def fit(rows: npt.NDArray[np.bool_]) -> dict[str, float] | None:
        m = mask[rows]
        return _standardised_fit(z[rows], m, _peptide_loglik(o_pad[rows], m, k), PEPTIDE_FEATURES)

    full, folds, per, reason = _cross_fit(tus, fit, cfg, rng)
    if full is None or per is None:
        return _unusable("peptide", k, n_fit, reason or "separated_or_not_converged", folds)
    beta = np.array([[c["intercept"], *(c[n] for n in PEPTIDE_FEATURES)] for c in per])
    fold = np.array([folds[t] for t in tus])
    eta = beta[fold, 0, None] + np.einsum("gjf,gf->gj", z, beta[fold, 1:])
    auc = _auc(o_pad[mask], eta[mask])
    return PMissModel("peptide", full, k, n_fit, auc, True, folds, tuple(per))


# ------------------------------------------------------------------------------ replicate model


def _runs_seen(
    genes: pd.DataFrame, runs: pd.DataFrame, hits: pd.DataFrame | None
) -> tuple[pd.DataFrame, int]:
    """(gene x run) presence from C16 joined to C11 assignments, and the number of runs R."""
    names = sorted({str(r) for r in runs["run"]})
    pairs = pd.DataFrame(sorted(_assigned(hits)), columns=["gene", "peptide"])
    pr = pd.DataFrame({"peptide": runs["peptide"].astype(str).str.replace("I", "L"),
                       "run": runs["run"].astype(str)})  # fmt: skip
    m = pairs.merge(pr, on="peptide").drop_duplicates(["gene", "run"])
    cnt = m.groupby("gene")["run"].nunique()
    return cnt.reindex(genes["locus_tag"].astype(str)).fillna(0).to_frame("s"), len(names)


def fit_replicate_model(
    genes: pd.DataFrame,
    peptide_runs: pd.DataFrame | None,
    cfg: GiapConfig,
    rng: np.random.Generator,
    *,
    peptide_hits: pd.DataFrame | None,
) -> PMissModel:
    """Per-run detection on abundance + ``detectability.covariates`` (spec 4), fitted by the
    doubly truncated binomial on genes seen in 1..R-1 runs; R = 2 is ``replicate_two_runs``.
    Covariates absent from ``genes`` are not used; zero-variance ones are dropped."""
    r = 0 if peptide_runs is None else peptide_runs["run"].astype(str).nunique()
    if peptide_runs is None or r < 2:
        return _unusable("replicate", 1, 0, "too_few_runs")
    if r == 2:  # decided from the run count: S = 1 of 2 for every training gene, no slope
        return _unusable("replicate", 1, 0, "replicate_two_runs")
    counts, r = _runs_seen(genes, peptide_runs, peptide_hits)
    s_all = counts["s"].to_numpy(float)
    a_all = _numeric(genes["abundance_log2"])
    keep = (s_all >= 1) & (s_all <= r - 1) & np.isfinite(a_all)
    g = genes[keep].assign(_s=s_all[keep]).sort_values("locus_tag", kind="stable")
    n_fit = len(g) * r
    if n_fit < MIN_ROWS:
        return _unusable("replicate", 1, n_fit, "too_few_rows")
    names = [c for c in cfg.detectability.covariates if c in g.columns]
    cov = transform_covariates(g, names)
    names = [c for c in names if float(np.std(cov[c].to_numpy())) > 0]
    z = np.column_stack([g["abundance_log2"].to_numpy(float), cov[names].to_numpy(float)])
    s = g["_s"].to_numpy(float)
    feats = ["abundance_log2", *names]
    valid = np.ones(len(g), dtype=bool)

    def fit(rows: npt.NDArray[np.bool_]) -> dict[str, float] | None:
        return _standardised_fit(z[rows], valid[rows], _replicate_loglik(s[rows], r), feats)

    tus = g["tu_id"].astype(str).to_numpy()
    full, folds, per, reason = _cross_fit(tus, fit, cfg, rng)
    if full is None or per is None:
        return _unusable("replicate", 1, n_fit, reason or "separated_or_not_converged", folds)
    beta = np.array([[c["intercept"], *(c[n] for n in feats)] for c in per])
    fold = np.array([folds[t] for t in tus])
    q = beta[fold, 0] + np.einsum("gf,gf->g", z, beta[fold, 1:])
    y = np.concatenate([np.ones(int(x)) for x in s] + [np.zeros(int(r - x)) for x in s])
    sc = np.concatenate([np.full(int(x), v) for x, v in zip(s, q, strict=True)]
                        + [np.full(int(r - x), v) for x, v in zip(s, q, strict=True)])  # fmt: skip
    return PMissModel("replicate", full, 1, n_fit, _auc(y, sc), True, folds, tuple(per))


# ------------------------------------------------------------------------------ dispatcher


def _row(group: str, metric: str, value: float) -> dict[str, object]:
    nan = float("nan")
    return dict(zip(COLUMNS, (SECTION, "model", group, metric, float(value), nan, nan),
                    strict=True))  # fmt: skip


def fit_detection_limit(
    genes: pd.DataFrame,
    peptide_hits: pd.DataFrame | None,
    peptide_runs: pd.DataFrame | None,
    cfg: GiapConfig,
    rng: np.random.Generator,
) -> DetectionLimitFit:
    """The peptide model whenever peptides exist (Q3), the replicate model when C16 is given;
    in force: the peptide model if usable, else the replicate model if usable, else the
    (unusable) peptide model, else the replicate model. Both AUCs go to ``summary``."""
    has_peps = peptide_hits is not None and len(peptide_hits) > 0
    pep = fit_peptide_model(genes, peptide_hits, cfg, rng) if has_peps else None
    rep = None
    if peptide_runs is not None:
        rep = fit_replicate_model(genes, peptide_runs, cfg, rng, peptide_hits=peptide_hits)
    fitted = [m for m in (pep, rep) if m is not None]
    usable = [m for m in fitted if m.usable]
    model = (usable or fitted or [_unusable("peptide", 1, 0, "no_input")])[0]
    rows = [_row(m.kind, met, v) for m in (pep, rep) if m is not None
            for met, v in (("auc", m.auc), ("n_fit", m.n_fit), ("usable", m.usable))]  # fmt: skip
    rows.append(_row(model.kind, "in_force", 1.0))
    summary = pd.DataFrame(rows, columns=COLUMNS).astype(
        {"value": "float64", "ci_low": "float64", "ci_high": "float64"}
    )
    return DetectionLimitFit(model, pep, rep, summary)
