"""UniOP / TerminatorNet cut-offs and the d0 / s refit on K-12 truth (issue M6-12).

GIAP_SPEC.md §15.4 steps 3-5, §7.1.3, §7.1.5, §7.4.6; spec docs/specs/M6-12-uniop-terminator-
cutoffs.md (lead decisions 1, 6, 7); methods docs/methods/operons.md "K-12 cut-off calibration
(M6-12)". Pure functions: no I/O; randomness only from the ``rng`` passed in (callers use
``giap.rng.for_step(STEP, seed)``); 0-based half-open coordinates. An item is called at t iff
score >= t.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import NamedTuple

import numpy as np
import numpy.typing as npt
import pandas as pd
from numpy.typing import ArrayLike
from scipy.special import expit
from scipy.stats import norm, rankdata

from giap.annotate.terminators import DOWNSTREAM_NO_NEXT_BP, _window
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene
from giap.stats.bootstrap import BootCI, tu_block_bootstrap

__all__ = [
    "D0_DEFAULT",
    "DEFAULT_B",
    "DEFAULT_SEED",
    "MIN_LABELLED",
    "STEP",
    "S_DEFAULT",
    "Cutoff",
    "DistanceFit",
    "Window",
    "auroc",
    "auroc_ci",
    "auroc_diff_ci",
    "choose_cutoff",
    "coincidence_rate",
    "distance_fit_ci",
    "distance_log_loss",
    "downstream_windows",
    "fit_distance_model",
    "fixed_cutoff_ci",
    "interval_gap",
    "precision_at",
    "precision_cv",
    "precision_lower_bound",
    "recall_at",
    "require_labelled",
    "terminator_match",
    "threshold_boot_q",
    "threshold_curve",
    "wilson_ci",
    "window_of",
]

STEP = "m6_12_cutoffs"
DEFAULT_SEED = 20261014
DEFAULT_B = 2000
MIN_LABELLED = 500
D0_DEFAULT = 50.0
S_DEFAULT = 20.0

Floats = npt.NDArray[np.float64]
Ints = npt.NDArray[np.int64]


def _xy(x: ArrayLike, y: ArrayLike) -> tuple[Floats, Ints]:
    """Validated (finite float, 0/1 int) arrays of equal, non-zero length."""
    xs = np.asarray(x, dtype=np.float64).ravel()
    ys = np.asarray(y).ravel()
    if xs.size != ys.size or xs.size == 0:
        raise ValueError(f"need equal non-empty lengths, got {xs.size} and {ys.size}")
    if not np.all(np.isfinite(xs)):
        raise ValueError("values must be finite")
    if not np.all((ys == 0) | (ys == 1)):
        raise ValueError("labels must be 0 or 1")
    return xs, ys.astype(np.int64)


def require_labelled(y: ArrayLike, min_n: int = MIN_LABELLED) -> tuple[int, int]:
    """(P, N); ``ValueError`` below ``min_n`` items or with one class (no calibration claimed)."""
    ys = np.asarray(y).ravel()
    p, n = int((ys == 1).sum()), int((ys == 0).sum())
    if ys.size < min_n or p == 0 or n == 0:
        raise ValueError(f"need >= {min_n} labelled items of both classes, got P={p}, N={n}")
    return p, n


def auroc(scores: ArrayLike, y: ArrayLike) -> float:
    """Mann-Whitney AUROC through average ranks (ties count 1/2); NaN when P = 0 or N = 0."""
    s, ys = _xy(scores, y)
    p = int(ys.sum())
    n = ys.size - p
    if p == 0 or n == 0:
        return math.nan
    r = rankdata(s, method="average")
    return float((r[ys == 1].sum() - p * (p + 1) / 2) / (p * n))


def threshold_curve(scores: ArrayLike, y: ArrayLike) -> pd.DataFrame:
    """Distinct scores descending with cumulative integer ``tp`` and ``fp`` at score >= t."""
    s, ys = _xy(scores, y)
    u, inv = np.unique(s, return_inverse=True)
    pos = np.bincount(inv[ys == 1], minlength=u.size)[::-1]
    tot = np.bincount(inv, minlength=u.size)[::-1]
    tp = np.cumsum(pos).astype(np.int64)
    return pd.DataFrame({"threshold": u[::-1], "tp": tp, "fp": np.cumsum(tot) - tp})


@dataclass(frozen=True, slots=True)
class Cutoff:
    """The chosen cut-off; ``literal_threshold`` is the smallest t reaching the target."""

    threshold: float
    rule: str
    precision: float
    recall: float
    f1: float
    tp: int
    fp: int
    literal_threshold: float | None


def choose_cutoff(curve: pd.DataFrame, precision: float = 0.90) -> Cutoff:
    """Stable rule (lead decision 1): the last candidate before the first one failing
    ``b*TP >= a*(TP+FP)`` (target a/b), descending; if the top fails, max F1 (exact), smallest t.
    """
    target = Fraction(str(precision))
    a, b = target.numerator, target.denominator
    th = [float(v) for v in curve["threshold"]]
    tp = [int(v) for v in curve["tp"]]
    fp = [int(v) for v in curve["fp"]]
    p = tp[-1]
    ok = [b * t >= a * (t + f) for t, f in zip(tp, fp, strict=True)]
    literal = next((th[k] for k in range(len(th) - 1, -1, -1) if ok[k]), None)
    n_ok = next((k for k, good in enumerate(ok) if not good), len(ok))
    if n_ok:
        k, rule = n_ok - 1, f"precision_{precision:.2f}"
    else:
        f1 = [Fraction(2 * t, t + f + p) for t, f in zip(tp, fp, strict=True)]
        best = max(f1)
        k, rule = max(j for j, v in enumerate(f1) if v == best), "max_f1"
    t, f = tp[k], fp[k]
    recall = t / p if p else math.nan
    return Cutoff(th[k], rule, t / (t + f), recall, 2 * t / (t + f + p), t, f, literal)


def precision_at(scores: ArrayLike, y: ArrayLike, threshold: float) -> float:
    """Share of positives among items with score >= threshold; NaN when none is called."""
    s, ys = _xy(scores, y)
    called = s >= threshold
    return float(ys[called].mean()) if called.any() else math.nan


def recall_at(scores: ArrayLike, y: ArrayLike, threshold: float) -> float:
    """Share of positives with score >= threshold; NaN when P = 0."""
    s, ys = _xy(scores, y)
    p = int(ys.sum())
    return int(ys[s >= threshold].sum()) / p if p else math.nan


# --- directon bootstrap (methods "Bootstrap streams": one tu_block_bootstrap call per CI) -----

_Stat = Callable[[pd.DataFrame], float]


def _frame(cluster: ArrayLike, y: ArrayLike, **cols: ArrayLike) -> pd.DataFrame:
    data = {k: np.asarray(v, dtype=np.float64).ravel() for k, v in cols.items()}
    return pd.DataFrame({**data, "y": np.asarray(y).ravel(), "cluster": np.asarray(cluster)})


def _boot(df: pd.DataFrame, stat: _Stat, B: int, rng: np.random.Generator, lv: float) -> BootCI:
    return tu_block_bootstrap(df, "cluster", stat, B, rng, lv)


def auroc_ci(
    scores: ArrayLike,
    y: ArrayLike,
    cluster: ArrayLike,
    *,
    rng: np.random.Generator,
    B: int = DEFAULT_B,
    level: float = 0.95,
) -> BootCI:
    """AUROC with its cluster (directon) percentile CI; one-class replicates are dropped."""
    df = _frame(cluster, y, x=scores)
    return _boot(df, lambda d: auroc(d["x"], d["y"]), B, rng, level)


def fixed_cutoff_ci(
    scores: ArrayLike,
    y: ArrayLike,
    cluster: ArrayLike,
    threshold: float,
    *,
    rng: np.random.Generator,
    B: int = DEFAULT_B,
    level: float = 0.95,
) -> tuple[BootCI, BootCI]:
    """(precision, recall) CIs at a fixed ``threshold`` (not re-selected); precision drawn first."""
    df = _frame(cluster, y, x=scores)
    prec = _boot(df, lambda d: precision_at(d["x"], d["y"], threshold), B, rng, level)
    rec = _boot(df, lambda d: recall_at(d["x"], d["y"], threshold), B, rng, level)
    return prec, rec


def auroc_diff_ci(
    score_a: ArrayLike,
    score_b: ArrayLike,
    y: ArrayLike,
    cluster: ArrayLike,
    *,
    rng: np.random.Generator,
    B: int = DEFAULT_B,
    level: float = 0.95,
) -> BootCI:
    """Paired AUROC(a) - AUROC(b), both on the same replicate (R4)."""
    df = _frame(cluster, y, a=score_a, b=score_b)
    return _boot(df, lambda d: auroc(d["a"], d["y"]) - auroc(d["b"], d["y"]), B, rng, level)


def threshold_boot_q(
    scores: ArrayLike,
    y: ArrayLike,
    cluster: ArrayLike,
    *,
    rng: np.random.Generator,
    B: int = DEFAULT_B,
    precision: float = 0.90,
) -> tuple[float, float, float]:
    """2.5 / 50 / 97.5 % linear quantiles of t* re-selected per replicate (R1); one-class
    replicates are NaN and dropped; fewer than B/2 valid gives NaNs."""
    seen: list[float] = []

    def stat(d: pd.DataFrame) -> float:
        ys = d["y"].to_numpy()
        one_class = bool(ys.min() == ys.max())
        cut = None if one_class else choose_cutoff(threshold_curve(d["x"], ys), precision)
        seen.append(math.nan if cut is None else cut.threshold)
        return seen[-1]

    _boot(_frame(cluster, y, x=scores), stat, B, rng, 0.95)
    vals = np.asarray(seen[1:])
    vals = vals[~np.isnan(vals)]
    if vals.size == 0 or vals.size < B / 2:
        return (math.nan, math.nan, math.nan)
    lo, mid, hi = np.quantile(vals, [0.025, 0.5, 0.975], method="linear")
    return (float(lo), float(mid), float(hi))


def precision_cv(
    scores: ArrayLike,
    y: ArrayLike,
    cluster: ArrayLike,
    *,
    rng: np.random.Generator,
    k: int = 10,
    precision: float = 0.90,
) -> float:
    """Pooled k-fold cluster CV precision of the whole rule (R1): t* chosen on k-1 folds, TP/FP
    counted on the held-out fold. ``ValueError`` with fewer than ``k`` clusters."""
    s, ys = _xy(scores, y)
    labels, inv = np.unique(np.asarray(cluster).ravel(), return_inverse=True)
    if labels.size < k or inv.size != s.size:
        raise ValueError(f"need >= {k} clusters, one label per item; got {labels.size}")
    fold = (rng.permutation(labels.size) % k)[inv]
    tp = fp = 0
    for j in range(k):
        train = fold != j
        t = choose_cutoff(threshold_curve(s[train], ys[train]), precision).threshold
        held = ys[~train & (s >= t)]
        tp, fp = tp + int(held.sum()), fp + int((held == 0).sum())
    return tp / (tp + fp) if tp + fp else math.nan


def wilson_ci(k: int, n: int, level: float = 0.95) -> tuple[float, float]:
    """Wilson score interval for k of n; exact 0 / 1 endpoints at k = 0 / n; NaNs for n = 0."""
    if n == 0:
        return (math.nan, math.nan)
    z = float(norm.ppf((1.0 + level) / 2.0))
    ph, z2 = k / n, z * z
    centre = (ph + z2 / (2 * n)) / (1 + z2 / n)
    half = z * math.sqrt(ph * (1 - ph) / n + z2 / (4 * n * n)) / (1 + z2 / n)
    lo = 0.0 if k == 0 else max(0.0, centre - half)
    return (lo, 1.0 if k == n else min(1.0, centre + half))


# --- the d0 / s refit (methods B1) --------------------------------------------------------


@dataclass(frozen=True, slots=True)
class DistanceFit:
    """p(operon | d) = 1/(1 + exp((d - d0)/s)); a failure keeps (50, 20) and ``reason``."""

    d0: float
    s: float
    beta0: float
    beta1: float
    converged: bool
    n_iter: int
    reason: str | None

    @property
    def calibrated(self) -> bool:
        return self.converged


def _failed(beta: Floats | None, n_iter: int, reason: str) -> DistanceFit:
    b0, b1 = (math.nan, math.nan) if beta is None else (float(beta[0]), float(beta[1]))
    return DistanceFit(D0_DEFAULT, S_DEFAULT, b0, b1, False, n_iter, reason)


def fit_distance_model(
    d: ArrayLike, y: ArrayLike, *, max_iter: int = 100, tol: float = 1e-10
) -> DistanceFit:
    """Unpenalised logistic MLE of y on d: exact separation check, start at (logit ybar, 0),
    Newton with step-halving (<= 30). Never raises on degenerate data; ``ValueError`` only on
    invalid input (lengths, non-finite d, labels outside {0, 1})."""
    x, ys = _xy(d, y)
    pos, neg = x[ys == 1], x[ys == 0]
    if pos.size == 0 or neg.size == 0:
        return _failed(None, 0, "one_class")
    if not (pos.max() > neg.min() and neg.max() > pos.min()):
        return _failed(None, 0, "separation")
    X = np.column_stack((np.ones_like(x), x))
    yf = ys.astype(np.float64)

    def loglik(b: Floats) -> float:
        eta = X @ b
        return float(np.sum(yf * eta - np.logaddexp(0.0, eta)))

    beta = np.array([math.log(pos.size / neg.size), 0.0])
    with np.errstate(all="ignore"):
        for it in range(1, max_iter + 1):
            p = expit(X @ beta)
            try:
                step = np.linalg.solve(X.T @ (X * (p * (1 - p))[:, None]), X.T @ (yf - p))
            except np.linalg.LinAlgError:
                return _failed(beta, it, "singular")
            cur, h = loglik(beta), 1.0
            while loglik(beta + h * step) < cur and h > 2.0**-30:
                h /= 2.0
            beta = beta + h * step
            if not np.all(np.isfinite(beta)):
                return _failed(beta, it, "non_finite")
            if float(np.max(np.abs(h * step))) < tol:
                break
        else:
            return _failed(beta, max_iter, "no_convergence")
        b0, b1 = float(beta[0]), float(beta[1])
        if b1 >= 0:
            return _failed(beta, it, "beta1_nonnegative")
        if abs(b1) > 1:  # s < 1 bp: saturated weights under quasi-separation (backstop)
            return _failed(beta, it, "steep")
        d0 = -b0 / b1
    if not math.isfinite(d0):
        return _failed(beta, it, "non_finite")
    return DistanceFit(d0, -1.0 / b1, b0, b1, True, it, None)


def distance_fit_ci(
    d: ArrayLike,
    y: ArrayLike,
    cluster: ArrayLike,
    *,
    rng: np.random.Generator,
    B: int = DEFAULT_B,
    level: float = 0.95,
) -> tuple[BootCI, BootCI]:
    """(d0, s) cluster-bootstrap CIs, d0 drawn first, separate draws; failed fits are NaN."""
    df = _frame(cluster, y, x=d)

    def refit(attr: str) -> _Stat:
        def stat(r: pd.DataFrame) -> float:
            fit = fit_distance_model(r["x"].to_numpy(), r["y"].to_numpy())
            return float(getattr(fit, attr)) if fit.converged else math.nan

        return stat

    return _boot(df, refit("d0"), B, rng, level), _boot(df, refit("s"), B, rng, level)


def distance_log_loss(d: ArrayLike, y: ArrayLike, d0: float, s: float) -> float:
    """Mean negative log-likelihood (nats) of the distance model, stable in the tails."""
    x, ys = _xy(d, y)
    z = (x - d0) / s  # -log p = log(1 + e^z); -log(1 - p) = log(1 + e^z) - z
    return float(np.mean(np.logaddexp(0.0, z) - (1 - ys) * z))


# --- TerminatorNet windows, matching and the coincidence null (methods B2) ----------------


class Window(NamedTuple):
    """A gene's downstream window ``[lo, hi)``, unrolled across the origin; empty if lo >= hi."""

    gene: str
    strand: str
    lo: int
    hi: int


def downstream_windows(
    genes: Sequence[Gene], replicon: Replicon, slack_bp: int = 50
) -> list[Window]:
    """One window per gene on ``replicon``, in genome order: from the stop moved ``slack_bp``
    into the gene to the genome-order neighbour on either strand (``_window``); a linear-end gene
    gets ``DOWNSTREAM_NO_NEXT_BP`` past its stop, clipped to the contig."""
    gs = sorted(
        (g for g in genes if g.interval.contig == replicon.id),
        key=lambda g: (g.interval.start, g.interval.end, g.locus_tag),
    )
    circ, n, length = replicon.topology == "circular", len(gs), replicon.length
    out = []
    for i, g in enumerate(gs):
        iv = g.interval
        j = i + 1 if iv.strand == "+" else i - 1
        if 0 <= j < n or circ:
            lo, hi = _window(g, gs[j % n], replicon, slack_bp)
            hi += length if n == 1 else 0  # a lone gene on a circle is its own neighbour
        elif iv.strand == "-":
            lo, hi = max(0, iv.start - DOWNSTREAM_NO_NEXT_BP), iv.start + slack_bp
        else:
            lo, hi = iv.end - slack_bp, min(length, iv.end + DOWNSTREAM_NO_NEXT_BP)
        out.append(Window(g.locus_tag, iv.strand, lo, hi))
    return out


def _shifts(replicon: Replicon) -> tuple[int, ...]:
    return (0, -replicon.length, replicon.length) if replicon.topology == "circular" else (0,)


def _containing(
    iv: Interval, windows: Sequence[Window], replicon: Replicon
) -> tuple[Window, int] | None:
    """The first same-strand window wholly holding ``iv`` shifted by k, with k."""
    for w in windows:
        if w.strand == iv.strand:
            for k in _shifts(replicon):
                if w.lo <= iv.start + k and iv.end + k <= w.hi:
                    return w, k
    return None


def window_of(iv: Interval, windows: Sequence[Window], replicon: Replicon) -> str | None:
    """Locus tag of the first same-strand window containing ``iv`` (any shift), else None."""
    hit = _containing(iv, windows, replicon)
    return None if hit is None else hit[0].gene


def interval_gap(a: Interval, b: Interval, replicon: Replicon) -> int:
    """min over shifts of max(0, b.start - a.end, a.start - b.end); overlap is 0."""
    return min(max(0, b.start + k - a.end, a.start - b.end - k) for k in _shifts(replicon))


def terminator_match(
    iv: Interval, truth: Sequence[Interval], replicon: Replicon, tol_bp: int = 50
) -> int | None:
    """Index of the nearest same-contig, same-strand truth with gap <= ``tol_bp`` (first on tie)."""
    best: tuple[int, int] | None = None
    for j, t in enumerate(truth):
        if t.contig == iv.contig and t.strand == iv.strand:
            gap = interval_gap(iv, t, replicon)
            if gap <= tol_bp and (best is None or gap < best[0]):
                best = (gap, j)
    return None if best is None else best[1]


def coincidence_rate(
    calls: Sequence[Interval],
    confidence: ArrayLike,
    windows: Sequence[Window],
    truth: Sequence[Interval],
    replicon: Replicon,
    threshold: float,
    *,
    rng: np.random.Generator,
    R: int = 200,
    tol_bp: int = 50,
) -> float:
    """B2 placement null: the mean over ``R`` placements of the matched share of eligible calls
    with confidence >= ``threshold``. Each call moves to a uniform same-strand window long enough
    to hold it, keeping its length and its offset from the window's 3' side (lo for +, hi for -).
    NaN when no call is at or above the threshold; ``ValueError`` on an ineligible call."""
    conf = np.asarray(confidence, dtype=np.float64).ravel()
    if conf.size != len(calls):
        raise ValueError(f"{len(calls)} calls but {conf.size} confidences")
    kept = [c for c, p in zip(calls, conf, strict=True) if p >= threshold]
    if not kept:
        return math.nan
    shifts = np.asarray(_shifts(replicon))[:, None, None]
    matched = np.zeros(R, dtype=np.int64)
    for strand in ("+", "-"):
        cs = [c for c in kept if c.strand == strand]
        if not cs:
            continue
        ws = sorted((w for w in windows if w.strand == strand), key=lambda w: w.lo - w.hi)
        lo, hi = np.array([w.lo for w in ws]), np.array([w.hi for w in ws])
        ln = np.array([c.end - c.start for c in cs])
        off = np.array([_offset(c, windows, replicon) for c in cs])
        n_fit = np.searchsorted(lo - hi, -(off + ln), side="right")  # windows longest first
        pick = rng.integers(0, n_fit, size=(R, len(cs)))
        starts = lo[pick] + off if strand == "+" else hi[pick] - off - ln
        tt = [t for t in truth if t.contig == cs[0].contig and t.strand == strand]
        if not tt:
            continue
        ts, te = np.array([t.start for t in tt]), np.array([t.end for t in tt])
        for r in range(R):
            s, e = starts[r][:, None], (starts[r] + ln)[:, None]
            near = (ts + shifts - e <= tol_bp) & (s - te - shifts <= tol_bp)
            matched[r] += int(near.any(axis=(0, 2)).sum())
    return float(np.mean(matched / len(kept)))


def _offset(c: Interval, windows: Sequence[Window], replicon: Replicon) -> int:
    hit = _containing(c, windows, replicon)
    if hit is None:
        raise ValueError(f"call {c.contig}:{c.start}-{c.end}{c.strand} lies in no window")
    w, k = hit
    return c.start + k - w.lo if c.strand == "+" else w.hi - (c.end + k)


def precision_lower_bound(obs: float, c: float) -> float:
    """Coincidence-corrected bound max(0, (obs - c)/(1 - c)); NaN when c >= 1 or non-finite."""
    if not (math.isfinite(obs) and math.isfinite(c)) or c >= 1.0:
        return math.nan
    return max(0.0, (obs - c) / (1.0 - c))
