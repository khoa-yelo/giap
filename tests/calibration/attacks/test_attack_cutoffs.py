"""Attack nulls for the M6-12 cut-off statistics (GIAP_SPEC.md §8.5; methods
docs/methods/operons.md "K-12 cut-off calibration (M6-12)" > "Attack nulls owed"; spec lead
decision 6).

1. *Directon-clustered null* (nightly): K-12-like directon sizes (1-32 pairs, median 2), a shared
   directon random effect on the scores (ICC 0.3) and directon-level label prevalence, drawn
   independently, so the true AUROC is 0.5. Over R = 500 the 95 % directon CI must cover 0.5 in
   at least 93 % of replicates (the methods' criterion). The single-pair CI's coverage is printed.
2. *Coincidence null* (B2): calls placed by the null itself give an observed match share equal to
   the coincidence rate c within Monte Carlo error, and a corrected bound near 0; calls placed
   exactly on truth give a small c and a bound equal to the observed precision.
3. *Top-group spike*: a top tie group at precision 0.85 above a broad region at 0.95 sends the
   stable rule to max-F1, and the literal value is reported beside it.
4. *Hard refit start* (B1), the separated and quasi-separated cases and random hostile inputs:
   none raises; failures keep (50, 20) with ``calibrated`` false.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
import pytest

from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene
from giap.rng import for_step
from giap.stats.cutoffs import (
    DEFAULT_SEED,
    STEP,
    auroc_ci,
    choose_cutoff,
    coincidence_rate,
    downstream_windows,
    fit_distance_model,
    precision_lower_bound,
    terminator_match,
    threshold_curve,
    window_of,
)

pytestmark = pytest.mark.calsuite("regions")

COVERAGE_MIN = 0.93  # methods "Attack nulls owed" item 1


# --- 1. directon-clustered null --------------------------------------------------------------


def k12_like_null(g: np.random.Generator, n_dir: int = 200) -> tuple[np.ndarray, ...]:
    sizes = np.minimum(g.geometric(0.4, n_dir), 32)  # median 2, max 32
    cl = np.repeat(np.arange(n_dir), sizes)
    pi = g.beta(2.0, 1.2, n_dir)[cl]  # directon-level prevalence, mean ~0.63
    y = (g.random(cl.size) < pi).astype(int)
    s = g.normal(0.0, math.sqrt(0.3), n_dir)[cl] + g.normal(0.0, math.sqrt(0.7), cl.size)
    return s, y, cl


@pytest.mark.calibration
@pytest.mark.calibration_nightly
def test_directon_clustered_null_auroc_ci_coverage(tolerances: dict[str, object]) -> None:
    seeds = tolerances["seeds"]
    assert isinstance(seeds, dict)
    gen = for_step("calibration.m6_12.attack_directon_null", int(seeds["calibration_base"]))
    r, b = 500, 500
    cov_dir = cov_iid = 0
    for _ in range(r):
        s, y, cl = k12_like_null(gen)
        a = auroc_ci(s, y, cl, rng=gen, B=b)
        i = auroc_ci(s, y, np.arange(y.size), rng=gen, B=b)
        cov_dir += a.ci_low <= 0.5 <= a.ci_high
        cov_iid += i.ci_low <= 0.5 <= i.ci_high
    print(f"M6-12 directon null (R={r}, B={b}): directon {cov_dir / r:.3f}, iid {cov_iid / r:.3f}")
    assert cov_dir / r >= COVERAGE_MIN


# --- 2. coincidence null (B2) ----------------------------------------------------------------

N_GENES = 150
SPACING = 3400  # gene length 900, so every window is 3400 - 900 + 50 = 2550 bp
RING = Replicon("ring", N_GENES * SPACING, "circular")


def ring_genes() -> list[Gene]:
    return [
        Gene(
            locus_tag=f"r{i:03d}",
            feature_type="cds",
            interval=Interval("ring", i * SPACING, i * SPACING + 900, "+"),
            aa="M" * 299,
            nt=None,
            start_codon="ATG",
        )
        for i in range(N_GENES)
    ]


def _at(lo: int, off: int, ln: int = 30) -> Interval:
    start = (lo + off) % RING.length
    return Interval("ring", start, start + ln, "+")


def test_coincidence_null_calls_placed_by_the_null() -> None:
    g = for_step("calibration.m6_12.attack_coincidence", DEFAULT_SEED)
    w = downstream_windows(ring_genes(), RING)
    assert {x.hi - x.lo for x in w} == {2550}
    truth = [_at(x.lo, int(min(abs(g.normal(0, 150)), 2000)), 20) for x in w]
    calls = []
    for _ in range(400):
        x = w[int(g.integers(0, len(w)))]
        calls.append(_at(x.lo, int(min(abs(g.normal(0, 150)), 2000))))
    assert all(window_of(c, w, RING) is not None for c in calls)
    obs = float(np.mean([terminator_match(c, truth, RING) is not None for c in calls]))
    c = coincidence_rate(calls, [0.9] * len(calls), w, truth, RING, 0.5, rng=g, R=100)
    bound = precision_lower_bound(obs, c)
    print(f"coincidence null: obs {obs:.3f}, c {c:.3f}, bound {bound:.3f}")
    assert 0.1 < c < 0.9  # offsets concentrate past the stop: coincidences are common
    assert abs(obs - c) <= 0.08
    assert bound <= 0.15


def test_coincidence_null_calls_on_truth() -> None:
    g = for_step("calibration.m6_12.attack_on_truth", DEFAULT_SEED)
    w = downstream_windows(ring_genes(), RING)
    truth = [_at(x.lo, int(g.integers(0, 2400))) for x in w]
    calls = list(truth)
    obs = float(np.mean([terminator_match(c, truth, RING) is not None for c in calls]))
    c = coincidence_rate(calls, [0.9] * len(calls), w, truth, RING, 0.5, rng=g)
    print(f"calls on truth: obs {obs:.3f}, c {c:.3f}")
    assert obs == 1.0
    assert c < 0.15
    assert precision_lower_bound(obs, c) == pytest.approx(obs)


# --- 3. top-group spike ----------------------------------------------------------------------


def test_top_group_spike_sends_the_rule_to_max_f1() -> None:
    g = np.random.default_rng(DEFAULT_SEED)
    scores = [1.0] * 20
    y = [1] * 17 + [0] * 3  # top tie group at 0.85
    broad = np.round(g.uniform(0.5, 0.99, 400), 4)
    yb = [1] * 380 + [0] * 20  # broad region at 0.95
    scores += broad.tolist() + np.round(g.uniform(0.0, 0.4, 300), 4).tolist()
    y += yb + [0] * 300
    curve = threshold_curve(scores, y)
    c = choose_cutoff(curve)
    assert c.rule == "max_f1"
    assert c.literal_threshold is not None and c.literal_threshold < 1.0
    p = sum(y)
    f1 = [Fraction(2 * t, t + f + p) for t, f in zip(curve["tp"], curve["fp"], strict=True)]
    best = max(f1)
    want = min(t for t, v in zip(curve["threshold"], f1, strict=True) if v == best)
    assert c.threshold == want


# --- 4. refit: hard start, separation, hostile inputs ----------------------------------------


def _grid(xs: range, n: int, d0: float, s: float) -> tuple[np.ndarray, np.ndarray]:
    d: list[int] = []
    y: list[int] = []
    for x in xs:
        k = round(n / (1.0 + math.exp((x - d0) / s)))
        d += [x] * n
        y += [1] * k + [0] * (n - k)
    return np.asarray(d), np.asarray(y)


def test_refit_hard_start_separated_and_quasi_never_raise() -> None:
    d, y = _grid(range(-20, 401, 20), 100, 100.0, 50.0)
    f = fit_distance_model(d, y)
    assert f.converged and abs(f.d0 - 100.14) <= 0.05 and abs(f.s - 50.06) <= 0.05
    d5, _ = _grid(range(-20, 201, 10), 1000, 40.0, 15.0)
    sep = (d5 < 60).astype(int)
    for dd, yy in ((d5, sep), (np.r_[d5, [60, 60, 60]], np.r_[sep, [1, 1, 1]])):
        f = fit_distance_model(dd, yy)
        assert not f.converged and not f.calibrated and (f.d0, f.s) == (50.0, 20.0)


def test_refit_random_hostile_inputs_never_raise() -> None:
    g = for_step("calibration.m6_12.attack_refit", DEFAULT_SEED)
    for _ in range(300):
        n = int(g.integers(2, 40))
        scale = float(10.0 ** g.integers(0, 6))
        d = np.round(g.normal(0.0, scale, n))
        y = g.integers(0, 2, n)
        f = fit_distance_model(d, y)
        if f.converged:
            assert math.isfinite(f.d0) and f.s >= 1.0 and f.beta1 < 0
        else:
            assert (f.d0, f.s) == (50.0, 20.0) and not f.calibrated


def test_rng_stream_name() -> None:
    assert STEP == "m6_12_cutoffs"
