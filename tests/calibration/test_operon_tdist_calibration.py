"""Calibration of the M6-01 Student-t engine (GIAP_SPEC.md §7.4.2, §7.4.3, §8.4; suite
``discordance``).

Engine-level claims only; the operon-level discordance level (the t fitted on the OP pairs of a
simulated genome) is M6-09 on M6-03's generators. R is ``replicates_fast/full.discordance``
(100 / 1,000). Replicate ``i`` of stream ``discordance.engine_tdist.<name>`` draws from
``_criteria.replicate_rng``; levels are gated by ``_criteria.calibrate`` through
``_gate.calibrate_levels`` (ADR-0015 two-block confirmation, sharding, jsonl log).

- ``test_two_sided_p_uniform_under_a_known_fit``: with the fit given (no estimation),
  ``2 min(F, 1 - F)`` is exactly U(0, 1) -- the definition's own correctness check (KS criterion).
- ``test_two_sided_p_pooled_level_under_a_known_fit``: the same claim as a pooled level check at
  every alpha in ``levels`` (m p-values per replicate, §8.4 pooled criterion).
- ``test_two_sided_p_level_on_a_held_out_draw``: the fit is estimated from ``N_FIT`` draws and the
  p is computed for an independent draw -- the honest version of the within-operon test.
- ``test_plugin_p_level_recorded``: the fit and the p use the *same* draws, which is what M6-05
  does. It is a plug-in p and cannot be exactly uniform, so the rate is **recorded**
  (``local/logs/calibration/discordance_plugin.jsonl``) with only a sanity bound; the lead decision
  that P(miss) stays out of this BH family does not change that this family is plug-in, and M6-09
  owns the verdict.
- ``test_fit_recovers_its_parameters``: the MLE is close to the drawing parameters on average
  (example 1 over R replicates instead of one seed).
"""

from __future__ import annotations

import json
import math
import time

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels

from giap.stats.tdist import StudentT, fit_student_t, student_t_two_sided

SUITE = "discordance"
FAST_SECONDS = 60.0
#: usable OP pairs on a medium genome (>= MIN_OP_PAIRS = 20 by a wide margin)
N_FIT = 400
#: the within-operon distribution used throughout: a mild location shift and heavy tails
TRUE = StudentT(loc=0.3, scale=0.5, df=5.0, n=N_FIT, converged=True)
#: p-values pooled per replicate in the pooled level check
M_POOLED = 50

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]


def _draw(rng: np.random.Generator, n: int) -> np.ndarray:
    return TRUE.loc + TRUE.scale * rng.standard_t(TRUE.df, size=n)


@pytest.mark.parametrize("full", MODES)
def test_two_sided_p_uniform_under_a_known_fit(full: bool) -> None:
    r = C.replicates(SUITE, full=full)

    def rep(rng: np.random.Generator) -> float:
        return float(student_t_two_sided(_draw(rng, 1), TRUE)[0])

    rec = C.calibrate(
        SUITE, "uniform", 0.0, rep, r, stream=f"{SUITE}.engine_tdist.uniform_known_fit"
    )
    assert rec["status"] in {"pass", "unconfirmed"}, rec


@pytest.mark.parametrize("full", MODES)
def test_two_sided_p_pooled_level_under_a_known_fit(full: bool) -> None:
    levels = C.load_tolerances().levels
    r = C.replicates(SUITE, full=full)

    def rep(rng: np.random.Generator) -> np.ndarray:
        p = student_t_two_sided(_draw(rng, M_POOLED), TRUE)
        return np.array([float((p <= a).mean()) for a in levels])

    gate = calibrate_levels(
        SUITE,
        "level",
        levels,
        rep,
        r,
        stream=f"{SUITE}.engine_tdist.pooled_known_fit",
        m=M_POOLED,
    )
    gate.assert_ok()
    if not full:
        assert gate.block_seconds.get(1, 0.0) < FAST_SECONDS, gate.block_seconds


@pytest.mark.parametrize("full", MODES)
def test_two_sided_p_level_on_a_held_out_draw(full: bool) -> None:
    """Fit on N_FIT draws, test one independent draw: the estimation error is the only slack."""
    levels = C.load_tolerances().levels
    r = C.replicates(SUITE, full=full)

    def rep(rng: np.random.Generator) -> float:
        fit = fit_student_t(_draw(rng, N_FIT))
        assert fit.converged and fit.n == N_FIT
        return float(student_t_two_sided(_draw(rng, 1), fit)[0])

    gate = calibrate_levels(
        SUITE, "level", levels, rep, r, stream=f"{SUITE}.engine_tdist.level_held_out"
    )
    gate.assert_ok()
    if not full:
        assert gate.block_seconds.get(1, 0.0) < FAST_SECONDS, gate.block_seconds


@pytest.mark.parametrize("full", MODES)
def test_plugin_p_level_recorded(full: bool) -> None:
    """The M6-05 situation: the same draws fit the t and are tested against it."""
    levels = C.load_tolerances().levels
    r = C.replicates(SUITE, full=full)
    t0 = time.perf_counter()

    def rep(rng: np.random.Generator) -> np.ndarray:
        x = _draw(rng, N_FIT)
        p = student_t_two_sided(x, fit_student_t(x))
        return np.array([float((p <= a).mean()) for a in levels])

    shares = np.vstack(C.run_replicates(f"{SUITE}.engine_tdist.plugin", rep, r))
    rates = shares.mean(axis=0)
    path = C.log_dir() / "discordance_plugin.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(
            json.dumps(
                {
                    "suite": SUITE,
                    "stream": f"{SUITE}.engine_tdist.plugin",
                    "claim": "plugin_two_sided_p_rejection_rate",
                    "levels": list(levels),
                    "rejection_rate": [float(v) for v in rates],
                    "n_fit": N_FIT,
                    "R": r,
                    "seconds": time.perf_counter() - t0,
                }
            )
            + "\n"
        )
    assert np.all(np.isfinite(rates))
    for alpha, rate in zip(levels, rates, strict=True):
        assert rate <= 2.0 * alpha + 0.01, f"plug-in rate {rate:.4f} at alpha {alpha}"


@pytest.mark.parametrize("full", MODES)
def test_fit_recovers_its_parameters(full: bool) -> None:
    """Example 1 over R replicates: the MLE is unbiased enough to be usable at N_FIT pairs."""
    r = C.replicates(SUITE, full=full)

    def rep(rng: np.random.Generator) -> np.ndarray:
        fit = fit_student_t(_draw(rng, N_FIT))
        return np.array([fit.loc, fit.scale, fit.df, float(fit.converged)])

    fits = np.vstack(C.run_replicates(f"{SUITE}.engine_tdist.recovery", rep, r))
    loc, scale, df, converged = (fits[:, j] for j in range(4))
    assert float(converged.mean()) > 0.95, "the Nelder-Mead fit must converge almost always"
    assert float(np.mean(loc)) == pytest.approx(TRUE.loc, abs=0.02)
    assert float(np.mean(scale)) == pytest.approx(TRUE.scale, abs=0.03)
    assert float(np.median(df)) == pytest.approx(TRUE.df, abs=1.5)
    assert not math.isnan(float(np.mean(df)))
