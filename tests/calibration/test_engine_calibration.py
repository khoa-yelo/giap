"""Engine calibration through the harness (issue M4-05; GIAP_SPEC.md §8.4).

Each test runs fast (``calibration``, ``replicates_fast``) and full (``calibration_full``,
``replicates_full``). Engine tests reuse the R of the closest consumer suite (M4-index lead
decision 6): ``gene_fdr`` for BH, ``regions`` for permutation, circular shifts and Besag-Clifford.
Replicate i of stream ``<suite>.<test>`` uses ``_criteria.replicate_rng`` (seed
``suite_seed(stream, i)``). Every level is gated by ``_criteria.calibrate`` (ADR-0015 two-block
confirmation, sharding) through ``_gate.calibrate_levels``, which logs to
``local/logs/calibration/<suite>.jsonl``.

- ``test_bh_fdr``: 1,000 p per replicate, 900 U(0, 1) nulls + 100 Beta(0.1, 1) alternatives;
  mean FDP of BH at q <= q + k SE for q in ``levels``.
- ``test_perm_p_level``: two groups of 30 iid N(0, 1), stat = difference of means, B = 199 label
  permutations; level check (discrete p, no KS).
- ``test_circular_shift_level_iid``: 1,000 iid N(0, 1) on a circular index, region = 20 consecutive
  indices, stat = sum, exhaustive shifts; level check.
- ``test_besag_clifford_level``: as ``test_perm_p_level`` with |difference|, h = 10, max_B = 1,000.
"""

from __future__ import annotations

import time
from collections.abc import Callable

import _criteria as C
import numpy as np
import pytest
from _gate import Gate, calibrate_levels

from giap.stats.multitest import bh
from giap.stats.permutation import (
    apply_shift,
    besag_clifford,
    choose_shifts,
    circular_shifts,
    perm_p,
)

FAST_SECONDS = 60.0
N_GROUP = 30

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]


def _diff(x: np.ndarray, labels: np.ndarray) -> float:
    return float(x[labels].mean() - x[~labels].mean())


def _level(
    suite: str, stream: str, full: bool, replicate: Callable[[np.random.Generator], float]
) -> None:
    r = C.replicates(suite, full=full)

    def checked(rng: np.random.Generator) -> float:
        p = replicate(rng)
        assert 0 < p <= 1, p
        return p

    gate = calibrate_levels(suite, "level", C.load_tolerances().levels, checked, r, stream=stream)
    _finish(gate, r, full)


def _finish(gate: Gate, r: int, full: bool) -> None:
    assert all(rec["R"] == r for rec in gate.records)
    gate.assert_ok()
    if not full:  # one simulation per replicate, timed in the first level's block-1 record
        duration = sum(rec["duration"] for rec in gate.records if rec["block"] == 1)
        assert duration < FAST_SECONDS, duration


@pytest.mark.parametrize("full", MODES)
def test_bh_fdr(full: bool) -> None:
    suite, stream = "gene_fdr", "gene_fdr.engine_bh"
    levels = C.load_tolerances().levels
    r = C.replicates(suite, full=full)

    def rep(rng: np.random.Generator) -> np.ndarray:
        p = np.concatenate([rng.uniform(size=900), rng.beta(0.1, 1.0, size=100)])
        null = np.arange(1000) < 900
        q = bh(p)
        out = np.zeros(len(levels))
        for j, a in enumerate(levels):
            called = q <= a
            n = int(called.sum())
            out[j] = (called & null).sum() / n if n else 0.0
        return out

    t0 = time.perf_counter()
    fdp = np.vstack(C.run_replicates(stream, rep, r))
    duration = time.perf_counter() - t0
    assert fdp.shape == (r, len(levels))
    failures = []
    for j, a in enumerate(levels):
        res = C.check_fdr(fdp[:, j], a)
        assert r == res.R
        C.log_result(suite, f"fdr:{stream}", a, res, C.suite_seed(stream, 0), duration)
        if not res.ok:
            failures.append(res.detail)
    assert not failures, failures
    if not full:
        assert duration < FAST_SECONDS, duration


@pytest.mark.parametrize("full", MODES)
def test_perm_p_level(full: bool) -> None:
    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=2 * N_GROUP)
        labels = np.arange(2 * N_GROUP) < N_GROUP
        null = np.array([_diff(x, rng.permutation(labels)) for _ in range(199)])
        return perm_p(_diff(x, labels), null)

    _level("regions", "regions.engine_perm_p", full, rep)


@pytest.mark.parametrize("full", MODES)
def test_circular_shift_level_iid(full: bool) -> None:
    n, width = 1000, 20

    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=n)
        idx = (int(rng.integers(0, n)) + np.arange(width)) % n
        shifts, exhaustive = choose_shifts(
            circular_shifts(n, idx, circular=True), 20_000, n, 1_000, rng
        )
        assert bool(exhaustive) and len(shifts) == n - 1
        null = np.array([x[apply_shift(idx, int(k), n)].sum() for k in shifts])
        return perm_p(float(x[idx].sum()), null)

    _level("regions", "regions.engine_circular_shift", full, rep)


@pytest.mark.parametrize("full", MODES)
def test_besag_clifford_level(full: bool) -> None:
    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=2 * N_GROUP)
        labels = np.arange(2 * N_GROUP) < N_GROUP

        def stat(g: np.random.Generator) -> float:
            return abs(_diff(x, g.permutation(labels)))

        p, b_used = besag_clifford(stat, abs(_diff(x, labels)), 10, 1_000, rng)
        assert 10 <= b_used <= 1_000
        return p

    _level("regions", "regions.engine_besag_clifford", full, rep)
