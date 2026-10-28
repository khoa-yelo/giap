"""Calibration of the M6-01 decoy-FDR engine (GIPPA_SPEC.md §7.4.4, §8.2, §8.4; suite
``missed_operons``).

Engine-level claims only: the operon-level gate (planted missed operons through the pair table, the
AR(1) attack, the detection/distance attack) is M6-09, and the generators there are M6-03's. R is
``replicates_fast/full.missed_operons`` (50 / 300), not the spec's "R = 200": ``tolerances.yaml`` is
protected and every suite draws its own R (M4-05). Replicate ``i`` of stream
``missed_operons.engine_decoy_fdr.<name>`` draws from ``_criteria.replicate_rng``; every level is
gated by ``_criteria.calibrate`` through ``_gate.calibrate_levels`` (ADR-0015 two-block
confirmation, sharding, jsonl log).

Scores are ``s = -|Delta|`` and "false" means a target whose truth is "not a missed operon".

- ``test_all_null_fdp``: targets and decoys from the same distribution, pi0 = 1 (the spec's
  acceptance criterion). Every call is false, so mean FDP <= q + 2 SE says the engine almost never
  calls under the complete null.
- ``test_planted_fdp``: 20 % of the targets are co-abundant (small \\|Delta\\|), so the check is not
  vacuous: calls happen at q <= 0.1 and the realised FDP is bounded.
- ``test_distance_confounded_fdp``: targets and decoys share the within-bin null but not the bin
  mix, the case the reweighting exists for. The reweighted FDR is gated; the single-bin FDR on the
  same replicates is recorded (not gated here: an anti-conservative procedure is a teeth test, and
  attacks belong to M6-09).
- ``test_lfdr_centred_under_a_shared_null``: example 9 over R replicates rather than one seed.
"""

from __future__ import annotations

import json
import math
import time
from collections.abc import Callable, Mapping, Sequence

import _criteria as C
import numpy as np
import numpy.typing as npt
import pytest
from _gate import Gate, calibrate_levels

from gippa.stats.decoy_fdr import local_fdr, max_decoy_fdr, reweighted_fdr

SUITE = "missed_operons"
FAST_SECONDS = 60.0
N_TARGET = 1200
N_DECOY = 1200
N_BINS = 3
#: a planted pair is co-abundant, so its |Delta| is far inside the independent null: anything
#: weaker is power, not calibration, and the power grid is M6-09's
PLANTED = 0.4
PLANTED_SCALE = 0.02
#: a bin mix that differs between targets and decoys (the confounding attack of §8.5 applied here)
BIN_SCALE = (0.3, 0.7, 1.5)

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]

Replicate = Callable[[np.random.Generator], npt.NDArray[np.float64]]


def _fdp(
    q: npt.NDArray[np.float64], false: npt.NDArray[np.bool_], levels: Sequence[float]
) -> npt.NDArray[np.float64]:
    out = np.zeros(len(levels))
    for j, a in enumerate(levels):
        called = ~np.isnan(q) & (q <= a)
        n = int(called.sum())
        out[j] = (false & called).sum() / n if n else 0.0
    return out


def _gate(name: str, full: bool, rep: Replicate) -> Gate:
    levels = C.load_tolerances().levels
    r = C.replicates(SUITE, full=full)
    gate = calibrate_levels(SUITE, "fdr", levels, rep, r, stream=f"{SUITE}.engine_decoy_fdr.{name}")
    assert all(rec["R"] == r for rec in gate.records)
    gate.assert_ok()
    if not full:
        assert gate.block_seconds.get(1, 0.0) < FAST_SECONDS, gate.block_seconds
    return gate


def _flat_scenario(
    rng: np.random.Generator, planted: float
) -> tuple[
    npt.NDArray[np.float64],
    npt.NDArray[np.int64],
    Mapping[str, tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]]],
    npt.NDArray[np.bool_],
]:
    """Targets and decoys with the same bin mix; a ``planted`` share of co-abundant targets."""
    is_planted = rng.random(N_TARGET) < planted
    delta = np.where(
        is_planted, rng.normal(0.0, PLANTED_SCALE, N_TARGET), rng.normal(0.0, 1.0, N_TARGET)
    )
    target_bin = rng.integers(0, N_BINS, N_TARGET).astype(np.int64)
    decoys = {
        name: (
            -np.abs(rng.normal(0.0, 1.0, N_DECOY)),
            rng.integers(0, N_BINS, N_DECOY).astype(np.int64),
        )
        for name in ("convergent", "sandwich")
    }
    return -np.abs(delta), target_bin, decoys, ~is_planted


def _binned(
    rng: np.random.Generator, n: int, weights: Sequence[float]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]]:
    """|Delta| ~ HalfNormal(BIN_SCALE[b]) with the bin drawn from ``weights``."""
    bins = rng.choice(N_BINS, size=n, p=np.asarray(weights) / np.sum(weights))
    scale = np.asarray(BIN_SCALE)[bins]
    return -np.abs(rng.normal(0.0, scale)), bins.astype(np.int64)


@pytest.mark.parametrize("full", MODES)
def test_all_null_fdp(full: bool) -> None:
    """The spec's acceptance criterion: same distribution, pi0 = 1, mean FDP <= q + 2 SE."""
    levels = C.load_tolerances().levels

    def rep(rng: np.random.Generator) -> npt.NDArray[np.float64]:
        s, b, decoys, false = _flat_scenario(rng, planted=0.0)
        return _fdp(max_decoy_fdr(s, b, decoys, n_bins=N_BINS).q, false, levels)

    _gate("all_null", full, rep)


@pytest.mark.parametrize("full", MODES)
def test_planted_fdp(full: bool) -> None:
    """40 % strongly co-abundant targets: calls happen, and the realised FDP is bounded.

    The planted |Delta| has to be far inside the null for the engine to call anything at q <= 0.1:
    with pi0 = 1 and N_NO,b counting every target, FDR(t) ~ N_total F_decoy(t) / T(t), so at a
    realistic effect size this family simply has little power (recorded in M6-09's power grid).
    """
    levels = C.load_tolerances().levels
    called_any: list[float] = []

    def rep(rng: np.random.Generator) -> npt.NDArray[np.float64]:
        s, b, decoys, false = _flat_scenario(rng, planted=PLANTED)
        q = max_decoy_fdr(s, b, decoys, n_bins=N_BINS).q
        called_any.append(float((~np.isnan(q) & (q <= 0.10)).sum()))
        return _fdp(q, false, levels)

    _gate("planted", full, rep)
    assert np.mean(called_any) > 10.0, f"vacuous: {np.mean(called_any)} mean calls at q <= 0.10"


@pytest.mark.parametrize("full", MODES)
def test_distance_confounded_fdp(full: bool) -> None:
    """Targets and decoys share the within-bin null but not the bin mix (all targets null).

    This is the case the reweighting exists for: pooled over bins the decoys look systematically
    worse than the targets, so a single-bin FDR underestimates the false count. The gate is that
    the **reweighted** FDR still controls the FDP; the single-bin FDP on the same replicates is
    *recorded* (``local/logs/calibration/missed_operons_reweighting.jsonl``) and only has to be no
    better than the reweighted one. Whether a realistic bin mix can push the single-bin FDR past
    its level at all is attack A, which is M6-09's with M6-03's generators.
    """
    levels = C.load_tolerances().levels
    naive: list[npt.NDArray[np.float64]] = []

    def rep(rng: np.random.Generator) -> npt.NDArray[np.float64]:
        s, b = _binned(rng, N_TARGET, (0.7, 0.2, 0.1))
        d_s, d_b = _binned(rng, N_DECOY, (0.1, 0.2, 0.7))
        false = np.ones(N_TARGET, dtype=bool)  # every target is null
        flat = reweighted_fdr(s, np.zeros_like(b), d_s, np.zeros_like(d_b), n_bins=1)
        order = np.argsort(-s, kind="stable")
        q_flat = np.empty_like(flat)
        q_flat[order] = np.minimum.accumulate(np.clip(flat[order], 0.0, 1.0)[::-1])[::-1]
        naive.append(_fdp(q_flat, false, levels))
        return _fdp(max_decoy_fdr(s, b, {"convergent": (d_s, d_b)}, n_bins=N_BINS).q, false, levels)

    gate = _gate("distance_confounded", full, rep)
    reweighted = gate.mean()
    assert reweighted is not None
    single_bin = np.vstack(naive).mean(axis=0)
    path = C.log_dir() / "missed_operons_reweighting.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(
            json.dumps(
                {
                    "suite": SUITE,
                    "stream": f"{SUITE}.engine_decoy_fdr.distance_confounded",
                    "claim": "reweighting_vs_single_bin",
                    "levels": list(levels),
                    "mean_fdp_reweighted": [float(v) for v in reweighted],
                    "mean_fdp_single_bin": [float(v) for v in single_bin],
                    "R": C.replicates(SUITE, full=full),
                }
            )
            + "\n"
        )
    assert np.all(single_bin >= reweighted - 1e-12), (
        f"single-bin {single_bin} beats reweighted {reweighted} under confounding"
    )
    check = np.random.default_rng(6101)
    assert np.mean(-_binned(check, 4000, (0.1, 0.2, 0.7))[0]) > 2.0 * np.mean(
        -_binned(check, 4000, (0.7, 0.2, 0.1))[0]
    ), "the generator is not confounded: decoy |Delta| must dominate target |Delta|"


@pytest.mark.parametrize("full", MODES)
def test_lfdr_centred_under_a_shared_null(full: bool) -> None:
    """Example 9 over R replicates: targets and decoys from one distribution give lfdr ~ 1."""
    r = C.replicates(SUITE, full=full)
    t0 = time.perf_counter()

    def rep(rng: np.random.Generator) -> float:
        target = -np.abs(rng.normal(size=500))
        decoy = -np.abs(rng.normal(size=500))
        return float(np.nanmedian(local_fdr(target, decoy, np.ones(500))))

    medians = C.run_replicates(f"{SUITE}.engine_decoy_fdr.lfdr_null", rep, r)
    mean_median = float(np.mean(medians))
    path = C.log_dir() / "missed_operons_lfdr.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(
            json.dumps(
                {
                    "suite": SUITE,
                    "stream": f"{SUITE}.engine_decoy_fdr.lfdr_null",
                    "claim": "lfdr_median_under_shared_null",
                    "mean_median_lfdr": mean_median,
                    "min": float(np.min(medians)),
                    "max": float(np.max(medians)),
                    "R": r,
                    "seconds": time.perf_counter() - t0,
                }
            )
            + "\n"
        )
    assert not math.isnan(mean_median)
    assert 0.8 <= mean_median <= 1.2, medians
