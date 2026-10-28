"""Attack null 2 for the M4-13 detectability model (GIPPA_SPEC.md §8.5 "short genes concentrated in
prophages"; docs/specs/M4-13-detectability.md, lead decision 5 after the statistician's review;
statistician's draft).

Cross-fitted ``p_detect`` must stay calibrated INSIDE the confounded subgroups, not only on
average: D1's marginal slope cannot see this (dropping the size covariates keeps the slope at
0.93-0.97 while the prophage O - E z is about -3.4; review sims).

Generator (one replicate): ``detectability_data.sim_geneset(rng, 2000, short_genes_in="prophage")``
(M4-07 layout with the short genes placed in prophages + M4-08 ``covariates_only`` detection);
``fit_detectability`` with the default ``DetectabilityConfig`` and
``gippa.rng.for_step("detectability", seed)``, seed drawn from the replicate generator. Subgroups of
the eligible genes: ``prophage`` (``extra["in_prophage"]``) and ``short_decile`` (length_aa at or
below the 10th percentile). Per subgroup z = (Σ d - Σ p) / sqrt(Σ p(1 - p)); replicate value =
two-sided normal p of z.

Checks (suite ``detectability``; R = ``replicates_{fast,full}.detectability``):
- level gate: ``calibrate(level)`` at every alpha in ``levels``, per subgroup. Level only: z is
  conservative under the null (SD ≈ 0.9 prophage, ≈ 0.6 short decile; review sims), so no
  uniformity check.
- teeth: a model without the size covariates (``gravy``, ``pi``, ``tm_segments`` only) is a
  confirmed ``fail`` of the prophage level check at alpha = 0.05 through ``calibrate`` (both
  blocks, PR #69 decision 3; review sims: |z| > 1.96 in ≈ 95 % of replicates).
"""

from __future__ import annotations

import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detectability_data as dd

from gippa.config.model import DetectabilityConfig
from gippa.rng import for_step

SUITE = "detectability"
MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
GROUPS = ("prophage", "short_decile")
NO_SIZE = ["gravy", "pi", "tm_segments"]
N_GENES = 2000
_MEMO: dict[tuple[str, int], dict[str, float]] = {}


def _subgroup_z(rng: np.random.Generator, covs: list[str] | None) -> dict[str, float]:
    from gippa.genetable.detectability import fit_detectability

    gs = dd.sim_geneset(rng, N_GENES, short_genes_in="prophage")
    cfg = DetectabilityConfig() if covs is None else DetectabilityConfig(covariates=covs)
    seed = int(rng.integers(0, 2**31 - 1))
    res = fit_detectability(
        gs.detected, gs.cov, gs.tu_of_gene, cfg, for_step("detectability", seed)
    )
    el = np.asarray(gs.extra["eligible"], bool)
    d_all = gs.detected[el].astype(bool).to_numpy()
    p_all = res.p_detect[el].to_numpy(float)
    assert np.isfinite(p_all).all()
    ln = gs.cov.loc[el, "length_aa"].to_numpy(float)
    masks = {
        "prophage": np.asarray(gs.extra["in_prophage"], bool)[el],
        "short_decile": ln <= np.quantile(ln, 0.1),
    }
    out: dict[str, float] = {}
    for group, m in masks.items():
        assert m.sum() >= 20, (group, int(m.sum()))
        d, p = d_all[m], p_all[m]
        out[group] = float((d.sum() - p.sum()) / np.sqrt(np.sum(p * (1 - p))))
    return out


def _z_fn(stream: str, covs: list[str] | None) -> Callable[[np.random.Generator], dict[str, float]]:
    def wrapped(rng: np.random.Generator) -> dict[str, float]:
        key = (stream, int(rng.bit_generator.seed_seq.entropy))  # type: ignore[attr-defined]
        if key not in _MEMO:
            _MEMO[key] = _subgroup_z(rng, covs)
        return _MEMO[key]

    return wrapped


def _p_two_sided(z: float) -> float:
    return float(2 * stats.norm.sf(abs(z)))


@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize("group", GROUPS)
def test_attack2_subgroup_oe_at_level(full: bool, group: str) -> None:
    stream = "detectability.attack2.default"
    R = C.replicates(SUITE, full)
    zf = _z_fn(stream, None)

    def fn(rng: np.random.Generator) -> float:
        return _p_two_sided(zf(rng)[group])

    statuses: list[Any] = [
        C.calibrate(SUITE, "level", a, fn, R, stream=f"{stream}.{group}")["status"]
        for a in C.load_tolerances().levels
    ]
    assert "fail" not in statuses, (group, statuses)


@pytest.mark.calibration
def test_attack2_no_size_covariates_is_confirmed_fail(teeth_log: Path) -> None:
    """Teeth: without length / n_observable_peptides the prophage O - E level check fails in both
    blocks (the subgroup gate has power against the §8.5 confounding)."""
    stream = "detectability.attack2.no_size"
    zf = _z_fn(stream, NO_SIZE)

    def fn(rng: np.random.Generator) -> float:
        return _p_two_sided(zf(rng)["prophage"])

    rec = C.calibrate(
        SUITE, "level", 0.05, fn, C.replicates(SUITE, False), stream=stream, log_dir=teeth_log
    )
    print(f"no-size detectability model, prophage O - E: {rec}")
    assert rec["status"] == "fail", f"no-size model not a confirmed fail: {rec}"
