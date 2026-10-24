"""M8-16c: the determinants sensitivity ``agreement`` separates signal from fold noise (spec §8.4:
every statistical claim gets a calibration test; GIAP_SPEC.md §8.6 ``agreement`` for the dominant
determinant block; PR #299 statistician review, finding 1; Khoa's option (a), 2026-10-24).

Claim. ``agreement`` between two determinants reruns that differ only in their TU-fold stream
measures how stable the dominant-block attribution is. If the features explain abundance, the
attribution is stable and ``agreement`` is high; if they explain nothing, it is fold noise. Under
M8-08 rule 3 (largest **positive** SHAP sum) over every explained gene, ``agreement`` floored at
about 0.7 in both cases (genes predicted below the mean pick a noise block just above 0). Under
option (a), the block with the largest **|SHAP|** sum, the two cases separate.

Probe (the statistician's, reproduced). ``sf.planted_world(n=400, seed=s, signal=x)``: log2
abundance = x * cai + N(0, 1), noise in the other four blocks. ``signal`` 0.35 is the no-signal
arm (mean LightGBM CV R² 0.01-0.08), 2.0 the strong arm (0.55-0.59). The primary rerun uses fold
stream 0 of ``sensitivity/determinants``, the comparison reruns streams 1-5 on identical data.

Gate (test author). Per world, **min over the five strong-arm agreements minus max over the
five null-arm agreements >= 0.05**: every strong-signal rerun pair must beat every no-signal pair.
How 0.05 was chosen, from pre-runs (session allocation, ``sf.det_oracle`` with the written-out
|SHAP| rule against the current implementation):

=====  ================  ================  ===============  ===================
seed   |SHAP| null       |SHAP| strong     margin |SHAP|    margin rule 3 (now)
=====  ================  ================  ===============  ===================
5      0.80-0.86         0.97-0.98         0.11             about 0.00
6      0.78-0.79         0.89-0.92         0.10             < 0
7      0.76-0.80         0.93-0.95         0.13             < 0
8      0.73-0.78         0.98-0.99         0.20             0.04
=====  ================  ================  ===============  ===================

The statistician's figures for seed 5 are 0.80-0.86 against 0.97-0.98. 0.05 sits about halfway
between the smallest |SHAP| margin (0.10) and the largest rule-3 margin (0.04). So the gate fails
for a definition with the noise floor and passes the fix with about 2x headroom. It is a
separation criterion, not a level: the absolute agreement at no signal (about 0.75-0.85) is a
fold-noise reference value for the methods note and is not gated here. The worlds and fold streams
are fixed, so the test is deterministic. Four worlds (seeds 5-8) guard against one lucky draw.
Runtime about 12 s per world (12 reruns at n = 400).
"""

from __future__ import annotations

from typing import Any

import pytest
from helpers import sensitivity_function_data as sf

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("determinants")]

STREAM = "sensitivity/determinants"
N_GENES = 400
NULL_SIGNAL, STRONG_SIGNAL = 0.35, 2.0
FOLD_STREAMS = range(1, 6)
MIN_SEPARATION = 0.05


def _arm(seed: int, signal: float) -> tuple[list[float], float]:
    """Agreements of fold streams 1-5 against stream 0, and the mean LightGBM CV R²."""
    from giap.context.sensitivity import agreement
    from giap.rng import for_step

    genes, feats = sf.planted_world(n=N_GENES, seed=seed, signal=signal)
    cfg = sf.base_config()
    runs: list[Any] = [sf.det_result(genes, feats, cfg, for_step(STREAM, s))
                       for s in (0, *FOLD_STREAMS)]  # fmt: skip
    assert all(r is not None for r in runs)
    vals = [agreement(runs[0], r) for r in runs[1:]]
    assert all(v is not None for v in vals), "agreement NA on a 400-gene world"
    return [float(v) for v in vals], sum(float(r.cv_r2) for r in runs) / len(runs)


@pytest.mark.parametrize("seed", [5, 6, 7, 8])
def test_agreement_separates_strong_signal_from_fold_noise(seed: int) -> None:
    null, r2_null = _arm(seed, NULL_SIGNAL)
    strong, r2_strong = _arm(seed, STRONG_SIGNAL)
    assert r2_null < 0.15, f"fixture: null arm has signal (CV R² {r2_null:.3f})"
    assert r2_strong > 0.45, f"fixture: strong arm too weak (CV R² {r2_strong:.3f})"
    margin = min(strong) - max(null)
    assert margin >= MIN_SEPARATION, (
        f"seed {seed}: agreement does not separate signal from fold noise: strong "
        f"{[round(v, 3) for v in strong]} (CV R² {r2_strong:.2f}) vs null "
        f"{[round(v, 3) for v in null]} (CV R² {r2_null:.2f}); margin {margin:.3f} < "
        f"{MIN_SEPARATION}"
    )
