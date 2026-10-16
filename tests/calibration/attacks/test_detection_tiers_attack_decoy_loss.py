"""Attack null for ADR-0009: picked gene FDR under random decoy loss (M4-12; GIAP_SPEC.md §7.3.2,
§8.5; M4-index lead decision 8).

Generator: the ``gene_fdr`` replicate of tests/helpers/detection_tiers_calsim.py (M4-07 -> M4-08 ->
M4-09, pi_wrong = 0.05, decoys emitted) with ``floor(L n)`` of the n decoys that pass
q <= ``q_value_max`` re-pointed to unpairable ``rev_LOST_<k>`` IDs, L in {0.10, 0.20, 0.30}. Before
ADR-0009 a 20 % loss inflated picked FDP to about 1.2 q (docs/methods/target-decoy.md).

Claims (criterion ``fdr``, ADR-0015 two-block confirmation, suite ``gene_fdr``; stream
``gene_fdr.attack_decoy_loss.L<pct>``), each at q in {0.01, 0.05, 0.10} with detection rerun at
``stats.alpha.identification`` = q: mean FDP <= q + k SE.

- L = 10 % and 20 %: the paired share is >= 0.80 by construction, so the tier must be
  ``decoy_picked`` (asserted per replicate) and the 1/f rescaling is what keeps the level.
- L = 30 %: the paired share may fall below 0.80 (MAYU fallback, the upper guard) or not (small
  n); the claim holds for whichever tier ``detect`` chooses.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detection_tiers_calsim as cs

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
SUITE = "gene_fdr"
PI = 0.05
LOSSES = [10, 20, 30]
_MEMO: dict[tuple[int, bool, int], dict[float, float]] = {}


def _replicate(pct: int, full: bool, rep: int) -> dict[float, float]:
    key = (pct, full, rep)
    if key not in _MEMO:
        n = cs.N_GENES_FULL if full else cs.N_GENES_FAST
        stream = f"gene_fdr.attack_decoy_loss.sim.L{pct}" + (".full" if full else "")
        rng = C.replicate_rng(stream, rep)
        lay, e = cs.simulate(rng, PI, n)
        p, h = cs.masked(e, "picked")
        p, k = cs.lose_decoys(p, pct / 100, 0.01, rng)
        fdps, tiers = cs.run_variant(lay, e, p, h)
        if pct < 30:
            assert tiers == {"decoy_picked"}, (pct, k, tiers)
        _MEMO[key] = fdps
    return _MEMO[key]


@pytest.mark.parametrize("pct", LOSSES)
@pytest.mark.parametrize("full", MODES)
def test_picked_fdr_under_random_decoy_loss(full: bool, pct: int) -> None:
    R = C.replicates(SUITE, full)
    stream = f"gene_fdr.attack_decoy_loss.L{pct}"
    statuses = {}
    for q in cs.LEVELS:

        def fn(rng: np.random.Generator, q: float = q) -> Any:
            return _replicate(pct, full, cs.rep_index(stream, rng, R))[q]

        statuses[q] = C.calibrate(SUITE, "fdr", q, fn, R, stream=stream)["status"]
    assert "fail" not in statuses.values(), statuses
