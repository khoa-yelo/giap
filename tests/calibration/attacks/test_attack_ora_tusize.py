"""Attack arms 4b and 4c on the TU-counted ORA (issue M8-04b; docs/specs/M8-04 lead decision 12,
ADR-0035 decision 3; GIAP_SPEC.md §8.4, §8.5).

Reference: the statistician's ``$GIAP_ROOT/stat/m8-04a-review/attack_tusize.py`` (Slurm
1867453/1867454). Before recalibration (raw q = 1 - prod(1 - p̂)), rejection rates of 3,000 p
values per cell against the bounds 0.0136 / 0.058 (alpha 0.01 / 0.05):

| scenario | categories | tail | 0.01 | 0.05 |
|---|---|---|---|---|
| c 1.5 | large TUs | under | 0.107 | 0.248 |
| c 1.5 | single-gene TUs | over | 0.019 | 0.079 |
| c 2.5 | large TUs | under | 0.278 | 0.490 |
| c 2.5 | single-gene TUs | over | 0.044 | 0.143 |
| p̂ intercept +0.5 | large TUs | under | 0.0137 | 0.067 |

The statistician's in-sample recalibration probe was at level in every cell (``*_recal.log``).

Arms (runnable here on a gene frame plus ``GeneSets``; M8-14/15 own the gated versions):

- **4b** TU-packed: gene logit = eta + c u_TU, u_TU ~ N(0, 1), c in {1.5, 2.5}; p̂ = the exact
  marginal over u (so p̂ is calibrated per gene, but q = 1 - prod(1 - p̂) is not, per TU);
  categories of 15 whole TUs from TUs of >= 4 genes (``large``) and from single-gene TUs
  (``single``);
- **4c** gene-level miscalibration: c = 0, p̂ = expit(+-0.5 + logit p_true); categories
  ``large``, ``single`` and ``random`` whole TUs.

Each (scenario, category type, direction) is a pooled level gate (m = 10 categories per
replicate, R = ``replicates_fast.ora``) on ``detection_enrichment(..., folds=)`` p-values, the
decision-12 recalibration on M4-13-style folds.

**Fallback (ADR-0034 pattern, decision 12).** If any gate of a scenario is a confirmed fail
after recalibration, the test passes only if every ORA row of that scenario carries the flag
``ora_uncalibrated``; the rates are logged either way (``local/logs/calibration/ora.jsonl``).

**Teeth.** The same 4b data at c = 2.5 with the recalibration off (``folds=None``, the raw q)
over-rejects ``under`` on large-TU categories beyond the level bound at alpha = 0.01.
"""

from __future__ import annotations

from collections.abc import Callable
from functools import cache
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels
from helpers.gsea_data import as_list
from helpers.ora_calsim import collection, simulate, whole_tu_sets

pytestmark = [
    pytest.mark.calibration,
    pytest.mark.calibration_nightly,
    pytest.mark.calsuite("ora"),
]

SUITE = "ora"
M = 10
SIZE = (5, 500)
FLAG = "ora_uncalibrated"
SCENARIOS: dict[str, tuple[dict[str, float], tuple[str, ...]]] = {
    "4b_c1.5": ({"c": 1.5}, ("large", "single")),
    "4b_c2.5": ({"c": 2.5}, ("large", "single")),
    "4c_a+0.5": ({"a": 0.5}, ("large", "single", "random")),
    "4c_a-0.5": ({"a": -0.5}, ("large", "single", "random")),
}


def _tables(scenario: str, i: int, recal: bool = True) -> dict[str, pd.DataFrame]:
    from giap.function import ora

    kw, kinds = SCENARIOS[scenario]
    rng = C.replicate_rng(f"{SUITE}.attack_tusize.{scenario}.sim", i)
    sim = simulate(rng, **kw)
    out = {}
    for kind in kinds:
        sets = whole_tu_sets(sim, rng, kind, M)
        assert len(sets) == M, (scenario, kind, len(sets))
        folds = sim.folds if recal else None
        out[kind] = ora.detection_enrichment(sim.genes, [collection(sets)], SIZE, 0.05, folds=folds)
    return out


@cache
def _pvalues(scenario: str, i: int) -> dict[tuple[str, str], np.ndarray]:
    out = {}
    for kind, tab in _tables(scenario, i).items():
        assert "tu_recal_skipped" not in tab.attrs["run_flags"]
        for direction in ("over", "under"):
            sub = tab[tab["direction"] == direction].sort_values("term")
            assert len(sub) == M
            out[(kind, direction)] = sub["p_value"].to_numpy(float)
    return out


def _fn(stream: str, R: int, scenario: str, cell: tuple[str, str], levels: list[float]
        ) -> Callable[[np.random.Generator], np.ndarray]:  # fmt: skip
    index = {C.suite_seed(stream, i): i for i in range(2 * R)}

    def rep(rng: np.random.Generator) -> np.ndarray:
        p = _pvalues(scenario, index[int(rng.bit_generator.seed_seq.entropy)])[cell]
        assert np.all(np.isfinite(p)) and np.all((p >= 0) & (p <= 1 + 1e-12))
        return np.array([np.count_nonzero(p <= a) / M for a in levels])

    return rep


def _all_flagged(tabs: dict[str, pd.DataFrame]) -> bool:
    return all(FLAG in as_list(f) for t in tabs.values() for f in t["flags"])


@pytest.mark.parametrize("scenario", sorted(SCENARIOS))
def test_attack_tusize_at_level_or_flagged(scenario: str) -> None:
    r = C.replicates(SUITE, full=False)
    levels = list(C.load_tolerances().levels)
    failed: list[Any] = []
    detail: list[str] = []
    for kind in SCENARIOS[scenario][1]:
        for direction in ("over", "under"):
            stream = f"{SUITE}.attack_tusize.{scenario}.{kind}.{direction}"
            gate = calibrate_levels(
                SUITE,
                "level",
                levels,
                _fn(stream, r, scenario, (kind, direction), levels),
                r,
                stream=stream,
                m=M,
            )
            failed += gate.failed
            detail += gate.detail()
    if failed:
        assert _all_flagged(_tables(scenario, 0)), (
            f"{scenario}: level fails after recalibration and the rows lack {FLAG}: {detail}"
        )


def test_teeth_raw_q_over_rejects_large_tus_under() -> None:
    r = 20
    alpha = 0.01
    rejected = 0
    for i in range(r):
        tab = _tables("4b_c2.5", i, recal=False)["large"]
        assert "tu_recal_skipped" in tab.attrs["run_flags"]
        p = tab.loc[tab["direction"] == "under", "p_value"].to_numpy(float)
        rejected += int(np.count_nonzero(p <= alpha))
    rate = rejected / (r * M)
    bound = C.level_bound(alpha, r * M)
    assert rate > bound, f"raw q rejected {rate:.3f} <= bound {bound:.4f}"
