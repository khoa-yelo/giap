"""Attack arms 4d-4f on the recalibrated TU-counted ORA, plus the skip path (statistician review
of PR #179; spec M8-04 decision 12, ADR-0035 decision 3; GIAP_SPEC.md §8.4, §8.5).

Proposed repo path: tests/calibration/attacks/test_attack_ora_recal_misfit.py

ADR-0035 claims that the cross-fitted recalibration d_TU ~ logistic(1, logit q_raw, log n)
restores level for TU-packed detection. Arms 4b/4c test that claim with 15-TU random whole-TU
categories, and pass. These arms keep the 4b world (``helpers.ora_calsim.simulate``) and change
only the categories or the detection depth, to the kinds real collections have:

- **4d, selected on detectability.** 15-TU categories drawn from the top third (by q_raw) of the
  >= 4-gene TUs (``large_hiq``: ribosome- or ATP-synthase-like operons). c = 1.5.
- **4e, large categories.** All >= 4-gene TUs in two halves (``big_large``, ~115 TUs each), and
  all >= 2-gene TUs in thirds (``big_multi``, ~200 TUs each), well inside ``ora_size`` [5, 500].
  c = 1.5.
- **4f, deep coverage.** 4b at c = 2.5 with eta0 = +1.5 (~75 % of genes detected), using
  4b's own ``large`` categories.

Review evidence (``$GIAP_ROOT/stat/m804b/``, Slurm 1867861/1867862/1867888, R = 3,000):
recalibrated rejection rates against the pooled bound (alpha 0.01 / 0.05 / 0.10). The oracle q
(exact TU marginal) is at level in every cell, so the misfit is in the recalibration model.

| arm | cell | recal | bound |
|---|---|---|---|
| 4d | large_hiq under | 0.030 / 0.108 / 0.193 | 0.0117 / 0.054 / 0.105 |
| 4e | big_large under | 0.023 / 0.114 / 0.212 | 0.0126 / 0.056 / 0.108 |
| 4e | big_multi over | 0.015 / 0.069 / 0.135 | 0.0121 / 0.055 / 0.106 |
| 4f | large under (fit ran) | 0.018 / 0.070 / 0.126 | 0.0112 / 0.053 / 0.104 |

**Skip path.** When a fold's training TUs have < 10 of either class (14 % of 4f replicates) or
``folds`` is None (M4-13's constant model), q is the raw q, which the 4b teeth show over-rejects
(4f skipped replicates: large under 0.45 at alpha 0.01). The rows must then carry ``FLAG``.

Fallback (ADR-0034 pattern, as test_attack_ora_tusize.py): a scenario whose gate fails passes
only if every ORA row carries ``FLAG``; rates are logged either way.
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
from helpers.ora_calsim import Sim, collection, simulate

pytestmark = [
    pytest.mark.calibration,
    pytest.mark.calibration_nightly,
    pytest.mark.calsuite("ora"),
]

SUITE = "ora"
SIZE = (5, 500)
FLAG = "ora_uncalibrated"
# scenario: (simulate kwargs, {kind: categories per replicate})
SCENARIOS: dict[str, tuple[dict[str, float], dict[str, int]]] = {
    "4d_c1.5": ({"c": 1.5}, {"large_hiq": 5}),
    "4e_c1.5": ({"c": 1.5}, {"big_large": 2, "big_multi": 3}),
    "4f_c2.5_deep": ({"c": 2.5, "eta0": 1.5}, {"large": 10}),
}


def _q_raw(sim: Sim) -> np.ndarray:
    p = sim.genes["p_detect"].to_numpy(float)
    return -np.expm1(np.bincount(sim.tu_of_gene, weights=np.log1p(-p), minlength=sim.sizes.size))


def _sets(sim: Sim, rng: np.random.Generator, kind: str, m: int) -> dict[str, frozenset[str]]:
    sizes = sim.sizes
    large = np.flatnonzero(sizes >= 4)
    if kind == "large":
        pool = rng.permutation(large)
        tus = [pool[i * 15 : (i + 1) * 15] for i in range(m)]
    elif kind == "large_hiq":
        top = large[np.argsort(_q_raw(sim)[large])][-(large.size // 3) :]
        pool = rng.permutation(top)
        tus = [pool[i * 15 : (i + 1) * 15] for i in range(m)]
    else:
        base = large if kind == "big_large" else np.flatnonzero(sizes >= 2)
        pool = rng.permutation(base)
        h = pool.size // m
        tus = [pool[i * h : (i + 1) * h] for i in range(m)]
    out = {f"{kind}{i:03d}": frozenset(g for t in c for g in sim.tu_genes(int(t)))
           for i, c in enumerate(tus) if len(c) >= SIZE[0]}  # fmt: skip
    assert len(out) == m, (kind, len(out))
    return out


def _tables(scenario: str, i: int) -> tuple[dict[str, pd.DataFrame], Sim]:
    from giap.function import ora

    kw, kinds = SCENARIOS[scenario]
    rng = C.replicate_rng(f"{SUITE}.attack_recal_misfit.{scenario}.sim", i)
    sim = simulate(rng, **kw)
    out = {}
    for kind, m in kinds.items():
        sets = _sets(sim, rng, kind, m)
        out[kind] = ora.detection_enrichment(
            sim.genes, [collection(sets)], SIZE, 0.05, folds=sim.folds
        )
    return out, sim


@cache
def _pvalues(scenario: str, i: int) -> dict[tuple[str, str], np.ndarray]:
    out = {}
    for kind, tab in _tables(scenario, i)[0].items():
        for direction in ("over", "under"):
            sub = tab[tab["direction"] == direction].sort_values("term")
            out[(kind, direction)] = np.minimum(sub["p_value"].to_numpy(float), 1.0)
    return out


def _fn(stream: str, R: int, scenario: str, cell: tuple[str, str], m: int, levels: list[float]
        ) -> Callable[[np.random.Generator], np.ndarray]:  # fmt: skip
    index = {C.suite_seed(stream, i): i for i in range(2 * R)}

    def rep(rng: np.random.Generator) -> np.ndarray:
        p = _pvalues(scenario, index[int(rng.bit_generator.seed_seq.entropy)])[cell]
        assert p.size == m and np.all(np.isfinite(p))
        return np.array([np.count_nonzero(p <= a) / m for a in levels])

    return rep


def _all_flagged(tabs: dict[str, pd.DataFrame]) -> bool:
    return all(FLAG in as_list(f) for t in tabs.values() for f in t["flags"])


@pytest.mark.parametrize("scenario", sorted(SCENARIOS))
def test_attack_recal_misfit_at_level_or_flagged(scenario: str) -> None:
    r = C.replicates(SUITE, full=False)
    levels = list(C.load_tolerances().levels)
    failed: list[Any] = []
    detail: list[str] = []
    for kind, m in SCENARIOS[scenario][1].items():
        for direction in ("over", "under"):
            stream = f"{SUITE}.attack_recal_misfit.{scenario}.{kind}.{direction}"
            gate = calibrate_levels(
                SUITE, "level", levels, _fn(stream, r, scenario, (kind, direction), m, levels),
                r, stream=stream, m=m,
            )  # fmt: skip
            failed += gate.failed
            detail += gate.detail()
    if failed:
        assert _all_flagged(_tables(scenario, 0)[0]), (
            f"{scenario}: level fails after recalibration and the rows lack {FLAG}: {detail}"
        )


def test_skipped_recalibration_flags_every_row() -> None:
    """folds=None (M4-13 constant model) uses the raw q: every row carries FLAG."""
    from giap.function import ora

    rng = np.random.default_rng(7)
    sim = simulate(rng, c=2.5)
    tab = ora.detection_enrichment(sim.genes, [collection(_sets(sim, rng, "large", 10))], SIZE,
                                   0.05, folds=None)  # fmt: skip
    assert "tu_recal_skipped" in tab.attrs["run_flags"]
    assert all(FLAG in as_list(f) for f in tab["flags"])


def test_deep_coverage_skip_flags_every_row() -> None:
    """At ~75 % detection some replicates have < 10 undetected training TUs in a fold; those runs
    use the raw q (large-TU under 0.45 at alpha 0.01 in the review) and must carry FLAG."""
    for i in range(200):
        tabs, _ = _tables("4f_c2.5_deep", i)
        tab = tabs["large"]
        if "tu_recal_skipped" in tab.attrs["run_flags"]:
            assert all(FLAG in as_list(f) for f in tab["flags"]), i
            return
    pytest.skip("no skipped replicate in 200 draws")
