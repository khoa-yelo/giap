"""Calibration of the M5-05 autocorrelation rows (suite ``autocorrelation``; §8.4, §8.1, §7.3.6).

Claims, each gated through ``_criteria.calibrate`` (ADR-0015 two-block confirmation),
R = ``replicates_fast.autocorrelation`` / ``replicates_full.autocorrelation``, B = 199
permutations per null, 1,000-gene chromosome fast / 4,500 full:

C1 **label null** (§8.4 "gene-level noise (label null)"): with no structure at all both nulls
   hold, so for each (statistic, null_model) and alpha in the tolerance levels the gates of
   ADR-0020 decision 5 / lead decision 18 apply:
   - lag-1 level, lag-20 level (``check_level``);
   - the BH family's any-rejection rate over the 20 lags (``check_fwer``, FWER-style bound);
   - the pooled 20-lag rejection share is **recorded**, not gated (lags are positively
     correlated; pooling hid a lag-1 failure).
C2 **TU null** (§8.4 "TU-structured noise without inter-TU structure"; §7.3.6 "autocorrelation
   tests at level under the TU-block null"): the same three gates for ``null_model = tu_block``
   and both statistics. ``detection = "abundance"``, so detection inherits the TU structure and
   the join count is tested where it is hardest.
C3 **teeth** (recorded, not gated; spec acceptance): the ``gene_label`` rejection rate at lag 1
   under the TU null -- the null the TU-block test is deaf to by construction -- is appended to
   ``local/logs/calibration/autocorrelation_teeth.jsonl``. The planning target is >= 0.8
   (spec example 3).

The scenarios are the M5-04 ``label_null`` / ``tu_null`` generators, written out in
``autocorr_calib`` until PR #119 merges (see that module's docstring). The ori-ter gradient
attack of ADR-0020 decision 19 (a, b) needs M5-04's ``oriter_gradient`` and is not written here.

Named ``test_autocorrelation_calibration.py``, not the spec's ``test_autocorrelation.py``:
pytest imports test modules by basename, and ``tests/unit/context/test_autocorrelation.py``
already owns that one (the repo keeps test basenames unique inside ``tests/``). The name also
follows the other suites (``test_detectability_calibration.py``, ...).
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import autocorr_calib as acal
import autocorr_data as ad

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
#: (statistic, null_model) pairs; the TU-null scenario gates ``tu_block`` only
COMBOS = [(s, m) for s in ad.STATISTICS for m in ad.NULLS]
FAST_GENES, FULL_GENES = 1000, 4500


def _scenario(name: str, full: bool) -> acal.Scenario:
    n = FULL_GENES if full else FAST_GENES
    return acal.label_null(n) if name == "label_null" else acal.tu_null(n)


def _levels() -> list[float]:
    return list(C.load_tolerances().levels)


def _level_gate(sc: acal.Scenario, statistic: str, null_model: str, lag: int, R: int) -> list[str]:
    levels = _levels()
    stream = f"{sc.stream}.{statistic}.{null_model}.lag{lag}"
    fn = acal.shared_fn(sc, stream, R, lambda r: acal.lag_p(r, statistic, null_model, lag))
    gate = calibrate_levels(acal.SUITE, "level", levels, fn, R, stream=stream)
    return gate.detail() if gate.failed else []


def _family_gate(sc: acal.Scenario, statistic: str, null_model: str, R: int) -> list[str]:
    levels = _levels()
    stream = f"{sc.stream}.{statistic}.{null_model}.bh_family"
    fn = acal.shared_fn(
        sc, stream, R, lambda r: acal.family_rejection(r, statistic, null_model, levels)
    )
    gate = calibrate_levels(acal.SUITE, "fwer", levels, fn, R, stream=stream)
    return gate.detail() if gate.failed else []


# --------------------------------------------------------------------------- C1 / C2 level


@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize(("statistic", "null_model"), COMBOS)
def test_label_null_lag_level(statistic: str, null_model: str, full: bool) -> None:
    """Lag-1 and lag-20 level under the complete null (both nulls, both statistics)."""
    R = C.replicates(acal.SUITE, full)
    sc = _scenario("label_null", full)
    failed = _level_gate(sc, statistic, null_model, 1, R)
    failed += _level_gate(sc, statistic, null_model, 20, R)
    assert not failed, failed


@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize(("statistic", "null_model"), COMBOS)
def test_label_null_bh_family_any_rejection(statistic: str, null_model: str, full: bool) -> None:
    """The BH family over the 20 lags rejects at most at the nominal rate (decision 18)."""
    R = C.replicates(acal.SUITE, full)
    failed = _family_gate(_scenario("label_null", full), statistic, null_model, R)
    assert not failed, failed


@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize("statistic", ad.STATISTICS)
def test_tu_null_lag_level(statistic: str, full: bool) -> None:
    """Under TU-structured noise the TU-block null holds at lags 1 and 20."""
    R = C.replicates(acal.SUITE, full)
    sc = _scenario("tu_null", full)
    failed = _level_gate(sc, statistic, "tu_block", 1, R)
    failed += _level_gate(sc, statistic, "tu_block", 20, R)
    assert not failed, failed


@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize("statistic", ad.STATISTICS)
def test_tu_null_bh_family_any_rejection(statistic: str, full: bool) -> None:
    R = C.replicates(acal.SUITE, full)
    failed = _family_gate(_scenario("tu_null", full), statistic, "tu_block", R)
    assert not failed, failed


def test_tu_null_scenario_gates_only_the_tu_block_null() -> None:
    """The gene-label null is false under the TU null, so it is never gated there."""
    tu, label = acal.tu_null(FAST_GENES), acal.label_null(FAST_GENES)
    assert tu.nulls == ("tu_block",)
    assert label.nulls == ("gene_label", "tu_block")
    assert tu.expr_cfg.ar1_rho == 0.0 and tu.expr_cfg.tu_sd == 2.0
    assert tu.expr_cfg.detection == "abundance"
    assert label.expr_cfg.tu_sd == 0.0 and label.expr_cfg.position_decay_log2 == 0.0
    # no gradient outside the attack arms, so no oriC either (M5-04 lead decision 16)
    assert tu.expr_cfg.oriter_log2 == 0.0 and label.layout_cfg.oric_frac is None


# --------------------------------------------------------------------------- recorded


def _record(name: str, payload: dict[str, Any]) -> Path:
    path = C.log_dir() / "autocorrelation_teeth.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps({"suite": acal.SUITE, "claim": name, **payload}) + "\n")
    return path


@pytest.mark.calibration
@pytest.mark.parametrize(("statistic", "null_model"), COMBOS)
def test_label_null_pooled_rate_recorded(statistic: str, null_model: str) -> None:
    """Decision 18: the pooled 20-lag rejection share is recorded, never gated."""
    R = C.replicates(acal.SUITE, False)
    sc = _scenario("label_null", False)
    levels = _levels()
    t0 = time.perf_counter()
    stream = f"{sc.stream}.{statistic}.{null_model}.pooled"
    fn = acal.shared_fn(
        sc, stream, R, lambda r: acal.pooled_share(r, statistic, null_model, levels)
    )
    shares = np.vstack(C.run_replicates(stream, fn, R)).mean(axis=0)
    path = _record(
        "pooled_20_lag",
        {
            "stream": stream,
            "statistic": statistic,
            "null_model": null_model,
            "R": R,
            "levels": levels,
            "rates": [float(v) for v in shares],
            "seconds": time.perf_counter() - t0,
        },
    )
    assert path.exists()
    assert all(0.0 <= float(v) <= 1.0 and not math.isnan(float(v)) for v in shares)


@pytest.mark.calibration
@pytest.mark.parametrize("statistic", ad.STATISTICS)
def test_teeth_gene_label_rejects_under_the_tu_null(statistic: str) -> None:
    """Spec example 3 / acceptance "teeth": under TU-structured noise the gene-label test must
    see the structure the TU-block test is built to ignore (planning target >= 0.8 at lag 1).

    Lead decision 28: the target is met by ``morans_i`` (rejection 1.00), not by
    ``join_count_bb`` (0.04), and that is expected, not a bug. giapsim draws detection as
    ``d ~ Bernoulli(p_detect_true)`` and the calibration feeds that same truth in as
    ``p_detect``, so the Pearson residuals (d - p)/sqrt(p(1 - p)) are exact Bernoulli
    residuals: independent across genes whatever structure p itself has, hence no power under
    any abundance-structured scenario. A join-count teeth gate would need a miscalibrated p
    (``p_source="constant"``, or the short-gene prophage of
    ``attacks/test_attack_autocorr.py``). Both rates are recorded and neither is gated; the
    statistician decides at review whether such an arm must gate the join count.
    """
    R = C.replicates(acal.SUITE, False)
    sc = _scenario("tu_null", False)
    t0 = time.perf_counter()
    stream = f"{sc.stream}.{statistic}.gene_label.teeth_lag1"
    fn = acal.shared_fn(sc, stream, R, lambda r: acal.lag_p(r, statistic, "gene_label", 1))
    p = np.array([float(v) for v in C.run_replicates(stream, fn, R)])
    rate = float(np.mean(p <= 0.05))
    path = _record(
        "gene_label_power_tu_null",
        {
            "stream": stream,
            "statistic": statistic,
            "null_model": "gene_label",
            "lag": 1,
            "alpha": 0.05,
            "R": R,
            "rejection_rate": rate,
            "target": 0.8,
            "seconds": time.perf_counter() - t0,
        },
    )
    assert path.exists()
    assert 0.0 <= rate <= 1.0


@pytest.mark.calibration
def test_replicate_frame_has_the_full_row_grid() -> None:
    """One replicate of each scenario: the table the gates read has the 80 rows of example 2."""
    for name in ("label_null", "tu_null"):
        res = acal.replicate_rows(_scenario(name, False), 0)
        assert len(res) == 2 * 2 * len(acal.LAGS)
        assert set(res["contig"]) == {acal.CONTIG}
        assert int(res["n_perm"].max()) == acal.PERMUTATIONS
