"""Attack nulls for the M4-17 region tests (GIAP_SPEC.md §8.5 attacks 2 and 4; statistician-owned).

A2 detectability confounded with region membership (spec acceptance "Attack nulls" (2); backlog:
   attack-2 suites gate on M4-13 ``p_detect``, PR #83): prophage genes drawn short
   (``short_genes_in="prophage"``, median length x 0.3), detection from covariates only
   (``detection = "covariates_only"``), p = M4-13 cross-fitted ``p_detect`` on the M4-10 covariates
   in fast and full mode. Claim: the prophage O - E p satisfies the level bound at alpha in
   {0.01, 0.05, 0.10}. Teeth: the same data with a constant p (the eligible detected share, i.e.
   detectability ignored) is a confirmed ``fail`` at alpha = 0.05 through ``calibrate``.
A4 high-copy plasmid (spec acceptance "Attack nulls" (4)): ``plasmid_copy_log2 = {"p1": log2(20)}``.
   Claim: the plasmid's own AMR set (O - E and HL against the plasmid's other genes) and every
   chromosome region test (prophage, bgc, crispr, amr x O - E, HL) stay at level. Fast: p =
   ``p_detect_true``; full: M4-13 ``p_detect``.

Layout as the ``regions`` suite (M4-07 default regions + 60-gene plasmid with its own AMR set,
origin crossing; 1,000 genes fast / 4,500 full), M4-08 AR(1) rho = 0.5.
R = ``replicates_*.regions``.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import _criteria as C
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import region_tests_calib as rc

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]


def _a2(full: bool, p_source: str = "m4_13") -> rc.Scenario:
    return rc.Scenario(
        f"attack2_short_prophage_{'full' if full else 'fast'}_{p_source}",
        4500 if full else 1000,
        {"ar1_rho": 0.5},
        p_source,
        {"short_genes_in": "prophage"},
    )


@pytest.mark.parametrize("full", MODES)
def test_a2_short_prophage_genes_o_minus_e_at_level(full: bool) -> None:
    R = C.replicates("regions", full)
    failed = rc.level_checks(_a2(full), [("prophage", "o_minus_e", "chr1")], R)
    assert not failed, failed


@pytest.mark.calibration
def test_a2_teeth_constant_p_detect_rejects(teeth_log: Path) -> None:
    """Ignoring detectability (constant p) on the attack-2 data is a confirmed fail at 0.05."""
    R = C.replicates("regions", False)
    sc = _a2(False, "constant")
    stream = f"{sc.stream}.prophage.o_minus_e"
    fn = rc.shared_fn(sc, stream, R, lambda res: rc.region_p(res, "prophage", "o_minus_e", "chr1"))
    rec = C.calibrate("regions", "level", 0.05, fn, R, stream=stream, log_dir=teeth_log)
    assert rec["status"] == "fail", f"constant p_detect not a confirmed fail: {rec}"


def _a4(full: bool) -> rc.Scenario:
    return rc.Scenario(
        f"attack4_copy20_{'full' if full else 'fast'}",
        4500 if full else 1000,
        {"ar1_rho": 0.5, "plasmid_copy_log2": {"p1": math.log2(20.0)}},
        "m4_13" if full else "true",
    )


@pytest.mark.parametrize("full", MODES)
def test_a4_high_copy_plasmid_at_level(full: bool) -> None:
    R = C.replicates("regions", full)
    failed = rc.level_checks(_a4(full), rc.PLASMID_AMR_CHECKS + rc.CHROMOSOME_CHECKS, R)
    assert not failed, failed
