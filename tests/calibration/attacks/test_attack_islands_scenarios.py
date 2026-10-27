"""M5-07 attack 1 (GIPPA_SPEC.md §8.5 "Autocorrelated abundance with no islands -> islands,
autocorrelation") on the gippasim scenarios; spec docs/specs/M5-07-island-calibration-power.md
(Examples 3, lead decisions 28-30, lead note 2026-10-21).

Every arm stays inside the island null (no island planted) and is gated with ``fwer`` through
``_gate.calibrate_levels`` (ADR-0015 two-block confirmation), full mode only (spec), R =
``replicates_full.islands`` = 500 (Example 3: each <= 0.0695 at alpha 0.05), B = PERM = 1000 as in
the gate (test_island_fwer.py).

- (a) ``ar1_background(4500, 0.95)``: strong background; l' stays capped at 20 (asserted in every
  replicate).
- (b) ``oriter_gradient`` 2 log2 cosine at rho 0.5 is ``test_attack_islands.py::
  test_attack_oriter_gradient_scenario`` (1,800 genes); here the same at the 4,500-gene gate layout.
- (c) ``tu_null`` with mean TU size 3 (large TUs, no inter-TU structure; detection from abundance,
  the scenario default).
- (d) ``ar1_background(4500, 0.5)`` with ``detection = "abundance"`` (censoring carries signal).
- (e) lead note 2026-10-21 / M5-06 decision 27: ``tu_null(1800, short_genes_in="prophage")``,
  clustered censoring in the short-gene prophage, gated.

These files' replicates reuse ``test_island_fwer.replicate`` (the step driver, l' <= 20 check).
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path
from typing import Any

import _criteria as C
import pytest
from _gate import calibrate_levels

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from test_island_fwer import replicate

SUITE = "islands"
pytestmark = [pytest.mark.calsuite(SUITE), pytest.mark.calibration_full]


def _abundance_detection(pair: Any) -> Any:
    from gippasim import scenarios as sc

    lcfg, ecfg = pair
    return lcfg, dataclasses.replace(ecfg, detection="abundance", beta=sc.ABUNDANCE_BETA)


def _arms() -> dict[str, Any]:
    from gippasim import scenarios as sc

    tu3_l, tu3_e = sc.tu_null(4500)
    return {
        "a_ar1_rho095": sc.ar1_background(4500, 0.95),
        "b_oriter_gradient_2log2_n4500": sc.oriter_gradient(4500, 0.5, 2.0),
        "c_tu_null_mean_tu3": (dataclasses.replace(tu3_l, mean_tu_size=3.0), tu3_e),
        "d_ar1_rho05_detection_abundance": _abundance_detection(sc.ar1_background(4500, 0.5)),
        "e_tu_null_short_gene_prophage": sc.tu_null(1800, short_genes_in="prophage"),
    }


ARMS = (
    "a_ar1_rho095",
    "b_oriter_gradient_2log2_n4500",
    "c_tu_null_mean_tu3",
    "d_ar1_rho05_detection_abundance",
    "e_tu_null_short_gene_prophage",
)


@pytest.mark.parametrize("arm", ARMS)
def test_attack1_island_fwer_at_level(arm: str) -> None:
    """Example 3: the share of replicates with an island stays within the FWER bound."""
    levels = C.load_tolerances().levels
    pair = _arms()[arm]
    gate = calibrate_levels(
        SUITE, "fwer", levels, lambda rng: replicate(rng, pair, levels),
        C.replicates(SUITE, full=True), stream=f"{SUITE}.attack1.{arm}",
    )  # fmt: skip
    print(f"attack 1 {arm}: {gate.detail()}")
    gate.assert_ok()


def test_arms_are_the_spec_scenarios() -> None:
    """The arms carry the knobs the spec names (no island, the stated rho / detection / TUs)."""
    arms = _arms()
    assert set(arms) == set(ARMS)
    for _lcfg, ecfg in arms.values():
        assert ecfg.islands == ()
    assert arms["a_ar1_rho095"][1].ar1_rho == 0.95
    assert arms["b_oriter_gradient_2log2_n4500"][1].oriter_log2 == 2.0
    assert arms["b_oriter_gradient_2log2_n4500"][1].oriter_shape == "cosine"
    assert arms["c_tu_null_mean_tu3"][0].mean_tu_size == 3.0
    assert arms["c_tu_null_mean_tu3"][1].ar1_rho == 0.0
    assert arms["d_ar1_rho05_detection_abundance"][1].detection == "abundance"
    assert arms["e_tu_null_short_gene_prophage"][0].short_genes_in == "prophage"
    assert arms["e_tu_null_short_gene_prophage"][0].n_genes_chromosome == 1800
