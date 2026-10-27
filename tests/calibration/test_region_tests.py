"""Calibration of the M4-17 region tests (suite ``regions``; GIPPA_SPEC.md §8.4, §7.3.6, §8.1).

Claims (each gated by ``_criteria.calibrate`` through ``_gate.calibrate_levels``, ADR-0015 two-block
confirmation; R = ``replicates_fast.regions`` / ``replicates_full.regions``):

C1 level (spec acceptance "Calibration"; §8.4 "Region tests | AR(1) abundance, detection from
   p_detect only | level"): with no region effects, ``detection = "covariates_only"`` and the M4-08
   AR(1) TU background, the circular-shift p of the region tests is at level. Gate structure (lead
   decision 7, after the statistician's review; ADR-0015 harness, no tolerance change):
   - chromosome: per statistic, one pooled stream over the first region of each of the four types
     (``check_level_pooled``, m = 4) at alpha in {0.01, 0.05, 0.10};
   - plasmid AMR set (K = 59 shifts): per statistic, ``check_level`` at alpha in {0.05, 0.10} only;
   - full mode: HL at rho in {0, 0.5, 0.9}; O - E (which does not depend on rho under
     covariate-only detection) at rho = 0.5 only.
   Discrete p (K = N - 1 shifts) -> level only, no uniformity check. Fast: 1,000-gene chromosome,
   rho = 0.5, p = ``p_detect_true``. Full: 4,500 genes, p = M4-13 cross-fitted ``p_detect`` on
   M4-10 covariates (xfail "needs M4-13 (#83)" until PR #83 is merged).
   The statistics are those of ADR-0018: studentised O - E, HL on residualised abundance.
C2 FDR across region families (§8.1 "BH within type | 0.05"; ADR-0018 decision 3: family = (region
   type, statistic)): on a genome with 4 prophages, 4 BGCs, 4 CRISPR arrays and 2 AMR sets, all
   null, the per-family FDP at q (= 1 if any q_value <= q, as every region is null) satisfies the
   ``check_fdr`` bound at q in {0.01, 0.05, 0.10}.
C3 power (§7.3.6 "regions with detection shifts of 10-50 %", recorded, not gated): BGC detection
   logit shifts that raise the BGC's mean detection probability by +10, 20, 30, 50 % (relative to
   its covariate-only mean, capped at 0.99), p = the covariate-only ``p_detect_true``; the O - E
   rejection rate at alpha = 0.05 over R = 100 replicates is appended to
   ``local/logs/calibration/regions_power.jsonl``.

Layout (M4-07, spec acceptance): default regions + a 60-gene plasmid ``p1`` with its own 6-gene
AMR set (``RegionSpec.contig``), ``origin_crossing=True``.
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
import pandas as pd
import pytest
from _gate import calibrate_levels

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import region_tests_calib as rc
import region_tests_data as rd

FAST = rc.Scenario("null_fast", 1000, {"ar1_rho": 0.5}, "true")
FULL = [rc.Scenario(f"null_full_rho{r}", 4500, {"ar1_rho": r}, "m4_13") for r in (0.0, 0.5, 0.9)]


# --------------------------------------------------------------------------- C1 level


@pytest.mark.calibration
@pytest.mark.parametrize("stat", rd.STATS)
def test_c1_level_fast(stat: str) -> None:
    R = C.replicates("regions", False)
    failed = rc.pooled_level_check(FAST, stat, R)
    failed += rc.level_checks(FAST, [("amr", stat, "p1")], R, rc.PLASMID_LEVELS)
    assert not failed, failed


#: full mode (lead decision 7): HL at every rho, O - E at rho = 0.5 only
FULL_CHECKS = [(sc, "hl_shift") for sc in FULL] + [(FULL[1], "o_minus_e")]


@pytest.mark.calibration_full
@rc.needs_m4_13
@pytest.mark.parametrize(("sc", "stat"), FULL_CHECKS, ids=lambda x: getattr(x, "name", x))
def test_c1_level_full(sc: rc.Scenario, stat: str) -> None:
    R = C.replicates("regions", True)
    failed = rc.pooled_level_check(sc, stat, R)
    failed += rc.level_checks(sc, [("amr", stat, "p1")], R, rc.PLASMID_LEVELS)
    assert not failed, failed


def test_full_checks_follow_lead_decision_7() -> None:
    assert FULL[1].expr_kw["ar1_rho"] == 0.5
    hl_rhos = sorted(sc.expr_kw["ar1_rho"] for sc, s in FULL_CHECKS if s == "hl_shift")
    assert hl_rhos == [0, 0.5, 0.9]
    assert [sc.expr_kw["ar1_rho"] for sc, s in FULL_CHECKS if s == "o_minus_e"] == [0.5]
    assert 0.01 not in rc.PLASMID_LEVELS


# --------------------------------------------------------------------------- C2 FDR


def _many_regions() -> tuple[Any, ...]:
    import gippasim.layout as la

    spec = la.RegionSpec
    return (
        *(spec("prophage", 40) for _ in range(4)),
        *(spec("bgc", 25) for _ in range(4)),
        *(spec("crispr", 6) for _ in range(4)),
        *(spec("amr", 10, scattered=True) for _ in range(2)),
    )


def _fdp(res: pd.DataFrame, typ: str, stat: str, levels: list[float]) -> np.ndarray:
    q = res.loc[(res["region_type"] == typ) & (res["stat"] == stat), "q_value"].to_numpy(float)
    q = np.where(np.isnan(q), 1.0, q)
    return np.array([float((q <= a).any()) for a in levels])


@pytest.mark.parametrize(
    "full",
    [
        pytest.param(False, marks=pytest.mark.calibration, id="fast"),
        pytest.param(True, marks=[pytest.mark.calibration_full, rc.needs_m4_13], id="full"),
    ],
)
def test_c2_fdr_within_type_and_stat(full: bool) -> None:
    R = C.replicates("regions", full)
    n = 4500 if full else 1000
    sc = rc.Scenario(
        f"fdr_{'full' if full else 'fast'}", n, {"ar1_rho": 0.5}, "m4_13" if full else "true",
        {"regions": _many_regions(), "plasmids": ()},
    )  # fmt: skip
    levels = list(C.load_tolerances().levels)
    failed: list[str] = []
    for typ in rd.REGION_TYPES:
        for stat in rd.STATS:
            stream = f"{sc.stream}.{typ}.{stat}"
            fn = rc.shared_fn(sc, stream, R, lambda r, t=typ, s=stat: _fdp(r, t, s, levels))
            gate = calibrate_levels("regions", "fdr", levels, fn, R, stream=stream)
            failed += gate.detail() if gate.failed else []
    assert not failed, failed


def test_fdr_layout_fits_and_numbers_regions() -> None:
    """The C2 genome has the 14 planted regions, numbered ``<type>_<n>`` (M4-07 §6.5)."""
    import gippasim.layout as la

    lay = la.simulate_layout(
        la.LayoutConfig(n_genes_chromosome=1000, regions=_many_regions()),
        np.random.default_rng(0),
    )
    counts = lay.regions["region_type"].value_counts().to_dict()
    assert counts == {"prophage": 4, "bgc": 4, "crispr": 4, "amr": 2}


# --------------------------------------------------------------------------- C3 power


POWER_SHIFTS = (0.10, 0.20, 0.30, 0.50)


def _logit_shift(p0: np.ndarray, rel: float) -> float:
    """Logit shift raising mean(p0) to min((1 + rel) mean(p0), 0.99)."""
    from scipy.optimize import brentq

    lg = np.log(p0) - np.log1p(-p0)
    target = min((1.0 + rel) * float(p0.mean()), 0.99)

    def f(delta: float) -> float:
        return float(np.mean(1.0 / (1.0 + np.exp(-(lg + delta))))) - target

    return float(brentq(f, 0.0, 50.0))


def _power_replicate(rng: np.random.Generator) -> dict[float, float]:
    import gippasim.expression as ex
    import gippasim.layout as la
    from gippa.context import regions as R
    from gippa.rng import for_step

    lay = la.simulate_layout(rd.calibration_layout_config(1000), rng)
    state = rng.bit_generator.state
    base = ex.ExpressionConfig(detection="covariates_only", ar1_rho=0.5)
    expr0 = ex.simulate_expression(lay, base, rng)
    seed = int(rng.integers(0, 2**31 - 1))
    p0 = expr0.genes["p_detect_true"].to_numpy(float)
    bgc = lay.region_genes.loc[lay.region_genes["region_id"] == "bgc_1", "locus_tag"].tolist()
    in_bgc = expr0.genes["locus_tag"].isin(bgc).to_numpy()
    out: dict[float, float] = {}
    for rel in POWER_SHIFTS:
        delta = _logit_shift(p0[in_bgc], rel)
        rng.bit_generator.state = state
        cfg = ex.ExpressionConfig(
            detection="covariates_only", ar1_rho=0.5, region_detect_shift={"bgc": delta}
        )
        expr = ex.simulate_expression(lay, cfg, rng)
        t = rd.tables_from(lay, expr, p0)
        rows = R.test_region(
            t.genes, bgc, "bgc", "chr1", rd.config(bootstrap=rc.CAL_BOOTSTRAP),
            for_step("context_stats", seed),
        )  # fmt: skip
        oe = next(r for r in rows if r.stat == "o_minus_e")
        out[rel] = rd.p_or_one(oe.p_value)
    return out


@pytest.mark.calibration
def test_c3_power_recorded() -> None:
    R = C.replicates("regions", False)
    t0 = time.perf_counter()
    reps = C.run_replicates("regions.region_tests.power_bgc", _power_replicate, R)
    path = C.log_dir() / "regions_power.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    rates = {rel: float(np.mean([r[rel] <= 0.05 for r in reps])) for rel in POWER_SHIFTS}
    with path.open("a") as fh:
        for rel, rate in rates.items():
            fh.write(
                json.dumps(
                    {
                        "suite": "regions",
                        "stream": "regions.region_tests.power_bgc",
                        "claim": "power",
                        "region_type": "bgc",
                        "stat": "o_minus_e",
                        "detection_shift_rel": rel,
                        "alpha": 0.05,
                        "R": R,
                        "rejection_rate": rate,
                        "seconds": time.perf_counter() - t0,
                    }
                )
                + "\n"
            )
    assert all(0.0 <= v <= 1.0 and not math.isnan(v) for v in rates.values())
    assert path.exists()


def test_logit_shift_hits_relative_target() -> None:
    p0 = np.array([0.2, 0.4, 0.5, 0.6, 0.7])
    for rel in POWER_SHIFTS:
        d = _logit_shift(p0, rel)
        lg = np.log(p0 / (1 - p0)) + d
        got = float(np.mean(1 / (1 + np.exp(-lg))))
        assert got == pytest.approx(min((1 + rel) * p0.mean(), 0.99), abs=1e-9)
