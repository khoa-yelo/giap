"""Attack nulls for the M4-17 region tests (GIAP_SPEC.md §8.5 attacks 2 and 4; statistician-owned).

Statistics as ADR-0018: studentised O - E (effect O - E), HL on log2 abundance residualised on the
M4-13 detectability covariates. Layout as the ``regions`` suite (M4-07 default regions + 60-gene
plasmid with its own AMR set, origin crossing; 1,000 genes fast / 4,500 full), M4-08 AR(1)
rho = 0.5; R = ``replicates_*.regions``; level at alpha in {0.01, 0.05, 0.10}.

A2 detectability confounded with region membership (prophage genes drawn short,
   ``short_genes_in="prophage"``, median length x 0.3), both directions (lead decision 8):
   - A2-low, O - E (spec acceptance "Attack nulls" (2); conservative direction for the plain sum):
     default detection intercept, so the short prophage genes sit at low p; p = M4-13 cross-fitted
     ``p_detect`` (backlog: attack-2 suites gate on M4-13, PR #83); xfail(strict=False)
     "needs M4-13 (#83)" while M4-13 is absent (lead decisions 4, 8).
   - A2-mid, O - E (anti-conservative direction, statistician limitation (i)): detection intercept
     +2.0, so the chromosome's p_detect averages ~0.87 and the short prophage genes ~0.47
     (p(1 - p) ~0.21 inside vs ~0.09 outside). The plain O - E sum rejects ~0.17 at alpha = 0.05
     here (test-author run, 300 replicates, p = p_detect_true); the studentised statistic ~0.057.
     Fast: p = ``p_detect_true`` (isolates the variance heterogeneity); full: M4-13 ``p_detect``.
   - A2-HL (anti-conservative direction, statistician limitation (ii)): ``detection = "abundance"``
     (beta_log2_abundance 0.5, intercept offset -0.5 x mu_log2 as in M4-08), so only the short
     prophage's abundant genes are detected. Prophage HL must stay at level (raw HL rejects
     0.21-0.25 at alpha = 0.05). HL does not use p_detect: p = ``p_detect_true`` in both modes.
   - Teeth (unchanged): constant p on the A2-low data is a confirmed ``fail`` at alpha = 0.05.
A4 high-copy plasmid (spec acceptance "Attack nulls" (4); lead decision 7): a deterministic
   invariance instead of a calibration stream. Every comparison is within one replicon, so results
   at ``plasmid_copy_log2 = {"p1": log2(20)}`` equal those at copy number 1 on the same seed:
   chromosome rows and all O - E rows bitwise; plasmid HL ``effect`` / CI within 1e-9 (log2 + c
   then a difference is exact only up to rounding) with identical ``p_value`` and ``n_perm``.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import _criteria as C
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import region_tests_calib as rc
import region_tests_data as rd

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
A2_MID_INTERCEPT_SHIFT = 2.0
A2_HL_BETA_ABUNDANCE = 0.5


def _beta(**shift: float) -> dict[str, float]:
    import giapsim.expression as ex

    b = dict(ex.DEFAULT_BETA)
    for k, v in shift.items():
        b[k] = b[k] + v
    return b


def _a2(full: bool, p_source: str = "m4_13") -> rc.Scenario:
    return rc.Scenario(
        f"attack2_short_prophage_{'full' if full else 'fast'}_{p_source}",
        4500 if full else 1000,
        {"ar1_rho": 0.5},
        p_source,
        {"short_genes_in": "prophage"},
    )


def _a2_mid(full: bool) -> rc.Scenario:
    return rc.Scenario(
        f"attack2_mid_prophage_{'full' if full else 'fast'}",
        4500 if full else 1000,
        {"ar1_rho": 0.5, "beta": _beta(intercept=A2_MID_INTERCEPT_SHIFT)},
        "m4_13" if full else "true",
        {"short_genes_in": "prophage"},
    )


def _a2_hl(full: bool) -> rc.Scenario:
    b4 = A2_HL_BETA_ABUNDANCE
    return rc.Scenario(
        f"attack2_hl_abundance_{'full' if full else 'fast'}",
        4500 if full else 1000,
        {
            "ar1_rho": 0.5,
            "detection": "abundance",
            "beta": _beta(log2_abundance=b4, intercept=-b4 * 20.0),
        },
        "true",
        {"short_genes_in": "prophage"},
    )


@pytest.mark.parametrize("full", MODES)
@rc.needs_m4_13
def test_a2_short_prophage_genes_o_minus_e_at_level(full: bool) -> None:
    R = C.replicates("regions", full)
    failed = rc.level_checks(_a2(full), [("prophage", "o_minus_e", "chr1")], R)
    assert not failed, failed


@pytest.mark.parametrize(
    "full",
    [
        pytest.param(False, marks=pytest.mark.calibration, id="fast"),
        pytest.param(True, marks=[pytest.mark.calibration_full, rc.needs_m4_13], id="full"),
    ],
)
def test_a2_mid_range_prophage_p_o_minus_e_at_level(full: bool) -> None:
    """Anti-conservative direction: region p_detect nearer 0.5 than the replicon's."""
    R = C.replicates("regions", full)
    failed = rc.level_checks(_a2_mid(full), [("prophage", "o_minus_e", "chr1")], R)
    assert not failed, failed


@pytest.mark.parametrize("full", MODES)
def test_a2_detection_by_abundance_prophage_hl_at_level(full: bool) -> None:
    """Anti-conservative direction of HL: detection selects the short prophage's abundant genes."""
    R = C.replicates("regions", full)
    failed = rc.level_checks(_a2_hl(full), [("prophage", "hl_shift", "chr1")], R)
    assert not failed, failed


@pytest.mark.parametrize("full", MODES)
def test_a2_short_prophage_hl_at_level_without_abundance_selection(full: bool) -> None:
    """A2-HL-MAR (statistician, review round 2): the A2-low data (short prophage genes, detection
    from covariates only, so abundance is independent of detection and of membership and there is
    no selection bias to residualise). Only ~6 of the prophage's 40 genes are detected against ~22
    at a shifted position, so a plain HL shift statistic has a wider sampling distribution at the
    observed position than across the shifts: plain or residualised HL rejects ~0.13-0.14 at
    alpha = 0.05 (R = 1,000-2,000, 1,000 and 4,500 genes). Isolates the variance mechanism of the
    A2-HL failure from the detection-selection bias; HL does not use p_detect (p = truth)."""
    R = C.replicates("regions", full)
    failed = rc.level_checks(_a2(full, "true"), [("prophage", "hl_shift", "chr1")], R)
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


def test_a2_scenarios_put_the_prophage_where_claimed() -> None:
    """Scenario sanity (one draw each): A2-mid puts the prophage's mean p_detect nearer 0.5 than
    the chromosome's (larger mean p (1 - p)); A2-low puts it below; A2-HL detects by abundance."""
    out = {}
    for name, sc in (("low", _a2(False, "true")), ("mid", _a2_mid(False))):
        t = rd.sim_tables(
            np.random.default_rng(31), rd.calibration_layout_config(1000, **dict(sc.layout_kw)),
            sc.expr_kw, "true",
        )  # fmt: skip
        g = t.genes[(t.genes["contig"] == "chr1") & t.genes["p_detect"].notna()]
        pro = t.region_genes.loc[t.region_genes["region_id"] == "prophage_1", "locus_tag"]
        inside = g["locus_tag"].isin(set(pro)).to_numpy()
        p = g["p_detect"].to_numpy(float)
        out[name] = (p[inside].mean(), p[~inside].mean(), (p * (1 - p))[inside].mean(),
                     (p * (1 - p))[~inside].mean())  # fmt: skip
    lo_in, lo_out, _, _ = out["low"]
    mid_in, mid_out, v_in, v_out = out["mid"]
    assert lo_in < lo_out
    assert abs(mid_in - 0.5) < abs(mid_out - 0.5) and mid_out > 0.8
    assert v_in > 1.5 * v_out
    assert _a2_hl(False).expr_kw["detection"] == "abundance"


# --------------------------------------------------------------------------- A4


def _a4_results(seed: int, copy_log2: float) -> pd.DataFrame:
    from giap.rng import for_step

    rng = np.random.default_rng(seed)
    kw = {"ar1_rho": 0.5, "plasmid_copy_log2": {"p1": copy_log2}}
    t = rd.sim_tables(rng, rd.calibration_layout_config(1000), kw, p_source="true")
    s = int(rng.integers(0, 2**31 - 1))
    res = rd.run_region_tests(
        t, rd.config(bootstrap=rc.CAL_BOOTSTRAP), for_step("context_stats", s)
    )
    out = res.sort_values(["region_id", "contig", "stat"]).reset_index(drop=True)
    out["flags"] = [";".join(sorted(rd.flag_set(f))) for f in out["flags"]]
    return out


@pytest.mark.parametrize("seed", [401, 402, 403])
def test_a4_high_copy_plasmid_invariance(seed: int) -> None:
    a = _a4_results(seed, 0.0)
    b = _a4_results(seed, math.log2(20.0))
    cols = list(rd.TABLE_COLUMNS)
    assert set(a["contig"]) == {"chr1", "p1"}
    loose = (b["contig"] == "p1") & (b["stat"] == "hl_shift")
    pd.testing.assert_frame_equal(a.loc[~loose, cols], b.loc[~loose, cols], check_exact=True)
    pd.testing.assert_frame_equal(
        a.loc[loose, cols].drop(columns=["effect", "ci_low", "ci_high", "observed"]),
        b.loc[loose, cols].drop(columns=["effect", "ci_low", "ci_high", "observed"]),
        check_exact=True,
    )
    for col in ("effect", "ci_low", "ci_high", "observed"):
        np.testing.assert_allclose(
            a.loc[loose, col].to_numpy(float), b.loc[loose, col].to_numpy(float),
            rtol=0.0, atol=1e-9, equal_nan=True,
        )  # fmt: skip
