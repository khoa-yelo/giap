"""Coverage of the M5-08 ori-ter slope CI (suite ``islands``; GIAP_SPEC.md §8.4, §7.3.5(e)).

The ori-ter regression is **descriptive** (no p-value, flag ``descriptive``), so the only
statistical claim is the circular block-bootstrap CI of the slope. ADR-0021 decision 1 fixes its
block length (l' = min(20, 2 l-hat + 1, T)) and Amendment 1 states the validity domain as a
dependence regime -- giapsim ``ar1_background`` at rho <= 0.5 -- not a replicon size; no row
carries ``uncalibrated``.

Claims:

C1 (spec example 7, **recorded**, not gated) ``scenarios.oriter_gradient(4500, 0.5, 2.0)``:
   the pooled-slope CI covers the OLS slope of the planted cosine curve in >= 90 % of R seeds
   (R = ``replicates_fast/full.islands``: 50 / 500). 2 log2 ori minus ter is a 4-fold
   copy-number ratio and, quoted in background SD (giapsim's AR(1) innovation SD is 2, so the
   stationary SD at rho = 0.5 is 2.31), a 0.87 SD peak-to-trough gradient -- the realistic
   amplitude the statistician measured at level. Rate appended to
   ``local/logs/calibration/oriter_coverage.jsonl``.
C2 (**recorded**, ``calibration_full``) the same at 4 log2 (1.73 SD, a 16-fold ratio): outside
   the regime the island null is calibrated for; recorded so the hardening pass can see it.
C3 (**teeth**, ``calibration_full``) resampling single TUs (``context.islands.block_len = 1``)
   instead of blocks must cover worse at rho = 0.5: without it the CI would not depend on the
   block length at all.

Each replicate draws its layout and expression from ``replicate_rng`` and runs the real
``giap.context.oriter.oriter_rows`` with ``giap.rng.for_step("context_stats.oriter", seed)``.
"""

from __future__ import annotations

import json
import math
import sys
import time
import warnings
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import oriter_data as od

SUITE = "islands"
B, PERM = 199, 199
#: spec example 7: 4,500 genes, rho 0.5, 2 log2 (0.87 background SD peak to trough)
FAST_GENES, FULL_GENES = 4500, 4500
RHO = 0.5
AMPLITUDE = 2.0
#: planning target of example 7 (recorded, not gated: the analysis is descriptive)
COVERAGE_TARGET = 0.90
MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]


def amplitude_sd(amplitude: float, rho: float, innovation_sd: float = 2.0) -> float:
    """``amplitude`` log2 quoted in stationary background SD (lead decision 49's convention)."""
    return amplitude * math.sqrt(1.0 - rho * rho) / innovation_sd


def _replicate(
    rng: np.random.Generator,
    n_genes: int,
    amplitude: float,
    rho: float = RHO,
    block_len: int | str = "auto",
) -> tuple[float, float]:
    """(covered?, CI width) of one replicate; NaN when the CI does not exist."""
    from giapsim import scenarios as sc

    genes, tus, replicons, feats = od.scenario_frames(
        rng, sc.oriter_gradient(n_genes, rho, amplitude)
    )
    cfg = od.config(bootstrap=B, permutations=PERM, block_len=block_len,
                    seed=int(rng.integers(1, 2**31 - 1)))  # fmt: skip
    rows = od.rows(genes, tus, replicons, feats, cfg)
    _, lo, hi = od.pooled_ci(rows)
    target = od.planted_slope(od.oracle_tus(genes, tus, replicons, feats, "chr1"), amplitude)
    if not (math.isfinite(lo) and math.isfinite(hi) and math.isfinite(target)):
        return math.nan, math.nan
    return float(lo <= target <= hi), float(hi - lo)


def _coverage(name: str, R: int, **kw: Any) -> dict[str, Any]:
    """Run R replicates, append the record to the suite log and return it."""
    t0 = time.perf_counter()
    out = C.run_replicates(f"{SUITE}.oriter.{name}", lambda r: _replicate(r, **kw), R)
    cov = np.array([c for c, _ in out], dtype=float)
    width = np.array([w for _, w in out], dtype=float)
    ok = np.isfinite(cov)
    record = {
        "suite": SUITE,
        "stream": f"{SUITE}.oriter.{name}",
        "claim": "ci_coverage",
        "statistic": "oriter_slope_pooled",
        "R": R,
        "estimable": float(ok.mean()) if R else math.nan,
        "coverage": float(cov[ok].mean()) if ok.any() else math.nan,
        "mean_ci_width": float(np.nanmean(width)) if ok.any() else math.nan,
        "target": COVERAGE_TARGET,
        "seconds": time.perf_counter() - t0,
        **{k: (v if isinstance(v, (int, float, str)) else str(v)) for k, v in kw.items()},
    }
    path = C.log_dir() / "oriter_coverage.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(record) + "\n")
    print(f"recorded oriter coverage {name}: {record}")
    return record


# --------------------------------------------------------------------------- C1


@pytest.mark.parametrize("full", MODES)
def test_pooled_slope_ci_covers_the_planted_curve(full: bool) -> None:
    """Example 7: coverage of the planted OLS slope at 2 log2 (0.87 background SD), rho 0.5."""
    R = C.replicates(SUITE, full=full)
    n_genes = FULL_GENES if full else FAST_GENES
    assert amplitude_sd(AMPLITUDE, RHO) == pytest.approx(0.866, abs=1e-3)
    rec = _coverage("gradient2log2", R, n_genes=n_genes, amplitude=AMPLITUDE)
    assert rec["estimable"] >= 0.95, rec
    assert 0.0 <= rec["coverage"] <= 1.0
    if rec["coverage"] < COVERAGE_TARGET:
        warnings.warn(
            f"ori-ter pooled CI coverage {rec['coverage']:.3f} < {COVERAGE_TARGET}"
            " (recorded, not gated: the analysis is descriptive)",
            stacklevel=1,
        )


# --------------------------------------------------------------------------- C2


@pytest.mark.calibration_full
def test_pooled_slope_ci_coverage_at_four_log2_recorded() -> None:
    """4 log2 (1.73 background SD): recorded only -- outside the measured regime."""
    R = C.replicates(SUITE, full=True)
    assert amplitude_sd(4.0, RHO) == pytest.approx(1.732, abs=1e-3)
    rec = _coverage("gradient4log2", R, n_genes=FULL_GENES, amplitude=4.0)
    assert 0.0 <= rec["coverage"] <= 1.0
    assert rec["estimable"] >= 0.95, rec


# --------------------------------------------------------------------------- C3 teeth


@pytest.mark.calibration_full
def test_teeth_single_tu_resampling_covers_worse() -> None:
    """Teeth: with ``block_len = 1`` the CI ignores the AR(1) dependence and under-covers."""
    R = C.replicates(SUITE, full=False)
    blocks = _coverage("teeth_blocks", R, n_genes=1800, amplitude=AMPLITUDE)
    iid = _coverage("teeth_iid", R, n_genes=1800, amplitude=AMPLITUDE, block_len=1)
    assert iid["coverage"] < blocks["coverage"], (iid, blocks)
    assert iid["mean_ci_width"] < blocks["mean_ci_width"], (iid, blocks)
