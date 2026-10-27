"""Attack null 4 for the M4-18 plasmid rows (GIPPA_SPEC.md §8.5 attack 4 "High-copy plasmids →
region
tests, reported as replicon-level effects only"; docs/specs/M4-18-regions-stats-table.md acceptance
"Attack null 4" and lead decisions 2 and 6).

Layout and expression as the M4-17 ``regions`` suite (tests/helpers/region_tests_data.py: M4-07
default regions + 60-gene plasmid ``p1`` with its own AMR set, origin crossing, 1,000 chromosome
genes; M4-08 AR(1) rho = 0.5, detection from covariates only, ``p = p_detect_true``).

- Gated (deterministic, like the M4-17 A4 invariance): at ``plasmid_copy_log2 = {"p1": log2(20)}``
  no plasmid row has a p-value (``null_model = none``, ``n_perm = 0``), and on the same seed the
  plasmid HL and its CI equal those at copy number 1 plus exactly log2(20), the O - E row is
  unchanged. A per-replicon residualised HL (HL ~ 0 whatever the dosage) or a plasmid tested inside
  a concatenated chromosome + plasmid index both fail. The chromosome region tests' level under
  this attack is M4-17's A4 test (tests/calibration/attacks/test_region_tests_attacks.py).
- Recorded, not gated (lead decisions 15-16): bias of the plasmid HL and coverage of its 95 %
  TU-block CI against the dosage log2(20) and against the realised dosage-plus-regulation shift (HL
  of the simulator's true abundance over eligible genes, the estimand), over R =
  ``replicates_fast.regions`` replicates; a warning when the realised-shift coverage is < 0.93.
  §8.4 and tolerances.yaml carry no criterion for this CI; the numbers are printed for the
  methods note.
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

import _criteria as C
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import region_tests_data as rd
import regions_stats_data as sd

pytestmark = pytest.mark.calibration
COPY_20 = math.log2(20.0)
STREAM = "regions.regions_stats.a4_plasmid"


def _sim(seed: int, copy_log2: float, B: int = 50) -> tuple[pd.DataFrame, float]:
    """(plasmid rows, realised shift) of one replicate. ``rd.simulate`` + ``rd.tables_from``
    consume the generator exactly as ``rd.sim_tables`` does. The realised shift (lead decision 16)
    is the HL of the simulator's ``abundance_log2_true`` over all eligible genes, plasmid against
    the reference replicons: the dosage-plus-regulation shift the CI estimates."""
    from gippa.context import regions_table as T
    from gippa.rng import for_step

    rng = np.random.default_rng(seed)
    kw = {"ar1_rho": 0.5, "plasmid_copy_log2": {"p1": copy_log2}}
    lay, expr = rd.simulate(rng, rd.calibration_layout_config(1000), kw)
    t = rd.tables_from(lay, expr, "true", rng)
    genes = t.genes.copy()
    genes["abundance_imputed"] = False
    s = int(rng.integers(0, 2**31 - 1))
    rows = T.plasmid_rows(genes, t.replicons, None, rd.config(bootstrap=B),
                          for_step("context_stats", s))  # fmt: skip
    truth = expr.genes.set_index("locus_tag")["abundance_log2_true"]
    a = truth.reindex(genes["locus_tag"].astype(str)).to_numpy(float)
    ref, _ = sd.reference_contigs(t.replicons, ["p1"])
    contig = genes["contig"].astype(str)
    el = sd.eligible(genes) & np.isfinite(a)
    x = a[el & (contig == "p1").to_numpy()]
    y = a[el & contig.isin(ref).to_numpy()]
    realised = float(np.median(np.subtract.outer(x, y)))
    return rows, realised


def _plasmid_rows(seed: int, copy_log2: float, B: int = 50) -> pd.DataFrame:
    return _sim(seed, copy_log2, B)[0]


@pytest.mark.parametrize("seed", [401, 402, 403])
def test_a4_plasmid_rows_untested_and_carry_the_dosage(seed: int) -> None:
    a = _plasmid_rows(seed, 0.0)
    b = _plasmid_rows(seed, COPY_20)
    for res in (a, b):
        assert set(res["region_id"]) == {"plasmid_1"}
        assert (res["null_model"] == "none").all()
        assert res["p_value"].isna().all() and res["q_value"].isna().all()
        assert (res["n_perm"].astype(int) == 0).all()
    ha, hb = sd.unadjusted(a, "plasmid_1", "hl_shift"), sd.unadjusted(b, "plasmid_1", "hl_shift")
    assert not sd.na(ha["effect"])
    for col in ("effect", "ci_low", "ci_high"):
        assert float(hb[col]) == pytest.approx(float(ha[col]) + COPY_20, abs=1e-9), col
    assert float(hb["effect"]) > 2.0  # the dosage is in the contrast, not residualised away
    oa, ob = sd.unadjusted(a, "plasmid_1", "o_minus_e"), sd.unadjusted(b, "plasmid_1", "o_minus_e")
    for col in ("effect", "ci_low", "ci_high", "observed", "expected"):
        assert float(oa[col]) == float(ob[col]), col


def test_recorded_plasmid_hl_bias_and_ci_coverage(full: bool) -> None:
    """Recorded, not gated (lead decisions 15-16): bias and 95 % CI coverage against log2(20)
    and against the realised shift; warns when the realised-shift coverage is below 0.93."""
    R = C.replicates("regions", full)
    err, err_real, cover, cover_real = [], [], 0, 0
    for i in range(R):
        seed = C.suite_seed(STREAM, i)
        rows, realised = _sim(seed, COPY_20, B=200)
        hl = sd.unadjusted(rows, "plasmid_1", "hl_shift")
        e, lo, hi = float(hl["effect"]), float(hl["ci_low"]), float(hl["ci_high"])
        err.append(e - COPY_20)
        err_real.append(e - realised)
        cover += lo <= COPY_20 <= hi
        cover_real += lo <= realised <= hi
    cov, cov_real = cover / R, cover_real / R
    msg = (
        f"recorded: plasmid HL at copy 20 (R={R}): mean HL - log2(20) = {np.mean(err):+.3f} "
        f"(sd {np.std(err, ddof=1):.3f}), coverage of log2(20) {cov:.3f}; "
        f"mean HL - realised shift = {np.mean(err_real):+.3f} "
        f"(sd {np.std(err_real, ddof=1):.3f}), coverage of the realised shift {cov_real:.3f}"
    )
    print(msg)
    if cov_real < 0.93:
        warnings.warn(msg, stacklevel=1)
    assert all(math.isfinite(x) for x in err + err_real)
