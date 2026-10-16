"""Attack null for the M4-16 basic O - E CI in a low-E group (PR #87 lead decision 4; §8.4).

The statistician's review (local/scratch/review-m4-16/rare.py) showed the TU-block percentile CI
of ``o_minus_e`` under-covers badly when the group expects few detections (or few
non-detections), even with many TUs: with 50 single-gene TUs at p = 0.02 (E = 1) the 95 % CI misses
the true value 0 in ~37 % of null replicates (O = 0 makes every TU sum equal and the CI
degenerate at -E). Decision 4 makes ``basic_rows`` flag such rows ``ci_few_units`` when
min(E, n_genes - E) < 10, in addition to < 10 TUs.

Attacks (all with >= 50 TUs, so only the E rule can set the flag; detections drawn
comonotonically inside a TU, so P(d_g = 1) = p_g and E[O - E] = 0 exactly; group = the
``sorf`` feature class):

- ``sorf_e1``: 50 single-gene TUs, p = 0.02 (E = 1).
- ``sorf_e2_operons``: 50 two-gene TUs, p = 0.02 (E = 2).
- ``saturated_n_minus_e1``: 50 single-gene TUs, p = 0.98 (n - E = 1).

Gates, per attack:
1. ``test_low_e_group_is_flagged``: every replicate's sORF ``o_minus_e`` row carries
   ``ci_few_units`` (deterministic; R = ``replicates_*.regions`` draws of the null).
2. ``test_low_e_unflagged_miss_rate_exceeds_level`` (teeth, PR #69 decision 3): the CI, read as if
   unflagged, misses 0 more often than the level allows: the ``fwer`` criterion on 0/1 miss
   indicators at 0.05 (bound 0.05 + k sqrt(0.05 * 0.95 / R)) is a confirmed ``fail`` through
   ``_criteria.calibrate`` (stream ``regions.attack_detection_summary.<attack>``). This is why the
   flag is required: the miss rate is far above the level, not a Monte Carlo fluctuation.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import _criteria as C
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detection_summary_fixture as fx

SUITE = "regions"
B = 1000
MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]


@dataclass(frozen=True)
class Attack:
    n_tu: int
    tu_size: int
    p: float


ATTACKS = {
    "sorf_e1": Attack(50, 1, 0.02),
    "sorf_e2_operons": Attack(50, 2, 0.02),
    "saturated_n_minus_e1": Attack(50, 1, 0.98),
}


def null_genes(a: Attack, rng: np.random.Generator) -> pd.DataFrame:
    """sORF group of ``n_tu`` TUs of ``tu_size`` genes; one uniform per TU (comonotone)."""
    n = a.n_tu * a.tu_size
    u = np.repeat(rng.random(a.n_tu), a.tu_size)
    d = u < a.p
    return pd.DataFrame(
        {
            "locus_tag": [f"S{i:04d}" for i in range(n)],
            "contig": "chr1",
            "feature_type": "sorf",
            "hypothetical": np.zeros(n, dtype=bool),
            "partial": np.zeros(n, dtype=bool),
            "tu_id": np.repeat([f"TU{j:04d}" for j in range(a.n_tu)], a.tu_size),
            "detected": d,
            "detection_status": np.where(d, "detected", "not_detected").astype(object),
            "p_detect": np.full(n, a.p, dtype=np.float64),
            "cog_category": pd.Series([None] * n, dtype=object),
        }
    )


def sorf_row(a: Attack, rng: np.random.Generator) -> pd.Series:
    from giap.context.summary import basic_rows

    genes = null_genes(a, rng)
    out = basic_rows(genes, B, np.random.default_rng(rng.integers(2**63)))
    assert int(fx.metric(out, "feature_class", "sorf", "n_tu")["value"]) == a.n_tu >= 10
    return fx.metric(out, "feature_class", "sorf", "o_minus_e")


@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize("attack", sorted(ATTACKS))
def test_low_e_group_is_flagged(attack: str, full: bool) -> None:
    a = ATTACKS[attack]
    R = C.replicates(SUITE, full)
    unflagged = []
    for i in range(R):
        r = sorf_row(a, C.replicate_rng(f"{SUITE}.attack_detection_summary.flag.{attack}", i))
        if fx.FEW_UNITS_FLAG not in fx.flag_set(r["flags"]):
            unflagged.append(i)
    assert not unflagged, f"{attack}: {len(unflagged)}/{R} low-E rows without ci_few_units"


@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize("attack", sorted(ATTACKS))
def test_low_e_unflagged_miss_rate_exceeds_level(attack: str, full: bool, teeth_log: Path) -> None:
    a = ATTACKS[attack]
    R = C.replicates(SUITE, full)

    def miss(rng: np.random.Generator) -> float:
        r = sorf_row(a, rng)
        lo, hi = float(r["ci_low"]), float(r["ci_high"])
        assert np.isfinite(lo) and np.isfinite(hi), r
        return 0.0 if lo <= 0.0 <= hi else 1.0

    rec = C.calibrate(SUITE, "fwer", 0.05, miss, R,
                      stream=f"{SUITE}.attack_detection_summary.{attack}",
                      log_dir=teeth_log)  # fmt: skip
    print(f"low-E attack {attack} (E = {a.n_tu * a.tu_size * a.p:g}): {rec}")
    assert rec["status"] == "fail", f"{attack}: unflagged miss rate not a confirmed fail: {rec}"
