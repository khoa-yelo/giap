"""Calibration of the gippasim operon generators (issue M6-03; GIPPA_SPEC.md §8.4 ``discordance``
and ``missed_operons``, §7.4.7 "Calibration"; docs/specs/M6-03-sim-operon-scenarios.md acceptance).

The generators must do what M6-09 assumes of them, pooled over R replicate 1,800-gene genomes
(R = 20 fast as the spec says; full R = the suite's ``replicates_full``), replicate ``i`` drawn by
``simulate_operons(cfgs, suite_seed(<stream>, i))``; Delta is the true ``b - a`` abundance
difference of pairs with both genes detected:

* **missed-operon null** (``operon_null``): no excess of small |Delta| among true NO pairs (class
  NO, not ``true_operon``) relative to the decoys (CV and SW pairs). Per replicate the cut-off ``c``
  is the decoys' 10 % quantile of |Delta|; pooled, the NO share below ``c`` minus the decoy share
  has a two-proportion z <= 3 (one-sided: "no excess");
* **missed-operon power** (``missed_operon_power(1800, 50, 0.1)``): the same excess is there, with
  z >= 5 and at least 0.03 in share;
* **discordance null** (``operon_null``): the share of true OP pairs with |Delta - loc| > 2 scale
  under ``truth["operon_t"]`` equals the t tail ``2 sf_t(2; df)`` within 3 binomial SEs;
* **discordance power** (``discordance_power(1800, 40, 2.0)``): that tail share exceeds the t tail
  with z >= 4 (this assumes a generating scale below 1 log2, so the planted 2.0 lies beyond
  2 scale; the methods note's within-operon scale is 0.5).

Interface as in tests/unit/sim/test_operon_scenarios.py (``simulate_operons``, ``pairs`` with
``true_operon``, ``truth["operon_t"]``). Plain pooled assertions in the style of
test_sim_layout_calibration.py (generator facts, not p-value criteria of a method).
"""

from __future__ import annotations

import math
from functools import cache
from typing import Any

import _criteria as C
import numpy as np
import pytest
from scipy import stats

R_FAST = 20
MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
N = 1800


def _R(suite: str, full: bool) -> int:
    return C.replicates(suite, True) if full else R_FAST


def _sc() -> Any:
    import gippasim.scenarios

    return gippasim.scenarios


@cache
def _sim(name: str, args: tuple[Any, ...], rep: int) -> Any:
    sc = _sc()
    return sc.simulate_operons(getattr(sc, name)(*args), C.suite_seed(f"m6_03_{name}", rep))


def _b(s: Any) -> np.ndarray:
    return np.asarray(s.astype(bool).to_numpy(), dtype=bool)


def _delta(s: Any, p: Any) -> np.ndarray:
    g = s.expression.genes
    ab = dict(zip(g["locus_tag"].astype(str), g["abundance_log2_true"].astype(float),
                  strict=True))  # fmt: skip
    return np.array([ab[str(b)] - ab[str(a)] for a, b in zip(p["gene_a"], p["gene_b"],
                                                             strict=True)])  # fmt: skip


def _both(p: Any) -> Any:
    return p[_b(p["detected_a"]) & _b(p["detected_b"])]


def _small_counts(s: Any) -> tuple[int, int, int, int]:
    """(NO below c, NO total, decoys below c, decoys total); c = the decoy 10 % |Delta| quantile."""
    p = _both(s.pairs)
    no = p[(p["class"] == "NO").to_numpy() & ~_b(p["true_operon"])]
    dec = p[p["class"].isin(["CV", "SW"]).to_numpy()]
    d_no, d_dec = np.abs(_delta(s, no)), np.abs(_delta(s, dec))
    c = float(np.quantile(d_dec, 0.10))
    return int((d_no < c).sum()), len(d_no), int((d_dec < c).sum()), len(d_dec)


def _excess(name: str, args: tuple[Any, ...], R: int) -> tuple[float, float]:
    """Pooled (NO share - decoy share) below c and its two-proportion z."""
    k1 = n1 = k2 = n2 = 0
    for i in range(R):
        a, b, c, d = _small_counts(_sim(name, args, i))
        k1, n1, k2, n2 = k1 + a, n1 + b, k2 + c, n2 + d
    p1, p2 = k1 / n1, k2 / n2
    pool = (k1 + k2) / (n1 + n2)
    se = math.sqrt(pool * (1 - pool) * (1 / n1 + 1 / n2))
    return p1 - p2, (p1 - p2) / se


@pytest.mark.calsuite("missed_operons")
@pytest.mark.parametrize("full", MODES)
def test_null_no_excess_of_small_delta_among_no_pairs(full: bool) -> None:
    diff, z = _excess("operon_null", (N,), _R("missed_operons", full))
    assert z <= 3.0, (diff, z)


@pytest.mark.calsuite("missed_operons")
@pytest.mark.parametrize("full", MODES)
def test_power_scenario_has_an_excess_of_small_delta(full: bool) -> None:
    diff, z = _excess("missed_operon_power", (N, 50, 0.1), _R("missed_operons", full))
    assert z >= 5.0 and diff >= 0.03, (diff, z)


def _tail(name: str, args: tuple[Any, ...], R: int) -> tuple[int, int, float]:
    """(true OP pairs beyond 2 scale, true OP pairs, expected t tail share) pooled over R."""
    k = n = 0
    tail = math.nan
    for i in range(R):
        s = _sim(name, args, i)
        t = s.truth["operon_t"]
        loc, scale, df = (float(t[c].iloc[0]) for c in ("loc", "scale", "df"))
        tail = float(2 * stats.t.sf(2.0, df))
        p = _both(s.pairs)
        op = p[(p["class"] == "OP").to_numpy() & _b(p["true_operon"])]
        d = _delta(s, op)
        k += int((np.abs(d - loc) > 2.0 * scale).sum())
        n += len(d)
    return k, n, tail


@pytest.mark.calsuite("discordance")
@pytest.mark.parametrize("full", MODES)
def test_null_op_tail_matches_the_generating_t(full: bool) -> None:
    k, n, tail = _tail("operon_null", (N,), _R("discordance", full))
    se = math.sqrt(tail * (1 - tail) / n)
    assert abs(k / n - tail) <= 3.0 * se, (k, n, tail)


@pytest.mark.calsuite("discordance")
@pytest.mark.parametrize("full", MODES)
def test_discordance_power_has_an_excess_tail(full: bool) -> None:
    k, n, tail = _tail("discordance_power", (N, 40, 2.0), _R("discordance", full))
    se = math.sqrt(tail * (1 - tail) / n)
    assert (k / n - tail) / se >= 4.0, (k, n, tail)
