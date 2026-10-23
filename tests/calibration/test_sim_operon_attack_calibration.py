"""Calibration sanity of the giapsim detection attack (issue M6-03, lead decisions 10 and 12 of
docs/specs/M6-03-sim-operon-scenarios.md; GIAP_SPEC.md §8.4 ``missed_operons``, §8.5, ADR-0018).

Replicate ``i`` of a scenario is ``simulate_operons(cfgs, suite_seed("m6_03_<name>", i))`` (the
stream names of tests/calibration/test_sim_operon_calibration.py), R = 20 fast, full R = the
``missed_operons`` suite's ``replicates_full``. Delta is the true ``b - a`` abundance difference
of pairs with both genes detected; targets are true NO pairs (class NO, not ``true_operon``), so
a UniOP false negative cannot fake the excess; decoys are the CV and SW pairs; distance bins are
the merged method's default edges (-inf, 0, 50, 150, 300, inf) on ``intergenic_bp``.

* **decision 12 (the attack bites)**: ``c`` is the 25th percentile of |Delta| of the null's
  (``operon_null(1800)``) targets pooled over its R replicates. In
  ``operon_attack_detection(1800, "prophage")`` the targets' share of |Delta| < c exceeds the
  decoys' share *in the same distance bin*: the stratified (Cochran-Mantel-Haenszel) one-sided z,
  pooled over the R replicates and the bins, is >= 2. Bins lacking targets or decoys contribute
  nothing; at least one bin below 50 bp must contribute (the short-gap decoys of decision 11);
* **decision 10**: pooled over the R replicates, >= 25 % of the block's adjacent (non-SW) pairs
  have both genes detected.

Interface as in tests/unit/sim/test_operon_scenarios_attacks.py.
"""

from __future__ import annotations

import math
from functools import cache
from typing import Any

import _criteria as C
import numpy as np
import pytest

R_FAST = 20
MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
N = 1800
EDGES = np.array([-np.inf, 0.0, 50.0, 150.0, 300.0, np.inf])
ATTACK = ("operon_attack_detection", (N, "prophage"))
NULL = ("operon_null", (N,))


def _R(full: bool) -> int:
    return C.replicates("missed_operons", True) if full else R_FAST


@cache
def _sim(name: str, args: tuple[Any, ...], rep: int) -> Any:
    import giapsim.scenarios as sc

    return sc.simulate_operons(getattr(sc, name)(*args), C.suite_seed(f"m6_03_{name}", rep))


def _b(s: Any) -> np.ndarray:
    return np.asarray(s.astype(bool).to_numpy(), dtype=bool)


def _both(p: Any) -> Any:
    return p[_b(p["detected_a"]) & _b(p["detected_b"])]


def _absdelta(s: Any, p: Any) -> np.ndarray:
    g = s.expression.genes
    ab = dict(zip(g["locus_tag"].astype(str), g["abundance_log2_true"].astype(float),
                  strict=True))  # fmt: skip
    return np.array([abs(ab[str(b)] - ab[str(a)]) for a, b in zip(p["gene_a"], p["gene_b"],
                                                                  strict=True)])  # fmt: skip


def _targets(p: Any) -> Any:
    return p[(p["class"] == "NO").to_numpy() & ~_b(p["true_operon"])]


def _decoys(p: Any) -> Any:
    return p[p["class"].isin(["CV", "SW"]).to_numpy()]


def _bin(p: Any) -> np.ndarray:
    d = p["intergenic_bp"].astype(float).to_numpy()
    return np.searchsorted(EDGES, d, side="right") - 1


def _null_cut(R: int) -> float:
    vals = [_absdelta(s, _targets(_both(s.pairs))) for s in (_sim(*NULL, i) for i in range(R))]
    return float(np.quantile(np.concatenate(vals), 0.25))


def _cmh(R: int, c: float) -> tuple[float, list[tuple[int, int, int, int]]]:
    """(stratified one-sided z, per-bin (targets small, targets, decoys small, decoys))."""
    nb = len(EDGES) - 1
    cells = np.zeros((nb, 4), dtype=np.int64)
    for i in range(R):
        s = _sim(*ATTACK, i)
        p = _both(s.pairs)
        for k, fr in ((0, _targets(p)), (2, _decoys(p))):
            small = _absdelta(s, fr) < c
            b = _bin(fr)
            cells[:, k] += np.bincount(b[small], minlength=nb)
            cells[:, k + 1] += np.bincount(b, minlength=nb)
    obs = exp = var = 0.0
    for k1, n1, k2, n2 in cells.tolist():
        tot, m = n1 + n2, k1 + k2
        if n1 == 0 or n2 == 0 or tot < 2:
            continue
        obs += k1
        exp += n1 * m / tot
        var += n1 * n2 * m * (tot - m) / (tot * tot * (tot - 1))
    z = (obs - exp) / math.sqrt(var) if var > 0 else -math.inf
    return z, [tuple(r) for r in cells.tolist()]


@pytest.mark.calsuite("missed_operons")
@pytest.mark.parametrize("full", MODES)
def test_detection_attack_creates_a_within_bin_small_delta_excess(full: bool) -> None:
    R = _R(full)
    c = _null_cut(R)
    z, cells = _cmh(R, c)
    short = [r for r in cells[:2] if r[1] and r[3]]
    assert short, cells  # some bin below 50 bp holds targets and decoys
    assert z >= 2.0, (z, c, cells)


@pytest.mark.calsuite("missed_operons")
@pytest.mark.parametrize("full", MODES)
def test_detection_attack_block_pairs_are_both_detected_at_least_a_quarter(full: bool) -> None:
    k = n = 0
    for i in range(_R(full)):
        s = _sim(*ATTACK, i)
        ids = set(s.layout.regions.loc[s.layout.regions["region_type"] == "prophage",
                                       "region_id"])  # fmt: skip
        rg = s.layout.region_genes
        block = set(rg.loc[rg["region_id"].isin(ids), "locus_tag"].astype(str))
        p = s.pairs[s.pairs["class"] != "SW"]
        inb = p["gene_a"].astype(str).isin(block) & p["gene_b"].astype(str).isin(block)
        bp = p[inb.to_numpy()]
        k += int((_b(bp["detected_a"]) & _b(bp["detected_b"])).sum())
        n += len(bp)
    assert n > 0 and k / n >= 0.25, (k, n)
