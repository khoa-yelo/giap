"""Calibration of the gippasim null scenarios (issue M5-04; GIPPA_SPEC.md §8.4 "Autocorrelation").

The two nulls the M5 autocorrelation tests are calibrated on must hold exactly in the generator:

* ``label_null``: gene abundance iid normal (no TU, no position, no AR(1) structure), so the
  Pearson correlation of disjoint neighbour pairs along chr1 (genes 2j, 2j+1 in genome order) and of
  disjoint within-TU pairs (the first two genes of each TU) has an exact t(m - 2) null;
* ``tu_null``: TU-structured noise without inter-TU structure (rho = 0), so after removing the
  known position decay (``abundance + decay * tu_position``), the last gene of TU 2k and the first
  gene of TU 2k+1 (genome order, each TU used once) are independent iid normals: exact t test.

Each replicate's p-value is U(0, 1) when the scenario follows the spec; gated through the M4-05
harness (KS uniformity and the level bound at every alpha, ADR-0015 two-block confirmation) in
suite ``autocorrelation`` (the generator of that suite's nulls). Replicate ``i`` draws layout and
expression from ``replicate_rng(stream, i)``.
"""

from __future__ import annotations

import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import expression_oracle as eo

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
SUITE = "autocorrelation"
pytestmark = pytest.mark.calsuite(SUITE)


def _run(cfgs: tuple[Any, Any], rng: np.random.Generator) -> tuple[Any, Any]:
    from gippasim.expression import simulate_expression
    from gippasim.layout import simulate_layout

    lay = simulate_layout(cfgs[0], rng)
    return lay, simulate_expression(lay, cfgs[1], rng)


def _r_pvalue(x: np.ndarray, y: np.ndarray) -> float:
    """Two-sided p of Pearson r = 0 (exact t(m - 2) for independent normal pairs)."""
    m = len(x)
    r = float(np.corrcoef(x, y)[0, 1])
    t = r * np.sqrt((m - 2) / max(1e-300, 1 - r * r))
    return float(2 * stats.t.sf(abs(t), m - 2))


def _chr1(lay: Any, ex: Any, decay: float = 0.0) -> tuple[list[str], list[float]]:
    """chr1 genes in genome order: TU ids and abundance + decay * TU position."""
    pos = eo.gene_tu_position(lay)
    ab = dict(zip(ex.genes["locus_tag"], ex.genes["abundance_log2_true"], strict=True))
    g = lay.genes[lay.genes["contig"] == "chr1"].sort_values("start", kind="mergesort")
    tags = g["locus_tag"].astype(str).tolist()
    return [pos[t][0] for t in tags], [float(ab[t]) + decay * pos[t][1] for t in tags]


def _gate(stream: str, full: bool, fn: Callable[[np.random.Generator], float]) -> None:
    memo: dict[int, float] = {}

    def once(rng: np.random.Generator) -> float:
        key = int(rng.bit_generator.seed_seq.entropy)  # type: ignore[attr-defined]
        if key not in memo:
            memo[key] = float(fn(rng))
        return memo[key]

    R = C.replicates(SUITE, full)
    statuses = [C.calibrate(SUITE, "uniform", 0.0, once, R, stream=stream)["status"]]
    for a in C.load_tolerances().levels:
        statuses.append(C.calibrate(SUITE, "level", a, once, R, stream=stream)["status"])
    assert "fail" not in statuses, (stream, statuses)


@pytest.mark.parametrize("full", MODES)
def test_label_null_neighbour_pairs(full: bool) -> None:
    def fn(rng: np.random.Generator) -> float:
        import gippasim.scenarios as sc

        lay, ex = _run(sc.label_null(1000), rng)
        _, x = _chr1(lay, ex)
        k = len(x) // 2 * 2
        a = np.asarray(x[:k], float)
        return _r_pvalue(a[0::2], a[1::2])

    _gate("sim_scenarios.label_null.neighbours", full, fn)


@pytest.mark.parametrize("full", MODES)
def test_label_null_within_tu_pairs(full: bool) -> None:
    def fn(rng: np.random.Generator) -> float:
        import gippasim.scenarios as sc

        lay, ex = _run(sc.label_null(1000), rng)
        tus, x = _chr1(lay, ex)
        first: dict[str, list[float]] = {}
        for t, v in zip(tus, x, strict=True):
            first.setdefault(t, []).append(v)
        pairs = np.array([v[:2] for v in first.values() if len(v) >= 2])
        return _r_pvalue(pairs[:, 0], pairs[:, 1])

    _gate("sim_scenarios.label_null.within_tu", full, fn)


@pytest.mark.parametrize("full", MODES)
def test_tu_null_between_tu_pairs(full: bool) -> None:
    def fn(rng: np.random.Generator) -> float:
        import gippasim.scenarios as sc

        cfgs = sc.tu_null(1000)
        lay, ex = _run(cfgs, rng)
        tus, x = _chr1(lay, ex, decay=cfgs[1].position_decay_log2)
        runs: list[list[float]] = []
        for i, (t, v) in enumerate(zip(tus, x, strict=True)):
            if i == 0 or t != tus[i - 1]:
                runs.append([])
            runs[-1].append(v)
        a = np.array([runs[k][-1] for k in range(0, len(runs) - 1, 2)])
        b = np.array([runs[k + 1][0] for k in range(0, len(runs) - 1, 2)])
        return _r_pvalue(a, b)

    _gate("sim_scenarios.tu_null.between_tu", full, fn)
