"""Performance budgets for M4-11 accepted evidence (docs/specs/M4-11-accepted-evidence.md: "runs
in < 1 s for 4,500 genes"); ``accepted_peptides`` + ``gene_evidence``.

- M4-09 default giapsim fixture (layout seed 1, expression seed 2, evidence seed 3): one run
  < 1 s per 4,500 genes (floor 1 s). From ``tests/unit/genetable/test_evidence.py::
  test_m409_default_fixture_detected_true_have_evidence``.
- 4,500-gene chromosome + plasmid + linear contig, origin crossing (``_sim_inputs(4500, 5,
  True)``): best of 3 < 1 s. From ``test_evidence.py::test_simulated_4500_genes_under_one_second``
  (now ``test_simulated_4500_genes_detected_have_evidence``).

The unit tests keep the correctness checks; the wall-clock asserts were flaky on loaded Slurm
nodes. Not part of the per-commit suite; run with
``pytest benchmarks --benchmark-only``. Recorded in ``benchmarks/baseline.json``.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from ._unit_helpers import load_test_module, timed

_E = load_test_module("unit/genetable/test_evidence.py")


def test_m409_default_fixture_evidence_budget(benchmark: Any) -> None:
    evmod = pytest.importorskip("giapsim.evidence")
    from giapsim.expression import ExpressionConfig, simulate_expression
    from giapsim.layout import LayoutConfig, simulate_layout

    layout = simulate_layout(LayoutConfig(), np.random.default_rng(1))
    expr = simulate_expression(layout, ExpressionConfig(), np.random.default_rng(2))
    sim = evmod.simulate_evidence(layout, expr, evmod.EvidenceConfig(), np.random.default_rng(3))
    ev = _E._ev()
    cfg = _E._cfg()

    def once() -> Any:
        acc = ev.accepted_peptides(sim.peptide_hits, sim.peptides, cfg)
        return ev.gene_evidence(layout.genes, sim.peptide_hits, acc, layout.replicons)

    def run() -> float:
        secs, out = timed(once)
        assert len(out) == len(layout.genes)
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    budget = 1.0 * max(1.0, len(layout.genes) / 4500)
    benchmark.extra_info.update(secs=secs, budget_s=budget, n_genes=len(layout.genes))
    assert secs < budget, f"{secs:.3f} s for {len(layout.genes)} genes (budget {budget:.3f} s)"


def test_simulated_4500_genes_evidence_budget(benchmark: Any) -> None:
    layout, _detected, peps, hits = _E._sim_inputs(4500, 5, True)
    ev = _E._ev()
    cfg = _E._cfg()

    def once() -> Any:
        acc = ev.accepted_peptides(hits, peps, cfg)
        return ev.gene_evidence(layout.genes, hits, acc, layout.replicons)

    def run() -> float:
        best = float("inf")
        for _ in range(3):
            secs, out = timed(once)
            assert len(out) == len(layout.genes) >= 4500
            best = min(best, secs)
        return best

    best = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info.update(best_secs=best, n_genes=len(layout.genes))
    assert best < 1.0, f"best of 3: {best:.3f} s for {len(layout.genes)} genes"
