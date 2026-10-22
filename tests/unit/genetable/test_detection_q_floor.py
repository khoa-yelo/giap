"""Khoa 2026-10-22: the picked-FDR floor 1/T (T = target-win genes) is reported as the ``tier`` row
``q_floor``; when it exceeds ``stats.alpha.identification`` the run warns ``q_floor_above_alpha``
and the gene q-values and statuses are left exactly as computed (no alpha change)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detection_tiers_fixtures as fx
from detection_tiers_fixtures import Builder


def _picked(n_genes: int, n_decoys: int) -> fx.Run:
    b = Builder(n_genes)
    t = b.tags or []
    for i, g in enumerate(t):
        b.target(g, 0.001, 0.01, 100.0 + i)
    for i in range(n_decoys):
        b.decoy([f"rev_{t[i]}"], 0.001, 0.01, 1.0 + i)
    return b.run()


@pytest.mark.parametrize(("n_genes", "flagged"), [(85, True), (100, False), (150, False)])
def test_q_floor_row_and_flag(n_genes: int, flagged: bool) -> None:
    r = _picked(n_genes, 8)
    assert r.det.tier == "decoy_picked"
    q = r.det.table["gene_q_value"].astype(float).to_numpy()
    n_win = int(np.isfinite(q).sum())
    assert fx.tier_metric(r.det, "q_floor") == pytest.approx(1.0 / n_win)
    assert np.nanmin(q) >= 1.0 / n_win - 1e-15
    assert ("q_floor_above_alpha" in r.det.warnings) is flagged


def test_flag_leaves_q_and_status_unchanged() -> None:
    """The 85-gene case: q stays at its floor 1/85, nothing is detected at alpha 0.01."""
    r = _picked(85, 8)
    assert r.det.tier == "decoy_picked"
    q = r.det.table["gene_q_value"].astype(float)
    assert q.min() == pytest.approx(1 / 85)
    assert "q_floor_above_alpha" in r.det.warnings
    assert not r.det.table["detected"].astype(bool).any()
    assert set(r.det.table["detection_status"]) == {"not_detected"}
