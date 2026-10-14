"""M4-06 contract conformance: simulator output now, real M2/M3 outputs later (M4-21).

Acceptance: "``tests/contract/test_m4_input_contract.py`` validates the M4-07 simulator output
(example 1) and, once M2/M3 merge, their real outputs (marked ``pending_sample`` until then)", and
"``TABLES`` is the single source of the column names; a test asserts [giapsim's] column lists equal
``TABLES``". The giapsim parts skip until M4-07 adds ``giapsim.layout``/``giapsim.contract_tables``.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "unit" / "genetable"))
from m4_contract_fixture import build_frames, empty_like, write_work

REAL_WORK = os.environ.get("GIAP_M4_CONFORMANCE_WORK", "")
M3_KEYS = ("peptide_hits", "protein_hits", "unmapped")


def _sim() -> tuple[object, object]:
    layout = pytest.importorskip("giapsim.layout", reason="M4-07 not merged yet")
    tables = pytest.importorskip("giapsim.contract_tables", reason="M4-07 not merged yet")
    return layout, tables


def test_sim_example1_200_genes_valid(tmp_path: Path) -> None:
    """Example 1: ``giapsim.write_contract_tables`` for a 200-gene layout loads with no errors."""
    import numpy as np

    layout_mod, tables_mod = _sim()
    from giap.io.work_tables import default_paths, load_m4_inputs

    work = tmp_path / "work"
    layout = layout_mod.simulate_layout(  # type: ignore[attr-defined]
        layout_mod.LayoutConfig(n_genes_chromosome=200),  # type: ignore[attr-defined]
        np.random.default_rng(1),
    )
    tables_mod.write_contract_tables(layout, work)  # type: ignore[attr-defined]
    paths = default_paths(work)
    # M4-07 writes the M2-side files (C1-C8); M3's tables are added header-only (nothing placed).
    stub = build_frames(20)
    stub = {k: empty_like(stub[k]) for k in M3_KEYS}
    write_work(stub, {k: paths[k] for k in M3_KEYS if not paths[k].exists()})
    inputs, issues = load_m4_inputs(paths)
    assert [i for i in issues if i.severity == "error"] == []
    assert len(inputs.genes) == 200


def test_sim_columns_equal_tables() -> None:
    """Drift check both ways: giapsim's literal ``CONTRACT_COLUMNS`` equal ``TABLES``."""
    _, tables_mod = _sim()
    from giap.genetable.contract import TABLES

    sim_cols = tables_mod.CONTRACT_COLUMNS  # type: ignore[attr-defined]
    assert set(sim_cols) <= set(TABLES)
    for key, cols in sim_cols.items():
        assert tuple(cols) == tuple(c.name for c in TABLES[key]), key


@pytest.mark.pending_sample
@pytest.mark.skipif(
    not REAL_WORK or not Path(REAL_WORK).is_dir(),
    reason="pending_sample: no real M2/M3 work directory (GIAP_M4_CONFORMANCE_WORK; M4-21)",
)
def test_real_m2_m3_outputs_conform() -> None:
    from giap.io.work_tables import default_paths, load_m4_inputs

    _, issues = load_m4_inputs(default_paths(Path(REAL_WORK)))
    assert [i for i in issues if i.severity == "error"] == []
