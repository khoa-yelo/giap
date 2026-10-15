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


def test_m1_12_peptides_conform_to_c9(tmp_path: Path) -> None:
    """C9 against the merged M1-12 producer: ``aggregate_peptides`` writes ``protein_ids`` as a
    native list (lead decision 1); its ``peptides.norm.parquet`` loads with no peptides error."""
    from giap.config.model import InputConfig
    from giap.io.peptides import aggregate_peptides, read_peptide_table, validate_peptides
    from giap.io.work_tables import default_paths, load_m4_inputs

    tsv = tmp_path / "peptides.tsv"
    tsv.write_text(
        "peptide\tprotein_ids\tq_value\trun\tintensity\n"
        "AAAAAAAK\tGSIM_00025;GSIM_00030\t0.001\tr1\t100\n"
        "AAAAAAAK\tGSIM_00030\t0.002\tr2\t50\n"
        "CCCCCCCR\tNA\t0.01\tr1\t1\n"
        "KAAAAAAA\trev_GSIM_00005\t0.3\tr1\t3\n",
        encoding="utf-8",
    )
    raw, _ = read_peptide_table(tsv)
    df, _ = validate_peptides(raw, InputConfig())
    agg, _ = aggregate_peptides(df)
    paths = default_paths(tmp_path / "work")
    write_work(build_frames(20), paths)
    agg.to_parquet(paths["peptides"], engine="pyarrow", index=False)  # as giap.preflight.run
    inputs, issues = load_m4_inputs(paths)
    bad = [i for i in issues if i.severity == "error" and i.code.startswith("contract.peptides.")]
    assert bad == []
    assert inputs.peptides is not None
    ids = dict(zip(inputs.peptides["peptide"], inputs.peptides["protein_ids"], strict=True))
    assert [str(x) for x in ids["AAAAAAAK"]] == ["GSIM_00025", "GSIM_00030"]
    assert len(ids["CCCCCCCR"]) == 0


@pytest.mark.pending_sample
@pytest.mark.skipif(
    not REAL_WORK or not Path(REAL_WORK).is_dir(),
    reason="pending_sample: no real M2/M3 work directory (GIAP_M4_CONFORMANCE_WORK; M4-21)",
)
def test_real_m2_m3_outputs_conform() -> None:
    from giap.io.work_tables import default_paths, load_m4_inputs

    _, issues = load_m4_inputs(default_paths(Path(REAL_WORK)))
    assert [i for i in issues if i.severity == "error"] == []
