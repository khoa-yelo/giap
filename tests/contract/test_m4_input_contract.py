"""M4-06 contract conformance: simulator output now, real M2/M3 outputs later (M4-21).

Acceptance: "``tests/contract/test_m4_input_contract.py`` validates the M4-07 simulator output
(example 1) and, once M2/M3 merge, their real outputs" (M4-21 unmarked ``pending_sample``), and
"``TABLES`` is the single source of the column names; a test asserts [gippasim's] column lists equal
``TABLES``". The gippasim parts skip until M4-07 adds
``gippasim.layout``/``gippasim.contract_tables``.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "unit" / "genetable"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
from m4_contract_fixture import build_frames, empty_like, write_work

M3_KEYS = ("peptide_hits", "protein_hits", "unmapped")
SEED = 20261014
#: the required M4-06 tables (C1-C8, C11-C13) and the M1-12 peptides the real run writes
REAL_KEYS = ("replicons", "genes", "xrefs", "annotation_meta", "tus", "layers", "regions",
             "region_genes", "peptides", "peptide_hits", "protein_hits", "unmapped")  # fmt: skip


def _sim() -> tuple[object, object]:
    layout = pytest.importorskip("gippasim.layout", reason="M4-07 not merged yet")
    tables = pytest.importorskip("gippasim.contract_tables", reason="M4-07 not merged yet")
    return layout, tables


def test_sim_example1_200_genes_valid(tmp_path: Path) -> None:
    """Example 1: ``gippasim.write_contract_tables`` for a 200-gene layout loads with no errors."""
    import numpy as np

    layout_mod, tables_mod = _sim()
    from gippa.io.work_tables import default_paths, load_m4_inputs

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
    """Drift check both ways: gippasim's literal ``CONTRACT_COLUMNS`` equal ``TABLES``."""
    _, tables_mod = _sim()
    from gippa.genetable.contract import TABLES

    sim_cols = tables_mod.CONTRACT_COLUMNS  # type: ignore[attr-defined]
    assert set(sim_cols) <= set(TABLES)
    for key, cols in sim_cols.items():
        assert tuple(cols) == tuple(c.name for c in TABLES[key]), key


def test_m1_12_peptides_conform_to_c9(tmp_path: Path) -> None:
    """C9 against the merged M1-12 producer: ``aggregate_peptides`` writes ``protein_ids`` as a
    native list (lead decision 1); its ``peptides.norm.parquet`` loads with no peptides error."""
    from gippa.config.model import InputConfig
    from gippa.io.peptides import aggregate_peptides, read_peptide_table, validate_peptides
    from gippa.io.work_tables import default_paths, load_m4_inputs

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
    agg.to_parquet(paths["peptides"], engine="pyarrow", index=False)  # as gippa.preflight.run
    inputs, issues = load_m4_inputs(paths)
    bad = [i for i in issues if i.severity == "error" and i.code.startswith("contract.peptides.")]
    assert bad == []
    assert inputs.peptides is not None
    ids = dict(zip(inputs.peptides["peptide"], inputs.peptides["protein_ids"], strict=True))
    assert [str(x) for x in ids["AAAAAAAK"]] == ["GSIM_00025", "GSIM_00030"]
    assert len(ids["CCCCCCCR"]) == 0


# --------------------------------------------------------------------------- M4-21 real path
# Examples 1-3 of docs/specs/M4-21-m2-m3-conformance.md on the outputs of the merged M2 and M3
# steps (``gippa run --executor local`` on gippasim SIM100K with its precomputed layers,
# tests/helpers/m4_real.py). The M4-06 contract is the arbiter (M4-21 lead decision 1): a failure
# here is a producer deviation unless the lead rules the contract wrong.


@pytest.fixture(scope="module")
def real(tmp_path_factory: pytest.TempPathFactory) -> tuple[Any, Path]:
    import m4_chain as mc
    import m4_real as mr

    root = tmp_path_factory.mktemp("m4_real")
    sim = mr.sim100k(root / "sim")
    return sim, mr.gippa_run(sim, root / "OUT", mc.test_db(root / "db"), SEED)


def _load(work: Path) -> tuple[Any, list[Any]]:
    import m4_real as mr

    from gippa.io.work_tables import load_m4_inputs

    return load_m4_inputs(mr.real_paths(work))


@pytest.mark.parametrize("key", REAL_KEYS)
def test_real_output_conforms(real: tuple[Any, Path], key: str) -> None:
    """Example 1 (and every other required table): the real producer's file loads with no
    ``contract.<key>.*`` error. ``regions``/``region_genes`` are C7/C8 of M2-19 ``region_layer``,
    required by ``context_stats`` and ``sensitivity`` even when the regions come from
    ``annotation.precomputed.regions``."""
    _, issues = _load(real[1])
    bad = [(i.code, i.message) for i in issues
           if i.severity == "error" and i.code.startswith(f"contract.{key}.")]  # fmt: skip
    assert bad == [], f"producer of {key} deviates from M4-06: {bad}"


def test_real_outputs_have_no_contract_error(real: tuple[Any, Path]) -> None:
    """Cross-table invariants (``validate_inputs``) included."""
    _, issues = _load(real[1])
    assert [(i.code, i.message) for i in issues if i.severity == "error"] == []


def test_example_1_real_genes_validate(real: tuple[Any, Path]) -> None:
    """Example 1: ``validate_table("genes", ...)`` on M2's ``genes`` returns no errors."""
    import pandas as pd

    from gippa.genetable.contract import validate_table

    genes = pd.read_parquet(real[1] / "genes")
    assert [i.code for i in validate_table("genes", genes) if i.severity == "error"] == []
    assert len(genes) > 50, len(genes)


def test_example_2_real_tus_partition_the_genes(real: tuple[Any, Path]) -> None:
    """Example 2 (M4-06 Q1): every gene is in exactly one TU, singletons included."""
    import pandas as pd

    work = real[1]
    genes = pd.read_parquet(work / "genes")["locus_tag"].astype(str)
    members = [str(g) for gs in pd.read_parquet(work / "tus")["genes"] for g in gs]
    counts = pd.Series(members).value_counts()
    assert sorted(counts.index) == sorted(genes), "TU members differ from the gene set"
    assert counts.max() == 1, f"genes in more than one TU: {list(counts[counts > 1].index)[:5]}"


def test_example_3_real_peptide_hits_assign_and_are_accepted(real: tuple[Any, Path]) -> None:
    """Example 3: M3's ``peptide_hits`` from the simulated proteome's tryptic peptides --
    assigning classes carry ``gene`` (no ``contract.peptide_hits.*`` error) and M4-11
    ``accepted_peptides`` accepts >= 99.9 % of the planted peptides.

    Planted = correct (``is_wrong`` false) gene peptides with exactly one exact genome locus and
    no alternative-start N-terminal Met (``n_loci_genome == 1``, ``alt_start_nterm`` false, as
    tests/helpers/placement_sims.py): the others are shared or unplaceable by design."""
    import pandas as pd

    from gippa.config import resolve_config
    from gippa.genetable.evidence import accepted_peptides

    sim, work = real
    inputs, issues = _load(work)
    assert not [i for i in issues if i.code.startswith("contract.peptide_hits.")
                and i.severity == "error"]  # fmt: skip
    truth = pd.read_csv(sim.truth["peptides"], sep="\t", dtype=str, keep_default_na=False)
    planted = truth[(truth["origin"] == "gene") & (truth["is_wrong"].str.lower() == "false")
                    & (truth["n_loci_genome"] == "1")
                    & (truth["alt_start_nterm"].str.lower() == "false")]  # fmt: skip
    assert len(planted) > 1000, len(planted)
    acc = accepted_peptides(inputs.peptide_hits, inputs.peptides, resolve_config(sim.config_yaml))
    pairs = set(zip(acc["gene"].astype(str), acc["peptide"].astype(str), strict=True))
    missed = [(g, p) for g, p in zip(planted["locus_tag"], planted["peptide"], strict=True)
              if (g, p) not in pairs]  # fmt: skip
    assert 1 - len(missed) / len(planted) >= 0.999, f"{len(missed)} missed, e.g. {missed[:5]}"
