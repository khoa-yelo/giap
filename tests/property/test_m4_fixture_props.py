"""Hypothesis properties of the M4-20 fixture builder ``giapsim.fixtures.build_m4_fixture``.

For any seed and any ``EvidenceConfig`` the fixture must be a usable M4 input: the C1-C14 files
exist in both layouts the steps read (flat for ``run_context_stats`` / ``run_sensitivity``,
``mapping/`` for ``gene_table.work_inputs``), ``load_m4_inputs`` reports no ``error`` issue, the
planted regions are there (``regions`` non-empty, AMR rows with NA coordinates per the M4-06 rule,
every member gene known), and the whole directory is a pure function of (size, seed, evidence) --
the fixture half of TE7.

``sim50k`` builds in about 0.2 s, so the properties run on it; ``sim1000`` is checked once per
evidence shape and ``sim4500`` only in the benchmark (benchmarks/test_m4_chain_bench.py).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m4_chain as mc

SETTINGS = settings(
    max_examples=10, deadline=None, suppress_health_check=[HealthCheck.too_slow]
)  # fmt: skip
EMITS = [
    frozenset({"q_value", "pep", "is_decoy", "intensity", "protein_ids"}),
    frozenset({"q_value", "pep"}),
    frozenset({"intensity"}),
    frozenset(),
]


def _inputs(work: Path) -> tuple[Any, list[Any]]:
    from giap.io.work_tables import default_paths, load_m4_inputs

    paths = {k: p for k, p in default_paths(work).items() if p.exists()}
    return load_m4_inputs(paths)


@given(seed=st.integers(1, 10_000), emit=st.sampled_from(EMITS))
@SETTINGS
def test_fixture_passes_the_m4_input_contract(
    tmp_path_factory: pytest.TempPathFactory, seed: int, emit: Any
) -> None:
    work = mc.build(
        tmp_path_factory.mktemp("m4fix") / "work",
        seed=seed,
        size="sim50k",
        evidence=mc.evidence_config(emit=emit),
    )
    missing = [rel for rel in mc.CONTRACT_FILES if not (work / rel).is_file()]
    assert not missing, missing
    _m4, issues = _inputs(work)
    errors = [i for i in issues if i.severity == "error"]
    assert not errors, [i.message for i in errors]


@given(seed=st.integers(1, 10_000))
@SETTINGS
def test_fixture_is_a_pure_function_of_size_and_seed(
    tmp_path_factory: pytest.TempPathFactory, seed: int
) -> None:
    root = tmp_path_factory.mktemp("m4pure")
    a = mc.build(root / "a", seed=seed, size="sim50k")
    b = mc.build(root / "b", seed=seed, size="sim50k")
    differ = [rel for rel in mc.CONTRACT_FILES if (a / rel).read_bytes() != (b / rel).read_bytes()]
    assert not differ, differ


@given(seed=st.integers(1, 10_000))
@SETTINGS
def test_planted_regions_are_present_and_well_formed(
    tmp_path_factory: pytest.TempPathFactory, seed: int
) -> None:
    work = mc.build(tmp_path_factory.mktemp("m4reg") / "work", seed=seed, size="sim50k")
    regions = pd.read_parquet(work / "regions.parquet")
    members = pd.read_parquet(work / "region_genes.parquet")
    genes = set(pd.read_parquet(work / "genes.parquet")["locus_tag"].astype(str))
    assert set(regions["region_type"]) >= {"prophage", "amr"}, sorted(set(regions["region_type"]))
    assert set(regions["region_id"]) == set(members["region_id"]), "a region with no member gene"
    assert set(members["locus_tag"].astype(str)) <= genes
    amr = regions[regions["region_type"] == "amr"]
    assert amr["start"].isna().all() and amr["end"].isna().all(), "AMR rows carry coordinates"
    other = regions[regions["region_type"] != "amr"]
    assert other["start"].notna().all() and (other["end"] > other["start"]).all()


@given(seed=st.integers(1, 10_000))
@SETTINGS
def test_the_m3_tsvs_are_written_in_both_layouts(
    tmp_path_factory: pytest.TempPathFactory, seed: int
) -> None:
    """``gene_table.work_inputs`` reads ``mapping/``, ``default_paths`` the flat copies."""
    work = mc.build(tmp_path_factory.mktemp("m4lay") / "work", seed=seed, size="sim50k")
    for name in ("peptide_hits.tsv", "protein_hits.tsv", "unmapped.tsv"):
        assert (work / name).read_bytes() == (work / "mapping" / name).read_bytes(), name


@given(seed=st.integers(1, 10_000))
@SETTINGS
def test_the_id_map_covers_every_gene(tmp_path_factory: pytest.TempPathFactory, seed: int) -> None:
    """C14 (``id_map.tsv``): one ``source_id`` -> ``locus_tag`` row per C2 gene."""
    work = mc.build(tmp_path_factory.mktemp("m4ids") / "work", seed=seed, size="sim50k")
    id_map = pd.read_csv(work / "id_map.tsv", sep="\t", dtype=str)
    genes = pd.read_parquet(work / "genes.parquet")
    assert list(id_map.columns) == ["source_id", "locus_tag"]
    assert set(id_map["locus_tag"]) == set(genes["locus_tag"].astype(str))
    assert id_map["source_id"].is_unique


def test_sim1000_passes_the_contract(tmp_path: Path) -> None:
    work = mc.build(tmp_path / "work", seed=11, size="sim1000")
    m4, issues = _inputs(work)
    errors = [i for i in issues if i.severity == "error"]
    assert not errors, [i.message for i in errors]
    assert len(m4.genes) >= 1000
    assert len(m4.peptide_hits) > 0 and len(m4.peptides) > 0
