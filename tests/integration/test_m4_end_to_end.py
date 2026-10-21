"""M4 end to end on giapsim fixtures (docs/specs/M4-20-integration-golden.md examples 1, 3-8).

``build_m4_fixture`` writes the M4-06 contract tables, then ``run_gene_table`` ->
``run_detectability`` -> ``run_context_stats`` -> ``run_sensitivity`` (the spec's call order, the
one MVP-01 wired into ``giap run``) produce ``tables/genes.tsv``, ``detection_summary.tsv``,
``regions_stats.tsv`` and ``sensitivity.tsv``. Every table validates with ``frictionless`` against
its §9.5 resource.

Example 1 runs on ``sim50k`` (the golden fixture, tests/golden/test_m4_golden.py compares its
tables with the committed ones); examples 4-8 on ``sim1000`` at ``stats.permutations = 199`` and
``stats.bootstrap = 200`` (the structure, not the numbers, is under test). Example 3 is TE7
(§11.6): the same seed in two directories gives byte-identical tables.

``genes.tsv`` is written by ``context_stats`` (M5-09 lead decision 3), which
``test_genes_tsv_is_written_by_context_stats`` pins together with its known limitation: a run
whose ``context_stats`` is skipped (no C7/C8) has no ``genes.tsv`` until §6.2's ``finalize``
owns the file.

``sim50k`` is 51 genes, 30 of them with peptide evidence, and the smallest attainable picked gene
q-value is 1/29 = 0.0345: at the default ``stats.alpha.identification`` = 0.01 nothing is ever
detected and every downstream table is vacuous. The fixture runs at ``alpha.identification`` =
0.05 (``GOLDEN_ALPHA``), where 29 of the 51 genes are detected; the vacuous default is pinned in
tests/golden/test_m4_golden.py.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m4_chain as mc

SEED = 20261021
#: the 50-kb genome cannot reach the default 0.01 (see the module docstring)
GOLDEN_ALPHA = 0.05
FAST_STATS = {"permutations": 199, "bootstrap": 200}
#: M4-20 lead decision 14: the writer still sits in ``run_detectability``; M5-09 moves it and
#: removes this mark (strict, so the move cannot land unnoticed)
XFAIL_OWNER = pytest.mark.xfail(
    strict=True, reason="genes.tsv owner moves to context_stats in M5-09 (#138)"
)


@pytest.fixture(scope="module")
def db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return mc.test_db(tmp_path_factory.mktemp("m4_db") / "db")


def _run(
    work: Path, db: Path, *, size: str, seed: int = SEED, cfg: Any | None = None, **build: Any
) -> dict[str, Path]:
    mc.build(work, seed=seed, size=size, **build)
    return mc.run_m4(work, cfg if cfg is not None else _cfg(seed), db)


def _cfg(seed: int = SEED, **over: Any) -> Any:
    stats = {**FAST_STATS, "alpha": {"identification": GOLDEN_ALPHA}, **over.pop("stats", {})}
    return mc.config(seed, stats=stats, **over)


@pytest.fixture(scope="module")
def sim50k(tmp_path_factory: pytest.TempPathFactory, db: Path) -> dict[str, Path]:
    """Example 1: ``sim50k`` with the default evidence (decoys emitted)."""
    work = tmp_path_factory.mktemp("m4_sim50k") / "work"
    return _run(work, db, size="sim50k", cfg=mc.config(SEED, stats={"alpha": {"identification":
                GOLDEN_ALPHA}}))  # fmt: skip


# --------------------------------------------------------------------------- example 1


def test_example_1_writes_the_four_tables(sim50k: dict[str, Path]) -> None:
    missing = [name for name, p in sim50k.items() if not p.is_file()]
    assert not missing, f"M4 chain wrote no {missing} (spec example 1)"
    for name, p in sim50k.items():
        assert p.stat().st_size > 0, name


@pytest.mark.parametrize("name", sorted(mc.OUTPUT_TABLES))
def test_example_1_tables_validate_against_their_resource(
    sim50k: dict[str, Path], name: str
) -> None:
    mc.validate(sim50k[name], mc.OUTPUT_TABLES[name])


def test_example_1_tier_is_decoy_picked_and_detection_is_reached(
    sim50k: dict[str, Path],
) -> None:
    """The fixture must exercise detection: a vacuous golden (no gene detected) is worthless."""
    work = sim50k["genes.tsv"].parent.parent
    g = mc.genes(work)
    assert set(g["detection_tier"].astype(str)) == {"decoy_picked"}
    n = int(g["detected"].astype(bool).sum())
    assert n > 0, "no gene detected: the golden tables would be statistically vacuous"
    assert 0 < n < len(g), (n, len(g))


def test_example_1_every_table_has_rows(sim50k: dict[str, Path]) -> None:
    """Spec "Interface" ("so every M4 output is non-empty"): the synthetic prophage block and AMR
    set give region rows, the detection summary its three sections, sensitivity its variants."""
    rs = mc.read(sim50k["regions_stats.tsv"])
    assert set(rs["region_type"]) >= {"prophage", "amr"}, sorted(set(rs["region_type"]))
    ds = mc.read(sim50k["detection_summary.tsv"])
    assert set(ds["section"]) >= {"tier", "basic", "detectability", "qc"}, sorted(
        set(ds["section"])
    )
    sv = mc.read(sim50k["sensitivity.tsv"])
    assert sorted(set(sv["variant"])) == ["missing=exclude", "missing=min"]
    genes = mc.read(sim50k["genes.tsv"])
    assert len(genes) == len(mc.genes(sim50k["genes.tsv"].parent.parent))


@XFAIL_OWNER
def test_genes_tsv_is_written_by_context_stats(tmp_path: Path, db: Path) -> None:
    """M5-09 lead decision 3: ``context_stats`` owns ``tables/genes.tsv`` -- it is the only step
    holding both the gene table and M5-09's ``island_id`` / ``replichore`` /
    ``replichore_position``, so an earlier write from ``detectability`` would write the file
    twice.

    The known limitation this pins, recorded and not worked around: after ``gene_table`` and
    ``detectability`` there is no ``genes.tsv``, so a run whose ``context_stats`` is skipped (no
    C7/C8) publishes none. §6.2 gives the file to ``finalize``, which is its proper owner.

    Strict xfail (M4-20 lead decision 14): ``run_detectability`` still writes the file here, and
    M5-09 (#138) moves the call. The mark is strict, so the day M5-09 lands this test turns the
    branch red until the mark is removed -- it cannot be forgotten.
    """
    from giap.steps.context_stats import run_context_stats
    from giap.steps.detectability import run_detectability
    from giap.steps.gene_table import run_gene_table

    work = mc.build(tmp_path / "work", seed=SEED, size="sim50k")
    cfg = _cfg()
    tsv = work / "tables" / "genes.tsv"
    with mc.db_env(db):
        run_gene_table(work, cfg)
        run_detectability(work, cfg)
        assert not tsv.exists(), "genes.tsv is written by context_stats, not before it"
        run_context_stats(work, cfg, None)
    assert tsv.is_file(), "context_stats did not write tables/genes.tsv"
    assert len(mc.read(tsv)) == len(mc.genes(work))


# --------------------------------------------------------------------------- example 3 (TE7)


def test_example_3_two_runs_same_seed_are_byte_identical(tmp_path: Path, db: Path) -> None:
    """§11.6 TE7. Two work directories, one seed, byte-identical tables."""
    a = _run(tmp_path / "a" / "work", db, size="sim50k")
    b = _run(tmp_path / "b" / "work", db, size="sim50k")
    da, dbg = mc.digests(a), mc.digests(b)
    assert sorted(da) == sorted(dbg), "the two runs wrote different tables"
    differ = [k for k in da if da[k] != dbg[k]]
    assert not differ, {k: mc.first_diff(a[k], b[k]) for k in differ}


def test_example_3_contract_tables_are_a_pure_function_of_the_seed(tmp_path: Path) -> None:
    """The fixture itself (not only the steps) is reproducible."""
    a = mc.build(tmp_path / "a", seed=SEED, size="sim50k")
    b = mc.build(tmp_path / "b", seed=SEED, size="sim50k")
    for rel in mc.CONTRACT_FILES:
        assert (a / rel).read_bytes() == (b / rel).read_bytes(), rel


# --------------------------------------------------------------------------- examples 4-7


@pytest.fixture(scope="module")
def mayu(tmp_path_factory: pytest.TempPathFactory, db: Path) -> dict[str, Path]:
    """Example 4: no decoys, q-values and PEPs kept."""
    work = tmp_path_factory.mktemp("m4_mayu") / "work"
    return _run(work, db, size="sim1000",
                evidence=mc.evidence_config(emit=frozenset({"q_value", "pep"})))  # fmt: skip


def test_example_4_mayu_tier_and_fdr_2pep_curve_rows(mayu: dict[str, Path]) -> None:
    work = mayu["genes.tsv"].parent.parent
    assert set(mc.genes(work)["detection_tier"].astype(str)) == {"mayu"}
    ds = mc.read(mayu["detection_summary.tsv"])
    tier = ds[ds["section"] == "tier"]
    assert not tier.empty, "detection_summary has no tier section"
    curve = tier[tier["metric"] == "fdr_2pep"]
    assert len(curve) >= 3, "MAYU q_threshold curve rows (metric fdr_2pep) missing"
    assert set(curve["group_type"]) == {"q_threshold"}


@pytest.fixture(scope="module")
def two_peptide(tmp_path_factory: pytest.TempPathFactory, db: Path) -> dict[str, Path]:
    """Example 5: no q-values, no PEPs, no decoys."""
    work = tmp_path_factory.mktemp("m4_2pep") / "work"
    return _run(work, db, size="sim1000", evidence=mc.evidence_config(emit=frozenset()))


def test_example_5_two_peptide_tier_single_peptide_genes_and_variant(
    two_peptide: dict[str, Path],
) -> None:
    work = two_peptide["genes.tsv"].parent.parent
    g = mc.genes(work)
    assert set(g["detection_tier"].astype(str)) == {"two_peptide"}
    status = g["detection_status"].astype(str).value_counts().to_dict()
    assert status.get("single_peptide_unestimated", 0) > 0, status
    sv = mc.read(two_peptide["sensitivity.tsv"])
    assert "single_peptide_detected" in set(sv["variant"]), sorted(set(sv["variant"]))


@pytest.fixture(scope="module")
def protein_list(tmp_path_factory: pytest.TempPathFactory, db: Path) -> dict[str, Path]:
    """Example 6: a protein list with a MaxLFQ abundance table."""
    work = tmp_path_factory.mktemp("m4_protein") / "work"
    ev = mc.evidence_config(proteins_only=True, abundance_table="maxlfq",
                            emit=frozenset({"q_value", "pep", "protein_ids"}))  # fmt: skip
    return _run(work, db, size="sim1000", evidence=ev,
                cfg=_cfg(abundance={"scale": "linear"}))  # fmt: skip


def test_example_6_protein_list_tier_and_relative_scale(protein_list: dict[str, Path]) -> None:
    work = protein_list["genes.tsv"].parent.parent
    g = mc.genes(work)
    assert set(g["detection_tier"].astype(str)) == {"protein_list"}
    scales = set(g["abundance_scale"].dropna().astype(str)) - {"NA", ""}
    assert scales == {"relative"}, scales
    assert (work / "abundance.norm.parquet").is_file(), "C10 not written for abundance_table"


@pytest.fixture(scope="module")
def draft(tmp_path_factory: pytest.TempPathFactory, db: Path) -> dict[str, Path]:
    """Example 7: a 30-gene linear contig and a plasmid."""
    work = tmp_path_factory.mktemp("m4_draft") / "work"
    return _run(work, db, size="sim1000", layout={"linear_contigs": (30,), "plasmids": (60,)})


def test_example_7_partial_genes_are_excluded(draft: dict[str, Path]) -> None:
    g = mc.genes(draft["genes.tsv"].parent.parent)
    partial = g[g["partial"].astype(bool)]
    assert len(partial) >= 2, "a linear contig has a partial gene at each end (§7.3)"
    assert set(partial["detection_status"].astype(str)) == {"partial_excluded"}
    assert not partial["detected"].astype(bool).any()


def test_example_7_plasmid_rows_are_dosage_plus_regulation(draft: dict[str, Path]) -> None:
    rs = mc.read(draft["regions_stats.tsv"])
    plasmid = rs[rs["region_type"] == "plasmid"]
    assert not plasmid.empty, "no plasmid row in regions_stats"
    assert set(plasmid["label"]) == {"dosage_plus_regulation"}
    assert set(plasmid["p_value"]) <= {"NA", ""}, "a plasmid row carries no p-value (M4-18)"


# --------------------------------------------------------------------------- example 8


def test_example_8_contract_violation_stops_the_step(tmp_path: Path, db: Path) -> None:
    """A gene in no TU is a contract ``error``: ``run_gene_table`` raises and writes nothing."""
    from giap.steps.gene_table import run_gene_table
    from giap.steps.registry import StepError

    work = mc.build(tmp_path / "work", seed=SEED, size="sim1000")
    tus = pd.read_parquet(work / "tus.parquet")
    members = list(tus.loc[0, "genes"])
    assert len(members) >= 1
    tus.at[0, "genes"] = members[1:]
    pq.write_table(pa.Table.from_pandas(tus, preserve_index=False), work / "tus.parquet")
    with pytest.raises(StepError) as err, mc.db_env(db):
        run_gene_table(work, _cfg())
    assert "genes in no TU" in str(err.value), str(err.value)
    assert not (work / "genetable" / "genes_pre.parquet").exists()
    assert not (work / "tables").exists()
