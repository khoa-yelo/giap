"""M4 end to end on the real M2/M3 outputs (docs/specs/M4-21-m2-m3-conformance.md work item 3,
example 4; GIAP_SPEC.md §6.4, §9.5, §11.6 TE7).

``giap run --executor local`` on giapsim SIM100K with its precomputed layers produces the M2/M3
tables (tests/helpers/m4_real.py); the M4 chain of M4-20 (``m4_chain.run_m4``: ``run_gene_table``
-> ``run_detectability`` -> ``run_context_stats`` -> ``run_sensitivity``) then runs on that very
``OUT/work`` directory -- the layout ``giap run`` writes, not the flat fixture layout -- and its
four ``tables/*.tsv`` validate with ``frictionless``. TE7: the whole real path twice, with the same
seed in two directories, gives byte-identical tables (one machine, so ``regions_stats.tsv`` is
compared byte for byte too; M4-20 decision 9 only concerns cross-machine goldens).

SIM100K has 87 genes and a picked-q floor of 1/85 = 0.0118 > 0.01 (MVP-01 run 1861972), so the
chain runs at ``alpha.identification`` = 0.05, as the M4-20 fixture tests do, with the fast
resampling counts (structure, not numbers, is under test).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m4_chain as mc
import m4_real as mr

SEED = 20261014
STATS = {"permutations": 199, "bootstrap": 200, "alpha": {"identification": 0.05}}


def real_m4(root: Path, db: Path) -> dict[str, Path]:
    sim = mr.sim100k(root / "sim")
    work = mr.giap_run(sim, root / "OUT", db, SEED)
    return mc.run_m4(work, mc.config(SEED, stats=STATS), db)


@pytest.fixture(scope="module")
def db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return mc.test_db(tmp_path_factory.mktemp("m4_real_db") / "db")


@pytest.fixture(scope="module")
def run_a(tmp_path_factory: pytest.TempPathFactory, db: Path) -> dict[str, Path]:
    return real_m4(tmp_path_factory.mktemp("real_a"), db)


def test_real_path_writes_the_four_tables(run_a: dict[str, Path]) -> None:
    missing = [n for n, p in run_a.items() if not p.is_file()]
    assert not missing, f"M4 on the real outputs wrote no {missing}"


@pytest.mark.parametrize(("name", "stem"), sorted(mc.OUTPUT_TABLES.items()))
def test_real_path_tables_validate(run_a: dict[str, Path], name: str, stem: str) -> None:
    """§9.5: every M4 table written from real M2/M3 outputs validates against its resource."""
    mc.validate(run_a[name], stem)


def test_real_path_is_not_vacuous(run_a: dict[str, Path]) -> None:
    """The real run must detect genes and test the planted prophage (SIM100K plants a prophage,
    a BGC, a CRISPR array and a plasmid through ``annotation.precomputed.regions``)."""
    genes = mc.read(run_a["genes.tsv"])
    det = int((genes["detection_status"] == "detected").sum())
    assert 0 < det < len(genes), f"{det} detected of {len(genes)}"
    rs = mc.read(run_a["regions_stats.tsv"])
    assert "prophage" in set(rs["region_type"]), sorted(set(rs["region_type"]))
    assert any(v not in ("NA", "") for v in rs["p_value"]), "no region test ran"


def test_example_4_real_path_twice_is_byte_identical(
    run_a: dict[str, Path], tmp_path: Path, db: Path
) -> None:
    """Example 4 (TE7): giapsim -> ``giap run`` -> M4, twice, gives the same bytes."""
    run_b = real_m4(tmp_path, db)
    a, b = mc.digests(run_a), mc.digests(run_b)
    assert set(a) == set(mc.OUTPUT_TABLES), f"run A wrote only {sorted(a)}"
    differ = [n for n in sorted(mc.OUTPUT_TABLES) if a.get(n) != b.get(n)]
    assert not differ, f"tables differ between two runs: {differ}"
