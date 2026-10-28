"""M6-08 examples 1-3 through ``gippa run`` on sim100k (docs/specs/M6-08-operons-step-and-run.md;
GIPPA_SPEC.md §6.4 "byte-identical tables", §9.1).

The inputs are the M2-21 executor-equality inputs
(tests/integration/test_m2_21_local_vs_nextflow.py):
the gippasim SIM100K genome with M3-19 proteomics and every precomputed layer, ``--modules core``
and one seed, with the test DB from tests/helpers/mapping_steps.py.

* Example 1: ``--executor local`` with ``--abundance`` writes ``OUT/tables/pairs.tsv``,
  ``operons.tsv`` and ``operon_flags.tsv`` and exits 0. ``pairs.tsv`` has one row per adjacent
  pair plus the SW rows. Since M4-24 (sim100k under ``mayu``, detection q estimated) pairs get
  ``both_detected_t`` and ``missed_operon_decoy`` tests, and untested (``none``) rows have p
  NA. Publishing ``tables/`` under the local executor relies on
  ``PUBLISHED_FILES`` (PR #145, lead decision 5).
* Example 3: the same run without ``--abundance`` writes ``pairs.tsv``, and ``operon_flags.tsv``
  is header-only.
* Example 2: ``--executor nextflow --profile farmshare_local`` gives byte-identical tables.

All three tables validate against their schemas (example 6).
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from click.testing import CliRunner

pytestmark = [pytest.mark.slow]

SEED = 7
TABLES = ("pairs.tsv", "operons.tsv", "operon_flags.tsv")
SCHEMA_DIR = Path(__file__).resolve().parents[2] / "src" / "gippa" / "schemas" / "tables"


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import gippasim
    from gippasim.proteomics import ProteomicsConfig

    cfg = dataclasses.replace(
        gippasim.SIM100K, proteomics=ProteomicsConfig(), write_precomputed=True
    )
    return gippasim.simulate(cfg, gippasim.SIM100K_SEED, tmp_path_factory.mktemp("sim100k") / "sim")


@pytest.fixture(scope="module")
def db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    from helpers.mapping_steps import write_test_db

    return Path(write_test_db(tmp_path_factory.mktemp("test_db_m6_08")))


def _run(sim: Any, out: Path, db: Path, *extra: str, abundance: bool = True) -> None:
    from gippa.cli import cli

    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides)]
    if abundance:
        assert sim.abundance is not None, "the fixture writes abundance.tsv"
        args += ["--abundance", str(sim.abundance)]
    args += ["--annotation", str(sim.bakta_json), "--config", str(sim.config_yaml)]
    args += ["--outdir", str(out), "--modules", "core", "--seed", str(SEED), *extra]
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIPPA_DB_DIR", str(db))
        mp.setenv("NXF_ANSI_LOG", "false")
        res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == 0, res.output


@pytest.fixture(scope="module")
def local(sim: Any, db: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("LOCAL") / "OUT"
    _run(sim, out, db, "--executor", "local")
    return out


@pytest.fixture(scope="module")
def local_noab(sim: Any, db: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("LOCAL_NOAB") / "OUT"
    _run(sim, out, db, "--executor", "local", abundance=False)
    return out


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _validate(path: Path, stem: str) -> None:
    frictionless = pytest.importorskip("frictionless")
    d = json.loads((SCHEMA_DIR / f"{stem}.schema.json").read_text(encoding="utf-8"))
    d = {k: v for k, v in d.items() if k != "coordinates"}
    res = frictionless.Resource(path=path.name, basepath=str(path.parent),
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    report = res.validate()
    assert report.valid, (stem, report.flatten(["rowNumber", "fieldName", "type", "note"]))


def _chain_counts(out: Path) -> tuple[int, int]:
    """Adjacent pairs and SW triples of the run's own C2 genes and C1 replicons."""
    from gippa.annotate.layers import genes_of, replicons_of
    from gippa.operons.pairs import adjacent_chain, sandwich_pairs

    genes = genes_of(pd.read_parquet(out / "work" / "genes"))
    reps = replicons_of(pd.read_parquet(out / "work" / "replicons"))
    return len(adjacent_chain(genes, reps)), len(sandwich_pairs(genes, reps))


def test_example_1_three_tables_under_tables(local: Path) -> None:
    for name in TABLES:
        assert (local / "tables" / name).is_file(), f"tables/{name} not written"
        _validate(local / "tables" / name, name.removesuffix(".tsv"))


def test_example_1_pairs_rows_and_tests(local: Path) -> None:
    df = _read(local / "tables" / "pairs.tsv")
    assert "which" not in df.columns
    n_adj, n_sw = _chain_counts(local)
    assert len(df) == n_adj + n_sw
    assert int((df["class"] == "SW").sum()) == n_sw
    tests = {"both_detected_t", "one_detected_pmiss", "missed_operon_decoy", "none"}
    assert set(df["test"]) <= tests
    assert set(df["test"]) - {"none"}, "no pair tested (M4-24: detection q is estimated)"
    assert set(df.loc[df["test"] == "none", "p_value"]) <= {"NA"}


def test_m6_10_sensitivity_has_operon_flags_rows(local: Path) -> None:
    """M6-10 acceptance: a real ``gippa run`` writes the ``operon_flags`` family into
    ``tables/sensitivity.tsv``; sim100k calls no flag, so n is 0 and J is NA (review #243 B1)."""
    df = _read(local / "tables" / "sensitivity.tsv")
    rows = df[df["family"] == "operon_flags"]
    assert len(rows) > 0, sorted(set(df["family"]))
    assert set(rows["n_primary"]) == {"0"} and set(rows["jaccard"]) == {"NA"}


def test_example_3_no_abundance(local_noab: Path) -> None:
    pairs = local_noab / "tables" / "pairs.tsv"
    assert pairs.is_file()
    assert len(_read(pairs)) > 0
    flags = local_noab / "tables" / "operon_flags.tsv"
    from gippa.operons.flags import FLAG_COLUMNS

    assert flags.read_text(encoding="utf-8") == "\t".join(FLAG_COLUMNS) + "\n"
    _validate(flags, "operon_flags")


@pytest.mark.needs_nextflow
def test_example_2_nextflow_gives_byte_identical_tables(
    sim: Any, db: Path, local: Path, tmp_path_factory: pytest.TempPathFactory
) -> None:
    nf = tmp_path_factory.mktemp("NF") / "OUT"
    _run(sim, nf, db, "--executor", "nextflow", "--profile", "farmshare_local")
    for name in TABLES:
        got = nf / "tables" / name
        assert got.is_file(), f"Nextflow run did not publish tables/{name}"
        assert got.read_bytes() == (local / "tables" / name).read_bytes(), name
