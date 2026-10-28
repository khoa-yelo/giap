"""``--executor local`` and ``--executor nextflow`` give byte-identical tables (M2-21 example 4;
GIPPA_SPEC.md §6.4 "an integration test checks that they produce byte-identical tables on the
synthetic genome", §11.2 "local and Nextflow executors give identical tables").

Example 4: ``test_local_vs_nextflow`` → ``cmp`` equal for ``operon_pairs.tsv``, ``tus.tsv``,
``regions.tsv`` and every BED.

Inputs: the gippasim SIM100K genome with M3-19 proteomics and every precomputed layer
(``config.yaml``, §6.4), peptides only (the local executor cannot run miniprot), ``--modules
core,mobile,bgc``, one seed. The same ``gippa run`` arguments run once with ``--executor local`` and
once with ``--executor nextflow --profile farmshare_local`` (real Nextflow, prebuilt GIPPA env under
the real ``$GIPPA_ROOT``; no annotation tool runs because every layer is precomputed). The database
directory is the test DB of tests/helpers/mapping_steps.py, never the real one.

Compared byte for byte: every ``*.tsv`` and ``*.bed`` under ``OUT/annotation/`` (the same file
set in both runs; it includes the three named tables) and every ``*.bed`` anywhere else under
``OUT`` outside ``work/`` (e.g. ``mapping/mappings.bed``).
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner

pytestmark = [pytest.mark.needs_nextflow, pytest.mark.slow]

SEED = 7
NAMED = ("operon_pairs.tsv", "tus.tsv", "regions.tsv")


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import gippasim
    from gippasim.proteomics import ProteomicsConfig

    cfg = dataclasses.replace(
        gippasim.SIM100K, proteomics=ProteomicsConfig(), write_precomputed=True
    )
    return gippasim.simulate(cfg, gippasim.SIM100K_SEED, tmp_path_factory.mktemp("sim100k") / "sim")


def _run(sim: Any, out: Path, *extra: str) -> None:
    from gippa.cli import cli

    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides)]
    args += ["--annotation", str(sim.bakta_json), "--config", str(sim.config_yaml)]
    args += ["--outdir", str(out), "--modules", "core,mobile,bgc", "--seed", str(SEED), *extra]
    res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == 0, res.output


def _tables(out: Path) -> dict[str, bytes]:
    ann = {f"annotation/{p.name}": p.read_bytes() for p in sorted((out / "annotation").glob("*"))
           if p.suffix in (".tsv", ".bed")}  # fmt: skip
    beds = {str(p.relative_to(out)): p.read_bytes() for p in sorted(out.rglob("*.bed"))
            if "work" not in p.relative_to(out).parts}  # fmt: skip
    return {**ann, **beds}


@pytest.fixture(scope="module")
def both(
    sim: Any, tmp_path_factory: pytest.TempPathFactory
) -> tuple[dict[str, bytes], dict[str, bytes]]:
    mp = pytest.MonkeyPatch()
    try:
        from helpers.mapping_steps import write_test_db

        mp.setenv("GIPPA_DB_DIR", str(write_test_db(tmp_path_factory.mktemp("test_db_mod"))))
        mp.setenv("NXF_ANSI_LOG", "false")
        local = tmp_path_factory.mktemp("LOCAL") / "OUT"
        nf = tmp_path_factory.mktemp("NF") / "OUT"
        _run(sim, local, "--executor", "local")
        _run(sim, nf, "--executor", "nextflow", "--profile", "farmshare_local")
    finally:
        mp.undo()
    return _tables(local), _tables(nf)


def test_example_4_same_table_set(both: tuple[dict[str, bytes], dict[str, bytes]]) -> None:
    local, nf = both
    for name in NAMED:
        assert f"annotation/{name}" in local, f"local run did not publish annotation/{name}"
    assert any(k.endswith(".bed") for k in local)
    assert sorted(nf) == sorted(local)


@pytest.mark.parametrize("name", NAMED)
def test_example_4_named_tables_byte_identical(
    both: tuple[dict[str, bytes], dict[str, bytes]], name: str
) -> None:
    local, nf = both
    key = f"annotation/{name}"
    assert key in nf, f"Nextflow run did not publish {key}"
    assert nf[key] == local[key], key


def test_example_4_every_tsv_and_bed_byte_identical(
    both: tuple[dict[str, bytes], dict[str, bytes]],
) -> None:
    local, nf = both
    differ = sorted(k for k in local if nf.get(k) != local[k])
    assert not differ, f"local and Nextflow differ: {differ}"
