"""``gippa run --executor local`` through the mapping steps (issue M3-17 examples 2, 3, 5;
GIPPA_SPEC.md §6.3, §6.4, §7.2.9).

Inputs: a 50-kb gippasim genome with M3-19 proteomics and every precomputed layer
(``config.yaml``, §6.4), ``--modules core``. Ambiguities resolved here:

* Example 2: the seven §7.2.9 files exist in one directory under ``OUT`` (``OUT/mapping`` per
  §9.1 is expected; any directory holding all seven is accepted); ``protein_hits.tsv`` is exactly
  its header line (``gippa.mapping.outputs.PROTEIN_HITS_COLUMNS``).
* Example 3: with ``--proteins`` the local executor cannot run miniprot (no precomputed protein
  layer exists, §6.4): exit 4 (§6.3 "gippa run exits 4 naming the failed step") and a message
  naming ``miniprot``; nothing is written to the mapping tables.
* Example 5: ``--executor nextflow --profile farmshare_local`` on the same inputs and seed gives
  byte-identical ``mapping/*.tsv`` (needs Nextflow; slow).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.mapping_steps import MAPPING_TSVS, MODULE2_FILES, module2_files, write_test_db

from gippa.mapping.outputs import PROTEIN_HITS_COLUMNS

SEED = 7


@pytest.fixture(autouse=True)
def _test_db(tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch) -> None:
    """``gippa run`` runs preflight: a test ``$GIPPA_DB_DIR``, never the real one (M3-17 rev. 12)."""
    monkeypatch.setenv("GIPPA_DB_DIR", str(write_test_db(tmp_path_factory.mktemp("test_db"))))


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    from gippasim.proteomics import ProteomicsConfig
    from gippasim.simulate import SimConfig, simulate

    out = tmp_path_factory.mktemp("sim50k")
    cfg = SimConfig(proteomics=ProteomicsConfig(), write_precomputed=True)
    return simulate(cfg, SEED, out / "sim")


def _args(sim: Any, out: Path, *extra: str) -> list[str]:
    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides)]
    args += ["--annotation", str(sim.bakta_json), "--config", str(sim.config_yaml)]
    args += ["--outdir", str(out), "--modules", "core", "--seed", str(SEED), *extra]
    return args


def _run(args: list[str]) -> Any:
    from gippa.cli import cli

    return CliRunner().invoke(cli, args)


def test_example_2_peptides_only_local_run_writes_all_seven_outputs(
    sim: Any, tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    r = _run(_args(sim, out, "--executor", "local"))
    assert r.exit_code == 0, r.output
    files = module2_files(out)
    assert set(files) == set(MODULE2_FILES)
    lines = files["protein_hits.tsv"].read_text().splitlines()
    assert lines == ["\t".join(PROTEIN_HITS_COLUMNS)]
    assert len(files["peptide_hits.tsv"].read_text().splitlines()) > 1


def test_example_3_proteins_without_miniprot_on_local_exit_4(sim: Any, tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    r = _run(_args(sim, out, "--proteins", str(sim.proteins_faa), "--executor", "local"))
    assert r.exit_code == 4, r.output
    assert "miniprot" in r.output
    assert not list(out.rglob("peptide_hits.tsv"))


def test_local_rerun_is_byte_identical(sim: Any, tmp_path: Path) -> None:
    a, b = tmp_path / "A", tmp_path / "B"
    for out in (a, b):
        r = _run(_args(sim, out, "--executor", "local"))
        assert r.exit_code == 0, r.output
    fa, fb = module2_files(a), module2_files(b)
    for f in MAPPING_TSVS:
        assert fa[f].read_bytes() == fb[f].read_bytes(), f


@pytest.mark.xfail(strict=False, reason="needs M2-20 context_layers")
@pytest.mark.needs_nextflow
@pytest.mark.slow
def test_example_5_local_and_nextflow_mapping_tsvs_are_byte_identical(
    sim: Any, tmp_path: Path
) -> None:
    local, nf = tmp_path / "LOCAL", tmp_path / "NF"
    r = _run(_args(sim, local, "--executor", "local"))
    assert r.exit_code == 0, r.output
    r = _run(_args(sim, nf, "--executor", "nextflow", "--profile", "farmshare_local"))
    assert r.exit_code == 0, r.output
    fl = module2_files(local)
    fn = nf / "mapping"
    for f in MAPPING_TSVS:
        assert (fn / f).is_file(), f"Nextflow run did not publish mapping/{f}"
        assert (fn / f).read_bytes() == fl[f].read_bytes(), f
