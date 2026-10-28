"""M11-14d, ADR-0061 decisions 1 and 4, through ``gippa run`` (both executors).

The ``tu_source_distance_model`` row flags come from ``operon_source()`` reading
``operon_meta.json``, an OPTIONAL input of the ``operons`` and ``context_stats`` steps. The unit
tests (tests/unit/test_m11_14d_fallback_flags.py) feed that file by hand; here the run's own
wiring must deliver it (workflow/main.nf OPERONS and CONTEXT calls, context.nf's filter, and the
local runner's registry plan). If the wiring drops it, the flags vanish silently (fail-open).

Inputs: the gippasim SIM100K genome with M3-19 proteomics, abundance and every precomputed layer
(as tests/integration/test_m6_08_operons_run.py), ``--modules core``, one seed.

* Distance-model run: ``annotation.precomputed.operon_pairs`` removed and
  ``annotation.uniop.enabled: false`` -> ``operon_source = distance_model``. Every
  ``pairs.tsv`` row with ``test = both_detected_t`` and every ``autocorrelation.tsv`` row with
  ``statistic = morans_i`` and ``null_model = tu_block`` carries ``tu_source_distance_model``; no
  other row does.
* Non-fallback run: the precomputed operon pairs (``operon_source = precomputed``; a UniOP-source
  run needs the UniOP tool, which test runs lack; the stub nf-test
  workflow/tests/main_operons_fallback.nf.test pins the run_uniop = true wiring) -> no row
  carries the flag.
* Decision 4 on both: every TU-block ``morans_i`` row carries ``tu_block_tu_error``.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
import yaml
from click.testing import CliRunner

pytestmark = [pytest.mark.slow]

SEED = 7
FLAG = "tu_source_distance_model"
TU_ERROR = "tu_block_tu_error"


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

    return Path(write_test_db(tmp_path_factory.mktemp("test_db_m11_14d")))


@pytest.fixture(scope="module")
def distance_config(sim: Any, tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The sim's config without precomputed operon pairs and with UniOP off."""
    doc = yaml.safe_load(Path(sim.config_yaml).read_text(encoding="utf-8"))
    pre = doc["annotation"]["precomputed"]
    assert "operon_pairs" in pre, "the sim config sets precomputed operon_pairs"
    del pre["operon_pairs"]
    doc["annotation"]["uniop"] = {"enabled": False}
    path = tmp_path_factory.mktemp("cfg_distance") / "config.yaml"
    path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")
    return path


def _run(sim: Any, out: Path, db: Path, config: Path, *extra: str) -> Path:
    from gippa.cli import cli

    assert sim.abundance is not None, "the fixture writes abundance.tsv"
    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides),
            "--abundance", str(sim.abundance), "--annotation", str(sim.bakta_json),
            "--config", str(config), "--outdir", str(out), "--modules", "core",
            "--seed", str(SEED), *extra]  # fmt: skip
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIPPA_DB_DIR", str(db))
        mp.setenv("NXF_ANSI_LOG", "false")
        res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == 0, res.output
    return out


NEXTFLOW = ("--executor", "nextflow", "--profile", "farmshare_local")
LOCAL = ("--executor", "local")


@pytest.fixture(scope="module")
def nf_distance(sim: Any, db: Path, distance_config: Path, tmp_path_factory: Any) -> Path:
    return _run(sim, tmp_path_factory.mktemp("NF_DIST") / "OUT", db, distance_config, *NEXTFLOW)


@pytest.fixture(scope="module")
def nf_precomputed(sim: Any, db: Path, tmp_path_factory: Any) -> Path:
    out = tmp_path_factory.mktemp("NF_PRE") / "OUT"
    return _run(sim, out, db, Path(sim.config_yaml), *NEXTFLOW)


@pytest.fixture(scope="module")
def local_distance(sim: Any, db: Path, distance_config: Path, tmp_path_factory: Any) -> Path:
    return _run(sim, tmp_path_factory.mktemp("LOCAL_DIST") / "OUT", db, distance_config, *LOCAL)


def _read(path: Path) -> pd.DataFrame:
    assert path.is_file(), f"{path} not published"
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _tokens(cell: str) -> set[str]:
    return {t for t in cell.split(";") if t}


def _has(df: pd.DataFrame, flag: str) -> pd.Series:
    return df["flags"].map(lambda c: flag in _tokens(c))


def _targets(pairs: pd.DataFrame, ac: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    on_pairs = pairs["test"] == "both_detected_t"
    on_ac = (ac["statistic"] == "morans_i") & (ac["null_model"] == "tu_block")
    assert on_pairs.any(), "no both_detected_t pair row: the check would be vacuous"
    assert on_ac.any(), "no TU-block morans_i row: the check would be vacuous"
    return on_pairs, on_ac


def _assert_fallback_flags(out: Path) -> None:
    pairs = _read(out / "tables" / "pairs.tsv")
    ac = _read(out / "tables" / "autocorrelation.tsv")
    on_pairs, on_ac = _targets(pairs, ac)
    got_p, got_a = _has(pairs, FLAG), _has(ac, FLAG)
    assert got_p[on_pairs].all(), (
        f"{int((~got_p[on_pairs]).sum())} of {int(on_pairs.sum())} both_detected_t rows lack "
        f"{FLAG}: operon_meta did not reach operons (ADR-0061 decision 1)"
    )
    assert not got_p[~on_pairs].any(), "the flag is on a non-both_detected_t pair row"
    assert got_a[on_ac].all(), (
        f"TU-block morans_i rows lack {FLAG}: operon_meta did not reach context_stats"
    )
    assert not got_a[~on_ac].any(), "the flag is on a row other than TU-block morans_i"
    assert _has(ac, TU_ERROR)[on_ac].all()


def test_nextflow_distance_model_run_flags_the_adr_0061_rows(nf_distance: Path) -> None:
    _assert_fallback_flags(nf_distance)


def test_local_distance_model_run_flags_the_adr_0061_rows(local_distance: Path) -> None:
    _assert_fallback_flags(local_distance)


def test_nextflow_precomputed_run_never_flags_the_distance_model(nf_precomputed: Path) -> None:
    pairs = _read(nf_precomputed / "tables" / "pairs.tsv")
    ac = _read(nf_precomputed / "tables" / "autocorrelation.tsv")
    _, on_ac = _targets(pairs, ac)
    assert not _has(pairs, FLAG).any()
    assert not _has(ac, FLAG).any()
    assert _has(ac, TU_ERROR)[on_ac].all(), "decision 4: every TU-block morans_i row, every run"
