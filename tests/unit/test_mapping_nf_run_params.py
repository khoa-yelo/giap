"""``gippa run`` writes MAPPING's inputs into ``OUT/work/params.json`` (M3-17 decision 17).

docs/specs/M3-17-mapping-steps.md, lead decision 17 after the code review of PR #96: the
reviewer's mutants H (``_NORM_INPUTS`` dropped: the raw peptide TSV reaches Nextflow) and I
(``inputs["contaminants"]`` dropped) survived because no test read ``params.json``. ``main.nf``
runs MAPPING from ``params.inputs.peptides`` = the preflight-normalised
``OUT/work/preflight/peptides.norm.parquet``, ``params.inputs.proteins_long`` =
``OUT/work/preflight/proteins.long.faa`` and ``params.inputs.contaminants`` =
``<db_dir>/contaminants/contaminants.fasta`` when the db manifest lists ``contaminants``
(``gippa.runner.local.contaminant_library``).

A fake ``nextflow`` on ``PATH`` exits 0; the database directory is a tmp copy of the real-format
``gippa setup`` manifest (tests/data/db_manifest); the real ``$GIPPA_DB_DIR`` is never read.
"""

from __future__ import annotations

import json
import os
import random
import stat
from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner
from helpers.fake_env import fake_gippa_root, write_setup_manifest

from gippa.exitcodes import ExitCode

FAKE_NF = "#!/bin/sh\nexit 0\n"
LONG_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"


@pytest.fixture(autouse=True)
def _env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    root = fake_gippa_root(tmp_path / "gippa_root")
    monkeypatch.setenv("GIPPA_ROOT", str(root))
    monkeypatch.delenv("GIPPA_DB_DIR", raising=False)
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture
def ins(tmp_path: Path) -> dict[str, Path]:
    rng = random.Random(11)
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    return {
        "genome": _write(tmp_path / "in" / "genome.fa", f">chr1 [topology=circular]\n{seq}\n"),
        "proteins": _write(tmp_path / "in" / "proteins.faa", f">P1\n{LONG_PROTEIN}\n"),
        "peptides": _write(
            tmp_path / "in" / "peptides.tsv",
            "peptide\tprotein_ids\trun\tintensity\n"
            "PEPTIDEK\tP1\tr1\t100\n"
            "AAGLIEVQAPILSR\tP1\tr1\t50\n",
        ),
    }


def _db(path: Path, *, library: bool = True, drop: tuple[str, ...] = ()) -> Path:
    write_setup_manifest(path, drop=drop)
    if library:
        _write(path / "contaminants" / "contaminants.fasta", ">c1\nMKVLAAGIVGRSTEQ\n")
    return path


def _run(ins: dict[str, Path], out: Path, db: Path, *, peptides: bool = True) -> dict[str, str]:
    from gippa.cli import cli

    argv = ["run", "--genome", str(ins["genome"]), "--proteins", str(ins["proteins"])]
    if peptides:
        argv += ["--peptides", str(ins["peptides"])]
    argv += ["--outdir", str(out), "--db-dir", str(db), "--modules", "core"]
    res = CliRunner().invoke(cli, argv)
    log = res.output + (res.stderr if res.stderr_bytes is not None else "")
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == ExitCode.SUCCESS, log
    params = json.loads((out / "work" / "params.json").read_text(encoding="utf-8"))
    inputs = params["inputs"]
    assert isinstance(inputs, dict), params
    return {k: v for k, v in inputs.items() if v}


def test_params_peptides_is_the_normalised_parquet(ins: dict[str, Path], tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    inputs = _run(ins, out, _db(tmp_path / "db"))
    norm = (out / "work" / "preflight" / "peptides.norm.parquet").absolute()
    assert "peptides" in inputs, inputs
    assert Path(inputs["peptides"]) == norm, f"raw peptide table passed to Nextflow: {inputs}"
    assert Path(inputs["peptides"]) != ins["peptides"]
    assert set(pd.read_parquet(inputs["peptides"])["peptide"]) == {"PEPTIDEK", "AAGLIEVQAPILSR"}


def test_params_proteins_long_is_the_preflight_long_fasta(
    ins: dict[str, Path], tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    inputs = _run(ins, out, _db(tmp_path / "db"))
    long = (out / "work" / "preflight" / "proteins.long.faa").absolute()
    assert "proteins_long" in inputs, inputs
    assert Path(inputs["proteins_long"]) == long, inputs
    assert ">P1" in long.read_text(encoding="utf-8")


def test_params_contaminants_from_a_setup_manifest(ins: dict[str, Path], tmp_path: Path) -> None:
    db = _db(tmp_path / "db")
    inputs = _run(ins, tmp_path / "OUT", db)
    assert "contaminants" in inputs, f"no contaminant library in params.json: {inputs}"
    assert (
        Path(inputs["contaminants"]).resolve()
        == (db / "contaminants" / "contaminants.fasta").resolve()
    )


@pytest.mark.parametrize(
    ("library", "drop"), [(False, ()), (True, ("contaminants",))], ids=["no_file", "not_listed"]
)
def test_params_no_contaminants_without_a_listed_library(
    ins: dict[str, Path], tmp_path: Path, library: bool, drop: tuple[str, ...]
) -> None:
    inputs = _run(ins, tmp_path / "OUT", _db(tmp_path / "db", library=library, drop=drop))
    assert "contaminants" not in inputs, inputs


def test_params_without_peptides_has_no_peptides_key(ins: dict[str, Path], tmp_path: Path) -> None:
    """main.nf runs MAPPING only when ``params.inputs.peptides`` is set."""
    inputs = _run(ins, tmp_path / "OUT", _db(tmp_path / "db"), peptides=False)
    assert "peptides" not in inputs, inputs
