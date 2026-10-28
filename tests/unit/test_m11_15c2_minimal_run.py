"""M11-15c2: ``gippa run`` without ``--modules`` on a DB directory set up with ``--minimal``.

docs/specs/M11-15-portable-clusters.md item 6 ("``--minimal``"), example 6.3, lead decision Q-L25
(auto ``core,function`` with a notice, not an error).

Interface pinned here: the DB directory's ``db_manifest.json`` carries ``{"setup": {"set":
"minimal"}}`` (written by ``gippa setup --minimal``). A run on it (``--db-dir`` or
``$GIPPA_DB_DIR``) with no ``--modules`` uses ``modules: [core, function]`` from preflight on
(the configuration ``_run_preflight`` receives) and prints the notice below; ``--modules`` given
is used as is, without the notice; a DB directory without that record keeps the default modules.

``_run_preflight`` is faked (it records the configuration) and ``nextflow`` is not on ``PATH``,
so the run stops before launching anything.
"""

from __future__ import annotations

import json
import tomllib
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner

REPO = Path(__file__).resolve().parents[2]
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
NOTICE = (
    "modules core,function (DB directory set up with --minimal; mobile and bgc need "
    f"{PKG} setup without --minimal)"
)
ALL = ["core", "mobile", "bgc", "function"]


@pytest.fixture
def seen(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> list[list[str]]:
    import gippa.cli as cli_mod

    got: list[list[str]] = []

    def fake_preflight(cfg: Any, outdir: Path, inputs: Any, *_a: Any, **_k: Any) -> Any:
        got.append(list(cfg.modules))
        return {}, cfg

    monkeypatch.setattr(cli_mod, "_run_preflight", fake_preflight)
    empty = tmp_path / "empty_bin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))  # no nextflow: the run stops before launching
    for var in ("GIPPA_DB_DIR", "GIPPA_ROOT", "SLURM_JOB_ID", "GIPPA_CONTAINER_RUNTIME"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg_config"))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg_data"))
    return got


def db_dir(tmp_path: Path, setup: dict[str, Any] | None) -> Path:
    db = tmp_path / "db"
    db.mkdir()
    doc: dict[str, Any] = {"schema_version": "1.0", "items": []}
    if setup is not None:
        doc["setup"] = setup
    (db / "db_manifest.json").write_text(json.dumps(doc), encoding="utf-8")
    return db


def run(tmp_path: Path, *args: str) -> Any:
    from gippa.cli import cli

    genome = tmp_path / "g.fa"
    genome.write_text(">c1\nACGTACGTACGT\n", encoding="utf-8")
    argv = ["run", "--genome", str(genome), "--outdir", str(tmp_path / "OUT"),
            "--profile", "docker", *args]  # fmt: skip
    return CliRunner().invoke(cli, argv, catch_exceptions=False)


@pytest.mark.parametrize("via", ["flag", "env"])
def test_minimal_db_dir_runs_core_function_with_the_notice(
    tmp_path: Path, seen: list[list[str]], monkeypatch: pytest.MonkeyPatch, via: str
) -> None:
    db = db_dir(tmp_path, {"set": "minimal"})
    if via == "env":
        monkeypatch.setenv("GIPPA_DB_DIR", str(db))
        r = run(tmp_path)
    else:
        r = run(tmp_path, "--db-dir", str(db))
    assert seen == [["core", "function"]], r.output
    assert NOTICE in r.output
    assert r.output.count(NOTICE) == 1


def test_explicit_modules_win_without_a_notice(tmp_path: Path, seen: list[list[str]]) -> None:
    db = db_dir(tmp_path, {"set": "minimal"})
    r = run(tmp_path, "--db-dir", str(db), "--modules", "core,mobile")
    assert seen == [["core", "mobile"]], r.output
    assert "set up with --minimal" not in r.output


@pytest.mark.parametrize("setup", [None, {"set": "default"}, {}])
def test_a_full_db_dir_keeps_the_default_modules(
    tmp_path: Path, seen: list[list[str]], setup: dict[str, Any] | None
) -> None:
    r = run(tmp_path, "--db-dir", str(db_dir(tmp_path, setup)))
    assert seen == [ALL], r.output
    assert "set up with --minimal" not in r.output
