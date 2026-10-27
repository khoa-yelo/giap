"""`gippa setup` on the command line (issue M2-08; GIPPA_SPEC.md §5.1, §7.9.1).

Only runs that need no download are exercised: ``--check`` and adoption of already-present
databases. Network access fails the test.
"""

from __future__ import annotations

import json
import socket
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from setup_guard import setup_guard  # noqa: F401  (autouse: no network, no tools; decision 3)

from gippa.cli import cli


@pytest.fixture(autouse=True)
def no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    def deny(*_a: Any, **_k: Any) -> None:
        raise AssertionError("gippa setup touched the network")

    monkeypatch.setattr(socket, "create_connection", deny)
    monkeypatch.setattr(socket.socket, "connect", deny)
    monkeypatch.setenv("SLURM_JOB_ID", "123")


def _bakta_light(db: Path) -> None:
    from gippa.setupdb import SOURCES

    d = db / SOURCES["bakta_light"].path
    d.mkdir(parents=True)
    (d / "version.json").write_text(json.dumps({"major": 6, "minor": 0, "type": "light"}))


def invoke(*args: str) -> Any:
    return CliRunner().invoke(cli, ["setup", *args], catch_exceptions=False)


NOTHING_ELSE = ("--no-genomad", "--no-promotech-models", "--no-antismash")


def test_check_with_nothing_installed_exits_2(tmp_path: Path) -> None:
    r = invoke("--db-dir", str(tmp_path), "--bakta", "light", *NOTHING_ELSE, "--check")
    assert r.exit_code == 2, r.output
    assert "bakta_light" in r.output + (r.stderr if r.stderr_bytes is not None else "")
    assert not (tmp_path / "bakta").exists()


def test_check_with_installed_db_exits_0(tmp_path: Path) -> None:
    _bakta_light(tmp_path)
    r = invoke("--db-dir", str(tmp_path), "--bakta", "light", *NOTHING_ELSE, "--check")
    assert r.exit_code == 0, r.output


def test_adopting_present_db_writes_manifest(tmp_path: Path) -> None:
    _bakta_light(tmp_path)
    r = invoke("--db-dir", str(tmp_path), "--bakta", "light", *NOTHING_ELSE)
    assert r.exit_code == 0, r.output
    data = json.loads((tmp_path / "db_manifest.json").read_text())
    assert data["schema_version"] == "1.0"
    assert [(i["name"], i["status"]) for i in data["items"]] == [("bakta_light", "ok")]


def test_nothing_planned_exits_0(tmp_path: Path) -> None:
    r = invoke("--db-dir", str(tmp_path), "--bakta", "skip", *NOTHING_ELSE)
    assert r.exit_code == 0, r.output


def test_setup_is_no_longer_a_stub(tmp_path: Path) -> None:
    r = invoke("--db-dir", str(tmp_path), "--bakta", "skip", *NOTHING_ELSE)
    assert "not implemented" not in r.output.lower()
