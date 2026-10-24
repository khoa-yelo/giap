"""``scripts/pilot/cascadia_offline.py``: the offline Unimod launcher (M10-08a, deliverable 4).

Spec: docs/specs/M10-08-a3-feasibility.md, deliverable 4 and acceptance criterion 6. Cascadia 0.0.7
downloads ``http://www.unimod.org/xml/unimod_tables.xml`` at the end of every run (pyteomics
``mass.unimod.Unimod()`` with no path); R6 forbids it, so the launcher

(a) refuses a ``UNIMOD_XML`` that is not an existing local file (exit 2, before importing
    ``cascadia``);
(b) sets ``pyteomics.mass.unimod._unimod_xml_download_url`` to its absolute path;
(c) hands the rest of ``argv`` to ``cascadia.cascadia.main()`` unchanged.

Stub ``cascadia`` and ``pyteomics`` packages (``helpers.a3_kit.write_stub_packages``) record what
the launcher did; nothing real is imported. The production call is ``python -I`` inside the
``giap-cascadia`` env, where ``-I`` ignores ``PYTHONPATH``: most tests use ``PYTHONPATH`` without
``-I`` and one test runs ``-I`` against a throw-away venv (no pip) holding the stubs.
No statistical claim here, so no calibration test (spec §8.4).
"""

from __future__ import annotations

import ast
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
from helpers.a3_kit import (
    DEFAULT_UNIMOD_URL,
    LAUNCHER,
    STUB_UNIMOD_XML,
    make_stub_venv,
    require,
    write_stub_packages,
)

ARGS = ["sequence", "/data/mzml/WT_1.mzML", "/db/cascadia/cascadia_0.0.7.ckpt", "-o",
        "/root/pilot/A3_cascadia/WT_1/WT_1", "-b", "32", "-t", "0.8"]  # fmt: skip


@pytest.fixture
def stubs(tmp_path: Path) -> Path:
    return write_stub_packages(tmp_path / "site")


@pytest.fixture
def unimod(tmp_path: Path) -> Path:
    p = tmp_path / "db" / "unimod_tables.xml"
    p.parent.mkdir(parents=True)
    p.write_bytes(STUB_UNIMOD_XML)
    return p


def _launch(
    tmp: Path,
    site: Path | None,
    *argv: str,
    python: str = sys.executable,
    isolated: bool = False,
    rc: int = 0,
    cwd: Path | None = None,
) -> tuple[subprocess.CompletedProcess[str], dict[str, object] | None, bool]:
    require(LAUNCHER, "deliverable 4")
    log, imported = tmp / "stub_main.json", tmp / "stub_import.log"
    for p in (log, imported):
        p.unlink(missing_ok=True)
    env = {"PATH": "/usr/bin:/bin", "LANG": "C.UTF-8", "STUB_LOG": str(log),
           "STUB_IMPORT_LOG": str(imported), "STUB_RC": str(rc)}  # fmt: skip
    if site is not None:
        env["PYTHONPATH"] = str(site)
    cmd = [python, *(["-I"] if isolated else []), str(LAUNCHER), *argv]
    res = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=60,
                         cwd=str(cwd or tmp), check=False)  # fmt: skip
    rec = json.loads(log.read_text()) if log.is_file() else None
    return res, rec, imported.is_file()


def _same_path(a: object, b: Path) -> bool:
    return isinstance(a, str) and os.path.isabs(a) and os.path.realpath(a) == os.path.realpath(b)


def test_sets_the_unimod_url_before_cascadia_runs(
    tmp_path: Path, stubs: Path, unimod: Path
) -> None:
    res, rec, imported = _launch(tmp_path, stubs, str(unimod), *ARGS)
    assert res.returncode == 0, res.stdout + res.stderr
    assert rec is not None, "cascadia.cascadia.main() was not called"
    assert imported
    assert rec["url_in_main"] != DEFAULT_UNIMOD_URL
    assert _same_path(rec["url_in_main"], unimod), rec
    # set before cascadia is even imported (deliverable 4: (b) before (c))
    assert _same_path(rec["url_at_import"], unimod), rec


def test_argv_is_passed_through_unchanged(tmp_path: Path, stubs: Path, unimod: Path) -> None:
    odd = [*ARGS, "--path with spaces", "", "-w", "2"]
    res, rec, _ = _launch(tmp_path, stubs, str(unimod), *odd)
    assert res.returncode == 0, res.stdout + res.stderr
    assert rec is not None and rec["argv"] == odd


def test_relative_unimod_path_becomes_absolute(tmp_path: Path, stubs: Path, unimod: Path) -> None:
    rel = os.path.relpath(unimod, tmp_path)
    res, rec, _ = _launch(tmp_path, stubs, rel, *ARGS, cwd=tmp_path)
    assert res.returncode == 0, res.stdout + res.stderr
    assert rec is not None and _same_path(rec["url_in_main"], unimod), rec


def test_cascadia_exit_code_is_propagated(tmp_path: Path, stubs: Path, unimod: Path) -> None:
    res, rec, _ = _launch(tmp_path, stubs, str(unimod), *ARGS, rc=5)
    assert rec is not None
    assert res.returncode == 5


@pytest.mark.parametrize(
    "bad",
    [
        "http://www.unimod.org/xml/unimod_tables.xml",
        "https://www.unimod.org/xml/unimod_tables.xml",
        "file:///etc/unimod_tables.xml",
        "missing",
        "dir",
    ],
)
def test_refuses_a_non_local_or_missing_table_without_importing_cascadia(
    tmp_path: Path, stubs: Path, bad: str
) -> None:
    if bad == "missing":
        bad = str(tmp_path / "absent" / "unimod_tables.xml")
    elif bad == "dir":
        (tmp_path / "a_dir").mkdir()
        bad = str(tmp_path / "a_dir")
    res, rec, imported = _launch(tmp_path, stubs, bad, *ARGS)
    assert res.returncode == 2, res.stdout + res.stderr
    assert rec is None, "cascadia main() ran"
    assert not imported, "cascadia was imported before the table was checked"


def test_no_arguments_exits_2(tmp_path: Path, stubs: Path) -> None:
    res, rec, imported = _launch(tmp_path, stubs)
    assert res.returncode == 2, res.stdout + res.stderr
    assert rec is None and not imported


def test_runs_isolated_inside_an_env_that_has_cascadia(tmp_path: Path, unimod: Path) -> None:
    # The production call: `micromamba run -n giap-cascadia python -I cascadia_offline.py …`.
    # No giap on the path, no PYTHONPATH: only the env's own site-packages.
    env = make_stub_venv(tmp_path / "giap-cascadia")
    res, rec, _ = _launch(tmp_path, None, str(unimod), *ARGS, python=str(env / "bin" / "python"),
                          isolated=True)  # fmt: skip
    assert res.returncode == 0, res.stdout + res.stderr
    assert rec is not None and rec["argv"] == ARGS
    assert _same_path(rec["url_in_main"], unimod)


def test_launcher_imports_no_giap_and_no_network_modules() -> None:
    tree = ast.parse(require(LAUNCHER, "deliverable 4").read_text())
    mods = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module:
            mods.add(node.module.split(".")[0])
    assert "giap" not in mods  # it runs inside giap-cascadia, where giap is not installed
    assert not mods & {"http", "requests", "socket", "ftplib", "urllib3"}
