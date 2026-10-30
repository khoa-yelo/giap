"""#503: ``gippa setup --db-dir <path>`` where the path cannot be created is a clean error.

GitHub issue khoa-yelo/gippa#503 (target 1.0.1). On Sherlock ``--db-dir /scratch/$USER/gippa-db``
(``/scratch/$USER`` missing, ``/scratch`` not writable) crashed in ``db_dir.mkdir`` with an uncaught
``PermissionError``: "internal error (a bug)", exit 5. The plan, with a free-space line for a path
that did not exist, had already been printed.

Interface pinned here (readings chosen where the issue leaves room):

- Exit code 1, ``ExitCode.INVALID_INPUT`` (GIPPA_SPEC.md §5.1 "invalid input"). This is the code
  the CLI already uses for a bad output path: ``gippa demo`` exits ``INVALID_INPUT`` with
  ``gippa demo: cannot create <outdir> (<reason>)`` when ``outdir.mkdir`` raises ``OSError``. A
  ``--db-dir`` is the user's input; exit 2 (``ENVIRONMENT``) means "a tool, database or model is
  missing", which is not the case here.
- stderr has one line ``gippa setup: cannot create --db-dir <path>: <reason>``. For ``EACCES``,
  ``<reason>`` contains ``permission denied`` (the ``strerror``, any case). For other ``OSError``s
  (a file in the way) the line keeps the ``cannot create --db-dir <path>`` prefix; the reason
  wording is not pinned. ``<path>`` is the directory as the user gave it. There is no
  traceback and no "internal error" text.
- The check comes before anything that implies success: no ``plan:`` and no ``free:`` line, no
  download, no pull, and no recorded DB-directory pointer (``$XDG_CONFIG_HOME/gippa/db_dir``). The
  licence lines may still come first; this file does not pin them.
- The free-space line (``free: <x> GiB on <where> (DB directory)``) names the directory that was
  actually measured. A ``--db-dir`` that does not exist yet is measured on its nearest existing
  ancestor, so the line names that ancestor as a whole path component. It may also name the
  DB directory, for example ``free: 9.8 GiB on /a (nearest existing directory of /a/b/db)``. The
  implementer may instead create the DB directory before the plan, so that the directory itself
  is measured. The test tells the two apart by answering "no" at the prompt: if the DB directory
  exists afterwards, the line may name it.

The tests call ``gippa.cli.main`` (the console script): it turns an uncaught exception into exit 5
as in the issue, while ``CliRunner`` would report exit 1 for any exception. No test touches the
network or a real tool. ``setup_guard`` blocks them; the fetcher is a fake that records the URL,
and ``subprocess.run`` answers only the conftest ``java``/``nextflow`` fakes, refusing the quota
commands (best effort, ignored by setup).
"""

from __future__ import annotations

import errno
import os
import re
import shutil
import subprocess
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any, NamedTuple

import pytest
from click.testing import CliRunner
from helpers.prereqs import run_fake_prereq
from helpers.setup_opt_outs import NO_RESOURCE_FLAGS, opt_outs
from setup_guard import setup_guard  # noqa: F401  (autouse: no network, no tools)

from gippa.cli import cli, main
from gippa.exitcodes import ExitCode
from gippa.setupdb import HttpsFetcher, SourceError

GIB = 1024**3
# one HTTPS item (geNOMAD), so a plan is printed and a download would follow; no images
ONE_DB = ("--bakta", "skip", "--genomad", "--no-promotech-models",
          *opt_outs(*NO_RESOURCE_FLAGS), *opt_outs("images", "uniop"))  # fmt: skip
ROOT = hasattr(os, "geteuid") and os.geteuid() == 0


class _Usage(NamedTuple):
    total: int
    used: int
    free: int


@pytest.fixture(autouse=True)
def env(
    setup_guard: list[str],  # noqa: F811  (patched after the guard, so the fakes win)
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> Iterator[dict[str, Any]]:
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USER", "alice")
    monkeypatch.setenv("LOGNAME", "alice")
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    for var in ("XDG_DATA_HOME", "GIPPA_DB_DIR", "GIPPA_ROOT", "NXF_APPTAINER_CACHEDIR",
                "APPTAINER_CACHEDIR", "APPTAINER_TMPDIR", "SINGULARITY_CACHEDIR"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(10_000 * GIB, 0, 9_000 * GIB))
    commands: list[list[str]] = []

    def run(argv: Sequence[Any], *a: Any, **k: Any) -> Any:
        argv = [str(x) for x in argv]
        if (done := run_fake_prereq(argv, *a, **k)) is not None:  # 15d: java -version
            return done
        commands.append(argv)
        if argv[0] in ("quota", "lfs"):
            return subprocess.CompletedProcess(argv, 1, "", "")
        raise AssertionError(f"unexpected subprocess.run {argv}")

    monkeypatch.setattr(subprocess, "run", run)
    fetched: list[str] = []

    def fetch(_self: Any, url: str, dest: Path) -> None:
        fetched.append(url)
        raise SourceError("fake.example", "refused by the test")

    monkeypatch.setattr(HttpsFetcher, "fetch", fetch)
    yield {"fetched": fetched, "commands": commands, "pointer": tmp_path / "cfg" / "gippa"}


def _main(args: Sequence[str], capsys: pytest.CaptureFixture[str]) -> tuple[int, str]:
    """Run the console script; (exit code, stdout + stderr)."""
    with pytest.raises(SystemExit) as info:
        main(list(args))
    out = capsys.readouterr()
    code = info.value.code
    return int(code) if isinstance(code, int) else 1, out.out + out.err


@pytest.fixture
def read_only(tmp_path: Path) -> Iterator[Path]:
    """A directory the test user can list but not write (``0o500``); restored afterwards."""
    if ROOT:
        pytest.skip("root writes into a read-only directory")
    ro = tmp_path / "ro"
    ro.mkdir()
    ro.chmod(0o500)
    try:
        yield ro
    finally:
        ro.chmod(0o700)


def _assert_clean_error(code: int, out: str, db: Path, env: dict[str, Any]) -> str:
    """The shared contract of an uncreatable --db-dir; returns the error line."""
    assert "internal error" not in out, out
    assert "Traceback" not in out and "uncaught exception" not in out, out
    assert code == ExitCode.INVALID_INPUT, f"exit {code}, want 1 (invalid input):\n{out}"
    lines = [ln for ln in out.splitlines() if f"cannot create --db-dir {db}" in ln]
    assert len(lines) == 1, f"no single 'cannot create --db-dir {db}' line:\n{out}"
    assert lines[0].startswith("gippa setup: cannot create --db-dir "), lines[0]
    # nothing that implies success: no plan, no free space, no download, no pointer
    assert not re.search(r"(?m)^gippa setup: (plan|free):", out), out
    assert env["fetched"] == [], env["fetched"]
    assert not (env["pointer"] / "db_dir").exists()
    assert not db.exists()
    return lines[0]


# ---------------------------------------------------------------- an uncreatable --db-dir


def test_db_dir_under_a_read_only_parent_is_a_clean_permission_error(
    read_only: Path, env: dict[str, Any], capsys: pytest.CaptureFixture[str]
) -> None:
    """The issue's case: the missing parent sits in a directory the user cannot write."""
    db = read_only / "khoang99" / "gippa-db"
    code, out = _main(["setup", "--db-dir", str(db), "--yes", *ONE_DB], capsys)
    line = _assert_clean_error(code, out, db, env)
    assert os.strerror(errno.EACCES).lower() in line.lower(), line
    assert "permission denied" in line.lower(), line


def test_db_dir_directly_in_a_read_only_directory(
    read_only: Path, env: dict[str, Any], capsys: pytest.CaptureFixture[str]
) -> None:
    """Only the last component is missing; its parent exists but is read-only."""
    db = read_only / "gippa-db"
    code, out = _main(["setup", "--db-dir", str(db), "--yes", *ONE_DB], capsys)
    line = _assert_clean_error(code, out, db, env)
    assert "permission denied" in line.lower(), line


def test_db_dir_from_gippa_db_dir_under_a_read_only_parent(
    read_only: Path,
    env: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """``$GIPPA_DB_DIR`` resolves to the same directory: the same clean error and exit 1. The
    line names the path; it may say ``GIPPA_DB_DIR`` instead of ``--db-dir``."""
    db = read_only / "missing" / "db"
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    code, out = _main(["setup", "--yes", *ONE_DB], capsys)
    assert "internal error" not in out and "Traceback" not in out, out
    assert code == ExitCode.INVALID_INPUT, f"exit {code}:\n{out}"
    lines = [ln for ln in out.splitlines() if "cannot create" in ln and str(db) in ln]
    assert len(lines) == 1, out
    assert lines[0].startswith("gippa setup: cannot create "), lines[0]
    assert "permission denied" in lines[0].lower(), lines[0]
    assert not re.search(r"(?m)^gippa setup: (plan|free):", out), out
    assert env["fetched"] == []


@pytest.mark.parametrize("where", ["at", "parent"], ids=["file-at-path", "file-as-parent"])
def test_db_dir_blocked_by_a_file_is_a_clean_error(
    tmp_path: Path, env: dict[str, Any], capsys: pytest.CaptureFixture[str], where: str
) -> None:
    """A regular file where the directory (or one of its parents) should be: ``FileExistsError``
    or ``NotADirectoryError`` from ``mkdir``, the same contract, with the OS reason."""
    blocker = tmp_path / "blocker"
    blocker.write_text("not a directory\n")
    db = blocker if where == "at" else blocker / "db"
    code, out = _main(["setup", "--db-dir", str(db), "--yes", *ONE_DB], capsys)
    assert "internal error" not in out and "Traceback" not in out, out
    assert code == ExitCode.INVALID_INPUT, f"exit {code}:\n{out}"
    lines = [ln for ln in out.splitlines() if f"cannot create --db-dir {db}" in ln]
    assert len(lines) == 1, out
    assert lines[0].startswith("gippa setup: cannot create --db-dir "), lines[0]
    reason = lines[0].split(f"cannot create --db-dir {db}", 1)[1]
    assert reason.strip(" :()"), f"no reason given: {lines[0]!r}"
    assert not re.search(r"(?m)^gippa setup: (plan|free):", out), out
    assert env["fetched"] == []
    assert blocker.read_text() == "not a directory\n"  # never replaced or removed


def test_check_with_an_uncreatable_db_dir_is_not_an_internal_error(
    read_only: Path, env: dict[str, Any], capsys: pytest.CaptureFixture[str]
) -> None:
    """``--check`` creates nothing; a missing DB directory is reported, never a crash."""
    db = read_only / "missing" / "db"
    code, out = _main(["setup", "--db-dir", str(db), "--check", *ONE_DB], capsys)
    assert "internal error" not in out and "Traceback" not in out, out
    assert code != ExitCode.INTERNAL, out
    assert env["fetched"] == []


def test_writable_db_dir_still_proceeds_to_the_downloads(
    tmp_path: Path, env: dict[str, Any], capsys: pytest.CaptureFixture[str]
) -> None:
    """Control: a creatable --db-dir (two missing levels) is created and setup goes on."""
    db = tmp_path / "new" / "deeper" / "db"
    code, out = _main(["setup", "--db-dir", str(db), "--yes", *ONE_DB], capsys)
    assert "cannot create" not in out, out
    assert db.is_dir()
    assert env["fetched"], f"setup did not reach the download:\n{out}"
    assert code == ExitCode.ENVIRONMENT, out  # the fake fetcher refuses: exit 2 by the item


# ---------------------------------------------------------------- the free-space line


def _tty(monkeypatch: pytest.MonkeyPatch) -> None:
    from gippa.setupdb import space, uniop_build

    monkeypatch.setattr(uniop_build, "stdin_is_tty", lambda: True)
    monkeypatch.setattr(space, "stdin_is_tty", lambda: True, raising=False)


def _names_path(line: str, path: Path) -> bool:
    """``path`` occurs in ``line`` as a whole path (not the prefix of a longer path)."""
    return re.search(re.escape(str(path)) + r"(?![\w./-])", line) is not None


@pytest.mark.parametrize("missing", [1, 3], ids=["one-missing-level", "three-missing-levels"])
def test_free_line_names_the_measured_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any], missing: int
) -> None:
    """The DB directory does not exist; ``disk_usage`` measures its nearest existing ancestor,
    so the free line names that ancestor (unless the directory was created before the plan)."""
    base = tmp_path / "proj"
    base.mkdir()
    db = base.joinpath(*[f"m{i}" for i in range(missing)])
    measured: list[str] = []

    def usage(p: Any) -> _Usage:
        measured.append(str(p))
        return _Usage(10_000 * GIB, 0, 9_000 * GIB)

    monkeypatch.setattr(shutil, "disk_usage", usage)
    _tty(monkeypatch)
    r = CliRunner().invoke(cli, ["setup", "--db-dir", str(db), *ONE_DB], input="n\n",
                           catch_exceptions=False)  # fmt: skip
    assert r.exit_code == 0, r.output  # "no" at the prompt: nothing downloaded, exit 0
    assert env["fetched"] == []
    free = [ln for ln in r.output.splitlines() if "free:" in ln and "DB directory" in ln]
    assert len(free) == 1, f"no free-space line for the DB directory:\n{r.output}"
    if db.exists():  # created before the plan: the directory itself was measured
        assert _names_path(free[0], db), free[0]
        return
    assert str(base) in measured, measured  # the ancestor is what was measured
    assert _names_path(free[0], base), (
        f"the free line names {db}, which does not exist; {base} was measured:\n{free[0]}"
    )
