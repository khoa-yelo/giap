"""M11-15c1: ``--db-dir`` optional; ``resolve_db_dir`` and the ``gippa setup`` pointer file.

docs/specs/M11-15-portable-clusters.md item 6 ("``--db-dir`` becomes optional (15c1)"), item 5
(``<data dir>``), examples 6.5 and 6.11, the edge cases of item 6 (a new ``--db-dir`` is
recorded and the old directory is not touched; ``GIPPA_DB_DIR`` beats the pointer file; an empty
or unreadable pointer file is ignored with a warning).

Interface pinned here (readings chosen where the spec leaves room):

- New module ``gippa.setupdb.dbdir`` with ``resolve_db_dir(flag, environ, config_home)``:
  ``flag`` is the ``--db-dir`` value or ``None``; ``environ`` a mapping (``GIPPA_DB_DIR``,
  ``XDG_DATA_HOME``, ``HOME``); ``config_home`` the XDG config home directory (what
  ``$XDG_CONFIG_HOME`` names, else ``~/.config``), so the pointer file is
  ``config_home / "gippa" / "db_dir"``. It returns an object with ``.path`` (a ``Path``) and
  ``.source``, one of ``"flag"``, ``"env"``, ``"pointer"``, ``"default"`` (a NamedTuple or a
  frozen dataclass both fit). Order: flag, ``GIPPA_DB_DIR``, pointer file, ``<data dir>/db``
  with ``<data dir>`` = ``$XDG_DATA_HOME/gippa``, else ``$HOME/.local/share/gippa``. An empty
  ``GIPPA_DB_DIR`` or ``XDG_DATA_HOME`` counts as unset (as ``gippa run`` treats
  ``GIPPA_DB_DIR`` today); the pointer file's one line is stripped of surrounding whitespace.
- ``gippa setup`` writes the absolute DB directory (one line) to
  ``$XDG_CONFIG_HOME/gippa/db_dir``, else ``$HOME/.config/gippa/db_dir``.
- ``gippa run``, ``gippa preflight --profile`` and ``gippa makedb`` resolve through the helper and,
  when the pointer file was used, print ``gippa <cmd>: DB directory <path> (from gippa setup)``.
  ``gippa run`` hands the resolved directory to ``_run_preflight`` (4th argument, ``db_dir``);
  ``gippa preflight --profile`` to ``gippa.preflight.run.run_preflight`` (6th positional
  argument, ``db_dir``), as today.
- An empty or unreadable pointer file: a warning naming the file on stderr, then the default.

Every test runs with ``HOME``, ``XDG_CONFIG_HOME`` and ``XDG_DATA_HOME`` inside ``tmp_path``; the
real ``$GIPPA_DB_DIR`` is unset, no network, no tool (``setup_guard``).
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from setup_guard import setup_guard  # noqa: F401  (autouse: no network, no tools)

from gippa.cli import cli

NOTHING = ("--bakta", "skip", "--no-genomad", "--no-promotech-models", "--no-antismash")
FROM_SETUP = "(from gippa setup)"


class _Stop(Exception):
    """Raised by the stubs once the resolved DB directory has been captured."""


@pytest.fixture(autouse=True)
def _home(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    for var in ("XDG_CONFIG_HOME", "XDG_DATA_HOME", "GIPPA_DB_DIR", "GIPPA_ROOT",
                "NXF_APPTAINER_CACHEDIR", "APPTAINER_CACHEDIR"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    return home


def _dbdir() -> Any:
    try:
        from gippa.setupdb import dbdir
    except ImportError as exc:
        pytest.fail(f"gippa.setupdb.dbdir is not implemented ({exc})")
    return dbdir


def _resolve(flag: Path | None, environ: dict[str, str], config_home: Path) -> Any:
    return _dbdir().resolve_db_dir(flag, environ, config_home)


def _pointer(config_home: Path, text: str) -> Path:
    p = config_home / "gippa" / "db_dir"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")
    return p


def _setup(*args: str) -> Any:
    return CliRunner().invoke(cli, ["setup", *args], catch_exceptions=False)


# ---------------------------------------------------------------- resolve_db_dir (pure)


def test_flag_wins_over_everything(tmp_path: Path) -> None:
    cfg = tmp_path / "cfg"
    _pointer(cfg, str(tmp_path / "recorded") + "\n")
    env = {"GIPPA_DB_DIR": str(tmp_path / "env"), "HOME": str(tmp_path / "home")}
    got = _resolve(tmp_path / "flag", env, cfg)
    assert (got.path, got.source) == (tmp_path / "flag", "flag")


def test_env_wins_over_pointer(tmp_path: Path) -> None:
    """Edge case: ``GIPPA_DB_DIR`` set and the file pointing elsewhere: the variable wins."""
    cfg = tmp_path / "cfg"
    _pointer(cfg, str(tmp_path / "recorded") + "\n")
    env = {"GIPPA_DB_DIR": str(tmp_path / "env"), "HOME": str(tmp_path / "home")}
    got = _resolve(None, env, cfg)
    assert (got.path, got.source) == (tmp_path / "env", "env")


def test_pointer_used_without_flag_or_env(tmp_path: Path) -> None:
    cfg = tmp_path / "cfg"
    _pointer(cfg, f"  {tmp_path / 'recorded'}  \n")
    got = _resolve(None, {"HOME": str(tmp_path / "home")}, cfg)
    assert (got.path, got.source) == (tmp_path / "recorded", "pointer")


def test_default_is_xdg_data_home(tmp_path: Path) -> None:
    env = {"XDG_DATA_HOME": str(tmp_path / "data"), "HOME": str(tmp_path / "home")}
    got = _resolve(None, env, tmp_path / "cfg")
    assert (got.path, got.source) == (tmp_path / "data" / "gippa" / "db", "default")


@pytest.mark.parametrize("xdg", [None, ""])
def test_default_falls_back_to_local_share(tmp_path: Path, xdg: str | None) -> None:
    env = {"HOME": str(tmp_path / "home"), "GIPPA_DB_DIR": ""}
    if xdg is not None:
        env["XDG_DATA_HOME"] = xdg
    got = _resolve(None, env, tmp_path / "cfg")
    assert got.path == tmp_path / "home" / ".local" / "share" / "gippa" / "db"
    assert got.source == "default"


@pytest.mark.parametrize("content", ["", "   \n", "\n\n"])
def test_empty_pointer_is_ignored(tmp_path: Path, content: str) -> None:
    cfg = tmp_path / "cfg"
    _pointer(cfg, content)
    got = _resolve(None, {"HOME": str(tmp_path / "home")}, cfg)
    assert got.source == "default"
    assert got.path == tmp_path / "home" / ".local" / "share" / "gippa" / "db"


def test_unreadable_pointer_is_ignored(tmp_path: Path) -> None:
    cfg = tmp_path / "cfg"
    (cfg / "gippa" / "db_dir").mkdir(parents=True)  # a directory where the file should be
    got = _resolve(None, {"HOME": str(tmp_path / "home")}, cfg)
    assert got.source == "default"


@settings(
    max_examples=60, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(have_flag=st.booleans(), have_env=st.booleans(), have_pointer=st.booleans())
def test_precedence_property(
    tmp_path: Path, have_flag: bool, have_env: bool, have_pointer: bool
) -> None:
    """The first present of flag, GIPPA_DB_DIR, pointer file, default is used, always."""
    case = tmp_path / f"case-{int(have_flag)}{int(have_env)}{int(have_pointer)}"
    cfg = case / "cfg"
    ptr = cfg / "gippa" / "db_dir"
    if have_pointer:
        _pointer(cfg, str(case / "recorded") + "\n")
    elif ptr.exists():
        ptr.unlink()
    env = {"HOME": str(case / "home")}
    if have_env:
        env["GIPPA_DB_DIR"] = str(case / "env")
    got = _resolve(case / "flag" if have_flag else None, env, cfg)
    order = [("flag", have_flag, case / "flag"), ("env", have_env, case / "env"),
             ("pointer", have_pointer, case / "recorded"),
             ("default", True, case / "home" / ".local" / "share" / "gippa" / "db")]  # fmt: skip
    want = next((src, path) for src, present, path in order if present)
    assert (got.source, got.path) == want


# ---------------------------------------------------------------- gippa setup writes the pointer


def test_setup_without_db_dir_is_not_a_usage_error(tmp_path: Path, _home: Path) -> None:
    r = _setup(*NOTHING)
    assert r.exit_code == 0, r.output
    assert "Missing option '--db-dir'" not in r.output


def test_setup_records_default_db_dir_under_dot_config(_home: Path) -> None:
    r = _setup(*NOTHING)
    assert r.exit_code == 0, r.output
    ptr = _home / ".config" / "gippa" / "db_dir"
    assert ptr.is_file(), "setup must record the DB directory in ~/.config/gippa/db_dir"
    want = _home / ".local" / "share" / "gippa" / "db"
    assert ptr.read_text(encoding="utf-8").strip() == str(want)
    assert len(ptr.read_text(encoding="utf-8").strip().splitlines()) == 1


def test_setup_records_given_db_dir_under_xdg_config_home(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    monkeypatch.chdir(tmp_path)
    r = _setup("--db-dir", "rel/db", *NOTHING)
    assert r.exit_code == 0, r.output
    ptr = tmp_path / "cfg" / "gippa" / "db_dir"
    assert ptr.read_text(encoding="utf-8").strip() == str((tmp_path / "rel" / "db").absolute())


def test_setup_uses_gippa_db_dir_and_records_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    db = tmp_path / "envdb"
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    r = _setup(*NOTHING)
    assert r.exit_code == 0, r.output
    assert (tmp_path / "cfg" / "gippa" / "db_dir").read_text().strip() == str(db)


def test_new_db_dir_is_recorded_and_old_left_alone(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Edge case: ``--db-dir`` different from the recorded one."""
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    old, new = tmp_path / "old", tmp_path / "new"
    assert _setup("--db-dir", str(old), *NOTHING).exit_code == 0
    (old / "keep.txt").parent.mkdir(parents=True, exist_ok=True)
    (old / "keep.txt").write_text("mine")
    before = sorted(p.relative_to(old) for p in old.rglob("*"))
    r = _setup("--db-dir", str(new), *NOTHING)
    assert r.exit_code == 0, r.output
    assert (tmp_path / "cfg" / "gippa" / "db_dir").read_text().strip() == str(new)
    assert sorted(p.relative_to(old) for p in old.rglob("*")) == before
    assert (old / "keep.txt").read_text() == "mine"


# ---------------------------------------------------------------- the readers


def _bakta_json(tmp_path: Path) -> Path:
    p = tmp_path / "in" / "bakta.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"sequences": [], "features": []}))
    return p


def test_makedb_uses_recorded_db_dir_and_says_so(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Example 6.11: ``gippa makedb`` with no ``--db-dir`` after setup."""
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    rec = tmp_path / "recorded"
    rec.mkdir()
    _pointer(tmp_path / "cfg", f"{rec}\n")
    r = CliRunner().invoke(
        cli, ["makedb", "--annotation", str(_bakta_json(tmp_path)), "--outdir", str(tmp_path / "o")]
    )
    assert f"gippa makedb: DB directory {rec} {FROM_SETUP}" in r.output
    # no contaminants there: the error names the recorded directory, not "no --db-dir"
    assert str(rec / "contaminants" / "contaminants.fasta") in r.output
    assert "no --db-dir" not in r.output


def test_makedb_without_anything_uses_default_db_dir(tmp_path: Path, _home: Path) -> None:
    r = CliRunner().invoke(
        cli, ["makedb", "--annotation", str(_bakta_json(tmp_path)), "--outdir", str(tmp_path / "o")]
    )
    default = _home / ".local" / "share" / "gippa" / "db"
    assert str(default / "contaminants" / "contaminants.fasta") in r.output, r.output
    assert FROM_SETUP not in r.output


@pytest.mark.parametrize("content", ["", "\n"])
def test_makedb_empty_pointer_warns_and_uses_default(
    tmp_path: Path, _home: Path, monkeypatch: pytest.MonkeyPatch, content: str
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    ptr = _pointer(tmp_path / "cfg", content)
    r = CliRunner().invoke(
        cli, ["makedb", "--annotation", str(_bakta_json(tmp_path)), "--outdir", str(tmp_path / "o")]
    )
    warn = [ln for ln in r.output.splitlines() if "warning" in ln.lower() and str(ptr) in ln]
    assert warn, f"no warning naming {ptr}:\n{r.output}"
    default = _home / ".local" / "share" / "gippa" / "db"
    assert str(default / "contaminants" / "contaminants.fasta") in r.output


def test_makedb_env_beats_pointer(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    _pointer(tmp_path / "cfg", f"{tmp_path / 'recorded'}\n")
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "env"))
    r = CliRunner().invoke(
        cli, ["makedb", "--annotation", str(_bakta_json(tmp_path)), "--outdir", str(tmp_path / "o")]
    )
    assert str(tmp_path / "env" / "contaminants" / "contaminants.fasta") in r.output
    assert FROM_SETUP not in r.output


def _run_inputs(tmp_path: Path) -> list[str]:
    d = tmp_path / "in"
    d.mkdir(parents=True, exist_ok=True)
    g, p = d / "g.fa", d / "p.faa"
    g.write_text(">chr1\n" + "ACGT" * 50 + "\n")
    p.write_text(">P1\nMKV\n")
    return ["--genome", str(g), "--proteins", str(p), "--outdir", str(tmp_path / "OUT")]


@pytest.fixture
def captured_run(monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    import gippa.cli

    seen: dict[str, Any] = {}

    def stub(cfg: Any, *a: Any, **k: Any) -> Any:
        seen["db_dir"] = k["db_dir"] if "db_dir" in k else a[2]
        raise _Stop

    monkeypatch.setattr(gippa.cli, "_run_preflight", stub)
    return seen


def test_run_uses_recorded_db_dir_and_says_so(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, captured_run: dict[str, Any]
) -> None:
    """Example 6.5: ``gippa setup`` then ``gippa run`` with no ``--db-dir``."""
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    rec = tmp_path / "recorded"
    assert _setup("--db-dir", str(rec), *NOTHING).exit_code == 0
    r = CliRunner().invoke(cli, ["run", *_run_inputs(tmp_path)])
    assert isinstance(r.exception, _Stop), r.output
    assert captured_run["db_dir"] is not None, "gippa run did not resolve a DB directory"
    assert Path(captured_run["db_dir"]) == rec
    assert f"gippa run: DB directory {rec} {FROM_SETUP}" in r.output


def test_run_default_db_dir_without_pointer(
    tmp_path: Path, _home: Path, captured_run: dict[str, Any]
) -> None:
    r = CliRunner().invoke(cli, ["run", *_run_inputs(tmp_path)])
    assert isinstance(r.exception, _Stop), r.output
    assert captured_run["db_dir"] is not None, "gippa run did not default the DB directory"
    assert Path(captured_run["db_dir"]) == _home / ".local" / "share" / "gippa" / "db"
    assert FROM_SETUP not in r.output


def test_run_env_beats_pointer(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, captured_run: dict[str, Any]
) -> None:
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    _pointer(tmp_path / "cfg", f"{tmp_path / 'recorded'}\n")
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "env"))
    r = CliRunner().invoke(cli, ["run", *_run_inputs(tmp_path)])
    assert isinstance(r.exception, _Stop), r.output
    assert Path(captured_run["db_dir"]) == tmp_path / "env"
    assert FROM_SETUP not in r.output


def test_preflight_profile_uses_recorded_db_dir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import gippa.preflight.run as pr

    seen: dict[str, Any] = {}

    def stub(*a: Any, **k: Any) -> Any:
        seen["db_dir"] = k["db_dir"] if "db_dir" in k else a[5]
        raise _Stop

    monkeypatch.setattr(pr, "run_preflight", stub)
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    rec = tmp_path / "recorded"
    _pointer(tmp_path / "cfg", f"{rec}\n")
    args = [a for a in _run_inputs(tmp_path) if a not in ("--outdir", str(tmp_path / "OUT"))]
    r = CliRunner().invoke(cli, ["preflight", *args, "--profile", "apptainer"])
    assert isinstance(r.exception, _Stop), r.output
    assert seen["db_dir"] is not None and Path(seen["db_dir"]) == rec
    assert f"gippa preflight: DB directory {rec} {FROM_SETUP}" in r.output


def test_helper_reads_the_same_pointer_setup_writes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Writer and readers agree on the file: the helper finds what setup recorded."""
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    rec = tmp_path / "recorded"
    assert _setup("--db-dir", str(rec), *NOTHING).exit_code == 0
    got = _resolve(None, {"HOME": os.environ["HOME"]}, tmp_path / "cfg")
    assert (got.path, got.source) == (rec, "pointer")
