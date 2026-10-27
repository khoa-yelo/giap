"""M11-15b review round 1 (PR #438): B2, B3 and nits 1, 3 and 5 (fast tests).

- **B2.** ``gzip.compress(..., mtime=0)`` writes the header OS byte (offset 9) as ``0x03`` on
  Python 3.12 and ``0xff`` on 3.13, so the bundled layers were not "reproducible byte for byte"
  (PROVENANCE.md). Choice pinned here: every gzip file ``make_demo_layers.py`` writes (the
  ``layers/`` directory) carries OS byte ``0xff``, RFC 1952 "unknown", which is what Python 3.13+
  writes by itself, so forcing it changes nothing on newer interpreters. ``--check`` treats a
  layer whose OS byte differs as out of date. ``sim100k.fna.gz`` and ``sim100k_peptides.tsv.gz``
  are byte-for-byte *copies* of committed files whose sha256 other provenance pins (the recorded
  Bakta outputs name the genome's sha256; docs/tutorials/data/PROVENANCE.md the peptides'), so
  the script never re-gzips them and they are left out of the OS-byte rule.
- **B3.** ``GIPPA_DB_DIR`` is unset for the local demo run (ADR-0066 decision 1) and restored
  afterwards, after a success and after a failure; ``--full`` keeps it.
- **Nit 1.** ``build_manifest`` writes ``demo`` only when true; ``demo_banner`` only for
  ``demo: true``.
- **Nit 3.** A failed ``gippa demo --full`` (example 2.3, missing prerequisites) does not leave the
  outdir it created, so the retry after ``gippa setup`` is not refused by example 2.2.
- **Nit 5.** ``--full`` with only part of the function resources installed leaves ``function``
  out.

No network, no Nextflow, no real ``$GIPPA_DB_DIR``: ``_run_local`` and ``_run_preflight`` are
replaced by recorders, as in tests/unit/test_m11_15b_demo.py.
"""

from __future__ import annotations

import gzip
import inspect
import os
import shutil
import subprocess
import sys
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from click.testing import CliRunner

REPO = Path(__file__).resolve().parents[2]
DEMO = REPO / "src" / "gippa" / "data" / "demo"
LAYERS = DEMO / "layers"
GENOME = REPO / "tests" / "data" / "genomes" / "sim100k.fna.gz"
PEPTIDES = REPO / "docs" / "tutorials" / "data" / "sim100k_peptides.tsv.gz"
SCRIPT = REPO / "scripts" / "docs" / "make_demo_layers.py"
OS_UNKNOWN = 0xFF  # RFC 1952 section 2.3.1: OS = 255, "unknown"
OS_OFFSET = 9
REFUSAL = "exists; use --outdir or remove it"
DEMO_MODULES = ["core", "mobile", "bgc"]
BANNER = "Synthetic demo data (SIM100K, gippasim seed 20261014)"


def _layers() -> list[Path]:
    return sorted(LAYERS.glob("*.gz"))


# ------------------------------------------------------------------ B2: fixed gzip OS byte


def test_there_are_bundled_layers() -> None:
    assert len(_layers()) == 8, [p.name for p in _layers()]


@pytest.mark.parametrize("path", _layers(), ids=lambda p: p.name)
def test_bundled_layer_has_os_byte_unknown(path: Path) -> None:
    head = path.read_bytes()[:10]
    assert head[:2] == b"\x1f\x8b", f"{path.name} is not gzip"
    assert head[OS_OFFSET] == OS_UNKNOWN, (
        f"{path.name}: gzip OS byte 0x{head[OS_OFFSET]:02x}, want 0xff (RFC 1952 'unknown'); "
        "Python 3.12 writes 0x03 and 3.13 writes 0xff, so it must be forced"
    )


def _script(script: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(script), *args], capture_output=True, text=True, timeout=300
    )


@pytest.fixture
def fake_checkout(tmp_path: Path) -> Iterator[Path]:
    """Checkout-shaped tree with copies of the script, the demo data, the genome and peptides
    (as in test_m11_15b_demo.py), so the script's checkout-relative paths stay in ``tmp_path``."""
    root = tmp_path / "checkout"
    (root / "scripts" / "docs").mkdir(parents=True)
    shutil.copy2(SCRIPT, root / "scripts" / "docs" / SCRIPT.name)
    shutil.copytree(DEMO, root / "src" / "gippa" / "data" / "demo")
    (root / "tests" / "data" / "genomes").mkdir(parents=True)
    shutil.copy2(GENOME, root / "tests" / "data" / "genomes" / GENOME.name)
    (root / "tests" / "sim").symlink_to(REPO / "tests" / "sim")
    (root / "docs" / "tutorials" / "data").mkdir(parents=True)
    shutil.copy2(PEPTIDES, root / "docs" / "tutorials" / "data" / PEPTIDES.name)
    yield root


@pytest.mark.parametrize("os_byte", [0x00, 0x03, 0x0B])
def test_check_reports_a_layer_whose_os_byte_differs(fake_checkout: Path, os_byte: int) -> None:
    """Same data, mtime 0, only the OS byte differs: ``--check`` must call it out of date."""
    demo = fake_checkout / "src" / "gippa" / "data" / "demo"
    for layer in sorted((demo / "layers").glob("*.gz")):  # the rest as the fixed writer gives
        data = bytearray(gzip.compress(gzip.decompress(layer.read_bytes()), 9, mtime=0))
        data[OS_OFFSET] = OS_UNKNOWN
        layer.write_bytes(bytes(data))
    target = demo / "layers" / "promoters.tsv.gz"
    data = bytearray(target.read_bytes())
    data[OS_OFFSET] = os_byte
    target.write_bytes(bytes(data))
    assert gzip.decompress(target.read_bytes()) == gzip.decompress(
        (LAYERS / target.name).read_bytes()
    )
    res = _script(fake_checkout / "scripts" / "docs" / SCRIPT.name, "--check")
    assert res.returncode == 1, res.stdout + res.stderr
    assert "promoters" in res.stdout + res.stderr


def test_regenerated_layers_have_os_byte_unknown(fake_checkout: Path) -> None:
    """The writer forces the byte on whatever Python runs it (3.12 would write 0x03)."""
    demo = fake_checkout / "src" / "gippa" / "data" / "demo"
    shutil.rmtree(demo / "layers")
    res = _script(fake_checkout / "scripts" / "docs" / SCRIPT.name)
    assert res.returncode == 0, res.stdout + res.stderr
    written = sorted((demo / "layers").glob("*.gz"))
    assert len(written) == 8
    for path in written:
        assert path.read_bytes()[OS_OFFSET] == OS_UNKNOWN, path.name


# ------------------------------------------------------------------ CLI helpers


@pytest.fixture
def clean_env(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    for name in ("GIPPA_DB_DIR", "GIPPA_CONTAINER_RUNTIME", "NUMBA_CACHE_DIR"):
        monkeypatch.delenv(name, raising=False)
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    return cwd


class Recorder:
    """Replaces a ``gippa.cli`` function; records its bound arguments and the environment."""

    def __init__(self, original: Callable[..., Any], then: Callable[[dict[str, Any]], Any]):
        self.sig = inspect.signature(original)
        self.then = then
        self.calls: list[dict[str, Any]] = []
        self.env: list[dict[str, str]] = []

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        bound = self.sig.bind(*args, **kwargs)
        bound.apply_defaults()
        self.calls.append(dict(bound.arguments))
        self.env.append(dict(os.environ))
        return self.then(dict(bound.arguments))


def _local_success(arguments: dict[str, Any]) -> None:
    report = Path(arguments["outdir"]) / "report" / "report.html"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(f"<html><body>{BANNER}</body></html>\n")
    sys.exit(0)


def _local_failure(_arguments: dict[str, Any]) -> None:
    sys.exit(4)


def _patch_local(
    monkeypatch: pytest.MonkeyPatch, then: Callable[[dict[str, Any]], Any]
) -> Recorder:
    import gippa.cli as gc

    rec = Recorder(gc._run_local, then)
    monkeypatch.setattr(gc, "_run_local", rec)

    def boom(*_a: Any, **_k: Any) -> Any:
        raise AssertionError("gippa demo (without --full) must not launch Nextflow")

    monkeypatch.setattr(gc, "_run_nextflow", boom)
    monkeypatch.setattr(gc, "_run_preflight", boom)
    return rec


MISSING_DBS = "environment.tools: missing: bakta_db, genomad_db; run: gippa setup"


def _preflight_failure(_arguments: dict[str, Any]) -> None:
    import click

    click.echo(MISSING_DBS, err=True)
    sys.exit(2)


@pytest.fixture
def preflight(monkeypatch: pytest.MonkeyPatch) -> Recorder:
    import gippa.cli as gc

    rec = Recorder(gc._run_preflight, _preflight_failure)
    monkeypatch.setattr(gc, "_run_preflight", rec)

    def boom(*_a: Any, **_k: Any) -> Any:
        raise AssertionError("Nextflow launched although preflight failed")

    def no_local(*_a: Any, **_k: Any) -> Any:
        raise AssertionError("--full must go through gippa run's Nextflow path")

    monkeypatch.setattr(gc, "_run_nextflow", boom)
    monkeypatch.setattr(gc, "_run_local", no_local)
    return rec


def _invoke(*args: str) -> Any:
    from gippa.cli import cli

    res = CliRunner().invoke(cli, ["demo", *args])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


@pytest.fixture
def user_db(clean_env: Path, monkeypatch: pytest.MonkeyPatch) -> str:
    """The caller's ``$GIPPA_DB_DIR`` (an existing directory, as after ``gippa setup``)."""
    db = clean_env.parent / "user-db"
    db.mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    return str(db)


# ------------------------------------------------------------------ B3: GIPPA_DB_DIR


def test_local_demo_gets_no_db_dir_although_the_caller_set_one(
    user_db: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    local = _patch_local(monkeypatch, _local_success)
    res = _invoke()
    assert res.exit_code == 0, res.output
    (call,) = local.calls
    assert call["inputs"].get("db_dir") is None, call["inputs"]
    assert call["info"].get("db_dir") is None, call["info"]
    (env,) = local.env
    assert "GIPPA_DB_DIR" not in env, "the local demo run sees the caller's GIPPA_DB_DIR"


def test_callers_db_dir_is_restored_after_a_successful_demo(
    user_db: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    _patch_local(monkeypatch, _local_success)
    res = _invoke()
    assert res.exit_code == 0, res.output
    assert os.environ.get("GIPPA_DB_DIR") == user_db


def test_callers_db_dir_is_restored_after_a_failed_demo(
    user_db: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    local = _patch_local(monkeypatch, _local_failure)
    res = _invoke()
    assert res.exit_code == 4, res.output
    assert "GIPPA_DB_DIR" not in local.env[0]
    assert os.environ.get("GIPPA_DB_DIR") == user_db


def test_full_keeps_the_callers_db_dir(user_db: str, preflight: Recorder) -> None:
    res = _invoke("--full")
    assert res.exit_code == 2, res.output
    (call,) = preflight.calls
    assert call["db_dir"] is not None and str(call["db_dir"]) == user_db
    (env,) = preflight.env
    assert env.get("GIPPA_DB_DIR") == user_db, "--full must keep GIPPA_DB_DIR for gippa run"
    assert os.environ.get("GIPPA_DB_DIR") == user_db


# ------------------------------------------------------------------ nit 1: manifest and banner


def test_build_manifest_omits_demo_for_an_ordinary_run() -> None:
    from helpers.manifest_data import build

    assert "demo" not in build()
    assert "demo" not in build(demo=False)


def test_build_manifest_writes_demo_true_for_a_demo_run() -> None:
    from helpers.manifest_data import build, schema_errors

    m = build(demo=True)
    assert m["demo"] is True
    assert schema_errors(m) == []


def test_demo_key_is_outside_the_provenance_core() -> None:
    """``demo`` does not change the provenance digest (a demo run's tables equal a plain run's)."""
    from helpers.manifest_data import build

    plain, demo = build(), build(demo=True)
    demo.pop("demo")
    assert demo == plain


def test_demo_banner_shown_only_for_demo_true() -> None:
    from gippa.report.content import demo_banner

    assert demo_banner({"demo": True}) == BANNER


@pytest.mark.parametrize(
    "doc", [{}, {"demo": False}, {"demo": "true"}, {"demo": 1}, {"demo": None}, None]
)
def test_demo_banner_absent_otherwise(doc: Any) -> None:
    from gippa.report.content import demo_banner

    assert demo_banner(doc) is None


# ------------------------------------------------------------------ nit 3: failed --full


def test_failed_full_demo_leaves_no_outdir(clean_env: Path, preflight: Recorder) -> None:
    res = _invoke("--full")
    assert res.exit_code == 2, res.output
    assert "run: gippa setup" in res.output
    assert not os.path.lexists(clean_env / "gippa-demo"), (
        "a failed gippa demo --full left gippa-demo/ behind; the retry would be refused (2.2)"
    )


def test_retry_after_a_failed_full_demo_is_not_refused(
    clean_env: Path, preflight: Recorder
) -> None:
    """Example 2.3, then (after ``gippa setup``) ``gippa demo --full`` again: preflight runs
    again instead of example 2.2's refusal."""
    first = _invoke("--full")
    assert first.exit_code == 2, first.output
    second = _invoke("--full")
    assert REFUSAL not in second.output, second.output
    assert second.exit_code == 2, second.output
    assert len(preflight.calls) == 2


def test_failed_full_demo_with_custom_outdir_leaves_no_outdir(
    clean_env: Path, preflight: Recorder
) -> None:
    res = _invoke("--full", "--outdir", "sub/try1")
    assert res.exit_code == 2, res.output
    assert not os.path.lexists(clean_env / "sub" / "try1")


# ------------------------------------------------------------------ nit 5: partial resources


def _function_db(db: Path) -> None:
    from helpers.function_db import write_function_db

    write_function_db(db, pd.DataFrame({"locus_tag": [], "db": [], "id": []}))


def _function_resources() -> tuple[str, ...]:
    from gippa.setupdb.installed import FUNCTION_RESOURCES

    return FUNCTION_RESOURCES


@pytest.mark.parametrize("missing", _function_resources())
def test_full_with_one_function_resource_missing_leaves_function_out(
    clean_env: Path, preflight: Recorder, monkeypatch: pytest.MonkeyPatch, missing: str
) -> None:
    from gippa.setupdb.installed import resource

    db = clean_env.parent / "db"
    _function_db(db)
    path = resource(db, missing)
    assert path is not None, missing
    path.unlink()
    assert resource(db, missing) is None
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    res = _invoke("--full")
    assert res.exit_code == 2, res.output
    (call,) = preflight.calls
    assert sorted(call["cfg"].modules) == sorted(DEMO_MODULES), missing


@pytest.mark.parametrize("kept", _function_resources())
def test_full_with_a_single_function_resource_leaves_function_out(
    clean_env: Path, preflight: Recorder, monkeypatch: pytest.MonkeyPatch, kept: str
) -> None:
    from gippa.setupdb.installed import resource

    db = clean_env.parent / "db"
    _function_db(db)
    for name in _function_resources():
        if name != kept:
            path = resource(db, name)
            assert path is not None
            path.unlink()
    assert resource(db, kept) is not None
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    res = _invoke("--full")
    assert res.exit_code == 2, res.output
    (call,) = preflight.calls
    assert sorted(call["cfg"].modules) == sorted(DEMO_MODULES), kept
