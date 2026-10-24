"""The GNU-userland probe that skips FarmShare script tests elsewhere (spec G19; macOS CI).

``helpers.gnu_userland`` probes ``sha256sum``, ``flock``, ``realpath -e`` and ``stat -c`` once per
session; marked tests skip at collection, and ``Sandbox.run`` / ``M10Sandbox.popen`` skip their
caller, so every test that runs a script through the shared helpers is covered.
"""

from __future__ import annotations

import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest
from helpers import gnu_userland
from helpers.gnu_userland import (
    GNU_USERLAND_REASON,
    have_gnu_userland,
    probe_gnu_userland,
    require_gnu_userland,
)
from helpers.slurm_sandbox import make_sandbox

REPO = Path(__file__).resolve().parents[2]
ALL = ("sha256sum", "flock", "realpath -e", "stat -c")


def _stub(d: Path, name: str, code: int) -> None:
    p = d / name
    p.write_text(f"#!/bin/sh\nexit {code}\n")
    p.chmod(p.stat().st_mode | stat.S_IXUSR)


def _stubs(tmp: Path, codes: dict[str, int]) -> str:
    d = tmp / "bin"
    d.mkdir()
    for name, code in codes.items():
        _stub(d, name, code)
    return str(d)


def test_reason_text() -> None:
    assert (
        GNU_USERLAND_REASON == "FarmShare/GNU userland only (sha256sum, flock, GNU realpath/stat)"
    )


def test_empty_path_misses_everything(tmp_path: Path) -> None:
    assert probe_gnu_userland(_stubs(tmp_path, {})) == ALL


def test_working_tools_are_present(tmp_path: Path) -> None:
    path = _stubs(tmp_path, {"sha256sum": 0, "flock": 0, "realpath": 0, "stat": 0})
    assert probe_gnu_userland(path) == ()


def test_bsd_realpath_and_stat_are_missing(tmp_path: Path) -> None:
    # macOS: BSD realpath/stat exist but reject -e / -c ("illegal option", exit 1)
    path = _stubs(tmp_path, {"sha256sum": 0, "flock": 0, "realpath": 1, "stat": 1})
    assert probe_gnu_userland(path) == ("realpath -e", "stat -c")


def test_failing_stub_counts_as_missing(tmp_path: Path) -> None:
    path = _stubs(tmp_path, {"sha256sum": 127, "flock": 0, "realpath": 0, "stat": 0})
    assert probe_gnu_userland(path) == ("sha256sum",)


@pytest.mark.skipif(sys.platform != "linux", reason="FarmShare (GNU/Linux) expectation")
def test_farmshare_host_has_the_gnu_userland() -> None:
    # on Linux CI and FarmShare the script tests must run, never skip silently
    assert gnu_userland.missing_gnu_userland() == ()
    assert have_gnu_userland()


def test_probe_is_cached() -> None:
    assert gnu_userland.missing_gnu_userland() is gnu_userland.missing_gnu_userland()


def test_require_skips_with_the_reason(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(gnu_userland, "missing_gnu_userland", lambda: ("flock",))
    with pytest.raises(pytest.skip.Exception, match=r"FarmShare/GNU userland only"):
        require_gnu_userland()


def test_require_passes_when_present(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(gnu_userland, "missing_gnu_userland", lambda: ())
    require_gnu_userland()


def test_sandbox_run_skips_without_gnu_userland(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sb = make_sandbox(tmp_path)
    monkeypatch.setattr(gnu_userland, "missing_gnu_userland", lambda: ("sha256sum",))
    with pytest.raises(pytest.skip.Exception, match=r"FarmShare/GNU userland only"):
        sb.run(["true"])
    monkeypatch.setattr(gnu_userland, "missing_gnu_userland", lambda: ())
    assert sb.run(["true"]).returncode == 0


@pytest.mark.gnu_userland
def test_marked_canary() -> None:
    """Runs only where the scripts can run; the next test checks it skips elsewhere."""
    assert have_gnu_userland()


def test_marker_skips_at_collection_when_tools_fail(tmp_path: Path) -> None:
    # failing stubs ahead of the real tools on the host PATH, as on macOS
    stubs = _stubs(tmp_path, {"sha256sum": 1, "flock": 1, "realpath": 1, "stat": 1})
    env = {**os.environ, "PATH": f"{stubs}{os.pathsep}{os.environ.get('PATH', '')}"}
    node = "tests/unit/test_gnu_userland_probe.py::test_marked_canary"
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "-p", "no:cacheprovider", "-p", "no:xdist", "-q",
         "-rs", node],
        cwd=REPO, env=env, capture_output=True, text=True, timeout=300, check=False,
    )  # fmt: skip
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "1 skipped" in proc.stdout, proc.stdout
    assert GNU_USERLAND_REASON in proc.stdout, proc.stdout


def test_x87_probe_reads_the_longdouble_mantissa(monkeypatch: pytest.MonkeyPatch) -> None:
    import numpy as np

    real = np.finfo

    class Float64Info:
        nmant = 52

    def fake(t: object) -> object:
        return Float64Info() if np.dtype(t) == np.dtype(np.longdouble) else real(t)  # type: ignore[arg-type]

    assert gnu_userland.have_x87_longdouble() == (real(np.longdouble).nmant >= 63)
    monkeypatch.setattr(np, "finfo", fake)
    assert not gnu_userland.have_x87_longdouble()
