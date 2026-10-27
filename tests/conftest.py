"""Capability markers, seeds and shared fixtures (GIPPA_SPEC.md §3.6, §11).

A test that needs a capability carries a marker (see pyproject.toml). When the capability is
missing the test is skipped with an explicit reason, unless ``--require-capabilities`` names it,
in which case the test fails instead, so silent skips can never hide failures in CI.
"""

from __future__ import annotations

import os
import shutil
import socket
import subprocess
from collections.abc import Callable
from functools import cache
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get("GIPPA_ROOT", REPO / "local"))
DB = Path(os.environ.get("GIPPA_DB_DIR", ROOT / "db"))
TOOLS = Path(os.environ.get("GIPPA_TOOLS_DIR", ROOT / "tools"))
TOLERANCES = REPO / "tests" / "calibration" / "tolerances.yaml"


def _probe_timeout() -> float:
    """Seconds per capability probe: ``GIPPA_PROBE_TIMEOUT`` (the nightly sets 180), default 30."""
    try:
        return float(os.environ.get("GIPPA_PROBE_TIMEOUT") or 30)  # GIPPA_PROBE_TIMEOUT default 30
    except ValueError:
        return 30.0


def _ok(cmd: list[str], timeout: float | None = None) -> bool:
    timeout = _probe_timeout() if timeout is None else timeout
    if shutil.which(cmd[0]) is None:
        return False
    try:
        return (
            subprocess.run(cmd, capture_output=True, timeout=timeout, check=False).returncode == 0
        )
    except (subprocess.TimeoutExpired, OSError):
        return False


@cache
def _host_reachable(host: str) -> bool:
    try:
        with socket.create_connection((host, 443), timeout=5):
            return True
    except OSError:
        return False


# M3-24 decision 1: the DIA-NN 2.3.1 Academia image Khoa supplied (never copied into the repo).
DIANN_SIF_DEFAULT = "/farmshare/home/classes/bios/270/jordan/diann_2.3.1.sif"


def _have_diann() -> bool:
    """DIA-NN: the Apptainer image of ``$DIANN_SIF`` (M3-24 decision 12, default above) or an
    unpacked 2.7.0 build under ``$GIPPA_TOOLS_DIR`` (GIPPA_SPEC.md §15.2 step 3)."""
    sif = Path(os.environ.get("DIANN_SIF") or DIANN_SIF_DEFAULT)
    return sif.is_file() or any(TOOLS.glob("diann-2.7.0*/**/diann*"))


CHECKS: dict[str, Callable[[], bool]] = {
    "needs_apptainer": lambda: _ok(["apptainer", "--version"]),
    "needs_docker": lambda: _ok(["docker", "info"]),
    "needs_conda": lambda: _ok(["micromamba", "--version"]),
    "needs_nextflow": lambda: _ok(["nextflow", "-version"]) and _ok(["nf-test", "version"]),
    "needs_slurm": lambda: _ok(
        ["sbatch", "--test-only", "--qos=dev", "-p", "normal", "-t", "1", "--wrap", "true"]
    ),
    "needs_gpu": lambda: bool(os.environ.get("CUDA_VISIBLE_DEVICES")),
    "needs_bakta_db": lambda: (
        (DB / "bakta" / "db-light").is_dir() or (DB / "bakta" / "db").is_dir()
    ),
    "needs_bakta_full": lambda: (DB / "bakta" / "db").is_dir(),
    "needs_genomad_db": lambda: (DB / "genomad_db").is_dir(),
    "needs_promotech_models": lambda: any(DB.glob("promotech/**/RF-HOT.model")),
    "needs_diann": _have_diann,
}

_cache: dict[str, bool] = {}


def _have(name: str) -> bool:
    if name not in _cache:
        _cache[name] = CHECKS[name]()
    return _cache[name]


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--require-capabilities",
        default="",
        help="comma-separated capability markers that must be satisfied (missing = failure)",
    )


def _skip_host_only(item: pytest.Item) -> None:
    """``gnu_userland`` / ``x87_longdouble``: FarmShare-only scripts (spec G19) skip elsewhere.

    Imports ``helpers`` only for a marked item: pytester runs of this conftest in a subprocess
    don't have ``tests/`` on ``sys.path``."""
    if not (item.get_closest_marker("gnu_userland") or item.get_closest_marker("x87_longdouble")):
        return
    from helpers.gnu_userland import (
        GNU_USERLAND_REASON,
        X87_LONGDOUBLE_REASON,
        have_gnu_userland,
        have_x87_longdouble,
    )

    if item.get_closest_marker("gnu_userland") and not have_gnu_userland():
        item.add_marker(pytest.mark.skip(reason=GNU_USERLAND_REASON))
    if item.get_closest_marker("x87_longdouble") and not have_x87_longdouble():
        item.add_marker(pytest.mark.skip(reason=X87_LONGDOUBLE_REASON))


@pytest.hookimpl(trylast=True)
def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    # trylast: runs after ``-m``/``-k`` deselection, so capability probes (including the
    # ``needs_net`` socket) only happen for selected items (PR #115 review)
    required = {c.strip() for c in config.getoption("--require-capabilities").split(",") if c}
    unknown = required - set(CHECKS) - {"needs_net"}
    if unknown:
        raise pytest.UsageError(f"unknown capabilities: {sorted(unknown)}")
    for item in items:
        _skip_host_only(item)
        missing: list[str] = []
        for mark in item.iter_markers():
            if mark.name in CHECKS and not _have(mark.name):
                missing.append(mark.name)
            elif mark.name == "needs_net":
                host = str(mark.args[0]) if mark.args else ""
                if not host or not _host_reachable(host):
                    missing.append(f"needs_net({host})")
        if not missing:
            continue
        base = {m.split("(")[0] for m in missing}
        reason = "missing capability: " + ", ".join(missing)
        if base & required:
            continue  # not skipped: pytest_runtest_setup fails it
        item.add_marker(pytest.mark.skip(reason=reason))


def pytest_runtest_setup(item: pytest.Item) -> None:
    required = {c.strip() for c in item.config.getoption("--require-capabilities").split(",") if c}
    for mark in item.iter_markers():
        if mark.name not in required:
            continue
        if mark.name in CHECKS and not _have(mark.name):
            pytest.fail(f"required capability {mark.name} is missing", pytrace=False)
        if mark.name == "needs_net" and not _host_reachable(str(mark.args[0]) if mark.args else ""):
            pytest.fail(f"required capability needs_net({mark.args}) is missing", pytrace=False)


def pytest_terminal_summary(terminalreporter: pytest.TerminalReporter) -> None:
    # Read the reports, not module state, so the summary also works under pytest-xdist
    lines = []
    for rep in terminalreporter.stats.get("skipped", []):
        longrepr = getattr(rep, "longrepr", None)
        reason = longrepr[2] if isinstance(longrepr, tuple) else str(longrepr)
        if "missing capability" in reason:
            lines.append(f"{rep.nodeid}: {reason.removeprefix('Skipped: ')}")
    if lines:
        terminalreporter.section("skipped for missing capabilities (GIPPA_SPEC §3.6)")
        for line in sorted(lines):
            terminalreporter.write_line(line)


@pytest.fixture(scope="session")
def tolerances() -> dict[str, object]:
    """Protected calibration constants and fixed seeds (tests/calibration/tolerances.yaml)."""
    with TOLERANCES.open() as fh:
        data: dict[str, object] = yaml.safe_load(fh)
    return data


@pytest.fixture(scope="session")
def golden_seed(tolerances: dict[str, object]) -> int:
    seeds = tolerances["seeds"]
    assert isinstance(seeds, dict)
    return int(seeds["GOLDEN_SEED"])


@pytest.fixture(autouse=True)
def _no_real_uniop_probe(request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch) -> None:
    """Preflight's UniOP probe (14d) reports success unless a test opts in with
    ``real_uniop_probe``: it would otherwise run the session's real UniOP env or image."""
    if request.node.get_closest_marker("real_uniop_probe"):
        return
    from gippa.preflight import environment

    monkeypatch.setattr(environment, "probe_uniop", lambda *a, **k: (0, ""), raising=False)


# M11-15a item 1 edge case (d): no test sees the real host's container runtimes, Slurm, cluster
# name, `gippa setup` pointer file or $GIPPA_DB_DIR; the probe is pinned to a FarmShare host
# (PINNED_PROBE). Opt out with ``@pytest.mark.real_host_probe``.
_REAL_DB_MARKERS = frozenset({
    "e2e", "needs_bakta_db", "needs_bakta_full", "needs_genomad_db", "needs_promotech_models",
    "needs_diann",
})  # fmt: skip
_REAL_DB_DIRS = ("contract", "e2e")
# The pinned host (lead decision on 15a): FarmShare, so `--profile auto` resolves to
# `farmshare,apptainer`, today's default. Defined in tests/helpers/probe_pin.py.


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers",
        "real_host_probe: keeps the real profile-detection probe, XDG dirs and $GIPPA_DB_DIR "
        "(tests/conftest.py pins them otherwise, M11-15a)",
    )


def _keeps_real_db(item: pytest.Item) -> bool:
    """Real-database tests (capability markers, tests/contract, tests/e2e) keep $GIPPA_DB_DIR."""
    if any(item.get_closest_marker(m) for m in _REAL_DB_MARKERS):
        return True
    try:
        parts = Path(item.path).resolve().relative_to(REPO / "tests").parts
    except ValueError:
        return False
    return bool(parts) and parts[0] in _REAL_DB_DIRS


@pytest.fixture(autouse=True)
def _pinned_host_probe(
    request: pytest.FixtureRequest, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Pin ``gippa.runner.profiles.probe_host`` to the FarmShare ``PINNED_PROBE``, point
    ``XDG_CONFIG_HOME``/``XDG_DATA_HOME`` at (uncreated) ``tmp_path`` subdirectories and unset
    ``GIPPA_DB_DIR``. A test that needs a probe patches ``probe_host`` itself."""
    if request.node.get_closest_marker("real_host_probe"):
        return
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / ".xdg" / "config"))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / ".xdg" / "data"))
    if not _keeps_real_db(request.node):
        monkeypatch.delenv("GIPPA_DB_DIR", raising=False)
    try:
        from gippa.runner import profiles
    except ImportError:  # before M11-15a: nothing to pin
        return
    if hasattr(profiles, "Probe"):
        from helpers.probe_pin import PINNED_PROBE

        pinned = profiles.Probe(**PINNED_PROBE)

        def probe_host(*_a: object, **_k: object) -> object:
            return pinned

        probe_host.gippa_pinned = True  # type: ignore[attr-defined]
        monkeypatch.setattr(profiles, "probe_host", probe_host, raising=False)
