"""Shared fixtures for unit tests.

``registry`` gives each test an empty step registry (M2-02 "temporary registry"): the module-level
``gippa.steps.registry._REGISTRY`` mapping (name -> StepSpec) is swapped for an empty dict, so toy
steps never collide with the real annotation steps or with other tests.

``_fake_prereqs`` (autouse) puts a fake ``java`` and ``nextflow`` first on ``PATH`` (M11-15d;
``tests/helpers/prereqs.py``).
"""

from __future__ import annotations

import os
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

# Imported now, before any test patches subprocess: the helper keeps the real subprocess.run for
# the conftest fakes (helpers.prereqs.run_fake_prereq).
import helpers.prereqs  # noqa: F401
import pytest


@pytest.fixture
def registry(monkeypatch: pytest.MonkeyPatch) -> Iterator[ModuleType]:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    # Import every real step module first: a module imported for the first time while the registry
    # is swapped would register into the temporary dict and, being cached, never into the real one.
    load_step_modules()
    monkeypatch.setattr(reg, "_REGISTRY", {})
    yield reg


@pytest.fixture
def no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    """No socket connection and no DNS lookup in this test (M11-14 AC 9; not autouse).

    ``socket.socket.connect``, ``socket.create_connection`` and ``socket.getaddrinfo`` raise
    ``OSError``, so a real download fails at once instead of reaching a host. Every 14u test
    module sets ``pytestmark = pytest.mark.usefixtures("no_network")``.
    """
    import socket

    def deny(*args: object, **_kw: object) -> None:
        raise OSError(f"no_network: network access is blocked in this test ({args!r})")

    monkeypatch.setattr(socket.socket, "connect", deny)
    monkeypatch.setattr(socket.socket, "connect_ex", deny)
    monkeypatch.setattr(socket, "create_connection", deny)
    monkeypatch.setattr(socket, "getaddrinfo", deny)


@pytest.fixture(scope="session")
def _fake_prereqs_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    from helpers.prereqs import write_prereqs

    return write_prereqs(tmp_path_factory.mktemp("fake_prereqs"))


@pytest.fixture(autouse=True)
def _fake_prereqs(request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch) -> None:
    """A fake ``java`` (21) and ``nextflow`` first on ``PATH``; ``JAVA_HOME``, ``JAVA_CMD`` and
    ``NXF_JAVA_HOME`` unset (M11-15d).

    ``gippa setup`` and ``gippa run`` check Java and look up Nextflow; no unit test may run the
    host's tools or download the launcher. Opt out with ``@pytest.mark.real_prereqs`` (the 15d
    tests build their own ``PATH``); ``needs_nextflow`` tests keep the real tools. A test's own
    fake, prepended to ``PATH`` after this autouse fixture ran, wins over these.
    """
    node = request.node
    if node.get_closest_marker("real_prereqs") or node.get_closest_marker("needs_nextflow"):
        return
    from helpers.prereqs import FAKE_PREREQS_ENV

    bindir = request.getfixturevalue("_fake_prereqs_dir")
    monkeypatch.setenv(FAKE_PREREQS_ENV, str(bindir))
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    for var in ("JAVA_HOME", "JAVA_CMD", "NXF_JAVA_HOME"):
        monkeypatch.delenv(var, raising=False)


#: the free space every unit test sees unless it fakes ``shutil.disk_usage`` itself
PLENTY_FREE = 100 * 1024**4


@pytest.fixture(autouse=True)
def _plenty_of_disk(monkeypatch: pytest.MonkeyPatch) -> None:
    """``shutil.disk_usage`` reports 100 TiB free on every path, so no unit test depends on how
    full the host's ``/tmp`` is (M11-15c2: the all-on default plan is ~16 GiB, and the 15c1 space
    check refuses it on a host with less free). The path must still exist, as for the real call.
    Tests about free space patch ``shutil.disk_usage`` themselves, which overrides this."""
    import shutil

    def disk_usage(path: object) -> Any:
        os.stat(path)  # type: ignore[arg-type]  # a missing path raises, as the real call does
        return shutil._ntuple_diskusage(PLENTY_FREE, 0, PLENTY_FREE)  # type: ignore[attr-defined]

    monkeypatch.setattr(shutil, "disk_usage", disk_usage)
