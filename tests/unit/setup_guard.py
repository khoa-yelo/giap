"""No-network, no-tool guard for the ``giap setup`` tests (M2-08 lead decision 3).

Import the autouse fixture into a test module to apply it to that module only::

    from setup_guard import setup_guard  # noqa: F401

While it is active, ``ToolRunner.run``, ``subprocess.run``/``Popen``, ``HttpsFetcher.fetch``,
socket connections and DNS lookups raise ``GuardError`` (and are recorded in ``BLOCKED``) unless
a test injects a fake explicitly. ``ToolRunner.run`` still calls the real ``ToolRunner.command``
first, so a missing environment/image surfaces as ``ToolUnavailable`` as in production; only the
start of the tool is blocked. Tests that exercise the real ``HttpsFetcher.fetch`` against fake
``http.client`` connections request the ``real_https_fetch`` fixture, which restores only that
method (sockets and DNS stay blocked).
"""

from __future__ import annotations

import socket
import subprocess
from collections.abc import Iterator, Sequence
from typing import Any

import pytest

BLOCKED: list[str] = []
_REAL: dict[str, Any] = {}


class GuardError(AssertionError):
    """A setup test tried to reach the network or start a tool."""


def _deny(what: str) -> Any:
    def deny(*args: Any, **_kw: Any) -> Any:
        BLOCKED.append(f"{what}: {args[1:] if args else ()}")
        raise GuardError(f"setup_guard: {what} is blocked in setup tests ({args!r})")

    return deny


@pytest.fixture(autouse=True)
def setup_guard(monkeypatch: pytest.MonkeyPatch) -> Iterator[list[str]]:
    from giap.setupdb import HttpsFetcher, ToolRunner

    BLOCKED.clear()
    _REAL.setdefault("fetch", HttpsFetcher.fetch)

    def guarded_run(self: Any, tool: str, argv: Sequence[str]) -> int:
        self.command(tool, argv)  # ToolUnavailable propagates as in production
        BLOCKED.append(f"ToolRunner.run: {tool} {list(argv)}")
        raise GuardError(f"setup_guard: ToolRunner.run({tool!r}) is blocked in setup tests")

    monkeypatch.setattr(ToolRunner, "run", guarded_run)
    monkeypatch.setattr(HttpsFetcher, "fetch", _deny("HttpsFetcher.fetch"))
    monkeypatch.setattr(subprocess, "run", _deny("subprocess.run"))
    monkeypatch.setattr(subprocess, "Popen", _deny("subprocess.Popen"))
    monkeypatch.setattr(socket, "create_connection", _deny("socket.create_connection"))
    monkeypatch.setattr(socket.socket, "connect", _deny("socket.connect"))
    monkeypatch.setattr(socket.socket, "connect_ex", _deny("socket.connect_ex"))
    for name in ("getaddrinfo", "gethostbyname", "gethostbyname_ex"):  # no DNS either
        monkeypatch.setattr(socket, name, _deny(f"socket.{name}"))
    yield BLOCKED


@pytest.fixture
def real_https_fetch(setup_guard: list[str], monkeypatch: pytest.MonkeyPatch) -> None:
    """Restore only the real ``HttpsFetcher.fetch``; sockets stay blocked (fake http.client)."""
    from giap.setupdb import HttpsFetcher

    for var in ("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY", "all_proxy",
                "ALL_PROXY", "no_proxy", "NO_PROXY"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setattr(HttpsFetcher, "fetch", _REAL["fetch"])
