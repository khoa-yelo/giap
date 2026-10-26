"""Shared fixtures for unit tests.

``registry`` gives each test an empty step registry (M2-02 "temporary registry"): the module-level
``giap.steps.registry._REGISTRY`` mapping (name -> StepSpec) is swapped for an empty dict, so toy
steps never collide with the real annotation steps or with other tests.
"""

from __future__ import annotations

from collections.abc import Iterator
from types import ModuleType

import pytest


@pytest.fixture
def registry(monkeypatch: pytest.MonkeyPatch) -> Iterator[ModuleType]:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

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
