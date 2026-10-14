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
    from giap.steps import registry as reg

    monkeypatch.setattr(reg, "_REGISTRY", {})
    yield reg
