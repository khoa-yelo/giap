"""Markers used only by the integration layer (registered here; pyproject's list is strict)."""

from __future__ import annotations

import pytest


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers", "integration: a whole `giap run` on sim100k (M9-13; tests/integration only)"
    )
