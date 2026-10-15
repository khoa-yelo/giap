"""Calibration harness fixtures and suite selection (issue M4-05; GIAP_SPEC.md §8.4, §11.2).

- ``full``: true for tests marked ``calibration_full`` (nightly R), false otherwise.
- Every calibration test gets the keyword ``calsuite_<suite>`` for each ``tolerances.yaml`` suite
  it names (in its own source, else in its module), so the job array selects one suite with
  ``-k calsuite_<suite>`` (``_criteria.pytest_args``).
- Results are appended by ``_criteria.log_result`` to ``local/logs/calibration/<suite>.jsonl``.
"""

from __future__ import annotations

import inspect
import re
from functools import cache
from types import ModuleType

import _criteria as C
import pytest


@pytest.fixture
def full(request: pytest.FixtureRequest) -> bool:
    """True under ``-m calibration_full`` (the test carries the marker)."""
    return request.node.get_closest_marker("calibration_full") is not None


def _named_suites(source: str) -> frozenset[str]:
    return frozenset(
        s for s in C.array_suites() if re.search(rf"[\"']{re.escape(s)}[\"'.]", source)
    )


@cache
def _module_suites(module: ModuleType) -> frozenset[str]:
    try:
        return _named_suites(inspect.getsource(module))
    except (OSError, TypeError):
        return frozenset()


def _item_suites(item: pytest.Function) -> frozenset[str]:
    try:
        own = _named_suites(inspect.getsource(item.function))
    except (OSError, TypeError):
        own = frozenset()
    return own or _module_suites(item.module)


@pytest.hookimpl(tryfirst=True)
def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    for item in items:
        if not isinstance(item, pytest.Function):
            continue
        if not any(item.iter_markers("calibration")) and not any(
            item.iter_markers("calibration_full")
        ):
            continue
        for suite in _item_suites(item):
            item.extra_keyword_matches.add(f"calsuite_{suite}")
