"""Load data-building helpers from the unit-test modules (``tests/`` is not a package).

The benchmarks for budgets moved out of the per-commit unit suite reuse the exact inputs the unit
tests build, so the two cannot drift apart.
"""

from __future__ import annotations

import importlib.util
import sys
import time
from collections.abc import Callable
from pathlib import Path
from types import ModuleType
from typing import Any

_ROOT = Path(__file__).resolve().parent.parent


def load_test_module(relpath: str) -> ModuleType:
    """Import ``tests/<relpath>`` under a private module name and return it."""
    name = "_bench_" + relpath.replace("/", "_").removesuffix(".py")
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _ROOT / "tests" / relpath)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def timed(fn: Callable[[], Any]) -> tuple[float, Any]:
    """Wall-clock seconds of one call of ``fn`` and its return value."""
    t0 = time.perf_counter()
    out = fn()
    return time.perf_counter() - t0, out
