"""Step registry: one function per pipeline step, files in to files out."""

from __future__ import annotations

import importlib
import pkgutil


def load_step_modules() -> list[str]:
    """Import every module under ``giap.steps`` so its steps register; returns their names."""
    names = sorted(m.name for m in pkgutil.iter_modules(__path__))
    for name in names:
        importlib.import_module(f"{__name__}.{name}")
    return names
