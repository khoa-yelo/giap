"""Plotly twins of the static plots (GIAP_SPEC.md §10.2 "Interactive in report"; M9-14).

A twin is ``@interactive("V17")`` over ``(data, title) -> plotly Figure``, where ``data`` is the
plot's data TSV (the exact drawn data, so the static and interactive views cannot disagree).
Twins live in the module-level ``_REGISTRY`` (read at call time). Pure: no file I/O.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import pandas as pd

InteractiveFn = Callable[[pd.DataFrame, str], Any]  # -> plotly.graph_objects.Figure

_REGISTRY: dict[str, InteractiveFn] = {}


def interactive(pid: str) -> Callable[[InteractiveFn], InteractiveFn]:
    """Register the decorated function as the twin of plot ``pid``; one twin per plot."""

    def register(fn: InteractiveFn) -> InteractiveFn:
        if pid in _REGISTRY:
            raise ValueError(f"plot {pid} already has an interactive twin")
        _REGISTRY[pid] = fn
        return fn

    return register


def has_twin(pid: str) -> bool:
    return pid in _REGISTRY


def figure_json(pid: str, data: pd.DataFrame, title: str) -> str | None:
    """The twin's figure as a JSON object string (``data``, ``layout``); ``None`` without a twin.
    A twin that raises is not caught (the ``report`` step fails, as for static plots)."""
    fn = _REGISTRY.get(pid)
    if fn is None:
        return None
    text: str = fn(data, title).to_json(validate=True, pretty=False, remove_uids=True)
    return text
