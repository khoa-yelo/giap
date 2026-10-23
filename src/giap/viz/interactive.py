"""Plotly twins of the static plots (GIAP_SPEC.md §10.2 "Interactive in report"; M9-14, M9-15).

A twin is ``@interactive("V17")`` over ``(data, title) -> plotly Figure`` or ``(data, title,
meta) -> Figure``, where ``data`` is the plot's data TSV (the exact drawn data, so the static and
interactive views cannot disagree) and ``meta`` the optional run context built by
``giap.steps.plots.twins`` (``{"alpha": {name: level}, "replicons": {id: length}}`` in FASTA
order; M9-15 lead decision 5). Twins live in the module-level ``_REGISTRY`` (read at call time).
Pure: no file I/O.
"""

from __future__ import annotations

import inspect
from collections.abc import Callable, Mapping, Sequence
from typing import Any

import pandas as pd

InteractiveFn = Callable[..., Any]  # (data, title[, meta]) -> plotly.graph_objects.Figure

#: §10.2 plots whose "Interactive in report" column is Plotly (V16: its concordance view); the
#: twins of V14, V15 and V16 land with those plots (M9-19, M9-20, M9-23)
INTERACTIVE_PLOTS: tuple[str, ...] = (
    *("V02", "V05", "V06", "V07", "V08", "V10", "V11", "V14", "V15", "V16", "V17", "V18"),
)

_REGISTRY: dict[str, InteractiveFn] = {}


def interactive(pid: str, columns: Sequence[str] = ()) -> Callable[[InteractiveFn], InteractiveFn]:
    """Register the decorated function as the twin of plot ``pid``; one twin per plot. A data
    TSV without every one of ``columns`` gets no twin (not that plot's TSV)."""

    def register(fn: InteractiveFn) -> InteractiveFn:
        if pid in _REGISTRY:
            raise ValueError(f"plot {pid} already has an interactive twin")
        fn.giap_columns = tuple(columns)  # type: ignore[attr-defined]
        _REGISTRY[pid] = fn
        return fn

    return register


def has_twin(pid: str) -> bool:
    return pid in _REGISTRY


def figure_json(
    pid: str, data: pd.DataFrame, title: str, meta: Mapping[str, Any] | None = None
) -> str | None:
    """The twin's figure as a JSON object string (``data``, ``layout``); ``None`` without a twin,
    for a header-only data TSV (the static empty state) or one without the twin's columns. A
    twin that raises is not caught (the ``report`` step fails, as for static plots)."""
    fn = _REGISTRY.get(pid)
    if fn is None or not len(data) or not set(getattr(fn, "giap_columns", ())) <= set(data):
        return None
    takes_meta = len(inspect.signature(fn).parameters) >= 3
    fig = fn(data, title, meta) if takes_meta else fn(data, title)
    text: str = fig.to_json(validate=True, pretty=False, remove_uids=True)
    return text
