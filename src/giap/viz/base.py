"""The plot registry, its skip wrapper and shared plot helpers (M9-01 behaviour 6, 8, 9;
M9-index decision 4). Pure: tables come from a ``RunView``, files are written by
``giap.io.plot_files``.

A plot is ``@plot("V07", "obs_vs_exp", needs=[...])`` over a function
``(view, cfg) -> list[Instance] | Skip``. The returned (wrapped) function checks, in this order,
``module_off:<m>``, ``layer_missing:<l>`` and ``input_absent:<key>`` before calling the plot; a
plot that raises is not caught (decision 4(c)).
"""

from __future__ import annotations

import functools
import re
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from giap.viz.style import new_figure

if TYPE_CHECKING:
    from giap.config import GiapConfig
    from giap.io.run_view import RunView


@dataclass(frozen=True)
class Instance:
    """One output file set: ``Vnn_<name>[_<suffix>].{svg,png,pdf}`` and its data TSV.

    ``data`` is 0-based half-open like ``RunView`` tables; ``coords`` names its start columns,
    which ``write_instance`` writes ``+ 1`` so the data TSV is 1-based closed (§6.5, §9.2)."""

    suffix: str | None
    figure: Figure
    data: pd.DataFrame
    title: str
    coords: tuple[str, ...] = ()


@dataclass(frozen=True)
class Skip:
    """A declared skip: ``module_off:<m>``, ``layer_missing:<l>`` or ``input_absent:<key>``."""

    reason: str


PlotFn = Callable[["RunView", "GiapConfig"], "list[Instance] | Skip"]


@dataclass(frozen=True)
class PlotSpec:
    pid: str
    name: str
    needs: tuple[str, ...]
    optional: tuple[str, ...]
    module: str
    layers: tuple[str, ...]
    fn: PlotFn


_REGISTRY: dict[str, PlotSpec] = {}
_PID = re.compile(r"^[Vv](\d{1,2})$")


def _pid(text: str) -> str | None:
    m = _PID.match(text.strip())
    return f"V{int(m.group(1)):02d}" if m else None


def _skip_reason(spec: PlotSpec, view: RunView) -> str | None:
    if spec.module != "core" and spec.module in view.modules_off():
        return f"module_off:{spec.module}"
    flags = view.run_flags()
    for layer in spec.layers:
        if f"layer_missing:{layer}" in flags:
            return f"layer_missing:{layer}"
    for key in spec.needs:
        if not view.has(key):
            return f"input_absent:{key}"
    return None


def plot(
    pid: str,
    name: str,
    *,
    needs: Sequence[str] = (),
    optional: Sequence[str] = (),
    module: str = "core",
    layers: Sequence[str] = (),
) -> Callable[[PlotFn], PlotFn]:
    """Register a plot under ``pid`` (``Vnn``); returns the skip-checking wrapper."""
    norm = _pid(pid)
    if norm != pid:
        raise ValueError(f"plot ID {pid!r} is not of the form Vnn")

    def deco(fn: PlotFn) -> PlotFn:
        if pid in _REGISTRY:
            raise ValueError(f"plot ID {pid} is already registered ({_REGISTRY[pid].name})")

        @functools.wraps(fn)
        def wrapped(view: RunView, cfg: GiapConfig) -> list[Instance] | Skip:
            reason = _skip_reason(spec, view)
            return Skip(reason) if reason is not None else fn(view, cfg)

        spec = PlotSpec(pid, name, tuple(needs), tuple(optional), module, tuple(layers), wrapped)
        _REGISTRY[pid] = spec
        return wrapped

    return deco


def registered() -> Mapping[str, PlotSpec]:
    """Every registered plot, sorted by ID."""
    return {k: _REGISTRY[k] for k in sorted(_REGISTRY)}


def normalise_ids(text: str) -> list[str]:
    """``"V1, v18,V03"`` -> ``["V01", "V03", "V18"]``; an unknown ID raises ``ValueError``."""
    out: set[str] = set()
    for token in (t.strip() for t in text.split(",")):
        if not token:
            continue
        norm = _pid(token)
        if norm is None or norm not in _REGISTRY:
            raise ValueError(f"unknown plot ID {norm or token}")
        out.add(norm)
    return sorted(out)


def empty_state(
    pid: str,
    name: str,
    message: str,
    columns: Sequence[str],
    width: Literal["single", "double"],
) -> Instance:
    """An axes-free figure with ``message`` centred and a header-only data frame (decision 4(b))."""
    fig = new_figure(width, 1.2)
    fig.set_label(f"{pid}_{name}")
    ax = fig.add_subplot()
    ax.set_axis_off()
    ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes, color="#555555")
    data = pd.DataFrame({c: pd.Series([], dtype=object) for c in columns}, columns=list(columns))
    return Instance(suffix=None, figure=fig, data=data, title=message)


def binned(
    starts: np.ndarray, ends: np.ndarray, length: int, max_elements: int = 20_000
) -> np.ndarray | None:
    """``None`` when a track has at most ``max_elements`` elements, otherwise the ``n + 1``
    integer edges (``int64``) of ``n = min(max_elements, length)`` bins over ``[0, length)``;
    every bin is at least 1 bp and widths differ by at most 1 bp (decision 12, PR #191 B2).

    Elements must lie in ``[0, length)``: the caller splits circular features with
    ``end > length`` at the origin first."""
    if len(starts) != len(ends):
        raise ValueError(f"{len(starts)} starts but {len(ends)} ends")
    if len(starts) <= max_elements:
        return None
    if length < 1:
        raise ValueError(f"replicon length {length} < 1")
    n = min(max_elements, length)
    return (np.arange(n + 1, dtype=np.int64) * np.int64(length)) // np.int64(n)


def title_pct(n: int, d: int) -> str:
    """``n / d`` as a whole percent for titles (``"12 %"``); ``d == 0`` -> ``"NA"``."""
    if d == 0:
        return "NA"
    return f"{int(100 * n / d + 0.5)} %"
