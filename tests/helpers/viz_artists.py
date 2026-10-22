"""Artist lookups for the M9 plot tests: markers, texts, rectangles (no image comparison).

A "marker" is one drawn point: a ``Line2D`` vertex with a visible marker (not the ``|``/``_``
errorbar caps) or one ``PathCollection`` offset (``scatter``). Colours are RGBA tuples; a hollow
marker has a face alpha of 0 (``facecolor="none"``).
"""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.colors as mcolors
import numpy as np
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import PathCollection, PolyCollection
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.text import Text

_NO_MARKER = (None, "", "None", " ", "none", "|", "_")
_CLEAR = (0.0, 0.0, 0.0, 0.0)


@dataclass(frozen=True)
class Marker:
    x: float
    y: float
    face: tuple[float, float, float, float]
    edge: tuple[float, float, float, float]

    @property
    def hollow(self) -> bool:
        return self.face[3] == 0.0

    def colours(self) -> set[str]:
        """The hex of the face and of the edge (whichever is visible)."""
        out = set()
        for rgba in (self.face, self.edge):
            if rgba[3] > 0:
                out.add(mcolors.to_hex(rgba).upper())
        return out


def _rgba(value: object) -> tuple[float, float, float, float]:
    if isinstance(value, str) and value.lower() == "none":
        return _CLEAR
    r, g, b, a = mcolors.to_rgba(value)  # type: ignore[arg-type]
    return (float(r), float(g), float(b), float(a))


def _pick(arr: np.ndarray, i: int) -> tuple[float, float, float, float]:
    if len(arr) == 0:
        return _CLEAR
    row = arr[i if len(arr) > 1 else 0]
    return (float(row[0]), float(row[1]), float(row[2]), float(row[3]))


def markers(ax: Axes) -> list[Marker]:
    out: list[Marker] = []
    for line in ax.get_lines():
        if line.get_marker() in _NO_MARKER or not line.get_visible():
            continue
        face = _rgba(line.get_markerfacecolor())
        edge = _rgba(line.get_markeredgecolor())
        xs = np.asarray(line.get_xdata(), dtype=float)
        ys = np.asarray(line.get_ydata(), dtype=float)
        out += [Marker(float(x), float(y), face, edge) for x, y in zip(xs, ys, strict=True)]
    for coll in ax.collections:
        if not isinstance(coll, PathCollection) or not coll.get_visible():
            continue
        offs = np.asarray(coll.get_offsets(), dtype=float).reshape(-1, 2)
        faces = np.asarray(coll.get_facecolors(), dtype=float).reshape(-1, 4)
        edges = np.asarray(coll.get_edgecolors(), dtype=float).reshape(-1, 4)
        for i, (x, y) in enumerate(offs):
            out.append(Marker(float(x), float(y), _pick(faces, i), _pick(edges, i)))
    return [m for m in out if np.isfinite(m.x) and np.isfinite(m.y)]


def markers_at_x(ax: Axes, x: float, tol: float = 1e-9) -> list[Marker]:
    return [m for m in markers(ax) if abs(m.x - x) <= tol]


def polys(ax: Axes) -> list[PolyCollection]:
    """Violin bodies (``violinplot``) and any other filled polygons."""
    return [c for c in ax.collections if isinstance(c, PolyCollection)]


def rects(ax: Axes) -> list[Rectangle]:
    return [
        p
        for p in ax.patches
        if isinstance(p, Rectangle) and p.get_width() > 0 and p.get_height() > 0
    ]


def texts(fig: Figure) -> list[str]:
    """Every non-empty text of the figure after a draw (tick labels included)."""
    draw(fig)
    return [t.get_text() for t in fig.findobj(Text) if t.get_text() and t.get_visible()]


def draw(fig: Figure) -> FigureCanvasAgg:
    canvas = fig.canvas if isinstance(fig.canvas, FigureCanvasAgg) else FigureCanvasAgg(fig)
    canvas.draw()
    return canvas


def texts_top_to_bottom(fig: Figure, contains: str) -> list[str]:
    """Visible texts containing ``contains``, ordered by their screen position, top first."""
    canvas = draw(fig)
    renderer = canvas.get_renderer()
    found = []
    for t in fig.findobj(Text):
        if contains in t.get_text() and t.get_visible():
            box = t.get_window_extent(renderer)
            found.append((-(box.y0 + box.y1) / 2, t.get_text()))
    return [s for _, s in sorted(found)]


def texts_left_to_right(fig: Figure, ax: Axes, labels: set[str]) -> list[str]:
    """The x tick label texts of ``ax`` (only those in ``labels`` when given), left to right."""
    canvas = draw(fig)
    renderer = canvas.get_renderer()
    found = []
    for t in ax.get_xticklabels():
        if t.get_text() and t.get_visible():
            box = t.get_window_extent(renderer)
            found.append(((box.x0 + box.x1) / 2, t.get_text()))
    return [s for _, s in sorted(found) if not labels or s in labels]
