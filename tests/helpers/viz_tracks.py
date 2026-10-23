"""Per-track artist lookups for the V1 ``circos`` tests (M9-08; no image comparison).

A track's artists carry ``Artist.get_gid()`` equal to the track's data-TSV enum value
(``gene_plus``, ``region``, ``island`` …). One drawn **element** is one ``Patch`` (a pyCirclize
``rect``/``arrow``/``bar`` piece), one path of a ``PatchCollection``/``PolyCollection``, one
offset of a ``PathCollection`` (``scatter``), one segment of a ``LineCollection``, one vertex of a
``Line2D`` with a marker, or one ``Line2D`` without a marker (a tick or a curve).

pyCirclize (1.10.1) merges every plain patch (``zorder == 1``, no hatch) into one untagged
``PatchCollection``; a renderer that keeps track gids must give its patches a hatch, a
non-default ``zorder`` or its own tagged collections.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import matplotlib.colors as mcolors
import numpy as np
from matplotlib.artist import Artist
from matplotlib.collections import Collection, LineCollection, PathCollection
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

_NO_MARKER = (None, "", "None", " ", "none")
RGBA = tuple[float, float, float, float]
_CLEAR: RGBA = (0.0, 0.0, 0.0, 0.0)


@dataclass(frozen=True)
class Element:
    face: RGBA
    edge: RGBA
    hatch: str | None

    @property
    def face_hex(self) -> str | None:
        return None if self.face[3] == 0 else mcolors.to_hex(self.face).upper()

    @property
    def hollow(self) -> bool:
        return self.face[3] == 0

    def colours(self) -> set[str]:
        """The hex of the face and of the edge, whichever is visible."""
        return {mcolors.to_hex(c).upper() for c in (self.face, self.edge) if c[3] > 0}


def _rgba(value: object) -> RGBA:
    if value is None or (isinstance(value, str) and value.lower() == "none"):
        return _CLEAR
    r, g, b, a = mcolors.to_rgba(value)  # type: ignore[arg-type]
    return (float(r), float(g), float(b), float(a))


def _row(arr: np.ndarray, i: int) -> RGBA:
    arr = np.asarray(arr, dtype=float).reshape(-1, 4)
    if len(arr) == 0:
        return _CLEAR
    r = arr[i if len(arr) > 1 else 0]
    return (float(r[0]), float(r[1]), float(r[2]), float(r[3]))


def _hatch(artist: Patch | Collection) -> str | None:
    h = artist.get_hatch()
    return h or None


def _of(artist: Artist) -> list[Element]:
    if not artist.get_visible():
        return []
    if isinstance(artist, Patch):
        return [
            Element(_rgba(artist.get_facecolor()), _rgba(artist.get_edgecolor()), _hatch(artist))
        ]
    if isinstance(artist, LineCollection):
        n = len(artist.get_segments())
        cols = artist.get_colors()
        return [Element(_CLEAR, _row(cols, i), None) for i in range(n)]
    if isinstance(artist, PathCollection):
        n = len(np.asarray(artist.get_offsets()).reshape(-1, 2))
        faces, edges = artist.get_facecolors(), artist.get_edgecolors()
        return [Element(_row(faces, i), _row(edges, i), _hatch(artist)) for i in range(n)]
    if isinstance(artist, Collection):
        n = len(artist.get_paths())
        faces, edges = artist.get_facecolors(), artist.get_edgecolors()
        return [Element(_row(faces, i), _row(edges, i), _hatch(artist)) for i in range(n)]
    if isinstance(artist, Line2D):
        if artist.get_marker() in _NO_MARKER:
            return [Element(_CLEAR, _rgba(artist.get_color()), None)]
        n = len(np.asarray(artist.get_xdata()))
        face = _rgba(artist.get_markerfacecolor())
        edge = _rgba(artist.get_markeredgecolor())
        return [Element(face, edge, None) for _ in range(n)]
    return []


def tagged(fig: Figure) -> dict[str, list[Artist]]:
    """gid -> the figure's artists carrying it."""
    out: dict[str, list[Artist]] = {}
    for a in fig.findobj(lambda x: isinstance(x, Artist) and x.get_gid() is not None):
        out.setdefault(str(a.get_gid()), []).append(a)
    return out


def elements(fig: Figure, track: str) -> list[Element]:
    """Every drawn element of ``track`` (artists whose gid is ``track``)."""
    found: list[Element] = []
    for a in tagged(fig).get(track, []):
        found += _of(a)
    return found


def face_hexes(els: Iterable[Element]) -> list[str | None]:
    return sorted((e.face_hex for e in els), key=lambda h: h or "")
