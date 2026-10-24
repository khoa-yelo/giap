"""The GIAP_SPEC.md §14.4 style checklist as code (M9-24; §10.1; M9-index decisions 30-31).

``check_instance`` returns one ``CheckResult`` per item, in the order of the M9-24 table; the
colour-blind item is recorded (``passed = None``), never gated (lead decision 1). Deuteranopia is
simulated with the Machado, Oliveira & Fernandes 2009 matrix (IEEE TVCG 15(6), 1291-1298, Table 1,
deuteranomaly at severity 1.0), applied in linear RGB. Chart furniture (spines, grid, ticks,
legends, ``axhline``/``axvline`` reference lines) is exempt from the palette item (decision 9).
File access goes through ``giap.io``: the data TSV is read back and the simulated PNG is written
next to the instance's PNG as ``deuteranopia_<stem>.png`` (a name the report never picks up).
"""

from __future__ import annotations

import io
import itertools
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np
from matplotlib.axes import Axes
from matplotlib.axis import Axis
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.cm import ScalarMappable
from matplotlib.colors import to_hex, to_rgba_array
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.text import Text

from giap.io.plot_files import read_png, write_png
from giap.io.tsv import read_tsv
from giap.viz import circos, detection, style
from giap.viz.base import Instance

ITEMS: tuple[str, ...] = (
    "title states the finding",
    "axes carry units",
    "palette as §10.1",
    "no overlapping labels",
    "legible at single-column width",
    "data TSV written",
    "colour-blind check",
    "≤ 20,000 drawn elements",
)
#: Machado et al. 2009, Table 1, deuteranomaly, severity 1.0 (rows act on linear RGB columns)
MACHADO_DEUTERANOPIA: np.ndarray = np.array(
    [
        [0.367322, 0.860646, -0.227968],
        [0.280085, 0.672501, 0.047413],
        [-0.011820, 0.042940, 0.968881],
    ]
)
#: axis labels (or label endings) that need no unit: unitless statistics, counts and categories
#: (lead decision 7; extend as plots need it). Anything else needs a unit in parentheses.
UNITLESS: tuple[str, ...] = (
    "log2 abundance", "q", "lfdr", "local fdr", "fdr", "\u03c1", "nes", "p-value", "p_detect",
    "counts", "count", "checks", "genes", "peptides", "placements", "density", "fraction",
    "rate", "category", "score bin", "tus",
    "max |s_w| over windows",  # V18 (c): the standardised scan statistic (M9-21 N2)
)  # fmt: skip
MIN_PT = 6.0  # index decision 31, lead decision 2
WIDTHS_IN = (style.SINGLE_IN, style.DOUBLE_IN)
MAX_ELEMENTS = 20_000
CONTINUOUS_MAPS = ("viridis", "cividis")
_NEUTRAL = ("#000000", "#FFFFFF", "#BBBBBB")
_PALETTES = (
    style.OKABE_ITO, style.DETECTION, style.STATUS_COLOUR, style.REGION, style.COG_SUPER_COLOUR,
)  # fmt: skip
_TABLE = {c.upper() for p in _PALETTES for c in p.values()}
ALLOWED: frozenset[str] = frozenset(
    _TABLE | {c.upper() for c, _ in style.ISLAND.values()} | set(_NEUTRAL)
)
#: colours drawn straight from a continuous map (``cmap(norm(x))``) are LUT entries
_LUT = frozenset(
    to_hex(c).upper() for m in CONTINUOUS_MAPS for c in matplotlib.colormaps[m](np.arange(256))
)
_DET = frozenset(c.upper() for c in style.DETECTION.values())
_COG = frozenset(c.upper() for c in style.COG_SUPER_COLOUR.values())
_TOL_PX = 0.5  # boxes that only touch (float round-off) do not overlap
_NA = "not applicable (empty state)"


@dataclass(frozen=True)
class CheckResult:
    plot: str
    item: str
    passed: bool | None  # None = recorded, not gated (or not applicable)
    detail: str


#: Titles that state a finding without a number (lead decision 10): V08 with no measured abundance
#: plots a detected fraction, so its title names what is shown instead.
DECLARED_TITLES: tuple[str, ...] = (detection.NO_ABUNDANCE,)


def is_empty_state(inst: Instance) -> bool:
    """What ``base.empty_state`` produces: a header-only frame, axes-free, the title drawn as
    the only text (lead decision 6: the phrases are derived from ``base``, not pinned)."""
    if len(inst.data) != 0:
        return False
    axes = inst.figure.axes
    texts = [t.get_text() for ax in axes for t in ax.texts]
    return all(not ax.axison for ax in axes) and texts == [inst.title]


def check_instance(inst: Instance, pid: str, files: Sequence[Path]) -> list[CheckResult]:
    """Run every item on one written instance (``files`` as ``write_instance`` returns them)."""
    images = [p for p in files if p.suffix in (".svg", ".png", ".pdf")]
    name = images[0].stem if images else pid
    fig = inst.figure
    canvas = FigureCanvasAgg(fig)
    canvas.draw()  # type: ignore[no-untyped-call]
    texts = _texts(fig, canvas)
    empty = is_empty_state(inst)
    per_axes, maps = _drawn_colours(fig)
    colours = set().union(*per_axes)
    checks: list[tuple[bool | None, str]] = [
        _title(inst, empty),
        (None, _NA) if empty else _units(fig),
        (None, _NA) if empty else _palette(per_axes, maps),
        _overlap(texts),
        _legible(fig, texts),
        _tsv(files, pid, name, empty),
        (None, _cvd(inst, files, name, colours)),
        (None, _NA) if empty else _elements(fig),
    ]
    return [CheckResult(name, item, ok, d) for item, (ok, d) in zip(ITEMS, checks, strict=True)]


def score(results: Sequence[CheckResult]) -> float:
    """Passed over gated (``passed`` not ``None``) items; 1.0 when nothing is gated."""
    gated = [r.passed for r in results if r.passed is not None]
    return sum(1 for g in gated if g) / len(gated) if gated else 1.0


# ------------------------------------------------------------------ items


def _title(inst: Instance, empty: bool) -> tuple[bool, str]:
    t = inst.title.strip()
    if not t:
        return False, "empty title"
    if "express" in t.lower():
        return False, f"{t!r} says 'express'"
    if re.search(r"\d", t) or empty or t.startswith(DECLARED_TITLES):
        return True, t
    return False, f"{t!r} states no number or empty-state message"


def _needs_unit(label: str) -> bool:
    low = label.strip().lower()
    if not low or re.search(r"\(.+\)", low):
        return False
    return not any(re.search(rf"(?:^|[\s|]){re.escape(u)}$", low) for u in UNITLESS)


def _units(fig: Figure) -> tuple[bool, str]:
    bad = [
        lab.get_text()
        for ax in _visible_axes(fig)
        for lab in (ax.xaxis.label, ax.yaxis.label)
        if lab.get_visible() and _needs_unit(lab.get_text())
    ]
    return (not bad, "; ".join(f"{b!r} needs a unit" for b in bad) or "ok")


def _palette(per_axes: list[set[str]], maps: set[str]) -> tuple[bool, str]:
    colours = set().union(*per_axes)
    problems = [f"{c} not in §10.1" for c in sorted(colours - ALLOWED - _LUT)]
    base_maps = {m.removesuffix("_r") for m in maps}
    problems += [f"colourmap {m} (use viridis or cividis)" for m in sorted(base_maps)
                 if m not in CONTINUOUS_MAPS]  # fmt: skip
    # M9-01 Q4: one axes draws the whole COG palette and one (maybe the same) the detection one
    if any(c >= _DET for c in per_axes) and any(c >= _COG for c in per_axes):
        problems.append("COG and detection palettes in one figure")
    return (not problems, "; ".join(problems) or "ok")


def _overlap(texts: list[tuple[Text, np.ndarray]]) -> tuple[bool, str]:
    pairs = []
    for (a, ba), (b, bb) in itertools.combinations(texts, 2):
        w = min(ba[2], bb[2]) - max(ba[0], bb[0])
        h = min(ba[3], bb[3]) - max(ba[1], bb[1])
        if w > _TOL_PX and h > _TOL_PX:
            pairs.append(f'"{a.get_text()}" ∩ "{b.get_text()}"')
    more = f" (+{len(pairs) - 5} more)" if len(pairs) > 5 else ""
    return (not pairs, "; ".join(pairs[:5]) + more if pairs else "ok")


def _legible(fig: Figure, texts: list[tuple[Text, np.ndarray]]) -> tuple[bool, str]:
    problems = [f'"{t.get_text()}" {t.get_fontsize():g} pt' for t, _ in texts
                if float(t.get_fontsize()) < MIN_PT]  # fmt: skip
    width = fig.get_figwidth()
    if not any(abs(width - w) < 1e-6 for w in WIDTHS_IN):
        problems.insert(0, f"width {width:g} in (not 3.5 or 7.2)")
    return (not problems, "; ".join(problems) or "ok")


def _tsv(files: Sequence[Path], pid: str, name: str, empty: bool) -> tuple[bool, str]:
    want = [p for p in files if p.suffix == ".tsv" and p.parent.name == "data"
            and p.stem == name and p.stem.startswith(f"{pid}_")]  # fmt: skip
    if not want or not want[0].is_file():
        return False, f"data/{name}.tsv missing"
    _, rows = read_tsv(want[0])
    if not rows and not empty:
        return False, f"data/{name}.tsv has no rows"
    return True, f"data/{name}.tsv: {len(rows)} rows"


def _elements(fig: Figure) -> tuple[bool, str]:
    """Counts per axes and, within one axes, per V1 track ``gid`` (``circos.TRACK_GIDS``; M9-24:
    V1's polar axes holds every track). Any other gid and untagged artists count together as one
    group per axes (M9-26 lead decision 4, B1); every group is gated at ``MAX_ELEMENTS``."""
    worst: list[str] = []
    for i, ax in enumerate(fig.axes):
        counts: dict[str | None, int] = {}

        def add(gid: str | None, n: int, counts: dict[str | None, int] = counts) -> None:
            key = gid if gid in circos.TRACK_GIDS else None
            counts[key] = counts.get(key, 0) + n

        for a in [*ax.patches, *ax.lines, *ax.texts]:
            add(a.get_gid(), 1)
        for c in ax.collections:
            add(c.get_gid(), max(len(np.asarray(c.get_offsets())), len(c.get_paths()), 1))
        for gid, n in counts.items():
            if n > MAX_ELEMENTS:
                where = f"axes {i}" if gid is None else f"axes {i} track {gid}"
                worst.append(f"{where}: {n:,} elements")
    return (not worst, "; ".join(worst) or "ok")


def _cvd(inst: Instance, files: Sequence[Path], name: str, colours: set[str]) -> str:
    """Writes the simulated PNG and returns ``min ΔE before X.X → after Y.Y`` (decision 8)."""
    png = next((p for p in files if p.suffix == ".png"), None)
    if png is not None:
        img = read_png(png)
    else:
        buf = io.BytesIO()
        inst.figure.savefig(buf, format="png", dpi=100)
        buf.seek(0)
        img = read_png(buf)
    sim = img.copy()
    sim[..., :3] = simulate_deuteranopia(img[..., :3])
    folder = png.parent if png is not None else next(iter(files), Path(".")).parent
    write_png(folder / f"deuteranopia_{name}.png", sim)
    cat = sorted(colours - {"#000000", "#FFFFFF"} - _LUT)
    if len(cat) < 2:
        return "min ΔE before NA → after NA (fewer than 2 categorical colours)"
    rgb = to_rgba_array(cat)[:, :3]
    before, after = _min_delta_e(rgb), _min_delta_e(simulate_deuteranopia(rgb))
    return f"min ΔE before {before:.1f} → after {after:.1f}"


# ------------------------------------------------------------------ colour science


def _to_linear(c: np.ndarray) -> np.ndarray:
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _to_srgb(c: np.ndarray) -> np.ndarray:
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * np.power(c, 1 / 2.4) - 0.055)


def simulate_deuteranopia(srgb: np.ndarray) -> np.ndarray:
    """Float sRGB in [0, 1], shape ``(..., 3)``: the Machado matrix in linear RGB, clipped."""
    lin = _to_linear(np.clip(np.asarray(srgb, dtype=float), 0.0, 1.0))
    out: np.ndarray = np.clip(_to_srgb(np.clip(lin @ MACHADO_DEUTERANOPIA.T, 0.0, 1.0)), 0.0, 1.0)
    return out


_XYZ = np.array(
    [[0.4124564, 0.3575761, 0.1804375],
     [0.2126729, 0.7151522, 0.0721750],
     [0.0193339, 0.1191920, 0.9503041]]
)  # fmt: skip
_WHITE = np.array([0.95047, 1.0, 1.08883])  # D65


def _lab(srgb: np.ndarray) -> np.ndarray:
    xyz = (_to_linear(srgb) @ _XYZ.T) / _WHITE
    f = np.where(xyz > (6 / 29) ** 3, np.cbrt(xyz), xyz / (3 * (6 / 29) ** 2) + 4 / 29)
    return np.stack([116 * f[:, 1] - 16, 500 * (f[:, 0] - f[:, 1]), 200 * (f[:, 1] - f[:, 2])], 1)


def _min_delta_e(srgb: np.ndarray) -> float:
    """Minimum pairwise CIE76 ΔE (Euclidean distance in CIELAB)."""
    lab = _lab(srgb)
    d = np.linalg.norm(lab[:, None, :] - lab[None, :, :], axis=-1)
    return float(d[np.triu_indices(len(lab), 1)].min())


# ------------------------------------------------------------------ figure walking


def _visible_axes(fig: Figure) -> list[Axes]:
    return [ax for ax in fig.axes if ax.get_visible() and ax.axison]


def _texts(fig: Figure, canvas: FigureCanvasAgg) -> list[tuple[Text, np.ndarray]]:
    """Every visible, non-empty text artist with its rendered (rotated) box in pixels."""
    found: list[Text] = list(fig.texts)
    found += [t for leg in fig.legends for t in (*leg.get_texts(), leg.get_title())]
    sup = getattr(fig, "_suptitle", None)
    if isinstance(sup, Text):
        found.append(sup)
    for ax in fig.axes:
        if not ax.get_visible():
            continue
        found += [ax.title, ax._left_title, ax._right_title, *ax.texts]  # type: ignore[attr-defined]
        leg = ax.get_legend()
        if leg is not None and leg.get_visible():
            found += [*leg.get_texts(), leg.get_title()]
        if ax.axison:
            found += [ax.xaxis.label, ax.yaxis.label]
            found += [t for axis in (ax.xaxis, ax.yaxis) for t in _tick_labels(axis)]
    renderer = canvas.get_renderer()  # type: ignore[no-untyped-call]
    out = []
    for t in dict.fromkeys(found):  # de-duplicate, keep order
        if t.get_visible() and t.get_text().strip():
            e = t.get_window_extent(renderer)
            out.append((t, np.array([e.x0, e.y0, e.x1, e.y1])))
    return out


def _tick_labels(axis: Axis) -> Iterable[Text]:
    # the ticks actually drawn (inside the view)
    ticks = axis._update_ticks()  # type: ignore[attr-defined]
    return [lab for tick in ticks for lab in (tick.label1, tick.label2) if lab.get_visible()]


def _is_reference_line(ax: Axes, line: Line2D) -> bool:
    tr = line.get_transform()
    return tr in (ax.get_xaxis_transform(), ax.get_yaxis_transform())


def _drawn_colours(fig: Figure) -> tuple[list[set[str]], set[str]]:
    """Per axes, the visible face/edge colours of its data artists (hex, upper case, alpha
    ignored, fully transparent dropped); and every colourmap drawn."""
    out: list[set[str]] = []
    maps: set[str] = set()
    for ax in fig.axes:
        rgba: list[np.ndarray] = []
        for p in ax.patches:
            if p.get_visible():
                rgba.append(to_rgba_array(p.get_facecolor()))
                if p.get_linewidth() > 0:
                    rgba.append(to_rgba_array(p.get_edgecolor()))
        for ln in ax.lines:
            if ln.get_visible() and not _is_reference_line(ax, ln):
                rgba += [to_rgba_array(ln.get_color())]
                if ln.get_marker() not in (None, "None", "", " "):
                    rgba += [to_rgba_array(ln.get_markerfacecolor()),
                             to_rgba_array(ln.get_markeredgecolor())]  # fmt: skip
        for c in ax.collections:
            if not c.get_visible():
                continue
            if c.get_array() is not None:
                maps.add(c.get_cmap().name)
            else:
                rgba.append(to_rgba_array(c.get_facecolor()))
            if np.any(np.asarray(c.get_linewidth()) > 0):
                rgba.append(to_rgba_array(c.get_edgecolor()))
        maps |= {im.get_cmap().name for im in ax.images if isinstance(im, ScalarMappable)}
        cols = np.concatenate(rgba) if rgba else np.empty((0, 4))
        out.append({to_hex((float(c[0]), float(c[1]), float(c[2]))).upper()
                    for c in cols if c[3] > 0})  # fmt: skip
    return out, maps
