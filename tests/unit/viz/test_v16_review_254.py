# ruff: noqa: RUF001, RUF002, RUF003  (ρ is the spec's own character)
"""V16 review fixes for PR #254 (M9-23 lead decision 5, B1-B3 and the plot half of B5).

- B1: the title's "detected in all k runs" count comes from the full intersection table, not the
  30 drawn bars (the all-runs set can rank below 30 and be folded into ``other``).
- B2: at 12 runs no two texts of panel (b) overlap and every cell text fits inside its cell
  (§14.4 "no overlapping labels"; the M9-24 checker is not on main yet, so the window extents are
  measured here after an Agg draw, as ``helpers.viz_landscape.text_overlaps`` does).
- B3: 12 runs x 400 shared genes (26,400 points, binned) draw at most 20,000 elements in the SVG,
  and the 20,001-point two-run case is still binned in the SVG too.
- B5: a negative ρ is hatched and clipped to the colour of ρ = 0; ρ cells sit in the lower
  triangle and the scatters and identity lines in the upper; every pair has its identity line.

Geometry is read in display coordinates against the diagonal run labels (the ``Text`` artists of
axes ``"b"`` whose text is a run label), so the tests do not depend on the axis direction.
"""

from __future__ import annotations

import io
import itertools
import math
import re

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import pytest
from helpers import viz_compare as vc
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib import colormaps
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import Collection
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text

from giap.viz import base, style
from giap.viz.base import Instance

MAX_ELEMENTS = 20_000
#: SVG elements that put ink on the page (matplotlib writes markers as one <use> each)
_DRAWABLE = re.compile(r"<(path|use|rect|circle|ellipse|line|polyline|polygon|text|image)\b")


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


def _plot(genes: dict[str, pd.DataFrame], abundance: pd.DataFrame | None = None,
          B: int = 20) -> Instance:  # fmt: skip
    from giap.viz.compare import compare_plot

    d, a = vc.inputs(genes, B=B)
    return compare_plot(d, abundance if abundance is not None else a, genes)


def _axes(fig: Figure, label: str) -> Axes:
    found = [ax for ax in fig.axes if ax.get_label() == label]
    assert len(found) == 1, [ax.get_label() for ax in fig.axes]
    return found[0]


def _draw(fig: Figure) -> object:
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    return canvas.get_renderer()


def _svg(fig: Figure) -> str:
    buf = io.StringIO()
    fig.savefig(buf, format="svg")
    return buf.getvalue()


def _drawable(svg: str) -> int:
    return len(_DRAWABLE.findall(svg))


def _diag(ax: Axes, labels: list[str]) -> dict[str, tuple[float, float]]:
    """Display centre of each diagonal run label of panel (b)."""
    out: dict[str, tuple[float, float]] = {}
    for t in ax.texts:
        if t.get_text() in labels:
            x, y = ax.transData.transform(t.get_position())
            out[t.get_text()] = (float(x), float(y))
    assert sorted(out) == sorted(labels), out
    return out


def _cell_px(ax: Axes) -> float:
    (x0, _), (x1, _) = ax.transData.transform([(0, 0), (1, 0)])
    return abs(float(x1 - x0))


def _column(x: float, diag: dict[str, tuple[float, float]]) -> tuple[str, float]:
    """The run whose diagonal label is nearest in x, and that label's display y."""
    lab = min(diag, key=lambda k: abs(diag[k][0] - x))
    return lab, diag[lab][1]


def _rho_rects(ax: Axes) -> list[Rectangle]:
    return [p for p in ax.patches if isinstance(p, Rectangle)]


def _centre(ax: Axes, r: Rectangle) -> tuple[float, float]:
    x, y = r.get_xy()
    cx, cy = ax.transData.transform((x + r.get_width() / 2, y + r.get_height() / 2))
    return float(cx), float(cy)


def _three_full() -> dict[str, pd.DataFrame]:
    return vc.detect_runs({"A1": set(range(8)), "A2": set(range(8)), "A3": set(range(8))}, 8)


# --- B1: the all-runs count comes from the full intersection table ----------------------------


def _all_runs_below_30(k: int = 6, per: int = 5, n_all: int = 1) -> dict[str, pd.DataFrame]:
    """41 one-, two- and three-run intersections of ``per`` genes each and ``n_all`` genes in
    all ``k`` runs: the all-runs set ranks below 30 and is folded into ``other``."""
    labels = [f"R{i}" for i in range(k)]
    sets: dict[str, set[int]] = {lab: set() for lab in labels}
    subsets = [s for m in (1, 2, 3) for s in itertools.combinations(range(k), m)]
    g = 0
    for s in subsets:
        for _ in range(per):
            for b in s:
                sets[labels[b]].add(g)
            g += 1
    for _ in range(n_all):
        for lab in labels:
            sets[lab].add(g)
        g += 1
    return vc.detect_runs(sets, g)


def test_b1_all_runs_count_when_the_all_runs_set_is_folded_into_other() -> None:
    genes = _all_runs_below_30()
    d, _ = vc.inputs(genes, B=5)
    inter = vc.intersection_rows(d)
    assert len(inter[(inter["intersection"] != "none") & (inter["n_genes"] > 0)]) == 42
    inst = _plot(genes, B=5)
    drawn = inst.data[(inst.data["panel"] == "a")]
    assert "R0&R1&R2&R3&R4&R5" not in set(drawn["intersection"].astype(str))  # in "other"
    assert inst.title.startswith("1 gene detected in all 6 runs; "), inst.title


@settings(max_examples=8, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(st.integers(0, 3))
def test_b1_all_runs_count_property(n_all: int) -> None:
    """Whatever the all-runs count, when it ranks below the 30 largest the title still states it."""
    inst = _plot(_all_runs_below_30(n_all=n_all), B=5)
    assert inst.title.startswith(f"{base.plural(n_all, 'gene')} detected in all 6 runs; "), (
        inst.title
    )


# --- B2: no cell-text overlap at 12 runs ------------------------------------------------------


def _overlaps(boxes: list[tuple[str, object]]) -> list[tuple[str, str]]:
    hits = []
    for (s1, b1), (s2, b2) in itertools.combinations(boxes, 2):
        w = min(b1.x1, b2.x1) - max(b1.x0, b2.x0)  # type: ignore[attr-defined]
        h = min(b1.y1, b2.y1) - max(b1.y0, b2.y0)  # type: ignore[attr-defined]
        if w > 1.0 and h > 1.0:
            hits.append((s1, s2))
    return hits


@pytest.fixture(scope="module")
def twelve() -> Instance:
    style.apply()
    genes = vc.many_intersections(40, 12)
    d, a = vc.inputs(genes, B=5)
    # one pair too few shared (ρ NA) and one negative ρ, the two widest cell texts
    from giap.viz.compare import compare_plot

    a = vc.set_rho(a, "R01", "R02", None, None, None, flags="too_few_shared")
    a = vc.set_rho(a, "R03", "R04", -0.25, -0.61, -0.013)
    return compare_plot(d, a, genes)


def test_b2_twelve_runs_no_two_panel_b_texts_overlap(twelve: Instance) -> None:
    renderer = _draw(twelve.figure)
    ax = _axes(twelve.figure, "b")
    texts = [t for t in ax.texts if t.get_visible() and t.get_text().strip()]
    boxes = [(t.get_text(), t.get_window_extent(renderer)) for t in texts]
    hits = _overlaps(boxes)
    assert hits == [], hits[:5]


def test_b2_twelve_runs_every_cell_text_fits_its_cell(twelve: Instance) -> None:
    renderer = _draw(twelve.figure)
    ax = _axes(twelve.figure, "b")
    bad = []
    for t in ax.texts:
        if not t.get_text().strip():
            continue
        x, y = t.get_position()
        col, row = math.floor(x), math.floor(y)
        (cx0, cy0), (cx1, cy1) = ax.transData.transform([(col, row), (col + 1, row + 1)])
        lo_x, hi_x = sorted((cx0, cx1))
        lo_y, hi_y = sorted((cy0, cy1))
        bb = t.get_window_extent(renderer)
        if bb.x0 < lo_x - 0.5 or bb.x1 > hi_x + 0.5 or bb.y0 < lo_y - 0.5 or bb.y1 > hi_y + 0.5:
            bad.append((t.get_text(), round(bb.width, 1), round(hi_x - lo_x, 1)))
    assert bad == [], bad[:5]


def test_b2_whole_figure_has_no_text_overlap_at_twelve_runs(twelve: Instance) -> None:
    from helpers.viz_landscape import text_overlaps

    assert text_overlaps(twelve.figure) == []


def test_b2_too_few_shared_reason_stays_in_the_figure(twelve: Instance) -> None:
    """Lead decision 5: the cell may shorten to "ρ NA" but the reason stays (a footnote)."""
    texts = " ".join(t.get_text() for t in twelve.figure.findobj(Text))
    assert "too few shared" in texts


# --- B3: at most 20,000 drawn elements when binned --------------------------------------------


@pytest.fixture(scope="module")
def binned_twelve() -> Instance:
    style.apply()
    n = 400
    genes = vc.detect_runs({f"R{i + 1:02d}": set(range(n)) for i in range(12)}, n)
    return _plot(genes, B=2)


def test_b3_twelve_by_400_is_binned(binned_twelve: Instance) -> None:
    pts = binned_twelve.data[binned_twelve.data["locus_tag"].notna()]
    assert len(pts) == 66 * 400 > MAX_ELEMENTS  # 26,400 points: binning applies
    ax = _axes(binned_twelve.figure, "b")
    marks = sum(len(c.get_offsets()) for c in ax.collections if c.get_offsets() is not None)
    assert marks < len(pts), "26,400 points drawn one by one: not binned"


def test_b3_twelve_by_400_svg_has_at_most_20000_elements(binned_twelve: Instance) -> None:
    svg = _svg(binned_twelve.figure)
    n = _drawable(svg)
    assert n <= MAX_ELEMENTS, f"{n} drawable SVG elements ({len(svg)} bytes)"


def test_b3_twelve_by_400_artists_draw_at_most_20000_cells(binned_twelve: Instance) -> None:
    """QuadMesh quads count too: every quad (masked or not) is a <path> in the SVG."""
    from matplotlib.collections import QuadMesh

    total = 0
    for c in binned_twelve.figure.findobj(Collection):
        if isinstance(c, QuadMesh):
            total += int(np.asarray(c.get_coordinates()).shape[0] - 1) * int(
                np.asarray(c.get_coordinates()).shape[1] - 1
            )
        else:
            total += max(len(c.get_offsets()), len(c.get_paths()))
    assert total <= MAX_ELEMENTS, total


def test_b3_two_runs_above_20000_points_svg_is_capped() -> None:
    n = 20_001
    genes = vc.detect_runs({"A1": set(range(n)), "A2": set(range(n))}, n)
    inst = _plot(genes, B=2)
    assert _drawable(_svg(inst.figure)) <= MAX_ELEMENTS


def test_b3_binned_cells_keep_the_identity_line(binned_twelve: Instance) -> None:
    ax = _axes(binned_twelve.figure, "b")
    assert len(_identity_lines(ax)) >= 66


# --- B5: negative ρ hatched and clipped -------------------------------------------------------


def _hex(c: object) -> str:
    return mcolors.to_hex(c).upper()  # type: ignore[arg-type]


def _cell_rect(ax: Axes, labels: list[str], ra: str, rb: str) -> Rectangle:
    """The ρ rectangle at the crossing of run ``ra``'s and run ``rb``'s diagonal labels."""
    diag = _diag(ax, labels)
    half = _cell_px(ax) / 2
    want = {(diag[ra][0], diag[rb][1]), (diag[rb][0], diag[ra][1])}
    found = [r for r in _rho_rects(ax)
             if any(abs(_centre(ax, r)[0] - x) < half and abs(_centre(ax, r)[1] - y) < half
                    for x, y in want)]  # fmt: skip
    assert len(found) == 1, (ra, rb, len(found))
    return found[0]


def _with_rho(cells: list[tuple[str, str, float, float, float]]) -> Instance:
    genes = _three_full()
    _, a = vc.inputs(genes)
    for ra, rb, rho, lo, hi in cells:
        a = vc.set_rho(a, ra, rb, rho, lo, hi)
    return _plot(genes, a)


def test_b5_negative_rho_is_hatched() -> None:
    inst = _with_rho([("A1", "A2", -0.5, -0.8, -0.1), ("A1", "A3", 0.6, 0.4, 0.8),
                      ("A2", "A3", 0.9, 0.8, 0.95)])  # fmt: skip
    ax = _axes(inst.figure, "b")
    labels = ["A1", "A2", "A3"]
    assert _cell_rect(ax, labels, "A1", "A2").get_hatch()
    assert not _cell_rect(ax, labels, "A1", "A3").get_hatch()
    assert not _cell_rect(ax, labels, "A2", "A3").get_hatch()


@pytest.mark.parametrize("rho", [-0.5, -0.05, -1.0])
def test_b5_negative_rho_is_clipped_to_the_bottom_colour(rho: float) -> None:
    inst = _with_rho([("A1", "A2", rho, rho - 0.1 if rho > -1 else -1.0, rho / 2)])
    ax = _axes(inst.figure, "b")
    face = _cell_rect(ax, ["A1", "A2", "A3"], "A1", "A2").get_facecolor()
    assert _hex(face) == _hex(colormaps["cividis"](0.0)), _hex(face)


def test_b5_positive_rho_is_coloured_by_its_value() -> None:
    inst = _with_rho([("A1", "A2", 0.5, 0.2, 0.7)])
    ax = _axes(inst.figure, "b")
    face = _cell_rect(ax, ["A1", "A2", "A3"], "A1", "A2").get_facecolor()
    assert _hex(face) == _hex(colormaps["cividis"](0.5))


# --- B5: lower triangle ρ, upper triangle scatters, identity line -----------------------------


def _identity_lines(ax: Axes) -> list[Line2D]:
    """Two-point lines rising at 45 degrees in display space (the identity, y = x)."""
    out = []
    for ln in ax.lines:
        xy = ax.transData.transform(np.column_stack([ln.get_xdata(), ln.get_ydata()]))
        if len(xy) != 2:
            continue
        dx, dy = xy[1] - xy[0]
        if dx == 0:
            continue
        if dx < 0:
            dx, dy = -dx, -dy
        if dy > 0 and math.isclose(dx, dy, rel_tol=0.02):
            out.append(ln)
    return out


@pytest.fixture(scope="module")
def four() -> Instance:
    style.apply()
    labels = ["A1", "A2", "A3", "A4"]
    return _plot(vc.detect_runs({lab: set(range(10)) for lab in labels}, 10))


LABELS4 = ["A1", "A2", "A3", "A4"]


def test_b5_rho_cells_are_in_the_lower_triangle(four: Instance) -> None:
    ax = _axes(four.figure, "b")
    diag = _diag(ax, LABELS4)
    rects = _rho_rects(ax)
    assert len(rects) == 6
    half = _cell_px(ax) / 2
    for r in rects:
        cx, cy = _centre(ax, r)
        _, ly = _column(cx, diag)
        assert cy < ly - half, ("ρ cell above the diagonal", r.get_xy())


def test_b5_each_pair_has_one_lower_cell(four: Instance) -> None:
    ax = _axes(four.figure, "b")
    for ra, rb in itertools.combinations(LABELS4, 2):
        _cell_rect(ax, LABELS4, ra, rb)


def test_b5_scatters_are_in_the_upper_triangle(four: Instance) -> None:
    ax = _axes(four.figure, "b")
    diag = _diag(ax, LABELS4)
    half = _cell_px(ax) / 2
    n = 0
    for c in ax.collections:
        offs = np.asarray(c.get_offsets())
        if not len(offs):
            continue
        for x, y in c.get_offset_transform().transform(offs):
            _, ly = _column(float(x), diag)
            assert y > ly + 0.5 * half, ("scatter point below the diagonal", x, y)
            n += 1
    assert n == 6 * 10


def test_b5_every_pair_has_an_identity_line_in_its_upper_cell(four: Instance) -> None:
    ax = _axes(four.figure, "b")
    diag = _diag(ax, LABELS4)
    half = _cell_px(ax) / 2
    lines = _identity_lines(ax)
    cells: set[tuple[str, str]] = set()
    for ln in lines:
        xy = ax.transData.transform(np.column_stack([ln.get_xdata(), ln.get_ydata()]))
        cx, cy = xy.mean(axis=0)
        col, ly = _column(float(cx), diag)
        assert cy > ly + half, ("identity line below the diagonal", cx, cy)
        row = min(diag, key=lambda k: abs(diag[k][1] - cy))
        cells.add((row, col))
    want = {(ra, rb) for ra, rb in itertools.product(LABELS4, LABELS4)
            if LABELS4.index(ra) < LABELS4.index(rb)}  # fmt: skip
    assert cells == want, sorted(cells)
