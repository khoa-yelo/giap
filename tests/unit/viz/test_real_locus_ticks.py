# ruff: noqa: RUF001  (thin / no-break spaces and U+2212 are what tick formatters emit)
"""K-12 plot survey, batch 1 item 6: locus (V03 / V13) x-axis ticks at single-column width.

The "genes" axes' position ticks ("100,000", "102,000", ...) collided at 3.5 in. Required: at most
5 labelled major ticks inside the view, kb formatting, and no two tick labels overlapping in the
rendered (Agg) boxes. On a circular window the labels still read the wrapped position.

Interface assumed: kb formatting means either each label ends in "kb" or the x label says
"(kb)"; a label's number times 1,000 is the position in bp (thousands separators, thin spaces and
U+2212 allowed).
"""

from __future__ import annotations

import itertools
from pathlib import Path

import pytest
from helpers import viz_locus as vl
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from giap.config import GiapConfig, resolve_config
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.base import Instance, Skip

MAX_TICKS = 5
#: the helper's chr1 / lin1 plus a K-12-sized chromosome (7-digit positions)
REPLICONS = (
    *vl.REPLICONS,
    {"id": "big", "kind": "chromosome", "length": 4_641_652, "topology": "circular"},
)


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


def _loci(tmp: Path, cfg: GiapConfig, contig: str, start: int, end: int) -> Instance:
    import giap.viz.locus  # noqa: F401  (registers V03)

    out = vl.build_out(
        tmp,
        replicons=REPLICONS,
        genes=[vl.gene("g1", contig, start + 100, start + 1_000, log2=20.0)],
        regions_stats=[vl.region_stat("prophage_1", "prophage", contig, start, end, 0.003)],
    )
    result = base.registered()["V03"].fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    (inst,) = [i for i in result if i.suffix == "prophage_1"]
    return inst


def _genes_ax(fig: Figure) -> Axes:
    (ax,) = [a for a in fig.axes if a.get_label() == "genes"]
    return ax


def _labels(fig: Figure) -> list[tuple[str, float, tuple[float, float, float, float]]]:
    """(text, tick location, Agg box) of the labelled major x ticks inside the view."""
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    renderer = canvas.get_renderer()
    ax = _genes_ax(fig)
    lo, hi = sorted(ax.get_xlim())
    out = []
    for tick in ax.xaxis._update_ticks():
        loc = float(tick.get_loc())
        if not lo <= loc <= hi:
            continue
        for lab in (tick.label1, tick.label2):
            if lab.get_visible() and lab.get_text().strip():
                e = lab.get_window_extent(renderer)
                out.append((lab.get_text(), loc, (e.x0, e.y0, e.x1, e.y1)))
    return out


def _kb(label: str) -> float:
    s = label.strip().lower().removesuffix("kb").strip()
    for ch in (",", " ", " ", " ", " "):
        s = s.replace(ch, "")
    return float(s.replace("−", "-")) * 1_000


def _is_kb(fig: Figure, labels: list[str]) -> bool:
    xlabel = _genes_ax(fig).get_xlabel().lower()
    return "(kb)" in xlabel or all(t.strip().lower().endswith("kb") for t in labels)


CASES = {
    # circular chr1 (100 kb): a 19 kb window near 40 kb, and one across the origin
    "circular": ("chr1", 40_001, 55_000, 100_000),
    "origin": ("chr1", 95_001, 110_000, 100_000),
    # K-12 scale: 7-digit positions ("4,600,000") at 3.5 in
    "k12": ("big", 4_600_001, 4_615_000, 4_641_652),
    # linear lin1 (30 kb)
    "linear": ("lin1", 8_001, 23_000, None),
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_locus_is_single_column(tmp_path: Path, cfg: GiapConfig, case: str) -> None:
    contig, s, e, _ = CASES[case]
    fig = _loci(tmp_path, cfg, contig, s, e).figure
    assert fig.get_figwidth() == pytest.approx(style.SINGLE_IN)


@pytest.mark.parametrize("case", sorted(CASES))
def test_at_most_five_major_ticks(tmp_path: Path, cfg: GiapConfig, case: str) -> None:
    contig, s, e, _ = CASES[case]
    labels = _labels(_loci(tmp_path, cfg, contig, s, e).figure)
    assert 2 <= len(labels) <= MAX_TICKS, [t for t, _, _ in labels]


@pytest.mark.parametrize("case", sorted(CASES))
def test_tick_labels_do_not_overlap(tmp_path: Path, cfg: GiapConfig, case: str) -> None:
    contig, s, e, _ = CASES[case]
    labels = _labels(_loci(tmp_path, cfg, contig, s, e).figure)
    bad = []
    for (ta, _, a), (tb, _, b) in itertools.combinations(labels, 2):
        w = min(a[2], b[2]) - max(a[0], b[0])
        h = min(a[3], b[3]) - max(a[1], b[1])
        if w > 0.5 and h > 0.5:
            bad.append((ta, tb))
    assert bad == [], bad


@pytest.mark.parametrize("case", sorted(CASES))
def test_ticks_are_in_kb_and_read_the_position(tmp_path: Path, cfg: GiapConfig, case: str) -> None:
    contig, s, e, length = CASES[case]
    fig = _loci(tmp_path, cfg, contig, s, e).figure
    labels = _labels(fig)
    assert _is_kb(fig, [t for t, _, _ in labels]), ([t for t, _, _ in labels],
                                                     _genes_ax(fig).get_xlabel())  # fmt: skip
    for text, loc, _ in labels:
        want = loc % length if length else loc
        assert _kb(text) == pytest.approx(want, abs=0.5), (text, loc)
        if length:
            assert 0 <= _kb(text) < length, text
