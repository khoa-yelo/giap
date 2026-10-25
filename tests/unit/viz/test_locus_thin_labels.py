"""V3/V13 gene labels do not overlap (#323 review S1; ``giap.viz.locus._thin_labels``).

Gene labels are the texts on the ``"genes"`` axes whose string is a gene name of the window. On
a dense window (short adjacent genes, long names) some labels are hidden; the visible ones have
pairwise disjoint laid-out boxes (Agg renderer, after layout). On a sparse window every label is
kept. Two renders hide the same labels.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from helpers import viz_locus as vl
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.text import Text

from giap.config import resolve_config
from giap.io.run_view import RunView
from giap.viz import style
from giap.viz.locus import Locus, render_locus

SCALE = (18.0, 24.0)


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


def _genes(n: int, length: int, gap: int, name: str) -> list[dict[str, object]]:
    rows = []
    for i in range(n):
        start = 10_001 + i * (length + gap)
        rows.append(vl.gene(f"G{i:03d}", "chr1", start, start + length - 1, log2=20.0,
                            name=f"{name}{i:02d}"))  # fmt: skip
    return rows


def _render(tmp: Path, genes: list[dict[str, object]]) -> tuple[Figure, list[str]]:
    out = vl.build_out(tmp, genes=genes)
    view = RunView.from_outdir(out)
    first, last = int(genes[0]["start"]) - 1, int(genes[-1]["end"])  # type: ignore[call-overload]
    lo = Locus("L1", "island", "chr1", first, last, None, "test locus")
    inst = render_locus(lo, view, resolve_config(), SCALE)
    return inst.figure, [str(g["gene"]) for g in genes]


def _gene_axes(fig: Figure) -> Axes:
    (ax,) = [a for a in fig.axes if a.get_label() == "genes"]
    return ax


def _labels(fig: Figure, names: list[str]) -> list[Text]:
    want = set(names)
    return [t for t in _gene_axes(fig).texts if t.get_text() in want]


def _visible_boxes(fig: Figure, names: list[str]) -> list[tuple[str, object]]:
    canvas = FigureCanvasAgg(fig)  # the public Agg renderer, after a full layout and draw
    canvas.draw()
    renderer = canvas.get_renderer()
    return [(t.get_text(), t.get_window_extent(renderer))
            for t in _labels(fig, names) if t.get_visible()]  # fmt: skip


def _overlapping(boxes: list[tuple[str, object]]) -> list[tuple[str, str]]:
    bad = []
    for i, (a, ba) in enumerate(boxes):
        for b, bb in boxes[i + 1 :]:
            if ba.overlaps(bb):  # type: ignore[attr-defined]
                bad.append((a, b))
    return bad


def test_dense_labels_are_thinned_and_the_kept_ones_do_not_overlap(tmp_path: Path) -> None:
    genes = _genes(10, 250, 10, "longgenename_")
    fig, names = _render(tmp_path, genes)
    labels = _labels(fig, names)
    assert len(labels) >= 2, [t.get_text() for t in _gene_axes(fig).texts]
    hidden = [t.get_text() for t in labels if not t.get_visible()]
    assert hidden, "overlapping gene labels are all still visible"
    boxes = _visible_boxes(fig, names)
    assert boxes, "every gene label is hidden"
    assert not _overlapping(boxes), _overlapping(boxes)


def test_sparse_labels_are_all_kept(tmp_path: Path) -> None:
    genes = _genes(3, 900, 3000, "g")
    fig, names = _render(tmp_path, genes)
    labels = _labels(fig, names)
    assert sorted(t.get_text() for t in labels) == sorted(names)
    assert all(t.get_visible() for t in labels), [
        t.get_text() for t in labels if not t.get_visible()
    ]
    assert not _overlapping(_visible_boxes(fig, names))


def test_thinning_is_deterministic(tmp_path: Path) -> None:
    genes = _genes(10, 250, 10, "longgenename_")
    runs = []
    for k in range(2):
        fig, names = _render(tmp_path / f"r{k}", genes)
        runs.append([(t.get_text(), t.get_visible()) for t in _labels(fig, names)])
    assert runs[0] == runs[1]
    assert any(not v for _, v in runs[0])
