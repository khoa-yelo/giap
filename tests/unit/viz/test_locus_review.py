"""V3 ``locus_<id>`` review round (PR #209 review; M9-10 lead decisions for the round, B1-B5).

- B1: an operon flag whose two genes are in different TUs and straddle the origin of a circular
  contig spans the shorter arc (``b`` unwrapped relative to ``a``), not the whole replicon.
- B2: the title wraps and stays inside the figure (window extent measured with
  ``FigureCanvasAgg`` after a draw; ``FigureCanvasBase`` gives false passes); every clause, with its
  q value, is still in the title text (whitespace-collapsed).
- B3: a colour key: a colorbar on an Axes labelled ``"log2 abundance"`` (the label chosen here)
  whose limits are the run-wide scale and whose cmap is cividis, plus legend swatches for
  undetected (``#BBBBBB`` fill), imputed (hatched) and single_peptide_unestimated (``#E69F00``
  outline); still >= 6 pt and no overlaps between the title, key and tracks.
- B4: strand drawing: a ``+`` arrow's head points right, a ``-`` arrow's left; the promoter bent
  arrow's stem sits at the low coordinate on ``+`` and the high coordinate on ``-`` (lead decision
  for this round) and points downstream; ``-`` strand marks are below the gene-track centre, ``+``
  above.
- B5: imputed genes keep their data row with ``value`` NA and ``category = imputed``.

Plus the reviewer's test-gap survivors: ``_pct`` half-up rounding, the missing-``gene_b`` drop,
the different-TU span when ``a`` ends last, ``run_scale`` with equal percentiles and with
imputed-but-detected rows, and the origin-line position to < 1 bp.
"""

from __future__ import annotations

import itertools
from pathlib import Path
from typing import Any

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import pytest
from helpers import viz_locus as vl
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.patches import Patch
from matplotlib.text import Text
from matplotlib.transforms import Bbox

from giap.config import GiapConfig, resolve_config
from giap.io.plot_files import write_instance
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.base import Instance, Skip
from giap.viz.locus import render_locus, run_scale, select_loci

SCALE = (18.0, 24.0)
KEY_LABEL = "log2 abundance"
PROMOTER = "#009E73"
TERMINATOR = "#D55E00"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


# --- helpers ------------------------------------------------------------------------------------


def _render(
    out: Path, cfg: GiapConfig, locus_id: str, scale: tuple[float, float] = SCALE
) -> Instance:
    view = RunView.from_outdir(out)
    (lo,) = [x for x in select_loci(view, cfg) if x.locus_id == locus_id]
    inst = render_locus(lo, view, cfg, scale)
    assert isinstance(inst, Instance)
    return inst


def _run_plot(out: Path, cfg: GiapConfig) -> list[Instance]:
    import giap.viz.locus  # noqa: F401  (registers V03)

    result = base.registered()["V03"].fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip)
    return list(result)


def _axes(fig: Figure, label: str) -> Axes:
    found = [ax for ax in fig.axes if ax.get_label() == label]
    assert len(found) == 1, [ax.get_label() for ax in fig.axes]
    return found[0]


def _agg(fig: Figure) -> Any:
    """Attach an Agg canvas, draw, and return its renderer (real text extents)."""
    canvas = FigureCanvasAgg(fig)
    renderer = canvas.get_renderer()
    fig.draw(renderer)
    return renderer


def _collapse(s: str) -> str:
    return " ".join(s.split())


def _hex(c: Any) -> str | None:
    try:
        return mcolors.to_hex(c).upper()
    except (ValueError, TypeError):
        return None


def _dashed_x(ax: Axes) -> list[float]:
    """Data x of every dashed vertical line on ``ax``."""
    xs = []
    for ln in ax.get_lines():
        x = np.asarray(ln.get_xdata(), dtype=float)
        if ln.get_linestyle() == "--" and len(x) >= 2 and np.allclose(x, x[0]):
            xs.append(float(x[0]))
    return xs


def _data_vertices(ax: Axes, p: Patch) -> np.ndarray:
    """The patch outline in ``ax`` data coordinates."""
    verts = p.get_transform().transform(p.get_path().vertices)
    return ax.transData.inverted().transform(verts)


def _opaque_patches(ax: Axes) -> list[Patch]:
    return [p for p in ax.patches if p.get_facecolor()[3] > 0.5]


def _gene_patch(ax: Axes, s: float, e: float) -> np.ndarray:
    """Vertices of the single opaque patch on ``ax`` spanning ``[s, e)`` in x."""
    width = abs(ax.get_xlim()[1] - ax.get_xlim()[0])
    tol = 0.002 * width
    found = []
    for p in _opaque_patches(ax):
        v = _data_vertices(ax, p)
        if abs(v[:, 0].min() - s) <= tol and abs(v[:, 0].max() - e) <= tol:
            found.append(v)
    assert len(found) == 1, (s, e, len(found))
    return found[0]


def _end_spread(v: np.ndarray, end: str, width: float) -> float:
    """The y spread of the outline's vertices at its ``"min"`` or ``"max"`` x end."""
    x = v[:, 0].min() if end == "min" else v[:, 0].max()
    at = v[np.abs(v[:, 0] - x) <= 1e-6 * width]
    return float(at[:, 1].max() - at[:, 1].min())


# --- B1: operon flag straddling the origin ------------------------------------------------------


def _straddle_out(tmp: Path, *, a_first: bool = True) -> Path:
    ga = vl.gene("ga", "chr1", 99_501, 100_000, tu_id="TUa", log2=20.0)
    gb = vl.gene("gb", "chr1", 1, 500, tu_id="TUb", log2=21.0)
    tus = [
        vl.tu("TUa", "chr1", 99_501, 100_000, "+", ["ga"]),
        vl.tu("TUb", "chr1", 1, 500, "+", ["gb"]),
    ]
    pair = ("ga", "gb") if a_first else ("gb", "ga")
    return vl.build_out(tmp, genes=[ga, gb], tus=tus, flags=[vl.flag("PAIR_1", *pair, 0.01)])


def test_b1_flag_across_the_origin_spans_the_shorter_arc(tmp_path: Path, cfg: GiapConfig) -> None:
    (lo,) = select_loci(RunView.from_outdir(_straddle_out(tmp_path)), cfg)
    assert (lo.start0, lo.end0) == (99_500, 100_500)


@pytest.mark.parametrize("a_first", [True, False], ids=["a_before_origin", "a_after_origin"])
def test_b1_flag_across_the_origin_span_is_1_kb(
    tmp_path: Path, cfg: GiapConfig, a_first: bool
) -> None:
    (lo,) = select_loci(RunView.from_outdir(_straddle_out(tmp_path, a_first=a_first)), cfg)
    assert lo.end0 - lo.start0 == 1_000
    assert lo.start0 % 100_000 in (0, 99_500)


@pytest.mark.parametrize("a_first", [True, False], ids=["a_before_origin", "a_after_origin"])
def test_b1_flag_across_the_origin_window_is_single_with_one_origin_line(
    tmp_path: Path, cfg: GiapConfig, a_first: bool
) -> None:
    inst = _render(_straddle_out(tmp_path, a_first=a_first), cfg, "PAIR_1")
    ax = _axes(inst.figure, "genes")
    lo, hi = ax.get_xlim()
    assert hi - lo == pytest.approx(1_000 + 2 * 2_000, abs=1)
    assert round(float(inst.figure.get_size_inches()[0]), 6) == style.SINGLE_IN
    xs = _dashed_x(ax)
    assert len(xs) == 1
    assert xs[0] % 100_000 == pytest.approx(0, abs=0.5)
    assert lo < xs[0] < hi
    genes = inst.data[inst.data["track"] == "gene"]
    assert set(genes["id"]) == {"ga", "gb"}
    assert (genes["end"] - genes["start"]).tolist() == [500, 500]
    assert genes["start"].max() - genes["start"].min() == 500  # adjacent in the window frame


def test_b1_index_row_is_not_the_whole_replicon(tmp_path: Path, cfg: GiapConfig) -> None:
    result = _run_plot(_straddle_out(tmp_path), cfg)
    (index,) = [i for i in result if i.suffix is None]
    row = index.data.set_index("locus_id").loc["PAIR_1"]
    assert (int(row["start"]), int(row["end"])) == (99_500, 100_500)
    write_instance(index, "V03", "locus", tmp_path / "plots", ["svg"], 100)
    tsv = (tmp_path / "plots" / "data" / "V03_locus.tsv").read_text(encoding="utf-8")
    line = next(ln for ln in tsv.splitlines() if "\tPAIR_1\t" in ln)
    assert "\tchr1\t99501\t100500\t" in line
    assert "\t1\t100000\t" not in line


def test_different_tu_span_when_a_ends_last(tmp_path: Path, cfg: GiapConfig) -> None:
    genes = [
        vl.gene("ga", "chr1", 12_001, 13_000, tu_id="TUa", log2=20.0),
        vl.gene("gb", "chr1", 10_001, 11_000, tu_id="TUb", log2=21.0),
    ]
    tus = [
        vl.tu("TUa", "chr1", 12_001, 13_000, "+", ["ga"]),
        vl.tu("TUb", "chr1", 10_001, 11_000, "+", ["gb"]),
    ]
    out = vl.build_out(tmp_path, genes=genes, tus=tus, flags=[vl.flag("PAIR_1", "ga", "gb", 0.01)])
    (lo,) = select_loci(RunView.from_outdir(out), cfg)
    assert (lo.start0, lo.end0) == (10_000, 13_000)


def test_flag_with_missing_gene_b_is_dropped(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vl.build_out(
        tmp_path,
        genes=[vl.gene("ga", "chr1", 12_001, 13_000, tu_id="TUa", log2=20.0)],
        tus=[vl.tu("TUa", "chr1", 12_001, 13_000, "+", ["ga"])],
        flags=[vl.flag("PAIR_1", "ga", "nosuch", 0.01)],
    )
    assert select_loci(RunView.from_outdir(out), cfg) == []


# --- B2: the title wraps inside the figure ------------------------------------------------------


def _single_flag_out(tmp: Path) -> tuple[Path, list[str]]:
    genes = [
        vl.gene("BSU_12340", "chr1", 40_001, 41_000, tu_id="TU0001", log2=20.0),
        vl.gene("BSU_12350", "chr1", 41_101, 42_000, tu_id="TU0001", log2=22.0),
    ]
    tus = [vl.tu("TU0001", "chr1", 40_001, 42_000, "+", ["BSU_12340", "BSU_12350"])]
    flags = [
        vl.flag("PAIR_1", "BSU_12340", "BSU_12350", 0.0123, flag="coabundant_despite_terminator")
    ]
    clauses = ["coabundant_despite_terminator between BSU_12340 and BSU_12350 (q = 0.0123)"]
    return vl.build_out(tmp, genes=genes, tus=tus, flags=flags), clauses


def _merged_flags_out(tmp: Path) -> tuple[Path, list[str]]:
    genes = [
        vl.gene("g3", "chr1", 40_001, 41_000, tu_id="TU0001", log2=20.0, name="lacZ"),
        vl.gene("g4", "chr1", 41_101, 42_000, tu_id="TU0001", log2=22.0, name="lacY"),
        vl.gene("g4b", "chr1", 42_101, 42_600, tu_id="TU0001", log2=23.0, name="lacA"),
    ]
    tus = [vl.tu("TU0001", "chr1", 40_001, 42_600, "+", ["g3", "g4", "g4b"])]
    flags = [
        vl.flag("PAIR_7", "g3", "g4", 0.02),
        vl.flag("PAIR_8", "g4", "g4b", 0.05, flag="unexplained_discordance"),
        vl.flag("PAIR_9", "g3", "g4b", 0.001, flag="conflict_terminator_in_operon"),
    ]
    clauses = [
        "conflict_terminator_in_operon between g3 and g4b (q = 0.001)",
        "internal_promoter_supported between g3 and g4 (q = 0.02)",
        "unexplained_discordance between g4 and g4b (q = 0.05)",
    ]
    return vl.build_out(tmp, genes=genes, tus=tus, flags=flags), clauses


def _title_text(fig: Figure, clauses: list[str]) -> Text:
    """The one drawn Text carrying every clause (whitespace-collapsed)."""
    found = [
        t
        for t in fig.findobj(Text)
        if t.get_visible() and all(c in _collapse(t.get_text()) for c in clauses)
    ]
    assert len(found) == 1, [t.get_text() for t in fig.findobj(Text) if t.get_text()]
    return found[0]


@pytest.mark.parametrize("builder", [_single_flag_out, _merged_flags_out], ids=["one", "three"])
def test_b2_title_keeps_every_clause(tmp_path: Path, cfg: GiapConfig, builder: Any) -> None:
    out, clauses = builder(tmp_path)
    (lo,) = select_loci(RunView.from_outdir(out), cfg)
    inst = _render(out, cfg, lo.locus_id)
    title = _title_text(inst.figure, clauses)
    for c in clauses:
        assert c in _collapse(title.get_text())


@pytest.mark.parametrize("builder", [_single_flag_out, _merged_flags_out], ids=["one", "three"])
def test_b2_title_stays_inside_the_figure(tmp_path: Path, cfg: GiapConfig, builder: Any) -> None:
    out, clauses = builder(tmp_path)
    (lo,) = select_loci(RunView.from_outdir(out), cfg)
    inst = _render(out, cfg, lo.locus_id)
    fig = inst.figure
    renderer = _agg(fig)
    title = _title_text(fig, clauses)
    bb = title.get_window_extent(renderer)
    fb = fig.bbox
    tol = 0.5  # px
    assert bb.x0 >= fb.x0 - tol, (bb, fb)
    assert bb.x1 <= fb.x1 + tol, (bb, fb)
    assert bb.y0 >= fb.y0 - tol, (bb, fb)
    assert bb.y1 <= fb.y1 + tol, (bb, fb)
    for label in ("pileup", "genes"):
        ab = _axes(fig, label).get_window_extent(renderer)
        assert not _overlaps(bb, ab), (label, bb, ab)


def test_b2_merged_title_is_wrapped(tmp_path: Path, cfg: GiapConfig) -> None:
    out, clauses = _merged_flags_out(tmp_path)
    (lo,) = select_loci(RunView.from_outdir(out), cfg)
    inst = _render(out, cfg, lo.locus_id)
    assert _title_text(inst.figure, clauses).get_text().count("\n") >= 2


# --- B3: colour key -----------------------------------------------------------------------------


def _key_axes(fig: Figure) -> Axes:
    found = [ax for ax in fig.axes if KEY_LABEL in (ax.get_xlabel(), ax.get_ylabel())]
    assert len(found) == 1, [(ax.get_label(), ax.get_xlabel(), ax.get_ylabel()) for ax in fig.axes]
    return found[0]


def _key_scale(cax: Axes) -> tuple[float, float]:
    """The colour key's value range: the long axis of the colorbar Axes."""
    x, y = sorted(cax.get_xlim()), sorted(cax.get_ylim())
    horizontal = cax.get_xlabel() == KEY_LABEL
    lo, hi = x if horizontal else y
    return float(lo), float(hi)


def _key_mappables(cax: Axes) -> list[Any]:
    return [a for a in [*cax.collections, *cax.images] if hasattr(a, "get_cmap")]


def _mixed_out(tmp: Path) -> Path:
    genes = [
        vl.gene("G1", "chr1", 10_101, 10_600, log2=21.0),
        vl.gene("G2", "chr1", 10_701, 11_200, strand="-", log2=22.0),
    ]
    return vl.build_out(
        tmp,
        genes=genes,
        regions_stats=[vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 12_000, 0.01)],
    )


def test_b3_colorbar_is_cividis_over_the_given_scale(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_mixed_out(tmp_path), cfg, "prophage_1", (18.0, 24.0))
    cax = _key_axes(inst.figure)
    assert _key_scale(cax) == pytest.approx((18.0, 24.0))
    maps = _key_mappables(cax)
    assert maps
    for m in maps:
        assert m.get_cmap().name == "cividis"
        assert (m.norm.vmin, m.norm.vmax) == pytest.approx((18.0, 24.0))


def _scale_run(tmp: Path) -> Path:
    """101 measured genes, log2 = 18 + 0.06 i; run scale (18.12, 23.88)."""
    genes = [
        vl.gene(f"m{i:03d}", "chr1", 900 * i + 1, 900 * i + 600, log2=18.0 + 0.06 * i)
        for i in range(101)
    ]
    genes += [
        vl.gene("impdet", "chr1", 95_001, 95_600, log2=40.0, imputed=True, status="detected"),
        vl.gene("imp_lo", "chr1", 96_001, 96_600, log2=5.0, imputed=True, status="not_detected"),
    ]
    return vl.build_out(
        tmp,
        genes=genes,
        regions_stats=[
            vl.region_stat("prophage_1", "prophage", "chr1", 88_201, 88_800, 0.01),
            vl.region_stat("prophage_2", "prophage", "chr1", 1_801, 2_400, 0.02),
        ],
    )


def test_b3_every_locus_key_shows_the_run_scale(tmp_path: Path, cfg: GiapConfig) -> None:
    out = _scale_run(tmp_path)
    expected = run_scale(RunView.from_outdir(out).table("genes_tsv"))
    assert expected == pytest.approx((18.12, 23.88))
    shown = [i for i in _run_plot(out, cfg) if i.suffix is not None]
    assert {i.suffix for i in shown} == {"prophage_1", "prophage_2"}
    for inst in shown:
        cax = _key_axes(inst.figure)
        assert _key_scale(cax) == pytest.approx(expected), inst.suffix
        assert all(m.get_cmap().name == "cividis" for m in _key_mappables(cax))


def _legend_handles(fig: Figure) -> list[Any]:
    legends = [lg for lg in fig.findobj(Legend)]
    assert legends, "no legend in the figure"
    return [h for lg in legends for h in lg.legend_handles]


def _legend_texts(fig: Figure) -> list[str]:
    return [t.get_text().lower() for lg in fig.findobj(Legend) for t in lg.get_texts()]


def test_b3_legend_swatches_for_undetected_imputed_unestimated(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    fig = _render(_mixed_out(tmp_path), cfg, "prophage_1").figure  # no such genes: key anyway
    handles = [h for h in _legend_handles(fig) if isinstance(h, Patch)]
    assert any(_hex(h.get_facecolor()) == "#BBBBBB" and not h.get_hatch() for h in handles)
    assert any(h.get_hatch() for h in handles)
    assert any(_hex(h.get_edgecolor()) == "#E69F00" for h in handles)
    texts = _legend_texts(fig)
    assert any("undetected" in t for t in texts), texts
    assert any("imputed" in t for t in texts), texts
    assert any("single" in t or "unestimated" in t for t in texts), texts


def test_b3_every_text_is_at_least_6_pt_after_drawing(tmp_path: Path, cfg: GiapConfig) -> None:
    fig = _render(_mixed_out(tmp_path), cfg, "prophage_1").figure
    _agg(fig)
    sizes = [t.get_fontsize() for t in fig.findobj(Text) if t.get_visible() and t.get_text()]
    assert sizes
    assert min(sizes) >= 6


def _overlaps(a: Bbox, b: Bbox, tol: float = 0.5) -> bool:
    w = min(a.x1, b.x1) - max(a.x0, b.x0)
    h = min(a.y1, b.y1) - max(a.y0, b.y0)
    return w > tol and h > tol


@pytest.mark.parametrize("which", ["mixed", "merged"])
def test_b3_key_title_and_tracks_do_not_overlap(
    tmp_path: Path, cfg: GiapConfig, which: str
) -> None:
    if which == "mixed":
        fig = _render(_mixed_out(tmp_path), cfg, "prophage_1").figure
        clauses = ["prophage_1:"]
    else:
        out, clauses = _merged_flags_out(tmp_path)
        (lo,) = select_loci(RunView.from_outdir(out), cfg)
        fig = _render(out, cfg, lo.locus_id).figure
    renderer = _agg(fig)
    boxes: dict[str, Bbox] = {
        "pileup": _axes(fig, "pileup").get_window_extent(renderer),
        "genes": _axes(fig, "genes").get_window_extent(renderer),
        "colorbar": _key_axes(fig).get_window_extent(renderer),
        "title": _title_text(fig, clauses).get_window_extent(renderer),
    }
    for i, lg in enumerate(fig.findobj(Legend)):
        boxes[f"legend{i}"] = lg.get_window_extent(renderer)
    fb = fig.bbox
    names = sorted(boxes)
    for n in names:
        b = boxes[n]
        assert b.x0 >= fb.x0 - 0.5 and b.x1 <= fb.x1 + 0.5, (n, b)
        assert b.y0 >= fb.y0 - 0.5 and b.y1 <= fb.y1 + 0.5, (n, b)
    for i, n1 in enumerate(names):
        for n2 in names[i + 1 :]:
            assert not _overlaps(boxes[n1], boxes[n2]), (n1, n2, boxes[n1], boxes[n2])


# --- B4: strand drawing -------------------------------------------------------------------------


def _strand_out(tmp: Path) -> Path:
    return vl.build_out(
        tmp,
        genes=[
            vl.gene("Gp", "chr1", 10_201, 11_000, strand="+", log2=20.0),
            vl.gene("Gm", "chr1", 12_201, 13_000, strand="-", log2=21.0),
        ],
        regions_stats=[vl.region_stat("bgc_1", "bgc", "chr1", 10_001, 14_000, 0.01)],
        promoters=[
            ("chr1", 10_000, 10_200, "Pp", 0.8, "+"),
            ("chr1", 13_000, 13_200, "Pm", 0.9, "-"),
        ],
        terminators=[
            ("chr1", 11_000, 11_050, "Tp", 0.7, "+"),
            ("chr1", 12_100, 12_150, "Tm", 0.6, "-"),
        ],
    )


def test_b4_arrow_heads_point_downstream(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = _axes(_render(_strand_out(tmp_path), cfg, "bgc_1").figure, "genes")
    width = abs(ax.get_xlim()[1] - ax.get_xlim()[0])
    plus = _gene_patch(ax, 10_200, 11_000)
    minus = _gene_patch(ax, 12_200, 13_000)
    # the pointed end is a single tip vertex (zero y spread); the blunt end has height
    assert _end_spread(plus, "max", width) == pytest.approx(0, abs=1e-9)
    assert _end_spread(plus, "min", width) > 0
    assert _end_spread(minus, "min", width) == pytest.approx(0, abs=1e-9)
    assert _end_spread(minus, "max", width) > 0


def _track_centre(ax: Axes) -> float:
    ys = np.concatenate([_data_vertices(ax, p)[:, 1] for p in _opaque_patches(ax)])
    return float((ys.min() + ys.max()) / 2)


def _mark_lines(ax: Axes, colour: str, s: float, e: float) -> list[np.ndarray]:
    """``(x, y)`` arrays of every ``colour`` line on ``ax`` lying within 400 bp of ``[s, e)``."""
    out = []
    for ln in ax.get_lines():
        if _hex(ln.get_color()) != colour:
            continue
        x = np.asarray(ln.get_xdata(), dtype=float)
        y = np.asarray(ln.get_ydata(), dtype=float)
        if len(x) and x.min() >= s - 400 and x.max() <= e + 400:
            out.append(np.column_stack([x, y]))
    return out


@pytest.mark.parametrize(
    ("colour", "s", "e", "side"),
    [
        (PROMOTER, 10_000, 10_200, +1),
        (PROMOTER, 13_000, 13_200, -1),
        (TERMINATOR, 11_000, 11_050, +1),
        (TERMINATOR, 12_100, 12_150, -1),
    ],
    ids=["promoter+", "promoter-", "terminator+", "terminator-"],
)
def test_b4_marks_are_drawn_on_their_strand_side(
    tmp_path: Path, cfg: GiapConfig, colour: str, s: int, e: int, side: int
) -> None:
    ax = _axes(_render(_strand_out(tmp_path), cfg, "bgc_1").figure, "genes")
    centre = _track_centre(ax)
    lines = _mark_lines(ax, colour, s, e)
    assert lines
    ys = np.concatenate([ln[:, 1] for ln in lines])
    if side > 0:
        assert (ys > centre).all(), (ys, centre)
    else:
        assert (ys < centre).all(), (ys, centre)


def _stem_and_head(lines: list[np.ndarray]) -> tuple[float, float]:
    """The bent arrow's stem x (its vertical segment) and its far x (the head)."""
    stems = []
    for ln in lines:
        for (x0, y0), (x1, y1) in itertools.pairwise(ln):
            if x0 == x1 and y0 != y1:
                stems.append(x0)
    assert len(set(stems)) == 1, stems
    stem = stems[0]
    xs = np.concatenate([ln[:, 0] for ln in lines])
    head = float(xs[np.argmax(np.abs(xs - stem))])
    return float(stem), head


@pytest.mark.parametrize(
    ("s", "e", "strand"), [(10_000, 10_200, "+"), (13_000, 13_200, "-")], ids=["plus", "minus"]
)
def test_b4_promoter_bent_arrow_sits_at_its_strand_start_and_points_downstream(
    tmp_path: Path, cfg: GiapConfig, s: int, e: int, strand: str
) -> None:
    ax = _axes(_render(_strand_out(tmp_path), cfg, "bgc_1").figure, "genes")
    stem, head = _stem_and_head(_mark_lines(ax, PROMOTER, s, e))
    if strand == "+":
        assert stem == pytest.approx(e, abs=1)  # TSS: high end of the upstream window on +
        assert head > stem
    else:
        assert stem == pytest.approx(s, abs=1)  # TSS: low end of the upstream window on -
        assert head < stem


def test_b4_promoter_bent_arrow_points_downstream(tmp_path: Path, cfg: GiapConfig) -> None:
    """Direction alone (independent of where the stem sits)."""
    ax = _axes(_render(_strand_out(tmp_path), cfg, "bgc_1").figure, "genes")
    stem_p, head_p = _stem_and_head(_mark_lines(ax, PROMOTER, 10_000, 10_200))
    stem_m, head_m = _stem_and_head(_mark_lines(ax, PROMOTER, 13_000, 13_200))
    assert head_p > stem_p
    assert head_m < stem_m


# --- B5: imputed values are not shown -----------------------------------------------------------


def test_b5_imputed_gene_row_has_na_value_and_imputed_category(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    genes = [
        vl.gene("G1", "chr1", 10_101, 10_600, log2=21.0),
        vl.gene("IMP", "chr1", 10_701, 11_200, log2=21.0, imputed=True, status="not_detected"),
        vl.gene("IMPD", "chr1", 11_301, 11_800, log2=19.0, imputed=True, status="detected"),
    ]
    out = vl.build_out(
        tmp_path / "run",
        genes=genes,
        regions_stats=[vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 12_000, 0.01)],
    )
    inst = _render(out, cfg, "prophage_1")
    g = inst.data[inst.data["track"] == "gene"].set_index("id")
    assert set(g.index) == {"G1", "IMP", "IMPD"}
    assert float(g.loc["G1", "value"]) == 21.0
    assert g.loc["G1", "category"] == "detected"
    for tag in ("IMP", "IMPD"):
        assert pd.isna(g.loc[tag, "value"]), tag
        assert g.loc[tag, "category"] == "imputed", tag
    write_instance(inst, "V03", "locus", tmp_path / "plots", ["svg"], 100)
    lines = (tmp_path / "plots" / "data" / "V03_locus_prophage_1.tsv").read_text().splitlines()
    assert "gene\tG1\tchr1\t10101\t10600\t+\t21\tdetected" in lines
    assert "gene\tIMP\tchr1\t10701\t11200\t+\tNA\timputed" in lines
    assert "gene\tIMPD\tchr1\t11301\t11800\t+\tNA\timputed" in lines


# --- test-gap survivors -------------------------------------------------------------------------


def test_pct_rounds_half_up(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vl.build_out(
        tmp_path,
        genes=[vl.gene("G", "chr1", 10_101, 11_000, log2=20.0)],
        regions_stats=[
            vl.region_stat(
                "prophage_1",
                "prophage",
                "chr1",
                10_001,
                12_000,
                0.003,
                observed=12.5,
                expected=30.5,
                n_eligible=100,
            )
        ],
    )
    (lo,) = select_loci(RunView.from_outdir(out), cfg)
    assert lo.label == "prophage_1: 13 % detected vs 31 % expected (q = 0.003)"


def test_run_scale_with_equal_percentiles_is_one_unit_wide() -> None:
    genes = pd.DataFrame(
        {
            "detected": [True] * 5,
            "abundance_imputed": [False] * 5,
            "abundance_log2": [20.0] * 5,
        }
    )
    assert run_scale(genes) == pytest.approx((19.5, 20.5))


def test_run_scale_ignores_imputed_rows_even_when_detected(tmp_path: Path) -> None:
    genes = RunView.from_outdir(_scale_run(tmp_path)).table("genes_tsv")
    assert genes is not None
    assert bool(genes.set_index("locus_tag").loc["impdet", "detected"])
    assert run_scale(genes) == pytest.approx((18.12, 23.88))


def test_origin_line_is_at_the_replicon_length_to_under_1_bp(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    genes = [
        vl.gene("G1", "chr1", 98_001, 98_500, log2=20.0),
        vl.gene("GX", "chr1", 99_801, 100_200, log2=20.0, crosses_origin=True),
    ]
    out = vl.build_out(
        tmp_path, genes=genes, islands=[vl.island("ISL_1", "chr1", 99_001, 101_500, 0.01)]
    )
    xs = _dashed_x(_axes(_render(out, cfg, "ISL_1").figure, "genes"))
    assert len(xs) == 1
    assert abs(xs[0] - 100_000) < 1
