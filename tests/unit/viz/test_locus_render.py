"""V3 ``locus_<id>`` rendering (M9-10 "Rendering", examples 3-6; GIAP_SPEC.md §10.2 V3, §10.1,
§6.5; M9-index decisions 9-13, 20).

Tests the plotted data (the ``Instance`` data frame = ``plots/data/V03_locus_<id>.tsv``), the
window, the colours the spec names, the figure structure, the title and determinism. No
``--mpl`` baselines (M9-25 owns them).

Interface assumed where the spec leaves it open (reported to the lead):

- ``render_locus(locus, view, cfg, scale)`` returns one ``Instance`` with ``suffix =
  locus_id``, ``title = locus.label`` (drawn in the figure) and ``"start"`` in ``coords``.
- Data columns are exactly ``track, id, contig, start, end, strand, value, category``; rows carry
  each feature's own coordinates (not clipped to the window), unwrapped into the window frame on
  circular contigs (a feature past the origin is shifted by ± the replicon length; ``start`` may
  be negative or ``end`` past the length). ``gene`` rows: ``id = locus_tag``, ``value`` =
  measured ``abundance_log2``, ``category = detection_status``. ``pileup`` rows: maximal runs of
  constant non-zero depth, ``value`` = depth. ``promoter`` / ``terminator`` rows: the BED
  ``name`` / ``score`` / ``strand``. ``region`` rows: ``category = region_type``; ``island`` rows:
  ``category = direction``. Every peptide_hits row counts as one accepted placement.
- The gene track is the axes labelled ``"genes"`` with x limits = the window (its width in bp is
  ``xlim[1] - xlim[0]``, whatever the offset); the pileup is the axes labelled ``"pileup"`` with
  y label ``"peptides"``. The origin is a dashed (``"--"``) vertical ``Line2D`` on the gene axes.
- Gene arrow fills are opaque (alpha > 0.5) patches or patch collections on the gene axes;
  region shading is drawn with alpha 0.2.
- ``scale`` maps to cividis with clipping (``Normalize(lo, hi)`` then the colormap, values outside
  clipped to the ends); the plot computes it once per run as the 2nd-98th percentile
  (``numpy.percentile``) of ``abundance_log2`` over genes with ``detected`` true,
  ``abundance_imputed`` false and a finite ``abundance_log2``.
"""

from __future__ import annotations

import itertools
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import pytest
from helpers import viz_locus as vl
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.collections import Collection
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.text import Text

from giap.config import GiapConfig, resolve_config
from giap.io.plot_files import write_instance
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.base import Instance, Skip
from giap.viz.locus import render_locus, select_loci

COLUMNS = ["track", "id", "contig", "start", "end", "strand", "value", "category"]
CIVIDIS = matplotlib.colormaps["cividis"]
SCALE = (18.0, 24.0)


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


def _axes(fig: Figure, label: str) -> Axes:
    found = [ax for ax in fig.axes if ax.get_label() == label]
    assert len(found) == 1, [ax.get_label() for ax in fig.axes]
    return found[0]


def _rows(inst: Instance, track: str) -> list[tuple[Any, ...]]:
    """``(id, start, end)`` of ``track``, sorted by start then id (0-based frame values)."""
    d = inst.data[inst.data["track"] == track]
    return sorted(
        ((str(r.id), int(r.start), int(r.end)) for r in d.itertuples()), key=lambda t: (t[1], t[0])
    )


def _width_bp(ax: Axes) -> float:
    lo, hi = ax.get_xlim()
    return abs(hi - lo)


def _origin_fracs(ax: Axes) -> list[float]:
    """Axes-fraction x of every dashed vertical line on ``ax``."""
    lo, hi = ax.get_xlim()
    fracs = []
    for ln in ax.get_lines():
        x = np.asarray(ln.get_xdata(), dtype=float)
        if ln.get_linestyle() == "--" and len(x) >= 2 and np.allclose(x, x[0]):
            fracs.append(float((x[0] - lo) / (hi - lo)))
    return fracs


def _fills(ax: Axes) -> list[tuple[tuple[float, ...], tuple[float, ...] | None, str | None]]:
    """(facecolour, edgecolour, hatch) of every patch and every collection member on ``ax``."""
    out: list[tuple[tuple[float, ...], tuple[float, ...] | None, str | None]] = []
    for p in ax.patches:
        out.append((tuple(p.get_facecolor()), tuple(p.get_edgecolor()), p.get_hatch()))
    for c in ax.collections:
        if not isinstance(c, Collection):
            continue
        fcs, ecs = c.get_facecolors(), c.get_edgecolors()
        for i, fc in enumerate(fcs):
            ec = tuple(ecs[i % len(ecs)]) if len(ecs) else None
            out.append((tuple(fc), ec, c.get_hatch()))
    return out


def _opaque(ax: Axes) -> list[tuple[tuple[float, ...], tuple[float, ...] | None, str | None]]:
    return [f for f in _fills(ax) if len(f[0]) == 4 and f[0][3] > 0.5]


def _rgb(c: Sequence[float]) -> tuple[float, ...]:
    return tuple(round(float(v), 4) for v in c[:3])


LUT = {_rgb(CIVIDIS(i)) for i in range(CIVIDIS.N)}


def _cividis(v: float, lo: float, hi: float) -> tuple[float, ...]:
    return _rgb(CIVIDIS(float(np.clip((v - lo) / (hi - lo), 0.0, 1.0))))


def _hex_colours(fig: Figure) -> set[str]:
    found: set[str] = set()
    for art in fig.findobj():
        cols: list[Any] = []
        if isinstance(art, Line2D):
            cols += [art.get_color(), art.get_markerfacecolor(), art.get_markeredgecolor()]
        elif isinstance(art, Patch):
            cols += [art.get_facecolor(), art.get_edgecolor()]
        elif isinstance(art, Collection):
            cols += list(art.get_facecolors()) + list(art.get_edgecolors())
        for c in cols:
            try:
                found.add(mcolors.to_hex(c).upper())
            except (ValueError, TypeError):
                continue
    return found


def _texts(fig: Figure) -> list[str]:
    return [t.get_text() for t in fig.findobj(Text) if t.get_text()]


# --- example 3: window, genes, width ---------------------------------------------------------


def _example3_out(tmp: Path) -> Path:
    genes = [
        vl.gene("A", "chr1", 7_501, 8_200, log2=20.0),  # overlaps the window start
        vl.gene("B", "chr1", 37_101, 38_000, log2=20.0),  # after the window
        vl.gene("C", "chr1", 7_001, 8_000, log2=20.0),  # ends exactly at the window start
        vl.gene("D", "chr1", 37_001, 37_500, log2=20.0),  # starts exactly at the window end
        vl.gene("E", "chr1", 36_901, 37_000, log2=20.0),  # last base inside
        vl.gene("F", "chr1", 20_001, 21_000, strand="-", log2=21.0),
    ]
    return vl.build_out(
        tmp,
        genes=genes,
        regions_stats=[vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 35_000, 0.003)],
        regions=[
            vl.region("prophage_1", "prophage", "chr1", 10_001, 35_000),
            vl.region("prophage_9", "prophage", "chr1", 80_001, 81_000),
        ],
    )


def test_example3_window_draws_overlapping_genes_only(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_example3_out(tmp_path), cfg, "prophage_1")
    assert _rows(inst, "gene") == [
        ("A", 7_500, 8_200),
        ("F", 20_000, 21_000),
        ("E", 36_900, 37_000),
    ]
    assert _width_bp(_axes(inst.figure, "genes")) == pytest.approx(29_000, abs=1)


def test_example3_region_rows_in_the_window(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_example3_out(tmp_path), cfg, "prophage_1")
    assert _rows(inst, "region") == [("prophage_1", 10_000, 35_000)]
    reg = inst.data[inst.data["track"] == "region"]
    assert list(reg["category"]) == ["prophage"]


def test_example3_is_double_width_and_has_no_origin_line(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_example3_out(tmp_path), cfg, "prophage_1")
    assert round(float(inst.figure.get_size_inches()[0]), 6) == style.DOUBLE_IN
    assert _origin_fracs(_axes(inst.figure, "genes")) == []


def test_data_columns_suffix_title_and_coords(tmp_path: Path, cfg: GiapConfig) -> None:
    out = _example3_out(tmp_path)
    view = RunView.from_outdir(out)
    (lo,) = select_loci(view, cfg)
    inst = render_locus(lo, view, cfg, SCALE)
    assert list(inst.data.columns) == COLUMNS
    assert inst.suffix == "prophage_1"
    assert inst.title == lo.label
    assert inst.title == "prophage_1: 12 % detected vs 31 % expected (q = 0.003)"
    assert inst.title in _texts(inst.figure)
    assert "start" in inst.coords
    assert "end" not in inst.coords
    assert set(inst.data["track"]) <= {"gene", "pileup", "promoter", "terminator", "region",
                                       "island"}  # fmt: skip


def test_gene_rows_carry_strand_abundance_and_status(tmp_path: Path, cfg: GiapConfig) -> None:
    genes = [
        vl.gene("G1", "chr1", 10_101, 11_000, log2=21.0),
        vl.gene("G2", "chr1", 11_101, 11_900, strand="-", status="not_detected"),
    ]
    out = vl.build_out(
        tmp_path / "run",
        genes=genes,
        regions_stats=[vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 12_000, 0.01)],
    )
    inst = _render(out, cfg, "prophage_1")
    g = inst.data[inst.data["track"] == "gene"].sort_values("start")
    assert list(g["id"]) == ["G1", "G2"]
    assert list(g["strand"]) == ["+", "-"]
    assert list(g["category"]) == ["detected", "not_detected"]
    assert float(g["value"].iloc[0]) == 21.0
    assert pd.isna(g["value"].iloc[1])
    write_instance(inst, "V03", "locus", tmp_path / "plots", ["svg"], 100)
    tsv = (tmp_path / "plots" / "data" / "V03_locus_prophage_1.tsv").read_text(encoding="utf-8")
    lines = tsv.splitlines()
    assert lines[0] == "\t".join(COLUMNS)
    assert "gene\tG1\tchr1\t10101\t11000\t+\t21\tdetected" in lines
    assert "gene\tG2\tchr1\t11101\t11900\t-\tNA\tnot_detected" in lines


# --- example 4: origin crossing ----------------------------------------------------------------


def _example4_out(tmp: Path) -> Path:
    genes = [
        vl.gene("G1", "chr1", 98_001, 98_500, log2=20.0),
        vl.gene("GX", "chr1", 99_801, 100_200, log2=20.0, crosses_origin=True),
        vl.gene("G2", "chr1", 501, 900, log2=20.0),  # past the origin: unwrapped to 100,500
        vl.gene("GF", "chr1", 50_001, 50_500, log2=20.0),  # far away
    ]
    return vl.build_out(
        tmp, genes=genes, islands=[vl.island("ISL_1", "chr1", 99_001, 101_500, 0.01)]
    )


def test_example4_window_is_unwrapped_across_the_origin(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_example4_out(tmp_path), cfg, "ISL_1")
    genes_ax = _axes(inst.figure, "genes")
    assert _width_bp(genes_ax) == pytest.approx(6_500, abs=1)
    assert _rows(inst, "gene") == [
        ("G1", 98_000, 98_500),
        ("GX", 99_800, 100_200),
        ("G2", 100_500, 100_900),
    ]
    assert _rows(inst, "island") == [("ISL_1", 99_000, 101_500)]
    isl = inst.data[inst.data["track"] == "island"]
    assert list(isl["category"]) == ["high"]


def test_example4_origin_line_at_100000(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_example4_out(tmp_path), cfg, "ISL_1")
    fracs = _origin_fracs(_axes(inst.figure, "genes"))
    assert len(fracs) == 1
    assert fracs[0] == pytest.approx((100_000 - 97_000) / 6_500, abs=1e-3)


def test_example4_is_single_width(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_example4_out(tmp_path), cfg, "ISL_1")
    assert round(float(inst.figure.get_size_inches()[0]), 6) == style.SINGLE_IN


def test_island_shading_uses_the_direction_hatch(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vl.build_out(
        tmp_path, islands=[vl.island("ISL_1", "chr1", 20_001, 24_000, 0.01, direction="low")]
    )
    inst = _render(out, cfg, "ISL_1")
    colour, hatch = style.ISLAND["low"]
    shaded = [f for f in _fills(_axes(inst.figure, "genes")) if f[2] == hatch]
    assert shaded
    assert any(_rgb(f[0]) == _rgb(mcolors.to_rgba(colour)) or
               (f[1] is not None and _rgb(f[1]) == _rgb(mcolors.to_rgba(colour)))
               for f in shaded)  # fmt: skip


# --- circular below 0 and linear clipping --------------------------------------------------------


def _edge_out(tmp: Path, contig: str) -> Path:
    length = 100_000 if contig == "chr1" else 30_000
    genes = [
        vl.gene("NEAR", contig, 1_001, 1_500, log2=20.0),
        vl.gene("TAIL", contig, length - 499, length - 200, log2=20.0),
    ]
    return vl.build_out(
        tmp,
        genes=genes,
        regions_stats=[
            vl.region_stat("prophage_1", "prophage", contig, 501, 5_000, 0.01),
            vl.region_stat("prophage_2", "prophage", contig, length - 2_999, length - 1_000, 0.02),
        ],
    )


def test_circular_window_extends_below_zero(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_edge_out(tmp_path, "chr1"), cfg, "prophage_1")
    genes_ax = _axes(inst.figure, "genes")
    assert _width_bp(genes_ax) == pytest.approx(8_500, abs=1)  # [-1,500, 7,000)
    assert _rows(inst, "gene") == [("TAIL", -500, -200), ("NEAR", 1_000, 1_500)]
    fracs = _origin_fracs(genes_ax)
    assert len(fracs) == 1
    assert fracs[0] == pytest.approx(1_500 / 8_500, abs=1e-3)


def test_linear_window_is_clipped_at_zero(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_edge_out(tmp_path, "lin1"), cfg, "prophage_1")
    genes_ax = _axes(inst.figure, "genes")
    assert _width_bp(genes_ax) == pytest.approx(7_000, abs=1)  # [0, 7,000)
    assert _rows(inst, "gene") == [("NEAR", 1_000, 1_500)]
    assert _origin_fracs(genes_ax) == []


def test_linear_window_is_clipped_at_the_contig_end(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_edge_out(tmp_path, "lin1"), cfg, "prophage_2")
    # span [27,000, 29,000) ± 2 kb -> [25,000, 31,000) clipped to [25,000, 30,000)
    assert _width_bp(_axes(inst.figure, "genes")) == pytest.approx(5_000, abs=1)
    assert _rows(inst, "gene") == [("TAIL", 29_500, 29_800)]


# --- widths and height --------------------------------------------------------------------------


@pytest.mark.parametrize(("end", "width"), [(18_000, style.SINGLE_IN), (18_001, style.DOUBLE_IN)])
def test_width_switches_above_20_kb(
    tmp_path: Path, cfg: GiapConfig, end: int, width: float
) -> None:
    # span [2,000, end) ± 2 kb -> a window of exactly 20,000 or 20,001 bp
    out = vl.build_out(
        tmp_path, regions_stats=[vl.region_stat("prophage_1", "prophage", "chr1", 2_001, end, 0.01)]
    )
    inst = _render(out, cfg, "prophage_1")
    assert round(float(inst.figure.get_size_inches()[0]), 6) == width


def test_height_is_fixed(tmp_path: Path, cfg: GiapConfig) -> None:
    small = _render(_example4_out(tmp_path / "a"), cfg, "ISL_1")
    large = _render(_example3_out(tmp_path / "b"), cfg, "prophage_1")
    assert small.figure.get_size_inches()[1] == pytest.approx(large.figure.get_size_inches()[1])


# --- example 5 and the gene colours -------------------------------------------------------------


def _one_gene_out(tmp: Path, **kw: Any) -> Path:
    return vl.build_out(
        tmp,
        genes=[vl.gene("G", "chr1", 10_101, 11_000, **kw)],
        regions_stats=[vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 12_000, 0.01)],
    )


def test_example5_abundance_above_the_scale_is_the_top_of_cividis(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    inst = _render(_one_gene_out(tmp_path, log2=25.0), cfg, "prophage_1", (18.0, 24.0))
    fills = {_rgb(f[0]) for f in _opaque(_axes(inst.figure, "genes"))}
    assert _rgb(CIVIDIS(1.0)) in fills


def test_abundance_below_the_scale_is_the_bottom_of_cividis(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    inst = _render(_one_gene_out(tmp_path, log2=10.0), cfg, "prophage_1", (18.0, 24.0))
    fills = {_rgb(f[0]) for f in _opaque(_axes(inst.figure, "genes"))}
    assert _rgb(CIVIDIS(0.0)) in fills


def test_the_scale_argument_sets_the_colour(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_one_gene_out(tmp_path, log2=25.0), cfg, "prophage_1", (18.0, 30.0))
    fills = {_rgb(f[0]) for f in _opaque(_axes(inst.figure, "genes"))}
    assert _cividis(25.0, 18.0, 30.0) in fills
    assert _rgb(CIVIDIS(1.0)) not in fills


def test_undetected_gene_is_grey(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_one_gene_out(tmp_path, status="not_detected"), cfg, "prophage_1")
    fills = {_rgb(f[0]) for f in _opaque(_axes(inst.figure, "genes"))}
    assert _rgb(mcolors.to_rgba("#BBBBBB")) in fills
    assert not fills & LUT


def test_single_peptide_unestimated_has_an_orange_outline(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_one_gene_out(tmp_path, status="single_peptide_unestimated"), cfg, "prophage_1")
    edges = {_rgb(f[1]) for f in _opaque(_axes(inst.figure, "genes")) if f[1] is not None}
    assert _rgb(mcolors.to_rgba("#E69F00")) in edges


def test_imputed_gene_is_hatched_grey_not_cividis(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(
        _one_gene_out(tmp_path, log2=21.0, imputed=True, status="not_detected"), cfg, "prophage_1"
    )
    opaque = _opaque(_axes(inst.figure, "genes"))
    assert any(f[2] for f in opaque)
    assert not {_rgb(f[0]) for f in opaque} & LUT


def _scale_run(tmp: Path) -> Path:
    """101 measured genes, log2 = 18 + 0.06 i at 0-based [900 i, 900 i + 600); 2nd-98th
    percentile = (18.12, 23.88). Two imputed outliers and an undetected gene with an abundance
    must not move the scale."""
    genes = [
        vl.gene(f"m{i:03d}", "chr1", 900 * i + 1, 900 * i + 600, log2=18.0 + 0.06 * i)
        for i in range(101)
    ]
    genes += [
        vl.gene("imp_lo", "chr1", 95_001, 95_600, log2=5.0, imputed=True, status="not_detected"),
        vl.gene("imp_hi", "chr1", 96_001, 96_600, log2=40.0, imputed=True, status="not_detected"),
        vl.gene("undet", "chr1", 97_001, 97_600, log2=39.0, status="not_detected"),
    ]
    return vl.build_out(
        tmp,
        genes=genes,
        regions_stats=[
            vl.region_stat("prophage_1", "prophage", "chr1", 88_201, 88_800, 0.01),  # gene 98
            vl.region_stat("prophage_2", "prophage", "chr1", 1_801, 2_400, 0.02),  # gene 2
        ],
    )


def _expected_colours(window: tuple[int, int]) -> set[tuple[float, ...]]:
    lo, hi = 18.12, 23.88
    return {
        _cividis(18.0 + 0.06 * i, lo, hi)
        for i in range(101)
        if 900 * i < window[1] and 900 * i + 600 > window[0]
    }


def test_one_colour_scale_is_shared_by_all_locus_plots(tmp_path: Path, cfg: GiapConfig) -> None:
    import giap.viz.locus  # noqa: F401  (registers V03)

    result = base.registered()["V03"].fn(RunView.from_outdir(_scale_run(tmp_path)), cfg)
    assert not isinstance(result, Skip)
    by_id = {i.suffix: i for i in result if i.suffix is not None}
    assert set(by_id) == {"prophage_1", "prophage_2"}
    windows = {"prophage_1": (86_200, 90_800), "prophage_2": (-200, 4_400)}
    for lid, window in windows.items():
        fills = {_rgb(f[0]) for f in _opaque(_axes(by_id[lid].figure, "genes"))}
        assert fills & LUT == _expected_colours(window), lid


# --- example 6: the peptide pileup --------------------------------------------------------------


def _pileup_out(tmp: Path, hits: Sequence[dict[str, object]] | None) -> Path:
    return vl.build_out(
        tmp,
        genes=[vl.gene("G", "chr1", 2_101, 2_900, log2=20.0)],
        regions_stats=[vl.region_stat("prophage_1", "prophage", "chr1", 2_001, 3_000, 0.01)],
        hits=hits,
    )


def _pileup_runs(inst: Instance) -> list[tuple[int, int, int]]:
    d = inst.data[inst.data["track"] == "pileup"]
    return sorted((int(r.start), int(r.end), int(r.value)) for r in d.itertuples())


def test_example6_pileup_depth(tmp_path: Path, cfg: GiapConfig) -> None:
    hits = [
        vl.hit("PEPTIDEK", "chr1", 101, 120),
        vl.hit("PEPTIDEK", "chr1", 111, 130),
        vl.hit("OUTSIDEK", "chr1", 60_001, 60_020),  # outside the window
        vl.hit("OTHERCK", "lin1", 101, 120),  # another contig
    ]
    inst = _render(_pileup_out(tmp_path, hits), cfg, "prophage_1")
    # 1-based 101-110 depth 1, 111-120 depth 2, 121-130 depth 1
    assert _pileup_runs(inst) == [(100, 110, 1), (110, 120, 2), (120, 130, 1)]
    pile = _axes(inst.figure, "pileup")
    assert pile.get_ylabel() == "peptides"


def test_no_peptide_hits_means_no_pileup_rows(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_pileup_out(tmp_path, None), cfg, "prophage_1")
    assert _pileup_runs(inst) == []
    assert len(inst.data[inst.data["track"] == "gene"]) == 1


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(st.lists(st.tuples(st.integers(1, 400), st.integers(1, 60)), min_size=1, max_size=12))
def test_property_pileup_runs_are_the_per_base_depth(placements: list[tuple[int, int]]) -> None:
    hits = [vl.hit("PEPTIDEK", "chr1", s, s + n - 1) for s, n in placements]
    depth = np.zeros(500, dtype=int)
    for s, n in placements:
        depth[s - 1 : s - 1 + n] += 1
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(_pileup_out(Path(tmp), hits), resolve_config(), "prophage_1")
    runs = _pileup_runs(inst)
    drawn = np.zeros(500, dtype=int)
    for a, b, v in runs:
        assert 0 <= a < b <= 500 and v >= 1
        assert not drawn[a:b].any()  # runs are disjoint
        drawn[a:b] = v
    assert (drawn == depth).all()
    for (_, b1, v1), (a2, _, v2) in itertools.pairwise(runs):
        assert not (b1 == a2 and v1 == v2)  # maximal runs


# --- promoters, terminators, regions ------------------------------------------------------------


def _marks_out(tmp: Path) -> Path:
    return vl.build_out(
        tmp,
        genes=[vl.gene("G", "chr1", 2_201, 2_800, log2=20.0)],
        regions_stats=[vl.region_stat("bgc_1", "bgc", "chr1", 2_001, 3_000, 0.01)],
        regions=[
            vl.region("bgc_1", "bgc", "chr1", 2_001, 3_000),
            vl.region("prophage_9", "prophage", "chr1", 80_001, 81_000),
        ],
        islands=[vl.island("ISL_1", "chr1", 2_501, 4_000, 0.5, direction="high")],
        promoters=[
            ("chr1", 2_100, 2_150, "P1", 0.8, "+"),
            ("chr1", 50_000, 50_050, "P2", 0.9, "+"),
            ("lin1", 2_100, 2_150, "P3", 0.9, "+"),
        ],
        terminators=[("chr1", 2_900, 2_950, "T1", 0.7, "-")],
    )


def test_promoter_and_terminator_rows_are_strand_aware(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_marks_out(tmp_path), cfg, "bgc_1")
    assert _rows(inst, "promoter") == [("P1", 2_100, 2_150)]
    assert _rows(inst, "terminator") == [("T1", 2_900, 2_950)]
    d = inst.data.set_index("id")
    assert d.loc["P1", "strand"] == "+"
    assert d.loc["T1", "strand"] == "-"
    assert float(d.loc["P1", "value"]) == pytest.approx(0.8)
    assert float(d.loc["T1", "value"]) == pytest.approx(0.7)


def test_promoter_and_terminator_colours(tmp_path: Path, cfg: GiapConfig) -> None:
    colours = _hex_colours(_render(_marks_out(tmp_path), cfg, "bgc_1").figure)
    assert "#009E73" in colours
    assert "#D55E00" in colours


def test_region_and_island_rows_and_region_shading(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_marks_out(tmp_path), cfg, "bgc_1")
    assert _rows(inst, "region") == [("bgc_1", 2_000, 3_000)]
    assert _rows(inst, "island") == [("ISL_1", 2_500, 4_000)]
    d = inst.data
    assert list(d[d["track"] == "region"]["category"]) == ["bgc"]
    assert list(d[d["track"] == "island"]["category"]) == ["high"]
    shade = [f for f in _fills(_axes(inst.figure, "genes")) if len(f[0]) == 4]
    want = _rgb(mcolors.to_rgba(style.REGION["bgc"]))
    assert any(_rgb(f[0]) == want and f[0][3] == pytest.approx(0.2) for f in shade)


# --- labels ---------------------------------------------------------------------------------------


def test_genes_are_labelled_by_name_else_locus_tag(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vl.build_out(
        tmp_path,
        genes=[
            vl.gene("g1", "chr1", 10_101, 10_900, name="dnaA", log2=20.0),
            vl.gene("g2", "chr1", 11_001, 11_800, log2=20.0),
        ],
        regions_stats=[vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 12_000, 0.01)],
    )
    texts = _texts(_render(out, cfg, "prophage_1").figure)
    assert "dnaA" in texts
    assert "g2" in texts
    assert "g1" not in texts


def test_large_bgc_draws_200_genes_with_at_most_12_labels(tmp_path: Path, cfg: GiapConfig) -> None:
    genes = [
        vl.gene(f"t{i:03d}", "chr1", 10_001 + 400 * i, 10_000 + 400 * i + 350, name=f"bg{i:03d}",
                log2=20.0)
        for i in range(200)
    ]  # fmt: skip
    out = vl.build_out(
        tmp_path,
        genes=genes,
        regions_stats=[vl.region_stat("bgc_1", "bgc", "chr1", 10_001, 90_000, 0.01)],
    )
    inst = _render(out, cfg, "bgc_1")
    assert len(inst.data[inst.data["track"] == "gene"]) == 200
    names = {f"bg{i:03d}" for i in range(200)} | {f"t{i:03d}" for i in range(200)}
    n_labels = sum(t in names for t in _texts(inst.figure))
    assert 1 <= n_labels <= 12


# --- house rules ----------------------------------------------------------------------------------


def test_every_text_is_at_least_6_pt(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_marks_out(tmp_path), cfg, "bgc_1")
    sizes = [t.get_fontsize() for t in inst.figure.findobj(Text) if t.get_text()]
    assert sizes
    assert min(sizes) >= 6


def test_render_leaves_no_pyplot_figure(tmp_path: Path, cfg: GiapConfig) -> None:
    import matplotlib.pyplot as plt

    before = plt.get_fignums()
    inst = _render(_example3_out(tmp_path), cfg, "prophage_1")
    assert isinstance(inst.figure, Figure)
    assert plt.get_fignums() == before


def test_files_are_byte_identical_across_renders(tmp_path: Path, cfg: GiapConfig) -> None:
    out = _marks_out(tmp_path / "run")
    written = []
    for name in ("x", "y"):
        inst = _render(out, cfg, "bgc_1")
        paths = write_instance(inst, "V03", "locus", tmp_path / name, ["svg", "png", "pdf"], 100)
        written.append({p.relative_to(tmp_path / name): p.read_bytes() for p in paths})
    assert set(written[0]) == {
        Path("V03_locus_bgc_1.svg"),
        Path("V03_locus_bgc_1.png"),
        Path("V03_locus_bgc_1.pdf"),
        Path("data/V03_locus_bgc_1.tsv"),
    }
    assert written[0] == written[1]


def test_render_does_not_mutate_the_view(tmp_path: Path, cfg: GiapConfig) -> None:
    view = RunView.from_outdir(_marks_out(tmp_path))
    (lo,) = [x for x in select_loci(view, cfg) if x.locus_id == "bgc_1"]
    genes_before = view.table("genes_tsv")
    first = render_locus(lo, view, cfg, SCALE)
    second = render_locus(lo, view, cfg, SCALE)
    pd.testing.assert_frame_equal(first.data, second.data)
    pd.testing.assert_frame_equal(view.table("genes_tsv"), genes_before)
