"""V2 ``abundance_landscape``: where along each replicon are proteins abundant? (M9-09; GIPPA_SPEC.md
§10.2 V2, §7.3.5(d) islands, M5-01 ``gippa.stats.smooth.loess_detrend``; lead decisions 2026-10-22).

V2 draws measured detected genes, a rug of undetected genes, a gene-level LOESS trend at the
islands' span fraction (``context.islands.detrend_span_frac``; legend "trend (LOESS, span <frac>)",
no claim) and shades islands with ``p_value ≤ stats.alpha.context``. oriC / terminus lines come from
the M9-03 support file, which is not merged: they are tested through their empty state only
(lead decision 3).

Oracles: the trend is pinned against statsmodels ``lowess(it=0, delta=0)`` at M5-01's span with
explicit circular padding (``helpers.viz_landscape.lowess_oracle``) and against the engine itself
(the spec's example 1 oracle, to 1e-12). Counts and spans come from the fixture rows.

Interface assumed where the spec leaves it open:

- The plot is registered as ``V02`` / ``abundance_landscape`` by some ``gippa.viz`` module (found by
  walking the package, as the ``plots`` step does); ``genes_tsv`` is in ``needs`` and ``islands``
  is not (a missing islands table is drawn without shading). It returns one double-width
  ``Instance`` with ``suffix=None``.
- Each replicon panel is an axes whose ``get_label()`` is the contig ID; panels run top to bottom
  in FASTA (preflight) order.
- Positions in the data TSV follow M9-03 example 1: ``position_bp = (start0 + end) / 2 mod L + 1``
  (fixture genes have integer midpoints). The tests read 1-based values from the written TSV only.
- Points are markers (Line2D or scatter) in ``#0072B2``; the trend is a black Line2D without
  markers; rug ticks carry their status colour; islands are hatched patches or collections.
"""

from __future__ import annotations

import math
import os
import tempfile
from pathlib import Path

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import pytest
from helpers.viz_landscape import (
    CHR1,
    correlogram_rows,
    gene,
    island,
    landscape_out,
    lowess_oracle,
    markers,
    na_correlogram_rows,
    plot_spec,
    text_overlaps,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.collections import LineCollection, PathCollection
from matplotlib.lines import Line2D
from matplotlib.text import Text

from gippa.config import GippaConfig, resolve_config
from gippa.io.plot_files import write_instance
from gippa.io.run_view import RunView
from gippa.stats.smooth import loess_detrend
from gippa.viz import style
from gippa.viz.base import Instance, Skip

REPO = Path(__file__).resolve().parents[3]
# product: M9-15 lead decision 1
COLUMNS = [
    *("contig", "kind", "locus_tag", "position_bp", "start", "end", "value"),
    "category",
    "product",
]
KINDS = ["point", "rug", "trend", "island", "oric", "ter"]
BLUE, GREY, ORANGE = "#0072B2", "#BBBBBB", "#E69F00"
EX1_MIDS = [10_000, 30_000, 50_000, 70_000, 90_000]
EX1_LOG2 = [1.0, 2.0, 3.0, 2.0, 1.0]
FRAC1 = "context.islands.detrend_span_frac=1.0"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


@pytest.fixture
def cfg1() -> GippaConfig:
    return resolve_config(overrides=[FRAC1])


def _render(out: Path, cfg: GippaConfig) -> Instance:
    result = plot_spec("V02").fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert isinstance(inst, Instance)
    return inst


def _ex1_genes() -> list[dict]:
    return [
        gene(f"g{i + 1}", mid, y) for i, (mid, y) in enumerate(zip(EX1_MIDS, EX1_LOG2, strict=True))
    ]


def _tsv(inst: Instance, tmp: Path) -> pd.DataFrame:
    """The written ``plots/data/V02_abundance_landscape.tsv`` (1-based closed)."""
    paths = write_instance(inst, "V02", "abundance_landscape", tmp, ["svg"], 72)
    tsv = next(p for p in paths if p.suffix == ".tsv")
    return pd.read_csv(tsv, sep="\t", na_values=["NA"], keep_default_na=False)


def _kind(frame: pd.DataFrame, kind: str, contig: str | None = None) -> pd.DataFrame:
    hit = frame[frame["kind"].astype(str) == kind]
    if contig is not None:
        hit = hit[hit["contig"].astype(str) == contig]
    return hit.reset_index(drop=True)


def _ax(inst: Instance, label: str) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_label() == label]
    assert len(found) == 1, f"axes {label!r}: {[a.get_label() for a in inst.figure.axes]}"
    return found[0]


def _hex(colour: object) -> str:
    return mcolors.to_hex(colour).upper()  # type: ignore[arg-type]


def _points(ax: Axes) -> list[tuple[float, float]]:
    return sorted((x, y) for x, y, hollow, c in markers(ax) if c == BLUE and not hollow)


def _trend_lines(ax: Axes) -> list[Line2D]:
    out = []
    for line in ax.get_lines():
        no_marker = line.get_marker() in (None, "", "None", " ", "none")
        if no_marker and _hex(line.get_color()) == "#000000" and len(line.get_xdata()) >= 2:
            out.append(line)
    return out


def _xs_with_colour(ax: Axes, colour: str) -> list[float]:
    """x of every line vertex, line-collection segment and scatter offset drawn in ``colour``."""
    xs: list[float] = []
    for line in ax.get_lines():
        cols = {_hex(line.get_color()), _hex(line.get_markeredgecolor())}
        if colour in cols:
            xs += [float(v) for v in line.get_xdata()]
    for coll in ax.collections:
        if isinstance(coll, LineCollection):
            cols = [_hex(c) for c in coll.get_colors()]
            for i, seg in enumerate(coll.get_segments()):
                if cols and cols[i % len(cols)] == colour:
                    xs += [float(p[0]) for p in seg]
        elif isinstance(coll, PathCollection):
            faces = [_hex(c) for c in coll.get_facecolors()]
            edges = [_hex(c) for c in coll.get_edgecolors()]
            for i, (x, _y) in enumerate(coll.get_offsets()):
                f = faces[i % len(faces)] if faces else None
                e = edges[i % len(edges)] if edges else None
                if colour in (f, e):
                    xs.append(float(x))
    return xs


def _spans(ax: Axes) -> list[tuple[str, float, float, str, float]]:
    """(hatch, x0, x1, face colour hex, alpha) of every hatched patch or collection path."""
    out = []
    for p in ax.patches:
        hatch = p.get_hatch()
        if not hatch:
            continue
        v = p.get_patch_transform().transform(p.get_path().vertices)
        rgba = mcolors.to_rgba(p.get_facecolor())
        alpha = p.get_alpha() if p.get_alpha() is not None else rgba[3]
        out.append((hatch, float(v[:, 0].min()), float(v[:, 0].max()), _hex(rgba), alpha))
    for coll in ax.collections:
        hatch = coll.get_hatch()
        if not hatch or isinstance(coll, PathCollection):
            continue
        faces = coll.get_facecolors()
        for i, path in enumerate(coll.get_paths()):
            v = path.vertices
            rgba = faces[i % len(faces)] if len(faces) else (0, 0, 0, 0)
            alpha = coll.get_alpha() if coll.get_alpha() is not None else float(rgba[3])
            out.append((hatch, float(v[:, 0].min()), float(v[:, 0].max()), _hex(rgba), alpha))
    return sorted(out, key=lambda s: (s[1], s[2]))


def _fig_texts(inst: Instance) -> list[str]:
    return [t.get_text() for t in inst.figure.findobj(Text) if t.get_text()]


def _isna(value: object) -> bool:
    return value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value))


# --- registration and skips -------------------------------------------------------------------


def test_v02_is_registered_as_abundance_landscape() -> None:
    spec = plot_spec("V02")
    assert spec.pid == "V02"
    assert spec.name == "abundance_landscape"
    assert "genes_tsv" in spec.needs
    assert "islands" not in spec.needs
    assert spec.module == "core"


def test_no_genes_table_is_skipped_as_input_absent(tmp_path: Path, cfg: GippaConfig) -> None:
    out = landscape_out(tmp_path, genes=None, islands=[])
    assert plot_spec("V02").fn(RunView.from_outdir(out), cfg) == Skip("input_absent:genes_tsv")


# --- examples 1-3 -----------------------------------------------------------------------------


def test_example1_points_and_trend_rows(tmp_path: Path, cfg1: GippaConfig) -> None:
    inst = _render(landscape_out(tmp_path, genes=_ex1_genes(), islands=[]), cfg1)
    assert inst.suffix is None
    assert list(inst.data.columns) == COLUMNS
    assert round(float(inst.figure.get_size_inches()[0]), 6) == style.DOUBLE_IN
    tsv = _tsv(inst, tmp_path / "plots")
    pts = _kind(tsv, "point")
    assert list(pts["locus_tag"]) == ["g1", "g2", "g3", "g4", "g5"]
    assert [int(v) for v in pts["position_bp"]] == [m + 1 for m in EX1_MIDS]
    assert list(pts["value"]) == pytest.approx(EX1_LOG2)
    assert set(pts["category"].astype(str)) == {"detected"}
    trend = _kind(tsv, "trend")
    assert [int(v) for v in trend["position_bp"]] == [m + 1 for m in EX1_MIDS]
    oracle = lowess_oracle(EX1_LOG2, EX1_MIDS, 1.0, True, 100_000)
    assert list(trend["value"]) == pytest.approx(list(oracle), rel=1e-5)


def test_trend_equals_the_engine_fit_to_1e_12(tmp_path: Path, cfg1: GippaConfig) -> None:
    inst = _render(landscape_out(tmp_path, genes=_ex1_genes(), islands=[]), cfg1)
    trend = _kind(inst.data, "trend")
    got = np.array([float(v) for v in trend["value"]])
    engine = loess_detrend(EX1_LOG2, EX1_MIDS, 1.0, circular=True, period=100_000).fit
    sm = lowess_oracle(EX1_LOG2, EX1_MIDS, 1.0, True, 100_000)
    np.testing.assert_allclose(got, engine, rtol=0, atol=1e-12)
    np.testing.assert_allclose(got, sm, rtol=0, atol=1e-9)
    line = _trend_lines(_ax(inst, "chr1"))
    assert len(line) == 1
    np.testing.assert_allclose([float(v) for v in line[0].get_ydata()], engine, atol=1e-12)
    xs = [float(v) for v in line[0].get_xdata()]
    assert xs == sorted(xs)  # position order, not joined across the origin
    assert xs[0] == pytest.approx(10.0, abs=1e-3)
    assert xs[-1] == pytest.approx(90.0, abs=1e-3)


def test_points_are_drawn_in_kb_and_blue(tmp_path: Path, cfg1: GippaConfig) -> None:
    inst = _render(landscape_out(tmp_path, genes=_ex1_genes(), islands=[]), cfg1)
    ax = _ax(inst, "chr1")
    pts = _points(ax)
    assert [x for x, _ in pts] == pytest.approx([m / 1000 for m in EX1_MIDS], abs=1e-3)
    assert [y for _, y in pts] == pytest.approx(EX1_LOG2)
    assert ax.get_xlabel() == "position (kb)"
    assert ax.get_ylabel() == "log2 abundance"


def test_legend_names_the_loess_span(tmp_path: Path, cfg: GippaConfig, cfg1: GippaConfig) -> None:
    for c, want in ((cfg1, "trend (LOESS, span 1)"), (cfg, "trend (LOESS, span 0.1)")):
        genes = [gene(f"g{i}", 2_000 * i + 1_000, float(i % 5)) for i in range(40)]
        inst = _render(landscape_out(tmp_path / want, genes=genes, islands=[]), c)
        assert want in _fig_texts(inst)
        assert "express" not in " ".join(_fig_texts(inst)).lower()


def test_example2_undetected_gene_is_one_rug_row(tmp_path: Path, cfg1: GippaConfig) -> None:
    genes = [*_ex1_genes(), gene("u1", 40_000, None, status="not_detected")]
    inst = _render(landscape_out(tmp_path, genes=genes, islands=[]), cfg1)
    tsv = _tsv(inst, tmp_path / "plots")
    rug = _kind(tsv, "rug")
    assert len(rug) == 1
    assert rug.loc[0, "locus_tag"] == "u1"
    assert rug.loc[0, "category"] == "not_detected"
    assert int(rug.loc[0, "position_bp"]) == 40_001
    assert len(_kind(tsv, "point")) == 5
    xs = _xs_with_colour(_ax(inst, "chr1"), GREY)
    assert any(abs(x - 40.0) < 1e-3 for x in xs)


def test_single_peptide_unestimated_rug_is_orange(tmp_path: Path, cfg1: GippaConfig) -> None:
    genes = [
        *_ex1_genes(),
        gene("u1", 40_000, None, status="not_detected"),
        gene("s1", 60_000, 4.0, status="single_peptide_unestimated"),
    ]
    inst = _render(landscape_out(tmp_path, genes=genes, islands=[]), cfg1)
    rug = _kind(_tsv(inst, tmp_path / "plots"), "rug")
    assert dict(zip(rug["locus_tag"], rug["category"], strict=True)) == {
        "u1": "not_detected",
        "s1": "single_peptide_unestimated",
    }
    ax = _ax(inst, "chr1")
    assert any(abs(x - 60.0) < 1e-3 for x in _xs_with_colour(ax, ORANGE))
    assert not any(abs(x - 60.0) < 1e-3 for x, _ in _points(ax))  # not a point


def test_example3_one_significant_island_and_title(tmp_path: Path, cfg1: GippaConfig) -> None:
    islands = [
        island("ISL_1", "high", 20_001, 40_000, 0.01),
        island("ISL_2", "low", 60_001, 80_000, 0.30),
    ]
    inst = _render(landscape_out(tmp_path, genes=_ex1_genes(), islands=islands), cfg1)
    tsv = _tsv(inst, tmp_path / "plots")
    isl = _kind(tsv, "island")
    assert len(isl) == 1
    assert (int(isl.loc[0, "start"]), int(isl.loc[0, "end"])) == (20_001, 40_000)
    assert isl.loc[0, "category"] == "high"
    assert inst.title == "1 high and 0 low abundance islands at FWER ≤ 0.05"
    assert inst.title in _fig_texts(inst)
    spans = _spans(_ax(inst, "chr1"))
    assert len(spans) == 1
    hatch, x0, x1, colour, alpha = spans[0]
    assert (hatch, colour) == (style.ISLAND["high"][1], style.ISLAND["high"][0])
    assert (x0, x1) == (pytest.approx(20.0, abs=1e-3), pytest.approx(40.0, abs=1e-3))
    assert alpha == pytest.approx(0.25)


def test_low_island_uses_the_low_hatch(tmp_path: Path, cfg1: GippaConfig) -> None:
    islands = [island("ISL_2", "low", 60_001, 80_000, 0.002)]
    inst = _render(landscape_out(tmp_path, genes=_ex1_genes(), islands=islands), cfg1)
    spans = _spans(_ax(inst, "chr1"))
    assert [(s[0], s[3]) for s in spans] == [(style.ISLAND["low"][1], style.ISLAND["low"][0])]
    assert inst.title == "0 high and 1 low abundance islands at FWER ≤ 0.05"


def test_island_alpha_follows_the_config_inclusively(tmp_path: Path) -> None:
    cfg = resolve_config(overrides=[FRAC1, "stats.alpha.context=0.01"])
    islands = [
        island("ISL_1", "high", 20_001, 40_000, 0.01),
        island("ISL_3", "high", 45_001, 55_000, 0.02),
    ]
    inst = _render(landscape_out(tmp_path, genes=_ex1_genes(), islands=islands), cfg)
    assert inst.title == "1 high and 0 low abundance islands at FWER ≤ 0.01"
    assert len(_kind(inst.data, "island")) == 1
    assert len(_spans(_ax(inst, "chr1"))) == 1


def test_title_adds_the_lag1_gene_label_morans_i(tmp_path: Path, cfg1: GippaConfig) -> None:
    rows = correlogram_rows(
        overrides={
            ("morans_i", "gene_label", 1): (0.31, -0.002, -0.05, 0.05, 0.002),
            ("morans_i", "tu_block", 1): (0.31, -0.004, -0.05, 0.12, 0.03),
        }
    )
    islands = [island("ISL_1", "high", 20_001, 40_000, 0.01)]
    out = landscape_out(tmp_path, genes=_ex1_genes(), islands=islands, autocorr=rows)
    inst = _render(out, cfg1)
    assert inst.title == (
        "1 high and 0 low abundance islands at FWER ≤ 0.05; lag-1 Moran's I 0.31 (q = 0.002)"
    )


def test_title_has_no_morans_clause_without_a_finite_q(tmp_path: Path, cfg1: GippaConfig) -> None:
    out = landscape_out(
        tmp_path, genes=_ex1_genes(), islands=[], autocorr=na_correlogram_rows("chr1")
    )
    assert _render(out, cfg1).title == "0 high and 0 low abundance islands at FWER ≤ 0.05"


# --- degenerate cases -------------------------------------------------------------------------


def test_origin_crossing_island_is_two_spans(tmp_path: Path, cfg1: GippaConfig) -> None:
    islands = [island("ISL_9", "high", 95_001, 105_000, 0.001)]
    inst = _render(landscape_out(tmp_path, genes=_ex1_genes(), islands=islands), cfg1)
    spans = [(round(s[1], 3), round(s[2], 3)) for s in _spans(_ax(inst, "chr1"))]
    assert spans == [(0.0, 5.0), (95.0, 100.0)]
    isl = _kind(_tsv(inst, tmp_path / "plots"), "island")
    assert len(isl) == 1  # decision 11: one row, end > length
    assert (int(isl.loc[0, "start"]), int(isl.loc[0, "end"])) == (95_001, 105_000)


def test_imputed_and_partial_values_are_never_drawn(tmp_path: Path, cfg1: GippaConfig) -> None:
    genes = [
        *_ex1_genes(),
        gene("imp", 60_000, 9.0, imputed=True),
        gene("par", 80_000, 8.0, partial=True),
        gene("nav", 85_000, None),  # detected, no value
    ]
    inst = _render(landscape_out(tmp_path, genes=genes, islands=[]), cfg1)
    pts = _kind(inst.data, "point")
    assert set(pts["locus_tag"]) == {"g1", "g2", "g3", "g4", "g5"}
    assert not ({"imp", "par", "nav"} & set(_kind(inst.data, "trend")["locus_tag"].dropna()))
    ys = [y for _, y in _points(_ax(inst, "chr1"))]
    assert 9.0 not in ys and 8.0 not in ys
    trend = [float(v) for v in _kind(inst.data, "trend")["value"]]
    engine = loess_detrend(EX1_LOG2, EX1_MIDS, 1.0, circular=True, period=100_000).fit
    np.testing.assert_allclose(trend, engine, atol=1e-12)


def test_fewer_than_three_points_draw_no_trend(tmp_path: Path, cfg1: GippaConfig) -> None:
    genes = [
        gene("g1", 10_000, 1.0),
        gene("g2", 30_000, 2.0),
        gene("u1", 50_000, None, status="not_detected"),
    ]
    inst = _render(landscape_out(tmp_path, genes=genes, islands=[]), cfg1)
    assert len(_kind(inst.data, "point")) == 2
    assert len(_kind(inst.data, "trend")) == 0
    assert len(_kind(inst.data, "rug")) == 1
    assert _trend_lines(_ax(inst, "chr1")) == []


def test_linear_contig_uses_an_unpadded_loess(tmp_path: Path, cfg1: GippaConfig) -> None:
    rep = ("lin1", 100_000, "linear", "chromosome")
    genes = [
        gene(f"g{i}", m, y, contig="lin1")
        for i, (m, y) in enumerate(zip(EX1_MIDS, EX1_LOG2, strict=True))
    ]
    inst = _render(landscape_out(tmp_path, replicons=[rep], genes=genes, islands=[]), cfg1)
    trend = [float(v) for v in _kind(inst.data, "trend")["value"]]
    linear = lowess_oracle(EX1_LOG2, EX1_MIDS, 1.0, False, 100_000)
    circular = lowess_oracle(EX1_LOG2, EX1_MIDS, 1.0, True, 100_000)
    assert np.max(np.abs(linear - circular)) > 0.05  # the two differ here
    np.testing.assert_allclose(trend, linear, atol=1e-9)


def test_origin_crossing_gene_midpoint_wraps(tmp_path: Path, cfg1: GippaConfig) -> None:
    """Midpoint unwrapped then mod L: a gene 99,501-100,500 (1-based) sits at position 1."""
    genes = [
        gene("ori", 100_000, 0.5, crosses_origin=True),
        *[gene(f"g{i}", m, y) for i, (m, y) in enumerate(zip(EX1_MIDS, EX1_LOG2, strict=True))],
    ]
    inst = _render(landscape_out(tmp_path, genes=genes, islands=[]), cfg1)
    pts = _kind(_tsv(inst, tmp_path / "plots"), "point")
    assert dict(zip(pts["locus_tag"], pts["position_bp"].astype(int), strict=True))["ori"] == 1
    pos = [0, *EX1_MIDS]
    ys = [0.5, *EX1_LOG2]
    trend = [float(v) for v in _kind(inst.data, "trend")["value"]]
    np.testing.assert_allclose(trend, lowess_oracle(ys, pos, 1.0, True, 100_000), atol=1e-9)
    xs = [x for x, _ in _points(_ax(inst, "chr1"))]
    assert min(xs) == pytest.approx(0.0, abs=1e-3)


def test_no_abundance_is_a_rug_for_every_gene(tmp_path: Path, cfg: GippaConfig) -> None:
    genes = [
        gene("d1", 10_000, None, abundance_type="none"),
        gene("d2", 30_000, None, abundance_type="none"),
        gene("u1", 50_000, None, status="not_detected", abundance_type="none"),
        gene("s1", 70_000, None, status="single_peptide_unestimated", abundance_type="none"),
    ]
    inst = _render(landscape_out(tmp_path, genes=genes, islands=[]), cfg)
    assert inst.title == "No measured abundance: detection only"
    assert set(_kind(inst.data, "rug")["locus_tag"]) == {"d1", "d2", "u1", "s1"}
    assert len(_kind(inst.data, "point")) == 0
    assert len(_kind(inst.data, "trend")) == 0
    ax = _ax(inst, "chr1")
    assert _points(ax) == []
    assert _trend_lines(ax) == []


def test_no_islands_table_draws_no_shading_and_says_so(tmp_path: Path, cfg1: GippaConfig) -> None:
    inst = _render(landscape_out(tmp_path, genes=_ex1_genes(), islands=None), cfg1)
    assert _spans(_ax(inst, "chr1")) == []
    assert len(_kind(inst.data, "island")) == 0
    notes = [t for t in _fig_texts(inst) if t != inst.title]
    assert any("island" in t.lower() for t in notes), notes


def test_without_support_files_there_are_no_landmark_lines(
    tmp_path: Path, cfg1: GippaConfig
) -> None:
    """M9-03 (#210) is not merged: no ``plots/support/replicon_landmarks.tsv`` -> no lines."""
    inst = _render(landscape_out(tmp_path, genes=_ex1_genes(), islands=[]), cfg1)
    kinds = set(inst.data["kind"].astype(str))
    assert not ({"oric", "ter"} & kinds)
    labels = [a.get_label() for ax in inst.figure.axes for a in ax.get_children()]
    assert not ({"oriC", "ter"} & set(labels))
    assert not ({"oriC", "ter"} & set(_fig_texts(inst)))


def test_one_panel_per_replicon_with_genes_in_fasta_order(tmp_path: Path, cfg1: GippaConfig) -> None:
    reps = [
        CHR1,
        ("plasmid_1", 8_000, "circular", "plasmid"),
        ("plasmid_2", 5_000, "circular", "plasmid"),
    ]
    genes = [
        *_ex1_genes(),
        *[
            gene(f"p{i}", m, 1.5, contig="plasmid_1", half=100)
            for i, m in enumerate((1_000, 3_000, 5_000, 7_000))
        ],
    ]
    inst = _render(landscape_out(tmp_path, replicons=reps, genes=genes, islands=[]), cfg1)
    labels = [ax.get_label() for ax in inst.figure.axes]
    assert "plasmid_2" not in labels
    top, bottom = _ax(inst, "chr1"), _ax(inst, "plasmid_1")
    assert top.get_position().y0 > bottom.get_position().y0
    assert len(_kind(inst.data, "point", "plasmid_1")) == 4
    assert len(_kind(inst.data, "trend", "plasmid_1")) == 4


def test_long_replicon_is_drawn_in_mb(tmp_path: Path, cfg1: GippaConfig) -> None:
    rep = ("big", 3_000_000, "circular", "chromosome")
    genes = [gene(f"g{i}", m, float(i)) for i, m in enumerate((500_000, 1_000_000, 2_000_000))]
    genes = [{**g, "contig": "big"} for g in genes]
    inst = _render(landscape_out(tmp_path, replicons=[rep], genes=genes, islands=[]), cfg1)
    ax = _ax(inst, "big")
    assert ax.get_xlabel() == "position (Mb)"
    assert [x for x, _ in _points(ax)] == pytest.approx([0.5, 1.0, 2.0], abs=1e-6)


def test_rows_sort_by_contig_kind_position_locus(tmp_path: Path, cfg1: GippaConfig) -> None:
    genes = [
        gene("g5", 90_000, 1.0),
        gene("u2", 20_000, None, status="not_detected"),
        *_ex1_genes()[:4],
        gene("u1", 20_000, None, status="not_detected", half=300),
    ]
    islands = [island("ISL_1", "high", 60_001, 80_000, 0.01)]
    inst = _render(landscape_out(tmp_path, genes=genes, islands=islands), cfg1)
    tsv = _tsv(inst, tmp_path / "plots")
    kinds = list(tsv["kind"].astype(str))
    assert kinds == sorted(kinds, key=KINDS.index)
    for kind in ("point", "rug", "trend"):
        rows = _kind(tsv, kind)
        keys = list(
            zip(rows["position_bp"].astype(int), rows["locus_tag"].astype(str), strict=True)
        )
        assert keys == sorted(keys), kind
    assert list(_kind(tsv, "rug")["locus_tag"]) == ["u1", "u2"]


def test_more_than_20000_points_are_binned(tmp_path: Path, cfg: GippaConfig) -> None:
    """20,001 measured genes on one replicon: at most 20,000 point rows, ``locus_tag`` NA, each
    the median of its equal-width bin (M9-01 ``binned`` edges)."""
    from gippa.viz.base import binned

    n, length = 20_001, 200_010
    mids = np.arange(n, dtype=np.int64) * 10 + 5
    ys = (np.arange(n) % 7).astype(float)
    rep = ("chrB", length, "circular", "chromosome")
    genes = [
        gene(f"b{i:05d}", int(m), float(y), contig="chrB", half=5)
        for i, (m, y) in enumerate(zip(mids, ys, strict=True))
    ]
    inst = _render(landscape_out(tmp_path, replicons=[rep], genes=genes, islands=[]), cfg)
    pts = _kind(inst.data, "point")
    assert 0 < len(pts) <= 20_000
    assert all(_isna(v) for v in pts["locus_tag"])
    edges = binned(mids, mids + 1, length)
    assert edges is not None
    which = np.searchsorted(edges, mids, side="right") - 1
    medians = sorted(float(np.median(ys[which == b])) for b in np.unique(which))
    assert sorted(float(v) for v in pts["value"]) == pytest.approx(medians)
    assert len(markers(_ax(inst, "chrB"))) <= 20_000


# --- property ---------------------------------------------------------------------------------

_STATUS = st.sampled_from(["detected", "detected", "not_detected", "single_peptide_unestimated"])


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    slots=st.lists(
        st.tuples(
            st.integers(min_value=1, max_value=98),
            _STATUS,
            st.booleans(),
            st.booleans(),
            st.integers(min_value=-40, max_value=40),
        ),
        min_size=1,
        max_size=30,
        unique_by=lambda t: t[0],
    ),
    frac=st.sampled_from([0.1, 0.3, 1.0]),
)
def test_property_points_rug_and_trend(
    slots: list[tuple[int, str, bool, bool, int]], frac: float
) -> None:
    """Points are exactly the detected, measured, non-imputed, non-partial genes; the rug is
    every ``detected = false`` gene, or every gene when no gene is a point (lead decisions 9 and
    24); the trend is the statsmodels lowess at the islands' span."""
    style.apply()
    cfg = resolve_config(overrides=[f"context.islands.detrend_span_frac={frac}"])
    genes, want_pts, want_rug = [], {}, set()
    for slot, status, imputed, partial, y10 in slots:
        mid, y = slot * 1_000, y10 / 10
        tag = f"g{slot:03d}"
        genes.append(gene(tag, mid, y, status=status, imputed=imputed, partial=partial, half=200))
        if status == "detected" and not imputed and not partial:
            want_pts[tag] = (mid, y)
        if status != "detected":
            want_rug.add(tag)
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(landscape_out(Path(tmp), genes=genes, islands=[]), cfg)
    pts = _kind(inst.data, "point")
    assert set(pts["locus_tag"]) == set(want_pts)
    if not want_pts:
        want_rug = {g["locus_tag"] for g in genes}
    assert set(_kind(inst.data, "rug")["locus_tag"]) == want_rug
    trend = [float(v) for v in _kind(inst.data, "trend")["value"]]
    order = sorted(want_pts.values())
    if len(order) < 3:
        assert trend == []
    else:
        fit = lowess_oracle([y for _, y in order], [m for m, _ in order], frac, True, 100_000)
        np.testing.assert_allclose(trend, fit, atol=1e-9)


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    islands=st.lists(
        st.tuples(
            st.sampled_from(["high", "low"]),
            st.sampled_from([0.001, 0.01, 0.04, 0.05, 0.06, 0.3, None]),
        ),
        max_size=6,
    ),
    alpha=st.sampled_from([0.01, 0.05, 0.1]),
)
def test_property_title_counts_islands_at_or_below_alpha(
    islands: list[tuple[str, float | None]], alpha: float
) -> None:
    style.apply()
    cfg = resolve_config(overrides=[FRAC1, f"stats.alpha.context={alpha}"])
    rows = [
        island(f"ISL_{i}", d, 10_001 + 12_000 * i, 20_000 + 12_000 * i, p)
        for i, (d, p) in enumerate(islands)
    ]
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(landscape_out(Path(tmp), genes=_ex1_genes(), islands=rows), cfg)
    sig = [d for d, p in islands if p is not None and p <= alpha]
    assert inst.title == (
        f"{sig.count('high')} high and {sig.count('low')} low abundance islands at FWER ≤ {alpha:g}"
    )
    assert len(_kind(inst.data, "island")) == len(sig)
    assert len(_spans(_ax(inst, "chr1"))) == len(sig)


# --- style, files and determinism -------------------------------------------------------------


def _rich(tmp: Path) -> Path:
    reps = [CHR1, ("plasmid_1", 8_000, "circular", "plasmid")]
    genes = [
        *[gene(f"g{i:02d}", 2_000 * i + 1_000, float(i % 6)) for i in range(48)],
        gene("u1", 97_000, None, status="not_detected"),
        gene("s1", 98_500, 2.0, status="single_peptide_unestimated", half=200),
        *[
            gene(f"p{i}", m, 1.0 + i, contig="plasmid_1", half=100)
            for i, m in enumerate((1_000, 3_000, 5_000, 7_000))
        ],
    ]
    islands = [
        island("ISL_1", "high", 20_001, 40_000, 0.01),
        island("ISL_2", "low", 95_001, 105_000, 0.02),
    ]
    rows = correlogram_rows(
        overrides={("morans_i", "gene_label", 1): (0.31, -0.002, -0.05, 0.05, 0.002)}
    )
    return landscape_out(tmp, replicons=reps, genes=genes, islands=islands, autocorr=rows)


def test_every_text_is_at_least_6_pt(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_rich(tmp_path), cfg)
    sizes = [t.get_fontsize() for t in inst.figure.findobj(Text) if t.get_text()]
    assert sizes
    assert min(sizes) >= 6


def test_no_text_overlaps_after_an_agg_draw(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_rich(tmp_path), cfg)
    assert text_overlaps(inst.figure) == []


def test_data_tsv_and_images_are_byte_identical(tmp_path: Path, cfg: GippaConfig) -> None:
    out = _rich(tmp_path / "run")
    written = []
    for name in ("x", "y"):
        inst = _render(out, cfg)
        paths = write_instance(
            inst, "V02", "abundance_landscape", tmp_path / name, ["svg", "png", "pdf"], 100
        )
        written.append({p.relative_to(tmp_path / name): p.read_bytes() for p in paths})
    assert set(written[0]) == {
        Path("V02_abundance_landscape.svg"),
        Path("V02_abundance_landscape.png"),
        Path("V02_abundance_landscape.pdf"),
        Path("data/V02_abundance_landscape.tsv"),
    }
    assert written[0] == written[1]
    tsv = written[0][Path("data/V02_abundance_landscape.tsv")].decode("utf-8").splitlines()
    assert tsv[0] == "\t".join(COLUMNS)


def test_render_does_not_mutate_between_calls(tmp_path: Path, cfg: GippaConfig) -> None:
    view = RunView.from_outdir(_rich(tmp_path))
    first = plot_spec("V02").fn(view, cfg)
    second = plot_spec("V02").fn(view, cfg)
    assert isinstance(first, list) and isinstance(second, list)
    pd.testing.assert_frame_equal(first[0].data, second[0].data)
    assert first[0].title == second[0].title


# --- example 6: the validation run ------------------------------------------------------------


def test_validation_like_run_is_rug_only(tmp_path: Path, cfg: GippaConfig) -> None:
    """Example 6 in miniature: measured values on genes none of which is detected (the
    2026-10-22 validation run), header-only islands, all-NA autocorrelation."""
    reps = [CHR1, ("plasmid_1", 8_000, "circular", "plasmid")]
    genes = [
        *[gene(f"g{i}", 4_000 * i + 2_000, 0.1 * i, status="not_detected") for i in range(20)],
        gene("p1", 4_000, -0.5, contig="plasmid_1", status="not_detected", half=100),
    ]
    rows = [*na_correlogram_rows("chr1"), *na_correlogram_rows("plasmid_1")]
    out = landscape_out(tmp_path, replicons=reps, genes=genes, islands=[], autocorr=rows)
    inst = _render(out, cfg)
    assert len(_kind(inst.data, "rug")) == 21
    assert len(_kind(inst.data, "point")) == 0
    assert len(_kind(inst.data, "trend")) == 0
    assert len(_kind(inst.data, "island")) == 0
    assert "Moran" not in inst.title


VAL_OUT = Path(os.environ.get("GIPPA_ROOT", REPO / "local")) / "runs" / "val-mvp-2026-10-22" / "OUT"


@pytest.mark.slow
@pytest.mark.skipif(not VAL_OUT.is_dir(), reason=f"validation run absent: {VAL_OUT}")
def test_example6_validation_run_is_rug_only(tmp_path: Path, cfg: GippaConfig) -> None:
    """87 genes (85 chr1, 2 plasmid_1), none detected; islands header-only."""
    inst = _render(VAL_OUT, cfg)
    assert len(_kind(inst.data, "rug")) == 87
    assert len(_kind(inst.data, "point")) == 0
    assert len(_kind(inst.data, "trend")) == 0
    assert {ax.get_label() for ax in inst.figure.axes} >= {"chr1", "plasmid_1"}
    paths = write_instance(inst, "V02", "abundance_landscape", tmp_path, ["svg", "png"], 100)
    assert all(p.stat().st_size > 0 for p in paths)
