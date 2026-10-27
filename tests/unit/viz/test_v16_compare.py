# ruff: noqa: RUF001, RUF002, RUF003  (ρ and − are the spec's own characters)
"""V16 ``compare``: UpSet of detected genes across runs and the pairwise abundance concordance
matrix (M9-23; GIPPA_SPEC.md §10.2 V16, §7.9.3; M9-index decision 29; lead decisions 1-2 of
2026-10-23).

Covers examples 1-3, every degenerate row (two runs, ρ NA, no detections, 12 runs), the data TSV,
the title (through ``gippa.viz.base.plural``, decision 1), the registry entry, the binning of large
scatters, determinism, the double width and the "no ``upsetplot`` import" AST rule. No ``--mpl``
baselines: M9-25 snapshots V16 on a three-run fixture (acceptance criterion 3), as M9-18 (#227)
left V12's baselines to M9-25.

Interface assumed where the spec leaves it open (M9-23 test author):

- ``gippa.viz.compare.compare_plot(detection, abundance, genes) -> Instance`` (pure):
  ``detection`` is ``compare_detection`` (both sections), ``abundance`` is ``compare_abundance``,
  ``genes`` maps each run label, in label order, to its ``genes.tsv`` frame (``RunView`` dtypes)
  for the upper-triangle scatters (the three compare tables hold no per-gene log2). One
  ``Instance`` with ``suffix=None``, ``coords=()``; the figure is double width (7.2 in) and draws
  ``Instance.title``.
- Importing ``gippa.viz.compare`` registers ``V16`` under the name ``compare`` (file stem
  ``V16_compare``); the registered ``(view, cfg)`` function returns a ``Skip`` on an ordinary run
  (V16 is written only by ``gippa compare``).
- Axes labels (``Axes.get_label()``): ``"a"`` the intersection bars (one ``Rectangle`` per drawn
  intersection, height = ``n_genes``), ``"a_sets"`` the set-size bars (one horizontal
  ``Rectangle`` per run, width = detected genes in that run), ``"b"`` the concordance matrix.
  Every other text (cell text, footnote, empty state) is a ``Text`` artist somewhere in the
  figure.
- Data columns ⊇ ``panel, intersection, n_genes, runs, run_a, run_b, spearman_rho, ci_low,
  ci_high, n_shared, locus_tag, x, y`` with ``panel`` first. Panel (a): one row per drawn
  intersection in drawing order, then (more than 30) one ``intersection = other`` row; ``runs``
  lists the intersection's labels (separated by ``,``, ``;`` or ``&``); a ``none`` row, if any,
  is ignored here. Panel (b): one row per unordered pair (label order) with ``locus_tag`` NA and
  the ρ columns, and one scatter row per shared measured gene of the pair (``locus_tag`` set,
  ``x`` = log2 in ``run_a``, ``y`` = log2 in ``run_b``).
- Numbers in the title and the cell text use the shortest ``%g`` form (``0.8``, ``0.93``);
  ``median pairwise ρ NA`` when no pair has a ρ.
"""

from __future__ import annotations

import ast
import math
import re
from collections.abc import Sequence
from pathlib import Path

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import pytest
from helpers import viz_compare as vc
from helpers.compare_data import gene, genes_frame
from helpers.viz_out_tree import make_out
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.collections import Collection
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text

from gippa.config import resolve_config
from gippa.io.plot_files import write_instance
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance, Skip

SRC = Path(__file__).resolve().parents[3] / "src"
COLUMNS = {
    "panel", "intersection", "n_genes", "runs", "run_a", "run_b", "spearman_rho", "ci_low",
    "ci_high", "n_shared", "locus_tag", "x", "y",
}  # fmt: skip
BLUE = "#0072B2"
GREY = "#BBBBBB"
EMPTY = "No gene detected in any run"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


# --- calling the plot -------------------------------------------------------------------------


def _plot(
    genes: dict[str, pd.DataFrame],
    detection: pd.DataFrame | None = None,
    abundance: pd.DataFrame | None = None,
) -> Instance:
    from gippa.viz.compare import compare_plot

    d, a = vc.inputs(genes)
    inst = compare_plot(detection if detection is not None else d,
                        abundance if abundance is not None else a, genes)  # fmt: skip
    assert isinstance(inst, Instance)
    assert inst.suffix is None
    return inst


def _isna(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


def _panel(inst: Instance, p: str) -> pd.DataFrame:
    return inst.data[inst.data["panel"].astype(str) == p].reset_index(drop=True)


def _bars(inst: Instance) -> pd.DataFrame:
    a = _panel(inst, "a")
    return a[a["intersection"].astype(str) != "none"].reset_index(drop=True)


def _stats(inst: Instance) -> pd.DataFrame:
    b = _panel(inst, "b")
    return b[[_isna(v) for v in b["locus_tag"]]].reset_index(drop=True)


def _scatter(inst: Instance, ra: str, rb: str) -> pd.DataFrame:
    b = _panel(inst, "b")
    keep = [not _isna(t) and x == ra and y == rb
            for t, x, y in zip(b["locus_tag"], b["run_a"], b["run_b"], strict=True)]  # fmt: skip
    return b[keep].reset_index(drop=True)


def _axes(fig: Figure, label: str) -> Axes:
    found = [ax for ax in fig.axes if ax.get_label() == label]
    assert len(found) == 1, [ax.get_label() for ax in fig.axes]
    return found[0]


def _rects(ax: Axes) -> list[Rectangle]:
    return [p for p in ax.patches if isinstance(p, Rectangle)]


def _texts(fig: Figure) -> list[str]:
    return [t.get_text() for t in fig.findobj(Text) if t.get_text()]


def _runs(value: object) -> set[str]:
    return {s for s in re.split(r"[,;&]\s*", str(value)) if s}


def _hex(c: object) -> str:
    return mcolors.to_hex(c).upper()  # type: ignore[arg-type]


def _colours(fig: Figure) -> set[str]:
    out: set[str] = set()
    for c in fig.findobj(Collection):
        out |= {_hex(x) for x in c.get_facecolors()}
    for ln in fig.findobj(Line2D):
        if ln.get_marker() not in (None, "", "None", " "):
            out.add(_hex(ln.get_markerfacecolor()))
    return out


def _with_rho(genes: dict[str, pd.DataFrame], cells: Sequence[tuple]) -> Instance:
    d, a = vc.inputs(genes)
    for ra, rb, rho, lo, hi in cells:
        a = vc.set_rho(a, ra, rb, rho, lo, hi)
    return _plot(genes, d, a)


def _three_full() -> dict[str, pd.DataFrame]:
    """Three runs detecting and measuring all eight genes (every pair has a ρ)."""
    return vc.detect_runs({"A1": set(range(8)), "A2": set(range(8)), "A3": set(range(8))}, 8)


# --- plural helper (lead decision 1) -----------------------------------------------------------


@pytest.mark.parametrize(
    ("n", "text"), [(0, "0 genes"), (1, "1 gene"), (2, "2 genes"), (30, "30 genes")]
)
def test_plural(n: int, text: str) -> None:
    assert base.plural(n, "gene") == text


def test_plural_other_words() -> None:
    assert base.plural(1, "run") == "1 run"
    assert base.plural(3, "run") == "3 runs"


@given(st.integers(0, 10_000), st.sampled_from(["gene", "run", "island", "peptide"]))
def test_plural_property(n: int, word: str) -> None:
    text = base.plural(n, word)
    assert text == (f"{n} {word}" if n == 1 else f"{n} {word}s")


# --- registry and the data TSV ---------------------------------------------------------------


def test_v16_is_registered_as_compare() -> None:
    import gippa.viz.compare  # noqa: F401  (registers V16)

    spec = base.registered()["V16"]
    assert spec.pid == "V16"
    assert spec.name == "compare"
    assert base.normalise_ids("v16") == ["V16"]


def test_v16_is_not_drawn_for_a_single_run(tmp_path: Path) -> None:
    """V16 is written only by ``gippa compare`` (§10.2): on an ordinary run it is a skip."""
    import gippa.viz.compare  # noqa: F401

    view = RunView.from_outdir(make_out(tmp_path))
    result = base.registered()["V16"].fn(view, resolve_config())
    assert isinstance(result, Skip), result


def test_data_tsv_columns_and_panel_order(tmp_path: Path) -> None:
    inst = _plot(vc.example_1())
    assert set(inst.data.columns) >= COLUMNS, sorted(COLUMNS - set(inst.data.columns))
    assert inst.data.columns[0] == "panel"
    panels = [str(p) for p in inst.data["panel"]]
    assert set(panels) == {"a", "b"}
    assert panels == sorted(panels)
    paths = write_instance(inst, "V16", "compare", tmp_path / "plots", ["svg"], 50)
    tsv = tmp_path / "plots" / "data" / "V16_compare.tsv"
    assert tsv in paths and (tmp_path / "plots" / "V16_compare.svg") in paths
    back = pd.read_csv(tsv, sep="\t", na_values=["NA"], keep_default_na=False)
    assert list(back.columns) == list(inst.data.columns)
    assert len(back) == len(inst.data)


def test_title_is_in_the_svg(tmp_path: Path) -> None:
    inst = _plot(vc.example_1())
    assert inst.title in _texts(inst.figure)
    write_instance(inst, "V16", "compare", tmp_path, ["svg"], 50)
    svg = (tmp_path / "V16_compare.svg").read_text(encoding="utf-8")
    assert inst.title in svg


def test_double_width() -> None:
    inst = _plot(vc.example_1())
    assert round(float(inst.figure.get_size_inches()[0]), 6) == 7.2


def test_two_renders_are_identical(tmp_path: Path) -> None:
    genes = _three_full()
    a, b = _plot(genes), _plot(genes)
    pd.testing.assert_frame_equal(a.data, b.data)
    pa = write_instance(a, "V16", "compare", tmp_path / "a", ["svg"], 50)
    pb = write_instance(b, "V16", "compare", tmp_path / "b", ["svg"], 50)
    assert pa[0].read_bytes() == pb[0].read_bytes()


# --- example 1 (M9-22 example 1) ---------------------------------------------------------------


def test_example_1_bars_and_order() -> None:
    inst = _plot(vc.example_1())
    bars = _bars(inst)
    assert [str(x) for x in bars["intersection"]] == ["A1", "A1&A2", "A1&A2&A3", "A2"]
    assert [int(n) for n in bars["n_genes"]] == [1, 1, 1, 1]
    assert [_runs(r) for r in bars["runs"]] == [{"A1"}, {"A1", "A2"}, {"A1", "A2", "A3"}, {"A2"}]
    heights = [r.get_height() for r in _rects(_axes(inst.figure, "a")) if r.get_height() > 0]
    assert heights == [1, 1, 1, 1]


def test_example_1_set_sizes() -> None:
    inst = _plot(vc.example_1())
    widths = [abs(r.get_width()) for r in _rects(_axes(inst.figure, "a_sets"))]
    assert sorted(w for w in widths if w > 0) == [1, 3, 3]


def test_example_1_none_is_a_footnote_not_a_bar() -> None:
    inst = _plot(vc.example_1())
    assert "none" not in [str(x) for x in _bars(inst)["intersection"]]
    assert len([r for r in _rects(_axes(inst.figure, "a")) if r.get_height() > 0]) == 4
    notes = [t for t in _texts(inst.figure) if re.search(r"none|no run", t, re.I)]
    assert notes, _texts(inst.figure)
    assert any(re.search(r"\b1\b", t) for t in notes), notes


def test_example_1_dot_matrix_colours() -> None:
    colours = _colours(_plot(vc.example_1()).figure)
    assert BLUE in colours and GREY in colours, colours


def test_example_1_title() -> None:
    """All three pairs share fewer than three measured genes, so every ρ is NA."""
    inst = _plot(vc.example_1())
    assert inst.title == "1 gene detected in all 3 runs; median pairwise ρ NA"


def test_title_median_rho_and_plural() -> None:
    genes = _three_full()
    inst = _with_rho(genes, [("A1", "A2", 0.5, 0.2, 0.7), ("A1", "A3", 0.8, 0.6, 0.93),
                             ("A2", "A3", 0.9, 0.8, 0.95)])  # fmt: skip
    assert inst.title == "8 genes detected in all 3 runs; median pairwise ρ 0.8"


def test_title_with_no_gene_in_all_runs() -> None:
    genes = vc.detect_runs({"A1": {0, 1, 2, 3}, "A2": {0, 1, 2, 3}, "A3": {4}}, 5)
    assert _plot(genes).title.startswith("0 genes detected in all 3 runs; median pairwise ρ ")


def test_title_has_no_expressed_wording() -> None:
    inst = _plot(_three_full())
    assert "express" not in " ".join([inst.title, *_texts(inst.figure)]).lower()


# --- example 2: cell text ---------------------------------------------------------------------


def test_example_2_cell_text() -> None:
    inst = _with_rho(_three_full(), [("A1", "A2", 0.8, 0.6, 0.93)])
    assert "ρ 0.8 [0.6, 0.93]" in _texts(inst.figure)


def test_every_pair_has_its_rho_and_ci_in_the_figure() -> None:
    inst = _with_rho(_three_full(), [("A1", "A2", 0.5, 0.2, 0.7), ("A1", "A3", 0.8, 0.6, 0.93),
                                     ("A2", "A3", 0.9, 0.8, 0.95)])  # fmt: skip
    texts = _texts(inst.figure)
    for want in ("ρ 0.5 [0.2, 0.7]", "ρ 0.8 [0.6, 0.93]", "ρ 0.9 [0.8, 0.95]"):
        assert want in texts, texts


def test_negative_rho_is_printed_not_hidden() -> None:
    inst = _with_rho(_three_full(), [("A1", "A2", -0.5, -0.8, -0.1)])
    assert any(re.fullmatch(r"ρ [-−]0\.5 \[[-−]0\.8, [-−]0\.1\]", t) for t in _texts(inst.figure))


def test_panel_b_stat_rows_copy_compare_abundance() -> None:
    genes = _three_full()
    d, a = vc.inputs(genes)
    inst = _plot(genes, d, a)
    rows = _stats(inst)
    assert list(zip(rows["run_a"], rows["run_b"], strict=True)) == [
        ("A1", "A2"), ("A1", "A3"), ("A2", "A3")]  # fmt: skip
    for col in ("spearman_rho", "ci_low", "ci_high"):
        np.testing.assert_allclose(rows[col].astype(float), a[col].astype(float), rtol=1e-12)
    assert [int(n) for n in rows["n_shared"]] == [int(n) for n in a["n_shared"]]


def test_panel_b_scatter_rows_are_shared_measured_genes() -> None:
    rows = [
        [gene(f"g{i}", start=1000 * i, log2=float(i), tu=f"TU_{i}") for i in range(6)],
        [gene("g0", start=0, log2=0.5, tu="TU_0"),
         gene("g1", start=1000, detected=False, tu="TU_1"),
         gene("g2", start=2000, log2=2.5, imputed=True, tu="TU_2"),
         *[gene(f"g{i}", start=1000 * i, log2=i + 0.5, tu=f"TU_{i}") for i in range(3, 6)]],
    ]  # fmt: skip
    genes = {"A1": genes_frame(rows[0]), "A2": genes_frame(rows[1])}
    sc = _scatter(_plot(genes), "A1", "A2")
    assert sorted(sc["locus_tag"]) == ["g0", "g3", "g4", "g5"]
    got = {t: (float(x), float(y)) for t, x, y in zip(sc["locus_tag"], sc["x"], sc["y"],
                                                      strict=True)}  # fmt: skip
    assert got == {"g0": (0.0, 0.5), "g3": (3.0, 3.5), "g4": (4.0, 4.5), "g5": (5.0, 5.5)}


# --- example 3 and the 12-run case: the 30 largest ---------------------------------------------


def test_example_3_forty_intersections() -> None:
    genes = vc.many_intersections(40, 6)
    labels = list(genes)
    inst = _plot(genes)
    bars = _bars(inst)
    want = [vc.labels_of(m, labels) for m in range(40, 10, -1)]  # n_genes 40 .. 11
    drawn = bars[bars["intersection"].astype(str) != "other"]
    assert [str(x) for x in drawn["intersection"]] == want
    assert [int(n) for n in drawn["n_genes"]] == list(range(40, 10, -1))
    other = bars[bars["intersection"].astype(str) == "other"]
    assert len(other) == 1
    assert int(other["n_genes"].iloc[0]) == sum(range(1, 11))
    assert bars.index[bars["intersection"].astype(str) == "other"][0] == len(bars) - 1
    assert len([r for r in _rects(_axes(inst.figure, "a")) if r.get_height() > 0]) == 30


def test_thirty_intersections_have_no_other_row() -> None:
    inst = _plot(vc.many_intersections(30, 5))
    bars = _bars(inst)
    assert len(bars) == 30
    assert "other" not in [str(x) for x in bars["intersection"]]


def test_twelve_runs_draw_the_thirty_largest() -> None:
    genes = vc.many_intersections(45, 12)
    inst = _plot(genes)
    bars = _bars(inst)
    drawn = bars[bars["intersection"].astype(str) != "other"]
    assert len(drawn) == 30
    assert [int(n) for n in drawn["n_genes"]] == list(range(45, 15, -1))
    assert int(bars["n_genes"].iloc[-1]) == sum(range(1, 16))
    assert len([r for r in _rects(_axes(inst.figure, "a")) if r.get_height() > 0]) == 30
    assert len(_stats(inst)) == 12 * 11 // 2


# --- degenerate cases -------------------------------------------------------------------------


def test_two_runs() -> None:
    genes = vc.detect_runs({"A1": {0, 1, 2, 3, 4}, "A2": {2, 3, 4, 5, 6}}, 8)
    inst = _plot(genes)
    assert len(_bars(inst)) <= 3
    assert [str(x) for x in _bars(inst)["intersection"]] == ["A1&A2", "A1", "A2"]
    rows = _stats(inst)
    assert list(zip(rows["run_a"], rows["run_b"], strict=True)) == [("A1", "A2")]
    assert inst.title.startswith("3 genes detected in all 2 runs; median pairwise ρ ")


def test_rho_na_too_few_shared_cell() -> None:
    genes = vc.detect_runs({"A1": {0, 1, 2, 3}, "A2": {3, 4, 5, 6}}, 8)  # 1 shared gene
    d, a = vc.inputs(genes)
    assert "too_few_shared" in str(a["flags"].iloc[0])
    inst = _plot(genes, d, a)
    assert "ρ NA" in _texts(inst.figure)  # M9-23 decision 5: the reason moves to a footnote
    assert any("too few shared" in t for t in _texts(inst.figure))
    assert _isna(_stats(inst)["spearman_rho"].iloc[0])


def test_no_detection_anywhere_is_the_empty_state() -> None:
    genes = vc.detect_runs({"A1": set(), "A2": set(), "A3": set()}, 5)
    inst = _plot(genes)
    assert EMPTY in _texts(inst.figure)
    assert inst.title == EMPTY
    assert inst.data.empty or len(_bars(inst)) == 0


def test_large_scatter_is_binned() -> None:
    """More than 20,000 shared points: no artist draws more than 20,000 elements (M9-01)."""
    n = 20_001
    genes = vc.detect_runs({"A1": set(range(n)), "A2": set(range(n))}, n)
    d, a = vc.inputs(genes, B=2)
    inst = _plot(genes, d, a)
    for c in inst.figure.findobj(Collection):
        assert len(c.get_offsets()) <= 20_000, c
    for ln in inst.figure.findobj(Line2D):
        assert len(ln.get_xdata()) <= 20_000, ln


# --- properties -------------------------------------------------------------------------------


@settings(max_examples=12, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    st.integers(2, 5).flatmap(
        lambda k: st.lists(st.lists(st.booleans(), min_size=k, max_size=k), min_size=1, max_size=40)
    )
)
def test_panel_a_counts_add_up(det: list[list[bool]]) -> None:
    k = len(det[0])
    labels = [f"A{i + 1}" for i in range(k)]
    sets = {lab: {g for g, row in enumerate(det) if row[i]} for i, lab in enumerate(labels)}
    genes = vc.detect_runs(sets, len(det))
    inst = _plot(genes)
    n_any = sum(any(row) for row in det)
    if n_any == 0:
        assert inst.title == EMPTY
        return
    bars = _bars(inst)
    assert sum(int(n) for n in bars["n_genes"]) == n_any
    drawn = bars[bars["intersection"].astype(str) != "other"]
    assert len(drawn) <= 30
    counts = [int(n) for n in drawn["n_genes"]]
    assert counts == sorted(counts, reverse=True)
    n_all = sum(all(row) for row in det)
    assert inst.title.startswith(f"{base.plural(n_all, 'gene')} detected in all {k} runs; ")


# --- no UpSetPlot (index decision 29) ---------------------------------------------------------


def test_no_upsetplot_import_in_src() -> None:
    bad: list[str] = []
    for path in sorted(SRC.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            names: list[str] = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                names = [node.module]
            elif isinstance(node, ast.Call):
                func = node.func
                fname = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", "")
                if fname in ("import_module", "__import__"):
                    names = [str(a.value) for a in node.args if isinstance(a, ast.Constant)]
            bad += [f"{path}:{n}" for n in names if n.split(".")[0] == "upsetplot"]
    assert bad == []
    assert (SRC / "gippa" / "viz" / "compare.py").is_file()
