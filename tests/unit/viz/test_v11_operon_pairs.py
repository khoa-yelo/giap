"""V11 ``operon_pairs``: do operon partners look alike? (M9-11; GIPPA_SPEC.md §10.2 V11, §7.4,
§7.4.5, §9.3 ``pairs.tsv`` / ``operon_flags.tsv``; M9-index decision 16).

Tests the plotted data (the ``Instance`` data frame = ``plots/data/V11_operon_pairs.tsv``), the
panel structure, the title, the degenerate cases, determinism and the colours the spec names. No
``--mpl`` baselines (Khoa approves those at M9-25).

Interface assumed where the spec leaves it open:

- The plot lives in ``gippa.viz.operon_pairs``; importing it registers ``V11`` (name
  ``operon_pairs``, module ``core``, ``needs=("pairs",)``, ``optional`` ⊇ {``operon_flags``}). It
  returns one ``Instance`` with ``suffix=None`` and ``coords=()``. The module exports
  ``FLAG_MARKERS``: the 7 §7.4.5 flag values -> 7 distinct matplotlib marker codes.
- The three panel axes carry ``Axes.get_label()`` ``"a"``, ``"b"``, ``"c"``.
- Data columns are exactly ``panel, pair_id, class, abs_delta_log2, lfdr, upstream_log2,
  downstream_log2, flag, q_value, bin, density, n`` (``n`` last). Panel (a): per class present
  among {OP, NO, CV, SW}, 30 bin rows (``pair_id`` NA, integer ``bin`` contiguous, ``density``);
  plus always five count rows, one per class OP, NO, CV, SW, DV (``bin``, ``density``,
  ``abs_delta_log2`` NA; ``n`` = pairs of that class with a finite ``delta_log2``). ``bin`` may
  start at 0 or 1. Panel (b): one row per pair with a finite ``lfdr``. Panel (c): one row per
  OP / NO pair with both ``abundance_log2`` finite; ``flag`` / ``q_value`` come from
  ``operon_flags.tsv`` (``flag`` NA for a pair that has no flag row).
- **Upstream** in panel (c) is ``gene_a`` on both strands: ``pairs.tsv`` already orients OP / NO
  pairs in transcription order (schema: "gene_a: Upstream gene in transcription order"; M6-04
  ``gippa.operons.pairs`` swaps the genes of a ``-`` pair). The spec's examples 2-3 swap again on
  ``-``; that conflict is reported to the lead (``test_upstream_is_gene_a_on_both_strands``).
- Panel (c) draws its points with ``scatter`` (``PathCollection``); flagged pairs at
  ``q_value <= stats.alpha.operon_flags`` are the points drawn with a marker size larger than the
  smallest size in the panel and a non-zero edge width; one marker shape per flag.
- The title's ``k`` counts ``operon_flags.tsv`` rows with ``q_value <= alpha``; the medians and
  ``alpha`` are printed with ``%g``-style shortest form (``0.4``, ``1.3``, ``0.1``).
- A header-only ``pairs.tsv`` is the framework empty state (zero rows, the same columns).
"""

from __future__ import annotations

import math
import os
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import pytest
from helpers import flags_data
from helpers.viz_out_tree import make_out, write_schema_table
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.collections import PathCollection
from matplotlib.lines import Line2D
from matplotlib.text import Text

from gippa.config import GippaConfig, resolve_config
from gippa.io.plot_files import write_instance
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance, Skip

REPO = Path(__file__).resolve().parents[3]
VAL_OUT = Path(os.environ.get("GIPPA_ROOT", REPO / "local")) / "runs" / "val-mvp-2026-10-22" / "OUT"
COLUMNS = [
    "panel",
    "pair_id",
    "class",
    "abs_delta_log2",
    "lfdr",
    "upstream_log2",
    "downstream_log2",
    "flag",
    "q_value",
    "bin",
    "density",
    "n",
]
CLASS_ORDER = ("OP", "NO", "CV", "SW", "DV")
CLASS_COLOUR = {"OP": "#0072B2", "NO": "#E69F00", "CV": "#009E73", "SW": "#CC79A7"}
FLAGS = (
    "internal_promoter_supported",
    "internal_terminator_supported",
    "unexplained_discordance",
    "conflict_terminator_in_operon",
    "novel_operon_candidate",
    "coabundant_despite_terminator",
    "other",
)
N_BINS = 30


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


# --- builders -------------------------------------------------------------------------------


def _pair(
    n: int,
    cls: str,
    a: float | None = None,
    b: float | None = None,
    *,
    strand: str | None = None,
    lfdr: float | None = None,
    q: float | None = None,
) -> dict[str, object]:
    """One ``pairs.tsv`` row; ``delta_log2 = b - a`` when both abundances are given."""
    measured = a is not None and b is not None
    test = "none"
    if lfdr is not None:
        test = "missed_operon_decoy"
    elif q is not None:
        test = "both_detected_t"
    return {
        "pair_id": f"PAIR_{n}",
        "class": cls,
        "gene_a": f"G{n}a",
        "gene_b": f"G{n}b",
        "strand": strand if strand is not None else ("." if cls in ("CV", "DV") else "+"),
        "intergenic_bp": 40,
        "detected_a": a is not None,
        "detected_b": b is not None,
        "abundance_log2_a": a,
        "abundance_log2_b": b,
        "delta_log2": (b - a) if measured else None,  # type: ignore[operator]
        "test": test,
        "stat": abs(b - a) if measured and lfdr is not None else None,  # type: ignore[operator]
        "q_value": q,
        "lfdr": lfdr,
    }


def _delta_pairs(cls: str, deltas: Sequence[float], start: int) -> list[dict[str, object]]:
    return [_pair(start + i, cls, 20.0, 20.0 + d) for i, d in enumerate(deltas)]


def _flag_row(
    pair: Mapping[str, object], flag: str, q: float, lfdr: float | None = None
) -> dict[str, object]:
    return {
        "pair_id": pair["pair_id"],
        "gene_a": pair["gene_a"],
        "gene_b": pair["gene_b"],
        "uniop_class": pair["class"],
        "promoter_between_yes": False,
        "terminator_between_yes": False,
        "abundance_pattern": "discordant",
        "flag": flag,
        "q_value": q,
        "lfdr": lfdr,
        "evidence_note": f"q={q:g}",
    }


def _out(
    tmp: Path,
    pairs: Sequence[Mapping[str, object]] | None,
    flags: Sequence[Mapping[str, object]] | None = (),
) -> Path:
    """``flags=None`` leaves ``operon_flags.tsv`` out; ``()`` writes it header-only."""
    out = make_out(tmp)
    if pairs is not None:
        write_schema_table(out / "tables" / "pairs.tsv", "pairs", pairs)
    if flags is not None:
        write_schema_table(out / "tables" / "operon_flags.tsv", "operon_flags", flags)
    return out


#: example 1: OP Δ -0.2, 0.4, 0.6; NO Δ 1.0, -1.6
EXAMPLE1 = [*_delta_pairs("OP", [-0.2, 0.4, 0.6], 1), *_delta_pairs("NO", [1.0, -1.6], 4)]


def _full_pairs() -> list[dict[str, object]]:
    """All five classes, NO pairs with lfdr, one minus-strand OP pair."""
    return [
        _pair(1, "OP", 20.0, 22.0, strand="+"),
        _pair(2, "OP", 20.0, 22.0, strand="-"),
        _pair(3, "OP", 18.0, 21.0, q=0.04),
        _pair(4, "OP", 19.0, 19.5, q=0.2),
        _pair(5, "NO", 20.0, 21.0, lfdr=0.9, q=0.6),
        _pair(6, "NO", 20.0, 20.2, lfdr=0.4, q=0.3),
        _pair(7, "CV", 20.0, 21.5),
        _pair(8, "SW", 20.0, 18.5),
        _pair(9, "DV", 20.0, 25.0),
    ]


def _full_flags(pairs: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    """Examples 4-5: PAIR_3 ``unexplained_discordance`` q 0.04; PAIR_4 q 0.2."""
    by_id = {p["pair_id"]: p for p in pairs}
    return [
        _flag_row(by_id["PAIR_3"], "unexplained_discordance", 0.04),
        _flag_row(by_id["PAIR_4"], "conflict_terminator_in_operon", 0.2),
    ]


# --- calling the plot -------------------------------------------------------------------------


def _v11() -> base.PlotSpec:
    import gippa.viz.operon_pairs  # noqa: F401  (registers V11)

    return base.registered()["V11"]


def _render(out: Path, cfg: GippaConfig) -> Instance:
    result = _v11().fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert isinstance(inst, Instance)
    return inst


def _isna(value: object) -> bool:
    return value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value))


def _num(value: object) -> float:
    return float(str(value))


def _rows(inst: Instance, panel: str) -> pd.DataFrame:
    return inst.data[inst.data["panel"] == panel].reset_index(drop=True)


def _bin_rows(inst: Instance, cls: str | None = None) -> pd.DataFrame:
    a = _rows(inst, "a")
    sel = a[[not _isna(v) for v in a["bin"]]]
    if cls is not None:
        sel = sel[sel["class"] == cls]
    return sel.reset_index(drop=True)


def _counts(inst: Instance) -> dict[str, int]:
    a = _rows(inst, "a")
    sel = a[[_isna(v) for v in a["bin"]]]
    assert len(sel) == len(set(sel["class"])), sel
    return {str(c): int(_num(n)) for c, n in zip(sel["class"], sel["n"], strict=True)}


def _by_pair(frame: pd.DataFrame) -> dict[str, pd.Series]:
    assert frame["pair_id"].is_unique, frame
    return {str(r["pair_id"]): r for _, r in frame.iterrows()}


def _ax(inst: Instance, label: str) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_label() == label]
    assert len(found) == 1, f"axes {label!r}: {[a.get_label() for a in inst.figure.axes]}"
    return found[0]


def _texts(ax: Axes) -> list[str]:
    return [t.get_text() for t in ax.findobj(Text) if t.get_text()]


def _legend_texts(ax: Axes) -> list[str]:
    leg = ax.get_legend()
    return [] if leg is None else [t.get_text() for t in leg.get_texts()]


def _drawn_colours(ax: Axes) -> set[str]:
    """Hex colours of every line and every patch edge/face drawn in ``ax`` (legend excluded)."""
    out: set[str] = set()
    for line in ax.get_lines():
        out.add(mcolors.to_hex(line.get_color()).upper())
    for patch in ax.patches:
        for rgba in (patch.get_edgecolor(), patch.get_facecolor()):
            if rgba[3] > 0:
                out.add(mcolors.to_hex(rgba[:3]).upper())
    return out


def _points(ax: Axes) -> list[tuple[float, float, float, float, PathCollection]]:
    """``(x, y, size, edge width, collection)`` of every scatter point in ``ax``."""
    pts = []
    for coll in ax.collections:
        if not isinstance(coll, PathCollection):
            continue
        offsets = np.asarray(coll.get_offsets(), dtype=float)
        sizes = np.broadcast_to(np.asarray(coll.get_sizes(), dtype=float), (len(offsets),))
        lws = np.asarray(coll.get_linewidths(), dtype=float)
        lws = np.broadcast_to(lws if lws.size else np.zeros(1), (len(offsets),))
        for (x, y), s, w in zip(offsets, sizes, lws, strict=True):
            pts.append((float(x), float(y), float(s), float(w), coll))
    return pts


def _outlined(ax: Axes) -> list[tuple[float, float, PathCollection]]:
    pts = _points(ax)
    if not pts:
        return []
    smallest = min(p[2] for p in pts)
    return [(x, y, c) for x, y, s, w, c in pts if s > smallest and w > 0]


def _point_colours(ax: Axes) -> set[str]:
    out: set[str] = set()
    for coll in ax.collections:
        if isinstance(coll, PathCollection):
            for rgba in np.atleast_2d(coll.get_facecolors()):
                if len(rgba) and rgba[3] > 0:
                    out.add(mcolors.to_hex(rgba[:3]).upper())
    return out


def _hlines(ax: Axes, linestyle: str) -> list[tuple[float, Line2D]]:
    found = []
    for line in ax.get_lines():
        y = [float(v) for v in line.get_ydata()]
        if line.get_linestyle() == linestyle and len(y) >= 2 and len(set(y)) == 1:
            found.append((y[0], line))
    return found


# --- registration and skips -------------------------------------------------------------------


def test_v11_is_registered_as_operon_pairs_needing_pairs() -> None:
    spec = _v11()
    assert spec.pid == "V11"
    assert spec.name == "operon_pairs"
    assert spec.needs == ("pairs",)
    assert "operon_flags" in spec.optional
    assert "pairs" not in spec.optional
    assert spec.module == "core"


def test_no_pairs_table_is_skipped_as_input_absent(tmp_path: Path, cfg: GippaConfig) -> None:
    out = _out(tmp_path, None, flags=[])
    assert _v11().fn(RunView.from_outdir(out), cfg) == Skip("input_absent:pairs")


def test_flag_markers_cover_the_seven_flags_with_distinct_shapes() -> None:
    import gippa.viz.operon_pairs as v11

    assert set(v11.FLAG_MARKERS) == set(FLAGS)
    assert len(set(v11.FLAG_MARKERS.values())) == 7


# --- structure, title, style ------------------------------------------------------------------


def test_one_instance_double_width_with_three_labelled_panels(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    assert inst.suffix is None
    assert inst.coords == ()
    assert round(float(inst.figure.get_size_inches()[0]), 6) == style.DOUBLE_IN
    assert {ax.get_label() for ax in inst.figure.axes} >= {"a", "b", "c"}
    assert list(inst.data.columns) == COLUMNS
    assert set(inst.data["panel"]) == {"a", "b", "c"}


def test_panel_c_axis_labels_name_upstream_and_downstream(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    c = _ax(inst, "c")
    assert "upstream" in c.get_xlabel().lower()
    assert "downstream" in c.get_ylabel().lower()
    for label in ("a", "b"):
        assert _ax(inst, label).get_xlabel()
        assert _ax(inst, label).get_ylabel()


def test_example1_title_states_medians_and_flag_count(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path, EXAMPLE1), cfg)
    expected = "Median |Δ| OP 0.4 vs NO 1.3 log2; 0 pairs flagged at q ≤ 0.1"
    assert inst.title == expected
    assert expected in [t.get_text() for t in inst.figure.findobj(Text)]
    assert "express" not in inst.title.lower()


def test_title_counts_flags_at_or_below_alpha(tmp_path: Path, cfg: GippaConfig) -> None:
    """Examples 4-5: PAIR_3 (q 0.04) counts, PAIR_4 (q 0.2) does not."""
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    assert inst.title.startswith("Median |Δ| OP ")
    assert inst.title.endswith("; 1 of 2 candidate flags at q ≤ 0.1")


def test_title_and_threshold_follow_the_configured_alpha(tmp_path: Path) -> None:
    cfg = resolve_config(overrides=["stats.alpha.operon_flags=0.05"])
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    assert inst.title.endswith("; 1 of 2 candidate flags at q ≤ 0.05")
    assert [y for y, _ in _hlines(_ax(inst, "b"), "--")] == [pytest.approx(0.05)]
    cfg_strict = resolve_config(overrides=["stats.alpha.operon_flags=0.01"])
    strict = _render(_out(tmp_path / "strict", pairs, _full_flags(pairs)), cfg_strict)
    assert strict.title.endswith("; 0 of 2 candidate flags at q ≤ 0.01")
    assert _outlined(_ax(strict, "c")) == []


def test_every_text_is_at_least_6_pt(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    sizes = [t.get_fontsize() for t in inst.figure.findobj(Text) if t.get_text()]
    assert sizes
    assert min(sizes) >= 6


def test_rows_are_sorted_by_panel_class_bin(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    data = inst.data
    panels = list(data["panel"])
    assert panels == sorted(panels)
    for panel, frame in data.groupby("panel", sort=False):
        ranks = [CLASS_ORDER.index(str(c)) for c in frame["class"]]
        assert ranks == sorted(ranks), panel
        for cls, sub in frame.groupby("class", sort=False):
            if panel == "a":
                bins = [_num(b) for b in sub["bin"] if not _isna(b)]
                assert bins == sorted(bins), cls
            else:  # single-digit IDs: string and number order agree
                ids = [str(p) for p in sub["pair_id"]]
                assert ids == sorted(ids), (panel, cls)


# --- panel (a): Δ densities -------------------------------------------------------------------


def test_example1_medians_bins_and_densities(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path, EXAMPLE1), cfg)
    width = 1.6 / N_BINS  # max |Δ| over the four classes
    for cls, n in (("OP", 3), ("NO", 2)):
        rows = _bin_rows(inst, cls)
        assert len(rows) == N_BINS
        bins = [int(_num(b)) for b in rows["bin"]]
        assert bins == list(range(bins[0], bins[0] + N_BINS))
        assert bins[0] in (0, 1)
        dens = [_num(d) for d in rows["density"]]
        assert sum(dens) * width == pytest.approx(1.0, rel=1e-5)
        assert all(_isna(p) for p in rows["pair_id"])
        nonzero = [b - bins[0] for b, d in zip(bins, dens, strict=True) if d > 0]
        expected = {"OP": [3, 7, 11], "NO": [18, 29]}[cls]
        assert nonzero == expected, cls
        for d in dens:
            assert d == pytest.approx(0.0) or d == pytest.approx(1 / (n * width), rel=1e-5)
    assert len(_bin_rows(inst)) == 2 * N_BINS  # CV and SW absent


def test_count_rows_for_all_five_classes(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path, EXAMPLE1), cfg)
    assert _counts(inst) == {"OP": 3, "NO": 2, "CV": 0, "SW": 0, "DV": 0}
    a = _rows(inst, "a")
    count_rows = a[[_isna(v) for v in a["bin"]]]
    for col in ("pair_id", "density", "abs_delta_log2"):
        assert all(_isna(v) for v in count_rows[col]), col


def test_example6_dv_is_counted_but_not_drawn(tmp_path: Path, cfg: GippaConfig) -> None:
    dv = [_pair(10 + i, "DV", 20.0, 20.0 + d) for i, d in enumerate([5.0, -4.0, 7.0])]
    with_dv = _render(_out(tmp_path / "dv", [*EXAMPLE1, *dv]), cfg)
    without = _render(_out(tmp_path / "no_dv", EXAMPLE1), cfg)
    assert _counts(with_dv)["DV"] == 3
    assert "DV" not in set(_bin_rows(with_dv)["class"])
    # shared bins over OP, NO, CV, SW only: DV's |Δ| = 7 does not widen them
    for cls in ("OP", "NO"):
        got = [_num(d) for d in _bin_rows(with_dv, cls)["density"]]
        ref = [_num(d) for d in _bin_rows(without, cls)["density"]]
        assert got == pytest.approx(ref)
    assert not any(t.startswith("DV") for t in _legend_texts(_ax(with_dv, "a")))


def test_shared_bins_span_the_largest_class(tmp_path: Path, cfg: GippaConfig) -> None:
    """CV's |Δ| = 3 sets the range: OP's |Δ| = 0.1 and 0.2 fall in the first two bins."""
    pairs = [
        *_delta_pairs("OP", [0.05, 0.15], 1),
        *_delta_pairs("NO", [1.5], 3),
        *_delta_pairs("CV", [3.0, -1.0], 4),
        *_delta_pairs("SW", [0.5], 6),
    ]
    inst = _render(_out(tmp_path, pairs), cfg)
    width = 3.0 / N_BINS
    bins_by_class = {}
    for cls in ("OP", "NO", "CV", "SW"):
        rows = _bin_rows(inst, cls)
        assert len(rows) == N_BINS, cls
        bins_by_class[cls] = [int(_num(b)) for b in rows["bin"]]
        dens = [_num(d) for d in rows["density"]]
        assert sum(dens) * width == pytest.approx(1.0, rel=1e-5), cls
    assert len({tuple(b) for b in bins_by_class.values()}) == 1
    op = [_num(d) for d in _bin_rows(inst, "OP")["density"]]
    assert op[0] == pytest.approx(0.5 / width, rel=1e-5)
    assert op[1] == pytest.approx(0.5 / width, rel=1e-5)
    assert sum(op[2:]) == pytest.approx(0.0)


def test_zero_max_delta_bins_over_unit_interval(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = [*_delta_pairs("OP", [0.0, 0.0], 1), *_delta_pairs("NO", [0.0], 3)]
    inst = _render(_out(tmp_path, pairs), cfg)
    width = 1.0 / N_BINS
    for cls in ("OP", "NO"):
        dens = [_num(d) for d in _bin_rows(inst, cls)["density"]]
        assert len(dens) == N_BINS
        assert dens[0] == pytest.approx(1 / width, rel=1e-5)
        assert sum(dens[1:]) == pytest.approx(0.0)
    assert inst.title.startswith("Median |Δ| OP 0 vs NO 0 log2;")


def test_single_pair_class_is_a_single_bin_histogram(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = [*_delta_pairs("OP", [0.3, 0.9, -1.2], 1), *_delta_pairs("NO", [0.6], 4)]
    inst = _render(_out(tmp_path, pairs), cfg)
    dens = [_num(d) for d in _bin_rows(inst, "NO")["density"]]
    assert sum(d > 0 for d in dens) == 1
    assert max(dens) == pytest.approx(N_BINS / 1.2, rel=1e-5)


def test_panel_a_draws_class_colours_and_a_legend_with_n(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    ax = _ax(inst, "a")
    assert set(CLASS_COLOUR.values()) <= _drawn_colours(ax)
    legend = _legend_texts(ax)
    counts = {"OP": 4, "NO": 2, "CV": 1, "SW": 1}
    for cls, n in counts.items():
        hits = [t for t in legend if t.startswith(cls)]
        assert len(hits) == 1, (cls, legend)
        assert str(n) in hits[0], (cls, hits)
    assert not any(t.startswith("DV") for t in legend)


def test_panel_a_omits_absent_classes(tmp_path: Path, cfg: GippaConfig) -> None:
    """Single-strand genome: no CV, no SW pairs."""
    inst = _render(_out(tmp_path, EXAMPLE1), cfg)
    ax = _ax(inst, "a")
    drawn = _drawn_colours(ax)
    assert {CLASS_COLOUR["OP"], CLASS_COLOUR["NO"]} <= drawn
    assert CLASS_COLOUR["CV"] not in drawn
    assert CLASS_COLOUR["SW"] not in drawn
    legend = _legend_texts(ax)
    assert not any(t.startswith(("CV", "SW")) for t in legend)


# --- panel (b): local fdr ---------------------------------------------------------------------


def test_panel_b_rows_are_the_pairs_with_an_lfdr(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    b = _by_pair(_rows(inst, "b"))
    assert set(b) == {"PAIR_5", "PAIR_6"}
    assert _num(b["PAIR_5"]["abs_delta_log2"]) == pytest.approx(1.0)
    assert _num(b["PAIR_5"]["lfdr"]) == pytest.approx(0.9)
    assert _num(b["PAIR_6"]["abs_delta_log2"]) == pytest.approx(0.2)
    assert _num(b["PAIR_6"]["lfdr"]) == pytest.approx(0.4)
    assert all(str(r["class"]) == "NO" for r in b.values())


def test_panel_b_curve_is_joined_in_x_order(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = [
        *_delta_pairs("OP", [0.1, 0.3], 1),
        _pair(3, "NO", 20.0, 21.5, lfdr=0.95),
        _pair(4, "NO", 20.0, 20.1, lfdr=0.2),
        _pair(5, "NO", 20.0, 19.2, lfdr=0.7),
    ]
    inst = _render(_out(tmp_path, pairs), cfg)
    ax = _ax(inst, "b")
    want_x, want_y = [0.1, 0.8, 1.5], [0.2, 0.7, 0.95]
    curves = [
        line
        for line in ax.get_lines()
        if len(line.get_xdata()) == 3
        and np.allclose(np.asarray(line.get_xdata(), dtype=float), want_x)
    ]
    assert len(curves) == 1, [list(line.get_xdata()) for line in ax.get_lines()]
    assert np.allclose(np.asarray(curves[0].get_ydata(), dtype=float), want_y)


def test_panel_b_threshold_line_and_note(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    ax = _ax(inst, "b")
    lines = _hlines(ax, "--")
    assert [y for y, _ in lines] == [pytest.approx(0.1)]
    texts = _texts(ax)
    assert lines[0][1].get_label() == "flag threshold (q)" or "flag threshold (q)" in texts
    assert "lfdr ranks pairs; flags use q" in texts


def test_panel_b_without_lfdr_says_so(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path, EXAMPLE1), cfg)
    assert len(_rows(inst, "b")) == 0
    assert "no lfdr values" in _texts(_ax(inst, "b"))


# --- panel (c): partner scatter ---------------------------------------------------------------


def test_example2_plus_strand_upstream_is_gene_a(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    row = _by_pair(_rows(inst, "c"))["PAIR_1"]
    assert _num(row["upstream_log2"]) == 20.0
    assert _num(row["downstream_log2"]) == 22.0


def test_upstream_is_gene_a_on_both_strands(tmp_path: Path, cfg: GippaConfig) -> None:
    """``pairs.tsv`` gene_a is already the upstream gene in transcription order (schema; M6-04
    swaps a ``-`` pair). Spec example 3 (upstream 22 for a ``-`` pair with a = 20, b = 22) would
    swap twice; this pins the schema reading pending the lead's decision."""
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    row = _by_pair(_rows(inst, "c"))["PAIR_2"]
    assert _num(row["upstream_log2"]) == 20.0
    assert _num(row["downstream_log2"]) == 22.0


def test_panel_c_rows_are_measured_same_strand_pairs(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = [*_full_pairs(), _pair(10, "OP", 20.0, None), _pair(11, "NO", None, None)]
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    c = _by_pair(_rows(inst, "c"))
    assert set(c) == {"PAIR_1", "PAIR_2", "PAIR_3", "PAIR_4", "PAIR_5", "PAIR_6"}
    assert {str(r["class"]) for r in c.values()} == {"OP", "NO"}
    assert _counts(inst)["OP"] == 4  # the half-measured OP pair has no finite Δ
    assert len(_points(_ax(inst, "c"))) >= 6


def test_example4_flagged_pair_is_outlined(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    c = _by_pair(_rows(inst, "c"))
    assert str(c["PAIR_3"]["flag"]) == "unexplained_discordance"
    assert _num(c["PAIR_3"]["q_value"]) == pytest.approx(0.04)
    outlined = _outlined(_ax(inst, "c"))
    assert [(x, y) for x, y, _ in outlined] == [(18.0, 21.0)]


def test_example5_flag_above_alpha_is_not_outlined(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    c = _by_pair(_rows(inst, "c"))
    assert _num(c["PAIR_4"]["q_value"]) == pytest.approx(0.2)
    assert (19.0, 19.5) not in [(x, y) for x, y, _ in _outlined(_ax(inst, "c"))]
    assert all(_isna(c[p]["flag"]) for p in ("PAIR_1", "PAIR_2", "PAIR_5", "PAIR_6"))


def test_different_flags_get_different_marker_shapes(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = _full_pairs()
    flags = [
        _flag_row(pairs[2], "unexplained_discordance", 0.04),
        _flag_row(pairs[0], "internal_promoter_supported", 0.01),
    ]
    inst = _render(_out(tmp_path, pairs, flags), cfg)
    outlined = _outlined(_ax(inst, "c"))
    assert sorted((x, y) for x, y, _ in outlined) == [(18.0, 21.0), (20.0, 22.0)]
    paths = {}
    for x, y, coll in outlined:
        idx = 0 if len(coll.get_paths()) == 1 else None
        assert idx is not None, "one marker path per flag collection"
        paths[(x, y)] = np.round(coll.get_paths()[0].vertices, 6).tobytes()
    assert paths[(18.0, 21.0)] != paths[(20.0, 22.0)]


def test_no_flags_table_draws_no_flag_markers(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, flags=None), cfg)
    assert _outlined(_ax(inst, "c")) == []
    assert all(_isna(v) for v in _rows(inst, "c")["flag"])
    assert len(_rows(inst, "c")) == 6


def test_panel_c_colours_by_class_and_draws_identity(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    ax = _ax(inst, "c")
    assert {CLASS_COLOUR["OP"], CLASS_COLOUR["NO"]} <= _point_colours(ax)
    identity = [
        line
        for line in ax.get_lines()
        if len(line.get_xdata()) >= 2
        and np.allclose(
            np.asarray(line.get_xdata(), dtype=float), np.asarray(line.get_ydata(), dtype=float)
        )
        and len(set(np.asarray(line.get_xdata(), dtype=float))) >= 2
    ]
    assert identity


# --- degenerate cases -------------------------------------------------------------------------


def test_validation_like_unmeasured_pairs_empty_states(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = [_pair(1, "OP"), _pair(2, "NO"), _pair(3, "CV"), _pair(4, "DV"), _pair(5, "SW")]
    inst = _render(_out(tmp_path, pairs, flags=[]), cfg)
    assert inst.title == "0 pairs flagged at q ≤ 0.1; too few measured pairs for densities"
    assert len(_bin_rows(inst)) == 0
    assert _counts(inst) == {c: 0 for c in CLASS_ORDER}
    assert len(_rows(inst, "b")) == 0
    assert len(_rows(inst, "c")) == 0
    assert _texts(_ax(inst, "a"))
    assert _texts(_ax(inst, "c"))
    assert "no lfdr values" in _texts(_ax(inst, "b"))
    assert _points(_ax(inst, "c")) == []


def test_no_no_pairs_uses_the_too_few_title(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = [*_delta_pairs("OP", [0.2, 0.5], 1), *_delta_pairs("CV", [1.0], 3)]
    pairs_flags = [_flag_row(pairs[0], "unexplained_discordance", 0.01)]
    inst = _render(_out(tmp_path, pairs, pairs_flags), cfg)
    assert inst.title == "1 pairs flagged at q ≤ 0.1; too few measured pairs for densities"


def test_header_only_pairs_is_the_framework_empty_state(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path, [], flags=[]), cfg)
    assert len(inst.data) == 0
    assert list(inst.data.columns) == COLUMNS
    assert inst.title


# --- an M6 fixture with flags (the real flag_pairs) ---------------------------------------------


def test_renders_on_an_m6_flag_fixture(tmp_path: Path, cfg: GippaConfig) -> None:
    from gippa.io.operon_tables import write_operon_flags_tsv, write_pairs_tsv
    from gippa.operons.flags import flag_pairs

    r = flags_data.row
    pairs = flags_data.frame(
        r("OP", delta=2.0, q=0.01, prom=0.99),  # internal_promoter_supported
        r("OP", delta=3.0, q=0.02),  # unexplained_discordance
        r("OP", delta=0.1, q=0.5, term=0.99),  # conflict_terminator_in_operon, q > alpha
        r("OP", delta=0.3),
        r("NO", delta=0.05, q=0.05, lfdr=0.3),  # novel_operon_candidate
        r("NO", delta=1.2, q=0.8, lfdr=0.95),
        r("CV", delta=1.0),
        r("SW", delta=-0.5),
        r("DV", delta=0.7),
    )
    flagged = flag_pairs(pairs, flags_data.cfg())
    assert len(flagged) == 4
    out = make_out(tmp_path)
    write_pairs_tsv(out / "tables" / "pairs.tsv", pairs)
    write_operon_flags_tsv(out / "tables" / "operon_flags.tsv", flagged)
    inst = _render(out, cfg)
    assert inst.title.endswith("; 3 of 4 candidate flags at q ≤ 0.1")
    c = _by_pair(_rows(inst, "c"))
    for _, f in flagged.iterrows():
        assert str(c[str(f["pair_id"])]["flag"]) == f["flag"]
    assert len(_outlined(_ax(inst, "c"))) == 3
    assert set(_by_pair(_rows(inst, "b"))) == {"PAIR_5", "PAIR_6"}
    assert _counts(inst) == {"OP": 4, "NO": 2, "CV": 1, "SW": 1, "DV": 1}


# --- properties -------------------------------------------------------------------------------

_DELTA = st.floats(-6.0, 6.0, allow_nan=False).map(lambda v: round(v, 3))


@settings(max_examples=12, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    deltas=st.fixed_dictionaries(
        {c: st.lists(_DELTA, max_size=6) for c in ("OP", "NO", "CV", "SW", "DV")}
    ).filter(lambda d: len(d["OP"]) + len(d["NO"]) + len(d["CV"]) + len(d["SW"]) > 0)
)
def test_property_shared_bins_integrate_to_one(deltas: dict[str, list[float]]) -> None:
    style.apply()
    pairs: list[dict[str, object]] = []
    for cls in CLASS_ORDER:
        pairs += _delta_pairs(cls, deltas[cls], len(pairs) + 1)
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(_out(Path(tmp), pairs), resolve_config())
    # |Δ| as the TSV writer stores it (b - a on 20 + d)
    absd = {c: [abs((20.0 + d) - 20.0) for d in deltas[c]] for c in CLASS_ORDER}
    top = max(v for c in ("OP", "NO", "CV", "SW") for v in absd[c])
    width = (top if top > 0 else 1.0) / N_BINS
    assert _counts(inst) == {c: len(deltas[c]) for c in CLASS_ORDER}
    seen = set()
    for cls in ("OP", "NO", "CV", "SW"):
        rows = _bin_rows(inst, cls)
        if not deltas[cls]:
            assert len(rows) == 0, cls
            continue
        assert len(rows) == N_BINS, cls
        seen.add(tuple(int(_num(b)) for b in rows["bin"]))
        dens = [_num(d) for d in rows["density"]]
        assert all(d >= 0 for d in dens)
        assert sum(dens) * width == pytest.approx(1.0, rel=1e-4), cls
    assert len(seen) == 1
    assert "DV" not in set(_bin_rows(inst)["class"])


@settings(max_examples=12, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    specs=st.lists(
        st.tuples(
            st.sampled_from(CLASS_ORDER),
            st.sampled_from(["+", "-"]),
            st.one_of(st.none(), st.integers(10, 30).map(float)),
            st.one_of(st.none(), st.integers(10, 30).map(float)),
            st.one_of(st.none(), st.floats(0.0, 1.0).map(lambda v: round(v, 3))),
        ),
        min_size=1,
        max_size=9,
    )
)
def test_property_panel_b_and_c_row_sets(
    specs: list[tuple[str, str, float | None, float | None, float | None]],
) -> None:
    style.apply()
    pairs = []
    for i, (cls, strand, a, b, lfdr) in enumerate(specs, start=1):
        measured = a is not None and b is not None
        pairs.append(
            _pair(
                i,
                cls,
                a,
                b,
                strand=strand if cls in ("OP", "NO", "SW") else ".",
                lfdr=lfdr if cls == "NO" and measured else None,
            )
        )
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(_out(Path(tmp), pairs), resolve_config())
    c = _by_pair(_rows(inst, "c")) if len(_rows(inst, "c")) else {}
    want_c = {
        str(p["pair_id"]): p
        for p in pairs
        if p["class"] in ("OP", "NO")
        and p["abundance_log2_a"] is not None
        and p["abundance_log2_b"] is not None
    }
    assert set(c) == set(want_c)
    for pid, p in want_c.items():
        assert _num(c[pid]["upstream_log2"]) == p["abundance_log2_a"]
        assert _num(c[pid]["downstream_log2"]) == p["abundance_log2_b"]
    b_rows = _by_pair(_rows(inst, "b")) if len(_rows(inst, "b")) else {}
    want_b = {str(p["pair_id"]): p for p in pairs if p["lfdr"] is not None}
    assert set(b_rows) == set(want_b)
    for pid, p in want_b.items():
        assert _num(b_rows[pid]["lfdr"]) == pytest.approx(float(p["lfdr"]))  # type: ignore[arg-type]
        assert _num(b_rows[pid]["abs_delta_log2"]) == pytest.approx(abs(float(p["delta_log2"])))  # type: ignore[arg-type]


# --- files and determinism --------------------------------------------------------------------


def test_data_tsv_and_images_are_byte_identical(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = _full_pairs()
    out = _out(tmp_path / "run", pairs, _full_flags(pairs))
    written = []
    for name in ("x", "y"):
        inst = _render(out, cfg)
        paths = write_instance(
            inst, "V11", "operon_pairs", tmp_path / name, ["svg", "png", "pdf"], 100
        )
        written.append({p.relative_to(tmp_path / name): p.read_bytes() for p in paths})
    assert set(written[0]) == {
        Path("V11_operon_pairs.svg"),
        Path("V11_operon_pairs.png"),
        Path("V11_operon_pairs.pdf"),
        Path("data/V11_operon_pairs.tsv"),
    }
    assert written[0] == written[1]
    tsv = written[0][Path("data/V11_operon_pairs.tsv")].decode("utf-8").splitlines()
    assert tsv[0] == "\t".join(COLUMNS)
    assert "c\tPAIR_3\tOP\t" in "\n".join(tsv)


def test_render_does_not_mutate_between_calls(tmp_path: Path, cfg: GippaConfig) -> None:
    pairs = _full_pairs()
    view = RunView.from_outdir(_out(tmp_path, pairs, _full_flags(pairs)))
    first = _v11().fn(view, cfg)
    second = _v11().fn(view, cfg)
    assert isinstance(first, list) and isinstance(second, list)
    pd.testing.assert_frame_equal(first[0].data, second[0].data)
    assert first[0].title == second[0].title


# --- the validation run (local/, slow) --------------------------------------------------------


@pytest.mark.slow
@pytest.mark.skipif(
    not (VAL_OUT / "tables" / "pairs.tsv").is_file(), reason=f"validation run absent: {VAL_OUT}"
)
def test_renders_on_the_validation_run(cfg: GippaConfig) -> None:
    """87 genes; every pair untested, ``operon_flags.tsv`` header-only; most pairs have a finite
    ``delta_log2`` (OP 65, NO 1, CV 8, SW 5, DV 8 on 2026-10-22)."""
    inst = _render(VAL_OUT, cfg)
    assert list(inst.data.columns) == COLUMNS
    assert inst.title.endswith("0 pairs flagged at q ≤ 0.1")
    assert _counts(inst) == {"OP": 65, "NO": 1, "CV": 8, "SW": 5, "DV": 8}
    assert len(_rows(inst, "b")) == 0
    assert "no lfdr values" in _texts(_ax(inst, "b"))
    assert _outlined(_ax(inst, "c")) == []
    assert len(_rows(inst, "c")) == 66
