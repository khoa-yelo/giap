"""V11 ``operon_pairs``, PR #208 review round (M9-11 lead decisions 3 and 6; review items B1-B5).

- B1: panel (b)'s legend, the "lfdr ranks pairs; flags use q" note and the dashed alpha line do
  not overlap, measured on the Agg canvas after a draw (``FigureCanvasBase`` gave false passes).
- B2 (decision 3): pairs whose ``pairs.tsv`` ``flags`` (semicolon-joined pair-scope tokens,
  schema ``pairs.flags``; M6-08 writes ``abundance_imputed``) carry ``abundance_imputed`` are left
  out of panels (a) and (c), the count rows and the title medians, but stay in panel (b).
- B3: the title medians are medians (not means) printed with ``%.2g``.
- B4: ``q == alpha`` counts in the title ``k`` and is outlined in panel (c).
- B5 (decision 6): bins are 0-based, pair rows keep ``pairs.tsv`` order, the empty-state text.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.text import Text
from matplotlib.transforms import Bbox
from test_v11_operon_pairs import (  # type: ignore[import-not-found]
    EXAMPLE1,
    N_BINS,
    _ax,
    _bin_rows,
    _by_pair,
    _counts,
    _delta_pairs,
    _flag_row,
    _full_flags,
    _full_pairs,
    _hlines,
    _num,
    _out,
    _outlined,
    _pair,
    _points,
    _render,
    _rows,
)

from giap.config import GiapConfig, resolve_config
from giap.viz import style

NOTE_B = "lfdr ranks pairs; flags use q"
IMPUTED = "abundance_imputed"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


def _with_flags(pair: dict[str, object], flags: str | None) -> dict[str, object]:
    return {**pair, "flags": flags}


# --- B1: panel (b) layout -----------------------------------------------------------------------


def _overlaps(a: Bbox, b: Bbox) -> bool:
    """Strict interior overlap of two display boxes (touching edges do not count)."""
    return a.x0 < b.x1 and b.x0 < a.x1 and a.y0 < b.y1 and b.y0 < a.y1


@pytest.mark.parametrize("alpha", [0.1, 0.05])
def test_b1_panel_b_legend_note_and_threshold_line_do_not_overlap(
    tmp_path: Path, alpha: float
) -> None:
    cfg = resolve_config(overrides=[f"stats.alpha.operon_flags={alpha}"])
    pairs = _full_pairs()
    inst = _render(_out(tmp_path, pairs, _full_flags(pairs)), cfg)
    fig = inst.figure
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    renderer = canvas.get_renderer()
    ax: Axes = _ax(inst, "b")

    leg = ax.get_legend()
    assert leg is not None, "panel (b) has no legend for the dashed threshold line"
    leg_box = leg.get_window_extent(renderer)
    notes = [t for t in ax.findobj(Text) if t.get_text() == NOTE_B]
    assert len(notes) == 1
    note_box = notes[0].get_window_extent(renderer)

    lines = _hlines(ax, "--")
    assert [y for y, _ in lines] == [pytest.approx(alpha)]
    line = lines[0][1]
    y_disp = float(ax.transData.transform((0.0, alpha))[1])
    half = float(line.get_linewidth()) * fig.dpi / 72.0 / 2.0 + 1.0
    ax_box = ax.get_window_extent(renderer)
    line_box = Bbox.from_extents(ax_box.x0, y_disp - half, ax_box.x1, y_disp + half)

    assert not _overlaps(leg_box, note_box), (leg_box, note_box)
    assert not _overlaps(leg_box, line_box), ("legend on the dashed line", leg_box, line_box)
    assert not _overlaps(note_box, line_box), ("note on the dashed line", note_box, line_box)
    # the legend's label and handle each, too (the window extent covers both, but be explicit)
    for t in leg.get_texts():
        box = t.get_window_extent(renderer)
        assert not _overlaps(box, note_box) and not _overlaps(box, line_box), t.get_text()
    for h in leg.legend_handles:
        assert h is not None
        box = h.get_window_extent(renderer)
        assert not _overlaps(box, note_box) and not _overlaps(box, line_box), h


# --- B2: imputed pairs (lead decision 3) ----------------------------------------------------------


def _imputed_case() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Example 1 plus a non-imputed ``rna_between`` OP pair (kept) and two imputed pairs (one OP
    with a large Δ that would move the bins and the median, one NO with an lfdr)."""
    kept = [_with_flags(p, None) for p in EXAMPLE1]
    kept.append(_with_flags(_pair(6, "OP", 20.0, 20.3), "rna_between"))
    imputed = [
        _with_flags(_pair(7, "OP", 10.0, 15.0), IMPUTED),
        _with_flags(_pair(8, "NO", 12.0, 16.0, lfdr=0.3), f"rna_between;{IMPUTED}"),
    ]
    return kept, imputed


def test_b2_imputed_pairs_leave_panel_a_counts_densities_and_medians(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    kept, imputed = _imputed_case()
    base = _render(_out(tmp_path / "base", kept), cfg)
    inst = _render(_out(tmp_path / "imp", [*kept, *imputed]), cfg)
    assert _counts(inst) == {"OP": 4, "NO": 2, "CV": 0, "SW": 0, "DV": 0}
    assert _counts(inst) == _counts(base)
    for cls in ("OP", "NO"):
        got = [_num(d) for d in _bin_rows(inst, cls)["density"]]
        want = [_num(d) for d in _bin_rows(base, cls)["density"]]
        assert got == pytest.approx(want), cls
    assert inst.title.startswith("Median |Δ| OP ")
    assert inst.title == base.title


def test_b2_imputed_pairs_leave_panel_c_rows_and_points(tmp_path: Path, cfg: GiapConfig) -> None:
    kept, imputed = _imputed_case()
    inst = _render(_out(tmp_path, [*kept, *imputed]), cfg)
    c = _by_pair(_rows(inst, "c"))
    assert set(c) == {str(p["pair_id"]) for p in kept}
    drawn = {(x, y) for x, y, _, _, _ in _points(_ax(inst, "c"))}
    assert (10.0, 15.0) not in drawn
    assert (12.0, 16.0) not in drawn
    assert (20.0, 20.3) in drawn  # rna_between alone is not imputed


def test_b2_imputed_pair_with_lfdr_stays_in_panel_b(tmp_path: Path, cfg: GiapConfig) -> None:
    kept, imputed = _imputed_case()
    inst = _render(_out(tmp_path, [*kept, *imputed]), cfg)
    b = _by_pair(_rows(inst, "b"))
    assert set(b) == {"PAIR_8"}
    assert _num(b["PAIR_8"]["lfdr"]) == pytest.approx(0.3)
    assert _num(b["PAIR_8"]["abs_delta_log2"]) == pytest.approx(4.0)


_ABUND = st.integers(10, 30).map(float)


@settings(max_examples=10, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    specs=st.lists(
        st.tuples(st.sampled_from(["OP", "NO"]), _ABUND, _ABUND, st.booleans()),
        min_size=1,
        max_size=8,
    )
)
def test_b2_property_adding_imputed_pairs_changes_neither_a_nor_c(
    specs: list[tuple[str, float, float, bool]],
) -> None:
    style.apply()
    pairs = [
        _with_flags(_pair(i, cls, a, b), IMPUTED if imp else None)
        for i, (cls, a, b, imp) in enumerate(specs, start=1)
    ]
    kept = [p for p in pairs if p["flags"] is None]
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(_out(Path(tmp) / "all", pairs), resolve_config())
        want = _render(_out(Path(tmp) / "kept", kept), resolve_config()) if kept else None
    want_ids = {str(p["pair_id"]) for p in kept}
    got_c = set(_by_pair(_rows(inst, "c"))) if len(_rows(inst, "c")) else set()
    assert got_c == want_ids
    for cls in ("OP", "NO"):
        assert _counts(inst)[cls] == sum(1 for p in kept if p["class"] == cls)
    if want is not None:
        for cls in ("OP", "NO"):
            got = [_num(d) for d in _bin_rows(inst, cls)["density"]]
            exp = [_num(d) for d in _bin_rows(want, cls)["density"]]
            assert got == pytest.approx(exp), cls
        assert inst.title == want.title


# --- B3: title medians --------------------------------------------------------------------------


def test_b3_title_medians_are_medians_printed_with_2g(tmp_path: Path, cfg: GiapConfig) -> None:
    """OP |Δ| 0.1, 0.123, 0.9 (median 0.123 -> "0.12"; mean 0.374 -> "0.37");
    NO |Δ| 1.0, 1.04, 2.5 (median 1.04 -> "1", %.3g "1.04"; mean 1.51 -> "1.5")."""
    pairs = [*_delta_pairs("OP", [0.1, -0.123, 0.9], 1), *_delta_pairs("NO", [1.0, 1.04, -2.5], 4)]
    inst = _render(_out(tmp_path, pairs), cfg)
    assert inst.title == "Median |Δ| OP 0.12 vs NO 1 log2; 0 pairs flagged at q ≤ 0.1"


# --- B4: q exactly alpha ------------------------------------------------------------------------


@pytest.mark.parametrize("alpha", [0.1, 0.05])
def test_b4_q_equal_to_alpha_is_counted_and_outlined(tmp_path: Path, alpha: float) -> None:
    cfg = resolve_config(overrides=[f"stats.alpha.operon_flags={alpha}"])
    pairs = _full_pairs()
    by_id = {p["pair_id"]: p for p in pairs}
    flags = [
        _flag_row(by_id["PAIR_3"], "unexplained_discordance", alpha),
        _flag_row(by_id["PAIR_4"], "conflict_terminator_in_operon", 0.2),
    ]
    inst = _render(_out(tmp_path, pairs, flags), cfg)
    assert inst.title.endswith(f"; 1 of 2 candidate flags at q ≤ {alpha:g}")
    assert [(x, y) for x, y, _ in _outlined(_ax(inst, "c"))] == [(18.0, 21.0)]


# --- B5: lead decision 6 ------------------------------------------------------------------------


def test_b5_bins_are_zero_based_and_the_smallest_delta_is_bin_0(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    pairs = [*_delta_pairs("OP", [0.0, 0.4, 0.6], 1), *_delta_pairs("NO", [1.0, -1.6], 4)]
    inst = _render(_out(tmp_path, pairs), cfg)
    for cls in ("OP", "NO"):
        rows = _bin_rows(inst, cls)
        assert [int(_num(b)) for b in rows["bin"]] == list(range(N_BINS)), cls
    op = {
        int(_num(b)): _num(d)
        for b, d in zip(*(_bin_rows(inst, "OP")[k] for k in ("bin", "density")), strict=True)
    }
    no = {
        int(_num(b)): _num(d)
        for b, d in zip(*(_bin_rows(inst, "NO")[k] for k in ("bin", "density")), strict=True)
    }
    assert op[0] > 0  # |Δ| = 0 falls in bin 0
    assert no[N_BINS - 1] > 0  # the largest |Δ| (1.6) falls in the last bin
    assert sorted(b for b, d in op.items() if d > 0) == [0, 7, 11]
    # Example 1 (smallest |Δ| 0.2 = bin 3 on 1.6 / 30 wide bins)
    ex1 = _render(_out(tmp_path / "ex1", EXAMPLE1), cfg)
    rows = _bin_rows(ex1, "OP")
    assert int(_num(rows["bin"][0])) == 0
    assert [
        int(_num(b)) for b, d in zip(rows["bin"], rows["density"], strict=True) if _num(d) > 0
    ] == [3, 7, 11]


def test_b5_pair_rows_keep_pairs_tsv_order_with_multi_digit_ids(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    pairs = [
        _pair(2, "NO", 20.0, 21.0, lfdr=0.5),
        _pair(10, "OP", 20.0, 20.5),
        _pair(1, "NO", 20.0, 22.0, lfdr=0.2),
        _pair(3, "OP", 20.0, 20.1),
        _pair(11, "NO", 20.0, 23.0, lfdr=0.1),
    ]
    inst = _render(_out(tmp_path, pairs), cfg)
    b = _rows(inst, "b")
    assert list(b["pair_id"]) == ["PAIR_2", "PAIR_1", "PAIR_11"]
    c = _rows(inst, "c")
    assert list(c["pair_id"]) == ["PAIR_10", "PAIR_3", "PAIR_2", "PAIR_1", "PAIR_11"]
    assert list(c["class"]) == ["OP", "OP", "NO", "NO", "NO"]


def test_b5_header_only_pairs_gives_the_exact_empty_state_text(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    inst = _render(_out(tmp_path, [], flags=[]), cfg)
    assert inst.title == "No gene pairs in this run"
    texts = [t.get_text() for t in inst.figure.findobj(Text) if t.get_text()]
    assert "No gene pairs in this run" in texts
    assert len(inst.data) == 0
