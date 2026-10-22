"""V6 ``mapping_fdr``: can I trust tolerant hits and de novo peptides? (M9-05; GIAP_SPEC.md
§10.2 V6, §9.3 ``mapping_fdr.tsv``, §5.2 ``stats.alpha.placement``, ``mapping.length_classes``).

V6 draws the M3 estimates; it computes none. Tests the data frame (=
``plots/data/V06_mapping_fdr.tsv``), the three panels, the alpha line, the title, the degenerate
cases and determinism. No ``--mpl`` baselines (M9-25).

Interface assumed where the spec leaves it open:

- The plot lives in ``giap.viz.mapping_fdr``; importing it registers ``V06`` (name
  ``mapping_fdr``, ``needs=("mapping_fdr",)``, module ``core``). It returns one double-width
  ``Instance``, ``suffix=None``, ``coords=()``.
- Panel axes carry ``Axes.get_label()`` ``"a"``, ``"b"`` and (only with de novo rows) ``"c"``.
- Data columns are exactly ``panel, length_class, n_diff, score_bin, metric, value, ci_low,
  ci_high, flags``: the spec's eight plus ``score_bin`` (panel (c) rows need their bin; the spec
  lists no column for it). Panel (a): ``metric`` ∈ {``decoys``, ``targets``}, one row per length
  class x n_diff x metric (42 with the default classes). Panel (b): ``length_class = all``,
  ``metric`` ∈ {``fdr_stratum``, ``n_accepted``} per n_diff (the title's accepted count comes from
  the data, index decision 5); NA strata keep ``value`` NA. Panel (c): ``metric = precision``,
  ``length_class`` and ``n_diff`` NA, ``flags`` as in ``mapping_fdr``.
- Rows sort by panel, length class (config order, ``all`` for (b)), n_diff, score_bin, metric
  (alphabetical).
- ``mapping_fdr`` carries one row set per decoy type; a decoy type whose rows are all zero (its
  peptide source is absent, as in the validation run) does not change V6.
- The title lists the strata with finite ``fdr_stratum`` ≤ alpha joined by ``", "``; alpha is
  printed with ``%g``; the NA clause is "stratum <m> unestimated (no targets)".
"""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

import matplotlib.colors as mcolors
import pandas as pd
import pytest
from helpers.viz_mapping import (
    EXAMPLE4_TOTALS,
    LENGTH_LABELS,
    STRATA,
    denovo_row,
    fdr_row,
    full_placement,
    mapping_out,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.collections import PathCollection
from matplotlib.patches import Rectangle
from matplotlib.text import Text

from giap.config import GiapConfig, resolve_config
from giap.io.plot_files import write_instance
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.base import Instance, Skip

REPO = Path(__file__).resolve().parents[3]
GOLDEN = REPO / "tests" / "data" / "golden" / "sim50k" / "mapping"
COLUMNS = [
    "panel",
    "length_class",
    "n_diff",
    "score_bin",
    "metric",
    "value",
    "ci_low",
    "ci_high",
    "flags",
]
EX4_TITLE = (
    "Placement FDR ≤ 0.01 in strata 0; 920 placements accepted; stratum 2 unestimated (no targets)"
)
#: per-class counts consistent with example 4 (example 5 is class 7, m = 0)
EX4_CLASSES = {("7", 0): (5, 0), ("12-14", 0): (895, 3), ("15-19", 1): (120, 2)}
EX6_DENOVO = [
    denovo_row(1, 0.8, 0.7, 0.88, "c_assumed"),
    denovo_row(2, 0.9, 0.82, 0.95, "c_assumed"),
    denovo_row(3, 0.95, 0.9, 0.98, "c_assumed"),
]
FOOTNOTE = "c assumed 0.99 (no reference peptides)"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


def _v06() -> base.PlotSpec:
    import giap.viz.mapping_fdr  # noqa: F401  (registers V06)

    return base.registered()["V06"]


def _render(out: Path, cfg: GiapConfig) -> Instance:
    result = _v06().fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert isinstance(inst, Instance)
    return inst


def _ex4_rows() -> list[dict]:
    return full_placement(EXAMPLE4_TOTALS, EX4_CLASSES)


def _ex4(tmp: Path, denovo: list[dict] | None = None) -> Path:
    return mapping_out(tmp, fdr=[*_ex4_rows(), *(denovo or [])])


def _rows(inst: Instance, panel: str) -> pd.DataFrame:
    return inst.data[inst.data["panel"] == panel].reset_index(drop=True)


def _isna(value: object) -> bool:
    return value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value))


def _num(value: object) -> float:
    return math.nan if _isna(value) else float(str(value))


def _cell(frame: pd.DataFrame, length_class: str, n_diff: int, metric: str) -> pd.Series:
    hit = frame[
        (frame["length_class"].astype(str) == length_class)
        & (frame["n_diff"].map(_num) == n_diff)
        & (frame["metric"].astype(str) == metric)
    ]
    assert len(hit) == 1, f"expected one row {length_class}/{n_diff}/{metric}:\n{frame}"
    return hit.iloc[0]


def _ax(inst: Instance, label: str) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_label() == label]
    assert len(found) == 1, f"axes {label!r}: {[a.get_label() for a in inst.figure.axes]}"
    return found[0]


def _notes(ax: Axes) -> list[str]:
    """Texts placed with ``ax.text`` / ``ax.annotate`` (not tick labels)."""
    return [t.get_text() for t in ax.texts if t.get_text()]


def _fig_texts(inst: Instance) -> list[str]:
    return [t.get_text() for t in inst.figure.findobj(Text) if t.get_text()]


def _hex(colour: object) -> str:
    return mcolors.to_hex(colour).upper()  # type: ignore[arg-type]


def _bars(ax: Axes) -> list[Rectangle]:
    return [p for p in ax.patches if isinstance(p, Rectangle)]


def _hlines(ax: Axes, linestyle: str) -> list[float]:
    ys = []
    for line in ax.get_lines():
        y = [float(v) for v in line.get_ydata()]
        if line.get_linestyle() == linestyle and len(y) >= 2 and len(set(y)) == 1:
            ys.append(y[0])
    return sorted(ys)


def _markers(ax: Axes) -> list[tuple[float, float, bool]]:
    """(x, y, hollow) of every marker drawn in ``ax`` (lines with a marker, scatter offsets)."""
    out = []
    for line in ax.get_lines():
        if line.get_marker() in (None, "", "None", " ", "none"):
            continue
        face = line.get_markerfacecolor()
        hollow = isinstance(face, str) and face.lower() == "none"
        for x, y in zip(line.get_xdata(), line.get_ydata(), strict=True):
            out.append((float(x), float(y), hollow))
    for coll in ax.collections:
        if not isinstance(coll, PathCollection):
            continue
        faces = coll.get_facecolors()
        for i, (x, y) in enumerate(coll.get_offsets()):
            hollow = len(faces) == 0 or faces[i % len(faces)][3] == 0
            out.append((float(x), float(y), hollow))
    return out


# --- registration and skips -------------------------------------------------------------------


def test_v06_is_registered_as_mapping_fdr() -> None:
    spec = _v06()
    assert spec.pid == "V06"
    assert spec.name == "mapping_fdr"
    assert spec.needs == ("mapping_fdr",)
    assert spec.module == "core"


def test_no_mapping_fdr_is_skipped_as_input_absent(tmp_path: Path, cfg: GiapConfig) -> None:
    out = mapping_out(tmp_path, hits=[], unmapped_rows=[])
    assert _v06().fn(RunView.from_outdir(out), cfg) == Skip("input_absent:mapping_fdr")


# --- structure --------------------------------------------------------------------------------


def test_double_width_two_panels_without_denovo(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_ex4(tmp_path), cfg)
    assert inst.suffix is None
    assert inst.coords == ()
    assert round(float(inst.figure.get_size_inches()[0]), 6) == style.DOUBLE_IN
    labels = {ax.get_label() for ax in inst.figure.axes}
    assert {"a", "b"} <= labels
    assert "c" not in labels
    assert list(inst.data.columns) == COLUMNS
    assert set(inst.data["panel"]) == {"a", "b"}


def test_rows_are_sorted_by_panel_class_ndiff_bin_metric(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_ex4(tmp_path, EX6_DENOVO), cfg)
    a = _rows(inst, "a")
    expected_a = [
        (label, m, metric)
        for label in LENGTH_LABELS
        for m in STRATA
        for metric in ("decoys", "targets")
    ]
    got_a = list(
        zip(a["length_class"].astype(str), a["n_diff"].map(_num), a["metric"], strict=True)
    )
    assert got_a == expected_a
    b = _rows(inst, "b")
    expected_b = [("all", m, metric) for m in STRATA for metric in ("fdr_stratum", "n_accepted")]
    got_b = list(
        zip(b["length_class"].astype(str), b["n_diff"].map(_num), b["metric"], strict=True)
    )
    assert got_b == expected_b
    assert list(inst.data["panel"]) == ["a"] * 42 + ["b"] * 6 + ["c"] * 3
    c = _rows(inst, "c")
    assert [int(_num(v)) for v in c["score_bin"]] == [1, 2, 3]


def test_every_text_is_at_least_6_pt(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_ex4(tmp_path, EX6_DENOVO), cfg)
    sizes = [t.get_fontsize() for t in inst.figure.findobj(Text) if t.get_text()]
    assert sizes
    assert min(sizes) >= 6


# --- panel (a) --------------------------------------------------------------------------------


def test_example5_targets_and_decoys_rows(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_ex4(tmp_path), cfg)
    a = _rows(inst, "a")
    assert len(a) == 7 * 3 * 2
    assert _num(_cell(a, "7", 0, "targets")["value"]) == 5
    assert _num(_cell(a, "7", 0, "decoys")["value"]) == 0
    assert _num(_cell(a, "12-14", 0, "decoys")["value"]) == 3
    assert _num(_cell(a, "15-19", 1, "targets")["value"]) == 120
    assert _num(_cell(a, "20-40", 2, "targets")["value"]) == 0
    assert not (a["length_class"].astype(str) == "all").any()
    assert all(_isna(v) for v in a["ci_low"])
    assert all(_isna(v) for v in a["score_bin"])


def test_panel_a_log_scale_and_colours(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_ex4(tmp_path), cfg)
    ax = _ax(inst, "a")
    assert ax.get_yscale() == "log"
    colours = [_hex(r.get_facecolor()) for r in _bars(ax)]
    assert {"#0072B2", "#D55E00"} <= set(colours)
    assert colours.count("#0072B2") <= 21
    assert colours.count("#D55E00") <= 21
    tops = sorted(
        r.get_y() + r.get_height() for r in _bars(ax) if _hex(r.get_facecolor()) == "#0072B2"
    )
    assert tops[-1] == pytest.approx(895)


def test_panel_a_follows_configured_length_classes(tmp_path: Path) -> None:
    cfg = resolve_config(overrides=["mapping.length_classes=[[7, 9], [10, 40]]"])
    rows = []
    for m in STRATA:
        rows += [
            fdr_row("7-9", m, 2 + m, 0, 0.0, 2),
            fdr_row("10-40", m, 10, 1, 0.1, 5),
            fdr_row("all", m, 12 + m, 1, 0.08, 7),
        ]
    inst = _render(mapping_out(tmp_path, fdr=rows), cfg)
    a = _rows(inst, "a")
    assert list(dict.fromkeys(a["length_class"].astype(str))) == ["7-9", "10-40"]
    assert len(a) == 2 * 3 * 2


# --- panel (b) --------------------------------------------------------------------------------


def test_example4_stratum_rows_and_title(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_ex4(tmp_path), cfg)
    b = _rows(inst, "b")
    fdr = b[b["metric"] == "fdr_stratum"]
    assert len(fdr) == 3
    assert _num(_cell(b, "all", 0, "fdr_stratum")["value"]) == pytest.approx(0.0033)
    assert _num(_cell(b, "all", 1, "fdr_stratum")["value"]) == pytest.approx(0.025)
    assert _isna(_cell(b, "all", 2, "fdr_stratum")["value"])
    accepted = {m: _num(_cell(b, "all", m, "n_accepted")["value"]) for m in STRATA}
    assert accepted == {0: 880, 1: 40, 2: 0}
    assert inst.title == EX4_TITLE
    assert EX4_TITLE in _fig_texts(inst)
    assert "express" not in inst.title.lower()


def test_alpha_line_is_dashed_at_the_default(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_ex4(tmp_path), cfg)
    assert _hlines(_ax(inst, "b"), "--") == [pytest.approx(cfg.stats.alpha.placement)]
    assert cfg.stats.alpha.placement == 0.01


def test_alpha_line_and_title_follow_the_config(tmp_path: Path) -> None:
    cfg = resolve_config(overrides=["stats.alpha.placement=0.05"])
    inst = _render(_ex4(tmp_path), cfg)
    assert _hlines(_ax(inst, "b"), "--") == [pytest.approx(0.05)]
    assert inst.title == (
        "Placement FDR ≤ 0.05 in strata 0, 1; 920 placements accepted; "
        "stratum 2 unestimated (no targets)"
    )
    assert inst.data.equals(_render(_ex4(tmp_path / "again"), resolve_config()).data)


def test_points_accepted_labels_and_na_marker(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_ex4(tmp_path), cfg)
    ax = _ax(inst, "b")
    pts = _markers(ax)
    filled = sorted((x, y) for x, y, hollow in pts if not hollow)
    assert (0.0, pytest.approx(0.0033)) in filled
    assert (1.0, pytest.approx(0.025)) in filled
    hollow = [(x, y) for x, y, h in pts if h]
    assert [x for x, _ in hollow] == [2.0]
    assert hollow[0][1] == pytest.approx(ax.get_ylim()[0], rel=1e-6, abs=1e-12)
    texts = _notes(ax)
    assert "NA (no targets)" in texts
    for n in ("880", "40", "0"):
        assert n in texts, n


def test_every_stratum_estimated_and_below_alpha(tmp_path: Path, cfg: GiapConfig) -> None:
    totals = {0: (500, 1, 0.002, 450), 1: (80, 0, 0.0, 70), 2: (30, 0, 0.005, 20)}
    inst = _render(mapping_out(tmp_path, fdr=full_placement(totals)), cfg)
    assert inst.title == "Placement FDR ≤ 0.01 in strata 0, 1, 2; 540 placements accepted"
    assert "NA (no targets)" not in _notes(_ax(inst, "b"))


def test_every_target_count_zero(tmp_path: Path, cfg: GiapConfig) -> None:
    """The 2026-10-22 validation run's first rows: (a) at the floor, (b) NA, title unestimated."""
    inst = _render(mapping_out(tmp_path, fdr=full_placement({})), cfg)
    a = _rows(inst, "a")
    assert len(a) == 42
    assert all(_num(v) == 0 for v in a["value"])
    b = _rows(inst, "b")
    assert all(_isna(v) for v in b[b["metric"] == "fdr_stratum"]["value"])
    ax = _ax(inst, "b")
    hollow = [x for x, _, h in _markers(ax) if h]
    assert sorted(hollow) == [0.0, 1.0, 2.0]
    assert "unestimated" in inst.title
    assert not inst.title.startswith("Placement FDR ≤")  # no stratum is estimated


def test_an_all_zero_decoy_type_does_not_change_v6(tmp_path: Path, cfg: GiapConfig) -> None:
    only = _render(_ex4(tmp_path / "one"), cfg)
    both_rows = [*_ex4_rows(), *full_placement({}, decoy_type="shuffle")]
    both = _render(mapping_out(tmp_path / "two", fdr=both_rows), cfg)
    pd.testing.assert_frame_equal(only.data, both.data)
    assert only.title == both.title


def test_header_only_mapping_fdr_is_no_peptide_placements(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(mapping_out(tmp_path, fdr=[]), cfg)
    assert list(inst.data.columns) == COLUMNS
    assert len(inst.data) == 0
    for label in ("a", "b"):
        ax = _ax(inst, label)
        assert "no peptide placements" in _notes(ax), label
        assert _markers(ax) == []
    assert "%" not in inst.title


# --- panel (c) --------------------------------------------------------------------------------


def test_example6_denovo_precision_with_footnote(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_ex4(tmp_path, EX6_DENOVO), cfg)
    c = _rows(inst, "c")
    assert len(c) == 3
    assert list(c["metric"]) == ["precision"] * 3
    assert [_num(v) for v in c["value"]] == pytest.approx([0.8, 0.9, 0.95])
    assert [_num(v) for v in c["ci_low"]] == pytest.approx([0.7, 0.82, 0.9])
    assert [_num(v) for v in c["ci_high"]] == pytest.approx([0.88, 0.95, 0.98])
    assert all("c_assumed" in str(f) for f in c["flags"])
    assert all(_isna(v) for v in c["length_class"])
    assert all(_isna(v) for v in c["n_diff"])
    assert FOOTNOTE in _fig_texts(inst)
    ax = _ax(inst, "c")
    pts = _markers(ax)
    assert sorted(y for _, y, _ in pts) == pytest.approx([0.8, 0.9, 0.95])
    assert all(h for _, _, h in pts)  # c assumed: hollow


def test_denovo_without_c_assumed_has_no_footnote(tmp_path: Path, cfg: GiapConfig) -> None:
    rows = [denovo_row(1, 0.85, 0.8, 0.9), denovo_row(2, 0.93, 0.9, 0.96)]
    inst = _render(_ex4(tmp_path, rows), cfg)
    assert FOOTNOTE not in _fig_texts(inst)
    pts = _markers(_ax(inst, "c"))
    assert len(pts) == 2
    assert not any(h for _, _, h in pts)


def test_denovo_error_bars_span_the_ci(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_ex4(tmp_path, EX6_DENOVO), cfg)
    ax = _ax(inst, "c")
    spans = []
    for coll in ax.collections:
        for seg in getattr(coll, "get_segments", lambda: [])():
            ys = sorted(float(p[1]) for p in seg)
            if len(ys) == 2 and ys[0] != ys[1]:
                spans.append((round(ys[0], 6), round(ys[1], 6)))
    for line in ax.get_lines():
        ys = sorted(float(v) for v in line.get_ydata())
        xs = {float(v) for v in line.get_xdata()}
        if len(ys) == 2 and len(xs) == 1 and ys[0] != ys[1]:
            spans.append((round(ys[0], 6), round(ys[1], 6)))
    assert {(0.7, 0.88), (0.82, 0.95), (0.9, 0.98)} <= set(spans)


# --- properties -------------------------------------------------------------------------------

_FDR = st.one_of(st.none(), st.floats(min_value=0.0, max_value=0.2, allow_nan=False))


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    fdrs=st.tuples(_FDR, _FDR, _FDR),
    accepted=st.tuples(*(st.integers(min_value=0, max_value=500) for _ in STRATA)),
    alpha=st.sampled_from([0.01, 0.05, 0.1]),
)
def test_property_title_lists_exactly_the_strata_at_or_below_alpha(
    fdrs: tuple[float | None, ...], accepted: tuple[int, ...], alpha: float
) -> None:
    style.apply()
    totals = {
        m: ((100 if f is not None else 0), 1, f, (accepted[m] if f is not None else 0))
        for m, f in zip(STRATA, fdrs, strict=True)
    }
    cfg = resolve_config(overrides=[f"stats.alpha.placement={alpha}"])
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(mapping_out(Path(tmp), fdr=full_placement(totals)), cfg)
    b = _rows(inst, "b")
    for m, f in zip(STRATA, fdrs, strict=True):
        v = _cell(b, "all", m, "fdr_stratum")["value"]
        assert _isna(v) if f is None else _num(v) == pytest.approx(f)
    listed = [str(m) for m, f in zip(STRATA, fdrs, strict=True) if f is not None and f <= alpha]
    total = sum(totals[m][3] for m in STRATA)
    if listed:
        assert inst.title.startswith(
            f"Placement FDR ≤ {alpha:g} in strata {', '.join(listed)}; {total} placements accepted"
        )
    na = [m for m, f in zip(STRATA, fdrs, strict=True) if f is None]
    if na:
        assert "unestimated (no targets)" in inst.title
    else:
        assert "unestimated" not in inst.title
    assert _hlines(_ax(inst, "b"), "--") == [pytest.approx(alpha)]


# --- the golden M3 tables, files and determinism ----------------------------------------------


def test_renders_on_the_sim50k_golden_tables(tmp_path: Path, cfg: GiapConfig) -> None:
    out = mapping_out(tmp_path)
    (out / "mapping").mkdir(parents=True, exist_ok=True)
    (out / "mapping" / "mapping_fdr.tsv").write_bytes((GOLDEN / "mapping_fdr.tsv").read_bytes())
    inst = _render(out, cfg)
    b = _rows(inst, "b")
    assert _num(_cell(b, "all", 1, "fdr_stratum")["value"]) == pytest.approx(0.25)
    assert _isna(_cell(b, "all", 0, "fdr_stratum")["value"])
    assert _num(_cell(_rows(inst, "a"), "15-19", 1, "targets")["value"]) == 2
    assert len(_rows(inst, "a")) == 42
    assert "unestimated" in inst.title
    assert "in strata 1" not in inst.title  # 0.25 > 0.01: no claim
    assert "c" not in {ax.get_label() for ax in inst.figure.axes}


def test_data_tsv_and_images_are_byte_identical(tmp_path: Path, cfg: GiapConfig) -> None:
    out = _ex4(tmp_path / "run", EX6_DENOVO)
    written = []
    for name in ("x", "y"):
        inst = _render(out, cfg)
        paths = write_instance(
            inst, "V06", "mapping_fdr", tmp_path / name, ["svg", "png", "pdf"], 100
        )
        written.append({p.relative_to(tmp_path / name): p.read_bytes() for p in paths})
    assert set(written[0]) == {
        Path("V06_mapping_fdr.svg"),
        Path("V06_mapping_fdr.png"),
        Path("V06_mapping_fdr.pdf"),
        Path("data/V06_mapping_fdr.tsv"),
    }
    assert written[0] == written[1]
    tsv = written[0][Path("data/V06_mapping_fdr.tsv")].decode("utf-8").splitlines()
    assert tsv[0] == "\t".join(COLUMNS)
    assert tsv[1].startswith("a\t7\t0\tNA\tdecoys\t0\t")
    assert any(line.startswith("b\tall\t2\tNA\tfdr_stratum\tNA\t") for line in tsv)
    assert any(line.startswith("c\tNA\tNA\t1\tprecision\t0.8\t0.7\t0.88\t") for line in tsv)


def test_render_does_not_mutate_between_calls(tmp_path: Path, cfg: GiapConfig) -> None:
    view = RunView.from_outdir(_ex4(tmp_path, EX6_DENOVO))
    first = _v06().fn(view, cfg)
    second = _v06().fn(view, cfg)
    assert isinstance(first, list) and isinstance(second, list)
    pd.testing.assert_frame_equal(first[0].data, second[0].data)
    assert first[0].title == second[0].title
