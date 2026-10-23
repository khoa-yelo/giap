"""``giap.viz.checklist``: the §14.4 style checklist as code (M9-24 examples 1-6, degenerate cases,
lead decisions 1-3; M9-index decisions 30-31; GIAP_SPEC.md §10.1, §14.4).

Crafted figures each violate one item (or none). They are written with
``giap.io.plot_files.write_instance`` so ``check_instance`` sees real ``files``.

Interface assumed beyond the spec text (the spec's public interface is ``CheckResult``,
``check_instance`` and ``score``):

- ``check_instance`` returns exactly one ``CheckResult`` per item, ``item`` spelled as in the spec's
  table (``ITEMS`` below), in that order; ``plot`` starts with the plot ID.
- the title item reads ``Instance.title``; an instance is an empty state when its ``data`` has no
  rows (``base.empty_state``).
- ``checklist.MACHADO_DEUTERANOPIA``: the 3x3 matrix (rows act on linear RGB column vectors);
  ``checklist.simulate_deuteranopia(srgb)``: float sRGB in [0, 1], shape ``(..., 3)`` -> same
  shape, the matrix applied in linear RGB and the result clipped to [0, 1] and re-encoded.
- the colour-blind item writes one extra PNG into the directory of the instance's PNG (any name).
- line colours count as drawn edge colours for the palette item; alpha is ignored.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.figure import Figure

from giap.io.plot_files import write_instance
from giap.viz import base, style
from giap.viz import checklist as ck
from giap.viz.base import Instance

ITEMS = (
    "title states the finding",
    "axes carry units",
    "palette as §10.1",
    "no overlapping labels",
    "legible at single-column width",
    "data TSV written",
    "colour-blind check",
    "≤ 20,000 drawn elements",
)
TITLE, UNITS, PALETTE, OVERLAP, LEGIBLE, TSV, CVD, ELEMENTS = ITEMS
GATED = tuple(i for i in ITEMS if i != CVD)
GOOD_TITLE = "Prophage genes: 12 % detected vs 31 % expected"
PID = "V90"

# Machado, Oliveira & Fernandes 2009, IEEE TVCG 15(6), Table 1, deuteranomaly, severity 1.0
MACHADO_D10 = np.array(
    [
        [0.367322, 0.860646, -0.227968],
        [0.280085, 0.672501, 0.047413],
        [-0.011820, 0.042940, 0.968881],
    ]
)


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


# ------------------------------------------------------------------ builders


def _figure(width_in: float = style.SINGLE_IN, height_in: float = 2.5) -> Figure:
    return Figure(figsize=(width_in, height_in))


def _good_axes(fig: Figure) -> None:
    """Three detection-coloured bars, labelled axes with units, readable ticks."""
    ax = fig.add_subplot()
    colours = list(style.DETECTION.values())
    ax.bar([0, 1, 2], [3.0, 1.0, 2.0], color=colours)
    ax.set_xticks([0, 1, 2], ["detected", "not detected", "unestimated"])
    ax.set_ylabel("log2 abundance")
    ax.set_title(GOOD_TITLE)
    fig.tight_layout()


def _data(n: int = 3) -> pd.DataFrame:
    return pd.DataFrame({"category": [f"c{i}" for i in range(n)], "value": np.arange(n) * 1.0})


def _inst(fig: Figure, title: str = GOOD_TITLE, data: pd.DataFrame | None = None) -> Instance:
    return Instance(suffix=None, figure=fig, data=_data() if data is None else data, title=title)


def _files(inst: Instance, tmp: Path, formats: Sequence[str] = ("svg", "png")) -> list[Path]:
    return write_instance(inst, PID, "fake", tmp / "plots", list(formats), 100)


def _check(
    inst: Instance, tmp: Path, formats: Sequence[str] = ("svg", "png")
) -> dict[str, ck.CheckResult]:
    results = ck.check_instance(inst, PID, _files(inst, tmp, formats))
    return {r.item: r for r in results}


def _good(tmp: Path) -> dict[str, ck.CheckResult]:
    fig = _figure()
    _good_axes(fig)
    return _check(_inst(fig), tmp)


# ------------------------------------------------------------------ shape of the result


def test_one_result_per_item_in_the_spec_order(tmp_path: Path) -> None:
    fig = _figure()
    _good_axes(fig)
    inst = _inst(fig)
    results = ck.check_instance(inst, PID, _files(inst, tmp_path))
    assert [r.item for r in results] == list(ITEMS)
    assert all(isinstance(r, ck.CheckResult) for r in results)
    assert all(r.plot.startswith(PID) for r in results)
    assert all(isinstance(r.detail, str) for r in results)


def test_check_result_is_frozen() -> None:
    r = ck.CheckResult(plot="V90_fake", item=TITLE, passed=True, detail="")
    with pytest.raises(AttributeError):
        r.passed = False  # type: ignore[misc]


def test_a_clean_figure_passes_every_gated_item_and_scores_one(tmp_path: Path) -> None:
    res = _good(tmp_path)
    failed = {k: r.detail for k, r in res.items() if k in GATED and r.passed is not True}
    assert failed == {}
    assert ck.score(list(res.values())) == 1.0


# ------------------------------------------------------------------ title


@pytest.mark.parametrize(
    ("title", "ok"),
    [
        (GOOD_TITLE, True),
        ("Operon partners: \u03c1 = 0.41", True),
        ("", False),
        ("Abundance by COG category", False),
        ("Prophage genes: 12 % expressed", False),
    ],
)
def test_title_states_the_finding(tmp_path: Path, title: str, ok: bool) -> None:
    fig = _figure()
    _good_axes(fig)
    assert _check(_inst(fig, title=title), tmp_path)[TITLE].passed is ok


def test_empty_state_message_counts_as_a_finding(tmp_path: Path) -> None:
    inst = base.empty_state(PID, "fake", "No islands at FWER ≤ 0.05 in this run", ["a"], "single")
    assert _check(inst, tmp_path)[TITLE].passed is True


# ------------------------------------------------------------------ units (example 2)


@pytest.mark.parametrize(
    ("label", "ok"),
    [
        ("length", False),  # example 2
        ("length (bp)", True),
        ("position (Mb)", True),
        ("log2 abundance", True),
        ("q", True),
        ("lfdr", True),
        ("NES", True),
        ("\u03c1", True),
        ("counts", True),
        ("", True),  # no label, no quantity
    ],
)
def test_axes_carry_units(tmp_path: Path, label: str, ok: bool) -> None:
    fig = _figure()
    _good_axes(fig)
    fig.axes[0].set_xlabel(label)
    res = _check(_inst(fig), tmp_path)[UNITS]
    assert res.passed is ok, res.detail
    if not ok:
        assert label in res.detail


def test_units_on_the_y_axis_too(tmp_path: Path) -> None:
    fig = _figure()
    _good_axes(fig)
    fig.axes[0].set_ylabel("length")
    res = _check(_inst(fig), tmp_path)[UNITS]
    assert res.passed is False
    assert "length" in res.detail


def test_hidden_axis_labels_are_not_checked(tmp_path: Path) -> None:
    fig = _figure()
    _good_axes(fig)
    ax = fig.add_axes((0.7, 0.7, 0.2, 0.2))
    ax.set_xlabel("length")
    ax.set_axis_off()
    assert _check(_inst(fig), tmp_path)[UNITS].passed is True


# ------------------------------------------------------------------ palette (example 3)


def test_example_3_a_red_bar_fails_naming_the_colour(tmp_path: Path) -> None:
    fig = _figure()
    _good_axes(fig)
    fig.axes[0].bar([3], [1.0], color="#FF0000")
    res = _check(_inst(fig), tmp_path)[PALETTE]
    assert res.passed is False
    assert "#FF0000" in res.detail.upper()


def test_a_red_edge_fails(tmp_path: Path) -> None:
    fig = _figure()
    _good_axes(fig)
    fig.axes[0].bar([3], [1.0], color="white", edgecolor="#FF0000", linewidth=1)
    res = _check(_inst(fig), tmp_path)[PALETTE]
    assert res.passed is False
    assert "#FF0000" in res.detail.upper()


def test_a_red_line_fails(tmp_path: Path) -> None:
    fig = _figure()
    _good_axes(fig)
    fig.axes[0].plot([0, 2], [1, 2], color="#FF0000")
    assert _check(_inst(fig), tmp_path)[PALETTE].passed is False


@pytest.mark.parametrize(
    "colours",
    [
        list(style.OKABE_ITO.values()),
        list(style.REGION.values()),
        list(style.DETECTION.values()),
        list(style.COG_SUPER_COLOUR.values()),
        [c for c, _ in style.ISLAND.values()] + ["#000000", "#FFFFFF", "#BBBBBB"],
    ],
    ids=["okabe_ito", "region", "detection", "cog", "island_and_neutrals"],
)
def test_section_10_1_colours_pass(tmp_path: Path, colours: list[str]) -> None:
    fig = _figure()
    ax = fig.add_subplot()
    ax.bar(range(len(colours)), [1.0] * len(colours), color=colours)
    ax.set_xticks([])
    ax.set_ylabel("q")
    res = _check(_inst(fig), tmp_path)[PALETTE]
    assert res.passed is True, res.detail


def test_alpha_is_ignored(tmp_path: Path) -> None:
    fig = _figure()
    _good_axes(fig)
    fig.axes[0].bar([3], [1.0], color="#0072B2", alpha=0.3)
    assert _check(_inst(fig), tmp_path)[PALETTE].passed is True


@pytest.mark.parametrize(("cmap", "ok"), [("viridis", True), ("cividis", True), ("jet", False)])
def test_continuous_mappables_use_viridis_or_cividis(tmp_path: Path, cmap: str, ok: bool) -> None:
    fig = _figure()
    ax = fig.add_subplot()
    ax.imshow(np.arange(12.0).reshape(3, 4), cmap=cmap)
    ax.set_xticks([])
    ax.set_yticks([])
    res = _check(_inst(fig), tmp_path)[PALETTE]
    assert res.passed is ok, res.detail
    if not ok:
        assert cmap in res.detail


@pytest.mark.parametrize(("cmap", "ok"), [("viridis", True), ("plasma", False)])
def test_scatter_colourmaps(tmp_path: Path, cmap: str, ok: bool) -> None:
    fig = _figure()
    ax = fig.add_subplot()
    ax.scatter([0, 1, 2], [0, 1, 2], c=[0.0, 0.5, 1.0], cmap=cmap)
    ax.set_xlabel("position (Mb)")
    ax.set_ylabel("log2 abundance")
    assert _check(_inst(fig), tmp_path)[PALETTE].passed is ok


def test_cog_and_detection_palettes_never_in_one_figure(tmp_path: Path) -> None:
    fig = _figure(style.DOUBLE_IN)
    left, right = fig.subplots(1, 2)
    cog = list(style.COG_SUPER_COLOUR.values())
    det = list(style.DETECTION.values())
    left.bar(range(len(cog)), [1.0] * len(cog), color=cog)
    right.bar(range(len(det)), [1.0] * len(det), color=det)
    for ax in (left, right):
        ax.set_xticks([])
    res = _check(_inst(fig), tmp_path)[PALETTE]
    assert res.passed is False
    assert "COG" in res.detail or "cog" in res.detail


# ------------------------------------------------------------------ overlaps (example 1)


def _standalone(token: str, text: str) -> bool:
    return re.search(rf"(?<![\w.]){re.escape(token)}(?![\w.])", text) is not None


def test_example_1_overlapping_tick_labels_fail_naming_both(tmp_path: Path) -> None:
    fig = _figure()
    ax = fig.add_subplot()
    ax.plot([0, 2], [0, 1], color="#0072B2")
    ax.set_xlim(0, 2)
    ax.set_xticks([1.0, 1.001], ["1", "10"])
    ax.set_xlabel("position (Mb)")
    ax.set_ylabel("q")
    res = _check(_inst(fig), tmp_path)[OVERLAP]
    assert res.passed is False
    assert _standalone("1", res.detail) and _standalone("10", res.detail), res.detail


def test_overlapping_annotations_fail(tmp_path: Path) -> None:
    fig = _figure()
    _good_axes(fig)
    ax = fig.axes[0]
    ax.text(1.0, 2.5, "alpha")
    ax.text(1.0, 2.5, "beta")
    res = _check(_inst(fig), tmp_path)[OVERLAP]
    assert res.passed is False
    assert "alpha" in res.detail and "beta" in res.detail


def test_invisible_and_empty_texts_are_ignored(tmp_path: Path) -> None:
    fig = _figure()
    _good_axes(fig)
    ax = fig.axes[0]
    ax.text(1.0, 2.5, "alpha")
    hidden = ax.text(1.0, 2.5, "beta")
    hidden.set_visible(False)
    ax.text(1.0, 2.5, "")
    assert _check(_inst(fig), tmp_path)[OVERLAP].passed is True


@pytest.mark.parametrize(("rotation", "ok"), [(90, True), (0, False)])
def test_rotated_labels_use_their_rotated_box(tmp_path: Path, rotation: int, ok: bool) -> None:
    """Twelve 8-character labels: side by side they collide, rotated 90° they do not."""
    fig = _figure(style.SINGLE_IN, 3.0)
    ax = fig.add_subplot()
    ax.bar(range(12), [1.0] * 12, color="#0072B2")
    ax.set_xticks(range(12), [f"label_{i:02d}" for i in range(12)], rotation=rotation)
    ax.set_ylabel("q")
    fig.tight_layout()
    res = _check(_inst(fig), tmp_path)[OVERLAP]
    assert res.passed is ok, res.detail


# ------------------------------------------------------------------ legibility (example 4)


def test_example_4_a_5pt_annotation_fails(tmp_path: Path) -> None:
    fig = _figure()
    _good_axes(fig)
    fig.axes[0].annotate("tiny note", (0.2, 0.9), xycoords="axes fraction", fontsize=5)
    res = _check(_inst(fig), tmp_path)[LEGIBLE]
    assert res.passed is False
    assert "tiny note" in res.detail or "5" in res.detail


@pytest.mark.parametrize(("width", "ok"), [(3.5, True), (7.2, True), (5.0, False), (3.0, False)])
def test_figure_width_is_single_or_double_column(tmp_path: Path, width: float, ok: bool) -> None:
    fig = _figure(width)
    _good_axes(fig)
    assert _check(_inst(fig), tmp_path)[LEGIBLE].passed is ok


@settings(
    max_examples=12,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
@given(size=st.sampled_from([3.0, 4.0, 5.0, 5.5, 5.9, 6.0, 6.5, 7.0, 8.0, 10.0]))
def test_legibility_threshold_is_6pt(tmp_path_factory: pytest.TempPathFactory, size: float) -> None:
    """Index decision 31 / lead decision 2: every text >= 6 pt passes, anything smaller fails."""
    fig = _figure()
    _good_axes(fig)
    fig.axes[0].text(0.1, 0.9, "x", transform=fig.axes[0].transAxes, fontsize=size)
    res = _check(_inst(fig), tmp_path_factory.mktemp("leg"))[LEGIBLE]
    assert res.passed is (size >= 6.0), (size, res.detail)


# ------------------------------------------------------------------ data TSV


def test_missing_tsv_fails(tmp_path: Path) -> None:
    fig = _figure()
    _good_axes(fig)
    inst = _inst(fig)
    files = [p for p in _files(inst, tmp_path) if p.suffix != ".tsv"]
    res = {r.item: r for r in ck.check_instance(inst, PID, files)}
    assert res[TSV].passed is False


def test_tsv_listed_but_deleted_fails(tmp_path: Path) -> None:
    fig = _figure()
    _good_axes(fig)
    inst = _inst(fig)
    files = _files(inst, tmp_path)
    next(p for p in files if p.suffix == ".tsv").unlink()
    res = {r.item: r for r in ck.check_instance(inst, PID, files)}
    assert res[TSV].passed is False


def test_tsv_must_sit_under_data_and_carry_the_plot_id(tmp_path: Path) -> None:
    fig = _figure()
    _good_axes(fig)
    inst = _inst(fig)
    files = _files(inst, tmp_path)
    tsv = next(p for p in files if p.suffix == ".tsv")
    stray = tmp_path / "plots" / "V91_other.tsv"
    stray.write_bytes(tsv.read_bytes())
    tsv.unlink()
    kept = [p for p in files if p.suffix != ".tsv"] + [stray]
    res = {r.item: r for r in ck.check_instance(inst, PID, kept)}
    assert res[TSV].passed is False


def test_header_only_tsv_is_accepted_for_an_empty_state(tmp_path: Path) -> None:
    inst = base.empty_state(PID, "fake", "No gene pairs in this run", ["a", "b"], "single")
    assert _check(inst, tmp_path)[TSV].passed is True


def test_tsv_rows_are_read_from_the_file(tmp_path: Path) -> None:
    """The data frame has rows but the written TSV was truncated to its header: fail."""
    fig = _figure()
    _good_axes(fig)
    inst = _inst(fig)
    files = _files(inst, tmp_path)
    tsv = next(p for p in files if p.suffix == ".tsv")
    tsv.write_text(tsv.read_text(encoding="utf-8").splitlines()[0] + "\n", encoding="utf-8")
    res = {r.item: r for r in ck.check_instance(inst, PID, files)}
    assert res[TSV].passed is False


# ------------------------------------------------------------------ elements (example 6)


def test_example_6_25000_points_on_one_axes_fail(tmp_path: Path) -> None:
    rng = np.random.default_rng(0)
    fig = _figure()
    ax = fig.add_subplot()
    ax.scatter(rng.random(25_000), rng.random(25_000), s=1, color="#0072B2")
    ax.set_xlabel("position (Mb)")
    ax.set_ylabel("log2 abundance")
    res = _check(_inst(fig), tmp_path)[ELEMENTS]
    assert res.passed is False


def test_elements_are_counted_per_axes(tmp_path: Path) -> None:
    rng = np.random.default_rng(1)
    fig = _figure(style.DOUBLE_IN)
    for ax in fig.subplots(1, 2):
        ax.scatter(rng.random(15_000), rng.random(15_000), s=1, color="#0072B2")
        ax.set_xlabel("position (Mb)")
        ax.set_ylabel("log2 abundance")
    assert _check(_inst(fig), tmp_path)[ELEMENTS].passed is True


def test_a_long_line_is_one_element(tmp_path: Path) -> None:
    fig = _figure()
    ax = fig.add_subplot()
    x = np.linspace(0, 1, 25_000)
    ax.plot(x, np.sin(x), color="#0072B2")
    ax.set_xlabel("position (Mb)")
    ax.set_ylabel("log2 abundance")
    assert _check(_inst(fig), tmp_path)[ELEMENTS].passed is True


def test_polar_track_axes_are_counted(tmp_path: Path) -> None:
    rng = np.random.default_rng(2)
    fig = _figure(style.SINGLE_IN, style.SINGLE_IN)
    ax = fig.add_subplot(projection="polar")
    ax.scatter(rng.random(25_000) * 6.28, rng.random(25_000), s=1, color="#0072B2")
    assert _check(_inst(fig), tmp_path)[ELEMENTS].passed is False


def test_many_patches_fail(tmp_path: Path) -> None:
    fig = _figure()
    ax = fig.add_subplot()
    ax.bar(np.arange(20_001), np.ones(20_001), color="#0072B2")
    ax.set_xticks([])
    ax.set_ylabel("q")
    assert _check(_inst(fig), tmp_path)[ELEMENTS].passed is False


# ------------------------------------------------------------------ empty state


def test_empty_state_marks_units_palette_elements_not_applicable(tmp_path: Path) -> None:
    inst = base.empty_state(PID, "fake", "No gene pairs in this run", ["a", "b"], "double")
    res = _check(inst, tmp_path)
    assert res[UNITS].passed is None
    assert res[PALETTE].passed is None
    assert res[ELEMENTS].passed is None
    for item in (TITLE, OVERLAP, LEGIBLE, TSV):
        assert res[item].passed is True, (item, res[item].detail)
    assert ck.score(list(res.values())) == 1.0


# ------------------------------------------------------------------ colour-blind (example 5)


def _numbers(text: str) -> list[float]:
    return [float(x) for x in re.findall(r"\d+(?:\.\d+)?", text)]


def test_example_5_detection_colours_are_recorded_not_gated(tmp_path: Path) -> None:
    res = _good(tmp_path)[CVD]
    assert res.passed is None  # lead decision 1, index decision 30
    assert "ΔE" in res.detail
    assert len(_numbers(res.detail)) >= 2  # before and after


def test_simulated_png_is_written_next_to_the_original(tmp_path: Path) -> None:
    fig = _figure()
    _good_axes(fig)
    inst = _inst(fig)
    files = _files(inst, tmp_path)
    png = next(p for p in files if p.suffix == ".png")
    before = set(png.parent.glob("*.png"))
    ck.check_instance(inst, PID, files)
    new = set(png.parent.glob("*.png")) - before
    assert len(new) == 1
    assert next(iter(new)).stat().st_size > 0


def test_svg_only_still_records_delta_e(tmp_path: Path) -> None:
    fig = _figure()
    _good_axes(fig)
    res = _check(_inst(fig), tmp_path, formats=("svg",))[CVD]
    assert res.passed is None
    assert len(_numbers(res.detail)) >= 2


def test_machado_matrix_is_the_published_deuteranomaly_severity_1() -> None:
    np.testing.assert_allclose(np.asarray(ck.MACHADO_DEUTERANOPIA), MACHADO_D10, atol=5e-7)


def test_module_docstring_cites_machado_2009() -> None:
    doc = ck.__doc__ or ""
    assert "Machado" in doc and "2009" in doc


def _to_linear(c: np.ndarray) -> np.ndarray:
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _to_srgb(c: np.ndarray) -> np.ndarray:
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * c ** (1 / 2.4) - 0.055)


def _oracle(srgb: np.ndarray) -> np.ndarray:
    lin = _to_linear(np.asarray(srgb, dtype=float))
    out = np.clip(lin @ MACHADO_D10.T, 0.0, 1.0)
    return _to_srgb(out)


def test_simulation_is_applied_in_linear_rgb() -> None:
    red = np.array([1.0, 0.0, 0.0])
    got = np.asarray(ck.simulate_deuteranopia(red))
    np.testing.assert_allclose(got, _oracle(red), atol=2e-3)
    # gamma-space application would give the raw matrix column (0.367, 0.280, 0)
    assert abs(got[0] - 0.367322) > 0.1


def test_simulation_of_the_detection_palette_matches_the_oracle() -> None:
    hexes = list(style.DETECTION.values())
    rgb = np.array([[int(h[i : i + 2], 16) / 255 for i in (1, 3, 5)] for h in hexes])
    np.testing.assert_allclose(np.asarray(ck.simulate_deuteranopia(rgb)), _oracle(rgb), atol=2e-3)


unit = st.floats(min_value=0.0, max_value=1.0, allow_nan=False)


@given(g=unit)
def test_greys_are_fixed_points(g: float) -> None:
    """Every Machado row sums to 1, so achromatic colours are unchanged."""
    got = np.asarray(ck.simulate_deuteranopia(np.array([g, g, g])))
    np.testing.assert_allclose(got, [g, g, g], atol=2e-3)


@given(
    arr=st.lists(st.tuples(unit, unit, unit), min_size=1, max_size=12).map(
        lambda xs: np.array(xs, dtype=float)
    )
)
def test_simulation_stays_in_gamut_and_keeps_shape(arr: np.ndarray) -> None:
    img = arr.reshape(1, -1, 3)
    got = np.asarray(ck.simulate_deuteranopia(img))
    assert got.shape == img.shape
    assert np.all(got >= 0.0) and np.all(got <= 1.0)
    np.testing.assert_allclose(got, _oracle(img), atol=2e-3)


# ------------------------------------------------------------------ score


def test_score_counts_passed_over_gated() -> None:
    rs = [
        ck.CheckResult("V90_fake", TITLE, True, ""),
        ck.CheckResult("V90_fake", UNITS, False, "length"),
        ck.CheckResult("V90_fake", CVD, None, "ΔE 10 → 5"),
        ck.CheckResult("V90_fake", TSV, True, ""),
    ]
    assert ck.score(rs) == pytest.approx(2 / 3)


@given(st.lists(st.sampled_from([True, False, None]), min_size=1, max_size=40))
def test_score_property(flags: list[bool | None]) -> None:
    rs = [ck.CheckResult("V90_fake", f"item{i}", f, "") for i, f in enumerate(flags)]
    gated = [f for f in flags if f is not None]
    if not gated:
        return
    s = ck.score(rs)
    assert 0.0 <= s <= 1.0
    assert s == pytest.approx(sum(f is True for f in gated) / len(gated))
