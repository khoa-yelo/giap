"""V10 ``correlogram``: is abundance (and detection) autocorrelated along the replicon? (M9-09;
GIPPA_SPEC.md §10.2 V10, §7.3.5(c); ADR-0020: raw ``value``, envelopes ``null_lo``/``null_hi`` "for
the V10 figure"; lead decision 2026-10-22: filled markers follow the TU-block null, both envelopes
are drawn).

V10 draws ``tables/autocorrelation.tsv`` and recomputes nothing: the tests patch every public
function of ``gippa.stats.autocorr`` to raise. The oracle for markers, bands, lines and the claim is
the fixture table itself (its ``null_lo``, ``null_hi``, ``null_mean`` and ``q_value`` columns).

Interface assumed where the spec leaves it open:

- The plot is registered as ``V10`` / ``correlogram`` by some ``gippa.viz`` module (found by
  walking the package); ``autocorrelation`` is in ``needs``. One double-width ``Instance``,
  ``suffix=None``.
- Panel axes carry ``get_label()`` ``"a:<contig>"`` (Moran's I) and ``"b:<contig>"`` (join
  count); columns run left to right in FASTA (preflight) order. An all-NA contig keeps its column
  (at least one axes labelled ``a:`` or ``b:<contig>``) with the "not estimable (<flags>)" note.
- The data TSV keeps every ``morans_i`` and ``join_count_bb`` row (NA rows included, so their
  ``flags`` reach the TSV); ``significant`` is the row's own ``q_value ≤ stats.alpha.context``
  (NA q -> false); the drawn marker uses the TU-block row's.
- Title numbers print with ``%g``; the title's q is the TU-block row's.
"""

from __future__ import annotations

import math
import os
import sys
import tempfile
from collections.abc import Callable
from pathlib import Path

import matplotlib.colors as mcolors
import pandas as pd
import pytest
from helpers.viz_landscape import (
    Cell,
    correlogram_rows,
    default_cell,
    gene,
    landscape_out,
    markers,
    na_correlogram_rows,
    plot_spec,
    text_overlaps,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.collections import LineCollection, PathCollection, PolyCollection
from matplotlib.text import Text

import gippa.stats.autocorr as autocorr
from gippa.config import GippaConfig, resolve_config
from gippa.io.plot_files import write_instance
from gippa.io.run_view import RunView
from gippa.viz import style
from gippa.viz.base import Instance, Skip

REPO = Path(__file__).resolve().parents[3]
COLUMNS = [
    "contig",
    "statistic",
    "null_model",
    "lag",
    "value",
    "null_mean",
    "null_lo",
    "null_hi",
    "effect",
    "q_value",
    "significant",
    "flags",
]
TU_BLUE = "#56B4E9"
EX4: dict[tuple[str, str, int], Cell] = {
    ("morans_i", "gene_label", 1): (0.31, -0.002, -0.04, 0.06, 0.001),
    ("morans_i", "tu_block", 1): (0.31, -0.004, -0.05, 0.12, 0.002),
    ("morans_i", "gene_label", 2): (0.08, -0.002, -0.05, 0.1, 0.4),
    ("morans_i", "tu_block", 2): (0.08, -0.004, -0.09, 0.11, 0.4),
}
EX4_TITLE = "Lag-1 Moran's I 0.31: above the TU-block envelope (q = 0.002)"
NOT_ESTIMABLE = "Autocorrelation not estimable on this genome"
CHR_A = ("chrA", 100_000, "circular", "chromosome")
CHR_B = ("chrB", 60_000, "circular", "chromosome")


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


@pytest.fixture
def no_autocorr(monkeypatch: pytest.MonkeyPatch) -> None:
    """Every public ``gippa.stats.autocorr`` function raises, also where a ``gippa.viz`` module
    imported it by name: V10 must not recompute a statistic or a null."""
    plot_spec("V10")  # import the viz modules first, so their bound names can be patched

    def boom(name: str) -> Callable[..., object]:
        def raiser(*_a: object, **_k: object) -> object:
            raise AssertionError(f"V10 called gippa.stats.autocorr.{name}")

        return raiser

    names = [n for n in getattr(autocorr, "__all__", dir(autocorr)) if not n.startswith("_")]
    for name in names:
        obj = getattr(autocorr, name, None)
        if not callable(obj) or isinstance(obj, type):
            continue
        monkeypatch.setattr(autocorr, name, boom(name))
        for mod_name, mod in list(sys.modules.items()):
            if mod_name.startswith("gippa.viz") and getattr(mod, name, None) is obj:
                monkeypatch.setattr(mod, name, boom(name))


def _render(out: Path, cfg: GippaConfig) -> Instance:
    result = plot_spec("V10").fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert isinstance(inst, Instance)
    return inst


def _out(tmp: Path, rows: list[dict], replicons: tuple = (CHR_A,)) -> Path:
    return landscape_out(tmp, replicons=replicons, genes=[gene("g1", 10_000, 1.0)], autocorr=rows)


def _ex4(tmp: Path) -> Path:
    return _out(tmp, correlogram_rows("chrA", overrides=EX4))


def _ax(inst: Instance, label: str) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_label() == label]
    assert len(found) == 1, f"axes {label!r}: {[a.get_label() for a in inst.figure.axes]}"
    return found[0]


def _contig_axes(inst: Instance, contig: str) -> list[Axes]:
    return [ax for ax in inst.figure.axes if ax.get_label() in (f"a:{contig}", f"b:{contig}")]


def _row(inst: Instance, contig: str, statistic: str, null_model: str, lag: int) -> pd.Series:
    d = inst.data
    hit = d[
        (d["contig"].astype(str) == contig)
        & (d["statistic"].astype(str) == statistic)
        & (d["null_model"].astype(str) == null_model)
        & (pd.to_numeric(d["lag"]) == lag)
    ]
    assert len(hit) == 1, f"{contig}/{statistic}/{null_model}/{lag}:\n{d}"
    return hit.iloc[0]


def _isna(value: object) -> bool:
    return value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value))


def _marker_at(ax: Axes, x: float, y: float) -> list[bool]:
    """``hollow`` of every marker drawn at (x, y)."""
    return [h for mx, my, h, _ in markers(ax) if abs(mx - x) < 1e-9 and abs(my - y) < 1e-9]


def _bands(ax: Axes) -> list[tuple[str, float, set[tuple[float, float]]]]:
    """(face colour hex, alpha, rounded vertices) of every filled polygon collection."""
    out = []
    for coll in ax.collections:
        if not isinstance(coll, PolyCollection) or isinstance(
            coll, PathCollection | LineCollection
        ):
            continue
        faces = coll.get_facecolors()
        if not len(faces):
            continue
        rgba = faces[0]
        alpha = coll.get_alpha() if coll.get_alpha() is not None else float(rgba[3])
        verts = {
            (round(float(x), 6), round(float(y), 6))
            for p in coll.get_paths()
            for x, y in p.vertices
        }
        out.append((mcolors.to_hex(rgba).upper(), alpha, verts))
    return out


def _is_light_grey(hex_colour: str) -> bool:
    r, g, b = mcolors.to_rgb(hex_colour)
    return max(r, g, b) - min(r, g, b) < 0.02 and r >= 0.6


def _dashed(ax: Axes) -> list[list[tuple[float, float]]]:
    out = []
    for line in ax.get_lines():
        if line.get_linestyle() == "--" and len(line.get_xdata()) >= 2:
            out.append(
                [
                    (round(float(x), 6), round(float(y), 6))
                    for x, y in zip(line.get_xdata(), line.get_ydata(), strict=True)
                ]
            )
    return out


def _notes(ax: Axes) -> list[str]:
    return [t.get_text() for t in ax.texts if t.get_text()]


# --- registration and skips -------------------------------------------------------------------


def test_v10_is_registered_as_correlogram() -> None:
    spec = plot_spec("V10")
    assert spec.pid == "V10"
    assert spec.name == "correlogram"
    assert "autocorrelation" in spec.needs
    assert spec.module == "core"


def test_no_autocorrelation_table_is_skipped(tmp_path: Path, cfg: GippaConfig) -> None:
    out = landscape_out(tmp_path, genes=[gene("g1", 10_000, 1.0)], autocorr=None)
    assert plot_spec("V10").fn(RunView.from_outdir(out), cfg) == Skip(
        "input_absent:autocorrelation"
    )


# --- structure --------------------------------------------------------------------------------


def test_double_width_two_panels_and_exact_columns(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_ex4(tmp_path), cfg)
    assert inst.suffix is None
    assert round(float(inst.figure.get_size_inches()[0]), 6) == style.DOUBLE_IN
    assert {"a:chrA", "b:chrA"} <= {ax.get_label() for ax in inst.figure.axes}
    assert list(inst.data.columns) == COLUMNS
    assert set(inst.data["statistic"].astype(str)) == {"morans_i", "join_count_bb"}  # no oriter
    assert len(inst.data) == 2 * 2 * 20


def test_rows_sort_by_contig_statistic_null_lag(tmp_path: Path, cfg: GippaConfig) -> None:
    rows = [*correlogram_rows("chrA", overrides=EX4), *correlogram_rows("chrB")]
    rows = rows[::-1]  # the table order must not leak through
    inst = _render(_out(tmp_path, rows, (CHR_A, CHR_B)), cfg)
    d = inst.data
    contigs = list(d["contig"].astype(str))
    assert contigs == sorted(contigs)
    keys = list(zip(contigs, d["statistic"].astype(str), d["null_model"].astype(str), strict=True))
    runs = [k for i, k in enumerate(keys) if i == 0 or keys[i - 1] != k]
    blocks = list(dict.fromkeys(keys))
    assert runs == blocks  # each (contig, statistic, null) is one contiguous block
    assert len(blocks) == 8
    for contig, statistic in {(c, s) for c, s, _ in blocks}:
        nulls = [n for c, s, n in blocks if (c, s) == (contig, statistic)]
        assert nulls == ["gene_label", "tu_block"]
    for block in blocks:
        mask = [k == block for k in keys]
        assert [int(v) for v in d.loc[mask, "lag"]] == list(range(1, 21)), block


def test_axes_labels_and_lag_ticks(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_ex4(tmp_path), cfg)
    for label in ("a:chrA", "b:chrA"):
        ax = _ax(inst, label)
        assert ax.get_xlabel() == "lag (genes)"
        ticks = [float(t) for t in ax.get_xticks() if 1 <= t <= 20]
        assert 1.0 in ticks and 20.0 in ticks
        assert all(t == int(t) for t in ticks)
    assert _ax(inst, "b:chrA").get_ylabel() == "BB join statistic (Pearson residual)"


def test_lag_ticks_follow_moran_lags(tmp_path: Path) -> None:
    cfg = resolve_config(overrides=["context.moran_lags=[1, 10]"])
    inst = _render(_out(tmp_path, correlogram_rows("chrA", lags=range(1, 11))), cfg)
    ax = _ax(inst, "a:chrA")
    lo, hi = ax.get_xlim()
    assert lo < 1 <= 10 < hi < 12
    assert len(inst.data) == 2 * 2 * 10


def test_columns_follow_fasta_order(tmp_path: Path, cfg: GippaConfig) -> None:
    rows = [*correlogram_rows("chrA"), *correlogram_rows("chrB")]
    inst = _render(_out(tmp_path, rows, (CHR_B, CHR_A)), cfg)
    assert _ax(inst, "a:chrB").get_position().x0 < _ax(inst, "a:chrA").get_position().x0
    assert _ax(inst, "a:chrB").get_position().y0 > _ax(inst, "b:chrB").get_position().y0


# --- examples 4 and 5 -------------------------------------------------------------------------


def test_example4_filled_marker_and_title(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_ex4(tmp_path), cfg)
    assert bool(_row(inst, "chrA", "morans_i", "tu_block", 1)["significant"]) is True
    assert _marker_at(_ax(inst, "a:chrA"), 1, 0.31) == [False]  # one filled marker
    assert inst.title == EX4_TITLE
    assert EX4_TITLE in [t.get_text() for t in inst.figure.findobj(Text)]
    assert "express" not in inst.title.lower()


def test_example5_inside_both_envelopes_is_open(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_ex4(tmp_path), cfg)
    for null in ("gene_label", "tu_block"):
        row = _row(inst, "chrA", "morans_i", null, 2)
        assert bool(row["significant"]) is False
        assert float(row["null_lo"]) < 0.08 < float(row["null_hi"])
    assert _marker_at(_ax(inst, "a:chrA"), 2, 0.08) == [True]


def test_one_marker_per_finite_lag(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_ex4(tmp_path), cfg)
    for label, statistic in (("a:chrA", "morans_i"), ("b:chrA", "join_count_bb")):
        xs = sorted(x for x, _, _, _ in markers(_ax(inst, label)))
        assert xs == [float(lag) for lag in range(1, 21)], statistic


@pytest.mark.parametrize(
    ("gene_q", "tu_q", "hollow"),
    [(0.001, 0.3, True), (0.2, 0.01, False), (0.05, 0.05, False), (0.001, None, True)],
)
def test_filled_markers_follow_the_tu_block_null(
    tmp_path: Path, cfg: GippaConfig, gene_q: float, tu_q: float | None, hollow: bool
) -> None:
    """Lead decision 2: the TU-block q decides (≤ alpha, inclusive); the gene-label q does not."""
    over = {
        ("morans_i", "gene_label", 3): (0.21, -0.002, -0.05, 0.05, gene_q),
        ("morans_i", "tu_block", 3): (0.21, -0.004, -0.094, 0.086, tu_q),
        ("join_count_bb", "gene_label", 3): (55.0, 38.0, 35.0, 41.0, gene_q),
        ("join_count_bb", "tu_block", 3): (55.0, 37.5, 32.5, 42.5, tu_q),
    }
    inst = _render(_out(tmp_path, correlogram_rows("chrA", overrides=over)), cfg)
    assert _marker_at(_ax(inst, "a:chrA"), 3, 0.21) == [hollow]
    assert _marker_at(_ax(inst, "b:chrA"), 3, 55.0) == [hollow]
    alpha = cfg.stats.alpha.context
    for null, q in (("gene_label", gene_q), ("tu_block", tu_q)):
        row = _row(inst, "chrA", "morans_i", null, 3)
        assert bool(row["significant"]) is (q is not None and q <= alpha)


def test_alpha_follows_the_config(tmp_path: Path) -> None:
    over = {("morans_i", "tu_block", 1): (0.31, -0.004, -0.05, 0.12, 0.03)}
    out = _out(tmp_path, correlogram_rows("chrA", overrides=over))
    strict = _render(out, resolve_config(overrides=["stats.alpha.context=0.01"]))
    loose = _render(out, resolve_config())
    assert _marker_at(_ax(strict, "a:chrA"), 1, 0.31) == [True]
    assert _marker_at(_ax(loose, "a:chrA"), 1, 0.31) == [False]


# --- envelopes and null means come from the table ---------------------------------------------


def test_both_envelopes_are_drawn_from_the_table(
    tmp_path: Path, cfg: GippaConfig, no_autocorr: None
) -> None:
    rows = correlogram_rows("chrA", overrides=EX4)
    inst = _render(_out(tmp_path, rows), cfg)
    for label, statistic in (("a:chrA", "morans_i"), ("b:chrA", "join_count_bb")):
        bands = _bands(_ax(inst, label))
        for null in ("gene_label", "tu_block"):
            want = set()
            for r in rows:
                if r["statistic"] == statistic and r["null_model"] == null:
                    want |= {
                        (float(r["lag"]), round(float(r["null_lo"]), 6)),  # type: ignore[arg-type]
                        (float(r["lag"]), round(float(r["null_hi"]), 6)),  # type: ignore[arg-type]
                    }
            match = [b for b in bands if want <= b[2]]
            assert match, f"{label}: no band holds the {null} envelope"
            colour, alpha = match[0][0], match[0][1]
            if null == "tu_block":
                assert colour == TU_BLUE
                assert alpha == pytest.approx(0.3)
            else:
                assert _is_light_grey(colour), colour


def test_null_means_are_dashed_lines(tmp_path: Path, cfg: GippaConfig, no_autocorr: None) -> None:
    rows = correlogram_rows("chrA", overrides=EX4)
    inst = _render(_out(tmp_path, rows), cfg)
    for label, statistic in (("a:chrA", "morans_i"), ("b:chrA", "join_count_bb")):
        dashed = [set(points) for points in _dashed(_ax(inst, label))]
        for null in ("gene_label", "tu_block"):
            want = {
                (float(r["lag"]), round(float(r["null_mean"]), 6))  # type: ignore[arg-type]
                for r in rows
                if r["statistic"] == statistic and r["null_model"] == null
            }
            assert any(want <= d for d in dashed), f"{label}/{null}: {dashed}"


def test_data_rows_copy_the_table(tmp_path: Path, cfg: GippaConfig, no_autocorr: None) -> None:
    inst = _render(_ex4(tmp_path), cfg)
    row = _row(inst, "chrA", "morans_i", "tu_block", 1)
    assert [float(row[c]) for c in ("value", "null_mean", "null_lo", "null_hi", "q_value")] == (
        pytest.approx([0.31, -0.004, -0.05, 0.12, 0.002])
    )
    jc = _row(inst, "chrA", "join_count_bb", "gene_label", 4)
    assert [float(jc[c]) for c in ("value", "null_mean", "null_lo", "null_hi")] == pytest.approx(
        [36.0, 38.0, 35.0, 41.0]
    )
    assert str(jc["flags"]) == "ci_unestimated"


# --- NA rows and the empty states -------------------------------------------------------------


def test_na_lags_draw_no_marker(tmp_path: Path, cfg: GippaConfig) -> None:
    na: Cell = (None, None, None, None, None)
    over = {("morans_i", n, lag): na for n in ("gene_label", "tu_block") for lag in (18, 19, 20)}
    rows = correlogram_rows("chrA", overrides=over)
    for r in rows:
        if r["statistic"] == "morans_i" and r["lag"] in (18, 19, 20):
            r["flags"] = "ci_unestimated;few_values;not_estimable"
    inst = _render(_out(tmp_path, rows), cfg)
    xs = sorted(x for x, _, _, _ in markers(_ax(inst, "a:chrA")))
    assert xs == [float(lag) for lag in range(1, 18)]
    row = _row(inst, "chrA", "morans_i", "tu_block", 19)
    assert _isna(row["value"])
    assert bool(row["significant"]) is False
    assert "few_values" in str(row["flags"])


def test_all_na_contig_shows_not_estimable(tmp_path: Path, cfg: GippaConfig) -> None:
    rows = [*correlogram_rows("chrA", overrides=EX4), *na_correlogram_rows("chrB")]
    inst = _render(_out(tmp_path, rows, (CHR_A, CHR_B)), cfg)
    axes = _contig_axes(inst, "chrB")
    assert axes, [a.get_label() for a in inst.figure.axes]
    notes = [t for ax in axes for t in _notes(ax)]
    assert any(t.startswith("not estimable (") and "few_values" in t for t in notes), notes
    assert all(markers(ax) == [] for ax in axes)
    assert inst.title == EX4_TITLE


def test_title_uses_the_first_fasta_contig_with_a_finite_row(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    rows = [*correlogram_rows("chrA", overrides=EX4), *correlogram_rows("chrB")]
    inst = _render(_out(tmp_path / "ba", rows, (CHR_B, CHR_A)), cfg)
    # chrB (first in FASTA order) lag 1: 0.16 above its TU envelope (-0.094, 0.086), q 0.5
    assert inst.title == "Lag-1 Moran's I 0.16: above the TU-block envelope (q = 0.5)"
    rows = [*correlogram_rows("chrA", overrides=EX4), *na_correlogram_rows("chrB")]
    inst = _render(_out(tmp_path / "na", rows, (CHR_B, CHR_A)), cfg)
    assert inst.title == EX4_TITLE


@pytest.mark.parametrize(("value", "word"), [(-0.2, "below"), (0.05, "within"), (0.31, "above")])
def test_title_places_lag1_against_the_tu_envelope(
    tmp_path: Path, cfg: GippaConfig, value: float, word: str
) -> None:
    over = {
        ("morans_i", "gene_label", 1): (value, -0.002, -0.01, 0.01, 0.001),
        ("morans_i", "tu_block", 1): (value, -0.004, -0.05, 0.12, 0.04),
    }
    inst = _render(_out(tmp_path, correlogram_rows("chrA", overrides=over)), cfg)
    assert inst.title == f"Lag-1 Moran's I {value:g}: {word} the TU-block envelope (q = 0.04)"


def test_validation_like_run_is_not_estimable(tmp_path: Path, cfg: GippaConfig) -> None:
    """Example 6 in miniature: every row NA (``few_values`` / ``not_estimable``) on two contigs."""
    rows = [*na_correlogram_rows("chrA"), *na_correlogram_rows("chrB")]
    inst = _render(_out(tmp_path, rows, (CHR_A, CHR_B)), cfg)
    assert inst.title == NOT_ESTIMABLE
    for contig in ("chrA", "chrB"):
        axes = _contig_axes(inst, contig)
        assert axes, contig
        notes = [t for ax in axes for t in _notes(ax)]
        assert any(t.startswith("not estimable (") for t in notes), (contig, notes)
        assert all(markers(ax) == [] for ax in axes)
    assert len(inst.data) == 2 * 2 * 2 * 20
    assert not inst.data["significant"].astype(bool).any()


def test_header_only_table_is_an_empty_state(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path, []), cfg)
    assert inst.title == NOT_ESTIMABLE
    assert list(inst.data.columns) == COLUMNS
    assert len(inst.data) == 0


# --- property ---------------------------------------------------------------------------------

_HUNDREDTHS = (
    st.integers(min_value=-40, max_value=40).filter(lambda k: k % 10).map(lambda k: k / 100)
)


@settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    value=_HUNDREDTHS,
    lo=st.integers(min_value=-30, max_value=0).map(lambda k: k / 100 - 0.005),
    width=st.integers(min_value=1, max_value=40).map(lambda k: k / 100),
    tu_q=st.one_of(st.none(), st.sampled_from([0.001, 0.01, 0.05, 0.051, 0.2, 0.9])),
    gene_q=st.sampled_from([0.001, 0.2]),
)
def test_property_marker_title_and_significance_follow_the_tu_row(
    value: float, lo: float, width: float, tu_q: float | None, gene_q: float
) -> None:
    style.apply()
    cfg = resolve_config()
    hi = round(lo + width, 4)
    lo = round(lo, 4)
    over = {
        ("morans_i", "gene_label", 1): (value, 0.0, -0.01, 0.01, gene_q),
        ("morans_i", "tu_block", 1): (value, round((lo + hi) / 2, 4), lo, hi, tu_q),
    }
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(_out(Path(tmp), correlogram_rows("chrA", overrides=over)), cfg)
    sig = tu_q is not None and tu_q <= cfg.stats.alpha.context
    assert _marker_at(_ax(inst, "a:chrA"), 1, value) == [not sig]
    assert bool(_row(inst, "chrA", "morans_i", "tu_block", 1)["significant"]) is sig
    word = "above" if value > hi else "below" if value < lo else "within"
    if tu_q is not None:
        assert inst.title == (
            f"Lag-1 Moran's I {value:g}: {word} the TU-block envelope (q = {tu_q:g})"
        )


# --- style, files and determinism -------------------------------------------------------------


def _rich(tmp: Path) -> Path:
    rows = [*correlogram_rows("chrA", overrides=EX4), *correlogram_rows("chrB", cell=default_cell)]
    return _out(tmp, rows, (CHR_A, CHR_B))


def test_every_text_is_at_least_6_pt(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_rich(tmp_path), cfg)
    sizes = [t.get_fontsize() for t in inst.figure.findobj(Text) if t.get_text()]
    assert sizes
    assert min(sizes) >= 6


def test_no_text_overlaps_after_an_agg_draw(tmp_path: Path, cfg: GippaConfig) -> None:
    assert text_overlaps(_render(_rich(tmp_path), cfg).figure) == []
    rows = [*na_correlogram_rows("chrA"), *na_correlogram_rows("chrB")]
    assert text_overlaps(_render(_out(tmp_path / "na", rows, (CHR_A, CHR_B)), cfg).figure) == []


def test_data_tsv_and_images_are_byte_identical(tmp_path: Path, cfg: GippaConfig) -> None:
    out = _rich(tmp_path / "run")
    written = []
    for name in ("x", "y"):
        inst = _render(out, cfg)
        paths = write_instance(
            inst, "V10", "correlogram", tmp_path / name, ["svg", "png", "pdf"], 100
        )
        written.append({p.relative_to(tmp_path / name): p.read_bytes() for p in paths})
    assert set(written[0]) == {
        Path("V10_correlogram.svg"),
        Path("V10_correlogram.png"),
        Path("V10_correlogram.pdf"),
        Path("data/V10_correlogram.tsv"),
    }
    assert written[0] == written[1]
    tsv = written[0][Path("data/V10_correlogram.tsv")].decode("utf-8").splitlines()
    assert tsv[0] == "\t".join(COLUMNS)
    assert (
        "chrA\tmorans_i\ttu_block\t1\t0.31\t-0.004\t-0.05\t0.12\t0.314\t0.002\ttrue\t"
        in "\n".join(tsv)
    )


def test_render_does_not_mutate_between_calls(tmp_path: Path, cfg: GippaConfig) -> None:
    view = RunView.from_outdir(_rich(tmp_path))
    first = plot_spec("V10").fn(view, cfg)
    second = plot_spec("V10").fn(view, cfg)
    assert isinstance(first, list) and isinstance(second, list)
    pd.testing.assert_frame_equal(first[0].data, second[0].data)
    assert first[0].title == second[0].title


VAL_OUT = Path(os.environ.get("GIPPA_ROOT", REPO / "local")) / "runs" / "val-mvp-2026-10-22" / "OUT"


@pytest.mark.slow
@pytest.mark.skipif(not VAL_OUT.is_dir(), reason=f"validation run absent: {VAL_OUT}")
def test_example6_validation_run_is_not_estimable(tmp_path: Path, cfg: GippaConfig) -> None:
    """chr1 and plasmid_1: 160 morans_i / join_count_bb rows, all NA; oriter_slope dropped."""
    inst = _render(VAL_OUT, cfg)
    assert inst.title == NOT_ESTIMABLE
    assert len(inst.data) == 160
    for contig in ("chr1", "plasmid_1"):
        notes = [t for ax in _contig_axes(inst, contig) for t in _notes(ax)]
        assert any(t.startswith("not estimable (") for t in notes), (contig, notes)
    paths = write_instance(inst, "V10", "correlogram", tmp_path, ["svg", "png"], 100)
    assert all(p.stat().st_size > 0 for p in paths)
