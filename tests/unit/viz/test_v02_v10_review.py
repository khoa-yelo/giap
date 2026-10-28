"""M9-09 review round (PR #218): V2 ``abundance_landscape`` and V10 ``correlogram``.

Pins the lead decisions on the review of 1f5fba6 (docs/specs/M9-09-v02-v10-landscape-correlogram.md)
(numbered as there):

- 21: an ``autocorrelation.tsv`` row flagged ``not_estimable`` or ``few_values`` is NA whatever its
  value (the producer writes the join count as 0 with zero nulls): it is not drawn and counts as NA
  for the "not estimable (<flags>)" empty state and the title;
- 22: V2's legend overlaps neither the title nor the data, with the oriC/ter entries from
  ``plots/support/replicon_landmarks.tsv`` (M9-03 columns), and everything sits inside the figure;
- 23: V10 wraps its contig columns into rows of at most 3, the height growing per row;
- 25: V10 rows sort ``statistic`` in the schema's enum order (lead decision 6; mutant C7);

plus the cheap NIT pins of the review (mutants L7, L11, L12, L14, C4, C13, C9).

Oracles are the fixture rows and ``gippa.stats.smooth.loess_detrend`` (the spec's example 1 oracle).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from helpers.viz_landscape import (
    CHR1,
    Cell,
    Replicon,
    correlogram_rows,
    gene,
    island,
    landscape_out,
    markers,
    na_correlogram_rows,
    plot_spec,
    text_overlaps,
)
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.text import Text

from gippa.config import GippaConfig, resolve_config
from gippa.io.plot_files import write_instance
from gippa.io.run_view import RunView
from gippa.stats.smooth import loess_detrend
from gippa.viz import style
from gippa.viz.base import Instance, Skip

NOT_ESTIMABLE = "Autocorrelation not estimable on this genome"
CHR_A: Replicon = ("chrA", 100_000, "circular", "chromosome")
CHR_B: Replicon = ("chrB", 60_000, "circular", "chromosome")
NE_FLAGS = "ci_unestimated;not_estimable"
FRAC1 = "context.islands.detrend_span_frac=1.0"
LANDMARK_HEADER = "contig\tlength_bp\ttopology\toric_position\tter_position\tstatus\n"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


@pytest.fixture
def cfg1() -> GippaConfig:
    return resolve_config(overrides=[FRAC1])


def _render(pid: str, out: Path, cfg: GippaConfig) -> Instance:
    result = plot_spec(pid).fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert isinstance(inst, Instance)
    return inst


def _ax(inst: Instance, label: str) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_label() == label]
    assert len(found) == 1, f"axes {label!r}: {[a.get_label() for a in inst.figure.axes]}"
    return found[0]


def _notes(ax: Axes) -> list[str]:
    return [t.get_text() for t in ax.texts if t.get_text()]


def _v10_out(tmp: Path, rows: list[dict], replicons: tuple = (CHR_A,)) -> Path:
    return landscape_out(tmp, replicons=replicons, genes=[gene("g1", 10_000, 1.0)], autocorr=rows)


def _v02_tsv(inst: Instance, tmp: Path) -> pd.DataFrame:
    """The written ``V02_abundance_landscape.tsv`` (1-based closed)."""
    paths = write_instance(inst, "V02", "abundance_landscape", tmp, ["svg"], 72)
    tsv = next(p for p in paths if p.suffix == ".tsv")
    return pd.read_csv(tsv, sep="\t", na_values=["NA"], keep_default_na=False)


def _kind(frame: pd.DataFrame, kind: str) -> pd.DataFrame:
    return frame[frame["kind"].astype(str) == kind].reset_index(drop=True)


def _drawn_outside_figure(fig: Figure) -> list[str]:
    """Visible texts and legends whose rendered box leaves the figure by more than 1 px."""
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    renderer = canvas.get_renderer()
    w, h = fig.bbox.width, fig.bbox.height
    out = []
    items: list = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    items += _legends(fig)
    for item in items:
        b = item.get_window_extent(renderer)
        if b.x0 < -1 or b.y0 < -1 or b.x1 > w + 1 or b.y1 > h + 1:
            out.append(item.get_text() if isinstance(item, Text) else repr(item))
    return out


# --- decision 21: not-estimable rows are NA whatever their value ------------------------------


def _flag_lags(rows: list[dict], statistic: str, lags: set[int], flags: str, zero: bool) -> None:
    """Rewrite the (statistic, lag ∈ lags) rows of both nulls as the producer writes a
    not-estimable row: the given flags, q NA and, when ``zero``, value and nulls 0."""
    for r in rows:
        if r["statistic"] == statistic and r["lag"] in lags:
            r["flags"] = flags
            r["q_value"] = r["p_value"] = None
            if zero:
                r.update(value=0.0, null_mean=0.0, null_lo=0.0, null_hi=0.0, effect=0.0)


def test_all_not_estimable_contig_with_zero_join_counts_is_the_empty_state(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    """chrB as the producer writes it (join count 0, zero nulls, q NA, ``not_estimable``): both
    panels carry the note and no marker, beside an estimable chrA (review B1)."""
    rows = [*correlogram_rows("chrA"), *na_correlogram_rows("chrB")]
    jc = [r for r in rows if r["contig"] == "chrB" and r["statistic"] == "join_count_bb"]
    assert jc and all(r["value"] == 0.0 and r["null_hi"] == 0.0 for r in jc)  # finite zeros
    inst = _render("V10", _v10_out(tmp_path, rows, (CHR_A, CHR_B)), cfg)
    for panel in ("a", "b"):
        ax = _ax(inst, f"{panel}:chrB")
        notes = _notes(ax)
        assert any(t.startswith("not estimable (") for t in notes), (panel, notes)
        assert markers(ax) == [], panel
    assert len(markers(_ax(inst, "b:chrA"))) == 20


def test_flagged_rows_are_na_even_when_every_value_is_zero(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    """Moran's I written as 0 too: flags decide, not finiteness; the title is not estimable."""
    rows = correlogram_rows("chrA")
    _flag_lags(rows, "morans_i", set(range(1, 21)), "ci_unestimated;few_values;not_estimable", True)
    _flag_lags(rows, "join_count_bb", set(range(1, 21)), NE_FLAGS, True)
    inst = _render("V10", _v10_out(tmp_path, rows), cfg)
    assert inst.title == NOT_ESTIMABLE
    for panel in ("a", "b"):
        ax = _ax(inst, f"{panel}:chrA")
        assert markers(ax) == [], panel
        assert any(t.startswith("not estimable (") for t in _notes(ax)), panel
    assert not inst.data["significant"].astype(bool).any()


@pytest.mark.parametrize("flags", [NE_FLAGS, "ci_unestimated;few_values"])
def test_mixed_contig_draws_only_the_estimable_lags(
    tmp_path: Path, cfg: GippaConfig, flags: str
) -> None:
    """Join count lags 11-20 written as 0 and flagged: markers at lags 1-10 only, no note; Moran's
    I lags 15-20 flagged with their finite values kept: markers at lags 1-14 only."""
    rows = correlogram_rows("chrA")
    _flag_lags(rows, "join_count_bb", set(range(11, 21)), flags, True)
    _flag_lags(rows, "morans_i", set(range(15, 21)), flags, False)
    inst = _render("V10", _v10_out(tmp_path, rows), cfg)
    b = _ax(inst, "b:chrA")
    assert sorted(x for x, _, _, _ in markers(b)) == [float(lag) for lag in range(1, 11)]
    assert all(y != 0.0 for _, y, _, _ in markers(b))
    assert not any(t.startswith("not estimable") for t in _notes(b))
    a = _ax(inst, "a:chrA")
    assert sorted(x for x, _, _, _ in markers(a)) == [float(lag) for lag in range(1, 15)]


def test_a_flagged_lag1_row_is_not_the_title_row(tmp_path: Path, cfg: GippaConfig) -> None:
    """chrA's Moran's I is all flagged (written as 0): the title falls through to chrB."""
    rows = [*correlogram_rows("chrA"), *correlogram_rows("chrB")]
    for r in rows:
        if r["contig"] == "chrA" and r["statistic"] == "morans_i":
            r.update(flags=NE_FLAGS, q_value=None, value=0.0, null_mean=0.0, null_lo=0.0)
            r.update(null_hi=0.0)
    inst = _render("V10", _v10_out(tmp_path, rows, (CHR_A, CHR_B)), cfg)
    # chrB lag 1: 0.16 above its TU envelope (-0.094, 0.086), q 0.5 (default_cell)
    assert inst.title == "Lag-1 Moran's I 0.16: above the TU-block envelope (q = 0.5)"


# --- decision 22: V2 legend with oriC / ter entries -------------------------------------------


def _landmark_out(tmp: Path) -> Path:
    """The reviewer's case in miniature: a 2.5 Mb circular chromosome with oriC/ter called, a
    60 kb linear plasmid, islands (one crossing the origin), all four legend statuses."""
    rng = np.random.default_rng(7)
    reps: list[Replicon] = [
        ("chr1", 2_500_000, "circular", "chromosome"),
        ("linear_plasmid_pXYZ_long", 60_000, "linear", "plasmid"),
    ]
    statuses = ["detected", "detected", "detected", "not_detected", "single_peptide_unestimated"]
    genes = []
    for i in range(120):
        mid = 10_000 + i * 20_800
        y = 2 + 1.5 * np.cos(2 * np.pi * mid / 2_500_000) + float(rng.normal(0, 0.6))
        genes.append(gene(f"c{i:04d}", mid, float(y), status=statuses[i % 5]))
    for i in range(30):
        mid = 1_000 + i * 1_900
        status = "not_detected" if i % 5 == 0 else "detected"
        genes.append(
            gene(f"p{i:02d}", mid, 1.0 + i / 10, contig=reps[1][0], status=status, half=300)
        )
    islands = [
        island("ISL_1", "high", 2_400_001, 2_600_000, 0.001),
        island("ISL_2", "low", 1_200_001, 1_300_000, 0.01),
        island("ISL_3", "high", 10_001, 20_000, 0.02, contig=reps[1][0]),
    ]
    over: dict[tuple[str, str, int], Cell] = {
        ("morans_i", "gene_label", 1): (0.31, -0.002, -0.05, 0.05, 0.001),
    }
    ac = [*correlogram_rows("chr1", overrides=over), *correlogram_rows(reps[1][0])]
    out = landscape_out(tmp, replicons=reps, genes=genes, islands=islands, autocorr=ac)
    sup = out / "plots" / "support"
    sup.mkdir(parents=True, exist_ok=True)
    (sup / "replicon_landmarks.tsv").write_text(
        LANDMARK_HEADER
        + "chr1\t2500000\tcircular\t1\t1250001\tcalled\n"
        + f"{reps[1][0]}\t60000\tlinear\tNA\tNA\tlinear\n",
        encoding="utf-8",
    )
    return out


def test_landmarks_are_drawn_from_the_support_file(tmp_path: Path, cfg: GippaConfig) -> None:
    """Positive control for the layout tests: the oriC and ter rows and legend entries exist."""
    inst = _render("V02", _landmark_out(tmp_path), cfg)
    kinds = set(inst.data["kind"].astype(str))
    assert {"oric", "ter"} <= kinds
    labels = [t.get_text() for lg in _legends(inst.figure) for t in lg.get_texts()]
    assert {"oriC", "ter"} <= set(labels), labels


def _legends(fig: Figure) -> list:  # type: ignore[type-arg]
    return [*fig.legends, *[ax.get_legend() for ax in fig.axes if ax.get_legend()]]


def test_legend_with_landmarks_overlaps_no_text(tmp_path: Path, cfg: GippaConfig) -> None:
    """Review B2: with the oriC/ter entries the legend ran into the suptitle."""
    inst = _render("V02", _landmark_out(tmp_path), cfg)
    assert text_overlaps(inst.figure) == []


def test_legend_with_landmarks_stays_off_the_title_and_the_data(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    """The legend box intersects neither the suptitle's box nor any panel's data area
    (decision 22: below the axes or in reserved space)."""
    inst = _render("V02", _landmark_out(tmp_path), cfg)
    fig = inst.figure
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    renderer = canvas.get_renderer()
    legends = _legends(fig)
    assert legends
    suptitle = fig._suptitle
    assert suptitle is not None and suptitle.get_text() == inst.title
    boxes = [("title", suptitle.get_window_extent(renderer))]
    boxes += [(ax.get_label(), ax.get_window_extent(renderer)) for ax in fig.axes]
    for lg in legends:
        lb = lg.get_window_extent(renderer)
        for name, b in boxes:
            w = min(lb.x1, b.x1) - max(lb.x0, b.x0)
            h = min(lb.y1, b.y1) - max(lb.y0, b.y0)
            assert not (w > 1.0 and h > 1.0), f"legend overlaps {name}"


def test_landmark_figure_keeps_every_text_inside(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render("V02", _landmark_out(tmp_path), cfg)
    assert _drawn_outside_figure(inst.figure) == []


# --- decision 23: V10 with many replicons -----------------------------------------------------


def _many(n: int) -> tuple[Replicon, ...]:
    return (
        ("chr1", 3_000_000, "circular", "chromosome"),
        *[(f"plasmid_{i}", 50_000, "circular", "plasmid") for i in range(1, n)],
    )


def _many_inst(tmp: Path, n: int, cfg: GippaConfig) -> Instance:
    reps = _many(n)
    rows = [r for rep in reps for r in correlogram_rows(rep[0])]
    return _render("V10", _v10_out(tmp, rows, reps), cfg)


@pytest.mark.parametrize("n", [4, 6])
def test_many_replicons_wrap_into_rows_of_at_most_3(
    tmp_path: Path, cfg: GippaConfig, n: int
) -> None:
    inst = _many_inst(tmp_path, n, cfg)
    reps = [r[0] for r in _many(n)]
    a_axes = [_ax(inst, f"a:{c}") for c in reps]
    for c in reps:
        _ax(inst, f"b:{c}")
    columns = {round(ax.get_position().x0, 3) for ax in a_axes}
    assert len(columns) <= 3, columns
    rows: dict[float, int] = {}
    for ax in a_axes:
        key = round(ax.get_position().y0, 3)
        rows[key] = rows.get(key, 0) + 1
    assert max(rows.values()) <= 3, rows
    assert len(rows) == -(-n // 3), rows
    # reading order (rows top to bottom, then left to right) is FASTA order
    reading = sorted(a_axes, key=lambda ax: (-ax.get_position().y0, ax.get_position().x0))
    assert [ax.get_label() for ax in reading] == [f"a:{c}" for c in reps]
    # each contig's (b) panel sits below its (a) panel in the same column
    for c in reps:
        a, b = _ax(inst, f"a:{c}").get_position(), _ax(inst, f"b:{c}").get_position()
        assert b.y0 < a.y0 and abs(b.x0 - a.x0) < 1e-6, c


@pytest.mark.parametrize("n", [4, 6])
def test_many_replicons_have_no_text_overlaps(tmp_path: Path, cfg: GippaConfig, n: int) -> None:
    inst = _many_inst(tmp_path, n, cfg)
    assert text_overlaps(inst.figure) == []
    assert _drawn_outside_figure(inst.figure) == []


def test_figure_height_grows_with_the_rows(tmp_path: Path, cfg: GippaConfig) -> None:
    three = _many_inst(tmp_path / "3", 3, cfg).figure.get_size_inches()[1]
    four = _many_inst(tmp_path / "4", 4, cfg).figure.get_size_inches()[1]
    six = _many_inst(tmp_path / "6", 6, cfg).figure.get_size_inches()[1]
    assert four > three
    assert six == pytest.approx(four)


# --- decision 25: statistic order -------------------------------------------------------------


def test_statistics_sort_in_enum_order_not_alphabetically(tmp_path: Path, cfg: GippaConfig) -> None:
    """Lead decision 6: morans_i before join_count_bb in each contig (mutant C7)."""
    rows = [*correlogram_rows("chrA"), *correlogram_rows("chrB")][::-1]
    inst = _render("V10", _v10_out(tmp_path, rows, (CHR_A, CHR_B)), cfg)
    paths = write_instance(inst, "V10", "correlogram", tmp_path / "plots", ["svg"], 72)
    tsv = pd.read_csv(next(p for p in paths if p.suffix == ".tsv"), sep="\t", na_values=["NA"])
    for frame in (inst.data, tsv):
        for contig in ("chrA", "chrB"):
            stats = list(frame.loc[frame["contig"].astype(str) == contig, "statistic"].astype(str))
            assert stats == ["morans_i"] * 40 + ["join_count_bb"] * 40, contig
    keys = list(zip(tsv["contig"], tsv["statistic"], tsv["null_model"], strict=True))
    blocks = list(dict.fromkeys(keys))
    assert [s for _, s, _ in blocks[:4]] == ["morans_i"] * 2 + ["join_count_bb"] * 2


# --- NIT pins ---------------------------------------------------------------------------------


def test_rows_sharing_a_midpoint_are_averaged_for_the_trend(
    tmp_path: Path, cfg1: GippaConfig
) -> None:
    """Two genes at midpoint 30,000 with log2 1 and 3 enter the LOESS once, at 2 (mutant L7)."""
    genes = [
        gene("g1", 10_000, 1.0),
        gene("g2a", 30_000, 1.0),
        gene("g2b", 30_000, 3.0, half=300),
        gene("g3", 50_000, 3.0),
        gene("g4", 70_000, 2.0),
        gene("g5", 90_000, 1.0),
    ]
    inst = _render("V02", landscape_out(tmp_path, genes=genes, islands=[]), cfg1)
    assert len(_kind(inst.data, "point")) == 6
    trend = _kind(inst.data, "trend")
    by_pos = trend.groupby("position_bp", sort=True)["value"].first()
    assert len(by_pos) == 5
    xs = [10_000, 30_000, 50_000, 70_000, 90_000]
    engine = loess_detrend([1.0, 2.0, 3.0, 2.0, 1.0], xs, 1.0, circular=True, period=100_000).fit
    np.testing.assert_allclose(by_pos.to_numpy(dtype=float), engine, atol=1e-12)


def test_origin_crossing_island_position_wraps_mod_length(
    tmp_path: Path, cfg1: GippaConfig
) -> None:
    """ISL 95,001-105,000 on L = 100,000: (start0 + end) / 2 mod L + 1 = 1 (mutant L14)."""
    genes = [gene(f"g{i}", m, 1.0) for i, m in enumerate((10_000, 30_000, 50_000))]
    islands = [island("ISL_9", "high", 95_001, 105_000, 0.001)]
    inst = _render("V02", landscape_out(tmp_path, genes=genes, islands=islands), cfg1)
    isl = _kind(_v02_tsv(inst, tmp_path / "plots"), "island")
    assert len(isl) == 1
    assert int(isl.loc[0, "position_bp"]) == 1
    assert 1 <= int(isl.loc[0, "position_bp"]) <= CHR1[1]


def test_bins_take_the_median_of_several_genes(tmp_path: Path, cfg: GippaConfig) -> None:
    """20,001 genes, three per 10 bp bin with log2 0, 0, 9: every bin's point is 0, the median,
    not 3, the mean (mutant L11)."""
    n_bins, length = 6_667, 200_000
    rep: Replicon = ("chrM", length, "circular", "chromosome")
    genes = []
    for b in range(n_bins):
        for k, y in ((2, 0.0), (4, 0.0), (6, 9.0)):
            genes.append(gene(f"m{b:05d}_{k}", 10 * b + k, y, contig="chrM", half=1))
    assert len(genes) == 20_001
    inst = _render("V02", landscape_out(tmp_path, replicons=[rep], genes=genes, islands=[]), cfg)
    pts = _kind(inst.data, "point")
    assert len(pts) == n_bins
    assert [float(v) for v in pts["value"]] == [0.0] * n_bins


def test_exactly_2_mb_stays_in_kb(tmp_path: Path, cfg1: GippaConfig) -> None:
    """Mb only when length > 2 Mb (mutant L12)."""
    rep: Replicon = ("two", 2_000_000, "circular", "chromosome")
    genes = [gene(f"g{i}", m, float(i), contig="two") for i, m in enumerate((500_000, 1_500_000))]
    inst = _render("V02", landscape_out(tmp_path, replicons=[rep], genes=genes, islands=[]), cfg1)
    assert _ax(inst, "two").get_xlabel() == "position (kb)"


@pytest.mark.parametrize(("value", "lo", "hi"), [(0.12, -0.05, 0.12), (-0.05, -0.05, 0.12)])
def test_value_on_the_envelope_edge_is_within(
    tmp_path: Path, cfg: GippaConfig, value: float, lo: float, hi: float
) -> None:
    """``value == null_hi`` and ``value == null_lo`` read "within" (mutants C4, C13)."""
    over: dict[tuple[str, str, int], Cell] = {
        ("morans_i", "gene_label", 1): (value, -0.002, -0.01, 0.01, 0.001),
        ("morans_i", "tu_block", 1): (value, -0.004, lo, hi, 0.04),
    }
    inst = _render("V10", _v10_out(tmp_path, correlogram_rows("chrA", overrides=over)), cfg)
    assert inst.title == f"Lag-1 Moran's I {value:g}: within the TU-block envelope (q = 0.04)"


def test_marker_falls_back_to_the_gene_label_value(tmp_path: Path, cfg: GippaConfig) -> None:
    """TU-block lag 3 has value NA (not flagged not-estimable); the gene-label row's 0.2 is drawn,
    open because the TU-block q is NA (mutant C9)."""
    na: Cell = (None, None, None, None, None)
    over: dict[tuple[str, str, int], Cell] = {
        ("morans_i", "tu_block", 3): na,
        ("morans_i", "gene_label", 3): (0.2, -0.002, -0.05, 0.05, 0.001),
    }
    inst = _render("V10", _v10_out(tmp_path, correlogram_rows("chrA", overrides=over)), cfg)
    at3 = [(y, hollow) for x, y, hollow, _ in markers(_ax(inst, "a:chrA")) if x == 3.0]
    assert at3 == [(pytest.approx(0.2), True)]
