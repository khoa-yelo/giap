"""V3 / V13 locus views with coordinate-less regions (regression: MVP-04 K-12 run, Slurm 1884567,
``REPORT:PLOTS`` failed with ``PlotFailed: plot V03 failed: TypeError: int() argument must be ...
not 'NAType'``).

Cause: ``annotation/regions.tsv`` writes AMR gene sets with ``NA`` start/end by design
(GIPPA_SPEC.md §7.1.7: ``amr`` is "a gene set rather than a contiguous region";
``gippa.annotate.regions`` writes ``NA`` for them). ``RunView`` reads them as nullable ``Int64``
``<NA>``. The region shading in ``gippa.viz.locus.render_locus`` passes every ``regions.tsv`` row of
the locus's contig through ``_in_window``, which calls ``int(r.start)`` before any overlap test, so
any locus on a contig that carries an AMR set raised. V13 renders its locus views with the same
``render_locus``, so it fails the same way.

Expected behaviour (M9-08 lead decision 5 for V01, applied to V03 / V13 shading; reversible
default, reported to the lead): a ``regions.tsv`` row with NA ``start`` or ``end`` is not shaded
and gets no ``region`` data row; the locus still renders, and regions with coordinates are shaded
as before. ``regions_stats`` prophage / BGC rows with NA coordinates (not expected upstream) are
left out of the V3 pool rather than failing the plot.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest
from helpers import viz_locus as vl
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.config import GippaConfig, resolve_config
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance
from gippa.viz.locus import Locus, OrfMark, render_locus, select_loci

SCALE = (18.0, 24.0)


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


def _amr_set(region_id: str = "amr_1", contig: str = "chr1") -> tuple:
    """An ``annotation/regions.tsv`` AMR gene set: NA start and end (§7.1.7)."""
    return (region_id, "amr", contig, None, None, "AMRFinderPlus", None, "g_amr1;g_amr2", "")


def _out(tmp: Path, regions: list[tuple]) -> Path:
    return vl.build_out(
        tmp,
        genes=[vl.gene("G1", "chr1", 10_101, 11_000, log2=20.0)],
        regions_stats=[vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 12_000, 0.01)],
        regions=regions,
    )


def _region_rows(inst: Instance) -> list[tuple[str, int, int, str]]:
    d = inst.data[inst.data["track"] == "region"]
    return sorted((str(r.id), int(r.start), int(r.end), str(r.category)) for r in d.itertuples())


def _render(out: Path, cfg: GippaConfig, locus_id: str) -> Instance:
    view = RunView.from_outdir(out)
    (lo,) = [x for x in select_loci(view, cfg) if x.locus_id == locus_id]
    return render_locus(lo, view, cfg, SCALE)


def test_amr_set_without_coordinates_is_not_shaded(tmp_path: Path, cfg: GippaConfig) -> None:
    out = _out(tmp_path, [vl.region("prophage_1", "prophage", "chr1", 10_001, 12_000), _amr_set()])
    inst = _render(out, cfg, "prophage_1")
    assert _region_rows(inst) == [("prophage_1", 10_000, 12_000, "prophage")]


def test_amr_set_alone_on_the_contig_renders(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path, [_amr_set()]), cfg, "prophage_1")
    assert _region_rows(inst) == []
    assert ("gene", "G1") in {(str(r.track), str(r.id)) for r in inst.data.itertuples()}


def test_amr_set_on_another_contig_is_ignored(tmp_path: Path, cfg: GippaConfig) -> None:
    out = _out(tmp_path, [_amr_set("amr_2", "lin1"), _amr_set("amr_1", "chr1")])
    assert _region_rows(_render(out, cfg, "prophage_1")) == []


def test_region_with_only_end_missing_is_not_shaded(tmp_path: Path, cfg: GippaConfig) -> None:
    half = ("crispr_1", "crispr", "chr1", 10_501, None, "Bakta", None, "", "")
    out = _out(tmp_path, [half, vl.region("prophage_1", "prophage", "chr1", 10_001, 12_000)])
    assert [r[0] for r in _region_rows(_render(out, cfg, "prophage_1"))] == ["prophage_1"]


def test_v03_plot_runs_with_an_amr_set(tmp_path: Path, cfg: GippaConfig) -> None:
    """The failing call of the MVP-04 run: the registered V03 function, not only the renderer."""
    import gippa.viz.locus  # noqa: F401  (registers V03)

    out = _out(tmp_path, [vl.region("prophage_1", "prophage", "chr1", 10_001, 12_000), _amr_set()])
    result = base.registered()["V03"].fn(RunView.from_outdir(out), cfg)
    assert isinstance(result, list)
    assert [i.suffix for i in result] == ["prophage_1", None]


def test_v13_locus_view_with_an_amr_set(tmp_path: Path, cfg: GippaConfig) -> None:
    """V13 draws its novel ORF views with ``render_locus`` (``Locus(kind="novel_orf")``)."""
    out = _out(tmp_path, [_amr_set()])
    lo = Locus("ORF_1", "novel_orf", "chr1", 10_400, 10_700, 0.001, "ORF_1")
    inst = render_locus(lo, RunView.from_outdir(out), cfg, SCALE, orf=OrfMark("+", 10_400))
    assert _region_rows(inst) == []
    assert "novel_orf" in set(inst.data["track"])


def test_prophage_stat_row_without_coordinates_is_left_out(tmp_path: Path, cfg: GippaConfig) -> None:
    """Proposed default: a prophage / BGC ``regions_stats`` row with NA start/end has no window,
    so it is not a V3 locus; the other loci are unaffected."""
    out = vl.build_out(
        tmp_path,
        genes=[vl.gene("G1", "chr1", 10_101, 11_000, log2=20.0)],
        regions_stats=[
            vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 12_000, 0.01),
            vl.region_stat("bgc_1", "bgc", "chr1", None, None, 0.02),  # type: ignore[arg-type]
            vl.region_stat("amr_1", "amr", "chr1", None, None, 0.3),  # type: ignore[arg-type]
        ],
    )
    loci = select_loci(RunView.from_outdir(out), cfg)
    assert [lo.locus_id for lo in loci] == ["prophage_1"]


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    spans=st.lists(
        st.one_of(
            st.none(),
            st.tuples(st.integers(1, 99_000), st.integers(1, 900)).map(
                lambda t: (t[0], t[0] + t[1])
            ),
        ),
        min_size=1,
        max_size=6,
    )
)
def test_region_rows_are_exactly_the_coordinate_bearing_overlaps(
    spans: list[tuple[int, int] | None],
) -> None:
    """Property: whatever mix of NA and real regions, the plot renders and its ``region`` rows are
    the regions with coordinates that overlap the window (8 001-14 000, 1-based). Every span
    lies inside chr1 (100 kb) and far from the origin, so no +/- length shift can reach the
    window."""
    cfg = resolve_config()
    regions = [
        _amr_set(f"amr_{i}") if s is None else vl.region(f"crispr_{i}", "crispr", "chr1", *s)
        for i, s in enumerate(spans)
    ]
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(_out(Path(tmp), regions), cfg, "prophage_1")
    lo, hi = 10_000 - 2_000, 12_000 + 2_000  # 0-based half-open window of prophage_1
    want = sorted(
        (f"crispr_{i}", s[0] - 1, s[1], "crispr")
        for i, s in enumerate(spans)
        if s is not None and s[0] - 1 < hi and s[1] > lo
    )
    assert _region_rows(inst) == want
