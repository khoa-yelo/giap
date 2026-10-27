"""M9-26 example 1 and work item 2: the 20-Mb stress fixture (25,000 genes on one contig) bins
V1's gene tracks and V2's points, and every track stays at or under 20,000 drawn elements
(docs/specs/M9-26-scale-and-determinism.md; GIPPA_SPEC.md §10.1 "Plots aggregate into bins when
a track would draw more than 20,000 elements"; M9-index decision 12; M9-24 checklist item
"≤ 20,000 drawn elements", which for V1's polar axes counts per track, M9-24 "V1 polar axes").

The fixture is ``helpers.scale_fixture.stress_out`` (gippasim layout coordinates, ``RunView``
frames only, not a pipeline run). Oracles are recomputed from the written ``genes.tsv``:

- V1 bins genes by their 0-based ``start`` over ``gippa.viz.base.bin_edges(L, 20000)`` (one
  replicon, so it gets every bin; M9-08's ``_bins`` contract); ``gene_plus``/``gene_minus``
  values are the detected share of that strand's genes in the bin, ``abundance`` the median of
  the measured ``abundance_log2`` of all genes in the bin (``None`` when none is measured);
- V2 bins the measured detected genes by their 0-based midpoint ``(start0 + end) / 2 mod L``
  over the same edges and draws one point per non-empty bin at the bin centre with the median
  ``abundance_log2``; undetected genes stay a rug (2,500 rows here, under the limit).
- Origin wrap (lead decision 4, review N7): the fixture's two origin-crossing genes are stored
  unwrapped (``length < end < 2 length``). V1 puts each in the bin of its 0-based start (the last
  bin); V2 puts each at the bin of its midpoint mod L: one before the origin (the last bin), one
  past it (the first bin).

The whole module renders V1 and V2 once at 25,000 genes (≈ 20-40 s): marked ``slow``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.scale_fixture import MAX_ELEMENTS, N_GENES, StressRun, stress_out
from helpers.viz_tracks import elements

from gippa.config import resolve_config
from gippa.io.plot_files import write_instance
from gippa.io.run_view import RunView
from gippa.viz import base

pytestmark = pytest.mark.slow

ELEMENT_ITEM = "≤ 20,000 drawn elements"
GENE_TRACKS = ("gene_plus", "gene_minus", "abundance")


@pytest.fixture(scope="module")
def stress(tmp_path_factory: pytest.TempPathFactory) -> StressRun:
    return stress_out(tmp_path_factory.mktemp("stress20mb"))


@pytest.fixture(scope="module")
def genes(stress: StressRun) -> pd.DataFrame:
    return RunView.from_outdir(stress.out).table("genes_tsv")  # 0-based starts


def _render(stress: StressRun, pid: str) -> base.Instance:
    from gippa.steps.plots import load_plot_modules
    from gippa.viz import style

    load_plot_modules()
    style.apply()
    result = base.registered()[pid].fn(RunView.from_outdir(stress.out), resolve_config(None))
    assert not isinstance(result, base.Skip), result
    (inst,) = [i for i in result if i.suffix is None]
    return inst


@pytest.fixture(scope="module")
def v01(stress: StressRun) -> base.Instance:
    return _render(stress, "V01")


@pytest.fixture(scope="module")
def v02(stress: StressRun) -> base.Instance:
    return _render(stress, "V02")


def _edges(stress: StressRun) -> np.ndarray:
    return base.bin_edges(stress.length, min(MAX_ELEMENTS, stress.length))


def _bin_of(edges: np.ndarray, pos: np.ndarray) -> np.ndarray:
    return np.searchsorted(edges, pos, side="right") - 1


def _na(v: Any) -> bool:
    return v is None or (isinstance(v, float) and np.isnan(v)) or v is pd.NA


# ---------------------------------------------------------------- the fixture


def test_stress_fixture_is_25000_genes_on_one_20_mb_contig(stress: StressRun) -> None:
    assert stress.n_genes == N_GENES
    assert 19_000_000 <= stress.length <= 22_000_000
    # V2 must have more than 20,000 points on its one panel for its points to bin
    assert stress.n_points > MAX_ELEMENTS
    # neither strand alone exceeds the limit: V1 bins on the gene count of the plot
    assert stress.n_plus <= MAX_ELEMENTS and stress.n_minus <= MAX_ELEMENTS


def _crossing(stress: StressRun, genes: pd.DataFrame) -> pd.DataFrame:
    g = genes[genes["locus_tag"].astype(str).isin(stress.crossing)]
    assert len(g) == 2, stress.crossing
    return g


def test_stress_fixture_has_two_origin_crossing_genes_stored_unwrapped(
    stress: StressRun, genes: pd.DataFrame
) -> None:
    g = _crossing(stress, genes)
    assert set(g["crosses_origin"].astype(bool)) == {True}
    s, e, L = g["start"].astype(int), g["end"].astype(int), stress.length
    assert ((s >= 0) & (s < L) & (e > L) & (e < 2 * L)).all(), g[["start", "end"]]
    assert (g["detection_status"] == "detected").all()
    assert int(genes["crosses_origin"].astype(bool).sum()) == 2
    # one midpoint before the origin, one past it (wrapped by mod L)
    mids = np.sort(np.mod((s + e) / 2.0, float(L)).to_numpy())
    edges = _edges(stress)
    assert _bin_of(edges, mids).tolist() == [0, len(edges) - 2], mids


# ---------------------------------------------------------------- V1


@pytest.mark.parametrize("track", GENE_TRACKS)
def test_example_1_v1_gene_tracks_binned_with_at_most_20000_rows(
    v01: base.Instance, track: str
) -> None:
    rows = v01.data[v01.data["track"] == track]
    assert len(rows) > 0, f"no {track} rows"
    assert len(rows) <= MAX_ELEMENTS, f"{track}: {len(rows)} rows"
    assert set(rows["binned"].astype(bool)) == {True}, f"{track} rows not all binned"


@pytest.mark.parametrize("track", GENE_TRACKS)
def test_v1_bins_are_disjoint_sorted_and_on_the_edge_grid(
    v01: base.Instance, stress: StressRun, track: str
) -> None:
    rows = v01.data[v01.data["track"] == track]
    s = rows["start"].astype(int).to_numpy()
    e = rows["end"].astype(int).to_numpy()
    assert np.all(np.diff(s) > 0) and np.all(e[:-1] <= s[1:]), f"{track} bins overlap"
    grid = set(_edges(stress).tolist())
    assert set(s.tolist()) <= grid and set(e.tolist()) <= grid, f"{track} bins off the grid"
    assert s.min() >= 0 and e.max() <= stress.length


@pytest.mark.parametrize(("track", "strand"), [("gene_plus", "+"), ("gene_minus", "-")])
def test_v1_gene_track_bins_hold_the_detected_share(
    v01: base.Instance, stress: StressRun, genes: pd.DataFrame, track: str, strand: str
) -> None:
    edges = _edges(stress)
    g = genes[genes["strand"] == strand]
    idx = _bin_of(edges, g["start"].astype(int).to_numpy())
    det = (g["detection_status"] == "detected").to_numpy()
    want = {int(edges[k]): float(det[idx == k].mean()) for k in np.unique(idx)}
    rows = v01.data[v01.data["track"] == track]
    got = dict(zip(rows["start"].astype(int), rows["value"].astype(float), strict=True))
    assert got.keys() == want.keys(), "non-empty bins differ from the genes' bins"
    np.testing.assert_allclose([got[k] for k in want], list(want.values()), atol=1e-12)
    assert set(rows["category"]) == {"detected_share"}


def test_v1_abundance_bins_hold_the_median_measured_log2(
    v01: base.Instance, stress: StressRun, genes: pd.DataFrame
) -> None:
    edges = _edges(stress)
    idx = _bin_of(edges, genes["start"].astype(int).to_numpy())
    y = pd.to_numeric(genes["abundance_log2"]).to_numpy(dtype=float)
    want: dict[int, float | None] = {}
    for k in np.unique(idx):
        v = y[(idx == k) & np.isfinite(y)]
        want[int(edges[k])] = float(np.median(v)) if len(v) else None
    rows = v01.data[v01.data["track"] == "abundance"]
    got = dict(zip(rows["start"].astype(int), rows["value"], strict=True))
    assert got.keys() == want.keys()
    for k, w in want.items():
        if w is None:
            assert _na(got[k]), (k, got[k])
        else:
            assert float(got[k]) == pytest.approx(w, abs=1e-12), k
    assert set(rows["category"]) == {"bin_median"}


def test_v1_origin_crossing_genes_land_in_the_bin_of_their_start(
    v01: base.Instance, stress: StressRun, genes: pd.DataFrame
) -> None:
    edges = _edges(stress)
    track = {"+": "gene_plus", "-": "gene_minus"}
    for g in _crossing(stress, genes).itertuples(index=False):
        k = int(_bin_of(edges, np.array([int(g.start)]))[0])
        assert k == len(edges) - 2, (g.locus_tag, k)  # the last bin, never wrapped to bin 0
        rows = v01.data[v01.data["track"] == track[str(g.strand)]]
        hit = rows[rows["start"].astype(int) == int(edges[k])]
        assert len(hit) == 1 and int(hit["end"].iloc[0]) == stress.length, (g.locus_tag, hit)
        assert int(rows["end"].astype(int).max()) <= stress.length


@pytest.mark.parametrize("track", GENE_TRACKS)
def test_v1_each_track_draws_at_most_20000_elements(v01: base.Instance, track: str) -> None:
    n = len(elements(v01.figure, track))
    assert 0 < n <= MAX_ELEMENTS, f"{track}: {n:,} drawn elements"


def test_v1_passes_the_checklist_element_item_at_20_mb(v01: base.Instance, tmp_path: Path) -> None:
    """M9-26 work item 1 / M9-24: every track ≤ 20,000 elements, as the checklist gates it."""
    from gippa.viz import checklist

    files = write_instance(v01, "V01", "circos", tmp_path / "plots", ["svg", "png"], 100)
    (res,) = [r for r in checklist.check_instance(v01, "V01", files) if r.item == ELEMENT_ITEM]
    assert res.passed is True, res.detail


# ---------------------------------------------------------------- V2


def test_example_1_v2_points_binned_to_bin_medians(
    v02: base.Instance, stress: StressRun, genes: pd.DataFrame
) -> None:
    edges = _edges(stress)
    g = genes[genes["contig"].astype(str) == stress.contig]
    y = pd.to_numeric(g["abundance_log2"]).to_numpy(dtype=float)
    keep = (g["detection_status"] == "detected").to_numpy() & np.isfinite(y)
    s = pd.to_numeric(g["start"]).to_numpy(dtype=float)
    e = pd.to_numeric(g["end"]).to_numpy(dtype=float)
    mids = np.mod((s + e) / 2.0, float(stress.length))[keep]
    idx = _bin_of(edges, mids)
    want = {
        int((edges[k] + edges[k + 1]) // 2): float(np.median(y[keep][idx == k]))
        for k in np.unique(idx)
    }
    pts = v02.data[v02.data["kind"] == "point"]
    assert len(pts) <= MAX_ELEMENTS, f"{len(pts)} V2 points"
    assert len(pts) < stress.n_points, "V2 drew one point per gene: not binned"
    got = dict(zip(pts["position_bp"].astype(int), pts["value"].astype(float), strict=True))
    assert got.keys() == want.keys()
    np.testing.assert_allclose([got[k] for k in want], list(want.values()), atol=1e-12)


def test_v2_origin_crossing_genes_sit_at_their_midpoint_bin_mod_l(
    v02: base.Instance, stress: StressRun, genes: pd.DataFrame
) -> None:
    edges = _edges(stress)
    g = _crossing(stress, genes)
    s, e = g["start"].astype(float).to_numpy(), g["end"].astype(float).to_numpy()
    mids = np.mod((s + e) / 2.0, float(stress.length))
    pts = v02.data[v02.data["kind"] == "point"]
    centres = set(pts["position_bp"].astype(int).tolist())
    for k in _bin_of(edges, mids).tolist():
        assert int((edges[k] + edges[k + 1]) // 2) in centres, k
    assert int(pts["position_bp"].astype(float).max()) < stress.length
    assert int(pts["position_bp"].astype(float).min()) >= 0


def test_v2_rug_stays_one_row_per_undetected_gene_under_the_limit(
    v02: base.Instance, stress: StressRun
) -> None:
    rug = v02.data[v02.data["kind"] == "rug"]
    assert len(rug) == stress.n_genes - stress.n_points
    assert len(rug) <= MAX_ELEMENTS


def test_v2_passes_the_checklist_element_item_at_20_mb(v02: base.Instance, tmp_path: Path) -> None:
    from gippa.viz import checklist

    files = write_instance(
        v02, "V02", "abundance_landscape", tmp_path / "plots", ["svg", "png"], 100
    )
    (res,) = [r for r in checklist.check_instance(v02, "V02", files) if r.item == ELEMENT_ITEM]
    assert res.passed is True, res.detail
