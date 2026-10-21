"""M5-06 islands: spec examples, the degenerate-case table and the lead decisions.

docs/specs/M5-06-islands.md (GIAP_SPEC.md §7.3.5(d) steps 1-10, §9.3 ``tables/islands.tsv``),
ADR-0021 with Amendment 1 (block length l' = min(20, 2 l-hat + 1, T); ``uncalibrated`` dropped,
the run report carries l-hat, l' and the lag-1 Moran's I instead) and lead decisions 21-27 and
40-51. Inputs are hand-built §9.3 gene tables with the
C5 ``tus`` and C1 ``replicons`` frames (tests/helpers/islands_data.py); the expectations come from
the spec text and from the oracle there, which composes the merged engines in the documented draw
order.

Interface assumed (tests/helpers/islands_data.py lists it in full): ``TuSeries`` carries the
replicon ``length`` and the per-TU ``start`` / ``end`` that the unwrapped coordinates of step 7
need; ``scan_replicons`` returns an ``IslandScan`` with ``series``, ``scanned``, ``block_len``,
``windows``, ``p_value``, ``clusters``, ``null_max``, ``n_perm``, per-replicon ``flags`` and the
run-report fields ``l_hat`` / ``morans_i_lag1`` (decision 48); ``low_direction_ceiling(q, w)`` is
the module-level ceiling of decision 50;
``islands.py`` holds the engine functions as module attributes (``scan_null_max`` is
monkeypatched below).

Adaptations of the examples, with the reason:

* examples 4-5 are tested on hand-built 600-TU replicons, where a planted island is found at
  ``stats.alpha.context`` for every seed tried; on giapsim they need ``detection = "abundance"``
  and a much larger shift (test_islands_sim.py), because ~45 % of genes are censored independently
  of abundance there (ADR-0021 "Power").
* example 6 keeps its structure (two overlapping significant high windows merge into one island
  reported by its best window) but not the literal S values 5.1 / 5.6, which depend on the Blom
  scores of the frame; the expected S comes from the window table.
* example 9 plants a 40-TU plasmid, not the 30-TU one of the spec table: at the 10 % span a 30-TU
  replicon gives k = 3 < 4, which decision 41 excludes. The 30-TU case is the decision-41 row.
* content tests that need islands to exist (coordinates, gene lists, wrap, flags) set
  ``stats.alpha.context = 1``: every non-zero window is then significant, the merged clusters are
  wide, and the reported span must still be the best window (decision 42).
"""

from __future__ import annotations

import dataclasses
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import islands_data as idat
from islands_data import Contig


def _isl() -> Any:
    from giap.context import islands

    return islands


def _rng(seed: int = 11) -> Any:
    from giap.rng import for_step

    return for_step(idat.STEP, seed)


def _ix(genes: pd.DataFrame, contig: str) -> pd.DataFrame:
    from giap.context.regions import replicon_index

    return replicon_index(genes, contig)


def _run(
    genes: pd.DataFrame, tus: pd.DataFrame, replicons: pd.DataFrame, cfg: Any, seed: int = 11
) -> tuple[pd.DataFrame, Any]:
    """The whole step on one random stream: (islands frame, IslandScan)."""
    return idat.step(genes, idat.build_series(genes, tus, replicons, cfg), cfg, _rng(seed))


def _bump(
    name: str,
    n_genes: int,
    bumps: dict[tuple[int, int], float],
    *,
    circular: bool = True,
    undetected: range | None = None,
) -> Contig:
    """A replicon of single-gene TUs with log2 abundance ``10 + i / n`` (a ramp the LOESS removes)
    plus the given bumps over half-open gene ranges; every gene detected except ``undetected``."""
    a = [10.0 + i / n_genes for i in range(n_genes)]
    for (lo, hi), v in bumps.items():
        for i in range(lo, hi):
            a[i] += v
    d = [1] * n_genes
    for i in undetected or ():
        d[i] = 0
    return Contig(name=name, d=d, p=[0.6] * n_genes, a=a, circular=circular)


def _window_tus(scan: Any, r: int, i: int) -> set[int]:
    res, n = scan.windows[r], len(scan.series[r].tu_ids)
    return {(int(res.start[i]) + k) % n for k in range(int(res.length[i]))}


# --------------------------------------------------------------------------- interface


def test_tu_series_dataclass() -> None:
    M = _isl()
    assert dataclasses.is_dataclass(M.TuSeries)
    assert {f.name for f in dataclasses.fields(M.TuSeries)} == {
        "contig", "circular", "length", "tu_ids", "genes", "start", "end", "score", "detrended",
    }  # fmt: skip
    assert M.TuSeries.__dataclass_params__.frozen
    assert hasattr(M.TuSeries, "__slots__")


def test_island_scan_dataclass() -> None:
    M = _isl()
    assert dataclasses.is_dataclass(M.IslandScan)
    assert {f.name for f in dataclasses.fields(M.IslandScan)} == {
        "series", "scanned", "l_hat", "block_len", "morans_i_lag1", "windows", "p_value",
        "clusters", "null_max", "n_perm", "flags",
    }  # fmt: skip


def test_islands_frame_columns_and_constants() -> None:
    genes, tus, reps = idat.frames([_bump("chr1", 200, {(40, 48): 2.0})])
    cfg = idat.config(alpha=1.0)
    out, _ = _run(genes, tus, reps, cfg)
    assert set(out.columns) == set(idat.ISLAND_COLUMNS)
    assert len(out) > 0
    assert (out["method"] == idat.METHOD).all()
    assert (out["null_model"] == idat.NULL_MODEL).all()
    assert out["q_value"].map(idat.is_na).all()
    assert np.array_equal(out["stat"].to_numpy(float), out["scan_s"].to_numpy(float))
    assert list(out["island_id"]) == [f"ISL_{i + 1}" for i in range(len(out))]
    assert (out["n_perm"].to_numpy(int) == cfg.stats.bootstrap).all()
    idat.assert_matches_oracle(out, idat.oracle_islands(genes, tus, reps, cfg, 11)[0])


# --------------------------------------------------------------------------- examples 1-3


def test_example1_gene_scores_left_censored() -> None:
    """Spec example 1: detected log2 [5, -, 3, -] -> [1.0491, -0.6289, 0.2993, -0.6289]."""
    genes, _, _ = idat.frames([Contig("chr1", d=[1, 0, 1, 0], p=[0.5] * 4, a=[5, 0, 3, 0])])
    s = _isl().gene_scores(_ix(genes, "chr1"))
    assert list(s.index) == [f"chr1_g{i:05d}" for i in range(4)]
    assert np.allclose(s.to_numpy(float), [1.0491, -0.6289, 0.2993, -0.6289], atol=5e-5)


def test_example2_tu_scores_are_gene_score_means() -> None:
    """Spec example 2: TU_1 = {G1, G2}, TU_2 = {G3}, TU_3 = {G4} gives
    z = [0.2101, 0.2993, -0.6289]."""
    genes, tus, _ = idat.frames(
        [
            Contig(
                "chr1",
                d=[1, 0, 1, 0],
                p=[0.5] * 4,
                a=[5, 0, 3, 0],
                tu=["TU_1", "TU_1", "TU_2", "TU_3"],
            )
        ]
    )
    ts = _isl().tu_series(_ix(genes, "chr1"), tus, contig="chr1", circular=True, length=1800,
                          span_frac=0.10)  # fmt: skip
    assert ts.tu_ids == ("TU_1", "TU_2", "TU_3")
    assert ts.genes == (("chr1_g00000", "chr1_g00001"), ("chr1_g00002",), ("chr1_g00003",))
    assert np.allclose(ts.score, [0.2101, 0.2993, -0.6289], atol=5e-5)
    assert ts.contig == "chr1" and ts.circular and ts.length == 1800
    assert list(ts.start) == [100, 900, 1300] and list(ts.end) == [800, 1200, 1600]


def test_example3_detected_gene_without_abundance_is_not_scored() -> None:
    """Spec example 3 / decision 22: under ``missing = exclude`` a detected gene with NaN
    abundance is dropped; its TU keeps its other genes, and a TU with none left is dropped."""
    genes, tus, reps = idat.frames(
        [
            Contig(
                "chr1",
                d=[1, 1, 1, 0],
                p=[0.5] * 4,
                a=[5, math.nan, math.nan, 0],
                tu=["TU_1", "TU_1", "TU_2", "TU_3"],
            )
        ]
    )
    cfg = idat.config(abundance_missing="exclude")
    s = _isl().gene_scores(_ix(genes, "chr1"))
    assert list(s.index) == ["chr1_g00000", "chr1_g00003"]
    ts = _isl().tu_series(_ix(genes, "chr1"), tus, contig="chr1", circular=True, length=1800,
                          span_frac=cfg.context.islands.detrend_span_frac)  # fmt: skip
    assert ts.tu_ids == ("TU_1", "TU_3")  # TU_2 kept only the unscored gene
    assert ts.genes == (("chr1_g00000",), ("chr1_g00003",))
    assert np.allclose(ts.score, idat.oracle_series(genes, tus, reps, cfg)[0].score)


# --------------------------------------------------------------------------- examples 4-5


def test_example4_planted_high_island_is_found() -> None:
    """Spec example 4: a planted high island gives one ``high`` island overlapping the planted
    TUs, ``p_value <= alpha``, ``method`` ``scan_max_stat`` and ``q_value`` NA."""
    genes, tus, reps = idat.frames([_bump("chr1", 600, {(100, 120): 10.0})])
    cfg = idat.config(alpha=0.05)
    out, _ = _run(genes, tus, reps, cfg)
    assert len(out) == 1
    r = out.iloc[0]
    assert r["direction"] == "high" and float(r["scan_s"]) > 0
    assert float(r["p_value"]) <= cfg.stats.alpha.context
    assert r["method"] == idat.METHOD and idat.is_na(r["q_value"])
    planted = {f"TU_chr1_{i}" for i in range(100, 120)}
    span = set(r["genes"]) & {f"chr1_g{i:05d}" for i in range(100, 120)}
    assert span, "the island must overlap the planted TUs"
    assert r["tu_first"] in planted and r["tu_last"] in planted
    idat.assert_matches_oracle(out, idat.oracle_islands(genes, tus, reps, cfg, 11)[0])


def test_example5_planted_low_island_is_found() -> None:
    """Spec example 5: the mirror case gives one ``low`` island with ``scan_s < 0`` and a
    ``detected_frac`` below the replicon's."""
    genes, tus, reps = idat.frames([_bump("chr1", 600, {}, undetected=range(500, 520))])
    cfg = idat.config(alpha=0.05)
    out, _ = _run(genes, tus, reps, cfg)
    assert len(out) == 1
    r = out.iloc[0]
    assert r["direction"] == "low" and float(r["scan_s"]) < 0
    assert float(r["p_value"]) <= cfg.stats.alpha.context
    assert r["tu_first"] == "TU_chr1_500" and r["tu_last"] == "TU_chr1_519"
    replicon_frac = float(genes["detected"].mean())
    assert float(r["detected_frac"]) == 0.0 < replicon_frac
    assert idat.is_na(r["mean_abundance_log2"])  # degenerate row 5: no detected gene


# --------------------------------------------------------------------------- examples 6-11


def test_example6_overlapping_windows_merge_into_one_island() -> None:
    """Spec example 6: overlapping significant high windows give one island, reported by its best
    window; its S is the largest |S| of the cluster."""
    genes, tus, reps = idat.frames([_bump("chr1", 200, {(10, 13): 1.0, (13, 21): 3.0})])
    cfg = idat.config(alpha=1.0)
    out, scan = _run(genes, tus, reps, cfg)
    high = out[out["direction"] == "high"].reset_index(drop=True)
    assert len(high) == 1
    r = high.iloc[0]
    cluster = next(c for c in scan.clusters[0] if c.direction == "high")
    assert len(cluster.members) >= 2, "the fixture is meant to merge several windows"
    res = scan.windows[0]
    best = cluster.best
    assert math.isclose(float(r["scan_s"]), float(res.s[best]), abs_tol=1e-12)
    assert abs(res.s[best]) == max(abs(res.s[i]) for i in cluster.members)
    tu_ids = scan.series[0].tu_ids
    assert int(r["n_tus"]) == int(res.length[best]) <= cfg.context.islands.max_tus
    assert r["tu_first"] == tu_ids[int(res.start[best])]
    planted = {f"TU_chr1_{i}" for i in range(10, 21)}
    assert {tu_ids[i] for i in _window_tus(scan, 0, best)} <= planted
    idat.assert_matches_oracle(out, idat.oracle_islands(genes, tus, reps, cfg, 11)[0])


def test_decision42_span_is_the_best_window_never_the_cluster_union() -> None:
    """Decision 42: the island's span is the best window's (start, length); the union of its
    cluster's members spills past it and must not be reported."""
    genes, tus, reps = idat.frames([_bump("chr1", 200, {(60, 80): 2.5})])
    cfg = idat.config(alpha=1.0)
    out, scan = _run(genes, tus, reps, cfg)
    assert len(out) >= 1
    wider = 0
    for row, cluster in zip(out.itertuples(), scan.clusters[0], strict=True):
        res, tu_ids = scan.windows[0], scan.series[0].tu_ids
        span = _window_tus(scan, 0, cluster.best)
        union = set().union(*(_window_tus(scan, 0, i) for i in cluster.members))
        assert span <= union
        wider += len(union) > len(span)
        assert int(row.n_tus) == len(span) <= cfg.context.islands.max_tus
        assert row.tu_first == tu_ids[int(res.start[cluster.best])]
        last = (int(res.start[cluster.best]) + int(res.length[cluster.best]) - 1) % len(tu_ids)
        assert row.tu_last == tu_ids[last]
        assert set(row.genes) == {f"chr1_g{i:05d}" for i in span}
    assert wider, "the fixture is meant to have at least one cluster wider than its best window"


def test_example7_wrapping_window_reports_unwrapped_coordinates() -> None:
    """Spec example 7: on a circular replicon the best window may wrap; ``tu_first`` is the TU at
    the window start (near the end of the replicon) and ``end`` passes the replicon length."""
    n = 200
    genes, tus, reps = idat.frames([_bump("chr1", n, {(0, 4): 3.0, (n - 4, n): 3.0})])
    cfg = idat.config(alpha=1.0)
    out, _ = _run(genes, tus, reps, cfg)
    length = int(reps.loc[reps["contig"] == "chr1", "length"].iloc[0])
    high = out[out["direction"] == "high"].reset_index(drop=True)
    assert len(high) == 1
    r = high.iloc[0]
    first, last = int(r["tu_first"].rsplit("_", 1)[1]), int(r["tu_last"].rsplit("_", 1)[1])
    assert first > last, "the best window is meant to wrap the origin"
    assert int(r["start"]) == 400 * first + 100
    assert int(r["end"]) == 400 * last + 400 + length > length
    assert int(r["n_tus"]) == (n - first) + last + 1
    idat.assert_matches_oracle(out, idat.oracle_islands(genes, tus, reps, cfg, 11)[0])


@pytest.mark.parametrize(
    ("bump", "first", "last", "start", "end"),
    [
        # T196..T199: ends at the crossing TU, whose unwrapped end is 400 + L (no extra wrap).
        ((196, 201), "T196", "T199", 400 * 196 + 100, 400 + 80_200),
        # T199, T1, T2: starts at the crossing TU near L and wraps; end is T2's end + L.
        ((199, 203), "T199", "T2", 400 * 199 + 100, 400 * 2 + 400 + 80_200),
    ],
    ids=["window-ends-at-crossing-tu", "window-starts-at-crossing-tu"],
)
def test_origin_crossing_tu_sorts_last_and_gives_valid_coordinates(
    bump: tuple[int, int], first: str, last: str, start: int, end: int
) -> None:
    """Regression (post-merge review of #131): a C5 TU holding the last and first genes of a
    circular replicon (start near the length, end past it) is last in TU order (C5 start), so a
    best window ending or starting at it has 0 <= start < length and start < end <= start + L."""
    n = 200
    a = [10.0 + 0.01 * ((i * 37) % 11) for i in range(n)]  # flat, no ramp to favour a neighbour
    for i in range(*bump):  # gene indices mod n: the bump covers both genes of the crossing TU
        a[i % n] += 4.0
    tu = [f"T{i}" for i in range(n)]
    tu[0] = f"T{n - 1}"  # genes n-1 and 0 form one TU across the origin
    c = Contig(name="chr1", d=[1] * n, p=[0.6] * n, a=a, circular=True, tu=tu)
    genes, tus, reps = idat.frames([c])
    tus = idat.cross_origin(genes, tus, reps, "chr1", f"T{n - 1}")
    length = int(reps["length"].iloc[0])
    cfg = idat.config(alpha=1.0)
    out, scan = _run(genes, tus, reps, cfg)
    ts = scan.series[0]
    assert ts.tu_ids[-1] == f"T{n - 1}" and ts.tu_ids[0] == "T1"
    assert int(ts.start[-1]) == 400 * (n - 1) + 100
    assert int(ts.end[-1]) == 400 + length
    for _, r in out.iterrows():
        assert 0 <= int(r["start"]) < length
        assert int(r["start"]) < int(r["end"]) <= int(r["start"]) + length
    high = out[out["direction"] == "high"].reset_index(drop=True)
    assert len(high) == 1
    r = high.iloc[0]
    assert (r["tu_first"], r["tu_last"]) == (first, last)
    assert length == 80_200
    assert (int(r["start"]), int(r["end"])) == (start, end)
    idat.assert_matches_oracle(out, idat.oracle_islands(genes, tus, reps, cfg, 11)[0])


def test_example8_integer_block_len_is_used_and_draws_no_permutations() -> None:
    """Spec example 8: ``context.islands.block_len = 5`` is reported in every row and no
    block-length permutation is drawn (the frame equals the oracle, which skips them)."""
    genes, tus, reps = idat.frames([_bump("chr1", 200, {(40, 48): 2.5})])
    cfg = idat.config(alpha=1.0, block_len=5)
    out, scan = _run(genes, tus, reps, cfg)
    assert scan.block_len == (5,)
    assert (out["block_len"].to_numpy(int) == 5).all()
    idat.assert_matches_oracle(out, idat.oracle_islands(genes, tus, reps, cfg, 11)[0])


def test_example9_plasmid_and_chromosome_share_one_null_maximum() -> None:
    """Spec example 9 (40 TUs, see the module docstring): the plasmid is scanned and enters the
    genome-wide maximum, so M_b is the row-wise maximum over both replicons."""
    genes, tus, reps = idat.frames(
        [_bump("chr1", 200, {(40, 48): 2.0}), _bump("p1", 40, {(5, 9): 2.0})]
    )
    cfg = idat.config(alpha=1.0)
    out, scan = _run(genes, tus, reps, cfg)
    want_frame, want = idat.oracle_islands(genes, tus, reps, cfg, 11)
    assert scan.scanned == (True, True)
    assert set(out["contig"]) == {"chr1", "p1"}
    assert np.allclose(scan.null_max, want.null_max)
    assert np.allclose(scan.null_max, np.maximum(*want.replicon_null))
    assert int(scan.n_perm) == cfg.stats.bootstrap
    idat.assert_matches_oracle(out, want_frame)


def test_example10_unscanned_plasmid_changes_nothing_bitwise() -> None:
    """Spec example 10: a plasmid with 3 scored TUs is not scanned and draws nothing, so the
    chromosome's islands are bitwise identical with and without it."""
    chrom = _bump("chr1", 200, {(40, 48): 2.0})
    genes1, tus1, reps1 = idat.frames([chrom])
    genes2, tus2, reps2 = idat.frames([chrom, _bump("p1", 3, {})])
    cfg = idat.config(alpha=1.0)
    a, scan_a = _run(genes1, tus1, reps1, cfg)
    b, scan_b = _run(genes2, tus2, reps2, cfg)
    assert scan_b.scanned == (True, False)
    assert idat.FEW_VALUES in set(scan_b.flags[1])
    assert np.array_equal(scan_a.null_max, scan_b.null_max)
    pd.testing.assert_frame_equal(a, b)


def test_example11_same_seed_same_output() -> None:
    """Spec example 11: identical islands frame and gene map for the same seed."""
    genes, tus, reps = idat.frames([_bump("chr1", 200, {(40, 48): 2.0})])
    cfg = idat.config(alpha=1.0)
    a, _ = _run(genes, tus, reps, cfg, seed=5)
    b, _ = _run(genes, tus, reps, cfg, seed=5)
    pd.testing.assert_frame_equal(a, b)
    pd.testing.assert_frame_equal(_isl().gene_island_map(a), _isl().gene_island_map(b))


# --------------------------------------------------------------------------- the family


def test_family_is_the_genome_wide_maximum_over_scanned_replicons() -> None:
    """Decisions 21 and 42: one family, not one per replicon. M_b is the row-wise maximum over
    every scanned replicon, so each replicon's p-values are at least the p-values its own maximum
    would give; with k scanned replicons (read from ``IslandScan.scanned``) a per-replicon family
    would leave the run-level rate at 1 - (1 - alpha)^k instead."""
    from giap.stats.scan import fwer_p

    genes, tus, reps = idat.frames(
        [_bump("chr1", 200, {(40, 48): 2.0}), _bump("p1", 60, {(5, 9): 2.0})]
    )
    cfg = idat.config(alpha=1.0)
    _, scan = _run(genes, tus, reps, cfg)
    _, want = idat.oracle_islands(genes, tus, reps, cfg, 11)
    assert sum(scan.scanned) == 2
    assert np.allclose(scan.null_max, np.maximum(*want.replicon_null))
    for r in range(2):
        assert np.allclose(scan.p_value[r], fwer_p(np.abs(scan.windows[r].s), scan.null_max))
        own = fwer_p(np.abs(scan.windows[r].s), want.replicon_null[r])
        assert (scan.p_value[r] >= own - 1e-12).all()
        assert (scan.p_value[r] > own + 1e-12).any(), "the family must bind on some window"


def test_negative_s_is_passed_as_abs_to_fwer_p() -> None:
    """Decision 42: ``fwer_p`` raises on a negative |S|, so the step passes ``np.abs(res.s)``;
    low windows get the p of their mirror."""
    from giap.stats.scan import fwer_p

    genes, tus, reps = idat.frames([_bump("chr1", 200, {(40, 48): -3.0})])
    cfg = idat.config(alpha=1.0)
    out, scan = _run(genes, tus, reps, cfg)
    assert (scan.windows[0].s < 0).any()
    assert np.allclose(scan.p_value[0], fwer_p(np.abs(scan.windows[0].s), scan.null_max))
    low = out[out["direction"] == "low"]
    assert len(low) >= 1
    assert (low["scan_s"].to_numpy(float) < 0).all()
    p = low["p_value"].to_numpy(float)
    assert ((p > 0) & (p <= 1)).all()


# --------------------------------------------------------------------------- degenerate table


@pytest.mark.parametrize("n_tus", [1, 3, 30])
def test_replicon_with_too_small_a_span_is_not_scanned(n_tus: int) -> None:
    """Degenerate row 1 with decisions 40/41: T < 2 * min_tus, or a span giving k < 4
    (T <= 30 at the 10 % span), is not scanned: no island, the run flag ``few_values``, no
    contribution to M_b."""
    genes, tus, reps = idat.frames([_bump("p1", n_tus, {})])
    out, scan = _run(genes, tus, reps, idat.config(alpha=1.0))
    assert scan.scanned == (False,)
    assert scan.flags[0] == (idat.FEW_VALUES,)
    assert len(out) == 0
    assert scan.n_perm == 0


def test_span_of_31_tus_is_scanned() -> None:
    """Decision 41 boundary: 31 TUs give k = 4 at the 10 % span, so the replicon is scanned."""
    genes, tus, reps = idat.frames([_bump("p1", 31, {})])
    _, scan = _run(genes, tus, reps, idat.config(alpha=1.0))
    assert scan.scanned == (True,)
    assert idat.FEW_VALUES not in set(scan.flags[0])


def test_no_replicon_scanned_gives_an_empty_table() -> None:
    """Degenerate row 2: ``islands.tsv`` header only, and no gene gets an island."""
    genes, tus, reps = idat.frames([_bump("p1", 3, {}), _bump("p2", 2, {})])
    out, scan = _run(genes, tus, reps, idat.config(alpha=1.0))
    assert len(out) == 0
    assert set(out.columns) == set(idat.ISLAND_COLUMNS)
    assert scan.scanned == (False, False)
    gmap = _isl().gene_island_map(out)
    assert len(gmap) == 0
    assert list(gmap.columns) == list(idat.GENE_MAP_COLUMNS)


def test_all_genes_undetected_gives_zero_scores_and_no_island() -> None:
    """Degenerate row 3: every Blom score is 0, so S = 0, p = 1 and no window has a direction."""
    n = 60
    genes, tus, reps = idat.frames([Contig("chr1", d=[0] * n, p=[0.3] * n, a=[math.nan] * n)])
    out, scan = _run(genes, tus, reps, idat.config(alpha=1.0))
    assert scan.scanned == (True,)
    assert np.allclose(scan.series[0].score, 0.0)
    assert np.allclose(scan.windows[0].s, 0.0)
    assert np.allclose(scan.p_value[0], 1.0)
    assert len(out) == 0


def test_nan_null_maximum_gives_no_island(monkeypatch: pytest.MonkeyPatch) -> None:
    """Degenerate row 4 with decision 42: with no valid M_b the p-values are NaN, and ``p <=
    alpha`` keeps every window out (a NaN p never counts as significant)."""
    M = _isl()
    genes, tus, reps = idat.frames([_bump("chr1", 200, {(40, 48): 2.0})])
    monkeypatch.setattr(
        M, "scan_null_max", lambda Z, *a, **k: np.full(np.shape(Z)[0], np.nan), raising=True
    )
    out, scan = _run(genes, tus, reps, idat.config(alpha=1.0))
    assert np.isnan(scan.p_value[0]).all()
    assert scan.n_perm == 0
    assert len(out) == 0


def test_best_window_of_two_tus_is_flagged_ci_few_units() -> None:
    """Degenerate row 6: an effect CI from fewer than 10 TUs carries the M4-16 ``ci_few_units``."""
    genes, tus, reps = idat.frames([_bump("chr1", 200, {(40, 42): 6.0})])
    out, _ = _run(genes, tus, reps, idat.config(alpha=1.0))
    row = out[out["tu_first"] == "TU_chr1_40"].iloc[0]
    assert int(row["n_tus"]) == 2
    assert idat.FEW_UNITS in idat.flag_set(row["flags"])
    assert float(row["ci_low"]) <= float(row["ci_high"])
    assert all(
        (idat.FEW_UNITS in idat.flag_set(f)) == (n < 10)
        for f, n in zip(out["flags"], out["n_tus"].to_numpy(int), strict=True)
    )


def test_block_len_above_t_is_bounded_to_t() -> None:
    """Degenerate row 7: an integer ``block_len`` larger than T is bounded to T."""
    genes, tus, reps = idat.frames([_bump("p1", 40, {(5, 9): 2.0})])
    out, scan = _run(genes, tus, reps, idat.config(alpha=1.0, block_len=200))
    assert len(scan.series[0].tu_ids) == 40
    assert scan.block_len == (40,)
    assert (out["block_len"].to_numpy(int) == 40).all()


def test_max_tus_above_t_is_capped_by_the_engine() -> None:
    """Degenerate row 8: on a circular replicon the windows stop at T - 1 TUs (M5-03)."""
    genes, tus, reps = idat.frames([_bump("p1", 40, {(5, 9): 2.0})])
    out, scan = _run(genes, tus, reps, idat.config(alpha=1.0, max_tus=100))
    t = len(scan.series[0].tu_ids)
    assert int(scan.windows[0].length.max()) == t - 1
    assert int(out["n_tus"].max()) <= t - 1


def test_block_length_follows_adr0021() -> None:
    """ADR-0021 decision 1: l' = min(20, 2 l-hat + 1, T), l-hat from ``choose_block_len`` on the
    detrended scores with the gene-label null at the literal 0.05 (decision 24)."""
    from giap.stats.autocorr import block_len_prime, choose_block_len

    genes, tus, reps = idat.frames([_bump("chr1", 200, {(40, 48): 2.0})])
    cfg = idat.config(alpha=1.0)
    M = _isl()
    series = idat.oracle_series(genes, tus, reps, cfg)
    ts = M.tu_series(_ix(genes, "chr1"), tus, contig="chr1", circular=True,
                     length=int(reps["length"].iloc[0]), span_frac=0.10)  # fmt: skip
    ell = M.block_length(ts, cfg, _rng(3))
    l_hat = choose_block_len(series[0].detrended, True, cfg.stats.permutations, _rng(3), 0.05, 20)
    assert ell == block_len_prime(l_hat, series[0].detrended.size) <= 20


# ------------------------------------------- ADR-0021 Amendment 1: run report, no uncalibrated


def test_no_replicon_or_island_is_flagged_uncalibrated() -> None:
    """ADR-0021 Amendment 1 / decision 48: ``uncalibrated`` is dropped; the only island flag on a
    small replicon is ``ci_few_units``."""
    genes, tus, reps = idat.frames([_bump("chr1", 200, {(40, 48): 2.0})])
    out, scan = _run(genes, tus, reps, idat.config(alpha=1.0))
    assert len(scan.series[0].tu_ids) < 1000 and len(out) >= 1
    assert all(idat.flag_set(f) <= {idat.FEW_UNITS} for f in out["flags"])
    assert all("uncalibrated" not in set(f) for f in scan.flags)


def test_run_report_carries_l_hat_block_len_and_lag1_moran() -> None:
    """Decision 48: the §9.4 run report needs l-hat, l' (``block_len``) and the lag-1 Moran's I of
    the detrended TU scores per replicon, so ``IslandScan`` exposes all three (0 / NaN for a
    replicon that is not scanned)."""
    from giap.stats.autocorr import block_len_prime, morans_i

    genes, tus, reps = idat.frames([_bump("chr1", 200, {(40, 48): 2.0}), _bump("p1", 3, {})])
    cfg = idat.config(alpha=1.0)
    _, scan = _run(genes, tus, reps, cfg)
    _, want = idat.oracle_islands(genes, tus, reps, cfg, 11)
    assert all(hasattr(scan, f) for f in idat.RUN_REPORT_FIELDS)
    assert scan.l_hat == want.l_hat and scan.block_len == want.block_len
    assert scan.block_len[0] == block_len_prime(scan.l_hat[0], len(scan.series[0].tu_ids))
    assert math.isclose(
        scan.morans_i_lag1[0], morans_i(scan.series[0].detrended, 1, True), abs_tol=1e-12
    )
    assert scan.l_hat[1] == scan.block_len[1] == 0  # the plasmid is not scanned
    assert math.isnan(scan.morans_i_lag1[1])


def test_low_direction_ceiling_formula_and_silenced_tus() -> None:
    """Decision 50: left censoring caps a low window of w TUs at sqrt(w) * |Phi^-1(q / 2)|, so
    silencing whole TUs cannot produce a low island once q is large. The ceiling is a public
    function of the step, because the run report prints it against the null maximum."""
    n, q = 200, 0.45
    ceiling = _isl().low_direction_ceiling(q, 20)
    assert math.isclose(ceiling, idat.ceiling_from_q(q, 20), rel_tol=0, abs_tol=1e-12)
    assert 3.3 < ceiling < 3.5  # the statistician's 3.4 at q = 0.45, w = 20
    rng = np.random.default_rng(50)
    d = [int(v) for v in rng.random(n) >= q]
    for i in range(100, 120):  # silence 20 consecutive TUs entirely
        d[i] = 0
    a = [10.0 + i / n for i in range(n)]
    genes, tus, reps = idat.frames([Contig("chr1", d=d, p=[0.55] * n, a=a)])
    cfg = idat.config(alpha=0.05)
    out, scan = _run(genes, tus, reps, cfg)
    w = idat.max_window(scan.series[0], cfg)
    got = _isl().low_direction_ceiling(idat.censored_share(genes, "chr1"), w)
    assert float(np.min(scan.windows[0].s)) >= -got
    assert "low" not in set(out["direction"])


def test_module_docstring_carries_the_low_direction_caveat() -> None:
    """Decision 50: M5-06 does not claim a two-sided procedure without the ceiling caveat."""
    doc = (_isl().__doc__ or "").lower()
    assert "ceiling" in doc
    assert "censor" in doc and "low" in doc


# --------------------------------------------------------------------------- gene map


def test_gene_island_map_assigns_every_island_gene() -> None:
    """Step 9: ``locus_tag`` -> ``island_id`` for every gene of every island."""
    genes, tus, reps = idat.frames([_bump("chr1", 200, {(40, 48): 2.0, (120, 128): -2.0})])
    out, _ = _run(genes, tus, reps, idat.config(alpha=1.0))
    gmap = _isl().gene_island_map(out)
    assert list(gmap.columns) == list(idat.GENE_MAP_COLUMNS)
    assert gmap["locus_tag"].is_unique
    assert set(gmap["locus_tag"]) == {g for row in out["genes"] for g in row}
    for _, row in out.iterrows():
        for g in row["genes"]:
            assert gmap.loc[gmap["locus_tag"] == g, "island_id"].iloc[0] in set(out["island_id"])


def _island_row(**kw: Any) -> dict[str, Any]:
    row: dict[str, Any] = dict.fromkeys(idat.ISLAND_COLUMNS)
    row.update(
        contig="chr1", direction="high", method=idat.METHOD, null_model=idat.NULL_MODEL,
        q_value=math.nan, flags=(),
    )  # fmt: skip
    row.update(kw)
    return row


def test_gene_in_two_islands_takes_smaller_p_then_larger_abs_s_then_smaller_id() -> None:
    """Step 9 tie-break on a hand-built islands frame (a high and a low island sharing genes)."""
    gmap = _isl().gene_island_map(
        pd.DataFrame(
            [
                _island_row(
                    island_id="ISL_1",
                    direction="high",
                    p_value=0.04,
                    scan_s=3.0,
                    genes=["g1", "g2"],
                ),
                _island_row(
                    island_id="ISL_2",
                    direction="low",
                    p_value=0.01,
                    scan_s=-2.0,
                    genes=["g2", "g3"],
                ),
                _island_row(
                    island_id="ISL_3",
                    direction="high",
                    p_value=0.01,
                    scan_s=5.0,
                    genes=["g3", "g4"],
                ),
                _island_row(
                    island_id="ISL_4",
                    direction="low",
                    p_value=0.01,
                    scan_s=-5.0,
                    genes=["g4", "g5"],
                ),
            ]
        )
    )
    got = dict(zip(gmap["locus_tag"], gmap["island_id"], strict=True))
    assert got == {
        "g1": "ISL_1",  # only island
        "g2": "ISL_2",  # smaller p wins
        "g3": "ISL_3",  # equal p, larger |S| wins
        "g4": "ISL_3",  # equal p and |S|, smaller ISL number wins
        "g5": "ISL_4",
    }
