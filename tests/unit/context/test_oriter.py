"""M5-08 ori-ter rows: spec examples 1-6 and 8 and the degenerate table.

docs/specs/M5-08-oriter.md with lead decisions 31-35, GIAP_SPEC.md §7.3.5(e) and §9.3, ADR-0021
decision 1 (the CI block length l') **with Amendment 1** (the validity domain is a dependence
regime; the `uncalibrated` flag is dropped -- no row here may carry it).

Decisions that supersede the spec text and are tested here:

* 31: several oriC calls are averaged (circular mean) only when every pair lies within 0.05 L of
  each other; farther apart the replicon is ``no_oric``. Both branches log
  ``context_stats.oriter.multiple_oric``.
* 32: the arm is the free ``stat`` column (``oriter_slope_right`` / ``_left`` / ``_pooled``),
  the genes.tsv column ``replichore`` is ``right`` / ``left``.
* 33: linear replicons and ``unknown``-kind contigs are skipped (no row, NA gene columns), as are
  plasmids.
* 34: the regression uses detected genes' abundance only.
* 35 / ADR-0021: l' = min(20, 2 l-hat + 1, T) unless ``context.islands.block_len`` is an integer.

Examples 3-5 read the committed Bakta 1.12.2 JSONs through ``parse_bakta_json`` (no tool run).
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import oriter_data as od
from oriter_data import flag_set, is_na

K12_LENGTH = 4_641_652
K12_MIDS = (3_884_032.0, 3_925_819.0)
NISSLE_ORIC = 5_010_344.0


def _oriter() -> Any:
    from giap.context import oriter

    return oriter


# --------------------------------------------------------------------------- example 1


def test_example1_replichore_of_four_midpoints() -> None:
    """L = 1,000, oriC midpoint 900 -> (right, 0.4), (right, 1.0), (left, 0.8), (left, 0.5)."""
    arm, pos = _oriter().replichore([100, 400, 500, 650], 900.0, 1000)
    assert [str(a) for a in arm] == ["right", "right", "left", "left"]
    assert list(np.asarray(pos, dtype=float)) == pytest.approx([0.4, 1.0, 0.8, 0.5])


def test_example1_through_the_gene_columns() -> None:
    """The same four genes as a gene table: ``gene_replichores`` gives the same labels."""
    genes, _, replicons = od.placed([od.Placed("chr1", 1000, [100, 400, 500, 650], width=50)])
    out = od.gene_columns(genes, replicons, od.oric_at("chr1", 900))
    assert list(out["locus_tag"]) == list(genes["locus_tag"])
    assert list(out["replichore"]) == ["right", "right", "left", "left"]
    assert list(out["replichore_position"]) == pytest.approx([0.4, 1.0, 0.8, 0.5])


# --------------------------------------------------------------------------- example 2


def test_example2_origin_crossing_gene_and_unwrapped_oric() -> None:
    """Gene start 990, end 1,010 (unwrapped) with oriC at 0: midpoint 0 -> (right, 0.0)."""
    genes, _, replicons = od.placed([od.Placed("chr1", 1000, [1000], width=20)])
    assert int(genes["start"].iloc[0]) == 990
    assert int(genes["end"].iloc[0]) == 1010
    feats = od.oric("chr1", (990, 1010))  # the oriC may be unwrapped too
    assert _oriter().oric_position(feats, "chr1", 1000) == pytest.approx(0.0)
    out = od.gene_columns(genes, replicons, feats)
    assert str(out["replichore"].iloc[0]) == "right"
    assert float(out["replichore_position"].iloc[0]) == pytest.approx(0.0)


# --------------------------------------------------------------------------- example 3


def test_example3_k12_two_oric_calls_give_their_circular_mean() -> None:
    """Recorded K-12: oriC midpoints 3,884,032 and 3,925,819 -> 3,904,925.5 (+- 1 bp)."""
    feats, reps = od.recorded("k12")
    oric_rows = feats[feats["type"].astype(str) == "oric"]
    assert len(oric_rows) == 2
    assert int(reps.loc[reps["contig"] == "U00096.3", "length"].iloc[0]) == K12_LENGTH
    got = _oriter().oric_position(feats, "U00096.3", K12_LENGTH)
    assert got is not None
    assert got == pytest.approx(sum(K12_MIDS) / 2.0, abs=1.0)


def test_example3_k12_logs_multiple_oric(caplog: pytest.LogCaptureFixture) -> None:
    """The two calls are listed in a ``context_stats.oriter.multiple_oric`` log event (Q1)."""
    feats, _ = od.recorded("k12")
    with caplog.at_level(logging.INFO, logger="giap"):
        _oriter().oric_position(feats, "U00096.3", K12_LENGTH)
    events = [r for r in caplog.records if getattr(r, "event", None) == od.MULTIPLE_ORIC_EVENT]
    assert events, [getattr(r, "event", r.getMessage()) for r in caplog.records]
    data = getattr(events[0], "data", {})
    listed = [float(v) for v in (data.get("midpoints") or data.get("oric") or [])]
    assert sorted(listed) == pytest.approx(list(K12_MIDS), abs=1.0)


def test_decision31_calls_farther_than_five_percent_apart_are_no_oric(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Lead decision 31: beyond 0.05 L apart the circular mean is refused (``no_oric``)."""
    feats = od.oric_at("chr1", 1_000, 300_000)  # 0.0643 of L = 4,641,652
    with caplog.at_level(logging.INFO, logger="giap"):
        assert _oriter().oric_position(feats, "chr1", K12_LENGTH) is None
    assert any(getattr(r, "event", None) == od.MULTIPLE_ORIC_EVENT for r in caplog.records)
    # just inside the window: averaged
    close = od.oric_at("chr1", 1_000, 200_000)  # 0.0429 of L
    assert _oriter().oric_position(close, "chr1", K12_LENGTH) == pytest.approx(100_500.0, abs=1.0)


def test_circular_mean_wraps_across_the_origin() -> None:
    """Two calls either side of coordinate 0 average to a point near 0, not to L/2."""
    import oric_oracle as oo

    length = 1_000_000
    got = _oriter().oric_position(od.oric_at("chr1", 999_900, 100, width=40), "chr1", length)
    assert got is not None
    assert oo.circ_dist(float(got), 0.0, float(length)) == pytest.approx(0.0, abs=1.0)


# --------------------------------------------------------------------------- example 4


def test_example4_nissle_chromosome_oric_and_plasmids_ignored() -> None:
    """Nissle: oriC midpoint 5,010,344 on CP171242.1; oriV / oriT on plasmids are not oriC."""
    feats, reps = od.recorded("nissle")
    lengths = dict(zip(reps["contig"].astype(str), reps["length"], strict=True))
    assert set(feats.loc[feats["type"].isin(["oriv", "orit"]), "contig"]) == {
        "CP171244.1",
        "CP171245.1",
    }
    got = _oriter().oric_position(feats, "CP171242.1", int(lengths["CP171242.1"]))
    assert got == pytest.approx(NISSLE_ORIC, abs=1.0)
    for contig in ("CP171243.1", "CP171244.1", "CP171245.1"):
        assert _oriter().oric_position(feats, contig, int(lengths[contig])) is None


def test_example4_plasmids_get_no_rows_and_na_gene_columns() -> None:
    """The three Nissle plasmids contribute no row and NA gene columns (Q3 / decision 33)."""
    feats, reps = od.recorded("nissle")
    lengths = dict(zip(reps["contig"].astype(str), reps["length"], strict=True))
    specs = [
        od.Placed(
            contig,
            int(lengths[contig]),
            od.right_arm_mids(
                NISSLE_ORIC if contig == "CP171242.1" else 0.0,
                int(lengths[contig]),
                [0.1, 0.3, 0.5, 0.7, 0.9],
            ),
            kind="chromosome" if contig == "CP171242.1" else "plasmid",
        )
        for contig in lengths
    ]
    genes, tus, replicons = od.placed(specs)
    out = od.rows(genes, tus, replicons, feats, od.config(bootstrap=20, permutations=9))
    assert set(out["contig"].astype(str)) == {"CP171242.1"}
    cols = od.gene_columns(genes, replicons, feats).set_index("locus_tag")
    plasmid_tags = genes.loc[genes["contig"] != "CP171242.1", "locus_tag"]
    assert all(is_na(cols.loc[t, "replichore"]) for t in plasmid_tags)
    assert all(is_na(cols.loc[t, "replichore_position"]) for t in plasmid_tags)
    chrom_tags = genes.loc[genes["contig"] == "CP171242.1", "locus_tag"]
    assert all(str(cols.loc[t, "replichore"]) == "right" for t in chrom_tags)


# --------------------------------------------------------------------------- example 5


def test_example5_sim100k_has_no_oric_row_only() -> None:
    """sim100k calls no oriC: one pooled row with ``no_oric``, gene columns NA (degenerate 1)."""
    feats, reps = od.recorded("sim100k")
    assert not len(feats[feats["type"].astype(str) == "oric"])
    lengths = dict(zip(reps["contig"].astype(str), reps["length"], strict=True))
    genes, tus, replicons = od.placed(
        [
            od.Placed("chr1", int(lengths["chr1"]), [10_000, 30_000, 60_000, 90_000]),
            od.Placed("plasmid_1", int(lengths["plasmid_1"]), [1_000, 4_000], kind="plasmid"),
        ]
    )
    out = od.rows(genes, tus, replicons, feats, od.config(bootstrap=20, permutations=9))
    assert len(out) == 1
    row = od.row_of(out, "pooled")
    assert str(row["contig"]) == "chr1"
    assert is_na(row["value"]) and is_na(row["ci_low"]) and is_na(row["ci_high"])
    assert flag_set(row["flags"]) == {od.DESCRIPTIVE, od.NO_ORIC}
    assert int(row["n_perm"]) == 0
    cols = od.gene_columns(genes, replicons, feats)
    assert cols["replichore"].map(is_na).all()
    assert cols["replichore_position"].map(is_na).all()


def test_no_features_table_at_all_is_no_oric() -> None:
    """``features = None`` (no ``features.parquet``) behaves like no oriC call."""
    genes, tus, replicons = od.placed([od.Placed("chr1", 10_000, [1_000, 3_000, 6_000, 9_000])])
    out = od.rows(genes, tus, replicons, None, od.config(bootstrap=20, permutations=9))
    assert len(out) == 1
    assert flag_set(od.row_of(out, "pooled")["flags"]) == {od.DESCRIPTIVE, od.NO_ORIC}
    assert _oriter().oric_position(None, "chr1", 10_000) is None


# --------------------------------------------------------------------------- example 6


def _example6_frames() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Five right-arm TUs at x = 0, 0.25, 0.5, 0.75, 1 with y = 10, 9.4, 9.1, 8.4, 8.1."""
    length, oric = 2000, 500.0
    mids = od.right_arm_mids(oric, length, [0.0, 0.25, 0.5, 0.75, 1.0])
    genes, tus, replicons = od.placed(
        [od.Placed("chr1", length, mids, values=[10.0, 9.4, 9.1, 8.4, 8.1])]
    )
    return genes, tus, replicons, od.oric_at("chr1", oric, width=40)


def test_example6_right_slope_is_minus_one_point_nine_two() -> None:
    """x = [0, .25, .5, .75, 1], y = [10, 9.4, 9.1, 8.4, 8.1] -> slope -1.92 (ter minus ori)."""
    genes, tus, replicons, feats = _example6_frames()
    ot = od.oracle_tus(genes, tus, replicons, feats, "chr1")
    assert list(ot.arm) == ["right"] * 5
    assert list(ot.x) == pytest.approx([0.0, 0.25, 0.5, 0.75, 1.0])
    out = od.rows(genes, tus, replicons, feats, od.config(bootstrap=40, permutations=9))
    right = od.row_of(out, "right")
    assert float(right["value"]) == pytest.approx(-1.92, abs=5e-3)
    assert float(right["effect"]) == pytest.approx(float(right["value"]))
    # with every TU on the right arm the pooled fit is the same line
    assert float(od.row_of(out, "pooled")["value"]) == pytest.approx(-1.92, abs=5e-3)
    assert is_na(od.row_of(out, "left")["value"])
    assert flag_set(od.row_of(out, "left")["flags"]) == {od.DESCRIPTIVE, od.NOT_ESTIMABLE}


def test_example6_ci_brackets_the_slope_and_matches_the_oracle() -> None:
    """The bootstrap CI is the 95 % percentile interval of the refitted slopes."""
    genes, tus, replicons, feats = _example6_frames()
    cfg = od.config(bootstrap=40, permutations=9, block_len=2)
    got = od.rows(genes, tus, replicons, feats, cfg)
    want = od.oracle_rows(genes, tus, replicons, feats, cfg)
    pd.testing.assert_frame_equal(
        got.reset_index(drop=True), want.reset_index(drop=True), check_dtype=False
    )
    row = od.row_of(got, "right")
    assert float(row["ci_low"]) <= float(row["value"]) <= float(row["ci_high"])


# --------------------------------------------------------------------------- example 8


def _mixed_genome() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """A circular chromosome with both arms populated, a plasmid and a linear contig."""
    length, oric = 20_000, 2_000.0
    rng = np.random.default_rng(3)
    pos = [0.05, 0.22, 0.35, 0.55, 0.7, 0.9, 0.98]
    mids = od.right_arm_mids(oric, length, pos) + [(oric - p * length / 2.0) % length for p in pos]
    values = list(20.0 - 2.0 * np.asarray(pos + pos) + 0.1 * rng.standard_normal(len(pos) * 2))
    specs = [
        od.Placed("chr1", length, mids, values=values),
        od.Placed("p1", 4_000, [500, 1_500, 2_500, 3_500], kind="plasmid"),
        od.Placed("c2", 3_000, [500, 1_500, 2_500], kind="unknown", topology="linear"),
    ]
    genes, tus, replicons = od.placed(specs)
    return genes, tus, replicons, od.oric_at("chr1", oric, width=40)


def test_example8_every_row_is_descriptive_without_a_p_value() -> None:
    """§7.3.5(e): ``p_value`` NA, ``q_value`` NA and the flag ``descriptive`` on every row."""
    genes, tus, replicons, feats = _mixed_genome()
    out = od.rows(genes, tus, replicons, feats, od.config(bootstrap=40, permutations=9))
    assert len(out) == 3
    for _, row in out.iterrows():
        assert is_na(row["p_value"]) and is_na(row["q_value"])
        assert od.DESCRIPTIVE in flag_set(row["flags"])
        assert "uncalibrated" not in flag_set(row["flags"])  # ADR-0021 Amendment 1
        assert str(row["statistic"]) == od.STATISTIC
        assert str(row["null_model"]) == od.NULL_MODEL
        assert str(row["method"]) == od.METHOD
        assert is_na(row["lag"])
        assert is_na(row["null_mean"]) and is_na(row["null_lo"]) and is_na(row["null_hi"])


def test_rows_use_the_autocorrelation_columns_in_order() -> None:
    """The frame is §9.3 ``autocorrelation.tsv`` shaped (M5-05 ``AUTOCORR_COLUMNS``)."""
    from giap.context.autocorrelation import AUTOCORR_COLUMNS

    genes, tus, replicons, feats = _mixed_genome()
    out = od.rows(genes, tus, replicons, feats, od.config(bootstrap=20, permutations=9))
    assert list(out.columns) == list(AUTOCORR_COLUMNS)
    assert list(AUTOCORR_COLUMNS) == list(od.COLUMNS)
    assert sorted(out["stat"].astype(str)) == [f"{od.STATISTIC}_{a}" for a in sorted(od.ARMS)]


def test_only_circular_chromosomes_are_analysed() -> None:
    """Plasmids, linear and ``unknown``-kind contigs: no row, NA gene columns (decision 33)."""
    genes, tus, replicons, feats = _mixed_genome()
    out = od.rows(genes, tus, replicons, feats, od.config(bootstrap=20, permutations=9))
    assert set(out["contig"].astype(str)) == {"chr1"}
    cols = od.gene_columns(genes, replicons, feats).set_index("locus_tag")
    off = genes.loc[genes["contig"] != "chr1", "locus_tag"]
    assert all(is_na(cols.loc[t, "replichore"]) for t in off)
    on = genes.loc[genes["contig"] == "chr1", "locus_tag"]
    assert all(0.0 <= float(cols.loc[t, "replichore_position"]) <= 1.0 for t in on)


def test_matches_the_oracle_on_a_mixed_genome() -> None:
    """Whole-frame parity with the spec oracle (values, CIs, n_perm, flags and draw order)."""
    genes, tus, replicons, feats = _mixed_genome()
    cfg = od.config(bootstrap=60, permutations=19)
    got = od.rows(genes, tus, replicons, feats, cfg)
    want = od.oracle_rows(genes, tus, replicons, feats, cfg)
    pd.testing.assert_frame_equal(
        got.reset_index(drop=True), want.reset_index(drop=True), check_dtype=False
    )


def test_same_seed_twice_is_identical() -> None:
    """§6.6 / R7: the rows depend on the seed only."""
    genes, tus, replicons, feats = _mixed_genome()
    cfg = od.config(bootstrap=40, permutations=9)
    a = od.rows(genes, tus, replicons, feats, cfg)
    b = od.rows(genes, tus, replicons, feats, cfg)
    pd.testing.assert_frame_equal(a, b)
    c = od.rows(genes, tus, replicons, feats, cfg, seed=cfg.seed + 1)
    assert not np.allclose(a["ci_low"].to_numpy(float), c["ci_low"].to_numpy(float), equal_nan=True)


# --------------------------------------------------------------------------- TU values


def test_tu_value_is_the_mean_of_its_detected_genes() -> None:
    """y_t = mean ``abundance_log2`` of the TU's detected genes with a value (decision 34)."""
    length, oric = 8_000, 0.0
    mids = od.right_arm_mids(oric, length, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    genes, tus, replicons = od.placed(
        [
            od.Placed(
                "chr1",
                length,
                mids,
                values=[10.0, 12.0, 9.0, 7.0, 8.0, 6.0],
                detected=[1, 1, 1, 0, 1, 1],
                tu=["T1", "T1", "T2", "T2", "T3", "T4"],
            )
        ]
    )
    vals = od.tu_values(genes, "chr1")
    assert vals == pytest.approx({"T1": 11.0, "T2": 9.0, "T3": 8.0, "T4": 6.0})
    ot = od.oracle_tus(genes, tus, replicons, od.oric_at("chr1", oric, width=40), "chr1")
    assert list(ot.tu_id) == ["T1", "T2", "T3", "T4"]
    assert list(ot.y) == pytest.approx([11.0, 9.0, 8.0, 6.0])


def test_tu_without_a_value_is_dropped() -> None:
    """A TU whose genes are all undetected contributes no point."""
    length = 8_000
    mids = od.right_arm_mids(0.0, length, [0.1, 0.3, 0.5, 0.7])
    genes, tus, replicons = od.placed([od.Placed("chr1", length, mids, detected=[1, 0, 1, 1])])
    ot = od.oracle_tus(genes, tus, replicons, od.oric_at("chr1", 0, width=40), "chr1")
    assert len(ot.tu_id) == 3
    out = od.rows(genes, tus, replicons, od.oric_at("chr1", 0, width=40),
                  od.config(bootstrap=20, permutations=9))  # fmt: skip
    assert float(od.row_of(out, "pooled")["value"]) == pytest.approx(od.ols_slope(ot.x, ot.y))


# --------------------------------------------------------------------------- degenerate table


def test_degenerate_fit_with_two_tus_is_not_estimable() -> None:
    """A fit with fewer than 3 TUs with values: ``value`` and CI NA, ``not_estimable``."""
    length, oric = 20_000, 5_000.0
    right = od.right_arm_mids(oric, length, [0.2, 0.6])
    left = [(oric - p * length / 2.0) % length for p in (0.1, 0.4, 0.7, 0.95)]
    genes, tus, replicons = od.placed([od.Placed("chr1", length, right + left)])
    feats = od.oric_at("chr1", oric, width=40)
    out = od.rows(genes, tus, replicons, feats, od.config(bootstrap=40, permutations=9))
    r = od.row_of(out, "right")
    assert is_na(r["value"]) and is_na(r["ci_low"]) and is_na(r["ci_high"])
    assert flag_set(r["flags"]) == {od.DESCRIPTIVE, od.NOT_ESTIMABLE}
    assert int(r["n_perm"]) == 0
    for arm in ("pooled", "left"):
        assert not is_na(od.row_of(out, arm)["value"])
        assert flag_set(od.row_of(out, arm)["flags"]) == {od.DESCRIPTIVE}


def test_degenerate_zero_variance_of_x_is_not_estimable() -> None:
    """Three concentric TUs share one midpoint, so x has no variance: NA with
    ``not_estimable``."""
    length, oric = 10_000, 0.0
    genes, tus, replicons = od.placed([od.Placed("chr1", length, [2_500] * 3, width=100)])
    genes.loc[1, ["start", "end"]] = [2_300, 2_700]  # concentric, so every TU midpoint is 2,500
    genes.loc[2, ["start", "end"]] = [2_100, 2_900]
    tus = od.tu_frame(genes)
    feats = od.oric_at("chr1", oric, width=40)
    ot = od.oracle_tus(genes, tus, replicons, feats, "chr1")
    assert float(np.var(ot.x)) == 0.0
    out = od.rows(genes, tus, replicons, feats, od.config(bootstrap=40, permutations=9))
    for arm in od.ARMS:
        row = od.row_of(out, arm)
        assert is_na(row["value"])
        assert flag_set(row["flags"]) == {od.DESCRIPTIVE, od.NOT_ESTIMABLE}


def test_degenerate_no_detected_gene_gives_three_na_rows() -> None:
    """No detected gene on the chromosome: all three rows NA with ``not_estimable``."""
    length, oric = 10_000, 0.0
    mids = od.right_arm_mids(oric, length, [0.1, 0.3, 0.5, 0.7])
    genes, tus, replicons = od.placed([od.Placed("chr1", length, mids, detected=[0, 0, 0, 0])])
    feats = od.oric_at("chr1", oric, width=40)
    out = od.rows(genes, tus, replicons, feats, od.config(bootstrap=40, permutations=9))
    assert len(out) == 3
    for arm in od.ARMS:
        row = od.row_of(out, arm)
        assert is_na(row["value"]) and int(row["n_perm"]) == 0
        assert flag_set(row["flags"]) == {od.DESCRIPTIVE, od.NOT_ESTIMABLE}
    # the gene columns still exist: the oriC is known
    cols = od.gene_columns(genes, replicons, feats)
    assert not cols["replichore"].map(is_na).any()


def _tie_frames(length: int = 10_000) -> tuple[Any, ...]:
    """Four right-arm TUs, three of them at the same position: many bootstrap replicates have
    no variance in x, so their slope is NaN and drops out of ``n_perm``."""
    mids = od.right_arm_mids(0.0, length, [0.2, 0.2, 0.2, 0.8])
    genes, _tus, replicons = od.placed(
        [od.Placed("chr1", length, mids, values=[10.0, 11.0, 9.0, 7.0], width=20)]
    )
    genes.loc[1, ["start", "end"]] = [int(mids[1]) - 40, int(mids[1]) + 40]
    genes.loc[2, ["start", "end"]] = [int(mids[2]) - 60, int(mids[2]) + 60]
    return genes, od.tu_frame(genes), replicons, od.oric_at("chr1", 0, width=40)


def test_degenerate_bootstrap_replicates_drop_out_of_n_perm() -> None:
    """A replicate whose resampled x has no variance gives a NaN slope and is not counted."""
    genes, tus, replicons, feats = _tie_frames()
    cfg = od.config(bootstrap=199, permutations=9, block_len=1)
    out = od.rows(genes, tus, replicons, feats, cfg)
    want = od.oracle_rows(genes, tus, replicons, feats, cfg)
    row, oracle = od.row_of(out, "pooled"), od.row_of(want, "pooled")
    assert int(row["n_perm"]) == int(oracle["n_perm"])
    assert 0 < int(row["n_perm"]) < int(cfg.stats.bootstrap)
    assert int(row["n_perm"]) >= cfg.stats.bootstrap / 2
    assert not is_na(row["ci_low"]) and od.NOT_ESTIMABLE not in flag_set(row["flags"])


def test_degenerate_fewer_than_half_valid_replicates_drops_the_ci() -> None:
    """Below B/2 valid replicates the CI is NA and ``not_estimable`` is added (``value`` stays)."""
    genes, tus, replicons, feats = _tie_frames()
    seeds = []
    for seed in range(200):
        cfg = od.config(bootstrap=4, permutations=9, block_len=1, seed=seed)
        want = od.row_of(od.oracle_rows(genes, tus, replicons, feats, cfg), "pooled")
        if int(want["n_perm"]) < 2:
            seeds.append(seed)
            break
    assert seeds, "no seed in 0..199 gave fewer than B/2 valid replicates"
    cfg = od.config(bootstrap=4, permutations=9, block_len=1, seed=seeds[0])
    row = od.row_of(od.rows(genes, tus, replicons, feats, cfg), "pooled")
    assert not is_na(row["value"])
    assert is_na(row["ci_low"]) and is_na(row["ci_high"])
    assert flag_set(row["flags"]) == {od.DESCRIPTIVE, od.NOT_ESTIMABLE}


# --------------------------------------------------------------------------- block length


def test_integer_block_len_overrides_the_rule() -> None:
    """ADR-0021 decision 1: an integer ``context.islands.block_len`` is used unchanged."""
    genes, tus, replicons, feats = _mixed_genome()
    cfg = od.config(bootstrap=40, permutations=9, block_len=3)
    got = od.rows(genes, tus, replicons, feats, cfg)
    want = od.oracle_rows(genes, tus, replicons, feats, cfg)
    pd.testing.assert_frame_equal(
        got.reset_index(drop=True), want.reset_index(drop=True), check_dtype=False
    )


def test_auto_block_len_follows_block_len_prime() -> None:
    """``auto``: l-hat from ``choose_block_len`` on the pooled residuals, then l' (lead 35)."""
    from giap.rng import for_step
    from giap.stats.autocorr import block_len_prime

    genes, tus, replicons, feats = _mixed_genome()
    cfg = od.config(bootstrap=40, permutations=19)
    ot = od.oracle_tus(genes, tus, replicons, feats, "chr1")
    _, diag = od.oracle_replicon(ot, cfg, for_step(od.STEP, cfg.seed))
    assert diag["l_hat"] is not None and diag["l_hat"] >= 1
    pooled = od.arm_order(ot, "pooled")
    assert diag["block_len"]["pooled"] == block_len_prime(diag["l_hat"], pooled.size)
    assert all(v <= 20 for v in diag["block_len"].values())
    got = od.rows(genes, tus, replicons, feats, cfg)
    want = od.oracle_rows(genes, tus, replicons, feats, cfg)
    pd.testing.assert_frame_equal(
        got.reset_index(drop=True), want.reset_index(drop=True), check_dtype=False
    )


def test_block_len_never_exceeds_the_fit_sequence() -> None:
    """Every fit's block length is at most its own number of TUs (the bootstrap needs it)."""
    genes, tus, replicons, feats = _mixed_genome()
    cfg = od.config(bootstrap=20, permutations=9, block_len=50)
    ot = od.oracle_tus(genes, tus, replicons, feats, "chr1")
    rows, diag = od.oracle_replicon(ot, cfg, np.random.default_rng(0))
    assert len(rows) == 3
    for arm, ell in diag["block_len"].items():
        assert 1 <= ell <= od.arm_order(ot, arm).size
    od.rows(genes, tus, replicons, feats, cfg)  # the step itself must not raise
