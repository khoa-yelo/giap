"""M5-05 autocorrelation rows: spec examples 1, 2, 4-8 and the degenerate table.

docs/specs/M5-05-autocorrelation.md with lead decisions 14-22, ADR-0020 (+ Amendment 1:
the Pearson residual join count) and docs/methods/autocorrelation.md "Statistics and nulls".

Decisions that supersede the spec text and are tested here:

* 21/22: a replicon whose 10 % span gives fewer than 4 points (n <= 30 values) has degenerate
  (all-zero) residuals, so its ``morans_i`` rows are NA with ``few_values`` and ``not_estimable``
  and do not enter a BH family. This supersedes example 1 (10 genes, I_1 = -1 -- kept here as the
  engine fact of Appendix B.4) and the Moran half of example 6 and of the "circular N = 30" row.
* 22: the values are detrended over the ordinals 0 ... n - 1 of the genes with values
  (``period`` = n when circular), after the wide 30 % pre-trend of ADR-0020 decision 1, and every
  permutation replicate is detrended again with the same smoother.
* ADR-0020 decision 4 / lead decision 15: ``ci_low`` / ``ci_high`` are NA on every row, flagged
  ``ci_unestimated``.
* ADR-0020 Amendment 1 / lead decision 17: ``value`` of a ``join_count_bb`` row is the raw BB_k;
  the test statistic (hence ``null_mean``, ``null_lo``, ``null_hi``, ``effect``) is the Pearson
  residual sum.
"""

from __future__ import annotations

import math
import sys
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import autocorr_data as ad
from autocorr_data import flag_set, is_na

LAGS = list(range(1, 21))


def _rows(contigs: list[Any], permutations: int = 9, seed: int = 7) -> pd.DataFrame:
    genes, replicons = ad.genome(contigs)
    return ad.rows(genes, replicons, ad.config(permutations), seed)


@cache
def _big_chromosome(n: int = 1000, n_detected: int = 700) -> Any:
    rng = np.random.default_rng(5)
    values = ad.ramp(n, rng)
    detected = [1] * n_detected + [0] * (n - n_detected)
    rng.shuffle(detected)
    return ad.series_contig(
        "chr1", values, detected=detected, p=list(rng.uniform(0.2, 0.9, n)), tu_size=3
    )


@cache
def _big_frame() -> pd.DataFrame:
    return _rows([_big_chromosome()])


# --------------------------------------------------------------------------- example 1


def test_example1_alternating_series_is_moran_minus_one_in_the_engine() -> None:
    """Appendix B.4 / §7.3.6: the alternating series has I_1 = -1 at lag 1 (circular, n = 10)."""
    from gippa.stats.autocorr import morans_i

    x = np.array([1.0, -1.0] * 5)
    assert morans_i(x, 1, True) == pytest.approx(-1.0)
    assert morans_i(x, 2, True) == pytest.approx(1.0)


def test_example1_ten_gene_replicon_moran_rows_are_few_values() -> None:
    """Decisions 21/22 supersede example 1: at n = 10 the 10 % span is 3 points, so the
    detrended values are identically 0 and the ``morans_i`` rows are NA."""
    out = _rows([ad.series_contig("chr1", [1.0, -1.0] * 5)])
    moran = out[out["statistic"] == "morans_i"]
    assert len(moran) == 2 * len(LAGS)
    for _, row in moran.iterrows():
        assert is_na(row["value"]) and is_na(row["p_value"]) and is_na(row["q_value"]), row
        assert {"few_values", "not_estimable"} <= flag_set(row["flags"]), row


def test_example1_lags_above_max_lag_are_not_estimable() -> None:
    """max_lag(10, circular) = 4: lags 5 ... 20 exist as rows but carry no statistic."""
    out = _rows([ad.series_contig("chr1", [1.0, -1.0] * 5)])
    jc = out[out["statistic"] == "join_count_bb"]
    for _, row in jc.iterrows():
        if int(row["lag"]) > 4:
            assert is_na(row["value"]) and is_na(row["p_value"]), row
            assert "not_estimable" in flag_set(row["flags"]), row


# --------------------------------------------------------------------------- example 2


def test_example2_row_grid_is_two_by_two_by_twenty() -> None:
    out = _big_frame()
    assert len(out) == 2 * 2 * len(LAGS) == 80
    keys = set(zip(out["statistic"], out["null_model"], out["lag"].astype(int), strict=True))
    assert keys == {(s, m, k) for s in ad.STATISTICS for m in ad.NULLS for k in LAGS}
    assert set(out["contig"]) == {"chr1"}


def test_example2_rows_are_sorted_contig_statistic_null_model_lag() -> None:
    """Example 2 "sorted contig, statistic, null_model, lag"; contigs in C1 order."""
    contigs = [_big_chromosome(), ad.series_contig("p1", ad.ramp(40), tu_size=2)]
    genes, replicons = ad.genome(contigs)
    out = ad.rows(genes, replicons, ad.config(5))
    assert list(dict.fromkeys(out["contig"])) == list(replicons["contig"])
    for contig, sub in out.groupby("contig", sort=False):
        key = list(zip(sub["statistic"], sub["null_model"], sub["lag"].astype(int), strict=True))
        assert key == sorted(key), contig


def test_columns_methods_and_test_column_conventions() -> None:
    from gippa.context.autocorrelation import AUTOCORR_COLUMNS

    out = _big_frame()
    assert tuple(AUTOCORR_COLUMNS) == ad.COLUMNS
    assert list(out.columns) == list(AUTOCORR_COLUMNS)
    assert (out["stat"] == out["statistic"]).all()
    assert set(out["null_model"]) == set(ad.NULLS)
    for _, row in out.iterrows():
        assert row["method"] == ad.METHOD[row["statistic"]], row
        assert is_na(row["ci_low"]) and is_na(row["ci_high"]), row
        assert "ci_unestimated" in flag_set(row["flags"]), row


# --------------------------------------------------------------------------- examples 4, 5, 7


def test_example4_q_value_is_above_p_value_and_na_together() -> None:
    out = _big_frame()
    p = out["p_value"].to_numpy(float)
    q = out["q_value"].to_numpy(float)
    assert np.array_equal(np.isnan(p), np.isnan(q))
    ok = ~np.isnan(p)
    assert np.all(q[ok] >= p[ok] - 1e-12)
    assert np.all((q[ok] >= 0.0) & (q[ok] <= 1.0))


def test_example4_q_value_is_bh_within_contig_statistic_null_model() -> None:
    from gippa.stats.multitest import bh

    out = _big_frame()
    for _, sub in out.groupby(["contig", "statistic", "null_model"], sort=False):
        want = bh(sub["p_value"].to_numpy(float))
        got = sub["q_value"].to_numpy(float)
        assert np.allclose(got, want, equal_nan=True, atol=1e-12)


def test_example5_n_perm_is_the_configured_permutation_count() -> None:
    out = _rows([ad.series_contig("chr1", ad.ramp(60), tu_size=3)], permutations=199)
    estimable = out[~out["p_value"].isna()]
    assert len(estimable) > 0
    assert set(estimable["n_perm"].astype(int)) == {199}
    assert (out["n_perm"].astype(int) <= 199).all()


def test_example7_same_seed_gives_the_same_frame() -> None:
    contigs = [
        ad.series_contig("chr1", ad.ramp(60), tu_size=3),
        ad.series_contig("p1", ad.ramp(44), tu_size=4),
    ]
    genes, replicons = ad.genome(contigs)
    cfg = ad.config(19)
    a = ad.rows(genes, replicons, cfg, seed=3)
    b = ad.rows(genes, replicons, cfg, seed=3)
    assert not ad.differences(a, b, tol=0.0)
    c = ad.rows(genes, replicons, cfg, seed=4)
    assert ad.differences(a, c, tol=0.0), "a different seed must redraw the permutations"


# --------------------------------------------------------------------------- example 6


def test_example6_small_plasmid_rows_and_bh_family() -> None:
    """30-gene circular plasmid next to a 1,000-gene chromosome: the plasmid's join-count lags
    15 ... 20 are NA (max_lag = 14) and its BH family has 14 tests; its Moran rows are
    ``few_values`` (decision 21) and enter no family."""
    from gippa.stats.multitest import bh

    rng = np.random.default_rng(11)
    plasmid = ad.series_contig(
        "p1",
        ad.ramp(30, rng),
        detected=[1, 1, 0] * 10,
        p=list(rng.uniform(0.2, 0.8, 30)),
        tu_size=2,
    )
    genes, replicons = ad.genome([_big_chromosome(), plasmid])
    out = ad.rows(genes, replicons, ad.config(49))

    for null_model in ad.NULLS:
        moran = ad.block(out, "p1", "morans_i", null_model)
        assert moran["value"].isna().all() and moran["q_value"].isna().all()
        assert all("few_values" in flag_set(f) for f in moran["flags"])

        jc = ad.block(out, "p1", "join_count_bb", null_model).sort_values("lag")
        high = jc[jc["lag"].astype(int) > 14]
        assert high["value"].isna().all() and high["p_value"].isna().all()
        assert all("not_estimable" in flag_set(f) for f in high["flags"])
        low = jc[jc["lag"].astype(int) <= 14]
        assert low["p_value"].notna().all(), low
        assert int(low["p_value"].notna().sum()) == 14
        assert np.allclose(
            jc["q_value"].to_numpy(float),
            bh(jc["p_value"].to_numpy(float)),
            equal_nan=True,
            atol=1e-12,
        )


def test_few_values_boundary_is_thirty_one_values() -> None:
    """k = ceil(0.1 n) >= 4 from n = 31 values on (``loess_span_points``)."""
    assert ad.few_values(30) and not ad.few_values(31)
    out = _rows([ad.series_contig("chr1", ad.ramp(31), tu_size=3)], permutations=19)
    moran = ad.block(out, "chr1", "morans_i", "gene_label")
    estimable = moran[moran["lag"].astype(int) <= 15]
    assert estimable["value"].notna().all()
    assert all("few_values" not in flag_set(f) for f in moran["flags"])


# --------------------------------------------------------------------------- example 8


def test_example8_no_detected_gene() -> None:
    """Moran rows NA (no values); join-count ``value`` 0 and, with a constant ``p_detect``, a
    constant null, so ``p_value`` is NA and the row is ``not_estimable``."""
    out = _rows([ad.series_contig("chr1", [0.0] * 40, detected=[0] * 40, p=0.4, tu_size=2)])
    moran = out[out["statistic"] == "morans_i"]
    assert moran["value"].isna().all() and moran["p_value"].isna().all()
    assert all("not_estimable" in flag_set(f) for f in moran["flags"])

    jc = out[(out["statistic"] == "join_count_bb") & (out["lag"].astype(int) <= 19)]
    assert (jc["value"].astype(float) == 0.0).all()
    assert jc["p_value"].isna().all()
    assert all("not_estimable" in flag_set(f) for f in jc["flags"])


# --------------------------------------------------------------------------- degenerate table


def test_degenerate_replicon_with_no_eligible_gene() -> None:
    """Every gene partial: no eligible gene, all rows NA with ``n_perm`` 0."""
    import region_tests_data as rd

    contig = rd.Contig(name="chr1", d=[0] * 12, p=[0.5] * 12, partial=[True] * 12)
    out = _rows([contig])
    assert len(out) == 80
    for col in ("value", "null_mean", "null_lo", "null_hi", "effect", "p_value", "q_value"):
        assert out[col].isna().all(), col
    assert set(out["n_perm"].astype(int)) == {0}
    assert all("not_estimable" in flag_set(f) for f in out["flags"])


def test_degenerate_every_eligible_gene_detected() -> None:
    """Moran rows normal; the join count is the number of lag-k pairs and its null is constant
    under a constant ``p_detect``, so ``p_value`` is NA."""
    out = _rows([ad.series_contig("chr1", ad.ramp(60), detected=[1] * 60, p=0.3, tu_size=3)])
    moran = ad.block(out, "chr1", "morans_i", "gene_label")
    assert moran[moran["lag"].astype(int) <= 29]["p_value"].notna().all()

    jc = ad.block(out, "chr1", "join_count_bb", "gene_label")
    pairs = {k: 60 for k in range(1, 30)}  # circular: n pairs at every valid lag
    for _, row in jc.iterrows():
        lag = int(row["lag"])
        assert float(row["value"]) == pairs[lag], row
        assert is_na(row["p_value"]) and "not_estimable" in flag_set(row["flags"]), row


def test_degenerate_one_gene_with_a_value() -> None:
    out = _rows([ad.series_contig("chr1", ad.ramp(40), detected=[1] + [0] * 39, p=0.5)])
    moran = out[out["statistic"] == "morans_i"]
    assert moran["value"].isna().all() and moran["p_value"].isna().all()
    assert all("not_estimable" in flag_set(f) for f in moran["flags"])


def test_degenerate_all_values_equal_after_detrending() -> None:
    """A linear series on a linear contig is reproduced exactly by the local-linear LOESS, so
    every residual is 0: Moran's I has no variance (and k = 4, so this is not ``few_values``)."""
    values = [10.0 + 0.5 * i for i in range(40)]
    contig = ad.series_contig("chr1", values, circular=False, tu_size=2)
    out = _rows([contig])
    moran = out[out["statistic"] == "morans_i"]
    assert moran["value"].isna().all(), moran[["lag", "null_model", "value"]]
    assert all("not_estimable" in flag_set(f) for f in moran["flags"])
    assert all("few_values" not in flag_set(f) for f in moran["flags"])


def test_degenerate_single_tu_replicon() -> None:
    """One TU: the TU-block null is the identity, so its rows are ``not_estimable``; the
    gene-label rows are computed normally."""
    rng = np.random.default_rng(3)
    values = ad.ramp(60, rng)
    contig = ad.series_contig("chr1", values, p=list(rng.uniform(0.3, 0.8, 60)), tu_size=60)
    out = _rows([contig], permutations=19)
    for statistic in ad.STATISTICS:
        tu = ad.block(out, "chr1", statistic, "tu_block")
        assert tu["p_value"].isna().all(), tu[["lag", "p_value"]]
        assert all("not_estimable" in flag_set(f) for f in tu["flags"])
        label = ad.block(out, "chr1", statistic, "gene_label")
        valid = label[label["lag"].astype(int) <= 29]
        assert valid["p_value"].notna().all(), valid[["lag", "p_value"]]


def test_degenerate_circular_replicon_with_thirty_genes() -> None:
    """N = 30 circular: lags 1-14 are computed for the join count, 15-20 are NA; the Moran rows
    are ``few_values`` (decisions 21/22 supersede this row of the table)."""
    rng = np.random.default_rng(7)
    contig = ad.series_contig(
        "chr1", ad.ramp(30, rng), detected=[1, 0] * 15, p=list(rng.uniform(0.2, 0.8, 30)), tu_size=2
    )
    out = _rows([contig], permutations=19)
    jc = ad.block(out, "chr1", "join_count_bb", "gene_label").sort_values("lag")
    assert jc[jc["lag"].astype(int) <= 14]["p_value"].notna().all()
    assert jc[jc["lag"].astype(int) >= 15]["value"].isna().all()
    moran = ad.block(out, "chr1", "morans_i", "gene_label")
    assert all({"few_values", "not_estimable"} <= flag_set(f) for f in moran["flags"])


# --------------------------------------------------------------------------- statistics


def test_moran_values_follow_decision_22() -> None:
    """``moran_values`` = the wide pre-trend then the 10 % detrend over the ordinals of the genes
    with values; NaN at genes without a value."""
    from gippa.context.autocorrelation import moran_values
    from gippa.context.regions import replicon_index

    rng = np.random.default_rng(13)
    values = ad.ramp(60, rng)
    detected = [1] * 45 + [0] * 15
    rng.shuffle(detected)
    genes, _ = ad.genome([ad.series_contig("chr1", values, detected=detected, tu_size=3)])
    ix = replicon_index(genes, "chr1")
    x = np.asarray(moran_values(ix, True), dtype=float)
    a = ix["a"].to_numpy(float)
    has = np.isfinite(a)
    assert x.shape == (len(ix),)
    assert np.array_equal(np.isfinite(x), has)
    assert x[has] == pytest.approx(ad.oracle_values(a[has], True), abs=1e-9)


def test_moran_values_are_degenerate_below_four_span_points() -> None:
    """M5-01 decision 4: with k = 3 every residual is 0 (consumers test k < 4, never x == 0)."""
    from gippa.context.autocorrelation import moran_values
    from gippa.context.regions import replicon_index

    genes, _ = ad.genome([ad.series_contig("chr1", ad.ramp(30, np.random.default_rng(2)))])
    x = np.asarray(moran_values(replicon_index(genes, "chr1"), True), dtype=float)
    assert np.allclose(x, 0.0, atol=1e-12)


def test_join_count_value_is_the_raw_count_and_effect_is_the_pearson_statistic() -> None:
    """ADR-0020 Amendment 1: ``value`` = BB_k, ``effect`` = Pearson statistic - null mean."""
    from gippa.context.regions import replicon_index
    from gippa.stats.autocorr import join_count

    rng = np.random.default_rng(17)
    detected = list(rng.integers(0, 2, 60))
    contig = ad.series_contig(
        "chr1", ad.ramp(60, rng), detected=detected, p=list(rng.uniform(0.2, 0.8, 60)), tu_size=3
    )
    genes, replicons = ad.genome([contig])
    out = ad.rows(genes, replicons, ad.config(19))
    ix = replicon_index(genes, "chr1")
    d, p = ix["d"].to_numpy(bool), ix["p"].to_numpy(float)
    for null_model in ad.NULLS:
        for lag in (1, 7, 20):
            row = ad.pick(out, "chr1", "join_count_bb", null_model, lag)
            assert float(row["value"]) == pytest.approx(join_count(d, lag, True))
            assert float(row["value"]) == float(int(row["value"]))  # a count
            pearson = join_count(d, lag, True, p=p)
            assert float(row["effect"]) == pytest.approx(
                pearson - float(row["null_mean"]), abs=1e-9
            )
            assert float(row["value"]) != pytest.approx(float(row["null_mean"]), abs=1e-6)


def test_moran_effect_is_the_value_minus_the_null_mean() -> None:
    out = _big_frame()
    moran = out[(out["statistic"] == "morans_i") & out["p_value"].notna()]
    assert len(moran) > 0
    assert np.allclose(
        moran["effect"].to_numpy(float),
        moran["value"].to_numpy(float) - moran["null_mean"].to_numpy(float),
        atol=1e-12,
    )


def test_null_envelope_brackets_the_null_mean() -> None:
    out = _big_frame()
    ok = out[out["p_value"].notna()]
    lo, mean, hi = (ok[c].to_numpy(float) for c in ("null_lo", "null_mean", "null_hi"))
    assert np.all(lo <= mean + 1e-12) and np.all(mean <= hi + 1e-12)


# --------------------------------------------------------------------------- oracle


def test_frame_matches_the_naive_oracle() -> None:
    """Every cell of a two-replicon genome, recomputed from the merged engines in the spec's
    draw order (label matrix for x, label matrix for d, TU-block matrix, replicons in C1 order)."""
    rng = np.random.default_rng(23)
    detected = [1] * 48 + [0] * 12
    rng.shuffle(detected)
    contigs = [
        ad.series_contig(
            "chr1",
            ad.ramp(60, rng),
            detected=detected,
            p=list(rng.uniform(0.2, 0.9, 60)),
            tu_size=3,
        ),
        ad.series_contig(
            "p1", ad.ramp(44, rng), circular=False, p=list(rng.uniform(0.3, 0.7, 44)), tu_size=4
        ),
    ]
    genes, replicons = ad.genome(contigs)
    cfg = ad.config(29)
    got = ad.rows(genes, replicons, cfg, seed=5)
    want = ad.oracle_rows(genes, replicons, cfg, seed=5)
    assert not ad.differences(got, want)


def test_replicon_autocorrelation_matches_the_oracle_for_one_replicon() -> None:
    from gippa.context.autocorrelation import replicon_autocorrelation
    from gippa.context.regions import replicon_index
    from gippa.rng import for_step

    rng = np.random.default_rng(31)
    contig = ad.series_contig(
        "chr1",
        ad.ramp(50, rng),
        detected=list(rng.integers(0, 2, 50)),
        p=list(rng.uniform(0.2, 0.9, 50)),
        tu_size=2,
    )
    genes, _ = ad.genome([contig])
    ix = replicon_index(genes, "chr1")
    got = replicon_autocorrelation(
        ix, "chr1", True, LAGS, 29, for_step("context_stats.autocorrelation", 9)
    )
    want = ad.oracle_replicon(
        genes, "chr1", True, LAGS, 29, for_step("context_stats.autocorrelation", 9)
    )
    assert set(ad.COLUMNS) - {"q_value"} <= set(got.columns)
    assert not ad.differences(got.assign(q_value=math.nan), want.assign(q_value=math.nan))
