"""Hypothesis properties of the M5-05 autocorrelation table.

Invariants from the acceptance list, §9.2 and the lead decisions:

* the row grid is the same for every replicon (spec "Lags": invalid lags are rows with NA
  statistics) and carries no duplicate (contig, statistic, null_model, lag);
* §9.2 test columns: ``stat`` = the statistic, ``method`` per statistic, ``ci_low`` / ``ci_high``
  NA with ``ci_unestimated`` (ADR-0020 decision 4), ``p`` in [0, 1] or NA, ``n_perm`` <= B, an NA
  statistic or p always flagged ``not_estimable``, ``few_values`` exactly on the Moran rows of a
  replicon whose 10 % span is below 4 points (decisions 21/22);
* BH within (contig, statistic, null_model): q >= p, q NA exactly where p is (example 4);
* Moran's I is invariant to affine maps of the abundances (docs/methods/autocorrelation.md), and
  the permutation draws depend only on the sizes, so the whole frame is unchanged by
  ``a -> c a + b`` with c > 0;
* Moran rows do not depend on ``p_detect`` (only the join count adjusts for detectability);
* the values are rotation-equivariant on a circular replicon (circular padding, M5-01), so the
  ``morans_i`` values are unchanged when the genome is rotated;
* the raw join count is an integer in [0, number of lag-k pairs].
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import autocorr_data as ad
from autocorr_data import flag_set, is_na

B = 5
LAGS = list(range(1, 21))
SETTINGS = settings(
    max_examples=25,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
PROBS = [0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0]


@st.composite
def genomes(draw: Any, max_contigs: int = 2) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One or two hand-built replicons with varied size, detection, TU runs and topology."""
    contigs = []
    for c in range(draw(st.integers(1, max_contigs))):
        n = draw(st.integers(4, 70))
        rng = np.random.default_rng(draw(st.integers(0, 2**31 - 1)))
        detected = draw(st.lists(st.booleans(), min_size=n, max_size=n))
        probs = draw(st.lists(st.sampled_from(PROBS), min_size=n, max_size=n))
        contigs.append(
            ad.series_contig(
                f"c{c}",
                ad.ramp(n, rng),
                circular=draw(st.booleans()),
                detected=[int(x) for x in detected],
                p=probs,
                tu_size=draw(st.integers(1, 5)),
            )
        )
    return ad.genome(contigs)


def _rows(genes: pd.DataFrame, replicons: pd.DataFrame, seed: int = 2) -> pd.DataFrame:
    return ad.rows(genes, replicons, ad.config(B), seed)


def _n_values(genes: pd.DataFrame, contig: str) -> int:
    from giap.context.regions import replicon_index

    ix = replicon_index(genes, contig)
    return int(np.isfinite(ix["a"].to_numpy(float)).sum()) if len(ix) else 0


@SETTINGS
@given(genomes())
def test_row_grid_is_complete_and_unique(g: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    genes, replicons = g
    out = _rows(genes, replicons)
    want = {
        (c, s, m, k)
        for c in replicons["contig"]
        for s in ad.STATISTICS
        for m in ad.NULLS
        for k in LAGS
    }
    keys = list(
        zip(out["contig"], out["statistic"], out["null_model"], out["lag"].astype(int), strict=True)
    )
    assert len(keys) == len(set(keys)) == len(want)
    assert set(keys) == want


@SETTINGS
@given(genomes())
def test_test_columns_and_flags(g: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    from giap.context.autocorrelation import AUTOCORR_COLUMNS

    genes, replicons = g
    out = _rows(genes, replicons)
    assert list(out.columns) == list(AUTOCORR_COLUMNS)
    few = {c: ad.few_values(_n_values(genes, str(c))) for c in replicons["contig"]}
    for _, row in out.iterrows():
        flags = flag_set(row["flags"])
        assert row["stat"] == row["statistic"]
        assert row["method"] == ad.METHOD[row["statistic"]]
        assert is_na(row["ci_low"]) and is_na(row["ci_high"]) and "ci_unestimated" in flags
        p = row["p_value"]
        assert is_na(p) or 0.0 <= float(p) <= 1.0
        assert 0 <= int(row["n_perm"]) <= B
        if is_na(p) or is_na(row["value"]):
            assert "not_estimable" in flags, row
        else:
            assert "not_estimable" not in flags, row
        want_few = few[row["contig"]] and row["statistic"] == "morans_i"
        assert ("few_values" in flags) == want_few, row
        if want_few:
            assert is_na(row["value"]) and is_na(p) and is_na(row["effect"]), row


@SETTINGS
@given(genomes())
def test_bh_within_contig_statistic_and_null(g: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    from giap.stats.multitest import bh

    genes, replicons = g
    out = _rows(genes, replicons)
    p = out["p_value"].to_numpy(float)
    q = out["q_value"].to_numpy(float)
    assert np.array_equal(np.isnan(p), np.isnan(q))
    ok = ~np.isnan(p)
    assert np.all(q[ok] + 1e-12 >= p[ok])
    for _, sub in out.groupby(["contig", "statistic", "null_model"], sort=False):
        assert np.allclose(
            sub["q_value"].to_numpy(float),
            bh(sub["p_value"].to_numpy(float)),
            equal_nan=True,
            atol=1e-12,
        )


@SETTINGS
@given(genomes())
def test_raw_join_count_is_an_integer_within_its_support(
    g: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    from giap.context.regions import replicon_index
    from giap.stats.autocorr import lag_pairs

    genes, replicons = g
    out = _rows(genes, replicons)
    jc = out[out["statistic"] == "join_count_bb"]
    for _, row in jc.iterrows():
        contig = str(row["contig"])
        circular = (
            str(replicons.loc[replicons["contig"] == contig, "topology"].iloc[0]) == "circular"
        )
        n_pairs = len(lag_pairs(len(replicon_index(genes, contig)), int(row["lag"]), circular)[0])
        if n_pairs == 0:
            assert is_na(row["value"]), row
        else:
            v = float(row["value"])
            assert v == float(int(v)) and 0 <= v <= n_pairs, row


@SETTINGS
@given(genomes(max_contigs=1), st.floats(0.5, 4.0), st.floats(-5.0, 5.0))
def test_affine_abundance_leaves_the_frame_unchanged(
    g: tuple[pd.DataFrame, pd.DataFrame], scale: float, shift: float
) -> None:
    genes, replicons = g
    other = genes.copy()
    other["abundance_log2"] = scale * genes["abundance_log2"].to_numpy(float) + shift
    assert not ad.differences(_rows(genes, replicons), _rows(other, replicons), tol=1e-7)


@SETTINGS
@given(genomes(max_contigs=1), st.lists(st.sampled_from(PROBS), min_size=1))
def test_moran_rows_do_not_depend_on_p_detect(
    g: tuple[pd.DataFrame, pd.DataFrame], probs: list[float]
) -> None:
    genes, replicons = g
    other = genes.copy()
    n = len(genes)
    other["p_detect"] = np.resize(np.asarray(probs, dtype=float), n)
    a = _rows(genes, replicons)
    b = _rows(other, replicons)
    moran_a = a[a["statistic"] == "morans_i"].reset_index(drop=True)
    moran_b = b[b["statistic"] == "morans_i"].reset_index(drop=True)
    assert not ad.differences(moran_a, moran_b, tol=0.0)


@SETTINGS
@given(st.integers(4, 70), st.integers(1, 5), st.integers(1, 69), st.integers(0, 2**31 - 1))
def test_moran_values_are_rotation_equivariant_on_a_circular_replicon(
    n: int, tu_size: int, shift: int, seed: int
) -> None:
    rng = np.random.default_rng(seed)
    values = ad.ramp(n, rng)
    detected = [int(x) for x in rng.integers(0, 2, n)]
    probs = list(rng.uniform(0.1, 0.9, n))
    k = shift % n
    base = ad.series_contig("chr1", values, detected=detected, p=probs, tu_size=tu_size)
    rolled = ad.series_contig(
        "chr1",
        np.roll(values, k),
        detected=list(np.roll(detected, k)),
        p=list(np.roll(probs, k)),
        tu_size=tu_size,
    )
    # the TU labels of the rotated genome are the rotated labels (tu_size runs rotate too)
    rolled.tu = list(np.roll(np.asarray(base.tu, dtype=object), k))
    a = _rows(*ad.genome([base]))
    b = _rows(*ad.genome([rolled]))
    va = a[a["statistic"] == "morans_i"]["value"].to_numpy(float)
    vb = b[b["statistic"] == "morans_i"]["value"].to_numpy(float)
    assert np.allclose(va, vb, atol=1e-8, equal_nan=True)


@SETTINGS
@given(genomes())
def test_same_seed_same_frame(g: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    genes, replicons = g
    assert not ad.differences(_rows(genes, replicons, 11), _rows(genes, replicons, 11), tol=0.0)


@SETTINGS
@given(genomes(max_contigs=1))
def test_values_match_the_decision_22_oracle(g: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    from giap.context.autocorrelation import moran_values
    from giap.context.regions import replicon_index

    genes, replicons = g
    contig = str(replicons["contig"].iloc[0])
    circular = str(replicons["topology"].iloc[0]) == "circular"
    ix = replicon_index(genes, contig)
    x = np.asarray(moran_values(ix, circular), dtype=float)
    a = ix["a"].to_numpy(float) if len(ix) else np.zeros(0)
    has = np.isfinite(a)
    assert x.shape == (len(ix),)
    assert np.array_equal(~np.isnan(x), has), "a value exactly at the genes with an abundance"
    assert np.allclose(x[has], ad.oracle_values(a[has], circular), atol=1e-9)
