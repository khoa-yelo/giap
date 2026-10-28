"""MAYU-style gene FDR, closed form (issue M4-02; GIPPA_SPEC.md §7.3.2, §8.3, Appendix B.3).

Examples 1-6 of docs/specs/M4-02-stats-mayu.md, the degenerate guards (p_b = 1, FDR clipped to 1,
no thresholds -> empty curve), PEP vs t*|P(t)|, rank binning, edge cases and a seeded Monte Carlo
cross-check of the closed form's model (false peptides land on genes with probability ∝ a_g).

Note: the 20 % unpaired-decoy fallback to MAYU (ADR-0009, upper guard) is tier logic in M4-12;
here we test only what ``mayu_gene_fdr`` itself guarantees.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from gippa.stats.mayu import MayuResult, mayu_gene_fdr
from gippa.stats.poisson_binomial import poisson_binomial_pmf
from gippa.stats.tda import monotone_q

NAN = float("nan")
CURVE_COLUMNS = ["t", "n_peptides", "F", "n_genes", "fp", "fdr", "fdr_2pep"]


def _row(res: MayuResult, t: float) -> pd.Series:
    rows = res.curve[np.isclose(res.curve["t"].to_numpy(), t, rtol=0, atol=1e-15)]
    assert len(rows) == 1, f"expected one curve row at t={t}, got {len(rows)}"
    return rows.iloc[0]


def _b3_inputs(n_bins: int = 1) -> tuple[np.ndarray, np.ndarray, None, np.ndarray, int]:
    """B.3: 5,000 genes of equal size, 3,000 detected at q 0.01, F = 0.01 * 5,000 = 50."""
    best_q = np.full(5000, NAN)
    best_q[:3000] = 0.01
    peptide_q = np.full(5000, 0.01)
    return best_q, peptide_q, None, np.full(5000, 7.0), n_bins


# ----------------------------------------------------------------------------- examples


@pytest.mark.parametrize("n_bins", [1, 7, 10])
def test_example1_b3_closed_form_any_n_bins(n_bins: int) -> None:
    res = mayu_gene_fdr(*_b3_inputs(n_bins))
    assert isinstance(res, MayuResult)
    assert list(res.curve.columns) == CURVE_COLUMNS
    row = _row(res, 0.01)
    assert row["n_peptides"] == 5000
    assert row["F"] == pytest.approx(50.0, rel=1e-12)
    assert row["n_genes"] == 3000
    assert row["fp"] == pytest.approx(20.1024, rel=1e-4)
    assert row["fdr"] == pytest.approx(0.0067008, rel=1e-4)
    # p_b recovered from the fixed point FP = (N - n + FP) p
    p_b = row["fp"] / (2000 + row["fp"])
    assert p_b == pytest.approx(0.0099512, rel=1e-4)


def test_example1_gene_q() -> None:
    res = mayu_gene_fdr(*_b3_inputs())
    assert res.gene_q.shape == (5000,)
    assert np.allclose(res.gene_q[:3000], 0.0067008, rtol=1e-4)
    assert np.isnan(res.gene_q[3000:]).all()


def test_example1_with_peps_same_answer() -> None:
    best_q, _, _, sizes, _ = _b3_inputs()
    peptide_q = np.full(100, 0.01)
    pep = np.full(100, 0.5)  # F = 50 from PEPs
    res = mayu_gene_fdr(best_q, peptide_q, pep, sizes, 1)
    assert _row(res, 0.01)["fdr"] == pytest.approx(0.0067008, rel=1e-4)


def test_example2_two_peptide_variant() -> None:
    res = mayu_gene_fdr(*_b3_inputs())
    row = _row(res, 0.01)
    lam = 0.01
    p2 = 1.0 - math.exp(-lam) * (1.0 + lam)
    assert p2 == pytest.approx(4.9668e-5, rel=1e-4)
    assert row["fdr_2pep"] == pytest.approx(3.3114e-5, rel=1e-4)
    assert row["fdr_2pep"] * 3000 == pytest.approx(0.099341, rel=1e-4)


@pytest.mark.parametrize("n_bins", [5, 10])
def test_example3_two_size_groups(n_bins: int) -> None:
    """3,000 genes a = 100 (1,000 hit), 2,000 genes a = 400 (1,500 hit), F = 40.

    Equal-count rank bins with 5 or 10 bins align with the two size groups, so w_b and the
    summed FP equal the two-bin answer.
    """
    sizes = np.r_[np.full(3000, 100.0), np.full(2000, 400.0)]
    best_q = np.full(5000, NAN)
    # spread hits through each group so every bin gets some
    best_q[0:3000:3] = 0.01
    best_q[3000:5000][np.arange(2000) % 4 != 0] = 0.01
    assert np.isfinite(best_q[:3000]).sum() == 1000
    assert np.isfinite(best_q[3000:]).sum() == 1500
    peptide_q = np.full(4000, 0.01)  # F = 0.01 * 4,000 = 40
    res = mayu_gene_fdr(best_q, peptide_q, None, sizes, n_bins)
    row = _row(res, 0.01)
    assert row["F"] == pytest.approx(40.0, rel=1e-12)
    assert row["n_genes"] == 2500
    assert row["fp"] == pytest.approx(7.2863 + 7.3272, rel=1e-4)
    assert row["fdr"] == pytest.approx(0.0058454, rel=1e-4)
    # bins never mix the two size groups here
    assert set(res.bin_of_gene[:3000]).isdisjoint(set(res.bin_of_gene[3000:]))


def test_example4_F_from_peps() -> None:
    best_q = np.array([0.01, 0.02, NAN, NAN])
    res = mayu_gene_fdr(best_q, [0.001, 0.005, 0.02], [0.01, 0.02, 0.5], np.ones(4), 1)
    assert _row(res, 0.02)["F"] == pytest.approx(0.53, rel=1e-12)
    assert _row(res, 0.02)["n_peptides"] == 3
    assert _row(res, 0.01)["F"] == pytest.approx(0.03, rel=1e-12)
    assert _row(res, 0.01)["n_peptides"] == 2


def test_example5_F_without_peps() -> None:
    best_q = np.r_[np.full(10, 0.01), np.full(90, NAN)]
    peptide_q = np.r_[np.linspace(0.0, 0.01, 200), np.full(50, 0.2)]  # 200 with q <= 0.01
    res = mayu_gene_fdr(best_q, peptide_q, None, np.ones(100), 1)
    row = _row(res, 0.01)
    assert row["n_peptides"] == 200
    assert row["F"] == pytest.approx(2.0, rel=1e-12)


def _example6() -> MayuResult:
    """FDR curve that is not monotone in t: about (0.183, 0.020, 0.102) at t = 0.001/0.005/0.01."""
    best_q = np.full(100, NAN)
    best_q[:10] = 0.001
    best_q[10:50] = 0.005
    best_q[50] = 0.01
    peptide_q = np.array([0.001, 0.001, 0.01])
    pep = np.array([1.0, 1.0, 8.0])  # F = 2, 2, 10
    return mayu_gene_fdr(best_q, peptide_q, pep, np.ones(100), 1)


def test_example6_gene_q_is_min_over_looser_thresholds() -> None:
    res = _example6()
    fdr = res.curve.set_index("t")["fdr"]
    f1, f2, f3 = fdr.loc[0.001], fdr.loc[0.005], fdr.loc[0.01]
    # hand values: N - n = 90, 50, 49; p = 1 - 0.99^F
    for f, n, F in ((f1, 10, 2), (f2, 50, 2), (f3, 51, 10)):
        p = 1 - 0.99**F
        assert f == pytest.approx((100 - n) * p / (1 - p) / n, rel=1e-9)
    assert f1 > f2 < f3  # non-monotone curve
    assert res.gene_q[0] == pytest.approx(f2, rel=1e-12)  # best_q 0.001 -> min(f1, f2, f3)
    assert res.gene_q[10] == pytest.approx(f2, rel=1e-12)  # best_q 0.005 -> min(f2, f3)
    assert res.gene_q[50] == pytest.approx(f3, rel=1e-12)  # best_q 0.01 -> f3
    assert np.isnan(res.gene_q[51:]).all()


def test_example6_matches_monotone_q() -> None:
    res = _example6()
    curve = res.curve.sort_values("t")
    q_at_t = dict(zip(curve["t"], monotone_q(curve["fdr"].to_numpy()), strict=True))
    for g in (0, 10, 50):
        best = float(_example6_best_q()[g])
        assert res.gene_q[g] == pytest.approx(q_at_t[best], rel=1e-12)


def _example6_best_q() -> np.ndarray:
    b = np.full(100, NAN)
    b[:10], b[10:50], b[50] = 0.001, 0.005, 0.01
    return b


# ----------------------------------------------------------------------------- curve shape


def test_curve_thresholds_are_distinct_finite_best_q_ascending() -> None:
    best_q = np.array([0.02, NAN, 0.001, 0.02, 0.005, NAN])
    res = mayu_gene_fdr(best_q, [0.001, 0.004, 0.02], None, np.ones(6), 2)
    assert res.curve["t"].tolist() == [0.001, 0.005, 0.02]
    assert res.curve["n_genes"].tolist() == [1, 2, 4]
    assert res.curve["n_peptides"].tolist() == [1, 2, 3]


def test_genes_without_peptide_count_in_N_b() -> None:
    """Adding undetected genes increases N_b - n_b and therefore FP."""
    base = mayu_gene_fdr([0.01] * 10 + [NAN] * 10, [0.01] * 100, None, np.ones(20), 1)
    more = mayu_gene_fdr([0.01] * 10 + [NAN] * 40, [0.01] * 100, None, np.ones(50), 1)
    assert _row(more, 0.01)["fp"] > _row(base, 0.01)["fp"]


# ----------------------------------------------------------------------------- degenerate cases


def test_all_detected_gives_zero_fp() -> None:
    """n_b = N_b -> FP_b = 0."""
    res = mayu_gene_fdr([0.01] * 20, [0.01] * 500, None, np.arange(1, 21), 4)
    row = _row(res, 0.01)
    assert row["fp"] == 0.0 and row["fdr"] == 0.0 and row["fdr_2pep"] == 0.0
    assert np.all(res.gene_q == 0.0)


def test_F_zero_gives_zero_fdr() -> None:
    """F(t) = 0 -> p_b = 0 -> FDR 0 (all PEPs zero)."""
    res = mayu_gene_fdr([0.01, NAN, NAN], [0.01, 0.01], [0.0, 0.0], np.ones(3), 1)
    row = _row(res, 0.01)
    assert row["F"] == 0.0 and row["fp"] == 0.0 and row["fdr"] == 0.0
    assert res.gene_q[0] == 0.0


def test_p_b_one_gives_fp_n_minus_n_and_fdr_clipped() -> None:
    """w = 1/3, F = 2,000 -> (2/3)^2000 underflows -> p_b = 1 -> FP = N - n = 2, FDR 2 -> 1."""
    res = mayu_gene_fdr([0.01, NAN, NAN], [0.01] * 2000, [1.0] * 2000, np.ones(3), 1)
    row = _row(res, 0.01)
    assert row["fp"] == pytest.approx(2.0)
    assert np.isfinite(row["fp"])
    assert row["fdr"] == 1.0
    assert row["fdr_2pep"] <= 1.0
    assert res.gene_q[0] == 1.0


def test_no_detected_gene_gives_empty_curve_and_nan_q() -> None:
    res = mayu_gene_fdr([NAN, NAN, NAN], [0.01, 0.02], None, [3, 4, 5], 2)
    assert list(res.curve.columns) == CURVE_COLUMNS
    assert len(res.curve) == 0
    assert res.gene_q.shape == (3,) and np.isnan(res.gene_q).all()
    assert res.bin_of_gene.shape == (3,)


def test_no_peptides_and_no_peps() -> None:
    """No accepted peptides (e.g. no decoy-qualified identifications): F = 0, FDR 0."""
    res = mayu_gene_fdr([0.01, NAN], [], None, [1, 1], 1)
    row = _row(res, 0.01)
    assert row["n_peptides"] == 0 and row["F"] == 0.0 and row["fdr"] == 0.0


def test_single_gene() -> None:
    res = mayu_gene_fdr([0.01], [0.01] * 10, None, [5], 10)
    assert res.bin_of_gene.tolist() == [0]
    row = _row(res, 0.01)
    assert row["fp"] == 0.0 and row["fdr"] == 0.0
    assert res.gene_q.tolist() == [0.0]


def test_single_gene_undetected() -> None:
    res = mayu_gene_fdr([NAN], [0.01], None, [5], 10)
    assert len(res.curve) == 0 and np.isnan(res.gene_q).all()


def test_all_single_hit_genes() -> None:
    """Each detected gene carries one accepted peptide: |P(t)| = n(t) at every threshold."""
    rng = np.random.default_rng(7)
    best_q = np.full(300, NAN)
    best_q[:120] = np.round(rng.uniform(0.0, 0.01, 120), 4)
    peptide_q = best_q[:120].copy()
    res = mayu_gene_fdr(best_q, peptide_q, None, rng.integers(1, 50, 300), 10)
    assert (res.curve["n_peptides"].to_numpy() == res.curve["n_genes"].to_numpy()).all()
    assert np.allclose(res.curve["F"], res.curve["t"] * res.curve["n_peptides"], rtol=1e-12)
    assert (res.curve["fdr_2pep"] <= res.curve["fdr"] + 1e-15).all()


# ----------------------------------------------------------------------------- binning


def test_sizes_floored_at_one() -> None:
    best_q = [0.01, 0.01, NAN, NAN, 0.02]
    a = mayu_gene_fdr(best_q, [0.01] * 300, None, [0, -3, 1, 0.5, 1], 2)
    b = mayu_gene_fdr(best_q, [0.01] * 300, None, [1, 1, 1, 1, 1], 2)
    pd.testing.assert_frame_equal(a.curve, b.curve)
    np.testing.assert_array_equal(a.gene_q, b.gene_q)


def test_fewer_genes_than_bins() -> None:
    res = mayu_gene_fdr([0.01, NAN, 0.02], [0.01], None, [5, 1, 9], 10)
    assert res.bin_of_gene.dtype == np.int64
    assert sorted(res.bin_of_gene.tolist()) == [0, 1, 2]
    assert res.bin_of_gene.tolist() == [1, 0, 2]  # ascending a_g


def test_equal_count_rank_bins() -> None:
    sizes = np.array([5, 1, 9, 9, 9, 9, 2, 3, 100, 1, 4, 9, 9])
    res = mayu_gene_fdr(np.full(13, 0.01), [0.01], None, sizes, 4)
    b = res.bin_of_gene
    assert b.dtype == np.int64 and b.shape == (13,)
    counts = np.bincount(b)
    assert len(counts) == 4 and counts.max() - counts.min() <= 1
    order = np.argsort(sizes, kind="stable")  # ties broken by input order
    assert (np.diff(b[order]) >= 0).all()


def test_ties_split_by_input_order() -> None:
    """All sizes equal: bins follow input order (method='first')."""
    res = mayu_gene_fdr(np.full(6, 0.01), [0.01], None, np.full(6, 4), 3)
    assert res.bin_of_gene.tolist() == [0, 0, 1, 1, 2, 2]


# ----------------------------------------------------------------------------- validation


def test_partly_nan_pep_raises() -> None:
    with pytest.raises(ValueError):
        mayu_gene_fdr([0.01, NAN], [0.01, 0.01], [0.1, NAN], [1, 1], 1)


@pytest.mark.parametrize(
    ("best_q", "peptide_q", "pep", "sizes"),
    [
        ([0.01, NAN], [0.01], None, [1, 1, 1]),  # sizes length != best_q length
        ([0.01, NAN], [0.01, 0.02], [0.1], [1, 1]),  # pep length != peptide_q length
    ],
)
def test_length_mismatch_raises(
    best_q: list[float], peptide_q: list[float], pep: list[float] | None, sizes: list[int]
) -> None:
    with pytest.raises(ValueError):
        mayu_gene_fdr(best_q, peptide_q, pep, sizes, 1)


@pytest.mark.parametrize("n_bins", [0, -1])
def test_bad_n_bins_raises(n_bins: int) -> None:
    with pytest.raises(ValueError):
        mayu_gene_fdr([0.01], [0.01], None, [1], n_bins)


# ----------------------------------------------------------------------------- determinism


def test_deterministic_and_inputs_not_mutated() -> None:
    rng = np.random.default_rng(3)
    best_q = np.where(rng.random(500) < 0.4, np.round(rng.uniform(0, 0.05, 500), 3), NAN)
    peptide_q = np.round(rng.uniform(0, 0.05, 800), 3)
    pep = rng.uniform(0, 0.2, 800)
    sizes = rng.integers(0, 80, 500).astype(float)
    copies = [x.copy() for x in (best_q, peptide_q, pep, sizes)]
    r1 = mayu_gene_fdr(best_q, peptide_q, pep, sizes, 10)
    r2 = mayu_gene_fdr(best_q, peptide_q, pep, sizes, 10)
    pd.testing.assert_frame_equal(r1.curve, r2.curve)
    np.testing.assert_array_equal(r1.gene_q, r2.gene_q)
    np.testing.assert_array_equal(r1.bin_of_gene, r2.bin_of_gene)
    for orig, now in zip(copies, (best_q, peptide_q, pep, sizes), strict=True):
        np.testing.assert_array_equal(orig, now)
    assert r1.gene_q.dtype == np.float64


def test_result_is_frozen() -> None:
    res = mayu_gene_fdr([0.01], [0.01], None, [1], 1)
    with pytest.raises(AttributeError):
        res.gene_q = np.zeros(1)  # type: ignore[misc]


# ----------------------------------------------------------------------------- Monte Carlo model


def test_monte_carlo_model_cross_check() -> None:
    """Drop F false peptides ∝ a_g (2,000 seeded draws); newly hit unhit genes ≈ Σ (N_b - n_b) p_b.

    Bins align with the two size groups so the binned w_b equals each gene's exact share.
    """
    sizes = np.r_[np.full(600, 3.0), np.full(400, 12.0)]
    n_genes = sizes.size
    rng_hits = np.random.default_rng(11)
    detected = rng_hits.random(n_genes) < 0.5
    best_q = np.where(detected, 0.01, NAN)
    F = 30
    peptide_q = np.full(F, 0.01)
    res = mayu_gene_fdr(best_q, peptide_q, np.ones(F), sizes, 10)
    fp_model = _row(res, 0.01)["fp"]

    # model expectation from the bins the function used
    a = sizes
    w = np.array([a[res.bin_of_gene == b].mean() for b in range(10)]) / a.sum()
    p_b = 1.0 - (1.0 - w) ** F
    unhit = ~detected
    N_minus_n = np.array([(unhit & (res.bin_of_gene == b)).sum() for b in range(10)])
    expected_new = float((N_minus_n * p_b).sum())
    # closed form is the fixed point FP_b = (N_b - n_b + FP_b) p_b
    assert fp_model == pytest.approx(float((N_minus_n * p_b / (1 - p_b)).sum()), rel=1e-9)
    # the Poisson-binomial mean of per-gene hit indicators is the same expectation
    p_gene = p_b[res.bin_of_gene][unhit]
    pmf = poisson_binomial_pmf(p_gene)
    assert float((np.arange(pmf.size) * pmf).sum()) == pytest.approx(expected_new, rel=1e-9)

    rng = np.random.default_rng(20261014)
    draws = 2000
    new_hits = np.empty(draws)
    for i in range(draws):
        g = rng.choice(n_genes, size=F, replace=True, p=a / a.sum())
        hit = np.zeros(n_genes, dtype=bool)
        hit[g] = True
        new_hits[i] = (hit & unhit).sum()
    assert new_hits.mean() == pytest.approx(expected_new, rel=0.10)
