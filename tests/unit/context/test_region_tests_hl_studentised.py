"""M4-17 lead decision 11 / ADR-0018 amendment 1: the TU-cluster-studentised HL shift test.

docs/adr/ADR-0018-region-tests-confounding.md "Amendment 1" (accepted-by-khoa 2026-10-16),
docs/specs/M4-17-region-tests.md decision 11, formula in docs/methods/region-tests.md "Review round
2" item 3 (statistician's scripts local/scratch/review-m4-17/v2/exp5.py):

- the circular-shift p-value of the HL row comes from T = HL / sqrt(V), recomputed at every shift;
- V = sigma_u^2 (sum_t m_t^2 / m^2 + sum_t m'_t^2 / m'^2) + sigma_e^2 (1/m + 1/m'), m_t / m'_t the
  detected genes with a finite residual in TU t inside / outside the (shifted) set, m / m' their
  totals;
- sigma_u^2 (between-TU, truncated at 0) and sigma_e^2 (within-TU) are one-way ANOVA moment
  estimates on the replicon's residualised detected genes, computed once per replicon;
- the reported HL effect and its CI are unchanged (residual scale, decision 10).

Public interface pinned here (test-author choice, ``giap.context.regions``):

- ``tu_variance_components(a, ok, tu) -> (sigma_u2, sigma_e2)``;
- ``hl_variance(ok, tu, idx, shifts, sigma_u2, sigma_e2) -> V per shift`` (NaN if a side is empty);
- ``studentised_hl_null(a, ok, tu, idx, shifts) -> T per shift`` (shift 0 = the observed T).

Test-author readings of the edge cases: one TU -> sigma_u^2 = 0 (no between-TU information; the
reviewer's script divides by n0 = 0); all TUs singletons -> sigma_e^2 = 0 (SSE / max(M - n_T, 1));
V = 0 only when every value is equal, so HL = 0 and T = 0 (p = 1, as the plain HL gives now);
no detected gene -> NA, never an exception.
"""

from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import region_tests_data as rd
from region_tests_data import Contig, gene_table, pick, region_frames, tags_of


def _regions() -> Any:
    from giap.context import regions as R

    return R


def _run(contigs: list[Contig], members: dict[str, tuple[str, list[str]]], seed: int = 1,
         cfg: Any = None) -> pd.DataFrame:  # fmt: skip
    genes, replicons = gene_table(contigs)
    genes, regions, rg = region_frames(genes, members)
    cfg = cfg if cfg is not None else rd.config()
    return _regions().test_all_regions(
        genes, regions, rg, replicons, cfg, np.random.default_rng(seed)
    )


# --------------------------------------------------------------------------- hand example
# Seven genes in genome order; TU t0 = {0, 2}, t1 = {1} (one detected gene), t2 = {3, 4, 5, 6}
# with gene 4 undetected. Detected values 1, 2, 3, 4, 6, 8 (M = 6, n_T = 3).
#   TU means 2, 2, 6; grand mean 4; SSE = 2 + 0 + 8 = 10 -> sigma_e^2 = 10 / (6 - 3) = 10/3.
#   SSB = 2 * 4 + 1 * 4 + 3 * 4 = 24 -> MSB = 12; n0 = (6 - 14/6) / 2 = 11/6;
#   sigma_u^2 = (12 - 10/3) / (11/6) = 52/11.
# Region I = {0, 1}:
#   k = 0: in t0:1, t1:1 (m 2, sum m_t^2 2); out t0:1, t2:3 (m' 4, 10)
#          V = 52/11 (2/4 + 10/16) + 10/3 (1/2 + 1/4) = 86/11;  HL({1,2} - {3,4,6,8}) = -3.5
#   k = 2: in t0:1, t2:1 (2, 2); out t0:1, t1:1, t2:2 (4, 6)
#          V = 52/11 (2/4 + 6/16) + 10/3 (3/4) = 73/11;     HL({3,4} - {1,2,6,8}) = -0.5
#   k = 3: in {3, 4} -> gene 4 undetected: t2:1 (1, 1); out t0:2, t1:1, t2:2 (5, 9)
#          V = 52/11 (1 + 9/25) + 10/3 (1 + 1/5) = 2868/275; HL({4} - {1,2,3,6,8}) = 1
#   k = 5: in {5, 6}: t2:2 (2, 4); out t0:2, t1:1, t2:1 (4, 6)
#          V = 52/11 (4/4 + 6/16) + 10/3 (1/2 + 1/4) = 9;     HL({6,8} - {1,2,3,4}) = 4.5

A = np.array([1.0, 2.0, 3.0, 4.0, np.nan, 6.0, 8.0])
OK = np.isfinite(A)
TU = np.array(["t0", "t1", "t0", "t2", "t2", "t2", "t2"], dtype=object)
IDX = np.array([0, 1])
SU2, SE2 = Fraction(52, 11), Fraction(10, 3)
HAND = {  # shift -> (V, HL)
    0: (Fraction(86, 11), -3.5),
    2: (Fraction(73, 11), -0.5),
    3: (Fraction(2868, 275), 1.0),
    5: (Fraction(9), 4.5),
}


def test_hand_example_is_consistent_with_oracle() -> None:
    """The hand arithmetic above agrees with the naive oracle (guards the oracle itself)."""
    su2, se2 = rd.naive_tu_components(A, OK, list(TU))
    assert su2 == pytest.approx(float(SU2)) and se2 == pytest.approx(float(SE2))
    ks = list(HAND)
    v = rd.naive_hl_variance(OK, list(TU), list(IDX), ks, su2, se2)
    np.testing.assert_allclose(v, [float(HAND[k][0]) for k in ks], rtol=1e-12)
    np.testing.assert_allclose(rd.naive_hl_null(A, OK, list(IDX), ks), [HAND[k][1] for k in ks])


def test_tu_variance_components_hand_example() -> None:
    su2, se2 = _regions().tu_variance_components(A, OK, TU)
    assert su2 == pytest.approx(float(SU2), rel=1e-12)
    assert se2 == pytest.approx(float(SE2), rel=1e-12)


def test_tu_variance_components_ignore_undetected_values() -> None:
    """A value at a gene with ok = False (here 1000 at the undetected gene 4) does not enter."""
    a = A.copy()
    a[4] = 1000.0
    su2, se2 = _regions().tu_variance_components(a, OK, TU)
    assert su2 == pytest.approx(float(SU2), rel=1e-12)
    assert se2 == pytest.approx(float(SE2), rel=1e-12)


def test_hl_variance_hand_example() -> None:
    ks = np.array(list(HAND))
    v = _regions().hl_variance(OK, TU, IDX, ks, float(SU2), float(SE2))
    np.testing.assert_allclose(v, [float(HAND[k][0]) for k in HAND], rtol=1e-12)


def test_studentised_hl_null_hand_example() -> None:
    ks = np.array(list(HAND))
    t = _regions().studentised_hl_null(A, OK, TU, IDX, ks)
    want = [h / math.sqrt(float(v)) for v, h in HAND.values()]
    np.testing.assert_allclose(t, want, rtol=1e-12)


def test_studentised_hl_null_equals_oracle_at_every_shift() -> None:
    ks = np.arange(7)
    t = _regions().studentised_hl_null(A, OK, TU, IDX, ks)
    want = rd.naive_studentised_hl_null(A, OK, list(TU), list(IDX), list(ks))
    np.testing.assert_allclose(t, want, rtol=1e-12, equal_nan=True)


def test_sigmas_are_estimated_once_per_replicon_not_per_shift() -> None:
    """V_k uses the replicon's sigma's at every shift: V / (su2-part, se2-part) is linear in the
    supplied sigma's (doubling both doubles V) and T does not depend on the region."""
    R = _regions()
    ks = np.array(list(HAND))
    v1 = R.hl_variance(OK, TU, IDX, ks, 1.0, 0.0)
    v2 = R.hl_variance(OK, TU, IDX, ks, 0.0, 1.0)
    v = R.hl_variance(OK, TU, IDX, ks, float(SU2), float(SE2))
    np.testing.assert_allclose(v, float(SU2) * v1 + float(SE2) * v2, rtol=1e-12)
    np.testing.assert_allclose(v1, [2 / 4 + 10 / 16, 2 / 4 + 6 / 16, 1 + 9 / 25, 4 / 4 + 6 / 16])
    np.testing.assert_allclose(v2, [3 / 4, 3 / 4, 6 / 5, 3 / 4])


# --------------------------------------------------------------------------- edge cases


def test_one_tu_has_no_between_tu_variance() -> None:
    """All detected genes in one TU: sigma_u^2 = 0 (finite, not NaN), sigma_e^2 = sample
    variance."""
    a = np.array([1.0, 2.0, 4.0, 7.0])
    su2, se2 = _regions().tu_variance_components(a, np.ones(4, bool), np.array(["t"] * 4))
    assert su2 == 0.0
    assert se2 == pytest.approx(float(np.var(a, ddof=1)))
    t = _regions().studentised_hl_null(a, np.ones(4, bool), np.array(["t"] * 4), np.array([0]),
                                       np.arange(4))  # fmt: skip
    assert np.all(np.isfinite(t))


def test_singleton_tus_put_everything_between() -> None:
    """Every TU one detected gene (the hand-built default): sigma_e^2 = 0, n0 = 1, sigma_u^2 =
    sample variance; V = sigma_u^2 (1/m + 1/m')."""
    a = np.array([3.0, 1.0, 4.0, 1.0, 5.0, 9.0])
    tu = np.array([f"t{i}" for i in range(6)])
    R = _regions()
    su2, se2 = R.tu_variance_components(a, np.ones(6, bool), tu)
    assert se2 == 0.0
    assert su2 == pytest.approx(float(np.var(a, ddof=1)))
    v = R.hl_variance(np.ones(6, bool), tu, np.array([0, 1]), np.array([0, 3]), su2, se2)
    np.testing.assert_allclose(v, su2 * (1 / 2 + 1 / 4))


def test_between_tu_variance_truncated_at_zero() -> None:
    """TU means closer than within-TU noise allows: MSB < sigma_e^2 -> sigma_u^2 = 0, not < 0."""
    a = np.array([0.0, 10.0, 1.0, 9.0, 2.0, 8.0])
    tu = np.array(["a", "a", "b", "b", "c", "c"])
    su2, se2 = _regions().tu_variance_components(a, np.ones(6, bool), tu)
    assert su2 == 0.0
    assert se2 == pytest.approx((50 + 32 + 18) / 3)


def test_v_zero_gives_t_zero() -> None:
    """All detected values equal: sigma's are 0, V = 0 and HL = 0 at every shift -> T = 0 (never
    NaN, never a ZeroDivisionError)."""
    a = np.full(8, 5.0)
    tu = np.array(["a", "a", "b", "b", "c", "c", "d", "d"])
    R = _regions()
    assert R.tu_variance_components(a, np.ones(8, bool), tu) == (0.0, 0.0)
    t = R.studentised_hl_null(a, np.ones(8, bool), tu, np.array([1, 2]), np.arange(8))
    np.testing.assert_array_equal(t, np.zeros(8))


def test_empty_side_gives_nan_t() -> None:
    """A shift whose set holds no detected gene (or every detected gene) has T = NaN."""
    a = np.array([1.0, np.nan, np.nan, 4.0, 6.0])
    ok = np.isfinite(a)
    tu = np.array(["a", "a", "b", "c", "c"])
    t = _regions().studentised_hl_null(a, ok, tu, np.array([0]), np.array([0, 1, 2, 3]))
    assert math.isnan(t[1]) and math.isnan(t[2])
    assert np.isfinite(t[0]) and np.isfinite(t[3])


def test_no_detected_gene_does_not_raise() -> None:
    a = np.full(5, np.nan)
    ok = np.zeros(5, bool)
    tu = np.array(["a", "a", "b", "c", "c"])
    R = _regions()
    R.tu_variance_components(a, ok, tu)  # no exception; values unused
    t = R.studentised_hl_null(a, ok, tu, np.array([0, 1]), np.arange(5))
    assert t.shape == (5,) and np.all(np.isnan(t))


# --------------------------------------------------------------------------- test_region rows


def _tu_contig(seed: int = 31, n: int = 48) -> tuple[Contig, list[int]]:
    """TUs of sizes 1, 2, 3, 4, ... in turn (single-gene TUs included), a TU random effect on
    abundance, ~40 % undetected, random detectability covariates; region = genes 12 ... 21."""
    rng = np.random.default_rng(seed)
    tus: list[str] = []
    size = 1
    while len(tus) < n:
        tus += [f"TU{len(tus)}"] * size
        size = size % 4 + 1
    tus = tus[:n]
    u = {t: rng.normal(0, 1.0) for t in sorted(set(tus))}
    d = (rng.random(n) < 0.6).astype(int).tolist()
    length = rng.integers(40, 800, n)
    n_obs = np.maximum(1, length // 25 + rng.integers(-3, 4, n))
    gravy = rng.normal(-0.3, 0.4, n)
    a = 18 + 1.5 * np.log(length) - 0.8 * gravy + np.array([u[t] for t in tus])
    a = a + rng.normal(0, 0.7, n)
    c = Contig("chr1", d, rng.uniform(0.3, 0.9, n).tolist(), a=a.tolist(), tu=tus,
               length_aa=length.tolist(), n_obs=n_obs.tolist(), gravy=gravy.tolist())  # fmt: skip
    return c, list(range(12, 22))


def _oracle(c: Contig, idx: list[int]) -> dict[str, Any]:
    det = np.array(c.d, bool)
    a = np.where(det, np.array(c.a, float), np.nan)
    n = len(c.d)
    L0, N0, G0 = rd.DEFAULT_COVARIATES
    ar = rd.naive_residualise(a, det, c.length_aa or [L0] * n, c.n_obs or [N0] * n,
                              c.gravy or [G0] * n)  # fmt: skip
    ok = det & np.isfinite(ar)
    tu = list(c.tu or [f"TU_{c.name}_{i}" for i in range(len(c.d))])
    shifts = list(range(1, len(ar)))
    x = [ar[i] for i in idx if ok[i]]
    y = [ar[i] for i in range(len(ar)) if ok[i] and i not in idx]
    t_obs = float(rd.naive_studentised_hl_null(ar, ok, tu, idx, [0])[0])
    t_null = rd.naive_studentised_hl_null(ar, ok, tu, idx, shifts)
    hl_null = rd.naive_hl_null(ar, ok, idx, shifts)
    eff = rd.naive_hl(x, y)
    return {
        "effect": eff,
        "p": rd.naive_perm_p(t_obs, t_null),
        "p_plain": rd.naive_perm_p(eff, hl_null),
        "k": int(np.count_nonzero(~np.isnan(t_null))),
    }


def test_hl_p_value_from_studentised_statistic() -> None:
    """HL row p = two-sided centred perm p of T over the circular shifts (oracle), which differs
    here from the plain-HL p of decision 10."""
    c, idx = _tu_contig()
    genes, _ = gene_table([c])
    row = pick(_run([c], {"prophage_1": ("prophage", tags_of(genes, "chr1", idx))}),
               "prophage_1", "hl_shift")  # fmt: skip
    o = _oracle(c, idx)
    assert o["p"] != pytest.approx(o["p_plain"]), "fixture does not separate T from plain HL"
    assert row["p_value"] == pytest.approx(o["p"])
    assert row["n_perm"] == o["k"]
    assert row["effect"] == pytest.approx(o["effect"], abs=1e-9)


# Reported HL effect and CI under decision 10 (current residual-scale implementation, commit
# 529bd14, seed 1, rd.config()): amendment 1 changes the p-value only.
GOLDEN_EFFECT = 0.07878994901260405
GOLDEN_CI = (-0.6003296150369698, 0.5595939753194148)


def test_hl_effect_and_ci_unchanged_by_amendment() -> None:
    c, idx = _tu_contig()
    genes, _ = gene_table([c])
    row = pick(_run([c], {"prophage_1": ("prophage", tags_of(genes, "chr1", idx))}),
               "prophage_1", "hl_shift")  # fmt: skip
    assert row["effect"] == pytest.approx(GOLDEN_EFFECT, abs=1e-12)
    assert row["ci_low"] == pytest.approx(GOLDEN_CI[0], abs=1e-12)
    assert row["ci_high"] == pytest.approx(GOLDEN_CI[1], abs=1e-12)


def test_constant_n_detected_gives_plain_hl_p() -> None:
    """Every gene detected and every TU a singleton: m and sum m_t^2 are the same at every shift,
    so V is constant and T orders the shifts as HL does: p equals the plain-HL p."""
    rng = np.random.default_rng(5)
    n = 40
    a = rng.normal(10, 1, n).tolist()
    c = Contig("chr1", [1] * n, [0.5] * n, a=a)
    genes, _ = gene_table([c])
    idx = list(range(4, 11))
    row = pick(_run([c], {"bgc_1": ("bgc", tags_of(genes, "chr1", idx))}), "bgc_1", "hl_shift")
    av = np.array(a)
    av = av - av.mean()
    null = rd.naive_hl_null(av, np.ones(n, bool), idx, range(1, n))
    eff = rd.naive_hl([av[i] for i in idx], [av[i] for i in range(n) if i not in idx])
    assert row["p_value"] == pytest.approx(rd.naive_perm_p(eff, null))


def test_constant_abundance_gives_p_one() -> None:
    """Constant abundance (V = 0 up to rounding of the residual fit): effect 0 and p = 1, as the
    plain HL gives; T must not turn rounding noise into signal."""
    n = 30
    tus = [f"TU{i // 3}" for i in range(n)]
    d = [0 if i % 4 == 0 else 1 for i in range(n)]
    c = Contig("chr1", d, [0.5] * n, a=[12.3] * n, tu=tus)
    genes, _ = gene_table([c])
    row = pick(_run([c], {"bgc_1": ("bgc", tags_of(genes, "chr1", range(6, 12)))}),
               "bgc_1", "hl_shift")  # fmt: skip
    assert abs(row["effect"]) < 1e-9
    assert row["p_value"] == pytest.approx(1.0)


def test_one_tu_replicon_runs() -> None:
    """Whole replicon one TU (sigma_u^2 = 0): the HL row has a finite p in (0, 1]."""
    c, idx = _tu_contig(33)
    c = Contig(c.name, c.d, c.p, a=c.a, tu=["TU_all"] * len(c.d), length_aa=c.length_aa,
               n_obs=c.n_obs, gravy=c.gravy)  # fmt: skip
    genes, _ = gene_table([c])
    row = pick(_run([c], {"bgc_1": ("bgc", tags_of(genes, "chr1", idx))}), "bgc_1", "hl_shift")
    o = _oracle(c, idx)
    assert 0.0 < row["p_value"] <= 1.0
    assert row["p_value"] == pytest.approx(o["p"])


def test_region_with_one_detected_gene() -> None:
    """A region whose only detected gene sits in a multi-gene TU: m = 1 at the observed position;
    p equals the oracle and n_perm counts the shifts with a detected gene inside."""
    n = 24
    tus = [f"TU{i // 3}" for i in range(n)]
    d = [1] * n
    d[7], d[8] = 0, 0  # region {6, 7, 8} = TU2, only gene 6 detected
    rng = np.random.default_rng(9)
    a = (10 + rng.normal(0, 1, n)).tolist()
    c = Contig("chr1", d, [0.5] * n, a=a, tu=tus)
    genes, _ = gene_table([c])
    idx = [6, 7, 8]
    row = pick(_run([c], {"crispr_1": ("crispr", tags_of(genes, "chr1", idx))}),
               "crispr_1", "hl_shift")  # fmt: skip
    o = _oracle(c, idx)
    assert row["p_value"] == pytest.approx(o["p"])
    assert row["n_perm"] == o["k"]


def test_all_genes_undetected_gives_na_row() -> None:
    """No detected gene on the replicon: HL row effect / p NA, no exception; O - E row normal."""
    n = 20
    c = Contig("chr1", [0] * n, [0.4] * n, tu=[f"TU{i // 2}" for i in range(n)])
    genes, _ = gene_table([c])
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "chr1", range(3, 7)))})
    hl = pick(res, "bgc_1", "hl_shift")
    assert rd.is_na(hl["effect"]) and rd.is_na(hl["p_value"]) and hl["n_perm"] == 0
    assert 0.0 < pick(res, "bgc_1", "o_minus_e")["p_value"] <= 1.0


def test_all_region_genes_detected_but_nothing_outside() -> None:
    """Only the region's genes are detected: the outside side is empty -> HL NA, no exception."""
    n = 20
    d = [1 if 3 <= i < 7 else 0 for i in range(n)]
    c = Contig("chr1", d, [0.4] * n, tu=[f"TU{i // 2}" for i in range(n)])
    genes, _ = gene_table([c])
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "chr1", range(3, 7)))})
    hl = pick(res, "bgc_1", "hl_shift")
    assert rd.is_na(hl["effect"]) and rd.is_na(hl["p_value"])


def test_hl_row_invariant_to_tu_renaming() -> None:
    """TU labels are only labels: renaming them (a bijection) leaves the HL row unchanged."""
    c, idx = _tu_contig(35)
    ren = Contig(c.name, c.d, c.p, a=c.a, tu=[f"z{t[::-1]}" for t in (c.tu or [])],
                 length_aa=c.length_aa, n_obs=c.n_obs, gravy=c.gravy)  # fmt: skip
    genes, _ = gene_table([c])
    members = {"bgc_1": ("bgc", tags_of(genes, "chr1", idx))}
    a = pick(_run([c], members), "bgc_1", "hl_shift")
    b = pick(_run([ren], members), "bgc_1", "hl_shift")
    for col in ("effect", "p_value", "ci_low", "ci_high", "n_perm"):
        assert a[col] == pytest.approx(b[col], abs=1e-12)
