"""M4-17 lead decisions 13 and 15 / ADR-0018 amendment 2: the label-free studentised HL shift test.

docs/adr/ADR-0018-region-tests-confounding.md "Amendment 2" (accepted-by-khoa 2026-10-16),
docs/specs/M4-17-region-tests.md decisions 13 and 15, formula in docs/methods/region-tests.md
"Review round 3" item 4 (statistician's prototype local/scratch/review-m4-17/v3/lagv.py):

- the circular-shift p-value of the HL row comes from T = HL / sqrt(V), recomputed at every shift
  (as amendment 1); the reported HL effect and its CI are unchanged (residual scale, decision 10);
- V = Var(mean_in - mean_out) = S_in,in / m^2 + S_out,out / m'^2 - 2 S_in,out / (m m'), S_A,B the
  sum of gamma(d(i, j)) over detected genes i in A, j in B with a finite residual, d the circular
  gene-index distance, m / m' the detected genes inside / outside the (shifted) set;
- gamma(h) is the replicon's autocovariance of the residualised detected abundances at circular
  gene-index lag h <= 30, estimated once per replicon and truncated at its first non-positive
  value (gamma(h') = 0 from there on; oracle ``rd.naive_lag_autocov``);
- V uses no TU label: relabelling, fragmenting or merging TUs leaves V, T and p unchanged
  (supersedes amendment 1's per-TU components, decision 12's ``tu_variance_components``);
- decision 12's guard stays, scaled to the data: V counts as 0 (T = 0) when gamma(0) <=
  (V_RTOL * scale)^2, ``scale`` the magnitude of the data before residualisation (default max |a|);
- decision 13: a missing ``tu_id`` raises ValueError (C6 makes ``tu_id`` non-null; no per-gene
  fallback). The bootstrap CI still resamples TU blocks.

Public interface pinned here (test-author choice, ``giap.context.regions``), replacing
``tu_variance_components`` and the TU arguments of ``hl_variance`` / ``studentised_hl_null``:

- ``LAG_MAX = 30``;
- ``lag_autocovariance(a, ok, max_lag=LAG_MAX) -> gamma[0 .. max_lag]`` (truncated);
- ``hl_variance(ok, idx, shifts, gamma) -> V per shift`` (NaN if a side is empty);
- ``studentised_hl_null(a, ok, idx, shifts, *, scale=None) -> T per shift`` (shift 0 = observed).

Not pinned: V <= 0 with gamma(0) above the guard (the truncated gamma need not be positive
definite); oracle comparisons skip such inputs (``rd.naive_v_positive``).
"""

from __future__ import annotations

import inspect
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


def _pad(g: list[Any], n: int = 31) -> np.ndarray:
    return np.array([float(x) for x in g] + [0.0] * (n - len(g)))


# --------------------------------------------------------------------------- hand example 1
# Eight genes on a circular replicon, gene 4 undetected: a = 1, 3, 2, 4, -, 6, 5, 7 (M = 7).
# Mean 4 -> x = -3, -1, -2, 0, (0), 2, 1, 3.
#   gamma(0) = (9 + 1 + 4 + 0 + 4 + 1 + 9) / 7 = 4
#   gamma(1): products of neighbours (0,1) 3, (1,2) 2, (2,3) 0, (5,6) 2, (6,7) 3, (7,0) -9 -> 1;
#             6 pairs with both detected -> 1/6
#   gamma(2): (0,2) 6, (1,3) 0, (3,5) 0, (5,7) 6, (6,0) -3, (7,1) -3 -> 6 / 6 = 1
#   gamma(3): (0,3) 0, (2,5) -4, (3,6) 0, (5,0) -6, (6,1) -1, (7,2) -6 -> -17/6 <= 0: truncated,
#             gamma(h) = 0 for every h >= 3 (gamma(4) = -8/6 would be negative too).
# Region I = {0, 1}:
#   k = 0: in {0, 1} (m 2), out {2, 3, 5, 6, 7} (m' 5).
#          S_in,in = 2 gamma0 + 2 gamma1 = 25/3; S_out,out = 5 gamma0 + 2 (3 gamma1 + 2 gamma2) = 25
#          (distance-1 pairs 23, 56, 67; distance-2 pairs 35, 57); S_in,out = 2 gamma1 + 4 gamma2
#          = 13/3 (distance 1: 12, 07; distance 2: 02, 06, 13, 17).
#          V = 25/12 + 25/25 - 2 (13/3) / 10 = 133/60;  HL({1, 3} - {2, 4, 6, 5, 7}) = -3
#   k = 3: in {3} (gene 4 undetected): HL({4} - {1, 3, 2, 6, 5, 7}) = 0, V = 455/108, T = 0
#   k = 5: in {5, 6}: HL({6, 5} - {1, 3, 2, 4, 7}) = 2.5, V = 763/300

A = np.array([1.0, 3.0, 2.0, 4.0, np.nan, 6.0, 5.0, 7.0])
OK = np.isfinite(A)
IDX = np.array([0, 1])
GAMMA = [Fraction(4), Fraction(1, 6), Fraction(1)]  # gamma(h) = 0 for h >= 3
V_HAND = [
    Fraction(133, 60),
    Fraction(749, 300),
    Fraction(763, 300),
    Fraction(455, 108),
    Fraction(455, 108),
    Fraction(763, 300),
    Fraction(749, 300),
    Fraction(133, 60),
]
HAND = {0: (V_HAND[0], -3.0), 3: (V_HAND[3], 0.0), 5: (V_HAND[5], 2.5)}  # shift -> (V, HL)

# Hand example 2 (truncation, not clipping): ten genes, all detected, a = 6 0 4 2 1 0 1 1 0 5,
# mean 2, x = 4 -2 2 0 -1 -2 -1 -1 -2 3. gamma(0) = 44/10 = 22/5; gamma(1) = (-8 - 4 + 0 + 0 + 2 + 2
# + 1 + 2 - 6 + 12) / 10 = 1/10; gamma(2) = -6/10 <= 0, so gamma(2 ...) = 0 although the raw
# gamma(3) = 3/5 is positive.
A2 = np.array([6.0, 0.0, 4.0, 2.0, 1.0, 0.0, 1.0, 1.0, 0.0, 5.0])
GAMMA2 = [Fraction(22, 5), Fraction(1, 10)]


def test_hand_examples_agree_with_oracle() -> None:
    """The hand arithmetic above agrees with the exact (Fraction) oracle (guards the oracle)."""
    af = [Fraction(int(x)) if np.isfinite(x) else math.nan for x in A]
    g = rd.naive_lag_autocov(af, OK)
    assert g[:3] == GAMMA and all(x == 0 for x in g[3:]) and len(g) == rd.LAG_MAX + 1
    assert rd.naive_lag_hl_variance(OK, list(IDX), range(8), g) == V_HAND
    np.testing.assert_allclose(rd.naive_hl_null(A, OK, list(IDX), list(HAND)),
                               [h for _, h in HAND.values()])  # fmt: skip
    g2 = rd.naive_lag_autocov([Fraction(int(x)) for x in A2], np.ones(10, bool))
    assert g2[:2] == GAMMA2 and all(x == 0 for x in g2[2:])


def test_lag_max_is_30() -> None:
    assert _regions().LAG_MAX == 30


def test_lag_autocovariance_hand_example() -> None:
    g = _regions().lag_autocovariance(A, OK)
    assert g.shape == (31,)
    np.testing.assert_allclose(g, _pad(GAMMA), rtol=1e-12, atol=1e-15)


def test_lag_autocovariance_truncates_at_first_non_positive_lag() -> None:
    g = _regions().lag_autocovariance(A2, np.ones(10, bool))
    np.testing.assert_allclose(g, _pad(GAMMA2), rtol=1e-12, atol=1e-15)


def test_lag_autocovariance_max_lag() -> None:
    """max_lag bounds the lags: gamma[0 .. max_lag]."""
    g = _regions().lag_autocovariance(A, OK, max_lag=1)
    np.testing.assert_allclose(g, [4.0, 1 / 6], rtol=1e-12)


def test_lag_autocovariance_ignores_undetected_values() -> None:
    """A value at a gene with ok = False (1000 at gene 4) enters neither the mean nor any lag."""
    a = A.copy()
    a[4] = 1000.0
    np.testing.assert_allclose(_regions().lag_autocovariance(a, OK), _pad(GAMMA), rtol=1e-12,
                               atol=1e-15)  # fmt: skip


def test_lag_autocovariance_is_circular() -> None:
    """Lags wrap around the origin: rotating the replicon leaves gamma unchanged."""
    R = _regions()
    np.testing.assert_allclose(R.lag_autocovariance(np.roll(A, 3), np.roll(OK, 3)),
                               R.lag_autocovariance(A, OK), rtol=1e-12, atol=1e-15)  # fmt: skip


def test_hl_variance_hand_example() -> None:
    v = _regions().hl_variance(OK, IDX, np.arange(8), _pad(GAMMA))
    np.testing.assert_allclose(v, [float(x) for x in V_HAND], rtol=1e-12)


def test_hl_variance_white_noise_is_two_sample_variance() -> None:
    """gamma = (g0, 0, 0, ...): V = g0 (1/m + 1/m') (k = 0: 4 (1/2 + 1/5) = 14/5; k = 3: m = 1,
    m' = 6: 4 (1 + 1/6) = 14/3)."""
    v = _regions().hl_variance(OK, IDX, np.array([0, 3]), _pad([4]))
    np.testing.assert_allclose(v, [14 / 5, 14 / 3], rtol=1e-12)


def test_hl_variance_linear_in_gamma() -> None:
    R = _regions()
    ks = np.arange(8)
    g1, g2 = _pad([1, 0.5, 0.25]), _pad([2, 0, 0, 1])
    np.testing.assert_allclose(R.hl_variance(OK, IDX, ks, 3 * g1 + g2),
                               3 * R.hl_variance(OK, IDX, ks, g1) + R.hl_variance(OK, IDX, ks, g2),
                               rtol=1e-12)  # fmt: skip


def test_studentised_hl_null_hand_example() -> None:
    ks = np.array(list(HAND))
    t = _regions().studentised_hl_null(A, OK, IDX, ks)
    np.testing.assert_allclose(t, [h / math.sqrt(float(v)) for v, h in HAND.values()], rtol=1e-12)


def test_studentised_hl_null_equals_oracle_at_every_shift() -> None:
    ks = np.arange(8)
    t = _regions().studentised_hl_null(A, OK, IDX, ks)
    want = rd.naive_studentised_hl_null(A, OK, list(IDX), list(ks))
    np.testing.assert_allclose(t, want, rtol=1e-12, equal_nan=True)


def test_gamma_is_estimated_once_per_replicon_not_per_shift() -> None:
    """T_k = HL_k / sqrt(V_k) with V_k from the replicon's one gamma (not re-estimated per set)."""
    R = _regions()
    ks = np.arange(8)
    g = R.lag_autocovariance(A, OK)
    want = R.hl_null(A, OK, IDX, ks) / np.sqrt(R.hl_variance(OK, IDX, ks, g))
    np.testing.assert_allclose(R.studentised_hl_null(A, OK, IDX, ks), want, rtol=1e-12)


def test_interface_takes_no_tu_labels() -> None:
    """Decision 15: V is label-free; the variance functions take no TU argument."""
    R = _regions()
    assert list(inspect.signature(R.hl_variance).parameters) == ["ok", "idx", "shifts", "gamma"]
    p = inspect.signature(R.studentised_hl_null).parameters
    assert list(p)[:4] == ["a", "ok", "idx", "shifts"] and "tu" not in p
    assert p["scale"].kind is inspect.Parameter.KEYWORD_ONLY


# --------------------------------------------------------------------------- edge cases and guard


def test_v_zero_gives_t_zero() -> None:
    """All detected values equal: gamma = 0, V = 0 and HL = 0 at every shift -> T = 0 (never NaN,
    never a ZeroDivisionError)."""
    a = np.full(8, 5.0)
    R = _regions()
    np.testing.assert_array_equal(R.lag_autocovariance(a, np.ones(8, bool)), np.zeros(31))
    t = R.studentised_hl_null(a, np.ones(8, bool), np.array([1, 2]), np.arange(8))
    np.testing.assert_array_equal(t, np.zeros(8))


def test_guard_treats_rounding_noise_as_zero() -> None:
    """Values equal up to rounding (relative 1e-13): gamma(0) <= (V_RTOL max|a|)^2, so V counts as
    0 and T = 0 at every shift; noise never becomes signal."""
    rng = np.random.default_rng(3)
    a = 12.3 + 12.3 * 1e-13 * rng.standard_normal(20)
    t = _regions().studentised_hl_null(a, np.ones(20, bool), np.arange(5, 9), np.arange(20))
    np.testing.assert_array_equal(t, np.zeros(20))


def test_guard_is_relative_to_the_data_scale() -> None:
    """Tiny but genuine variation (the whole data scaled by 1e-9) is not zeroed: T is the same as
    for the unscaled data (HL scales by c, V by c^2)."""
    R = _regions()
    rng = np.random.default_rng(4)
    a = rng.normal(0, 1, 24)
    ok = np.ones(24, bool)
    ks, idx = np.arange(24), np.arange(3, 9)
    t0 = R.studentised_hl_null(a, ok, idx, ks)
    assert np.any(t0 != 0)
    np.testing.assert_allclose(R.studentised_hl_null(1e-9 * a, ok, idx, ks), t0, rtol=1e-9)


def test_guard_uses_the_supplied_scale() -> None:
    """``scale`` (magnitude before residualisation) sets the guard: residuals of size 1e-3 count as
    rounding against scale 1e5 ((1e-7 * 1e5)^2 = 1e-4 > gamma(0) ~ 1e-6) but not by default."""
    R = _regions()
    rng = np.random.default_rng(5)
    a = 1e-3 * rng.standard_normal(16)
    ok, idx, ks = np.ones(16, bool), np.arange(2, 6), np.arange(16)
    np.testing.assert_array_equal(R.studentised_hl_null(a, ok, idx, ks, scale=1e5), np.zeros(16))
    assert np.any(R.studentised_hl_null(a, ok, idx, ks) != 0)


def test_empty_side_gives_nan_t() -> None:
    """A shift whose set holds no detected gene (or every detected gene) has T = NaN."""
    a = np.array([1.0, np.nan, np.nan, 4.0, 6.0])
    ok = np.isfinite(a)
    t = _regions().studentised_hl_null(a, ok, np.array([0]), np.array([0, 1, 2, 3]))
    assert math.isnan(t[1]) and math.isnan(t[2])
    assert np.isfinite(t[0]) and np.isfinite(t[3])
    v = _regions().hl_variance(ok, np.array([0]), np.array([0, 1, 2, 3]), _pad([1.0]))
    assert math.isnan(v[1]) and math.isnan(v[2]) and np.isfinite(v[0])


def test_no_detected_gene_does_not_raise() -> None:
    a = np.full(5, np.nan)
    ok = np.zeros(5, bool)
    R = _regions()
    g = R.lag_autocovariance(a, ok)  # no exception; values unused
    assert g.shape == (31,)
    t = R.studentised_hl_null(a, ok, np.array([0, 1]), np.arange(5))
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
    shifts = list(range(1, len(ar)))
    x = [ar[i] for i in idx if ok[i]]
    y = [ar[i] for i in range(len(ar)) if ok[i] and i not in idx]
    t_obs = float(rd.naive_studentised_hl_null(ar, ok, idx, [0])[0])
    t_null = rd.naive_studentised_hl_null(ar, ok, idx, shifts)
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
# 529bd14, seed 1, rd.config()): amendments 1 and 2 change the p-value only.
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
    """Every gene detected on a circular replicon: the lag covariance matrix is circulant, so V is
    the same at every shift and T orders the shifts as HL does: p equals the plain-HL p."""
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
    """Whole replicon one TU: the HL row has a finite p in (0, 1] and equals the oracle."""
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


# --------------------------------------------------------------------------- decision 15: no labels


def _with_tus(c: Contig, tus: list[str]) -> Contig:
    return Contig(c.name, c.d, c.p, a=c.a, tu=tus, length_aa=c.length_aa, n_obs=c.n_obs,
                  gravy=c.gravy)  # fmt: skip


@pytest.mark.parametrize("variant", ["singletons", "fragment_region", "merge_pairs", "one_tu"])
def test_hl_p_unchanged_by_tu_relabelling_or_fragmenting(variant: str) -> None:
    """A-TU in miniature: TU labels do not enter V, so fragmenting the region's TUs (or all TUs)
    or merging TUs leaves the HL effect, p and n_perm unchanged (the CI may change: its bootstrap
    resamples TU blocks)."""
    c, idx = _tu_contig(37)
    tus = list(c.tu or [])
    n = len(tus)
    new = {
        "singletons": [f"S{i}" for i in range(n)],
        "fragment_region": [f"F{i}" if i in idx else t for i, t in enumerate(tus)],
        "merge_pairs": [f"M{i // 8}" for i in range(n)],
        "one_tu": ["TU_all"] * n,
    }[variant]
    genes, _ = gene_table([c])
    members = {"prophage_1": ("prophage", tags_of(genes, "chr1", idx))}
    a = pick(_run([c], members), "prophage_1", "hl_shift")
    b = pick(_run([_with_tus(c, new)], members), "prophage_1", "hl_shift")
    assert b["effect"] == pytest.approx(a["effect"], abs=1e-12)
    assert b["p_value"] == pytest.approx(a["p_value"], abs=1e-12)
    assert b["n_perm"] == a["n_perm"]
    assert a["p_value"] == pytest.approx(_oracle(c, idx)["p"])


# --------------------------------------------------------------------------- decision 13


@pytest.mark.parametrize("where", ["inside", "outside"])
def test_missing_tu_id_raises_value_error(where: str) -> None:
    """Decision 13: a gene on a tested replicon without ``tu_id`` (C6 makes it non-null) raises
    ValueError; there is no per-gene fallback. Inside or outside the region alike."""
    c, idx = _tu_contig(41)
    tus: list[Any] = list(c.tu or [])
    tus[idx[2] if where == "inside" else 40] = None
    genes, replicons = gene_table([_with_tus(c, tus)])
    assert genes["tu_id"].isna().sum() == 1
    genes, regions, rg = region_frames(genes, {"bgc_1": ("bgc", tags_of(genes, "chr1", idx))})
    with pytest.raises(ValueError, match="tu_id"):
        _regions().test_all_regions(genes, regions, rg, replicons, rd.config(),
                                    np.random.default_rng(1))  # fmt: skip


def test_missing_tu_id_raises_on_nan_label() -> None:
    """NaN (as read from a TSV with an empty cell) counts as missing too."""
    c, idx = _tu_contig(43)
    genes, replicons = gene_table([c])
    genes["tu_id"] = genes["tu_id"].astype(object)
    genes.loc[genes.index[5], "tu_id"] = np.nan
    genes, regions, rg = region_frames(genes, {"bgc_1": ("bgc", tags_of(genes, "chr1", idx))})
    with pytest.raises(ValueError, match="tu_id"):
        _regions().test_all_regions(genes, regions, rg, replicons, rd.config(),
                                    np.random.default_rng(1))  # fmt: skip
