"""Properties of the M4-17 region tests (docs/specs/M4-17-region-tests.md acceptance "Property";
GIAP_SPEC.md §7.3.6).

- P1 O - E (observed, expected, effect) is unchanged when the values of genes outside the region
  are permuted among the outside indices;
- P2 every shifted index set has |I| members (o_minus_e_null of ones is |I| at every shift);
- P3 on linear contigs no shifted set wraps: every shift used satisfies max(I) + k < N or
  min(I) + k >= N, and the count equals the spec's rule (``expected_shifts``);
- P4 the vectorised nulls equal a naive per-shift loop (O - E and HL, NaN where a side is empty);
- P5 end to end on circular contigs below the threshold: n_perm = N - 1 (O - E), and p_value
  equals the two-sided centred p of the naive null (both statistics; NaN shifts dropped). O - E:
  the null of the ADR-0018 studentised statistic sum(d - p) / sqrt(sum p (1 - p)), effect O - E;
  HL with the default constant covariates (residualising only centres a);
- P6 (ADR-0018 decision 1) with random detectability covariates, HL effect and p equal the naive
  HL on ``rd.naive_residualise`` residuals and its naive shift null;
- ADR-0018 amendments 1 and 2 (lead decisions 11, 15): in P5 / P6 the HL p comes from the
  studentised T = HL / sqrt(V), V label-free from the replicon's lag autocovariance
  (``rd.naive_studentised_hl_null``; the effect stays the HL);
- P8 (decision 15) ``lag_autocovariance`` and ``hl_variance`` equal the oracles, and the HL row's
  effect, p and n_perm do not depend on the TU partition (any relabelling / fragmenting / merging);
- P7 the HL p (and T at every shift) is invariant to adding a constant to log2 abundance and to
  scaling it by c > 0; the effect and CI scale by c.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import region_tests_data as rd
from region_tests_data import Contig, gene_table, pick, region_frames, tags_of

SETTINGS = settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])


def _regions() -> Any:
    from giap.context import regions as R

    return R


@st.composite
def replicon(draw: st.DrawFn, min_n: int = 3, max_n: int = 30) -> dict[str, Any]:
    n = draw(st.integers(min_n, max_n))
    idx = sorted(draw(st.sets(st.integers(0, n - 1), min_size=1, max_size=max(1, n - 1))))
    d = draw(st.lists(st.integers(0, 1), min_size=n, max_size=n))
    p = draw(st.lists(st.sampled_from([0.1, 0.25, 0.3, 0.5, 0.7, 0.9]), min_size=n, max_size=n))
    a = draw(st.lists(st.integers(5, 12).map(float), min_size=n, max_size=n))
    circular = draw(st.booleans())
    return {"n": n, "idx": idx, "d": d, "p": p, "a": a, "circular": circular}


def _run(c: Contig, idx: list[int], seed: int = 0) -> Any:
    genes, replicons = gene_table([c])
    genes, regions, rg = region_frames(genes, {"bgc_1": ("bgc", tags_of(genes, c.name, idx))})
    return _regions().test_all_regions(
        genes, regions, rg, replicons, rd.config(bootstrap=5), np.random.default_rng(seed)
    )


@SETTINGS
@given(rep=replicon(), perm_seed=st.integers(0, 2**31 - 1))
def test_p1_o_minus_e_invariant_to_outside_relabelling(rep: dict[str, Any], perm_seed: int) -> None:
    n, idx = rep["n"], rep["idx"]
    out = [i for i in range(n) if i not in idx]
    order = np.random.default_rng(perm_seed).permutation(len(out))
    d2, p2 = list(rep["d"]), list(rep["p"])
    for j, src in zip(out, [out[k] for k in order], strict=True):
        d2[j], p2[j] = rep["d"][src], rep["p"][src]
    a = _run(Contig("chr1", rep["d"], rep["p"], circular=rep["circular"]), idx)
    b = _run(Contig("chr1", d2, p2, circular=rep["circular"]), idx)
    ra, rb = pick(a, "bgc_1", "o_minus_e"), pick(b, "bgc_1", "o_minus_e")
    for col in ("observed", "expected", "effect"):
        assert ra[col] == pytest.approx(rb[col])
    assert ra["n_perm"] == rb["n_perm"]


@SETTINGS
@given(rep=replicon())
def test_p2_shifted_sets_keep_size(rep: dict[str, Any]) -> None:
    from giap.stats.permutation import circular_shifts

    n, idx = rep["n"], np.array(rep["idx"])
    shifts = circular_shifts(n, idx, rep["circular"])
    out = _regions().o_minus_e_null(np.ones(n), idx, shifts)
    assert out.shape == shifts.shape
    np.testing.assert_allclose(out, np.full(shifts.size, float(idx.size)))


@SETTINGS
@given(rep=replicon(min_n=3, max_n=25))
def test_p3_linear_shifts_never_wrap(rep: dict[str, Any], monkeypatch_shifts: Any) -> None:
    calls = monkeypatch_shifts
    calls.clear()
    n, idx = rep["n"], rep["idx"]
    res = _run(Contig("ctg1", rep["d"], rep["p"], circular=False), idx)
    want = rd.expected_shifts(n, idx, circular=False)
    assert pick(res, "bgc_1", "o_minus_e")["n_perm"] == len(want)
    for shifts in calls:
        for k in shifts.tolist():
            assert max(idx) + k < n or min(idx) + k >= n
        assert sorted(shifts.tolist()) == want


@SETTINGS
@given(rep=replicon(min_n=2, max_n=25))
def test_p4_vectorised_nulls_equal_naive_loop(rep: dict[str, Any]) -> None:
    from giap.stats.permutation import circular_shifts

    R = _regions()
    n, idx = rep["n"], rep["idx"]
    shifts = circular_shifts(n, np.array(idx), rep["circular"])
    r = np.array(rep["d"], float) - np.array(rep["p"], float)
    np.testing.assert_allclose(
        R.o_minus_e_null(r, np.array(idx), shifts), rd.naive_o_minus_e_null(r, idx, shifts),
        atol=1e-12,
    )  # fmt: skip
    det = np.array(rep["d"], bool)
    a = np.where(det, np.array(rep["a"], float), np.nan)
    np.testing.assert_allclose(
        R.hl_null(a, det, np.array(idx), shifts), rd.naive_hl_null(a, det, idx, shifts),
        equal_nan=True, atol=1e-12,
    )  # fmt: skip


@SETTINGS
@given(rep=replicon(min_n=4, max_n=30))
def test_p5_p_values_equal_naive_exhaustive(rep: dict[str, Any]) -> None:
    n, idx = rep["n"], rep["idx"]
    c = Contig("chr1", rep["d"], rep["p"], a=rep["a"], circular=True)
    res = _run(c, idx)
    shifts = list(range(1, n))
    r = np.array(rep["d"], float) - np.array(rep["p"], float)
    v = np.array(rep["p"], float) * (1 - np.array(rep["p"], float))
    oe = pick(res, "bgc_1", "o_minus_e")
    assert oe["n_perm"] == n - 1
    obs = float(r[idx].sum())
    assert oe["effect"] == pytest.approx(obs)
    want = rd.naive_perm_p(rd.studentised(r, v, idx), rd.naive_studentised_null(r, v, idx, shifts))
    assert oe["p_value"] == pytest.approx(want)

    det = np.array(rep["d"], bool)
    a = np.where(det, np.array(rep["a"], float), np.nan)
    hl = pick(res, "bgc_1", "hl_shift")
    x = [a[i] for i in idx if det[i]]
    y = [a[i] for i in range(n) if det[i] and i not in idx]
    hl_obs = rd.naive_hl(x, y)
    assume(rd.naive_v_positive(a, det, idx, [0, *shifts]))
    null = rd.naive_studentised_hl_null(a, det, idx, shifts)
    t_obs = float(rd.naive_studentised_hl_null(a, det, idx, [0])[0])
    if math.isnan(hl_obs):
        assert rd.is_na(hl["effect"]) and rd.is_na(hl["p_value"])
    else:
        assert hl["effect"] == pytest.approx(hl_obs)
        assert hl["n_perm"] == int(np.count_nonzero(~np.isnan(null)))
        want_hl = rd.naive_perm_p(t_obs, null)
        if math.isnan(want_hl):
            assert rd.is_na(hl["p_value"])
        else:
            assert hl["p_value"] == pytest.approx(want_hl)


@SETTINGS
@given(
    rep=replicon(min_n=10, max_n=30),
    cov=st.lists(
        st.tuples(st.integers(10, 1200), st.integers(0, 60), st.floats(-1.5, 1.5)),
        min_size=30, max_size=30,
    ),
)  # fmt: skip
def test_p6_hl_on_residualised_abundance(rep: dict[str, Any], cov: list[Any]) -> None:
    n, idx = rep["n"], rep["idx"]
    L = [c[0] for c in cov[:n]]
    k = [c[1] for c in cov[:n]]
    g = [round(c[2], 3) for c in cov[:n]]
    c = Contig("chr1", rep["d"], rep["p"], a=rep["a"], length_aa=L, n_obs=k, gravy=g)
    res = _run(c, idx)
    det = np.array(rep["d"], bool)
    assume(int(det.sum()) >= 8)  # an identified fit: 4 coefficients, a few residual dof
    a = np.where(det, np.array(rep["a"], float), np.nan)
    ar = rd.naive_residualise(a, det, L, k, g)
    ok = det & np.isfinite(ar)
    x = [ar[i] for i in idx if ok[i]]
    y = [ar[i] for i in range(n) if ok[i] and i not in idx]
    hl_obs = rd.naive_hl(x, y)
    hl = pick(res, "bgc_1", "hl_shift")
    if math.isnan(hl_obs):
        assert rd.is_na(hl["effect"]) and rd.is_na(hl["p_value"])
        return
    assert hl["effect"] == pytest.approx(hl_obs, abs=1e-8)
    scale = float(np.nanmax(np.abs(a[ok]))) if ok.any() else None  # decision 12: raw magnitude
    null = rd.naive_studentised_hl_null(ar, ok, idx, list(range(1, n)), scale=scale)
    obs = rd.naive_studentised_hl_null(ar, ok, idx, [0], scale=scale)
    if not np.all((null == 0) | np.isnan(null)):  # guard not triggered: V must be positive
        assume(rd.naive_v_positive(ar, ok, idx, list(range(n))))
    want = rd.naive_perm_p(float(obs[0]), null)
    if math.isnan(want):
        assert rd.is_na(hl["p_value"])
    else:
        assert hl["p_value"] == pytest.approx(want)


@st.composite
def tu_replicon(draw: st.DrawFn) -> dict[str, Any]:
    """Circular replicon with multi-gene TUs, continuous abundance and random covariates."""
    n = draw(st.integers(10, 28))
    sizes = draw(st.lists(st.integers(1, 4), min_size=n, max_size=n))
    tu: list[str] = []
    for j, s in enumerate(sizes):
        tu += [f"TU{j}"] * s
    tu = tu[:n]
    idx = sorted(draw(st.sets(st.integers(0, n - 1), min_size=1, max_size=n // 2)))
    d = draw(st.lists(st.integers(0, 1), min_size=n, max_size=n))
    a = draw(st.lists(st.floats(5, 15, allow_nan=False), min_size=n, max_size=n))
    L = draw(st.lists(st.integers(10, 1200), min_size=n, max_size=n))
    k = draw(st.lists(st.integers(0, 60), min_size=n, max_size=n))
    g = draw(st.lists(st.floats(-1.5, 1.5, allow_nan=False), min_size=n, max_size=n))
    return {"n": n, "idx": idx, "d": d, "a": a, "tu": tu, "L": L, "k": k,
            "g": [round(x, 3) for x in g]}  # fmt: skip


@SETTINGS
@given(rep=tu_replicon(), c=st.sampled_from([0.125, 0.5, 2.0, 3.0, 10.0]),
       b=st.floats(-20.0, 20.0, allow_nan=False))  # fmt: skip
def test_p7_hl_p_invariant_to_affine_abundance(rep: dict[str, Any], c: float, b: float) -> None:
    n, idx = rep["n"], rep["idx"]
    assume(sum(rep["d"]) >= 8)
    a2 = [c * x + b for x in rep["a"]]
    mk = [Contig("chr1", rep["d"], [0.5] * n, a=av, tu=rep["tu"], length_aa=rep["L"],
                 n_obs=rep["k"], gravy=rep["g"]) for av in (rep["a"], a2)]  # fmt: skip
    r0, r1 = (pick(_run(m, idx), "bgc_1", "hl_shift") for m in mk)
    if rd.is_na(r0["p_value"]):
        assert rd.is_na(r1["p_value"])
        return
    assert r1["p_value"] == pytest.approx(r0["p_value"], abs=1e-12)
    assert r1["n_perm"] == r0["n_perm"]
    assert r1["effect"] == pytest.approx(c * r0["effect"], abs=1e-7)
    for col in ("ci_low", "ci_high"):  # NA CI (too few TU blocks) stays NA
        assert r1[col] == pytest.approx(c * r0[col], abs=1e-7, nan_ok=True)


@SETTINGS
@given(rep=tu_replicon(), c=st.sampled_from([0.125, 0.5, 2.0, 3.0, 10.0]),
       b=st.floats(-20.0, 20.0, allow_nan=False))  # fmt: skip
def test_p7_studentised_t_invariant_to_affine_values(
    rep: dict[str, Any], c: float, b: float
) -> None:
    """T(c a + b) = T(a) at every shift (HL scales by c, V by c^2; constants cancel)."""
    n, idx = rep["n"], np.array(rep["idx"])
    ok = np.array(rep["d"], bool)
    a = np.where(ok, np.array(rep["a"], float), np.nan)
    assume(ok.sum() >= 2 and np.ptp(a[ok]) > 1e-6)
    ks = np.arange(n)
    assume(rd.naive_v_positive(a, ok, list(idx), list(ks)))
    R = _regions()
    t0 = R.studentised_hl_null(a, ok, idx, ks)
    t1 = R.studentised_hl_null(c * a + b, ok, idx, ks)
    np.testing.assert_allclose(t1, t0, rtol=1e-9, atol=1e-9, equal_nan=True)
    want = rd.naive_studentised_hl_null(a, ok, list(idx), list(ks))
    np.testing.assert_allclose(t0, want, rtol=1e-9, atol=1e-12, equal_nan=True)


@pytest.fixture(scope="module")
def monkeypatch_shifts() -> Any:
    """Module-scoped spy on ``choose_shifts`` (hypothesis forbids function-scoped fixtures)."""
    import giap.stats.permutation as perm

    R = _regions()
    mp = pytest.MonkeyPatch()
    orig = perm.choose_shifts
    calls: list[np.ndarray] = []

    def spy(*args: Any, **kw: Any) -> Any:
        out = orig(*args, **kw)
        calls.append(np.asarray(out[0]).copy())
        return out

    mp.setattr(perm, "choose_shifts", spy)
    if hasattr(R, "choose_shifts"):
        mp.setattr(R, "choose_shifts", spy)
    yield calls
    mp.undo()


# --------------------------------------------------------------------------- P8 (decision 15)


@st.composite
def lag_replicon(draw: st.DrawFn) -> dict[str, Any]:
    """Detected flags and values (integers: exact oracle arithmetic is meaningful), a region set,
    and a lag bound (small bounds exercise max_lag; 30 is the default)."""
    n = draw(st.integers(2, 40))
    ok = draw(st.lists(st.booleans(), min_size=n, max_size=n))
    a = draw(st.lists(st.integers(-8, 8).map(float), min_size=n, max_size=n))
    idx = sorted(draw(st.sets(st.integers(0, n - 1), min_size=1, max_size=max(1, n - 1))))
    return {"n": n, "ok": ok, "a": a, "idx": idx, "max_lag": draw(st.sampled_from([0, 1, 3, 30]))}


@SETTINGS
@given(rep=lag_replicon())
def test_p8_lag_autocovariance_equals_oracle(rep: dict[str, Any]) -> None:
    ok = np.array(rep["ok"], bool)
    a = np.where(ok, np.array(rep["a"]), np.nan)
    g = _regions().lag_autocovariance(a, ok, max_lag=rep["max_lag"])
    want = rd.naive_lag_autocov(list(a), ok, rep["max_lag"])
    np.testing.assert_allclose(g, np.array(want, float), rtol=1e-9, atol=1e-12)


@SETTINGS
@given(rep=lag_replicon(), g=st.lists(st.floats(-2.0, 5.0, allow_nan=False), min_size=31,
                                      max_size=31))  # fmt: skip
def test_p8_hl_variance_equals_oracle(rep: dict[str, Any], g: list[float]) -> None:
    """Any gamma vector (the quadratic form, not only estimated gammas); NaN where a side is
    empty."""
    ok = np.array(rep["ok"], bool)
    ks = np.arange(rep["n"])
    v = _regions().hl_variance(ok, np.array(rep["idx"]), ks, np.array(g))
    want = rd.naive_lag_hl_variance(ok, rep["idx"], list(ks), g)
    np.testing.assert_allclose(v, np.array(want, float), rtol=1e-9, atol=1e-9, equal_nan=True)


@SETTINGS
@given(rep=tu_replicon(), data=st.data())
def test_p8_hl_row_independent_of_tu_partition(rep: dict[str, Any], data: Any) -> None:
    """Any other TU partition of the same genes (random contiguous blocks, labels drawn at random,
    so relabelling, fragmenting and merging all occur): effect, p and n_perm unchanged."""
    n, idx = rep["n"], rep["idx"]
    assume(sum(rep["d"]) >= 8)
    cuts = data.draw(st.lists(st.booleans(), min_size=n, max_size=n))
    names = data.draw(st.lists(st.integers(0, 10**6), min_size=n, max_size=n, unique=True))
    blk, tu2 = 0, []
    for i in range(n):
        blk += int(cuts[i])
        tu2.append(f"Q{names[blk]}")
    mk = [Contig("chr1", rep["d"], [0.5] * n, a=rep["a"], tu=t, length_aa=rep["L"],
                 n_obs=rep["k"], gravy=rep["g"]) for t in (rep["tu"], tu2)]  # fmt: skip
    r0, r1 = (pick(_run(m, idx), "bgc_1", "hl_shift") for m in mk)
    if rd.is_na(r0["p_value"]):
        assert rd.is_na(r1["p_value"])
        return
    assert r1["effect"] == pytest.approx(r0["effect"], abs=1e-12)
    assert r1["p_value"] == pytest.approx(r0["p_value"], abs=1e-12)
    assert r1["n_perm"] == r0["n_perm"]
