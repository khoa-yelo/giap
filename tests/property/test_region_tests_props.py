"""Properties of the M4-17 region tests (docs/specs/M4-17-region-tests.md acceptance "Property";
GIAP_SPEC.md §7.3.6).

- P1 O - E (observed, expected, effect) is unchanged when the values of genes outside the region
  are permuted among the outside indices;
- P2 every shifted index set has |I| members (o_minus_e_null of ones is |I| at every shift);
- P3 on linear contigs no shifted set wraps: every shift used satisfies max(I) + k < N or
  min(I) + k >= N, and the count equals the spec's rule (``expected_shifts``);
- P4 the vectorised nulls equal a naive per-shift loop (O - E and HL, NaN where a side is empty);
- P5 end to end on circular contigs below the threshold: n_perm = N - 1 (O - E), and p_value
  equals the two-sided centred p of the naive null (both statistics; NaN shifts dropped).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings
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
    oe = pick(res, "bgc_1", "o_minus_e")
    assert oe["n_perm"] == n - 1
    obs = float(r[idx].sum())
    assert oe["effect"] == pytest.approx(obs)
    want = rd.naive_perm_p(obs, rd.naive_o_minus_e_null(r, idx, shifts))
    assert oe["p_value"] == pytest.approx(want)

    det = np.array(rep["d"], bool)
    a = np.where(det, np.array(rep["a"], float), np.nan)
    hl = pick(res, "bgc_1", "hl_shift")
    x = [a[i] for i in idx if det[i]]
    y = [a[i] for i in range(n) if det[i] and i not in idx]
    hl_obs = rd.naive_hl(x, y)
    null = rd.naive_hl_null(a, det, idx, shifts)
    if math.isnan(hl_obs):
        assert rd.is_na(hl["effect"]) and rd.is_na(hl["p_value"])
    else:
        assert hl["effect"] == pytest.approx(hl_obs)
        assert hl["n_perm"] == int(np.count_nonzero(~np.isnan(null)))
        want_hl = rd.naive_perm_p(hl_obs, null)
        if math.isnan(want_hl):
            assert rd.is_na(hl["p_value"])
        else:
            assert hl["p_value"] == pytest.approx(want_hl)


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
