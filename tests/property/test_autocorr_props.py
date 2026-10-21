"""Hypothesis properties of the autocorrelation engine (issue M5-02; GIAP_SPEC.md §7.3.5(c), §7.3.6,
§8.3; docs/specs/M5-02-stats-moran-joincount.md acceptance criteria, lead decision 7).

- every TU-block permutation row is a permutation of 0 ... n - 1, keeps each block contiguous
  (cyclically) in its internal order (a wrapped block: end-part then start-part), and preserves
  the multiset of block sizes;
- label permutations with ``movable`` never move fixed positions;
- ``morans_i_matrix`` / ``join_count_matrix`` (raw and Pearson residual, decision 17; p drawn
  from [0, 1] with 0 and 1 likely, i.e. missing genes) equal scalar loops and an independent
  oracle on random inputs, valid and invalid lags alike (the join-count comparison is made at a
  tolerance scaled to the sum's conditioning, ``_bb_tol``: near p = 0 or 1 the oracle's naive
  summation, not the engine, loses digits);
- Moran's I is invariant to affine maps x -> a x + b (a != 0); I_k in [-1, 1] is NOT asserted;
- ``lag_test`` p equals ``perm_p`` unless the null is constant (then NaN, decision 6);
- ``choose_block_len`` lies in [1, min(max_len, max_lag + 1, n)] and its draws depend on B and n
  only; ``block_len_prime`` = min(20, 2 l-hat + 1, T) (decision 10).
"""

from __future__ import annotations

import math

import numpy as np
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from giap.stats.autocorr import (
    block_len_prime,
    choose_block_len,
    join_count,
    join_count_matrix,
    label_permutations,
    lag_pairs,
    lag_test,
    max_lag,
    morans_i,
    morans_i_matrix,
    tu_block_permutations,
    tu_blocks,
)
from giap.stats.permutation import perm_p

FAST = settings(max_examples=150, deadline=None)
SEED = st.integers(0, 2**32 - 1)


# ----------------------------------------------------------------------------- oracles


def _oracle_pairs(n: int, k: int, circular: bool) -> list[tuple[int, int]]:
    if circular:
        if not 1 <= k <= (n - 1) // 2:
            return []
        return [(i, (i + k) % n) for i in range(n)]
    if not 1 <= k <= n - 1:
        return []
    return [(i, i + k) for i in range(n - k)]


def _oracle_moran(x: np.ndarray, k: int, circular: bool) -> float:
    pairs = _oracle_pairs(len(x), k, circular)
    has = ~np.isnan(x)
    n = int(has.sum())
    if not pairs or n < 2:
        return math.nan
    z = np.where(has, x - x[has].mean(), 0.0)
    ss = float((z[has] ** 2).sum())
    if ss <= 1e-24 * max(1.0, float(np.max(x[has] ** 2))):
        return math.nan
    good = [(i, j) for i, j in pairs if has[i] and has[j]]
    if not good:
        return math.nan
    num = sum(z[i] * z[j] for i, j in good) / len(good)
    return float(num / (ss / n))


def _oracle_bb(d: np.ndarray, k: int, circular: bool, p: np.ndarray | None) -> tuple[float, float]:
    """Raw Σ d_i d_j, or (decision 17) Pearson Σ r_i r_j, r = (d - p)/sqrt(p(1 - p)), over pairs
    where both genes have a residual (p in (0, 1)).

    Returns
    -------
    tuple
        The sum and ``scale`` = Σ |r_i r_j| over the same pairs, the conditioning of the sum
        (``_bb_tol``). With p near 0 or 1 (unclipped, decision 17) the residuals reach ~3e4, so
        the terms reach ~1e8 while the sum is O(1): a naive left-to-right ``sum`` of m <= 25 such
        terms carries an error of up to ~5.5e-15 * scale, and it is this oracle, not the engine
        (correctly rounded, checked with ``math.fsum``), that loses the digits.
    """
    pairs = _oracle_pairs(len(d), k, circular)
    if not pairs:
        return math.nan, 0.0
    if p is None:
        r = d.astype(float)
        ok = np.ones(len(d), dtype=bool)
    else:
        ok = (p > 0) & (p < 1)
        r = np.zeros(len(d))
        r[ok] = (d[ok] - p[ok]) / np.sqrt(p[ok] * (1 - p[ok]))
    terms = [r[i] * r[j] for i, j in pairs if ok[i] and ok[j]]
    return float(sum(terms)), float(sum(abs(t) for t in terms))


def _bb_tol(scale: float) -> float:
    """Absolute tolerance for a join-count sum of conditioning ``scale`` = Σ |terms|."""
    return 1e-9 + 1e-12 * scale


# ----------------------------------------------------------------------------- strategies


@st.composite
def tu_labels(draw: st.DrawFn) -> list[str]:
    """Runs of TU labels from a small alphabet: adjacent runs differ, a label may recur later
    (interrupted TU), and the first and last runs often share a label (origin-crossing TU)."""
    n_runs = draw(st.integers(0, 8))
    labels: list[str] = []
    prev = None
    for _ in range(n_runs):
        lab = draw(st.sampled_from("ABCDE").filter(lambda c, p=prev: c != p))
        labels.extend([lab] * draw(st.integers(1, 4)))
        prev = lab
    return labels


def _expected_blocks(blocks: np.ndarray) -> dict[int, list[int]]:
    """Positions of each block in emission order (wrapped block: end-part then start-part)."""
    n = blocks.size
    out: dict[int, list[int]] = {}
    for i in range(n):
        out.setdefault(int(blocks[i]), []).append(i)
    for b, pos in out.items():
        gaps = [j for j in range(1, len(pos)) if pos[j] != pos[j - 1] + 1]
        if gaps:  # a wrapped block: the start-part [0 ...] and the end-part [... n-1]
            assert len(gaps) == 1 and pos[0] == 0 and pos[-1] == n - 1, (b, pos)
            g = gaps[0]
            out[b] = pos[g:] + pos[:g]
    return out


# ----------------------------------------------------------------------------- TU blocks


@FAST
@given(tu_labels(), st.booleans())
def test_tu_blocks_are_maximal_runs(labels: list[str], circular: bool) -> None:
    blocks = tu_blocks(labels, circular)
    assert blocks.shape == (len(labels),) and blocks.dtype == np.int64
    n = len(labels)
    for i in range(1, n):
        assert (blocks[i] == blocks[i - 1]) == (labels[i] == labels[i - 1])
    if n:
        assert blocks[0] == 0
        wrap = circular and n > 1 and labels[0] == labels[-1]
        assert (blocks[-1] == blocks[0]) == (wrap or len(set(blocks.tolist())) == 1)
    _expected_blocks(blocks)  # each block a run, or one wrapped run


@FAST
@given(tu_labels(), st.booleans(), st.integers(1, 12), SEED)
def test_tu_block_permutations_keep_blocks(
    labels: list[str], circular: bool, B: int, seed: int
) -> None:
    blocks = tu_blocks(labels, circular)
    n = blocks.size
    out = tu_block_permutations(blocks, B, np.random.default_rng(seed))
    assert out.shape == (B, n) and out.dtype == np.int64
    expected = _expected_blocks(blocks)
    sizes = sorted(len(v) for v in expected.values())
    for row in out:
        assert sorted(row.tolist()) == list(range(n))
        rows = row.tolist()
        for seq in expected.values():
            start = rows.index(seq[0])
            assert [rows[(start + t) % n] for t in range(len(seq))] == seq
        # cyclic runs of block ids along the row = the blocks, sizes preserved
        ids = blocks[row]
        if len(expected) > 1:
            cuts = [t for t in range(n) if ids[t] != ids[t - 1]]
            run_sizes = sorted(
                (cuts[(m + 1) % len(cuts)] - cuts[m]) % n or n for m in range(len(cuts))
            )
            assert run_sizes == sizes
        elif n:
            assert rows == expected[int(blocks[0])]


@FAST
@given(st.integers(0, 30), st.integers(1, 10), SEED, st.data())
def test_label_permutations_movable_never_moves_fixed(
    n: int, B: int, seed: int, data: st.DataObject
) -> None:
    movable = np.array(data.draw(st.lists(st.booleans(), min_size=n, max_size=n)), dtype=bool)
    out = label_permutations(n, B, np.random.default_rng(seed), movable=movable)
    assert out.shape == (B, n) and out.dtype == np.int64
    fixed = np.flatnonzero(~movable)
    mov = np.flatnonzero(movable)
    for row in out:
        assert sorted(row.tolist()) == list(range(n))
        np.testing.assert_array_equal(row[fixed], fixed)
        assert sorted(row[mov].tolist()) == mov.tolist()


# ----------------------------------------------------------------------------- matrices


VALUES = st.one_of(st.floats(-50, 50, allow_nan=False), st.just(math.nan))


@FAST
@given(
    st.integers(0, 25).flatmap(
        lambda n: st.tuples(st.just(n), st.lists(VALUES, min_size=n, max_size=n))
    ),
    st.lists(st.integers(-1, 28), min_size=1, max_size=6),
    st.booleans(),
    st.integers(1, 6),
    SEED,
)
def test_morans_i_matrix_equals_scalar_and_oracle(
    nx: tuple[int, list[float]], lags: list[int], circular: bool, B: int, seed: int
) -> None:
    n, vals = nx
    x = np.array(vals, dtype=float)
    orders = label_permutations(n, B, np.random.default_rng(seed))
    X = x[orders] if n else np.zeros((B, 0))
    got = morans_i_matrix(X, lags, circular)
    assert got.shape == (B, len(lags))
    for b in range(B):
        for c, k in enumerate(lags):
            want = _oracle_moran(X[b], k, circular)
            scalar = morans_i(X[b], k, circular)
            if math.isnan(want):
                assert math.isnan(scalar) and math.isnan(got[b, c]), (b, k)
            else:
                assert math.isclose(scalar, want, rel_tol=1e-9, abs_tol=1e-9)
                assert math.isclose(got[b, c], scalar, rel_tol=1e-9, abs_tol=1e-9)


@FAST
@given(
    st.integers(0, 25).flatmap(
        lambda n: st.tuples(
            st.lists(st.booleans(), min_size=n, max_size=n),
            st.lists(
                st.one_of(
                    st.floats(0, 1),
                    st.sampled_from([0.0, 1.0, 1e-9, 1 - 1e-9, math.nan]),
                ),  # in [0, 1] or NaN: out-of-range p raises (decision 21)
                min_size=n,
                max_size=n,
            ),
        )
    ),
    st.lists(st.integers(-1, 28), min_size=1, max_size=6),
    st.booleans(),
    st.booleans(),
    st.integers(1, 6),
    SEED,
)
def test_join_count_matrix_equals_scalar_and_oracle(
    dp: tuple[list[bool], list[float]],
    lags: list[int],
    circular: bool,
    residual: bool,
    B: int,
    seed: int,
) -> None:
    d = np.array(dp[0], dtype=bool)
    p = np.array(dp[1], dtype=float)
    n = d.size
    orders = label_permutations(n, B, np.random.default_rng(seed))
    D = d[orders] if n else np.zeros((B, 0), dtype=bool)
    P = (p[orders] if n else np.zeros((B, 0))) if residual else None
    got = join_count_matrix(D, lags, circular, P=P)
    assert got.shape == (B, len(lags))
    for b in range(B):
        pb = None if P is None else P[b]
        for c, k in enumerate(lags):
            want, scale = _oracle_bb(D[b], k, circular, pb)
            scalar = join_count(D[b], k, circular, p=pb)
            if math.isnan(want):
                assert math.isnan(scalar) and math.isnan(got[b, c]), (b, k)
            else:
                tol = _bb_tol(scale)
                assert math.isclose(scalar, want, rel_tol=1e-9, abs_tol=tol), (b, k, scale)
                assert math.isclose(got[b, c], scalar, rel_tol=1e-9, abs_tol=tol), (b, k, scale)


@FAST
@given(st.integers(0, 40), st.integers(-2, 45), st.booleans())
def test_lag_pairs_match_oracle(n: int, k: int, circular: bool) -> None:
    i, j = lag_pairs(n, k, circular)
    assert sorted(zip(i.tolist(), j.tolist(), strict=True)) == sorted(_oracle_pairs(n, k, circular))
    assert (i.size > 0) == (1 <= k <= max_lag(n, circular))


# ----------------------------------------------------------------------------- Moran invariance


@FAST
@given(
    st.lists(st.one_of(st.floats(-10, 10, allow_nan=False), st.just(math.nan)), min_size=3,
             max_size=30),
    st.integers(1, 14),
    st.booleans(),
    st.one_of(st.floats(-100, -0.01), st.floats(0.01, 100)),
    st.floats(-100, 100),
)  # fmt: skip
def test_morans_i_affine_invariant(
    vals: list[float], k: int, circular: bool, a: float, b: float
) -> None:
    x = np.array(vals, dtype=float)
    has = ~np.isnan(x)
    assume(has.sum() >= 2 and np.ptp(x[has]) > 0.1)
    i0 = morans_i(x, k, circular)
    i1 = morans_i(a * x + b, k, circular)
    if math.isnan(i0):
        assert math.isnan(i1)
    else:
        assert math.isclose(i0, i1, rel_tol=1e-7, abs_tol=1e-9)


# ----------------------------------------------------------------------------- lag_test


@FAST
@given(
    st.floats(-20, 20, allow_nan=False),
    st.lists(st.one_of(st.integers(-3, 3).map(float), st.just(math.nan)), max_size=40),
)
def test_lag_test_p_is_perm_p_unless_null_constant(obs: float, null: list[float]) -> None:
    res = lag_test(obs, null)
    valid = np.array([v for v in null if not math.isnan(v)])
    assert res.n_valid == valid.size
    assert res.value == obs
    if valid.size == 0:
        assert math.isnan(res.p_value) and math.isnan(res.null_mean)
    elif np.ptp(valid) == 0:
        assert math.isnan(res.p_value)
        assert res.null_mean == valid[0] == res.null_lo == res.null_hi
    else:
        assert res.p_value == perm_p(obs, null, "two-sided")
        # lead decision 15: lo <= mean <= hi is not a theorem (a 0.975 quantile can fall below
        # the mean); the mean lies within the range and the envelope is ordered
        assert valid.min() <= res.null_mean <= valid.max()
        assert res.null_lo <= res.null_hi


# ----------------------------------------------------------------------------- choose_block_len


@settings(max_examples=40, deadline=None)
@given(
    st.integers(0, 60),
    st.booleans(),
    st.integers(1, 25),
    SEED,
    st.sampled_from(["iid", "ar", "alternating", "constant"]),
)
def test_choose_block_len_bounds(
    n: int, circular: bool, max_len: int, seed: int, kind: str
) -> None:
    rng = np.random.default_rng(seed)
    if kind == "iid":
        z = rng.normal(size=n)
    elif kind == "ar":
        z = np.cumsum(rng.normal(size=n))
    elif kind == "alternating":
        z = np.array([1.0, -1.0] * (n // 2 + 1))[:n]
    else:
        z = np.full(n, 3.0)
    ell = choose_block_len(z, circular, 39, rng, max_len=max_len)
    assert isinstance(ell, int)
    if n < 3 or kind == "constant":
        assert ell == 1
    else:
        assert 1 <= ell <= min(max_len, max_lag(n, circular) + 1, n)


@settings(max_examples=25, deadline=None)
@given(st.integers(3, 80), st.integers(1, 30), SEED, SEED)
def test_choose_block_len_draws_depend_on_b_and_n_only(
    n: int, B: int, seed_z: int, seed_rng: int
) -> None:
    z1 = np.random.default_rng(seed_z).normal(size=n)
    z2 = np.array([1.0, -1.0] * n)[:n]
    after = []
    for z in (z1, z2):
        rng = np.random.default_rng(seed_rng)
        choose_block_len(z, True, B, rng)
        after.append(rng.bit_generator.state)
    assert after[0] == after[1]


@FAST
@given(st.integers(1, 20), st.integers(1, 5_000))
def test_block_len_prime_formula(l_hat: int, T: int) -> None:
    got = block_len_prime(l_hat, T)
    assert got == min(20, 2 * l_hat + 1, T)
    assert 1 <= got <= min(20, T)
