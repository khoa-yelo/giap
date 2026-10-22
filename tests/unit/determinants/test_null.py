"""M8-09 TU-size-preserving permutation null: examples 1-4 and 8, the degenerate rows and the
§7.8.3 units "TU-size-preserving permutation blocks" and "the p formula"
(docs/specs/M8-09-determinants-attribution-null.md rules 5-7, decision 16; GIAP_SPEC.md §7.8.2
step 4, §8.2).

Contract choices the interface leaves open are in ``tests/helpers/determinants_attrib_data.py``.
"""

from __future__ import annotations

import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_attrib_data as a
import determinants_models_data as d

from giap.rng import for_step
from giap.stats.permutation import perm_p

SEED = 20261021


def _rng(seed: int = SEED) -> np.random.Generator:
    return for_step(a.NULL_STREAM, seed)


def _world(n: int = 120, seed: int = 3) -> tuple[pd.DataFrame, np.ndarray, list[str], np.ndarray]:
    """n rows, TUs of 1-4 genes, two informative columns; repeat-1 folds from ``tu_folds``."""
    rng = np.random.default_rng(seed)
    tu = d.tu_ids(d.mixed_sizes(rng, n, 1, 4))
    X = pd.DataFrame(
        {"cai": rng.normal(size=n), "gravy": rng.normal(size=n), "sd_dg": rng.normal(size=n)},
        index=pd.Index([f"G{i:05d}" for i in range(n)], name="locus_tag"),
    )
    y = 1.0 * X["cai"].to_numpy() + rng.normal(0.0, 1.0, n)
    folds = d.mod().tu_folds(tu, 5, 1, for_step("determinants.cv", seed))
    return X, np.asarray(y, dtype=float), tu, np.asarray(folds[0], dtype=np.int64)


def _cfg(**lgb: Any) -> Any:
    return d.config(fast=True, n_estimators=20, **lgb)


def _null(
    X: pd.DataFrame, y: np.ndarray, tu: list[str], f1: np.ndarray, obs: dict[str, float], B: int,
    seed: int = SEED, threads: int = 1, cfg: Any = None,
) -> dict[str, Any]:  # fmt: skip
    cfg = _cfg() if cfg is None else cfg
    return dict(a.null_mod().null_p(X, y, tu, f1, obs, cfg, B, seed, threads))


def _flags(res: Any) -> set[str]:
    f = res.flags
    if isinstance(f, str):
        return {x for x in f.split(";") if x}
    return set(f)


# ---------------------------------------------------------------------------------------------
# rule 5 / example 1: permutation blocks

EX1_TU = a.tu_from_sizes([2, 2, 3, 1, 2], ["A", "B", "C", "D", "E"])  # rows A 0-1, B 2-3, C 4-6,
# D 7, E 8-9


def test_example_1_tu_size_blocks() -> None:
    got = a.null_mod().tu_size_blocks(EX1_TU)
    assert isinstance(got, list) and len(got) == 3
    for arr in got:
        assert np.asarray(arr).dtype == np.int64
    expect = [np.array([[7]]), np.array([[0, 1], [2, 3], [8, 9]]), np.array([[4, 5, 6]])]
    for g, e in zip(got, expect, strict=True):
        np.testing.assert_array_equal(np.asarray(g), e)


def test_example_1_permutation_never_moves_a_two_gene_block_onto_c_or_d() -> None:
    y = np.arange(10, dtype=float) * 10.0
    rng = _rng()
    blocks_ab_e = {(0.0, 10.0), (20.0, 30.0), (80.0, 90.0)}
    seen: set[tuple[tuple[float, float], ...]] = set()
    for _ in range(300):
        p = np.asarray(a.null_mod().permute_response(y, EX1_TU, rng))
        np.testing.assert_array_equal(p[4:8], y[4:8])  # C and D fixed
        pairs = (tuple(p[0:2]), tuple(p[2:4]), tuple(p[8:10]))
        assert set(pairs) == blocks_ab_e
        seen.add(pairs)
    assert len(seen) == 6  # every permutation of {A, B, E} occurs


def test_example_2_block_swap_keeps_gene_order() -> None:
    y = np.array([1.0, 2.0, 11.0, 12.0, 21.0, 22.0, 23.0])  # a1 a2 | b1 b2 | c1 c2 c3
    tu = a.tu_from_sizes([2, 2, 3], ["A", "B", "C"])
    swapped = np.array([11.0, 12.0, 1.0, 2.0, 21.0, 22.0, 23.0])
    rng = _rng(1)
    outcomes: Counter[str] = Counter()
    for _ in range(200):
        p = np.asarray(a.null_mod().permute_response(y, tu, rng), dtype=float)
        if np.array_equal(p, swapped):
            outcomes["swap"] += 1
        else:
            np.testing.assert_array_equal(p, y)
            outcomes["identity"] += 1
    assert outcomes["swap"] > 0 and outcomes["identity"] > 0


def test_permute_response_is_uniform_within_a_size_group() -> None:
    """Exchangeability needs a uniform permutation: 3 TUs of size 2 -> each of the 6 orders at
    rate 1/6 (6000 draws, ±5 binomial SD)."""
    y = np.array([0.0, 1.0, 10.0, 11.0, 20.0, 21.0])
    tu = a.tu_from_sizes([2, 2, 2], ["A", "B", "C"])
    rng = _rng(2)
    counts = Counter(
        tuple(np.asarray(a.null_mod().permute_response(y, tu, rng))[[0, 2, 4]].tolist())
        for _ in range(6000)
    )
    assert len(counts) == 6
    sd = math.sqrt(6000 * (1 / 6) * (5 / 6))
    for c in counts.values():
        assert abs(c - 1000) <= 5 * sd


def test_permute_response_non_contiguous_tu_rows() -> None:
    """A TU's genes need not be adjacent: block A = rows (0, 2), B = rows (1, 3)."""
    y = np.array([1.0, 11.0, 2.0, 12.0])
    tu = ["A", "B", "A", "B"]
    rng = _rng(3)
    seen = set()
    for _ in range(100):
        p = tuple(np.asarray(a.null_mod().permute_response(y, tu, rng)).tolist())
        assert p in {(1.0, 11.0, 2.0, 12.0), (11.0, 1.0, 12.0, 2.0)}
        seen.add(p)
    assert len(seen) == 2


def test_permute_response_pure_and_reproducible() -> None:
    y = np.arange(10, dtype=float)
    y0 = y.copy()
    p1 = np.asarray(a.null_mod().permute_response(y, EX1_TU, _rng(5)))
    p2 = np.asarray(a.null_mod().permute_response(y, EX1_TU, _rng(5)))
    np.testing.assert_array_equal(y, y0)
    np.testing.assert_array_equal(p1, p2)
    assert p1.dtype == np.float64 and p1.shape == y.shape


def test_permute_response_all_sizes_unique_is_identity() -> None:
    tu = a.tu_from_sizes([1, 2, 3, 4])
    y = np.arange(10, dtype=float)
    for s in range(20):
        np.testing.assert_array_equal(a.null_mod().permute_response(y, tu, _rng(s)), y)


# ---------------------------------------------------------------------------------------------
# rule 6 / examples 3-4: the p formula


def test_example_3_p_formula() -> None:
    """B = 9, two null R² >= 0.45 -> p = (1 + 2)/(1 + 9) = 0.3 (the formula null_p must use)."""
    null = np.array([0.1, 0.2, 0.3, 0.5, 0.05, 0.6, 0.0, 0.1, 0.2])
    assert perm_p(0.45, null, alternative="greater") == pytest.approx(0.3, abs=1e-15)


def test_null_p_uses_the_formula_on_its_null_vector() -> None:
    X, y, tu, f1 = _world()
    obs = {"lightgbm": 0.0, "elastic_net": 0.0}  # below most nulls: b > 0 for both models
    B = 9
    res = _null(X, y, tu, f1, obs, B)
    assert set(res) == set(obs)
    for name, r in res.items():
        null = np.asarray(r.null_r2, dtype=float)
        assert null.shape == (B,) and np.all(np.isfinite(null)), name
        b = int(np.count_nonzero(null >= obs[name]))
        assert r.b_count == b, name
        assert r.n_perm == B, name
        assert r.p == pytest.approx((1 + b) / (1 + B), abs=1e-15), name
        assert r.p == pytest.approx(perm_p(obs[name], null, alternative="greater"), abs=1e-15)
        assert r.observed == pytest.approx(obs[name]), name
        assert "null_degenerate" not in _flags(r), name


def test_example_4_no_null_exceeds_gives_floor() -> None:
    """Observed R² = 1 (no null refit reaches it) -> b = 0, p = 1/(1 + B); fast form B = 19."""
    X, y, tu, f1 = _world()
    res = _null(X, y, tu, f1, {"lightgbm": 1.0}, 19)
    assert set(res) == {"lightgbm"}
    assert res["lightgbm"].b_count == 0
    assert res["lightgbm"].p == pytest.approx(1 / 20, abs=1e-15)


@pytest.mark.slow
def test_example_4_B_999() -> None:
    X, y, tu, f1 = _world()
    res = _null(X, y, tu, f1, {"lightgbm": 1.0}, 999, threads=4)
    assert res["lightgbm"].b_count == 0
    assert res["lightgbm"].p == pytest.approx(1 / 1000, abs=1e-15)


def test_only_requested_models_run() -> None:
    X, y, tu, f1 = _world()
    res = _null(X, y, tu, f1, {"elastic_net": 0.3}, 3)
    assert set(res) == {"elastic_net"}


# ---------------------------------------------------------------------------------------------
# rule 6: stream, child seeds, repeat-1 folds (decision 16)


@pytest.mark.parametrize("model", ["elastic_net", "lightgbm"])
def test_null_refits_are_one_cv_repeat_on_repeat1_folds(model: str) -> None:
    """Null b = the pooled-OOF R² of ``cross_validate`` on y permuted with child seed b, one
    repeat on ``folds_rep1``. LightGBM without bagging or feature sampling is seed-free, so the
    refit is reproducible outside ``null_p`` (the elastic net has no randomness)."""
    X, y, tu, f1 = _world()
    cfg = _cfg(subsample=1.0, colsample_bytree=1.0)
    B = 4
    res = _null(X, y, tu, f1, {model: 0.2}, B, cfg=cfg)[model]
    expect = []
    for child in a.null_children(SEED, B):
        yb = np.asarray(a.null_mod().permute_response(y, tu, np.random.default_rng(child)))
        cv = d.mod().cross_validate(X, yb, f1[None, :], {}, cfg, SEED, 1, tu=tu)
        expect.append(float(cv[model].r2[0]))
    np.testing.assert_allclose(np.asarray(res.null_r2, dtype=float), expect, rtol=0, atol=1e-9)


def test_null_vector_prefix_stable_in_B() -> None:
    """Child seed b does not depend on B (SeedSequence.spawn): B = 3 is the head of B = 6."""
    X, y, tu, f1 = _world()
    obs = {"lightgbm": 0.2, "elastic_net": 0.2}
    r3, r6 = _null(X, y, tu, f1, obs, 3), _null(X, y, tu, f1, obs, 6)
    for m in obs:
        np.testing.assert_array_equal(np.asarray(r6[m].null_r2)[:3], np.asarray(r3[m].null_r2))


def test_seed_changes_the_null_and_reproduces() -> None:
    X, y, tu, f1 = _world()
    obs = {"elastic_net": 0.2}
    r1, r1b = _null(X, y, tu, f1, obs, 4, seed=1), _null(X, y, tu, f1, obs, 4, seed=1)
    r2 = _null(X, y, tu, f1, obs, 4, seed=2)
    np.testing.assert_array_equal(r1["elastic_net"].null_r2, r1b["elastic_net"].null_r2)
    assert not np.array_equal(r1["elastic_net"].null_r2, r2["elastic_net"].null_r2)


def test_example_8_threads_identical() -> None:
    X, y, tu, f1 = _world()
    obs = {"lightgbm": 0.2, "elastic_net": 0.2}
    r1 = _null(X, y, tu, f1, obs, 8, threads=1)
    r8 = _null(X, y, tu, f1, obs, 8, threads=8)
    for m in obs:
        np.testing.assert_array_equal(r1[m].null_r2, r8[m].null_r2)
        assert r1[m].p == r8[m].p and r1[m].b_count == r8[m].b_count


# ---------------------------------------------------------------------------------------------
# degenerate rows


def test_every_tu_size_unique_is_degenerate_p_one() -> None:
    """Sizes 1..14 (105 rows): no movable block -> b = B, p = 1, ``null_degenerate``."""
    sizes = list(range(1, 15))
    tu = a.tu_from_sizes(sizes)
    n = len(tu)
    rng = np.random.default_rng(4)
    X = pd.DataFrame(
        {"cai": rng.normal(size=n), "gravy": rng.normal(size=n)},
        index=pd.Index([f"G{i:05d}" for i in range(n)], name="locus_tag"),
    )
    y = X["cai"].to_numpy() + rng.normal(size=n)
    f1 = np.asarray(d.mod().tu_folds(tu, 5, 1, for_step("determinants.cv", 4))[0])
    B = 7
    res = _null(X, y, tu, f1, {"lightgbm": 0.9, "elastic_net": 0.9}, B)
    for name, r in res.items():
        assert r.b_count == B, name
        assert r.p == 1.0, name
        assert r.n_perm == B, name
        assert "null_degenerate" in _flags(r), name


def test_partly_unique_sizes_not_degenerate() -> None:
    X, y, tu, f1 = _world()
    res = _null(X, y, tu, f1, {"elastic_net": 0.0}, 2)
    assert "null_degenerate" not in _flags(res["elastic_net"])
