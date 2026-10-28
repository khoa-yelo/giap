"""Target-decoy q-values, monotone q and picked FDR (issue M4-01; GIPPA_SPEC.md §7.2.5, §7.3.2,
§8.3, Appendix B.6; Savitski 2015 picked rule: ties go to the decoy).

Conventions pinned here (M3-10 reuses ``tda_qvalues``): higher score = better; NaN and -inf =
absent; thresholds are the distinct observed scores and counts include every score >= t, so a decoy
tied with a target counts at that threshold (M4-01 Q1); output is one q per target in input order.
"""

from __future__ import annotations

import dataclasses
import inspect

import numpy as np
import pytest

from gippa.stats.tda import PickedResult, monotone_q, picked_fdr, tda_qvalues

NAN = float("nan")
INF = float("inf")


def close(actual: object, expected: object, atol: float = 1e-12) -> None:
    np.testing.assert_allclose(
        np.asarray(actual, dtype=float), expected, rtol=0, atol=atol, equal_nan=True
    )


# ---------------------------------------------------------------------------- interface (pinned)


def test_tda_qvalues_signature_is_pinned_for_m3_10() -> None:
    sig = inspect.signature(tda_qvalues)
    assert list(sig.parameters) == ["target_scores", "decoy_scores", "pi", "plus_one"]
    assert sig.parameters["pi"].default == 1.0
    assert sig.parameters["plus_one"].default is True


def test_picked_fdr_signature_is_pinned() -> None:
    assert list(inspect.signature(picked_fdr).parameters) == ["target_best", "decoy_best"]


def test_tda_qvalues_returns_float64_one_per_target() -> None:
    q = tda_qvalues(target_scores=[3.0, 2.0, 1.0], decoy_scores=[2.0], pi=1.0, plus_one=True)
    assert isinstance(q, np.ndarray)
    assert q.dtype == np.float64
    assert q.shape == (3,)


# ---------------------------------------------------------------------------- monotone_q


def test_monotone_q_running_min_from_loose_end() -> None:
    # FDR of example 8 at thresholds 5, 4, 3.6, 3, 2, 1 (strictest -> loosest)
    fdr = [1.0, 0.5, 1.0, 2 / 3, 1.0, 0.75]
    close(monotone_q(fdr), [0.5, 0.5, 2 / 3, 2 / 3, 0.75, 0.75])


def test_monotone_q_clips_to_unit_interval() -> None:
    close(monotone_q([2.0, 1.5]), [1.0, 1.0])
    close(monotone_q([3.0, 0.2, 4.0]), [0.2, 0.2, 1.0])


def test_monotone_q_empty() -> None:
    q = monotone_q([])
    assert q.shape == (0,)
    assert q.dtype == np.float64


# ---------------------------------------------------------------------------- tda_qvalues: examples


def test_example4_b6_150_targets_before_first_decoy() -> None:
    targets = np.arange(1000, 850, -1, dtype=float)
    assert targets.size == 150
    q = tda_qvalues(targets, [850.0])
    close(q, np.full(150, 1 / 150))
    assert np.all(q <= 0.01)  # all 150 accepted at 0.01
    assert round(float(q[0]), 4) == 0.0067


def test_example5_b6_50_targets_no_decoys() -> None:
    q = tda_qvalues(np.arange(50, 0, -1, dtype=float), [])
    close(q, np.full(50, 0.02))
    assert not np.any(q <= 0.01)  # none accepted at 0.01


def test_example6_plus_one_false() -> None:
    close(tda_qvalues([3, 2, 1], [2], plus_one=False), [0.0, 1 / 3, 1 / 3])


def test_example7_pi_half() -> None:
    close(tda_qvalues([3, 2, 1], [2], pi=0.5), [1 / 3, 1 / 3, 1 / 3])


def test_default_is_plus_one_with_pi_one() -> None:
    # FDR = [1/1, 2/2, 2/3] -> q = [2/3, 2/3, 2/3]
    close(tda_qvalues([3, 2, 1], [2]), [2 / 3, 2 / 3, 2 / 3])


def test_output_is_in_input_order_never_sorted() -> None:
    close(tda_qvalues([1, 3, 2], [2], plus_one=False), [1 / 3, 0.0, 1 / 3])


def test_tied_targets_get_equal_q() -> None:
    # t=2: (0+1)/2; t=1.5: (1+1)/2; t=1: (1+1)/3
    q = tda_qvalues([2, 1, 2], [1.5])
    close(q, [0.5, 2 / 3, 0.5])
    assert q[0] == q[2]


def test_decoy_tied_with_target_is_counted_at_that_threshold() -> None:
    # Q1 default: >= on both, so the decoy at 2 counts against the target at 2
    close(tda_qvalues([2], [2], plus_one=False), [1.0])
    close(tda_qvalues([3, 2], [2], plus_one=False), [0.0, 0.5])


def test_q_clipped_to_one() -> None:
    close(tda_qvalues([1], [5, 4, 3]), [1.0])


# ---------------------------------------------------------------------------- tda_qvalues: edges


def test_nan_targets_get_nan_and_do_not_count() -> None:
    # valid targets [3, 2], decoy 2.5: FDR(3)=0, FDR(2.5)=1/1, FDR(2)=1/2
    close(tda_qvalues([3, NAN, 2], [NAN, 2.5], plus_one=False), [0.0, NAN, 0.5])


def test_minus_inf_target_is_absent() -> None:
    close(tda_qvalues([3, -INF, 2], [], plus_one=False), [0.0, NAN, 0.0])
    close(tda_qvalues([-INF, 2], []), [NAN, 1.0])


def test_minus_inf_decoy_is_absent() -> None:
    # if -inf were a decoy at some threshold it could only matter below every real target
    close(tda_qvalues([1, 2], [-INF, NAN], plus_one=False), [0.0, 0.0])


def test_plus_inf_scores_allowed() -> None:
    # t=inf: (0+1)/1; t=1: (0+1)/2; t=0: (1+1)/2
    close(tda_qvalues([INF, 1], [0]), [0.5, 0.5])
    # +inf decoy tied with +inf target counts
    close(tda_qvalues([INF, 1], [INF], plus_one=False), [0.5, 0.5])


def test_empty_targets() -> None:
    q = tda_qvalues([], [1.0, 2.0])
    assert q.shape == (0,)
    assert q.dtype == np.float64


def test_all_nan_targets() -> None:
    q = tda_qvalues([NAN, NAN], [1.0])
    assert q.shape == (2,)
    assert np.isnan(q).all()


@pytest.mark.parametrize("pi", [0.0, -0.1, 1.0000001, 2.0, NAN])
def test_pi_outside_zero_one_raises(pi: float) -> None:
    with pytest.raises(ValueError):
        tda_qvalues([3, 2, 1], [2], pi=pi)


def test_pi_one_is_allowed_and_small_pi_scales() -> None:
    q1 = tda_qvalues([3, 2, 1], [2], pi=1.0)
    qs = tda_qvalues([3, 2, 1], [2], pi=1e-3)
    close(qs, q1 * 1e-3)


def test_inputs_not_mutated() -> None:
    t = np.array([3.0, NAN, 1.0])
    d = np.array([2.0, NAN])
    t0, d0 = t.copy(), d.copy()
    tda_qvalues(t, d)
    np.testing.assert_array_equal(t, t0)
    np.testing.assert_array_equal(d, d0)


# ---------------------------------------------------------------------------- picked_fdr


def test_example8_picked_competition() -> None:
    res = picked_fdr([5, 4, 3, 3, 2, 1], [1, NAN, 3.6, NAN, 2, NAN])
    assert isinstance(res, PickedResult)
    close(res.q, [0.5, 0.5, NAN, 2 / 3, NAN, 0.75])
    np.testing.assert_array_equal(res.decoy_wins, [False, False, True, False, True, False])
    assert res.decoy_wins.dtype == np.bool_
    assert res.q.dtype == np.float64
    assert res.n_target_wins == 4
    assert res.n_decoy_wins == 2
    assert isinstance(res.n_target_wins, int)
    assert isinstance(res.n_decoy_wins, int)


def test_example9_decoy_only_gene_and_clipping() -> None:
    res = picked_fdr([NAN, 2], [1, NAN])
    np.testing.assert_array_equal(res.decoy_wins, [True, False])
    close(res.q, [NAN, 1.0])
    assert res.n_target_wins == 1
    assert res.n_decoy_wins == 1


def test_picked_tie_goes_to_decoy() -> None:
    res = picked_fdr([2.0, 3.0], [2.0, 1.0])
    np.testing.assert_array_equal(res.decoy_wins, [True, False])
    assert np.isnan(res.q[0])
    # picked list 3T, 2D: FDR(3) = 1/1, FDR(2) = 2/1 -> q = 1
    close(res.q[1:], [1.0])


def test_picked_plus_inf_tie_goes_to_decoy() -> None:
    res = picked_fdr([INF], [INF])
    np.testing.assert_array_equal(res.decoy_wins, [True])
    assert np.isnan(res.q[0])


def test_picked_neither_exists() -> None:
    res = picked_fdr([NAN, -INF, 3.0], [-INF, NAN, NAN])
    np.testing.assert_array_equal(res.decoy_wins, [False, False, False])
    close(res.q, [NAN, NAN, 1.0])
    assert res.n_target_wins == 1
    assert res.n_decoy_wins == 0


def test_picked_all_nan() -> None:
    res = picked_fdr([NAN, NAN], [NAN, NAN])
    assert np.isnan(res.q).all()
    assert not res.decoy_wins.any()
    assert res.n_target_wins == 0
    assert res.n_decoy_wins == 0


def test_picked_empty() -> None:
    res = picked_fdr([], [])
    assert res.q.shape == (0,)
    assert res.decoy_wins.shape == (0,)
    assert res.n_target_wins == 0
    assert res.n_decoy_wins == 0


def test_picked_mismatched_lengths_raise() -> None:
    with pytest.raises(ValueError):
        picked_fdr([1.0, 2.0], [1.0])


def test_picked_output_in_input_order() -> None:
    perm = [3, 0, 5, 1, 4, 2]
    t = np.array([5, 4, 3, 3, 2, 1], dtype=float)
    d = np.array([1, NAN, 3.6, NAN, 2, NAN])
    res = picked_fdr(t[perm], d[perm])
    close(res.q, np.array([0.5, 0.5, NAN, 2 / 3, NAN, 0.75])[perm])


def test_picked_result_is_frozen_with_slots() -> None:
    res = picked_fdr([1.0], [NAN])
    assert dataclasses.is_dataclass(res)
    with pytest.raises(dataclasses.FrozenInstanceError):
        res.n_target_wins = 5  # type: ignore[misc]
    assert hasattr(PickedResult, "__slots__")
    assert {f.name for f in dataclasses.fields(PickedResult)} == {
        "decoy_wins",
        "q",
        "n_target_wins",
        "n_decoy_wins",
    }
