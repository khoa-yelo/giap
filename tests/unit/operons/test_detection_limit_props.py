"""Hypothesis properties of the M6-02 detection-limit model
(docs/specs/M6-02-detection-limit-model.md, acceptance criteria "Property: ...").

- p_detect equals the closed-form Poisson-binomial (2^n enumeration) of the per-peptide logits;
- p_detect is non-decreasing in a (abundance coefficient >= 0) and non-increasing in k_min;
- p_miss = mean_r [1 - p_detect(a_v + r)], lies in [0, 1] and is non-increasing in a_v;
- p_detect / p_miss are deterministic and do not mutate their inputs;
- observable_peptides agrees with M4-10's n_observable_peptides count.
"""

from __future__ import annotations

import dataclasses
import itertools
import math
from typing import Any

import numpy as np
import pytest
from helpers import detection_limit_data as D
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.config import GippaConfig

SETTINGS = settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])


def _dl() -> Any:
    from gippa.operons import detection_limit

    return detection_limit


coefs = st.fixed_dictionaries(
    {
        "intercept": st.floats(-20, 5),
        "abundance_log2": st.floats(0.0, 1.5),
        "length": st.floats(-0.2, 0.2),
        "gravy": st.floats(-2.0, 2.0),
    }
)
peps_st = st.lists(st.tuples(st.integers(7, 30), st.floats(-3.0, 3.0)), min_size=1, max_size=8)
abund = st.floats(0.0, 35.0)
ratios_st = st.lists(st.floats(-4.0, 4.0), min_size=1, max_size=6)


def _model(coef: dict[str, float], k_min: int = 1) -> Any:
    return _dl().PMissModel(kind="peptide", coef=coef, k_min=k_min, n_fit=500, auc=0.8, usable=True)


@SETTINGS
@given(coefs, abund, peps_st, st.integers(1, 3))
def test_p_detect_is_the_poisson_binomial_of_the_logits(
    coef: dict[str, float], a: float, peps: list[tuple[int, float]], k: int
) -> None:
    got = _dl().p_detect(_model(coef, k), a, peps)
    assert got == pytest.approx(D.p_detect_oracle(coef, a, peps, k), rel=1e-9, abs=1e-12)
    assert 0.0 <= got <= 1.0


@SETTINGS
@given(coefs, abund, abund, peps_st, st.integers(1, 3))
def test_p_detect_non_decreasing_in_a(
    coef: dict[str, float], a1: float, a2: float, peps: list[tuple[int, float]], k: int
) -> None:
    lo, hi = sorted((a1, a2))
    m = _model(coef, k)
    assert _dl().p_detect(m, lo, peps) <= _dl().p_detect(m, hi, peps) + 1e-12


@SETTINGS
@given(coefs, abund, peps_st)
def test_p_detect_non_increasing_in_k_min(
    coef: dict[str, float], a: float, peps: list[tuple[int, float]]
) -> None:
    m = _model(coef, 1)
    vals = [
        _dl().p_detect(dataclasses.replace(m, k_min=k), a, peps) for k in range(1, len(peps) + 2)
    ]
    assert all(x >= y - 1e-12 for x, y in itertools.pairwise(vals))
    assert vals[-1] == 0.0  # k_min > number of observable peptides: never called


@SETTINGS
@given(coefs, abund, ratios_st, peps_st, st.integers(1, 2))
def test_p_miss_is_mean_of_one_minus_p_detect_in_unit_interval(
    coef: dict[str, float], a_v: float, ratios: list[float], peps: list[tuple[int, float]], k: int
) -> None:
    dl = _dl()
    m = _model(coef, k)
    pm = dl.p_miss(m, a_v, ratios, peps)
    want = float(np.mean([1.0 - dl.p_detect(m, a_v + r, peps) for r in ratios]))
    assert pm == pytest.approx(want, rel=1e-12, abs=1e-15)
    assert 0.0 <= pm <= 1.0


@SETTINGS
@given(coefs, abund, abund, ratios_st, peps_st, st.integers(1, 2))
def test_p_miss_non_increasing_in_a_v(
    coef: dict[str, float],
    a1: float,
    a2: float,
    ratios: list[float],
    peps: list[tuple[int, float]],
    k: int,
) -> None:
    lo, hi = sorted((a1, a2))
    m = _model(coef, k)
    assert _dl().p_miss(m, hi, ratios, peps) <= _dl().p_miss(m, lo, ratios, peps) + 1e-12


@SETTINGS
@given(coefs, abund, ratios_st, peps_st)
def test_p_miss_deterministic_and_inputs_untouched(
    coef: dict[str, float], a_v: float, ratios: list[float], peps: list[tuple[int, float]]
) -> None:
    c0, r0, p0 = dict(coef), list(ratios), list(peps)
    m = _model(coef)
    arr = np.asarray(ratios, dtype=float)
    x = _dl().p_miss(m, a_v, arr, peps)
    y = _dl().p_miss(m, a_v, arr, peps)
    assert x == y or (math.isnan(x) and math.isnan(y))
    assert dict(m.coef) == c0 and coef == c0 and ratios == r0 and peps == p0
    np.testing.assert_array_equal(arr, np.asarray(r0, dtype=float))


@SETTINGS
@given(coefs, abund, ratios_st)
def test_p_miss_order_of_ratios_irrelevant(
    coef: dict[str, float], a_v: float, ratios: list[float]
) -> None:
    m = _model(coef)
    peps = list(D.TOY_PEPS)
    assert _dl().p_miss(m, a_v, ratios, peps) == pytest.approx(
        _dl().p_miss(m, a_v, ratios[::-1], peps), rel=1e-12, abs=1e-15
    )


aa_st = st.text(alphabet="ACDEFGHKMNPQRSTVWY", min_size=0, max_size=300)


@SETTINGS
@given(aa_st, st.sampled_from(["trypsin_p", "trypsin"]))
def test_observable_peptides_count_matches_m4_10(aa: str, enzyme: str) -> None:
    from gippa.genetable.covariates import n_observable_peptides

    cfg = GippaConfig(abundance={"enzyme": enzyme})
    got = _dl().observable_peptides(aa, cfg)
    assert len(got) == n_observable_peptides(aa, enzyme, (7, 30))  # type: ignore[arg-type]
    assert all(7 <= n <= 30 for n, _ in got)
    assert all(-4.5 <= g <= 4.5 for _, g in got)
