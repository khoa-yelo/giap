"""Hypothesis properties of placement q-values (issue M3-10; GIAP_SPEC.md §7.2.5, §8.3).

- ``stratum_qvalues``: q in [0, 1]; q non-decreasing in E; adding a decoy never lowers any q;
  invariant to input order; equals M4-01 ``tda_qvalues(-E_T, -E_D, pi)`` for constant pi;
  invariant to a strictly increasing transform of E (E is a ranking score only);
- ``placement``: q in [0, 1] (0 for ``same_db`` exact), adding a decoy-only peptide never lowers
  any q, row order never changes the result.
"""

from __future__ import annotations

import dataclasses
from typing import Any

import numpy as np
import pandas as pd
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.config import GiapConfig
from giap.mapping.placement_fdr import PeptideCosts, placement, stratum_qvalues
from giap.stats.tda import tda_qvalues

FAST = settings(max_examples=150, deadline=None)
EPS = 1e-12
# small grid makes ties between and among targets and decoys frequent
E = st.one_of(st.integers(0, 8).map(lambda k: k / 4), st.floats(1e-9, 10.0, allow_nan=False))
PI = st.one_of(st.just(1.0), st.floats(0.001, 1.0, allow_nan=False))
ALPHABET19 = "ACDEFGHKLMNPQRSTVWY"


@st.composite
def stratum(draw: st.DrawFn) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    et = draw(st.lists(E, min_size=1, max_size=40))
    ed = draw(st.lists(E, max_size=40))
    pi = draw(st.lists(PI, min_size=len(et), max_size=len(et)))
    return np.array(et), np.array(ed, dtype=float), np.array(pi)


@FAST
@given(stratum())
def test_q_in_unit_interval(s: tuple[np.ndarray, np.ndarray, np.ndarray]) -> None:
    q = stratum_qvalues(*s)
    assert q.shape == s[0].shape
    assert np.all((q >= 0.0) & (q <= 1.0))


@FAST
@given(stratum())
def test_q_monotone_in_e(s: tuple[np.ndarray, np.ndarray, np.ndarray]) -> None:
    et = s[0]
    q = stratum_qvalues(*s)
    order = np.argsort(et, kind="stable")
    assert np.all(np.diff(q[order]) >= -EPS)
    # tied E -> identical q
    for v in np.unique(et):
        assert np.ptp(q[et == v]) <= EPS


@FAST
@given(stratum(), E)
def test_adding_a_decoy_never_lowers_q(
    s: tuple[np.ndarray, np.ndarray, np.ndarray], extra: float
) -> None:
    et, ed, pi = s
    before = stratum_qvalues(et, ed, pi)
    after = stratum_qvalues(et, np.append(ed, extra), pi)
    assert np.all(after >= before - EPS)


@FAST
@given(stratum(), st.randoms(use_true_random=False))
def test_input_order_invariant(s: tuple[np.ndarray, np.ndarray, np.ndarray], rnd: Any) -> None:
    et, ed, pi = s
    pt = list(range(et.size))
    pd_ = list(range(ed.size))
    rnd.shuffle(pt)
    rnd.shuffle(pd_)
    q = stratum_qvalues(et, ed, pi)
    q2 = stratum_qvalues(et[pt], ed[pd_], pi[pt])
    np.testing.assert_allclose(q2, q[pt], rtol=1e-12, atol=EPS)


@FAST
@given(stratum(), PI)
def test_constant_pi_equals_tda_qvalues(
    s: tuple[np.ndarray, np.ndarray, np.ndarray], pi: float
) -> None:
    et, ed, _ = s
    q = stratum_qvalues(et, ed, np.full(et.size, pi))
    np.testing.assert_allclose(q, tda_qvalues(-et, -ed, pi=pi, plus_one=True), rtol=1e-12)


@FAST
@given(stratum(), st.sampled_from([2.0**-7, 2.0**10]))
def test_e_is_a_ranking_score_only(
    s: tuple[np.ndarray, np.ndarray, np.ndarray], scale: float
) -> None:
    """Lead decision (M3-09 review): E only ranks; rescaling target and decoy E alike (exactly,
    by a power of two) leaves every q unchanged (no absolute reading of E)."""
    et, ed, pi = s

    def f(x: np.ndarray) -> np.ndarray:
        return scale * x

    np.testing.assert_allclose(
        stratum_qvalues(f(et), f(ed), pi), stratum_qvalues(et, ed, pi), rtol=1e-12, atol=EPS
    )


# ----------------------------------------------------------------------------------- placement


def _pep(i: int, length: int) -> str:
    digits = "".join(ALPHABET19[(i // 19**k) % 19] for k in (2, 1, 0))
    return ("MDEFGHKLNPQRSTVWYAC" * 3)[: length - 3] + digits


@st.composite
def costs(draw: st.DrawFn) -> list[dict[str, Any]]:
    n = draw(st.integers(0, 30))
    rows = []
    for i in range(n):
        src = draw(st.sampled_from(["database", "denovo"]))
        has_t = draw(st.booleans())
        has_d = draw(st.booleans())
        rows.append(
            {
                "peptide": _pep(i, draw(st.integers(7, 40))),
                "source": src,
                "search_db": draw(st.sampled_from(["same", "other"]))
                if src == "database"
                else "other",
                "m_target": draw(st.integers(0, 2)) if has_t else None,
                "e_target": draw(E) if has_t else None,
                "m_decoy": draw(st.integers(0, 2)) if has_d else None,
                "e_decoy": draw(E) if has_d else None,
                "pep": draw(st.one_of(st.none(), st.floats(1e-4, 0.05))),
                "q_value": draw(st.one_of(st.none(), st.floats(1e-4, 0.01))),
                "pi_denovo": draw(st.one_of(st.none(), PI)) if src == "denovo" else None,
            }
        )
    return rows


def _frame(rows: list[dict[str, Any]]) -> pd.DataFrame:
    out = [
        dataclasses.asdict(PeptideCosts(length=len(r["peptide"]), **r))  # type: ignore[arg-type]
        for r in rows
    ]
    return pd.DataFrame(out, columns=[f.name for f in dataclasses.fields(PeptideCosts)])


def _q(rows: list[dict[str, Any]]) -> dict[str, float]:
    per, _ = placement(_frame(rows), GiapConfig())
    return {
        str(p): float(q)
        for p, q in zip(per["peptide"], per["placement_q"], strict=True)
        if not pd.isna(q)
    }


@FAST
@given(costs())
def test_placement_q_in_unit_interval(rows: list[dict[str, Any]]) -> None:
    q = _q(rows)
    targets = {r["peptide"] for r in rows if r["m_target"] is not None}
    assert set(q) == targets
    for r in rows:
        if r["peptide"] in q:
            assert 0.0 <= q[r["peptide"]] <= 1.0
            if r["search_db"] == "same" and r["m_target"] == 0:
                assert q[r["peptide"]] == 0.0


@FAST
@given(costs(), st.integers(0, 2), E, st.sampled_from(["database", "denovo"]))
def test_placement_adding_a_decoy_never_lowers_q(
    rows: list[dict[str, Any]], m: int, e: float, src: str
) -> None:
    before = _q(rows)
    extra = {
        "peptide": _pep(5000, 9), "source": src, "search_db": "other", "m_target": None,
        "e_target": None, "m_decoy": m, "e_decoy": e, "pep": None, "q_value": None,
        "pi_denovo": None,
    }  # fmt: skip
    after = _q([*rows, extra])
    for p, q in before.items():
        assert after[p] >= q - EPS


@FAST
@given(costs(), st.randoms(use_true_random=False))
def test_placement_row_order_invariant(rows: list[dict[str, Any]], rnd: Any) -> None:
    shuffled = list(rows)
    rnd.shuffle(shuffled)
    assert _q(rows) == _q(shuffled)
