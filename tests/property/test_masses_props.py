"""Hypothesis properties of the generated mass-equivalence list (issue M3-05; GIPPA_SPEC.md
Appendix B.5, §7.2.4).

- brute-force oracle: at any tolerance the list equals the irreducible brute-force equivalences
  (1<->1, 1<->2, 2<->2 with no shared residue, residue+mod <-> residue) up to the tolerance
  boundary; 2<->2 compositions sharing a residue are never listed (lead decision 2, PR #35);
- monotone in the tolerance: a smaller tolerance gives a sub-list;
- every entry's delta is the recomputed |m(left) - m(right)| and lies within the tolerance;
- the list is sorted and identical on repeated calls.
"""

from __future__ import annotations

from itertools import combinations, combinations_with_replacement

from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.mapping.masses import MOD_DELTA, RESIDUE_MASS, mass_equivalents

FAST = settings(max_examples=60, deadline=None)
TOL = st.floats(0.0, 0.06, allow_nan=False)
STANDARD = "ACDEFGHKLMNPQRSTVWY"
EPS = 1e-6  # slack at the tolerance boundary: delta_da is rounded to 6 decimals (decision 4)


def _canon(side: str) -> str:
    return side if "+" in side else "".join(sorted(side))


def _mass(side: str) -> float:
    if "+" in side:
        res, mod = side.split("+", 1)
        return RESIDUE_MASS[res] + MOD_DELTA[mod]
    return sum(RESIDUE_MASS[c] for c in side)


def _keys(tol: float) -> set[frozenset[str]]:
    return {frozenset({_canon(e.left), _canon(e.right)}) for e in mass_equivalents(tol)}


def _brute(tol: float) -> tuple[set[frozenset[str]], set[frozenset[str]], set[frozenset[str]]]:
    """(certain, borderline, excluded) brute-force keys at ``tol``."""
    comps = list(STANDARD) + ["".join(c) for c in combinations_with_replacement(STANDARD, 2)]
    certain: set[frozenset[str]] = set()
    border: set[frozenset[str]] = set()
    excluded: set[frozenset[str]] = set()

    def put(a: str, b: str, *, shared: bool = False) -> None:
        d = abs(_mass(a) - _mass(b))
        key = frozenset({a, b})
        if shared:
            excluded.add(key)
        elif abs(d - tol) <= EPS:
            border.add(key)
        elif d < tol:
            certain.add(key)

    for a, b in combinations(comps, 2):
        if len(a) == 1 or len(b) == 1:
            put(a, b)
        else:
            put(a, b, shared=bool(set(a) & set(b)))
    for r in STANDARD:
        for mod in MOD_DELTA:
            for s in STANDARD:
                put(f"{r}+{mod}", s)
    return certain, border, excluded


@FAST
@given(TOL)
def test_matches_brute_force_oracle(tol: float) -> None:
    got = _keys(tol)
    certain, border, excluded = _brute(tol)
    missing = certain - got
    assert not missing, sorted(sorted(k) for k in missing)[:5]
    shared = got & excluded
    assert not shared, sorted(sorted(k) for k in shared)[:5]
    extra = got - certain - border
    assert not extra, sorted(sorted(k) for k in extra)[:5]


@FAST
@given(TOL, TOL)
def test_monotone_in_tolerance(t1: float, t2: float) -> None:
    lo, hi = sorted((t1, t2))
    assert _keys(lo) <= _keys(hi)


@FAST
@given(TOL)
def test_entries_valid_sorted_and_repeatable(tol: float) -> None:
    entries = mass_equivalents(tol)
    for e in entries:
        assert e.delta_da == abs(_mass(e.left) - _mass(e.right)) or (
            abs(e.delta_da - abs(_mass(e.left) - _mass(e.right))) <= EPS
        )
        assert e.delta_da <= tol + EPS
    keys = [(e.kind, e.delta_da, e.left, e.right) for e in entries]
    assert keys == sorted(keys)
    assert mass_equivalents(tol) == entries
