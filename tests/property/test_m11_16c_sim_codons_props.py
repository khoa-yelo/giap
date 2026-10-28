"""Hypothesis properties of gippasim's per-table codon usage (M11-16c; docs/specs/
M11-16-translation-tables-4-25.md item 8).

For every ``tga_sense_share`` in (0, 1) and table 4 or 25: ``sense_usage`` has no stop of the
table, every weight is positive, TGA's share of its family is exactly the share, every amino
acid keeps its table-11 total, and the total weight is unchanged. ``SimConfig`` accepts every
such share and refuses every share outside (0, 1).
"""

from __future__ import annotations

import math
from collections import Counter

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

FAMILY = {4: "W", 25: "G"}
tables = st.sampled_from([4, 25])
shares = st.floats(min_value=1e-6, max_value=1 - 1e-6, allow_nan=False)
outside = st.one_of(
    st.floats(max_value=0.0, allow_nan=False, allow_infinity=False),
    st.floats(min_value=1.0, allow_nan=False, allow_infinity=False),
)


def _totals(usage: dict[str, float], code: dict[str, str]) -> Counter[str]:
    out: Counter[str] = Counter()
    for c, w in usage.items():
        out[code[c]] += w
    return out


@settings(deadline=None)
@given(tables, shares)
def test_sense_usage_invariants(table: int, share: float) -> None:
    from gippasim.codons import CODE_11, SENSE_USAGE, code, sense_usage, stops

    u = sense_usage(table, tga_sense_share=share)
    t = code(table)
    assert not set(u) & stops(table)
    assert "TGA" in u and all(w > 0 for w in u.values())
    fam = sum(w for c, w in u.items() if t[c] == FAMILY[table])
    assert math.isclose(u["TGA"] / fam, share, rel_tol=1e-9)
    want, got = _totals(SENSE_USAGE, CODE_11), _totals(u, t)
    assert set(got) == set(want)
    for a, w in want.items():
        assert math.isclose(got[a], w, rel_tol=1e-9)
    assert math.isclose(sum(u.values()), sum(SENSE_USAGE.values()), rel_tol=1e-9)


@settings(deadline=None)
@given(tables, shares)
def test_simconfig_accepts_open_interval(table: int, share: float) -> None:
    from gippasim.simulate import SimConfig

    cfg = SimConfig(genetic_code=table, tga_sense_share=share)
    assert cfg.genetic_code == table and cfg.tga_sense_share == share


@settings(deadline=None)
@given(tables, outside)
def test_simconfig_refuses_outside(table: int, share: float) -> None:
    from gippasim.simulate import SimConfig

    with pytest.raises(ValueError):
        SimConfig(genetic_code=table, tga_sense_share=share)
