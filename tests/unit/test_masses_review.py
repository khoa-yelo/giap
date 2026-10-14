"""Pins the lead decisions after code review of M3-05 (PR #35;
docs/specs/M3-05-masses-and-diff-ops.md, items 1-4): inclusive tolerance, shared-residue 2<->2 pairs
excluded, 255 padding in
``BlockTables``, and ``delta_da`` rounded to 6 decimals before sorting and output.
"""

from __future__ import annotations

import numpy as np
import pytest

from giap.mapping.diffops import block_tables
from giap.mapping.masses import RESIDUE_MASS, mass_equivalents


def _pairs(tol: float) -> set[frozenset[str]]:
    return {frozenset({e.left, e.right}) for e in mass_equivalents(tol)}


# (1) tolerance is inclusive: Δ(K,Q) exactly equal to tol_da is listed.
@pytest.mark.parametrize("rounded", [False, True], ids=["raw", "rounded6"])
def test_tolerance_inclusive_at_delta_k_q(rounded: bool) -> None:
    delta = abs(RESIDUE_MASS["K"] - RESIDUE_MASS["Q"])
    tol = round(delta, 6) if rounded else delta
    assert frozenset({"K", "Q"}) in _pairs(tol)


def test_tolerance_just_below_delta_k_q_excludes() -> None:
    delta = round(abs(RESIDUE_MASS["K"] - RESIDUE_MASS["Q"]), 6)
    assert frozenset({"K", "Q"}) not in _pairs(delta - 1e-6)


# (2) 2<->2 pairs sharing a residue are never listed.
@pytest.mark.parametrize("tol", [0.05, 0.001])
def test_shared_residue_two_two_never_listed(tol: float) -> None:
    entries = mass_equivalents(tol)
    assert frozenset({"AK", "AQ"}) not in {frozenset({e.left, e.right}) for e in entries}
    shared = [
        (e.left, e.right)
        for e in entries
        if "+" not in e.left + e.right
        and len(e.left) == len(e.right) == 2
        and set(e.left) & set(e.right)
    ]
    assert not shared, shared[:5]


# (3) padded slots of BlockTables hold 255, never a residue code.
def test_block_tables_padding_is_255() -> None:
    t = block_tables()
    n_pad = 0
    for arr, lens in ((t.left, t.left_len), (t.right, t.right_len)):
        slots = np.arange(arr.shape[1])[None, :]
        pad = slots >= lens[:, None]
        n_pad += int(pad.sum())
        assert np.all(arr[pad] == 255)
        assert np.all(arr[~pad] != 255)
    assert n_pad > 0  # length-1 sides exist (N, Q, W, R)


# (4) delta_da rounded to 6 decimals before sorting and output.
@pytest.mark.parametrize("tol", [0.05, 0.001, 0.06])
def test_delta_rounded_and_order_stable(tol: float) -> None:
    entries = mass_equivalents(tol)
    assert entries
    for e in entries:
        assert e.delta_da == round(e.delta_da, 6), e
    keys = [(e.kind, round(e.delta_da, 6), e.left, e.right) for e in entries]
    assert keys == sorted(keys)
