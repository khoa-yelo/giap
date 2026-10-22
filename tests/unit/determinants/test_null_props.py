"""M8-09 Hypothesis properties of the TU-size-preserving permutation (spec acceptance criterion 3):
``permute_response`` preserves the multiset of TU sizes, each TU's internal order, and fixes every
TU whose size is unique; ``tu_size_blocks`` partitions the rows by TU size."""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_attrib_data as a


@st.composite
def tu_layouts(draw: st.DrawFn) -> tuple[list[str], int]:
    """TU sizes 1-5, rows of each TU possibly interleaved (shuffled row order), and a seed."""
    sizes = draw(st.lists(st.integers(1, 5), min_size=1, max_size=25))
    tu = a.tu_from_sizes(sizes, [f"T{i:03d}" for i in range(len(sizes))])
    if draw(st.booleans()):
        order = draw(st.permutations(range(len(tu))))
        tu = [tu[i] for i in order]
    return tu, draw(st.integers(0, 2**32 - 1))


def _rows_by_tu(tu: list[str]) -> dict[str, list[int]]:
    out: dict[str, list[int]] = {}
    for i, t in enumerate(tu):
        out.setdefault(t, []).append(i)
    return out


@settings(max_examples=200, deadline=None)
@given(tu_layouts())
def test_permute_response_moves_whole_blocks_among_equal_sizes(
    layout: tuple[list[str], int],
) -> None:
    tu, seed = layout
    y = np.arange(len(tu), dtype=float) + 0.5  # distinct values identify the source gene
    out = np.asarray(a.null_mod().permute_response(y, tu, np.random.default_rng(seed)))
    assert out.shape == y.shape
    assert sorted(out.tolist()) == sorted(y.tolist())  # a permutation of y
    rows = _rows_by_tu(tu)
    size_count = Counter(len(r) for r in rows.values())
    source_of: dict[str, str] = {}
    for t, r in rows.items():
        vals = out[r].tolist()
        src = [s for s, rr in rows.items() if y[rr].tolist() == vals]
        assert len(src) == 1, (t, vals)  # the block is one TU's vector, in its gene order
        assert len(rows[src[0]]) == len(r)  # same size
        source_of[t] = src[0]
        if size_count[len(r)] == 1:
            assert src[0] == t  # unique size: fixed
    assert sorted(source_of.values()) == sorted(rows)  # bijection among TUs
    # the multiset of TU sizes seen by each response block is unchanged
    assert Counter(len(rows[s]) for s in source_of.values()) == size_count


@settings(max_examples=200, deadline=None)
@given(tu_layouts())
def test_tu_size_blocks_partition(layout: tuple[list[str], int]) -> None:
    tu, _ = layout
    blocks = a.null_mod().tu_size_blocks(tu)
    rows = _rows_by_tu(tu)
    sizes = [np.asarray(b).shape[1] for b in blocks]
    assert sizes == sorted(set(len(r) for r in rows.values()))  # one array per size, ascending
    flat = np.concatenate([np.asarray(b).ravel() for b in blocks])
    assert sorted(flat.tolist()) == list(range(len(tu)))
    for b in blocks:
        b = np.asarray(b)
        assert b.dtype == np.int64
        tus = [tu[int(r[0])] for r in b]
        assert tus == sorted(tus)  # TUs sorted by ID
        for r in b:
            assert r.tolist() == rows[tu[int(r[0])]]  # the TU's rows, in input order
