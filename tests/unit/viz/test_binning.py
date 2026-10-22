"""``giap.viz.base.binned``: the 20,000-element rule per track (M9-01 behaviour 8, examples 5-6,
index decision 12).

``binned`` returns ``None`` at or under ``max_elements`` elements, otherwise the
``min(max_elements, length) + 1`` integer edges of near-equal bins (widths differ by at most 1 bp,
PR #191 B2) over ``[0, length)``.
"""

from __future__ import annotations

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.viz.base import binned


def _intervals(n: int, length: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    starts = rng.integers(0, length, size=n)
    ends = np.minimum(starts + rng.integers(1, 3000, size=n), length)
    return starts, ends


def test_example_5_a_5_mb_genome_with_4500_genes_is_not_binned() -> None:
    starts, ends = _intervals(4_500, 5_000_000, 1)
    assert binned(starts, ends, 5_000_000) is None


def test_exactly_max_elements_is_not_binned() -> None:
    starts, ends = _intervals(20_000, 1_000_000, 2)
    assert binned(starts, ends, 1_000_000) is None


def test_example_6_25000_intervals_give_20000_bins_of_50_bp() -> None:
    starts, ends = _intervals(25_000, 1_000_000, 3)
    edges = binned(starts, ends, 1_000_000)
    assert edges is not None
    assert len(edges) == 20_001
    assert edges[0] == 0
    assert edges[-1] == 1_000_000
    np.testing.assert_allclose(np.diff(edges), 50.0)


def test_custom_max_elements() -> None:
    starts, ends = _intervals(11, 1_000, 4)
    edges = binned(starts, ends, 1_000, max_elements=10)
    assert edges is not None
    np.testing.assert_allclose(edges, np.linspace(0, 1_000, 11))


@settings(max_examples=40, deadline=None)
@given(
    max_elements=st.integers(1, 20_000),
    extra=st.integers(-50, 50),
    length=st.integers(1_000, 10_000_000),
    seed=st.integers(0, 2**32 - 1),
)
def test_property_bins_at_most_max_and_cover_every_element(
    max_elements: int, extra: int, length: int, seed: int
) -> None:
    n = max(0, max_elements + extra)
    starts, ends = _intervals(n, length, seed)
    edges = binned(starts, ends, length, max_elements=max_elements)
    if n <= max_elements:
        assert edges is None
        return
    assert edges is not None
    assert len(edges) - 1 <= 20_000
    assert len(edges) - 1 == min(max_elements, length)
    assert edges[0] == 0 and edges[-1] == length
    assert np.all(np.diff(edges) > 0)
    w = np.diff(edges)
    assert w.max() - w.min() <= 1  # integer edges (PR #191 B2)
    # every element falls into a bin: its start and its last base are inside [0, length)
    for pos in (starts, ends - 1):
        idx = np.searchsorted(edges, pos, side="right") - 1
        assert np.all((idx >= 0) & (idx < len(edges) - 1))
