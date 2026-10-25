"""``scripts/pilot/cascadia_shards.py``: the shard blocks of the CPU array (lead, 2026-10-24).

The blocks must partition the batches in order, so merging them reproduces the unsplit run's batch
sequence (and row numbers in the .ssl). Cascadia itself is not imported here.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "pilot" / "cascadia_shards.py"


def _load():  # type: ignore[no-untyped-def]
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location("cascadia_shards", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S = _load()


@given(st.integers(0, 20_000), st.integers(1, 64))
def test_blocks_partition_the_batches_in_order(n_batches: int, n_shards: int) -> None:
    blocks = [S.shard_batches(n_batches, k, n_shards) for k in range(n_shards)]
    assert [b for blk in blocks for b in blk] == list(range(n_batches))
    sizes = [len(b) for b in blocks]
    assert max(sizes) - min(sizes) <= 1


def test_example_4471_batches_in_16_shards() -> None:
    assert S.shard_batches(4471, 0, 16) == range(0, 280)
    assert S.shard_batches(4471, 15, 16) == range(4192, 4471)


@pytest.mark.parametrize("shard", [-1, 16])
def test_shard_out_of_range_raises(shard: int) -> None:
    with pytest.raises(ValueError):
        S.shard_batches(4471, shard, 16)


@given(st.integers(1, 5000), st.integers(1, 64), st.data())
def test_chunk_plan_tiles_everything_outside_the_reserved_blocks(
    n_batches: int, size: int, data: st.DataObject
) -> None:
    cuts = sorted(data.draw(st.lists(st.integers(0, n_batches), max_size=4)))
    reserved = [(cuts[i], cuts[i + 1]) for i in range(0, len(cuts) - 1, 2)]
    plan = S.chunk_plan(n_batches, reserved, size)
    got = [b for blk in plan for b in blk]
    want = [b for b in range(n_batches) if not any(a <= b < c for a, c in reserved)]
    assert got == want
    assert all(0 < len(blk) <= size for blk in plan)


def test_chunk_plan_example() -> None:
    plan = S.chunk_plan(100, [(0, 46)], 20)
    assert [(b.start, b.stop) for b in plan] == [(46, 66), (66, 86), (86, 100)]
