"""scripts/ci/shard_tests.py: the CI coverage shards partition the pytest files."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("shard_tests", ROOT / "scripts/ci/shard_tests.py")
assert _spec and _spec.loader
st = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(st)


@pytest.mark.parametrize("total", [1, 4, 6, 7])
def test_shards_partition_the_real_test_files(total: int) -> None:
    files = st.test_files(ROOT)
    parts = [st.shard(files, i, total) for i in range(1, total + 1)]
    flat = [f for p in parts for f in p]
    assert sorted(flat) == files and len(flat) == len(set(flat))
    assert max(map(len, parts)) - min(map(len, parts)) <= 1


def test_this_file_and_calibration_files_are_included() -> None:
    files = st.test_files(ROOT)
    assert "tests/unit/test_ci_shard_tests.py" in files
    assert any(f.startswith("tests/calibration/") for f in files)
    assert not any("__pycache__" in f for f in files)


def test_bad_index_raises() -> None:
    with pytest.raises(ValueError):
        st.shard(["a"], 0, 4)
    with pytest.raises(ValueError):
        st.shard(["a"], 5, 4)
