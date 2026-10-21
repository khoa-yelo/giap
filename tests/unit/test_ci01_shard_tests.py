"""``scripts/ci/shard_tests.py``: a stable partition of a pytest selection by test id (CI-01).

Lead decision 1: ``scripts/slurm/calibration_full.sbatch`` shards *by test id*, not by file, with
"hash of the node id modulo the shard count". The partition must therefore be

* complete and disjoint — every collected test lands in exactly one shard, so the shards together
  run the whole selection;
* stable — the same test id gives the same shard in every process, on every node and in every run
  (:func:`hash` is salted per process and would not do), because each array task computes the
  partition on its own, with no shared plan file.
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "ci" / "shard_tests.py"

_spec = importlib.util.spec_from_file_location("ci01_shard_tests", SCRIPT)
assert _spec is not None and _spec.loader is not None, f"{SCRIPT} does not exist (CI-01)"
shard_tests = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(shard_tests)

COLLECTED = """tests/calibration/test_engine_calibration.py::test_perm_p_level[full]
tests/calibration/test_engine_calibration.py::test_moran_lag1_level_tu_null[full]
tests/calibration/attacks/test_attack_mayu.py::test_mayu_vs_bh[decoy loss-full]
tests/calibration/test_region_tests.py::test_region_level[full]

161/9019 tests collected (8858 deselected) in 25.15s
"""


def _ids(n: int) -> list[str]:
    return [f"tests/calibration/test_mod{i % 7}.py::test_case{i}[full]" for i in range(n)]


def test_node_ids_keeps_only_test_ids_in_collection_order() -> None:
    ids = shard_tests.node_ids(COLLECTED)
    assert ids == COLLECTED.splitlines()[:4], ids


def test_node_ids_accepts_parameter_ids_with_spaces() -> None:
    assert shard_tests.node_ids("tests/a.py::test_x[decoy loss-full]\n") == [
        "tests/a.py::test_x[decoy loss-full]"
    ]


def test_node_ids_deduplicates() -> None:
    twice = "tests/a.py::test_x\ntests/a.py::test_x\n"
    assert shard_tests.node_ids(twice) == ["tests/a.py::test_x"]


@pytest.mark.parametrize("shards", [1, 2, 7, 8, 16, 40])
def test_the_shards_are_a_partition_of_the_selection(shards: int) -> None:
    ids = _ids(161)
    parts = [shard_tests.shard(ids, k, shards) for k in range(shards)]
    flat = [n for part in parts for n in part]
    assert sorted(flat) == sorted(ids), "a test is missing from the shards or ran twice"
    assert len(set(flat)) == len(flat)
    for k, part in enumerate(parts):
        assert part == [n for n in ids if n in set(part)], f"shard {k} is not in collection order"


def test_shard_of_is_stable_across_processes() -> None:
    """A fresh interpreter (new hash salt) must assign the same shards."""
    ids = _ids(50)
    code = (
        f"import sys; sys.path.insert(0, {str(SCRIPT.parent)!r}); import shard_tests;"
        "print(' '.join(str(shard_tests.shard_of(n, 8)) for n in sys.argv[1:]))"
    )
    env_runs = [
        subprocess.run(
            [sys.executable, "-c", code, *ids], capture_output=True, text=True, check=True
        ).stdout.strip()
        for _ in range(2)
    ]
    here = " ".join(str(shard_tests.shard_of(n, 8)) for n in ids)
    assert env_runs[0] == env_runs[1] == here, env_runs


def test_shards_are_roughly_even() -> None:
    """Sizing the shard count from a total only makes sense if no shard draws most of the tests."""
    counts = shard_tests.counts(_ids(161), 8)
    assert sum(counts) == 161
    assert max(counts) <= 2 * (161 / 8), counts


@pytest.mark.parametrize(("index", "shards"), [(0, 0), (-1, 4), (4, 4), (0, -2)])
def test_out_of_range_shards_are_refused(index: int, shards: int) -> None:
    with pytest.raises(ValueError):
        shard_tests.shard(_ids(3), index, shards)


def _run(args: list[str], stdin: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args], input=stdin, capture_output=True, text=True
    )


def test_cli_reads_stdin_and_prints_one_id_per_line() -> None:
    res = _run(["--shard", "0", "--of", "1"], COLLECTED)
    assert res.returncode == 0, res.stderr
    assert res.stdout.splitlines() == COLLECTED.splitlines()[:4]


def test_cli_counts_every_shard() -> None:
    res = _run(["--of", "4", "--count"], COLLECTED)
    assert res.returncode == 0, res.stderr
    rows = [ln.split("\t") for ln in res.stdout.splitlines()]
    assert [r[0] for r in rows] == ["0", "1", "2", "3"]
    assert sum(int(r[1]) for r in rows) == 4


def test_cli_fails_on_an_empty_collection() -> None:
    res = _run(["--shard", "0", "--of", "4"], "no tests ran\n")
    assert res.returncode == 1 and "no test ids" in res.stderr


def test_cli_refuses_a_shard_outside_the_range() -> None:
    res = _run(["--shard", "4", "--of", "4"], COLLECTED)
    assert res.returncode == 2 and "shard 4" in res.stderr
