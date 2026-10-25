"""``--merge`` folds in the shard exits (review of PR #327, B1 and B2) and the xfail bookkeeping
of ADR-0057 has no gaps a mutant can slip through (should-fix 1, mutants e-g).

B1/B2: ``main --shard K/N`` appends ``{suite, shard, n, block, pytest_exit}`` to
``<suite>.shards.jsonl``. A shard's exit carries every verdict that is not a harness record: shard
0's direct assertions (D-1 quartiles, the mvp04 non-vacuity check, the placement V0 floor ...), a
strict XPASS of a shard-0-only test, the shard-local power checks. So, for every suite, with or
without calibrate records:

- every recorded shard exit 0: the merge verdict is the records' (as before);
- some shard exit 1 (a test failed): ``MERGE_FAIL``; any other non-zero exit (2 interrupted /
  collection error, 3 internal, 4 usage, 5 nothing collected, a signal): ``MERGE_ERROR``. A set
  with both a 1 and another non-zero exit is non-zero either way (not pinned which);
- once ``<suite>.shards.jsonl`` exists, a block-1 shard set that is not exactly 0..N-1 of one N
  is ``MERGE_ERROR`` (a shard died before logging);
- no ``shards.jsonl`` at all (records produced in-process, as test_criteria_merge_xfail.py does):
  unchanged.

The shard records here are produced in-process (``calibrate`` under the shard variables) and the
exits are logged with ``_log_shard_done``, exactly what ``main`` writes. Unmarked: every commit.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import _criteria as C
import pytest

SUITE = "discordance"
STRICT = {"strict": True, "reason": "pending ADR"}
N = 3


@pytest.fixture
def log(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    for var in (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv(C.LOG_DIR_ENV, str(tmp_path))
    return tmp_path


def _records(suite: str, stream: str, value: float, log: Path, R: int = 6) -> None:
    """Block-1 and block-2 shard records (N shards) of a level check whose every p is ``value``."""
    with pytest.MonkeyPatch.context() as mp:
        for block in (1, 2):
            mp.setenv(C.BLOCK_ENV, str(block))
            for k in range(N):
                r = C.shard_replicates(R, N, k)
                mp.setenv(C.SHARD_ENV, f"{r.start}:{r.stop}@{k}/{N}")
                C.calibrate(suite, "level", 0.05, lambda _rng: value, R, stream=stream,
                            log_dir=log)  # fmt: skip


def _exits(suite: str, rcs: list[int], block: int = 1) -> None:
    for k, rc in enumerate(rcs):
        C._log_shard_done(suite, (k, len(rcs)), block, rc)


PASS, FAIL = 0.5, 0.0

#: (shard exits, merge exit); a set holding both 1 and another non-zero exit may be either
CASES: list[tuple[list[int], set[int]]] = [
    ([0, 0, 0], {0}),
    ([1, 1, 1], {C.MERGE_FAIL}),
    ([1, 0, 0], {C.MERGE_FAIL}),  # shard 0's direct assertion failed (B2)
    ([0, 1, 0], {C.MERGE_FAIL}),  # a shard-local power check failed
    ([4, 4, 4], {C.MERGE_ERROR}),  # usage error (PYTEST_ADDOPTS=--no-such-option)
    ([5, 5, 5], {C.MERGE_ERROR}),  # -k selecting nothing
    ([2, 2, 2], {C.MERGE_ERROR}),  # collection / import error
    ([0, 0, 3], {C.MERGE_ERROR}),  # internal error in one shard
    ([0, -9, 0], {C.MERGE_ERROR}),  # killed by a signal
    ([1, 4, 0], {C.MERGE_FAIL, C.MERGE_ERROR}),
]
IDS = ["-".join(map(str, rcs)) for rcs, _ in CASES]


# --------------------------------------------------------------------------- B1: no records


@pytest.mark.parametrize(("rcs", "want"), CASES, ids=IDS)
def test_suite_without_records_merges_its_shard_exits(
    log: Path, rcs: list[int], want: set[int]
) -> None:
    """The "no calibrate records; every shard finished" branch (denovo_precision-like suite)."""
    _exits("denovo_precision", rcs)
    assert C._merge_main("denovo_precision", True) in want


@pytest.mark.parametrize("suite", C.array_suites())
def test_every_suite_without_records_fails_closed(log: Path, suite: str) -> None:
    _exits(suite, [4, 4, 4])
    assert C._merge_main(suite, True) == C.MERGE_ERROR
    (log / f"{suite}.shards.jsonl").unlink()
    _exits(suite, [1, 0, 0])
    assert C._merge_main(suite, True) == C.MERGE_FAIL


# --------------------------------------------------------------------------- B2: with records


@pytest.mark.parametrize(("rcs", "want"), CASES, ids=IDS)
def test_passing_records_do_not_hide_a_shard_exit(
    log: Path, rcs: list[int], want: set[int]
) -> None:
    _records(SUITE, "discordance.rr.pass", PASS, log)
    _exits(SUITE, rcs)
    assert C._merge_main(SUITE, True) in want


@pytest.mark.parametrize("suite", C.array_suites())
def test_every_suite_folds_in_the_shard_exits(log: Path, suite: str) -> None:
    _records(suite, f"{suite}.rr.pass", PASS, log)
    _exits(suite, [1, 0, 0])
    assert C._merge_main(suite, True) == C.MERGE_FAIL
    (log / f"{suite}.shards.jsonl").unlink()
    _exits(suite, [0, 0, 4])
    assert C._merge_main(suite, True) == C.MERGE_ERROR
    (log / f"{suite}.shards.jsonl").unlink()
    _exits(suite, [0, 0, 0])
    assert C._merge_main(suite, True) == 0


def test_a_confirmed_fail_stays_a_fail_whatever_the_shard_exits(log: Path) -> None:
    _records(SUITE, "discordance.rr.fail", FAIL, log)
    _exits(SUITE, [0, 0, 0])
    assert C._merge_main(SUITE, True) == C.MERGE_FAIL


def test_a_block_two_shard_exit_counts_too(log: Path) -> None:
    """Every recorded shard exit is folded in, block 2 included."""
    _records(SUITE, "discordance.rr.pass", PASS, log)
    _exits(SUITE, [0, 0, 0])
    _exits(SUITE, [0, 1, 0], block=2)
    assert C._merge_main(SUITE, True) == C.MERGE_FAIL


@pytest.mark.parametrize("logged", [[0, 1], [1, 2], [0], [0, 1, 2, 3]])
def test_incomplete_shard_set_is_a_merge_error(log: Path, logged: list[int]) -> None:
    """``shards.jsonl`` exists but shard k of 3 never logged (died, timed out): no verdict, even
    though the records happen to cover every replicate."""
    _records(SUITE, "discordance.rr.pass", PASS, log)
    for k in logged:
        C._log_shard_done(SUITE, (k, N if k < N else N + 1), 1, 0)
    assert C._merge_main(SUITE, True) == C.MERGE_ERROR


def test_without_shards_log_the_records_decide(log: Path) -> None:
    """No ``shards.jsonl`` (in-process records): the merge is the records' verdict, as before."""
    _records(SUITE, "discordance.rr.pass", PASS, log)
    assert C._merge_main(SUITE, True) == 0


# --------------------------------------------------------------------------- mutants e-g


def test_xfail_check_written_by_two_tests_fails_closed(log: Path) -> None:
    """Mutant (e): ADR-0057 says an xfail check owned by two tests is a merge error."""
    _records(SUITE, "discordance.rr.two", FAIL, log)
    path = log / f"{SUITE}.jsonl"
    recs = C.read_results(path)
    for i, rec in enumerate(recs):
        rec["test"], rec["xfail"] = ("t.py::test_a" if i % 2 else "t.py::test_b"), STRICT
    path.write_text("".join(json.dumps(r) + "\n" for r in recs))
    with pytest.raises(ValueError, match="several tests"):
        C.merge_shards(SUITE, log)
    assert C._merge_main(SUITE, True) == C.MERGE_ERROR


def _rec(test: str, status: str) -> dict[str, Any]:
    return {"test": test, "xfail": dict(STRICT), "status": status}


@pytest.mark.parametrize("statuses", [["vacuous"], ["pass", "vacuous"], ["vacuous", "unconfirmed"]])
def test_a_vacuous_check_is_never_an_xpass(statuses: list[str]) -> None:
    """Mutant (f): rule 7, vacuous is neither a pass nor a fail, so a strict-xfail test with a
    vacuous check is not a strict XPASS and its outcome stays None."""
    merged = [_rec("t.py::test_v", s) for s in statuses]
    assert C.apply_xfail(merged, frozenset()) == []
    assert all(r["outcome"] is None for r in merged)


def _deferred(test: str, shards: list[tuple[int, int, int]]) -> list[dict[str, Any]]:
    return [{"suite": SUITE, "test": test, "outcome": C.DEFERRED, "xfail": dict(STRICT),
             "shard": list(s)} for s in shards]  # fmt: skip


@pytest.mark.parametrize(
    "shards",
    [[(0, 3, 2), (1, 3, 2), (2, 3, 2)],  # block 2 only
     [(0, 3, 1), (1, 3, 1), (2, 3, 2)],  # block-1 shard 2 missing, block 2 stands in
     [(0, 3, 2), (1, 3, 1), (2, 3, 1)]],
)  # fmt: skip
def test_direct_xpasses_ignore_block_two_shards(shards: list[tuple[int, int, int]]) -> None:
    """Mutant (g): only block-1 shards 0..N-1 make a clean run."""
    assert C.direct_xpasses(_deferred("t.py::test_d", shards), set()) == []


def test_direct_xpasses_need_every_block_one_shard() -> None:
    full = [(0, 3, 1), (1, 3, 1), (2, 3, 1)]
    assert C.direct_xpasses(_deferred("t.py::test_d", full), set()) == [("t.py::test_d", True)]
    with_b2 = [*full, (0, 3, 2)]
    assert C.direct_xpasses(_deferred("t.py::test_d", with_b2), set()) == [("t.py::test_d", True)]
