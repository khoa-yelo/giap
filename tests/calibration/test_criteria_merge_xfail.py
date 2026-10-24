"""``--merge`` honours the xfail mark of the test that owns a check (ADR-0057, accepted by Khoa
2026-10-24). Unsharded, pytest itself applies the mark; merged, the verdict must be the same:

- strict ``xfail(raises=AssertionError)`` and a confirmed ``fail``: an expected failure; reported
  ``xfail``, exit unchanged (0 when nothing else fails);
- strict xfail and every check of the test passes: strict XPASS, exit ``MERGE_FAIL`` (as pytest);
- non-strict xfail: reported (``xfail`` / ``xpass``), never changes the exit code;
- unmarked: a confirmed ``fail`` is ``MERGE_FAIL`` as before.

The owner comes from the records (``test`` = node ID, ``xfail`` = the evaluated mark, written by
``calibrate*`` while conftest's ``pytest_runtest_call`` names the running test), never a hand list.
The verdict is per *test*, like pytest's: a strict-xfail test with one failing and one passing
check is an expected failure, and a direct assertion that failed in shard 0 (an ``xfailed``
outcome record, ``log_test_outcome``) makes the test an expected failure even if its checks pass.
Here the records are produced in-process (block-1 and block-2 shards, so no confirm run is needed)
and their owner fields set as conftest would set them. Unmarked: runs on every commit.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import _criteria as C
import pytest

SUITE = "discordance"
STRICT = {"strict": True, "reason": "pending ADR"}
LOOSE = {"strict": False, "reason": "known weakness"}


@pytest.fixture
def log(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    for var in (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv(C.LOG_DIR_ENV, str(tmp_path))
    return tmp_path


def _shards(stream: str, value: float, test: str, xfail: dict[str, Any] | None, log: Path,
            R: int = 6, n: int = 3) -> None:  # fmt: skip
    """Block-1 and block-2 shard records of a level check whose every p is ``value``."""
    with pytest.MonkeyPatch.context() as mp:
        for block in (1, 2):
            mp.setenv(C.BLOCK_ENV, str(block))
            for k in range(n):
                r = C.shard_replicates(R, n, k)
                mp.setenv(C.SHARD_ENV, f"{r.start}:{r.stop}@{k}/{n}")
                C.calibrate(SUITE, "level", 0.05, lambda _rng: value, R, stream=stream,
                            log_dir=log)  # fmt: skip
    path = log / f"{SUITE}.jsonl"
    recs = C.read_results(path)
    for rec in recs:
        if rec.get("stream") == stream:
            rec["test"], rec["xfail"] = test, xfail
    path.write_text("".join(json.dumps(r) + "\n" for r in recs))


def _merged(log: Path) -> dict[str, dict[str, Any]]:
    return {r["stream"]: r for r in C.read_results(log / f"{SUITE}.merged.jsonl")}


FAIL, PASS = 0.0, 0.5


def test_strict_xfail_confirmed_fail_is_expected(log: Path) -> None:
    _shards("discordance.x.strict_fail", FAIL, "t.py::test_a", STRICT, log)
    assert C._merge_main(SUITE, True) == 0
    rec = _merged(log)["discordance.x.strict_fail"]
    assert rec["status"] == "fail" and rec["outcome"] == "xfail"


def test_strict_xfail_pass_is_a_strict_xpass(log: Path, capsys: pytest.CaptureFixture[str]) -> None:
    _shards("discordance.x.strict_pass", PASS, "t.py::test_b", STRICT, log)
    assert C._merge_main(SUITE, True) == C.MERGE_FAIL
    assert _merged(log)["discordance.x.strict_pass"]["outcome"] == "xpass"
    assert "XPASS(strict)" in capsys.readouterr().out


@pytest.mark.parametrize(("value", "outcome"), [(FAIL, "xfail"), (PASS, "xpass")])
def test_non_strict_xfail_is_reported_only(log: Path, value: float, outcome: str) -> None:
    _shards("discordance.x.loose", value, "t.py::test_c", LOOSE, log)
    assert C._merge_main(SUITE, True) == 0
    assert _merged(log)["discordance.x.loose"]["outcome"] == outcome


def test_unmarked_confirmed_fail_still_fails(log: Path) -> None:
    _shards("discordance.x.plain", FAIL, "t.py::test_d", None, log)
    _shards("discordance.x.strict_ok", FAIL, "t.py::test_e", STRICT, log)
    assert C._merge_main(SUITE, True) == C.MERGE_FAIL
    merged = _merged(log)
    assert merged["discordance.x.plain"].get("outcome") is None
    assert merged["discordance.x.strict_ok"]["outcome"] == "xfail"


def test_verdict_is_per_test_not_per_check(log: Path) -> None:
    """One strict-xfail test, two checks: one confirmed fail makes the whole test xfail (pytest
    raises at the first failing assertion), so the passing check is not an XPASS."""
    _shards("discordance.x.two_a", FAIL, "t.py::test_f", STRICT, log)
    _shards("discordance.x.two_b", PASS, "t.py::test_f", STRICT, log)
    assert C._merge_main(SUITE, True) == 0
    merged = _merged(log)
    assert merged["discordance.x.two_a"]["outcome"] == "xfail"
    assert merged["discordance.x.two_b"]["outcome"] == "xfail"


def test_direct_failure_in_shard_zero_makes_the_test_xfail(log: Path) -> None:
    """The check passes but the test's direct assertion failed in shard 0 (an ``xfailed``
    outcome record): unsharded pytest reports XFAIL, so the merge must not call it an XPASS."""
    _shards("discordance.x.direct", PASS, "t.py::test_g", STRICT, log)
    C.log_test_outcome(SUITE, "t.py::test_g", "xfailed", STRICT, log_dir=log)
    assert C._merge_main(SUITE, True) == 0
    assert _merged(log)["discordance.x.direct"]["outcome"] == "xfail"


def test_disagreeing_xfail_states_of_one_check_fail_closed(log: Path) -> None:
    _shards("discordance.x.mixed", FAIL, "t.py::test_h", STRICT, log)
    path = log / f"{SUITE}.jsonl"
    recs = C.read_results(path)
    recs[0]["xfail"] = None
    path.write_text("".join(json.dumps(r) + "\n" for r in recs))
    assert C._merge_main(SUITE, True) == C.MERGE_ERROR


def test_summary_names_the_expected_failure(log: Path, capsys: pytest.CaptureFixture[str]) -> None:
    _shards("discordance.x.summary", FAIL, "t.py::test_i", STRICT, log)
    assert C._merge_main(SUITE, True) == 0
    out = capsys.readouterr().out
    line = next(
        x for x in out.splitlines() if "discordance.x.summary" in x or "| discordance |" in x
    )
    assert "xfail" in line


def test_records_without_owner_fields_merge_as_before(log: Path) -> None:
    """Logs written before ADR-0057 (no ``test`` / ``xfail``) merge unchanged."""
    _shards("discordance.x.old", FAIL, "t.py::test_j", None, log)
    path = log / f"{SUITE}.jsonl"
    recs = [
        {k: v for k, v in r.items() if k not in ("test", "xfail")} for r in C.read_results(path)
    ]
    path.write_text("".join(json.dumps(r) + "\n" for r in recs))
    assert C._merge_main(SUITE, True) == C.MERGE_FAIL
