"""Shards of any R (``--shard K/N``, N up to 10) and the shard-0 owner of direct checks.

The discordance suite failed under ``--shard`` (Slurm 1885506-09): a check whose R differs from
the suite R (attack C's ``R_C = 300`` under ``replicates_full.discordance = 1000``) got the suite's
replicate range and raised; tests asserting ``status in {pass, unconfirmed}`` failed on a
``shard`` record; tests computing over every replicate repeated it in every shard. Pinned here:

- ``--shard K/N`` exports ``GIAP_CALIBRATION_REPLICATES=start:stop@K/N``; ``replicate_range(R)``
  is then ``shard_replicates(R, N, K)`` for *every* R (a check with R < N runs nothing in some
  shards; such an empty shard logs nothing). Without the ``@K/N`` suffix the variable means what
  it meant before (start:stop of the given R).
- Sharded and unsharded runs give the same records: the merged block-1 values equal the
  unsharded ones and the verdicts agree, for ``calibrate``, ``calibrate_pooled`` and
  ``calibrate_fdr_counts`` (R3: no bound, alpha or R changes; only where replicates run).
- ``owns_direct_checks()``: true unsharded and in shard 0 of block 1, false in every other shard,
  any block-2 shard and the ``--merge`` confirm run. Direct (non-``calibrate``) parts of a test run
  where it is true, so exactly one process evaluates them.
- ``skip_without_verdict(records)`` skips (``SHARD_SKIP``) when any record has no verdict.
- Records carry the ``test`` that wrote them (node ID) and its ``xfail`` state (ADR-0057).

Unmarked: runs on every commit.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import _criteria as C
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

FAST = settings(max_examples=40, deadline=None)
SUITE = "discordance"
SHARD_VARS = (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV)


@pytest.fixture(autouse=True)
def _clean(monkeypatch: pytest.MonkeyPatch) -> None:
    for var in SHARD_VARS:
        monkeypatch.delenv(var, raising=False)


def _spec(k: int, n: int, R_suite: int = 1000) -> str:
    r = C.shard_replicates(R_suite, n, k)
    return f"{r.start}:{r.stop}@{k}/{n}"


# --------------------------------------------------------------------------- replicate ranges


def test_shard_suffix_shards_a_check_r_other_than_the_suite_r(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The observed failure: shard 1/4 of R = 1000 is 250:500; a check with R = 300 runs 75:150."""
    monkeypatch.setenv(C.SHARD_ENV, _spec(1, 4))
    assert C.replicate_range(1000) == range(250, 500)
    assert C.replicate_range(300) == range(75, 150)
    assert C.shard_spec() == (1, 4)


def test_without_suffix_the_variable_keeps_its_meaning(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(C.SHARD_ENV, "2:5")
    assert C.replicate_range(10) == range(2, 5)
    assert C.shard_spec() is None
    with pytest.raises(ValueError):
        C.replicate_range(4)


@pytest.mark.parametrize("bad", ["0:3@3/3", "0:3@-1/3", "0:3@0/0", "0:3@a/b", "0:3@1", "x@0/2"])
def test_bad_suffix_fails_closed(monkeypatch: pytest.MonkeyPatch, bad: str) -> None:
    monkeypatch.setenv(C.SHARD_ENV, bad)
    with pytest.raises(ValueError):
        C.replicate_range(10)


@FAST
@given(R=st.integers(1, 60), n=st.integers(1, 10))
def test_suffix_shards_partition_every_r(R: int, n: int) -> None:
    with pytest.MonkeyPatch.context() as mp:
        got: list[int] = []
        for k in range(n):
            mp.setenv(C.SHARD_ENV, _spec(k, n))
            got += list(C.replicate_range(R))
    assert got == list(range(R))


# --------------------------------------------------------------------------- equivalence


def _uniform(rng: Any) -> float:
    return float(rng.uniform())


def _counts(i: int) -> dict[str, Any]:
    rng = C.replicate_rng("discordance.shard_spec.counts", i)
    n_pred = int(rng.integers(5, 15))
    n_false = int(rng.integers(0, 2))
    tp = n_pred - n_false
    return {"f1": {"tp_p": tp, "n_pred": n_pred, "tp_r": tp, "n_truth": tp + 1},
            "fdr": {"n_called": n_pred, "n_false": n_false}}  # fmt: skip


def _run_checks(R: int, log_dir: Path) -> list[dict[str, Any]]:
    return [
        C.calibrate(SUITE, "level", 0.05, _uniform, R, stream="discordance.shard_spec.level",
                    log_dir=log_dir),
        C.calibrate(SUITE, "fdr", 0.05, lambda rng: 0.04 * rng.uniform(), max(R, 2),
                    stream="discordance.shard_spec.fdr", log_dir=log_dir),
        C.calibrate_pooled(SUITE, "discordance.shard_spec.f1", "f1", 0.5, _counts, R,
                           log_dir=log_dir),
        C.calibrate_fdr_counts(SUITE, "discordance.shard_spec.fdrc", "fdr", 0.2, _counts,
                               max(R, 2), log_dir=log_dir),
    ]  # fmt: skip


def _essence(rec: dict[str, Any]) -> tuple[Any, ...]:
    return (rec["stream"], rec["criterion"], rec["level"], rec["R"], rec["block"], rec["status"],
            rec["observed"], rec["bound"], rec["values"])  # fmt: skip


@settings(max_examples=25, deadline=None)
@given(R=st.integers(1, 23), n=st.integers(1, 10), suite_r=st.sampled_from([100, 1000, 7]))
def test_sharded_merge_equals_unsharded(
    R: int, n: int, suite_r: int, tmp_path_factory: Any
) -> None:
    one, many = tmp_path_factory.mktemp("one"), tmp_path_factory.mktemp("many")
    with pytest.MonkeyPatch.context() as mp:
        for var in SHARD_VARS:
            mp.delenv(var, raising=False)
        _run_checks(R, one)
        for k in range(n):
            mp.setenv(C.SHARD_ENV, _spec(k, n, suite_r))
            for rec in _run_checks(R, many):  # a shard that is the check's whole R decides it
                whole = C.replicate_range(rec["R"]) == range(rec["R"])
                assert whole or rec["status"] == C.SHARD_STATUS, rec
        mp.delenv(C.SHARD_ENV)
        ref = {_essence(r)[:4]: _essence(r) for r in C.merge_shards(SUITE, one)}
        got = {_essence(r)[:4]: _essence(r) for r in C.merge_shards(SUITE, many)}
    assert got == ref


def test_empty_shard_logs_nothing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """R = 2 under 3 shards: shard 2 runs no replicate and appends no (malformed) record."""
    monkeypatch.setenv(C.SHARD_ENV, _spec(2, 3))
    assert C.replicate_range(2) == range(2, 2)
    rec = C.calibrate(SUITE, "level", 0.05, _uniform, 2, stream="discordance.shard_spec.empty",
                      log_dir=tmp_path)  # fmt: skip
    assert rec["status"] == C.SHARD_STATUS and rec["values"] == []
    pooled = C.calibrate_pooled(SUITE, "discordance.shard_spec.f1", "f1", 0.5, _counts, 2,
                                log_dir=tmp_path)  # fmt: skip
    assert pooled["status"] == C.SHARD_STATUS and pooled["values"] == []
    assert not (tmp_path / f"{SUITE}.jsonl").exists()


# --------------------------------------------------------------------------- owner and skip


@pytest.mark.parametrize(
    ("env", "owns"),
    [
        ({}, True),
        ({C.SHARD_ENV: "0:334@0/3"}, True),
        ({C.SHARD_ENV: "334:667@1/3"}, False),
        ({C.SHARD_ENV: "0:334@0/3", C.BLOCK_ENV: "2"}, False),
        ({C.SHARD_ENV: "0:5"}, True),  # legacy range: shard 0 is the one starting at 0
        ({C.SHARD_ENV: "5:10"}, False),
        ({C.CONFIRM_ENV: "[]"}, False),
    ],
)
def test_owns_direct_checks(
    monkeypatch: pytest.MonkeyPatch, env: dict[str, str], owns: bool
) -> None:
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    assert C.owns_direct_checks() is owns


@pytest.mark.parametrize("status", sorted(C.NO_VERDICT_STATUSES))
def test_skip_without_verdict_skips(status: str) -> None:
    with pytest.raises(pytest.skip.Exception, match="verdict comes from --merge"):
        C.skip_without_verdict([{"status": "pass"}, {"status": status}])
    with pytest.raises(pytest.skip.Exception):
        C.skip_without_verdict({"status": status})


@pytest.mark.parametrize("status", ["pass", "unconfirmed", "fail", C.VACUOUS_STATUS])
def test_skip_without_verdict_keeps_verdicts(status: str) -> None:
    C.skip_without_verdict([{"status": status}])
    C.skip_without_verdict({"status": status})


# --------------------------------------------------------------------------- CLI


def test_cli_shard_exports_the_suffix(capsys: pytest.CaptureFixture[str]) -> None:
    assert C.main(["--suite", SUITE, "--shard", "2/10", "--dry-run"]) == 0
    r = C.shard_replicates(C.replicates(SUITE, full=True), 10, 2)
    assert f"{C.SHARD_ENV}={r.start}:{r.stop}@2/10 " in capsys.readouterr().out


# --------------------------------------------------------------------------- owning test


def test_records_name_the_test_that_wrote_them(
    tmp_path: Path, request: pytest.FixtureRequest
) -> None:
    rec = C.calibrate(SUITE, "level", 0.05, _uniform, 3, stream="discordance.shard_spec.owner",
                      log_dir=tmp_path)  # fmt: skip
    logged = C.read_results(tmp_path / f"{SUITE}.jsonl")
    assert {r["test"] for r in logged} == {rec["test"]} == {request.node.nodeid}
    assert rec["xfail"] is None


@pytest.mark.xfail(strict=True, raises=AssertionError, reason="ADR-0057 owner record")
def test_records_carry_the_strict_xfail_of_their_test(tmp_path: Path) -> None:
    rec = C.calibrate(SUITE, "level", 0.05, _uniform, 3, stream="discordance.shard_spec.owner",
                      log_dir=tmp_path)  # fmt: skip
    if rec["xfail"] == {"strict": True, "reason": "ADR-0057 owner record"}:
        raise AssertionError("expected: the record names the strict xfail")


# --------------------------------------------------------------------------- suites without checks


def _fake_pytest(rc: int) -> Any:
    class Done:
        returncode = rc

    return lambda *_a, **_k: Done()


def test_shard_run_records_its_completion(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """``main --shard K/N`` appends ``{shard, n, block, pytest_exit}`` to ``<suite>.shards.jsonl``
    after its pytest run, so ``--merge`` can tell a suite without harness checks from dead shards."""
    monkeypatch.setenv(C.LOG_DIR_ENV, str(tmp_path))
    monkeypatch.setattr(C.subprocess, "run", _fake_pytest(1))
    assert C.main(["--suite", SUITE, "--shard", "1/3"]) == 1
    (rec,) = C.read_results(tmp_path / f"{SUITE}.shards.jsonl")
    assert rec == {"suite": SUITE, "shard": 1, "n": 3, "block": 1, "pytest_exit": 1}


@pytest.mark.parametrize(
    ("shards", "want"), [([0, 1, 2], 0), ([0, 2], C.MERGE_ERROR), ([], C.MERGE_ERROR)]
)
def test_merge_of_a_suite_without_harness_checks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, shards: list[int], want: int
) -> None:
    """denovo_precision has no ``calibrate`` check: its verdicts are the shard exit codes. Its merge
    passes when every block-1 shard finished, and stays a merge error when a shard is missing."""
    monkeypatch.setenv(C.LOG_DIR_ENV, str(tmp_path))
    monkeypatch.setattr(C.subprocess, "run", _fake_pytest(0))
    for k in shards:
        C.main(["--suite", SUITE, "--shard", f"{k}/3"])
    assert C._merge_main(SUITE, True) == want
