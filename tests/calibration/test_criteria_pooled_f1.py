"""Pooled F1 criterion of ``_criteria`` (docs/specs/M7-13-calibration-gate.md rule 2, rule 7,
§ Array "The harness gains a pooled criterion", Q4 / lead decision 2: point estimate, two-block
confirmation; ADR-0015). The harness extension is the implementer's; pinned here:

* ``C.pooled_f1(counts) -> (precision, recall, f1)`` over count dicts ``{n_pred, tp_p, n_truth,
  tp_r}``: P = sum tp_p / sum n_pred, Rc = sum tp_r / sum n_truth, F1 = 2 P Rc / (P + Rc); F1 = 0
  when sum n_pred = 0 or sum n_truth = 0. Fails closed (ValueError) on no replicates, negative or
  non-integer counts, tp_p > n_pred or tp_r > n_truth.
* ``C.jackknife_se(counts)``: leave-one-replicate-out pooled F1 values F_(i), SE =
  sqrt((R - 1) / R * sum (F_(i) - mean F_(.))^2); recorded only.
* ``C.check_pooled_min(counts, bound, *, guard=None) -> CheckResult``: ``observed`` = pooled F1,
  ``ok`` iff F1 >= bound and, for every ``guard`` item ``name: m``, sum of ``name`` >= m (a guard
  miss is a failed block, rule 7); ``R`` = number of replicates.
* ``C.calibrate_pooled(suite, stream, key, bound, records, R, *, guard=None, vacuous=None,
  log_dir=None) -> dict``: replicate i's counts are ``records(i)[key]`` (``records`` takes the
  absolute replicate index; it computes or caches one count record per replicate, rule 1).
  ADR-0015 blocks as ``calibrate``; ``vacuous`` (``{name: m}``) makes block 1 ``vacuous`` (no
  block 2, not a pass) when the block mean of ``name`` is below ``m`` (rule 6). Records carry
  ``criterion = "pooled_f1"``, ``key``, ``level = bound``, ``observed`` (F1), ``precision``,
  ``recall``, ``se``, ``block``, ``status``, ``start``, ``stop``, ``values`` (the count dicts).
  Sharded and ``GIPPA_CALIBRATION_CONFIRM`` runs as ``calibrate``; ``merge_shards`` gives pooled
  checks their verdict (pass | unconfirmed | fail | vacuous | pending), one per (stream, key).

Unmarked tests run on every commit; the operating-characteristic and jackknife checks are marked
``calibration`` (suite ``variants``).
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

SUITE, STREAM = "variants", "variants.selftest.pooled"


@pytest.fixture(autouse=True)
def _unsharded(monkeypatch: pytest.MonkeyPatch) -> None:
    for var in (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV):
        monkeypatch.delenv(var, raising=False)


def cnt(tp_p: int, n_pred: int, tp_r: int, n_truth: int) -> dict[str, int]:
    return {"tp_p": tp_p, "n_pred": n_pred, "tp_r": tp_r, "n_truth": n_truth}


def table(good: set[int] | None = None, bad: tuple[int, int, int, int] = (0, 10, 0, 10)) -> Any:
    """``records(i)``: perfect counts for i in ``good`` (all when None), else ``bad``."""
    calls: list[int] = []

    def records(i: int) -> dict[str, Any]:
        calls.append(i)
        ok = good is None or i in good
        return {"k": cnt(10, 10, 10, 10) if ok else cnt(*bad), "other": cnt(10, 10, 10, 10)}

    records.calls = calls  # type: ignore[attr-defined]
    return records


def test_example1_single_replicate() -> None:
    p, r, f = C.pooled_f1([cnt(8, 10, 7, 9)])
    assert (p, r, f) == pytest.approx((0.8, 7 / 9, 0.789), abs=5e-4)


def test_example2_pooling_counts_truth_of_silent_replicate() -> None:
    p, r, f = C.pooled_f1([cnt(8, 10, 7, 9), cnt(0, 0, 0, 3)])
    assert (p, r, f) == pytest.approx((0.8, 7 / 12, 0.675), abs=5e-4)


def test_degenerate_zero_sums_give_zero_f1() -> None:
    assert C.pooled_f1([cnt(0, 0, 0, 5)])[2] == 0.0
    assert C.pooled_f1([cnt(0, 4, 0, 0)])[2] == 0.0
    assert C.pooled_f1([cnt(0, 0, 0, 0), cnt(0, 0, 0, 0)])[2] == 0.0


@pytest.mark.parametrize("bad", [[], [cnt(3, 2, 0, 1)], [cnt(0, 1, 2, 1)], [cnt(-1, 1, 0, 1)],
                                 [{"tp_p": 0.5, "n_pred": 1, "tp_r": 0, "n_truth": 1}],
                                 [{"tp_p": 0, "n_pred": 1}]])  # fmt: skip
def test_fail_closed(bad: list[dict[str, Any]]) -> None:
    with pytest.raises((ValueError, KeyError)):
        C.pooled_f1(bad)


def test_jackknife_small_example() -> None:
    reps = [cnt(8, 10, 7, 9), cnt(9, 10, 9, 10), cnt(5, 10, 6, 10)]
    loo = [C.pooled_f1([r for j, r in enumerate(reps) if j != i])[2] for i in range(3)]
    want = math.sqrt(2 / 3 * sum((x - np.mean(loo)) ** 2 for x in loo))
    assert C.jackknife_se(reps) == pytest.approx(want)
    assert C.jackknife_se([cnt(5, 5, 5, 5)] * 4) == pytest.approx(0.0)


COUNTS = st.lists(st.tuples(st.integers(0, 30), st.integers(0, 30), st.integers(0, 30),
                            st.integers(0, 30)).map(lambda t: cnt(min(t[0], t[1]), t[1],
                                                                  min(t[2], t[3]), t[3])),
                  min_size=1, max_size=12)  # fmt: skip


@settings(max_examples=100, deadline=None)
@given(reps=COUNTS, perm=st.randoms(use_true_random=False))
def test_pooled_f1_is_a_function_of_the_sums(reps: list[dict[str, int]], perm: Any) -> None:
    total = {k: sum(r[k] for r in reps) for k in reps[0]}
    shuffled = list(reps)
    perm.shuffle(shuffled)
    f = C.pooled_f1(reps)[2]
    assert f == pytest.approx(C.pooled_f1([total])[2]) == pytest.approx(C.pooled_f1(shuffled)[2])
    assert 0.0 <= f <= 1.0
    p, r = total["tp_p"], total["tp_r"]
    assert (f == 0.0) == (p == 0 or r == 0)


def test_check_pooled_min_bound_and_guard() -> None:
    reps = [cnt(9, 10, 9, 10)] * 6  # F1 0.9, 60 predictions, 60 recall-truth rows
    res = C.check_pooled_min(reps, 0.9)
    assert res.ok and res.observed == pytest.approx(0.9) and res.bound == 0.9 and res.R == 6
    assert not C.check_pooled_min(reps, 0.91).ok
    assert C.check_pooled_min(reps, 0.9, guard={"n_pred": 60, "n_truth": 60}).ok
    miss = C.check_pooled_min(reps, 0.9, guard={"n_pred": 50, "n_truth": 61})
    assert not miss.ok and "guard" in miss.detail


def _run(records: Any, tmp: Path, R: int = 4, **kw: Any) -> dict[str, Any]:
    return C.calibrate_pooled(SUITE, STREAM, "k", 0.9, records, R, log_dir=tmp, **kw)


def test_pass_runs_block_one_only(tmp_path: Path) -> None:
    rec_fn = table()
    rec = _run(rec_fn, tmp_path)
    assert rec["status"] == "pass" and rec["block"] == 1 and rec_fn.calls == [0, 1, 2, 3]
    assert rec["criterion"] == "pooled_f1" and rec["key"] == "k" and rec["level"] == 0.9
    assert rec["observed"] == pytest.approx(1.0) and rec["values"] == [cnt(10, 10, 10, 10)] * 4
    assert {"precision", "recall", "se"} <= set(rec)
    logged = C.read_results(tmp_path / f"{SUITE}.jsonl")
    assert [r["status"] for r in logged] == ["pass"]


def test_unconfirmed_when_only_block_one_fails(tmp_path: Path) -> None:
    rec_fn = table(good=set(range(4, 8)))
    rec = _run(rec_fn, tmp_path)
    assert rec["status"] == "unconfirmed" and rec["block"] == 2
    assert rec_fn.calls == list(range(8))  # block 2 is replicates R..2R-1


def test_fail_when_both_blocks_fail(tmp_path: Path) -> None:
    rec = _run(table(good=set()), tmp_path)
    assert rec["status"] == "fail" and rec["observed"] == 0.0
    blocks = [r["block"] for r in C.read_results(tmp_path / f"{SUITE}.jsonl")]
    assert blocks == [1, 2]


def test_guard_miss_fails_like_a_failed_block(tmp_path: Path) -> None:
    rec = _run(table(), tmp_path, guard={"n_truth": 50})  # 4 x 10 recall-truth rows < 50
    assert rec["status"] == "fail"


def test_vacuous_block_is_neither_pass_nor_fail(tmp_path: Path) -> None:
    rec_fn = table()
    rec = _run(rec_fn, tmp_path, vacuous={"n_pred": 11.0})  # block mean 10 < 11
    assert rec["status"] == "vacuous" and rec_fn.calls == [0, 1, 2, 3]
    assert _run(table(), tmp_path, vacuous={"n_pred": 10.0})["status"] == "pass"


def test_unknown_suite_refused(tmp_path: Path) -> None:
    with pytest.raises(KeyError):
        C.calibrate_pooled("no_such_suite", STREAM, "k", 0.9, table(), 4, log_dir=tmp_path)


def _shards(monkeypatch: pytest.MonkeyPatch, tmp: Path, rec_fn: Any, parts: list[str],
            block: str = "1", R: int = 6, **kw: Any) -> None:  # fmt: skip
    for p in parts:
        monkeypatch.setenv(C.SHARD_ENV, p)
        monkeypatch.setenv(C.BLOCK_ENV, block)
        rec = C.calibrate_pooled(SUITE, STREAM, "k", 0.9, rec_fn, R, log_dir=tmp, **kw)
        assert rec["status"] == C.SHARD_STATUS
    monkeypatch.delenv(C.SHARD_ENV)
    monkeypatch.delenv(C.BLOCK_ENV)


def test_shards_merge_to_the_unsharded_verdict(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    good = {0, 1, 2, 4, 6, 7, 8, 9, 10, 11}  # block 1 F1 < 0.9, block 2 passes
    whole = C.calibrate_pooled(SUITE, STREAM, "k", 0.9, table(good), 6, log_dir=tmp_path / "a")
    _shards(monkeypatch, tmp_path / "b", table(good), ["0:2", "2:5", "5:6"])
    pend = [r for r in C.merge_shards(SUITE, tmp_path / "b") if r.get("key") == "k"]
    assert [r["status"] for r in pend] == [C.PENDING_STATUS]
    _shards(monkeypatch, tmp_path / "b", table(good), ["0:3", "3:6"], block="2")
    merged = [r for r in C.merge_shards(SUITE, tmp_path / "b") if r.get("key") == "k"]
    assert len(merged) == 1 and merged[0]["status"] == whole["status"] == "unconfirmed"
    assert merged[0]["observed"] == pytest.approx(whole["observed"])
    assert merged[0]["values"] == whole["values"]


def test_merge_refuses_missing_replicates(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _shards(monkeypatch, tmp_path, table(), ["0:2", "3:6"])
    with pytest.raises(ValueError):
        C.merge_shards(SUITE, tmp_path)


def test_confirm_runs_block_two_of_the_listed_pooled_check_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _shards(monkeypatch, tmp_path, table(good={6, 7, 8, 9, 10, 11}), ["0:3", "3:6"])
    other = C.calibrate_pooled(SUITE, STREAM, "other", 0.9, table(), 6, log_dir=tmp_path)
    assert other["status"] == "pass"
    (pending,) = [r for r in C.merge_shards(SUITE, tmp_path) if r["status"] == C.PENDING_STATUS]
    monkeypatch.setenv(C.CONFIRM_ENV, json.dumps([C._check_key(pending)]))
    rec_fn = table(good={6, 7, 8, 9, 10, 11})
    assert C.calibrate_pooled(SUITE, STREAM, "k", 0.9, rec_fn, 6, log_dir=tmp_path)["block"] == 2
    assert rec_fn.calls == list(range(6, 12))
    skip_fn = table()
    skipped = C.calibrate_pooled(SUITE, STREAM, "other", 0.9, skip_fn, 6, log_dir=tmp_path)
    assert skipped["status"] == C.SKIPPED_STATUS and skip_fn.calls == []
    monkeypatch.delenv(C.CONFIRM_ENV)
    final = {r["key"]: r["status"] for r in C.merge_shards(SUITE, tmp_path) if "key" in r}
    assert final == {"k": "unconfirmed", "other": "pass"}


# --------------------------------------------------------------------------- merge path
# Lead decision 7 (code review B3, N2-N4): the merge must reach the unsharded verdict on guards,
# vacuity and settings, fail closed on stale or non-finite records, and the job's exit code must
# carry the confirm run's pytest exit code.

GUARD = {"n_truth": 100}  # 6 x 10 recall-truth rows < 100: a guard miss in each block


def _keyed(tmp: Path) -> list[dict[str, Any]]:
    return [r for r in C.merge_shards(SUITE, tmp) if r.get("key") == "k"]


def test_merge_guard_miss_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    whole = C.calibrate_pooled(SUITE, STREAM, "k", 0.9, table(), 6, log_dir=tmp_path / "a",
                               guard=GUARD)  # fmt: skip
    assert whole["status"] == "fail"
    b = tmp_path / "b"
    _shards(monkeypatch, b, table(), ["0:3", "3:6"], guard=GUARD)
    assert [r["status"] for r in _keyed(b)] == [C.PENDING_STATUS]  # F1 1.0, yet block 1 fails
    _shards(monkeypatch, b, table(), ["0:2", "2:6"], block="2", guard=GUARD)
    assert [r["status"] for r in _keyed(b)] == ["fail"]


def test_merge_guard_miss_in_block_one_only_is_unconfirmed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rec_fn = table(good=set(range(12)))  # F1 1.0 everywhere
    big = {"n_truth": 60}
    _shards(monkeypatch, tmp_path, rec_fn, ["0:3", "3:6"], guard=big)
    assert [r["status"] for r in _keyed(tmp_path)] == ["pass"]  # 60 >= 60
    c = tmp_path / "c"
    _shards(monkeypatch, c, rec_fn, ["0:3", "3:6"], guard={"n_truth": 61})
    assert [r["status"] for r in _keyed(c)] == [C.PENDING_STATUS]


def test_merge_vacuous_block_one_is_vacuous(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    vac = {"n_pred": 11.0}  # block mean 10 < 11
    whole = C.calibrate_pooled(SUITE, STREAM, "k", 0.9, table(), 6, log_dir=tmp_path / "a",
                               vacuous=vac)  # fmt: skip
    _shards(monkeypatch, tmp_path / "b", table(), ["0:3", "3:6"], vacuous=vac)
    merged = _keyed(tmp_path / "b")
    assert [r["status"] for r in merged] == [whole["status"]] == [C.VACUOUS_STATUS]
    assert merged[0]["block"] == 1


@pytest.mark.parametrize(
    ("first", "second"),
    [({"guard": {"n_truth": 10}}, {"guard": {"n_truth": 20}}),
     ({"guard": {"n_truth": 10}}, {}),
     ({"vacuous": {"n_pred": 1.0}}, {"vacuous": {"n_pred": 2.0}}),
     ({}, {"vacuous": {"n_pred": 1.0}})],
)  # fmt: skip
def test_merge_refuses_shards_that_disagree_on_settings(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, first: dict[str, Any], second: dict[str, Any]
) -> None:
    _shards(monkeypatch, tmp_path, table(), ["0:3"], **first)
    _shards(monkeypatch, tmp_path, table(), ["3:6"], **second)
    with pytest.raises(ValueError, match="disagree"):
        C.merge_shards(SUITE, tmp_path)
    monkeypatch.setenv(C.LOG_DIR_ENV, str(tmp_path))
    assert C._merge_main(SUITE, False) == C.MERGE_ERROR


def test_merge_refuses_a_replicate_logged_twice_with_different_values(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``_collect``'s stale-merge guard (code review N2), pooled records."""
    _shards(monkeypatch, tmp_path, table(), ["0:3", "3:6"])
    _shards(monkeypatch, tmp_path, table(good=set()), ["0:3"])  # a stale rerun, other counts
    with pytest.raises(ValueError, match="twice"):
        C.merge_shards(SUITE, tmp_path)


def test_merge_accepts_a_replicate_logged_twice_with_the_same_values(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _shards(monkeypatch, tmp_path, table(), ["0:3", "3:6", "0:3"])
    assert [r["status"] for r in _keyed(tmp_path)] == ["pass"]


def test_merge_refuses_a_stale_calibrate_record(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The same guard for scalar ``calibrate`` records (predates the PR; N2)."""
    for part, value in (("0:5", 0.0), ("5:10", 0.0), ("0:5", 1.0)):
        monkeypatch.setenv(C.SHARD_ENV, part)
        monkeypatch.setenv(C.BLOCK_ENV, "1")
        C.calibrate(SUITE, "fwer", 0.05, lambda rng, v=value: v, 10, stream=f"{STREAM}.fwer",
                    log_dir=tmp_path)  # fmt: skip
    with pytest.raises(ValueError, match="twice"):
        C.merge_shards(SUITE, tmp_path)


@pytest.mark.parametrize("bad", [float("inf"), float("-inf"), float("nan")])
def test_merge_non_finite_count_is_a_merge_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, bad: float
) -> None:
    """A non-finite count in a shard record (json round-trips ``Infinity``) fails closed as a
    ValueError, so ``--merge`` exits MERGE_ERROR, never an uncaught OverflowError (N3)."""
    _shards(monkeypatch, tmp_path, table(), ["0:3", "3:6"])
    path = tmp_path / f"{SUITE}.jsonl"
    recs = C.read_results(path)
    recs[0]["values"][1]["n_pred"] = bad
    path.write_text("".join(json.dumps(r) + "\n" for r in recs))
    with pytest.raises(ValueError):
        C.merge_shards(SUITE, tmp_path)
    monkeypatch.setenv(C.LOG_DIR_ENV, str(tmp_path))
    assert C._merge_main(SUITE, False) == C.MERGE_ERROR


def _pending_then_block2(
    tmp: Path, monkeypatch: pytest.MonkeyPatch, code: int, good2: bool
) -> list[Any]:
    """Block 1 fails in shards; the fake confirm run logs block 2 and exits ``code``."""
    _shards(monkeypatch, tmp, table(good=set(range(6, 12))), ["0:3", "3:6"])
    monkeypatch.setenv(C.LOG_DIR_ENV, str(tmp))
    seen: list[Any] = []

    def fake_block2(suite: str, full: bool, pending: Any, lanes: int | None = None) -> int:
        seen.append([C._check_key(r) for r in pending])
        rec_fn = table(good=set(range(6, 12)) if good2 else set())
        _shards(monkeypatch, tmp, rec_fn, ["0:3", "3:6"], block="2")
        monkeypatch.setenv(C.LOG_DIR_ENV, str(tmp))
        return code

    monkeypatch.setattr(C, "_run_block2", fake_block2)
    return seen


@pytest.mark.parametrize("code", [1, 2, 4])
def test_merge_fails_when_the_confirm_run_pytest_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, code: int
) -> None:
    """Lead decision 7 (statistician B2, code review B1): block 2 passes (``unconfirmed``) but the
    confirm run's pytest exited non-zero (a test outside the harness failed): not exit 0."""
    seen = _pending_then_block2(tmp_path, monkeypatch, code, good2=True)
    assert C._merge_main(SUITE, False) != 0
    assert seen == [[[STREAM, C.POOLED_F1, 0.9, 6, "k"]]]


@pytest.mark.parametrize(("good2", "want"), [(True, 0), (False, C.MERGE_FAIL)])
def test_merge_exit_with_a_clean_confirm_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, good2: bool, want: int
) -> None:
    _pending_then_block2(tmp_path, monkeypatch, 0, good2=good2)
    assert C._merge_main(SUITE, False) == want


def _capture_child_env(monkeypatch: pytest.MonkeyPatch) -> list[dict[str, str]]:
    envs: list[dict[str, str]] = []

    class Done:
        returncode = 0

    def fake_run(cmd: Any, env: dict[str, str] | None = None, **kw: Any) -> Done:
        envs.append(dict(env or {}))
        return Done()

    monkeypatch.setattr(C.subprocess, "run", fake_run)
    return envs


@pytest.mark.parametrize("argv", [["--shard", "0/2"], ["--shard=1/2"], []])
def test_main_strips_a_leaked_confirm_variable(
    monkeypatch: pytest.MonkeyPatch, argv: list[str]
) -> None:
    """A ``GIPPA_CALIBRATION_CONFIRM`` left in the environment would make every gate test skip
    and exit 0; ``main()`` removes it from the shard (and unsharded) pytest environment."""
    monkeypatch.setenv(C.CONFIRM_ENV, json.dumps([[STREAM, "fwer", 0.05, 6]]))
    envs = _capture_child_env(monkeypatch)
    assert C.main(["--suite", SUITE, *argv]) == 0
    assert len(envs) == 1 and C.CONFIRM_ENV not in envs[0], envs


@pytest.mark.parametrize(("argv", "want"), [(["--shard=0/2"], range(0, 50)),
                                            (["--shard=1/2"], range(50, 100)),
                                            (["--shard", "1/2"], range(50, 100))])  # fmt: skip
def test_parse_args_shard_forms(argv: list[str], want: range) -> None:
    ns = C.parse_args(["--suite", SUITE, *argv])
    assert ns.replicates == want and not ns.merge
    with pytest.raises(SystemExit):
        C.parse_args(["--suite", SUITE, "--shard=2/2"])  # 0-based: K < N


def test_parse_args_merge_form() -> None:
    assert C.parse_args(["--suite", SUITE, "--merge"]).merge


# --------------------------------------------------------------------------- calibration


def _binomial_records(f_p: float, f_r: float, n: int = 40) -> Any:
    def records(i: int) -> dict[str, Any]:
        g = C.replicate_rng(f"{STREAM}.oc{f_p}", i)
        return {"k": cnt(int(g.binomial(n, f_p)), n, int(g.binomial(n, f_r)), n)}

    return records


@pytest.mark.calibration
@pytest.mark.calsuite("variants")
@pytest.mark.parametrize(("f", "want"), [(0.95, "pass"), (0.84, "fail")])
def test_operating_characteristic(f: float, want: str, teeth_log: Path) -> None:
    """At R = replicates_fast.variants and 40 predictions / truth rows per replicate, a method at
    true F1 0.95 passes and one at 0.84 (risk 2's preview) is a confirmed fail."""
    R = C.replicates(SUITE, False)
    rec = C.calibrate_pooled(SUITE, f"{STREAM}.oc{f}", "k", 0.9, _binomial_records(f, f), R,
                             log_dir=teeth_log)  # fmt: skip
    assert rec["status"] == want, rec


@pytest.mark.calibration
@pytest.mark.calsuite("variants")
def test_jackknife_se_tracks_the_spread_of_pooled_f1() -> None:
    """Recorded SE (rule 2): over 200 blocks of R = 20 binomial replicates, the mean jackknife SE
    is within 25 % of the empirical SD of the pooled F1."""
    R, B = C.replicates(SUITE, False), 200
    rec = _binomial_records(0.85, 0.8, n=30)
    f1, se = [], []
    for b in range(B):
        block = [rec(b * R + i)["k"] for i in range(R)]
        f1.append(C.pooled_f1(block)[2])
        se.append(C.jackknife_se(block))
    ratio = float(np.mean(se)) / float(np.std(f1, ddof=1))
    assert 0.8 <= ratio <= 1.25, ratio
