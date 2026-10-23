"""Count-record FDR criterion of ``_criteria`` for the M7-13b novel-ORF gate (docs/specs/
M7-13-calibration-gate.md rules 1, 6, G4 / G4-k3 rows, § Array "shard-mergeable, ADR-0015
blocks, statuses pass / unconfirmed / fail / vacuous"; §8.4 FDR criterion "Mean FDP over R
replicates <= q + 2 SE(FDP)"). The harness extension is the implementer's; pinned here:

* ``C.FDR_COUNTS == "fdr_counts"`` (the record ``criterion``).
* ``C.fdp_of(counts) -> float``: ``n_false / max(n_called, 1)`` (rule 6: 0 when nothing is
  called); fails closed (ValueError / KeyError) on a missing key, a negative or non-integer
  count, or ``n_false > n_called``.
* ``C.check_fdr_counts(values, q) -> CheckResult``: ``check_fdr`` on the per-replicate FDPs.
* ``C.calibrate_fdr_counts(suite, stream, key, q, records, R, *, vacuous=None, log_dir=None)
  -> dict``: replicate i's counts are ``records(i)[key]`` (absolute index, block 2 is
  R..2R-1), mirroring ``calibrate_pooled``: ADR-0015 blocks; ``vacuous`` (``{name: m}``) makes
  block 1 ``vacuous`` (no block 2, not a pass) when the block mean of ``name`` is below ``m``
  (rule 6: mean list targets below 1/q). Records carry ``criterion = "fdr_counts"``, ``key``,
  ``level = q``, ``observed`` (mean FDP), ``bound``, ``block``, ``status``, ``start``, ``stop``,
  ``values`` (the count dicts) and ``vacuous``. Sharded and ``GIAP_CALIBRATION_CONFIRM`` runs
  as ``calibrate_pooled`` (confirm key ``C._check_key(record)`` includes ``key``);
  ``merge_shards`` gives one verdict per (stream, key, level, R) (pass | unconfirmed | fail |
  vacuous | pending).

Unmarked tests run on every commit; the operating characteristic is marked ``calibration``
(suite ``novel_orfs``).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

SUITE, STREAM = "novel_orfs", "novel_orfs.selftest.fdr_counts"


@pytest.fixture(autouse=True)
def _unsharded(monkeypatch: pytest.MonkeyPatch) -> None:
    for var in (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV):
        monkeypatch.delenv(var, raising=False)


def cnt(n_called: int, n_false: int, n_targets: int = 500) -> dict[str, int]:
    return {"n_called": n_called, "n_false": n_false, "n_targets": n_targets}


def table(good: set[int] | None = None, bad: tuple[int, int] = (100, 50),
          n_targets: int = 500) -> Any:  # fmt: skip
    """``records(i)``: FDP 0 for i in ``good`` (all when None), else ``bad`` = (called, false)."""
    calls: list[int] = []

    def records(i: int) -> dict[str, Any]:
        calls.append(i)
        ok = good is None or i in good
        k = cnt(100, 0, n_targets) if ok else cnt(*bad, n_targets)
        return {"k": k, "other": cnt(100, 0)}

    records.calls = calls  # type: ignore[attr-defined]
    return records


def test_constant() -> None:
    assert C.FDR_COUNTS == "fdr_counts"


def test_example9_fdp() -> None:
    """Spec example 9: cds_outframe at q 0.05, 120 called, 4 not in truth -> FDP 0.0333."""
    assert C.fdp_of(cnt(120, 4)) == pytest.approx(0.0333, abs=5e-5)


def test_no_call_gives_fdp_zero() -> None:
    assert C.fdp_of(cnt(0, 0)) == 0.0


@pytest.mark.parametrize("bad", [{"n_called": 2, "n_false": 3}, {"n_called": -1, "n_false": 0},
                                 {"n_called": 2.5, "n_false": 0}, {"n_called": True,
                                                                   "n_false": 0},
                                 {"n_called": 2}])  # fmt: skip
def test_fdp_fails_closed(bad: dict[str, Any]) -> None:
    with pytest.raises((ValueError, KeyError)):
        C.fdp_of(bad)


@settings(max_examples=100, deadline=None)
@given(st.lists(st.tuples(st.integers(0, 50), st.integers(0, 50)), min_size=2, max_size=15),
       st.sampled_from([0.01, 0.05, 0.10]))  # fmt: skip
def test_check_fdr_counts_is_check_fdr_on_the_fdps(pairs: list[tuple[int, int]], q: float) -> None:
    values = [cnt(max(a, b), min(a, b)) for a, b in pairs]
    fdps = [min(a, b) / max(max(a, b), 1) for a, b in pairs]
    got, want = C.check_fdr_counts(values, q), C.check_fdr(fdps, q)
    assert (got.ok, got.R) == (want.ok, want.R)
    assert got.observed == pytest.approx(want.observed)
    assert got.bound == pytest.approx(want.bound)


def _run(records: Any, tmp: Path, R: int = 4, q: float = 0.05, **kw: Any) -> dict[str, Any]:
    return C.calibrate_fdr_counts(SUITE, STREAM, "k", q, records, R, log_dir=tmp, **kw)


def test_pass_runs_block_one_only(tmp_path: Path) -> None:
    rec_fn = table()
    rec = _run(rec_fn, tmp_path)
    assert rec["status"] == "pass" and rec["block"] == 1 and rec_fn.calls == [0, 1, 2, 3]
    assert rec["criterion"] == "fdr_counts" and rec["key"] == "k" and rec["level"] == 0.05
    assert rec["observed"] == 0.0 and rec["values"] == [cnt(100, 0)] * 4
    assert [r["status"] for r in C.read_results(tmp_path / f"{SUITE}.jsonl")] == ["pass"]


def test_unconfirmed_when_only_block_one_fails(tmp_path: Path) -> None:
    rec_fn = table(good=set(range(4, 8)))
    rec = _run(rec_fn, tmp_path)
    assert rec["status"] == "unconfirmed" and rec["block"] == 2
    assert rec_fn.calls == list(range(8))


def test_fail_when_both_blocks_fail(tmp_path: Path) -> None:
    rec = _run(table(good=set()), tmp_path)
    assert rec["status"] == "fail" and rec["observed"] == pytest.approx(0.5)
    assert [r["block"] for r in C.read_results(tmp_path / f"{SUITE}.jsonl")] == [1, 2]


def test_vacuous_block_is_neither_pass_nor_fail(tmp_path: Path) -> None:
    """Rule 6: mean list targets 50 < 1/q = 100 -> vacuous, block 2 not run (even when block 1
    would fail); at exactly 100 the check is judged."""
    rec_fn = table(good=set(), n_targets=50)
    rec = _run(rec_fn, tmp_path, q=0.01, vacuous={"n_targets": 100.0})
    assert rec["status"] == C.VACUOUS_STATUS and rec_fn.calls == [0, 1, 2, 3]
    judged = _run(table(n_targets=100), tmp_path, q=0.01, vacuous={"n_targets": 100.0})
    assert judged["status"] == "pass"


def test_fail_closed_before_logging(tmp_path: Path) -> None:
    def broken(i: int) -> dict[str, Any]:
        return {"k": {"n_called": 1, "n_false": 2, "n_targets": 1}}

    with pytest.raises(ValueError):
        _run(broken, tmp_path)
    assert not (tmp_path / f"{SUITE}.jsonl").exists()


def test_unknown_suite_refused(tmp_path: Path) -> None:
    with pytest.raises(KeyError):
        C.calibrate_fdr_counts("no_such_suite", STREAM, "k", 0.05, table(), 4, log_dir=tmp_path)


def _shards(monkeypatch: pytest.MonkeyPatch, tmp: Path, rec_fn: Any, parts: list[str],
            block: str = "1", R: int = 6, **kw: Any) -> None:  # fmt: skip
    for p in parts:
        monkeypatch.setenv(C.SHARD_ENV, p)
        monkeypatch.setenv(C.BLOCK_ENV, block)
        rec = C.calibrate_fdr_counts(SUITE, STREAM, "k", 0.05, rec_fn, R, log_dir=tmp, **kw)
        assert rec["status"] == C.SHARD_STATUS
    monkeypatch.delenv(C.SHARD_ENV)
    monkeypatch.delenv(C.BLOCK_ENV)


def test_shards_merge_to_the_unsharded_verdict(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    good = {0, 1, 2, 4, 6, 7, 8, 9, 10, 11}  # block 1 fails (two FDP 0.5), block 2 passes
    whole = C.calibrate_fdr_counts(SUITE, STREAM, "k", 0.05, table(good), 6,
                                   log_dir=tmp_path / "a")  # fmt: skip
    _shards(monkeypatch, tmp_path / "b", table(good), ["0:2", "2:5", "5:6"])
    pend = [r for r in C.merge_shards(SUITE, tmp_path / "b") if r.get("key") == "k"]
    assert [r["status"] for r in pend] == [C.PENDING_STATUS]
    _shards(monkeypatch, tmp_path / "b", table(good), ["0:3", "3:6"], block="2")
    merged = [r for r in C.merge_shards(SUITE, tmp_path / "b") if r.get("key") == "k"]
    assert len(merged) == 1 and merged[0]["status"] == whole["status"] == "unconfirmed"
    assert merged[0]["criterion"] == C.FDR_COUNTS
    assert merged[0]["observed"] == pytest.approx(whole["observed"])
    assert merged[0]["values"] == whole["values"]


def test_merge_keeps_vacuity(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A vacuous check stays vacuous after the merge (rule 6: never counted as a pass)."""
    _shards(monkeypatch, tmp_path, table(n_targets=20), ["0:3", "3:6"],
            vacuous={"n_targets": 100.0})  # fmt: skip
    (merged,) = [r for r in C.merge_shards(SUITE, tmp_path) if r.get("key") == "k"]
    assert merged["status"] == C.VACUOUS_STATUS


def test_merge_refuses_missing_replicates(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _shards(monkeypatch, tmp_path, table(), ["0:2", "3:6"])
    with pytest.raises(ValueError):
        C.merge_shards(SUITE, tmp_path)


def test_levels_are_separate_checks(tmp_path: Path) -> None:
    """One stream and key at q 0.05 and 0.10 are two checks (G4 runs three levels per class)."""
    _run(table(), tmp_path, q=0.05)
    _run(table(), tmp_path, q=0.10)
    merged = [r for r in C.merge_shards(SUITE, tmp_path) if r.get("key") == "k"]
    assert sorted(r["level"] for r in merged) == [0.05, 0.10]


def test_confirm_runs_block_two_of_the_listed_check_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _shards(monkeypatch, tmp_path, table(good={6, 7, 8, 9, 10, 11}), ["0:3", "3:6"])
    other = C.calibrate_fdr_counts(SUITE, STREAM, "other", 0.05, table(), 6, log_dir=tmp_path)
    assert other["status"] == "pass"
    (pending,) = [r for r in C.merge_shards(SUITE, tmp_path) if r["status"] == C.PENDING_STATUS]
    key = C._check_key(pending)
    assert key[-1] == "k" and key[1] == C.FDR_COUNTS
    monkeypatch.setenv(C.CONFIRM_ENV, json.dumps([key]))
    rec_fn = table(good={6, 7, 8, 9, 10, 11})
    rec = C.calibrate_fdr_counts(SUITE, STREAM, "k", 0.05, rec_fn, 6, log_dir=tmp_path)
    assert rec["block"] == 2 and rec_fn.calls == list(range(6, 12))
    skip_fn = table()
    skipped = C.calibrate_fdr_counts(SUITE, STREAM, "other", 0.05, skip_fn, 6, log_dir=tmp_path)
    assert skipped["status"] == C.SKIPPED_STATUS and skip_fn.calls == []
    monkeypatch.delenv(C.CONFIRM_ENV)
    final = {r["key"]: r["status"] for r in C.merge_shards(SUITE, tmp_path) if "key" in r}
    assert final == {"k": "unconfirmed", "other": "pass"}


# --------------------------------------------------------------------------- calibration


def _binomial_records(fdr: float, n: int = 120) -> Any:
    def records(i: int) -> dict[str, Any]:
        g = C.replicate_rng(f"{STREAM}.oc{fdr}", i)
        return {"k": cnt(n, int(g.binomial(n, fdr)))}

    return records


@pytest.mark.calibration
@pytest.mark.calsuite("novel_orfs")
@pytest.mark.parametrize(("fdr", "want"), [(0.03, "pass"), (0.15, "fail")])
def test_operating_characteristic(fdr: float, want: str, teeth_log: Path) -> None:
    """At R = 10 (the fast gate R) and 120 calls per replicate, a procedure at true FDR 0.03
    passes q = 0.05 and one at 0.15 (attack a's regime) is a confirmed fail."""
    rec = C.calibrate_fdr_counts(SUITE, f"{STREAM}.oc{fdr}", "k", 0.05, _binomial_records(fdr),
                                 10, log_dir=teeth_log)  # fmt: skip
    assert rec["status"] == want, rec
    assert np.isfinite(rec["bound"])
