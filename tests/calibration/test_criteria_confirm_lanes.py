"""M11-03c PR 2: block 2 of ``--merge`` split over L lanes gives the sequential confirm run's
records, verdicts and exit (R3: scheduling only; docs/specs/M11-03c-calibration-full-sharding.md,
"PR 2" and acceptance criterion 7).

``GIPPA_CALIBRATION_CONFIRM_SLICE=k/L`` is read only in the confirm branches of ``calibrate`` and
``_calibrate_counts`` (``calibrate_pooled``, ``calibrate_fdr_counts``): there a listed check runs
``shard_replicates(R, L, k)`` of block 2 (its own R) and logs it as a block-2 ``shard`` record; an
empty slice (R < L) logs nothing. ``GIPPA_CALIBRATION_CONFIRM_LANES=L`` makes ``--merge`` start L
block-2 pytests, each with ``--junitxml=<GIPPA_CALIBRATION_CONFIRM_JUNIT>.lane<k>.junit.xml`` on
the command line, and write ``<log>/<suite>.confirm_parts.json`` (the JUnit names it launched).
A bad ``CONFIRM_LANES`` is ``MERGE_ERROR``, never a silent fallback.

Part 1 (in process, fast): the slice in each confirm branch, the empty slice, unchanged readers
(``replicate_range``, ``owns_direct_checks``, ``shard_spec``), the slice ignored outside the
confirm run, and a Hypothesis property: the L slices partition block 2 exactly.

Part 2 (subprocess, the ``test_criteria_shard_equivalence.py`` pattern on
``_confirm_lanes_fixture_suite.py``): block 1 as 2 shards, then ``--merge`` with
``CONFIRM_LANES`` unset (today's sequential confirm run, the reference) and L = 1, 2, 3, plus the
unsharded run. Unmarked: runs on every commit (tiny R, synthetic checks).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import _criteria as C
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

SLICE_ENV = "GIPPA_CALIBRATION_CONFIRM_SLICE"
LANES_ENV = "GIPPA_CALIBRATION_CONFIRM_LANES"
JUNIT_ENV = "GIPPA_CALIBRATION_CONFIRM_JUNIT"
NEW_VARS = (SLICE_ENV, LANES_ENV, JUNIT_ENV)

HERE = Path(__file__).resolve().parent
FIXTURE = HERE / "_confirm_lanes_fixture_suite.py"
SUITE = "discordance"
FX = "discordance.confirm_lanes_fixture"
N_SHARDS = 2  # block-1 shards; L differs from it on purpose (L = 3)
JUNIT_BASE = "task-5.merge"  # <GIPPA_CALIBRATION_CONFIRM_JUNIT> = <log>/task-5.merge


def _clean(mp: pytest.MonkeyPatch) -> None:
    for var in (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV, *NEW_VARS, "PYTEST_ADDOPTS"):
        mp.delenv(var, raising=False)


# =========================================================================== part 1: in process


def _uniform(rng: Any) -> float:
    """p in [0, 0.04]: the level check fails in both blocks (so block 2 decides the record)."""
    return 0.04 * float(rng.uniform())


def _counts(stream: str) -> Any:
    def rec(i: int) -> dict[str, Any]:
        rng = C.replicate_rng(stream, i)
        n = int(rng.integers(5, 12))
        f = int(rng.integers(n // 2 + 1, n))  # F1 < 0.5, FDP > 0.5: both blocks fail
        return {"f1": {"tp_p": n - f, "n_pred": n, "tp_r": n - f, "n_truth": n},
                "fdr": {"n_called": n, "n_false": f}}  # fmt: skip

    return rec


def _run_check(kind: str, stream: str, R: int, log: Path) -> dict[str, Any]:
    if kind == "level":
        return C.calibrate(SUITE, "level", 0.05, _uniform, R, stream=stream, log_dir=log)
    if kind == "pooled":
        return C.calibrate_pooled(SUITE, stream, "f1", 0.5, _counts(stream), R, log_dir=log)
    return C.calibrate_fdr_counts(SUITE, stream, "fdr", 0.3, _counts(stream), R, log_dir=log)


def _confirm_key(kind: str, stream: str, R: int) -> list[Any]:
    if kind == "level":
        return [stream, "level", 0.05, R]
    if kind == "pooled":
        return [stream, C.POOLED_F1, 0.5, R, "f1"]
    return [stream, C.FDR_COUNTS, 0.3, R, "fdr"]


def _logged(log: Path) -> list[dict[str, Any]]:
    path = log / f"{SUITE}.jsonl"
    return C.read_results(path) if path.exists() else []


def _lanes(kind: str, stream: str, R: int, L: int, root: Path) -> tuple[list[Any], list[Any]]:
    """(sequential block-2 record, records of L sliced confirm lanes) under CONFIRM."""
    with pytest.MonkeyPatch.context() as mp:
        _clean(mp)
        mp.setenv(C.CONFIRM_ENV, json.dumps([_confirm_key(kind, stream, R)]))
        _run_check(kind, stream, R, root / "seq")
        for k in range(L):
            mp.setenv(SLICE_ENV, f"{k}/{L}")
            _run_check(kind, stream, R, root / "lanes")
    return _logged(root / "seq"), _logged(root / "lanes")


def _expected_slices(R: int, L: int) -> list[tuple[int, int]]:
    out = []
    for k in range(L):
        s = C.shard_replicates(R, L, k)
        if len(s):
            out.append((R + s.start, R + s.stop))
    return out


KINDS = ("level", "pooled", "fdr_counts")


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize(("R", "L"), [(5, 2), (5, 3), (7, 3), (4, 1), (2, 3), (1, 2)])
def test_confirm_slice_runs_only_its_share_of_block_2(
    kind: str, R: int, L: int, tmp_path: Path
) -> None:
    """Each lane logs one block-2 ``shard`` record of exactly ``shard_replicates(R, L, k)``
    (offset by R), whose values are byte-identical to that slice of the sequential confirm run;
    empty slices (R < L) log nothing (no ``start == stop`` record for ``_collect`` to reject)."""
    stream = f"{SUITE}.confirm_slice.{kind}.R{R}"
    (seq,), lanes = _lanes(kind, stream, R, L, tmp_path)
    assert (seq["block"], seq["start"], seq["stop"]) == (2, R, 2 * R)
    assert [(r["start"], r["stop"]) for r in lanes] == _expected_slices(R, L)
    assert all(r["block"] == 2 and r["status"] == C.SHARD_STATUS for r in lanes), lanes
    joined = [v for r in lanes for v in r["values"]]
    assert json.dumps(joined) == json.dumps(seq["values"])  # byte for byte, in index order
    assert all(r["seed"] == seq["seed"] for r in lanes)  # seed of the block's first replicate


@pytest.mark.parametrize("kind", KINDS)
def test_lane_records_merge_to_the_sequential_verdict(kind: str, tmp_path: Path) -> None:
    """Block 1 logged whole, block 2 as 3 lanes: ``merge_shards`` gives the sequential record."""
    R, L = 5, 3
    stream = f"{SUITE}.confirm_merge.{kind}"
    with pytest.MonkeyPatch.context() as mp:
        _clean(mp)
        mp.setenv(C.SHARD_ENV, f"0:{R}")  # block 1 as one shard record (no verdict, no block 2)
        for d in ("seq", "lanes"):
            _run_check(kind, stream, R, tmp_path / d)
    _seq, _lane = _lanes(kind, stream, R, L, tmp_path)
    ref = C.merge_shards(SUITE, tmp_path / "seq")
    got = C.merge_shards(SUITE, tmp_path / "lanes")
    assert [(r["status"], r["block"]) for r in ref] == [("fail", 2)]  # block 2 decided it
    fields = ("status", "observed", "bound", "seed", "ok", "R", "block", "start", "stop")
    assert [{f: r[f] for f in fields} for r in got] == [{f: r[f] for f in fields} for r in ref]
    assert [json.dumps(r["values"]) for r in got] == [json.dumps(r["values"]) for r in ref]


@pytest.mark.parametrize("kind", KINDS)
def test_unlisted_check_is_still_skipped_under_a_slice(kind: str, tmp_path: Path) -> None:
    with pytest.MonkeyPatch.context() as mp:
        _clean(mp)
        mp.setenv(C.CONFIRM_ENV, json.dumps([["other.stream", "level", 0.05, 5]]))
        mp.setenv(SLICE_ENV, "1/3")
        rec = _run_check(kind, f"{SUITE}.confirm_unlisted.{kind}", 5, tmp_path)
    assert rec["status"] == C.SKIPPED_STATUS
    assert _logged(tmp_path) == []


def test_slice_leaves_every_other_reader_unchanged() -> None:
    """The slice is read only in the two confirm branches: ``replicate_range`` (``_gc``'s
    vacuity guard reads replicate 0), ``owns_direct_checks``, ``shard_spec``, ``replicate_block``
    see what the sequential confirm run sees."""
    with pytest.MonkeyPatch.context() as mp:
        _clean(mp)
        mp.setenv(C.CONFIRM_ENV, json.dumps([["s", "level", 0.05, 7]]))
        ref = (C.replicate_range(7), C.owns_direct_checks(), C.shard_spec(), C.replicate_block())
        for k in range(3):
            mp.setenv(SLICE_ENV, f"{k}/3")
            got = (C.replicate_range(7), C.owns_direct_checks(), C.shard_spec(),
                   C.replicate_block())  # fmt: skip
            assert got == ref == (range(7), False, None, 1)
            assert next(iter(C.replicate_range(7))) == 0


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("shard", [None, "0:5@1/2"])
def test_slice_outside_a_confirm_run_changes_nothing(
    kind: str, shard: str | None, tmp_path: Path
) -> None:
    """Without ``GIPPA_CALIBRATION_CONFIRM`` the slice variable is ignored (unsharded: both
    blocks; block-1 shard: its own range), so behaviour is byte-unchanged."""
    stream = f"{SUITE}.confirm_ignored.{kind}"
    out = {}
    for name, slice_ in (("without", None), ("with", "2/3")):
        with pytest.MonkeyPatch.context() as mp:
            _clean(mp)
            if shard is not None:
                mp.setenv(C.SHARD_ENV, shard)
            if slice_ is not None:
                mp.setenv(SLICE_ENV, slice_)
            _run_check(kind, stream, 5, tmp_path / name)
        recs = _logged(tmp_path / name)
        for r in recs:
            r.pop("duration")
        out[name] = recs
    assert out["with"] == out["without"]
    assert out["without"], "the check logged nothing"


@settings(max_examples=40, deadline=None)
@given(R=st.integers(1, 12), L=st.integers(1, 6), kind=st.sampled_from(KINDS))
def test_property_lane_slices_partition_block_2(R: int, L: int, kind: str) -> None:
    """For every R and L: the non-empty lane records cover R..2R-1 exactly once, contiguous, in
    lane order, with the sequential run's values (R3: the same replicates, the same seeds)."""
    stream = f"{SUITE}.confirm_prop.{kind}.R{R}"
    with tempfile.TemporaryDirectory() as d:
        (seq,), lanes = _lanes(kind, stream, R, L, Path(d))
    covered = [i for r in lanes for i in range(r["start"], r["stop"])]
    assert covered == list(range(R, 2 * R))
    assert all(r["stop"] > r["start"] for r in lanes)
    assert len(lanes) == min(R, L)
    assert json.dumps([v for r in lanes for v in r["values"]]) == json.dumps(seq["values"])


# =========================================================================== part 2: --merge

SCENARIOS = {
    # every pending kind, both strict-xfail confirmed fails, a vacuous block 1, R = 2 checks
    "ok": "pass_b1 or cal_unconfirmed or pooled_unconfirmed or fdrc_unconfirmed or strict_fail"
    " or vacuous or r2_",
    "fail": "pass_b1 or unmarked_fail or unmarked_fdrc_fail",  # unmarked confirmed fails
    "xpass": "pass_b1 or xpass_strict",  # block 1 fails, block 2 passes under a strict xfail
    "vacuous": "pass_b1 or vacuous",  # nothing pending: no block-2 lane at all
}
WANT = {"ok": 0, "fail": 1, "xpass": 1, "vacuous": 0}  # unsharded pytest exit (0 / non-zero)
WANT_MERGE = {"ok": 0, "fail": C.MERGE_FAIL, "xpass": C.MERGE_FAIL, "vacuous": 0}
PENDING = {"ok": True, "fail": True, "xpass": True, "vacuous": False}
#: scenarios outside the equivalence matrix: ``raise`` = a pending check that logs its block-2
#: slice and then raises RuntimeError in the process whose slice holds the block's last replicate
#: (``_gc``'s reach check; review S2): one lane of several, the sequential run, the unsharded run
EXTRA = {"raise": "pass_b1 or raise_after_log or cal_unconfirmed"}
#: fields compared between a lane merge and the reference (``duration`` aside); ``values`` apart
FIELDS = ("status", "observed", "bound", "seed", "ok", "R", "block", "start", "stop", "m", "key",
          "test", "xfail")  # fmt: skip


def _args(selection: str) -> list[str]:
    return [str(FIXTURE), "-p", "no:cacheprovider", "-q", "-m", "calibration", "-k", selection]


def _env(log: Path) -> dict[str, str]:
    drop = (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV, *NEW_VARS, "PYTEST_ADDOPTS",
            "PYTEST_XDIST_WORKER")  # fmt: skip
    keep = {k: v for k, v in os.environ.items() if k not in drop}
    return {**keep, C.LOG_DIR_ENV: str(log), "GIPPA_CONFIRM_LANES_FIXTURE": "1"}


def _setup(mp: pytest.MonkeyPatch, log: Path, selection: str) -> None:
    mp.setattr(C, "pytest_args", lambda _suite, _full: _args(selection))
    _clean(mp)
    for k, v in _env(log).items():
        mp.setenv(k, v)
    mp.chdir(C.REPO)


_BASE: dict[str, tuple[Path, int, Path]] = {}


@pytest.fixture(scope="module")
def base_root(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return tmp_path_factory.mktemp("confirm_lanes_base")


def _base(scenario: str, root: Path) -> tuple[Path, int, Path]:
    """(unsharded log, unsharded pytest exit, block-1 shard log) of ``scenario``, built once per
    process; every merge works on a copy of the shard log."""
    if scenario not in _BASE:
        sel = {**SCENARIOS, **EXTRA}[scenario]
        one, shards = root / scenario / "one", root / scenario / "shards"
        cmd = [sys.executable, "-m", "pytest", *_args(sel)]
        rc = subprocess.run(cmd, env=_env(one), cwd=C.REPO, check=False).returncode
        with pytest.MonkeyPatch.context() as mp:
            _setup(mp, shards, sel)
            rcs = [C.main(["--suite", SUITE, "--fast", "--shard", f"{k}/{N_SHARDS}"])
                   for k in range(N_SHARDS)]  # fmt: skip
        assert rcs == [0] * N_SHARDS, rcs  # a shard has no verdict of its own
        _BASE[scenario] = (one, rc, shards)
    return _BASE[scenario]


def _merge(scenario: str, root: Path, dest: Path, env: dict[str, str]) -> int:
    _one, _rc, shards = _base(scenario, root)
    shutil.copytree(shards, dest)
    with pytest.MonkeyPatch.context() as mp:
        _setup(mp, dest, {**SCENARIOS, **EXTRA}[scenario])
        for k, v in env.items():
            mp.setenv(k, v)
        return C.main(["--suite", SUITE, "--fast", "--merge"])


def _lane_env(log: Path, L: int | None) -> dict[str, str]:
    """The job's merge environment (spec step 4): the merge's own JUnit in PYTEST_ADDOPTS, plus
    the lanes; ``None`` = CONFIRM_LANES unset (today)."""
    env = {"PYTEST_ADDOPTS": f"--junitxml={log / 'task-5.merge.junit.xml'}"}
    if L is not None:
        env |= {LANES_ENV: str(L), JUNIT_ENV: str(log / JUNIT_BASE)}
    return env


def _by_check(records: list[dict[str, Any]]) -> dict[tuple[Any, ...], dict[str, Any]]:
    out = {}
    for r in records:
        key = (r["stream"], r["criterion"], r["level"], r["R"], str(r.get("key")))
        assert key not in out, key
        out[key] = r
    return out


def _block2(log: Path) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {}
    for r in _logged(log):
        if int(r.get("block", 1)) == 2:
            out.setdefault(r["stream"], []).append((int(r["start"]), int(r["stop"])))
    return {k: sorted(v) for k, v in out.items()}


def _lines(path: Path) -> list[str]:
    return path.read_text().splitlines() if path.exists() else []


@pytest.fixture(scope="module")
def references(base_root: Path) -> dict[str, tuple[int, Path]]:
    """Merge exit and log of today's sequential confirm run (CONFIRM_LANES unset), per scenario."""
    out = {}
    for scenario in SCENARIOS:
        dest = base_root / scenario / "ref"
        out[scenario] = (_merge(scenario, base_root, dest, _lane_env(dest, None)), dest)
    return out


@pytest.mark.parametrize("L", [1, 2, 3])
@pytest.mark.parametrize("scenario", sorted(SCENARIOS))
def test_block_2_on_lanes_equals_the_sequential_and_unsharded_runs(
    scenario: str,
    L: int,
    tmp_path: Path,
    base_root: Path,
    references: dict[str, tuple[int, Path]],
) -> None:
    rc_ref, ref_log = references[scenario]
    one, rc_one, shards = _base(scenario, base_root)
    log = tmp_path / "log"
    rc = _merge(scenario, base_root, log, _lane_env(log, L))

    # exits: the unsharded run, today's merge and the lane merge agree
    assert (rc_one != 0) == (WANT[scenario] != 0), rc_one
    assert rc_ref == WANT_MERGE[scenario], rc_ref
    assert rc == rc_ref, (rc, rc_ref)

    # per-check verdict records: byte-identical values, same fields, same outcome
    ref = _by_check(C.read_results(ref_log / f"{SUITE}.merged.jsonl"))
    got = _by_check(C.read_results(log / f"{SUITE}.merged.jsonl"))
    assert got.keys() == ref.keys() and ref
    for key, r in ref.items():
        g = got[key]
        assert json.dumps(g["values"]) == json.dumps(r["values"]), key
        assert {f: g.get(f) for f in (*FIELDS, "outcome")} == {
            f: r.get(f) for f in (*FIELDS, "outcome")
        }, key
    unsharded = _by_check(C.merge_shards(SUITE, one))
    assert unsharded.keys() == got.keys()
    for key, u in unsharded.items():
        assert json.dumps(got[key]["values"]) == json.dumps(u["values"]), key
        assert {f: got[key].get(f) for f in FIELDS} == {f: u.get(f) for f in FIELDS}, key

    # where block 2 ran: L slices per pending check (R < L: fewer), never the whole block twice
    want_b2 = {s: [(r["start"], r["stop"])] for s, r in _by_check_stream(ref_log).items()}
    lanes_b2 = _block2(log)
    assert lanes_b2.keys() == _block2(ref_log).keys() == want_b2.keys()
    for stream, (whole,) in want_b2.items():
        R = whole[1] - whole[0]
        want = _expected_slices(R, L)
        assert lanes_b2[stream] == want, (stream, lanes_b2[stream], want)

    # block-2 lanes write no shard exit and no ``deferred`` outcome
    assert _lines(log / f"{SUITE}.shards.jsonl") == _lines(shards / f"{SUITE}.shards.jsonl")
    new_outcomes = _lines(log / f"{SUITE}.outcomes.jsonl")[
        len(_lines(shards / f"{SUITE}.outcomes.jsonl")) :
    ]
    assert not [o for o in new_outcomes if json.loads(o).get("outcome") == C.DEFERRED]

    # the L block-2 pytests write L distinct JUnit files, listed in confirm_parts.json; L = 1 too
    # (lead decision on review S1: one process with its own lane-0 JUnit, so the job's ``parts``
    # follows one rule whenever CONFIRM_LANES is set)
    names = [f"{JUNIT_BASE}.lane{k}.junit.xml" for k in range(L)]
    parts = json.loads((log / f"{SUITE}.confirm_parts.json").read_text())
    assert isinstance(parts, list)
    if PENDING[scenario]:
        assert sorted(parts) == sorted(names), parts
        assert len(set(parts)) == L
        for name in names:
            assert (log / name).is_file(), name
        # the per-lane --junitxml (command line) wins over the inherited PYTEST_ADDOPTS one
        assert not (log / "task-5.merge.junit.xml").exists()
    else:
        assert parts == [], parts
        assert not list(log.glob(f"{JUNIT_BASE}.lane*.junit.xml"))


@pytest.mark.parametrize("scenario", sorted(SCENARIOS))
def test_unset_lanes_keeps_todays_junit_and_writes_no_parts(
    scenario: str, references: dict[str, tuple[int, Path]]
) -> None:
    """CONFIRM_LANES unset: today's behaviour (one process, the PYTEST_ADDOPTS JUnit, no
    ``confirm_parts.json``, no lane JUnit)."""
    _rc, ref_log = references[scenario]
    assert not (ref_log / f"{SUITE}.confirm_parts.json").exists()
    assert not list(ref_log.glob(f"{JUNIT_BASE}.lane*.junit.xml"))
    assert (ref_log / "task-5.merge.junit.xml").is_file() == PENDING[scenario]


def _by_check_stream(log: Path) -> dict[str, dict[str, Any]]:
    """Block-2 records of the sequential confirm run, by stream (one per pending check)."""
    out: dict[str, dict[str, Any]] = {}
    for r in _logged(log):
        if int(r.get("block", 1)) == 2:
            assert r["stream"] not in out, r["stream"]
            out[r["stream"]] = r
    return out


def test_ok_scenario_outcomes_and_statuses(
    base_root: Path, references: dict[str, tuple[int, Path]], tmp_path: Path
) -> None:
    """Both strict-xfail directions and vacuity under L = 3 lanes (ADR-0057, rule 6)."""
    log = tmp_path / "ok"
    assert _merge("ok", base_root, log, _lane_env(log, 3)) == 0
    merged = {r["stream"]: r for r in C.read_results(log / f"{SUITE}.merged.jsonl")}
    for name in ("strict_fail_cal", "strict_fail_pooled", "strict_fail_fdrc",
                 "r2_fdrc_strict_fail"):  # fmt: skip
        assert merged[f"{FX}.{name}"]["status"] == "fail", name
        assert merged[f"{FX}.{name}"]["outcome"] == "xfail", name
    for name in ("cal_unconfirmed", "pooled_unconfirmed", "fdrc_unconfirmed", "r2_cal",
                 "r2_pooled"):  # fmt: skip
        assert merged[f"{FX}.{name}"]["status"] == "unconfirmed", name
        assert merged[f"{FX}.{name}"]["block"] == 2, name
    vac = merged[f"{FX}.vacuous"]
    assert (vac["status"], vac["block"]) == (C.VACUOUS_STATUS, 1)
    assert f"{FX}.vacuous" not in _block2(log)  # a vacuous block 1 never reaches block 2
    assert merged[f"{FX}.pass_b1"]["status"] == "pass"
    assert f"{FX}.pass_b1" not in _block2(log)


def test_xpass_scenario_is_a_strict_xpass_on_lanes(base_root: Path, tmp_path: Path) -> None:
    """Block 1 fails, block 2 passes on the lanes (``unconfirmed``) under a strict xfail:
    XPASS(strict), ``--merge`` exit 1, as unsharded (a lane's partial block-2 record must never
    log an ``xfailed`` outcome that would turn it into XFAIL)."""
    log = tmp_path / "xp"
    assert _merge("xpass", base_root, log, _lane_env(log, 3)) == C.MERGE_FAIL
    merged = {r["stream"]: r for r in C.read_results(log / f"{SUITE}.merged.jsonl")}
    for name in ("xpass_strict_cal", "xpass_strict_fdrc"):
        assert merged[f"{FX}.{name}"]["status"] == "unconfirmed"
        assert merged[f"{FX}.{name}"]["outcome"] == "xpass"
    outcomes = [json.loads(x) for x in _lines(log / f"{SUITE}.outcomes.jsonl")]
    assert not [o for o in outcomes if o.get("outcome") == C.XFAILED], outcomes
    assert len(_block2(log)[f"{FX}.xpass_strict_cal"]) == 3


def test_r2_with_three_lanes_logs_no_empty_record(base_root: Path, tmp_path: Path) -> None:
    """R = 2, L = 3: lane 2's slice is empty and logs nothing (else ``_collect`` raises
    "malformed record" and the merge is a spurious MERGE_ERROR)."""
    log = tmp_path / "r2"
    assert _merge("ok", base_root, log, _lane_env(log, 3)) == 0
    for name in ("r2_cal", "r2_pooled", "r2_fdrc_strict_fail"):
        assert _block2(log)[f"{FX}.{name}"] == [(2, 3), (3, 4)], name
    assert all(r["stop"] > r["start"] for r in _logged(log))


@pytest.mark.parametrize(
    "env",
    [{LANES_ENV: "x"}, {LANES_ENV: "0"}, {LANES_ENV: "-1"}, {LANES_ENV: "2.5"}],
    ids=["non-integer", "zero", "negative", "float"],
)
def test_bad_confirm_lanes_is_a_merge_error(
    env: dict[str, str], base_root: Path, tmp_path: Path
) -> None:
    log = tmp_path / "bad"
    rc = _merge("ok", base_root, log, {**env, JUNIT_ENV: str(log / JUNIT_BASE)})
    assert rc == C.MERGE_ERROR, rc


def test_confirm_lanes_without_confirm_junit_is_a_merge_error(
    base_root: Path, tmp_path: Path
) -> None:
    """Never a silent fallback to one lane (statistician S5)."""
    log = tmp_path / "nojunit"
    assert _merge("ok", base_root, log, {LANES_ENV: "2"}) == C.MERGE_ERROR


def test_a_dead_block_2_lane_is_a_merge_error(base_root: Path, tmp_path: Path) -> None:
    """A lane that dies before logging leaves block 2 incomplete: MERGE_ERROR, never a verdict
    (fail closed; here every lane's pytest is refused by a bad option)."""
    log = tmp_path / "dead"
    env = {**_lane_env(log, 2), "PYTEST_ADDOPTS": "--no-such-option"}
    assert _merge("ok", base_root, log, env) == C.MERGE_ERROR
    parts = json.loads((log / f"{SUITE}.confirm_parts.json").read_text())
    assert sorted(parts) == [f"{JUNIT_BASE}.lane{k}.junit.xml" for k in range(2)]  # launched


@pytest.mark.parametrize("L", [2, 3])
def test_a_lane_that_raises_after_logging_is_a_merge_error(
    L: int, base_root: Path, tmp_path: Path
) -> None:
    """Review S2: one lane logs its slice of a pending check, then raises RuntimeError (outside
    ``raises=AssertionError``); the other lanes succeed. Every replicate is logged and the check
    is ``unconfirmed``, so the records alone say exit 0: only the folded lane exit carries the
    failure. ``--merge`` must give MERGE_ERROR, as the sequential confirm run does."""
    _one, rc_one, _shards = _base("raise", base_root)
    assert rc_one != 0  # the unsharded run fails too (the raise is in its block 2)
    ref = tmp_path / "ref"
    assert _merge("raise", base_root, ref, _lane_env(ref, None)) == C.MERGE_ERROR
    log = tmp_path / "log"
    assert _merge("raise", base_root, log, _lane_env(log, L)) == C.MERGE_ERROR
    stream = f"{FX}.raise_after_log"
    R = 5
    assert _block2(log)[stream] == _expected_slices(R, L)  # every slice logged: block 2 complete
    raised = [k for k in range(L) if R + C.shard_replicates(R, L, k).stop == 2 * R]
    assert raised == [L - 1]  # exactly one lane raised; the others passed
    parts = json.loads((log / f"{SUITE}.confirm_parts.json").read_text())
    assert sorted(parts) == [f"{JUNIT_BASE}.lane{k}.junit.xml" for k in range(L)]


@pytest.mark.parametrize("L", [2, 3])
def test_one_dead_lane_of_several_is_a_merge_error(L: int, base_root: Path, tmp_path: Path) -> None:
    """Statistician's gap: only lane 1 dies (killed before it logs, exit 137 as after an OOM
    kill); the other lanes finish and log. Block 2 is incomplete: MERGE_ERROR, never a verdict."""
    log = tmp_path / "onedead"
    env = {**_lane_env(log, L), "GIPPA_CONFIRM_LANES_FIXTURE_DIE": "1"}
    assert _merge("ok", base_root, log, env) == C.MERGE_ERROR
    parts = json.loads((log / f"{SUITE}.confirm_parts.json").read_text())
    assert sorted(parts) == [f"{JUNIT_BASE}.lane{k}.junit.xml" for k in range(L)]
    assert not (log / f"{JUNIT_BASE}.lane1.junit.xml").exists()  # lane 1 died
    for k in set(range(L)) - {1}:
        assert (log / f"{JUNIT_BASE}.lane{k}.junit.xml").is_file(), k  # the others ran
    lane0 = C.shard_replicates(5, L, 0)
    assert (5 + lane0.start, 5 + lane0.stop) in _block2(log)[f"{FX}.cal_unconfirmed"]


@pytest.mark.parametrize("scenario", sorted(SCENARIOS))
def test_a_stray_confirm_slice_without_lanes_gives_the_reference(
    scenario: str, base_root: Path, tmp_path: Path, references: dict[str, tuple[int, Path]]
) -> None:
    """Review S4: CONFIRM_LANES unset and a stray ``GIPPA_CALIBRATION_CONFIRM_SLICE=1/3`` in the
    merge's environment: the sequential confirm child must not inherit it, so block 2 runs whole
    and the merge gives the reference exit and records."""
    rc_ref, ref_log = references[scenario]
    log = tmp_path / "stray"
    assert _merge(scenario, base_root, log, {**_lane_env(log, None), SLICE_ENV: "1/3"}) == rc_ref
    ref = _by_check(C.read_results(ref_log / f"{SUITE}.merged.jsonl"))
    got = _by_check(C.read_results(log / f"{SUITE}.merged.jsonl"))
    assert got.keys() == ref.keys() and ref
    for key, r in ref.items():
        assert json.dumps(got[key]["values"]) == json.dumps(r["values"]), key
        assert {f: got[key].get(f) for f in (*FIELDS, "outcome")} == {
            f: r.get(f) for f in (*FIELDS, "outcome")
        }, key
    assert _block2(log) == _block2(ref_log)
