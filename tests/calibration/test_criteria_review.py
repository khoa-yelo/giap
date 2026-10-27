"""Lead decisions after the statistician review of the calibration harness (PR #54, BLOCK).

Pins items 1-4 of "Lead decisions after statistician review" in
docs/specs/M4-05-calibration-harness.md and ADR-0015 (two-block confirmation):

1. Fail closed: ``check_level`` / ``check_fwer`` raise ValueError on non-finite or out-of-[0, 1]
   input; ``check_fdr`` raises on FDP outside [0, 1] and on R < 2. (The all-NaN and p > 1 cases of
   ``check_level`` live in attacks/test_attack_calibration_harness.py and are not repeated here.)
2. Sharding: suites read ``GIPPA_CALIBRATION_REPLICATES=start:stop`` through the harness and run only
   that range; ``merge_shards`` (CLI ``--merge``) combines shard records and evaluates the criterion
   on the full R. Without the variable behaviour is unchanged.
3. ADR-0015: a check failing on replicates 0..R-1 is re-run on R..2R-1 with the same seed formula
   ``SeedSequence([calibration_base, crc32(stream), rep])``; the gate fails only if both blocks
   fail. Records carry ``block`` (1 | 2) and ``status`` (pass | unconfirmed | fail).
4. Suite assignment is explicit: every calibration test runs in exactly one array task.

Harness entry points pinned here (test-author's reading of the decisions):

- ``replicate_range(R) -> range``: ``range(start, stop)`` from the variable, else ``range(R)``.
- ``calibrate(suite, criterion, level, fn, R, *, stream=None, log_dir=None) -> Mapping``: runs
  ``fn(replicate_rng(stream or suite, i))`` (one scalar per replicate) for the replicates of the
  current shard, applies ``criterion`` in {"level", "uniform", "fdr", "fwer"} at ``level``, applies
  the two-block confirmation and logs to ``<log_dir>/<suite>.jsonl``; returns the final record.
- ``merge_shards(suite, log_dir=None) -> list[Mapping]``: one verdict record per check.

Unmarked: these run on every commit.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
from collections.abc import Callable
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

FAST = settings(max_examples=50, deadline=None)
BAD = [float("nan"), float("inf"), float("-inf"), -0.1, -1e-12, 1.0 + 1e-12, 1.5]

# --------------------------------------------------------------------------- 1. fail closed


@pytest.mark.parametrize("bad", [*BAD[1:5], 1.5], ids=repr)
def test_check_level_rejects_bad_p_among_valid(bad: float) -> None:
    p = [0.5] * 99 + [bad]
    with pytest.raises(ValueError):
        C.check_level(p, 0.05)


def test_check_level_rejects_single_nan_among_valid() -> None:
    with pytest.raises(ValueError):
        C.check_level([0.5] * 99 + [float("nan")], 0.05)


def test_check_level_accepts_closed_unit_interval() -> None:
    res = C.check_level([0.0] * 2 + [1.0] * 98, 0.05)
    assert res.ok and res.observed == pytest.approx(0.02)


@pytest.mark.parametrize("bad", [*BAD, 2.0], ids=repr)
def test_check_fwer_rejects_bad_input(bad: float) -> None:
    with pytest.raises(ValueError):
        C.check_fwer([0.0] * 99 + [bad])


def test_check_fwer_accepts_bools_and_zero_one() -> None:
    assert C.check_fwer([True] * 3 + [False] * 97).observed == pytest.approx(0.03)
    assert C.check_fwer([1.0] * 3 + [0.0] * 97).observed == pytest.approx(0.03)


@pytest.mark.parametrize("bad", BAD, ids=repr)
def test_check_fdr_rejects_fdp_outside_unit_interval(bad: float) -> None:
    with pytest.raises(ValueError):
        C.check_fdr([0.01] * 49 + [bad], 0.05)


def test_check_fdr_rejects_single_replicate() -> None:
    with pytest.raises(ValueError):
        C.check_fdr([0.0], 0.05)


def test_check_fdr_accepts_two_replicates_and_unit_bounds() -> None:
    res = C.check_fdr([0.0, 0.0], 0.05)
    assert res.ok and res.R == 2
    res = C.check_fdr([0.0, 1.0], 0.05)
    assert res.R == 2 and res.observed == pytest.approx(0.5)


@FAST
@given(
    st.lists(st.floats(0.0, 1.0), min_size=2, max_size=50),
    st.sampled_from(BAD),
    st.integers(0, 49),
)
def test_any_bad_entry_raises_property(good: list[float], bad: float, pos: int) -> None:
    x = list(good)
    x.insert(pos % (len(x) + 1), bad)
    for check in (
        lambda v: C.check_level(v, 0.05),
        lambda v: C.check_fwer(v),
        lambda v: C.check_fdr(v, 0.05),
    ):
        with pytest.raises(ValueError):
            check(x)


# --------------------------------------------------------------------------- helpers


class Recorder:
    """Replicate function returning ``values(call_index, rng)`` and recording one draw per call."""

    def __init__(self, value: Callable[[int, float], float]) -> None:
        self.value = value
        self.draws: list[float] = []

    def __call__(self, rng: np.random.Generator) -> float:
        u = float(rng.uniform())
        out = self.value(len(self.draws), u)
        self.draws.append(u)
        return out


def _expected_draws(stream: str, reps: range) -> list[float]:
    return [float(C.replicate_rng(stream, i).uniform()) for i in reps]


def _records(log_dir: Path, suite: str) -> list[dict[str, Any]]:
    return C.read_results(log_dir / f"{suite}.jsonl")


# --------------------------------------------------------------------------- 2. sharding


def test_replicate_range_without_variable_is_full(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(C.SHARD_ENV, raising=False)
    assert C.replicate_range(10) == range(10)


def test_replicate_range_reads_variable(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(C.SHARD_ENV, "2:5")
    assert C.replicate_range(10) == range(2, 5)


@pytest.mark.parametrize("bad", ["5:2", "0:11", "-1:3", "a:b", "3", ""])
def test_replicate_range_rejects_bad_variable(monkeypatch: pytest.MonkeyPatch, bad: str) -> None:
    monkeypatch.setenv(C.SHARD_ENV, bad)
    with pytest.raises(ValueError):
        C.replicate_range(10)


def test_variable_matches_cli_shard_blocks(monkeypatch: pytest.MonkeyPatch) -> None:
    """``--shard K/N`` exports exactly the block that ``replicate_range`` then runs."""
    r = C.replicates("regions", full=True)
    for k in range(3):
        ns = C.parse_args(["--suite", "regions", "--shard", f"{k}/3"])
        monkeypatch.setenv(C.SHARD_ENV, f"{ns.replicates.start}:{ns.replicates.stop}")
        assert C.replicate_range(r) == C.shard_replicates(r, 3, k)


STREAM_FDR = "regions.review_shard"
R_SHARD = 10


def _fdr_fn() -> Recorder:
    return Recorder(lambda _i, u: 0.04 * u)  # FDP ~ U(0, 0.04): mean 0.02, always passes q = 0.05


def _run(log_dir: Path, fn: Recorder) -> Any:
    return C.calibrate("regions", "fdr", 0.05, fn, R_SHARD, stream=STREAM_FDR, log_dir=log_dir)


def test_shard_runs_only_its_range(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(C.SHARD_ENV, "4:7")
    fn = _fdr_fn()
    _run(tmp_path, fn)
    assert fn.draws == _expected_draws(STREAM_FDR, range(4, 7))


def test_unsharded_calibrate_runs_all_replicates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv(C.SHARD_ENV, raising=False)
    fn = _fdr_fn()
    rec = _run(tmp_path, fn)
    assert fn.draws == _expected_draws(STREAM_FDR, range(R_SHARD))
    assert rec["R"] == R_SHARD and rec["status"] == "pass" and rec["block"] == 1
    want = C.check_fdr([0.04 * u for u in fn.draws], 0.05)
    assert rec["observed"] == pytest.approx(want.observed)
    assert rec["bound"] == pytest.approx(want.bound)


@pytest.mark.parametrize("blocks", [["0:4", "4:10"], ["4:10", "0:4"], ["0:3", "3:6", "6:10"]])
def test_merged_shards_equal_unsharded_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, blocks: list[str]
) -> None:
    whole, sharded = tmp_path / "whole", tmp_path / "sharded"
    monkeypatch.delenv(C.SHARD_ENV, raising=False)
    ref = _run(whole, _fdr_fn())
    for b in blocks:
        monkeypatch.setenv(C.SHARD_ENV, b)
        _run(sharded, _fdr_fn())
    monkeypatch.delenv(C.SHARD_ENV, raising=False)
    merged = C.merge_shards("regions", log_dir=sharded)
    assert len(merged) == 1
    m = merged[0]
    for key in ("ok", "R", "seed", "status", "block", "criterion"):
        assert m[key] == ref[key], key
    assert m["level"] == pytest.approx(ref["level"])
    assert m["observed"] == pytest.approx(ref["observed"], rel=1e-12)
    assert m["bound"] == pytest.approx(ref["bound"], rel=1e-12)


def test_merge_refuses_missing_replicates(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A verdict needs the full R: shards 0:4 and 6:10 leave 4, 5 missing."""
    for b in ["0:4", "6:10"]:
        monkeypatch.setenv(C.SHARD_ENV, b)
        _run(tmp_path, _fdr_fn())
    monkeypatch.delenv(C.SHARD_ENV, raising=False)
    with pytest.raises(ValueError):
        C.merge_shards("regions", log_dir=tmp_path)


def test_cli_merge_flag(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    for b in ["0:5", "5:10"]:
        monkeypatch.setenv(C.SHARD_ENV, b)
        _run(tmp_path, _fdr_fn())
    monkeypatch.delenv(C.SHARD_ENV, raising=False)
    monkeypatch.setenv(C.LOG_DIR_ENV, str(tmp_path))
    ns = C.parse_args(["--suite", "regions", "--merge"])
    assert ns.merge is True
    assert C.parse_args(["--suite", "regions"]).merge is False

    def no_pytest(*_a: object, **_k: object) -> object:
        raise AssertionError("--merge must not run the suite")

    monkeypatch.setattr(C.subprocess, "run", no_pytest)
    assert C.main(["--suite", "regions", "--merge"]) == 0


# ------------------------------------------------------------------- 3. two-block (ADR-0015)

STREAM_2B = "regions.review_two_block"
R_2B = 20


def _calibrate_level(log_dir: Path, fn: Recorder) -> Any:
    return C.calibrate("regions", "level", 0.05, fn, R_2B, stream=STREAM_2B, log_dir=log_dir)


@pytest.fixture(autouse=True)
def _no_shard(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(C.SHARD_ENV, raising=False)


def test_passing_check_runs_one_block(tmp_path: Path) -> None:
    fn = Recorder(lambda _i, _u: 0.5)
    rec = _calibrate_level(tmp_path, fn)
    assert rec["status"] == "pass" and rec["block"] == 1 and rec["ok"] is True
    assert len(fn.draws) == R_2B
    recs = _records(tmp_path, "regions")
    assert all("block" in r and "status" in r for r in recs)
    assert {r["block"] for r in recs} == {1}


def test_failure_only_in_block_one_is_unconfirmed(tmp_path: Path) -> None:
    fn = Recorder(lambda i, _u: 0.0 if i < R_2B else 0.5)  # rejects every block-1 replicate
    rec = _calibrate_level(tmp_path, fn)
    assert rec["status"] == "unconfirmed"  # gate passes, visible in the summary
    assert rec["status"] != "fail"
    assert rec["block"] == 2 and rec["R"] == R_2B
    recs = _records(tmp_path, "regions")
    assert all("block" in r and "status" in r for r in recs)
    assert sorted({r["block"] for r in recs}) == [1, 2]
    b1 = [r for r in recs if r["block"] == 1]
    assert b1 and all(not r["ok"] for r in b1)


def test_failure_in_both_blocks_fails(tmp_path: Path) -> None:
    fn = Recorder(lambda _i, _u: 0.0)
    rec = _calibrate_level(tmp_path, fn)
    assert rec["status"] == "fail" and rec["ok"] is False and rec["block"] == 2
    recs = _records(tmp_path, "regions")
    assert sorted({r["block"] for r in recs}) == [1, 2]
    assert recs[-1]["status"] == "fail"


def test_second_block_uses_replicates_r_to_2r_with_same_seed_formula(tmp_path: Path) -> None:
    fn = Recorder(lambda _i, _u: 0.0)
    _calibrate_level(tmp_path, fn)
    assert len(fn.draws) == 2 * R_2B
    assert fn.draws == _expected_draws(STREAM_2B, range(2 * R_2B))
    seeds = {r["block"]: r["seed"] for r in _records(tmp_path, "regions")}
    assert seeds[1] != seeds[2]


def test_two_block_status_is_deterministic(tmp_path: Path) -> None:
    a = _calibrate_level(tmp_path / "a", Recorder(lambda i, _u: 0.0 if i < R_2B else 0.5))
    b = _calibrate_level(tmp_path / "b", Recorder(lambda i, _u: 0.0 if i < R_2B else 0.5))
    assert dict(a)["status"] == dict(b)["status"] and a["seed"] == b["seed"]


@FAST
@given(st.booleans(), st.booleans())
def test_status_truth_table(fail1: bool, fail2: bool) -> None:
    import tempfile

    def val(i: int, _u: float) -> float:
        bad = fail1 if i < R_2B else fail2
        return 0.0 if bad else 0.5

    with tempfile.TemporaryDirectory() as d:
        rec = _calibrate_level(Path(d), Recorder(val))
    want = "fail" if (fail1 and fail2) else "unconfirmed" if fail1 else "pass"
    assert rec["status"] == want


def test_summary_shows_status(tmp_path: Path) -> None:
    _calibrate_level(tmp_path, Recorder(lambda i, _u: 0.0 if i < R_2B else 0.5))
    md = C.summary_markdown(_records(tmp_path, "regions"))
    assert "unconfirmed" in md


# --------------------------------------------------------------------------- 4. suite assignment

_COLLECTOR = textwrap.dedent(
    """
    import json, sys
    import pytest

    class Rec:
        def __init__(self):
            self.ids = []
        def pytest_collection_finish(self, session):
            self.ids = [i.nodeid for i in session.items]

    out = {}
    for key, args in json.loads(sys.argv[1]).items():
        rec = Rec()
        code = pytest.main(
            ["--collect-only", "-q", "-p", "no:cacheprovider", *args],
            plugins=[rec],
        )
        out[key] = {"code": int(code), "ids": rec.ids}
    json.dump(out, open(sys.argv[2], "w"))
    """
)


@pytest.fixture(scope="module")
def selections(tmp_path_factory: pytest.TempPathFactory) -> dict[str, dict[str, Any]]:
    cal = str(C.REPO / "tests" / "calibration")
    jobs: dict[str, list[str]] = {
        "all|full": [cal, "-m", "calibration_full"],
        "all|fast": [cal, "-m", "calibration and not calibration_full"],
    }
    for suite in C.array_suites():
        for full in (True, False):
            jobs[f"{suite}|{'full' if full else 'fast'}"] = C.pytest_args(suite, full)
    d = tmp_path_factory.mktemp("collect")
    script, out = d / "collect.py", d / "out.json"
    script.write_text(_COLLECTOR)
    proc = subprocess.run(
        [sys.executable, str(script), json.dumps(jobs), str(out)],
        cwd=C.REPO,
        capture_output=True,
        text=True,
        timeout=600,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    result: dict[str, dict[str, Any]] = json.loads(out.read_text())
    return result


def _check_one_owner(selections: dict[str, dict[str, Any]], mode: str) -> None:
    everything = set(selections[f"all|{mode}"]["ids"])
    assert everything, f"no calibration tests collected ({mode})"
    owners: dict[str, list[str]] = {}
    for suite in C.array_suites():
        for nodeid in selections[f"{suite}|{mode}"]["ids"]:
            owners.setdefault(nodeid, []).append(suite)
    twice = {n: s for n, s in owners.items() if len(s) > 1}
    assert not twice, (
        f"{mode}: {len(twice)} tests in several array tasks, e.g. {sorted(twice.items())[:5]}"
    )
    orphans = sorted(everything - set(owners))
    assert not orphans, f"{mode}: {len(orphans)} tests in no array task, e.g. {orphans[:5]}"
    assert set(owners) <= everything


def test_every_calibration_test_runs_in_exactly_one_suite(
    selections: dict[str, dict[str, Any]],
) -> None:
    """Both modes and the harness attack file in ONE test: under xdist each test using the
    module fixture may land on another worker and rerun the collector subprocess (backlog)."""
    for mode in ("full", "fast"):
        _check_one_owner(selections, mode)
    # attacks/test_attack_calibration_harness.py names both gene_fdr and regions.
    hits = {
        suite
        for suite in C.array_suites()
        for nodeid in selections[f"{suite}|full"]["ids"]
        if "test_fdr_check_passes_exact_level_bh" in nodeid
    }
    assert len(hits) == 1, hits


# ------------------------------------------- 5. sharded runs honour ADR-0015 (lead decision 5)
#
# ``--merge`` evaluates block 1 from the shard records and, for every failing check, runs block 2
# itself (replicates R..2R-1, sequentially in the merge task) before the verdict; it never reports
# ``fail`` without block 2. The merge task selects a suite's tests through ``pytest_args`` (the one
# hook that maps a suite to pytest arguments), so these tests point it at a tiny probe test file
# whose replicate function is cheap and deterministic: replicate i returns p = 0.0 (rejected) if
# i is in block 1 (or in mode "both"), else 0.5. The probe recognises block 1 by its uniform draw.

STREAM_M = "regions.review_merge_block2"
R_M = 20
PROBE_MODE_ENV = "GIPPA_TEST_MERGE_PROBE_MODE"  # "b1": only block 1 fails; "both": both fail
PROBE_DRAWS_ENV = "GIPPA_TEST_MERGE_PROBE_DRAWS"  # file the probe appends its draws to

_PROBE = textwrap.dedent(
    f"""
    import os

    import _criteria as C

    STREAM, R = {STREAM_M!r}, {R_M}
    B1 = {{float(C.replicate_rng(STREAM, i).uniform()) for i in range(R)}}


    def fn(rng):
        u = float(rng.uniform())
        with open(os.environ[{PROBE_DRAWS_ENV!r}], "a") as fh:
            fh.write(repr(u) + "\\n")
        bad = os.environ[{PROBE_MODE_ENV!r}] == "both" or u in B1
        return 0.0 if bad else 0.5


    def test_probe():
        C.calibrate("regions", "level", 0.05, fn, R, stream=STREAM)
    """
)


def _probe_value(mode: str) -> Callable[[int, float], float]:
    b1 = set(_expected_draws(STREAM_M, range(R_M)))
    return lambda _i, u: 0.0 if (mode == "both" or u in b1) else 0.5


def _block1_shards(log_dir: Path, mode: str, monkeypatch: pytest.MonkeyPatch) -> None:
    """Shard records of block 1 only (0:7, 7:20), as the job array writes them."""
    monkeypatch.delenv(C.BLOCK_ENV, raising=False)
    for b in ["0:7", "7:20"]:
        monkeypatch.setenv(C.SHARD_ENV, b)
        C.calibrate(
            "regions", "level", 0.05, Recorder(_probe_value(mode)), R_M,
            stream=STREAM_M, log_dir=log_dir,
        )  # fmt: skip
    monkeypatch.delenv(C.SHARD_ENV, raising=False)


def _setup_probe(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str, body: str = _PROBE
) -> tuple[Path, Path]:
    log_dir, draws = tmp_path / "logs", tmp_path / "draws.txt"
    probe_dir = tmp_path / "probe"
    probe_dir.mkdir()
    probe = probe_dir / "test_merge_probe.py"
    probe.write_text(body)
    draws.write_text("")
    _block1_shards(log_dir, mode, monkeypatch)
    cal = str(Path(C.__file__).resolve().parent)
    monkeypatch.setenv("PYTHONPATH", f"{cal}:{os.environ.get('PYTHONPATH', '')}")
    monkeypatch.setenv(C.LOG_DIR_ENV, str(log_dir))
    monkeypatch.setenv(PROBE_MODE_ENV, mode)
    monkeypatch.setenv(PROBE_DRAWS_ENV, str(draws))
    monkeypatch.setattr(
        C, "pytest_args", lambda _suite, _full: [str(probe), "-q", "-p", "no:cacheprovider"]
    )
    return log_dir, draws


def _merged(log_dir: Path) -> list[dict[str, Any]]:
    path = log_dir / "regions.merged.jsonl"
    return C.read_results(path) if path.exists() else []


def _reference(tmp_path: Path, mode: str) -> dict[str, Any]:
    """Unsharded ``calibrate`` of the same check: block 1 fails, so block 2 runs."""
    return C.calibrate(
        "regions", "level", 0.05, Recorder(_probe_value(mode)), R_M,
        stream=STREAM_M, log_dir=tmp_path / "reference",
    )  # fmt: skip


@pytest.mark.parametrize(
    ("mode", "status", "code"), [("b1", "unconfirmed", 0), ("both", "fail", 1)]
)
def test_merge_runs_block_two_itself(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str, status: str, code: int
) -> None:
    log_dir, draws = _setup_probe(tmp_path, monkeypatch, mode)
    assert C.main(["--suite", "regions", "--merge"]) == code
    got = [float(x) for x in draws.read_text().split()]
    assert got == _expected_draws(STREAM_M, range(R_M, 2 * R_M)), "block 2 = replicates R..2R-1"
    merged = [r for r in _merged(log_dir) if r.get("stream") == STREAM_M]
    assert len(merged) == 1, merged
    m = merged[0]
    assert m["status"] == status and m["block"] == 2 and m["R"] == R_M
    assert m["ok"] is (mode == "b1")
    ref = _reference(tmp_path, mode)
    assert ref["status"] == status and ref["block"] == 2
    assert m["seed"] == ref["seed"] == C.suite_seed(STREAM_M, R_M)
    assert m["observed"] == pytest.approx(ref["observed"])
    assert m["bound"] == pytest.approx(ref["bound"])


def test_merge_does_not_rerun_passing_checks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Block 1 passes: no block 2 is run and the merged record is block 1, ``pass``."""
    log_dir, draws = tmp_path / "logs", tmp_path / "draws.txt"
    monkeypatch.delenv(C.BLOCK_ENV, raising=False)
    for b in ["0:7", "7:20"]:
        monkeypatch.setenv(C.SHARD_ENV, b)
        C.calibrate(
            "regions", "level", 0.05, Recorder(lambda _i, _u: 0.5), R_M,
            stream=STREAM_M, log_dir=log_dir,
        )  # fmt: skip
    monkeypatch.delenv(C.SHARD_ENV, raising=False)
    monkeypatch.setenv(C.LOG_DIR_ENV, str(log_dir))

    def no_pytest(*_a: object, **_k: object) -> object:
        raise AssertionError("a passing block 1 needs no block 2")

    monkeypatch.setattr(C.subprocess, "run", no_pytest)
    monkeypatch.setattr(C, "pytest_args", no_pytest)
    assert C.main(["--suite", "regions", "--merge"]) == 0
    assert not draws.exists()
    (m,) = _merged(log_dir)
    assert m["status"] == "pass" and m["block"] == 1


def test_merge_never_fails_without_block_two(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """If block 2 cannot be run (its tests log nothing), --merge errors instead of saying fail."""
    empty = "def test_probe():\n    pass\n"
    log_dir, _draws = _setup_probe(tmp_path, monkeypatch, "both", body=empty)
    code = C.main(["--suite", "regions", "--merge"])
    assert code not in (0, 1), "exit 1 means a confirmed fail; 0 means pass"
    assert all(r.get("status") != "fail" for r in _merged(log_dir))
    assert all(r.get("status") != "pass" for r in _merged(log_dir))


def test_merge_shards_never_reports_fail_from_block_one_alone(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The record-only merge cannot run block 2: it must not turn block 1 alone into ``fail``."""
    _block1_shards(tmp_path, "both", monkeypatch)
    try:
        merged = C.merge_shards("regions", log_dir=tmp_path)
    except (ValueError, LookupError):
        return  # refusing to decide is fail closed
    (m,) = merged
    assert m["status"] not in ("fail", "pass", "unconfirmed"), m["status"]
