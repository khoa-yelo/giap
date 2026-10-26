"""``scripts/slurm/nightly.sbatch`` robustness fixes from the 2026-10-23/26 nightly failures.

Run for real in the M2-21 sandbox (``tests/helpers/slurm_sandbox.py``), with this module's own
``python``/``git``/``nf-test`` wrappers in the sandbox ``bin/``. Pinned:

1. Slow-node guard: at task start the script times ``python -c '…import giap…'``; when that takes
   longer than ``NIGHTLY_SLOW_NODE_SECS`` (default 60) seconds the task resubmits itself
   (``sbatch --array=<task> --exclude=<previous exclusions>,<this host>``, the same commit and
   worktree, ``NIGHTLY_SLOW_HOPS`` + 1, ``NIGHTLY_EXCLUDE`` = the new list), extends the pending
   summary job and exits 0 without writing its own nightly.tsv row. Bounded by
   ``NIGHTLY_SLOW_MAX_HOPS`` (default 2); ``NIGHTLY_SLOW_NODE_GUARD=0`` disables it.
2. Rows on kill: ``#SBATCH --signal=B:TERM@600``, and the long pipelines run in the background
   (``cmd & wait $!``) so a SIGTERM sent to the batch shell alone runs the TERM trap at once and
   the task's row (exit 143) is written well before Slurm's KillWait.
3. Stale git lock: an ``index.lock`` left in the task worktree's git dir is removed before the
   checkout; a worktree whose checkout still fails is removed (``git worktree remove --force``),
   pruned and re-added; ``checkout -f`` is skipped when HEAD already is the target commit.
4. ``tests/conftest.py``'s capability probe timeout comes from ``GIAP_PROBE_TIMEOUT`` (default
   30 s), and nightly.sbatch exports ``GIAP_PROBE_TIMEOUT`` (default ≥ 180 s, an outer value wins).
"""

from __future__ import annotations

import contextlib
import importlib.util
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path
from types import ModuleType

import pytest
from helpers.gnu_userland import require_gnu_userland
from helpers.slurm_sandbox import (
    _PYTHON,
    REPO,
    Sandbox,
    assert_honours_giap_repo,
    export_pairs,
    make_sandbox,
    option,
    write_exe,
)

NIGHTLY = REPO / "scripts" / "slurm" / "nightly.sbatch"
CONFTEST = REPO / "tests" / "conftest.py"
JOB = "4242"
HOST = "rye-03"

# Inserted after the shebang of the sandbox python wrapper:
# * `python -c '…import giap…'` (the slow-node probe) sleeps $FAKE_IMPORT_SLEEP s, logs, exits 0;
# * `python -m pytest` logs GIAP_PROBE_TIMEOUT; FAKE_PYTEST=sleep marks its start and sleeps 60 s.
_HOOK = r"""
if [ "${1:-}" = "-c" ] && printf '%s' "${2:-}" | grep -Eq 'import giap([^a-zA-Z0-9_]|$)'; then
  printf 'import-probe\n' >> "$FAKE_SLURM_DIR/calls.log"
  sleep "${FAKE_IMPORT_SLEEP:-0}"; exit 0
fi
if [ "${1:-}" = "-m" ] && [ "${2:-}" = "pytest" ]; then
  printf 'pytest GIAP_PROBE_TIMEOUT=%s\n' "${GIAP_PROBE_TIMEOUT:-<unset>}" \
    >> "$FAKE_SLURM_DIR/calls.log"
  if [ "${FAKE_PYTEST:-}" = sleep ]; then
    : > "$FAKE_SLURM_DIR/long_started"; sleep 60; echo "1 passed in 60.00s"; exit 0
  fi
fi
"""

_GIT = r"""#!/bin/bash
printf '%s\n' "$*" >> "$FAKE_SLURM_DIR/git.log"
exec {git} "$@"
"""

_NF_TEST = r"""#!/bin/bash
printf '%s\n' "nf-test $*" >> "$FAKE_SLURM_DIR/calls.log"
if [ "${FAKE_NF_TEST:-}" = sleep ]; then : > "$FAKE_SLURM_DIR/long_started"; sleep 60; fi
echo "SUCCESS: Executed 1 tests"
"""


@pytest.fixture
def sb(tmp_path: Path) -> Sandbox:
    require_gnu_userland()
    assert_honours_giap_repo(NIGHTLY.read_text(), "nightly.sbatch")
    box = make_sandbox(tmp_path)
    head, _, rest = _PYTHON.replace("{python}", sys.executable).partition("\n")
    for name in ("python", "python3"):
        write_exe(box.bin / name, f"{head}\n{_HOOK}{rest}")
    git = shutil.which("git")
    assert git, "git is not on PATH"
    write_exe(box.bin / "git", _GIT.replace("{git}", git))
    write_exe(box.bin / "nf-test", _NF_TEST)
    return box


def _env(task: int, job: str = JOB, **extra: str) -> dict[str, str]:
    return {
        "SLURM_ARRAY_TASK_ID": str(task),
        "SLURM_ARRAY_JOB_ID": job,
        "SLURM_JOB_ID": str(int(job) + 1000),
        "SLURM_ARRAY_TASK_COUNT": "4",
        "SLURM_CPUS_PER_TASK": "2",
        "SLURM_JOB_NAME": "giap-nightly",
        "NIGHTLY_POST": "dry",
        "FAKE_HOSTNAME": HOST,
        **extra,
    }


def _nightly(sb: Sandbox, task: int, job: str = JOB, **extra: str):  # type: ignore[no-untyped-def]
    return sb.run(["bash", str(sb.repo / "scripts" / "slurm" / "nightly.sbatch")],
                  env=_env(task, job, **extra))  # fmt: skip


def _rows(sb: Sandbox, job_task: str) -> list[list[str]]:
    return [r for r in sb.tsv_rows() if len(r) >= 5 and r[1] == job_task]


def _self_resubmissions(sb: Sandbox) -> list[list[str]]:
    return [c["argv"] for c in sb.sbatch_calls()
            if c["argv"] and c["argv"][-1].endswith("nightly.sbatch")]  # fmt: skip


def _sbatch_excludes(text: str) -> list[str]:
    out: list[str] = []
    for m in re.finditer(r"^#SBATCH\s+(?:--exclude[= ]|-x\s*)(\S+)", text, re.M):
        out += [h for h in m[1].split(",") if h]
    return out


def _diag(res: subprocess.CompletedProcess[str], sb: Sandbox) -> str:
    return f"\nstdout:\n{res.stdout}\nstderr:\n{res.stderr}\nrows: {sb.tsv_rows()}"


# ---------------------------------------------------------------- 1. slow-node guard


@pytest.mark.parametrize("task", [0, 2])
def test_slow_node_resubmits_with_exclude_and_writes_no_row(sb: Sandbox, task: int) -> None:
    sb.set_next_job_id(7001)
    sb.set_squeue("giap-nightly-summary", [f"9100 nightly:{JOB}"])
    res = _nightly(sb, task, NIGHTLY_SLOW_NODE_SECS="1", FAKE_IMPORT_SLEEP="3",
                   NIGHTLY_EXCLUDE="oat-07")  # fmt: skip
    diag = _diag(res, sb)
    assert "import-probe" in sb.calls_log(), (
        "the slow-node probe (python -c 'import giap') never ran" + diag
    )
    subs = _self_resubmissions(sb)
    assert len(subs) == 1, f"expected one resubmission of nightly.sbatch, got {subs}" + diag
    argv = subs[0]
    subs_env = [c["env"] for c in sb.sbatch_calls()
                if c["argv"] and c["argv"][-1].endswith("nightly.sbatch")]  # fmt: skip
    excl = (option(argv, "exclude", "x") or "").split(",")
    assert HOST in excl, f"the slow node {HOST} is not excluded: {argv}"
    assert "oat-07" in excl, f"the earlier exclusion (NIGHTLY_EXCLUDE) is dropped: {argv}"
    for h in _sbatch_excludes(NIGHTLY.read_text()):  # e.g. the rye-02 stopgap (#402)
        assert h in excl, (
            f"the #SBATCH --exclude host {h} is dropped by --exclude on the command line"
        )
    assert option(argv, "array") == str(task), argv
    ex = export_pairs(argv)
    assert ex.get("NIGHTLY_ARRAY_JOB") == JOB, ex
    assert ex.get("NIGHTLY_EXPECT") == "4", ex
    assert ex.get("NIGHTLY_SLOW_HOPS") == "1", ex
    # a comma list cannot travel inside --export (sbatch splits it at commas): it rides in the
    # environment, which --export=ALL forwards (lead fix after the implementer's report, #403)
    env_excl = subs_env[0].get("NIGHTLY_EXCLUDE", "")
    assert set(env_excl.split(",")) >= {HOST, "oat-07"}, (env_excl, ex)
    ref = ex.get("NIGHTLY_REF", "")
    assert ref and sb.head.startswith(ref), f"the resubmission does not pin the commit: {ex}"
    assert res.returncode == 0, diag
    assert not _rows(sb, f"{JOB}_{task}"), "the slow-node hop wrote a row for the task" + diag
    assert not (sb.state / "pytest.argv").exists(), "the tests still ran on the slow node"
    assert any(re.search(rf"scontrol .*afterany:{JOB}:7001\b", ln) for ln in sb.calls_log()), (
        f"the pending summary job is not extended to the resubmitted job\n{sb.calls_log()}"
    )


def test_slow_node_hops_are_counted(sb: Sandbox) -> None:
    res = _nightly(sb, 2, NIGHTLY_SLOW_NODE_SECS="1", FAKE_IMPORT_SLEEP="3", NIGHTLY_SLOW_HOPS="1",
                   NIGHTLY_EXCLUDE="rye-01")  # fmt: skip
    subs = _self_resubmissions(sb)
    assert len(subs) == 1, f"one more hop is allowed (default max 2): {subs}" + _diag(res, sb)
    ex = export_pairs(subs[0])
    assert ex.get("NIGHTLY_SLOW_HOPS") == "2", ex
    assert set((option(subs[0], "exclude", "x") or "").split(",")) >= {"rye-01", HOST}, subs[0]


@pytest.mark.parametrize(
    ("hops", "extra"),
    [("2", {}), ("1", {"NIGHTLY_SLOW_MAX_HOPS": "1"}), ("0", {"NIGHTLY_SLOW_NODE_GUARD": "0"})],
    ids=["default-max-reached", "max-overridden", "guard-disabled"],
)
def test_slow_node_guard_is_bounded_and_can_be_disabled(
    sb: Sandbox, hops: str, extra: dict[str, str]
) -> None:
    res = _nightly(sb, 2, NIGHTLY_SLOW_NODE_SECS="1", FAKE_IMPORT_SLEEP="3",
                   NIGHTLY_SLOW_HOPS=hops, **extra)  # fmt: skip
    assert not _self_resubmissions(sb), "resubmitted beyond the hop limit / with the guard off"
    rows = _rows(sb, f"{JOB}_2")
    assert len(rows) == 1 and rows[0][3] == "0", "the task did not run and report" + _diag(res, sb)
    assert (sb.state / "pytest.argv").exists()


def test_fast_node_is_not_resubmitted(sb: Sandbox) -> None:
    res = _nightly(sb, 2, FAKE_IMPORT_SLEEP="0")
    assert not _self_resubmissions(sb)
    rows = _rows(sb, f"{JOB}_2")
    assert len(rows) == 1 and rows[0][3] == "0", _diag(res, sb)


def test_failed_slow_node_resubmission_still_writes_a_row(sb: Sandbox) -> None:
    res = _nightly(sb, 2, NIGHTLY_SLOW_NODE_SECS="1", FAKE_IMPORT_SLEEP="3", FAKE_SBATCH_FAIL="all")
    assert "import-probe" in sb.calls_log(), "the slow-node probe never ran" + _diag(res, sb)
    assert len(_rows(sb, f"{JOB}_2")) == 1, "no row after a refused resubmission" + _diag(res, sb)


def test_slow_node_default_threshold_is_60_seconds() -> None:
    text = NIGHTLY.read_text()
    m = re.search(r"\$\{NIGHTLY_SLOW_NODE_SECS:-([0-9]+)\}", text)
    assert m, "NIGHTLY_SLOW_NODE_SECS (default threshold) is not read by nightly.sbatch"
    assert 30 <= int(m[1]) <= 120, m[0]


# ---------------------------------------------------------------- 2. rows on SIGTERM


def test_sbatch_asks_for_an_early_term_to_the_batch_shell() -> None:
    assert re.search(r"^#SBATCH\s+--signal=B:TERM@600\s*$", NIGHTLY.read_text(), re.M), (
        "nightly.sbatch lacks `#SBATCH --signal=B:TERM@600`"
    )


@pytest.mark.parametrize(
    ("task", "fake"), [(2, "FAKE_PYTEST"), (0, "FAKE_NF_TEST")], ids=["pytest", "nf-test"]
)
def test_term_to_batch_shell_writes_row_promptly(sb: Sandbox, task: int, fake: str) -> None:
    """SIGTERM to the batch shell only (``--signal=B:``), mid-pipeline: the trap writes the row."""
    require_gnu_userland()
    env = {**sb.env, **_env(task, **{fake: "sleep"})}
    out = (sb.tmp / "term.out").open("w")
    proc = subprocess.Popen(
        ["bash", str(sb.repo / "scripts" / "slurm" / "nightly.sbatch")],
        env=env, cwd=str(sb.tmp), stdout=out, stderr=subprocess.STDOUT,
        start_new_session=True,
    )  # fmt: skip
    try:
        deadline = time.monotonic() + 60
        while not (sb.state / "long_started").exists():
            assert proc.poll() is None, (sb.tmp / "term.out").read_text()
            assert time.monotonic() < deadline, "the long-running step never started"
            time.sleep(0.1)
        time.sleep(0.5)
        os.kill(proc.pid, signal.SIGTERM)
        t0 = time.monotonic()
        try:
            rc = proc.wait(timeout=20)
        except subprocess.TimeoutExpired:
            pytest.fail("the batch shell did not act on SIGTERM while the pipeline ran "
                        "(run it as `cmd & wait $!`)")  # fmt: skip
        took = time.monotonic() - t0
    finally:
        with contextlib.suppress(ProcessLookupError):
            os.killpg(proc.pid, signal.SIGKILL)
        out.close()
    log = (sb.tmp / "term.out").read_text()
    assert took < 20, took
    assert rc == 143, f"exit {rc}\n{log}"
    rows = _rows(sb, f"{JOB}_{task}")
    assert len(rows) == 1, f"expected one row for the killed task, got {sb.tsv_rows()}\n{log}"
    assert rows[0][3] == "143", rows[0]
    assert re.search(r"SIGTERM|TIMEOUT|terminated", rows[0][4], re.I), rows[0]


# ---------------------------------------------------------------- 3. stale git lock


def _worktree(sb: Sandbox, task: int = 2) -> Path:
    return sb.root / "ci" / f"main-t{task}"


def _gitdir(w: Path) -> Path:
    out = subprocess.run(["git", "-C", str(w), "rev-parse", "--absolute-git-dir"],
                         capture_output=True, text=True, check=True).stdout.strip()  # fmt: skip
    return Path(out)


def _advance_origin(sb: Sandbox) -> str:
    env = {**os.environ, **sb.env}
    run = {"cwd": str(sb.repo), "env": env, "check": True, "capture_output": True}
    subprocess.run(["git", "commit", "-q", "--allow-empty", "-m", "next"], **run)  # type: ignore[call-overload]
    subprocess.run(["git", "push", "-q", "origin", "HEAD:main"], **run)  # type: ignore[call-overload]
    return sb.head


def _head(w: Path) -> str:
    return subprocess.run(["git", "-C", str(w), "rev-parse", "HEAD"],
                          capture_output=True, text=True, check=True).stdout.strip()  # fmt: skip


def _first_run(sb: Sandbox) -> Path:
    res = _nightly(sb, 2, job="4000")
    w = _worktree(sb)
    assert w.is_dir() and _rows(sb, "4000_2"), "the first (setup) run failed" + _diag(res, sb)
    (sb.state / "git.log").unlink(missing_ok=True)
    return w


def test_stale_index_lock_is_removed_before_checkout(sb: Sandbox) -> None:
    w = _first_run(sb)
    lock = _gitdir(w) / "index.lock"
    lock.write_text("")
    target = _advance_origin(sb)
    res = _nightly(sb, 2)
    diag = _diag(res, sb)
    assert not lock.exists(), "the stale index.lock is still there" + diag
    assert _head(w) == target, "the worktree is not at the new origin/main" + diag
    rows = _rows(sb, f"{JOB}_2")
    assert len(rows) == 1 and rows[0][3] == "0", diag


def test_checkout_is_skipped_when_head_is_already_the_target(sb: Sandbox) -> None:
    w = _first_run(sb)
    (_gitdir(w) / "index.lock").write_text("")
    res = _nightly(sb, 2)
    diag = _diag(res, sb)
    git_log = (sb.state / "git.log").read_text().splitlines()
    assert not [ln for ln in git_log if re.search(r"\bcheckout\b", ln)], (
        f"checkout ran although HEAD already is the target\n{git_log}" + diag
    )
    rows = _rows(sb, f"{JOB}_2")
    assert len(rows) == 1 and rows[0][3] == "0", diag


def test_broken_worktree_is_removed_pruned_and_re_added(sb: Sandbox) -> None:
    w = _first_run(sb)
    (_gitdir(w) / "index").write_text(
        "garbage\n"
    )  # checkout -f fails: "index file smaller than expected"
    target = _advance_origin(sb)
    res = _nightly(sb, 2)
    diag = _diag(res, sb)
    git_log = (sb.state / "git.log").read_text().splitlines()
    remove = [
        ln for ln in git_log if re.search(r"\bworktree remove\b", ln) and "--force" in ln.split()
    ]
    prune = [ln for ln in git_log if re.search(r"\bworktree prune\b", ln)]
    assert remove, f"no `git worktree remove --force`\n{git_log}" + diag
    assert prune, f"no `git worktree prune`\n{git_log}" + diag
    assert _head(w) == target, "the re-added worktree is not at the new origin/main" + diag
    st = subprocess.run(
        ["git", "-C", str(w), "status", "--porcelain"], capture_output=True, text=True
    )
    assert st.returncode == 0, st.stderr
    rows = _rows(sb, f"{JOB}_2")
    assert len(rows) == 1 and rows[0][3] == "0", diag


# ---------------------------------------------------------------- 4. probe timeout


def _load_conftest(monkeypatch: pytest.MonkeyPatch, value: str | None) -> ModuleType:
    if value is None:
        monkeypatch.delenv("GIAP_PROBE_TIMEOUT", raising=False)
    else:
        monkeypatch.setenv("GIAP_PROBE_TIMEOUT", value)
    spec = importlib.util.spec_from_file_location("giap_conftest_probe_copy", CONFTEST)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_probe_timeout_honours_env(monkeypatch: pytest.MonkeyPatch) -> None:
    mod = _load_conftest(monkeypatch, "0.5")
    t0 = time.monotonic()
    ok = mod._ok(["sleep", "3"])
    took = time.monotonic() - t0
    assert not ok, "a probe slower than GIAP_PROBE_TIMEOUT=0.5 still counts as available"
    assert took < 2.5, f"the probe ran {took:.1f}s despite GIAP_PROBE_TIMEOUT=0.5"


@pytest.mark.parametrize("value", [None, "", "not-a-number", "120"])
def test_probe_timeout_default_and_bad_values(
    monkeypatch: pytest.MonkeyPatch, value: str | None
) -> None:
    mod = _load_conftest(monkeypatch, value)
    assert mod._ok(["sleep", "0.1"]), f"a fast probe fails with GIAP_PROBE_TIMEOUT={value!r}"


def test_probe_timeout_default_is_30_seconds() -> None:
    text = CONFTEST.read_text()
    assert "GIAP_PROBE_TIMEOUT" in text, "tests/conftest.py does not read GIAP_PROBE_TIMEOUT"
    assert re.search(r"GIAP_PROBE_TIMEOUT[^\n]*\b30\b", text), "the default is no longer 30 s"


def _pytest_probe_timeouts(sb: Sandbox) -> list[str]:
    return [
        m[1] for ln in sb.calls_log() if (m := re.fullmatch(r"pytest GIAP_PROBE_TIMEOUT=(.*)", ln))
    ]


@pytest.mark.parametrize("task", [0, 2])
def test_nightly_exports_a_long_probe_timeout(sb: Sandbox, task: int) -> None:
    res = _nightly(sb, task)
    seen = _pytest_probe_timeouts(sb)
    assert seen, "pytest never ran" + _diag(res, sb)
    for v in seen:
        assert re.fullmatch(r"[0-9]+(\.[0-9]+)?", v), f"GIAP_PROBE_TIMEOUT={v} for pytest"
        assert float(v) >= 180, f"GIAP_PROBE_TIMEOUT={v} < 180 s"


def test_nightly_keeps_an_outer_probe_timeout(sb: Sandbox) -> None:
    res = _nightly(sb, 2, GIAP_PROBE_TIMEOUT="300")
    assert _pytest_probe_timeouts(sb) == ["300"], _diag(res, sb)
