"""``scripts/slurm/calibration_sharded.sh SUITE N [extra]``: N shard jobs plus one dependent
``--merge`` job of ``calibration_full.sbatch`` sharing one ``GIAP_CALIBRATION_LOG_DIR``.

Runs under bash with a stub ``sbatch`` on PATH (no job is submitted): the stub logs its arguments
(one call per line, fields separated by ``\\x1f``) and prints an increasing job id.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests" / "calibration"))
import _criteria as C  # noqa: E402

SCRIPT = REPO / "scripts" / "slurm" / "calibration_sharded.sh"

STUB_SBATCH = r"""#!/bin/bash
# stub sbatch: log the call; print id (1000 + call number); STUB_FAIL_AT=<call number> fails it
n=$(( $(cat "$STUB_COUNT" 2>/dev/null || echo 0) + 1 )); echo "$n" > "$STUB_COUNT"
if [ "${STUB_FAIL_AT:-0}" = "$n" ]; then echo "sbatch: error" >&2; exit 1; fi
(IFS=$'\x1f'; echo "$*") >> "$STUB_LOG"
echo "$((1000 + n));cluster"
"""


def run(tmp: Path, *args: str, env: dict[str, str] | None = None) -> tuple[int, list[list[str]], str]:
    stub_repo = tmp / "giap"
    (stub_repo / "scripts").mkdir(parents=True, exist_ok=True)
    (stub_repo / "scripts" / "farmshare_env.sh").write_text("# stub\n")
    bindir = tmp / "bin"
    bindir.mkdir(exist_ok=True)
    (bindir / "sbatch").write_text(STUB_SBATCH)
    (bindir / "sbatch").chmod(0o755)
    log = tmp / "calls.txt"
    for f in (log, tmp / "count"):
        f.unlink(missing_ok=True)
    full = {
        "PATH": f"{bindir}{os.pathsep}{os.environ['PATH']}",
        "HOME": str(tmp),
        "GIAP_REPO": str(stub_repo),
        "GIAP_ROOT": str(tmp / "root"),
        "GIAP_PYTHON": sys.executable,
        "STUB_LOG": str(log),
        "STUB_COUNT": str(tmp / "count"),
        **(env or {}),
    }
    bash = shutil.which("bash") or "/bin/bash"
    proc = subprocess.run([bash, str(SCRIPT), *args], env=full, capture_output=True, text=True,
                          timeout=60, check=False)  # fmt: skip
    calls = [line.split("\x1f") for line in log.read_text().splitlines()] if log.exists() else []
    return proc.returncode, calls, proc.stdout + proc.stderr


def opt(call: list[str], name: str) -> str | None:
    hits = [a.split("=", 1)[1] for a in call if a.startswith(f"--{name}=")]
    return hits[-1] if hits else None


def harness_args(call: list[str]) -> list[str]:
    i = next(i for i, a in enumerate(call) if a.endswith("calibration_full.sbatch"))
    return call[i + 1 :]


@pytest.mark.parametrize("n", [1, 3, 10])
def test_submits_n_shards_then_a_dependent_merge(tmp_path: Path, n: int) -> None:
    logdir = tmp_path / "shared"
    rc, calls, out = run(tmp_path, "discordance", str(n), "--fast",
                         env={"GIAP_CALIBRATION_LOG_DIR": str(logdir)})  # fmt: skip
    assert rc == 0, out
    assert len(calls) == n + 1, calls
    task = str(C.array_suites().index("discordance"))
    ids = [str(1001 + k) for k in range(n)]
    for k, call in enumerate(calls[:n]):
        assert opt(call, "array") == task and opt(call, "dependency") is None, call
        assert harness_args(call) == ["--shard", f"{k}/{n}", "--fast"], call
        assert opt(call, "job-name") == f"giap-cal-discordance-s{k}"
    merge = calls[-1]
    assert opt(merge, "dependency") == "afterany:" + ":".join(ids), merge
    assert harness_args(merge) == ["--merge", "--fast"] and opt(merge, "array") == task
    assert str(opt(merge, "job-name")).startswith("giap-")
    for call in calls:  # one shared log dir, and the job runs this checkout
        exp = str(opt(call, "export"))
        assert f"GIAP_CALIBRATION_LOG_DIR={logdir}" in exp.split(","), exp
        assert f"GIAP_WORKTREE={REPO}" in exp.split(","), exp
        assert call[-len(harness_args(call)) - 1] == str(REPO / "scripts/slurm/calibration_full.sbatch")
    assert f"shards={','.join(ids)} merge={1001 + n} log={logdir}" in out


def test_shards_cover_every_replicate_of_the_suite(tmp_path: Path) -> None:
    rc, calls, out = run(tmp_path, "regions", "7")
    assert rc == 0, out
    R = C.replicates("regions", True)
    got: list[int] = []
    for call in calls[:-1]:
        k, n = (int(x) for x in harness_args(call)[1].split("/"))
        got.extend(C.shard_replicates(R, n, k))
    assert got == list(range(R))


def test_default_log_dir_is_fresh_under_giap_root(tmp_path: Path) -> None:
    rc, calls, out = run(tmp_path, "discordance", "2")
    assert rc == 0, out
    dirs = {re.search(r"GIAP_CALIBRATION_LOG_DIR=([^,]+)", str(opt(c, "export"))).group(1)  # type: ignore[union-attr]
            for c in calls}  # fmt: skip
    (d,) = dirs
    assert d.startswith(str(tmp_path / "root" / "logs" / "calibration" / "sharded-discordance-"))


def test_sbatch_extra_options_reach_every_job(tmp_path: Path) -> None:
    rc, calls, out = run(tmp_path, "discordance", "2", env={"GIAP_SBATCH_EXTRA": "--time=01:00:00 --qos=normal"})
    assert rc == 0, out
    assert all(opt(c, "time") == "01:00:00" and opt(c, "qos") == "normal" for c in calls), calls


@pytest.mark.parametrize(
    "args",
    [("discordance",), ("discordance", "0"), ("discordance", "11"), ("discordance", "x"),
     ("no_such_suite", "3"), ("discordance", "3", "--shard", "1/3"),
     ("discordance", "3", "--merge"), ("discordance", "3", "--dry-run"),
     ("discordance", "3", "--task-id=4")],
)  # fmt: skip
def test_bad_arguments_submit_nothing(tmp_path: Path, args: tuple[str, ...]) -> None:
    rc, calls, out = run(tmp_path, *args)
    assert rc == 2 and calls == [], (rc, calls, out)


def test_existing_suite_log_is_refused(tmp_path: Path) -> None:
    logdir = tmp_path / "old"
    logdir.mkdir()
    (logdir / "discordance.jsonl").write_text("{}\n")
    rc, calls, out = run(tmp_path, "discordance", "3", env={"GIAP_CALIBRATION_LOG_DIR": str(logdir)})
    assert rc == 2 and calls == [], out


@pytest.mark.parametrize("fail_at", [1, 2, 4])
def test_failed_submission_stops_and_never_merges_a_partial_set(tmp_path: Path, fail_at: int) -> None:
    """3 shards: a failed shard submission (call 1 or 2) submits no merge; a failed merge
    submission (call 4) exits non-zero and names the shards left without one."""
    rc, calls, out = run(tmp_path, "discordance", "3", env={"STUB_FAIL_AT": str(fail_at)})
    assert rc == 1, out
    assert not any("--merge" in c for c in calls), calls
    assert len(calls) == fail_at - 1, calls
