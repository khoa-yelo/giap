"""``scripts/slurm/calibration_sharded.sh SUITE N [extra]``: N shard jobs plus one dependent
``--merge`` job of ``calibration_full.sbatch`` sharing one ``GIPPA_CALIBRATION_LOG_DIR``.

Runs under bash with a stub ``sbatch`` on PATH (no job is submitted): the stub logs its arguments
(one call per line, fields separated by ``\\x1f``), logs the ``GIPPA_CALIBRATION_LOG_DIR`` it was
called with, and prints an increasing job id. A stub ``squeue`` (review of PR #327, R14 guard)
prints the jobs of ``STUB_SQUEUE_ROWS`` (JSON) as real squeue does: ``-o``/``--format`` (``%i
%A %j %C %T %t %u %D %P %M %N %R``, with widths and truncation), ``-O``/``--Format``, ``-h``,
``-u``/``--user``/``--me`` (without one, another user's ``gippa-*`` jobs are listed too),
``-t``/``--states`` and ``-n``/``--name``; ``STUB_SQUEUE_FAIL=1`` makes it exit 1.
"""

from __future__ import annotations

import json
import os
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
printf '%s\n' "${GIPPA_CALIBRATION_LOG_DIR-<unset>}" >> "$STUB_ENV_LOG"
echo "$((1000 + n));cluster"
"""

STUB_SQUEUE = r"""#!PYTHON
import json, os, re, sys
argv = sys.argv[1:]
with open(os.environ["STUB_SQUEUE_LOG"], "a") as fh:
    fh.write(json.dumps(argv) + "\n")
if os.environ.get("STUB_SQUEUE_FAIL") == "1":
    print("squeue: error: fake failure", file=sys.stderr); sys.exit(1)
path = os.environ.get("STUB_SQUEUE_ROWS", "")
rows = json.load(open(path)) if path and os.path.exists(path) else []
fmt = "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R"
big, header, users, states, names = None, True, None, None, None
def val(i, a, short, long):
    if a == short or a == long:
        return argv[i + 1], i + 1
    if a.startswith(long + "="):
        return a.split("=", 1)[1], i
    return a[len(short):], i
i = 0
while i < len(argv):
    a = argv[i]
    if a in ("-h", "--noheader"):
        header = False
    elif a == "--me":
        users = [os.environ["USER"]]
    elif a.startswith(("-u", "--user")):
        v, i = val(i, a, "-u", "--user"); users = v.split(",")
    elif a.startswith(("-t", "--states")):
        v, i = val(i, a, "-t", "--states"); states = [x.upper() for x in v.split(",")]
    elif a.startswith(("-n", "--name")) and not a.startswith("--noheader"):
        v, i = val(i, a, "-n", "--name"); names = v.split(",")
    elif a.startswith(("-o", "--format")):
        v, i = val(i, a, "-o", "--format"); fmt = v
    elif a.startswith(("-O", "--Format")):
        v, i = val(i, a, "-O", "--Format"); big = v
    i += 1
SHORT = {"RUNNING": "R", "PENDING": "PD", "COMPLETING": "CG"}
def keep(r):
    if users is not None and r["user"] not in users:
        return False
    if states is not None and "ALL" not in states and r["state"] not in states \
            and SHORT.get(r["state"]) not in states:
        return False
    return names is None or r["name"] in names
def field(r, c):
    return {"i": r["id"], "A": r["id"], "j": r["name"], "C": r["cpus"], "T": r["state"],
            "t": SHORT.get(r["state"], r["state"]), "u": r["user"], "D": 1, "P": "normal",
            "M": "0:00", "N": "wheat-01", "R": "wheat-01"}.get(c, "")
HEAD = {"i": "JOBID", "A": "JOBID", "j": "NAME", "C": "CPUS", "T": "STATE", "t": "ST",
        "u": "USER", "D": "NODES", "P": "PARTITION", "M": "TIME", "N": "NODELIST",
        "R": "NODELIST(REASON)"}
def render(get):
    def sub(m):
        text = str(get(m.group(3)))
        if m.group(2):
            w = int(m.group(2)); text = text[:w]
            text = text.rjust(w) if m.group(1) else text.ljust(w)
        return text
    return re.sub(r"%(\.)?(\d+)?([A-Za-z])", sub, fmt)
BIG = {"JOBID": "i", "JOBARRAYID": "i", "NAME": "j", "NUMCPUS": "C", "STATE": "T",
       "STATECOMPACT": "t", "USERNAME": "u", "NUMNODES": "D"}
def render_big(get):
    out = []
    for spec in big.split(","):
        name, _, w = spec.partition(":")
        w = int(w.rstrip("|") or 20) if w else 20
        out.append(str(get(BIG.get(name.upper(), "?")))[:w].ljust(w))
    return "".join(out)
line = render_big if big is not None else render
if header:
    print(line(lambda c: HEAD.get(c, "")))
for r in rows:
    if keep(r):
        print(line(lambda c, r=r: field(r, c)))
"""


def run(
    tmp: Path, *args: str, env: dict[str, str] | None = None
) -> tuple[int, list[list[str]], str]:
    stub_repo = tmp / "gippa"
    (stub_repo / "scripts").mkdir(parents=True, exist_ok=True)
    (stub_repo / "scripts" / "farmshare_env.sh").write_text("# stub\n")
    bindir = tmp / "bin"
    bindir.mkdir(exist_ok=True)
    (bindir / "sbatch").write_text(STUB_SBATCH)
    (bindir / "sbatch").chmod(0o755)
    (bindir / "squeue").write_text(STUB_SQUEUE.replace("PYTHON", sys.executable, 1))
    (bindir / "squeue").chmod(0o755)
    log = tmp / "calls.txt"
    for f in (log, tmp / "count", tmp / "env.txt", tmp / "squeue.jsonl"):
        f.unlink(missing_ok=True)
    full = {
        "PATH": f"{bindir}{os.pathsep}{os.environ['PATH']}",
        "HOME": str(tmp),
        "GIPPA_REPO": str(stub_repo),
        "GIPPA_ROOT": str(tmp / "root"),
        "GIPPA_PYTHON": sys.executable,
        "STUB_LOG": str(log),
        "STUB_COUNT": str(tmp / "count"),
        "STUB_ENV_LOG": str(tmp / "env.txt"),
        "STUB_SQUEUE_LOG": str(tmp / "squeue.jsonl"),
        "USER": "tester",
        **(env or {}),
    }
    bash = shutil.which("bash") or "/bin/bash"
    proc = subprocess.run([bash, str(SCRIPT), *args], env=full, capture_output=True, text=True,
                          timeout=60, check=False)  # fmt: skip
    calls = [line.split("\x1f") for line in log.read_text().splitlines()] if log.exists() else []
    return proc.returncode, calls, proc.stdout + proc.stderr


def submit_envs(tmp: Path) -> list[str]:
    """``GIPPA_CALIBRATION_LOG_DIR`` in the environment of each sbatch call."""
    path = tmp / "env.txt"
    return path.read_text().splitlines() if path.exists() else []


def exported_log_dir(call: list[str], submit_env: str) -> str | None:
    """The log dir the job gets, as sbatch parses ``--export`` (comma-separated, ``ALL`` keeps
    the submitting environment): an explicit assignment wins, else the inherited value."""
    entries = str(opt(call, "export")).split(",")
    assigned = [e.split("=", 1)[1] for e in entries if e.startswith("GIPPA_CALIBRATION_LOG_DIR=")]
    if assigned:
        return assigned[-1]
    return submit_env if "ALL" in entries and submit_env != "<unset>" else None


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
                         env={"GIPPA_CALIBRATION_LOG_DIR": str(logdir)})  # fmt: skip
    assert rc == 0, out
    assert len(calls) == n + 1, calls
    task = str(C.array_suites().index("discordance"))
    ids = [str(1001 + k) for k in range(n)]
    for k, call in enumerate(calls[:n]):
        assert opt(call, "array") == task and opt(call, "dependency") is None, call
        assert harness_args(call) == ["--shard", f"{k}/{n}", "--fast"], call
        assert opt(call, "job-name") == f"gippa-cal-discordance-s{k}"
    merge = calls[-1]
    assert opt(merge, "dependency") == "afterany:" + ":".join(ids), merge
    assert harness_args(merge) == ["--merge", "--fast"] and opt(merge, "array") == task
    assert str(opt(merge, "job-name")).startswith("gippa-")
    for call, env in zip(calls, submit_envs(tmp_path), strict=True):  # one shared log dir
        exp = str(opt(call, "export"))
        assert exported_log_dir(call, env) == str(logdir), (exp, env)
        assert f"GIPPA_WORKTREE={REPO}" in exp.split(","), exp
        assert call[-len(harness_args(call)) - 1] == str(
            REPO / "scripts/slurm/calibration_full.sbatch"
        )
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


def test_default_log_dir_is_fresh_under_gippa_root(tmp_path: Path) -> None:
    rc, calls, out = run(tmp_path, "discordance", "2")
    assert rc == 0, out
    dirs = {exported_log_dir(c, e) for c, e in zip(calls, submit_envs(tmp_path), strict=True)}
    (d,) = dirs
    assert d.startswith(str(tmp_path / "root" / "logs" / "calibration" / "sharded-discordance-"))


def test_sbatch_extra_options_reach_every_job(tmp_path: Path) -> None:
    rc, calls, out = run(
        tmp_path, "discordance", "2", env={"GIPPA_SBATCH_EXTRA": "--time=01:00:00 --qos=normal"}
    )
    assert rc == 0, out
    assert all(opt(c, "time") == "01:00:00" and opt(c, "qos") == "normal" for c in calls), calls


@pytest.mark.parametrize(
    "args",
    [("discordance",), ("discordance", "0"), ("discordance", "11"), ("discordance", "x"),
     ("no_such_suite", "3"), ("discordance", "3", "--shard", "1/3"),
     ("discordance", "3", "--merge"), ("discordance", "3", "--dry-run"),
     ("discordance", "3", "--task-id=4"),
     # review of PR #327: a --block 2 run never owns the direct checks and merges to an error
     ("discordance", "3", "--block", "2"), ("discordance", "3", "--block=2"),
     ("discordance", "3", "--fast", "--block", "1"), ("discordance", "3", "--block")],
)  # fmt: skip
def test_bad_arguments_submit_nothing(tmp_path: Path, args: tuple[str, ...]) -> None:
    rc, calls, out = run(tmp_path, *args)
    assert rc == 2 and calls == [], (rc, calls, out)


def test_existing_suite_log_is_refused(tmp_path: Path) -> None:
    logdir = tmp_path / "old"
    logdir.mkdir()
    (logdir / "discordance.jsonl").write_text("{}\n")
    rc, calls, out = run(
        tmp_path, "discordance", "3", env={"GIPPA_CALIBRATION_LOG_DIR": str(logdir)}
    )
    assert rc == 2 and calls == [], out


@pytest.mark.parametrize("fail_at", [1, 2, 4])
def test_failed_submission_stops_and_never_merges_a_partial_set(
    tmp_path: Path, fail_at: int
) -> None:
    """3 shards: a failed shard submission (call 1 or 2) submits no merge; a failed merge
    submission (call 4) exits non-zero and names the shards left without one."""
    rc, calls, out = run(tmp_path, "discordance", "3", env={"STUB_FAIL_AT": str(fail_at)})
    assert rc == 1, out
    assert not any("--merge" in c for c in calls), calls
    assert len(calls) == fail_at - 1, calls


# --------------------------------------------------------------------------- review of PR #327


@pytest.mark.parametrize(
    "stale", ["discordance.jsonl", "discordance.outcomes.jsonl", "discordance.shards.jsonl",
              "discordance.merged.jsonl", "discordance.anything"],
)  # fmt: skip
def test_any_existing_suite_file_is_refused(tmp_path: Path, stale: str) -> None:
    """Should-fix 4: a stale ``xfailed`` outcome or shard-exit line would be merged too; a suite
    without checks never writes ``<suite>.jsonl``, so every ``<suite>.*`` file is refused."""
    logdir = tmp_path / "old"
    logdir.mkdir()
    (logdir / stale).write_text("{}\n")
    rc, calls, out = run(
        tmp_path, "discordance", "3", env={"GIPPA_CALIBRATION_LOG_DIR": str(logdir)}
    )
    assert rc == 2 and calls == [], out


def test_other_suites_files_are_not_stale(tmp_path: Path) -> None:
    logdir = tmp_path / "shared"
    logdir.mkdir()
    (logdir / "regions.jsonl").write_text("{}\n")
    (logdir / "discordance_notes.txt").write_text("x\n")
    rc, calls, out = run(
        tmp_path, "discordance", "2", env={"GIPPA_CALIBRATION_LOG_DIR": str(logdir)}
    )
    assert rc == 0 and len(calls) == 3, out


def test_a_comma_in_the_log_dir_reaches_the_jobs_whole(tmp_path: Path) -> None:
    """Nit: ``--export=ALL,GIPPA_CALIBRATION_LOG_DIR=/x/a,b`` hands the job ``/x/a``. Either refuse
    the path (exit 2, nothing submitted) or pass it so that sbatch gives every job all of it."""
    logdir = tmp_path / "a,b"
    rc, calls, out = run(tmp_path, "discordance", "2",
                         env={"GIPPA_CALIBRATION_LOG_DIR": str(logdir)})  # fmt: skip
    if rc == 2:
        assert calls == [], out
        return
    assert rc == 0 and len(calls) == 3, out
    for call, env in zip(calls, submit_envs(tmp_path), strict=True):
        assert exported_log_dir(call, env) == str(logdir), (opt(call, "export"), env)


def _job(i: int, name: str, cpus: int, state: str = "RUNNING", user: str = "tester") -> dict:
    return {"id": str(500 + i), "name": name, "cpus": cpus, "state": state, "user": user}


#: jobs that never count: another user's gippa-* job, the user's own non-gippa job
NOISE = [_job(90, "gippa-someone-else", 200, "PENDING", user="other"), _job(91, "jupyter", 100)]


def _budget(tmp_path: Path, rows: list[dict], n: int, **env: str) -> tuple[int, list, str]:
    path = tmp_path / "squeue_rows.json"
    path.write_text(json.dumps([*NOISE, *rows]))
    return run(tmp_path, "discordance", str(n), env={"STUB_SQUEUE_ROWS": str(path), **env})


@pytest.mark.parametrize(
    ("rows", "n", "ok"),
    [
        ([_job(0, "gippa-run-x", 96)], 7, True),  # 96 + 8 jobs x 4 CPUs = 128
        ([_job(0, "gippa-run-x", 96)], 8, False),  # 96 + 36 = 132 CPUs
        ([_job(0, "gippa-run-x", 96, "PENDING")], 8, False),  # pending CPUs count too
        ([_job(0, "gippa-a", 60), _job(1, "gippa-b", 40, "PENDING")], 6, True),  # 100 + 28
        ([_job(0, "gippa-a", 60), _job(1, "gippa-b", 40, "PENDING")], 7, False),  # 100 + 32
        ([_job(i, f"gippa-t{i}", 1) for i in range(15)], 4, True),  # 15 + 5 = 20 jobs
        ([_job(i, f"gippa-t{i}", 1) for i in range(15)], 5, False),  # 21 jobs
        ([_job(i, f"gippa-t{i}", 1, "PENDING") for i in range(15)], 5, False),
        ([_job(i, f"gippa-t{i}", 1, "PENDING" if i % 2 else "RUNNING") for i in range(19)], 1,
         False),  # 19 + 2 = 21 jobs
        ([], 10, True),  # only noise: 0 + 11 jobs x 4 = 44 CPUs
    ],
)  # fmt: skip
def test_r14_budget_is_checked_before_any_submission(
    tmp_path: Path, rows: list[dict], n: int, ok: bool
) -> None:
    """Should-fix 5 (R14): sum the user's running + pending ``gippa-*`` jobs and CPUs (squeue) and
    refuse when the N + 1 new jobs (4 CPUs each, calibration_full.sbatch) would pass 128 CPUs or
    20 jobs. Another user's jobs and the user's non-gippa jobs do not count."""
    rc, calls, out = _budget(tmp_path, rows, n)
    assert (tmp_path / "squeue.jsonl").exists(), "the script never queried squeue"
    if ok:
        assert rc == 0 and len(calls) == n + 1, (rc, out)
    else:
        assert rc != 0 and calls == [], (rc, calls, out)


def test_r14_budget_fails_closed_when_squeue_fails(tmp_path: Path) -> None:
    rc, calls, out = _budget(tmp_path, [], 2, STUB_SQUEUE_FAIL="1")
    assert rc != 0 and calls == [], (rc, calls, out)
