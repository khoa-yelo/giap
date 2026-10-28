"""``scripts/slurm/calibration_suite.sbatch``: one Slurm job per calibration suite, run as lanes
(issue M11-03c, PR 3).

Spec: docs/specs/M11-03c-calibration-full-sharding.md, "``calibration_suite.sbatch`` (new)" steps
1-5, "Thread limits per lane (S1)", "In-job partition check (S2)", "Exit folding", "Record",
"R4 under a sharded suite (S3)"; examples 2-5 and 10; acceptance criteria 4, 5 and the job side of
3 (a lane whose collected set or JUnit count differs from ``expected`` gives exit 2).

Static checks of the script text, then runs under bash in a sandbox (the pattern of
tests/unit/test_m11_03b_sbatch.py). Nothing is submitted: ``sbatch``, ``srun``, ``squeue``,
``scancel``, ``uv`` and a bare ``pytest`` on ``PATH`` are fakes that log a forbidden call.

The sandbox:

- ``GIPPA_REPO`` = ``GIPPA_WORKTREE`` = a git checkout in a temporary directory holding a stub
  ``scripts/farmshare_env.sh``, a copy of ``scripts/ci/calibration_lanes.py`` (PR 1a/1b), a plan
  ``scripts/slurm/calibration_plan.tsv`` (the committed one with the determinants row of the
  spec's example 2, ``0 5 *``: the committed row is 12 shards + 3 direct groups since dry run
  45bc3599 and the example's all-direct collection would leave no group-only lanes to check;
  a test may edit a row), a
  durations file ``scripts/slurm/calibration_durations.tsv`` written by the test,
  ``tests/calibration/tolerances.yaml`` and a placeholder ``tests/calibration/_criteria.py``.
  The job must use the checkout's files (in Slurm the batch script runs from the spool directory).
- ``UV_PROJECT_ENVIRONMENT/bin/python`` is a stub interpreter (spec step 1: every Python step runs
  with the venv interpreter, never ``uv run``). A call naming ``tests/calibration/_criteria.py``
  acts as the harness wrapper; ``-m pytest`` acts as pytest; anything else (``-c``, ``-``,
  ``scripts/ci/calibration_lanes.py``, ``-m <other>``) is ``exec``'d with the real interpreter.
  The wrapper runs its pytest as a child with ``subprocess.run``, as ``_criteria.py`` does, so on
  SIGINT the wrapper kills its child 0.25 s after its own ``KeyboardInterrupt`` (spec: the job
  signals pytest, not the wrapper). Both are Python, so SIGINT left ignored by bash for an
  asynchronous command never raises ``KeyboardInterrupt`` (spec: ``trap - INT``).
- The stub pytest reads its collection from ``STUB_COLLECTION`` (``{nodeid: shard0_only}``) and
  applies ``--deselect`` from ``PYTEST_ADDOPTS`` (``shlex.split``) and the command line as pytest
  does: a prefix match. ``--collect-only`` prints the kept nodeids, a blank line and ``N tests
  collected``, writes ``--calibration-lane-report`` when given (exit 5 when nothing is kept); a run
  writes one JUnit testcase per kept nodeid at the last ``--junitxml`` (command line after
  ``PYTEST_ADDOPTS``). Its role is ``collect``, ``merge`` (wrapper with ``--merge``), or the part
  name of its JUnit: ``shard<k>``, ``group<g>``, ``lane<k>`` (a block-2 pytest of ``--merge``).
  Per role: ``STUB_RC_<role>`` exit, ``STUB_SLEEP_<role>`` seconds, ``STUB_FAIL_<role>`` a nodeid
  that fails, ``STUB_DROP_<role>`` one testcase fewer, ``STUB_NO_JUNIT_<role>`` no JUnit.
- The stub ``--merge`` (PR 2 behaviour): with ``GIPPA_CALIBRATION_CONFIRM_LANES=L`` it writes
  ``<log>/<suite>.confirm_parts.json`` (unless ``STUB_NO_CONFIRM_PARTS``) and, with
  ``STUB_PENDING``, starts ``L`` block-2 pytests with
  ``--junitxml=<CONFIRM_JUNIT>.lane<k>.junit.xml`` on their command line; their names are the
  JSON list.
- Review round (#389): ``STUB_DIE_ONCE=<word>`` makes the first stub-interpreter call whose
  argv holds ``<word>`` SIGKILL itself (the job's record helper dying, S1);
  ``STUB_SPAWN_DELAY_<role>`` delays a wrapper's pytest child (the interrupt rescan). The fold of
  a hand-built job state (S2) calls the job's Python helper (the ``<<'PY'`` here-document of the
  script) directly, using its state files (``meta.json``, ``plan.json``, ``lanes.json``,
  ``launched``, ``exit.<pos>``, ``merge.launched``, ``exit.merge``): implementation-specific by
  design (reviewer's suggestion).
- ``taskset`` is a fake that accepts ``-c|--cpu-list LIST cmd…``, exports ``STUB_TASKSET=LIST`` and
  ``exec``s the command.

Conventions settled by the test author where the spec leaves them open (the implementer must
follow them or say why not):

- Test seams of the job: ``GIPPA_CALIBRATION_CPUSET`` (a Linux CPU list, e.g. ``0-3,8-11``)
  replaces the job's ``sched_getaffinity`` set, and ``GIPPA_CALIBRATION_SYSFS`` replaces
  ``/sys/devices/system/cpu`` (``cpu<n>/topology/thread_siblings_list``).
  ``GIPPA_CALIBRATION_GRACE``
  is the spec's grace in seconds. The job reads ``GIPPA_SUITE``, ``GIPPA_WORKTREE``,
  ``GIPPA_CALIBRATION_LOG_DIR`` and ``SLURM_JOB_ID``; there is no ``SLURM_ARRAY_*``.
- A lane is pinned to ``cpus_per_lane / 2`` whole cores: its ``taskset`` list holds both siblings
  of each of its cores, all inside the job's set, disjoint from the other lanes'. Too few whole
  cores for the plan row's ``shards + groups`` lanes ends the job with exit 2 and a record before
  any lane starts.
- The lane processes (wrapper and pytest) see ``OMP_NUM_THREADS``, ``OPENBLAS_NUM_THREADS``,
  ``MKL_NUM_THREADS``, ``NUMBA_NUM_THREADS`` and ``SLURM_CPUS_PER_TASK`` = ``cpus_per_lane``.
- A shard lane's wrapper gets ``--shard k/N``; a group lane's wrapper gets neither ``--shard`` nor
  ``--merge``. The merge starts after every lane has ended (step 4, "Wait for all lanes").
- The record's ``parts`` equal ``calibration_lanes.record_parts(lanes, confirm_parts)``; the merge
  part names come from ``<suite>.confirm_parts.json``. A missing file after a merge exit 0 or 1 is
  a check problem (exit 2); a JUnit part that is missing or has another testcase count than its
  lane's ``expected`` is a count problem (exit 2). ``checks`` is non-empty exactly when the job
  found a partition, count or confirm-parts problem.
- On USR1 the job's lane pytests get SIGINT (a pytest grandchild of the batch shell, not its
  ``_criteria.py`` wrapper) within 3 s; the record says ``exit 143``; the job's own exit is
  non-zero. SIGTERM also gives ``exit 143``.

No statistical claim is made here (scheduling only, R3), so no calibration test applies. The fold
rule itself is property-tested in tests/unit/test_m11_03c_plan.py; here a matrix of lane and merge
exits pins that the job applies it (each case is a full job run, too slow for Hypothesis).
"""

from __future__ import annotations

import contextlib
import importlib.util
import json
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from helpers.slurm_sandbox import assert_honours_gippa_repo

REPO = Path(__file__).resolve().parents[2]
SBATCH = REPO / "scripts" / "slurm" / "calibration_suite.sbatch"
LANES_PY = REPO / "scripts" / "ci" / "calibration_lanes.py"
ROWS_PY = REPO / "scripts" / "ci" / "calibration_rows.py"
PLAN = REPO / "scripts" / "slurm" / "calibration_plan.tsv"
SUITES = ["placement", "gene_fdr", "detectability", "regions", "autocorrelation", "islands",
          "discordance", "missed_operons", "ora", "gsea", "determinants", "novel_orfs",
          "denovo_precision", "variants", "compare"]  # fmt: skip
JOB = "1900100"
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS",
               "SLURM_CPUS_PER_TASK")  # fmt: skip
CONFIRM_VARS = ("GIPPA_CALIBRATION_CONFIRM_LANES", "GIPPA_CALIBRATION_CONFIRM_JUNIT",
                "GIPPA_CALIBRATION_CONFIRM")  # fmt: skip

DET = "tests/calibration/attacks/test_attack_determinants_null.py::test_null_level_under_attack"
EXAMPLE2_ROW = "determinants\t0\t5\t*\t2\t4000\t08:00:00\t30"
ATTACK5 = ("tests/calibration/attacks/test_attack5_reversed_stretches.py::"
           "test_attack5_reversed_stretches_shuffle_decoys_full")  # fmt: skip
VAR = "tests/calibration/test_variants_calibration.py"
VAR_ATT = "tests/calibration/attacks/test_attack_variants.py"

# ------------------------------------------------------------------ stubs

STUB_PY = r'''#!@REAL@
"""stub venv python (tests/unit/test_m11_03c_suite_sbatch.py)."""
import json, os, re, shlex, signal, subprocess, sys, time
from pathlib import Path

REAL = "@REAL@"
argv = sys.argv[1:]
i = 0
while i < len(argv) and argv[i].startswith("-") and argv[i] not in ("-m", "-c", "-"):
    i += 1
_die = os.environ.get("STUB_DIE_ONCE")
if _die and _die in argv and not os.path.exists(os.environ["STUB_LOG"] + ".died"):
    open(os.environ["STUB_LOG"] + ".died", "w").close()
    os.kill(os.getpid(), signal.SIGKILL)  # the first matching call dies (review S1)
mode = None
if i + 1 < len(argv) and argv[i] == "-m" and argv[i + 1] == "pytest":
    mode, rest = "pytest", argv[i + 2:]
elif i < len(argv) and argv[i].endswith("_criteria.py"):
    mode, rest = "criteria", argv[i + 1:]
if mode is None:
    os.execv(REAL, [REAL, *argv])

ENV_KEYS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS",
            "SLURM_CPUS_PER_TASK", "STUB_TASKSET", "GIPPA_CALIBRATION_LOG_DIR",
            "GIPPA_CALIBRATION_CONFIRM_LANES", "GIPPA_CALIBRATION_CONFIRM_JUNIT",
            "GIPPA_CALIBRATION_CONFIRM", "GIPPA_CALIBRATION_REPLICATES", "PYTEST_ADDOPTS",
            "GIPPA_CALIBRATION_CONFIRM_SLICE")
addopts = shlex.split(os.environ.get("PYTEST_ADDOPTS", ""))
opts = addopts + rest


def values(name):
    out, j = [], 0
    while j < len(opts):
        a = opts[j]
        if a == name and j + 1 < len(opts):
            out.append(opts[j + 1]); j += 2; continue
        if a.startswith(name + "="):
            out.append(a[len(name) + 1:])
        j += 1
    return out


deselect = values("--deselect")
junits = values("--junitxml") + values("--junit-xml")
junit = junits[-1] if junits else None
collect = "--collect-only" in opts or "--co" in opts
report = (values("--calibration-lane-report") or [None])[-1]


def role():
    if collect:
        return "collect"
    if mode == "criteria" and "--merge" in rest:
        return "merge"
    m = re.search(r"\.(shard\d+|group\d+|lane\d+)\.junit\.xml$", junit or "")
    return m.group(1) if m else "unknown"


ROLE = role()


def env(name, default=""):
    return os.environ.get(f"STUB_{name}_{ROLE}", default)


def log(**kw):
    kw.update(pid=os.getpid(), ppid=os.getppid(), t=time.time(), role=ROLE, mode=mode,
              argv=rest, deselect=deselect, junit=junit, collect=collect,
              env={k: os.environ.get(k) for k in ENV_KEYS})
    fd = os.open(os.environ["STUB_LOG"], os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
    os.write(fd, (json.dumps(kw) + "\n").encode())
    os.close(fd)


collection = json.loads(Path(os.environ["STUB_COLLECTION"]).read_text())
pref = tuple(deselect)
kept = [n for n in sorted(collection) if not (pref and n.startswith(pref))]


def write_junit(path, cases, fails):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    body = []
    for n in cases:
        mod, _, name = n.partition("::")
        cls = mod[:-3].replace("/", ".") if mod.endswith(".py") else mod
        name = name.replace("&", "&amp;").replace('"', "&quot;").replace("<", "&lt;")
        inner = '<failure message="stub">x</failure>' if n in fails else ""
        body.append(f'<testcase classname="{cls}" name="{name}" time="0.1">{inner}</testcase>')
    Path(path).write_text('<?xml version="1.0" encoding="utf-8"?><testsuites>'
                          f'<testsuite name="pytest" tests="{len(cases)}">' + "".join(body)
                          + "</testsuite></testsuites>")


def nap(seconds):
    end = time.time() + seconds
    while time.time() < end:
        time.sleep(0.05)


def pytest_main():
    if os.environ.get("STUB_IGNORE_INT"):
        signal.signal(signal.SIGINT, signal.SIG_IGN)
    log(event="start", proc="pytest")
    fails = {f for f in env("FAIL").split(",") if f}
    cases = list(kept)
    if env("DROP"):
        cases = cases[:-1]
    try:
        nap(float(env("SLEEP", "0")))
    except KeyboardInterrupt:
        log(event="int", proc="pytest")
        if junit and not collect:
            write_junit(junit, cases[:1], fails)
        return 2
    if collect:
        hide = os.environ.get("STUB_COLLECT_HIDE")
        shown = [n for n in kept if not (hide and report is None and n == hide)]
        if report is not None:
            Path(report).write_text(json.dumps({n: {"shard0_only": collection[n]}
                                                for n in kept}))
        print("\n".join(shown))
        print()
        print(f"{len(shown)} tests collected in 0.01s")
        rc = int(env("RC", "0")) or (0 if shown else 5)
        log(event="end", proc="pytest", rc=rc)
        return rc
    if junit and not env("NO_JUNIT"):
        write_junit(junit, cases, fails)
    rc = int(env("RC", "1" if fails else "0"))
    log(event="end", proc="pytest", rc=rc)
    return rc


def child(extra):
    return [sys.argv[0], "-m", "pytest", "tests/calibration", "-m", "calibration_full", *extra]


def merge_main():
    log(event="start", proc="criteria")
    lanes = os.environ.get("GIPPA_CALIBRATION_CONFIRM_LANES")
    procs = []
    if lanes is not None:
        stem = os.environ["GIPPA_CALIBRATION_CONFIRM_JUNIT"]
        paths = [f"{stem}.lane{k}.junit.xml" for k in range(int(lanes))]
        if not os.environ.get("STUB_PENDING"):
            paths = []
        logdir = Path(os.environ["GIPPA_CALIBRATION_LOG_DIR"])
        if not os.environ.get("STUB_NO_CONFIRM_PARTS"):
            logdir.mkdir(parents=True, exist_ok=True)
            (logdir / f"{os.environ['STUB_SUITE']}.confirm_parts.json").write_text(
                json.dumps([Path(p).name for p in paths]) + "\n")
        procs = [subprocess.Popen(child([f"--junitxml={p}"])) for p in paths]
    try:
        codes = [p.wait() for p in procs]
    except KeyboardInterrupt:
        log(event="int", proc="criteria")
        raise
    rc = int(os.environ.get("STUB_RC_merge", "0"))
    if any(codes):
        rc = 2
    log(event="end", proc="criteria", rc=rc)
    return rc


def criteria_main():
    if ROLE == "merge":
        return merge_main()
    log(event="start", proc="criteria")
    try:
        nap(float(os.environ.get(f"STUB_SPAWN_DELAY_{ROLE}", "0")))
        rc = subprocess.run(child([]), check=False).returncode
    except KeyboardInterrupt:
        log(event="int", proc="criteria")
        raise
    log(event="end", proc="criteria", rc=rc)
    return rc


sys.exit(pytest_main() if mode == "pytest" else criteria_main())
'''

FAKE_TASKSET = r"""#!/bin/bash
# fake taskset: -c|--cpu-list LIST cmd...
case "${1:-}" in
  -c|--cpu-list) export STUB_TASKSET="$2"; shift 2 ;;
  --cpu-list=*) export STUB_TASKSET="${1#--cpu-list=}"; shift ;;
  *) echo "fake taskset: unsupported arguments: $*" >&2; exit 64 ;;
esac
exec "$@"
"""

FORBIDDEN = r"""#!/bin/bash
printf '{"event": "forbidden", "cmd": "%s", "args": "%s"}\n' "$(basename "$0")" "$*" >> "$STUB_LOG"
exit 1
"""

_GIT_ENV = {
    "GIT_AUTHOR_NAME": "Fixture", "GIT_AUTHOR_EMAIL": "fixture@example.invalid",
    "GIT_COMMITTER_NAME": "Fixture", "GIT_COMMITTER_EMAIL": "fixture@example.invalid",
    "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1",
}  # fmt: skip


def _git(repo: Path, *args: str) -> str:
    out = subprocess.run(["git", "-C", str(repo), *args], env={**os.environ, **_GIT_ENV},
                         check=True, capture_output=True, text=True)  # fmt: skip
    return out.stdout.strip()


def _load(path: Path, name: str) -> ModuleType:
    assert path.is_file(), f"{path} does not exist"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def lanes_mod() -> ModuleType:
    return _load(LANES_PY, "gippa_calibration_lanes_pr3")


@pytest.fixture(scope="module")
def rows_mod() -> ModuleType:
    return _load(ROWS_PY, "gippa_calibration_rows_pr3")


def _code() -> str:
    assert SBATCH.is_file(), f"{SBATCH} does not exist (M11-03c PR 3)"
    return "\n".join(ln for ln in SBATCH.read_text().splitlines()
                     if not ln.lstrip().startswith("#"))  # fmt: skip


# ------------------------------------------------------------------ sandbox


@dataclass
class Run:
    rc: int
    stdout: str
    stderr: str
    events: list[dict[str, Any]]
    log: Path
    elapsed: float = 0.0

    def of(self, *, proc: str | None = None, event: str = "start",
           role: str | None = None) -> list[dict[str, Any]]:  # fmt: skip
        return [e for e in self.events if e.get("event") == event
                and (proc is None or e.get("proc") == proc)
                and (role is None or e.get("role") == role)]  # fmt: skip

    def lane_wrappers(self) -> list[dict[str, Any]]:
        return [e for e in self.of(proc="criteria")
                if re.fullmatch(r"shard\d+|group\d+", e["role"])]  # fmt: skip

    def collects(self) -> list[dict[str, Any]]:
        return self.of(proc="pytest", role="collect")

    def forbidden(self) -> list[dict[str, Any]]:
        return [e for e in self.events if e.get("event") == "forbidden"]


@dataclass
class Box:
    base: Path
    repo: Path
    head: str
    log: Path
    bin: Path
    venv: Path
    sysfs: Path
    stub_log: Path
    collection: Path
    procs: list[subprocess.Popen[bytes]] = field(default_factory=list)

    def plan_row(self, suite: str, cells: str) -> None:
        """Replace the plan row of ``suite`` (tab-separated cells after the suite name)."""
        plan = self.repo / "scripts" / "slurm" / "calibration_plan.tsv"
        lines = plan.read_text().splitlines()
        lines = [f"{suite}\t{cells}" if ln.split("\t")[0] == suite else ln for ln in lines]
        plan.write_text("\n".join(lines) + "\n")

    def durations(self, secs: Mapping[str, float]) -> None:
        body = "".join(f"{n}\t{s:g}\tprobe\n" for n, s in sorted(secs.items()))
        (self.repo / "scripts" / "slurm" / "calibration_durations.tsv").write_text(
            "nodeid\tseconds\tsource\n" + body
        )

    def topology(self, ncpu: int, style: str = "split") -> None:
        """``split``: cpu n and n + ncpu/2 are siblings (``0,10``); ``pairs``: 2m, 2m+1."""
        shutil.rmtree(self.sysfs, ignore_errors=True)
        half = ncpu // 2
        for n in range(ncpu):
            if style == "split":
                sib = f"{n % half},{n % half + half}"
            else:
                sib = f"{n - n % 2}-{n - n % 2 + 1}"
            d = self.sysfs / f"cpu{n}" / "topology"
            d.mkdir(parents=True)
            (d / "thread_siblings_list").write_text(sib + "\n")

    def siblings(self, cpu: int) -> set[int]:
        text = (self.sysfs / f"cpu{cpu}" / "topology" / "thread_siblings_list").read_text()
        return _cpulist(text.strip())

    def env(self, suite: str, report: Mapping[str, bool], cpuset: str,
            extra: Mapping[str, str] | None = None) -> dict[str, str]:  # fmt: skip
        self.collection.write_text(json.dumps(dict(report)))
        return {
            "PATH": f"{self.bin}{os.pathsep}{os.environ['PATH']}",
            "HOME": str(self.base),
            "GIPPA_REPO": str(self.repo),
            "GIPPA_WORKTREE": str(self.repo),
            "GIPPA_ROOT": str(self.base),
            "GIPPA_SUITE": suite,
            "GIPPA_CALIBRATION_LOG_DIR": str(self.log),
            "GIPPA_CALIBRATION_CPUSET": cpuset,
            "GIPPA_CALIBRATION_SYSFS": str(self.sysfs),
            "GIPPA_CALIBRATION_GRACE": "3",
            "UV_PROJECT_ENVIRONMENT": str(self.venv),
            "SLURM_JOB_ID": JOB,
            "SLURM_JOB_NAME": f"gippa-cal-{suite}",
            "SLURM_CPUS_PER_TASK": str(len(_cpulist(cpuset))),
            "STUB_LOG": str(self.stub_log),
            "STUB_COLLECTION": str(self.collection),
            "STUB_SUITE": suite,
            **_GIT_ENV,
            **(extra or {}),
        }

    def start(self, env: Mapping[str, str]) -> tuple[subprocess.Popen[bytes], Path, Path]:
        out, err = self.base / "job.out", self.base / "job.err"
        bash = shutil.which("bash") or "/bin/bash"
        with out.open("wb") as fo, err.open("wb") as fe:  # files, not pipes: lingering stubs
            proc = subprocess.Popen([bash, str(SBATCH)], env=dict(env), stdout=fo, stderr=fe,
                                    cwd=self.base)  # fmt: skip
        self.procs.append(proc)
        return proc, out, err

    def events(self) -> list[dict[str, Any]]:
        if not self.stub_log.exists():
            return []
        return [json.loads(ln) for ln in self.stub_log.read_text().splitlines() if ln.strip()]

    def finish(self, proc: subprocess.Popen[bytes], out: Path, err: Path, t0: float,
               timeout: float = 120) -> Run:  # fmt: skip
        try:
            rc = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            pytest.fail(f"job did not finish in {timeout} s; stderr:\n{err.read_text()}")
        return Run(rc, out.read_text(), err.read_text(), self.events(), self.log,
                   time.monotonic() - t0)  # fmt: skip

    def run(self, suite: str, report: Mapping[str, bool], cpuset: str,
            extra: Mapping[str, str] | None = None) -> Run:  # fmt: skip
        assert SBATCH.is_file(), f"{SBATCH} does not exist (M11-03c PR 3)"
        t0 = time.monotonic()
        proc, out, err = self.start(self.env(suite, report, cpuset, extra))
        return self.finish(proc, out, err, t0)


def _cpulist(text: str) -> set[int]:
    out: set[int] = set()
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-")
            out.update(range(int(a), int(b) + 1))
        else:
            out.add(int(part))
    return out


def _make_box(base: Path) -> Box:
    repo = base / "repo"
    for d in ("scripts/ci", "scripts/slurm", "tests/calibration", "tests/sim"):
        (repo / d).mkdir(parents=True)
    (repo / "scripts" / "farmshare_env.sh").write_text("# stub\n")
    shutil.copyfile(LANES_PY, repo / "scripts" / "ci" / "calibration_lanes.py")
    plan = PLAN.read_text().splitlines()
    for k, ln in enumerate(plan):  # example 2's row: 5 group lanes of the 19 direct tests
        if ln.startswith("determinants\t"):
            plan[k] = "\t".join(["determinants", "0", "5", "*", *ln.split("\t")[4:]])
    assert EXAMPLE2_ROW in plan, plan
    (repo / "scripts" / "slurm" / "calibration_plan.tsv").write_text("\n".join(plan) + "\n")
    (repo / "scripts" / "slurm" / "calibration_durations.tsv").write_text(
        "nodeid\tseconds\tsource\n"
    )
    shutil.copyfile(REPO / "tests" / "calibration" / "tolerances.yaml",
                    repo / "tests" / "calibration" / "tolerances.yaml")  # fmt: skip
    (repo / "tests" / "calibration" / "_criteria.py").write_text("# placeholder: stub python\n")
    (repo / "tests" / "sim" / "keep.py").write_text("# sim\n")
    _git(repo, "init", "-q", "-b", "main")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "fixture")
    bindir, venv = base / "bin", base / "venv"
    (venv / "bin").mkdir(parents=True)
    bindir.mkdir()
    stub = venv / "bin" / "python"
    stub.write_text(STUB_PY.replace("@REAL@", sys.executable))
    stub.chmod(0o755)
    (bindir / "taskset").write_text(FAKE_TASKSET)
    (bindir / "taskset").chmod(0o755)
    for name in ("sbatch", "srun", "squeue", "scancel", "uv", "pytest"):
        (bindir / name).write_text(FORBIDDEN)
        (bindir / name).chmod(0o755)
    box = Box(base, repo, _git(repo, "rev-parse", "HEAD"), base / "log", bindir, venv,
              base / "sysfs", base / "stub.jsonl", base / "collection.json")  # fmt: skip
    box.topology(20)
    return box


def _kill_stubs(box: Box) -> None:
    for e in box.events():
        pid = e.get("pid")
        if not isinstance(pid, int):
            continue
        with contextlib.suppress(OSError):
            cmd = Path(f"/proc/{pid}/cmdline").read_bytes().decode(errors="replace")
            if str(box.base) in cmd:
                os.kill(pid, signal.SIGKILL)
    for p in box.procs:
        if p.poll() is None:
            p.kill()
            p.wait()


@pytest.fixture
def box() -> Iterator[Box]:
    # not tmp_path: its path contains "pytest", which would make every stub look like pytest
    base = Path(tempfile.mkdtemp(prefix="gippa-m1103c-"))
    b = _make_box(base)
    try:
        yield b
    finally:
        _kill_stubs(b)
        shutil.rmtree(base, ignore_errors=True)


# ------------------------------------------------------------------ collections


def det_report(n: int = 19) -> tuple[dict[str, bool], dict[str, float]]:
    """Example 2: 19 direct determinants tests; two long ones."""
    params = [f"{p}-{m}" for p in ("autocorr_pos", "detect", "size", "clustered", "gc", "len",
                                   "codon", "phase", "tss", "gene")
              for m in ("lightgbm", "elastic_net")]  # fmt: skip
    ids = [f"{DET}[full-{p}]" for p in params][:n]
    durations = {nid: 600.0 + 37.0 * j for j, nid in enumerate(sorted(ids))}
    durations[f"{DET}[full-autocorr_pos-lightgbm]"] = 8364.0
    durations[f"{DET}[full-detect-lightgbm]"] = 7782.0
    return dict.fromkeys(ids, True), durations


def placement_report() -> tuple[dict[str, bool], dict[str, float]]:
    """Example 3: 28 shard-0-only tests, 30 harness tests (attack5 among them)."""
    direct = [f"tests/calibration/test_placement_fdr_sim.py::test_direct_{j:02d}"
              for j in range(28)]  # fmt: skip
    harness = [f"tests/calibration/test_placement_harness.py::test_check_{j:02d}_full"
               for j in range(29)] + [ATTACK5]  # fmt: skip
    durations = {nid: 4140.0 + 60.0 * (j % 11) for j, nid in enumerate(direct)}
    durations[ATTACK5] = 12096.0
    return {**dict.fromkeys(direct, True), **dict.fromkeys(harness, False)}, durations


VAR_DIRECT = [f"{VAR}::test_direct_a", f"{VAR}::test_direct_b", f"{VAR}::test_direct_c[case one]"]
VAR_HARNESS = [f"{VAR_ATT}::test_g1_pooled_full", f"{VAR_ATT}::test_a6b_full",
               f"{VAR}::test_level_full"]  # fmt: skip


def var_report() -> dict[str, bool]:
    """variants (2 shards, 1 group, direct); one nodeid with a space (shlex quoting)."""
    return {**dict.fromkeys(VAR_DIRECT, True), **dict.fromkeys(VAR_HARNESS, False)}


def expected_lanes(lanes_mod: ModuleType, box: Box, suite: str, report: Mapping[str, bool]
                   ) -> list[Any]:  # fmt: skip
    rows = {r.suite: r for r in lanes_mod.read_plan(
        box.repo / "scripts" / "slurm" / "calibration_plan.tsv", SUITES)}  # fmt: skip
    durations = lanes_mod.read_durations(box.repo / "scripts" / "slurm" /
                                         "calibration_durations.tsv")  # fmt: skip
    return list(lanes_mod.lanes(rows[suite], dict(report), durations))


def record(log: Path, task: int) -> dict[str, Any]:
    recs = sorted(log.glob("task-*.json"))
    assert [p.name for p in recs] == [f"task-{task}.json"], recs
    assert not [p for p in log.iterdir() if ".tmp" in p.name], list(log.iterdir())
    doc: dict[str, Any] = json.loads(recs[0].read_text())
    return doc


def junit_count(path: Path) -> int:
    return sum(1 for _ in ET.parse(path).getroot().iter("testcase"))


def _role_of_part(part: str) -> str:
    m = re.search(r"\.(shard\d+|group\d+)\.junit\.xml$", part)
    assert m, part
    return m.group(1)


def nightly(rows_mod: ModuleType, box: Box, task: int, suite: str) -> tuple[Any, list[str]]:
    res = rows_mod._read_task(box.log, task, box.head, suite, {})
    return res, rows_mod.full_rows([res], JOB, box.head, "2026-10-25", expect=1)


def _no_lane_started(run: Run) -> None:
    assert run.lane_wrappers() == [], [e["role"] for e in run.lane_wrappers()]
    assert run.of(proc="criteria", role="merge") == []
    assert [e for e in run.of(proc="pytest") if e["role"] != "collect"] == []


# ------------------------------------------------------------------ static


def test_script_exists_and_parses() -> None:
    assert SBATCH.is_file(), f"{SBATCH} does not exist (M11-03c PR 3)"
    res = subprocess.run(["bash", "-n", str(SBATCH)], capture_output=True, text=True, check=False)
    assert res.returncode == 0, res.stderr


def test_signal_directive_and_no_array() -> None:
    _code()
    text = SBATCH.read_text()
    assert re.findall(r"^#SBATCH\s+--signal=(\S+)\s*$", text, re.M) == ["B:USR1@600"]
    assert not re.search(r"^#SBATCH\s+--array", text, re.M)
    names = re.findall(r"^#SBATCH\s+--job-name=(\S+)", text, re.M)
    assert all(n.startswith("gippa-") for n in names), names


def test_venv_interpreter_never_uv_run() -> None:
    code = _code()
    assert not re.search(r"\buv\s+run\b", code), "spec step 1: no uv between the shell and pytest"
    assert "UV_PROJECT_ENVIRONMENT" in code and "/bin/python" in code


def test_traps_usr1_term_exit_and_restores_sigint() -> None:
    code = _code()
    traps = re.findall(r"^\s*(?:.*[;&(]\s*)?trap\s+(.+)$", code, re.M)
    joined = " ".join(traps)
    for sig in ("USR1", "TERM", "EXIT"):
        assert re.search(rf"\b{sig}\b", joined), (sig, traps)
    assert re.search(r"trap\s+-\s+INT\b", code), "an asynchronous lane must get SIGINT back"


def test_signals_only_its_own_tree() -> None:
    code = _code()
    assert not re.search(r"\b(pkill|killall)\b", code)
    assert not re.search(r"-u\s+\"?\$\{?USER\b", code)
    assert not re.search(r"\bkill\b[^;\n]*\s0(\s|;|$)", code, re.M), "no process-group kill"
    assert not re.search(r"\bkill\b[^;\n]*\s--\s+-", code), "no process-group kill"
    assert "/proc/" in code or "pgrep -P" in code, "walk the batch shell's descendants"


def test_honours_gippa_repo_and_bash32() -> None:
    assert_honours_gippa_repo(SBATCH.read_text(), "calibration_suite.sbatch")
    res = subprocess.run([sys.executable, str(REPO / "scripts" / "ci" / "check_bash32.py"),
                          str(SBATCH)], capture_output=True, text=True, check=False)  # fmt: skip
    assert res.returncode == 0, res.stdout + res.stderr


def test_shellcheck_when_available() -> None:
    sc = shutil.which("shellcheck")
    if sc is None:
        pytest.skip("shellcheck not installed")
    res = subprocess.run([sc, "-S", "warning", str(SBATCH)], capture_output=True, text=True,
                         check=False)  # fmt: skip
    assert res.returncode == 0, res.stdout


# ------------------------------------------------------------------ example 2 (group lanes only)


def test_example_2_determinants_group_lanes(box: Box, lanes_mod: ModuleType,
                                            rows_mod: ModuleType) -> None:  # fmt: skip
    report, durations = det_report()
    box.durations(durations)
    want = expected_lanes(lanes_mod, box, "determinants", report)
    assert [ln.kind for ln in want] == ["group"] * 5
    run = box.run("determinants", report, "0-19")
    assert run.rc == 0, run.stderr
    assert run.forbidden() == [], run.forbidden()
    log = run.log

    # collection: one lane-report run, then one --collect-only per lane, all before any lane
    collects = run.collects()
    assert len(collects) == 1 + len(want), [c["deselect"] for c in collects]
    assert sum(1 for c in collects if not c["deselect"]) >= 1
    for ln in want:
        assert any(set(c["deselect"]) == set(ln.deselect) for c in collects), ln.index
    ends = [e["t"] for e in run.of(proc="pytest", role="collect", event="end")]
    wrappers = run.lane_wrappers()
    assert len(ends) == len(collects) and wrappers
    assert max(ends) <= min(w["t"] for w in wrappers), "collect-only checks run before any lane"

    # every lane: its own JUnit, deselect list, log dir, core, thread limits
    assert sorted(w["role"] for w in wrappers) == [f"group{g}" for g in range(5)]
    cores: list[set[int]] = []
    for ln in want:
        (w,) = [w for w in wrappers if w["role"] == f"group{ln.index}"]
        assert w["junit"] == f"{log}/task-10.group{ln.index}.junit.xml", w
        assert set(w["deselect"]) == set(ln.deselect)
        assert not {"--shard", "--merge"} & {a.split("=")[0] for a in w["argv"]}, w["argv"]
        assert w["env"]["GIPPA_CALIBRATION_LOG_DIR"] == f"{log}/lanes/determinants-group{ln.index}"
        (inner,) = run.of(proc="pytest", role=f"group{ln.index}")
        assert inner["ppid"] == w["pid"]
        for e in (w, inner):
            assert {k: e["env"][k] for k in THREAD_VARS} == dict.fromkeys(THREAD_VARS, "2"), e
            assert not any(e["env"][k] for k in CONFIRM_VARS), e["env"]
        cpus = _cpulist(w["env"]["STUB_TASKSET"] or "")
        assert len(cpus) == 2 and cpus <= set(range(20)), cpus
        assert all(box.siblings(c) <= cpus for c in cpus), f"not a whole core: {cpus}"
        cores.append(cpus)
        assert junit_count(log / f"task-10.group{ln.index}.junit.xml") == len(ln.expected)
    assert all(a.isdisjoint(b) for i, a in enumerate(cores) for b in cores[i + 1 :]), cores

    # no merge without shards
    assert run.of(proc="criteria", role="merge") == []

    rec = record(log, 10)
    assert {k: rec[k] for k in ("task", "suite", "job", "commit", "dirty", "exit", "args",
                                "sharded")} == {
        "task": 10, "suite": "determinants", "job": JOB, "commit": box.head, "dirty": False,
        "exit": 0, "args": [], "sharded": False,
    }  # fmt: skip
    assert rec["parts"] == [f"task-10.group{g}.junit.xml" for g in range(5)]
    assert rec["parts"] == lanes_mod.record_parts(want, [])
    assert re.fullmatch(r"\d{4}-\d\d-\d\d", rec["date"]), rec
    assert [(e["kind"], e["index"], e["exit"], e["tests"], e["junit_tests"])
            for e in rec["lanes"]] == [("group", ln.index, 0, len(ln.expected),
                                        len(ln.expected)) for ln in want]  # fmt: skip
    assert {k: str(rec["plan"][k]) for k in ("suite", "shards", "groups", "group_select")} == {
        "suite": "determinants",
        "shards": "0",
        "groups": "5",
        "group_select": "*",
    }
    assert rec["checks"] == []

    res, rows = nightly(rows_mod, box, 10, "determinants")
    assert res.problems == [] and not res.fail_closed, res
    assert [r.split("\t")[3] for r in rows] == ["0", "0"], rows


def test_empty_group_lanes_are_recorded_with_zero_tests(box: Box, lanes_mod: ModuleType) -> None:
    report, durations = det_report(3)
    box.durations(durations)
    want = expected_lanes(lanes_mod, box, "determinants", report)
    run = box.run("determinants", report, "0-19")
    assert run.rc == 0, run.stderr
    non_empty = [ln for ln in want if ln.expected]
    assert len(non_empty) == 3
    assert sorted(w["role"] for w in run.lane_wrappers()) == sorted(
        f"group{ln.index}" for ln in non_empty
    )
    rec = record(run.log, 10)
    assert rec["parts"] == [ln.part for ln in non_empty]
    assert len(rec["lanes"]) == 5
    assert sorted(e["tests"] for e in rec["lanes"]) == [0, 0, 1, 1, 1]


@pytest.mark.parametrize("style", ["split", "pairs"])
def test_four_cpu_lanes_get_two_whole_cores(box: Box, lanes_mod: ModuleType, style: str) -> None:
    box.topology(8, style)
    report = {f"tests/calibration/test_ora_calibration.py::test_ora_{c}_full": False
              for c in "abcd"}  # fmt: skip
    run = box.run("ora", report, "0-7")
    assert run.rc == 0, run.stderr
    wrappers = run.lane_wrappers()
    assert sorted(w["role"] for w in wrappers) == ["group0", "group1"]
    sets = []
    for w in wrappers:
        assert {k: w["env"][k] for k in THREAD_VARS} == dict.fromkeys(THREAD_VARS, "4"), w["env"]
        cpus = _cpulist(w["env"]["STUB_TASKSET"] or "")
        assert len(cpus) == 4 and all(box.siblings(c) <= cpus for c in cpus), (style, cpus)
        sets.append(cpus)
    assert sets[0].isdisjoint(sets[1])


@pytest.mark.parametrize("cpuset", ["0-9", "0-3,10-13"], ids=["no-whole-core", "four-cores"])
def test_too_few_whole_cores_exit_2_before_any_lane(box: Box, cpuset: str) -> None:
    # determinants needs 5 cores; "0-9" has 10 CPUs but no whole core (siblings n, n + 10)
    report, durations = det_report()
    box.durations(durations)
    run = box.run("determinants", report, cpuset)
    assert run.rc == 2, run.stderr
    _no_lane_started(run)
    rec = record(run.log, 10)
    assert rec["exit"] == 2 and rec["job"] == JOB


# ------------------------------------------------------------------ example 3 (shards + merge)


def test_example_3_placement_shards_groups_and_merge(box: Box, lanes_mod: ModuleType,
                                                     rows_mod: ModuleType) -> None:  # fmt: skip
    report, durations = placement_report()
    box.durations(durations)
    want = expected_lanes(lanes_mod, box, "placement", report)
    direct = {n for n, s in report.items() if s}
    run = box.run("placement", report, "0-19", {"STUB_PENDING": "1"})
    assert run.rc == 0, run.stderr
    assert run.forbidden() == []
    log = run.log
    wrappers = run.lane_wrappers()
    assert sorted(w["role"] for w in wrappers) == sorted(
        [f"shard{k}" for k in range(2)] + [f"group{g}" for g in range(8)]
    )
    for k in range(2):
        (w,) = [w for w in wrappers if w["role"] == f"shard{k}"]
        assert w["junit"] == f"{log}/task-0.shard{k}.junit.xml"
        assert set(w["deselect"]) == direct
        args = " ".join(w["argv"])
        assert re.search(rf"--shard[= ]{k}/2\b", args), w["argv"]
        assert w["env"]["GIPPA_CALIBRATION_LOG_DIR"] == str(log)
        assert junit_count(log / f"task-0.shard{k}.junit.xml") == 30
    for ln in [ln for ln in want if ln.kind == "group"]:
        (w,) = [w for w in wrappers if w["role"] == f"group{ln.index}"]
        assert set(w["deselect"]) == set(ln.deselect)
        assert w["env"]["GIPPA_CALIBRATION_LOG_DIR"] == f"{log}/lanes/placement-group{ln.index}"
    cores = [_cpulist(w["env"]["STUB_TASKSET"] or "") for w in wrappers]
    assert all(len(c) == 2 for c in cores)
    assert len(set().union(*cores)) == 20, "10 lanes on 10 disjoint whole cores"

    (merge,) = run.of(proc="criteria", role="merge")
    lane_ends = [e["t"] for e in run.of(proc="criteria", event="end")
                 if re.fullmatch(r"shard\d+|group\d+", e["role"])]  # fmt: skip
    assert len(lane_ends) == 10 and merge["t"] >= max(lane_ends), "merge after every lane"
    assert "--merge" in merge["argv"]
    assert set(merge["deselect"]) == direct
    assert merge["env"]["GIPPA_CALIBRATION_CONFIRM_LANES"] == "2"
    assert merge["env"]["GIPPA_CALIBRATION_CONFIRM_JUNIT"] == f"{log}/task-0.merge"
    assert merge["env"]["GIPPA_CALIBRATION_LOG_DIR"] == str(log)
    # step 4: the merge runs on the shard lanes' cores with the lanes' thread limits
    shard_cpus = set().union(*(_cpulist(w["env"]["STUB_TASKSET"] or "") for w in wrappers
                               if w["role"].startswith("shard")))  # fmt: skip
    assert len(shard_cpus) == 4
    assert _cpulist(merge["env"]["STUB_TASKSET"] or "") == shard_cpus, merge["env"]
    assert {k: merge["env"][k] for k in THREAD_VARS} == dict.fromkeys(THREAD_VARS, "2")
    assert sorted(e["role"] for e in run.of(proc="pytest") if e["role"].startswith("lane")) == [
        "lane0",
        "lane1",
    ]

    rec = record(log, 0)
    assert rec["sharded"] is True and rec["exit"] == 0 and rec["job"] == JOB, rec
    merge_parts = ["task-0.merge.lane0.junit.xml", "task-0.merge.lane1.junit.xml"]
    assert rec["parts"] == lanes_mod.record_parts(want, merge_parts)
    assert rec["parts"][:2] == ["task-0.shard0.junit.xml", "task-0.shard1.junit.xml"]
    assert rec["parts"][-2:] == merge_parts
    assert rec["checks"] == []
    res, _rows = nightly(rows_mod, box, 0, "placement")
    assert res.problems == [] and not res.fail_closed, res


# ------------------------------------------------------------------ errors before any lane


def test_example_4_dash_with_shard0_only_tests_exits_2(box: Box) -> None:
    box.plan_row("novel_orfs", "6\t0\t-\t2\t4000\t06:00:00\t30")
    report = {"tests/calibration/test_novel_orfs_cal.py::test_planted_tir_calibration_300k_full"
              "[a]": True,
              "tests/calibration/test_novel_orfs_cal.py::test_gc_full": False}  # fmt: skip
    run = box.run("novel_orfs", report, "0-19")
    assert run.rc == 2, run.stderr
    _no_lane_started(run)
    assert record(run.log, 11)["exit"] == 2


@pytest.mark.parametrize("case", ["collect-error", "empty"])
def test_collection_error_exits_2_with_a_record(box: Box, case: str) -> None:
    report, durations = det_report()
    box.durations(durations)
    extra = {"STUB_RC_collect": "2"} if case == "collect-error" else {}
    run = box.run("determinants", report if case != "empty" else {}, "0-19", extra)
    assert run.rc == 2, run.stderr
    _no_lane_started(run)
    assert record(run.log, 10)["exit"] == 2


def test_partition_problem_exits_2_before_any_lane(box: Box) -> None:
    report, durations = det_report()
    box.durations(durations)
    hide = f"{DET}[full-gc-lightgbm]"
    assert hide in report
    run = box.run("determinants", report, "0-19", {"STUB_COLLECT_HIDE": hide})
    assert run.rc == 2, run.stderr
    _no_lane_started(run)
    rec = record(run.log, 10)
    assert rec["exit"] == 2 and rec["checks"], rec


def test_unknown_suite_starts_nothing(box: Box) -> None:
    report, _ = det_report()
    run = box.run("no_such_suite", report, "0-19")
    assert run.rc == 2, run.stderr
    _no_lane_started(run)
    assert run.forbidden() == []


# ------------------------------------------------------------------ counts and confirm parts


def test_short_junit_is_a_count_problem(box: Box, lanes_mod: ModuleType) -> None:
    report, durations = det_report()
    box.durations(durations)
    want = expected_lanes(lanes_mod, box, "determinants", report)
    run = box.run("determinants", report, "0-19", {"STUB_DROP_group2": "1"})
    assert run.rc == 2, run.stderr
    rec = record(run.log, 10)
    assert rec["exit"] == 2 and rec["checks"], rec
    (lane,) = [e for e in rec["lanes"] if e["index"] == 2]
    assert lane["exit"] == 0
    assert lane["junit_tests"] == len(want[2].expected) - 1


def test_missing_junit_is_recorded_and_fails_the_rows(box: Box, rows_mod: ModuleType) -> None:
    report, durations = det_report()
    box.durations(durations)
    run = box.run("determinants", report, "0-19", {"STUB_NO_JUNIT_group3": "1"})
    assert run.rc == 2, run.stderr
    rec = record(run.log, 10)
    assert "task-10.group3.junit.xml" in rec["parts"], "parts are the parts launched"
    assert not (run.log / "task-10.group3.junit.xml").exists()
    assert rec["exit"] == 2 and rec["checks"]
    res, rows = nightly(rows_mod, box, 10, "determinants")
    assert "junit missing" in res.problems and res.fail_closed, res
    assert [r.split("\t")[3] for r in rows] == ["1", "1"]


@pytest.mark.parametrize("merge_rc", [0, 1])
def test_missing_confirm_parts_is_a_check_problem(box: Box, merge_rc: int) -> None:
    run = box.run("variants", var_report(), "0-19",
                  {"STUB_NO_CONFIRM_PARTS": "1", "STUB_RC_merge": str(merge_rc)})  # fmt: skip
    assert run.of(proc="criteria", role="merge"), "variants has shards: the merge runs"
    assert run.rc == 2, run.stderr
    rec = record(run.log, 13)
    assert rec["exit"] == 2 and rec["checks"], rec


def test_confirm_parts_empty_list_adds_no_parts(box: Box, lanes_mod: ModuleType) -> None:
    want = expected_lanes(lanes_mod, box, "variants", var_report())
    run = box.run("variants", var_report(), "0-19")
    assert run.rc == 0, run.stderr
    (merge,) = run.of(proc="criteria", role="merge")
    assert merge["env"]["GIPPA_CALIBRATION_CONFIRM_LANES"] == "2"
    rec = record(run.log, 13)
    assert rec["parts"] == lanes_mod.record_parts(want, [])
    assert rec["parts"] == ["task-13.shard0.junit.xml", "task-13.shard1.junit.xml",
                            "task-13.group0.junit.xml"]  # fmt: skip
    assert rec["sharded"] is True and rec["checks"] == []


def test_deselect_list_survives_shell_quoting(box: Box) -> None:
    run = box.run("variants", var_report(), "0-19")
    assert run.rc == 0, run.stderr
    for k in range(2):
        (w,) = [w for w in run.lane_wrappers() if w["role"] == f"shard{k}"]
        assert set(w["deselect"]) == set(VAR_DIRECT), w["deselect"]
        assert junit_count(run.log / f"task-13.shard{k}.junit.xml") == len(VAR_HARNESS)
    assert junit_count(run.log / "task-13.group0.junit.xml") == len(VAR_DIRECT)


# ------------------------------------------------------------------ exit folding (example 5)


FOLD_CASES = [
    ({}, 0),
    ({"STUB_RC_group0": "1"}, 1),
    ({"STUB_RC_merge": "1"}, 1),
    ({"STUB_RC_shard0": "1", "STUB_RC_merge": "1"}, 1),
    ({"STUB_RC_merge": "2"}, 2),
    ({"STUB_RC_shard1": "3", "STUB_RC_merge": "2"}, 3),
    ({"STUB_RC_shard1": "3", "STUB_RC_group0": "4"}, 3),
    ({"STUB_RC_group0": "4", "STUB_RC_merge": "2"}, 4),
    ({"STUB_RC_group0": "1", "STUB_RC_merge": "1", "STUB_DROP_shard0": "1"}, 2),
]


@pytest.mark.parametrize(("extra", "want"), FOLD_CASES,
                         ids=[",".join(f"{k[5:]}={v}" for k, v in e.items()) or "clean"
                              for e, _ in FOLD_CASES])  # fmt: skip
def test_exit_folding(box: Box, lanes_mod: ModuleType, rows_mod: ModuleType,
                      extra: dict[str, str], want: int) -> None:  # fmt: skip
    run = box.run("variants", var_report(), "0-19", extra)
    assert run.rc == want, (run.rc, run.stderr)
    rec = record(run.log, 13)
    assert rec["exit"] == want
    lane_rc = {f"{e['kind']}{e['index']}": e["exit"] for e in rec["lanes"]}
    for role in ("shard0", "shard1", "group0"):
        assert lane_rc[role] == int(extra.get(f"STUB_RC_{role}", "0")), lane_rc
    merge_rc = int(extra.get("STUB_RC_merge", "0"))
    assert lanes_mod.fold_exit([lane_rc[r] for r in ("shard0", "shard1", "group0")], merge_rc) == (
        want if "STUB_DROP_shard0" not in extra else 1
    )
    res, rows = nightly(rows_mod, box, 13, "variants")
    assert res.exit == want
    assert rows[0].split("\t")[3] == ("0" if want == 0 else "1"), rows[0]


def test_r4_strictness_group_lane_non_attack_failure_fails_both_rows(
    box: Box, rows_mod: ModuleType
) -> None:
    bad = f"{VAR}::test_direct_a"
    run = box.run("variants", var_report(), "0-19", {"STUB_FAIL_group0": bad})
    assert run.rc == 1, run.stderr
    rec = record(run.log, 13)
    assert rec["exit"] == 1 and rec["sharded"] is True
    res, rows = nightly(rows_mod, box, 13, "variants")
    assert res.attacks.total > 0 and res.attacks.bad == 0, res
    assert res.full_failed and res.attacks_failed, res
    assert [r.split("\t")[3] for r in rows] == ["1", "1"], rows


# ------------------------------------------------------------------ signals (example 10)


def _wait_for(box: Box, pred: Any, timeout: float = 60) -> list[dict[str, Any]]:
    end = time.monotonic() + timeout
    while time.monotonic() < end:
        ev = box.events()
        if pred(ev):
            return ev
        time.sleep(0.1)
    pytest.fail(f"timed out waiting; events: {[(e.get('event'), e.get('role')) for e in ev]}")


@pytest.fixture
def decoy(box: Box) -> Iterator[subprocess.Popen[bytes]]:
    """A ``sleep`` named pytest outside the job's process tree."""
    d = box.base / "decoy"
    d.mkdir()
    sleep = shutil.which("sleep")
    assert sleep is not None
    (d / "pytest").symlink_to(sleep)
    proc = subprocess.Popen([str(d / "pytest"), "300"])
    try:
        yield proc
    finally:
        proc.kill()
        proc.wait()


def _signal_run(box: Box, suite: str, report: Mapping[str, bool], extra: Mapping[str, str],
                ready: Any, sig: int = signal.SIGUSR1,
                timeout: float = 60) -> tuple[Run, float]:  # fmt: skip
    code = _code()  # a pkill-style implementation would hit this test session: refuse first
    assert not re.search(r"\b(pkill|killall)\b", code)
    t0 = time.monotonic()
    proc, out, err = box.start(box.env(suite, report, "0-19", extra))
    _wait_for(box, ready)
    sent = time.time()
    proc.send_signal(sig)
    return box.finish(proc, out, err, t0, timeout=timeout), sent


def _inner_started(roles: set[str]) -> Any:
    def pred(ev: list[dict[str, Any]]) -> bool:
        got = {e["role"] for e in ev if e.get("event") == "start" and e.get("proc") == "pytest"}
        return roles <= got

    return pred


# walks the job's process tree through Linux ``/proc`` (FarmShare only, spec G19)
@pytest.mark.gnu_userland
def test_usr1_during_lanes_interrupts_pytest_and_records_143(
    box: Box, rows_mod: ModuleType, decoy: subprocess.Popen[bytes]
) -> None:
    report, durations = det_report()
    box.durations(durations)
    roles = {f"group{g}" for g in range(5)}
    extra = {f"STUB_SLEEP_{r}": "120" for r in roles}
    run, sent = _signal_run(box, "determinants", report, extra, _inner_started(roles))
    ints = [e for e in run.events if e.get("event") == "int"]
    assert {e["role"] for e in ints if e["proc"] == "pytest"} == roles, ints
    assert not [e for e in ints if e["proc"] == "criteria"], "signal pytest, not the wrapper"
    assert max(e["t"] for e in ints) - sent <= 3.0, "USR1 handled at once (& wait)"
    assert run.rc != 0
    rec = record(run.log, 10)
    assert rec["exit"] == 143, rec
    for r in roles:
        assert junit_count(run.log / f"task-10.{r}.junit.xml") == 1, "partial JUnit"
    assert decoy.poll() is None, "a pytest outside the job's tree was signalled"
    res, rows = nightly(rows_mod, box, 10, "determinants")
    assert res.fail_closed and "exit 143" in res.problems, res
    assert "determinants exit 143" in rows[0] and rows[0].split("\t")[3] == "1", rows[0]


# walks the job's process tree through Linux ``/proc`` (FarmShare only, spec G19)
@pytest.mark.gnu_userland
def test_usr1_during_merge_interrupts_block2_pytests(
    box: Box, rows_mod: ModuleType, decoy: subprocess.Popen[bytes]
) -> None:
    roles = {"lane0", "lane1"}
    extra = {"STUB_PENDING": "1", **{f"STUB_SLEEP_{r}": "120" for r in roles}}
    run, sent = _signal_run(box, "variants", var_report(), extra, _inner_started(roles))
    ints = [e for e in run.events if e.get("event") == "int" and e["proc"] == "pytest"]
    assert {e["role"] for e in ints} == roles, ints
    assert max(e["t"] for e in ints) - sent <= 3.0
    rec = record(run.log, 13)
    assert rec["exit"] == 143, rec
    assert {"task-13.merge.lane0.junit.xml", "task-13.merge.lane1.junit.xml"} <= set(rec["parts"])
    for k in range(2):
        assert junit_count(run.log / f"task-13.merge.lane{k}.junit.xml") == 1
    assert decoy.poll() is None
    res, _rows = nightly(rows_mod, box, 13, "variants")
    assert res.fail_closed and "exit 143" in res.problems, res


# walks the job's process tree through Linux ``/proc`` (FarmShare only, spec G19)
@pytest.mark.gnu_userland
def test_usr1_during_collection_records_143(box: Box) -> None:
    report, durations = det_report()
    box.durations(durations)

    def ready(ev: list[dict[str, Any]]) -> bool:
        return any(e.get("event") == "start" and e.get("role") == "collect" for e in ev)

    run, sent = _signal_run(box, "determinants", report, {"STUB_SLEEP_collect": "120"}, ready)
    ints = [e for e in run.events if e.get("event") == "int" and e["proc"] == "pytest"]
    assert ints and min(e["t"] for e in ints) - sent <= 3.0, ints
    _no_lane_started(run)
    assert record(run.log, 10)["exit"] == 143


# walks the job's process tree through Linux ``/proc`` (FarmShare only, spec G19)
@pytest.mark.gnu_userland
def test_pytest_ignoring_sigint_still_gets_a_record_after_the_grace(box: Box) -> None:
    report, durations = det_report()
    box.durations(durations)
    roles = {f"group{g}" for g in range(5)}
    extra = {"STUB_IGNORE_INT": "1", "GIPPA_CALIBRATION_GRACE": "2",
             **{f"STUB_SLEEP_{r}": "120" for r in roles}}  # fmt: skip
    run, _sent = _signal_run(box, "determinants", report, extra, _inner_started(roles),
                             timeout=30)  # fmt: skip
    assert run.elapsed < 60
    rec = record(run.log, 10)
    assert rec["exit"] == 143, rec
    assert run.rc != 0


# walks the job's process tree through Linux ``/proc`` (FarmShare only, spec G19)
@pytest.mark.gnu_userland
def test_sigterm_during_lanes_records_143_once(box: Box) -> None:
    report, durations = det_report()
    box.durations(durations)
    roles = {f"group{g}" for g in range(5)}
    extra = {f"STUB_SLEEP_{r}": "120" for r in roles}
    run, _sent = _signal_run(box, "determinants", report, extra, _inner_started(roles),
                             sig=signal.SIGTERM)  # fmt: skip
    assert run.rc != 0
    assert record(run.log, 10)["exit"] == 143


# ------------------------------------------------------------------ review round (#389)


def test_stale_confirm_parts_is_removed_before_the_merge(box: Box) -> None:
    box.log.mkdir(parents=True)
    (box.log / "variants.confirm_parts.json").write_text('["task-13.bogus.junit.xml"]\n')
    run = box.run("variants", var_report(), "0-19", {"STUB_NO_CONFIRM_PARTS": "1"})
    assert run.rc == 2, run.stderr
    rec = record(run.log, 13)
    assert rec["exit"] == 2 and rec["checks"], rec
    assert "task-13.bogus.junit.xml" not in rec["parts"]


def test_leaked_confirm_and_replicate_variables_are_unset(box: Box) -> None:
    leak = {"GIPPA_CALIBRATION_CONFIRM": '[["x", "c", 0.05, 10]]',
            "GIPPA_CALIBRATION_CONFIRM_LANES": "5", "GIPPA_CALIBRATION_CONFIRM_JUNIT": "/nowhere",
            "GIPPA_CALIBRATION_CONFIRM_SLICE": "0/2",
            "GIPPA_CALIBRATION_REPLICATES": "0:1"}  # fmt: skip
    run = box.run("variants", var_report(), "0-19", leak)
    assert run.rc == 0, run.stderr
    for e in run.lane_wrappers() + [e for e in run.of(proc="pytest") if e["role"] != "collect"
                                    and not e["role"].startswith("lane")]:  # fmt: skip
        assert not any(e["env"].get(k) for k in leak), (e["role"], e["env"])
    (merge,) = run.of(proc="criteria", role="merge")
    assert merge["env"]["GIPPA_CALIBRATION_CONFIRM_LANES"] == "2"
    assert merge["env"]["GIPPA_CALIBRATION_CONFIRM_JUNIT"] == f"{run.log}/task-13.merge"
    for k in ("GIPPA_CALIBRATION_CONFIRM", "GIPPA_CALIBRATION_CONFIRM_SLICE",
              "GIPPA_CALIBRATION_REPLICATES"):  # fmt: skip
        assert not merge["env"].get(k), (k, merge["env"])


def test_record_is_written_when_the_record_helper_dies(box: Box) -> None:
    """S1: the traps must still write ``task-<i>.json`` when the first attempt is killed."""
    run = box.run("variants", var_report(), "0-19", {"STUB_DIE_ONCE": "record"})
    assert Path(f"{box.stub_log}.died").exists(), "no interpreter call with 'record' in argv"
    rec = record(run.log, 13)
    assert rec["exit"] == run.rc, (rec["exit"], run.rc)
    assert rec["exit"] in (0, 2), rec


def _job_helper(box: Box) -> Path:
    m = re.search(r"<<'PY'\n(.*?)\nPY\n", SBATCH.read_text(), re.S)
    assert m, "the job's Python helper here-document (<<'PY' ... PY)"
    path = box.base / "suite_job.py"
    path.write_text(m.group(1) + "\n")
    return path


def _hand_state(box: Box, lanes_mod: ModuleType, *, lane_exits: Mapping[int, str | None],
                merge_exit: str | None, merged: bool = True) -> Path:  # fmt: skip
    """A finished variants job's state: every lane's JUnit present and complete."""
    import dataclasses

    report = var_report()
    rows = {r.suite: r for r in lanes_mod.read_plan(
        box.repo / "scripts" / "slurm" / "calibration_plan.tsv", SUITES)}  # fmt: skip
    lanes = lanes_mod.lanes(rows["variants"], report, {})
    st = box.base / "state"
    st.mkdir()
    box.log.mkdir(parents=True, exist_ok=True)
    (st / "meta.json").write_text(
        json.dumps(
            {
                "suite": "variants",
                "task": 13,
                "job": JOB,
                "checkout": str(box.repo),
                "log": str(box.log),
            }
        )
    )
    (st / "plan.json").write_text(json.dumps(dataclasses.asdict(rows["variants"])))
    (st / "lanes.json").write_text(json.dumps([dataclasses.asdict(ln) for ln in lanes]))
    (st / "launched").write_text("")
    for pos, ln in enumerate(lanes):
        cases = "".join(f'<testcase classname="c" name="t{j}" time="0.1"/>'
                        for j in range(len(ln.expected)))  # fmt: skip
        (box.log / ln.part).write_text(f"<testsuites><testsuite>{cases}</testsuite></testsuites>")
        rc = lane_exits.get(pos, "0")
        if rc is not None:
            (st / f"exit.{pos}").write_text(f"{rc}\n")
    if merged:
        (st / "merge.launched").write_text("")
        (box.log / "variants.confirm_parts.json").write_text("[]\n")
        if merge_exit is not None:
            (st / "exit.merge").write_text(merge_exit)
    return st


@pytest.mark.parametrize(
    ("lane_exits", "merge_exit", "want"),
    [
        ({}, "0\n", 0),
        ({}, None, 2),  # S2: launched merge, exit lost
        ({}, "garbage\n", 2),  # S2: unparsable merge exit
        ({}, "", 2),
        ({1: None}, "0\n", 2),  # a lane with no exit file
        ({2: None}, "1\n", 2),
    ],
    ids=[
        "clean",
        "merge-exit-missing",
        "merge-exit-garbage",
        "merge-exit-empty",
        "lane-exit-missing",
        "group-exit-missing",
    ],
)
def test_fold_of_a_lost_exit_is_2(box: Box, lanes_mod: ModuleType,
                                  lane_exits: dict[int, str | None], merge_exit: str | None,
                                  want: int) -> None:  # fmt: skip
    helper = _job_helper(box)
    st = _hand_state(box, lanes_mod, lane_exits=lane_exits, merge_exit=merge_exit)
    res = subprocess.run([sys.executable, str(helper), "record", str(st), "-"],
                         capture_output=True, text=True, check=False, timeout=60)  # fmt: skip
    assert res.returncode == 0, res.stderr
    assert res.stdout.strip() == str(want), (res.stdout, res.stderr)
    assert record(box.log, 13)["exit"] == want


def _alive(pid: int) -> bool:
    try:
        stat = Path(f"/proc/{pid}/stat").read_text()
    except OSError:
        return False
    return stat.rsplit(")", 1)[1].split()[0] != "Z"


# walks the job's process tree through Linux ``/proc`` (FarmShare only, spec G19)
@pytest.mark.gnu_userland
def test_usr1_keeps_the_grace_then_kills_survivors_in_its_tree(
    box: Box, decoy: subprocess.Popen[bytes]
) -> None:
    report, durations = det_report()
    box.durations(durations)
    roles = {f"group{g}" for g in range(5)}
    grace = 4
    extra = {"STUB_IGNORE_INT": "1", "GIPPA_CALIBRATION_GRACE": str(grace),
             **{f"STUB_SLEEP_{r}": "120" for r in roles}}  # fmt: skip
    run, sent = _signal_run(box, "determinants", report, extra, _inner_started(roles),
                            timeout=60)  # fmt: skip
    done = time.time()
    assert done - sent >= grace - 0.5, "USR1 waits the grace for pytests that ignore SIGINT"
    assert record(run.log, 10)["exit"] == 143
    lane_pids = [e["pid"] for e in run.of(proc="pytest") if e["role"] in roles]
    assert len(lane_pids) == 5
    deadline = time.monotonic() + 3
    while any(_alive(p) for p in lane_pids) and time.monotonic() < deadline:
        time.sleep(0.1)
    assert not [p for p in lane_pids if _alive(p)], "survivors after the grace get SIGKILL"
    assert decoy.poll() is None, "only the job's own tree is killed"


def test_sigterm_writes_the_record_without_waiting_the_grace(box: Box) -> None:
    """S3: KillWait is 30 s on FarmShare; a TERM record must not wait a 120 s grace."""
    report, durations = det_report()
    box.durations(durations)
    roles = {f"group{g}" for g in range(5)}
    extra = {"STUB_IGNORE_INT": "1", "GIPPA_CALIBRATION_GRACE": "60",
             **{f"STUB_SLEEP_{r}": "120" for r in roles}}  # fmt: skip
    proc, _out, err = box.start(box.env("determinants", report, "0-19", extra))
    _wait_for(box, _inner_started(roles))
    proc.send_signal(signal.SIGTERM)
    sent = time.monotonic()
    rec_path = box.log / "task-10.json"
    while not rec_path.exists() and time.monotonic() - sent < 20:
        time.sleep(0.1)
    took = time.monotonic() - sent
    assert rec_path.exists(), (
        f"no record {took:.0f} s after SIGTERM (grace 60 s): {err.read_text()}"
    )
    time.sleep(0.3)  # the rename is atomic; let a second (wrong) write show up
    assert record(box.log, 10)["exit"] == 143


# walks the job's process tree through Linux ``/proc`` (FarmShare only, spec G19)
@pytest.mark.gnu_userland
def test_interrupt_rescans_for_pytests_started_after_usr1(box: Box) -> None:
    report, durations = det_report()
    box.durations(durations)
    early = {f"group{g}" for g in range(4)}
    extra = {"GIPPA_CALIBRATION_GRACE": "15", "STUB_SPAWN_DELAY_group4": "2.5",
             **{f"STUB_SLEEP_group{g}": "120" for g in range(5)}}  # fmt: skip

    def ready(ev: list[dict[str, Any]]) -> bool:
        wrappers = {e["role"] for e in ev if e.get("event") == "start"
                    and e.get("proc") == "criteria"}  # fmt: skip
        return _inner_started(early)(ev) and "group4" in wrappers

    run, sent = _signal_run(box, "determinants", report, extra, ready, timeout=60)
    (late,) = run.of(proc="pytest", role="group4")
    assert late["t"] > sent, "the group4 pytest must start after USR1 for this test"
    ints = {e["role"] for e in run.events if e.get("event") == "int" and e["proc"] == "pytest"}
    assert ints == early | {"group4"}, ints
    assert record(run.log, 10)["exit"] == 143


# ------------------------------------------------------------------ the stub itself


def test_stub_pytest_models_prefix_deselection(box: Box) -> None:
    """Guard on the fake: ``--deselect`` is a prefix match, quoted items survive ``shlex``."""
    box.collection.write_text(json.dumps({"t.py::a": False, "t.py::a_full": False,
                                          "t.py::b[x y]": True}))  # fmt: skip
    env = {
        **os.environ,
        "STUB_LOG": str(box.stub_log),
        "STUB_COLLECTION": str(box.collection),
        "PYTEST_ADDOPTS": f"--deselect {shlex.quote('t.py::b[x y]')} --deselect=t.py::a",
    }
    res = subprocess.run([str(box.venv / "bin" / "python"), "-m", "pytest", "--collect-only",
                          "-q"], env=env, capture_output=True, text=True, check=False)  # fmt: skip
    assert res.returncode == 5 and res.stdout.split("\n\n")[0].strip() == "", res.stdout
    env["PYTEST_ADDOPTS"] = "--deselect=t.py::b"
    res = subprocess.run([str(box.venv / "bin" / "python"), "-m", "pytest", "--collect-only",
                          "-q"], env=env, capture_output=True, text=True, check=False)  # fmt: skip
    assert res.returncode == 0
    assert res.stdout.split("\n\n")[0].split() == ["t.py::a", "t.py::a_full"]
    res = subprocess.run([str(box.venv / "bin" / "python"), "-c", "print(6 * 7)"], env=env,
                         capture_output=True, text=True, check=False)  # fmt: skip
    assert res.stdout.strip() == "42"
