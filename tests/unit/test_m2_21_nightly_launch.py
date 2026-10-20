"""Nightly scheduling and the R14 worst case (M2-21 lead decisions 12, 14 and 18).

``scripts/slurm/nightly_launch.sh`` runs in the sandbox of ``tests/helpers/slurm_sandbox.py``
(fake ``sbatch``/``sacct``/``gh`` first on ``PATH``, sandbox ``GIAP_REPO``/``GIAP_ROOT``); the
scrontab and the sbatch scripts are read as text. Nothing is submitted or installed.

Conventions settled by the test author:

* Decision 18: the launcher submits ``nightly.sbatch`` with ``--array=0-3`` (``--array 0-3`` and
  ``-a 0-3`` accepted) and then one summary job with ``--dependency`` containing
  ``afterany:<array job ID>``; the summary job asks for 1 CPU (``--cpus-per-task``/``-c`` absent
  or 1, also in its script's ``#SBATCH`` lines), is named ``giap-*`` and runs
  ``scripts/ci/post_summary.py`` (``--wrap`` or a script). The launcher exits 0; when the array
  submission fails it exits non-zero and submits no summary job.
* Decision 14: the summary job, run in the sandbox after one task of the array was killed (no row;
  ``sacct`` says ``TIMEOUT``), posts exactly one comment for the array job; it is not green and has
  a table row for the missing task. When the array already posted itself it posts nothing.
* Decision 18: the scrontab has one cron line running ``scripts/slurm/nightly_launch.sh``; its
  ``#SCRON`` options carry no ``--array``, a ``giap-*`` job name and 1 CPU.
* Decision 12: worst case = array tasks x ``--cpus-per-task`` + Nextflow child jobs + the 1-CPU
  summary job <= 128 (R14). Child jobs exist only when the default ``*_PROFILE`` uses the Slurm
  executor (``farmshare``, not ``farmshare_local``): then each task that runs ``giap run`` adds
  ``queueSize`` x ``min(16, cpus-per-task)`` CPUs, ``queueSize`` read from the script if it sets
  one, else from ``workflow/conf/farmshare.config``.
"""

from __future__ import annotations

import re
import shlex
import subprocess
from pathlib import Path
from typing import Any

import pytest
from helpers.slurm_sandbox import (
    REPO,
    Sandbox,
    assert_honours_giap_repo,
    export_pairs,
    make_sandbox,  # fmt: skip
    option,
)

SLURM = REPO / "scripts" / "slurm"
LAUNCH = SLURM / "nightly_launch.sh"
NIGHTLY = SLURM / "nightly.sbatch"
WEEKLY = SLURM / "weekly.sbatch"
SCRONTAB = SLURM / "nightly.scrontab"
ARRAY_ID = "4242"  # the fake sbatch numbers its jobs 4242, 4243, …


def _launcher_text() -> str:
    assert LAUNCH.is_file(), "scripts/slurm/nightly_launch.sh does not exist (M2-21 decision 18)"
    return LAUNCH.read_text()


@pytest.fixture
def sb(tmp_path: Path) -> Sandbox:
    assert_honours_giap_repo(_launcher_text(), "nightly_launch.sh")
    return make_sandbox(tmp_path)


def _launch(sb: Sandbox, **env: str):  # type: ignore[no-untyped-def]
    return sb.run(["bash", str(sb.repo / "scripts" / "slurm" / "nightly_launch.sh")], env=env)


def _is_array(call: dict[str, Any]) -> bool:
    return option(call["argv"], "array", "a") == "0-3" and any(
        a.endswith("nightly.sbatch") for a in call["argv"]
    )


def _summary_calls(calls: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [c for c in calls if "afterany" in (option(c["argv"], "dependency", "d") or "")]


def _sbatch_directive(script: str | None, long: str, short: str) -> str | None:
    if not script:
        return None
    m = re.search(rf"^#SBATCH\s+(?:--{long}[= ]|-{short}\s*)(\S+)", script, re.M)
    return m.group(1) if m else None


def test_launcher_parses() -> None:
    _launcher_text()
    res = subprocess.run(["bash", "-n", str(LAUNCH)], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr


def test_launcher_submits_the_array_then_an_afterany_summary(sb: Sandbox) -> None:
    res = _launch(sb)
    assert res.returncode == 0, f"{res.stdout}\n{res.stderr}"
    calls = sb.sbatch_calls()
    arrays = [c for c in calls if _is_array(c)]
    assert len(arrays) == 1, f"expected one `sbatch --array=0-3 … nightly.sbatch`: {calls}"
    array_id = str(arrays[0]["id"])
    summaries = _summary_calls(calls)
    assert len(summaries) == 1, f"expected one afterany summary job: {calls}"
    s = summaries[0]
    dep = option(s["argv"], "dependency", "d") or ""
    assert re.search(rf"afterany:(?:[0-9_:]*:)?{array_id}\b", dep), (dep, array_id)
    assert calls.index(s) > calls.index(arrays[0])
    assert option(s["argv"], "array", "a") is None
    cpus = option(s["argv"], "cpus-per-task", "c") or _sbatch_directive(
        s["script"], "cpus-per-task", "c"
    )
    assert cpus in (None, "1"), f"the summary job takes 1 CPU, not {cpus}"
    name = option(s["argv"], "job-name", "J") or _sbatch_directive(s["script"], "job-name", "J")
    assert name and name.startswith("giap-"), f"summary job name {name!r} (R14: giap-*)"
    body = (option(s["argv"], "wrap") or "") + (s["script"] or "")
    assert "post_summary.py" in body, f"the summary job does not run post_summary.py: {s}"


def test_launcher_fails_without_a_summary_when_the_array_is_refused(sb: Sandbox) -> None:
    res = _launch(sb, FAKE_SBATCH_FAIL="0")
    assert res.returncode != 0, "the launcher hid a failed array submission"
    later = [c for c in sb.sbatch_calls() if not c["failed"]]
    assert not _summary_calls(later), f"summary submitted for a job that does not exist: {later}"


def _summary_command(sb: Sandbox, call: dict[str, Any]) -> tuple[list[str], dict[str, str]]:
    env = {k: v for k, v in export_pairs(call["argv"]).items() if k != "ALL"}
    env.update({"SLURM_JOB_ID": str(call["id"]), "SLURM_JOB_NAME": "giap-nightly-summary",
                "SLURM_CPUS_PER_TASK": "1"})  # fmt: skip
    wrap = option(call["argv"], "wrap")
    if wrap is not None:
        assert_honours_giap_repo(wrap, "summary --wrap")
        return ["bash", "-c", wrap], env
    argv = call["argv"]
    cwd = Path(call["cwd"])
    idx = next(i for i, a in enumerate(argv) if not a.startswith("-") and (cwd / a).is_file())
    script = cwd / argv[idx]
    assert_honours_giap_repo(script.read_text(), str(script))
    return ["bash", str(script), *argv[idx + 1 :]], env


def _array_with_rows(sb: Sandbox, tasks: list[int]) -> dict[str, Any]:
    res = _launch(sb)
    assert res.returncode == 0, res.stderr
    (summary,) = _summary_calls(sb.sbatch_calls())
    commit = sb.head[:7]
    tsv = sb.root / "logs" / "nightly.tsv"
    tsv.write_text("".join(
        f"2026-10-20\t{ARRAY_ID}_{t}\t{commit}\t0\ttask {t}: 1 passed\n" for t in tasks
    ))  # fmt: skip
    states = {t: "COMPLETED" for t in range(4)} | {t: "TIMEOUT" for t in range(4) if t not in tasks}
    sb.set_sacct([
        {"JobID": f"{ARRAY_ID}_{t}", "JobName": "giap-nightly", "State": st,
         "ExitCode": "0:0" if st == "COMPLETED" else "0:15", "NodeList": "barley-01",
         "Elapsed": "06:00:00" if st == "TIMEOUT" else "00:10:00"}
        for t, st in states.items()
    ])  # fmt: skip
    return summary


def _comments(sb: Sandbox) -> list[dict[str, Any]]:
    return [c for c in sb.gh_calls() if c["argv"][:2] == ["issue", "comment"]]


def test_summary_job_posts_when_a_task_was_killed(sb: Sandbox) -> None:
    """Decision 14: task 2 hit its time limit and wrote nothing; the night is still posted."""
    summary = _array_with_rows(sb, [0, 1, 3])
    cmd, env = _summary_command(sb, summary)
    res = sb.run(cmd, env=env, cwd=Path(summary["cwd"]) if summary["cwd"].startswith(
        str(sb.tmp)) else sb.repo)  # fmt: skip
    comments = _comments(sb)
    assert len(comments) == 1, (
        f"expected one issue comment, got {comments}\n{shlex.join(cmd)}\n{res.stdout}\n{res.stderr}"
    )
    body = comments[0]["body"] or ""
    assert ARRAY_ID in body, body
    assert "**green**" not in body and "green" not in body.split("\n", 1)[0].lower(), body
    assert re.search(r"^\|\s*2\s*\|", body, re.M), f"no row for the killed task 2:\n{body}"


def test_summary_job_does_not_post_twice(sb: Sandbox) -> None:
    summary = _array_with_rows(sb, [0, 1, 2, 3])
    (sb.root / "logs" / "nightly.tsv.posted").write_text(f"{ARRAY_ID}\n")
    cmd, env = _summary_command(sb, summary)
    res = sb.run(cmd, env=env, cwd=sb.repo)
    assert not _comments(sb), f"posted {ARRAY_ID} a second time:\n{res.stdout}\n{res.stderr}"


def _scron() -> tuple[list[str], list[str]]:
    lines = SCRONTAB.read_text().splitlines()
    opts = [ln for ln in lines if ln.startswith("#SCRON")]
    cron = [ln for ln in lines if ln.strip() and not ln.lstrip().startswith("#")]
    return opts, cron


def test_scrontab_runs_the_one_cpu_launcher() -> None:
    opts, cron = _scron()
    assert len(cron) == 1, cron
    assert cron[0].split()[-1].endswith("scripts/slurm/nightly_launch.sh"), cron
    joined = " ".join(opts)
    assert "--array" not in joined and not re.search(r"(^|\s)-a\s", joined), opts
    name = re.search(r"--job-name[= ](\S+)|-J\s*(\S+)", joined)
    assert name and (name.group(1) or name.group(2)).startswith("giap-"), opts
    cpus = re.search(r"--cpus-per-task[= ](\d+)|-c\s*(\d+)", joined)
    assert cpus is None or (cpus.group(1) or cpus.group(2)) == "1", opts


# --- decision 12: R14 worst case -----------------------------------------------------------------


def _code(path: Path) -> str:
    return "\n".join(ln for ln in path.read_text().splitlines()
                     if not ln.lstrip().startswith("#") or ln.startswith("#SBATCH"))  # fmt: skip


def _directive(path: Path, key: str) -> str | None:
    m = re.search(rf"^#SBATCH\s+--{key}[= ](\S+)", path.read_text(), re.M)
    return m.group(1) if m else None


def _array_size(spec: str) -> int:
    ids: set[int] = set()
    for part in spec.split("%")[0].split(","):
        lo, _, hi = part.partition("-")
        ids |= set(range(int(lo), int(hi or lo) + 1))
    if "%" in spec:
        return min(len(ids), int(spec.split("%")[1]))
    return len(ids)


def _default_profile(path: Path, var: str) -> str:
    m = re.search(rf"\$\{{{var}:-([^}}]*)\}}", _code(path))
    assert m, f"{path.name}: no ${{{var}:-…}} default"
    return m.group(1)


def _queue_size(path: Path) -> int:
    m = re.search(r"queueSize\D{0,6}(\d+)", _code(path))
    if m:
        return int(m.group(1))
    conf = (REPO / "workflow" / "conf" / "farmshare.config").read_text()
    m = re.search(r"queueSize\s*=\s*(\d+)", conf)
    assert m
    return int(m.group(1))


def _slurm_children(profile: str) -> bool:
    return "farmshare" in profile.split(",")


def _worst_case(path: Path, var: str, heads: int) -> tuple[int, str]:
    cpus = int(_directive(path, "cpus-per-task") or 1)
    array = _directive(path, "array") or "0"
    if path == NIGHTLY and LAUNCH.is_file():
        m = re.search(r"--array[= ](\S+)", LAUNCH.read_text())
        array = m.group(1).strip("\"'") if m else array
    tasks = _array_size(array)
    profile = _default_profile(path, var)
    children = heads * _queue_size(path) * min(16, cpus) if _slurm_children(profile) else 0
    total = tasks * cpus + children + 1
    return total, f"{tasks} x {cpus} + children {children} + summary 1 (profile {profile})"


def _nightly_nextflow_heads() -> int:
    """Array arms of nightly.sbatch whose text (with the functions they call) uses $PROFILE."""
    code = _code(NIGHTLY)
    funcs = {m.group(1): m.end() for m in re.finditer(r"^(\w+)\(\)\s*\{", code, re.M)}
    bodies = {}
    for name, start in funcs.items():
        end = re.search(r"^\}", code[start:], re.M)
        bodies[name] = code[start : start + (end.start() if end else len(code))]
    m = re.search(r'case\s+"?\$\{?SLURM_ARRAY_TASK_ID[^\n]*\bin\b(.*?)\besac\b', code, re.S)
    assert m
    heads = 0
    for arm in re.findall(r"^\s*\d+\)(.*?);;", m.group(1), re.S | re.M):
        text = arm + "".join(b for n, b in bodies.items() if re.search(rf"\b{n}\b", arm))
        heads += bool(re.search(r"\$\{?PROFILE\b", text))
    return max(heads, 1)


def test_nightly_worst_case_cpus_within_r14() -> None:
    total, how = _worst_case(NIGHTLY, "NIGHTLY_PROFILE", _nightly_nextflow_heads())
    assert total <= 128, f"nightly worst case {total} CPUs > 128 (R14): {how}"


def test_weekly_worst_case_cpus_within_r14() -> None:
    tasks = _array_size(_directive(WEEKLY, "array") or "0")
    total, how = _worst_case(WEEKLY, "WEEKLY_PROFILE", tasks)
    assert total <= 128, f"weekly worst case {total} CPUs > 128 (R14): {how}"


@pytest.mark.parametrize(
    ("path", "var"), [(NIGHTLY, "NIGHTLY_PROFILE"), (WEEKLY, "WEEKLY_PROFILE")]
)
def test_default_profile_is_local_or_queue_capped(path: Path, var: str) -> None:
    """Decision 12: default ``farmshare_local,apptainer``, or ``farmshare`` with a queueSize cap"""
    profile = _default_profile(path, var)
    if _slurm_children(profile):
        assert re.search(r"queueSize", _code(path)), (
            f"{path.name}: {var} defaults to {profile} without capping queueSize"
        )
    else:
        assert "farmshare_local" in profile.split(","), profile
