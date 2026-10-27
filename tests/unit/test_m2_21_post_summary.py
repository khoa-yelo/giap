"""``scripts/ci/post_summary.py`` behaviour with a fake ``gh`` (M2-21 lead decisions 4 and 17).

``gh`` is the fake of ``tests/helpers/slurm_sandbox.py`` first on ``PATH``: ``issue list`` answers
from a state file, ``issue create`` adds issue #77 to it, ``issue comment`` records its body; every
call is logged. Nothing reaches GitHub.

Pinned (decision 4, 17; reviewer's notes on PR #107):

* the "Nightly (FarmShare)" issue is created once (only when no open issue has exactly that
  title) and every array job gets exactly one comment on it, also when the script runs again
  (two tasks finishing together, or the afterany summary job after the last task);
* fewer rows than ``--expect``: nothing is posted and the exit is 0;
* ``--dry-run`` prints the comment, calls no ``gh`` and leaves no ``.posted`` mark;
* the comment is green only when every task exited 0 and there is at least one row; it names the
  job, the commit and every task; the last row of a task wins; ``|`` in a summary cannot break the
  table;
* a failed ``gh issue comment`` exits non-zero and does not mark the job posted (a later run
  retries).
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
from helpers.slurm_sandbox import FAKE_GH, FAKE_SACCT, write_exe

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "ci" / "post_summary.py"
TITLE = "Nightly (FarmShare)"


class Env:
    def __init__(self, tmp: Path) -> None:
        self.tmp = tmp
        self.state = tmp / "state"
        self.state.mkdir()
        bindir = tmp / "bin"
        bindir.mkdir()
        write_exe(bindir / "gh", FAKE_GH.replace("{python}", sys.executable))
        write_exe(bindir / "sacct", FAKE_SACCT.replace("{python}", sys.executable))
        self.tsv = tmp / "logs" / "nightly.tsv"
        self.tsv.parent.mkdir()
        self.env = {
            "PATH": f"{bindir}:/usr/bin:/bin",
            "HOME": str(tmp),
            "FAKE_SLURM_DIR": str(self.state),
            "GH_TOKEN": "fake-token-never-used",
            "GH_HOST": "github.invalid",
        }

    def rows(self, *rows: tuple[str, str, str, str]) -> None:
        with self.tsv.open("a") as fh:
            for job_task, commit, rc, summary in rows:
                fh.write(f"2026-10-20\t{job_task}\t{commit}\t{rc}\t{summary}\n")

    def run(self, *args: str, **env: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(SCRIPT), "--tsv", str(self.tsv), *args],
            env={**self.env, **env}, cwd=str(self.tmp), capture_output=True, text=True,
            timeout=60, check=False,
        )  # fmt: skip

    def sacct(self, *rows: tuple[str, str]) -> None:
        (self.state / "sacct_rows.json").write_text(
            json.dumps([{"JobID": j, "State": st, "JobName": "gippa-nightly"} for j, st in rows])
        )

    def calls(self) -> list[dict[str, Any]]:
        p = self.state / "gh.jsonl"
        return [json.loads(ln) for ln in p.read_text().splitlines()] if p.is_file() else []

    def comments(self) -> list[dict[str, Any]]:
        return [c for c in self.calls() if c["argv"][:2] == ["issue", "comment"]]

    def creates(self) -> list[dict[str, Any]]:
        return [c for c in self.calls() if c["argv"][:2] == ["issue", "create"]]

    def issues(self, *issues: tuple[int, str]) -> None:
        (self.state / "gh_state.json").write_text(json.dumps(
            {"issues": [{"number": n, "title": t} for n, t in issues], "next": 77}))  # fmt: skip


@pytest.fixture
def e(tmp_path: Path) -> Env:
    return Env(tmp_path)


def _four(e: Env, job: str, rcs: str = "0000", commit: str = "abc1234") -> None:
    e.rows(*((f"{job}_{t}", commit, rcs[t], f"task {t}: 3 passed") for t in range(4)))


def test_issue_created_once_and_one_comment_per_job(e: Env) -> None:
    _four(e, "100")
    r1 = e.run("--job", "100", "--expect", "4")
    _four(e, "200")
    r2 = e.run("--job", "200", "--expect", "4")
    assert r1.returncode == 0 and r2.returncode == 0, (r1.stderr, r2.stderr)
    assert len(e.creates()) == 1, e.calls()
    create = e.creates()[0]["argv"]
    assert create[create.index("--title") + 1] == TITLE
    comments = e.comments()
    assert [c["argv"][2] for c in comments] == ["77", "77"], comments
    assert "100" in comments[0]["body"] and "200" in comments[1]["body"]
    assert "200" not in comments[0]["body"]


def test_existing_issue_is_reused_only_on_an_exact_title(e: Env) -> None:
    e.issues((12, "Nightly (FarmShare) archive"), (31, TITLE))
    _four(e, "100")
    assert e.run("--job", "100", "--expect", "4").returncode == 0
    assert not e.creates(), e.calls()
    assert [c["argv"][2] for c in e.comments()] == ["31"]


def test_near_miss_title_creates_the_issue(e: Env) -> None:
    e.issues((12, "Nightly (FarmShare) archive"))
    _four(e, "100")
    assert e.run("--job", "100", "--expect", "4").returncode == 0
    assert len(e.creates()) == 1
    assert [c["argv"][2] for c in e.comments()] == ["77"]


def test_same_job_is_posted_once(e: Env) -> None:
    _four(e, "100")
    for _ in range(3):
        assert e.run("--job", "100", "--expect", "4").returncode == 0
    assert len(e.comments()) == 1, e.comments()


def test_incomplete_array_posts_nothing(e: Env) -> None:
    e.rows(("100_0", "abc1234", "0", "a"), ("100_1", "abc1234", "0", "b"),
           ("999_2", "abc1234", "0", "other job"))  # fmt: skip
    res = e.run("--job", "100", "--expect", "4")
    assert res.returncode == 0, res.stderr
    assert not e.calls(), e.calls()


def test_dry_run_is_offline_and_leaves_no_mark(e: Env) -> None:
    _four(e, "100", rcs="0010")
    res = e.run("--job", "100", "--expect", "4", "--dry-run")
    assert res.returncode == 0, res.stderr
    assert not e.calls(), "--dry-run called gh"
    assert "100" in res.stdout and "abc1234" in res.stdout
    for t in range(4):
        assert re.search(rf"^\|\s*{t}\s*\|", res.stdout, re.M), res.stdout
    assert not any(p.name.endswith(".posted") for p in e.tsv.parent.iterdir())
    assert e.run("--job", "100", "--expect", "4").returncode == 0
    assert len(e.comments()) == 1, "a dry run must not count as posted"


@pytest.mark.parametrize(("rcs", "green"), [("0000", True), ("0100", False), ("0002", False)])
def test_green_only_when_every_task_passed(e: Env, rcs: str, green: bool) -> None:
    _four(e, "100", rcs=rcs)
    res = e.run("--job", "100", "--expect", "4", "--dry-run")
    head = res.stdout.strip().splitlines()[0].lower()
    assert ("green" in head) is green, res.stdout
    if not green:
        assert "red" in head


def test_an_empty_night_is_never_green(e: Env) -> None:
    e.rows(("999_0", "abc1234", "0", "another job"))
    res = e.run("--job", "100", "--expect", "0", "--dry-run")
    assert "green" not in res.stdout.lower(), res.stdout


def test_last_row_of_a_task_wins(e: Env) -> None:
    _four(e, "100", rcs="0300")
    e.rows(("100_1", "abc1234", "0", "task 1: rerun passed"))
    res = e.run("--job", "100", "--expect", "4", "--dry-run")
    assert "green" in res.stdout.splitlines()[0].lower(), res.stdout
    assert "rerun passed" in res.stdout
    assert len(re.findall(r"^\|\s*1\s*\|", res.stdout, re.M)) == 1, res.stdout


def test_pipes_in_a_summary_do_not_add_columns(e: Env) -> None:
    e.rows(("100_0", "abc1234", "0", "a | b | c"))
    res = e.run("--job", "100", "--expect", "1", "--dry-run")
    row = next(ln for ln in res.stdout.splitlines() if re.match(r"^\|\s*0\s*\|", ln))
    header = next(ln for ln in res.stdout.splitlines() if ln.startswith("|") and "task" in ln)
    assert row.count("|") == header.count("|"), (header, row)


def test_failed_comment_is_an_error_and_is_retried(e: Env) -> None:
    _four(e, "100")
    res = e.run("--job", "100", "--expect", "4", FAKE_GH_FAIL="comment")
    assert res.returncode != 0
    assert e.run("--job", "100", "--expect", "4").returncode == 0
    assert len(e.comments()) == 2, "the failed post must not mark the job as posted"


def test_missing_tsv_is_an_error(e: Env) -> None:
    res = e.run("--job", "100", "--expect", "4")
    assert res.returncode != 0
    assert "nightly.tsv" in res.stderr
    assert not e.calls()


def test_token_never_printed(e: Env) -> None:
    _four(e, "100")
    res = e.run("--job", "100", "--expect", "4")
    assert "fake-token-never-used" not in res.stdout + res.stderr
    assert all("fake-token-never-used" not in json.dumps(c) for c in e.calls())


# --- second review of PR #107: decisions 25, 26 -------------------------------------------------


@pytest.mark.parametrize("done", [0, 1, 2, 3])
def test_fewer_than_expect_posts_nothing(e: Env, done: int) -> None:
    """Decision 26 (``expect - 1``): 3 of 4 tasks is not yet a night."""
    e.rows(*((f"100_{t}", "abc1234", "0", f"task {t}") for t in range(done)))
    res = e.run("--job", "100", "--expect", "4")
    assert res.returncode == 0, res.stderr
    assert not e.comments(), f"posted with {done}/4 tasks: {e.comments()}"


def test_note_rows_are_not_tasks(e: Env) -> None:
    """Decision 26: a ``<job>_note`` row neither completes the night nor shows as a task."""
    e.rows(*((f"100_{t}", "abc1234", "0", f"task {t}") for t in (0, 2, 3)))
    e.rows(("100_note", "abc1234", "0", "benchmarks resubmitted as job 7001"))
    assert e.run("--job", "100", "--expect", "4").returncode == 0
    assert not e.comments(), "a note row counted as task 1"
    e.rows(("100_1", "abc1234", "0", "task 1"))
    res = e.run("--job", "100", "--expect", "4", "--dry-run")
    assert res.returncode == 0, res.stderr
    assert "green" in res.stdout.splitlines()[0].lower(), res.stdout
    assert not re.search(r"^\|\s*note\s*\|", res.stdout, re.M), res.stdout
    assert "resubmitted as job 7001" in res.stdout, "the note is not shown"
    table = [ln for ln in res.stdout.splitlines() if re.match(r"^\|\s*\d+\s*\|", ln)]
    assert len(table) == 4, table


def test_final_shows_the_sacct_state_of_a_missing_task(e: Env) -> None:
    """Decision 26: the killed task's row carries its own sacct state (other jobs ignored)."""
    e.rows(*((f"100_{t}", "abc1234", "0", f"task {t}") for t in (0, 1, 3)))
    e.sacct(("100_2", "TIMEOUT"), ("1000_2", "OUT_OF_MEMORY"), ("9999_2", "COMPLETED"))
    res = e.run("--job", "100", "--expect", "4", "--final", FAKE_SACCT_NOFILTER="1")
    assert res.returncode == 0, res.stderr
    (c,) = e.comments()
    body = c["body"]
    assert "green" not in body.splitlines()[0].lower(), body
    row = re.search(r"^\|\s*2\s*\|.*$", body, re.M)
    assert row and "TIMEOUT" in row.group(0), body
    assert "COMPLETED" not in row.group(0) and "OUT_OF_MEMORY" not in row.group(0), row.group(0)


def test_final_without_a_tsv_posts_a_red_night(e: Env) -> None:
    """Decision 25: every task was killed before writing; nightly.tsv does not exist."""
    e.sacct(*((f"100_{t}", "TIMEOUT") for t in range(4)))
    assert not e.tsv.exists()
    res = e.run("--job", "100", "--expect", "4", "--final")
    assert res.returncode == 0, res.stderr
    comments = e.comments()
    assert len(comments) == 1, f"{comments}\n{res.stderr}"
    body = comments[0]["body"]
    assert "red" in body.splitlines()[0].lower() and "100" in body, body
    for t in range(4):
        assert re.search(rf"^\|\s*{t}\s*\|.*TIMEOUT", body, re.M), body
