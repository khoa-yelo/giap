"""``scripts/slurm/pilot_cleanup.sbatch``: dry-run default, the DELETED.tsv ledger, the
preconditions, the plan steps and the K2 default (docs/specs/M10-14-pilot-reproducibility.md,
"The clean-up", example 10, edge case "Partial deletes and symlinks", acceptance criteria 3-4).

Lead decision (2026-10-24) on K2: until Khoa answers, the DIA-NN reports and libraries (plan
step 2) are **kept**: the script deletes only the Nextflow work directories (plan step 1) and the
mzML files (plan step 4), and lists the reports and libraries in the dry run as "kept pending K2".
So criterion 3's "CLEAN=1 deletes exactly plan steps 1, 2 and 4" is tested as "steps 1 and 4".

Everything runs in a sandbox (``helpers.m10_14_kit.cleanup_sandbox``): a throw-away checkout, a
fake ``$GIPPA_ROOT`` in ``tmp_path``, and fake ``squeue``, ``sacct``, ``nextflow`` and ``rm`` first
on a minimal ``PATH`` (``rm`` records whether DELETED.tsv already named its target, then runs
``/bin/rm``). No Slurm, no network. No statistical claim: no calibration test (spec §8.4).
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
from helpers.m10_14_kit import (
    CLEANUP,
    DELETED_COLUMNS,
    F7,
    CleanupSandbox,
    check_house_rules,
    cleanup_runbook,
    cleanup_sandbox,
    commit_runbook,
    sha256_file,
)
from helpers.slurm_sandbox import REPO, assert_honours_gippa_repo

pytestmark = pytest.mark.skipif(sys.platform == "win32", reason="bash script")


@pytest.fixture
def cs(tmp_path: Path) -> CleanupSandbox:
    assert CLEANUP.is_file(), "scripts/slurm/pilot_cleanup.sbatch does not exist (M10-14)"
    assert_honours_gippa_repo(CLEANUP.read_text(), "pilot_cleanup.sbatch")
    return cleanup_sandbox(tmp_path)


def _out(res: subprocess.CompletedProcess[str]) -> str:
    return f"exit {res.returncode}\n--- stdout\n{res.stdout}\n--- stderr\n{res.stderr}"


def _plan_text(cs: CleanupSandbox) -> str:
    p = cs.cleanup_dir / "plan.tsv"
    assert p.is_file(), f"no plan.tsv in {cs.cleanup_dir}"
    return p.read_text()


def _missing(before: dict[str, str], after: dict[str, str]) -> set[str]:
    return set(before) - set(after)


def _changed(before: dict[str, str], after: dict[str, str]) -> set[str]:
    return {p for p in before if p in after and before[p] != after[p]}


# ---------------------------------------------------------------- house rules (criterion 4)


def test_cleanup_sbatch_house_rules() -> None:
    check_house_rules(CLEANUP, name="gippa-pilot-cleanup", cpus=2, mem_gb=8, minutes=120)


def test_cleanup_sbatch_bash32() -> None:
    assert CLEANUP.is_file(), "scripts/slurm/pilot_cleanup.sbatch does not exist"
    res = subprocess.run(
        [sys.executable, str(REPO / "scripts" / "ci" / "check_bash32.py"), str(CLEANUP)],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stdout + res.stderr


def test_cleanup_sbatch_static_guards() -> None:
    """``set -f`` (no glob spans directories), ``realpath -e``, ``nextflow clean -f -k``."""
    assert CLEANUP.is_file(), "scripts/slurm/pilot_cleanup.sbatch does not exist"
    code = "\n".join(
        ln for ln in CLEANUP.read_text().splitlines() if not ln.lstrip().startswith("#")
    )
    assert re.search(r"^\s*set -f\b|set -o noglob", code, re.M), "set -f before any rm"
    assert "realpath -e" in code
    assert re.search(r"nextflow clean -f -k|nextflow clean -k -f", code)
    assert "DRY_RUN" in code and "CLEAN" in code
    assert not re.search(r"\brm\s+-[a-z]*r[a-z]*\s+[^\n]*\*", code), "rm with a glob"


# ---------------------------------------------------------------- the dry run (default)


def test_dry_run_is_the_default_and_deletes_nothing(cs: CleanupSandbox) -> None:
    before = cs.snapshot()
    res = cs.run()
    assert res.returncode == 0, _out(res)
    after = cs.snapshot()
    assert not _missing(before, after), sorted(_missing(before, after))
    assert not _changed(before, after), sorted(_changed(before, after))
    assert not any("-f" in c["argv"] for c in cs.nextflow_calls()), cs.nextflow_calls()
    assert cs.rm_log() == [], cs.rm_log()
    assert cs.repo_status() == "", f"wrote into $GIPPA_REPO: {cs.repo_status()}"


def test_dry_run_plan_lists_the_targets(cs: CleanupSandbox) -> None:
    res = cs.run(DRY_RUN="1")
    assert res.returncode == 0, _out(res)
    plan = _plan_text(cs)
    rows = [ln for ln in plan.splitlines() if ln.strip()]
    assert len(rows) >= 2 and "\t" in rows[0], "plan.tsv: a header line and one row per target"
    for work in cs.work_dirs:
        assert str(work) in plan, f"plan step 1 target {work} not in plan.tsv"
    assert str(cs.mzml_dir) in plan or all(
        str(f) in plan for f in cs.delete_files if f.suffix == ".mzML"
    ), "plan step 4 (mzML)"
    notes = cs.root / "runs" / "pilot" / "pass1" / "notes" / "work"
    assert not any(str(notes) == c.split("\t")[0] for c in rows), "work/ without .nextflow/"


def test_dry_run_lists_reports_and_libraries_as_kept_pending_k2(cs: CleanupSandbox) -> None:
    res = cs.run()
    assert res.returncode == 0, _out(res)
    lines = (_plan_text(cs) + "\n" + res.stdout + "\n" + res.stderr).splitlines()
    for name in ("A1.speclib", "A1.predicted.speclib", "A1_report.parquet"):
        hits = [ln for ln in lines if name in ln]
        assert hits, f"{name} is not listed in the dry run"
        assert all("K2" in ln for ln in hits), f"{name} not marked as kept pending K2: {hits}"


# ---------------------------------------------------------------- CLEAN=1


def test_clean_deletes_work_dirs_and_mzml_only(cs: CleanupSandbox) -> None:
    before = cs.snapshot()
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    after = cs.snapshot()
    expected = {str(p) for p in cs.delete_files}
    assert _missing(before, after) == expected, (
        f"unexpected deletions: {sorted(_missing(before, after) - expected)}; "
        f"not deleted: {sorted(expected - _missing(before, after))}"
    )
    assert not _changed(before, after), sorted(_changed(before, after))
    for p in cs.keep_files:
        assert p.is_file() and before[str(p)] == sha256_file(p), f"kept file changed: {p}"
    for p in cs.k2_files:
        assert p.is_file(), f"{p} deleted although K2 is unanswered (lead decision: keep)"
    for w in cs.work_dirs:
        assert not w.exists(), f"{w} still exists"
    for d in cs.run_dirs:
        assert (d / ".nextflow").is_dir(), "only OUT/work goes; .nextflow/ stays"
    assert cs.repo_status() == "", f"wrote into $GIPPA_REPO: {cs.repo_status()}"


def test_clean_calls_nextflow_clean_in_each_run_dir(cs: CleanupSandbox) -> None:
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    calls = cs.nextflow_calls()
    cleaned = {
        os.path.realpath(str(c["cwd"]))
        for c in calls
        if c["argv"] and c["argv"][0] == "clean" and {"-f", "-k"} <= set(c["argv"])
    }
    for d in cs.run_dirs:
        assert os.path.realpath(d) in cleaned, f"no `nextflow clean -f -k` from {d}: {calls}"
    notes = cs.root / "runs" / "pilot" / "pass1" / "notes"
    assert os.path.realpath(notes) not in cleaned


def test_ledger_rows_precede_each_delete(cs: CleanupSandbox) -> None:
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    rm = cs.rm_log()
    assert rm, "no rm call under $GIPPA_ROOT (the plan's explicit paths go to rm)"
    unledgered = [p for state, p in rm if state != "ledgered"]
    assert not unledgered, f"rm before the DELETED.tsv row was written: {unledgered}"
    rows = cs.ledger()
    for w in cs.work_dirs:
        mine = [r for r in rows if r["path"] == str(w)]
        assert len(mine) >= 2, f"{w}: a row before the delete and one completing it: {mine}"
        assert mine[0]["action"] != "deleted" and mine[-1]["action"] == "deleted", mine
        assert mine[-1]["job"] == "4242", mine[-1]
        assert re.fullmatch(r"\d+", mine[-1]["bytes"]), mine[-1]
        assert mine[-1]["sha256"] == "NA", "a directory has no single sha256"
        assert re.match(r"\d{4}-\d{2}-\d{2}T", mine[-1]["utc"]), mine[-1]


def test_ledger_hashes_deleted_regular_files(cs: CleanupSandbox) -> None:
    before = cs.snapshot()
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    for r in cs.ledger():
        assert set(r) == set(DELETED_COLUMNS), r
        if r["sha256"] != "NA":
            assert r["path"] in before, f"hashed path {r['path']} was not a regular file"
            assert r["sha256"] == before[r["path"]], r


def test_example_10_resume_after_kill(cs: CleanupSandbox) -> None:
    """Killed after a delete but before its row was completed: the re-run lists it ``absent``."""
    gone = cs.work_dirs[0]
    cs.cleanup_dir.mkdir(parents=True, exist_ok=True)
    (cs.cleanup_dir / "DELETED.tsv").write_text(
        "\t".join(DELETED_COLUMNS)
        + "\n"
        + "\t".join([str(gone), "dir", "123", "NA", "deleting", "2026-10-24T01:00:00Z", "4241"])
        + "\n"
    )
    subprocess.run(["/bin/rm", "-rf", str(gone)], check=True)
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    rows = cs.ledger()
    assert rows[0]["job"] == "4241", "the ledger is appended to, never rewritten"
    assert any(r["path"] == str(gone) and r["action"] == "absent" for r in rows), rows
    for w in cs.work_dirs[1:]:
        assert not w.exists(), f"the re-run did not continue with {w}"
    for p in cs.keep_files:
        assert p.is_file(), p


def test_symlink_target_is_skipped(cs: CleanupSandbox) -> None:
    real = cs.outside / "work_real"
    (real / "ab").mkdir(parents=True)
    (real / "ab" / "task.out").write_text("someone else's data\n")
    run = cs.root / "runs" / "pilot" / "pass1" / "out_link"
    (run / ".nextflow").mkdir(parents=True)
    (run / "work").symlink_to(real)
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    assert (run / "work").is_symlink(), "the link itself was removed"
    assert (real / "ab" / "task.out").is_file(), "followed a symlink out of $GIPPA_ROOT"
    text = (cs.cleanup_dir / "DELETED.tsv").read_text() + _plan_text(cs) + res.stdout
    assert any(str(run / "work") in ln and "skipped_symlink" in ln for ln in text.splitlines()), (
        "the symlink target is not listed as skipped_symlink"
    )


def test_empty_plan_is_exit_0_with_an_empty_ledger(cs: CleanupSandbox) -> None:
    for p in cs.delete_files:
        p.unlink()
    for w in cs.work_dirs:
        subprocess.run(["/bin/rm", "-rf", str(w)], check=True)
    subprocess.run(["/bin/rm", "-rf", str(cs.mzml_dir)], check=True)
    for d in cs.run_dirs:
        subprocess.run(["/bin/rm", "-rf", str(d / ".nextflow")], check=True)
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    assert [r for r in cs.ledger() if r["action"] == "deleted"] == []


# ---------------------------------------------------------------- preconditions: exit 2


def _refuses(cs: CleanupSandbox, **env: str) -> None:
    before = cs.snapshot()
    res = cs.run(CLEAN="1", **env)
    assert res.returncode == 2, _out(res)
    after = cs.snapshot()
    assert not _missing(before, after), f"deleted despite exit 2: {sorted(_missing(before, after))}"
    assert cs.rm_log() == [], cs.rm_log()
    assert not any("-f" in c["argv"] for c in cs.nextflow_calls())


def test_refuses_on_failing_evidence(cs: CleanupSandbox) -> None:
    ev = cs.root / "runs" / "pilot" / "freeze" / F7 / "evidence" / "A1.json"
    ev.write_text('{"arm": "A1", "failures": ["commit_mismatch:fedcba9"]}')
    _refuses(cs)


def test_refuses_on_failing_compare_evidence(cs: CleanupSandbox) -> None:
    ev = cs.root / "runs" / "pilot" / "freeze" / F7 / "evidence" / "compare.json"
    ev.write_text('{"failures": ["exit:1"]}')
    _refuses(cs)


def test_refuses_without_the_verify_json(cs: CleanupSandbox) -> None:
    (cs.root / "pilot" / "reproduce" / f"verify-{F7}.json").unlink()
    _refuses(cs)


def test_refuses_on_verify_failures(cs: CleanupSandbox) -> None:
    """Example 5: an altered DIA-NN copy fails verify, and S15 refuses to run."""
    (cs.root / "pilot" / "reproduce" / f"verify-{F7}.json").write_text(
        '{"failures": ["mismatch:diann_image"], "reconvert": {"a1": "identical"}}'
    )
    _refuses(cs)


def test_refuses_without_the_freeze_metrics(cs: CleanupSandbox) -> None:
    (cs.root / "pilot" / "metrics" / f"freeze-{F7}" / "metrics.json").unlink()
    _refuses(cs)


@pytest.mark.parametrize(
    "name",
    [
        "gippa-run-freeze-a1",
        "gippa-diann-A1",
        "gippa-entrapment",
        "gippa-pilot-compare",
        "gippa-bench",
    ],
)
def test_refuses_while_pilot_jobs_run(cs: CleanupSandbox, name: str) -> None:
    cs.set_jobs(
        [
            {"id": 4242, "name": "gippa-pilot-cleanup", "state": "RUNNING"},
            {"id": 5001, "name": name, "state": "PENDING" if "bench" in name else "RUNNING"},
        ]
    )
    _refuses(cs)


def test_own_job_and_unrelated_jobs_do_not_block(cs: CleanupSandbox) -> None:
    cs.set_jobs(
        [
            {"id": 4242, "name": "gippa-pilot-cleanup", "state": "RUNNING"},
            {"id": 5002, "name": "gippa-test-branch", "state": "RUNNING"},
            {"id": 5003, "name": "gippa-nightly", "state": "PENDING"},
        ]
    )
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    assert not cs.work_dirs[0].exists()


def test_refuses_a_target_resolving_outside_gippa_root(cs: CleanupSandbox) -> None:
    """``realpath -e`` of every target must lie inside $GIPPA_ROOT (an ancestor link leads out)."""
    te7 = cs.root / "runs" / "pilot" / "te7"
    moved = cs.outside / "te7"
    subprocess.run(["/bin/mv", str(te7), str(moved)], check=True)
    te7.symlink_to(moved)
    before_out = sorted(str(p) for p in moved.rglob("*") if p.is_file())
    _refuses(cs)
    assert sorted(str(p) for p in moved.rglob("*") if p.is_file()) == before_out


def test_refuses_a_target_that_is_a_kept_output(cs: CleanupSandbox) -> None:
    """No target is an Inputs row or an Outputs row with kept = yes."""
    work = cs.work_dirs[2]
    row = f"| pass1_work | $GIPPA_ROOT/{work.relative_to(cs.root)} | {'a' * 64} | - | yes |"
    commit_runbook(cs, cleanup_runbook(extra_outputs=[row]))
    _refuses(cs)


def test_refuses_a_target_that_is_an_input(cs: CleanupSandbox) -> None:
    mz = next(f for f in cs.delete_files if f.suffix == ".mzML")
    row = (
        f"| wt1_mzml | $GIPPA_DATA_DIR/{mz.relative_to(cs.root / 'data')} | {sha256_file(mz)} | "
        "$GIPPA_ROOT/pilot/A1/provenance.json |"
    )
    commit_runbook(cs, cleanup_runbook(extra_inputs=[row]))
    _refuses(cs)
