"""``scripts/slurm/diann_convert.sbatch <arm>`` in a sandbox (M10-03 criteria 5 and 7).

Spec: docs/specs/M10-03-diann-convert.md, "The conversion job", examples 8-10, edge case
"Partial outputs". The sandbox (``helpers.diann_convert_m10.convert_sandbox``) runs the **real**
converter of this checkout on a hand-made 10-row report; only ``-m pytest`` (step 7) is faked,
and ``GIAP_PYTHON`` is a recording wrapper. The determinism failure (Q8) uses the wrapper's
``FAKE_CONVERT_NONDET`` (it appends a byte to the second conversion's output), so the script
needs no test hook.

Fixed here: ``$OUT_DIR`` overrides ``$GIAP_ROOT/pilot``; ``SLURM_JOB_ID`` names the temporary
files; the default ``GIAP_PYTHON`` is ``$UV_PROJECT_ENVIRONMENT/bin/python``.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner
from helpers.diann_convert_m10 import (
    ARMS,
    JOB_ID,
    SCRIPT,
    ConvertSandbox,
    convert_sandbox,
    job_rows,
    sha256,
    table,
)

from giap.cli import cli


@pytest.fixture
def box(tmp_path: Path) -> ConvertSandbox:
    return convert_sandbox(tmp_path)


def _outputs(stem: str) -> set[str]:
    return {
        f"{stem}.tsv.gz",
        f"{stem}_abundance.tsv.gz",
        f"{stem}.convert.log",
        f"{stem}.convert.json",
    }


def _tables(stem: str) -> set[str]:
    """The final names M10-04 reads (``convert.json`` may be written at step 6)."""
    return {f"{stem}.tsv.gz", f"{stem}_abundance.tsv.gz"}


def _flag(argv: list[str], flag: str) -> str:
    assert flag in argv, (flag, argv)
    return argv[argv.index(flag) + 1]


def _nothing_written(box: ConvertSandbox, res) -> None:  # type: ignore[no-untyped-def]
    assert res.returncode == 2, (res.returncode, res.stdout, res.stderr)
    assert box.inputs_files() == [], f"wrote {box.inputs_files()} on a precondition failure"
    assert not box.inputs.exists(), "$PILOT/inputs was created on a precondition failure (M13)"
    assert box.calls("convert") == [], "giap convert ran despite a failed precondition"
    assert box.repo_is_clean() == ""


# ------------------------------------------------------------------ static house rules


def _text() -> str:
    assert SCRIPT.is_file(), f"{SCRIPT} does not exist (M10-03)"
    return SCRIPT.read_text()


def _directive(key: str) -> str | None:
    m = re.search(rf"^#SBATCH\s+--{key}(?:[= ](\S+))?", _text(), re.M)
    return (m.group(1) or "") if m else None


def test_sbatch_house_rules() -> None:
    """Criterion 7 (R14, §3.8)."""
    job = _directive("job-name")
    assert job and job.startswith("giap-"), f"job name {job!r}"
    assert _directive("partition") == "normal"
    assert _directive("qos") == "normal"
    assert _directive("cpus-per-task") == "2"
    assert _directive("mem") == "8G"
    assert _directive("time") == "00:30:00"
    out = _directive("output")
    assert out and "/local/logs/" in out, f"--output {out}"
    code = "\n".join(ln for ln in _text().splitlines() if not ln.lstrip().startswith("#"))
    assert "scripts/farmshare_env.sh" in code
    assert re.search(r"^set -euo pipefail$", code, re.M)
    assert not re.search(r"\buv\s+(run|sync)\b|\bpip3?\s+install\b|\bpip\b", code), (
        "no uv run, uv sync or pip"
    )


# ------------------------------------------------------------------ arguments and preconditions


@pytest.mark.parametrize("args", [(), ("A3",), ("A1", "A2"), ("a1",)], ids=repr)
def test_bad_arguments_exit_2(box: ConvertSandbox, args: tuple[str, ...]) -> None:
    box.stage("A1")
    _nothing_written(box, box.run(*args))


def test_example8_check_not_ok_exits_2(box: ConvertSandbox) -> None:
    box.stage("A1", ok=False)
    res = box.run("A1")
    _nothing_written(box, res)
    assert "check.json" in res.stdout + res.stderr


def test_missing_check_json_exits_2(box: ConvertSandbox) -> None:
    box.stage("A1")
    (box.arm_dir("A1") / "check.json").unlink()
    res = box.run("A1")
    _nothing_written(box, res)
    assert "check.json" in res.stdout + res.stderr


def test_example9_report_rewritten_after_the_search_exits_2(box: ConvertSandbox) -> None:
    rep = box.stage("A2", stale_sum=True)
    res = box.run("A2")
    _nothing_written(box, res)
    msg = res.stdout + res.stderr
    assert sha256(rep) in msg and "0" * 64 in msg, "the message names both hashes"


def test_missing_outputs_sha256_exits_2(box: ConvertSandbox) -> None:
    box.stage("A1", sums=False)
    _nothing_written(box, box.run("A1"))


@pytest.mark.parametrize("state", ["missing", "empty"])
def test_missing_or_empty_report_exits_2(box: ConvertSandbox, state: str) -> None:
    rep = box.stage("A1")
    if state == "missing":
        rep.unlink()
    else:
        rep.write_bytes(b"")
    res = box.run("A1")
    _nothing_written(box, res)
    assert rep.name in res.stdout + res.stderr


def test_failed_marker_exits_2(box: ConvertSandbox) -> None:
    box.stage("A1", failed=True)
    res = box.run("A1")
    _nothing_written(box, res)
    assert "FAILED" in res.stdout + res.stderr


# ------------------------------------------------------------------ success


@pytest.mark.parametrize("arm", sorted(ARMS))
def test_success_writes_the_checked_outputs(box: ConvertSandbox, arm: str) -> None:
    stem, search_db = ARMS[arm]
    rep = box.stage(arm)
    res = box.run(arm)
    assert res.returncode == 0, (res.stdout, res.stderr)
    assert set(box.inputs_files()) == _outputs(stem), box.inputs_files()
    assert box.inputs.stat().st_mode & 0o777 == 0o755
    assert box.repo_is_clean() == "", "the job wrote into the checkout"

    for name in (f"{stem}.tsv.gz", f"{stem}_abundance.tsv.gz"):
        assert (box.inputs / name).read_bytes()[:2] == b"\x1f\x8b", f"{name} is not gzip"
    # the tables are this checkout's conversion of the report, byte for byte
    ref = box.sb.tmp / "ref"
    ref.mkdir()
    cli_res = CliRunner().invoke(
        cli,
        [
            "convert",
            "--from",
            "diann",
            str(rep),
            "--out",
            str(ref / "a.tsv.gz"),
            "--abundance-out",
            str(ref / "ab.tsv.gz"),
            "--search-db",
            search_db,
        ],
    )
    assert cli_res.exit_code == 0, cli_res.output
    assert (box.inputs / f"{stem}.tsv.gz").read_bytes() == (ref / "a.tsv.gz").read_bytes()
    assert (box.inputs / f"{stem}_abundance.tsv.gz").read_bytes() == (
        ref / "ab.tsv.gz"
    ).read_bytes()

    doc = json.loads((box.inputs / f"{stem}.convert.json").read_text())
    assert doc["failures"] == []
    assert doc["arm"] == arm and doc["search_db"] == search_db
    assert doc["report"]["sha256"] == sha256(rep)
    assert doc["counts"]["rows"] == 10 and doc["counts"]["decoy_rows"] == 2


def test_the_convert_calls(box: ConvertSandbox) -> None:
    """Steps 3-5: branch code, temporary names, a second conversion outside ``$OUT``."""
    box.stage("A1")
    res = box.run("A1")
    assert res.returncode == 0, (res.stdout, res.stderr)
    calls = box.calls("convert")
    assert len(calls) == 2, [c["argv"] for c in calls]
    for c in calls:
        assert str(box.repo / "src") in c["pythonpath"].split(":"), c["pythonpath"]
        assert _flag(c["argv"], "--from") == "diann"
        assert _flag(c["argv"], "--search-db") == "other"
        assert str(box.report("A1")) in c["argv"]
        assert c["rc"] == 0
    first, second = calls
    for flag in ("--out", "--abundance-out"):
        p = Path(_flag(first["argv"], flag))
        assert p.parent == box.inputs, p
        # ``.<name>.tmp.$SLURM_JOB_ID``; a trailing ``.gz`` may be kept so the CLI still gzips
        assert p.name.startswith(".") and f"tmp.{JOB_ID}" in p.name, p.name
        q = Path(_flag(second["argv"], flag))
        assert box.inputs not in q.parents and q.parent != box.inputs, (
            f"the second conversion goes to $TMPDIR, not $OUT: {q}"
        )
    for c in calls:  # no final table exists while converting
        assert not set(c["present"] or []) & _tables("a1"), c["present"]


@pytest.mark.parametrize("arm", sorted(ARMS))
def test_outputs_are_moved_only_after_the_checks(box: ConvertSandbox, arm: str) -> None:
    """Steps 6-8: summary and contract checks see temporary files only, then the move.

    Review B1: step 7 checks this job's report (``GIAP_M10_03_REPORT``) and this arm only
    (``GIAP_M10_03_ARM``, ``-k "full_report and <arm>"``).
    """
    stem = ARMS[arm][0]
    box.stage(arm)
    res = box.run(arm)
    assert res.returncode == 0, (res.stdout, res.stderr)
    summaries, pytests = box.calls("summary"), box.calls("pytest")
    assert len(summaries) == 1, summaries
    assert "-I" in summaries[0]["argv"]
    assert _flag(summaries[0]["argv"], "--arm") == arm
    assert len(pytests) == 1, pytests
    argv = pytests[0]["argv"]
    assert "slow" in argv, argv
    assert _flag(argv, "-k") == f"full_report and {arm}", argv
    assert any(a.endswith("test_m10_03_diann_contract.py") for a in argv), argv
    assert pytests[0]["env"] == {
        "GIAP_M10_03_REPORT": str(box.report(arm)),
        "GIAP_M10_03_ARM": arm,
    }, pytests[0]["env"]
    for c in (*summaries, *pytests):
        assert not set(c["present"] or []) & _tables(stem), (c["kind"], c["present"])


def test_contract_step_checks_the_out_dir_report(box: ConvertSandbox) -> None:
    """Review B1: with ``OUT_DIR`` the contract step gets ``$OUT_DIR/<arm>/<arm>_report.parquet``,
    not ``$GIAP_ROOT/pilot/…``."""
    box.stage("A1")
    alt = box.sb.tmp / "alt"
    (alt / "A1").mkdir(parents=True)
    for p in box.arm_dir("A1").iterdir():
        p.rename(alt / "A1" / p.name)
    res = box.run("A1", OUT_DIR=str(alt))
    assert res.returncode == 0, (res.stdout, res.stderr)
    (call,) = box.calls("pytest")
    assert call["env"]["GIAP_M10_03_REPORT"] == str(alt / "A1" / "A1_report.parquet"), call
    assert call["env"]["GIAP_M10_03_ARM"] == "A1", call


def test_all_skipped_contract_run_fails_the_job(box: ConvertSandbox) -> None:
    """Review B1: pytest exits 0 when every test is skipped; the job counts that as a failure."""
    box.stage("A1")
    res = box.run("A1", FAKE_PYTEST="skip")
    assert len(box.calls("pytest")) == 1
    assert res.returncode == 1, (res.returncode, res.stdout, res.stderr)
    _failed_cleanly(box)


def test_example10_two_runs_same_bytes(box: ConvertSandbox) -> None:
    """Two runs under different ``SLURM_JOB_ID``s (review S1: the temporary names then differ,
    so a temporary name leaking into ``convert.json`` breaks the equality)."""
    box.stage("A1")
    res = box.run("A1", SLURM_JOB_ID="4242")
    assert res.returncode == 0, (res.stdout, res.stderr)
    first = {n: sha256(box.inputs / n) for n in ("a1.tsv.gz", "a1_abundance.tsv.gz")}
    doc1 = json.loads((box.inputs / "a1.convert.json").read_text())
    box.reset_calls()
    res = box.run("A1", SLURM_JOB_ID="5151")
    assert res.returncode == 0, (res.stdout, res.stderr)
    second = {n: sha256(box.inputs / n) for n in ("a1.tsv.gz", "a1_abundance.tsv.gz")}
    doc2 = json.loads((box.inputs / "a1.convert.json").read_text())
    assert first == second
    assert (doc1["job"]["id"], doc2["job"]["id"]) == ("4242", "5151")
    assert {k: v for k, v in doc1.items() if k != "job"} == {
        k: v for k, v in doc2.items() if k != "job"
    }
    assert set(box.inputs_files()) == _outputs("a1"), box.inputs_files()


@pytest.mark.parametrize("arm", sorted(ARMS))
def test_convert_json_names_the_final_paths(box: ConvertSandbox, arm: str) -> None:
    """Review S1 (mutants M1, M2, M15): ``command`` and ``outputs.*.path`` hold the final names,
    never the temporary ones; the provenance is that of the sandbox checkout."""
    stem, search_db = ARMS[arm]
    rep = box.stage(arm)
    res = box.run(arm)
    assert res.returncode == 0, (res.stdout, res.stderr)
    doc = json.loads((box.inputs / f"{stem}.convert.json").read_text())
    assert doc["command"] == [
        "giap",
        "convert",
        "--from",
        "diann",
        str(rep),
        "--out",
        str(box.inputs / f"{stem}.tsv.gz"),
        "--abundance-out",
        str(box.inputs / f"{stem}_abundance.tsv.gz"),
        "--search-db",
        search_db,
    ], doc["command"]
    assert doc["report"]["path"] == str(rep)
    assert {k: v["path"] for k, v in doc["outputs"].items()} == {
        "peptides": str(box.inputs / f"{stem}.tsv.gz"),
        "abundance": str(box.inputs / f"{stem}_abundance.tsv.gz"),
        "log": str(box.inputs / f"{stem}.convert.log"),
    }
    for k, name in (("peptides", f"{stem}.tsv.gz"), ("abundance", f"{stem}_abundance.tsv.gz")):
        assert doc["outputs"][k]["sha256"] == sha256(box.inputs / name), k
    # review B3: the code that ran is the sandbox checkout, committed; bytecode is not "dirty"
    assert doc["giap"]["commit"] == box.sb.head
    assert doc["giap"]["dirty"] is False, doc["giap"]


def test_out_dir_overrides_the_pilot_root(box: ConvertSandbox) -> None:
    box.stage("A1")
    alt = box.sb.tmp / "alt"
    (alt / "A1").mkdir(parents=True)
    for p in box.arm_dir("A1").iterdir():
        p.rename(alt / "A1" / p.name)
    res = box.run("A1", OUT_DIR=str(alt))
    assert res.returncode == 0, (res.stdout, res.stderr)
    assert {p.name for p in (alt / "inputs").iterdir()} == _outputs("a1")
    assert not box.inputs.exists()


# ------------------------------------------------------------- failures after the preconditions


def _no_tables(box: ConvertSandbox) -> None:
    left = [n for n in box.inputs_files() if ".tsv" in n]
    assert not left, f"tables left behind: {left}"


def _failed_cleanly(box: ConvertSandbox, stem: str = "a1") -> str:
    """Exit-1 aftermath (review S3, S5): ``<a>.FAILED`` and nothing else under ``$OUT``, no
    hidden temporary file, nothing in the checkout. Returns the marker's text."""
    _no_tables(box)
    hidden = [n for n in box.inputs_files() if n.startswith(".") or ".tmp." in n]
    assert not hidden, f"temporary files left behind: {hidden}"
    assert box.inputs_files() == [f"{stem}.FAILED"], box.inputs_files()
    assert box.repo_is_clean() == ""
    return (box.inputs / f"{stem}.FAILED").read_text()


def test_converter_failure_exits_1_with_failed_marker(box: ConvertSandbox) -> None:
    rows = [{k: v for k, v in r.items() if k != "Stripped.Sequence"} for r in job_rows()]
    box.stage("A1", table(rows))
    res = box.run("A1")
    assert res.returncode == 1, (res.returncode, res.stdout, res.stderr)
    _no_tables(box)
    text = _failed_cleanly(box)
    assert re.search(r"\b1\b", text), "the marker records the exit code"


def test_stale_failed_marker_is_removed_on_success(box: ConvertSandbox) -> None:
    box.stage("A1")
    box.inputs.mkdir(parents=True)
    (box.inputs / "a1.FAILED").write_text("exit 1\n")
    assert box.run("A1").returncode == 0
    assert set(box.inputs_files()) == _outputs("a1")


@pytest.mark.parametrize("flag", ["--out", "--abundance-out"])
def test_nondeterministic_second_output_exits_1(box: ConvertSandbox, flag: str) -> None:
    """Either table differing between the two conversions fails the job (review S5, M9)."""
    box.stage("A1")
    res = box.run("A1", FAKE_CONVERT_NONDET=flag)
    assert res.returncode == 1, (res.returncode, res.stdout, res.stderr)
    _failed_cleanly(box)


def test_summary_failure_exits_1_and_moves_nothing(box: ConvertSandbox) -> None:
    """A report with no decoy converts, but summary rule 2 fails it."""
    box.stage("A1", table(job_rows(decoys=False)))
    res = box.run("A1")
    assert res.returncode == 1, (res.returncode, res.stdout, res.stderr)
    _failed_cleanly(box)


def test_contract_failure_exits_1_and_moves_nothing(box: ConvertSandbox) -> None:
    box.stage("A1")
    res = box.run("A1", FAKE_PYTEST="fail")
    assert res.returncode == 1, (res.returncode, res.stdout, res.stderr)
    _failed_cleanly(box)


def test_failed_marker_holds_the_converter_exit_code(box: ConvertSandbox) -> None:
    """Review S5 (M14): the marker carries the converter's own exit code, not just ``exit 1``."""
    box.stage("A1")
    res = box.run("A1", FAKE_CONVERT_RC="7")
    assert res.returncode == 1, (res.returncode, res.stdout, res.stderr)
    text = _failed_cleanly(box)
    assert re.search(r"\b7\b", text), f"the converter's exit code 7 is not in a1.FAILED: {text!r}"


# ------------------------------------------------------------- re-runs (review S3)


def _stale_rerun_setup(box: ConvertSandbox) -> None:
    box.stage("A1")
    res = box.run("A1", SLURM_JOB_ID="1000")
    assert res.returncode == 0, (res.stdout, res.stderr)
    assert set(box.inputs_files()) == _outputs("a1")
    box.reset_calls()


@pytest.mark.parametrize("failure", ["converter", "determinism", "summary", "contract"])
def test_failed_rerun_leaves_no_stale_outputs(box: ConvertSandbox, failure: str) -> None:
    """Review S3: once the preconditions pass, the previous outputs are moved aside, and deleted
    when the re-run fails; M10-04 never reads a table of an older report next to ``FAILED``."""
    _stale_rerun_setup(box)
    env = {"SLURM_JOB_ID": "2000"}
    if failure == "converter":
        rows = [{k: v for k, v in r.items() if k != "Stripped.Sequence"} for r in job_rows()]
        box.stage("A1", table(rows))
    elif failure == "summary":
        box.stage("A1", table(job_rows(decoys=False)))
    elif failure == "determinism":
        env["FAKE_CONVERT_NONDET"] = "--out"
    else:
        env["FAKE_PYTEST"] = "fail"
    res = box.run("A1", **env)
    assert res.returncode == 1, (res.returncode, res.stdout, res.stderr)
    _failed_cleanly(box)


def test_successful_rerun_leaves_only_the_new_outputs(box: ConvertSandbox) -> None:
    """Review S3: nothing moved aside is left behind after a successful re-run."""
    _stale_rerun_setup(box)
    rep = box.stage("A1", table(job_rows()[:9]))
    res = box.run("A1", SLURM_JOB_ID="2000")
    assert res.returncode == 0, (res.stdout, res.stderr)
    assert set(box.inputs_files()) == _outputs("a1"), box.inputs_files()
    doc = json.loads((box.inputs / "a1.convert.json").read_text())
    assert doc["report"]["sha256"] == sha256(rep) and doc["counts"]["rows"] == 9


def test_precondition_failure_on_rerun_keeps_the_previous_outputs(box: ConvertSandbox) -> None:
    """Review S3: outputs are moved aside only *after* the preconditions pass."""
    _stale_rerun_setup(box)
    before = {n: sha256(box.inputs / n) for n in _outputs("a1")}
    (box.arm_dir("A1") / "check.json").write_text(json.dumps({"ok": False}))
    res = box.run("A1", SLURM_JOB_ID="2000")
    assert res.returncode == 2, (res.returncode, res.stdout, res.stderr)
    assert {n: sha256(box.inputs / n) for n in box.inputs_files()} == before
    assert box.calls("convert") == []


# ------------------------------------------------------------- set -e exits (review S4)

_ROOT = pytest.mark.skipif(
    hasattr(os, "geteuid") and os.geteuid() == 0, reason="root reads unreadable files"
)


@_ROOT
@pytest.mark.parametrize("what", ["report", "outputs.sha256"])
def test_precondition_read_error_exits_2(box: ConvertSandbox, what: str) -> None:
    """Review S4: an I/O error while checking a precondition is exit 2, not ``set -e``'s 1."""
    rep = box.stage("A1")
    target = rep if what == "report" else box.arm_dir("A1") / "outputs.sha256"
    target.chmod(0)
    try:
        res = box.run("A1")
    finally:
        target.chmod(0o644)
    _nothing_written(box, res)


def test_later_command_failure_writes_failed(box: ConvertSandbox) -> None:
    """Review S4: a failing command after the preconditions (here the publishing ``mv``) is
    exit 1 with ``<a>.FAILED`` (an ERR trap), not a bare ``set -e`` exit."""
    box.stage("A1")
    res = box.run("A1", FAKE_MV_FAIL="1")
    assert res.returncode == 1, (res.returncode, res.stdout, res.stderr)
    assert (box.inputs / "a1.FAILED").is_file(), box.inputs_files()
    hidden = [n for n in box.inputs_files() if ".tmp." in n]
    assert not hidden, f"temporary files left behind: {hidden}"


# ------------------------------------------------------------- $OUT inside the checkout (S2)


def _move_arm(box: ConvertSandbox, pilot: Path, arm: str = "A1") -> None:
    (pilot / arm).mkdir(parents=True)
    for p in box.arm_dir(arm).iterdir():
        p.rename(pilot / arm / p.name)


def test_out_dir_inside_the_checkout_and_not_ignored_exits_2(box: ConvertSandbox) -> None:
    """Review S2: an ``OUT_DIR`` such as ``$GIAP_REPO/validation/pilot`` would write derived
    tables where ``git add -A`` stages them; it is a precondition failure. (``validation/pilot``
    itself already holds tracked M10-01 inputs, so an untracked, unignored sibling is used.)"""
    box.stage("A1")
    pilot = box.repo / "derived" / "pilot"
    _move_arm(box, pilot)
    res = box.run("A1", OUT_DIR=str(pilot))
    assert res.returncode == 2, (res.returncode, res.stdout, res.stderr)
    assert not (pilot / "inputs").exists(), sorted(p.name for p in (pilot / "inputs").iterdir())
    assert box.calls("convert") == []


def test_out_dir_under_an_ignored_directory_of_the_checkout_runs(box: ConvertSandbox) -> None:
    """Review S2 / B2: the FarmShare default ``GIAP_ROOT=$GIAP_REPO/local`` (git-ignored)."""
    (box.repo / ".gitignore").write_text("/local/\n")
    for a in (["add", ".gitignore"], ["commit", "-q", "-m", "ignore local/"]):
        subprocess.run(["git", "-C", str(box.repo), *a], check=True, capture_output=True,
                       env=box.sb.env)  # fmt: skip
    box.stage("A1")
    pilot = box.repo / "local" / "pilot"
    _move_arm(box, pilot)
    res = box.run("A1", OUT_DIR=str(pilot))
    assert res.returncode == 0, (res.returncode, res.stdout, res.stderr)
    assert {p.name for p in (pilot / "inputs").iterdir()} == _outputs("a1")
    doc = json.loads((pilot / "inputs" / "a1.convert.json").read_text())
    assert doc["giap"]["dirty"] is False, doc["giap"]
