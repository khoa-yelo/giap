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
import re
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


def test_outputs_are_moved_only_after_the_checks(box: ConvertSandbox) -> None:
    """Steps 6-8: summary and contract checks see temporary files only, then the move."""
    box.stage("A1")
    assert box.run("A1").returncode == 0
    summaries, pytests = box.calls("summary"), box.calls("pytest")
    assert len(summaries) == 1, summaries
    assert "-I" in summaries[0]["argv"]
    assert _flag(summaries[0]["argv"], "--arm") == "A1"
    assert len(pytests) == 1, pytests
    argv = pytests[0]["argv"]
    assert "slow" in argv and "full_report" in argv, argv
    assert any(a.endswith("test_m10_03_diann_contract.py") for a in argv), argv
    for c in (*summaries, *pytests):
        assert not set(c["present"] or []) & _tables("a1"), (c["kind"], c["present"])


def test_example10_two_runs_same_bytes(box: ConvertSandbox) -> None:
    box.stage("A1")
    assert box.run("A1").returncode == 0
    first = {n: sha256(box.inputs / n) for n in ("a1.tsv.gz", "a1_abundance.tsv.gz")}
    doc1 = json.loads((box.inputs / "a1.convert.json").read_text())
    box.reset_calls()
    assert box.run("A1").returncode == 0
    second = {n: sha256(box.inputs / n) for n in ("a1.tsv.gz", "a1_abundance.tsv.gz")}
    doc2 = json.loads((box.inputs / "a1.convert.json").read_text())
    assert first == second
    assert {k: v for k, v in doc1.items() if k != "job"} == {
        k: v for k, v in doc2.items() if k != "job"
    }


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


def test_converter_failure_exits_1_with_failed_marker(box: ConvertSandbox) -> None:
    rows = [{k: v for k, v in r.items() if k != "Stripped.Sequence"} for r in job_rows()]
    box.stage("A1", table(rows))
    res = box.run("A1")
    assert res.returncode == 1, (res.returncode, res.stdout, res.stderr)
    _no_tables(box)
    marker = box.inputs / "a1.FAILED"
    assert marker.is_file(), box.inputs_files()
    assert re.search(r"\b1\b", marker.read_text()), "the marker records the exit code"
    assert box.repo_is_clean() == ""


def test_stale_failed_marker_is_removed_on_success(box: ConvertSandbox) -> None:
    box.stage("A1")
    box.inputs.mkdir(parents=True)
    (box.inputs / "a1.FAILED").write_text("exit 1\n")
    assert box.run("A1").returncode == 0
    assert set(box.inputs_files()) == _outputs("a1")


def test_nondeterministic_second_output_exits_1(box: ConvertSandbox) -> None:
    box.stage("A1")
    res = box.run("A1", FAKE_CONVERT_NONDET="1")
    assert res.returncode == 1, (res.returncode, res.stdout, res.stderr)
    _no_tables(box)


def test_summary_failure_exits_1_and_moves_nothing(box: ConvertSandbox) -> None:
    """A report with no decoy converts, but summary rule 2 fails it."""
    box.stage("A1", table(job_rows(decoys=False)))
    res = box.run("A1")
    assert res.returncode == 1, (res.returncode, res.stdout, res.stderr)
    _no_tables(box)


def test_contract_failure_exits_1_and_moves_nothing(box: ConvertSandbox) -> None:
    box.stage("A1")
    res = box.run("A1", FAKE_PYTEST="fail")
    assert res.returncode == 1, (res.returncode, res.stdout, res.stderr)
    _no_tables(box)
