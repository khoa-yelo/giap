"""``scripts/slurm/cascadia_convert.sbatch`` in a sandbox (M10-10a criteria 3 and 12).

Spec: docs/specs/M10-10-a3-run.md, "``cascadia_convert.sbatch`` (step 3)", examples 4-6, edge case
"Empty input". The sandbox (``helpers.a3_run_m10_10.convert_sandbox``) runs the **real** M10-09
converter and the real ``cascadia_convert.py`` of this checkout on three tiny ``.ssl`` written by
the test; only ``a3_feasibility.py evaluate`` (M10-08b) is faked, by exit code (``FAKE_L6``).
Nothing is submitted to Slurm and nothing touches the real ``$GIAP_ROOT``.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest
from click.testing import CliRunner
from helpers.a3_run_m10_10 import (
    CONVERT_JOB,
    JOB_SSL,
    RUNS,
    SEQ_LIKE,
    ConvertSandbox,
    convert_sandbox,
    feasibility_doc,
    sha256,
    strings,
)

from giap.cli import cli

OUTPUTS = {"a3.tsv.gz", "a3.convert.json", "a3.convert.log"}


@pytest.fixture
def box(tmp_path: Path) -> ConvertSandbox:
    return convert_sandbox(tmp_path)


def _msg(res) -> str:  # type: ignore[no-untyped-def]
    return (
        f"rc={res.returncode}\n--- stdout\n{res.stdout[-3000:]}\n--- stderr\n{res.stderr[-3000:]}"
    )


# ------------------------------------------------------------------ criterion 12: house rules


def _text() -> str:
    assert CONVERT_JOB.is_file(), f"{CONVERT_JOB.name} does not exist yet (M10-10a deliverable 2)"
    return CONVERT_JOB.read_text()


def _directive(key: str) -> str | None:
    m = re.search(rf"^#SBATCH\s+--{key}(?:[= ](\S+))?", _text(), re.M)
    return (m.group(1) or "") if m else None


def test_sbatch_house_rules() -> None:
    """Criterion 12 (R14, §3.8)."""
    job = _directive("job-name")
    assert job and job.startswith("giap-"), f"job name {job!r}"
    assert _directive("partition") == "normal"
    assert _directive("qos") == "normal"
    assert _directive("cpus-per-task") == "2"
    assert _directive("mem") == "8000M"
    assert _directive("time") == "00:30:00"
    out = _directive("output")
    assert out and "/local/logs/" in out, f"--output {out}"
    code = "\n".join(ln for ln in _text().splitlines() if not ln.lstrip().startswith("#"))
    assert "scripts/farmshare_env.sh" in code
    assert re.search(r"^set -euo pipefail$", code, re.M)
    assert re.search(r"^set -f$", code, re.M), "set -f (spec: header)"
    assert not re.search(r"\buv\s+(run|sync)\b|\bpip3?\b", code), "no uv run, uv sync or pip"
    assert "A3_OUT_DIR" in code, "the A3 root honours A3_OUT_DIR"


# ------------------------------------------------------------------ success (example 4)


def test_example4_success(box: ConvertSandbox) -> None:
    ssl = box.stage()
    res = box.run()
    assert res.returncode == 0, _msg(res)
    assert set(box.inputs_files()) == OUTPUTS, box.inputs_files()
    assert box.repo_is_clean() == "", "the job wrote into the checkout"

    # ssl.sha256: sha256sum format, one line per run, sorted by run
    lines = (box.a3 / "ssl.sha256").read_text().splitlines()
    assert len(lines) == 3
    parsed = [re.fullmatch(r"([0-9a-f]{64}) [ *](\S+)", ln) for ln in lines]
    assert all(parsed), lines
    names = [Path(m.group(2)).name for m in parsed if m]
    assert names == ["WT_1.ssl", "WT_2.ssl", "WT_3.ssl"]
    assert [m.group(1) for m in parsed if m] == [sha256(p) for p in ssl]

    for r in RUNS:  # the per-run L6 check left its aggregates next to each .ssl
        feas = box.a3 / r / "feasibility.json"
        assert feas.is_file(), f"no {feas}"
    assert sorted(c["run"] for c in box.calls("evaluate")) == list(RUNS)

    # the table is this checkout's conversion of the three .ssl, in run order, byte for byte
    ref = box.sb.tmp / "ref.tsv.gz"
    argv = ["convert", "--from", "cascadia", *map(str, ssl), "--out", str(ref), "--source",
            "denovo"]  # fmt: skip
    cli_res = CliRunner().invoke(cli, argv)
    assert cli_res.exit_code == 0, cli_res.output
    table = box.inputs / "a3.tsv.gz"
    assert table.read_bytes()[:2] == b"\x1f\x8b", "a3.tsv.gz is gzip"
    assert table.read_bytes() == ref.read_bytes()
    assert len(box.calls("convert")) == 2, "two conversions (step 5, determinism)"

    doc = json.loads((box.inputs / "a3.convert.json").read_text())
    assert doc["failures"] == []
    assert doc["counts"]["rows"] == sum(len(v) for v in JOB_SSL.values())
    assert doc["counts"]["rows_per_run"] == {r: len(JOB_SSL[r]) for r in RUNS}
    assert doc["outputs"]["peptides"]["sha256"] == sha256(table)
    assert doc["outputs"]["peptides"]["path"] == str(table), "the final path, not a temporary"
    bad = [s for s in strings(doc) if re.fullmatch(SEQ_LIKE, s)]
    assert bad == [], bad


def test_a3_out_dir_overrides_the_a3_root(box: ConvertSandbox) -> None:
    other = box.sb.root / "elsewhere" / "A3"
    box.stage(a3=other)
    res = box.run(A3_OUT_DIR=str(other))
    assert res.returncode == 0, _msg(res)
    assert (other / "ssl.sha256").is_file()
    assert not box.a3.exists(), "the default $GIAP_ROOT/pilot/A3 was used despite A3_OUT_DIR"
    assert set(box.inputs_files()) == OUTPUTS


def test_stale_failed_marker_is_removed(box: ConvertSandbox) -> None:
    box.stage()
    box.inputs.mkdir(parents=True)
    (box.inputs / "a3.FAILED").write_text("exit 1\n")
    res = box.run()
    assert res.returncode == 0, _msg(res)
    assert set(box.inputs_files()) == OUTPUTS


# ------------------------------------------------------------------ the gate and preconditions


def _wrote_nothing(box: ConvertSandbox) -> None:
    assert box.inputs_files() == [], box.inputs_files()
    assert not (box.a3 / "ssl.sha256").exists()
    assert not list(box.a3.glob("*/feasibility.json"))
    assert box.calls("convert") == [] and box.calls("evaluate") == []
    assert box.repo_is_clean() == ""


@pytest.mark.parametrize("reason", ["low_genome_match", "route_unavailable"])
def test_gate_fail_exits_3_and_writes_nothing(box: ConvertSandbox, reason: str) -> None:
    box.stage()
    box.set_gate(feasibility_doc(verdict="fail", reason=reason))
    res = box.run()
    assert res.returncode == 3, _msg(res)
    assert "M10-11a" in res.stdout + res.stderr
    _wrote_nothing(box)


def test_gate_missing_feasibility_exits_2(box: ConvertSandbox) -> None:
    box.stage()
    box.set_gate(None)
    res = box.run()
    assert res.returncode == 2, _msg(res)
    _wrote_nothing(box)


@pytest.mark.parametrize("state", ["missing", "empty"])
def test_missing_or_empty_ssl_exits_2(box: ConvertSandbox, state: str) -> None:
    ssl = box.stage()
    if state == "missing":
        ssl[1].unlink()
    else:
        ssl[1].write_bytes(b"")
    res = box.run()
    assert res.returncode == 2, _msg(res)
    assert "WT_2" in res.stdout + res.stderr, "the message names the run"
    assert box.inputs_files() == [] and box.calls("convert") == []


def test_arguments_are_rejected(box: ConvertSandbox) -> None:
    box.stage()
    res = box.run("WT_1")
    assert res.returncode == 2, _msg(res)
    assert box.inputs_files() == []


# ------------------------------------------------------------------ per-run L6 (examples 5, 6)


def test_example5_one_run_fails_l6_and_all_three_are_converted(box: ConvertSandbox) -> None:
    box.stage()
    res = box.run(FAKE_L6="WT_3=3")
    assert res.returncode == 0, _msg(res)
    assert set(box.inputs_files()) == OUTPUTS
    feas = json.loads((box.a3 / "WT_3" / "feasibility.json").read_text())
    assert feas["verdict"] == "fail"
    doc = json.loads((box.inputs / "a3.convert.json").read_text())
    assert sorted(doc["counts"]["rows_per_run"]) == list(RUNS), "never a data-chosen subset"


@pytest.mark.parametrize("plan", ["WT_2=3,WT_3=3", "WT_1=3,WT_2=3,WT_3=3"])
def test_example6_l6_fails_on_two_or_more_runs_exits_3(box: ConvertSandbox, plan: str) -> None:
    box.stage()
    res = box.run(FAKE_L6=plan)
    assert res.returncode == 3, _msg(res)
    n = plan.count("=3")
    assert f"L6 fails on {n} of 3 runs: lead decision before A3 runs" in res.stdout + res.stderr
    assert "a3.tsv.gz" not in box.inputs_files()
    assert box.calls("convert") == [], "no conversion after the L6 stop"
    assert box.repo_is_clean() == ""


@pytest.mark.parametrize("rc", ["1", "2"])
def test_l6_evaluation_error_exits_1(box: ConvertSandbox, rc: str) -> None:
    box.stage()
    res = box.run(FAKE_L6=f"WT_2={rc}")
    assert res.returncode == 1, _msg(res)
    assert "a3.tsv.gz" not in box.inputs_files()


# ------------------------------------------------------------------ conversion failures


def test_converter_error_exits_1_with_failed_marker(box: ConvertSandbox) -> None:
    box.stage(score_type={"WT_2": "PERCOLATOR QVALUE"})
    res = box.run()
    assert res.returncode == 1, _msg(res)
    files = box.inputs_files()
    assert "a3.tsv.gz" not in files and "a3.convert.json" not in files
    assert "a3.FAILED" in files, files
    assert not [f for f in files if f.startswith(".") or ".tmp" in f], f"temporaries left: {files}"
    assert box.repo_is_clean() == ""


def test_nondeterministic_conversion_exits_1(box: ConvertSandbox) -> None:
    box.stage()
    res = box.run(FAKE_CONVERT_NONDET="1")
    assert res.returncode == 1, _msg(res)
    files = box.inputs_files()
    assert "a3.tsv.gz" not in files and "a3.convert.json" not in files, files


def test_summary_rule_failure_exits_1(box: ConvertSandbox) -> None:
    """A run label that is not in acquisition.json fails summary rule 2 (D1 stems ≠ the runs)."""
    ssl = box.stage()
    text = ssl[2].read_text().replace("/WT_3.mzML", "/WT_3b.mzML")
    ssl[2].write_text(text)
    res = box.run()
    assert res.returncode == 1, _msg(res)
    assert "a3.tsv.gz" not in box.inputs_files()
    assert "a3.convert.json" not in box.inputs_files()
