"""``scripts/slurm/pilot_compare.sbatch OUT_A1 OUT_A2 CMP`` (docs/specs/M10-04-pilot-runs.md,
"The compare job", acceptance criteria 5 and 6, examples 7, 8 and 9).

The job runs in a sandbox checkout (``helpers.pilot_runs_m10_04.compare_job_sandbox``: fake Slurm
tools, a logging ``gippa`` that must never be called, ``src/`` and ``validation/`` of this branch)
on two tiny arms written by ``helpers.pilot_runs_m10_04.write_arm``, with the real ``gippa compare``.

Pinned by the test-author where the spec leaves it open: the evidence goes to
``$(dirname CMP)/evidence/compare.json`` (the spec's ``<R>/evidence/compare.json`` with
``CMP = <R>/cmp_a1_a2``); when the first ``gippa compare`` exits non-zero, ``run_evidence.py`` is not
run (its exit is the job's).
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest
from helpers.pilot_runs_m10_04 import (
    COMPARE_JOB,
    TABLES,
    CompareJobSandbox,
    compare_job_sandbox,
    two_arms,
)
from helpers.slurm_sandbox import assert_honours_gippa_repo, write_exe


def _text() -> str:
    assert COMPARE_JOB.is_file(), "scripts/slurm/pilot_compare.sbatch does not exist (M10-04)"
    return COMPARE_JOB.read_text(encoding="utf-8")


@pytest.fixture
def job(tmp_path: Path) -> CompareJobSandbox:
    assert_honours_gippa_repo(_text(), "pilot_compare.sbatch")
    return compare_job_sandbox(tmp_path)


def _evidence(cmp: Path) -> Path:
    return cmp.parent / "evidence" / "compare.json"


def _out(res) -> str:  # type: ignore[no-untyped-def]
    return res.stdout + res.stderr


# ------------------------------------------------------------------ criterion 5


def test_two_compares_then_evidence_exit_0(job: CompareJobSandbox) -> None:
    """Example 7: three tables, V16, compare.html, a byte-identical second call, failures []."""
    a1, a2 = two_arms(job.runs)
    cmp = job.runs / "cmp_a1_a2"
    res = job.run(str(a1), str(a2), str(cmp))
    assert res.returncode == 0, _out(res)[-3000:]
    for t in TABLES:
        assert (cmp / "tables" / f"{t}.tsv").is_file(), t
    assert (cmp / "report" / "compare.html").stat().st_size > 0
    assert (cmp / "plots" / "V16_compare.svg").stat().st_size > 0
    doc = json.loads(_evidence(cmp).read_text())
    assert doc["failures"] == []
    assert job.repo_is_clean(), "the job wrote into the checkout"
    assert not any(ln.startswith("gippa ") for ln in job.sb.calls_log()), "gippa runs via GIPPA_PYTHON"


def test_gippa_python_runs_compare_twice_and_the_evidence(job: CompareJobSandbox) -> None:
    a1, a2 = two_arms(job.runs)
    cmp = job.runs / "cmp_a1_a2"
    res = job.run(str(a1), str(a2), str(cmp), env={"GIPPA_PYTHON": str(job.recording_python())})
    assert res.returncode == 0, _out(res)[-3000:]
    calls = job.python_calls()
    compares = [c for c in calls if "compare" in c and "-c" in c and "--outdir" in c]
    assert len(compares) == 2, calls
    for c in compares:
        assert "from gippa.cli import cli" in c[c.index("-c") + 1]
        assert c[c.index("--labels") + 1] == "A1,A2"
        assert c[c.index("compare") + 1 : c.index("compare") + 3] == [str(a1), str(a2)]
    outdirs = [c[c.index("--outdir") + 1] for c in compares]
    assert outdirs[0] == str(cmp) and outdirs[1] != str(cmp)
    assert not outdirs[1].startswith(str(job.sb.repo)), "the second compare goes to $TMPDIR"
    ev = [c for c in calls if any(a.endswith("run_evidence.py") for a in c)]
    assert len(ev) == 1, calls
    e = ev[0]
    assert "compare" in e and e[e.index("compare") + 1] == str(cmp)
    assert e[e.index("--second") + 1] == outdirs[1]
    assert e[e.index("--labels") + 1] == "A1,A2"
    assert Path(e[e.index("--out") + 1]) == _evidence(cmp)
    assert "--genes" in e, "review B4: G comes from OUT_A1/tables/genes.tsv"
    assert e[e.index("--genes") + 1] == str(a1 / "tables" / "genes.tsv")


def test_wrong_argument_count_exits_2(job: CompareJobSandbox) -> None:
    a1, a2 = two_arms(job.runs)
    for args in ([], [str(a1), str(a2)], [str(a1), str(a2), str(job.runs / "c"), "x"]):
        res = job.run(*args)
        assert res.returncode == 2, (args, _out(res))
    assert not (job.runs / "c").exists() and not (job.runs / "evidence").exists()


def test_missing_run_directory_exits_2(job: CompareJobSandbox) -> None:
    a1, _ = two_arms(job.runs)
    cmp = job.runs / "cmp_a1_a2"
    res = job.run(str(a1), str(job.runs / "out_missing"), str(cmp))
    assert res.returncode == 2, _out(res)
    assert not cmp.exists() and not (job.runs / "evidence").exists()


def test_existing_cmp_exits_2_and_writes_nothing(job: CompareJobSandbox) -> None:
    """Example 9: no merging into stale output."""
    a1, a2 = two_arms(job.runs)
    cmp = job.runs / "cmp_a1_a2"
    cmp.mkdir()
    (cmp / "old.txt").write_text("stale\n")
    res = job.run(str(a1), str(a2), str(cmp))
    assert res.returncode == 2, _out(res)
    assert "exists" in _out(res).lower()
    assert sorted(p.name for p in cmp.iterdir()) == ["old.txt"]
    assert not (job.runs / "evidence").exists()


def test_different_genome_digests_exit_1(job: CompareJobSandbox) -> None:
    """Example 8: ``gippa compare`` exit 1 (CompareInputError) is the job's exit code."""
    a1, a2 = two_arms(job.runs, sha_b="9" * 64)
    cmp = job.runs / "cmp_a1_a2"
    res = job.run(str(a1), str(a2), str(cmp))
    assert res.returncode == 1, _out(res)[-3000:]
    assert "different genomes" in _out(res)
    assert not _evidence(cmp).exists()
    assert job.repo_is_clean()


# ------------------------------------------------------------------ review S4: job exit codes

# A GIPPA_PYTHON that runs the real interpreter, except: ``FAULT=evidence_fails`` deletes
# ``CMP/plots/V16_compare.svg`` before ``run_evidence.py`` (rule C5 then fails, exit 1);
# ``FAULT=second_fails`` exits 3 for the second ``gippa compare`` (its ``--outdir`` is not CMP).
_FAULTY_PYTHON = r"""#!/bin/bash
args=("$@")
if [ "${FAULT:-}" = evidence_fails ]; then
  for a in "$@"; do
    case "$a" in *run_evidence.py) rm -f "$FAULT_CMP/plots/V16_compare.svg" ;; esac
  done
fi
if [ "${FAULT:-}" = second_fails ]; then
  for i in "${!args[@]}"; do
    if [ "${args[$i]}" = --outdir ] && [ "${args[$((i + 1))]}" != "$FAULT_CMP" ]; then
      echo "fake: second compare fails" >&2; exit 3
    fi
  done
fi
exec {python} "$@"
"""


def _faulty(job: CompareJobSandbox, fault: str, cmp: Path) -> dict[str, str]:
    exe = write_exe(
        job.sb.bin / "faulty_python", _FAULTY_PYTHON.replace("{python}", sys.executable)
    )
    return {"GIPPA_PYTHON": str(exe), "FAULT": fault, "FAULT_CMP": str(cmp)}


def test_evidence_failure_is_the_job_exit_code(job: CompareJobSandbox) -> None:
    """run_evidence.py exits 1 (a rule of set C fails): the job exits 1, the evidence is kept."""
    a1, a2 = two_arms(job.runs)
    cmp = job.runs / "cmp_a1_a2"
    res = job.run(str(a1), str(a2), str(cmp), env=_faulty(job, "evidence_fails", cmp))
    assert res.returncode == 1, _out(res)[-3000:]
    assert json.loads(_evidence(cmp).read_text())["failures"] == ["no_compare_report"]


def test_second_compare_failure_fails_the_job(job: CompareJobSandbox) -> None:
    """The second ``gippa compare`` fails: the job fails and no evidence is written."""
    a1, a2 = two_arms(job.runs)
    cmp = job.runs / "cmp_a1_a2"
    res = job.run(str(a1), str(a2), str(cmp), env=_faulty(job, "second_fails", cmp))
    assert res.returncode != 0, _out(res)[-3000:]
    assert (cmp / "tables" / "compare_detection.tsv").is_file(), "the first compare ran"
    assert not _evidence(cmp).exists()


# ------------------------------------------------------------------ criterion 6: house rules


def _directives(text: str) -> dict[str, str]:
    out = {}
    for m in re.finditer(r"^#SBATCH\s+--([\w-]+)(?:[= ](\S+))?", text, re.M):
        out[m.group(1)] = m.group(2) or ""
    return out


def test_sbatch_header() -> None:
    d = _directives(_text())
    assert d.get("job-name", "").startswith("gippa-"), d
    assert d.get("job-name") == "gippa-pilot-compare"
    assert d.get("partition") == "normal" and d.get("qos") == "normal"
    assert d.get("cpus-per-task") == "2"
    assert d.get("mem") == "8G"
    assert d.get("time") == "01:00:00"
    assert "/local/logs/" in d.get("output", ""), d


def test_sbatch_shell_rules() -> None:
    text = _text()
    code = [ln for ln in text.splitlines() if not ln.lstrip().startswith("#")]
    body = "\n".join(code)
    assert re.search(r"^set -euo pipefail\b", body, re.M)
    assert re.search(r"^set -f\b", body, re.M) or re.search(r"^set -[a-z]*f", body, re.M)
    assert re.search(r"source\s+\S*scripts/farmshare_env\.sh", body)
    assert not re.search(r"\buv\s+(run|sync)\b", body), "no uv run / uv sync"
    assert not re.search(r"(^|[\s;|&])(pip3?|python3? -m pip)\s", body), "no pip"
    assert "GIPPA_PYTHON" in body and "UV_PROJECT_ENVIRONMENT" in body
    assert re.search(r"PYTHONPATH=\S*\$\{?GIPPA_REPO\}?/src", body)
