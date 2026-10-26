"""Static checks of ``scripts/slurm/nightly.sbatch`` and ``weekly.sbatch`` (M2-21; GIAP_SPEC.md §12,
§3.6, R14; docs/specs/M2-21-integration-nightly-runtimes.md and its lead decisions 3-5).

Nothing is submitted: the scripts are parsed (``bash -n``) and read as text.

Conventions settled by the test author:

* Array dispatch is a ``case "$SLURM_ARRAY_TASK_ID" in`` (any quoting / ``${…:-…}`` default) with
  arms ``0)`` … ``3)``; an arm's text includes the bodies of shell functions it calls (one level).
* Decision 3: ``0`` = ``tests/integration`` + ``tests/contract`` with ``--require-capabilities``
  and the ``farmshare_nightly`` list of ``ci/gates.yaml`` (named literally or read from that file);
  ``1`` = benchmarks pinned to wheat; ``2`` = the ``calibration_nightly`` tests; ``3`` = a real
  ``giap run`` on sim100k with Bakta light and with the recorded full-DB JSON.
* ``MODE`` defaults to ``test``; ``record`` runs the ``record_tools.sbatch`` logic.
* Never pushes or commits (decision 5: ``main`` only; §12 worktrees under ``$GIAP_ROOT/ci/``).
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
SLURM = REPO / "scripts" / "slurm"
NIGHTLY = SLURM / "nightly.sbatch"
WEEKLY = SLURM / "weekly.sbatch"
GATES = yaml.safe_load((REPO / "ci" / "gates.yaml").read_text())["required_capabilities"]
SCRIPTS = {"nightly": NIGHTLY, "weekly": WEEKLY}


def _text(path: Path) -> str:
    assert path.is_file(), f"{path.relative_to(REPO)} does not exist (M2-21)"
    return path.read_text()


def _code(path: Path) -> str:
    """The script without comment lines (``#SBATCH`` directives kept)."""
    keep = [ln for ln in _text(path).splitlines()
            if not ln.lstrip().startswith("#") or ln.startswith("#SBATCH")]  # fmt: skip
    return "\n".join(keep)


def _sbatch(path: Path, key: str) -> str | None:
    m = re.search(rf"^#SBATCH\s+--{key}[= ](\S+)", _text(path), re.M)
    return m.group(1) if m else None


def _functions(code: str) -> dict[str, str]:
    out: dict[str, str] = {}
    pat = re.compile(r"^(?:function\s+)?([A-Za-z_][\w-]*)\s*(?:\(\s*\))?\s*\{\s*$", re.M)
    for m in pat.finditer(code):
        end = re.search(r"^\}\s*$", code[m.end() :], re.M)
        out[m.group(1)] = code[m.end() : m.end() + end.start()] if end else ""
    return out


def _array_arms(path: Path) -> dict[str, str]:
    code = _code(path)
    m = re.search(r"case\s+\"?\$\{?SLURM_ARRAY_TASK_ID[^\s]*\s+in\b(.*?)\besac\b", code, re.S)
    assert m, f'{path.name}: no `case "$SLURM_ARRAY_TASK_ID" in … esac` dispatch'
    arms: dict[str, str] = {}
    for am in re.finditer(r"^\s*([0-9|\s\"']+)\)(.*?)(?=;;)", m.group(1), re.S | re.M):
        for label in re.split(r"[|\s\"']+", am.group(1)):
            if label:
                arms[label] = am.group(2)
    funcs = _functions(code)
    for k, body in arms.items():
        arms[k] = body + "".join(f for name, f in funcs.items() if re.search(rf"\b{name}\b", body))
    return arms


@pytest.mark.parametrize("name", sorted(SCRIPTS))
def test_script_parses(name: str) -> None:
    bash = shutil.which("bash")
    assert bash is not None
    res = subprocess.run(
        [bash, "-n", str(SCRIPTS[name])], capture_output=True, text=True, timeout=30
    )
    assert res.returncode == 0, res.stderr


@pytest.mark.parametrize("name", sorted(SCRIPTS))
def test_job_name_is_giap_prefixed(name: str) -> None:
    assert _sbatch(SCRIPTS[name], "job-name") == f"giap-{name}"


@pytest.mark.parametrize("name", sorted(SCRIPTS))
def test_mode_defaults_to_test_and_accepts_record(name: str) -> None:
    code = _code(SCRIPTS[name])
    assert re.search(r"\$\{MODE:-test\}", code), "MODE must default to test"
    assert re.search(r"\brecord\b", code) and re.search(r"\btest\b", code)


@pytest.mark.parametrize("name", sorted(SCRIPTS))
def test_never_pushes_or_commits(name: str) -> None:
    code = _code(SCRIPTS[name])
    assert not re.search(r"\bgit\b[^\n;|&]*\b(push|commit)\b", code)
    assert not re.search(r"\bgh\s+pr\s+(merge|create)\b", code)
    assert "waitmerge" not in code and "replant" not in code


@pytest.mark.parametrize("name", sorted(SCRIPTS))
def test_requires_its_gate_capabilities(name: str) -> None:
    """Example 2: the gate list is required, so a missing DB fails instead of skipping."""
    code = _code(SCRIPTS[name])
    assert "--require-capabilities" in code
    gate = f"farmshare_{name}"
    if gate not in code:
        missing = [c for c in GATES[gate] if c not in code]
        assert not missing, f"{name}: capabilities of {gate} not required: {missing}"


@pytest.mark.parametrize("name", sorted(SCRIPTS))
def test_tests_main_only_in_a_ci_worktree(name: str) -> None:
    code = _code(SCRIPTS[name])
    assert "origin/main" in code
    assert re.search(r"\$\{?GIAP_ROOT\}?/ci/", code), "worktrees live under $GIAP_ROOT/ci/ (§12)"


def test_nightly_array_matches_decision_3() -> None:
    array = _sbatch(NIGHTLY, "array")
    assert array is not None, "nightly.sbatch is a job array (#SBATCH --array)"
    ids = set()
    for part in array.split("%")[0].split(","):
        lo, _, hi = part.partition("-")
        ids |= set(range(int(lo), int(hi or lo) + 1))
    assert ids == {0, 1, 2, 3}, array
    arms = _array_arms(NIGHTLY)
    assert {"0", "1", "2", "3"} <= set(arms), sorted(arms)
    a0, a1, a2, a3 = (arms[k] for k in "0123")
    assert "tests/integration" in a0 and "tests/contract" in a0 and "--require-capabilities" in a0
    assert "benchmark" in a1 and "wheat" in a1
    assert "calibration_nightly" in a2
    assert re.search(r"\bgiap\s+run\b", a3) and "sim100k" in a3
    assert "light" in a3, "arm 3 runs Bakta with the light DB"
    assert re.search(r"--annotation|bakta_json|bakta_full", a3), (
        "and with the recorded full-DB JSON"
    )


def test_nightly_r14_cpu_cap() -> None:
    cpus = int(_sbatch(NIGHTLY, "cpus-per-task") or 1)
    array = _sbatch(NIGHTLY, "array") or "0"
    throttle = int(array.split("%")[1]) if "%" in array else 4
    assert cpus * throttle <= 128, (cpus, array)


def test_nightly_record_mode_runs_record_tools_logic() -> None:
    assert "record_tools" in _code(NIGHTLY)


def test_nightly_runs_nf_test_and_reports() -> None:
    code = _code(NIGHTLY)
    # the suites of workflow/tests, all of them or found there minus the docker-only ones (M11-14p)
    assert re.search(r"nf-test\s+test\s+(workflow/tests|\$\(find workflow/tests)", code)
    assert "nightly.tsv" in code, "decision 4: one line per sub-task in local/logs/nightly.tsv"
    assert "post_summary.py" in code
    assert (REPO / "scripts" / "ci" / "post_summary.py").is_file()


def test_scrontab_points_at_the_launcher() -> None:
    """Decision 18 supersedes decision 3's "the scrontab line keeps pointing at nightly.sbatch":
    the cron line runs ``nightly_launch.sh``, which submits nightly.sbatch as the array
    (details in ``test_m2_21_nightly_launch.py``)."""
    lines = [ln for ln in (SLURM / "nightly.scrontab").read_text().splitlines()
             if ln.strip() and not ln.lstrip().startswith("#")]  # fmt: skip
    assert len(lines) == 1 and lines[0].split()[-1].endswith("scripts/slurm/nightly_launch.sh")
    assert "nightly.sbatch" in (SLURM / "nightly_launch.sh").read_text()


def test_weekly_covers_k12_nissle_and_antismash() -> None:
    code = _code(WEEKLY).lower()
    assert "k12" in code and "nissle" in code
    assert "antismash" in code
    assert re.search(r"\bgiap\s+run\b", code)
