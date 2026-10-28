"""Static checks of the MVP-04 Slurm scripts (docs/specs/MVP-04-ecoli-evidence-run.md "The run",
"Where the run sits", lead decision 18; GIPPA_SPEC.md §3.8, R14).

Nothing is submitted: the scripts are parsed (``bash -n``) and read as text, as in
tests/unit/test_m2_21_nightly_weekly_sbatch.py.

* ``scripts/slurm/mvp04_inputs.sbatch`` (new): the gippasim input step, job ``gippa-mvp04-inputs``,
  <= 4 CPU, <= 16 GB, <= 30 min. It runs ``python -m gippasim.real_genome`` (lead decision 6.4)
  on the committed genome and recorded JSON, with RegulonDB ``TUSet.txt`` (decision 1), into
  ``$GIPPA_ROOT/mvp04/<seed>/`` (rule 2), seed ``K12_EVIDENCE_SEED`` unless overridden.
  Conventions settled by the test author: gate 9 (TE7 on the inputs) needs a second gippasim call
  with the same seed, which the same job writes to ``$GIPPA_ROOT/mvp04/<seed>_te7/``.
* ``scripts/slurm/gippa_run.sbatch`` (existing): the acceptance run's head job, 2 CPU / 8 GB,
  ``--profile farmshare,apptainer``. These checks pass today; they guard the spec's run shape.
"""

from __future__ import annotations

import importlib
import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SLURM = REPO / "scripts" / "slurm"
INPUTS = SLURM / "mvp04_inputs.sbatch"
RUN = SLURM / "gippa_run.sbatch"


def _text(path: Path) -> str:
    assert path.is_file(), f"{path.relative_to(REPO)} does not exist (MVP-04 part 2)"
    return path.read_text()


def _code(path: Path) -> str:
    """The script without comment lines (``#SBATCH`` directives kept)."""
    keep = [ln for ln in _text(path).splitlines()
            if not ln.lstrip().startswith("#") or ln.startswith("#SBATCH")]  # fmt: skip
    return "\n".join(keep)


def _sbatch(path: Path, key: str) -> str | None:
    m = re.search(rf"^#SBATCH\s+--{key}[= ](\S+)", _text(path), re.M)
    return m.group(1) if m else None


def _mem_gb(text: str) -> float:
    m = re.fullmatch(r"(\d+(?:\.\d+)?)([KMGT]?)B?", text.upper())
    assert m, text
    scale = {"K": 1 / 1024**2, "M": 1 / 1024, "": 1 / 1024, "G": 1.0, "T": 1024.0}
    return float(m.group(1)) * scale[m.group(2)]


def _minutes(text: str) -> float:
    days, _, rest = text.rpartition("-")
    parts = [int(p) for p in rest.split(":")]
    if len(parts) == 1:
        h, mnt, s = 0, parts[0], 0
    elif len(parts) == 2:
        h, mnt, s = 0, parts[0], parts[1]
    else:
        h, mnt, s = parts
    return (int(days) if days else 0) * 1440 + h * 60 + mnt + s / 60


@pytest.mark.parametrize("path", [INPUTS, RUN], ids=lambda p: p.name)
def test_script_parses(path: Path) -> None:
    bash = shutil.which("bash")
    assert bash is not None
    _text(path)
    proc = subprocess.run([bash, "-n", str(path)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr


@pytest.mark.parametrize("path", [INPUTS, RUN], ids=lambda p: p.name)
def test_sources_farmshare_env_and_strict_mode(path: Path) -> None:
    code = _code(path)
    assert re.search(r"^\s*(source|\.)\s+\S*scripts/farmshare_env\.sh", code, re.M), path.name
    assert re.search(r"^\s*set -euo pipefail", code, re.M), path.name


@pytest.mark.parametrize("path", [INPUTS, RUN], ids=lambda p: p.name)
def test_logs_go_to_local_logs(path: Path) -> None:
    out = _sbatch(path, "output")
    assert out is not None and re.fullmatch(r"\S*/local/logs/%x\.%j\.out", out), out


# ---------------------------------------------------------------- mvp04_inputs.sbatch


def test_inputs_job_name_and_caps() -> None:
    """Spec "The run": 1 job, <= 4 CPU, <= 16 GB, <= 30 min (the gippasim step budget)."""
    assert _sbatch(INPUTS, "job-name") == "gippa-mvp04-inputs"
    cpus = _sbatch(INPUTS, "cpus-per-task")
    assert cpus is not None and 1 <= int(cpus) <= 4, cpus
    mem = _sbatch(INPUTS, "mem")
    assert mem is not None and _mem_gb(mem) <= 16, mem
    wall = _sbatch(INPUTS, "time")
    assert wall is not None and _minutes(wall) <= 30, wall
    assert _sbatch(INPUTS, "array") is None, "one job, not an array"
    assert _sbatch(INPUTS, "gres") is None, "no GPU"


def test_inputs_runs_the_gippasim_input_step_on_the_committed_inputs() -> None:
    code = _code(INPUTS)
    assert re.search(r"python\S*\s+-m\s+gippasim\.real_genome\b", code), code
    assert "tests/data/genomes/k12_mg1655.fna.gz" in code
    assert "tests/data/recorded/bakta/1.12.2-full/k12/k12_mg1655.json.gz" in code
    for flag in ("--genome", "--annotation", "--seed", "--outdir", "--tu-set"):
        assert flag in code, flag
    assert "regulondb/14.5.0/TUSet.txt" in code, "RegulonDB 14.5.0 TUs (lead decision 1)"
    assert "tests/sim" in code, "gippasim is imported from tests/sim (PYTHONPATH)"


def test_inputs_writes_under_gippa_root_mvp04() -> None:
    """Rule 2: everything under ``$GIPPA_ROOT/mvp04/<seed>/``; nothing into the checkout."""
    code = _code(INPUTS)
    outs = re.findall(r"--outdir[= ]+(\S+)", code)
    assert outs, code
    assert re.search(r"\$\{?GIPPA_ROOT\}?/mvp04/", code), "outputs under $GIPPA_ROOT/mvp04/"
    for o in outs:  # a variable (set from $GIPPA_ROOT above) or a $GIPPA_ROOT path, never the repo
        o = o.strip("\"'")
        assert o.startswith("$"), o
        assert not re.match(r"\$\{?(GIPPA_REPO|PWD)\b", o), o


def test_inputs_seed_defaults_to_k12_evidence_seed() -> None:
    seed = importlib.import_module("gippasim.scenarios").K12_EVIDENCE_SEED
    code = _code(INPUTS)
    assert str(seed) in code or "K12_EVIDENCE_SEED" in code, "the default seed (Q11)"
    assert re.search(r"\$\{(?:1|SEED):-", code), "the seed can be overridden (positional or SEED)"


def test_inputs_makes_the_second_te7_call() -> None:
    """Gate 9: two gippasim calls with the same seed; the second into ``<seed>_te7``."""
    code = _code(INPUTS)
    calls = len(re.findall(r"-m\s+gippasim\.real_genome\b", code))
    loop = re.search(r"^\s*for\s+\w+\s+in\b[^\n]*_te7", code, re.M)
    assert calls >= 2 or loop, "a second gippasim call (gate 9)"
    assert "_te7" in code


def test_inputs_never_syncs_the_shared_venv() -> None:
    """The shared venv (MEMORY: no plain ``uv run`` / ``uv sync`` / pip in jobs)."""
    code = _code(INPUTS)
    assert not re.search(r"\buv\s+sync\b|\bpip\s+install\b", code), code
    for m in re.finditer(r"\buv\s+run\b[^\n]*", code):
        assert "--no-sync" in m.group(0), m.group(0)


def test_inputs_runs_no_gippa_and_no_container() -> None:
    """The input step is gippasim only (rule 7 tuning and step 8 outputs); gippa runs elsewhere."""
    code = _code(INPUTS)
    assert not re.search(r"\bgippa\s+run\b|\bnextflow\b|\bapptainer\b", code), code


# ---------------------------------------------------------------- gippa_run.sbatch (acceptance run)


def test_run_head_job_shape() -> None:
    """Spec "Runtime budget": the head job takes 2 CPU / 8 GB; jobs named gippa-*."""
    name = _sbatch(RUN, "job-name")
    assert name is not None and name.startswith("gippa-"), name
    assert _sbatch(RUN, "cpus-per-task") == "2"
    mem = _sbatch(RUN, "mem")
    assert mem is not None and _mem_gb(mem) == 8, mem


def test_run_uses_the_farmshare_apptainer_profile() -> None:
    code = _code(RUN)
    assert re.search(r"gippa run\s+--?profile[= ]farmshare,apptainer\b", code), code
    assert '"$@"' in code
