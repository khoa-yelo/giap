"""The two M8-12b Slurm jobs, run against stubs (docs/specs/M8-12b-golden-median-exact.md §2.1,
§3.6, §4.2, example 6, B4, B8; GIPPA_SPEC.md §3.8, R14).

Nothing reaches Slurm, micromamba, R or the real checkout: each job runs in a
``helpers.slurm_sandbox`` checkout with a fake ``micromamba`` (it answers the version check and
the lock export, and "runs" ``make_golden_r.R`` by copying the committed M8-12 outputs) and, for
the Wallenius job, a stand-in ``make_wallenius_exact.py``.

Interface fixed here (spec §2.1, §3.6):

- ``golden_r.sbatch``: exits 1 before running R when ``GOLDEN_R_REASON`` is unset; exits 1
  before touching ``tests/data/golden/r/`` when the new ``ora_genes``, ``ora_categories`` or
  ``outputs/goseq.tsv.gz`` differ from the committed ones, unless
  ``GOLDEN_R_ALLOW_ORA_CHANGE=1``; a first run (nothing committed) passes. PROVENANCE.md gains
  ``- Reason: $GOLDEN_R_REASON`` and ``- Supersedes: Slurm job <old id>`` (or ``none``).
- ``golden_wallenius.sbatch [checkout]``: requires ``SLURM_JOB_ID``; runs
  ``$UV_PROJECT_ENVIRONMENT/bin/python scripts/golden/make_wallenius_exact.py`` twice with
  ``--workers 4`` and ``--out`` outside the checkout, ``cmp``s the two files and only then
  writes ``tests/data/golden/wallenius/{exact_tails.tsv.gz,PROVENANCE.md}``.
"""

from __future__ import annotations

import gzip
import hashlib
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from helpers.slurm_sandbox import Sandbox, assert_honours_gippa_repo, make_sandbox, write_exe

REPO = Path(__file__).resolve().parents[2]
GOLDEN_R = REPO / "scripts" / "slurm" / "golden_r.sbatch"
WALLENIUS = REPO / "scripts" / "slurm" / "golden_wallenius.sbatch"
REAL_GENERATOR = REPO / "scripts" / "golden" / "make_wallenius_exact.py"
COMMITTED = REPO / "tests" / "data" / "golden" / "r"
DEST_R = Path("tests/data/golden/r")
DEST_W = Path("tests/data/golden/wallenius")
JOB = "5550001"
REASON = "M8-12b: median-exact GSEA inputs (test)"

_MICROMAMBA = r"""#!{python}
import json, os, shutil, sys
argv = sys.argv[1:]
state = os.environ["FAKE_SLURM_DIR"]
with open(os.path.join(state, "micromamba.jsonl"), "a") as fh:
    fh.write(json.dumps(argv) + "\n")
if argv[:1] == ["create"]:
    sys.exit("fake micromamba: no env builds in tests")
if argv[:2] == ["env", "export"]:  # the committed lock, so the re-export is byte-identical
    sys.stdout.write(open(os.environ["FAKE_LOCK"]).read())
    sys.exit(0)
if argv[:1] == ["run"] and "-e" in argv:
    code = argv[argv.index("-e") + 1]
    if "packageVersion" in code:
        print("versions match the pins")
    else:
        print("4.5.3", end="")
    sys.exit(0)
if argv[:1] == ["run"] and any(a.endswith("make_golden_r.R") for a in argv):
    i = [j for j, a in enumerate(argv) if a.endswith("make_golden_r.R")][0]
    _inp, out, versions = argv[i + 1 : i + 4]
    shutil.copytree(os.environ["FAKE_R_OUTPUTS"], out)
    if os.environ.get("FAKE_GOSEQ_CHANGE"):
        with open(os.path.join(out, "goseq.tsv.gz"), "ab") as fh:
            fh.write(b"\0")
    with open(versions, "w") as fh:
        fh.write("fake VERSIONS\n")
    sys.exit(0)
sys.exit("fake micromamba: unexpected call " + " ".join(argv))
"""


def _git(sb: Sandbox, *args: str) -> None:
    subprocess.run(["git", *args], cwd=sb.repo, env=sb.env, check=True, capture_output=True)


def _commit_all(sb: Sandbox, msg: str) -> None:
    _git(sb, "add", "-A")
    _git(sb, "commit", "-q", "--allow-empty", "-m", msg)


def _sandbox(tmp: Path, *, committed: bool = True) -> Sandbox:
    sb = make_sandbox(tmp)
    (sb.repo / "tests" / "calibration").mkdir(parents=True)
    shutil.copy(REPO / "tests/calibration/tolerances.yaml", sb.repo / "tests/calibration")
    if committed:
        shutil.copytree(COMMITTED, sb.repo / DEST_R)
    _commit_all(sb, "golden data")
    mamba = tmp / "mamba"
    prefix = mamba / "envs" / "gippa-golden-r"
    prefix.mkdir(parents=True)
    yml = sb.repo / "scripts" / "golden" / "golden_r.yml"
    (prefix / ".gippa-built").write_text(hashlib.sha256(yml.read_bytes()).hexdigest() + "\n")
    shutil.copytree(COMMITTED / "outputs", tmp / "fake_r_outputs")
    write_exe(sb.bin / "micromamba", _MICROMAMBA.replace("{python}", sys.executable))
    sb.env.update(
        {
            "MAMBA_ROOT_PREFIX": str(mamba),
            "FAKE_R_OUTPUTS": str(tmp / "fake_r_outputs"),
            "SLURM_JOB_ID": JOB,
            "FAKE_LOCK": str(REPO / "scripts" / "golden" / "golden_r.explicit.txt"),
        }
    )
    return sb


def _run_golden_r(sb: Sandbox, **env: str) -> subprocess.CompletedProcess[str]:
    script = sb.repo / "scripts" / "slurm" / "golden_r.sbatch"
    assert_honours_gippa_repo(script.read_text(), "golden_r.sbatch")
    return sb.run(["bash", str(script), str(sb.repo)], env=env)


def _r_ran(sb: Sandbox) -> bool:
    return any(any(a.endswith("make_golden_r.R") for a in c) for c in sb.jsonl("micromamba.jsonl"))


def _snapshot(root: Path) -> dict[str, bytes]:
    return {
        str(p.relative_to(root)): p.read_bytes() for p in sorted(root.rglob("*")) if p.is_file()
    }


# --- golden_r.sbatch ----------------------------------------------------------------------------


def test_golden_r_requires_a_reason(tmp_path: Path) -> None:
    sb = _sandbox(tmp_path)
    before = _snapshot(sb.repo / DEST_R)
    proc = _run_golden_r(sb)
    assert proc.returncode != 0, proc.stdout
    assert "GOLDEN_R_REASON" in proc.stdout + proc.stderr
    assert not _r_ran(sb), "exits before running R"
    assert _snapshot(sb.repo / DEST_R) == before


def test_golden_r_records_reason_and_supersedes(tmp_path: Path) -> None:
    sb = _sandbox(tmp_path)
    proc = _run_golden_r(sb, GOLDEN_R_REASON=REASON)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    text = (sb.repo / DEST_R / "PROVENANCE.md").read_text()
    assert f"- Reason: {REASON}" in text.splitlines()
    prior = re.search(r"^- Slurm job: (\d+)", (COMMITTED / "PROVENANCE.md").read_text(), re.M)
    assert prior, "the committed PROVENANCE names its Slurm job"
    assert f"- Supersedes: Slurm job {prior.group(1)}" in text.splitlines(), "the committed job"
    assert re.search(rf"Slurm job: {JOB} on ", text), "the new job"
    # the M8-12 template is kept (tests/golden/test_golden_r_files.py::test_provenance_...)
    assert "GOLDEN_SEED" in text and "50000" in text
    assert re.search(r"commit\W{0,3}[0-9a-f]{40}\b", text)
    assert "(scripts modified after this commit)" not in text


def test_golden_r_first_run_supersedes_none(tmp_path: Path) -> None:
    sb = _sandbox(tmp_path, committed=False)
    proc = _run_golden_r(sb, GOLDEN_R_REASON=REASON)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    text = (sb.repo / DEST_R / "PROVENANCE.md").read_text()
    assert "- Supersedes: none" in text.splitlines()
    assert f"- Reason: {REASON}" in text.splitlines()


def _disturb_committed_ora(sb: Sandbox) -> bytes:
    path = sb.repo / DEST_R / "inputs" / "ora_genes.tsv.gz"
    path.write_bytes(path.read_bytes() + b"\0")  # the committed bytes no longer match a run
    _commit_all(sb, "disturb ORA")
    return path.read_bytes()


def test_golden_r_ora_guard_exits_before_writing(tmp_path: Path) -> None:
    """Example 6: changed ORA bytes stop the job before tests/data/golden/r/ is touched."""
    sb = _sandbox(tmp_path)
    _disturb_committed_ora(sb)
    before = _snapshot(sb.repo / DEST_R)
    proc = _run_golden_r(sb, GOLDEN_R_REASON=REASON)
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert "ora_genes" in proc.stdout + proc.stderr
    assert _snapshot(sb.repo / DEST_R) == before


def test_golden_r_goseq_guard_exits_before_writing(tmp_path: Path) -> None:
    sb = _sandbox(tmp_path)
    before = _snapshot(sb.repo / DEST_R)
    proc = _run_golden_r(sb, GOLDEN_R_REASON=REASON, FAKE_GOSEQ_CHANGE="1")
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert "goseq" in proc.stdout + proc.stderr
    assert _snapshot(sb.repo / DEST_R) == before


def test_golden_r_ora_guard_override(tmp_path: Path) -> None:
    sb = _sandbox(tmp_path)
    disturbed = _disturb_committed_ora(sb)
    proc = _run_golden_r(sb, GOLDEN_R_REASON=REASON, GOLDEN_R_ALLOW_ORA_CHANGE="1")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    new = (sb.repo / DEST_R / "inputs" / "ora_genes.tsv.gz").read_bytes()
    assert new != disturbed
    assert new == (COMMITTED / "inputs" / "ora_genes.tsv.gz").read_bytes()


# --- golden_wallenius.sbatch: the house rules ----------------------------------------------------


def _wtext() -> str:
    assert WALLENIUS.is_file(), f"missing {WALLENIUS.relative_to(REPO)}"
    return WALLENIUS.read_text()


def test_wallenius_job_header_and_house_rules() -> None:
    text = _wtext()
    proc = subprocess.run(["bash", "-n", str(WALLENIUS)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert text.splitlines()[0] == "#!/bin/bash"
    assert re.search(r"^#SBATCH --job-name=gippa-golden-wallenius$", text, re.M)
    assert re.search(r"^#SBATCH --output=\S*/local/logs/%x\.%j\.out$", text, re.M)
    cpus = re.search(r"^#SBATCH --cpus-per-task=(\d+)$", text, re.M)
    assert cpus and int(cpus.group(1)) <= 4, "ADR-0044: <= 4 CPUs"
    mem = re.search(r"^#SBATCH --mem=(\d+)G$", text, re.M)
    assert mem and int(mem.group(1)) <= 4
    hours = re.search(r"^#SBATCH --time=(\d\d):(\d\d):(\d\d)$", text, re.M)
    assert hours and (int(hours.group(1)), int(hours.group(2))) <= (1, 0), "ADR-0044: <= 1 h"
    assert "set -euo pipefail" in text
    assert "scripts/farmshare_env.sh" in text
    assert "SLURM_JOB_ID:?" in text
    assert re.search(r"\$\{?UV_PROJECT_ENVIRONMENT\}?/bin/python", text)
    assert not re.search(r"\buv\s+(run|sync|pip)\b", text)
    assert "make_wallenius_exact.py" in text
    assert re.search(r"--workers[ =]4\b", text)
    assert re.search(r"\bcmp\b", text), "the run-twice byte check"
    assert "tests/data/golden/wallenius" in text
    assert_honours_gippa_repo(text, "golden_wallenius.sbatch")


# --- golden_wallenius.sbatch: runs against a stand-in generator ----------------------------------

_FAKE_NONDETERMINISTIC = r"""#!/usr/bin/env python3
# Stand-in make_wallenius_exact.py: every run writes different bytes.
import gzip, json, os, sys, time
argv = sys.argv[1:]
state = os.environ["FAKE_SLURM_DIR"]
with open(os.path.join(state, "wallenius.jsonl"), "a") as fh:
    fh.write(json.dumps(argv) + "\n")
out = argv[argv.index("--out") + 1]
n = sum(1 for _ in open(os.path.join(state, "wallenius.jsonl")))
text = "category\ttail\tp_exact\nCAT_001\tover\t0.5\nCAT_001\tunder\t0.%d\n" % (5 + n)
open(out, "wb").write(gzip.compress(text.encode(), mtime=0))
if "--log" in argv:
    open(argv[argv.index("--log") + 1], "w").write("{}\n")
"""

_FAKE_REAL_SUBSET = r"""#!/usr/bin/env python3
# Stand-in make_wallenius_exact.py: the real generator restricted to CAT_010 (one R2 point).
import json, os, subprocess, sys
argv = sys.argv[1:]
with open(os.path.join(os.environ["FAKE_SLURM_DIR"], "wallenius.jsonl"), "a") as fh:
    fh.write(json.dumps(argv) + "\n")
real = os.environ["FAKE_REAL_GENERATOR"]
sys.exit(subprocess.call([sys.executable, real, *argv, "--categories", "CAT_010"]))
"""


def _wallenius_sandbox(tmp: Path, fake: str) -> Sandbox:
    sb = _sandbox(tmp)
    assert (sb.repo / "scripts" / "slurm" / "golden_wallenius.sbatch").is_file(), (
        "missing scripts/slurm/golden_wallenius.sbatch"
    )
    gen = sb.repo / "scripts" / "golden" / "make_wallenius_exact.py"
    gen.write_text(fake)
    _commit_all(sb, "stand-in generator")
    sb.env["FAKE_REAL_GENERATOR"] = str(REAL_GENERATOR)
    return sb


def _run_wallenius(
    sb: Sandbox, *, unset: tuple[str, ...] = (), timeout: float = 180
) -> subprocess.CompletedProcess[str]:
    script = sb.repo / "scripts" / "slurm" / "golden_wallenius.sbatch"
    assert_honours_gippa_repo(script.read_text(), "golden_wallenius.sbatch")
    return sb.run(["bash", str(script), str(sb.repo)], unset=unset, timeout=timeout)


def test_wallenius_job_requires_slurm_job_id(tmp_path: Path) -> None:
    sb = _wallenius_sandbox(tmp_path, _FAKE_NONDETERMINISTIC)
    proc = _run_wallenius(sb, unset=("SLURM_JOB_ID",))
    assert proc.returncode != 0
    assert not sb.jsonl("wallenius.jsonl"), "the generator never ran"
    assert not (sb.repo / DEST_W).exists()


def test_wallenius_job_stops_when_the_two_runs_differ(tmp_path: Path) -> None:
    sb = _wallenius_sandbox(tmp_path, _FAKE_NONDETERMINISTIC)
    proc = _run_wallenius(sb)
    assert proc.returncode != 0, proc.stdout + proc.stderr
    calls = sb.jsonl("wallenius.jsonl")
    assert len(calls) == 2, "two runs, then the byte check"
    assert not (sb.repo / DEST_W / "exact_tails.tsv.gz").exists()
    assert not (sb.repo / DEST_W / "PROVENANCE.md").exists()


@pytest.mark.slow
def test_wallenius_job_writes_the_file_and_provenance(tmp_path: Path) -> None:
    sb = _wallenius_sandbox(tmp_path, _FAKE_REAL_SUBSET)
    proc = _run_wallenius(sb, timeout=600)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    calls = sb.jsonl("wallenius.jsonl")
    assert len(calls) == 2
    outs = []
    for argv in calls:
        assert argv[argv.index("--workers") + 1] == "4"
        out = Path(argv[argv.index("--out") + 1])
        assert sb.repo not in out.parents, "runs write outside the checkout"
        outs.append(out)
    assert outs[0] != outs[1]
    dest = sb.repo / DEST_W
    text = gzip.decompress((dest / "exact_tails.tsv.gz").read_bytes()).decode()
    rows = [ln.split("\t") for ln in text.splitlines()]
    assert rows[0] == ["category", "tail", "p_exact"]
    assert rows[1] == ["CAT_010", "over", "6.0681368749300294e-46"]
    text = (dest / "PROVENANCE.md").read_text()
    assert JOB in text
    assert "golden_wallenius.sbatch" in text
    assert re.search(rf"commit\W{{0,3}}{sb.head}\b", text)
    assert "byte-identical" in text
    assert "ADR-0044" in text
    assert re.search(r"nmant\W{0,3}63", text)
    for name in ("ora_genes.tsv.gz", "ora_categories.tsv.gz"):
        digest = hashlib.sha256((COMMITTED / "inputs" / name).read_bytes()).hexdigest()
        assert digest in text, name
