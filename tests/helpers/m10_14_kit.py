"""Fixtures for the M10-14 tests (docs/specs/M10-14-pilot-reproducibility.md).

Nothing here reaches Slurm, the network, the real ``$GIAP_ROOT`` or any DIA-NN output:

* ``synthetic_repo`` builds a small checkout in ``tmp_path``: three ``scripts/slurm/*.sbatch``
  with known ``#SBATCH`` headers, a copy of the committed Nissle genome (so example 4's hash
  holds), a stand-in ``makedb_manifest.json`` that records that hash, the DIA-NN pin file, and,
  on request, ``validation/pilot/metrics/a3_feasibility.json``;
* ``runbook`` renders a ``reproduce.md`` in the format of the spec section "The runbook
  format"; the default text passes every static check against ``synthetic_repo``;
* ``cleanup_sandbox`` builds on ``helpers.slurm_sandbox.make_sandbox`` (a throw-away checkout of
  ``scripts/`` and ``ci/`` plus ``validation/``, fake Slurm tools first on a minimal ``PATH``) and
  lays out a fake ``$GIAP_ROOT`` with every target class of the clean-up plan, the freeze
  evidence, and fakes for ``squeue``, ``nextflow`` and ``rm`` (bash 3.2 or Python).

Conventions settled by the test author (the spec leaves them open):

* ``reproduce_check.py`` is imported by path (``load_reproduce_check``) and run as
  ``python -I validation/pilot/reproduce_check.py static RUNBOOK [--repo DIR]`` and
  ``… verify RUNBOOK --out J [--repo DIR]``, with ``$GIAP_ROOT`` from the environment; ``--repo``
  defaults to the checkout holding the script. Exit 0 when there is no failure, 1 otherwise.
* ``parse(text)`` returns an object with ``status``, ``freeze_commit``, ``a3_branch`` (strings),
  ``inputs``, ``evidence``, ``outputs`` (lists of ``dict[str, str]`` keyed by the column names),
  ``a3`` (``{"cascadia": status, "fallback": status}``) and ``steps``: objects with ``id``
  (``"S7"``), ``title``, ``issue`` (``"M10-03"``), ``job`` (``dict[str, str]`` of the ``Job:``
  fields, or ``None`` for ``Job: none (session)``), ``commands`` (the bash block bodies),
  ``reads`` and ``writes`` (lists of strings).
* Step ids are ``S<n>`` or ``S-<word>`` (the steps table names ``S-check`` and ``S-bench``).
* ``static_check`` returns failure strings; the only warning is ``duplicate_hash:<id>,<id>``,
  which the CLI prints without failing.
* ``verify`` returns ``{"failures": [...], "inputs": {id: {"status": ...}, ...}, ...}``.
* ``a3_branch=pending`` needs no ``active`` subsection and does not consult the verdict.
* ``pilot_cleanup.sbatch`` reads the freeze commit F from the runbook header, the runbook being
  ``$RUNBOOK`` (default ``$GIAP_REPO/validation/pilot/reproduce.md``), and the compare evidence
  is a JSON under ``$R/evidence/`` like the arms'.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType

from helpers.slurm_sandbox import REPO, Sandbox, make_sandbox, write_exe

PILOT = REPO / "validation" / "pilot"
REPRODUCE_MD = PILOT / "reproduce.md"
REPRODUCE_CHECK = PILOT / "reproduce_check.py"
SLURM = REPO / "scripts" / "slurm"
CLEANUP = SLURM / "pilot_cleanup.sbatch"
REPRO_SBATCH = SLURM / "pilot_reproduce_check.sbatch"

GENOME_REL = "tests/data/genomes/nissle1917.fna.gz"
GENOME_SHA = "fe74801661e688a525eef7a7a6be6210659ebc9a2e61691efb1f10b65ea623ee"
MANIFEST_REL = "validation/pilot/db_a2/makedb_manifest.json"
DIANN_SHA = "1946452d56d2fb0f66e05ced473df88e53484e28040f1ecbd83c4a43405a6381"
PIN_REL = "validation/pilot/diann_image.sha256"

F = "0123456789abcdef0123456789abcdef01234567"
F7 = F[:7]
OTHER_COMMIT = "fedcba9876543210fedcba9876543210fedcba98"

DELETED_COLUMNS = ["path", "kind", "bytes", "sha256", "action", "utc", "job"]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def load_reproduce_check() -> ModuleType:
    """Import ``validation/pilot/reproduce_check.py`` (stdlib only) as a module."""
    assert REPRODUCE_CHECK.is_file(), (
        f"{REPRODUCE_CHECK.relative_to(REPO)} does not exist (M10-14 deliverable 2)"
    )
    spec = importlib.util.spec_from_file_location("giap_pilot_reproduce_check", REPRODUCE_CHECK)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # dataclasses need the module registered
    spec.loader.exec_module(mod)
    return mod


def run_cli(*args: str, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    assert REPRODUCE_CHECK.is_file(), f"{REPRODUCE_CHECK.relative_to(REPO)} does not exist"
    full = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"), "LANG": "C.UTF-8", **(env or {})}
    return subprocess.run(
        [sys.executable, "-I", str(REPRODUCE_CHECK), *args],
        env=full,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )


# ------------------------------------------------------------------ the synthetic checkout

_HEADER = """#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
set -euo pipefail
source "${{GIAP_REPO:-/home/classes/bios/270/khoa/giap}}/scripts/farmshare_env.sh"
echo stand-in
"""

SCRIPTS = {
    "diann_pilot.sbatch": dict(name="giap-diann", cpus=16, mem="64000M", time="24:00:00"),
    "diann_convert.sbatch": dict(name="giap-diann-convert", cpus=2, mem="8G", time="00:30:00"),
    "giap_run.sbatch": dict(name="giap-run", cpus=2, mem="8G", time="2-00:00:00"),
}


def synthetic_repo(tmp: Path, *, verdict: str | None = None) -> Path:
    repo = tmp / "repo"
    slurm = repo / "scripts" / "slurm"
    slurm.mkdir(parents=True)
    for name, fields in SCRIPTS.items():
        (slurm / name).write_text(_HEADER.format(**fields))
    genome = repo / GENOME_REL
    genome.parent.mkdir(parents=True)
    shutil.copyfile(REPO / GENOME_REL, genome)
    manifest = repo / MANIFEST_REL
    manifest.parent.mkdir(parents=True)
    manifest.write_text(
        json.dumps({"inputs": {"genome": {"path": GENOME_REL, "sha256": GENOME_SHA}}})
    )
    (repo / PIN_REL).write_text(f"{DIANN_SHA}  diann_2.3.1.sif\n")
    (repo / "validation" / "pilot" / "diann_A1.cfg").write_text("[step1]\n--threads 16\n")
    if verdict is not None:
        set_verdict(repo, verdict)
    return repo


def set_verdict(repo: Path, verdict: str) -> None:
    p = repo / "validation" / "pilot" / "metrics" / "a3_feasibility.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"verdict": verdict, "reason": None}))


# ------------------------------------------------------------------ the runbook renderer

INPUT_GENOME = f"| genome | {GENOME_REL} | {GENOME_SHA} | {MANIFEST_REL} |"
INPUT_IMAGE = f"| diann_image | $GIAP_ROOT/images/diann_2.3.1.sif | {DIANN_SHA} | {PIN_REL} |"

STEP_S0 = """### S0 checkout (M10-14)

Job: none (session)

```bash
git -C "$GIAP_REPO" rev-parse HEAD
```

Reads: -
Writes: -
"""

STEP_S5 = """### S5 searches (M10-02)

Job: script=scripts/slurm/diann_pilot.sbatch name=giap-diann-A1 partition=normal qos=normal \
cpus=16 gpus=0 mem={mem} time=24:00:00 constraint=- array=-

```bash
sbatch -J {jname} {extra}"$GIAP_REPO/{script}" A1
```

Reads: genome, diann_image
Writes: $GIAP_ROOT/pilot/A1/
"""

# Spec example 1, verbatim.
STEP_S7 = """### S7 conversion (M10-03)

Job: script=scripts/slurm/diann_convert.sbatch name=giap-diann-convert-A1 partition=normal \
qos=normal cpus=2 gpus=0 mem=8G time=00:30:00 constraint=- array=-

```bash
sbatch -J giap-diann-convert-A1 "$GIAP_REPO/scripts/slurm/diann_convert.sbatch" A1
```

Reads: genome
Writes: $GIAP_ROOT/pilot/inputs/a1.tsv.gz
"""

STEP_S8 = """### S8 freeze A2 (M10-14)

Job: script=scripts/slurm/giap_run.sbatch name=giap-run-freeze-a2 partition=normal qos=normal \
cpus=2 gpus=0 mem=8G time=2-00:00:00 constraint=- array=- children=nextflow:normalx14cpu

```bash
G=$GIAP_REPO/tests/data/genomes/nissle1917.fna.gz
sbatch -J giap-run-freeze-a2 scripts/slurm/giap_run.sbatch --genome $G --outdir $R/out_a2
```

Reads: genome
Writes: $R/out_a2
"""

EVIDENCE_HEAD = (
    "| step | label | job | state | exit | elapsed | maxrss | node | commit |\n"
    "|---|---|---|---|---|---|---|---|---|"
)
EVIDENCE_S5 = f"| S5 | A1 | 1879913 | COMPLETED | 0:0 | 05:12:00 | 30G | barley-07 | {'c' * 40} |"
EVIDENCE_S8_PENDING = (
    "| S8 | A2 freeze | pending | pending | pending | pending | pending | pending | pending |"
)


def evidence_s8(commit: str = F) -> str:
    return f"| S8 | A2 freeze | 1900001 | COMPLETED | 0:0 | 06:00:00 | 7G | barley-03 | {commit} |"


def runbook(
    *,
    status: str = "draft",
    freeze: str = "pending",
    a3_branch: str = "pending",
    cascadia: str = "inactive",
    fallback: str = "inactive",
    inputs: list[str] | None = None,
    steps: list[str] | None = None,
    evidence: list[str] | None = None,
    outputs: list[str] | None = None,
    prose: str = "",
) -> str:
    """A ``reproduce.md`` in the spec's format; the defaults pass ``static`` on synthetic_repo."""
    inputs = [INPUT_GENOME, INPUT_IMAGE] if inputs is None else inputs
    steps = [STEP_S0, s5(), STEP_S7, STEP_S8] if steps is None else steps
    evidence = [EVIDENCE_S5, EVIDENCE_S8_PENDING] if evidence is None else evidence
    outputs = [] if outputs is None else outputs
    parts = [
        f"<!-- reproduce: status={status} freeze_commit={freeze} a3_branch={a3_branch} -->",
        "# Pilot reproduction runbook (synthetic test fixture)",
        "",
        "Every command runs with scripts/farmshare_env.sh sourced.",
        prose,
        "",
        "## Inputs",
        "",
        "| id | path | sha256 | recorded_in |",
        "|---|---|---|---|",
        *inputs,
        "",
        "## Steps",
        "",
        "\n".join(steps),
        "## A3",
        "",
        f"### A3-cascadia (status: {cascadia})",
        "",
        "Pending M10-10.",
        "",
        f"### A3-fallback (status: {fallback})",
        "",
        "Pending M10-11b.",
        "",
        "## Evidence",
        "",
        EVIDENCE_HEAD,
        *evidence,
        "",
        "## Outputs",
        "",
        "| id | path | sha256 | recorded_in | kept |",
        "|---|---|---|---|---|",
        *outputs,
        "",
    ]
    return "\n".join(parts)


def s5(
    *,
    mem: str = "64000M",
    jname: str = "giap-diann-A1",
    extra: str = "",
    script: str = "scripts/slurm/diann_pilot.sbatch",
) -> str:
    return STEP_S5.format(mem=mem, jname=jname, extra=extra, script=script)


# ------------------------------------------------------------------ the clean-up sandbox

_SQUEUE = r"""#!{python}
'''Fake squeue: prints the jobs of $FAKE_SLURM_DIR/squeue_jobs.json in the -o/--format given.

Each job is {"id": ..., "name": ..., "state": ...}. Supports %i %A %j %T %t %u (widths such as
%.18i are accepted); -h/--noheader drops the header line. Every call is logged.'''
import json, os, re, sys
state = os.environ["FAKE_SLURM_DIR"]
argv = sys.argv[1:]
with open(os.path.join(state, "squeue.jsonl"), "a") as fh:
    fh.write(json.dumps(argv) + "\n")
p = os.path.join(state, "squeue_jobs.json")
jobs = json.load(open(p)) if os.path.exists(p) else []
fmt, header = "%i %j %T", True
i = 0
while i < len(argv):
    a = argv[i]
    if a in ("-o", "--format", "-O") and i + 1 < len(argv):
        fmt = argv[i + 1]; i += 1
    elif a.startswith("--format="):
        fmt = a.split("=", 1)[1]
    elif a.startswith("-o") and len(a) > 2:
        fmt = a[2:]
    elif a in ("-h", "--noheader"):
        header = False
    i += 1
keys = {"i": "id", "A": "id", "j": "name", "T": "state", "t": "state", "u": "user"}
heads = {"i": "JOBID", "A": "JOBID", "j": "NAME", "T": "STATE", "t": "ST", "u": "USER"}
pat = re.compile(r"%\.?\d*([iAjTtu])")
if header:
    print(pat.sub(lambda m: heads[m.group(1)], fmt))
for job in jobs:
    job = dict(job, user=os.environ.get("USER", "giap"))
    print(pat.sub(lambda m: str(job[keys[m.group(1)]]), fmt))
"""

_NEXTFLOW = r"""#!{python}
'''Fake nextflow: logs argv and cwd; deletes nothing.'''
import json, os, sys
with open(os.path.join(os.environ["FAKE_SLURM_DIR"], "nextflow.jsonl"), "a") as fh:
    fh.write(json.dumps({"argv": sys.argv[1:], "cwd": os.getcwd()}) + "\n")
"""

# bash 3.2: no arrays beyond "$@", no [[ =~ ]] needed. Records, for every argument under the
# sandbox $GIAP_ROOT, whether DELETED.tsv already named it when rm was called; then runs rm.
_RM = r"""#!/bin/bash
ledger="$GIAP_ROOT/pilot/cleanup/DELETED.tsv"
for a in "$@"; do
  case "$a" in
    "$FAKE_ROOT"/*)
      if [ -f "$ledger" ] && grep -F -q -- "$a" "$ledger"; then
        printf 'ledgered\t%s\n' "$a" >> "$FAKE_SLURM_DIR/rm.log"
      else
        printf 'unledgered\t%s\n' "$a" >> "$FAKE_SLURM_DIR/rm.log"
      fi ;;
  esac
done
exec /bin/rm "$@"
"""


@dataclass
class CleanupSandbox:
    sb: Sandbox
    root: Path
    outside: Path
    runbook: Path
    work_dirs: list[Path]
    run_dirs: list[Path]
    delete_files: list[Path]  # every regular file CLEAN=1 must remove
    keep_files: list[Path]  # every regular file that must survive byte for byte
    k2_files: list[Path]  # DIA-NN reports/libraries: kept pending K2 (lead decision)
    mzml_dir: Path

    @property
    def cleanup_dir(self) -> Path:
        return self.root / "pilot" / "cleanup"

    def run(self, **env: str) -> subprocess.CompletedProcess[str]:
        script = self.sb.repo / "scripts" / "slurm" / "pilot_cleanup.sbatch"
        assert script.is_file(), "scripts/slurm/pilot_cleanup.sbatch does not exist (M10-14)"
        return self.sb.run(["bash", str(script)], env=env, cwd=self.sb.tmp, timeout=120)

    def set_jobs(self, jobs: list[dict[str, object]]) -> None:
        (self.sb.state / "squeue_jobs.json").write_text(json.dumps(jobs))

    def nextflow_calls(self) -> list[dict[str, object]]:
        return self.sb.jsonl("nextflow.jsonl")

    def rm_log(self) -> list[tuple[str, str]]:
        p = self.sb.state / "rm.log"
        if not p.is_file():
            return []
        return [tuple(ln.split("\t", 1)) for ln in p.read_text().splitlines()]  # type: ignore[misc]

    def snapshot(self) -> dict[str, str]:
        """sha256 of each regular file under $GIAP_ROOT and the outside tree (cleanup/ skipped)."""
        out: dict[str, str] = {}
        for base in (self.root, self.outside):
            for dirpath, _dirnames, filenames in os.walk(base):
                d = Path(dirpath)
                if d == self.cleanup_dir or self.cleanup_dir in d.parents:
                    continue
                if d == self.root / "venv" or (self.root / "venv") in d.parents:
                    continue
                if d == self.root / "logs" or (self.root / "logs") in d.parents:
                    continue
                for f in filenames:
                    p = d / f
                    if p.is_symlink() or not p.is_file():
                        continue
                    out[str(p)] = sha256_file(p)
        return out

    def ledger(self) -> list[dict[str, str]]:
        p = self.cleanup_dir / "DELETED.tsv"
        assert p.is_file(), f"no ledger at {p}"
        lines = [ln for ln in p.read_text().splitlines() if ln.strip()]
        assert lines, "DELETED.tsv is empty (it must have at least its header)"
        head = lines[0].split("\t")
        assert head == DELETED_COLUMNS, f"DELETED.tsv header {head} (spec: {DELETED_COLUMNS})"
        return [dict(zip(head, ln.split("\t"), strict=False)) for ln in lines[1:]]

    def repo_status(self) -> str:
        """``git status --porcelain`` of the sandbox checkout; ``__pycache__`` is ignored there
        as in the real repository's .gitignore (the sandbox has none)."""
        cmd = ["git", "-C", str(self.sb.repo), "status", "--porcelain"]
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        return "\n".join(ln for ln in out.splitlines() if "__pycache__" not in ln).strip()


def _w(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def _run_dir(base: Path, *, nextflow: bool = True) -> tuple[Path, list[Path], list[Path]]:
    """A GIAP outdir: (.nextflow, work/ to delete, tables/logs/manifest/report to keep)."""
    keep = [
        _w(base / "tables" / "proteins.tsv", "locus_tag\tclass\nX_0001\tdetected\n"),
        _w(base / "mapping" / "peptides.tsv", "peptide\nAAAK\n"),
        _w(base / "logs" / "giap.log", "ok\n"),
        _w(base / "run_manifest.json", json.dumps({"git_commit": F, "git_dirty": False})),
        _w(base / "report" / "report.html", "<html></html>\n"),
        _w(base / ".nextflow.log", "nextflow log\n"),
    ]
    delete: list[Path] = []
    if nextflow:
        keep.append(_w(base / ".nextflow" / "history", "session\n"))
        delete = [
            _w(base / "work" / "ab" / "cdef01" / ".command.sh", "echo task\n"),
            _w(base / "work" / "ab" / "cdef01" / "proteins.tsv", "intermediate\n"),
            _w(base / "work" / "12" / "345678" / ".exitcode", "0"),
        ]
    else:
        keep.append(_w(base / "work" / "notes.txt", "not a nextflow run\n"))
    return base, delete, keep


def cleanup_runbook(
    extra_inputs: list[str] | None = None, extra_outputs: list[str] | None = None
) -> str:
    return runbook(
        status="final",
        freeze=F,
        a3_branch="fallback",
        fallback="active",
        inputs=[INPUT_GENOME, *(extra_inputs or [])],
        steps=[STEP_S0],
        evidence=[evidence_s8()],
        outputs=extra_outputs or [],
    )


def cleanup_sandbox(tmp: Path) -> CleanupSandbox:
    tmp = tmp.resolve()
    sb = make_sandbox(tmp)
    shutil.copytree(
        REPO / "validation",
        sb.repo / "validation",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    rb = sb.repo / "validation" / "pilot" / "reproduce.md"
    rb.write_text(cleanup_runbook())
    for args in (["add", "-A"], ["commit", "-q", "-m", "validation + runbook"]):
        subprocess.run(
            ["git", "-C", str(sb.repo), *args], check=True, capture_output=True, env=sb.env
        )
    py = sys.executable
    write_exe(sb.bin / "squeue", _SQUEUE.replace("{python}", py))
    write_exe(sb.bin / "nextflow", _NEXTFLOW.replace("{python}", py))
    write_exe(sb.bin / "rm", _RM)
    root = sb.root
    outside = tmp / "outside"
    outside.mkdir()
    sb.env.update(
        {
            "GIAP_PYTHON": py,
            "SLURM_JOB_ID": "4242",
            "SLURM_JOB_NAME": "giap-pilot-cleanup",
            "FAKE_ROOT": str(root),
            "TMPDIR": str(tmp / "tmpdir"),
        }
    )
    (tmp / "tmpdir").mkdir()

    work_dirs: list[Path] = []
    run_dirs: list[Path] = []
    delete: list[Path] = []
    keep: list[Path] = []
    R = root / "runs" / "pilot" / "freeze" / F7
    for base in (
        R / "out_a2",
        R / "out_a1",
        root / "runs" / "pilot" / "pass1" / "out_a1",
        root / "runs" / "pilot" / "te7" / "rep1",
        root / "bench" / "b1" / "r1" / "standard" / "out",
    ):
        d, dl, kp = _run_dir(base)
        run_dirs.append(d)
        work_dirs.append(d / "work")
        delete += dl
        keep += kp
    _, _, kp = _run_dir(root / "runs" / "pilot" / "pass1" / "notes", nextflow=False)
    keep += kp

    # freeze evidence (the preconditions)
    keep += [
        _w(R / "evidence" / "A2.json", json.dumps({"arm": "A2", "failures": []})),
        _w(R / "evidence" / "A1.json", json.dumps({"arm": "A1", "failures": []})),
        _w(R / "evidence" / "compare.json", json.dumps({"failures": []})),
        _w(R / "ATTEMPTS.tsv", "arm\tattempt\tjob\nA2\t1\t1900001\n"),
        _w(root / "pilot" / "metrics" / f"freeze-{F7}" / "metrics.json", "{}\n"),
        _w(
            root / "pilot" / "reproduce" / f"verify-{F7}.json",
            json.dumps(
                {
                    "freeze_commit": F,
                    "static_ok": True,
                    "failures": [],
                    "reconvert": {"a1": "identical", "a2": "identical"},
                }
            ),
        ),
    ]

    # DIA-NN directories: plan step 2 (kept pending K2) and step 3 (kept)
    k2: list[Path] = []
    for arm, d in (
        ("A1", root / "pilot" / "A1"),
        ("A2", root / "pilot" / "A2"),
        ("A2", root / "pilot" / "te7" / "A2"),
    ):
        k2 += [
            _w(d / f"{arm}.speclib", "library\n"),
            _w(d / f"{arm}.predicted.speclib", "predicted library\n"),
            _w(d / f"{arm}.parquet", "step-1 library\n"),
            _w(d / f"{arm}_empirical.parquet", "empirical\n"),
            _w(d / f"{arm}_report.parquet", "report\n"),
            _w(d / f"{arm}_report.pg_matrix.tsv", "matrix\n"),
            _w(d / "quant" / "WT_1.quant", "quant\n"),
        ]
        keep += [
            _w(d / "provenance.json", json.dumps({"arm": arm})),
            _w(d / "outputs.sha256", "x  y\n"),
            _w(d / "command_step1.txt", "diann --lib\n"),
            _w(d / "command_step2.txt", "diann --f\n"),
            _w(d / "check.json", json.dumps({"failures": []})),
            _w(d / "step1.done", "done\n"),
            _w(d / "diann_version.txt", "DIA-NN 2.3.1\n"),
            _w(d / "fasta" / "db.fasta", ">p\nMAAAK\n"),
        ]
    # plan step 4: mzML (deleted); plan steps 5 and 6: kept
    mzml = root / "data" / "pxd078459" / "mzml"
    delete += [_w(mzml / "WT_1.mzML", "<mzML/>\n"), _w(mzml / "WT_2.mzML", "<mzML/>\n")]
    keep += [
        _w(root / "pilot" / "recorded" / "diann" / "2.3.1" / "A1" / "report.sample.parquet", "s\n"),
        _w(root / "pilot" / "inputs" / "a1.tsv.gz", "converted\n"),
        _w(root / "pilot" / "inputs" / "a1.convert.json", "{}\n"),
        _w(root / "data" / "pxd078459" / "WT_1.d" / "analysis.tdf", "tdf\n"),
        _w(root / "data" / "pxd078459" / "WT_1.d" / "analysis.tdf_bin", "tdf_bin\n"),
        _w(root / "data" / "pxd078459" / "WT_1.zip", "archive\n"),
    ]
    return CleanupSandbox(
        sb=sb,
        root=root,
        outside=outside,
        runbook=rb,
        work_dirs=work_dirs,
        run_dirs=run_dirs,
        delete_files=delete,
        keep_files=keep + k2,
        k2_files=k2,
        mzml_dir=mzml,
    )


def commit_runbook(cs: CleanupSandbox, text: str) -> None:
    cs.runbook.write_text(text)
    for args in (["add", "-A"], ["commit", "-q", "-m", "runbook"]):
        subprocess.run(
            ["git", "-C", str(cs.sb.repo), *args], check=True, capture_output=True, env=cs.sb.env
        )


# ------------------------------------------------------------------ sbatch house rules


def sbatch_header(path: Path) -> dict[str, str]:
    out = {}
    for m in re.finditer(r"^#SBATCH\s+--([\w-]+)(?:[= ](\S+))?", path.read_text(), re.M):
        out[m.group(1)] = m.group(2) or ""
    return out


def sbatch_code(path: Path) -> str:
    assert path.is_file(), f"{path.relative_to(REPO)} does not exist (M10-14 deliverable)"
    return "\n".join(
        ln
        for ln in path.read_text().splitlines()
        if not ln.lstrip().startswith("#") or ln.startswith("#SBATCH")
    )


def mem_mb(value: str) -> float:
    m = re.fullmatch(r"(\d+(?:\.\d+)?)([KMGT]?)B?", value.strip(), re.I)
    assert m, value
    return (
        float(m.group(1))
        * {"": 1, "K": 1 / 1024, "M": 1, "G": 1024, "T": 1024**2}[m.group(2).upper()]
    )


def time_minutes(value: str) -> int:
    """Slurm --time (MM, MM:SS, HH:MM:SS, D-HH, D-HH:MM, D-HH:MM:SS) in whole minutes."""
    days, rest = value.split("-", 1) if "-" in value else ("0", value)
    parts = [int(x) for x in rest.split(":")]
    if "-" in value:
        parts += [0] * (3 - len(parts))
        h, mi, sec = parts
    elif len(parts) == 3:
        h, mi, sec = parts
    elif len(parts) == 2:
        h, (mi, sec) = 0, parts
    else:
        h, mi, sec = 0, parts[0], 0
    return int(days) * 1440 + h * 60 + mi + (sec > 0)


def check_house_rules(path: Path, *, name: str, cpus: int, mem_gb: int, minutes: int) -> None:
    """Criterion 4 (as test_sbatch_house_rules checks) plus the steps-table job class."""
    code = sbatch_code(path)
    hdr = sbatch_header(path)
    assert hdr.get("job-name") == name, hdr
    assert hdr.get("partition") == "normal" and hdr.get("qos") == "normal", hdr
    assert int(hdr.get("cpus-per-task", "0")) == cpus, hdr
    assert 1000 * mem_gb <= mem_mb(hdr.get("mem", "0")) <= 1024 * mem_gb, hdr
    assert time_minutes(hdr.get("time", "0")) == minutes, hdr
    out = hdr.get("output", "")
    assert "/local/logs/" in out, f"--output {out!r} not under local/logs/"
    assert "gres" not in hdr or "gpu" not in hdr["gres"], "R14: no GPU for these jobs"
    assert "scripts/farmshare_env.sh" in code, "the environment is not sourced"
    assert re.search(r"^set -euo pipefail$", code, re.M), "no `set -euo pipefail`"
    assert not re.search(r"\buv\s+(run|sync)\b|\bpip3?\s+install\b|-m pip\b", code), "uv/pip"
    into_repo = re.compile(
        r"(>>?|\btee\b|\bcp\b.*|\bmv\b.*|--out(dir)?[= ]|mkdir -p)\s*\"?\$\{?GIAP_REPO\b"
    )
    bad = [ln for ln in code.splitlines() if into_repo.search(ln)]
    assert not bad, f"writes under $GIAP_REPO: {bad}"
