"""M10-02 sandbox additions on top of ``helpers.pilot_kit`` (docs/specs/M10-02-diann-searches.md).

``helpers.pilot_kit`` is the M3-24 sandbox and is amended by the lead (M10-02 lead decision 1);
this module spells out the additions the M10-02 spec lists under "Helper changes" without editing
it, so the M10-02 tests can run before and after the lead's amendment:

* the stand-in image (``stub sif\\n``) is made mode 0444, and the sandbox checkout's
  ``validation/pilot/diann_image.sha256`` holds its sha256 (so the sandbox pin matches it);
* ``validation/pilot/inputs/acquisition.json`` (M10-01's shape: ``dataset_mode`` and one object
  per run with ``run``, ``mode`` and ``notes``) names ``wt_rep1``-``wt_rep3``, and every run gets
  an ``analysis.tdf_bin`` stand-in;
* ``GIAP_PYTHON`` is ``sys.executable`` (it has pyarrow);
* ``DIANN_SIF_SHARED`` is a second stand-in with the same bytes, mode 0777 like Khoa's image;
* the fake ``apptainer`` is replaced by an extended one (below), and ``cp`` is wrapped so the
  install tests can see which files were copied.

Both edits to the sandbox repository are committed there, so ``repo_is_clean`` keeps meaning
"the script wrote nothing into the checkout".

The extended fake ``apptainer``/``singularity`` never executes anything. Environment knobs:

* ``FAKE_DIANN_BANNER``: the banner line, printed after a blank line as the real binary does
  (default ``DIA-NN 2.3.1 Academia (…)``);
* ``FAKE_DIANN_RC``: exit code of every real (non-smoke) DIA-NN call;
* ``FAKE_DIANN_RC_STEP1`` / ``FAKE_DIANN_RC_STEP2``: exit code of the predict call (it has
  ``--predictor``) or of the search call (it has ``--lib`` or ``--f``) only;
* ``FAKE_DIANN_NO_REPORT=1``: the search call exits 0 but writes no ``--out`` file;
* ``FAKE_DIANN_NO_LIB=1``: the predict call exits 0 but writes no ``.predicted.speclib``;
* ``FAKE_DIANN_SLEEP`` / ``FAKE_DIANN_SLEEP_STEP1`` / ``FAKE_DIANN_SLEEP_STEP2``: seconds a real
  call (or only the predict / the search call) sleeps before it writes anything; while it sleeps
  ``$FAKE_SLURM_DIR/diann.sleeping`` holds ``step1``/``step2`` (``wait_for_sleep``), so a test can
  signal the job mid-search (Slurm TIMEOUT, scancel).

The image of an ``exec`` is the first non-option token after ``exec`` (``-B``/``--bind`` take a
value), whatever its name, so a banner check on a temporary copy is seen; ``image_sha256`` is
the sha256 of that file when the call starts.

A search call with ``--f`` runs writes ``<basename of run>.quant`` into its ``--temp`` directory,
or next to the run when no ``--temp`` is given (DIA-NN's default); ``temp`` and
``temp_listing`` (the ``--temp`` directory's file names when the call starts) are logged.

Each call is logged to ``$FAKE_SLURM_DIR/diann.jsonl`` with the fields of the M3-24 stub plus
``inspect`` (an ``apptainer inspect``) and ``preexisting``: for each of ``--out``, ``--out-lib``,
``--lib`` and the ``.predicted.speclib`` and ``.log.txt`` siblings of ``--out-lib``, whether it
existed when the call started (so a test can tell that a leftover file was deleted *before*
the call).
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from helpers.gnu_userland import require_gnu_userland
from helpers.pilot_kit import (
    DIANN_BIN,
    DIANN_SIF,
    PILOT,
    README,
    SLURM,
    PilotSandbox,
    pilot_sandbox,
)
from helpers.slurm_sandbox import REPO, write_exe

# ------------------------------------------------------------------ paths and constants

PIN_FILE = PILOT / "diann_image.sha256"
PIN_NAME = "diann_2.3.1.sif"
DIANN_REPORT = PILOT / "diann_report.py"
ACQUISITION = PILOT / "inputs" / "acquisition.json"
RECORDED_DIANN = REPO / "tests" / "data" / "recorded" / "diann"
M10_SCRIPTS = {
    "diann_image": SLURM / "diann_image.sbatch",
    "diann_pilot": SLURM / "diann_pilot.sbatch",
    "diann_report": SLURM / "diann_report.sbatch",
}
DIANN_SIF_SHARED = DIANN_SIF  # Khoa's image; only diann_image.sbatch may read it
COPY_REL = Path("images") / PIN_NAME  # under $GIAP_ROOT (index lead decision 2)
STAND_IN = b"stub sif\n"
STAND_IN_SHA256 = hashlib.sha256(STAND_IN).hexdigest()
RUNS = ("wt_rep1", "wt_rep2", "wt_rep3")
BANNER = "DIA-NN 2.3.1 Academia"
HEX64 = r"[0-9a-f]{64}"

# The columns the pilot depends on (check rule 2) and the converter's optional ones.
REPORT_REQUIRED = (
    "Stripped.Sequence",
    "Modified.Sequence",
    "Protein.Ids",
    "Protein.Group",
    "Run",
    "Precursor.Charge",
    "Global.Q.Value",
    "Decoy",
    "PG.MaxLFQ",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ------------------------------------------------------------------ acquisition.json fixtures


def acquisition_doc(
    runs: Sequence[str] = RUNS,
    *,
    mode: str = "diaPASEF",
    notes: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """An ``acquisition.json`` in M10-01's shape (``dataset_mode`` plus one object per run)."""
    notes = dict(notes or {})
    return {
        "dataset_mode": mode,
        "runs": [
            {
                "run": r,
                "mode": mode
                if mode in {"diaPASEF", "ddaPASEF", "ms1_only", "mixed"}
                else "unknown",
                "incomplete_run": False,
                "notes": notes.get(r, ""),
                "source_folder": f"{r}.d",
            }
            for r in runs
        ],
    }


def write_json(path: Path, doc: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, indent=2) + "\n")
    return path


# ------------------------------------------------------------------ the extended fake container

_APPTAINER = r"""#!{python}
'''Extended fake apptainer/singularity (tests/helpers/pilot_kit_m10.py): never executes anything.'''
import hashlib, json, os, sys, time

state = os.environ["FAKE_SLURM_DIR"]
argv = sys.argv[1:]
rec = {"tool": os.path.basename(sys.argv[0]), "argv": argv, "cwd": os.getcwd(),
       "binds": [], "image": None, "inner": [], "created": [], "smoke": False,
       "inspect": False, "preexisting": {}, "image_sha256": None, "temp": None,
       "temp_listing": None}

def log(rc):
    rec["rc"] = rc
    with open(os.path.join(state, "diann.jsonl"), "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    sys.exit(rc)

banner = os.environ.get("FAKE_DIANN_BANNER",
                        "DIA-NN 2.3.1 Academia (Data-Independent Acquisition by Neural Networks)")

if argv and argv[0] == "inspect":
    rec["inspect"] = rec["smoke"] = True
    rec["image"] = next((a for a in argv[1:] if not a.startswith("-")), None)
    if rec["image"] is None or not os.path.exists(rec["image"]):
        print("FATAL: image not found", file=sys.stderr)
        log(255)
    labels = {"org.label-schema.build-date": "Tuesday_6_January_2026_19:18:54_PST",
              "org.label-schema.usage.singularity.deffile.bootstrap": "docker",
              "org.label-schema.usage.singularity.deffile.from": "khoayelo189/diann:2.3.1",
              "org.label-schema.usage.apptainer.version": "1.4.4"}
    if "--json" in argv:
        print(json.dumps({"data": {"attributes": {"labels": labels}}}))
    else:
        for k, v in labels.items():
            print(k + ": " + v)
    log(0)

i = 1 if argv and argv[0] in ("exec", "run") else 0
while i < len(argv):
    a = argv[i]
    if a in ("-B", "--bind"):
        rec["binds"] += argv[i + 1].split(","); i += 2; continue
    if a.startswith("--bind="):
        rec["binds"] += a.split("=", 1)[1].split(","); i += 1; continue
    if a.startswith("-"):
        i += 1; continue
    rec["image"] = a
    rec["inner"] = argv[i + 1:]
    break

if rec["image"] is None:
    rec["smoke"] = True
    print("apptainer version 1.4.4")
    log(0)

if not os.path.exists(rec["image"]):
    print("FATAL: image not found: " + rec["image"], file=sys.stderr)
    log(255)
try:
    with open(rec["image"], "rb") as fh:
        rec["image_sha256"] = hashlib.sha256(fh.read()).hexdigest()
except OSError:
    pass

inner = rec["inner"][1:]
if not inner or set(inner) <= {"--help", "--version", "--dry-run"}:
    rec["smoke"] = True
    print()  # the real binary prints a blank line first (job 1879883)
    print(banner)
    print("Compiled on Dec  5 2025")
    if "--version" in inner:
        print("WARNING: unrecognised option [--version]")
    log(0)

def value(flag):
    return inner[inner.index(flag) + 1] if flag in inner and inner.index(flag) + 1 < len(inner) \
        else None

out, out_lib, lib = value("--out"), value("--out-lib"), value("--lib")
watch = [p for p in (out, out_lib, lib) if p]
if out_lib:
    watch.append(os.path.splitext(out_lib)[0] + ".predicted.speclib")
    watch.append(os.path.splitext(out_lib)[0] + ".log.txt")
rec["preexisting"] = {p: os.path.exists(p) for p in watch}
temp = value("--temp")
rec["temp"] = temp
if temp is not None:
    rec["temp_listing"] = sorted(os.listdir(temp)) if os.path.isdir(temp) else None
runs = [inner[k + 1] for k, a in enumerate(inner) if a == "--f" and k + 1 < len(inner)]

def touch(path):
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("stub DIA-NN output for " + os.path.basename(path) + "\n")
    rec["created"].append(path)

if lib is not None and not os.path.exists(lib):
    print("FATAL: library not found: " + lib, file=sys.stderr)
    log(1)

search = lib is not None or "--f" in inner
rc = int(os.environ.get("FAKE_DIANN_RC", "0"))
step_rc = os.environ.get("FAKE_DIANN_RC_STEP2" if search else "FAKE_DIANN_RC_STEP1")
if step_rc:
    rc = int(step_rc)
nap = os.environ.get("FAKE_DIANN_SLEEP_STEP2" if search else "FAKE_DIANN_SLEEP_STEP1") \
    or os.environ.get("FAKE_DIANN_SLEEP")
if nap:
    with open(os.path.join(state, "diann.sleeping"), "w") as fh:
        fh.write("step2" if search else "step1")
    time.sleep(float(nap))
if rc:
    print("DIA-NN: fake failure", file=sys.stderr)
    log(rc)

for r in runs:
    q = os.path.basename(r.rstrip("/")) + ".quant"
    touch(os.path.join(temp, q) if temp is not None else r.rstrip("/") + ".quant")
if out_lib is not None:
    touch(out_lib)
    if search or not os.environ.get("FAKE_DIANN_NO_LIB"):
        touch(os.path.splitext(out_lib)[0] + ".predicted.speclib")
if out is not None and not os.environ.get("FAKE_DIANN_NO_REPORT"):
    touch(out)
    touch(os.path.splitext(out)[0] + ".stats.tsv")
log(0)
"""

_CP = r"""#!/bin/bash
# Logs every cp call (one line, tab-separated argv) and runs the real cp.
( IFS=$'\t'; printf '%s\n' "$*" ) >> "$FAKE_SLURM_DIR/cp.log"
exec /bin/cp "$@"
"""

_GZIP = r"""#!/bin/bash
# gzip/gunzip wrapper (tests/helpers/pilot_kit_m10.py): FAKE_GZIP_FAIL=1 is a gunzip killed partway.
printf '%s\n' "$(basename "$0") $*" >> "$FAKE_SLURM_DIR/gzip.log"
if [ -n "${FAKE_GZIP_FAIL:-}" ]; then
  printf '>truncated_protein\nMKK\n'
  exit 1
fi
exec "{real}" "$@"
"""


# ------------------------------------------------------------------ the sandbox


@dataclass
class M10Sandbox:
    pilot: PilotSandbox
    shared: Path  # stand-in for Khoa's world-writable image (DIANN_SIF_SHARED)
    stand_in: Path  # the read-only stand-in DIANN_SIF points at (mode 0444)

    @property
    def sb(self):  # type: ignore[no-untyped-def]
        return self.pilot.sb

    @property
    def root(self) -> Path:
        return self.pilot.sb.root

    @property
    def repo(self) -> Path:
        return self.pilot.sb.repo

    @property
    def copy(self) -> Path:
        """Where ``diann_image.sbatch`` installs the copy by default."""
        return self.root / COPY_REL

    @property
    def raw_dir(self) -> Path:
        return self.pilot.raw_dir

    def out(self, arm: str) -> Path:
        return self.pilot.out_dir / arm

    @property
    def pin_file(self) -> Path:
        return self.repo / "validation" / "pilot" / "diann_image.sha256"

    @property
    def acquisition(self) -> Path:
        return self.repo / "validation" / "pilot" / "inputs" / "acquisition.json"

    def run(self, script: str, *args: str, unset: Iterable[str] = (), **env: str):  # type: ignore[no-untyped-def]
        path = self.repo / "scripts" / "slurm" / f"{script}.sbatch"
        return self.sb.run(
            ["bash", str(path), *args], env=env, cwd=self.sb.tmp, unset=tuple(unset), timeout=300
        )

    def popen(self, script: str, *args: str, **env: str) -> subprocess.Popen[str]:
        """Start ``scripts/slurm/<script>.sbatch`` in its own session (a process group, like a
        Slurm job step), so a test can signal the whole job as ``scancel``/TIMEOUT does."""
        require_gnu_userland()  # the scripts are FarmShare-only (spec G19)
        path = self.repo / "scripts" / "slurm" / f"{script}.sbatch"
        return subprocess.Popen(
            ["bash", str(path), *args],
            env={**self.sb.env, **env},
            cwd=str(self.sb.tmp),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )

    def wait_for_sleep(self, proc: subprocess.Popen[str], timeout: float = 120) -> str:
        """Wait until the fake DIA-NN sleeps (``FAKE_DIANN_SLEEP*``); return ``step1``/``step2``."""
        flag = self.sb.state / "diann.sleeping"
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            if flag.is_file() and flag.read_text().strip():
                return flag.read_text().strip()
            if proc.poll() is not None:
                out, err = proc.communicate()
                raise AssertionError(f"the job ended (exit {proc.returncode}) before DIA-NN "
                                     f"slept\n{out}\n{err}")  # fmt: skip
            time.sleep(0.05)
        proc.kill()
        raise AssertionError("the fake DIA-NN never started sleeping")

    def install_gzip_wrappers(self) -> None:
        """``gzip``/``gunzip`` on the sandbox PATH that run the real tools, except that with
        ``FAKE_GZIP_FAIL=1`` they print a truncated FASTA and exit 1 (a gunzip killed partway)."""
        path = os.pathsep.join(d for d in self.sb.env["PATH"].split(os.pathsep)
                               if d != str(self.sb.bin))  # fmt: skip
        for name in ("gzip", "gunzip"):
            real = shutil.which(name, path=path)
            assert real, f"no real {name} on {path}"
            write_exe(self.sb.bin / name, _GZIP.replace("{real}", real))

    def calls(self) -> list[dict[str, Any]]:
        """Every container call, smoke checks and ``inspect`` included."""
        return self.sb.jsonl("diann.jsonl")

    def searches(self) -> list[dict[str, Any]]:
        return [c for c in self.calls() if not c["smoke"]]

    def reset_calls(self) -> None:
        p = self.sb.state / "diann.jsonl"
        if p.exists():
            p.unlink()

    def cp_calls(self) -> list[list[str]]:
        p = self.sb.state / "cp.log"
        return [ln.split("\t") for ln in p.read_text().splitlines()] if p.is_file() else []

    def provenance(self, arm: str) -> dict[str, Any]:
        p = self.out(arm) / "provenance.json"
        assert p.is_file(), f"{p} was not written (provenance on every exit, spec rule 7)"
        doc: dict[str, Any] = json.loads(p.read_text())
        return doc

    def repo_is_clean(self) -> str:
        return self.pilot.repo_is_clean()

    def commit(self, message: str) -> None:
        for args in (["add", "-A"], ["commit", "-q", "--allow-empty", "-m", message]):
            subprocess.run(["git", "-C", str(self.repo), *args], check=True, capture_output=True,
                           env=self.sb.env)  # fmt: skip


def m10_sandbox(tmp: Path, *, raws: int = 3, tdf_bin: bool = True) -> M10Sandbox:
    """``pilot_sandbox`` plus the M10-02 helper changes (module docstring)."""
    pilot = pilot_sandbox(tmp, raws=raws, m10_02=False)
    sb = pilot.sb
    py = sys.executable
    for name in ("apptainer", "singularity"):
        write_exe(sb.bin / name, _APPTAINER.replace("{python}", py))
    write_exe(sb.bin / "cp", _CP)

    stand_in = pilot.sif
    stand_in.write_bytes(STAND_IN)
    stand_in.chmod(0o444)
    shared = sb.tmp / "shared" / PIN_NAME
    shared.parent.mkdir(parents=True, exist_ok=True)
    shared.write_bytes(STAND_IN)
    shared.chmod(0o777)

    if tdf_bin:
        for d in pilot.raw_dir.glob("*.d"):
            (d / "analysis.tdf_bin").write_text("stub timsTOF binary\n")

    m10 = M10Sandbox(pilot=pilot, shared=shared, stand_in=stand_in)
    m10.pin_file.write_text(f"{STAND_IN_SHA256}  {PIN_NAME}\n")
    write_json(m10.acquisition, acquisition_doc(RUNS[:raws] if raws else RUNS))
    m10.commit("M10-02 sandbox: pin and acquisition.json")
    sb.env.update({"GIAP_PYTHON": py, "DIANN_SIF_SHARED": str(shared)})
    return m10


# ------------------------------------------------------------------ provenance accessors


def step_state(prov: Mapping[str, Any], step: str) -> tuple[str, int | None]:
    """``(ran|skipped|not_reached, exit code)`` of ``step1``/``step2`` in ``provenance.json``.

    Accepted shapes: ``"step1": {"status"|"state": …, "exit_code": …}`` or
    ``"step1": "ran"`` with ``"step1_exit_code": …`` next to it.
    """
    v = prov[step]
    if isinstance(v, Mapping):
        state = v.get("status", v.get("state"))
        code = v.get("exit_code", v.get("rc"))
    else:
        state, code = v, prov.get(f"{step}_exit_code")
    assert state in {"ran", "skipped", "not_reached"}, (step, v)
    return str(state), (None if code is None else int(code))


def commands(prov: Mapping[str, Any]) -> list[list[str]]:
    """The argv arrays of ``provenance.json`` ``commands`` (a list or a step-keyed object)."""
    c = prov["commands"]
    if isinstance(c, Mapping):
        return [list(c[k]) for k in sorted(c) if c[k] is not None]
    return [list(x) for x in c if x is not None]


def strings(doc: Any) -> list[str]:
    """Every string anywhere inside a JSON value."""
    if isinstance(doc, str):
        return [doc]
    if isinstance(doc, Mapping):
        return [s for v in doc.values() for s in strings(v)]
    if isinstance(doc, list):
        return [s for v in doc for s in strings(v)]
    return []


def ends_with(tokens: Sequence[str], tail: Sequence[str]) -> bool:
    return len(tokens) >= len(tail) and list(tokens[len(tokens) - len(tail) :]) == list(tail)


def diann_args(call: Mapping[str, Any]) -> list[str]:
    inner = list(call["inner"])
    assert inner and inner[0] == DIANN_BIN, inner[:1]
    return inner[1:]


def flag_values(argv: Sequence[str], flag: str) -> list[str]:
    return [argv[i + 1] for i, a in enumerate(argv) if a == flag and i + 1 < len(argv)]


# ------------------------------------------------------------------ DIA-NN report fixtures


def report_rows(
    targets: Mapping[str, int], decoys: Mapping[str, int], *, interleave: bool = True
) -> list[dict[str, Any]]:
    """Rows of a DIA-NN-like report: ``targets[run]`` target and ``decoys[run]`` decoy rows.

    With ``interleave`` the runs alternate in file order (as when DIA-NN writes by precursor).
    ``Fixture.Row`` is the row's file-order index, so a sample can be traced back.
    """
    per_run: dict[str, list[tuple[int, int]]] = {}
    for run in dict.fromkeys([*targets, *decoys]):
        nt, nd = targets.get(run, 0), decoys.get(run, 0)
        per_run[run] = [(j, 0) for j in range(nt)] + [(nt + j, 1) for j in range(nd)]
    order: list[tuple[str, int, int]] = []
    if interleave:
        width = max((len(v) for v in per_run.values()), default=0)
        for j in range(width):
            for run, items in per_run.items():
                if j < len(items):
                    order.append((run, *items[j]))
    else:
        order = [(run, j, d) for run, items in per_run.items() for j, d in items]
    rows = []
    for i, (run, j, decoy) in enumerate(order):
        seq = "PEPTIDE" + "ACDEFGHIKLMNPQRSTVWY"[j % 20] * (1 + j % 5) + "K"
        q = round(((j * 7919) % 1000) / 1000 * (0.05 if not decoy else 1.0), 6)
        rows.append({
            "Run": run,
            "Precursor.Id": f"{seq}{2 + j % 3}#{j}",
            "Stripped.Sequence": seq,
            "Modified.Sequence": seq,
            "Precursor.Charge": 2 + j % 3,
            "Protein.Ids": f"P{j % 50:05d}",
            "Protein.Group": f"P{j % 50:05d}",
            "Q.Value": q,
            "Global.Q.Value": q,
            "PEP": q,
            "Decoy": decoy,
            "PG.MaxLFQ": float(1000 + j),
            "Precursor.Normalised": float(500 + j),
            "Fixture.Row": i,
        })  # fmt: skip
    return rows


def write_report(
    path: Path, rows: Sequence[Mapping[str, Any]], *, drop: Iterable[str] = ()
) -> Path:
    import pyarrow as pa
    import pyarrow.parquet as pq

    drop = set(drop)
    cols = [c for c in rows[0] if c not in drop] if rows else list(REPORT_REQUIRED)
    table = (
        pa.table({c: [r[c] for r in rows] for c in cols})
        if rows
        else pa.table({c: pa.array([], pa.string()) for c in cols})
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, path)
    return path


def write_sums(sums: Path, *files: Path) -> Path:
    sums.parent.mkdir(parents=True, exist_ok=True)
    sums.write_text("".join(f"{sha256(f)}  {f.name}\n" for f in sorted(files)))
    return sums


def run_report(*args: str, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    """``python -I validation/pilot/diann_report.py …`` from this checkout (no sandbox needed)."""
    assert DIANN_REPORT.is_file(), f"{DIANN_REPORT.relative_to(REPO)} does not exist (M10-02)"
    env = {k: v for k, v in os.environ.items() if not k.startswith(("http", "HTTP"))}
    return subprocess.run(
        [sys.executable, "-I", str(DIANN_REPORT), *args],
        capture_output=True,
        text=True,
        timeout=300,
        cwd=str(cwd or REPO),
        env=env,
        check=False,
    )


def load_diann_report_module():  # type: ignore[no-untyped-def]
    """Import ``validation/pilot/diann_report.py`` as a module (for its pure functions)."""
    import importlib.util

    assert DIANN_REPORT.is_file(), f"{DIANN_REPORT.relative_to(REPO)} does not exist (M10-02)"
    spec = importlib.util.spec_from_file_location("giap_pilot_diann_report", DIANN_REPORT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def readme_text() -> str:
    return README.read_text()
