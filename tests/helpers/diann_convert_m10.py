"""Hand-made DIA-NN 2.3.1-shaped reports and a job sandbox for the M10-03 tests.

Spec: docs/specs/M10-03-diann-convert.md. Khoa's rule (2026-10-23): nothing derived from the real
DIA-NN reports is committed. Every row below is invented from the DIA-NN 2.x documentation and
the spec's examples; a fixture names only the columns it is about plus the converter's required
ones (criterion 8, Q-K2). The column types are set explicitly (``pa.float32()``, ``pa.int8()``,
``pa.bool_()``, ``pa.dictionary(pa.int32(), pa.string())``), as the documentation allows them.

The job sandbox (``convert_sandbox``) builds on ``helpers.slurm_sandbox.make_sandbox``, as the
M10-02 helper ``pilot_kit_m10`` does, and adds:

* ``src/``, ``validation/`` and ``tests/contract/test_m10_03_diann_contract.py`` (when present)
  in the sandbox checkout, committed, so ``PYTHONPATH=$GIAP_REPO/src`` runs this branch's
  converter and ``repo_is_clean`` means "the script wrote nothing into the checkout";
* ``$UV_PROJECT_ENVIRONMENT/bin/python`` (the default ``GIAP_PYTHON``) replaced by a recording
  wrapper that runs the real interpreter. It logs one JSON line per call to
  ``$FAKE_SLURM_DIR/python.jsonl`` (``argv``, ``kind``, ``pythonpath``, ``rc`` and ``present``:
  the names under ``$GIAP_ROOT/pilot/inputs`` when the call started). ``kind`` is ``convert``
  (``-c … giap.cli … convert``), ``summary`` (``diann_convert.py``), ``pytest`` (``-m pytest``,
  faked: ``FAKE_PYTEST=fail`` exits 1; ``FAKE_PYTEST=skip`` exits 0 with every test skipped,
  printing ``SKIPPED [n] …`` lines and ``10 skipped in …`` and writing a junit XML with
  ``skipped == tests`` when ``--junitxml`` is given; the default passes, 5 tests) or ``other``.
  A pytest record also holds ``env``: ``GIAP_M10_03_REPORT`` and ``GIAP_M10_03_ARM`` as the call
  saw them. ``FAKE_CONVERT_NONDET=<flag>`` (``--out`` or ``--abundance-out``; ``1`` means
  ``--out``) appends one byte to that output of the **second** successful convert call (the
  determinism stub, Q8). ``FAKE_CONVERT_RC=<n>`` makes the **first** convert call exit ``n``
  without converting;
* ``$UV_PROJECT_ENVIRONMENT/bin/mv``: ``/bin/mv``, except that ``FAKE_MV_FAIL=1`` makes every
  ``mv`` with an argument containing ``.tmp.`` (the publishing move) fail with exit 1;
* ``$UV_PROJECT_ENVIRONMENT/bin/giap`` runs this branch's CLI (for ``giap version --json``).

Nothing here touches the network, Slurm or the real ``$GIAP_ROOT``.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from helpers.slurm_sandbox import REPO, Sandbox, make_sandbox, write_exe

SCRIPT = REPO / "scripts" / "slurm" / "diann_convert.sbatch"
SUMMARY = REPO / "validation" / "pilot" / "diann_convert.py"
CONTRACT = REPO / "tests" / "contract" / "test_m10_03_diann_contract.py"
ARMS = {"A1": ("a1", "other"), "A2": ("a2", "same")}  # arm → (output stem, --search-db)
JOB_ID = "4242"

DICT = pa.dictionary(pa.int32(), pa.string())

# Arrow type of each documented column, as the 2.3.1 documentation allows it (F2).
TYPES: dict[str, pa.DataType] = {
    "Run": DICT,
    "Protein.Ids": DICT,
    "Protein.Group": DICT,
    "Stripped.Sequence": pa.string(),
    "Modified.Sequence": pa.string(),
    "Precursor.Charge": pa.int8(),
    "Global.Q.Value": pa.float32(),
    "Q.Value": pa.float32(),
    "PEP": pa.float32(),
    "Decoy": pa.int8(),
    "Precursor.Normalised": pa.float32(),
    "Precursor.Quantity": pa.float32(),
    "PG.MaxLFQ": pa.float32(),
    # unmapped extra 2.x columns named in the documentation; the converter ignores them
    "Precursor.Id": pa.string(),
    "RT": pa.float32(),
    "IM": pa.float32(),
}


def row(**cells: Any) -> dict[str, Any]:
    """One report row; keyword names use ``_`` for ``.`` (``Global_Q_Value``)."""
    return {k.replace("_", "."): v for k, v in cells.items()}


def table(
    rows: Sequence[Mapping[str, Any]], types: Mapping[str, pa.DataType] | None = None
) -> pa.Table:
    """A pyarrow table of ``rows`` (every row has the same keys), typed by ``TYPES``/``types``."""
    types = {**TYPES, **(types or {})}
    cols = list(rows[0])
    arrays = {}
    for c in cols:
        t = types.get(c, pa.string())
        vals = [r[c] for r in rows]
        if pa.types.is_dictionary(t):
            arrays[c] = pa.array(vals, pa.string()).dictionary_encode().cast(t)
        else:
            arrays[c] = pa.array(vals, t)
    return pa.table(arrays)


def write_parquet(path: Path, tbl: pa.Table) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(tbl, path)
    return path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ------------------------------------------------------------------ the 10-row job report


def job_rows(decoys: bool = True) -> list[dict[str, Any]]:
    """A 10-row report in the 2.3.1 shape: 3 runs, 2 decoys (when ``decoys``), Cont_ IDs."""
    peptides = ["PEPTMIDEK", "ACDLMTGK", "LLGNVLVR", "TTPSYVAFTDTER", "IGGIGTVPVGR"]
    out = []
    for i in range(10):
        decoy = int(decoys and i in (3, 8))
        seq = peptides[i % 5]
        prot = "Cont_P00761" if i == 6 else ("Cont_P00722;P00722" if i == 7 else f"P0A7V{i % 4}")
        out.append(
            row(
                Run=f"WT_{1 + i % 3}",
                Protein_Ids=prot,
                Protein_Group=prot,
                Stripped_Sequence=seq,
                Modified_Sequence=seq.replace("M", "M(UniMod:35)").replace("C", "C(UniMod:4)"),
                Precursor_Charge=2 + i % 2,
                Global_Q_Value=None if decoy and i == 8 else 0.001 * (i + 1) + (0.4 * decoy),
                PEP=0.002 * (i + 1),
                Decoy=decoy,
                Precursor_Normalised=1000.0 + 37 * i,
                PG_MaxLFQ=5000.0 + 11 * (i % 4),
            )
        )
    return out


# ------------------------------------------------------------------ the job sandbox

_PYTHON = r"""#!{python}
'''Recording GIAP_PYTHON wrapper (tests/helpers/diann_convert_m10.py).'''
import json, os, subprocess, sys
argv = sys.argv[1:]
state = os.environ["FAKE_SLURM_DIR"]
inputs = os.path.join(os.environ.get("OUT_DIR") or os.path.join(os.environ["GIAP_ROOT"], "pilot"),
                      "inputs")
code = argv[argv.index("-c") + 1] if "-c" in argv and argv.index("-c") + 1 < len(argv) else ""
rest = argv[argv.index("-c") + 2:] if code else []
if argv[:2] == ["-m", "pytest"]:
    kind = "pytest"
elif "giap.cli" in code and "convert" in rest:
    kind = "convert"
elif any(a.endswith("diann_convert.py") for a in argv):
    kind = "summary"
else:
    kind = "other"
rec = {"argv": argv, "kind": kind, "pythonpath": os.environ.get("PYTHONPATH", ""),
       "present": sorted(os.listdir(inputs)) if os.path.isdir(inputs) else None,
       "env": {k: os.environ.get(k) for k in ("GIAP_M10_03_REPORT", "GIAP_M10_03_ARM")}}
def log(rc):
    rec["rc"] = rc
    with open(os.path.join(state, "python.jsonl"), "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    sys.exit(rc)
def junit(tests, skipped, failures):
    path = None
    for i, a in enumerate(argv):
        if a.startswith(("--junitxml=", "--junit-xml=")):
            path = a.split("=", 1)[1]
        elif a in ("--junitxml", "--junit-xml") and i + 1 < len(argv):
            path = argv[i + 1]
    if not path:
        return
    cases = "".join(
        f'<testcase classname="test_m10_03_diann_contract" name="test_full_report[A1-f{i}]">'
        + ('<skipped message="fake skip"/>' if i < skipped else "")
        + ('<failure message="fake failure"/>' if skipped <= i < skipped + failures else "")
        + "</testcase>" for i in range(tests))
    with open(path, "w") as fh:
        fh.write('<?xml version="1.0" encoding="utf-8"?><testsuites><testsuite name="pytest" '
                 f'errors="0" failures="{failures}" skipped="{skipped}" tests="{tests}">'
                 f'{cases}</testsuite></testsuites>')
if kind == "pytest":
    mode = os.environ.get("FAKE_PYTEST", "ok")
    if mode == "fail":
        junit(5, 0, 1)
        print("FAILED tests/contract/test_m10_03_diann_contract.py::test_full_report[A1-f1]")
        print("1 failed, 4 passed in 0.01s"); log(1)
    if mode == "skip":
        junit(5, 5, 0)
        print("SKIPPED [5] tests/contract/test_m10_03_diann_contract.py:88: GIAP_ROOT is not set")
        print("5 skipped in 0.01s"); log(0)
    junit(5, 0, 0)
    print("5 passed in 0.01s"); log(0)
if kind == "convert" and os.environ.get("FAKE_CONVERT_RC"):
    c_path = os.path.join(state, "convert_calls")
    c = int(open(c_path).read()) + 1 if os.path.exists(c_path) else 1
    open(c_path, "w").write(str(c))
    if c == 1:
        print("giap convert: error: fake_failure: stub (0 rows)", file=sys.stderr)
        log(int(os.environ["FAKE_CONVERT_RC"]))
rc = subprocess.run([{python!r}, *argv]).returncode
nondet = os.environ.get("FAKE_CONVERT_NONDET")
if kind == "convert" and rc == 0 and nondet:
    flag = "--out" if nondet == "1" else nondet
    n_path = os.path.join(state, "convert_ok_calls")
    n = int(open(n_path).read()) + 1 if os.path.exists(n_path) else 1
    open(n_path, "w").write(str(n))
    if n == 2 and flag in rest:
        with open(rest[rest.index(flag) + 1], "ab") as fh:
            fh.write(b"\0")
log(rc)
"""

_MV = """#!/bin/bash
if [ -n "${FAKE_MV_FAIL:-}" ]; then
  for a in "$@"; do
    case "$a" in *.tmp.*) echo "mv: fake failure: $a" >&2; exit 1 ;; esac
  done
fi
exec {mv} "$@"
"""

_GIAP = """#!/bin/bash
PYTHONPATH="$GIAP_REPO/src" exec {python} -c 'from giap.cli import cli; cli()' "$@"
"""


@dataclass
class ConvertSandbox:
    sb: Sandbox

    @property
    def repo(self) -> Path:
        return self.sb.repo

    @property
    def pilot(self) -> Path:
        return self.sb.root / "pilot"

    @property
    def inputs(self) -> Path:
        return self.pilot / "inputs"

    def arm_dir(self, arm: str) -> Path:
        return self.pilot / arm

    def report(self, arm: str) -> Path:
        return self.arm_dir(arm) / f"{arm}_report.parquet"

    def stage(
        self,
        arm: str,
        tbl: pa.Table | None = None,
        *,
        ok: bool = True,
        sums: bool = True,
        stale_sum: bool = False,
        failed: bool = False,
    ) -> Path:
        """``$PILOT/<arm>/``: the report, ``check.json``, ``outputs.sha256`` (and ``FAILED``)."""
        d = self.arm_dir(arm)
        rep = write_parquet(self.report(arm), tbl if tbl is not None else table(job_rows()))
        (d / "check.json").write_text(json.dumps({"ok": ok, "failures": [] if ok else ["x"]}))
        digest = "0" * 64 if stale_sum else sha256(rep)
        if sums:
            (d / "outputs.sha256").write_text(f"{digest}  {rep.name}\n")
        if failed:
            (d / "FAILED").write_text("exit 1\n")
        return rep

    def run(self, *args: str, **env: str) -> subprocess.CompletedProcess[str]:
        assert SCRIPT.is_file(), f"{SCRIPT.relative_to(REPO)} does not exist (M10-03)"
        script = self.repo / "scripts" / "slurm" / "diann_convert.sbatch"
        if not script.is_file():  # created after the sandbox was copied: add and commit it
            shutil.copy2(SCRIPT, script)
            for a in (["add", "-A"], ["commit", "-q", "-m", "M10-03 sandbox: the job script"]):
                subprocess.run(
                    ["git", "-C", str(self.repo), *a],
                    check=True,
                    capture_output=True,
                    env=self.sb.env,
                )
        return self.sb.run(["bash", str(script), *args], env=env, cwd=self.sb.tmp, timeout=300)

    def calls(self, kind: str | None = None) -> list[dict[str, Any]]:
        recs = self.sb.jsonl("python.jsonl")
        return [r for r in recs if kind is None or r["kind"] == kind]

    def reset_calls(self) -> None:
        for name in ("python.jsonl", "convert_ok_calls", "convert_calls"):
            p = self.sb.state / name
            if p.exists():
                p.unlink()

    def inputs_files(self) -> list[str]:
        return sorted(p.name for p in self.inputs.iterdir()) if self.inputs.is_dir() else []

    def repo_is_clean(self) -> str:
        """``git status --porcelain --ignored`` of the sandbox checkout, bytecode caches excepted
        (``python -I`` ignores ``PYTHONDONTWRITEBYTECODE``)."""
        cmd = ["git", "-C", str(self.repo), "status", "--porcelain", "--ignored"]
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        return "\n".join(ln for ln in out.splitlines() if "__pycache__" not in ln).strip()


def convert_sandbox(tmp: Path) -> ConvertSandbox:
    sb = make_sandbox(tmp)
    ignore = shutil.ignore_patterns("__pycache__", "*.pyc")
    shutil.copytree(REPO / "src", sb.repo / "src", ignore=ignore)
    shutil.copytree(REPO / "validation", sb.repo / "validation", ignore=ignore)
    if CONTRACT.is_file():
        (sb.repo / "tests" / "contract").mkdir(parents=True)
        shutil.copy2(CONTRACT, sb.repo / "tests" / "contract" / CONTRACT.name)
    for args in (["add", "-A"], ["commit", "-q", "-m", "M10-03 sandbox: src, validation"]):
        subprocess.run(
            ["git", "-C", str(sb.repo), *args], check=True, capture_output=True, env=sb.env
        )
    py = sys.executable
    for name in ("python", "python3"):
        write_exe(sb.bin / name, _PYTHON.replace("{python!r}", repr(py)).replace("{python}", py))
    write_exe(sb.bin / "giap", _GIAP.replace("{python}", py))
    write_exe(sb.bin / "mv", _MV.replace("{mv}", shutil.which("mv") or "/bin/mv"))
    sb.env.update({"SLURM_JOB_ID": JOB_ID, "PYTHONDONTWRITEBYTECODE": "1"})
    sb.env.pop("GIAP_PYTHON", None)
    return ConvertSandbox(sb=sb)


# ------------------------------------------------------------------ running the converter


def giap_convert(
    report: Path, out: Path, abundance: Path, search_db: str
) -> subprocess.CompletedProcess[str]:
    """``giap convert --from diann`` of this checkout in a subprocess (stderr = the issue log)."""
    env = {**os.environ, "PYTHONPATH": str(REPO / "src")}
    cmd = [
        sys.executable,
        "-c",
        "from giap.cli import cli; cli()",
        "convert",
        "--from",
        "diann",
        str(report),
        "--out",
        str(out),
        "--abundance-out",
        str(abundance),
        "--search-db",
        search_db,
    ]
    return subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=300, check=False)


def load_summary_module():  # type: ignore[no-untyped-def]
    """Import ``validation/pilot/diann_convert.py`` (for ``summarize``)."""
    import importlib.util

    assert SUMMARY.is_file(), f"{SUMMARY.relative_to(REPO)} does not exist (M10-03)"
    spec = importlib.util.spec_from_file_location("giap_pilot_diann_convert", SUMMARY)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
