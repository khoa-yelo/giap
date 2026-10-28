"""Fixtures and sandboxes for the M10-10a tests (the A3 run's jobs and evidence).

Spec: docs/specs/M10-10-a3-run.md, piece 10a (entry gate, ``cascadia_convert.sbatch``,
``validation/pilot/cascadia_convert.py``, the ``A3_OUT_DIR`` override of ``denovo_pilot.sbatch``,
``run_evidence.py`` rules A10/A11 and rule set C for N labels, the 4-argument
``pilot_compare.sbatch``). Nothing here reaches the network, a GPU, Cascadia, real data or the
real ``$GIPPA_ROOT``; every ``.ssl`` row, hash and count is invented in ``tmp_path``.

* ``write_feasibility`` / ``write_assets``: an M10-08 ``a3_feasibility.json`` (aggregates only, the
  M10-08 field list) and an ``a3_assets.tsv`` whose ``cascadia_weights`` row carries a sha256.
* ``write_ssl``: a Cascadia 0.0.7-shaped ``.ssl`` (M10-09's header, ``score-type = UNKNOWN``).
* ``convert_sandbox``: ``helpers.slurm_sandbox.make_sandbox`` plus this checkout's ``src/`` and
  ``validation/`` (committed in the sandbox checkout), a passing gate, and a recording
  ``GIPPA_PYTHON`` (``$UV_PROJECT_ENVIRONMENT/bin/python``) that runs the **real** converter and
  summary, but fakes ``a3_feasibility.py evaluate`` (08b) by exit code: ``FAKE_L6="WT_2=3,WT_3=1"``
  (default 0 = pass). The fake writes a small aggregates JSON to the ``--out`` it is given.
  ``FAKE_CONVERT_NONDET=1`` appends one byte to the second conversion's ``--out`` (the determinism
  check of step 5 then has to fail).
* ``m10_04()``: ``helpers.pilot_runs_m10_04`` (M10-04, PR #296, not on main on 2026-10-24); the
  evidence tests build their run trees and the compare-job sandbox with it.

Conventions settled by the test author (the spec leaves them open; listed in the report):

* ``cascadia_convert.py`` has ``main(argv) -> int``; ``gate`` takes ``--feasibility`` and
  ``--assets`` (defaults: the committed files of the checkout it runs from).
* ``summarize(ssl_rows, peptides, issues, acquisition_runs) -> dict`` takes ``ssl_rows`` as
  ``{run: data rows}``, ``peptides`` as the converted table (a DataFrame with the file's columns),
  ``issues`` as ``gippa.io.issues.Issue`` objects (the parsed log) and the acquisition run labels in
  file order; it returns at least ``counts``, ``issues`` and ``failures``. A failure string starts
  with ``rule <n>`` (as M10-03's ``diann_convert.summarize``).
* ``a3.convert.json``: ``ssl`` is keyed by run, each ``{path, sha256, bytes, rows}``; ``outputs``
  has ``peptides`` and ``log``, each ``{path, sha256, bytes, rows}``.
* A11 needs the accepted ``out_a2``: ``arm_evidence``'s ``expect["a2"]`` (a path) and the CLI's
  ``--a2 OUT_A2``. C7 needs ``pass1/cmp_a1_a2``: ``compare_evidence(..., pair_cmp=PATH)`` and the
  CLI's ``--pair-cmp PATH``. Expected (not failed) findings go to ``expected_warnings`` (the C2
  entries name the table, e.g. ``compare_abundance:A1,A3:missing_in:A3``), the C7 finding to
  ``warnings``.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from helpers.slurm_sandbox import REPO, Sandbox, make_sandbox, write_exe

PILOT = REPO / "validation" / "pilot"
CASCADIA_CONVERT = PILOT / "cascadia_convert.py"
RUN_EVIDENCE = PILOT / "run_evidence.py"
ACQUISITION = PILOT / "inputs" / "acquisition.json"
SLURM = REPO / "scripts" / "slurm"
CONVERT_JOB = SLURM / "cascadia_convert.sbatch"
DENOVO_JOB = SLURM / "denovo_pilot.sbatch"
COMPARE_JOB = SLURM / "pilot_compare.sbatch"

RUNS = ("WT_1", "WT_2", "WT_3")
WEIGHTS_SHA = "7" * 64
FAIL_REASONS = (
    "too_few_peptides", "low_genome_match", "both", "route_unavailable",
    "structurally_incompatible", "tool_incompatible",
)  # fmt: skip
ASSET_COLUMNS = ("asset", "version", "url", "sha256", "bytes", "licence", "path", "retrieved")
SSL_HEADER = (
    "file", "scan", "charge", "sequence", "score-type", "score", "retention-time", "start-time",
    "end-time",
)  # fmt: skip
SEQ_LIKE = r"^[A-Z]{7,}$"  # criterion 2: no string value of the JSON looks like a sequence


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strings(doc: Any) -> list[str]:
    """Every string value (and key) of a JSON document."""
    if isinstance(doc, str):
        return [doc]
    if isinstance(doc, Mapping):
        return [s for k, v in doc.items() for s in (*strings(k), *strings(v))]
    if isinstance(doc, list | tuple):
        return [s for v in doc for s in strings(v)]
    return []


def load(path: Path, name: str, what: str) -> ModuleType:
    assert path.is_file(), f"{path.relative_to(REPO)} does not exist yet ({what})"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def load_cascadia_convert() -> ModuleType:
    return load(CASCADIA_CONVERT, "gippa_pilot_cascadia_convert", "M10-10a deliverable 3")


def m10_04() -> ModuleType:
    """``helpers.pilot_runs_m10_04``; a clear failure while M10-04 (PR #296) is not merged."""
    try:
        import helpers.pilot_runs_m10_04 as kit
    except ModuleNotFoundError as exc:  # pragma: no cover - depends on the merge state
        pytest.fail(f"M10-04 (PR #296: run_evidence.py and its test kit) is not merged: {exc}")
    return kit


# ------------------------------------------------------------------ the entry gate


def feasibility_doc(
    *, verdict: str = "pass", reason: str | None = None, weights: str = WEIGHTS_SHA, **over: Any
) -> dict[str, Any]:
    """An M10-08 ``a3_feasibility.json`` (every field of its list; aggregates only)."""
    counted = reason not in ("route_unavailable", "structurally_incompatible", "tool_incompatible")
    n = 2400 if counted else None
    doc: dict[str, Any] = {
        "schema_version": 1, "run": "WT_1", "cascadia_version": "0.0.7",
        "weights_sha256": weights, "route": "B2", "route_version": "3cb96aef",
        "mzml_sha256": "8" * 64, "mzml_check": {f"M{i}": True for i in range(1, 9)},
        "genome_sha256": "9" * 64, "cutoff": 0.8, "cutoff_source": "cli_default",
        "n_rows": n, "n_below_cutoff": 0 if counted else None, "n_len_7_30": n, "n_pass": n,
        "n_top": 1000 if counted else None, "n_match": 640 if counted else None,
        "n_nonstandard": 0 if counted else None,
        "share": 0.64 if counted else None,
        "share_ci95": [0.6098, 0.6692] if counted else None,
        "decoy_share": 0.009 if counted else None, "contaminant_share": 0.0 if counted else None,
        "matches_by_replicon": {"CP171242.1": 640} if counted else {},
        "length_hist": {"9": 1000} if counted else {},
        "rule": {"min_peptides": 500, "min_share": 0.2, "top_n": 1000},
        "fails": {"too_few_peptides": reason in ("too_few_peptides", "both"),
                  "low_genome_match": reason in ("low_genome_match", "both")}
        if counted else {"too_few_peptides": None, "low_genome_match": None},
        "verdict": verdict, "reason": reason, "evidence": None, "job_ids": ["4242"],
        "date": "2026-10-24",
    }  # fmt: skip
    doc.update(over)
    return doc


def write_feasibility(path: Path, doc: Mapping[str, Any] | str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = doc if isinstance(doc, str) else json.dumps(doc, indent=2, sort_keys=True) + "\n"
    path.write_text(text, encoding="utf-8")
    return path


def write_assets(path: Path, weights_sha: str = WEIGHTS_SHA) -> Path:
    """``a3_assets.tsv`` in M10-08's layout; the ``cascadia_weights`` row holds ``weights_sha``."""
    rows = [
        {"asset": "cascadia_weights", "version": "0.0.7 (Drive file cascadia.ckpt)",
         "url": "https://drive.usercontent.google.com/download?id=x", "sha256": weights_sha,
         "bytes": "1000", "licence": "Apache-2.0",
         "path": "$GIPPA_DB_DIR/cascadia/cascadia_0.0.7.ckpt", "retrieved": "2026-10-24"},
        {"asset": "unimod_table", "version": "psims 1.4.0", "url": "https://anaconda.org/x",
         "sha256": "6" * 64, "bytes": "500", "licence": "public domain",
         "path": "$GIPPA_DB_DIR/cascadia/unimod_tables.xml", "retrieved": "2026-10-24"},
    ]  # fmt: skip
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(ASSET_COLUMNS)] + ["\t".join(r[c] for c in ASSET_COLUMNS) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# ------------------------------------------------------------------ .ssl fixtures


def ssl_line(run: str, scan: int, charge: int, seq: str, score: float, st: str = "UNKNOWN") -> str:
    return "\t".join([f"/scratch/x/mzml/{run}.mzML", str(scan), str(charge), seq, st, repr(score),
                      "4.5", "4.4", "4.6"])  # fmt: skip


def write_ssl(path: Path, run: str, rows: Sequence[tuple[str, int, float]], **kw: str) -> Path:
    """``rows``: (sequence, charge, score); scans count from 0."""
    path.parent.mkdir(parents=True, exist_ok=True)
    body = [ssl_line(run, i, c, s, sc, **kw) for i, (s, c, sc) in enumerate(rows)]
    path.write_text("\n".join(["\t".join(SSL_HEADER), *body]) + "\n", encoding="utf-8")
    return path


#: The 3-run, 12-row fixture of criterion 2 (counts worked out by hand in the test).
TWELVE: dict[str, list[tuple[str, int, float]]] = {
    "WT_1": [("PEPTLDEK", 2, 0.95), ("M[+15.995]PEPTLDEKR", 3, 0.88), ("SAMPLLER", 2, 0.91),
             ("PEPTLDEK", 3, 0.97)],
    "WT_2": [("ACDEFGHK", 2, 0.85), ("LLLKR", 2, 0.99), ("PEPTLDEK", 2, 0.93),
             ("C[+57.021]AMPLLERK", 3, 0.90)],
    "WT_3": [("SAMPLLER", 4, 0.86), ("ACDEFGHK", 3, 0.92), ("V" * 31, 2, 0.89),
             ("GGGGGGG", 2, 0.94)],
}  # fmt: skip

#: Small per-run fixtures for the job sandbox (3 / 2 / 4 data rows).
JOB_SSL: dict[str, list[tuple[str, int, float]]] = {
    "WT_1": [("PEPTLDEK", 2, 0.95), ("M[+15.995]PEPTLDEKR", 3, 0.88), ("SAMPLLER", 2, 0.91)],
    "WT_2": [("ACDEFGHK", 2, 0.85), ("PEPTLDEK", 2, 0.93)],
    "WT_3": [("SAMPLLER", 4, 0.86), ("ACDEFGHK", 3, 0.92), ("GGGGGGGK", 2, 0.94),
             ("LLLLLLLK", 2, 0.81)],
}  # fmt: skip


# ------------------------------------------------------------------ the conversion job sandbox

_PYTHON = r"""#!{python}
'''Recording GIPPA_PYTHON for the M10-10 conversion job (tests/helpers/a3_run_m10_10.py).'''
import json, os, subprocess, sys
argv = sys.argv[1:]
state = os.environ["FAKE_SLURM_DIR"]
code = argv[argv.index("-c") + 1] if "-c" in argv and argv.index("-c") + 1 < len(argv) else ""
rest = argv[argv.index("-c") + 2:] if code else []
script = next((a for a in argv if a.endswith(".py")), "")
if script.endswith("a3_feasibility.py") and "evaluate" in argv:
    kind = "evaluate"
elif "gippa.cli" in code and "convert" in rest:
    kind = "convert"
elif script.endswith("cascadia_convert.py"):
    kind = "gate" if "gate" in argv else "summary" if "summary" in argv else "cascadia_convert"
else:
    kind = "other"
rec = {"argv": argv, "kind": kind}
def log(rc):
    rec["rc"] = rc
    with open(os.path.join(state, "python.jsonl"), "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    sys.exit(rc)
if kind == "evaluate":
    ssl = next((a for a in argv if a.endswith(".ssl")), "")
    run = os.path.basename(ssl)[: -len(".ssl")] if ssl else ""
    plan = dict(p.split("=", 1) for p in os.environ.get("FAKE_L6", "").split(",") if "=" in p)
    rc = int(plan.get(run, "0"))
    rec["run"] = run
    if rc in (0, 3) and "--out" in argv:
        out = argv[argv.index("--out") + 1]
        os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
        doc = {"schema_version": 1, "run": run, "verdict": "pass" if rc == 0 else "fail",
               "reason": None if rc == 0 else "too_few_peptides", "n_pass": 0 if rc else 600}
        with open(out, "w") as fh:
            fh.write(json.dumps(doc, indent=2, sort_keys=True) + "\n")
    log(rc)
rc = subprocess.run([{python!r}, *argv]).returncode
if kind == "convert" and rc == 0 and os.environ.get("FAKE_CONVERT_NONDET") and "--out" in rest:
    n_path = os.path.join(state, "convert_ok_calls")
    n = int(open(n_path).read()) + 1 if os.path.exists(n_path) else 1
    open(n_path, "w").write(str(n))
    if n == 2:
        with open(rest[rest.index("--out") + 1], "ab") as fh:
            fh.write(b"\0")
log(rc)
"""


@dataclass
class ConvertSandbox:
    sb: Sandbox

    @property
    def repo(self) -> Path:
        return self.sb.repo

    @property
    def a3(self) -> Path:
        return self.sb.root / "pilot" / "A3"

    @property
    def inputs(self) -> Path:
        return self.sb.root / "pilot" / "inputs"

    def ssl(self, run: str, a3: Path | None = None) -> Path:
        return (a3 or self.a3) / run / f"{run}.ssl"

    def stage(
        self, rows: Mapping[str, Sequence[tuple[str, int, float]]] | None = None,
        *, a3: Path | None = None, score_type: Mapping[str, str] | None = None,
    ) -> list[Path]:  # fmt: skip
        """The three ``.ssl`` under ``$A3/<RUN>/`` (``JOB_SSL`` by default)."""
        rows = JOB_SSL if rows is None else rows
        st = score_type or {}
        return [write_ssl(self.ssl(r, a3), r, rows[r], st=st.get(r, "UNKNOWN")) for r in rows]

    def set_gate(self, doc: Mapping[str, Any] | str | None) -> None:
        """Rewrite (``None``: delete) the sandbox checkout's ``a3_feasibility.json``, committed."""
        path = self.repo / "validation" / "pilot" / "metrics" / "a3_feasibility.json"
        if doc is None:
            path.unlink()
        else:
            write_feasibility(path, doc)
        self._commit("M10-10 sandbox: gate")

    def _commit(self, msg: str) -> None:
        for a in (["add", "-A"], ["commit", "-q", "--allow-empty", "-m", msg]):
            subprocess.run(["git", "-C", str(self.repo), *a], check=True, capture_output=True,
                           env=self.sb.env)  # fmt: skip

    def run(self, *args: str, **env: str) -> subprocess.CompletedProcess[str]:
        assert CONVERT_JOB.is_file(), (
            f"{CONVERT_JOB.relative_to(REPO)} does not exist yet (M10-10a deliverable 2)"
        )
        script = self.repo / "scripts" / "slurm" / CONVERT_JOB.name
        if not script.is_file():
            shutil.copy2(CONVERT_JOB, script)
            self._commit("M10-10 sandbox: the job script")
        from helpers.slurm_sandbox import assert_honours_gippa_repo

        assert_honours_gippa_repo(script.read_text(), script.name)
        return self.sb.run(["bash", str(script), *args], env=env, cwd=self.sb.tmp, timeout=300)

    def calls(self, kind: str | None = None) -> list[dict[str, Any]]:
        return [r for r in self.sb.jsonl("python.jsonl") if kind is None or r["kind"] == kind]

    def inputs_files(self) -> list[str]:
        return sorted(p.name for p in self.inputs.iterdir()) if self.inputs.is_dir() else []

    def repo_is_clean(self) -> str:
        cmd = ["git", "-C", str(self.repo), "status", "--porcelain", "--ignored"]
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        return "\n".join(ln for ln in out.splitlines() if "__pycache__" not in ln).strip()


def convert_sandbox(tmp: Path) -> ConvertSandbox:
    sb = make_sandbox(tmp)
    ignore = shutil.ignore_patterns("__pycache__", "*.pyc")
    shutil.copytree(REPO / "src", sb.repo / "src", ignore=ignore)
    shutil.copytree(REPO / "validation", sb.repo / "validation", ignore=ignore,
                    dirs_exist_ok=True)  # fmt: skip
    pilot = sb.repo / "validation" / "pilot"
    write_feasibility(pilot / "metrics" / "a3_feasibility.json", feasibility_doc())
    write_assets(pilot / "a3_assets.tsv")
    box = ConvertSandbox(sb=sb)
    box._commit("M10-10 sandbox: src, validation, a passing gate")
    py = sys.executable
    for name in ("python", "python3"):
        write_exe(sb.bin / name, _PYTHON.replace("{python!r}", repr(py)).replace("{python}", py))
    sb.env.update({"SLURM_JOB_ID": "5151", "PYTHONDONTWRITEBYTECODE": "1"})
    sb.env.pop("GIPPA_PYTHON", None)
    return box


# ------------------------------------------------------------------ evidence fixtures (M10-04 kit)

SHA_A3_PEPTIDES = "4" * 64  # sha256 of a3.tsv.gz (invented)
REFERENCE_ROLE = "reference_peptides"
SPECTRUM_PAIRING_WARNING = "preflight:peptides.spectrum_pairing"  # C17 in the manifest


def write_mapping_fdr(out: Path, *, n_bins: int = 10, flags: str | None = None,
                      placement: bool = True) -> Path:  # fmt: skip
    """``OUT/mapping/mapping_fdr.tsv``: placement rows and ``n_bins`` ``denovo_precision`` rows."""
    from helpers.viz_out_tree import write_schema_table

    rows: list[dict[str, object]] = []
    if placement:
        rows += [{"section": "placement", "length_class": lc, "n_diff": 0, "targets": 40,
                  "decoys": 1, "pi_wrong": 1.0, "fdr_stratum": 0.025, "n_accepted": 40,
                  "decoy_type": "shuffle"} for lc in ("short", "long")]  # fmt: skip
    rows += [{"section": "denovo_precision", "score_bin": b, "precision": 0.5 + 0.04 * b,
              "ci_low": 0.4 + 0.04 * b, "ci_high": 0.6 + 0.035 * b, "c_value": 0.95,
              "flags": flags} for b in range(n_bins)]  # fmt: skip
    return write_schema_table(out / "mapping" / "mapping_fdr.tsv", "mapping_fdr", rows)


def add_reference_input(out: Path, sha: str) -> None:
    """The manifest row ``launch_sources`` adds for ``--reference-peptides`` (self-hashed)."""
    path = out / "run_manifest.json"
    doc = json.loads(path.read_text(encoding="utf-8"))
    row = {"role": REFERENCE_ROLE, "file": "a2.tsv.gz", "path": "/pilot/inputs/a2.tsv.gz",
           "sha256": sha, "size": 2048}  # fmt: skip
    doc["inputs"].append(row)
    path.write_text(json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def a3_genes() -> list[dict[str, object]]:
    """A3's genes: detected as A1/A2's, no abundance (de novo rows carry no intensity)."""
    rows = m10_04().arm_genes(shift=0.25)
    for r in rows:
        r.update(abundance_log2=None, abundance_raw=None, abundance_type="none")
    return rows


def write_a3_arm(
    out: Path,
    *,
    reference: str | None,
    fdr_bins: int = 10,
    fdr_flags: str | None = None,
    fdr: bool = True,
    spectrum_warning: bool = True,
    **manifest_kw: Any,
) -> Path:
    """An A3 ``OUT/`` tree (example 7 by default: reference = A2's peptides table)."""
    kit = m10_04()
    inputs = {"peptides": SHA_A3_PEPTIDES, **manifest_kw.pop("inputs", {})}
    kit.write_arm(out, genes=a3_genes(), inputs=inputs, **manifest_kw)
    if (out / "run_manifest.json").is_file():
        if reference is not None:
            add_reference_input(out, reference)
        if spectrum_warning:
            doc = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            kit.edit_manifest(out, warnings=sorted({*doc.get("warnings", []),
                                                    SPECTRUM_PAIRING_WARNING}))  # fmt: skip
    if fdr:
        write_mapping_fdr(out, n_bins=fdr_bins, flags=fdr_flags)
    return out


def a3_expect(out_a2: Path, **over: Any) -> dict[str, Any]:
    kit = m10_04()
    inputs = {**kit.EXPECT_INPUTS, "peptides": SHA_A3_PEPTIDES, REFERENCE_ROLE: kit.SHA_PEPTIDES}
    return kit.expect(arm="A3", inputs=over.pop("inputs", inputs), a2=out_a2, **over)


def a3_cli_args(out: Path, dest: Path, out_a2: Path, **kw: str) -> list[str]:
    kit = m10_04()
    args = kit.arm_cli_args(out, dest, arm="A3", peptides=kw.pop("peptides", SHA_A3_PEPTIDES), **kw)
    return [*args, "--expect-input", f"{REFERENCE_ROLE}={kit.SHA_PEPTIDES}", "--a2", str(out_a2)]


def gippa_compare_n(
    runs: Sequence[Path], labels: Sequence[str], cmp: Path, *, env_root: Path, bootstrap: int = 200
) -> subprocess.CompletedProcess[str]:
    """The real ``gippa compare RUN… --outdir CMP --labels L1,L2,…`` of this checkout."""
    env_root.mkdir(parents=True, exist_ok=True)
    env = {"PATH": "/usr/bin:/bin", "HOME": str(env_root), "LANG": "C.UTF-8",
           "PYTHONPATH": str(REPO / "src"), "GIPPA_DB_DIR": str(env_root),
           "GIPPA_ROOT": str(env_root), "PYTHONDONTWRITEBYTECODE": "1",
           "MPLCONFIGDIR": str(env_root / "mpl")}  # fmt: skip
    cmd = [sys.executable, "-c", "from gippa.cli import cli; cli()", "compare",
           *(str(r) for r in runs), "--outdir", str(cmp), "--labels", ",".join(labels),
           "--set", f"stats.bootstrap={bootstrap}"]  # fmt: skip
    return subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=600, check=False)


def three_arms(tmp: Path) -> tuple[Path, Path, Path]:
    """``out_a1``, ``out_a2`` (pass 1) and ``out_a3`` (no abundance) on one genome."""
    kit = m10_04()
    a1, a2 = kit.two_arms(tmp)
    a3 = write_a3_arm(tmp / "out_a3", reference=kit.SHA_PEPTIDES)
    return a1, a2, a3
