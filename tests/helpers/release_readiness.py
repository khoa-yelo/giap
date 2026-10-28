"""Fixture repository for ``scripts/release/readiness.py`` (issue M11-03, acceptance criterion 1).

``Fixture.build(tmp_path)`` writes a git repository in ``tmp_path/repo`` with two commits (X, then
R = HEAD) and every evidence file of a release that passes: fake ``PROGRESS.md`` (gate blocks
M1-M10, M10 with a Khoa review line), ``validation/pilot/metrics/metrics.json``,
``validation/pilot/REPORT.md``, ``benchmarks/baseline.json``, version 1.0.0 in
``src/gippa/__init__.py`` and ``CITATION.cff``, a built ``CHANGELOG.md``, an empty ``changelog.d/``,
and copies of the checkout's ``workflow/conf/containers.config``, ``src/gippa/annotate/tools.py``,
``NOTICE.md`` and ``LICENSE`` (H4/H5 read them). ``tmp_path/logs`` holds ``nightly.tsv`` and
``weekly.tsv`` rows for R (date, job_task, commit (8-hex short), exit, summary: the
``scripts/slurm/nightly.sbatch`` layout), and ``tmp_path/check_runs.json`` a ``gh api
.../check-runs`` response with every run ``success``.

A test changes one input (``write_json``/``write_text``/``commit``), runs ``compute`` with
``--offline`` (never the network) and reads ``readiness.json``.
"""

from __future__ import annotations

import functools
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType
from typing import Any

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "release" / "readiness.py"
VERSION = "1.0.0"

ROW_IDS = [f"R{i}" for i in range(1, 11)] + [f"H{i}" for i in range(1, 7)]
SOURCES = {
    "R1": "repo", "R2": "farmshare", "R3": "farmshare", "R4": "farmshare", "R5": "both",
    "R6": "both", "R7": "repo", "R8": "repo", "R9": "repo", "R10": "repo",
    "H1": "repo", "H2": "repo", "H3": "repo", "H4": "repo", "H5": "repo", "H6": "farmshare",
}  # fmt: skip
STATUSES = {"pass", "fail", "pending", "deferred", "waived"}

_GIT_ENV = {
    "GIT_AUTHOR_NAME": "Fixture", "GIT_AUTHOR_EMAIL": "fixture@example.invalid",
    "GIT_COMMITTER_NAME": "Fixture", "GIT_COMMITTER_EMAIL": "fixture@example.invalid",
    "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_SYSTEM": os.devnull,
}  # fmt: skip

COPIED = (
    "workflow/conf/containers.config",
    "src/gippa/annotate/tools.py",
    "NOTICE.md",
    "LICENSE",
)


@functools.cache
def load_module() -> ModuleType:
    """``scripts/release/readiness.py`` as a module (its CLI sits behind ``__main__``)."""
    spec = importlib.util.spec_from_file_location("gippa_release_readiness", SCRIPT)
    assert spec is not None and spec.loader is not None, f"{SCRIPT} missing"
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def git(repo: Path, *args: str) -> str:
    env = {**os.environ, **_GIT_ENV}
    out = subprocess.run(
        ["git", "-C", str(repo), *args], env=env, check=True, capture_output=True, text=True
    )
    return out.stdout.strip()


def progress_md(gates: range | list[int] = range(1, 11), khoa_m10: bool = True) -> str:
    out = ["# GIPPA progress (fixture)", "", "## Gates", ""]
    for k in gates:
        out += [f"### Gate M{k}", f"- Evidence: fixture evidence for M{k}."]
        if k == 10 and khoa_m10:
            out.append("- Khoa review: approved 2026-11-01")
        out.append("")
    return "\n".join(out)


def report_md(statuses: tuple[str, ...] = ("pass",) * 5) -> str:
    rows = "\n".join(f"| bullet {i} | {s} | fixture |" for i, s in enumerate(statuses))
    return (
        "# Pilot report (fixture)\n\n## 1. Status\n\n| §15.6 bullet | Status | Evidence |\n"
        f"|---|---|---|\n{rows}\n\nNotes.\n"
    )


def metrics_doc(share: float = 0.9971, passed: bool | None = None) -> dict[str, Any]:
    p = {"status": "pass"}
    return {
        "schema": "gippa-pilot-metrics/1",
        "te1": {
            "n_peptides": 10000,
            "exact_inframe": round(share * 10000),
            "share_exact_inframe": share,
            "threshold": 0.995,
            "passed": share >= 0.995 if passed is None else passed,
        },
        "mapping_fdr": {"value": 0.0083, "threshold": 0.01},
        "te_real": {"schema": "gippa-pilot-te-real/1",
                    "te2": dict(p), "te3": dict(p), "te4": dict(p), "te5": dict(p)},
        "te7": {"schema": "gippa-pilot-te7/1", "gippa": dict(p),
                "diann": {"status": "identical_report"}, "end_to_end": True},
    }  # fmt: skip


def baseline_doc(standard_s: float | None = 2280.0, core_s: float | None = 1500.0) -> dict:
    cells: dict[str, Any] = {}
    for cell, s in (("k12_standard", standard_s), ("k12_core", core_s)):
        if s is not None:
            cells[cell] = {"median_s": s, "threads": 16, "mem_mb": 64000,
                           "nodes": ["wheat-05"] * 3, "cpu": "Intel Xeon Gold 6148"}  # fmt: skip
    doc: dict[str, Any] = {"annotation_steps": {"wheat": {"median_s": 12.5}}}
    if cells:
        doc["standard_run"] = {"wheat": cells}
    return doc


def adr_text(number: int, status: str, row: str = "R7") -> str:
    return (
        f"# ADR-{number:04d}: release runtime waiver (fixture)\n- Status: {status}\n"
        f"- Date: 2026-11-01\n\n## Decision\nWaives readiness {row}: fixture.\n"
    )


def nightly_row(date: str, job: str, commit: str, exit_: int, summary: str) -> str:
    return f"{date}\t{job}\t{commit}\t{exit_}\t{summary}\n"


def passing_nightly(commit: str, date: str = "2026-11-02", job: str = "2000") -> str:
    return "".join(
        nightly_row(date, f"{job}_{t}", commit, 0, s)
        for t, s in enumerate(
            (
                "integration: 812 passed, 92 skipped in 1438.86s (0:23:58)",
                "benchmarks: 57 passed in 519.92s (0:08:39) (wheat wheat-02)",
                "calibration_nightly: 17290 passed, 108 skipped, 13 xfailed in 20556.28s",
                "calibration_full: 0 failed, 3120 passed in 30111.00s",
                "calibration_attacks: 240 passed in 5100.00s",
            )
        )
    )


def passing_weekly(commit: str, date: str = "2026-11-01", job: str = "3000") -> str:
    return nightly_row(
        date, f"{job}_0", commit, 0,
        "k12_mg1655: gippa run exit 0; TE8 /scratch/ci/validation/2026-11-01_te8_3000.json;",
    )  # fmt: skip


@dataclass
class Fixture:
    root: Path
    repo: Path
    logs: Path
    out: Path
    check_runs: Path
    x: str = ""  # first commit (full hash)
    r: str = ""  # HEAD after build (full hash)
    extra: dict[str, str] = field(default_factory=dict)

    @classmethod
    def build(cls, tmp: Path, *, nightly: bool = True, weekly: bool = True) -> Fixture:
        repo, logs = tmp / "repo", tmp / "logs"
        fx = cls(tmp, repo, logs, repo / "validation" / "release" / VERSION, tmp / "cr.json")
        repo.mkdir()
        logs.mkdir()
        git(repo, "init", "-q", "-b", "main")
        for rel in COPIED:
            dst = repo / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(REPO / rel, dst)
        fx.write_text("src/gippa/__init__.py", f'"""GIPPA."""\n\n__version__ = "{VERSION}"\n')
        fx.write_text("src/gippa/cli.py", '"""CLI (fixture)."""\n')
        fx.write_text("src/gippa/stats/bh.py", '"""BH (fixture)."""\n')
        fx.write_text(
            "CITATION.cff",
            f'cff-version: 1.2.0\ntitle: "GIPPA"\nversion: {VERSION}\nlicense: MIT\n',
        )
        fx.write_text(
            "CHANGELOG.md",
            "# Changelog\n\n<!-- towncrier release notes start -->\n\n"
            f"## GIPPA {VERSION} (2026-11-01)\n\n### Features\n\n- First release.\n",
        )
        fx.write_text("changelog.d/.gitkeep", "")
        fx.write_text("PROGRESS.md", progress_md())
        fx.write_text("validation/pilot/REPORT.md", report_md())
        fx.write_json("validation/pilot/metrics/metrics.json", metrics_doc())
        fx.write_json("benchmarks/baseline.json", baseline_doc())
        fx.write_text("docs/user-guide/install.md", "# Install\n")
        fx.write_text("docs/adr/ADR-0001-fixture.md", adr_text(1, "accepted-by-khoa", "none"))
        fx.x = fx.commit("X: fixture")
        fx.write_text("README.md", "# GIPPA\n")
        fx.r = fx.commit("R: fixture")
        if nightly:
            fx.set_nightly(passing_nightly(fx.short(fx.r)))
        if weekly:
            fx.set_weekly(passing_weekly(fx.short(fx.r)))
        fx.set_check_runs(["lint", "tests (3.11)", "package"])
        return fx

    # ---------------------------------------------------------------- inputs
    def write_text(self, rel: str, text: str) -> Path:
        p = self.repo / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
        return p

    def write_json(self, rel: str, doc: Any) -> Path:
        return self.write_text(rel, json.dumps(doc, indent=2) + "\n")

    def commit(self, msg: str) -> str:
        git(self.repo, "add", "-A")
        git(self.repo, "commit", "-q", "-m", msg)
        return git(self.repo, "rev-parse", "HEAD")

    @staticmethod
    def short(sha: str) -> str:
        return sha[:8]

    def set_nightly(self, text: str) -> None:
        (self.logs / "nightly.tsv").write_text(text, encoding="utf-8")

    def set_weekly(self, text: str) -> None:
        (self.logs / "weekly.tsv").write_text(text, encoding="utf-8")

    def set_check_runs(
        self, names: list[str], conclusion: dict[str, str] | None = None, sha: str | None = None
    ) -> None:
        conclusion = conclusion or {}
        head = sha or git(self.repo, "rev-parse", "HEAD").strip()
        runs = [
            {
                "name": n,
                "status": "completed",
                "conclusion": conclusion.get(n, "success"),
                "head_sha": head,
            }
            for n in names
        ]
        self.check_runs.write_text(
            json.dumps({"total_count": len(runs), "check_runs": runs}), encoding="utf-8"
        )

    # ---------------------------------------------------------------- runs
    def compute(
        self, *extra: str, check_runs: bool = True, logs: bool = True
    ) -> tuple[subprocess.CompletedProcess[str], dict[str, Any] | None]:
        args = ["compute", "--version", VERSION, "--repo", str(self.repo)]
        if logs:
            args += ["--logs", str(self.logs)]
        args += ["--out", str(self.out), "--offline"]
        if check_runs:
            args += ["--check-runs", str(self.check_runs)]
        proc = run_cli(*args, *extra)
        path = self.out / "readiness.json"
        doc = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None
        return proc, doc

    def verify(self, readiness: Path | None = None) -> subprocess.CompletedProcess[str]:
        path = readiness or self.out / "readiness.json"
        return run_cli("verify", str(path), "--repo", str(self.repo))


def run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    env = {**os.environ, **_GIT_ENV}
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        env=env, capture_output=True, text=True, timeout=120, check=False,
    )  # fmt: skip


def row(doc: dict[str, Any] | None, rid: str) -> dict[str, Any]:
    assert doc is not None, "compute wrote no readiness.json"
    hits = [r for r in doc["rows"] if r["id"] == rid]
    assert len(hits) == 1, f"row {rid}: {len(hits)} entries"
    return hits[0]


def text_of(r: dict[str, Any]) -> str:
    """Every string of a row (status text, value, threshold, evidence) for substring checks."""
    return json.dumps(r, ensure_ascii=False)
