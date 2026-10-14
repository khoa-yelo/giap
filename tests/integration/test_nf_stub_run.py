"""Stub run of the skeleton pipeline and ``-resume`` (M2-06 example 8; GIAP_SPEC.md §3.8 "Resume").

Runs ``workflow/main.nf`` with ``-profile test,stub -stub-run`` (the profiles ``ci.yml`` gives
nf-test) twice in the same work directory; the second run with ``-resume`` must take every task
from the cache. Local executor only, inside the session allocation.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.needs_nextflow

REPO = Path(__file__).resolve().parents[2]
MAIN = REPO / "workflow" / "main.nf"


def _run(tmp: Path, trace: str, resume: bool) -> list[dict[str, str]]:
    out = tmp / "OUT"
    (tmp / "params.json").write_text(json.dumps({"outdir": str(out)}))
    cmd = ["nextflow", "-q", "run", str(MAIN), "-profile", "test,stub", "-stub-run"]
    cmd += ["-params-file", str(tmp / "params.json"), "-work-dir", str(out / "work")]
    cmd += ["-with-trace", str(out / "logs" / trace)]
    if resume:
        cmd.append("-resume")
    env = {**os.environ, "NXF_ANSI_LOG": "false"}
    res = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True, env=env, timeout=900)
    assert res.returncode == 0, f"{' '.join(cmd)}\n{res.stdout}\n{res.stderr}"
    with (out / "logs" / trace).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_stub_run_then_resume_is_fully_cached(tmp_path: Path) -> None:
    first = _run(tmp_path, "trace1.tsv", resume=False)
    assert first, "stub run executed no task"
    assert any("GIAP_STEP" in r["name"] for r in first), [r["name"] for r in first]
    assert all(r["status"] == "COMPLETED" for r in first), first
    versions = list((tmp_path / "OUT" / "work").rglob("versions.yml"))
    assert versions, "no task emitted versions.yml"

    second = _run(tmp_path, "trace2.tsv", resume=True)
    assert sorted(r["name"] for r in second) == sorted(r["name"] for r in first)
    assert all(r["status"] == "CACHED" for r in second), second


def test_stub_run_publishes_into_outdir(tmp_path: Path) -> None:
    _run(tmp_path, "trace.tsv", resume=False)
    out = tmp_path / "OUT"
    published = [
        p
        for p in out.rglob("*")
        if p.is_file() and not p.is_relative_to(out / "work") and not p.is_relative_to(out / "logs")
    ]
    assert published, "publishDir copied nothing into OUT/"
    assert not any(p.is_symlink() for p in published), "publishDir mode must be copy"
