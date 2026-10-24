"""The M7-15 real-tool arm, checked on its recorded run (docs/specs/M7-15-integration-hardening.md
item 4, example 5; lead decision 2, 2026-10-22: K-12 with the in-silico TE8 peptides, Slurm only).

The run itself is one ``sbatch scripts/slurm/giap_run.sbatch`` (``-profile farmshare,apptainer``;
Promotech and OSTIR from their conda envs, ADR-0008) on ``tests/data/genomes/k12_mg1655.fna.gz``
with the TE8 in-silico tryptic peptide set and ``--genome-relation related``; its job ID goes to
PROGRESS.md. This test never starts it: it reads the finished output directory named by
``$GIAP_M7_15_K12_OUT`` (skipped when unset) and checks example 5:

* the run completed (``run_manifest.json``, the profile and ``related`` on its command line);
* the second round ran on real candidates: PROMOTECH_NOVEL and OSTIR_NOVEL tasks ``COMPLETED`` in
  ``logs/trace.tsv`` (none FAILED, none skipped as precomputed);
* ``tables/novel_orfs.tsv`` has rows and validates; the seven Module 5/6 files exist;
* their runtimes are in ``benchmarks/baseline.json`` (tests/unit/test_m7_15_baseline.py pins the
  entry; here: the recorded Slurm job is this run's).
"""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m7_chain as m7

pytestmark = [pytest.mark.e2e, pytest.mark.slow]

REPO = Path(__file__).resolve().parents[2]
ENV = "GIAP_M7_15_K12_OUT"
BASELINE_KEY = "novel_second_round_k12"


@pytest.fixture(scope="module")
def out() -> Path:
    path = os.environ.get(ENV)
    if not path:
        pytest.skip(f"${ENV} not set: no recorded K-12 real-tool run")
    p = Path(path)
    assert p.is_dir(), p
    return p


def _trace(out: Path) -> list[dict[str, str]]:
    with (out / "logs" / "trace.tsv").open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_example_5_run_completed_with_the_real_profile(out: Path) -> None:
    man = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    argv = " ".join(man["command_line"])
    assert "farmshare,apptainer" in argv, argv
    assert man["config"]["genome_relation"] == "related"
    assert "k12" in argv.lower(), argv


def test_example_5_second_round_tools_ran(out: Path) -> None:
    rows = _trace(out)
    for proc in ("PROMOTECH_NOVEL", "OSTIR_NOVEL"):
        mine = [r for r in rows if r["name"].split(" (")[0].split(":")[-1] == proc]
        assert mine, f"{proc} never ran (precomputed or no candidate start)"
        assert all(r["status"] in ("COMPLETED", "CACHED") for r in mine), mine
        assert any(r["status"] == "COMPLETED" for r in mine), mine
    assert not [r for r in rows if r["status"] == "FAILED"]


def test_example_5_novel_orfs_non_empty_and_valid(out: Path) -> None:
    for rel in m7.SEVEN:
        assert (out / rel).is_file(), rel
    lines = (out / "tables" / "novel_orfs.tsv").read_text(encoding="utf-8").splitlines()
    assert len(lines) > 1, "tables/novel_orfs.tsv is header-only on K-12"
    m7.validate_tsv(out / "tables" / "novel_orfs.tsv", "novel_orfs")


def test_example_5_runtimes_recorded_for_this_run(out: Path) -> None:
    base = json.loads((REPO / "benchmarks" / "baseline.json").read_text(encoding="utf-8"))
    entry = base.get(BASELINE_KEY)
    assert entry, f"benchmarks/baseline.json has no {BASELINE_KEY}"
    man = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    job = os.environ.get("GIAP_M7_15_K12_JOB") or str(man.get("slurm_job_id", ""))
    if job:
        assert str(entry["slurm_job"]) == job, (entry["slurm_job"], job)
