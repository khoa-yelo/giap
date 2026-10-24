"""M2-14b example 5 and item 4: ``-resume`` reruns only the Promoter Calculator batch whose region
changed (docs/specs/M2-14b-promoter-calculator-batching.md; GIAP_SPEC.md §6.3 ``-resume``).

Stub runs of ``workflow/main.nf`` (ANNOTATE, ``modules = ['core']``) with ``-profile test,stub
-stub-run`` and ``params.stub.tu_regions = 250``: the stubbed ``operon_pairs`` step emits
``tu_upstream/tu000.fa`` ... ``tu249.fa`` (annotate.nf ``scatterSpec``), which batch into
``pc_tu000_tu099``, ``pc_tu100_tu199`` and ``pc_tu200_tu249`` (tag form pinned by the test author,
spec item 4). All runs of one test share one work directory; each run's ``-with-trace`` gives the
task status (``CACHED`` / ``COMPLETED``) per PROMOTER_CALCULATOR task.

"One region changes" (test author): after the first run, the bytes of ``tu150.fa`` in the cached
``operon_pairs`` task directory are rewritten (a one-record FASTA instead of the empty stub file).
On ``-resume`` operon_pairs stays cached (its inputs are unchanged) and emits the changed file;
Nextflow's standard input hash (path, size, mtime) then changes for exactly the batch that holds
tu150, as a rerun upstream step would deliver one new region file. The pattern follows
tests/integration/test_report_resume_stub.py (M9-26). Local executor inside the allocation.
"""

from __future__ import annotations

import csv
import json
import os
import re
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.needs_nextflow

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"
DATA = WF / "tests" / "data" / "annotate"
N_REGIONS = 250
BATCHES = ["pc_tu000_tu099", "pc_tu100_tu199", "pc_tu200_tu249"]


def _params(tmp: Path) -> Path:
    data = {
        "outdir": str(tmp / "OUT"),
        "modules": ["core"],
        "inputs": {"genome": str(DATA / "genome.fa"), "replicons": str(DATA / "replicons.tsv")},
        "annotation": json.loads((DATA / "annotation.json").read_text()),
        "stub": {"windows_chunks": 1, "tu_regions": N_REGIONS},
    }
    path = tmp / "params.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def _run(tmp: Path, n: int, resume: bool) -> dict[str, str]:
    """One stub run; returns ``{PROMOTER_CALCULATOR tag: status}``."""
    trace = tmp / f"trace{n}.tsv"
    cmd = ["nextflow", "-q", "run", str(WF / "main.nf"), "-profile", "test,stub", "-stub-run"]
    cmd += ["-params-file", str(_params(tmp)), "-work-dir", str(tmp / "work")]
    cmd += ["-with-trace", str(trace)]
    if resume:
        cmd.append("-resume")
    drop = {"SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE"}
    env = {k: v for k, v in os.environ.items() if k not in drop}
    env["NXF_ANSI_LOG"] = "false"
    res = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True, env=env, timeout=1200)
    assert res.returncode == 0, f"{' '.join(cmd)}\n{res.stdout}\n{res.stderr}"
    with trace.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    pc: dict[str, str] = {}
    for row in rows:
        m = re.fullmatch(r"(?P<proc>[^(]+?)\s*\((?P<tag>.*)\)\s*", row["name"])
        if m and m["proc"].split(":")[-1] == "PROMOTER_CALCULATOR":
            assert m["tag"] not in pc, f"two PROMOTER_CALCULATOR tasks tagged {m['tag']}"
            pc[m["tag"]] = row["status"]
    return pc


def _operon_pairs_region(work: Path, name: str) -> Path:
    for run in work.glob("*/*/.command.run"):
        head = run.read_text(errors="replace").splitlines()
        line = next((x for x in head if x.startswith("### name:")), "")
        if "OPERON_PAIRS" in line or "(operon_pairs)" in line:
            region = run.parent / "tu_upstream" / name
            assert region.is_file(), f"operon_pairs task {run.parent} has no {region.name}"
            return region
    raise AssertionError(f"no operon_pairs task directory under {work}")


def test_example_5_resume_reruns_only_the_changed_regions_batch(tmp_path: Path) -> None:
    first = _run(tmp_path, 1, resume=False)
    assert sorted(first) == BATCHES, f"{len(first)} PROMOTER_CALCULATOR tasks: {sorted(first)[:5]}"
    assert set(first.values()) == {"COMPLETED"}, first

    second = _run(tmp_path, 2, resume=True)
    assert second == dict.fromkeys(BATCHES, "CACHED"), (
        f"unchanged regions: every batch must be cached (membership depends only on TU ids): "
        f"{second}"
    )

    region = _operon_pairs_region(tmp_path / "work", "tu150.fa")
    region.write_text(">tu150\nACGTTGACAATTAATCATCGGCTCGTATAATGTGTGG\n", encoding="utf-8")
    third = _run(tmp_path, 3, resume=True)
    assert third == {
        "pc_tu000_tu099": "CACHED",
        "pc_tu100_tu199": "COMPLETED",
        "pc_tu200_tu249": "CACHED",
    }, third
