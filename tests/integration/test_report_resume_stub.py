"""M9-26 example 3 and work item 5: ``-resume`` caches ``PLOTS`` and ``REPORT_HTML`` when no table
changed and reruns both when one table changed (docs/specs/M9-26-scale-and-determinism.md;
GIAP_SPEC.md §6.3 ``-resume``, §3.8 "Resume").

Stub runs of the REPORT subworkflow (``workflow/subworkflows/local/report.nf``) through a small
driver script written into ``tmp_path`` that feeds it ``[key, path]`` pairs (the same keys as
``workflow/tests/report.nf.test``: preflight, pairs, regions_tsv, events) from copies of the
``workflow/tests/data`` files, with ``-profile test,stub -stub-run`` and the repository's
``workflow/nextflow.config``. Every run of one test shares one work directory; the trace
(``-with-trace``) of each run gives the task status (``CACHED`` / ``COMPLETED``) per process.

Interface choices (M9-26 test author): the spec asks for an nf-test; this test author may not
write under ``workflow/``, so the same check is a pytest that drives Nextflow directly (as
tests/integration/test_nf_stub_run.py does for the whole pipeline). "A table changed" is a
change of the bytes of ``regions.tsv`` (one appended row), the way an upstream rerun would
deliver a new file. Marked ``needs_nextflow``; local executor inside the session allocation.
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.needs_nextflow

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"
DATA = WF / "tests" / "data"
INPUTS = {
    "preflight": DATA / "genetable" / "preflight.json",
    "pairs": DATA / "annotate" / "operon_pairs.tsv",
    "regions_tsv": DATA / "annotate" / "regions.tsv",
    "events": DATA / "annotate" / "annotation.json",
}
DRIVER = """\
nextflow.enable.dsl = 2
include {{ REPORT }} from '{report}'

workflow {{
    REPORT(channel.of({pairs}))
}}
"""


def _driver(tmp: Path) -> tuple[Path, dict[str, Path]]:
    data = tmp / "data"
    data.mkdir()
    copies = {}
    for key, src in INPUTS.items():
        dst = data / src.name
        shutil.copyfile(src, dst)
        copies[key] = dst
    pairs = ", ".join(f"['{k}', file('{p}')]" for k, p in copies.items())
    script = tmp / "report_driver.nf"
    report = WF / "subworkflows" / "local" / "report.nf"
    script.write_text(DRIVER.format(report=report, pairs=pairs), encoding="utf-8")
    return script, copies


def _run(tmp: Path, script: Path, n: int, resume: bool) -> dict[str, str]:
    """One stub run; returns ``{process: status}`` for PLOTS and REPORT_HTML."""
    out = tmp / "OUT"
    trace = tmp / f"trace{n}.tsv"
    (tmp / "params.json").write_text(json.dumps({"outdir": str(out)}), encoding="utf-8")
    cmd = ["nextflow", "-q", "-c", str(WF / "nextflow.config"), "run", str(script)]
    cmd += ["-profile", "test,stub", "-stub-run", "-params-file", str(tmp / "params.json")]
    cmd += ["-work-dir", str(tmp / "work"), "-with-trace", str(trace)]
    if resume:
        cmd.append("-resume")
    env = {**os.environ, "NXF_ANSI_LOG": "false"}
    res = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True, env=env, timeout=900)
    assert res.returncode == 0, f"{' '.join(cmd)}\n{res.stdout}\n{res.stderr}"
    with trace.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    status: dict[str, str] = {}
    for row in rows:
        proc = row["name"].split(" (")[0].rsplit(":", 1)[-1]
        status[proc] = row["status"]
    assert {"PLOTS", "REPORT_HTML"} <= status.keys(), [r["name"] for r in rows]
    assert (out / "report" / "report.html").is_file(), "report/report.html not published"
    return status


def test_example_3_resume_without_change_caches_plots_and_report(tmp_path: Path) -> None:
    script, _ = _driver(tmp_path)
    first = _run(tmp_path, script, 1, resume=False)
    assert first["PLOTS"] == "COMPLETED" and first["REPORT_HTML"] == "COMPLETED", first
    second = _run(tmp_path, script, 2, resume=True)
    assert second["PLOTS"] == "CACHED", second
    assert second["REPORT_HTML"] == "CACHED", second


def test_resume_after_one_table_changed_reruns_plots_and_report(tmp_path: Path) -> None:
    script, copies = _driver(tmp_path)
    _run(tmp_path, script, 1, resume=False)
    regions = copies["regions_tsv"]
    lines = regions.read_text(encoding="utf-8").splitlines()
    assert len(lines) >= 2, "regions.tsv fixture has no data row to duplicate"
    regions.write_text("\n".join([*lines, lines[-1]]) + "\n", encoding="utf-8")
    second = _run(tmp_path, script, 2, resume=True)
    assert second["PLOTS"] == "COMPLETED", second
    assert second["REPORT_HTML"] == "COMPLETED", second
