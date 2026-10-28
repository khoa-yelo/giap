"""M10-06b: the TE7 half of the ``metrics.py`` hook (``--te7``), the TE7 diagnosis rule of
``render``, the TE7 row of the template and the committed TE7 pieces (M10-06 acceptance
criteria 9 and 11, TE7 half). The TE2-TE5 half is ``test_m10_06_hook.py`` (M10-06a).

Spec: docs/specs/M10-06-te-real-and-te7.md ("Outputs" > "metrics.py hook": "TE7 is shown in its
own report section (it is not a §15.6 bullet)"; "Investigation of failures": "for TE7 when it
fails"). ``te_real.json`` comes from ``te_real.py compute`` on the synthetic pilot, ``te7.json``
from ``te_real.py te7`` on a hand-made pair; scenarios edit only status fields.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path
from typing import Any

import pytest
from helpers.pilot_metrics_m10 import TEMPLATE, run_cli, sim_world
from helpers.te7_m10 import base_tree, load_te7, make_repo, write_run
from helpers.te_real_m10 import REPO, build_te_pilot, run_te_real, te_compute_args

SECTION = "True examples (TE2\u2013TE5) and determinism (TE7)"


@pytest.fixture(scope="module")
def inputs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    root = tmp_path_factory.mktemp("m10_06b_hook")
    world = sim_world(root / "sim")
    p, contam = build_te_pilot(root / "pilot", world)
    te_out = root / "te_metrics"
    res = run_te_real(*te_compute_args(p, contam, te_out), cwd=p.root)
    assert res.returncode == 0, res.stdout + res.stderr
    shas = make_repo(root / "repo")
    ref = write_run(
        root / "runs" / "out_a2", base_tree(), commit=shas["base"], outdir_label="out_a2"
    )
    rep = write_run(root / "runs" / "rep1", base_tree(), commit=shas["base"], outdir_label="rep1")
    res = run_te_real(
        "te7", "--ref", str(ref), "--rep", str(rep), "--repo", str(root / "repo"),
        "--out", str(te_out),
    )  # fmt: skip
    assert res.returncode == 0, res.stdout + res.stderr
    te_real = json.loads((te_out / "te_real.json").read_text(encoding="utf-8"))
    te7 = json.loads((te_out / "te7.json").read_text(encoding="utf-8"))
    assert te7["gippa"]["status"] == "pass"
    return {"pilot": p, "te_real": te_real, "te7": te7}


def _files(inputs: dict[str, Any], tmp: Path, *, te7_fail: bool = False) -> tuple[Path, Path]:
    te7 = json.loads(json.dumps(inputs["te7"]))
    if te7_fail:
        te7["gippa"]["status"] = "fail"
    a, b = tmp / "te_real.json", tmp / "te7.json"
    a.write_text(json.dumps(inputs["te_real"], sort_keys=True, indent=2) + "\n", encoding="utf-8")
    b.write_text(json.dumps(te7, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return a, b


def _metrics(inputs: dict[str, Any], tmp: Path, **kw: Any) -> tuple[Path, dict[str, Any]]:
    a, b = _files(inputs, tmp, **kw)
    out = tmp / "m"
    p = inputs["pilot"]
    res = run_cli(*p.compute_args(out, "--te-real", str(a), "--te7", str(b)), cwd=p.root)
    assert res.returncode == 0, res.stdout + res.stderr
    return out / "metrics.json", json.loads((out / "metrics.json").read_text(encoding="utf-8"))


def _render(metrics: Path, out: Path, notes: str | None) -> Any:
    args = ["render", "--metrics", str(metrics), "--out", str(out)]
    if notes is not None:
        nf = out.with_suffix(".notes.md")
        nf.write_text(notes, encoding="utf-8")
        args += ["--notes", str(nf)]
    return run_cli(*args)


def test_te7_json_is_copied_unchanged(inputs: dict[str, Any], tmp_path: Path) -> None:
    _, doc = _metrics(inputs, tmp_path)
    assert doc["te7"] == inputs["te7"]


def test_te7_is_not_a_bullet(inputs: dict[str, Any], tmp_path: Path) -> None:
    _, doc = _metrics(inputs, tmp_path, te7_fail=True)
    e = doc["checklist"][4]
    assert (e["status"] if isinstance(e, dict) else e) == "pass"  # TE7 has its own section
    assert len(doc["checklist"]) == 5


def test_render_exit_1_for_te7_fail_without_diagnosis(
    inputs: dict[str, Any], tmp_path: Path
) -> None:
    metrics, _ = _metrics(inputs, tmp_path, te7_fail=True)
    res = _render(metrics, tmp_path / "REPORT.md", None)
    assert res.returncode == 1 and "TE7" in res.stderr, res.stdout + res.stderr
    notes = "## Diagnosis: TE7\n\nCause class: pipeline_bug. Linked backlog entry.\n"
    res = _render(metrics, tmp_path / "REPORT.md", notes)
    assert res.returncode == 0, res.stdout + res.stderr


def test_render_te7_pass_needs_no_diagnosis(inputs: dict[str, Any], tmp_path: Path) -> None:
    metrics, _ = _metrics(inputs, tmp_path)
    report = tmp_path / "REPORT.md"
    res = _render(metrics, report, None)
    assert res.returncode == 0, res.stdout + res.stderr
    assert SECTION in report.read_text(encoding="utf-8")


def test_template_has_the_te7_row() -> None:
    text = TEMPLATE.read_text(encoding="utf-8")
    assert SECTION in text
    assert "TE7" in text.split(SECTION, 1)[1]
    paths = re.findall(r"\{\{\s*([^}\s]+)\s*\}\}", text)
    assert any(p.startswith("te7.") for p in paths)


def test_template_pending_section_drops_te7() -> None:
    """M10-15 removed the draft's "Pending in pass 2" section; TE7 renders in its own."""
    text = TEMPLATE.read_text(encoding="utf-8")
    assert "pending in pass 2" not in text.lower()
    assert "{{te7.gippa.status}}" in text


def test_committed_te7_pieces_are_tracked() -> None:
    load_te7()  # te7.py exists (M10-06b)
    res = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "--", "validation/pilot/"],
        capture_output=True, text=True, check=True,
    )  # fmt: skip
    tracked = set(res.stdout.split())
    assert "validation/pilot/te7.py" in tracked
    assert "validation/pilot/te7_known_columns.tsv" in tracked
