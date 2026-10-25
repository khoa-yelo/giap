"""M10-06a: the ``metrics.py`` hook for TE2-TE5 (``--te-real``), the ``render`` diagnosis rule,
the template section and the guard (M10-06 acceptance criteria 9 and 11; §15.6 bullet 5; review
round 1 N2). The TE7 half (``--te7``, the TE7 row, a TE7 diagnosis) is
``test_m10_06_te7_hook.py`` (M10-06b).

Spec: docs/specs/M10-06-te-real-and-te7.md ("Investigation of failures", "Outputs" > "metrics.py
hook", "Committed"); Khoa's rule of 2026-10-24 (aggregates and named example genes may be
committed, row-level DIA-NN content may not: the two TSVs are never tracked; a tracked
``te_real.json`` / ``te7.json`` must hold no sequence).

``te_real.json`` comes from ``te_real.py compute`` on the synthetic pilot; scenarios edit only its
``status`` fields. Interface fixed here: ``metrics.py compute`` copies it under ``te_real`` of
``metrics.json``; checklist bullet 5 is ``checklist.4.status``; a TE is "not passed" when its
``status`` is not ``pass`` (TE5 by its combined ``te5.status``); ``render`` exit 1 names the TE
(``TE3``) on stderr; a diagnosis section names exactly one cause class (N2).
"""

from __future__ import annotations

import fnmatch
import json
import re
import subprocess
from pathlib import Path, PurePosixPath
from typing import Any

import pytest
from helpers.pilot_metrics_m10 import TEMPLATE, run_cli, sim_world
from helpers.te_real_m10 import (
    CAUSE_CLASSES,
    REPO,
    build_te_pilot,
    run_te_real,
    te_compute_args,
)

#: the 06a heading may stop here; 06b completes it with " and determinism (TE7)"
SECTION_PREFIX = "True examples (TE2\u2013TE5)"


@pytest.fixture(scope="module")
def inputs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    root = tmp_path_factory.mktemp("m10_06_hook")
    world = sim_world(root / "sim")
    p, contam = build_te_pilot(root / "pilot", world)
    te_out = root / "te_metrics"
    res = run_te_real(*te_compute_args(p, contam, te_out), cwd=p.root)
    assert res.returncode == 0, res.stdout + res.stderr
    te_real = json.loads((te_out / "te_real.json").read_text(encoding="utf-8"))
    assert te_real["te2"]["status"] == te_real["te3"]["status"] == "pass"
    return {"root": root, "pilot": p, "te_real": te_real}


def _variant(inputs: dict[str, Any], tmp: Path, *, fail: tuple[str, ...] = ()) -> Path:
    """A copy of te_real.json with the named TEs set to ``fail``."""
    te_real = json.loads(json.dumps(inputs["te_real"]))
    for te in fail:
        te_real[te]["status"] = "fail"
    a = tmp / "te_real.json"
    a.write_text(json.dumps(te_real, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return a


def _metrics(inputs: dict[str, Any], out: Path, *extra: str) -> dict[str, Any]:
    p = inputs["pilot"]
    res = run_cli(*p.compute_args(out, *extra), cwd=p.root)
    assert res.returncode == 0, res.stdout + res.stderr
    doc: dict[str, Any] = json.loads((out / "metrics.json").read_text(encoding="utf-8"))
    return doc


def _bullet5(doc: dict[str, Any]) -> str:
    e = doc["checklist"][4]
    return str(e["status"] if isinstance(e, dict) else e)


def _render(metrics: Path, out: Path, notes: str | None) -> Any:
    args = ["render", "--metrics", str(metrics), "--out", str(out)]
    if notes is not None:
        nf = out.with_suffix(".notes.md")
        nf.write_text(notes, encoding="utf-8")
        args += ["--notes", str(nf)]
    return run_cli(*args)


# ------------------------------------------------------------------ criterion 9: bullet 5


def test_bullet_5_pending_without_te_real(inputs: dict[str, Any], tmp_path: Path) -> None:
    assert _bullet5(_metrics(inputs, tmp_path / "m")) == "pending (M10-06)"


def test_bullet_5_pass_with_four_passes(inputs: dict[str, Any], tmp_path: Path) -> None:
    a = _variant(inputs, tmp_path)
    doc = _metrics(inputs, tmp_path / "m", "--te-real", str(a))
    assert _bullet5(doc) == "pass"
    assert doc["te_real"] == inputs["te_real"]  # copied unchanged


@pytest.mark.parametrize(
    ("fail", "k"), [(("te3",), 1), (("te2", "te5"), 2), (("te2", "te3", "te4", "te5"), 4)]
)
def test_bullet_5_reported_k_not_passed(
    inputs: dict[str, Any], tmp_path: Path, fail: tuple[str, ...], k: int
) -> None:
    a = _variant(inputs, tmp_path, fail=fail)
    doc = _metrics(inputs, tmp_path / "m", "--te-real", str(a))
    assert _bullet5(doc) == f"reported ({k} not passed, explained)"
    assert len(doc["checklist"]) == 5


def test_not_evaluable_counts_as_not_passed(inputs: dict[str, Any], tmp_path: Path) -> None:
    a = _variant(inputs, tmp_path)
    doc = json.loads(a.read_text(encoding="utf-8"))
    doc["te4"]["status"] = "not_evaluable"
    a.write_text(json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    m = _metrics(inputs, tmp_path / "m", "--te-real", str(a))
    assert _bullet5(m) == "reported (1 not passed, explained)"


# ------------------------------------------------------------------ criterion 9: render rule


def _computed(inputs: dict[str, Any], tmp: Path, **kw: Any) -> Path:
    a = _variant(inputs, tmp, **kw)
    out = tmp / "m"
    _metrics(inputs, out, "--te-real", str(a))
    return out / "metrics.json"


def test_render_all_pass_without_notes(inputs: dict[str, Any], tmp_path: Path) -> None:
    metrics = _computed(inputs, tmp_path)
    report = tmp_path / "REPORT.md"
    res = _render(metrics, report, None)
    assert res.returncode == 0, res.stdout + res.stderr
    text = report.read_text(encoding="utf-8")
    assert SECTION_PREFIX in text
    assert "**missing:" not in text


def test_render_exit_1_without_diagnosis(inputs: dict[str, Any], tmp_path: Path) -> None:
    metrics = _computed(inputs, tmp_path, fail=("te3",))
    res = _render(metrics, tmp_path / "REPORT.md", None)
    assert res.returncode == 1 and "TE3" in res.stderr, res.stdout + res.stderr
    res = _render(metrics, tmp_path / "REPORT.md", "Nothing to say about the pilot.\n")
    assert res.returncode == 1 and "TE3" in res.stderr, res.stdout + res.stderr


@pytest.mark.parametrize("cause", CAUSE_CLASSES)
def test_render_ok_with_diagnosis_and_cause_class(
    inputs: dict[str, Any], tmp_path: Path, cause: str
) -> None:
    metrics = _computed(inputs, tmp_path, fail=("te3",))
    notes = f"## Diagnosis: TE3\n\nCause class: {cause}. Evidence: the aggregates above.\n"
    res = _render(metrics, tmp_path / "REPORT.md", notes)
    assert res.returncode == 0, res.stdout + res.stderr
    assert "## Diagnosis: TE3" in (tmp_path / "REPORT.md").read_text(encoding="utf-8")


def test_render_exit_1_when_the_section_names_no_cause_class(
    inputs: dict[str, Any], tmp_path: Path
) -> None:
    metrics = _computed(inputs, tmp_path, fail=("te2", "te3"))
    notes = (
        "## Diagnosis: TE2\n\nCause class: annotation. A start the model misses.\n\n"
        "## Diagnosis: TE3\n\nWe looked at it.\n"
    )
    res = _render(metrics, tmp_path / "REPORT.md", notes)
    assert res.returncode == 1 and "TE3" in res.stderr, res.stdout + res.stderr


@pytest.mark.parametrize(
    "line",
    ["Cause class: annotation, biology.", "Cause class: search or definition."],
    ids=["two_classes", "either_or"],
)
def test_render_exit_1_when_the_section_names_two_cause_classes(
    inputs: dict[str, Any], tmp_path: Path, line: str
) -> None:
    """N2: "names one cause class" means exactly one of the five, not one or more."""
    metrics = _computed(inputs, tmp_path, fail=("te3",))
    notes = f"## Diagnosis: TE3\n\n{line} Evidence: the aggregates above.\n"
    res = _render(metrics, tmp_path / "REPORT.md", notes)
    assert res.returncode == 1 and "TE3" in res.stderr, res.stdout + res.stderr


def test_vacuous_te2_pass_needs_no_diagnosis(inputs: dict[str, Any], tmp_path: Path) -> None:
    a = _variant(inputs, tmp_path)
    doc = json.loads(a.read_text(encoding="utf-8"))
    doc["te2"]["vacuous"] = True
    a.write_text(json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    out = tmp_path / "m"
    m = _metrics(inputs, out, "--te-real", str(a))
    assert _bullet5(m) == "pass"
    res = _render(out / "metrics.json", tmp_path / "REPORT.md", None)
    assert res.returncode == 0, res.stdout + res.stderr


# ------------------------------------------------------------------ the template


def test_template_has_the_te_section() -> None:
    text = TEMPLATE.read_text(encoding="utf-8")
    assert SECTION_PREFIX in text
    for te in ("TE2", "TE3", "TE4", "TE5"):
        assert te in text.split(SECTION_PREFIX, 1)[1], te
    paths = re.findall(r"\{\{\s*([^}\s]+)\s*\}\}", text)
    assert any(p.startswith("te_real.") for p in paths)


def test_template_pending_section_drops_te2_to_te5() -> None:
    """M10-15 removed the draft's "Pending in pass 2" section; TE2-TE5 render in their own."""
    text = TEMPLATE.read_text(encoding="utf-8")
    assert "pending in pass 2" not in text.lower()
    for te in ("te2", "te3", "te4", "te5"):
        assert f"{{{{te_real.{te}.status}}}}" in text, te


# ------------------------------------------------------------------ criterion 11: guard

FORBIDDEN_NAMES = ("te_real_exceptions.tsv", "te7_files.tsv")
JSON_NAMES = ("te_real.json", "te7.json")
PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")


def _tracked(*paths: str) -> list[str]:
    res = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-z", "--", *paths],
        capture_output=True,
        text=True,
        check=True,
    )
    return [p for p in res.stdout.split("\0") if p]


def _match(path: str, names: tuple[str, ...]) -> bool:
    return any(fnmatch.fnmatchcase(PurePosixPath(path).name, n) for n in names)


def test_row_level_outputs_are_never_tracked() -> None:
    bad = [p for p in _tracked() if _match(p, FORBIDDEN_NAMES)]
    assert not bad, f"row-level pilot outputs are never committed (Khoa 2026-10-24): {bad}"


def test_tracked_te_json_holds_no_sequence() -> None:
    """Aggregates may be committed (2026-10-24); a tracked te_real.json / te7.json holds no
    peptide-like string."""
    for p in (p for p in _tracked() if _match(p, JSON_NAMES)):
        text = (REPO / p).read_text(encoding="utf-8")
        assert not PEPTIDE_LIKE.search(text), p


def test_guard_patterns_catch_the_outputs() -> None:
    for name in FORBIDDEN_NAMES:
        assert _match(f"x/{name}", FORBIDDEN_NAMES)
    assert not _match("validation/pilot/te_real.py", FORBIDDEN_NAMES)
    assert not _match("validation/pilot/te7_known_columns.tsv", FORBIDDEN_NAMES)


def test_committed_pieces_are_tracked() -> None:
    tracked = set(_tracked("validation/pilot/"))
    assert "validation/pilot/te_real.py" in tracked
