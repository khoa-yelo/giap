"""Guard: no pilot metric output is tracked (M10-05 acceptance criterion 9).

Spec: docs/specs/M10-05-pilot-metrics.md, "Outputs" and "REPORT.md draft / What it may contain";
Khoa's rule of 2026-10-23 (nothing derived from the DIA-NN reports is committed), read by default
(Q-K1, lead decision 3) as covering ``metrics.json``, ``te1.json``, ``te1_exceptions.tsv`` and
the rendered ``REPORT.md``: they live under ``$GIAP_ROOT/pilot/metrics/`` only. Extends M10-03's
guard on the converted tables.

The ``git ls-files`` tests pass while nothing forbidden is tracked; they exist to fail the commit
that would add one. The template tests fail until ``validation/pilot/report_template.md`` (the
committed, number-free template) exists.
"""

from __future__ import annotations

import fnmatch
import re
import subprocess
from pathlib import PurePosixPath

from helpers.pilot_metrics_m10 import REPO, TEMPLATE

FORBIDDEN_NAMES = ("metrics.json", "te1.json", "te1_exceptions.tsv")
#: top-level keys of metrics.json (spec "Outputs")
TOP_KEYS = {"schema", "inputs", "spec", "m1", "m2", "m3", "m4", "m5", "m7", "te1",
            "te1_a1_diagnostic", "h1", "h3", "checklist", "warnings"}  # fmt: skip


def _tracked(*paths: str) -> list[str]:
    res = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-z", "--", *paths],
        capture_output=True,
        text=True,
        check=True,
    )
    return [p for p in res.stdout.split("\0") if p]


def test_no_rendered_report_is_tracked() -> None:
    assert not _tracked("validation/pilot/REPORT.md"), "REPORT.md stays on FarmShare (Q-K1)"


def test_nothing_under_validation_pilot_metrics() -> None:
    assert not _tracked("validation/pilot/metrics/"), _tracked("validation/pilot/metrics/")


def test_no_metric_file_is_tracked_anywhere() -> None:
    bad = [
        p
        for p in _tracked()
        if any(fnmatch.fnmatchcase(PurePosixPath(p).name, pat) for pat in FORBIDDEN_NAMES)
    ]
    assert not bad, f"pilot metric outputs are never committed (Khoa 2026-10-23): {bad}"


def test_the_patterns_catch_the_compute_outputs() -> None:
    """Self-check of the patterns against the names ``metrics.py compute`` writes."""
    for name in ("metrics.json", "te1.json", "te1_exceptions.tsv"):
        assert any(fnmatch.fnmatchcase(name, pat) for pat in FORBIDDEN_NAMES), name
    assert not any(fnmatch.fnmatchcase("metrics.py", pat) for pat in FORBIDDEN_NAMES)


def test_the_template_is_committed() -> None:
    assert "validation/pilot/report_template.md" in _tracked("validation/pilot/")


def test_template_placeholders_are_dotted_paths_into_metrics_json() -> None:
    text = TEMPLATE.read_text(encoding="utf-8")
    paths = re.findall(r"\{\{\s*([^}\s]+)\s*\}\}", text)
    assert paths, "the template renders no value"
    unknown = sorted({p for p in paths if p.split(".")[0] not in TOP_KEYS})
    assert not unknown, f"placeholders outside metrics.json: {unknown}"


def test_template_states_the_te1_limits_and_the_h1_wording() -> None:
    """S6: "says nothing about identification error"; H1 "not supported by the
    pre-registered rule", never "refuted"; the CI is "nominal"."""
    low = TEMPLATE.read_text(encoding="utf-8").lower()
    assert "identification error" in low
    assert "not supported by the pre-registered rule" in low
    assert "refuted" not in low
    assert "nominal" in low
