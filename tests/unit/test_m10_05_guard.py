"""Guard: no row-level pilot output is tracked; aggregates may be (M10-05 acceptance criterion 9).

Spec: docs/specs/M10-05-pilot-metrics.md, "Outputs" and "REPORT.md draft / What it may contain".
Khoa's answer to Q-K1 (2026-10-24): aggregates and example gene names may be committed. The guard
forbids only row-level files: peptide tables (the converted DIA-NN tables, any ``.parquet``),
the TE1 exception lists (``te1_exceptions*.tsv``, anywhere) and GIAP run tables (``genes.tsv``,
``peptide_hits.tsv``, ``compare_*.tsv``, ...) under ``validation/pilot/``; under
``validation/pilot/metrics/`` only the aggregates. An aggregate ``metrics.json``, ``te1.json``
or ``REPORT.md`` is allowed,
provided no peptide-like string (7+ amino-acid letters, ``metrics.PEPTIDE_LIKE``) appears in it.
Extends M10-03's guard on the converted tables.

The ``git ls-files`` tests pass while nothing forbidden is tracked; they exist to fail the commit
that would add one.
"""

from __future__ import annotations

import fnmatch
import re
import subprocess
from pathlib import PurePosixPath

from helpers.pilot_metrics_m10 import REPO, TEMPLATE

#: row-level outputs, forbidden anywhere in the checkout
FORBIDDEN_NAMES = ("te1_exceptions*.tsv", "te1_exceptions*.tsv.gz")
#: row-level tables, forbidden under validation/pilot/ (peptide tables and GIAP run tables)
FORBIDDEN_UNDER_PILOT = (
    "a1.tsv*", "a2.tsv*", "a3.tsv*", "*peptide*.tsv*", "report*.tsv*", "*.parquet",
    "genes.tsv*", "peptide_hits.tsv*", "unmapped.tsv*", "detection_summary.tsv*",
    "proteins.tsv*", "compare_*.tsv*",
)  # fmt: skip
#: aggregate outputs that may be committed under validation/pilot/ (Q-K1, 2026-10-24;
#: a3_feasibility.json: M10-08 deliverable 10, aggregates only)
AGGREGATES = ("metrics.json", "te1.json", "REPORT.md", "a3_feasibility.json")
#: the peptide-like pattern of ``metrics.check_notes`` (spec: 7 or more amino-acid letters)
PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")
#: top-level keys of metrics.json (spec "Outputs"; giap_commit: review S4; te_real and te7:
#: M10-06, m8: M10-07, both stacked on this branch)
TOP_KEYS = {"schema", "inputs", "spec", "giap_commit", "m1", "m2", "m3", "m4", "m5", "m6", "m7",
            "m8", "te1", "te1_a1_diagnostic", "te_real", "te7", "h1", "h3", "h4", "checklist",
            "evidence", "warnings"}  # fmt: skip


def _tracked(*paths: str) -> list[str]:
    res = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-z", "--", *paths],
        capture_output=True,
        text=True,
        check=True,
    )
    return [p for p in res.stdout.split("\0") if p]


def _forbidden(path: str) -> bool:
    pp = PurePosixPath(path)
    if any(fnmatch.fnmatchcase(pp.name, pat) for pat in FORBIDDEN_NAMES):
        return True
    if pp.parts[:2] != ("validation", "pilot"):
        return False
    if pp.parts[:3] == ("validation", "pilot", "metrics") and pp.name not in AGGREGATES:
        return True
    return any(fnmatch.fnmatchcase(pp.name, pat) for pat in FORBIDDEN_UNDER_PILOT)


def test_no_row_level_output_is_tracked() -> None:
    bad = [p for p in _tracked() if _forbidden(p)]
    assert not bad, f"row-level pilot outputs are never committed (Q-K1): {bad}"


def test_tracked_aggregates_hold_no_peptide_like_string() -> None:
    hits = []
    for p in _tracked("validation/pilot/"):
        if PurePosixPath(p).name in AGGREGATES or fnmatch.fnmatchcase(p, "validation/pilot/*/*.md"):
            text = (REPO / p).read_text(encoding="utf-8")
            if m := PEPTIDE_LIKE.search(text):
                hits.append((p, m.start()))
    assert not hits, f"peptide-like strings in committed aggregates: {hits}"


def test_the_patterns_forbid_rows_and_allow_aggregates() -> None:
    """Self-check of the patterns against the names ``metrics.py compute`` / ``render`` write."""
    assert _forbidden("validation/pilot/metrics/te1_exceptions.tsv")
    assert _forbidden("local/review/te1_exceptions.tsv")
    row_level = (
        "validation/pilot/metrics/peptides_a2.tsv.gz",
        "validation/pilot/metrics/m1_per_run.json",
        "validation/pilot/a2.tsv.gz",
        "validation/pilot/runs/genes.tsv",
        "validation/pilot/peptide_hits.tsv",
        "validation/pilot/compare_detection.tsv",
        "validation/pilot/a2.parquet",
    )
    for path in row_level:
        assert _forbidden(path), path
    for aggregate in ("validation/pilot/metrics/metrics.json", "validation/pilot/metrics/te1.json",
                      "validation/pilot/REPORT.md", "validation/pilot/metrics.py",
                      "validation/pilot/report_template.md",
                      "validation/pilot/db_a2/db_id_map.tsv"):  # fmt: skip
        assert not _forbidden(aggregate), aggregate
    assert not _forbidden("tests/data/recorded/genes.tsv")  # fixtures elsewhere are not pilot rows


def test_the_peptide_pattern_matches_metrics_py() -> None:
    from helpers.pilot_metrics_m10 import load_metrics

    assert load_metrics().PEPTIDE_LIKE.pattern == PEPTIDE_LIKE.pattern
    assert PEPTIDE_LIKE.search("x QQWWRAK y") and not PEPTIDE_LIKE.search("x QQWWRA y")


def test_the_template_holds_no_peptide_like_string() -> None:
    m = PEPTIDE_LIKE.search(TEMPLATE.read_text(encoding="utf-8"))
    assert m is None, m


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
