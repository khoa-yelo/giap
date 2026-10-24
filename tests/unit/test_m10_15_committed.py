"""The committed pilot report and its inputs (M10-15 acceptance criteria 2, 4, 6 and 8).

Spec: docs/specs/M10-15-pilot-gate.md, "The final report" (sections 1-10), "Acceptance
criteria". In 15a the committed ``REPORT.md``, ``metrics.json``, ``te1.json`` and
``report_notes.md`` do not exist yet: the re-render (2) and consistency (8) tests skip, and the
re-render runs on a fixture report instead (``test_m10_15_render_final.py``). The final template
(criterion 4) exists in 15a. The data-rule checkers of criterion 6 are defined here and tested
on hand-made documents, then applied to whichever committed files exist.

"No number" in the template (criterion 4, M10-05's rule for notes): after removing placeholders,
section references (``§15.6``), dates (``2026-10-24``), versions (``2.3.1``) and identifiers
(``ADR-0051``, ``PXD078459``, ``UP000000625``, ``GCA_000013265.1``), no run of 4 or more
digits remains.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path, PurePosixPath
from typing import Any

import pytest
from helpers.pilot_gate_m10_15 import render_final
from helpers.pilot_metrics_m10 import REPO, TEMPLATE

PILOT = REPO / "validation" / "pilot"
REPORT = PILOT / "REPORT.md"
METRICS = PILOT / "metrics" / "metrics.json"
TE1 = PILOT / "metrics" / "te1.json"
NOTES = PILOT / "report_notes.md"
RUNBOOK = PILOT / "reproduce.md"

PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")
LOCUS_TAG = re.compile(r"[A-Z][A-Z0-9]+_\d{3,}")
UNIPROT = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}")
MAX_IDS = 20
PLACEHOLDER = re.compile(r"\{\{\s*([^}\s]+)\s*\}\}")
SECTIONS = (
    "Status",
    "Data and arms",
    "Metrics",
    "TE1",
    "TE2",
    "Hypothes",
    "Not evaluated",
    "Runtime",
    "Reproducib",
    "Deviation",
)


def _tracked(*paths: str) -> list[str]:
    res = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-z", "--", *paths],
        capture_output=True,
        text=True,
        check=True,
    )
    return [p for p in res.stdout.split("\0") if p]


# ------------------------------------------------------------------ criterion 6 checkers


def _is_identifier(s: str) -> bool:
    return bool(LOCUS_TAG.fullmatch(s) or UNIPROT.fullmatch(s))


def identifier_lists(doc: Any, path: str = "") -> list[str]:
    """Paths of arrays holding more than ``MAX_IDS`` locus-tag or accession strings."""
    out: list[str] = []
    if isinstance(doc, dict):
        for k, v in doc.items():
            out += identifier_lists(v, f"{path}.{k}" if path else str(k))
    elif isinstance(doc, list):
        if sum(isinstance(x, str) and _is_identifier(x) for x in doc) > MAX_IDS:
            out.append(path or "<root>")
        for i, v in enumerate(doc):
            out += identifier_lists(v, f"{path}.{i}")
    return out


def test_identifier_list_checker_on_hand_made_documents() -> None:
    tags = [f"ECOLIN_{i:05d}" for i in range(21)]
    accs = ["P0A7V8", "Q9XYZ1", "A0A0H3JKH7"] * 7
    assert identifier_lists({"warnings": tags}) == ["warnings"]
    assert identifier_lists({"te7": {"differing": [{"ids": accs}]}}) == ["te7.differing.0.ids"]
    assert identifier_lists({"warnings": tags[:20]}) == []  # 20 is allowed
    paths = [f"tables/genes_{i}.tsv" for i in range(40)]  # table paths and contig names stay
    assert identifier_lists({"te7": {"differing": paths}, "contigs": ["chromosome"] * 30}) == []
    assert _is_identifier("ECOLIN_00012") and _is_identifier("P0A7V8")
    assert not _is_identifier("chromosome") and not _is_identifier("pending (M10-12)")


def _committed_texts() -> dict[str, str]:
    return {
        p.name: p.read_text(encoding="utf-8") for p in (REPORT, METRICS, TE1, NOTES) if p.is_file()
    }


def test_committed_report_files_hold_no_peptide_like_string() -> None:
    hits = {n: m.start() for n, t in _committed_texts().items() if (m := PEPTIDE_LIKE.search(t))}
    assert not hits, f"7+ amino-acid strings in committed report files: {hits}"


def test_committed_metrics_json_holds_no_identifier_list() -> None:
    if not METRICS.is_file():
        pytest.skip("validation/pilot/metrics/metrics.json is committed in 15b")
    bad = identifier_lists(json.loads(METRICS.read_text(encoding="utf-8")))
    assert not bad, f"arrays of more than {MAX_IDS} locus tags or accessions: {bad}"


def test_no_row_level_te_or_entrapment_output_is_tracked() -> None:
    names = {
        "te_real.json",
        "te7.json",
        "entrapment.json",
        "te7_files.tsv",
        "entrapment_discoveries.tsv",
        "ATTEMPTS.tsv",
    }
    bad = []
    for p in _tracked():
        pp = PurePosixPath(p)
        if pp.name in names or (
            pp.parts[:2] == ("validation", "pilot")
            and ("exceptions" in pp.name or (pp.suffix == ".tsv" and pp.parts[2:3] == ("metrics",)))
        ):
            bad.append(p)
    assert not bad, f"row-level pilot outputs are never committed: {bad}"


# ------------------------------------------------------------------ criterion 2: re-render


def test_committed_report_reproduces_byte_for_byte(tmp_path: Path) -> None:
    if not REPORT.is_file():
        pytest.skip("validation/pilot/REPORT.md is committed in 15b")
    out = tmp_path / "REPORT.md"
    res = render_final(METRICS, out, notes=NOTES if NOTES.is_file() else None)
    assert res.returncode == 0, res.stderr
    assert out.read_bytes() == REPORT.read_bytes(), "re-render differs from the committed report"


# ------------------------------------------------------------------ criterion 8: consistency


def _header(text: str) -> dict[str, str]:
    m = re.fullmatch(r"<!-- reproduce: (.*) -->", text.splitlines()[0].strip())
    assert m, "runbook header line"
    return dict(kv.split("=", 1) for kv in m.group(1).split() if "=" in kv)


def test_committed_metrics_commit_equals_the_freeze_commit() -> None:
    if not METRICS.is_file():
        pytest.skip("validation/pilot/metrics/metrics.json is committed in 15b")
    doc = json.loads(METRICS.read_text(encoding="utf-8"))
    head = _header(RUNBOOK.read_text(encoding="utf-8"))
    assert head.get("status") == "final"
    assert doc["giap_commit"] == head["freeze_commit"]


def test_committed_metrics_inputs_hash_to_committed_files() -> None:
    if not METRICS.is_file():
        pytest.skip("validation/pilot/metrics/metrics.json is committed in 15b")
    import hashlib

    doc = json.loads(METRICS.read_text(encoding="utf-8"))
    tracked = _tracked()
    bad = []
    for key, entry in doc.get("inputs", {}).items():
        path = str(entry.get("path", ""))
        rel = next((t for t in tracked if path.endswith("/" + t) or path == t), None)
        if rel is not None:
            digest = hashlib.sha256((REPO / rel).read_bytes()).hexdigest()
            if digest != entry["sha256"]:
                bad.append((key, rel))
    assert not bad, f"inputs whose committed file has another hash: {bad}"


# ------------------------------------------------------------------ criterion 4: final template


def _template() -> str:
    return TEMPLATE.read_text(encoding="utf-8")


def test_final_template_has_sections_1_to_10_in_order() -> None:
    heads = re.findall(r"^## (\d+)\. (.*)$", _template(), re.M)
    assert [int(n) for n, _ in heads] == list(range(1, 11)), heads
    for (_, title), word in zip(heads, SECTIONS, strict=True):
        assert word.lower() in title.lower(), (title, word)


def test_final_template_drops_pending_in_pass_2() -> None:
    assert "pending in pass 2" not in _template().lower()


def test_final_template_holds_no_number() -> None:
    text = PLACEHOLDER.sub("", _template())
    text = re.sub(r"§[\d.]+", "", text)
    text = re.sub(r"\b\d{4}-\d{2}-\d{2}\b", "", text)  # dates of decisions are not values
    text = re.sub(r"\b\d+\.\d+\.\d+\b", "", text)
    text = re.sub(r"\b[A-Z][A-Za-z]*[-_]?\d+(?:\.\d+)?\b", "", text)
    long_numbers = re.findall(r"\d{4,}", text)
    assert not long_numbers, f"numbers typed into the template: {long_numbers}"


@pytest.mark.parametrize("prefix", ["m6.", "m8.", "h4."])
def test_final_template_has_placeholders_for(prefix: str) -> None:
    assert any(p.startswith(prefix) for p in PLACEHOLDER.findall(_template())), prefix


def test_final_template_h2_deferred_beside_fdp_and_genome_inconsistent_share() -> None:
    t = _template()
    paths = PLACEHOLDER.findall(t)
    assert "m8.peptide.fdp_combined" in paths
    assert any(p.startswith("m2.A1.genome_inconsistent") for p in paths)
    h2 = [ln for ln in t.splitlines() if "H2" in ln or "Hypothesis 2" in ln]
    assert h2 and any("deferred" in ln for ln in h2), h2


def test_final_template_status_table_has_evidence_and_review_statement() -> None:
    t = _template()
    paths = PLACEHOLDER.findall(t)
    for i in range(5):
        assert f"checklist.{i}.status" in paths and f"checklist.{i}.evidence" in paths, i
    assert "Submitted for Khoa's review" in t
    low = t.lower()
    assert "reviewed by khoa" not in low  # the report never states that it was reviewed
