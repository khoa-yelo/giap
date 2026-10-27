"""M11-14w surfaces: AC 6 (``run_manifest.json`` ``warnings``) and AC 7 (the report summary
alert), plus the optional slow real-genome validation (docs/specs/M11-14w-wrong-genetic-code-
warning.md "Where it runs and where it surfaces" items 3-4, "Test plan"; GIAP_SPEC.md §9.2, §10.3).

Interface pinned here (beyond the spec text): the alert is the first element of the ``summary``
section after its heading, ``<div class="alert" data-warning="genome.genetic_code">``, whose text
(whitespace-normalised) is the check's message, HTML-escaped; it is shown only at level ``warn``.
"""

from __future__ import annotations

import dataclasses
import json
import os
from pathlib import Path
from typing import Any

import pytest
from helpers.genetic_code_fixture import synthetic_genome
from helpers.manifest_data import build, preflight_doc
from helpers.report_content import html, page, words
from helpers.report_fixtures import run_out

from giap.preflight.checks import Check

RULE_ID = "genome.genetic_code"
ALERT = f'data-warning="{RULE_ID}"'


def _check_dict(c: Check) -> dict[str, Any]:
    return dataclasses.asdict(c)


def _real(code: int) -> Check:
    from giap.preflight.genetic_code import genetic_code_check

    return genetic_code_check([synthetic_genome(7, 0.28, code)], 11)


# ---------------------------------------------------------------- AC 6: manifest


def test_ac6_manifest_lists_the_warning() -> None:
    c = _real(4)
    assert c.level == "warn"
    m = build(preflight=preflight_doc(checks=[_check_dict(c)]))
    assert f"preflight:{RULE_ID}" in m["warnings"]


def test_ac6_no_warning_without_it() -> None:
    c = _real(11)
    assert c.level == "ok"
    m = build(preflight=preflight_doc(checks=[_check_dict(c)]))
    assert f"preflight:{RULE_ID}" not in m["warnings"]


# ---------------------------------------------------------------- AC 7: report


def _with_checks(out: Path, checks: list[dict[str, Any]] | None) -> Path:
    path = out / "preflight" / "preflight.json"
    body = json.loads(path.read_text(encoding="utf-8"))
    if checks is None:
        body.pop("checks", None)
    else:
        body["checks"] = checks
    path.write_text(json.dumps(body, sort_keys=True, indent=1), encoding="utf-8")
    return out


def _row(level: str, message: str) -> dict[str, Any]:
    return {"rule_id": RULE_ID, "input": "genome", "level": level, "rule": "rule",
            "value": None, "message": message, "exit_code": 0}  # fmt: skip


def test_ac7_warning_alert_in_summary(tmp_path: Path) -> None:
    c = _real(4)
    out = _with_checks(run_out(tmp_path), [_check_dict(c)])
    p = page(out)
    summary = p.section("summary")
    alerts = summary.by_attr("data-warning", RULE_ID)
    assert len(alerts) == 1
    assert alerts[0].tag == "div" and alerts[0].has_class("alert")
    # styled inline so reports without the warning stay byte-identical (#423 review N1)
    assert "border-left" in (alerts[0].attrs.get("style") or "")
    assert words(alerts[0].text()) == words(c.message)
    assert len(p.root.by_attr("data-warning")) == 1


def test_ac7_alert_starts_the_summary_section(tmp_path: Path) -> None:
    out = _with_checks(run_out(tmp_path), [_check_dict(_real(4))])
    text = html(out)
    start = text.index('<section id="summary">')
    end = text.index("</section>", start)
    at = text.find(ALERT, start)
    assert start < at < end
    head = text.index("</h2>", start)
    between = text[head + len("</h2>") : at]
    assert between.strip().startswith("<div") and between.count("<") == 1, between


def test_ac7_message_is_escaped(tmp_path: Path) -> None:
    msg = "This genome probably uses genetic code 4 <script>alert(1)</script> & more"
    out = _with_checks(run_out(tmp_path), [_row("warn", msg)])
    text = html(out)
    assert "<script>alert(1)</script>" not in text
    assert "&lt;script&gt;" in text
    p = page(out)
    (alert,) = p.section("summary").by_attr("data-warning", RULE_ID)
    assert words(alert.text()) == words(msg)


@pytest.mark.parametrize("level", ["ok", "info"])
def test_ac7_no_alert_below_warn(tmp_path: Path, level: str) -> None:
    out = _with_checks(run_out(tmp_path), [_row(level, "not shown")])
    text = html(out)
    assert "data-warning" not in text
    assert not page(out).root.by_attr("data-warning")


def test_ac7_html_byte_identical_without_the_warning(tmp_path: Path) -> None:
    a = _with_checks(run_out(tmp_path / "a"), None)
    b = _with_checks(run_out(tmp_path / "b"), [_check_dict(_real(11))])
    c = _with_checks(run_out(tmp_path / "c"), [_row("info", "x")])
    assert html(a) == html(b) == html(c)
    w = _with_checks(run_out(tmp_path / "w"), [_check_dict(_real(4))])
    assert html(w) != html(a)


def test_ac7_other_warnings_are_not_shown(tmp_path: Path) -> None:
    """Q3: only this warning is shown in 1.0 (the general list is backlog)."""
    row = _row("warn", "other warning text")
    row["rule_id"] = "genome.contiguity"
    out = _with_checks(run_out(tmp_path), [row])
    assert "data-warning" not in html(out)


# ---------------------------------------------------------------- optional real-genome validation

REAL = [
    ("more_datasets/mhyo/raw/genome/GCF_000313635.1_ASM31363v1_genomic.fna.gz", 0.227, "warn"),
    ("data/genomes/k12_mg1655/GCA_000005845.2_ASM584v2_genomic.fna.gz", 0.881, "ok"),
    ("more_datasets/gsulf/raw/genome/GCF_000007985.2_ASM798v2_genomic.fna.gz", 0.823, "ok"),
    ("more_datasets/msmeg/raw/genome/GCF_000015005.1_ASM1500v1_genomic.fna.gz", 0.864, "ok"),
]


@pytest.mark.slow
@pytest.mark.parametrize(("rel", "ratio", "level"), REAL, ids=["mhyo", "k12", "gsulf", "msmeg"])
def test_real_genomes_match_the_spec_table(rel: str, ratio: float, level: str) -> None:
    """Spec "Thresholds" table, ratios to +-0.01; skipped when the local copies are absent."""
    root = os.environ.get("GIAP_ROOT")
    if not root or not (Path(root) / rel).is_file():
        pytest.skip(f"{rel} not under $GIAP_ROOT")
    from giap.io.fasta import read_genome_fasta
    from giap.preflight.genetic_code import genetic_code_check

    records = read_genome_fasta(Path(root) / rel)
    c = genetic_code_check([r.sequence for r in records], 11)
    assert c.level == level
    assert isinstance(c.value, dict)
    assert c.value["ratio"] == pytest.approx(ratio, abs=0.01)
