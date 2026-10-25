"""What M10-08b commits: the verdict JSON, the README section and the pass-only provenance note.

Spec: docs/specs/M10-08-a3-feasibility.md, deliverables 10, 12 and 13, the verdict JSON fields
("aggregates only"; "No sequence and no per-peptide value appears in the file"), the decision rule
L6 and its exit/verdict semantics, and acceptance criteria 3 and 12. Reads tracked files only.
The ``.ssl`` guard is ``tests/unit/test_m10_09_guard.py`` and the binary-artefact guard is
``tests/unit/test_m10_08_committed.py``; they are not duplicated here.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pytest
from helpers.a3_eval_kit import (
    COUNT_FIELDS,
    JSON_FIELDS,
    REASONS,
    RULE,
    SSL_COLUMNS,
    VERDICT_REASONS,
    sequence_like_strings,
)
from helpers.a3_kit import PILOT, is_canonical_json, require

from giap.stats.cutoffs import wilson_ci

METRICS = PILOT / "metrics" / "a3_feasibility.json"
README = PILOT / "README.md"
RECORDED = PILOT.parents[1] / "tests" / "data" / "recorded" / "cascadia"
PROVENANCE = RECORDED / "0.0.7" / "PROVENANCE.md"
HEX64 = re.compile(r"^[0-9a-f]{64}$")


@pytest.fixture(scope="module")
def verdict() -> dict[str, Any]:
    assert METRICS.is_file(), (
        "validation/pilot/metrics/a3_feasibility.json does not exist yet "
        "(M10-08b: deliverable 10, the committed verdict)"
    )
    text = METRICS.read_text()
    assert is_canonical_json(text), "sorted keys, 2-space indent, trailing newline"
    j = json.loads(text)
    assert isinstance(j, dict)
    return j


def test_every_field_and_aggregates_only(verdict: dict[str, Any]) -> None:
    assert set(verdict) >= JSON_FIELDS, sorted(JSON_FIELDS - set(verdict))
    assert verdict["schema_version"] == 1 and verdict["run"] == "WT_1"
    assert verdict["cascadia_version"] == "0.0.7"
    assert verdict["rule"] == RULE
    assert not sequence_like_strings(verdict), sequence_like_strings(verdict)
    for k in COUNT_FIELDS:
        assert verdict[k] is None or (isinstance(verdict[k], int) and verdict[k] >= 0), k
    for k in ("matches_by_replicon", "length_hist", "mzml_check"):
        v = verdict[k]
        assert v is None or isinstance(v, dict), k
    for v in (verdict["matches_by_replicon"] or {}).values():
        assert isinstance(v, int)
    for k, v in (verdict["length_hist"] or {}).items():
        assert 7 <= int(k) <= 30 and isinstance(v, int)
    for k, v in (verdict["mzml_check"] or {}).items():
        assert re.fullmatch(r"M[1-8]", k) and isinstance(v, bool)
    assert isinstance(verdict["job_ids"], list)
    assert all(isinstance(x, int | str) for x in verdict["job_ids"])
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", str(verdict["date"]))
    for k in ("weights_sha256", "mzml_sha256", "genome_sha256"):
        assert verdict[k] is None or HEX64.match(verdict[k]), k
    assert verdict["route"] in ("B2", "B1", None)
    # every list is a short aggregate (the CI pair or job IDs), never a per-peptide table
    assert verdict["share_ci95"] is None or len(verdict["share_ci95"]) == 2
    assert len(verdict["job_ids"]) <= 20


def test_verdict_is_consistent_with_rule_l6(verdict: dict[str, Any]) -> None:
    assert verdict["verdict"] in ("pass", "fail")
    assert verdict["reason"] in VERDICT_REASONS
    if verdict["reason"] in REASONS:  # early fail: no counts
        assert verdict["verdict"] == "fail" and verdict["evidence"]
        assert all(verdict[k] is None for k in COUNT_FIELDS)
        assert verdict["fails"] == {"too_few_peptides": None, "low_genome_match": None}
        return
    n_pass, n_top, n_match = verdict["n_pass"], verdict["n_top"], verdict["n_match"]
    assert n_top == min(RULE["top_n"], n_pass)
    assert verdict["route"] in ("B2", "B1") and verdict["mzml_check"]
    assert all(verdict["mzml_check"].values()), "a verdict needs an mzML that passed the check"
    for k in ("weights_sha256", "genome_sha256"):
        assert HEX64.match(verdict[k] or ""), k
    few = n_pass < RULE["min_peptides"]
    if n_top == 0:
        assert verdict["share"] is None and verdict["fails"]["low_genome_match"] is None
        low = False
    else:
        assert abs(verdict["share"] - n_match / n_top) <= 5e-5
        lo, hi = wilson_ci(n_match, n_top)
        assert abs(verdict["share_ci95"][0] - lo) <= 5e-5
        assert abs(verdict["share_ci95"][1] - hi) <= 5e-5
        low = n_match / n_top < RULE["min_share"]
        assert verdict["fails"]["low_genome_match"] is low
    assert verdict["fails"]["too_few_peptides"] is few
    want = {(False, False): None, (True, False): "too_few_peptides",
            (False, True): "low_genome_match", (True, True): "both"}[(few, low)]  # fmt: skip
    assert verdict["reason"] == want
    assert verdict["verdict"] == ("pass" if want is None else "fail")
    assert verdict["evidence"] is None
    assert sum(verdict["length_hist"].values()) == n_pass


def test_provenance_note_exists_only_on_a_pass(verdict: dict[str, Any]) -> None:
    if verdict["verdict"] != "pass":
        assert not PROVENANCE.exists(), "PROVENANCE.md is committed on a pass only (deliverable 13)"
        return
    text = require(PROVENANCE, "deliverable 13").read_text()
    assert "0.0.7" in text and verdict["weights_sha256"] in text
    for col in SSL_COLUMNS:
        assert col in text, f"the .ssl header column {col!r} (from the Cascadia source)"
    assert re.search(r"\b[0-9a-f]{64}\b", text.replace(verdict["weights_sha256"], "")), (
        "the env lock's and the input archive's sha256"
    )
    assert "WT_1.ssl.gz" in text and "$GIAP_ROOT" in text, "where the sample sits (not committed)"


def test_readme_has_the_a3_de_novo_section() -> None:
    text = require(README, "deliverable 12").read_text()
    m = re.search(r"^## A3 de novo\b.*?(?=^## |\Z)", text, re.M | re.S)
    assert m, "validation/pilot/README.md needs an '## A3 de novo' section"
    sec = m.group(0)
    for needle in ("Cascadia 0.0.7", "Apache-2.0", "TIMSCONVERT", "TDF-SDK", "a3_assets.tsv",
                   "a3_feasibility.py evaluate", "denovo_pilot.sbatch", "verdict"):  # fmt: skip
        assert needle in sec, needle


def test_metrics_dir_holds_no_per_peptide_file() -> None:
    d = PILOT / "metrics"
    names = sorted(p.name for p in d.iterdir()) if d.is_dir() else []
    assert not [n for n in names if "top1000" in n or n.endswith((".tsv", ".ssl", ".gz"))], names


def test_recorded_cascadia_dir_is_absent_or_provenance_only() -> None:
    d: Path = PROVENANCE.parent.parent
    if d.exists():
        files = sorted(p.relative_to(d).as_posix() for p in d.rglob("*") if p.is_file())
        assert files in ([], ["0.0.7/PROVENANCE.md"]), files
