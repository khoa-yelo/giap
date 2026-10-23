"""Gate M7 evidence ledger (GIAP_SPEC.md §13 row M7; docs/specs/M7-13c-gate-evidence.md rules 3-5,
lead decisions 1-2), modelled on tests/integration/test_m5_gate.py.

Reads the merged harness records of one full-tier run (``$GIAP_CALIBRATION_LOG_DIR``, rule 2:
``$GIAP_ROOT/logs/calibration/m7-gate-<C>``): ``variants.merged.jsonl`` (array task 13) and
``<NOVEL_SUITE>.merged.jsonl`` (task 11, ``novel_orfs``; ``GIAP_M7_NOVEL_SUITE=novel_orfs_gate``
for task 15 once ADR-0048 is accepted). Rule 3: a record counts only at the full R (variants 100,
novel gate 300) with its replicates covering ``0..R-1`` (and ``R..2R-1`` when block 2 decided) in
the raw ``<suite>.jsonl``, and a status other than ``shard``/``pending``/``skipped``. Missing or
short evidence raises ``RuntimeError``, never the strict xfail of a known failure (A6c, ADR-0049
decision 5, lead decision 2). G2 and A6a bypass the harness (M7-13a decision 7): their evidence is
the per-replicate ``variants_gate_replicates.jsonl`` plus a ``full`` diagnostics record. G4 cells
may be ``vacuous`` only where M7-13 rule 6 predicts it (picked_decoy intergenic q 0.01, or the de
novo regime). Attacks a/b/d/GC (fast tier, own log) and TE2 (pipeline run) are outside the ledger.
Skips unless ``GIAP_M7_GATE=1``.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests" / "calibration"))
sys.path.insert(0, str(REPO / "tests"))

FULL_R = {"variants": 100, "novel": 300}
LEVELS = (0.01, 0.05, 0.1)
CLASSES = ("intergenic", "cds_outframe", "antisense", "nterm_extension")
SIX_FRAME = CLASSES[:3]
REGIMES = ("picked_decoy", "denovo_subgroup")
PASSED = frozenset({"pass", "unconfirmed"})
NO_VERDICT = frozenset({"shard", "pending", "skipped"})
#: (suite, stream) -> the ADR accepting the failure; a strict xfail, so a fix XPASSes loudly
KNOWN_FAILURES = {("variants", "variants.a6c"): "ADR-0049 decision 5 (accepted known failure)"}
Check = tuple[str, str, str, str | None, float]  # (suite, stream, criterion, key, level)


def novel_suite() -> str:
    return os.environ.get("GIAP_M7_NOVEL_SUITE", "novel_orfs")


def _k(q: float) -> str:
    return f"{q:g}"


def required() -> dict[str, Check]:
    v, n = "variants", novel_suite()
    out: dict[str, Check] = {
        "G1": (v, "variants.g1", "pooled_f1", "g1", 0.90),
        "A6b": (v, "variants.a6b", "fwer", None, 0.05),
        "A6c": (v, "variants.a6c", "fwer", None, 0.05),
        "G3": (n, "novel_orfs.g3", "pooled_f1", "g3", 0.90),
    }
    for r in REGIMES:
        for c in CLASSES:
            for q in LEVELS:
                key = f"{r}.{c}.{_k(q)}"
                out[f"G4:{key}"] = (n, f"novel_orfs.g4.{r}.{c}", "fdr_counts", key, q)
    for c in SIX_FRAME:
        for q in LEVELS:
            key = f"denovo_subgroup.{c}.{_k(q)}"
            out[f"G4-k3:{c}.{_k(q)}"] = (n, f"novel_orfs.g4k3.{c}", "fdr_counts", key, q)
    return out


def predicted_vacuous(name: str) -> bool:
    """M7-13 rule 6: intergenic at q 0.01 (gap supply) or a de novo class; never G4-k3."""
    return name.startswith("G4:") and (
        name == "G4:picked_decoy.intergenic.0.01" or name.startswith("G4:denovo_subgroup.")
    )


def demanded() -> bool:
    return os.environ.get("GIAP_M7_GATE") == "1"


def log_dir() -> Path:
    return Path(os.environ.get("GIAP_CALIBRATION_LOG_DIR", REPO / "local" / "logs" / "calibration"))


def _read(name: str, d: Path | None = None) -> list[dict[str, Any]]:
    path = (d or log_dir()) / name
    if not path.is_file():
        raise RuntimeError(f"no evidence: {path} is missing (rule 3)")
    return [json.loads(x) for x in path.read_text(encoding="utf-8").splitlines() if x.strip()]


def _same(rec: dict[str, Any], check: Check, R: int) -> bool:
    _s, stream, crit, key, level = check
    return (rec.get("stream"), rec.get("criterion"), int(rec.get("R", -1))) == (stream, crit, R) \
        and abs(float(rec.get("level", -1)) - level) < 1e-9 \
        and (key is None or str(rec.get("key")) == key)  # fmt: skip


def evidence(check: Check, d: Path | None = None) -> dict[str, Any]:
    """The full-R merged verdict record of ``check`` (rule 3), or ``RuntimeError``."""
    suite = check[0]
    R = FULL_R["variants" if suite == "variants" else "novel"]
    merged = [r for r in _read(f"{suite}.merged.jsonl", d) if _same(r, check, R)]
    if len(merged) != 1:
        raise RuntimeError(f"{check}: {len(merged)} merged records at R = {R} (rule 3)")
    rec = merged[0]
    if rec["status"] in NO_VERDICT:
        raise RuntimeError(f"{check}: status {rec['status']} is no verdict (rule 3)")
    seen = {i for r in _read(f"{suite}.jsonl", d) if _same(r, check, R)
            for i in range(int(r["start"]), int(r["stop"]))}  # fmt: skip
    want = set(range(R * (2 if int(rec.get("block", 1)) == 2 else 1)))
    if want - seen:
        raise RuntimeError(f"{check}: replicates {min(want - seen)}.. not covered (rule 3)")
    return rec


def assert_passed(name: str, rec: dict[str, Any]) -> None:
    if rec["status"] == "vacuous" and predicted_vacuous(name):
        return
    assert rec["status"] in PASSED, f"{name}: {rec['status']} {rec.get('observed')} {rec}"


GATE = pytest.mark.skipif(not demanded(), reason="set GIAP_M7_GATE=1 to read the M7 gate run")
GATED = [n for n, c in required().items() if (c[0], c[1]) not in KNOWN_FAILURES]
KNOWN = [n for n, c in required().items() if (c[0], c[1]) in KNOWN_FAILURES]


@GATE
def test_ledger_matches_the_gate_suites() -> None:
    import _criteria as C
    import _m7_gate as G
    from helpers import m7_gate_novel_records as N
    from helpers import m7_gate_records as M

    assert int(G.GATE_R["full"]) == FULL_R["novel"]
    assert C.replicates("variants", True) == FULL_R["variants"]
    assert tuple(C.load_tolerances().levels) == LEVELS
    assert (C.POOLED_F1, C.FDR_COUNTS) == ("pooled_f1", "fdr_counts")
    assert (M.G1_BASE, M.A6C_BASE) == ("variants.g1", "variants.a6c")
    assert tuple(N.NOVEL_CLASSES) == CLASSES and tuple(N.SIX_CLASSES) == SIX_FRAME
    assert tuple(N.REGIMES) == REGIMES
    assert {N.g4_key(r, c, q) for r in REGIMES for c in CLASSES for q in LEVELS} == {
        c[3] for n, c in required().items() if n.startswith("G4:")
    }
    assert len(GATED) == 3 + 24 + 9 and KNOWN == ["A6c"]


@GATE
@pytest.mark.parametrize("name", GATED)
def test_gated_check_has_full_r_evidence_and_passed(name: str) -> None:
    assert_passed(name, evidence(required()[name]))


@GATE
@pytest.mark.parametrize(
    "name",
    [pytest.param(n, marks=pytest.mark.xfail(strict=True, raises=AssertionError,
                                             reason=KNOWN_FAILURES[required()[n][:2]]))
     for n in KNOWN],
)  # fmt: skip
def test_known_failure_has_full_r_evidence(name: str) -> None:
    rec = evidence(required()[name])  # missing evidence is a RuntimeError, never the xfail
    assert rec["status"] in PASSED, f"{name}: {rec['status']} (known failure)"


def _replicates(stream: str) -> dict[int, dict[str, Any]]:
    rows = [r for r in _read("variants_gate_replicates.jsonl") if r.get("stream") == stream]
    got = {int(r["replicate"]): r for r in rows}
    missing = sorted(set(range(FULL_R["variants"])) - set(got))
    check = {"variants.g2": "G2", "variants.a6": "A6a"}[stream]
    diag = _read("variants_gate_diagnostics.jsonl")
    if missing or not any(r.get("check") == check and r.get("full") is True for r in diag):
        raise RuntimeError(f"{stream}: replicates {missing[:5]} or the full diagnostics missing")
    return got


@GATE
def test_g2_same_strain_has_no_internal_high_row_in_any_replicate() -> None:
    bad = {i: r["violations"] for i, r in _replicates("variants.g2").items() if r["violations"]}
    assert not bad, f"G2 violations per replicate: {bad}"


@GATE
def test_a6a_every_isobaric_row_is_flagged_in_every_replicate() -> None:
    got = _replicates("variants.a6")
    if sum(int(r["a6a"]["n_rows"]) for r in got.values()) < 1:
        raise RuntimeError("A6a vacuous: no isobaric-only row (never evidence)")
    bad = {i: r["a6a"] for i, r in got.items() if r["a6a"]["n_unflagged"]}
    assert not bad, f"A6a unflagged isobaric rows: {bad}"


def test_the_evidence_rule_on_synthetic_records(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Switch, and Example 8 / rule 3: a fast record is ignored, a coverage gap is no evidence."""
    monkeypatch.delenv("GIAP_M7_GATE", raising=False)
    assert not demanded()
    monkeypatch.setenv("GIAP_M7_GATE", "1")
    assert demanded()
    check = required()["G3"]
    base = {"stream": "novel_orfs.g3", "criterion": "pooled_f1", "key": "g3", "level": 0.9}

    def write(name: str, recs: list[dict[str, Any]]) -> None:
        (tmp_path / name).write_text("".join(json.dumps({**base, **r}) + "\n" for r in recs))

    write("novel_orfs.merged.jsonl", [{"R": 10, "status": "pass", "start": 0, "stop": 10}])
    write("novel_orfs.jsonl", [{"R": 10, "start": 0, "stop": 10}])
    with pytest.raises(RuntimeError, match="0 merged records"):
        evidence(check, tmp_path)
    write("novel_orfs.merged.jsonl", [{"R": 300, "status": "pass", "block": 1}])
    write("novel_orfs.jsonl", [{"R": 300, "start": 0, "stop": 150}])
    with pytest.raises(RuntimeError, match=r"150\.\. not covered"):
        evidence(check, tmp_path)
    write("novel_orfs.jsonl", [{"R": 300, "start": 0, "stop": 150}, {"R": 300, "start": 150,
                                                                    "stop": 300}])  # fmt: skip
    assert evidence(check, tmp_path)["status"] == "pass"
    assert predicted_vacuous("G4:picked_decoy.intergenic.0.01")
    assert not predicted_vacuous("G4:picked_decoy.antisense.0.05")  # Example 4
    assert not predicted_vacuous("G4-k3:antisense.0.01")
