"""Gate M7 evidence ledger (GIAP_SPEC.md §13 row M7; docs/specs/M7-13c-gate-evidence.md rules 2-5,
lead decisions 1-2, 9-10), modelled on tests/integration/test_m5_gate.py.

Reads the records of one full-tier run in ``$GIAP_CALIBRATION_LOG_DIR`` (required, rule 2:
``$GIAP_ROOT/logs/calibration/m7-gate-<C>``, never the shared fast-tier ``local/logs``):
``variants.merged.jsonl`` (array task 13) and ``<NOVEL_SUITE>.merged.jsonl`` (task 11,
``novel_orfs``; ``GIAP_M7_NOVEL_SUITE=novel_orfs_gate`` for task 15 once ADR-0048 is accepted).
Rule 3: a record counts only once per check at the full R (variants 100, novel gate 300), with its
replicates covering ``0..R-1`` (and ``R..2R-1`` when block 2 decided) in the raw ``<suite>.jsonl``
and a status other than ``shard``/``pending``/``skipped``. Missing or short evidence raises
``RuntimeError``, never the strict xfail of a known failure (A6c, ADR-0049 decision 5), which
xfails on status ``fail`` only. G2 and A6a bypass the harness (M7-13a decision 7): their evidence
is replicates ``0..99`` of ``variants_gate_replicates.jsonl`` plus a ``full`` diagnostics record.
G4 cells may be ``vacuous`` only where M7-13 rule 6 predicts it, and fail on any
``UNESTIMATED_BUGS`` reason in their full diagnostics (spec O3). Not checked here: attacks
a/b/d/GC, TE2, the reference arm, strand and origin facts (lead decision 9). The log-reading tests
skip unless ``GIAP_M7_GATE=1``; the synthetic-record tests always run.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml

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
#: test_novel_orfs_gate.UNESTIMATED_BUGS (pinned in the ledger test)
UNESTIMATED_BUGS = frozenset({"paired_share_low", "decoys_unscored", "no_makedb_pair",
                              "decoys_not_makedb", "denovo_inputs_missing"})  # fmt: skip
#: (suite, stream) -> the ADR accepting the failure; a strict xfail, so a fix XPASSes loudly
KNOWN_FAILURES = {("variants", "variants.a6c"): "ADR-0049 decision 5 (accepted known failure)"}
Check = tuple[str, str, str, str | None, float]  # (suite, stream, criterion, key, level)


def novel_suite() -> str:
    return os.environ.get("GIAP_M7_NOVEL_SUITE", "novel_orfs")


def f1_min() -> dict[str, float]:
    tol = yaml.safe_load((REPO / "tests" / "calibration" / "tolerances.yaml").read_text())
    return {k: float(v) for k, v in tol["f1_min"].items()}


def _k(q: float) -> str:
    return f"{q:g}"


def required() -> dict[str, Check]:
    v, n, f1 = "variants", novel_suite(), f1_min()
    out: dict[str, Check] = {
        "G1": (v, "variants.g1", "pooled_f1", "g1", f1["variants_high"]),
        "A6b": (v, "variants.a6b", "fwer", None, 0.05),
        "A6c": (v, "variants.a6c", "fwer", None, 0.05),
        "G3": (n, "novel_orfs.g3", "pooled_f1", "g3", f1["novel_orfs"]),
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
    raw = os.environ.get("GIAP_CALIBRATION_LOG_DIR")
    if not raw:
        raise RuntimeError("set GIAP_CALIBRATION_LOG_DIR to the gate run's own directory (rule 2)")
    return Path(raw)


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


def known_failure_verdict(name: str, rec: dict[str, Any]) -> None:
    """Pass (strict XPASS) on a passed status, AssertionError (the xfail) only on ``fail``."""
    if rec["status"] in PASSED:
        return
    if rec["status"] != "fail":
        pytest.fail(f"{name}: status {rec['status']} is not the known failure (never an xfail)")
    raise AssertionError(f"{name}: fail (known failure)")


def g4_bugs(name: str, d: Path | None = None) -> None:
    """Spec O3: an ``UNESTIMATED_BUGS`` reason in a G4 cell's full diagnostics is a failure."""
    _s, _stream, _c, key, q = required()[name]
    regime, cls = str(key).split(".")[:2]
    rows = [r for r in _read("novel_orfs_gate_diagnostics.jsonl", d)
            if r.get("check") == "G4" and r.get("full") is True and r.get("regime") == regime
            and r.get("class") == cls and abs(float(r.get("q", -1)) - q) < 1e-9]  # fmt: skip
    if not rows:
        raise RuntimeError(f"{name}: no full G4 diagnostics record")
    bugs = sorted({k for r in rows for k in r.get("reasons", {}) if k in UNESTIMATED_BUGS})
    assert not bugs, f"{name}: unestimated reasons the scenario rules out: {bugs}"


def _replicates(stream: str, d: Path | None = None) -> dict[int, dict[str, Any]]:
    """Replicates ``0..R-1`` of ``stream`` (later ones are other checks' block 2, n2)."""
    R = FULL_R["variants"]
    got = {int(r["replicate"]): r for r in _read("variants_gate_replicates.jsonl", d)
           if r.get("stream") == stream and int(r["replicate"]) < R}  # fmt: skip
    missing = sorted(set(range(R)) - set(got))
    check = {"variants.g2": "G2", "variants.a6": "A6a"}[stream]
    diag = _read("variants_gate_diagnostics.jsonl", d)
    if missing or not any(r.get("check") == check and r.get("full") is True for r in diag):
        raise RuntimeError(f"{stream}: replicates {missing[:5]} or the full diagnostics missing")
    return got


def g2_verdict(d: Path | None = None) -> None:
    bad = {i: r["violations"] for i, r in _replicates("variants.g2", d).items() if r["violations"]}
    assert not bad, f"G2 violations per replicate: {bad}"


def a6a_verdict(d: Path | None = None) -> None:
    got = _replicates("variants.a6", d)
    if sum(int(r["a6a"]["n_rows"]) for r in got.values()) < 1:
        raise RuntimeError("A6a vacuous: no isobaric-only row (never evidence)")
    bad = {i: r["a6a"] for i, r in got.items() if r["a6a"]["n_unflagged"]}
    assert not bad, f"A6a unflagged isobaric rows: {bad}"


GATE = pytest.mark.skipif(not demanded(), reason="set GIAP_M7_GATE=1 to read the M7 gate run")
GATED = [n for n, c in required().items() if (c[0], c[1]) not in KNOWN_FAILURES]
KNOWN = [n for n, c in required().items() if (c[0], c[1]) in KNOWN_FAILURES]
G4_CELLS = [n for n in GATED if n.startswith("G4:")]


def known_xfail(name: str) -> pytest.MarkDecorator:
    reason = KNOWN_FAILURES[required()[name][:2]]
    return pytest.mark.xfail(strict=True, raises=AssertionError, reason=reason)


def _gate_test_constants() -> dict[str, Any]:
    """Literal constants of test_novel_orfs_gate.py, read with ``ast`` (importing it outside
    tests/calibration meets its conftest-registered markers, test_m5_gate.py)."""
    import ast

    path = REPO / "tests" / "calibration" / "test_novel_orfs_gate.py"
    out: dict[str, Any] = {}
    for node in ast.parse(path.read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
            name, v = node.targets[0].id, node.value
            if name == "UNESTIMATED_BUGS" and isinstance(v, ast.Call) and v.args:
                out[name] = frozenset(ast.literal_eval(v.args[0]))
            elif isinstance(v, ast.Constant):
                out[name] = v.value
    return out


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
    t = _gate_test_constants()
    assert (t["SUITE"], t["G3_STREAM"]) == ("novel_orfs", "novel_orfs.g3")
    assert t["UNESTIMATED_BUGS"] == UNESTIMATED_BUGS
    assert tuple(N.NOVEL_CLASSES) == CLASSES and tuple(N.SIX_CLASSES) == SIX_FRAME
    assert tuple(N.REGIMES) == REGIMES and N.K3_BASE == "novel_orfs.g4k3"
    assert {N.g4_key(r, c, q) for r in REGIMES for c in CLASSES for q in LEVELS} == {
        c[3] for n, c in required().items() if n.startswith("G4:")
    }
    assert {N.g4_key("denovo_subgroup", c, q) for c in SIX_FRAME for q in LEVELS} == {
        c[3] for n, c in required().items() if n.startswith("G4-k3:")
    }
    assert len(GATED) == 3 + 24 + 9 and KNOWN == ["A6c"]


@GATE
@pytest.mark.parametrize("name", GATED)
def test_gated_check_has_full_r_evidence_and_passed(name: str) -> None:
    assert_passed(name, evidence(required()[name]))


@GATE
@pytest.mark.parametrize("name", G4_CELLS)
def test_g4_cell_has_no_unestimated_bug_reason(name: str) -> None:
    g4_bugs(name)


@GATE
@pytest.mark.parametrize("name", [pytest.param(n, marks=known_xfail(n)) for n in KNOWN])
def test_known_failure_has_full_r_evidence(name: str) -> None:
    known_failure_verdict(name, evidence(required()[name]))  # no evidence: RuntimeError


@GATE
def test_g2_same_strain_has_no_internal_high_row_in_any_replicate() -> None:
    g2_verdict()


@GATE
def test_a6a_every_isobaric_row_is_flagged_in_every_replicate() -> None:
    a6a_verdict()


# --------------------------------------------------------------------------- always on


def _write(d: Path, name: str, recs: list[dict[str, Any]], **base: Any) -> None:
    (d / name).write_text("".join(json.dumps({**base, **r}) + "\n" for r in recs))


G3_BASE = {"stream": "novel_orfs.g3", "criterion": "pooled_f1", "key": "g3", "level": 0.9}


def test_the_switch_and_the_log_dir_are_required(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("GIAP_M7_GATE", raising=False)
    assert not demanded()
    monkeypatch.setenv("GIAP_M7_GATE", "1")
    assert demanded()
    monkeypatch.delenv("GIAP_CALIBRATION_LOG_DIR", raising=False)
    with pytest.raises(RuntimeError, match="rule 2"):
        log_dir()


@pytest.mark.parametrize(
    ("merged", "raw", "match"),
    [
        ([{"R": 10, "status": "pass"}], [{"R": 10, "start": 0, "stop": 10}], "0 merged"),
        ([{"R": 300, "status": "pass"}] * 2, [{"R": 300, "start": 0, "stop": 300}], "2 merged"),
        ([{"R": 300, "status": "pending"}], [{"R": 300, "start": 0, "stop": 300}], "no verdict"),
        ([{"R": 300, "status": "pass", "block": 1}], [{"R": 300, "start": 0, "stop": 150}],
         r"150\.\. not covered"),
        ([{"R": 300, "status": "unconfirmed", "block": 2}], [{"R": 300, "start": 0, "stop": 300}],
         r"300\.\. not covered"),
    ],
)  # fmt: skip
def test_evidence_rule_rejects(
    tmp_path: Path, merged: list[dict[str, Any]], raw: list[dict[str, Any]], match: str
) -> None:
    _write(tmp_path, "novel_orfs.merged.jsonl", merged, **G3_BASE)
    _write(tmp_path, "novel_orfs.jsonl", raw, **G3_BASE)
    with pytest.raises(RuntimeError, match=match):
        evidence(required()["G3"], tmp_path)


def test_evidence_rule_accepts_shards_covering_both_blocks(tmp_path: Path) -> None:
    _write(tmp_path, "novel_orfs.merged.jsonl", [{"R": 300, "status": "unconfirmed", "block": 2}],
           **G3_BASE)  # fmt: skip
    _write(tmp_path, "novel_orfs.jsonl", [{"R": 300, "start": a, "stop": a + 150}
                                          for a in (0, 150, 300, 450)], **G3_BASE)  # fmt: skip
    assert evidence(required()["G3"], tmp_path)["status"] == "unconfirmed"


def test_vacuity_counts_only_where_rule_6_predicts_it() -> None:
    for name in ("G4:picked_decoy.intergenic.0.01", "G4:denovo_subgroup.antisense.0.05"):
        assert_passed(name, {"status": "vacuous"})
    for name in ("G4:picked_decoy.antisense.0.05", "G4-k3:antisense.0.01", "G1", "G3"):
        with pytest.raises(AssertionError):  # Example 4; decision 9
            assert_passed(name, {"status": "vacuous"})


def test_known_failure_xfails_only_on_fail() -> None:
    (param,) = test_known_failure_has_full_r_evidence.pytestmark[0].args[1]  # type: ignore[attr-defined]
    (mark,) = param.marks
    assert mark.name == "xfail" and mark.kwargs["strict"] is True
    assert mark.kwargs["raises"] is AssertionError and "ADR-0049" in mark.kwargs["reason"]
    known_failure_verdict("A6c", {"status": "pass"})  # strict XPASS
    with pytest.raises(AssertionError):
        known_failure_verdict("A6c", {"status": "fail"})
    with pytest.raises(pytest.fail.Exception) as exc:  # n1: vacuous is an error, not the xfail
        known_failure_verdict("A6c", {"status": "vacuous"})
    assert not isinstance(exc.value, AssertionError)


def test_g4_cell_fails_on_an_unestimated_bug_reason(tmp_path: Path) -> None:
    name, f = "G4:picked_decoy.antisense.0.05", "novel_orfs_gate_diagnostics.jsonl"
    with pytest.raises(RuntimeError, match="missing"):
        g4_bugs(name, tmp_path)
    base = {"check": "G4", "regime": "picked_decoy", "class": "antisense", "q": 0.05}
    _write(tmp_path, f, [{"full": False, "reasons": {}}], **base)
    with pytest.raises(RuntimeError, match="no full G4"):
        g4_bugs(name, tmp_path)
    _write(tmp_path, f, [{"full": True, "reasons": {"lost_pick": 3}}], **base)
    g4_bugs(name, tmp_path)
    _write(tmp_path, f, [{"full": True, "reasons": {}},
                         {"full": True, "reasons": {"decoys_unscored": 1}}], **base)  # fmt: skip
    with pytest.raises(AssertionError, match="decoys_unscored"):
        g4_bugs(name, tmp_path)


def _replicate_logs(d: Path, stream: str, rows: dict[int, dict[str, Any]], check: str) -> None:
    _write(d, "variants_gate_replicates.jsonl", [{"replicate": i, **r} for i, r in rows.items()],
           stream=stream)  # fmt: skip
    _write(d, "variants_gate_diagnostics.jsonl", [{"check": check, "full": True}])


def test_g2_needs_every_replicate(tmp_path: Path) -> None:
    _replicate_logs(tmp_path, "variants.g2", {i: {"violations": 0} for i in range(99)}, "G2")
    with pytest.raises(RuntimeError, match=r"\[99\]"):
        g2_verdict(tmp_path)
    _replicate_logs(tmp_path, "variants.g2", {i: {"violations": 0} for i in range(100)}, "G2")
    g2_verdict(tmp_path)


def test_a6a_empty_arm_is_an_error_and_block_2_rows_are_ignored(tmp_path: Path) -> None:
    def rows(n_rows: int) -> dict[int, dict[str, Any]]:
        out = {i: {"a6a": {"n_rows": n_rows, "n_unflagged": 0}} for i in range(100)}
        return {**out, 150: {"a6a": {"n_rows": 1, "n_unflagged": 1}}}  # n2: A6b/A6c block 2

    _replicate_logs(tmp_path, "variants.a6", rows(0), "A6a")
    with pytest.raises(RuntimeError, match="A6a vacuous"):
        a6a_verdict(tmp_path)
    _replicate_logs(tmp_path, "variants.a6", rows(1), "A6a")
    a6a_verdict(tmp_path)
