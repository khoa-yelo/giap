"""M7 gate, novel-ORF half (M7-13b): G3 novel-ORF F1, G4 class FDR per regime and class, G4-k3
clustered false ORFs (docs/specs/M7-13-calibration-gate.md § Checks, rules 1, 2, 5-8, 10-12,
risks 3-5, lead decisions 1-2; GIPPA_SPEC.md §7.6.4, §7.6.6, §8.4 row "Novel-ORF class FDR", §13
row M7; ADR-0015; M7-11b rules 6 and 11; M7-11b statistician review, PR #256 comment 6365221592).

R is the module constant ``_m7_gate.GATE_R`` (10 fast / 300 full; lead decision 1, ADR-0048
pending Khoa), the checks logged under suite ``novel_orfs`` (array task 11). One simulation
(``novel_orf_power(seed, "1M", n_per_class=NOVEL_PER_CLASS_GATE)``) per replicate feeds G3 and
both G4 regimes (``helpers.m7_gate_novel_records.gate_record``; seeds there). Fast checks are
nightly; the per-PR subset is the one-replicate smoke test at 100 kb.

* G3: pooled F1 over the four classes of the ``called = true`` rows (database input, regime
  ``picked_decoy``) >= ``f1_min.novel_orfs`` through ``calibrate_pooled`` (point estimate, Q4),
  guards: summed recall truth >= 100 per class, >= 50 for intergenic (rule 7). Recorded: F1 per
  class, the ``denovo_subgroup`` and ``unestimated`` (decoys stripped) arms' F1, start accuracy,
  the share of planted RBS starts passing the threshold, realised and false-candidate counts.
* G4: per regime (database -> ``picked_decoy``, ``peptides_denovo.tsv`` -> ``denovo_subgroup``),
  class and q in ``levels``, the ORFs with side-table ``q <= q`` and ``class_rule_met`` (NaN q,
  i.e. pick losers and the >= 2-peptide fallback calls, never enter) have mean FDP <= q + 2 SE
  (``calibrate_fdr_counts``). Vacuous (rule 6: block-1 mean list targets < 1/q) is logged and
  the test skips: neither a pass nor a fail. ``unestimated`` ORFs of a class for a reason the
  scenario rules out (scored makedb decoys and every de novo input are provided) fail the check
  (spec "Degenerate cases"). The fallback calls per class are recorded, not gated.
* G4-k3: G4's de novo regime with 20 injected false ORFs per six-frame class carrying 3-5 wrong
  de novo reads each (rule 8; ``inject_false_orfs``). Risk 3, confirmed by the full tier (R =
  300, Slurm 1875441 / 1876296 / 1877124): cds_outframe and antisense at q 0.01 and 0.05 and
  intergenic at 0.01 are strict xfails naming ADR-0050 (``raises=AssertionError``, lead decision
  10), so a fix shows as XPASS; the other cells stay hard checks.
* G4-k3 raises ``RuntimeError`` (so no xfail can absorb it) when the injected keys reaching
  GIPPA's table, summed over the replicates run, fall below ``REACH_MIN_SHARE`` (0.90) of those
  injected; when its ORFs are ``unestimated`` for an ``UNESTIMATED_BUGS`` reason (as G4); and
  when an ADR-0050 xfail cell turns vacuous (review of #271, B2 and nits 2-3).
* Lead decision 11: a ``no_makedb_pair`` ORF whose every peptide row names another planted
  ORF's makedb entry (a chance shared peptide) is exempt from ``UNESTIMATED_BUGS`` and recorded
  as ``shared_peptide_no_pair``; any other ``no_makedb_pair`` still fails.

Under ``--shard`` / ``GIPPA_CALIBRATION_CONFIRM`` a check has no verdict and the test skips after
writing its record (lead decision 6). No threshold is widened (R3).
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from typing import Any

import _criteria as C
import _m7_gate as G
import numpy as np
import pytest
import yaml
from helpers import m7_gate_novel_records as M

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("novel_orfs")]
NIGHTLY = pytest.mark.calibration_nightly  # fast tests > 60 s (rule 11)

SUITE = "novel_orfs"
G3_STREAM = "novel_orfs.g3"
NOT_FAILED = frozenset({"pass", "unconfirmed", C.SHARD_STATUS, C.SKIPPED_STATUS})
NO_VERDICT = frozenset({C.SHARD_STATUS, C.PENDING_STATUS, C.SKIPPED_STATUS})
G3_GUARD = {"n_truth_intergenic": 50, "n_truth_cds_outframe": 100, "n_truth_antisense": 100,
            "n_truth_nterm_extension": 100}  # fmt: skip
#: reasons an ORF of the four classes cannot have on this scenario ("Degenerate cases": the
#: scenario provides scored makedb-tagged decoys, and the step gets every de novo input)
UNESTIMATED_BUGS = frozenset({"paired_share_low", "decoys_unscored", "no_makedb_pair",
                              "decoys_not_makedb", "denovo_inputs_missing"})  # fmt: skip
LEVELS = tuple(C.load_tolerances().levels)


def gate_R(full: bool) -> int:
    return int(G.GATE_R["full" if full else "fast"])


def f1_min() -> float:
    return float(yaml.safe_load(C.TOLERANCES_PATH.read_text())["f1_min"]["novel_orfs"])


def diagnostics(record: Mapping[str, Any]) -> None:
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "novel_orfs_gate_diagnostics.jsonl").open("a") as fh:
        fh.write(json.dumps(dict(record), default=float) + "\n")


def ran_range(rec: Mapping[str, Any]) -> range:
    return range(0) if rec["status"] == C.SKIPPED_STATUS else range(rec["start"], rec["stop"])


def verdict(rec: Mapping[str, Any]) -> None:
    """Skip without a verdict (shard / confirm runs) or when vacuous (rule 6), else assert."""
    if rec["status"] in NO_VERDICT:
        pytest.skip("shard run: the verdict comes from --merge")
    if rec["status"] == C.VACUOUS_STATUS:
        pytest.skip(f"vacuous (rule 6: mean list targets < 1/q): {rec['stream']} {rec['key']}")
    assert rec["status"] in NOT_FAILED, rec


def _f1(c: Mapping[str, int]) -> float:
    return float(C.pooled_f1([c])[2])


def _sum(rows: list[Mapping[str, int]]) -> dict[str, int]:
    return {k: int(sum(r[k] for r in rows)) for k in rows[0]} if rows else {}


# --------------------------------------------------------------------------- per-PR smoke


def test_smoke_one_replicate_record_is_complete() -> None:
    rec = M.gate_record(0, "100k")
    for key in ("g3", "g3_denovo", "g3_unestimated"):
        assert set(G.NOVEL_COUNT_KEYS) <= set(rec[key]), (key, rec[key])
        assert all(isinstance(v, int) and v >= 0 for v in rec[key].values())
        assert rec[key]["tp_p"] <= rec[key]["n_pred"] and rec[key]["tp_r"] <= rec[key]["n_truth"]
    assert set(rec["g3_by_class"]) == set(M.NOVEL_CLASSES)
    want = {M.g4_key(r, c, q) for r in M.REGIMES for c in M.NOVEL_CLASSES for q in LEVELS}
    assert set(rec["g4"]) == want
    for v in rec["g4"].values():
        assert set(v) == set(M.FDR_KEYS) and 0 <= v["n_false"] <= v["n_called"]
    assert set(rec["attack_a"]) == set(M.NOVEL_CLASSES)
    assert set(rec["attack_d"]) == {M.g4_key("picked_decoy", c, q) for c in M.NOVEL_CLASSES
                                    for q in LEVELS}  # fmt: skip
    assert rec["g3"]["n_truth"] > 0
    assert sum(v["n_targets"] for k, v in rec["g4"].items() if k.startswith("picked")) > 0, (
        "no ORF entered a picked class list: the class FDR is not wired"
    )


# --------------------------------------------------------------------------- G3


def _g3(full: bool) -> None:
    R = gate_R(full)
    rec = C.calibrate_pooled(SUITE, G3_STREAM, "g3", f1_min(), M.gate_record, R,
                             guard=G3_GUARD)  # fmt: skip
    ran = [M.gate_record(i) for i in ran_range(rec)]
    if ran:
        by = {c: _sum([r["g3_by_class"][c] for r in ran]) for c in M.NOVEL_CLASSES}
        tp = sum(b["tp_p"] for b in by.values())
        diagnostics({
            "check": "G3", "full": full, "status": rec["status"], "f1": rec["observed"],
            "precision": rec["precision"], "recall": rec["recall"], "se": rec["se"],
            "f1_by_class": {c: _f1(b) for c, b in by.items()}, "by_class": by,
            "f1_denovo_subgroup": float(C.pooled_f1([r["g3_denovo"] for r in ran])[2]),
            "f1_unestimated": float(C.pooled_f1([r["g3_unestimated"] for r in ran])[2]),
            "start_accuracy": sum(b["start_ok"] for b in by.values()) / tp if tp else None,
            "rbs_pass_share": {c: b["rbs_pass"] / b["rbs_planted"] if b["rbs_planted"] else None
                               for c, b in by.items()},
            "mismatch": sum(b["mismatch"] for b in by.values()),
            "realised": _sum([r["realised"] for r in ran]),
            "false_candidates": {g: _sum([r["false_candidates"][g] for r in ran])
                                 for g in M.REGIMES},
            "seconds_per_replicate": float(np.mean([r["seconds"] for r in ran])),
        })  # fmt: skip
    verdict(rec)


@NIGHTLY
def test_g3_novel_orf_f1() -> None:
    _g3(False)


@pytest.mark.calibration_full
def test_g3_novel_orf_f1_full() -> None:
    _g3(True)


# --------------------------------------------------------------------------- G4


G4_CASES = [(r, c, q) for r in M.REGIMES for c in M.NOVEL_CLASSES for q in LEVELS]


def _g4(full: bool, regime: str, cls: str, q: float) -> None:
    R, key = gate_R(full), M.g4_key(regime, cls, q)
    rec = C.calibrate_fdr_counts(SUITE, f"novel_orfs.g4.{regime}.{cls}", key, q,
                                 lambda i: M.gate_record(i)["g4"], R,
                                 vacuous={"n_targets": 1.0 / q})  # fmt: skip
    ran = [M.gate_record(i) for i in ran_range(rec)]
    reasons: dict[str, int] = {}
    for r in ran:
        for why, n in r["diag"][regime]["reasons"][cls].items():
            reasons[why] = reasons.get(why, 0) + int(n)
    fallback = sum(int(r["diag"][regime]["fallback_calls"][cls]) for r in ran)
    shared = sum(int(r["diag"][regime].get(M.SHARED_NO_PAIR, {}).get(cls, 0)) for r in ran)
    vals = [r["g4"][key] for r in ran]
    diagnostics({"check": "G4", "full": full, "regime": regime, "class": cls, "q": q,
                 "status": rec["status"], "mean_fdp": rec["observed"], "bound": rec["bound"],
                 "mean_called": float(np.mean([v["n_called"] for v in vals])) if vals else None,
                 "mean_targets": float(np.mean([v["n_targets"] for v in vals]))
                 if vals else None,
                 "reasons": reasons, "fallback_calls": fallback,
                 M.SHARED_NO_PAIR: shared})  # fmt: skip
    bugs = M.unestimated_bugs(reasons, shared, UNESTIMATED_BUGS)  # lead decision 11
    assert not bugs, f"G4 {regime}/{cls}: unestimated ORFs on a scenario that rules them out {bugs}"
    verdict(rec)


@NIGHTLY
@pytest.mark.parametrize(("regime", "cls", "q"), G4_CASES)
def test_g4_class_fdr(regime: str, cls: str, q: float) -> None:
    _g4(False, regime, cls, q)


@pytest.mark.calibration_full
@pytest.mark.parametrize(("regime", "cls", "q"), G4_CASES)
def test_g4_class_fdr_full(regime: str, cls: str, q: float) -> None:
    _g4(True, regime, cls, q)


# --------------------------------------------------------------------------- G4-k3 (risk 3)


ADR_0050 = "ADR-0050 (proposed, pending Khoa; M7-13b lead decision 10)"
#: G4-k3 cells confirmed failing in the full tier (ADR-0050: exact systematic de novo errors and
#: the per-ORF minimum q); every other cell stays a hard check
K3_XFAIL = frozenset({("cds_outframe", 0.01), ("cds_outframe", 0.05), ("antisense", 0.01),
                      ("antisense", 0.05), ("intergenic", 0.01)})  # fmt: skip


def xfail_0050(why: str) -> pytest.MarkDecorator:
    """House strict xfail: only an AssertionError counts, so missing evidence still errors."""
    return pytest.mark.xfail(strict=True, raises=AssertionError, reason=f"{ADR_0050}: {why}")


K3_CASES = [
    pytest.param(c, q, marks=xfail_0050("G4-k3 FDP above q + 2 SE") if (c, q) in K3_XFAIL else ())
    for c in M.SIX_CLASSES
    for q in LEVELS
]


def _g4k3(full: bool, cls: str, q: float) -> None:
    R, key = gate_R(full), M.g4_key("denovo_subgroup", cls, q)
    first = M.k3_record(next(iter(C.replicate_range(R))))
    if first["injected"].get(cls, 0) < 1:  # not an AssertionError: a broken transform
        raise RuntimeError(f"G4-k3 vacuous: no {cls} false ORF injected {first['injected']}")
    rec = C.calibrate_fdr_counts(SUITE, f"novel_orfs.g4k3.{cls}", key, q,
                                 lambda i: M.k3_record(i)["g4k3"], R,
                                 vacuous={"n_targets": 1.0 / q})  # fmt: skip
    ran = [M.k3_record(i) for i in ran_range(rec)]
    injected, entered = _sum([r["injected"] for r in ran]), _sum([r["entered"] for r in ran])
    reasons: dict[str, int] = {}
    for r in ran:
        for why, n in r["diag"]["reasons"][cls].items():
            reasons[why] = reasons.get(why, 0) + int(n)
    shared = sum(int(r["diag"].get(M.SHARED_NO_PAIR, {}).get(cls, 0)) for r in ran)
    diagnostics({"check": "G4-k3", "full": full, "class": cls, "q": q, "status": rec["status"],
                 "mean_fdp": rec["observed"], "bound": rec["bound"], "injected": injected,
                 "entered": entered, "reasons": reasons,
                 M.SHARED_NO_PAIR: shared})  # fmt: skip
    # the checks below raise RuntimeError, never AssertionError: the ADR-0050 strict xfails
    # (raises=AssertionError) must not absorb a broken injection or a degenerate scenario
    if ran:
        M.check_reach(entered.get("rows", 0), sum(injected.values()), f"G4-k3 {cls}")
    bugs = M.unestimated_bugs(reasons, shared, UNESTIMATED_BUGS)  # as G4 (review nit 2)
    if bugs:
        raise RuntimeError(
            f"G4-k3 {cls}: unestimated ORFs on a scenario that rules them out {bugs}"
        )
    if (cls, q) in K3_XFAIL and rec["status"] == C.VACUOUS_STATUS:  # review nit 3
        raise RuntimeError(f"G4-k3 {cls} q {q}: an ADR-0050 xfail cell turned vacuous {rec}")
    verdict(rec)


@NIGHTLY
@pytest.mark.parametrize(("cls", "q"), K3_CASES)
def test_g4k3_clustered_false_orfs(cls: str, q: float) -> None:
    _g4k3(False, cls, q)


@pytest.mark.calibration_full
@pytest.mark.parametrize(("cls", "q"), K3_CASES)
def test_g4k3_clustered_false_orfs_full(cls: str, q: float) -> None:
    _g4k3(True, cls, q)
