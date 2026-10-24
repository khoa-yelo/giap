"""G4-k3 companion: the exact-error ORFs carry ``denovo_only_unsupported``
(docs/specs/M7-17-denovo-support-flag.md acceptance criterion 6, lead decision 1; ADR-0050 option
3; M7-13 rule 8).

Needs the M7-13b harness ``helpers.m7_gate_novel_records`` and the novel half of
``_m7_gate`` (PR #271, ``m7/calibration-gate-novel``); both are imported inside the record so the
module collects without them.

One record per replicate (stream ``novel_orfs.g4k3.flag``): replicate i simulates
``novel_orf_power`` with the gate seed, injects the G4-k3 false ORFs with the G4-k3 stream
(``replicate_rng("novel_orfs.g4k3", i)``, so the injected reads are those of G4-k3 replicate i)
and runs GIAP through ``novel_orfs`` on the injected de novo file and on the database file.

* *Asserted, every replicate:* each injected false ORF (a row of ``novel_orfs.tsv`` at an
  ``inject_false_orfs`` key; its peptides are all de novo, so ``has_denovo and not db_support``)
  whose side row has ``q <= 0.05`` and ``class_rule_met`` carries the flag (side cell and
  ``flags``) unless its TSV class is ``cds_outframe``/``antisense`` with ``rbs_supported`` true;
  such an exempt row does **not** carry it (rule 4 both ways). ``intergenic`` rows are never
  exempt. A violation is an ``AssertionError``. Zero injected rows reaching ``q <= 0.05`` with
  the class rule in a class over the block, or the injection not reaching GIAP
  (``check_reach``), is a ``RuntimeError``.
* *Recorded* (``novel_orfs_gate_diagnostics.jsonl``): per class the reached rows, the flagged
  share and the escape share (exempt rows / reached rows). Lead decision 1: an escape share above
  0.2 in ``cds_outframe`` or ``antisense`` goes to Khoa; it is not asserted.
* *Recorded:* per regime and class, the share of true ORFs (truth key and class match) with
  ``q <= 0.01`` and the class rule that carry the flag: ``picked_decoy`` from the database run,
  ``denovo_subgroup`` from the injected de novo run (the injected reads touch non-coding segments
  only, so its true ORFs are those of the plain de novo arm).

Reading pinned by the test author: the record runs its own G4-k3 arm (the M7-13b ``k3_record``
keeps only counts, not the tables), so the fast and full tiers each simulate once more per
replicate; R is ``GATE_R`` (10 fast / 300 full).
"""

from __future__ import annotations

import json
import math
import shutil
from collections.abc import Mapping
from functools import cache
from pathlib import Path
from typing import Any

import _criteria as C
import pandas as pd
import pytest

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("novel_orfs")]
NIGHTLY = pytest.mark.calibration_nightly

FLAG = "denovo_only_unsupported"
STREAM = "novel_orfs.g4k3.flag"
SIX_CLASSES = ("intergenic", "cds_outframe", "antisense")
EXEMPT_CLASSES = frozenset({"cds_outframe", "antisense"})
GATED_CLASSES = ("intergenic", "cds_outframe", "antisense", "nterm_extension")
REGIMES = ("picked_decoy", "denovo_subgroup")
Q_INJECTED, Q_TRUE = 0.05, 0.01
ESCAPE_REPORT = 0.2  # lead decision 1: above this the lead takes the RBS exemption to Khoa


def _M() -> Any:
    from helpers import m7_gate_novel_records as M

    return M


def _G() -> Any:
    import _m7_gate as G

    return G


def _true(v: object) -> bool:
    if isinstance(v, str):
        return v.strip().lower() == "true"
    return False if v is None or pd.isna(v) else bool(v)  # type: ignore[call-overload]


def _flag_cells(run: Any) -> tuple[dict[str, bool | None], dict[str, bool]]:
    """``(side cell, FLAG in flags)`` by ``orf_id``; a side cell is None when the column is
    absent."""
    side = run.side
    cells: dict[str, bool | None] = {}
    for r in side.to_dict("records"):
        v = r.get(FLAG)
        cells[str(r["orf_id"])] = None if FLAG not in side.columns or pd.isna(v) else bool(v)
    tsv = {str(r["orf_id"]): FLAG in str(r["flags"]).split(";")
           for r in run.novel_orfs.to_dict("records")}  # fmt: skip
    return cells, tsv


def injected_rows(run: Any, keys: Mapping[str, list[list[Any]]],
                  lengths: Mapping[str, int]) -> list[dict[str, Any]]:  # fmt: skip
    """One dict per ``novel_orfs.tsv`` row at an injected key: the injected class, the row's
    class, ``q``, ``class_rule_met``, ``rbs_supported`` and both flag readings."""
    G = _G()
    want = {(str(c), str(s), int(p)): cls for cls, v in keys.items() for c, s, p in v}
    side = {str(r["orf_id"]): r for r in run.side.to_dict("records")}
    cells, tsv = _flag_cells(run)
    out = []
    for r in run.novel_orfs.to_dict("records"):
        c, s = str(r["contig"]), str(r["strand"])
        key = (c, s, G.novel_stop(s, int(r["end"]), lengths.get(c)))
        if key not in want:
            continue
        oid = str(r["orf_id"])
        sr = side[oid]
        out.append({"injected_class": want[key], "class": str(r["class"]), "orf_id": oid,
                    "q": float(sr["q"]), "rule": bool(sr["class_rule_met"]),
                    "rbs": _true(r["rbs_supported"]), "side_flag": cells[oid],
                    "tsv_flag": tsv[oid]})  # fmt: skip
    return out


def true_flag_counts(run: Any, truth: pd.DataFrame, lengths: Mapping[str, int],
                     regime: str) -> dict[str, dict[str, int]]:  # fmt: skip
    """Per class: true ORFs (truth key and class) at ``q <= 0.01`` with the class rule in
    ``regime``, and how many carry the flag in ``flags``."""
    G = _G()
    truth_keys = {(str(t["contig"]), str(t["strand"]),
                   int(t["stop"]) if lengths.get(str(t["contig"])) is None else
                   (int(t["stop"]) - 1) % lengths[str(t["contig"])] + 1): str(t["class"])
                  for t in truth.to_dict("records")}  # fmt: skip
    _, tsv = _flag_cells(run)
    by_id = {str(r["orf_id"]): r for r in run.novel_orfs.to_dict("records")}
    out = {c: {"n": 0, "flagged": 0} for c in GATED_CLASSES}
    s = run.side
    for r in s.to_dict("records"):
        cls = str(r["class"])
        if cls not in out or r["regime"] != regime or not bool(r["class_rule_met"]):
            continue
        q = float(r["q"])
        if not (math.isfinite(q) and q <= Q_TRUE):
            continue
        t = by_id[str(r["orf_id"])]
        c, st = str(t["contig"]), str(t["strand"])
        if truth_keys.get((c, st, G.novel_stop(st, int(t["end"]), lengths.get(c)))) != cls:
            continue
        out[cls]["n"] += 1
        out[cls]["flagged"] += int(tsv[str(r["orf_id"])])
    return out


@cache
def flag_record(i: int, size: str = "1M") -> dict[str, Any]:
    M = _M()
    root = M._root("flag")
    try:
        res = M.simulate(i, size, root)
        seed, lengths, truth = M.sim_seed(i), M._lengths(res), M._truth(res)
        path, info = M.inject_false_orfs(res, C.replicate_rng(M.K3_BASE, i))
        dn = M.run_giap(res, path, root / "k3" / "OUT", seed)
        db = M.run_giap(res, Path(res.outdir) / "peptides.tsv", root / "db" / "OUT", seed)
        rows = injected_rows(dn, info["keys"], lengths)
        rec = {"stream": STREAM, "replicate": i, "seed": seed, "size": size,
               "injected": {c: len(v) for c, v in info["keys"].items()},
               "entered": M.entered(dn, info["keys"], lengths), "rows": rows,
               "side_has_flag": FLAG in dn.side.columns and FLAG in db.side.columns,
               "true": {"picked_decoy": true_flag_counts(db, truth, lengths, "picked_decoy"),
                        "denovo_subgroup": true_flag_counts(dn, truth, lengths,
                                                            "denovo_subgroup")}}  # fmt: skip
    finally:
        shutil.rmtree(root, ignore_errors=True)
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "novel_orfs_gate_replicates.jsonl").open("a") as fh:
        fh.write(json.dumps(rec, default=float) + "\n")
    return rec


def _reached(r: Mapping[str, Any]) -> bool:
    return r["rule"] and math.isfinite(r["q"]) and r["q"] <= Q_INJECTED


def _exempt(r: Mapping[str, Any]) -> bool:
    return r["rbs"] and r["class"] in EXEMPT_CLASSES


def summarise(recs: list[Mapping[str, Any]]) -> dict[str, Any]:
    """Per injected class: reached, flagged, exempt, escape share; per regime/class the true-ORF
    flagged share; the violations (row dicts with the replicate)."""
    per = {c: {"reached": 0, "flagged": 0, "exempt": 0} for c in SIX_CLASSES}
    violations = []
    for rec in recs:
        for r in rec["rows"]:
            if not _reached(r) or r["injected_class"] not in per:
                continue
            p = per[r["injected_class"]]
            p["reached"] += 1
            p["exempt"] += int(_exempt(r))
            flagged = bool(r["side_flag"]) and bool(r["tsv_flag"])
            p["flagged"] += int(flagged)
            if flagged is _exempt(r) or r["side_flag"] is None:
                violations.append({**r, "replicate": rec["replicate"]})
    for p in per.values():
        n = p["reached"]
        p["flagged_share"] = p["flagged"] / n if n else math.nan
        p["escape_share"] = p["exempt"] / n if n else math.nan
    true: dict[str, dict[str, Any]] = {}
    for g in REGIMES:
        for c in GATED_CLASSES:
            n = sum(rec["true"][g][c]["n"] for rec in recs)
            f = sum(rec["true"][g][c]["flagged"] for rec in recs)
            true[f"{g}.{c}"] = {"n": n, "flagged": f, "share": f / n if n else math.nan}
    return {"per_class": per, "true": true, "violations": violations}


def _check(full: bool) -> None:
    M = _M()
    R = int(_G().GATE_R["full" if full else "fast"])
    recs = [flag_record(i) for i in C.replicate_range(R)]
    s = summarise(recs)
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    escape = {c: s["per_class"][c]["escape_share"] for c in sorted(EXEMPT_CLASSES)}
    reached = [r for rec in recs for r in rec["rows"] if _reached(r)]
    moved: dict[str, int] = {}  # injected -> TSV class of reached rows (stat note 7f8e4ed0)
    for r in reached:
        k = f"{r['injected_class']}->{r['class']}"
        moved[k] = moved.get(k, 0) + 1
    diag = {"check": "G4-k3 flag", "full": full, "R": len(recs), "per_class": s["per_class"],
            "true": s["true"], "n_violations": len(s["violations"]),
            "escape_above_report": {c: v > ESCAPE_REPORT for c, v in escape.items()},
            "injected_to_tsv_class": dict(sorted(moved.items()))}  # fmt: skip
    with (d / "novel_orfs_gate_diagnostics.jsonl").open("a") as fh:
        fh.write(json.dumps(diag, default=float) + "\n")
    injected = {c: sum(r["injected"].get(c, 0) for r in recs) for c in SIX_CLASSES}
    entered = sum(r["entered"]["rows"] for r in recs)
    M.check_reach(entered, sum(injected.values()), "G4-k3 flag")
    for c in SIX_CLASSES:
        if s["per_class"][c]["reached"] < 1:
            raise RuntimeError(f"G4-k3 flag vacuous: no injected {c} row at q <= {Q_INJECTED} "
                               f"with the class rule over {len(recs)} replicates")  # fmt: skip
    assert all(r["side_has_flag"] for r in recs), f"novel_class_fdr has no {FLAG} column"
    assert not s["violations"], s["violations"][:10]
    # Every reached row whose TSV class is intergenic carries the flag (rule 4 has no intergenic
    # exemption). Binned by TSV class, not the injected class: an injected segment can share its
    # stop with a planted antisense ORF and then *is* that ORF (§7.6.2 grouping, §7.2.6 priority);
    # statistician verdict on stat/m7-17-gate 7f8e4ed0. Rule 4 is still checked both ways on every
    # reached row above, with zero tolerance.
    tsv_intergenic = [r for r in reached if r["class"] == "intergenic"]
    assert tsv_intergenic, "G4-k3 flag vacuous: no reached row of TSV class intergenic"
    unflagged = [r for r in tsv_intergenic if not (r["side_flag"] and r["tsv_flag"])]
    assert not unflagged, unflagged[:10]


@NIGHTLY
def test_g4k3_exact_error_orfs_carry_the_flag() -> None:
    _check(False)


@pytest.mark.calibration_full
def test_g4k3_exact_error_orfs_carry_the_flag_full() -> None:
    _check(True)
