"""Teeth of the M7-13a variants gate: each check has a mutant that must fail
(docs/specs/M7-13-calibration-gate.md § Attacks: "Each attack has a mutant that must fail";
attack 6 row "A6c doubles as the teeth of A6b"; GIAP_SPEC.md §8.5 "A failing attack blocks the
merge"). The mutants reuse the gate's own replicates (``helpers.m7_gate_records``: one GIAP run
feeds the check and its mutant, rule 1) and are evaluated through the same criteria, logged to a
private directory (``teeth_log``), never the nightly suite log.

* G1 with the confidence filter removed (every substitution row predicted): a confirmed ``fail``
  of the pooled F1 criterion.
* G2's check on the isobaric attack's de novo reads (error-laden, no strain): violations in at
  least one replicate.
* A6a with ``mass_equivalent`` stripped from every row: unflagged rows in at least one replicate.
* A6b with every swap row taken as ``high``: a confirmed FWER ``fail``.

M7-13b (novel-ORF gate, suite ``novel_orfs``, ``helpers.m7_gate_novel_records``; R = the fast
``_m7_gate.GATE_R``; these tests carry their own ``calsuite`` marker):

* a: G4's picked regime with the class-FDR filter removed (every candidate meeting the class
  rule called): a confirmed FDR ``fail`` at q = 0.10 in at least one class.
* a and d run on the clustered-null arm (``cluster_record``: 400 null ``antisense`` makedb
  entries with 3-5 products each, lead decision 7; statistician Slurm 1874863, R = 200: a
  0.356 / 0.364 against 0.104, d 0.131 / 0.136 against 0.105 at q 0.10). On the plain gate
  replicates nearly every false ORF has one peptide and fails the class rule, so neither has
  teeth there. The arm's filtered G4 is recorded in ``novel_orfs_gate_diagnostics.jsonl``, not
  gated (statistician: at level, 0.004 / 0.037 / 0.078 against 0.011 / 0.053 / 0.105).
* b: the matched-TIR arm (``planted_tir = MATCHED_TIR``) at ``rbs_min_quantile`` 0.10 (G3's own
  database run) against a mutant run at 0.50: intergenic + ``nterm_extension`` recall (recorded)
  drops by more than 2 SE (paired over replicates).
* d: one picked list over the six-frame classes (``pooled_side``): a confirmed FDR ``fail`` in at
  least one (class, q); the per-class ``nterm_extension`` FDP is G4's own check.
* GC: false de novo peptides from GC-poor intergenic segments; intergenic FDP <= q + 2 SE in
  the de novo regime (a gate, logged to the suite log; no mutant). Fast tier: every level a hard
  check. Full tier (R = ``GATE_R["full"]``): q 0.01 only, a strict xfail naming ADR-0050 (lead
  decisions 10 and 13).
* 5i (reversed stretches in intergenic space) is covered by G4-k3.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import _criteria as C
import _m7_gate as G
import numpy as np
import pytest
import yaml
from helpers import m7_gate_novel_records as N
from helpers import m7_gate_records as M

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("variants")]
NIGHTLY = pytest.mark.calibration_nightly
SUITE = "variants"


def test_smoke_one_attack_replicate_record_is_complete() -> None:
    rec = M.a6_record(0)
    assert set(rec["a6a"]) == set(rec["a6a_mutant"]) == {"n_rows", "n_unflagged"}
    assert set(rec["a6b"]) == set(rec["a6b_mutant"]) == {"folded", "unfolded"}
    assert rec["a6a_mutant"]["n_rows"] == rec["a6a"]["n_rows"]
    assert rec["a6b_mutant"]["folded"] >= rec["a6b"]["folded"]
    assert isinstance(rec["g2_mutant_violations"], int)


@NIGHTLY
def test_g1_mutant_without_confidence_filter_fails(teeth_log: Path) -> None:
    R = C.replicates(SUITE, False)
    bound = float(yaml.safe_load(C.TOLERANCES_PATH.read_text())["f1_min"]["variants_high"])
    rec = C.calibrate_pooled(SUITE, M.G1_BASE, "mutant_any_confidence", bound,
                             lambda i: M.g1_record(i, "100k"), R, log_dir=teeth_log,
                             guard={"n_truth": 50, "n_pred": 50})  # fmt: skip
    assert rec["status"] == "fail", rec


def _a6(R: int) -> list[dict[str, Any]]:
    return [M.a6_record(i) for i in range(R)]


@NIGHTLY
def test_g2_check_finds_internal_high_rows_on_erroneous_reads() -> None:
    got = _a6(C.replicates(SUITE, False))
    assert any(r["g2_mutant_violations"] > 0 for r in got), "G2's check has no teeth"


@NIGHTLY
def test_a6a_mutant_without_mass_equivalent_fails() -> None:
    got = _a6(C.replicates(SUITE, False))
    assert any(r["a6a_mutant"]["n_unflagged"] > 0 for r in got), "A6a's check has no teeth"


@NIGHTLY
def test_a6b_mutant_every_swap_high_fails(teeth_log: Path) -> None:
    R, stream = C.replicates(SUITE, False), "variants.a6b.mutant"
    index = {C.suite_seed(stream, i): i for i in range(2 * R)}

    def fn(rng: Any) -> float:
        i = index[int(rng.bit_generator.seed_seq.entropy)]
        return float(M.a6_record(i)["a6b_mutant"]["folded"] > 0)

    rec = C.calibrate(SUITE, "fwer", 0.05, fn, R, stream=stream, log_dir=teeth_log)
    assert rec["status"] == "fail", rec


# --------------------------------------------------------------------------- M7-13b (novel ORFs)

NOVEL = pytest.mark.calsuite("novel_orfs")
NSUITE = "novel_orfs"
NLEVELS = tuple(C.load_tolerances().levels)


def _nR() -> int:
    return int(G.GATE_R["fast"])


@NOVEL
def test_smoke_one_novel_attack_replicate_record_is_complete() -> None:
    rec = N.gc_record(0, "100k")
    assert set(rec["gc"]) == {N.g4_key("denovo_subgroup", "intergenic", q) for q in NLEVELS}
    assert rec["injected"] >= 1 and rec["n_peptides"] >= rec["injected"]
    assert set(rec["entered"]) == {"rows", "finite_q", "called"}


@NOVEL
@NIGHTLY
def test_attack_a_without_class_fdr_filter_fails(teeth_log: Path) -> None:
    got = {c: C.calibrate_fdr_counts(NSUITE, "novel_orfs.attack_a", c, 0.10,
                                     lambda i: N.cluster_record(i)["attack_a"], _nR(),
                                     log_dir=teeth_log)["status"]
           for c in N.NOVEL_CLASSES}  # fmt: skip
    assert "fail" in got.values(), f"attack a has no teeth: {got}"


@NOVEL
@NIGHTLY
def test_attack_b_rbs_threshold_mutant_loses_recall() -> None:
    classes = ("intergenic", "nterm_extension")
    R = _nR()
    base = [N.gate_record(i)["g3_by_class"] for i in range(R)]
    mut = [N.rbs_record(i)["g3_by_class"] for i in range(R)]
    n = [sum(b[c]["n_truth"] for c in classes) for b in base]
    assert n == [sum(m[c]["n_truth"] for c in classes) for m in mut], "arms are not paired"
    tb = np.array([sum(b[c]["tp_r"] for c in classes) for b in base], dtype=float)
    tm = np.array([sum(m[c]["tp_r"] for c in classes) for m in mut], dtype=float)
    total = float(sum(n))
    assert total >= 50, f"attack b vacuous: {total} recall-truth rows"
    d = tb - tm
    drop = float(d.sum()) / total
    se = float(d.std(ddof=1)) * R**0.5 / total
    bound = C.load_tolerances().se_multiplier * se
    print(
        json.dumps(
            {
                "check": "attack_b",
                "recall_q010": float(tb.sum()) / total,
                "recall_q050": float(tm.sum()) / total,
                "drop": drop,
                "se": se,
            }
        )
    )
    assert drop > bound, f"attack b: recall drop {drop:.4g} <= 2 SE {bound:.4g}"


@NOVEL
@NIGHTLY
def test_attack_d_pooled_list_fails_in_a_class(teeth_log: Path) -> None:
    got = {}
    for c in N.NOVEL_CLASSES:
        for q in NLEVELS:
            k = N.g4_key("picked_decoy", c, q)
            got[k] = C.calibrate_fdr_counts(NSUITE, "novel_orfs.attack_d", k, q,
                                            lambda i: N.cluster_record(i)["attack_d"], _nR(),
                                            vacuous={"n_targets": 1.0 / q},
                                            log_dir=teeth_log)["status"]  # fmt: skip
    assert "fail" in got.values(), f"attack d has no teeth: {got}"


def _diagnostics(record: dict[str, Any]) -> None:
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "novel_orfs_gate_diagnostics.jsonl").open("a") as fh:
        fh.write(json.dumps(record, default=float) + "\n")


@NOVEL
def test_smoke_cluster_null_entries_avoid_planted_orfs_and_sim_hits() -> None:
    """The clustered-null transform at 100 kb: eligible entries are null by construction (no
    planted key, no overlap with a planted ORF on either strand, no simulated hit) and each has
    3-5 distinct products of its own segment; decoys are reversed under ``rev_<id>``."""
    import shutil
    import tempfile

    root = Path(tempfile.mkdtemp(prefix="m7gate-novel-cltest-"))
    try:
        res = N.simulate(0, "100k", root)
        pep = N._tsv(Path(res.outdir) / "peptides.tsv")
        seen = {p.removeprefix("rev_") for p in pep["protein_ids"] if "GIAPORF_" in p}
        truth = N._truth(res)
        keys = set(zip(truth["contig"].astype(str), truth["strand"].astype(str),
                       truth["stop"].astype(int), strict=True))  # fmt: skip
        segs = {N.makedb_id(N.CL_CLASS, s): s for s in N._segments(res, (N.CL_CLASS,))[N.CL_CLASS]}
        path, info = N.inject_null_entries(res, C.replicate_rng(N.CL_BASE, 0), m=60)
        out = N._tsv(path)
    finally:
        shutil.rmtree(root, ignore_errors=True)
    entries = info["entries"]
    assert len(entries) == min(60, info["n_eligible"]) and entries, info["n_eligible"]
    assert len({e["id"] for e in entries}) == len(entries)
    assert seen, "the simulation hit no GIAPORF_ entry"
    for e in entries:
        assert e["id"].startswith("GIAPORF_antisense_") and e["id"] not in seen, e["id"]
        assert tuple(e["key"]) not in keys, e
        lo, hi = e["span"]
        for c, a, b in zip(truth["contig"], truth["start"], truth["end"], strict=True):
            s0, e0 = sorted((int(a), int(b)))
            assert c != e["contig"] or not (s0 - 1 < hi and lo < e0), (e["id"], a, b)
        seg = segs[e["id"]]
        assert N.CL_PEPTIDES[0] <= len(e["products"]) <= N.CL_PEPTIDES[1]
        assert len(set(e["products"])) == len(e["products"])
        assert set(e["products"]) <= set(seg["products"])
        prot = ("rev_" if e["decoy"] else "") + e["id"]
        rows = out[out["protein_ids"] == prot]
        want = sorted(p[::-1] if e["decoy"] else p for p in e["products"])
        assert sorted(rows["peptide"]) == want, e["id"]
        assert set(rows["is_decoy"]) == {"true" if e["decoy"] else "false"}
    assert len(out) == len(pep) + info["n_rows"]
    assert out.iloc[: len(pep)].equals(pep), "the simulated rows must be kept unchanged"
    null = pep[pep["protein_ids"].str.startswith("rev_GIAPORF_")]
    trip = set(zip(null["q_value"], null["pep"], null["score"], strict=True))
    add = out.iloc[len(pep) :]
    assert set(zip(add["q_value"], add["pep"], add["score"], strict=True)) <= trip


@NOVEL
@NIGHTLY
def test_cluster_arm_filtered_g4_recorded(teeth_log: Path) -> None:
    """The clustered-null arm's own G4 (picked regime, q filter on): recorded, not gated (lead
    decision 7); the arm must give a non-vacuous antisense list at q = 0.10."""
    got = {}
    for c in N.NOVEL_CLASSES:
        for q in NLEVELS:
            k = N.g4_key("picked_decoy", c, q)
            rec = C.calibrate_fdr_counts(NSUITE, "novel_orfs.cluster.g4", k, q,
                                         lambda i: N.cluster_record(i)["g4"], _nR(),
                                         vacuous={"n_targets": 1.0 / q},
                                         log_dir=teeth_log)  # fmt: skip
            ran = [N.cluster_record(i) for i in range(rec["stop"])]  # every block run
            vals = [r["g4"][k] for r in ran]
            got[k] = rec["status"]
            _diagnostics({"check": "cluster_g4", "regime": "picked_decoy", "class": c, "q": q,
                          "status": rec["status"], "mean_fdp": rec["observed"],
                          "bound": rec["bound"], "gated": False, "replicates": len(ran),
                          "mean_called": float(np.mean([v["n_called"] for v in vals])),
                          "mean_targets": float(np.mean([v["n_targets"] for v in vals])),
                          "n_target_entries": float(np.mean([r["n_target_entries"]
                                                             for r in ran])),
                          "entered": {e: int(sum(r["entered"][e] for r in ran))
                                      for e in ("rows", "finite_q", "called")}})  # fmt: skip
    assert got[N.g4_key("picked_decoy", N.CL_CLASS, 0.1)] != C.VACUOUS_STATUS, got


def _gc(R: int, q: float) -> None:
    first = N.gc_record(next(iter(C.replicate_range(R))))
    if first["injected"] < 1:  # not an AssertionError: a broken transform
        raise RuntimeError(f"attack GC vacuous: nothing injected {first}")
    k = N.g4_key("denovo_subgroup", "intergenic", q)
    rec = C.calibrate_fdr_counts(NSUITE, "novel_orfs.gc", k, q, lambda i: N.gc_record(i)["gc"],
                                 R, vacuous={"n_targets": 1.0 / q})  # fmt: skip
    if rec["status"] in (C.SHARD_STATUS, C.PENDING_STATUS, C.SKIPPED_STATUS):
        pytest.skip("shard run: the verdict comes from --merge")
    if rec["status"] == C.VACUOUS_STATUS:
        pytest.skip(f"vacuous (rule 6): {k}")
    assert rec["status"] in ("pass", "unconfirmed"), rec


@NOVEL
@NIGHTLY
@pytest.mark.parametrize("q", NLEVELS)
def test_attack_gc_poor_intergenic_fdp_at_level(q: float) -> None:
    """Fast tier, every level a hard check (lead decision 13: its fixed seeds pass at 0.01)."""
    _gc(_nR(), q)


#: lead decision 13 (amending 10): GC at q 0.01 fails in the full tier (R = 300: 0.0122 / 0.0126
#: against 0.0112, ADR-0050); a strict xfail there, so a fix shows as XPASS
GC_FULL_XFAIL = pytest.mark.xfail(
    strict=True, raises=AssertionError,
    reason="ADR-0050 (proposed, pending Khoa; M7-13b lead decisions 10 and 13): exact GC-poor "
           "intergenic de novo reads above q + 2 SE")  # fmt: skip


@NOVEL
@pytest.mark.calibration_full
@GC_FULL_XFAIL
def test_attack_gc_poor_intergenic_fdp_at_level_full() -> None:
    _gc(int(G.GATE_R["full"]), 0.01)
