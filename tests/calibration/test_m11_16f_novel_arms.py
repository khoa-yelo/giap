"""Genetic-code arms 4 and 25 of the M7-13b novel-ORF gate (M11-16f; docs/specs/
M11-16-translation-tables-4-25.md item 9, D14; statistician review of #481, B2).

Arm 11 is ``test_novel_orfs_gate.py`` (G3, G4, G4-k3), unchanged. An arm here is G4 of that suite,
the novel-ORF class FDR per regime (database input -> ``picked_decoy``, ``peptides_denovo.tsv``
-> ``denovo_subgroup``), class and q in ``levels``, with the simulation at genetic code ``t``:
``dataclasses.replace(novel_orf_power(seed, "1M", n_per_class=NOVEL_PER_CLASS_GATE),
genetic_code=t)`` and ``annotation.translation_table: t`` in the run config
(``helpers.genetic_code_arms.stamp_config``; checked through ``resolve_config`` on every
replicate). Everything else is arm 11's:

* the criterion and path: ``calibrate_fdr_counts`` (mean FDP <= q + k SE, two blocks), vacuous
  when the block-1 mean list targets < 1/q (rule 6; logged, the test skips), the counts from
  ``helpers.m7_gate_novel_records.g4_counts``, the ``UNESTIMATED_BUGS`` check with lead decision
  11's shared-peptide exemption and its M11-16f amendment (``nonatg_nterm_recurrence``: a
  GTG / TTG gene's initiator-M N-terminal tryptic peptide recurring at a chance ATG elsewhere,
  recorded separately; diagnosis at code 25, replicates 14, 83, 179 and 213);
* R = ``_m7_gate.GATE_R["full"]`` (``calibration_full`` only, D14) and the seeds: replicate ``i``
  simulates with arm 11's seed for index ``i`` (``arm_seeds`` maps the arm stream's draw back to
  stream ``novel_orfs.gate``), so the arms are paired with arm 11;
* records under ``novel_orfs.<check>.code<t>`` (``novel_orfs.g4.<regime>.<class>.code<t>``; the
  replicate records go to ``novel_orfs_gate_replicates.jsonl`` with stream
  ``novel_orfs.gate.code<t>``).

ADR-0050: arm 11's strict xfail cells (``K3_XFAIL``) are G4-k3 cells only; arm 11's G4 has no
xfail cell, so neither has this arm: every G4 cell is a hard check, as in arm 11.

G3 (pooled novel-ORF F1 of the database run, ``called = true``) is recorded per arm, with the
per-class F1 and the planted ORFs holding an in-frame TGA (sense at 4 and 25: the ORFs a table-11
reading would cut), and not gated. The arm's premise is checked: some planted ORF holds an
in-frame TGA.

Arms 11, 4 and 25 use common random numbers (the same seeds per replicate index), so they are
highly correlated and not independent confirmations (review Q3).

Cost: two ``gippa run`` per replicate (database, de novo; arm 11 also runs the decoys-stripped
arm). The block's replicates are computed once per process in a forked pool
(``helpers.function_gate_chain.prefetch``; ``GIPPA_M815_WORKERS``, at most the process's CPU
share) and every check reads them; tests of one table share an ``xdist_group`` (run with
``--dist loadgroup`` so one worker owns a table).
"""

from __future__ import annotations

import dataclasses
import json
import shutil
import tempfile
import time
from functools import cache
from pathlib import Path
from typing import Any

import _criteria as C
import _m7_gate as G
import numpy as np
import pytest
from helpers import function_gate_chain as F
from helpers import m7_gate_novel_records as M
from helpers.genetic_code_arms import SENSE_TABLES, arm_seeds, in_frame_tga, stamp_config
from test_novel_orfs_gate import LEVELS, SUITE, UNESTIMATED_BUGS, diagnostics, ran_range, verdict

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("novel_orfs")]

SIZE = "1M"
_CACHE: dict[tuple[int, int, str], dict[str, Any]] = {}
_RC = str.maketrans("ACGTN", "TGCAN")


def gate_R() -> int:
    return int(G.GATE_R["full"])


def arm_stream(table: int) -> str:
    return f"{M.GATE_BASE}.code{table}"


@cache
def _seed_of(table: int) -> dict[int, int]:
    return arm_seeds(M.GATE_BASE, arm_stream(table), 2 * gate_R())


def sim_seed(table: int, i: int) -> int:
    """Arm 11's simulation seed of replicate ``i`` (``M.sim_seed(i)``), through ``arm_seeds``."""
    arm = int(C.replicate_rng(arm_stream(table), i).integers(0, 2**31 - 1))
    return _seed_of(table)[arm]


def _tga_orfs(res: Any) -> dict[str, int]:
    """Planted novel ORFs per class whose coding stretch holds an in-frame TGA."""
    genome, truth = M._genome(res), M._truth(res)
    out = dict.fromkeys(M.NOVEL_CLASSES, 0)
    for cls, c, s, e, st in zip(truth["class"], truth["contig"], truth["start"], truth["end"],
                                truth["strand"], strict=True):  # fmt: skip
        s, e = int(s), int(e)
        if s > e or str(cls) not in out:  # origin-crossing ORFs are not unwrapped
            continue
        nt = genome[str(c)][s - 1 : e]
        nt = nt if str(st) == "+" else nt.translate(_RC)[::-1]
        out[str(cls)] += int(in_frame_tga(nt) > 0)
    return out


def _record(key: tuple[int, int, str]) -> dict[str, Any]:
    """G3 and G4 counts of replicate ``i`` at genetic code ``table`` (pure in ``key``)."""
    from gippa.config import resolve_config
    from gippasim import scenarios as sc
    from gippasim.simulate import simulate

    table, i, size = key
    root, t0 = Path(tempfile.mkdtemp(prefix=f"m11-16f-novel-{table}-")), time.perf_counter()
    try:
        seed = sim_seed(table, i)
        n = sc.NOVEL_PER_CLASS_GATE if size == "1M" else None
        cfg = dataclasses.replace(sc.novel_orf_power(seed, size, n_per_class=n), genetic_code=table)
        res = simulate(cfg, seed, root / "sim")
        stamp_config(res.config_yaml, table)
        got = resolve_config(Path(res.config_yaml)).annotation.translation_table
        if got != table:  # not an AssertionError: a broken arm, never a calibration verdict
            raise RuntimeError(f"run config has translation_table {got}, want {table}")
        lengths, truth = M._lengths(res), M._truth(res)
        out = Path(res.outdir)
        db = M.run_gippa(res, out / "peptides.tsv", root / "db" / "OUT", seed)
        dn = M.run_gippa(res, out / "peptides_denovo.tsv", root / "dn" / "OUT", seed)
        g3 = M.g3_counts(db, truth, lengths)
        genes = M.gene_table(res.bakta_json)
        rec = {
            "stream": arm_stream(table), "replicate": i, "seed": seed, "size": size,
            "table": table, "g3": M._pool_counts(g3), "g3_by_class": g3,
            "g3_denovo": M._pool_counts(M.g3_counts(dn, truth, lengths)),
            "g4": {**M.g4_counts(db, truth, "picked_decoy", lengths),
                   **M.g4_counts(dn, truth, "denovo_subgroup", lengths)},
            "diag": {"picked_decoy": M._diag(db, M._tsv(out / "peptides.tsv"), truth, lengths,
                                             genes),
                     "denovo_subgroup": M._diag(dn, M._tsv(out / "peptides_denovo.tsv"), truth,
                                                lengths, genes)},
            "realised": json.loads((out / "truth" / "novel_orfs.json").read_text())["realised"],
            "tga_orfs": _tga_orfs(res),
            "seconds": time.perf_counter() - t0,
        }  # fmt: skip
    finally:
        shutil.rmtree(root, ignore_errors=True)
    M._log(rec)
    return rec


def record(table: int, i: int, size: str = SIZE) -> dict[str, Any]:
    """Cached ``_record``; a miss fills this process's block of replicates in a forked pool."""
    key = (table, i, size)
    if key not in _CACHE and size != SIZE:  # the smoke size: one replicate, no pool
        _CACHE[key] = _record(key)
    if key not in _CACHE:
        R = gate_R()
        idx = F.prefetch_indices(R)
        if i not in idx:  # block 2 of an unsharded run, or a confirm run
            idx = [] if F.sharded() else list(range((i // R) * R, (i // R + 1) * R))
        F.prefetch(_CACHE, _record, sorted({key, *((table, j, size) for j in idx)}))
    return _CACHE[key]


def _sum(rows: list[dict[str, int]]) -> dict[str, int]:
    return {k: int(sum(r[k] for r in rows)) for k in rows[0]} if rows else {}


# --------------------------------------------------------------------------- smoke


@pytest.mark.calibration_full
def test_smoke_code_record_is_complete() -> None:
    """One 100-kb replicate at code 4: the record is complete and the class lists are wired."""
    rec = record(4, 0, "100k")
    assert rec["table"] == 4 and rec["seed"] == M.sim_seed(0)
    want = {M.g4_key(r, c, q) for r in M.REGIMES for c in M.NOVEL_CLASSES for q in LEVELS}
    assert set(rec["g4"]) == want
    for v in rec["g4"].values():
        assert set(v) == set(M.FDR_KEYS) and 0 <= v["n_false"] <= v["n_called"]
    assert rec["g3"]["n_truth"] > 0
    assert sum(v["n_targets"] for k, v in rec["g4"].items() if k.startswith("picked")) > 0


# --------------------------------------------------------------------------- G3 (recorded)


@pytest.mark.calibration_full
@pytest.mark.parametrize(
    "table",
    [pytest.param(t, marks=pytest.mark.xdist_group(f"novel_code{t}")) for t in SENSE_TABLES],
)
def test_g3_novel_orf_f1_code_arm_recorded(table: int) -> None:
    """G3 at code ``table``: recorded, not gated (review B2). Skips under shard / confirm runs."""
    if F.sharded():
        pytest.skip("recorded-only arm: runs once per suite run, not per shard")
    ran = [record(table, i) for i in range(gate_R())]
    by = {c: _sum([r["g3_by_class"][c] for r in ran]) for c in M.NOVEL_CLASSES}
    p, r, f1 = C.pooled_f1([x["g3"] for x in ran])
    tga = _sum([x["tga_orfs"] for x in ran])
    diagnostics({
        "check": "G3", "full": True, "table": table, "stream": f"novel_orfs.g3.code{table}",
        "status": "recorded", "R": len(ran), "f1": f1, "precision": p, "recall": r,
        "f1_by_class": {c: float(C.pooled_f1([b])[2]) for c, b in by.items()}, "by_class": by,
        "f1_denovo_subgroup": float(C.pooled_f1([x["g3_denovo"] for x in ran])[2]),
        "tga_orfs": tga, "realised": _sum([x["realised"] for x in ran]),
        "seconds_per_replicate": float(np.mean([x["seconds"] for x in ran])),
    })  # fmt: skip
    print(f"G3 code {table}: F1 {f1:.4f} (P {p:.4f}, R {r:.4f}), TGA ORFs {tga}")
    assert sum(tga.values()) > 0, f"code {table}: no planted ORF holds an in-frame TGA"


# --------------------------------------------------------------------------- G4


G4_CASES = [
    pytest.param(t, r, c, q, marks=pytest.mark.xdist_group(f"novel_code{t}"))
    for t in SENSE_TABLES
    for r in M.REGIMES
    for c in M.NOVEL_CLASSES
    for q in LEVELS
]


@pytest.mark.calibration_full
@pytest.mark.parametrize(("table", "regime", "cls", "q"), G4_CASES)
def test_g4_class_fdr_code_arm_full(table: int, regime: str, cls: str, q: float) -> None:
    R, key = gate_R(), M.g4_key(regime, cls, q)
    rec = C.calibrate_fdr_counts(SUITE, f"novel_orfs.g4.{regime}.{cls}.code{table}", key, q,
                                 lambda i: record(table, i)["g4"], R,
                                 vacuous={"n_targets": 1.0 / q})  # fmt: skip
    ran = [record(table, i) for i in ran_range(rec)]
    reasons: dict[str, int] = {}
    for r in ran:
        for why, n in r["diag"][regime]["reasons"][cls].items():
            reasons[why] = reasons.get(why, 0) + int(n)
    fallback = sum(int(r["diag"][regime]["fallback_calls"][cls]) for r in ran)
    shared = sum(int(r["diag"][regime].get(M.SHARED_NO_PAIR, {}).get(cls, 0)) for r in ran)
    nonatg = sum(int(r["diag"][regime].get(M.NONATG_NTERM, {}).get(cls, 0)) for r in ran)
    vals = [r["g4"][key] for r in ran]
    diagnostics({"check": "G4", "full": True, "table": table, "regime": regime, "class": cls,
                 "q": q, "status": rec["status"], "mean_fdp": rec["observed"],
                 "bound": rec["bound"],
                 "mean_called": float(np.mean([v["n_called"] for v in vals])) if vals else None,
                 "mean_targets": float(np.mean([v["n_targets"] for v in vals]))
                 if vals else None,
                 "reasons": reasons, "fallback_calls": fallback,
                 M.SHARED_NO_PAIR: shared, M.NONATG_NTERM: nonatg})  # fmt: skip
    print(f"G4 code {table} {regime}/{cls} q {q}: {rec['status']} FDP {rec['observed']}"
          f" bound {rec['bound']} {M.SHARED_NO_PAIR} {shared} {M.NONATG_NTERM} {nonatg}"
          f" reasons {reasons}")  # fmt: skip
    # lead decision 11 and its M11-16f amendment (disjoint records)
    bugs = M.unestimated_bugs(reasons, shared + nonatg, UNESTIMATED_BUGS)
    assert not bugs, f"G4 {regime}/{cls}: unestimated ORFs on a scenario that rules them out {bugs}"
    verdict(rec)
