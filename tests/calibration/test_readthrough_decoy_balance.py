"""Decoy-balance check of the annotated-stop cells (M7-16a; spec test plan "Decoy-balance check";
statistician review change 2). Suite ``placement`` (``replicates_fast.placement`` /
``replicates_full.placement``); existing key ``se_multiplier`` only, no new tolerance (R3).

giapsim 1-Mb null genomes (no read-through planted) through ``run_chain`` (preflight and the
registered mapping steps), with the M3-19 database wrong-peptide arm (``reverse`` decoys,
``pi_wrong`` 0.05) and the M3-20 de novo arm (``shuffle`` decoys). Per replicate and decoy type:

- T_rt: best-cost target rows of the *false* target placements (no best-cost locus overlaps a
  true origin; helpers.placement_truth) whose window holds a marked target cell;
- D_rt: best-cost decoy rows of that decoy type whose window holds a mirrored cell;
- T: all best-cost target rows of the false placements; the read-through share T_rt / T.

Criterion, one-sided (fails only in the anti-conservative direction), the ADR-0019 estimate-check
form: mean(D_rt) >= mean(T_rt) - se_multiplier * SE(D_rt - T_rt). It is reported as
``exercised`` (mean T_rt >= 1 per replicate) or ``vacuous``; vacuous is an allowed outcome for
this check only, and expected at 1 Mb (T_rt ~ 0). The evidence for symmetry is the stressed arm
(test_readthrough_stressed.py) and the statistician's attack. All three numbers go to
``<log_dir>/placement_readthrough_balance.jsonl`` per decoy type.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from helpers.placement_sims import hits_of, run_chain
from helpers.placement_truth import best_loci
from helpers.readthrough_counts import best_cost_rows, holds_cell

pytestmark = [
    pytest.mark.calibration,
    pytest.mark.calsuite("placement"),
    pytest.mark.calibration_unshardable,
]

SUITE = "placement"
LENGTH = 1_000_000


def _counts(ch, kind: str) -> dict[str, float]:  # type: ignore[no-untyped-def]
    from giap.mapping.sixframe import load_stop_cells

    cells = load_stop_cells(ch.avail["sixframe"])
    hits = best_cost_rows(hits_of(ch.avail))
    hits = hits.assign(rt=holds_cell(hits, cells))
    per = pd.read_parquet(ch.avail["placement"]).drop_duplicates("peptide")
    per = per[per["decoy_type"] == kind]
    m_t = dict(zip(per["peptide"], per["m_target"].astype(int), strict=True))
    loci = best_loci(hits, m_t)
    false = {p for p in m_t if not any(ch.truth.overlaps(iv, p) for iv in loci.get(p, ()))}
    tgt = hits[(hits["search"] == "target") & hits["peptide"].isin(false)]
    tgt = tgt[tgt["peptide"].map(m_t).to_numpy() == tgt["n_diff"].to_numpy()]
    dec = hits[(hits["search"] == "decoy") & (hits["decoy_type"] == kind)]
    return {"T_rt": float(tgt["rt"].sum()), "D_rt": float(dec["rt"].sum()), "T": float(len(tgt)),
            "cells": float(sum(len(c) for c in cells.target))}  # fmt: skip


def _replicate(kind: str, seed: int, d: Path) -> dict[str, float]:
    from giapsim.proteomics import ProteomicsConfig

    pro = ProteomicsConfig(search_db="other", pi_wrong=0.05)  # type: ignore[arg-type]
    if kind == "reverse":
        ch = run_chain("database", seed, d, length=LENGTH, proteomics=pro)
    else:
        from giapsim.denovo import DenovoConfig

        ch = run_chain("denovo", seed, d, length=LENGTH, proteomics=pro, denovo=DenovoConfig())
    out = _counts(ch, kind)
    shutil.rmtree(d, ignore_errors=True)
    return out


def run_balance(kind: str, full: bool, tmp: Path) -> dict[str, object]:
    stream = f"placement.readthrough_balance.{kind}"
    R = C.replicates(SUITE, full)
    k = C.load_tolerances().se_multiplier
    reps = []
    for i in range(R):
        seed = int(C.replicate_rng(stream, i).integers(0, 2**31 - 1))
        reps.append(_replicate(kind, seed, tmp / f"r{i}"))
    t = np.array([r["T_rt"] for r in reps])
    dd = np.array([r["D_rt"] for r in reps])
    tot = np.array([r["T"] for r in reps])
    se = float((dd - t).std(ddof=1) / np.sqrt(len(t)))
    status = "exercised" if t.mean() >= 1.0 else "vacuous"
    rec = {
        "suite": stream, "R": R, "full": full, "status": status, "mean_T_rt": float(t.mean()),
        "mean_D_rt": float(dd.mean()), "mean_T": float(tot.mean()),
        "rt_share": float(t.sum() / tot.sum()) if tot.sum() else 0.0, "se": se,
        "mean_cells": float(np.mean([r["cells"] for r in reps])),
    }  # fmt: skip
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "placement_readthrough_balance.jsonl").open("a") as fh:
        fh.write(json.dumps(rec) + "\n")
    assert dd.mean() >= t.mean() - k * se, (
        f"{kind} ({status}): mean D_rt {dd.mean():.3f} < mean T_rt {t.mean():.3f} - {k} SE"
        f" ({se:.3f}), R = {R}"
    )
    return rec


@pytest.mark.parametrize("kind", ["reverse", "shuffle"])
def test_readthrough_decoy_balance(tmp_path: Path, kind: str) -> None:
    run_balance(kind, False, tmp_path)


@pytest.mark.calibration_full
@pytest.mark.parametrize("kind", ["reverse", "shuffle"])
def test_readthrough_decoy_balance_full(tmp_path: Path, kind: str) -> None:
    run_balance(kind, True, tmp_path)
