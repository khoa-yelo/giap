"""Genetic-code arms 4 and 25 of the placement-FDR and read-through decoy-balance suites
(M11-16f; docs/specs/M11-16-translation-tables-4-25.md item 9: "target-decoy balance with fewer
stops", "cells only at TAA/TAG"; ``calibration_full`` only, D14; the statistician reviews these
arms before 16f merges).

Arm 11 is ``test_placement_fdr_sim.py`` / ``test_readthrough_decoy_balance.py``, unchanged. An arm
here is the same replicate function and criterion with every gippasim genome at genetic code
``t`` and the run's ``annotation.translation_table = t`` (``helpers.genetic_code_arms.sim_at``).
The arm logs under its own stream (``<arm-11 stream>.code<t>``) and maps each replicate back to
arm 11's seed for the same index (``arm_seeds``): same seeds, same ``replicates_full.placement``,
same ``levels`` and ``se_multiplier``.

The read-through arm also checks that every annotated-stop cell of the code-``t`` index is at a
TAA/TAG codon (TGA is sense there), on every replicate.

The 1-Mb stressed stratum-0 estimate check
(``test_database_stressed_stratum0_estimate_covers_false_placements``, R = ``STRESS_EST_R`` = 100,
floor ``STRESS_EST_V0_FLOOR``) runs at code 4 too (statistician review of #481, N2), with arm 11's
seeds (stream ``placement.est_db_stressed``) and its record under ``.code4``.

Arms 11, 4 and 25 use common random numbers (the same replicate seeds), so they are highly
correlated and not independent confirmations of one another (review Q3).
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _placement_suite import run_fdr_suite, write_diagnostics
from helpers.genetic_code_arms import SENSE_TABLES, arm_seeds, sim_at
from helpers.placement_sims import run_chain, score_chain
from test_placement_fdr_sim import (
    STRESS_EST_R,
    STRESS_EST_V0_FLOOR,
    _class_diag,
    _collect_class_diag,
)

pytestmark = [
    pytest.mark.calibration_full,
    pytest.mark.calsuite("placement"),
    pytest.mark.calibration_unshardable,
]

LEVELS = C.load_tolerances().levels
LENGTH = 1_000_000


def _seeds(base: str, table: int) -> tuple[str, dict[int, int]]:
    arm = f"{base}.code{table}"
    return arm, arm_seeds(base, arm, 2 * C.replicates("placement", True))


def _proteomics(**kw: object) -> object:
    from gippasim.proteomics import ProteomicsConfig

    return ProteomicsConfig(search_db="other", **kw)  # type: ignore[arg-type]


@pytest.mark.parametrize("table", SENSE_TABLES)
@pytest.mark.parametrize("pi_wrong", [0.01, 0.05])
def test_database_shuffled_placement_fdr_code_arm_full(
    tmp_path: Path, pi_wrong: float, table: int
) -> None:
    arm, seed_of = _seeds(f"placement.sim_db_pi{pi_wrong}", table)

    def replicate(seed: int, d: Path):  # type: ignore[no-untyped-def]
        with sim_at(table):
            ch = run_chain("database", seed_of[seed], d, length=LENGTH,
                           proteomics=_proteomics(pi_wrong=pi_wrong), paralogs=20)  # fmt: skip
        assert ch.cfg.annotation.translation_table == table
        return score_chain(ch, LEVELS, d)

    run_fdr_suite(arm, replicate, tmp_path, (), full=True, mutant=True)


@pytest.mark.parametrize("table", SENSE_TABLES)
def test_database_stressed_stratum0_placement_fdr_code_arm_full(tmp_path: Path, table: int) -> None:
    arm, seed_of = _seeds("placement.sim_db_stressed", table)

    def replicate(seed: int, d: Path):  # type: ignore[no-untyped-def]
        with sim_at(table):
            ch = run_chain("database", seed_of[seed], d, length=LENGTH,
                           proteomics=_proteomics(pi_wrong=0.05, n_runs=1),
                           inject_short=400)  # fmt: skip
        assert ch.cfg.annotation.translation_table == table
        return score_chain(ch, LEVELS, d)

    run_fdr_suite(arm, replicate, tmp_path, [("reverse", 0)], full=True)


@pytest.mark.parametrize("table", (4,))
def test_database_stressed_stratum0_estimate_code_arm_full(tmp_path: Path, table: int) -> None:
    """Arm ``table`` of ``test_database_stressed_stratum0_estimate_covers_false_placements``: the
    same replicates (arm 11's seeds), floor and rule mean(W0) >= mean(V0) - k SE."""
    stream = "placement.est_db_stressed"  # arm 11's seeds
    w0, v0, cs = [], [], []
    by_class: dict[str, list[list[float]]] = {}
    for i in range(STRESS_EST_R):
        seed = int(C.replicate_rng(stream, i).integers(0, 2**31 - 1))
        d = tmp_path / f"r{i}"
        with sim_at(table):
            ch = run_chain("database", seed, d, length=LENGTH,
                           proteomics=_proteomics(pi_wrong=0.05, n_runs=1),
                           inject_short=400)  # fmt: skip
        assert ch.cfg.annotation.translation_table == table
        rep = score_chain(ch, LEVELS, d)
        w0.append(float(rep.extra["W0_reverse"]))
        v0.append(float(rep.extra["V0_reverse"]))
        _collect_class_diag(rep.extra, cs, by_class)
        shutil.rmtree(d, ignore_errors=True)
    w, v = np.array(w0), np.array(v0)
    k = C.load_tolerances().se_multiplier
    se = float((v - w).std(ddof=1) / np.sqrt(len(v)))
    rec = {"suite": f"{stream}.code{table}", "n": len(v), "se": se,
           **_class_diag(cs, by_class)}  # fmt: skip
    write_diagnostics({**rec, "mean_W0": float(w.mean()), "mean_V0": float(v.mean())})
    print(f"code {table} stressed estimate: {rec | {'W0': w.mean(), 'V0': v.mean()}}")
    assert v.sum() >= STRESS_EST_V0_FLOOR, f"stratum 0 under-exercised: total V0 {v.sum():.0f}"
    assert w.mean() >= v.mean() - k * se, (
        f"code {table}: estimated stratum-0 false placements {w.mean():.3f} < realised"
        f" {v.mean():.3f} - {k} SE ({se:.3f}), {len(v)} replicates"
    )


# ------------------------------------------------------------------ read-through decoy balance


def _counts(ch, kind: str, table: int) -> dict[str, float]:  # type: ignore[no-untyped-def]
    """``test_readthrough_decoy_balance._counts`` with the stop cells read at ``table``."""
    from helpers.placement_sims import hits_of
    from helpers.placement_truth import best_loci
    from helpers.readthrough_counts import best_cost_rows, holds_cell

    from gippa.io.fasta import read_genome_fasta
    from gippa.mapping.sixframe import load_stop_cells

    cells = load_stop_cells(ch.avail["sixframe"], table=table)
    genome = {r.id: r.sequence for r in read_genome_fasta(Path(ch.res.genome_fasta))}
    for s in cells.sites:
        n = len(genome[s.contig])
        nt = "".join(genome[s.contig][(s.codon_start + j) % n] for j in range(3))
        if s.strand == "-":
            nt = nt.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        assert nt in {"TAA", "TAG"}, (table, s, nt)
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
    return {
        "T_rt": float(tgt["rt"].sum()),
        "D_rt": float(dec["rt"].sum()),
        "T": float(len(tgt)),
        "cells": float(sum(len(c) for c in cells.target)),
        "sites": float(len(cells.sites)),
    }


def _replicate(kind: str, seed: int, d: Path, table: int) -> dict[str, float]:
    pro = _proteomics(pi_wrong=0.05)
    with sim_at(table):
        if kind == "reverse":
            ch = run_chain("database", seed, d, length=LENGTH, proteomics=pro)
        else:
            from gippasim.denovo import DenovoConfig

            ch = run_chain("denovo", seed, d, length=LENGTH, proteomics=pro, denovo=DenovoConfig())
    assert ch.cfg.annotation.translation_table == table
    out = _counts(ch, kind, table)
    shutil.rmtree(d, ignore_errors=True)
    return out


@pytest.mark.parametrize("table", SENSE_TABLES)
@pytest.mark.parametrize("kind", ["reverse", "shuffle"])
def test_readthrough_decoy_balance_code_arm_full(tmp_path: Path, kind: str, table: int) -> None:
    stream = f"placement.readthrough_balance.{kind}"  # arm 11's seeds
    R = C.replicates("placement", True)
    k = C.load_tolerances().se_multiplier
    reps = []
    for i in range(R):
        seed = int(C.replicate_rng(stream, i).integers(0, 2**31 - 1))
        reps.append(_replicate(kind, seed, tmp_path / f"r{i}", table))
    t = np.array([r["T_rt"] for r in reps])
    dd = np.array([r["D_rt"] for r in reps])
    tot = np.array([r["T"] for r in reps])
    se = float((dd - t).std(ddof=1) / np.sqrt(len(t)))
    status = "exercised" if t.mean() >= 1.0 else "vacuous"
    rec = {
        "suite": f"{stream}.code{table}", "R": R, "full": True, "status": status,
        "mean_T_rt": float(t.mean()), "mean_D_rt": float(dd.mean()), "mean_T": float(tot.mean()),
        "rt_share": float(t.sum() / tot.sum()) if tot.sum() else 0.0, "se": se,
        "mean_cells": float(np.mean([r["cells"] for r in reps])),
        "mean_sites": float(np.mean([r["sites"] for r in reps])),
    }  # fmt: skip
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "placement_readthrough_balance.jsonl").open("a") as fh:
        fh.write(json.dumps(rec) + "\n")
    assert min(r["sites"] for r in reps) > 0, "no annotated-stop site: the arm checks nothing"
    assert dd.mean() >= t.mean() - k * se, (
        f"code {table} {kind} ({status}): mean D_rt {dd.mean():.3f} < mean T_rt {t.mean():.3f}"
        f" - {k} SE ({se:.3f}), R = {R}"
    )
