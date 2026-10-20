"""Genome-level placement FDR on giapsim through the mapping steps (issue M3-21; GIAP_SPEC.md
§7.2.9 "Calibration: placement FDR <= q + Monte Carlo error on simulated wrong peptides (database
and de novo error models)", §8.4 "Placement FDR" row; docs/specs/M3-21 lead decisions 1-12;
ADR-0019).

Every replicate: a giapsim genome (100 kb fast / 1 Mb full), M3-19 proteomics, ``giap
preflight`` and the registered mapping steps (helpers.placement_sims). FDP per list (decoy type x
m) and pooled, false placement = no accepted locus overlaps a true origin (helpers.placement_truth),
criterion mean FDP <= q + 2 SE at every level (``_placement_suite.run_fdr_suite``).

"Exercised" (decision 4, PR #113 review B2) = mean >= 1 per replicate of the list's peptides (target
hit at cost m, any q) that are *false*: no best-cost locus overlaps a true origin. A wrong string
placed at its true origin (de novo read, other-strain peptide) is not a null candidate.

**What the genome level shows.** At 100 kb no suite below except the stressed arm exercises any
list (measured 2026-10-21): tryptic >= 7-aa peptides that pass the q filter almost never hit a
100-kb-1-Mb six-frame translation by chance, and tolerant search needs >= 10 aa (``max_diff``
10-11 -> 1, 12-40 -> 2), so chance non-origin hits at m >= 1 are ~0 at <= 1 Mb. The unclaimed gates
are regression smoke checks (coordinates, oracle, wiring), not calibration evidence. The FDR claim
for **m >= 1 rests on the list-level suites** (test_placement_fdr_lists.py,
attacks/test_attack_placement_pi.py); stratum 0 at genome level on the stressed arm and the
estimate checks.

* ``database_shuffled`` (pi_wrong 0.01 / 0.05, 20 planted near-paralogs for the unique guard and
  its mutant): smoke check; the stratum-0 π̂ is checked directly
  (``test_database_stratum0_pi_hat_covers_the_wrong_share``).
* ``database_stressed`` (review B4, b2sim ``--mode short``): M3-19 table plus injected shuffled
  7-8-mer wrong database peptides, PEP = the true wrong share (delta = 0) -> claims reverse m = 0
  (1.6 false per replicate at 100 kb, K = 2,000; 2.85 at 1 Mb, K = 400). Exercised but not
  binding (FDR estimate < 0.002 at the end of the list). The FDP gate passes under the current
  rule (Slurm 1860995, 1 Mb, R = 300: FDP 0.0002 at every q). The 100-kb arm also asserts the
  estimate check mean(Ŵ0) >= mean(V0) - 2 SE: Ŵ0 2.44 vs V0 1.6, because the injected share
  dominates the pooled π̂. At 1 Mb the estimate check is a separate test
  (``test_database_stressed_stratum0_estimate_covers_false_placements``, ``calibration_full``,
  strict xfail "ADR-0019 implementation pending"). Mean PEP gives Ŵ0 0.57 vs V0 2.85 there. The
  no-origin share is length-correlated: 75 % of the injected peptides are 7-mers, and the 7-mer
  class has a no-origin share of 0.22 against 0.031 pooled. ADR-0019's per-class π̂ corrects
  this; the statistician's reference rule gives Ŵ0 3.65 vs V0 3.50 (review of #113, round 2).
* ``related_other_strain`` (5 SNV/kb, other-strain wrong peptides at their lifted ortholog):
  smoke check (all of them are placed at their origin).
* ``denovo``: smoke check (errors are placed at their spectrum's origin).
* ``related_dbstrain`` (decision 8, ADR-0019): correct peptides reported as the database strain's
  sequence; FDP smoke check at 100 kb, and the 1-Mb stratum-0 estimate check at 5 SNV/kb
  (``calibration_full``, strict xfail "ADR-0019 implementation pending").
"""

from __future__ import annotations

import shutil
from pathlib import Path

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _placement_suite import run_fdr_suite, write_diagnostics
from helpers.placement_sims import null_share_stratum0, run_chain, score_chain

pytestmark = pytest.mark.calibration

LEVELS = C.load_tolerances().levels
ADR0019 = pytest.mark.xfail(strict=True, reason="ADR-0019 implementation pending")


def _length(full: bool) -> int:
    return 1_000_000 if full else 100_000


def _proteomics(**kw: object) -> object:
    from giapsim.proteomics import ProteomicsConfig

    return ProteomicsConfig(search_db="other", **kw)  # type: ignore[arg-type]


def _strain(snv: float = 5.0) -> object:
    from giapsim.strain import StrainConfig

    return StrainConfig(snv_rate_per_kb=snv)


def _suite(full: bool, pi_wrong: float, tmp: Path) -> None:
    def replicate(seed: int, d: Path):  # type: ignore[no-untyped-def]
        ch = run_chain(
            "database", seed, d, length=_length(full), proteomics=_proteomics(pi_wrong=pi_wrong),
            paralogs=20,
        )  # fmt: skip
        return score_chain(ch, LEVELS, d)

    run_fdr_suite(f"placement.sim_db_pi{pi_wrong}", replicate, tmp, (), full=full, mutant=True)


@pytest.mark.parametrize("pi_wrong", [0.01, 0.05])
def test_database_shuffled_placement_fdr(tmp_path: Path, pi_wrong: float) -> None:
    _suite(False, pi_wrong, tmp_path)


@pytest.mark.calibration_full
@pytest.mark.parametrize("pi_wrong", [0.01, 0.05])
def test_database_shuffled_placement_fdr_full(tmp_path: Path, pi_wrong: float) -> None:
    _suite(True, pi_wrong, tmp_path)


def _stressed_replicate(full: bool, seed: int, d: Path):  # type: ignore[no-untyped-def]
    ch = run_chain(
        "database", seed, d, length=_length(full),
        proteomics=_proteomics(pi_wrong=0.05, n_runs=1), inject_short=400 if full else 2000,
    )  # fmt: skip
    return score_chain(ch, LEVELS, d)


def _stressed(full: bool, tmp: Path) -> None:
    def replicate(seed: int, d: Path):  # type: ignore[no-untyped-def]
        return _stressed_replicate(full, seed, d)

    # The 100-kb arm also carries the estimate check (it passes under mean PEP, docstring);
    # at 1 Mb that check is the separate strict xfail below.
    run_fdr_suite(
        "placement.sim_db_stressed", replicate, tmp, [("reverse", 0)], full=full,
        estimate=() if full else ("reverse",),
    )  # fmt: skip


def test_database_stressed_stratum0_placement_fdr(tmp_path: Path) -> None:
    """Review B4: the stratum-0 database list exercised at genome level (see module docstring)."""
    _stressed(False, tmp_path)


@pytest.mark.calibration_full
def test_database_stressed_stratum0_placement_fdr_full(tmp_path: Path) -> None:
    _stressed(True, tmp_path)


STRESS_EST_R = 100
STRESS_EST_V0_FLOOR = 100


@pytest.mark.calibration_full
@ADR0019
def test_database_stressed_stratum0_estimate_covers_false_placements(tmp_path: Path) -> None:
    """Statistician review of #113, round 2: the 1-Mb stressed arm (K = 400 injected short wrong
    database peptides, PEP = the true wrong share, delta = 0). Estimate rule as in the related
    db-strain check: mean(Ŵ0) >= mean(V0) - k SE, where Ŵ0 = fdr_stratum x targets of the reverse
    m = 0 ``mapping_fdr`` row (all length classes) and V0 = its false peptides, at any q. The floor
    is ΣV0 >= 100 (about 285 expected). The no-origin share is length-correlated (7-mer class
    0.22 against 0.031 pooled), so mean PEP under-counts: Ŵ0 0.57 vs V0 2.85, more than 10 SE at
    R = 100. ADR-0019's per-class π̂ corrects this (reference Ŵ0 3.65 vs V0 3.50)."""
    w0, v0 = [], []
    for i in range(STRESS_EST_R):
        seed = int(C.replicate_rng("placement.est_db_stressed", i).integers(0, 2**31 - 1))
        d = tmp_path / f"r{i}"
        rep = _stressed_replicate(True, seed, d)
        w0.append(float(rep.extra["W0_reverse"]))
        v0.append(float(rep.extra["V0_reverse"]))
        shutil.rmtree(d, ignore_errors=True)
    w, v = np.array(w0), np.array(v0)
    k = C.load_tolerances().se_multiplier
    se = float((v - w).std(ddof=1) / np.sqrt(len(v)))
    rec = {"suite": "placement.est_db_stressed", "n": len(v), "se": se}
    write_diagnostics({**rec, "mean_W0": float(w.mean()), "mean_V0": float(v.mean())})
    assert v.sum() >= STRESS_EST_V0_FLOOR, f"stratum 0 under-exercised: total V0 {v.sum():.0f}"
    assert w.mean() >= v.mean() - k * se, (
        f"estimated stratum-0 false placements {w.mean():.3f} < realised {v.mean():.3f} - {k} SE"
        f" ({se:.3f}), {len(v)} replicates"
    )


def _related(full: bool, tmp: Path) -> None:
    def replicate(seed: int, d: Path):  # type: ignore[no-untyped-def]
        ch = run_chain(
            "database", seed, d, length=_length(full), strain=_strain(),
            proteomics=_proteomics(pi_wrong=0.05, wrong_source="other_strain"),
            overrides=["genome_relation=related"],
        )  # fmt: skip
        return score_chain(ch, LEVELS, d)

    run_fdr_suite("placement.sim_related", replicate, tmp, (), full=full)


def test_related_other_strain_placement_fdr(tmp_path: Path) -> None:
    _related(False, tmp_path)


@pytest.mark.calibration_full
def test_related_other_strain_placement_fdr_full(tmp_path: Path) -> None:
    _related(True, tmp_path)


def _denovo(full: bool, tmp: Path) -> None:
    from giapsim.denovo import DenovoConfig

    def replicate(seed: int, d: Path):  # type: ignore[no-untyped-def]
        ch = run_chain(
            "denovo", seed, d, length=_length(full), proteomics=_proteomics(pi_wrong=0.05),
            denovo=DenovoConfig(),
        )  # fmt: skip
        return score_chain(ch, LEVELS, d)

    run_fdr_suite("placement.sim_denovo", replicate, tmp, (), full=full)


def test_denovo_error_model_placement_fdr(tmp_path: Path) -> None:
    _denovo(False, tmp_path)


@pytest.mark.calibration_full
def test_denovo_error_model_placement_fdr_full(tmp_path: Path) -> None:
    _denovo(True, tmp_path)


# ---------------------------------------------------------------- stratum-0 π̂ (database input)

PI_R = 20  # the shortfall is several-fold; 20 replicates resolve it (SE ~ 0.001)


def _pi_check(kind: str, stream: str, tmp: Path, **kw: object) -> tuple[float, float, float]:
    """(mean π̂, mean true null share, SE of the difference) over PI_R replicates; null share =
    stratum-0 searched database peptides that are wrong or absent from the sample."""
    pi, null = [], []
    for i in range(PI_R):
        seed = int(C.replicate_rng(stream, i).integers(0, 2**31 - 1))
        ch = run_chain(kind, seed, tmp / f"r{i}", **kw)  # type: ignore[arg-type]
        a, b = null_share_stratum0(ch)
        shutil.rmtree(tmp / f"r{i}", ignore_errors=True)
        pi.append(a)
        null.append(b)
    d = np.array(null) - np.array(pi)
    return float(np.mean(pi)), float(np.mean(null)), float(d.std(ddof=1) / np.sqrt(PI_R))


def _assert_pi(mean_pi: float, mean_null: float, se: float) -> None:
    k = C.load_tolerances().se_multiplier
    assert mean_pi >= mean_null - k * se, (
        f"stratum-0 π̂ {mean_pi:.5f} below the true share of wrong or absent searched peptides"
        f" {mean_null:.5f} (k SE = {k * se:.5f}); validity needs π̂ >= π (placement-fdr.md)"
    )


@pytest.mark.parametrize("n_runs", [1, pytest.param(3, marks=ADR0019)])
def test_database_stratum0_pi_hat_covers_the_wrong_share(tmp_path: Path, n_runs: int) -> None:
    """π̂ (mean input PEP over the stratum-0 database queries) must not understate their wrong
    share (methods note: "validity needs π̂ >= π"). M3-19 row PEPs are calibrated; with 3 runs
    preflight collapses a peptide's runs with ``pep = min`` and the mean of the minima is not a
    peptide-level wrong share (ADR-0019 "Run collapse"; decision 9: [3] is the regression, strict
    xfail until ADR-0019 replaces mean PEP; [1] is the control)."""
    stats = _pi_check(
        "database", f"placement.pi_db_runs{n_runs}", tmp_path,
        proteomics=_proteomics(pi_wrong=0.05, n_runs=n_runs),
    )  # fmt: skip
    _assert_pi(*stats)


def _dbstrain(length: int, snv: float) -> dict[str, object]:
    return {"length": length, "strain": _strain(snv),
            "proteomics": _proteomics(pi_wrong=0.01, n_runs=1)}  # fmt: skip


def _related_dbstrain(full: bool, tmp: Path) -> None:
    kw = _dbstrain(_length(full), 5.0)

    def replicate(seed: int, d: Path):  # type: ignore[no-untyped-def]
        ch = run_chain("related_dbstrain", seed, d, **kw)  # type: ignore[arg-type]
        return score_chain(ch, LEVELS, d)

    run_fdr_suite("placement.sim_related_dbstrain", replicate, tmp, (), full=full)


def test_related_dbstrain_placement_fdr(tmp_path: Path) -> None:
    """Decision 8 (replaces the decision-7 xfail): related genome, ``search_db = other``, correct
    peptides reported as the database (M3-18 reference) strain's sequence at the same protein
    offset (SNV-only strain, 5 SNV/kb; one run). Peptides spanning a strain difference have no exact
    origin; their chance exact hits elsewhere are the stratum-0 null, about 0 per replicate at
    100 kb (review B1: L/I SNV pairs keep their origin), so this is a smoke check; the 1-Mb
    estimate check below is the evidence."""
    _related_dbstrain(False, tmp_path)


@pytest.mark.calibration_full
def test_related_dbstrain_placement_fdr_full(tmp_path: Path) -> None:
    _related_dbstrain(True, tmp_path)


EST_R = 100  # review B3: 1-Mb genomes, about 13 s each; SE ~ 0.09 at V0 ~ 0.75
EST_SNV = 5.0  # review B3 (decision 8 allows >= 2/kb): V0 ~ 0.75, D0 ~ 9 per replicate
EST_V0_FLOOR = 30


@pytest.mark.calibration_full
@ADR0019
def test_related_dbstrain_stratum0_estimate_covers_false_placements(tmp_path: Path) -> None:
    """Decision 8 estimate check (statistician, B2; review B3): 1-Mb genome, 5 SNV/kb, correct
    peptides from the reference strain as above. Per replicate and decoy type at m = 0 (database
    input: reverse; the shuffle list is empty), Ŵ0 = ``fdr_stratum`` x ``targets`` of the
    ``mapping_fdr`` row (all length classes) is the estimated number of false placements among
    all stratum-0 targets, V0 the realised number (no best-cost locus overlaps a true origin).
    Assert mean(Ŵ0) >= mean(V0) - k SE(V0 - Ŵ0) over EST_R replicates and total V0 >=
    EST_V0_FLOOR. Mean PEP counts only identification errors, so Ŵ0 under-states V0
    (statistician: Ŵ0 ~ 0.1 vs V0 ~ 0.75); strict xfail until ADR-0019 is implemented (measured
    ratio 0.96-2.4 there)."""
    from helpers.placement_sims import hits_of
    from helpers.placement_truth import best_loci

    kw = _dbstrain(1_000_000, EST_SNV)
    w0, v0 = [], []
    for i in range(EST_R):
        seed = int(C.replicate_rng("placement.est_related_dbstrain", i).integers(0, 2**31 - 1))
        d = tmp_path / f"r{i}"
        ch = run_chain("related_dbstrain", seed, d, **kw)  # type: ignore[arg-type]
        per = pd.read_parquet(ch.avail["placement"]).drop_duplicates("peptide")
        loci = best_loci(hits_of(ch.avail), dict(zip(per["peptide"], per["m_target"].astype(int),
                                                     strict=True)))  # fmt: skip
        f = pd.read_parquet(ch.avail["mapping_fdr"])
        for dt in ("reverse", "shuffle"):
            row = f[(f["section"] == "placement") & (f["decoy_type"] == dt) & (f["n_diff"] == 0)
                    & (f["length_class"] == "all")]  # fmt: skip
            if row.empty or int(row["targets"].iloc[0]) == 0:
                continue
            t0 = per[(per["m_target"] == 0) & (per["decoy_type"] == dt)]["peptide"]
            v0.append(sum(not any(ch.truth.overlaps(iv, p) for iv in loci.get(p, ())) for p in t0))
            w0.append(float(row["fdr_stratum"].iloc[0]) * int(row["targets"].iloc[0]))
        shutil.rmtree(d, ignore_errors=True)
    w, v = np.array(w0), np.array(v0, dtype=float)
    k = C.load_tolerances().se_multiplier
    se = float((v - w).std(ddof=1) / np.sqrt(len(v)))
    rec = {"suite": "placement.est_related_dbstrain", "n": len(v), "se": se}
    write_diagnostics({**rec, "mean_W0": float(w.mean()), "mean_V0": float(v.mean())})
    assert v.sum() >= EST_V0_FLOOR, f"stratum 0 under-exercised: total V0 {v.sum():.0f}"
    assert w.mean() >= v.mean() - k * se, (
        f"estimated stratum-0 false placements {w.mean():.3f} < realised {v.mean():.3f} - {k} SE"
        f" ({se:.3f}), {len(v)} lists"
    )
