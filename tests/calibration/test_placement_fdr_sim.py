"""Genome-level placement FDR on giapsim through the mapping steps (issue M3-21; GIAP_SPEC.md
§7.2.9 "Calibration: placement FDR <= q + Monte Carlo error on simulated wrong peptides (database
and de novo error models)", §8.4 "Placement FDR" row; docs/specs/M3-21 lead decisions 1-7).

Every replicate: a giapsim genome (100 kb fast / 1 Mb full), M3-19 proteomics, ``giap
preflight`` and the registered mapping steps (helpers.placement_sims). FDP per list (decoy type x
m) and pooled, false placement = no accepted locus overlaps a true origin (helpers.placement_truth),
criterion mean FDP <= q + 2 SE at every level (``_placement_suite.run_fdr_suite``).

Suites and what they exercise (mean wrong or non-origin peptides with a target hit per replicate,
measured at 100 kb, 2026-10-20):

* ``database_shuffled`` (pi_wrong 0.01 / 0.05, ``search_db = other``, 20 planted near-paralogs for
  the unique guard and its mutant): **no list is exercised** -- after the q <= q_value_max filter
  about 1-2 % of the mapped queries are wrong and none hits the 100-kb target by chance. Nothing is
  claimed; every list and the pooled FDP are still gated. The stratum-0 risk is checked on π̂
  directly (``test_database_stratum0_pi_hat_covers_the_wrong_share``).
* ``related_other_strain``: M3-18 strain (5 SNV/kb), other-strain wrong peptides lifted to their
  ortholog (decision 2) -> list reverse m = 1 (about 30 per replicate).
* ``denovo``: M3-20 reads, shuffle decoys -> shuffle m = 1 (about 1,500) and m = 2 (about 60).
* ``related_dbstrain`` (decision 7): correct peptides reported as the database strain's sequence;
  strict xfail "ADR pending: stratum-0 π̂ for search_db=other".
"""

from __future__ import annotations

import shutil
from pathlib import Path

import _criteria as C
import numpy as np
import pytest
from _placement_suite import run_fdr_suite
from helpers.placement_sims import null_share_stratum0, run_chain, score_chain

pytestmark = pytest.mark.calibration

LEVELS = C.load_tolerances().levels


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


def _related(full: bool, tmp: Path) -> None:
    def replicate(seed: int, d: Path):  # type: ignore[no-untyped-def]
        ch = run_chain(
            "database", seed, d, length=_length(full), strain=_strain(),
            proteomics=_proteomics(pi_wrong=0.05, wrong_source="other_strain"),
            overrides=["genome_relation=related"],
        )  # fmt: skip
        return score_chain(ch, LEVELS, d)

    run_fdr_suite("placement.sim_related", replicate, tmp, [("reverse", 1)], full=full)


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

    claimed = [("shuffle", 1), ("shuffle", 2)]
    run_fdr_suite("placement.sim_denovo", replicate, tmp, claimed, full=full)


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


@pytest.mark.parametrize("n_runs", [1, 3])
def test_database_stratum0_pi_hat_covers_the_wrong_share(tmp_path: Path, n_runs: int) -> None:
    """π̂ (mean input PEP over the stratum-0 database queries) must not understate their wrong
    share (methods note: "validity needs π̂ >= π"). M3-19 row PEPs are calibrated; with 3 runs
    the per-peptide input PEP is aggregated over the runs."""
    stats = _pi_check(
        "database", f"placement.pi_db_runs{n_runs}", tmp_path,
        proteomics=_proteomics(pi_wrong=0.05, n_runs=n_runs),
    )  # fmt: skip
    _assert_pi(*stats)


@pytest.mark.xfail(strict=True, reason="ADR pending: stratum-0 π̂ for search_db=other")
def test_related_dbstrain_correct_peptides_stratum0(tmp_path: Path) -> None:
    """Decision 7 (statistician B2): related genome, ``search_db = other``, correct peptides drawn
    from the database strain (one run, so input PEPs are the calibrated row PEPs). Correct
    peptides crossing a strain difference have no exact origin in the sample, so for stratum 0
    they are null peptides that the mean PEP does not count. Checks: the placement FDR gates
    (all lists and pooled), then π̂ >= the true null share of the stratum-0 database queries."""

    def replicate(seed: int, d: Path):  # type: ignore[no-untyped-def]
        ch = run_chain(
            "related_dbstrain", seed, d, strain=_strain(),
            proteomics=_proteomics(pi_wrong=0.01, n_runs=1),
        )  # fmt: skip
        return score_chain(ch, LEVELS, d)

    run_fdr_suite("placement.sim_related_dbstrain", replicate, tmp_path / "fdr", ())
    stats = _pi_check(
        "related_dbstrain", "placement.pi_related_dbstrain", tmp_path / "pi", strain=_strain(),
        proteomics=_proteomics(pi_wrong=0.01, n_runs=1),
    )  # fmt: skip
    _assert_pi(*stats)
