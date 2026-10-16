"""Calibration harness fixtures and suite selection (issue M4-05; GIAP_SPEC.md §8.4, §11.2).

- ``full``: true for tests marked ``calibration_full`` (nightly R), false otherwise.
- Suite assignment is explicit (PR #54 lead decision 4): every calibration test belongs to exactly
  one ``tolerances.yaml`` suite, taken from (first match wins) a
  ``@pytest.mark.calsuite("<suite>")`` marker, ``SUITE_REGISTRY["<file>::<function>"]`` or
  ``SUITE_REGISTRY["<file>"]`` (paths relative to tests/calibration). The test gets the keyword
  ``calsuite_<suite>`` so the job array selects one suite with ``-k calsuite_<suite>``
  (``_criteria.pytest_args``). A calibration test in no entry
  runs in no array task; ``test_criteria_review.py`` fails on such orphans, so register new files.
- Results are appended by ``_criteria`` to ``local/logs/calibration/<suite>.jsonl``.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from types import MappingProxyType

import _criteria as C
import pytest

HERE = Path(__file__).resolve().parent

# Engine and attack tests use the suite whose R they draw (replicates_fast/full.<suite>); tests
# with no R of their own sit with their consumer: chance hits and exact recall with placement, the
# giapsim v2 knobs (promoters, terminators, TUs) with missed_operons, harness self-tests with
# regions.
SUITE_REGISTRY: Mapping[str, str] = MappingProxyType(
    {
        "test_abundance_calibration.py": "regions",
        "test_abundance_calibration.py::test_ibaq_denominator_beats_raw_sum": "regions",
        "test_abundance_calibration.py::test_teeth_permuted_denominator_fails": "regions",
        "test_chance_hits.py": "placement",
        "test_criteria.py": "regions",
        "test_criteria_review.py": "regions",
        "test_criteria_pooled.py": "regions",
        "test_denovo_precision_calibration.py": "denovo_precision",
        "test_detection_tiers_gene_fdr.py": "gene_fdr",
        "test_detectability_calibration.py": "detectability",
        "test_engine_calibration.py": "regions",
        "test_engine_calibration.py::test_bh_fdr": "gene_fdr",
        "test_exact_recall.py": "placement",
        "test_tolerant_recall.py": "placement",
        "test_mapping_steps_gate.py": "placement",
        "test_mayu_calibration.py": "gene_fdr",
        "test_placement_fdr_genome.py": "placement",
        "test_placement_fdr_lists.py": "placement",
        "test_sim_v2_calibration.py": "missed_operons",
        "test_sim_layout_calibration.py": "missed_operons",
        "test_sim_denovo_calibration.py": "denovo_precision",
        "test_sim_proteomics_calibration.py": "gene_fdr",
        "test_sim_evidence_calibration.py": "gene_fdr",
        "test_sim_expression_calibration.py": "regions",
        **{
            f"test_sim_expression_calibration.py::{name}": "detectability"
            for name in (
                "test_detection_bernoulli_total",
                "test_detection_bernoulli_by_decile",
                "test_covariates_only_detection_independent_of_abundance",
                "test_bgc_detection_shift",
                "test_default_detected_share_every_replicate",
            )
        },
        "test_stats_bootstrap_calibration.py": "regions",
        "test_detection_summary_calibration.py": "regions",
        "test_stats_fdr_calibration.py": "gene_fdr",
        "test_stats_permutation_calibration.py": "regions",
        "attacks/test_attack_bootstrap.py": "regions",
        "attacks/test_detection_summary_attack_low_e.py": "regions",
        "attacks/test_attack_calibration_harness.py": "regions",
        "attacks/test_attack_calibration_harness.py::test_fdr_check_passes_exact_level_bh": (
            "gene_fdr"
        ),
        "attacks/test_attack_calibration_harness.py::test_fdr_check_flags_inflated_bh": (
            "gene_fdr"
        ),
        "attacks/test_attack_chance_hits.py": "placement",
        "attacks/test_attack_denovo_precision.py": "placement",
        "attacks/test_detection_tiers_attack_decoy_loss.py": "gene_fdr",
        "attacks/test_attack_detectability.py": "detectability",
        "attacks/test_attack_fdr_engine.py": "gene_fdr",
        "attacks/test_attack_genome_decoys.py": "placement",
        "attacks/test_attack_mayu.py": "gene_fdr",
        "attacks/test_attack_permutation.py": "regions",
        "attacks/test_attack_placement_pi.py": "placement",
        "attacks/test_attack_placement_same_db.py": "placement",
    }
)


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers", "calsuite(name): tolerances.yaml suite (array task) of a calibration test"
    )


@pytest.fixture
def full(request: pytest.FixtureRequest) -> bool:
    """True under ``-m calibration_full`` (the test carries the marker)."""
    return request.node.get_closest_marker("calibration_full") is not None


@pytest.fixture
def teeth_log(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    """Log dir for teeth tests (PR #69 decision 3): an anti-conservative procedure must be a
    confirmed ``fail`` through ``calibrate``, so it runs unsharded (both blocks, whatever
    ``--shard`` / ``--merge`` set) and logs to a private directory, never the nightly suite log."""
    for var in (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV):
        monkeypatch.delenv(var, raising=False)
    return tmp_path


def item_suite(item: pytest.Function) -> str | None:
    """The one suite of ``item`` (marker, then function entry, then file entry); None if unset."""
    marker = item.get_closest_marker("calsuite")
    if marker is not None:
        suite = str(marker.args[0])
    else:
        try:
            rel = Path(str(item.path)).resolve().relative_to(HERE).as_posix()
        except ValueError:
            return None
        name = item.originalname
        suite_or_none = SUITE_REGISTRY.get(f"{rel}::{name}", SUITE_REGISTRY.get(rel))
        if suite_or_none is None:
            return None
        suite = suite_or_none
    if suite not in C.array_suites():
        raise pytest.UsageError(f"{item.nodeid}: calsuite {suite!r} is not in tolerances.yaml")
    return suite


# Fast-mode calibration tests slower than 60 s (single core; Slurm 1851461, 2026-10-16): marked
# ``calibration_nightly``, so GitHub CI (4-core runners, budget 15-30 min) skips them and
# scripts/slurm/nightly.sbatch runs them every night on main. Add a test here when it passes 60 s.
NIGHTLY: frozenset[str] = frozenset(
    {
        "attacks/test_attack_bootstrap.py::test_recorded_few_tus",
        "attacks/test_attack_bootstrap.py::test_tu_block_ci_coverage_under_attack",
        "attacks/test_attack_detectability.py::test_attack2_no_size_covariates_is_confirmed_fail",
        "attacks/test_attack_detectability.py::test_attack2_subgroup_oe_at_level",
        "attacks/test_attack_genome_decoys.py::test_composition_null_decoy_hits_equal_target_hits",
        "attacks/test_detection_summary_attack_low_e.py::test_low_e_group_is_flagged",
        "attacks/test_detection_summary_attack_low_e.py::test_low_e_unflagged_miss_rate_exceeds_level",
        "attacks/test_detection_tiers_attack_decoy_loss.py::test_picked_fdr_under_random_decoy_loss",
        "attacks/test_detection_tiers_attack_inputs.py::test_gene_fdr_under_input_attack",
        "attacks/test_region_tests_attacks.py::test_a2_detection_by_abundance_prophage_hl_at_level",
        "attacks/test_region_tests_attacks.py::test_a2_mid_range_prophage_p_o_minus_e_at_level",
        "attacks/test_region_tests_attacks.py::test_a2_short_prophage_genes_o_minus_e_at_level",
        "attacks/test_region_tests_attacks.py::test_a2_short_prophage_hl_at_level_without_abundance_selection",
        "attacks/test_region_tests_attacks.py::test_a2_teeth_constant_p_detect_rejects",
        "attacks/test_region_tests_attacks.py::test_alin_linear_replicon_hl_at_level",
        "attacks/test_region_tests_attacks.py::test_atu_region_specific_tu_fragmentation_hl_at_level",
        "test_detectability_calibration.py::test_d1_mean_cross_fitted_slope_in_tolerance",
        "test_detection_summary_calibration.py::test_o_minus_e_ci_covers_zero_under_tu_correlated_null",
        "test_detection_tiers_gene_fdr.py::test_gene_fdr_controlled",
        "test_detection_tiers_gene_fdr.py::test_teeth_no_q_filter_fails",
        "test_mapping_steps_gate.py::test_placement_fdr_through_the_steps",
        "test_mapping_steps_gate.py::test_recall_of_planted_peptides_through_the_steps",
        "test_placement_fdr_genome.py::test_realised_placement_fdr_on_giapsim",
        "test_region_tests.py::test_c1_level_fast",
        "test_region_tests.py::test_c2_fdr_within_type_and_stat",
        "test_region_tests.py::test_c3_power_recorded",
        "test_sim_evidence_calibration.py::test_false_targets_size_proportional",
        "test_sim_evidence_calibration.py::test_intensity_lognormal",
        "test_sim_evidence_calibration.py::test_pep_sum_matches_false_count",
        "test_sim_evidence_calibration.py::test_q_values_control_fdp",
        "test_sim_evidence_calibration.py::test_score_components",
        "test_stats_bootstrap_calibration.py::test_recorded_coverage_twenty_tus",
        "test_tolerant_recall.py::test_planted_tolerant_peptide_recall",
    }
)


@pytest.hookimpl(tryfirst=True)
def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    for item in items:
        if not isinstance(item, pytest.Function):
            continue
        if not any(item.iter_markers("calibration")) and not any(
            item.iter_markers("calibration_full")
        ):
            continue
        suite = item_suite(item)
        if suite is not None:
            item.extra_keyword_matches.add(f"calsuite_{suite}")
        try:
            rel = Path(str(item.path)).resolve().relative_to(HERE).as_posix()
        except ValueError:
            continue
        if f"{rel}::{item.originalname}" in NIGHTLY and any(item.iter_markers("calibration")):
            item.add_marker(pytest.mark.calibration_nightly)
