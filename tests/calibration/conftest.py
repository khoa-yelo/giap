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
- Confirm runs (``GIAP_CALIBRATION_CONFIRM``, set only by ``--merge``'s block-2 run; M7-13 lead
  decisions 7-8): a calibration test that cannot reach the harness skips before setup, and a
  harness test that stops before its first harness call skips (its block-1 path was judged in
  the shards). See ``reaches_harness``.
"""

from __future__ import annotations

import os
import sys
import types
from collections.abc import Callable, Generator, Iterator, Mapping
from pathlib import Path
from types import MappingProxyType
from typing import Any

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
        "test_denovo_precision_sim.py": "denovo_precision",
        "test_detection_tiers_gene_fdr.py": "gene_fdr",
        "test_gene_fdr_real_mapping.py": "gene_fdr",
        "test_detectability_calibration.py": "detectability",
        "test_engine_calibration.py": "regions",
        "test_engine_calibration.py::test_bh_fdr": "gene_fdr",
        "test_exact_recall.py": "placement",
        "test_tolerant_recall.py": "placement",
        "test_mapping_steps_gate.py": "placement",
        "test_mayu_calibration.py": "gene_fdr",
        "test_placement_fdr_genome.py": "placement",
        "test_placement_fdr_lists.py": "placement",
        "test_placement_fdr_sim.py": "placement",
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
        "test_stats_scan_calibration.py": "islands",
        "test_islands_calibration.py": "islands",
        "attacks/test_attack_islands.py": "islands",
        "test_region_tests.py": "regions",
        "test_autocorrelation_calibration.py": "autocorrelation",
        "test_oriter_calibration.py": "islands",
        "test_operon_fdr_calibration.py": "missed_operons",
        "test_operon_tdist_calibration.py": "discordance",
        "attacks/test_region_tests_attacks.py": "regions",
        "attacks/test_regions_stats_attacks.py": "regions",
        "attacks/test_attack_bootstrap.py": "regions",
        "attacks/test_attack_autocorr.py": "autocorrelation",
        "attacks/test_attack_autocorr_engine.py": "autocorrelation",
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
        "attacks/test_attack_placement_contaminants.py": "placement",
        "attacks/test_attack_placement_pi.py": "placement",
        "attacks/test_attack_placement_same_db.py": "placement",
        "attacks/test_attack5_reversed_stretches.py": "placement",
        "attacks/test_attack_placement_strain_span.py": "placement",
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
        "attacks/test_attack5_reversed_stretches.py::test_attack5_reversed_stretches_shuffle_decoys",
        "attacks/test_attack_autocorr.py::test_attack_oriter_gradient_gene_label_at_level",
        "attacks/test_attack_autocorr.py::test_attack_oriter_gradient_tu_block_at_level",
        "attacks/test_attack_autocorr.py::test_attack_short_prophage_join_count_at_level",
        "attacks/test_attack_autocorr.py::test_teeth_permuting_residuals_without_redetrending_fails",
        "attacks/test_attack_bootstrap.py::test_recorded_few_tus",
        "attacks/test_attack_bootstrap.py::test_tu_block_ci_coverage_under_attack",
        "attacks/test_attack_detectability.py::test_attack2_no_size_covariates_is_confirmed_fail",
        "attacks/test_attack_detectability.py::test_attack2_subgroup_oe_at_level",
        "attacks/test_regions_stats_attacks.py::test_recorded_plasmid_hl_bias_and_ci_coverage",
        "attacks/test_attack_genome_decoys.py::test_composition_null_decoy_hits_equal_target_hits",
        "attacks/test_attack_islands.py::test_attack_draft_assembly",
        "attacks/test_attack_islands.py::test_attack_many_replicons",
        "attacks/test_attack_islands.py::test_attack_oriter_gradient_detection_abundance",
        "attacks/test_attack_islands.py::test_attack_oriter_gradient_scenario",
        "attacks/test_attack_islands.py::test_recorded_autocorrelation_rho09_sizes",
        "attacks/test_attack_islands.py::test_recorded_oriter_gradient_extreme",
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
        "test_autocorrelation_calibration.py::test_label_null_lag_level",
        "test_autocorrelation_calibration.py::test_label_null_bh_family_any_rejection",
        "test_autocorrelation_calibration.py::test_label_null_pooled_rate_recorded",
        "test_autocorrelation_calibration.py::test_tu_null_lag_level",
        "test_autocorrelation_calibration.py::test_tu_null_bh_family_any_rejection",
        "test_denovo_precision_sim.py::test_denovo_precision_on_giapsim",
        "test_detectability_calibration.py::test_d1_mean_cross_fitted_slope_in_tolerance",
        "test_detection_summary_calibration.py::test_o_minus_e_ci_covers_zero_under_tu_correlated_null",
        "test_detection_tiers_gene_fdr.py::test_gene_fdr_controlled",
        "test_detection_tiers_gene_fdr.py::test_teeth_no_q_filter_fails",
        "test_islands_calibration.py::test_genome_wide_family_over_two_replicons",
        "test_islands_calibration.py::test_island_fwer_at_level",
        "test_islands_calibration.py::test_teeth_spec_block_length_is_anticonservative",
        "test_mapping_steps_gate.py::test_placement_fdr_through_the_steps",
        "test_mapping_steps_gate.py::test_recall_of_planted_peptides_through_the_steps",
        "test_oriter_calibration.py::test_pooled_slope_ci_covers_the_planted_curve",
        "test_placement_fdr_genome.py::test_realised_placement_fdr_on_giapsim",
        "test_placement_fdr_sim.py::test_database_shuffled_placement_fdr",
        "test_placement_fdr_sim.py::test_database_stressed_stratum0_placement_fdr",
        "test_placement_fdr_sim.py::test_database_stratum0_pi_hat_covers_the_wrong_share",
        "test_placement_fdr_sim.py::test_denovo_error_model_placement_fdr",
        "test_placement_fdr_sim.py::test_related_dbstrain_placement_fdr",
        "test_placement_fdr_sim.py::test_related_other_strain_placement_fdr",
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
        "test_variant_class_error.py::test_class_error_rate_unbiased_with_planted_variants",
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


# --------------------------------------------------------------------------- confirm runs
# M7-13 lead decisions 7-8. The block-2 confirm run reruns only the harness checks listed in
# GIAP_CALIBRATION_CONFIRM. A test that asserts directly (G2, A6a, the reference arm, M7-10's
# class-error tests) has no block 2: its verdict is its shard exit code, so rerunning it there
# repeats block 1 and its result would be thrown away. Such tests skip before setup. Harness
# tests are found statically (no marker to maintain): the test function, its fixtures and the
# test-side functions they reach (tests/ only, never _criteria itself) name ``calibrate`` or
# ``calibrate_pooled``. A harness test still runs; if it stops before its first harness call
# (for example A6c's block-1 vacuity ``RuntimeError``, already raised in the shard), it skips.
# A miss here fails closed: a harness test skipped by mistake leaves its block 2 unlogged, and
# ``--merge`` reports MERGE_ERROR.

HARNESS_NAMES: frozenset[str] = frozenset({"calibrate", "calibrate_pooled"})
TESTS_ROOT = HERE.parent
CONFIRM_SKIP = "confirm run: no harness check (its verdict is the shard exit code)"
CONFIRM_SKIP_EARLY = "confirm run: stopped before the harness (block 1 is judged in the shards)"


def _codes(code: types.CodeType) -> Iterator[types.CodeType]:
    yield code
    for const in code.co_consts:
        if isinstance(const, types.CodeType):
            yield from _codes(const)


def _test_side(obj: object) -> types.FunctionType | None:
    """``obj`` as a plain function defined under tests/ (not in _criteria), else None."""
    obj = getattr(obj, "__wrapped__", obj)
    if not isinstance(obj, types.FunctionType):
        return None
    path = Path(obj.__code__.co_filename).resolve()
    if path == Path(C.__file__).resolve() or not path.is_relative_to(TESTS_ROOT):
        return None
    return obj


def _is_test_side_module(mod: types.ModuleType) -> bool:
    path = getattr(mod, "__file__", None)
    return path is not None and Path(path).resolve().is_relative_to(TESTS_ROOT)


def reaches_harness(roots: list[Callable[..., Any]]) -> bool:
    """True if any function reachable from ``roots`` names a harness entry point."""
    todo = [f for f in map(_test_side, roots) if f is not None]
    seen: set[types.CodeType] = set()
    while todo:
        fn = todo.pop()
        if fn.__code__ in seen:
            continue
        seen.add(fn.__code__)
        names = {n for code in _codes(fn.__code__) for n in code.co_names}
        if names & HARNESS_NAMES:
            return True
        glob = fn.__globals__
        mods = [glob[n] for n in names if isinstance(glob.get(n), types.ModuleType)]
        mods += [sys.modules[n] for n in names if n in sys.modules]  # function-local imports
        mods = [m for m in mods if _is_test_side_module(m)]  # never getattr on numpy & co.
        cands = [glob.get(n) for n in names]
        cands += [getattr(m, n, None) for m in mods for n in names]
        todo += [f for f in map(_test_side, cands) if f is not None]
    return False


def _is_calibration(item: pytest.Item) -> bool:
    return isinstance(item, pytest.Function) and (
        any(item.iter_markers("calibration")) or any(item.iter_markers("calibration_full"))
    )


def _item_reaches_harness(item: pytest.Function) -> bool:
    fixtures = [d.func for defs in item._fixtureinfo.name2fixturedefs.values() for d in defs]
    return reaches_harness([item.function, *fixtures])


def pytest_runtest_setup(item: pytest.Item) -> None:
    if C.CONFIRM_ENV not in os.environ or not _is_calibration(item):
        return
    assert isinstance(item, pytest.Function)
    if not _item_reaches_harness(item):
        pytest.skip(CONFIRM_SKIP)


@pytest.hookimpl(wrapper=True)
def pytest_pyfunc_call(pyfuncitem: pytest.Function) -> Generator[None, object, object]:
    if C.CONFIRM_ENV not in os.environ or not _is_calibration(pyfuncitem):
        return (yield)
    calls = [0]
    saved = {n: getattr(C, n) for n in HARNESS_NAMES}

    def counted(f: Callable[..., Any]) -> Callable[..., Any]:
        def call(*a: Any, **k: Any) -> Any:
            calls[0] += 1
            return f(*a, **k)

        return call

    for n, f in saved.items():
        setattr(C, n, counted(f))
    try:
        return (yield)
    except Exception:
        if calls[0] == 0 and C.CONFIRM_ENV in os.environ:
            pytest.skip(CONFIRM_SKIP_EARLY)
        raise
    finally:
        for n, f in saved.items():
            setattr(C, n, f)
