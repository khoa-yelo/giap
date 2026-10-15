"""Shared replicate machinery for the M4-17 region-test calibration and attack suites.

One simulated genome per replicate feeds every (region, statistic) check: each check has its
own ``_criteria`` stream (so its log record and block-2 confirmation are its own), and
replicate ``i`` of any stream maps back to the scenario's replicate ``i`` (seed
``suite_seed(<scenario stream>, i)``), simulated once per process and cached. Region tests run with
``giap.rng.for_step("context_stats", seed)`` (spec "RNG"), ``seed`` drawn from the replicate.

A check's replicate value is the p-value of one region row; NA (no valid shift, or no detected
gene for HL) counts as p = 1 (``p_or_one``: a test that cannot be computed does not reject).
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
import region_tests_data as rd

#: bootstrap B for calibration runs: the CI is not under test here, the p-values are
CAL_BOOTSTRAP = 20


@dataclass(frozen=True)
class Scenario:
    """A null generator: layout size, expression knobs and the p_detect source."""

    name: str
    n_genes: int
    expr_kw: Mapping[str, Any] = field(default_factory=dict)
    p_source: str = "true"
    layout_kw: Mapping[str, Any] = field(default_factory=dict)

    @property
    def stream(self) -> str:
        return f"regions.region_tests.{self.name}"


_TABLES: dict[tuple[str, int], pd.DataFrame] = {}


def _C() -> Any:
    import _criteria as C

    return C


def replicate_results(sc: Scenario, i: int) -> pd.DataFrame:
    """``test_all_regions`` output of replicate ``i`` of scenario ``sc`` (cached)."""
    key = (sc.stream, i)
    if key not in _TABLES:
        from giap.rng import for_step

        rng = _C().replicate_rng(sc.stream, i)
        lc = rd.calibration_layout_config(sc.n_genes, **dict(sc.layout_kw))
        tables = rd.sim_tables(rng, lc, sc.expr_kw, p_source=sc.p_source)
        seed = int(rng.integers(0, 2**31 - 1))
        cfg = rd.config(bootstrap=CAL_BOOTSTRAP)
        _TABLES[key] = rd.run_region_tests(tables, cfg, for_step("context_stats", seed))
    return _TABLES[key]


def region_p(res: pd.DataFrame, region_type: str, stat: str, contig: str) -> float:
    """p of the first region (sorted region_id) of ``region_type`` on ``contig``, NA -> 1."""
    rows = res[(res["region_type"] == region_type) & (res["stat"] == stat)]
    rows = rows[rows["contig"] == contig].sort_values("region_id")
    if rows.empty:
        raise AssertionError(f"no {region_type}/{stat} row on {contig}")
    return rd.p_or_one(rows.iloc[0]["p_value"])


def shared_fn(
    sc: Scenario, stream: str, R: int, value: Callable[[pd.DataFrame], Any]
) -> Callable[[np.random.Generator], Any]:
    """Replicate function for ``calibrate(..., stream=stream)``: maps the generator handed in by
    ``_criteria`` (seed ``suite_seed(stream, i)``) to replicate ``i`` of ``sc``."""
    C = _C()
    index = {C.suite_seed(stream, i): i for i in range(2 * R)}

    def fn(rng: np.random.Generator) -> Any:
        key = int(rng.bit_generator.seed_seq.entropy)  # type: ignore[attr-defined]
        return value(replicate_results(sc, index[key]))

    return fn


def level_checks(sc: Scenario, checks: list[tuple[str, str, str]], R: int) -> list[str]:
    """Level at every ``levels`` alpha for each (region_type, stat, contig) check of ``sc``
    (stream ``<scenario stream>.<type>[@contig].<stat>``); returns failure details."""
    from _gate import calibrate_levels

    C = _C()
    levels = list(C.load_tolerances().levels)
    failed: list[str] = []
    for typ, stat, contig in checks:
        tag = typ if contig == "chr1" else f"{typ}@{contig}"
        stream = f"{sc.stream}.{tag}.{stat}"
        fn = shared_fn(sc, stream, R, lambda res, t=typ, s=stat, c=contig: region_p(res, t, s, c))
        gate = calibrate_levels("regions", "level", levels, fn, R, stream=stream)
        failed += gate.detail() if gate.failed else []
    return failed


CHROMOSOME_CHECKS = [(t, s, "chr1") for t in rd.REGION_TYPES for s in rd.STATS]
PLASMID_AMR_CHECKS = [("amr", s, "p1") for s in rd.STATS]
