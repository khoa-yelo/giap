"""ADR-0019 decision 4 through the registered ``placement`` step: the run warning
``input_prefiltered`` (an input already filtered to genome matches cannot be calibrated by any
genome-based π̂; ADR-0019 Evidence, "Limitation").

Convention (as for ``c_assumed``, tests/unit/test_mapping_steps.py): a run warning is a WARNING
record on a ``giap.steps.*`` logger whose message names the flag. Inputs: a 50-kb giapsim genome,
database input (``search_db = other``, one run, so the row-level mean PEP and the peptide-level one
agree). Control: the stressed arm (``inject_short`` = 200 short wrong peptides absent from the
genome, PEP = the true wrong share): k (stratum-0 database peptides without an exact target hit)
is near n * mean PEP * (1 - r), no warning is due. Pre-filtered (giapsim's wrong identifications
dropped, every PEP = 0.05, ``helpers.placement_sims.prefilter_to_genome``): k is ~0 against an
expected ~35, far below the 0.001 binomial quantile (~20), so the warning must fire.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from helpers.placement_sims import Chain, run_chain

SEED = 20261022


class _Sink(logging.Handler):
    def __init__(self) -> None:
        super().__init__(logging.INFO)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


def _chain(tmp: Path, prefilter: bool) -> tuple[Chain, list[logging.LogRecord]]:
    from giapsim.proteomics import ProteomicsConfig

    kw = {"prefilter": 0.05} if prefilter else {"inject_short": 200}

    sink = _Sink()
    log = logging.getLogger("giap.steps")
    old = log.level
    log.addHandler(sink)
    log.setLevel(logging.INFO)
    try:
        ch = run_chain(
            "database", SEED, tmp, length=50_000,
            proteomics=ProteomicsConfig(pi_wrong=0.05, search_db="other", n_runs=1), **kw,
        )  # fmt: skip
    finally:
        log.removeHandler(sink)
        log.setLevel(old)
    return ch, sink.records


@pytest.fixture(scope="module")
def chains(tmp_path_factory: pytest.TempPathFactory) -> Iterator[dict[str, tuple]]:
    yield {
        "raw": _chain(tmp_path_factory.mktemp("raw"), prefilter=False),
        "prefiltered": _chain(tmp_path_factory.mktemp("prefiltered"), prefilter=True),
    }


def _warned(records: list[logging.LogRecord]) -> bool:
    return any(
        r.levelno >= logging.WARNING and "input_prefiltered" in r.getMessage()
        for r in records
        if r.name.startswith("giap.steps")
    )


def _stratum0(ch: Chain) -> tuple[int, int, float]:
    """(n, k, mean PEP) of the database stratum-0 population, from the step's own outputs."""
    from giap.steps.mapping import select_queries

    q = select_queries(pd.read_parquet(ch.avail["peptides"]), ch.cfg)
    q = q[q["source"] == "database"]
    per = pd.read_parquet(ch.avail["placement"])
    exact = set(per.loc[per["m_target"] == 0, "peptide"])
    peps = q.drop_duplicates("peptide")
    n = len(peps)
    k = int((~peps["peptide"].isin(exact)).sum())
    pep = pd.to_numeric(peps["pep"], errors="coerce").to_numpy(dtype=float)
    return n, k, float(np.nanmean(pep))


def test_unfiltered_input_does_not_warn(chains: dict[str, tuple]) -> None:
    ch, records = chains["raw"]
    n, k, mean_pep = _stratum0(ch)
    assert n > 200 and k >= 0.5 * n * mean_pep, (n, k, mean_pep)  # the design: k ~ its expectation
    assert not _warned(records), (n, k, mean_pep)


def test_prefiltered_input_warns_input_prefiltered(chains: dict[str, tuple]) -> None:
    ch, records = chains["prefiltered"]
    n, k, mean_pep = _stratum0(ch)
    assert n > 200 and mean_pep == pytest.approx(0.05), (n, mean_pep)
    assert k <= 0.2 * n * mean_pep, (n, k, mean_pep)  # far below the expected unmapped count
    assert _warned(records), (n, k, mean_pep)
