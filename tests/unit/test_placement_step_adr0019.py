"""ADR-0019 decisions 3 and 4 through the registered ``placement`` step (M3-25 decisions 3, 4):
the ``stratum0_classes`` output and the run warning ``input_prefiltered`` (an input already
filtered to genome matches cannot be calibrated by any genome-based π̂; ADR-0019 Evidence,
"Limitation").

M3-25 decision 3: ``placement`` writes the class table as a new output ``stratum0_classes``
(``mapping/stratum0_classes.parquet`` in the run layout, the columns of ``stratum0_pi_classes``),
following M4-18's ``detectability_summary`` pattern; ``mapping_fdr.tsv`` columns are unchanged.

Convention (as for ``c_assumed``, tests/unit/test_mapping_steps.py): a run warning is a WARNING
record on a ``gippa.steps.*`` logger whose message names the flag. Inputs: a 50-kb gippasim genome,
database input (``search_db = other``, one run, so the row-level mean PEP and the peptide-level one
agree). Control: the stressed arm (``inject_short`` = 200 short wrong peptides absent from the
genome, PEP = the true wrong share): k (stratum-0 database peptides without an exact target hit)
is near n * mean PEP * (1 - r), no warning is due. Pre-filtered (gippasim's wrong identifications
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
from helpers.mapping_steps import fresh_steps, key_like
from helpers.placement_sims import Chain, run_chain

SEED = 20261022


class _Sink(logging.Handler):
    def __init__(self) -> None:
        super().__init__(logging.INFO)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


def _chain(tmp: Path, prefilter: bool) -> tuple[Chain, list[logging.LogRecord]]:
    from gippasim.proteomics import ProteomicsConfig

    kw = {"prefilter": 0.05} if prefilter else {"inject_short": 200}

    sink = _Sink()
    log = logging.getLogger("gippa.steps")
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
        if r.name.startswith("gippa.steps")
    )


def _stratum0(ch: Chain) -> tuple[int, int, float]:
    """(n, k, mean PEP) of the database stratum-0 population, from the step's own outputs."""
    from gippa.steps.mapping import select_queries

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


# ----------------------------------------------------------- decision 3: stratum0_classes output

STRATUM0_COLUMNS = ("length_class", "n", "k", "r", "upper", "pi_wrong")


def test_placement_declares_the_stratum0_classes_output() -> None:
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        spec = reg.get("placement")
    assert key_like(spec.outputs, "stratum0_classes")
    assert key_like(spec.outputs, "mapping_fdr")  # unchanged alongside it


def test_placement_writes_the_class_table(chains: dict[str, tuple]) -> None:
    """The output holds ``stratum0_pi_classes`` of the step's own stratum-0 population: one row
    per ``mapping.length_classes`` class, in config order."""
    from gippa.mapping.placement_fdr import length_class_labels, stratum0_pi_classes

    ch, _records = chains["raw"]
    tab = pd.read_parquet(ch.avail[key_like(list(ch.avail), "stratum0_classes")])
    labels = length_class_labels([list(c) for c in ch.cfg.mapping.length_classes])
    assert list(tab.columns) == list(STRATUM0_COLUMNS)
    assert tab["length_class"].astype(str).tolist() == labels
    assert int(tab["n"].sum()) > 200
    want = stratum0_pi_classes(_population(ch), [list(c) for c in ch.cfg.mapping.length_classes])
    pd.testing.assert_frame_equal(tab.reset_index(drop=True), want.reset_index(drop=True))


def _population(ch: Chain) -> pd.DataFrame:
    """P: the database queries except ``same_db`` exact peptides and contaminants (decision 7),
    at their best costs, as ``placement`` sees them."""
    from helpers.placement_sims import hits_of

    from gippa.mapping.chance import residue_freqs
    from gippa.mapping.placement_fdr import best_costs
    from gippa.mapping.sixframe import load_index
    from gippa.steps.mapping import select_queries

    q = select_queries(pd.read_parquet(ch.avail["peptides"]), ch.cfg)
    idx, _decoys = load_index(ch.avail[key_like(list(ch.avail), "sixframe")], table=11)
    q = q.copy()
    q["pi_denovo"] = np.nan
    costs = best_costs(hits_of(ch.avail), q, residue_freqs(idx))
    keep = costs["source"].astype(str) == "database"
    keep &= ~((costs["search_db"].astype(str) == "same") & (costs["m_target"] == 0))
    if "is_contaminant" in costs:
        keep &= ~costs["is_contaminant"].fillna(False).astype(bool)
    return costs[keep].reset_index(drop=True)
