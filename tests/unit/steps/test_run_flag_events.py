"""MVP-03a: every D-E "add an event" site emits a §6.7 event carrying ``data.flag``
(docs/specs/MVP-03-manifest-finalize.md D-E and acceptance item
``tests/unit/steps/test_run_flag_events.py``; GIAP_SPEC.md §6.7, §9.4 Run row).

The D-E event shape is the existing one: the step logs through ``ctx.log`` with
``extra={"event": <name>, "data": {"flag": <token>, ...}}``; ``giap.io.events.EventsHandler``
turns such a record into an ``events.jsonl`` line. The event *name* is not pinned
(``<step>.run_flag`` or ``<step>.layer_missing``); a record counts when its ``event`` is a string
and its ``data`` maps ``flag`` to the token. Sites (step: trigger -> token):

* ``region_layer``: no ``genomad`` / ``antismash`` input -> ``layer_missing:genomad`` /
  ``layer_missing:antismash``;
* ``operons``: no ``promoter_windows`` / ``terminators`` / ``tx_rates`` ->
  ``layer_missing:promoters`` / ``layer_missing:terminators`` / ``layer_missing:tx_rates``;
* ``operons``: pass 1 (always) -> ``no_detection_limit_model``;
* ``context_stats``: no ``features`` input (M5 outputs wired) -> ``layer_missing:features``;
* ``placement``: de novo peptides, no reference peptides -> ``c_assumed``;
* ``operon_pairs``: no ``uniop_pred`` input -> ``operon_source_distance_model``;
* ``codon_usage``: no ribosomal-protein CAI reference -> ``cai_reference_missing``.

Inputs: ``helpers.layers_chain.build_chain`` (merged producers on giapsim SIM100K), the
``context_tables_data`` fixture without features (M5-09 example 5) and the M3-17 mapping chain on
de novo peptides (``test_mapping_steps.py`` example 6). The raw simulator Bakta JSON has no
``ribosomal protein`` product, so ``parse_bakta`` on it gives a gene set without a CAI reference.
"""

from __future__ import annotations

import logging
import math
import sys
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.layers_chain import PARSE_OUT, Chain, build_chain, registry

from giap.config import GiapConfig

HELPERS = Path(__file__).resolve().parents[2] / "helpers"


class _Recorder(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


@contextmanager
def recording() -> Iterator[list[logging.LogRecord]]:
    """Every record of the ``giap.steps`` loggers (any level) while the block runs."""
    rec = _Recorder()
    log = logging.getLogger("giap.steps")
    old = log.level
    log.addHandler(rec)
    log.setLevel(logging.DEBUG)
    try:
        yield rec.records
    finally:
        log.removeHandler(rec)
        log.setLevel(old)


def flag_events(records: list[logging.LogRecord], token: str) -> list[logging.LogRecord]:
    out = []
    for r in records:
        event, data = getattr(r, "event", None), getattr(r, "data", None)
        if isinstance(event, str) and isinstance(data, Mapping) and data.get("flag") == token:
            out.append(r)
    return out


def assert_flag_event(records: list[logging.LogRecord], token: str, step: str) -> None:
    hits = flag_events(records, token)
    named = [r.name for r in records if token in r.getMessage()]
    assert hits, f"no event with data.flag == {token!r}; log-only records: {named}"
    assert any(r.name == f"giap.steps.{step}" for r in hits), [r.name for r in hits]


# ----------------------------------------------------------------------------- fixtures


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(Path(tmp_path_factory.mktemp("run_flag_chain")))


def _work(chain: Chain) -> Path:
    return chain.inputs["genes"].parent


# ----------------------------------------------------------------------------- region_layer


@pytest.mark.parametrize("layer", ["genomad", "antismash"])
def test_region_layer_layer_missing_event(chain: Chain, tmp_path: Path, layer: str) -> None:
    ins = {k: chain.inputs[k] for k in ("genes", "features", "replicons")}
    outs = {"regions": tmp_path / "regions.parquet", "region_genes": tmp_path / "rg.parquet"}
    with registry() as reg, recording() as recs:
        reg.run_step("region_layer", GiapConfig(seed=1), ins, outs)
    assert_flag_event(recs, f"layer_missing:{layer}", "region_layer")


# ----------------------------------------------------------------------------- operon_pairs


def test_operon_pairs_distance_model_event(chain: Chain, tmp_path: Path) -> None:
    w = _work(chain)
    ins = {"genes": chain.inputs["genes"], "replicons": chain.inputs["replicons"],
           "genome": chain.inputs["genome"], "uniop_cds": w / "uniop_cds.gff",
           "uniop_map": w / "uniop_map.tsv"}  # fmt: skip
    outs = {"pairs_raw": tmp_path / "pairs_raw.parquet", "tus": tmp_path / "tus.parquet",
            "tu_upstream": tmp_path / "tu_upstream",
            "operon_meta": tmp_path / "operon_meta.json"}  # fmt: skip
    with registry() as reg, recording() as recs:
        reg.run_step("operon_pairs", GiapConfig(seed=1), ins, outs)
    assert_flag_event(recs, "operon_source_distance_model", "operon_pairs")


# ----------------------------------------------------------------------------- codon_usage


def test_codon_usage_cai_reference_missing_event(chain: Chain, tmp_path: Path) -> None:
    w = tmp_path / "raw"
    with registry() as reg:
        reg.run_step("parse_bakta", GiapConfig(seed=1),
                     {"bakta_json": Path(chain.sim.bakta_json), "genome": chain.inputs["genome"]},
                     {k: w / v for k, v in PARSE_OUT.items()})  # fmt: skip
        with recording() as recs:
            reg.run_step("codon_usage", GiapConfig(seed=1),
                         {"genes": w / "genes.parquet", "trna_counts": w / "trna_counts.tsv"},
                         {"codon": tmp_path / "codon.parquet"})  # fmt: skip
    assert_flag_event(recs, "cai_reference_missing", "codon_usage")


def test_codon_usage_no_event_with_a_reference(chain: Chain, tmp_path: Path) -> None:
    """The chain's Bakta JSON names ribosomal proteins: no ``cai_reference_missing`` event."""
    w = _work(chain)
    with registry() as reg, recording() as recs:
        reg.run_step("codon_usage", GiapConfig(seed=1),
                     {"genes": w / "genes.parquet", "trna_counts": w / "trna_counts.tsv"},
                     {"codon": tmp_path / "codon.parquet"})  # fmt: skip
    assert not flag_events(recs, "cai_reference_missing")


# ----------------------------------------------------------------------------- operons


def _gene_table(chain: Chain, path: Path) -> Path:
    """A §9.3 gene table over the C2 genes (as tests/unit/steps/test_operons_step.py builds it)."""
    from giap.genetable.assemble import GENES_COLUMNS, GENES_TYPES
    from giap.io.genes_tsv import write_genes_parquet

    c2 = chain.frame("genes")
    n = len(c2)
    default = {"list": lambda: [], "bool": lambda: False, "int": lambda: None,
               "float": lambda: math.nan, "str": lambda: None}  # fmt: skip
    df = pd.DataFrame({c: [default[GENES_TYPES[c]]() for _ in range(n)] for c in GENES_COLUMNS})
    for c in ("locus_tag", "gene", "product", "feature_type", "contig", "start", "end", "strand",
              "partial", "crosses_origin", "hypothetical"):  # fmt: skip
        df[c] = c2[c].to_numpy()
    det = np.arange(n) % 4 != 3
    df["detected"] = det
    df["detection_status"] = np.where(det, "detected", "not_detected")
    df["p_detect"] = np.linspace(0.2, 0.9, n)
    df["abundance_imputed"] = False
    path.parent.mkdir(parents=True, exist_ok=True)
    write_genes_parquet(path, df)
    return path


@pytest.fixture(scope="module")
def operons_bare(chain: Chain, tmp_path_factory: pytest.TempPathFactory) -> list[Any]:
    """The ``operons`` step without its three optional layers."""
    out = tmp_path_factory.mktemp("operons_bare")
    ins = {k: chain.inputs[k] for k in ("genes", "replicons", "tus", "pairs_raw", "features")}
    ins["gene_table"] = _gene_table(chain, out / "genes.parquet")
    outs = {"pairs": out / "pairs.tsv", "operons": out / "operons.tsv",
            "operon_flags": out / "operon_flags.tsv"}  # fmt: skip
    with registry() as reg, recording() as recs:
        reg.run_step("operons", GiapConfig(seed=1), ins, outs)
    return recs


@pytest.mark.parametrize("layer", ["promoters", "terminators", "tx_rates"])
def test_operons_layer_missing_event(operons_bare: list[Any], layer: str) -> None:
    assert_flag_event(operons_bare, f"layer_missing:{layer}", "operons")


def test_operons_no_detection_limit_model_event(operons_bare: list[Any]) -> None:
    assert_flag_event(operons_bare, "no_detection_limit_model", "operons")


def test_operons_missed_operons_event_unchanged(operons_bare: list[Any]) -> None:
    """The M6-08 event stays (03c reads it); it is not a flag event."""
    names = [getattr(r, "event", None) for r in operons_bare]
    assert "operons.missed_operons" in names


# ----------------------------------------------------------------------------- context_stats


def test_context_stats_layer_missing_features_event(tmp_path: Path) -> None:
    sys.path.insert(0, str(HELPERS))
    try:
        import context_tables_data as cd
    finally:
        sys.path.remove(str(HELPERS))
    from giap.steps.context_stats import run_context_stats

    work = tmp_path / "work"
    cd.build(work, features=False)
    with recording() as recs:
        run_context_stats(work, cd.config(), None)
    assert_flag_event(recs, "layer_missing:features", "context_stats")


# ----------------------------------------------------------------------------- placement


def test_placement_c_assumed_event(tmp_path: Path) -> None:
    from helpers.mapping_steps import SimRun, drive, roots_of, simulate_with_proteomics

    from giap.config import resolve_config
    from giapsim.proteomics import ProteomicsConfig

    seed = 11
    run = simulate_with_proteomics(
        tmp_path / "in", seed, proteomics=ProteomicsConfig(pi_wrong=0.05, search_db="other")
    )
    df = pd.read_parquet(run.peptides)
    df["source"] = pd.Series(["denovo"] * len(df), dtype="string")
    df["q_value"] = np.nan
    df["pep"] = np.nan
    peps = tmp_path / "peptides.denovo.parquet"
    df.to_parquet(peps, index=False)
    run = SimRun(run.res, run.genome, peps, run.proteins_long, run.config_yaml)
    cfg = resolve_config(run.config_yaml, [f"seed={seed}"])
    with registry() as reg, recording() as recs:
        drive(reg, cfg, tmp_path / "work", roots_of(run))
    assert_flag_event(recs, "c_assumed", "placement")
