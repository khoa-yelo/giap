"""Inputs, oracles and frame builders for the M5-09 context tables
(docs/specs/M5-09-context-tables-step.md).

The fixture is the M4-18 integration fixture (``detection_summary_fixture.sim_fixture``: M4-07
layout + M4-08 expression + perfect C9/C11-C14 evidence) with two M5 additions:

* an **oriC** (C15 ``features``): ``LayoutConfig.oric_frac`` places a 232-bp ``oric`` feature and
  draws nothing, so the layout is byte-identical to the one ``sim_fixture`` wrote and
  ``write_contract_tables`` can simply rewrite C1-C8 + C15 over it;
* a **planted island** (``ExpressionConfig.islands``), which only changes the gene table this
  module builds (as M4-18's integration test does, through ``region_tests_data.tables_from``).

M4-20's ``giapsim.fixtures.build_m4_fixture`` (the builder the spec names) is not merged; when it
is, :func:`build` should become ``m4_chain.build(work, seed=, size=, layout={"oric_frac": ...},
evidence=...)`` plus the same gene-table override. Nothing else in these tests depends on how the
work directory is filled.

Oracles call the **merged** producers ``giap.context.autocorrelation`` (M5-05) and
``giap.context.islands`` (M5-06) and the **unmerged** ``giap.context.oriter`` (M5-08, PR #132)
directly on the step's own random streams (§6.6: ``giap.rng.for_step(<stream>, cfg.seed)``), so a
test can compare what the step wrote with what the producers produce without re-deriving any
statistic.
"""

from __future__ import annotations

import math
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import detection_summary_fixture as fx
import region_tests_data as rd

REPO = Path(__file__).resolve().parents[2]
TABLE_DIR = REPO / "src" / "giap" / "schemas" / "tables"

#: §9.3 ``tables/autocorrelation.tsv`` columns (= ``giap.context.autocorrelation``
#: ``AUTOCORR_COLUMNS``; a literal copy, so a silent reorder fails a test)
AUTOCORR_COLUMNS: tuple[str, ...] = (
    "contig", "statistic", "lag", "null_model", "value", "null_mean", "null_lo", "null_hi",
    "stat", "effect", "ci_low", "ci_high", "p_value", "q_value", "method", "n_perm", "flags",
)  # fmt: skip
#: §9.3 ``tables/islands.tsv`` columns (= ``giap.context.islands.ISLAND_COLUMNS``)
ISLAND_COLUMNS: tuple[str, ...] = (
    "island_id", "direction", "contig", "start", "end", "tu_first", "tu_last", "n_tus",
    "n_genes", "genes", "mean_abundance_log2", "detected_frac", "scan_s", "block_len",
    "stat", "effect", "ci_low", "ci_high", "p_value", "q_value", "method", "null_model",
    "n_perm", "flags",
)  # fmt: skip
#: §9.3 genes.tsv context group, as ``work/context/gene_context.parquet`` holds it
GENE_CONTEXT_COLUMNS: tuple[str, ...] = (
    "locus_tag", "island_id", "replichore", "replichore_position",
)  # fmt: skip

#: integer columns of each table (Frictionless ``integer``; ``lag`` carries NA, edge case 2)
AUTOCORR_INTS: tuple[str, ...] = ("lag", "n_perm")
ISLAND_INTS: tuple[str, ...] = ("start", "end", "n_tus", "n_genes", "block_len", "n_perm")
AUTOCORR_FLOATS: tuple[str, ...] = ("value", "null_mean", "null_lo", "null_hi", "effect",
                                    "ci_low", "ci_high", "p_value", "q_value")  # fmt: skip
ISLAND_FLOATS: tuple[str, ...] = ("mean_abundance_log2", "detected_frac", "scan_s", "effect",
                                  "ci_low", "ci_high", "p_value", "q_value")  # fmt: skip
#: list-valued cell (``;``-joined in the TSV)
ISLAND_LISTS: tuple[str, ...] = ("genes",)

#: §6.6 random streams of the M5 parts of ``context_stats`` (the producers' docstrings)
STREAMS: Mapping[str, str] = {
    "autocorrelation": "context_stats.autocorrelation",
    "islands": "context_stats.islands",
    "oriter": "context_stats.oriter",
}
#: §9.4 "Region and context" flags that this issue's rows can carry
DESCRIPTIVE = "descriptive"
NO_ORIC = "no_oric"
NOT_ESTIMABLE = "not_estimable"
FEW_VALUES = "few_values"
CI_UNESTIMATED = "ci_unestimated"
CI_FEW_UNITS = "ci_few_units"
#: ADR-0021 Amendment 1 / M5-06 decision 48: never emitted by M5
UNCALIBRATED = "uncalibrated"

SEED = 11
N_GENES = 260
ORIC_FRAC = 0.3


# --------------------------------------------------------------------------- configuration


def config(
    *,
    seed: int = SEED,
    permutations: int = 39,
    bootstrap: int = 39,
    alpha_context: float = 0.05,
    moran_lags: Sequence[int] = (1, 2),
    min_tus: int = 2,
    max_tus: int = 20,
    span_frac: float = 0.10,
    block_len: int | str = "auto",
) -> Any:
    """``GiapConfig`` with the §5.2 knobs the M5 parts of ``context_stats`` read."""
    from giap.config.model import GiapConfig

    return GiapConfig.model_validate(
        {
            "seed": seed,
            "stats": {
                "permutations": permutations,
                "bootstrap": bootstrap,
                "alpha": {"context": alpha_context},
            },
            "context": {
                "moran_lags": list(moran_lags),
                "islands": {
                    "min_tus": min_tus,
                    "max_tus": max_tus,
                    "detrend_span_frac": span_frac,
                    "block_len": block_len,
                },
            },
        }
    )


# --------------------------------------------------------------------------- the work directory


@dataclass(frozen=True)
class Fixture:
    """A work directory ready for ``run_context_stats`` and the frames behind it."""

    work: Path
    genes: pd.DataFrame  # the gene table the step reads (work/genetable/genes.parquet)
    replicons: pd.DataFrame  # C1
    tus: pd.DataFrame  # C5
    features: pd.DataFrame | None  # C15 (None when features.parquet was removed)


def _layout(seed: int, n_genes: int, oric_frac: float | None) -> Any:
    """The ``sim_fixture`` layout, optionally with an oriC (placing one draws nothing)."""
    from giapsim.layout import LayoutConfig, simulate_layout

    return simulate_layout(
        LayoutConfig(
            n_genes_chromosome=n_genes,
            plasmids=(25,),
            linear_contigs=(15,),
            oric_frac=oric_frac,
        ),
        np.random.default_rng(seed),
    )


def _beta(coefficient: float, mu_log2: float = 20.0) -> dict[str, float]:
    """M4-08 detection coefficients with an abundance term, centred on ``mu_log2`` so the overall
    detected share stays where the fixture had it and only the planted island moves: a gene in a
    +12 log2 island is then detected with probability near 1, which is what makes an island
    reachable at all (M5-06 decision 50: left censoring caps the high direction too)."""
    return {
        "intercept": -8.65 - coefficient * mu_log2,
        "log_length": 1.2,
        "log1p_n_obs": 0.8,
        "gravy": -0.6,
        "log2_abundance": coefficient,
    }


def _expression(
    layout: Any,
    seed: int,
    islands: Sequence[Any],
    oriter_log2: float,
    detection_beta: float,
) -> Any:
    from giapsim.expression import ExpressionConfig, simulate_expression

    cfg = ExpressionConfig(
        detection="abundance",
        ar1_rho=0.5,
        beta=_beta(detection_beta),
        islands=tuple(islands),
        oriter_log2=oriter_log2,
    )
    return simulate_expression(layout, cfg, np.random.default_rng(seed + 1))


def island_spec(first_tu: int, n_tus: int, effect: float, replicon: str = "chr1") -> Any:
    from giapsim.expression import IslandSpec

    return IslandSpec(replicon=replicon, first_tu=first_tu, n_tus=n_tus, effect=effect)


def build(
    work: Path,
    *,
    seed: int = SEED,
    n_genes: int = N_GENES,
    oric_frac: float | None = ORIC_FRAC,
    islands: Sequence[Any] = (),
    oriter_log2: float = 0.0,
    detection_beta: float = 0.5,
    features: bool = True,
) -> Fixture:
    """Write a work directory: C1-C15 contract tables and ``genetable/genes.parquet``.

    ``features=False`` deletes ``features.parquet`` (example 5); ``oric_frac=None`` keeps it but
    header-only (no oriC call).
    """
    from giap.io.work_tables import default_paths, load_m4_inputs
    from giapsim.contract_tables import write_contract_tables

    fx.sim_fixture(work, seed=seed, n_genes=n_genes)
    layout = _layout(seed, n_genes, oric_frac)
    write_contract_tables(layout, work)  # same C1-C8 (nothing drawn), C15 with the oriC
    paths = default_paths(work)
    inputs, issues = load_m4_inputs(paths)
    assert not [i for i in issues if i.severity == "error"], issues[:5]

    expr = _expression(layout, seed, islands, oriter_log2, detection_beta)
    tables = rd.tables_from(layout, expr, "true", genes_c2=inputs.genes, layers=inputs.layers)
    genes = tables.genes.copy()
    c2 = inputs.genes.set_index("locus_tag")
    tags = genes["locus_tag"].astype(str)
    genes["hypothetical"] = tags.map(c2["hypothetical"].astype(bool)).to_numpy(bool)
    genes["cog_category"] = fx.cog_categories(inputs.genes, inputs.xrefs).to_numpy()
    genes["abundance_imputed"] = False
    (work / "genetable").mkdir(parents=True, exist_ok=True)
    genes.to_parquet(work / "genetable" / "genes.parquet", index=False)

    feature_path = work / "features.parquet"
    if not features:
        feature_path.unlink()
    feats = pd.read_parquet(feature_path) if feature_path.exists() else None
    return Fixture(work, genes, inputs.replicons, inputs.tus, feats)


# --------------------------------------------------------------------------- oracles


def oracle_autocorrelation(f: Fixture, cfg: Any) -> pd.DataFrame:
    """M5-05 rows plus M5-08 rows, each on its own stream, in the spec's "Row assembly" order."""
    from giap.context.autocorrelation import autocorrelation_rows
    from giap.context.oriter import oriter_rows
    from giap.rng import for_step

    ac = autocorrelation_rows(
        f.genes, f.replicons, cfg, for_step(STREAMS["autocorrelation"], cfg.seed)
    )
    ot = oriter_rows(
        f.genes, f.tus, f.replicons, f.features, cfg, for_step(STREAMS["oriter"], cfg.seed)
    )
    return sort_autocorrelation(pd.concat([ac, ot], ignore_index=True), f.replicons)


def oracle_islands(f: Fixture, cfg: Any) -> tuple[pd.DataFrame, Any]:
    """The M5-06 island frame and its ``IslandScan`` on the step's island stream."""
    import islands_data as isd

    from giap.rng import for_step

    series = isd.build_series(f.genes, f.tus, f.replicons, cfg)
    rng = for_step(STREAMS["islands"], cfg.seed)
    frame, scan = isd.step(f.genes, series, cfg, rng)
    return frame, scan


def oracle_gene_context(f: Fixture, islands: pd.DataFrame) -> pd.DataFrame:
    """``gene_context``: the M5-06 gene→island map and the M5-08 replichore columns, left-joined
    onto the gene table (one row per gene, NA where a gene has neither)."""
    from giap.context.islands import gene_island_map
    from giap.context.oriter import gene_replichores

    tags = f.genes["locus_tag"].astype(str)
    isl = gene_island_map(islands).set_index("locus_tag")["island_id"]
    rep = gene_replichores(f.genes, f.replicons, f.features).set_index("locus_tag")
    return pd.DataFrame(
        {
            "locus_tag": tags.to_numpy(),
            "island_id": tags.map(isl).to_numpy(object),
            "replichore": tags.map(rep["replichore"]).to_numpy(object),
            "replichore_position": tags.map(rep["replichore_position"]).to_numpy(float),
        },
        columns=list(GENE_CONTEXT_COLUMNS),
    )


#: §9.3 ``statistic`` enum order (the spec's "Row assembly" sort key)
STATISTIC_ORDER: tuple[str, ...] = ("morans_i", "join_count_bb", "oriter_slope")
#: §9.3 ``null_model`` enum order
NULL_ORDER: tuple[str, ...] = ("gene_label", "tu_block", "bootstrap")


def sort_autocorrelation(df: pd.DataFrame, replicons: pd.DataFrame) -> pd.DataFrame:
    """Spec "Row assembly": contig (C1 order), ``statistic`` and ``null_model`` in §9.3 enum
    order, ``lag`` (NA last), ``stat``."""
    order = {str(c): i for i, c in enumerate(replicons["contig"])}
    key = pd.DataFrame(
        {
            "c": [order[str(v)] for v in df["contig"]],
            "s": [STATISTIC_ORDER.index(str(v)) for v in df["statistic"]],
            "n": [NULL_ORDER.index(str(v)) for v in df["null_model"]],
            "l": [math.inf if pd.isna(v) else float(v) for v in df["lag"]],
            "t": [str(v) for v in df["stat"]],
        }
    )
    out = df.loc[key.sort_values(["c", "s", "n", "l", "t"], kind="mergesort").index]
    return out.reset_index(drop=True)


# --------------------------------------------------------------------------- reading the output


def read_table(path: Path) -> pd.DataFrame:
    """A written TSV as strings (``NA`` → None), for the schema-rule checks."""
    from giap.io.tsv import read_tsv

    header, rows = read_tsv(path)
    return pd.DataFrame(rows, columns=header, dtype=object)


def flag_set(cell: Any) -> set[str]:
    """Flags of a table cell, whether it is a tuple (frame) or a ``;``-joined string (TSV)."""
    if cell is None or (isinstance(cell, float) and math.isnan(cell)):
        return set()
    if isinstance(cell, str):
        return {f for f in cell.split(";") if f and f != "NA"}
    return {str(f) for f in cell}


def is_na(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, str):
        return value == "NA"
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def row_key(row: Mapping[str, Any], columns: Sequence[str]) -> tuple[str, ...]:
    """Identity of a table row: its key columns as strings ("NA" for missing)."""
    return tuple("NA" if is_na(row[c]) else str(row[c]) for c in columns)
