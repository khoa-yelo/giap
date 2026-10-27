"""Step ``determinants`` (GIAP_SPEC.md §6.2, §7.8; docs/specs/M8-10-determinants-step-and-run.md):
gene table + C1/C2/C5/C6 (+ optional C7/C8 and M5-09's gene context) →
``tables/determinants_{features,model,attribution}.tsv``.

Order: features (M8-07) → univariate, repeated grouped CV and drop-one-block ΔR² (M8-08) →
LightGBM SHAP and the TU-size-preserving null (M8-09) → writers. Streams from ``config.seed``:
``determinants.bootstrap`` (univariate CIs, ``stats.bootstrap`` draws), ``determinants.cv`` (the
outer folds) and ``determinants.null`` (inside ``null_p``, B = ``determinants.permutations``).
``threads`` is ``config.threads`` capped by the CPUs this process may use; it never changes the
tables (fixed per-fit ``num_threads``, the null collected in b order).

Without abundance (every ``abundance_log2`` NA) the features table and the univariate AUROC rows
are written, the attribution table is header-only and the run flag is
``layer_missing:abundance``. Below the 100-row minimum each model gets one ``too_few_rows`` row.
The §7.8.2 caveats and the features dropped from the models go to the step's events
(``determinants.caveats``, ``determinants.features_dropped``) for the M9 report.
"""

from __future__ import annotations

import logging
import math
import os
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq  # type: ignore[import-untyped]

from giap.annotate.operons import Tu
from giap.config.model import GiapConfig
from giap.determinants import CAVEATS
from giap.determinants.attribution import (
    ATTRIBUTION_COLUMNS,
    attribution_table,
    shap_blocks,
    shap_features,
)
from giap.determinants.features import FEATURE_BLOCKS, feature_frame, model_matrix
from giap.determinants.models import (
    MODELS,
    TABLE_COLUMNS,
    CvResult,
    cross_validate,
    enough_data,
    model_rows,
    model_rows_table,
    tu_folds,
    univariate,
)
from giap.determinants.null import NullResult, null_p
from giap.determinants.support import null_rows, pd_rows, shap_rows
from giap.io.determinants_tables import (
    write_determinants_attribution,
    write_determinants_features,
    write_determinants_model,
)
from giap.io.plot_support import write_support
from giap.model.coords import Interval
from giap.rng import for_step
from giap.steps.codon import _read_genes
from giap.steps.common import one_path as _path
from giap.steps.common import read_replicons as _replicons
from giap.steps.registry import StepContext, step

__all__ = ["determinants_step"]

#: ``determinants_model.tsv`` row order (spec step 3): model, metric, then ``feature_or_block``
MODEL_ORDER: tuple[str, ...] = ("lightgbm", "elastic_net", "univariate")
METRIC_ORDER: tuple[str, ...] = (
    "cv_r2_mean", "cv_r2_sd", "null_p", "delta_r2", "spearman_rho", "auroc",
)  # fmt: skip
NO_ABUNDANCE_FLAG = "layer_missing:abundance"
#: M9-17: optional directory output; the SHAP, PD and null files go to its ``support/``
SUPPORT = "determinants_support"
#: the null p's effect has no interval (lead decision 4, ADR-0031 pattern)
CI_UNESTIMATED = "ci_unestimated"


def _threads(cfg: GiapConfig) -> int:
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()
    return max(1, min(cfg.threads, cpus or 1))


def _tus(path: Path) -> list[Tu]:
    """C5 rows -> ``Tu``. The contract table has no upstream window and ``feature_frame`` does
    not use one, so ``upstream`` is a 1-bp stand-in at the TU start."""
    out = []
    for r in pq.read_table(path, columns=["tu_id", "contig", "start", "end", "strand",
                                          "genes"]).to_pylist():  # fmt: skip
        iv = Interval(r["contig"], r["start"], r["end"], r["strand"])
        up = Interval(r["contig"], r["start"], r["start"] + 1, r["strand"])
        out.append(Tu(r["tu_id"], r["contig"], iv, r["strand"], tuple(r["genes"]), up))
    return out


def _regions(inputs: Mapping[str, Path]) -> pd.DataFrame | None:
    """C8 with the C7 ``region_type`` joined on (``feature_frame``'s ``regions``), or None."""
    if "regions" not in inputs or "region_genes" not in inputs:
        return None
    c7 = pd.read_parquet(inputs["regions"])[["region_id", "region_type"]]
    c8 = pd.read_parquet(inputs["region_genes"])
    return c8.merge(c7.drop_duplicates("region_id"), on="region_id", how="left",
                    validate="many_to_one")  # fmt: skip


def _features(inputs: Mapping[str, Path]) -> pd.DataFrame:
    ctx = inputs.get("gene_context")
    return feature_frame(
        pd.read_parquet(inputs["gene_table"]),
        _read_genes(inputs["genes"]),
        _tus(inputs["tus"]),
        pd.read_parquet(inputs["layers"]),
        _regions(inputs),
        None if ctx is None else pd.read_parquet(ctx),
        _replicons(inputs["replicons"]),
    )


def _add_flag(flags: object, flag: str) -> str:
    have = [f for f in str(flags).split(";") if f] if isinstance(flags, str) else []
    return ";".join(have if flag in have else [*have, flag])


def _null_rows(nulls: Mapping[str, NullResult], n: int) -> pd.DataFrame:
    """One ``null_p`` row per model; ``n_perm`` counts the non-NaN null values (decision 4)."""
    rows = [
        {"model": m, "feature_or_block": "all", "metric": "null_p", "value": r.p,
         "ci_low": math.nan, "ci_high": math.nan, "n": n,
         "n_perm": int(np.count_nonzero(~np.isnan(r.null_r2))), "flags": r.flags}
        for m, r in nulls.items()
    ]  # fmt: skip
    df = pd.DataFrame(rows, columns=list(TABLE_COLUMNS))
    df[["value", "ci_low", "ci_high"]] = df[["value", "ci_low", "ci_high"]].astype(float)
    df[["n", "n_perm"]] = df[["n", "n_perm"]].astype("Int64")
    return df


def _model_table(
    cv: Mapping[str, CvResult], uni: pd.DataFrame, nulls: Mapping[str, NullResult], n: int
) -> pd.DataFrame:
    """Spec step 3: CV, ΔR² and univariate rows plus the null rows, sorted; every
    ``cv_r2_mean`` row carries ``ci_unestimated``."""
    tbl = model_rows_table(cv, uni)
    if nulls:
        tbl = pd.concat([tbl, _null_rows(nulls, n)], ignore_index=True)
    mean = tbl["metric"] == "cv_r2_mean"
    tbl.loc[mean, "flags"] = [_add_flag(f, CI_UNESTIMATED) for f in tbl.loc[mean, "flags"]]
    return _sorted(tbl)


def _sorted(tbl: pd.DataFrame) -> pd.DataFrame:
    key = pd.DataFrame({
        "m": tbl["model"].map(MODEL_ORDER.index), "k": tbl["metric"].map(METRIC_ORDER.index),
        "f": tbl["feature_or_block"].astype(str),
    })  # fmt: skip
    order = key.sort_values(["m", "k", "f"], kind="stable").index
    return tbl.loc[order].reset_index(drop=True)


def _events(features: pd.DataFrame, log: logging.Logger) -> None:
    caveats = {"event": "determinants.caveats", "data": {"caveats": list(CAVEATS)}}
    log.info("determinants caveats (§7.8.2)", extra=caveats)
    dropped = [f for fs in FEATURE_BLOCKS.values() for f in fs if bool(features[f].isna().all())]
    event = {"event": "determinants.features_dropped", "data": {"features": dropped}}
    log.info("determinants: features NA in every row, dropped from the models: %s",
             ", ".join(dropped) or "none", extra=event)  # fmt: skip


def _determinants(
    inputs: Mapping[str, Path], outputs: Mapping[str, Path], cfg: GiapConfig, log: logging.Logger
) -> None:
    if SUPPORT in outputs:
        outputs[SUPPORT].mkdir(parents=True)
    features = _features(inputs)
    write_determinants_features(outputs["determinants_features"], features)
    _events(features, log)
    seed = cfg.seed
    uni = univariate(features, cfg.stats.bootstrap, for_step("determinants.bootstrap", seed))
    if bool(features["abundance_log2"].isna().all()):
        log.warning("%s: no abundance; determinants_model.tsv has the AUROC rows only and "
                    "determinants_attribution.tsv is header-only", NO_ABUNDANCE_FLAG)  # fmt: skip
        write_determinants_model(outputs["determinants_model"],
                                 _sorted(uni[uni["metric"] == "auroc"]))  # fmt: skip
        write_determinants_attribution(outputs["determinants_attribution"],
                                       pd.DataFrame(columns=list(ATTRIBUTION_COLUMNS)))  # fmt: skip
        return

    rows = model_rows(features)
    X_all, blocks = model_matrix(features)
    X = X_all[rows]
    y = pd.to_numeric(features.loc[rows, "abundance_log2"]).to_numpy(float)
    tu = features.loc[rows, "tu_id"].astype(str).tolist()
    cv: dict[str, CvResult] = {}
    nulls: dict[str, NullResult] = {}
    shap = None
    if enough_data(tu):
        d, threads = cfg.determinants, _threads(cfg)
        folds = tu_folds(tu, d.outer_folds, d.repeats, for_step("determinants.cv", seed))
        cv = cross_validate(X, y, folds, blocks, d, seed, threads, tu=tu)
        observed: dict[str, float] = {m: float(np.mean(cv[m].r2)) for m in MODELS}
        nulls = null_p(X, y, tu, folds[0], observed, d, d.permutations, seed, threads)
        flag = features[["ribosomal_protein", "detected"]].astype("boolean").fillna(False)
        reference = flag.all(axis=1).to_numpy(bool) & ~rows  # the CAI reference genes (Q3)
        shap = shap_blocks(cv["lightgbm"], X, blocks, extra=X_all[reference])
        shap.index = pd.Index(features.loc[shap.index, "locus_tag"].astype(str))
    table = _model_table(cv, uni, nulls, len(y))
    write_determinants_model(outputs["determinants_model"], table)
    att = attribution_table(features, rows, cv.get("lightgbm"), shap)
    write_determinants_attribution(outputs["determinants_attribution"], att)
    if SUPPORT in outputs and cv:
        _support(
            outputs[SUPPORT], cv["lightgbm"], X, blocks, features.loc[rows], nulls, table,
            _threads(cfg),
        )  # fmt: skip


def _support(
    out: Path, cv: CvResult, X: pd.DataFrame, blocks: Mapping[str, str], model: pd.DataFrame,
    nulls: Mapping[str, NullResult], table: pd.DataFrame, threads: int = 1,
) -> None:  # fmt: skip
    """``determinants_shap``, ``determinants_pd`` and ``null_determinants`` (M9-17)."""
    tags = model["locus_tag"].astype(str).tolist()
    write_support("determinants_shap", shap_rows(shap_features(cv, X, blocks), X, tags), out)
    write_support("determinants_pd", pd_rows(cv, X, model, threads), out)
    mean = table[table["metric"] == "cv_r2_mean"]
    observed = dict(zip(mean["model"].astype(str), mean["value"].astype(float), strict=True))
    write_support("null_determinants", null_rows(nulls, observed), out)


@step(
    "determinants",
    inputs=["gene_table", "genes", "tus", "layers", "replicons"],
    optional_inputs=["regions", "region_genes", "gene_context"],
    outputs=["determinants_features", "determinants_model", "determinants_attribution", SUPPORT],
    optional_outputs=[SUPPORT],
)
def determinants_step(ctx: StepContext) -> None:
    """Gene table + contract tables → the three determinants tables."""
    inputs = {k: _path(v) for k, v in ctx.inputs.items()}
    _determinants(inputs, ctx.outputs, ctx.cfg, ctx.log)
