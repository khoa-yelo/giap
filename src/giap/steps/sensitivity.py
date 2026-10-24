"""Sensitivity reruns → ``tables/sensitivity.tsv`` (GIAP_SPEC.md §8.6, §2.4 G9; issue M4-19).

Reads the M4-15 gene table (``genetable/genes.parquet``) and C1 / C5 / C7 / C8 from
``work_dir``. Each variant's gene table comes from ``variant_genes``; ``single_peptide_detected``
also refits ``p_detect`` (M4-13) on the primary's folds, ``for_step("detectability", seed)``
(lead decision 6). Every family runs on the primary and on each variant at
``stats.permutations = 199`` (``islands`` also sets ``stats.bootstrap = 199``, M5-10), each
call from the same family stream ``for_step("sensitivity/<family>", seed)`` (common random
numbers, lead decision 2). Tables a family declares as keyword parameters (``replicons``,
``regions``, ``region_genes``, ``tus``) are bound with ``functools.partial`` (lead decision 5).
``operon_flags`` (M6-10) also takes C2 ``genes`` (as ``gene_models``), ``pairs_raw``, C3
``features`` and the optional evidence layers; without C2, ``pairs_raw`` or C3 it is skipped
(lead decision 4). A variant's config is ``base`` with its ``overrides`` and ``config_patch``.

M8-16: ``enriched_categories`` takes the collections built from the function resources
(``FUNCTION_RESOURCES``, as the ``function`` step builds them on the primary gene table), the
optional ``tu_folds`` as ``folds`` and the lazy R fgsea engine; it is skipped without
``function`` in ``modules`` or without any resource. ``determinants`` takes the M8-10
``determinants_features`` table; it is skipped without it, and writes no row when the primary
rerun is below the 100-row minimum. Each skip is an event naming the family; the determinants
CV R² of every run is the event ``sensitivity.determinants_cv_r2``.

The registered step ``sensitivity`` (M4-19 decision 16, wired by MVP-01) takes those files as
explicit inputs instead of a work directory, so ``giap run`` can call it like any other step;
``run_sensitivity`` keeps the work-directory form for the local executor and the tests.
"""

from __future__ import annotations

import functools
import inspect
import logging
import math
from collections.abc import Mapping
from pathlib import Path

import pandas as pd

from giap.annotate.layers import replicons_of
from giap.config import GiapConfig
from giap.context import sensitivity as S
from giap.determinants.sensitivity import DeterminantsSensitivity
from giap.function.collections import build_collections
from giap.genetable.covariates import COVARIATE_COLUMNS
from giap.genetable.detectability import fit_detectability
from giap.genetable.eligibility import eligible
from giap.io.determinants_tables import read_determinants_features
from giap.io.function_resources import read_resources
from giap.io.genes_tsv import read_genes_parquet
from giap.io.sensitivity_tsv import write_sensitivity_tsv
from giap.io.work_tables import load_tables
from giap.rng import for_step
from giap.setupdb.installed import FUNCTION_RESOURCES
from giap.steps import function as _function
from giap.steps.gene_table import check_issues, contract_paths
from giap.steps.operons import LAYERS, read_layer
from giap.steps.registry import InputPath, StepContext, step

__all__ = ["run_sensitivity", "sensitivity", "with_overrides"]

_log = logging.getLogger("giap.steps.sensitivity")
#: C1 / C7 / C8 / C5, bound to the families that take them (lead decisions 5, 12, 13; M5-10)
_TABLES = ("replicons", "regions", "region_genes", "tus")
#: the ``operon_flags`` tables (M6-10 lead decision 4): step key -> family keyword
_PAIRS = {"genes": "gene_models", "pairs_raw": "pairs_raw", "features": "features"}
_OPERONS = "operon_flags"
_ENRICHED = "enriched_categories"
_DETERMINANTS = "determinants"
#: the M8-16 inputs: function resources, TU folds (ORA) and the determinants features table
_M8 = (*FUNCTION_RESOURCES, "tu_folds", "determinants_features")


def with_overrides(cfg: GiapConfig, overrides: Mapping[str, object]) -> GiapConfig:
    """``cfg`` with dotted keys (``"abundance.missing"``) replaced, re-validated."""
    raw = cfg.model_dump()
    for key, value in overrides.items():
        *path, last = key.split(".")
        node = raw
        for p in path:
            node = node[p]
        node[last] = value
    return GiapConfig.model_validate(raw)


def _refit(genes: pd.DataFrame, cfg: GiapConfig) -> pd.DataFrame:
    """``genes`` with ``p_detect`` refitted on its detection calls (primary folds)."""
    idx = pd.Index(genes["locus_tag"].astype(str))
    det = pd.Series(genes["detected"].astype(bool).to_numpy(), index=idx, dtype="boolean")
    det[genes["partial"].astype(bool).to_numpy()] = pd.NA
    cols = list(COVARIATE_COLUMNS)
    cov = genes[["locus_tag", *cols]].astype(dict.fromkeys(cols, "float64"))
    tus = pd.Series(genes["tu_id"].to_numpy(), index=idx)
    res = fit_detectability(det, cov, tus, cfg.detectability, for_step("detectability", cfg.seed))
    out = genes.copy()
    p = res.p_detect.set_axis(res.p_detect.index.astype(str))
    out["p_detect"] = p.reindex(idx).to_numpy(dtype="float64")
    return out


def _bind(family: S.Family, tables: Mapping[str, object]) -> S.Family:
    """``family`` with the tables it names as parameters bound (a family taking only
    ``(genes, cfg, rng)`` is returned as is)."""
    params = inspect.signature(family).parameters
    kw = {k: v for k, v in tables.items() if k in params}
    return functools.partial(family, **kw) if kw else family


def _sensitivity(gene_table: Path, tables: Mapping[str, Path], cfg: GiapConfig, out: Path) -> None:
    """Rerun every registered family under every variant and write the §8.6 table to ``out``.

    ``gene_table`` is M4-15's ``genes.parquet``; ``tables`` maps the C1 / C5 / C7 / C8 contract keys
    (``_TABLES``) and any of the pair keys (``_PAIRS``, ``LAYERS``) to their files.
    """
    genes = read_genes_parquet(gene_table)
    tier = str(genes["detection_tier"].iloc[0]) if len(genes) else ""
    variants = S.variants(cfg, tier)
    if any(v.name == "top3" for v in variants) and not genes["top3_log2"].notna().any():
        _log.info("sensitivity: top3 skipped: no peptide intensities")
        variants = [v for v in variants if v.name != "top3"]
    loaded, issues = load_tables(tables, _TABLES)
    check_issues(issues, _log, "sensitivity")
    frames: dict[str, object] = dict(loaded)
    families = dict(S.FAMILIES)
    # per-family tables: ``determinants`` and ``operon_flags`` both name a ``features`` keyword
    # (M8-07's TSV vs C3), so the M8-16 tables are bound to their own family only
    own: dict[str, dict[str, object]] = {}
    for name in (_ENRICHED, _DETERMINANTS):
        if name not in families:
            continue
        extra = (_function_inputs(genes, tables, cfg) if name == _ENRICHED
                 else _determinants_inputs(tables))  # fmt: skip
        if extra is None:
            families.pop(name)
        else:
            own[name] = extra
    pairs = _pair_tables(tables, loaded, cfg)
    if pairs is None:  # the UniOP-threshold patch only means something with the pair tables
        families.pop(_OPERONS, None)
    else:
        frames = {**frames, **pairs}
        if "relative" in set(genes["abundance_scale"].dropna().astype(str)):
            _log.info("sensitivity: %s: abundance scale is relative; no pair is tested, every "
                      "%s set is empty (jaccard NA)", _OPERONS, _OPERONS)  # fmt: skip
    # a variant rerun for no remaining family (``pair_threshold=*`` without the pair tables) is
    # dropped before its gene table is built (M6-10 lead decision 8)
    variants = [v for v in variants if any(v.runs_for(f) for f in families)]
    base = with_overrides(cfg, {"stats.permutations": S.SENSITIVITY_B})
    runs: list[tuple[S.Variant, pd.DataFrame, GiapConfig]] = []
    for v in variants:
        g = S.variant_genes(genes, v, cfg)
        runs.append((v, _refit(g, cfg) if v.single_peptide_detected else g,
                     with_overrides(base, {**v.overrides, **v.config_patch})))  # fmt: skip
    primary: dict[str, S.Flagged] = {}
    reruns: dict[tuple[str, str], S.Flagged] = {}
    for name, family in families.items():
        fam = _bind(family, {**frames, **own.get(name, {})})
        p = fam(genes, base, for_step(f"sensitivity/{name}", cfg.seed))
        if p is None:
            _skip(name, _below_minimum_cause(genes))
            continue
        primary[name] = p
        for v, g, vcfg in runs:
            if v.runs_for(name):
                r = fam(g, vcfg, for_step(f"sensitivity/{name}", cfg.seed))
                if r is None:  # a variant below the minimum has no row
                    _variant_skip(name, v.name, _below_minimum_cause(g))
                else:
                    reruns[(name, v.name)] = r
        _cv_event(name, primary, reruns)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_sensitivity_tsv(out, S.sensitivity_rows(primary, reruns))
    _log.info("sensitivity: %d rows (%d variants)", len(reruns), len(variants))


def _skip(family: str, why: str) -> None:
    event = {"event": f"sensitivity.{family}_skipped", "data": {"reason": why}}
    _log.info("sensitivity: %s skipped: %s", family, why, extra=event)


def _variant_skip(family: str, variant: str, why: str) -> None:
    event = {"event": f"sensitivity.{family}_variant_skipped",
             "data": {"variant": variant, "reason": why}}  # fmt: skip
    _log.info("sensitivity: %s variant %s skipped (no row): %s", family, variant, why, extra=event)


def _below_minimum_cause(genes: pd.DataFrame) -> str:
    """Why a family returned ``None``: no abundance at all, or too few model rows."""
    if (
        "abundance_log2" in genes
        and pd.to_numeric(genes["abundance_log2"], errors="coerce").notna().any()
    ):
        return "too few model rows (below the 100-row / 25-TU minimum)"
    return "no abundance (every abundance_log2 is NA)"


def _function_inputs(
    genes: pd.DataFrame, paths: Mapping[str, Path], cfg: GiapConfig
) -> dict[str, object] | None:
    """``collections``, ``folds`` and ``engine`` of ``enriched_categories``, or ``None`` (logged)
    when ``function`` is not in ``modules`` or no function resource is given."""
    if "function" not in cfg.modules:
        _skip(_ENRICHED, "the function module is off")
        return None
    res_paths = {k: paths[k] for k in FUNCTION_RESOURCES if k in paths}
    if not res_paths:
        _skip(_ENRICHED, "no function resource")
        return None
    colls, _ = build_collections(genes, eligible(genes), read_resources(res_paths))
    folds = _function.read_folds(paths["tu_folds"]) if "tu_folds" in paths else {}
    return {"collections": colls, "folds": folds or None,
            "engine": _function._LazyRFgsea(_log)}  # fmt: skip


def _determinants_inputs(paths: Mapping[str, Path]) -> dict[str, object] | None:
    """``features`` of ``determinants`` (M8-10's TSV), or ``None`` (logged) when absent."""
    if "determinants_features" not in paths:
        _skip(_DETERMINANTS, "no determinants_features table")
        return None
    return {"features": read_determinants_features(paths["determinants_features"])}


def _cv_event(
    name: str, primary: Mapping[str, S.Flagged], reruns: Mapping[tuple[str, str], S.Flagged]
) -> None:
    """The CV R² of every determinants run (§8.6 has no R² column; spec item 2)."""
    p = primary[name]
    if not isinstance(p, DeterminantsSensitivity):
        return
    runs = {"primary": p, **{v: r for (f, v), r in reruns.items() if f == name}}
    data = {k: (float(r.cv_r2) if math.isfinite(r.cv_r2) else None)
            for k, r in runs.items() if isinstance(r, DeterminantsSensitivity)}  # fmt: skip
    _log.info("sensitivity: %s CV R² per run", name,
              extra={"event": "sensitivity.determinants_cv_r2", "data": data})  # fmt: skip


def _pair_tables(
    paths: Mapping[str, Path], frames: Mapping[str, pd.DataFrame], cfg: GiapConfig
) -> dict[str, pd.DataFrame] | None:
    """The ``operon_flags`` tables under their family keywords, or ``None`` (logged) when C2,
    ``pairs_raw`` or C3 is absent: the family and the ``pair_threshold`` variants are then
    skipped. Evidence layers are read as the ``operons`` step reads them."""
    absent = sorted(k for k in _PAIRS if k not in paths)
    if absent:
        _log.info("sensitivity: %s and the pair_threshold variants skipped: no %s", _OPERONS,
                  ", ".join(absent),
                  extra={"event": "sensitivity.operon_flags_skipped",
                         "data": {"missing": absent}})  # fmt: skip
        return None
    out = {kw: pd.read_parquet(paths[k]) for k, kw in _PAIRS.items()}
    reps = replicons_of(frames["replicons"])
    for key in LAYERS:
        if key in paths:
            out[key] = read_layer(cfg, key, paths[key], reps, frames["tus"])
    return out


def run_sensitivity(work_dir: Path, cfg: GiapConfig) -> Path | None:
    """Rerun every registered family under every variant; writes and returns
    ``work_dir/tables/sensitivity.tsv`` (``None`` when ``output.sensitivity = off``)."""
    out = work_dir / "tables" / "sensitivity.tsv"
    out.unlink(missing_ok=True)  # a failed or ``off`` run never leaves a stale table
    if cfg.output.sensitivity == "off":
        _log.info(
            "sensitivity: output.sensitivity = off; families %s not rerun", sorted(S.FAMILIES)
        )
        return None
    tables = contract_paths(work_dir)  # either work layout (M4-21 lead decision 3)
    _sensitivity(work_dir / "genetable" / "genes.parquet", tables, cfg, out)
    return out


def _one(p: InputPath) -> Path:
    if isinstance(p, list):
        raise ValueError(f"expected one path, got {len(p)}")
    return p


@step("sensitivity", inputs=["gene_table", *_TABLES], optional_inputs=[*_PAIRS, *LAYERS, *_M8],
      outputs=["sensitivity"])  # fmt: skip
def sensitivity(ctx: StepContext) -> None:
    """Gene table + C1 / C5 / C7 / C8 (+ the M6-10 pair tables) → the §8.6 sensitivity table
    (M4-19 decision 16).

    ``output.sensitivity = off`` writes the table with its header and no row: a declared step
    output must exist, so the pipeline cannot hand a half-written run to the next step.
    """
    out = ctx.outputs["sensitivity"]
    if ctx.cfg.output.sensitivity == "off":
        ctx.log.info(
            "sensitivity: output.sensitivity = off; families %s not rerun, %s has no rows",
            sorted(S.FAMILIES), out.name,
        )  # fmt: skip
        out.parent.mkdir(parents=True, exist_ok=True)
        write_sensitivity_tsv(out, S.sensitivity_rows({}, {}))
        return
    tables = {k: _one(v) for k, v in ctx.inputs.items() if k in (*_TABLES, *_PAIRS, *LAYERS, *_M8)}
    _sensitivity(_one(ctx.inputs["gene_table"]), tables, ctx.cfg, out)
