"""Sensitivity reruns → ``tables/sensitivity.tsv`` (GIAP_SPEC.md §8.6, §2.4 G9; issue M4-19).

Reads the M4-15 gene table (``genetable/genes.parquet``) and C1 / C7 / C8 from ``work_dir``. Each
variant's gene table comes from ``variant_genes``; ``single_peptide_detected`` also refits
``p_detect`` (M4-13) on the primary's folds, ``for_step("detectability", seed)`` (lead decision
6). Every family runs on the primary and on each variant at ``stats.permutations = 199``, each
call from the same family stream ``for_step("sensitivity/<family>", seed)`` (common random
numbers, lead decision 2). Tables a family declares as keyword parameters (``replicons``,
``regions``, ``region_genes``) are bound with ``functools.partial`` (lead decision 5).
"""

from __future__ import annotations

import functools
import inspect
import logging
from collections.abc import Mapping
from pathlib import Path

import pandas as pd

from giap.config import GiapConfig
from giap.context import sensitivity as S
from giap.genetable.covariates import COVARIATE_COLUMNS
from giap.genetable.detectability import fit_detectability
from giap.io.genes_tsv import read_genes_parquet
from giap.io.sensitivity_tsv import write_sensitivity_tsv
from giap.io.work_tables import default_paths, load_m4_inputs
from giap.rng import for_step

__all__ = ["run_sensitivity", "with_overrides"]

_log = logging.getLogger("giap.steps.sensitivity")


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


def _bind(family: S.Family, tables: Mapping[str, pd.DataFrame]) -> S.Family:
    params = inspect.signature(family).parameters
    kw = {k: v for k, v in tables.items() if k in params}
    return functools.partial(family, **kw) if kw else family


def run_sensitivity(work_dir: Path, cfg: GiapConfig) -> Path | None:
    """Rerun every registered family under every variant; writes and returns
    ``work_dir/tables/sensitivity.tsv`` (``None`` when ``output.sensitivity = off``)."""
    if cfg.output.sensitivity == "off":
        _log.info(
            "sensitivity: output.sensitivity = off; families %s not rerun", sorted(S.FAMILIES)
        )
        return None
    genes = read_genes_parquet(work_dir / "genetable" / "genes.parquet")
    tier = str(genes["detection_tier"].iloc[0]) if len(genes) else ""
    variants = S.variants(cfg, tier)
    if any(v.name == "top3" for v in variants) and not genes["top3_log2"].notna().any():
        _log.info("sensitivity: top3 skipped: no peptide intensities")
        variants = [v for v in variants if v.name != "top3"]
    paths = {k: p for k, p in default_paths(work_dir).items() if p.exists()}
    m4, _ = load_m4_inputs(paths)
    tables = {"replicons": m4.replicons, "regions": m4.regions, "region_genes": m4.region_genes}
    base = with_overrides(cfg, {"stats.permutations": S.SENSITIVITY_B})
    runs: list[tuple[str, pd.DataFrame, GiapConfig]] = []
    for v in variants:
        g = S.variant_genes(genes, v, cfg)
        runs.append((v.name, _refit(g, cfg) if v.single_peptide_detected else g,
                     with_overrides(base, v.overrides)))  # fmt: skip
    primary: dict[str, frozenset[str]] = {}
    reruns: dict[tuple[str, str], frozenset[str]] = {}
    for name, family in S.FAMILIES.items():
        fam = _bind(family, tables)
        primary[name] = fam(genes, base, for_step(f"sensitivity/{name}", cfg.seed))
        for vname, g, vcfg in runs:
            reruns[(name, vname)] = fam(g, vcfg, for_step(f"sensitivity/{name}", cfg.seed))
    out = work_dir / "tables" / "sensitivity.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    write_sensitivity_tsv(out, S.sensitivity_rows(primary, reruns))
    _log.info("sensitivity: %d rows (%d variants)", len(reruns), len(variants))
    return out
