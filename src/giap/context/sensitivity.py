"""Sensitivity framework and the ``regions`` family (issue M4-19; GIAP_SPEC.md §8.6, §2.4 G9,
§7.3.2, §5.2; docs/specs/M4-19-sensitivity-regions.md).

A *variant* changes one analysis choice (missing-abundance policy, single-peptide genes counted as
detected, Top3 instead of iBAQ). A *family* maps a gene table to its flagged IDs; the step reruns
every family under every variant at ``stats.permutations = 199`` and reports the Jaccard index of
each variant's flagged set against the primary's (recomputed at the same B). Families receive
``(genes, cfg, rng)``; extra tables are keyword-only parameters that the step binds by name
(lead decision 5). Pure functions: randomness only from the generator passed in.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import cast

import numpy as np
import pandas as pd

from giap.config.model import AbundanceMissing, GiapConfig
from giap.context.regions import test_all_regions
from giap.genetable.abundance import impute

__all__ = ["COLUMNS", "FAMILIES", "SENSITIVITY_B", "Family", "Variant", "flagged_ids", "jaccard",
           "regions_family", "sensitivity_rows", "variant_genes", "variants"]  # fmt: skip

#: §8.6 ``sensitivity.tsv`` columns (``agreement`` is NA except for the determinants family)
COLUMNS = ("family", "variant", "n_primary", "n_variant", "n_shared", "jaccard", "agreement")
#: §8.6 "Sensitivity reruns use B = 199"
SENSITIVITY_B = 199


@dataclass(frozen=True, slots=True)
class Variant:
    """One sensitivity rerun: its name, the config keys it changes (dotted) and whether
    single-peptide genes count as detected (§7.3.2)."""

    name: str
    overrides: Mapping[str, object]
    single_peptide_detected: bool = False


Family = Callable[[pd.DataFrame, GiapConfig, np.random.Generator], frozenset[str]]


def variants(cfg: GiapConfig, tier: str) -> list[Variant]:
    """Variants of ``output.sensitivity`` for detection ``tier`` (the primary choice is never
    rerun): ``missing=<p>`` per ``abundance.sensitivity`` policy, ``single_peptide_detected`` for
    ``two_peptide``, and ``top3`` under ``full``."""
    level = cfg.output.sensitivity
    if level == "off":
        return []
    out = [
        Variant(f"missing={p}", {"abundance.missing": p})
        for p in dict.fromkeys(cfg.abundance.sensitivity)
        if p != cfg.abundance.missing
    ]
    if tier == "two_peptide":
        out.append(Variant("single_peptide_detected", {}, single_peptide_detected=True))
    if level == "full" and cfg.abundance.method != "top3":
        out.append(Variant("top3", {"abundance.method": "top3"}))
    return out


def jaccard(a: frozenset[str], b: frozenset[str]) -> float | None:
    """|A and B| / |A or B| (intersection over union); ``None`` when both sets are empty."""
    union = len(a | b)
    return None if union == 0 else len(a & b) / union


def sensitivity_rows(
    primary: Mapping[str, frozenset[str]], reruns: Mapping[tuple[str, str], frozenset[str]]
) -> pd.DataFrame:
    """One row per (family, variant) rerun, in ``reruns`` order (columns ``COLUMNS``)."""
    rows = []
    for (fam, var), v in reruns.items():
        p = primary[fam]
        rows.append((fam, var, len(p), len(v), len(p & v), jaccard(p, v), None))
    df = pd.DataFrame(rows, columns=list(COLUMNS))
    return df.astype({"n_primary": "int64", "n_variant": "int64", "n_shared": "int64",
                      "jaccard": "Float64", "agreement": "Float64"})  # fmt: skip


def variant_genes(genes: pd.DataFrame, variant: Variant, cfg: GiapConfig) -> pd.DataFrame:
    """The primary gene table under ``variant``: detection widened to single-peptide genes when
    asked, then the M4-14 imputation redone with the variant's policy on the observed values
    (``top3_log2`` for ``top3``). ``p_detect`` is left as is (the step refits it)."""
    g = genes.copy()
    tags = pd.Index(g["locus_tag"].astype(str))
    det = pd.Series(g["detected"].astype(bool).to_numpy(), index=tags)
    if variant.single_peptide_detected:
        single = (g["detection_status"].astype(str) == "single_peptide_unestimated").to_numpy()
        det |= single
        g.loc[single, "detection_status"] = "detected"
        g["detected"] = det.to_numpy()
    policy = cast(
        AbundanceMissing, variant.overrides.get("abundance.missing", cfg.abundance.missing)
    )
    if variant.overrides.get("abundance.method") == "top3":
        base = g["top3_log2"].astype("float64")
    else:
        base = g["abundance_log2"].astype("float64").where(~g["abundance_imputed"].astype(bool))
    eligible = pd.Series(~g["partial"].astype(bool).to_numpy(), index=tags)
    filled, mask, _ = impute(base.set_axis(tags), det, policy, eligible)
    g["abundance_log2"] = filled.to_numpy("float64")
    g["abundance_imputed"] = mask.to_numpy(bool)
    g["imputation_policy"] = policy
    return g


def flagged_ids(res: pd.DataFrame, alpha: float) -> frozenset[str]:
    """Lead decisions 3-4: ``<region_id>:<contig>:<stat>`` of the rows with
    ``null_model = circular_shift``, p not NA and q <= ``alpha``."""
    q = res["q_value"].astype("float64")
    keep = (
        (res["null_model"].astype(str) == "circular_shift")
        & res["p_value"].astype("float64").notna()
        & q.notna()
        & (q <= alpha)
    )
    sub = res[keep]
    return frozenset(
        f"{r}:{c}:{s}" for r, c, s in zip(sub["region_id"], sub["contig"], sub["stat"], strict=True)
    )


def _tables_from_genes(genes: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """C7 / C8 / C1 stand-ins from the gene table alone: members from ``region_ids``, the region
    type from its §6.5 ID prefix ``<type>_<n>``, contig of its first member; every replicon
    circular (C1 topology is not in the gene table)."""
    members = [
        (str(r), str(t), str(c))
        for t, c, ids in zip(genes["locus_tag"], genes["contig"], genes["region_ids"], strict=True)
        for r in (ids if ids is not None else [])
    ]
    rg = pd.DataFrame([(r, t) for r, t, _ in members], columns=["region_id", "locus_tag"])
    first: dict[str, str] = {}
    for r, _, c in members:
        first.setdefault(r, c)
    regions = pd.DataFrame(
        [(r, r.rsplit("_", 1)[0], c) for r, c in first.items()],
        columns=["region_id", "region_type", "contig"],
    )
    contigs = list(dict.fromkeys(genes["contig"].astype(str)))
    replicons = pd.DataFrame({"contig": contigs, "topology": ["circular"] * len(contigs)})
    return regions, rg, replicons


def regions_family(
    genes: pd.DataFrame,
    cfg: GiapConfig,
    rng: np.random.Generator,
    *,
    regions: pd.DataFrame | None = None,
    region_genes: pd.DataFrame | None = None,
    replicons: pd.DataFrame | None = None,
) -> frozenset[str]:
    """Flagged M4-17 region tests (``flagged_ids`` at ``alpha.context``). The step binds C7
    ``regions``, C8 ``region_genes`` and C1 ``replicons``; without them they are derived from
    ``genes`` (``region_ids``; replicons circular)."""
    if regions is None or region_genes is None or replicons is None:
        d_reg, d_rg, d_rep = _tables_from_genes(genes)
        regions = d_reg if regions is None else regions
        region_genes = d_rg if region_genes is None else region_genes
        replicons = d_rep if replicons is None else replicons
    res = test_all_regions(genes, regions, region_genes, replicons, cfg, rng)
    return flagged_ids(res, cfg.stats.alpha.context)


#: registered families (M4: ``regions``; M5-M8 add theirs)
FAMILIES: dict[str, Family] = {"regions": regions_family}
