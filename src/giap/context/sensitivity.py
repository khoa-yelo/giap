"""Sensitivity framework and the ``regions`` family (issue M4-19; GIAP_SPEC.md §8.6, §2.4 G9,
§7.3.2, §5.2; docs/specs/M4-19-sensitivity-regions.md).

A *variant* changes one analysis choice (missing-abundance policy, single-peptide genes counted as
detected, Top3 instead of iBAQ). A *family* maps a gene table to its flagged IDs; the step reruns
every family under every variant at ``stats.permutations = 199`` and reports the Jaccard index of
each variant's flagged set against the primary's (recomputed at the same B). Families receive
``(genes, cfg, rng)``; the tables a family needs are keyword-only parameters that the step binds
by name (lead decisions 5, 13). Pure functions: randomness only from the generator passed in.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from typing import cast

import numpy as np
import pandas as pd

from giap.config.model import AbundanceMissing, GiapConfig
from giap.context.regions import test_all_regions
from giap.genetable.abundance import impute
from giap.genetable.assemble import GENE_FLAGS

__all__ = ["COLUMNS", "FAMILIES", "SENSITIVITY_B", "Family", "Variant", "flagged_ids", "jaccard",
           "regions_family", "sensitivity_rows", "variant_genes", "variants"]  # fmt: skip

#: §8.6 ``sensitivity.tsv`` columns (``agreement`` is NA except for the determinants family)
COLUMNS = ("family", "variant", "n_primary", "n_variant", "n_shared", "jaccard", "agreement")
#: §8.6 "Sensitivity reruns use B = 199"
SENSITIVITY_B = 199
_RANK = {f: i for i, f in enumerate(GENE_FLAGS)}


@dataclass(frozen=True, slots=True)
class Variant:
    """One sensitivity rerun: its name, the config keys it changes (dotted) and whether
    single-peptide genes count as detected (§7.3.2)."""

    name: str
    overrides: Mapping[str, object]
    single_peptide_detected: bool = False


#: ``(genes, cfg, rng, **tables) -> flagged IDs``; ``...`` because families add keyword-only tables
Family = Callable[..., frozenset[str]]


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
    (``top3_log2`` for ``top3``, with the Top3 raw values and ``abundance_type``), and
    ``imputed_abundance`` flags matching the new imputation. ``p_detect`` is left as is (the
    step refits it)."""
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
    top3 = variant.overrides.get("abundance.method") == "top3"
    if top3:
        base = g["top3_log2"].astype("float64")
        g["abundance_raw"] = np.exp2(base.to_numpy("float64"))
        g["abundance_type"] = "top3"
    else:
        base = g["abundance_log2"].astype("float64").where(~g["abundance_imputed"].astype(bool))
    eligible = pd.Series(~g["partial"].astype(bool).to_numpy(), index=tags)
    filled, mask, _ = impute(base.set_axis(tags), det, policy, eligible)
    g["abundance_log2"] = filled.to_numpy("float64")
    g["abundance_imputed"] = mask.to_numpy(bool)
    g["imputation_policy"] = policy
    g["flags"] = [
        _imputed_flag(f, m) for f, m in zip(g["flags"], g["abundance_imputed"], strict=True)
    ]
    return g


def _imputed_flag(flags: Iterable[str], imputed: bool) -> list[str]:
    """``flags`` with ``imputed_abundance`` iff ``imputed`` (§9.4 vocabulary order)."""
    out = [f for f in flags if f != "imputed_abundance"] + ["imputed_abundance"] * bool(imputed)
    return sorted(out, key=lambda f: _RANK.get(f, len(_RANK)))


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


def regions_family(
    genes: pd.DataFrame,
    cfg: GiapConfig,
    rng: np.random.Generator,
    *,
    regions: pd.DataFrame,
    region_genes: pd.DataFrame,
    replicons: pd.DataFrame,
) -> frozenset[str]:
    """Flagged M4-17 region tests (``flagged_ids`` at ``alpha.context``) on C7 ``regions``,
    C8 ``region_genes`` and C1 ``replicons`` (bound by the step; topology from C1)."""
    res = test_all_regions(genes, regions, region_genes, replicons, cfg, rng)
    return flagged_ids(res, cfg.stats.alpha.context)


#: registered families (M4: ``regions``; M5-M8 add theirs)
FAMILIES: dict[str, Family] = {"regions": regions_family}
