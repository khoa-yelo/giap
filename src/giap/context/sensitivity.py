"""Sensitivity framework and the ``regions`` family (issue M4-19; GIAP_SPEC.md §8.6, §2.4 G9,
§7.3.2, §5.2; docs/specs/M4-19-sensitivity-regions.md).

A *variant* changes one analysis choice (missing-abundance policy, single-peptide genes counted as
detected, Top3 instead of iBAQ). A *family* maps a gene table to its flagged IDs; the step reruns
every family under every variant at ``stats.permutations = 199`` and reports the Jaccard index of
each variant's flagged set against the primary's (recomputed at the same B). Families receive
``(genes, cfg, rng)``; the tables a family needs are keyword-only parameters that the step binds
by name (lead decisions 5, 13). Pure functions: randomness only from the generator passed in.

The ``islands`` family (M5-10) flags ``<direction>:<locus_tag>`` for every gene of every M5-06
island, rerunning the scan at ``stats.bootstrap = stats.permutations = 199``.

The ``operon_flags`` family (M6-10) flags the ``pair_id`` values of the M6-08 chain's
``operon_flags`` frame. Its ``full`` variants ``pair_threshold=<t>`` carry a ``config_patch``
(lead decision 1): they change the UniOP threshold, hence the pair classes, not the abundances.

M8-16 adds ``enriched_categories`` (``giap.function.sensitivity``) and ``determinants``
(``giap.determinants.sensitivity``). A determinants result is a ``DeterminantsSensitivity``:
``n_*`` / ``jaccard`` use its ``flagged`` set and ``agreement`` the dominant blocks of the genes
the features explain in both runs (M8-16c, lead decision 6).
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from typing import cast

import numpy as np
import pandas as pd

from giap.config.model import AbundanceMissing, GiapConfig
from giap.context import islands as _islands
from giap.context.regions import replicon_index, test_all_regions
from giap.determinants.sensitivity import DeterminantsSensitivity, determinants_family
from giap.function.sensitivity import enriched_categories_family
from giap.genetable.abundance import impute
from giap.genetable.assemble import GENE_FLAGS
from giap.steps.operons import run_operons

__all__ = ["COLUMNS", "FAMILIES", "PAIR_THRESHOLD", "SENSITIVITY_B", "UNIOP_THRESHOLDS", "Family",
           "Flagged", "Variant", "agreement", "flagged_ids", "islands_family", "jaccard",
           "operon_flags_family", "regions_family", "sensitivity_rows", "variant_genes",
           "variants"]  # fmt: skip

#: §8.6 ``sensitivity.tsv`` columns (``agreement`` is NA except for the determinants family)
COLUMNS = ("family", "variant", "n_primary", "n_variant", "n_shared", "jaccard", "agreement")
#: §8.6 "Sensitivity reruns use B = 199"
SENSITIVITY_B = 199
#: §8.6 ``full``: "UniOP thresholds 0.3 and 0.7" (dotted config key of the patch)
UNIOP_THRESHOLDS = (0.3, 0.7)
PAIR_THRESHOLD = "annotation.uniop.pair_threshold"
_RANK = {f: i for i, f in enumerate(GENE_FLAGS)}


@dataclass(frozen=True, slots=True)
class Variant:
    """One sensitivity rerun: its name, the config keys it changes (dotted), whether
    single-peptide genes count as detected (§7.3.2), a config patch (dotted keys) that changes
    the analysis but not the abundances (M6-10 lead decision 1; empty for M4-19's variants) and
    the families it is rerun for (``None`` = every family; M6-10 lead decision 8)."""

    name: str
    overrides: Mapping[str, object]
    single_peptide_detected: bool = False
    config_patch: Mapping[str, object] = field(default_factory=dict)
    families: tuple[str, ...] | None = None

    @property
    def patch_only(self) -> bool:
        """True when the variant only patches config: no overrides, no single-peptide switch
        (its gene table is the primary's)."""
        return bool(self.config_patch) and not self.overrides and not self.single_peptide_detected

    def runs_for(self, family: str) -> bool:
        """Whether this variant is rerun for ``family`` (``families`` is ``None`` or names it)."""
        return self.families is None or family in self.families


#: a family's result: flagged IDs, or the determinants rerun (flagged genes and dominant blocks)
Flagged = frozenset[str] | DeterminantsSensitivity
#: ``(genes, cfg, rng, **tables) -> Flagged`` (``None``: too little data, no row); ``...`` because
#: families add keyword-only tables
Family = Callable[..., Flagged | None]


def variants(cfg: GiapConfig, tier: str) -> list[Variant]:
    """Variants of ``output.sensitivity`` for detection ``tier`` (the primary choice is never
    rerun): ``missing=<p>`` per ``abundance.sensitivity`` policy, ``single_peptide_detected`` for
    ``two_peptide``, and ``top3`` then ``pair_threshold=<t>`` (``UNIOP_THRESHOLDS``) under
    ``full``."""
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
    if level == "full":
        out += [
            Variant(
                f"pair_threshold={t:g}",
                {},
                config_patch={PAIR_THRESHOLD: t},
                families=("operon_flags",),
            )
            for t in UNIOP_THRESHOLDS
            if t != cfg.annotation.uniop.pair_threshold
        ]
    return out


def jaccard(a: frozenset[str], b: frozenset[str]) -> float | None:
    """|A and B| / |A or B| (intersection over union); ``None`` when both sets are empty."""
    union = len(a | b)
    return None if union == 0 else len(a & b) / union


def agreement(a: DeterminantsSensitivity, b: DeterminantsSensitivity) -> float | None:
    """Share of the genes explained in both runs whose dominant block is the same.

    A gene is explained in a run when it is in that run's ``dominant`` map and not in its
    ``flagged`` set (``abundance_not_explained_by_sequence``); a gene flagged in either run is
    left out. Genes with a block of ``None`` in both runs are left out; ``None`` in one run only
    is a disagreement. ``None`` (NA) when no gene is left (M8-16 lead decisions 5-6)."""
    shared = [g for g in a.dominant.keys() & b.dominant.keys()
              if g not in a.flagged and g not in b.flagged
              and not (a.dominant[g] is None and b.dominant[g] is None)]  # fmt: skip
    if not shared:
        return None
    same = sum(a.dominant[g] is not None and a.dominant[g] == b.dominant[g] for g in shared)
    return same / len(shared)


def sensitivity_rows(
    primary: Mapping[str, Flagged], reruns: Mapping[tuple[str, str], Flagged]
) -> pd.DataFrame:
    """One row per (family, variant) rerun, in ``reruns`` order (columns ``COLUMNS``);
    ``agreement`` only for ``DeterminantsSensitivity`` results."""
    rows = []
    for (fam, var), v in reruns.items():
        p, agree = primary[fam], None
        if isinstance(p, DeterminantsSensitivity) and isinstance(v, DeterminantsSensitivity):
            agree = agreement(p, v)
        ps = p.flagged if isinstance(p, DeterminantsSensitivity) else p
        vs = v.flagged if isinstance(v, DeterminantsSensitivity) else v
        rows.append((fam, var, len(ps), len(vs), len(ps & vs), jaccard(ps, vs), agree))
    df = pd.DataFrame(rows, columns=list(COLUMNS))
    return df.astype({"n_primary": "int64", "n_variant": "int64", "n_shared": "int64",
                      "jaccard": "Float64", "agreement": "Float64"})  # fmt: skip


def variant_genes(genes: pd.DataFrame, variant: Variant, cfg: GiapConfig) -> pd.DataFrame:
    """The primary gene table under ``variant``: detection widened to single-peptide genes when
    asked, then the M4-14 imputation redone with the variant's policy on the observed values
    (``top3_log2`` for ``top3``, with the Top3 raw values and ``abundance_type``), and
    ``imputed_abundance`` flags matching the new imputation. ``p_detect`` is left as is (the
    step refits it). A patch-only variant gives the primary table back (a copy)."""
    g = genes.copy()
    if variant.patch_only:
        return g
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


def islands_family(
    genes: pd.DataFrame,
    cfg: GiapConfig,
    rng: np.random.Generator,
    *,
    tus: pd.DataFrame,
    replicons: pd.DataFrame,
) -> frozenset[str]:
    """``<direction>:<locus_tag>`` over the best-window genes of every M5-06 island (M5-10).

    The M5-06 pipeline as ``context_stats`` runs it (``tu_series`` per C1 replicon,
    ``scan_replicons`` then ``islands_frame`` on ``rng``) on C5 ``tus`` and C1 ``replicons``,
    with ``stats.bootstrap`` and ``stats.permutations`` at ``SENSITIVITY_B`` (§8.6). The
    direction is part of the ID, so a gene moving from a high to a low island is a change.
    """
    stats = cfg.stats.model_copy(update={"bootstrap": SENSITIVITY_B, "permutations": SENSITIVITY_B})
    run_cfg = cfg.model_copy(update={"stats": stats})
    span = float(cfg.context.islands.detrend_span_frac)
    series = [
        _islands.tu_series(
            replicon_index(genes, str(contig)),
            tus,
            contig=str(contig),
            circular=str(topology) == "circular",
            length=int(length),
            span_frac=span,
        )
        for contig, topology, length in zip(
            replicons["contig"], replicons["topology"], replicons["length"], strict=True
        )
    ]
    scan = _islands.scan_replicons(series, run_cfg, rng)
    frame = _islands.islands_frame(scan, genes, run_cfg, rng)
    return frozenset(
        f"{d}:{g}"
        for d, members in zip(frame["direction"], frame["genes"], strict=True)
        for g in members
    )


def operon_flags_family(
    genes: pd.DataFrame,
    cfg: GiapConfig,
    rng: np.random.Generator,
    *,
    pairs_raw: pd.DataFrame,
    tus: pd.DataFrame,
    replicons: pd.DataFrame,
    gene_models: pd.DataFrame,
    features: pd.DataFrame,
    promoter_windows: pd.DataFrame | None = None,
    terminators: pd.DataFrame | None = None,
    tx_rates: pd.DataFrame | None = None,
) -> frozenset[str]:
    """The ``pair_id`` values of the M6-08 ``operon_flags`` frame (``run_operons``) on the gene
    table ``genes``, C2 ``gene_models``, C1, C5, ``pairs_raw`` and C3 ``features`` (M6-10).

    An evidence layer left at ``None`` is a missing layer (``layer_missing:<layer>``), as in the
    ``operons`` step. Module 4 draws no random number, so ``rng`` is unused. Pure.
    """
    given = {"promoter_windows": promoter_windows, "terminators": terminators, "tx_rates": tx_rates}
    layers = {k: v for k, v in given.items() if v is not None}
    flagged = run_operons(genes, gene_models, replicons, tus, pairs_raw, features, layers, cfg)[2]
    return frozenset(str(p) for p in flagged["pair_id"])


#: registered families (M4: ``regions``; M5: ``islands``; M6: ``operon_flags``; M8-16:
#: ``enriched_categories``, ``determinants``)
FAMILIES: dict[str, Family] = {
    "regions": regions_family, "islands": islands_family, "operon_flags": operon_flags_family,
    "enriched_categories": enriched_categories_family, "determinants": determinants_family,
}  # fmt: skip
