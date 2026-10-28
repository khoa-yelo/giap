"""Missed operons: NO-pair scores against CV and SW decoys, q-values and local fdr (GIPPA_SPEC.md
§7.4.4, §8.1 row "Missed operon", §9.3; docs/specs/M6-06-missed-operons.md, pass 1).

Targets are the NO pairs with both genes eligible (not ``partial``), both detected and both
``abundance_log2`` finite, scored ``s = -|Delta|``; imputed abundances and ``rna_between`` pairs
stay (lead decision 2). Decoys are the CV and SW pairs with both genes detected and finite
abundances (``partial`` is a target rule only). Both are binned by ``intergenic_bp`` over
``operons.distance_bins_bp`` and handed to the M6-01 engine ``gippa.stats.decoy_fdr.max_decoy_fdr``,
which reweights each decoy set to the targets' distance mix, takes the larger FDR, monotonises it
into ``q_value`` and computes the local fdr (ranking only, never a call). A decoy FDR has no p, so
``p_value`` is NA (lead decision 1).

**Power.** The family controls the FDR but is near-powerless at realistic effects (M6-01 lead
decision 12): with ``pi0 = 1`` a call at q <= 0.1 needs a missed operon ~25-50 times more
concentrated near Delta = 0 than a null pair, so effects of 0.25-1.0 SD give essentially no calls.
An empty call list is the expected outcome, not evidence that a genome has no missed operons.

Pure and deterministic: no I/O, no RNG; logging only.
"""

from __future__ import annotations

import logging
import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pandas as pd

from gippa.config import GippaConfig
from gippa.operons.pairs import TEST_COLUMNS
from gippa.stats.decoy_fdr import bin_index, max_decoy_fdr

__all__ = ["MIN_NO_PAIRS", "MissedOperons", "missed_operons"]

log = logging.getLogger(__name__)

#: fewest scored NO pairs for the decoy FDR (lead decision 3, as M6-05)
MIN_NO_PAIRS: int = 20
TEST = "missed_operon_decoy"
METHOD = "decoy_fdr"
TOO_FEW = "too_few_no_pairs"
RELATIVE = "relative_scale"
UNESTIMATED = "unestimated"
#: scored rows have a q but no CI for a single pair (as M6-05, #147 review F1)
CI_UNESTIMATED = "ci_unestimated"
IMPUTED = "abundance_imputed"
ROW_COLUMNS: tuple[str, ...] = ("pair_id", "test", *TEST_COLUMNS, "which", "flags")
#: ``operons.decoys`` name -> pair class and ``bin_counts`` column
DECOY_CLASS: dict[str, tuple[str, str]] = {
    "convergent": ("CV", "N_cv"),
    "sandwich": ("SW", "N_sw"),
}


@dataclass(frozen=True, slots=True)
class MissedOperons:
    """Result of :func:`missed_operons`.

    Attributes:
        rows: one row per NO pair: ``pair_id``, ``test``, the §9.2 test columns, ``which`` (the
            decoy set that decided the FDR, "" when unscored) and ``flags`` (a list per row).
            Scored targets have ``test = "missed_operon_decoy"``; the rest keep ``"none"``.
        n_targets: scored NO pairs (also when too few to test).
        n_decoys: usable decoys per configured decoy set, keyed by config name.
        bin_counts: one row per distance bin, in bin order: ``bin``, ``lower_bp``, ``upper_bp``,
            ``N_NO``, ``N_cv``, ``N_sw`` (usable pairs; 0 for an unconfigured set).
    """

    rows: pd.DataFrame
    n_targets: int
    n_decoys: dict[str, int]
    bin_counts: pd.DataFrame


def _bool(s: pd.Series) -> npt.NDArray[np.bool_]:
    return np.asarray(s.astype("boolean").fillna(False).to_numpy(dtype=bool), dtype=np.bool_)


def _flags(v: object) -> list[str]:
    if isinstance(v, str):
        return [f for f in v.replace(",", ";").split(";") if f]
    if isinstance(v, list | tuple | np.ndarray):
        return [str(f) for f in v]
    return []


def _per_gene(genes: pd.DataFrame | None, col: str, tags: pd.Series) -> npt.NDArray[np.bool_]:
    if genes is None or col not in genes.columns:
        return np.zeros(len(tags), dtype=np.bool_)
    m = genes.drop_duplicates("locus_tag").set_index("locus_tag")[col]
    return _bool(tags.map(m))


def _is_relative(genes: pd.DataFrame | None) -> bool:
    """The run's abundance scale is relative (gene table ``abundance_scale``, §7.3.4)."""
    if genes is None or "abundance_scale" not in genes.columns:
        return False
    return "relative" in set(genes["abundance_scale"].dropna().astype(str))


def _measured(pairs: pd.DataFrame) -> npt.NDArray[np.bool_]:
    """Both detected, both abundances and Delta finite."""
    a = pairs["abundance_log2_a"].astype(float).to_numpy()
    b = pairs["abundance_log2_b"].astype(float).to_numpy()
    d = pairs["delta_log2"].astype(float).to_numpy()
    ok = _bool(pairs["detected_a"]) & _bool(pairs["detected_b"])
    ok &= np.isfinite(a) & np.isfinite(b) & np.isfinite(d)
    return np.asarray(ok, dtype=np.bool_)


def _empty_rows(no: pd.DataFrame) -> pd.DataFrame:
    n = len(no)
    rows = pd.DataFrame({"pair_id": no["pair_id"].to_numpy(), "test": ["none"] * n})
    for col in TEST_COLUMNS:
        rows[col] = None if col in ("method", "null_model") else math.nan
    rows["n_perm"] = pd.array([pd.NA] * n, dtype="Int64")
    rows["which"] = ""
    return rows


def _bin_counts(edges: Sequence[float], counts: dict[str, npt.NDArray[np.int64]]) -> pd.DataFrame:
    n_bins = len(edges) - 1
    zero = np.zeros(n_bins, dtype=np.int64)
    return pd.DataFrame(
        {
            "bin": np.arange(n_bins, dtype=np.int64),
            "lower_bp": np.asarray(edges[:-1], dtype=np.float64),
            "upper_bp": np.asarray(edges[1:], dtype=np.float64),
            **{c: counts.get(c, zero) for c in ("N_NO", "N_cv", "N_sw")},
        }
    )


def missed_operons(
    pairs: pd.DataFrame, cfg: GippaConfig, *, genes: pd.DataFrame | None = None
) -> MissedOperons:
    """Score every NO pair as a candidate missed operon against the CV/SW decoys (§7.4.4).

    Args:
        pairs: the M6-04 pair frame (``gippa.operons.pairs.pair_frame``); not modified.
        cfg: run configuration (``operons.decoys``, ``operons.distance_bins_bp``,
            ``operons.pi0``).
        genes: the gene table (``locus_tag``, ``partial``, ``abundance_scale``); ``None`` means
            every gene is eligible and the scale is molar.

    Returns:
        The :class:`MissedOperons` result. ``q_value <= stats.alpha.operon_flags`` marks a pair
        co-abundant (M6-07 names it); ``lfdr`` ranks and never calls. A relative abundance scale
        flags every NO row ``relative_scale`` and ``unestimated``; fewer than ``MIN_NO_PAIRS``
        scored targets or no usable decoy set flags them ``too_few_no_pairs``; neither path scores
        a row. Scored rows carry ``ci_unestimated`` (a q, no CI); a NO row with an imputed gene
        carries ``abundance_imputed``. Eligible pairs outside every distance bin are not scored
        and are counted in a WARNING.

    Raises:
        ValueError: if ``operons.decoys`` is empty.
    """
    oc = cfg.operons
    if not oc.decoys:
        raise ValueError("operons.decoys is empty: at least one decoy set is required")
    if oc.pi0 != 1.0:
        log.warning(
            "missed operons: operons.pi0 = %g scales every q and lfdr with no estimate behind it; "
            "the default 1 is conservative",
            oc.pi0,
        )
    edges = [float(e) for e in oc.distance_bins_bp]
    n_bins = len(edges) - 1
    cls = pairs["class"].astype(str).to_numpy()
    bins = bin_index(pairs["intergenic_bp"].astype(float).to_numpy(), edges)
    measured = _measured(pairs)
    in_bin = (bins >= 0) & (bins < n_bins)
    usable = measured & in_bin
    score = -np.abs(pairs["delta_log2"].astype(float).to_numpy())

    is_no = cls == "NO"
    partial = _per_gene(genes, "partial", pairs["gene_a"]) | _per_gene(
        genes, "partial", pairs["gene_b"]
    )
    target = is_no & usable & ~partial
    imputed = _per_gene(genes, "abundance_imputed", pairs["gene_a"]) | _per_gene(
        genes, "abundance_imputed", pairs["gene_b"]
    )
    is_decoy = np.isin(cls, [DECOY_CLASS[n][0] for n in oc.decoys])
    out_no = int((is_no & measured & ~partial & ~in_bin).sum())
    out_d = int((is_decoy & measured & ~in_bin).sum())
    if out_no or out_d:
        log.warning(
            "missed operons: %d eligible NO pair(s) and %d decoy pair(s) fall outside every "
            "operons.distance_bins_bp bin and are not scored",
            out_no,
            out_d,
        )
    counts: dict[str, npt.NDArray[np.int64]] = {"N_NO": np.bincount(bins[target], minlength=n_bins)}
    decoys: dict[str, tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]]] = {}
    n_decoys: dict[str, int] = {}
    for name in dict.fromkeys(oc.decoys):
        klass, col = DECOY_CLASS[name]
        keep = (cls == klass) & usable
        counts[col] = np.bincount(bins[keep], minlength=n_bins)
        n_decoys[name] = int(keep.sum())
        if keep.any():
            decoys[name] = (score[keep], bins[keep])
        else:
            log.info("missed operons: no usable %s decoy pair; the set is dropped", name)
    bin_counts = _bin_counts(edges, counts)

    no = pairs.loc[is_no]
    rows = _empty_rows(no)
    raw = no["flags"] if "flags" in no.columns else [None] * len(no)
    flags = [_flags(f) for f in raw]
    flags = [[*f, IMPUTED] if imp and IMPUTED not in f else f
             for f, imp in zip(flags, imputed[is_no], strict=True)]  # fmt: skip
    n_targets = int(target.sum())

    extra: list[str] = []
    if _is_relative(genes):
        extra = [RELATIVE, UNESTIMATED]
    elif n_targets < MIN_NO_PAIRS or not decoys:
        log.info(
            "missed operons: %d scored NO pairs (minimum %d), %d usable decoy set(s); untested",
            n_targets,
            MIN_NO_PAIRS,
            len(decoys),
        )
        extra = [TOO_FEW]
    if extra:
        rows["flags"] = [[*f, *(x for x in extra if x not in f)] for f in flags]
        return MissedOperons(rows[list(ROW_COLUMNS)], n_targets, n_decoys, bin_counts)

    for used, (_, d_bin) in decoys.items():
        n_d = np.bincount(d_bin, minlength=n_bins)
        for b in np.flatnonzero((counts["N_NO"] > 0) & (n_d == 0)):
            log.warning(
                "missed operons: distance bin %d [%g, %g) has %d targets and no %s decoys; it "
                "contributes its whole target count to the FDR (conservative)",
                b,
                edges[b],
                edges[b + 1],
                counts["N_NO"][b],
                used,
            )

    t_score, t_bin = score[target], bins[target]
    res = max_decoy_fdr(t_score, t_bin, decoys, n_bins, pi0=oc.pi0)
    idx = np.flatnonzero(target[is_no])
    rows.loc[idx, "test"] = TEST
    rows.loc[idx, "stat"] = t_score
    rows.loc[idx, "effect"] = -t_score
    rows.loc[idx, "q_value"] = res.q
    rows.loc[idx, "lfdr"] = res.lfdr
    rows.loc[idx, "method"] = METHOD
    rows.loc[idx, "null_model"] = "+".join(decoys) + "_reweighted"
    rows.loc[idx, "which"] = res.which.astype(str)
    scored = set(idx.tolist())
    rows["flags"] = [
        [*f, CI_UNESTIMATED] if i in scored and CI_UNESTIMATED not in f else f
        for i, f in enumerate(flags)
    ]
    return MissedOperons(rows[list(ROW_COLUMNS)], n_targets, n_decoys, bin_counts)
