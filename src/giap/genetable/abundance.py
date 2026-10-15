"""Gene-level abundance: iBAQ, Top3, protein-level route, log2, missing values (M4-14).

Spec §2.1 decisions 2-3, §7.3.4 (iBAQ, Top3, protein-level route, transform, missing values),
§7.2.6 (two proteins on one gene, cross-strain proteins), §9.3 abundance columns, §9.4 flags,
§5.2 ``abundance`` config; docs/specs/M4-14-abundance.md.

Lead decisions (spec questions 1-3, proposed defaults): Top3 is the linear mean of the (up to) 3
largest intensities, then log2; the fill statistic uses observed values of eligible (non-partial)
detected genes only; undetected genes keep an observed value (reporting, not inference).

Peptide route: only ``abundance_eligible`` accepted rows with a finite intensity > 0 (the C9
median-across-runs value) are quantified. Protein route: C10 values joined through the M4-12
protein -> gene map, the larger log2 value per gene. Undetected genes are never imputed.
Everything here is pure; the only side effect is a ``no_observed_abundance`` log record.
"""

from __future__ import annotations

import logging
import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from giap.config import GiapConfig
    from giap.config.model import AbundanceMissing

_log = logging.getLogger(__name__)

ABUNDANCE_COLUMNS = (
    "abundance_raw",
    "abundance_log2",
    "abundance_imputed",
    "imputation_policy",
    "abundance_type",
    "abundance_scale",
    "top3_log2",
)
FLAGS_COLUMN = "abundance_flags"
MOLAR_TYPES = frozenset({"ibaq", "top3"})
POLICIES = ("median", "mean", "min", "zero", "exclude")
Source = Literal["peptide_intensity", "protein_table", "none"]


@dataclass(frozen=True, slots=True)
class AbundanceResult:
    """``table``: ``locus_tag`` + §9.3 abundance columns + ``abundance_flags`` (list per gene)."""

    table: pd.DataFrame
    source: Source
    fill_value: float | None


def _positive(values: pd.Series) -> pd.Series:
    """Finite values > 0 (log2 undefined otherwise); everything else becomes NA."""
    v = pd.to_numeric(values, errors="coerce").astype("float64")
    return v.where(np.isfinite(v) & (v > 0))


def _quantified(accepted: pd.DataFrame) -> pd.DataFrame:
    """Eligible accepted rows with a usable intensity: (gene, intensity)."""
    if accepted.empty:
        return pd.DataFrame({"gene": pd.Series([], dtype=object), "intensity": []})
    keep = accepted["abundance_eligible"].astype(bool).to_numpy()
    sub = accepted.loc[keep, ["gene", "intensity"]].copy()
    sub["intensity"] = _positive(sub["intensity"])
    return sub.dropna(subset=["intensity"]).astype({"gene": str})


def peptide_abundance(accepted: pd.DataFrame, n_observable: pd.Series) -> pd.DataFrame:
    """Per gene with >= 1 quantified peptide: ``locus_tag``, ``ibaq``, ``top3`` (both linear) and
    ``n_quant_peptides``. iBAQ = sum / max(1, n_observable); Top3 = mean of the 3 largest values.
    Sums are taken over sorted values so the result does not depend on the row order."""
    tags: list[str] = []
    ibaq: list[float] = []
    top3: list[float] = []
    count: list[int] = []
    for gene, vals in _quantified(accepted).groupby("gene", sort=True)["intensity"]:
        v = sorted(float(x) for x in vals)
        n_obs = n_observable.get(gene, 0)
        n_obs = 0 if pd.isna(n_obs) else int(n_obs)
        top = v[-3:]
        tags.append(str(gene))
        ibaq.append(math.fsum(v) / max(1, n_obs))
        top3.append(math.fsum(top) / len(top))
        count.append(len(v))
    return pd.DataFrame(
        {
            "locus_tag": pd.Series(tags, dtype=object),
            "ibaq": pd.Series(ibaq, dtype="float64"),
            "top3": pd.Series(top3, dtype="float64"),
            "n_quant_peptides": pd.Series(count, dtype="int64"),
        }
    )


def _to_log2(values: pd.Series, scale: str) -> pd.Series:
    if scale == "log2":
        v = pd.to_numeric(values, errors="coerce").astype("float64")
        return v.where(np.isfinite(v))
    return pd.Series(np.log2(_positive(values).to_numpy()), index=values.index)


def resolve_c10_scale(
    abundance: pd.DataFrame, configured: str | None, c10_scale: str | None
) -> str:
    """Measurement scale of the C10 values (PR #82 lead decision 1): an explicit
    ``abundance.scale`` (not ``auto``) first, then the preflight scale (``c10_scale``,
    ``M4Inputs.abundance_scale``), else ValueError. C10 has no ``scale`` column; an optional one
    carrying more than one value raises when ``abundance.scale`` is ``auto``."""
    if configured is not None and configured != "auto":
        return configured
    if "scale" in abundance.columns:
        scales = set(abundance["scale"].dropna().astype(str))
        if len(scales) > 1:
            raise ValueError(f"C10 abundance: one measurement scale expected, got {sorted(scales)}")
    if c10_scale is not None and c10_scale != "auto":
        return c10_scale
    raise ValueError("C10 abundance: scale unknown (abundance.scale=auto and no preflight scale)")


def _flags_of(cell: Any) -> list[str]:
    if cell is None or (isinstance(cell, float) and math.isnan(cell)):
        return []
    if isinstance(cell, str):
        return [f for f in cell.split(";") if f]
    return [str(f) for f in cell]


def _cross_strain_proteins(protein_hits: pd.DataFrame | None) -> frozenset[str]:
    """Proteins with ``cross_strain_protein`` on any of their C12 hit rows (M3-16 flags)."""
    if protein_hits is None or protein_hits.empty or "flags" not in protein_hits.columns:
        return frozenset()
    return frozenset(
        str(p)
        for p, f in zip(protein_hits["protein_id"], protein_hits["flags"], strict=True)
        if "cross_strain_protein" in _flags_of(f)
    )


def protein_abundance(
    abundance: pd.DataFrame,
    pmap: pd.DataFrame,
    declared_type: str,
    scale: str,
    protein_hits: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """C10 values per gene through the protein -> gene map: ``locus_tag``, ``value_log2`` (the
    largest over the gene's proteins), ``value_raw`` (that protein's input value), ``n_proteins``
    and ``flags`` (``two_proteins_one_gene``, ``cross_strain_protein``). ``scale`` is the resolved
    measurement scale (``resolve_c10_scale``); ``cross_strain_protein`` is set when a protein mapped
    by ``protein_hit`` carries that flag in the C12 ``protein_hits``. +-inf values become NA."""
    if declared_type not in (*MOLAR_TYPES, "maxlfq", "intensity", "unknown"):
        raise ValueError(f"unknown abundance type {declared_type!r}")
    if scale not in ("log2", "linear"):
        raise ValueError(f"unknown C10 measurement scale {scale!r}")
    cross = _cross_strain_proteins(protein_hits)
    prot = pd.DataFrame(
        {
            "protein_id": abundance["protein_id"].astype(str).to_numpy(),
            "value_raw": pd.to_numeric(abundance["value"], errors="coerce").to_numpy("float64"),
        }
    )
    prot["value_log2"] = _to_log2(prot["value_raw"], scale).to_numpy()
    pm = pmap.loc[pmap["gene"].notna()]
    pm = pd.DataFrame(
        {
            "protein_id": pm["protein_id"].astype(str).to_numpy(),
            "locus_tag": pm["gene"].astype(str).to_numpy(),
            "cross": [
                str(m) == "protein_hit" and str(p) in cross
                for p, m in zip(pm["protein_id"], pm["method"], strict=True)
            ],
        }
    )
    joined = pm.merge(prot, on="protein_id", how="inner")
    rows: list[dict[str, Any]] = []
    for tag, grp in joined.groupby("locus_tag", sort=True):
        g = grp.sort_values(["value_log2", "protein_id"], ascending=[False, True])
        best = g.iloc[0]
        fl = ["two_proteins_one_gene"] if len(g) >= 2 else []
        if bool(g["cross"].any()):
            fl.append("cross_strain_protein")
        rows.append(
            {
                "locus_tag": str(tag),
                "value_log2": float(best["value_log2"]),
                "value_raw": float(best["value_raw"]),
                "n_proteins": len(g),
                "flags": fl,
            }
        )
    cols = ["locus_tag", "value_log2", "value_raw", "n_proteins", "flags"]
    out = pd.DataFrame(rows, columns=cols)
    return out.astype({"value_log2": "float64", "value_raw": "float64", "n_proteins": "int64"})


def _fill_statistic(observed: np.ndarray, policy: str) -> float | None:
    if policy not in POLICIES:
        raise ValueError(f"unknown missing-value policy {policy!r}")
    if policy == "exclude" or observed.size == 0:
        return None
    if policy == "median":
        return float(np.median(observed))
    if policy == "mean":
        return math.fsum(observed.tolist()) / observed.size
    lo = float(observed.min())
    return lo - 1.0 if policy == "zero" else lo


def impute(
    log2: pd.Series,
    detected: pd.Series,
    policy: AbundanceMissing,
    eligible: pd.Series | None = None,
) -> tuple[pd.Series, pd.Series, float | None]:
    """Fill missing ``log2`` of ``detected`` genes with the policy statistic over the observed
    values of detected ``eligible`` (non-partial; default all) genes (median, mean, min,
    zero = min - 1, exclude = no fill). Detected partial genes are filled too (PR #82 lead
    decision 5). Returns (filled log2, imputed mask, fill value); undetected genes are never
    filled. +-inf counts as missing."""
    values = pd.to_numeric(log2, errors="coerce").astype("float64")
    values = values.where(np.isfinite(values))
    det = detected.reindex(values.index).fillna(False).astype(bool)
    elig = (
        pd.Series(True, index=values.index)
        if eligible is None
        else eligible.reindex(values.index).fillna(True).astype(bool)
    )
    missing = values.isna()
    fill = _fill_statistic(values[det & elig & ~missing].to_numpy(), policy)
    mask = det & missing if fill is not None else pd.Series(False, index=values.index)
    filled = values.where(~mask, fill if fill is not None else math.nan)
    return filled, mask.astype(bool), fill


def _none_table(tags: Sequence[str]) -> pd.DataFrame:
    n = len(tags)
    return pd.DataFrame(
        {
            "locus_tag": pd.Series(list(tags), dtype=object),
            "abundance_raw": np.full(n, np.nan),
            "abundance_log2": np.full(n, np.nan),
            "abundance_imputed": np.zeros(n, dtype=bool),
            "imputation_policy": pd.Series([None] * n, dtype=object),
            "abundance_type": pd.Series(["none"] * n, dtype=object),
            "abundance_scale": pd.Series([None] * n, dtype=object),
            "top3_log2": np.full(n, np.nan),
            FLAGS_COLUMN: [[] for _ in range(n)],
        }
    )


def _has_intensity(accepted: pd.DataFrame) -> bool:
    if accepted.empty:
        return False
    keep = accepted["abundance_eligible"].astype(bool).to_numpy()
    return bool(pd.to_numeric(accepted.loc[keep, "intensity"], errors="coerce").notna().any())


def gene_abundance(
    genes: pd.DataFrame,
    accepted: pd.DataFrame,
    covariates: pd.DataFrame,
    detection: Any,
    abundance: pd.DataFrame | None,
    pmap: pd.DataFrame,
    cfg: GiapConfig,
    policy: AbundanceMissing | None = None,
    *,
    c10_scale: str | None = None,
    protein_hits: pd.DataFrame | None = None,
) -> AbundanceResult:
    """§9.3 abundance columns, one row per C2 gene in ``genes`` order. ``detection`` is the M4-12
    ``Detection`` (only ``.table``: locus_tag, detected is read); ``policy`` overrides
    ``abundance.missing`` (sensitivity reruns); ``c10_scale`` is the preflight C10 scale
    (``M4Inputs.abundance_scale``) and ``protein_hits`` the C12 hits (``cross_strain_protein``).
    Partial genes stay out of the fill statistic but, when detected, are filled."""
    conf = cfg.abundance
    pol: AbundanceMissing = policy if policy is not None else conf.missing
    tags = genes["locus_tag"].astype(str)
    idx = pd.Index(tags.tolist())
    n = len(idx)
    if _has_intensity(accepted):
        source: Source = "peptide_intensity"
        n_obs = covariates.set_index(covariates["locus_tag"].astype(str))["n_observable_peptides"]
        pa = peptide_abundance(accepted, n_obs).set_index("locus_tag").reindex(idx)
        raw = pa[conf.method].astype("float64")
        log2 = pd.Series(np.log2(raw.to_numpy()), index=idx)
        top3_log2 = pd.Series(np.log2(pa["top3"].to_numpy("float64")), index=idx)
        nq = pa["n_quant_peptides"].fillna(0).to_numpy()
        gene_flags: list[list[str]] = [["top3_lt3"] if 1 <= k <= 2 else [] for k in nq]
        a_type, a_scale = str(conf.method), "molar"
    elif abundance is not None:
        source = "protein_table"
        sc = resolve_c10_scale(abundance, conf.scale, c10_scale)
        pr = protein_abundance(abundance, pmap, conf.type, sc, protein_hits)
        pr = pr.set_index("locus_tag").reindex(idx)
        raw = pr["value_raw"].astype("float64")
        log2 = pr["value_log2"].astype("float64")
        top3_log2 = pd.Series(np.full(n, np.nan), index=idx)
        gene_flags = [_flags_of(f) for f in pr["flags"]]
        a_type = str(conf.type)
        a_scale = "molar" if a_type in MOLAR_TYPES else "relative"
    else:
        return AbundanceResult(_none_table(tags.tolist()), "none", None)

    det_tab = detection.table
    det = pd.Series(det_tab["detected"].astype(bool).to_numpy(), index=det_tab["locus_tag"])
    det = det.reindex(idx).fillna(False).astype(bool)
    if "partial" in genes.columns:
        elig = pd.Series(~genes["partial"].astype(bool).to_numpy(), index=idx)
    else:
        elig = pd.Series(True, index=idx)
    filled, mask, fill = impute(log2, det, pol, elig)
    if fill is None and pol != "exclude" and bool(det.any()):
        _log.warning("no_observed_abundance: no detected gene has an observed abundance")
    imputed = mask.to_numpy(bool)
    for fl, imp in zip(gene_flags, imputed, strict=True):
        if imp:
            fl.append("imputed_abundance")
    table = pd.DataFrame(
        {
            "locus_tag": pd.Series(tags.tolist(), dtype=object),
            "abundance_raw": raw.to_numpy("float64"),
            "abundance_log2": filled.to_numpy("float64"),
            "abundance_imputed": imputed,
            "imputation_policy": pd.Series([pol] * n, dtype=object),
            "abundance_type": pd.Series([a_type] * n, dtype=object),
            "abundance_scale": pd.Series([a_scale] * n, dtype=object),
            "top3_log2": top3_log2.to_numpy("float64"),
            FLAGS_COLUMN: gene_flags,
        }
    )
    return AbundanceResult(table, source, fill)
