"""Inputs for the M4-13 detectability tests: small synthetic gene sets and giapsim (M4-07 layout +
M4-08 ``covariates_only`` detection) data sets with M4-10 covariates.

Nothing here imports ``giap.genetable.detectability``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from scipy.special import expit

NAMES: tuple[str, ...] = ("n_observable_peptides", "length_aa", "gravy", "pi", "tm_segments")
COV_COLUMNS = ("length_aa", "n_observable_peptides", "gravy", "pi", "mw_da", "tm_segments")


@dataclass
class GeneSet:
    detected: pd.Series  # "boolean" dtype, index locus_tag; <NA> = not eligible
    cov: pd.DataFrame  # index locus_tag
    tu_of_gene: pd.Series  # index locus_tag -> TU id
    p_true: pd.Series = field(default_factory=lambda: pd.Series(dtype=float))
    extra: dict[str, Any] = field(default_factory=dict)


def synthetic(
    n_tus: int,
    seed: int,
    *,
    intercept: float = -7.0,
    tm_values: tuple[int, ...] = (0, 1, 2),
    detected: np.ndarray | None = None,
) -> GeneSet:
    """1-3 genes per TU, covariates with realistic ranges, detection logistic in log length,
    log1p n_obs and GRAVY (``detected`` overrides the Bernoulli draw)."""
    r = np.random.default_rng(seed)
    rows = []
    i = 0
    for j in range(n_tus):
        for _ in range(int(r.integers(1, 4))):
            i += 1
            length = int(np.exp(r.normal(5.5, 0.6))) + 10
            rows.append(
                {
                    "locus_tag": f"G{i:05d}",
                    "tu": f"TU_{j + 1}",
                    "length_aa": length,
                    "n_observable_peptides": int(r.poisson(length / 25)),
                    "gravy": float(r.normal(-0.3, 0.4)),
                    "pi": float(r.uniform(4, 11)),
                    "mw_da": float(length * 110.0),
                    "tm_segments": int(r.choice(tm_values)),
                }
            )
    df = pd.DataFrame(rows).set_index("locus_tag")
    eta = (
        intercept
        + 1.1 * np.log(df["length_aa"])
        + 0.6 * np.log1p(df["n_observable_peptides"])
        - 0.8 * df["gravy"]
    )
    p = expit(eta.to_numpy())
    d = r.random(len(df)) < p if detected is None else np.asarray(detected, bool)
    return GeneSet(
        detected=pd.Series(pd.array(d, dtype="boolean"), index=df.index, name="detected"),
        cov=df[list(COV_COLUMNS)].copy(),
        tu_of_gene=df["tu"].rename("tu_id"),
        p_true=pd.Series(p, index=df.index),
    )


def with_counts(n_det: int, n_undet: int, seed: int, genes_per_tu: int = 2) -> GeneSet:
    """Exactly ``n_det`` detected and ``n_undet`` undetected genes, ``genes_per_tu`` per TU."""
    n = n_det + n_undet
    r = np.random.default_rng(seed)
    d = np.zeros(n, bool)
    d[r.choice(n, n_det, replace=False)] = True
    tags = [f"G{i + 1:05d}" for i in range(n)]
    length = (np.exp(r.normal(5.5, 0.6, n)) + 10).astype(int)
    cov = pd.DataFrame(
        {
            "length_aa": length,
            "n_observable_peptides": r.poisson(length / 25),
            "gravy": r.normal(-0.3, 0.4, n),
            "pi": r.uniform(4, 11, n),
            "mw_da": length * 110.0,
            "tm_segments": r.choice([0, 1, 2], n),
        },
        index=tags,
    )
    tu = pd.Series([f"TU_{i // genes_per_tu + 1}" for i in range(n)], index=tags, name="tu_id")
    return GeneSet(pd.Series(pd.array(d, dtype="boolean"), index=tags), cov, tu)


def sim_geneset(
    rng: np.random.Generator,
    n_genes: int,
    *,
    short_genes_in: str | None = None,
) -> GeneSet:
    """M4-07 layout (``n_genes_chromosome = n_genes``) + M4-08 ``covariates_only`` detection.

    Eligible = non-partial ``cds``/``sorf`` rows (giapsim gives pseudogenes an intercept-only
    detection probability, so they are marked not eligible here); covariates from M4-10
    (``covariates(aa, "trypsin_p", (7, 30))``, the giapsim ``n_obs`` definition).
    """
    import giapsim.expression as ex
    import giapsim.layout as la
    from giap.genetable.covariates import covariates

    lay = la.simulate_layout(
        la.LayoutConfig(n_genes_chromosome=n_genes, short_genes_in=short_genes_in), rng
    )
    expr = ex.simulate_expression(lay, ex.ExpressionConfig(detection="covariates_only"), rng)
    e = expr.genes.set_index("locus_tag")
    g = lay.genes.set_index("locus_tag")
    eligible = g["feature_type"].isin(["cds", "sorf"]).to_numpy() & ~g["partial"].to_numpy()
    rows = []
    for aa, ok in zip(g["aa"], eligible, strict=True):
        if ok and isinstance(aa, str) and aa:
            c = covariates(aa, "trypsin_p", (7, 30))
            rows.append({k: getattr(c, k) for k in COV_COLUMNS})
        else:
            rows.append({k: np.nan for k in COV_COLUMNS})
    cov = pd.DataFrame(rows, index=g.index)
    d = pd.array(e.loc[g.index, "detected_true"].to_numpy(bool), dtype="boolean")
    d[~eligible] = pd.NA
    tu = lay.truth.set_index("locus_tag")["tu_id"].loc[g.index]
    in_prophage = set(
        lay.region_genes.merge(lay.regions[["region_id", "region_type"]], on="region_id").query(
            "region_type == 'prophage'"
        )["locus_tag"]
    )
    return GeneSet(
        detected=pd.Series(d, index=g.index, name="detected"),
        cov=cov,
        tu_of_gene=tu,
        p_true=e.loc[g.index, "p_detect_true"].astype(float),
        extra={"in_prophage": g.index.isin(list(in_prophage)), "eligible": eligible},
    )
