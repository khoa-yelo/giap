"""Gene-frame generators for the M8-04b ORA calibration and attack arms (spec M8-04 decision 12).

The reference is the statistician's ``$GIAP_ROOT/stat/m8-04a-review/attack_tusize.py`` (Slurm
1867453/1867454; ADR-0035), rewritten to emit what ``detection_enrichment`` consumes: a gene frame,
``GeneSets`` and M4-13-style TU folds.

- ``simulate(rng, c=, a=, b=, oracle=)``: G genes in TUs of size min(Geometric(0.45), 8); gene
  logit = eta + c * u_TU with eta = -0.6 + N(0, 1) per gene and u_TU ~ N(0, 1) shared by the TU
  (the TU-packed effect); detection ~ Bernoulli(p_true). ``p_detect`` is the exact marginal
  E_u[expit(eta + c u)] (40-point Gauss-Hermite, as the reference), then miscalibrated as
  expit(a + b * logit(.)) (arm 4c: a = +-0.5); ``oracle=True`` gives p_true itself.
- ``whole_tu_sets(sim, rng, kind, n, size)``: ``n`` disjoint categories of ``size`` whole TUs
  (every gene of the TU is a member) drawn from ``kind`` in {random, large (>= 4 genes),
  single (1 gene)}; fewer when the pool is short (as the reference).
- ``random_gene_sets(sim, rng, n, size)``: ``n`` disjoint categories of ``size`` random genes
  (arm 1/6 style; TU membership then follows the any-eligible-gene rule).
- ``folds``: ``giap.genetable.detectability.tu_folds`` over the TU IDs, k = 5 (M4-13's default),
  drawn from the replicate generator.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.polynomial.hermite_e import hermegauss
from scipy.special import expit, logit

from giap.function.collections import GeneSets
from giap.genetable.detectability import tu_folds

G = 3000
SET_TUS = 15
K_FOLDS = 5
_GH_X, _GH_W = hermegauss(40)
_GH_W = _GH_W / _GH_W.sum()


@dataclass(frozen=True)
class Sim:
    genes: pd.DataFrame
    sizes: np.ndarray  # genes per TU, TU i = TU_<i+1>
    tu_of_gene: np.ndarray  # TU index per gene, gene j = g<j+1>
    d_tu: np.ndarray  # TU detected
    folds: Mapping[str, int]

    @property
    def tu_ids(self) -> list[str]:
        return [f"TU_{i + 1}" for i in range(self.sizes.size)]

    def tu_genes(self, i: int) -> list[str]:
        return [f"g{j + 1}" for j in np.flatnonzero(self.tu_of_gene == i)]


def simulate(
    rng: np.random.Generator,
    *,
    n_genes: int = G,
    c: float = 0.0,
    a: float = 0.0,
    b: float = 1.0,
    oracle: bool = False,
    eta0: float = -0.6,
) -> Sim:
    sizes_l: list[int] = []
    while sum(sizes_l) < n_genes:
        sizes_l.append(min(int(rng.geometric(0.45)), 8))
    sizes = np.array(sizes_l)
    T = sizes.size
    tu = np.repeat(np.arange(T), sizes)
    g = tu.size
    eta = eta0 + rng.normal(0.0, 1.0, g)
    u = rng.normal(0.0, 1.0, T)
    p_true = expit(eta + c * u[tu])
    if oracle:
        p_hat = p_true
    else:
        p_marg = (expit(eta[:, None] + c * _GH_X[None, :]) * _GH_W).sum(1)
        p_hat = expit(a + b * logit(p_marg))
    det = rng.random(g) < p_true
    pos = np.concatenate([np.arange(1, s + 1) for s in sizes])
    genes = pd.DataFrame(
        {
            "locus_tag": [f"g{j + 1}" for j in range(g)],
            "tu_id": [f"TU_{i + 1}" for i in tu],
            "tu_position": pos,
            "detected": det,
            "p_detect": p_hat.astype(float),
            "partial": np.zeros(g, dtype=bool),
            "abundance_log2": np.where(det, 20.0, np.nan),
            "abundance_imputed": np.zeros(g, dtype=bool),
            "abundance_scale": "molar",
        }
    )
    d_tu = np.bincount(tu, weights=det.astype(float), minlength=T) > 0
    folds = tu_folds([f"TU_{i + 1}" for i in range(T)], K_FOLDS, rng)
    return Sim(genes, sizes, tu, d_tu, folds)


def _pool(sim: Sim, kind: str) -> np.ndarray:
    if kind == "random":
        return np.arange(sim.sizes.size)
    if kind == "large":
        return np.flatnonzero(sim.sizes >= 4)
    if kind == "single":
        return np.flatnonzero(sim.sizes == 1)
    raise ValueError(kind)


def whole_tu_sets(
    sim: Sim, rng: np.random.Generator, kind: str, n: int, size: int = SET_TUS
) -> dict[str, frozenset[str]]:
    pool = rng.permutation(_pool(sim, kind))
    n = min(n, pool.size // size)
    out: dict[str, frozenset[str]] = {}
    for i in range(n):
        tus = pool[i * size : (i + 1) * size]
        out[f"{kind}{i:03d}"] = frozenset(t for x in tus for t in sim.tu_genes(int(x)))
    return out


def random_gene_sets(
    sim: Sim, rng: np.random.Generator, n: int, size: int = 40
) -> dict[str, frozenset[str]]:
    perm = rng.permutation(len(sim.genes))
    tags = sim.genes["locus_tag"].to_numpy()
    return {f"genes{i:03d}": frozenset(tags[perm[i * size : (i + 1) * size]]) for i in range(n)}


def collection(sets: Mapping[str, frozenset[str]], name: str = "kegg") -> GeneSets:
    s = {t: frozenset(g) for t, g in sorted(sets.items())}
    return GeneSets(name, s, {t: t for t in s}, {}, 1.0, 0)  # type: ignore[arg-type]
