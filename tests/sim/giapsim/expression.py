"""Expression and detection truth for sequence-free layouts (issue M4-08; GIAP_SPEC.md §11.3, §8.4).

TU level (log2): ``tu_log2 = mu + replicon copy + region effect + ori-ter gradient + island +
AR(1)``, one value per TU; the AR(1) runs over the TUs of each replicon in ``TU_<n>`` order with a
stationary start (``x_0 ~ N(0, s^2/(1 - rho^2))``, ``x_t = rho x_{t-1} + e_t``, ``s = tu_sd``).
Gene level: ``abundance_log2_true = tu_log2 - decay * position + stoichiometry + noise`` (+ the
gene's own region effect when its TU straddles a region boundary). Detection:
``p = sigmoid(b0 + b1 log len + b2 log1p n_obs + b3 GRAVY [+ b4 abundance] + region shift)``,
``detected_true ~ Bernoulli(p)``. Planted shifts and islands are deterministic, so runs with the
same seed and different shifts differ by exactly the shift. Nothing is imported from ``giap``.

Replichores (M5-04, GIAP_SPEC.md §7.3.5(e)): with a chr1 ``oric`` feature at midpoint o, a point m
has clockwise distance d = (m - o) mod L; ``right`` with position d / (L/2) when d <= L/2, else
``left`` with (L - d) / (L/2) (``NA`` off chr1 or without an oriC). The ori-ter gradient is the TU
shift g(p) at the TU midpoint's position p: cosine (a/2) cos(pi p), linear a (0.5 - p), with
``a = oriter_log2``; it is deterministic, and ``oriter_log2 = 0`` adds nothing at all.
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from itertools import repeat
from typing import Any, Literal

import numpy as np
import pandas as pd

__all__ = [
    "DEFAULT_BETA",
    "Expression",
    "ExpressionConfig",
    "FrozenMap",
    "IslandSpec",
    "detection_terms",
    "expression_for_simulation",
    "simple_covariates",
    "simulate_expression",
]

#: Kyte & Doolittle (1982) hydropathy; any other letter scores 0.
KD: Mapping[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}  # fmt: skip
_TRYPSIN_P = re.compile(r"(?<=[KR])")
_PEPTIDE_LEN = (7, 30)
_TU_NUM = re.compile(r"(\d+)$")
ISLAND_COLUMNS = (
    "island_id", "replicon", "first_tu", "n_tus", "effect_sd", "effect_log2", "first_tu_id",
    "last_tu_id",
)  # fmt: skip


class FrozenMap(Mapping[str, float]):
    """An immutable, hashable ``str -> float`` mapping (keeps the config frozen and hashable)."""

    def __init__(self, items: Mapping[str, float] | Sequence[tuple[str, float]] = ()) -> None:
        pairs = items.items() if isinstance(items, Mapping) else items
        self._items = tuple(sorted((str(k), float(v)) for k, v in pairs))
        self._d = dict(self._items)

    def __getitem__(self, key: str) -> float:
        return self._d[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._d)

    def __len__(self) -> int:
        return len(self._d)

    def __hash__(self) -> int:
        return hash(self._items)

    def __repr__(self) -> str:
        return f"FrozenMap({self._d!r})"


#: Detection logit coefficients (lead Q1). The intercept gives ~55 % detected on the default
#: layout; ``log2_abundance`` is used only with ``detection = "abundance"`` (set it together with
#: an intercept offset of ``-b4 * mu_log2``).
DEFAULT_BETA: Mapping[str, float] = FrozenMap(
    {"intercept": -8.65, "log_length": 1.2, "log1p_n_obs": 0.8, "gravy": -0.6,
     "log2_abundance": 0.0}
)  # fmt: skip


@dataclass(frozen=True)
class IslandSpec:
    """A planted island: ``n_tus`` TUs from 0-based ``first_tu`` (replicon TU order), shifted by
    ``effect`` stationary SDs of the AR(1) background (``tu_sd`` at rho = 0)."""

    replicon: str
    first_tu: int
    n_tus: int
    effect: float


@dataclass(frozen=True)
class ExpressionConfig:
    """Knobs of the expression model (spec M4-08 "Interface"); all shifts are log2."""

    mu_log2: float = 20.0
    tu_sd: float = 2.0  # AR(1) innovation SD
    gene_sd: float = 0.5
    stoich_sd: float = 0.7
    ar1_rho: float = 0.5
    position_decay_log2: float = 0.2
    plasmid_copy_log2: Mapping[str, float] = field(default_factory=FrozenMap)
    region_effect_log2: Mapping[str, float] = field(default_factory=FrozenMap)
    region_detect_shift: Mapping[str, float] = field(default_factory=FrozenMap)
    islands: tuple[IslandSpec, ...] = ()
    detection: Literal["covariates_only", "abundance"] = "covariates_only"
    beta: Mapping[str, float] = DEFAULT_BETA
    oriter_log2: float = 0.0  # ori minus ter difference of TU-level log2 abundance (M5-04)
    oriter_shape: Literal["cosine", "linear"] = "cosine"

    def __post_init__(self) -> None:
        for name in ("plasmid_copy_log2", "region_effect_log2", "region_detect_shift", "beta"):
            object.__setattr__(self, name, FrozenMap(getattr(self, name)))
        object.__setattr__(self, "islands", tuple(self.islands))
        if not 0.0 <= self.ar1_rho < 1.0:
            raise ValueError(f"ar1_rho must be in [0, 1), got {self.ar1_rho}")
        if self.oriter_shape not in ("cosine", "linear"):
            raise ValueError(f"unknown oriter_shape {self.oriter_shape!r}")
        if not math.isfinite(self.oriter_log2):
            raise ValueError(f"oriter_log2 must be finite, got {self.oriter_log2}")
        if self.detection not in ("covariates_only", "abundance"):
            raise ValueError(f"unknown detection mode {self.detection!r}")
        missing = set(DEFAULT_BETA) - set(self.beta)
        if missing:
            raise ValueError(f"beta lacks {sorted(missing)}")
        if min(self.tu_sd, self.gene_sd, self.stoich_sd) < 0:
            raise ValueError("standard deviations must be >= 0")


@dataclass(frozen=True)
class Expression:
    """``genes``: one row per layout gene (layout order); ``islands``: planted-island truth."""

    genes: pd.DataFrame
    islands: pd.DataFrame


def simple_covariates(aa: str) -> tuple[int, float]:
    """(distinct Trypsin/P peptides of 7-30 aa, no missed cleavage; Kyte-Doolittle GRAVY)."""
    lo, hi = _PEPTIDE_LEN
    n = len({p for p in _TRYPSIN_P.split(aa) if lo <= len(p) <= hi})
    gravy = sum(map(KD.get, aa, repeat(0.0))) / len(aa) if aa else math.nan
    return n, gravy


def _tu_order(tus: pd.DataFrame) -> pd.DataFrame:
    """TUs in ``TU_<n>`` order (contig order, then start); frame order when ids are not numbered."""
    nums = [_TU_NUM.search(str(t)) for t in tus["tu_id"]]
    if all(nums):
        key = [int(m.group(1)) for m in nums if m]
        return tus.iloc[np.argsort(key, kind="stable")].reset_index(drop=True)
    return tus.reset_index(drop=True)


def _ar1(e: np.ndarray, rho: float, sd: float) -> np.ndarray:
    x = np.empty_like(e)
    if len(e):
        x[0] = sd * e[0] / math.sqrt(1.0 - rho * rho)
        for t in range(1, len(e)):
            x[t] = rho * x[t - 1] + sd * e[t]
    return x


def _replichore(mid: np.ndarray, oric_mid: float, length: float) -> tuple[np.ndarray, np.ndarray]:
    """(right arm?, replichore position in [0, 1]) of midpoints ``mid`` (any real; taken mod L)."""
    d = np.mod(np.mod(mid, length) - oric_mid, length)
    half = length / 2.0
    right = d <= half
    return right, np.where(right, d, length - d) / half


def _gradient(p: np.ndarray, cfg: ExpressionConfig) -> np.ndarray:
    a = cfg.oriter_log2
    if cfg.oriter_shape == "cosine":
        out: np.ndarray = a / 2.0 * np.cos(np.pi * p)
    else:
        out = a * (0.5 - p)
    return out


def _abundance(
    genes: pd.DataFrame,
    tus: pd.DataFrame,
    gene_types: Mapping[str, frozenset[str]],
    cfg: ExpressionConfig,
    rng: np.random.Generator,
    oric: tuple[float, float] | None = None,
) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    """TU- and gene-level abundance in ``genes`` order, and the island truth frame; ``oric`` is
    (chr1 oriC midpoint, chr1 length), needed when ``cfg.oriter_log2 != 0``."""
    tus = _tu_order(tus)
    contigs = tus["contig"].astype(str).to_numpy()
    members = [[str(g) for g in gs] for gs in tus["genes"]]
    unknown = set(cfg.plasmid_copy_log2) - set(contigs)
    if unknown:
        raise ValueError(f"plasmid_copy_log2 names unknown replicons {sorted(unknown)}")
    eps = rng.standard_normal(len(tus))
    x = np.zeros(len(tus))
    pos_in: dict[str, list[int]] = {}
    for i, c in enumerate(contigs):
        pos_in.setdefault(c, []).append(i)
    for idx in pos_in.values():
        x[idx] = _ar1(eps[idx], cfg.ar1_rho, cfg.tu_sd)
    shift = np.array([cfg.plasmid_copy_log2.get(c, 0.0) for c in contigs])
    gene_eff = {
        g: sum(cfg.region_effect_log2.get(t, 0.0) for t in gene_types.get(g, frozenset()))
        for g in genes["locus_tag"].astype(str)
    }
    region_tu = np.zeros(len(tus))  # effect of the region types shared by every TU member
    for i, gs in enumerate(members):
        if gs:
            shared = frozenset.intersection(*(gene_types.get(g, frozenset()) for g in gs))
            region_tu[i] = sum(cfg.region_effect_log2.get(t, 0.0) for t in shared)
    shift += region_tu
    if cfg.oriter_log2 != 0.0:
        if oric is None:
            raise ValueError("oriter_log2 != 0 needs an oriC on chr1 (LayoutConfig.oric_frac)")
        on = contigs == "chr1"
        mid = (tus["start"].to_numpy(float)[on] + tus["end"].to_numpy(float)[on]) / 2.0
        shift[on] += _gradient(_replichore(mid, *oric)[1], cfg)
    island_tu = np.zeros(len(tus), dtype=bool)
    island_of = np.full(len(tus), None, dtype=object)  # island id per TU, None outside islands
    rows = []
    stationary = cfg.tu_sd / math.sqrt(1.0 - cfg.ar1_rho**2)
    for k, isl in enumerate(cfg.islands):
        idx = pos_in.get(isl.replicon, [])
        if isl.n_tus < 1 or isl.first_tu < 0 or isl.first_tu + isl.n_tus > len(idx):
            raise ValueError(f"island {isl} outside replicon {isl.replicon!r} ({len(idx)} TUs)")
        sel = idx[isl.first_tu : isl.first_tu + isl.n_tus]
        shift[sel] += isl.effect * stationary
        island_tu[sel] = True
        island_of[sel] = f"island_{k + 1}"
        rows.append((f"island_{k + 1}", isl.replicon, isl.first_tu, isl.n_tus, float(isl.effect),
                     float(isl.effect * stationary), str(tus["tu_id"].iloc[sel[0]]),
                     str(tus["tu_id"].iloc[sel[-1]])))  # fmt: skip
    tu_val = cfg.mu_log2 + shift + x
    where = {g: (i, k) for i, gs in enumerate(members) for k, g in enumerate(gs)}
    tags = genes["locus_tag"].astype(str).tolist()
    missing = [g for g in tags if g not in where]
    if missing:
        raise ValueError(f"{len(missing)} genes are in no TU (first: {missing[0]})")
    tu_i = np.array([where[g][0] for g in tags], dtype=int)
    pos = np.array([where[g][1] for g in tags], dtype=float)
    tu_g = tu_val[tu_i]
    eff = np.array([gene_eff[g] for g in tags])
    own = eff - region_tu[tu_i]  # nonzero only in TUs straddling a region boundary
    stoich = cfg.stoich_sd * rng.standard_normal(len(tags))
    noise = cfg.gene_sd * rng.standard_normal(len(tags))
    ab = tu_g - cfg.position_decay_log2 * pos + stoich + noise + own
    islands = pd.DataFrame(rows, columns=list(ISLAND_COLUMNS))
    return {"tu_log2": tu_g, "abundance": ab, "in_island": island_tu[tu_i],
            "island_id": island_of[tu_i], "region_effect": eff}, islands  # fmt: skip


def _gene_region_types(layout: Any) -> dict[str, frozenset[str]]:
    typ = dict(zip(layout.regions["region_id"], layout.regions["region_type"], strict=True))
    out: dict[str, set[str]] = {}
    for rid, g in zip(layout.region_genes["region_id"], layout.region_genes["locus_tag"],
                      strict=True):  # fmt: skip
        out.setdefault(str(g), set()).add(str(typ[rid]))
    return {g: frozenset(t) for g, t in out.items()}


def _sigmoid(x: np.ndarray) -> np.ndarray:
    out: np.ndarray = 1.0 / (1.0 + np.exp(-x))
    return out


def _chr1_oric(layout: Any) -> tuple[float, float] | None:
    """(oriC midpoint mod L, L) of the layout's chr1 ``oric`` feature; None without one."""
    f = getattr(layout, "features", None)
    if f is None or len(f) == 0:
        return None
    f = f[(f["type"] == "oric") & (f["contig"] == "chr1")]
    if len(f) == 0:
        return None
    reps = layout.replicons
    length = float(reps.loc[reps["contig"] == "chr1", "length"].iloc[0])
    mid = (float(f["start"].iloc[0]) + float(f["end"].iloc[0])) / 2.0
    return mid % length, length


def detection_terms(
    layout: Any, cfg: ExpressionConfig
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(covariate logit ``b0 + b1 log len + b2 log1p n_obs + b3 GRAVY``, region shift, n_obs,
    GRAVY) per gene in ``layout.genes`` order; the abundance term is added by the caller."""
    genes = layout.genes
    types = _gene_region_types(layout)
    b = cfg.beta
    n_obs = np.zeros(len(genes), dtype=np.int64)
    gravy = np.full(len(genes), math.nan)
    logit = np.full(len(genes), b["intercept"])
    for i, aa in enumerate(genes["aa"]):
        if isinstance(aa, str) and aa:
            n_obs[i], gravy[i] = simple_covariates(aa)
            logit[i] += (b["log_length"] * math.log(len(aa)) + b["log1p_n_obs"]
                         * math.log1p(n_obs[i]) + b["gravy"] * gravy[i])  # fmt: skip
    tags = genes["locus_tag"].astype(str).tolist()
    shift = np.array([sum(cfg.region_detect_shift.get(t, 0.0) for t in types.get(g, ()))
                      for g in tags], dtype=float)  # fmt: skip
    return logit, shift, n_obs, gravy


def simulate_expression(layout: Any, cfg: ExpressionConfig, rng: np.random.Generator) -> Expression:
    """Expression and detection truth for a ``giapsim.layout.Layout`` (not mutated)."""
    genes = layout.genes
    types = _gene_region_types(layout)
    oric = _chr1_oric(layout)
    parts, islands = _abundance(genes, layout.tus, types, cfg, rng, oric)
    logit, shift, n_obs, gravy = detection_terms(layout, cfg)
    if cfg.detection == "abundance":
        logit += cfg.beta["log2_abundance"] * parts["abundance"]
    tags = genes["locus_tag"].astype(str).tolist()
    logit += shift
    p = _sigmoid(logit)
    detected = rng.random(len(tags)) < p
    frame = pd.DataFrame(
        {
            "locus_tag": tags,
            "tu_log2": parts["tu_log2"],
            "abundance_log2_true": parts["abundance"],
            "p_detect_true": p,
            "detected_true": detected,
            "n_obs_true": n_obs,
            "gravy_true": gravy,
            "in_island": parts["in_island"],
            "island_id": pd.array(parts["island_id"], dtype="string"),  # NA outside islands
            "region_effect_log2": parts["region_effect"].astype(float),
        }
    )
    arm = np.full(len(tags), None, dtype=object)
    pos = np.full(len(tags), math.nan)
    if oric is not None:
        on = (genes["contig"].astype(str) == "chr1").to_numpy()
        mid = (genes["start"].to_numpy(float)[on] + genes["end"].to_numpy(float)[on]) / 2.0
        right, p_on = _replichore(mid, *oric)
        arm[on] = np.where(right, "right", "left")
        pos[on] = p_on
    frame["replichore_true"] = pd.array(arm, dtype="string")  # NA off chr1 / without an oriC
    frame["replichore_position_true"] = pos
    return Expression(frame, islands)


def expression_for_simulation(
    genes: pd.DataFrame, tus: pd.DataFrame, cfg: ExpressionConfig, rng: np.random.Generator
) -> np.ndarray:
    """Per-gene natural-log abundance (``genes`` order) for the nucleotide simulator: the M3-19
    ``gene_log_abundance`` replacement (opt-in; region effects and the ori-ter gradient need a
    layout, so none here: ``oriter_log2 != 0`` is a ``ValueError``)."""
    parts, _ = _abundance(genes, tus, {}, cfg, rng)
    out: np.ndarray = parts["abundance"] * math.log(2.0)
    return out
