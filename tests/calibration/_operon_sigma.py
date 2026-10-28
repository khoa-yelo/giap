"""Shared helpers of the M6-13 calibration (docs/specs/M6-13-discordance-abundance-sigma.md
"Calibration and the gate"; the M6-09 harness in ``_operon_gate.py``).

Every operon_null arm (D-1 with the log-linear null, its teeth, D-1b-step, D-1b-floor, D-1c and
the null quartile check) reads **genome i of the D-1 stream**
(``suite_seed("discordance.attack_d1.noise_by_abundance", i)``), the genomes the statistician
measured (Slurm 1885387 / 1885429 / 1885466), whatever stream name its harness records carry. A
replicate rng of stream s is mapped back to its index i through ``index_of(s, R)``, so a check's
log key (stream, criterion, level, R) stays unique per arm while the genomes are shared.

- ``noise_by_shape(pairs, shape)``: Delta of the both-detected OP pairs rescaled to
  ``loc + (Delta - loc) * exp(h(z))``, z the standardised pair-mean abundance, with
  ``h = -0.5 z`` (D-1), the centred ratio-2 step at the median (D-1b-step) or the centred floor
  ``0.5 log(0.15 + e^-z)`` (D-1b-floor). The existing ``_noise_by_abundance`` and M-2 are
  untouched.
- ``gene_noise(sim, seed)``: D-1c, per-gene measurement noise
  ``obs_g = true_g + 0.35 exp(-0.5 z_g) eps_g``, z_g the standardised true abundance of the
  detected genes, eps from ``default_rng([seed, 7])``; the pair columns are rebuilt from the genes.
- ``quartile_rates``: per-quartile hits / n pooled over replicates (quartiles cut per genome,
  group 0 = ``<= q25``, as the D-1 test).
"""

from __future__ import annotations

import json
import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from typing import Any

import _criteria as C
import _operon_gate as G
import numpy as np
import pandas as pd

D1_STREAM = "discordance.attack_d1.noise_by_abundance"
DISC_STREAM = "discordance.gate.operon_null_1800"
DI = "discordance"
#: a quartile rate above this multiple of alpha fails the test (S-9; the attack-C pattern). Not a
#: level bound and not a tolerance
ESCALATE_QUARTILE = 2.0
M_POOLED = 200
D1C_KAPPA, D1C_SLOPE = 0.35, -0.5
SHAPES = ("d1", "step", "floor", "none")


@dataclass(frozen=True)
class Base:
    """The parts of one simulated genome the arms read (the full simulation is not cached)."""

    seed: int
    pairs: pd.DataFrame
    genes: pd.DataFrame  # locus_tag, abundance_log2_true, detected_true
    truth_discordant: frozenset[str]


@cache
def base(seed: int, scenario: tuple[Any, ...] = ("operon_null", G.N_GENES)) -> Base:
    sim = G.simulate(G.scenario(*scenario), int(seed))
    g = sim.expression.genes[["locus_tag", "abundance_log2_true", "detected_true"]].copy()
    disc = frozenset(sim.truth["discordant_pairs"]["pair_id"].astype(str))
    return Base(int(seed), sim.pairs, g, disc)


@cache
def index_of(stream: str, R: int) -> Mapping[int, int]:
    """replicate seed -> index, both blocks (0 .. 2R - 1)."""
    return {C.suite_seed(stream, i): i for i in range(2 * R)}


def d1_seed(stream: str, R: int, rng: np.random.Generator) -> int:
    """The D-1 genome seed of the replicate whose rng is ``rng`` on ``stream``."""
    i = index_of(stream, R)[int(rng.bit_generator.seed_seq.entropy)]
    return C.suite_seed(D1_STREAM, i)


def _h(shape: str, z: np.ndarray) -> np.ndarray:
    if shape == "none":
        return np.zeros_like(z)
    if shape == "d1":
        return -0.5 * z
    if shape == "step":
        h = np.where(z < np.median(z), math.log(2) / 2, -math.log(2) / 2)
    elif shape == "floor":
        h = 0.5 * np.log(0.15 + np.exp(-z))
    else:
        raise ValueError(f"unknown shape {shape!r}")
    return h - h.mean()


def noise_by_shape(pairs: pd.DataFrame, shape: str) -> pd.DataFrame:
    """D-1 / D-1b: OP pairs with both abundances finite, z on their pair mean."""
    from gippasim.scenarios import OPERON_T

    loc = OPERON_T[0]
    out = pairs.copy()
    a = out["abundance_log2_a"].astype(float).to_numpy()
    b = out["abundance_log2_b"].astype(float).to_numpy()
    keep = (out["class"] == "OP").to_numpy() & np.isfinite(a) & np.isfinite(b)
    m = (a + b) / 2
    h = np.zeros(len(out))
    z = (m[keep] - m[keep].mean()) / m[keep].std(ddof=1)
    h[keep] = _h(shape, z)
    d = b - a
    new = np.where(keep, loc + (d - loc) * np.exp(h), d)
    out["abundance_log2_b"] = np.where(keep, a + new, b)
    out["delta_log2"] = np.where(keep, new, out["delta_log2"].astype(float))
    return out


def gene_noise(bs: Base) -> tuple[pd.DataFrame, np.ndarray]:
    """D-1c pairs and the true upstream abundance of every pair (row order of the pairs)."""
    rng = np.random.default_rng([bs.seed, 7])
    g = bs.genes
    tags = g["locus_tag"].astype(str).to_numpy()
    tru = g["abundance_log2_true"].astype(float).to_numpy()
    det = g["detected_true"].astype(bool).to_numpy()
    zg = np.zeros(len(g))
    zg[det] = (tru[det] - tru[det].mean()) / tru[det].std(ddof=1)
    eps = rng.standard_normal(len(g))
    obs = np.where(det, tru + D1C_KAPPA * np.exp(D1C_SLOPE * zg) * eps, np.nan)
    pos = {t: k for k, t in enumerate(tags)}
    p = bs.pairs.copy()
    ia = np.array([pos[str(t)] for t in p["gene_a"]])
    ib = np.array([pos[str(t)] for t in p["gene_b"]])
    p["abundance_log2_a"], p["abundance_log2_b"] = obs[ia], obs[ib]
    p["delta_log2"] = obs[ib] - obs[ia]
    return p, tru[ia]


@dataclass(frozen=True)
class SigmaRun:
    p: np.ndarray  # both_detected_t p, row order
    q: np.ndarray
    pair_id: np.ndarray
    strata: Mapping[str, np.ndarray]  # covariates of the tested pairs, by name
    scaled: Any  # Discordance.scaled
    fit: Any  # Discordance.fit (constant)


def sigma_run(pairs: pd.DataFrame, genes: pd.DataFrame, sigma_model: str,
              strata: Mapping[str, pd.Series]) -> SigmaRun:  # fmt: skip
    """``discordance(..., sigma_model=sigma_model)``; ``strata`` maps a name to a Series indexed
    by ``pair_id`` (read at the tested pairs)."""
    from gippa.config import GippaConfig
    from gippa.operons.discordance import discordance

    res = discordance(pairs, genes[["locus_tag"]], GippaConfig(), None, sigma_model=sigma_model)
    t = res.rows[res.rows["test"] == "both_detected_t"]
    ids = t["pair_id"].astype(str).to_numpy()
    return SigmaRun(
        p=t["p_value"].astype(float).to_numpy(),
        q=t["q_value"].astype(float).to_numpy(),
        pair_id=ids,
        strata={k: v.loc[ids].to_numpy(float) for k, v in strata.items()},
        scaled=res.scaled,
        fit=res.fit,
    )


def _by_id(pairs: pd.DataFrame, values: np.ndarray) -> pd.Series:
    return pd.Series(values, index=pairs["pair_id"].astype(str).to_numpy())


def planted_strata(pairs: pd.DataFrame) -> dict[str, pd.Series]:
    """m: the planted (pre-transform) pair mean; a: the upstream gene's abundance."""
    a = pairs["abundance_log2_a"].astype(float).to_numpy()
    b = pairs["abundance_log2_b"].astype(float).to_numpy()
    return {"m": _by_id(pairs, (a + b) / 2), "a": _by_id(pairs, a)}


@cache
def shape_run(seed: int, shape: str, sigma_model: str) -> SigmaRun:
    bs = base(seed)
    noisy = bs.pairs if shape == "none" else noise_by_shape(bs.pairs, shape)
    return sigma_run(noisy, bs.genes, sigma_model, planted_strata(bs.pairs))


@cache
def d1c_run(seed: int, sigma_model: str) -> SigmaRun:
    bs = base(seed)
    pairs, a_true = gene_noise(bs)
    return sigma_run(pairs, bs.genes, sigma_model, {"a_true": _by_id(pairs, a_true)})


def quartile_rates(runs: Sequence[SigmaRun], stratum: str,
                   levels: Sequence[float]) -> dict[str, Any]:  # fmt: skip
    """Rates per quartile of ``stratum`` (cut per genome), pooled over ``runs``."""
    hits = np.zeros((4, len(levels)))
    n = np.zeros(4, dtype=int)
    for r in runs:
        s = r.strata[stratum]
        grp = np.searchsorted(np.quantile(s, [0.25, 0.5, 0.75]), s, side="left")
        for q in range(4):
            sel = grp == q
            n[q] += int(sel.sum())
            hits[q] += [(r.p[sel] <= a).sum() for a in levels]
    rate = hits / n[:, None]
    return {"stratum": stratum, "R": len(runs), "n": n.tolist(), "rate": rate.tolist(),
            "ratio": (rate / np.asarray(levels)).tolist(),
            "bounds": [[C.level_bound(a, int(k)) for a in levels] for k in n]}  # fmt: skip


def escalation(rec: Mapping[str, Any], levels: Sequence[float]) -> list[str]:
    """Quartile x level cells above ``ESCALATE_QUARTILE * alpha``."""
    return [f"Q{q + 1} alpha={a:g}: {rec['rate'][q][j]:.4f} > {ESCALATE_QUARTILE:g} alpha"
            for q in range(4) for j, a in enumerate(levels)
            if rec["rate"][q][j] > ESCALATE_QUARTILE * a]  # fmt: skip


def gamma_record(runs: Sequence[SigmaRun]) -> dict[str, Any]:
    """gamma-hat, gamma-hat / s, the LR and the share of each model / reason (recorded)."""
    s = [r.scaled for r in runs]
    g = np.array([x.gamma for x in s], dtype=float)
    per = np.array([x.gamma / x.sd if x.sd else math.nan for x in s], dtype=float)
    lr = np.array([x.lr_stat for x in s], dtype=float)
    reasons: dict[str, int] = {}
    for x in s:
        key = f"{x.model}:{x.reason}"
        reasons[key] = reasons.get(key, 0) + 1
    return {"R": len(s), "gamma_mean": float(np.nanmean(g)),
            "gamma_sd": float(np.nanstd(g, ddof=1)), "gamma_per_log2_mean": float(np.nanmean(per)),
            "lr_mean": float(np.nanmean(lr)), "models": reasons}  # fmt: skip


def pooled_rep(get_run: Callable[[np.random.Generator], SigmaRun],
               levels: Sequence[float]) -> Callable[[np.random.Generator], np.ndarray]:  # fmt: skip
    """The M6-09 pooled replicate: ``M_POOLED`` tested pairs drawn at random per genome."""

    def rep(rng: np.random.Generator) -> np.ndarray:
        run = get_run(rng)
        assert run.p.size >= M_POOLED, run.p.size
        p = rng.choice(run.p, M_POOLED, replace=False)
        return np.array([float((p <= a).mean()) for a in levels])

    return rep


def log(name: str, row: Mapping[str, Any]) -> None:
    path = C.log_dir() / f"operons_sigma_{name}.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(dict(row)) + "\n")
