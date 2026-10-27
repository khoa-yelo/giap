"""Attack nulls for the MAYU closed-form gene FDR (M4-02; GIPPA_SPEC.md §7.3.2, §8.5; statistician).

Generator (one replicate): G = 3,000 eligible genes with protein lengths log-normal (median
300 aa, log-SD ``sd``) and ``n_observable_peptides`` a_g = round(length / 25) floored at 1. A gene
is present with a logistic probability rising in log a_g (detectability tied to size, as in real
data, so absent genes are short). Present genes receive Poisson(0.3·a_g) true peptides (at least
one). Calibrated PEPs: as many decoy-like candidates as true peptides; every peptide is false with
probability equal to its PEP; false peptides land on genes with probability ∝ hit weight
(∝ a_g^gamma, or ∝ protein length in the ``noisy_size`` scenario). Peptide q = running mean of
sorted PEPs; accepted = q ≤ 0.01.
A gene is false when it is not present (all its accepted peptides are false).

Attacks inside the method's assumptions (false hits ∝ a_g or more concentrated); each must stay at
level, mean FDP ≤ q + k·SE (§8.4) at q in ``levels``:

- ``detectability_bias``: the base generator (absent genes are short);
- ``heavy_tail``: log-SD 1.2 (wide within-bin size spread; Jensen gap of the binned w_b);
- ``one_hit_wonders``: present genes carry ~1 true peptide (most genes are single-peptide genes);
- ``long_gene_concentration``: false hits ∝ a_g² (wrong hits pile up on long genes) — conservative;
- ``no_pep``: F(t) = t·|P(t)| instead of Σ PEP.

Assumption violations, direction documented (docs/methods/mayu-gene-fdr.md), not level gates:

- ``uniform_hits``: false hits uniform over genes (gamma = 0) → realised FDP > estimate
  (ratio ≈ 1.24);
- ``noisy_size``: hits ∝ protein length while a_g = Poisson(length/25) → mild underestimate
  (ratio ≈ 1.05).

R is ``replicates_fast.gene_fdr`` (marker ``calibration``) or ``replicates_full.gene_fdr``
(``calibration_full``); seed ``seeds.calibration_base``; k = ``se_multiplier``. The level gates go
through ``_criteria.calibrate`` (ADR-0015 two-block confirmation, sharding; stream
``gene_fdr.attack_mayu.<attack>``); the violation directions keep their single-stream estimate.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest
from _gate import calibrate_levels

from gippa.stats.mayu import mayu_gene_fdr

G = 3000
N_BINS = 10
ACCEPT_Q = 0.01


@dataclass(frozen=True)
class Scenario:
    sd: float = 0.6
    gamma: float = 1.0
    noisy_size: bool = False
    one_hit: bool = False
    use_pep: bool = True


@dataclass(frozen=True)
class Replicate:
    fdp: np.ndarray  # realised FDP at gene q <= each level
    fdr_loose: float  # estimated FDR at the loosest threshold
    fdp_loose: float  # realised FDP among all genes with an accepted peptide


def _replicate(rng: np.random.Generator, s: Scenario, levels: list[float]) -> Replicate:
    length = np.exp(rng.normal(np.log(300.0), s.sd, G))
    if s.noisy_size:
        a = rng.poisson(length / 25.0).astype(np.float64)
        hit_w = length
    else:
        a = np.maximum(1.0, np.round(length / 25.0))
        hit_w = a**s.gamma
    af = np.maximum(a, 1.0)
    logit = -0.5 + 1.5 * (np.log(af) - np.log(np.median(af)))
    present = rng.random(G) < 1.0 / (1.0 + np.exp(-logit))
    lam = np.full(G, 0.3) if s.one_hit else 0.3 * af
    n_true = np.where(present, np.maximum(rng.poisson(lam), 1), 0)
    true_gene = np.repeat(np.arange(G), n_true)
    n_t = true_gene.size
    pep = np.r_[rng.beta(0.3, 30.0, n_t), rng.beta(2.0, 2.0, n_t)]
    false = rng.random(pep.size) < pep
    gene = np.r_[true_gene, np.zeros(n_t, dtype=np.int64)]
    candidate_true = ~false & (np.arange(pep.size) >= n_t)
    gene[candidate_true] = rng.choice(np.flatnonzero(present), size=int(candidate_true.sum()))
    gene[false] = rng.choice(G, size=int(false.sum()), p=hit_w / hit_w.sum())

    order = np.argsort(pep, kind="stable")
    q = np.empty(pep.size)
    q[order] = np.minimum.accumulate((np.cumsum(pep[order]) / np.arange(1, pep.size + 1))[::-1])[
        ::-1
    ]
    acc = q <= ACCEPT_Q
    best = np.full(G, np.nan)
    np.fmin.at(best, gene[acc], q[acc])
    res = mayu_gene_fdr(best, q[acc], pep[acc] if s.use_pep else None, a, N_BINS)

    absent = ~present
    fdp = np.zeros(len(levels))
    for j, alpha in enumerate(levels):
        called = res.gene_q <= alpha
        n = int(called.sum())
        fdp[j] = (absent & called).sum() / n if n else 0.0
    any_hit = ~np.isnan(best)
    return Replicate(
        fdp=fdp,
        fdr_loose=float(res.curve["fdr"].iloc[-1]),
        fdp_loose=float((absent & any_hit).sum() / any_hit.sum()),
    )


AT_LEVEL: dict[str, Scenario] = {
    "detectability_bias": Scenario(),
    "heavy_tail": Scenario(sd=1.2),
    "one_hit_wonders": Scenario(one_hit=True),
    "long_gene_concentration": Scenario(gamma=2.0),
    "no_pep": Scenario(use_pep=False),
}

# documented direction: ratio band for mean realised FDP / mean estimated FDR (loosest threshold)
VIOLATIONS: dict[str, tuple[Scenario, float, float]] = {
    "uniform_hits": (Scenario(gamma=0.0), 1.10, 1.40),
    "noisy_size": (Scenario(noisy_size=True), 1.0, 1.15),
}

SUITES = [
    pytest.param("replicates_fast", marks=pytest.mark.calibration, id="fast"),
    pytest.param("replicates_full", marks=pytest.mark.calibration_full, id="full"),
]


def _settings(tolerances: dict[str, object], suite: str) -> tuple[int, int, float, list[float]]:
    reps, seeds = tolerances[suite], tolerances["seeds"]
    assert isinstance(reps, dict) and isinstance(seeds, dict)
    levels = [float(x) for x in tolerances["levels"]]  # type: ignore[union-attr]
    return (
        int(reps["gene_fdr"]),
        int(seeds["calibration_base"]),
        float(tolerances["se_multiplier"]),  # type: ignore[arg-type]
        levels,
    )


@pytest.mark.parametrize("suite", SUITES)
@pytest.mark.parametrize("attack", sorted(AT_LEVEL))
def test_mayu_attack_stays_at_level(tolerances: dict[str, object], attack: str, suite: str) -> None:
    r, _, _, levels = _settings(tolerances, suite)
    loose: list[tuple[float, float]] = []  # (estimated FDR, realised FDP), replicates run here

    def fdp(rng: np.random.Generator) -> np.ndarray:
        x = _replicate(rng, AT_LEVEL[attack], levels)
        loose.append((x.fdr_loose, x.fdp_loose))
        return x.fdp

    stream = f"gene_fdr.attack_mayu.{attack}"
    calibrate_levels("gene_fdr", "fdr", levels, fdp, r, stream=stream).assert_ok()
    if attack == "long_gene_concentration" and loose:  # documented direction: conservative
        est, real = np.mean(loose, axis=0)
        assert real < est, f"expected conservative: realised {real:.4f} >= estimated {est:.4f}"


@pytest.mark.parametrize("suite", SUITES)
@pytest.mark.parametrize("violation", sorted(VIOLATIONS))
def test_mayu_violation_direction(
    tolerances: dict[str, object], violation: str, suite: str
) -> None:
    r, seed, _, levels = _settings(tolerances, suite)
    scenario, lo, hi = VIOLATIONS[violation]
    rng = np.random.default_rng(seed)
    reps = [_replicate(rng, scenario, levels) for _ in range(r)]
    est = np.mean([x.fdr_loose for x in reps])
    real = np.mean([x.fdp_loose for x in reps])
    ratio = real / est
    assert lo < ratio < hi, (
        f"{violation}: realised/estimated = {ratio:.3f}, documented ({lo}, {hi})"
    )
