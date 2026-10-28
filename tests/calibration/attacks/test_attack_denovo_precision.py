"""Attack nulls for the de novo π̂ handed to the placement FDR (M3-11 → M3-10; GIPPA_SPEC.md
§7.2.5, §8.5; ADR-0010; statistician). Method: docs/methods/denovo-precision.md, sections
"π̂ for M3-10" and "Calibration and attacks".

Generator (one replicate, de novo, placement stratum m = 0): N peptides in 10 score bins of equal
size with true precision ``P_BIN`` (π_b = 1 - p_b), one length class. A correct peptide matches the
target exactly with probability ``c_true``; a wrong one with probability ``q_mult`` · r_b (chance
hit); every peptide matches the shuffled decoy with probability r_b (r_b = 0.1, or 0.25 → 0.05 from
the worst to the best bin when ``r_by_bin``: low-score bins hold short peptides with more chance
hits). ``cluster`` > 1 makes runs of peptides share correctness and both match draws (nested or
overlapping peptides of one locus: dependent rows). Target E: correct 10^N(-3.5, 1), chance
10^N(-1.5, 1); decoy E 10^N(-1.5, 1).
The estimator sees the same exact-match calls (as in the pipeline), with c_l = ``c_ref`` and
π̂ = ``pi_wrong_denovo`` (1 - ci_low). The FDR is ADR-0010 written out directly:
FDR(e) = (Σ π̂_j over decoy hits with E ≤ e + max π̂ over searched peptides) / T(e), q = running min.

At level (mean FDP ≤ q + k·SE at every q in ``levels``):
- ``base`` (N = 4000), ``small_bins`` (N = 200, 20 per bin), ``tiny_bins`` (N = 50, 5 per bin);
- ``chance_by_bin``: r_b tied to π_b (ADR-0010's hardest heterogeneity);
- ``clustered``: rows dependent in runs of 5 (bootstrap treats them as independent);
- ``easy_reference``: reference peptides match more often than correct de novo peptides
  (c_ref 0.99 > c_true 0.90) → precision biased DOWN, π̂ up: conservative;
- ``c_true_one``: c_ref 0.99 (the ``c_assumed`` default) while correct peptides always match.

Assumption violations, direction documented (not level gates):
- ``hard_reference``: c_ref 0.85 < c_true 0.99 (reference set harder than the searched set, or
  c_l from a handful of reference peptides) → precision biased UP, π̂ too low → FDP above q at
  q = 0.01 (≈ 0.012);
- ``near_correct_wrong``: wrong peptides match the target 2x as often as the decoy (r_b
  underestimates the wrong-peptide target rate) → FDP ≈ 2q; fails with the oracle π too, so it is
  a decoy-model violation owned by M3-09/M3-10/M3-21, not by the π̂ bound;
- ``small_bin_bound``: per-bin coverage of the lower bound P(1 - ci_low ≥ π) in a 20-peptide
  bin with p = 0.98 is ≈ 0.5, far below 0.975 (percentile bootstrap at the clip): the bound alone
  is not conservative there; the placement FDR holds only through the + max π̂ term (``tiny_bins``,
  ``small_bins`` above).

R is ``replicates_fast.placement`` (marker ``calibration``) or ``replicates_full.placement``
(``calibration_full``); seed ``seeds.calibration_base``; k = ``se_multiplier``. The at-level gates
go through ``_criteria.calibrate`` (ADR-0015 two-block confirmation, sharding; stream
``placement.attack_denovo_pi.<attack>``, one generator per replicate for the simulation and the
bootstrap); the documented violations keep their single-stream estimates.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels

from gippa.mapping.denovo_precision import denovo_precision, pi_wrong_denovo, score_bins
from gippa.rng import for_step

P_BIN = np.array([0.3, 0.45, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98, 0.99])
N_BIN = P_BIN.size
CLASS = "8"


@dataclass(frozen=True)
class Scenario:
    n: int = 4000
    c_ref: float = 0.99
    c_true: float = 0.99
    r_by_bin: bool = False
    q_mult: float = 1.0
    cluster: int = 1
    oracle: bool = False


AT_LEVEL = {
    "base": Scenario(),
    "small_bins": Scenario(n=200),
    "tiny_bins": Scenario(n=50),
    "chance_by_bin": Scenario(r_by_bin=True),
    "clustered": Scenario(cluster=5),
    "easy_reference": Scenario(c_true=0.90),
    "c_true_one": Scenario(c_true=1.0),
}


def _simulate(rng: np.random.Generator, s: Scenario) -> pd.DataFrame:
    n = s.n // N_BIN
    n_cl = -(-n // s.cluster)
    parts = []
    for b, p in enumerate(P_BIN):
        r = 0.25 - 0.2 * b / (N_BIN - 1) if s.r_by_bin else 0.1

        def rep(x: np.ndarray, n: int = n) -> np.ndarray:
            return np.repeat(x, s.cluster)[:n]

        correct = rep(rng.random(n_cl) < p)
        u_t, u_d = rep(rng.random(n_cl)), rep(rng.random(n_cl))
        target = np.where(correct, u_t < s.c_true, u_t < min(1.0, s.q_mult * r))
        e_t = np.where(correct, 10 ** rng.normal(-3.5, 1.0, n), 10 ** rng.normal(-1.5, 1.0, n))
        parts.append(
            pd.DataFrame(
                {
                    "peptide": "PEPTIDEK",
                    "score": b + rng.uniform(0.01, 0.99, n),
                    "length_class": CLASS,
                    "target_exact": target,
                    "decoy_exact": u_d < r,
                    "correct": correct,
                    "e_t": e_t,
                    "e_d": 10 ** rng.normal(-1.5, 1.0, n),
                    "pi_true": 1.0 - p,
                }
            )
        )
    return pd.concat(parts, ignore_index=True)


def _pi_hat(df: pd.DataFrame, s: Scenario, rng: np.random.Generator) -> np.ndarray:
    if s.oracle:
        return np.asarray(df["pi_true"].to_numpy(), dtype=float)
    bins = score_bins(df["score"].to_numpy())
    rows = denovo_precision(df, {CLASS: s.c_ref}, False, rng)
    return pi_wrong_denovo(rows, bins)


def _fdp(df: pd.DataFrame, pi: np.ndarray, levels: list[float]) -> np.ndarray:
    """Realised FDP at each q under ADR-0010: FDR(e) = (W(e) + max π̂) / T(e)."""
    t = df["target_exact"].to_numpy()
    d = df["decoy_exact"].to_numpy()
    e_t = df["e_t"].to_numpy()[t]
    e_d = df["e_d"].to_numpy()[d]
    order = np.argsort(e_t, kind="stable")
    s = e_t[order]
    n_t = np.searchsorted(s, s, side="right").astype(float)
    od = np.argsort(e_d, kind="stable")
    w_cum = np.concatenate([[0.0], np.cumsum(pi[d][od])])
    w = w_cum[np.searchsorted(e_d[od], s, side="right")]
    fdr = (w + float(pi.max())) / n_t
    q = np.minimum(np.minimum.accumulate(fdr[::-1])[::-1], 1.0)
    wrong = ~df["correct"].to_numpy()[t][order]
    out = np.zeros(len(levels))
    for j, alpha in enumerate(levels):
        acc = q <= alpha
        out[j] = wrong[acc].sum() / acc.sum() if acc.any() else 0.0
    return out


def _run(s: Scenario, r: int, seed: int, levels: list[float]) -> tuple[np.ndarray, np.ndarray]:
    sim = for_step("denovo_precision_attack", seed)
    fdp = np.zeros((r, len(levels)))
    for i in range(r):
        df = _simulate(sim, s)
        fdp[i] = _fdp(df, _pi_hat(df, s, for_step("placement", seed + i)), levels)
    return fdp.mean(axis=0), fdp.std(axis=0, ddof=1) / np.sqrt(r)


SUITES = [
    pytest.param("replicates_fast", marks=pytest.mark.calibration, id="fast"),
    pytest.param("replicates_full", marks=pytest.mark.calibration_full, id="full"),
]


def _settings(tolerances: dict[str, object], suite: str) -> tuple[int, int, float, list[float]]:
    reps, seeds = tolerances[suite], tolerances["seeds"]
    assert isinstance(reps, dict) and isinstance(seeds, dict)
    lv = tolerances["levels"]
    assert isinstance(lv, list)
    levels = [float(x) for x in lv]
    k = float(tolerances["se_multiplier"])  # type: ignore[arg-type]
    return int(reps["placement"]), int(seeds["calibration_base"]), k, levels


@pytest.mark.parametrize("suite", SUITES)
@pytest.mark.parametrize("attack", sorted(AT_LEVEL))
def test_denovo_pi_attack_stays_at_level(
    tolerances: dict[str, object], attack: str, suite: str
) -> None:
    r, _, _, levels = _settings(tolerances, suite)
    s = AT_LEVEL[attack]

    def fdp(rng: np.random.Generator) -> np.ndarray:
        df = _simulate(rng, s)
        return _fdp(df, _pi_hat(df, s, rng), levels)

    stream = f"placement.attack_denovo_pi.{attack}"
    calibrate_levels("placement", "fdr", levels, fdp, r, stream=stream).assert_ok()


@pytest.mark.parametrize("suite", SUITES)
def test_hard_reference_is_anticonservative(tolerances: dict[str, object], suite: str) -> None:
    """c_ref < c_true biases precision up: documented failure direction (FDP > q at q = 0.01)."""
    r, seed, k, levels = _settings(tolerances, suite)
    mean, se = _run(Scenario(c_ref=0.85), r, seed, levels)
    j = levels.index(0.01)
    assert mean[j] > 0.01 + k * se[j], f"expected FDP > 0.01, got {mean[j]:.4f}"


@pytest.mark.parametrize("suite", SUITES)
def test_near_correct_wrong_breaks_even_oracle(tolerances: dict[str, object], suite: str) -> None:
    """Wrong peptides hitting the target 2x the decoy rate: FDP ≈ 2q, with π̂ and with oracle π."""
    r, seed, k, levels = _settings(tolerances, suite)
    for oracle in (False, True):
        mean, se = _run(Scenario(q_mult=2.0, oracle=oracle), r, seed, levels)
        j = levels.index(0.05)
        assert mean[j] > 0.05 + k * se[j], f"oracle={oracle}: FDP {mean[j]:.4f} not above 0.05"


@pytest.mark.parametrize("suite", SUITES)
def test_small_bin_lower_bound_undercovers(tolerances: dict[str, object], suite: str) -> None:
    """Per-bin coverage of 1 - ci_low ≥ π in 20-peptide bins at p = 0.98 is far below 0.975."""
    r, seed, _, _ = _settings(tolerances, suite)
    sim = for_step("denovo_precision_attack", seed)
    s = Scenario(n=200)
    top = int(np.flatnonzero(P_BIN == 0.98)[0])
    hit = []
    for i in range(r):
        df = _simulate(sim, s)
        pi = _pi_hat(df, s, for_step("placement", seed + i))
        sel = np.floor(df["score"].to_numpy()).astype(int) == top
        hit.append(bool(pi[sel].min() >= 1.0 - P_BIN[top]))
    assert np.mean(hit) < 0.9, f"coverage {np.mean(hit):.3f}: bound unexpectedly conservative"
