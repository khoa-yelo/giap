"""Monte Carlo reproduction of the B.2 chance-hit table (issue M3-09; GIAP_SPEC.md §7.2.5,
§7.2.9, Appendix B.2: "A Monte Carlo simulation on random sequences with these frequencies MUST
reproduce the <= 1 column for L = 8-10 within 20 %").

A random target of N = 10^6 residues is drawn from a 19-letter frequency vector with
sum f^2 = 0.069 (E. coli composition with I merged into L, shrunk toward uniform to hit 0.069
exactly); random peptides of length L are drawn from the same vector. Hits with <= 1
substitution per peptide are counted exactly (one sorted key array per wildcard position:
count(<= 1) = sum_j count(match except j) - (L - 1) * count(exact)), and the mean, scaled to
N = 10^7, must be within 20 % of the B.2 <= 1 column and of ``expected_hits_uniform``, with a
Monte Carlo SE <= 5 %. Seeds: ``giap.rng.for_step(..., seeds.calibration_base)``.

The full suite (L = 8, 9, 10) takes ~10 s; a fast L = 8 version runs on every commit.
"""

from __future__ import annotations

import numpy as np
import pytest

from giap.mapping.chance import ResidueFreqs, expected_hits, expected_hits_uniform
from giap.mapping.sixframe import ALPHABET
from giap.rng import for_step

S = 0.069
B2_LE1 = {8: 0.56, 9: 0.043, 10: 0.0033}
N_PEPTIDES = {8: 200_000, 9: 1_000_000, 10: 4_000_000}
RES = ALPHABET[:19]
# E. coli K-12 amino-acid composition (%), I merged into L.
ECOLI = {
    "A": 9.5,
    "C": 1.1,
    "D": 5.1,
    "E": 5.8,
    "F": 3.9,
    "G": 7.4,
    "H": 2.3,
    "K": 4.4,
    "L": 16.7,
    "M": 2.8,
    "N": 3.9,
    "P": 4.4,
    "Q": 4.4,
    "R": 5.5,
    "S": 5.8,
    "T": 5.4,
    "V": 7.1,
    "W": 1.5,
    "Y": 2.8,
}


def freqs_with_s(s: float = S) -> np.ndarray:
    """(1 - t) * uniform + t * E. coli with t chosen so that sum f^2 = s exactly."""
    b = np.array([ECOLI[a] for a in RES])
    b /= b.sum()
    u = np.full(19, 1 / 19)
    t = np.sqrt((s - 1 / 19) / ((b**2).sum() - 1 / 19))
    f = (1 - t) * u + t * b
    assert abs((f**2).sum() - s) < 1e-12 and f.min() > 0
    return f


def _count(sorted_keys: np.ndarray, q: np.ndarray) -> np.ndarray:
    return np.searchsorted(sorted_keys, q, "right") - np.searchsorted(sorted_keys, q, "left")


def hits_le1(target: np.ndarray, peptides: np.ndarray) -> np.ndarray:
    """Exact number of target windows within <= 1 substitution of each peptide row."""
    L = peptides.shape[1]
    n = len(target) - L + 1
    w = 19 ** np.arange(L, dtype=np.int64)
    tk = np.zeros(n, np.int64)
    for i in range(L):
        tk += target[i : i + n].astype(np.int64) * w[i]
    pk = (peptides.astype(np.int64) * w).sum(axis=1)
    total = -(L - 1) * _count(np.sort(tk), pk)
    for j in range(L):
        masked_t = np.sort(tk - (tk // w[j] % 19) * w[j])
        total += _count(masked_t, pk - peptides[:, j].astype(np.int64) * w[j])
    return total


def simulate(L: int, n_peptides: int, n_target: int, seed: int) -> tuple[float, float]:
    """(mean hits per peptide scaled to N = 1e7, its SE), deterministic in ``seed``."""
    rng = for_step(f"calibration.chance_hits.L{L}", seed)
    f = freqs_with_s()
    target = rng.choice(19, size=n_target, p=f).astype(np.uint8)
    peptides = rng.choice(19, size=(n_peptides, L), p=f).astype(np.uint8)
    h = hits_le1(target, peptides)
    scale = 1e7 / (n_target - L + 1)
    return float(h.mean() * scale), float(h.std(ddof=1) / np.sqrt(n_peptides) * scale)


def _seed(tolerances: dict[str, object]) -> int:
    seeds = tolerances["seeds"]
    assert isinstance(seeds, dict)
    return int(seeds["calibration_base"])


def test_counting_kernel_matches_brute_force() -> None:
    rng = np.random.default_rng(0)
    target = rng.integers(0, 3, size=400).astype(np.uint8)  # small alphabet: many hits
    peptides = rng.integers(0, 3, size=(50, 5)).astype(np.uint8)
    win = np.lib.stride_tricks.sliding_window_view(target, 5)
    brute = np.array([((win != p).sum(axis=1) <= 1).sum() for p in peptides])
    np.testing.assert_array_equal(hits_le1(target, peptides), brute)


def test_simulation_is_deterministic() -> None:
    assert simulate(8, 2000, 20_000, 1) == simulate(8, 2000, 20_000, 1)


def test_fast_b2_le1_l8(tolerances: dict[str, object]) -> None:
    """Fast version (every commit): L = 8, N = 2e5 target, 1e5 peptides."""
    mean, se = simulate(8, 100_000, 200_000, _seed(tolerances))
    model = expected_hits_uniform(8, 1, S, 1e7)
    assert se <= 0.08 * mean
    assert abs(mean - B2_LE1[8]) <= 0.2 * B2_LE1[8], (mean, se)
    assert abs(mean - model) <= 0.2 * model, (mean, model)


@pytest.mark.calibration
@pytest.mark.parametrize("L", [8, 9, 10])
def test_b2_le1_column_within_20pct(L: int, tolerances: dict[str, object]) -> None:
    mean, se = simulate(L, N_PEPTIDES[L], 10**6, _seed(tolerances))
    model = expected_hits_uniform(L, 1, S, 1e7)
    assert se <= 0.05 * mean, f"L={L}: MC SE {se:.3g} > 5 % of {mean:.3g}"
    assert abs(mean - B2_LE1[L]) <= 0.2 * B2_LE1[L], f"L={L}: MC {mean:.4g} vs B.2 {B2_LE1[L]}"
    assert abs(mean - model) <= 0.2 * model, f"L={L}: MC {mean:.4g} vs E(L,1) {model:.4g}"


@pytest.mark.calibration
def test_ranking_score_mean_matches_simulation(tolerances: dict[str, object]) -> None:
    """Per-peptide E_i averaged over random peptides matches the simulated hit rate (L = 8)."""
    L, n = 8, 200_000
    mean, se = simulate(L, n, 10**6, _seed(tolerances))
    rng = for_step(f"calibration.chance_hits.L{L}", _seed(tolerances))
    f = freqs_with_s()
    rng.choice(19, size=10**6, p=f)  # same target draw as simulate()
    peptides = rng.choice(19, size=(n, L), p=f)
    fr = ResidueFreqs(f=f, n_positions=10**7)
    sub = peptides[:5000]
    e_mean = float(np.mean([expected_hits("".join(RES[c] for c in row), 1, fr) for row in sub]))
    assert abs(e_mean - mean) <= 0.2 * mean, (e_mean, mean, se)
