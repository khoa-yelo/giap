"""Attack nulls for the M4-03 permutation engine (GIAP_SPEC.md §8.5; statistician-owned).

Each generator is built to push ``perm_p``, the circular/linear shift null or ``besag_clifford``
above level while staying inside the assumptions in docs/methods/permutation-and-shifts.md:

- ``ar1_circular_exhaustive`` / ``ar1_circular_sampled``: a circular-stationary AR(1)-like process
  (rho = 0.9) on a circular replicon; one contiguous region; exhaustive shifts, or B = 199 shifts
  drawn without replacement (``n > n_max_exhaustive``). Shifts must stay at level.
  ``test_naive_permutation_is_anticonservative`` shows the attack has teeth: permuting gene labels
  on the same data ignores the autocorrelation and is far above level.
- ``tiny_circular_n10`` / ``tiny_circular_n20`` / ``tiny_linear_n12``: small plasmids and short
  linear contigs, where K = n - 1 or n - span is 9-19 and the two-sided centring matters. Each
  replicate is a draft genome of ``CONTIGS`` independent replicons, so the level check counts
  R * CONTIGS independent p-values.
- ``heavy_ties``: discrete values (70 % at one imputed minimum) on non-dyadic decimals, so equal
  multisets summed in different orders differ by float noise; the tie tolerance must keep the
  one-sided p at least as large as the exact rational p (``test_tie_tolerance_never_lowers_p``).
- ``besag_clifford_ties``: two groups of Bernoulli(0.1) indicators (heavily tied difference in
  means), sequential p with h = 10.

The AR(1) and heavy-ties attacks simulate ``REPLICONS`` independent replicons per replicate, so
their level check also counts R * REPLICONS independent p-values.

Criterion (§8.4 level): for alpha in ``levels``, rejection rate <= alpha + k sqrt(alpha(1-alpha)/M),
M independent p-values. R = ``replicates_fast.regions`` (marker ``calibration``) or
``replicates_full.regions`` (``calibration_full``); seed ``seeds.calibration_base``. Attacks with
one p per replicate (``ONE_P``) are gated by ``_criteria.calibrate`` (ADR-0015 two-block
confirmation, sharding; stream ``regions.attack_permutation.<attack>``); the pooled ones keep
their own check over M (calibrate takes one scalar per replicate).
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pytest
from _gate import calibrate_levels

from giap.rng import for_step
from giap.stats.permutation import (
    apply_shift,
    besag_clifford,
    choose_shifts,
    circular_shifts,
    perm_p,
)

Replicate = Callable[[np.random.Generator], np.ndarray]
CONTIGS = 50
REPLICONS = 10
N_MAX_EXHAUSTIVE = 20_000


def _circular_ar1(rng: np.random.Generator, n: int, rho: float = 0.9) -> np.ndarray:
    """Stationary AR(1)-like process on a circle: white noise circularly filtered by rho^j."""
    kernel = np.zeros(n)
    j = np.arange(min(n, 400))
    kernel[j] = rho**j
    e = rng.normal(size=n)
    return np.asarray(np.fft.irfft(np.fft.rfft(e) * np.fft.rfft(kernel), n), dtype=np.float64)


def _window_sums(x: np.ndarray, w: int) -> np.ndarray:
    """Sum of x over the circular window starting at each position."""
    n = x.size
    cs = np.concatenate([[0.0], np.cumsum(np.concatenate([x, x[: w - 1]]))])
    return cs[w : w + n] - cs[:n]


def _shift_p(
    x: np.ndarray, start: int, w: int, circular: bool, n_max: int, b: int, rng: np.random.Generator
) -> float:
    """Two-sided shift p for the contiguous region [start, start + w) (circular indices)."""
    n = x.size
    idx = (start + np.arange(w)) % n
    shifts, _ = choose_shifts(circular_shifts(n, idx, circular), n_max, n, b, rng)
    sums = _window_sums(x, w)
    starts = np.array([apply_shift([start], int(k), n)[0] for k in shifts], dtype=np.int64)
    return float(perm_p(float(x[idx].sum()), sums[starts]))


def _ar1_circular(n: int, n_max: int) -> Replicate:
    def rep(rng: np.random.Generator) -> np.ndarray:
        out = np.empty(REPLICONS)
        for c in range(REPLICONS):
            x = _circular_ar1(rng, n)
            out[c] = _shift_p(x, int(rng.integers(0, n)), 20, True, n_max, 199, rng)
        return out

    return rep


def _tiny(n: int, w: int, circular: bool) -> Replicate:
    def rep(rng: np.random.Generator) -> np.ndarray:
        out = np.empty(CONTIGS)
        for c in range(CONTIGS):
            x = rng.normal(size=n)
            start = int(rng.integers(0, n if circular else n - w + 1))
            out[c] = _shift_p(x, start, w, circular, N_MAX_EXHAUSTIVE, 1_000, rng)
        return out

    return rep


def _tied_values(rng: np.random.Generator, n: int) -> np.ndarray:
    levels = np.array([201, 237, 263, 299], dtype=np.int64)  # tenths: 20.1 ... 29.9 (log2-like)
    return rng.choice(levels, size=n, p=[0.7, 0.1, 0.1, 0.1])


def _heavy_ties(rng: np.random.Generator) -> np.ndarray:
    n, w = 300, 10
    out = np.empty(REPLICONS)
    for c in range(REPLICONS):
        x = _tied_values(rng, n) / 10.0
        out[c] = _shift_p(x, int(rng.integers(0, n)), w, True, N_MAX_EXHAUSTIVE, 1_000, rng)
    return out


def _bc_ties(rng: np.random.Generator) -> np.ndarray:
    labels = np.arange(60) < 30
    x = (rng.random(60) < 0.1).astype(float)

    def diff(lab: np.ndarray) -> float:
        return abs(float(x[lab].mean() - x[~lab].mean()))  # permutation mean is exactly 0

    p, _ = besag_clifford(lambda g: diff(g.permutation(labels)), diff(labels), 10, 1_000, rng)
    return np.array([p])


ATTACKS: dict[str, Replicate] = {
    "ar1_circular_exhaustive": _ar1_circular(1000, N_MAX_EXHAUSTIVE),
    "ar1_circular_sampled": _ar1_circular(3000, 1000),
    "tiny_circular_n10": _tiny(10, 2, circular=True),
    "tiny_circular_n20": _tiny(20, 3, circular=True),
    "tiny_linear_n12": _tiny(12, 2, circular=False),
    "heavy_ties": _heavy_ties,
    "besag_clifford_ties": _bc_ties,
}

#: attacks with one p per replicate. The others pool CONTIGS / REPLICONS independent p per
#: replicate (M = R x CONTIGS p-values per check); ``calibrate`` takes one scalar per replicate
#: with R = the suite's R (shard ranges), so they keep their own level check over M.
ONE_P = frozenset({"besag_clifford_ties"})


def _one_p(rep: Replicate) -> Callable[[np.random.Generator], float]:
    def fn(rng: np.random.Generator) -> float:
        (p,) = rep(rng)
        assert np.isfinite(p) and 0 < p <= 1, p
        return float(p)

    return fn


SUITES = [
    pytest.param("replicates_fast", marks=pytest.mark.calibration, id="fast"),
    pytest.param("replicates_full", marks=pytest.mark.calibration_full, id="full"),
]


def _settings(tolerances: dict[str, object], suite: str) -> tuple[list[float], int, int, float]:
    lv, reps, seeds, se_k = (tolerances[key] for key in ("levels", suite, "seeds", "se_multiplier"))
    assert isinstance(lv, list) and isinstance(reps, dict) and isinstance(seeds, dict)
    assert isinstance(se_k, (int, float))
    levels, k = [float(x) for x in lv], float(se_k)
    return levels, int(reps["regions"]), int(seeds["calibration_base"]), k


@pytest.mark.parametrize("suite", SUITES)
@pytest.mark.parametrize("attack", sorted(ATTACKS))
def test_attack_stays_at_level(tolerances: dict[str, object], attack: str, suite: str) -> None:
    levels, r, seed, k = _settings(tolerances, suite)
    if attack in ONE_P:  # one p per replicate: gated by calibrate (ADR-0015, sharding)
        stream = f"regions.attack_permutation.{attack}"
        fn = _one_p(ATTACKS[attack])
        calibrate_levels("regions", "level", levels, fn, r, stream=stream).assert_ok()
        return
    rng = for_step(f"attack.permutation.{attack}.{suite}", seed)
    p = np.concatenate([ATTACKS[attack](rng) for _ in range(r)])
    assert np.isfinite(p).all() and ((p > 0) & (p <= 1)).all()
    m = p.size
    for a in levels:
        rate = float((p <= a).mean())
        bound = a + k * np.sqrt(a * (1 - a) / m)
        assert rate <= bound, f"{attack} alpha={a}: rate {rate:.4f} > bound {bound:.4f}, M={m}"


@pytest.mark.parametrize("suite", SUITES)
def test_naive_permutation_is_anticonservative(tolerances: dict[str, object], suite: str) -> None:
    """Teeth check: gene-label permutation on the AR(1) attack data is far above level."""
    _, r, seed, k = _settings(tolerances, suite)
    rng = for_step(f"attack.permutation.naive_ar1.{suite}", seed)
    n, w, a = 1000, 20, 0.05
    p = np.empty(r)
    for i in range(r):
        x = _circular_ar1(rng, n)
        idx = (int(rng.integers(0, n)) + np.arange(w)) % n
        null = np.array([x[rng.choice(n, w, replace=False)].sum() for _ in range(199)])
        p[i] = perm_p(float(x[idx].sum()), null)
    rate = float((p <= a).mean())
    assert rate > a + k * np.sqrt(a * (1 - a) / r), f"naive permutation rate {rate:.3f} at level"


@pytest.mark.parametrize("suite", SUITES)
def test_tie_tolerance_never_lowers_p(tolerances: dict[str, object], suite: str) -> None:
    """With exact ties summed in different orders, float p >= exact rational p ("greater")."""
    _, r, seed, _ = _settings(tolerances, suite)
    rng = for_step(f"attack.permutation.tie_exact.{suite}", seed)
    n, w = 300, 50
    for _ in range(r):
        xi = _tied_values(rng, n)
        x = xi / 10.0
        idx = (int(rng.integers(0, n)) + np.arange(w)) % n
        shifts = circular_shifts(n, idx, circular=True)
        moved = [apply_shift(idx, int(s), n) for s in shifts]
        null = np.array([x[j].sum() for j in moved])
        exact_obs = int(xi[idx].sum())
        b = sum(int(xi[j].sum()) >= exact_obs for j in moved)
        exact = (1 + b) / (1 + len(shifts))
        assert perm_p(float(x[idx].sum()), null, "greater") >= exact - 1e-15
