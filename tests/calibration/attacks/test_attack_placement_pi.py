"""Attack nulls for per-peptide pi_wrong in the placement FDR (M3-10; GIAP_SPEC.md §7.2.5, §8.5;
statistician). Decision and derivation: docs/methods/placement-fdr.md, section
"pi_wrong in the stratum FDR".

Signed rule (de novo stratum 0, where pi_wrong varies by score bin): in one stratum list ordered by
E, FDR(e) = (sum of pi_j over decoy hits j with E <= e + c) / T(e), c = max pi over all searched
peptides of the decoy type. With constant pi this is the spec formula pi (D(e) + 1) / T(e).

Generator (one replicate, de novo, stratum 0): N peptides in score bins with true wrong rates
``rates``; ``pi_denovo`` = rate x ``mult`` (capped at 1). Every peptide is searched against the
decoy and gets a chance decoy hit with probability p_chance(pi) and E from the chance law; a wrong
peptide gets a chance target hit independently from the same law (decoy hit <=> chance target hit
of the same peptide, exchangeable); a correct peptide always hits its target with low E.
p_chance is 0.25, or 0.05 + 0.4 pi ("corr": bad bins hit more often by chance, e.g. short peptides).

Attacks the signed rule must survive (mean FDP <= q + k SE, q in ``levels``):
- ``two_bins``: pi = 0.02 / 0.9 (extreme heterogeneity);
- ``two_bins_corr``, ``ten_bins_corr``: chance-hit rate correlated with pi;
- ``two_bins_small``: N = 300 (B.6-sized lists, where the +c term matters);
- ``ten_bins_over``: pi_hat = 1.5 x true pi (conservative misestimate);
- ``constant``: pi = 0.3 everywhere (must equal the spec formula).

Power check (no implementation needed): the rejected rules exceed the bound on ``two_bins_corr``:
mean pi over targets with E <= e (lead decision M3-10 Q2 as first written) and mean pi over all
searched peptides times (D + 1).

R is ``replicates_fast.placement`` (marker ``calibration``) or ``replicates_full.placement``
(``calibration_full``); seed ``seeds.calibration_base``; k = ``se_multiplier``.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import pandas as pd
import pytest

ALPHABET19 = "ACDEFGHKLMNPQRSTVWY"
COLUMNS = [
    "peptide", "length", "source", "search_db", "m_target", "m_decoy", "e_target", "e_decoy",
    "pep", "q_value", "pi_denovo",
]  # fmt: skip


@dataclass(frozen=True)
class Scenario:
    rates: tuple[float, ...]
    n: int = 2000
    mult: float = 1.0
    corr: bool = False


TEN = tuple(float(x) for x in np.linspace(0.02, 0.9, 10))
SCENARIOS = {
    "two_bins": Scenario((0.02, 0.9)),
    "two_bins_corr": Scenario((0.02, 0.9), corr=True),
    "ten_bins_corr": Scenario(TEN, corr=True),
    "two_bins_small": Scenario((0.02, 0.9), n=300),
    "ten_bins_over": Scenario(TEN, mult=1.5),
    "constant": Scenario((0.3,)),
}


@dataclass(frozen=True)
class Replicate:
    costs: pd.DataFrame
    wrong: np.ndarray
    pi_hat: np.ndarray
    has_t: np.ndarray
    has_d: np.ndarray
    e_t: np.ndarray
    e_d: np.ndarray


def _names(lengths: np.ndarray) -> list[str]:
    out = []
    for i, L in enumerate(lengths):
        digits = "".join(ALPHABET19[(i // 19**k) % 19] for k in (3, 2, 1, 0))
        out.append(("MDEFGHKLNPQRSTVWYAC" * 2)[: int(L) - 4] + digits)
    return out


def _replicate(rng: np.random.Generator, s: Scenario) -> Replicate:
    n = s.n
    rates = np.asarray(s.rates)
    pi = rates[rng.integers(0, rates.size, n)]
    pi_hat = np.minimum(1.0, s.mult * pi)
    p_chance = 0.05 + 0.4 * pi if s.corr else np.full(n, 0.25)
    wrong = rng.random(n) < pi
    chance_t = rng.random(n) < p_chance
    has_d = rng.random(n) < p_chance
    e_chance_t = 10 ** rng.normal(-1.5, 1.0, n)
    e_d = 10 ** rng.normal(-1.5, 1.0, n)
    e_correct = 10 ** rng.normal(-3.5, 1.0, n)
    has_t = ~wrong | chance_t
    e_t = np.where(wrong, e_chance_t, e_correct)
    lengths = rng.integers(7, 21, n)
    nan = np.full(n, np.nan)
    costs = pd.DataFrame(
        {
            "peptide": _names(lengths),
            "length": lengths,
            "source": "denovo",
            "search_db": "other",
            "m_target": np.where(has_t, 0.0, np.nan),
            "m_decoy": np.where(has_d, 0.0, np.nan),
            "e_target": np.where(has_t, e_t, np.nan),
            "e_decoy": np.where(has_d, e_d, np.nan),
            "pep": nan,
            "q_value": nan,
            "pi_denovo": pi_hat,
        },
        columns=COLUMNS,
    )
    return Replicate(costs, wrong, pi_hat, has_t, has_d, e_t, e_d)


def _ref_q(et: np.ndarray, ed: np.ndarray, fdr_of: Callable[..., np.ndarray]) -> np.ndarray:
    """q over one list: FDR at every target E (decoys with E <= e count), running min from the
    bottom, clipped to 1; ``fdr_of(order, n_t, n_d)`` gives FDR at the sorted target thresholds."""
    order = np.argsort(et, kind="stable")
    s = et[order]
    n_t = np.searchsorted(s, s, side="right").astype(float)
    n_d = np.searchsorted(np.sort(ed), s, side="right").astype(float)
    fdr = fdr_of(order, n_t, n_d)
    q = np.minimum(np.minimum.accumulate(fdr[::-1])[::-1], 1.0)
    out = np.empty_like(q)
    out[order] = q
    return out


def _q_rejected(rep: Replicate, rule: str) -> np.ndarray:
    t = np.flatnonzero(rep.has_t)
    et, ed, pit = rep.e_t[t], rep.e_d[rep.has_d], rep.pi_hat[t]
    if rule == "mean_over_targets":

        def fdr(order: np.ndarray, n_t: np.ndarray, n_d: np.ndarray) -> np.ndarray:
            mean_pi = np.cumsum(pit[order])[n_t.astype(int) - 1] / n_t
            return np.asarray(mean_pi * (n_d + 1.0) / n_t)

    else:  # mean over all searched peptides

        def fdr(order: np.ndarray, n_t: np.ndarray, n_d: np.ndarray) -> np.ndarray:
            return np.asarray(float(rep.pi_hat.mean()) * (n_d + 1.0) / n_t)

    q = np.full(rep.wrong.size, np.nan)
    q[t] = _ref_q(et, ed, fdr)
    return q


def _q_placement(rep: Replicate) -> np.ndarray:
    from giap.config import GiapConfig
    from giap.mapping.placement_fdr import placement

    per, _ = placement(rep.costs, GiapConfig())
    q = per.set_index("peptide")["placement_q"].astype(float)
    return np.asarray(q.reindex(rep.costs["peptide"]).to_numpy(), dtype=float)


def _settings(tolerances: dict[str, object], suite: str) -> tuple[list[float], int, int, float]:
    levels = [float(x) for x in tolerances["levels"]]  # type: ignore[union-attr]
    reps, seeds = tolerances[suite], tolerances["seeds"]
    assert isinstance(reps, dict) and isinstance(seeds, dict)
    k = float(tolerances["se_multiplier"])  # type: ignore[arg-type]
    return levels, int(reps["placement"]), int(seeds["calibration_base"]), k


def _fdp(
    tolerances: dict[str, object],
    suite: str,
    scenario: str,
    salt: int,
    q_of: Callable[[Replicate], np.ndarray],
) -> tuple[np.ndarray, np.ndarray, list[float], float]:
    levels, r, seed, k = _settings(tolerances, suite)
    rng = np.random.default_rng([seed, 310, salt])
    fdp = np.zeros((r, len(levels)))
    for i in range(r):
        rep = _replicate(rng, SCENARIOS[scenario])
        q = q_of(rep)
        for j, a in enumerate(levels):
            called = ~np.isnan(q) & (q <= a)
            n = int(called.sum())
            fdp[i, j] = (rep.wrong & called).sum() / n if n else 0.0
    return fdp.mean(axis=0), fdp.std(axis=0, ddof=1) / np.sqrt(r), levels, k


SUITES = [
    pytest.param("replicates_fast", marks=pytest.mark.calibration, id="fast"),
    pytest.param("replicates_full", marks=pytest.mark.calibration_full, id="full"),
]


@pytest.mark.parametrize("suite", SUITES)
@pytest.mark.parametrize("scenario", sorted(SCENARIOS))
def test_attack_signed_rule_stays_at_level(
    tolerances: dict[str, object], scenario: str, suite: str
) -> None:
    salt = sorted(SCENARIOS).index(scenario)
    mean, se, levels, k = _fdp(tolerances, suite, scenario, salt, _q_placement)
    for j, a in enumerate(levels):
        assert mean[j] <= a + k * se[j], (
            f"{scenario} q={a}: mean FDP {mean[j]:.4f} > {a} + {k}*{se[j]:.4f}"
        )


@pytest.mark.parametrize("suite", SUITES)
@pytest.mark.parametrize("rule", ["mean_over_targets", "mean_over_searched"])
def test_attack_has_power_against_rejected_rules(
    tolerances: dict[str, object], rule: str, suite: str
) -> None:
    """The attack detects the rejected rules: both exceed q + k SE at q = 0.01 and 0.05."""
    mean, se, levels, k = _fdp(
        tolerances, suite, "two_bins_corr", 99, lambda rep: _q_rejected(rep, rule)
    )
    for j, a in enumerate(levels[:2]):
        assert mean[j] > a + k * se[j], (
            f"{rule} q={a}: mean FDP {mean[j]:.4f} not above {a} + {k}*{se[j]:.4f}"
        )
