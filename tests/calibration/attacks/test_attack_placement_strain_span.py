"""Attack nulls for stratum-0 pi_wrong with cross-strain database input (ADR-0019; GIAP_SPEC.md
§7.2.5 step 4, §8.5; statistician). Methods: docs/methods/placement-fdr.md, "Cross-strain input".

Status: committed with M3-21 (PR #113); the ``placement`` test's strict xfail was removed with the
ADR-0019 implementation (branch m3/adr0019-pi-hat), so it gates ``placement`` itself. With
``search_db = other`` some *correct* identifications are database-strain sequences spanning a strain
difference: they have no exact true locus, so an exact chance hit elsewhere is a false stratum-0
placement, while their PEP (identification error only) is small. Mean PEP then counts only the
identification errors among the peptides without an exact origin, and the FDR is under-estimated by
about (pi_wrong + delta) / pi_wrong.

Generator (one replicate = one stratum-0 database list, ``search_db = other``, every PEP = pi):
origin peptides (exact true locus, low-E target hit), wrong peptides (share pi) and spanning
peptides (share delta of the correct ones): wrong and spanning peptides get a chance target hit
w.p. p_j; every peptide a chance decoy hit w.p. p_j (A1-A2 of the methods note).
Scenarios: ``flat_d05`` / ``flat_d20`` (p_j = 0.25, delta = 0.05 / 0.20); ``short_d20`` (wrong
and spanning peptides 8-10 aa, p_j falls with length: the length-correlated case); ``pad_d20``
(flat_d20 plus N long exact peptides with PEP 0: padding attack on any single pooled pi_hat).
Power at delta = 0 relative to mean PEP is not tested here; it was measured for ADR-0019
(Evidence: 0.984-1.0 at N = 4,000, Slurm 1860138 / 1860180).

- ``test_attack_strain_span_stays_at_level``: through ``placement`` (fails under mean PEP; a
  normal test since the ADR-0019 implementation).
- ``test_attack_strain_span_reference_passes``: reference ADR-0019 rule (per length class pi_hat
  from the exact-match shares, 95 % upper bound, ADR-0010 decoy weights, c = max) on the same
  generator.
- ``test_attack_has_power_against_mean_pep``: teeth: the mean-PEP rule is a confirmed ``fail`` at
  q = 0.01 and 0.05 on ``flat_d20`` (PR #69 decision 3).
- ``test_estimator_edge_cases``: the reference rule's edge cases (r = 1, k = n, r = 0, empty
  classes, c = max).
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels

ALPHABET19 = "ACDEFGHKLMNPQRSTVWY"
N = 4000
PI = 0.01
P_CHANCE = 0.25

SUITES = [
    pytest.param("replicates_fast", marks=pytest.mark.calibration, id="fast"),
    pytest.param("replicates_full", marks=pytest.mark.calibration_full, id="full"),
]


def _names(lengths: np.ndarray) -> list[str]:
    out = []
    for i, length in enumerate(lengths):
        digits = "".join(ALPHABET19[(i // 19**k) % 19] for k in (4, 3, 2, 1, 0))
        out.append(("MDEFGHKLNPQRSTVWYAC" * 2)[: int(length) - 5] + digits)
    return out


def _list(rng: np.random.Generator, delta: float, short: bool) -> tuple[pd.DataFrame, np.ndarray]:
    u = rng.random(N)
    d1 = delta * (1 - PI)
    no_origin = u < PI + d1  # wrong (u < PI) or spanning
    lengths = rng.integers(8, 21, N)
    if short:
        lengths = np.where(no_origin, rng.integers(8, 11, N), lengths)
        p = np.clip(P_CHANCE * np.exp(-(lengths - 8) / 3.0), 0.0, 1.0)
    else:
        p = np.full(N, P_CHANCE)
    chance_t = rng.random(N) < p
    chance_d = rng.random(N) < p
    e_ct = 10 ** rng.normal(-1.5, 1.0, N) * p / P_CHANCE
    e_cd = 10 ** rng.normal(-1.5, 1.0, N) * p / P_CHANCE
    e_ok = 10 ** rng.normal(-3.5, 1.0, N)
    has_t = ~no_origin | chance_t
    df = pd.DataFrame(
        {
            "peptide": _names(lengths),
            "length": lengths,
            "source": "database",
            "search_db": "other",
            "m_target": np.where(has_t, 0.0, np.nan),
            "m_decoy": np.where(chance_d, 0.0, np.nan),
            "e_target": np.where(has_t, np.where(no_origin, e_ct, e_ok), np.nan),
            "e_decoy": np.where(chance_d, e_cd, np.nan),
            "pep": PI,
            "q_value": np.nan,
            "pi_denovo": np.nan,
        }
    )
    return df, no_origin


def _padded(rng: np.random.Generator) -> tuple[pd.DataFrame, np.ndarray]:
    """flat_d20 plus N long exact peptides with PEP 0 and no chance hits (padding attack)."""
    df, no_origin = _list(rng, 0.20, False)
    pad = pd.DataFrame(
        {
            "peptide": [x + "W" for x in _names(np.full(N, 29))], "length": 30,
            "source": "database", "search_db": "other", "m_target": 0.0, "m_decoy": np.nan,
            "e_target": 10 ** rng.normal(-8.0, 1.0, N), "e_decoy": np.nan, "pep": 0.0,
            "q_value": np.nan, "pi_denovo": np.nan,
        }
    )  # fmt: skip
    return pd.concat([df, pad], ignore_index=True), np.r_[no_origin, np.zeros(N, bool)]


SCENARIOS: dict[str, Callable[[np.random.Generator], tuple[pd.DataFrame, np.ndarray]]] = {
    "flat_d05": lambda g: _list(g, 0.05, False),
    "flat_d20": lambda g: _list(g, 0.20, False),
    "short_d20": lambda g: _list(g, 0.20, True),
    "pad_d20": _padded,
}


CLASSES = [(7, 7), (8, 8), (9, 9), (10, 11), (12, 14), (15, 19), (20, 40)]  # mapping.length_classes


def _pi_adr0019(costs: pd.DataFrame) -> tuple[np.ndarray, float]:
    """ADR-0019 per-peptide pi_hat (its length class) and c = max over classes present.

    Class l: n peptides of the stratum-0 database population, k without an exact target hit, r the
    share with an exact decoy hit; pi_l = min(1, U95(k, n) / (1 - r)), U95 the one-sided 95 %
    Clopper-Pearson upper bound (1 when k = n); pi_l = 1 when r = 1. No PEP is used.
    """
    from scipy.stats import beta

    L = costs["length"].to_numpy()
    mt0 = (costs["m_target"] == 0).to_numpy()
    md0 = (costs["m_decoy"] == 0).to_numpy()
    pi = np.ones(len(costs))
    present = []
    for lo, hi in CLASSES:
        s = (lo <= L) & (hi >= L)
        n = int(s.sum())
        if n == 0:
            continue
        k, r = int(n - mt0[s].sum()), float(md0[s].mean())
        u = 1.0 if k >= n else float(beta.ppf(0.95, k + 1, n - k))
        pi[s] = 1.0 if r >= 1.0 else min(1.0, u / (1.0 - r))
        present.append(pi[s][0])
    return pi, (max(present) if present else 1.0)


def _q_stratum0(costs: pd.DataFrame, pi: np.ndarray, c: float) -> np.ndarray:
    from giap.mapping.placement_fdr import stratum_qvalues

    m_t, m_d = costs["m_target"].to_numpy(), costs["m_decoy"].to_numpy()
    is_t, is_d = m_t == 0, m_d == 0
    q = np.full(len(costs), np.nan)
    q[is_t] = stratum_qvalues(
        costs["e_target"].to_numpy()[is_t], costs["e_decoy"].to_numpy()[is_d], pi[is_d], c
    )
    return q


def _q_placement(costs: pd.DataFrame) -> np.ndarray:
    from giap.config import GiapConfig
    from giap.mapping.placement_fdr import placement

    per, _ = placement(costs, GiapConfig())
    q = per.set_index("peptide")["placement_q"].astype(float)
    return np.asarray(q.reindex(costs["peptide"]).to_numpy(), dtype=float)


RULES: dict[str, Callable[[pd.DataFrame], np.ndarray]] = {
    "placement": _q_placement,
    "reference": lambda c: _q_stratum0(c, *_pi_adr0019(c)),
    "mean_pep": lambda c: _q_stratum0(
        c, np.full(len(c), float(c["pep"].mean())), float(c["pep"].mean())
    ),
}


def _fdp(
    scenario: str, rule: str, levels: list[float]
) -> Callable[[np.random.Generator], np.ndarray]:
    def fn(rng: np.random.Generator) -> np.ndarray:
        costs, no_origin = SCENARIOS[scenario](rng)
        q = RULES[rule](costs)
        out = np.zeros(len(levels))
        for j, a in enumerate(levels):
            called = ~np.isnan(q) & (q <= a)
            n = int(called.sum())
            out[j] = (no_origin & called).sum() / n if n else 0.0
        return out

    return fn


@pytest.mark.parametrize("suite", SUITES)
@pytest.mark.parametrize("scenario", sorted(SCENARIOS))
def test_attack_strain_span_stays_at_level(scenario: str, suite: str) -> None:
    levels = list(C.load_tolerances().levels)
    r = C.replicates("placement", full=suite == "replicates_full")
    stream = f"placement.attack_strain_span.placement.{scenario}"
    calibrate_levels("placement", "fdr", levels, _fdp(scenario, "placement", levels), r,
                     stream=stream).assert_ok()  # fmt: skip


@pytest.mark.parametrize("suite", SUITES)
@pytest.mark.parametrize("scenario", sorted(SCENARIOS))
def test_attack_strain_span_reference_passes(scenario: str, suite: str) -> None:
    levels = list(C.load_tolerances().levels)
    r = C.replicates("placement", full=suite == "replicates_full")
    stream = f"placement.attack_strain_span.reference.{scenario}"
    gate = calibrate_levels("placement", "fdr", levels, _fdp(scenario, "reference", levels), r,
                            stream=stream)  # fmt: skip
    gate.assert_ok()
    mean = gate.mean()
    assert mean is None or mean[0] > 0.5 * levels[0]  # not vacuous: the list binds at q = 0.01


@pytest.mark.parametrize("suite", SUITES)
def test_attack_has_power_against_mean_pep(suite: str, teeth_log: Path) -> None:
    levels = list(C.load_tolerances().levels)
    r = C.replicates("placement", full=suite == "replicates_full")
    stream = "placement.attack_strain_span.rejected.mean_pep"
    gate = calibrate_levels("placement", "fdr", levels, _fdp("flat_d20", "mean_pep", levels), r,
                            stream=stream, log_dir=teeth_log)  # fmt: skip
    assert [rec["status"] for rec in gate.records[:2]] == ["fail", "fail"], gate.detail()


@pytest.mark.calibration
def test_estimator_edge_cases() -> None:
    """r = 1 -> 1; k = n -> 1; r = 0 -> U95(k, n); empty classes ignored; c = max class pi."""
    from scipy.stats import beta

    def costs(mt: list[float], md: list[float], length: int = 9) -> pd.DataFrame:
        return pd.DataFrame({"m_target": mt, "m_decoy": md, "length": length})

    nan = np.nan
    assert _pi_adr0019(costs([0, 0], [0, 0]))[1] == 1.0  # r = 1
    assert _pi_adr0019(costs([nan, nan], [nan, nan]))[1] == 1.0  # k = n
    pi, c = _pi_adr0019(costs([0] * 99 + [nan], [nan] * 100))  # r = 0, k = 1 of 100
    assert c == pytest.approx(float(beta.ppf(0.95, 2, 99)))
    two = pd.concat([costs([0] * 100, [nan] * 100, 9), costs([nan] * 10, [nan] * 10, 25)])
    pi, c = _pi_adr0019(two.reset_index(drop=True))
    assert c == 1.0 and pi[0] < 0.05  # per-class labels; c = max
