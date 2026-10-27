"""Attack null: contaminant peptides in the stratum-0 database list (GIPPA_SPEC.md §8.5;
ADR-0019 decision 1, M3-25 lead decision 7; statistician review of PR #124).

ADR-0019 estimates π̂_l over the population P, and M3-25 decision 7 takes the contaminants out of
P: they are not genome-derived, so they carry no information about the chance-hit rate. That part
is right (see the review of #124). What decision 10 left open is the other half: a contaminant
peptide stays in the stratum-0 lists, where its exact decoy hit weighs its length class's π̂ — an
estimate of the *genome-derived* peptides' no-origin share. A contaminant has no genome origin at
all, so its own π_wrong is 1, like a peptide outside every length class (ADR-0019 decision 2) or
a peptide of stratum m >= 1. Weighing it π̂_l << 1 under-counts exactly the chance hits that are
false placements by construction, and its exact target hits inflate T_0 at the same time.

Generator (one replicate = one stratum-0 database list, ``search_db = other``): N genome-derived
peptides, a share PI of them without an exact true origin, plus M_CONT contaminant peptides, never
with a true origin. A peptide without a true origin gets an exact target hit with the chance
probability p of its length; every peptide draws an exact decoy hit with the same p (A1-A2 of
docs/methods/placement-fdr.md). False = accepted without a true origin, so every accepted
contaminant counts.

Scenarios: ``flat_c200`` (p = 0.25 for every length, the convention of the sibling attacks) and
``decay_c200`` (p falls by e per 2 residues, as it does at genome scale). Both carry M_CONT = 200
contaminants against N = 4,000 peptides (5 %, an ordinary contaminant load).

- ``test_attack_contaminants_reference_passes``: the reference rule (ADR-0019 π̂ over P, a
  contaminant's own weight 1, c = max over the classes and the contaminants) holds level on both,
  so the attack is passable.
- ``test_attack_contaminants_stays_at_level``: ``placement`` as shipped. ``decay_c200`` holds
  (measured mean FDP 0.0029 at q = 0.01, R = 60): at realistic chance rates too few contaminants
  are placed to matter. ``flat_c200`` is a **strict xfail**: measured mean FDP 0.0152 (SE 0.0002,
  R = 60) against q = 0.01, about 26 SE over the level, while the reference rule gives 0.0086.
  The hole is not a regression of ADR-0019 — the superseded mean-PEP rule gives 0.0149 on the same
  replicates, because it too weighs a contaminant far below 1 — so it does not block #124; it
  waits for a decision on whether a contaminant's weight is 1 or whether contaminants leave the
  stratum lists with the ``same_db`` exact peptides.
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
N = 4000  # genome-derived database peptides
M_CONT = 200  # contaminant peptides in the same input
PI = 0.01  # share of the genome-derived peptides without an exact true origin
P_CHANCE = 0.25
CLASSES = [(7, 7), (8, 8), (9, 9), (10, 11), (12, 14), (15, 19), (20, 40)]  # mapping.length_classes

SUITES = [
    pytest.param("replicates_fast", marks=pytest.mark.calibration, id="fast"),
    pytest.param("replicates_full", marks=pytest.mark.calibration_full, id="full"),
]
XFAIL_WEIGHT = pytest.mark.xfail(
    strict=True, reason="contaminant weight pending (M3-25 decision 10 left it open)"
)


def _names(lengths: np.ndarray, tag: str) -> list[str]:
    out = []
    for i, length in enumerate(lengths):
        digits = "".join(ALPHABET19[(i // 19**k) % 19] for k in (4, 3, 2, 1, 0))
        out.append((tag + "DEFGHKLNPQRSTVWYAC" * 2)[: int(length) - 5] + digits)
    return out


def _part(
    rng: np.random.Generator, lengths: np.ndarray, p: np.ndarray, origin: np.ndarray,
    contaminant: bool, tag: str,
) -> pd.DataFrame:  # fmt: skip
    """One block of the input: ``origin`` marks the peptides with an exact true locus."""
    n = len(lengths)
    chance_t = rng.random(n) < p
    chance_d = rng.random(n) < p
    has_t = origin | chance_t
    e_ok = 10 ** rng.normal(-3.5, 1.0, n)
    e_ch = 10 ** rng.normal(-1.5, 1.0, n) * p / P_CHANCE
    return pd.DataFrame(
        {
            "peptide": _names(lengths, tag),
            "length": lengths,
            "source": "database",
            "search_db": "other",
            "m_target": np.where(has_t, 0.0, np.nan),
            "m_decoy": np.where(chance_d, 0.0, np.nan),
            "e_target": np.where(has_t, np.where(origin, e_ok, e_ch), np.nan),
            "e_decoy": np.where(chance_d, e_ch, np.nan),
            "pep": PI,
            "q_value": np.nan,
            "pi_denovo": np.nan,
            "is_contaminant": contaminant,
        }
    )


def _list(rng: np.random.Generator, decay: bool) -> tuple[pd.DataFrame, np.ndarray]:
    def rate(lengths: np.ndarray) -> np.ndarray:
        return P_CHANCE * np.exp(-(lengths - 8) / 2.0) if decay else np.full(len(lengths), P_CHANCE)

    lengths = rng.integers(8, 21, N)
    no_origin = rng.random(N) < PI
    cont_len = rng.integers(8, 21, M_CONT)
    genome = _part(rng, lengths, rate(lengths), ~no_origin, False, "M")
    cont = _part(rng, cont_len, rate(cont_len), np.zeros(M_CONT, bool), True, "K")
    return pd.concat([genome, cont], ignore_index=True), np.r_[no_origin, np.ones(M_CONT, bool)]


SCENARIOS: dict[str, Callable[[np.random.Generator], tuple[pd.DataFrame, np.ndarray]]] = {
    "flat_c200": lambda g: _list(g, False),
    "decay_c200": lambda g: _list(g, True),
}


def _pi_adr0019(costs: pd.DataFrame) -> tuple[np.ndarray, float]:
    """ADR-0019 π̂ per length class over P (the non-contaminant rows), spread over every row of
    the list, and c = max over the classes with n_l > 0; a peptide outside every class keeps 1."""
    from scipy.stats import beta

    L = costs["length"].to_numpy()
    pop = ~costs["is_contaminant"].to_numpy(dtype=bool)
    mt0 = (costs["m_target"] == 0).to_numpy()
    md0 = (costs["m_decoy"] == 0).to_numpy()
    pi = np.ones(len(costs))
    present = []
    for lo, hi in CLASSES:
        s = (lo <= L) & (hi >= L)
        sp = s & pop
        n = int(sp.sum())
        if n == 0:
            continue
        k, r = int(n - mt0[sp].sum()), float(md0[sp].mean())
        u = 1.0 if k >= n else float(beta.ppf(0.95, k + 1, n - k))
        pi[s] = 1.0 if r >= 1.0 else min(1.0, u / (1.0 - r))
        present.append(pi[s][0])
    return pi, (max(present) if present else 1.0)


def _q_stratum0(costs: pd.DataFrame, pi: np.ndarray, c: float) -> np.ndarray:
    from gippa.mapping.placement_fdr import stratum_qvalues

    m_t, m_d = costs["m_target"].to_numpy(), costs["m_decoy"].to_numpy()
    is_t, is_d = m_t == 0, m_d == 0
    q = np.full(len(costs), np.nan)
    q[is_t] = stratum_qvalues(
        costs["e_target"].to_numpy()[is_t], costs["e_decoy"].to_numpy()[is_d], pi[is_d], c
    )
    return q


def _q_reference(costs: pd.DataFrame) -> np.ndarray:
    """ADR-0019 with a contaminant's own weight 1 (it has no genome origin, so π_wrong = 1)."""
    pi, c = _pi_adr0019(costs)
    cont = costs["is_contaminant"].to_numpy(dtype=bool)
    pi = np.where(cont, 1.0, pi)
    return _q_stratum0(costs, pi, max(c, 1.0) if cont.any() else c)


def _q_placement(costs: pd.DataFrame) -> np.ndarray:
    from gippa.config import GippaConfig
    from gippa.mapping.placement_fdr import placement

    per, _ = placement(costs, GippaConfig())
    q = per.set_index("peptide")["placement_q"].astype(float)
    return np.asarray(q.reindex(costs["peptide"]).to_numpy(), dtype=float)


RULES: dict[str, Callable[[pd.DataFrame], np.ndarray]] = {
    "placement": _q_placement,
    "reference": _q_reference,
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
@pytest.mark.parametrize("scenario", [pytest.param("flat_c200", marks=XFAIL_WEIGHT), "decay_c200"])
def test_attack_contaminants_stays_at_level(scenario: str, suite: str) -> None:
    levels = list(C.load_tolerances().levels)
    r = C.replicates("placement", full=suite == "replicates_full")
    stream = f"placement.attack_contaminants.placement.{scenario}"
    calibrate_levels("placement", "fdr", levels, _fdp(scenario, "placement", levels), r,
                     stream=stream).assert_ok()  # fmt: skip


@pytest.mark.parametrize("suite", SUITES)
@pytest.mark.parametrize("scenario", sorted(SCENARIOS))
def test_attack_contaminants_reference_passes(scenario: str, suite: str) -> None:
    levels = list(C.load_tolerances().levels)
    r = C.replicates("placement", full=suite == "replicates_full")
    stream = f"placement.attack_contaminants.reference.{scenario}"
    gate = calibrate_levels("placement", "fdr", levels, _fdp(scenario, "reference", levels), r,
                            stream=stream)  # fmt: skip
    gate.assert_ok()
    mean = gate.mean()
    assert mean is None or mean[0] > 0.1 * levels[0]  # not vacuous: false placements do occur


@pytest.mark.parametrize("suite", SUITES)
def test_attack_has_power_against_the_class_weight(suite: str, teeth_log: Path) -> None:
    """Teeth: on ``flat_c200`` the shipped weight is a confirmed ``fail`` at q = 0.01."""
    levels = list(C.load_tolerances().levels)
    r = C.replicates("placement", full=suite == "replicates_full")
    stream = "placement.attack_contaminants.rejected.class_weight"
    gate = calibrate_levels("placement", "fdr", levels, _fdp("flat_c200", "placement", levels), r,
                            stream=stream, log_dir=teeth_log)  # fmt: skip
    assert gate.records[0]["status"] == "fail", gate.detail()
