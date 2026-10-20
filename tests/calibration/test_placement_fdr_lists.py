"""Statistics-only calibration of placement q-values (issue M3-10; GIAP_SPEC.md §7.2.5, §8.4 "FDR").

Criterion: mean FDP over R replicates <= q + k SE(FDP) for q in ``levels``; R =
``replicates_fast.placement``, k = ``se_multiplier`` (tests/calibration/tolerances.yaml). Gated by
``_criteria.calibrate`` (ADR-0015 two-block confirmation, sharding): replicate i of stream
``placement.lists.<scenario>`` draws from ``replicate_rng(stream, i)``.

Generator (one replicate, one stratum m): N peptides, each wrong with a known probability.
- Every peptide (correct or wrong) is searched against the decoy (lead decision M3-10 Q1) and gets a
  chance decoy hit with probability ``P_CHANCE``, E drawn from the chance law.
- A wrong peptide gets a chance target hit independently with the same probability and the same E
  law (target and decoy exchangeable for wrong peptides); otherwise it has no target hit.
- A correct peptide always has its target hit, with E from a low-E law.
So D estimates the chance hits of the whole input and pi_wrong * D the false targets: the factor
pi_wrong is the share of wrong peptides among all searched peptides. (The issue text gives correct
peptides "no decoy hit"; with pi_wrong < 1 that would make the spec formula anti-conservative by the
factor pi_wrong, so correct peptides draw chance decoy hits like every other searched peptide.)

For database input, stratum 0, ADR-0019 replaces §7.2.5 step 4: π̂ is the per-length-class
exact-share estimate min(1, U95(k, n)/(1 - r)) and PEP and q-values are not read. In this
generator a correct peptide always has its exact target hit and a wrong one has one with
probability P_CHANCE, so k/n estimates pi_wrong (1 - P_CHANCE) and r = P_CHANCE: the estimator
is unbiased up to the upper bound. The scenarios that used to exercise each step-4 route stay as
regression scenarios for the class estimator under different wrong-share laws:
- ``known_pi``: every PEP equals the true wrong probability (database, stratum 0 and 1);
- ``calibrated_pep``: per-peptide PEP ~ Beta, wrong ~ Bernoulli(PEP);
- ``max_q``: no PEP, q-values up to ``input.q_value_max`` (formerly the max-q route);
- ``deflated_pep`` (ADR-0019 "Run collapse"): PEP = a tenth of the true wrong probability, as the
  §5.3 min over runs produces; mean PEP fails this scenario (FDP ~ 10 q), the class estimator must
  not (it does not read PEP);
- ``denovo_bins``: de novo peptides in 10 score bins with different wrong rates, ``pi_denovo`` = the
  true per-bin wrong rate; per-peptide pi enters as decoy weights, FDR(e) = (W(e) + c)/T(e) with
  c = max pi (ADR-0010; docs/methods/placement-fdr.md "pi_wrong in the stratum FDR", where this
  scenario gave mean FDP 0.0095 / 0.050 / 0.100 at R = 300).
With constant pi (the first three scenarios) ADR-0010 reduces to the §7.2.5 formula.
"""

from __future__ import annotations

from collections.abc import Callable

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels

from giap.config import GiapConfig
from giap.mapping.placement_fdr import placement

pytestmark = pytest.mark.calibration

SUITE = "placement"

N = 4000
P_CHANCE = 0.25
ALPHABET19 = "ACDEFGHKLMNPQRSTVWY"
COLUMNS = [
    "peptide", "length", "source", "search_db", "m_target", "m_decoy", "e_target", "e_decoy",
    "pep", "q_value", "pi_denovo",
]  # fmt: skip


def _names(lengths: np.ndarray) -> list[str]:
    out = []
    for i, L in enumerate(lengths):
        digits = "".join(ALPHABET19[(i // 19**k) % 19] for k in (3, 2, 1, 0))
        out.append(("MDEFGHKLNPQRSTVWYAC" * 2)[: int(L) - 4] + digits)
    return out


def _replicate(
    rng: np.random.Generator,
    m: int,
    wrong_prob: np.ndarray,
    *,
    source: str = "database",
    pep: np.ndarray | None = None,
    q_value: np.ndarray | None = None,
    pi_denovo: np.ndarray | None = None,
) -> tuple[pd.DataFrame, np.ndarray]:
    n = wrong_prob.size
    wrong = rng.random(n) < wrong_prob
    chance_t = rng.random(n) < P_CHANCE
    chance_d = rng.random(n) < P_CHANCE
    e_chance_t = 10 ** rng.normal(-1.5, 1.0, n)
    e_chance_d = 10 ** rng.normal(-1.5, 1.0, n)
    e_correct = 10 ** rng.normal(-3.5, 1.0, n)
    has_t = ~wrong | chance_t
    e_t = np.where(wrong, e_chance_t, e_correct)
    lengths = rng.integers(7, 21, n)
    nan = np.full(n, np.nan)
    df = pd.DataFrame(
        {
            "peptide": _names(lengths),
            "length": lengths,
            "source": source,
            "search_db": "other",
            "m_target": np.where(has_t, float(m), np.nan),
            "m_decoy": np.where(chance_d, float(m), np.nan),
            "e_target": np.where(has_t, e_t, np.nan),
            "e_decoy": np.where(chance_d, e_chance_d, np.nan),
            "pep": nan if pep is None else pep,
            "q_value": nan if q_value is None else q_value,
            "pi_denovo": nan if pi_denovo is None else pi_denovo,
        },
        columns=COLUMNS,
    )
    return df, wrong


def _check_fdr(
    name: str, make: Callable[[np.random.Generator], tuple[pd.DataFrame, np.ndarray]]
) -> None:
    levels = C.load_tolerances().levels
    cfg = GiapConfig()
    n_acc: list[float] = []  # accepted at the loosest level, per replicate run here

    def fdp(rng: np.random.Generator) -> np.ndarray:
        costs, wrong = make(rng)
        per, _ = placement(costs, cfg)
        q = per.set_index("peptide")["placement_q"].astype(float)
        q = q.reindex(costs["peptide"]).to_numpy()
        out = np.zeros(len(levels))
        for j, a in enumerate(levels):
            called = ~np.isnan(q) & (q <= a)
            n = int(called.sum())
            out[j] = (wrong & called).sum() / n if n else 0.0
            if j == len(levels) - 1:
                n_acc.append(n)
        return out

    r = C.replicates(SUITE, full=False)
    calibrate_levels(SUITE, "fdr", levels, fdp, r, stream=f"{SUITE}.lists.{name}").assert_ok()
    # non-trivial: something is accepted at the loosest level
    assert not n_acc or np.mean(n_acc) > 0


@pytest.mark.parametrize("m", [0, 1])
def test_known_pi_controls_fdr(m: int) -> None:
    pi = 0.3

    def make(rng: np.random.Generator) -> tuple[pd.DataFrame, np.ndarray]:
        return _replicate(rng, m, np.full(N, pi), pep=np.full(N, pi))

    _check_fdr(f"known_pi.m{m}", make)


def test_calibrated_pep_controls_fdr() -> None:
    def make(rng: np.random.Generator) -> tuple[pd.DataFrame, np.ndarray]:
        pep = rng.beta(0.6, 1.4, N)
        return _replicate(rng, 0, pep, pep=pep)

    _check_fdr("calibrated_pep", make)


def test_max_q_route_controls_fdr() -> None:
    def make(rng: np.random.Generator) -> tuple[pd.DataFrame, np.ndarray]:
        q = rng.uniform(0.0, 0.01, N)  # max q ~ 0.01 >= the true wrong share 0.008
        return _replicate(rng, 0, np.full(N, 0.008), q_value=q)

    _check_fdr("max_q", make)


def test_deflated_pep_controls_fdr() -> None:
    """ADR-0019 Evidence row "PEP deflated 10x (min-collapse)": 14 failing cells under mean PEP,
    0 under the class estimator."""
    pi = 0.3

    def make(rng: np.random.Generator) -> tuple[pd.DataFrame, np.ndarray]:
        return _replicate(rng, 0, np.full(N, pi), pep=np.full(N, pi / 10))

    _check_fdr("deflated_pep", make)


def test_denovo_per_bin_pi_controls_fdr() -> None:
    rates = np.linspace(0.02, 0.9, 10)

    def make(rng: np.random.Generator) -> tuple[pd.DataFrame, np.ndarray]:
        w = rates[rng.integers(0, rates.size, N)]
        return _replicate(rng, 0, w, source="denovo", pi_denovo=w)

    _check_fdr("denovo_bins", make)
