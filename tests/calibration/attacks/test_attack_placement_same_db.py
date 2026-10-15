"""Attack nulls for the same-strain shortcut in the placement FDR (M3-10; GIAP_SPEC.md §7.2.5
step 5, §8.5; statistician). Methods: docs/methods/placement-fdr.md, section "Same-strain input".

Step 5 lets peptides from a search against this genome's own proteins (``search_db = same``) skip
the FDR for exact hits (``placement_q = 0``, flag ``same_db``) and leaves them out of both
stratum-0 lists. Two ways this can leak into the FDR of the *other* peptides:

- ``denovo_labelled_same``: de novo peptides are not from a protein search, but ``search_db`` is
  filled from ``input.search_db`` for every row (``io/peptides.py``), so a de novo input run with
  ``input.search_db: same`` carries ``same`` on every row. The shortcut must apply to database
  input only; otherwise every exact de novo hit, right or wrong, gets q = 0.
- ``mixed_same_other``: one input with database peptides from an own-protein search (exact hits,
  correct, small PEP) and from a cross-strain search (``other``, PEP = true wrong rate 0.3).
  pi_wrong (mean PEP) must be estimated on the peptides that enter the stratum-0 lists, the same
  population that defines c in ADR-0010; averaging in the same_db exact peptides (which have left
  both lists) halves pi_wrong and doubles the FDP of the ``other`` list.

FDP is measured over the peptides accepted through a stratum list (``same_db`` exact peptides
excluded: their q = 0 is a label, not an estimate). Criterion as everywhere: mean FDP over R
replicates <= q + k SE at every level in ``levels``. A reference implementation of the corrected
rule (``stratum_qvalues`` with pi_wrong over the searched population) is checked on the same
generator, so the attack is known to be passable. R = ``replicates_fast.placement`` (marker
``calibration``) or ``replicates_full.placement`` (``calibration_full``); k = ``se_multiplier``.
Gated by ``_criteria.calibrate`` (ADR-0015 two-block confirmation, sharding): replicate i of
stream ``placement.attack_same_db.<rule>.<scenario>`` draws from ``replicate_rng(stream, i)``.
"""

from __future__ import annotations

from collections.abc import Callable

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import Gate, calibrate_levels

ALPHABET19 = "ACDEFGHKLMNPQRSTVWY"
N = 2000
SAME_PEP = 0.001
OTHER_PI = 0.3
P_CHANCE = 0.25


def _names(lengths: np.ndarray, offset: int) -> list[str]:
    out = []
    for i, length in enumerate(lengths):
        j = i + offset
        digits = "".join(ALPHABET19[(j // 19**k) % 19] for k in (4, 3, 2, 1, 0))
        out.append(("MDEFGHKLNPQRSTVWYAC" * 2)[: int(length) - 5] + digits)
    return out


def _block(
    rng: np.random.Generator, n: int, pi: np.ndarray, offset: int
) -> tuple[dict[str, object], np.ndarray]:
    """Stratum-0 lists: chance decoy hit for all, chance target hit for wrong, low-E for correct."""
    wrong = rng.random(n) < pi
    chance_t = rng.random(n) < P_CHANCE
    has_d = rng.random(n) < P_CHANCE
    e_chance = 10 ** rng.normal(-1.5, 1.0, n)
    e_d = 10 ** rng.normal(-1.5, 1.0, n)
    e_ok = 10 ** rng.normal(-3.5, 1.0, n)
    has_t = ~wrong | chance_t
    e_t = np.where(wrong, e_chance, e_ok)
    lengths = rng.integers(8, 21, n)
    cols: dict[str, object] = {
        "peptide": _names(lengths, offset),
        "length": lengths,
        "m_target": np.where(has_t, 0.0, np.nan),
        "m_decoy": np.where(has_d, 0.0, np.nan),
        "e_target": np.where(has_t, e_t, np.nan),
        "e_decoy": np.where(has_d, e_d, np.nan),
        "q_value": np.full(n, np.nan),
    }
    return cols, wrong


def _denovo_labelled_same(rng: np.random.Generator) -> tuple[pd.DataFrame, np.ndarray]:
    pi = np.array([0.02, 0.9])[rng.integers(0, 2, N)]
    cols, wrong = _block(rng, N, pi, 0)
    df = pd.DataFrame(cols)
    df["source"] = "denovo"
    df["search_db"] = "same"  # filled from input.search_db
    df["pep"] = np.nan
    df["pi_denovo"] = pi
    return df, wrong


def _mixed_same_other(rng: np.random.Generator) -> tuple[pd.DataFrame, np.ndarray]:
    other, wrong_o = _block(rng, N, np.full(N, OTHER_PI), 0)
    d_o = pd.DataFrame(other)
    d_o["search_db"] = "other"
    d_o["pep"] = OTHER_PI
    lengths = rng.integers(8, 21, N)
    has_d = rng.random(N) < P_CHANCE
    d_s = pd.DataFrame(
        {
            "peptide": _names(lengths, N),
            "length": lengths,
            "m_target": 0.0,
            "m_decoy": np.where(has_d, 0.0, np.nan),
            "e_target": 10 ** rng.normal(-3.5, 1.0, N),
            "e_decoy": np.where(has_d, 10 ** rng.normal(-1.5, 1.0, N), np.nan),
            "q_value": np.nan,
            "search_db": "same",
            "pep": SAME_PEP,
        }
    )
    df = pd.concat([d_o, d_s], ignore_index=True)
    df["source"] = "database"
    df["pi_denovo"] = np.nan
    wrong = np.concatenate([wrong_o, np.zeros(N, dtype=bool)])
    return df, wrong


SCENARIOS: dict[str, Callable[[np.random.Generator], tuple[pd.DataFrame, np.ndarray]]] = {
    "denovo_labelled_same": _denovo_labelled_same,
    "mixed_same_other": _mixed_same_other,
}


def _q_placement(costs: pd.DataFrame) -> np.ndarray:
    from giap.config import GiapConfig
    from giap.mapping.placement_fdr import placement

    per, _ = placement(costs, GiapConfig())
    q = per.set_index("peptide")["placement_q"].astype(float)
    return np.asarray(q.reindex(costs["peptide"]).to_numpy(), dtype=float)


def _q_reference(costs: pd.DataFrame) -> np.ndarray:
    """Corrected rule, stratum 0: shortcut for database same_db exact only; pi over the rest."""
    from giap.mapping.placement_fdr import stratum_qvalues

    db = (costs["source"] == "database").to_numpy()
    m_t = costs["m_target"].to_numpy()
    m_d = costs["m_decoy"].to_numpy()
    shortcut = db & (costs["search_db"] == "same").to_numpy() & (m_t == 0)
    searched = ~shortcut
    pi_dn = costs["pi_denovo"].to_numpy(dtype=float)
    pep = costs["pep"].to_numpy(dtype=float)
    pi_db = float(np.mean(pep[searched & db])) if (searched & db).any() else 1.0
    pi = np.where(db, pi_db, np.where(np.isnan(pi_dn), 1.0, pi_dn))
    q = np.full(len(costs), np.nan)
    q[shortcut] = 0.0
    for src in (db, ~db):
        s = searched & src
        if not s.any():
            continue
        is_t, is_d = s & (m_t == 0), s & (m_d == 0)
        e_t, e_d = costs["e_target"].to_numpy(), costs["e_decoy"].to_numpy()
        q[is_t] = stratum_qvalues(e_t[is_t], e_d[is_d], pi[is_d], float(pi[s].max()))
    return q


SUITES = [
    pytest.param("replicates_fast", marks=pytest.mark.calibration, id="fast"),
    pytest.param("replicates_full", marks=pytest.mark.calibration_full, id="full"),
]


def _gate(suite: str, scenario: str, q_of: Callable[..., np.ndarray], rule: str) -> Gate:
    levels = C.load_tolerances().levels

    def fdp(rng: np.random.Generator) -> np.ndarray:
        costs, wrong = SCENARIOS[scenario](rng)
        q = q_of(costs)
        listed = ~(
            (costs["source"] == "database").to_numpy()
            & (costs["search_db"] == "same").to_numpy()
            & (costs["m_target"].to_numpy() == 0)
        )
        out = np.zeros(len(levels))
        for j, a in enumerate(levels):
            called = listed & ~np.isnan(q) & (q <= a)
            n = int(called.sum())
            out[j] = (wrong & called).sum() / n if n else 0.0
        return out

    r = C.replicates("placement", full=suite == "replicates_full")
    stream = f"placement.attack_same_db.{rule}.{scenario}"
    return calibrate_levels("placement", "fdr", levels, fdp, r, stream=stream)


@pytest.mark.parametrize("suite", SUITES)
@pytest.mark.parametrize("scenario", sorted(SCENARIOS))
def test_attack_same_db_shortcut_stays_at_level(scenario: str, suite: str) -> None:
    _gate(suite, scenario, _q_placement, "placement").assert_ok()


@pytest.mark.parametrize("suite", SUITES)
@pytest.mark.parametrize("scenario", sorted(SCENARIOS))
def test_attack_same_db_reference_rule_passes(scenario: str, suite: str) -> None:
    """The corrected rule survives the same attack, and has power (accepts something)."""
    gate = _gate(suite, scenario, _q_reference, "reference")
    gate.assert_ok()
    mean = gate.mean()  # over the replicates run here (block 1 unsharded)
    levels = C.load_tolerances().levels
    assert mean is None or mean[-1] > 0.5 * levels[-1]
