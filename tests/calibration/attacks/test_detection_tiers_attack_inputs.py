"""Attack nulls for the M4-12 review fixes (PR #92 lead decisions 5 and 11; GIPPA_SPEC.md §7.3.2,
§8.5). Added by the statistician's re-review.

Generator: the ``gene_fdr`` replicate of tests/helpers/detection_tiers_calsim.py (M4-07 -> M4-08 ->
M4-09, pi_wrong = 0.05, 2,000 genes fast / 4,500 full). FDP as in the helper.

1. ``contaminant_groups``: a random share X = 30 % of the genes share every peptide with a
   contaminant protein: their target rows get ``CON__<g>`` and their decoy rows ``rev_CON__<g>``
   appended to ``protein_ids`` (a protein group ``P;CON__P``, as an engine reports a peptide
   shared with a contaminant entry). Targets are placed by sequence (C11), so they keep counting
   in T; only the decoy side sees the contaminant member. Whatever rule decides which decoys are
   "contaminant decoys" (decision 11), it must not lift the run into ``decoy_picked`` with fewer
   decoys than false targets: the claim is gene-FDR control at q in {0.01, 0.05, 0.10} for
   whichever tier ``detect`` chooses. Under the "any stripped ID is a contaminant" rule these
   decoys leave both pairing counts, f stays ~1 and D loses X of the gene decoys, so
   FDP ~ q / (1 - X) (review run, pi 0.01, 2,000 genes, q = 0.01: 0.0118 at X = 0.2 and 0.0160
   at X = 0.4, tier always ``decoy_picked``). Pairing on the non-contaminant members (or counting
   mixed decoys as unpaired) keeps the level.
2. ``pep_only``: decoys dropped and every q-value NA (C9 ``q_value``, C11 ``best_q_input``), PEP
   kept (spec §2.1 "PEPs only -> MAYU"; decision 5: peptide q = running mean of the sorted
   accepted PEPs). All target peptides pass the q filter, so unlike the gated MAYU claims the
   q = 0.05 / 0.10 levels are reachable here. Claim: gene-FDR control at each level, tier ``mayu``.

Criterion ``fdr`` (mean FDP <= q + k SE, ADR-0015 two-block confirmation), suite ``gene_fdr``,
streams ``gene_fdr.attack_inputs.<attack>``.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detection_tiers_calsim as cs

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
SUITE = "gene_fdr"
PI = 0.05
SHARE = 0.30
ATTACKS = ("contaminant_groups", "pep_only")
_MEMO: dict[tuple[str, bool, int], tuple[dict[float, float], set[str]]] = {}


def contaminant_groups(p: pd.DataFrame, genes: np.ndarray, share: float,
                       rng: np.random.Generator) -> pd.DataFrame:  # fmt: skip
    """Append ``CON__<g>`` / ``rev_CON__<g>`` to rows of a random ``share`` of the genes."""
    p = p.copy()
    chosen = set(rng.choice(genes, size=int(share * len(genes)), replace=False).tolist())
    out = []
    decoy = p["is_decoy"].astype(bool).tolist()
    for cell, dec in zip(p["protein_ids"].tolist(), decoy, strict=True):
        ids = [str(c) for c in cell]
        base = [c[len("rev_") :] if dec and c.startswith("rev_") else c for c in ids]
        if any(b in chosen for b in base):
            ids.append(("rev_CON__" if dec else "CON__") + base[0])
        out.append(ids)
    p["protein_ids"] = pd.Series(out, index=p.index, dtype=object)
    return p


def _replicate(attack: str, full: bool, rep: int) -> tuple[dict[float, float], set[str]]:
    key = (attack, full, rep)
    if key not in _MEMO:
        n = cs.N_GENES_FULL if full else cs.N_GENES_FAST
        sim = f"gene_fdr.attack_inputs.sim.{attack}" + (".full" if full else "")
        rng = C.replicate_rng(sim, rep)
        lay, e = cs.simulate(rng, PI, n)
        if attack == "contaminant_groups":
            p, h = cs.masked(e, "picked")
            p = contaminant_groups(p, lay.genes["locus_tag"].astype(str).to_numpy(), SHARE, rng)
        else:
            p, h = cs.masked(e, "mayu_pep")
            p["q_value"] = np.nan
            h["best_q_input"] = np.nan
        _MEMO[key] = cs.run_variant(lay, e, p, h)
    return _MEMO[key]


@pytest.mark.calsuite(SUITE)
@pytest.mark.parametrize("attack", ATTACKS)
@pytest.mark.parametrize("full", MODES)
def test_gene_fdr_under_input_attack(full: bool, attack: str) -> None:
    R = C.replicates(SUITE, full)
    stream = f"gene_fdr.attack_inputs.{attack}"
    statuses = {}
    tiers: set[str] = set()
    for q in cs.LEVELS:

        def fn(rng: np.random.Generator, q: float = q) -> Any:
            fdps, t = _replicate(attack, full, cs.rep_index(stream, rng, R))
            tiers.update(t)
            return fdps[q]

        statuses[q] = C.calibrate(SUITE, "fdr", q, fn, R, stream=stream)["status"]
    if attack == "pep_only":  # tiers of the replicates run here (none in a confirm run)
        assert tiers <= {"mayu"}, tiers
        if C.owns_direct_checks():
            assert tiers == {"mayu"}, tiers
    C.skip_without_verdict([{"status": s} for s in statuses.values()])
    assert "fail" not in statuses.values(), (attack, sorted(tiers), statuses)
