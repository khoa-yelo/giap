"""Gene-level FDR calibration of M4-12 detection tiers (docs/specs/M4-12-detection-tiers.md,
acceptance criteria 3-5; GIAP_SPEC.md §7.3.2, §8.4 suite "Gene FDR (picked, MAYU) | false peptides
on random genes, size-proportional | FDR | 50 | 300").

Suite ``gene_fdr`` (SUITE_REGISTRY file entry; the spec's ``test_gene_fdr.py`` carries the
``detection_tiers_`` prefix for unique test basenames). R = ``replicates_fast.gene_fdr`` (50) /
``replicates_full.gene_fdr`` (300); replicate seeds from the harness: the simulation of replicate i
draws from ``replicate_rng("gene_fdr.detect.sim.pi<pi>[.full]", i)`` and is shared by the claims
below (each claim has its own harness stream, which only names the check). Generator and FDP:
tests/helpers/detection_tiers_calsim.py (M4-07 layout, 2,000 genes fast / 4,500 full -> M4-08 ->
M4-09 with pi_wrong in {0.01, 0.05} -> M4-10/M4-11 -> ``detect``).

Statistical claims (criterion ``fdr``: mean FDP <= q + k SE, ADR-0015 two-block confirmation),
each at q in {0.01, 0.05, 0.10} with detection rerun at ``stats.alpha.identification`` = q:

1. ``gene_fdr.picked.pi<pi>``: picked gene-level target-decoy FDR (decoys emitted; tier
   ``decoy_picked`` asserted per replicate) controls gene FDR.
2. ``gene_fdr.mayu_pep.pi<pi>``: decoys dropped, PEPs present: the MAYU tier (F = sum PEP)
   controls gene FDR.
3. ``gene_fdr.mayu_nopep.pi<pi>``: decoys and PEPs dropped: the MAYU tier with F = t |P(t)|.
4. Recorded, not gated (spec: the tier is ``unestimated`` by design): mean FDP of the
   ``two_peptide`` tier on the same replicates, appended to ``<log_dir>/gene_fdr_two_peptide.jsonl``
   (``_criteria.log_result``; one record per pi, ``observed`` = mean FDP).
5. Teeth (private log): detection with no q filter at all (``input.q_value_max`` = 1,
   ``stats.alpha.identification`` = 1) is a confirmed ``fail`` of the FDR criterion at q = 0.01.

ADR-0009 (decoy loss 10/20/30 %) is the attack null
tests/calibration/attacks/test_detection_tiers_attack_decoy_loss.py.
"""

from __future__ import annotations

import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import detection_tiers_calsim as cs

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
PIS = [0.01, 0.05]
SUITE = "gene_fdr"
GATED = ("picked", "mayu_pep", "mayu_nopep")
_MEMO: dict[tuple[float, bool, int], dict[str, Any]] = {}


def _sim_stream(pi: float, full: bool) -> str:
    return f"gene_fdr.detect.sim.pi{pi}" + (".full" if full else "")


def _replicate(pi: float, full: bool, rep: int) -> dict[str, Any]:
    """FDP per level for every variant of replicate ``rep`` (memoised; one simulation)."""
    key = (pi, full, rep)
    if key not in _MEMO:
        n = cs.N_GENES_FULL if full else cs.N_GENES_FAST
        lay, e = cs.simulate(C.replicate_rng(_sim_stream(pi, full), rep), pi, n)
        out: dict[str, Any] = {}
        for variant in (*GATED, "two_peptide"):
            p, h = cs.masked(e, variant)
            levels = cs.LEVELS if variant != "two_peptide" else (cs.LEVELS[0],)
            fdps, tiers = cs.run_variant(lay, e, p, h, levels)
            assert tiers == {cs.TIER_OF[variant]}, (variant, tiers)
            out[variant] = fdps
        _MEMO[key] = out
    return _MEMO[key]


def _claim(pi: float, full: bool, variant: str, q: float, stream: str,
           R: int) -> Callable[[np.random.Generator], float]:  # fmt: skip
    def fn(rng: np.random.Generator) -> float:
        rep = cs.rep_index(stream, rng, R)
        return float(_replicate(pi, full, rep)[variant][q])

    return fn


@pytest.mark.parametrize("pi", PIS)
@pytest.mark.parametrize("full", MODES)
def test_gene_fdr_controlled(full: bool, pi: float) -> None:
    R = C.replicates(SUITE, full)
    statuses = {}
    for variant in GATED:
        stream = f"gene_fdr.{variant}.pi{pi}"
        for q in cs.LEVELS:
            fn = _claim(pi, full, variant, q, stream, R)
            rec = C.calibrate(SUITE, "fdr", q, fn, R, stream=stream)
            statuses[(variant, q)] = rec["status"]
    assert "fail" not in statuses.values(), statuses
    if C.replicate_range(R) == range(R) and C.confirm_checks() is None:
        values = [_replicate(pi, full, i)["two_peptide"][cs.LEVELS[0]] for i in range(R)]
        res = C.check_fdr(values, cs.LEVELS[0])
        seed = C.suite_seed(_sim_stream(pi, full), 0)
        path = C.log_result("gene_fdr_two_peptide", "fdr", cs.LEVELS[0], res, seed, 0.0)
        assert path.exists()
        assert 0.0 <= res.observed <= 1.0


@pytest.mark.calibration
def test_teeth_no_q_filter_fails(teeth_log: Path) -> None:
    pi = 0.05

    def fn(rng: np.random.Generator) -> float:
        lay, e = cs.simulate(rng, pi, cs.N_GENES_FAST)
        p, h = cs.masked(e, "picked")
        fdps, _ = cs.run_variant(lay, e, p, h, (1.0,), input={"q_value_max": 1.0})
        return fdps[1.0]

    rec = C.calibrate(SUITE, "fdr", 0.01, fn, 20, stream="gene_fdr.detect.teeth_no_q_filter",
                      log_dir=teeth_log)  # fmt: skip
    assert rec["status"] == "fail", rec


@pytest.mark.calibration
def test_replicate_sanity() -> None:
    """The generator gives the claims something to bite on: false genes exist among genes with
    accepted evidence, decoys pass q <= 0.01 and all pair (rev_<locus_tag>)."""
    lay, e = cs.simulate(C.replicate_rng("gene_fdr.detect.sanity", 0), 0.05, cs.N_GENES_FAST)
    d = e.peptides[e.peptides["is_decoy"].astype(bool)]
    assert int((d["q_value"] <= 0.01).sum()) >= 5
    tags = set(lay.genes["locus_tag"].astype(str))
    assert all(ids[0][len("rev_") :] in tags for ids in d["protein_ids"])
    absent = ~e.gene_truth["present"].astype(bool)
    assert absent.sum() >= 200
