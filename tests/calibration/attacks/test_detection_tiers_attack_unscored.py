"""Attack null for M4-24: picked gene FDR with unscored and score-only decoys
(docs/specs/M4-24-decoy-unscored.md acceptance criterion "Calibration attack"; GIAP_SPEC.md §7.3.2
picked steps 1-3, §8 preamble / R4, §8.4 suite "Gene FDR (picked, MAYU)", §8.5).

Generator: the ``picked`` input of tests/helpers/detection_tiers_calsim.py (M4-07 layout, 2,000
genes fast / 4,500 full -> M4-08 -> M4-09 with pi_wrong in {0.01, 0.05}, matched decoys
``rev_<locus_tag>``, q / PEP / score emitted on both sides), then the decoy rows are masked:

- ``decoys_score_only``: decoy ``q_value`` and ``pep`` set to NA, ``score`` kept on both sides.
  S = (score,) != C_T: the tier must be ``decoy_picked`` with ``decoy_unscored`` in every
  replicate, and the picked ranking on the score alone must keep the level. Teeth (recorded in
  the M4-24 test-author hand-off): on main, before the fix, every NA-q decoy ranks below every
  target, D(t) = 0 and this variant is a confirmed ``fail`` at q = 0.01.
- ``decoys_unscored``: decoy q, PEP and score NA (giapsim M3-19 shape). S is empty: the tier must
  be ``mayu`` (rule 3) in every replicate, and the MAYU estimate must keep the level.

Claims (criterion ``fdr``, mean FDP <= q + k SE, ADR-0015 two-block confirmation, suite
``gene_fdr``; stream ``gene_fdr.attack_unscored.<variant>.pi<pi>``; R = 50 fast / 300 full), each at
q in {0.01, 0.05, 0.10} with detection rerun at ``stats.alpha.identification`` = q. The replicate
of index i draws from ``replicate_rng("gene_fdr.attack_unscored.sim.pi<pi>[.full]", i)`` and is
shared by both variants.

Teeth (private log, stays valid after the fix): decoys with q / PEP NA and a score below every
target (-1e9) are what main did to every NA-q decoy (they never win); the picked tier then gives
every target-won gene q = 1 / n_won and is a confirmed ``fail`` at q = 0.01.
"""

from __future__ import annotations

import sys
from collections.abc import Sequence
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
PIS = [0.01, 0.05]
UNSCORED = "decoy_unscored"
#: variant -> (decoy columns set to NA, required tier); ``decoy_unscored`` is required in both
VARIANTS: dict[str, tuple[tuple[str, ...], str]] = {
    "decoys_score_only": (("q_value", "pep"), "decoy_picked"),
    "decoys_unscored": (("q_value", "pep", "score"), "mayu"),
}
_MEMO: dict[tuple[float, bool, int], dict[str, dict[float, float]]] = {}


def mask_decoys(
    peptides: pd.DataFrame, cols: Sequence[str], score: float | None = None
) -> pd.DataFrame:
    """C9 copy with ``cols`` NA on the decoy rows (and their score set to ``score`` if given)."""
    p = peptides.copy()
    dec = p["is_decoy"].astype(bool).to_numpy()
    for c in cols:
        p[c] = p[c].astype("float64").mask(dec)
    if score is not None:
        p.loc[dec, "score"] = score
    return p


def run(
    lay: Any, e: Any, peptides: pd.DataFrame, hits: pd.DataFrame, levels: Sequence[float]
) -> tuple[dict[float, float], set[str], list[tuple[str, ...]]]:
    """``cs.run_variant`` that also returns the warnings of each ``detect`` call."""
    from giap.config import GiapConfig
    from giap.genetable.covariates import covariate_table
    from giap.genetable.detection import detect
    from giap.genetable.evidence import accepted_peptides, gene_evidence

    present = e.gene_truth.set_index(e.gene_truth["locus_tag"].astype(str))["present"]
    base = GiapConfig()
    acc = accepted_peptides(hits, peptides, base)
    ev = gene_evidence(lay.genes, hits, acc, lay.replicons,
                       min_peptides=base.detection.min_peptides_sequences_only)  # fmt: skip
    cov = covariate_table(lay.genes, base)
    out: dict[float, float] = {}
    tiers: set[str] = set()
    warns: list[tuple[str, ...]] = []
    for q in levels:
        c = GiapConfig(stats={"alpha": {"identification": q}})
        d = detect(lay.genes, ev, acc, peptides, hits, cs.protein_hits(e), cov, c)
        out[q] = cs.fdp(d.table, present)
        tiers.add(d.tier)
        warns.append(tuple(d.warnings))
    return out, tiers, warns


def _sim_stream(pi: float, full: bool) -> str:
    return f"gene_fdr.attack_unscored.sim.pi{pi}" + (".full" if full else "")


def _replicate(pi: float, full: bool, rep: int) -> dict[str, dict[float, float]]:
    key = (pi, full, rep)
    if key not in _MEMO:
        n = cs.N_GENES_FULL if full else cs.N_GENES_FAST
        lay, e = cs.simulate(C.replicate_rng(_sim_stream(pi, full), rep), pi, n)
        p0, h = cs.masked(e, "picked")
        assert int(p0["is_decoy"].astype(bool).sum()) > 0, "the picked input emits decoys"
        res: dict[str, dict[float, float]] = {}
        for variant, (cols, tier) in VARIANTS.items():
            fdps, tiers, warns = run(lay, e, mask_decoys(p0, cols), h, cs.LEVELS)
            assert tiers == {tier}, (variant, rep, tiers)
            assert all(UNSCORED in w for w in warns), (variant, rep, warns)
            res[variant] = fdps
        _MEMO[key] = res
    return _MEMO[key]


@pytest.mark.parametrize("variant", sorted(VARIANTS))
@pytest.mark.parametrize("pi", PIS)
@pytest.mark.parametrize("full", MODES)
def test_gene_fdr_with_unscored_decoys(full: bool, pi: float, variant: str) -> None:
    R = C.replicates(SUITE, full)
    stream = f"gene_fdr.attack_unscored.{variant}.pi{pi}"
    statuses = {}
    for q in cs.LEVELS:

        def fn(rng: np.random.Generator, q: float = q) -> float:
            return _replicate(pi, full, cs.rep_index(stream, rng, R))[variant][q]

        statuses[q] = C.calibrate(SUITE, "fdr", q, fn, R, stream=stream)["status"]
    assert "fail" not in statuses.values(), statuses


@pytest.mark.calibration
def test_teeth_decoys_that_never_win_fail(teeth_log: Path) -> None:
    pi = 0.05

    def fn(rng: np.random.Generator) -> float:
        lay, e = cs.simulate(rng, pi, cs.N_GENES_FAST)
        p, h = cs.masked(e, "picked")
        fdps, tiers, _ = run(lay, e, mask_decoys(p, ("q_value", "pep"), score=-1e9), h, (0.01,))
        assert tiers == {"decoy_picked"}, tiers
        return fdps[0.01]

    rec = C.calibrate(SUITE, "fdr", 0.01, fn, 20, stream="gene_fdr.attack_unscored.teeth",
                      log_dir=teeth_log)  # fmt: skip
    assert rec["status"] == "fail", rec


@pytest.mark.calibration
def test_masking_touches_decoy_rows_only() -> None:
    _lay, e = cs.simulate(C.replicate_rng("gene_fdr.attack_unscored.sanity", 0), 0.05, 500)
    p0, _ = cs.masked(e, "picked")
    p = mask_decoys(p0, ("q_value", "pep"))
    dec = p0["is_decoy"].astype(bool).to_numpy()
    assert dec.any() and (~dec).any()
    pd.testing.assert_frame_equal(p[~dec], p0[~dec])
    assert p.loc[dec, ["q_value", "pep"]].isna().all(axis=None)
    np.testing.assert_array_equal(p.loc[dec, "score"].to_numpy(), p0.loc[dec, "score"].to_numpy())
    assert p0.loc[dec, "score"].notna().all(), "decoys carry a score in the picked input"
