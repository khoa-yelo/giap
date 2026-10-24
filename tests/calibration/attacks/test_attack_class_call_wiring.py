"""Attack nulls on the M7-11b wired call path (GIAP_SPEC.md §7.6.4, §8.1, §8.5; ADR-0009;
docs/specs/M7-11b-class-fdr-wiring.md rules 2, 5, 6; statistician, PR #256 review).

``tests/calibration/test_class_call_calibration.py`` measures ORFs whose class-FDR *cell* is
numeric and <= a. Two parts of the wiring are outside that claim, and these attacks target them.
Criterion ``fdr`` (mean FDP <= a + 2 SE, ADR-0015 two-block confirmation), suite ``novel_orfs``,
a in ``levels``, streams ``novel_orfs.attack.class_call.<name>``.

- ``lost_pick_called``: FDP of the **``called``** column from ``apply_evidence`` (rule 6), run
  with ``stats.alpha.novel_orfs = a``. The world is ``cds_outframe`` ORFs with 2-3 peptides each,
  so every ORF meets its class rule (>= 2 peptides) and would pass the ``unestimated`` fallback.
  25 % of entries are true. A false entry's target and its decoy are exchangeable (best of 1-3
  N(0, 1)). About half of the false targets lose the pick. Rule 5 (``lost_pick`` -> 1, Q1) is
  what keeps them uncalled. Teeth: with ``class_fdr_cell`` mapping ``lost_pick`` to
  ``"unestimated"`` (the rejected reading of Q1), the fallback calls every pick loser and FDP is
  far above level (checked by ``test_lost_pick_attack_has_teeth``).
- ``repointed_decoys``: the ADR-0009 decoy-loss attack, routed through the pipeline's own
  ``makedb_id_map`` (rule 2, no ``db_id_map.tsv``). 20 % or 50 % of decoy rows are re-pointed to
  class-tagged IDs whose targets do not exist. With a real map these decoys would be unpaired and
  the 80 % guard would fire at 50 %. In v1 every observed tagged ID becomes an entry, so the paired
  share is 1 (Q5). This attack checks that a share of 1 is not anti-conservative: the re-pointed
  decoys enter as decoy-only entries, which are decoy wins, so D is not lost. FDP of the numeric
  cells is measured within the class.
"""

from __future__ import annotations

import math
from collections.abc import Callable

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels
from helpers.class_call_fixture import config, gene_model, gene_table, world_genome
from helpers.novel_fdr_fixture import PickedWorld

pytestmark = pytest.mark.calsuite("novel_orfs")

SUITE = "novel_orfs"
MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
GENOME = world_genome()
GM = gene_model(GENOME)
CONTIGS = tuple(GENOME)
GT = gene_table()


def _cc():  # type: ignore[no-untyped-def]
    import giap.novel.class_call as m

    return m


def _q_in(score: float) -> float:
    return 0.01 / (1.0 + math.exp(score))


def _best(rng: np.random.Generator, mu: float) -> float:
    return float(rng.normal(mu, 1.0, rng.integers(1, 4)).max())


# ---------------------------------------------------------------- lost_pick_called

N_ENTRIES = 400  # 100 per contig: entries at 1000 + 600 i stay inside the 200-kb contigs
P_TRUE = 0.25
P_HIT = 0.5


def _cds_world(rng: np.random.Generator) -> tuple[PickedWorld, set[str]]:
    w = PickedWorld()
    false_ids: set[str] = set()
    for k in range(N_ENTRIES):
        contig = CONTIGS[k % len(CONTIGS)]
        true = rng.random() < P_TRUE
        has_t = true or rng.random() < P_HIT
        has_d = rng.random() < P_HIT
        if not (has_t or has_d):
            continue
        entry = w.entry("cds_outframe", contig, "+")
        if has_t:
            s = _best(rng, 3.0 if true else 0.0)
            n_pep = int(rng.integers(2, 4))
            orf, _ = w.target("cds_outframe", contig, [0] * n_pep, entry=entry, q=_q_in(s),
                              score=s)  # fmt: skip
            if not true:
                false_ids.add(orf.orf_id)
        if has_d:
            s = _best(rng, 0.0)
            w.decoy(entry[0], 0, q=_q_in(s), score=s)
    return w, false_ids


def _called_fdp(
    cell: Callable[[object], float | str],
) -> Callable[[np.random.Generator], np.ndarray]:
    """FDP of ``called`` at every level: class FDR once, then ``apply_evidence`` per alpha, with the
    step's float ``class_q`` (lead decision 6, gap 2)."""
    from giap.novel.candidates import minimal_call
    from giap.novel.evidence import apply_evidence
    from giap.novel.outputs import orf_table

    levels = C.load_tolerances().levels

    def fdp(rng: np.random.Generator) -> np.ndarray:
        w, false_ids = _cds_world(rng)
        res = _cc().class_fdr(w.orfs, w.peptides, gm=GM, genome=GENOME, cfg=config(),
                              denovo=None)  # fmt: skip
        scored = minimal_call(w.orfs, config())
        scored["class_fdr"] = pd.Series([cell(res.orfs[i]) for i in scored["orf_id"]],
                                        dtype="object")  # fmt: skip
        pre = orf_table(w.orfs, scored, genome=GENOME, replicons=GM.replicons)
        class_q = {k: r.q for k, r in res.orfs.items() if math.isfinite(r.q)}
        out = np.zeros(len(levels))
        for j, a in enumerate(levels):
            t = apply_evidence(pre, w.orfs, gene_table=GT, terminators=None, gm=GM,
                               genome=GENOME, cfg=config(a), peptide_intervals={},
                               class_q=class_q)  # fmt: skip
            called = set(t.loc[t["called"].astype(bool), "orf_id"].astype(str))
            out[j] = len(called & false_ids) / len(called) if called else 0.0
        return out

    return fdp


def _lost_pick_as_unestimated(r: object) -> float | str:
    """The rejected Q1 reading: a pick loser gets ``unestimated`` (and so the fallback)."""
    v = _cc().class_fdr_cell(r)
    return "unestimated" if getattr(r, "reason", None) == "lost_pick" else v


@pytest.mark.parametrize("full", MODES)
def test_lost_pick_called_fdp_at_level(full: bool) -> None:
    levels = C.load_tolerances().levels
    fn = _called_fdp(lambda r: _cc().class_fdr_cell(r))
    calibrate_levels(SUITE, "fdr", levels, fn, C.replicates(SUITE, full=full),
                     stream=f"{SUITE}.attack.class_call.lost_pick_called").assert_ok()  # fmt: skip


@pytest.mark.calibration
def test_lost_pick_attack_has_teeth() -> None:
    """Under the rejected reading the fallback calls every pick loser: FDP fails at every level.
    Five replicates suffice: the effect is not marginal (every false pick loser is called)."""
    levels = C.load_tolerances().levels
    fn = _called_fdp(_lost_pick_as_unestimated)
    rngs = [np.random.default_rng(C.suite_seed(f"{SUITE}.attack.class_call.teeth", i))
            for i in range(5)]  # fmt: skip
    mean = np.mean([fn(g) for g in rngs], axis=0)
    assert all(m > 2 * a for m, a in zip(mean, levels, strict=True)), mean


# ---------------------------------------------------------------- repointed_decoys

N_INTERGENIC = 1500
P_TRUE_IG = 0.4


def _repointed(loss: float) -> Callable[[np.random.Generator], np.ndarray]:
    levels = C.load_tolerances().levels

    def fdp(rng: np.random.Generator) -> np.ndarray:
        w = PickedWorld()
        false_ids: set[str] = set()
        n_lost = 0
        for k in range(N_INTERGENIC):
            contig = CONTIGS[k % len(CONTIGS)]
            true = rng.random() < P_TRUE_IG
            has_t = true or rng.random() < P_HIT
            has_d = rng.random() < P_HIT
            if not (has_t or has_d):
                continue
            entry = w.entry("intergenic", contig, "+", in_map=False)
            if has_t:
                s = _best(rng, 3.0 if true else 0.0)
                orf, _ = w.target("intergenic", contig, 0, entry=entry, q=_q_in(s), score=s)
                if not true:
                    false_ids.add(orf.orf_id)
            if has_d:
                s = _best(rng, 0.0)
                did = entry[0]
                if rng.random() < loss:  # tagged, parseable, no target anywhere
                    n_lost += 1
                    did = f"GIAPORF_intergenic_ghost_{600 * n_lost + 1}_{600 * n_lost + 300}_p"
                w.decoy(did, 0, q=_q_in(s), score=s)
        res = _cc().class_fdr(w.orfs, w.peptides, gm=GM, genome=GENOME, cfg=config(),
                              denovo=None)  # fmt: skip
        cells = {k: _cc().class_fdr_cell(r) for k, r in res.orfs.items()}
        out = np.zeros(len(levels))
        for j, a in enumerate(levels):
            called = {o for o, v in cells.items() if not isinstance(v, str) and v <= a}
            out[j] = len(called & false_ids) / len(called) if called else 0.0
        return out

    return fdp


@pytest.mark.parametrize("loss", [0.2, 0.5], ids=["loss20", "loss50"])
@pytest.mark.parametrize("full", MODES)
def test_repointed_decoys_paired_share_one_at_level(full: bool, loss: float) -> None:
    levels = C.load_tolerances().levels
    calibrate_levels(SUITE, "fdr", levels, _repointed(loss), C.replicates(SUITE, full=full),
                     stream=f"{SUITE}.attack.class_call.repointed_{int(loss * 100)}",
                     ).assert_ok()  # fmt: skip
