"""M7-11 picked ranking key per class (code review #237 B2; lead decisions 11, 12, 13).

Decision 11: a class is ranked on (q, score) only when every target and decoy row has a finite q;
otherwise on score alone when every row has a score; otherwise it is unestimated with reason
``decoys_unscored``. Decision 12: inside the all-q branch, score breaks q ties only when every
decoy row has a score (else q alone, ties to the decoy); an unscored target then ranks last on
score. Each world below has two makedb entries of one class, every row at the same input q
(a single q plateau) unless the case needs a missing q, so the outcome is decided by the
tie-break rule under test. Review mutants killed: (a) C4, (b) C3 and C5, (c) C9,
(d) C8, C13 and C14.
"""

from __future__ import annotations

import math

from helpers.novel_fdr_fixture import PickedWorld

from giap.novel.fdr import PickedClassFdr, picked_class_fdr

NAN = math.nan
Q = 1e-3  # the shared q plateau
CLS = "intergenic"


def _run(w: PickedWorld) -> PickedClassFdr:
    return picked_class_fdr(w.orfs, w.peptides, w.id_map)


def test_a_all_q_finite_some_decoys_unscored_ranks_on_q_alone_ties_to_decoy() -> None:
    """(a) Every q is equal and finite; one decoy has no score. The class ranks on q alone, so
    both targets tie their decoys and lose, although each target outscores a decoy that has a
    score (score would break the tie toward the targets)."""
    w = PickedWorld()
    t1, d1 = w.target(CLS, "c1", 1, q=Q, score=900.0)
    w.decoy(d1, 1, q=Q, score=NAN)  # the unscored decoy
    t2, d2 = w.target(CLS, "c1", 1, q=Q, score=800.0)
    w.decoy(d2, 1, q=Q, score=500.0)
    res = _run(w)
    assert res.reason == {t1.orf_id: "lost_pick", t2.orf_id: "lost_pick"}
    assert math.isnan(res.q[t1.orf_id]) and math.isnan(res.q[t2.orf_id])
    st = res.classes[CLS]
    assert (st.n_target_wins, st.n_decoy_wins, st.estimated) == (0, 2, True)


def test_b_all_decoys_scored_unscored_target_ranks_last_class_uses_score() -> None:
    """(b) Every q is equal and finite, every decoy has a score, one target has none. The class
    ranks on (q, score): the scored target beats its lower-scored decoy, and the unscored target
    ranks last in the q plateau, so it loses to its decoy."""
    w = PickedWorld()
    t1, d1 = w.target(CLS, "c1", 1, q=Q, score=NAN)  # the unscored target
    w.decoy(d1, 1, q=Q, score=100.0)
    t2, d2 = w.target(CLS, "c1", 1, q=Q, score=900.0)
    w.decoy(d2, 1, q=Q, score=100.0)
    res = _run(w)
    assert res.reason == {t1.orf_id: "lost_pick"}
    # one target win above one decoy win (ranks: t2 > both decoys > t1): q = (0 + 1) / 1
    assert res.q[t2.orf_id] == 1.0
    st = res.classes[CLS]
    assert (st.n_target_wins, st.n_decoy_wins, st.estimated) == (1, 1, True)


def test_b_unscored_target_never_beats_a_scored_decoy_in_its_q_plateau() -> None:
    """(b) The conservative direction alone: a lone unscored target against a scored decoy at the
    same q loses (an unscored target ranked first would win)."""
    w = PickedWorld()
    t1, d1 = w.target(CLS, "c1", 1, q=Q, score=NAN)
    w.decoy(d1, 1, q=Q, score=1.0)
    t2, d2 = w.target(CLS, "c1", 1, q=Q, score=2.0)
    w.decoy(d2, 1, q=Q, score=1.0)
    res = _run(w)
    assert res.reason.get(t1.orf_id) == "lost_pick"
    assert res.reason.get(t2.orf_id) is None


def test_c_some_q_missing_every_score_present_higher_score_wins() -> None:
    """(c) One target has no q, every row has a score: the class ranks on score alone, and the
    higher score wins in both directions."""
    w = PickedWorld()
    t1, d1 = w.target(CLS, "c1", 1, q=NAN, score=900.0)  # no q, high score: wins
    w.decoy(d1, 1, q=Q, score=100.0)
    t2, d2 = w.target(CLS, "c1", 1, q=Q, score=100.0)  # low score: loses
    w.decoy(d2, 1, q=Q, score=900.0)
    res = _run(w)
    assert res.reason == {t2.orf_id: "lost_pick"}
    # ranks: t1 = d2 > t2 = d1; the d2 win ties t1 and counts against it: q = min(1, 2 / 1)
    assert res.q[t1.orf_id] == 1.0
    st = res.classes[CLS]
    assert (st.n_target_wins, st.n_decoy_wins, st.estimated) == (1, 1, True)


def test_d_some_q_missing_some_scores_missing_is_decoys_unscored() -> None:
    """(d) A target has no q and a decoy has no score: no key covers every row, so the class is
    unestimated (``decoys_unscored``) and so are its ORFs."""
    w = PickedWorld()
    t1, d1 = w.target(CLS, "c1", 1, q=NAN, score=900.0)
    w.decoy(d1, 1, q=Q, score=NAN)
    t2, d2 = w.target(CLS, "c1", 1, q=Q, score=900.0)
    w.decoy(d2, 1, q=Q, score=100.0)
    res = _run(w)
    st = res.classes[CLS]
    assert st.estimated is False
    assert st.reason == "decoys_unscored"
    assert (st.n_target_wins, st.n_decoy_wins) == (0, 0)
    assert res.reason == {t1.orf_id: "decoys_unscored", t2.orf_id: "decoys_unscored"}
    assert all(math.isnan(v) for v in res.q.values())
