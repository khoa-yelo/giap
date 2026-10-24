"""Hypothesis properties of the M7-11b wiring (docs/specs/M7-11b-class-fdr-wiring.md rules 3, 5,
6 and 8):

* ``denovo_class_fdr(...).q`` equals ``denovo_class_q`` (its results do not change), and its
  per-class stats are the list sizes: distinct target peptides with a hit, best-cost decoy
  peptides of the class, their pi-hat sum (missing = 1) and c = the largest searched pi-hat;
* picked: every class q is >= its ``q_floor = 1 / n_target_wins``; ``floor_above_alpha`` is the
  float comparison ``q_floor > alpha``; a class with fewer than 1/alpha target wins calls nothing;
* de novo: every class q is >= ``c / n_targets``;
* ``reason is None`` iff ``q`` is finite, and ``class_fdr_cell`` is q, 1.0 (``lost_pick``) or
  ``"unestimated"``;
* the call rule: ``called = rule and q <= alpha`` for a numeric cell, ``rule and n_peptides >=
  2`` otherwise;
* ``class_fdr`` does not depend on ORF or peptide-row order.
"""

from __future__ import annotations

import math
import random
from typing import Any

import pandas as pd
from helpers.class_call_fixture import assemble, cc, config, gene_model, world_genome
from helpers.novel_fdr_fixture import ContigClassifier, PickedWorld, hit_rows, make_orf, pep_name
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

from giap.model.coords import Interval

SETTINGS = settings(max_examples=40, deadline=None,
                    suppress_health_check=[HealthCheck.too_slow])  # fmt: skip
SLOW = settings(max_examples=12, deadline=None,
                suppress_health_check=[HealthCheck.too_slow])  # fmt: skip
CLASSES = ("intergenic", "cds_outframe", "antisense")
GENOME = world_genome()
GM = gene_model(GENOME)


# ---------------------------------------------------------------- de novo


@st.composite
def denovo_case(draw: st.DrawFn) -> dict[str, Any]:
    n_t = draw(st.integers(1, 30))
    n_d = draw(st.integers(0, 15))
    e = st.floats(1e-6, 10.0, allow_nan=False)
    pi_v = st.one_of(st.none(), st.floats(0.01, 1.0, allow_nan=False))
    return {
        "t_cls": draw(st.lists(st.sampled_from(CLASSES), min_size=n_t, max_size=n_t)),
        "d_cls": draw(st.lists(st.sampled_from(CLASSES), min_size=n_d, max_size=n_d)),
        "e_t": draw(st.lists(e, min_size=n_t, max_size=n_t, unique=True)),
        "e_d": draw(st.lists(e, min_size=n_d, max_size=n_d)),
        "pi_t": draw(st.lists(pi_v, min_size=n_t, max_size=n_t)),
        "pi_d": draw(st.lists(pi_v, min_size=n_d, max_size=n_d)),
        "group": draw(st.lists(st.integers(0, 3), min_size=n_t, max_size=n_t)),
    }


def _dn_frames(case: dict[str, Any]) -> tuple[list[Any], pd.DataFrame, pd.DataFrame, Any, dict]:
    tp = [pep_name(k, "DNVO") for k in range(len(case["t_cls"]))]
    dp = [pep_name(k, "DCOY") for k in range(len(case["d_cls"]))]
    th = hit_rows(tp, [Interval(c, 100, 130, "-") for c in case["t_cls"]], case["e_t"],
                  search="target")  # fmt: skip
    dh = hit_rows(dp, [Interval(c, 100, 130, "-") for c in case["d_cls"]], case["e_d"],
                  search="decoy")  # fmt: skip
    pi = {p: v for p, v in zip([*tp, *dp], [*case["pi_t"], *case["pi_d"]], strict=True)
          if v is not None}  # fmt: skip
    slots: dict[tuple[str, int], list[str]] = {}
    for i, p in enumerate(tp):
        slots.setdefault((case["t_cls"][i], case["group"][i]), []).append(p)
    orfs = [make_orf(f"N_{c}_{g}", c, c, "-", 97, ps) for (c, g), ps in slots.items()]
    return orfs, th, dh, ContigClassifier({c: c for c in CLASSES}), pi


@SETTINGS
@given(denovo_case())
def test_denovo_class_fdr_q_equals_denovo_class_q(case: dict[str, Any]) -> None:
    from giap.novel import fdr

    orfs, th, dh, clf, pi = _dn_frames(case)
    want = fdr.denovo_class_q(orfs, th, dh, clf, pi)
    got = fdr.denovo_class_fdr(orfs, th, dh, clf, pi).q
    assert set(got) == set(want)
    for k, v in want.items():
        assert (math.isnan(v) and math.isnan(got[k])) or got[k] == v, k


@SETTINGS
@given(denovo_case())
def test_denovo_class_stats_are_the_list_sizes(case: dict[str, Any]) -> None:
    from giap.novel import fdr

    orfs, th, dh, clf, pi = _dn_frames(case)
    out = fdr.denovo_class_fdr(orfs, th, dh, clf, pi)
    w = [1.0 if v is None else v for v in [*case["pi_t"], *case["pi_d"]]]
    c = max(w)
    for cls in set(case["t_cls"]):
        s = out.classes[cls]
        assert s.n_targets == case["t_cls"].count(cls)
        di = [i for i, k in enumerate(case["d_cls"]) if k == cls]
        assert s.n_decoys == len(di)
        wd = [1.0 if case["pi_d"][i] is None else case["pi_d"][i] for i in di]
        assert math.isclose(s.w_decoys, sum(wd), rel_tol=1e-9, abs_tol=1e-12)
        assert math.isclose(s.c, c, rel_tol=1e-12)
        floor = s.c / s.n_targets
        for o in orfs:
            if o.cls == cls:
                assert out.q[o.orf_id] >= floor * (1 - 1e-12), (o.orf_id, out.q[o.orf_id], floor)


@SETTINGS
@given(denovo_case(), st.floats(0.001, 0.2))
def test_denovo_summary_floor_through_class_fdr(case: dict[str, Any], alpha: float) -> None:
    orfs, th, dh, clf, pi = _dn_frames(case)
    peps = _c9_denovo([p for o in orfs for p in o.peptides])
    d = cc().DenovoInputs(th[["peptide", "e_target"]], dh, pi, clf)
    res = cc().class_fdr(orfs, peps, gm=GM, genome=GENOME, cfg=config(alpha), denovo=d)
    for s in res.classes:
        assert s.regime == "denovo_subgroup"
        if s.n_targets == 0:
            assert s.q_floor is None
        else:
            assert s.q_floor is not None and math.isclose(s.q_floor, s.c / s.n_targets)
        assert s.floor_above_alpha is (s.q_floor is not None and s.q_floor > alpha)


def _c9_denovo(peptides: list[str]) -> pd.DataFrame:
    from helpers.novel_fdr_fixture import pep_row, peptides_frame

    return peptides_frame([pep_row(p, [], 0, is_decoy=False, source="denovo", q=math.nan)
                           for p in peptides])  # fmt: skip


# ---------------------------------------------------------------- picked


@st.composite
def picked_world(draw: st.DrawFn) -> tuple[PickedWorld, dict[str, str]]:
    """Up to three classes on their own contigs; each entry a target ORF (1-3 peptides) with an
    optional decoy, plus decoy-only entries; all ranks distinct. Returns the world and the
    ORF -> class map."""
    w = PickedWorld()
    n_cls = draw(st.integers(1, 3))
    specs = []
    for _ in range(n_cls):
        n = draw(st.integers(1, 25))
        specs.append([(draw(st.integers(1, 3)), draw(st.booleans())) for _ in range(n)])
    n_only = draw(st.integers(0, 4))
    total = sum(p + d for sp in specs for p, d in sp) + n_only + 1
    ranks = iter(draw(st.permutations(list(range(1, total + 1)))))
    cls_of = {}
    for k, sp in enumerate(specs):
        cls, contig = CLASSES[k], f"c{k + 1}"
        for n_p, has_d in sp:
            o, did = w.target(cls, contig, [next(ranks) for _ in range(n_p)])
            cls_of[o.orf_id] = cls
            if has_d:
                w.decoy(did, next(ranks))
    for _ in range(n_only):
        did, _ = w.entry("antisense", "c4", in_map=False)
        w.decoy(did, next(ranks))
    did, _ = w.entry("intergenic", "c4", in_map=False)  # always a makedb run
    w.decoy(did, next(ranks))
    return w, cls_of


@SETTINGS
@given(picked_world(), st.floats(0.001, 0.3))
def test_picked_q_never_below_the_floor(data: tuple[PickedWorld, dict], alpha: float) -> None:
    w, _ = data
    res = cc().class_fdr(w.orfs, w.peptides, gm=GM, genome=GENOME, cfg=config(alpha),
                         denovo=None)  # fmt: skip
    by = {(s.regime, s.cls): s for s in res.classes}
    for r in res.orfs.values():
        assert r.regime == "picked_decoy"
        if not math.isfinite(r.q):
            continue
        s = by[("picked_decoy", r.fdr_class)]
        assert s.n_target_wins and math.isclose(s.q_floor, 1 / s.n_target_wins)
        assert r.q >= s.q_floor * (1 - 1e-12)
        if s.n_target_wins < 1 / alpha * (1 - 1e-9):
            assert r.q > alpha  # fewer than about 1/alpha candidates: none callable
    for s in res.classes:
        if not s.estimated or not s.n_target_wins:
            assert s.q_floor is None
        else:
            assert s.q_floor is not None and math.isclose(s.q_floor, 1 / s.n_target_wins)
        assert s.floor_above_alpha is (s.q_floor is not None and s.q_floor > alpha)


@SETTINGS
@given(picked_world())
def test_reason_iff_not_finite_and_cell(data: tuple[PickedWorld, dict]) -> None:
    w, _ = data
    m = cc()
    res = m.class_fdr(w.orfs, w.peptides, gm=GM, genome=GENOME, cfg=config(), denovo=None)
    for r in res.orfs.values():
        assert (r.reason is None) is math.isfinite(r.q)
        cell = m.class_fdr_cell(r)
        if math.isfinite(r.q):
            assert cell == r.q and 0.0 < r.q <= 1.0
        elif r.reason == "lost_pick":
            assert cell == 1.0
        else:
            assert cell == "unestimated"


@SETTINGS
@given(picked_world(), st.randoms(use_true_random=False))
def test_class_fdr_is_order_free(data: tuple[PickedWorld, dict], rnd: random.Random) -> None:
    w, _ = data
    m = cc()
    a = m.class_fdr(w.orfs, w.peptides, gm=GM, genome=GENOME, cfg=config(), denovo=None)
    orfs = list(w.orfs)
    rnd.shuffle(orfs)
    idx = list(range(len(w.pep_rows)))
    rnd.shuffle(idx)
    peps = w.peptides.iloc[idx].reset_index(drop=True)
    b = m.class_fdr(orfs, peps, gm=GM, genome=GENOME, cfg=config(), denovo=None)
    assert a.classes == b.classes
    assert set(a.orfs) == set(b.orfs)
    for k, r in a.orfs.items():
        s = b.orfs[k]
        assert (r.regime, r.fdr_class, r.reason, r.class_mismatch) == (
            s.regime, s.fdr_class, s.reason, s.class_mismatch)  # fmt: skip
        assert (math.isnan(r.q) and math.isnan(s.q)) or r.q == s.q


@SLOW
@given(picked_world(), st.sampled_from([0.005, 0.01, 0.05, 0.2]))
def test_call_rule(data: tuple[PickedWorld, dict], alpha: float) -> None:
    w, _ = data
    assume(len(w.orfs) <= 40)
    a = assemble(w.orfs, w.peptides, gm=GM, genome=GENOME, cfg=config(alpha))
    for o in w.orfs:
        r = a.res.orfs[o.orf_id]
        rule = a.rules[o.orf_id]
        if math.isfinite(r.q):
            want = rule and r.q <= alpha
        elif r.reason == "lost_pick":
            want = False
        else:
            want = rule and len(o.peptides) >= 2
        assert a.called(o.orf_id) is want, (o.orf_id, o.cls, r, rule)
