"""Hypothesis properties of the M7-11 class FDR (docs/specs/M7-11-novel-class-fdr.md rules 2-3).

* picked, one class, every decoy paired: the class q equals ``gippa.stats.tda.picked_fdr`` on the
  per-entry best scores (oracle; ADR-0009 rescaling is the identity at f = 1);
* picked classes are never pooled: a class's q does not change when another class is added;
* q lies in (0, 1], is monotone in the ORF score, and does not depend on row or ORF order;
* de novo, per class: single-peptide ORF q equals ``stratum_qvalues`` on that class's list with
  c = the largest searched π̂ (oracle), and a multi-peptide ORF gets the minimum (decision 3).
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pandas as pd
from helpers.novel_fdr_fixture import (
    ContigClassifier,
    PickedWorld,
    hit_rows,
    make_orf,
    pep_name,
    q_of,
)
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

from gippa.mapping.placement_fdr import stratum_qvalues
from gippa.model.coords import Interval
from gippa.stats.tda import picked_fdr

SETTINGS = settings(max_examples=40, deadline=None,
                    suppress_health_check=[HealthCheck.too_slow])  # fmt: skip


def _fdr():  # type: ignore[no-untyped-def]
    import gippa.novel.fdr as m  # fails (ModuleNotFoundError) until M7-11 lands

    return m


@st.composite
def picked_class(draw: st.DrawFn, max_n: int = 30) -> list[tuple[int, int | None]]:
    """Per entry: (target rank, decoy rank or None); all ranks distinct, every entry has a target
    ORF, at least one decoy."""
    n = draw(st.integers(1, max_n))
    has_d = draw(st.lists(st.booleans(), min_size=n, max_size=n))
    assume(any(has_d))
    ranks = draw(st.permutations(list(range(1, n + sum(has_d) + 1))))
    it = iter(ranks)
    return [(next(it), next(it) if d else None) for d in has_d]


def _build(w: PickedWorld, cls: str, contig: str, spec: list[tuple[int, int | None]]) -> list[str]:
    ids = []
    for t, d in spec:
        orf, did = w.target(cls, contig, t)
        if d is not None:
            w.decoy(did, d)
        ids.append(orf.orf_id)
    return ids


@SETTINGS
@given(picked_class())
def test_picked_matches_oracle(spec: list[tuple[int, int | None]]) -> None:
    w = PickedWorld()
    ids = _build(w, "intergenic", "c1", spec)
    qs = _fdr().picked_class_q(w.orfs, w.peptides, w.id_map)
    t = np.array([-r for r, _ in spec], dtype=float)
    d = np.array([np.nan if r is None else -r for _, r in spec], dtype=float)
    res = picked_fdr(t, d)
    for o, q, lost in zip(ids, res.q, res.decoy_wins, strict=True):
        if lost:
            assert math.isnan(q_of(qs, o))  # lost the pick: unestimated, never a class q
        else:
            assert math.isclose(q_of(qs, o), q, rel_tol=1e-12)


@SETTINGS
@given(picked_class(), picked_class())
def test_picked_classes_never_pooled(
    a: list[tuple[int, int | None]], b: list[tuple[int, int | None]]
) -> None:
    w1 = PickedWorld()
    ids = _build(w1, "intergenic", "c1", a)
    alone = _fdr().picked_class_q(w1.orfs, w1.peptides, w1.id_map)
    _build(w1, "antisense", "c3", b)
    both = _fdr().picked_class_q(w1.orfs, w1.peptides, w1.id_map)
    for o in ids:
        qa, qb = q_of(alone, o), q_of(both, o)
        assert (math.isnan(qa) and math.isnan(qb)) or qa == qb


@SETTINGS
@given(picked_class(), st.randoms(use_true_random=False))
def test_picked_range_monotone_and_order_free(
    spec: list[tuple[int, int | None]], rnd: object
) -> None:
    w = PickedWorld()
    ids = _build(w, "cds_outframe", "c2", spec)
    peptides, idm = w.peptides, w.id_map
    qs = _fdr().picked_class_q(w.orfs, peptides, idm)
    rank = {o: t for o, (t, _) in zip(ids, spec, strict=True)}
    got = [(rank[o], q_of(qs, o)) for o in ids if not math.isnan(q_of(qs, o))]
    for _, q in got:
        assert 0.0 < q <= 1.0
    got.sort()
    assert all(x[1] <= y[1] for x, y in itertools.pairwise(got))  # better rank, q <=
    perm = list(range(len(peptides)))
    rnd.shuffle(perm)  # type: ignore[attr-defined]
    orfs = list(w.orfs)
    rnd.shuffle(orfs)  # type: ignore[attr-defined]
    qs2 = _fdr().picked_class_q(orfs, peptides.iloc[perm].reset_index(drop=True),
                                idm.iloc[::-1].reset_index(drop=True))  # fmt: skip
    for o in ids:
        q1, q2 = q_of(qs, o), q_of(qs2, o)
        assert (math.isnan(q1) and math.isnan(q2)) or q1 == q2


# ---------------------------------------------------------------- de novo

CLASSES = ("intergenic", "cds_outframe", "antisense")


@st.composite
def denovo_case(draw: st.DrawFn) -> dict[str, object]:
    n_t = draw(st.integers(1, 40))
    n_d = draw(st.integers(0, 20))
    t_cls = draw(st.lists(st.sampled_from(CLASSES), min_size=n_t, max_size=n_t))
    d_cls = draw(st.lists(st.sampled_from(CLASSES), min_size=n_d, max_size=n_d))
    e = st.floats(1e-6, 10.0, allow_nan=False)
    e_t = draw(st.lists(e, min_size=n_t, max_size=n_t, unique=True))
    e_d = draw(st.lists(e, min_size=n_d, max_size=n_d))
    pi_v = st.one_of(st.none(), st.floats(0.01, 1.0, allow_nan=False))
    pi_t = draw(st.lists(pi_v, min_size=n_t, max_size=n_t))
    pi_d = draw(st.lists(pi_v, min_size=n_d, max_size=n_d))
    group = draw(st.lists(st.integers(0, 3), min_size=n_t, max_size=n_t))  # ORF slot per peptide
    return {"t_cls": t_cls, "d_cls": d_cls, "e_t": e_t, "e_d": e_d, "pi_t": pi_t,
            "pi_d": pi_d, "group": group}  # fmt: skip


def _dn(case: dict[str, object], grouped: bool) -> tuple[dict[str, float], dict[str, float]]:
    """(ORF q, oracle peptide q)."""
    t_cls, d_cls = case["t_cls"], case["d_cls"]
    e_t, e_d = case["e_t"], case["e_d"]
    pi_t, pi_d, group = case["pi_t"], case["pi_d"], case["group"]
    assert isinstance(t_cls, list) and isinstance(d_cls, list) and isinstance(group, list)
    assert isinstance(e_t, list) and isinstance(e_d, list)
    assert isinstance(pi_t, list) and isinstance(pi_d, list)
    clf = ContigClassifier({c: c for c in CLASSES})
    tp = [pep_name(k, "DNVO") for k in range(len(t_cls))]
    dp = [pep_name(k, "DCOY") for k in range(len(d_cls))]
    targets = hit_rows(tp, [Interval(c, 100, 130, "-") for c in t_cls], e_t, search="target")
    decoys = hit_rows(dp, [Interval(c, 100, 130, "-") for c in d_cls], e_d, search="decoy")
    pi = {p: v for p, v in zip([*tp, *dp], [*pi_t, *pi_d], strict=True) if v is not None}
    w = [1.0 if v is None else v for v in [*pi_t, *pi_d]]
    c = max(w)
    oracle: dict[str, float] = {}
    for cls in CLASSES:
        ti = [i for i, k in enumerate(t_cls) if k == cls]
        di = [i for i, k in enumerate(d_cls) if k == cls]
        q = stratum_qvalues([e_t[i] for i in ti], [e_d[i] for i in di],
                            [w[len(tp) + i] for i in di], c)  # fmt: skip
        oracle.update({tp[i]: float(v) for i, v in zip(ti, q, strict=True)})
    slots: dict[tuple[str, int], list[str]] = {}
    for i, p in enumerate(tp):
        slots.setdefault((t_cls[i], group[i] if grouped else i), []).append(p)
    orfs = [make_orf(f"N_{cls}_{g}", cls, cls, "-", 97, ps) for (cls, g), ps in slots.items()]
    qs = _fdr().denovo_class_q(orfs, targets, decoys, clf, pi)
    by_orf = {o.orf_id: min(oracle[p] for p in o.peptides) for o in orfs}
    return qs, by_orf


@SETTINGS
@given(denovo_case())
def test_denovo_single_peptide_matches_stratum_oracle(case: dict[str, object]) -> None:
    qs, oracle = _dn(case, grouped=False)
    assert set(oracle) <= set(qs)
    for o, q in oracle.items():
        assert math.isclose(q_of(qs, o), q, rel_tol=1e-9, abs_tol=1e-12)
        assert 0.0 < q_of(qs, o) <= 1.0


@SETTINGS
@given(denovo_case())
def test_denovo_orf_q_is_min_over_peptides(case: dict[str, object]) -> None:
    qs, oracle = _dn(case, grouped=True)
    for o, q in oracle.items():
        assert math.isclose(q_of(qs, o), q, rel_tol=1e-9, abs_tol=1e-12)


def test_property_helpers_build_frames() -> None:
    """Guard: the builders the properties use produce the pinned layouts."""
    w = PickedWorld()
    _build(w, "intergenic", "c1", [(1, 2)])
    assert isinstance(w.peptides, pd.DataFrame) and w.peptides["is_decoy"].tolist() == [False, True]
    assert w.id_map["db_id"].str.startswith("GIPPAORF_intergenic_c1_").all()
