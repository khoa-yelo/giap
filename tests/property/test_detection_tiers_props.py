"""Hypothesis properties of M4-12 detection tiers (docs/specs/M4-12-detection-tiers.md, acceptance
criterion 2; GIPPA_SPEC.md §7.3.2).

- adding a decoy never turns a ``not_detected`` gene into ``detected``. Scope (ADR-0009): every
  decoy is paired before and after, so the paired share f stays 1 and the tier stays
  ``decoy_picked``; with unpaired decoys an extra paired decoy raises f and may lawfully lower q;
- ``detected`` implies ``n_accepted >= 1`` (M4-11 evidence), in every tier;
- the tier and the detection table are functions of the input content only: permuting the rows of
  C9, C11 and C12 changes neither;
- statuses: partial genes are ``partial_excluded``; q is NA outside the q-value tiers and in
  [0, 1] inside them; ``detected`` iff status ``detected``.

Scenarios are drawn small (<= 8 genes) and run through the merged producers
(tests/helpers/detection_tiers_fixtures.py).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import detection_tiers_fixtures as fx
from detection_tiers_fixtures import Builder, cfg

FAST = settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
Q = st.sampled_from([0.0, 0.0005, 0.001, 0.005, 0.01, 0.05])  # 0.05 > q_value_max: rejected
PEP = st.sampled_from([0.001, 0.01, 0.05, 0.2])
SCORE = st.sampled_from([10.0, 20.0, 30.0, 40.0])
ALPHA = st.sampled_from([0.01, 0.1, 0.5, 1.0])


@st.composite
def scenario(draw: st.DrawFn, mode: str | None = None) -> dict[str, Any]:
    """Gene count, partial genes, targets (gene, q, pep, score), paired decoys (gene, q, pep,
    score), unpaired decoys, C12 rank-1 hits; ``mode`` picks the evidence kind."""
    n = draw(st.integers(1, 8))
    kind = mode or draw(st.sampled_from(["scored", "pep_only", "sequences", "proteins"]))
    partial = draw(st.sets(st.integers(0, n - 1), max_size=2))
    targets = draw(st.lists(st.tuples(st.integers(0, n - 1), Q, PEP, SCORE), max_size=12))
    decoys = draw(st.lists(st.tuples(st.integers(0, n - 1), Q, PEP, SCORE), max_size=6))
    unpaired = draw(st.integers(0, 2)) if kind == "scored" else 0
    prots = draw(st.sets(st.integers(0, n - 1), max_size=n))
    return {"n": n, "kind": kind, "partial": partial, "targets": targets, "decoys": decoys,
            "unpaired": unpaired, "prots": prots, "alpha": draw(ALPHA)}  # fmt: skip


def build(
    sc: dict[str, Any], extra_decoy: tuple[int, float, float, float] | None = None
) -> Builder:
    b = Builder(sc["n"], partial={fx.tag(i) for i in sc["partial"]})
    kind = sc["kind"]
    for g, q, pep, score in sc["targets"]:
        if kind == "scored":
            b.target(fx.tag(g), q, pep, score)
        elif kind == "pep_only":
            b.target(fx.tag(g), None, pep, score)
        elif kind == "sequences":
            b.target(fx.tag(g), None, None, None)
    if kind == "scored":
        for g, q, pep, score in [*sc["decoys"], *([extra_decoy] if extra_decoy else [])]:
            b.decoy([f"rev_{fx.tag(g)}"], q, pep, score)
        for i in range(sc["unpaired"]):
            b.decoy([f"rev_UNKNOWN{i}"], 0.001, 0.01, 10.0)
    for g in sorted(sc["prots"]):
        b.protein_hit(f"UP_{g}", fx.tag(g))
    return b


def _run(sc: dict[str, Any], **kw: Any) -> fx.Run:
    return build(sc, kw.pop("extra_decoy", None)).run(
        cfg(stats={"alpha": {"identification": sc["alpha"]}}), no_peptides=sc["kind"] == "proteins",
        **kw,
    )  # fmt: skip


@FAST
@given(sc=scenario("scored"), extra=st.tuples(st.integers(0, 7), Q, PEP, SCORE))
def test_adding_a_decoy_never_creates_a_detection(
    sc: dict[str, Any], extra: tuple[int, float, float, float]
) -> None:
    sc = {**sc, "unpaired": 0}
    extra = (extra[0] % sc["n"], *extra[1:])
    before = _run(sc)
    after = _run(sc, extra_decoy=extra)
    if before.det.tier != "decoy_picked":
        return
    assert after.det.tier == "decoy_picked"
    b = before.det.table.set_index("locus_tag")
    a = after.det.table.set_index("locus_tag")
    newly = a["detected"].astype(bool) & ~b["detected"].astype(bool)
    assert not newly.any(), list(newly[newly].index)
    qb = b["gene_q_value"].astype(float)
    qa = a["gene_q_value"].astype(float)
    both = qb.notna() & qa.notna()
    assert (qa[both] >= qb[both] - 1e-12).all(), "an extra paired decoy never lowers a gene q"


@FAST
@given(sc=scenario())
def test_detected_implies_accepted_evidence_and_status_rules(sc: dict[str, Any]) -> None:
    r = _run(sc)
    t = r.det.table
    assert t["locus_tag"].astype(str).tolist() == r.genes["locus_tag"].astype(str).tolist()
    ev = r.evidence.set_index("locus_tag")
    det = t["detected"].astype(bool).to_numpy()
    status = t["detection_status"].astype(str).to_numpy()
    assert (det == (status == "detected")).all()
    assert set(status) <= set(fx.STATUSES)
    partial = r.genes["partial"].astype(bool).to_numpy()
    assert (status[partial] == "partial_excluded").all()
    assert (status[~partial] != "partial_excluded").all()
    q = t["gene_q_value"].astype(float).to_numpy()
    assert np.isnan(q[partial]).all()
    if r.det.tier in ("two_peptide", "protein_list"):
        assert np.isnan(q).all(), "unestimated tiers carry no gene q"
    else:
        ok = ~np.isnan(q)
        assert ((q[ok] >= 0) & (q[ok] <= 1)).all()
        assert (det == (ok & (q <= sc["alpha"]) & ~partial)).all()
    if r.det.tier != "protein_list":
        n_acc = ev.loc[t["locus_tag"].astype(str), "n_accepted"].to_numpy()
        assert (n_acc[det] >= 1).all()
    if r.det.tier == "two_peptide":
        single = status == "single_peptide_unestimated"
        n_acc = ev.loc[t["locus_tag"].astype(str), "n_accepted"].to_numpy()
        assert (n_acc[single] == 1).all() and not det[single].any()


@FAST
@given(sc=scenario(), seed=st.integers(0, 2**32 - 1))
def test_row_permutation_invariant(sc: dict[str, Any], seed: int) -> None:
    a = _run(sc)
    b = _run(sc, shuffle=np.random.default_rng(seed))
    assert a.det.tier == b.det.tier
    assert sorted(a.det.warnings) == sorted(b.det.warnings)
    ta, tb = a.det.table, b.det.table
    assert ta["detection_status"].tolist() == tb["detection_status"].tolist()
    assert ta["detected"].astype(bool).tolist() == tb["detected"].astype(bool).tolist()
    qa = ta["gene_q_value"].astype(float).to_numpy()
    qb = tb["gene_q_value"].astype(float).to_numpy()
    np.testing.assert_allclose(qa, qb, rtol=1e-12, equal_nan=True)


@FAST
@given(sc=scenario("scored"))
def test_tier_depends_on_pairing_share(sc: dict[str, Any]) -> None:
    """The tier follows the §7.3.2 table from the input columns: decoy_picked iff >= 1 decoy with
    q <= q_value_max carries protein_ids and >= 80 % of those pair. Without any target q the
    decoy q filter is off (M4-24 rule 6, ADR-0047)."""
    r = _run(sc)
    passing = [d for d in sc["decoys"] if not sc["targets"] or d[1] <= 0.01]
    n_pass = len(passing) + sc["unpaired"]
    expect_picked = n_pass > 0 and len(passing) / n_pass >= 0.8
    assert (r.det.tier == "decoy_picked") == expect_picked, (r.det.tier, len(passing), n_pass)
    if not expect_picked:
        assert r.det.tier == "mayu" or not any(t[1] <= 0.01 for t in sc["targets"])
    assert not math.isnan(fx.tier_metric(r.det, "n_decoys_paired"))
