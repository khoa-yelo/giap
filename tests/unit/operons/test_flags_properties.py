"""Hypothesis properties of M6-07 (docs/specs/M6-07-flags-and-operons-table.md, acceptance
criteria): the flag is a pure function of (class, promoter_yes, terminator_yes, pattern); every
``operon_flags.tsv`` row exists in ``pairs.tsv`` with the same ``q_value`` and never has pattern
``none``; ``sum(n_flags) <= flagged OP pairs``; one ``operons.tsv`` row per TU; no mutation,
determinism and invariance to the input row order.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pandas as pd
from helpers.flags_data import NAN, cfg, chain_pairs, frame, gene_table, row, tu, tus
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.operons import flags as fl
from gippa.operons import tu_table as tt

SETTINGS = settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
EVIDENCE = st.sampled_from([NAN, 0.1, 0.4999, 0.5, 0.81])
Q = st.sampled_from([NAN, 0.02, 0.1, 0.10001, 0.6])
ROW = st.builds(
    lambda cls, prom, term, det, delta, q, flags: row(
        cls, prom=prom, term=term, detected=det, delta=delta, q=q, lfdr=0.3, flags=flags
    ),
    st.sampled_from(["OP", "NO", "CV", "DV", "SW"]),
    EVIDENCE,
    EVIDENCE,
    st.sampled_from([(True, True), (True, False), (False, True), (False, False)]),
    st.sampled_from([-1.4, 0.0, 1.4]),
    Q,
    st.sampled_from(
        [[], ["unestimated"], ["ci_unestimated"], ["abundance_imputed"], ["rna_between"],
         ["layer_missing:terminators"]]
    ),
)  # fmt: skip
FRAMES = st.lists(ROW, min_size=0, max_size=25).map(lambda rs: frame(*rs))


def _same(a: float, b: float) -> bool:
    return (math.isnan(a) and math.isnan(b)) or a == b


@SETTINGS
@given(FRAMES)
def test_flag_rows_come_from_pairs_with_the_same_q(pairs: pd.DataFrame) -> None:
    out = fl.flag_pairs(pairs, cfg())
    assert tuple(out.columns) == tuple(fl.FLAG_COLUMNS)
    by = pairs.set_index("pair_id")
    for r in out.itertuples(index=False):
        src = by.loc[r.pair_id]
        assert _same(float(src["q_value"]), float(r.q_value))
        assert _same(float(src["lfdr"]), float(r.lfdr))
        assert (src["gene_a"], src["gene_b"], src["class"]) == (r.gene_a, r.gene_b, r.uniop_class)
    assert "none" not in set(out["abundance_pattern"])
    assert set(out["flag"]) <= set(fl.FLAGS)
    assert set(out["uniop_class"]) <= {"OP", "NO"}
    assert out["pair_id"].is_unique


@SETTINGS
@given(FRAMES)
def test_flag_is_a_function_of_class_evidence_and_pattern(pairs: pd.DataFrame) -> None:
    c = cfg()
    ev = fl.evidence_yes(pairs, c)
    pat = fl.abundance_pattern(pairs, c.stats.alpha.operon_flags)
    expected = {}
    for i, r in pairs.iterrows():
        f = fl.flag_for(
            r["class"],
            bool(ev.loc[i, "promoter_between_yes"]),
            bool(ev.loc[i, "terminator_between_yes"]),
            pat.loc[i],
        )
        if f is not None:
            expected[r["pair_id"]] = (f, pat.loc[i])
    out = fl.flag_pairs(pairs, c)
    got = {
        p: (f, a)
        for p, f, a in zip(out["pair_id"], out["flag"], out["abundance_pattern"], strict=True)
    }
    assert got == expected


@SETTINGS
@given(FRAMES)
def test_unestimated_rows_are_never_flagged(pairs: pd.DataFrame) -> None:
    out = fl.flag_pairs(pairs, cfg())
    unest = {p for p, f in zip(pairs["pair_id"], pairs["flags"], strict=True) if "unestimated" in f}
    assert not unest & set(out["pair_id"])


@SETTINGS
@given(FRAMES, st.randoms(use_true_random=False))
def test_pure_deterministic_and_order_free(pairs: pd.DataFrame, rnd) -> None:
    before = pairs.copy(deep=True)
    a = fl.flag_pairs(pairs, cfg())
    pd.testing.assert_frame_equal(pairs, before)
    perm = list(range(len(pairs)))
    rnd.shuffle(perm)
    b = fl.flag_pairs(pairs.iloc[perm].reset_index(drop=True), cfg())
    pd.testing.assert_frame_equal(a.reset_index(drop=True), b.reset_index(drop=True))
    ids = [int(p.split("_")[1]) for p in a["pair_id"]]
    assert ids == sorted(ids)


# ------------------------------------------------------------------------------ operons.tsv


@st.composite
def layouts(draw):
    sizes = draw(st.lists(st.integers(1, 5), min_size=1, max_size=8))
    tags, k = [], 0
    for s in sizes:
        tags.append([f"g{k + j}" for j in range(s)])
        k += s
    flat = [g for gs in tags for g in gs]
    ab = {g: draw(st.floats(5.0, 25.0, allow_nan=False)) for g in flat if draw(st.booleans())}
    names = list(itertools.pairwise(flat))
    q = {p: draw(st.sampled_from([NAN, 0.02, 0.6])) for p in names}
    prom = {p: draw(st.sampled_from([0.1, 0.81])) for p in names}
    term = {p: draw(st.sampled_from([0.1, 0.8])) for p in names}
    return tags, ab, q, prom, term


@SETTINGS
@given(layouts())
def test_tu_table_properties(layout) -> None:
    tu_genes, ab, q, prom, term = layout
    rows, pos = [], 0
    for k, gs in enumerate(tu_genes):
        rows.append(tu(f"TU_{k + 1}", gs, pos, pos + 300 * len(gs)))
        pos += 300 * len(gs) + 50
    t = tus(*rows)
    pairs = chain_pairs(tu_genes, ab, q=q, prom=prom, term=term)
    genes = gene_table(ab, [g for gs in tu_genes for g in gs])
    flagged = fl.flag_pairs(pairs, cfg())
    out = tt.operons_table(t, pairs, genes, flagged)

    assert list(out["tu_id"]) == [f"TU_{k + 1}" for k in range(len(tu_genes))]
    assert [list(g) for g in out["genes"]] == [list(g) for g in tu_genes]
    assert list(out["n_genes"]) == [len(g) for g in tu_genes]
    n_op = int((flagged["uniop_class"] == "OP").sum())
    assert int(out["n_flags"].sum()) <= n_op
    # in this layout every OP pair lies inside one TU, so the bound is reached
    assert int(out["n_flags"].sum()) == n_op
    for r, gs in zip(out.itertuples(index=False), tu_genes, strict=True):
        vals = [ab[g] for g in gs if g in ab]
        assert r.detected_frac == len(vals) / len(gs)
        if vals:
            assert r.median_abundance_log2 == np.median(vals)
        else:
            assert math.isnan(r.median_abundance_log2)
        deltas = [abs(ab[b] - ab[a]) for a, b in itertools.pairwise(gs) if a in ab and b in ab]
        if deltas:
            assert math.isclose(r.max_abs_delta_log2, max(deltas), rel_tol=1e-12, abs_tol=1e-12)
        else:
            assert math.isnan(r.max_abs_delta_log2)
