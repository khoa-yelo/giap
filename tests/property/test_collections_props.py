"""Properties of the M8-03 gene-set collections (docs/specs/M8-03-function-collections.md).

Random gene tables over the synthetic vocabulary of tests/helpers/function_resources_data.py
(live GO terms, an alt ID, an obsolete and an absent ID; KOs with pathways, overview-only and
unknown KOs; COG letter strings with ``-``/blank/NA), random partial flags and ``p_detect``:

- P1 propagation is closed: every gene of a GO set is in the set of each of its non-root parents;
- P2 GO, slim, KEGG and COG sets equal independent oracles (is_a + part_of ancestors, roots
  removed; ``path:map`` rows only, overview maps excluded; one set per category letter);
- P3 no set contains a non-eligible gene and no set is empty;
- P4 coverage = (# eligible genes in >= 1 set) / (# eligible genes);
- P5 the result does not depend on the gene-table row order.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.function_resources_data import (
    ABSENT,
    B2_ALT,
    COG_FUN,
    GO_COLLECTION,
    KEGG_LINKS,
    KEGG_OVERVIEW,
    OBSOLETE,
    ROOTS,
    SLIM_TERMS,
    TERMS,
    eligible_mask,
    genes_frame,
    propagated,
    write_resources,
)
from hypothesis import given, settings
from hypothesis import strategies as st

GO_IDS = [*sorted(TERMS), B2_ALT, OBSOLETE, ABSENT]
KOS = [*sorted({k for k, _ in KEGG_LINKS}), "K99999"]
COG_VALUES = ["K", "KL", "J", "R", "S", "M", "LM", "-", "", None]


@pytest.fixture(scope="module")
def res(tmp_path_factory: pytest.TempPathFactory) -> Any:
    from giap.io.function_resources import read_resources

    root: Path = tmp_path_factory.mktemp("function_resources")
    return read_resources(write_resources(root))


_row = st.fixed_dictionaries(
    {
        "partial": st.booleans(),
        "p_detect": st.sampled_from([0.1, 0.5, 0.9, float("nan")]),
        "go_terms": st.lists(st.sampled_from(GO_IDS), max_size=4, unique=True),
        "kegg_ko": st.lists(st.sampled_from(KOS), max_size=3, unique=True),
        "cog_category": st.sampled_from(COG_VALUES),
    }
)
_tables = st.lists(_row, min_size=0, max_size=12)


def _frame(rows: list[dict[str, Any]]) -> pd.DataFrame:
    return genes_frame([{**r, "locus_tag": f"t{i:02d}"} for i, r in enumerate(rows)])


def _build(genes: pd.DataFrame, res: Any) -> dict[str, Any]:
    from giap.function.collections import build_collections

    out, _ = build_collections(genes, eligible_mask(genes), res)
    return {g.collection: g for g in out}


def _oracle(genes: pd.DataFrame) -> dict[str, dict[str, set[str]]]:
    el = eligible_mask(genes)
    links: dict[str, set[str]] = {}
    for ko, m in KEGG_LINKS:
        if m not in KEGG_OVERVIEW:
            links.setdefault(ko, set()).add(m)
    out: dict[str, dict[str, set[str]]] = {
        c: {} for c in ("cog", "kegg", "go_bp", "go_mf", "go_cc", "go_slim")
    }
    slim = set(SLIM_TERMS) - ROOTS
    for i, row in enumerate(genes.itertuples(index=False)):
        if not el[i]:
            continue
        tag = row.locus_tag
        prop = propagated(row.go_terms)
        for t in prop:
            out[GO_COLLECTION[TERMS[t][1]]].setdefault(t, set()).add(tag)
        for t in prop & slim:
            out["go_slim"].setdefault(t, set()).add(tag)
        for ko in row.kegg_ko:
            for m in links.get(ko, ()):
                out["kegg"].setdefault(m, set()).add(tag)
        if isinstance(row.cog_category, str):
            for ch in row.cog_category:
                if ch in COG_FUN:
                    out["cog"].setdefault(ch, set()).add(tag)
    return out


@settings(max_examples=60, deadline=None)
@given(rows=_tables)
def test_p1_propagation_closed(res: Any, rows: list[dict[str, Any]]) -> None:
    by = _build(_frame(rows), res)
    for c in GO_COLLECTION.values():
        gs = by[c]
        for t, members in gs.sets.items():
            for p in gs.parents[t] - ROOTS:
                assert members <= gs.sets[p], (c, t, p)


@settings(max_examples=60, deadline=None)
@given(rows=_tables)
def test_p2_p3_sets_equal_oracle(res: Any, rows: list[dict[str, Any]]) -> None:
    g = _frame(rows)
    by = _build(g, res)
    want = _oracle(g)
    assert set(by) == set(want)
    ineligible = set(g["locus_tag"][~eligible_mask(g)])
    for c, gs in by.items():
        got = {k: set(v) for k, v in gs.sets.items()}
        assert got == want[c], c
        for v in gs.sets.values():
            assert v and not (v & ineligible)


@settings(max_examples=60, deadline=None)
@given(rows=_tables)
def test_p4_coverage(res: Any, rows: list[dict[str, Any]]) -> None:
    g = _frame(rows)
    by = _build(g, res)
    want = _oracle(g)
    n = int(eligible_mask(g).sum())
    for c, gs in by.items():
        if n == 0:
            assert math.isnan(gs.coverage)
            continue
        covered = set().union(*want[c].values()) if want[c] else set()
        assert gs.coverage == pytest.approx(len(covered) / n), c
        assert 0.0 <= gs.coverage <= 1.0


@settings(max_examples=40, deadline=None)
@given(rows=_tables, data=st.data())
def test_p5_row_order_invariant(res: Any, rows: list[dict[str, Any]], data: Any) -> None:
    g = _frame(rows)
    order = data.draw(st.permutations(list(range(len(g)))))
    h = g.iloc[list(order)].reset_index(drop=True)
    a, b = _build(g, res), _build(h, res)
    assert set(a) == set(b)
    for c in a:
        assert dict(a[c].sets) == dict(b[c].sets)
        assert list(a[c].sets) == list(b[c].sets)
        assert a[c].n_unknown_ids == b[c].n_unknown_ids
        ca, cb = a[c].coverage, b[c].coverage
        assert (math.isnan(ca) and math.isnan(cb)) or ca == cb
