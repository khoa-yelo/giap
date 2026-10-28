"""M8-03: function resources -> gene-set collections (docs/specs/M8-03-function-collections.md).

Examples 1-10 and every degenerate row of the spec, on synthetic resources written by the tests
(tests/helpers/function_resources_data.py): a toy go-basic OBO (is_a + part_of, an alt ID, an
obsolete term), a toy slim OBO, a synthetic two-pathway KEGG table (``path:map`` rows, their
``path:ko`` duplicates, overview maps) and a toy COG fun-20 table. Never real KEGG content
(M8 index, lead decision 4).

Pinned interface (spec "Public interface"): ``gippa.function.collections.build_collections(genes,
eligible, res) -> (list[GeneSets], run_flags)``, ``GeneSets`` fields ``collection``, ``sets``,
``names``, ``parents``, ``coverage``, ``n_unknown_ids``; ``gippa.io.function_resources.
read_resources(paths) -> Resources`` keyed by the db_sources.yaml item names (``go_basic``,
``goslim_generic``, ``kegg_ko_pathway``, ``kegg_pathway_names``, ``cog2020_def``,
``cog2020_fun``); ``gippa.genetable.eligible(genes)`` (the shared §7.3 helper).
"""

from __future__ import annotations

import copy
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.function_resources_data import (
    ABSENT,
    B1,
    B2,
    B2_ALT,
    B3,
    B4,
    B5,
    BP_ROOT,
    C1,
    C2,
    CC_ROOT,
    DIRECT_PARENTS,
    KEGG_OVERVIEW,
    KEGG_PATHWAYS,
    M1,
    M2,
    MF_ROOT,
    OBSOLETE,
    ROOTS,
    TERMS,
    eligible_mask,
    genes_frame,
    write_resources,
)

COLLECTIONS = ("cog", "kegg", "go_bp", "go_mf", "go_cc", "go_slim")
GO = ("go_bp", "go_mf", "go_cc")


def _mods() -> tuple[Any, Any]:
    from gippa.function import collections as col
    from gippa.io import function_resources as fr

    return col, fr


def _res(tmp: Path, omit: tuple[str, ...] = ()) -> Any:
    _, fr = _mods()
    return fr.read_resources(write_resources(tmp, omit))


def _build(genes: pd.DataFrame, res: Any, eligible: np.ndarray | None = None
           ) -> tuple[dict[str, Any], list[str]]:  # fmt: skip
    col, _ = _mods()
    el = eligible_mask(genes) if eligible is None else eligible
    out, flags = col.build_collections(genes, el, res)
    by = {g.collection: g for g in out}
    assert len(by) == len(out), "one GeneSets per collection"
    return by, list(flags)


def _has_flag(flags: list[str], word: str) -> bool:
    return any(f == word or f.startswith(word + ":") for f in flags)


def _members(gs: Any, tag: str) -> set[str]:
    return {t for t, genes in gs.sets.items() if tag in genes}


EXAMPLE_ROWS: list[dict[str, Any]] = [
    {"locus_tag": "g1", "go_terms": [B2]},
    {"locus_tag": "g2", "go_terms": [B3]},
    {"locus_tag": "g3", "go_terms": [M1, B2]},
    {"locus_tag": "g4"},
]


@pytest.fixture
def res(tmp_path: Path) -> Any:
    return _res(tmp_path)


# --- GO (examples 1-5) --------------------------------------------------------------------------


def test_example1_is_a_propagation_and_parents(res: Any) -> None:
    by, _ = _build(genes_frame(EXAMPLE_ROWS), res)
    bp = by["go_bp"]
    assert _members(bp, "g1") == {B2, B1}
    assert BP_ROOT not in bp.sets
    assert bp.parents[B2] == frozenset({B1})
    assert bp.names[B2] == "toy process two"


def test_example2_part_of_propagates(res: Any) -> None:
    by, _ = _build(genes_frame(EXAMPLE_ROWS), res)
    assert _members(by["go_bp"], "g2") == {B3, B1}  # B1 only through part_of (decision 3)


def test_example3_split_by_namespace(res: Any) -> None:
    by, _ = _build(genes_frame(EXAMPLE_ROWS), res)
    assert _members(by["go_mf"], "g3") == {M1}
    assert _members(by["go_bp"], "g3") == {B2, B1}
    assert _members(by["go_cc"], "g3") == set()
    assert MF_ROOT not in by["go_mf"].sets


def test_example4_coverage_per_namespace(res: Any) -> None:
    by, _ = _build(genes_frame(EXAMPLE_ROWS), res)
    assert by["go_bp"].coverage == pytest.approx(0.75)
    assert by["go_mf"].coverage == pytest.approx(0.25)
    assert by["go_cc"].coverage == pytest.approx(0.0)


def test_example5_absent_id_dropped_and_counted(res: Any) -> None:
    rows = [*EXAMPLE_ROWS, {"locus_tag": "g5", "go_terms": [ABSENT]}]
    by, _ = _build(genes_frame(rows), res)
    for c in GO:
        assert _members(by[c], "g5") == set()
    counts = [by[c].n_unknown_ids for c in GO]
    assert max(counts) == 1 and min(counts) >= 0
    assert by["go_bp"].coverage == pytest.approx(3 / 5)  # g5 eligible, unannotated


def test_obsolete_term_dropped_and_counted(res: Any) -> None:
    rows = [{"locus_tag": "o1", "go_terms": [OBSOLETE]}, {"locus_tag": "o2", "go_terms": [B1]}]
    by, _ = _build(genes_frame(rows), res)
    assert OBSOLETE not in by["go_bp"].sets
    assert _members(by["go_bp"], "o1") == set()
    assert max(by[c].n_unknown_ids for c in GO) == 1
    assert by["go_bp"].coverage == pytest.approx(0.5)


def test_unknown_ids_count_every_dropped_id(res: Any) -> None:
    rows = [{"locus_tag": "u1", "go_terms": [OBSOLETE, ABSENT]}, {"locus_tag": "u2",
             "go_terms": [B1, "GO:8888888"]}]  # fmt: skip
    by, _ = _build(genes_frame(rows), res)
    assert max(by[c].n_unknown_ids for c in GO) == 3
    assert _members(by["go_bp"], "u2") == {B1}


def test_alt_id_maps_to_primary(res: Any) -> None:
    by, _ = _build(genes_frame([{"locus_tag": "a1", "go_terms": [B2_ALT]}]), res)
    bp = by["go_bp"]
    assert _members(bp, "a1") == {B2, B1}
    assert B2_ALT not in bp.sets
    assert all(c.n_unknown_ids == 0 for c in by.values() if c.collection in GO)


def test_diamond_and_mixed_relations(res: Any) -> None:
    rows = [{"locus_tag": "d4", "go_terms": [B4]}, {"locus_tag": "d5", "go_terms": [B5]}]
    by, _ = _build(genes_frame(rows), res)
    bp = by["go_bp"]
    assert _members(bp, "d4") == {B4, B2, B3, B1}
    assert _members(bp, "d5") == {B5, B4, B2, B3, B1}  # part_of then is_a then part_of


def test_cc_part_of(res: Any) -> None:
    by, _ = _build(genes_frame([{"locus_tag": "c", "go_terms": [C2]}]), res)
    assert _members(by["go_cc"], "c") == {C2, C1}
    assert CC_ROOT not in by["go_cc"].sets


def test_parents_direct_within_namespace(res: Any) -> None:
    rows = [{"locus_tag": f"p{i}", "go_terms": [t]} for i, t in enumerate(sorted(TERMS))]
    by, _ = _build(genes_frame(rows), res)
    for c in GO:
        gs = by[c]
        for t in gs.sets:
            assert gs.parents[t] - ROOTS == DIRECT_PARENTS[t] - ROOTS, (c, t)
            assert all(p in gs.sets or p in ROOTS for p in gs.parents[t])
    assert by["go_bp"].parents[B4] == frozenset({B2, B3})
    assert B1 in by["go_bp"].parents[B3]  # the part_of parent is a parent


def test_go_names_cover_sets(res: Any) -> None:
    rows = [{"locus_tag": f"p{i}", "go_terms": [t]} for i, t in enumerate(sorted(TERMS))]
    by, _ = _build(genes_frame(rows), res)
    for c in (*GO, "go_slim"):
        assert set(by[c].sets) <= set(by[c].names)
        for t in by[c].sets:
            assert by[c].names[t] == TERMS[t][0]


# --- GO slim ------------------------------------------------------------------------------------


def test_slim_membership_through_all_ancestors(res: Any) -> None:
    rows = [
        {"locus_tag": "s2", "go_terms": [B2]},
        {"locus_tag": "s3", "go_terms": [B3]},
        {"locus_tag": "s5", "go_terms": [B5]},
        {"locus_tag": "sm", "go_terms": [M2]},
        {"locus_tag": "sc", "go_terms": [C2]},
        {"locus_tag": "sx"},
    ]
    by, _ = _build(genes_frame(rows), res)
    slim = by["go_slim"]
    assert _members(slim, "s2") == {B1}
    assert _members(slim, "s3") == {B3, B1}
    assert _members(slim, "s5") == {B3, B1}  # via part_of B4 then the full DAG
    assert _members(slim, "sm") == {M1}  # one collection across namespaces
    assert _members(slim, "sc") == set()  # C1 is not a slim term; the CC root is removed
    assert not set(slim.sets) & ROOTS
    assert len(slim.parents) == 0  # classic tests only (decision 5)
    assert slim.coverage == pytest.approx(4 / 6)


# --- KEGG (example 6) ---------------------------------------------------------------------------


def test_example6_kegg_dedup_and_overview_excluded(res: Any) -> None:
    rows = [{"locus_tag": "g6", "kegg_ko": ["K00001"]}, {"locus_tag": "g7", "kegg_ko": ["K00002"]}]
    by, _ = _build(genes_frame(rows), res)
    kegg = by["kegg"]
    assert _members(kegg, "g6") == {"map00010"}
    assert _members(kegg, "g7") == {"map00010", "map00020"}
    assert set(kegg.sets) <= set(KEGG_PATHWAYS)
    assert not any(k.startswith(("ko", "path:")) for k in kegg.sets)
    assert not set(kegg.sets) & set(KEGG_OVERVIEW)
    assert kegg.sets["map00010"] == frozenset({"g6", "g7"})


def test_kegg_coverage_counts_pathways_not_kos(res: Any) -> None:
    rows = [
        {"locus_tag": "k1", "kegg_ko": ["K00003"]},
        {"locus_tag": "k2", "kegg_ko": ["K00004"]},  # overview-only KO
        {"locus_tag": "k3", "kegg_ko": ["K99999"]},  # KO with no pathway
        {"locus_tag": "k4"},
    ]
    by, _ = _build(genes_frame(rows), res)
    kegg = by["kegg"]
    assert _members(kegg, "k1") == {"map00020"}
    assert _members(kegg, "k2") == _members(kegg, "k3") == set()
    assert kegg.coverage == pytest.approx(0.25)
    assert len(kegg.parents) == 0


def test_kegg_names_and_missing_name(res: Any) -> None:
    rows = [{"locus_tag": "g7", "kegg_ko": ["K00002"]}]
    by, _ = _build(genes_frame(rows), res)
    kegg = by["kegg"]
    assert kegg.names["map00010"] == "Toy pathway alpha"
    assert "map00020" in kegg.sets  # no name row: kept
    nm = kegg.names.get("map00020", "NA")
    assert nm == "NA" or nm is None or nm is pd.NA or (isinstance(nm, float) and math.isnan(nm))


# --- COG (example 7) ----------------------------------------------------------------------------


def test_example7_cog_letters(res: Any) -> None:
    rows = [
        {"locus_tag": "c1", "cog_id": "COG9001", "cog_category": "KL"},
        {"locus_tag": "c2", "cog_id": "COG9002", "cog_category": "J"},
        {"locus_tag": "c3", "cog_category": "R"},
        {"locus_tag": "c4", "cog_category": "S"},
        {"locus_tag": "c5", "cog_category": "-"},
        {"locus_tag": "c6", "cog_category": ""},
        {"locus_tag": "c7", "cog_category": None},
    ]
    by, _ = _build(genes_frame(rows), res)
    cog = by["cog"]
    assert _members(cog, "c1") == {"K", "L"}
    assert _members(cog, "c2") == {"J"}
    assert _members(cog, "c3") == {"R"} and _members(cog, "c4") == {"S"}  # real categories
    for t in ("c5", "c6", "c7"):
        assert _members(cog, t) == set()
    assert not {"-", "", "NA"} & set(cog.sets)
    assert cog.names["K"] == "Toy transcription" and cog.names["L"] == "Toy replication"
    assert cog.coverage == pytest.approx(4 / 7)
    assert len(cog.parents) == 0


def test_cog_category_pandas_na(res: Any) -> None:
    g = genes_frame([{"locus_tag": "n1", "cog_category": "K"}, {"locus_tag": "n2"}])
    g["cog_category"] = pd.Series(["K", pd.NA], dtype=object)
    by, _ = _build(g, res)
    assert by["cog"].sets == {"K": frozenset({"n1"})}
    assert by["cog"].coverage == pytest.approx(0.5)


# --- eligibility (example 8) --------------------------------------------------------------------


def test_example8_partial_gene_excluded(res: Any) -> None:
    rows = [
        *EXAMPLE_ROWS,
        {"locus_tag": "gp", "partial": True, "go_terms": [B2], "kegg_ko": ["K00001"],
         "cog_category": "K"},
        {"locus_tag": "gn", "p_detect": float("nan"), "go_terms": [M1], "cog_category": "J"},
    ]  # fmt: skip
    by, _ = _build(genes_frame(rows), res)
    for gs in by.values():
        assert _members(gs, "gp") == set() and _members(gs, "gn") == set()
    assert by["go_bp"].coverage == pytest.approx(0.75)  # denominator: g1..g4
    assert by["go_mf"].coverage == pytest.approx(0.25)


def test_eligible_argument_is_used(res: Any) -> None:
    g = genes_frame(EXAMPLE_ROWS)
    el = np.array([True, True, False, True])
    by, _ = _build(g, res, eligible=el)
    assert _members(by["go_bp"], "g3") == set()
    assert by["go_bp"].coverage == pytest.approx(2 / 3)
    assert by["go_mf"].coverage == pytest.approx(0.0)


def test_shared_eligible_helper() -> None:
    from gippa.genetable import eligible

    g = pd.DataFrame({
        "locus_tag": ["a", "b", "c", "d", "e"],
        "partial": pd.Series([False, True, None, False, False], dtype="boolean"),
        "p_detect": [0.3, 0.4, 0.5, float("nan"), float("inf")],
    })  # fmt: skip
    m = eligible(g)
    assert m.dtype == np.bool_
    assert m.tolist() == [True, False, True, False, False]


# --- resources missing (example 9), light DB (example 10) ---------------------------------------


def test_all_installed_gives_six_collections_in_order(res: Any) -> None:
    col, _ = _mods()
    out, flags = col.build_collections(genes_frame(EXAMPLE_ROWS),
                                       eligible_mask(genes_frame(EXAMPLE_ROWS)), res)  # fmt: skip
    assert [g.collection for g in out] == list(COLLECTIONS)
    assert not any(f.startswith("layer_missing") for f in flags)
    assert not _has_flag(list(flags), "bakta_db_light")


def test_example9_go_missing(tmp_path: Path) -> None:
    r = _res(tmp_path, omit=("go_basic",))
    assert r.go is None
    by, flags = _build(genes_frame(EXAMPLE_ROWS), r)
    assert not set(by) & {*GO, "go_slim"}
    assert set(by) == {"cog", "kegg"}
    assert _has_flag(flags, "layer_missing:go")
    assert not _has_flag(flags, "layer_missing:kegg") and not _has_flag(flags, "layer_missing:cog")


def test_go_slim_missing_keeps_go(tmp_path: Path) -> None:
    r = _res(tmp_path, omit=("goslim_generic",))
    assert r.slim is None and r.go is not None
    by, _ = _build(genes_frame(EXAMPLE_ROWS), r)
    assert "go_slim" not in by
    assert set(GO) <= set(by)


def test_kegg_missing(tmp_path: Path) -> None:
    r = _res(tmp_path, omit=("kegg_ko_pathway", "kegg_pathway_names"))
    assert r.kegg_links is None and r.kegg_names is None
    by, flags = _build(genes_frame(EXAMPLE_ROWS), r)
    assert "kegg" not in by and _has_flag(flags, "layer_missing:kegg")
    assert not _has_flag(flags, "layer_missing:go")


def test_cog_missing(tmp_path: Path) -> None:
    r = _res(tmp_path, omit=("cog2020_def", "cog2020_fun"))
    assert r.cog_names is None
    by, flags = _build(genes_frame(EXAMPLE_ROWS), r)
    assert "cog" not in by and _has_flag(flags, "layer_missing:cog")


def test_nothing_installed() -> None:
    _, fr = _mods()
    r = fr.read_resources(dict.fromkeys(
        ("go_basic", "goslim_generic", "kegg_ko_pathway", "kegg_pathway_names", "cog2020_def",
         "cog2020_fun")))  # fmt: skip
    assert (r.go, r.slim, r.kegg_links, r.kegg_names, r.cog_names) == (None,) * 5
    by, flags = _build(genes_frame(EXAMPLE_ROWS), r)
    assert by == {}
    for w in ("layer_missing:go", "layer_missing:kegg", "layer_missing:cog"):
        assert _has_flag(flags, w)


def test_example10_bakta_light_flag(res: Any) -> None:
    rows = [*EXAMPLE_ROWS, {"locus_tag": "g6", "kegg_ko": ["K00001"]}]
    by, flags = _build(genes_frame(rows, annotation_source="bakta_light"), res)
    assert _has_flag(flags, "bakta_db_light")
    assert set(by) == set(COLLECTIONS)
    assert _members(by["kegg"], "g6") == {"map00010"}


# --- degenerate tables --------------------------------------------------------------------------


def test_no_eligible_gene(res: Any) -> None:
    rows = [{**r, "partial": True} for r in EXAMPLE_ROWS]
    by, _ = _build(genes_frame(rows), res)
    assert set(by) == set(COLLECTIONS)
    for gs in by.values():
        assert len(gs.sets) == 0
        assert math.isnan(gs.coverage)


def test_empty_gene_table(res: Any) -> None:
    by, _ = _build(genes_frame([]), res)
    assert set(by) == set(COLLECTIONS)
    assert all(len(gs.sets) == 0 and math.isnan(gs.coverage) for gs in by.values())


def test_no_function_columns(res: Any) -> None:
    g = genes_frame(EXAMPLE_ROWS, function_columns=False)
    by, flags = _build(g, res)
    assert set(by) == set(COLLECTIONS)
    for gs in by.values():
        assert len(gs.sets) == 0
        assert gs.coverage == 0.0
    assert not _has_flag(flags, "bakta_db_light")


def test_no_empty_sets_and_eligible_only(res: Any) -> None:
    rows = [*EXAMPLE_ROWS, {"locus_tag": "gp", "partial": True, "go_terms": [C2]}]
    by, _ = _build(genes_frame(rows), res)
    assert C2 not in by["go_cc"].sets and C1 not in by["go_cc"].sets
    for gs in by.values():
        assert all(len(v) > 0 for v in gs.sets.values())
        assert all(isinstance(v, frozenset) for v in gs.sets.values())


# --- determinism, order, no mutation ------------------------------------------------------------


def _snapshot(by: dict[str, Any]) -> dict[str, Any]:
    return {c: (dict(g.sets), dict(g.names), dict(g.parents),
                None if math.isnan(g.coverage) else g.coverage, g.n_unknown_ids)
            for c, g in by.items()}  # fmt: skip


RICH_ROWS: list[dict[str, Any]] = [
    {"locus_tag": "r1", "go_terms": [B5, M2], "kegg_ko": ["K00002"], "cog_category": "KL"},
    {"locus_tag": "r2", "go_terms": [C2, ABSENT], "kegg_ko": ["K00003"], "cog_category": "J"},
    {"locus_tag": "r3", "go_terms": [B2_ALT], "kegg_ko": ["K00001"], "cog_category": "S"},
    {"locus_tag": "r4", "go_terms": [B3, OBSOLETE], "cog_category": "-"},
    {"locus_tag": "r5", "partial": True, "go_terms": [B1], "cog_category": "M"},
    {"locus_tag": "r6"},
]


def test_deterministic_and_sorted_keys(res: Any) -> None:
    g = genes_frame(RICH_ROWS)
    a, fa = _build(g, res)
    b, fb = _build(g, res)
    assert _snapshot(a) == _snapshot(b) and fa == fb
    for gs in a.values():
        assert list(gs.sets) == sorted(gs.sets)


def test_row_order_invariant(res: Any) -> None:
    g = genes_frame(RICH_ROWS)
    perm = g.iloc[[5, 2, 0, 4, 1, 3]].reset_index(drop=True)
    a, fa = _build(g, res)
    b, fb = _build(perm, res)
    assert _snapshot(a) == _snapshot(b) and sorted(fa) == sorted(fb)


def test_inputs_not_mutated(res: Any) -> None:
    col, _ = _mods()
    g = genes_frame(RICH_ROWS)
    before = g.copy(deep=True)
    lists = copy.deepcopy(g["go_terms"].tolist())
    el = eligible_mask(g)
    el_before = el.copy()
    links = None if res.kegg_links is None else res.kegg_links.copy(deep=True)
    col.build_collections(g, el, res)
    pd.testing.assert_frame_equal(g, before)
    assert g["go_terms"].tolist() == lists
    np.testing.assert_array_equal(el, el_before)
    pd.testing.assert_frame_equal(res.kegg_links, links)


def test_go_terms_as_numpy_arrays(res: Any) -> None:
    """Parquet round-trips list columns as numpy arrays; the result must not change."""
    g = genes_frame(RICH_ROWS)
    a, _ = _build(g, res)
    h = g.copy()
    for c in ("go_terms", "kegg_ko"):
        h[c] = pd.Series([np.array(v, dtype=object) for v in g[c]], dtype=object)
    b, _ = _build(h, res)
    assert _snapshot(a) == _snapshot(b)
