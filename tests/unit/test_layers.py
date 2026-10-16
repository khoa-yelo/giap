"""M2-20 ``giap.annotate.layers``: ``LAYER_COLUMNS``, ``gene_layers``, ``layer_flags``
(docs/specs/M2-20-context-layers.md examples 2-6, edge cases; GIAP_SPEC.md §9.3, §9.4, §6.3).

Inputs are the merged producers' outputs on giapsim SIM100K (``helpers.layers_chain``): every
``gene_layers`` argument is the DataFrame read with ``pd.read_parquet`` from the producer's file
(``genes``, ``tus``, ``pairs_raw``, ``promoter_windows``, ``tx_rates``, ``terminators``, ``rbs``,
``codon``, ``region_genes``, ``replicons``), ``None`` for a missing optional layer.

Pinned where the spec is silent (test author):

* One row per gene of ``genes.parquet``, in that order; columns exactly ``LAYER_COLUMNS``.
* ``region_ids`` and ``flags`` are native list cells (C6 contract, ``giap.genetable.contract``);
  ``flags`` is exactly the set {partial, crosses_origin, pseudogene, cai_reference} that applies.
* Promoter / tx / terminator columns equal the merged per-gene functions (M2-13
  ``gene_promoter_features`` with ``promotech_yes``; M2-14 ``tu_tx_features`` of the gene's TU;
  M2-15 ``terminator_conf_downstream``, checked where the next same-strand gene is a pair
  partner); rbs / codon columns equal the ``rbs`` / ``codon`` rows of the gene (NA without one).
* ``layer_flags(cfg, present)``: ``present`` is the set of ``context_layers`` input keys given;
  result sorted, unique, every item ``layer_missing:<layer>``; ``genomad`` / ``antismash`` follow
  ``cfg.modules`` (``mobile`` / ``bgc``), ``tx_rates`` follows the key (and
  ``promoter_calculator.enabled``).
"""

from __future__ import annotations

import itertools
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd
import pytest
from helpers import layers_oracle as orc
from helpers.layers_chain import (
    NO_RBS_EVERY,
    PARTIAL_LOCUS,
    Chain,
    build_chain,
    call_gene_layers,
    fresh_steps,
    precomputed_operon_pairs,
    read_layers_frames,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.config import GiapConfig

SPEC_COLUMNS = (
    "locus_tag", "tu_id", "tu_position", "tu_size", "uniop_p_prev", "uniop_p_next", "region_ids",
    "promoter_score", "promoter_score_tu", "n_promoter_windows_tu",
    "tx_rate_max_tu", "tx_rate_sum_tu", "tss_distance_bp", "terminator_conf_downstream",
    "tir_ostir", "sd_dg", "mrna_dg", "spacing_dg", "standby_dg", "spacer_bp", "start_codon",
    "cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "rare_codons_2_10", "flags",
)  # fmt: skip
RBS_COLS = ("tir_ostir", "sd_dg", "mrna_dg", "spacing_dg", "standby_dg", "spacer_bp", "start_codon")
CODON_COLS = ("cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "rare_codons_2_10")
TX_COLS = ("tx_rate_max_tu", "tx_rate_sum_tu", "tss_distance_bp")
ALL_KEYS = {
    "genes", "features", "replicons", "pairs_raw", "tus", "promoter_windows", "tx_rates",
    "terminators", "terminators_crosscheck", "rbs", "codon", "regions", "region_genes",
}  # fmt: skip


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(tmp_path_factory.mktemp("layers_chain"))


@pytest.fixture(scope="module")
def layers(chain: Chain) -> pd.DataFrame:
    return call_gene_layers(read_layers_frames(chain))


def _row(df: pd.DataFrame, tag: str) -> dict[str, Any]:
    rows = df[df["locus_tag"] == tag]
    assert len(rows) == 1, f"{tag}: {len(rows)} rows"
    return rows.iloc[0].to_dict()


def _list(cell: Any) -> list[str]:
    assert cell is not None and not isinstance(cell, str), f"expected a list cell, got {cell!r}"
    return [str(x) for x in cell]


def _gene_runs(chain: Chain, strand: str) -> list[tuple[str, str, str]]:
    """(A, B, C) consecutive same-strand genes in genome order linked by pairs (A,B), (B,C)."""
    p = chain.frame("pairs_raw")
    g = chain.frame("genes").set_index("locus_tag")
    nxt = dict(zip(p["gene_a"], p["gene_b"], strict=True))
    return [
        (a, b, nxt[b]) for a, b in nxt.items()
        if b in nxt and g.loc[a, "strand"] == strand and g.loc[a, "contig"] == "chr1"
        and g.loc[a, "start"] < g.loc[b, "start"] < g.loc[nxt[b], "start"]
    ]  # fmt: skip


# ----------------------------------------------------------------------------- columns / shape


def test_layer_columns_are_the_spec_tuple() -> None:
    from giap.annotate.layers import LAYER_COLUMNS

    assert tuple(LAYER_COLUMNS) == SPEC_COLUMNS


def test_layer_columns_equal_the_c6_contract_columns() -> None:
    """M4-06 C6 consumes layers.parquet with exactly these columns."""
    from giap.annotate.layers import LAYER_COLUMNS
    from giap.genetable.contract import CONTRACT_COLUMNS

    assert tuple(LAYER_COLUMNS) == CONTRACT_COLUMNS["layers"]


def test_one_row_per_gene_in_gene_table_order(chain: Chain, layers: pd.DataFrame) -> None:
    assert tuple(layers.columns) == SPEC_COLUMNS
    assert layers["locus_tag"].tolist() == chain.frame("genes")["locus_tag"].tolist()


def test_layers_frame_passes_the_c6_contract(chain: Chain, layers: pd.DataFrame) -> None:
    from giap.genetable.contract import validate_inputs, validate_table

    errors = [i for i in validate_table("layers", layers) if i.severity == "error"]
    assert not errors, errors
    empty = pd.DataFrame()
    tables = SimpleNamespace(
        replicons=chain.frame("replicons"), genes=chain.frame("genes"), tus=chain.frame("tus"),
        layers=layers, regions=chain.frame("regions"), region_genes=chain.frame("region_genes"),
        peptides=None, peptide_hits=empty, protein_hits=empty,
    )  # fmt: skip
    bad = [i for i in validate_inputs(tables) if i.code.startswith("contract.layers")]
    assert not bad, bad


# ----------------------------------------------------------------------------- TU context


def test_tu_columns_follow_tus_parquet(chain: Chain, layers: pd.DataFrame) -> None:
    want = orc.tu_positions(chain.frame("tus"))
    for r in layers.itertuples(index=False):
        assert (r.tu_id, int(r.tu_position), int(r.tu_size)) == want[r.locus_tag], r.locus_tag


def test_every_gene_has_a_tu(layers: pd.DataFrame) -> None:
    """Edge case: genes without a TU cannot exist (singletons are TUs)."""
    assert layers["tu_id"].notna().all()
    assert (layers["tu_size"] >= 1).all() and (layers["tu_position"] >= 1).all()


def test_uniop_p_prev_next_in_transcription_order(chain: Chain, layers: pd.DataFrame) -> None:
    want = orc.prev_next(chain.frame("genes"), chain.frame("pairs_raw"))
    for r in layers.itertuples(index=False):
        wp, wn = want[r.locus_tag]
        assert orc.same(r.uniop_p_prev, wp), (r.locus_tag, r.uniop_p_prev, wp)
        assert orc.same(r.uniop_p_next, wn), (r.locus_tag, r.uniop_p_next, wn)


def test_example_2_plus_tu_middle_gene(chain: Chain, tmp_path: Path) -> None:
    """TU (A, B, C) on +, p(A,B)=0.9, p(B,C)=0.7 -> B: position 2, size 3, prev 0.9, next 0.7."""
    runs = _gene_runs(chain, "+")
    assert runs, "SIM100K has three linked + genes"
    a, b, c = runs[0]
    p = chain.frame("pairs_raw")
    over = {(a, b): 0.9, (b, c): 0.7}
    over |= {(x, a): 0.1 for x in p.loc[p["gene_b"] == a, "gene_a"]}
    over |= {(c, y): 0.1 for y in p.loc[p["gene_a"] == c, "gene_b"]}
    pairs, tus = precomputed_operon_pairs(chain, tmp_path, over)
    tu = [list(g) for g in pd.read_parquet(tus)["genes"] if a in list(g)]
    assert tu == [[a, b, c]], "the precomputed pairs make TU (A, B, C)"
    df = call_gene_layers(read_layers_frames(chain, pairs_raw=pairs, tus=tus))
    rb = _row(df, b)
    assert (int(rb["tu_position"]), int(rb["tu_size"])) == (2, 3)
    assert rb["uniop_p_prev"] == pytest.approx(0.9) and rb["uniop_p_next"] == pytest.approx(0.7)
    ra = _row(df, a)
    assert int(ra["tu_position"]) == 1 and ra["uniop_p_next"] == pytest.approx(0.9)
    if (p["gene_b"] == a).any():  # lead Q2: the first gene keeps p to its upstream neighbour
        assert ra["uniop_p_prev"] == pytest.approx(0.1)


def test_example_3_minus_tu_last_gene(chain: Chain, tmp_path: Path) -> None:
    """- TU (C, B, A) (genome order A, B, C), the gene before A in genome order is on + ->
    A: position 3, uniop_p_prev = p(A,B), uniop_p_next = NA."""
    g = chain.frame("genes")
    order = g["locus_tag"].tolist()
    strand = dict(zip(g["locus_tag"], g["strand"], strict=True))
    runs = [r for r in _gene_runs(chain, "-")
            if order.index(r[0]) > 0 and strand[order[order.index(r[0]) - 1]] == "+"]  # fmt: skip
    assert runs, "SIM100K has a - run whose genome predecessor is on +"
    a, b, c = runs[0]
    p = chain.frame("pairs_raw")
    over = {(a, b): 0.8, (b, c): 0.85}
    over |= {(c, y): 0.1 for y in p.loc[p["gene_a"] == c, "gene_b"]}
    pairs, tus = precomputed_operon_pairs(chain, tmp_path, over)
    tu = [list(x) for x in pd.read_parquet(tus)["genes"] if a in list(x)]
    assert tu == [[c, b, a]], "transcription order of a - TU"
    df = call_gene_layers(read_layers_frames(chain, pairs_raw=pairs, tus=tus))
    ra = _row(df, a)
    assert (int(ra["tu_position"]), int(ra["tu_size"])) == (3, 3)
    assert ra["uniop_p_prev"] == pytest.approx(0.8)
    assert orc.isna(ra["uniop_p_next"])
    rc = _row(df, c)
    assert int(rc["tu_position"]) == 1 and rc["uniop_p_next"] == pytest.approx(0.85)
    if (p["gene_a"] == c).any():  # upstream same-strand neighbour below the threshold (Q2)
        assert rc["uniop_p_prev"] == pytest.approx(0.1)


@settings(
    max_examples=12, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(data=st.data())
def test_property_tu_positions_and_shared_pair_probabilities(
    chain: Chain, tmp_path_factory: pytest.TempPathFactory, data: st.DataObject
) -> None:
    """For any pair probabilities (wrap pairs unlinked): positions are 1..size in each TU;
    consecutive genes of a TU share their pair p (next of k == prev of k+1, >= threshold); a
    TU's first gene has prev NA or < threshold and its last gene next NA or < threshold."""
    p = chain.frame("pairs_raw")
    probs = data.draw(st.lists(st.sampled_from([0.05, 0.3, 0.49, 0.5, 0.51, 0.8, 0.99]),
                               min_size=len(p), max_size=len(p)))  # fmt: skip
    start = dict(zip(chain.frame("genes")["locus_tag"], chain.frame("genes")["start"], strict=True))
    # wrap pairs (last, first) stay unlinked: a fully linked ring breaks at the origin (M2-12)
    over = {(a, b): (x if start[a] < start[b] else 0.05)
            for a, b, x in zip(p["gene_a"], p["gene_b"], probs, strict=True)}  # fmt: skip
    pairs, tus = precomputed_operon_pairs(chain, tmp_path_factory.mktemp("prop"), over)
    df = call_gene_layers(read_layers_frames(chain, pairs_raw=pairs, tus=tus)).set_index(
        "locus_tag"
    )
    thr = GiapConfig().annotation.uniop.pair_threshold
    for gs in pd.read_parquet(tus)["genes"]:
        gs = list(gs)
        assert [int(df.loc[x, "tu_position"]) for x in gs] == list(range(1, len(gs) + 1))
        assert {int(df.loc[x, "tu_size"]) for x in gs} == {len(gs)}
        for x, y in itertools.pairwise(gs):
            assert df.loc[x, "uniop_p_next"] == df.loc[y, "uniop_p_prev"] >= thr
        first, last = df.loc[gs[0], "uniop_p_prev"], df.loc[gs[-1], "uniop_p_next"]
        assert orc.isna(first) or first < thr
        assert orc.isna(last) or last < thr


# ----------------------------------------------------------------------------- regions


def test_region_ids_sorted_from_region_genes(chain: Chain, layers: pd.DataFrame) -> None:
    want = orc.region_ids(chain.frame("region_genes"))
    for r in layers.itertuples(index=False):
        assert _list(r.region_ids) == want.get(r.locus_tag, []), r.locus_tag


def test_example_4_gene_in_prophage_and_plasmid(chain: Chain, layers: pd.DataFrame) -> None:
    rg = chain.frame("region_genes")
    both = set(rg.loc[rg["region_id"] == "prophage_1", "locus_tag"]) & set(
        rg.loc[rg["region_id"] == "plasmid_1", "locus_tag"]
    )
    assert both, "the chain plants a provirus on plasmid_1"
    for tag in sorted(both):
        assert _list(_row(layers, tag)["region_ids"]) == ["plasmid_1", "prophage_1"]


def test_example_5_core_only_has_no_prophage_or_bgc_ids(chain: Chain, tmp_path: Path) -> None:
    """modules=[core]: region_layer runs without geNOMAD / antiSMASH (annotate.nf), so no
    prophage / bgc ids; the run flags name both missing tools."""
    from giap.annotate.layers import layer_flags

    outs = {"regions": tmp_path / "regions.parquet", "region_genes": tmp_path / "rg.parquet"}
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        reg.run_step(
            "region_layer", GiapConfig(seed=1, modules=["core"]),
            {k: chain.inputs[k] for k in ("genes", "features", "replicons")}, outs,
        )  # fmt: skip
    df = call_gene_layers(read_layers_frames(chain, region_genes=outs["region_genes"]))
    ids = {x for cell in df["region_ids"] for x in _list(cell)}
    assert not {x for x in ids if x.startswith(("prophage_", "bgc_"))}
    assert "plasmid_1" in ids  # the Bakta replicon table still gives the plasmid region
    flags = layer_flags(GiapConfig(seed=1, modules=["core"]), ALL_KEYS)
    assert {"layer_missing:genomad", "layer_missing:antismash"} <= set(flags)


# ----------------------------------------------------------------------- tx, promoters, terminators


def test_example_6_no_tx_rates_gives_na_and_flag(chain: Chain) -> None:
    from giap.annotate.layers import layer_flags

    df = call_gene_layers(read_layers_frames(chain, tx_rates=None))
    for c in TX_COLS:
        assert df[c].isna().all(), c
    cfg = GiapConfig.model_validate({"annotation": {"promoter_calculator": {"enabled": False}}})
    assert "layer_missing:tx_rates" in layer_flags(cfg, ALL_KEYS - {"tx_rates"})


def test_tx_columns_are_the_tu_features(chain: Chain, layers: pd.DataFrame) -> None:
    want = orc.tx_features(read_layers_frames(chain))
    assert (want["tx_rate_max"] > 0).any(), "some TUs have Promoter Calculator hits"
    for r in layers.itertuples(index=False):
        w = want.loc[r.tu_id]
        assert orc.same(r.tx_rate_max_tu, float(w["tx_rate_max"])), r.locus_tag
        assert orc.same(r.tx_rate_sum_tu, float(w["tx_rate_sum"])), r.locus_tag
        wd = w["tss_distance_bp"]
        assert orc.same(r.tss_distance_bp, None if pd.isna(wd) else int(wd)), r.locus_tag


def test_promoter_columns_are_the_gene_promoter_features(
    chain: Chain, layers: pd.DataFrame
) -> None:
    yes = GiapConfig().annotation.promoters.promotech_yes
    want = orc.promoter_features(read_layers_frames(chain), yes)
    for r in layers.itertuples(index=False):
        w = want.loc[r.locus_tag]
        for c in ("promoter_score", "promoter_score_tu", "n_promoter_windows_tu"):
            assert orc.same(getattr(r, c), w[c]), (r.locus_tag, c)


def test_no_promoter_layer_gives_na(chain: Chain) -> None:
    df = call_gene_layers(read_layers_frames(chain, promoter_windows=None))
    for c in ("promoter_score", "promoter_score_tu", "n_promoter_windows_tu"):
        assert df[c].isna().all(), c


def test_terminator_conf_downstream(chain: Chain, layers: pd.DataFrame) -> None:
    want = orc.terminator_downstream_with_partner(read_layers_frames(chain))
    assert any(not math.isnan(v) for v in want.values()), "recorded TerminatorNet calls hit"
    got = dict(zip(layers["locus_tag"], layers["terminator_conf_downstream"], strict=True))
    for tag, v in want.items():
        assert orc.same(got[tag], v), tag
    conf = set(chain.frame("terminators")["confidence"].astype(float))
    assert all(orc.isna(v) or float(v) in conf for v in got.values())


def test_no_terminator_layer_gives_na(chain: Chain) -> None:
    df = call_gene_layers(read_layers_frames(chain, terminators=None))
    assert df["terminator_conf_downstream"].isna().all()


# ----------------------------------------------------------------------------- rbs / codon / flags


def test_rbs_columns_and_missing_rbs_rows_are_na(chain: Chain, layers: pd.DataFrame) -> None:
    rbs = chain.frame("rbs").set_index("locus_tag")
    n_missing = 0
    for r in layers.itertuples(index=False):
        for c in RBS_COLS:
            w = rbs.loc[r.locus_tag, c] if r.locus_tag in rbs.index else None
            assert orc.same(getattr(r, c), w), (r.locus_tag, c)
        n_missing += r.locus_tag not in rbs.index
    assert n_missing >= len(layers) // NO_RBS_EVERY  # edge case exercised


def test_codon_columns_and_gene_missing_from_codon_is_na(chain: Chain) -> None:
    codon = chain.frame("codon")
    dropped = set(codon["locus_tag"].iloc[::5])
    frames = read_layers_frames(chain)
    frames["codon"] = codon[~codon["locus_tag"].isin(dropped)].reset_index(drop=True)
    df = call_gene_layers(frames)
    full = codon.set_index("locus_tag")
    for r in df.itertuples(index=False):
        for c in CODON_COLS:
            w = None if r.locus_tag in dropped else full.loc[r.locus_tag, c]
            assert orc.same(getattr(r, c), w), (r.locus_tag, c)


def test_no_rbs_or_codon_layer_gives_na(chain: Chain) -> None:
    df = call_gene_layers(read_layers_frames(chain, rbs=None, codon=None))
    for c in (*RBS_COLS, *CODON_COLS):
        assert df[c].isna().all(), c


def test_gene_flags(chain: Chain, layers: pd.DataFrame) -> None:
    want = orc.gene_flags(chain.frame("genes"), chain.frame("codon"))
    assert any("partial" in f for f in want.values())
    assert any("pseudogene" in f for f in want.values())
    assert any("cai_reference" in f for f in want.values())
    for r in layers.itertuples(index=False):
        assert set(_list(r.flags)) == want[r.locus_tag], r.locus_tag


def test_partial_gene_keeps_its_layer_values(chain: Chain, layers: pd.DataFrame) -> None:
    """Edge case: partial genes keep layer values (statistics exclude them later, §6.5)."""
    r = _row(layers, PARTIAL_LOCUS)
    assert "partial" in _list(r["flags"])
    assert not orc.isna(r["promoter_score"])
    assert not orc.isna(r["tu_id"]) and not orc.isna(r["at_2_10"])


# ----------------------------------------------------------------------------- layer_flags


def test_layer_flags_all_present_all_modules_is_empty() -> None:
    from giap.annotate.layers import layer_flags

    assert tuple(layer_flags(GiapConfig(), ALL_KEYS)) == ()


@pytest.mark.parametrize(
    ("modules", "present", "want"),
    [
        (["core"], ALL_KEYS, {"layer_missing:genomad", "layer_missing:antismash"}),
        (["core", "mobile"], ALL_KEYS, {"layer_missing:antismash"}),
        (["core", "bgc"], ALL_KEYS, {"layer_missing:genomad"}),
        (["core", "mobile", "bgc"], ALL_KEYS - {"tx_rates"}, {"layer_missing:tx_rates"}),
    ],
)
def test_layer_flags(modules: list[str], present: set[str], want: set[str]) -> None:
    from giap.annotate.layers import layer_flags

    got = tuple(layer_flags(GiapConfig(modules=modules), set(present)))
    assert want <= set(got)
    assert list(got) == sorted(set(got))
    assert all(f.startswith("layer_missing:") for f in got)
    if modules == ["core", "mobile"]:
        assert "layer_missing:genomad" not in got
    if modules == ["core", "bgc"]:
        assert "layer_missing:antismash" not in got
