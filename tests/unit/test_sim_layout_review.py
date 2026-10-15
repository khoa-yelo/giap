"""giapsim layout: lead decisions after code review of PR #62 (spec M4-07-sim-layout-tables.md).

1. Circular replicons carry the wrap pair (last gene, first gene) in ``operon_pairs`` (layout and M0
   adapter; M2-12 convention, two-gene circular contigs included); linear contigs never do.
2. ``LayoutConfig(origin_crossing=True)`` rotates ``chr1`` so one gene, one TU and one block region
   straddle the origin (unwrapped ``end > length``, §6.5), all contract-valid; the default False
   leaves the outputs unchanged.
3. One test per surviving mutant of the review (adapter p, singleton TUs, region numbering,
   ``_MIN_AA``, sorted complete ``region_ids``, linear contigs never ``chromosome``).
4. ``genes.amr`` = membership in an AMR region (written); ``ribosomal_protein`` stays unwritten; the
   ``"p"`` fallback of the pair lookup is gone.
5. M0 ``simulate`` with ``real_rrna_trna=True, n_trna=0`` runs.

Checks come from tests/helpers/layout_oracle.py (written from the specs, not from the simulator).
"""

from __future__ import annotations

import csv
import hashlib
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import layout_oracle as lo

BLOCK_TYPES = ("prophage", "bgc", "crispr")


def _layout_mod() -> Any:
    import giapsim.layout

    return giapsim.layout


def _tables_mod() -> Any:
    import giapsim.contract_tables

    return giapsim.contract_tables


def simulate(seed: int = 1, **kw: Any) -> Any:
    m = _layout_mod()
    return m.simulate_layout(m.LayoutConfig(**kw), np.random.default_rng(seed))


def write(layout: Any, work: Path) -> dict[str, pd.DataFrame]:
    _tables_mod().write_contract_tables(layout, work)
    return dict(lo.read_parquets(work))


def _ordered(genes: pd.DataFrame) -> dict[str, list[tuple[str, str]]]:
    """(locus_tag, strand) per contig in genome order."""
    return {
        c: list(zip(g["locus_tag"], g["strand"], strict=True))
        for c, g in lo.genes_by_contig(genes).items()
    }


def _pair_set(pairs: pd.DataFrame) -> set[tuple[str, str]]:
    return set(zip(pairs["gene_a"], pairs["gene_b"], strict=True))


def _wrap_expectation(layout: Any) -> tuple[set[tuple[str, str]], set[tuple[str, str]]]:
    """(wrap pairs that must be listed, wrap pairs that must not be listed)."""
    topo = dict(zip(layout.replicons["contig"], layout.replicons["topology"], strict=True))
    want, forbid = set(), set()
    for c, gs in _ordered(layout.genes).items():
        if len(gs) < 2:
            continue
        (first, s0), (last, s1) = gs[0], gs[-1]
        if topo[c] == "circular" and s0 == s1:
            want.add((last, first))
        elif topo[c] == "linear":
            forbid.add((last, first))
    return want, forbid


# ------------------------------------------------------------------- 1. wrap pair on circular


@pytest.fixture(scope="module")
def mixed() -> Any:
    return simulate(7, n_genes_chromosome=400, plasmids=(25, 2), linear_contigs=(12, 3))


def test_wrap_pair_listed_for_circular_replicons_layout() -> None:
    hits = 0
    for seed in range(12):
        # p_same_strand = 1 makes every contig single-strand, so every circular contig wraps
        for ps in (1.0, 0.6):
            layout = simulate(seed, n_genes_chromosome=120, plasmids=(15, 2),
                              linear_contigs=(9,), regions=(), p_same_strand=ps)  # fmt: skip
            want, forbid = _wrap_expectation(layout)
            got = _pair_set(layout.operon_pairs)
            assert want <= got, sorted(want - got)
            assert not (forbid & got), sorted(forbid & got)
            hits += len(want)
    assert hits > 0


def test_wrap_pair_only_on_circular_same_strand_layout(mixed: Any) -> None:
    """operon_pairs = adjacent same-strand pairs plus the circular wrap pairs, nothing else."""
    want = set(lo.adjacent_same_strand(mixed.genes, mixed.replicons, wrap=False))
    want |= _wrap_expectation(mixed)[0]
    assert _pair_set(mixed.operon_pairs) == want


def test_two_gene_plasmid_single_strand_has_forward_and_wrap_pair() -> None:
    layout = simulate(4, n_genes_chromosome=50, regions=(), plasmids=(2,), p_same_strand=1.0)
    g = lo.genes_by_contig(layout.genes)["p1"]["locus_tag"].tolist()
    got = _pair_set(layout.operon_pairs)
    assert (g[0], g[1]) in got and (g[1], g[0]) in got


def test_wrap_pair_gives_first_and_last_gene_uniop_like_production(tmp_path: Path) -> None:
    layout = simulate(3, n_genes_chromosome=150, regions=(), plasmids=(20,), p_same_strand=1.0,
                      linear_contigs=(10,))  # fmt: skip
    t = write(layout, tmp_path / "w")
    want = lo.expected_uniop(layout.genes, layout.replicons, layout.operon_pairs)
    lay = t["layers"].set_index("locus_tag")
    for c, gs in _ordered(layout.genes).items():
        circ = c != "ctg1"
        for tag in (gs[0][0], gs[-1][0]):
            prev, nxt = want[tag]
            gp, gn = lay.loc[tag, "uniop_p_prev"], lay.loc[tag, "uniop_p_next"]
            assert (prev is None) == lo.is_na(gp) and (nxt is None) == lo.is_na(gn), (c, tag)
            if prev is not None:
                assert gp == pytest.approx(prev)
            if nxt is not None:
                assert gn == pytest.approx(nxt)
            # single-strand circular contig: both neighbours exist, so both p are set
            assert circ == (not lo.is_na(gp) and not lo.is_na(gn)), (c, tag)


@pytest.fixture(scope="module")
def m0_rich(tmp_path_factory: pytest.TempPathFactory) -> tuple[Any, Any]:
    import giapsim

    cfg = giapsim.SimConfig(genome_length=60000, n_plasmids=2, linear_contig_length=20000,
                            origin_crossing_gene=True, n_pseudogenes=5, n_sorfs=5,
                            prophage_bp=10000, bgc_bp=8000)  # fmt: skip
    sim = giapsim.simulate(cfg, 0, tmp_path_factory.mktemp("m0rich"))
    return sim, _layout_mod().layout_from_simulation(sim)


def test_wrap_pair_listed_for_circular_replicons_adapter(m0_rich: Any) -> None:
    _, layout = m0_rich
    want, forbid = _wrap_expectation(layout)
    got = _pair_set(layout.operon_pairs)
    assert want <= got, sorted(want - got)
    assert not (forbid & got)
    full = set(lo.adjacent_same_strand(layout.genes, layout.replicons, wrap=False)) | want
    assert got == full


def test_wrap_pair_listed_adapter_default_m0(tmp_path: Path) -> None:
    import giapsim

    for seed in range(1, 6):
        sim = giapsim.simulate(giapsim.SimConfig(), seed, tmp_path / str(seed))
        layout = _layout_mod().layout_from_simulation(sim)
        want, _ = _wrap_expectation(layout)
        assert want <= _pair_set(layout.operon_pairs), (seed, sorted(want))


# ------------------------------------------------------------------- 2. origin_crossing


def _crossing_layout(seed: int, **kw: Any) -> Any:
    return simulate(seed, n_genes_chromosome=kw.pop("n", 300), origin_crossing=True, **kw)


def _chr1_length(layout: Any) -> int:
    return int(layout.replicons.set_index("contig").loc["chr1", "length"])


def test_origin_crossing_defaults_to_false() -> None:
    assert _layout_mod().LayoutConfig().origin_crossing is False


@pytest.mark.parametrize("seed", range(6))
def test_origin_crossing_gene_tu_and_block_region(seed: int, tmp_path: Path) -> None:
    layout = _crossing_layout(seed, plasmids=(20,), linear_contigs=(10,))
    n = _chr1_length(layout)
    g = layout.genes
    cross = g[g["crosses_origin"].astype(bool)]
    assert len(cross) >= 1
    assert set(cross["contig"]) == {"chr1"}
    assert ((cross["start"] < n) & (cross["end"] > n) & (cross["end"] < 2 * n)).all()
    # TUs and block regions straddling the origin are stored unwrapped (end > length)
    tus = layout.tus[(layout.tus["contig"] == "chr1") & (layout.tus["end"] > n)]
    assert len(tus) >= 1
    assert ((tus["start"] < n) & (tus["end"] < 2 * n)).all()
    reg = layout.regions
    reg = reg[reg["region_type"].isin(BLOCK_TYPES) & (reg["contig"] == "chr1")]
    xreg = reg[reg["end"].astype("float") > n]
    assert len(xreg) >= 1
    assert ((xreg["start"].astype("float") < n) & (xreg["end"].astype("float") < 2 * n)).all()
    # the written tables pass the M4-06 contract
    t = write(layout, tmp_path / "w")
    _, issues = lo.load_work(tmp_path / "w")
    assert lo.errors(issues) == []
    from giap.genetable.contract import validate_table

    for key in lo.WRITTEN_TABLES:
        assert lo.errors(validate_table(key, t[key])) == [], key
    assert t["genes"]["crosses_origin"].sum() == len(cross)


def _covers(gs: int, ge: int, s: int, e: int, n: int) -> bool:
    return (s <= gs and ge <= e) or (s <= gs + n and ge + n <= e)


@pytest.mark.parametrize("seed", range(4))
def test_origin_crossing_tu_and_region_members_consistent(seed: int) -> None:
    layout = _crossing_layout(seed)
    n = _chr1_length(layout)
    order = lo.genes_by_contig(layout.genes)["chr1"]
    pos = {t: i for i, t in enumerate(order["locus_tag"])}
    g = layout.genes.set_index("locus_tag")
    for r in layout.tus[(layout.tus["contig"] == "chr1") & (layout.tus["end"] > n)].itertuples():
        members = lo.as_list(r.genes)
        assert lo.circular_consecutive([pos[m] for m in members], len(pos)), r.tu_id
        assert {g.loc[m, "strand"] for m in members} == {r.strand}
        for m in members:
            assert _covers(int(g.loc[m, "start"]), int(g.loc[m, "end"]), r.start, r.end, n), m
    reg = layout.regions
    for r in reg[reg["region_type"].isin(BLOCK_TYPES)].itertuples():
        members = lo.members(layout.region_genes, r.region_id)
        assert lo.circular_consecutive([pos[m] for m in members], len(pos)), r.region_id
        for m in members:
            assert _covers(int(g.loc[m, "start"]), int(g.loc[m, "end"]), int(r.start),
                           int(r.end), n), (r.region_id, m)  # fmt: skip


def test_origin_crossing_keeps_other_invariants() -> None:
    layout = _crossing_layout(2, plasmids=(15,))
    assert lo.check_gene_coordinates(layout.genes, layout.replicons) == []
    assert lo.check_locus_tags(layout.genes, layout.replicons) == []
    assert lo.check_gaps(layout.genes, layout.tus) == []
    # every gene in exactly one TU; the crossing gene's TU holds it once
    seen = [t for cell in layout.tus["genes"] for t in lo.as_list(cell)]
    assert sorted(seen) == sorted(layout.genes["locus_tag"])
    # the wrap pair is listed when the first and last chr1 genes share a strand
    want, _ = _wrap_expectation(layout)
    assert want <= _pair_set(layout.operon_pairs)


# Digests of LayoutConfig(n_genes_chromosome=300, plasmids=(20,), linear_contigs=(10,)), seed 3,
# taken on the PR #62 head before the review fixes. Tables untouched by decisions 1 and 4 (and the
# gene frame without ``amr``) must stay byte-identical with origin_crossing off.
_DEFAULT_KW: dict[str, Any] = {
    "n_genes_chromosome": 300,
    "plasmids": (20,),
    "linear_contigs": (10,),
}
_PINNED = {
    "replicons.parquet": "344ec075ddc22b024b0f5bf71e148e9457d7d596c816919898af39339c4f3e66",
    "tus.parquet": "52851f936b3d5a5acd0fce8a50e82743c948743a1805e2a6cc40535056fb1dfe",
    "regions.parquet": "cd7687bcc7e1e09baec42eab7f39428fada5b61a27713bfb12fe0ec27407f025",
}
_PINNED_GENES_CSV = "c5f362b1475398f9020bb9e23f2dddd5fba671dc2b112ebfd15a7f3340122e96"


def _genes_digest(genes: pd.DataFrame) -> str:
    g = genes.drop(columns=[c for c in ("amr", "ribosomal_protein") if c in genes.columns])
    return hashlib.sha256(g.to_csv(index=False).encode()).hexdigest()


def _digests(work: Path) -> dict[str, str]:
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(work.iterdir())}


def test_origin_crossing_false_leaves_outputs_unchanged(tmp_path: Path) -> None:
    off = simulate(3, origin_crossing=False, **_DEFAULT_KW)
    _tables_mod().write_contract_tables(off, tmp_path / "off")
    _tables_mod().write_contract_tables(simulate(3, **_DEFAULT_KW), tmp_path / "dflt")
    assert _digests(tmp_path / "off") == _digests(tmp_path / "dflt")
    d = _digests(tmp_path / "off")
    for name, digest in _PINNED.items():
        assert d[name] == digest, name
    assert _genes_digest(off.genes) == _PINNED_GENES_CSV
    assert not off.genes["crosses_origin"].any()


def test_origin_crossing_deterministic(tmp_path: Path) -> None:
    _tables_mod().write_contract_tables(_crossing_layout(5), tmp_path / "a")
    _tables_mod().write_contract_tables(_crossing_layout(5), tmp_path / "b")
    assert _digests(tmp_path / "a") == _digests(tmp_path / "b")


# ------------------------------------------------------------------- 3. surviving mutants


def _truth_tu(sim: Any) -> dict[str, str]:
    with Path(sim.genes_truth).open(newline="") as fh:
        return {r["locus_tag"]: r.get("tu_id") or "" for r in csv.DictReader(fh, delimiter="\t")}


def test_mutant_adapter_pair_p_is_one_within_zero_across_tus(m0_rich: Any) -> None:
    _, layout = m0_rich
    tu = lo.tu_of(layout.tus)
    p = layout.operon_pairs
    same = [tu[a] == tu[b] for a, b in zip(p["gene_a"], p["gene_b"], strict=True)]
    assert p["uniop_p"].tolist() == [1.0 if s else 0.0 for s in same]
    assert any(same) and not all(same)


def test_mutant_adapter_genes_without_truth_tu_are_singletons(m0_rich: Any) -> None:
    sim, layout = m0_rich
    truth = _truth_tu(sim)
    orphans = [g for g in layout.genes["locus_tag"] if not truth.get(g)]
    assert len(orphans) >= 2  # pseudogenes, sORFs, prophage/BGC genes have no truth TU
    size = {str(t): len(lo.as_list(c)) for t, c in zip(layout.tus["tu_id"], layout.tus["genes"],
                                                        strict=True)}  # fmt: skip
    tu = lo.tu_of(layout.tus)
    assert [g for g in orphans if size[tu[g]] != 1] == []


def test_mutant_region_numbering_two_regions_of_same_type() -> None:
    spec = _layout_mod().RegionSpec
    orders = set()
    for seed in range(12):
        layout = simulate(seed, n_genes_chromosome=200,
                          regions=(spec("prophage", 20), spec("prophage", 30)))  # fmt: skip
        r = layout.regions.set_index("region_id")
        assert sorted(r.index) == ["prophage_1", "prophage_2"]
        assert int(r.loc["prophage_1", "start"]) < int(r.loc["prophage_2", "start"]), seed
        orders.add(len(lo.members(layout.region_genes, "prophage_1")))
    assert orders == {20, 30}  # both spec orders were placed first at least once


def test_mutant_min_aa_no_sorf_sized_cds_or_pseudogene() -> None:
    layout = simulate(2, n_genes_chromosome=600, regions=(), gene_len_median_aa=25.0,
                      frac_sorf=0.0, frac_pseudogene=0.2)  # fmt: skip
    n_aa = layout.genes["nt"].str.len() // 3 - 1
    assert (n_aa >= 30).all(), int(n_aa.min())
    assert (n_aa == 30).sum() > 50  # the clip is binding at this median


def test_mutant_layers_region_ids_sorted_and_complete(tmp_path: Path) -> None:
    """A gene in two regions (overlap, as M4-17 will plant) lists both IDs, sorted."""
    layout = simulate(1, n_genes_chromosome=200)
    rg = layout.region_genes
    tag = lo.members(rg, "amr_1")[0]
    extra = pd.DataFrame([("bgc_1", tag, 1.0, None, False)], columns=rg.columns).astype(rg.dtypes)
    # rows in sorted order: a reversal (or truncation) of the per-gene list is caught
    rg2 = pd.concat([extra, rg], ignore_index=True)
    lay = _tables_mod().layers_frame(replace(layout, region_genes=rg2)).set_index("locus_tag")
    assert lo.as_list(lay.loc[tag, "region_ids"]) == ["amr_1", "bgc_1"]
    rg3 = pd.concat([rg, extra], ignore_index=True)
    lay = _tables_mod().layers_frame(replace(layout, region_genes=rg3)).set_index("locus_tag")
    assert lo.as_list(lay.loc[tag, "region_ids"]) == ["amr_1", "bgc_1"]


def test_mutant_truth_region_ids_sorted_and_complete() -> None:
    # ``_truth`` is the layout's truth-frame builder (tooling); overlapping regions arrive in M4-17
    truth_fn = _layout_mod()._truth
    rg = pd.DataFrame({"region_id": ["prophage_1", "amr_1", "bgc_1"], "locus_tag": ["a", "a", "b"]})
    t = truth_fn(["a", "b", "c"], ["TU_1", "TU_1", "TU_2"], rg).set_index("locus_tag")
    assert [lo.as_list(t.loc[x, "region_ids"]) for x in "abc"] == [
        ["amr_1", "prophage_1"],
        ["bgc_1"],
        [],
    ]


def test_mutant_layout_truth_region_ids_match_region_genes(mixed: Any) -> None:
    want: dict[str, set[str]] = {}
    for rid, g in zip(mixed.region_genes["region_id"], mixed.region_genes["locus_tag"],
                      strict=True):  # fmt: skip
        want.setdefault(g, set()).add(rid)
    for tag, ids in zip(mixed.truth["locus_tag"], mixed.truth["region_ids"], strict=True):
        ids = lo.as_list(ids)
        assert ids == sorted(want.get(tag, set())), tag


def test_mutant_linear_contig_never_kind_chromosome(mixed: Any, m0_rich: Any) -> None:
    for layout in (mixed, m0_rich[1]):
        rep = layout.replicons
        lin = rep[rep["topology"] == "linear"]
        assert len(lin) >= 1
        assert (lin["kind"] != "chromosome").all()
    rep = mixed.replicons.set_index("contig")
    assert rep.loc[["ctg1", "ctg2"], "kind"].tolist() == ["unknown", "unknown"]


# ------------------------------------------------------------------- 4. amr column


def _amr_members(layout: Any) -> set[str]:
    reg = layout.regions
    ids = set(reg.loc[reg["region_type"] == "amr", "region_id"])
    rg = layout.region_genes
    return set(rg.loc[rg["region_id"].isin(ids), "locus_tag"])


@pytest.mark.parametrize("kw", [{}, {"plasmids": (30,), "linear_contigs": (10,)},
                                {"regions": ()}])  # fmt: skip
def test_amr_column_is_membership_in_an_amr_region(kw: dict[str, Any], tmp_path: Path) -> None:
    layout = simulate(5, n_genes_chromosome=300, **kw)
    assert "amr" in layout.genes.columns
    want = layout.genes["locus_tag"].isin(_amr_members(layout))
    assert layout.genes["amr"].astype(bool).tolist() == want.tolist()
    t = write(layout, tmp_path / "w")
    assert "amr" in t["genes"].columns
    assert t["genes"]["amr"].astype(bool).tolist() == want.tolist()
    assert "ribosomal_protein" not in t["genes"].columns
    assert "ribosomal_protein" not in layout.genes.columns


def test_amr_column_counts_scattered_amr_genes() -> None:
    spec = _layout_mod().RegionSpec
    layout = simulate(
        6, n_genes_chromosome=200, regions=(spec("amr", 7, True), spec("amr", 5, True))
    )
    assert int(layout.genes["amr"].sum()) == 12


def test_amr_column_adapter_all_false(m0_rich: Any, tmp_path: Path) -> None:
    _, layout = m0_rich
    assert "amr" in layout.genes.columns and not layout.genes["amr"].astype(bool).any()
    t = write(layout, tmp_path / "w")
    assert "amr" in t["genes"].columns and "ribosomal_protein" not in t["genes"].columns


def test_pair_lookup_has_no_p_fallback() -> None:
    layout = simulate(1, n_genes_chromosome=100, regions=())
    renamed = layout.operon_pairs.rename(columns={"uniop_p": "p"})
    with pytest.raises(KeyError):
        _tables_mod().layers_frame(replace(layout, operon_pairs=renamed))


# ------------------------------------------------------------------- 5. M0 real tRNA with n_trna=0


def test_m0_real_rrna_trna_with_zero_trna_runs(tmp_path: Path) -> None:
    import giapsim

    cfg = giapsim.SimConfig(real_rrna_trna=True, n_trna=0, n_rrna_operons=1)
    sim = giapsim.simulate(cfg, 1, tmp_path / "sim")
    with Path(sim.genes_truth).open(newline="") as fh:
        assert list(csv.DictReader(fh, delimiter="\t"))
    layout = _layout_mod().layout_from_simulation(sim)
    assert lo.errors(lo.load_work(_write_dir(layout, tmp_path))[1]) == []


def test_m0_real_rrna_trna_with_zero_trna_and_no_rrna_runs(tmp_path: Path) -> None:
    import giapsim

    giapsim.simulate(giapsim.SimConfig(real_rrna_trna=True, n_trna=0), 2, tmp_path / "sim")


def _write_dir(layout: Any, tmp_path: Path) -> Path:
    work = tmp_path / "w"
    _tables_mod().write_contract_tables(layout, work)
    return work
