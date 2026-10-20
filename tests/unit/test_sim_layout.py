"""giapsim sequence-free layout and M4-06 contract tables (issue M4-07; GIAP_SPEC.md §11.3, §6.5).

Examples 1-9, the edge cases and acceptance criteria of docs/specs/M4-07-sim-layout-tables.md:
``simulate_layout`` → ``write_contract_tables`` passes ``giap.io.work_tables.load_m4_inputs`` with
no error (and every written table passes ``validate_table`` with no issue at all),
``CONTRACT_COLUMNS`` equals ``giap.genetable.contract.TABLES``, the layout agrees with itself
(TUs, pairs, layers, regions) and with the M0 simulator truth (adapter), runs are deterministic,
the M0 outputs are unchanged, and the new modules import nothing from ``giap``.

Checks come from tests/helpers/layout_oracle.py (written from the spec, not from the simulator).
This file lives in tests/unit (like test_sim_basic.py) because tests/sim is implementer tooling.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
import time
from dataclasses import replace
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import layout_oracle as lo

REPO = Path(__file__).resolve().parents[2]
GIAPSIM = REPO / "tests" / "sim" / "giapsim"

# M0 digests (tests/unit/test_sim_features.py, main at 4c72c2a): unchanged by M4-07.
M0_SHA256 = {
    "genome.fasta": "698373a45c787caa64c01d7050dd1d2a7b31e03771e4e389056ca510c6df1971",
    "bakta_like.json": "d312420eb3d404cb1abbb6d07844cccbc83c0ce5e9c912d3a02aa17da6c49064",
    "truth/genes.tsv": "9c2c98152eed218d3855ccddf478bf49ea42bc69727fce466111f704e0387b0f",
}
NUMERIC_LAYER_COLUMNS = (
    "promoter_score", "promoter_score_tu", "n_promoter_windows_tu", "tx_rate_max_tu",
    "tx_rate_sum_tu", "tss_distance_bp", "terminator_conf_downstream", "tir_ostir", "sd_dg",
    "mrna_dg", "spacing_dg", "standby_dg", "spacer_bp", "cai", "tai", "enc", "fop", "gc3", "cpb",
    "at_2_10", "rare_codons_2_10",
)  # fmt: skip


def _layout_mod() -> ModuleType:
    import giapsim.layout

    return giapsim.layout


def _tables_mod() -> ModuleType:
    import giapsim.contract_tables

    return giapsim.contract_tables


def simulate(seed: int = 1, **kw: Any) -> Any:
    m = _layout_mod()
    return m.simulate_layout(m.LayoutConfig(**kw), np.random.default_rng(seed))


def write(layout: Any, work: Path) -> dict[str, pd.DataFrame]:
    _tables_mod().write_contract_tables(layout, work)
    return dict(lo.read_parquets(work))


@pytest.fixture(scope="module")
def default(tmp_path_factory: pytest.TempPathFactory) -> tuple[Any, dict[str, pd.DataFrame], Path]:
    """Example 1: ``simulate_layout(LayoutConfig(), default_rng(1))`` and its written tables."""
    layout = simulate(1)
    work = tmp_path_factory.mktemp("default") / "work"
    return layout, write(layout, work), work


@pytest.fixture(scope="module")
def mixed(tmp_path_factory: pytest.TempPathFactory) -> tuple[Any, dict[str, pd.DataFrame], Path]:
    """Chromosome + two plasmids + two linear contigs (every replicon kind at once)."""
    layout = simulate(7, n_genes_chromosome=400, plasmids=(25, 1), linear_contigs=(12, 1))
    work = tmp_path_factory.mktemp("mixed") / "work"
    return layout, write(layout, work), work


# ---------------------------------------------------------------------------- examples 1 and 2


def test_example1_default_layout(default: Any) -> None:
    layout, t, _ = default
    rep, genes, tus = t["replicons"], t["genes"], t["tus"]
    assert rep[["contig", "topology", "kind"]].values.tolist() == [
        ["chr1", "circular", "chromosome"]
    ]
    assert len(genes) == 1000
    assert set(genes["contig"]) == {"chr1"}
    counts: dict[str, int] = {}
    for cell in tus["genes"]:
        for g in lo.as_list(cell):
            counts[g] = counts.get(g, 0) + 1
    assert set(counts) == set(genes["locus_tag"]) and set(counts.values()) == {1}
    sizes = [len(lo.as_list(c)) for c in tus["genes"]]
    assert abs(np.mean(sizes) - 1.8) <= 0.2, np.mean(sizes)
    # the in-memory layout carries the same genes and TUs
    assert sorted(layout.genes["locus_tag"]) == sorted(genes["locus_tag"])
    assert sorted(map(str, layout.tus["tu_id"])) == sorted(map(str, tus["tu_id"]))


def test_example2_contract_tables_load_without_errors(default: Any) -> None:
    _, _, work = default
    inputs, issues = lo.load_work(work)
    assert lo.errors(issues) == []
    assert len(inputs.genes) == 1000
    assert inputs.annotation_source


@pytest.mark.parametrize("key", lo.WRITTEN_TABLES)
def test_example2_each_written_table_has_no_issue_at_all(default: Any, key: str) -> None:
    """No error and no warning (e.g. extra columns) for any written table."""
    from giap.genetable.contract import validate_table

    _, t, _ = default
    assert validate_table(key, t[key]) == []


def test_example2_cross_table_invariants_hold_in_memory(default: Any) -> None:
    """``validate_inputs`` on the loaded frames (genes in TUs, layers, region subsets, bounds)."""
    from giap.genetable.contract import validate_inputs

    _, _, work = default
    inputs, _ = lo.load_work(work)
    assert validate_inputs(inputs) == []


def test_writes_exactly_the_m2_side_files_and_creates_work_dir(tmp_path: Path) -> None:
    work = tmp_path / "a" / "work"  # does not exist yet
    _tables_mod().write_contract_tables(simulate(2, n_genes_chromosome=200), work)
    assert sorted(p.name for p in work.iterdir()) == sorted(lo.M2_FILES)
    meta = json.loads((work / "annotation_meta.json").read_text(encoding="utf-8"))
    assert isinstance(meta["annotation_source"], str) and meta["annotation_source"]


# ---------------------------------------------------------------------------- drift check


def test_contract_columns_equal_tables_for_every_written_table() -> None:
    from giap.genetable.contract import TABLES

    cc = _tables_mod().CONTRACT_COLUMNS
    assert set(lo.WRITTEN_TABLES) <= set(cc)
    assert set(cc) <= set(TABLES)
    want = {k: tuple(c.name for c in v) for k, v in TABLES.items()}
    for key, cols in cc.items():
        assert tuple(cols) == want[key], key


def test_contract_columns_is_a_literal_copy() -> None:
    """The simulator keeps its own literal copy (independence): not the giap object itself."""
    from giap.genetable import contract

    assert _tables_mod().CONTRACT_COLUMNS is not contract.CONTRACT_COLUMNS


def test_written_columns_follow_contract_order(default: Any) -> None:
    from giap.genetable.contract import TABLES

    _, t, _ = default
    for key in lo.WRITTEN_TABLES:
        required = [c.name for c in TABLES[key] if not c.optional]
        allowed = [c.name for c in TABLES[key]]
        cols = list(t[key].columns)
        assert [c for c in cols if c in required] == required, key
        assert [c for c in allowed if c in cols] == cols, key


# ---------------------------------------------------------------------------- examples 3 and 4


def test_example3_plasmid_replicon_and_region(tmp_path: Path) -> None:
    t = write(simulate(1, plasmids=(60,)), tmp_path / "w")
    rep = t["replicons"].set_index("contig")
    assert list(t["replicons"]["contig"]) == ["chr1", "p1"]
    assert rep.loc["p1", "topology"] == "circular" and rep.loc["p1", "kind"] == "plasmid"
    p1 = t["genes"][t["genes"]["contig"] == "p1"]
    assert len(p1) == 60
    reg = t["regions"].set_index("region_id")
    assert reg.loc["plasmid_1", "region_type"] == "plasmid"
    assert reg.loc["plasmid_1", "contig"] == "p1"
    assert (
        reg.loc["plasmid_1", "start"] == 0
        and reg.loc["plasmid_1", "end"] == rep.loc["p1", "length"]
    )
    assert sorted(lo.members(t["region_genes"], "plasmid_1")) == sorted(p1["locus_tag"])
    assert (t["genes"]["contig"] == "chr1").sum() == 1000


def test_example4_linear_contigs_partial_ends(tmp_path: Path) -> None:
    t = write(simulate(1, linear_contigs=(30, 30)), tmp_path / "w")
    rep = t["replicons"].set_index("contig")
    assert list(t["replicons"]["contig"]) == ["chr1", "ctg1", "ctg2"]
    assert list(rep.loc[["ctg1", "ctg2"], "topology"]) == ["linear", "linear"]
    g = t["genes"]
    assert int(g["partial"].sum()) == 4
    for c, gg in lo.genes_by_contig(g).items():
        want = [False] * len(gg) if c == "chr1" else [True] + [False] * (len(gg) - 2) + [True]
        assert gg["partial"].astype(bool).tolist() == want, c
        if c != "chr1":
            assert len(gg) == 30


def test_mixed_replicons_contig_order_and_counts(mixed: Any) -> None:
    _, t, work = mixed
    assert list(t["replicons"]["contig"]) == ["chr1", "p1", "p2", "ctg1", "ctg2"]
    n = t["genes"].groupby("contig").size().to_dict()
    assert n == {"chr1": 400, "p1": 25, "p2": 1, "ctg1": 12, "ctg2": 1}
    kinds = dict(zip(t["replicons"]["contig"], t["replicons"]["kind"], strict=True))
    assert kinds["chr1"] == "chromosome" and kinds["p1"] == kinds["p2"] == "plasmid"
    # a one-gene linear contig: that gene is both the first and the last (partial)
    one = t["genes"][t["genes"]["contig"] == "ctg2"]
    assert one["partial"].astype(bool).tolist() == [True]
    assert lo.errors(lo.load_work(work)[1]) == []
    assert {"plasmid_1", "plasmid_2"} <= set(t["regions"]["region_id"])


# ---------------------------------------------------------------------------- example 5: regions


def test_example5_default_region_ids(default: Any) -> None:
    _, t, _ = default
    assert sorted(t["regions"]["region_id"]) == sorted(["prophage_1", "bgc_1", "crispr_1", "amr_1"])
    types = dict(zip(t["regions"]["region_id"], t["regions"]["region_type"], strict=True))
    assert types == {"prophage_1": "prophage", "bgc_1": "bgc", "crispr_1": "crispr", "amr_1": "amr"}


@pytest.mark.parametrize(("rid", "n"), [("prophage_1", 40), ("bgc_1", 25), ("crispr_1", 6)])
def test_example5_blocks_are_consecutive_genes(default: Any, rid: str, n: int) -> None:
    _, t, _ = default
    chrom = lo.genes_by_contig(t["genes"])["chr1"]
    index = {g: i for i, g in enumerate(chrom["locus_tag"])}
    mem = lo.members(t["region_genes"], rid)
    assert len(mem) == n == len(set(mem))
    assert lo.circular_consecutive([index[g] for g in mem], len(chrom)), rid


def test_example5_blocks_do_not_overlap(default: Any) -> None:
    _, t, _ = default
    sets = [set(lo.members(t["region_genes"], r)) for r in ("prophage_1", "bgc_1", "crispr_1")]
    assert not (sets[0] & sets[1] or sets[0] & sets[2] or sets[1] & sets[2])


def test_example5_prophage_subsets_and_core_biosynthetic(default: Any) -> None:
    _, t, _ = default
    rg = t["region_genes"]
    pro = rg[rg["region_id"] == "prophage_1"]
    assert pro["subset"].notna().all()
    assert set(pro["subset"]) <= set(lo.SUBSETS)
    assert rg.loc[rg["region_id"] != "prophage_1", "subset"].isna().all()
    assert not rg.loc[rg["region_id"] != "bgc_1", "core_biosynthetic"].astype(bool).any()
    assert rg.loc[rg["region_id"] == "bgc_1", "core_biosynthetic"].astype(bool).any()


def test_example5_amr_scattered_with_na_coordinates(default: Any) -> None:
    _, t, _ = default
    amr = t["regions"].set_index("region_id").loc["amr_1"]
    assert amr["contig"] == "chr1"
    assert pd.isna(amr["start"]) and pd.isna(amr["end"])
    mem = lo.members(t["region_genes"], "amr_1")
    assert len(mem) == len(set(mem)) == 10
    chrom = lo.genes_by_contig(t["genes"])["chr1"]
    index = {g: i for i, g in enumerate(chrom["locus_tag"])}
    assert not lo.circular_consecutive([index[g] for g in mem], len(chrom))  # scattered


def test_example5_region_coordinates_cover_members(default: Any) -> None:
    _, t, _ = default
    g = t["genes"].set_index("locus_tag")
    rep = t["replicons"].set_index("contig")
    for r in t["regions"].itertuples(index=False):
        if r.region_type == "amr":
            continue
        length = int(rep.loc[r.contig, "length"])
        assert 0 <= r.start < r.end, r.region_id
        for row in t["region_genes"][t["region_genes"]["region_id"] == r.region_id].itertuples():
            gg = g.loc[row.locus_tag]
            assert gg["contig"] == r.contig
            frac = lo.overlap_fraction(
                int(gg["start"]), int(gg["end"]), int(r.start), int(r.end), length,
                rep.loc[r.contig, "topology"] == "circular",
            )  # fmt: skip
            assert frac >= 0.5, (r.region_id, row.locus_tag, frac)
            assert 0.5 <= row.overlap_frac <= 1.0
            assert row.overlap_frac == pytest.approx(frac, abs=1e-9)


def test_example5_region_ids_numbered_by_contig_order_then_start(mixed: Any) -> None:
    _, t, _ = mixed
    order = lo.contig_order(t["replicons"])
    for typ, sub in t["regions"].groupby("region_type"):
        sub = sub.assign(_o=sub["contig"].map(order), _s=sub["start"].fillna(-1))
        ids = sub.sort_values(["_o", "_s"], kind="mergesort")["region_id"].tolist()
        assert ids == [f"{typ}_{i + 1}" for i in range(len(ids))], typ


# ---------------------------------------------------------------------------- structure (spec text)


def test_gene_coordinates_sequences_and_lengths(default: Any, mixed: Any) -> None:
    for _, t, _ in (default, mixed):
        assert lo.check_gene_coordinates(t["genes"], t["replicons"]) == []


def test_gaps_within_and_between_tus(default: Any, mixed: Any) -> None:
    for _, t, _ in (default, mixed):
        assert lo.check_gaps(t["genes"], t["tus"]) == []


def test_tus_partition_order_and_ids(default: Any, mixed: Any) -> None:
    for _, t, _ in (default, mixed):
        assert lo.check_tus(t["genes"], t["tus"], t["replicons"]) == []


def test_locus_tags_in_contig_order(default: Any, mixed: Any) -> None:
    for _, t, _ in (default, mixed):
        assert lo.check_locus_tags(t["genes"], t["replicons"]) == []


def test_gene_frame_sorted_by_contig_then_start(mixed: Any) -> None:
    _, t, _ = mixed
    g = t["genes"]
    order = lo.contig_order(t["replicons"])
    keys = list(zip(g["contig"].map(order), g["start"], strict=True))
    assert keys == sorted(keys)


def test_feature_types_present_by_default(default: Any) -> None:
    _, t, _ = default
    g = t["genes"]
    ft = g["feature_type"].value_counts().to_dict()
    assert ft.get("sorf", 0) > 0 and ft.get("pseudogene", 0) > 0
    assert g["hypothetical"].astype(bool).any() and not g["hypothetical"].astype(bool).all()


def test_operon_pairs_every_adjacent_same_strand_pair(default: Any, mixed: Any) -> None:
    for layout, t, _ in (default, mixed):
        pairs = layout.operon_pairs
        got = {frozenset(p) for p in zip(pairs["gene_a"], pairs["gene_b"], strict=True)}
        need = {frozenset(p) for p in lo.adjacent_same_strand(t["genes"], t["replicons"])}
        allowed = {
            frozenset(p) for p in lo.adjacent_same_strand(t["genes"], t["replicons"], wrap=True)
        }
        assert need <= got
        assert got <= allowed, sorted(map(sorted, got - allowed))[:3]
        assert len(pairs) == len(got)  # no duplicate pair
        p = pairs[lo.pair_prob_column(pairs)].to_numpy(dtype=float)
        assert ((p >= 0) & (p <= 1)).all()


def test_operon_pairs_intra_tu_higher_than_inter(default: Any) -> None:
    layout, t, _ = default
    tmap = lo.tu_of(t["tus"])
    pairs = layout.operon_pairs
    col = lo.pair_prob_column(pairs)
    intra = [p for a, b, p in zip(pairs["gene_a"], pairs["gene_b"], pairs[col], strict=True)
             if tmap[a] == tmap[b]]  # fmt: skip
    inter = [p for a, b, p in zip(pairs["gene_a"], pairs["gene_b"], pairs[col], strict=True)
             if tmap[a] != tmap[b]]  # fmt: skip
    assert len(intra) > 100 and len(inter) > 50
    assert np.mean(intra) > 0.7 > 0.3 > np.mean(inter)


def test_tus_are_the_planted_ones_not_rederived(tmp_path: Path) -> None:
    """``tus`` = planted TUs: with Beta(8, 2) both ways, pairs say "linked" but TUs stay planted."""
    t = write(
        simulate(3, n_genes_chromosome=300, regions=(), uniop_beta_out=(8, 2)), tmp_path / "w"
    )
    sizes = [len(lo.as_list(c)) for c in t["tus"]["genes"]]
    assert np.mean(sizes) < 3.0  # re-deriving from p would merge most same-strand runs
    assert lo.check_gaps(t["genes"], t["tus"]) == []


def test_layers_tu_columns(default: Any, mixed: Any) -> None:
    for _, t, _ in (default, mixed):
        lay = t["layers"].set_index("locus_tag")
        assert sorted(lay.index) == sorted(t["genes"]["locus_tag"])
        for r in t["tus"].itertuples(index=False):
            members = lo.as_list(r.genes)
            for k, gtag in enumerate(members):
                row = lay.loc[gtag]
                assert row["tu_id"] == r.tu_id
                assert int(row["tu_position"]) == k + 1
                assert int(row["tu_size"]) == len(members)


def test_layers_uniop_prev_next_in_transcription_order(default: Any, mixed: Any) -> None:
    for layout, t, _ in (default, mixed):
        want = lo.expected_uniop(t["genes"], t["replicons"], layout.operon_pairs)
        lay = t["layers"]
        for tag, prev, nxt in zip(
            lay["locus_tag"], lay["uniop_p_prev"], lay["uniop_p_next"], strict=True
        ):
            for got, exp, what in ((prev, want[tag][0], "prev"), (nxt, want[tag][1], "next")):
                if exp is None:
                    assert pd.isna(got), (tag, what, got)
                else:
                    assert got == pytest.approx(exp, abs=1e-12), (tag, what)


def test_layers_region_ids_and_flags(default: Any, mixed: Any) -> None:
    for _, t, _ in (default, mixed):
        reg_of: dict[str, list[str]] = {}
        for rid, tag in zip(
            t["region_genes"]["region_id"], t["region_genes"]["locus_tag"], strict=True
        ):
            reg_of.setdefault(tag, []).append(rid)
        g = t["genes"].set_index("locus_tag")
        for tag, ids, flags in zip(
            t["layers"]["locus_tag"], t["layers"]["region_ids"], t["layers"]["flags"], strict=True
        ):
            assert lo.as_list(ids) == sorted(reg_of.get(tag, [])), tag
            fl = set(lo.as_list(flags))
            assert ("partial" in fl) == bool(g.loc[tag, "partial"]), tag
            assert ("crosses_origin" in fl) == bool(g.loc[tag, "crosses_origin"]), tag
            assert ("pseudogene" in fl) == (g.loc[tag, "feature_type"] == "pseudogene"), tag
            assert fl <= {"partial", "crosses_origin", "pseudogene"}, fl


@pytest.mark.parametrize("col", NUMERIC_LAYER_COLUMNS)
def test_layers_other_columns_all_na(default: Any, col: str) -> None:
    _, t, _ = default
    assert t["layers"][col].isna().all(), col


def test_xrefs_formats_and_genes(default: Any) -> None:
    _, t, _ = default
    x = t["xrefs"]
    assert set(x["locus_tag"]) <= set(t["genes"]["locus_tag"])
    assert set(x["db"]) <= set(lo.XREF_RE)
    assert {"COG", "GO", "KEGG"} <= set(x["db"])
    for db, i in zip(x["db"], x["id"], strict=True):
        assert lo.XREF_RE[db].match(i), (db, i)
    frac = x["locus_tag"].nunique() / len(t["genes"])
    assert 0.6 <= frac <= 0.8


def test_truth_frame_has_tu_ids(default: Any) -> None:
    layout, t, _ = default
    truth = layout.truth
    assert isinstance(truth, pd.DataFrame) and "tu_id" in truth.columns
    assert set(map(str, truth["tu_id"])) <= set(map(str, t["tus"]["tu_id"]))
    if "locus_tag" in truth.columns:
        tmap = lo.tu_of(t["tus"])
        assert all(tmap[g] == tu for g, tu in zip(truth["locus_tag"], truth["tu_id"], strict=True))


def test_written_tables_equal_layout_frames(default: Any) -> None:
    layout, t, _ = default
    for key in ("replicons", "genes", "tus", "regions", "region_genes"):
        mem: pd.DataFrame = getattr(layout, key)
        disk = t[key]
        common = [c for c in disk.columns if c in mem.columns and disk[c].map(np.isscalar).all()]
        assert {"contig", "locus_tag", "tu_id", "region_id"} & set(common) or key == "tus"
        a = (
            mem[common]
            .reset_index(drop=True)
            .astype(object)
            .where(mem[common].notna().values, None)
        )
        b = (
            disk[common]
            .reset_index(drop=True)
            .astype(object)
            .where(disk[common].notna().values, None)
        )
        assert a.values.tolist() == b.values.tolist(), key


# ---------------------------------------------------------------------------- knobs


def test_p_same_strand_one_gives_single_strand_contigs(tmp_path: Path) -> None:
    t = write(simulate(4, n_genes_chromosome=200, regions=(), p_same_strand=1.0, plasmids=(20,)),
              tmp_path / "w")  # fmt: skip
    for c, gg in lo.genes_by_contig(t["genes"]).items():
        assert gg["strand"].nunique() == 1, c


def test_mean_tu_size_one_gives_singletons(tmp_path: Path) -> None:
    layout = simulate(5, n_genes_chromosome=200, regions=(), mean_tu_size=1.0)
    t = write(layout, tmp_path / "w")
    assert {len(lo.as_list(c)) for c in t["tus"]["genes"]} == {1}
    assert lo.check_gaps(t["genes"], t["tus"]) == []


def test_zero_fractions(tmp_path: Path) -> None:
    t = write(simulate(6, n_genes_chromosome=300, regions=(), frac_sorf=0.0, frac_pseudogene=0.0,
                       frac_hypothetical=0.0), tmp_path / "w")  # fmt: skip
    assert set(t["genes"]["feature_type"]) == {"cds"}
    assert not t["genes"]["hypothetical"].astype(bool).any()
    t = write(
        simulate(6, n_genes_chromosome=300, regions=(), frac_hypothetical=1.0), tmp_path / "x"
    )
    assert t["genes"]["hypothetical"].astype(bool).all()


def test_no_regions_gives_empty_region_tables(tmp_path: Path) -> None:
    work = tmp_path / "w"
    t = write(simulate(1, n_genes_chromosome=100, regions=()), work)
    assert len(t["regions"]) == 0 and len(t["region_genes"]) == 0
    assert all(lo.as_list(v) == [] for v in t["layers"]["region_ids"])
    assert lo.errors(lo.load_work(work)[1]) == []


def test_example9_short_genes_in_prophage(tmp_path: Path) -> None:
    t = write(simulate(1, short_genes_in="prophage"), tmp_path / "w")
    g = t["genes"].set_index("locus_tag")
    pro = set(lo.members(t["region_genes"], "prophage_1"))
    span = (g["end"] - g["start"]).astype(float)
    inside = span[span.index.isin(pro)].median()
    outside = span[~span.index.isin(pro)].median()
    assert inside <= 0.5 * outside, (inside, outside)


def test_example9_default_prophage_not_short(default: Any) -> None:
    _, t, _ = default
    g = t["genes"].set_index("locus_tag")
    pro = set(lo.members(t["region_genes"], "prophage_1"))
    span = (g["end"] - g["start"]).astype(float)
    assert span[span.index.isin(pro)].median() > 0.5 * span[~span.index.isin(pro)].median()


# ---------------------------------------------------------------------------- edge cases


def test_regions_that_cannot_fit_raise() -> None:
    with pytest.raises(ValueError):
        simulate(1, n_genes_chromosome=50)  # 40 + 25 + 6 consecutive genes do not fit


@pytest.mark.parametrize("n", [0, 1, -3])
def test_fewer_than_two_chromosome_genes_raise(n: int) -> None:
    with pytest.raises(ValueError):
        simulate(1, n_genes_chromosome=n, regions=())


def test_two_gene_chromosome_is_valid(tmp_path: Path) -> None:
    work = tmp_path / "w"
    t = write(simulate(1, n_genes_chromosome=2, regions=()), work)
    assert len(t["genes"]) == 2
    assert lo.errors(lo.load_work(work)[1]) == []


def test_one_gene_plasmid_one_tu_no_pairs(tmp_path: Path) -> None:
    layout = simulate(2, n_genes_chromosome=200, regions=(), plasmids=(1,))
    t = write(layout, tmp_path / "w")
    p1 = t["genes"][t["genes"]["contig"] == "p1"]["locus_tag"].tolist()
    assert len(p1) == 1
    tus = t["tus"][t["tus"]["contig"] == "p1"]
    assert len(tus) == 1 and lo.as_list(tus["genes"].iloc[0]) == p1
    pairs = layout.operon_pairs
    assert not (pairs["gene_a"].isin(p1) | pairs["gene_b"].isin(p1)).any()
    lay = t["layers"].set_index("locus_tag").loc[p1[0]]
    assert pd.isna(lay["uniop_p_prev"]) and pd.isna(lay["uniop_p_next"])
    assert lo.members(t["region_genes"], "plasmid_1") == p1


def test_locus_tags_unique_across_replicons(mixed: Any) -> None:
    _, t, _ = mixed
    tags = t["genes"]["locus_tag"]
    assert tags.is_unique
    assert all(lo.LOCUS_RE.match(x) for x in tags)


# ------------------------------------------------------------------------ example 6: determinism


def _digests(work: Path) -> dict[str, str]:
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(work.iterdir())}


def test_example6_same_seed_byte_identical(tmp_path: Path) -> None:
    kw: dict[str, Any] = dict(n_genes_chromosome=300, plasmids=(20,), linear_contigs=(10,))
    np.random.seed(0)  # noqa: NPY002
    _tables_mod().write_contract_tables(simulate(11, **kw), tmp_path / "a")
    np.random.seed(12345)  # noqa: NPY002 (the global numpy state must not matter)
    _ = np.random.random(100)  # noqa: NPY002
    _tables_mod().write_contract_tables(simulate(11, **kw), tmp_path / "b")
    da, db = _digests(tmp_path / "a"), _digests(tmp_path / "b")
    assert sorted(da) == sorted(lo.M2_FILES)
    assert da == db


def test_example6_different_seed_differs(tmp_path: Path) -> None:
    _tables_mod().write_contract_tables(simulate(11, n_genes_chromosome=300), tmp_path / "a")
    _tables_mod().write_contract_tables(simulate(12, n_genes_chromosome=300), tmp_path / "b")
    da, db = _digests(tmp_path / "a"), _digests(tmp_path / "b")
    assert da["genes.parquet"] != db["genes.parquet"]
    assert da["tus.parquet"] != db["tus.parquet"]


def test_layout_config_is_frozen_and_hashable() -> None:
    m = _layout_mod()
    cfg = m.LayoutConfig()
    hash(cfg)
    with pytest.raises(AttributeError):
        cfg.n_genes_chromosome = 5  # type: ignore[misc]
    assert replace(cfg, n_genes_chromosome=10).n_genes_chromosome == 10
    assert cfg.n_genes_chromosome == 1000 and cfg.mean_tu_size == 1.8
    assert cfg.p_same_strand == 0.6 and cfg.gene_len_median_aa == 280.0
    assert tuple(cfg.uniop_beta_in) == (8, 2) and tuple(cfg.uniop_beta_out) == (2, 8)
    assert cfg.short_genes_in is None and cfg.plasmids == () and cfg.linear_contigs == ()


# ---------------------------------------------------------------------------- example 7: speed


def _best_of(n: int, fn: Any) -> float:
    best = float("inf")
    for _ in range(n):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def test_example7_speed_4500_genes() -> None:
    assert _best_of(3, lambda: simulate(1, n_genes_chromosome=4500)) < 2.0


# The 1,000-gene target (best of 3 < 0.2 s, M4-07 acceptance) is checked in
# benchmarks/test_sim_layout_bench.py: wall clock is flaky on loaded nodes.


# ---------------------------------------------------------------------------- example 8: M0 adapter


def _truth_rows(path: Path) -> list[dict[str, str]]:
    import csv

    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def m0(tmp_path_factory: pytest.TempPathFactory) -> tuple[Any, Any]:
    import giapsim

    sim = giapsim.simulate(giapsim.SimConfig(), 1, tmp_path_factory.mktemp("m0"))
    return sim, _layout_mod().layout_from_simulation(sim)


def test_example8_adapter_genes_match_truth(m0: Any) -> None:
    sim, layout = m0
    truth = _truth_rows(sim.genes_truth)
    g = layout.genes.set_index("locus_tag")
    assert sorted(g.index) == sorted(r["locus_tag"] for r in truth)
    for r in truth:
        row = g.loc[r["locus_tag"]]
        assert int(row["start"]) == int(r["start"]) - 1  # 0-based
        assert int(row["end"]) == int(r["end"])
        assert row["strand"] == r["strand"]
        assert row["contig"] == "chr1"
        if "start_codon" in g.columns:
            assert row["start_codon"] == r["start_codon"]


def test_example8_adapter_tus_match_truth(m0: Any) -> None:
    sim, layout = m0
    want: dict[str, set[str]] = {}
    for r in _truth_rows(sim.genes_truth):
        want.setdefault(r["tu_id"], set()).add(r["locus_tag"])
    got = {
        str(t): set(lo.as_list(c))
        for t, c in zip(layout.tus["tu_id"], layout.tus["genes"], strict=True)
    }
    assert got == want


def test_example8_adapter_sequences_and_replicons_from_json(m0: Any) -> None:
    sim, layout = m0
    data = json.loads(Path(sim.bakta_json).read_text())
    cds = {f["locus"]: f for f in data["features"] if f["type"] == "cds"}
    g = layout.genes.set_index("locus_tag")
    for tag, f in cds.items():
        assert g.loc[tag, "nt"] == f["nt"]
        assert g.loc[tag, "aa"] == f["aa"].rstrip("*")
    rep = layout.replicons
    assert rep["contig"].tolist() == [s["id"] for s in data["sequences"]]
    assert rep["length"].tolist() == [s["length"] for s in data["sequences"]]
    assert rep["topology"].tolist() == ["circular"]


def test_example8_adapter_has_no_regions_and_valid_tables(m0: Any, tmp_path: Path) -> None:
    _, layout = m0
    assert len(layout.regions) == 0 and len(layout.region_genes) == 0
    work = tmp_path / "w"
    t = write(layout, work)
    assert lo.errors(lo.load_work(work)[1]) == []
    assert lo.check_tus(t["genes"], t["tus"], t["replicons"]) == []
    from giap.genetable.contract import validate_table

    for key in lo.WRITTEN_TABLES:
        assert lo.errors(validate_table(key, t[key])) == [], key


def test_adapter_feature_rich_m0_output_is_contract_valid(tmp_path: Path) -> None:
    """Pseudogenes, sORFs and a plasmid (M3-18 features) still give a valid gene partition."""
    import giapsim

    cfg = giapsim.SimConfig(n_plasmids=1, n_pseudogenes=2, n_sorfs=2, plasmid_length=(8000, 8000))
    sim = giapsim.simulate(cfg, 3, tmp_path / "sim")
    layout = _layout_mod().layout_from_simulation(sim)
    work = tmp_path / "w"
    t = write(layout, work)
    assert lo.errors(lo.load_work(work)[1]) == []
    assert {"sorf", "pseudogene"} <= set(t["genes"]["feature_type"])
    assert set(t["replicons"]["kind"]) == {"chromosome", "plasmid"}


def test_adapter_does_not_modify_simulator_output(tmp_path: Path) -> None:
    import giapsim

    sim = giapsim.simulate(giapsim.SimConfig(), 1, tmp_path / "sim")
    before = {
        p: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted((tmp_path / "sim").rglob("*"))
        if p.is_file()
    }
    _layout_mod().layout_from_simulation(sim)
    after = {
        p: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted((tmp_path / "sim").rglob("*"))
        if p.is_file()
    }
    assert before == after


# ------------------------------------------------------------------------ backward compatibility


def test_m0_outputs_unchanged_with_layout_modules_loaded(tmp_path: Path) -> None:
    import giapsim

    _layout_mod(), _tables_mod()
    giapsim.simulate(giapsim.SimConfig(), seed=1, outdir=tmp_path)
    for name, digest in M0_SHA256.items():
        assert hashlib.sha256((tmp_path / name).read_bytes()).hexdigest() == digest, name


def test_existing_public_api_unchanged() -> None:
    import inspect

    import giapsim

    _layout_mod()
    assert list(inspect.signature(giapsim.simulate).parameters) == ["config", "seed", "outdir"]
    for name in ("SIM100K", "SIM100K_SEED", "SimConfig", "SimResult", "StrainConfig", "simulate"):
        assert hasattr(giapsim, name), name


# ---------------------------------------------------------------------------- independence


_GIAP_IMPORT = re.compile(
    r"^\s*(?:from\s+giap(?:\.\w+)*\s+import|import\s+giap(?:\.|\s|,|$))", re.M
)
# Modules that predate M4-07 and already use giap's writers (M0, M2-05); every other giapsim
# module, in particular the M4-07 ones, must import nothing from ``giap``.
_GRANDFATHERED = frozenset({"simulate.py", "layers.py"})


def test_new_giapsim_modules_import_nothing_from_giap() -> None:
    _layout_mod(), _tables_mod()
    files = sorted(p for p in GIAPSIM.rglob("*.py") if p.name not in _GRANDFATHERED)
    assert {"layout.py", "contract_tables.py"} <= {p.name for p in files}
    bad = [str(p.relative_to(REPO)) for p in files if _GIAP_IMPORT.search(p.read_text())]
    assert bad == []


def test_import_regex_sanity() -> None:
    assert _GIAP_IMPORT.search("from giap.genetable.contract import TABLES")
    assert _GIAP_IMPORT.search("import giap.io")
    assert _GIAP_IMPORT.search("    from giap import rng")
    assert not _GIAP_IMPORT.search("from giapsim.codons import STOPS")
    assert not _GIAP_IMPORT.search("import giapsim")
