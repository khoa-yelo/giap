"""M4-15 ``genes.tsv`` assembly: examples 2-9 and the edge cases of
docs/specs/M4-15-gene-table-assembly.md (GIAP_SPEC.md §6.5, §7.3, §9.2-§9.4).

Inputs are the hand-made M4-06 contract frames of ``m4_contract_fixture`` run through the merged
M4-10 ... M4-14 producers (tests/helpers/gene_table_fixture.py), so ``assemble_genes`` sees exactly
what its producers return.

Pinned by the test author (spec gaps, see the hand-back report):
- ``assemble_genes`` returns the internal representation: 0-based half-open ``start``/``end``,
  Python lists (or arrays) in list columns, NaN/None for missing; columns exactly
  ``GENES_COLUMNS``. The TSV conversion is only in ``giap.io.genes_tsv.write_genes_tsv``.
- ``function_columns(xrefs, cog_def)`` takes the COG table as a frame with columns ``cog_id`` and
  ``cog_category`` (extra columns ignored); list cells (``go_terms``, ``kegg_ko``, ``ec``) are
  sorted and de-duplicated; xrefs of other databases are ignored.
- Lead decisions 4, 5, 9: ``gc_local_dev`` is copied from C6 when present; ``assemble_genes(...,
  cog_def=None)`` (keyword) sets ``cog_category``; several COG IDs -> the smallest and its
  category, bare category xrefs (``COG:L``) ignored; no id map -> ``source_id = locus_tag``;
  ``island_id`` / ``replichore`` string and ``replichore_position`` float, all NA in M4.
- ``detectability_constant`` goes on eligible (non-partial) genes only (M4-13 spec "M4-15 sets
  it"); ``source_id`` comes from C14 when the id map is given.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import gene_table_fixture as gt
from gene_table_fixture import as_list, is_na, row, tag

ASSEMBLE = Path(__file__).resolve().parents[3] / "src" / "giap" / "genetable" / "assemble.py"


@pytest.fixture(scope="module")
def world() -> gt.World:
    return gt.build_world()


@pytest.fixture(scope="module")
def table(world: gt.World) -> pd.DataFrame:
    return world.assemble()


def _tsv(tmp_path: Path, df: pd.DataFrame) -> pd.DataFrame:
    """Write with ``write_genes_tsv`` and read back the raw cells (strings, ``NA`` kept)."""
    from giap.io.genes_tsv import write_genes_tsv

    path = tmp_path / "genes.tsv"
    write_genes_tsv(path, df)
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


# ----------------------------------------------------------------------------- vocabularies


def test_genes_columns_follow_section_9_3() -> None:
    from giap.genetable.assemble import GENES_COLUMNS

    assert tuple(GENES_COLUMNS) == gt.GENES_COLUMNS
    assert "aa" not in GENES_COLUMNS and "nt" not in GENES_COLUMNS


def test_gene_flags_vocabulary_order() -> None:
    from giap.genetable.assemble import GENE_FLAGS

    assert tuple(GENE_FLAGS) == gt.GENE_FLAGS


def test_no_column_says_expressed() -> None:
    """§1.2 rule 6: say "detected", never "expressed"."""
    from giap.genetable.assemble import GENE_FLAGS, GENES_COLUMNS

    assert not [c for c in (*GENES_COLUMNS, *GENE_FLAGS) if "express" in c.lower()]


def test_gene_flags_union_in_vocabulary_order() -> None:
    from giap.genetable.assemble import gene_flags

    parts = {
        "layers": ["pseudogene"],
        "evidence": ["strain_difference", "tolerant_evidence", "multi_identical"],
        "abundance": ["imputed_abundance", "multi_identical"],
    }
    assert gene_flags(parts) == [
        "pseudogene", "multi_identical", "imputed_abundance", "tolerant_evidence",
        "strain_difference",
    ]  # fmt: skip
    assert gene_flags({}) == []
    assert gene_flags({"a": []}) == []


@pytest.mark.parametrize("bad", ["expressed", "same_db", "layer_missing:cog", ""])
def test_gene_flags_vocabulary_is_closed(bad: str) -> None:
    from giap.genetable.assemble import gene_flags

    with pytest.raises(ValueError, match="flag"):
        gene_flags({"layers": ["partial", bad]})


# ----------------------------------------------------------------------------- whole table


def test_one_row_per_gene_columns_and_order(world: gt.World, table: pd.DataFrame) -> None:
    assert tuple(table.columns) == gt.GENES_COLUMNS
    assert len(table) == len(world.inputs.genes) == 20
    assert table["locus_tag"].tolist() == [tag(i) for i in range(20)]  # already in §6.5 order


def test_identity_columns(world: gt.World, table: pd.DataFrame) -> None:
    g = world.inputs.genes.set_index("locus_tag")
    cov = world.cov.set_index("locus_tag")
    for t in table["locus_tag"]:
        r = row(table, t)
        for c in ("feature_type", "contig", "start", "end", "strand"):
            assert r[c] == g.loc[t, c], (t, c)
        for c in ("hypothetical", "partial", "crosses_origin"):
            assert bool(r[c]) is bool(g.loc[t, c]), (t, c)
        for c in ("gene", "product"):
            assert (is_na(r[c]) and is_na(g.loc[t, c])) or r[c] == g.loc[t, c], (t, c)
        assert int(r["length_aa"]) == int(cov.loc[t, "length_aa"])
        assert r["replicon_type"] == "chromosome"
        assert r["annotation_source"] == "bakta_full"
    assert row(table, tag(3))["source_id"] == "orig_3"  # C14 id map
    assert row(table, tag(7))["feature_type"] == "pseudogene"


def test_protein_evidence_abundance_columns_come_from_their_producers(
    world: gt.World, table: pd.DataFrame
) -> None:
    sources: list[tuple[pd.DataFrame, tuple[str, ...]]] = [
        (world.cov, ("mw_da", "pi", "gravy", "tm_segments", "n_observable_peptides")),
        (world.detection.table, ("detected", "detection_tier", "detection_status",
                                 "gene_q_value")),
        (world.evidence, ("detected_2pep", "n_peptides_unique", "n_peptides_multi_identical",
                          "n_peptides_shared", "n_peptides_tolerant", "coverage_frac",
                          "best_peptide_q")),
        (world.abundance.table, gt.ABUNDANCE),
    ]  # fmt: skip
    for src, cols in sources:
        s = src.set_index(src["locus_tag"].astype(str))
        for t in table["locus_tag"]:
            r = row(table, t)
            for c in cols:
                want, got = s.loc[t, c], r[c]
                if is_na(want):
                    assert is_na(got), (t, c, got)
                elif isinstance(want, float):
                    assert math.isclose(float(got), want, rel_tol=1e-12), (t, c)
                else:
                    assert got == want, (t, c, got, want)
    p = world.detectability.p_detect
    for t in table["locus_tag"]:
        got = row(table, t)["p_detect"]
        assert (is_na(got) and is_na(p[t])) or math.isclose(float(got), float(p[t])), t


def test_context_copied_from_layers(world: gt.World, table: pd.DataFrame) -> None:
    """Example 4 (and: no operon neighbour -> ``uniop_p_prev/next`` NA)."""
    lay = world.inputs.layers.set_index("locus_tag")
    for t in table["locus_tag"]:
        r = row(table, t)
        assert r["tu_id"] == lay.loc[t, "tu_id"]
        assert int(r["tu_position"]) == int(lay.loc[t, "tu_position"])
        assert int(r["tu_size"]) == int(lay.loc[t, "tu_size"])
        for c in ("uniop_p_prev", "uniop_p_next"):
            want = lay.loc[t, c]
            assert (is_na(r[c]) and is_na(want)) or float(r[c]) == float(want), (t, c)
        assert as_list(r["region_ids"]) == as_list(lay.loc[t, "region_ids"])
    g8 = row(table, tag(8))  # TU_5 = [gene 7, gene 8] on "+"
    assert (g8["tu_id"], int(g8["tu_position"]), int(g8["tu_size"])) == ("TU_5", 2, 2)
    assert float(g8["uniop_p_prev"]) == 0.9
    assert is_na(g8["uniop_p_next"])  # last gene of its TU: no next neighbour
    g0 = row(table, tag(0))  # single-gene TU: no neighbour either side
    assert is_na(g0["uniop_p_prev"]) and is_na(g0["uniop_p_next"])


def test_regulatory_and_codon_values_copied_from_layers() -> None:
    fr = gt.frames_with()
    lay = fr["layers"]
    lay.loc[0, ["cai", "promoter_score", "sd_dg"]] = [0.75, 0.3, -4.5]
    lay.loc[1, "spacer_bp"] = 7
    lay.loc[1, "start_codon"] = "GTG"
    t = gt.build_world(frames=fr).assemble()
    g0, g1 = row(t, tag(0)), row(t, tag(1))
    assert (float(g0["cai"]), float(g0["promoter_score"]), float(g0["sd_dg"])) == (0.75, 0.3, -4.5)
    assert int(g1["spacer_bp"]) == 7 and g1["start_codon"] == "GTG"
    for c in ("tx_rate_max_tu", "terminator_conf_downstream", "tai", "rare_codons_2_10"):
        assert t[c].map(is_na).all(), c  # missing layers stay NA


def test_m5_columns_are_na(table: pd.DataFrame) -> None:
    """Example 9: ``island_id``, ``replichore``, ``replichore_position`` are filled by M5."""
    for c in ("island_id", "replichore", "replichore_position"):
        assert table[c].map(is_na).all(), c
    assert pd.api.types.is_float_dtype(table["replichore_position"].dtype)
    for c in ("island_id", "replichore"):
        dt = table[c].dtype
        assert pd.api.types.is_object_dtype(dt) or pd.api.types.is_string_dtype(dt), (c, dt)


def test_gc_local_dev_copied_from_layers() -> None:
    """Lead decision 4: M2-20 writes ``gc_local_dev``; it is copied unchanged."""
    fr = gt.frames_with()
    fr["layers"].insert(7, "gc_local_dev", [0.01 * i - 0.05 for i in range(20)])
    fr["layers"].loc[3, "gc_local_dev"] = float("nan")
    t = gt.build_world(frames=fr).assemble()
    for i in range(20):
        got = row(t, tag(i))["gc_local_dev"]
        if i == 3:
            assert is_na(got)
        else:
            assert math.isclose(float(got), 0.01 * i - 0.05, abs_tol=1e-15), i


def test_source_id_defaults_to_locus_tag_without_id_map(world: gt.World) -> None:
    """Lead decision 9: no C14 id map -> ``source_id = locus_tag``."""
    t = world.assemble(inputs=world.with_inputs(id_map=None))
    assert (t["source_id"] == t["locus_tag"]).all()


def test_gc_local_dev_na_without_layer_column() -> None:
    fr = gt.frames_with()
    assert "gc_local_dev" not in fr["layers"].columns  # optional C6 column (M2-20 decision 7)
    t = gt.build_world(frames=fr).assemble()
    assert t["gc_local_dev"].map(is_na).all()


def test_assemble_cog_def_keyword(world: gt.World, table: pd.DataFrame) -> None:
    """Lead decision 5: ``cog_def`` sets ``cog_category``; the default (None) leaves it NA."""
    assert table["cog_category"].map(is_na).all()
    t = world.assemble(cog_def=gt.cog_def_frame())
    assert row(t, tag(0))["cog_category"] == "J"  # COG0000
    assert row(t, tag(2))["cog_category"] == "KT"  # COG0002
    assert is_na(row(t, tag(4))["cog_category"])  # COG0004 not in the table
    assert is_na(row(t, tag(1))["cog_category"])  # no COG xref


def test_function_columns_in_table(table: pd.DataFrame) -> None:
    g0, g1 = row(table, tag(0)), row(table, tag(1))
    assert g0["cog_id"] == "COG0000" and as_list(g0["go_terms"]) == ["GO:0000000"]
    assert is_na(g1["cog_id"]) and as_list(g1["go_terms"]) == []
    assert as_list(g1["kegg_ko"]) == [] and as_list(g1["ec"]) == []


# ----------------------------------------------------------------------------- spec examples


def test_example_2_internal_gene_written_one_based(tmp_path: Path) -> None:
    fr = gt.frames_with()
    fr["genes"].loc[0, ["start", "end"]] = [0, 300]
    fr["tus"].loc[0, "start"] = 0
    t = gt.build_world(frames=fr).assemble()
    assert (int(row(t, tag(0))["start"]), int(row(t, tag(0))["end"])) == (0, 300)  # internal
    cells = _tsv(tmp_path, t).set_index("locus_tag")
    assert (cells.loc[tag(0), "start"], cells.loc[tag(0), "end"]) == ("1", "300")


def test_example_3_origin_crossing_gene(tmp_path: Path) -> None:
    fr = gt.frames_with()
    fr["replicons"].loc[0, "length"] = 50_000
    g = fr["genes"]
    g.loc[19, ["start", "end", "crosses_origin"]] = [49_900, 50_200, True]
    fr["layers"].at[19, "flags"] = ["crosses_origin"]
    t = gt.build_world(frames=fr).assemble()
    r = row(t, tag(19))
    assert (int(r["start"]), int(r["end"]), bool(r["crosses_origin"])) == (49_900, 50_200, True)
    assert "crosses_origin" in as_list(r["flags"])
    cells = _tsv(tmp_path, t).set_index("locus_tag")
    assert tuple(cells.loc[tag(19), ["start", "end", "crosses_origin"]]) == (
        "49901", "50200", "true",
    )  # fmt: skip
    assert "crosses_origin" in cells.loc[tag(19), "flags"].split(";")


def test_example_5_region_ids_list(tmp_path: Path, table: pd.DataFrame) -> None:
    assert as_list(row(table, tag(19))["region_ids"]) == ["amr_1", "prophage_1"]
    cells = _tsv(tmp_path, table).set_index("locus_tag")
    assert cells.loc[tag(19), "region_ids"] == "amr_1;prophage_1"
    assert cells.loc[tag(0), "region_ids"] == "NA"  # empty list


def test_example_6_function_pivot() -> None:
    from giap.genetable.assemble import function_columns

    xrefs = pd.DataFrame(
        [("G1", "COG", "COG0593"), ("G1", "GO", "GO:0005524"), ("G1", "EC", "3.6.4.12"),
         ("G1", "KEGG", "K02313"), ("G1", "UniRef", "UniRef90_P03004"),
         ("G2", "GO", "GO:0003677"), ("G2", "GO", "GO:0003677"), ("G2", "GO", "GO:0000001"),
         ("G3", "COG", "COG9999"), ("G4", "COG", "COG0002"),
         ("G5", "COG", "COG0593"), ("G5", "COG", "COG0002"), ("G5", "COG", "KT"),
         ("G6", "COG", "L")],
        columns=["locus_tag", "db", "id"],
    )  # fmt: skip
    cols = ["locus_tag", "cog_id", "cog_category", "go_terms", "kegg_ko", "ec"]
    with_cog = function_columns(xrefs, gt.cog_def_frame())
    assert list(with_cog.columns) == cols
    r = row(with_cog, "G1")
    assert (r["cog_id"], r["cog_category"]) == ("COG0593", "L")
    assert as_list(r["go_terms"]) == ["GO:0005524"]
    assert as_list(r["ec"]) == ["3.6.4.12"] and as_list(r["kegg_ko"]) == ["K02313"]
    r2 = row(with_cog, "G2")
    assert is_na(r2["cog_id"]) and is_na(r2["cog_category"])
    assert as_list(r2["go_terms"]) == ["GO:0000001", "GO:0003677"]  # sorted, unique
    assert row(with_cog, "G3")["cog_id"] == "COG9999" and is_na(row(with_cog, "G3")["cog_category"])
    assert row(with_cog, "G4")["cog_category"] == "KT"  # multi-letter category kept
    # decision 9: several COG IDs -> smallest ID and its category; bare categories ignored
    assert tuple(row(with_cog, "G5")[["cog_id", "cog_category"]]) == ("COG0002", "KT")
    g6 = with_cog[with_cog["locus_tag"] == "G6"]
    assert g6.empty or (is_na(g6.iloc[0]["cog_id"]) and is_na(g6.iloc[0]["cog_category"]))
    without = function_columns(xrefs, None)
    r = row(without, "G1")
    assert r["cog_id"] == "COG0593" and is_na(r["cog_category"])


def test_example_7_flags_in_vocabulary_order(world: gt.World, tmp_path: Path) -> None:
    t4 = tag(4)
    assert bool(row(world.detection.table, t4)["detected"])
    ev = world.evidence.copy()
    ev["evidence_flags"] = [
        ["tolerant_evidence", "strain_difference", "multi_identical"] if t == t4 else []
        for t in ev["locus_tag"]
    ]
    ab = world.abundance.table.copy()
    ab["abundance_flags"] = [["imputed_abundance"] if t == t4 else [] for t in ab["locus_tag"]]
    import dataclasses

    out = world.assemble(evidence=ev, abundance=dataclasses.replace(world.abundance, table=ab))
    want = ["multi_identical", "imputed_abundance", "tolerant_evidence", "strain_difference"]
    assert as_list(row(out, t4)["flags"]) == want
    cells = _tsv(tmp_path, out).set_index("locus_tag")
    assert cells.loc[t4, "flags"] == ";".join(want)
    assert cells.loc[tag(0), "flags"] == "NA"


def test_layer_flags_reach_the_table(table: pd.DataFrame) -> None:
    # gene 7: C6 ``pseudogene`` plus M4-14 ``top3_lt3`` (one quantified peptide)
    assert as_list(row(table, tag(7))["flags"]) == ["pseudogene", "top3_lt3"]


def test_example_8_partial_gene() -> None:
    fr = gt.frames_with()
    fr["genes"].loc[5, "partial"] = True
    fr["layers"].at[5, "flags"] = ["partial"]
    t = gt.build_world(frames=fr).assemble()
    r = row(t, tag(5))
    assert r["detection_status"] == "partial_excluded"
    assert bool(r["detected"]) is False and bool(r["partial"]) is True
    assert is_na(r["p_detect"]) and is_na(r["gene_q_value"])
    assert as_list(r["flags"]) == ["partial"]


def test_detectability_constant_on_eligible_genes_only() -> None:
    fr = gt.frames_with()
    fr["genes"].loc[5, "partial"] = True
    fr["layers"].at[5, "flags"] = ["partial"]
    t = gt.build_world(frames=fr, constant=True).assemble()
    for tg in t["locus_tag"]:
        flags = as_list(row(t, tg)["flags"])
        assert ("detectability_constant" in flags) is (tg != tag(5)), (tg, flags)
    assert as_list(row(t, tag(7))["flags"]) == ["pseudogene", "top3_lt3", "detectability_constant"]


# ----------------------------------------------------------------------------- sort order


def _two_contigs() -> dict[str, Any]:
    """Genes 15-19 on ``a_plasmid`` (listed after ``chr1`` in C1, so FASTA order != name order)."""
    fr = gt.frames_with()
    fr["replicons"] = pd.concat(
        [fr["replicons"], pd.DataFrame({"contig": ["a_plasmid"], "length": [5000],
                                        "topology": ["circular"], "kind": ["plasmid"]})],
        ignore_index=True,
    ).astype({"length": "int64"})  # fmt: skip
    g = fr["genes"]
    for k, i in enumerate(range(15, 20)):
        start = 100 + 350 * (4 - k)
        g.loc[i, ["contig", "start", "end"]] = ["a_plasmid", start, start + 300]
    fr["tus"].loc[fr["tus"]["tu_id"].isin(["TU_9", "TU_10", "TU_11"]), "contig"] = "a_plasmid"
    return fr


def test_sorted_by_contig_order_then_coordinates() -> None:
    w = gt.build_world(frames=_two_contigs())
    t = w.assemble()
    keys = gt.sort_key(t, ["chr1", "a_plasmid"])
    assert keys == sorted(keys)
    assert t["contig"].tolist() == ["chr1"] * 15 + ["a_plasmid"] * 5
    assert t["locus_tag"].tolist()[15:] == [tag(i) for i in (19, 18, 17, 16, 15)]
    assert row(t, tag(16))["replicon_type"] == "plasmid"


def test_sort_does_not_depend_on_input_order() -> None:
    w = gt.build_world(frames=_two_contigs())
    base = w.assemble()
    perm = np.random.default_rng(7).permutation(20)

    def shuffle(df: pd.DataFrame) -> pd.DataFrame:
        return df.iloc[perm].reset_index(drop=True)

    import dataclasses

    inputs = w.with_inputs(genes=shuffle(w.inputs.genes), layers=shuffle(w.inputs.layers),
                           id_map=shuffle(w.inputs.id_map))  # fmt: skip
    det = dataclasses.replace(w.detection, table=shuffle(w.detection.table))
    ab = dataclasses.replace(w.abundance, table=shuffle(w.abundance.table))
    out = w.assemble(inputs=inputs, cov=shuffle(w.cov), evidence=shuffle(w.evidence),
                     detection=det, abundance=ab)  # fmt: skip
    pd.testing.assert_frame_equal(out.reset_index(drop=True), base.reset_index(drop=True))


# ----------------------------------------------------------------------------- purity


def test_assemble_module_is_pure() -> None:
    """AGENTS.md: file I/O only in io/ and steps/; ``assemble.py`` is pure (acceptance)."""
    import ast

    tree = ast.parse(ASSEMBLE.read_text(encoding="utf-8"))
    banned_mod = {"os", "pathlib", "shutil", "io", "pyarrow.parquet", "giap.steps"}
    banned_call = {"open", "read_parquet", "to_parquet", "read_csv", "to_csv", "read_text",
                   "write_text", "read_bytes", "write_bytes", "read_tsv", "write_tsv"}  # fmt: skip
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            assert not {a.name for a in node.names} & banned_mod, ast.dump(node)
        elif isinstance(node, ast.ImportFrom):
            assert node.module not in banned_mod, node.module
        elif isinstance(node, ast.Call):
            f = node.func
            name = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", "")
            assert name not in banned_call, f"assemble.py calls {name}"
