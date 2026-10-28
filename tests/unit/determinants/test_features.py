"""M8-07 ``determinants_features.tsv``: examples 1-10 and every degenerate row of
docs/specs/M8-07-determinants-features.md (GIPPA_SPEC.md §7.8.1, §7.8.2 step 1, §7.8.3).

The world is ``determinants_features_data.toy()`` (13 genes; see that module for the layout and
for the argument conventions this file assumes). Examples 1-4 are the §7.8.3 "Unit: the
Met-processing rule" test. Run flags (``layer_missing:<layer>``) are the step's business (M8-10):
the pure functions are only checked to keep the TSV column and drop it from the model matrix.
"""

from __future__ import annotations

import copy
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_features_data as d
from determinants_features_data import is_na, mod

# ---------------------------------------------------------------------------------------------
# nterm_residue: the Met-processing rule (examples 1-4, decision 19)


@pytest.mark.parametrize(
    ("aa", "want"),
    [
        ("MAKVLE", "A"),  # example 1
        ("MKLVAE", "M"),  # example 2
        ("MP", "P"),  # example 3
        ("VKLAE", "V"),  # example 4
    ],
)
def test_met_processing_rule_examples_1_to_4(aa: str, want: str) -> None:
    assert mod().nterm_residue(aa) == want


@pytest.mark.parametrize("second", list("ACGPSTV"))
def test_met_is_removed_for_every_small_second_residue(second: str) -> None:
    assert mod().nterm_residue("M" + second + "KLE") == second


@pytest.mark.parametrize("second", [r for r in d.RESIDUES if r not in "ACGPSTV"])
def test_met_is_kept_for_every_other_second_residue(second: str) -> None:
    assert mod().nterm_residue("M" + second + "KLE") == "M"


@pytest.mark.parametrize("aa", ["VA", "LS", "AG", "TP"])
def test_no_initiator_met_gives_residue_1_even_when_residue_2_is_small(aa: str) -> None:
    """Degenerate row: residue 1 not M -> residue 1 (the rule needs an initiator Met)."""
    assert mod().nterm_residue(aa) == aa[0]


@pytest.mark.parametrize("aa", ["M", "V", ""])
def test_proteins_shorter_than_two_residues_are_na(aa: str) -> None:
    """Degenerate rows: ``aa`` = "M" -> NA; shorter than 2 aa -> NA; empty (pseudogene) -> NA."""
    assert is_na(mod().nterm_residue(aa))


# ---------------------------------------------------------------------------------------------
# dist_to_tu_start (examples 5-7)


def _dist(genes: list, tus: list, reps: dict) -> dict:
    return dict(mod().dist_to_tu_start(genes, tus, reps))


def test_example_5_plus_strand_tu() -> None:
    reps = {"c": d.replicon("c", 5000, "linear")}
    a, b = d.gene("a", "c", 1000, 1300, "+"), d.gene("b", "c", 1450, 1750, "+")
    got = _dist([a, b], [d.tu("TU_1", [a, b])], reps)
    assert got == {"a": 0, "b": 450}


def test_example_6_minus_strand_tu_counts_from_the_right_hand_start_codon() -> None:
    reps = {"c": d.replicon("c", 10_000, "linear")}
    first, second = d.gene("r", "c", 5000, 5600, "-"), d.gene("l", "c", 4300, 4950, "-")
    got = _dist([second, first], [d.tu("TU_1", [first, second])], reps)
    assert got == {"r": 0, "l": 650}


def test_example_7_circular_tu_across_the_origin_wraps() -> None:
    reps = {"c": d.replicon("c", 10_000, "circular")}
    a, b = d.gene("a", "c", 9800, 9950, "+"), d.gene("b", "c", 200, 500, "+")
    got = _dist([b, a], [d.tu("TU_1", [a, b], 10_000)], reps)
    assert got == {"a": 0, "b": 400}


def test_minus_strand_tu_across_the_origin_wraps() -> None:
    """Start codons at 399 and 9899 on a 10-kb ring, transcribed leftwards: (399 - 9899) mod L."""
    reps = {"c": d.replicon("c", 10_000, "circular")}
    first, second = d.gene("a", "c", 100, 400, "-"), d.gene("b", "c", 9500, 9900, "-")
    got = _dist([first, second], [d.tu("TU_1", [first, second], 10_000)], reps)
    assert got == {"a": 0, "b": 500}


def test_first_gene_of_every_tu_is_zero_and_single_gene_tus_are_zero() -> None:
    w = d.toy()
    got = _dist(w.genes, w.tus, w.replicons)
    for t in w.tus:
        assert got[t.genes[0]] == 0, t.tu_id
    assert got == d.EXPECTED_DIST


def test_distances_are_plain_ints() -> None:
    w = d.toy()
    got = _dist(w.genes, w.tus, w.replicons)
    assert all(type(v) is int for v in got.values()), {k: type(v) for k, v in got.items()}


def test_gene_in_no_tu_has_no_distance() -> None:
    """Degenerate: a gene in no TU (should not occur after M2-12) -> NA."""
    w = d.toy()
    tus = [t for t in w.tus if t.tu_id != "TU_1"]
    got = mod().dist_to_tu_start(w.genes, tus, w.replicons)
    for tag in ("G01", "G02"):
        assert is_na(got.get(tag)), (tag, got.get(tag))
    assert got["G04"] == 650


# ---------------------------------------------------------------------------------------------
# feature_frame: rows, columns, dtypes


def test_feature_columns_constant_is_the_tsv_layout() -> None:
    """Identity columns first (spec step 1), then the §7.8.1 features in table order."""
    assert tuple(mod().FEATURE_COLUMNS) == d.TSV_COLUMNS


def test_feature_frame_columns_in_order() -> None:
    assert list(d.toy().features().columns) == list(d.TSV_COLUMNS)


def test_one_row_per_eligible_gene_in_section_6_5_order() -> None:
    """§7.3: partial genes and genes without a finite ``p_detect`` have no row."""
    df = d.toy().features()
    assert df["locus_tag"].astype(str).tolist() == list(d.ELIGIBLE_ORDER)
    assert not set(d.NOT_ELIGIBLE) & set(df["locus_tag"].astype(str))


def test_row_order_does_not_depend_on_the_gene_table_row_order() -> None:
    w = d.toy()
    want = w.features()
    w2 = w.copy()
    w2.genes_tbl = w2.genes_tbl.iloc[::-1].reset_index(drop=True)
    w2.layers = w2.layers.sample(frac=1.0, random_state=3).reset_index(drop=True)
    w2.genes = list(reversed(w2.genes))
    w2.tus = list(reversed(w2.tus))
    got = w2.features()
    pd.testing.assert_frame_equal(
        got.reset_index(drop=True), want.reset_index(drop=True), check_exact=True
    )


def test_no_eligible_gene_gives_an_empty_frame_with_every_column() -> None:
    w = d.toy()
    w.genes_tbl["p_detect"] = np.nan
    df = w.features()
    assert len(df) == 0
    assert list(df.columns) == list(d.TSV_COLUMNS)


@pytest.mark.parametrize("col", d.INT_FEATURES)
def test_integer_features_have_an_integer_dtype(col: str) -> None:
    df = d.toy().features()
    assert pd.api.types.is_integer_dtype(df[col]), (col, df[col].dtype)


@pytest.mark.parametrize("col", d.FLOAT_FEATURES)
def test_float_features_have_a_float_dtype(col: str) -> None:
    df = d.toy().features()
    assert pd.api.types.is_float_dtype(df[col]), (col, df[col].dtype)


@pytest.mark.parametrize("col", d.BOOL_COLUMNS)
def test_boolean_columns_have_a_boolean_dtype(col: str) -> None:
    df = d.toy().features()
    assert pd.api.types.is_bool_dtype(df[col]), (col, df[col].dtype)


@pytest.mark.parametrize("col", ("locus_tag", "tu_id", *d.STR_FEATURES))
def test_string_columns_hold_strings_or_na(col: str) -> None:
    df = d.toy().features()
    assert not pd.api.types.is_numeric_dtype(df[col]), (col, df[col].dtype)
    assert all(is_na(v) or isinstance(v, str) for v in df[col]), (col, df[col].tolist())


def test_integer_dtype_survives_na_cells() -> None:
    """``dist_to_tu_start_bp`` NA for a gene in no TU keeps the column an integer column."""
    w = d.toy()
    w.tus = [t for t in w.tus if t.tu_id != "TU_1"]
    df = w.features()
    assert pd.api.types.is_integer_dtype(df["dist_to_tu_start_bp"])
    rows = d.by_tag(df)
    assert is_na(rows["G01"]["dist_to_tu_start_bp"]) and is_na(rows["G02"]["dist_to_tu_start_bp"])
    assert rows["G04"]["dist_to_tu_start_bp"] == 650


# ---------------------------------------------------------------------------------------------
# feature_frame: values


@pytest.mark.parametrize(("feature", "source"), sorted(d.GENE_TABLE_SOURCES.items()))
def test_feature_is_copied_from_its_gene_table_column(feature: str, source: str) -> None:
    """Spec step 2 name map; the TSV holds raw values (no transform)."""
    w = d.toy()
    rows = d.by_tag(w.features())
    for tag, r in rows.items():
        want = d.gene_table_row(w, tag)[source]
        assert d.same(r[feature], want), (feature, tag, r[feature], want)


@pytest.mark.parametrize(("col", "source"), sorted(d.IDENTITY_SOURCES.items()))
def test_identity_columns_are_copied_from_the_gene_table(col: str, source: str) -> None:
    w = d.toy()
    for tag, r in d.by_tag(w.features()).items():
        want = d.gene_table_row(w, tag)[source]
        assert d.same(r[col], want), (col, tag, r[col], want)


def test_n_promoters_tu_to_gene_is_c6_n_promoter_windows_tu() -> None:
    """Lead decision 1 (question 2)."""
    w = d.toy()
    for tag, r in d.by_tag(w.features()).items():
        want = d.layers_row(w, tag)["n_promoter_windows_tu"]
        assert d.same(r["n_promoters_tu_to_gene"], want), (tag, r["n_promoters_tu_to_gene"], want)


def test_dist_to_tu_start_bp_column_examples_5_to_7() -> None:
    rows = d.by_tag(d.toy().features())
    for tag, r in rows.items():
        assert r["dist_to_tu_start_bp"] == d.EXPECTED_DIST[tag], tag


def test_nterm_residue_column_uses_c2_aa() -> None:
    """Examples 1-4 on real rows, ``aa`` = "M" (G05) and the empty pseudogene ``aa`` (G09) NA."""
    rows = d.by_tag(d.toy().features())
    for tag, r in rows.items():
        want = d.EXPECTED_NTERM[tag]
        if want is None:
            assert is_na(r["nterm_residue"]), (tag, r["nterm_residue"])
        else:
            assert r["nterm_residue"] == want, (tag, r["nterm_residue"])


def test_ribosomal_protein_row_is_present_and_flagged() -> None:
    """Degenerate row: the CAI reference set stays in the TSV with ``ribosomal_protein`` true."""
    rows = d.by_tag(d.toy().features())
    assert "G07" in rows
    assert bool(rows["G07"]["ribosomal_protein"]) is True
    others = [t for t in rows if t != "G07"]
    assert all(bool(rows[t]["ribosomal_protein"]) is False for t in others)


def test_example_8_prophage_and_bgc_membership() -> None:
    rows = d.by_tag(d.toy().features())
    assert bool(rows["G08"]["in_prophage"]) is True and bool(rows["G08"]["in_bgc"]) is False
    assert bool(rows["G09"]["in_bgc"]) is True and bool(rows["G09"]["in_prophage"]) is False
    for tag in set(rows) - {"G08", "G09"}:
        assert bool(rows[tag]["in_prophage"]) is False, tag
        assert bool(rows[tag]["in_bgc"]) is False, tag


def test_other_region_types_do_not_set_either_flag() -> None:
    w = d.toy()
    extra = pd.DataFrame(
        {"region_id": ["crispr_1", "plasmid_1", "amr_1"], "locus_tag": ["G01", "G12", "G13"],
         "overlap_frac": [1.0, 1.0, 1.0], "subset": [None, None, None],
         "core_biosynthetic": [False, False, False],
         "region_type": ["crispr", "plasmid", "amr"]}
    )  # fmt: skip
    assert w.regions is not None
    w.regions = pd.concat([w.regions, extra], ignore_index=True)
    rows = d.by_tag(w.features())
    for tag in ("G01", "G12", "G13"):
        assert bool(rows[tag]["in_prophage"]) is False and bool(rows[tag]["in_bgc"]) is False


def test_no_regions_layer_gives_na_membership() -> None:
    """Degenerate row: no C7/C8 -> ``in_prophage``, ``in_bgc`` NA (not false)."""
    w = d.toy()
    w.regions = None
    df = w.features()
    assert df["in_prophage"].isna().all() and df["in_bgc"].isna().all()
    assert len(df) == len(d.ELIGIBLE_ORDER)


def test_replichore_position_comes_from_gene_context() -> None:
    w = d.toy()
    assert w.gene_context is not None
    ctx = w.gene_context.set_index("locus_tag")["replichore_position"]
    for tag, r in d.by_tag(w.features()).items():
        assert d.same(r["replichore_position"], ctx[tag]), tag


def test_linear_replicon_without_oric_has_na_replichore_position() -> None:
    rows = d.by_tag(d.toy().features())
    assert is_na(rows["G12"]["replichore_position"])
    assert rows["G12"]["replicon_type"] == "plasmid"


def test_no_gene_context_gives_na_replichore_position() -> None:
    w = d.toy()
    w.gene_context = None
    assert w.features()["replichore_position"].isna().all()


def test_example_10_tir_ostir_zero_is_kept_raw_in_the_tsv() -> None:
    rows = d.by_tag(d.toy().features())
    assert rows["G02"]["tir_ostir"] == 0.0
    assert rows["G03"]["tir_ostir"] == -2.0


def test_tsv_has_no_transforms() -> None:
    """``tir_ostir`` and ``length_aa`` raw; ``start_codon`` / ``replicon_type`` not one-hot."""
    w = d.toy()
    rows = d.by_tag(w.features())
    assert rows["G06"]["start_codon"] == "CTG"
    assert rows["G04"]["start_codon"] == "GTG"
    assert rows["G01"]["length_aa"] == d.gene_table_row(w, "G01")["length_aa"]
    assert not any(c.startswith(("start_codon_", "replicon_type_", "nterm_residue_"))
                   for c in w.features().columns)  # fmt: skip


def test_example_9_missing_promotech_keeps_the_columns_all_na() -> None:
    w = d.toy().blank(d.LAYER_GROUPS["promoters"])
    df = w.features()
    for c in d.LAYER_FEATURES["promoters"]:
        assert c in df.columns and df[c].isna().all(), c
    assert df["tx_rate_max_tu"].notna().all()


# ---------------------------------------------------------------------------------------------
# BLOCKS


def test_blocks_are_the_five_attribution_blocks_in_order() -> None:
    """``determinants_attribution.tsv``: ``shap_transcription`` ... ``shap_context``."""
    assert tuple(mod().BLOCKS) == d.BLOCK_NAMES


@pytest.mark.parametrize("block", d.BLOCK_NAMES)
def test_block_model_columns(block: str) -> None:
    cols = tuple(mod().BLOCKS[block])
    assert len(cols) == len(set(cols)), cols
    assert set(cols) == d.EXPECTED_BLOCKS[block]


def test_blocks_are_disjoint() -> None:
    seen: dict[str, str] = {}
    for block, cols in mod().BLOCKS.items():
        for c in cols:
            assert c not in seen, f"{c} in {seen.get(c)} and {block}"
            seen[c] = block


# ---------------------------------------------------------------------------------------------
# model_matrix


def _mm(df: pd.DataFrame | None = None) -> tuple[pd.DataFrame, dict[str, str]]:
    df = d.toy().features() if df is None else df
    x, blocks = mod().model_matrix(df)
    return x, dict(blocks)


def test_model_matrix_one_row_per_feature_row_in_order() -> None:
    df = d.toy().features()
    x, _ = _mm(df)
    assert len(x) == len(df)
    if not isinstance(x.index, pd.RangeIndex):
        assert x.index.astype(str).tolist() == df["locus_tag"].astype(str).tolist()


def test_model_matrix_has_no_identity_columns() -> None:
    x, _ = _mm()
    assert not set(x.columns) & set(d.IDENTITY)


def test_model_matrix_is_numeric() -> None:
    x, _ = _mm()
    for c in x.columns:
        assert pd.api.types.is_numeric_dtype(x[c]), (c, x[c].dtype)
    np.asarray(x, dtype=float)


def test_model_matrix_metadata_maps_every_column_to_its_block() -> None:
    x, blocks = _mm()
    assert set(blocks) == set(x.columns)
    allowed = mod().BLOCKS
    for c, b in blocks.items():
        assert c in allowed[b], (c, b)


def test_model_matrix_full_world_keeps_every_present_column() -> None:
    """With every layer present (toy world) all block columns appear."""
    x, _ = _mm()
    want = set().union(*d.EXPECTED_BLOCKS.values())
    assert set(x.columns) == want


def _col(x: pd.DataFrame, c: str) -> np.ndarray:
    return np.asarray(x[c], dtype=float)


def test_example_10_tir_ostir_is_logged_and_non_positive_is_na() -> None:
    df = d.toy().features()
    x, _ = _mm(df)
    raw = df["tir_ostir"].to_numpy(float)
    got = _col(x, "tir_ostir")
    for i, tag in enumerate(df["locus_tag"]):
        if raw[i] > 0:
            assert got[i] == pytest.approx(math.log(raw[i]), rel=1e-12), tag
        else:
            assert np.isnan(got[i]), (tag, raw[i], got[i])


def test_length_aa_is_natural_log() -> None:
    """Lead decision 1 (question 3): natural logs."""
    df = d.toy().features()
    x, _ = _mm(df)
    np.testing.assert_allclose(_col(x, "length_aa"), np.log(df["length_aa"].to_numpy(float)),
                               rtol=1e-12)  # fmt: skip


@pytest.mark.parametrize(
    "c",
    [c for c in d.FEATURES if c not in ("tir_ostir", "length_aa", "start_codon", "replicon_type",
                                        "nterm_residue", "in_prophage", "in_bgc")],
)  # fmt: skip
def test_other_numeric_features_enter_the_model_untransformed(c: str) -> None:
    df = d.toy().features()
    x, _ = _mm(df)
    np.testing.assert_array_equal(_col(x, c), df[c].to_numpy(float, na_value=np.nan))


@pytest.mark.parametrize("c", ["in_prophage", "in_bgc"])
def test_region_membership_is_zero_one(c: str) -> None:
    df = d.toy().features()
    x, _ = _mm(df)
    np.testing.assert_array_equal(_col(x, c), df[c].astype(float).to_numpy())


def _one_hot(x: pd.DataFrame, feature: str, levels: tuple[str, ...]) -> pd.DataFrame:
    return pd.DataFrame({lv: _col(x, f"{feature}_{lv}") for lv in levels})


def test_start_codon_one_hot_with_other() -> None:
    """Degenerate row: a start codon outside {ATG, GTG, TTG} -> ``other`` (CTG on G06)."""
    df = d.toy().features()
    x, _ = _mm(df)
    oh = _one_hot(x, "start_codon", d.START_LEVELS)
    for i, codon in enumerate(df["start_codon"]):
        level = codon if codon in ("ATG", "GTG", "TTG") else "other"
        assert oh.iloc[i].tolist() == [1.0 if lv == level else 0.0 for lv in d.START_LEVELS], (
            df["locus_tag"].iloc[i],
            codon,
        )


def test_replicon_type_one_hot() -> None:
    df = d.toy().features()
    x, _ = _mm(df)
    oh = _one_hot(x, "replicon_type", d.REPLICON_LEVELS)
    for i, kind in enumerate(df["replicon_type"]):
        assert oh.iloc[i].tolist() == [1.0 if lv == kind else 0.0 for lv in d.REPLICON_LEVELS]


def test_nterm_residue_one_hot_over_twenty_residues() -> None:
    """Rare levels kept: every residue has a column even when no gene has it."""
    df = d.toy().features()
    x, _ = _mm(df)
    oh = _one_hot(x, "nterm_residue", tuple(d.RESIDUES))
    for i, res in enumerate(df["nterm_residue"]):
        row = oh.iloc[i]
        if is_na(res):
            assert not (row == 1.0).any(), (df["locus_tag"].iloc[i], row[row == 1.0])
        else:
            assert row.tolist() == [1.0 if r == res else 0.0 for r in d.RESIDUES]


@pytest.mark.parametrize("layer", sorted(d.LAYER_FEATURES))
def test_all_na_layer_columns_are_dropped_from_the_model(layer: str) -> None:
    """Spec step 4 / example 9: kept in the TSV, absent from ``model_matrix`` and its metadata."""
    df = d.toy().blank(d.LAYER_GROUPS[layer]).features()
    x, blocks = _mm(df)
    for c in d.LAYER_FEATURES[layer]:
        assert c in df.columns and df[c].isna().all(), c
        assert c not in x.columns and c not in blocks, c
    assert set(x.columns) == set(blocks)


def test_example_9_promotech_missing_drops_three_transcription_columns_only() -> None:
    df = d.toy().blank(d.LAYER_GROUPS["promoters"]).features()
    _, blocks = _mm(df)
    assert d.model_columns_of(blocks, "transcription") == {
        "tx_rate_max_tu", "tx_rate_sum_tu", "tu_position", "dist_to_tu_start_bp",
    }  # fmt: skip


def test_a_block_with_no_column_left_is_absent_from_the_metadata() -> None:
    df = d.toy().blank(d.LAYER_GROUPS["codon"]).features()
    _, blocks = _mm(df)
    assert "codon" not in set(blocks.values())


def test_all_na_start_codon_drops_its_one_hot_columns() -> None:
    df = d.toy().blank(["start_codon"]).features()
    x, blocks = _mm(df)
    assert not [c for c in x.columns if c.startswith("start_codon")]
    assert not [c for c in blocks if c.startswith("start_codon")]


def test_no_regions_drops_membership_from_the_model() -> None:
    w = d.toy()
    w.regions = None
    x, blocks = _mm(w.features())
    for c in ("in_prophage", "in_bgc"):
        assert c not in x.columns and c not in blocks


def test_no_gene_context_drops_replichore_position_from_the_model() -> None:
    w = d.toy()
    w.gene_context = None
    x, _ = _mm(w.features())
    assert "replichore_position" not in x.columns


def test_partly_na_column_is_kept() -> None:
    """Only an all-NA column is dropped; ``replichore_position`` is NA on the plasmid gene only."""
    x, _ = _mm()
    assert "replichore_position" in x.columns
    assert np.isnan(_col(x, "replichore_position")).sum() == 1


# ---------------------------------------------------------------------------------------------
# determinism, purity


def test_feature_frame_is_deterministic() -> None:
    a, b = d.toy().features(), d.toy().features()
    pd.testing.assert_frame_equal(a, b, check_exact=True)


def test_model_matrix_is_deterministic() -> None:
    df = d.toy().features()
    (x1, b1), (x2, b2) = _mm(df), _mm(df.copy())
    pd.testing.assert_frame_equal(x1, x2, check_exact=True)
    assert b1 == b2 and list(b1) == list(b2)


def test_feature_frame_does_not_mutate_its_inputs() -> None:
    w = d.toy()
    before = copy.deepcopy(w)
    w.features()
    pd.testing.assert_frame_equal(w.genes_tbl, before.genes_tbl)
    pd.testing.assert_frame_equal(w.layers, before.layers)
    assert w.regions is not None and before.regions is not None
    pd.testing.assert_frame_equal(w.regions, before.regions)
    assert w.gene_context is not None and before.gene_context is not None
    pd.testing.assert_frame_equal(w.gene_context, before.gene_context)
    assert w.genes == before.genes and w.tus == before.tus
    assert w.replicons == before.replicons


def test_model_matrix_does_not_mutate_the_feature_frame() -> None:
    df = d.toy().features()
    before = df.copy(deep=True)
    _mm(df)
    pd.testing.assert_frame_equal(df, before, check_exact=True)
