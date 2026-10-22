"""M8-07 conformance on the merged M4 hand world (docs/specs/M8-07, acceptance "Conformance").

``gene_table_fixture.build_world()`` runs the merged M4-10 ... M4-14 producers on the M4-06
contract frames (20 genes on circular ``chr1``, a sORF, a pseudogene, ``prophage_1`` on genes
10-19, every regulatory layer NA) and ``assemble_genes`` builds the §9.3 gene table. The M8-07
arguments come straight from those frames (C1, C2, C5, C6, C7 + C8); there is no M5-09 gene
context, so ``replichore_position`` is NA. The real-path sim100k version is
tests/integration/test_determinants_features_sim100k.py (slow).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_features_data as d
import gene_table_fixture as gt
from determinants_features_data import is_na, mod


@pytest.fixture(scope="module")
def world() -> dict[str, Any]:
    w = gt.build_world()
    fr = w.frames
    return {
        "genes_tbl": w.assemble(),
        "genes": d.genes_from_frame(fr["genes"]),
        "tus": d.tus_from_frame(fr["tus"]),
        "layers": fr["layers"],
        "regions": d.regions_arg(fr["regions"], fr["region_genes"]),
        "c2": fr["genes"],
        "replicons": d.replicons_from_frame(fr["replicons"]),
    }


def _features(w: dict[str, Any]) -> Any:
    return mod().feature_frame(w["genes_tbl"], w["genes"], w["tus"], w["layers"], w["regions"],
                               None, w["replicons"])  # fmt: skip


def test_every_eligible_gene_has_one_row_in_gene_table_order(world: dict[str, Any]) -> None:
    df = _features(world)
    assert df["locus_tag"].astype(str).tolist() == d.eligible_tags(world["genes_tbl"])
    assert list(df.columns) == list(d.TSV_COLUMNS)


@pytest.mark.parametrize(("feature", "source"),
                         d.CONFORMANCE_SOURCES)  # fmt: skip
def test_columns_equal_their_gene_table_sources(
    world: dict[str, Any], feature: str, source: str
) -> None:
    df = _features(world)
    tbl = world["genes_tbl"].set_index("locus_tag")
    for tag, r in d.by_tag(df).items():
        assert d.same(r[feature], tbl.loc[tag, source]), (feature, tag, r[feature])


def test_prophage_membership_matches_c8(world: dict[str, Any]) -> None:
    rg = world["regions"]
    members = set(rg.loc[rg["region_type"] == "prophage", "locus_tag"].astype(str))
    assert members, "the hand world plants prophage_1"
    for tag, r in d.by_tag(_features(world)).items():
        assert bool(r["in_prophage"]) is (tag in members), tag
        assert bool(r["in_bgc"]) is False, tag


def test_nterm_residue_and_pseudogene(world: dict[str, Any]) -> None:
    c2 = world["c2"].set_index("locus_tag")
    for tag, r in d.by_tag(_features(world)).items():
        aa = c2.loc[tag, "aa"]
        if not isinstance(aa, str) or len(aa) < 2:
            assert is_na(r["nterm_residue"]), tag
        else:
            assert r["nterm_residue"] == "A", (tag, aa[:3])  # "MAAA..." -> Met removed


def test_distances_follow_c5_transcription_order(world: dict[str, Any]) -> None:
    rows = d.by_tag(_features(world))
    by_tag = {g.locus_tag: g for g in world["genes"]}
    for t in world["tus"]:
        p0 = d.start_codon_pos(by_tag[t.genes[0]])
        for tag in t.genes:
            if tag not in rows:
                continue
            want = abs(d.start_codon_pos(by_tag[tag]) - p0)
            assert rows[tag]["dist_to_tu_start_bp"] == want, (t.tu_id, tag)


def test_no_gene_context_means_na_replichore_position(world: dict[str, Any]) -> None:
    assert _features(world)["replichore_position"].isna().all()


def test_model_matrix_drops_the_all_na_layers_of_the_hand_world(world: dict[str, Any]) -> None:
    """Every regulatory layer is NA in the hand world: none of them reaches the model."""
    df = _features(world)
    x, blocks = mod().model_matrix(df)
    assert set(blocks) == set(x.columns)
    for c in ("promotech_max", "promotech_tu", "tir_ostir", "sd_dg", "cai", "replichore_position"):
        assert df[c].isna().all() and c not in x.columns, c
    assert {"tu_position", "dist_to_tu_start_bp", "length_aa"} <= set(x.columns)
