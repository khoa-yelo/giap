"""Hypothesis properties of M8-07 (docs/specs/M8-07-determinants-features.md, acceptance list).

* **Met rule** -- ``nterm_residue`` is residue 2 exactly when residue 1 is M and residue 2 is in
  {A, C, G, P, S, T, V}, residue 1 otherwise, NA below 2 aa; it only reads the first two residues;
* **distance oracle** -- genes laid out at transcription offsets u_0 = 0 < u_1 < ... from a TU
  start (either strand, linear or circular, anywhere on the ring including across the origin)
  get ``dist_to_tu_start_bp`` = u_i;
* **mirror** -- reflecting the genome (x -> L - x, strands flipped) leaves every distance unchanged;
* **blocks** -- whatever layers are missing, every ``model_matrix`` column maps to exactly one
  block, the metadata keys equal the columns, and exactly the all-NA feature columns are dropped.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_features_data as d
from determinants_features_data import is_na, mod

SETTINGS = settings(
    max_examples=60,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
SMALL = set("ACGPSTV")


@SETTINGS
@given(st.text(alphabet=d.RESIDUES, max_size=12))
def test_met_rule(aa: str) -> None:
    got = mod().nterm_residue(aa)
    if len(aa) < 2:
        assert is_na(got)
    elif aa[0] == "M" and aa[1] in SMALL:
        assert got == aa[1]
    else:
        assert got == aa[0]


@SETTINGS
@given(st.text(alphabet=d.RESIDUES, min_size=2, max_size=4),
       st.text(alphabet=d.RESIDUES, max_size=20))  # fmt: skip
def test_met_rule_reads_only_the_first_two_residues(head: str, tail: str) -> None:
    f = mod().nterm_residue
    assert f(head[:2] + tail) == f(head[:2])


@st.composite
def tu_layouts(draw: Any) -> dict[str, Any]:
    """One TU of k genes on a replicon of length L, with its expected distances."""
    circular = draw(st.booleans())
    strand = draw(st.sampled_from("+-"))
    k = draw(st.integers(1, 5))
    lens = [3 * draw(st.integers(20, 200)) for _ in range(k)]
    gaps = [draw(st.integers(-20, 300)) for _ in range(k - 1)]
    u = [0]
    for i in range(k - 1):
        u.append(u[-1] + lens[i] + max(gaps[i], -min(lens[i], lens[i + 1]) // 2 + 1))
    span = u[-1] + lens[-1]
    length = span + draw(st.integers(1, 20_000))
    if circular:
        o = draw(st.integers(0, length - 1))
    else:
        o = draw(st.integers(0, length - span)) if strand == "+" else (
            draw(st.integers(span - 1, length - 1)))  # fmt: skip
    genes = []
    for i in range(k):
        if strand == "+":
            s = (o + u[i]) % length if circular else o + u[i]
            e = s + lens[i]
        else:
            sc = (o - u[i]) % length if circular else o - u[i]  # start codon
            s = sc - lens[i] + 1
            if s < 0:  # circular only: store unwrapped across the origin
                s += length
            e = s + lens[i]
        genes.append(d.gene(f"g{i}", "c", s, e, strand))
    rep = d.replicon("c", length, "circular" if circular else "linear")
    return {"genes": genes, "length": length, "circular": circular, "rep": rep,
            "want": {f"g{i}": u[i] for i in range(k)}}  # fmt: skip


def _distances(genes: list[Any], length: int, rep: Any) -> dict[str, Any]:
    return dict(mod().dist_to_tu_start(genes, [d.tu("TU_1", genes, length)], {"c": rep}))


@SETTINGS
@given(tu_layouts())
def test_distance_oracle(lay: dict[str, Any]) -> None:
    got = _distances(lay["genes"], lay["length"], lay["rep"])
    assert got == lay["want"]


@SETTINGS
@given(tu_layouts())
def test_mirroring_the_genome_leaves_distances_unchanged(lay: dict[str, Any]) -> None:
    before = _distances(lay["genes"], lay["length"], lay["rep"])
    mirrored = [d.mirror_gene(g, lay["length"], lay["circular"]) for g in lay["genes"]]
    after = _distances(mirrored, lay["length"], lay["rep"])
    assert after == before


@SETTINGS
@given(tu_layouts(), st.randoms(use_true_random=False))
def test_distances_do_not_depend_on_the_gene_list_order(lay: dict[str, Any], rnd: Any) -> None:
    genes = list(lay["genes"])
    tu = d.tu("TU_1", genes, lay["length"])
    rnd.shuffle(genes)
    got = dict(mod().dist_to_tu_start(genes, [tu], {"c": lay["rep"]}))
    assert got == lay["want"]


_BLANKABLE = (*d.LAYER_GROUPS, "start_codon", "gc_local_dev", "regions", "gene_context")
_WORLD = d.toy()


def _world_without(missing: frozenset[str]) -> d.World:
    cols: list[str] = []
    for m in missing:
        if m in d.LAYER_GROUPS:
            cols += list(d.LAYER_GROUPS[m])
        elif m in ("start_codon", "gc_local_dev"):
            cols.append(m)
    w = _WORLD.blank(cols)
    if "regions" in missing:
        w.regions = None
    if "gene_context" in missing:
        w.gene_context = None
    return w


@SETTINGS
@given(st.frozensets(st.sampled_from(_BLANKABLE)))
def test_every_model_column_maps_to_exactly_one_block(missing: frozenset[str]) -> None:
    f = mod()
    df = _world_without(missing).features()
    x, blocks = f.model_matrix(df)
    blocks = dict(blocks)
    assert list(df.columns) == list(d.TSV_COLUMNS)
    assert set(blocks) == set(x.columns)
    for c in x.columns:
        owners = [b for b, cols in f.BLOCKS.items() if c in cols]
        assert owners == [blocks[c]], (c, owners, blocks[c])


@SETTINGS
@given(st.frozensets(st.sampled_from(_BLANKABLE)))
def test_exactly_the_all_na_features_are_dropped(missing: frozenset[str]) -> None:
    df = _world_without(missing).features()
    x, _ = mod().model_matrix(df)
    cols = set(x.columns)
    for feat in d.FEATURES:
        encoded = {c for c in d.EXPECTED_BLOCKS_UNION if c == feat or c.startswith(feat + "_")}
        all_na = bool(df[feat].isna().all())
        if all_na:
            assert not encoded & cols, (feat, encoded & cols)
        else:
            assert encoded <= cols, (feat, encoded - cols)


@SETTINGS
@given(st.frozensets(st.sampled_from(_BLANKABLE)))
def test_missing_layers_never_change_rows_or_other_columns(missing: frozenset[str]) -> None:
    full = _WORLD.features()
    df = _world_without(missing).features()
    assert df["locus_tag"].tolist() == full["locus_tag"].tolist()
    untouched = ["tu_position", "dist_to_tu_start_bp", "length_aa", "mw", "pi", "gravy",
                 "nterm_residue", "tm_segments", "replicon_type", "ribosomal_protein"]  # fmt: skip
    pd.testing.assert_frame_equal(df[untouched], full[untouched], check_exact=True)
