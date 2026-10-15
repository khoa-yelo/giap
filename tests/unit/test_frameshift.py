"""Tolerant matching, ±1-nt frameshift splits (issue M3-08; GIAP_SPEC.md §7.2.4).

Examples 1-7 of docs/specs/M3-08-frameshift-splits.md, its edge cases (stop at the split point,
δ = -1 overlap, origin-crossing splits on circular replicons, equal-cost ordinary tolerant hits),
sides of ``min_side_aa`` (4) residues — which M3-08 checks with its own anchored check, never with
``tolerant_dp.verify`` (M3-06 lead decision 1: ``verify`` raises for queries of ≤ 4 aa) —, decoys
searched identically (each peptide against its own decoy type), the M3-10 Q1 target/decoy split,
config handling and determinism.

Conventions resolved by the test author (spec silent), shared with
tests/helpers/oracle_frameshift.py:

- rows follow the M3-07 driver: ``HIT_COLUMNS``, ``route = "tolerant"``, ``query`` the I→L
  normalised text, ``decoy_type`` NA on target rows, sorted by peptide, search (target first),
  decoy type, contig order, start, end, strand; ``diff_types`` ``;``-joined, sorted, with repeats;
- per peptide and search: every split locus at the minimal split cost (≤ k);
- ``blocks`` are in the same unwrapped coordinates as ``start``/``end`` (``s1 = start``,
  ``e2 = end``, ``start`` in ``[0, length)``); ``frame`` is the label of the prefix block;
- the target search runs only for ``target_queries``; the decoy search for every searched peptide;
- reverse decoys are read in their own direction, so the reverse-decoy splits of ``Q`` are the
  target splits of ``reversed(Q)`` (same loci and blocks); shuffle decoys keep the target codons;
- "at the same locus the non-frameshift explanation wins": a split is not reported when the
  whole query aligns without a frameshift (M3-06 operations) at a cost ≤ the split's cost on an
  interval overlapping the split locus (same contig and strand).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))

from fs_genomes import background, codons_for, join_split, rc
from oracle_frameshift import TYPE_ORDER, Result, frameshift_oracle

from giap.config.model import TolerantConfig
from giap.mapping.decoys import DecoyType, build_decoys
from giap.mapping.frameshift import match_frameshift
from giap.mapping.hits import HIT_COLUMNS
from giap.mapping.sixframe import Frame, SixFrameIndex, build_sixframe
from giap.mapping.tolerant import match_tolerant
from giap.model import Replicon
from giap.rng import for_step

CFG = TolerantConfig()
FS1 = "GGCCGGAACCGACCATTGATAGAAAAACTGATGAACGCGTAA"  # spec example sequence (42 nt)
INT_COLS = ("frame_id", "aa_start", "aa_end", "start", "end", "frame", "n_diff")
KEY = ["contig", "start", "end", "strand", "frame", "n_diff", "diff_types", "blocks"]


def _cfg(**fs: object) -> TolerantConfig:
    return TolerantConfig(frameshift={"enabled": True, "min_side_aa": 4, **fs})  # type: ignore[arg-type]


def _index(genome: dict[str, str], topology: str = "linear") -> SixFrameIndex:
    reps = {c: Replicon(c, len(s), topology) for c, s in genome.items()}  # type: ignore[arg-type]
    return build_sixframe(genome, reps)


def _decoys(idx: SixFrameIndex, *kinds: DecoyType) -> dict[DecoyType, tuple[Frame, ...]]:
    return {k: build_decoys(idx, k, 50, for_step("sixframe_index", 0)) for k in kinds}


def _run(
    idx: SixFrameIndex,
    queries: list[str],
    *,
    targets: set[str] | None = None,
    kind: DecoyType = "reverse",
    cfg: TolerantConfig = CFG,
    decoys: dict[DecoyType, tuple[Frame, ...]] | None = None,
) -> pd.DataFrame:
    dec = _decoys(idx, kind) if decoys is None else decoys
    tq = set(queries) if targets is None else targets
    return match_frameshift(queries, idx, dec, dict.fromkeys(queries, kind), tq, cfg)


def _target(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["search"] == "target"].reset_index(drop=True)


def _decoy(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["search"] == "decoy"].reset_index(drop=True)


def _rows(df: pd.DataFrame) -> list[tuple[object, ...]]:
    return [tuple(r) for r in df[KEY].itertuples(index=False, name=None)]


def _expected(res: Result) -> list[tuple[object, ...]]:
    return sorted(
        (p.contig, p.start, p.end, p.strand, p.frame, p.cost, p.diff_types, p.blocks_str)
        for p in res.placements
    )


def _k(peptide: str) -> int:
    return 1 if len(peptide) <= 11 else 2


def _precedence(n_diff: object, types: object) -> tuple[int, tuple[int, ...]]:
    """(cost, diff types ranked substitution < swap < block < near_block < frameshift)."""
    return int(n_diff), tuple(sorted(TYPE_ORDER.index(t) for t in str(types).split(";")))  # type: ignore[call-overload]


def _assert_oracle(
    df: pd.DataFrame, peptide: str, idx: SixFrameIndex, min_side: int = 4
) -> list[tuple[object, ...]]:
    """Target rows of ``peptide`` agree with the brute-force oracle.

    One row per oracle locus interval; when several splits give the same interval (the junction
    can move by one codon when both frames read the same residue there), the row carries the
    lowest (cost, diff types by M3-05 precedence) and the blocks of one of the splits reaching it.
    """
    res = frameshift_oracle(peptide, idx, idx.frames, None, _k(peptide), min_side)
    groups: dict[tuple[object, ...], list[tuple[object, ...]]] = {}
    for row in _expected(res):
        groups.setdefault(row[:4], []).append(row)
    got = sorted(_rows(_target(df[df["peptide"] == peptide])))
    assert [r[:4] for r in got] == sorted(groups)
    for r in got:
        options = groups[r[:4]]
        best = min(_precedence(o[5], o[6]) for o in options)
        assert _precedence(r[5], r[6]) == best, (r, options)
        assert r in options, (r, options)
    return got


def _check_schema(df: pd.DataFrame) -> None:
    assert list(df.columns) == list(HIT_COLUMNS)
    for col in INT_COLS:
        assert df[col].dtype == np.int64, col


# --- examples 1-7 ------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def fs1() -> SixFrameIndex:
    return _index({"fs1": FS1})


def test_example_sequence_frames(fs1: SixFrameIndex) -> None:
    from giap.mapping.sixframe import decode

    assert decode(fs1.frames[2].aa) == "PEPTLDRKTDERV"  # +3 (I -> L)
    assert decode(fs1.frames[0].aa) == "GRNRPLLEKLMNA*"  # +1


def test_example_1_one_target_row(fs1: SixFrameIndex) -> None:
    df = _run(fs1, ["PEPTIDEKLMNA"])
    _check_schema(df)
    tgt = _target(df)
    assert len(tgt) == 1
    row = tgt.iloc[0]
    assert row["peptide"] == "PEPTIDEKLMNA"
    assert row["query"] == "PEPTLDEKLMNA"
    assert pd.isna(row["decoy_type"])
    assert (row["contig"], row["start"], row["end"], row["strand"]) == ("fs1", 2, 39, "+")
    assert row["frame"] == 3
    assert row["route"] == "tolerant"
    assert row["n_diff"] == 1
    assert row["diff_types"] == "frameshift"
    assert row["blocks"] == "2-20;21-39"
    _assert_oracle(df, "PEPTIDEKLMNA", fs1)


def test_example_2_min_side_7_gives_no_row(fs1: SixFrameIndex) -> None:
    df = _run(fs1, ["PEPTIDEKLMNA"], cfg=_cfg(min_side_aa=7))
    _check_schema(df)
    assert df.empty


def test_example_3_frameshift_plus_substitution(fs1: SixFrameIndex) -> None:
    tgt = _target(_run(fs1, ["PEPTVDEKLMNA"]))
    assert len(tgt) == 1
    row = tgt.iloc[0]
    assert (row["start"], row["end"], row["strand"], row["frame"]) == (2, 39, "+", 3)
    assert row["n_diff"] == 2
    assert row["diff_types"] == "frameshift;substitution"
    assert row["blocks"] == "2-20;21-39"


def test_example_4_ten_aa_frameshift_and_substitution_exceeds_k(fs1: SixFrameIndex) -> None:
    df = _run(fs1, ["PEPTVDEKLM"])  # 10 aa: k = 1; split (6 | 4) + substitution = cost 2
    assert _target(df).empty


def test_example_4_control_ten_aa_with_four_aa_side(fs1: SixFrameIndex) -> None:
    """Same 10-aa query without the substitution: one row; the suffix side is 4 aa (anchored)."""
    tgt = _target(_run(fs1, ["PEPTIDEKLM"]))
    assert _rows(tgt) == [("fs1", 2, 33, "+", 3, 1, "frameshift", "2-20;21-33")]


def test_four_aa_prefix_side(fs1: SixFrameIndex) -> None:
    tgt = _target(_run(fs1, ["PTIDEKLMNA"]))  # prefix PTID (4 aa) in +3, suffix EKLMNA in +1
    assert _rows(tgt) == [("fs1", 8, 39, "+", 3, 1, "frameshift", "8-20;21-39")]


def test_example_5_reverse_complement(fs1: SixFrameIndex) -> None:
    idx = _index({"fs1": rc(FS1)})
    df = _run(idx, ["PEPTIDEKLMNA"])
    tgt = _target(df)
    assert len(tgt) == 1
    row = tgt.iloc[0]
    assert (row["contig"], row["start"], row["end"], row["strand"]) == ("fs1", 3, 40, "-")
    assert row["blocks"] == "3-21;22-40"
    assert (row["n_diff"], row["diff_types"]) == (1, "frameshift")
    assert row["frame"] == -3  # label of the prefix block [22, 40) on the minus strand
    _assert_oracle(df, "PEPTIDEKLMNA", idx)


# Example 6: one nucleotide read twice (δ = -1). PEPTID ends with GAT; the suffix S K L M N A
# starts with TCT, whose T is the last base of GAT.
LEFT6, RIGHT6 = background(30, 6), background(30, 66)
PLANT6 = "CCAGAACCAACAATTGAT" + "CT" + "AAACTGATGAACGCG"


def test_example_6_delta_minus_one_overlapping_blocks() -> None:
    idx = _index({"g6": LEFT6 + PLANT6 + RIGHT6})
    df = _run(idx, ["PEPTIDSKLMNA"])
    tgt = _target(df)
    assert _rows(tgt) == [("g6", 30, 65, "+", 1, 1, "frameshift", "30-48;47-65")]
    s1e1, s2e2 = tgt.iloc[0]["blocks"].split(";")
    e1, s2 = int(s1e1.split("-")[1]), int(s2e2.split("-")[0])
    assert s2 == e1 - 1  # the second block starts 1 nt before the first ends
    _assert_oracle(df, "PEPTIDSKLMNA", idx)


def test_example_6_minus_strand() -> None:
    n = len(LEFT6 + PLANT6 + RIGHT6)
    idx = _index({"g6": rc(LEFT6 + PLANT6 + RIGHT6)})
    df = _run(idx, ["PEPTIDSKLMNA"])
    assert _rows(_target(df)) == [
        ("g6", n - 65, n - 30, "-", -1, 1, "frameshift", f"{n - 65}-{n - 47};{n - 48}-{n - 30}")
    ]


def test_example_7_disabled_gives_empty_frame(fs1: SixFrameIndex) -> None:
    df = _run(fs1, ["PEPTIDEKLMNA", "PEPTVDEKLMNA"], cfg=_cfg(enabled=False))
    _check_schema(df)
    assert df.empty


def test_empty_queries_give_empty_frame(fs1: SixFrameIndex) -> None:
    df = match_frameshift([], fs1, _decoys(fs1, "reverse"), {}, set(), CFG)
    _check_schema(df)
    assert df.empty


# --- edge cases: stops, NOMATCH, side lengths ---------------------------------------------------

LEFT, RIGHT = background(30, 11), background(30, 12)
PEPTID = "CCAGAACCAACAATTGAT"  # P E P T I D
KLMNA = "AAACTGATGAACGCG"


@pytest.mark.parametrize(
    ("first_suffix_codon", "query", "expected"),
    [
        # sense control: W (TGG) is the first suffix codon -> split at p = 6, cost 1
        ("TGG", "PEPTIDWKLMNA", [("g", 30, 67, "+", 1, 1, "frameshift", "30-48;49-67")]),
        # a stop there: the split may not consume it; the cheapest explanation moves the split one
        # residue right (W vs ATA = I/L in the prefix frame: a substitution) -> cost 2; the split
        # may also sit one codon later (AAA = K in both frames): blocks 30-54;55-67
        (
            "TAA",
            "PEPTIDWKLMNA",
            [("g", 30, 67, "+", 1, 2, "frameshift;substitution", "30-51;52-67")],
        ),
        ("TGG", "PEPTIDWKLM", [("g", 30, 61, "+", 1, 1, "frameshift", "30-48;49-61")]),
        ("TAA", "PEPTIDWKLM", []),  # k = 1: the stop forbids the cost-1 split
    ],
)
def test_split_points_on_a_stop_are_rejected(
    first_suffix_codon: str, query: str, expected: list[tuple[object, ...]]
) -> None:
    idx = _index({"g": LEFT + PEPTID + "A" + first_suffix_codon + KLMNA + RIGHT})
    df = _run(idx, [query])
    got = _rows(_target(df))
    if first_suffix_codon == "TAA" and expected:
        assert [r[:7] for r in got] == [e[:7] for e in expected]
        assert got[0][7] in ("30-51;52-67", "30-54;55-67")  # never 30-48;49-67 (on the stop)
    else:
        assert got == expected
    _assert_oracle(df, query, idx)


def test_nomatch_query_residue_costs_a_substitution(fs1: SixFrameIndex) -> None:
    tgt = _target(_run(fs1, ["PEPTBDEKLMNA"]))  # B is outside the residue alphabet
    assert _rows(tgt) == [("fs1", 2, 39, "+", 3, 2, "frameshift;substitution", "2-20;21-39")]


# 12-aa queries split 8 | 4 or 4 | 8 with the one remaining difference on the 4-aa side
# (k = 2: frameshift + one difference); the genome carries the true sides, the query the change.
SHORT_SIDE = [
    # (genome prefix, genome suffix, query prefix, query suffix, diff type)
    ("MKWHPEPT", "LMNA", "MKWHPEPT", "LWNA", "substitution"),
    ("MKWHPEPT", "LMNA", "MKWHPEPT", "MLNA", "swap"),
    ("MKWHPEPT", "GGKAC", "MKWHPEPT", "NKAC", "block"),
    ("MKWHPEPT", "ADKHC", "MKWHPEPT", "WKHC", "near_block"),
    ("WHPE", "TLMNAKDF", "WHPF", "TLMNAKDF", "substitution"),
    ("WHPE", "TLMNAKDF", "WHEP", "TLMNAKDF", "swap"),
    ("WGGPE", "TLMNAKDF", "WNPE", "TLMNAKDF", "block"),
    ("WGVPE", "TLMNAKDF", "WRPE", "TLMNAKDF", "near_block"),
]


@pytest.mark.parametrize("delta", [1, -1])
@pytest.mark.parametrize(("gp", "gs", "qp", "qs", "dtype"), SHORT_SIDE)
def test_four_aa_side_with_one_difference(
    gp: str, gs: str, qp: str, qs: str, dtype: str, delta: int
) -> None:
    nt = join_split(codons_for(gp), codons_for(gs), delta)
    assert nt is not None
    idx = _index({"g": LEFT + nt + RIGHT})
    query = qp + qs
    df = _run(idx, [query])
    got = _assert_oracle(df, query, idx)
    assert len(got) == 1
    _, start, end, strand, _, n_diff, types, blocks = got[0]
    assert (start, strand, n_diff, types) == (30, "+", 2, ";".join(sorted(["frameshift", dtype])))
    assert end == 30 + len(nt)
    (_s1, e1), (s2, _e2) = (tuple(map(int, b.split("-"))) for b in str(blocks).split(";"))
    assert s2 - e1 == delta


def test_min_side_counts_query_residues(fs1: SixFrameIndex) -> None:
    """``min_side_aa: 5`` keeps the 6 | 6 split of example 1 and drops the 6 | 4 split."""
    cfg = _cfg(min_side_aa=5)
    assert len(_target(_run(fs1, ["PEPTIDEKLMNA"], cfg=cfg))) == 1
    assert _target(_run(fs1, ["PEPTIDEKLM"], cfg=cfg)).empty
    assert _target(_run(fs1, ["PTIDEKLMNA"], cfg=cfg)).empty


def test_short_queries_are_not_searched(fs1: SixFrameIndex) -> None:
    """Below ``min_len`` (10) nothing is searched (and no decoy type is needed)."""
    df = match_frameshift(["PEPTIDEKL"], fs1, _decoys(fs1, "reverse"), {}, {"PEPTIDEKL"}, CFG)
    assert df.empty


def test_lowered_min_len_searches_shorter_queries(fs1: SixFrameIndex) -> None:
    cfg = TolerantConfig(min_len=8, max_diff={"8-11": 1, "12-40": 2})
    tgt = _target(_run(fs1, ["PTIDEKLM"], cfg=cfg))  # 4 | 4
    assert _rows(tgt) == [("fs1", 8, 33, "+", 3, 1, "frameshift", "8-20;21-33")]


def test_queries_of_at_most_four_aa_give_no_rows_and_no_error(fs1: SixFrameIndex) -> None:
    """M3-07 lead decision 2, carried over: even with ``min_len`` and ``min_side_aa`` lowered."""
    cfg = TolerantConfig(
        min_len=1, max_diff={"1-11": 1, "12-40": 2}, frameshift={"enabled": True, "min_side_aa": 2}
    )
    df = _run(fs1, ["IDEK", "DEK"], cfg=cfg)
    _check_schema(df)
    assert df.empty


def test_two_aa_side_with_lowered_min_side() -> None:
    """Sides shorter than 4 aa (config ``min_side_aa: 2``) use the anchored check too."""
    cfg = _cfg(min_side_aa=2)
    nt = join_split(codons_for("MKWHPEPT"), codons_for("QR"), 1)
    assert nt is not None
    idx = _index({"g": LEFT + nt + RIGHT})
    q = "MKWHPEPTQR"  # 8 | 2, k = 1
    df = _run(idx, [q], cfg=cfg)
    got = _assert_oracle(df, q, idx, min_side=2)
    assert got == [("g", 30, 61, "+", 1, 1, "frameshift", "30-54;55-61")]
    assert _target(_run(idx, [q])).empty  # default min_side_aa: 4


def test_k_above_two_is_rejected(fs1: SixFrameIndex) -> None:
    cfg = TolerantConfig(max_diff={"10-11": 1, "12-40": 3})
    with pytest.raises(ValueError):
        _run(fs1, ["PEPTIDEKLMNA"], cfg=cfg)


def test_decoy_type_without_frames_raises(fs1: SixFrameIndex) -> None:
    with pytest.raises(ValueError):
        match_frameshift(
            ["PEPTIDEKLMNA"],
            fs1,
            _decoys(fs1, "reverse"),
            {"PEPTIDEKLMNA": "shuffle"},
            {"PEPTIDEKLMNA"},
            CFG,
        )


def test_searched_peptide_without_decoy_type_raises(fs1: SixFrameIndex) -> None:
    with pytest.raises(ValueError):
        match_frameshift(["PEPTIDEKLMNA"], fs1, _decoys(fs1, "reverse"), {}, {"PEPTIDEKLMNA"}, CFG)


# --- circular replicons: splits across the origin ----------------------------------------------


def _origin_case(strand: str, delta: int, cut: int, peptide: str) -> tuple[SixFrameIndex, int]:
    """Circular replicon of 300 nt whose origin falls ``cut`` nt into the planted split."""
    p = len(peptide) // 2
    nt = join_split(codons_for(peptide[:p]), codons_for(peptide[p:]), delta)
    assert nt is not None
    body = background(150, 40 + cut) + nt
    body += background(300 - len(body), 41 + cut)
    plus = body[150 + cut :] + body[: 150 + cut]  # origin `cut` nt into the planted split
    seq = plus if strand == "+" else rc(plus)
    return _index({"c": seq}, "circular"), len(nt)


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("delta", [1, -1])
@pytest.mark.parametrize("cut", [1, 2, 3, 17, 18, 19, 20, 33])
def test_split_across_the_origin(strand: str, delta: int, cut: int) -> None:
    peptide = "MKWHPEPTDFLMNAQR"  # 16 aa, 8 | 8
    idx, span = _origin_case(strand, delta, cut, peptide)
    df = _run(idx, [peptide])
    got = _assert_oracle(df, peptide, idx)
    assert len(got) == 1
    _, start, end, st, _, n_diff, types, _ = got[0]
    assert (st, n_diff, types) == (strand, 1, "frameshift")
    assert end - start == span
    assert 0 <= start < 300 < end  # unwrapped across the origin (§6.5)


@pytest.mark.parametrize("back", [0, 1, 2, 5, 12])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_longest_split_across_the_origin_fits_the_wrap(back: int, strand: str) -> None:
    """40 aa, k = 2: one GG -> N block in the prefix (41 genome residues) + 1 nt = 124 nt."""
    gp = "MKWHPEPTDFLMGGAQRSTVWYCH"  # 24 aa
    gs = "EKLMNAQRSTVWYPDFK"  # 17 aa
    peptide = gp.replace("GG", "N") + gs  # 23 + 17 = 40 aa
    nt = join_split(codons_for(gp), codons_for(gs), 1)
    assert nt is not None and len(nt) == 124 and len(peptide) == 40
    body = background(200, 7) + nt
    body += background(400 - len(body), 8)
    r = 200 + back  # the split starts `back` nt before the origin of the rotated replicon
    plus = body[r:] + body[:r]
    seq = plus if strand == "+" else rc(plus)
    idx = _index({"c": seq}, "circular")
    df = _run(idx, [peptide])
    got = _assert_oracle(df, peptide, idx)
    assert len(got) == 1
    _, start, end, _, _, n_diff, types, _ = got[0]
    assert end - start == 124 and (n_diff, types) == (2, "block;frameshift")
    if back:
        assert 0 <= start < 400 < end


# --- decoys: searched identically, each peptide against its own decoy type ---------------------


def _planted(seed: int = 3) -> tuple[SixFrameIndex, str]:
    gp, gs = "MKWHPEPTDF", "EKLMNAQRCH"
    nt = join_split(codons_for(gp), codons_for(gs), 1)
    assert nt is not None
    genome = background(240, seed) + nt + background(240, seed + 1)
    return _index({"g": genome}), gp + gs


def test_shuffle_decoys_keep_target_codons() -> None:
    """With the target frames passed as the 'shuffle' decoy frames, the decoy search finds the
    same split rows as the target search (shuffle mirror = identity)."""
    idx, pep = _planted()
    df = _run(idx, [pep], kind="shuffle", decoys={"shuffle": idx.frames})
    tgt, dec = _target(df), _decoy(df)
    assert len(tgt) == 1
    assert _rows(dec) == _rows(tgt)
    assert set(dec["decoy_type"]) == {"shuffle"}
    assert (dec["route"] == "tolerant").all()


def test_reverse_decoy_split_of_reversed_peptide_is_the_target_split() -> None:
    """Reverse decoys are read in their own direction: the reverse-decoy splits of reversed(P)
    sit on the target split of P (same locus, blocks, cost and types)."""
    idx, pep = _planted()
    rev = pep[::-1]
    df = match_frameshift(
        [pep, rev], idx, _decoys(idx, "reverse"), {pep: "reverse", rev: "reverse"}, {pep}, CFG
    )
    tgt = _target(df[df["peptide"] == pep])
    dec = _decoy(df[df["peptide"] == rev])
    cols = ["contig", "start", "end", "strand", "n_diff", "diff_types", "blocks"]
    assert len(tgt) == 1
    assert dec[cols].values.tolist() == tgt[cols].values.tolist()
    res = frameshift_oracle(rev, idx, _decoys(idx, "reverse")["reverse"], "reverse", 2, 4)
    assert sorted(_rows(dec)) == _expected(res)
    assert _decoy(df[df["peptide"] == pep]).empty  # reversed(P) is not planted


def test_reverse_decoy_example_on_fs1(fs1: SixFrameIndex) -> None:
    rev = "PEPTIDEKLMNA"[::-1]
    df = _run(fs1, [rev], targets=set())
    dec = _decoy(df)
    assert len(dec) == 1
    row = dec.iloc[0]
    assert (row["start"], row["end"], row["strand"], row["blocks"]) == (2, 39, "+", "2-20;21-39")
    assert (row["n_diff"], row["diff_types"], row["decoy_type"]) == (1, "frameshift", "reverse")
    assert row["frame"] == 1  # the prefix of the decoy split mirrors the block [21, 39)


def test_each_peptide_searches_its_own_decoy_type() -> None:
    idx, pep = _planted()
    rev = pep[::-1]
    decoys: dict[DecoyType, tuple[Frame, ...]] = {
        "reverse": _decoys(idx, "reverse")["reverse"],
        "shuffle": idx.frames,  # identity "shuffle" so that a planted split is findable
    }
    # pep against shuffle (finds the planted split), rev against reverse (finds it as well)
    df = match_frameshift([pep, rev], idx, decoys, {pep: "shuffle", rev: "reverse"}, set(), CFG)
    assert _target(df).empty
    dec = _decoy(df)
    assert dict(zip(dec["peptide"], dec["decoy_type"], strict=True)) == {
        pep: "shuffle",
        rev: "reverse",
    }
    # swapped types: neither peptide's own decoy frames hold a split
    df2 = match_frameshift([pep, rev], idx, decoys, {pep: "reverse", rev: "shuffle"}, set(), CFG)
    assert df2.empty


def test_il_twins_keep_their_own_decoy_type_and_target_membership() -> None:
    idx, pep = _planted()
    twin = pep.replace("L", "I")
    assert twin != pep
    decoys: dict[DecoyType, tuple[Frame, ...]] = {
        "reverse": _decoys(idx, "reverse")["reverse"],
        "shuffle": idx.frames,
    }
    df = match_frameshift([pep, twin], idx, decoys, {pep: "shuffle", twin: "reverse"}, {twin}, CFG)
    assert set(_target(df)["peptide"]) == {twin}
    assert set(_decoy(df)["peptide"]) == {pep}
    assert set(df["query"]) == {pep}  # both normalise to the L form


def test_target_search_only_for_target_queries_decoys_for_all() -> None:
    """M3-10 Q1: a peptide with an exact hit is no leftover; its decoy split search still runs."""
    idx, pep = _planted()
    df = _run(idx, [pep], targets=set(), kind="shuffle", decoys={"shuffle": idx.frames})
    assert _target(df).empty
    assert len(_decoy(df)) == 1


# --- interaction with ordinary tolerant hits -----------------------------------------------------


def test_equal_cost_tolerant_and_split_at_different_loci_are_both_kept() -> None:
    """Edge case: a 1-substitution in-frame copy and a cost-1 split elsewhere: both loci stay —
    the tolerant driver reports the first, the split matcher the second; no row is shared."""
    gp, gs = "MKWHPEPTDF", "EKLMNAQRCH"
    pep = gp + gs
    split = join_split(codons_for(gp), codons_for(gs), 1)
    assert split is not None
    inframe = "".join(codons_for(pep[:5] + "W" + pep[6:]))  # one substitution
    genome = background(150, 21) + split + background(150, 22) + inframe + background(150, 23)
    idx = _index({"g": genome})
    decoys = _decoys(idx, "reverse")
    args = ([pep], idx, decoys, {pep: "reverse"}, {pep}, CFG)
    tol = _target(match_tolerant(*args))
    fs = _target(match_frameshift(*args))
    assert tol[["start", "n_diff", "diff_types"]].values.tolist() == [
        [150 + len(split) + 150, 1, "substitution"]
    ]
    assert fs[["start", "end", "n_diff", "diff_types"]].values.tolist() == [
        [150, 150 + len(split), 1, "frameshift"]
    ]
    both = pd.concat([match_tolerant(*args), match_frameshift(*args)])
    key = ["peptide", "search", "decoy_type", "contig", "start", "end", "strand"]
    assert not both.astype({"decoy_type": str}).duplicated(subset=key).any()


def test_same_locus_non_frameshift_explanation_wins() -> None:
    """Edge case: at one locus a 1-substitution in-frame alignment and a cost-1 split explain the
    query equally well; the split is not reported (precedence, M3-05: frameshift last).

    Genome: prefix codons, then G + A*12. In the prefix frame: GAA AAA AAA AAA = E K K K; one
    nucleotide later (δ = +1): AAA AAA AAA AAA = K K K K. Query = prefix + KKKK.
    """
    gp = "MWHPEPTF"
    tail = "G" + "A" * 12 + "CC"
    genome = background(60, 31) + "".join(codons_for(gp)) + tail + background(60, 32)
    idx = _index({"g": genome})
    pep = gp + "KKKK"
    decoys = _decoys(idx, "reverse")
    args = ([pep], idx, decoys, {pep: "reverse"}, {pep}, CFG)
    tol = _target(match_tolerant(*args))
    assert tol[["start", "n_diff", "diff_types"]].values.tolist() == [[60, 1, "substitution"]]
    fs = _target(match_frameshift(*args))
    assert not ((fs["start"] < 60 + 36) & (fs["end"] > 60)).any()
    _assert_oracle(match_frameshift(*args), pep, idx)


def test_rows_never_carry_non_frameshift_explanations() -> None:
    idx, pep = _planted()
    df = _run(idx, [pep, pep[:5] + "W" + pep[6:]])
    assert len(df) >= 2
    assert all("frameshift" in t.split(";") for t in df["diff_types"])
    assert all(len(str(b).split(";")) == 2 for b in df["blocks"])
    assert (df["n_diff"] == [len(t.split(";")) for t in df["diff_types"]]).all()


# --- order and determinism ----------------------------------------------------------------------


def _multi() -> tuple[SixFrameIndex, list[str]]:
    contigs: dict[str, str] = {}
    peps: list[str] = []
    for i, (gp, gs, delta) in enumerate(
        [
            ("MKWHPEPTDF", "EKLMNAQRCH", 1),
            ("PEPTIDSAKR", "WWHHKLMN", -1),
            ("GHKLMNPQ", "RSTVWYAC", 1),
        ]
    ):
        nt = join_split(codons_for(gp), codons_for(gs), delta)
        assert nt is not None
        body = background(90, 50 + i) + nt + background(90, 60 + i)
        contigs[f"c{i}"] = body if i != 1 else rc(body)
        peps.append(gp + gs)
    return _index(contigs), peps


def test_rows_are_sorted_and_input_order_does_not_matter() -> None:
    idx, peps = _multi()
    queries = [*peps, peps[0][:4] + "Y" + peps[0][5:], peps[2][::-1]]
    decoys = _decoys(idx, "reverse")
    decoy_of = dict.fromkeys(queries, "reverse")
    a = match_frameshift(queries, idx, decoys, decoy_of, set(queries), CFG)  # type: ignore[arg-type]
    b = match_frameshift(queries[::-1], idx, decoys, decoy_of, set(queries), CFG)  # type: ignore[arg-type]
    c = match_frameshift([*queries, *queries], idx, decoys, decoy_of, set(queries), CFG)  # type: ignore[arg-type]
    pd.testing.assert_frame_equal(a, b)
    pd.testing.assert_frame_equal(a, c)
    assert set(_target(a)["peptide"]) >= set(peps)
    order = {"c0": 0, "c1": 1, "c2": 2}
    keys = [
        (
            r.peptide,
            r.search != "target",
            "" if pd.isna(r.decoy_type) else str(r.decoy_type),
            order[r.contig],
            r.start,
            r.end,
            r.strand,
        )
        for r in a.itertuples(index=False)
    ]
    assert keys == sorted(keys)


def test_repeated_calls_are_identical() -> None:
    idx, peps = _multi()
    decoys = _decoys(idx, "shuffle")
    decoy_of = dict.fromkeys(peps, "shuffle")
    a = match_frameshift(peps, idx, decoys, decoy_of, set(peps), CFG)  # type: ignore[arg-type]
    b = match_frameshift(peps, idx, decoys, decoy_of, set(peps), CFG)  # type: ignore[arg-type]
    pd.testing.assert_frame_equal(a, b)
    for pep in peps:
        _assert_oracle(a, pep, idx)
