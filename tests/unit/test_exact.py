"""Route 1: exact Aho-Corasick peptide matching (issue M3-04; GIPPA_SPEC.md §7.2.2, §7.2.9).

Genomes ``G`` (linear ``c1``) and ``R`` (circular ``r``, ``wrap_nt=12``) are those of M3-02.
``match_exact`` takes the length window as keyword arguments ``min_len``/``max_len`` (defaults 7
and 40, config ``mapping.exact``); the spec's examples run with ``min_len=4`` (or 3).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gippa.io.fasta import read_genome_fasta
from gippa.mapping.decoys import DecoyType, build_decoys, mirror_aa
from gippa.mapping.exact import match_exact
from gippa.mapping.hits import HIT_COLUMNS
from gippa.mapping.sixframe import SixFrameIndex, build_sixframe, decode
from gippa.model import Interval, Replicon, to_one_based_closed
from gippa.model.genes import gene_model_from_bakta_json
from gippa.rng import for_step
from gippasim.codons import CODE_11
from gippasim.simulate import SimConfig, simulate

G = "ATGGCTAAAATTCTGTGGTAACGTNNA"  # frames: MAKLLW*RX WLKFCGNX G*NSVVTX XTLPQNFSH XRYHRLLA XVTTEF*P
R = "AAAGAACTGTAAGGCTCCATGCCG"  # circular, wrap 12: +0 = KEL*GSMPKEL*
_RC = str.maketrans("ACGTN", "TGCAN")
INT_COLS = ("frame_id", "aa_start", "aa_end", "start", "end", "frame", "n_diff")
SPEC_COLUMNS = (
    "peptide",
    "query",
    "search",
    "decoy_type",
    "contig",
    "frame_id",
    "aa_start",
    "aa_end",
    "start",
    "end",
    "strand",
    "frame",
    "route",
    "n_diff",
    "diff_types",
    "blocks",
)
SENSE = sorted(c for c, a in CODE_11.items() if a != "*")
BACK = {a: next(c for c in SENSE if CODE_11[c] == a) for a in set(CODE_11.values()) - {"*"}}


def _rc(s: str) -> str:
    return s.translate(_RC)[::-1]


def _index(genome: dict[str, str], topology: str = "linear", wrap_nt: int = 120) -> SixFrameIndex:
    reps = {c: Replicon(c, len(s), topology) for c, s in genome.items()}  # type: ignore[arg-type]
    return build_sixframe(genome, reps, wrap_nt=wrap_nt)


def _match(
    idx: SixFrameIndex,
    queries: list[str],
    kind: DecoyType = "reverse",
    **kw: int,
) -> pd.DataFrame:
    decoys = {kind: build_decoys(idx, kind, 50, for_step("sixframe_index", 0))}
    return match_exact(queries, idx, decoys, {q: kind for q in queries}, **kw)


def _target(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["search"] == "target"].reset_index(drop=True)


def _decoy(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["search"] == "decoy"].reset_index(drop=True)


def _ivs(df: pd.DataFrame) -> list[Interval]:
    return [
        Interval(r.contig, int(r.start), int(r.end), r.strand) for r in df.itertuples(index=False)
    ]


def _nt(peptide: str) -> str:
    return "".join(BACK[a] for a in peptide)


def _background(n: int, seed: int) -> str:
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list("ACGT"), size=n))


# --- schema ------------------------------------------------------------------------------------


def test_hit_columns_match_spec() -> None:
    assert tuple(HIT_COLUMNS) == SPEC_COLUMNS


def test_empty_query_list_returns_typed_empty_frame() -> None:
    idx = _index({"c1": G})
    df = _match(idx, [])
    assert list(df.columns) == list(HIT_COLUMNS)
    assert len(df) == 0
    full = _match(idx, ["MAKLLW"], min_len=4)
    for col in INT_COLS:
        assert pd.api.types.is_integer_dtype(df[col]), col
        assert df[col].dtype == full[col].dtype, col


# --- examples ----------------------------------------------------------------------------------


def test_example_1_single_target_locus() -> None:
    idx = _index({"c1": G})
    df = _target(_match(idx, ["MAKILW"], min_len=4))
    assert len(df) == 1
    row = df.iloc[0]
    assert _ivs(df) == [Interval("c1", 0, 18, "+")]
    assert row["peptide"] == "MAKILW"
    assert row["query"] == "MAKLLW"
    assert int(row["frame"]) == 1
    assert int(row["frame_id"]) == 0
    assert (int(row["aa_start"]), int(row["aa_end"])) == (0, 6)
    assert row["route"] == "exact"
    assert int(row["n_diff"]) == 0
    assert row["diff_types"] == ""
    assert pd.isna(row["blocks"])
    assert pd.isna(row["decoy_type"])


def test_example_2_i_equals_l() -> None:
    idx = _index({"c1": G})
    df = _target(_match(idx, ["MAKLLW"], min_len=4))
    assert _ivs(df) == [Interval("c1", 0, 18, "+")]
    assert df.loc[0, "query"] == "MAKLLW"
    assert df.loc[0, "peptide"] == "MAKLLW"


def test_duplicate_queries_share_a_key_but_each_peptide_gets_rows() -> None:
    pep = "PEPTLDEK"
    seq = _background(60, 1) + _nt(pep) + _background(60, 2)
    idx = _index({"c1": seq})
    df = _target(_match(idx, ["PEPTIDEK", "PEPTLDEK"]))
    assert sorted(df["peptide"]) == ["PEPTIDEK", "PEPTLDEK"]
    assert set(df["query"]) == {"PEPTLDEK"}
    assert set(_ivs(df)) == {Interval("c1", 60, 84, "+")}


@pytest.fixture(scope="module")
def sim50(tmp_path_factory: pytest.TempPathFactory) -> tuple[SixFrameIndex, object]:
    res = simulate(SimConfig(genome_length=50_000), seed=11, outdir=tmp_path_factory.mktemp("s"))
    (rec,) = read_genome_fasta(res.genome_fasta)
    reps = {rec.id: Replicon(rec.id, len(rec.sequence), "circular")}
    idx = build_sixframe({rec.id: rec.sequence}, reps)
    return idx, gene_model_from_bakta_json(res.bakta_json, reps)


def test_example_3_peptide_from_sim_gene_maps_inside_gene(sim50: tuple) -> None:  # type: ignore[type-arg]
    idx, model = sim50
    gene = model.by_locus("GSIM_00005")
    pep = gene.aa[10:25]
    df = _target(_match(idx, [pep]))
    rep = idx.replicons[gene.interval.contig]
    g = gene.interval
    inside = [
        iv
        for iv in _ivs(df)
        if iv.strand == g.strand
        and iv.contig == g.contig
        and (iv.start - g.start) % rep.length + iv.length_bp <= g.length_bp
    ]
    assert inside
    assert all(idx.frame_label(iv) == idx.frame_label(g) for iv in inside)
    labels = {int(f) for f, iv in zip(df["frame"], _ivs(df), strict=True) if iv in inside}
    assert labels == {idx.frame_label(g)}


def test_example_4_planted_duplicate_keeps_every_locus() -> None:
    pep = "WMHKCYFENQ"  # 10 aa, 30 nt
    nt = _nt(pep)
    seq = _background(300, 3) + nt + _background(900, 4) + nt + _background(300, 5)
    idx = _index({"c1": seq})
    df = _target(_match(idx, [pep]))
    assert _ivs(df) == [Interval("c1", 300, 330, "+"), Interval("c1", 1230, 1260, "+")]


def test_example_5_origin_crossing_hit_is_one_row() -> None:
    idx = _index({"r": R}, "circular", wrap_nt=12)
    df = _target(_match(idx, ["MPKEL"], min_len=4))
    assert _ivs(df) == [Interval("r", 18, 33, "+")]
    assert int(df.loc[0, "frame"]) == 1
    assert to_one_based_closed(_ivs(df)[0]) == (19, 33)  # TSV coordinates


def test_origin_copies_are_deduplicated() -> None:
    idx = _index({"r": R}, "circular", wrap_nt=12)
    df = _target(_match(idx, ["KEL"], min_len=3))
    assert _ivs(df) == [Interval("r", 0, 9, "+")]


def test_example_6_linear_contig_never_wraps() -> None:
    idx = _index({"r": R}, "linear")
    df = _target(_match(idx, ["MPKEL"], min_len=4))
    assert len(df) == 0


def test_example_7_nomatch_letters_never_match() -> None:
    seq = _background(30, 6) + _nt("SECKLLAAR") + _background(30, 7)
    idx = _index({"c1": seq})
    assert len(_match(idx, ["SEUKLLAAR"])) == 0
    assert len(_target(_match(idx, ["SECKLLAAR"]))) == 1


def test_query_x_never_matches_target_x() -> None:
    seq = "ATGGCTAAAATTCTGTGGNNNAAA"  # +0 = MAKLLWXK
    idx = _index({"c1": seq})
    assert decode(idx.frames[0].aa) == "MAKLLWXK"
    assert len(_match(idx, ["MAKLLWXK"])) == 0
    assert len(_target(_match(idx, ["MAKLLW"], min_len=4))) == 1


def test_example_8_queries_never_match_across_stop() -> None:
    idx = _index({"c1": G})  # +0 = MAKLLW*RX
    assert len(_match(idx, ["LLWR"], min_len=4)) == 0
    assert len(_match(idx, ["KLLWR"], min_len=4)) == 0


def test_example_9_reverse_decoy_hit_has_mirrored_interval() -> None:
    idx = _index({"c1": G})  # reversed +0 = XR*WLLKAM
    df = _match(idx, ["WLLKAM", "WILKAM"], min_len=4)
    assert len(_target(df)) == 0
    dec = _decoy(df)
    assert sorted(dec["peptide"]) == ["WILKAM", "WLLKAM"]
    for r in dec.itertuples(index=False):
        assert r.decoy_type == "reverse"
        assert (r.contig, int(r.frame_id), int(r.aa_start), int(r.aa_end)) == ("c1", 0, 3, 9)
        assert (int(r.start), int(r.end), r.strand, int(r.frame)) == (0, 18, "+", 1)
        assert mirror_aa("reverse", 9, 3, 9) == (0, 6)
        assert r.route == "exact" and int(r.n_diff) == 0


def test_each_peptide_scans_only_its_own_decoy_type() -> None:
    idx = _index({"c1": G})
    rev = build_decoys(idx, "reverse", 50, for_step("sixframe_index", 0))
    shu = build_decoys(idx, "shuffle", 50, for_step("sixframe_index", 0))
    assert all("WLLKAM" not in decode(f.aa) for f in shu)  # precondition
    decoys = {"reverse": rev, "shuffle": shu}
    df = match_exact(["WLLKAM"], idx, decoys, {"WLLKAM": "shuffle"}, min_len=4)
    assert len(df) == 0
    df = match_exact(["WLLKAM"], idx, decoys, {"WLLKAM": "reverse"}, min_len=4)
    assert list(df["decoy_type"]) == ["reverse"]


def test_shuffle_decoy_rows_mirror_the_same_window() -> None:
    seq = _background(600, 8)
    idx = _index({"c1": seq})
    shu = build_decoys(idx, "shuffle", 50, for_step("sixframe_index", 3))
    fid = 1
    aa = decode(shu[fid].aa)
    start = next(i for i in range(len(aa) - 8) if "*" not in aa[i : i + 8])
    pep = aa[start : start + 8]
    df = match_exact([pep], idx, {"shuffle": shu}, {pep: "shuffle"})
    want = idx.to_interval(fid, start, start + 8)
    rows = _decoy(df)
    assert want in _ivs(rows)
    assert set(rows["decoy_type"]) == {"shuffle"}


# --- edge cases --------------------------------------------------------------------------------


def test_overlapping_occurrences_each_start_is_a_locus() -> None:
    idx = _index({"c1": "GCT" * 10})  # +0 = AAAAAAAAAA
    df = _target(_match(idx, ["AAAAAAAA"]))
    plus = [iv for iv in _ivs(df) if iv.strand == "+"]
    assert plus == [Interval("c1", s, s + 24, "+") for s in (0, 3, 6)]


def test_peptide_on_both_strands_gives_two_loci() -> None:
    pep = "WMHKCYFENQ"
    nt = _nt(pep)
    seq = _background(90, 9) + nt + _background(150, 10) + _rc(nt) + _background(90, 11)
    idx = _index({"c1": seq})
    df = _target(_match(idx, [pep]))
    assert _ivs(df) == [Interval("c1", 90, 120, "+"), Interval("c1", 270, 300, "-")]
    assert [int(f) for f in df["frame"]] == [1, idx.frame_label(Interval("c1", 270, 300, "-"))]


def test_minus_strand_origin_crossing_hit() -> None:
    pep = "WMHKCYFENQ"
    nt = _rc(_nt(pep))  # 30 nt on the minus strand
    body = _background(240, 12)
    seq = nt[15:] + body + nt[:15]  # the minus-strand locus spans the origin
    idx = _index({"c": seq}, "circular")
    df = _target(_match(idx, [pep]))
    L = len(seq)
    assert _ivs(df) == [Interval("c", L - 15, L + 15, "-")]


def test_length_window_excludes_short_and_long_queries() -> None:
    long = "WMHKCYFENQ" * 5  # 50 aa
    seq = _background(30, 13) + _nt(long) + _background(30, 14)
    idx = _index({"c1": seq})
    assert len(_match(idx, ["WMHKCY"])) == 0  # 6 aa < min_len 7
    assert len(_target(_match(idx, ["WMHKCY"], min_len=4))) >= 1
    assert len(_match(idx, [long[:41]])) == 0  # 41 aa > max_len 40
    assert len(_target(_match(idx, [long[:40]]))) == 2  # periodic: offsets 0 and 10


def test_multi_contig_rows_follow_contig_order() -> None:
    pep = "WMHKCYFENQ"
    nt = _nt(pep)
    genome = {"z": _background(60, 15) + nt, "a": nt + _background(60, 16)}
    idx = _index(genome)
    df = _target(_match(idx, [pep]))
    assert _ivs(df) == [Interval("z", 60, 90, "+"), Interval("a", 0, 30, "+")]
    assert [int(f) for f in df["frame_id"]] == [0, 6]


def test_rows_sorted_and_deterministic() -> None:
    seq = _background(3000, 17)
    idx = _index({"c1": seq, "c2": _background(1500, 18)})
    frames = [decode(f.aa) for f in idx.frames]
    peps = sorted(
        {
            fr[i : i + 7]
            for fr in frames
            for i in range(0, len(fr) - 7, 37)
            if "*" not in fr[i : i + 7] and "X" not in fr[i : i + 7]
        }
    )[:40]
    a = _match(idx, peps)
    b = _match(idx, list(reversed(peps)))
    pd.testing.assert_frame_equal(a.reset_index(drop=True), b.reset_index(drop=True))
    order = {c: i for i, c in enumerate(idx.replicons)}
    assert list(a["peptide"]) == sorted(a["peptide"])
    for _, grp in a.groupby(["peptide", "search"], sort=False):
        keys = [(order[r.contig], r.start, r.end, r.strand) for r in grp.itertuples(index=False)]
        assert keys == sorted(keys)
        assert len(set(keys)) == len(keys)
    blocks = [(p, s) for p, s in zip(a["peptide"], a["search"], strict=True)]
    runs = [k for i, k in enumerate(blocks) if i == 0 or blocks[i - 1] != k]
    assert len(runs) == len(set(runs))  # each (peptide, search) block is contiguous


def test_written_coordinates_are_one_based() -> None:
    idx = _index({"c1": G})
    df = _target(_match(idx, ["MAKILW"], min_len=4))
    assert [to_one_based_closed(iv) for iv in _ivs(df)] == [(1, 18)]
