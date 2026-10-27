"""Tolerant matching: q-gram seeds and the target/decoy search driver (issue M3-07; §7.2.4).

Examples 1-10 of docs/specs/M3-07-tolerant-seeds.md, the q-gram index and seed rules, the
edge cases (frame boundaries, circular extension, low complexity, NOMATCH), config handling,
the decoy-type and ``decoy_of`` validation rules inherited from M3-04 (lead decision 1) and the
M3-10 Q1 split (target search only for ``target_queries``, decoy search for every query).

Conventions resolved by the test author (spec silent): rows use ``route = "tolerant"``,
``blocks`` NA and ``decoy_type`` NA on target rows (M3-04 schema); ``diff_types`` is the
``;``-joined, *sorted* list with repeats (M3-04 schema "sorted", M3-08 example
``"frameshift;substitution"``); candidate diagonals are global (``p - i`` with ``p`` a position
of ``QGramIndex.positions``) and ``frame_starts[f]`` is the global position of residue 0 of frame
``f``.
"""

from __future__ import annotations

import re
from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gippa.config.model import TolerantConfig
from gippa.io.fasta import read_genome_fasta
from gippa.mapping.decoys import DecoyType, build_decoys
from gippa.mapping.exact import match_exact
from gippa.mapping.hits import HIT_COLUMNS
from gippa.mapping.qgram import QGramIndex, build_qgram_index, candidate_diagonals, choose_q
from gippa.mapping.sixframe import (
    Frame,
    SixFrameIndex,
    build_sixframe,
    decode,
    encode_query,
    encode_target,
)
from gippa.mapping.tolerant import match_tolerant
from gippa.model import Interval, Replicon
from gippa.model.genes import gene_model_from_bakta_json
from gippa.rng import for_step
from gippasim.codons import CODE_11
from gippasim.simulate import SimConfig, simulate

CFG = TolerantConfig()
SENSE = sorted(c for c, a in CODE_11.items() if a != "*")
BACK = {a: next(c for c in SENSE if CODE_11[c] == a) for a in set(CODE_11.values()) - {"*"}}
INT_COLS = ("frame_id", "aa_start", "aa_end", "start", "end", "frame", "n_diff")
_RC = str.maketrans("ACGTN", "TGCAN")


def _rc(s: str) -> str:
    return s.translate(_RC)[::-1]


def _nt(peptide: str) -> str:
    return "".join(BACK[a] for a in peptide)


def _background(n: int, seed: int) -> str:
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list("ACGT"), size=n))


def _index(genome: dict[str, str], topology: str = "linear", wrap_nt: int = 120) -> SixFrameIndex:
    reps = {c: Replicon(c, len(s), topology) for c, s in genome.items()}  # type: ignore[arg-type]
    return build_sixframe(genome, reps, wrap_nt=wrap_nt)


def _decoys(idx: SixFrameIndex, *kinds: DecoyType) -> dict[DecoyType, tuple[Frame, ...]]:
    return {k: build_decoys(idx, k, 50, for_step("sixframe_index", 0)) for k in kinds}


def _run(
    idx: SixFrameIndex,
    queries: list[str],
    *,
    targets: set[str] | None = None,
    kind: DecoyType = "reverse",
    cfg: TolerantConfig = CFG,
) -> pd.DataFrame:
    decoys = _decoys(idx, kind)
    tq = set(queries) if targets is None else targets
    return match_tolerant(queries, idx, decoys, dict.fromkeys(queries, kind), tq, cfg)


def _target(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["search"] == "target"].reset_index(drop=True)


def _decoy(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["search"] == "decoy"].reset_index(drop=True)


def _ivs(df: pd.DataFrame) -> list[Interval]:
    return [
        Interval(r.contig, int(r.start), int(r.end), r.strand) for r in df.itertuples(index=False)
    ]


def _frame(aa: str, contig: str = "c", strand: str = "+") -> Frame:
    codes = encode_target(aa)
    return Frame(contig, strand, 0, codes, 3 * len(codes))  # type: ignore[arg-type]


# --- choose_q (examples 1-5) -------------------------------------------------------------------


@pytest.mark.parametrize(
    ("L", "k", "expected"),
    [
        (10, 1, (4, 2)),  # example 1: q = 5 gives t = 0
        (11, 1, (5, 1)),  # example 2
        (12, 2, (3, 2)),  # example 3: q = 4 gives -1
        (13, 2, (3, 3)),  # example 4
        (14, 2, (4, 1)),  # example 5
        (20, 2, (5, 4)),
        (40, 2, (5, 24)),
    ],
)
def test_examples_1_to_5_choose_q(L: int, k: int, expected: tuple[int, int]) -> None:
    assert choose_q(L, k) == expected


@pytest.mark.parametrize(("L", "k"), [(10, 2), (9, 2), (5, 1), (4, 1), (2, 0), (1, 1)])
def test_choose_q_raises_without_lossless_q(L: int, k: int) -> None:
    with pytest.raises(ValueError):
        choose_q(L, k)


def test_choose_q_matches_formula_everywhere() -> None:
    for L in range(1, 61):
        for k in range(4):
            ok = [(q, (L - q + 1) - k * (q + 1)) for q in (5, 4, 3)]
            ok = [(q, t) for q, t in ok if t >= 1]
            if ok:
                assert choose_q(L, k) == ok[0], (L, k)  # largest q first
            else:
                with pytest.raises(ValueError):
                    choose_q(L, k)


def test_every_configured_length_has_a_lossless_q() -> None:
    for L in range(10, 41):
        k = 1 if L <= 11 else 2
        q, t = choose_q(L, k)
        assert q in (3, 4, 5) and t >= 1


# --- q-gram index ------------------------------------------------------------------------------

FRAMES3 = (_frame("ACDEFGH"), _frame("*KLMN"), _frame("PQXRSTV"))


def _decode_positions(idx: QGramIndex, frames: tuple[Frame, ...]) -> list[tuple[int, int]]:
    starts = np.asarray(idx.frame_starts)[: len(frames)]
    out = []
    for p in np.asarray(idx.positions):
        fid = int(np.searchsorted(starts, p, side="right")) - 1
        out.append((fid, int(p) - int(starts[fid])))
    return out


def test_qgram_index_skips_stops_x_and_frame_boundaries() -> None:
    idx = build_qgram_index(FRAMES3, 3)
    assert idx.q == 3
    assert np.asarray(idx.positions).dtype == np.int32
    got = sorted(_decode_positions(idx, FRAMES3))
    # ACD CDE DEF EFG FGH | KLM LMN | RST STV ; no q-gram spans two frames, none holds * or X
    expected = [(0, i) for i in range(5)] + [(1, 1), (1, 2)] + [(2, 3), (2, 4)]
    assert got == expected


def test_qgram_index_is_csr_grouped_by_qgram() -> None:
    frames = (_frame("ACDACDACDEF"), _frame("CDACDK"))
    idx = build_qgram_index(frames, 3)
    offsets = np.asarray(idx.offsets)
    assert offsets[0] == 0
    assert np.all(np.diff(offsets) >= 0)
    assert offsets[-1] == len(idx.positions)
    loci = _decode_positions(idx, frames)
    text = [decode(f.aa) for f in frames]
    grams = [text[f][s : s + 3] for f, s in loci]
    seen: dict[str, int] = {}
    for b in range(len(offsets) - 1):
        bucket = {grams[i] for i in range(int(offsets[b]), int(offsets[b + 1]))}
        assert len(bucket) <= 1, bucket  # one q-gram per CSR row
        for g in bucket:
            assert g not in seen, g  # each q-gram in one row only
            seen[g] = b
    assert grams.count("ACD") == 4 and grams.count("CDA") == 3


def test_qgram_index_folds_i_into_l() -> None:
    frames = (_frame("LLLL"), _frame("IIII"))
    idx = build_qgram_index(frames, 4)
    offsets = np.asarray(idx.offsets)
    sizes = np.diff(offsets)
    assert sorted(sizes[sizes > 0].tolist()) == [2]


@pytest.mark.parametrize("q", [3, 4, 5])
def test_qgram_index_counts_every_valid_qgram(q: int) -> None:
    rng = np.random.default_rng(q)
    frames = tuple(
        Frame("c", "+", 0, rng.integers(0, 21, size=n).astype(np.uint8), 3 * n)
        for n in (0, 2, 50, 301)
    )
    idx = build_qgram_index(frames, q)
    want = sum(
        int(np.all(f.aa[i : i + q] < 19)) for f in frames for i in range(max(0, len(f.aa) - q + 1))
    )
    assert len(idx.positions) == want
    for fid, s in _decode_positions(idx, frames):
        assert s >= 0 and s + q <= len(frames[fid].aa)
        assert np.all(frames[fid].aa[s : s + q] < 19)


# --- candidate diagonals -----------------------------------------------------------------------


def _random_frames(seed: int, sizes: tuple[int, ...] = (300, 300, 300)) -> tuple[Frame, ...]:
    rng = np.random.default_rng(seed)
    return tuple(
        Frame("c", "+", 0, rng.integers(0, 19, size=n).astype(np.uint8), 3 * n) for n in sizes
    )


def test_candidate_diagonals_exact_copy() -> None:
    frames = _random_frames(1)
    q, t = choose_q(12, 2)
    idx = build_qgram_index(frames, q)
    query = np.ascontiguousarray(frames[1].aa[100:112])
    d = int(np.asarray(idx.frame_starts)[1]) + 100
    cands = np.asarray(candidate_diagonals(query, idx, t))
    assert cands.ndim == 1 and np.issubdtype(cands.dtype, np.integer)
    assert np.any(np.abs(cands - d) <= 2), (cands, d)


def test_candidate_diagonals_none_when_threshold_unreachable() -> None:
    frames = _random_frames(2)
    idx = build_qgram_index(frames, 3)
    query = np.ascontiguousarray(frames[0].aa[10:22])
    assert len(candidate_diagonals(query, idx, 12 - 3 + 2)) == 0


def test_candidate_diagonals_count_each_query_qgram_once_per_band() -> None:
    # query q-gram AAA (offset 0) hits target positions 0, 1, 2: diagonals 0, 1, 2. The band of
    # centre 1 covers all three, yet AAA counts once: 1 shared q-gram, never 3.
    frames = (_frame("AAAAA"),)
    idx = build_qgram_index(frames, 3)
    query = encode_query("AAACDEFGHK")
    assert len(candidate_diagonals(query, idx, 2)) == 0
    assert len(candidate_diagonals(query, idx, 1)) > 0


def test_candidate_diagonals_band_absorbs_one_shift() -> None:
    # 1 -> 2 exchange (N -> GG) in the middle: q-grams left of it sit on diagonal d, right of it
    # on d + 1; only the band sum reaches the threshold.
    target = "CDEFHKMPQRSTVWY" + "GG" + "ALSDTE"
    query = "CDEFHKMPQRSTVWY" + "N" + "ALSDTE"
    frames = (_frame(target),)
    idx = build_qgram_index(frames, 3)
    q_codes = encode_query(query)
    # 13 q-grams on diagonal 0 (query offsets 0..12), 4 on diagonal 1 (offsets 16..19)
    cands = np.asarray(candidate_diagonals(q_codes, idx, 17))
    assert len(cands) > 0 and np.any(np.abs(cands - 0) <= 2) and np.any(np.abs(cands - 1) <= 2)
    assert len(candidate_diagonals(q_codes, idx, 18)) == 0


def test_candidate_diagonals_merged_at_least_three_apart() -> None:
    frames = (_frame("A" * 60), _frame("AG" * 30))
    idx = build_qgram_index(frames, 3)
    for query in ("AAAAAAAAAAAA", "AGAGAGAGAGAG", "AAAAAAGAGAGA"):
        cands = np.sort(np.asarray(candidate_diagonals(encode_query(query), idx, 2)))
        assert len(cands) > 0
        assert np.all(np.diff(cands) >= 3), cands


def test_candidate_diagonals_skip_nomatch_qgrams() -> None:
    frames = (_frame("ACDEFGHKLMNPQRSTVWY"),)
    idx = build_qgram_index(frames, 3)
    clean = encode_query("CDEFGHKLMN")  # 8 q-grams on one diagonal
    assert len(candidate_diagonals(clean, idx, 8)) > 0
    nomatch = encode_query("CDEFBHKLMN")  # B destroys 3 of them
    assert len(candidate_diagonals(nomatch, idx, 5)) > 0
    assert len(candidate_diagonals(nomatch, idx, 6)) == 0


# --- examples 6-8 on a 50-kb simulated genome --------------------------------------------------


@cache
def _sim50(tmp: str) -> tuple[SixFrameIndex, list[tuple[str, Interval]]]:
    """50-kb gippasim genome and its genes: (protein aa from offset 1, gene-origin interval)."""
    res = simulate(SimConfig(genome_length=50_000), seed=7, outdir=Path(tmp))
    (rec,) = read_genome_fasta(res.genome_fasta)
    rep = Replicon(rec.id, len(rec.sequence), "circular")
    model = gene_model_from_bakta_json(res.bakta_json, {rec.id: rep})
    idx = build_sixframe({rec.id: rec.sequence}, {rec.id: rep})
    genes = [(g.aa, g.interval) for g in model.genes if "*" not in g.aa[:-1]]
    return idx, genes


@pytest.fixture(scope="module")
def sim50(tmp_path_factory: pytest.TempPathFactory) -> tuple[SixFrameIndex, list]:  # type: ignore[type-arg]
    return _sim50(str(tmp_path_factory.mktemp("sim50")))


def _window_iv(gene_iv: Interval, off: int, n: int, length: int) -> Interval:
    """Genomic interval of protein residues ``[off, off + n)`` of a gene."""
    s = gene_iv.start + 3 * off if gene_iv.strand == "+" else gene_iv.end - 3 * (off + n)
    s %= length
    return Interval(gene_iv.contig, s, s + 3 * n, gene_iv.strand)


def _frame_texts(idx: SixFrameIndex) -> list[str]:
    return [decode(f.aa) for f in idx.frames]


def _count_exact(texts: list[str], window: str) -> int:
    w = window.replace("I", "L")
    return sum(len(re.findall(f"(?={re.escape(w)})", t)) for t in texts)


def _pick_window(
    sim: tuple[SixFrameIndex, list],
    n: int,
    need: str | None,  # type: ignore[type-arg]
) -> tuple[str, Interval]:
    """First unique ``n``-aa gene window (offset >= 1) holding ``need`` once, at 4..n-5."""
    idx, genes = sim
    texts = _frame_texts(idx)
    length = next(iter(idx.replicons.values())).length
    for aa, iv in genes:
        for off in range(1, len(aa) - n - 1, 7):
            w = aa[off : off + n]
            if "*" in w or "I" in w or "L" in w:
                continue
            if need is not None and (w.count(need) != 1 or not 4 <= w.index(need) <= n - 5):
                continue
            if need is not None and ("G" in w or "N" in w.replace(need, "")):
                continue
            if _count_exact(texts, w) == 1:
                return w, _window_iv(iv, off, n, length)
    raise AssertionError("no suitable gene window in the simulated genome")


def _one_target(df: pd.DataFrame, peptide: str) -> pd.Series:
    tgt = _target(df)
    tgt = tgt[tgt["peptide"] == peptide]
    assert len(tgt) == 1, tgt
    return tgt.iloc[0]


def test_example_6_block_n_to_gg(sim50: tuple[SixFrameIndex, list]) -> None:  # type: ignore[type-arg]
    idx, _ = sim50
    window, want = _pick_window(sim50, 12, "N")
    query = window.replace("N", "GG")
    df = _run(idx, [query])
    row = _one_target(df, query)
    assert _ivs(pd.DataFrame([row])) == [want]
    assert row["n_diff"] == 1
    assert row["diff_types"] == "block"
    assert row["route"] == "tolerant"
    assert row["query"] == query
    assert pd.isna(row["decoy_type"]) and pd.isna(row["blocks"])
    assert row["frame"] == idx.frame_label(want)


def _plant_sub_swap(window: str) -> str:
    s = list(window)
    j = next(j for j in range(2, 5) if s[j] != s[j + 1])
    s[j], s[j + 1] = s[j + 1], s[j]
    p = 8
    s[p] = "W" if s[p] != "W" else "Y"
    return "".join(s)


def test_example_7_substitution_and_swap(sim50: tuple[SixFrameIndex, list]) -> None:  # type: ignore[type-arg]
    idx, _ = sim50
    window, want = _pick_window(sim50, 12, None)
    query = _plant_sub_swap(window)
    df = _run(idx, [query])
    row = _one_target(df, query)
    assert _ivs(pd.DataFrame([row])) == [want]
    assert row["n_diff"] == 2
    assert row["diff_types"] == "substitution;swap"


def test_example_8_ten_aa_two_substitutions_not_placed(
    sim50: tuple[SixFrameIndex, list],  # type: ignore[type-arg]
) -> None:
    idx, _ = sim50
    window, _ = _pick_window(sim50, 10, None)
    s = list(window)
    for p in (2, 7):
        s[p] = "W" if s[p] != "W" else "Y"
    query = "".join(s)
    df = _run(idx, [query, window])
    assert len(_target(df)[lambda d: d["peptide"] == query]) == 0  # k = 1 at 10 aa
    # control: the unmodified 10-aa window is found at cost 0
    assert _one_target(df, window)["n_diff"] == 0


# --- example 9: below min_len ------------------------------------------------------------------


def test_example_9_nine_aa_not_searched() -> None:
    pep9 = "MKTWQEPHY"
    genome = _background(300, 1) + _nt(pep9) + _background(300, 2)
    idx = _index({"c1": genome})
    df = _run(idx, [pep9])
    assert list(df.columns) == list(HIT_COLUMNS)
    assert len(df) == 0


def test_min_len_from_config() -> None:
    pep = "MKTWQEPHYDR"  # 11 aa
    genome = _background(300, 1) + _nt(pep) + _background(300, 2)
    idx = _index({"c1": genome})
    assert len(_target(_run(idx, [pep]))) == 1
    assert len(_run(idx, [pep], cfg=TolerantConfig(min_len=12))) == 0


def test_max_diff_from_config() -> None:
    pep = "MKTWQEPHYDRF"  # 12 aa
    genome = _background(300, 1) + _nt(pep) + _background(300, 2)
    idx = _index({"c1": genome})
    two = "MKTCQEPHYDAF"  # two substitutions
    one = "MKTCQEPHYDRF"
    cfg = TolerantConfig(max_diff={"10-11": 1, "12-40": 1})
    assert len(_target(_run(idx, [two]))) == 1
    df = _run(idx, [two, one], cfg=cfg)
    assert set(_target(df)["peptide"]) == {one}


# --- example 10: exact hit -> no tolerant target search; decoys still searched ------------------

W = "MKTWQEPHYDRF"
W_SUB = "MKTWQECHYDRF"  # the second genomic copy: one substitution (P -> C)
W_DEC = "MKTYQEPHYDRF"  # planted reversed: a reverse-decoy hit at cost 1 (W -> Y)


def _example10_genome() -> tuple[dict[str, str], Interval, Interval, Interval]:
    parts = [_background(600, 11), _nt(W), _background(300, 12), _nt(W_SUB)]
    parts += [_background(300, 13), _nt(W_DEC[::-1]), _background(300, 14)]
    genome = "".join(parts)
    a = 600
    b = a + 36 + 300
    c = b + 36 + 300
    return (
        {"c1": genome},
        Interval("c1", a, a + 36, "+"),
        Interval("c1", b, b + 36, "+"),
        Interval("c1", c, c + 36, "+"),
    )


def test_example_10_exact_peptide_skips_tolerant_target_search() -> None:
    genome, exact_iv, _, decoy_iv = _example10_genome()
    idx = _index(genome)
    decoys = _decoys(idx, "reverse")
    ex = match_exact([W], idx, decoys, {W: "reverse"})
    assert _ivs(_target(ex)) == [exact_iv]
    leftovers = {W} - set(_target(ex)["peptide"])
    df = match_tolerant([W], idx, decoys, {W: "reverse"}, leftovers, CFG)
    assert len(_target(df)) == 0
    dec = _decoy(df)
    assert _ivs(dec) == [decoy_iv]  # mirrored target interval of the reversed copy
    assert dec.loc[0, "decoy_type"] == "reverse"
    assert dec.loc[0, "n_diff"] == 1
    assert dec.loc[0, "diff_types"] == "substitution"
    assert dec.loc[0, "route"] == "tolerant"


def test_example_10_target_search_keeps_only_minimal_cost() -> None:
    genome, exact_iv, sub_iv, _ = _example10_genome()
    idx = _index(genome)
    df = _run(idx, [W], targets={W})
    tgt = _target(df)
    assert _ivs(tgt) == [exact_iv]  # the cost-1 copy at sub_iv is not minimal
    assert sub_iv not in _ivs(tgt)
    assert tgt.loc[0, "n_diff"] == 0 and tgt.loc[0, "diff_types"] == ""


def test_all_minimal_cost_loci_are_reported() -> None:
    genome = _background(300, 21) + _nt(W_SUB) + _background(300, 22)
    genome += _rc(_nt("MKTWQECHYDRW")) + _background(300, 23) + _nt(W_SUB) + _background(99, 24)
    idx = _index({"c1": genome})
    tgt = _target(_run(idx, [W]))
    assert _ivs(tgt) == [
        Interval("c1", 300, 336, "+"),
        Interval("c1", 972, 1008, "+"),
    ]
    assert tgt["n_diff"].tolist() == [1, 1]


# --- decoy types and decoy_of validation (M3-04 lead decision 1) -------------------------------


def test_each_peptide_searched_against_its_own_decoy_type() -> None:
    genome, _, _, decoy_iv = _example10_genome()
    idx = _index(genome)
    decoys = _decoys(idx, "reverse", "shuffle")
    other = "PEPTLDEKAGRW"
    queries = [W, other]
    df = match_tolerant(queries, idx, decoys, {W: "reverse", other: "shuffle"}, set(), CFG)
    dec = _decoy(df)
    assert set(dec[dec["peptide"] == W]["decoy_type"]) == {"reverse"}
    assert set(dec[dec["peptide"] == other]["decoy_type"]) <= {"shuffle"}
    assert decoy_iv in _ivs(dec[dec["peptide"] == W])
    swapped = match_tolerant(queries, idx, decoys, {W: "shuffle", other: "reverse"}, set(), CFG)
    sdec = _decoy(swapped)
    assert set(sdec[sdec["peptide"] == W]["decoy_type"]) <= {"shuffle"}
    assert decoy_iv not in _ivs(sdec[sdec["peptide"] == W])


def test_decoy_of_value_without_frames_raises() -> None:
    idx = _index({"c1": _background(600, 3)})
    decoys = _decoys(idx, "reverse")
    with pytest.raises(ValueError, match="shuffle"):
        match_tolerant([W], idx, decoys, {W: "shuffle"}, {W}, CFG)


def test_searched_peptide_missing_from_decoy_of_raises() -> None:
    idx = _index({"c1": _background(600, 3)})
    decoys = _decoys(idx, "reverse")
    with pytest.raises(ValueError, match=W):
        match_tolerant([W, W_SUB], idx, decoys, {W_SUB: "reverse"}, {W}, CFG)


def test_unsearched_short_peptide_needs_no_decoy_type() -> None:
    idx = _index({"c1": _background(600, 3)})
    decoys = _decoys(idx, "reverse")
    df = match_tolerant(["MKTWQEPHY"], idx, decoys, {}, {"MKTWQEPHY"}, CFG)
    assert len(df) == 0


# --- schema, determinism -----------------------------------------------------------------------


def test_empty_result_is_typed() -> None:
    idx = _index({"c1": _background(600, 3)})
    df = _run(idx, [])
    assert list(df.columns) == list(HIT_COLUMNS)
    full = _run(_index(_example10_genome()[0]), [W])
    assert len(full) > 0
    assert list(full.columns) == list(HIT_COLUMNS)
    for col in INT_COLS:
        assert pd.api.types.is_integer_dtype(df[col]), col
        assert pd.api.types.is_integer_dtype(full[col]), col


def test_peptide_keeps_input_letters_and_query_is_i_to_l() -> None:
    pep = "MKTWQEPHYDIF"
    genome = _background(300, 1) + _nt("MKTWQEPHYDLW") + _background(300, 2)
    idx = _index({"c1": genome})
    row = _one_target(_run(idx, [pep]), pep)
    assert row["query"] == "MKTWQEPHYDLF"
    assert row["n_diff"] == 1  # I = L is free; F -> W is the only difference


def test_deterministic_and_independent_of_query_order() -> None:
    genome, *_ = _example10_genome()
    idx = _index(genome)
    queries = [W, W_SUB, W_DEC, "MKTWQEPHYDRFAG", "PEPTLDEKAGRW", "MKTWQEPHY"]
    a = _run(idx, queries)
    b = _run(idx, list(reversed(queries)))
    c = _run(idx, queries + queries)
    pd.testing.assert_frame_equal(a, b)
    pd.testing.assert_frame_equal(a, c)
    pd.testing.assert_frame_equal(a, _run(idx, queries))


# --- edge cases --------------------------------------------------------------------------------


def test_circular_origin_spanning_window_reported_once() -> None:
    pep_nt = _nt(W_SUB)
    bg = _background(900, 31)
    genome = pep_nt[15:] + bg + pep_nt[:15]  # the window crosses the origin
    L = len(genome)
    idx = _index({"r": genome}, topology="circular")
    tgt = _target(_run(idx, [W]))
    assert _ivs(tgt) == [Interval("r", L - 15, L + 21, "+")]
    assert tgt.loc[0, "n_diff"] == 1


def test_circular_extension_copy_deduplicated() -> None:
    bg = _background(900, 32)
    genome = bg[:30] + _nt(W_SUB) + bg[30:]  # inside the first 120 nt: also in the extension
    idx = _index({"r": genome}, topology="circular")
    df = _run(idx, [W])
    tgt = _target(df)
    assert _ivs(tgt) == [Interval("r", 30, 66, "+")]
    dec = _decoy(df)
    assert len(dec) == len(set(_ivs(dec)))


def test_minus_strand_window() -> None:
    genome = _background(300, 41) + _rc(_nt(W_SUB)) + _background(300, 42)
    idx = _index({"c1": genome})
    tgt = _target(_run(idx, [W]))
    assert _ivs(tgt) == [Interval("c1", 300, 336, "-")]
    assert tgt.loc[0, "n_diff"] == 1


def test_windows_at_contig_ends() -> None:
    first = "MKTWQEPHYDRF"
    last = "CHWPKDYQTMEF"
    genome = _nt(first) + _background(600, 51) + _nt(last)
    idx = _index({"c1": genome})
    q1 = "WKTWQEPHYDRF"  # substitution at residue 0
    q2 = "CHWPKDYQTMEW"  # substitution at the last residue
    df = _target(_run(idx, [q1, q2]))
    got = {r.peptide: (int(r.start), int(r.end), int(r.n_diff)) for r in df.itertuples()}
    n = len(genome)
    assert got == {q1: (0, 36, 1), q2: (n - 36, n, 1)}


def test_low_complexity_queries_stay_correct() -> None:
    run = "A" * 30 + "GK"
    genome = _background(300, 61) + _nt(run) + _background(300, 62)
    idx = _index({"c1": genome})
    tail = Interval("c1", 300 + 3 * 20, 300 + 3 * 32, "+")
    exact_tail = "AAAAAAAAAAGK"  # cost 0 only at the end of the run
    sub_tail = "AAAAAAAAAWGK"  # cost 1 only at the end of the run
    poly = "AAAAAAAAAAAA"  # cost 0 at many overlapping windows
    df = _target(_run(idx, [exact_tail, sub_tail, poly]))
    by = {p: g for p, g in df.groupby("peptide")}
    assert _ivs(by[exact_tail]) == [tail] and by[exact_tail]["n_diff"].tolist() == [0]
    assert _ivs(by[sub_tail]) == [tail] and by[sub_tail]["n_diff"].tolist() == [1]
    assert len(by[poly]) >= 1 and set(by[poly]["n_diff"]) == {0}
    for iv in _ivs(by[poly]):  # GCA repeats read as poly-A on + and (GCT) on -
        assert iv.start < 390 and iv.end > 300, iv


def test_nomatch_residue_costs_a_substitution() -> None:
    genome = _background(300, 71) + _nt(W) + _background(300, 72)
    idx = _index({"c1": genome})
    q = "MKTWQBPHYDRF"
    row = _one_target(_run(idx, [q]), q)
    assert row["n_diff"] == 1 and row["diff_types"] == "substitution"
    assert _ivs(pd.DataFrame([row])) == [Interval("c1", 300, 336, "+")]
    many = "BKTBQEBHYBRF"  # too many NOMATCH residues: may find nothing, never an error
    df = _run(idx, [many])
    assert set(df["n_diff"]) <= {0, 1, 2}


def test_no_duplicate_target_loci_with_exact_matcher(
    sim50: tuple[SixFrameIndex, list],  # type: ignore[type-arg]
) -> None:
    idx, genes = sim50
    rng = np.random.default_rng(5)
    queries = []
    for aa, _ in genes[:40]:
        off = int(rng.integers(1, len(aa) - 20))
        w = aa[off : off + int(rng.integers(10, 19))]
        if "*" in w:
            continue
        queries.append(w)
        s = list(w)
        s[3] = "W" if s[3] != "W" else "Y"
        queries.append("".join(s))
    decoys = _decoys(idx, "reverse")
    decoy_of = dict.fromkeys(queries, "reverse")
    ex = match_exact(queries, idx, decoys, decoy_of)
    leftovers = set(queries) - set(_target(ex)["peptide"])
    assert 0 < len(leftovers) < len(set(queries))
    tol = match_tolerant(queries, idx, decoys, decoy_of, leftovers, CFG)  # type: ignore[arg-type]
    assert set(_target(tol)["peptide"]) <= leftovers
    both = pd.concat([_target(ex), _target(tol)])
    keys = list(
        zip(
            both["peptide"], both["contig"], both["start"], both["end"], both["strand"], strict=True
        )
    )
    assert len(keys) == len(set(keys))
    assert set(_decoy(tol)["peptide"]) <= set(queries)


# --- short queries: verify raises for L <= 4 (M3-06 lead decision 1) ---------------------------


def test_queries_of_at_most_four_residues_are_never_verified() -> None:
    # A config that lowers min_len below the verifier's limit must not crash the driver: queries
    # of <= 2 * BAND (4) residues are not tolerant-searched; 5 residues are (k = 0, q = 5, t = 1).
    pep5, pep4 = "MKTWQ", "HYDR"
    genome = _background(300, 81) + _nt(pep5 + "EP" + pep4) + _background(300, 82)
    idx = _index({"c1": genome})
    cfg = TolerantConfig(min_len=3, max_diff={"3-40": 0})
    df = _run(idx, [pep5, pep4, "KTW"], cfg=cfg)
    assert set(df["peptide"]) == {pep5}
    assert _ivs(_target(df)[lambda d: d["peptide"] == pep5]) == [Interval("c1", 300, 315, "+")]


# --- input normalisation (M3-06 lead decision 4: one uint8 C-contiguous specialisation) -------


def test_candidate_diagonals_accept_any_integer_dtype_and_stride() -> None:
    frames = _random_frames(9)
    q, t = choose_q(14, 2)
    idx = build_qgram_index(frames, q)
    base = np.ascontiguousarray(frames[2].aa[40:54])
    want = np.asarray(candidate_diagonals(base, idx, t))
    assert len(want) > 0
    strided = np.repeat(base, 2)[::2]
    assert not strided.flags["C_CONTIGUOUS"] or strided.strides != base.strides
    for variant in (base.astype(np.int64), base.astype(np.int32), strided):
        np.testing.assert_array_equal(np.asarray(candidate_diagonals(variant, idx, t)), want)


# --- public njit batch entry point over the verifier (M3-06 lead decision 6) -------------------
#
# Interface pinned by these tests (the M3-07 spec defers it here):
#   tolerant_dp.exchange_lookup(tables) -> int8 array   (the kernel's exchange table)
#   tolerant_dp.verify_batch(query, frame, d0s, k, ex)  numba-jitted; query/frame uint8 C-contig,
#       d0s int64; returns int64 (len(d0s), 3) rows (key, t_start, t_end), key -1 when cost > k
#   tolerant_dp.decode_key(key) -> (cost, diff_types)   (the same decoding as verify)


def _batch_api():  # type: ignore[no-untyped-def]
    from gippa.mapping import tolerant_dp

    for name in ("verify_batch", "exchange_lookup", "decode_key"):
        assert name in tolerant_dp.__all__, name
    return tolerant_dp


def test_verify_batch_is_a_public_numba_dispatcher() -> None:
    from numba.core.registry import CPUDispatcher

    td = _batch_api()
    assert isinstance(td.verify_batch, CPUDispatcher)


def test_verify_batch_agrees_with_verify_row_by_row() -> None:
    from gippa.mapping.diffops import block_tables
    from gippa.mapping.tolerant_dp import verify

    td = _batch_api()
    tables = block_tables()
    ex = td.exchange_lookup(tables)
    rng = np.random.default_rng(13)
    n_hits = 0
    for trial in range(40):
        frame = rng.integers(0, 21, size=int(rng.integers(0, 200))).astype(np.uint8)
        frame[rng.random(len(frame)) < 0.85] %= 19
        L = int(rng.integers(5, 25))
        if len(frame) > L and trial % 2 == 0:
            s = int(rng.integers(0, len(frame) - L))
            query = frame[s : s + L].copy() % 19
            query[int(rng.integers(L))] = int(rng.integers(19))
        else:
            query = rng.integers(0, 19, size=L).astype(np.uint8)
        query = np.ascontiguousarray(query, dtype=np.uint8)
        d0s = np.arange(-L - 3, len(frame) + 4, dtype=np.int64)
        k = int(rng.integers(0, 3))
        out = np.asarray(td.verify_batch(query, frame, d0s, k, ex))
        assert out.shape == (len(d0s), 3) and out.dtype == np.int64
        for d0, (key, s, e) in zip(d0s, out, strict=True):
            want = verify(query, frame, int(d0), k, tables)
            if not want:
                assert key == -1, (trial, d0)
                continue
            n_hits += 1
            cost, types = td.decode_key(int(key))
            assert (cost, int(s), int(e), tuple(types)) == (
                want[0].cost,
                want[0].t_start,
                want[0].t_end,
                want[0].diff_types,
            ), (trial, d0)
    assert n_hits > 50


def test_verify_batch_is_callable_from_jitted_code() -> None:
    from numba import njit, prange

    from gippa.mapping.diffops import block_tables

    td = _batch_api()
    vb = td.verify_batch
    ex = td.exchange_lookup(block_tables())

    @njit(parallel=True)
    def many(queries, frame, d0s, k, ex):  # type: ignore[no-untyped-def]
        out = np.empty((queries.shape[0], d0s.shape[0], 3), dtype=np.int64)
        for i in prange(queries.shape[0]):
            out[i] = vb(queries[i], frame, d0s, k, ex)
        return out

    frame = encode_target("ACDEFGHKLMNPQRSTVWY" * 3)
    queries = np.stack([frame[i : i + 12] for i in (0, 5, 19)]).astype(np.uint8)
    d0s = np.array([0, 5, 19, 30], dtype=np.int64)
    got = many(queries, frame, d0s, 0, ex)
    for i, d in enumerate((0, 5, 19)):
        direct = np.asarray(vb(np.ascontiguousarray(queries[i]), frame, d0s, 0, ex))
        np.testing.assert_array_equal(got[i], direct)
        row = int(np.nonzero(d0s == d)[0][0])
        assert tuple(got[i, row, 1:]) == (d, d + 12)
        assert td.decode_key(int(got[i, row, 0])) == (0, ())


def test_verify_batch_rejects_queries_of_at_most_four_residues() -> None:
    from gippa.mapping.diffops import block_tables

    td = _batch_api()
    ex = td.exchange_lookup(block_tables())
    frame = encode_target("AKGAKAKGAK")
    with pytest.raises(ValueError):
        td.verify_batch(encode_query("AKGA"), frame, np.array([0], dtype=np.int64), 1, ex)


def _assert_codes(arr: object, what: str) -> None:
    assert isinstance(arr, np.ndarray), (what, type(arr))
    assert arr.dtype == np.uint8 and arr.ndim == 1 and arr.flags.c_contiguous, (what, arr.dtype)


def test_driver_passes_uint8_contiguous_codes_to_the_kernels(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Lead decision 3 after review: argument types seen at the Python boundary of
    ``_search_many`` / ``_search_one`` (independent of a cold or warm numba cache)."""
    from gippa.mapping import tolerant

    genome, *_ = _example10_genome()
    idx = _index(genome)
    queries = [W, W_SUB, "MKTWQEPHYDRF", "AAAAAAAAAAAA"]
    _run(idx, queries)  # compile (or load) the kernels before they are wrapped
    calls: dict[str, list[tuple[object, ...]]] = {"many": [], "one": []}
    many, one = tolerant._search_many, tolerant._search_one

    def spy_many(*args: object) -> object:
        calls["many"].append(args)
        return many(*args)

    def spy_one(*args: object) -> object:
        calls["one"].append(args)
        return one(*args)

    monkeypatch.setattr(tolerant, "_search_many", spy_many)
    monkeypatch.setattr(tolerant, "_search_one", spy_one)
    many_copies = _background(150, 7) + "".join(
        _nt(W_SUB) + _background(150, 8 + i) for i in range(20)
    )
    _run(idx, queries)
    _run(_index({"m": many_copies}), [W])  # > 16 loci: the serial _search_one path
    assert calls["many"], "match_tolerant does not call _search_many"
    assert calls["one"], "more than 16 loci did not take the serial _search_one path"
    for args in calls["many"]:
        _assert_codes(args[0], "qcat")
        _assert_codes(args[7], "text")
    for args in calls["one"]:
        _assert_codes(args[0], "query codes")
        _assert_codes(args[6], "text")
