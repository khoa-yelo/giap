"""Lead decisions after code review of M3-07 (PR #48; docs/specs/M3-07-tolerant-seeds.md 1-5).

1. ``k > BAND`` (2) in ``max_diff`` is a configuration error (``ValueError``), not a silent loss.
2. The circular wrap of the six-frame index holds a 40-aa query with ``k = 2`` blocks whose window
   starts up to 12 nt before the origin (``WRAP_NT >= 3 * (40 + 2) = 126``).
4. Tests that kill the mutants that survived review: the > 16-loci serial recompute, the merge
   cover ``c + 2``, the target and decoy-type filters (I/L twins), the replicon-length check
   (a window of exactly the replicon length is kept) and ``k < 0`` in ``verify_batch``.
5. ``max_diff`` length classes beyond ``MAX_QUERY`` (512) are rejected, not skipped silently.

The repro cases are the reviewer's (k3.py, wrap.py, mut.py).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from giap.config.model import TolerantConfig
from giap.mapping import sixframe
from giap.mapping.decoys import DecoyType, build_decoys
from giap.mapping.diffops import block_tables
from giap.mapping.qgram import MAX_QUERY, build_qgram_index, candidate_diagonals
from giap.mapping.sixframe import Frame, SixFrameIndex, build_sixframe, encode_query, encode_target
from giap.mapping.tolerant import match_tolerant
from giap.mapping.tolerant_dp import BAND, exchange_lookup, verify_batch
from giap.model import Interval, Replicon
from giap.rng import for_step
from giapsim.codons import CODE_11

SENSE = sorted(c for c, a in CODE_11.items() if a != "*")
BACK = {a: next(c for c in SENSE if CODE_11[c] == a) for a in set(CODE_11.values()) - {"*"}}


def _nt(peptide: str) -> str:
    return "".join(BACK[a] for a in peptide)


def _background(n: int, seed: int) -> str:
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list("ACGT"), size=n))


def _index(genome: dict[str, str], topology: str = "linear") -> SixFrameIndex:
    reps = {c: Replicon(c, len(s), topology) for c, s in genome.items()}  # type: ignore[arg-type]
    return build_sixframe(genome, reps)


def _decoys(idx: SixFrameIndex, *kinds: DecoyType) -> dict[DecoyType, tuple[Frame, ...]]:
    return {k: build_decoys(idx, k, 50, for_step("sixframe_index", 0)) for k in kinds}


def _target(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["search"] == "target"].reset_index(drop=True)


def _decoy(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["search"] == "decoy"].reset_index(drop=True)


def _ivs(df: pd.DataFrame) -> list[Interval]:
    return [
        Interval(r.contig, int(r.start), int(r.end), r.strand) for r in df.itertuples(index=False)
    ]


# --- item 1: k > BAND is a configuration error ---------------------------------------------------

K3_WINDOW = "MKTWYEPHFDGGRCMKHPWYEGGFDKCMHEPWTGGFMKHCYDE"  # reviewer k3.py
K3_QUERY = K3_WINDOW.replace("GG", "N")  # 29 aa, three N <-> GG blocks: path spans 4 diagonals


def _k3_world() -> tuple[SixFrameIndex, dict[DecoyType, tuple[Frame, ...]]]:
    genome = _background(600, 1) + _nt(K3_WINDOW) + _background(600, 2)
    idx = _index({"g": genome})
    return idx, _decoys(idx, "reverse")


@pytest.mark.parametrize("max_diff", [{"10-40": 3}, {"10-11": 1, "12-40": 3}, {"10-40": 7}])
def test_k_above_band_is_rejected(max_diff: dict[str, int]) -> None:
    assert BAND == 2
    idx, decoys = _k3_world()
    with pytest.raises(ValueError):
        cfg = TolerantConfig(max_diff=max_diff)
        match_tolerant([K3_QUERY], idx, decoys, {K3_QUERY: "reverse"}, {K3_QUERY}, cfg)


def test_k3_repro_is_rejected_not_silently_wrong() -> None:
    """Reviewer k3.py: k = 3 reported a cost-3 hit only by chance or nothing at all."""
    idx, decoys = _k3_world()
    with pytest.raises(ValueError):
        cfg = TolerantConfig(max_diff={"10-40": 3})
        match_tolerant([K3_QUERY], idx, decoys, {K3_QUERY: "reverse"}, {K3_QUERY}, cfg)


@pytest.mark.parametrize("max_diff", [{"10-40": 0}, {"10-40": 1}, {"10-11": 1, "12-40": 2}])
def test_k_up_to_band_is_accepted(max_diff: dict[str, int]) -> None:
    idx, decoys = _k3_world()
    cfg = TolerantConfig(max_diff=max_diff)
    match_tolerant([K3_QUERY], idx, decoys, {K3_QUERY: "reverse"}, {K3_QUERY}, cfg)


# --- item 5: length classes beyond MAX_QUERY are rejected ----------------------------------------


@pytest.mark.parametrize(
    "max_diff", [{"10-513": 2}, {"10-40": 2, "41-1000": 2}, {f"{MAX_QUERY + 1}-{MAX_QUERY + 9}": 1}]
)
def test_max_len_above_max_query_is_rejected(max_diff: dict[str, int]) -> None:
    assert MAX_QUERY == 512
    idx, decoys = _k3_world()
    with pytest.raises(ValueError):
        cfg = TolerantConfig(max_diff=max_diff)
        match_tolerant([K3_QUERY], idx, decoys, {K3_QUERY: "reverse"}, {K3_QUERY}, cfg)


def test_max_len_equal_to_max_query_is_accepted() -> None:
    idx, decoys = _k3_world()
    cfg = TolerantConfig(max_diff={f"10-{MAX_QUERY}": 2})
    df = match_tolerant([K3_QUERY[:20]], idx, decoys, {K3_QUERY[:20]: "reverse"}, set(), cfg)
    assert set(df["search"]) <= {"decoy"}


# --- item 2: circular wrap holds 40-aa windows crossing the origin -------------------------------

WRAP_WINDOW = "MKTWYEPHFDGGRCMKHPWYEGGFDKCMHEPWTMKHCYDEFH"  # reviewer wrap.py: 42 aa = 126 nt
WRAP_QUERY = WRAP_WINDOW.replace("GG", "N")  # 40 aa, two blocks


def test_wrap_nt_covers_longest_query_with_max_k() -> None:
    assert len(WRAP_QUERY) == 40
    assert sixframe.WRAP_NT >= 3 * (40 + 2)


@pytest.mark.parametrize("back", list(range(13)))
def test_forty_aa_two_blocks_across_origin(back: int) -> None:
    """Window starting ``back`` nt before the origin (reviewer wrap.py, back = 3 was missed)."""
    seg = _nt(WRAP_WINDOW)
    g = seg + _background(1200 - len(seg), 5)
    n = len(g)
    g2 = g[back:] + g[:back]  # the window now starts at n - back
    idx2 = build_sixframe({"g": g2}, {"g": Replicon("g", n, "circular")})
    q = WRAP_QUERY
    df2 = match_tolerant([q], idx2, _decoys(idx2, "reverse"), {q: "reverse"}, {q}, TolerantConfig())
    t2 = _target(df2)
    got = [
        (int(r.start) % n, int(r.end) - int(r.start), int(r.n_diff), r.strand)
        for r in t2.itertuples(index=False)
    ]
    assert ((n - back) % n, 126, 2, "+") in got, (back, got)
    assert all(c == 2 for _, _, c, _ in got), got


# --- item 4: mutants that survived review --------------------------------------------------------

W = "MKTWQEPHYDRF"
W_SUB = "MKTWQECHYDRF"


def test_more_than_sixteen_minimal_loci_all_reported() -> None:
    """``_CAP`` = 16 rows per query in the parallel pass; larger results are recomputed."""
    n_copies = 23
    parts = [_background(150, 100)]
    for i in range(n_copies):
        parts += [_nt(W_SUB), _background(150, 101 + i)]
    genome = "".join(parts)
    idx = _index({"c1": genome})
    df = match_tolerant([W], idx, _decoys(idx, "reverse"), {W: "reverse"}, {W}, TolerantConfig())
    tgt = _target(df)
    want = [Interval("c1", 150 + i * 186, 150 + i * 186 + 36, "+") for i in range(n_copies)]
    assert _ivs(tgt) == want
    assert tgt["n_diff"].tolist() == [1] * n_copies


def test_merged_centres_cover_only_two_more_raw_diagonals() -> None:
    """Spec: the smallest uncovered raw ``r`` gives centre ``r + 1`` covering ``r .. r + 2``.

    Query ``A + B`` against target ``A + XXX + B``: ``A`` hits diagonal ``D`` (raw ``D-1..D+1``),
    ``B`` diagonal ``D + 3`` (raw ``D+2..D+4``). Centres ``D`` (covers ``D-1..D+1``) and ``D+3``.
    """
    a, b = "MKTWQEPHYD", "FRCWNMHPEK"
    pre = "ACDEFGHKLMNPQRSTVWY"[::-1]
    target = pre + a + "GSV" + b + pre
    frames = (Frame("c", "+", 0, encode_target(target), 3 * len(target)),)
    idx = build_qgram_index(frames, 3)
    got = candidate_diagonals(encode_query(a + b), idx, 8)
    d = len(pre)
    assert got.tolist() == [d, d + 3]


def _twin_genome() -> tuple[dict[str, str], Interval, Interval]:
    pep = "MKTWQELHYDRF"
    genome = _background(600, 41) + _nt(pep) + _background(300, 42)
    genome += _nt("MKTWQELHYDRW"[::-1]) + _background(300, 43)  # reversed: a reverse-decoy hit
    return (
        {"c1": genome},
        Interval("c1", 600, 636, "+"),
        Interval("c1", 936, 972, "+"),
    )


TWIN_L = "MKTWQELHYDRW"  # cost 1 against the genome (F -> W)
TWIN_I = "MKTWQEIHYDRW"  # the same query text after I -> L


def test_target_filter_on_i_l_twins() -> None:
    genome, tgt_iv, _ = _twin_genome()
    idx = _index(genome)
    decoys = _decoys(idx, "reverse")
    df = match_tolerant(
        [TWIN_L, TWIN_I],
        idx,
        decoys,
        {TWIN_L: "reverse", TWIN_I: "reverse"},
        {TWIN_L},
        TolerantConfig(),
    )
    tgt = _target(df)
    assert set(tgt["peptide"]) == {TWIN_L}
    assert _ivs(tgt) == [tgt_iv]
    assert set(_decoy(df)["peptide"]) == {TWIN_L, TWIN_I}  # decoys searched for every query


def test_decoy_type_filter_on_i_l_twins() -> None:
    genome, _, dec_iv = _twin_genome()
    idx = _index(genome)
    decoys = _decoys(idx, "reverse", "shuffle")
    df = match_tolerant(
        [TWIN_L, TWIN_I],
        idx,
        decoys,
        {TWIN_L: "reverse", TWIN_I: "shuffle"},
        set(),
        TolerantConfig(),
    )
    dec = _decoy(df)
    rows_l = dec[dec["peptide"] == TWIN_L]
    rows_i = dec[dec["peptide"] == TWIN_I]
    assert set(rows_l["decoy_type"]) == {"reverse"}
    assert dec_iv in _ivs(rows_l)
    assert set(rows_i["decoy_type"]) <= {"shuffle"}
    assert dec_iv not in _ivs(rows_i)


def test_window_of_exactly_the_replicon_length_is_kept() -> None:
    """``3 * (b - a) == replicon length`` is a valid locus (only longer windows are dropped)."""
    pep = "MKTWQEPHYDRF"  # 12 aa = 36 nt: the whole circular replicon
    genome = _nt(pep)
    idx = _index({"r": genome}, topology="circular")
    q = "MKTWQECHYDRF"
    df = match_tolerant([q], idx, _decoys(idx, "reverse"), {q: "reverse"}, {q}, TolerantConfig())
    tgt = _target(df)
    hits = {
        (int(r.start) % 36, int(r.end) - int(r.start), r.strand, int(r.n_diff))
        for r in tgt.itertuples()
    }
    assert (0, 36, "+", 1) in hits, hits
    assert all(length <= 36 for _, length, _, _ in hits)


@pytest.mark.parametrize("k", [-1, -2, -5])
def test_verify_batch_negative_k_finds_nothing(k: int) -> None:
    ex = exchange_lookup(block_tables())
    frame = encode_target("ACDEFGHKLMNPQRSTVWY" * 3)
    query = np.ascontiguousarray(frame[5:17])
    d0s = np.array([5, 0, 19], dtype=np.int64)
    assert verify_batch(query, frame, d0s, 0, ex)[0, 0] == 0  # exact copy at k = 0
    out = verify_batch(query, frame, d0s, k, ex)
    assert out.shape == (3, 3)
    assert (out == -1).all(), out
