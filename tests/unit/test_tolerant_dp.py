"""Unit tests for the banded tolerant verifier (issue M3-06; GIAP_SPEC.md §7.2.4, G4).

Spec examples 1-10, the edge cases of docs/specs/M3-06-tolerant-verifier.md, both orientations of
every ``BlockTables`` row (rows are directional; the DP tries both, lead decision 5 of PR #35),
stop / ``X`` / query NOMATCH (decodes to ``#``) never matched, frame ends, an origin-crossing hit
on a circular replicon, the band limits, ``k = 0`` as exact matching, and determinism.

``diff_types`` is a tuple sorted by the precedence of ``DIFF_TYPES``
(substitution < swap < block < near_block), with repeats.
"""

from __future__ import annotations

import dataclasses
from itertools import permutations

import numpy as np
import pytest

from giap.mapping.diffops import BLOCKS, DIFF_TYPES, NEAR_BLOCKS, block_tables
from giap.mapping.sixframe import build_sixframe, decode, encode_query, encode_target
from giap.mapping.tolerant_dp import Alignment, verify
from giap.model import Replicon

TABLES = block_tables()
ORDER = {name: i for i, name in enumerate(DIFF_TYPES)}


def run(query: str, frame: str, d0: int, k: int) -> list[Alignment]:
    return verify(encode_query(query), encode_target(frame), d0, k, TABLES)


def one(query: str, frame: str, d0: int, k: int) -> Alignment:
    out = run(query, frame, d0, k)
    assert len(out) == 1, out  # every window in a +-2 band overlaps: one locus
    return out[0]


# (query, k, target window, expected cost, expected diff types); None = no alignment
EXAMPLES = [
    ("AKLMNPQRSTVW", 2, "AKLMPNQRSTVW", 1, ("swap",)),
    ("PEPTNDEKAR", 1, "PEPTGGDEKAR", 1, ("block",)),
    ("VTEGTAPKLR", 1, "VTEASAPKLR", 1, ("block",)),
    ("SWAPKELLTR", 1, "SADAPKELLTR", 1, ("near_block",)),
    ("MKRLLAEQTR", 1, "MKGVLLAEQTR", 1, ("near_block",)),
    ("PEPTIDEKAR", 1, "PEPTLDEKAR", 0, ()),
    ("PEPTIDEKAR", 1, "PAPTIDEKGR", None, None),
    ("AKLMNPQRSTVW", 2, "AKLMPNQRSTVF", 2, ("substitution", "swap")),
    ("PEPTIDEKAR", 1, "PEPT*DEKAR", None, None),
    ("AGLLKDEPTR", 1, "QLLKDEPTR", 1, ("block",)),
]


@pytest.mark.parametrize("flank", ["", "CHY", "**"], ids=["bare", "flanked", "stops"])
@pytest.mark.parametrize(
    ("query", "k", "window", "cost", "types"), EXAMPLES, ids=[f"ex{i}" for i in range(1, 11)]
)
def test_spec_examples(
    query: str,
    k: int,
    window: str,
    cost: int | None,
    types: tuple[str, ...] | None,
    flank: str,
) -> None:
    frame = flank + window + flank[::-1]
    out = run(query, frame, len(flank), k)
    if cost is None:
        assert out == []
        return
    assert len(out) == 1, out
    aln = out[0]
    assert aln.cost == cost
    assert aln.diff_types == types
    assert (aln.t_start, aln.t_end) == (len(flank), len(flank) + len(window))


def test_example2_window_length_11_and_example10_length_9() -> None:
    a2 = one("PEPTNDEKAR", "PEPTGGDEKAR", 0, 1)
    a10 = one("AGLLKDEPTR", "QLLKDEPTR", 0, 1)
    assert a2.t_end - a2.t_start == 11
    assert a10.t_end - a10.t_start == 9


def test_example3_is_one_block_not_two_substitutions() -> None:
    assert one("VTEGTAPKLR", "VTEASAPKLR", 0, 2).diff_types == ("block",)


def test_example7_needs_k2() -> None:
    aln = one("PEPTIDEKAR", "PAPTIDEKGR", 0, 2)
    assert (aln.cost, aln.diff_types) == (2, ("substitution", "substitution"))


def test_alignment_is_frozen_dataclass() -> None:
    aln = one("PEPTIDEKAR", "PEPTLDEKAR", 0, 1)
    assert isinstance(aln, Alignment)
    assert [f.name for f in dataclasses.fields(Alignment)] == [
        "cost",
        "t_start",
        "t_end",
        "diff_types",
    ]
    with pytest.raises(dataclasses.FrozenInstanceError):
        aln.cost = 5  # type: ignore[misc]
    assert isinstance(aln.diff_types, tuple)


# ---- every block / near-block pair, every residue order, both orientations ---------------------

PREFIX, SUFFIX = "MKPHC", "FYHPKM"  # residues outside every exchange pair


def _exchange_cases() -> list[tuple[str, str, str]]:
    cases: list[tuple[str, str, str]] = []
    for pairs, kind in ((BLOCKS, "block"), (NEAR_BLOCKS, "near_block")):
        for a, b in pairs:
            for x in sorted({"".join(p) for p in permutations(a)}):
                for y in sorted({"".join(p) for p in permutations(b)}):
                    cases.append((x, y, kind))  # query side a, target side b
                    cases.append((y, x, kind))  # the reverse orientation
    return cases


EXCHANGES = _exchange_cases()


@pytest.mark.parametrize(
    ("qseg", "tseg", "kind"), EXCHANGES, ids=[f"{q}>{t}" for q, t, _ in EXCHANGES]
)
def test_every_exchange_in_both_orientations_costs_one(qseg: str, tseg: str, kind: str) -> None:
    query, window = PREFIX + qseg + SUFFIX, PREFIX + tseg + SUFFIX
    aln = one(query, "CHY" + window + "YHC", 3, 1)
    assert (aln.cost, aln.diff_types) == (1, (kind,))
    assert (aln.t_start, aln.t_end) == (3, 3 + len(window))


@pytest.mark.parametrize(
    ("qseg", "tseg", "kind"), EXCHANGES, ids=[f"{q}>{t}" for q, t, _ in EXCHANGES]
)
def test_exchange_at_query_ends(qseg: str, tseg: str, kind: str) -> None:
    """Free target ends: when the target side is longer, the shorter window explains the end
    residue by one substitution at the same cost, and the precedence keeps the substitution."""
    body = "MKPHCFYHPKM"
    start = one(qseg + body, tseg + body, 0, 1)
    end = one(body + qseg, body + tseg, 0, 1)
    if len(tseg) > len(qseg):
        assert (start.cost, start.diff_types) == (1, ("substitution",))
        assert (end.cost, end.diff_types) == (1, ("substitution",))
        assert (start.t_start, start.t_end) == (1, len(tseg) + len(body))
        assert (end.t_start, end.t_end) == (0, len(body) + 1)
    else:
        assert (start.cost, start.diff_types) == (1, (kind,))
        assert (end.cost, end.diff_types) == (1, (kind,))
        assert (start.t_start, start.t_end) == (0, len(tseg) + len(body))
        assert (end.t_start, end.t_end) == (0, len(tseg) + len(body))


def test_block_sides_with_isoleucine_fold_to_leucine() -> None:
    # {A,V} <-> {G,L}: a query "GI" is "GL" after the I -> L fold
    aln = one("PEPTLDEKGIR", "PEPTLDEKAVR", 0, 1)
    assert (aln.cost, aln.diff_types) == (1, ("block",))
    aln = one("PEPTLDEKAVR", "PEPTLDEKIGR", 0, 1)
    assert (aln.cost, aln.diff_types) == (1, ("block",))


def test_il_free_everywhere() -> None:
    aln = one("IIKLLIIKLLMR", "LLKIILLKIIMR", 0, 0)
    assert (aln.cost, aln.diff_types) == (0, ())


def test_swap_needs_distinct_residues_and_aa_is_a_match() -> None:
    aln = one("AAKLMNPQRSTV", "AAKLMNPQRSTV", 0, 2)
    assert (aln.cost, aln.diff_types) == (0, ())
    aln = one("PEPTAADEKAR", "PEPTAADEKAR", 0, 1)
    assert (aln.cost, aln.diff_types) == (0, ())


def test_two_length_changing_blocks_move_the_diagonal_by_two() -> None:
    longer = one("AKNLMPQNSTVW", "CH" + "AKGGLMPQGGSTVW" + "HC", 2, 2)
    assert (longer.cost, longer.diff_types) == (2, ("block", "block"))
    assert (longer.t_start, longer.t_end) == (2, 16)
    shorter = one("AKGGLMPQGGSTVW", "CH" + "AKNLMPQNSTVW" + "HC", 2, 2)
    assert (shorter.cost, shorter.diff_types) == (2, ("block", "block"))
    assert (shorter.t_start, shorter.t_end) == (2, 14)


def test_two_near_blocks_in_one_direction() -> None:
    aln = one("MKWLHPCRYHPK", "MKADLHPCGVYHPK", 0, 2)
    assert (aln.cost, aln.diff_types) == (2, ("near_block", "near_block"))
    assert aln.t_end == 14


def test_mixed_types_sorted_by_precedence_not_alphabet() -> None:
    aln = one("PEPTNDEKARW", "PEPTGGDEKARF", 0, 2)
    assert (aln.cost, aln.diff_types) == (2, ("substitution", "block"))
    assert list(aln.diff_types) == sorted(aln.diff_types, key=ORDER.__getitem__)


# ---- k limit ----------------------------------------------------------------------------------


def test_cost_above_k_returns_empty() -> None:
    assert run("AKLMNPQRSTVW", "AKLMPNQRSTVF", 0, 1) == []
    assert run("AKLMNPQRSTVW", "AKLMPNQRSTVW", 0, 0) == []
    assert run("AKNLMPQNSTVW", "AKGGLMPQGGSTVW", 0, 1) == []


def test_reported_cost_is_minimal_and_at_most_k() -> None:
    for k in range(4):
        out = run("AKLMNPQRSTVW", "AKLMPNQRSTVF", 0, k)
        if k < 2:
            assert out == []
        else:
            assert [a.cost for a in out] == [2]


# ---- stop, X and query NOMATCH ----------------------------------------------------------------


@pytest.mark.parametrize("k", [1, 2])
def test_target_x_is_never_consumed(k: int) -> None:
    assert run("PEPTIDEKAR", "PEPTXDEKAR", 0, k) == []


@pytest.mark.parametrize("k", [1, 2])
def test_stop_inside_a_block_window_is_never_consumed(k: int) -> None:
    assert run("PEPTNDEKAR", "PEPTG*DEKAR", 0, k) == []  # N <-> G* is no block
    assert run("PEPTWDEKAR", "PEPTA*DEKAR", 0, k) == []  # W <-> A* is no near-block
    assert run("PEPTNDEKAR", "PEPT*NDEKAR", 0, k) == []


def test_swap_never_consumes_a_stop_or_x() -> None:
    assert run("PEPTLDEKARMK", "PEPTLDE*ARMK", 0, 2) == []
    assert run("PEPTLDEKARMK", "PEPTLDXKARMK", 0, 2) == []


def test_query_nomatch_only_by_substitution() -> None:
    assert decode(encode_query("PEPTBDEKAR"))[4] == "#"
    aln = one("PEPTBDEKAR", "PEPTLDEKAR", 0, 1)
    assert (aln.cost, aln.diff_types) == (1, ("substitution",))
    assert run("PEPTBDEKAR", "PEPTLDEKAR", 0, 0) == []


@pytest.mark.parametrize("letter", ["X", "B", "*", "U"])
@pytest.mark.parametrize("target", ["X", "*"])
def test_query_nomatch_never_matches_target_x_or_stop(letter: str, target: str) -> None:
    for k in (0, 1, 2):
        assert run(f"PEPT{letter}DEKAR", f"PEPT{target}DEKAR", 0, k) == []


def test_nomatch_query_never_in_a_block() -> None:
    # "#" + G is never {A,G}; the query J (NOMATCH) vs Q must be a substitution, not a block
    aln = one("PEPTJGDEKAR", "PEPTQGDEKAR", 0, 1)
    assert (aln.cost, aln.diff_types) == (1, ("substitution",))


# ---- frame ends, band, ties --------------------------------------------------------------------


def test_window_never_runs_past_the_frame_end() -> None:
    aln = one("PEPTLDEKGG", "CHYPEPTLDEKN", 3, 1)  # GG <-> N with the window at the frame end
    assert (aln.cost, aln.diff_types, aln.t_end) == (1, ("block",), 12)
    aln = one("PEPTLDEKARN", "PEPTLDEKARG", 0, 1)  # N <-> GG would need a residue past the end
    assert (aln.cost, aln.diff_types) == (1, ("substitution",))
    assert (aln.t_start, aln.t_end) == (0, 11)


def test_window_never_starts_before_the_frame() -> None:
    aln = one("NPEPTLDEKR", "GPEPTLDEKR", 0, 1)
    assert (aln.cost, aln.diff_types, aln.t_start, aln.t_end) == (1, ("substitution",), 0, 10)


@pytest.mark.parametrize("d0", [-2, -1, 0, 1, 2])
def test_peptide_at_frame_start_with_band_past_the_start(d0: int) -> None:
    aln = one("PEPTLDEKAR", "PEPTLDEKARMMMMM", d0, 0)
    assert (aln.cost, aln.t_start, aln.t_end) == (0, 0, 10)


@pytest.mark.parametrize("shift", [-2, -1, 0, 1, 2])
def test_peptide_at_frame_end_with_band_past_the_end(shift: int) -> None:
    frame = "MMMMM" + "PEPTLDEKAR"
    aln = one("PEPTLDEKAR", frame, 5 + shift, 0)
    assert (aln.cost, aln.t_start, aln.t_end) == (0, 5, 15)


@pytest.mark.parametrize("shift", [-3, 3, 5])
def test_band_is_plus_minus_two(shift: int) -> None:
    frame = "MMMMMMMM" + "PEPTLDEKAR" + "MMMMMMMM"
    assert run("PEPTLDEKAR", frame, 8 + shift, 0) == []


def test_equal_cost_overlapping_windows_keep_the_lowest_precedence() -> None:
    # window ...G (substitution N->G) and window ...GG (block N->GG) both cost 1
    aln = one("PEPTLDEKARN", "PEPTLDEKARGG", 0, 1)
    assert (aln.cost, aln.diff_types, aln.t_end) == (1, ("substitution",), 11)
    aln = one("PEPTLDEKARQ", "PEPTLDEKARAG", 0, 1)
    assert (aln.cost, aln.diff_types, aln.t_end) == (1, ("substitution",), 11)


def test_equal_windows_keep_the_smallest_t_start() -> None:
    aln = one("GGGGGGGGGG", "GGGGGGGGGGGG", 1, 0)
    assert (aln.cost, aln.diff_types, aln.t_start, aln.t_end) == (0, (), 0, 10)


# ---- k = 0 is exact matching ------------------------------------------------------------------


def test_k0_equals_exact_matching_on_random_frames() -> None:
    rng = np.random.default_rng(6)
    for _ in range(300):
        frame = rng.integers(0, 21, size=40).astype(np.uint8)  # residues, stops and X
        L = int(rng.integers(7, 13))
        s = int(rng.integers(0, 40 - L + 1))
        query = frame[s : s + L].copy()
        if rng.random() < 0.5:
            query[int(rng.integers(L))] = rng.integers(0, 19)
        query[query >= 19] = 21  # a query cannot hold a stop or X: NOMATCH
        out = verify(query, frame, s, 0, TABLES)
        starts = [
            t
            for t in range(max(0, s - 2), min(40 - L, s + 2) + 1)
            if np.array_equal(frame[t : t + L], query)
        ]
        if starts:
            assert [(a.cost, a.t_start, a.t_end, a.diff_types) for a in out] == [
                (0, starts[0], starts[0] + L, ())
            ]
        else:
            assert out == []


# ---- origin-crossing hits on circular replicons -----------------------------------------------

_CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
    "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCG", "Q": "CAG", "R": "CGT",
    "S": "TCT", "T": "ACC", "V": "GTT", "W": "TGG", "Y": "TAT",
}  # fmt: skip


def _revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _protein(n: int, seed: int) -> str:
    rng = np.random.default_rng(seed)
    letters = "ACDEFGHKLMNPQRSTVWY"
    return "".join(letters[i] for i in rng.integers(0, 19, size=n))


@pytest.mark.parametrize("strand", ["+", "-"])
def test_origin_crossing_hit_on_circular_replicon(strand: str) -> None:
    protein = _protein(100, 3)
    cds = "".join(_CODON[a] for a in protein)  # 300 nt, read from the origin
    genome = cds if strand == "+" else _revcomp(cds)
    rep = Replicon("c", 300, "circular")
    idx = build_sixframe({"c": genome}, {"c": rep})
    true = protein[95:] + protein[:7]  # 12 aa across the origin
    query = true[:5] + ("W" if true[5] != "W" else "F") + true[6:]  # one substitution
    hits = []
    for fid, fr in enumerate(idx.frames):
        if fr.strand != strand:
            continue
        pos = decode(fr.aa).find(true)
        if pos < 0:
            continue
        aln = one(query, decode(fr.aa), pos, 1)
        assert (aln.cost, aln.diff_types) == (1, ("substitution",))
        assert (aln.t_start, aln.t_end) == (pos, pos + 12)
        iv = idx.to_interval(fid, aln.t_start, aln.t_end)
        hits.append((iv, idx.frame_label(iv)))
    assert hits, "the origin-crossing window must exist in a wrapped frame"
    for iv, label in hits:
        assert iv.length_bp == 36
        assert iv.end > 300  # unwrapped origin-crossing interval (spec §6.5)
        assert label == (1 if strand == "+" else -1)


def test_hit_absent_from_unwrapped_linear_replicon() -> None:
    protein = _protein(100, 3)
    cds = "".join(_CODON[a] for a in protein)
    idx = build_sixframe({"c": cds}, {"c": Replicon("c", 300, "linear")})
    frame = decode(idx.frames[0].aa)
    query = protein[95:] + protein[:7]
    assert run(query, frame, 95, 2) == []


# ---- determinism and purity -------------------------------------------------------------------


def test_deterministic_and_does_not_mutate_inputs() -> None:
    q = encode_query("AKLMNPQRSTVW")
    f = encode_target("CHYAKLMPNQRSTVFYHC")
    q0, f0 = q.copy(), f.copy()
    first = verify(q, f, 3, 2, TABLES)
    for tables in (TABLES, block_tables()):
        assert verify(q, f, 3, 2, tables) == first
    assert np.array_equal(q, q0) and np.array_equal(f, f0)
    assert first == [Alignment(2, 3, 15, ("substitution", "swap"))]


def test_accepts_read_only_frame_slices() -> None:
    f = encode_target("CHYPEPTGGDEKARYHC")
    f.setflags(write=False)
    aln = verify(encode_query("PEPTNDEKAR"), f, 3, 1, TABLES)
    assert [(a.cost, a.diff_types, a.t_start, a.t_end) for a in aln] == [(1, ("block",), 3, 14)]
