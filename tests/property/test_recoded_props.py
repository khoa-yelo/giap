"""Decoy symmetry and completeness with recoded Sec/Pyl cells (M7-16c; ADR-0058 item 6; spec
"16c amendment (2026-10-25)" rule 21, test plan C6; statistician pre-code check Q1).

The 16a properties (``test_readthrough_props.py``) re-run on random genomes (linear or circular,
both strands, genes crossing the origin) whose annotated CDS carry random ``U``-over-TGA (and
``O``-over-TAG) residues, plus deliberate mismatches (a ``U`` over a sense codon, which rule 21 (e)
skips) and pseudogenes with ``U`` (never scanned). The expected sites are derived here from the
final genome, independently of ``gippa.mapping.readthrough``: a terminal key per scanned feature
whose last codon is a stop, a recoded key per ``U``/``O`` at ``i >= 1`` whose codon matches. The
oracle cells are every ``*`` whose codon is one of those keys (origin copies included).

- ``frame_cells(annotated_stops(..., recoded=True))`` equals the oracle cells, and the kinds are
  right (a key shared by a terminal and a recoded codon is ``recoded``);
- (a) per frame the target and every decoy have equal cell counts, every decoy cell holds ``*``;
- (b) ``mirror_cells`` is an involution for ``reverse`` and the identity for ``shuffle``;
- zero false negatives against the brute force (target and decoy search) for queries over a
  recoded cell, with every cell (terminal and recoded) consumable;
- (c) for ``reverse``, the decoy rows of ``Q`` are the target rows of ``reversed(Q)`` after
  interval mapping.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import pandas as pd
from brute_tolerant import overlaps_mod
from helpers.readthrough_toys import gene, model
from hypothesis import HealthCheck, assume, example, given, settings
from hypothesis import strategies as st
from test_readthrough_props import CFG, Plan, _query, check_search, oracle_cells, oracle_mirror

from gippa.mapping.decoys import DecoyType, build_decoys
from gippa.mapping.sixframe import STOP, Frame, SixFrameIndex, build_sixframe
from gippa.mapping.tolerant import match_tolerant
from gippa.model import Replicon
from gippa.model.genes import GeneModel
from gippa.rng import for_step

SETTINGS = settings(max_examples=40, deadline=None, suppress_health_check=list(HealthCheck))
STOPS = ("TAA", "TAG", "TGA")
SENSE = tuple(a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS)


def _rt():  # type: ignore[no-untyped-def]
    import gippa.mapping.readthrough as rt

    return rt


@dataclass(frozen=True)
class RPlan:
    seed: int
    size: int
    topology: str
    n_genes: int
    n_pseudo: int
    p_mismatch: float  # chance per gene of one extra U over a sense codon
    letter: str  # query residue over the recoded cell
    ops: tuple[str, ...]
    adjacent: bool
    kind: str


@dataclass
class RWorld:
    plan: RPlan
    genome: dict[str, str]
    gm: GeneModel
    index: SixFrameIndex
    decoys: dict[DecoyType, tuple[Frame, ...]]
    terminal_keys: set[tuple[str, str, int]]
    recoded_keys: set[tuple[str, str, int]]
    target: list[np.ndarray]  # oracle cells, every key
    recoded_cells: list[np.ndarray]  # oracle cells, recoded keys only
    mirrored: list[np.ndarray]


def _pos(strand: str, s: int, e: int, j: int, n: int) -> int:
    """Genome position (mod n) of gene-strand base ``j`` of the gene ``[s, e)``."""
    return (s + j) % n if strand == "+" else (e - 1 - j) % n


def _read(seq: list[str], strand: str, s: int, e: int, i: int, n: int) -> str:
    bases = "".join(seq[_pos(strand, s, e, 3 * i + k, n)] for k in range(3))
    return bases if strand == "+" else bases.translate(str.maketrans("ACGT", "TGCA"))


def _codon_start(strand: str, s: int, e: int, i: int, n: int) -> int:
    return (s + 3 * i) % n if strand == "+" else (e - 3 * i - 3) % n


def build(plan: RPlan) -> RWorld:
    rng = np.random.default_rng(plan.seed)
    n = plan.size
    circular = plan.topology == "circular"
    seq = list("".join(rng.choice(list("ACGT"), size=n)))
    specs = []  # (tag, ftype, strand, s, e, n_codons, U residues, O residues, mismatch residue)
    for k in range(plan.n_genes + plan.n_pseudo):
        c = int(rng.integers(8, 41))
        strand = "+" if rng.random() < 0.5 else "-"
        s = int(rng.integers(0, n)) if circular else int(rng.integers(0, n - 3 * c + 1))
        e = s + 3 * c
        inner = rng.permutation(np.arange(1, c - 1))
        n_u = int(rng.integers(1, 3))
        us = sorted(int(x) for x in inner[:n_u])
        os_ = sorted(int(x) for x in inner[n_u : n_u + 1]) if rng.random() < 0.3 else []
        mm = int(inner[n_u + 1]) if rng.random() < plan.p_mismatch and c > 5 else None
        codons = ["ATG"] + [str(rng.choice(SENSE)) for _ in range(c - 2)] + [str(rng.choice(STOPS))]
        for i in us:
            codons[i] = "TGA"
        for i in os_:
            codons[i] = "TAG"
        nt = "".join(codons)
        for j, b in enumerate(nt):  # later genes may overwrite earlier ones: rule 21 skips those
            pos = _pos(strand, s, e, j, n)
            seq[pos] = b if strand == "+" else b.translate(str.maketrans("ACGT", "TGCA"))
        ftype = "cds" if k < plan.n_genes else "pseudogene"
        specs.append((f"T{k}", ftype, strand, s, e, c, us, os_, mm))
    genome = {"g": "".join(seq)}
    reps = {"g": Replicon("g", n, plan.topology)}  # type: ignore[arg-type]
    genes, terminal, recoded = [], set(), set()
    for tag, ftype, strand, s, e, c, us, os_, mm in specs:
        aa = ["M"] + ["A"] * (c - 2)
        for i in us:
            aa[i] = "U"
        for i in os_:
            aa[i] = "O"
        if mm is not None:
            aa[mm] = "U"  # over whatever codon is there; a mismatch unless it reads TGA
        genes.append(replace(gene(tag, "g", s, e, strand, ftype), aa="".join(aa)))
        if ftype != "cds" or _read(seq, strand, s, e, c - 1, n) not in STOPS:
            continue
        terminal.add(("g", strand, _codon_start(strand, s, e, c - 1, n)))
        for i, letter in enumerate(aa):
            want = {"U": "TGA", "O": "TAG"}.get(letter)
            if i >= 1 and want and _read(seq, strand, s, e, i, n) == want:
                recoded.add(("g", strand, _codon_start(strand, s, e, i, n)))
    gm = model(reps, genes)
    index = build_sixframe(genome, reps)
    kind: DecoyType = plan.kind  # type: ignore[assignment]
    decoys = {kind: build_decoys(index, kind, 50, for_step("sixframe_index", plan.seed))}
    target = oracle_cells(index, terminal | recoded)
    return RWorld(
        plan,
        genome,
        gm,
        index,
        decoys,
        terminal,
        recoded,
        target,
        oracle_cells(index, recoded),
        oracle_mirror(kind, decoys[kind], target),
    )


def impl(w: RWorld):  # type: ignore[no-untyped-def]
    rt = _rt()
    sites, _ = rt.annotated_stops(w.gm, w.genome, recoded=True)
    target = rt.frame_cells(w.index, sites)
    decoys = {k: rt.mirror_cells(k, fr, target) for k, fr in w.decoys.items()}
    return rt.StopCells(sites=sites, target=target, decoys=decoys)


plans = st.builds(
    RPlan,
    seed=st.integers(0, 2**31 - 1),
    size=st.sampled_from([240, 450, 900, 1800]),
    topology=st.sampled_from(["linear", "circular"]),
    n_genes=st.integers(1, 8),
    n_pseudo=st.integers(0, 2),
    p_mismatch=st.sampled_from([0.0, 0.5]),
    letter=st.sampled_from(["U", "W", "Q", "K", "X"]),
    ops=st.lists(st.sampled_from(["substitution", "swap", "block", "near_block"]), max_size=1).map(
        tuple
    ),
    adjacent=st.booleans(),
    kind=st.sampled_from(["reverse", "shuffle"]),
)

_EX = [
    RPlan(21, 900, "linear", 4, 1, 0.5, "U", (), False, "reverse"),
    RPlan(22, 450, "circular", 5, 1, 0.0, "U", ("substitution",), True, "shuffle"),
    RPlan(23, 240, "circular", 6, 0, 0.5, "W", ("swap",), True, "reverse"),
    RPlan(24, 1800, "linear", 8, 2, 0.5, "U", ("block",), True, "shuffle"),
]


def _with_examples(f):  # type: ignore[no-untyped-def]
    for p in _EX:
        f = example(plan=p)(f)
    return f


def _q(w: RWorld) -> str:
    p = w.plan
    qp = Plan(
        p.seed, p.size, p.topology, p.n_genes, "marked", p.ops, p.adjacent, p.letter, p.kind, False
    )
    return _query(qp, w.index, w.recoded_cells, np.random.default_rng(p.seed + 1))


def _lists(cells) -> list[list[int]]:  # type: ignore[no-untyped-def]
    return [[int(x) for x in c] for c in cells]


# ---------------------------------------------------------------- the cells and their kinds


@SETTINGS
@given(plan=plans)
@_with_examples
def test_recoded_frame_cells_equal_the_independent_cells(plan: RPlan) -> None:
    w = build(plan)
    c = impl(w)
    assert _lists(c.target) == _lists(w.target)
    keys = {(s.contig, s.strand, s.codon_start): s.kind for s in c.sites}
    assert set(keys) == w.terminal_keys | w.recoded_keys
    assert {k for k, v in keys.items() if v == "recoded"} == w.recoded_keys
    assert len(c.sites) == len(keys)  # merged: one site per key
    for s in c.sites:
        assert s.codon in (("TGA", "TAG") if s.kind == "recoded" else STOPS)


@SETTINGS
@given(plan=plans)
@_with_examples
def test_recoded_decoy_counts_equal_hold_stops_and_mirror(plan: RPlan) -> None:
    """(a) and (b)."""
    w = build(plan)
    rt = _rt()
    c = impl(w)
    for kind in ("reverse", "shuffle"):
        frames = build_decoys(w.index, kind, 50, for_step("sixframe_index", plan.seed))
        dc = rt.mirror_cells(kind, frames, c.target)
        assert [len(x) for x in dc] == [len(x) for x in c.target]
        assert all(np.all(f.aa[x] == STOP) for f, x in zip(frames, dc, strict=True))
        assert _lists(dc) == _lists(oracle_mirror(kind, frames, w.target))
        back = rt.mirror_cells(kind, w.index.frames if kind == "reverse" else frames, dc)
        assert _lists(back) == _lists(c.target)
        if kind == "shuffle":
            assert _lists(dc) == _lists(c.target)


# ---------------------------------------------------------------- completeness


def _search(df: pd.DataFrame, search: str) -> pd.DataFrame:
    return df[df["search"] == search].reset_index(drop=True)


@SETTINGS
@given(plan=plans)
@_with_examples
def test_recoded_queries_equal_the_brute_force(plan: RPlan) -> None:
    w = build(plan)
    assume(any(len(x) for x in w.recoded_cells))
    pep = _q(w)
    assume(10 <= len(pep) <= 40)
    kind: DecoyType = plan.kind  # type: ignore[assignment]
    df = match_tolerant([pep], w.index, w.decoys, {pep: kind}, {pep}, CFG, cells=impl(w))
    check_search(_search(df, "target"), pep, w.index, w.index.frames, w.target, None)
    check_search(_search(df, "decoy"), pep, w.index, w.decoys[kind], w.mirrored, kind)


@SETTINGS
@given(plan=plans)
@_with_examples
def test_recoded_reverse_decoy_rows_of_q_are_target_rows_of_reversed_q(plan: RPlan) -> None:
    """(c) after interval mapping."""
    w = build(replace(plan, kind="reverse"))
    assume(any(len(x) for x in w.recoded_cells))
    q = _q(w)
    assume(10 <= len(q) <= 40)
    rq = q[::-1]
    df = match_tolerant(
        [q, rq], w.index, w.decoys, {q: "reverse", rq: "reverse"}, {rq}, CFG, cells=impl(w)
    )
    dec = df[(df["peptide"] == q) & (df["search"] == "decoy")]
    tgt = df[(df["peptide"] == rq) & (df["search"] == "target")]
    n = plan.size
    assert set(dec["n_diff"]) == set(tgt["n_diff"])

    def covered(xs: pd.DataFrame, ys: pd.DataFrame) -> bool:
        return all(
            any(
                x.strand == y.strand
                and overlaps_mod(
                    int(x.start) % n,
                    int(x.end) - int(x.start),
                    int(y.start) % n,
                    int(y.end) - int(y.start),
                    n,
                )
                for y in ys.itertuples(index=False)
            )
            for x in xs.itertuples(index=False)
        )

    assert covered(dec, tgt) and covered(tgt, dec)


def test_the_examples_carry_recoded_cells_and_skips() -> None:
    """Non-vacuity: every @example world has recoded cells; the mismatch worlds have a skip."""
    rt = _rt()
    for p in _EX:
        w = build(p)
        assert sum(len(x) for x in w.recoded_cells) >= 1, p
        if p.p_mismatch:
            got = rt.recoded_skips(w.gm, w.genome)
            assert set(got) == {"length_mismatch", "codon_mismatch"}
            assert got["length_mismatch"] == 0
