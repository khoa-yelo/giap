"""Properties of the annotated-stop rule (M7-16a; ADR-0058 items 1-2; spec test plan, statistician
review changes 7 and 9; GIAP_SPEC.md §11.2 "zero false negatives against a brute-force matcher
on random genomes <= 50 kb, for every difference type", now **with** the consumable-stop rule).

Random genomes (linear or circular, both strands, occasional origin-near stops) carry random
annotated CDS ending at random stop codons. The consumable cells are derived here independently
of ``giap.mapping.readthrough`` (every ``*`` whose ``to_interval`` codon is an annotated stop,
origin copies included) and must equal ``frame_cells``. Queries plant 0-2 differences
(substitution, swap, block, near-block; often next to the cell) into windows spanning a marked
cell, take windows spanning an unmarked stop, or are random. ``brute_readthrough.py`` scans every
window of every frame with the cells:

- zero false negatives / soundness (exact, no tolerance): the minimal cost equals the brute force,
  every brute-force cluster is covered by a reported locus, every reported window costs the
  minimum, and no reported window holds an unmarked ``*``; for the target and the decoy search
  (decoy cells mirrored by the oracle: ``n - 1 - p`` / ``p``);
- decoy symmetry: (a) per-frame cell counts equal, every decoy cell holds ``*``; (b)
  ``mirror_cells`` is an involution for ``reverse``, the identity for ``shuffle``; (c) the target
  search of ``Q`` on the decoy frames (cells = the mirrored cells) equals the decoy search of
  ``Q``, and for ``reverse`` the target rows of ``reversed(Q)`` equal the decoy rows of ``Q``
  (cost and locus); (d) example 5's broken build (decoy cells removed) breaks (c); (e) the
  frameshift route: the reverse-decoy rows of ``reversed(Q)`` are the target rows of ``Q`` with
  cells;
- no-op: empty cells give the frames of the call without the keyword (the M3-07 and M3-08
  property strategies).
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations

import numpy as np
import pandas as pd
import pytest
import test_frameshift_props as fsp
from brute_readthrough import brute_force_rt, masks, window_cost_rt
from brute_tolerant import overlaps_mod, window_interval
from helpers.readthrough_toys import gene, model, rc, toy_c1
from hypothesis import HealthCheck, assume, example, given, settings
from hypothesis import strategies as st
from test_tolerant_seeds_props import worlds as seed_worlds

from giap.config.model import TolerantConfig
from giap.mapping.decoys import DecoyType, build_decoys
from giap.mapping.frameshift import match_frameshift
from giap.mapping.sixframe import STOP, Frame, SixFrameIndex, build_sixframe, decode, encode_query
from giap.mapping.tolerant import match_tolerant
from giap.model import Replicon
from giap.model.genes import GeneModel
from giap.rng import for_step

CFG = TolerantConfig()
SETTINGS = settings(max_examples=60, deadline=None, suppress_health_check=list(HealthCheck))
RESIDUES = "ACDEFGHKLMNPQRSTVWY"
DIFF_VOCAB = {"substitution", "swap", "block", "near_block"}
_BLOCK_RAW = (("N", "GG"), ("Q", "AG"), ("GT", "AS"), ("AV", "GL"), ("AN", "GQ"))
_NEAR_RAW = (("W", "AD"), ("W", "EG"), ("R", "GV"))


def _expand(raw: tuple[tuple[str, str], ...]) -> list[tuple[str, str]]:
    out = []
    for a, b in raw:
        for x in sorted({"".join(p) for p in permutations(a)}):
            for y in sorted({"".join(p) for p in permutations(b)}):
                out += [(x, y), (y, x)]
    return out


SEGS = {"block": _expand(_BLOCK_RAW), "near_block": _expand(_NEAR_RAW)}  # (target, query)


def _k(L: int) -> int:
    return 1 if L <= 11 else 2


def _rt():  # type: ignore[no-untyped-def]
    import giap.mapping.readthrough as rt

    return rt


# ---------------------------------------------------------------- worlds


@dataclass(frozen=True)
class Plan:
    seed: int
    size: int
    topology: str
    n_genes: int
    query: str  # "marked" | "unmarked" | "random"
    ops: tuple[str, ...]
    adjacent: bool  # plant the first op next to the cell
    letter: str  # query residue over the cell
    kind: str
    origin: bool  # plant a stop codon at or across the origin


@dataclass
class World:
    plan: Plan
    genome: dict[str, str]
    reps: dict[str, Replicon]
    gm: GeneModel
    index: SixFrameIndex
    decoys: dict[DecoyType, tuple[Frame, ...]]
    target: list[np.ndarray]  # oracle cells per frame (positions)
    mirrored: list[np.ndarray]
    query: str


def oracle_cells(index: SixFrameIndex, keys: set[tuple[str, str, int]]) -> list[np.ndarray]:
    """Every ``*`` whose codon (``to_interval``, start mod length) is an annotated stop."""
    out = []
    for fid, f in enumerate(index.frames):
        n = index.replicons[f.contig].length
        ps = []
        for p in np.nonzero(f.aa == STOP)[0]:
            iv = index.to_interval(fid, int(p), int(p) + 1)
            if (iv.contig, iv.strand, iv.start % n) in keys:
                ps.append(int(p))
        out.append(np.array(ps, dtype=np.int64))
    return out


def oracle_mirror(
    kind: str, frames: tuple[Frame, ...], cells: list[np.ndarray]
) -> list[np.ndarray]:
    if kind == "shuffle":
        return [c.copy() for c in cells]
    return [np.sort(len(f.aa) - 1 - c).astype(np.int64) for f, c in zip(frames, cells, strict=True)]


def _stops(index: SixFrameIndex) -> list[tuple[int, int, str, int]]:
    """(frame_id, p, strand, codon start mod n) of every '*'."""
    out = []
    for fid, f in enumerate(index.frames):
        n = index.replicons[f.contig].length
        for p in np.nonzero(f.aa == STOP)[0]:
            iv = index.to_interval(fid, int(p), int(p) + 1)
            out.append((fid, int(p), iv.strand, iv.start % n))
    return out


def _gene_for(tag: str, strand: str, cs: int, n: int, circular: bool, codons: int):  # type: ignore[no-untyped-def]
    if strand == "+":
        end = cs + 3
        start = end - 3 * codons
        if start < 0:
            if circular:
                start, end = start + n, end + n
            else:
                start = end % 3
        if circular and start >= n:
            start, end = start - n, end - n
        return gene(tag, "g", start, end, "+")
    end = cs + 3 * codons
    if end > n and not circular:
        end = cs + 3 * ((n - cs) // 3)
    return gene(tag, "g", cs, end, "-")


def build_world(plan: Plan) -> World:
    rng = np.random.default_rng(plan.seed)
    n = plan.size
    seq = list("".join(rng.choice(list("ACGT"), size=n)))
    circular = plan.topology == "circular"
    if plan.origin:
        at = int(rng.choice([n - 2, n - 1, 0, 1, 2, 3, 30]))
        codon = str(rng.choice(["TAA", "TAG", "TGA"]))
        if rng.random() < 0.5:
            codon = rc(codon)
        for i, c in enumerate(codon):
            if circular or at + i < n:
                seq[(at + i) % n] = c
    genome = {"g": "".join(seq)}
    reps = {"g": Replicon("g", n, plan.topology)}  # type: ignore[arg-type]
    index = build_sixframe(genome, reps)
    stops = sorted({(s, c) for _, _, s, c in _stops(index)})
    if plan.origin:
        near = [x for x in stops if x[1] >= n - 2 or x[1] <= 3][:2]
        rest = [x for x in stops if x not in near]
        pick = near + [rest[int(i)] for i in rng.permutation(len(rest))[: plan.n_genes - len(near)]]
    else:
        idx = rng.permutation(len(stops))[: plan.n_genes]
        pick = [stops[int(i)] for i in idx]
    genes = [
        _gene_for(f"T{i}", s, c, n, circular, int(rng.integers(5, 40)))
        for i, (s, c) in enumerate(pick)
    ]
    gm = model(reps, genes)
    keys = {("g", s, c) for s, c in pick}
    target = oracle_cells(index, keys)
    kind: DecoyType = plan.kind  # type: ignore[assignment]
    decoys = {kind: build_decoys(index, kind, 50, for_step("sixframe_index", plan.seed))}
    mirrored = oracle_mirror(kind, decoys[kind], target)
    query = _query(plan, index, target, rng)
    return World(plan, genome, reps, gm, index, decoys, target, mirrored, query)


def _plant(
    q: list[str], cell: int | None, op: str, near: bool, rng: np.random.Generator
) -> list[str]:
    L = len(q)
    spots = [j for j in range(L - 1) if cell is None or cell not in (j, j + 1)]
    if near and cell is not None:
        spots = [j for j in (cell - 2, cell + 1) if 0 <= j < L - 1] or spots
    if not spots:
        return q
    if op == "substitution":
        j = int(rng.choice(spots))
        q[j] = (
            RESIDUES[(RESIDUES.index(q[j]) + 1 + int(rng.integers(18))) % 19]
            if q[j] in RESIDUES
            else "A"
        )
        return q
    if op == "swap":
        for j in rng.permutation(spots):
            j = int(j)
            if q[j] != q[j + 1]:
                q[j], q[j + 1] = q[j + 1], q[j]
                return q
        return q
    for j in rng.permutation(spots):  # block / near-block: a target segment already in q
        j = int(j)
        for t_seg, q_seg in SEGS[op]:
            if "".join(q[j : j + len(t_seg)]) == t_seg and (
                cell is None or not j <= cell < j + len(t_seg)
            ):
                return q[:j] + list(q_seg) + q[j + len(t_seg) :]
    # none present: write one in (target and query differ by one exchange relative to the frame)
    t_seg, q_seg = SEGS[op][int(rng.integers(len(SEGS[op])))]
    j = int(rng.choice(spots))
    return q[:j] + list(q_seg) + q[j + 1 :]


def _query(
    plan: Plan, index: SixFrameIndex, target: list[np.ndarray], rng: np.random.Generator
) -> str:
    if plan.query == "random":
        return "".join(rng.choice(list(RESIDUES), size=int(rng.integers(10, 21))))
    marked = {(fid, int(p)) for fid, c in enumerate(target) for p in c}
    pool = (
        sorted(marked)
        if plan.query == "marked"
        else sorted((fid, p) for fid, p, _, _ in _stops(index) if (fid, p) not in marked)
    )
    if not pool:
        return "".join(rng.choice(list(RESIDUES), size=12))
    fid, p = pool[int(rng.integers(len(pool)))]
    aa = decode(index.frames[fid].aa)
    L = min(int(rng.integers(10, 21)), len(aa))
    lo = max(0, p - L + 1)
    hi = max(lo + 1, min(p, len(aa) - L) + 1)
    best = None
    for _ in range(12):
        a = int(rng.integers(lo, hi))
        w = aa[a : a + L]
        others = sum(ch in "*X" for i, ch in enumerate(w) if a + i != p)
        if best is None or others < best[0]:
            best = (others, a, w)
    _, a, w = best  # type: ignore[misc]
    q = [ch if ch not in "*X" else str(rng.choice(list(RESIDUES))) for ch in w]
    cell = p - a if 0 <= p - a < len(q) else None
    if cell is not None:
        q[cell] = plan.letter
    for i, op in enumerate(plan.ops):
        q = _plant(q, cell if cell is not None and cell < len(q) else None, op,
                   plan.adjacent and i == 0, rng)  # fmt: skip
    return "".join(q)[:40]


plans = st.builds(
    Plan,
    seed=st.integers(0, 2**31 - 1),
    size=st.sampled_from([240, 450, 900, 1800, 3000]),
    topology=st.sampled_from(["linear", "circular"]),
    n_genes=st.integers(1, 10),
    query=st.sampled_from(["marked", "marked", "marked", "unmarked", "random"]),
    ops=st.lists(st.sampled_from(sorted(DIFF_VOCAB)), max_size=2).map(tuple),
    adjacent=st.booleans(),
    letter=st.sampled_from([*RESIDUES, "U", "X"]),
    kind=st.sampled_from(["reverse", "shuffle"]),
    origin=st.sampled_from([False, False, True]),
)

# Hypothesis @examples: each difference type next to a marked cell, and the wrap
_EX = [
    Plan(11, 900, "linear", 4, "marked", ("substitution",), True, "W", "reverse", False),
    Plan(12, 900, "linear", 4, "marked", ("swap",), True, "Q", "shuffle", False),
    Plan(13, 1800, "linear", 6, "marked", ("block",), True, "W", "reverse", False),
    Plan(14, 1800, "circular", 6, "marked", ("near_block",), True, "Y", "shuffle", False),
    Plan(15, 450, "circular", 3, "marked", ("substitution",), False, "U", "reverse", True),
    Plan(16, 450, "circular", 3, "marked", (), False, "W", "shuffle", True),
    Plan(17, 900, "linear", 5, "unmarked", (), False, "W", "reverse", False),
    Plan(18, 240, "circular", 2, "marked", ("swap", "substitution"), True, "K", "reverse", True),
]


def _with_examples(f):  # type: ignore[no-untyped-def]
    for p in _EX:
        f = example(plan=p)(f)
    return f


def impl_cells(w: World):  # type: ignore[no-untyped-def]
    rt = _rt()
    sites, _ = rt.annotated_stops(w.gm, w.genome)
    target = rt.frame_cells(w.index, sites)
    decoys = {k: rt.mirror_cells(k, fr, target) for k, fr in w.decoys.items()}
    return rt.StopCells(sites=sites, target=target, decoys=decoys)


def _search(df: pd.DataFrame, search: str) -> pd.DataFrame:
    return df[df["search"] == search].reset_index(drop=True)


def check_search(
    df: pd.DataFrame,
    peptide: str,
    index: SixFrameIndex,
    frames: tuple[Frame, ...],
    cells: list[np.ndarray],
    kind: str | None,
) -> None:
    k = _k(len(peptide))
    m = masks(frames, cells)
    brute = brute_force_rt(peptide, frames, m, k)
    q = encode_query(peptide)
    if brute.cost is None:
        assert len(df) == 0, df
        return
    assert len(df) > 0, ("false negative", peptide, brute)
    assert set(df["n_diff"].tolist()) == {brute.cost}, (df, brute)
    seen = set()
    for r in df.itertuples(index=False):
        fid, a, b, n = int(r.frame_id), int(r.aa_start), int(r.aa_end), int(r.n_diff)
        assert window_cost_rt(q, frames[fid].aa, m[fid], a, b, n) == n, ("false positive", r)
        win = frames[fid].aa[a:b]
        for j in np.nonzero(win == STOP)[0]:
            assert m[fid][a + int(j)], ("an unmarked '*' consumed", r)
        types = [] if r.diff_types == "" else str(r.diff_types).split(";")
        assert len(types) == n and set(types) <= DIFF_VOCAB, r
        contig, start, length, strand, label = window_interval(index, frames, kind, fid, a, b)
        rep_len = index.replicons[contig].length
        assert (r.contig, r.strand, int(r.frame)) == (contig, strand, label), r
        assert int(r.start) % rep_len == start and int(r.end) - int(r.start) == length, r
        key = (r.contig, int(r.start) % rep_len, int(r.end), r.strand)
        assert key not in seen, ("duplicate locus", r)
        seen.add(key)
    for cluster in brute.clusters:
        covered = False
        for fid, s, e in cluster:
            try:
                contig, start, length, strand, label = window_interval(
                    index, frames, kind, fid, s, e
                )
            except ValueError:
                continue
            rep_len = index.replicons[contig].length
            for r in df.itertuples(index=False):
                if (r.contig, r.strand, int(r.frame)) != (contig, strand, label):
                    continue
                if overlaps_mod(start, length, int(r.start) % rep_len, int(r.end) - int(r.start),
                                rep_len):  # fmt: skip
                    covered = True
                    break
            if covered:
                break
        assert covered, ("false negative: missed locus", peptide, cluster, df)


# ---------------------------------------------------------------- the oracle and the cells


def test_brute_force_reproduces_the_spec_examples() -> None:
    genome, reps, _ = toy_c1()
    index = build_sixframe(genome, reps)
    cells = [np.array([18])] + [np.zeros(0, dtype=np.int64)] * 5
    m = masks(index.frames, cells)
    assert brute_force_rt("WQLNEFKWVPAGHTL", index.frames, m, 2).clusters == (((0, 11, 26),),)
    assert brute_force_rt("HTLDMRKQAAGLKE", index.frames, m, 2).cost is None
    assert brute_force_rt("NEFKWVPAGHA", index.frames, m, 1).cost is None
    assert brute_force_rt("NEFKWVPAGHAL", index.frames, m, 2).cost == 2
    plain = masks(index.frames, [np.zeros(0, dtype=np.int64)] * 6)
    assert brute_force_rt("WQLNEFKWVPAGHTL", index.frames, plain, 2).cost is None


@SETTINGS
@given(plan=plans)
@_with_examples
def test_frame_cells_equal_the_independent_cells(plan: Plan) -> None:
    w = build_world(plan)
    c = impl_cells(w)
    assert [list(x) for x in c.target] == [list(x) for x in w.target]
    assert all(x.dtype == np.int64 for x in c.target)
    for fid, x in enumerate(c.target):  # rule 2's round trip
        n = w.index.replicons[w.index.frames[fid].contig].length
        for p in x:
            assert w.index.frames[fid].aa[p] == STOP
            iv = w.index.to_interval(fid, int(p), int(p) + 1)
            assert ("g", iv.strand, iv.start % n) in _rt().site_keys(c.sites)


# ---------------------------------------------------------------- zero false negatives


@SETTINGS
@given(plan=plans)
@_with_examples
def test_match_tolerant_equals_the_brute_force_with_cells(plan: Plan) -> None:
    w = build_world(plan)
    pep = w.query
    assume(10 <= len(pep) <= 40)
    kind: DecoyType = plan.kind  # type: ignore[assignment]
    df = match_tolerant([pep], w.index, w.decoys, {pep: kind}, {pep}, CFG, cells=impl_cells(w))
    check_search(_search(df, "target"), pep, w.index, w.index.frames, w.target, None)
    check_search(_search(df, "decoy"), pep, w.index, w.decoys[kind], w.mirrored, kind)


@pytest.mark.parametrize("topology", ["linear", "circular"])
def test_fifty_kb_genomes_equal_the_brute_force(topology: str) -> None:
    """Spec size: 50 kb, 300 annotated genes, 16 queries over marked cells (0-1 more difference)."""
    plan = Plan(50 if topology == "linear" else 51, 50_000, topology, 300, "marked", (), False,
                "W", "reverse", topology == "circular")  # fmt: skip
    w = build_world(plan)
    cells = impl_cells(w)
    queries = []
    for i in range(16):
        p2 = Plan(plan.seed * 100 + i, plan.size, topology, plan.n_genes, "marked",
                  ("substitution",) if i % 2 else (), bool(i % 3), "WQKU"[i % 4], "reverse",
                  plan.origin)  # fmt: skip
        queries.append(_query(p2, w.index, w.target, np.random.default_rng(p2.seed)))
    queries = sorted(set(q for q in queries if 10 <= len(q) <= 40))
    assert len(queries) >= 12
    df = match_tolerant(queries, w.index, w.decoys, dict.fromkeys(queries, "reverse"),
                        set(queries), CFG, cells=cells)  # fmt: skip
    found = 0
    for q in queries:
        rows = df[df["peptide"] == q]
        tgt = _search(rows, "target")
        check_search(tgt, q, w.index, w.index.frames, w.target, None)
        check_search(_search(rows, "decoy"), q, w.index, w.decoys["reverse"], w.mirrored,
                     "reverse")  # fmt: skip
        found += int(len(tgt) > 0)
    assert found >= len(queries) // 2  # most planted read-through windows are placed


# ---------------------------------------------------------------- decoy symmetry


@SETTINGS
@given(plan=plans)
@_with_examples
def test_decoy_cell_counts_equal_and_hold_stops(plan: Plan) -> None:
    """(a) and (b)."""
    w = build_world(plan)
    rt = _rt()
    c = impl_cells(w)
    for kind in ("reverse", "shuffle"):
        frames = build_decoys(w.index, kind, 50, for_step("sixframe_index", plan.seed))
        dc = rt.mirror_cells(kind, frames, c.target)
        assert [len(x) for x in dc] == [len(x) for x in c.target]
        assert all(np.all(f.aa[x] == STOP) for f, x in zip(frames, dc, strict=True))
        assert [list(x) for x in dc] == [list(x) for x in oracle_mirror(kind, frames, w.target)]
        back = rt.mirror_cells(kind, w.index.frames if kind == "reverse" else frames, dc)
        assert [list(x) for x in back] == [list(x) for x in c.target]


def _cols(df: pd.DataFrame) -> pd.DataFrame:
    cols = ["frame_id", "aa_start", "aa_end", "n_diff", "diff_types"]
    return df[cols].sort_values(cols).reset_index(drop=True)


def decoy_equals_target_on_decoy_frames(w: World, pep: str, cells) -> bool:  # type: ignore[no-untyped-def]
    """(c): the decoy search of ``pep`` equals the target search of ``pep`` on the decoy frames
    whose cells are the decoy cells."""
    rt = _rt()
    kind: DecoyType = w.plan.kind  # type: ignore[assignment]
    a = match_tolerant([pep], w.index, w.decoys, {pep: kind}, set(), CFG, cells=cells)
    flipped = SixFrameIndex(replicons=w.index.replicons, frames=w.decoys[kind])
    swapped = rt.StopCells(sites=cells.sites, target=cells.decoys[kind],
                           decoys={kind: cells.target})  # fmt: skip
    b = match_tolerant([pep], flipped, {kind: w.index.frames}, {pep: kind}, {pep}, CFG,
                       cells=swapped)  # fmt: skip
    return _cols(_search(a, "decoy")).equals(_cols(_search(b, "target")))


@SETTINGS
@given(plan=plans)
@_with_examples
def test_target_search_on_decoy_frames_equals_decoy_search_with_cells(plan: Plan) -> None:
    w = build_world(plan)
    assume(plan.topology == "linear")  # no origin extension (as the M3-07 analogue)
    assume(10 <= len(w.query) <= 40)
    assert decoy_equals_target_on_decoy_frames(w, w.query, impl_cells(w))


@SETTINGS
@given(plan=plans)
@_with_examples
def test_reverse_decoy_rows_of_q_are_target_rows_of_reversed_q(plan: Plan) -> None:
    w = build_world(Plan(**{**plan.__dict__, "kind": "reverse"}))
    q = w.query
    assume(10 <= len(q) <= 40)
    rq = q[::-1]
    df = match_tolerant([q, rq], w.index, w.decoys, {q: "reverse", rq: "reverse"}, {rq}, CFG,
                        cells=impl_cells(w))  # fmt: skip
    dec = df[(df["peptide"] == q) & (df["search"] == "decoy")]
    tgt = df[(df["peptide"] == rq) & (df["search"] == "target")]
    n = plan.size
    assert set(dec["n_diff"]) == set(tgt["n_diff"])

    def covered(xs: pd.DataFrame, ys: pd.DataFrame) -> bool:
        return all(
            any(x.strand == y.strand and overlaps_mod(int(x.start) % n, int(x.end) - int(x.start),
                                                      int(y.start) % n, int(y.end) - int(y.start),
                                                      n)
                for y in ys.itertuples(index=False))
            for x in xs.itertuples(index=False)
        )  # fmt: skip

    assert covered(dec, tgt) and covered(tgt, dec)


def test_example_5_decoy_row_and_the_broken_build_breaks_symmetry() -> None:
    """(d): spec example 5 through ``match_tolerant``; removing the decoy cells loses the decoy
    row, which property (c) detects."""
    genome, reps, gm = toy_c1()
    index = build_sixframe(genome, reps)
    decoys = {"reverse": build_decoys(index, "reverse", 50, for_step("sixframe_index", 1))}
    plan = Plan(0, 126, "linear", 1, "marked", (), False, "W", "reverse", False)
    w = World(plan, genome, reps, gm, index, decoys,  # type: ignore[arg-type]
              oracle_cells(index, {("c1", "+", 54)}), [], "LTHGAPVWKFENLQW")  # fmt: skip
    rt = _rt()
    good = impl_cells(w)
    pep = "LTHGAPVWKFENLQW"
    df = match_tolerant([pep], index, decoys, {pep: "reverse"}, {pep}, CFG, cells=good)
    assert len(_search(df, "target")) == 0
    d = _search(df, "decoy")
    assert len(d) == 1
    r = d.iloc[0]
    assert (r["aa_start"], r["aa_end"], r["n_diff"], r["start"], r["end"], r["strand"]) == (
        16, 31, 1, 33, 78, "+",
    )  # fmt: skip
    assert decoy_equals_target_on_decoy_frames(w, pep, good)
    none = tuple(np.zeros(0, dtype=np.int64) for _ in range(6))
    broken = rt.StopCells(sites=good.sites, target=good.target, decoys={"reverse": none})
    df_b = match_tolerant([pep], index, decoys, {pep: "reverse"}, {pep}, CFG, cells=broken)
    assert len(_search(df_b, "decoy")) == 0
    flipped = SixFrameIndex(replicons=index.replicons, frames=decoys["reverse"])
    swapped = rt.StopCells(sites=good.sites, target=good.decoys["reverse"],
                           decoys={"reverse": good.target})  # fmt: skip
    b = match_tolerant([pep], flipped, {"reverse": index.frames}, {pep: "reverse"}, {pep}, CFG,
                       cells=swapped)  # fmt: skip
    assert not _cols(_search(df_b, "decoy")).equals(_cols(_search(b, "target")))


@SETTINGS
@given(plan=plans)
@_with_examples
def test_frameshift_reverse_decoys_of_reversed_query_are_the_target_rows_with_cells(
    plan: Plan,
) -> None:
    """(e), extending M3-08's ``test_reverse_decoys_of_reversed_query_are_the_target_rows``: with
    the cells, the reverse-decoy split rows of ``reversed(Q)`` are the target split rows of
    ``Q`` (the decoy search uses the reverse type's mirrored cells, so precedence suppresses the
    same splits)."""
    w = build_world(Plan(**{**plan.__dict__, "kind": "reverse"}))
    q = w.query
    assume(10 <= len(q) <= 40 and q[::-1] != q)
    rq = q[::-1]
    df = match_frameshift([q, rq], w.index, w.decoys, {q: "reverse", rq: "reverse"}, {q}, CFG,
                          cells=impl_cells(w))  # fmt: skip
    tgt = [r[:4] + r[5:] for r in fsp._rows(df[(df["peptide"] == q) & (df["search"] == "target")])]
    dec = [r[:4] + r[5:] for r in fsp._rows(df[(df["peptide"] == rq) & (df["search"] == "decoy")])]
    assert {r[4] for r in tgt} == {r[4] for r in dec}
    assert fsp._covers(tgt, dec, plan.size) and fsp._covers(dec, tgt, plan.size)


# ---------------------------------------------------------------- no-op


def _empty_cells(index: SixFrameIndex, decoys) -> object:  # type: ignore[no-untyped-def]
    rt = _rt()
    none = tuple(np.zeros(0, dtype=np.int64) for _ in index.frames)
    return rt.StopCells(sites=(), target=none, decoys=dict.fromkeys(decoys, none))


@settings(max_examples=60, deadline=None, suppress_health_check=list(HealthCheck))
@given(w=seed_worlds())
def test_empty_cells_are_a_no_op_for_match_tolerant(w) -> None:  # type: ignore[no-untyped-def]
    pep = w.plant.query
    a = match_tolerant([pep], w.index, w.decoys, {pep: w.kind}, {pep}, CFG)
    for cells in (None, _empty_cells(w.index, w.decoys)):
        b = match_tolerant([pep], w.index, w.decoys, {pep: w.kind}, {pep}, CFG, cells=cells)
        pd.testing.assert_frame_equal(a, b)


@settings(max_examples=40, deadline=None, suppress_health_check=list(HealthCheck))
@given(case=fsp.cases())
def test_empty_cells_are_a_no_op_for_match_frameshift(case) -> None:  # type: ignore[no-untyped-def]
    index = fsp._index({"g": case.genome}, case.topology)
    decoys = {case.kind: build_decoys(index, case.kind, 50, for_step("sixframe_index", case.seed))}
    q = case.query
    cfg = fsp._cfg(case.min_side)
    tq = {q} if case.in_target else set()
    a = match_frameshift([q], index, decoys, {q: case.kind}, tq, cfg)
    for cells in (None, _empty_cells(index, decoys)):
        b = match_frameshift([q], index, decoys, {q: case.kind}, tq, cfg, cells=cells)
        pd.testing.assert_frame_equal(a, b)


def test_empty_gene_model_is_a_no_op_on_a_world() -> None:
    plan = Plan(3, 1800, "circular", 1, "marked", (), False, "W", "shuffle", True)
    w = build_world(plan)
    none = [np.zeros(0, dtype=np.int64)] * 6
    w0 = World(plan, w.genome, w.reps, model(w.reps, []), w.index, w.decoys, none, none, w.query)
    c = impl_cells(w0)
    assert c.sites == () and all(len(x) == 0 for x in c.target)
    pep = w.query
    a = match_tolerant([pep], w.index, w.decoys, {pep: "shuffle"}, {pep}, CFG)
    b = match_tolerant([pep], w.index, w.decoys, {pep: "shuffle"}, {pep}, CFG, cells=c)
    pd.testing.assert_frame_equal(a, b)
