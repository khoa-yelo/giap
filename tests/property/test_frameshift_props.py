"""Completeness and symmetry of the ±1-nt frameshift split matcher (issue M3-08; §7.2.4, §11.2).

"Property tests MUST show zero false negatives against a brute-force matcher on random genomes
<= 50 kb, for every difference type." Random genomes (linear or circular, rotated so that the
planted split may cross the origin, occasional ``N``) carry a split built from two adjacent
frames (``δ = ±1``, either strand); the query is the two sides with 0..k-1 further planted
operations (substitution, swap, every order of every ``BLOCKS`` / ``NEAR_BLOCKS`` pair) on
either side. The brute force (tests/helpers/oracle_frameshift.py) scans every split point, every
junction of every frame pair and every window, for the target and the decoy search:

- the reported cost equals the brute-force minimum (nothing when there is none);
- every brute-force locus at that cost is overlapped by a reported row (zero false negatives);
- every reported row is a brute-force locus: interval, blocks, cost, diff types and prefix frame
  (zero false positives); no interval is reported twice;
- strand symmetry: the reverse complement genome gives the mirrored rows (example 5 generalised);
- reverse decoys: the reverse-decoy rows of reversed(Q) are the target rows of Q;
- thread count 1 vs 4 gives identical frames.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
from dataclasses import dataclass
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))

from fs_genomes import SENSE_RESIDUES, background, codons_for, join_split, rc
from oracle_frameshift import TYPE_ORDER, Result, frameshift_oracle, overlaps_mod

import giap
from giap.config.model import TolerantConfig
from giap.mapping.decoys import DecoyType, build_decoys
from giap.mapping.frameshift import match_frameshift
from giap.mapping.sixframe import SixFrameIndex, build_sixframe
from giap.model import Replicon
from giap.rng import for_step

SETTINGS = settings(max_examples=60, deadline=None, suppress_health_check=list(HealthCheck))
RESIDUES = "ACDEFGHKLMNPQRSTVWY"
OPS = ("substitution", "swap", "block", "near_block")
_BLOCK_RAW = (("N", "GG"), ("Q", "AG"), ("GT", "AS"), ("AV", "GL"), ("AN", "GQ"))
_NEAR_RAW = (("W", "AD"), ("W", "EG"), ("R", "GV"))


def _expand(raw: tuple[tuple[str, str], ...]) -> list[tuple[str, str]]:
    out = []
    for a, b in raw:
        for x in sorted({"".join(p) for p in permutations(a)}):
            for y in sorted({"".join(p) for p in permutations(b)}):
                out += [(x, y), (y, x)]
    return out


BLOCK_SEGS = _expand(_BLOCK_RAW)  # (genome segment, query segment)
NEAR_SEGS = _expand(_NEAR_RAW)


def _k(L: int) -> int:
    return 1 if L <= 11 else 2


def _cfg(min_side: int) -> TolerantConfig:
    return TolerantConfig(frameshift={"enabled": True, "min_side_aa": min_side})  # type: ignore[arg-type]


@dataclass(frozen=True)
class Case:
    genome: str
    topology: str
    query: str
    kind: DecoyType
    in_target: bool
    min_side: int
    seed: int


@st.composite
def op_piece(draw: st.DrawFn, kind: str) -> tuple[str, str]:
    """(genome segment, query segment) of one planted operation."""
    if kind == "substitution":
        a = draw(st.sampled_from(RESIDUES))
        return a, draw(st.sampled_from(RESIDUES.replace(a, "")))
    if kind == "swap":
        a = draw(st.sampled_from(RESIDUES))
        b = draw(st.sampled_from(RESIDUES.replace(a, "")))
        return a + b, b + a
    return draw(st.sampled_from(BLOCK_SEGS if kind == "block" else NEAR_SEGS))


@st.composite
def split_sides(draw: st.DrawFn, min_side: int) -> tuple[str, str, str, str]:
    """(genome prefix, genome suffix, query prefix, query suffix) with 0..k-1 planted ops."""
    side = st.text(alphabet=SENSE_RESIDUES, min_size=min_side, max_size=19)
    gp, gs = draw(side), draw(side)
    qp, qs = gp, gs
    if draw(st.booleans()):
        t_seg, q_seg = draw(op_piece(draw(st.sampled_from(OPS))))
        on_prefix = draw(st.booleans())
        host = gp if on_prefix else gs
        i = draw(st.integers(0, len(host)))
        g2, q2 = host[:i] + t_seg + host[i:], host[:i] + q_seg + host[i:]
        if on_prefix:
            gp, qp = g2, q2
        else:
            gs, qs = g2, q2
        assume(len(qp) + len(qs) >= 12)  # one more difference needs k = 2
    assume(len(qp) >= min_side and len(qs) >= min_side and 10 <= len(qp) + len(qs) <= 40)
    return gp, gs, qp, qs


@st.composite
def cases(draw: st.DrawFn) -> Case:
    topology = draw(st.sampled_from(["linear", "circular"]))
    size = draw(st.sampled_from([240, 600, 1500, 4000]))
    seed = draw(st.integers(0, 2**31 - 1))
    min_side = draw(st.sampled_from([4, 4, 4, 3, 5]))
    kind = draw(st.sampled_from(["reverse", "shuffle"]))
    in_target = draw(st.sampled_from([True, True, False]))
    rng = np.random.default_rng(seed)
    genome = background(size, seed, n_rate=draw(st.sampled_from([0.0, 0.0, 0.003])))
    if draw(st.integers(0, 9)) == 7:  # about 1 in 10: an unplanted random query
        q = "".join(draw(st.lists(st.sampled_from(RESIDUES), min_size=10, max_size=30)))
        return Case(genome, topology, q, kind, in_target, min_side, seed)  # type: ignore[arg-type]
    gp, gs, qp, qs = draw(split_sides(min_side))
    delta = draw(st.sampled_from([1, -1]))
    nt = join_split(codons_for(gp, rng), codons_for(gs, rng), delta, extra="ACGT"[seed % 4])
    assume(nt is not None)
    assert nt is not None
    if draw(st.booleans()):
        nt = rc(nt)  # planted on the minus strand
    pos = draw(st.integers(0, size - len(nt)))
    genome = genome[:pos] + nt + genome[pos + len(nt) :]
    if topology == "circular":
        cut = draw(st.one_of(st.integers(-12, len(nt) + 12), st.integers(0, size - 1)))
        r = (pos + cut) % size
        genome = genome[r:] + genome[:r]
    return Case(genome, topology, qp + qs, kind, in_target, min_side, seed)  # type: ignore[arg-type]


def _index(genome: dict[str, str], topology: str) -> SixFrameIndex:
    reps = {c: Replicon(c, len(s), topology) for c, s in genome.items()}  # type: ignore[arg-type]
    return build_sixframe(genome, reps)


def _rows(df: pd.DataFrame) -> list[tuple[object, ...]]:
    cols = ["contig", "start", "end", "strand", "frame", "n_diff", "diff_types", "blocks"]
    return sorted(tuple(r) for r in df[cols].itertuples(index=False, name=None))


def precedence(n_diff: object, types: object) -> tuple[int, tuple[int, ...]]:
    """(cost, diff types ranked substitution < swap < block < near_block < frameshift)."""
    return int(n_diff), tuple(sorted(TYPE_ORDER.index(t) for t in str(types).split(";")))  # type: ignore[call-overload]


def check_search(df: pd.DataFrame, res: Result, index: SixFrameIndex) -> None:
    """Rows of one peptide and search agree with the brute force (see module docstring)."""
    if res.cost is None:
        assert df.empty, df
        return
    assert not df.empty, res
    rows = _rows(df)
    assert {r[5] for r in rows} == {res.cost}
    groups: dict[tuple[object, ...], list[tuple[object, ...]]] = {}
    for p in res.placements:
        groups.setdefault((p.contig, p.start, p.end, p.strand), []).append(
            (p.contig, p.start, p.end, p.strand, p.frame, p.cost, p.diff_types, p.blocks_str)
        )
    keys = [r[:4] for r in rows]
    assert len(keys) == len(set(keys)), "an interval is reported twice"
    for r in rows:
        assert r[:4] in groups, (r, sorted(groups))  # zero false positives
        options = groups[r[:4]]
        assert precedence(r[5], r[6]) == min(precedence(o[5], o[6]) for o in options), (r, options)
        assert r in options, (r, options)
    for key in groups:  # zero false negatives: every brute-force locus is covered
        contig, start, end, strand = key
        length = index.replicons[str(contig)].length
        assert any(
            r[0] == contig
            and r[3] == strand
            and overlaps_mod(
                int(start), int(end) - int(start), int(r[1]), int(r[2]) - int(r[1]), length
            )
            for r in rows
        ), (key, rows)


def _search(df: pd.DataFrame, peptide: str, search: str) -> pd.DataFrame:
    return df[(df["peptide"] == peptide) & (df["search"] == search)]


@SETTINGS
@given(cases())
def test_split_matcher_equals_brute_force(case: Case) -> None:
    index = _index({"g": case.genome}, case.topology)
    decoys = {case.kind: build_decoys(index, case.kind, 50, for_step("sixframe_index", case.seed))}
    q = case.query
    targets = {q} if case.in_target else set()
    df = match_frameshift([q], index, decoys, {q: case.kind}, targets, _cfg(case.min_side))
    k = _k(len(q))
    assert set(df["route"]) <= {"tolerant"}
    if case.in_target:
        res = frameshift_oracle(q, index, index.frames, None, k, case.min_side)
        check_search(_search(df, q, "target"), res, index)
    else:
        assert _search(df, q, "target").empty
    dec = _search(df, q, "decoy")
    assert set(dec["decoy_type"]) <= {case.kind}
    res = frameshift_oracle(q, index, decoys[case.kind], case.kind, k, case.min_side)
    check_search(dec, res, index)


# --- every operation shape on either side, both δ, both strands --------------------------------

ALL_SEGS = [("substitution", ("E", "W")), ("swap", ("KD", "DK"))]
ALL_SEGS += [("block", seg) for seg in BLOCK_SEGS] + [("near_block", seg) for seg in NEAR_SEGS]


@pytest.mark.parametrize(("i", "op", "seg"), [(i, op, seg) for i, (op, seg) in enumerate(ALL_SEGS)])
def test_every_operation_order_next_to_a_split(i: int, op: str, seg: tuple[str, str]) -> None:
    """Each planted operation (every order of every BLOCKS / NEAR_BLOCKS pair) beside a split."""
    t_seg, q_seg = seg
    gp, gs = "MKHPEPTF", "YCRDLMQH"
    if i % 2:  # on the suffix, right after the junction
        gs_g, gs_q, gp_g, gp_q = t_seg + gs, q_seg + gs, gp, gp
    else:  # on the prefix, right before the junction
        gp_g, gp_q, gs_g, gs_q = gp + t_seg, gp + q_seg, gs, gs
    rng = np.random.default_rng(i)
    nt = join_split(codons_for(gp_g, rng), codons_for(gs_g, rng), (1, -1)[i // 2 % 2])
    if nt is None:
        nt = join_split(codons_for(gp_g), codons_for(gs_g), 1)
    assert nt is not None
    nt = nt if i // 4 % 2 else rc(nt)
    genome = background(250, 500 + i) + nt + background(250, 600 + i)
    index = _index({"g": genome}, "linear")
    q = gp_q + gs_q
    decoys = {"reverse": build_decoys(index, "reverse", 50, for_step("sixframe_index", i))}
    df = match_frameshift([q], index, decoys, {q: "reverse"}, {q}, _cfg(4))
    res = frameshift_oracle(q, index, index.frames, None, _k(len(q)), 4)
    # the planted operation explains some minimal split (the junction may also move instead)
    assert res.cost == 2 and any(op in p.types for p in res.placements)
    check_search(_search(df, q, "target"), res, index)
    res = frameshift_oracle(q, index, decoys["reverse"], "reverse", _k(len(q)), 4)
    check_search(_search(df, q, "decoy"), res, index)


# --- strand symmetry (example 5 generalised) ---------------------------------------------------


def _mirror(row: tuple[object, ...], length: int) -> tuple[object, ...]:
    contig, start, end, strand, frame, n_diff, types, blocks = row
    s, e = int(start), int(end)
    ms = (length - e) % length
    shift = ms - (length - e)
    (a1, b1), (a2, b2) = (tuple(map(int, b.split("-"))) for b in str(blocks).split(";"))
    mb = sorted(
        [(length - b2 + shift, length - a2 + shift), (length - b1 + shift, length - a1 + shift)]
    )
    flipped = "-" if strand == "+" else "+"
    return (
        contig,
        ms,
        ms + e - s,
        flipped,
        -int(frame),
        n_diff,
        types,
        ";".join(f"{a}-{b}" for a, b in mb),
    )  # type: ignore[arg-type]


def _covers(a: list[tuple[object, ...]], b: list[tuple[object, ...]], length: int) -> bool:
    return all(
        any(
            x[0] == y[0]
            and x[3] == y[3]
            and overlaps_mod(
                int(x[1]), int(x[2]) - int(x[1]), int(y[1]), int(y[2]) - int(y[1]), length
            )
            for y in b
        )
        for x in a
    )


@SETTINGS
@given(cases())
def test_reverse_complement_mirrors_the_rows(case: Case) -> None:
    q = case.query
    out = []
    for seq in (case.genome, rc(case.genome)):
        index = _index({"g": seq}, case.topology)
        decoys = {"reverse": build_decoys(index, "reverse", 50, for_step("sixframe_index", 0))}
        df = match_frameshift([q], index, decoys, {q: "reverse"}, {q}, _cfg(case.min_side))
        out.append(_rows(_search(df, q, "target")))
    fwd, bwd = out
    length = len(case.genome)
    mirrored = sorted(_mirror(r, length) for r in fwd)
    assert {r[5] for r in mirrored} == {r[5] for r in bwd}
    assert _covers(mirrored, bwd, length) and _covers(bwd, mirrored, length)
    if len(fwd) == 1 and len(bwd) == 1:
        assert mirrored[0][:7] == bwd[0][:7]


@SETTINGS
@given(cases())
def test_reverse_decoys_of_reversed_query_are_the_target_rows(case: Case) -> None:
    q = case.query
    rq = q[::-1]
    assume(rq != q)
    index = _index({"g": case.genome}, case.topology)
    decoys = {"reverse": build_decoys(index, "reverse", 50, for_step("sixframe_index", 0))}
    df = match_frameshift(
        [q, rq], index, decoys, {q: "reverse", rq: "reverse"}, {q}, _cfg(case.min_side)
    )
    tgt = [r[:4] + r[5:] for r in _rows(_search(df, q, "target"))]
    dec = [r[:4] + r[5:] for r in _rows(_search(df, rq, "decoy"))]
    length = len(case.genome)
    assert {r[4] for r in tgt} == {r[4] for r in dec}
    assert _covers(tgt, dec, length) and _covers(dec, tgt, length)
    if len(tgt) == 1 and len(dec) == 1:
        assert tgt[0][:6] == dec[0][:6]


# --- 50-kb genomes (spec size), several planted splits ------------------------------------------


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_fifty_kb_genomes_equal_brute_force(seed: int) -> None:
    rng = np.random.default_rng(1000 + seed)
    size = 50_000
    genome = background(size, 2000 + seed)
    queries: list[str] = []
    for j, pos in enumerate(range(3_000, 48_000, 9_000)):
        while True:
            gp = "".join(rng.choice(list(SENSE_RESIDUES), size=int(rng.integers(4, 20))))
            gs = "".join(rng.choice(list(SENSE_RESIDUES), size=int(rng.integers(4, 20))))
            if 12 <= len(gp) + len(gs) <= 40:
                break
        qs = gs
        if j % 2:  # one planted substitution on the suffix side
            i = int(rng.integers(len(gs)))
            qs = gs[:i] + ("W" if gs[i] != "W" else "Y") + gs[i + 1 :]
        nt = join_split(codons_for(gp, rng), codons_for(gs, rng), (1, -1)[j % 2])
        if nt is None:
            continue
        nt = nt if j % 3 else rc(nt)
        genome = genome[:pos] + nt + genome[pos + len(nt) :]
        queries.append(gp + qs)
    topology = "circular" if seed % 2 else "linear"
    if topology == "circular":
        genome = genome[3_010:] + genome[:3_010]  # the first split crosses the origin
    index = _index({"g": genome}, topology)
    decoys = {"shuffle": build_decoys(index, "shuffle", 50, for_step("sixframe_index", seed))}
    df = match_frameshift(
        queries, index, decoys, dict.fromkeys(queries, "shuffle"), set(queries), _cfg(4)
    )  # type: ignore[arg-type]
    for q in queries:
        k = _k(len(q))
        res = frameshift_oracle(q, index, index.frames, None, k, 4)
        assert res.cost is not None
        check_search(_search(df, q, "target"), res, index)
        res = frameshift_oracle(q, index, decoys["shuffle"], "shuffle", k, 4)
        check_search(_search(df, q, "decoy"), res, index)


# --- thread count ------------------------------------------------------------------------------

THREAD_SCRIPT = textwrap.dedent(
    """
    import sys
    sys.path.insert(0, sys.argv[1])
    import numpy as np
    from fs_genomes import background, codons_for, join_split, rc
    from giap.config.model import TolerantConfig
    from giap.mapping.decoys import build_decoys
    from giap.mapping.frameshift import match_frameshift
    from giap.mapping.sixframe import build_sixframe
    from giap.model import Replicon
    from giap.rng import for_step

    out = []
    for seed in range(2):
        rng = np.random.default_rng(seed)
        genome = background(20_000, seed)
        queries = []
        for j in range(60):
            gp = "".join(rng.choice(list("ACDEFGHKLMNPQRSTVWY"), size=4 + j % 9))
            gs = "".join(rng.choice(list("ACDEFGHKLMNPQRSTVWY"), size=12 - j % 9))
            nt = join_split(codons_for(gp, rng), codons_for(gs, rng), (1, -1)[j % 2])
            if nt is None:
                continue
            nt = nt if j % 3 else rc(nt)
            pos = 100 + 300 * j
            genome = genome[:pos] + nt + genome[pos + len(nt):]
            queries.append(gp + gs if j % 4 else gp + gs[:-1] + "W")
        queries += ["A" * 12, "AG" * 8]
        topo = "circular" if seed else "linear"
        index = build_sixframe({"g": genome}, {"g": Replicon("g", len(genome), topo)})
        decoys = {"reverse": build_decoys(index, "reverse", 50, for_step("sixframe_index", seed))}
        df = match_frameshift(queries, index, decoys, dict.fromkeys(queries, "reverse"),
                              set(queries[::2]), TolerantConfig())
        out.append(df.to_json(orient="split"))
    print("\\n".join(out))
    """
)


def _run_threads(n: int) -> str:
    root = Path(__file__).resolve().parents[2]
    env = dict(os.environ, NUMBA_NUM_THREADS=str(n), OMP_NUM_THREADS=str(n))
    giap_src = str(Path(giap.__file__).resolve().parents[1])  # the giap under test
    env["PYTHONPATH"] = os.pathsep.join(
        p for p in (giap_src, str(root / "tests" / "sim"), env.get("PYTHONPATH", "")) if p
    )
    proc = subprocess.run(
        [sys.executable, "-c", THREAD_SCRIPT, str(root / "tests" / "helpers")],
        capture_output=True,
        text=True,
        env=env,
        timeout=900,
    )
    assert proc.returncode == 0, proc.stderr
    return proc.stdout


def test_thread_count_does_not_change_output() -> None:
    one, four = _run_threads(1), _run_threads(4)
    first = json.loads(one.splitlines()[0])
    assert first["columns"] and first["data"]
    assert one == four
