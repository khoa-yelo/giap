"""Completeness of the q-gram seeds and the tolerant search driver (issue M3-07; §7.2.4, §11.2).

"Property tests MUST show zero false negatives against a brute-force matcher on random genomes
<= 50 kb, for every difference type." Random genomes (linear or circular, rotated so that planted
windows may cross the origin, with occasional ``N``) carry a designed target window; the query is
that window with 0..k planted operations (substitution, swap, every residue order of every
``BLOCKS`` / ``NEAR_BLOCKS`` pair in both directions). The brute force (``brute_tolerant.py``)
scans every window of every frame:

- target and decoy search: the reported minimal cost equals the brute-force minimum (nothing
  above k is reported); every cluster of overlapping minimal-cost windows is covered by a
  reported locus (zero false negatives); every reported window really has that cost and maps to
  the reported (mirrored) interval (zero false positives); no interval is reported twice;
- seeds: for a planted query, :func:`candidate_diagonals` returns a centre ``c`` whose ``±2``
  band holds the whole planted alignment path (the verifier can then find it);
- the target search on the decoy frames equals the decoy search (decoys searched identically);
- thread count 1 vs 4 gives identical frames (subprocesses with ``NUMBA_NUM_THREADS``).
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
from brute_tolerant import (
    INF,
    brute_force,
    overlaps_mod,
    semi_global_costs,
    window_interval,
)
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st
from oracle_tolerant import align_window

import giap
from giap.config.model import TolerantConfig
from giap.mapping.decoys import DecoyType, build_decoys
from giap.mapping.qgram import build_qgram_index, candidate_diagonals, choose_q
from giap.mapping.sixframe import (
    Frame,
    SixFrameIndex,
    build_sixframe,
    decode,
    encode_query,
    encode_target,
)
from giap.mapping.tolerant import match_tolerant
from giap.model import Replicon
from giap.rng import for_step
from giapsim.codons import CODE_11

CFG = TolerantConfig()
RESIDUES = "ACDEFGHKLMNPQRSTVWY"
SETTINGS = settings(max_examples=80, deadline=None, suppress_health_check=[HealthCheck.too_slow])
CODONS: dict[str, list[str]] = {}
for _c, _a in sorted(CODE_11.items()):
    CODONS.setdefault(_a, []).append(_c)
_RC = str.maketrans("ACGTN", "TGCAN")
DIFF_VOCAB = {"substitution", "swap", "block", "near_block"}

_BLOCK_RAW = (("N", "GG"), ("Q", "AG"), ("GT", "AS"), ("AV", "GL"), ("AN", "GQ"))
_NEAR_RAW = (("W", "AD"), ("W", "EG"), ("R", "GV"))


def _orders(s: str) -> list[str]:
    return sorted({"".join(p) for p in permutations(s)})


def _expand(raw: tuple[tuple[str, str], ...]) -> list[tuple[str, str]]:
    out = []
    for a, b in raw:
        for x in _orders(a):
            for y in _orders(b):
                out += [(x, y), (y, x)]
    return out


BLOCK_SEGS = _expand(_BLOCK_RAW)  # (target segment, query segment)
NEAR_SEGS = _expand(_NEAR_RAW)


def _k(L: int) -> int:
    return 1 if L <= 11 else 2


def _rc(s: str) -> str:
    return s.translate(_RC)[::-1]


@dataclass(frozen=True)
class Plant:
    window: str  # target residues
    query: str
    n_ops: int
    lo: int  # min / max diagonal offset (target - query position) along the planted path
    hi: int


@st.composite
def op_piece(draw: st.DrawFn, kind: str) -> tuple[str, str]:
    """(target segment, query segment) of one planted operation."""
    if kind == "substitution":
        a = draw(st.sampled_from(RESIDUES))
        b = draw(st.sampled_from(RESIDUES.replace(a, "")))
        return a, b
    if kind == "swap":
        a = draw(st.sampled_from(RESIDUES))
        b = draw(st.sampled_from(RESIDUES.replace(a, "")))
        return a + b, b + a
    return draw(st.sampled_from(BLOCK_SEGS if kind == "block" else NEAR_SEGS))


@st.composite
def plants(draw: st.DrawFn) -> Plant:
    n_ops = draw(st.integers(0, 2))
    chunk = st.text(alphabet=RESIDUES, min_size=1, max_size=9)
    first = draw(chunk)
    window, query = [first], [first]
    off = lo = hi = 0
    for _ in range(n_ops):
        kind = draw(st.sampled_from(["substitution", "swap", "block", "near_block"]))
        t_seg, q_seg = draw(op_piece(kind))
        c = draw(chunk)
        window += [t_seg, c]
        query += [q_seg, c]
        off += len(t_seg) - len(q_seg)
        lo, hi = min(lo, off), max(hi, off)
    w, q = "".join(window), "".join(query)
    if len(q) < 10:
        pad = draw(st.text(alphabet=RESIDUES, min_size=10 - len(q), max_size=14 - len(q)))
        w, q = w + pad, q + pad
    assume(len(q) <= 30 and n_ops <= _k(len(q)))
    if draw(st.booleans()) and "L" in q:  # I = L is free
        i = q.index("L")
        q = q[:i] + "I" + q[i + 1 :]
    return Plant(w, q, n_ops, lo, hi)


def _back(aa: str, rng: np.random.Generator) -> str:
    return "".join(CODONS[a][int(rng.integers(len(CODONS[a])))] for a in aa)


@dataclass(frozen=True)
class World:
    index: SixFrameIndex
    decoys: dict[DecoyType, tuple[Frame, ...]]
    kind: DecoyType
    plant: Plant


@st.composite
def worlds(draw: st.DrawFn) -> World:
    plant = draw(plants())
    seed = draw(st.integers(0, 2**32 - 1))
    rng = np.random.default_rng(seed)
    n_bg = draw(st.integers(150, 1500))  # longer than any window
    bg = "".join(rng.choice(list("ACGT"), size=n_bg))
    if n_bg and draw(st.integers(0, 4)) == 0:  # some N -> X codons
        pos = rng.integers(0, n_bg, size=3)
        bg = "".join("N" if i in set(pos.tolist()) else c for i, c in enumerate(bg))
    seg = _back(plant.window, rng)
    if draw(st.booleans()):
        seg = _rc(seg)
    at = draw(st.integers(0, n_bg))
    genome = bg[:at] + seg + bg[at:]
    if draw(st.booleans()):  # a reversed copy: a positive reverse-decoy hit
        rev = _back(plant.window[::-1], rng)
        at2 = draw(st.integers(0, len(genome) - len(seg)))
        at2 += len(seg) if at2 > at else 0  # never inside the planted window
        genome = genome[:at2] + (rev if draw(st.booleans()) else _rc(rev)) + genome[at2:]
    topology = draw(st.sampled_from(["linear", "circular"]))
    if topology == "circular":
        r = draw(st.integers(0, len(genome) - 1))
        genome = genome[r:] + genome[:r]
    rep = {"g": Replicon("g", len(genome), topology)}
    index = build_sixframe({"g": genome}, rep)  # type: ignore[arg-type]
    kind: DecoyType = draw(st.sampled_from(["reverse", "shuffle"]))
    decoys = {kind: build_decoys(index, kind, 50, for_step("sixframe_index", seed))}
    return World(index, decoys, kind, plant)


def _check_search(
    df: pd.DataFrame,
    peptide: str,
    index: SixFrameIndex,
    frames: tuple[Frame, ...],
    kind: DecoyType | None,
    *,
    complete: bool = True,
) -> None:
    k = _k(len(peptide))
    brute = brute_force(peptide, frames, k)
    qletters = decode(encode_query(peptide))
    if brute.cost is None:
        assert len(df) == 0, df
        return
    if complete:
        assert len(df) > 0, ("false negative", peptide, brute)
        assert set(df["n_diff"].tolist()) == {brute.cost}, (df, brute)
    else:
        assert set(df["n_diff"].tolist()) <= set(range(brute.cost, k + 1))
    seen = set()
    for r in df.itertuples(index=False):
        n = int(r.n_diff)
        text = decode(frames[int(r.frame_id)].aa)
        got = align_window(qletters, text, int(r.aa_start), int(r.aa_end), None, n)
        assert got is not None and got[0] == n, ("false positive", r)
        types = [] if r.diff_types == "" else str(r.diff_types).split(";")
        assert len(types) == n and set(types) <= DIFF_VOCAB, r
        contig, start, length, strand, label = window_interval(
            index, frames, kind, int(r.frame_id), int(r.aa_start), int(r.aa_end)
        )
        rep_len = index.replicons[contig].length
        assert (r.contig, r.strand, int(r.frame)) == (contig, strand, label), r
        assert int(r.start) % rep_len == start and int(r.end) - int(r.start) == length, r
        key = (r.contig, int(r.start) % rep_len, int(r.end), r.strand)
        assert key not in seen, ("duplicate locus", r)
        seen.add(key)
    if not complete:
        return
    for cluster in brute.clusters:
        covered = False
        for fid, s, e in cluster:
            try:
                contig, start, length, strand, label = window_interval(
                    index, frames, kind, fid, s, e
                )
            except ValueError:  # longer than the replicon: never reported
                continue
            rep_len = index.replicons[contig].length
            for r in df.itertuples(index=False):
                if (r.contig, r.strand, int(r.frame)) != (contig, strand, label):
                    continue
                r_start = int(r.start) % rep_len
                r_len = int(r.end) - int(r.start)
                if overlaps_mod(start, length, r_start, r_len, rep_len):
                    covered = True
                    break
            if covered:
                break
        assert covered, ("false negative: missed locus", peptide, cluster, df)


def _search(df: pd.DataFrame, search: str) -> pd.DataFrame:
    return df[df["search"] == search].reset_index(drop=True)


# --- the brute force itself --------------------------------------------------------------------


@SETTINGS
@given(
    q=st.text(alphabet=RESIDUES + "B", min_size=5, max_size=8),
    t=st.text(alphabet=RESIDUES + "*X", min_size=0, max_size=12),
    k=st.integers(0, 2),
)
def test_brute_force_dp_equals_window_oracle(q: str, t: str, k: int) -> None:
    qc, tc = encode_query(q), encode_target(t)
    qletters, tletters = decode(qc), decode(tc)
    costs = semi_global_costs(qc, tc, k)
    for e in range(len(t) + 1):
        best = None
        for s in range(e + 1):
            got = align_window(qletters, tletters, s, e, None, k)
            if got is not None and (best is None or got[0] < best):
                best = got[0]
        assert (INF if best is None else best) == int(costs[e]), (q, t, e)


# --- completeness and soundness of match_tolerant ----------------------------------------------


@SETTINGS
@given(w=worlds())
def test_match_tolerant_equals_brute_force(w: World) -> None:
    pep = w.plant.query
    df = match_tolerant([pep], w.index, w.decoys, {pep: w.kind}, {pep}, CFG)
    assert set(df["peptide"]) <= {pep}
    tgt = _search(df, "target")
    dec = _search(df, "decoy")
    assert len(tgt) > 0, "planted window not found"
    assert int(tgt["n_diff"].min()) <= w.plant.n_ops
    assert tgt["decoy_type"].isna().all()
    assert (dec["decoy_type"] == w.kind).all()
    assert set(df["route"]) <= {"tolerant"}
    _check_search(tgt, pep, w.index, w.index.frames, None)
    _check_search(dec, pep, w.index, w.decoys[w.kind], w.kind)


@SETTINGS
@given(w=worlds(), data=st.data())
def test_nomatch_queries_never_report_false_positives(w: World, data: st.DataObject) -> None:
    q = w.plant.query
    pos = data.draw(st.integers(0, len(q) - 1))
    pep = q[:pos] + data.draw(st.sampled_from("BJOUXZ")) + q[pos + 1 :]
    df = match_tolerant([pep], w.index, w.decoys, {pep: w.kind}, {pep}, CFG)
    _check_search(_search(df, "target"), pep, w.index, w.index.frames, None, complete=False)
    dec = _search(df, "decoy")
    _check_search(dec, pep, w.index, w.decoys[w.kind], w.kind, complete=False)


@SETTINGS
@given(w=worlds())
def test_target_search_on_decoy_frames_equals_decoy_search(w: World) -> None:
    linear = {c: Replicon(r.id, r.length, "linear") for c, r in w.index.replicons.items()}
    assume(all(f.ext_nt == linear[f.contig].length for f in w.index.frames))  # no extension
    pep = w.plant.query
    frames = w.decoys[w.kind]
    flipped = SixFrameIndex(replicons=w.index.replicons, frames=frames)
    a = match_tolerant([pep], w.index, w.decoys, {pep: w.kind}, set(), CFG)
    b = match_tolerant([pep], flipped, {w.kind: w.index.frames}, {pep: w.kind}, {pep}, CFG)
    cols = ["frame_id", "aa_start", "aa_end", "n_diff", "diff_types"]
    dec = _search(a, "decoy")[cols].sort_values(cols).reset_index(drop=True)
    tgt = _search(b, "target")[cols].sort_values(cols).reset_index(drop=True)
    pd.testing.assert_frame_equal(dec, tgt)


# --- seeds -------------------------------------------------------------------------------------


@SETTINGS
@given(p=plants(), seed=st.integers(0, 2**32 - 1), data=st.data())
def test_seeds_cover_the_planted_path(p: Plant, seed: int, data: st.DataObject) -> None:
    rng = np.random.default_rng(seed)
    sizes = [int(n) for n in rng.integers(0, 400, size=data.draw(st.integers(1, 6)))]
    frames = []
    for n in sizes:
        aa = rng.integers(0, 19, size=n).astype(np.uint8)
        aa[rng.random(n) < 0.05] = 19
        aa[rng.random(n) < 0.01] = 20
        frames.append(aa)
    fid = data.draw(st.integers(0, len(frames) - 1))
    w = encode_target(p.window)
    host = frames[fid]
    s = data.draw(st.integers(0, len(host)))
    frames[fid] = np.concatenate([host[:s], w, host[s:]]).astype(np.uint8)
    frs = tuple(Frame("c", "+", 0, a, 3 * len(a)) for a in frames)
    L = len(p.query)
    q, t = choose_q(L, _k(L))
    idx = build_qgram_index(frs, q)
    d = int(np.asarray(idx.frame_starts)[fid]) + s
    cands = np.asarray(candidate_diagonals(encode_query(p.query), idx, t))
    ok = (cands >= d + p.hi - 2) & (cands <= d + p.lo + 2)
    assert ok.any(), (cands, d, p)


# --- fixed 50-kb genomes (the spec's upper size), both topologies ------------------------------


@pytest.mark.parametrize("topology", ["linear", "circular"])
def test_fifty_kb_genomes_equal_brute_force(topology: str) -> None:
    rng = np.random.default_rng(50 if topology == "linear" else 51)
    genome = "".join(rng.choice(list("ACGT"), size=50_000))
    rep = {"g": Replicon("g", len(genome), topology)}  # type: ignore[arg-type]
    index = build_sixframe({"g": genome}, rep)
    decoys = {
        "reverse": build_decoys(index, "reverse", 50, for_step("sixframe_index", 1)),
        "shuffle": build_decoys(index, "shuffle", 50, for_step("sixframe_index", 1)),
    }
    queries: list[str] = []
    decoy_of: dict[str, DecoyType] = {}
    ops = ["substitution", "swap"]
    for i in range(24):
        f = index.frames[int(rng.integers(6))]
        text = decode(f.aa)
        L = int(rng.integers(10, 25))
        while True:
            s = int(rng.integers(0, len(text) - L))
            if "*" not in text[s : s + L] and "X" not in text[s : s + L]:
                break
        pep = list(text[s : s + L])
        for _ in range(int(rng.integers(0, _k(L) + 1))):
            j = int(rng.integers(0, L - 1))
            if ops[i % 2] == "swap" and pep[j] != pep[j + 1]:
                pep[j], pep[j + 1] = pep[j + 1], pep[j]
            else:
                pep[j] = RESIDUES[(RESIDUES.index(pep[j]) + 3) % 19]
        p = "".join(pep)
        queries.append(p)
        decoy_of[p] = "reverse" if i % 3 else "shuffle"
    df = match_tolerant(queries, index, decoys, decoy_of, set(queries), CFG)  # type: ignore[arg-type]
    for p in sorted(set(queries)):
        rows = df[df["peptide"] == p]
        tgt = _search(rows, "target")
        assert len(tgt) > 0
        _check_search(tgt, p, index, index.frames, None)
        kind = decoy_of[p]
        _check_search(_search(rows, "decoy"), p, index, decoys[kind], kind)


# --- thread count ------------------------------------------------------------------------------

THREAD_SCRIPT = textwrap.dedent(
    """
    import numpy as np
    from giap.config.model import TolerantConfig
    from giap.mapping.decoys import build_decoys
    from giap.mapping.sixframe import build_sixframe, decode
    from giap.mapping.tolerant import match_tolerant
    from giap.model import Replicon
    from giap.rng import for_step

    out = []
    for seed in range(3):
        rng = np.random.default_rng(seed)
        genome = "".join(rng.choice(list("ACGT"), size=20_000))
        topo = "circular" if seed % 2 else "linear"
        index = build_sixframe({"g": genome}, {"g": Replicon("g", len(genome), topo)})
        decoys = {"shuffle": build_decoys(index, "shuffle", 50, for_step("sixframe_index", seed))}
        queries = []
        for i in range(200):
            text = decode(index.frames[i % 6].aa)
            s = int(rng.integers(0, len(text) - 30))
            q = text[s : s + 10 + i % 20].replace("*", "A").replace("X", "G")
            if i % 3:
                q = q[:4] + "W" + q[5:]
            queries.append(q)
        queries += ["A" * 12, "AG" * 8]
        df = match_tolerant(queries, index, decoys, dict.fromkeys(queries, "shuffle"),
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
        [sys.executable, "-c", THREAD_SCRIPT], capture_output=True, text=True, env=env, timeout=900
    )
    assert proc.returncode == 0, proc.stderr
    return proc.stdout


def test_thread_count_does_not_change_output() -> None:
    one, four = _run_threads(1), _run_threads(4)
    assert one.strip() and json.loads(one.splitlines()[0])["columns"]
    assert one == four
