"""Hypothesis properties of the exact matcher (issue M3-04; GIAP_SPEC.md §7.2.2, §7.2.9, §11.2).

- "A peptide cut from a random genome maps back to its origin": random genomes <= 50 kb (both
  topologies), a window of 7-40 aa from any strand and phase (origin-crossing on circular
  contigs), some L written as I -> the origin interval is among the target rows.
- The result equals a brute-force ``str.find`` loop over every target and decoy frame.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from giap.mapping.decoys import DecoyType, build_decoys, mirror_aa
from giap.mapping.exact import match_exact
from giap.mapping.hits import HIT_COLUMNS
from giap.mapping.sixframe import SixFrameIndex, build_sixframe, decode, dedupe_loci, encode_query
from giap.model import Interval, Replicon
from giap.rng import for_step
from giapsim.codons import CODE_11

_RC = str.maketrans("ACGT", "TGCA")


def _rc(s: str) -> str:
    return s.translate(_RC)[::-1]


def _genome(seed: int, n: int) -> str:
    codes = np.random.default_rng(seed).integers(0, 4, size=n, dtype=np.uint8)
    return np.frombuffer(b"ACGT", dtype=np.uint8)[codes].tobytes().decode("ascii")


def _translate(nt: str) -> str:
    return "".join(CODE_11[nt[i : i + 3]] for i in range(0, len(nt) - 2, 3))


@settings(max_examples=60, deadline=None)
@given(
    seed=st.integers(0, 2**32 - 1),
    n=st.integers(150, 50_000),
    circular=st.booleans(),
    strand=st.sampled_from("+-"),
    k=st.integers(7, 40),
    pos=st.floats(0, 1, exclude_max=True),
    near_origin=st.booleans(),
    i_mask=st.integers(0, 2**40 - 1),
)
def test_peptide_cut_from_random_genome_maps_back_to_origin(
    seed: int,
    n: int,
    circular: bool,
    strand: str,
    k: int,
    pos: float,
    near_origin: bool,
    i_mask: int,
) -> None:
    seq = _genome(seed, n)
    topo = "circular" if circular else "linear"
    rep = Replicon("c", n, topo)
    span = 3 * k
    if circular:
        start = n - span // 2 - (seed % 3) if near_origin else int(pos * n)
    else:
        start = int(pos * (n - span + 1))
    ext = seq + seq
    nt = ext[start : start + span]
    if strand == "-":
        nt = _rc(nt)
    aa = _translate(nt)
    assume("*" not in aa)
    pep = "".join("I" if c == "L" and (i_mask >> i) & 1 else c for i, c in enumerate(aa))
    idx = build_sixframe({"c": seq}, {"c": rep})
    decoys = {"reverse": build_decoys(idx, "reverse", 50, for_step("sixframe_index", 0))}
    df = match_exact([pep], idx, decoys, {pep: "reverse"})
    tgt = df[df["search"] == "target"]
    got = {Interval(r.contig, int(r.start), int(r.end), r.strand) for r in tgt.itertuples()}
    want = Interval("c", start, start + span, strand)  # type: ignore[arg-type]
    assert want in got
    assert (tgt["peptide"] == pep).all()
    assert (tgt["query"] == pep.replace("I", "L")).all()


def _occurrences(text: str, q: str) -> list[int]:
    out, i = [], text.find(q)
    while i >= 0:
        out.append(i)
        i = text.find(q, i + 1)
    return out


def _brute(
    queries: list[str],
    idx: SixFrameIndex,
    decoys: dict[DecoyType, tuple],  # type: ignore[type-arg]
    decoy_of: dict[str, DecoyType],
    min_len: int = 7,
    max_len: int = 40,
) -> set[tuple[object, ...]]:
    rows: set[tuple[object, ...]] = set()
    targets = [decode(f.aa) for f in idx.frames]
    for pep in set(queries):
        if not min_len <= len(pep) <= max_len:
            continue
        q = decode(encode_query(pep))
        if "#" in q:
            continue
        hits = [
            idx.to_interval(fid, a, a + len(q))
            for fid, text in enumerate(targets)
            for a in _occurrences(text, q)
        ]
        for iv in dedupe_loci(hits, idx.replicons):
            rows.add((pep, q, "target", None, iv.contig, iv.start, iv.end, iv.strand))
        kind = decoy_of[pep]
        dhits = []
        for fid, frame in enumerate(decoys[kind]):
            text = decode(frame.aa)
            for a in _occurrences(text, q):
                ma, mb = mirror_aa(kind, len(frame.aa), a, a + len(q))
                dhits.append(idx.to_interval(fid, ma, mb))
        for iv in dedupe_loci(dhits, idx.replicons):
            rows.add((pep, q, "decoy", kind, iv.contig, iv.start, iv.end, iv.strand))
    return rows


def _row_set(df: pd.DataFrame) -> set[tuple[object, ...]]:
    return {
        (
            r.peptide,
            r.query,
            r.search,
            None if pd.isna(r.decoy_type) else r.decoy_type,
            r.contig,
            int(r.start),
            int(r.end),
            r.strand,
        )
        for r in df.itertuples(index=False)
    }


@settings(max_examples=40, deadline=None)
@given(
    seed=st.integers(0, 2**32 - 1),
    lengths=st.lists(st.integers(150, 3000), min_size=1, max_size=3),
    circular=st.lists(st.booleans(), min_size=3, max_size=3),
    data=st.data(),
)
def test_matches_equal_brute_force_over_every_frame(
    seed: int, lengths: list[int], circular: list[bool], data: st.DataObject
) -> None:
    genome = {f"c{i}": _genome(seed + i, n) for i, n in enumerate(lengths)}
    reps = {
        c: Replicon(c, len(s), "circular" if circular[i] else "linear")
        for i, (c, s) in enumerate(genome.items())
    }
    # Short random repeats raise the chance of multi-locus and overlapping hits.
    idx = build_sixframe(genome, reps)
    rng = for_step("sixframe_index", seed % 1000)
    decoys = {
        "reverse": build_decoys(idx, "reverse", 50, rng),
        "shuffle": build_decoys(idx, "shuffle", 50, rng),
    }
    frames = [decode(f.aa) for f in idx.frames] + [decode(f.aa) for f in decoys["reverse"]]
    n_q = data.draw(st.integers(1, 25))
    queries: list[str] = []
    for _ in range(n_q):
        src = data.draw(st.sampled_from(frames))
        k = data.draw(st.integers(5, 42))
        a = data.draw(st.integers(0, max(0, len(src) - k)))
        pep = src[a : a + k].replace("*", "W")
        if data.draw(st.booleans()):
            pep = pep.replace("L", "I", 1)
        if data.draw(st.integers(0, 9)) == 0:
            pep = pep[:2] + "U" + pep[3:]
        if pep:
            queries.append(pep)
    decoy_of: dict[str, DecoyType] = {
        q: data.draw(st.sampled_from(["reverse", "shuffle"])) for q in queries
    }
    df = match_exact(queries, idx, decoys, decoy_of)
    assert list(df.columns) == list(HIT_COLUMNS)
    got = _row_set(df)
    assert len(got) == len(df)  # one row per peptide x locus x search
    assert got == _brute(queries, idx, decoys, decoy_of)
    for r in df.itertuples(index=False):
        iv = Interval(r.contig, int(r.start), int(r.end), r.strand)
        assert int(r.frame) == idx.frame_label(iv)
        assert r.route == "exact" and int(r.n_diff) == 0 and r.diff_types == ""
        search_frames = idx.frames if r.search == "target" else decoys[r.decoy_type]
        window = decode(search_frames[int(r.frame_id)].aa)[int(r.aa_start) : int(r.aa_end)]
        assert window == r.query
