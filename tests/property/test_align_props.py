"""Properties of the M7-01 residue aligner (docs/specs/M7-01-residue-alignment.md, acceptance).

Random circular or linear genomes (<= 5 kb, random background) carry a locus written on either
strand, possibly across the origin of a circular contig. The genome residues ``G`` and the
peptide ``P`` are built together from segments: matched residues (an I/L flip now and then) and
at most two planted differences of every §7.2.4 type (substitution, adjacent swap, every order
of every ``BLOCKS`` and ``NEAR_BLOCKS`` pair, both directions). ``n_diff`` is the oracle cost of
the locus (``helpers.align_oracle``: exhaustive banded enumeration, independent of ``giap``).

- the aligner's ops equal the oracle's rule-4 optimum (cost, diff-type multiset, leftmost);
- its cost (non-``match``/``il`` ops) equals ``n_diff`` and the M3-06 verifier's cost whenever
  the verifier's best window is the locus itself;
- replaying the ops gives back the peptide (I/L as given) and the genome translation, and the
  codons are the locus codons in peptide order (rule 8);
- ``n_diff`` off by one raises ``ValueError``;
- mirroring the genome (reverse complement) and the hit keeps every op's peptide positions,
  kinds and letters and mirrors its codon coordinates;
- frameshift splits (±1 nt, both strands, across the origin): the junction event, both sides'
  ops and the shared budget.
"""

from __future__ import annotations

import importlib
import sys
from itertools import permutations
from pathlib import Path
from types import ModuleType

import numpy as np
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))

from align_oracle import (
    _key,
    all_alignments,
    best_path,
    codon_at,
    cost_of,
    diff_types,
    expected_ops,
    fold,
    locus_aa,
    op_tuples,
    place,
    rc,
    to_ops,
)
from fs_genomes import background, codons_for, join_split

from giap.mapping.diffops import block_tables
from giap.mapping.sixframe import encode_query, encode_target
from giap.mapping.tolerant_dp import verify
from giap.model import Replicon

SETTINGS = settings(max_examples=60, deadline=None, suppress_health_check=list(HealthCheck))
SENSE = "ACDEFGHIKLMNPQRSTVWY"
_PAIRS = (("N", "GG"), ("Q", "AG"), ("GT", "AS"), ("AV", "GL"), ("AN", "GQ"))
_NEAR = (("W", "AD"), ("W", "EG"), ("R", "GV"))


def _orders(s: str) -> list[str]:
    return sorted({"".join(p) for p in permutations(s)})


#: (genome side, peptide side) for every order of every exchange, both directions
EXCHANGES: dict[str, list[tuple[str, str]]] = {
    kind: [
        (g, p)
        for a, b in pairs
        for x, y in ((a, b), (b, a))
        for g in _orders(x)
        for p in _orders(y)
    ]
    for kind, pairs in (("block", _PAIRS), ("near_block", _NEAR))
}


def align() -> ModuleType:
    return importlib.import_module("giap.variants.align")


@st.composite
def segments(draw: st.DrawFn, lo: int, hi: int, max_diffs: int) -> tuple[str, str]:
    """``(G, P)``: genome residues and peptide with <= ``max_diffs`` planted differences."""
    n = draw(st.integers(lo, hi))
    segs = [(a, a) for a in draw(st.lists(st.sampled_from(SENSE), min_size=n, max_size=n))]
    for i in draw(st.lists(st.integers(0, n - 1), max_size=3)):
        g = segs[i][0]
        if g in "IL":
            segs[i] = (g, "L" if g == "I" else "I")
    for _ in range(draw(st.integers(0, max_diffs))):
        i = draw(st.integers(0, len(segs) - 1))
        kind = draw(st.sampled_from(("substitution", "swap", "block", "near_block")))
        if kind == "substitution":
            g = draw(st.sampled_from(SENSE))
            p = draw(st.sampled_from([c for c in SENSE if fold(c) != fold(g)]))
        elif kind == "swap":
            x = draw(st.sampled_from(SENSE))
            y = draw(st.sampled_from([c for c in SENSE if fold(c) != fold(x)]))
            g, p = x + y, y + x
        else:
            g, p = draw(st.sampled_from(EXCHANGES[kind]))
        segs[i] = (g, p)
    return "".join(g for g, _ in segs), "".join(p for _, p in segs)


@st.composite
def loci(draw: st.DrawFn) -> dict:
    gres, pep = draw(segments(10, 30, 2))
    span = 3 * len(gres)
    circular = draw(st.booleans())
    length = draw(st.integers(max(300, span + 30), 5000))
    strand = draw(st.sampled_from("+-"))
    if circular:
        start = draw(
            st.one_of(st.integers(length - span + 1, length - 1), st.integers(0, length - 1))
        )
    else:
        start = draw(st.integers(0, length - span))
    seed = draw(st.integers(0, 2**31 - 1))
    body = "".join(codons_for(gres, np.random.default_rng(seed)))
    seq = place(background(length, seed), start, body, strand)
    rep = Replicon("c1", length, "circular" if circular else "linear")
    return {"g": gres, "p": pep, "seq": seq, "rep": rep, "start": start,
            "end": start + span, "strand": strand}  # fmt: skip


def _hit(mod: ModuleType, c: dict, n_diff: int, seq: str | None = None):
    return mod.align_hit(c["p"], "c1", c["start"], c["end"], c["strand"], n_diff, None,
                         c["seq"] if seq is None else seq, c["rep"])  # fmt: skip


@SETTINGS
@given(loci())
def test_ops_equal_the_oracle_and_replay(c: dict) -> None:
    mod = align()
    assert locus_aa(c["seq"], c["start"], c["end"], c["strand"]) == c["g"]
    path = best_path(c["p"], c["g"], 2)
    assume(path is not None)
    n_diff = _key(path)[0]
    a = _hit(mod, c, n_diff)
    ops = op_tuples(a)
    assert ops == expected_ops(c["p"], c["seq"], c["start"], c["end"], c["strand"])
    assert cost_of(ops) == n_diff
    assert "".join(o[4] for o in ops) == c["p"]
    assert "".join(o[3] for o in ops) == c["g"]
    want_codons = [codon_at(j, c["start"], c["end"], c["strand"]) for j in range(len(c["g"]))]
    assert [x for o in ops for x in o[2]] == want_codons
    assert [o[1] for o in ops] == sorted(o[1] for o in ops)
    assert all(o[1][0] < o[1][1] for o in ops)
    assert all(o[1][1] == ops[k + 1][1][0] for k, o in enumerate(ops[:-1]))
    assert a.frameshift is None and a.peptide == c["p"] and a.strand == c["strand"]
    for kind, _, _, ref, alt in ops:
        if kind == "match":
            assert ref == alt
        elif kind == "il":
            assert ref != alt and fold(ref) == fold(alt)
    # M3-06: the verifier's optimum on the locus translation, when its window is the locus
    got = verify(encode_query(c["p"]), encode_target(c["g"]), 0, 2, block_tables())
    assert got, "M3-06 finds nothing where the oracle aligns"
    assert got[0].cost <= n_diff
    if (got[0].t_start, got[0].t_end) == (0, len(c["g"])):
        assert got[0].cost == n_diff
        assert got[0].diff_types == diff_types(path)


@SETTINGS
@given(loci())
def test_n_diff_off_by_one_raises(c: dict) -> None:
    mod = align()
    path = best_path(c["p"], c["g"], 2)
    assume(path is not None)
    n_diff = _key(path)[0]
    try:
        _hit(mod, c, n_diff + 1)
    except ValueError as exc:
        assert c["p"] in str(exc)
    else:
        raise AssertionError(f"n_diff {n_diff + 1} accepted for a cost-{n_diff} locus")
    if n_diff >= 1:
        try:
            _hit(mod, c, n_diff - 1)
        except ValueError as exc:
            assert c["p"] in str(exc)
        else:
            raise AssertionError(f"n_diff {n_diff - 1} accepted for a cost-{n_diff} locus")


@SETTINGS
@given(loci())
def test_mirroring_maps_ops_to_the_same_peptide_positions(c: dict) -> None:
    mod = align()
    path = best_path(c["p"], c["g"], 2)
    assume(path is not None)
    n_diff = _key(path)[0]
    length = c["rep"].length
    a = _hit(mod, c, n_diff)
    ms, me = length - c["end"], length - c["start"]
    shift = length if ms < 0 else 0
    mirrored = dict(c, seq=rc(c["seq"]), start=ms + shift, end=me + shift,
                    strand="-" if c["strand"] == "+" else "+")  # fmt: skip
    b = _hit(mod, mirrored, n_diff)
    oa, ob = op_tuples(a), op_tuples(b)
    assert [(k, p, r, x) for k, p, _, r, x in oa] == [(k, p, r, x) for k, p, _, r, x in ob]
    assert [tuple(length - x - 3 + shift for x in o[2]) for o in oa] == [o[2] for o in ob]


@st.composite
def splits(draw: st.DrawFn) -> dict:
    (gp, pp), (gs, ps) = draw(segments(4, 15, 0)), draw(segments(4, 15, 0))
    if draw(st.booleans()):  # one further difference, on either side
        side = draw(st.sampled_from("ps"))
        g, p = draw(segments(4, 15, 1))
        gp, pp, gs, ps = (g, p, gs, ps) if side == "p" else (gp, pp, g, p)
    delta = draw(st.sampled_from((1, -1)))
    seed = draw(st.integers(0, 2**31 - 1))
    rng = np.random.default_rng(seed)
    body = join_split(codons_for(gp, rng), codons_for(gs, rng), delta)
    assume(body is not None)
    assert body is not None
    circular = draw(st.booleans())
    length = draw(st.integers(max(300, len(body) + 30), 5000))
    if circular:
        start = draw(st.one_of(st.integers(length - len(body) + 1, length - 1),
                               st.integers(0, length - 1)))  # fmt: skip
    else:
        start = draw(st.integers(0, length - len(body)))
    strand = draw(st.sampled_from("+-"))
    seq = place(background(length, seed), start, body, strand)
    end = start + len(body)
    if strand == "+":
        pre = (start, start + 3 * len(gp))
        suf = (pre[1] + delta, end)
    else:
        pre = (end - 3 * len(gp), end)
        suf = (start, start + 3 * len(gs))
    lo, hi = sorted((pre, suf))
    return {"gp": gp, "pp": pp, "gs": gs, "ps": ps, "delta": delta, "seq": seq, "start": start,
            "end": end, "strand": strand, "pre": pre, "suf": suf,
            "blocks": f"{lo[0]}-{lo[1]};{hi[0]}-{hi[1]}",
            "rep": Replicon("c1", length, "circular" if circular else "linear")}  # fmt: skip


def _codons(block: tuple[int, int], n: int, strand: str) -> list[int]:
    return [codon_at(j, block[0], block[1], strand) for j in range(n)]


@SETTINGS
@given(splits())
def test_frameshift_split_event_ops_and_budget(c: dict) -> None:
    mod = align()
    pep, gp, gs, strand = c["pp"] + c["ps"], c["gp"], c["gs"], c["strand"]
    assert locus_aa(c["seq"], *c["pre"], strand) == gp
    assert locus_aa(c["seq"], *c["suf"], strand) == gs
    # oracle: every peptide split point p, each side aligned on its own block (rule 7)
    cands = []
    for p in range(1, len(pep)):
        a, b = best_path(pep[:p], gp, 1), best_path(pep[p:], gs, 1)
        if a is None or b is None or _key(a)[0] + _key(b)[0] > 1:
            continue
        ops = to_ops(a, pep[:p], gp, _codons(c["pre"], len(gp), strand))
        ops += to_ops(b, pep[p:], gs, _codons(c["suf"], len(gs), strand), pep_offset=p)
        k = (cost_of(ops), tuple(sorted(o[0] for o in ops if o[0] not in ("match", "il"))),
             tuple(o[1][0] for o in ops if o[0] not in ("match", "il")))  # fmt: skip
        cands.append((k, p, ops))
    assume(cands)
    best = min(k for k, _, _ in cands)
    n_diff = 1 + best[0]
    got = mod.align_hit(pep, "c1", c["start"], c["end"], strand, n_diff, c["blocks"], c["seq"],
                        c["rep"])  # fmt: skip
    ops = op_tuples(got)
    assert 1 + cost_of(ops) == n_diff
    assert "".join(o[4] for o in ops) == pep
    assert "".join(o[3] for o in ops) == gp + gs
    want_codons = _codons(c["pre"], len(gp), strand) + _codons(c["suf"], len(gs), strand)
    assert [x for o in ops for x in o[2]] == want_codons
    fs = got.frameshift
    assert fs is not None
    assert fs.shift == c["delta"]
    assert fs.last_prefix_codon == want_codons[len(gp) - 1]
    assert fs.first_suffix_codon == want_codons[len(gp)]
    winners = [(p, o) for k, p, o in cands if k == best]
    if len({p for p, _ in winners}) == 1:
        assert fs.pep_junction == winners[0][0]
        assert ops == winners[0][1]
    else:
        assert fs.pep_junction in {p for p, _ in winners}
    try:
        mod.align_hit(pep, "c1", c["start"], c["end"], strand, n_diff + 1, c["blocks"],
                      c["seq"], c["rep"])  # fmt: skip
    except ValueError as exc:
        assert pep in str(exc)
    else:
        raise AssertionError("frameshift row with n_diff one too high accepted")


def test_oracle_enumeration_is_global_and_banded() -> None:
    """Guards the oracle: a 1->2 block shifts the diagonal by one and costs one."""
    paths = all_alignments("AGGGEFKLLR", "ANGEFKLLR", 1)
    assert paths and all(len([s for s in p if s[0] != "match"]) == 1 for p in paths)
    assert best_path("AAAAAAAAAA", "KKKKKKKKKK", 2) is None
