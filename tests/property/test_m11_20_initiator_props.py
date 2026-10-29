"""Properties P1-P5 of the initiator-M credit (docs/specs/M11-20-initiator-methionine.md § Tests,
rules 1-3, S1). Random genomes (300-3,000 nt, linear or circular, tables 4/11/25) carry 1-4
GTG/TTG (and ATG) genes on either strand, sometimes across the origin; ``aa`` is ``M`` plus the
table's translation. Cells are derived independently (``helpers.initiator_toys.ref_cells``:
every frame position mapping onto a start key) and the brute force of rule 3 is
``helpers.initiator_toys.credited_rows``.
"""

from __future__ import annotations

import importlib
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
from helpers.initiator_toys import (
    STARTS,
    cds,
    credited_rows,
    mirrored,
    ref_cells,
    row_set,
    start_keys,
    usable,
    windows,
    with_decoys,
)
from helpers.readthrough_toys import model, rc
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.mapping.decoys import mirror_aa
from gippa.mapping.exact import match_exact
from gippa.mapping.sixframe import decode, translate_text
from gippa.model import Replicon
from gippa.model.genetic_code import stop_codons

SETTINGS = settings(max_examples=40, deadline=None, suppress_health_check=list(HealthCheck))
KS = range(7, 41)


def ini() -> ModuleType:
    return importlib.import_module("gippa.mapping.initiator")


def _sense(table: int) -> list[str]:
    stops = stop_codons(table)
    return [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in stops]


@st.composite
def worlds(draw: Any) -> dict[str, Any]:
    """A random contig with 1-4 non-overlapping genes (start, body of sense codons, stop)."""
    seed = draw(st.integers(0, 2**31 - 1))
    table = draw(st.sampled_from((4, 11, 25)))
    circular = draw(st.booleans())
    n = draw(st.integers(300, 3000))
    rng = np.random.default_rng(seed)
    seq = list("".join(rng.choice(list("ACGT"), size=n)))
    sense, stop = _sense(table), sorted(stop_codons(table))[0]
    genes, used = [], 0
    pos = int(rng.integers(0, n)) if circular else int(rng.integers(0, 30))
    for i in range(draw(st.integers(1, 4))):
        n_aa = int(rng.integers(8, 50))
        span = 3 * (n_aa + 1)
        if used + span + 10 > n or (not circular and pos + span > n):
            break
        first = draw(st.sampled_from(("GTG", "TTG", "GTG", "ATG")))
        body = "".join(rng.choice(sense, size=n_aa - 1))
        gs = first + body + stop  # gene strand
        strand = draw(st.sampled_from("+-"))
        for j, b in enumerate(gs if strand == "+" else rc(gs)):
            seq[(pos + j) % n] = b
        genes.append(
            cds(f"G{i}", "c", pos, pos + span, strand, "M" + translate_text(body, table), first)
        )
        used += span + 10
        pos = (pos + span + 10) % n if circular else pos + span + 10
    reps = {"c": Replicon("c", n, "circular" if circular else "linear")}
    return {"genome": {"c": "".join(seq)}, "reps": reps, "gm": model(reps, genes), "table": table}


def _built(w: dict[str, Any]) -> tuple[Any, Any, Any, dict[str, Any]]:
    m = ini()
    index, decoys = with_decoys(w["genome"], w["reps"], w["table"])
    sites, _ = m.annotated_starts(w["gm"], w["genome"], STARTS, table=w["table"])
    target = m.start_frame_cells(index, sites)
    cells = m.StartCells(
        sites, target, {k: m.mirror_start_cells(k, fr, target) for k, fr in decoys.items()}
    )
    ref_t = ref_cells(index, start_keys(w["gm"], w["genome"]))
    ref = {
        "target": ref_t,
        "reverse": mirrored("reverse", decoys["reverse"], ref_t),
        "shuffle": ref_t,
    }
    return index, decoys, cells, ref


def _prefix_ivs(g: Any, k: int) -> tuple[int, int]:
    iv = g.interval
    return (iv.start, iv.start + 3 * k) if iv.strand == "+" else (iv.end - 3 * k, iv.end)


@SETTINGS
@given(worlds())
def test_p1_completeness(w: dict[str, Any]) -> None:
    index, decoys, cells, _ = _built(w)
    reps = w["reps"]
    rows, want = [], []
    for g in w["gm"].genes:
        if g.start_codon == "ATG":
            continue
        for k in range(7, min(40, len(g.aa)) + 1):
            if "*" in g.aa[:k] or "X" in g.aa[:k]:
                continue
            rows.append(g.aa[:k])
            want.append((g.aa[:k], g.interval.strand, *_prefix_ivs(g, k)))
    if not rows:
        return
    df = match_exact(rows, index, decoys, dict.fromkeys(rows, "reverse"), starts=cells)
    tgt = df[df["search"] == "target"].reset_index(drop=True)
    flags = ini().initiator_rows(
        tgt, cells.sites, {"c": reps["c"].length}, {"c": reps["c"].topology == "circular"}
    )
    got = {
        (p, s, int(a), int(b))
        for p, s, a, b, f in zip(
            tgt["peptide"], tgt["strand"], tgt["start"], tgt["end"], flags, strict=True
        )
        if f
    }
    n = reps["c"].length
    norm = {(p, s, a % n, a % n + (b - a)) for p, s, a, b in want}
    assert norm <= got, sorted(norm - got)[:5]


def _queries(
    w: dict[str, Any], index: Any, decoys: Any, ref: dict[str, Any], seed: int
) -> list[str]:
    """True prefixes, decoy-cell windows (both rules), V/L-initial and random M-initial strings."""
    rng = np.random.default_rng(seed)
    out = {g.aa[:k] for g in w["gm"].genes for k in (7, 9, 12) if len(g.aa) >= k}
    for kind, frames in (
        ("target", index.frames),
        ("reverse", decoys["reverse"]),
        ("shuffle", decoys["shuffle"]),
    ):
        for fid, f in enumerate(frames):
            for k in (7, 8, 11):
                for _, a, b in windows(kind, len(f.aa), ref[kind][fid], k):
                    rest = decode(f.aa[a : b - 1] if kind == "reverse" else f.aa[a + 1 : b])
                    out |= {"M" + rest, "V" + rest, "M" + rest[::-1]}
    out |= {"M" + "".join(rng.choice(list("ACDEFGHKLMNPQRSTVWY"), size=7)) for _ in range(20)}
    return sorted(q for q in out if "*" not in q and "X" not in q)


@SETTINGS
@given(worlds(), st.integers(0, 1000))
def test_p2_brute_force(w: dict[str, Any], seed: int) -> None:
    index, decoys, cells, ref = _built(w)
    qs = _queries(w, index, decoys, ref, seed)
    rng = np.random.default_rng(seed)
    of = {q: ("reverse", "shuffle")[int(rng.integers(2))] for q in qs}
    got = match_exact(qs, index, decoys, of, starts=cells)
    plain = match_exact(qs, index, decoys, of)
    want = row_set(plain) | credited_rows(qs, index, decoys, of, ref)
    assert row_set(got) == want
    assert set(got["route"]) <= {"exact"} and set(got["n_diff"]) <= {0}
    assert len(row_set(got)) == len(got)  # one row per (mirrored) interval


@SETTINGS
@given(worlds())
def test_p3_no_non_atg_sites_is_a_no_op(w: dict[str, Any]) -> None:
    from dataclasses import replace

    genes = []
    for g in w["gm"].genes:
        x = g.interval.start if g.interval.strand == "+" else g.interval.end - 3
        nt = "ATG" if g.interval.strand == "+" else rc("ATG")
        s = list(w["genome"]["c"])
        for j, b in enumerate(nt):
            s[(x + j) % len(s)] = b
        w["genome"]["c"] = "".join(s)
        genes.append(replace(g, start_codon="ATG"))
    w["gm"] = model(w["reps"], genes)
    index, decoys, cells, ref = _built(w)
    assert cells.sites == ()
    qs = _queries(w, index, decoys, ref, 1)
    of = dict.fromkeys(qs, "reverse")
    pd.testing.assert_frame_equal(
        match_exact(qs, index, decoys, of, starts=cells), match_exact(qs, index, decoys, of)
    )


@SETTINGS
@given(worlds())
def test_p4_decoy_mirror(w: dict[str, Any]) -> None:
    index, decoys, cells, ref = _built(w)
    assert [list(c) for c in cells.target] == [list(c) for c in ref["target"]]
    for kind in ("reverse", "shuffle"):
        assert [list(c) for c in cells.decoys[kind]] == [list(c) for c in ref[kind]]
    rep_len = w["reps"]["c"].length
    for fid, f in enumerate(index.frames):
        n, t = len(f.aa), cells.target[fid]
        for kind in ("reverse", "shuffle"):
            assert len(cells.decoys[kind][fid]) == len(t)
            for k in KS:
                assert usable(kind, decoys[kind][fid], cells.decoys[kind][fid], k, rep_len) == (
                    usable("target", f, t, k, rep_len)
                ), (kind, fid, k)
        for k in KS:
            tw = {(a, b) for _, a, b in windows("target", n, t, k)}
            for _, a, b in windows("reverse", n, cells.decoys["reverse"][fid], k):
                assert mirror_aa("reverse", n, a, b) in tw


def _rc_world(w: dict[str, Any]) -> dict[str, Any]:
    from dataclasses import replace

    from gippa.model.coords import Interval

    n = w["reps"]["c"].length
    genes = []
    for g in w["gm"].genes:
        iv = g.interval
        s, e = n - iv.end, n - iv.start
        if s < 0:
            s, e = s + n, e + n
        genes.append(replace(g, interval=Interval("c", s, e, "-" if iv.strand == "+" else "+")))
    return {**w, "genome": {"c": rc(w["genome"]["c"])}, "gm": model(w["reps"], genes)}


@SETTINGS
@given(worlds())
def test_p5_strand_symmetry(w: dict[str, Any]) -> None:
    w2 = _rc_world(w)
    n = w["reps"]["c"].length

    def credited(world: dict[str, Any]) -> set[tuple[str, str, int, int]]:
        index, decoys, cells, _ = _built(world)
        qs = sorted(
            {
                g.aa[:k]
                for g in world["gm"].genes
                for k in (7, 10)
                if len(g.aa) >= k and "*" not in g.aa[:k] and "X" not in g.aa[:k]
            }
        )
        if not qs:
            return set()
        df = match_exact(qs, index, decoys, dict.fromkeys(qs, "reverse"), starts=cells)
        tgt = df[df["search"] == "target"].reset_index(drop=True)
        fl = ini().initiator_rows(
            tgt, cells.sites, {"c": n}, {"c": world["reps"]["c"].topology == "circular"}
        )
        return {
            (p, s, int(a), int(b))
            for p, s, a, b, f in zip(
                tgt["peptide"], tgt["strand"], tgt["start"], tgt["end"], fl, strict=True
            )
            if f
        }

    a, b = credited(w), credited(w2)

    def flip(r: tuple[str, str, int, int]) -> tuple[str, str, int, int]:
        p, s, x, y = r
        u, v = n - y, n - x
        if u < 0:
            u, v = u + n, v + n
        return (p, "-" if s == "+" else "+", u, v)

    assert {flip(r) for r in a} == b
    assert {r[0] for r in a} == {r[0] for r in b}
