"""Annotated stop sites, frame cells, decoy cells and index format 2 (M7-16a; ADR-0058 items 1-2;
docs/specs/M7-16-stop-codon-readthrough.md rules 1-3 and 8, examples 5, 8, edge cases, "Determinism
of the decoy cells" = statistician review change 7).

Interface pinned here (spec "Public interface")::

    giap.mapping.readthrough.StopSite(contig, strand, codon_start, locus_tags, codon, kind)
    giap.mapping.readthrough.StopCells(sites, target, decoys)
    annotated_stops(gm, genome) -> (sites, n_not_a_stop)
    frame_cells(index, sites) -> tuple[np.ndarray, ...]          # per frame_id, sorted int64
    mirror_cells(kind, frames, cells) -> tuple[np.ndarray, ...]   # ``frames``: that decoy's frames
    site_keys(sites) -> frozenset[(contig, strand, codon_start)]
    giap.mapping.sixframe.READTHROUGH == 22
    save_index(path, index, decoys, cells=None)  -> format 2;  load_stop_cells(path) -> StopCells

``mirror_cells`` is given the decoy frames of ``kind`` (rule 3: "Every mirrored cell must hold
STOP, or the function raises", which only the decoy frames can show); for ``reverse`` the decoy
frame of the decoy is the target frame, so the involution is ``mirror(rev, target,
mirror(rev, decoy, c))``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from helpers.readthrough_index import FIXTURE, SEED, world
from helpers.readthrough_toys import gene, model, rc, toy_c1, toy_c2

from giap.mapping.decoys import build_decoys
from giap.mapping.sixframe import (
    STOP,
    Frame,
    SixFrameIndex,
    build_sixframe,
    load_index,
    save_index,
)
from giap.model import Replicon
from giap.rng import for_step


def _rt():  # type: ignore[no-untyped-def]
    import giap.mapping.readthrough as rt

    return rt


def _cells_of(index: SixFrameIndex, gm, genome):  # type: ignore[no-untyped-def]
    rt = _rt()
    sites, n_not = rt.annotated_stops(gm, genome)
    return sites, n_not, rt.frame_cells(index, sites)


def _as_lists(cells: tuple[np.ndarray, ...]) -> list[list[int]]:
    return [[int(x) for x in c] for c in cells]


def _site(contig: str, strand: str, start: int, tags: tuple[str, ...], codon: str):  # type: ignore[no-untyped-def]
    return _rt().StopSite(
        contig=contig, strand=strand, codon_start=start, locus_tags=tags, codon=codon,
        kind="terminal",
    )  # fmt: skip


def _circular(n: int, edits: dict[int, str], seed: int = 1) -> str:
    rng = np.random.default_rng(seed)
    seq = list("".join(rng.choice(list("ACGT"), size=n)))
    for at, codon in edits.items():
        for i, c in enumerate(codon):
            seq[(at + i) % n] = c
    return "".join(seq)


# ---------------------------------------------------------------- rule 1: sites


def test_readthrough_code_is_22() -> None:
    from giap.mapping.sixframe import AMBIG, NOMATCH, READTHROUGH

    assert (STOP, AMBIG, NOMATCH, READTHROUGH) == (19, 20, 21, 22)


def test_example_1_plus_strand_site_and_cell() -> None:
    genome, reps, gm = toy_c1()
    index = build_sixframe(genome, reps)
    sites, n_not, cells = _cells_of(index, gm, genome)
    assert sites == (_site("c1", "+", 54, ("G1",), "TGA"),)
    assert n_not == 0
    assert _as_lists(cells) == [[18], [], [], [], [], []]
    assert all(c.dtype == np.int64 for c in cells)
    assert index.frames[0].aa[30] == STOP and index.frames[0].aa[41] == STOP  # not marked


def test_example_8_minus_strand_cell() -> None:
    genome, reps, gm = toy_c2()
    index = build_sixframe(genome, reps)
    sites, n_not, cells = _cells_of(index, gm, genome)
    assert sites == (_site("c2", "-", 120, ("G2",), "TAA"),)
    assert n_not == 0
    assert _as_lists(cells) == [[], [], [], [59], [], []]  # frame -0: (300 - 0 - 3 - 120) / 3
    iv = index.to_interval(3, 59, 60)
    assert (iv.contig, iv.start, iv.end, iv.strand) == ("c2", 120, 123, "-")


def test_wrap_plus_stop_crossing_the_origin_has_one_cell() -> None:
    """Edge case: length 300, ext_nt 426; + gene [240, 302): codon_start 299, frame +2, p = 99."""
    seq = _circular(300, {299: "TAA"})
    reps = {"r": Replicon("r", 300, "circular")}
    gm = model(reps, [gene("W1", "r", 240, 302, "+")])
    index = build_sixframe({"r": seq}, reps)
    assert index.frames[2].ext_nt == 426
    sites, _, cells = _cells_of(index, gm, {"r": seq})
    assert sites == (_site("r", "+", 299, ("W1",), "TAA"),)
    assert _as_lists(cells) == [[], [], [99], [], [], []]


def test_wrap_plus_stop_at_zero_has_two_cells() -> None:
    seq = _circular(300, {0: "TGA"})
    reps = {"r": Replicon("r", 300, "circular")}
    gm = model(reps, [gene("W2", "r", 240, 303, "+")])
    index = build_sixframe({"r": seq}, reps)
    sites, _, cells = _cells_of(index, gm, {"r": seq})
    assert sites == (_site("r", "+", 0, ("W2",), "TGA"),)
    assert _as_lists(cells) == [[0, 100], [], [], [], [], []]


def test_wrap_minus_stop_at_zero_has_two_cells() -> None:
    seq = _circular(300, {0: rc("TAG")})
    reps = {"r": Replicon("r", 300, "circular")}
    gm = model(reps, [gene("W3", "r", 0, 90, "-")])
    index = build_sixframe({"r": seq}, reps)
    sites, _, cells = _cells_of(index, gm, {"r": seq})
    assert sites == (_site("r", "-", 0, ("W3",), "TAG"),)
    assert _as_lists(cells) == [[], [], [], [41, 141], [], []]
    for p in (41, 141):  # rule 2's round trip, start wrapped
        iv = index.to_interval(3, p, p + 1)
        assert (iv.start % 300, iv.end - iv.start, iv.strand) == (0, 3, "-")


def test_shared_stop_is_one_site_with_sorted_tags_and_one_cell() -> None:
    nt = "ATG" + "GCT" * 28 + "TAA" + "GCT" * 20
    reps = {"s": Replicon("s", len(nt), "linear")}
    gm = model(reps, [gene("G5", "s", 30, 90, "+"), gene("G4", "s", 0, 90, "+")])
    index = build_sixframe({"s": nt}, reps)
    sites, _, cells = _cells_of(index, gm, {"s": nt})
    assert sites == (_site("s", "+", 87, ("G4", "G5"), "TAA"),)
    assert _as_lists(cells)[0] == [29]
    assert sum(len(c) for c in cells) == 1


def test_world_sites_sorted_sorf_included_pseudogene_excluded_not_a_stop_counted() -> None:
    genome, reps, gm = world()
    index = build_sixframe(genome, reps)
    sites, n_not, cells = _cells_of(index, gm, genome)
    assert sites == (
        _site("a", "+", 297, ("A1",), "TAA"),
        _site("a", "+", 897, ("A2", "A3"), "TGA"),
        _site("a", "+", 2097, ("A5",), "TAA"),  # sORF
        _site("a", "-", 1500, ("A4",), "TAG"),
        _site("b", "+", 1199, ("B1",), "TAG"),  # crosses the origin
        _site("b", "-", 3, ("B2",), "TAA"),
    )
    assert n_not == 1  # A7 ends in GCT; the pseudogene A6 is neither a site nor counted
    want = [[] for _ in range(12)]
    want[0] = [99, 299, 699]
    want[3] = [499]
    want[8] = [399]
    want[9] = [40, 440]
    assert _as_lists(cells) == want
    assert _rt().site_keys(sites) == frozenset(
        {("a", "+", 297), ("a", "+", 897), ("a", "+", 2097), ("a", "-", 1500), ("b", "+", 1199),
         ("b", "-", 3)}
    )  # fmt: skip


def test_partial_gene_included_iff_its_last_codon_is_a_stop() -> None:
    nt = "GCT" * 10 + "TGA" + "GCT" * 10 + "AAA"
    reps = {"p": Replicon("p", len(nt), "linear")}
    gm = model(reps, [gene("P1", "p", 0, 33, "+", partial=True),
                      gene("P2", "p", 33, len(nt), "+", partial=True)])  # fmt: skip
    index = build_sixframe({"p": nt}, reps)
    sites, n_not, _ = _cells_of(index, gm, {"p": nt})
    assert sites == (_site("p", "+", 30, ("P1",), "TGA"),)
    assert n_not == 1


def test_empty_annotation_gives_no_sites() -> None:
    genome, reps, _ = toy_c1()
    gm = model(reps, [])
    index = build_sixframe(genome, reps)
    sites, n_not, cells = _cells_of(index, gm, genome)
    assert sites == () and n_not == 0
    assert all(len(c) == 0 for c in cells) and len(cells) == 6


# ---------------------------------------------------------------- rule 2: the round trip raises


def test_frame_cells_raises_on_a_site_that_is_not_a_stop_in_its_frame() -> None:
    genome, reps, _ = toy_c1()
    index = build_sixframe(genome, reps)
    bad = _site("c1", "+", 51, ("G1",), "TGA")  # p = 17 is K, not '*'
    with pytest.raises(ValueError, match="c1"):
        _rt().frame_cells(index, (bad,))


# ---------------------------------------------------------------- rule 3: decoy cells


def test_example_5_reverse_and_shuffle_cells() -> None:
    genome, reps, gm = toy_c1()
    index = build_sixframe(genome, reps)
    _, _, cells = _cells_of(index, gm, genome)
    rt = _rt()
    rev = build_decoys(index, "reverse", 50, for_step("sixframe_index", 1))
    shf = build_decoys(index, "shuffle", 50, for_step("sixframe_index", 1))
    rc_cells = rt.mirror_cells("reverse", rev, cells)
    sh_cells = rt.mirror_cells("shuffle", shf, cells)
    assert _as_lists(rc_cells)[0] == [23]  # 42 - 1 - 18
    assert rev[0].aa[23] == STOP and rev[0].aa[11] == STOP  # 30 -> 11 holds '*' but is unmarked
    assert _as_lists(sh_cells)[0] == [18]
    assert all(shf[0].aa[p] == STOP for p in (18, 30, 41))
    assert [len(c) for c in rc_cells] == [len(c) for c in cells] == [len(c) for c in sh_cells]


def test_mirror_cells_reverse_is_an_involution_and_shuffle_the_identity() -> None:
    genome, reps, gm = world()
    index = build_sixframe(genome, reps)
    _, _, cells = _cells_of(index, gm, genome)
    rt = _rt()
    rev = build_decoys(index, "reverse", 50, for_step("sixframe_index", 2))
    shf = build_decoys(index, "shuffle", 50, for_step("sixframe_index", 2))
    once = rt.mirror_cells("reverse", rev, cells)
    assert _as_lists(rt.mirror_cells("reverse", index.frames, once)) == _as_lists(cells)
    assert _as_lists(rt.mirror_cells("shuffle", shf, cells)) == _as_lists(cells)
    for f, c in zip(rev, once, strict=True):
        assert np.all(f.aa[c] == STOP) and list(c) == sorted(c)


def test_mirror_cells_raises_when_a_mirrored_cell_is_not_a_stop() -> None:
    genome, reps, gm = toy_c1()
    index = build_sixframe(genome, reps)
    _, _, cells = _cells_of(index, gm, genome)
    broken = list(build_decoys(index, "shuffle", 50, for_step("sixframe_index", 1)))
    aa = broken[0].aa.copy()
    aa[18] = 0
    broken[0] = Frame(broken[0].contig, broken[0].strand, broken[0].offset, aa, broken[0].ext_nt)
    with pytest.raises(ValueError):
        _rt().mirror_cells("shuffle", tuple(broken), cells)


def test_mirror_cells_draws_no_randomness(monkeypatch: pytest.MonkeyPatch) -> None:
    """Review change 7 / §6.6: no call to ``giap.rng.for_step`` or ``numpy.random``, the legacy
    global state is untouched, and the step's shuffle stream gives the same next draw with and
    without the call."""
    import giap.rng

    genome, reps, gm = world()
    index = build_sixframe(genome, reps)
    _, _, cells = _cells_of(index, gm, genome)
    rt = _rt()
    calls: list[str] = []
    real_for_step = giap.rng.for_step
    real_default_rng = np.random.default_rng

    def spy_for_step(*a, **k):  # type: ignore[no-untyped-def]
        calls.append("for_step")
        return real_for_step(*a, **k)

    def spy_default_rng(*a, **k):  # type: ignore[no-untyped-def]
        calls.append("default_rng")
        return real_default_rng(*a, **k)

    g1 = for_step("sixframe_index", SEED)
    shf1 = build_decoys(index, "shuffle", 50, g1)
    monkeypatch.setattr(giap.rng, "for_step", spy_for_step)
    monkeypatch.setattr(rt, "for_step", spy_for_step, raising=False)
    monkeypatch.setattr(np.random, "default_rng", spy_default_rng)
    legacy = np.random.get_state()  # noqa: NPY002 - the legacy global state must not move
    rt.mirror_cells("shuffle", shf1, cells)
    rt.mirror_cells("reverse", build_decoys(index, "reverse", 50, g1), cells)
    after = np.random.get_state()  # noqa: NPY002
    monkeypatch.undo()
    assert calls == []
    assert legacy[0] == after[0] and np.array_equal(legacy[1], after[1]) and legacy[2:] == after[2:]
    g2 = for_step("sixframe_index", SEED)
    shf2 = build_decoys(index, "shuffle", 50, g2)
    assert all(np.array_equal(a.aa, b.aa) for a, b in zip(shf1, shf2, strict=True))
    assert int(g1.integers(0, 2**62)) == int(g2.integers(0, 2**62))


# ---------------------------------------------------------------- rule 8: index format 2


def _stop_cells(index, decoys, gm, genome):  # type: ignore[no-untyped-def]
    rt = _rt()
    sites, _ = rt.annotated_stops(gm, genome)
    target = rt.frame_cells(index, sites)
    return rt.StopCells(
        sites=sites, target=target,
        decoys={k: rt.mirror_cells(k, fr, target) for k, fr in decoys.items()},
    )  # fmt: skip


def test_save_index_writes_format_2_and_load_stop_cells_round_trips(tmp_path: Path) -> None:
    from giap.mapping.sixframe import load_stop_cells

    genome, reps, gm = world()
    index = build_sixframe(genome, reps)
    decoys = {
        "reverse": build_decoys(index, "reverse", 50, for_step("sixframe_index", 3)),
        "shuffle": build_decoys(index, "shuffle", 50, for_step("sixframe_index", 3)),
    }
    cells = _stop_cells(index, decoys, gm, genome)
    path = tmp_path / "six.npz"
    save_index(path, index, decoys, cells)
    with np.load(path) as npz:
        meta = json.loads(str(npz["meta"]))
        keys = set(npz.keys())
    assert meta["format"] == 2
    assert len(meta["stop_sites"]) == 6
    want = {"target_rt", "target_rt_bounds", "reverse_rt", "reverse_rt_bounds", "shuffle_rt",
            "shuffle_rt_bounds"}  # fmt: skip
    assert want <= keys
    got = load_stop_cells(path)
    assert got.sites == cells.sites
    assert _as_lists(got.target) == _as_lists(cells.target)
    assert set(got.decoys) == {"reverse", "shuffle"}
    for k in decoys:
        assert _as_lists(got.decoys[k]) == _as_lists(cells.decoys[k])
    index2, decoys2 = load_index(path)  # signature unchanged
    assert all(np.array_equal(a.aa, b.aa) for a, b in zip(index2.frames, index.frames, strict=True))
    assert set(decoys2) == {"reverse", "shuffle"}


def test_save_index_without_cells_writes_empty_cells(tmp_path: Path) -> None:
    from giap.mapping.sixframe import load_stop_cells

    genome, reps, _ = toy_c1()
    index = build_sixframe(genome, reps)
    decoys = {"reverse": build_decoys(index, "reverse", 50, for_step("sixframe_index", 1))}
    path = tmp_path / "six.npz"
    save_index(path, index, decoys)
    got = load_stop_cells(path)
    assert got.sites == ()
    assert [len(c) for c in got.target] == [0] * 6
    assert [len(c) for c in got.decoys["reverse"]] == [0] * 6


def _format1_file(path: Path) -> Path:
    genome, reps, _ = toy_c1()
    index = build_sixframe(genome, reps)
    rev = build_decoys(index, "reverse", 50, for_step("sixframe_index", 1))
    meta = {"format": 1, "replicons": [{"id": "c1", "length": 126, "topology": "linear",
                                        "kind": "unknown"}],
            "frames": [[f.contig, f.strand, f.offset, f.ext_nt] for f in index.frames],
            "decoys": ["reverse"]}  # fmt: skip

    def pack(frames: tuple[Frame, ...]) -> tuple[np.ndarray, np.ndarray]:
        bounds = np.zeros(len(frames) + 1, dtype=np.int64)
        np.cumsum([len(f.aa) for f in frames], out=bounds[1:])
        return np.concatenate([f.aa for f in frames]).astype(np.uint8), bounds

    ta, tb = pack(index.frames)
    ra, rb = pack(rev)
    with path.open("wb") as fh:
        np.savez_compressed(fh, meta=np.array(json.dumps(meta)), target_aa=ta, target_bounds=tb,
                            reverse_aa=ra, reverse_bounds=rb)  # fmt: skip
    return path


@pytest.mark.parametrize("which", ["crafted", "fixture"])
def test_load_stop_cells_rejects_format_1(tmp_path: Path, which: str) -> None:
    from giap.mapping.sixframe import load_stop_cells

    path = _format1_file(tmp_path / "old.npz") if which == "crafted" else FIXTURE
    with pytest.raises(ValueError, match="rebuild the six-frame index"):
        load_stop_cells(path)


def test_step_decoy_frames_are_byte_identical_to_format_1(tmp_path: Path) -> None:
    """Review change 7: for the same seed the ``sixframe_index`` step's format-2 file holds the
    target and decoy frames of the committed format-1 file (written on main), byte for byte."""
    from helpers.readthrough_index import run_sixframe_index

    out = run_sixframe_index(tmp_path)
    with np.load(FIXTURE) as old, np.load(out) as new:
        m_old, m_new = json.loads(str(old["meta"])), json.loads(str(new["meta"]))
        assert (m_old["format"], m_new["format"]) == (1, 2)
        assert m_new["frames"] == m_old["frames"] and m_new["decoys"] == m_old["decoys"]
        for key in ("target_aa", "target_bounds", "reverse_aa", "reverse_bounds", "shuffle_aa",
                    "shuffle_bounds"):  # fmt: skip
            a, b = old[key], new[key]
            assert a.dtype == b.dtype and a.tobytes() == b.tobytes(), key


# ---------------------------------------------------------------- review round 1 (PR #331)


def test_linear_contig_stop_in_the_last_codon_is_a_site() -> None:
    """Edge case (item 3): a ``+`` gene ``[96, 126)`` on the linear 126-nt ``c1`` ends in the
    contig's last codon (TAG, p = 41): ``codon_start + 3 == n`` is in range, so it is a site."""
    genome, reps, gm0 = toy_c1()
    gm = model(reps, [*gm0.genes, gene("G9", "c1", 96, 126, "+")])
    index = build_sixframe(genome, reps)
    sites, n_not, cells = _cells_of(index, gm, genome)
    assert sites == (_site("c1", "+", 54, ("G1",), "TGA"), _site("c1", "+", 123, ("G9",), "TAG"))
    assert n_not == 0
    assert _as_lists(cells) == [[18, 41], [], [], [], [], []]


def test_linear_contig_minus_gene_at_zero_is_a_site() -> None:
    """Edge case (item 3): a ``-`` gene ``[0, 30)`` on a linear contig: its stop is ``[0, 3)``."""
    gene_strand = "ATG" + "GCT" * 8 + "TAA"
    nt = rc(gene_strand) + "GCT" * 10
    reps = {"m": Replicon("m", len(nt), "linear")}
    gm = model(reps, [gene("M1", "m", 0, 30, "-")])
    index = build_sixframe({"m": nt}, reps)
    sites, n_not, cells = _cells_of(index, gm, {"m": nt})
    assert sites == (_site("m", "-", 0, ("M1",), "TAA"),)
    assert n_not == 0
    assert sum(len(c) for c in cells) == 1
    fid = next(i for i, c in enumerate(cells) if len(c))
    p = int(cells[fid][0])
    assert index.frames[fid].strand == "-" and index.frames[fid].aa[p] == STOP
    iv = index.to_interval(fid, p, p + 1)
    assert (iv.contig, iv.start, iv.end, iv.strand) == ("m", 0, 3, "-")


def test_shared_stop_tags_are_sorted_when_the_features_arrive_unsorted() -> None:
    """Rule 1 (item 6): ``Z9`` (start 0) comes before ``A1`` (start 30) in the gene model; the
    site's ``locus_tags`` are still sorted."""
    nt = "ATG" + "GCT" * 28 + "TAA" + "GCT" * 20
    reps = {"s": Replicon("s", len(nt), "linear")}
    gm = model(reps, [gene("A1", "s", 30, 90, "+"), gene("Z9", "s", 0, 90, "+")])
    assert [g.locus_tag for g in gm.genes] == ["Z9", "A1"]
    index = build_sixframe({"s": nt}, reps)
    sites, _, cells = _cells_of(index, gm, {"s": nt})
    assert sites == (_site("s", "+", 87, ("A1", "Z9"), "TAA"),)
    assert _as_lists(cells)[0] == [29]


@pytest.mark.parametrize("which", ["crafted", "fixture"])
def test_load_index_rejects_format_1(tmp_path: Path, which: str) -> None:
    """Item 4: ``load_index`` itself checks the format (not only ``load_stop_cells``)."""
    path = _format1_file(tmp_path / "old.npz") if which == "crafted" else FIXTURE
    with pytest.raises(ValueError, match="rebuild the six-frame index"):
        load_index(path)
