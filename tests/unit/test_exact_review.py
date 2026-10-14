"""Lead decisions after code review of PR #27 (docs/specs/M3-04-exact-matcher.md).

1. ``match_exact`` validates ``decoy_of`` up front (``ValueError``).
2. I/L twins sharing one automaton key keep their own decoy types; the replicon-length guard on a
   short circular replicon; a ``U`` query is never added to the automaton.
3. ``benchmarks`` is a package, so ``tests/unit/test_exact.py`` and ``benchmarks/test_exact.py``
   are collected together without an import-file mismatch.
"""

from __future__ import annotations

import os
import subprocess
import sys
import types
from pathlib import Path
from typing import ClassVar

import ahocorasick
import pandas as pd
import pytest

import giap.mapping.exact as exact_mod
from giap.mapping.decoys import build_decoys
from giap.mapping.exact import match_exact
from giap.mapping.sixframe import Frame, SixFrameIndex, build_sixframe, decode
from giap.model import Interval, Replicon
from giap.rng import for_step

REPO = Path(__file__).resolve().parents[2]
G = "ATGGCTAAAATTCTGTGGTAACGTNNA"  # +0 = MAKLLW*RX (M3-02 genome)
# +0 = MAKLLW*WLLKAM*, so the reversed +0 (*MAKLLW*WLLKAM) contains MAKLLW as well.
TWIN = "ATGGCTAAACTGCTGTGG" + "TAA" + "TGGCTGCTGAAAGCTATG" + "TAA"
# 30-nt circular replicon whose +0 frame has no stop: MAKLLWGSPE, repeated by the wrap.
SHORT = "ATGGCTAAACTGCTGTGGGGTTCTCCGGAA"


def _index(genome: dict[str, str], topology: str = "linear", wrap_nt: int = 120) -> SixFrameIndex:
    reps = {c: Replicon(c, len(s), topology) for c, s in genome.items()}  # type: ignore[arg-type]
    return build_sixframe(genome, reps, wrap_nt=wrap_nt)


def _identity_decoys(idx: SixFrameIndex) -> tuple[Frame, ...]:
    """A "shuffle" decoy set that leaves every frame unchanged, so its hits are known exactly."""
    return tuple(Frame(f.contig, f.strand, f.offset, f.aa.copy(), f.ext_nt) for f in idx.frames)


def _rev(idx: SixFrameIndex) -> tuple[Frame, ...]:
    return build_decoys(idx, "reverse", 50, for_step("sixframe_index", 0))


# --- 1. decoy_of validation --------------------------------------------------------------------


def test_decoy_type_without_decoy_frames_raises() -> None:
    idx = _index({"c1": G})
    with pytest.raises(ValueError, match="shuffle"):
        match_exact(["MAKLLW"], idx, {"reverse": _rev(idx)}, {"MAKLLW": "shuffle"}, min_len=4)


def test_decoy_type_without_frames_raises_even_with_no_target_hit() -> None:
    idx = _index({"c1": G})
    with pytest.raises(ValueError, match="shuffle"):
        match_exact(["PPPPPP"], idx, {"reverse": _rev(idx)}, {"PPPPPP": "shuffle"}, min_len=4)


def test_peptide_missing_from_decoy_of_raises_value_error() -> None:
    idx = _index({"c1": G})
    with pytest.raises(ValueError, match="MAKILW"):
        match_exact(
            ["MAKLLW", "MAKILW"], idx, {"reverse": _rev(idx)}, {"MAKLLW": "reverse"}, min_len=4
        )


def test_peptide_missing_from_decoy_of_raises_even_without_hits() -> None:
    idx = _index({"c1": G})
    with pytest.raises(ValueError, match="PPPPPP"):
        match_exact(["PPPPPP"], idx, {"reverse": _rev(idx)}, {}, min_len=4)


def test_valid_decoy_of_does_not_raise() -> None:
    idx = _index({"c1": G})
    df = match_exact(["MAKLLW"], idx, {"reverse": _rev(idx)}, {"MAKLLW": "reverse"}, min_len=4)
    assert list(df["search"]) == ["target"]


# --- 2a. I/L twins sharing one automaton key ----------------------------------------------------


def test_il_twins_each_get_only_their_own_decoy_type() -> None:
    idx = _index({"c1": TWIN})
    assert decode(idx.frames[0].aa).startswith("MAKLLW*WLLKAM*")  # precondition
    decoys = {"reverse": _rev(idx), "shuffle": _identity_decoys(idx)}
    df = match_exact(
        ["MAKILW", "MAKLLW"],
        idx,
        decoys,
        {"MAKILW": "reverse", "MAKLLW": "shuffle"},
        min_len=4,
    )
    assert set(df["query"]) == {"MAKLLW"}
    for pep, kind in (("MAKILW", "reverse"), ("MAKLLW", "shuffle")):
        rows = df[df["peptide"] == pep]
        tgt = rows[rows["search"] == "target"]
        dec = rows[rows["search"] == "decoy"]
        assert len(tgt) == 1, pep
        assert (int(tgt["start"].iloc[0]), int(tgt["end"].iloc[0])) == (0, 18)
        assert len(dec) >= 1, pep
        assert set(dec["decoy_type"]) == {kind}, pep
    # reverse: MAKLLW at reversed +0 [1, 7) mirrors to target window [7, 13) = WLLKAM, nt 21..39
    rev = df[(df["peptide"] == "MAKILW") & (df["search"] == "decoy")]
    assert ("c1", 21, 39, "+") in {
        (r.contig, int(r.start), int(r.end), r.strand) for r in rev.itertuples(index=False)
    }
    # identity "shuffle": the decoy window is the target window itself
    shu = df[(df["peptide"] == "MAKLLW") & (df["search"] == "decoy")]
    assert ("c1", 0, 18, "+") in {
        (r.contig, int(r.start), int(r.end), r.strand) for r in shu.itertuples(index=False)
    }


def test_il_twins_target_rows_are_identical_but_for_peptide() -> None:
    idx = _index({"c1": TWIN})
    decoys = {"reverse": _rev(idx), "shuffle": _identity_decoys(idx)}
    df = match_exact(
        ["MAKILW", "MAKLLW"], idx, decoys, {"MAKILW": "reverse", "MAKLLW": "shuffle"}, min_len=4
    )
    tgt = df[df["search"] == "target"].drop(columns="peptide").reset_index(drop=True)
    pd.testing.assert_frame_equal(
        tgt.iloc[[0]].reset_index(drop=True), tgt.iloc[[1]].reset_index(drop=True)
    )


# --- 2b. replicon-length guard on a ~30-nt circular replicon ----------------------------------


def test_short_circular_replicon_precondition() -> None:
    idx = _index({"r": SHORT}, "circular")
    assert idx.replicons["r"].length == 30
    assert "MAKLLWGSPEMAKLLWGSPE" in decode(idx.frames[0].aa)
    with pytest.raises(ValueError):
        idx.to_interval(0, 0, 11)


@pytest.mark.parametrize("kind", ["reverse", "shuffle"])
def test_peptide_longer_than_replicon_never_raises_nor_hits(kind: str) -> None:
    idx = _index({"r": SHORT}, "circular")
    longer = ["MAKLLWGSPEM", "MAKLLWGSPEMA", "AKLLWGSPEMAKLLWGSPEM", "EPSGWLLKAMEP"]
    decoys = {"reverse": _rev(idx), "shuffle": _identity_decoys(idx)}
    df = match_exact(longer, idx, decoys, dict.fromkeys(longer, kind), min_len=4)  # type: ignore[arg-type]
    assert len(df) == 0


def test_peptide_as_long_as_replicon_maps_once_including_origin_rotation() -> None:
    idx = _index({"r": SHORT}, "circular")
    peps = ["MAKLLWGSPE", "SPEMAKLLWG"]  # the second crosses the origin
    df = match_exact(peps, idx, {"reverse": _rev(idx)}, dict.fromkeys(peps, "reverse"), min_len=4)
    tgt = df[df["search"] == "target"]
    got = {
        r.peptide: Interval(r.contig, int(r.start), int(r.end), r.strand)
        for r in tgt.itertuples(index=False)
    }
    assert len(tgt) == 2
    assert got["MAKLLWGSPE"] == Interval("r", 0, 30, "+")
    assert got["SPEMAKLLWG"] == Interval("r", 21, 51, "+")


def test_long_and_fitting_peptides_mixed_in_one_call() -> None:
    idx = _index({"r": SHORT}, "circular")
    peps = ["MAKLLWGSPE", "MAKLLWGSPEMA", "AKLLWGS"]
    df = match_exact(peps, idx, {"reverse": _rev(idx)}, dict.fromkeys(peps, "reverse"), min_len=4)
    tgt = df[df["search"] == "target"]
    assert sorted(tgt["peptide"]) == ["AKLLWGS", "MAKLLWGSPE"]
    assert "MAKLLWGSPEMA" not in set(df["peptide"])


# --- 2c. a U query is never added to the automaton --------------------------------------------


class _Recording:
    """Delegating wrapper around ``ahocorasick.Automaton`` (not subclassable) that logs words."""

    words: ClassVar[list[str]] = []

    def __init__(self, *args: object) -> None:
        self._a = ahocorasick.Automaton(*args)

    def add_word(self, word: str, value: object) -> bool:
        _Recording.words.append(word)
        return bool(self._a.add_word(word, value))

    def __getattr__(self, name: str) -> object:
        return getattr(self._a, name)


def test_u_query_is_never_added_to_the_automaton(monkeypatch: pytest.MonkeyPatch) -> None:
    _Recording.words = []
    fake = types.SimpleNamespace(Automaton=_Recording, STORE_INTS=ahocorasick.STORE_INTS)
    monkeypatch.setattr(exact_mod, "ahocorasick", fake)
    idx = _index({"c1": G})
    peps = ["MAKLLW", "MAKUUW", "UAKLLW"]
    df = match_exact(peps, idx, {"reverse": _rev(idx)}, dict.fromkeys(peps, "reverse"), min_len=4)
    assert _Recording.words == ["MAKLLW"]
    assert set(df["peptide"]) == {"MAKLLW"}


def test_u_only_queries_return_empty_frame() -> None:
    idx = _index({"c1": G})
    peps = ["MAKLUW", "UUUUUU"]
    df = match_exact(peps, idx, {"reverse": _rev(idx)}, dict.fromkeys(peps, "reverse"), min_len=4)
    assert len(df) == 0


# --- 3. benchmarks is a package ------------------------------------------------------------------


def test_benchmarks_is_a_package() -> None:
    assert (REPO / "benchmarks" / "__init__.py").is_file()


def test_unit_and_benchmark_test_exact_collect_together() -> None:
    env = {k: v for k, v in os.environ.items() if not k.startswith("PYTEST_")}
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "tests/unit/test_exact.py",
            "benchmarks/test_exact.py",
            "--collect-only",
            "-q",
            "-p",
            "no:cacheprovider",
            "-p",
            "no:xdist",
        ],
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
    )
    out = proc.stdout + proc.stderr
    assert proc.returncode == 0, out
    assert "import file mismatch" not in out
    assert "tests/unit/test_exact.py::" in out
    assert "benchmarks/test_exact.py::" in out
