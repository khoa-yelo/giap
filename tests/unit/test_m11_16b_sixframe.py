"""M11-16b six-frame translation per table and index format 3 (docs/specs/M11-16-translation-
tables-4-25.md item 4, examples 1, 2 and 4, AC 7, edge cases; decision D3 and D13).

Interface (spec item 4)::

    translate(nt: bytes, table: int = 11) -> np.ndarray          # codes, frame 0
    translate_text(nt: str, table: int = 11, offset: int = 0) -> str
    translate11 / translate11_text                                # aliases, table 11 (D13)
    build_sixframe(genome, replicons, wrap_nt=WRAP_NT, table=11) -> SixFrameIndex
    SixFrameIndex.translation_table: int = 11
    save_index(...)            # format 3, meta["translation_table"] = index.translation_table
    load_index(path, *, table) / load_stop_cells(path, *, table)  # table is a required keyword

Table 11 byte identity: ``tests/data/mapping/readthrough/sixframe_format2.npz`` is the work file
that ``sixframe_index`` wrote for ``helpers.readthrough_index.world()`` on main before M11-16
(``make_format2_fixture.py``). A table-11 index of today's code must hold the same arrays, byte for
byte, and the same meta apart from ``format`` and ``translation_table``.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pytest
from helpers.readthrough_index import SEED, run_sixframe_index, world

from gippa.mapping import sixframe as sf
from gippa.mapping.decoys import build_decoys
from gippa.mapping.sixframe import STOP, build_sixframe, decode, save_index
from gippa.model import Replicon
from gippa.rng import for_step

TABLES = (4, 11, 25)
FORMAT2 = Path(__file__).resolve().parents[1] / "data/mapping/readthrough/sixframe_format2.npz"
FORMAT1 = Path(__file__).resolve().parents[1] / "data/mapping/readthrough/sixframe_format1.npz"
EX2 = "ATGGCCTGAGCCTAA"
EX2_RC = "TTAGGCTCAGGCCAT"
P0, M0 = 0, 3


def _mismatch(path: Path, built: int, run: int) -> str:
    return re.escape(
        f"{path}: six-frame index built for translation table {built}; this run uses {run}; "
        "rebuild the six-frame index"
    )


def _linear(seq: str, table: int, cid: str = "c") -> sf.SixFrameIndex:
    return build_sixframe({cid: seq}, {cid: Replicon(cid, len(seq), "linear")}, table=table)


def _save(path: Path, index: sf.SixFrameIndex) -> Path:
    save_index(path, index, {"reverse": build_decoys(index, "reverse", 50, for_step("t", 1))})
    return path


def _meta(path: Path) -> dict[str, object]:
    with np.load(path) as npz:
        meta: dict[str, object] = json.loads(str(npz["meta"]))
    return meta


# ---------------------------------------------------------------- example 1: translate_text


@pytest.mark.parametrize(("table", "want"), [(11, "MW**"), (4, "MWW*"), (25, "MWG*")])
def test_example_1_translate_text(table: int, want: str) -> None:
    assert sf.translate_text("ATGTGGTGATAA", table) == want
    assert sf.translate_text("ATGTGGTGATAA", table=table) == want


@pytest.mark.parametrize("table", TABLES)
def test_example_1_ambiguous_codon_and_partial_codon(table: int) -> None:
    assert sf.translate_text("ATGTGN", table) == "MX"
    assert sf.translate_text("ATGTGATG", table) == sf.translate_text("ATGTGA", table)
    assert sf.translate_text("AT", table) == ""


def test_translate_text_offset_is_the_third_argument() -> None:
    assert sf.translate_text("CATGTGATAA", 4, 1) == "MW*"
    assert sf.translate_text("CATGTGATAA", table=25, offset=1) == "MG*"


@pytest.mark.parametrize(("table", "tga"), [(11, STOP), (4, 17), (25, 5)])
def test_translate_codes_tga(table: int, tga: int) -> None:
    """Item 4: TGA is code 19 (``*``) under 11, ``W`` (17) under 4 and ``G`` (5) under 25."""
    assert sf.ALPHABET[17] == "W" and sf.ALPHABET[5] == "G"
    assert sf.translate(b"TGATAATAG", table).tolist() == [tga, STOP, STOP]


def test_d13_aliases_are_table_11() -> None:
    nt = world()[0]["a"]
    assert sf.translate(nt.encode()).tobytes() == sf.translate11(nt.encode()).tobytes()
    assert sf.translate(nt.encode(), 11).tobytes() == sf.translate11(nt.encode()).tobytes()
    assert sf.translate_text(nt) == sf.translate11_text(nt) == sf.translate_text(nt, 11)
    assert sf.translate_text(nt, 11, 2) == sf.translate11_text(nt, 2)


@pytest.mark.parametrize("table", [1, 2, 12, 26])
def test_unsupported_table_raises(table: int) -> None:
    with pytest.raises(ValueError, match="unsupported translation table"):
        sf.translate_text("ATG", table)
    with pytest.raises(ValueError, match="unsupported translation table"):
        _linear(EX2, table)


# ---------------------------------------------------------------- example 2: six frames


@pytest.mark.parametrize(("table", "want"), [(11, "MA*A*"), (4, "MAWA*"), (25, "MAGA*")])
def test_example_2_six_frames_both_strands(table: int, want: str) -> None:
    plus = _linear(EX2, table)
    assert decode(plus.frames[P0].aa) == want
    minus = _linear(EX2_RC, table)
    assert decode(minus.frames[M0].aa) == want
    iv = minus.to_interval(M0, 2, 3)
    assert (iv.contig, iv.start, iv.end, iv.strand) == ("c", 6, 9, "-")
    assert EX2_RC[6:9] == "TCA"
    assert decode(minus.frames[M0].aa[2:3]) == ("*" if table == 11 else want[2])


@pytest.mark.parametrize("table", TABLES)
def test_index_records_its_table(table: int) -> None:
    assert _linear(EX2, table).translation_table == table


def test_index_default_table_is_11() -> None:
    idx = build_sixframe({"c": EX2}, {"c": Replicon("c", 15, "linear")})
    assert idx.translation_table == 11
    assert sf.SixFrameIndex(replicons={}, frames=()).translation_table == 11


@pytest.mark.parametrize(("table", "letter"), [(4, "W"), (25, "G"), (11, "*")])
def test_circular_tga_across_the_origin(table: int, letter: str) -> None:
    """Edge case: a TGA over bases n-2, n-1, 0 of a circular replicon reads W / G in the wrapped
    frame (``+``, codon at 27 of a 29-nt contig with frame offset 0)."""
    seq = "A" + "GCT" * 8 + "AATG"  # 29 nt: [27:29] = "TG", [0] = "A" -> TGA across the origin
    assert len(seq) == 29 and (seq[27:] + seq[0]) == "TGA"
    idx = build_sixframe({"r": seq}, {"r": Replicon("r", 29, "circular")}, wrap_nt=12, table=table)
    f = next(i for i, fr in enumerate(idx.frames) if fr.strand == "+" and fr.offset == 0)
    assert decode(idx.frames[f].aa[9:10]) == letter
    iv = idx.to_interval(f, 9, 10)
    assert (iv.start, iv.end) == (27, 30)


@pytest.mark.parametrize("table", (4, 25))
def test_linear_partial_tga_is_dropped(table: int) -> None:
    """Edge case: a linear contig ending in ``...TG`` gives no codon there, in every table."""
    idx = _linear("ATGGCTTG", table)
    assert decode(idx.frames[P0].aa) == "MA"


# ---------------------------------------------------------------- example 4 and AC 7: format 3


@pytest.mark.parametrize("table", TABLES)
def test_ac7_format_3_round_trips_at_each_table(tmp_path: Path, table: int) -> None:
    index = _linear(EX2 * 3, table)
    path = _save(tmp_path / "six.npz", index)
    meta = _meta(path)
    assert meta["format"] == 3 and meta["translation_table"] == table
    got, decoys = sf.load_index(path, table=table)
    assert got.translation_table == table
    assert [decode(f.aa) for f in got.frames] == [decode(f.aa) for f in index.frames]
    assert set(decoys) == {"reverse"}
    cells = sf.load_stop_cells(path, table=table)
    assert cells.sites == ()


@pytest.mark.parametrize(("built", "run"), [(4, 11), (11, 4), (25, 4), (4, 25), (11, 25)])
def test_example_4_cross_table_load_is_refused(tmp_path: Path, built: int, run: int) -> None:
    path = _save(tmp_path / "six.npz", _linear(EX2, built))
    with pytest.raises(ValueError, match=_mismatch(path, built, run)):
        sf.load_index(path, table=run)
    with pytest.raises(ValueError, match=_mismatch(path, built, run)):
        sf.load_stop_cells(path, table=run)


def test_table_is_a_required_keyword(tmp_path: Path) -> None:
    path = _save(tmp_path / "six.npz", _linear(EX2, 11))
    for fn in (sf.load_index, sf.load_stop_cells):
        with pytest.raises(TypeError):
            fn(path)  # type: ignore[call-arg]
        with pytest.raises(TypeError):
            fn(path, 11)  # type: ignore[misc]


def test_d3_format_2_loads_only_at_table_11() -> None:
    assert _meta(FORMAT2)["format"] == 2 and "translation_table" not in _meta(FORMAT2)
    index, decoys = sf.load_index(FORMAT2, table=11)
    assert index.translation_table == 11 and set(decoys) == {"reverse", "shuffle"}
    assert len(sf.load_stop_cells(FORMAT2, table=11).sites) == 6
    for run in (4, 25):
        with pytest.raises(ValueError, match=_mismatch(FORMAT2, 11, run)):
            sf.load_index(FORMAT2, table=run)
        with pytest.raises(ValueError, match=_mismatch(FORMAT2, 11, run)):
            sf.load_stop_cells(FORMAT2, table=run)


@pytest.mark.parametrize("table", TABLES)
def test_format_1_is_still_refused_at_every_table(table: int) -> None:
    for fn in (sf.load_index, sf.load_stop_cells):
        with pytest.raises(ValueError, match="format 1 has no stop cells; rebuild"):
            fn(FORMAT1, table=table)


def test_format_3_without_a_table_key_is_unsupported(tmp_path: Path) -> None:
    """Edge case: the key is required in format 3 (a hand-edited file)."""
    path = _save(tmp_path / "six.npz", _linear(EX2, 4))
    with np.load(path) as npz:
        arrays = {k: npz[k] for k in npz.files}
    meta = json.loads(str(arrays["meta"]))
    del meta["translation_table"]
    arrays["meta"] = np.array(json.dumps(meta))
    bad = tmp_path / "bad.npz"
    with bad.open("wb") as fh:
        np.savez_compressed(fh, **arrays)  # type: ignore[arg-type]
    for fn in (sf.load_index, sf.load_stop_cells):
        with pytest.raises(ValueError, match="unsupported sixframe_index format"):
            fn(bad, table=4)
        with pytest.raises(ValueError, match="unsupported sixframe_index format"):
            fn(bad, table=11)


def test_empty_genome_index_keeps_its_table(tmp_path: Path) -> None:
    """Edge case: ``build_sixframe({}, {}, table=4)`` has no frames; saved and loaded with 4."""
    index = build_sixframe({}, {}, table=4)
    assert index.frames == () and index.translation_table == 4
    path = tmp_path / "empty.npz"
    save_index(path, index, {})
    got, decoys = sf.load_index(path, table=4)
    assert got.frames == () and got.translation_table == 4 and decoys == {}
    with pytest.raises(ValueError, match="rebuild the six-frame index"):
        sf.load_index(path, table=11)


# ---------------------------------------------------------------- table 11 byte identity


def test_table_11_frames_are_byte_identical_to_main() -> None:
    """Default and explicit ``table=11`` give the target frames main wrote (format-2 fixture)."""
    genome, reps, _ = world()
    explicit = build_sixframe(genome, reps, table=11)
    default = build_sixframe(genome, reps)
    with np.load(FORMAT2) as old:
        want, bounds = old["target_aa"].tobytes(), old["target_bounds"].tolist()
    for idx in (explicit, default):
        assert b"".join(f.aa.tobytes() for f in idx.frames) == want
        assert np.cumsum([0] + [len(f.aa) for f in idx.frames]).tolist() == bounds
        assert idx.translation_table == 11


def test_table_11_step_index_is_byte_identical_to_format_2(tmp_path: Path) -> None:
    """AC 1 for the index: every array of the step's table-11 work file equals main's format-2
    file, and the meta differs only in ``format`` (2 -> 3) and the new ``translation_table``."""
    out = run_sixframe_index(tmp_path)
    assert SEED == 7
    with np.load(FORMAT2) as old, np.load(out) as new:
        assert set(new.files) == set(old.files)
        for key in old.files:
            if key == "meta":
                continue
            a, b = old[key], new[key]
            assert a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes(), key
        m_old, m_new = json.loads(str(old["meta"])), json.loads(str(new["meta"]))
    assert (m_old["format"], m_new["format"]) == (2, 3)
    assert m_new["translation_table"] == 11
    assert set(m_new) == set(m_old) | {"translation_table"}
    assert {k: v for k, v in m_new.items() if k not in ("format", "translation_table")} == {
        k: v for k, v in m_old.items() if k != "format"
    }
    loaded, _ = sf.load_index(out, table=11)
    old_idx, _ = sf.load_index(FORMAT2, table=11)
    assert [f.aa.tobytes() for f in loaded.frames] == [f.aa.tobytes() for f in old_idx.frames]
    assert sf.load_stop_cells(out, table=11).sites == sf.load_stop_cells(FORMAT2, table=11).sites


def test_table_4_changes_only_the_tga_cells() -> None:
    """Under 4 the frames differ from 11 exactly where 11 has ``*`` over a TGA, and there hold W."""
    genome, reps, _ = world()
    i11, i4 = build_sixframe(genome, reps, table=11), build_sixframe(genome, reps, table=4)
    n_changed = 0
    for f11, f4 in zip(i11.frames, i4.frames, strict=True):
        diff = np.flatnonzero(f11.aa != f4.aa)
        assert np.all(f11.aa[diff] == STOP) and np.all(f4.aa[diff] == 17)
        n_changed += len(diff)
    assert n_changed > 0
