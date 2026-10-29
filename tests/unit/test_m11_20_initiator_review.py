"""M11-20 review findings on PR #490 (docs/specs/M11-20-initiator-methionine.md § Definition rules
1-3 and 7, § Examples 6-7).

- N1: ``match_exact`` checks the start cells of every decoy kind the queries use whenever
  ``starts`` has sites, independent of the data (also when no query starts with M); ``save_index``
  refuses start cells without a decoy kind's cells.
- N5: ``annotation.start_codons`` is validated in the config: each entry 3 letters of ACGT (any
  case, rule 2 upper-cases) and not a stop under ``annotation.translation_table`` (TGA is a stop
  at 11 only; TAA and TAG at 4, 11 and 25).
- m04: no credit when ``3k`` exceeds the replicon length (rule 3), even where the window fits the
  extended circular frame.
- m17: the ``mapping.start_sites`` event's ``usable_k7`` equals the independent count of
  ``helpers.initiator_toys.usable``.
- A minus-strand gene whose first codon lies past the origin (``end - 3 >= L``) is credited.

Readings taken (report): N1's raise is the existing ``ValueError`` text ``start cells have no
<kind> decoy cells``; a kind the queries do not use may lack cells; ``starts`` without sites needs
no cells. N5's error is a ``ValueError`` (pydantic ``ValidationError`` included, as every config
test) naming ``start_codons`` and the offending codon.
"""

from __future__ import annotations

import importlib
import logging
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.initiator_toys import (
    G3_NT,
    STARTS,
    cds,
    credited_rows,
    mirrored,
    ref_cells,
    row_set,
    start_keys,
    toy_c1,
    toy_c5,
    usable,
    with_decoys,
)
from helpers.readthrough_toys import model, rc, write_annotation, write_genome

from gippa.mapping.exact import match_exact
from gippa.mapping.sixframe import save_index
from gippa.model import Replicon


def ini() -> ModuleType:
    return importlib.import_module("gippa.mapping.initiator")


def starts_of(world: Any, table: int = 11) -> tuple[Any, Any, Any]:
    genome, reps, gm = world
    m = ini()
    index, decoys = with_decoys(genome, reps, table)
    sites, _ = m.annotated_starts(gm, genome, STARTS, table=table)
    target = m.start_frame_cells(index, sites)
    dc = {k: m.mirror_start_cells(k, fr, target) for k, fr in decoys.items()}
    return index, decoys, m.StartCells(sites, target, dc)


def model_rows(world: Any, queries: list[str], kind: str) -> set[Any]:
    """The independent rule-3 rows (helpers.initiator_toys) plus nothing else."""
    genome, reps, gm = world
    index, decoys = with_decoys(genome, reps)
    tc = ref_cells(index, start_keys(gm, genome))
    cells = {"target": tc, **{k: mirrored(k, decoys[k], tc) for k in decoys}}
    return credited_rows(queries, index, decoys, dict.fromkeys(queries, kind), cells)


# ---------------------------------------------------------------- N1: decoy cells, data-free


@pytest.mark.parametrize("missing", ["reverse", "shuffle"])
def test_match_exact_raises_on_missing_wanted_decoy_cells_without_m_queries(missing: str) -> None:
    index, decoys, cells = starts_of(toy_c1())
    assert cells.sites
    bad = ini().StartCells(
        cells.sites, cells.target, {k: v for k, v in cells.decoys.items() if k != missing}
    )
    queries = ["PPVAELWHK"]  # V-initial: no query starts with M
    with pytest.raises(ValueError, match=f"start cells have no {missing} decoy cells"):
        match_exact(queries, index, decoys, dict.fromkeys(queries, missing), starts=bad)


def test_match_exact_raises_on_missing_decoy_cells_with_m_queries() -> None:
    index, decoys, cells = starts_of(toy_c1())
    bad = ini().StartCells(cells.sites, cells.target, {"shuffle": cells.decoys["shuffle"]})
    with pytest.raises(ValueError, match="start cells have no reverse decoy cells"):
        match_exact(["MAELWHK"], index, decoys, {"MAELWHK": "reverse"}, starts=bad)


def test_match_exact_accepts_missing_cells_of_an_unused_kind() -> None:
    index, decoys, cells = starts_of(toy_c1())
    part = ini().StartCells(cells.sites, cells.target, {"reverse": cells.decoys["reverse"]})
    queries = ["MAELWHK", "PPVAELWHK"]
    got = match_exact(queries, index, decoys, dict.fromkeys(queries, "reverse"), starts=part)
    full = match_exact(queries, index, decoys, dict.fromkeys(queries, "reverse"), starts=cells)
    assert row_set(got) == row_set(full)


def test_match_exact_without_sites_needs_no_decoy_cells() -> None:
    index, decoys, cells = starts_of(toy_c1())
    empty = ini().StartCells((), cells.target, {})
    queries = ["MAELWHK", "PPVAELWHK"]
    got = match_exact(queries, index, decoys, dict.fromkeys(queries, "reverse"), starts=empty)
    old = match_exact(queries, index, decoys, dict.fromkeys(queries, "reverse"))
    assert row_set(got) == row_set(old)


@pytest.mark.parametrize("missing", ["reverse", "shuffle"])
def test_save_index_raises_on_start_cells_without_a_decoy_kind(
    tmp_path: Path, missing: str
) -> None:
    index, decoys, cells = starts_of(toy_c1())
    bad = ini().StartCells(
        cells.sites, cells.target, {k: v for k, v in cells.decoys.items() if k != missing}
    )
    with pytest.raises(ValueError, match=f"start cells have no {missing} decoy cells"):
        save_index(tmp_path / "six.npz", index, decoys, None, starts=bad)


# ---------------------------------------------------------------- N5: config validation


def _cfg(*overrides: str) -> Any:
    from gippa.config import resolve_config

    return resolve_config(overrides=list(overrides))


@pytest.mark.parametrize(
    ("table", "codon"),
    [(11, "TGA"), (11, "TAA"), (11, "TAG"), (4, "TAA"), (4, "TAG"), (25, "TAA"), (25, "TAG")],
)
def test_start_codon_that_is_a_stop_under_the_table_is_rejected(table: int, codon: str) -> None:
    with pytest.raises(ValueError, match=rf"(?s)start_codons.*{codon}|{codon}.*start_codons"):
        _cfg(f"annotation.translation_table={table}", f"annotation.start_codons=[ATG, {codon}]")


def test_lower_case_stop_start_codon_is_rejected() -> None:
    with pytest.raises(ValueError, match="start_codons"):
        _cfg("annotation.start_codons=[ATG, tga]")


def test_stop_message_names_the_table() -> None:
    with pytest.raises(ValueError, match="stop") as err:
        _cfg("annotation.start_codons=[ATG, TGA]")
    assert "11" in str(err.value)


@pytest.mark.parametrize("table", [4, 25])
def test_tga_start_codon_is_allowed_where_tga_is_not_a_stop(table: int) -> None:
    cfg = _cfg(f"annotation.translation_table={table}", "annotation.start_codons=[ATG, TGA]")
    assert list(cfg.annotation.start_codons) == ["ATG", "TGA"]


@pytest.mark.parametrize("codon", ["ANG", "AT", "ATGA", "AUG", "", "A-G"])
def test_start_codon_not_three_acgt_letters_is_rejected(codon: str) -> None:
    with pytest.raises(ValueError, match="start_codons"):
        _cfg(f"annotation.start_codons=[ATG, '{codon}']")


def test_default_and_lower_case_start_codons_are_valid() -> None:
    assert list(_cfg().annotation.start_codons) == ["ATG", "GTG", "TTG"]
    got = _cfg("annotation.start_codons=[atg, gtg]").annotation.start_codons
    assert [c.upper() for c in got] == ["ATG", "GTG"]


# ---------------------------------------------------------------- m04: 3k > L


def toy_short_circle():  # type: ignore[no-untyped-def]
    """Circular 21 nt, all one GTG gene ``+ [0, 21)`` (aa ``MAELWHK``, no stop): frame 0 of
    the 42-nt extended frame is ``VAELWHKVAELWHK``, so an 8-aa window fits it but 3k = 24 > 21."""
    nt = "GTGGCTGAACTGTGGCATAAA"
    reps = {"s": Replicon("s", 21, "circular")}
    return {"s": nt}, reps, model(reps, [cds("GS", "s", 0, 21, "+", "MAELWHK", "GTG")])


@pytest.mark.parametrize("kind", ["reverse", "shuffle"])
def test_m04_no_credit_when_3k_exceeds_the_replicon(kind: str) -> None:
    world = toy_short_circle()
    index, decoys, cells = starts_of(world)
    assert len(cells.sites) == 1
    longq = ["MAELWHKV"]  # k = 8, 3k = 24 > L = 21
    got = match_exact(longq, index, decoys, dict.fromkeys(longq, kind), starts=cells)
    old = match_exact(longq, index, decoys, dict.fromkeys(longq, kind))
    assert row_set(got) == row_set(old)
    assert got[got["search"] == "target"].empty
    assert not model_rows(world, longq, kind)


def test_m04_credit_when_3k_equals_the_replicon() -> None:
    world = toy_short_circle()
    index, decoys, cells = starts_of(world)
    q = ["MAELWHK"]  # k = 7, 3k = 21 = L: credited
    got = match_exact(q, index, decoys, {"MAELWHK": "reverse"}, starts=cells)
    tgt = got[got["search"] == "target"]
    assert list(zip(tgt["start"], tgt["end"], tgt["strand"], strict=True)) == [(0, 21, "+")]
    assert row_set(got) == model_rows(world, q, "reverse")


# ---------------------------------------------------------------- minus strand past the origin


def toy_c6():  # type: ignore[no-untyped-def]
    """Circular 90 nt, G6 ``- [69, 93)`` reading ``TTG GCT GAA CTG TGG CAT AAA TAA`` (aa
    ``MAELWHK``): its first codon sits at plus bases 90..92, i.e. 0..2 (``end - 3 = 90 >= L``)."""
    seq = list("C" * 90)
    for i, b in enumerate(rc(G3_NT)):
        seq[(69 + i) % 90] = b
    assert rc("".join(seq[0:3])) == "TTG"
    reps = {"c6": Replicon("c6", 90, "circular")}
    return (
        {"c6": "".join(seq)},
        reps,
        model(reps, [cds("G6", "c6", 69, 93, "-", "MAELWHK", "TTG")]),
    )


def test_minus_strand_first_codon_past_the_origin_site() -> None:
    genome, _reps, gm = toy_c6()
    sites, skipped = ini().annotated_starts(gm, genome, STARTS, table=11)
    assert [(s.contig, s.strand, s.codon_start, tuple(s.locus_tags), s.codon) for s in sites] == [
        ("c6", "-", 0, ("G6",), "TTG")
    ]
    assert sum(skipped.values()) == 0


def test_minus_strand_first_codon_past_the_origin_is_credited() -> None:
    world = toy_c6()
    index, decoys, cells = starts_of(world)
    df = match_exact(["MAELWHK"], index, decoys, {"MAELWHK": "reverse"}, starts=cells)
    tgt = df[df["search"] == "target"]
    assert list(zip(tgt["start"], tgt["end"], tgt["strand"], strict=True)) == [(72, 93, "-")]
    assert row_set(df) == model_rows(world, ["MAELWHK"], "reverse")
    assert match_exact(["MAELWHK"], index, decoys, {"MAELWHK": "reverse"}).empty


# ---------------------------------------------------------------- m17: usable_k7 of the event


def _step_event(tmp: Path, world: Any, caplog: pytest.LogCaptureFixture) -> dict[str, Any]:
    from gippa.config import resolve_config
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    genome, _reps, gm = world
    ann = write_annotation(gm, tmp / "ann")
    pq = tmp / "peptides.parquet"
    pd.DataFrame(
        {
            "peptide": ["MAELWHK", "SAMPLEWLDK"],
            "source": ["database", "denovo"],  # both decoy kinds are built
            "q_value": [0.001, np.nan],
            "pep": 0.01,
            "score": 0.9,
            "is_decoy": False,
            "search_db": "other",
        }
    ).to_parquet(pq, index=False)
    ins = {
        "genome": write_genome(genome, tmp / "genome.fa"),
        "replicons": ann["replicons"],
        "genes": ann["genes"],
        "peptides": pq,
    }
    load_step_modules()
    caplog.set_level(logging.INFO)
    reg.run_step(
        "sixframe_index",
        resolve_config(overrides=["seed=7"]),
        {k: ins[k] for k in reg.get("sixframe_index").inputs},
        {"sixframe": tmp / "sixframe.npz"},
    )
    ev = [r for r in caplog.records if getattr(r, "event", None) == "mapping.start_sites"]
    assert len(ev) == 1
    return ev[0].data  # type: ignore[attr-defined,no-any-return]


def _model_usable(world: Any, k: int = 7) -> dict[str, int]:
    genome, reps, gm = world
    index, decoys = with_decoys(genome, reps)
    tc = ref_cells(index, start_keys(gm, genome))
    out = {
        "target": sum(
            usable("target", f, c, k, reps[f.contig].length)
            for f, c in zip(index.frames, tc, strict=True)
        )
    }
    for kind, frames in decoys.items():
        dc = mirrored(kind, frames, tc)
        out[kind] = sum(
            usable(kind, f, c, k, reps[f.contig].length) for f, c in zip(frames, dc, strict=True)
        )
    return out


@pytest.mark.parametrize(
    ("name", "world", "target"),
    [
        ("c1", toy_c1, 2),  # the site and its origin copy in the extended frame, both stop-free
        ("c5", lambda: toy_c5(4), 0),  # at table 11 the TGA in codon 4 is a stop: no window
        ("short", toy_short_circle, 2),  # 3k = 21 = L: cells 0 and 7 of the 14-aa frame fit
    ],
)
def test_m17_usable_k7_equals_the_toy_model(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, name: str, world: Any, target: int
) -> None:
    w = world()
    data = _step_event(tmp_path, w, caplog)
    want = _model_usable(w)
    assert want["target"] == target, name  # the model against a hand count
    assert data["usable_k7"] == want, name
