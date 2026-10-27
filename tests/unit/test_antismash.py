"""antiSMASH command builder and JSON parser (M2-18 examples 6-8; GIPPA_SPEC.md §7.1.1, §6.5).

Rules (spec M2-18): ``records[].id`` = contig; ``records[].areas[]`` ``start``/``end`` are 0-based,
end exclusive; ``start > end`` on a circular record → origin-crossing, stored unwrapped
(``end + length``, ``gippa.model.coords``); ``products`` in antiSMASH order with duplicates
collapsed; core genes = ``locus_tag`` of CDS features in the area whose ``gene_kind`` qualifier
contains ``biosynthetic``.

Ambiguities resolved here (test author):

* "``gene_kind`` contains ``biosynthetic``" means the qualifier **list** contains the exact value
  ``biosynthetic``; ``biosynthetic-additional`` (also in the recordings) is not a core gene.
* "CDS in the area": the CDS lies inside the area (for an origin-crossing area, either side of
  the origin); core genes are listed in genome (feature) order, without duplicates.
* CDS locations are antiSMASH/Biopython strings ``[a:b](±)``, possibly fuzzy (``[<a:b]``,
  ``[a:>b]``, seen in the Nissle recording) or ``join{…}`` across the origin.
* Topology is taken from ``replicons`` (the GIPPA replicon table); ``start > end`` on a linear
  replicon is corrupt input → ``ValueError``.
* A record id absent from ``replicons`` → ``ValueError`` naming it, even when it has no areas.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from gippa.model.coords import Interval, Replicon, check_interval, split_at_origin

REPS = {
    "chr1": Replicon("chr1", 100_000, "circular", "chromosome"),
    "lin": Replicon("lin", 50_000, "linear", "unknown"),
}


def _a() -> ModuleType:
    from gippa.annotate import antismash

    return antismash


def _cds(loc: str, tag: str, kinds: list[str] | None = None) -> dict[str, Any]:
    q: dict[str, Any] = {"locus_tag": [tag], "product": ["hypothetical protein"]}
    if kinds is not None:
        q["gene_kind"] = kinds
    return {"location": loc, "type": "CDS", "qualifiers": q}


def _record(
    rid: str,
    areas: list[dict[str, Any]],
    features: list[dict[str, Any]],
    topology: str = "circular",
) -> dict[str, Any]:
    return {
        "id": rid,
        "name": rid,
        "annotations": {"topology": topology},
        "features": features,
        "areas": areas,
    }


def _area(start: int, end: int, products: list[str]) -> dict[str, Any]:
    return {
        "start": start,
        "end": end,
        "products": products,
        "protoclusters": {},
        "candidates": [],
        "subregions": [],
    }


def _write(tmp: Path, records: list[dict[str, Any]]) -> Path:
    p = tmp / "genome.json"
    p.write_text(json.dumps({"version": "8.0.4", "records": records, "schema": 4}))
    return p


# ---------------------------------------------------------------- example 6: command
# (amended by M2-18 lead decision 4: ``databases`` is required)


def test_example_6_command() -> None:
    from gippa.annotate.antismash_cmd import antismash_command

    got = antismash_command(Path("bakta.gbff"), Path("out"), 8, Path("/db/antismash"))
    assert got == [
        "antismash",
        "--genefinding-tool",
        "none",
        "--minimal",
        "--databases",
        "/db/antismash",
        "--cpus",
        "8",
        "--output-dir",
        "out",
        "bakta.gbff",
    ]
    assert all(isinstance(a, str) for a in got)


@pytest.mark.parametrize("threads", [1, 17])
def test_command_threads(threads: int) -> None:
    from gippa.annotate.antismash_cmd import antismash_command

    got = antismash_command(Path("x/in.gbff"), Path("o"), threads, Path("d/antismash"))
    assert got[got.index("--cpus") + 1] == str(threads)
    assert got[-1] == "x/in.gbff"


# ---------------------------------------------------------------- example 7: one region


def test_example_7_region_with_core_gene(tmp_path: Path) -> None:
    a = _a()
    rec = _record(
        "chr1",
        [_area(1000, 25000, ["NRPS"])],
        [
            _cds("[500:900](+)", "ESC_1000", ["biosynthetic"]),  # outside the area
            _cds("[1200:2400](+)", "ESC_1233"),  # no gene_kind
            _cds("[3000:9000](+)", "ESC_1234", ["biosynthetic"]),
            _cds("[9100:9900](-)", "ESC_1235", ["biosynthetic-additional"]),
            _cds("[26000:27000](+)", "ESC_1240", ["biosynthetic"]),  # outside the area
        ],
    )
    got = a.parse_antismash_json(_write(tmp_path, [rec]), REPS)
    assert got == [
        a.BgcRegion("chr1", Interval("chr1", 1000, 25000, "."), ("NRPS",), ("ESC_1234",))
    ]


# ---------------------------------------------------------------- example 8: empty


def test_example_8_no_records(tmp_path: Path) -> None:
    a = _a()
    assert a.parse_antismash_json(_write(tmp_path, []), REPS) == []


def test_records_without_areas_give_no_regions(tmp_path: Path) -> None:
    a = _a()
    recs = [
        _record("chr1", [], [_cds("[0:300](+)", "T_1", ["biosynthetic"])]),
        _record("lin", [], []),
    ]
    assert a.parse_antismash_json(_write(tmp_path, recs), REPS) == []


# ---------------------------------------------------------------- edge cases


def test_unknown_record_id_raises_naming_it(tmp_path: Path) -> None:
    a = _a()
    rec = _record("chrZ", [_area(0, 100, ["NRPS"])], [])
    with pytest.raises(ValueError, match="chrZ"):
        a.parse_antismash_json(_write(tmp_path, [rec]), REPS)


def test_multiple_products_kept_in_order_duplicates_collapsed(tmp_path: Path) -> None:
    a = _a()
    rec = _record(
        "chr1", [_area(10, 5000, ["T1PKS", "NRPS", "T1PKS", "NRP-metallophore", "NRPS"])], []
    )
    (r,) = a.parse_antismash_json(_write(tmp_path, [rec]), REPS)
    assert r.products == ("T1PKS", "NRPS", "NRP-metallophore")
    assert r.core_genes == ()


def test_origin_crossing_area_on_circular_record(tmp_path: Path) -> None:
    """``start > end`` on a circular record: ``Interval(start, end + length)``; core genes on both
    sides of the origin (incl. a ``join`` CDS spanning it) are found, in genome order."""
    a = _a()
    rec = _record(
        "chr1",
        [_area(95_000, 4_000, ["terpene"])],
        [
            _cds("[1000:2000](+)", "T_0002", ["biosynthetic"]),
            _cds("[5000:6000](+)", "T_0005", ["biosynthetic"]),  # outside
            _cds("[90000:94000](-)", "T_0090", ["biosynthetic"]),  # outside
            _cds("[96000:97000](-)", "T_0096", ["biosynthetic"]),
            _cds("join{[99800:100000](+), [0:400](+)}", "T_0099", ["biosynthetic"]),
        ],
    )
    (r,) = a.parse_antismash_json(_write(tmp_path, [rec]), REPS)
    assert r.interval == Interval("chr1", 95_000, 104_000, ".")
    assert r.interval.crosses_origin(REPS["chr1"])
    assert set(r.core_genes) == {"T_0096", "T_0099", "T_0002"}
    assert len(r.core_genes) == 3


def test_start_after_end_on_linear_record_raises(tmp_path: Path) -> None:
    a = _a()
    rec = _record("lin", [_area(40_000, 1_000, ["NRPS"])], [], topology="linear")
    with pytest.raises(ValueError):
        a.parse_antismash_json(_write(tmp_path, [rec]), REPS)


def test_fuzzy_cds_locations_are_parsed(tmp_path: Path) -> None:
    """Partial CDS at contig/feature edges (``[<a:b]``, ``[a:>b]``, as in the Nissle run)."""
    a = _a()
    rec = _record(
        "lin",
        [_area(0, 9_000, ["arylpolyene"])],
        [
            _cds("[<0:300](+)", "L_0001", ["biosynthetic"]),
            _cds("[1000:>1500](-)", "L_0002", ["biosynthetic"]),
            _cds("[2000:2600](+)", "L_0003", ["biosynthetic-additional"]),
        ],
        topology="linear",
    )
    (r,) = a.parse_antismash_json(_write(tmp_path, [rec]), REPS)
    assert r.interval == Interval("lin", 0, 9_000, ".")
    assert r.core_genes == ("L_0001", "L_0002")


def test_area_ending_at_contig_end(tmp_path: Path) -> None:
    a = _a()
    rec = _record(
        "lin",
        [_area(45_000, 50_000, ["RiPP-like"])],
        [_cds("[49000:50000](+)", "L_0049", ["biosynthetic"])],
        topology="linear",
    )
    (r,) = a.parse_antismash_json(_write(tmp_path, [rec]), REPS)
    assert r.interval == Interval("lin", 45_000, 50_000, ".")
    assert r.core_genes == ("L_0049",)


def test_regions_across_records_keep_record_and_area_order(tmp_path: Path) -> None:
    a = _a()
    recs = [
        _record("chr1", [_area(100, 200, ["A"]), _area(5000, 6000, ["B"])], []),
        _record("lin", [_area(10, 20, ["C"])], [], topology="linear"),
    ]
    got = a.parse_antismash_json(_write(tmp_path, recs), REPS)
    assert [(r.contig, r.products) for r in got] == [
        ("chr1", ("A",)),
        ("chr1", ("B",)),
        ("lin", ("C",)),
    ]


def test_core_gene_lies_in_one_of_two_areas(tmp_path: Path) -> None:
    a = _a()
    rec = _record(
        "chr1",
        [_area(1000, 3000, ["A"]), _area(5000, 8000, ["B"])],
        [
            _cds("[1500:2000](+)", "X_1", ["biosynthetic"]),
            _cds("[6000:7000](-)", "X_2", ["biosynthetic"]),
        ],
    )
    r1, r2 = a.parse_antismash_json(_write(tmp_path, [rec]), REPS)
    assert r1.core_genes == ("X_1",) and r2.core_genes == ("X_2",)


def test_bgc_region_is_frozen_slotted_with_tuples(tmp_path: Path) -> None:
    a = _a()
    rec = _record(
        "chr1", [_area(1000, 25000, ["NRPS"])], [_cds("[3000:9000](+)", "E_1", ["biosynthetic"])]
    )
    path = _write(tmp_path, [rec])
    (r,) = a.parse_antismash_json(path, REPS)
    assert dataclasses.is_dataclass(r) and not hasattr(r, "__dict__")
    with pytest.raises(dataclasses.FrozenInstanceError):
        r.contig = "x"  # type: ignore[misc]
    assert isinstance(r.products, tuple) and isinstance(r.core_genes, tuple)
    assert a.parse_antismash_json(path, REPS) == a.parse_antismash_json(path, REPS)


# ---------------------------------------------------------------- Hypothesis properties


@st.composite
def _circular_area(draw: st.DrawFn) -> tuple[int, int, int]:
    length = draw(st.integers(100, 1_000_000))
    start = draw(st.integers(0, length - 1))
    end = draw(st.integers(1, length))
    return length, start, end


@settings(max_examples=80, deadline=None)
@given(_circular_area())
def test_property_area_coordinates(
    tmp_path_factory: pytest.TempPathFactory, case: tuple[int, int, int]
) -> None:
    """0-based half-open in = out; ``start > end`` unwraps by ``+length``; always a valid interval
    on the replicon whose origin split covers exactly the area's bases."""
    a = _a()
    length, start, end = case
    assume(start != end)
    reps = {"c": Replicon("c", length, "circular")}
    rec = _record("c", [_area(start, end, ["NRPS"])], [])
    (r,) = a.parse_antismash_json(_write(tmp_path_factory.mktemp("as"), [rec]), reps)
    want_end = end if start < end else end + length
    assert r.interval == Interval("c", start, want_end, ".")
    check_interval(r.interval, reps["c"])
    parts = split_at_origin(r.interval, reps["c"])
    covered = sum(p.length_bp for p in parts)
    assert covered == (end - start if start < end else length - start + end)


@settings(max_examples=60, deadline=None)
@given(
    st.lists(
        st.sampled_from(["NRPS", "T1PKS", "terpene", "RiPP-like", "NRP-metallophore"]),
        min_size=1,
        max_size=8,
    ),
    st.lists(
        st.tuples(
            st.integers(0, 9_000),
            st.integers(1, 900),
            st.sampled_from(
                [
                    ["biosynthetic"],
                    ["biosynthetic-additional"],
                    ["transport"],
                    ["regulatory"],
                    None,
                    ["biosynthetic", "transport"],
                ]
            ),
        ),
        max_size=12,
    ),
)
def test_property_products_and_core_genes(
    tmp_path_factory: pytest.TempPathFactory,
    products: list[str],
    cds: list[tuple[int, int, list[str] | None]],
) -> None:
    a = _a()
    area = (2_000, 7_000)
    # whole CDS inside or outside the area (antiSMASH areas end on CDS boundaries)
    whole = [
        (s, n, k)
        for s, n, k in cds
        if s + n <= area[0] or s >= area[1] or (area[0] <= s and s + n <= area[1])
    ]
    feats = sorted(
        (_cds(f"[{s}:{s + n}](+)", f"G_{i:04d}", kinds) for i, (s, n, kinds) in enumerate(whole)),
        key=lambda f: int(f["location"][1:].split(":")[0]),
    )
    rec = _record("lin", [_area(*area, products)], feats, topology="linear")
    reps = {"lin": Replicon("lin", 10_000, "linear")}
    (r,) = a.parse_antismash_json(_write(tmp_path_factory.mktemp("as"), [rec]), reps)
    assert r.products == tuple(dict.fromkeys(products))
    want = [
        f["qualifiers"]["locus_tag"][0]
        for f in feats
        if "biosynthetic" in (f["qualifiers"].get("gene_kind") or [])
        and area[0] <= int(f["location"][1:].split(":")[0])
        and int(f["location"].split(":")[1].split("]")[0]) <= area[1]
    ]
    assert r.core_genes == tuple(want)
