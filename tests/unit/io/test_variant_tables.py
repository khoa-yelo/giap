"""M7-03 table writers and schemas: ``pileup.tsv.gz`` and ``variants.tsv``
(docs/specs/M7-03-variant-tables-and-vcf.md rule 10; GIAP_SPEC.md §9.3, §9.4, §9.5, §6.5;
ADR-0006).

Interface reading (test author; see the report to the lead):

* ``write_pileup(path, pile)`` takes whatever ``pileup`` returned; ``write_variants(path, table)``
  takes the table returned by ``vcf_records``; ``read_variants`` returns the §9.3 columns with
  0-based ``codon_start`` (the inverse of the writer) and list cells for ``flags``.
* ``class_error_rate`` is written ``unestimated`` and ``class_error_n`` ``NA`` until M7-10.
* The Frictionless descriptors are ``src/giap/schemas/tables/{pileup,variants}.schema.json``
  and ``src/giap/schemas/datapackage.json`` lists both resources (§9.5).
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

import pytest
from helpers.variant_tables_fixture import (
    GM_,
    GP,
    PILEUP_COLUMNS,
    VARIANT_FLAGS,
    VARIANTS_COLUMNS,
    aa_at,
    al,
    flag_set,
    match_op,
    peptides_at,
    raw_hit,
    records,
    ref_peptides,
    run,
    sub,
    two,
    vt,
)

ROOT = Path(__file__).resolve().parents[3]
SCHEMAS = ROOT / "src" / "giap" / "schemas"
TABLES = SCHEMAS / "tables"


def _schema(name: str) -> dict[str, Any]:
    data = json.loads((TABLES / f"{name}.schema.json").read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def _fields(name: str) -> dict[str, dict[str, Any]]:
    return {f["name"]: f for f in _schema(name)["fields"]}


def _enum(name: str, col: str) -> set[str]:
    return set(_fields(name)[col]["constraints"]["enum"])


def _validate(path: Path, name: str) -> Any:
    frictionless = pytest.importorskip("frictionless")
    d = {k: v for k, v in _schema(name).items() if k != "coordinates"}
    res = frictionless.Resource(
        path=path.name, basepath=str(path.parent), format="tsv",
        schema=frictionless.Schema.from_descriptor(d),
    )  # fmt: skip
    return res.validate()


def _assert_valid(path: Path, name: str) -> None:
    report = _validate(path, name)
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def _c2_sub(codon: int) -> list[Any]:
    ref = aa_at("c2", codon, "+")
    alt = "W" if ref != "W" else "C"
    ops = [match_op("c2", codon - 3, "+"), ("substitution", (codon,), ref, alt),
           match_op("c2", codon + 3, "+")]  # fmt: skip
    return [raw_hit("c2", "+", ops)]


def _frameshift() -> list[Any]:
    ev = al().FrameshiftEvent(4, GP.codon(91), GP.codon(92) + 1, 1)
    ops = [match_op("chr", GP.codon(i), "+") for i in range(88, 92)]
    return [raw_hit("chr", "+", ops, frameshift=ev)]


def _bag() -> list[Any]:
    return (
        peptides_at(GP, 44, sub(GP, 44, "S"), 3)
        + peptides_at(GP, 50, sub(GP, 50, "D"), 2)
        + ref_peptides(GP, 50, 1)
        + peptides_at(GP, 60, two(GP, 60, "swap", "DS"), 1)
        + peptides_at(GP, 70, two(GP, 70, "block", "AS"), 2)
        + peptides_at(GP, 6, ("block", (GP.codon(6),), aa_at("chr", GP.codon(6), "+"), "GG"), 1)
        + peptides_at(GM_, 10, sub(GM_, 10, "S"), 1)
        + _frameshift()
        + _c2_sub(999)
    )


def _lines(path: Path) -> list[list[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as fh:  # type: ignore[operator]
        return [ln.rstrip("\n").split("\t") for ln in fh]


# --------------------------------------------------------------------------- schemas


def test_pileup_schema_fields_and_enums() -> None:
    s = _schema("pileup")
    assert tuple(f["name"] for f in s["fields"]) == PILEUP_COLUMNS
    assert s["coordinates"] == "1-based-closed"
    assert "NA" in s["missingValues"]
    assert _enum("pileup", "call") == {"REF", "ALT", "mixed", "no_call"}
    assert _enum("pileup", "strand") == {"+", "-"}
    f = _fields("pileup")
    for col in ("codon_start", "n_ref", "n_il"):
        assert f[col]["type"] == "integer"
    assert f["codon_start"]["constraints"]["minimum"] == 1


def test_variants_schema_fields_and_enums() -> None:
    s = _schema("variants")
    assert tuple(f["name"] for f in s["fields"]) == VARIANTS_COLUMNS
    assert s["coordinates"] == "1-based-closed"
    assert "NA" in s["missingValues"]
    assert {"ALT", "mixed"} <= _enum("variants", "call") <= {"REF", "ALT", "mixed", "no_call"}
    assert _enum("variants", "confidence") == {"high", "candidate"}
    assert _enum("variants", "alt_type") == {"substitution", "swap", "block", "frameshift"}
    f = _fields("variants")
    for col in ("codon_start", "n_ref", "n_alt", "nt_changes"):
        assert f[col]["type"] == "integer"
    assert f["codon_start"]["constraints"]["minimum"] == 1


@pytest.mark.parametrize("name", ["pileup", "variants"])
def test_schemas_are_valid_frictionless(name: str) -> None:
    frictionless = pytest.importorskip("frictionless")
    d = {k: v for k, v in _schema(name).items() if k != "coordinates"}
    frictionless.Schema.from_descriptor(d)


def test_datapackage_lists_both_resources() -> None:
    pkg = json.loads((SCHEMAS / "datapackage.json").read_text(encoding="utf-8"))
    paths = {r["name"]: r["path"] for r in pkg["resources"]}
    assert paths["pileup"].endswith("tables/pileup.tsv.gz")
    assert paths["variants"].endswith("tables/variants.tsv")


# --------------------------------------------------------------------------- pileup.tsv.gz


def test_write_pileup_header_and_b7_row(tmp_path: Path) -> None:
    pile, _ = run(peptides_at(GP, 44, sub(GP, 44, "S"), 3))
    path = tmp_path / "pileup.tsv.gz"
    vt().write_pileup(path, pile)
    with path.open("rb") as fh:
        assert fh.read(2) == b"\x1f\x8b"
    lines = _lines(path)
    assert tuple(lines[0]) == PILEUP_COLUMNS
    rows = {(r[0], r[1], r[2]): r for r in lines[1:]}
    assert rows[("chr", "1001", "+")] == [
        "chr", "1001", "+", "NISS_00012", "45", "A", "0", "0", "S:3", "ALT",
    ]  # fmt: skip
    assert all(not c.startswith("_") for c in lines[0])
    assert len(lines) - 1 == len(pile)


def test_write_pileup_rows_in_section_6_5_order(tmp_path: Path) -> None:
    pile, _ = run(_bag())
    path = tmp_path / "pileup.tsv.gz"
    vt().write_pileup(path, pile)
    body = _lines(path)[1:]
    order = {"chr": 0, "c2": 1}
    keys = [(order[r[0]], int(r[1]), r[2] == "-") for r in body]
    assert keys == sorted(keys)
    assert all(1 <= int(r[1]) <= (2400 if r[0] == "chr" else 1000) for r in body)
    assert {"NA"} <= {r[3] for r in body if r[0] == "c2"}  # off-gene: gene NA
    _assert_valid(path, "pileup")


def test_pileup_origin_codon_is_one_based_within_length(tmp_path: Path) -> None:
    pile, _ = run(_c2_sub(999))
    path = tmp_path / "pileup.tsv.gz"
    vt().write_pileup(path, pile)
    starts = sorted(int(r[1]) for r in _lines(path)[1:])
    assert starts == [3, 997, 1000]  # codons 996, 999 and 1002 (wrapped to 2)


def test_pileup_gz_is_byte_identical_across_runs(tmp_path: Path) -> None:
    """TE7: two independent runs give identical bytes (gzip ``mtime = 0``, no file name)."""
    a, b = tmp_path / "a" / "pileup.tsv.gz", tmp_path / "b" / "pileup.tsv.gz"
    a.parent.mkdir()
    b.parent.mkdir()
    vt().write_pileup(a, run(_bag())[0])
    vt().write_pileup(b, run(_bag())[0])
    raw = a.read_bytes()
    assert raw == b.read_bytes()
    assert raw[4:8] == b"\x00\x00\x00\x00"  # MTIME
    assert not raw[3] & 0x08  # FLG.FNAME unset


def test_empty_pileup_is_header_only_gzip(tmp_path: Path) -> None:
    pile, _ = run([])
    path = tmp_path / "pileup.tsv.gz"
    vt().write_pileup(path, pile)
    assert _lines(path) == [list(PILEUP_COLUMNS)]
    _assert_valid(path, "pileup")


# --------------------------------------------------------------------------- variants.tsv


def test_write_variants_columns_exactly_section_9_3(tmp_path: Path) -> None:
    _, table, _, _ = records(_bag())
    path = tmp_path / "variants.tsv"
    vt().write_variants(path, table)
    lines = _lines(path)
    assert tuple(lines[0]) == VARIANTS_COLUMNS
    assert len(lines) - 1 == len(table)
    _assert_valid(path, "variants")


def test_write_variants_b7_row(tmp_path: Path) -> None:
    _, table, _, _ = records(
        peptides_at(GP, 44, sub(GP, 44, "S"), 3), {"flags": [("mass_equivalent", "snv")]}
    )
    path = tmp_path / "variants.tsv"
    vt().write_variants(path, table)
    (row,) = _lines(path)[1:]
    assert row == [
        "VAR_1", "chr", "1001", "+", "NISS_00012", "45", "GCT", "A", "S", "substitution",
        "ALT", "high", "0", "3", "1", "TCT", "unestimated", "NA", "denovo",
        "mass_equivalent;snv",
    ]  # fmt: skip


def test_write_variants_minus_and_origin_rows(tmp_path: Path) -> None:
    _, table, _, _ = records(_bag())
    path = tmp_path / "variants.tsv"
    vt().write_variants(path, table)
    rows = [dict(zip(VARIANTS_COLUMNS, r, strict=True)) for r in _lines(path)[1:]]
    (m,) = [r for r in rows if r["contig"] == "chr" and r["codon_start"] == "1568"]
    assert (m["strand"], m["ref_codon"], m["alt_codon_vcf"], m["call"]) == ("-", "GCT", "AGA",
                                                                             "ALT")  # fmt: skip
    (o,) = [r for r in rows if r["contig"] == "c2"]
    assert o["codon_start"] == "1000" and o["alt_codon_vcf"] == "NA"
    assert "vcf_unrepresentable" in o["flags"].split(";")
    (fs,) = [r for r in rows if r["alt_type"] == "frameshift"]
    assert fs["alt_codon_vcf"] == "NA"
    (b12,) = [r for r in rows if r["alt_type"] == "block" and r["alt_aa"] == "GG"]
    assert b12["alt_codon_vcf"] == "NA"
    (sw,) = [r for r in rows if r["alt_type"] == "swap"]
    assert (sw["alt_codon_vcf"], sw["nt_changes"]) == ("GATAGT", "4")
    assert {r["class_error_rate"] for r in rows} == {"unestimated"}
    assert {r["class_error_n"] for r in rows} == {"NA"}


def test_variants_rows_in_section_6_5_order(tmp_path: Path) -> None:
    _, table, _, _ = records(_bag())
    path = tmp_path / "variants.tsv"
    vt().write_variants(path, table)
    body = _lines(path)[1:]
    order = {"chr": 0, "c2": 1}
    keys = [(order[r[1]], int(r[2]), r[3] == "-") for r in body]
    assert keys == sorted(keys)


def test_no_sites_header_only_variants(tmp_path: Path) -> None:
    _, table, _, _ = records(ref_peptides(GP, 44, 2))
    path = tmp_path / "variants.tsv"
    vt().write_variants(path, table)
    assert _lines(path) == [list(VARIANTS_COLUMNS)]
    _assert_valid(path, "variants")


def test_read_variants_round_trip(tmp_path: Path) -> None:
    _, table, _, _ = records(_bag(), None)
    path = tmp_path / "variants.tsv"
    m = vt()
    m.write_variants(path, table)
    back = m.read_variants(path)
    assert tuple(back.columns) == VARIANTS_COLUMNS
    assert back["variant_id"].tolist() == table["variant_id"].tolist()
    assert [int(x) for x in back["codon_start"]] == [int(x) for x in table["codon_start"]]
    for col in ("n_ref", "n_alt"):
        assert [int(x) for x in back[col]] == [int(x) for x in table[col]]
    assert [flag_set(x) for x in back["flags"]] == [flag_set(x) for x in table["flags"]]
    m.write_variants(tmp_path / "again.tsv", back)
    assert (tmp_path / "again.tsv").read_bytes() == path.read_bytes()


# --------------------------------------------------------------------------- flag vocabulary


def _swap_flags(path: Path, value: str) -> Path:
    lines = path.read_text(encoding="utf-8").splitlines()
    head = lines[0].split("\t")
    i = head.index("flags")
    cells = lines[1].split("\t")
    cells[i] = value
    out = path.with_name("edited.tsv")
    out.write_text("\n".join([lines[0], "\t".join(cells), *lines[2:]]) + "\n", encoding="utf-8")
    return out


@pytest.mark.parametrize("flags", [(f,) for f in VARIANT_FLAGS] + [tuple(VARIANT_FLAGS)])
def test_variants_schema_accepts_the_variants_flag_vocabulary(
    tmp_path: Path, flags: tuple[str, ...]
) -> None:
    _, table, _, _ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 1))
    path = tmp_path / "variants.tsv"
    vt().write_variants(path, table)
    _assert_valid(_swap_flags(path, ";".join(flags)), "variants")


@pytest.mark.parametrize("bad", ["not_a_flag", "snv;fs", "crosses_origin"])
def test_variants_schema_rejects_other_flags(tmp_path: Path, bad: str) -> None:
    _, table, _, _ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 1))
    path = tmp_path / "variants.tsv"
    vt().write_variants(path, table)
    assert not _validate(_swap_flags(path, bad), "variants").valid


def test_unknown_flags_are_rejected(tmp_path: Path) -> None:
    """An item outside the Variants vocabulary raises (in ``vcf_records`` or the writer)."""
    aligns = peptides_at(GP, 44, sub(GP, 44, "S"), 1)
    with pytest.raises(ValueError, match="bogus"):
        _, table, _, _ = records(aligns, {"flags": [("bogus",)]})
        vt().write_variants(tmp_path / "variants.tsv", table)
