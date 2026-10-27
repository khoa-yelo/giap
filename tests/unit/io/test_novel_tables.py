"""``tables/novel_orfs.tsv``: the Frictionless schema, the writer and reader, and the Data Package
resource (issue M7-07; docs/specs/M7-07-novel-orf-outputs.md rules 1-4, lead decision 1;
GIPPA_SPEC.md §9.3, §9.4, §9.5; ADR-0012 decision 3; M7-index lead decision 2 (new flag values
through the schema)).

Interface reading (test author): ``write_novel_orfs(path, table)`` writes the table returned by
``gippa.novel.outputs.orf_table``; ``read_novel_orfs(path)`` returns a frame with the same columns
that ``write_novel_orfs`` writes back to identical bytes. The descriptor is
``src/gippa/schemas/tables/novel_orfs.schema.json`` and ``src/gippa/schemas/datapackage.json`` lists
it as the resource ``novel_orfs`` at ``tables/novel_orfs.tsv``.
"""

from __future__ import annotations

import json
import re
import shutil
from pathlib import Path
from typing import Any

import pytest
from helpers.novel_outputs_fixture import (
    NOVEL_ORFS_COLUMNS,
    Case,
    ex1,
    ex1_ex6,
    mirror_case,
    novel_io,
    open3_case,
    orfs_of,
    scored,
    table,
    two_contigs,
    wrap_case,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

ROOT = Path(__file__).resolve().parents[3]
SCHEMAS = ROOT / "src" / "gippa" / "schemas"
TABLES = SCHEMAS / "tables"
NOVEL_CLASSES = {"nterm_extension", "cterm_extension", "intergenic", "cds_outframe", "antisense",
                 "pseudogene", "sorf"}  # fmt: skip
NOVEL_FLAGS = (
    "junction_spanned", "homopolymer_junction", "class_fdr_unestimated", "rbs_supported",
    "no_start_codon", "partial", "ostir_missing", "rbs_threshold_unestimated", "junction_unknown",
    "promoter_supported", "terminator_supported", "readthrough_candidate",
)  # fmt: skip
ID1 = "NOVEL_c1_43_p"


def _schema() -> dict[str, Any]:
    data = json.loads((TABLES / "novel_orfs.schema.json").read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def _fields() -> dict[str, dict[str, Any]]:
    return {f["name"]: f for f in _schema()["fields"]}


def _validate(path: Path) -> Any:
    frictionless = pytest.importorskip("frictionless")
    d = {k: v for k, v in _schema().items() if k != "coordinates"}
    res = frictionless.Resource(
        path=path.name, basepath=str(path.parent), format="tsv",
        schema=frictionless.Schema.from_descriptor(d),
    )  # fmt: skip
    return res.validate()


def _assert_valid(path: Path) -> None:
    report = _validate(path)
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def _write(
    tmp_path: Path, case: Case, sc_kw: dict[str, Any] | None = None, name: str = "t"
) -> Path:
    orfs = orfs_of(case)
    t = table(case, orfs, scored(orfs, **(sc_kw or {})))
    path = tmp_path / name / "novel_orfs.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    novel_io().write_novel_orfs(path, t)
    return path


# ---------------------------------------------------------------- schema


def test_schema_fields_in_column_order() -> None:
    s = _schema()
    assert tuple(f["name"] for f in s["fields"]) == NOVEL_ORFS_COLUMNS
    assert s["name"] == "novel_orfs"
    assert s["coordinates"] == "1-based-closed"
    assert "NA" in s["missingValues"]
    assert all(f.get("description") for f in s["fields"])


def test_schema_types_and_enums() -> None:
    f = _fields()
    assert set(f["class"]["constraints"]["enum"]) == NOVEL_CLASSES
    assert set(f["strand"]["constraints"]["enum"]) == {"+", "-"}
    for col in ("frame", "start_longest", "start_best_ostir", "end", "length_aa_longest",
                "length_aa_best", "n_peptides", "rank", "terminator_distance_nt"):  # fmt: skip
        assert f[col]["type"] == "integer", col
    for col in ("best_peptide_q", "ostir_tir", "promoter_score", "terminator_conf"):
        assert f[col]["type"] == "number", col
    for col in ("called", "promoter_supported", "terminator_downstream", "rbs_supported"):
        assert f[col]["type"] == "boolean", col
    for col in ("start_longest", "start_best_ostir", "end"):
        assert f[col]["constraints"]["minimum"] == 1, col
    for col in ("orf_id", "class", "contig", "strand", "end", "n_peptides", "class_fdr", "called"):
        assert f[col]["constraints"]["required"] is True, col
    assert re.fullmatch(f["orf_id"]["constraints"]["pattern"], "NOVEL_c1_43_p")
    assert re.fullmatch(f["orf_id"]["constraints"]["pattern"], "NOVEL_my.contig-2_878_m")
    assert not re.fullmatch(f["orf_id"]["constraints"]["pattern"], "NOVEL_c1_43_x")


def test_minus_strand_convention_is_described() -> None:
    """Lead decision 1: every other table is lowest/highest, so the descriptions must say so."""
    f = _fields()
    for col in ("start_longest", "start_best_ostir", "end"):
        assert "transcription" in f[col]["description"].lower(), col


def test_class_fdr_is_a_string_float_or_unestimated() -> None:
    c = _fields()["class_fdr"]
    assert c["type"] == "string"
    pat = c["constraints"]["pattern"]
    for ok in ("unestimated", "0.05", "1e-07", "0", "1", "0.0123457"):
        assert re.fullmatch(pat, ok), ok
    for bad in ("foo", "NaN-ish", "unestimated;0.1"):
        assert not re.fullmatch(pat, bad), bad


def test_flags_vocabulary_through_the_schema() -> None:
    pat = _fields()["flags"]["constraints"]["pattern"]
    for flag in NOVEL_FLAGS:
        assert re.fullmatch(pat, flag), flag
    assert re.fullmatch(pat, "class_fdr_unestimated;no_start_codon;ostir_missing;partial")
    for bad in ("bogus", "snv", "class_fdr_unestimated;bogus", "partial;"):
        assert not re.fullmatch(pat, bad), bad


def test_schema_is_a_valid_frictionless_schema() -> None:
    frictionless = pytest.importorskip("frictionless")
    d = {k: v for k, v in _schema().items() if k != "coordinates"}
    frictionless.Schema.from_descriptor(d)


def test_datapackage_lists_the_resource() -> None:
    pkg = json.loads((SCHEMAS / "datapackage.json").read_text(encoding="utf-8"))
    res = {r["name"]: r for r in pkg["resources"]}
    r = res["novel_orfs"]
    assert r["path"] == "tables/novel_orfs.tsv"
    assert r["schema"] == "tables/novel_orfs.schema.json"
    assert r["profile"] == "tabular-data-resource"
    assert r["format"] == "tsv"
    assert r["coordinates"] == "1-based-closed"


def test_datapackage_validates_with_the_written_table(tmp_path: Path) -> None:
    frictionless = pytest.importorskip("frictionless")
    pkg = json.loads((SCHEMAS / "datapackage.json").read_text(encoding="utf-8"))
    pkg["resources"] = [r for r in pkg["resources"] if r["name"] == "novel_orfs"]
    (tmp_path / "tables").mkdir()
    src = _write(tmp_path, ex1_ex6(), {"best": {ID1: 60}})
    shutil.copy(src, tmp_path / "tables" / "novel_orfs.tsv")
    shutil.copy(TABLES / "novel_orfs.schema.json", tmp_path / "tables")
    (tmp_path / "datapackage.json").write_text(json.dumps(pkg), encoding="utf-8")
    report = frictionless.validate(str(tmp_path / "datapackage.json"))
    assert report.valid, report.flatten(["type", "note"])


# ---------------------------------------------------------------- frictionless on written tables


SCENARIOS = {
    "ex1_ex6": lambda: (ex1_ex6(), {"best": {ID1: 60}}),
    "minus": lambda: (mirror_case(ex1()), {}),
    "origin": lambda: (wrap_case(), {}),
    "open3": lambda: (open3_case(), {}),
    "no_start": lambda: (ex1(atg=False, gtg=False), {}),
    "two_contigs": lambda: (two_contigs(), {"class_fdr": {ID1: 0.05}}),
    "scored": lambda: (ex1(), {"best": {ID1: 60}, "promoter_score": {ID1: 0.9},
                               "promoter_supported": {ID1: True}, "terminator_conf": {ID1: 0.5},
                               "terminator_downstream": {ID1: True},
                               "terminator_distance_nt": {ID1: 12}, "rank": {ID1: 1},
                               "rbs_supported": {ID1: False}}),
}  # fmt: skip


@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_written_table_validates(tmp_path: Path, name: str) -> None:
    case, kw = SCENARIOS[name]()
    _assert_valid(_write(tmp_path, case, kw))


def test_header_only_table_validates(tmp_path: Path) -> None:
    case = ex1()
    path = tmp_path / "novel_orfs.tsv"
    novel_io().write_novel_orfs(path, table(case, [], None))
    assert path.read_text(encoding="utf-8") == "\t".join(NOVEL_ORFS_COLUMNS) + "\n"
    _assert_valid(path)


@pytest.mark.parametrize(
    ("col", "bad"), [("class", "cds_inframe"), ("class_fdr", "lots"), ("flags", "bogus")]
)
def test_schema_rejects_bad_cells(tmp_path: Path, col: str, bad: str) -> None:
    path = _write(tmp_path, ex1())
    lines = path.read_text(encoding="utf-8").splitlines()
    cells = lines[1].split("\t")
    cells[NOVEL_ORFS_COLUMNS.index(col)] = bad
    path.write_text("\n".join([lines[0], "\t".join(cells)]) + "\n", encoding="utf-8")
    assert not _validate(path).valid


# ---------------------------------------------------------------- reader


def test_read_write_round_trip_is_byte_identical(tmp_path: Path) -> None:
    path = _write(tmp_path, two_contigs(), {"best": {ID1: 60}, "class_fdr": {ID1: 0.05}})
    back = novel_io().read_novel_orfs(path)
    assert tuple(back.columns) == NOVEL_ORFS_COLUMNS
    assert len(back) == 3
    again = tmp_path / "again.tsv"
    novel_io().write_novel_orfs(again, back)
    assert again.read_bytes() == path.read_bytes()


def test_read_header_only(tmp_path: Path) -> None:
    path = tmp_path / "novel_orfs.tsv"
    novel_io().write_novel_orfs(path, table(ex1(), [], None))
    back = novel_io().read_novel_orfs(path)
    assert tuple(back.columns) == NOVEL_ORFS_COLUMNS and len(back) == 0


_maybe = st.one_of(st.none(), st.floats(0, 1e4, allow_nan=False, width=32))


@settings(
    max_examples=30, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    best=st.sampled_from([None, 42, 60]),
    tir=_maybe,
    prom=_maybe,
    psup=st.one_of(st.none(), st.booleans()),
    tdist=st.one_of(st.none(), st.integers(-200, 500)),
    rank=st.one_of(st.none(), st.integers(1, 50)),
    fdr=st.one_of(st.just("unestimated"), st.floats(0, 1, allow_nan=False)),
    called=st.booleans(),
)
def test_round_trip_property(
    tmp_path_factory: Any, best: int | None, tir: float | None, prom: float | None,
    psup: bool | None, tdist: int | None, rank: int | None, fdr: object, called: bool,
) -> None:  # fmt: skip
    case = ex1()
    orfs = orfs_of(case)
    cols = {"ostir_tir": tir, "promoter_score": prom, "promoter_supported": psup,
            "terminator_distance_nt": tdist, "rank": rank, "class_fdr": fdr}  # fmt: skip
    sc = scored(orfs, best={ID1: best}, called={ID1: called},
                **{k: {ID1: v} for k, v in cols.items()})  # fmt: skip
    d = Path(tmp_path_factory.mktemp("rt"))
    path = d / "novel_orfs.tsv"
    novel_io().write_novel_orfs(path, table(case, orfs, sc))
    _assert_valid(path)
    again = d / "again.tsv"
    novel_io().write_novel_orfs(again, novel_io().read_novel_orfs(path))
    assert again.read_bytes() == path.read_bytes()
