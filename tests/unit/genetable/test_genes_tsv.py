"""M4-15 ``write_genes_tsv`` / ``read_genes_tsv`` and the ``genes`` table schema
(docs/specs/M4-15-gene-table-assembly.md; GIPPA_SPEC.md §6.5, §6.6, §9.2, §9.3, §9.5; ADR-0006).

Pinned by the test author:
- The schema file is ``src/gippa/schemas/tables/genes.schema.json`` (lead decision 16, the
  ``<table>.schema.json`` convention); ``genes.json`` must not exist. Like the M3-16 / M4-16 table
  schemas it is a Frictionless *schema* descriptor with a top-level ``coordinates`` key, which is
  removed before ``frictionless.Schema.from_descriptor``.
- The ``flags`` vocabulary is enforced by the descriptor: a TSV whose ``flags`` cell holds a word
  outside ``GENE_FLAGS`` fails ``frictionless validate`` (pattern or enum, the encoding is free).
- ``read_genes_tsv`` returns the internal representation: ``start - 1``, list columns as lists
  (``NA`` -> empty list), booleans as bool, integers as an integer dtype, floats as float
  (6 significant digits survive, §6.6 ``%.6g``), ``NA`` as missing.
- The writer emits the columns in ``GENES_COLUMNS`` order whatever the frame's column order.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import gene_table_fixture as gt
from gene_table_fixture import as_list, is_na, tag

TABLE_DIR = Path(__file__).resolve().parents[3] / "src" / "gippa" / "schemas" / "tables"


@pytest.fixture(scope="module")
def table() -> pd.DataFrame:
    fr = gt.frames_with()
    fr["genes"].loc[5, "partial"] = True
    fr["layers"].at[5, "flags"] = ["partial"]
    fr["layers"].loc[0, ["cai", "spacer_bp"]] = [0.123456789, 7]
    return gt.build_world(frames=fr).assemble()


def schema_path() -> Path:
    path = TABLE_DIR / "genes.schema.json"
    assert path.is_file(), f"expected {path} (lead decision 16)"
    return path


def test_schema_file_name() -> None:
    """Lead decision 16: one ``genes.schema.json``, no ``genes.json`` left behind."""
    assert (TABLE_DIR / "genes.schema.json").is_file()
    assert not (TABLE_DIR / "genes.json").exists()


def descriptor() -> dict[str, Any]:
    return json.loads(schema_path().read_text(encoding="utf-8"))  # type: ignore[no-any-return]


def frictionless_report(path: Path) -> Any:
    frictionless = pytest.importorskip("frictionless")
    d = {k: v for k, v in descriptor().items() if k != "coordinates"}
    res = frictionless.Resource(
        path=path.name, basepath=str(path.parent), schema=frictionless.Schema.from_descriptor(d)
    )
    return res.validate()


def write(tmp_path: Path, df: pd.DataFrame, name: str = "genes.tsv") -> Path:
    from gippa.io.genes_tsv import write_genes_tsv

    path = tmp_path / name
    write_genes_tsv(path, df)
    return path


def raw(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


# ----------------------------------------------------------------------------- writer


def test_header_and_cell_conventions(tmp_path: Path, table: pd.DataFrame) -> None:
    path = write(tmp_path, table)
    text = path.read_text(encoding="utf-8")
    assert text.endswith("\n") and "\r" not in text
    assert tuple(text.split("\n", 1)[0].split("\t")) == gt.GENES_COLUMNS
    cells = raw(path).set_index("locus_tag")
    assert len(cells) == len(table)
    for c in gt.BOOL_COLUMNS:
        assert set(cells[c]) <= {"true", "false"}, c
    assert cells.loc[tag(5), "partial"] == "true" and cells.loc[tag(0), "partial"] == "false"
    assert cells.loc[tag(0), "cai"] == "0.123457"  # %.6g
    assert cells.loc[tag(0), "spacer_bp"] == "7"
    assert cells.loc[tag(5), "p_detect"] == "NA"
    assert cells.loc[tag(1), "gene"] == "NA"  # missing string
    assert cells.loc[tag(5), "flags"] == "partial"
    for i in range(20):
        assert int(cells.loc[tag(i), "start"]) == int(table.loc[i, "start"]) + 1
        assert int(cells.loc[tag(i), "end"]) == int(table.loc[i, "end"])


def test_writer_uses_column_order_of_genes_columns(tmp_path: Path, table: pd.DataFrame) -> None:
    a = write(tmp_path, table, "a.tsv")
    b = write(tmp_path, table[list(reversed(table.columns))], "b.tsv")
    assert a.read_bytes() == b.read_bytes()


def test_list_item_with_separator_is_rejected(tmp_path: Path, table: pd.DataFrame) -> None:
    bad = table.copy()
    bad.at[0, "go_terms"] = ["GO:1;GO:2"]
    with pytest.raises(ValueError):
        write(tmp_path, bad)


# ----------------------------------------------------------------------------- round trip


def _equal_cells(c: str, a: Any, b: Any) -> bool:
    if c in gt.LIST_COLUMNS:
        return as_list(a) == as_list(b)
    if is_na(a) or is_na(b):
        return is_na(a) and is_na(b)
    if c in gt.FLOAT_COLUMNS:
        return math.isclose(float(a), float(b), rel_tol=1e-5, abs_tol=0.0)
    if c in gt.INT_COLUMNS:
        return int(a) == int(b)
    if c in gt.BOOL_COLUMNS:
        return bool(a) is bool(b)
    return str(a) == str(b)


def test_round_trip_dtype_aware(tmp_path: Path, table: pd.DataFrame) -> None:
    """Acceptance: ``read_genes_tsv(write_genes_tsv(df)) == df`` on the fixture."""
    from gippa.io.genes_tsv import read_genes_tsv

    back = read_genes_tsv(write(tmp_path, table))
    assert tuple(back.columns) == gt.GENES_COLUMNS
    assert len(back) == len(table)
    for c in gt.INT_COLUMNS:
        assert pd.api.types.is_integer_dtype(back[c].dtype), (c, back[c].dtype)
    for c in gt.BOOL_COLUMNS:
        assert pd.api.types.is_bool_dtype(back[c].dtype), (c, back[c].dtype)
    for c in ("coverage_frac", "gene_q_value", "p_detect", "uniop_p_prev", "mw_da"):
        assert pd.api.types.is_float_dtype(back[c].dtype), (c, back[c].dtype)
    for c in gt.LIST_COLUMNS:
        assert all(isinstance(v, list) for v in back[c]), c
    for c in gt.GENES_COLUMNS:
        for i in range(len(table)):
            assert _equal_cells(c, back.at[i, c], table.at[i, c]), (c, i, back.at[i, c],
                                                                    table.at[i, c])  # fmt: skip
    assert int(back.loc[back["locus_tag"] == tag(0), "start"].iloc[0]) == 100  # 0-based again


def test_round_trip_is_idempotent_on_bytes(tmp_path: Path, table: pd.DataFrame) -> None:
    from gippa.io.genes_tsv import read_genes_tsv

    first = write(tmp_path, table, "a.tsv")
    second = write(tmp_path, read_genes_tsv(first), "b.tsv")
    assert first.read_bytes() == second.read_bytes()


def test_parquet_round_trip_keeps_missing_values(tmp_path: Path, table: pd.DataFrame) -> None:
    """Lead decision 11 (blocking): missing strings stay null through the gene-table Parquet
    (never the text ``"nan"``), also after a second write as the ``detectability`` step does."""
    import pyarrow.parquet as pq

    from gippa.io.genes_tsv import read_genes_parquet, write_genes_parquet

    first = tmp_path / "pre.parquet"
    write_genes_parquet(first, table)
    back = read_genes_parquet(first)
    second = tmp_path / "genes.parquet"
    write_genes_parquet(second, back)
    assert tuple(back.columns) == gt.GENES_COLUMNS
    for c in gt.GENES_COLUMNS:
        for i in range(len(table)):
            assert _equal_cells(c, back.at[i, c], table.at[i, c]), (c, i, back.at[i, c],
                                                                    table.at[i, c])  # fmt: skip
    for c in gt.STRING_COLUMNS:
        assert "nan" not in [str(v).lower() for v in back[c] if not is_na(v)], c
    t1, t2 = pq.read_table(first), pq.read_table(second)
    for c in (*gt.STRING_COLUMNS, *gt.FLOAT_COLUMNS, *gt.INT_COLUMNS):
        n_missing = int(table[c].map(is_na).sum())
        assert t1.column(c).null_count == t2.column(c).null_count == n_missing, c
    for c in ("island_id", "replichore", "cog_category", "gene"):
        assert t2.column(c).null_count > 0, c  # the fixture has missing values there


def test_read_cog_def_falls_back_from_utf8(tmp_path: Path) -> None:
    """Lead decision 13: the NCBI file holds cp1252 bytes (0x92); ``read_cog_def`` reads it."""
    from gippa.io.genes_tsv import read_cog_def

    path = tmp_path / "cog-20.def.tab"
    rows = [b"COG0593\tL\tDnaA\x92s ATPase\tDnaA\t\t\t\n", b"COG0002\tKT\tRegulator\tRegX\t\t\t\n"]
    path.write_bytes(b"".join(rows))
    with pytest.raises(UnicodeDecodeError):
        path.read_bytes().decode("utf-8")
    df = read_cog_def(path)
    assert df["cog_id"].tolist() == ["COG0593", "COG0002"]
    assert df["cog_category"].tolist() == ["L", "KT"]


# ----------------------------------------------------------------------------- schema


def test_schema_lists_every_column_with_type_and_description() -> None:
    d = descriptor()
    fields = d["fields"]
    assert tuple(f["name"] for f in fields) == gt.GENES_COLUMNS
    for f in fields:
        assert f.get("type"), f["name"]
        assert isinstance(f.get("description"), str) and f["description"].strip(), f["name"]
    assert d["coordinates"] == "1-based-closed"
    assert "NA" in d["missingValues"]


def test_schema_types_follow_section_9_3() -> None:
    types = {f["name"]: f["type"] for f in descriptor()["fields"]}
    want = (
        {c: "integer" for c in gt.INT_COLUMNS}
        | {c: "boolean" for c in gt.BOOL_COLUMNS}
        | {c: "number" for c in gt.FLOAT_COLUMNS}
        | {c: "string" for c in (*gt.STRING_COLUMNS, *gt.LIST_COLUMNS)}
    )
    bad = {c: (types.get(c), t) for c, t in want.items() if types.get(c) != t}
    assert not bad, bad


def test_schema_enums() -> None:
    fields = {f["name"]: f for f in descriptor()["fields"]}
    for c, values in gt.ENUMS.items():
        assert fields[c].get("constraints", {}).get("enum") == list(values), c
    assert fields["locus_tag"].get("constraints", {}).get("required") is True


def test_schema_never_says_expressed() -> None:
    text = schema_path().read_text(encoding="utf-8").lower()
    for f in descriptor()["fields"]:
        assert "express" not in f["name"]
        for v in f.get("constraints", {}).get("enum", []):
            assert "express" not in str(v).lower()
    assert "expressed" not in text


def test_written_fixture_validates(tmp_path: Path, table: pd.DataFrame) -> None:
    report = frictionless_report(write(tmp_path, table))
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def _with_flags(tmp_path: Path, table: pd.DataFrame, cell: str) -> Path:
    path = write(tmp_path, table)
    lines = path.read_text(encoding="utf-8").split("\n")
    header = lines[0].split("\t")
    cells = lines[1].split("\t")
    cells[header.index("flags")] = cell
    lines[1] = "\t".join(cells)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def test_schema_accepts_every_gene_flag(tmp_path: Path, table: pd.DataFrame) -> None:
    report = frictionless_report(_with_flags(tmp_path, table, ";".join(gt.GENE_FLAGS)))
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


@pytest.mark.parametrize("cell", ["partial;bogus", "expressed", "same_db", "partial;;pseudogene"])
def test_schema_rejects_flags_outside_the_vocabulary(
    tmp_path: Path, table: pd.DataFrame, cell: str
) -> None:
    assert not frictionless_report(_with_flags(tmp_path, table, cell)).valid
