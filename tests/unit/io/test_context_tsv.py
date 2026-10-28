"""M5-09 I/O: ``write_autocorrelation`` / ``write_islands``, their readers and the two table
schemas (docs/specs/M5-09-context-tables-step.md "Public interface", examples 2, 3, 9, the edge
cases and the acceptance item "examples 2, 3, 9 and round trips (dtype-aware) of both tables").

GIPPA_SPEC.md §9.2 (table conventions and test columns), §9.3 (the two column dictionaries), §9.5
(Frictionless resources), §6.5 (0-based half-open internally, 1-based closed in the writer; an
origin-crossing span keeps ``end > length``), ADR-0006 (``NA``, ``;`` lists, ``%.6g``).

Interface assumed (spec "Public interface"; the writers do not exist yet)::

    gippa.io.context_tsv.write_autocorrelation(path, df)   # AUTOCORR_COLUMNS, no coordinates
    gippa.io.context_tsv.read_autocorrelation(path)        # inverse
    gippa.io.context_tsv.write_islands(path, df)           # ISLAND_COLUMNS, start + 1
    gippa.io.context_tsv.read_islands(path)                # inverse, start - 1

Reader dtypes (the "dtype-aware" round trip): every ``integer`` field of the schema comes back as
pandas nullable ``Int64`` (``lag`` carries NA, edge case 2, and so may ``start`` ... ``n_perm`` in
a header-only table), every ``number`` field as ``float64``, string fields as ``object`` holding
``str``, and the two list-valued cells (``flags``, ``genes``) as tuples of ``str``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import context_tables_data as cd

from gippa.context.autocorrelation import AUTOCORR_COLUMNS
from gippa.context.islands import ISLAND_COLUMNS

TABLE_DIR = Path(__file__).resolve().parents[3] / "src" / "gippa" / "schemas" / "tables"
NA = float("nan")


# --------------------------------------------------------------------------- frames


def _autocorr_row(**over: Any) -> dict[str, Any]:
    """One estimable Moran's I row (every column of §9.3 filled)."""
    row: dict[str, Any] = {
        "contig": "chr1", "statistic": "morans_i", "lag": 1, "null_model": "gene_label",
        "value": 0.25, "null_mean": -0.004, "null_lo": -0.12, "null_hi": 0.11,
        "stat": "morans_i", "effect": 0.254, "ci_low": NA, "ci_high": NA,
        "p_value": 0.01, "q_value": 0.02, "method": "moran_i_permutation", "n_perm": 999,
        "flags": ("ci_unestimated",),
    }  # fmt: skip
    row.update(over)
    return row


def _island_row(**over: Any) -> dict[str, Any]:
    """One island row with the 0-based half-open span of example 2."""
    row: dict[str, Any] = {
        "island_id": "ISL_1", "direction": "high", "contig": "chr1", "start": 1000, "end": 9500,
        "tu_first": "TU_3", "tu_last": "TU_9", "n_tus": 7, "n_genes": 12,
        "genes": ("G_00003", "G_00004", "G_00005"), "mean_abundance_log2": 22.5,
        "detected_frac": 0.75, "scan_s": 5.125, "block_len": 3, "stat": "scan_max",
        "effect": 1.5, "ci_low": 0.9, "ci_high": 2.1, "p_value": 0.004, "q_value": NA,
        "method": "scan_max_stat", "null_model": "circular_block_bootstrap", "n_perm": 999,
        "flags": ("ci_few_units",),
    }  # fmt: skip
    row.update(over)
    return row


def _frame(rows: list[dict[str, Any]], columns: tuple[str, ...]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=list(columns))


def _cells(path: Path) -> tuple[list[str], list[list[str]]]:
    text = path.read_text(encoding="utf-8")
    assert text.endswith("\n")
    lines = text.rstrip("\n").split("\n")
    return lines[0].split("\t"), [line.split("\t") for line in lines[1:]]


def _cell(path: Path, column: str, row: int = 0) -> str:
    header, rows = _cells(path)
    return rows[row][header.index(column)]


# --------------------------------------------------------------------------- headers


def test_autocorrelation_header_is_the_module_constant(tmp_path: Path) -> None:
    from gippa.io.context_tsv import write_autocorrelation

    out = tmp_path / "autocorrelation.tsv"
    write_autocorrelation(out, _frame([_autocorr_row()], AUTOCORR_COLUMNS))
    assert _cells(out)[0] == list(AUTOCORR_COLUMNS) == list(cd.AUTOCORR_COLUMNS)


def test_islands_header_is_the_module_constant(tmp_path: Path) -> None:
    from gippa.io.context_tsv import write_islands

    out = tmp_path / "islands.tsv"
    write_islands(out, _frame([_island_row()], ISLAND_COLUMNS))
    assert _cells(out)[0] == list(ISLAND_COLUMNS) == list(cd.ISLAND_COLUMNS)


@pytest.mark.parametrize("table", ["autocorrelation", "islands"])
def test_empty_frame_writes_the_header_line_only(tmp_path: Path, table: str) -> None:
    """Example 4 at the writer: no island → ``islands.tsv`` is the header line only."""
    from gippa.io import context_tsv

    columns = AUTOCORR_COLUMNS if table == "autocorrelation" else ISLAND_COLUMNS
    out = tmp_path / f"{table}.tsv"
    getattr(context_tsv, f"write_{table}")(out, _frame([], columns))
    assert out.read_text(encoding="utf-8") == "\t".join(columns) + "\n"


# --------------------------------------------------------------------------- coordinates


def test_example_2_island_coordinates_are_one_based_closed(tmp_path: Path) -> None:
    """Example 2: 0-based half-open [1000, 9500) → TSV ``start`` 1001, ``end`` 9500."""
    from gippa.io.context_tsv import write_islands

    out = tmp_path / "islands.tsv"
    write_islands(out, _frame([_island_row(start=1000, end=9500)], ISLAND_COLUMNS))
    assert (_cell(out, "start"), _cell(out, "end")) == ("1001", "9500")


def test_example_3_origin_crossing_island_keeps_an_unwrapped_end(tmp_path: Path) -> None:
    """Example 3: an island crossing the origin of a 50,000-bp replicon (0-based 49,000 - 51,200)
    is written ``start`` 49001, ``end`` 51200 — ``end > length`` is kept (§6.5)."""
    from gippa.io.context_tsv import write_islands

    out = tmp_path / "islands.tsv"
    write_islands(out, _frame([_island_row(start=49_000, end=51_200)], ISLAND_COLUMNS))
    assert (_cell(out, "start"), _cell(out, "end")) == ("49001", "51200")


def test_islands_round_trip_restores_zero_based_coordinates(tmp_path: Path) -> None:
    from gippa.io.context_tsv import read_islands, write_islands

    df = _frame([_island_row(start=49_000, end=51_200)], ISLAND_COLUMNS)
    out = tmp_path / "islands.tsv"
    write_islands(out, df)
    back = read_islands(out)
    assert int(back.loc[0, "start"]) == 49_000
    assert int(back.loc[0, "end"]) == 51_200


def test_autocorrelation_has_no_coordinates_to_convert(tmp_path: Path) -> None:
    """The autocorrelation schema declares no ``coordinates`` key (spec "Public interface")."""
    assert "coordinates" not in _schema("autocorrelation")
    assert _schema("islands")["coordinates"] == "1-based-closed"


# --------------------------------------------------------------------------- cell formats


def test_lag_is_an_integer_cell_and_na_when_absent(tmp_path: Path) -> None:
    """Edge case: ``lag`` is an integer column with NA (an ``oriter_slope`` row has no lag)."""
    from gippa.io.context_tsv import write_autocorrelation

    oriter = _autocorr_row(
        statistic="oriter_slope", stat="oriter_slope_pooled", lag=NA, null_model="bootstrap",
        p_value=NA, q_value=NA, null_mean=NA, null_lo=NA, null_hi=NA,
        method="ols_replichore", flags=("descriptive",),
    )  # fmt: skip
    rows = [_autocorr_row(lag=1), oriter]
    out = tmp_path / "autocorrelation.tsv"
    write_autocorrelation(out, _frame(rows, AUTOCORR_COLUMNS))
    assert _cell(out, "lag", 0) == "1"
    assert _cell(out, "lag", 1) == "NA"
    assert _cell(out, "p_value", 1) == "NA"


def test_very_small_floats_use_six_significant_digits(tmp_path: Path) -> None:
    """Edge case: ``%.6g`` as everywhere (M1-04 ``write_tsv``)."""
    from gippa.io.context_tsv import write_autocorrelation

    out = tmp_path / "autocorrelation.tsv"
    write_autocorrelation(out, _frame([_autocorr_row(p_value=1e-7, value=1 / 3)],
                                      AUTOCORR_COLUMNS))  # fmt: skip
    assert _cell(out, "p_value") == "1e-07"
    assert _cell(out, "value") == "0.333333"


def test_island_gene_list_is_semicolon_joined(tmp_path: Path) -> None:
    from gippa.io.context_tsv import write_islands

    out = tmp_path / "islands.tsv"
    write_islands(out, _frame([_island_row(genes=("G_1", "G_2", "G_3"))], ISLAND_COLUMNS))
    assert _cell(out, "genes") == "G_1;G_2;G_3"


@pytest.mark.parametrize("empty", [(), []])
def test_empty_flags_are_na(tmp_path: Path, empty: Any) -> None:
    """§9.2: a row with no flag writes ``NA`` (``regions_stats.tsv`` convention), never ``''``."""
    from gippa.io.context_tsv import write_islands

    out = tmp_path / "islands.tsv"
    write_islands(out, _frame([_island_row(flags=empty)], ISLAND_COLUMNS))
    assert _cell(out, "flags") == "NA"


def test_several_flags_are_semicolon_joined(tmp_path: Path) -> None:
    from gippa.io.context_tsv import write_autocorrelation

    out = tmp_path / "autocorrelation.tsv"
    flags = ("ci_unestimated", "few_values", "not_estimable")
    write_autocorrelation(out, _frame([_autocorr_row(flags=flags)], AUTOCORR_COLUMNS))
    assert _cell(out, "flags") == "ci_unestimated;few_values;not_estimable"


# --------------------------------------------------------------------------- the flags trap


@pytest.mark.parametrize("table", ["autocorrelation", "islands"])
def test_flags_column_is_not_the_pandas_flags_attribute(tmp_path: Path, table: str) -> None:
    """``pandas.Series.flags`` shadows the ``flags`` *column*: ``row.flags`` on an ``iterrows``
    row (or an ``itertuples`` row built with ``name=`` renaming) silently returns a ``Flags``
    object, so a writer that reaches for the attribute writes
    ``<Flags(allows_duplicate_labels=True)>`` into every cell. The writer must use
    ``row["flags"]``."""
    from gippa.io import context_tsv

    if table == "autocorrelation":
        columns, rows = AUTOCORR_COLUMNS, [_autocorr_row(flags=("not_estimable",))]
    else:
        columns, rows = ISLAND_COLUMNS, [_island_row(flags=("ci_few_units",))]
    df = _frame(rows, columns)
    # the trap itself: the attribute is not the column
    assert not isinstance(df.iloc[0].flags, tuple | list | str)
    out = tmp_path / f"{table}.tsv"
    getattr(context_tsv, f"write_{table}")(out, df)
    assert _cell(out, "flags") == rows[0]["flags"][0]
    assert "Flags" not in out.read_text(encoding="utf-8")


@pytest.mark.parametrize("table", ["autocorrelation", "islands"])
def test_every_row_keeps_its_own_flags(tmp_path: Path, table: str) -> None:
    """No row loses or inherits another row's error label on the way through the writer."""
    from gippa.io import context_tsv

    if table == "autocorrelation":
        columns = AUTOCORR_COLUMNS
        want = [("ci_unestimated",), (), ("ci_unestimated", "not_estimable")]
        rows = [_autocorr_row(lag=i + 1, flags=f) for i, f in enumerate(want)]
    else:
        columns = ISLAND_COLUMNS
        want = [("ci_few_units",), (), ("ci_few_units",)]
        rows = [_island_row(island_id=f"ISL_{i + 1}", flags=f) for i, f in enumerate(want)]
    out = tmp_path / f"{table}.tsv"
    getattr(context_tsv, f"write_{table}")(out, _frame(rows, columns))
    header, cells = _cells(out)
    got = [cd.flag_set(r[header.index("flags")]) for r in cells]
    assert got == [set(f) for f in want]


# --------------------------------------------------------------------------- closed vocabulary


@pytest.mark.parametrize("table", ["autocorrelation", "islands"])
def test_example_9_unknown_flag_raises(tmp_path: Path, table: str) -> None:
    """Example 9: an unknown flag reaching the writer is a ``ValueError`` (§9.4 is closed)."""
    from gippa.io import context_tsv

    if table == "autocorrelation":
        columns, rows = AUTOCORR_COLUMNS, [_autocorr_row(flags=("wobbly",))]
    else:
        columns, rows = ISLAND_COLUMNS, [_island_row(flags=("wobbly",))]
    with pytest.raises(ValueError, match="wobbly"):
        getattr(context_tsv, f"write_{table}")(tmp_path / "t.tsv", _frame(rows, columns))


@pytest.mark.parametrize("table", ["autocorrelation", "islands"])
def test_every_vocabulary_flag_is_accepted(tmp_path: Path, table: str) -> None:
    from gippa.context.vocab import CONTEXT_FLAGS
    from gippa.io import context_tsv

    if table == "autocorrelation":
        columns, make = AUTOCORR_COLUMNS, _autocorr_row
    else:
        columns, make = ISLAND_COLUMNS, _island_row
    rows = [make(flags=(flag,)) for flag in sorted(CONTEXT_FLAGS)]
    out = tmp_path / f"{table}.tsv"
    getattr(context_tsv, f"write_{table}")(out, _frame(rows, columns))
    header, cells = _cells(out)
    assert [r[header.index("flags")] for r in cells] == sorted(CONTEXT_FLAGS)


# --------------------------------------------------------------------------- round trips


def _canonical(df: pd.DataFrame, lists: tuple[str, ...]) -> list[tuple[Any, ...]]:
    """Rows as comparable Python values: lists → tuples of str, NA → None, numbers → the float
    the TSV can carry (``%.6g``, ADR-0006: the round trip is exact to six significant digits)."""
    out = []
    for row in df.itertuples(index=False, name=None):
        cells: list[Any] = []
        for column, value in zip(df.columns, row, strict=True):
            if column in lists:
                cells.append(tuple(str(v) for v in (value or ())))
            elif isinstance(value, str):
                cells.append(value)
            elif cd.is_na(value):
                cells.append(None)
            else:
                cells.append(float(f"{float(value):.6g}"))
        out.append(tuple(cells))
    return out


def test_autocorrelation_round_trip(tmp_path: Path) -> None:
    from gippa.io.context_tsv import read_autocorrelation, write_autocorrelation

    rows = [
        _autocorr_row(lag=1),
        _autocorr_row(lag=2, null_model="tu_block", value=-0.0125, flags=()),
        _autocorr_row(statistic="join_count_bb", stat="join_count_bb", lag=3, value=17.0,
                      method="join_count_permutation",
                      flags=("ci_unestimated", "not_estimable"), p_value=NA, q_value=NA),
        _autocorr_row(statistic="oriter_slope", stat="oriter_slope_pooled", lag=NA,
                      null_model="bootstrap", value=-2.4317, effect=-2.4317, ci_low=-2.9,
                      ci_high=-1.95, null_mean=NA, null_lo=NA, null_hi=NA, p_value=NA,
                      q_value=NA, method="ols_replichore", n_perm=1000,
                      flags=("descriptive",)),
    ]  # fmt: skip
    df = _frame(rows, AUTOCORR_COLUMNS)
    out = tmp_path / "autocorrelation.tsv"
    write_autocorrelation(out, df)
    back = read_autocorrelation(out)
    assert list(back.columns) == list(AUTOCORR_COLUMNS)
    assert _canonical(back, ("flags",)) == _canonical(df, ("flags",))


def test_islands_round_trip(tmp_path: Path) -> None:
    from gippa.io.context_tsv import read_islands, write_islands

    low = _island_row(
        island_id="ISL_2", direction="low", start=49_000, end=51_200, genes=("G_1",), n_tus=2,
        n_genes=1, scan_s=-4.5, effect=-1.2, ci_low=NA, ci_high=NA,
        flags=("ci_few_units", "not_estimable"),
    )  # fmt: skip
    rows = [_island_row(), low]
    df = _frame(rows, ISLAND_COLUMNS)
    out = tmp_path / "islands.tsv"
    write_islands(out, df)
    back = read_islands(out)
    assert list(back.columns) == list(ISLAND_COLUMNS)
    assert _canonical(back, ("genes", "flags")) == _canonical(df, ("genes", "flags"))


@pytest.mark.parametrize(
    ("table", "columns"),
    [("autocorrelation", AUTOCORR_COLUMNS), ("islands", ISLAND_COLUMNS)],
)
def test_reader_dtypes_are_the_schema_types(tmp_path: Path, table: str,
                                            columns: tuple[str, ...]) -> None:  # fmt: skip
    """Dtype-aware round trip: ``integer`` → nullable ``Int64``, ``number`` → ``float64``,
    ``string`` → ``object``, and the list cells hold tuples of ``str``."""
    from gippa.io import context_tsv

    make = _autocorr_row if table == "autocorrelation" else _island_row
    out = tmp_path / f"{table}.tsv"
    getattr(context_tsv, f"write_{table}")(out, _frame([make()], columns))
    back = getattr(context_tsv, f"read_{table}")(out)
    types = {f["name"]: f["type"] for f in _schema(table)["fields"]}
    lists = ("flags", "genes")
    for column in columns:
        dtype = back[column].dtype
        if column in lists:
            assert isinstance(back.loc[0, column], tuple), column
        elif types[column] == "integer":
            assert str(dtype) == "Int64", (column, dtype)
        elif types[column] == "number":
            assert str(dtype) == "float64", (column, dtype)
        else:
            assert pd.api.types.is_object_dtype(dtype), (column, dtype)
            assert isinstance(back.loc[0, column], str), column


@pytest.mark.parametrize("table", ["autocorrelation", "islands"])
def test_reader_of_a_header_only_table(tmp_path: Path, table: str) -> None:
    from gippa.io import context_tsv

    columns = AUTOCORR_COLUMNS if table == "autocorrelation" else ISLAND_COLUMNS
    out = tmp_path / f"{table}.tsv"
    getattr(context_tsv, f"write_{table}")(out, _frame([], columns))
    back = getattr(context_tsv, f"read_{table}")(out)
    assert list(back.columns) == list(columns) and len(back) == 0


def test_producer_frames_round_trip(tmp_path: Path) -> None:
    """The merged producers' own frames (M5-05 ``AUTOCORR_COLUMNS`` dtypes, M5-06
    ``ISLAND_COLUMNS`` dtypes) survive the writer unchanged in value."""
    import islands_data as isd

    from gippa.context.autocorrelation import autocorrelation_rows
    from gippa.io.context_tsv import (
        read_autocorrelation,
        read_islands,
        write_autocorrelation,
        write_islands,
    )
    from gippa.rng import for_step

    cfg = isd.config(bootstrap=39, permutations=19, block_len=1, alpha=0.5, span_frac=0.30)
    rng = np.random.default_rng(3)
    genes, tus, replicons = isd.frames(
        [isd.sim_contig(rng, "chr1", 80, rho=0.3, island=(20, 5, 3.0))]
    )
    ac = autocorrelation_rows(genes, replicons, cfg, for_step("x", 1))
    islands, _ = isd.step(genes, isd.build_series(genes, tus, replicons, cfg), cfg,
                          for_step("y", 1))  # fmt: skip
    for name, df, lists in (("autocorrelation", ac, ("flags",)),
                            ("islands", islands, ("genes", "flags"))):  # fmt: skip
        out = tmp_path / f"{name}.tsv"
        write = write_autocorrelation if name == "autocorrelation" else write_islands
        read = read_autocorrelation if name == "autocorrelation" else read_islands
        write(out, df)
        assert _canonical(read(out), lists) == _canonical(df, lists)


# --------------------------------------------------------------------------- schemas (§9.5)


def _schema(name: str) -> dict[str, Any]:
    found = [p for p in (TABLE_DIR / f"{name}.json", TABLE_DIR / f"{name}.schema.json")
             if p.exists()]  # fmt: skip
    assert len(found) == 1, found
    data = json.loads(found[0].read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


@pytest.mark.parametrize(
    ("name", "columns"),
    [("autocorrelation", AUTOCORR_COLUMNS), ("islands", ISLAND_COLUMNS)],
)
def test_schema_fields_are_the_table_columns(name: str, columns: tuple[str, ...]) -> None:
    assert tuple(f["name"] for f in _schema(name)["fields"]) == columns


@pytest.mark.parametrize("name", ["autocorrelation", "islands"])
def test_schema_is_a_valid_frictionless_descriptor(name: str) -> None:
    frictionless = pytest.importorskip("frictionless")
    descriptor = {k: v for k, v in _schema(name).items() if k != "coordinates"}
    frictionless.Schema.from_descriptor(descriptor)


@pytest.mark.parametrize("name", ["autocorrelation", "islands"])
def test_schema_treats_na_as_missing(name: str) -> None:
    assert "NA" in _schema(name)["missingValues"]


@pytest.mark.parametrize("name", ["autocorrelation", "islands"])
def test_schema_describes_every_field(name: str) -> None:
    """§9.2: every column has a description and carries its unit in the name."""
    for field in _schema(name)["fields"]:
        assert field.get("description"), field["name"]


def test_autocorrelation_schema_enums_and_types() -> None:
    fields = {f["name"]: f for f in _schema("autocorrelation")["fields"]}
    assert set(fields["statistic"]["constraints"]["enum"]) == {
        "morans_i", "join_count_bb", "oriter_slope",
    }  # fmt: skip
    assert set(fields["null_model"]["constraints"]["enum"]) == {
        "gene_label", "tu_block", "bootstrap",
    }  # fmt: skip
    assert fields["lag"]["type"] == "integer"  # edge case: integer with missing values
    assert fields["n_perm"]["type"] == "integer"
    for column in cd.AUTOCORR_FLOATS:
        assert fields[column]["type"] == "number", column
    for column in ("p_value", "q_value"):
        assert fields[column]["constraints"]["minimum"] == 0
        assert fields[column]["constraints"]["maximum"] == 1


def test_islands_schema_enums_and_types() -> None:
    fields = {f["name"]: f for f in _schema("islands")["fields"]}
    assert set(fields["direction"]["constraints"]["enum"]) == {"high", "low"}
    assert fields["island_id"]["constraints"]["pattern"] == "ISL_[1-9][0-9]*"
    for column in cd.ISLAND_INTS:
        assert fields[column]["type"] == "integer", column
    for column in cd.ISLAND_FLOATS:
        assert fields[column]["type"] == "number", column
    assert fields["start"]["constraints"]["minimum"] == 1  # 1-based
    assert fields["genes"]["type"] in ("string", "array")


@pytest.mark.parametrize("name", ["autocorrelation", "islands"])
def test_written_table_validates_as_a_frictionless_resource(tmp_path: Path, name: str) -> None:
    """§9.5: the shipped schema validates the writer's own output."""
    frictionless = pytest.importorskip("frictionless")
    from gippa.io import context_tsv

    if name == "autocorrelation":
        columns = AUTOCORR_COLUMNS
        oriter = _autocorr_row(
            statistic="oriter_slope", stat="oriter_slope_pooled", lag=NA,
            null_model="bootstrap", null_mean=NA, null_lo=NA, null_hi=NA, p_value=NA,
            q_value=NA, method="ols_replichore", flags=("descriptive",),
        )  # fmt: skip
        join = _autocorr_row(statistic="join_count_bb", stat="join_count_bb", lag=2,
                             method="join_count_permutation")  # fmt: skip
        rows = [_autocorr_row(), join, oriter]
    else:
        columns = ISLAND_COLUMNS
        rows = [_island_row(), _island_row(island_id="ISL_2", direction="low", flags=())]
    out = tmp_path / f"{name}.tsv"
    getattr(context_tsv, f"write_{name}")(out, _frame(rows, columns))
    schema = frictionless.Schema.from_descriptor(
        {k: v for k, v in _schema(name).items() if k != "coordinates"}
    )
    resource = frictionless.Resource(path=out.name, basepath=str(out.parent), schema=schema)
    report = resource.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_autocorrelation_schema_states_the_q_value_rule() -> None:
    """Schema rule (R4): ``p_value`` NA ⇒ ``q_value`` NA. Frictionless cannot express a
    cross-field rule, so the schema states it in the description (as ``regions_stats`` does) and
    the data check lives in the integration test."""
    fields = {f["name"]: f for f in _schema("autocorrelation")["fields"]}
    assert "p_value" in fields["q_value"]["description"]


def test_islands_schema_pins_the_constant_columns() -> None:
    """Schema rules: ``method = scan_max_stat``, ``null_model`` the M5-06 constant (§9.3)."""
    from gippa.context.islands import METHOD, NULL_MODEL

    fields = {f["name"]: f for f in _schema("islands")["fields"]}
    assert fields["method"]["constraints"]["enum"] == [METHOD]
    assert NULL_MODEL in fields["null_model"]["constraints"]["enum"]
