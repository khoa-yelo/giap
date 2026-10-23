"""M4-06 input contract: ``giap.genetable.contract`` (docs/specs/M4-06-input-contract.md).

Examples 1-10 one test each, a parametrised missing-column test over every ``Column`` in ``TABLES``,
dtype/enum/NA rules, cross-table invariants and the edge cases of the spec. C15 ``features``
(M5-04, docs/specs/M5-04-sim-oric-gradient.md) is at the end: optional in ``validate_inputs``,
``contig`` known, ``0 <= start < end <= 2 length``.
"""

from __future__ import annotations

import ast
import dataclasses
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from m4_contract_fixture import (
    AGG_COLUMNS,
    FEATURE_COLUMNS,
    LAYER_COLUMNS,
    PEPTIDE_HITS_COLUMNS,
    PEPTIDE_RUN_COLUMNS,
    PROMOTER_COLUMNS,
    PROTEIN_HITS_COLUMNS,
    UNMAPPED_COLUMNS,
    _hit,
    build_frames,
    empty_like,
    make_inputs,
    tag,
)

from giap.io.issues import Issue

REPO = Path(__file__).resolve().parents[3]
CONTRACT_PY = REPO / "src" / "giap" / "genetable" / "contract.py"

TABLE_KEYS = (
    "replicons", "genes", "xrefs", "tus", "layers", "regions", "region_genes", "peptides",
    "abundance", "peptide_hits", "protein_hits", "unmapped", "id_map", "features",
    "peptide_runs", "peptide_spectra",
)  # fmt: skip

# Required (non-optional) columns per table, from the spec table "The contract".
REQUIRED: dict[str, tuple[str, ...]] = {
    "replicons": ("contig", "length", "topology", "kind"),
    "genes": (
        "locus_tag", "gene", "product", "feature_type", "contig", "start", "end", "strand",
        "partial", "crosses_origin", "hypothetical", "aa", "nt", "start_codon",
    ),
    "xrefs": ("locus_tag", "db", "id"),
    "tus": ("tu_id", "contig", "start", "end", "strand", "genes"),
    "layers": LAYER_COLUMNS,
    "regions": ("region_id", "region_type", "contig", "start", "end", "source_tool", "score"),
    "region_genes": ("region_id", "locus_tag", "overlap_frac", "subset", "core_biosynthetic"),
    "peptides": AGG_COLUMNS,
    "abundance": ("protein_id", "value"),  # lead decision 1 after review: no scale column
    "peptide_hits": PEPTIDE_HITS_COLUMNS,
    "protein_hits": PROTEIN_HITS_COLUMNS,
    "unmapped": UNMAPPED_COLUMNS,
    "id_map": ("source_id", "locus_tag"),
    "features": FEATURE_COLUMNS,  # C15 (M5-04): parse_bakta's columns, in its order
    "peptide_runs": PEPTIDE_RUN_COLUMNS,  # C16 (M6-02): M1-06's columns, in its order
}  # fmt: skip
OPTIONAL: dict[str, tuple[str, ...]] = {
    "genes": ("ribosomal_protein", "amr"),
    "abundance": ("raw_id",),  # merged M1-12 writes it (lead decision 1 after review)
}

# (table, column) -> (dtype, nullable or None = unpinned, enum or None)
PINNED: dict[tuple[str, str], tuple[str, bool | None, tuple[str, ...] | None]] = {
    ("replicons", "contig"): ("str", False, None),
    ("replicons", "length"): ("int", False, None),
    ("replicons", "topology"): ("str", False, ("circular", "linear")),
    ("replicons", "kind"): ("str", False, ("chromosome", "plasmid", "unknown")),
    ("genes", "locus_tag"): ("str", False, None),
    ("genes", "gene"): ("str", True, None),
    ("genes", "product"): ("str", True, None),
    ("genes", "feature_type"): ("str", False, ("cds", "sorf", "pseudogene")),
    ("genes", "contig"): ("str", False, None),
    ("genes", "start"): ("int", False, None),
    ("genes", "end"): ("int", False, None),
    ("genes", "strand"): ("str", False, ("+", "-")),
    ("genes", "partial"): ("bool", False, None),
    ("genes", "crosses_origin"): ("bool", False, None),
    ("genes", "hypothetical"): ("bool", False, None),
    ("genes", "aa"): ("str", True, None),
    ("genes", "nt"): ("str", False, None),
    ("genes", "start_codon"): ("str", None, None),
    ("genes", "ribosomal_protein"): ("bool", None, None),
    ("genes", "amr"): ("bool", None, None),
    ("xrefs", "locus_tag"): ("str", False, None),
    ("tus", "tu_id"): ("str", False, None),
    ("tus", "start"): ("int", False, None),
    ("tus", "end"): ("int", False, None),
    ("tus", "strand"): ("str", False, ("+", "-")),
    ("tus", "genes"): ("list", False, None),
    ("layers", "locus_tag"): ("str", False, None),
    ("layers", "tu_id"): ("str", False, None),
    ("layers", "tu_position"): ("int", None, None),
    ("layers", "tu_size"): ("int", None, None),
    ("layers", "uniop_p_prev"): ("float", True, None),
    ("layers", "uniop_p_next"): ("float", True, None),
    ("layers", "region_ids"): ("list", None, None),
    ("layers", "promoter_score"): ("float", True, None),
    ("layers", "cai"): ("float", True, None),
    ("layers", "flags"): ("list", None, None),
    ("regions", "region_id"): ("str", False, None),
    ("regions", "region_type"): (
        "str",
        False,
        ("prophage", "plasmid", "bgc", "crispr", "amr"),
    ),
    ("regions", "start"): ("int", True, None),
    ("regions", "end"): ("int", True, None),
    ("regions", "score"): ("float", True, None),
    ("region_genes", "overlap_frac"): ("float", False, None),
    ("region_genes", "subset"): (
        "str",
        True,
        ("integrase", "repressor", "structural", "lysis", "cargo"),
    ),
    ("region_genes", "core_biosynthetic"): ("bool", False, None),
    ("peptides", "peptide"): ("str", False, None),
    ("peptides", "q_value"): ("float", None, None),
    ("peptides", "is_decoy"): ("bool", False, None),
    ("peptides", "protein_ids"): ("list", None, None),
    ("peptides", "n_rows"): ("int", None, None),
    ("peptides", "intensity"): ("float", True, None),
    ("abundance", "protein_id"): ("str", False, None),
    ("abundance", "value"): ("float", True, None),
    ("abundance", "raw_id"): ("str", None, None),
    ("peptide_hits", "start"): ("int", False, None),
    ("peptide_hits", "end"): ("int", False, None),
    ("peptide_hits", "strand"): ("str", False, ("+", "-")),
    ("peptide_hits", "crosses_origin"): ("bool", False, None),
    ("peptide_hits", "route"): ("str", False, ("exact", "tolerant")),
    ("peptide_hits", "diff_types"): ("list", None, None),
    ("peptide_hits", "placement_q"): ("float", None, None),
    ("peptide_hits", "locus_class"): (
        "str",
        False,
        (
            "rna_overlap",
            "pseudogene",
            "sorf",
            "cds_inframe",
            "nterm_extension",
            "cterm_extension",  # ADR-0013 (lead decision 2 after review)
            "cds_outframe",
            "antisense",
            "intergenic",
        ),
    ),
    ("peptide_hits", "uniqueness"): ("str", False, ("unique", "multi_identical", "shared")),
    ("peptide_hits", "gene"): ("str", True, None),
    ("peptide_hits", "novel_candidate_id"): ("str", True, None),
    ("peptide_hits", "flags"): ("list", None, None),
    ("protein_hits", "rank"): ("int", None, None),
    ("protein_hits", "class"): ("str", False, ("identical", "near", "divergent")),
    ("protein_hits", "gene"): ("str", True, None),
    ("protein_hits", "flags"): ("list", None, None),
    # C15 (M5-04): ``?`` in the spec marks the nullable columns
    ("features", "feature_id"): ("str", False, None),
    ("features", "type"): (
        "str",
        False,
        ("rrna", "trna", "tmrna", "ncrna", "ncrna_region", "crispr", "oric", "orit", "oriv", "gap"),
    ),
    ("features", "contig"): ("str", False, None),
    ("features", "start"): ("int", False, None),
    ("features", "end"): ("int", False, None),
    ("features", "strand"): ("str", False, ("+", "-", ".")),
    ("features", "gene"): ("str", True, None),
    ("features", "product"): ("str", True, None),
    ("features", "amino_acid"): ("str", True, None),
    ("features", "anticodon"): ("str", True, None),
    # C16 (M6-02 lead decision 1): M1-06 ``peptide_runs``
    ("peptide_runs", "peptide"): ("str", False, None),
    ("peptide_runs", "run"): ("str", False, None),
    ("peptide_runs", "intensity"): ("float", True, None),
    ("peptide_runs", "n_rows"): ("int", False, None),
}


# --------------------------------------------------------------------------- helpers


def _contract() -> Any:
    from giap.genetable import contract

    return contract


def _all_columns() -> list[tuple[str, str]]:
    """Every (table, column) of ``TABLES`` (spec list until the module exists)."""
    try:
        from giap.genetable.contract import TABLES

        return [(t, c.name) for t, cols in TABLES.items() for c in cols]
    except Exception:
        return [(t, c) for t in TABLE_KEYS for c in (*REQUIRED[t], *OPTIONAL.get(t, ()))]


def errors(issues: list[Issue]) -> list[Issue]:
    return [i for i in issues if i.severity == "error"]


def warnings_(issues: list[Issue]) -> list[Issue]:
    return [i for i in issues if i.severity == "warn"]


def text(i: Issue) -> str:
    return i.message + " " + " ".join(i.examples)


def with_code(issues: list[Issue], code: str) -> list[Issue]:
    return [i for i in issues if i.code == code]


def validate_all(frames: dict[str, Any]) -> list[Issue]:
    """Per-table and cross-table validation, as ``load_m4_inputs`` would report."""
    c = _contract()
    out: list[Issue] = []
    for key in TABLE_KEYS:
        if frames.get(key) is not None:
            out += c.validate_table(key, frames[key])
    out += c.validate_inputs(with_features(make_inputs(frames), frames.get("features")))
    return out


class with_features:
    """``M4Inputs`` plus a ``features`` attribute (C15), whatever ``M4Inputs`` itself carries."""

    def __init__(self, base: Any, features: pd.DataFrame | None) -> None:
        self._base = base
        self.features = features

    def __getattr__(self, name: str) -> Any:
        return getattr(self._base, name)


@pytest.fixture
def frames() -> dict[str, Any]:
    return build_frames(20)


# --------------------------------------------------------------------------- public interface


def test_column_dataclass_shape() -> None:
    c = _contract()
    col = c.Column("x", "int")
    assert dataclasses.is_dataclass(col)
    assert (col.nullable, col.enum, col.optional) == (False, None, False)
    with pytest.raises(dataclasses.FrozenInstanceError):
        col.name = "y"  # type: ignore[misc]
    assert {f.name for f in dataclasses.fields(c.Column)} == {
        "name", "dtype", "nullable", "enum", "optional",
    }  # fmt: skip


def test_tables_keys() -> None:
    assert set(_contract().TABLES) == set(TABLE_KEYS)


@pytest.mark.parametrize("table", TABLE_KEYS)
def test_tables_columns_pinned(table: str) -> None:
    cols = _contract().TABLES[table]
    names = [c.name for c in cols]
    assert len(names) == len(set(names)), f"duplicate column in {table}"
    required = tuple(c.name for c in cols if not c.optional)
    if table in (
        "layers",
        "peptides",
        "peptide_hits",
        "protein_hits",
        "unmapped",
        "features",
        "peptide_runs",
    ):
        assert required == REQUIRED[table]  # producer column order (M2-20, M1-06, §9.3)
    else:
        assert set(required) == set(REQUIRED[table])
    optional = {c.name for c in cols if c.optional}
    assert set(OPTIONAL.get(table, ())) <= optional
    for c in cols:
        assert c.dtype in ("str", "int", "float", "bool", "list")


def test_gc_local_dev_if_present_is_optional() -> None:
    """Lead decision 2 adds ``gc_local_dev`` to M2-20's layers; it is not in M4-06's C6, so it may
    only be optional (absent in today's M2-20 LAYER_COLUMNS)."""
    cols = {c.name: c for c in _contract().TABLES["layers"]}
    if "gc_local_dev" in cols:
        assert cols["gc_local_dev"].optional and cols["gc_local_dev"].dtype == "float"


@pytest.mark.parametrize(("key", "spec"), list(PINNED.items()), ids=lambda x: str(x))
def test_column_types_pinned(key: tuple[str, str], spec: tuple[Any, ...]) -> None:
    table, name = key
    col = next(c for c in _contract().TABLES[table] if c.name == name)
    dtype, nullable, enum = spec
    assert col.dtype == dtype
    if nullable is not None:
        assert col.nullable is nullable
    if enum is not None:
        assert col.enum is not None and set(col.enum) == set(enum)


def test_layer_value_columns_nullable() -> None:
    """Example 9 context: every regulatory/codon layer column may be all-NA (layer missing)."""
    cols = {c.name: c for c in _contract().TABLES["layers"]}
    for name in LAYER_COLUMNS[7:]:
        if name not in ("start_codon", "flags"):
            assert cols[name].nullable, name


def test_valid_fixture_has_no_errors(frames: dict[str, Any]) -> None:
    issues = validate_all(frames)
    assert errors(issues) == []
    assert all(isinstance(i, Issue) for i in issues)


def test_issue_codes_are_namespaced(frames: dict[str, Any]) -> None:
    frames["genes"].loc[0, "feature_type"] = "CDS"
    frames["genes"].loc[1, "locus_tag"] = frames["genes"].loc[2, "locus_tag"]
    issues = validate_all(frames)
    assert issues
    for i in issues:
        assert isinstance(i, Issue)
        assert i.code.startswith("contract."), i.code
        assert i.code.split(".")[1] in TABLE_KEYS
        assert i.severity in ("error", "warn", "info")


# --------------------------------------------------------------------------- examples 1-10


def test_example1_200_gene_layout_loads_without_errors(tmp_path: Path) -> None:
    """Example 1 on a hand-built 200-gene layout (the giapsim variant is in tests/contract)."""
    from m4_contract_fixture import write_work

    from giap.io.work_tables import default_paths, load_m4_inputs

    f = build_frames(200)
    paths = default_paths(tmp_path / "work")
    write_work(f, paths)
    inputs, issues = load_m4_inputs(paths)
    assert errors(issues) == []
    assert len(inputs.genes) == 200


def test_example2_duplicate_locus_tag() -> None:
    f = build_frames(20)
    g = f["genes"]
    i5 = g.index[g["locus_tag"] == "GSIM_00005"][0]
    dup = g.loc[[i5]].copy()
    dup["start"] = dup["start"] + 7000
    dup["end"] = dup["end"] + 7000
    f["genes"] = pd.concat([g, dup], ignore_index=True)
    c = _contract()
    issues = c.validate_table("genes", f["genes"]) + c.validate_inputs(make_inputs(f))
    hit = with_code(issues, "contract.genes.locus_tag_unique")
    assert hit and all(i.severity == "error" for i in hit)
    assert any("GSIM_00005" in text(i) for i in hit)


def test_example3_feature_type_uppercase() -> None:
    f = build_frames(20)
    f["genes"].loc[0, "feature_type"] = "CDS"
    issues = _contract().validate_table("genes", f["genes"])
    hit = with_code(issues, "contract.genes.enum")
    assert hit and hit[0].severity == "error"
    assert any("feature_type" in text(i) for i in hit)


def test_example4_gene_missing_from_tus() -> None:
    f = build_frames(20)
    tus = f["tus"].copy()
    tus["genes"] = [[g for g in gs if g != "GSIM_00010"] for gs in tus["genes"]]
    tus = tus[[len(gs) > 0 for gs in tus["genes"]]].reset_index(drop=True)
    f["tus"] = tus
    issues = _contract().validate_inputs(make_inputs(f))
    hit = with_code(issues, "contract.tus.gene_coverage")
    assert hit and all(i.severity == "error" for i in hit)
    assert any("GSIM_00010" in text(i) for i in hit)


def test_example5_gene_in_two_tus() -> None:
    f = build_frames(20)
    tus = f["tus"].copy()
    first, last = tus.index[0], tus.index[-1]
    extra = tus.at[first, "genes"][0]
    tus.at[last, "genes"] = [*tus.at[last, "genes"], extra]
    f["tus"] = tus
    issues = _contract().validate_inputs(make_inputs(f))
    hit = with_code(issues, "contract.tus.gene_coverage")
    assert hit and all(i.severity == "error" for i in hit)
    assert any(extra in text(i) for i in hit)


def _origin_frames(topology: str) -> dict[str, Any]:
    f = build_frames(20)
    f["replicons"].loc[0, "length"] = 50_000
    f["replicons"].loc[0, "topology"] = topology
    g = f["genes"]
    last = g.index[-1]
    g.loc[last, "start"] = 49_900
    g.loc[last, "end"] = 50_200
    g.loc[last, "crosses_origin"] = True
    tus = f["tus"]
    t = next(i for i, gs in tus["genes"].items() if g.loc[last, "locus_tag"] in gs)
    tus.loc[t, "end"] = 50_200
    if len(tus.at[t, "genes"]) == 1:
        tus.loc[t, "start"] = 49_900
    return f


def test_example6_origin_crossing_circular_valid() -> None:
    issues = validate_all(_origin_frames("circular"))
    assert errors(issues) == []


def test_example6_origin_crossing_linear_error() -> None:
    issues = validate_all(_origin_frames("linear"))
    hit = with_code(issues, "contract.genes.origin")
    assert hit and all(i.severity == "error" for i in hit)


def test_example7_assigning_class_without_gene() -> None:
    f = build_frames(20)
    h = f["peptide_hits"]
    row = h.index[h["locus_class"] == "cds_inframe"][0]
    assert h.loc[row, "uniqueness"] == "unique"
    h.loc[row, "gene"] = None
    issues = _contract().validate_inputs(make_inputs(f))
    hit = with_code(issues, "contract.peptide_hits.assignment")
    assert hit and all(i.severity == "error" for i in hit)


def test_example8_unknown_gene() -> None:
    f = build_frames(20)
    h = f["peptide_hits"]
    row = h.index[h["locus_class"] == "cds_inframe"][0]
    h.loc[row, "gene"] = "GSIM_99999"
    issues = _contract().validate_inputs(make_inputs(f))
    hit = with_code(issues, "contract.peptide_hits.unknown_gene")
    assert hit and all(i.severity == "error" for i in hit)
    assert any("GSIM_99999" in text(i) for i in hit)


@pytest.mark.parametrize("encoding", ["nan_float", "none_object", "pd_na_nullable"])
def test_example9_all_promoter_columns_na_valid(encoding: str) -> None:
    f = build_frames(20)
    lay = f["layers"]
    n = len(lay)
    for c in PROMOTER_COLUMNS:
        if encoding == "nan_float":
            lay[c] = pd.Series([math.nan] * n, dtype="float64")
        elif encoding == "none_object":
            lay[c] = pd.Series([None] * n, dtype=object)
        else:
            lay[c] = pd.array([pd.NA] * n, dtype="Float64")
    c_ = _contract()
    assert errors(c_.validate_table("layers", lay)) == []
    assert errors(c_.validate_inputs(make_inputs(f))) == []


def test_example10_protein_ids_joined_string() -> None:
    """Example 10 as changed by lead decision 1 (spec "Lead decisions"): in the Parquet work tables
    ``protein_ids`` is a native list of strings; a ``;``-joined string (the TSV-only encoding of
    ADR-0006) is error ``contract.peptides.dtype``."""
    f = build_frames(20)
    p = f["peptides"]
    p["protein_ids"] = [";".join(x) if len(x) else None for x in p["protein_ids"]]
    p.loc[0, "protein_ids"] = "GSIM_00025;GSIM_00030"
    issues = _contract().validate_table("peptides", p)
    hit = with_code(issues, "contract.peptides.dtype")
    assert hit and all(i.severity == "error" for i in hit)
    assert any("protein_ids" in text(i) for i in hit)


def test_protein_ids_native_list_valid() -> None:
    f = build_frames(20)
    p = f["peptides"]
    assert all(isinstance(x, list) for x in p["protein_ids"])
    p.at[0, "protein_ids"] = ["GSIM_00025", "GSIM_00030"]
    assert errors(_contract().validate_table("peptides", p)) == []
    p["protein_ids"] = [np.array(x, dtype=object) for x in p["protein_ids"]]  # as read back
    assert errors(_contract().validate_table("peptides", p)) == []


def test_protein_ids_one_joined_cell_is_error() -> None:
    """A single ``;``-joined string among native lists is still a dtype error naming the column."""
    f = build_frames(20)
    p = f["peptides"]
    p.at[0, "protein_ids"] = "GSIM_00025;GSIM_00030"
    hit = with_code(_contract().validate_table("peptides", p), "contract.peptides.dtype")
    assert hit and all(i.severity == "error" for i in hit)
    assert any("protein_ids" in text(i) for i in hit)


def test_protein_ids_empty_list_valid() -> None:
    """M1-12 writes ``[]`` for a peptide without protein IDs (fixture row ``CCCCCCCR``)."""
    p = build_frames(20)["peptides"]
    assert list(p.loc[p["peptide"] == "CCCCCCCR", "protein_ids"].iloc[0]) == []
    assert errors(_contract().validate_table("peptides", p)) == []


@pytest.mark.parametrize(
    "bad",
    [
        ["GSIM_00025", ""],
        [""],
        ["GSIM_00025", None],
        [None],
        ["GSIM_00025", math.nan],
        [pd.NA],
        np.array(["GSIM_00025", ""], dtype=object),
        np.array(["GSIM_00025", None], dtype=object),
    ],
    ids=["empty_str", "only_empty", "none", "only_none", "nan", "pd_na", "np_empty", "np_none"],
)
def test_protein_ids_empty_or_na_element_is_error(bad: Any) -> None:
    """Lead decision 1: lists containing empty strings or NA elements are errors."""
    f = build_frames(20)
    p = f["peptides"]
    p.at[1, "protein_ids"] = bad
    errs = [
        i
        for i in errors(_contract().validate_table("peptides", p))
        if i.code.startswith("contract.peptides.")
    ]
    assert errs and any("protein_ids" in text(i) for i in errs)


# --------------------------------------------------------------------------- missing columns


@pytest.mark.parametrize(("table", "column"), _all_columns(), ids=lambda x: str(x))
def test_missing_column(table: str, column: str) -> None:
    f = build_frames(20)
    df = f[table]
    cols = {c.name: c for c in _contract().TABLES[table]}
    dropped = df.drop(columns=[column]) if column in df.columns else df
    issues = _contract().validate_table(table, dropped)
    named = [i for i in errors(issues) if column in text(i)]
    if cols[column].optional:
        assert named == [], f"optional {table}.{column} missing must not be an error"
    else:
        assert named, f"missing required {table}.{column} gave no error naming it"
        assert all(i.code.startswith(f"contract.{table}.") for i in named)


# --------------------------------------------------------------------------- dtypes, enums, NA


def test_extra_column_is_warning_only(frames: dict[str, Any]) -> None:
    g = frames["genes"].assign(my_extra=1)
    issues = _contract().validate_table("genes", g)
    assert errors(issues) == []
    w = warnings_(issues)
    assert w and any("my_extra" in text(i) for i in w)


def test_wrong_dtype_is_error(frames: dict[str, Any]) -> None:
    g = frames["genes"].copy()
    g["start"] = g["start"].astype(str) + "x"
    hit = with_code(_contract().validate_table("genes", g), "contract.genes.dtype")
    assert hit and hit[0].severity == "error" and any("start" in text(i) for i in hit)


def test_float_in_bool_column_is_dtype_error(frames: dict[str, Any]) -> None:
    g = frames["genes"].copy()
    g["partial"] = 0.5
    hit = with_code(_contract().validate_table("genes", g), "contract.genes.dtype")
    assert hit and any("partial" in text(i) for i in hit)


def test_na_in_non_nullable_column(frames: dict[str, Any]) -> None:
    g = frames["genes"].copy()
    g["locus_tag"] = g["locus_tag"].astype(object)
    g.loc[3, "locus_tag"] = None
    errs = errors(_contract().validate_table("genes", g))
    assert errs and any("locus_tag" in text(i) for i in errs)
    assert all(i.code.startswith("contract.genes.") for i in errs)


def test_na_in_nullable_columns_ok(frames: dict[str, Any]) -> None:
    g = frames["genes"].copy()
    g["gene"] = None
    g["product"] = None
    assert errors(_contract().validate_table("genes", g)) == []


def test_bad_strand_enum(frames: dict[str, Any]) -> None:
    g = frames["genes"].copy()
    g.loc[0, "strand"] = "."
    hit = with_code(_contract().validate_table("genes", g), "contract.genes.enum")
    assert hit and any("strand" in text(i) for i in hit)


@pytest.mark.parametrize(
    ("table", "column", "bad"),
    [
        ("replicons", "topology", "Circular"),
        ("replicons", "kind", "chromid"),
        ("regions", "region_type", "island"),
        ("region_genes", "subset", "tail"),
        ("peptide_hits", "route", "miniprot"),
        ("peptide_hits", "locus_class", "cds"),
        ("peptide_hits", "uniqueness", "razor"),
        ("protein_hits", "class", "exact"),
    ],
)
def test_bad_enum_values(table: str, column: str, bad: str) -> None:
    f = build_frames(20)
    df = f[table].copy()
    df[column] = df[column].astype(object)
    df.loc[df.index[0], column] = bad
    hit = with_code(_contract().validate_table(table, df), f"contract.{table}.enum")
    assert hit and hit[0].severity == "error"
    assert any(column in text(i) for i in hit)


def test_enum_ignores_na_in_nullable_enum(frames: dict[str, Any]) -> None:
    rg = frames["region_genes"]
    assert rg["subset"].isna().any()
    assert (
        with_code(_contract().validate_table("region_genes", rg), "contract.region_genes.enum")
        == []
    )


def test_pandas_nullable_dtypes_accepted(frames: dict[str, Any]) -> None:
    """Edge case: ``Int64``, ``boolean`` and ``string`` are accepted like their numpy twins."""
    c = _contract()
    for key in ("replicons", "genes", "xrefs", "regions", "region_genes", "abundance", "unmapped"):
        df = frames[key].convert_dtypes()
        assert errors(c.validate_table(key, df)) == [], key
    g = frames["genes"].astype(
        {"locus_tag": "string", "start": "Int64", "end": "Int64", "partial": "boolean"}
    )
    assert errors(c.validate_table("genes", g)) == []


def test_numpy_int32_and_object_strings_accepted(frames: dict[str, Any]) -> None:
    g = frames["genes"].astype({"start": np.int32, "end": np.int32, "locus_tag": object})
    assert errors(_contract().validate_table("genes", g)) == []


def test_numpy_array_list_cells_accepted(frames: dict[str, Any]) -> None:
    """Parquet ``list<string>`` read by pandas gives numpy arrays per cell."""
    tus = frames["tus"].copy()
    tus["genes"] = [np.array(g, dtype=object) for g in tus["genes"]]
    assert errors(_contract().validate_table("tus", tus)) == []


def test_validate_table_does_not_mutate(frames: dict[str, Any]) -> None:
    g = frames["genes"]
    before = g.copy(deep=True)
    _contract().validate_table("genes", g)
    pd.testing.assert_frame_equal(g, before)


def test_header_only_tables_valid(frames: dict[str, Any]) -> None:
    """Edge case: header-only C11 (no peptide placed) and C12 (no proteins) are valid."""
    frames["peptide_hits"] = empty_like(frames["peptide_hits"])
    frames["protein_hits"] = empty_like(frames["protein_hits"])
    frames["unmapped"] = empty_like(frames["unmapped"])
    assert errors(validate_all(frames)) == []


def test_header_only_object_dtype_valid() -> None:
    """A header-only TSV read gives object columns: no dtype errors on zero rows."""
    df = pd.DataFrame({c: pd.Series([], dtype=object) for c in PEPTIDE_HITS_COLUMNS})
    assert errors(_contract().validate_table("peptide_hits", df)) == []


# --------------------------------------------------------------------------- C1/C2 invariants


def test_zero_gene_replicon_valid(frames: dict[str, Any]) -> None:
    """Edge case: a tiny plasmid with no genes is valid and has no TUs."""
    frames["replicons"] = pd.concat(
        [
            frames["replicons"],
            pd.DataFrame(
                {"contig": ["p1"], "length": [800], "topology": ["circular"], "kind": ["plasmid"]}
            ).astype({"length": "int64"}),
        ],
        ignore_index=True,
    )
    assert errors(validate_all(frames)) == []


def test_duplicate_replicon_contig(frames: dict[str, Any]) -> None:
    frames["replicons"] = pd.concat([frames["replicons"]] * 2, ignore_index=True)
    errs = errors(validate_all(frames))
    assert any(i.code.startswith("contract.replicons.") and "chr1" in text(i) for i in errs)


def _gene_error(frames: dict[str, Any]) -> list[Issue]:
    return [i for i in errors(validate_all(frames)) if i.code.startswith("contract.genes.")]


def test_start_not_below_end(frames: dict[str, Any]) -> None:
    frames["genes"].loc[4, "end"] = frames["genes"].loc[4, "start"]
    assert _gene_error(frames)


def test_negative_start(frames: dict[str, Any]) -> None:
    frames["genes"].loc[0, "start"] = -3
    assert _gene_error(frames)


def test_end_beyond_length_without_origin_flag(frames: dict[str, Any]) -> None:
    length = int(frames["replicons"].loc[0, "length"])
    last = frames["genes"].index[-1]
    frames["genes"].loc[last, "start"] = length - 100
    frames["genes"].loc[last, "end"] = length + 200
    assert _gene_error(frames)


def test_end_equal_length_valid(frames: dict[str, Any]) -> None:
    """Half-open: ``end == length`` is the last base, not an overflow."""
    last = frames["genes"].index[-1]
    end = int(frames["genes"].loc[last, "end"])
    frames["replicons"].loc[0, "length"] = end
    frames["peptide_hits"] = frames["peptide_hits"][
        frames["peptide_hits"]["locus_class"] != "intergenic"
    ].reset_index(drop=True)
    # the C15 fixture rows sit in the free stretch the shrunk replicon no longer has (C15 is
    # optional, and its own bounds are tested below)
    frames["features"] = empty_like(frames["features"])
    assert errors(validate_all(frames)) == []


@pytest.mark.parametrize("end", [100_000, 100_050, 50_000])
def test_origin_end_out_of_range(end: int) -> None:
    """``crosses_origin`` needs ``length < end < 2·length`` (L = 50,000)."""
    f = _origin_frames("circular")
    last = f["genes"].index[-1]
    f["genes"].loc[last, "end"] = end
    hit = [i for i in errors(validate_all(f)) if i.code.startswith("contract.genes.")]
    assert hit


def test_gene_on_unknown_contig(frames: dict[str, Any]) -> None:
    frames["genes"].loc[0, "contig"] = "chrX"
    assert _gene_error(frames)


def test_cds_with_na_aa_is_error(frames: dict[str, Any]) -> None:
    frames["genes"].loc[4, "aa"] = None
    assert _gene_error(frames)


def test_pseudogene_na_aa_with_nt_valid(frames: dict[str, Any]) -> None:
    g = frames["genes"]
    assert g.loc[7, "feature_type"] == "pseudogene" and pd.isna(g.loc[7, "aa"])
    assert errors(validate_all(frames)) == []


def test_pseudogene_without_nt_is_error(frames: dict[str, Any]) -> None:
    """Edge case: a pseudogene with ``aa`` NA must have ``nt`` (M4-10 translates it)."""
    frames["genes"].loc[7, "nt"] = None
    errs = [i for i in errors(validate_all(frames)) if i.code.startswith("contract.genes.")]
    assert errs and any("nt" in text(i) for i in errs)


def test_aa_with_trailing_stop_is_error(frames: dict[str, Any]) -> None:
    frames["genes"].loc[4, "aa"] = frames["genes"].loc[4, "aa"] + "*"
    assert _gene_error(frames)


# --------------------------------------------------------------------------- C5/C6 invariants


def test_tu_listing_unknown_gene(frames: dict[str, Any]) -> None:
    tus = frames["tus"]
    tus.at[tus.index[0], "genes"] = [*tus.at[tus.index[0], "genes"], "GSIM_77777"]
    errs = errors(validate_all(frames))
    assert any(i.code.startswith("contract.tus.") and "GSIM_77777" in text(i) for i in errs)


def test_singletons_are_tus() -> None:
    f = build_frames(12, tu_sizes=(1,))
    assert errors(validate_all(f)) == []


def test_layers_missing_gene_row(frames: dict[str, Any]) -> None:
    frames["layers"] = frames["layers"].iloc[1:].reset_index(drop=True)
    errs = errors(validate_all(frames))
    assert any(i.code.startswith("contract.layers.") and tag(0) in text(i) for i in errs)


def test_layers_duplicate_gene_row(frames: dict[str, Any]) -> None:
    lay = frames["layers"]
    frames["layers"] = pd.concat([lay, lay.iloc[[2]]], ignore_index=True)
    errs = errors(validate_all(frames))
    assert any(i.code.startswith("contract.layers.") for i in errs)


def test_layers_tu_id_disagrees_with_tus(frames: dict[str, Any]) -> None:
    lay = frames["layers"]
    i = lay.index[lay["tu_id"] != "TU_1"][0]
    lay.loc[i, "tu_id"] = "TU_1"
    errs = errors(validate_all(frames))
    assert any(i.code.startswith("contract.layers.") for i in errs)


# --------------------------------------------------------------------------- C7/C8 invariants


def test_amr_region_with_coordinates_is_error(frames: dict[str, Any]) -> None:
    r = frames["regions"]
    i = r.index[r["region_type"] == "amr"][0]
    r.loc[i, "start"] = 0
    r.loc[i, "end"] = 1000
    errs = errors(validate_all(frames))
    assert any(i.code.startswith("contract.regions.") for i in errs)


def test_subset_on_non_prophage_is_error(frames: dict[str, Any]) -> None:
    rg = frames["region_genes"]
    i = rg.index[rg["region_id"] == "amr_1"][0]
    rg.loc[i, "subset"] = "cargo"
    errs = errors(validate_all(frames))
    assert any(i.code.startswith("contract.region_genes.") for i in errs)


def test_prophage_member_without_subset_is_error(frames: dict[str, Any]) -> None:
    rg = frames["region_genes"]
    i = rg.index[rg["region_id"] == "prophage_1"][0]
    rg.loc[i, "subset"] = None
    errs = errors(validate_all(frames))
    assert any(i.code.startswith("contract.region_genes.") for i in errs)


# --------------------------------------------------------------------------- C9/C11/C12 invariants


def test_duplicate_peptide_in_c9(frames: dict[str, Any]) -> None:
    p = frames["peptides"]
    frames["peptides"] = pd.concat([p, p.iloc[[0]]], ignore_index=True)
    errs = errors(validate_all(frames))
    assert any(i.code.startswith("contract.peptides.") and "AAAAAAAK" in text(i) for i in errs)


def test_decoys_included_in_c9_valid(frames: dict[str, Any]) -> None:
    assert frames["peptides"]["is_decoy"].any()
    assert errors(validate_all(frames)) == []


def test_decoy_in_peptide_hits_is_error(frames: dict[str, Any]) -> None:
    hit = _hit("KAAAAAAA", 5000, "+", "intergenic", "unique", None)
    h = frames["peptide_hits"]
    frames["peptide_hits"] = pd.concat(
        [h, pd.DataFrame([hit]).astype(h.dtypes.to_dict())], ignore_index=True
    )
    errs = errors(validate_all(frames))
    assert any(i.code.startswith("contract.peptide_hits.") and "KAAAAAAA" in text(i) for i in errs)


def test_non_assigning_class_with_gene_is_error(frames: dict[str, Any]) -> None:
    """``gene`` is non-NA exactly for assigning classes."""
    h = frames["peptide_hits"]
    i = h.index[h["locus_class"] == "intergenic"][0]
    h.loc[i, "gene"] = tag(0)
    hit = with_code(validate_all(frames), "contract.peptide_hits.assignment")
    assert hit and hit[0].severity == "error"


@pytest.mark.parametrize("locus_class", ["sorf", "nterm_extension", "cds_inframe"])
def test_assigning_classes_with_gene_valid(frames: dict[str, Any], locus_class: str) -> None:
    h = frames["peptide_hits"]
    i = h.index[h["locus_class"] == "cds_inframe"][0]
    h.loc[i, "locus_class"] = locus_class
    assert errors(validate_all(frames)) == []


def test_pseudogene_placement_assigned_valid(frames: dict[str, Any]) -> None:
    """Lead decision 1 (M4-06 Q3 stands): pseudogene placements carry the pseudogene's tag."""
    h = frames["peptide_hits"]
    row = h[h["locus_class"] == "pseudogene"].iloc[0]
    assert row["gene"] == tag(7)
    assert errors(validate_all(frames)) == []


@pytest.mark.parametrize("uniqueness", ["multi_identical", "shared"])
def test_non_unique_assigned_rows_valid(frames: dict[str, Any], uniqueness: str) -> None:
    h = frames["peptide_hits"]
    i = h.index[h["locus_class"] == "cds_inframe"][0]
    h.loc[i, "uniqueness"] = uniqueness
    assert errors(validate_all(frames)) == []


def test_protein_hit_unknown_gene(frames: dict[str, Any]) -> None:
    frames["protein_hits"].loc[0, "gene"] = "GSIM_88888"
    errs = errors(validate_all(frames))
    assert any(
        i.code.startswith("contract.protein_hits.") and "GSIM_88888" in text(i) for i in errs
    )


def test_protein_hit_without_gene_valid(frames: dict[str, Any]) -> None:
    frames["protein_hits"].loc[0, "gene"] = None
    assert errors(validate_all(frames)) == []


def test_optional_tables_none_valid(frames: dict[str, Any]) -> None:
    """C10 and C14 absent (``None``) are valid; C14 absent → ``source_id`` NA downstream."""
    frames["abundance"] = None
    frames["id_map"] = None
    assert errors(validate_all(frames)) == []


def test_inf_score_input_valid(frames: dict[str, Any]) -> None:
    """ADR-0006: ``INF`` decodes to +inf; an infinite engine score is a float, not a dtype error."""
    frames["peptide_hits"].loc[0, "score_input"] = math.inf
    assert errors(validate_all(frames)) == []


# --------------------------------------------------------------------------- purity


def test_contract_module_does_no_io() -> None:
    """Acceptance: ``contract.py`` imports no I/O (file I/O only in ``io/work_tables.py``)."""
    src = CONTRACT_PY.read_text(encoding="utf-8")
    tree = ast.parse(src)
    banned_mods = ("pyarrow", "pathlib", "json", "io", "os", "shutil", "csv", "gzip")
    banned_giap = ("giap.io.tsv", "giap.io.parquet", "giap.io._open", "giap.io.work_tables")
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            names = [node.module or ""]
        else:
            continue
        for n in names:
            assert n.split(".")[0] not in banned_mods, n
            assert not any(n == b or n.startswith(b + ".") for b in banned_giap), n
        if isinstance(node, ast.ImportFrom) and node.module == "giap.io":
            assert {a.name for a in node.names} <= {"issues"}
    calls = {
        n.func.attr if isinstance(n.func, ast.Attribute) else getattr(n.func, "id", "")
        for n in ast.walk(tree)
        if isinstance(n, ast.Call)
    }
    assert not calls & {"open", "read_parquet", "to_parquet", "read_csv", "read_table", "read_tsv"}


def test_contract_uses_shared_issue_type() -> None:
    """``Issue`` is M1-04's shared type, not a local copy."""
    src = CONTRACT_PY.read_text(encoding="utf-8")
    assert "giap.io.issues" in src


# --------------------------------------------------------------------------- C15 features (M5-04)


def _features_issues(issues: list[Issue]) -> list[Issue]:
    return [i for i in issues if i.code.startswith("contract.features.")]


def _inputs_with(frames: dict[str, Any], features: pd.DataFrame | None) -> Any:
    return with_features(make_inputs(frames), features)


def test_c15_columns_follow_parse_bakta() -> None:
    """C15 is exactly the columns ``parse_bakta`` writes (``giap.steps.annotation._FEATURE_COLS``),
    in that order, all required (the table itself is optional)."""
    from giap.steps.annotation import _FEATURE_COLS

    cols = _contract().TABLES["features"]
    assert tuple(c.name for c in cols) == FEATURE_COLUMNS
    assert tuple(name for name, _ in _FEATURE_COLS) == FEATURE_COLUMNS
    assert not any(c.optional for c in cols)
    assert _contract().CONTRACT_COLUMNS["features"] == FEATURE_COLUMNS


def test_c15_type_enum_is_nctype() -> None:
    from typing import get_args

    from giap.model.genes import NcType

    col = next(c for c in _contract().TABLES["features"] if c.name == "type")
    assert col.enum is not None and set(col.enum) == set(get_args(NcType))


def test_c15_fixture_valid(frames: dict[str, Any]) -> None:
    assert _contract().validate_table("features", frames["features"]) == []
    assert _features_issues(_contract().validate_inputs(_inputs_with(frames, frames["features"])))\
        == []  # fmt: skip


def test_c15_header_only_valid(frames: dict[str, Any]) -> None:
    """No oriC and no other feature: an empty table is valid (M5-08 then reports ``no_oric``)."""
    empty = empty_like(frames["features"])
    assert _contract().validate_table("features", empty) == []
    assert _features_issues(_contract().validate_inputs(_inputs_with(frames, empty))) == []


def test_c15_absent_is_valid(frames: dict[str, Any]) -> None:
    """Example 7 (optional): ``features = None`` or no attribute at all gives no issue."""
    c = _contract()
    assert _features_issues(c.validate_inputs(_inputs_with(frames, None))) == []
    assert _features_issues(c.validate_inputs(make_inputs(frames))) == []


@pytest.mark.parametrize(
    ("column", "bad"),
    [("type", "oriC"), ("type", "cds"), ("strand", "?"), ("strand", "")],
)
def test_c15_enum_errors(frames: dict[str, Any], column: str, bad: str) -> None:
    """Bakta's raw ``oriC`` / ``?`` must be normalised by ``parse_bakta`` (``oric`` / ``.``)."""
    df = frames["features"].copy()
    df.loc[0, column] = bad
    errs = errors(_contract().validate_table("features", df))
    assert errs and all(i.code.startswith("contract.features.") for i in errs)
    assert any(column in text(i) for i in errs)


@pytest.mark.parametrize("column", ["feature_id", "type", "contig", "start", "end", "strand"])
def test_c15_na_in_required_value_is_error(frames: dict[str, Any], column: str) -> None:
    df = frames["features"].copy()
    df[column] = df[column].astype(object)
    df.loc[0, column] = None
    errs = errors(_contract().validate_table("features", df))
    assert errs and any(column in text(i) for i in errs)


@pytest.mark.parametrize("column", ["gene", "product", "amino_acid", "anticodon"])
def test_c15_nullable_columns_accept_na(frames: dict[str, Any], column: str) -> None:
    df = frames["features"].copy()
    df[column] = None
    assert errors(_contract().validate_table("features", df)) == []


def _feature_rows(frames: dict[str, Any], **row: Any) -> pd.DataFrame:
    df = frames["features"].iloc[:1].copy()
    for k, v in row.items():
        df[k] = df[k].astype(object) if isinstance(v, str) else df[k]
        df.loc[df.index[0], k] = v
    return df.astype({"start": "int64", "end": "int64"})


def _length(frames: dict[str, Any]) -> int:
    return int(frames["replicons"].loc[0, "length"])


@pytest.mark.parametrize(
    "row",
    [
        {"contig": "chrX"},  # unknown contig
        {"start": -1, "end": 231},  # start < 0
        {"start": 500, "end": 500},  # start == end
        {"start": 600, "end": 500},  # start > end
    ],
    ids=["unknown_contig", "negative_start", "empty", "reversed"],
)
def test_c15_out_of_range_rows_rejected(frames: dict[str, Any], row: dict[str, Any]) -> None:
    df = _feature_rows(frames, **row)
    issues = _contract().validate_table("features", df)
    issues += _contract().validate_inputs(_inputs_with(frames, df))
    errs = [i for i in errors(issues) if i.code.startswith("contract.features.")]
    assert errs, row


def test_c15_end_beyond_twice_length_rejected(frames: dict[str, Any]) -> None:
    n = _length(frames)
    df = _feature_rows(frames, start=2 * n - 100, end=2 * n + 1)
    issues = _contract().validate_table("features", df)
    issues += _contract().validate_inputs(_inputs_with(frames, df))
    assert [i for i in errors(issues) if i.code.startswith("contract.features.")]


@pytest.mark.parametrize("offset", [1, 116, None])
def test_c15_origin_crossing_feature_valid(frames: dict[str, Any], offset: int | None) -> None:
    """On a circular replicon ``end`` may pass ``length`` (unwrapped across the origin, §6.5) up
    to ``2 length``; ``start`` stays inside the replicon (lead decision 23)."""
    n = _length(frames)
    end = 2 * n if offset is None else n + offset
    start = n - 1 if offset is None else end - 232
    df = _feature_rows(frames, start=start, end=end)
    issues = _contract().validate_table("features", df)
    issues += _contract().validate_inputs(_inputs_with(frames, df))
    assert _features_issues(errors(issues)) == []


def test_c15_feature_on_plasmid_and_linear_contig_valid(frames: dict[str, Any]) -> None:
    reps = pd.concat(
        [frames["replicons"],
         pd.DataFrame({"contig": ["p1", "ctg1"], "length": [5000, 4000],
                       "topology": ["circular", "linear"], "kind": ["plasmid", "unknown"]})],
        ignore_index=True,
    ).astype({"length": "int64"})  # fmt: skip
    frames["replicons"] = reps
    df = pd.concat(
        [frames["features"], _feature_rows(frames, contig="p1", start=10, end=242),
         _feature_rows(frames, contig="ctg1", start=0, end=4000)],
        ignore_index=True,
    )  # fmt: skip
    df["feature_id"] = [f"f{i}" for i in range(len(df))]
    issues = _contract().validate_table("features", df)
    issues += _contract().validate_inputs(_inputs_with(frames, df))
    assert _features_issues(errors(issues)) == []


@pytest.mark.parametrize("topology", ["circular", "linear"])
@pytest.mark.parametrize("delta", [0, 1, 1000])
def test_c15_start_at_or_beyond_length_rejected(
    frames: dict[str, Any], topology: str, delta: int
) -> None:
    """Lead decision 23: ``start >= length`` is an error on every contig, circular or linear,
    with code ``contract.features.start_bounds`` (``start`` is never unwrapped)."""
    n = 4000
    frames["replicons"] = pd.DataFrame(
        {"contig": ["c1"], "length": pd.Series([n], dtype="int64"), "topology": [topology],
         "kind": ["chromosome" if topology == "circular" else "unknown"]}
    )  # fmt: skip
    frames["genes"] = empty_like(frames["genes"])
    df = _feature_rows(frames, contig="c1", start=n + delta, end=n + delta + 232)
    issues = _contract().validate_table("features", df)
    issues += _contract().validate_inputs(_inputs_with(frames, df))
    codes = [i.code for i in errors(issues)]
    assert "contract.features.start_bounds" in codes, codes


@pytest.mark.parametrize("topology", ["circular", "linear"])
def test_c15_start_just_below_length_valid(frames: dict[str, Any], topology: str) -> None:
    n = 4000
    frames["replicons"] = pd.DataFrame(
        {"contig": ["c1"], "length": pd.Series([n], dtype="int64"), "topology": [topology],
         "kind": ["chromosome" if topology == "circular" else "unknown"]}
    )  # fmt: skip
    frames["genes"] = empty_like(frames["genes"])
    df = _feature_rows(frames, contig="c1", start=n - 1, end=n)
    issues = _contract().validate_table("features", df)
    issues += _contract().validate_inputs(_inputs_with(frames, df))
    assert _features_issues(errors(issues)) == []


@pytest.mark.parametrize("end_delta", [1, 100, 4000])
def test_c15_end_beyond_length_on_a_linear_contig_rejected(
    frames: dict[str, Any], end_delta: int
) -> None:
    """Lead decision 23: nothing crosses the origin of a linear contig, so ``end > length`` there
    is ``contract.features.bounds`` (on a circular replicon the same row is valid)."""
    n = 4000
    base = pd.DataFrame(
        {"contig": ["lin", "circ"], "length": pd.Series([n, n], dtype="int64"),
         "topology": ["linear", "circular"], "kind": ["unknown", "chromosome"]}
    )  # fmt: skip
    frames["replicons"] = base
    frames["genes"] = empty_like(frames["genes"])
    bad = _feature_rows(frames, contig="lin", start=n - 10, end=n + end_delta)
    issues = _contract().validate_inputs(_inputs_with(frames, bad))
    assert "contract.features.bounds" in [i.code for i in errors(issues)]
    ok = _feature_rows(frames, contig="circ", start=n - 10, end=n + end_delta)
    assert _features_issues(errors(_contract().validate_inputs(_inputs_with(frames, ok)))) == []


def test_c15_end_equal_length_on_a_linear_contig_valid(frames: dict[str, Any]) -> None:
    n = 4000
    frames["replicons"] = pd.DataFrame(
        {"contig": ["lin"], "length": pd.Series([n], dtype="int64"), "topology": ["linear"],
         "kind": ["unknown"]}
    )  # fmt: skip
    frames["genes"] = empty_like(frames["genes"])
    df = _feature_rows(frames, contig="lin", start=n - 232, end=n)
    assert _features_issues(errors(_contract().validate_inputs(_inputs_with(frames, df)))) == []


def test_c15_validate_inputs_does_not_mutate(frames: dict[str, Any]) -> None:
    df = frames["features"].copy()
    _contract().validate_inputs(_inputs_with(frames, df))
    pd.testing.assert_frame_equal(df, frames["features"])


def _bakta_with_oric(tmp_path: Path) -> tuple[Path, Path]:
    """A minimal Bakta JSON whose oriC has Bakta's strand ``?`` (as K-12 and Nissle 1917 do)."""
    import json

    rng = np.random.default_rng(15)
    chrom = "".join(rng.choice(list("ACGT"), 3000))
    data = {
        "version": {"bakta": "1.12.2", "db": {"version": "6.0", "type": "light"}},
        "sequences": [
            {"id": "chr1", "type": "chromosome", "topology": "circular", "length": 3000,
             "nt": chrom},
        ],
        "features": [
            {"type": "cds", "sequence": "chr1", "start": 101, "stop": 400, "strand": "+",
             "locus": "T_00010", "aa": "MKV", "product": "protein"},
            {"type": "oriC", "sequence": "chr1", "start": 1001, "stop": 1232, "strand": "?",
             "id": "T_oriC"},
            {"type": "tRNA", "sequence": "chr1", "start": 2101, "stop": 2176, "strand": "-",
             "id": "T_t1", "product": "tRNA-Leu", "amino_acid": "Leu", "anti_codon": "cag"},
        ],
    }  # fmt: skip
    tmp_path.mkdir(parents=True, exist_ok=True)
    js, fa = tmp_path / "bakta.json", tmp_path / "genome.norm.fa"
    js.write_text(json.dumps(data), encoding="utf-8")
    fa.write_text(f">chr1\n{chrom}\n", encoding="utf-8")
    return js, fa


def test_c15_parse_bakta_output_conforms(tmp_path: Path) -> None:
    """Lead decision 11: Bakta's oriC strand ``?`` arrives as ``.`` through ``parse_bakta``, and
    its ``features.parquet`` passes C15 (table and cross-table with its ``replicons.parquet``)."""
    import importlib

    from giap.config import GiapConfig
    from giap.steps import registry

    mod = importlib.import_module("giap.steps.annotation")
    if "parse_bakta" not in registry._REGISTRY:
        importlib.reload(mod)
    js, fa = _bakta_with_oric(tmp_path / "in")
    names = {
        "genes": "genes.parquet", "features": "features.parquet",
        "replicons": "replicons.parquet", "rrna16s": "rrna16s.fa",
        "trna_counts": "trna_counts.tsv", "xrefs": "xrefs.parquet",
        "annotation_meta": "annotation_meta.json",
    }  # fmt: skip
    outs = {k: tmp_path / "work" / v for k, v in names.items()}
    registry.run_step("parse_bakta", GiapConfig(seed=1), {"bakta_json": js, "genome": fa}, outs)
    feats = pd.read_parquet(outs["features"])
    oric = feats[feats["type"] == "oric"]
    assert len(oric) == 1
    assert oric["strand"].tolist() == ["."]
    assert (int(oric["start"].iloc[0]), int(oric["end"].iloc[0])) == (1000, 1232)
    assert _contract().validate_table("features", feats) == []
    frames = build_frames(5)
    frames["replicons"] = pd.read_parquet(outs["replicons"])
    issues = _contract().validate_inputs(_inputs_with(frames, feats))
    assert _features_issues(issues) == []


def test_c15_giapsim_literal_copy_equals_tables() -> None:
    """Acceptance: giapsim's literal C15 column list equals ``TABLES`` (drift check)."""
    from giapsim.contract_tables import CONTRACT_COLUMNS

    assert "features" in CONTRACT_COLUMNS
    assert tuple(CONTRACT_COLUMNS["features"]) == tuple(
        c.name for c in _contract().TABLES["features"]
    )
