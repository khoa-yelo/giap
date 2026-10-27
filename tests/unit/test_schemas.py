"""Bundled JSON Schemas (GIPPA_SPEC.md §6.7) and the committed config schema."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import jsonschema
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCHEMA_DIR = ROOT / "src" / "gippa" / "schemas"
SCHEMAS = ["run_manifest.schema.json", "config.schema.json", "preflight.schema.json"]


def _load(name: str) -> dict[str, Any]:
    data = json.loads((SCHEMA_DIR / name).read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


@pytest.mark.parametrize("name", SCHEMAS)
def test_schema_is_valid_json_object(name: str) -> None:
    schema = _load(name)
    assert schema.get("type") == "object"
    assert "properties" in schema


@pytest.mark.parametrize("name", SCHEMAS)
def test_schema_is_valid_against_its_metaschema(name: str) -> None:
    schema = _load(name)
    # Draft 2020-12 unless the schema declares another draft.
    validator = jsonschema.validators.validator_for(schema, default=jsonschema.Draft202012Validator)
    validator.check_schema(schema)


def test_run_manifest_declares_draft_2020_12() -> None:
    assert _load("run_manifest.schema.json")["$schema"] == (
        "https://json-schema.org/draft/2020-12/schema"
    )


def test_config_schema_validates_default_config() -> None:
    from gippa.config import resolve_config

    schema = _load("config.schema.json")
    instance = resolve_config().model_dump(mode="json")
    jsonschema.Draft202012Validator(schema).validate(instance)


def test_committed_config_schema_is_up_to_date() -> None:
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "write_config_schema.py"), "--check"],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=False,
    )
    assert result.returncode == 0, result.stderr


# --- M3-16: Frictionless table schemas of the Module 2 outputs (§9.3, §9.5) ---------------------

TABLE_DIR = SCHEMA_DIR / "tables"
MAPPING_TABLES = ("peptide_hits", "protein_hits", "unmapped", "mapping_fdr")


def _table_schema(name: str) -> dict[str, Any]:
    data = json.loads((TABLE_DIR / f"{name}.schema.json").read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def _fields(name: str) -> dict[str, dict[str, Any]]:
    return {f["name"]: f for f in _table_schema(name)["fields"]}


def _enum(name: str, column: str) -> set[str]:
    return set(_fields(name)[column]["constraints"]["enum"])


def _columns(name: str) -> tuple[str, ...]:
    from gippa.mapping import outputs

    return getattr(outputs, f"{name.upper()}_COLUMNS")  # type: ignore[no-any-return]


@pytest.mark.parametrize("name", MAPPING_TABLES)
def test_mapping_table_schema_fields_match_columns(name: str) -> None:
    assert tuple(f["name"] for f in _table_schema(name)["fields"]) == _columns(name)


@pytest.mark.parametrize("name", MAPPING_TABLES)
def test_mapping_table_schema_is_valid_frictionless(name: str) -> None:
    frictionless = pytest.importorskip("frictionless")
    descriptor = {k: v for k, v in _table_schema(name).items() if k != "coordinates"}
    frictionless.Schema.from_descriptor(descriptor)  # raises on an invalid descriptor


@pytest.mark.parametrize("name", MAPPING_TABLES)
def test_mapping_table_schema_na_is_missing(name: str) -> None:
    assert "NA" in _table_schema(name)["missingValues"]


@pytest.mark.parametrize("name", ["peptide_hits", "protein_hits"])
def test_mapping_table_schema_coordinates(name: str) -> None:
    assert _table_schema(name)["coordinates"] == "1-based-closed"


def test_peptide_hits_schema_enums_and_types() -> None:
    from gippa.mapping.loci import LOCUS_CLASSES

    f = _fields("peptide_hits")
    assert _enum("peptide_hits", "route") == {"exact", "tolerant"}
    assert _enum("peptide_hits", "locus_class") == set(LOCUS_CLASSES)
    assert "cterm_extension" in _enum("peptide_hits", "locus_class")  # ADR-0013
    assert _enum("peptide_hits", "uniqueness") == {"unique", "multi_identical", "shared"}
    assert _enum("peptide_hits", "strand") >= {"+", "-"}
    for col in ("peptide_len", "start", "end", "frame", "n_diff"):
        assert f[col]["type"] == "integer", col
    for col in ("best_q_input", "pep_input", "score_input", "placement_q"):
        assert f[col]["type"] == "number", col
    assert f["crosses_origin"]["type"] == "boolean"


def test_protein_hits_schema_enums_and_types() -> None:
    f = _fields("protein_hits")
    assert _enum("protein_hits", "route") == {"miniprot", "tblastn"}
    assert _enum("protein_hits", "class") == {"identical", "near", "divergent"}
    for col in ("protein_len_aa", "start", "end", "rank"):
        assert f[col]["type"] == "integer", col
    for col in ("identity", "query_coverage", "reciprocal_overlap"):
        assert f[col]["type"] == "number", col


def test_unmapped_schema_enums() -> None:
    assert _enum("unmapped", "kind") == {"peptide", "protein"}
    assert _enum("unmapped", "label") == {"contaminant", "unmapped"}
    assert _enum("unmapped", "reason") == {
        "no_hit",
        "rejected_placement_fdr",
        "below_identity",
        "too_short",
        "contaminant_id",
    }


def test_mapping_fdr_schema_enums_and_types() -> None:
    f = _fields("mapping_fdr")
    assert _enum("mapping_fdr", "decoy_type") == {"reverse", "shuffle"}
    assert _enum("mapping_fdr", "section") == {"placement", "denovo_precision"}
    for col in ("n_diff", "targets", "decoys", "n_accepted", "score_bin"):
        assert f[col]["type"] == "integer", col
    for col in ("pi_wrong", "fdr_stratum", "precision", "ci_low", "ci_high", "c_value"):
        assert f[col]["type"] == "number", col


def _rich_outputs(outdir: Path) -> None:
    """Write outputs that use every enum value of §9.3 (stand-in for the 50-kb golden fixture)."""
    import pandas as pd

    from gippa.mapping.contaminants import UnmappedRow
    from gippa.mapping.loci import LOCUS_CLASSES
    from gippa.mapping.outputs import PEPTIDE_HITS_COLUMNS, PROTEIN_HITS_COLUMNS
    from gippa.model import Replicon

    reps = {"chr1": Replicon("chr1", 50_000, "circular"), "p1": Replicon("p1", 3000, "linear")}
    pep_rows = []
    uniq = ("unique", "multi_identical", "shared")
    for i, lc in enumerate(LOCUS_CLASSES):
        exact = i % 2 == 0
        start = 49_990 if i == 0 else 300 * i  # i == 0 crosses the origin
        pep_rows.append(
            {
                "peptide": f"PEPTIDE{'K' * (i + 1)}",
                "peptide_len": 8 + i,
                "source": "database" if exact else "denovo",
                "best_q_input": 0.001 * i if exact else None,
                "pep_input": 0.01,
                "score_input": 10.0 + i,
                "contig": "chr1",
                "start": start,
                "end": start + 3 * (8 + i),
                "strand": "+" if i % 3 else "-",
                "frame": 1 + start % 3 if i % 3 else -1,
                "route": "exact" if exact else "tolerant",
                "n_diff": 0 if exact else 1,
                "diff_types": None if exact else ("substitution",),
                "placement_q": 0.0001 * (i + 1),
                "length_class": "7-9" if i < 2 else "10-11",
                "locus_class": lc,
                "uniqueness": uniq[i % 3],
                "gene": None if lc in ("intergenic", "antisense") else f"G_{i:05d}",
                "novel_candidate_id": None,
                "flags": "same_db" if i == 1 else "",
                "blocks": None,
            }
        )
    peps = pd.DataFrame(pep_rows)
    assert set(PEPTIDE_HITS_COLUMNS) - {"crosses_origin"} <= set(peps.columns)
    prot_rows = [
        {
            "protein_id": f"P{i}",
            "protein_len_aa": 100,
            "route": route,
            "contig": "p1",
            "start": 100 * i,
            "end": 100 * i + 300,
            "strand": "+",
            "rank": 1,
            "identity": ident,
            "query_coverage": 0.9,
            "reciprocal_overlap": None if cls == "divergent" else 0.9,
            "class": cls,
            "gene": None if cls == "divergent" else f"G_P{i}",
            "flags": "fs" if i == 1 else "",
        }
        for i, (route, cls, ident) in enumerate(
            [
                ("miniprot", "identical", 1.0),
                ("tblastn", "near", 0.95),
                ("miniprot", "divergent", 0.5),
            ]
        )
    ]
    prots = pd.DataFrame(prot_rows, columns=list(PROTEIN_HITS_COLUMNS))
    unmapped = [
        UnmappedRow("u1", "AAAK", "peptide", "unmapped", "no_hit", None),
        UnmappedRow("u2", "AAAAK", "peptide", "unmapped", "rejected_placement_fdr", None),
        UnmappedRow("u3", "AK", "peptide", "unmapped", "too_short", None),
        UnmappedRow("u4", "MKKL", "protein", "unmapped", "below_identity", None),
        UnmappedRow("CON_1", "MSTV", "protein", "contaminant", "contaminant_id", "CON_1"),
    ]
    placement = pd.DataFrame(
        [
            {
                "section": "placement",
                "decoy_type": dt,
                "n_diff": m,
                "length_class": lc,
                "targets": 10,
                "decoys": 1,
                "pi_wrong": 0.4,
                "fdr_stratum": 0.05,
                "n_accepted": 8,
            }
            for dt in ("reverse", "shuffle")
            for m in (0, 1, 2)
            for lc in ("7-9", "10-11", "all")
        ]
    )
    denovo = pd.DataFrame(
        [
            {"section": "denovo_precision", "score_bin": 0, "precision": None, "ci_low": None,
             "ci_high": None, "c_value": 0.99, "flags": "unreported"},
            {"section": "denovo_precision", "score_bin": 1, "precision": 0.8, "ci_low": 0.7,
             "ci_high": 0.9, "c_value": 0.99, "flags": ""},
        ]
    )  # fmt: skip
    fdr = pd.concat([placement, denovo], ignore_index=True)
    from gippa.mapping.outputs import write_mapping_outputs

    write_mapping_outputs(outdir, peps, prots, unmapped, fdr, reps, "cd" * 32)


def test_example_09_frictionless_validate_mapping_tables(tmp_path: Path) -> None:
    frictionless = pytest.importorskip("frictionless")
    out = tmp_path / "mapping"
    out.mkdir()
    _rich_outputs(out)
    for name in MAPPING_TABLES:
        descriptor = {k: v for k, v in _table_schema(name).items() if k != "coordinates"}
        resource = frictionless.Resource(
            path=f"{name}.tsv",
            basepath=str(out),
            schema=frictionless.Schema.from_descriptor(descriptor),
        )
        report = resource.validate()
        assert report.valid, (name, report.flatten(["rowNumber", "fieldName", "type", "note"]))


def test_example_09_frictionless_cli_accepts_schema_files(tmp_path: Path) -> None:
    exe = Path(sys.executable).with_name("frictionless")
    if not exe.exists():
        pytest.skip("frictionless CLI not installed next to the interpreter")
    out = tmp_path / "mapping"
    out.mkdir()
    _rich_outputs(out)
    for name in MAPPING_TABLES:
        result = subprocess.run(
            [str(exe), "validate", "--trusted", "--schema", str(TABLE_DIR / f"{name}.schema.json"),
             f"{name}.tsv"],
            capture_output=True, text=True, cwd=out, check=False,
        )  # fmt: skip
        assert result.returncode == 0, (name, result.stdout[-2000:], result.stderr[-2000:])
