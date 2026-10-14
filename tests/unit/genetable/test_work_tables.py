"""M4-06 loader: ``giap.io.work_tables`` (docs/specs/M4-06-input-contract.md, "Public interface").

``load_m4_inputs`` reads M2's Parquet intermediates (0-based) and M3's final ``mapping/*.tsv``
(1-based closed, ``;`` lists, ``NA``/``INF`` per ADR-0006), converts on read and never raises on
content: it returns issues.
"""

from __future__ import annotations

import dataclasses
import json
import math
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from m4_contract_fixture import FILE_NAMES, build_frames, empty_like, write_work

from giap.io.issues import Issue

KEYS = set(FILE_NAMES)
OPTIONAL_FILES = ("abundance", "id_map")
REQUIRED_FILES = (
    "replicons", "genes", "xrefs", "annotation_meta", "tus", "layers", "regions",
    "region_genes", "peptide_hits", "protein_hits", "unmapped",
)  # fmt: skip


def _wt() -> Any:
    from giap.io import work_tables

    return work_tables


def errors(issues: list[Issue]) -> list[Issue]:
    return [i for i in issues if i.severity == "error"]


def _load(tmp_path: Path, frames: dict[str, Any] | None = None, skip: tuple[str, ...] = ()) -> Any:
    wt = _wt()
    f = frames if frames is not None else build_frames(20)
    paths = wt.default_paths(tmp_path / "work")
    write_work(f, paths, skip=skip)
    return wt.load_m4_inputs(paths)


def test_input_keys() -> None:
    assert set(_wt().M4_INPUT_KEYS) == KEYS


def test_default_paths_layout(tmp_path: Path) -> None:
    work = tmp_path / "work"
    paths = _wt().default_paths(work)
    assert set(paths) == KEYS
    for key, name in FILE_NAMES.items():
        assert paths[key].name == name, key
        assert work in paths[key].parents, key


def test_m4inputs_fields() -> None:
    names = {f.name for f in dataclasses.fields(_wt().M4Inputs)}
    assert names == {
        "replicons", "genes", "xrefs", "annotation_source", "tus", "layers", "regions",
        "region_genes", "peptides", "abundance", "peptide_hits", "protein_hits", "unmapped",
        "id_map",
    }  # fmt: skip


def test_roundtrip_no_errors(tmp_path: Path) -> None:
    inputs, issues = _load(tmp_path)
    assert errors(issues) == [], issues
    assert len(inputs.genes) == 20
    assert inputs.annotation_source == "bakta_full"
    assert list(inputs.genes["locus_tag"]) == list(build_frames(20)["genes"]["locus_tag"])


def test_tsv_coordinates_converted_to_zero_based(tmp_path: Path) -> None:
    """M3 TSVs are 1-based closed; on read ``start - 1``, ``end`` unchanged (spec edge cases)."""
    f = build_frames(20)
    inputs, _ = _load(tmp_path, f)
    exp = f["peptide_hits"].sort_values(["start", "peptide"]).reset_index(drop=True)
    got = inputs.peptide_hits.sort_values(["start", "peptide"]).reset_index(drop=True)
    assert [int(x) for x in got["start"]] == [int(x) for x in exp["start"]]
    assert [int(x) for x in got["end"]] == [int(x) for x in exp["end"]]
    assert int(inputs.protein_hits["start"].iloc[0]) == int(f["protein_hits"]["start"].iloc[0])
    assert int(inputs.protein_hits["end"].iloc[0]) == int(f["protein_hits"]["end"].iloc[0])


def test_tsv_raw_file_is_one_based(tmp_path: Path) -> None:
    """Guard on the fixture: the file on disk holds ``start + 1`` (so the test above is real)."""
    f = build_frames(20)
    paths = _wt().default_paths(tmp_path / "work")
    write_work(f, paths)
    raw = pd.read_csv(paths["peptide_hits"], sep="\t")
    assert sorted(raw["start"]) == sorted(int(s) + 1 for s in f["peptide_hits"]["start"])


def test_tsv_lists_split_and_values_typed(tmp_path: Path) -> None:
    f = build_frames(20)
    h = f["peptide_hits"]
    h.at[0, "diff_types"] = ["substitution", "swap"]
    h.at[0, "flags"] = ["same_db", "strain_difference"]
    h.at[0, "route"] = "tolerant"
    h.at[0, "n_diff"] = 2
    h.loc[0, "score_input"] = math.inf
    h.loc[0, "pep_input"] = -math.inf
    h.loc[0, "best_q_input"] = math.nan
    inputs, issues = _load(tmp_path, f)
    assert errors(issues) == []
    row = inputs.peptide_hits[inputs.peptide_hits["peptide"] == h.at[0, "peptide"]].iloc[0]
    assert list(row["diff_types"]) == ["substitution", "swap"]
    assert list(row["flags"]) == ["same_db", "strain_difference"]
    assert row["score_input"] == math.inf  # ADR-0006: INF -> +inf
    assert row["pep_input"] == -math.inf  # -INF -> -inf
    assert pd.isna(row["best_q_input"])  # NA -> missing
    assert int(row["n_diff"]) == 2
    assert bool(row["crosses_origin"]) is False
    assert pd.api.types.is_bool_dtype(inputs.peptide_hits["crosses_origin"])
    assert pd.api.types.is_integer_dtype(inputs.peptide_hits["start"])
    assert pd.api.types.is_float_dtype(inputs.peptide_hits["placement_q"])


def test_tsv_na_gene_is_missing(tmp_path: Path) -> None:
    inputs, _ = _load(tmp_path)
    h = inputs.peptide_hits
    inter = h[h["locus_class"] == "intergenic"]
    assert len(inter) == 1 and pd.isna(inter["gene"].iloc[0])
    assert pd.isna(inputs.unmapped["contaminant_match"].iloc[0])


def test_header_only_tsvs(tmp_path: Path) -> None:
    """Header-only C11 (no peptide placed) and C12 (no proteins) load empty without errors."""
    f = build_frames(20)
    f["peptide_hits"] = empty_like(f["peptide_hits"])
    f["protein_hits"] = empty_like(f["protein_hits"])
    inputs, issues = _load(tmp_path, f)
    assert errors(issues) == []
    assert len(inputs.peptide_hits) == 0 and len(inputs.protein_hits) == 0
    assert "placement_q" in inputs.peptide_hits.columns
    assert "rank" in inputs.protein_hits.columns


@pytest.mark.parametrize("key", OPTIONAL_FILES)
def test_missing_optional_file_warns(tmp_path: Path, key: str) -> None:
    inputs, issues = _load(tmp_path, skip=(key,))
    assert errors(issues) == []
    assert getattr(inputs, key) is None
    assert any(i.severity == "warn" and i.code.startswith("contract.") for i in issues)


@pytest.mark.parametrize("key", REQUIRED_FILES)
def test_missing_required_file_is_error_not_raise(tmp_path: Path, key: str) -> None:
    _, issues = _load(tmp_path, skip=(key,))
    errs = errors(issues)
    assert errs, f"missing {key} gave no error"
    assert all(i.code.startswith("contract.") for i in errs)
    assert any(FILE_NAMES[key] in i.message or key in i.code for i in errs)


def test_annotation_source_from_meta(tmp_path: Path) -> None:
    f = build_frames(20)
    f["annotation_source"] = "user"
    inputs, issues = _load(tmp_path, f)
    assert inputs.annotation_source == "user"
    assert errors(issues) == []


def test_meta_without_source_is_error(tmp_path: Path) -> None:
    wt = _wt()
    paths = wt.default_paths(tmp_path / "work")
    write_work(build_frames(20), paths)
    paths["annotation_meta"].write_text(json.dumps({"flags": []}), encoding="utf-8")
    _, issues = wt.load_m4_inputs(paths)
    assert errors(issues)


def test_bad_content_returns_issues_not_raise(tmp_path: Path) -> None:
    f = build_frames(20)
    f["genes"].loc[0, "feature_type"] = "CDS"
    f["genes"].loc[1, "locus_tag"] = f["genes"].loc[2, "locus_tag"]
    _, issues = _load(tmp_path, f)
    codes = {i.code for i in errors(issues)}
    assert "contract.genes.enum" in codes
    assert "contract.genes.locus_tag_unique" in codes


def test_extra_tsv_column_is_warning(tmp_path: Path) -> None:
    f = build_frames(20)
    f["peptide_hits"]["my_extra"] = "x"
    _, issues = _load(tmp_path, f)
    assert errors(issues) == []
    assert any(i.severity == "warn" and "my_extra" in i.message for i in issues)


def test_missing_tsv_column_is_error(tmp_path: Path) -> None:
    f = build_frames(20)
    f["peptide_hits"] = f["peptide_hits"].drop(columns=["placement_q"])
    _, issues = _load(tmp_path, f)
    assert any(
        i.code.startswith("contract.peptide_hits.") and "placement_q" in i.message
        for i in errors(issues)
    )


def test_parquet_list_columns_load_as_lists(tmp_path: Path) -> None:
    inputs, _ = _load(tmp_path)
    tu_genes = inputs.tus["genes"].iloc[1]
    assert [str(x) for x in tu_genes] == list(build_frames(20)["tus"]["genes"].iloc[1])


def test_id_map_absent_leaves_none(tmp_path: Path) -> None:
    """Edge case: C14 absent → ``id_map`` None (``source_id`` NA downstream)."""
    inputs, _ = _load(tmp_path, skip=("id_map",))
    assert inputs.id_map is None
