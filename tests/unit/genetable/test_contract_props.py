"""Hypothesis properties of the M4-06 validators (TU partition, column checks, no raising)."""

from __future__ import annotations

from typing import Any

import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from m4_contract_fixture import build_frames, make_inputs, tag

TABLE_KEYS = (
    "replicons", "genes", "xrefs", "tus", "layers", "regions", "region_genes", "peptides",
    "abundance", "peptide_hits", "protein_hits", "unmapped", "id_map",
)  # fmt: skip
ORDER_FREE = ("xrefs", "region_genes", "peptides", "abundance", "unmapped", "id_map")
SETTINGS = settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])


def _validate(frames: dict[str, Any]) -> list[Any]:
    from giap.genetable.contract import validate_inputs, validate_table

    out = []
    for k in TABLE_KEYS:
        if frames.get(k) is not None:
            out += validate_table(k, frames[k])
    return out + validate_inputs(make_inputs(frames))


layouts = st.tuples(
    st.integers(min_value=1, max_value=60),
    st.lists(st.integers(min_value=1, max_value=6), min_size=1, max_size=5),
)


@SETTINGS
@given(layouts)
def test_valid_layouts_have_no_errors(layout: tuple[int, list[int]]) -> None:
    n, sizes = layout
    issues = _validate(build_frames(n, sizes))
    assert [i for i in issues if i.severity == "error"] == []


@SETTINGS
@given(layouts, st.data())
def test_dropping_a_gene_from_tus_is_reported(layout: tuple[int, list[int]], data: Any) -> None:
    n, sizes = layout
    f = build_frames(n, sizes)
    victim = tag(data.draw(st.integers(min_value=0, max_value=n - 1)))
    tus = f["tus"].copy()
    tus["genes"] = [[g for g in gs if g != victim] for gs in tus["genes"]]
    f["tus"] = tus[[len(g) > 0 for g in tus["genes"]]].reset_index(drop=True)
    from giap.genetable.contract import validate_inputs

    hit = [i for i in validate_inputs(make_inputs(f)) if i.code == "contract.tus.gene_coverage"]
    assert hit and any(victim in i.message + " ".join(i.examples) for i in hit)


@SETTINGS
@given(layouts, st.data())
def test_gene_in_two_tus_is_reported(layout: tuple[int, list[int]], data: Any) -> None:
    n, sizes = layout
    f = build_frames(n, sizes)
    tus = f["tus"].copy()
    if len(tus) < 2:
        tus = pd.concat([tus, tus.iloc[[0]].assign(tu_id="TU_999")], ignore_index=True)
    i = data.draw(st.integers(min_value=0, max_value=len(tus) - 1))
    j = data.draw(st.integers(min_value=0, max_value=len(tus) - 1).filter(lambda x: x != i))
    g = tus.at[i, "genes"][0]
    tus.at[j, "genes"] = [*tus.at[j, "genes"], g]
    f["tus"] = tus
    from giap.genetable.contract import validate_inputs

    hit = [i for i in validate_inputs(make_inputs(f)) if i.code == "contract.tus.gene_coverage"]
    assert hit and all(x.severity == "error" for x in hit)


@SETTINGS
@given(st.sampled_from(TABLE_KEYS), st.data())
def test_any_column_subset_never_raises_and_reports_each_missing(table: str, data: Any) -> None:
    from giap.genetable.contract import TABLES, validate_table

    df = build_frames(20)[table]
    cols = list(df.columns)
    drop = data.draw(st.lists(st.sampled_from(cols), unique=True, max_size=len(cols)))
    issues = validate_table(table, df.drop(columns=drop))
    required = {c.name for c in TABLES[table] if not c.optional}
    errs = [i for i in issues if i.severity == "error"]
    for c in set(drop) & required:
        assert any(c in i.message + " ".join(i.examples) for i in errs), c


@SETTINGS
@given(st.sampled_from(ORDER_FREE), st.integers(min_value=0, max_value=10_000))
def test_row_permutation_does_not_change_issue_codes(table: str, seed: int) -> None:
    """Per-table checks are row-order free (codes and severities, not row numbers)."""
    from giap.genetable.contract import validate_table

    df = build_frames(20)[table]
    perm = df.sample(frac=1.0, random_state=seed).reset_index(drop=True)
    a = sorted((i.code, i.severity) for i in validate_table(table, df))
    b = sorted((i.code, i.severity) for i in validate_table(table, perm))
    assert a == b


# C9 ``protein_ids`` (lead decision 1): native list of non-empty, non-NA strings in Parquet.
_ids = st.text(
    alphabet=st.characters(min_codepoint=33, max_codepoint=126, exclude_characters=";"),
    min_size=1,
    max_size=12,
)
_id_lists = st.lists(_ids, max_size=5)


def _peptide_issues(protein_ids: list[Any]) -> list[Any]:
    from giap.genetable.contract import validate_table

    p = build_frames(20)["peptides"]
    p["protein_ids"] = pd.Series(protein_ids, dtype=object)
    return [
        i
        for i in validate_table("peptides", p)
        if i.severity == "error" and "protein_ids" in i.message + " ".join(i.examples)
    ]


@SETTINGS
@given(st.lists(_id_lists, min_size=4, max_size=4), st.booleans())
def test_protein_ids_lists_of_nonempty_strings_valid(lists: list[list[str]], arrays: bool) -> None:
    import numpy as np

    cells: list[Any] = [np.array(x, dtype=object) if arrays else x for x in lists]
    assert _peptide_issues(cells) == []


@SETTINGS
@given(
    st.lists(_id_lists, min_size=4, max_size=4),
    st.integers(min_value=0, max_value=3),
    st.sampled_from(["", None, float("nan"), pd.NA]),
    st.data(),
)
def test_protein_ids_bad_element_anywhere_is_error(
    lists: list[list[str]], row: int, bad: Any, data: Any
) -> None:
    pos = data.draw(st.integers(min_value=0, max_value=len(lists[row])))
    lists[row] = [*lists[row][:pos], bad, *lists[row][pos:]]
    assert _peptide_issues(list(lists))


@SETTINGS
@given(st.lists(_id_lists, min_size=4, max_size=4), st.integers(min_value=0, max_value=3))
def test_protein_ids_joined_string_cell_is_dtype_error(lists: list[list[str]], row: int) -> None:
    from giap.genetable.contract import validate_table

    cells: list[Any] = list(lists)
    cells[row] = ";".join(lists[row]) or "GSIM_00005"
    p = build_frames(20)["peptides"]
    p["protein_ids"] = pd.Series(cells, dtype=object)
    hit = [i for i in validate_table("peptides", p) if i.code == "contract.peptides.dtype"]
    assert hit and all(i.severity == "error" for i in hit)
