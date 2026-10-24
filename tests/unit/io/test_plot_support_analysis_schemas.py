"""M9-17 support schemas and descriptor (docs/specs/M9-17-plot-support-analysis.md "Support
files"; ADR-0043 rules 3 and 6).

Each of the six files has ``src/giap/schemas/plot_support/<name>.json`` with the spec's columns in
order and ``NA`` as a missing value, and one ``datapackage.json`` resource with ``role:
plot_support`` and a ``path`` under ``plots/support/``. Contract: tests/helpers/
support_analysis_data.py.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from helpers import support_analysis_data as sa

import giap.schemas

SCHEMAS = Path(giap.schemas.__file__).resolve().parent
NAMES = tuple(sa.COLUMNS)


def _schema(name: str) -> dict[str, Any]:
    path = SCHEMAS / "plot_support" / f"{name}.json"
    assert path.is_file(), f"missing schema {path}"
    data: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return data


def _fields(name: str) -> dict[str, dict[str, Any]]:
    return {f["name"]: f for f in _schema(name)["fields"]}


@pytest.mark.parametrize("name", NAMES)
def test_schema_fields_are_the_spec_columns(name: str) -> None:
    schema = _schema(name)
    assert [f["name"] for f in schema["fields"]] == list(sa.COLUMNS[name])
    assert "NA" in schema.get("missingValues", [])
    assert schema.get("name") == name


@pytest.mark.parametrize("name", NAMES)
def test_datapackage_resource_has_role_plot_support(name: str) -> None:
    pkg = json.loads((SCHEMAS / "datapackage.json").read_text(encoding="utf-8"))
    res = [r for r in pkg["resources"] if r["name"] == name]
    assert len(res) == 1, [r["name"] for r in pkg["resources"]]
    r = res[0]
    assert r["role"] == "plot_support"
    assert r["path"].startswith(f"plots/support/{name}.tsv")
    assert r.get("profile") == "tabular-data-resource"
    assert (SCHEMAS / r["schema"]).resolve() == (SCHEMAS / "plot_support" / f"{name}.json")


def test_field_types() -> None:
    """Integers for positions, replicates and permutations; numbers for every value; ``in_set``
    a boolean."""
    want = {
        "gsea_running": {"rank_position": "integer", "score": "number", "in_set": "boolean",
                         "running_es": "number"},
        "determinants_shap": {"feature_value": "number", "shap": "number"},
        "determinants_pd": {"grid_value": "number", "pd_mean": "number", "pd_low": "number",
                            "pd_high": "number"},
        "null_regions": {"shift_k": "integer", "null_value": "number", "observed": "number"},
        "null_islands": {"replicate": "integer", "max_scan_s": "number",
                         "observed_max": "number"},
        "null_determinants": {"permutation": "integer", "cv_r2": "number",
                              "observed_cv_r2_mean": "number"},
    }  # fmt: skip
    for name, types in want.items():
        fields = _fields(name)
        for col, t in types.items():
            assert fields[col]["type"] == t, (name, col, fields[col]["type"])


def test_enums() -> None:
    assert set(_fields("determinants_pd")["feature"]["constraints"]["enum"]) == set(sa.PD_FEATURES)
    assert set(_fields("null_determinants")["model"]["constraints"]["enum"]) == {
        "lightgbm", "elastic_net",
    }  # fmt: skip
    assert set(_fields("null_regions")["stat"]["constraints"]["enum"]) == {"o_minus_e", "hl_shift"}


def test_null_regions_schema_says_which_statistic_and_the_grid_header() -> None:
    """The null stores the statistic each test compares (ADR-0018 studentised O - E; the
    amendments' studentised HL), and Q2's grid header is documented."""
    schema = _schema("null_regions")
    text = json.dumps(schema).lower()
    assert "studenti" in text
    assert "percentile" in text, "Q2: the schema must describe the percentile-grid header"


def test_pd_schema_says_spread_not_ci() -> None:
    """pd_low / pd_high are the min / max over boosters: a spread, labelled so, not a CI."""
    fields = _fields("determinants_pd")
    for col in ("pd_low", "pd_high"):
        desc = fields[col].get("description", "").lower()
        assert "spread" in desc and "booster" in desc, (col, desc)
        assert "confidence" not in desc or "not a confidence" in desc, (col, desc)


# ----------------------------------------------------------------------------- review round
# Lead decision 8 (PR #272): the grid form's own schema and the wording fixes.


def test_grid_schema_fields_are_the_grid_columns() -> None:
    """``plot_support/null_regions_grid.json``: ``sa.GRID_COLUMNS`` in order, ``NA`` missing,
    ``percentile`` an integer 1..99 whose description says it is the 1-based rank for a null of
    at most 99 finite values (statistician N1)."""
    schema = _schema("null_regions_grid")
    assert [f["name"] for f in schema["fields"]] == list(sa.GRID_COLUMNS)
    assert "NA" in schema.get("missingValues", [])
    assert schema.get("name") == "null_regions_grid"
    f = _fields("null_regions_grid")
    assert f["percentile"]["type"] == "integer"
    c = f["percentile"].get("constraints", {})
    assert c.get("minimum") == 1 and c.get("maximum") == 99, c
    assert "rank" in f["percentile"].get("description", "").lower()
    for col in ("null_value", "observed"):
        assert f[col]["type"] == "number", col
    assert set(f["stat"]["constraints"]["enum"]) == {"o_minus_e", "hl_shift"}


def test_null_regions_schema_names_the_fourth_column() -> None:
    """Review N1: the grid's ``percentile`` replaces ``shift_k``, the fourth column (``contig``
    is the second), not the third."""
    text = _schema("null_regions")["description"].lower()
    assert "fourth column" in text, text
    assert "third column" not in text, text


def test_null_determinants_cv_r2_says_one_repeat_and_conservative() -> None:
    """Statistician N3: ``cv_r2`` is the pooled out-of-fold R² of one CV repeat; the observed
    value is the mean over repeats (conservative)."""
    desc = _fields("null_determinants")["cv_r2"].get("description", "").lower()
    for word in ("pooled", "out-of-fold", "one cv repeat", "conservative"):
        assert word in desc, (word, desc)
