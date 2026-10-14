"""Bundled JSON Schemas (GIAP_SPEC.md §6.7) and the committed config schema."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import jsonschema
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCHEMA_DIR = ROOT / "src" / "giap" / "schemas"
SCHEMAS = ["run_manifest.schema.json", "config.schema.json"]


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
    from giap.config import resolve_config

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
