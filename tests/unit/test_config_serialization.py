"""Config serialisation and ``--set`` typing (GIAP_SPEC.md §5.1, §6.7).

The resolved config is written to JSON (run manifest), so its JSON-mode dump must be strict JSON
(no NaN/Infinity tokens) and must round-trip, including ±inf in ``operons.distance_bins_bp``.
"""

from __future__ import annotations

import json
import math
from importlib import resources
from typing import Any

import jsonschema

from giap.config import GiapConfig, resolve_config


def _strict_json_dump() -> str:
    dump = resolve_config().model_dump(mode="json")
    return json.dumps(dump, allow_nan=False)


def test_default_config_dump_is_strict_json() -> None:
    _strict_json_dump()


def test_default_config_round_trips_through_json() -> None:
    original = resolve_config()
    restored = GiapConfig.model_validate(json.loads(_strict_json_dump()))
    assert restored == original


def test_distance_bins_infinities_survive_round_trip() -> None:
    restored = GiapConfig.model_validate(json.loads(_strict_json_dump()))
    bins = restored.operons.distance_bins_bp
    assert bins[0] == -math.inf
    assert bins[-1] == math.inf


def _config_schema() -> dict[str, Any]:
    text = resources.files("giap.schemas").joinpath("config.schema.json").read_text("utf-8")
    schema: dict[str, Any] = json.loads(text)
    return schema


def test_schema_parses_as_strict_json() -> None:
    text = resources.files("giap.schemas").joinpath("config.schema.json").read_text("utf-8")

    def reject(token: str) -> None:
        raise ValueError(f"non-standard JSON constant {token}")

    json.loads(text, parse_constant=reject)


def test_default_dump_validates_against_config_schema() -> None:
    dump = json.loads(_strict_json_dump())
    jsonschema.validate(dump, _config_schema())


def test_schema_default_for_distance_bins_validates_against_itself() -> None:
    schema = _config_schema()
    validator = jsonschema.Draft202012Validator(schema)
    operons_ref = schema["properties"]["operons"]
    if "$ref" in operons_ref:
        operons = validator.resolver.resolve(operons_ref["$ref"])[1]
    else:
        operons = operons_ref
    prop = operons["properties"]["distance_bins_bp"]
    assert "default" in prop
    sub = {**schema, **prop}
    jsonschema.validate(prop["default"], sub)


def test_set_numeric_looking_sample_name_stays_string() -> None:
    config = resolve_config(overrides=["sample_name=2024"])
    assert config.sample_name == "2024"
    assert isinstance(config.sample_name, str)


def test_set_numeric_looking_schema_version_stays_string() -> None:
    config = resolve_config(overrides=["schema_version=1.0"])
    assert config.schema_version == "1.0"
    assert isinstance(config.schema_version, str)
