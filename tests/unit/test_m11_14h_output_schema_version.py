"""M11-14h: ``output_schema_version`` (tables + manifest), apart from the config's
``schema_version``.

docs/specs/M11-14-container-install.md §4 "Output schema" ("14h adds ``output_schema_version``
(tables + manifest), ``MAJOR.MINOR``, starting at ``1.0`` with the 1.0.0 release ... The config's
``schema_version`` keeps its meaning, the config format (Q-L13)"); ADR-0060 §6.

Interface pinned here (test-author reading, beyond the spec text):

- One source: ``gippa.schemas.OUTPUT_SCHEMA_VERSION``, a ``"MAJOR.MINOR"`` string.
- The manifest: ``build_manifest`` writes the top-level ``output_schema_version`` (that
  constant); ``run_manifest.schema.json`` requires it with the pattern ``^[0-9]+\\.[0-9]+$``.
  ``schema_version`` is still the config's.
- The tables: ``src/gippa/schemas/datapackage.json`` (copied to ``OUT/schemas/``) carries the same
  top-level ``output_schema_version``; the schema diff test reads it from there.
- The value before the first release is not pinned beyond major ``1``: the spec starts the
  counter at ``1.0`` with 1.0.0, and the 14d manifest keys (``operon_source*``) come before that
  release (see the PR body).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import jsonschema
import pytest
from helpers import manifest_data as md

REPO = Path(__file__).resolve().parents[2]
SCHEMAS = REPO / "src" / "gippa" / "schemas"
PATTERN = re.compile(r"^[0-9]+\.[0-9]+$")


def _version() -> str:
    from gippa.schemas import OUTPUT_SCHEMA_VERSION

    assert isinstance(OUTPUT_SCHEMA_VERSION, str)
    return OUTPUT_SCHEMA_VERSION


def test_constant_is_major_minor_with_major_1() -> None:
    value = _version()
    assert PATTERN.match(value), value
    assert value.split(".")[0] == "1"


def test_manifest_carries_output_schema_version() -> None:
    m = md.build()
    assert m["output_schema_version"] == _version()


def test_config_schema_version_keeps_its_meaning() -> None:
    """Q-L13: ``schema_version`` stays the config format, copied from the config."""
    from gippa.config import GippaConfig

    cfg = GippaConfig(seed=7)
    m = md.build(cfg=cfg)
    assert m["schema_version"] == cfg.schema_version


def test_manifest_schema_requires_the_key() -> None:
    schema = md.schema()
    assert "output_schema_version" in schema["required"]
    m = md.build()
    md.validate(m)
    m.pop("output_schema_version")
    with pytest.raises(jsonschema.ValidationError):
        md.validate(m)


@pytest.mark.parametrize("bad", ["1", "v1.0", "1.0.0", "1.x", ""])
def test_manifest_schema_rejects_a_value_that_is_not_major_minor(bad: str) -> None:
    m = md.build()
    m["output_schema_version"] = bad
    with pytest.raises(jsonschema.ValidationError):
        md.validate(m)


def test_datapackage_carries_the_same_version() -> None:
    doc = json.loads((SCHEMAS / "datapackage.json").read_text(encoding="utf-8"))
    assert doc.get("output_schema_version") == _version()
