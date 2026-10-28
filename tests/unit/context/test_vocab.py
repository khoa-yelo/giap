"""M5-09 ``gippa.context.vocab``: the controlled flag vocabulary of the region-and-context tables
(docs/specs/M5-09-context-tables-step.md "Public interface" and the acceptance item
"``CONTEXT_FLAGS`` ⊇ §9.4 Region-and-context flags; schemas list the same set").

GIPPA_SPEC.md §9.4 "Flag vocabulary (controlled; extend only through the schema)", lead decision 36
(M5 adds ``not_estimable``, ``few_values`` and ``ci_unestimated``) and ADR-0021 Amendment 1 /
M5-06 decision 48 (``uncalibrated`` is **dropped**: no per-replicon statistic certifies the
validity domain, so the flag would be false reassurance — "M5-09's schema must not add
``uncalibrated``").

Interface assumed: ``gippa.context.vocab.CONTEXT_FLAGS``, an iterable of ``str``.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[3]
SPEC = REPO / "GIPPA_SPEC.md"
TABLE_DIR = REPO / "src" / "gippa" / "schemas" / "tables"

#: M4 additions the spec lists (the merged ``regions_stats`` / ``detection_summary`` producers)
M4_FLAGS = ("few_shifts", "split_contigs", "hl_variance_nonpositive", "no_reference_replicon",
            "ci_few_units", "depth_adjusted")  # fmt: skip
#: M5 additions (lead decision 36)
M5_FLAGS = ("not_estimable", "few_values", "ci_unestimated")
#: ADR-0021 Amendment 1: never emitted, never in the M5 schemas
UNCALIBRATED = "uncalibrated"


def _flags() -> set[str]:
    from gippa.context.vocab import CONTEXT_FLAGS

    return {str(f) for f in CONTEXT_FLAGS}


def _schema(name: str) -> dict[str, Any]:
    found = [p for p in (TABLE_DIR / f"{name}.json", TABLE_DIR / f"{name}.schema.json")
             if p.exists()]  # fmt: skip
    assert len(found) == 1, found
    data = json.loads(found[0].read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def _schema_flags(name: str) -> set[str]:
    """The flag vocabulary a table schema constrains ``flags`` to (an ``enum`` of the items, or
    the alternatives of a ``;``-separated pattern as ``regions_stats`` writes it)."""
    field = next(f for f in _schema(name)["fields"] if f["name"] == "flags")
    constraints = field.get("constraints", {})
    if "enum" in constraints:
        return set(constraints["enum"])
    items = field.get("items") or {}
    if "enum" in items.get("constraints", {}):
        return set(items["constraints"]["enum"])
    pattern = constraints["pattern"]
    return {m for m in re.findall(r"[a-z][a-z0-9_]*", pattern)}


def _spec_region_and_context_flags() -> set[str]:
    """The "Region and context" row of the §9.4 table."""
    text = SPEC.read_text(encoding="utf-8")
    row = re.search(r"^\| Region and context \| (.+?) \|$", text, re.MULTILINE)
    assert row, "§9.4 'Region and context' row not found"
    return set(re.findall(r"`([a-z][a-z0-9_]*)`", row.group(1)))


# --------------------------------------------------------------------------- the vocabulary


def test_vocabulary_covers_the_spec_row() -> None:
    """⊇ §9.4 "Region and context", minus the dropped ``uncalibrated``."""
    assert _spec_region_and_context_flags() - {UNCALIBRATED} <= _flags()


def test_vocabulary_covers_the_m4_and_m5_additions() -> None:
    assert set(M4_FLAGS) | set(M5_FLAGS) <= _flags()


def test_uncalibrated_is_not_in_the_vocabulary() -> None:
    """ADR-0021 Amendment 1 / decision 48: the flag is dropped, not merely unused."""
    assert UNCALIBRATED not in _flags()


def test_flag_names_are_snake_case_and_unique() -> None:
    from gippa.context.vocab import CONTEXT_FLAGS

    names = list(CONTEXT_FLAGS)
    assert len(names) == len(set(names))
    assert all(re.fullmatch(r"[a-z][a-z0-9_]*", str(n)) for n in names), names


@pytest.mark.parametrize("name", ["autocorrelation", "islands"])
def test_schemas_list_the_same_set(name: str) -> None:
    """Spec §9.4 "extend only through the schema": the two new schemas constrain ``flags`` to
    exactly ``CONTEXT_FLAGS``."""
    assert _schema_flags(name) == _flags()


def test_regions_stats_schema_is_inside_the_vocabulary() -> None:
    """The merged M4 table draws from the same vocabulary (plus the legacy ``uncalibrated``,
    which ADR-0021 Amendment 1 removed from M5 but which M4's schema still lists)."""
    assert _schema_flags("regions_stats") <= _flags() | {UNCALIBRATED}


# --------------------------------------------------------------------------- the producers


def test_every_flag_the_m5_producers_can_emit_is_in_the_vocabulary() -> None:
    """The flag constants of M5-05, M5-06 and M5-08 (the only sources of ``flags`` in the two
    new tables) are all vocabulary members."""
    from gippa.context import autocorrelation as ac
    from gippa.context import islands as isl
    from gippa.context import oriter as ot

    emitted = {
        ac.NOT_ESTIMABLE, ac.FEW_VALUES, ac.CI_UNESTIMATED,
        isl.FEW_VALUES, isl.CI_FEW_UNITS,
        ot.DESCRIPTIVE, ot.NO_ORIC, ot.NOT_ESTIMABLE,
    }  # fmt: skip
    assert emitted <= _flags()
    assert UNCALIBRATED not in emitted
