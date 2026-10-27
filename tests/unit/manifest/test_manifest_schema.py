"""MVP-03a ``run_manifest.schema.json`` extension (docs/specs/MVP-03-manifest-finalize.md
acceptance item 1, lead decisions 13 and 14; GIPPA_SPEC.md §6.7, §9.5).

Pinned: every D-C field of 03a is a declared property; ``additionalProperties: false`` at the top
level; ``flags`` items match ``^[a-z0-9_]+(:[a-z0-9_.-]+)?$``; the existing ``sensitivity``
property stays (03c fills it). The trace column ``container`` is in ``workflow/nextflow.config``
(D-C ``tools`` row: the digest comes from the trace).
"""

from __future__ import annotations

import copy
import re
from pathlib import Path

import pytest
from helpers.manifest_data import build, schema, schema_errors

REPO = Path(__file__).resolve().parents[3]
FIELDS_03A = [
    "schema_version", "gippa_version", "git_commit", "git_dirty", "command_line", "executor",
    "profile", "resume", "outdir", "started_at", "finished_at", "status", "exit_code",
    "failed_step", "failed_log", "publication", "config", "inputs", "inputs_absent", "tools",
    "databases", "annotation", "detection_tier", "seeds", "steps", "skipped", "modules_off",
    "flags", "layers_missing", "warnings", "volatile", "provenance_digest",
]  # fmt: skip


@pytest.mark.parametrize("field", FIELDS_03A)
def test_every_03a_field_is_declared(field: str) -> None:
    assert field in schema()["properties"], field


def test_sensitivity_property_kept_for_03c() -> None:
    assert "sensitivity" in schema()["properties"]


def test_top_level_additional_properties_false() -> None:
    assert schema().get("additionalProperties") is False
    m = build()
    m["not_a_field"] = 1
    assert schema_errors(m), "an unknown top-level key must fail"


@pytest.mark.parametrize("bad", ["Layer_Missing:cog", "layer missing", "a:b:c", "", "x:"])
def test_bad_flag_rejected(bad: str) -> None:
    m = build()
    m["flags"] = [bad]
    assert schema_errors(m), bad


@pytest.mark.parametrize("good", ["layer_missing:cog", "c_assumed", "uncalibrated:operon_distance",
                                  "layer_missing:rbs_novel"])  # fmt: skip
def test_good_flag_accepted(good: str) -> None:
    m = build()
    m["flags"] = [good]
    assert not schema_errors(m), good


def test_status_enum() -> None:
    m = build()
    m["status"] = "running"
    assert schema_errors(m)


@pytest.mark.parametrize("field", ["provenance_digest"])
def test_digest_is_hex_sha256(field: str) -> None:
    m = build()
    assert re.fullmatch(r"[0-9a-f]{64}", m[field])
    bad = copy.deepcopy(m)
    bad[field] = "xyz"
    assert schema_errors(bad)


def test_seeds_values_are_integers() -> None:
    m = build()
    m["seeds"]["gene_table"] = 1.5
    assert schema_errors(m)


def test_skipped_reason_pattern() -> None:
    m = build()
    m["skipped"] = [{"name": "x", "reason": "because"}]
    assert schema_errors(m)
    for reason in ("module_off:mobile", "opt_in_off", "precomputed:terminators", "not_run"):
        m["skipped"] = [{"name": "x", "reason": reason}]
        assert not schema_errors(m), reason


def test_trace_fields_include_container() -> None:
    text = (REPO / "workflow" / "nextflow.config").read_text(encoding="utf-8")
    m = re.search(r"trace\s*\{[^}]*fields\s*=\s*'([^']*)'", text)
    assert m, "no trace.fields in workflow/nextflow.config"
    fields = m.group(1).split(",")
    assert "container" in fields
    for kept in ("name", "status", "exit", "attempt", "realtime", "%cpu", "peak_rss", "workdir"):
        assert kept in fields, kept
