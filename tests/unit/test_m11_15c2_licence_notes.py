"""M11-15c2: ``licence_notes`` in ``run_manifest.json`` and the report's Methods section.

docs/specs/M11-15-portable-clusters.md item 6 ("Record of use (15c2)"), example 6.6, item 13
(licence notices recorded), Khoa's decision Q-K35 ("document the license stuff").

Interface pinned here (readings chosen where the spec leaves room):

- ``build_manifest(ManifestSources)`` returns a ``licence_notes`` list of strings, always present
  (``[]`` when nothing licence-restricted was used), in the order geNOMAD, KEGG, then any others;
  the two notes are exactly ``"geNOMAD: non-commercial"`` and ``"KEGG: academic use"``.
- "Used" = what the run executed: geNOMAD when module ``mobile`` is on and the GENOMAD process ran
  (trace), KEGG when module ``function`` is on and its step ran with the KEGG resources of the DB
  directory. A run with ``modules: [core]`` records neither. The UniOP "no licence" note (14u2)
  joins the list as an entry starting with ``UniOP`` when the operons came from UniOP.
- ``licence_notes`` is in ``run_manifest.schema.json`` (the manifest still validates).
- The report's Methods section shows every note of the manifest (one line), and none without.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from helpers.manifest_data import build, spec, trace_row, validate

GENOMAD_NOTE = "geNOMAD: non-commercial"
KEGG_NOTE = "KEGG: academic use"
KEGG_DATE = "2026-10-15"
DB_MANIFEST = {
    "schema_version": "1.0",
    "items": [
        {"name": "genomad", "version": "1.7", "url": "https://zenodo.org/x", "checksum":
         "sha256:" + "1" * 64, "size_bytes": 10, "retrieved": KEGG_DATE, "path": "genomad_db",
         "status": "ok", "message": ""},
        {"name": "kegg_ko_pathway", "version": KEGG_DATE, "url": "https://rest.kegg.jp/x",
         "checksum": "sha256:" + "8" * 64, "size_bytes": 2000, "retrieved": KEGG_DATE,
         "path": "function/kegg/ko_pathway.tsv", "status": "ok", "message": ""},
        {"name": "kegg_pathway_names", "version": KEGG_DATE, "url": "https://rest.kegg.jp/y",
         "checksum": "sha256:" + "7" * 64, "size_bytes": 300, "retrieved": KEGG_DATE,
         "path": "function/kegg/pathway_names.tsv", "status": "ok", "message": ""},
    ],
}  # fmt: skip


def _cfg(modules: list[str]) -> Any:
    from gippa.config import GippaConfig

    return GippaConfig(seed=7, modules=modules)  # type: ignore[arg-type]


def _run(modules: list[str], processes: list[str], **over: Any) -> dict[str, Any]:
    """A Nextflow run with ``modules`` on, whose trace completed ``processes``."""
    registered = [spec("detectability"), spec("gene_table"), spec("regions", module="mobile"),
                  spec("function", module="function")]  # fmt: skip
    base: dict[str, Any] = {
        "cfg": _cfg(modules),
        "registered": registered,
        "executor": "nextflow",
        "trace": [trace_row(p) for p in processes],
        "db_manifest": DB_MANIFEST,
        "tool_reports": [{"process": "GENOMAD", "name": "genomad", "version": "1.7.4",
                          "container_digest": None, "conda_lock_hash": None}]
        if any("GENOMAD" in p for p in processes) else [],
    }  # fmt: skip
    base.update(over)
    return build(**base)


EVERYTHING = ["GIPPA:MOBILE:GENOMAD (sample)", "GIPPA:MOBILE:GIPPA_STEP (regions)",
              "GIPPA:FUNCTION:FUNCTION_STEP (function)"]  # fmt: skip


# ---------------------------------------------------------------- manifest (example 6.6)


def test_run_with_genomad_and_kegg_lists_both() -> None:
    m = _run(["core", "mobile", "bgc", "function"], EVERYTHING)
    assert m["licence_notes"] == [GENOMAD_NOTE, KEGG_NOTE]


def test_core_only_run_lists_nothing() -> None:
    m = _run(["core"], [], db_manifest=DB_MANIFEST)
    assert m["licence_notes"] == []


def test_installed_but_unused_items_are_not_listed() -> None:
    """The DB directory has geNOMAD and KEGG; the run had mobile and function off."""
    m = _run(["core", "bgc"], [])
    assert GENOMAD_NOTE not in m["licence_notes"] and KEGG_NOTE not in m["licence_notes"]


def test_mobile_without_function_lists_only_genomad() -> None:
    m = _run(["core", "mobile"], EVERYTHING[:2])
    assert m["licence_notes"] == [GENOMAD_NOTE]


def test_function_without_mobile_lists_only_kegg() -> None:
    m = _run(["core", "function"], EVERYTHING[2:])
    assert m["licence_notes"] == [KEGG_NOTE]


def test_uniop_operons_add_the_uniop_note() -> None:
    m = _run(["core"], [], operon_meta={"operon_source": "uniop", "reason": "UniOP"})
    assert any(n.startswith("UniOP") for n in m["licence_notes"]), m["licence_notes"]


def test_distance_model_operons_add_no_uniop_note() -> None:
    m = _run(["core"], [], operon_meta={"operon_source": "distance_model", "reason": "x"})
    assert not any(n.startswith("UniOP") for n in m["licence_notes"])


def test_manifest_with_licence_notes_validates() -> None:
    m = _run(["core", "mobile", "bgc", "function"], EVERYTHING)
    assert m["licence_notes"]
    validate(m)


# ---------------------------------------------------------------- report Methods section


def _methods(tmp_path: Path, doc: dict[str, Any] | None) -> str:
    from helpers.report_content import page, section_text
    from helpers.report_fixtures import run_out

    out = run_out(tmp_path)
    if doc is not None:
        (out / "run_manifest.json").write_text(json.dumps(doc), encoding="utf-8")
    return section_text(page(out), "methods")


def test_report_methods_shows_the_licence_notes(tmp_path: Path) -> None:
    text = _methods(tmp_path, {"licence_notes": [GENOMAD_NOTE, KEGG_NOTE]})
    assert GENOMAD_NOTE in text and KEGG_NOTE in text
    assert text.index(GENOMAD_NOTE) < text.index(KEGG_NOTE)


@pytest.mark.parametrize("doc", [None, {"licence_notes": []}, {}])
def test_report_methods_without_notes_shows_none(tmp_path: Path, doc: Any) -> None:
    text = _methods(tmp_path, doc)
    assert GENOMAD_NOTE not in text and KEGG_NOTE not in text
