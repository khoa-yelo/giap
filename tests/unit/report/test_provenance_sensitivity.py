"""M9-16 behaviour 5 (sensitivity section), behaviour 6 (provenance, with and without
``run_manifest.json``; lead decision 1), example 6, the "databases block missing" degenerate row
and determinism (docs/specs/M9-16-report-content.md; GIAP_SPEC.md §6.7, §8.6, §10.3 section 13;
M9-index decision 26).

Interface pinned here (M9-16 test author, beyond the spec text):
- the manifest is ``OUT/run_manifest.json``; the "manifest digest" is the SHA-256 of its bytes;
- the command line is the manifest's ``command_line`` joined by single spaces;
- an empty ``databases`` list (what the MVP-03a writer emits until MVP-03c) reads as "not
  recorded", like a missing block; an unreadable manifest is a note naming ``run_manifest.json``,
  never an exception, and the page falls back to the pass-1 subset;
- each sensitivity family is its own ``<table data-sensitivity-family="<family>">`` in section
  ``sensitivity``, holding every row of that family (no top-20 cut) in file order, columns as
  ``sensitivity.tsv``; ``jaccard`` is ``%.2f``, ``NA``, or "no flagged set in either run" when
  ``n_primary`` and ``n_variant`` are both 0.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from helpers.report_content import (
    COMMIT,
    GO_RELEASE,
    KEGG_DATE,
    MANIFEST_INPUTS,
    SEED,
    TOOL_DIGEST,
    det_rows,
    html,
    manifest,
    page,
    section_text,
    sha256_file,
    write_category_summary,
    write_det_model,
    write_genes,
    write_islands,
    write_manifest,
    write_missed_event,
    write_pairs,
    write_sensitivity,
    write_tier_fdr,
)
from helpers.report_fixtures import INPUTS, Node, run_out, sha256, table_rows
from helpers.viz_out_tree import schema_columns, write_resolved_config

import giap

NOT_WRITTEN = "run_manifest.json not written by this GIAP version"


# ---------------------------------------------------------------- provenance with a manifest


def test_provenance_from_manifest(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    path = write_manifest(out)
    text = section_text(page(out), "provenance")
    doc = manifest()
    assert doc["giap_version"] in text and COMMIT in text, text
    assert " ".join(doc["command_line"]) in text, text
    for row in MANIFEST_INPUTS:
        assert str(row["file"]) in text and str(row["sha256"]) in text, row
    assert "diamond" in text and "2.1.9" in text and TOOL_DIGEST in text, text
    for db in doc["databases"]:
        assert db["name"] in text and db["checksum"] in text, db
    assert KEGG_DATE in text and GO_RELEASE in text, text
    assert str(SEED) in text, text
    assert sha256_file(path) in text, "manifest digest (SHA-256 of run_manifest.json)"
    assert NOT_WRITTEN not in text


def test_manifest_digest_changes_with_the_manifest(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    first = write_manifest(out)
    d1 = sha256_file(first)
    assert d1 in section_text(page(out), "provenance")
    second = write_manifest(out, manifest(command_line=["giap", "run", "--seed", "1"]))
    d2 = sha256_file(second)
    text = section_text(page(out), "provenance")
    assert d1 != d2 and d2 in text and d1 not in text


def test_manifest_missing_databases_block_says_not_recorded(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    doc = manifest()
    del doc["databases"]
    write_manifest(out, doc)
    text = section_text(page(out), "provenance")
    assert "not recorded" in text, text
    assert COMMIT in text, "the rest of the manifest is still shown"


def test_manifest_empty_databases_says_not_recorded(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_manifest(out, manifest(databases=[]))
    assert "not recorded" in section_text(page(out), "provenance")


def test_unreadable_manifest_is_a_note(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    (out / "run_manifest.json").write_text("{not json", encoding="utf-8")
    text = section_text(page(out), "provenance")  # must not raise
    assert "run_manifest.json" in text, text
    assert str(INPUTS[0]["sha256"]) in text, "falls back to the preflight SHA-256s"


# ---------------------------------------------------------------- example 6: no manifest


def test_example_6_no_manifest(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    resolved = write_resolved_config(out / "resolved_config.yaml", ["core"])
    text = section_text(page(out), "provenance")
    assert NOT_WRITTEN in text, text
    for item in INPUTS:
        assert str(item["sha256"]) in text, item
    assert giap.__version__ in text and sha256(resolved) in text


# ---------------------------------------------------------------- sensitivity (behaviour 5)


def _families(p_node: Node) -> dict[str, Node]:
    return {
        t.attrs["data-sensitivity-family"]: t for t in p_node.by_attr("data-sensitivity-family")
    }


def test_sensitivity_grouped_by_family_whole(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    rows = [("islands", f"block_len_{i}", 10, 12, 9, 9 / 13) for i in range(25)]
    rows += [("regions", "alpha_0.01", 4, 3, 3, 0.75), ("regions", "pi0_1", 0, 0, 0, None)]
    write_sensitivity(out, rows)
    p = page(out)
    fams = _families(p.section("sensitivity"))
    assert sorted(fams) == ["islands", "regions"], sorted(fams)
    header, body = table_rows(fams["islands"])
    assert len(body) == 25, "rendered whole, not the top 20"
    assert [r[header.index("variant")] for r in body] == [f"block_len_{i}" for i in range(25)]
    assert set(header) <= set(schema_columns("sensitivity")) and "jaccard" in header
    assert {r[header.index("jaccard")] for r in body} == {"0.69"}
    rh, rbody = table_rows(fams["regions"])
    jac = {r[rh.index("variant")]: r[rh.index("jaccard")] for r in rbody}
    assert jac == {"alpha_0.01": "0.75", "pi0_1": "no flagged set in either run"}, jac


def test_sensitivity_na_jaccard_with_nonzero_n(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_sensitivity(out, [("operon_flags", "v1", 2, 0, 0, None)])
    header, body = table_rows(_families(page(out).section("sensitivity"))["operon_flags"])
    assert body[0][header.index("jaccard")] == "NA", body


def test_sensitivity_absent_says_not_produced(tmp_path: Path) -> None:
    p = page(run_out(tmp_path))
    assert not _families(p.section("sensitivity"))
    assert "not produced in this run" in section_text(p, "sensitivity")


# ---------------------------------------------------------------- determinism


def _rich_out(tmp: Path) -> Path:
    out = run_out(tmp)
    write_genes(out, 12, 3, n_2pep=8)
    write_tier_fdr(out, {"0.01": 0.0098})
    write_islands(out, [0.001, 0.3])
    write_pairs(out)
    write_missed_event(out)
    write_category_summary(out, {"cog": 0.62, "kegg": 0.41})
    write_det_model(out, det_rows("elastic_net", 0.21, 0.03, 0.001, "null_uncalibrated"))
    write_sensitivity(out, [("islands", "v1", 1, 1, 1, 1.0)])
    write_manifest(out)
    return out


def test_report_is_deterministic(tmp_path: Path) -> None:
    out = _rich_out(tmp_path / "a")
    assert html(out) == html(out)


def test_report_does_not_depend_on_the_out_location(tmp_path: Path) -> None:
    a = _rich_out(tmp_path / "a")
    b = tmp_path / "elsewhere" / "OUT"
    shutil.copytree(a, b)
    ha, hb = html(a), html(b)
    assert ha == hb
    assert str(tmp_path) not in ha
