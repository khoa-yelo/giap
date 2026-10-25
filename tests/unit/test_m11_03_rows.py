"""``readiness.py compute``: every row on a fixture repository (issue M11-03, examples 1-10, 13;
edge cases; acceptance criteria 1-2).

Spec: docs/specs/M11-03-release-readiness.md. Fixture: ``helpers.release_readiness`` (a two-commit
git repository in ``tmp_path`` whose evidence passes every row). Every run passes ``--offline``;
the H6 check runs come from ``--check-runs FILE`` (a saved ``gh api .../check-runs`` response), the
offline input that lets example 1 reach ``overall = pass`` without the network.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from helpers.release_readiness import (
    ROW_IDS,
    SOURCES,
    STATUSES,
    Fixture,
    adr_text,
    baseline_doc,
    git,
    metrics_doc,
    nightly_row,
    passing_nightly,
    progress_md,
    report_md,
    row,
    text_of,
)


@pytest.fixture
def fx(tmp_path: Path) -> Fixture:
    return Fixture.build(tmp_path)


def _metrics(fx: Fixture, **kw: object) -> dict:
    doc = metrics_doc(**kw)  # type: ignore[arg-type]
    fx.write_json("validation/pilot/metrics/metrics.json", doc)
    return doc


# ------------------------------------------------------------------ example 1 and the JSON layout


def test_example1_all_pass(fx: Fixture) -> None:
    proc, doc = fx.compute()
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert doc is not None and doc["overall"] == "pass"
    r1, r7 = row(doc, "R1"), row(doc, "R7")
    assert r1["status"] == "pass" and "99.71 %" in r1["value"]
    assert r7["status"] == "pass" and "0:38" in r7["value"]
    assert "wheat" in text_of(r7)  # node class shown
    bad = {r["id"]: r["status"] for r in doc["rows"] if r["status"] not in ("pass", "deferred")}
    assert bad == {}


def test_json_layout(fx: Fixture) -> None:
    _, doc = fx.compute()
    assert doc is not None
    assert doc["version"] == "1.0.0"
    assert doc["commit"] == fx.r and re.fullmatch(r"[0-9a-f]{40}", doc["commit"])
    assert re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(\.\d+)?(Z|\+00:00)", doc["computed"])
    assert [r["id"] for r in doc["rows"]] == ROW_IDS
    for r in doc["rows"]:
        assert set(r) >= {"id", "criterion", "status", "value", "threshold", "evidence",
                          "waiver", "source"}  # fmt: skip
        assert r["status"] in STATUSES, r
        assert r["source"] == SOURCES[r["id"]], r
        assert isinstance(r["evidence"], list)
        assert all(isinstance(e, str) for e in r["evidence"])
        assert r["waiver"] is None or re.fullmatch(r"ADR-\d{4}", r["waiver"])
    assert (fx.out / "READINESS.md").is_file()


def test_r10_always_deferred(fx: Fixture) -> None:
    _, doc = fx.compute()
    assert row(doc, "R10")["status"] == "deferred"
    assert "co-abundance" in text_of(row(doc, "R10"))


# ------------------------------------------------------------------ R1 (examples 2, 3)


def test_example2_r1_below_bound_fails(fx: Fixture) -> None:
    _metrics(fx, share=0.9949)
    proc, doc = fx.compute()
    assert proc.returncode == 1
    r1 = row(doc, "R1")
    assert r1["status"] == "fail" and "99.49 % < 99.5 %" in text_of(r1)
    assert doc is not None and doc["overall"] == "fail"


def test_example3_r1_bound_inclusive(fx: Fixture) -> None:
    _metrics(fx, share=0.995)
    proc, doc = fx.compute()
    assert row(doc, "R1")["status"] == "pass" and proc.returncode == 0


def test_r1_needs_te1_passed(fx: Fixture) -> None:
    _metrics(fx, share=0.9971, passed=False)
    _, doc = fx.compute()
    assert row(doc, "R1")["status"] == "fail"


def test_missing_metrics_is_pending_not_crash(fx: Fixture) -> None:
    (fx.repo / "validation/pilot/metrics/metrics.json").unlink()
    fx.commit("drop metrics")
    proc, doc = fx.compute()
    assert proc.returncode == 3, proc.stdout + proc.stderr
    r1 = row(doc, "R1")
    assert r1["status"] == "pending"
    assert "missing: validation/pilot/metrics/metrics.json" in text_of(r1)
    assert doc is not None and doc["overall"] == "pending"


# ------------------------------------------------------------------ R2


@pytest.mark.parametrize(("fdr", "status"), [(0.0083, "pass"), (0.01, "pass"), (0.0101, "fail")])
def test_r2_mapping_fdr_key(fx: Fixture, fdr: float, status: str) -> None:
    doc = metrics_doc()
    doc["mapping_fdr"]["value"] = fdr
    fx.write_json("validation/pilot/metrics/metrics.json", doc)
    _, out = fx.compute()
    assert row(out, "R2")["status"] == status


def test_r2_no_key_no_freeze_run_pending(fx: Fixture) -> None:
    doc = metrics_doc()
    del doc["mapping_fdr"]
    fx.write_json("validation/pilot/metrics/metrics.json", doc)
    proc, out = fx.compute()
    assert row(out, "R2")["status"] == "pending"
    assert proc.returncode == 3


def test_r2_freeze_run_without_table_pending(fx: Fixture, tmp_path: Path) -> None:
    doc = metrics_doc()
    del doc["mapping_fdr"]
    fx.write_json("validation/pilot/metrics/metrics.json", doc)
    (tmp_path / "freeze").mkdir()
    _, out = fx.compute("--freeze-run", str(tmp_path / "freeze"))
    r2 = row(out, "R2")
    assert r2["status"] == "pending" and "mapping_fdr.tsv" in text_of(r2)


# ------------------------------------------------------------------ R3, R4 (examples 9, 10)


def test_example9_code_change_after_logged_commit_pending(fx: Fixture) -> None:
    fx.set_nightly(passing_nightly(fx.short(fx.r)))
    x = fx.r
    fx.write_text("src/giap/stats/bh.py", '"""BH (fixture), changed."""\n')
    fx.commit("change bh")
    _, doc = fx.compute()
    r3 = row(doc, "R3")
    assert r3["status"] == "pending", (x, r3)
    assert "no calibration run for R" in text_of(r3)


@pytest.mark.parametrize("path", ["tests/calibration/test_x.py", "ci/gates.yaml"])
def test_r3_other_guarded_paths_pending(fx: Fixture, path: str) -> None:
    fx.write_text(path, "changed\n")
    fx.commit("guarded change")
    _, doc = fx.compute()
    assert row(doc, "R3")["status"] == "pending"


def test_example10_docs_only_change_uses_ancestor_row(fx: Fixture) -> None:
    x = fx.short(fx.r)
    fx.set_nightly(passing_nightly(x, date="2026-11-02", job="2000"))
    fx.write_text("docs/user-guide/install.md", "# Install\n\nMore.\n")
    fx.commit("docs only")
    _, doc = fx.compute()
    r3 = row(doc, "R3")
    assert r3["status"] == "pass"
    assert any(
        "nightly.tsv" in e and "2026-11-02" in e and "2000_" in e and x in e for e in r3["evidence"]
    ), r3["evidence"]


def test_r3_failed_count_fails_even_with_exit0(fx: Fixture) -> None:
    c = fx.short(fx.r)
    rows = passing_nightly(c).replace(
        "calibration_full: 0 failed, 3120 passed", "calibration_full: 2 failed, 3118 passed"
    )
    fx.set_nightly(rows)
    _, doc = fx.compute()
    assert row(doc, "R3")["status"] == "fail"


def test_r3_nonzero_exit_fails(fx: Fixture) -> None:
    c = fx.short(fx.r)
    fx.set_nightly(
        passing_nightly(c) + nightly_row("2026-11-03", "2100_2", c, 1,
                                         "calibration_nightly: 17290 passed in 20000.00s")
    )  # fmt: skip
    _, doc = fx.compute()
    assert row(doc, "R3")["status"] == "fail"


def test_r3_later_failure_beats_earlier_pass(fx: Fixture) -> None:
    c = fx.short(fx.r)
    fx.set_nightly(
        passing_nightly(c, date="2026-11-02", job="2000")
        + nightly_row("2026-11-03", "2001_3", c, 1, "calibration_full: 1 failed, 3119 passed")
    )
    _, doc = fx.compute()
    assert row(doc, "R3")["status"] == "fail"


def test_r3_same_date_higher_job_is_newest(fx: Fixture) -> None:
    c = fx.short(fx.r)
    fx.set_nightly(
        nightly_row("2026-11-02", "2005_3", c, 1, "calibration_full: 1 failed, 3119 passed")
        + passing_nightly(c, date="2026-11-02", job="2009")
    )
    _, doc = fx.compute()
    assert row(doc, "R3")["status"] == "pass"


def test_r3_no_full_row_pending(fx: Fixture) -> None:
    c = fx.short(fx.r)
    fx.set_nightly(
        "".join(ln + "\n" for ln in passing_nightly(c).splitlines() if "calibration_full" not in ln)
    )
    _, doc = fx.compute()
    assert row(doc, "R3")["status"] == "pending"


def test_malformed_lines_skipped_never_pass(fx: Fixture) -> None:
    c = fx.short(fx.r)
    good = [ln for ln in passing_nightly(c).splitlines() if "calibration_full" not in ln]
    bad = ["2026-11-02\t2000_3\t" + c + "\t0", "garbage"]  # 4 fields (killed task), 1 field
    fx.set_nightly("\n".join(good + bad) + "\n")
    proc, doc = fx.compute()
    assert proc.returncode != 0
    assert row(doc, "R3")["status"] == "pending"
    md = (fx.out / "READINESS.md").read_text(encoding="utf-8")
    notes = md.split("## Notes", 1)
    assert len(notes) == 2, "READINESS.md has no Notes section"
    assert re.search(r"\b2\b[^\n]*malformed|malformed[^\n]*\b2\b", notes[1]), notes[1]


def test_r4_attack_row(fx: Fixture) -> None:
    c = fx.short(fx.r)
    fx.set_nightly(passing_nightly(c) + nightly_row(
        "2026-11-03", "2100_4", c, 1, "calibration_attacks: 3 failed, 237 passed"))  # fmt: skip
    _, doc = fx.compute()
    assert row(doc, "R4")["status"] == "fail"


def test_r4_no_attack_row_pending(fx: Fixture) -> None:
    c = fx.short(fx.r)
    fx.set_nightly(
        "".join(ln + "\n" for ln in passing_nightly(c).splitlines() if "attacks" not in ln)
    )
    _, doc = fx.compute()
    r4 = row(doc, "R4")
    assert r4["status"] == "pending" and "no attack row" in text_of(r4)
    assert row(doc, "R3")["status"] == "pass"


def test_no_logs_dir_pending(fx: Fixture) -> None:
    proc, doc = fx.compute(logs=False)
    assert row(doc, "R3")["status"] == "pending"
    assert proc.returncode == 3


# ------------------------------------------------------------------ R5, R6 (example 7)


@pytest.mark.parametrize("status", ["fail", "reported (not passed, explained)"])
def test_example7_reported_te_fails_r5(fx: Fixture, status: str) -> None:
    doc = metrics_doc()
    doc["te_real"]["te3"]["status"] = status
    fx.write_json("validation/pilot/metrics/metrics.json", doc)
    proc, out = fx.compute()
    r5 = row(out, "R5")
    assert r5["status"] == "fail" and "TE3" in text_of(r5)
    assert proc.returncode == 1


def test_r5_not_evaluable_te_deferred_with_reason(fx: Fixture) -> None:
    doc = metrics_doc()
    doc["te_real"]["te4"] = {"status": "not_evaluable", "reason": ["no_ribosomal_genes"]}
    fx.write_json("validation/pilot/metrics/metrics.json", doc)
    proc, out = fx.compute()
    r5 = row(out, "R5")
    assert r5["status"] in ("pass", "deferred")
    assert "TE4" in text_of(r5) and "no_ribosomal_genes" in text_of(r5)
    assert proc.returncode == 0


def test_r5_no_weekly_te8_row_pending(fx: Fixture) -> None:
    fx.set_weekly("")
    proc, out = fx.compute()
    r5 = row(out, "R5")
    assert r5["status"] == "pending" and "TE8" in text_of(r5)
    assert proc.returncode == 3


def test_r6_te7_giap_fail_fails(fx: Fixture) -> None:
    doc = metrics_doc()
    doc["te7"]["giap"]["status"] = "fail"
    fx.write_json("validation/pilot/metrics/metrics.json", doc)
    _, out = fx.compute()
    assert row(out, "R6")["status"] == "fail"


# ------------------------------------------------------------------ R7, R8 (examples 4-6)


def test_example4_runtime_above_range_fails(fx: Fixture) -> None:
    fx.write_json("benchmarks/baseline.json", baseline_doc(standard_s=3400.0))
    proc, doc = fx.compute()
    r7 = row(doc, "R7")
    assert r7["status"] == "fail" and "0:57 > 0:55" in text_of(r7)
    assert proc.returncode == 1


@pytest.mark.parametrize(
    ("secs", "status"), [(960.0, "pass"), (3300.0, "pass"), (959.0, "fail"), (3301.0, "fail")]
)
def test_r7_bounds_inclusive(fx: Fixture, secs: float, status: str) -> None:
    fx.write_json("benchmarks/baseline.json", baseline_doc(standard_s=secs))
    _, doc = fx.compute()
    assert row(doc, "R7")["status"] == status


def test_example5_accepted_adr_waives(fx: Fixture) -> None:
    fx.write_json("benchmarks/baseline.json", baseline_doc(standard_s=3400.0))
    fx.write_text("docs/adr/ADR-0060-runtime-waiver.md", adr_text(60, "accepted-by-khoa"))
    proc, doc = fx.compute()
    r7 = row(doc, "R7")
    assert r7["status"] == "waived" and r7["waiver"] == "ADR-0060"
    assert any("ADR-0060" in e for e in r7["evidence"])
    assert doc is not None and doc["overall"] == "pass" and proc.returncode == 0


def test_example6_proposed_adr_does_not_waive(fx: Fixture) -> None:
    fx.write_json("benchmarks/baseline.json", baseline_doc(standard_s=3400.0))
    fx.write_text("docs/adr/ADR-0060-runtime-waiver.md", adr_text(60, "proposed"))
    proc, doc = fx.compute()
    r7 = row(doc, "R7")
    assert r7["status"] == "fail" and r7["waiver"] is None
    assert any("ADR-0060" in e and "not accepted" in e for e in r7["evidence"])
    assert proc.returncode == 1


def test_waiver_only_for_named_row(fx: Fixture) -> None:
    _metrics(fx, share=0.9949)
    fx.write_text("docs/adr/ADR-0060-runtime-waiver.md", adr_text(60, "accepted-by-khoa", "R7"))
    _, doc = fx.compute()
    assert row(doc, "R1")["status"] == "fail"


def test_r7_missing_baseline_entry_pending(fx: Fixture) -> None:
    fx.write_json("benchmarks/baseline.json", baseline_doc(standard_s=None, core_s=None))
    _, doc = fx.compute()
    assert row(doc, "R7")["status"] == "pending"
    assert row(doc, "R8")["status"] == "pending"


def test_r8_reported_never_fails(fx: Fixture) -> None:
    fx.write_json("benchmarks/baseline.json", baseline_doc(core_s=99999.0))
    _, doc = fx.compute()
    assert row(doc, "R8")["status"] == "pass"


# ------------------------------------------------------------------ R9, H1 (example 8)


def test_r9_needs_khoa_review_line(fx: Fixture) -> None:
    fx.write_text("PROGRESS.md", progress_md(khoa_m10=False))
    _, doc = fx.compute()
    assert row(doc, "R9")["status"] != "pass"


@pytest.mark.parametrize("bad", ["fail", "pending (M10-06)"])
def test_r9_report_table_not_clean(fx: Fixture, bad: str) -> None:
    fx.write_text("validation/pilot/REPORT.md", report_md(("pass", "pass", bad, "pass", "pass")))
    _, doc = fx.compute()
    assert row(doc, "R9")["status"] != "pass"


def test_r9_missing_report_pending(fx: Fixture) -> None:
    (fx.repo / "validation/pilot/REPORT.md").unlink()
    _, doc = fx.compute()
    r9 = row(doc, "R9")
    assert r9["status"] == "pending" and "missing: validation/pilot/REPORT.md" in text_of(r9)


def test_example8_missing_gate_block(fx: Fixture) -> None:
    fx.write_text("PROGRESS.md", progress_md([1, 2, 3, 4, 5, 6, 7, 8, 10]))
    proc, doc = fx.compute()
    h1 = row(doc, "H1")
    assert h1["status"] == "fail" and "missing: M9" in text_of(h1)
    assert proc.returncode == 1


def test_h1_gate_m1_not_matched_by_m10(fx: Fixture) -> None:
    fx.write_text("PROGRESS.md", progress_md(range(2, 11)))
    _, doc = fx.compute()
    assert "missing: M1" in text_of(row(doc, "H1"))


# ------------------------------------------------------------------ H2, H3


@pytest.mark.parametrize("where", ["init", "citation"])
def test_h2_version_disagreement(fx: Fixture, where: str) -> None:
    if where == "init":
        fx.write_text("src/giap/__init__.py", '__version__ = "0.9.0"\n')
    else:
        fx.write_text("CITATION.cff", 'cff-version: 1.2.0\ntitle: "GIAP"\nversion: 0.9.0\n')
    _, doc = fx.compute()
    assert row(doc, "H2")["status"] == "fail"


def test_h2_tag_exists_fails_and_is_kept(fx: Fixture) -> None:
    git(fx.repo, "tag", "v1.0.0")
    _, doc = fx.compute()
    h2 = row(doc, "H2")
    assert h2["status"] == "fail" and "tag v1.0.0 exists" in text_of(h2)
    assert git(fx.repo, "tag", "-l") == "v1.0.0"  # never deletes tags
    assert row(doc, "R1")["status"] == "pass"


def test_compute_never_creates_tags(fx: Fixture) -> None:
    _, doc = fx.compute()
    assert row(doc, "H2")["status"] == "pass"
    assert git(fx.repo, "tag", "-l") == ""


def test_h3_unbuilt_fragment_fails(fx: Fixture) -> None:
    fx.write_text("changelog.d/+m11-release.feature.md", "Release tooling.\n")
    _, doc = fx.compute()
    assert row(doc, "H3")["status"] == "fail"
    assert "+m11-release.feature.md" in text_of(row(doc, "H3"))


def test_h3_no_release_section_fails(fx: Fixture) -> None:
    fx.write_text("CHANGELOG.md", "# Changelog\n\n<!-- towncrier release notes start -->\n")
    _, doc = fx.compute()
    assert row(doc, "H3")["status"] == "fail"


# ------------------------------------------------------------------ H4, H5 (example 13)


def _containers(fx: Fixture) -> Path:
    return fx.repo / "workflow/conf/containers.config"


def test_example13_tag_not_digest(fx: Fixture) -> None:
    p = _containers(fx)
    text = re.sub(
        r"container = 'docker://oschwengers/bakta@sha256:[0-9a-f]{64}'",
        "container = 'docker://oschwengers/bakta:v1.12.2'",
        p.read_text(encoding="utf-8"),
    )
    assert "bakta:v1.12.2'" in text
    p.write_text(text, encoding="utf-8")
    _, doc = fx.compute()
    h4 = row(doc, "H4")
    assert h4["status"] == "fail" and "tag, not digest: BAKTA" in text_of(h4)


def test_h4_registry_tool_neither_pinned_nor_conda_only(fx: Fixture) -> None:
    p = _containers(fx)
    lines = [ln for ln in p.read_text(encoding="utf-8").splitlines() if "'GENOMAD'" not in ln]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    _, doc = fx.compute()
    h4 = row(doc, "H4")
    assert h4["status"] == "fail" and "genomad" in text_of(h4).lower()


def test_h5_tool_missing_from_notice(fx: Fixture) -> None:
    p = fx.repo / "NOTICE.md"
    lines = [
        ln for ln in p.read_text(encoding="utf-8").splitlines() if not ln.startswith("| geNOMAD")
    ]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    _, doc = fx.compute()
    h5 = row(doc, "H5")
    assert h5["status"] == "fail" and "genomad" in text_of(h5).lower()


def test_h5_licence_not_mit(fx: Fixture) -> None:
    fx.write_text("LICENSE", "Apache License\nVersion 2.0, January 2004\n")
    _, doc = fx.compute()
    assert row(doc, "H5")["status"] == "fail"


# ------------------------------------------------------------------ H6 (acceptance criterion 2)


def test_h6_offline_without_check_runs_pending(fx: Fixture) -> None:
    proc, doc = fx.compute(check_runs=False)
    h6 = row(doc, "H6")
    assert h6["status"] == "pending" and "offline" in text_of(h6)
    assert doc is not None and doc["overall"] == "pending" and proc.returncode == 3


def test_h6_failed_check_run(fx: Fixture) -> None:
    fx.set_check_runs(["lint", "tests (3.11)"], {"tests (3.11)": "failure"})
    proc, doc = fx.compute()
    h6 = row(doc, "H6")
    assert h6["status"] == "fail" and "tests (3.11)" in text_of(h6)
    assert proc.returncode == 1


def test_offline_touches_no_network(fx: Fixture, monkeypatch: pytest.MonkeyPatch) -> None:
    # a proxy that cannot exist: any HTTP(S) attempt through urllib/requests/gh would fail loudly
    for var in ("HTTPS_PROXY", "HTTP_PROXY", "https_proxy", "http_proxy"):
        monkeypatch.setenv(var, "http://127.0.0.1:9")
    monkeypatch.setenv("GH_HOST", "invalid.invalid")
    proc, doc = fx.compute()
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert row(doc, "H6")["status"] == "pass"  # from --check-runs, not from GitHub
