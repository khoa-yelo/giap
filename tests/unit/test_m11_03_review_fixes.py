"""Regression tests for the PR #337 review (review 5644150619, items 1-10, 12, 13) of
``scripts/release/readiness.py`` (issue M11-03; spec docs/specs/M11-03-release-readiness.md).

The lead's decisions pinned here:

(a) a waiver applies to ``fail`` and ``pending`` rows, and the row keeps its ``missing:`` value;
(b) H6: ``success``/``skipped`` pass; ``neutral``, ``cancelled``, ``timed_out``, ``failure`` and
    ``action_required`` fail; a null conclusion (still running) is ``pending``; every run of a
    ``--check-runs`` file must carry ``head_sha`` = the commit under check, else ``fail``;
(c) TE7 reads the real ``metrics.json`` shape, ``te7.gippa.status`` (``validation/pilot/te7.py``
    writes no ``te7.passed``); a ``te7`` or ``te_real`` without ``schema`` is a placeholder
    (``metrics.te_placeholders()``, status ``pending (M10-06)``) and gives ``pending``;
(d) the peptide guard: no lowercasing; a peptide-like run not on a fixed allow-list is redacted,
    so a peptide never reaches ``READINESS.md`` in any case and ``CHANGELOG.md`` keeps its case;
(e) renames count at both paths (``--no-renames``); an ancestor is code-equal when
    ``git diff --quiet <c> HEAD -- <paths>`` holds (merge-history fixture);
(f) R3 stays ``pending`` (never ``fail``) until a ``calibration_full`` nightly row exists.

Every ``compute`` here writes a ``--check-runs`` file whose runs carry ``head_sha`` = HEAD.
"""

from __future__ import annotations

import functools
import importlib.util
import json
import math
import re
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from helpers.release_readiness import (
    REPO,
    Fixture,
    adr_text,
    git,
    load_module,
    metrics_doc,
    nightly_row,
    passing_nightly,
    passing_weekly,
    progress_md,
    report_md,
    row,
    text_of,
)
from hypothesis import assume, given, settings
from hypothesis import strategies as st

METRICS = "validation/pilot/metrics/metrics.json"
AA = "ACDEFGHIKLMNPQRSTVWY"
PEPTIDE_LIKE = re.compile(rf"[{AA}]{{7,}}")  # the M10-05 guard rule
TE_KEYS = ("te2", "te3", "te4", "te5")
ALL_SUCCESS = {"lint": "success", "tests (3.11)": "success", "package": "success"}


# ------------------------------------------------------------------ helpers


@functools.cache
def pilot_metrics() -> ModuleType:
    """``validation/pilot/metrics.py`` (``te_placeholders``, ``share``: the real shapes)."""
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))  # it imports validation.proxies.te_proxies
    spec = importlib.util.spec_from_file_location(
        "gippa_pilot_metrics_m11_03", REPO / "validation/pilot/metrics.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def te7_doc(status: str = "pass") -> dict[str, Any]:
    """``te7`` as ``validation/pilot/te7.py::te7`` writes it (schema gippa-pilot-te7/1): the verdict
    is ``gippa.status`` from ``classify_tree``; there is no ``passed`` key."""
    return {
        "schema": "gippa-pilot-te7/1",
        "spec": "0" * 64,
        "pair": {},
        "gippa": {"status": status, "differing": [], "external_differing": [],
                 "figures_differing": [], "n_gate": 41, "n_external": 2, "n_figures": 9,
                 "n_excluded": 4},
        "diann": {"status": "identical_report", "l0": True, "l3": None},
        "end_to_end": status == "pass",
    }  # fmt: skip


def te_real_doc(statuses: dict[str, str] | None = None) -> dict[str, Any]:
    """``te_real`` as ``validation/pilot/te_real.py`` writes it (schema gippa-pilot-te-real/1)."""
    st_ = {k: "pass" for k in TE_KEYS} | (statuses or {})
    return {"schema": "gippa-pilot-te-real/1", "spec": "0" * 64,
            **{k: {"status": s} for k, s in st_.items()}}  # fmt: skip


def real_metrics(te7: str = "pass", te_real: dict[str, str] | None = None) -> dict[str, Any]:
    doc = metrics_doc()
    doc["te7"] = te7_doc(te7)
    doc["te_real"] = te_real_doc(te_real)
    return doc


def write_runs(
    fx: Fixture, conclusions: dict[str, str | None] | None = None, sha: str | None = None
) -> None:
    """A saved ``gh api .../commits/<sha>/check-runs`` response; ``head_sha`` defaults to HEAD."""
    head = sha or git(fx.repo, "rev-parse", "HEAD")
    runs = [
        {"name": n, "status": "in_progress" if c is None else "completed", "conclusion": c,
         "head_sha": head}
        for n, c in (ALL_SUCCESS if conclusions is None else conclusions).items()
    ]  # fmt: skip
    fx.check_runs.write_text(
        json.dumps({"total_count": len(runs), "check_runs": runs}), encoding="utf-8"
    )


def compute(fx: Fixture, *extra: str, runs: bool = True) -> tuple[Any, dict[str, Any] | None]:
    if runs:
        write_runs(fx)
    return fx.compute(*extra)


def out_of(p: Any) -> str:
    return str(p.stdout) + str(p.stderr)


def md_of(fx: Fixture) -> str:
    return (fx.out / "READINESS.md").read_text(encoding="utf-8")


@pytest.fixture
def fx(tmp_path: Path) -> Fixture:
    return Fixture.build(tmp_path)


def test_real_shape_helpers_match_pilot_code() -> None:
    """The placeholder shape comes from the pilot code itself, not from a hand-made fixture."""
    ph = pilot_metrics().te_placeholders()
    assert "schema" not in ph["te7"] and "passed" not in ph["te7"]
    assert ph["te7"]["gippa"]["status"].startswith("pending")
    assert "schema" not in ph["te_real"]
    assert all(ph["te_real"][k]["status"].startswith("pending") for k in TE_KEYS)
    assert "passed" not in te7_doc()


# ------------------------------------------------------------------ item 1: TE7 = te7.gippa.status


def test_item1_te7_real_shape_passes(fx: Fixture) -> None:
    fx.write_json(METRICS, real_metrics())
    _, doc = compute(fx)
    r5, r6 = row(doc, "R5"), row(doc, "R6")
    assert r6["status"] == "pass", r6
    assert r5["status"] == "pass", r5
    assert re.search(r"TE7 pass", r5["value"]), r5["value"]


def test_item1_te7_gippa_fail_fails(fx: Fixture) -> None:
    fx.write_json(METRICS, real_metrics(te7="fail"))
    _, doc = compute(fx)
    assert row(doc, "R6")["status"] == "fail"
    r5 = row(doc, "R5")
    assert r5["status"] == "fail" and "TE7" in text_of(r5)


def test_item1_invented_passed_key_is_not_read(fx: Fixture) -> None:
    doc = real_metrics(te7="fail")
    doc["te7"]["passed"] = True  # a key the pilot never writes; the verdict is gippa.status
    fx.write_json(METRICS, doc)
    _, out = compute(fx)
    assert row(out, "R6")["status"] == "fail"
    assert row(out, "R5")["status"] == "fail"


def test_item1_te7_placeholder_is_pending(fx: Fixture) -> None:
    doc = real_metrics()
    doc["te7"] = pilot_metrics().te_placeholders()["te7"]
    fx.write_json(METRICS, doc)
    proc, out = compute(fx)
    assert row(out, "R6")["status"] == "pending", row(out, "R6")
    r5 = row(out, "R5")
    assert r5["status"] == "pending", r5
    assert re.search(r"TE7 pending", r5["value"]), r5["value"]
    assert proc.returncode == 3, out_of(proc)


# ------------------------------------------------------------------ item 2: te_real placeholders


def test_item2_te_real_placeholders_pending_not_fail(fx: Fixture) -> None:
    doc = real_metrics()
    doc["te_real"] = pilot_metrics().te_placeholders()["te_real"]
    fx.write_json(METRICS, doc)
    proc, out = compute(fx)
    r5 = row(out, "R5")
    assert r5["status"] == "pending", r5
    assert "fail" not in r5["value"], r5["value"]
    assert proc.returncode == 3, out_of(proc)


def test_item2_pending_status_with_schema_is_pending(fx: Fixture) -> None:
    fx.write_json(METRICS, real_metrics(te_real={"te3": "pending (M10-06)"}))
    _, out = compute(fx)
    r5 = row(out, "R5")
    assert r5["status"] == "pending", r5
    assert re.search(r"TE3 pending", r5["value"]), r5["value"]


# ------------------------------------------------------------------ item 3: R1 never crashes


def test_item3_empty_denominator_fails_with_reason(fx: Fixture) -> None:
    doc = real_metrics()
    s = pilot_metrics().share(0, 0)  # value None, reason empty_denominator
    doc["te1"] = {"n_peptides": 0, "exact_inframe": 0, "share_exact_inframe": s["value"],
                  "ci_low": None, "ci_high": None, "threshold": 0.995, "passed": False,
                  "reason": s["reason"]}  # fmt: skip
    fx.write_json(METRICS, doc)
    proc, out = compute(fx)
    assert "Traceback" not in proc.stderr, proc.stderr
    assert out is not None, out_of(proc)
    r1 = row(out, "R1")
    assert r1["status"] == "fail" and "empty_denominator" in text_of(r1), r1
    assert proc.returncode == 1


@pytest.mark.parametrize("share", [math.nan, math.inf])
def test_item3_non_finite_share_fails(fx: Fixture, share: float) -> None:
    doc = real_metrics()
    doc["te1"]["share_exact_inframe"] = share
    doc["te1"]["passed"] = True
    fx.write_json(METRICS, doc)  # json.dumps writes NaN / Infinity, json.loads reads them back
    proc, out = compute(fx)
    assert "Traceback" not in proc.stderr, proc.stderr
    assert row(out, "R1")["status"] == "fail"


# ------------------------------------------------------------------ item 4: R9 table parsing

TE_TABLE = (
    "\n## 5. True examples\n\n| TE | Status | Headline |\n|---|---|---|\n"
    "| TE3: ribosomal proteins | fail | detected 0.71 |\n"
    "| TE7: byte-identical tables | pending (M10-06) | gate files compared 0 |\n"
)


def test_item4_later_tables_not_read(fx: Fixture) -> None:
    rep = report_md(("pass", "pass", "pass", "pass", "reported (1 not passed, explained)"))
    fx.write_text("validation/pilot/REPORT.md", rep + TE_TABLE)
    _, doc = compute(fx)
    r9 = row(doc, "R9")
    assert r9["status"] == "pass", r9


@pytest.mark.parametrize(("cell", "status"), [("pass", "pass"), ("fail", "fail")])
def test_item4_status_header_with_suffix(fx: Fixture, cell: str, status: str) -> None:
    fx.write_text(
        "validation/pilot/REPORT.md",
        "# r\n\n| §15.6 bullet | Status (§15.6) | Evidence |\n|---|---|---|\n"
        f"| a | pass | x |\n| b | {cell} | y |\n",
    )
    proc, doc = compute(fx)
    assert "Traceback" not in proc.stderr, proc.stderr
    assert row(doc, "R9")["status"] == status


def test_item4_empty_status_table_pending(fx: Fixture) -> None:
    fx.write_text(
        "validation/pilot/REPORT.md", "# r\n\n| §15.6 bullet | Status | Evidence |\n|---|---|---|\n"
    )
    _, doc = compute(fx)
    assert row(doc, "R9")["status"] == "pending"


# ------------------------------------------------------------------ item 5: renames (verify)


def test_item5_rename_out_of_src_is_stale(fx: Fixture) -> None:
    compute(fx)
    fx.commit("readiness paperwork")
    assert "commit rule: pass" in out_of(fx.verify())  # control: paperwork only
    git(fx.repo, "mv", "src/gippa/cli.py", "docs/cli.py")
    git(fx.repo, "commit", "-q", "-m", "move cli into docs")
    p = fx.verify()
    assert p.returncode == 1
    assert "commit rule: fail (stale: src/gippa/cli.py)" in out_of(p), out_of(p)


# ------------------------------------------------------------------ item 6: ancestry with merges


def _merge_history(fx: Fixture) -> tuple[str, str]:
    """main: R - Z (docs) - M (merge); side: R - Y (src/gippa/stats/bh.py). Returns (Y, Z)."""
    git(fx.repo, "checkout", "-q", "-b", "side")
    fx.write_text("src/gippa/stats/bh.py", '"""BH (fixture), changed on side."""\n')
    y = fx.commit("Y: src change on side")
    git(fx.repo, "checkout", "-q", "main")
    fx.write_text("docs/a.md", "# a\n")
    z = fx.commit("Z: docs on main")
    git(fx.repo, "merge", "-q", "--no-ff", "side", "-m", "M: merge side")
    assert git(fx.repo, "diff", "--name-only", z, "HEAD", "--", "src") == "src/gippa/stats/bh.py"
    assert git(fx.repo, "diff", "--name-only", y, "HEAD", "--", "src") == ""
    return y, z


def test_item6_pre_merge_main_commit_not_code_equal(fx: Fixture) -> None:
    _, z = _merge_history(fx)
    fx.set_nightly(passing_nightly(fx.short(z)))
    _, doc = compute(fx)
    r3 = row(doc, "R3")
    assert r3["status"] == "pending", r3


def test_item6_code_equal_side_commit_counts(fx: Fixture) -> None:
    y, _ = _merge_history(fx)
    fx.set_nightly(passing_nightly(fx.short(y)))
    _, doc = compute(fx)
    r3 = row(doc, "R3")
    assert r3["status"] == "pass", r3
    assert any(fx.short(y) in e for e in r3["evidence"]), r3["evidence"]


def test_item6_rename_out_of_src_invalidates_row(fx: Fixture) -> None:
    git(fx.repo, "mv", "src/gippa/stats/bh.py", "docs/bh.py")
    git(fx.repo, "commit", "-q", "-m", "move bh into docs")
    _, doc = compute(fx)
    assert row(doc, "R3")["status"] == "pending"


# ------------------------------------------------------------------ item 7: peptide guard

PEP_STATUS, PEP_REASON, PEP_LOG = "PEPTIDEKLMNQR", "MKVLAAGIIR", "WQERTYHGFD"


def test_item7_peptides_never_reach_readiness_md(fx: Fixture) -> None:
    doc = real_metrics(te_real={"te3": PEP_STATUS})
    doc["te_real"]["te4"] = {
        "status": "not_evaluable",
        "reason": ["no_ribosomal_genes", PEP_REASON],
    }
    fx.write_json(METRICS, doc)
    c = fx.short(fx.r)
    fx.set_nightly(
        passing_nightly(c).replace(
            "calibration_nightly: 17290 passed", f"calibration_nightly: 17290 passed {PEP_LOG}"
        )
    )
    compute(fx)
    md = md_of(fx)
    for pep in (PEP_STATUS, PEP_REASON, PEP_LOG):
        assert pep.lower() not in md.lower(), (pep, [ln for ln in md.splitlines() if pep.lower()
                                                     in ln.lower()])  # fmt: skip
    assert "redacted" in md.lower()


def test_item7_file_names_keep_their_case(fx: Fixture) -> None:
    compute(fx)
    md = md_of(fx)
    assert "CHANGELOG.md" in md, [ln for ln in md.splitlines() if "changel" in ln.lower()]
    assert "LICENSE" in md
    assert "changelOG" not in md


def test_item7_render_redacts_notes_and_values(fx: Fixture) -> None:
    compute(fx)
    doc = json.loads((fx.out / "readiness.json").read_text(encoding="utf-8"))
    doc["notes"] = [f"gh api failed: {PEP_STATUS}"]
    for r in doc["rows"]:
        if r["id"] == "R5":
            r["value"] = f"TE3 fail ({PEP_REASON})"
    md = load_module().render(doc)
    assert PEP_STATUS.lower() not in md.lower() and PEP_REASON.lower() not in md.lower()
    assert "CHANGELOG.md" in md


def _mini_doc(value: str, note: str) -> dict[str, Any]:
    r = {"id": "R5", "criterion": "TE", "status": "fail", "value": value, "threshold": "",
         "evidence": ["CHANGELOG.md", "LICENSE"], "waiver": None, "source": "both"}  # fmt: skip
    return {"version": "1.0.0", "commit": "a" * 40, "computed": "2026-11-01T00:00:00+00:00",
            "overall": "fail", "rows": [r], "notes": [note]}  # fmt: skip


@settings(max_examples=60, deadline=None)
@given(st.text(alphabet=AA, min_size=7, max_size=30))
def test_item7_property_any_peptide_redacted(pep: str) -> None:
    assume(not any(w in pep for w in ("CHANGELOG", "LICENSE", "READINESS")))
    render = load_module().render
    base = render(_mini_doc("TE3 fail ()", "x")).lower()
    assume(pep.lower() not in base)
    md = render(_mini_doc(f"TE3 fail ({pep})", f"summary {pep};"))
    assert pep.lower() not in md.lower()
    assert "CHANGELOG.md" in md and "LICENSE" in md


# ------------------------------------------------------------------ item 8: H6 conclusions


def test_item8_skipped_counts_as_pass_and_is_listed(fx: Fixture) -> None:
    write_runs(fx, {**ALL_SUCCESS, "guards": "skipped", "summary": "skipped"})
    _, doc = fx.compute()
    h6 = row(doc, "H6")
    assert h6["status"] == "pass", h6
    assert "guards" in h6["value"] and "summary" in h6["value"], h6["value"]


@pytest.mark.parametrize(
    "conclusion", ["neutral", "cancelled", "timed_out", "failure", "action_required"]
)
def test_item8_bad_conclusions_fail(fx: Fixture, conclusion: str) -> None:
    write_runs(fx, {**ALL_SUCCESS, "tests (3.12)": conclusion})
    proc, doc = fx.compute()
    h6 = row(doc, "H6")
    assert h6["status"] == "fail" and "tests (3.12)" in text_of(h6), h6
    assert proc.returncode == 1


def test_item8_running_check_is_pending(fx: Fixture) -> None:
    write_runs(fx, {**ALL_SUCCESS, "tests (3.12)": None})
    proc, doc = fx.compute()
    h6 = row(doc, "H6")
    assert h6["status"] == "pending", h6
    assert proc.returncode == 3, out_of(proc)


def test_item8_failure_beats_running(fx: Fixture) -> None:
    write_runs(fx, {**ALL_SUCCESS, "tests (3.12)": None, "lint": "failure"})
    _, doc = fx.compute()
    assert row(doc, "H6")["status"] == "fail"


def test_item8_check_runs_for_another_commit_fail(fx: Fixture) -> None:
    write_runs(fx, sha=fx.x)  # a saved response for X, not for HEAD = R
    proc, doc = fx.compute()
    h6 = row(doc, "H6")
    assert h6["status"] == "fail", h6
    assert proc.returncode == 1


def test_item8_check_runs_stale_after_new_commit(fx: Fixture) -> None:
    write_runs(fx)  # for R
    fx.write_text("README.md", "# GIPPA 1.0\n")
    fx.commit("readme")
    _, doc = fx.compute()  # the file is not rewritten: its head_sha is R, HEAD is newer
    assert row(doc, "H6")["status"] == "fail"


# ------------------------------------------------------------------ item 9: every container =

DIGEST = "sha256:" + "0" * 64


@pytest.mark.parametrize(
    "extra",
    [
        "process.container = 'docker://ubuntu:latest'",
        "process {\n    container = 'docker://ubuntu:latest'\n}",
    ],
)
def test_item9_container_outside_with_name_checked(fx: Fixture, extra: str) -> None:
    p = fx.repo / "workflow/conf/containers.config"
    p.write_text(p.read_text(encoding="utf-8") + "\n" + extra + "\n", encoding="utf-8")
    _, doc = compute(fx)
    h4 = row(doc, "H4")
    assert h4["status"] == "fail" and "ubuntu" in text_of(h4), h4


@pytest.mark.parametrize(
    "extra",
    [
        "// process.container = 'docker://ubuntu:latest'",
        f"process.container = 'ubuntu@{DIGEST}'",
    ],
)
def test_item9_commented_or_digest_container_passes(fx: Fixture, extra: str) -> None:
    p = fx.repo / "workflow/conf/containers.config"
    p.write_text(p.read_text(encoding="utf-8") + "\n" + extra + "\n", encoding="utf-8")
    _, doc = compute(fx)
    assert row(doc, "H4")["status"] == "pass", row(doc, "H4")


# ------------------------------------------------------------------ item 10: R6 paths, TE8 rows


@pytest.mark.parametrize(
    "path", ["workflow/main.nf", "tests/integration/test_x.py", "tests/contract/test_x.py"]
)
def test_item10_r6_integration_paths_invalidate(fx: Fixture, path: str) -> None:
    fx.write_text(path, "// changed\n")
    fx.commit("integration-tree change")
    _, doc = compute(fx)
    r6 = row(doc, "R6")
    assert r6["status"] == "pending", r6


K12_DIED = "k12_mg1655: gate 30 passed; gippa run exit 1 ended during gippa run (exit 1)"


@pytest.mark.parametrize(("exit_", "summary"), [(1, K12_DIED), (0, "k12_mg1655: gate 30 passed;")])
def test_item10_newest_k12_row_without_te8_beats_older_pass(
    fx: Fixture, exit_: int, summary: str
) -> None:
    c = fx.short(fx.r)
    fx.set_weekly(
        passing_weekly(c, date="2026-11-01", job="3000")
        + nightly_row("2026-11-08", "3100_0", c, exit_, summary)
    )
    _, doc = compute(fx)
    r5 = row(doc, "R5")
    assert r5["status"] in ("fail", "pending"), r5
    assert "TE8 pass" not in r5["value"], r5["value"]


def test_item10_nissle_row_does_not_decide_te8(fx: Fixture) -> None:
    c = fx.short(fx.r)
    fx.set_weekly(
        passing_weekly(c, date="2026-11-01", job="3000")
        + nightly_row("2026-11-08", "3100_1", c, 1, "nissle1917: gate 30 passed; TE2 x.json;")
    )
    _, doc = compute(fx)
    assert "TE8 pass" in row(doc, "R5")["value"]


def test_item10_te8_uses_ok_not_exit_only(fx: Fixture) -> None:
    c = fx.short(fx.r)
    fx.set_weekly(
        nightly_row("2026-11-01", "3000_0", c, 0,
                    "k12_mg1655: gate 2 failed, 28 passed; gippa run exit 0; TE8 x.json;")
    )  # fmt: skip
    _, doc = compute(fx)
    r5 = row(doc, "R5")
    assert r5["status"] == "fail" and "TE8 fail" in r5["value"], r5


# ------------------------------------------------------------------ item 11 / decision (f): R3


def test_decision_f_no_calibration_full_row_is_pending(fx: Fixture) -> None:
    c = fx.short(fx.r)
    fx.set_nightly(
        nightly_row("2026-11-02", "2000_0", c, 0,
                    "integration: 812 passed, 92 skipped in 1438.86s (0:23:58); nf-test "
                    "farmshare_local: SUCCESS: Executed 40 tests")
        + nightly_row("2026-11-02", "2000_1", c, 0,
                      "benchmarks: 57 passed in 519.92s (0:08:39) (wheat wheat-02)")
        + nightly_row("2026-11-02", "2000_2", c, 0,
                      "calibration_nightly: 17290 passed, 108 skipped, 13 xfailed in 20556.28s")
        + nightly_row("2026-11-02", "2000_3", c, 0, "gippa run sim100k: exit 0")
    )  # fmt: skip
    _, doc = compute(fx)
    r3, r4 = row(doc, "R3"), row(doc, "R4")
    assert r3["status"] == "pending" and "calibration_full" in r3["value"], r3
    assert r4["status"] == "pending", r4


# ------------------------------------------------------------------ item 12: surviving mutants


def test_item12_waiver_applies_to_pending_and_keeps_missing_value(fx: Fixture) -> None:
    fx.write_json("benchmarks/baseline.json", {})
    fx.write_text("docs/adr/ADR-0060-runtime-waiver.md", adr_text(60, "accepted-by-khoa", "R7"))
    _, doc = compute(fx)
    r7 = row(doc, "R7")
    assert r7["status"] == "waived" and r7["waiver"] == "ADR-0060", r7
    assert r7["value"].startswith("missing:"), r7["value"]


@pytest.mark.parametrize(("rid", "status"), [("R1", "pass"), ("R10", "deferred")])
def test_item12_waiver_never_touches_pass_or_deferred(fx: Fixture, rid: str, status: str) -> None:
    fx.write_text("docs/adr/ADR-0061-waiver.md", adr_text(61, "accepted-by-khoa", rid))
    _, doc = compute(fx)
    r = row(doc, rid)
    assert r["status"] == status and r["waiver"] is None, r


def test_item12_waiver_for_r10_does_not_waive_r1(fx: Fixture) -> None:
    fx.write_json(METRICS, metrics_doc(share=0.9949))
    fx.write_text("docs/adr/ADR-0062-waiver.md", adr_text(62, "accepted-by-khoa", "R10"))
    _, doc = compute(fx)
    r1 = row(doc, "R1")
    assert r1["status"] == "fail" and r1["waiver"] is None, r1


def test_item12_newest_uses_job_id_not_file_order(fx: Fixture) -> None:
    c = fx.short(fx.r)
    fx.set_nightly(
        passing_nightly(c, date="2026-11-02", job="2009")
        + nightly_row("2026-11-02", "2005_3", c, 1, "calibration_full: 1 failed, 3119 passed")
    )
    _, doc = compute(fx)
    assert row(doc, "R3")["status"] == "pass"


def test_item12_h1_needs_gate_m10(fx: Fixture) -> None:
    fx.write_text("PROGRESS.md", progress_md(range(1, 10)))
    _, doc = compute(fx)
    h1 = row(doc, "H1")
    assert h1["status"] == "fail" and "M10" in h1["value"], h1


@pytest.fixture
def computed(fx: Fixture) -> Fixture:
    proc, _ = compute(fx)
    assert proc.returncode == 0, out_of(proc)
    assert fx.verify().returncode == 0
    return fx


def _tamper(fx: Fixture, edit: Any) -> Path:
    p = fx.out / "readiness.json"
    doc = json.loads(p.read_text(encoding="utf-8"))
    doc["rows"] = edit(doc["rows"])
    p.write_text(json.dumps(doc), encoding="utf-8")
    return p


def test_item12_verify_rejects_a_dropped_row(computed: Fixture) -> None:
    p = _tamper(computed, lambda rows: [r for r in rows if r["id"] != "H6"])
    assert computed.verify(p).returncode == 1


def test_item12_verify_rejects_reordered_rows(computed: Fixture) -> None:
    def swap(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
        i, j = (next(k for k, r in enumerate(rows) if r["id"] == x) for x in ("R7", "R8"))
        rows[i], rows[j] = rows[j], rows[i]
        return rows

    assert computed.verify(_tamper(computed, swap)).returncode == 1


# ------------------------------------------------------------------ item 13: sub-issue gates


def test_item13_sub_issue_gate_is_not_milestone_gate(fx: Fixture) -> None:
    prog = progress_md([k for k in range(1, 11) if k != 9]) + "\n### Gate M9-03\n- Evidence: x\n"
    fx.write_text("PROGRESS.md", prog)
    _, doc = compute(fx)
    h1 = row(doc, "H1")
    assert h1["status"] == "fail" and "missing: M9" in h1["value"], h1


def test_item13_r9_needs_gate_m10_not_m10_15(fx: Fixture) -> None:
    prog = progress_md(range(1, 10)) + "\n### Gate M10-15\n- Khoa review: approved 2026-11-01\n"
    fx.write_text("PROGRESS.md", prog)
    _, doc = compute(fx)
    assert row(doc, "R9")["status"] == "pending", row(doc, "R9")
    assert "missing: M10" in row(doc, "H1")["value"]
