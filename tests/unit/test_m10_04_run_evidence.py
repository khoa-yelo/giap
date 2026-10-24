"""``validation/pilot/run_evidence.py``: rule sets A (arm) and C (compare) on hand-made run
directories (docs/specs/M10-04-pilot-runs.md, acceptance criteria 1-4, examples 3, 6 and 7, edge
cases "Empty input", "Head job killed before the tail", "Detection tier differs").

Fixtures are built in ``tmp_path`` by ``helpers.pilot_runs_m10_04``: manifests come from the real
``build_manifest`` and validate against ``run_manifest.schema.json``; ``sacct`` rows are passed as
data to the pure functions and served by the fake ``sacct`` (``helpers.slurm_sandbox``) to the CLI.
The compare fixtures are the output of the real ``giap compare`` on two tiny arms (criterion 4).

Pinned by the test-author (spec Q6, lead decision 1): ``giap compare`` (M9-22) writes one
``compare_stability`` row per schema family even when neither run has that table (``n_a``/``n_b``
NA, flags ``missing_in:A1;missing_in:A2``), so C4 is "one row per enum family". Failure strings
whose suffix the spec leaves open (``manifest_invalid:<path>``, ``genes_table:<message>``,
``table:<name>:<message>``) are matched by prefix; ``commit_mismatch:<sha>`` carries the
manifest's (observed) commit.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import pytest
from helpers.manifest_data import validate
from helpers.pilot_runs_m10_04 import (
    COMMIT,
    HEAD_JOB,
    RUN_EVIDENCE,
    SHA_ANNOTATION,
    SHA_GENOME,
    SHA_PEPTIDES,
    TABLES,
    TASK_JOBS,
    arm_genes,
    bakta_task,
    default_trace,
    edit_manifest,
    expect,
    giap_compare,
    load_run_evidence,
    sacct_rows,
    sha256,
    trace_task,
    two_arms,
    write_arm,
)
from helpers.slurm_sandbox import FAKE_SACCT, write_exe
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

# ------------------------------------------------------------------ helpers


def _arm(out: Path, sacct: list[dict[str, str]] | None = None, **exp: Any) -> dict[str, Any]:
    mod = load_run_evidence()
    doc = mod.arm_evidence(out, sacct_rows() if sacct is None else sacct, expect(**exp))
    assert isinstance(doc, dict)
    json.dumps(doc)  # JSON-serialisable
    return doc


def _fake_sacct(tmp: Path, rows: list[dict[str, str]]) -> dict[str, str]:
    bindir, state = tmp / "fakebin", tmp / "fakestate"
    bindir.mkdir(exist_ok=True)
    state.mkdir(exist_ok=True)
    write_exe(bindir / "sacct", FAKE_SACCT.replace("{python}", sys.executable))
    (state / "sacct_rows.json").write_text(json.dumps(rows))
    return {"PATH": f"{bindir}:/usr/bin:/bin", "HOME": str(tmp), "LANG": "C.UTF-8",
            "FAKE_SLURM_DIR": str(state), "PYTHONPATH": str(RUN_EVIDENCE.parents[2] / "src"),
            "GIAP_ROOT": str(tmp / "no_root"), "PYTHONDONTWRITEBYTECODE": "1"}  # fmt: skip


def _cli(tmp: Path, *args: str, rows: list[dict[str, str]] | None = None):  # type: ignore[no-untyped-def]
    assert RUN_EVIDENCE.is_file(), "validation/pilot/run_evidence.py does not exist (M10-04)"
    env = _fake_sacct(tmp, sacct_rows() if rows is None else rows)
    return subprocess.run([sys.executable, str(RUN_EVIDENCE), *args], capture_output=True,
                          text=True, env=env, cwd=str(tmp), timeout=300,
                          check=False)  # fmt: skip


def _arm_cli(tmp: Path, out: Path, rows: list[dict[str, str]] | None = None, **kw: str):  # type: ignore[no-untyped-def]
    dest = tmp / "evidence" / "a2.json"
    args = ["arm", str(out), "--arm", kw.get("arm", "A2"), "--job", HEAD_JOB,
            "--expect-commit", kw.get("commit", COMMIT),
            "--expect-input", f"genome={SHA_GENOME}",
            "--expect-input", f"annotation={SHA_ANNOTATION}",
            "--expect-input", f"peptides={kw.get('peptides', SHA_PEPTIDES)}",
            "--out", str(dest)]  # fmt: skip
    return _cli(tmp, *args, rows=rows), dest


# ------------------------------------------------------------------ criterion 1: a passing arm


def test_fixture_manifest_is_schema_valid(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2")
    validate(json.loads((out / "run_manifest.json").read_text()))
    assert _arm(out)["failures"] == []


def test_passing_arm_lists_the_fields(tmp_path: Path) -> None:
    """Example 3 / criterion 1: ``failures: []`` and every listed field."""
    out = write_arm(tmp_path / "out_a2")
    doc = _arm(out)
    assert doc["failures"] == []
    assert doc["arm"] == "A2"
    assert Path(doc["outdir"]) == out
    for key in ("job", "manifest", "inputs", "report", "tables", "tasks"):
        assert key in doc, key
    assert str(doc["job"]["id"]) == HEAD_JOB
    assert doc["job"]["state"] == "COMPLETED"
    assert doc["job"]["exit_code"] == "0:0"
    man = doc["manifest"]
    for key in ("status", "exit_code", "failed_step", "failed_log", "git_commit", "git_dirty",
                "profile", "resume", "started_at", "finished_at", "seeds", "modules_off",
                "detection_tier", "annotation", "databases", "warnings"):  # fmt: skip
        assert key in man, key
    assert man["status"] == "completed" and man["git_commit"] == COMMIT
    assert man["seeds"]["master"] == 20261014
    assert man["modules_off"] == []
    assert doc["inputs"] == {"genome": SHA_GENOME, "annotation": SHA_ANNOTATION,
                             "peptides": SHA_PEPTIDES}  # fmt: skip


def test_report_and_tables_are_hashed_not_read(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2")
    doc = _arm(out)
    html = out / "report" / "report.html"
    assert doc["report"]["sha256"] == sha256(html)
    assert doc["report"]["bytes"] == html.stat().st_size
    assert Path(doc["report"]["path"]) in (html, html.relative_to(out))
    genes = out / "tables" / "genes.tsv"
    entry = doc["tables"].get("genes.tsv", doc["tables"].get("genes"))
    assert entry is not None, doc["tables"]
    assert entry["sha256"] == sha256(genes) and entry["bytes"] == genes.stat().st_size
    text = json.dumps(doc)
    assert "TST_00001" not in text, "no table content in the evidence (Khoa's rule)"


def test_task_summary_from_trace_and_sacct(tmp_path: Path) -> None:
    trace = [*default_trace(), trace_task("GIAP:CORE:DETECTION (detection)", "1900004",
                                          attempt=2, task_id=3)]  # fmt: skip
    out = write_arm(tmp_path / "out_a2", trace=trace)
    rows = sacct_rows(tasks={TASK_JOBS[0]: 14, TASK_JOBS[1]: 8, "1900004": 12})
    tasks = _arm(out, rows)["tasks"]
    assert tasks["count"] == 3
    assert tasks["failed"] == 0
    assert tasks["retried"] == 1
    assert tasks["max_alloc_cpus"] == 14


def test_warnings_are_listed_never_failed(tmp_path: Path) -> None:
    ev = {"ts": "2026-10-24T09:00:00+00:00", "step": "operons", "level": "WARNING",
          "event": "operons.layer_missing",
          "data": {"flag": "layer_missing:terminators"}}  # fmt: skip
    out = write_arm(tmp_path / "out_a2", events=[ev])
    doc = _arm(out)
    assert doc["failures"] == []
    assert doc["manifest"]["warnings"], "the manifest warnings are copied"


def test_detection_tier_is_recorded_not_failed(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a1", tier="sequences_only")
    doc = _arm(out, arm="A1")
    assert doc["failures"] == []
    assert doc["manifest"]["detection_tier"] == "sequences_only"


def test_cli_passing_arm_writes_json_and_prints_a_row(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2")
    res, dest = _arm_cli(tmp_path, out)
    assert res.returncode == 0, res.stdout + res.stderr
    doc = json.loads(dest.read_text())
    assert doc["failures"] == [] and doc["arm"] == "A2"
    assert str(doc["job"]["id"]) == HEAD_JOB and doc["tasks"]["max_alloc_cpus"] == 14
    rows = [ln for ln in res.stdout.splitlines() if ln.startswith("|")]
    assert len(rows) == 1 and "A2" in rows[0], res.stdout


# ------------------------------------------------------------------ criterion 2: rules A1-A9


def _break(case: str, tmp: Path) -> tuple[Path, list[dict[str, str]], dict[str, Any]]:
    out = tmp / "out_a2"
    rows: list[dict[str, str]] = sacct_rows()
    exp: dict[str, Any] = {}
    if case == "no_manifest":
        write_arm(out, manifest=False)
        rows = sacct_rows(state="TIMEOUT", exit_code="0:15")
    elif case == "run_failed":
        trace = [*default_trace(), trace_task("GIAP:OPERONS:OPERONS_STEP (operons)", "1900005",
                                              status="FAILED", task_id=5)]  # fmt: skip
        write_arm(out, trace=trace, exit_code=4, failed_step="operons")
        rows = sacct_rows(state="FAILED", exit_code="4:0")
    elif case == "commit_mismatch":
        write_arm(out, commit="deadbee")
    elif case == "dirty_checkout":
        write_arm(out, dirty=True)
    elif case == "live_bakta":
        write_arm(out, trace=[*default_trace(), bakta_task()])
    elif case == "input_mismatch":
        write_arm(out)
        exp = {"inputs": {"genome": SHA_GENOME, "annotation": SHA_ANNOTATION,
                          "peptides": "f" * 64}}  # fmt: skip
    elif case == "cpu_cap":
        write_arm(out)
        rows = sacct_rows(tasks={TASK_JOBS[0]: 14, TASK_JOBS[1]: 16})
    elif case == "no_report":
        write_arm(out, report=None)
    elif case == "genes_table":
        write_arm(out)
        from helpers.viz_out_tree import write_schema_table

        write_schema_table(out / "tables" / "genes.tsv", "genes", arm_genes(), drop=["strand"])
    elif case == "module_off":
        write_arm(out, modules=("core", "mobile", "bgc"))
    else:  # pragma: no cover
        raise AssertionError(case)
    return out, rows, exp


#: rule → (fixture case, exact failure string or prefix ending in ":")
RULES = {
    "A1": ("no_manifest", "no_manifest"),
    "A2": ("run_failed", "run_failed:operons"),
    "A3-commit": ("commit_mismatch", "commit_mismatch:deadbee"),
    "A3-dirty": ("dirty_checkout", "dirty_checkout"),
    "A4": ("live_bakta", "live_bakta"),
    "A5": ("input_mismatch", "input_mismatch:peptides"),
    "A6": ("cpu_cap", "cpu_cap:16"),
    "A7": ("no_report", "no_report"),
    "A8": ("genes_table", "genes_table:"),
    "A9": ("module_off", "module_off:function"),
}


def _has(failures: list[str], want: str) -> bool:
    return any(f.startswith(want) for f in failures) if want.endswith(":") else want in failures


@pytest.mark.parametrize("rule", sorted(RULES))
def test_each_rule_fails_on_its_fixture(tmp_path: Path, rule: str) -> None:
    case, want = RULES[rule]
    out, rows, exp = _break(case, tmp_path)
    failures = _arm(out, rows, **exp)["failures"]
    assert _has(failures, want), (rule, failures)
    if case != "no_manifest":  # a killed head job also lacks the tail's outputs
        assert len(failures) == 1, (rule, failures)


@pytest.mark.parametrize("rule", ["A2", "A3-commit", "A4", "A6", "A9"])
def test_cli_exits_1_on_a_broken_arm(tmp_path: Path, rule: str) -> None:
    case, want = RULES[rule]
    out, rows, exp = _break(case, tmp_path)
    assert not exp
    res, dest = _arm_cli(tmp_path, out, rows=rows)
    assert res.returncode == 1, res.stdout + res.stderr
    assert _has(json.loads(dest.read_text())["failures"], want)


def test_cli_exits_1_on_a_wrong_peptides_sha(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2")
    res, dest = _arm_cli(tmp_path, out, peptides="f" * 64)
    assert res.returncode == 1, res.stdout + res.stderr
    assert json.loads(dest.read_text())["failures"] == ["input_mismatch:peptides"]


def test_cli_exits_1_without_a_manifest(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2", manifest=False)
    res, dest = _arm_cli(tmp_path, out, rows=sacct_rows(state="TIMEOUT", exit_code="0:15"))
    assert res.returncode == 1, res.stdout + res.stderr
    assert "no_manifest" in json.loads(dest.read_text())["failures"]


def test_manifest_invalid_names_the_problem(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2")
    edit_manifest(out, status=None)
    failures = _arm(out)["failures"]
    assert any(f.startswith("manifest_invalid:") for f in failures), failures


def test_head_job_not_completed_fails_a2(tmp_path: Path) -> None:
    """A2 also needs the head job ``COMPLETED`` with ``0:0`` (not only the manifest)."""
    out = write_arm(tmp_path / "out_a2")
    failures = _arm(out, sacct_rows(state="FAILED", exit_code="1:0"))["failures"]
    assert any(f.startswith("run_failed:") for f in failures), failures


def test_light_bakta_db_is_a_live_annotation(tmp_path: Path) -> None:
    meta = {"annotation_source": "bakta", "bakta_version": "1.12.2", "db_version": "6.0",
            "db_type": "light"}  # fmt: skip
    out = write_arm(tmp_path / "out_a2", annotation_meta=meta)
    assert _arm(out)["failures"] == ["live_bakta"]


def test_missing_genes_table_fails_a8(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2")
    (out / "tables" / "genes.tsv").unlink()
    failures = _arm(out)["failures"]
    assert any(f.startswith("genes_table:") for f in failures), failures


def test_empty_report_fails_a7(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2", report="")
    assert _arm(out)["failures"] == ["no_report"]


_INDEPENDENT = {
    "commit": "commit_mismatch:deadbee",
    "dirty": "dirty_checkout",
    "cpu": "cpu_cap:16",
    "report": "no_report",
    "module": "module_off:function",
    "peptides": "input_mismatch:peptides",
}


@settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(st.sets(st.sampled_from(sorted(_INDEPENDENT))))
def test_rules_are_independent(breaks: set[str]) -> None:
    """Property: the failures are exactly the strings of the rules broken, one per rule."""
    with tempfile.TemporaryDirectory() as d:
        kw: dict[str, Any] = {}
        if "commit" in breaks:
            kw["commit"] = "deadbee"
        if "dirty" in breaks:
            kw["dirty"] = True
        if "report" in breaks:
            kw["report"] = None
        if "module" in breaks:
            kw["modules"] = ("core", "mobile", "bgc")
        out = write_arm(Path(d) / "out_a2", **kw)
        rows = sacct_rows(tasks={TASK_JOBS[0]: 16 if "cpu" in breaks else 14, TASK_JOBS[1]: 2})
        exp = {}
        if "peptides" in breaks:
            exp["inputs"] = {"genome": SHA_GENOME, "annotation": SHA_ANNOTATION,
                             "peptides": "f" * 64}  # fmt: skip
        failures = _arm(out, rows, **exp)["failures"]
    assert sorted(failures) == sorted(_INDEPENDENT[b] for b in breaks)


# ------------------------------------------------------------------ criterion 3: empty results


def test_empty_results_still_pass(tmp_path: Path) -> None:
    """All genes ``not_detected`` and a header-only ``novel_orfs.tsv``: the pipeline ran."""
    out = write_arm(tmp_path / "out_a1", genes=arm_genes(detected=False), novel_orfs_empty=True)
    assert _arm(out, arm="A1")["failures"] == []


# ------------------------------------------------------------------ criterion 4: rule set C


@pytest.fixture(scope="module")
def compare_pair(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    """``CMP`` and ``CMP2`` from the real ``giap compare`` on two tiny arms (example 7)."""
    tmp = tmp_path_factory.mktemp("m10_04_cmp")
    a1, a2 = two_arms(tmp)
    for name in ("cmp_a1_a2", "cmp2"):
        res = giap_compare(a1, a2, tmp / name, env_root=tmp / "env")
        assert res.returncode == 0, res.stdout[-2000:] + res.stderr[-2000:]
    return tmp / "cmp_a1_a2", tmp / "cmp2"


@pytest.fixture
def pair(compare_pair: tuple[Path, Path], tmp_path: Path) -> tuple[Path, Path]:
    cmp, cmp2 = compare_pair
    shutil.copytree(cmp, tmp_path / "cmp_a1_a2")
    shutil.copytree(cmp2, tmp_path / "cmp2")
    return tmp_path / "cmp_a1_a2", tmp_path / "cmp2"


def _compare(cmp: Path, cmp2: Path) -> dict[str, Any]:
    mod = load_run_evidence()
    doc = mod.compare_evidence(cmp, cmp2, ("A1", "A2"))
    assert isinstance(doc, dict)
    json.dumps(doc)
    return doc


def _rows(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    return lines[0].split("\t"), [ln.split("\t") for ln in lines[1:]]


def _write_rows(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.write_text("\n".join("\t".join(r) for r in [header, *rows]) + "\n", encoding="utf-8")


def test_compare_passes_and_holds_hashes_only(pair: tuple[Path, Path]) -> None:
    cmp, cmp2 = pair
    doc = _compare(cmp, cmp2)
    assert doc["failures"] == []
    text = json.dumps(doc)
    for t in TABLES:
        assert sha256(cmp / "tables" / f"{t}.tsv") in text, t
    assert sha256(cmp / "report" / "compare.html") in text
    assert str((cmp / "report" / "compare.html").stat().st_size) in text
    header, rows = _rows(cmp / "tables" / "compare_abundance.tsv")
    rho = rows[0][header.index("spearman_rho")]
    assert rho not in text, "no value from the tables (no rho) in compare.json"


def test_c1_a_missing_column(pair: tuple[Path, Path]) -> None:
    cmp, cmp2 = pair
    path = cmp / "tables" / "compare_abundance.tsv"
    header, rows = _rows(path)
    j = header.index("method")
    _write_rows(path, header[:j] + header[j + 1 :], [r[:j] + r[j + 1 :] for r in rows])
    failures = _compare(cmp, cmp2)["failures"]
    assert any(f.startswith("table:compare_abundance") for f in failures), failures


def test_c2_two_abundance_rows(pair: tuple[Path, Path]) -> None:
    cmp, cmp2 = pair
    path = cmp / "tables" / "compare_abundance.tsv"
    header, rows = _rows(path)
    _write_rows(path, header, [rows[0], rows[0]])
    assert "abundance_rows:2" in _compare(cmp, cmp2)["failures"]


def test_c3_a_matrix_row_short(pair: tuple[Path, Path]) -> None:
    cmp, cmp2 = pair
    path = cmp / "tables" / "compare_detection.tsv"
    header, rows = _rows(path)
    sec = header.index("section")
    matrix = [i for i, r in enumerate(rows) if r[sec] == "matrix"]
    drop = matrix[-1]
    _write_rows(path, header, [r for i, r in enumerate(rows) if i != drop])
    failures = _compare(cmp, cmp2)["failures"]
    assert f"matrix_rows:{len(matrix) - 1}" in failures, failures


def _pin_m9_22_stability(cmp: Path) -> None:
    """Q6: a family with no table in either run is a row (NA counts, ``missing_in`` flags)."""
    with (cmp / "tables" / "compare_stability.tsv").open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert [r["family"] for r in rows] == [
        "islands", "operon_flags", "enriched_categories", "novel_orfs", "alt_sites",
    ]  # fmt: skip
    for r in rows:
        assert (r["run_a"], r["run_b"]) == ("A1", "A2")
        assert set(r["flags"].split(";")) == {"missing_in:A1", "missing_in:A2"}


def test_c4_a_family_missing(pair: tuple[Path, Path]) -> None:
    """C4 is "one row per enum family", as M9-22 writes it (pinned first, spec Q6)."""
    cmp, cmp2 = pair
    _pin_m9_22_stability(cmp)
    path = cmp / "tables" / "compare_stability.tsv"
    header, rows = _rows(path)
    fam = header.index("family")
    _write_rows(path, header, [r for r in rows if r[fam] != "alt_sites"])
    assert "stability_family:alt_sites" in _compare(cmp, cmp2)["failures"]


@pytest.mark.parametrize("victim", ["report/compare.html", "plots/V16_compare.svg"])
def test_c5_no_compare_report(pair: tuple[Path, Path], victim: str) -> None:
    cmp, cmp2 = pair
    (cmp / victim).unlink()
    assert _compare(cmp, cmp2)["failures"] == ["no_compare_report"]


def test_c5_empty_compare_report(pair: tuple[Path, Path]) -> None:
    cmp, cmp2 = pair
    (cmp / "report" / "compare.html").write_bytes(b"")
    assert _compare(cmp, cmp2)["failures"] == ["no_compare_report"]


def test_c6_one_byte_changed_in_cmp2(pair: tuple[Path, Path]) -> None:
    cmp, cmp2 = pair
    path = cmp2 / "tables" / "compare_detection.tsv"
    data = bytearray(path.read_bytes())
    i = data.index(b"detected\t")  # a value cell, not the header
    i = data.index(b"detected", i + 1)
    data[i] = ord("D")
    path.write_bytes(bytes(data))
    failures = _compare(cmp, cmp2)["failures"]
    assert len(failures) == 1 and failures[0].startswith("nondeterministic:compare_detection")


def test_compare_on_empty_arms_passes(tmp_path: Path) -> None:
    """Edge case "Empty input": all ``not_detected``, the rho row NA plus its flag; C1 holds."""
    a1, a2 = two_arms(tmp_path, detected=False)
    for name in ("cmp_a1_a2", "cmp2"):
        res = giap_compare(a1, a2, tmp_path / name, env_root=tmp_path / "env")
        assert res.returncode == 0, res.stderr[-2000:]
    assert _compare(tmp_path / "cmp_a1_a2", tmp_path / "cmp2")["failures"] == []


def test_compare_cli_exit_codes(pair: tuple[Path, Path], tmp_path: Path) -> None:
    cmp, cmp2 = pair
    dest = tmp_path / "evidence" / "compare.json"
    args = ["compare", str(cmp), "--second", str(cmp2), "--labels", "A1,A2", "--out", str(dest)]
    res = _cli(tmp_path, *args)
    assert res.returncode == 0, res.stdout + res.stderr
    assert json.loads(dest.read_text())["failures"] == []
    (cmp / "report" / "compare.html").unlink()
    res = _cli(tmp_path, *args)
    assert res.returncode == 1, res.stdout + res.stderr
    assert json.loads(dest.read_text())["failures"] == ["no_compare_report"]
