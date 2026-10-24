"""``validation/pilot/metrics.py compute`` / ``render`` on synthetic run directories (M10-05
acceptance criteria 5-8; example 3 and the edge cases "missing table" and "--out inside the
repository").

Spec: docs/specs/M10-05-pilot-metrics.md ("Inputs", "Outputs", "REPORT.md draft", "§15.6 checklist
in pass 1", "Interface"). Inputs are built in ``tmp`` by ``tests/helpers/pilot_metrics_m10.py``
(a giapsim SIM100K genome and Bakta JSON; every table under its schema; a schema-valid
``run_manifest.json``); nothing comes from ``$GIAP_ROOT`` or a DIA-NN report.

Interface fixed here (beyond the spec's command lines and exit codes):

* ``metrics.json`` per-arm objects are keyed by the arm label: ``m1.<arm>`` (the
  ``identifications`` dict, so the template path ``m1.A2.union.peptides.non_contaminant``
  exists), ``m2.<arm>`` (``mapping_shares``; ``m2.A1.genome_inconsistent`` adds ``consistent``),
  ``m3.<arm>``, ``m7.<arm>``;
* a ``warnings`` entry is a string or a dict with ``code``; the metric-1 cross-check code is
  ``raw_vs_written_q01_out_of_range``;
* a metric whose table is missing is ``null`` or a dict with ``reason = "missing:<file>"``;
* ``checklist`` is a list in the spec's table order (bullet 1 = runs, 2 = TE1, 3 = compare,
  4 = runtimes, 5 = TE2-TE5); each entry a status string or a dict with ``status``;
* review round (lead decisions): exit 1 only for ``Inconsistent`` (a cross-check or a rejected
  notes file; a template placeholder missing under a metric that was computed); every read,
  schema, JSON, missing-key and config error is a precondition, exit 2, nothing written (B1). A
  metric with a ``reason`` renders its placeholders as ``**missing: <reason>**`` and ``render``
  exits 0 (B2). ``metrics.json`` ``giap_commit`` is the runs' ``run_manifest.json``
  ``git_commit`` (S4). Warning codes: ``unknown_source``, ``target_q_na``,
  ``sacct_ignored_no_attempts``, and one naming the compare table that fails
  ``giap.io.output_checks.check_tsv`` (B3).
"""

from __future__ import annotations

import json
import re
import shutil
from pathlib import Path
from typing import Any

import pytest
from helpers.pilot_metrics_m10 import (
    GIAP_COMMIT,
    REPO,
    TEMPLATE,
    Pilot,
    World,
    build_pilot,
    load_metrics,
    output_files,
    pep_row,
    run_cli,
    sha256,
    sim_world,
    write_a3,
    write_peptides,
)

PLACEHOLDER = re.compile(r"\{\{\s*([^}\s]+)\s*\}\}")


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> World:
    return sim_world(tmp_path_factory.mktemp("m10_05_sim"))


@pytest.fixture
def pilot(tmp_path: Path, world: World) -> Pilot:
    return build_pilot(tmp_path / "pilot", world)


def _compute(p: Pilot, out: Path, *extra: str, **kw: Any) -> Any:
    return run_cli(*p.compute_args(out, *extra, **kw), cwd=p.root)


def _ok(res: Any) -> None:
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


def _metrics(out: Path) -> dict[str, Any]:
    doc: dict[str, Any] = json.loads((out / "metrics.json").read_text(encoding="utf-8"))
    return doc


def _codes(metrics: dict[str, Any]) -> list[str]:
    return [w if isinstance(w, str) else str(w.get("code")) for w in metrics.get("warnings", [])]


def _cl(entry: Any) -> str:
    return str(entry["status"] if isinstance(entry, dict) else entry)


def _nothing_written(out: Path) -> bool:
    return not out.exists() or not any(out.iterdir())


# ------------------------------------------------------------------ criterion 5


def test_compute_writes_the_three_files_and_is_byte_identical(pilot: Pilot, tmp_path: Path) -> None:
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out))
    first = output_files(out)
    assert {"metrics.json", "te1.json", "te1_exceptions.tsv"} <= set(first)
    _ok(_compute(pilot, out))
    assert output_files(out) == first  # same inputs, same bytes (no timestamp)


def test_metrics_json_layout(pilot: Pilot, tmp_path: Path) -> None:
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out))
    raw = (out / "metrics.json").read_text(encoding="utf-8")
    assert raw.endswith("\n")
    doc = json.loads(raw)
    assert raw == json.dumps(doc, sort_keys=True, indent=2) + "\n"
    assert doc["schema"] == "giap-pilot-metrics/1"
    for key in ("inputs", "spec", "giap_commit", "m1", "m2", "m3", "m4", "m5", "m7", "te1",
                "te1_a1_diagnostic", "h1", "h3", "checklist", "warnings"):  # fmt: skip
        assert key in doc, key
    assert re.fullmatch(r"[0-9a-f]{40}", str(doc["spec"])), doc["spec"]
    for name, rec in doc["inputs"].items():
        assert set(rec) >= {"path", "sha256"}, name
    te1 = json.loads((out / "te1.json").read_text(encoding="utf-8"))
    assert te1 == doc["te1"] and "share_exact_inframe" in te1


def test_compute_values_on_the_fixture(pilot: Pilot, tmp_path: Path) -> None:
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out))
    doc = _metrics(out)
    assert doc["m1"]["A2"]["union"]["peptides"]["all"] == 21  # 20 strain peptides + 1 Cont_
    assert doc["m1"]["A2"]["union"]["peptides"]["contaminant_only"] == 1
    m2 = doc["m2"]["A1"]
    assert (m2["n_universe"], m2["n_contaminant_only"], m2["n_denominator"]) == (21, 1, 20)
    assert m2["counts"] == {"exact": 17, "tolerant": 2, "unmapped": 1}
    gi = m2["genome_inconsistent"]
    assert gi["count"] == 3 and gi["consistent"] is True
    te1 = doc["te1"]
    assert te1["n_peptides"] == 20 and te1["exact_inframe"] == 20 and te1["passed"] is True
    assert "passed" not in doc["te1_a1_diagnostic"]  # a diagnostic, no pass flag
    assert doc["m4"]["spearman_rho"] == pytest.approx(0.83)
    assert doc["m5"]["islands"]["jaccard"] == pytest.approx(0.5)
    assert doc["h3"]["verdict"] == "supported"
    assert [s["name"] for s in doc["m7"]["A2"]["steps"]] == ["mapping", "bakta", "gene_table"]
    assert doc["m7"]["A2"]["wall_time_s"] == pytest.approx(2550.0)
    assert doc["m7"]["A2"]["hosts"] == ["sh04-01n01"]
    assert doc["m7"]["A2"]["runtime_complete"] is False


def test_no_sequence_in_metrics_json_and_exceptions_listed(tmp_path: Path, world: World) -> None:
    p = build_pilot(tmp_path / "pilot", world, te1_misses=2)
    out = tmp_path / "metrics"
    _ok(_compute(p, out))
    text = (out / "metrics.json").read_text(encoding="utf-8")
    assert not [pep for _, pep in p.peptides if pep in text]
    exc = (out / "te1_exceptions.tsv").read_text(encoding="utf-8").splitlines()
    header = exc[0].split("\t")
    assert header[:4] == ["kind", "locus_tag", "peptide", "hits"]
    assert len(exc) - 1 == 2
    assert {line.split("\t")[2] for line in exc[1:]} == {pep for _, pep in p.peptides[:2]}


# ------------------------------------------------------------------ criterion 6 / example 3


def test_example_3a_difference_within_crossings_no_warning(tmp_path: Path, world: World) -> None:
    p = build_pilot(tmp_path / "pilot", world, q_crossings=4, raw_vs_written=2)
    out = tmp_path / "metrics"
    _ok(_compute(p, out))
    assert "raw_vs_written_q01_out_of_range" not in _codes(_metrics(out))


def test_example_3b_difference_above_crossings_warns(tmp_path: Path, world: World) -> None:
    p = build_pilot(tmp_path / "pilot", world, q_crossings=1, raw_vs_written=2)
    out = tmp_path / "metrics"
    _ok(_compute(p, out))  # a warning, not an exit code
    assert "raw_vs_written_q01_out_of_range" in _codes(_metrics(out))


def test_zero_crossings_need_zero_difference(tmp_path: Path, world: World) -> None:
    p = build_pilot(tmp_path / "pilot", world, q_crossings=0, raw_vs_written=1)
    out = tmp_path / "metrics"
    _ok(_compute(p, out))
    assert "raw_vs_written_q01_out_of_range" in _codes(_metrics(out))


def test_example_3c_convert_json_disagrees_exit_1(tmp_path: Path, world: World) -> None:
    p = build_pilot(tmp_path / "pilot", world, conv_offset=-1)
    out = tmp_path / "metrics"
    res = _compute(p, out)
    assert res.returncode == 1, res.stdout + res.stderr


def test_run_exit_code_not_zero_exit_2(tmp_path: Path, world: World) -> None:
    p = build_pilot(tmp_path / "pilot", world, exit_a1=1)
    out = tmp_path / "metrics"
    res = _compute(p, out)
    assert res.returncode == 2 and _nothing_written(out), res.stdout + res.stderr
    assert "run_manifest.json" in res.stdout + res.stderr or "out_a1" in res.stdout + res.stderr


def test_genome_hashes_differ_exit_2(tmp_path: Path, world: World) -> None:
    p = build_pilot(tmp_path / "pilot", world, genome_sha_a1="e" * 64)
    out = tmp_path / "metrics"
    res = _compute(p, out)
    assert res.returncode == 2 and _nothing_written(out), res.stdout + res.stderr


def test_peptide_table_hash_differs_from_manifest_exit_2(pilot: Pilot, tmp_path: Path) -> None:
    before = sha256(pilot.pep_a2)
    write_peptides(pilot.pep_a2, [pep_row(pep, tag) for tag, pep in pilot.peptides[:5]])
    assert sha256(pilot.pep_a2) != before
    out = tmp_path / "metrics"
    res = _compute(pilot, out)
    assert res.returncode == 2 and _nothing_written(out), res.stdout + res.stderr
    assert "a2.tsv.gz" in res.stdout + res.stderr


def test_missing_input_exit_2(pilot: Pilot, tmp_path: Path) -> None:
    pilot.k12.unlink()
    out = tmp_path / "metrics"
    res = _compute(pilot, out)
    assert res.returncode == 2 and _nothing_written(out), res.stdout + res.stderr


def test_out_inside_the_repository_exit_2(pilot: Pilot) -> None:
    out = REPO / "local" / "m10_05_cli_guard_out"
    shutil.rmtree(out, ignore_errors=True)
    try:
        res = _compute(pilot, out)
        assert res.returncode == 2, res.stdout + res.stderr
        assert _nothing_written(out)
    finally:
        shutil.rmtree(out, ignore_errors=True)


def test_missing_compare_table_nulls_its_metric_only(pilot: Pilot, tmp_path: Path) -> None:
    (pilot.cmp / "tables" / "compare_stability.tsv").unlink()
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out))
    doc = _metrics(out)
    m5 = doc["m5"]
    reason = "missing:compare_stability.tsv"
    assert (m5 is None and any(reason in c for c in map(str, doc["warnings"]))) or (
        isinstance(m5, dict) and m5.get("reason") == reason
    ), m5
    assert doc["m4"]["spearman_rho"] == pytest.approx(0.83)  # the others are computed
    assert _cl(doc["checklist"][2]) == "fail"


# ------------------------------------------------------------------ criterion 8


def test_checklist_a1_a2_only(pilot: Pilot, tmp_path: Path) -> None:
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out))
    cl = _metrics(out)["checklist"]
    assert len(cl) == 5
    assert _cl(cl[0]).startswith("partial") and "A3 pending" in _cl(cl[0])
    assert _cl(cl[1]) == "pass"
    assert _cl(cl[2]) == "pass"
    assert _cl(cl[3]).startswith("pending") and "M10-12" in _cl(cl[3])
    assert _cl(cl[4]).startswith("pending") and "M10-06" in _cl(cl[4])


def test_checklist_with_a3_at_exit_0(pilot: Pilot, tmp_path: Path) -> None:
    a3 = write_a3(pilot.root, sha256(pilot.genome))
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out, "--a3", str(a3)))
    cl = _metrics(out)["checklist"]
    assert _cl(cl[0]) == "pass"
    assert _cl(cl[3]).startswith("pending") and _cl(cl[4]).startswith("pending")


def test_checklist_te1_fail(tmp_path: Path, world: World) -> None:
    p = build_pilot(tmp_path / "pilot", world, te1_misses=2)  # 18 / 20 = 0.9
    out = tmp_path / "metrics"
    _ok(_compute(p, out))
    doc = _metrics(out)
    assert doc["te1"]["passed"] is False
    assert _cl(doc["checklist"][1]) == "fail"


# ------------------------------------------------------------------ criterion 7: render


def _render(metrics: Path, out: Path, notes: Path | None = None) -> Any:
    args = ["render", "--metrics", str(metrics), "--out", str(out)]
    if notes is not None:
        args += ["--notes", str(notes)]
    return run_cli(*args, cwd=out.parent)


def _resolve(doc: Any, path: str) -> Any:
    cur = doc
    for part in path.split("."):
        if isinstance(cur, dict) and part in cur:
            cur = cur[part]
        elif isinstance(cur, list) and part.isdigit() and int(part) < len(cur):
            cur = cur[int(part)]
        else:
            raise KeyError(path)
    return cur


def _numeric(v: Any) -> bool:
    if isinstance(v, bool):
        return False
    if isinstance(v, int | float):
        return True
    return isinstance(v, dict) and isinstance(v.get("value"), int | float)


@pytest.fixture
def computed(pilot: Pilot, tmp_path: Path) -> tuple[Pilot, Path]:
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out))
    return pilot, out


def test_template_placeholders_resolve_and_numbers_cite_sources(
    computed: tuple[Pilot, Path],
) -> None:
    p, out = computed
    doc = _metrics(out)
    # M10-07: <!-- if/ifnot PATH --> blocks are applied first (metric 8 only with --entrapment)
    template = load_metrics().sections(TEMPLATE.read_text(encoding="utf-8"), doc)
    paths = PLACEHOLDER.findall(template)
    assert paths, "the template has no placeholder"
    for path in paths:
        _resolve(doc, path)  # every placeholder of the committed template resolves
    report = out / "REPORT.md"
    _ok(_render(out / "metrics.json", report))
    text = report.read_text(encoding="utf-8")
    assert "**missing:" not in text and "{{" not in text
    numeric = [path for path in paths if _numeric(_resolve(doc, path))]
    assert numeric, "the template renders no number"
    assert text.count("(metrics.json") >= len(numeric)
    t_lines, r_lines = template.splitlines(), text.splitlines()
    if len(t_lines) == len(r_lines):
        for tl, rl in zip(t_lines, r_lines, strict=True):
            n = sum(_numeric(_resolve(doc, x)) for x in PLACEHOLDER.findall(tl))
            assert rl.count("(metrics.json") >= n, rl
    # aggregates only: no peptide sequence of the fixture in the report
    assert not [pep for _, pep in p.peptides if pep in text]


def test_report_wording(computed: tuple[Pilot, Path]) -> None:
    _, out = computed
    report = out / "REPORT.md"
    _ok(_render(out / "metrics.json", report))
    text = report.read_text(encoding="utf-8")
    low = text.lower()
    assert "refuted" not in low
    assert "assumes independent genes" in low  # p_gene's label
    assert "nominal" in low  # the TE1 CI
    assert "no protein-level fdr control" in low  # metric 1 groups
    # example 10: islands 0.5; operon flags both empty; enriched categories missing in A1
    assert "none in either arm" in low
    assert "function module off in a1" in low
    for section in ("Status", "TE1", "Hypotheses", "Pending in pass 2", "Deviations"):
        assert section in text, section


def test_render_missing_key_exit_1_naming_the_path(computed: tuple[Pilot, Path]) -> None:
    _, out = computed
    template = TEMPLATE.read_text(encoding="utf-8")
    path = next(x for x in PLACEHOLDER.findall(template) if "." in x)
    doc = _metrics(out)
    parent = _resolve(doc, path.rsplit(".", 1)[0])
    del parent[path.rsplit(".", 1)[1]]
    bad = out / "broken.json"
    bad.write_text(json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    res = _render(bad, out / "REPORT_bad.md")
    assert res.returncode == 1, res.stdout + res.stderr
    assert path in res.stdout + res.stderr


@pytest.mark.parametrize(
    "notes",
    [
        "The misses are explained by QQQWWWRAK sequences.",  # a 7+ letter amino-acid string
        "Of 31415926 peptides most mapped.",  # a 4+ digit number not in metrics.json
    ],
)
def test_render_rejects_notes(computed: tuple[Pilot, Path], notes: str) -> None:
    _, out = computed
    nf = out / "notes.md"
    nf.write_text(notes + "\n", encoding="utf-8")
    res = _render(out / "metrics.json", out / "REPORT.md", nf)
    assert res.returncode == 1, res.stdout + res.stderr


def test_render_includes_valid_notes_verbatim(computed: tuple[Pilot, Path]) -> None:
    _, out = computed
    notes = "The TE1 result is explained by the annotation; nothing else is pending."
    nf = out / "notes.md"
    nf.write_text(notes + "\n", encoding="utf-8")
    report = out / "REPORT.md"
    _ok(_render(out / "metrics.json", report, nf))
    text = report.read_text(encoding="utf-8")
    assert notes in text
    assert text.index("Notes") < text.index(notes)


# ================================================================== review-round fixes
# Code review rev-1005 (B1-B3, S1-S10) and statistician review stat-1005 of 83ba1a99, with the
# lead's decisions of the review round.

SEP = " \u203a "  # metrics.SOURCE_SEP


def _edit_tsv(path: Path, col: str, fn: Any, *, drop: bool = False) -> None:
    """Rewrite a TSV: ``fn(row_index, value) -> value`` on column ``col``, or drop the column."""
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    j = header.index(col)
    out = []
    for i, line in enumerate(lines):
        cells = line.split("\t")
        if drop:
            del cells[j]
        elif i > 0:
            cells[j] = fn(i - 1, cells[j])
        out.append("\t".join(cells))
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def _edit_json(path: Path, fn: Any) -> None:
    doc = json.loads(path.read_text(encoding="utf-8"))
    fn(doc)
    path.write_text(json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def _precondition(res: Any, out: Path) -> None:
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert _nothing_written(out)
    assert "Traceback" not in res.stderr, res.stderr


# ---------------------------------------------------------------- B1: exit 2 for bad inputs


def test_b1_compare_table_with_a_bad_type_exit_2(pilot: Pilot, tmp_path: Path) -> None:
    _edit_tsv(pilot.cmp / "tables" / "compare_abundance.tsv", "n_shared", lambda i, v: "abc")
    out = tmp_path / "metrics"
    _precondition(_compute(pilot, out), out)


def test_b1_genes_tsv_missing_a_schema_column_exit_2(pilot: Pilot, tmp_path: Path) -> None:
    _edit_tsv(pilot.out_a2 / "tables" / "genes.tsv", "region_ids", None, drop=True)
    out = tmp_path / "metrics"
    res = _compute(pilot, out)
    _precondition(res, out)
    assert "genes.tsv" in res.stderr


def test_b1_malformed_convert_json_exit_2(pilot: Pilot, tmp_path: Path) -> None:
    pilot.conv_a2.write_text('{"counts": {', encoding="utf-8")
    out = tmp_path / "metrics"
    res = _compute(pilot, out)
    _precondition(res, out)
    assert pilot.conv_a2.name in res.stderr


def test_b1_summary_json_without_distinct_stripped_q01_exit_2(pilot: Pilot, tmp_path: Path) -> None:
    _edit_json(pilot.sum_a1, lambda d: d.pop("distinct_stripped_q01"))
    out = tmp_path / "metrics"
    _precondition(_compute(pilot, out), out)


def test_b1_invalid_run_config_exit_2(pilot: Pilot, tmp_path: Path) -> None:
    def bad(doc: dict[str, Any]) -> None:
        doc["config"]["input"]["q_value_max"] = "not-a-number"

    _edit_json(pilot.out_a2 / "run_manifest.json", bad)
    out = tmp_path / "metrics"
    _precondition(_compute(pilot, out), out)


def test_b1_cross_check_failure_stays_exit_1(tmp_path: Path, world: World) -> None:
    p = build_pilot(tmp_path / "pilot", world, conv_offset=1)
    res = _compute(p, tmp_path / "metrics")
    assert res.returncode == 1, res.stdout + res.stderr
    assert "Traceback" not in res.stderr


# ---------------------------------------------------------------- S7: no config


def test_s7_manifest_without_config_exit_2(pilot: Pilot, tmp_path: Path) -> None:
    _edit_json(pilot.out_a1 / "run_manifest.json", lambda d: d.pop("config"))
    out = tmp_path / "metrics"
    res = _compute(pilot, out)
    _precondition(res, out)
    assert "config" in res.stderr


# ---------------------------------------------------------------- S4: commits


def test_s4_giap_commit_from_the_manifest(computed: tuple[Pilot, Path]) -> None:
    _, out = computed
    doc = _metrics(out)
    assert doc["giap_commit"] == GIAP_COMMIT
    report = out / "REPORT.md"
    _ok(_render(out / "metrics.json", report))
    assert GIAP_COMMIT in report.read_text(encoding="utf-8")


def test_s4_git_failure_is_a_precondition(pilot: Pilot, tmp_path: Path) -> None:
    fake = tmp_path / "bin"
    fake.mkdir()
    git = fake / "git"
    git.write_text(
        "#!/bin/sh\necho 'fatal: not a git repository' >&2\nexit 128\n", encoding="utf-8"
    )
    git.chmod(0o755)
    out = tmp_path / "metrics"
    res = run_cli(*pilot.compute_args(out), cwd=pilot.root, path_prefix=fake)
    _precondition(res, out)


# ---------------------------------------------------------------- S5: render --out guard


def test_s5_render_out_inside_the_repository_exit_2(computed: tuple[Pilot, Path]) -> None:
    _, out = computed
    target = REPO / "local" / "m10_05_render_guard" / "REPORT_probe.md"
    shutil.rmtree(target.parent, ignore_errors=True)
    target.parent.mkdir(parents=True)
    try:
        res = _render(out / "metrics.json", target)
        assert res.returncode == 2, res.stdout + res.stderr
        assert not target.exists()
    finally:
        shutil.rmtree(target.parent, ignore_errors=True)


# ---------------------------------------------------------------- B2: render a missing table

MISSING_MARK = re.compile(r"\*\*missing: ([^*]+)\*\*")


def _render_ok(out: Path) -> str:
    report = out / "REPORT.md"
    _ok(_render(out / "metrics.json", report))
    text = report.read_text(encoding="utf-8")
    assert "{{" not in text
    return text


def test_b2_missing_compare_stability_renders_markers(pilot: Pilot, tmp_path: Path) -> None:
    (pilot.cmp / "tables" / "compare_stability.tsv").unlink()
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out))
    text = _render_ok(out)
    reasons = MISSING_MARK.findall(text)
    assert any("compare_stability.tsv" in r for r in reasons), reasons
    assert f"(metrics.json{SEP}m5" not in text  # every metric-5 placeholder is a marker
    assert f"(metrics.json{SEP}m4)" in text  # the others are rendered


def test_b2_missing_a2_peptide_hits_renders_markers(pilot: Pilot, tmp_path: Path) -> None:
    (pilot.out_a2 / "tables" / "peptide_hits.tsv").unlink()
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out))
    doc = _metrics(out)
    assert "peptide_hits.tsv" in str(doc["m2"]["A2"].get("reason"))
    assert "peptide_hits.tsv" in str(doc["te1"].get("reason"))
    text = _render_ok(out)
    reasons = MISSING_MARK.findall(text)
    assert any("peptide_hits.tsv" in r for r in reasons), reasons
    assert f"(metrics.json{SEP}m2.A2" not in text
    assert f"(metrics.json{SEP}m2.A1" in text and f"(metrics.json{SEP}h1)" in text
    assert _cl(doc["checklist"][1]) == "fail (no data)"


# ---------------------------------------------------------------- B3: full schema validation


def test_b3_enum_violation_in_a_compare_table_fails_bullet_3(pilot: Pilot, tmp_path: Path) -> None:
    _edit_tsv(pilot.cmp / "tables" / "compare_detection.tsv", "detection_status",
              lambda i, v: "bogus_status" if i == 0 else v)  # fmt: skip
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out))  # a warning and a fail, not an exit code
    doc = _metrics(out)
    assert _cl(doc["checklist"][2]) == "fail"
    assert any("compare_detection" in json.dumps(w) for w in doc["warnings"]), doc["warnings"]


def test_b3_no_private_schema_import_and_check_tsv_used() -> None:
    src = (REPO / "validation" / "pilot" / "metrics.py").read_text(encoding="utf-8")
    assert not re.search(r"import[^\n]*\b_schema\b", src), "private run_view._schema import"
    assert "check_tsv" in src


# ---------------------------------------------------------------- S1, should-fix 2: warnings


def test_s1_unknown_source_count_and_warning(tmp_path: Path, world: World) -> None:
    p = build_pilot(tmp_path / "pilot", world, unknown_source=2)
    out = tmp_path / "metrics"
    _ok(_compute(p, out))
    doc = _metrics(out)
    assert doc["te1"]["unknown_source"] == 2
    assert "unknown_source" in _codes(doc)


def test_s1_no_unknown_source_no_warning(computed: tuple[Pilot, Path]) -> None:
    doc = _metrics(computed[1])
    assert doc["te1"]["unknown_source"] == 0
    assert "unknown_source" not in _codes(doc) and "target_q_na" not in _codes(doc)


def test_target_q_na_warns(tmp_path: Path, world: World) -> None:
    p = build_pilot(tmp_path / "pilot", world, na_q_a2=1)
    out = tmp_path / "metrics"
    _ok(_compute(p, out))
    doc = _metrics(out)
    assert doc["te1"]["target_q_na"] == 1 and doc["te1"]["n_peptides"] == 21
    assert "target_q_na" in _codes(doc)


# ---------------------------------------------------------------- S3: ATTEMPTS.tsv and sacct

ATTEMPT_COLS = ("arm", "attempt", "outdir", "job", "resume", "git_commit", "state", "exit_code")


def _attempts(runs: Path, rows: list[tuple[str, ...]]) -> Path:
    path = runs / "ATTEMPTS.tsv"
    path.write_text("\n".join("\t".join(r) for r in [ATTEMPT_COLS, *rows]) + "\n", encoding="utf-8")
    return path


def _attempt_rows(runs: Path) -> list[tuple[str, ...]]:
    c = GIAP_COMMIT[:7]
    return [
        ("A2", "1", str(runs / "out_a2"), "1001", "fresh", c, "superseded", "1"),
        ("A1", "1", str(runs / "out_a1"), "2001", "fresh", c, "accepted", "0"),
        ("A2", "2", str(runs / "out_a2.r2"), "1002", "fresh", c, "accepted", "4"),
        ("A2", "2", str(runs / "out_a2.r2"), "1003", "resume-infra", c, "accepted", "0"),
        ("A2", "3", str(runs / "out_a2.r3"), "1004", "fresh", c, "superseded", "1"),
    ]


def test_s3_resolve_attempt_last_non_superseded(tmp_path: Path) -> None:
    from helpers.pilot_metrics_m10 import load_metrics

    m = load_metrics()
    runs = tmp_path / "runs"
    runs.mkdir()
    assert m.resolve_attempt(runs / "out_a2") == (runs / "out_a2", None)  # no ATTEMPTS.tsv
    _attempts(runs, _attempt_rows(runs))
    assert m.resolve_attempt(runs / "out_a2") == (runs / "out_a2.r2", "1003")
    assert m.resolve_attempt(runs / "out_a1") == (runs / "out_a1", "2001")
    _attempts(
        runs,
        [r for r in _attempt_rows(runs) if r[0] == "A1"]
        + [("A2", "1", str(runs / "out_a2"), "1001", "fresh", "x", "superseded", "1")],
    )
    assert m.resolve_attempt(runs / "out_a2")[1] is None  # only superseded attempts


def _sacct(path: Path) -> Path:
    cols = "JobID|JobName|State|Elapsed|MaxRSS|NodeList|AllocCPUS"
    rows = ["1002|giap-run-a2|FAILED|00:10:00||sh04-01n02|114",
            "1003|giap-run-a2|COMPLETED|01:42:00||sh04-01n01|114",
            "1003.batch|batch|COMPLETED|01:42:00|3000000K|sh04-01n01|2",
            "1003.extern|extern|COMPLETED|01:42:01|0|sh04-01n01|114",
            "2001|giap-run-a1|COMPLETED|02:10:00||sh04-01n03|114",
            "2001.batch|batch|COMPLETED|02:10:00|2000000K|sh04-01n03|2"]  # fmt: skip
    path.write_text("\n".join([cols, *rows]) + "\n", encoding="utf-8")
    return path


def test_s3_compute_follows_attempts_and_picks_the_head_sacct_row(
    pilot: Pilot, tmp_path: Path
) -> None:
    runs = pilot.out_a2.parent
    final = runs / "out_a2.r2"
    shutil.copytree(pilot.out_a2, final)
    _edit_json(pilot.out_a2 / "run_manifest.json", lambda d: d.update(exit_code=1))  # superseded
    _attempts(runs, _attempt_rows(runs))
    sacct = _sacct(tmp_path / "sacct.tsv")
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out, "--sacct", str(sacct)))  # out_a2 itself would be exit 2
    doc = _metrics(out)
    assert doc["m7"]["A2"]["sacct"]["JobID"] == "1003"
    assert doc["m7"]["A2"]["sacct"]["Elapsed"] == "01:42:00"
    assert doc["m7"]["A1"]["sacct"]["JobID"] == "2001"
    assert doc["inputs"]["manifest_a2"]["path"] == str((final / "run_manifest.json").resolve())
    assert not [c for c in _codes(doc) if c.startswith("sacct")]


def test_s3_sacct_without_attempts_warns(pilot: Pilot, tmp_path: Path) -> None:
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out, "--sacct", str(_sacct(tmp_path / "sacct.tsv"))))
    doc = _metrics(out)
    assert "sacct_ignored_no_attempts" in _codes(doc)
    assert doc["m7"]["A2"]["sacct"] is None


# ---------------------------------------------------------------- S10: CLI-level mutants


def test_s10_locus_tags_differ_from_bakta_exit_2(pilot: Pilot, tmp_path: Path) -> None:  # 70
    genes = pilot.out_a2 / "tables" / "genes.tsv"
    lines = genes.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    extra = lines[1].split("\t")
    extra[header.index("locus_tag")] = "EXTRA_99999"
    genes.write_text("\n".join([*lines, "\t".join(extra)]) + "\n", encoding="utf-8")
    out = tmp_path / "metrics"
    _precondition(_compute(pilot, out), out)


def test_s10_negative_raw_vs_written_warns(tmp_path: Path, world: World) -> None:  # 72
    p = build_pilot(tmp_path / "pilot", world, q_crossings=4, raw_vs_written=-1)
    out = tmp_path / "metrics"
    _ok(_compute(p, out))
    assert "raw_vs_written_q01_out_of_range" in _codes(_metrics(out))


def test_s10_q_crossings_reach_te1(tmp_path: Path, world: World) -> None:  # 80
    p = build_pilot(tmp_path / "pilot", world, q_crossings=2, raw_vs_written=1)
    out = tmp_path / "metrics"
    _ok(_compute(p, out))
    sq = _metrics(out)["te1"]["share_without_q_at_cutoff"]
    assert sq["q_crossings"] == 2 and sq.get("reason") != "no_crossings"
    assert sq["removed"] == 0 and sq["value"] == pytest.approx(1.0)


def test_s10_checklist_a3_at_exit_1_fails(pilot: Pilot, tmp_path: Path) -> None:  # 56
    a3 = write_a3(pilot.root, sha256(pilot.genome), exit_code=1)
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out, "--a3", str(a3)))
    assert _cl(_metrics(out)["checklist"][0]) == "fail"


# ---------------------------------------------------------------- S8, statistician B1/B2/SF3/SF4


def _template() -> str:
    return TEMPLATE.read_text(encoding="utf-8")


def _paths() -> list[str]:
    return PLACEHOLDER.findall(_template())


@pytest.mark.parametrize(
    "path",
    [
        "giap_commit", "spec",  # S4, section 1
        "te1.misses.exact_other_class", "te1.misses.unplaced",  # S8, SF3
        "te1.non_atg_nterm.peptides",  # SF3
        "te1.share_without_q_at_cutoff.removed",  # statistician B2
        "te1.share_without_q_at_cutoff.q_crossings",
        "te1.share_without_q_at_cutoff.pass_becomes_fail",
        "h1.strain_specific_share", "h1.b", "h1.c",  # S8, SF4
    ],
)  # fmt: skip
def test_s8_template_has_the_placeholder(path: str) -> None:
    assert path in _paths(), path


def _flat(text: str) -> str:
    return " ".join(text.split())


def test_stat_b1_wilson_sentence_after_the_h1_intervals() -> None:
    t = _flat(_template())
    sentence = (
        "These are nominal 95 % Wilson intervals; they describe the effects and are not the "
        "test (a two-sided 95 % interval and the one-sided 0.05 test can disagree near the "
        "boundary)."
    )
    assert sentence in t
    assert t.index("{{h1.prop_a2_only}}") < t.index(sentence)


def test_stat_sf4_not_evaluable_stated_apart_from_not_supported() -> None:
    t = _flat(_template()).lower()
    assert "`not_evaluable` means the rule could not be applied" in t
    assert "`not_supported` means not supported by the pre-registered rule" in t
    assert 'a verdict other than "supported" means not supported' not in t


def test_s8_rendered_breakdowns_are_k_v_pairs(tmp_path: Path, world: World) -> None:
    p = build_pilot(tmp_path / "pilot", world, te1_misses=2)
    out = tmp_path / "metrics"
    _ok(_compute(p, out))
    text = _render_ok(out)
    assert "no_hit: 2" in text  # te1.misses.unplaced
    assert "{'" not in text and "{}" not in text


def test_s8_every_metric_3_row_is_rendered_for_both_arms(computed: tuple[Pilot, Path]) -> None:
    _, out = computed
    doc = _metrics(out)
    text = _render_ok(out)
    for arm in ("A1", "A2"):
        names = doc["m3"][arm]["order"]
        assert {"unknown", "no_peptides", "pseudogene"} <= set(names)
        for name in names:
            assert f"(metrics.json{SEP}m3.{arm}.subsets.{name})" in text, (arm, name)


# ---------------------------------------------------------------- lead decision 6: H3 numbers


FAMILY_KEYS = ("n_a", "n_b", "n_shared", "jaccard")


def test_decision6_template_shows_rho_ci_and_every_family_count() -> None:
    paths = set(_paths())
    rho = {"m4.spearman_rho", "m4.ci_low", "m4.ci_high"}
    alt = {"h3.criteria.rho.rho", "h3.criteria.rho.ci_low", "h3.criteria.rho.ci_high"}
    assert rho <= paths or alt <= paths
    for f in ("islands", "operon_flags", "enriched_categories"):
        assert f"h3.criteria.{f}.union" in paths, f
        for k in FAMILY_KEYS:
            assert {f"m5.{f}.{k}", f"h3.criteria.{f}.{k}"} & paths, (f, k)


def test_decision6_render_shows_h3_numbers_when_not_evaluable(tmp_path: Path, world: World) -> None:
    p = build_pilot(tmp_path / "pilot", world, abundance_over={"flags": "too_few_shared"})
    out = tmp_path / "metrics"
    _ok(_compute(p, out))
    doc = _metrics(out)
    assert doc["h3"]["verdict"] == "not_evaluable"
    text = _render_ok(out)
    assert "**missing:" not in text
    assert "0.83" in text and "0.76" in text and "0.88" in text  # rho [ci_low, ci_high]
    for f in ("islands", "operon_flags", "enriched_categories"):
        assert f"(metrics.json{SEP}h3.criteria.{f})" in text, f
