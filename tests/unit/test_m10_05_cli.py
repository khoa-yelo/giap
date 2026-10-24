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
  4 = runtimes, 5 = TE2-TE5); each entry a status string or a dict with ``status``.
"""

from __future__ import annotations

import json
import re
import shutil
from pathlib import Path
from typing import Any

import pytest
from helpers.pilot_metrics_m10 import (
    REPO,
    TEMPLATE,
    Pilot,
    World,
    build_pilot,
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
    for key in ("inputs", "spec", "m1", "m2", "m3", "m4", "m5", "m7", "te1",
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
    (pilot.cmp / "compare_stability.tsv").unlink()
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
    template = TEMPLATE.read_text(encoding="utf-8")
    paths = PLACEHOLDER.findall(template)
    assert paths, "the template has no placeholder"
    doc = _metrics(out)
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
