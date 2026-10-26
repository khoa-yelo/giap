"""``metrics.py compute --a3 …`` and ``render`` for metric 6 and hypothesis 4 (M10-10b;
docs/specs/M10-10-a3-run.md, "``metrics.py`` changes (10b)", acceptance criterion 10, example 13,
the lead decision after the M10-15a review).

Inputs: M10-05's synthetic pilot (``helpers.pilot_metrics_m10.build_pilot``) plus the A3 inputs
of ``helpers.a3_metrics_m10_10b.build_a3`` (synthetic de novo peptides, examples 9 and 10 as the
A2/A3 tables). Nothing comes from ``$GIAP_ROOT``.

Interface fixed here (test-author choices, reversible):

* the new options are ``--peptides-a3``, ``--convert-a3``, ``--cmp3`` (spec) and the optional
  ``--strain-db`` (lead decision S4); ``CMP3`` is a ``giap compare`` directory
  (``CMP3/tables/<stem>.tsv``);
* the A3 manifest records A2's table under the input role ``reference_peptides``;
  ``CONV_A3`` is 10a's ``a3.convert.json`` (``failures``, ``outputs.peptides.sha256``);
* without ``--a3`` the top-level keys are M10-05's (with M10-15's ``evidence``/``m8``/``te_*``)
  and nothing else; ``m6.overlap_a2_a3`` and ``h4`` carry ``reason = "a3_not_run"``;
* the fixture has no placement log line, so ``h4`` warns ``denovo_count_unchecked`` (S6);
  warning codes of ``m6`` and ``h4`` are also listed in the top-level ``warnings``;
* ``render`` writes one H4 table row per bin, citing ``h4.bins.<i>``, and every caveat string
  of ``h4.caveats``; the per-arm metric-6 tables cite ``m6.<arm>.novel_orfs.by_class.<class>``
  and ``m6.<arm>.variants.alt.<type>``; the overlap table cites both sets of each family and
  ``m6.overlap_a2_a3.novel_orfs.orf_id_jaccard``.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pytest
from helpers.a3_metrics_m10_10b import (
    A3,
    ALT_TYPES,
    NOVEL_CLASSES,
    VERDICTS,
    build_a3,
    example_stability,
)
from helpers.pilot_metrics_m10 import (
    TEMPLATE,
    Pilot,
    World,
    build_pilot,
    output_files,
    run_cli,
    sha256,
    sim_world,
)

SEP = " \u203a "
NEEDS = "--a3 needs --peptides-a3, --convert-a3 and --cmp3"
#: M10-05's top-level keys as main writes them before M10-10b (M10-06/07/15 hooks included)
M10_05_KEYS = {
    "schema",
    "spec",
    "warnings",
    "giap_commit",
    "m1",
    "m2",
    "m3",
    "m4",
    "m5",
    "m7",
    "te1",
    "te1_a1_diagnostic",
    "h1",
    "h3",
    "te_real",
    "te7",
    "m6",
    "h4",
    "m8",
    "evidence",
    "checklist",
    "inputs",
}
PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> World:
    return sim_world(tmp_path_factory.mktemp("m10_10b_sim"))


@pytest.fixture
def pilot(tmp_path: Path, world: World) -> Pilot:
    return build_pilot(tmp_path / "pilot", world)


def _compute(p: Pilot, out: Path, *extra: str) -> Any:
    return run_cli(*p.compute_args(out, *extra), cwd=p.root)


def _ok(res: Any) -> None:
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


def _metrics(out: Path) -> dict[str, Any]:
    doc: dict[str, Any] = json.loads((out / "metrics.json").read_text(encoding="utf-8"))
    return doc


def _codes(doc: dict[str, Any]) -> list[str]:
    return [w if isinstance(w, str) else str(w.get("code")) for w in doc.get("warnings", [])]


def _nothing_written(out: Path) -> bool:
    return not out.exists() or not any(out.iterdir())


def _with_a3(p: Pilot, world: World, tmp: Path, *extra: str, **kw: Any) -> tuple[Any, Path, A3]:
    a3 = build_a3(p, world, **kw)
    out = tmp / "metrics"
    return _compute(p, out, *a3.args(), *extra), out, a3


# ------------------------------------------------------------------ without --a3


def test_without_a3_keys_and_a3_not_run_stubs(pilot: Pilot, tmp_path: Path) -> None:
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out))
    doc = _metrics(out)
    assert set(doc) == M10_05_KEYS  # no new top-level key
    assert doc["m6"]["overlap_a2_a3"]["reason"] == "a3_not_run"
    assert doc["h4"]["reason"] == "a3_not_run"
    # M10-05's fixture has no novel_orfs.tsv / variants.tsv: per-arm, per-table reasons
    for arm in ("A1", "A2"):
        assert doc["m6"][arm]["novel_orfs"]["reason"] == "missing:novel_orfs.tsv"
        assert doc["m6"][arm]["variants"]["reason"] == "missing:variants.tsv"
    assert not {"peptides_a3", "convert_a3", "manifest_a3"} & set(doc["inputs"])
    assert "pending (M10-10)" not in (out / "metrics.json").read_text(encoding="utf-8")
    first = output_files(out)
    _ok(_compute(pilot, out))
    assert output_files(out) == first


def test_without_a3_per_arm_summaries_computed(pilot: Pilot, world: World, tmp_path: Path) -> None:
    build_a3(pilot, world)  # writes A1/A2's novel_orfs.tsv and variants.tsv; --a3 not passed
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out))
    m6 = _metrics(out)["m6"]
    assert m6["A2"]["novel_orfs"]["n_rows"] == 2 and m6["A1"]["novel_orfs"]["n_rows"] == 1
    assert m6["A2"]["variants"]["alt"]["substitution"]["high"] == 1
    assert m6["overlap_a2_a3"]["reason"] == "a3_not_run"


# ------------------------------------------------------------------ example 13 and preconditions


@pytest.mark.parametrize("skip", ["--peptides-a3", "--convert-a3", "--cmp3"])
def test_example_13_a3_without_a_companion_exit_2(
    pilot: Pilot, world: World, tmp_path: Path, skip: str
) -> None:
    a3 = build_a3(pilot, world)
    out = tmp_path / "metrics"
    res = _compute(pilot, out, *a3.args(skip=[skip]))
    assert res.returncode == 2, res.stdout + res.stderr
    assert NEEDS in res.stderr and _nothing_written(out), res.stderr


@pytest.mark.parametrize("keep", ["--peptides-a3", "--convert-a3", "--cmp3"])
def test_companion_without_a3_exit_2(pilot: Pilot, world: World, tmp_path: Path, keep: str) -> None:
    a3 = build_a3(pilot, world)
    out = tmp_path / "metrics"
    skip = [o for o in ("--a3", "--peptides-a3", "--convert-a3", "--cmp3") if o != keep]
    res = _compute(pilot, out, *a3.args(skip=skip))
    assert res.returncode == 2 and _nothing_written(out), res.stdout + res.stderr
    assert res.stderr.startswith("precondition:"), res.stderr  # not argparse's usage error


PRECONDITIONS = {
    "a3_exit_code": {"exit_code": 1},
    "genome_sha": {"genome_sha": "0" * 64},
    "peptides_sha": {"peptides_sha": "1" * 64},
    "reference_not_a2": {"reference_sha": "2" * 64},
    "convert_failures": {"conv_failures": ["rule 1: 3 table rows, 4 .ssl data rows"]},
    "convert_peptides_sha": {"conv_peptides_sha": "3" * 64},
    "cmp3_missing_table": {"drop_cmp_table": "compare_stability"},
}


@pytest.mark.parametrize("name", sorted(PRECONDITIONS))
def test_a3_preconditions_exit_2(pilot: Pilot, world: World, tmp_path: Path, name: str) -> None:
    res, out, _ = _with_a3(pilot, world, tmp_path, **PRECONDITIONS[name])
    assert res.returncode == 2, f"{name}: exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert _nothing_written(out), name
    assert res.stderr.startswith("precondition:"), res.stderr  # not argparse's usage error


# ------------------------------------------------------------------ with --a3


def test_with_a3_writes_m6_h4_inputs_and_checklist(
    pilot: Pilot, world: World, tmp_path: Path
) -> None:
    res, out, a3 = _with_a3(pilot, world, tmp_path)
    _ok(res)
    doc = _metrics(out)
    assert set(doc) == M10_05_KEYS
    m6 = doc["m6"]
    for arm in ("A1", "A2", "A3"):
        assert set(m6[arm]["novel_orfs"]["by_class"]) == set(NOVEL_CLASSES), arm
        assert set(m6[arm]["variants"]["alt"]) == set(ALT_TYPES), arm
    ov = m6["overlap_a2_a3"]
    assert "reason" not in ov
    assert ov["novel_orfs"]["called"]["jaccard"] == pytest.approx(0.25)  # example 9
    assert ov["novel_orfs"]["orf_id_jaccard"] == pytest.approx(0.0)
    assert ov["alt_sites"]["alt"]["jaccard"] == pytest.approx(0.5)  # example 10
    codes = _codes(doc)
    assert "novel_orf_key_differs" in codes and "denovo_count_unchecked" in codes
    h4 = doc["h4"]
    assert h4["verdict"] in VERDICTS and h4["verdict"] == "not_evaluable"  # n_b = 3 < 20
    assert sum(b["n_b"] for b in h4["bins"]) == len(a3.denovo) == 30
    assert sum(b["k_b"] for b in h4["bins"]) == a3.n_confirmed  # I/L folded on both sides
    assert [b["score_bin"] for b in h4["bins"]] == list(range(10))
    for b in h4["bins"]:
        s = b["L_b_searchable"]
        assert s == "not_evaluable" or s.get("reason") == "not_evaluable"  # no --strain-db
    inputs = doc["inputs"]
    assert inputs["peptides_a3"]["sha256"] == sha256(a3.pep_a3)
    assert inputs["convert_a3"]["sha256"] == sha256(a3.conv_a3)
    assert inputs["manifest_a3"]["sha256"] == sha256(a3.out_a3 / "run_manifest.json")
    status = doc["checklist"][0]
    assert (status["status"] if isinstance(status, dict) else status) == "pass"
    text = (out / "metrics.json").read_text(encoding="utf-8")
    assert not [p for p in a3.denovo if p in text]  # no sequence in metrics.json
    first = output_files(out)
    _ok(_compute(pilot, out, *a3.args()))
    assert output_files(out) == first


def test_strain_db_fills_l_b_searchable(pilot: Pilot, world: World, tmp_path: Path) -> None:
    res, out, a3 = _with_a3(pilot, world, tmp_path)
    _ok(res)
    _ok(_compute(pilot, out, *a3.args(), "--strain-db", str(a3.strain_db)))
    bins = _metrics(out)["h4"]["bins"]
    sh = [b["L_b_searchable"] for b in bins]
    assert all(isinstance(s, dict) and "k" in s for s in sh), sh
    assert sum(s["k"] for s in sh) == a3.n_confirmed
    assert a3.n_confirmed <= sum(s["n"] for s in sh) <= len(a3.denovo)


def test_doctored_cmp3_alt_row_exit_1(pilot: Pilot, world: World, tmp_path: Path) -> None:
    """Example 10: compare says 1.0 where metric 6 counts 0.5 → ``Inconsistent``, exit 1."""
    res, out, _ = _with_a3(
        pilot, world, tmp_path, stability=example_stability(alt=[1, 2, 1, 1.0, None])
    )
    assert res.returncode == 1, res.stdout + res.stderr
    assert _nothing_written(out)


def test_doctored_mapping_fdr_bins_exit_1(pilot: Pilot, world: World, tmp_path: Path) -> None:
    """Definition 2: a ``score_bin`` label GIAP's bins never give → ``Inconsistent``, exit 1."""
    res, out, _ = _with_a3(pilot, world, tmp_path, fdr_extra_bin=True)
    assert res.returncode == 1, res.stdout + res.stderr
    assert _nothing_written(out)


# ------------------------------------------------------------------ render


def _render(metrics: Path, dest: Path, cwd: Path) -> Any:
    return run_cli("render", "--metrics", str(metrics), "--out", str(dest), cwd=cwd)


def test_template_sections_10b() -> None:
    t = TEMPLATE.read_text(encoding="utf-8")
    assert "pending M10-10" not in t
    assert "{{m6.overlap_a2_a3.novel_orfs.orf_id_jaccard}}" in t
    assert "{{h4.verdict}}" in t
    assert not PEPTIDE_LIKE.search(re.sub(r"\{\{[^}]*\}\}", "", t))


def test_render_resolves_metric6_and_h4(pilot: Pilot, world: World, tmp_path: Path) -> None:
    res, out, _ = _with_a3(pilot, world, tmp_path)
    _ok(res)
    dest = tmp_path / "REPORT.md"
    r = _render(out / "metrics.json", dest, pilot.root)
    _ok(r)
    text = dest.read_text(encoding="utf-8")
    doc = _metrics(out)
    for arm in ("A1", "A2", "A3"):
        for cls in NOVEL_CLASSES:
            assert f"metrics.json{SEP}m6.{arm}.novel_orfs.by_class.{cls})" in text, (arm, cls)
        for t in ALT_TYPES:
            assert f"metrics.json{SEP}m6.{arm}.variants.alt.{t})" in text, (arm, t)
    for fam, sets in (("novel_orfs", ("called", "any_row")), ("alt_sites", ("alt", "high"))):
        for s in sets:
            assert f"metrics.json{SEP}m6.overlap_a2_a3.{fam}.{s})" in text, (fam, s)
    for i in range(len(doc["h4"]["bins"])):
        assert f"metrics.json{SEP}h4.bins.{i}" in text, i
    for cav in doc["h4"]["caveats"]:
        assert cav in text, cav
    m6_section = text.split("### Metric 6", 1)[1].split("## 4.", 1)[0]
    h4_section = text.split("**H4**", 1)[1].split("## 7.", 1)[0]
    assert "descriptive, unestimated" in m6_section.lower()
    assert "**missing:" not in m6_section and "**missing:" not in h4_section
    assert not PEPTIDE_LIKE.search(text)


def test_render_without_a3_exit_0(pilot: Pilot, tmp_path: Path) -> None:
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out))
    dest = tmp_path / "REPORT.md"
    _ok(_render(out / "metrics.json", dest, pilot.root))
    assert "a3_not_run" in dest.read_text(encoding="utf-8")
