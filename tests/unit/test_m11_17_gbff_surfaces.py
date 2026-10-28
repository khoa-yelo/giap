"""M11-17 ``inputs.bakta_gbff`` surfaces: AC 2 (``run_preflight`` wiring and placement), AC 3
(never changes results; ``--strict`` exit 3), AC 4 (schema), AC 5 (``run_manifest.json``
``warnings``) and AC 6 (``gippa preflight`` CLI) (docs/specs/M11-17-preflight-gbff-warning.md
"Where it runs and where it surfaces"; GIPPA_SPEC.md §5.1, §7.0, §9.2; M1-12).
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Any

import jsonschema
import pandas as pd
import pytest
from click.testing import CliRunner
from helpers.manifest_data import build, preflight_doc

from gippa.config import GippaConfig, resolve_config
from gippa.exitcodes import ExitCode
from gippa.preflight.checks import Check

RULE_ID = "inputs.bakta_gbff"
ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = ROOT / "src" / "gippa" / "schemas" / "preflight.schema.json"
LOW = [
    "preflight.min_n50_warn=1",
    "preflight.min_largest_contig_warn=1",
    "preflight.min_total_bp_warn=1",
]
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"
GENOME = ">c1\n" + ("ACGTTGCAAGGCTTAACCGT" * 300) + "\n"


@pytest.fixture(autouse=True)
def _isolated_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    db = tmp_path / "db_empty"
    db.mkdir()
    repo = tmp_path / "repo_empty"
    repo.mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    monkeypatch.setenv("GIPPA_REPO", str(repo))
    monkeypatch.setenv("GIPPA_ROOT", str(repo))
    monkeypatch.delenv("NXF_APPTAINER_CACHEDIR", raising=False)


def _cfg(*extra: str) -> GippaConfig:
    return resolve_config(overrides=[*LOW, *extra])


def _probe() -> Any:
    from gippa.preflight.environment import GIB, EnvProbe, requirements

    found: dict[str, str | None] = {
        r.name: r.min_version or "1.0"
        for r in requirements(["core", "mobile", "bgc", "function"], "genome")
    }
    return EnvProbe(found=found, free_disk_bytes=500 * GIB, total_ram_bytes=64 * GIB)


def _files(tmp_path: Path) -> tuple[Path, Path, Path]:
    g = tmp_path / "g.fa"
    g.write_text(GENOME)
    p = tmp_path / "proteins.faa"
    p.write_text(f">P1\n{PROTEIN}\n")
    pep = tmp_path / "peptides.tsv"
    pep.write_text("peptide\tprotein_ids\nPEPTIDEK\tP1\n")
    return g, p, pep


def _bakta(tmp_path: Path, name: str = "sample.json", gbff: str | None = None) -> Path:
    d = tmp_path / "bakta"
    d.mkdir(exist_ok=True)
    js = d / name
    js.write_text("{}")
    if gbff is not None:
        (d / gbff).write_bytes(b"LOCUS\n")
    return js


def _run(tmp_path: Path, annotation: Path | None, *extra: str) -> Any:
    from gippa.preflight.run import PreflightInputs, run_preflight

    g, p, pep = _files(tmp_path)
    inp = PreflightInputs(genome=g, proteins=p, peptides=pep, annotation=annotation)
    return run_preflight(inp, _cfg(*extra), _probe(), "preflight")


def _found(report: Any) -> list[Check]:
    return [c for c in report.checks if c.rule_id == RULE_ID]


# ---------------------------------------------------------------- AC 2: wiring


@pytest.mark.parametrize(("gbff", "want"), [(None, "warn"), ("sample.gbff", "ok"),
                                            ("sample.gbff.gz", "ok")])  # fmt: skip
def test_ac2_annotation_json_gives_one_check(tmp_path: Path, gbff: str | None, want: str) -> None:
    js = _bakta(tmp_path, gbff=gbff)
    report, _ = _run(tmp_path, js)
    found = _found(report)
    assert len(found) == 1, [c.rule_id for c in report.checks]
    assert found[0].level == want


def test_ac2_configured_json_gives_one_check(tmp_path: Path) -> None:
    js = _bakta(tmp_path)
    report, _ = _run(tmp_path, None, f"annotation.precomputed.bakta_json={js}")
    found = _found(report)
    assert len(found) == 1
    assert found[0].level == "warn"
    assert str(tmp_path / "bakta" / "sample.gbff") in found[0].message


def test_ac2_bgc_off_is_info(tmp_path: Path) -> None:
    js = _bakta(tmp_path)
    report, _ = _run(tmp_path, js, "modules=[core, function]")
    found = _found(report)
    assert len(found) == 1 and found[0].level == "info"


def test_ac2_equals_the_rule_function(tmp_path: Path) -> None:
    from gippa.preflight.gbff import bgc_gbff_check

    js = _bakta(tmp_path)
    report, _ = _run(tmp_path, js)
    assert _found(report) == [bgc_gbff_check(_cfg(), js)]


@pytest.mark.parametrize("name", [None, "a.gff3", "a.gbff"])
def test_ac2_no_bakta_json_no_check(tmp_path: Path, name: str | None) -> None:
    ann = None
    if name is not None:
        ann = tmp_path / name
        ann.write_text("")
    report, _ = _run(tmp_path, ann)
    assert _found(report) == []


def test_ac2_placed_after_inputs_before_translation_table(tmp_path: Path) -> None:
    js = _bakta(tmp_path)
    report, _ = _run(tmp_path, js)
    ids = [c.rule_id for c in report.checks]
    i = ids.index(RULE_ID)
    assert ids[i - 1] == "inputs.annotation_kind", ids
    assert ids[i + 1] == "config.translation_table", ids
    assert all(r.startswith("inputs.") for r in ids[:i]), ids


# ---------------------------------------------------------------- AC 3: never changes results


def _norm_equal(a: Any, b: Any) -> None:
    for f in dataclasses.fields(a):
        x, y = getattr(a, f.name), getattr(b, f.name)
        if isinstance(x, pd.DataFrame) or isinstance(y, pd.DataFrame):
            assert (x is None) == (y is None), f.name
            if x is not None:
                pd.testing.assert_frame_equal(x, y)
        else:
            assert x == y, f.name


def test_ac3_warning_changes_nothing_but_checks(tmp_path: Path) -> None:
    from gippa.preflight.run import exit_code

    warn_dir, ok_dir = tmp_path / "w", tmp_path / "o"
    warn_dir.mkdir()
    ok_dir.mkdir()
    rep_w, norm_w = _run(warn_dir, _bakta(warn_dir))
    rep_o, norm_o = _run(ok_dir, _bakta(ok_dir, gbff="sample.gbff"))
    assert [c.level for c in _found(rep_w)] == ["warn"]
    assert [c.level for c in _found(rep_o)] == ["ok"]
    assert exit_code(rep_w, False) == exit_code(rep_o, False) == ExitCode.SUCCESS
    _norm_equal(norm_w, norm_o)
    others_w = [(c.rule_id, c.level) for c in rep_w.checks if c.rule_id != RULE_ID]
    others_o = [(c.rule_id, c.level) for c in rep_o.checks if c.rule_id != RULE_ID]
    assert others_w == others_o


def test_ac3_strict_exits_3(tmp_path: Path) -> None:
    from gippa.preflight.run import exit_code

    report, _ = _run(tmp_path, _bakta(tmp_path))
    assert [c.level for c in _found(report)] == ["warn"]
    others = [c for c in report.checks if c.rule_id != RULE_ID and c.level in {"warn", "error"}]
    assert others == [], [(c.rule_id, c.level) for c in others]
    assert exit_code(report, True) == ExitCode.WARNINGS_STRICT


# ---------------------------------------------------------------- AC 4: schema


@pytest.mark.parametrize("gbff", [None, "sample.gbff"])
def test_ac4_preflight_json_validates(tmp_path: Path, gbff: str | None) -> None:
    from gippa.preflight.run import write_outputs

    out = tmp_path / "out"
    report, norm = _run(tmp_path, _bakta(tmp_path, gbff=gbff))
    write_outputs(out, report, norm)
    doc = json.loads((out / "preflight" / "preflight.json").read_text())
    jsonschema.validate(doc, json.loads(SCHEMA_PATH.read_text()))
    rows = [c for c in doc["checks"] if c["rule_id"] == RULE_ID]
    assert len(rows) == 1
    md = (out / "preflight" / "preflight.md").read_text()
    assert RULE_ID in md


# ---------------------------------------------------------------- AC 5: manifest


def test_ac5_manifest_lists_the_warning(tmp_path: Path) -> None:
    from gippa.preflight.gbff import bgc_gbff_check

    c = bgc_gbff_check(_cfg(), _bakta(tmp_path))
    assert c is not None and c.level == "warn"
    m = build(preflight=preflight_doc(checks=[dataclasses.asdict(c)]))
    assert f"preflight:{RULE_ID}" in m["warnings"]


@pytest.mark.parametrize(("extra", "gbff"), [((), "sample.gbff"), (("modules=[core]",), None)])
def test_ac5_no_manifest_warning_at_ok_or_info(
    tmp_path: Path, extra: tuple[str, ...], gbff: str | None
) -> None:
    from gippa.preflight.gbff import bgc_gbff_check

    c = bgc_gbff_check(_cfg(*extra), _bakta(tmp_path, gbff=gbff))
    assert c is not None and c.level in {"ok", "info"}
    m = build(preflight=preflight_doc(checks=[dataclasses.asdict(c)]))
    assert f"preflight:{RULE_ID}" not in m["warnings"]


# ---------------------------------------------------------------- AC 6: CLI


def _cli(tmp_path: Path, js: Path, *flags: str) -> Any:
    from gippa.cli import cli

    g, p, pep = _files(tmp_path)
    out = tmp_path / "out"
    args = ["preflight", "--genome", str(g), "--proteins", str(p), "--peptides", str(pep),
            "--annotation", str(js), "--outdir", str(out), "--profile", "none"]  # fmt: skip
    for s in LOW:
        args += ["--set", s]
    return CliRunner().invoke(cli, [*args, *flags], catch_exceptions=False), out


def test_ac6_cli_warns_and_exits_0(tmp_path: Path) -> None:
    js = _bakta(tmp_path)
    res, out = _cli(tmp_path, js)
    assert res.exit_code == ExitCode.SUCCESS, res.output
    md = (out / "preflight" / "preflight.md").read_text()
    assert RULE_ID in md
    assert str(tmp_path / "bakta" / "sample.gbff") in md
    doc = json.loads((out / "preflight" / "preflight.json").read_text())
    (row,) = [c for c in doc["checks"] if c["rule_id"] == RULE_ID]
    assert row["level"] == "warn"


def test_ac6_cli_strict_exits_3(tmp_path: Path) -> None:
    res, _ = _cli(tmp_path, _bakta(tmp_path), "--strict")
    assert res.exit_code == ExitCode.WARNINGS_STRICT, res.output


def test_ac6_cli_with_gbff_has_no_warning(tmp_path: Path) -> None:
    res, out = _cli(tmp_path, _bakta(tmp_path, gbff="sample.gbff"), "--strict")
    doc = json.loads((out / "preflight" / "preflight.json").read_text())
    (row,) = [c for c in doc["checks"] if c["rule_id"] == RULE_ID]
    assert row["level"] == "ok"
    assert res.exit_code == ExitCode.SUCCESS, res.output
