"""1.0 supports genetic code 11 only (Q-K30 default (b); table 4 and others: backlog, 1.1).

Bakta, the six-frame index and codon usage all assume code 11, so any other
``annotation.translation_table`` must stop the run in preflight with an input error (exit 1),
never a warning and never a silent code-11 run. Rule id ``config.translation_table``; the
message names the value, the setting and "only genetic code 11 is supported in this release".
``gippa run`` must not start Nextflow (a fake ``nextflow`` on PATH records any launch).
"""

from __future__ import annotations

import json
import os
import stat
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner

from gippa.config import resolve_config
from gippa.exitcodes import ExitCode

ROOT = Path(__file__).resolve().parents[3]
SCHEMA_PATH = ROOT / "src" / "gippa" / "schemas" / "preflight.schema.json"
RULE = "config.translation_table"
LOW = ["preflight.min_n50_warn=1", "preflight.min_largest_contig_warn=1"]
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    for name in ("db", "repo"):
        (tmp_path / name).mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db"))
    monkeypatch.setenv("GIPPA_REPO", str(tmp_path / "repo"))
    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "repo"))  # never the real $GIPPA_ROOT
    monkeypatch.delenv("NXF_APPTAINER_CACHEDIR", raising=False)


def _files(tmp_path: Path) -> tuple[Path, Path]:
    g = tmp_path / "g.fa"
    g.write_text(">chr1 circular=true\n" + ("ACGGTCAT" * 18_750) + "\n")
    p = tmp_path / "proteins.faa"
    p.write_text(f">P1\n{PROTEIN}\n")
    return g, p


def _report(tmp_path: Path, *sets: str) -> Any:
    from gippa.preflight.run import PreflightInputs, run_preflight

    g, p = _files(tmp_path)
    cfg = resolve_config(overrides=[*LOW, *sets])
    return run_preflight(PreflightInputs(genome=g, proteins=p), cfg, None, "preflight")[0]


@pytest.mark.parametrize("table", [4, 25])
def test_non_11_translation_table_is_a_preflight_error(tmp_path: Path, table: int) -> None:
    report = _report(tmp_path, f"annotation.translation_table={table}")
    hits = [c for c in report.checks if c.rule_id == RULE]
    assert len(hits) == 1, [c.rule_id for c in report.checks]
    c = hits[0]
    assert (c.level, c.exit_code) == ("error", int(ExitCode.INVALID_INPUT))
    assert c.input != "environment"
    assert str(table) in c.message
    assert "annotation.translation_table" in c.message
    assert "only genetic code 11 is supported in this release" in c.message.lower()
    assert report.exit_code == int(ExitCode.INVALID_INPUT)


def test_default_translation_table_11_passes(tmp_path: Path) -> None:
    report = _report(tmp_path)
    hits = [c for c in report.checks if c.rule_id == RULE]
    assert [c.level for c in hits] == ["ok"]
    assert report.exit_code == int(ExitCode.SUCCESS)


def test_gippa_preflight_exits_invalid_input(tmp_path: Path) -> None:
    import jsonschema

    from gippa.cli import cli

    g, p = _files(tmp_path)
    out = tmp_path / "OUT"
    base = ["preflight", "--profile", "none", "--genome", str(g), "--proteins", str(p)]
    base += ["--outdir", str(out)]
    base += [a for s in LOW for a in ("--set", s)]
    ok = CliRunner().invoke(cli, base)
    assert ok.exit_code == ExitCode.SUCCESS, ok.output
    res = CliRunner().invoke(cli, [*base, "--set", "annotation.translation_table=4"])
    assert res.exit_code == ExitCode.INVALID_INPUT, res.output
    assert RULE in res.output
    data = json.loads((out / "preflight" / "preflight.json").read_text(encoding="utf-8"))
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema).validate(data)
    assert any(c["rule_id"] == RULE and c["level"] == "error" for c in data["checks"])


def test_gippa_run_exits_invalid_input_without_launching_nextflow(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.cli import cli

    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text('#!/bin/sh\ntouch "$FAKE_NF_MARK"\nexit 0\n')
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    mark = tmp_path / "nextflow_launched"
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_MARK", str(mark))
    g, p = _files(tmp_path)
    argv = ["run", "--genome", str(g), "--proteins", str(p), "--outdir", str(tmp_path / "OUT"),
            "--profile", "apptainer", "--set", "annotation.translation_table=4"]  # fmt: skip
    res = CliRunner().invoke(cli, [*argv, *(a for s in LOW for a in ("--set", s))])
    assert res.exit_code == ExitCode.INVALID_INPUT, res.output
    assert RULE in res.output
    assert not mark.exists(), "gippa run launched nextflow despite a preflight error"


def test_inputs_page_states_genetic_code_11() -> None:
    text = (ROOT / "docs" / "user-guide" / "inputs.md").read_text(encoding="utf-8")
    assert "genetic code 11" in text.lower()
    assert "annotation.translation_table" in text
