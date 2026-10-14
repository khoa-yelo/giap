"""``run_preflight``, ``exit_code`` and ``write_outputs`` (issue M1-12; GIAP_SPEC.md §5.1, §6.2,
§7.0) plus the real environment probe left over from M1-11 (question 1 and the "probe raises"
edge case; ADR-0005 item 3).

The real probe is selected with ``probe=None``; it reads ``.giap/capabilities.json`` and
``$GIAP_DB_DIR/db_manifest.json``. A missing manifest is info under ``mode="preflight"`` and an
environment error under ``mode="run"``; a probe failure (here: an unreadable manifest) is exactly
one ``environment.tools`` error (exit 2) in either mode.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import jsonschema
import pytest

import giap
from giap.config import GiapConfig, resolve_config
from giap.exitcodes import ExitCode
from giap.preflight.checks import Check
from giap.preflight.environment import GIB, EnvProbe, requirements
from giap.preflight.report import PreflightReport

ROOT = Path(__file__).resolve().parents[3]
SCHEMA_PATH = ROOT / "src" / "giap" / "schemas" / "preflight.schema.json"
LOW = ["preflight.min_n50_warn=1", "preflight.min_largest_contig_warn=1"]
ALL_MODULES = ["core", "mobile", "bgc", "function"]
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"


# ---------------------------------------------------------------- helpers


@pytest.fixture(autouse=True)
def _empty_db_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    db = tmp_path / "db"
    db.mkdir()
    monkeypatch.setenv("GIAP_DB_DIR", str(db))
    return db


def cfg(*extra: str) -> GiapConfig:
    return resolve_config(overrides=[*LOW, *extra])


def acgt(n: int) -> str:
    return ("ACGGTCAT" * (n // 8 + 1))[:n]


def fake_probe(**override: str | None) -> EnvProbe:
    found: dict[str, str | None] = {
        r.name: r.min_version or "1.0" for r in requirements(ALL_MODULES, "genome")
    }
    found.update(override)
    return EnvProbe(found=found, free_disk_bytes=500 * GIB, total_ram_bytes=64 * GIB)


def inputs(tmp_path: Path, genome: str | None = None, with_peptides: bool = True) -> Any:
    from giap.preflight.run import PreflightInputs

    g = tmp_path / "g.fa"
    g.write_text(genome if genome is not None else f">chr1 circular=true\n{acgt(150_000)}\n")
    p = tmp_path / "proteins.faa"
    p.write_text(f">P1\n{PROTEIN}\n>S1\nMKVLAAGIVGR\n")
    pep = None
    if with_peptides:
        pep = tmp_path / "peptides.tsv"
        pep.write_text("peptide\tprotein_ids\nPEPTIDEK\tP1\n")
    return PreflightInputs(genome=g, proteins=p, peptides=pep)


def run(tmp_path: Path, probe: EnvProbe | None, mode: str = "preflight", **kw: Any) -> Any:
    from giap.preflight.run import run_preflight

    return run_preflight(inputs(tmp_path, **kw), cfg(), probe, mode)  # type: ignore[arg-type]


def env_checks(report: PreflightReport) -> list[Check]:
    return [c for c in report.checks if c.input == "environment"]


def chk(rule_id: str, level: str, exit_code: int = 0) -> Check:
    input_name = rule_id.split(".")[0]
    return Check(rule_id, input_name, level, "rule", 0, "", exit_code)  # type: ignore[arg-type]


def report_of(checks: Sequence[Check]) -> PreflightReport:
    return PreflightReport(
        giap_version=giap.__version__,
        inputs=[],
        checks=list(checks),
        routing={
            "n_proteins_long": 0,
            "n_proteins_peptide_path": 0,
            "n_proteins_dropped": 0,
            "n_peptides": 0,
            "abundance_used": False,
            "abundance_scale": "linear",
        },
        evidence_tier="protein_list",
        replicons=[],
        contig_map=[],
    )


def validate(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema).validate(data)
    assert isinstance(data, dict)
    return data


# ---------------------------------------------------------------- PreflightInputs


def test_preflight_inputs_defaults_and_frozen(tmp_path: Path) -> None:
    from dataclasses import FrozenInstanceError

    from giap.preflight.run import PreflightInputs

    inp = PreflightInputs(genome=tmp_path / "g.fa")
    assert (inp.proteins, inp.peptides, inp.abundance, inp.annotation) == (None,) * 4
    with pytest.raises(FrozenInstanceError):
        inp.genome = tmp_path  # type: ignore[misc]


# ---------------------------------------------------------------- exit_code precedence


@pytest.mark.parametrize(
    ("checks", "strict", "expected"),
    [
        ([chk("genome.parse", "ok")], False, 0),
        ([chk("genome.parse", "ok")], True, 0),
        ([chk("genome.topology", "info"), chk("environment.resources", "info")], True, 0),
        ([chk("genome.contiguity", "warn")], False, 0),
        ([chk("genome.contiguity", "warn")], True, 3),
        ([chk("environment.resources", "warn")], True, 3),
        ([chk("genome.parse", "error")], False, 1),
        ([chk("genome.parse", "error"), chk("genome.gc", "warn")], True, 1),
        ([chk("environment.tools", "error", 2)], False, 2),
        ([chk("environment.tools", "error", 2), chk("genome.gc", "warn")], True, 2),
        ([chk("environment.tools", "error", 2), chk("peptides.schema", "error")], True, 1),
        ([chk("inputs.evidence_present", "error"), chk("environment.tools", "error", 2)], False, 1),
    ],
)
def test_exit_code_precedence_1_over_2_over_3(
    checks: list[Check], strict: bool, expected: int
) -> None:
    from giap.preflight.run import exit_code

    code = exit_code(report_of(checks), strict)
    assert code == expected
    assert isinstance(code, ExitCode)


# ------------------------------------------------ run_preflight with an injected probe


def test_run_preflight_clean_with_fake_probe(tmp_path: Path) -> None:
    from giap.preflight.run import exit_code

    report, _norm = run(tmp_path, fake_probe())
    assert isinstance(report, PreflightReport)
    assert [c for c in report.checks if c.level == "error"] == []
    tools = [c for c in env_checks(report) if c.rule_id == "environment.tools"]
    assert len(tools) == 1
    assert tools[0].level == "ok"
    assert exit_code(report, strict=False) == ExitCode.SUCCESS
    ids = {c.rule_id for c in report.checks}
    assert {"inputs.evidence_present", "genome.parse", "proteins.ids", "peptides.schema"} <= ids
    assert report.contig_map == [("chr1", "chr1")]
    assert [r.kind for r in report.replicons] == ["chromosome"]


def test_run_preflight_missing_tool_is_exit_2(tmp_path: Path) -> None:
    from giap.preflight.run import exit_code

    report, _ = run(tmp_path, fake_probe(bakta_db=None))
    tools = [c for c in env_checks(report) if c.rule_id == "environment.tools"]
    assert len(tools) == 1
    assert tools[0].level == "error"
    assert tools[0].exit_code == 2
    assert exit_code(report, strict=False) == ExitCode.ENVIRONMENT
    assert exit_code(report, strict=True) == ExitCode.ENVIRONMENT


def test_run_preflight_uses_enabled_modules_only(tmp_path: Path) -> None:
    from giap.preflight.run import PreflightInputs, run_preflight

    inp = inputs(tmp_path)
    assert isinstance(inp, PreflightInputs)
    probe = fake_probe(genomad=None, antismash=None, go=None)
    report, _ = run_preflight(inp, cfg("modules=[core]"), probe, "preflight")
    tools = [c for c in env_checks(report) if c.rule_id == "environment.tools"]
    assert [c.level for c in tools] == ["ok"]


def test_run_preflight_routing_counts_short_proteins_as_peptides(tmp_path: Path) -> None:
    report, _ = run(tmp_path, fake_probe())
    assert report.routing["n_proteins_long"] == 1
    assert report.routing["n_proteins_peptide_path"] == 1
    assert report.routing["n_proteins_dropped"] == 0
    assert report.routing["n_peptides"] == 2  # PEPTIDEK + MKVLAAGIVGR from S1


def test_run_preflight_genome_error_is_exit_1(tmp_path: Path) -> None:
    from giap.preflight.run import exit_code

    report, _ = run(tmp_path, fake_probe(), genome=">c\nACGU\n")
    parse = [c for c in report.checks if c.rule_id == "genome.parse"]
    assert parse[0].level == "error"
    assert exit_code(report, strict=True) == ExitCode.INVALID_INPUT


# ---------------------------------------------------------------- the real probe (probe=None)


def test_real_probe_missing_manifest_is_info_under_preflight(tmp_path: Path) -> None:
    from giap.preflight.run import exit_code

    report, _ = run(tmp_path, None, mode="preflight")
    env = env_checks(report)
    assert env, "environment checks must be recorded"
    assert all(c.level in ("ok", "info") for c in env), env
    assert any("db_manifest" in c.message + json.dumps(c.value) for c in env)
    assert exit_code(report, strict=True) == ExitCode.SUCCESS


def test_real_probe_missing_manifest_is_error_under_run(tmp_path: Path) -> None:
    from giap.preflight.run import exit_code

    report, _ = run(tmp_path, None, mode="run")
    tools = [c for c in env_checks(report) if c.rule_id == "environment.tools"]
    assert len(tools) == 1
    assert tools[0].level == "error"
    assert tools[0].exit_code == 2
    assert "db_manifest" in tools[0].message + json.dumps(tools[0].value)
    assert exit_code(report, strict=False) == ExitCode.ENVIRONMENT


@pytest.mark.parametrize("mode", ["preflight", "run"])
@pytest.mark.parametrize("content", ["{not json", "[1, 2, 3]\n"], ids=["invalid", "not_object"])
def test_real_probe_failure_is_one_environment_tools_error(
    tmp_path: Path, _empty_db_dir: Path, mode: str, content: str
) -> None:
    from giap.preflight.run import exit_code

    (_empty_db_dir / "db_manifest.json").write_text(content)
    report, _ = run(tmp_path, None, mode=mode)
    tools = [c for c in env_checks(report) if c.rule_id == "environment.tools"]
    assert len(tools) == 1
    assert tools[0].level == "error"
    assert tools[0].exit_code == 2
    assert tools[0].message
    assert [c for c in report.checks if c.level == "error"] == tools
    assert exit_code(report, strict=False) == ExitCode.ENVIRONMENT


# ---------------------------------------------------------------- write_outputs


def test_write_outputs_clean_report_writes_report_and_work_files(tmp_path: Path) -> None:
    from giap.preflight.run import write_outputs

    report, norm = run(tmp_path, fake_probe())
    out = tmp_path / "out"
    write_outputs(out, report, norm)
    data = validate(out / "preflight" / "preflight.json")
    assert data["summary"]["n_error"] == 0
    assert (out / "preflight" / "preflight.md").read_text().startswith("# GIAP preflight")
    work = sorted(p.name for p in (out / "work" / "preflight").iterdir())
    assert work == sorted(
        [
            "contig_map.tsv",
            "genome.norm.fa",
            "peptide_runs.parquet",
            "peptides.norm.parquet",
            "proteins.long.faa",
            "replicons.tsv",
        ]
    )
    assert sorted(p.name for p in (out / "preflight").iterdir()) == [
        "preflight.json",
        "preflight.md",
    ]


def test_write_outputs_error_report_writes_no_work_files(tmp_path: Path) -> None:
    from giap.preflight.run import write_outputs

    report, norm = run(tmp_path, fake_probe(), genome=">c\nACGU\n")
    out = tmp_path / "out"
    write_outputs(out, report, norm)
    data = validate(out / "preflight" / "preflight.json")
    assert data["summary"]["n_error"] >= 1
    assert (out / "preflight" / "preflight.md").is_file()
    assert not (out / "work" / "preflight").exists()


def test_write_outputs_environment_error_still_writes_report(tmp_path: Path) -> None:
    from giap.preflight.run import write_outputs

    report, norm = run(tmp_path, fake_probe(bakta=None))
    out = tmp_path / "out"
    write_outputs(out, report, norm)
    data = validate(out / "preflight" / "preflight.json")
    assert data["summary"]["n_error"] == 1
    assert not (out / "work" / "preflight").exists()


def test_write_outputs_twice_is_byte_identical(tmp_path: Path) -> None:
    from giap.preflight.run import write_outputs

    report, norm = run(tmp_path, fake_probe())
    a, b = tmp_path / "a", tmp_path / "b"
    write_outputs(a, report, norm)
    write_outputs(b, report, norm)
    files = sorted(p.relative_to(a) for p in a.rglob("*") if p.is_file())
    assert files
    for rel in files:
        assert (a / rel).read_bytes() == (b / rel).read_bytes(), rel
