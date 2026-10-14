"""Review follow-ups for the preflight report and environment rules (issue M1-11;
GIAP_SPEC.md §7.0, §9.5; spec ``docs/specs/M1-11-preflight-env-report.md`` "Lead decisions after
code review").

Pins: tool versions match the pinned major.minor with patch ≥ pin and pre-releases below the
release; databases match the pin to its written precision; numpy scalars in every field; colliding
dict keys raise; per-check exit code 1/2/0; ``routing`` closed in the schema; Markdown escaping of
``|``, backticks, ``<`` and whitespace; non-finite floats → null; floats to 6 significant digits.
"""

from __future__ import annotations

import copy
import json
import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import jsonschema
import numpy as np
import pytest

from giap.config import resolve_config
from giap.io.fasta import read_genome_fasta
from giap.model.coords import Replicon
from giap.preflight.checks import Check
from giap.preflight.environment import EnvProbe, Requirement, environment_checks, requirements
from giap.preflight.evidence import inputs_checks
from giap.preflight.genome import genome_checks, summarize_genome
from giap.preflight.report import (
    InputRef,
    PreflightReport,
    render_preflight_md,
    to_json_dict,
    write_preflight_json,
)

ROOT = Path(__file__).resolve().parents[3]
SCHEMA_PATH = ROOT / "src" / "giap" / "schemas" / "preflight.schema.json"
GIB = 1024**3
ROUTING: dict[str, Any] = {
    "n_proteins_long": 1,
    "n_proteins_peptide_path": 0,
    "n_proteins_dropped": 0,
    "n_peptides": 0,
    "abundance_used": False,
    "abundance_scale": "linear",
}


# ---------------------------------------------------------------- helpers


def validate(instance: object) -> None:
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema).validate(instance)


def chk(
    rule_id: str, level: str, value: object = 0, message: str = "", exit_code: int = 0
) -> Check:
    return Check(
        rule_id, rule_id.split(".")[0], level, f"rule of {rule_id}", value, message, exit_code
    )  # type: ignore[arg-type]


def report(
    checks: Sequence[Check],
    routing: Mapping[str, Any] | None = None,
    inputs: Sequence[InputRef] | None = None,
    replicons: Sequence[Replicon] | None = None,
) -> PreflightReport:
    return PreflightReport(
        giap_version="0.0.0",
        inputs=list(
            inputs if inputs is not None else [InputRef("genome", "genome.fa", "a" * 64, 10)]
        ),
        checks=list(checks),
        routing=dict(routing if routing is not None else ROUTING),
        evidence_tier="protein_list",
        replicons=list(
            replicons
            if replicons is not None
            else [Replicon("chr1", 400, "circular", "chromosome")]
        ),
        contig_map=[("chr1", "chr1")],
    )


def written(tmp_path: Path, r: PreflightReport) -> dict[str, Any]:
    out = tmp_path / "preflight.json"
    write_preflight_json(out, r)
    d = json.loads(out.read_text(encoding="utf-8"))
    assert isinstance(d, dict)
    return d


def value_of(d: Mapping[str, Any], rule_id: str) -> Any:
    return next(c["value"] for c in d["checks"] if c["rule_id"] == rule_id)


def found_all(reqs: Sequence[Requirement], **override: str | None) -> dict[str, str | None]:
    found: dict[str, str | None] = {r.name: r.min_version or "1.0" for r in reqs}
    found.update(override)
    return found


def tools(modules: Sequence[str], **override: str | None) -> Check:
    reqs = requirements(list(modules), "targeted")
    p = EnvProbe(
        found=found_all(reqs, **override), free_disk_bytes=500 * GIB, total_ram_bytes=64 * GIB
    )
    out = {c.rule_id: c for c in environment_checks(reqs, p, 0, 0)}
    return out["environment.tools"]


def _cells(line: str) -> list[str]:
    return [c.strip() for c in re.split(r"(?<!\\)\|", line.strip().strip("|"))]


def check_rows(md: str) -> list[list[str]]:
    lines = md.splitlines()
    head = next(
        i for i, ln in enumerate(lines) if _cells(ln) == ["Level", "Rule", "Value", "Message"]
    )
    rows = []
    for ln in lines[head + 2 :]:
        if not ln.strip().startswith("|"):
            break
        rows.append(_cells(ln))
    return rows


def md_row(c: Check) -> list[str]:
    rows = check_rows(render_preflight_md(report([c])))
    assert len(rows) == 1, rows
    return rows[0]


# ---------------------------------------------------------------- tool versions


@pytest.mark.parametrize("bakta", ["1.12.2", "1.12.3", "1.12.10"])
def test_bakta_patch_at_or_above_pin_accepted(bakta: str) -> None:
    c = tools(["core"], bakta=bakta, bakta_db="6.0")
    assert c.level == "ok", c
    assert c.exit_code == 0


@pytest.mark.parametrize(
    "bakta",
    [
        "2.0.0",  # newer major
        "1.13.0",  # newer minor
        "1.11.9",
        "1.12.1",  # patch below 1.12.2
        "1.12.0",
        "1.12.2rc1",  # pre-releases rank below the release
        "1.12.2-pre",
        "1.12.2a1",
        "1.12.2b2",
        "1.12.2.dev0",
    ],
)
def test_bakta_version_rejected(bakta: str) -> None:
    c = tools(["core"], bakta=bakta, bakta_db="6.0")
    assert c.level == "error", c
    assert c.exit_code == 2
    assert isinstance(c.value, dict)
    assert "bakta" in c.value.get("wrong_version", {}), c.value
    assert "bakta" in c.message


@pytest.mark.parametrize("genomad", ["2.0.0", "1.13.0", "1.11.5"])
def test_genomad_other_major_minor_rejected(genomad: str) -> None:
    c = tools(["core", "mobile"], bakta="1.12.2", bakta_db="6.0", genomad=genomad)
    assert c.level == "error", c
    assert isinstance(c.value, dict)
    assert "genomad" in c.value.get("wrong_version", {}), c.value


@pytest.mark.parametrize("db", ["6", "6.0", "6.1", "6.10"])
def test_bakta_db_major_6_accepted(db: str) -> None:
    c = tools(["core"], bakta="1.12.2", bakta_db=db)
    assert c.level == "ok", c


@pytest.mark.parametrize("db", ["5.1", "7.0", "60"])
def test_bakta_db_other_major_rejected(db: str) -> None:
    c = tools(["core"], bakta="1.12.2", bakta_db=db)
    assert c.level == "error", c
    assert isinstance(c.value, dict)
    assert "bakta_db" in c.value.get("wrong_version", {}), c.value


def test_genomad_db_pin_matched_to_minor() -> None:
    ok = tools(["core", "mobile"], bakta="1.12.2", bakta_db="6.0", genomad_db="1.9")
    assert ok.level == "ok", ok
    for bad in ("1.0", "1.8", "1.10", "2.9"):
        c = tools(["core", "mobile"], bakta="1.12.2", bakta_db="6.0", genomad_db=bad)
        assert c.level == "error", (bad, c)
        assert c.exit_code == 2
        assert isinstance(c.value, dict)
        assert "genomad_db" in c.value.get("wrong_version", {}), (bad, c.value)


# ---------------------------------------------------------------- numpy scalars in every field


def test_numpy_ints_in_replicon_length_and_input_size_write(tmp_path: Path) -> None:
    r = report(
        [chk("genome.parse", "ok", 1)],
        inputs=[InputRef("genome", "genome.fa", "a" * 64, np.int64(1234))],  # type: ignore[arg-type]
        replicons=[Replicon("chr1", np.int64(400), "circular", "chromosome")],  # type: ignore[arg-type]
    )
    d = written(tmp_path, r)
    validate(d)
    assert d["inputs"][0]["size"] == 1234
    assert type(d["inputs"][0]["size"]) is int
    assert d["replicons"][0]["length"] == 400
    assert type(d["replicons"][0]["length"]) is int


def test_numpy_scalars_in_routing_write(tmp_path: Path) -> None:
    routing = dict(ROUTING, n_peptides=np.int32(5), abundance_used=np.bool_(True))
    d = written(tmp_path, report([chk("genome.parse", "ok", 1)], routing=routing))
    validate(d)
    assert d["routing"]["n_peptides"] == 5
    assert d["routing"]["abundance_used"] is True


# ---------------------------------------------------------------- colliding keys


def test_value_keys_colliding_as_strings_raise(tmp_path: Path) -> None:
    r = report([chk("genome.gc", "info", {1: "a", "1": "b"})])
    with pytest.raises(ValueError):
        to_json_dict(r)
    with pytest.raises(ValueError):
        write_preflight_json(tmp_path / "preflight.json", r)


def test_nested_colliding_keys_raise() -> None:
    r = report([chk("genome.gc", "info", {"x": {np.int64(2): 0, "2": 1}})])
    with pytest.raises(ValueError):
        to_json_dict(r)


def test_routing_colliding_keys_raise() -> None:
    routing: dict[Any, Any] = dict(ROUTING)
    routing[7] = 1
    routing["7"] = 2
    with pytest.raises(ValueError):
        to_json_dict(report([chk("genome.parse", "ok", 1)], routing=routing))


def test_distinct_int_keys_still_written_as_strings() -> None:
    d = to_json_dict(report([chk("genome.gc", "info", {1: "a", 2: "b"})]))
    assert value_of(d, "genome.gc") == {"1": "a", "2": "b"}


# ---------------------------------------------------------------- per-check exit codes


def test_input_errors_carry_exit_code_1_from_real_rules(tmp_path: Path) -> None:
    g = tmp_path / "genome.fa"
    g.write_text(">a\n" + "ACGT" * 10 + "\n>a\n" + "ACGT" * 10 + "\n", encoding="utf-8")
    cfg = resolve_config()
    s = summarize_genome(read_genome_fasta(g), cfg.preflight.n_run_min_bp)
    checks = [*genome_checks(s, cfg), *inputs_checks(has_peptides=False, has_proteins=False)]
    errors = [c for c in checks if c.level == "error"]
    assert {c.rule_id for c in errors} >= {"genome.ids_unique", "inputs.evidence_present"}
    for c in checks:
        assert c.exit_code == (1 if c.level == "error" else 0), c
    env = tools(["core"], bakta="1.12.2", bakta_db=None)
    r = report([*checks, env])
    d = to_json_dict(r)
    validate(d)
    codes = {c["rule_id"]: c["exit_code"] for c in d["checks"]}
    assert codes["genome.ids_unique"] == 1
    assert codes["inputs.evidence_present"] == 1
    assert codes["environment.tools"] == 2
    assert r.exit_code == 1


def test_genome_parse_error_carries_exit_code_1() -> None:
    cfg = resolve_config()
    s = summarize_genome([], cfg.preflight.n_run_min_bp)
    parse = next(
        c for c in genome_checks(s, cfg, parse_error="bad FASTA") if c.rule_id == "genome.parse"
    )
    assert parse.level == "error"
    assert parse.exit_code == 1


@pytest.mark.parametrize("code", [0, 1, 2])
def test_schema_allows_exit_codes_0_1_2(code: int) -> None:
    validate(to_json_dict(report([chk("genome.parse", "error", 0, "x", code)])))


@pytest.mark.parametrize("code", [3, -1, 64])
def test_schema_rejects_other_exit_codes(code: int) -> None:
    d = copy.deepcopy(to_json_dict(report([chk("genome.parse", "ok", 1)])))
    d["checks"][0]["exit_code"] = code
    with pytest.raises(jsonschema.ValidationError):
        validate(d)


# ---------------------------------------------------------------- schema: routing closed


def test_schema_rejects_extra_routing_key() -> None:
    d = copy.deepcopy(to_json_dict(report([chk("genome.parse", "ok", 1)])))
    validate(d)
    d["routing"]["n_proteins_mystery"] = 3
    with pytest.raises(jsonschema.ValidationError):
        validate(d)


# ---------------------------------------------------------------- Markdown escaping


def test_md_pipe_in_message_is_escaped() -> None:
    row = md_row(chk("genome.gc", "warn", 0.5, "left|right"))
    assert len(row) == 4, row
    assert row[3] == "left\\|right"


def test_md_pipe_in_string_value_is_escaped() -> None:
    row = md_row(chk("genome.gc", "warn", "a|b|c", "m"))
    assert len(row) == 4, row
    assert row[2] == "a\\|b\\|c"


def test_md_whitespace_and_newlines_collapse() -> None:
    md = render_preflight_md(
        report([chk("genome.gc", "warn", 0.5, "line one\nline\t two\r\n  three ")])
    )
    rows = check_rows(md)
    assert len(rows) == 1
    assert rows[0][3] == "line one line two three"
    line = next(ln for ln in md.splitlines() if ln.startswith("| warn"))
    assert "line one line two three" in line


def test_md_backticks_escaped() -> None:
    row = md_row(chk("genome.gc", "warn", "`v`", "run `rm` now"))
    for cell in row[2:]:
        assert "`" in cell
        assert not re.search(r"(?<!\\)`", cell), cell


def test_md_angle_bracket_escaped() -> None:
    row = md_row(chk("genome.gc", "warn", "<b>", "see <script>alert(1)</script>"))
    for cell in row[2:]:
        assert "<b>" not in cell and "<script>" not in cell, cell
        assert not re.search(r"(?<!\\)<", cell), cell


# ---------------------------------------------------------------- non-finite floats


@pytest.mark.parametrize(
    "bad", [math.nan, math.inf, -math.inf, np.float64("nan"), np.float32("inf")]
)
def test_non_finite_floats_become_null_and_write(tmp_path: Path, bad: float) -> None:
    value = {"x": bad, "xs": [1.0, bad]}
    r = report([chk("genome.gc", "info", value), chk("genome.parse", "ok", bad)])
    d = written(tmp_path, r)
    validate(d)
    assert value_of(d, "genome.gc") == {"x": None, "xs": [1.0, None]}
    assert value_of(d, "genome.parse") is None
    md = render_preflight_md(r)
    assert "NaN" not in md and "Infinity" not in md


def test_non_finite_float_in_routing_becomes_null(tmp_path: Path) -> None:
    routing = dict(ROUTING, extra_ratio=math.nan)
    d = to_json_dict(report([chk("genome.parse", "ok", 1)], routing=routing))
    assert d["routing"]["extra_ratio"] is None
    json.dumps(d, allow_nan=False)


# ---------------------------------------------------------------- 6 significant digits


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        (0.123456789, 0.123457),
        (1 / 3, 0.333333),
        (1234567.0, 1234570.0),
        (2.0000004, 2.0),
        (np.float64(98.7654321), 98.7654),
        (np.float32(0.1), 0.1),
        (-0.000123456789, -0.000123457),
    ],
)
def test_floats_rounded_to_6_significant_digits(
    tmp_path: Path, raw: float, expected: float
) -> None:
    d = written(tmp_path, report([chk("genome.gc", "info", {"gc": raw, "list": [raw]})]))
    assert value_of(d, "genome.gc") == {"gc": expected, "list": [expected]}


def test_routing_floats_rounded_to_6_significant_digits() -> None:
    d = to_json_dict(report([chk("genome.parse", "ok", 1)], routing=dict(ROUTING, ratio=2 / 3)))
    assert d["routing"]["ratio"] == 0.666667


def test_md_value_uses_rounded_float() -> None:
    row = md_row(chk("genome.gc", "info", 0.123456789))
    assert row[2] == "0.123457"
