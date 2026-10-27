"""``preflight.json`` schema and writer, ``preflight.md`` summary, report exit code
(issue M1-11; GIPPA_SPEC.md §7.0 Outputs, §9.5, §6.6, §6.7; ADR-0005 item 3). T-PRE-25.

Reports are built from real ``Check`` lists produced by the genome, evidence and environment
rules on tiny in-memory or ``tmp_path`` inputs. ``jsonschema`` is a test-only dependency here.
"""

from __future__ import annotations

import copy
import json
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import jsonschema
import numpy as np
import pytest

import gippa
from gippa.config import GippaConfig, resolve_config
from gippa.io.fasta import read_genome_fasta
from gippa.io.protein_fasta import read_protein_fasta
from gippa.preflight.checks import Check
from gippa.preflight.environment import EnvProbe, environment_checks, requirements
from gippa.preflight.evidence import candidate_tier, inputs_checks, protein_checks, tier_check
from gippa.preflight.genome import (
    classify_replicons,
    genome_checks,
    resolve_topology,
    summarize_genome,
)
from gippa.preflight.report import (
    InputRef,
    PreflightReport,
    render_preflight_md,
    to_json_dict,
    write_preflight_json,
)

ROOT = Path(__file__).resolve().parents[3]
SCHEMA_PATH = ROOT / "src" / "gippa" / "schemas" / "preflight.schema.json"
GIB = 1024**3
TIERS = ("decoy_picked", "mayu", "two_peptide", "protein_list")
TOP_KEYS = {
    "schema_version",
    "gippa_version",
    "inputs",
    "checks",
    "summary",
    "evidence_tier",
    "routing",
    "replicons",
    "contig_map",
}
SMALL = [
    "preflight.min_total_bp_warn=1",
    "preflight.min_n50_warn=1",
    "preflight.min_largest_contig_warn=1",
    "preflight.small_replicon_bp=1",
]
RANK = {"error": 0, "warn": 1, "info": 2, "ok": 3}
TIMESTAMP = re.compile(r"\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}")


# ---------------------------------------------------------------- helpers


def schema() -> dict[str, Any]:
    data = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def validate(instance: object) -> None:
    jsonschema.Draft202012Validator(schema()).validate(instance)


def cfg() -> GippaConfig:
    return resolve_config(overrides=SMALL)


def acgt(n: int) -> str:
    return ("ACGT" * (n // 4 + 1))[:n]


def env_ok() -> list[Check]:
    reqs = requirements(["core"], "targeted")
    found: dict[str, str | None] = {r.name: r.min_version or "1.0" for r in reqs}
    found.update(bakta="1.12.2", bakta_db="6.0")
    p = EnvProbe(found=found, free_disk_bytes=500 * GIB, total_ram_bytes=64 * GIB)
    return environment_checks(reqs, p, 0, 0)


def env_missing_db() -> list[Check]:
    reqs = requirements(["core"], "targeted")
    found: dict[str, str | None] = {r.name: r.min_version or "1.0" for r in reqs}
    found.update(bakta="1.12.2", bakta_db=None)
    p = EnvProbe(found=found, free_disk_bytes=500 * GIB, total_ram_bytes=64 * GIB)
    return environment_checks(reqs, p, 0, 0)


def build_report(
    tmp_path: Path, genome_text: str | None = None, env: Sequence[Check] | None = None
) -> PreflightReport:
    """A report from a tiny genome + protein file through the real rules."""
    config = cfg()
    gpath = tmp_path / "genome.fa"
    gpath.write_text(
        genome_text
        if genome_text is not None
        else f">chr1 circular=true\n{acgt(400)}\n>p|1\n{acgt(80)}\n",
        encoding="utf-8",
    )
    ppath = tmp_path / "proteins.faa"
    ppath.write_text(">P1\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA\n")
    records = read_genome_fasta(gpath)
    s = summarize_genome(records, config.preflight.n_run_min_bp)
    topo = resolve_topology(s, config.annotation.bakta.complete) if s.ids else ()
    replicons = classify_replicons(s, topo, config.preflight.small_replicon_bp) if s.ids else []
    proteins = read_protein_fasta(ppath)
    tier = candidate_tier(None, has_proteins=True)
    checks = [
        *inputs_checks(has_peptides=False, has_proteins=True),
        *genome_checks(s, config),
        *protein_checks(proteins, None, (7, 40)),
        tier_check(tier),
        *(env if env is not None else env_ok()),
    ]
    inputs = [
        InputRef("genome", "genome.fa", "a" * 64, gpath.stat().st_size),
        InputRef("proteins", "proteins.faa", "b" * 64, ppath.stat().st_size),
    ]
    return PreflightReport(
        gippa_version=gippa.__version__,
        inputs=inputs,
        checks=checks,
        routing={
            "n_proteins_long": 1,
            "n_proteins_peptide_path": 0,
            "n_proteins_dropped": 0,
            "n_peptides": 0,
            "abundance_used": False,
            "abundance_scale": "linear",
        },
        evidence_tier=tier,
        replicons=replicons,
        contig_map=list(zip(s.ids, s.original_ids, strict=True)),
    )


def chk(
    rule_id: str, level: str, value: object = 0, message: str = "", exit_code: int = 0
) -> Check:
    input_name = rule_id.split(".")[0]
    return Check(rule_id, input_name, level, f"rule of {rule_id}", value, message, exit_code)  # type: ignore[arg-type]


# ---------------------------------------------------------------- schema itself


def test_preflight_schema_is_valid_draft_2020_12() -> None:
    s = schema()
    assert s["$schema"] == "https://json-schema.org/draft/2020-12/schema"
    assert s.get("type") == "object"
    jsonschema.Draft202012Validator.check_schema(s)


def test_preflight_schema_top_level_contract() -> None:
    s = schema()
    assert set(s["required"]) >= TOP_KEYS
    assert s["additionalProperties"] is False


# ---------------------------------------------------------------- T-PRE-25 validation


def test_T_PRE_25_pass(tmp_path: Path) -> None:
    r = build_report(tmp_path)
    d = to_json_dict(r)
    validate(d)
    assert set(d) == TOP_KEYS
    assert d["summary"]["n_error"] == 0
    assert d["gippa_version"] == gippa.__version__
    assert d["evidence_tier"] == "protein_list"


def test_T_PRE_25_fail_unknown_level(tmp_path: Path) -> None:
    d = to_json_dict(build_report(tmp_path))
    bad = copy.deepcopy(d)
    bad["checks"][0]["level"] = "fatal"
    with pytest.raises(jsonschema.ValidationError):
        validate(bad)


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("checks", 0, "rule_id"), "Genome.Parse"),
        (("checks", 0, "rule_id"), "genome"),
        (("checks", 0, "input"), "metadata"),
        (("checks", 0, "exit_code"), 3),
        (("evidence_tier",), "spectral_count"),
        (("summary", "n_error"), "0"),
        (("routing", "abundance_used"), "no"),
    ],
)
def test_T_PRE_25_fail_field_constraints(
    tmp_path: Path, path: tuple[str | int, ...], value: object
) -> None:
    bad = copy.deepcopy(to_json_dict(build_report(tmp_path)))
    target: Any = bad
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises(jsonschema.ValidationError):
        validate(bad)


@pytest.mark.parametrize("key", sorted(TOP_KEYS))
def test_T_PRE_25_fail_missing_top_level_key(tmp_path: Path, key: str) -> None:
    bad = copy.deepcopy(to_json_dict(build_report(tmp_path)))
    del bad[key]
    with pytest.raises(jsonschema.ValidationError):
        validate(bad)


def test_T_PRE_25_fail_extra_top_level_key(tmp_path: Path) -> None:
    bad = copy.deepcopy(to_json_dict(build_report(tmp_path)))
    bad["timestamp"] = "2026-10-14T00:00:00"
    with pytest.raises(jsonschema.ValidationError):
        validate(bad)


@pytest.mark.parametrize(
    "key",
    [
        "n_proteins_long",
        "n_proteins_peptide_path",
        "n_proteins_dropped",
        "n_peptides",
        "abundance_used",
        "abundance_scale",
    ],
)
def test_T_PRE_25_fail_missing_routing_key(tmp_path: Path, key: str) -> None:
    bad = copy.deepcopy(to_json_dict(build_report(tmp_path)))
    del bad["routing"][key]
    with pytest.raises(jsonschema.ValidationError):
        validate(bad)


def test_T_PRE_25_report_with_errors_validates(tmp_path: Path) -> None:
    """Duplicate contig IDs (input error) plus a missing DB (environment error)."""
    r = build_report(tmp_path, f">a\n{acgt(40)}\n>a\n{acgt(40)}\n", env=env_missing_db())
    d = to_json_dict(r)
    validate(d)
    assert d["summary"]["n_error"] == 2


# ---------------------------------------------------------------- JSON content


def test_every_evaluated_rule_present(tmp_path: Path) -> None:
    r = build_report(tmp_path)
    d = to_json_dict(r)
    assert sorted(c["rule_id"] for c in d["checks"]) == sorted(c.rule_id for c in r.checks)
    assert any(c["level"] == "ok" for c in d["checks"])
    ids = {c["rule_id"] for c in d["checks"]}
    assert {"genome.parse", "environment.tools", "environment.resources"} <= ids


def test_check_fields_round_trip(tmp_path: Path) -> None:
    r = build_report(tmp_path)
    d = to_json_dict(r)
    by_id = {c["rule_id"]: c for c in d["checks"]}
    for c in r.checks:
        j = by_id[c.rule_id]
        assert j["input"] == c.input
        assert j["level"] == c.level
        assert j["rule"] == c.rule
        assert j["message"] == c.message
        assert j["exit_code"] == c.exit_code
        assert j["value"] == json.loads(json.dumps(c.value))


def test_summary_counts_levels(tmp_path: Path) -> None:
    r = build_report(tmp_path, f">a\n{acgt(40)}\n>a\n{acgt(40)}\n", env=env_missing_db())
    d = to_json_dict(r)
    levels = [c.level for c in r.checks]
    assert d["summary"] == {
        "n_error": levels.count("error"),
        "n_warn": levels.count("warn"),
        "n_info": levels.count("info"),
        "n_ok": levels.count("ok"),
    }


def test_replicons_and_contig_map(tmp_path: Path) -> None:
    r = build_report(tmp_path)
    d = to_json_dict(r)
    assert d["replicons"] == [
        {"id": x.id, "length": x.length, "topology": x.topology, "kind": x.kind}
        for x in r.replicons
    ]
    assert d["replicons"][0] == {
        "id": "chr1",
        "length": 400,
        "topology": "circular",
        "kind": "chromosome",
    }
    pairs = [
        tuple(p) if isinstance(p, list) else (p["id"], p["original_id"]) for p in d["contig_map"]
    ]
    assert pairs == [("chr1", "chr1"), ("p_1", "p|1")]


def test_inputs_hold_file_names_only(tmp_path: Path) -> None:
    d = to_json_dict(build_report(tmp_path))
    assert len(d["inputs"]) == 2
    genome = next(i for i in d["inputs"] if i["role"] == "genome")
    assert "genome.fa" in genome.values()
    assert genome["sha256"] == "a" * 64
    assert genome["size"] == (tmp_path / "genome.fa").stat().st_size
    for entry in d["inputs"]:
        for v in entry.values():
            if isinstance(v, str):
                assert "/" not in v, entry


def test_routing_copied(tmp_path: Path) -> None:
    r = build_report(tmp_path)
    assert to_json_dict(r)["routing"] == dict(r.routing)


# ---------------------------------------------------------------- writer


def test_write_twice_byte_identical_sorted_no_paths_or_timestamps(tmp_path: Path) -> None:
    r = build_report(tmp_path)
    a, b = tmp_path / "out_a" / "preflight.json", tmp_path / "out_b" / "preflight.json"
    a.parent.mkdir()
    b.parent.mkdir()
    write_preflight_json(a, r)
    write_preflight_json(b, r)
    raw = a.read_bytes()
    assert raw == b.read_bytes()
    text = raw.decode("utf-8")
    assert text.endswith("}\n")
    assert not text.endswith("\n\n")
    assert text.splitlines()[1].startswith('  "')  # indent 2
    assert str(tmp_path) not in text
    assert str(Path.home()) not in text
    assert not TIMESTAMP.search(text)
    assert json.loads(text) == to_json_dict(r)

    def hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        keys = [k for k, _ in pairs]
        assert keys == sorted(keys), keys
        return dict(pairs)

    validate(json.loads(text, object_pairs_hook=hook))


def test_written_bytes_independent_of_construction(tmp_path: Path) -> None:
    """Two reports built separately from the same inputs write the same bytes."""
    (tmp_path / "x").mkdir()
    (tmp_path / "y").mkdir()
    r1 = build_report(tmp_path / "x")
    r2 = build_report(tmp_path / "y")
    write_preflight_json(tmp_path / "x" / "preflight.json", r1)
    write_preflight_json(tmp_path / "y" / "preflight.json", r2)
    assert (tmp_path / "x" / "preflight.json").read_bytes() == (
        tmp_path / "y" / "preflight.json"
    ).read_bytes()


def test_numpy_scalars_in_values_are_converted(tmp_path: Path) -> None:
    r = build_report(tmp_path)
    extra = Check(
        "genome.gc",
        "genome",
        "ok",
        "gc",
        {"gc": np.float64(0.5), "n": np.int64(3), "flag": np.bool_(True)},  # type: ignore[dict-item]
    )
    r2 = PreflightReport(
        gippa_version=r.gippa_version,
        inputs=r.inputs,
        checks=[c for c in r.checks if c.rule_id != "genome.gc"] + [extra],
        routing=r.routing,
        evidence_tier=r.evidence_tier,
        replicons=r.replicons,
        contig_map=r.contig_map,
    )
    out = tmp_path / "preflight.json"
    write_preflight_json(out, r2)
    d = json.loads(out.read_text(encoding="utf-8"))
    validate(d)
    gc = next(c for c in d["checks"] if c["rule_id"] == "genome.gc")
    assert gc["value"] == {"gc": 0.5, "n": 3, "flag": True}


# ---------------------------------------------------------------- exit code


def test_exit_code_clean_is_0(tmp_path: Path) -> None:
    assert build_report(tmp_path).exit_code == 0


def test_exit_code_warnings_only_is_0(tmp_path: Path) -> None:
    r = build_report(tmp_path, f">a\n{acgt(40)}NNNNNNNNNNNN{acgt(40)}\n")
    assert any(c.level == "warn" for c in r.checks)
    assert not any(c.level == "error" for c in r.checks)
    assert r.exit_code == 0


def test_exit_code_environment_error_is_2(tmp_path: Path) -> None:
    assert build_report(tmp_path, env=env_missing_db()).exit_code == 2


def test_exit_code_input_error_beats_environment_error(tmp_path: Path) -> None:
    r = build_report(tmp_path, f">a\n{acgt(40)}\n>a\n{acgt(40)}\n", env=env_missing_db())
    assert r.exit_code == 1
    r_input_only = build_report(tmp_path, f">a\n{acgt(40)}\n>a\n{acgt(40)}\n")
    assert r_input_only.exit_code == 1


# ---------------------------------------------------------------- preflight.md


def _md_report(tmp_path: Path) -> PreflightReport:
    r = build_report(tmp_path)
    checks = [
        chk("genome.gc", "info", 0.5),
        chk("genome.contig_count", "warn", 12, "more than 10 contigs"),
        chk("genome.parse", "ok", 1),
        chk("environment.tools", "error", {"missing": ["bakta_db"]}, "missing: bakta_db", 2),
        chk("genome.topology", "ok", []),
    ]
    return PreflightReport(
        gippa_version=r.gippa_version,
        inputs=r.inputs,
        checks=checks,
        routing=r.routing,
        evidence_tier=r.evidence_tier,
        replicons=r.replicons,
        contig_map=r.contig_map,
    )


def _cells(line: str) -> list[str]:
    """Cells of a Markdown table row; ``\\|`` is an escaped pipe inside a cell."""
    return [c.strip() for c in re.split(r"(?<!\\)\|", line.strip().strip("|"))]


def _table_rows(md: str) -> list[list[str]]:
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


def test_md_summary_line_counts(tmp_path: Path) -> None:
    md = render_preflight_md(_md_report(tmp_path))
    assert "1 error, 1 warning" in md
    assert md.index("1 error, 1 warning") < md.index("| Level")


def test_md_table_ordered_error_warn_info_ok(tmp_path: Path) -> None:
    rows = _table_rows(render_preflight_md(_md_report(tmp_path)))
    assert len(rows) == 5
    levels = [r[0].lower() for r in rows]
    assert levels == sorted(levels, key=RANK.__getitem__)
    assert levels[:2] == ["error", "warn"]
    assert "environment.tools" in rows[0][1]
    assert "genome.contig_count" in rows[1][1]
    assert "bakta_db" in rows[0][2]
    assert "more than 10 contigs" in rows[1][3]


def test_md_from_real_report_has_no_paths_and_is_deterministic(tmp_path: Path) -> None:
    r = build_report(tmp_path, f">a\n{acgt(40)}\n>a\n{acgt(40)}\n", env=env_missing_db())
    md = render_preflight_md(r)
    assert md == render_preflight_md(r)
    assert md.endswith("\n")
    assert str(tmp_path) not in md
    assert "2 errors" in md
    assert len(_table_rows(md)) == len(r.checks)
