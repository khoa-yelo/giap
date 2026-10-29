"""M11-16d: tool plumbing for genetic codes 4 and 25 (docs/specs/M11-16-translation-tables-4-25.md
items 3.1-3.4 and 7; examples 7, 8 and 13; AC 9, 10 (parse part) and 13).

Interface pinned here (test-author reading, beyond the spec text):

- ``bakta_command(..., translation_table=t)`` (keyword, as in the spec signature). At 11 the argv
  is byte-identical to today (D2, Khoa 2026-10-28); at 4/25 ``--translation-table t`` sits right
  before ``--replicons`` (example 8), after ``--complete`` when that is passed.
- ``miniprot_index_cmd(genome, mpi, threads, table=t)`` gives ``-T t``; the align command never
  carries ``-T`` (§7.1.1 gotcha).
- ``effective_tblastn_args(args, table)`` in ``gippa.mapping.tblastn`` (spec 3.3);
  ``tblastn_cmd(db, proteins, threads, args, table=t)`` splices the effective string. The parity
  case table lives in ``workflow/tests/data/tblastn/db_gencode_cases.json``; the same cases run
  through the Nextflow closure in ``workflow/tests/tblastn.nf.test``.
- D1 option A (Khoa, 2026-10-28): the shipped default becomes ``"-seg no -evalue 1e-5"``, and the
  16a bridge (the shipped-default exemption in ``gippa.preflight.tool_genetic_code``) is gone, so
  the old default string at table 4 is refused like any other mismatching ``-db_gencode``.
- ``parse_bakta_json(path, ..., translation_table=t)`` (keyword) raises
  ``gippa.annotate.bakta.AnnotationTableMismatch`` (a ``ValueError``) when the JSON's
  ``genome.translation_table`` differs from ``t``; a JSON without the key is read at ``t`` and the
  parser logs ``annotation.table_unrecorded`` at INFO (as a record's ``event`` attribute or in
  its message). ``gippa _step parse_bakta`` then exits ``INVALID_INPUT`` (3.4).
- ``annotation_meta.json`` gains ``translation_table`` (the config value); the manifest's
  ``annotation.translation_table`` is read from it, or from the resolved config when the meta has
  no such key. ``OUTPUT_SCHEMA_VERSION`` goes from ``1.0`` to ``1.1`` (additive, M11-14h).
- The report reads the table from the run (the tests set the config and the manifest alike): one
  Methods sentence on every run; at 4/25 the summary opens with the ``data-genetic-code`` note,
  after any ``SUMMARY_WARNINGS`` alert (D8).
"""

from __future__ import annotations

import json
import logging
import re
import shlex
from pathlib import Path
from typing import Any

import jsonschema
import pytest
import yaml
from helpers import manifest_data as md
from helpers.report_content import html, page, section_text, write_manifest
from helpers.report_fixtures import run_out
from hypothesis import given, settings
from hypothesis import strategies as st
from unit.test_m11_16a_preflight import SMALL, TOOL_RULE, _only, _run

from gippa.exitcodes import ExitCode

REPO = Path(__file__).resolve().parents[2]
CASES = REPO / "workflow" / "tests" / "data" / "tblastn" / "db_gencode_cases.json"
DEFAULTS_YAML = REPO / "src" / "gippa" / "config" / "defaults.yaml"
CONFIG_SCHEMA = REPO / "src" / "gippa" / "schemas" / "config.schema.json"
MODULES = REPO / "workflow" / "modules" / "local"
BASE_CONFIG = REPO / "workflow" / "conf" / "base.config"

TABLES = (4, 11, 25)
NEW_DEFAULT = "-seg no -evalue 1e-5"
OLD_DEFAULT = "-db_gencode 11 -seg no -evalue 1e-5"
NAMES = {
    11: "Bacterial, Archaeal and Plant Plastid",
    4: "Mycoplasma and Spiroplasma",
    25: "Candidate Division SR1 and Gracilibacteria",
}
TGA = {4: "W", 25: "G"}


def _cases() -> list[dict[str, Any]]:
    data = json.loads(CASES.read_text(encoding="utf-8"))
    cases: list[dict[str, Any]] = data["cases"]
    return cases


def _strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"(?m)^\s*//.*$", "", text)


def _script(tool: str) -> str:
    text = _strip_comments((MODULES / tool / "main.nf").read_text(encoding="utf-8"))
    m = re.search(r"\bscript\s*:(.*?)\bstub\s*:", text, flags=re.S)
    assert m, f"{tool}: no script: block followed by stub:"
    return m.group(1)


# ================================================================ 3.1 Bakta argv

GENOME = Path("/w/genome.norm.fa")
OUT = Path("/w/bakta")
DB = Path("/db/bakta/db")
REPLICONS = Path("/w/replicons.tsv")
REGIONS = Path("/w/regions.gff3")
HEAD = ["bakta", "--db", str(DB), "--output", str(OUT), "--prefix", "P", "--threads", "16",
        "--keep-contig-headers"]  # fmt: skip


def _bakta(table: int, **kw: Any) -> list[str]:
    from gippa.annotate.bakta_cmd import bakta_command

    args: dict[str, Any] = {
        "genome": GENOME,
        "outdir": OUT,
        "prefix": "P",
        "db": DB,
        "threads": 16,
        "all_circular": False,
        "complete": "auto",
        "replicons_tsv": REPLICONS,
        "regions_gff3": None,
        "extra_args": [],
        "translation_table": table,
    }
    args.update(kw)
    cmd = bakta_command(**args)
    assert all(isinstance(a, str) for a in cmd)
    return cmd


def test_bakta_example8_table_25_right_before_replicons() -> None:
    cmd = _bakta(25)
    assert cmd == [*HEAD, "--translation-table", "25", "--replicons", str(REPLICONS), str(GENOME)]


def test_bakta_table_11_argv_is_todays_argv_d2() -> None:
    """D2: no flag at 11 (Bakta's own default), so the code-11 Bakta caches stay valid."""
    assert _bakta(11) == [*HEAD, "--replicons", str(REPLICONS), str(GENOME)]
    assert "--translation-table" not in _bakta(11, all_circular=True, regions_gff3=REGIONS)


@pytest.mark.parametrize("table", [4, 25])
@pytest.mark.parametrize("complete", [True, False, "auto"])
@pytest.mark.parametrize("all_circular", [True, False])
@pytest.mark.parametrize("regions", [None, REGIONS])
def test_bakta_flag_only_off_11_and_otherwise_identical(
    table: int, complete: Any, all_circular: bool, regions: Path | None
) -> None:
    kw = {"complete": complete, "all_circular": all_circular, "regions_gff3": regions,
          "extra_args": ["--skip-plot"]}  # fmt: skip
    cmd, base = _bakta(table, **kw), _bakta(11, **kw)
    assert cmd.count("--translation-table") == 1, cmd
    i = cmd.index("--translation-table")
    assert cmd[i + 1] == str(table)
    assert cmd[i + 2] == "--replicons", "example 8: the flag sits right before --replicons"
    assert cmd[:i] + cmd[i + 2 :] == base, "nothing else changes with the table"
    if "--complete" in cmd:
        assert cmd.index("--complete") < i


def test_bakta_module_passes_the_flag_only_off_11() -> None:
    """The BAKTA module builds the same flag from ``params.annotation.translation_table``."""
    script = _script("bakta")
    assert "--translation-table" in script
    assert "translation_table" in script, "the module must read params.annotation.translation_table"
    assert re.search(r"\b11\b", script), "the module must leave the flag out at table 11 (D2)"
    line = next(ln for ln in script.splitlines() if re.match(r"\s*mode=", ln))
    assert "--replicons" in line


# ================================================================ 3.2 miniprot -T


@pytest.mark.parametrize("table", TABLES)
def test_miniprot_index_takes_the_table(table: int) -> None:
    from gippa.mapping.miniprot import miniprot_index_cmd

    cmd = miniprot_index_cmd(Path("g.fa"), Path("g.mpi"), 8, table=table)
    assert cmd == ["miniprot", "-t", "8", "-T", str(table), "-d", "g.mpi", "g.fa"]


def test_miniprot_example8_table_25() -> None:
    from gippa.mapping.miniprot import miniprot_index_cmd

    cmd = miniprot_index_cmd(Path("/a/genome.norm.fa"), Path("/b/genome.mpi"), 2, table=25)
    i = cmd.index("-T")
    assert cmd[i : i + 2] == ["-T", "25"]


def test_miniprot_module_index_uses_the_table() -> None:
    script = _script("miniprot")
    index = next(ln for ln in script.splitlines() if re.search(r"(?<![\w-])-d\s", ln))
    m = re.search(r"-T\s+\$\{([^}]*)\}", index)
    assert m, f"miniprot index line must give -T ${{...}}: {index.strip()}"
    assert "translation_table" in m.group(1) and "11" in m.group(1), m.group(1)
    assert not re.search(r"-T\s+11\b", script), "-T 11 hard-coded"


# ================================================================ 3.3 tblastn -db_gencode


def _eff(args: str, table: int) -> str:
    from gippa.mapping.tblastn import effective_tblastn_args

    out = effective_tblastn_args(args, table)
    assert isinstance(out, str)
    return out


@pytest.mark.parametrize("case", _cases(), ids=lambda c: c["id"])
def test_effective_args_parity_table(case: dict[str, Any]) -> None:
    if case["expected"] is None:
        with pytest.raises(ValueError):
            _eff(case["args"], case["table"])
    else:
        assert _eff(case["args"], case["table"]) == case["expected"]


@pytest.mark.parametrize(
    ("args", "table"),
    [("-db_gencode=11", 4), ("-seg no -db_gencode 25", 4), ("-db_gencode 4", 11)],
)
def test_effective_args_mismatch_raises(args: str, table: int) -> None:
    with pytest.raises(ValueError):
        _eff(args, table)


_TOKENS = st.lists(
    st.sampled_from(
        ["-seg", "no", "yes", "-evalue", "1e-5", "-max_target_seqs", "5", "-word_size"]
    ),
    max_size=8,
)


@settings(max_examples=200, deadline=None)
@given(tokens=_TOKENS, table=st.sampled_from(TABLES))
def test_property_append_once_and_idempotent(tokens: list[str], table: int) -> None:
    """Without ``-db_gencode`` the table is appended after the user's tokens; applying the
    function again changes nothing (no duplicate: BLAST+ rejects a repeated argument)."""
    args = " ".join(tokens)
    out = _eff(args, table)
    assert shlex.split(out) == [*tokens, "-db_gencode", str(table)]
    assert out == out.strip()
    assert _eff(out, table) == out


@settings(max_examples=200, deadline=None)
@given(tokens=_TOKENS, table=st.sampled_from(TABLES), at=st.integers(0, 8), data=st.data())
def test_property_own_value_unchanged_other_value_raises(
    tokens: list[str], table: int, at: int, data: st.DataObject
) -> None:
    at = min(at, len(tokens))
    same = " ".join([*tokens[:at], "-db_gencode", str(table), *tokens[at:]])
    assert _eff(same, table) == same
    other = data.draw(st.sampled_from([t for t in TABLES if t != table]))
    with pytest.raises(ValueError):
        _eff(" ".join([*tokens[:at], "-db_gencode", str(other), *tokens[at:]]), table)


@pytest.mark.parametrize("table", TABLES)
def test_tblastn_cmd_uses_the_effective_args(table: int) -> None:
    from gippa.mapping.tblastn import OUTFMT, tblastn_cmd

    cmd = tblastn_cmd(Path("db"), Path("p.faa"), 2, NEW_DEFAULT, table=table)
    assert cmd == [
        "tblastn", "-db", "db", "-query", "p.faa", "-num_threads", "2",
        "-seg", "no", "-evalue", "1e-5", "-db_gencode", str(table), "-outfmt", OUTFMT,
    ]  # fmt: skip


def test_tblastn_cmd_mismatch_raises() -> None:
    from gippa.mapping.tblastn import tblastn_cmd

    with pytest.raises(ValueError):
        tblastn_cmd(Path("db"), Path("p.faa"), 2, "-db_gencode 11 -seg no", table=4)


def _tblastn_block(conf: str) -> str:
    """The body of the ``withName: 'TBLASTN'`` selector, braces matched (closures nest)."""
    m = re.search(r"withName\s*:\s*['\"]TBLASTN['\"]\s*\{", conf)
    assert m, "conf/base.config: no withName: 'TBLASTN' selector"
    depth, i = 1, m.end()
    while depth:
        depth += {"{": 1, "}": -1}.get(conf[i], 0)
        i += 1
    return conf[m.end() : i - 1]


def test_base_config_closure_reads_args_and_table() -> None:
    body = _tblastn_block(_strip_comments(BASE_CONFIG.read_text(encoding="utf-8")))
    assert "ext.args" in body
    assert "tblastn_args" in body and "translation_table" in body, body
    assert "db_gencode" in body, body


# ---------------------------------------------------------------- D1: the shipped default


def test_shipped_tblastn_default_has_no_db_gencode() -> None:
    from gippa.config import resolve_config
    from gippa.config.model import ProteinsConfig

    assert ProteinsConfig().tblastn_args == NEW_DEFAULT
    assert resolve_config(overrides=[]).mapping.proteins.tblastn_args == NEW_DEFAULT
    schema = json.loads(CONFIG_SCHEMA.read_text(encoding="utf-8"))
    found = [d["properties"]["tblastn_args"]["default"] for d in schema.get("$defs", {}).values()
             if "tblastn_args" in d.get("properties", {})]  # fmt: skip
    assert found == [NEW_DEFAULT], found


def test_defaults_yaml_value_and_comment() -> None:
    text = DEFAULTS_YAML.read_text(encoding="utf-8")
    assert yaml.safe_load(text)["mapping"]["proteins"]["tblastn_args"] == NEW_DEFAULT
    assert "-T 11 is applied at indexing time" not in text
    assert "# -T <annotation.translation_table> is applied at indexing time" in text


def test_bridge_removed_old_default_string_at_table_4_is_an_error(tmp_path: Path) -> None:
    """16a exempted exactly the old shipped default; that exemption goes with D1."""
    report, _ = _run(
        tmp_path,
        SMALL,
        "annotation.translation_table=4",
        f"mapping.proteins.tblastn_args={OLD_DEFAULT}",
    )
    check = _only(report, TOOL_RULE)
    assert check.level == "error"
    assert check.value == ["mapping.proteins.tblastn_args"]
    assert report.exit_code == int(ExitCode.INVALID_INPUT)


# ================================================================ 3.4 precomputed Bakta JSON

CONTIG = "ATG" + "GCT" * 18 + "TAA"  # 60 nt


def _json(tmp: Path, table: int | None) -> tuple[Path, Path]:
    tmp.mkdir(parents=True, exist_ok=True)
    data: dict[str, Any] = {
        "sequences": [{"id": "c1", "type": "chromosome", "topology": "linear", "length": 60,
                       "nt": CONTIG}],
        "features": [{"type": "cds", "sequence": "c1", "start": 1, "stop": 60, "strand": "+",
                      "locus": "G1", "aa": "M" + "A" * 18}],
        "version": {"bakta": "1.12.2", "db": {"version": "6.0", "type": "light"}},
    }  # fmt: skip
    if table is not None:
        data["genome"] = {"genus": "Mycoplasma", "species": "hyorhinis", "translation_table": table}
    js, fa = tmp / "bakta.json", tmp / "genome.norm.fa"
    js.write_text(json.dumps(data), encoding="utf-8")
    fa.write_text(f">c1\n{CONTIG}\n", encoding="utf-8")
    return js, fa


def _mismatch_message(j: int, t: int) -> str:
    return (
        f"annotation JSON was made with translation table {j}; this run uses {t}: re-run Bakta "
        f"with --translation-table {t} or set annotation.translation_table: {j}"
    )


def _unrecorded(records: list[logging.LogRecord]) -> list[logging.LogRecord]:
    return [
        r
        for r in records
        if r.levelno == logging.INFO
        and (
            getattr(r, "event", None) == "annotation.table_unrecorded"
            or "annotation.table_unrecorded" in r.getMessage()
        )
    ]


def test_example13_json_at_11_run_at_4_raises(tmp_path: Path) -> None:
    from gippa.annotate.bakta import AnnotationTableMismatch, parse_bakta_json

    assert issubclass(AnnotationTableMismatch, ValueError)
    js, _ = _json(tmp_path, 11)
    with pytest.raises(AnnotationTableMismatch) as info:
        parse_bakta_json(js, translation_table=4)
    assert _mismatch_message(11, 4) in str(info.value)


@pytest.mark.parametrize(("j", "t"), [(4, 11), (25, 4), (4, 25)])
def test_any_recorded_mismatch_raises(tmp_path: Path, j: int, t: int) -> None:
    from gippa.annotate.bakta import AnnotationTableMismatch, parse_bakta_json

    js, _ = _json(tmp_path, j)
    with pytest.raises(AnnotationTableMismatch) as info:
        parse_bakta_json(js, translation_table=t)
    assert _mismatch_message(j, t) in str(info.value)


@pytest.mark.parametrize("table", TABLES)
def test_recorded_table_equal_parses(tmp_path: Path, table: int, caplog: Any) -> None:
    from gippa.annotate.bakta import parse_bakta_json

    js, _ = _json(tmp_path, table)
    with caplog.at_level(logging.DEBUG):
        ann = parse_bakta_json(js, translation_table=table)
    assert [g.locus_tag for g in ann.model.genes] == ["G1"]
    assert not _unrecorded(caplog.records)


@pytest.mark.parametrize("table", [4, 11])
def test_unrecorded_table_parses_and_logs(tmp_path: Path, table: int, caplog: Any) -> None:
    from gippa.annotate.bakta import parse_bakta_json

    js, _ = _json(tmp_path, None)
    with caplog.at_level(logging.DEBUG):
        ann = parse_bakta_json(js, translation_table=table)
    assert [g.locus_tag for g in ann.model.genes] == ["G1"]
    assert _unrecorded(caplog.records), [r.getMessage() for r in caplog.records]


OUT_FILES = {
    "genes": "genes.parquet",
    "features": "features.parquet",
    "replicons": "replicons.parquet",
    "rrna16s": "rrna16s.fa",
    "trna_counts": "trna_counts.tsv",
    "xrefs": "xrefs.parquet",
    "annotation_meta": "annotation_meta.json",
}


def _cfg(table: int) -> Any:
    from gippa.config import resolve_config

    return resolve_config(overrides=["seed=7", f"annotation.translation_table={table}"])


def _parse_step(tmp: Path, js: Path, fa: Path, table: int) -> dict[str, Path]:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    outs = {k: tmp / "work" / v for k, v in OUT_FILES.items()}
    reg.run_step("parse_bakta", _cfg(table), {"bakta_json": js, "genome": fa}, outs)
    return outs


@pytest.mark.parametrize(("j", "t"), [(11, 11), (None, 4), (4, 4), (25, 25), (None, 11)])
def test_annotation_meta_records_the_table(tmp_path: Path, j: int | None, t: int) -> None:
    js, fa = _json(tmp_path / "in", j)
    outs = _parse_step(tmp_path, js, fa, t)
    meta = json.loads(outs["annotation_meta"].read_text(encoding="utf-8"))
    assert meta["translation_table"] == t
    assert isinstance(meta["translation_table"], int)


def test_parse_bakta_step_refuses_a_mismatch(tmp_path: Path) -> None:
    from gippa.annotate.bakta import AnnotationTableMismatch
    from gippa.steps.registry import StepError

    js, fa = _json(tmp_path / "in", 11)
    with pytest.raises(StepError) as info:
        _parse_step(tmp_path, js, fa, 4)
    assert isinstance(info.value.__cause__, AnnotationTableMismatch)


def test_example13_step_cli_exits_invalid_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: Any
) -> None:
    from gippa.cli import main

    js, fa = _json(tmp_path / "in", 11)
    resolved = tmp_path / "resolved.yaml"
    resolved.write_text(yaml.safe_dump(_cfg(4).model_dump(mode="json"), sort_keys=True))
    monkeypatch.chdir(tmp_path)
    argv = ["_step", "parse_bakta", "--config", str(resolved),
            "--in", f"bakta_json={js}", "--in", f"genome={fa}"]  # fmt: skip
    for k, v in OUT_FILES.items():
        argv += ["--out", f"{k}={tmp_path / 'work' / v}"]
    with pytest.raises(SystemExit) as info:
        main(argv)
    assert info.value.code == int(ExitCode.INVALID_INPUT)
    err = " ".join(capsys.readouterr().err.split())
    assert _mismatch_message(11, 4) in err


# ================================================================ 7 manifest


def _meta(table: int | None) -> dict[str, Any]:
    meta: dict[str, Any] = {"db_type": "light", "flags": [], "bakta_version": "1.12.2",
                            "db_version": "6.0", "annotation_source": "bakta_light",
                            "anti_sd": "CCTCCT", "n_genes": 1, "n_features": 0,
                            "skipped_types": {}}  # fmt: skip
    if table is not None:
        meta["translation_table"] = table
    return meta


@pytest.mark.parametrize("table", TABLES)
def test_manifest_annotation_translation_table_from_meta(table: int) -> None:
    m = md.build(cfg=_cfg(table), annotation_meta=_meta(table))
    assert m["annotation"]["translation_table"] == table
    assert type(m["annotation"]["translation_table"]) is int
    md.validate(m)


@pytest.mark.parametrize("table", TABLES)
def test_manifest_annotation_translation_table_falls_back_to_config(table: int) -> None:
    m = md.build(cfg=_cfg(table), annotation_meta=_meta(None))
    assert m["annotation"]["translation_table"] == table
    md.validate(m)


def test_manifest_schema_declares_the_integer_property() -> None:
    ann = md.schema()["properties"]["annotation"]
    assert ann["properties"]["translation_table"] == {"type": "integer", "enum": [4, 11, 25]}


@pytest.mark.parametrize("bad", [12, 1, "4", 4.5, None])
def test_manifest_schema_rejects_a_bad_table(bad: Any) -> None:
    m = md.build(cfg=_cfg(4), annotation_meta=_meta(4))
    assert m["annotation"]["translation_table"] == 4
    md.validate(m)
    m["annotation"]["translation_table"] = bad
    with pytest.raises(jsonschema.ValidationError):
        md.validate(m)


def test_output_schema_version_minor_bump() -> None:
    from gippa.schemas import OUTPUT_SCHEMA_VERSION

    assert OUTPUT_SCHEMA_VERSION == "1.2"  # 1.1: M11-16; 1.2: flag initiator_m, M11-20
    assert md.build()["output_schema_version"] == "1.2"


# ================================================================ 7 report


ALERT_RULE = "genome.genetic_code"


def _report_out(tmp: Path, table: int, warn: str | None = None) -> Path:
    out = run_out(tmp)
    doc = md.build(cfg=_cfg(table), annotation_meta=_meta(table))
    write_manifest(out, doc)
    if warn is not None:
        path = out / "preflight" / "preflight.json"
        body = json.loads(path.read_text(encoding="utf-8"))
        row = {"rule_id": ALERT_RULE, "input": "genome", "level": "warn", "rule": "rule",
               "value": None, "message": warn, "exit_code": 0}  # fmt: skip
        body["checks"] = [row]
        path.write_text(json.dumps(body, sort_keys=True, indent=1), encoding="utf-8")
    return out


def _methods(t: int) -> str:
    return (
        f"Genes were called and all six frames translated with NCBI genetic code {t} ({NAMES[t]})."
    )


def _note(t: int) -> str:
    return (
        f'<p class="note" data-genetic-code="{t}">Genetic code {t} ({NAMES[t]}): '
        f"TGA is read as {TGA[t]}.</p>"
    )


@pytest.mark.parametrize("table", TABLES)
def test_report_methods_sentence_every_run(tmp_path: Path, table: int) -> None:
    out = _report_out(tmp_path, table)
    p = page(out, f"annotation.translation_table={table}")
    text = section_text(p, "methods")
    assert text.count(_methods(table)) == 1, text[:400]


@pytest.mark.parametrize("table", [4, 25])
def test_report_note_opens_the_summary_off_11(tmp_path: Path, table: int) -> None:
    out = _report_out(tmp_path, table)
    text = html(out, f"annotation.translation_table={table}")
    assert text.count(_note(table)) == 1
    start = text.index('<section id="summary">')
    head = text.index("</h2>", start)
    assert text[head + len("</h2>") :].lstrip().startswith(_note(table))
    (node,) = page(out, f"annotation.translation_table={table}").root.by_attr("data-genetic-code")
    assert node.attrs["data-genetic-code"] == str(table)


def test_report_note_comes_after_the_genetic_code_alert(tmp_path: Path) -> None:
    out = _report_out(tmp_path, 4, warn="This genome probably uses genetic code 11 (TGA = stop).")
    text = html(out, "annotation.translation_table=4")
    start = text.index('<section id="summary">')
    end = text.index("</section>", start)
    alert = text.index(f'data-warning="{ALERT_RULE}"', start)
    note = text.index(_note(4), start)
    assert start < alert < note < end


def test_report_no_note_at_11(tmp_path: Path) -> None:
    out = _report_out(tmp_path, 11)
    text = html(out, "annotation.translation_table=11")
    assert "data-genetic-code" not in text
    assert _methods(11) in " ".join(text.split())
