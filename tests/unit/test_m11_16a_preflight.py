"""M11-16a preflight: AC 2-4, examples 7c, 8 and 9, the tie rows, and the schema description
(docs/specs/M11-16-translation-tables-4-25.md, Definition items 2.1-2.3).

Interface assumed (spec):
- ``gippa.preflight.run._translation_table_check(table) -> Check`` (same rule id as #413);
- a new check ``config.tool_genetic_code`` (input ``inputs``) in every ``run_preflight`` report;
- ``gippa.preflight.genetic_code._decide(bp11, bp4, n4, table=t)`` (example 9) returning a level
  or a ``Check``; ``genetic_code_check(seqs, table)`` looks ``long_orf_cover`` up at call time
  (as in the M11-14w tests), so the wording tests patch it.

Example 9 names M. hyorhinis (ratio 0.227) and E. coli K-12 (ratio 0.88). Those genomes are not
available offline, so their rows use patched coverages with the same ratios; the synthetic
fixtures from ``helpers.genetic_code_fixture`` cover the scan itself. A "code-25 genome" is the
generator's TGA-sense genome (``code=4``): the scan reads nucleotides only, and code 4 and code 25
genomes are the same to it (spec item 2.3; the generator gains an explicit code 25 in 16c).

End-to-end runs at tables 4 and 25 set ``mapping.proteins.tblastn_args`` without
``-db_gencode``: the shipped default still carries ``-db_gencode 11`` until 16d (spec 3.3, D1).
Item 2.2 exempts exactly that default until 16d (the bridge); any other ``-db_gencode`` that
differs from the table is refused.
"""

from __future__ import annotations

import json
import re
from collections.abc import Sequence
from functools import cache
from pathlib import Path
from typing import Any

import jsonschema
import pytest
from helpers.genetic_code_fixture import synthetic_genome
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.config import GippaConfig, resolve_config
from gippa.exitcodes import ExitCode
from gippa.preflight.checks import Check

ROOT = Path(__file__).resolve().parents[2]
CONFIG_SCHEMA = ROOT / "src" / "gippa" / "schemas" / "config.schema.json"
PREFLIGHT_SCHEMA = ROOT / "src" / "gippa" / "schemas" / "preflight.schema.json"
LOW = ["preflight.min_n50_warn=1", "preflight.min_largest_contig_warn=1"]
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"
NO_GENCODE = "mapping.proteins.tblastn_args=-seg no -evalue 1e-5"

TT_RULE = "config.translation_table"
TT_RULE_TEXT = "annotation.translation_table in {4, 11, 25} → ok; else error"
TOOL_RULE = "config.tool_genetic_code"
GC_RULE = "genome.genetic_code"
GC_RULE_TEXT = (
    "long ORFs (>= 300 codons, six frames): bp with TGA as stop / bp with TGA as sense; "
    "table 11: < 0.55 -> warn; tables 4, 25: >= 0.55 -> warn (needs >= 20 such ORFs)"
)
MSG_UNSUPPORTED = "not evaluated: unsupported annotation.translation_table"
MSG_FEW = "not evaluated: fewer than 20 open reading frames of at least 300 codons"
VALUE_KEYS = {
    "min_codons", "bp_tga_stop", "bp_tga_sense", "n_orfs_tga_sense", "ratio", "threshold",
    "translation_table",
}  # fmt: skip
NAMES = {4: "Mycoplasma and Spiroplasma", 25: "Candidate Division SR1 and Gracilibacteria"}
BAKTA, MINIPROT, TBLASTN = (
    "annotation.bakta.extra_args",
    "mapping.proteins.miniprot_args",
    "mapping.proteins.tblastn_args",
)
SEEDS = [1, 2, 3]
TGA_SENSE = [(0.28, 0.01), (0.40, 0.005)]  # (gc, TGA rate): code-4/25-like genomes
CODE11 = [0.28, 0.50, 0.70]


def w11(q4: int, q11: int, ratio: str) -> str:
    return (
        "This genome probably uses a genetic code in which TGA is not a stop: code 4 (TGA = Trp; "
        "Mycoplasma, Spiroplasma, Mesoplasma, Ureaplasma and other Mollicutes) or code 25 "
        "(TGA = Gly; candidate division SR1, Gracilibacteria). Long open reading frames cover "
        f"{q4} % of the genome when TGA is read as a codon but only {q11} % when it is a stop "
        f"(ratio {ratio} < 0.55). This run uses annotation.translation_table: 11, so genes will "
        "be cut at TGA, peptides with W or G at TGA codons will not map, and strain_difference, "
        "readthrough and nterm_extension_evidence flags will be artefacts. Set "
        "annotation.translation_table to 4 or 25 from the organism's taxonomy; the genome "
        "sequence alone cannot tell the two apart. This warning does not change any result."
    )


def w4(table: int, q4: int, q11: int, ratio: str) -> str:
    return (
        "This genome probably uses genetic code 11 (TGA = stop). Long open reading frames cover "
        f"{q11} % of the genome when TGA is a stop and {q4} % when it is read as a codon (ratio "
        f"{ratio} >= 0.55; genomes that use code 4 or 25 are below 0.55). This run uses "
        f"annotation.translation_table: {table} ({NAMES[table]}), so genes that end in TGA will "
        "be read through into the next gene and peptides may be placed in the wrong frame. Most "
        "bacteria and archaea use code 11; check the organism's genetic code. This warning does "
        "not change any result."
    )


def norm(text: str) -> str:
    return " ".join(text.split())


def level_of(result: Any) -> str:
    return result if isinstance(result, str) else str(result.level)


def gcc(seqs: Sequence[str], table: int) -> Check:
    from gippa.preflight.genetic_code import genetic_code_check

    c = genetic_code_check(seqs, table)
    assert isinstance(c, Check)
    return c


@cache
def tga_sense(seed: int = 7, gc: float = 0.28, rate: float = 0.01) -> str:
    return synthetic_genome(seed, gc, 4, tga_rate=rate)


@cache
def code11(seed: int = 7, gc: float = 0.50) -> str:
    return synthetic_genome(seed, gc, 11)


def assert_gc_shape(c: Check, table: int) -> dict[str, Any]:
    assert c.rule_id == GC_RULE
    assert c.input == "genome"
    assert c.level in {"ok", "info", "warn"}
    assert c.exit_code == 0
    assert c.rule == GC_RULE_TEXT
    v = c.value
    assert isinstance(v, dict)
    assert set(v) == VALUE_KEYS
    assert v["translation_table"] == table
    assert v["min_codons"] == 300 and v["threshold"] == 0.55
    return v


def fake_cover(bp11: int, bp4: int, n: int = 20) -> Any:
    def fake(seq: str, stops: Any, min_codons: int) -> tuple[int, int]:
        assert min_codons == 300
        return (bp11, n) if "TGA" in set(stops) else (bp4, n)

    return fake


GENOME_2KB = ["ACGT" * 500]  # with patched coverages: q = round(100 * bp / 2000)


# ---------------------------------------------------------------- 2.1 config.translation_table


@pytest.mark.parametrize("table", [4, 11, 25])
def test_ac2_supported_tables_are_ok(table: int) -> None:
    from gippa.preflight.run import _translation_table_check

    c = _translation_table_check(table)
    assert (c.rule_id, c.input, c.level, c.message, c.exit_code) == (
        TT_RULE, "inputs", "ok", "", 0,
    )  # fmt: skip
    assert c.rule == TT_RULE_TEXT
    assert c.value == table


@pytest.mark.parametrize("table", [1, 2, 12, 0, -1, 26])
def test_ac2_other_tables_are_an_input_error(table: int) -> None:
    from gippa.preflight.run import _translation_table_check

    c = _translation_table_check(table)
    assert (c.rule_id, c.level, c.exit_code) == (TT_RULE, "error", int(ExitCode.INVALID_INPUT))
    assert c.rule == TT_RULE_TEXT
    assert c.value == table
    assert c.message == (
        f"annotation.translation_table={table}: GIPPA supports genetic codes 4, 11 and 25 "
        "(the tables Bakta accepts)"
    )


def test_string_table_is_coerced_by_pydantic() -> None:
    """Edge cases: ``--set annotation.translation_table="4"`` gives 4."""
    cfg = resolve_config(overrides=['annotation.translation_table="4"'])
    assert cfg.annotation.translation_table == 4


def test_schema_lists_the_tables_in_the_description_only() -> None:
    """Item 2.1 / D7: the values appear in the description, not as an ``enum`` constraint, so a
    bad value reaches the preflight table."""
    from gippa.config import json_schema

    for schema in (json.loads(CONFIG_SCHEMA.read_text(encoding="utf-8")), json_schema()):
        prop = schema["$defs"]["AnnotationConfig"]["properties"]["translation_table"]
        assert prop["type"] == "integer" and prop["default"] == 11
        assert "enum" not in prop
        desc = prop.get("description", "")
        for t in ("4", "11", "25"):
            assert re.search(rf"\b{t}\b", desc), desc
    schema = json.loads(CONFIG_SCHEMA.read_text(encoding="utf-8"))
    bad = resolve_config(overrides=["annotation.translation_table=26"])
    dump = bad.model_dump(mode="json")
    jsonschema.validate(dump, schema)


# ---------------------------------------------------------------- preflight runs


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    for name in ("db", "repo"):
        (tmp_path / name).mkdir(exist_ok=True)
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db"))
    monkeypatch.setenv("GIPPA_REPO", str(tmp_path / "repo"))
    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "repo"))  # never the real $GIPPA_ROOT
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


def _run(tmp_path: Path, seqs: Sequence[str], *extra: str) -> Any:
    from gippa.preflight.run import PreflightInputs, run_preflight

    g = tmp_path / "g.fa"
    g.write_text("".join(f">c{i + 1}\n{s}\n" for i, s in enumerate(seqs)))
    p = tmp_path / "proteins.faa"
    p.write_text(f">P1\n{PROTEIN}\n")
    pep = tmp_path / "peptides.tsv"
    pep.write_text("peptide\tprotein_ids\nPEPTIDEK\tP1\n")
    inputs = PreflightInputs(genome=g, proteins=p, peptides=pep)
    return run_preflight(inputs, _cfg(*extra), _probe(), "preflight")


def _only(report: Any, rule_id: str) -> Check:
    hits = [c for c in report.checks if c.rule_id == rule_id]
    assert len(hits) == 1, [c.rule_id for c in report.checks]
    return hits[0]


SMALL = ["ACGGTCAT" * 18_750]  # 150 kb; genome.genetic_code is info


@pytest.mark.parametrize("table", [4, 25])
def test_ac2_tables_4_and_25_pass_preflight(tmp_path: Path, table: int) -> None:
    report, norm_inputs = _run(tmp_path, SMALL, f"annotation.translation_table={table}", NO_GENCODE)
    assert _only(report, TT_RULE).level == "ok"
    assert _only(report, TOOL_RULE).level == "ok"
    assert report.exit_code == int(ExitCode.SUCCESS)
    assert norm_inputs.cfg.annotation.translation_table == table


# ---------------------------------------------------------------- 2.2 config.tool_genetic_code


TOOL_CASES = [
    pytest.param(11, [], [], id="default-11"),
    pytest.param(4, [NO_GENCODE], [], id="no-gencode-4"),
    pytest.param(
        25, [f"{TBLASTN}=-db_gencode 25 -seg yes -evalue 1e-7"], [], id="ex7b-same-gencode"
    ),
    pytest.param(11, [f"{TBLASTN}=-db_gencode=11 -seg no"], [], id="gencode-eq-same"),
    pytest.param(11, [f"{BAKTA}=['--keep-contig-headers']"], [], id="bakta-other-arg"),
    pytest.param(11, [f"{MINIPROT}=-S -t 4 -N 5"], [], id="miniprot-threads-not-T"),
    pytest.param(11, [f"{BAKTA}=['--translation-table', '4']"], [BAKTA], id="ex8-bakta"),
    pytest.param(11, [f"{BAKTA}=['--translation-table', '11']"], [BAKTA], id="bakta-same-value"),
    pytest.param(4, [NO_GENCODE, f"{BAKTA}=['--translation-table=4']"], [BAKTA], id="bakta-eq"),
    pytest.param(11, [f"{MINIPROT}=-S -T 4"], [MINIPROT], id="miniprot-T-token"),
    pytest.param(11, [f"{MINIPROT}=-S -T11 -N 5"], [MINIPROT], id="miniprot-T-attached"),
    pytest.param(4, [NO_GENCODE, f"{MINIPROT}=-T 4"], [MINIPROT], id="miniprot-T-same"),
    pytest.param(4, [f"{TBLASTN}=-db_gencode 11 -seg no"], [TBLASTN], id="ex7c-tblastn"),
    pytest.param(11, [f"{TBLASTN}=-db_gencode=4"], [TBLASTN], id="tblastn-eq-differs"),
    pytest.param(
        11, [f"{TBLASTN}=-db_gencode 11 -seg no -db_gencode 11"], [TBLASTN], id="tblastn-twice"
    ),
    pytest.param(
        25,
        [
            f"{BAKTA}=['--translation-table', '25']",
            f"{MINIPROT}=-S -T25",
            f"{TBLASTN}=-db_gencode 4",
        ],
        [BAKTA, MINIPROT, TBLASTN],
        id="all-three",
    ),
]


@pytest.mark.parametrize(("table", "sets", "keys"), TOOL_CASES)
def test_ac3_tool_genetic_code(
    tmp_path: Path, table: int, sets: list[str], keys: list[str]
) -> None:
    report, _ = _run(tmp_path, SMALL, f"annotation.translation_table={table}", *sets)
    c = _only(report, TOOL_RULE)
    assert c.input == "inputs"
    assert isinstance(c.value, list)
    assert sorted(c.value) == sorted(keys) and len(c.value) == len(set(c.value))
    if not keys:
        assert (c.level, c.message, c.exit_code) == ("ok", "", 0)
        assert report.exit_code == int(ExitCode.SUCCESS)
        return
    assert (c.level, c.exit_code) == ("error", int(ExitCode.INVALID_INPUT))
    assert report.exit_code == int(ExitCode.INVALID_INPUT)
    for key in keys:
        assert key in c.message
    for key in {BAKTA, MINIPROT, TBLASTN} - set(keys):
        assert key not in c.message
    assert "annotation.translation_table" in c.message
    assert "remove" in c.message.lower()


def test_tool_genetic_code_in_preflight_json(tmp_path: Path) -> None:
    from gippa.preflight.run import write_outputs

    report, norm_inputs = _run(tmp_path, SMALL, f"{BAKTA}=['--translation-table', '4']")
    out = tmp_path / "OUT"
    out.mkdir()
    write_outputs(out, report, norm_inputs)
    data = json.loads((out / "preflight" / "preflight.json").read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(json.loads(PREFLIGHT_SCHEMA.read_text("utf-8"))).validate(data)
    rows = [r for r in data["checks"] if r["rule_id"] == TOOL_RULE]
    assert len(rows) == 1
    assert rows[0]["level"] == "error" and rows[0]["value"] == [BAKTA]


# ---------------------------------------------------------------- 2.3 _decide


def test_example9_decide_rows() -> None:
    from gippa.preflight.genetic_code import _decide

    assert level_of(_decide(550, 1000, 20, table=4)) == "warn"
    assert level_of(_decide(549, 1000, 20, table=4)) == "ok"
    assert level_of(_decide(550, 1000, 20, table=11)) == "ok"
    assert level_of(_decide(549, 1000, 20, table=11)) == "warn"
    assert level_of(_decide(550, 1000, 20, table=25)) == "warn"
    assert level_of(_decide(549, 1000, 20, table=25)) == "ok"
    assert level_of(_decide(1000, 1000, 20, table=25)) == "warn"
    for table in (4, 11, 25):
        assert level_of(_decide(550, 1000, 19, table=table)) == "info"
        assert level_of(_decide(0, 0, 0, table=table)) == "info"


@settings(max_examples=300, deadline=None)
@given(bp4=st.integers(1, 10**9), frac=st.integers(0, 10**6), n4=st.integers(20, 10**5))
def test_ac4_decide_at_4_and_25_is_the_complement_of_11(bp4: int, frac: int, n4: int) -> None:
    from gippa.preflight.genetic_code import _decide

    bp11 = bp4 * frac // 10**6
    l11 = level_of(_decide(bp11, bp4, n4, table=11))
    l4 = level_of(_decide(bp11, bp4, n4, table=4))
    l25 = level_of(_decide(bp11, bp4, n4, table=25))
    assert l11 == ("warn" if 20 * bp11 < 11 * bp4 else "ok")
    assert l4 == l25 == ("warn" if 20 * bp11 >= 11 * bp4 else "ok")
    assert {l11, l4} == {"ok", "warn"}


@settings(max_examples=100, deadline=None)
@given(bp4=st.integers(0, 10**9), frac=st.integers(0, 10**6), n4=st.integers(0, 19))
def test_ac4_decide_below_20_orfs_is_info_at_every_table(bp4: int, frac: int, n4: int) -> None:
    from gippa.preflight.genetic_code import _decide

    bp11 = bp4 * frac // 10**6
    for table in (4, 11, 25):
        assert level_of(_decide(bp11, bp4, n4, table=table)) == "info"


# ---------------------------------------------------------------- 2.3 example 9 and the texts


def test_example9_mhyorhinis_at_11_is_w11(monkeypatch: pytest.MonkeyPatch) -> None:
    import gippa.preflight.genetic_code as gc_mod

    monkeypatch.setattr(gc_mod, "long_orf_cover", fake_cover(227, 1000))  # ratio 0.227
    c = gcc(GENOME_2KB, 11)
    v = assert_gc_shape(c, 11)
    assert c.level == "warn"
    assert (v["bp_tga_stop"], v["bp_tga_sense"]) == (227, 1000)
    assert norm(c.message) == w11(50, 11, "0.22")


def test_example9_mhyorhinis_at_4_is_ok(monkeypatch: pytest.MonkeyPatch) -> None:
    import gippa.preflight.genetic_code as gc_mod

    monkeypatch.setattr(gc_mod, "long_orf_cover", fake_cover(227, 1000))
    c = gcc(GENOME_2KB, 4)
    assert_gc_shape(c, 4)
    assert (c.level, c.message) == ("ok", "")


@pytest.mark.parametrize("table", [4, 25])
def test_example9_ecoli_at_4_or_25_is_w4(monkeypatch: pytest.MonkeyPatch, table: int) -> None:
    import gippa.preflight.genetic_code as gc_mod

    monkeypatch.setattr(gc_mod, "long_orf_cover", fake_cover(880, 1000))  # ratio 0.88
    c = gcc(GENOME_2KB, table)
    assert_gc_shape(c, table)
    assert c.level == "warn"
    assert norm(c.message) == w4(table, 50, 44, "0.88")
    assert f"({NAMES[table]})" in c.message


def test_example9_ecoli_at_11_is_ok(monkeypatch: pytest.MonkeyPatch) -> None:
    import gippa.preflight.genetic_code as gc_mod

    monkeypatch.setattr(gc_mod, "long_orf_cover", fake_cover(880, 1000))
    c = gcc(GENOME_2KB, 11)
    assert_gc_shape(c, 11)
    assert (c.level, c.message) == ("ok", "")


@pytest.mark.parametrize(
    ("bp11", "table", "level", "ratio"),
    [
        (550, 4, "warn", "0.55"),  # tie: ratio exactly 0.55
        (550, 25, "warn", "0.55"),
        (550, 11, "ok", None),
        (549, 11, "warn", "0.54"),  # truncated, never "0.55 < 0.55"
        (549, 4, "ok", None),
        (1000, 25, "warn", "1.00"),  # bp11 == bp4: not "0.100"
        (1000, 4, "warn", "1.00"),
        (999, 4, "warn", "0.99"),
    ],
)
def test_tie_rows_and_ratio_printing(
    monkeypatch: pytest.MonkeyPatch, bp11: int, table: int, level: str, ratio: str | None
) -> None:
    import gippa.preflight.genetic_code as gc_mod

    monkeypatch.setattr(gc_mod, "long_orf_cover", fake_cover(bp11, 1000))
    c = gcc(GENOME_2KB, table)
    assert_gc_shape(c, table)
    assert c.level == level
    if ratio is None:
        assert c.message == ""
        return
    text = norm(c.message)
    op = "<" if table == 11 else ">="
    assert f"(ratio {ratio} {op} 0.55" in text
    assert "0.100" not in text and "0.55 < 0.55" not in text
    q11 = (200 * bp11 + 2000) // 4000
    want = w11(50, q11, ratio) if table == 11 else w4(table, 50, q11, ratio)
    assert text == want


@pytest.mark.parametrize("table", [1, 2, 12, 26, 0])
def test_unsupported_table_is_info(table: int) -> None:
    for seqs in ([tga_sense()], [code11()], []):
        c = gcc(seqs, table)
        assert c.rule_id == GC_RULE and c.rule == GC_RULE_TEXT
        assert (c.level, c.exit_code) == ("info", 0)
        assert c.message == MSG_UNSUPPORTED


@pytest.mark.parametrize("table", [4, 11, 25])
def test_info_paths_at_every_table(table: int) -> None:
    c = gcc([], table)
    v = assert_gc_shape(c, table)
    assert (c.level, c.message) == ("info", MSG_FEW)
    assert v["ratio"] is None
    c = gcc(["A" * 10_000], table)
    assert_gc_shape(c, table)
    assert (c.level, c.message) == ("info", MSG_FEW)


# ---------------------------------------------------------------- fixtures at 4, 11 and 25


@pytest.mark.parametrize(("gc", "rate"), TGA_SENSE)
@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize(("table", "want"), [(11, "warn"), (4, "ok"), (25, "ok")])
def test_tga_sense_fixture(seed: int, gc: float, rate: float, table: int, want: str) -> None:
    c = gcc([tga_sense(seed, gc, rate)], table)
    v = assert_gc_shape(c, table)
    assert c.level == want, v
    assert v["n_orfs_tga_sense"] >= 20
    if want == "warn":
        assert c.message.startswith("This genome probably uses a genetic code in which TGA is")
    else:
        assert c.message == ""


@pytest.mark.parametrize("gc", CODE11)
@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize(("table", "want"), [(11, "ok"), (4, "warn"), (25, "warn")])
def test_code11_fixture(seed: int, gc: float, table: int, want: str) -> None:
    c = gcc([code11(seed, gc)], table)
    v = assert_gc_shape(c, table)
    assert c.level == want, v
    if want == "warn":
        assert c.message.startswith("This genome probably uses genetic code 11 (TGA = stop).")
        assert f"annotation.translation_table: {table} ({NAMES[table]})" in c.message
    else:
        assert c.message == ""


def test_same_genome_value_differs_only_in_the_table() -> None:
    seqs = [code11(2)]
    vals = {t: assert_gc_shape(gcc(seqs, t), t) for t in (4, 11, 25)}
    strip = [{k: x for k, x in v.items() if k != "translation_table"} for v in vals.values()]
    assert strip[0] == strip[1] == strip[2]


# ---------------------------------------------------------------- wiring at table 4


def test_preflight_at_4_with_code11_genome_warns_and_strict_exits_3(tmp_path: Path) -> None:
    from gippa.preflight.run import exit_code

    report, _ = _run(tmp_path, [code11()], "annotation.translation_table=4", NO_GENCODE)
    c = _only(report, GC_RULE)
    assert_gc_shape(c, 4)
    assert c.level == "warn"
    others = [x.rule_id for x in report.checks if x.level in {"warn", "error"} and x is not c]
    assert not others
    assert exit_code(report, False) == ExitCode.SUCCESS
    assert exit_code(report, True) == ExitCode.WARNINGS_STRICT


@pytest.mark.parametrize("table", [4, 25])
def test_preflight_at_4_or_25_with_tga_sense_genome_is_ok(tmp_path: Path, table: int) -> None:
    from gippa.preflight.run import exit_code

    extra = (f"annotation.translation_table={table}", NO_GENCODE)
    report, _ = _run(tmp_path, [tga_sense()], *extra)
    c = _only(report, GC_RULE)
    assert_gc_shape(c, table)
    assert c.level == "ok"
    assert c == gcc([tga_sense()], table)
    assert exit_code(report, True) == ExitCode.SUCCESS
