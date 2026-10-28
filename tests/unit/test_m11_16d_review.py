"""M11-16d review round (PR #479): tblastn ``-db_gencode`` parity, quoting, the BAKTA module flag
and the type of a recorded Bakta ``translation_table`` (docs/specs/M11-16-translation-tables-4-25.md
items 2.2, 3.1, 3.3 and 3.4).

Lead decisions pinned here (2026-10-28, #479 review):

- N3: ``effective_tblastn_args`` and the TBLASTN ``ext.args`` closure
  (``workflow/conf/base.config``) agree on a repeated equal ``-db_gencode`` (refused), a bare
  trailing ``-db_gencode`` (refused with a ``ValueError``, never an ``IndexError``) and an equal
  ``-db_gencode`` behind leading spaces (returned unchanged, untrimmed). The cases live in
  ``workflow/tests/data/tblastn/db_gencode_cases.json`` (read by
  ``test_m11_16d_tools.test_effective_args_parity_table``) and one nf-test per case in
  ``workflow/tests/tblastn.nf.test``; ``test_parity_cases_match_the_nf_tests_one_to_one`` keeps
  the two in sync.
- N1: a quoted ``-db_gencode`` value (``"4"``, ``'4'``, ``-db_gencode="4"``) is refused everywhere:
  ``effective_tblastn_args`` raises ``ValueError``, ``config.tool_genetic_code`` reports ``error``
  and the closure throws (it compares raw whitespace tokens). A shell would pass ``4`` to BLAST+,
  but the closure cannot unquote, so the Python side follows it.
- N8: ``parse_bakta_json(..., translation_table=t)`` raises ``BaktaParseError`` naming the field
  when the recorded ``genome.translation_table`` is not an int (a string, a float, a bool, a
  list); it is not an ``AnnotationTableMismatch``.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from unit.test_m11_16d_tools import CASES, TABLES, _cases, _json, _script

from gippa.config import GippaConfig, resolve_config
from gippa.exitcodes import ExitCode

REPO = Path(__file__).resolve().parents[2]
NF_TEST = REPO / "workflow" / "tests" / "tblastn.nf.test"
TBLASTN = "mapping.proteins.tblastn_args"


def _eff(args: str, table: int) -> str:
    from gippa.mapping.tblastn import effective_tblastn_args

    return effective_tblastn_args(args, table)


def _cfg(args: str, table: int) -> GippaConfig:
    cfg = resolve_config(overrides=[f"annotation.translation_table={table}"])
    proteins = cfg.mapping.proteins.model_copy(update={"tblastn_args": args})
    mapping = cfg.mapping.model_copy(update={"proteins": proteins})
    cfg = cfg.model_copy(update={"mapping": mapping})
    assert cfg.mapping.proteins.tblastn_args == args
    return cfg


def _tool_check(args: str, table: int) -> Any:
    from gippa.preflight.tool_genetic_code import tool_genetic_code_check

    return tool_genetic_code_check(_cfg(args, table))


# ================================================================ N3 parity cases


def _case(case_id: str) -> dict[str, Any]:
    (case,) = [c for c in _cases() if c["id"] == case_id]
    return case


def test_repeated_equal_db_gencode_raises() -> None:
    """BLAST+ rejects a repeated argument even when both values agree (kills M1)."""
    case = _case("repeat-same")
    assert case == {"id": "repeat-same", "args": "-db_gencode 4 -db_gencode 4", "table": 4,
                    "expected": None}  # fmt: skip
    with pytest.raises(ValueError, match="-db_gencode"):
        _eff("-db_gencode 4 -db_gencode 4", 4)
    for t in TABLES:
        with pytest.raises(ValueError):
            _eff(f"-db_gencode={t} -seg no -db_gencode {t}", t)


def test_bare_trailing_db_gencode_raises_value_error_not_index_error() -> None:
    """A trailing ``-db_gencode`` without a value is refused as a ``ValueError`` (kills M3, M7)."""
    case = _case("bare-trailing")
    assert case["args"] == "-seg no -db_gencode" and case["expected"] is None
    for t in TABLES:
        try:
            _eff("-seg no -db_gencode", t)
        except ValueError:
            pass
        except Exception as exc:  # an IndexError is a crash, not a refusal
            pytest.fail(f"bare trailing -db_gencode at {t}: {type(exc).__name__}: {exc}")
        else:
            pytest.fail(f"bare trailing -db_gencode at {t} was accepted")


def test_equal_db_gencode_behind_leading_spaces_is_returned_untrimmed() -> None:
    """One equal ``-db_gencode`` passes through byte for byte, as the closure returns ``args``
    (kills M5)."""
    case = _case("leading-spaces-same")
    assert case["expected"] == case["args"] == "  -db_gencode 4 -seg no"
    assert _eff("  -db_gencode 4 -seg no", 4) == "  -db_gencode 4 -seg no"
    assert _eff("-db_gencode 25 -seg no  ", 25) == "-db_gencode 25 -seg no  "


@pytest.mark.parametrize(
    "args", ["-db_gencode 4 -db_gencode 4", "-seg no -db_gencode", "-db_gencode=4 -db_gencode 4"]
)
def test_preflight_refuses_the_new_raise_cases(args: str) -> None:
    check = _tool_check(args, 4)
    assert check.level == "error"
    assert check.value == [TBLASTN]
    assert check.exit_code == int(ExitCode.INVALID_INPUT)


# ================================================================ N1 quoting

QUOTED = [
    pytest.param('-db_gencode "4" -seg no', 4, id="double"),
    pytest.param("-db_gencode '4' -seg no", 4, id="single"),
    pytest.param('-db_gencode="25" -seg no', 25, id="equals-double"),
    pytest.param("-seg no -db_gencode '11'", 11, id="single-trailing"),
    pytest.param("-db_gencode='25'", 25, id="equals-single-only"),
]


@pytest.mark.parametrize(("args", "table"), QUOTED)
def test_quoted_db_gencode_value_is_refused_by_effective_args(args: str, table: int) -> None:
    """The closure compares raw tokens (``"4"`` is not ``4``) and throws; Python refuses the same
    string (kills M6's twin: ``shlex.split`` unquotes and so accepts it today)."""
    with pytest.raises(ValueError):
        _eff(args, table)


@pytest.mark.parametrize(("args", "table"), QUOTED)
def test_quoted_db_gencode_value_is_refused_by_tblastn_cmd(args: str, table: int) -> None:
    from gippa.mapping.tblastn import tblastn_cmd

    with pytest.raises(ValueError):
        tblastn_cmd(Path("db"), Path("p.faa"), 2, args, table=table)


@pytest.mark.parametrize(("args", "table"), QUOTED)
def test_quoted_db_gencode_value_is_a_preflight_error(args: str, table: int) -> None:
    check = _tool_check(args, table)
    assert check.level == "error", check
    assert check.value == [TBLASTN]
    assert check.exit_code == int(ExitCode.INVALID_INPUT)


def test_quoted_cases_are_in_the_parity_table() -> None:
    quoted = {c["id"]: c for c in _cases() if c["id"].startswith("quoted-")}
    assert set(quoted) == {"quoted-double", "quoted-single", "quoted-equals"}
    assert all(c["expected"] is None for c in quoted.values())
    assert all(("'" in c["args"]) or ('"' in c["args"]) for c in quoted.values())


@pytest.mark.parametrize("table", TABLES)
def test_quoting_elsewhere_in_args_is_not_refused(table: int) -> None:
    """Only a quoted ``-db_gencode`` value is refused; quoting another option is still fine."""
    args = "-seg 'no' -evalue \"1e-5\""
    assert _tool_check(args, table).level == "ok"
    assert _eff(args, table) == f"{args} -db_gencode {table}"


_PLAIN = st.sampled_from(["-seg", "no", "yes", "-evalue", "1e-5", "-max_target_seqs", "5"])


@settings(max_examples=150, deadline=None)
@given(
    before=st.lists(_PLAIN, max_size=4),
    after=st.lists(_PLAIN, max_size=4),
    table=st.sampled_from(TABLES),
    value=st.sampled_from(TABLES),
    quote=st.sampled_from(['"', "'"]),
    equals=st.booleans(),
)
def test_property_a_quoted_value_is_always_refused(
    before: list[str],
    after: list[str],
    table: int,
    value: int,
    quote: str,
    equals: bool,
) -> None:
    """Whatever surrounds it and whatever the value, a quoted ``-db_gencode`` value is refused by
    both Python sides, as by the closure."""
    opt = f"-db_gencode={quote}{value}{quote}" if equals else f"-db_gencode {quote}{value}{quote}"
    args = " ".join([*before, opt, *after])
    with pytest.raises(ValueError):
        _eff(args, table)
    assert _tool_check(args, table).level == "error"


_TOKENS = st.lists(
    st.sampled_from(
        ["-seg", "no", "-evalue", "1e-5", "-db_gencode", "4", "11", "25", "-db_gencode=4",
         "-db_gencode=11", '"4"', "'25'", '-db_gencode="4"']
    ),
    max_size=7,
)  # fmt: skip


@settings(max_examples=300, deadline=None)
@given(tokens=_TOKENS, table=st.sampled_from(TABLES))
def test_property_effective_args_raises_iff_preflight_refuses(
    tokens: list[str], table: int
) -> None:
    """Preflight refuses first: ``effective_tblastn_args`` raises exactly when
    ``config.tool_genetic_code`` reports the tblastn key (one rule, two places)."""
    args = " ".join(tokens)
    try:
        _eff(args, table)
        raised = False
    except ValueError:
        raised = True
    check = _tool_check(args, table)
    assert raised == (TBLASTN in check.value), (args, table, check)


# ================================================================ sync: JSON <-> nf-test

_NF_TEST = re.compile(r'^\s*test\("M11-16d db_gencode parity ([^:]+): ', re.M)


def _nf_blocks() -> dict[str, str]:
    text = NF_TEST.read_text(encoding="utf-8")
    starts = [(m.group(1), m.start()) for m in _NF_TEST.finditer(text)]
    ids = [i for i, _ in starts]
    assert len(ids) == len(set(ids)), f"duplicate nf-test parity names: {ids}"
    ends = [s for _, s in starts[1:]] + [len(text)]
    return {i: text[s:e] for (i, s), e in zip(starts, ends, strict=True)}


def _groovy(s: str) -> str:
    return re.sub(r"\\(.)", r"\1", s)


def test_parity_cases_match_the_nf_tests_one_to_one() -> None:
    cases = _cases()
    ids = [c["id"] for c in cases]
    assert len(ids) == len(set(ids)), f"duplicate case ids in {CASES.name}: {ids}"
    blocks = _nf_blocks()
    assert sorted(blocks) == sorted(ids), {
        "only in json": sorted(set(ids) - set(blocks)),
        "only in nf-test": sorted(set(blocks) - set(ids)),
    }


@pytest.mark.parametrize("case", _cases(), ids=lambda c: c["id"])
def test_each_nf_test_runs_its_case(case: dict[str, Any]) -> None:
    """The nf-test for a case sets the same args and table and expects the same outcome."""
    block = _nf_blocks()[case["id"]]
    m = re.search(r'tblastn_args:\s*"((?:[^"\\]|\\.)*)"', block)
    assert m, block
    assert _groovy(m.group(1)) == case["args"]
    t = re.search(r"translation_table:\s*(\d+)", block)
    assert t and int(t.group(1)) == case["table"], block
    if case["expected"] is None:
        assert "assert process.failed" in block
        assert "sets -db_gencode" in block, "a failing case must check the closure's message"
    else:
        assert "assert process.success" in block
        g = re.search(r'assert got == "((?:[^"\\]|\\.)*)"', block)
        assert g, block
        assert _groovy(g.group(1)) == case["expected"]


def test_7c_nf_test_checks_the_failure_message() -> None:
    block = _nf_blocks()["7c"]
    assert "process.failed" in block
    assert re.search(r"process\.(stdout|stderr)", block), "7c must read the run's output"
    assert "this run uses annotation.translation_table 4" in block


# ================================================================ 3.1 BAKTA module flag


def test_bakta_module_flag_sits_in_a_conditional_on_the_table() -> None:
    """``--translation-table`` is built only on the not-11 branch of a test of the table, and the
    script uses it once, at the start of ``mode=`` (before ``--replicons``). An unconditional flag
    fails here."""
    script = _script("bakta")
    assert script.count("--translation-table") == 1, script
    (line,) = [ln for ln in script.splitlines() if "--translation-table" in ln]
    tbl = re.search(r"def\s+(\w+)\s*=\s*params\.annotation\?\.translation_table\s*\?:\s*11", script)
    assert tbl, "the module must read params.annotation.translation_table (default 11)"
    t = re.escape(tbl.group(1))
    flag = rf"\"--translation-table \$\{{{t}\}} ?\""
    eq11 = rf"def\s+(\w+)\s*=\s*{t}\s*==\s*11\s*\?\s*''\s*:\s*{flag}"
    ne11 = rf"def\s+(\w+)\s*=\s*{t}\s*!=\s*11\s*\?\s*{flag}\s*:\s*''"
    m = re.search(eq11, line) or re.search(ne11, line)
    assert m, f"the flag must be '' at 11 and --translation-table <t> otherwise: {line.strip()}"
    mode = next(ln for ln in script.splitlines() if re.match(r"\s*mode=", ln))
    assert re.match(rf'\s*mode="\$\{{{re.escape(m.group(1))}\}}--replicons ', mode), mode


# ================================================================ N8 Bakta table type


def _json_with(tmp: Path, recorded: Any) -> Path:
    js, _ = _json(tmp, 11)
    data = json.loads(js.read_text(encoding="utf-8"))
    data["genome"]["translation_table"] = recorded
    js.write_text(json.dumps(data), encoding="utf-8")
    return js


@pytest.mark.parametrize("table", TABLES)
@pytest.mark.parametrize(
    "recorded",
    [
        pytest.param("11", id="str-11"),
        pytest.param("4", id="str-4"),
        pytest.param(11.0, id="float-11"),
        pytest.param(4.0, id="float-4"),
        pytest.param(25.5, id="float-25.5"),
        pytest.param(True, id="bool-true"),
        pytest.param(False, id="bool-false"),
        pytest.param([11], id="list"),
        pytest.param({"id": 11}, id="dict"),
    ],
)
def test_non_int_recorded_table_is_a_parse_error(tmp_path: Path, recorded: Any, table: int) -> None:
    from gippa.annotate.bakta import AnnotationTableMismatch, BaktaParseError, parse_bakta_json

    js = _json_with(tmp_path, recorded)
    with pytest.raises(BaktaParseError) as info:
        parse_bakta_json(js, translation_table=table)
    assert not isinstance(info.value, AnnotationTableMismatch)
    assert "translation_table" in str(info.value), str(info.value)


@pytest.mark.parametrize("table", TABLES)
def test_int_recorded_table_still_compares(tmp_path: Path, table: int) -> None:
    from gippa.annotate.bakta import AnnotationTableMismatch, parse_bakta_json

    other = next(t for t in TABLES if t != table)
    assert parse_bakta_json(_json_with(tmp_path / "eq", table), translation_table=table)
    with pytest.raises(AnnotationTableMismatch):
        parse_bakta_json(_json_with(tmp_path / "ne", other), translation_table=table)
