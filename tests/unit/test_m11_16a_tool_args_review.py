"""M11-16a #473 review: ``config.tool_genetic_code`` parses options as the tools and the shell do
(docs/specs/M11-16-translation-tables-4-25.md item 2.2, "Parsing follows ...").

- option strings are ``shlex.split``; an ``extra_args`` element may hold several words;
- a Bakta token whose name (before ``=``) is an argparse prefix of ``--translation-table`` of at
  least 4 characters is a conflict;
- a miniprot single-dash token with ``T`` after the dash is a conflict (clusters ``-ST4``);
- a trailing bare ``-db_gencode`` is a conflict.

Tested through ``gippa.preflight.tool_genetic_code.tool_genetic_code_check(cfg)`` and once
through ``run_preflight`` (the check is in every report).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from gippa.config import GippaConfig, resolve_config
from gippa.exitcodes import ExitCode
from gippa.preflight.checks import Check

BAKTA = "annotation.bakta.extra_args"
MINIPROT = "mapping.proteins.miniprot_args"
TBLASTN = "mapping.proteins.tblastn_args"
NO_GENCODE = "-seg no -evalue 1e-5"


def _cfg(table: int, **fields: Any) -> GippaConfig:
    """The default config at ``table`` with dotted ``fields`` set (no YAML parsing of values,
    so shell quotes reach the check verbatim)."""
    data = resolve_config(overrides=[f"annotation.translation_table={table}"]).model_dump()
    if table != 11:
        data["mapping"]["proteins"]["tblastn_args"] = NO_GENCODE
    for dotted, value in fields.items():
        node = data
        *path, last = dotted.split("__")
        for part in path:
            node = node[part]
        node[last] = value
    return GippaConfig.model_validate(data)


def _check(table: int, **fields: Any) -> Check:
    from gippa.preflight.tool_genetic_code import tool_genetic_code_check

    c = tool_genetic_code_check(_cfg(table, **fields))
    assert c.rule_id == "config.tool_genetic_code" and c.input == "inputs"
    return c


def _assert_error(c: Check, key: str) -> None:
    assert (c.level, c.exit_code) == ("error", int(ExitCode.INVALID_INPUT)), c
    assert c.value == [key]
    assert key in c.message


def _assert_ok(c: Check) -> None:
    assert (c.level, c.value, c.message, c.exit_code) == ("ok", [], "", 0), c


# ---------------------------------------------------------------- 1. Bakta argparse prefixes

BAKTA_FIELD = "annotation__bakta__extra_args"


@pytest.mark.parametrize(
    "args",
    [
        ["--transl", "4"],
        ["--tr=4"],
        ["--translation-tab", "4"],
        ["--tr", "4"],
        ["--translation-table 4"],  # one element, several words: the module joins with spaces
        ["--keep-contig-headers --transl=4"],
    ],
    ids=lambda a: repr(a),
)
def test_bakta_prefix_or_joined_element_is_an_error(args: list[str]) -> None:
    _assert_error(_check(11, **{BAKTA_FIELD: args}), BAKTA)


@pytest.mark.parametrize(
    "args",
    [
        ["--threads", "4"],
        ["--tmp-dir", "x"],
        ["--threads 4"],
        ["--keep-contig-headers"],
        ["--t", "4"],  # shorter than 4 characters: not a conflict
    ],
    ids=lambda a: repr(a),
)
def test_bakta_unrelated_flags_are_ok(args: list[str]) -> None:
    _assert_ok(_check(11, **{BAKTA_FIELD: args}))


# ---------------------------------------------------------------- 2. miniprot clusters

MINIPROT_FIELD = "mapping__proteins__miniprot_args"


@pytest.mark.parametrize("args", ["-ST4", "-uT4", "-uT 4", "-S -N 5 -uT4"])
def test_miniprot_cluster_with_T_is_an_error(args: str) -> None:
    _assert_error(_check(11, **{MINIPROT_FIELD: args}), MINIPROT)


@pytest.mark.parametrize("args", ["-t 4", "-u", "--outn=1", "-S -N 5 --outs=0.95"])
def test_miniprot_other_options_are_ok(args: str) -> None:
    _assert_ok(_check(11, **{MINIPROT_FIELD: args}))


# ---------------------------------------------------------------- 3. tblastn and the shell

TBLASTN_FIELD = "mapping__proteins__tblastn_args"


def test_tblastn_quoted_option_name_is_split_like_the_shell() -> None:
    _assert_error(_check(11, **{TBLASTN_FIELD: "'-db_gencode' 4"}), TBLASTN)


def test_tblastn_quoted_value_equal_to_the_table_is_ok() -> None:
    _assert_ok(_check(4, **{TBLASTN_FIELD: "-db_gencode '4'"}))


def test_tblastn_trailing_bare_db_gencode_is_an_error() -> None:
    _assert_error(_check(11, **{TBLASTN_FIELD: "-seg no -db_gencode"}), TBLASTN)


# ---------------------------------------------------------------- wiring


@pytest.fixture
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    for name in ("db", "repo"):
        (tmp_path / name).mkdir(exist_ok=True)
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db"))
    monkeypatch.setenv("GIPPA_REPO", str(tmp_path / "repo"))
    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "repo"))  # never the real $GIPPA_ROOT
    monkeypatch.delenv("NXF_APPTAINER_CACHEDIR", raising=False)


@pytest.mark.usefixtures("_isolated")
def test_run_preflight_refuses_a_bakta_prefix(tmp_path: Path) -> None:
    from gippa.preflight.run import PreflightInputs, run_preflight

    g = tmp_path / "g.fa"
    g.write_text(">chr1\n" + ("ACGGTCAT" * 18_750) + "\n")
    p = tmp_path / "proteins.faa"
    p.write_text(">P1\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA\n")
    cfg = _cfg(11, **{BAKTA_FIELD: ["--transl", "4"]})
    report = run_preflight(PreflightInputs(genome=g, proteins=p), cfg, None, "preflight")[0]
    hits = [c for c in report.checks if c.rule_id == "config.tool_genetic_code"]
    assert len(hits) == 1
    _assert_error(hits[0], BAKTA)
    assert report.exit_code == int(ExitCode.INVALID_INPUT)
