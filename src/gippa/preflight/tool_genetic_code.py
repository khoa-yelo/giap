"""``config.tool_genetic_code``: no tool option may set a genetic code itself (M11-16 item 2.2).

GIPPA sets every tool's genetic code from ``annotation.translation_table``, so a Bakta
``--translation-table``, a miniprot ``-T`` or a tblastn ``-db_gencode`` that disagrees with the
table (or is repeated) is an input error. Options are parsed as the shell and the tools see them
(spec 2.2, "Parsing follows ..."): every string is ``shlex.split``; a Bakta argparse prefix of
``--translation-table`` of at least 4 characters counts; a miniprot single-dash cluster with ``T``
counts. A string the shell cannot split is reported as unparseable. Pure: reads the config only.
"""

from __future__ import annotations

import shlex
from collections.abc import Callable, Sequence

from gippa.config import GippaConfig
from gippa.exitcodes import ExitCode
from gippa.mapping.tblastn import db_gencode_values
from gippa.preflight.checks import Check

__all__ = ["RULE_ID", "tool_genetic_code_check"]

RULE_ID = "config.tool_genetic_code"
RULE_TEXT = (
    "annotation.bakta.extra_args without --translation-table (or an argparse prefix of it of at "
    "least 4 characters), mapping.proteins.miniprot_args without -T (also inside a single-dash "
    "cluster), mapping.proteins.tblastn_args with at most one unquoted -db_gencode equal "
    "to annotation.translation_table, every string split by shlex → ok; else error"
)
BAKTA = "annotation.bakta.extra_args"
MINIPROT = "mapping.proteins.miniprot_args"
TBLASTN = "mapping.proteins.tblastn_args"

_BAKTA_OPTION = "--translation-table"
_BAKTA_MIN_PREFIX = 4


class _Unparseable(ValueError):
    """An option string ``shlex.split`` rejects (for example an unbalanced quote)."""


def _split(args: str) -> list[str]:
    try:
        return shlex.split(args)
    except ValueError as exc:
        raise _Unparseable(str(exc)) from exc


def _bakta_sets_code(args: Sequence[str]) -> bool:
    """A token whose name (before ``=``) is an argparse prefix of ``--translation-table``."""
    for element in args:
        for tok in _split(element):
            name = tok.split("=", 1)[0]
            if len(name) >= _BAKTA_MIN_PREFIX and _BAKTA_OPTION.startswith(name):
                return True
    return False


def _miniprot_sets_code(args: str) -> bool:
    """A single-dash token (a short option or a cluster such as ``-uT4``) holding ``T``."""
    return any(
        tok.startswith("-") and not tok.startswith("--") and "T" in tok[1:] for tok in _split(args)
    )


def _tblastn_conflicts(args: str, table: int) -> bool:
    """``-db_gencode`` repeated, bare at the end, quoted, or with a value other than ``table``
    (``gippa.mapping.tblastn.db_gencode_values``, the rule ``effective_tblastn_args`` applies)."""
    _split(args)
    return db_gencode_values(args, table)[1]


def tool_genetic_code_check(cfg: GippaConfig) -> Check:
    """One ``config.tool_genetic_code`` check; ``value`` is the list of offending config keys."""
    table = cfg.annotation.translation_table
    proteins = cfg.mapping.proteins
    tests: tuple[tuple[str, Callable[[], bool]], ...] = (
        (BAKTA, lambda: _bakta_sets_code(cfg.annotation.bakta.extra_args)),
        (MINIPROT, lambda: _miniprot_sets_code(proteins.miniprot_args)),
        (TBLASTN, lambda: _tblastn_conflicts(proteins.tblastn_args, table)),
    )
    conflicts: list[str] = []
    unparseable: list[str] = []
    for key, test in tests:
        try:
            if test():
                conflicts.append(key)
        except _Unparseable:
            unparseable.append(key)
    keys = [key for key, _ in tests if key in conflicts or key in unparseable]
    if not keys:
        return Check(RULE_ID, "inputs", "ok", RULE_TEXT, [], "")
    parts: list[str] = []
    if conflicts:
        parts.append(
            f"{', '.join(conflicts)} set a genetic code; GIPPA sets every tool's code from "
            f"annotation.translation_table ({table}): remove the option and set "
            "annotation.translation_table instead"
        )
    if unparseable:
        parts.append(
            f"{', '.join(unparseable)} could not be parsed (shell quoting, for example an "
            "unbalanced quote): fix the quoting so the genetic-code options can be checked"
        )
    return Check(
        RULE_ID,
        "inputs",
        "error",
        RULE_TEXT,
        list(keys),
        "; ".join(parts),
        int(ExitCode.INVALID_INPUT),
    )
