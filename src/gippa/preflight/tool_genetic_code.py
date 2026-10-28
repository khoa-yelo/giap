"""``config.tool_genetic_code``: no tool option may set a genetic code itself (M11-16 item 2.2).

GIPPA sets every tool's genetic code from ``annotation.translation_table``, so a Bakta
``--translation-table``, a miniprot ``-T`` or a tblastn ``-db_gencode`` that disagrees with the
table (or is repeated) is an input error. Pure: reads the config only.
"""

from __future__ import annotations

from collections.abc import Sequence

from gippa.config import GippaConfig
from gippa.config.model import ProteinsConfig
from gippa.exitcodes import ExitCode
from gippa.preflight.checks import Check

__all__ = ["RULE_ID", "tool_genetic_code_check"]

RULE_ID = "config.tool_genetic_code"
RULE_TEXT = (
    "annotation.bakta.extra_args without --translation-table, mapping.proteins.miniprot_args "
    "without -T, mapping.proteins.tblastn_args with at most one -db_gencode equal to "
    "annotation.translation_table → ok; else error"
)
BAKTA = "annotation.bakta.extra_args"
MINIPROT = "mapping.proteins.miniprot_args"
TBLASTN = "mapping.proteins.tblastn_args"

#: The shipped ``tblastn_args`` default (it still carries ``-db_gencode 11``); never flagged.
_TBLASTN_DEFAULT: str = ProteinsConfig.model_fields["tblastn_args"].default


def _bakta_sets_code(args: Sequence[str]) -> bool:
    return any(a == "--translation-table" or a.startswith("--translation-table=") for a in args)


def _miniprot_sets_code(args: str) -> bool:
    return any(a.startswith("-T") for a in args.split())


def _tblastn_conflicts(args: str, table: int) -> bool:
    """``-db_gencode`` repeated, or one whose value is not ``table``.

    The shipped default is exempt at every table: M11-16d derives ``-db_gencode`` from the table
    for it (spec 3.3, D1 fallback), and that change removes this exemption.
    """
    if args == _TBLASTN_DEFAULT:
        return False
    tokens = args.split()
    values: list[str | None] = []
    for i, tok in enumerate(tokens):
        if tok == "-db_gencode":
            values.append(tokens[i + 1] if i + 1 < len(tokens) else None)
        elif tok.startswith("-db_gencode="):
            values.append(tok.split("=", 1)[1])
    return len(values) > 1 or any(v != str(table) for v in values)


def tool_genetic_code_check(cfg: GippaConfig) -> Check:
    """One ``config.tool_genetic_code`` check; ``value`` is the list of offending config keys."""
    table = cfg.annotation.translation_table
    proteins = cfg.mapping.proteins
    keys = [
        key
        for key, bad in (
            (BAKTA, _bakta_sets_code(cfg.annotation.bakta.extra_args)),
            (MINIPROT, _miniprot_sets_code(proteins.miniprot_args)),
            (TBLASTN, _tblastn_conflicts(proteins.tblastn_args, table)),
        )
        if bad
    ]
    if not keys:
        return Check(RULE_ID, "inputs", "ok", RULE_TEXT, [], "")
    message = (
        f"{', '.join(keys)} set a genetic code; GIPPA sets every tool's code from "
        f"annotation.translation_table ({table}): remove the option and set "
        "annotation.translation_table instead"
    )
    return Check(
        RULE_ID, "inputs", "error", RULE_TEXT, list(keys), message, int(ExitCode.INVALID_INPUT)
    )
