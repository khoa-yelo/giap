"""Preflight check records (GIPPA_SPEC.md §7.0; ADR-0005 item 3).

Every evaluated rule yields one ``Check``; passing rules are recorded with ``level="ok"``.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal, TypeAlias

Level = Literal["ok", "info", "warn", "error"]
InputName = Literal["inputs", "genome", "proteins", "peptides", "abundance", "environment"]
JsonValue: TypeAlias = "str | int | float | bool | list[JsonValue] | dict[str, JsonValue] | None"

LEVEL_ORDER: tuple[Level, ...] = ("ok", "info", "warn", "error")
_RANK: dict[str, int] = {lv: i for i, lv in enumerate(LEVEL_ORDER)}

__all__ = ["LEVEL_ORDER", "LIST_CAP", "Check", "InputName", "JsonValue", "Level", "capped", "worst"]

LIST_CAP = 1000


@dataclass(frozen=True, slots=True)
class Check:
    """One evaluated preflight rule.

    ``exit_code`` is 1 for an input error, 2 for an environment error and 0 otherwise; an error
    created with the default 0 gets the code its ``input`` implies.
    """

    rule_id: str
    input: InputName
    level: Level
    rule: str
    value: JsonValue
    message: str = ""
    exit_code: int = 0

    def __post_init__(self) -> None:
        # An error left at the default 0 takes the code its input implies (ADR-0005 item 3).
        if self.level == "error" and self.exit_code == 0:
            object.__setattr__(self, "exit_code", 2 if self.input == "environment" else 1)


def worst(checks: Iterable[Check]) -> Level:
    """The most severe level among ``checks`` (``"ok"`` when empty)."""
    level: Level = "ok"
    for c in checks:
        if _RANK[c.level] > _RANK[level]:
            level = c.level
    return level


def capped(items: list[JsonValue]) -> tuple[list[JsonValue], str]:
    """At most ``LIST_CAP`` items, plus a message giving the total when the list was cut."""
    if len(items) <= LIST_CAP:
        return items, ""
    return items[:LIST_CAP], f"{len(items)} entries in total; first {LIST_CAP} listed"
