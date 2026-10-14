"""Validation findings shared by the io validators and preflight (GIAP_SPEC.md §7.0)."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

Severity = Literal["error", "warn", "info"]

MAX_EXAMPLES = 5

__all__ = ["MAX_EXAMPLES", "Issue", "Severity", "row_issue"]


@dataclass(frozen=True, slots=True)
class Issue:
    """One finding; ``examples`` holds at most five 1-based data-row numbers (as strings)."""

    code: str
    severity: Severity
    message: str
    n_rows: int = 0
    examples: tuple[str, ...] = ()


def row_issue(code: str, severity: Severity, message: str, rows: Iterable[int]) -> Issue | None:
    """Issue for the given 1-based ``rows`` (``None`` if there are none)."""
    rows = list(rows)
    if not rows:
        return None
    examples = tuple(str(r) for r in rows[:MAX_EXAMPLES])
    return Issue(code, severity, f"{message} ({len(rows)} rows)", len(rows), examples)
