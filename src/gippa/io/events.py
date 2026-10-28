"""``OUT/logs/events.jsonl`` (GIPPA_SPEC.md §6.7): one JSON object per event (``ts``, ``step``,
``level``, ``event``, ``data``). Steps record events through ``ctx.log`` with
``extra={"event": name, "data": {...}}``; ``EventsHandler`` writes them."""

from __future__ import annotations

import json
import logging
from collections.abc import Mapping
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


def event_record(
    step: str, event: str, data: Mapping[str, Any], level: str = "INFO", ts: datetime | None = None
) -> dict[str, Any]:
    """One §6.7 event object."""
    when = (ts or datetime.now(UTC)).isoformat(timespec="seconds")
    return {"ts": when, "step": step, "level": level, "event": event, "data": dict(data)}


def append_event(logs_dir: Path, record: Mapping[str, Any]) -> None:
    """Append ``record`` as one line to ``logs_dir/events.jsonl``."""
    logs_dir.mkdir(parents=True, exist_ok=True)
    with (logs_dir / "events.jsonl").open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(record, sort_keys=True, default=str) + "\n")


class EventsHandler(logging.Handler):
    """Log records with an ``event`` extra → ``logs_dir/events.jsonl`` (step = logger leaf)."""

    def __init__(self, logs_dir: Path) -> None:
        super().__init__(logging.INFO)
        self.logs_dir = logs_dir

    def emit(self, record: logging.LogRecord) -> None:
        event, data = getattr(record, "event", None), getattr(record, "data", None)
        if isinstance(event, str):
            ts = datetime.fromtimestamp(record.created, UTC)
            step = record.name.rsplit(".", 1)[-1]
            data = data if isinstance(data, Mapping) else {}
            append_event(self.logs_dir, event_record(step, event, data, record.levelname, ts))
