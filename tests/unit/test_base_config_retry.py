"""base.config retry policy (M7-15, Khoa 2026-10-23): exit 137/143 retries twice, memory grows
with ``task.attempt``, and no process selector overrides the policy."""

from __future__ import annotations

import re
from pathlib import Path

BASE = Path(__file__).resolve().parents[2] / "workflow" / "conf" / "base.config"


def _text() -> str:
    return BASE.read_text(encoding="utf-8")


def test_retries_137_and_143_twice() -> None:
    text = _text()
    assert re.search(r"^\s*maxRetries\s*=\s*2\s*$", text, re.M)
    strategy = re.search(r"^\s*errorStrategy\s*=\s*\{(.*)\}\s*$", text, re.M)
    assert strategy, "no default errorStrategy closure"
    body = strategy.group(1)
    assert re.search(r"\[\s*137\s*,\s*143\s*\]", body) and "'retry'" in body
    assert "'terminate'" in body, "every other exit terminates"


def test_default_memory_grows_with_the_attempt() -> None:
    assert re.search(r"^\s*memory\s*=\s*\{\s*\d+\.GB\s*\*\s*task\.attempt\s*\}", _text(), re.M)


def test_no_selector_overrides_the_retry_policy() -> None:
    selectors = re.findall(r"with(?:Name|Label):\s*'[^']+'\s*\{([^}]*)\}", _text())
    assert selectors
    for body in selectors:
        assert "errorStrategy" not in body and "maxRetries" not in body, body
