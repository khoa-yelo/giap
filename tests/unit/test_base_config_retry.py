"""base.config retry policy (M7-15, Khoa 2026-10-23): exit 137/143 retries twice, memory grows
with ``task.attempt``, and no process selector overrides the policy.

Review #278 B3 (lead decision 11): the selector check parses balanced braces, so a selector whose
body holds a closure (``memory = { ... }``) is read to its end, and it reads every workflow
config (base, farmshare, farmshare_local, the profiles in nextflow.config), not base.config
only; module directives are checked too. The resolved per-profile view is in
tests/integration/test_m7_15_retry_policy.py (``nextflow config -flat``).
"""

from __future__ import annotations

import re
from pathlib import Path

WF = Path(__file__).resolve().parents[2] / "workflow"
BASE = WF / "conf" / "base.config"
CONFIGS = sorted([WF / "nextflow.config", *(WF / "conf").glob("*.config")])
POLICY = re.compile(r"\b(errorStrategy|maxRetries)\b")
_SELECTOR = re.compile(r"with(?:Name|Label)\s*:\s*(['\"])(.*?)\1\s*\{")


def _text() -> str:
    return BASE.read_text(encoding="utf-8")


def _strip_comments(text: str) -> str:
    return re.sub(r"(^|\s)//[^\n]*", r"\1", text)  # not the // of docker:// URIs


def _block(text: str, open_at: int) -> str:
    """The text between the ``{`` at ``open_at`` and its matching ``}`` (braces balanced)."""
    assert text[open_at] == "{"
    depth = 0
    for i in range(open_at, len(text)):
        depth += {"{": 1, "}": -1}.get(text[i], 0)
        if depth == 0:
            return text[open_at + 1 : i]
    raise AssertionError(f"unbalanced braces from offset {open_at}")


def selector_bodies(text: str) -> list[tuple[str, str]]:
    """(selector, whole body) for every withName/withLabel block, closures included."""
    text = _strip_comments(text)
    return [(m.group(2), _block(text, m.end() - 1)) for m in _SELECTOR.finditer(text)]


def test_retries_137_and_143_twice() -> None:
    text = _strip_comments(_text())
    assert re.search(r"^\s*maxRetries\s*=\s*2\s*$", text, re.M)
    strategy = re.search(r"^\s*errorStrategy\s*=\s*\{", text, re.M)
    assert strategy, "no default errorStrategy closure"
    body = _block(text, strategy.end() - 1)  # a multi-line closure is read to its end
    assert re.search(r"\[\s*137\s*,\s*143\s*\]", body) and "'retry'" in body
    assert "'terminate'" in body, "every other exit terminates"


def test_default_memory_grows_with_the_attempt() -> None:
    assert re.search(r"^\s*memory\s*=\s*\{\s*\d+\.GB\s*\*\s*task\.attempt\s*\}", _text(), re.M)


def test_selector_parser_reads_past_a_closure() -> None:
    """The review's probe: the old ``\\{([^}]*)\\}`` regex stopped inside the memory closure."""
    probe = "withName: 'VARIANTS_STEP' { memory = { 8.GB * task.attempt }; maxRetries = 0 }"
    ((name, body),) = selector_bodies(probe)
    assert name == "VARIANTS_STEP" and POLICY.search(body)
    nested = (
        'profiles { x { process { withLabel: "process_low" { cpus = { 2 }\n'
        "  errorStrategy = { task.attempt < 9 ? 'retry' : 'ignore' } } } } }"
    )
    ((name, body),) = selector_bodies(nested)
    assert name == "process_low" and POLICY.search(body)


def test_no_selector_overrides_the_retry_policy() -> None:
    seen = 0
    for path in CONFIGS:
        for name, body in selector_bodies(path.read_text(encoding="utf-8")):
            seen += 1
            assert not POLICY.search(body), f"{path.name}: {name} overrides the retry policy"
    assert seen >= 10, "selector parser found too few blocks"


def test_no_module_directive_overrides_the_retry_policy() -> None:
    nf = [p for p in WF.rglob("*.nf") if "tests" not in p.relative_to(WF).parts]
    assert nf
    for path in nf:
        text = _strip_comments(path.read_text(encoding="utf-8"))
        assert not POLICY.search(text), f"{path.relative_to(WF)} sets errorStrategy/maxRetries"
