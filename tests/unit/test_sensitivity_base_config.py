"""M8-16 review round, code review B1 (lead decision 5): ``workflow/conf/base.config`` sizes the
``sensitivity`` step (``SENSITIVITY_RERUNS:SENSITIVITY``) at least like ``FUNCTION_STEP``: it runs
R fgsea and the LightGBM CV reruns, so the default 1 CPU / 2 GB would be killed (exit 137).
Decision 5: 6 GB, 4 h, process_low CPUs (2). Static check, as ``test_m8_17b_wiring`` does for
``FUNCTION_STEP`` and ``test_determinants_step`` for ``DETERMINANTS:DETERMINANTS_STEP``.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CONF = REPO / "workflow" / "conf" / "base.config"
NF = REPO / "workflow" / "subworkflows" / "local" / "sensitivity.nf"
QUALIFIED = "SENSITIVITY_RERUNS:SENSITIVITY"


def _strip_comments(text: str) -> str:
    return "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("//"))


def _selectors(text: str) -> list[tuple[str, str]]:
    """``(pattern, body)`` of every ``withName: '<pattern>' { ... }`` (braces balanced)."""
    out = []
    for m in re.finditer(r"withName\s*:\s*['\"]([^'\"]+)['\"]\s*\{", text):
        depth, i = 1, m.end()
        while depth:
            depth += text[i] == "{"
            depth -= text[i] == "}"
            i += 1
        out.append((m.group(1), text[m.end() : i - 1]))
    return out


def _matches(pat: str, simple: str, qualified: str) -> bool:
    return bool(re.fullmatch(pat, simple) or re.fullmatch(pat, qualified))


def test_the_step_runs_as_sensitivity_inside_sensitivity_reruns() -> None:
    text = NF.read_text(encoding="utf-8")
    assert re.search(r"GIPPA_STEP\s+as\s+SENSITIVITY\b", text)
    assert re.search(r"workflow\s+SENSITIVITY_RERUNS\b", text)


def test_base_config_sizes_sensitivity_like_function_step() -> None:
    sel = _selectors(_strip_comments(CONF.read_text(encoding="utf-8")))
    hits = [(p, b) for p, b in sel if _matches(p, "SENSITIVITY", QUALIFIED)]
    assert hits, f"no withName selector for {QUALIFIED} in base.config: {[p for p, _ in sel]}"
    for pat, _ in hits:
        for other in ("GIPPA_STEP", "FUNCTION_STEP", "DETERMINANTS_STEP", "CONTEXT_STATS"):
            assert not _matches(pat, other, f"CONTEXT:{other}"), f"{pat} also matches {other}"
    body = hits[-1][1]  # the last matching selector wins in Nextflow
    mem = re.search(r"\bmemory\s*=\s*\{?\s*(\d+(?:\.\d+)?)\s*\.\s*GB", body)
    assert mem and float(mem.group(1)) >= 6, body
    assert re.search(r"task\.attempt", body), f"memory does not grow on retry: {body}"
    time = re.search(r"\btime\s*=\s*\{?\s*(\d+(?:\.\d+)?)\s*\.\s*h\b", body)
    assert time and float(time.group(1)) >= 4, body
    cpus = re.search(r"\bcpus\s*=\s*(\d+)\b", body)
    assert cpus and int(cpus.group(1)) >= 2, f"process_low CPUs (2) expected: {body}"
