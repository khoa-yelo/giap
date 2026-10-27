"""TerminatorNet and Promotech resource requests (ADR-0062; #410 review B2).

Static pins on ``workflow/conf/base.config`` and ``farmshare.config``, as for the other resource
``withName`` blocks: TerminatorNet gets the process_high CPUs at 16 GB x attempt, Promotech (and
its novel-ORF alias) 2 CPUs at 10 GB x attempt, and FarmShare keeps its max(cpus, ceil(memory /
4000 MB)) rule (MaxMemPerCPU, #416 D4) for the Promotech selector (base.config's withName CPUs
outrank the label rule), i.e. 3 / 6 / 8 CPUs at attempts 1 / 2 / 3.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

CONF = Path(__file__).resolve().parents[2] / "workflow" / "conf"


def _block(text: str, selector: str) -> str:
    m = re.search(rf"withName:\s*'{re.escape(selector)}'\s*\{{(.*)\}}\s*$", text, re.M)
    assert m, f"no withName: '{selector}' block"
    return m[1]


def test_terminatornet_request() -> None:
    body = _block((CONF / "base.config").read_text(), "TERMINATORNET")
    assert re.search(r"\bcpus\s*=\s*16\b", body), body
    assert re.search(r"memory\s*=\s*\{\s*16\.GB\s*\*\s*task\.attempt\s*\}", body), body


def test_promotech_request() -> None:
    body = _block((CONF / "base.config").read_text(), "PROMOTECH|PROMOTECH_NOVEL")
    assert re.search(r"\bcpus\s*=\s*2\b", body), body
    assert re.search(r"memory\s*=\s*\{\s*10\.GB\s*\*\s*task\.attempt\s*\}", body), body


def test_farmshare_keeps_four_gb_per_cpu_for_promotech() -> None:
    body = _block((CONF / "farmshare.config").read_text(), "PROMOTECH|PROMOTECH_NOVEL")
    m = re.search(
        r"Math\.max\(\s*(\d+)\s*,\s*Math\.ceil\(task\.memory\.toBytes\(\)\s*/\s*4000\.MB", body
    )
    assert m and int(m[1]) == 2, body
    assert [max(2, math.ceil(10 * 1024 * a / 4000)) for a in (1, 2, 3)] == [3, 6, 8]
