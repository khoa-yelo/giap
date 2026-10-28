"""FarmShare CPUs from memory use MaxMemPerCPU = 4000 MB, not 4 GB (#416 D4).

farmshare.config asked for max(cpus, ceil(memory / 4 GB)) CPUs, but Slurm's MaxMemPerCPU is
4000 MB (MiB), so a 12 CPU / 48 GB request (49152 MB / 4000 = 12.3) got a 13th CPU from Slurm.
Every memory-derived CPU rule in farmshare.config divides by ``4000.MB``; the derived count for
48 GB is 13, the CPUs Slurm allocates.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pytest

CONF = Path(__file__).resolve().parents[2] / "workflow" / "conf"
FARMSHARE = CONF / "farmshare.config"
UNITS = {"KB": 2**10, "MB": 2**20, "GB": 2**30}
RULE = re.compile(
    r"Math\.max\(\s*(\d+|params\.terminatornet_chunk_cpus)\s*,\s*Math\.ceil\(\s*task\.memory\.toBytes\(\)\s*/\s*"
    r"(\d+)\.(KB|MB|GB)\.toBytes\(\)\s*\)\s*as\s+int\s*\)"
)
#: a param base reads as its default (M2-15b: TERMINATORNET_CHUNK)
BASES = {"params.terminatornet_chunk_cpus": 4}
MAX_MEM_PER_CPU = 4000 * 2**20  # Slurm MaxMemPerCPU=4000 (MB = MiB)


def _rules() -> list[tuple[int, int]]:
    text = re.sub(r"(^|\s)//[^\n]*", r"\1", FARMSHARE.read_text(encoding="utf-8"))
    rules = [(BASES.get(m[1]) or int(m[1]), int(m[2]) * UNITS[m[3]]) for m in RULE.finditer(text)]
    assert len(rules) >= 6, f"too few memory-derived CPU rules: {rules}"
    assert len(rules) == text.count("Math.ceil("), "a CPU rule the pattern does not read"
    return rules


def _cpus(base: int, divisor: int, memory: int) -> int:
    return max(base, math.ceil(memory / divisor))


def test_every_rule_divides_by_4000_mb() -> None:
    assert {d for _, d in _rules()} == {MAX_MEM_PER_CPU}


def test_no_four_gb_divisor_left() -> None:
    text = re.sub(r"(^|\s)//[^\n]*", r"\1", FARMSHARE.read_text(encoding="utf-8"))
    assert not re.search(r"/\s*4\.GB", text)


@pytest.mark.parametrize("base", [1, 2, 8, 12])
def test_48_gb_derives_13_cpus(base: int) -> None:
    (divisor,) = {d for _, d in _rules()}
    assert _cpus(base, divisor, 48 * 2**30) == 13  # ceil(49152 / 4000)
    assert math.ceil(48 * 1024 / 4000) == 13


def test_derived_cpus_never_exceed_max_mem_per_cpu() -> None:
    """For every rule and memory 1..64 GB, memory / cpus <= 4000 MB (Slurm adds no CPU)."""
    for base, divisor in _rules():
        for gb in range(1, 65):
            mem = gb * 2**30
            assert mem / _cpus(base, divisor, mem) <= MAX_MEM_PER_CPU, (base, gb)
