"""M2 gate: every §7.1.1 annotation tool has an active contract test (M2-21 acceptance item 1).

"Every tool in §7.1.1 except miniprot/BLAST+ (M3) has an active contract test on a committed
recorded fixture (no ``pending_sample`` left in ``tests/contract/`` for annotation tools)."

Still pending, with the reason (checked so the list cannot silently grow):

* ``blastp``, ``tblastn``, ``miniprot`` (M3), ``diann`` (M7), ``m4_input`` (M4-21): not
  annotation tools.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
CONTRACT = REPO / "tests" / "contract"
RECORDED = REPO / "tests" / "data" / "recorded"
ANNOTATION_TOOLS = ("bakta", "genomad", "antismash", "terminatornet", "transtermhp", "ostir",
                    "promoter_calculator", "promotech", "uniop")  # fmt: skip
STILL_PENDING = {"blastp", "tblastn", "miniprot", "diann", "m4_input"}
PENDING = re.compile(r"pytest\.mark\.pending_sample\b")


def _pending_files() -> set[str]:
    return {p.name.removeprefix("test_").removesuffix("_contract.py")
            for p in CONTRACT.glob("test_*.py") if PENDING.search(p.read_text())}  # fmt: skip


@pytest.mark.parametrize("tool", ANNOTATION_TOOLS)
def test_annotation_tool_contract_is_active(tool: str) -> None:
    path = CONTRACT / f"test_{tool}_contract.py"
    assert path.is_file(), path
    text = path.read_text()
    assert not PENDING.search(text), f"{path.name} still carries pending_sample"
    assert not re.search(r"skipif\([^)]*pending_sample", text, re.S), path.name


def test_only_the_listed_contracts_remain_pending() -> None:
    assert _pending_files() <= STILL_PENDING, sorted(_pending_files() - STILL_PENDING)


def test_bactermfinder_is_no_longer_pending() -> None:
    """A BacTermFinder sample is recorded (fix-bactermfinder-parser): the marker is gone."""
    text = (CONTRACT / "test_bactermfinder_contract.py").read_text()
    assert not PENDING.search(text)


@pytest.mark.parametrize(
    "fixture",
    [
        "promotech/56251ad9/k12_mg1655_excerpt",
        "uniop/3ed14dfc/k12",
        "uniop/3ed14dfc/sim100k",
    ],
)
def test_new_excerpts_carry_provenance(fixture: str) -> None:
    import hashlib

    d = RECORDED / fixture
    prov = (d / "PROVENANCE.md").read_text()
    for f in sorted(d.iterdir()):
        if f.name == "PROVENANCE.md":
            continue
        assert hashlib.sha256(f.read_bytes()).hexdigest() in prov, f"{fixture}/{f.name}"
