"""``scripts/dev/progress_bars.py``: specs superseded by an accepted-by-khoa ADR (dev tooling)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts/dev/progress_bars.py"


@pytest.fixture(scope="module")
def pb() -> ModuleType:
    spec = importlib.util.spec_from_file_location("progress_bars", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["progress_bars"] = mod
    spec.loader.exec_module(mod)
    return mod


ACCEPTED = "# ADR-0038: x\n- Status: accepted-by-khoa (2026-10-22)\n"
PROPOSED = "# ADR-0099: y\n- Status: proposed\n"


def test_superseded_by_reads_the_head_only(pb: ModuleType) -> None:
    assert pb.superseded_by("# M8-11: t\n\nStatus: superseded by ADR-0038 (why)\n") == "ADR-0038"
    assert pb.superseded_by("# t\n- **Status:** Superseded by adr-0038\n") == "ADR-0038"
    assert pb.superseded_by("# t\n" + "x\n" * 12 + "Status: superseded by ADR-0038\n") is None
    assert pb.superseded_by("# t\nThis is not superseded by ADR-0038.\n") is None


def test_accepted_by_khoa(pb: ModuleType) -> None:
    assert pb.accepted_by_khoa(ACCEPTED)
    assert not pb.accepted_by_khoa(PROPOSED)
    assert not pb.accepted_by_khoa("# ADR\n- Status: accepted (lead)\n")


def test_superseded_requires_an_accepted_adr(pb: ModuleType) -> None:
    specs = [
        ("M8-11", "# M8-11\nStatus: superseded by ADR-0038\n"),
        ("M8-12", "# M8-12\nStatus: superseded by ADR-0099\n"),
        ("M8-13", "# M8-13\nStatus: superseded by ADR-0777\n"),
        ("M8-14", "# M8-14\nplain spec\n"),
    ]
    out, warnings = pb.superseded(specs, {"ADR-0038": ACCEPTED, "ADR-0099": PROPOSED})
    assert out == {"M8-11"}
    assert len(warnings) == 2
    assert any("M8-12" in w and "ADR-0099" in w for w in warnings)
    assert any("M8-13" in w and "not found" in w for w in warnings)


def test_render_table_counts_superseded_as_done(pb: ModuleType) -> None:
    issues = {8: ["M8-01", "M8-02", "M8-03"]}
    block = pb.render_table(issues, {"M8-01", "M8-02"}, set(), {"M8-03"})
    row = next(line for line in block.splitlines() if line.startswith("| M8 "))
    assert row == (
        f"| M8 Function and determinants | {pb.DONE * 2}{pb.SUPERSEDED} 100 % "
        "| 3 (1 superseded) | 0 | 0 |"
    )
    assert f"{pb.SUPERSEDED} superseded" in block  # legend
    assert "designed issues: 4/4 done (100 %)" in block  # M0 + three M8 issues


def test_render_table_superseded_beats_running(pb: ModuleType) -> None:
    block = pb.render_table({8: ["M8-01", "M8-02"]}, set(), {"M8-02"}, {"M8-02"})
    row = next(line for line in block.splitlines() if line.startswith("| M8 "))
    assert row.endswith("| 1 (1 superseded) | 0 | 1 |")
