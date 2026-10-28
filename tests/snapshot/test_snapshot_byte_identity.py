"""M9-26 work item 4 on the M9-25 snapshot fixture: every registered plot rendered twice (two
fresh interpreters, different ``PYTHONHASHSEED``s) gives byte-identical ``.svg``, ``.pdf``,
``.png`` and data TSVs (docs/specs/M9-26-scale-and-determinism.md example 2; M9-01 item 7;
GIPPA_SPEC.md R7/TE7). Within one environment only.

The fixtures are M9-25's frozen ``OUT/`` subsets under ``tests/data/snapshot/`` (as
tests/snapshot/test_checklist_all.py sweeps them). Until M9-25 commits them the
parametrisation is a single explicit skip; tests/unit/viz/test_plot_byte_identity.py runs the same
comparison on the committed sim50k golden tables meanwhile.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from helpers.byte_identity import differing, render_in_subprocess, svg_diff_report, tree_bytes

SNAPSHOT = Path(__file__).resolve().parents[1] / "data" / "snapshot"


def _fixtures() -> list[Path]:
    if not SNAPSHOT.is_dir():
        return []
    outs = []
    for d in sorted(p for p in SNAPSHOT.iterdir() if p.is_dir()):
        out = d / "OUT" if (d / "OUT").is_dir() else d
        if (out / "tables").is_dir() or (out / "preflight").is_dir():
            outs.append(out)
    return outs


FIXTURES = _fixtures()


@pytest.mark.parametrize(
    "out",
    FIXTURES
    or [pytest.param(None, marks=pytest.mark.skip(reason="M9-25 snapshot fixtures not committed"))],
    ids=[p.parent.name if p.name == "OUT" else p.name for p in FIXTURES] or ["absent"],
)
def test_every_plot_byte_identical_on_the_snapshot_fixture(out: Path, tmp_path: Path) -> None:
    a = render_in_subprocess(out, tmp_path / "a", hashseed=3)
    b = render_in_subprocess(out, tmp_path / "b", hashseed=4)
    assert a == b
    assert any(v == "rendered" for v in a.values()), a
    ta, tb = tree_bytes(tmp_path / "a"), tree_bytes(tmp_path / "b")
    assert differing(ta, tb) == [], svg_diff_report(ta, tb)
