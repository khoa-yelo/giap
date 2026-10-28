"""M9-24 acceptance: every registered plot on the snapshot fixtures passes every gated item of the
§14.4 style checklist (the score M9-25's approval packet reports; GIPPA_SPEC.md §14.4).

The fixtures are M9-25's frozen ``OUT/`` subsets under ``tests/data/snapshot/`` (the rich
``sim_plots`` and the header-only empty-state tree); each sub-directory that is a GIPPA output
(``tables/`` or ``preflight/``, directly or under ``OUT/``) is swept. Until M9-25 commits them the
parametrisation is a single explicit skip; ``tests/integration/test_checklist_sim100k.py`` runs the
same sweep on a real sim100k run meanwhile.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from helpers.checklist_sweep import sweep

from gippa.viz import checklist  # noqa: F401  (fails collection until M9-24 lands)

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
def test_every_registered_plot_passes_every_gated_item(out: Path, tmp_path: Path) -> None:
    failures, _ = sweep(out, tmp_path / "plots")
    assert failures == {}
