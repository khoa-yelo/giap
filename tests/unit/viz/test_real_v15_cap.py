"""K-12 plot survey, batch 1 item 4: V15's SHAP axes stay within the 20,000-element cap.

On K-12 the checklist counted 20,001 elements on the SHAP axes: ``MAX_DRAWN`` points plus the
zero line. Required: the drawn total on the axes, extras (reference lines, texts, patches)
included, is at most 20,000, counted as ``giap.viz.checklist`` counts it (patches + lines + texts
+ each collection's offsets / paths, untagged artists together).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from helpers import viz_determinants as vd
from helpers.viz_landscape import plot_spec

from giap.config import resolve_config
from giap.io.run_view import RunView
from giap.viz import checklist as ck
from giap.viz import style
from giap.viz.base import Skip

CAP = 20_000


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


def _count(ax) -> int:  # type: ignore[no-untyped-def]
    n = len(ax.patches) + len(ax.lines) + len(ax.texts)
    for c in ax.collections:
        n += max(len(np.asarray(c.get_offsets())), len(c.get_paths()), 1)
    return n


@pytest.mark.parametrize("n_genes", [1_001, 1_500])
def test_shap_axes_total_is_within_the_cap(tmp_path: Path, n_genes: int) -> None:
    """20 features x >= 1,001 genes: the per-feature share fills the cap exactly."""
    shap = vd.shap_frame(tags=[f"T{i:05d}" for i in range(n_genes)])
    out = vd.build_out(tmp_path, shap=shap)
    result = plot_spec("V15").fn(RunView.from_outdir(out), resolve_config())
    assert not isinstance(result, Skip), result
    fig = result[0].figure
    (ax,) = [a for a in fig.axes if a.get_label() == "c"]
    assert _count(ax) <= CAP, _count(ax)
    ok, detail = ck._elements(fig)
    assert ok, detail
