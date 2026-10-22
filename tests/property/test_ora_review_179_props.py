"""Properties for the PR #179 review blockers (docs/specs/M8-04 lead decision 13).

- Any small universe run with ``folds=None`` flags every row ``tu_recal_skipped`` and
  ``ora_uncalibrated``.
- Any eligible gene whose ``p_detect`` lies outside [0, 1] makes ``detection_enrichment`` raise
  ``ValueError`` naming that locus tag; the same value on a partial gene is ignored.
"""

from __future__ import annotations

import re
from typing import Any

import pytest
from helpers.gsea_data import as_list
from helpers.ora_data import Gene, frame, gene_sets, ora
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

SETTINGS = settings(max_examples=40, deadline=None, suppress_health_check=list(HealthCheck))
SIZE = (1, 500)
ALPHA = 0.05

_p = st.floats(min_value=0.05, max_value=0.95, allow_nan=False)
_bad = st.one_of(
    st.floats(min_value=1.0, max_value=10.0, exclude_min=True, allow_nan=False),
    st.floats(min_value=-10.0, max_value=0.0, exclude_max=True, allow_nan=False),
)


@st.composite
def universes(draw: Any) -> tuple[dict[str, list[Gene]], list[str]]:
    tus: dict[str, list[Gene]] = {}
    tags: list[str] = []
    for i in range(draw(st.integers(min_value=3, max_value=10))):
        genes = []
        for _ in range(draw(st.integers(min_value=1, max_value=3))):
            tag = f"g{len(tags) + 1}"
            det = draw(st.booleans())
            genes.append(Gene(tag, 10.0 if det else None, det, p_detect=draw(_p)))
            tags.append(tag)
        tus[f"TU_{i + 1}"] = genes
    return tus, tags


def _colls(draw: Any, tags: list[str]) -> list[Any]:
    sets = {
        f"T{j}": set(draw(st.lists(st.sampled_from(tags), min_size=1, max_size=5)))
        for j in range(draw(st.integers(min_value=1, max_value=3)))
    }
    return [gene_sets("kegg", sets)]


@SETTINGS
@given(data=st.data())
def test_folds_none_flags_every_row_skipped_and_uncalibrated(data: Any) -> None:
    tus, tags = data.draw(universes())
    out = ora().detection_enrichment(frame(tus), _colls(data.draw, tags), SIZE, ALPHA, folds=None)
    for f in out["flags"]:
        assert {"tu_recal_skipped", "ora_uncalibrated"} <= set(as_list(f)), as_list(f)


@SETTINGS
@given(data=st.data())
def test_out_of_range_p_detect_fails_closed_naming_the_gene(data: Any) -> None:
    tus, tags = data.draw(universes())
    tag = data.draw(st.sampled_from(tags))
    g = frame(tus)
    g.loc[g["locus_tag"] == tag, "p_detect"] = data.draw(_bad)
    colls = _colls(data.draw, tags)
    with pytest.raises(ValueError, match=rf"\b{re.escape(tag)}\b"):
        ora().detection_enrichment(g, colls, SIZE, ALPHA, folds=None)
    g.loc[g["locus_tag"] == tag, "partial"] = True
    if (~g["partial"].astype(bool)).any():
        ora().detection_enrichment(g, colls, SIZE, ALPHA, folds=None)
