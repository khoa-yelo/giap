"""M8-03 review (#159), ADR-0023: every KEGG "Global and overview maps" entry is excluded,
including the nitrogen (map01310) and sulfur (map01320) cycle maps outside map011xx/map012xx."""

from __future__ import annotations

import pandas as pd
import pytest

from giap.function.collections import kegg_pathways


@pytest.mark.parametrize("overview", ["map01100", "map01250", "map01310", "map01320"])
def test_overview_maps_excluded(overview: str) -> None:
    links = pd.DataFrame(
        {"ko": ["ko:K00001", "ko:K00001"], "pathway": [f"path:{overview}", "path:map00010"]}
    )
    assert kegg_pathways(links) == {"K00001": frozenset({"map00010"})}


@pytest.mark.parametrize("kept", ["map01010", "map01300", "map01330", "map00910"])
def test_neighbouring_maps_kept(kept: str) -> None:
    links = pd.DataFrame({"ko": ["ko:K00002"], "pathway": [f"path:{kept}"]})
    assert kegg_pathways(links) == {"K00002": frozenset({kept})}
