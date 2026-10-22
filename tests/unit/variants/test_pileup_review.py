"""M7-02 review (#158): ``sites`` refuses a pileup frame that lost its ``attrs``, since the
frameshift events travel only there and would otherwise vanish without an error."""

from __future__ import annotations

import pandas as pd
import pytest
from unit.variants.test_pileup import CFG, G1, al, gene_hit, pl, run


@pytest.mark.parametrize("drop", ["frameshifts", "contig_order"])
def test_sites_raises_when_pileup_attrs_are_missing(drop: str) -> None:
    pile, _ = run([gene_hit(al(), G1, 2, 8, {4: "D"})])
    rebuilt = pd.DataFrame(pile.to_dict("records"))
    with pytest.raises(ValueError, match="attrs"):
        pl().sites(rebuilt, CFG)
    pile.attrs.pop(drop)
    with pytest.raises(ValueError, match="attrs"):
        pl().sites(pile, CFG)
