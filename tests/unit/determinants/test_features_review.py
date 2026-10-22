"""M8-07 review (#162): a duplicated locus tag in a per-gene source frame keeps its first row
instead of raising a length mismatch."""

from __future__ import annotations

import pandas as pd

from giap.determinants import features


def test_lookup_with_duplicate_tags_keeps_first() -> None:
    frame = pd.DataFrame({"locus_tag": ["a", "a", "b"], "x": [1.0, 2.0, 3.0]})
    out = features._lookup(frame, "x", pd.Series(["b", "a", "c"]))
    assert out.tolist()[:2] == [3.0, 1.0]
    assert pd.isna(out.iloc[2])
