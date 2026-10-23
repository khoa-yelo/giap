"""#214 review round: lead decision 6a across collections, and complete engine rows."""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest
from helpers.fgsea_switch_data import RecordingEngine
from helpers.gsea_data import gene_sets, singleton_frame
from test_gsea_engine_switch import _genes, _run, _values


def _two_collections() -> tuple[pd.DataFrame, list[Any]]:
    tags = _genes(40)
    cog = gene_sets("cog", {"K": tags[::4]})
    kegg = gene_sets("kegg", {"WHOLE": tags})
    return singleton_frame(_values(40)), [cog, kegg]


def test_whole_ranking_set_of_a_collection_without_a_job_adds_no_row() -> None:
    genes, colls = _two_collections()
    eng = RecordingEngine()
    out = _run(genes, colls, eng)
    assert [j.job_id for j in eng.calls[0].jobs] == ["main:cog"]
    assert list(out["term"]) == ["K"]
    alone = _run(genes, colls[1:], RecordingEngine())
    assert alone.empty  # the same kegg collection alone: header-only


def test_an_engine_that_drops_a_sent_row_fails_the_call() -> None:
    tags = _genes(40)
    coll = gene_sets("cog", {"J": tags[::4], "K": tags[1::4]})
    eng = RecordingEngine(drop={("main:cog", "K")})
    with pytest.raises(ValueError, match="no row for main:cog K"):
        _run(singleton_frame(_values(40)), [coll], eng)
