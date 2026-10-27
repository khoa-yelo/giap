"""M6-04 conformance on the merged sim100k producer chain (``helpers.layers_chain``): the OP pairs
of the pair table are exactly the consecutive gene pairs of ``tus.parquet`` (same threshold,
``annotation.uniop.pair_threshold`` default 0.5, as M2-12 ``operon_pairs`` used), with all three
evidence layers present.
"""

from __future__ import annotations

import itertools
from pathlib import Path

import pandas as pd
import pytest
from helpers.layers_chain import Chain, build_chain

from gippa.annotate.layers import genes_of, replicons_of, terminators_of
from gippa.config import GippaConfig

RNA_TYPES = ("rrna", "trna", "ncrna", "tmrna")


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(Path(tmp_path_factory.mktemp("pairs_chain")))


@pytest.fixture(scope="module")
def pairs(chain: Chain) -> pd.DataFrame:
    from gippa.operons.pairs import Evidence, pair_frame

    gdf = chain.frame("genes")
    feats = chain.frame("features")
    ev = Evidence(
        windows=chain.frame("promoter_windows"),
        terminators=tuple(terminators_of(chain.frame("terminators"))),
        tx=chain.frame("tx_rates"),
        rna=feats[feats["type"].isin(RNA_TYPES)].reset_index(drop=True),
    )
    tbl = pd.DataFrame(
        {"locus_tag": gdf["locus_tag"], "detected": False, "p_detect": float("nan"),
         "abundance_log2": float("nan")}
    )  # fmt: skip
    return pair_frame(
        tbl, genes_of(gdf), replicons_of(chain.frame("replicons")), chain.frame("pairs_raw"),
        ev, GippaConfig(seed=1),
    )  # fmt: skip


def test_op_pairs_reproduce_tu_membership(chain: Chain, pairs: pd.DataFrame) -> None:
    tus = chain.frame("tus")
    want = {(a, b) for genes in tus["genes"] for a, b in itertools.pairwise(genes)}
    op = pairs[pairs["class"] == "OP"]
    assert set(zip(op["gene_a"], op["gene_b"], strict=True)) == want
    assert len(want) > 0


def test_every_gene_pair_once_and_evidence_present(chain: Chain, pairs: pd.DataFrame) -> None:
    adj = pairs[pairs["class"] != "SW"]
    keys = [frozenset(k) for k in zip(adj["gene_a"], adj["gene_b"], strict=True)]
    assert len(keys) == len(set(keys))
    assert pairs["pair_id"].is_unique
    tx_ok = pairs["class"].isin(["OP", "NO", "SW"])
    for col in ("promoter_between", "tx_rate_between", "terminator_between"):
        assert pairs.loc[~tx_ok, col].isna().all()
        assert pairs.loc[tx_ok, col].notna().any(), col
