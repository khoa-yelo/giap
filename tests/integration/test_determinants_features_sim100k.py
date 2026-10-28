"""M8-07 conformance on the merged sim100k real path (docs/specs/M8-07, acceptance
"Conformance": every eligible gene has a row and the feature columns equal the gene-table columns
they come from).

``gippa run --executor local`` on gippasim SIM100K (tests/helpers/m4_real.py) then the M4 chain
(``m4_chain.run_m4``) give the real ``work/`` tables: C1 ``replicons``, C2 ``genes``, C5 ``tus``,
C6 ``layers``, C7/C8 ``regions`` / ``region_genes``, the gene table ``genetable/genes.parquet``
and M5-09's ``context/gene_context.parquet``. About 40 s (a whole local run), like
tests/integration/test_m4_real_path.py.
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import determinants_features_data as d

SEED = 20261014
STATS = {"permutations": 99, "bootstrap": 100, "alpha": {"identification": 0.05}}


@pytest.fixture(scope="module")
def features_mod() -> ModuleType:
    """First fixture, so a missing module fails before the minutes-long run."""
    return d.mod()


@pytest.fixture(scope="module")
def work(features_mod: ModuleType, tmp_path_factory: pytest.TempPathFactory) -> Path:
    import m4_chain as mc
    import m4_real as mr

    root = tmp_path_factory.mktemp("m8_features_sim100k")
    db = mc.test_db(root / "db")
    w = mr.gippa_run(mr.sim100k(root / "sim"), root / "OUT", db, SEED)
    mc.run_m4(w, mc.config(SEED, stats=STATS), db)
    return w


@pytest.fixture(scope="module")
def args(work: Path) -> dict[str, Any]:
    import m4_real as mr

    p = mr.real_paths(work)
    c2 = pd.read_parquet(p["genes"])
    return {
        "genes_tbl": pd.read_parquet(work / "genetable" / "genes.parquet"),
        "genes": d.genes_from_frame(c2),
        "tus": d.tus_from_frame(pd.read_parquet(p["tus"])),
        "layers": pd.read_parquet(p["layers"]),
        "regions": d.regions_arg(pd.read_parquet(p["regions"]), pd.read_parquet(p["region_genes"])),
        "gene_context": pd.read_parquet(work / "context" / "gene_context.parquet"),
        "replicons": d.replicons_from_frame(pd.read_parquet(p["replicons"])),
        "c2": c2,
    }


@pytest.fixture(scope="module")
def features(features_mod: ModuleType, args: dict[str, Any]) -> pd.DataFrame:
    a = args
    return features_mod.feature_frame(a["genes_tbl"], a["genes"], a["tus"], a["layers"],
                                      a["regions"], a["gene_context"], a["replicons"])  # fmt: skip


def test_every_eligible_gene_has_a_row(features: pd.DataFrame, args: dict[str, Any]) -> None:
    want = d.eligible_tags(args["genes_tbl"])
    assert len(want) > 50, f"only {len(want)} eligible genes: a vacuous fixture"
    assert sorted(features["locus_tag"].astype(str)) == sorted(want)
    assert list(features.columns) == list(d.TSV_COLUMNS)


@pytest.mark.parametrize(("feature", "source"),
                         d.CONFORMANCE_SOURCES)  # fmt: skip
def test_feature_columns_equal_the_gene_table(
    features: pd.DataFrame, args: dict[str, Any], feature: str, source: str
) -> None:
    tbl = args["genes_tbl"].set_index("locus_tag")
    bad = [(t, r[feature], tbl.loc[t, source]) for t, r in d.by_tag(features).items()
           if not d.same(r[feature], tbl.loc[t, source])]  # fmt: skip
    assert not bad, bad[:5]


def test_c6_and_c2_derived_columns(features: pd.DataFrame, args: dict[str, Any]) -> None:
    lay = args["layers"].set_index("locus_tag")
    c2 = args["c2"].set_index("locus_tag")
    ctx = args["gene_context"].set_index("locus_tag")
    assert "ribosomal_protein" in c2.columns, "parse_bakta writes the C2 column"
    rows = d.by_tag(features)
    for t, r in rows.items():
        assert d.same(r["n_promoters_tu_to_gene"], lay.loc[t, "n_promoter_windows_tu"]), t
        assert d.same(r["replichore_position"], ctx.loc[t, "replichore_position"]), t
        if not d.is_na(c2.loc[t, "ribosomal_protein"]):  # SIM100K plants none: all false
            assert bool(r["ribosomal_protein"]) is bool(c2.loc[t, "ribosomal_protein"]), t


def test_planted_prophage_and_bgc_reach_the_features(features: pd.DataFrame) -> None:
    """SIM100K plants a prophage and a BGC (M4-21)."""
    assert features["in_prophage"].astype(bool).any()
    assert features["in_bgc"].astype(bool).any()


def test_model_matrix_blocks(features_mod: ModuleType, features: pd.DataFrame) -> None:
    x, blocks = features_mod.model_matrix(features)
    assert len(x) == len(features) and set(blocks) == set(x.columns)
    for c, b in blocks.items():
        assert c in features_mod.BLOCKS[b], (c, b)
