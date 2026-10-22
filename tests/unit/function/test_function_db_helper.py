"""The synthetic function DB helper (``tests/helpers/function_db.py``; M8-06 Q6, lead decision 1).

The helper is tooling for M8-06 examples 1 and 2: on the simulator's own GO/KEGG/COG IDs it must
give parseable resources whose collections are non-empty, with GO parents of >= 10 genes, and it
must carry no KEGG content (made-up ``map9xxxx`` pathways only).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from helpers.function_db import (
    COG_FUN,
    ITEMS,
    PATHS,
    go_ancestors,
    ko_pathway,
    write_function_db,
)


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory):  # type: ignore[no-untyped-def]
    import giapsim.layout as la

    lay = la.simulate_layout(la.LayoutConfig(n_genes_chromosome=400), np.random.default_rng(11))
    root = tmp_path_factory.mktemp("fdb")
    return lay, write_function_db(root / "db", lay.xrefs)


def _genes(lay, db_root: Path) -> pd.DataFrame:  # type: ignore[no-untyped-def]
    from giap.genetable.assemble import function_columns
    from giap.io.genes_tsv import read_cog_def

    tags = lay.genes["locus_tag"].astype(str).tolist()
    fn = function_columns(lay.xrefs, read_cog_def(db_root / PATHS["cog2020_def"]))
    fn = fn.set_index("locus_tag").reindex(tags)
    return pd.DataFrame({
        "locus_tag": tags,
        "partial": [False] * len(tags),
        "p_detect": [0.5] * len(tags),
        "cog_category": [v if isinstance(v, str) else None for v in fn["cog_category"]],
        "go_terms": [v if isinstance(v, list) else [] for v in fn["go_terms"]],
        "kegg_ko": [v if isinstance(v, list) else [] for v in fn["kegg_ko"]],
        "annotation_source": ["bakta_full"] * len(tags),
    })  # fmt: skip


def test_every_item_written_at_its_db_sources_path(sim) -> None:  # type: ignore[no-untyped-def]
    from giap.setupdb.items import SOURCES

    _, db = sim
    for name in ITEMS:
        assert db.paths[name] == db.root / SOURCES[name].path, name
        assert db.paths[name].is_file(), name
    man = json.loads((db.root / "db_manifest.json").read_text(encoding="utf-8"))
    assert sorted(i["name"] for i in man["items"] if i["status"] == "ok") == sorted(ITEMS)


def test_manifest_read_by_preflight(sim) -> None:  # type: ignore[no-untyped-def]
    from giap.preflight.environment import read_manifest

    _, db = sim
    found = read_manifest(str(db.root)) or {}
    assert set(ITEMS) <= set(found)


def test_resources_parse_and_collections_are_populated(sim) -> None:  # type: ignore[no-untyped-def]
    from giap.function.collections import COLLECTIONS, build_collections
    from giap.genetable.eligibility import eligible
    from giap.io.function_resources import read_resources

    lay, db = sim
    genes = _genes(lay, db.root)
    colls, flags = build_collections(genes, eligible(genes), read_resources(db.inputs()))
    assert flags == []
    assert [c.collection for c in colls] == list(COLLECTIONS)
    for c in colls:
        assert c.sets, c.collection
        assert c.n_unknown_ids == 0, c.collection  # every sim ID is in the DB
    for c in colls:
        if c.collection in ("go_bp", "go_mf", "go_cc", "go_slim"):
            assert max(len(s) for s in c.sets.values()) >= 10, c.collection


def test_go_hierarchy_parent_and_group(sim) -> None:  # type: ignore[no-untyped-def]
    from giap.io.function_resources import read_go

    lay, db = sim
    dag = read_go(db.paths["go_basic"])
    ids = sorted(set(lay.xrefs.loc[lay.xrefs["db"] == "GO", "id"]))
    assert ids
    for t in ids[:50]:
        parent, group = go_ancestors(t)
        assert dag.parents[t] == frozenset({parent})
        assert dag.parents[parent] == frozenset({group})
        assert dag.namespace[t] == dag.namespace[parent] == dag.namespace[group]


def test_cog_fun_has_the_ncbi_categories(sim) -> None:  # type: ignore[no-untyped-def]
    from giap.io.function_resources import read_resources

    _, db = sim
    res = read_resources(db.inputs(["cog2020_fun"]))
    assert res.cog_names is not None
    names = dict(zip(res.cog_names["category"], res.cog_names["name"], strict=True))
    assert names == COG_FUN
    assert names["J"] == "Translation, ribosomal structure and biogenesis"


def test_kegg_is_synthetic_only(sim) -> None:  # type: ignore[no-untyped-def]
    """Lead decision 1 / index decision 4: no real KEGG pathway IDs or names."""
    lay, db = sim
    rows = [line.split("\t") for line in db.paths["kegg_ko_pathway"].read_text().splitlines()]
    assert rows
    maps = {p.removeprefix("path:") for _, p in rows}
    assert all(m.startswith(("map9", "ko9")) for m in maps), maps
    kos = sorted(set(lay.xrefs.loc[lay.xrefs["db"] == "KEGG", "id"]))
    assert {k.removeprefix("ko:") for k, _ in rows} == set(kos)
    for ko in kos:
        assert [f"ko:{ko}", f"path:{ko_pathway(ko)}"] in rows
    names = db.paths["kegg_pathway_names"].read_text().splitlines()
    assert all(n.startswith("map9") and "Synthetic pathway" in n for n in names)


def test_deterministic(sim, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    lay, db = sim
    again = write_function_db(tmp_path / "db", lay.xrefs)
    for name in ITEMS:
        assert again.paths[name].read_bytes() == db.paths[name].read_bytes(), name
