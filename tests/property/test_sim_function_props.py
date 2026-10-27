"""Hypothesis properties of gippasim's synthetic ontology and function worlds (issue M8-14;
docs/specs/M8-14-sim-function-determinants.md §1, lead decision 1; interface as assumed in
tests/unit/test_sim_function.py).

For any term count, depth, gene list and seed, ``synthetic_ontology``: parses with M8-03's reader
(goatools); has exactly ``n_terms`` live terms, all ``GO:9xxxxxx``; is acyclic with every parent
defined; annotates only genes of the list with terms of the OBO; is a pure function of the seed.
For any small ORA/GSEA world: every GO, KO and pathway ID is in its reserved range.
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

SETTINGS = settings(max_examples=25, deadline=None,
                    suppress_health_check=[HealthCheck.too_slow])  # fmt: skip
RE_GO = re.compile(r"^GO:9\d{6}$")
RE_KO = re.compile(r"^K9\d{4}$")
RE_MAP = re.compile(r"^map9\d{4}$")


def _fs() -> Any:
    import gippasim.function_sim

    return gippasim.function_sim


def _read_go(obo: str) -> Any:
    from gippa.io.function_resources import read_go

    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "s.obo"
        path.write_text(obo)
        return read_go(path)


@SETTINGS
@given(
    n_terms=st.integers(3, 80),
    depth=st.integers(1, 6),
    n_genes=st.integers(0, 300),
    seed=st.integers(0, 2**32 - 1),
)
def test_synthetic_ontology_invariants(n_terms: int, depth: int, n_genes: int, seed: int) -> None:
    genes = [f"G{i:04d}" for i in range(n_genes)]
    obo, table = _fs().synthetic_ontology(n_terms, depth, np.random.default_rng(seed),
                                          genes=genes)  # fmt: skip
    dag = _read_go(obo)
    assert len(dag.names) == n_terms
    assert all(RE_GO.match(t) for t in dag.names)
    assert not dag.alt or all(RE_GO.match(a) for a in dag.alt)
    for t, ps in dag.parents.items():
        assert t not in ps
        assert ps <= set(dag.names)
    state: dict[str, int] = {}

    def visit(t: str) -> None:  # 1 = on stack, 2 = done
        state[t] = 1
        for p in dag.parents[t]:
            assert state.get(p) != 1, f"cycle through {p}"
            if p not in state:
                visit(p)
        state[t] = 2

    for t in dag.names:
        if t not in state:
            visit(t)
    assert list(table.columns) == ["locus_tag", "db", "id"]
    assert set(table["id"]) <= set(dag.names)
    assert set(table["locus_tag"]) <= set(genes)
    assert not table.duplicated().any()
    again = _fs().synthetic_ontology(n_terms, depth, np.random.default_rng(seed), genes=genes)
    assert again[0] == obo
    pd.testing.assert_frame_equal(again[1], table)


@SETTINGS
@given(
    which=st.sampled_from(["ora_null", "ora_attack", "gsea_null", "gsea_packed"]),
    n_genes=st.integers(400, 900),
    seed=st.integers(0, 2**32 - 1),
)
def test_world_ids_always_reserved(which: str, n_genes: int, seed: int) -> None:
    import gippasim.scenarios as sc

    pair = {
        "ora_null": lambda: sc.ora_null(n_genes),
        "ora_attack": lambda: sc.ora_attack_confounded(n_genes, -1.5),
        "gsea_null": lambda: sc.gsea_null(n_genes),
        "gsea_packed": lambda: sc.gsea_attack_operon_packed(n_genes, 3),
    }[which]()
    w = _fs().draw(pair, np.random.default_rng(seed))
    x = w.layout.xrefs
    assert all(RE_GO.match(t) for t in x.loc[x["db"] == "GO", "id"])
    assert all(RE_KO.match(t) for t in x.loc[x["db"] == "KEGG", "id"])
    assert all(RE_GO.match(t) for t in re.findall(r"^id:\s*(\S+)", w.obo, flags=re.M))
    assert all(RE_KO.match(str(k).removeprefix("ko:")) for k in w.kegg_links["ko"])
    assert all(RE_MAP.match(str(p).removeprefix("path:")) for p in w.kegg_links["pathway"])
    assert set(w.sets["locus_tag"]) <= set(w.layout.genes["locus_tag"])
