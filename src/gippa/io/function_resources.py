"""Readers for the M8 function resources (``db_sources.yaml`` items; spec §7.7.1, M8-03).

``go_basic`` / ``goslim_generic`` (OBO, read with goatools over ``is_a`` and ``part_of``),
``kegg_ko_pathway`` / ``kegg_pathway_names`` (KEGG REST ``link/pathway/ko`` and ``list/pathway``)
and ``cog2020_fun`` (COG 2020 ``fun-20.tab``). An item mapped to None is not installed.
"""

from __future__ import annotations

import csv
from collections.abc import Mapping, Sequence
from pathlib import Path

import pandas as pd
from goatools.obo_parser import GODag

from gippa.function.collections import GoDag, Resources


def read_go(path: Path) -> GoDag:
    """Live terms of an OBO file; alt IDs map to their primary, obsolete terms are dropped."""
    dag = GODag(str(path), optional_attrs={"relationship"}, load_obsolete=False, prt=None)
    names: dict[str, str] = {}
    namespace: dict[str, str] = {}
    parents: dict[str, frozenset[str]] = {}
    alt: dict[str, str] = {}
    for key, term in dag.items():
        if key != term.id:
            alt[key] = term.id
            continue
        names[key] = str(term.name)
        namespace[key] = str(term.namespace)
        part_of = getattr(term, "relationship", {}).get("part_of", ())
        parents[key] = frozenset(p.id for p in (*term.parents, *part_of))
    return GoDag(names, namespace, parents, alt)


def _tsv(path: Path, names: Sequence[str], usecols: Sequence[int]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, usecols=list(usecols), dtype=str,
                     keep_default_na=False, quoting=csv.QUOTE_NONE,
                     encoding_errors="replace")  # fmt: skip
    df.columns = list(names)
    return df


def read_resources(paths: Mapping[str, Path | None]) -> Resources:
    """Resources keyed by ``db_sources.yaml`` item name; absent or None items give None."""

    def get(item: str) -> Path | None:
        p = paths.get(item)
        return Path(p) if p is not None else None

    go, slim = get("go_basic"), get("goslim_generic")
    links, kegg_names, fun = get("kegg_ko_pathway"), get("kegg_pathway_names"), get("cog2020_fun")
    return Resources(
        go=read_go(go) if go is not None else None,
        slim=read_go(slim) if slim is not None else None,
        kegg_links=_tsv(links, ("ko", "pathway"), (0, 1)) if links is not None else None,
        kegg_names=_tsv(kegg_names, ("pathway", "name"), (0, 1)) if kegg_names else None,
        cog_names=_tsv(fun, ("category", "name"), (0, 2)) if fun is not None else None,
    )
