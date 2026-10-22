"""Gene-set collections from GO, GO slim, KEGG and COG (spec §7.7.1, §7.7.4; M8-03).

Pure: the resources are read by :func:`giap.io.function_resources.read_resources`. Gene
annotations come from the gene table (``cog_category``, ``go_terms``, ``kegg_ko``). GO terms are
mapped to their primary ID, propagated over ``is_a`` and ``part_of`` and split by namespace; the
three roots are removed. KEGG keeps ``path:map`` links and drops the global/overview maps
(``map011xx``, ``map012xx``, ``map01310``, ``map01320``; ADR-0023). Sets hold eligible genes
only and are never empty.
"""

from __future__ import annotations

import math
import re
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import numpy.typing as npt
import pandas as pd

B8 = npt.NDArray[np.bool_]
Collection = Literal["cog", "kegg", "go_bp", "go_mf", "go_cc", "go_slim"]
COLLECTIONS: tuple[Collection, ...] = ("cog", "kegg", "go_bp", "go_mf", "go_cc", "go_slim")
GO_ROOTS: Mapping[str, Collection] = {
    "GO:0008150": "go_bp",
    "GO:0003674": "go_mf",
    "GO:0005575": "go_cc",
}
NAMESPACES: Mapping[str, Collection] = {
    "biological_process": "go_bp",
    "molecular_function": "go_mf",
    "cellular_component": "go_cc",
}
# KEGG "Global and overview maps" (ADR-0023)
_OVERVIEW = re.compile(r"^map01(?:[12]\d\d|310|320)$")
_MISSING_NAME = "NA"


@dataclass(frozen=True, slots=True)
class GoDag:
    """A GO OBO as plain mappings over live (non-obsolete) primary IDs."""

    names: Mapping[str, str]
    namespace: Mapping[str, str]
    parents: Mapping[str, frozenset[str]]  # direct is_a + part_of parents
    alt: Mapping[str, str]  # alt_id -> primary ID


@dataclass(frozen=True, slots=True)
class Resources:
    """Installed function resources; every field is None when not installed."""

    go: GoDag | None
    slim: GoDag | None
    kegg_links: pd.DataFrame | None  # columns ko, pathway (raw link/pathway/ko rows)
    kegg_names: pd.DataFrame | None  # columns pathway, name
    cog_names: pd.DataFrame | None  # columns category, name (fun-20.tab)


@dataclass(frozen=True, slots=True)
class GeneSets:
    collection: Collection
    sets: Mapping[str, frozenset[str]]  # term -> locus tags (eligible genes only)
    names: Mapping[str, str]  # term -> human name
    parents: Mapping[str, frozenset[str]]  # GO only: term -> direct non-root parents
    coverage: float  # share of eligible genes with >= 1 term here
    n_unknown_ids: int  # annotation IDs not found in the resource


def _ids(v: Any) -> list[str]:
    if isinstance(v, list | tuple | np.ndarray):
        return [str(x) for x in v if isinstance(x, str) and x]
    return []


def _gene_sets(
    c: Collection,
    members: Mapping[str, set[str]],
    names: Mapping[str, str],
    n_eligible: int,
    *,
    parents: Mapping[str, frozenset[str]] | None = None,
    n_unknown: int = 0,
) -> GeneSets:
    sets = {t: frozenset(members[t]) for t in sorted(members) if members[t]}
    covered: set[str] = set().union(*sets.values()) if sets else set()
    cov = len(covered) / n_eligible if n_eligible else math.nan
    nm = {t: names.get(t, _MISSING_NAME) for t in sets}
    par = {t: parents.get(t, frozenset()) for t in sets} if parents is not None else {}
    return GeneSets(c, sets, nm, par, cov, n_unknown)


def _closure(dag: GoDag) -> Callable[[str], frozenset[str]]:
    memo: dict[str, frozenset[str]] = {}

    def up(t: str) -> frozenset[str]:
        if t not in memo:
            out = {t}
            for p in dag.parents.get(t, ()):
                out |= up(p)
            memo[t] = frozenset(out)
        return memo[t]

    return up


def _go(tags: list[str], terms: list[list[str]], res: Resources, n: int) -> list[GeneSets]:
    dag = res.go
    assert dag is not None
    up = _closure(dag)
    members: dict[Collection, dict[str, set[str]]] = {c: {} for c in COLLECTIONS}
    slim: frozenset[str] = frozenset()
    if res.slim is not None:
        slim = frozenset(dag.alt.get(t, t) for t in res.slim.names) - set(GO_ROOTS)
    unknown = 0
    for tag, ids in zip(tags, terms, strict=True):
        prop: set[str] = set()
        for raw in ids:
            t = dag.alt.get(raw, raw)
            if t not in dag.namespace:
                unknown += 1
                continue
            prop |= up(t)
        prop -= set(GO_ROOTS)
        for t in prop:
            c = NAMESPACES.get(dag.namespace[t])
            if c is not None:
                members[c].setdefault(t, set()).add(tag)
        for t in prop & slim:
            members["go_slim"].setdefault(t, set()).add(tag)
    seen = {t for c in ("go_bp", "go_mf", "go_cc") for t in members[c]}
    ns = dag.namespace
    parents = {
        t: frozenset(p for p in dag.parents[t] if p not in GO_ROOTS and ns.get(p) == ns[t])
        for t in seen
    }
    out = [
        _gene_sets(c, members[c], dag.names, n, parents=parents, n_unknown=unknown)
        for c in ("go_bp", "go_mf", "go_cc")
    ]
    if res.slim is not None:
        names = {**dag.names, **res.slim.names}
        out.append(_gene_sets("go_slim", members["go_slim"], names, n, n_unknown=unknown))
    return out


def kegg_pathways(links: pd.DataFrame) -> dict[str, frozenset[str]]:
    """KO -> ``mapNNNNN`` pathways from ``link/pathway/ko`` rows: ``path:map`` rows only (the
    ``path:ko`` rows duplicate them); the global/overview maps (``map011xx``, ``map012xx``,
    ``map01310``, ``map01320``; ADR-0023) are excluded."""
    out: dict[str, set[str]] = {}
    for ko, pw in zip(links["ko"].astype(str), links["pathway"].astype(str), strict=True):
        m = pw.removeprefix("path:")
        if m.startswith("map") and not _OVERVIEW.match(m):
            out.setdefault(ko.removeprefix("ko:"), set()).add(m)
    return {k: frozenset(v) for k, v in out.items()}


def _kegg(tags: list[str], kos: list[list[str]], res: Resources, n: int) -> GeneSets:
    assert res.kegg_links is not None
    path = kegg_pathways(res.kegg_links)
    known = {str(k).removeprefix("ko:") for k in res.kegg_links["ko"]}
    members: dict[str, set[str]] = {}
    unknown = 0
    for tag, ids in zip(tags, kos, strict=True):
        for ko in ids:
            ko = ko.removeprefix("KEGG:").removeprefix("ko:")
            unknown += ko not in known
            for m in path.get(ko, ()):
                members.setdefault(m, set()).add(tag)
    names: dict[str, str] = {}
    if res.kegg_names is not None:
        for k, v in zip(res.kegg_names["pathway"], res.kegg_names["name"], strict=True):
            names.setdefault(str(k).removeprefix("path:"), str(v))
    return _gene_sets("kegg", members, names, n, n_unknown=unknown)


def _cog(tags: list[str], cats: list[Any], res: Resources, n: int) -> GeneSets:
    assert res.cog_names is not None
    names = dict(zip(res.cog_names["category"].astype(str), res.cog_names["name"].astype(str),
                     strict=True))  # fmt: skip
    members: dict[str, set[str]] = {}
    unknown = 0
    for tag, cat in zip(tags, cats, strict=True):
        if not isinstance(cat, str):
            continue
        for ch in cat.strip():
            if ch in names:
                members.setdefault(ch, set()).add(tag)
            elif ch != "-":
                unknown += 1
    return _gene_sets("cog", members, names, n, n_unknown=unknown)


def _column(genes: pd.DataFrame, el: B8, col: str) -> list[Any]:
    if col not in genes:
        return [None] * int(el.sum())
    return list(genes[col].to_numpy(object)[el])


def build_collections(
    genes: pd.DataFrame, eligible: B8, res: Resources
) -> tuple[list[GeneSets], list[str]]:
    """The collections of the installed resources, in ``COLLECTIONS`` order, and run flags
    (``layer_missing:go|kegg|cog``; ``bakta_db_light`` when the gene table says so)."""
    el = np.asarray(eligible, dtype=bool)
    n = int(el.sum())
    tags = [str(t) for t in genes["locus_tag"].to_numpy(object)[el]]
    flags: list[str] = []
    out: list[GeneSets] = []
    if res.cog_names is not None:
        out.append(_cog(tags, _column(genes, el, "cog_category"), res, n))
    else:
        flags.append("layer_missing:cog")
    if res.kegg_links is not None:
        kos = [_ids(v) for v in _column(genes, el, "kegg_ko")]
        out.append(_kegg(tags, kos, res, n))
    else:
        flags.append("layer_missing:kegg")
    if res.go is not None:
        out.extend(_go(tags, [_ids(v) for v in _column(genes, el, "go_terms")], res, n))
    else:
        flags.append("layer_missing:go")
    if "annotation_source" in genes and (genes["annotation_source"] == "bakta_light").any():
        flags.append("bakta_db_light")
    return out, flags
