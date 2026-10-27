"""A synthetic ``$GIPPA_DB_DIR`` with the five M8 function resources (M8-06 Q6, lead decision 1).

Built from the IDs a run actually carries (the simulator's own GO, KEGG and COG IDs, or any C3
``xrefs`` frame), so the collections are non-empty. Nothing here is real GO or KEGG content: the
KEGG pathways are made-up ``map9xxxx`` IDs with made-up names (M8 index lead decision 4: no KEGG
file or KEGG-derived table is ever committed or generated). The COG 2020 functional categories in
``fun-20.tab`` are the NCBI category letters and names (public domain, as M8-06 example 4 needs
"Translation, ribosomal structure and biogenesis").

Layout (``db_sources.yaml`` paths; ``db_manifest.json`` is a ``gippa setup`` manifest)::

    go/go-basic.obo          the three roots; per namespace 2 group terms and 4 parent terms;
                             every input GO ID is a leaf under one parent (ID mod 3 -> namespace,
                             ID mod 4 -> parent; parents 0-1 under group 0, 2-3 under group 1)
    go/goslim_generic.obo    the roots and the 12 parent terms
    kegg/ko_pathway.tsv      ko:Knnnnn -> path:map9000{1..5} (KO number mod 5), + path:ko rows
    kegg/pathway_names.tsv   map9000{1..5} -> "Synthetic pathway <i>"
    cog/fun-20.tab           the 26 COG 2020 functional categories
    cog/cog-20.def.tab       every input COG ID -> ``COG_LETTERS[number mod 9]``

Nothing here imports ``gippa``.
"""

from __future__ import annotations

import json
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

FUNCTION_RESOURCES = ("go_basic", "goslim_generic", "kegg_ko_pathway", "kegg_pathway_names",
                      "cog2020_fun")  # fmt: skip
ITEMS = (*FUNCTION_RESOURCES, "cog2020_def")
PATHS = {
    "go_basic": "go/go-basic.obo",
    "goslim_generic": "go/goslim_generic.obo",
    "kegg_ko_pathway": "kegg/ko_pathway.tsv",
    "kegg_pathway_names": "kegg/pathway_names.tsv",
    "cog2020_fun": "cog/fun-20.tab",
    "cog2020_def": "cog/cog-20.def.tab",
}

ROOTS = {"biological_process": "GO:0008150", "molecular_function": "GO:0003674",
         "cellular_component": "GO:0005575"}  # fmt: skip
NAMESPACES = tuple(ROOTS)
COLLECTION = {"biological_process": "go_bp", "molecular_function": "go_mf",
              "cellular_component": "go_cc"}  # fmt: skip

#: NCBI COG 2020 functional categories (fun-20.tab: letter, colour, description)
COG_FUN: dict[str, str] = {
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "X": "Mobilome: prophages, transposons",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
}
#: cog-20.def.tab category of COG<n> is COG_LETTERS[n % 9] (one two-letter category)
COG_LETTERS = ("J", "K", "L", "M", "E", "C", "G", "S", "KT")
N_PATHWAYS = 5
N_PARENTS = 4
N_GROUPS = 2


def group_id(ns: int, j: int) -> str:
    return f"GO:{990100 + 10 * ns + j:07d}"


def parent_id(ns: int, j: int) -> str:
    return f"GO:{990000 + 10 * ns + j:07d}"


def pathway_id(i: int) -> str:
    return f"map9{i:04d}"


def _num(s: str) -> int:
    m = re.search(r"(\d+)$", s)
    if m is None:
        raise ValueError(f"no number in ID {s!r}")
    return int(m.group(1))


def go_namespace(go_id: str) -> str:
    """The namespace a leaf GO ID is given (ID number mod 3)."""
    return NAMESPACES[_num(go_id) % 3]


def go_parent(go_id: str) -> str:
    """The parent term of a leaf GO ID (ID number mod 4, within its namespace)."""
    return parent_id(NAMESPACES.index(go_namespace(go_id)), _num(go_id) % N_PARENTS)


def go_ancestors(go_id: str) -> tuple[str, str]:
    """(parent, group) of a leaf GO ID; the root is not listed."""
    ns = NAMESPACES.index(go_namespace(go_id))
    j = _num(go_id) % N_PARENTS
    return parent_id(ns, j), group_id(ns, j // 2)


def ko_pathway(ko: str) -> str:
    return pathway_id(_num(ko) % N_PATHWAYS + 1)


def cog_letters(cog_id: str) -> str:
    return COG_LETTERS[_num(cog_id) % len(COG_LETTERS)]


def _block(t: str, name: str, ns: str, isa: Iterable[str], subset: str | None = None) -> str:
    lines = ["[Term]", f"id: {t}", f"name: {name}", f"namespace: {ns}",
             f'def: "Synthetic test term {t}." [GOC:gippa]']  # fmt: skip
    if subset:
        lines.append(f"subset: {subset}")
    lines += [f"is_a: {p}" for p in isa]
    return "\n".join(lines) + "\n"


_HEADER = "format-version: 1.2\ndata-version: releases/2026-05-19\nontology: go\n"
_TYPEDEF = "[Typedef]\nid: part_of\nname: part of\nxref: BFO:0000050\nis_transitive: true\n"


def _upper_terms(subset: str | None, with_groups: bool) -> list[str]:
    blocks = [_block(r, ns, ns, (), subset) for ns, r in ROOTS.items()]
    for i, ns in enumerate(NAMESPACES):
        for g in range(N_GROUPS):
            if with_groups:
                blocks.append(_block(group_id(i, g), f"synthetic {ns} group {g}", ns,
                                     (ROOTS[ns],)))  # fmt: skip
        for j in range(N_PARENTS):
            up = group_id(i, j // 2) if with_groups else ROOTS[ns]
            blocks.append(_block(parent_id(i, j), f"synthetic {ns} parent {j}", ns, (up,),
                                 subset))  # fmt: skip
    return blocks


def go_obo_text(go_ids: Iterable[str]) -> str:
    blocks = _upper_terms(None, with_groups=True)
    for t in sorted(set(go_ids)):
        blocks.append(_block(t, f"synthetic leaf {t}", go_namespace(t), (go_parent(t),)))
    return _HEADER + "\n" + "\n".join(blocks) + "\n" + _TYPEDEF


def slim_obo_text() -> str:
    head = _HEADER + 'subsetdef: goslim_generic "Generic GO slim"\n'
    return (
        head
        + "\n"
        + "\n".join(_upper_terms("goslim_generic", with_groups=False))
        + "\n"
        + (_TYPEDEF)
    )


def kegg_links_text(kos: Iterable[str]) -> str:
    rows = []
    for ko in sorted(set(kos)):
        m = ko_pathway(ko)
        rows += [f"ko:{ko}\tpath:{m}", f"ko:{ko}\tpath:ko{m[3:]}"]
    return "".join(r + "\n" for r in rows)


def kegg_names_text() -> str:
    return "".join(f"{pathway_id(i)}\tSynthetic pathway {i}\n" for i in range(1, N_PATHWAYS + 1))


def cog_fun_text() -> str:
    return "".join(f"{k}\tE0E0E0\t{v}\n" for k, v in COG_FUN.items())


def cog_def_text(cog_ids: Iterable[str]) -> str:
    return "".join(
        f"{c}\t{cog_letters(c)}\tSynthetic {c}\tSYN\t\t\t\n" for c in sorted(set(cog_ids))
    )


def manifest_doc(items: Iterable[str]) -> dict[str, object]:
    """A ``gippa setup`` manifest listing ``items`` as ``ok``."""
    rows = [{"name": n, "version": "synthetic", "url": "", "checksum": "", "size_bytes": 0,
             "retrieved": "", "path": PATHS[n], "status": "ok", "message": ""}
            for n in items]  # fmt: skip
    return {"schema_version": "1.0", "items": rows}


@dataclass(frozen=True)
class FunctionDb:
    root: Path
    paths: Mapping[str, Path | None]  # every ITEMS name; None when not written

    def inputs(self, items: Iterable[str] = FUNCTION_RESOURCES) -> dict[str, Path]:
        """The ``function`` step's optional inputs for the written ``items``."""
        return {k: p for k in items if (p := self.paths.get(k)) is not None}


def ids_of(xrefs: pd.DataFrame, db: str) -> list[str]:
    """Distinct IDs of database ``db`` (C3 ``xrefs``: ``locus_tag``, ``db``, ``id``)."""
    ids = xrefs.loc[xrefs["db"].astype(str) == db, "id"].astype(str)
    if db == "COG":
        ids = ids[ids.str.match(r"^COG\d+$")]
    return sorted(set(ids))


def write_function_db(
    root: Path,
    xrefs: pd.DataFrame,
    *,
    items: Iterable[str] = ITEMS,
    listed: Iterable[str] | None = None,
) -> FunctionDb:
    """Write the synthetic DB under ``root`` for the IDs in ``xrefs``; ``items`` are the files
    written, ``listed`` the manifest entries (default: ``items``)."""
    write = set(items)
    texts = {
        "go_basic": lambda: go_obo_text(ids_of(xrefs, "GO")),
        "goslim_generic": slim_obo_text,
        "kegg_ko_pathway": lambda: kegg_links_text(ids_of(xrefs, "KEGG")),
        "kegg_pathway_names": kegg_names_text,
        "cog2020_fun": cog_fun_text,
        "cog2020_def": lambda: cog_def_text(ids_of(xrefs, "COG")),
    }
    root.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path | None] = {}
    for name in ITEMS:
        if name not in write:
            paths[name] = None
            continue
        p = root / PATHS[name]
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(texts[name](), encoding="utf-8")
        paths[name] = p
    names = [n for n in ITEMS if n in (write if listed is None else set(listed))]
    (root / "db_manifest.json").write_text(json.dumps(manifest_doc(names), indent=2) + "\n",
                                           encoding="utf-8")  # fmt: skip
    return FunctionDb(root, paths)
