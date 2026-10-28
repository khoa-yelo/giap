"""Synthetic function resources for the M8-03 tests (docs/specs/M8-03-function-collections.md).

Everything here is written by the tests themselves: a toy ``go-basic.obo`` (is_a and part_of
edges, an alt ID, an obsolete term), a toy ``goslim_generic.obo``, a two-pathway KEGG KO ->
pathway link table in the ``rest.kegg.jp/link/pathway/ko`` format (``path:map`` rows, their
``path:ko`` duplicates and global/overview maps) with a ``list/pathway`` name table, and a toy COG
2020 ``fun-20.tab`` / ``cog-20.def.tab``. No real GO, KEGG or COG content (M8 index, lead
decision 4: no KEGG file or KEGG-derived table is ever committed).

Toy GO DAG (``<-`` is ``is_a``, ``<~`` is ``part_of``)::

    BP root GO:0008150 <- B1 GO:0000101 <- B2 GO:0000102 (alt_id GO:0000112)
    BP root            <- B3 GO:0000103,  B3 <~ B1
    B2 <- B4 GO:0000104,  B3 <- B4        (diamond)
    BP root            <- B5 GO:0000105,  B5 <~ B4
    MF root GO:0003674 <- M1 GO:0000201 <- M2 GO:0000202
    CC root GO:0005575 <- C1 GO:0000301,  CC root <- C2 GO:0000302,  C2 <~ C1
    GO:0000199 obsolete (BP)
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

BP_ROOT, MF_ROOT, CC_ROOT = "GO:0008150", "GO:0003674", "GO:0005575"
ROOTS = frozenset({BP_ROOT, MF_ROOT, CC_ROOT})
B1, B2, B3, B4, B5 = (f"GO:000010{i}" for i in range(1, 6))
B2_ALT = "GO:0000112"
M1, M2 = "GO:0000201", "GO:0000202"
C1, C2 = "GO:0000301", "GO:0000302"
OBSOLETE = "GO:0000199"
ABSENT = "GO:9999999"

NAMESPACE = {BP_ROOT: "biological_process", MF_ROOT: "molecular_function",
             CC_ROOT: "cellular_component"}  # fmt: skip

#: term -> (name, namespace, is_a parents, part_of parents)
TERMS: dict[str, tuple[str, str, tuple[str, ...], tuple[str, ...]]] = {
    BP_ROOT: ("biological_process", "biological_process", (), ()),
    MF_ROOT: ("molecular_function", "molecular_function", (), ()),
    CC_ROOT: ("cellular_component", "cellular_component", (), ()),
    B1: ("toy process one", "biological_process", (BP_ROOT,), ()),
    B2: ("toy process two", "biological_process", (B1,), ()),
    B3: ("toy process three", "biological_process", (BP_ROOT,), (B1,)),
    B4: ("toy process four", "biological_process", (B2, B3), ()),
    B5: ("toy process five", "biological_process", (BP_ROOT,), (B4,)),
    M1: ("toy function one", "molecular_function", (MF_ROOT,), ()),
    M2: ("toy function two", "molecular_function", (M1,), ()),
    C1: ("toy component one", "cellular_component", (CC_ROOT,), ()),
    C2: ("toy component two", "cellular_component", (CC_ROOT,), (C1,)),
}
ALT_IDS = {B2: (B2_ALT,)}

#: every live term -> its ancestors over is_a and part_of (excluding itself)
ANCESTORS: dict[str, frozenset[str]] = {}


def _anc(t: str) -> frozenset[str]:
    _, _, isa, part = TERMS[t]
    out: set[str] = set()
    for p in (*isa, *part):
        out |= {p, *_anc(p)}
    return frozenset(out)


for _t in TERMS:
    ANCESTORS[_t] = _anc(_t)

#: direct is_a / part_of parents (what ``GeneSets.parents`` must hold)
DIRECT_PARENTS = {t: frozenset((*v[2], *v[3])) for t, v in TERMS.items()}
GO_COLLECTION = {"biological_process": "go_bp", "molecular_function": "go_mf",
                 "cellular_component": "go_cc"}  # fmt: skip

#: slim terms (the roots plus three toy terms)
SLIM_TERMS = (BP_ROOT, MF_ROOT, CC_ROOT, B1, B3, M1)


def _term_block(t: str, name: str, ns: str, isa: Sequence[str], part: Sequence[str],
                alt: Sequence[str] = (), subset: str | None = None) -> str:  # fmt: skip
    lines = ["[Term]", f"id: {t}", f"name: {name}", f"namespace: {ns}"]
    lines += [f"alt_id: {a}" for a in alt]
    lines.append(f'def: "Synthetic test term {t}." [GOC:gippa]')
    if subset:
        lines.append(f"subset: {subset}")
    lines += [f"is_a: {p} ! {TERMS[p][0]}" for p in isa]
    lines += [f"relationship: part_of {p} ! {TERMS[p][0]}" for p in part]
    return "\n".join(lines) + "\n"


_HEADER = "format-version: 1.2\ndata-version: releases/2026-05-19\nontology: go\n\n"
_TYPEDEF = "[Typedef]\nid: part_of\nname: part of\nxref: BFO:0000050\nis_transitive: true\n"


def go_obo_text() -> str:
    blocks = [_term_block(t, *v, alt=ALT_IDS.get(t, ())) for t, v in TERMS.items()]
    blocks.append(
        "[Term]\n"
        f"id: {OBSOLETE}\n"
        "name: obsolete toy process\n"
        "namespace: biological_process\n"
        'def: "OBSOLETE. Synthetic." [GOC:gippa]\n'
        "is_obsolete: true\n"
        f"replaced_by: {B1}\n"
    )
    return _HEADER + "\n".join(blocks) + "\n" + _TYPEDEF


def slim_obo_text() -> str:
    """The slim: its own terms with is_a / part_of among slim terms only (as goslim files do)."""
    keep = set(SLIM_TERMS)
    blocks = []
    for t in SLIM_TERMS:
        name, ns, isa, part = TERMS[t]
        isa_k = [p for p in isa if p in keep]
        part_k = [p for p in part if p in keep]
        blocks.append(_term_block(t, name, ns, isa_k, part_k, subset="goslim_generic"))
    return (
        _HEADER
        + 'subsetdef: goslim_generic "Generic GO slim"\n\n'
        + "\n".join(blocks)
        + ("\n" + _TYPEDEF)
    )


#: synthetic KEGG: two real pathways plus global/overview maps (both ``path:map`` and ``path:ko``)
KEGG_PATHWAYS = ("map00010", "map00020")
KEGG_OVERVIEW = ("map01100", "map01110", "map01200", "map01230")
KEGG_LINKS: tuple[tuple[str, str], ...] = (
    ("K00001", "map00010"),
    ("K00001", "map01100"),
    ("K00002", "map00010"),
    ("K00002", "map00020"),
    ("K00002", "map01110"),
    ("K00003", "map00020"),
    ("K00003", "map01200"),
    ("K00004", "map01230"),  # a KO only in an overview map
)
#: map00020 has no name row (degenerate case: name NA, kept)
KEGG_NAMES = {
    "map00010": "Toy pathway alpha",
    "map01100": "Toy global map",
    "map01110": "Toy overview one",
    "map01200": "Toy overview two",
    "map01230": "Toy overview three",
}


def kegg_links_text() -> str:
    rows = []
    for ko, m in KEGG_LINKS:  # same layout as link/pathway/ko: both rows per link
        rows.append(f"ko:{ko}\tpath:{m}")
        rows.append(f"ko:{ko}\tpath:ko{m[3:]}")
    return "\n".join(rows) + "\n"


def kegg_names_text() -> str:
    return "".join(f"{k}\t{v}\n" for k, v in KEGG_NAMES.items())


COG_FUN = {
    "J": "Toy translation",
    "K": "Toy transcription",
    "L": "Toy replication",
    "M": "Toy envelope",
    "R": "Toy general function prediction only",
    "S": "Toy function unknown",
}


def cog_fun_text() -> str:
    return "".join(f"{k}\tFFFFFF\t{v}\n" for k, v in COG_FUN.items())


def cog_def_text() -> str:
    return "COG9001\tKL\tToy COG\tTOY\t\t\t\nCOG9002\tJ\tToy COG two\tTOY2\t\t\t\n"


ITEM_FILES = {
    "go_basic": ("go/go-basic.obo", go_obo_text),
    "goslim_generic": ("go/goslim_generic.obo", slim_obo_text),
    "kegg_ko_pathway": ("kegg/ko_pathway.tsv", kegg_links_text),
    "kegg_pathway_names": ("kegg/pathway_names.tsv", kegg_names_text),
    "cog2020_def": ("cog/cog-20.def.tab", cog_def_text),
    "cog2020_fun": ("cog/fun-20.tab", cog_fun_text),
}


def write_resources(root: Path, omit: Iterable[str] = ()) -> dict[str, Path | None]:
    """Write the toy items under ``root`` (db_sources.yaml layout); ``omit`` items map to None."""
    skip = set(omit)
    out: dict[str, Path | None] = {}
    for item, (rel, text) in ITEM_FILES.items():
        if item in skip:
            out[item] = None
            continue
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text(), encoding="utf-8")
        out[item] = p
    return out


def genes_frame(rows: Sequence[Mapping[str, Any]], *, function_columns: bool = True,
                annotation_source: str = "bakta_full") -> pd.DataFrame:  # fmt: skip
    """A minimal gene table: ``locus_tag``, ``partial``, ``p_detect``, the function columns."""
    data: dict[str, list[Any]] = {"locus_tag": [], "partial": [], "p_detect": []}
    fn: dict[str, list[Any]] = {"cog_id": [], "cog_category": [], "go_terms": [], "kegg_ko": [],
                                "ec": []}  # fmt: skip
    for r in rows:
        data["locus_tag"].append(r["locus_tag"])
        data["partial"].append(bool(r.get("partial", False)))
        data["p_detect"].append(float(r.get("p_detect", 0.5)))
        fn["cog_id"].append(r.get("cog_id"))
        fn["cog_category"].append(r.get("cog_category"))
        fn["go_terms"].append(list(r.get("go_terms", [])))
        fn["kegg_ko"].append(list(r.get("kegg_ko", [])))
        fn["ec"].append(list(r.get("ec", [])))
    df = pd.DataFrame(data)
    if function_columns:
        for k, v in fn.items():
            df[k] = pd.Series(v, dtype=object)
        df["annotation_source"] = pd.Series([annotation_source] * len(df), dtype=object)
    return df


def eligible_mask(genes: pd.DataFrame) -> np.ndarray:
    """§7.3 eligibility, computed independently: non-partial with finite ``p_detect``."""
    p = pd.to_numeric(genes["p_detect"], errors="coerce").to_numpy(float)
    part = genes["partial"].astype("boolean").fillna(False).to_numpy(bool)
    return np.asarray(~part & np.isfinite(p), dtype=bool)


def propagated(terms: Iterable[str]) -> set[str]:
    """Oracle: live terms (alt IDs mapped) plus all is_a/part_of ancestors, roots removed."""
    out: set[str] = set()
    for t in terms:
        t = B2 if t == B2_ALT else t
        if t in TERMS:
            out |= {t, *ANCESTORS[t]}
    return out - ROOTS
