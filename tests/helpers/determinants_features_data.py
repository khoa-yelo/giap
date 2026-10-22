"""Inputs for the M8-07 feature-matrix tests (docs/specs/M8-07-determinants-features.md).

``toy()`` is a hand-made 13-gene genome whose genes carry the spec's examples:

* ``chr1`` circular, 10,000 bp, chromosome; ``p1`` linear, 3,000 bp, plasmid (no oriC);
* TU_1 ``+`` G01 [1000, 1300), G02 [1450, 1750): example 5 (distance 450);
* TU_2 ``-`` G03 [5000, 5600), G04 [4300, 4950): example 6 (distance 650);
* TU_3 ``+`` across the origin, G05 [9800, 9950), G06 [200, 500): example 7 (distance 400);
* G07 a ribosomal protein (``rplB``, "50S ribosomal protein L2", flag ``cai_reference``);
* G08 in ``prophage_1``; G09 a pseudogene (empty ``aa``) in ``bgc_1``;
* G10 partial and G11 without ``p_detect``: not eligible (§7.3), so no row;
* G12 the plasmid gene (``replichore_position`` NA); G13 a plain single-gene TU.

``aa`` carries examples 1-4 and the degenerate rows (G05 = "M" -> NA, G09 "" -> NA); start codons
ATG / GTG (G04) / TTG (G05) / CTG (G06 -> ``other``); ``tir_ostir`` = 0 on G02 (example 10) and
< 0 on G03. Every other feature source holds a distinct value per cell, so a column copied from
the wrong source or the wrong row shows up.

Argument conventions resolved by the test author (the spec's interface leaves them open):

* ``replicons`` is ``Mapping[str, Replicon]`` (as ``giap.annotate.operons``), C1 order = key order;
* ``regions`` is the C8 ``region_genes`` frame with the C7 ``region_type`` of each region joined
  on (``region_id``, ``locus_tag``, ``overlap_frac``, ``subset``, ``core_biosynthetic``,
  ``region_type``); region ids follow M2-19's ``<type>_<n>``;
* ``gene_context`` is M5-09's ``context/gene_context.parquet`` (``GENE_CONTEXT_COLUMNS``).
"""

from __future__ import annotations

import copy
import importlib
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd

NAN = math.nan
MODULE = "giap.determinants.features"

IDENTITY = (
    "locus_tag", "tu_id", "detected", "abundance_log2", "abundance_imputed", "ribosomal_protein",
)  # fmt: skip
TRANSCRIPTION = (
    "promotech_max", "promotech_tu", "tx_rate_max_tu", "tx_rate_sum_tu", "n_promoters_tu_to_gene",
    "tu_position", "dist_to_tu_start_bp",
)  # fmt: skip
TRANSLATION = ("tir_ostir", "sd_dg", "spacer_bp", "start_codon", "mrna_dg")
CODON = ("cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "rare_codons_2_10")
PROTEIN = ("length_aa", "mw", "pi", "gravy", "nterm_residue", "tm_segments")
CONTEXT = ("replicon_type", "replichore_position", "in_prophage", "in_bgc", "gc_local_dev")
FEATURES = (*TRANSCRIPTION, *TRANSLATION, *CODON, *PROTEIN, *CONTEXT)
#: the TSV: identity columns, then the §7.8.1 features in table order
TSV_COLUMNS = (*IDENTITY, *FEATURES)

BLOCK_NAMES = ("transcription", "translation", "codon", "protein", "context")
RESIDUES = "ACDEFGHIKLMNPQRSTVWY"
START_LEVELS = ("ATG", "GTG", "TTG", "other")
REPLICON_LEVELS = ("chromosome", "plasmid", "unknown")
#: model columns per block after encoding (one-hot names ``<feature>_<level>``)
EXPECTED_BLOCKS: dict[str, frozenset[str]] = {
    "transcription": frozenset(TRANSCRIPTION),
    "translation": frozenset(
        {"tir_ostir", "sd_dg", "spacer_bp", "mrna_dg"} | {f"start_codon_{x}" for x in START_LEVELS}
    ),
    "codon": frozenset(CODON),
    "protein": frozenset(
        {"length_aa", "mw", "pi", "gravy", "tm_segments"} | {f"nterm_residue_{r}" for r in RESIDUES}
    ),
    "context": frozenset(
        {"replichore_position", "in_prophage", "in_bgc", "gc_local_dev"}
        | {f"replicon_type_{x}" for x in REPLICON_LEVELS}
    ),
}

#: feature -> gene-table column it is copied from (spec step 2; ``n_promoters_tu_to_gene`` comes
#: from C6 and is checked separately)
GENE_TABLE_SOURCES: dict[str, str] = {
    "promotech_max": "promoter_score", "promotech_tu": "promoter_score_tu",
    "tx_rate_max_tu": "tx_rate_max_tu", "tx_rate_sum_tu": "tx_rate_sum_tu",
    "tu_position": "tu_position", "tir_ostir": "tir_ostir", "sd_dg": "sd_dg",
    "spacer_bp": "spacer_bp", "start_codon": "start_codon", "mrna_dg": "mrna_dg",
    **{c: c for c in CODON}, "length_aa": "length_aa", "mw": "mw_da", "pi": "pi",
    "gravy": "gravy", "tm_segments": "tm_segments", "replicon_type": "replicon_type",
    "gc_local_dev": "gc_local_dev",
}  # fmt: skip
#: identity column -> gene-table column
IDENTITY_SOURCES = {c: c for c in ("locus_tag", "tu_id", "detected", "abundance_log2",
                                   "abundance_imputed")}  # fmt: skip
INT_FEATURES = ("n_promoters_tu_to_gene", "tu_position", "dist_to_tu_start_bp", "spacer_bp",
                "rare_codons_2_10", "length_aa", "tm_segments")  # fmt: skip
FLOAT_FEATURES = ("promotech_max", "promotech_tu", "tx_rate_max_tu", "tx_rate_sum_tu",
                  "tir_ostir", "sd_dg", "mrna_dg", "cai", "tai", "enc", "fop", "gc3", "cpb",
                  "at_2_10", "mw", "pi", "gravy", "replichore_position",
                  "gc_local_dev")  # fmt: skip
STR_FEATURES = ("start_codon", "nterm_residue", "replicon_type")
BOOL_COLUMNS = ("detected", "abundance_imputed", "ribosomal_protein", "in_prophage", "in_bgc")
#: per layer: the gene-table / C6 columns it fills (NA everywhere when the layer is missing)
LAYER_GROUPS: dict[str, tuple[str, ...]] = {
    "promoters": ("promoter_score", "promoter_score_tu", "n_promoter_windows_tu"),
    "tx_rates": ("tx_rate_max_tu", "tx_rate_sum_tu"),
    "ostir": ("tir_ostir", "sd_dg", "mrna_dg", "spacer_bp"),
    "codon": CODON,
}
#: the features those layers feed
LAYER_FEATURES: dict[str, tuple[str, ...]] = {
    "promoters": ("promotech_max", "promotech_tu", "n_promoters_tu_to_gene"),
    "tx_rates": ("tx_rate_max_tu", "tx_rate_sum_tu"),
    "ostir": ("tir_ostir", "sd_dg", "mrna_dg", "spacer_bp"),
    "codon": CODON,
}


def mod() -> ModuleType:
    """``giap.determinants.features`` (imported per test so each test fails on its own)."""
    return importlib.import_module(MODULE)


def is_na(x: Any) -> bool:
    if x is None or x is pd.NA:
        return True
    try:
        return bool(pd.isna(x))
    except (TypeError, ValueError):
        return False


# ---------------------------------------------------------------------------------------------
# Model objects


def replicon(cid: str, length: int, topology: str = "circular", kind: str = "chromosome") -> Any:
    from giap.model.coords import Replicon

    return Replicon(cid, length, topology, kind)  # type: ignore[arg-type]


def gene(
    tag: str,
    contig: str,
    start: int,
    end: int,
    strand: str,
    aa: str = "MKV",
    *,
    start_codon: str | None = "ATG",
    feature_type: str = "cds",
    partial: bool = False,
    name: str | None = None,
    product: str = "",
) -> Any:
    from giap.model.coords import Interval
    from giap.model.genes import Gene

    return Gene(
        tag, feature_type, Interval(contig, start, end, strand), aa,  # type: ignore[arg-type]
        "ATG" * max(1, (end - start) // 3), start_codon, name, product, partial,
    )  # fmt: skip


def start_codon_pos(g: Any) -> int:
    """0-based start-codon coordinate (first base in transcription direction)."""
    return int(g.interval.start) if g.interval.strand == "+" else int(g.interval.end) - 1


def tu(tu_id: str, members: Sequence[Any], length: int | None = None) -> Any:
    """A ``Tu`` over ``members`` (transcription order); interval unwrapped across the origin."""
    from giap.annotate.operons import Tu
    from giap.model.coords import Interval

    first = members[0]
    contig, strand = first.interval.contig, first.interval.strand
    left = members[0] if strand == "+" else members[-1]
    right = members[-1] if strand == "+" else members[0]
    start, end = int(left.interval.start), int(right.interval.end)
    if end <= start:  # wraps: unwrap
        assert length is not None
        end += length
    up = Interval(contig, max(0, start - 100), max(1, start), strand)
    return Tu(tu_id, contig, Interval(contig, start, end, strand), strand,
              tuple(g.locus_tag for g in members), up)  # fmt: skip


def mirror_gene(g: Any, length: int, circular: bool) -> Any:
    """Reflect ``g`` through the genome (x -> L - x) and flip its strand."""
    from giap.model.coords import Interval

    s, e = int(g.interval.start), int(g.interval.end)
    n = e - s
    s2 = (length - e) % length if circular else length - e
    flip = "-" if g.interval.strand == "+" else "+"
    return replace(g, interval=Interval(g.interval.contig, s2, s2 + n, flip))


# ---------------------------------------------------------------------------------------------
# The toy world

AA = {
    "G01": "MAKVLE", "G02": "MKLVAE", "G03": "MPLKE", "G04": "VKLAE", "G05": "M",
    "G06": "MSTVK", "G07": "MAKGL", "G08": "MGKLE", "G09": "", "G10": "MAAAA", "G11": "MDDDK",
    "G12": "MCVKE", "G13": "MTWKE",
}  # fmt: skip
START = {"G04": "GTG", "G05": "TTG", "G06": "CTG"}
#: tag -> (contig, start, end, strand, tu_id, tu_position)
LAYOUT: dict[str, tuple[str, int, int, str, str, int]] = {
    "G01": ("chr1", 1000, 1300, "+", "TU_1", 1), "G02": ("chr1", 1450, 1750, "+", "TU_1", 2),
    "G03": ("chr1", 5000, 5600, "-", "TU_2", 1), "G04": ("chr1", 4300, 4950, "-", "TU_2", 2),
    "G05": ("chr1", 9800, 9950, "+", "TU_3", 1), "G06": ("chr1", 200, 500, "+", "TU_3", 2),
    "G07": ("chr1", 6000, 6300, "+", "TU_4", 1), "G08": ("chr1", 7000, 7300, "+", "TU_5", 1),
    "G09": ("chr1", 7500, 7800, "-", "TU_6", 1), "G10": ("chr1", 8000, 8300, "+", "TU_7", 1),
    "G11": ("chr1", 8500, 8800, "+", "TU_8", 1), "G12": ("p1", 100, 400, "+", "TU_9", 1),
    "G13": ("chr1", 3000, 3300, "+", "TU_10", 1),
}  # fmt: skip
TAGS = tuple(sorted(LAYOUT))
NOT_ELIGIBLE = ("G10", "G11")
#: eligible genes in §6.5 order (C1 contig order, start, end, strand, locus_tag)
ELIGIBLE_ORDER = ("G06", "G01", "G02", "G13", "G04", "G03", "G07", "G08", "G09", "G05", "G12")
EXPECTED_DIST = {"G01": 0, "G02": 450, "G03": 0, "G04": 650, "G05": 0, "G06": 400, "G07": 0,
                 "G08": 0, "G09": 0, "G10": 0, "G11": 0, "G12": 0, "G13": 0}  # fmt: skip
EXPECTED_NTERM: dict[str, str | None] = {
    "G01": "A", "G02": "M", "G03": "P", "G04": "V", "G05": None, "G06": "S", "G07": "A",
    "G08": "G", "G09": None, "G12": "C", "G13": "T",
}  # fmt: skip
TIR = {"G02": 0.0, "G03": -2.0}
REGION_OF = {"G08": ("prophage_1", "prophage"), "G09": ("bgc_1", "bgc")}

_FLOAT_SOURCES = ("promoter_score", "promoter_score_tu", "tx_rate_max_tu", "tx_rate_sum_tu",
                  "terminator_conf_downstream", "tir_ostir", "sd_dg", "mrna_dg", "cai", "tai",
                  "enc", "fop", "gc3", "cpb", "at_2_10", "mw_da", "pi", "gravy",
                  "gc_local_dev", "uniop_p_prev", "uniop_p_next")  # fmt: skip
_INT_SOURCES = ("spacer_bp", "rare_codons_2_10", "tm_segments", "length_aa",
                "n_observable_peptides")  # fmt: skip


def _float_value(j: int, i: int) -> float:
    return float(j + 1) + (i + 1) / 100.0


def _int_value(j: int, i: int) -> int:
    return 10 * (j + 1) + i


@dataclass
class World:
    """The seven ``feature_frame`` arguments."""

    genes_tbl: pd.DataFrame
    genes: list[Any]
    tus: list[Any]
    layers: pd.DataFrame
    regions: pd.DataFrame | None
    gene_context: pd.DataFrame | None
    replicons: dict[str, Any]

    def args(self) -> tuple[Any, ...]:
        return (self.genes_tbl, self.genes, self.tus, self.layers, self.regions,
                self.gene_context, self.replicons)  # fmt: skip

    def features(self) -> pd.DataFrame:
        return mod().feature_frame(*self.args())  # type: ignore[no-any-return]

    def copy(self) -> World:
        return copy.deepcopy(self)

    def blank(self, columns: Sequence[str]) -> World:
        """A copy with ``columns`` NA in the gene table and in C6 (a missing layer)."""
        w = self.copy()
        for frame in (w.genes_tbl, w.layers):
            for c in columns:
                if c not in frame.columns:
                    continue
                if pd.api.types.is_integer_dtype(frame[c]):
                    frame[c] = pd.array([pd.NA] * len(frame), dtype="Int64")
                elif pd.api.types.is_float_dtype(frame[c]):
                    frame[c] = np.nan
                else:
                    frame[c] = pd.Series([None] * len(frame), dtype=object, index=frame.index)
        return w


def _typed(values: list[Any], kind: str) -> pd.Series:
    if kind == "int":
        return pd.Series(pd.array(values, dtype="Int64"))
    if kind == "bool":
        return pd.Series(pd.array(values, dtype="boolean"))
    if kind == "float":
        return pd.Series(values, dtype="float64")
    return pd.Series(values, dtype=object)


def toy() -> World:
    """The 13-gene world of the module docstring (fresh objects on every call)."""
    from giap.genetable.assemble import GENES_COLUMNS, GENES_TYPES

    reps = {"chr1": replicon("chr1", 10_000), "p1": replicon("p1", 3_000, "linear", "plasmid")}
    genes = []
    for t in TAGS:
        contig, s, e, strand, _, _ = LAYOUT[t]
        kw: dict[str, Any] = {"start_codon": START.get(t, "ATG")}
        if t == "G07":
            kw |= {"name": "rplB", "product": "50S ribosomal protein L2"}
        if t == "G09":
            kw["feature_type"] = "pseudogene"
        if t == "G10":
            kw["partial"] = True
        genes.append(gene(t, contig, s, e, strand, AA[t], **kw))
    by_tag = {g.locus_tag: g for g in genes}
    tu_members: dict[str, list[Any]] = {}
    for t in TAGS:
        tu_members.setdefault(LAYOUT[t][4], []).append(t)
    tus = []
    for tid in sorted(tu_members, key=lambda x: int(x.split("_")[1])):
        members = sorted(tu_members[tid], key=lambda t: LAYOUT[t][5])
        tus.append(tu(tid, [by_tag[t] for t in members], 10_000))

    n = len(TAGS)
    tu_size = {tid: len(m) for tid, m in tu_members.items()}
    cols: dict[str, list[Any]] = {c: [None] * n for c in GENES_COLUMNS}
    for i, t in enumerate(TAGS):
        contig, s, e, strand, tid, pos = LAYOUT[t]
        g = by_tag[t]
        r = cols
        r["locus_tag"][i], r["source_id"][i] = t, t
        r["feature_type"][i] = g.feature_type
        r["product"][i] = g.product or None
        r["gene"][i] = g.gene
        r["hypothetical"][i] = False
        r["contig"][i], r["start"][i], r["end"][i], r["strand"][i] = contig, s, e, strand
        r["partial"][i], r["crosses_origin"][i] = g.partial, False
        r["replicon_type"][i] = reps[contig].kind
        r["tu_id"][i], r["tu_position"][i], r["tu_size"][i] = tid, pos, tu_size[tid]
        r["region_ids"][i] = [REGION_OF[t][0]] if t in REGION_OF else []
        r["replichore"][i] = None
        r["replichore_position"][i] = None
        r["start_codon"][i] = g.start_codon
        for j, c in enumerate(_FLOAT_SOURCES):
            r[c][i] = _float_value(j, i)
        for j, c in enumerate(_INT_SOURCES):
            r[c][i] = _int_value(j, i)
        r["tir_ostir"][i] = TIR.get(t, r["tir_ostir"][i])
        det = i % 2 == 0
        r["detected"][i] = det
        r["detection_tier"][i] = "detected" if det else "not_detected"
        r["detection_status"][i] = "detected" if det else "not_detected"
        r["abundance_log2"][i] = 20.0 + i / 4 if det else NAN
        r["abundance_imputed"][i] = t == "G13"
        r["p_detect"][i] = NAN if t == "G11" else 0.3 + i / 50
        for c in ("go_terms", "kegg_ko", "ec"):
            r[c][i] = []
        flags = []
        if t == "G07":
            flags.append("cai_reference")
        if t == "G09":
            flags.append("pseudogene")
        if t == "G10":
            flags.append("partial")
        r["flags"][i] = flags
    tbl = pd.DataFrame({c: _typed(cols[c], GENES_TYPES[c]) for c in GENES_COLUMNS})

    from giap.genetable.contract import CONTRACT_COLUMNS

    lay: dict[str, Any] = {}
    for c in CONTRACT_COLUMNS["layers"]:
        if c in tbl.columns:
            lay[c] = tbl[c].copy()
    lay["n_promoter_windows_tu"] = pd.Series(pd.array([i % 3 for i in range(n)], dtype="Int64"))
    lay["tss_distance_bp"] = pd.Series(pd.array([50 + i for i in range(n)], dtype="Int64"))
    for j, c in enumerate(("spacing_dg", "standby_dg")):
        lay[c] = pd.Series([-1.0 - j - i / 10 for i in range(n)], dtype="float64")
    lay["flags"] = pd.Series([[] for _ in range(n)], dtype=object)
    layers = pd.DataFrame({c: lay[c] for c in CONTRACT_COLUMNS["layers"]})

    regions = pd.DataFrame(
        {
            "region_id": [REGION_OF[t][0] for t in sorted(REGION_OF)],
            "locus_tag": sorted(REGION_OF),
            "overlap_frac": [1.0] * len(REGION_OF),
            "subset": [
                "cargo" if REGION_OF[t][1] == "prophage" else None for t in sorted(REGION_OF)
            ],
            "core_biosynthetic": [False] * len(REGION_OF),
            "region_type": [REGION_OF[t][1] for t in sorted(REGION_OF)],
        }
    )
    chr_tags = [t for t in TAGS if LAYOUT[t][0] == "chr1"]
    ctx = pd.DataFrame(
        {
            "locus_tag": list(TAGS),
            "island_id": [None] * n,
            "replichore": [
                "right"
                if t in chr_tags[: len(chr_tags) // 2]
                else ("left" if t in chr_tags else None)
                for t in TAGS
            ],
            "replichore_position": [
                0.05 * (i + 1) if t in chr_tags else NAN for i, t in enumerate(TAGS)
            ],
        }
    )
    return World(tbl, genes, tus, layers, regions, ctx, reps)


def by_tag(df: pd.DataFrame) -> dict[str, pd.Series]:
    return {str(r["locus_tag"]): r for _, r in df.iterrows()}


def gene_table_row(w: World, tag: str) -> pd.Series:
    sel = w.genes_tbl[w.genes_tbl["locus_tag"] == tag]
    assert len(sel) == 1
    return sel.iloc[0]


def layers_row(w: World, tag: str) -> pd.Series:
    sel = w.layers[w.layers["locus_tag"] == tag]
    assert len(sel) == 1
    return sel.iloc[0]


def same(a: Any, b: Any) -> bool:
    """Cell equality with NA == NA (floats compared exactly: values are copied, not computed)."""
    if is_na(a) or is_na(b):
        return is_na(a) and is_na(b)
    if isinstance(a, bool | np.bool_) or isinstance(b, bool | np.bool_):
        return bool(a) == bool(b)
    if isinstance(a, str) or isinstance(b, str):
        return str(a) == str(b)
    return float(a) == float(b)


def model_columns_of(mapping: Mapping[str, str], block: str) -> set[str]:
    return {c for c, b in mapping.items() if b == block}


#: every encoded model column of the five blocks
EXPECTED_BLOCKS_UNION: frozenset[str] = frozenset().union(*EXPECTED_BLOCKS.values())


# ---------------------------------------------------------------------------------------------
# Conformance: arguments from contract frames (C1, C2, C5, C6, C7/C8) and an assembled gene table


def genes_from_frame(df: pd.DataFrame) -> list[Any]:
    """C2 rows -> ``Gene`` (as ``giap.steps.mapping._gene_model`` builds them)."""
    from giap.model.coords import Interval
    from giap.model.genes import Gene

    return [
        Gene(
            r["locus_tag"], r["feature_type"],
            Interval(r["contig"], int(r["start"]), int(r["end"]), r["strand"]),
            r["aa"] if isinstance(r["aa"], str) else "", r["nt"],
            r["start_codon"] if isinstance(r["start_codon"], str) else None,
            r["gene"] if isinstance(r["gene"], str) else None,
            r["product"] if isinstance(r["product"], str) else "", bool(r["partial"]),
            bool(r["hypothetical"]),
        )
        for r in df.to_dict("records")
    ]  # fmt: skip


def tus_from_frame(df: pd.DataFrame) -> list[Any]:
    """C5 rows -> ``Tu`` (``upstream`` from the frame when present, else a 1-bp stand-in)."""
    from giap.annotate.operons import Tu
    from giap.model.coords import Interval

    out = []
    for r in df.to_dict("records"):
        s, e = int(r["start"]), int(r["end"])
        if "upstream_start" in r and not is_na(r["upstream_start"]):
            up = Interval(
                r["contig"], int(r["upstream_start"]), int(r["upstream_end"]), r["strand"]
            )
        else:
            up = Interval(r["contig"], s, s + 1, r["strand"])
        out.append(Tu(r["tu_id"], r["contig"], Interval(r["contig"], s, e, r["strand"]),
                      r["strand"], tuple(r["genes"]), up))  # fmt: skip
    return out


def replicons_from_frame(df: pd.DataFrame) -> dict[str, Any]:
    return {
        str(r["contig"]): replicon(str(r["contig"]), int(r["length"]), r["topology"], r["kind"])
        for r in df.to_dict("records")
    }


def regions_arg(regions: pd.DataFrame, region_genes: pd.DataFrame) -> pd.DataFrame:
    """C8 with the C7 ``region_type`` joined on (the ``regions`` argument convention)."""
    types = regions[["region_id", "region_type"]].drop_duplicates("region_id")
    return region_genes.merge(types, on="region_id", how="left", validate="many_to_one")


def eligible_tags(genes_tbl: pd.DataFrame) -> list[str]:
    """§7.3 (``regions_table._eligible``): non-partial with a finite ``p_detect``, table order."""
    p = pd.to_numeric(genes_tbl["p_detect"], errors="coerce").to_numpy(float)
    part = genes_tbl["partial"].astype("boolean").fillna(False).to_numpy(bool)
    keep = ~part & np.isfinite(p)
    return genes_tbl.loc[keep, "locus_tag"].astype(str).tolist()


#: (feature, gene-table column) pairs checked on real tables (``locus_tag`` is the join key)
CONFORMANCE_SOURCES: list[tuple[str, str]] = sorted(
    {
        **GENE_TABLE_SOURCES,
        **{k: v for k, v in IDENTITY_SOURCES.items() if k != "locus_tag"},
    }.items()
)
