"""Writer of the M4-06 contract tables C1-C8 and C15 for a gippasim layout (issues M4-07, M5-04).

The column names and logical types below are a literal copy of the M2-side tables of
``gippa.genetable.contract.TABLES`` (a test compares them, catching drift both ways); gippasim imports
nothing from ``gippa`` so the truth stays independent of the code under test. Parquet files hold
0-based half-open coordinates and Arrow ``list<string>`` list columns; every file is a pure function
of the layout (byte-identical for the same layout).
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

if TYPE_CHECKING:
    from gippasim.layout import Layout

__all__ = ["ANNOTATION_SOURCE", "CONTRACT_COLUMNS", "CONTRACT_TYPES", "OPERON_TRUTH_TYPES",
           "REQUIRED_COLUMNS", "layers_frame", "read_operon_truth", "write_contract_tables",
           "write_operon_truth"]  # fmt: skip

ANNOTATION_SOURCE = "bakta_full"  # §9.3 enum (M4-15 lead decision 7)

_LAYER_NA: tuple[tuple[str, str], ...] = (
    ("promoter_score", "float"), ("promoter_score_tu", "float"), ("n_promoter_windows_tu", "int"),
    ("tx_rate_max_tu", "float"), ("tx_rate_sum_tu", "float"), ("tss_distance_bp", "int"),
    ("terminator_conf_downstream", "float"), ("tir_ostir", "float"), ("sd_dg", "float"),
    ("mrna_dg", "float"), ("spacing_dg", "float"), ("standby_dg", "float"), ("spacer_bp", "int"),
)  # fmt: skip
_CODON_NA: tuple[tuple[str, str], ...] = (
    ("cai", "float"), ("tai", "float"), ("enc", "float"), ("fop", "float"), ("gc3", "float"),
    ("cpb", "float"), ("at_2_10", "float"), ("rare_codons_2_10", "int"),
)  # fmt: skip

#: (column, logical dtype) per written table, in contract order; ``?`` marks an optional column
#: (written only when the layout frame has it).
CONTRACT_TYPES: Mapping[str, tuple[tuple[str, str], ...]] = {
    "replicons": (("contig", "str"), ("length", "int"), ("topology", "str"), ("kind", "str")),
    "genes": (
        ("locus_tag", "str"), ("gene", "str"), ("product", "str"), ("feature_type", "str"),
        ("contig", "str"), ("start", "int"), ("end", "int"), ("strand", "str"),
        ("partial", "bool"), ("crosses_origin", "bool"), ("hypothetical", "bool"), ("aa", "str"),
        ("nt", "str"), ("start_codon", "str"), ("ribosomal_protein", "bool?"), ("amr", "bool?"),
    ),
    "xrefs": (("locus_tag", "str"), ("db", "str"), ("id", "str")),
    "tus": (
        ("tu_id", "str"), ("contig", "str"), ("start", "int"), ("end", "int"), ("strand", "str"),
        ("genes", "list"),
    ),
    "layers": (
        ("locus_tag", "str"), ("tu_id", "str"), ("tu_position", "int"), ("tu_size", "int"),
        ("uniop_p_prev", "float"), ("uniop_p_next", "float"), ("region_ids", "list"),
        ("gc_local_dev", "float?"), *_LAYER_NA, ("start_codon", "str"), *_CODON_NA,
        ("flags", "list"),
    ),
    "regions": (
        ("region_id", "str"), ("region_type", "str"), ("contig", "str"), ("start", "int"),
        ("end", "int"), ("source_tool", "str"), ("score", "float"), ("genes", "list?"),
        ("notes", "str?"),
    ),
    "region_genes": (
        ("region_id", "str"), ("locus_tag", "str"), ("overlap_frac", "float"), ("subset", "str"),
        ("core_biosynthetic", "bool"),
    ),
    # C15 (M5-04): parse_bakta's columns; gene, product, amino_acid, anticodon nullable
    "features": (
        ("feature_id", "str"), ("type", "str"), ("contig", "str"), ("start", "int"),
        ("end", "int"), ("strand", "str"), ("gene", "str"), ("product", "str"),
        ("amino_acid", "str"), ("anticodon", "str"),
    ),
}  # fmt: skip

#: Column names per written table (equal to ``gippa.genetable.contract.TABLES``; drift test).
CONTRACT_COLUMNS: Mapping[str, tuple[str, ...]] = {
    k: tuple(name for name, _ in cols) for k, cols in CONTRACT_TYPES.items()
}

#: Required (always written) columns per table.
REQUIRED_COLUMNS: Mapping[str, tuple[str, ...]] = {
    k: tuple(name for name, dtype in cols if not dtype.endswith("?"))
    for k, cols in CONTRACT_TYPES.items()
}

_ARROW: dict[str, pa.DataType] = {
    "str": pa.string(),
    "int": pa.int64(),
    "float": pa.float64(),
    "bool": pa.bool_(),
    "list": pa.list_(pa.string()),
}
_PANDAS_NA: dict[str, tuple[str, object]] = {"int": ("Int64", pd.NA), "float": ("float64", np.nan)}
_FILES: dict[str, str] = {k: f"{k}.parquet" for k in CONTRACT_TYPES}


def _uniop_lookup(pairs: pd.DataFrame) -> dict[tuple[str, str], float]:
    out: dict[tuple[str, str], float] = {}
    for a, b, p in zip(pairs["gene_a"], pairs["gene_b"], pairs["uniop_p"], strict=True):
        out[(a, b)] = out[(b, a)] = float(p)
    return out


def layers_frame(layout: Layout) -> pd.DataFrame:
    """C6 rows (one per gene, gene order): TU columns, UniOP p in transcription order, regions,
    start codon and flags; every other M2-20 layer column all-NA."""
    genes = layout.genes
    tags: list[str] = genes["locus_tag"].tolist()
    tu_id: dict[str, str] = {}
    tu_pos: dict[str, int] = {}
    tu_size: dict[str, int] = {}
    for tid, members in zip(layout.tus["tu_id"], layout.tus["genes"], strict=True):
        for k, g in enumerate(members):
            tu_id[g], tu_pos[g], tu_size[g] = str(tid), k + 1, len(members)
    regions: dict[str, list[str]] = {}
    for rid, g in zip(layout.region_genes["region_id"], layout.region_genes["locus_tag"],
                      strict=True):  # fmt: skip
        regions.setdefault(g, []).append(str(rid))
    pm = _uniop_lookup(layout.operon_pairs)
    contig, strand = genes["contig"].tolist(), genes["strand"].tolist()
    circular = set(layout.replicons.loc[layout.replicons["topology"] == "circular", "contig"])
    n = len(tags)
    # genome-order neighbours on the same contig; on a circular contig of > 2 genes the first and
    # last genes are neighbours too (the wrap pair, M2-12)
    up_of: list[int | None] = [i - 1 if i > 0 and contig[i - 1] == contig[i] else None
                               for i in range(n)]  # fmt: skip
    down_of: list[int | None] = [i + 1 if i + 1 < n and contig[i + 1] == contig[i] else None
                                 for i in range(n)]  # fmt: skip
    lo = 0
    while lo < n:
        hi = lo
        while hi + 1 < n and contig[hi + 1] == contig[lo]:
            hi += 1
        if contig[lo] in circular and hi - lo + 1 > 2:
            up_of[lo], down_of[hi] = hi, lo
        lo = hi + 1
    prev = np.full(n, np.nan)
    nxt = np.full(n, np.nan)
    for i in range(n):
        u, d = up_of[i], down_of[i]
        up = None if u is None else pm.get((tags[u], tags[i]))
        down = None if d is None else pm.get((tags[i], tags[d]))
        first, second = (up, down) if strand[i] == "+" else (down, up)
        prev[i] = np.nan if first is None else first
        nxt[i] = np.nan if second is None else second
    flags = []
    for part, cross, ft in zip(genes["partial"], genes["crosses_origin"], genes["feature_type"],
                               strict=True):  # fmt: skip
        f = ["partial"] if part else []
        f += ["crosses_origin"] if cross else []
        f += ["pseudogene"] if ft == "pseudogene" else []
        flags.append(f)
    data: dict[str, Any] = {
        "locus_tag": tags,
        "tu_id": [tu_id[g] for g in tags],
        "tu_position": np.array([tu_pos[g] for g in tags], dtype=np.int64),
        "tu_size": np.array([tu_size[g] for g in tags], dtype=np.int64),
        "uniop_p_prev": prev,
        "uniop_p_next": nxt,
        "region_ids": [sorted(regions.get(g, [])) for g in tags],
    }
    for name, dtype in (*_LAYER_NA, *_CODON_NA):
        pd_dtype, na = _PANDAS_NA[dtype]
        data[name] = pd.Series([na] * n, dtype=pd_dtype)
    data["start_codon"] = genes["start_codon"].tolist()
    data["flags"] = flags
    cols = REQUIRED_COLUMNS["layers"]
    return pd.DataFrame({c: data[c] for c in cols}, columns=list(cols))


def _arrow_table(key: str, df: pd.DataFrame) -> pa.Table:
    cols = tuple(
        (name, dtype.rstrip("?"))
        for name, dtype in CONTRACT_TYPES[key]
        if not dtype.endswith("?") or name in df.columns
    )
    schema = pa.schema([pa.field(name, _ARROW[dtype]) for name, dtype in cols])
    arrays = []
    for name, dtype in cols:
        values = df[name].tolist()
        if dtype == "list":
            values = [[str(x) for x in v] for v in values]
        else:
            values = [None if not isinstance(v, list) and pd.isna(v) else v for v in values]
            if dtype == "int":
                values = [None if v is None else int(v) for v in values]
            elif dtype == "float":
                values = [None if v is None else float(v) for v in values]
            elif dtype == "bool":
                values = [None if v is None else bool(v) for v in values]
        arrays.append(pa.array(values, type=_ARROW[dtype]))
    return pa.Table.from_arrays(arrays, schema=schema)


def write_contract_tables(layout: Layout, work_dir: Path) -> None:
    """Write C1-C8 (``replicons``, ``genes``, ``xrefs``, ``annotation_meta.json``, ``tus``,
    ``layers``, ``regions``, ``region_genes``) and C15 ``features`` (header-only without an oriC)
    directly under ``work_dir`` (created if needed)."""
    work_dir.mkdir(parents=True, exist_ok=True)
    frames: dict[str, pd.DataFrame] = {
        "replicons": layout.replicons,
        "genes": layout.genes,
        "xrefs": layout.xrefs,
        "tus": layout.tus,
        "layers": layers_frame(layout),
        "regions": layout.regions,
        "region_genes": layout.region_genes,
        "features": layout.features,
    }
    for key, df in frames.items():
        pq.write_table(_arrow_table(key, df), work_dir / _FILES[key], compression="zstd")
    meta = {"annotation_source": ANNOTATION_SOURCE, "flags": []}
    (work_dir / "annotation_meta.json").write_text(
        json.dumps(meta, sort_keys=True, indent=1) + "\n", encoding="utf-8"
    )


#: Operon-scenario truth files ``truth/<key>.tsv`` (M6-03; gippasim truth, not gippa contract
#: tables): (column, logical dtype) in file order.
OPERON_TRUTH_TYPES: Mapping[str, tuple[tuple[str, str], ...]] = {
    "pairs": (
        ("pair_id", "str"), ("class", "str"), ("gene_a", "str"), ("gene_b", "str"),
        ("true_operon", "bool"), ("promoter_between_true", "bool"),
        ("terminator_between_true", "bool"),
    ),
    "missed_operons": (("pair_id", "str"), ("gene_a", "str"), ("gene_b", "str")),
    "discordant_pairs": (
        ("pair_id", "str"), ("gene_a", "str"), ("gene_b", "str"), ("direction", "str"),
        ("element", "str"), ("effect_log2", "float"),
    ),
    "operon_t": (("loc", "float"), ("scale", "float"), ("df", "float")),
}  # fmt: skip


def _tsv_cell(v: Any, dtype: str) -> str:
    if v is None or (not isinstance(v, str) and pd.isna(v)):
        return "NA"
    if dtype == "bool":
        return "true" if bool(v) else "false"
    if dtype == "float":
        return repr(float(v))
    if dtype == "int":
        return str(int(v))
    return str(v)


def write_operon_truth(sim: Any, out_dir: Path) -> dict[str, Path]:
    """Write ``sim.truth`` as ``out_dir/truth/<key>.tsv`` (gippa TSV: ``true``/``false``, ``NA``);
    byte-identical for equal frames. Returns key -> path."""
    (out_dir / "truth").mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for key, cols in OPERON_TRUTH_TYPES.items():
        df = sim.truth[key]
        lines = ["\t".join(c for c, _ in cols)]
        data = [df[c].tolist() for c, _ in cols]
        for row in zip(*data, strict=True):
            lines.append("\t".join(_tsv_cell(v, d) for v, (_, d) in zip(row, cols, strict=True)))
        paths[key] = out_dir / "truth" / f"{key}.tsv"
        paths[key].write_text("\n".join(lines) + "\n", encoding="utf-8")
    return paths


def read_operon_truth(out_dir: Path) -> dict[str, pd.DataFrame]:
    """Read the ``write_operon_truth`` files back with their declared dtypes."""
    out: dict[str, pd.DataFrame] = {}
    for key, cols in OPERON_TRUTH_TYPES.items():
        df = pd.read_csv(out_dir / "truth" / f"{key}.tsv", sep="\t", dtype=str,
                         keep_default_na=False)  # fmt: skip
        for c, d in cols:
            v = df[c]
            if d == "bool":
                df[c] = v.map({"true": True, "false": False}).astype(bool)
            elif d == "float":
                df[c] = pd.to_numeric(v.replace("NA", "nan"), errors="raise").astype(float)
            elif d == "int":
                df[c] = v.map(lambda x: pd.NA if x == "NA" else int(x)).astype("Int64")
            else:
                df[c] = v.astype(object)
        out[key] = df[[c for c, _ in cols]]
    return out
