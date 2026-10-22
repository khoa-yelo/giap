"""Inputs for the M8-06a ``function`` step tests (docs/specs/M8-06-function-step-and-run.md).

* :func:`toy4` is example 4's 4-gene gene table (g1-g4 in TU_1-TU_4).
* :func:`world_pre` is a ~600-gene ``genes_pre`` frame (the M4-13 ``detectability_data.synthetic``
  layout, 300 TUs, with covariates, so ``detectability`` fits a non-constant model) carrying
  function columns from a synthetic C3 ``xrefs`` frame, detection and abundance. Its function DB
  comes from :mod:`helpers.function_db`.
* Writers for the step's files: the gene table (``write_genes_parquet``), C1 replicons and
  ``tu_folds.parquet`` (``tu_id`` string, ``fold`` int64).

Nothing here imports ``giap.steps.function`` or ``giap.io.function_tables``.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]

from helpers.detectability_data import synthetic

SEED = 7
CONTIG = "chr1"
GENOME_LENGTH = 1_000_000


def natural_key(s: str) -> tuple[tuple[int, int | str], ...]:
    return tuple((0, int(p)) if p.isdigit() else (1, p) for p in re.split(r"(\d+)", s) if p)


def full_frame(cols: Mapping[str, Sequence[Any]]) -> pd.DataFrame:
    """A §9.3 gene table: every ``GENES_COLUMNS`` column, defaults where ``cols`` gives none."""
    from giap.genetable.assemble import GENES_COLUMNS, GENES_TYPES

    n = len(next(iter(cols.values())))
    default = {"list": lambda: [], "bool": lambda: False, "int": lambda: None,
               "float": lambda: math.nan, "str": lambda: None}  # fmt: skip
    data: dict[str, Any] = {}
    for c in GENES_COLUMNS:
        if c in cols:
            data[c] = list(cols[c])
        else:
            data[c] = [default[GENES_TYPES[c]]() for _ in range(n)]
    df = pd.DataFrame({c: pd.Series(v, dtype=object) for c, v in data.items()})
    return df


def write_gene_table(path: Path, df: pd.DataFrame) -> Path:
    from giap.io.genes_tsv import write_genes_parquet

    path.parent.mkdir(parents=True, exist_ok=True)
    write_genes_parquet(path, df)
    return path


def write_replicons(path: Path, topology: str = "circular", length: int = GENOME_LENGTH) -> Path:
    from giap.io.parquet import write_parquet_deterministic

    path.parent.mkdir(parents=True, exist_ok=True)
    cols = [("contig", pa.string()), ("length", pa.int64()), ("topology", pa.string()),
            ("kind", pa.string())]  # fmt: skip
    write_parquet_deterministic(path, cols, [(CONTIG, length, topology, "chromosome")])
    return path


def write_tu_folds(path: Path, folds: Mapping[str, int]) -> Path:
    """``genetable/tu_folds.parquet`` as item 1 describes it (header-only for ``{}``)."""
    from giap.io.parquet import write_parquet_deterministic

    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [(t, int(folds[t])) for t in sorted(folds, key=natural_key)]
    write_parquet_deterministic(path, [("tu_id", pa.string()), ("fold", pa.int64())], rows)
    return path


def read_tu_folds(path: Path) -> dict[str, int]:
    df = pd.read_parquet(path)
    return {str(t): int(f) for t, f in zip(df["tu_id"], df["fold"], strict=True)}


# ----------------------------------------------------------------------------- example 4

TOY_DETECTED = (True, True, False, False)
TOY_P = (0.9, 0.5, 0.4, 0.3)
TOY_COG = ("J", "J", "J", None)


def toy4(**over: Sequence[Any]) -> pd.DataFrame:
    """Example 4: g1-g4 in TU_1-TU_4, non-partial, COG J, J, J, NA, detected T, T, F, F,
    ``p_detect`` 0.9, 0.5, 0.4, 0.3, no abundance. ``over`` replaces whole columns."""
    n = 4
    cols: dict[str, Sequence[Any]] = {
        "locus_tag": [f"g{i + 1}" for i in range(n)],
        "feature_type": ["cds"] * n,
        "contig": [CONTIG] * n,
        "start": [1000 * i for i in range(n)],
        "end": [1000 * i + 900 for i in range(n)],
        "strand": ["+"] * n,
        "partial": [False] * n,
        "crosses_origin": [False] * n,
        "hypothetical": [False] * n,
        "tu_id": [f"TU_{i + 1}" for i in range(n)],
        "tu_position": [1] * n,
        "tu_size": [1] * n,
        "detected": list(TOY_DETECTED),
        "detection_status": ["detected" if d else "not_detected" for d in TOY_DETECTED],
        "p_detect": list(TOY_P),
        "abundance_imputed": [False] * n,
        "cog_category": list(TOY_COG),
        "cog_id": ["COG0001", "COG0002", "COG0003", None],
        "annotation_source": ["bakta_full"] * n,
    }
    cols.update(over)
    return full_frame(cols)


# ----------------------------------------------------------------------------- the 600-gene world

N_TUS = 300
COG_POOL = tuple(f"COG{i:04d}" for i in range(1, 46))
GO_POOL = tuple(f"GO:{n:07d}" for n in range(1201, 100_000, 1667))  # 60 sim-like leaf IDs
KO_POOL = tuple(f"K{n:05d}" for n in range(101, 30_000, 997))  # 31 sim-like KOs


def world_xrefs(tags: Sequence[str], seed: int = SEED) -> pd.DataFrame:
    """A C3 ``xrefs`` frame in the giapsim shape: ~80 % of genes get a COG, ~60 % of those GO
    terms and a KO."""
    r = np.random.default_rng(seed + 101)
    rows: list[tuple[str, str, str]] = []
    for t in tags:
        if r.random() >= 0.8:
            continue
        rows.append((t, "COG", str(r.choice(COG_POOL))))
        if r.random() < 0.6:
            for g in sorted(set(r.choice(GO_POOL, size=int(r.integers(1, 3))).tolist())):
                rows.append((t, "GO", str(g)))
        if r.random() < 0.6:
            rows.append((t, "KEGG", str(r.choice(KO_POOL))))
    return pd.DataFrame(rows, columns=["locus_tag", "db", "id"])


def world_pre(xrefs_db: Path | None = None, seed: int = SEED) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(``genes_pre`` frame, its C3 ``xrefs``). ``xrefs_db`` (a function DB root with
    ``cog/cog-20.def.tab``) fills ``cog_category``; without it the category is derived from the
    same rule (:func:`helpers.function_db.cog_letters`)."""
    from giap.genetable.assemble import function_columns
    from giap.io.genes_tsv import read_cog_def
    from helpers.function_db import cog_letters

    gs = synthetic(N_TUS, seed)
    tags = [str(t) for t in gs.detected.index]
    n = len(tags)
    tu = gs.tu_of_gene.reindex(tags).astype(str).tolist()
    pos: list[int] = []
    size: dict[str, int] = {}
    for t in tu:
        size[t] = size.get(t, 0) + 1
        pos.append(size[t])
    det = gs.detected.reindex(tags).astype(bool).to_numpy()
    r = np.random.default_rng(seed + 202)
    abundance = np.where(det, 20.0 + r.normal(0.0, 2.0, n), np.nan)
    xrefs = world_xrefs(tags, seed)
    if xrefs_db is not None:
        cog_def = read_cog_def(xrefs_db / "cog" / "cog-20.def.tab")
    else:
        ids = sorted(set(xrefs.loc[xrefs["db"] == "COG", "id"]))
        cog_def = pd.DataFrame({"cog_id": ids, "cog_category": [cog_letters(c) for c in ids]})
    fn = function_columns(xrefs, cog_def).set_index("locus_tag").reindex(tags)
    cov = gs.cov.reindex(tags)
    cols: dict[str, Sequence[Any]] = {
        "locus_tag": tags,
        "feature_type": ["cds"] * n,
        "contig": [CONTIG] * n,
        "start": [1000 * i for i in range(n)],
        "end": [1000 * i + 900 for i in range(n)],
        "strand": ["+"] * n,
        "length_aa": [int(v) for v in cov["length_aa"]],
        "partial": [False] * n,
        "crosses_origin": [False] * n,
        "hypothetical": [False] * n,
        "tu_id": tu,
        "tu_position": pos,
        "tu_size": [size[t] for t in tu],
        "mw_da": [float(v) for v in cov["mw_da"]],
        "pi": [float(v) for v in cov["pi"]],
        "gravy": [float(v) for v in cov["gravy"]],
        "tm_segments": [int(v) for v in cov["tm_segments"]],
        "n_observable_peptides": [int(v) for v in cov["n_observable_peptides"]],
        "detected": det.tolist(),
        "detection_status": ["detected" if d else "not_detected" for d in det],
        "abundance_log2": abundance.tolist(),
        "abundance_imputed": [False] * n,
        "abundance_type": ["ibaq" if d else None for d in det],
        "abundance_scale": ["absolute"] * n,
        "cog_id": [None if pd.isna(v) else v for v in fn["cog_id"]],
        "cog_category": [None if pd.isna(v) else v for v in fn["cog_category"]],
        "go_terms": [list(v) if isinstance(v, list) else [] for v in fn["go_terms"]],
        "kegg_ko": [list(v) if isinstance(v, list) else [] for v in fn["kegg_ko"]],
        "ec": [[] for _ in range(n)],
        "annotation_source": ["bakta_full"] * n,
    }
    return full_frame(cols), xrefs


def oracle_tables(
    gene_table: Path,
    replicons: Path,
    folds: Mapping[str, int] | None,
    resources: Mapping[str, Path],
    cfg: Any,
    out: Path,
) -> tuple[dict[str, Path], list[str]]:
    """Spec item 2's chain on the step's own files, written with the M8-06a writers: the three
    table paths and the run-flag union (first-seen order). ``folds`` empty means None
    (lead decision 3)."""
    from giap.function.collections import build_collections
    from giap.function.gsea import gsea_table
    from giap.function.ora import detection_enrichment
    from giap.function.summary import category_summary
    from giap.genetable.eligibility import eligible
    from giap.io.function_resources import read_resources
    from giap.io.function_tables import (
        write_category_summary_tsv,
        write_enrichment_detection_tsv,
        write_enrichment_gsea_tsv,
    )
    from giap.io.genes_tsv import read_genes_parquet
    from giap.rng import for_step

    genes = read_genes_parquet(gene_table)
    reps = pd.read_parquet(replicons)
    res = read_resources(dict(resources))
    el = eligible(genes)
    colls, flags = build_collections(genes, el, res)
    alpha = cfg.stats.alpha.enrichment
    summary = category_summary(genes, colls, cfg.stats.bootstrap,
                               for_step("function.bootstrap", cfg.seed))  # fmt: skip
    det = detection_enrichment(genes, colls, tuple(cfg.function.ora_size), alpha,
                               folds=dict(folds) if folds else None)  # fmt: skip
    gsea = gsea_table(genes, colls, tuple(cfg.function.gsea_size), cfg.function.gsea, alpha,
                      for_step("function.gsea", cfg.seed), reps)  # fmt: skip
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        "category_summary": out / "category_summary.tsv",
        "enrichment_detection": out / "enrichment_detection.tsv",
        "enrichment_gsea": out / "enrichment_gsea.tsv",
    }
    write_category_summary_tsv(paths["category_summary"], summary)
    write_enrichment_detection_tsv(paths["enrichment_detection"], det)
    write_enrichment_gsea_tsv(paths["enrichment_gsea"], gsea)
    union = list(dict.fromkeys([*flags, *det.attrs["run_flags"], *gsea.attrs["run_flags"]]))
    return paths, union
