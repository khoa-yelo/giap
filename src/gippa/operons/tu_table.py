"""``operons.tsv``: one row per C5 transcription unit (GIPPA_SPEC.md §7.4.7, §9.3;
docs/specs/M6-07-flags-and-operons-table.md step 5).

Coordinates stay 0-based half-open and unwrapped across the origin (§6.5); only the M6-08
writer converts them. Pure and deterministic: no I/O, no RNG.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pandas as pd

__all__ = ["OPERON_COLUMNS", "operons_table"]

#: §9.3 ``operons.tsv`` columns, in table order
OPERON_COLUMNS: tuple[str, ...] = (
    "tu_id", "contig", "start", "end", "strand", "genes", "n_genes", "detected_frac",
    "median_abundance_log2", "max_abs_delta_log2", "n_flags",
)  # fmt: skip
_INT = ("start", "end", "n_genes", "n_flags")


def _gene_list(v: object) -> list[str]:
    return [str(g) for g in v] if isinstance(v, list | tuple | np.ndarray) else []


def operons_table(
    tus: pd.DataFrame, pairs: pd.DataFrame, genes: pd.DataFrame, flagged: pd.DataFrame
) -> pd.DataFrame:
    """Per TU: genes (transcription order), size, detected fraction of the non-partial genes,
    median log2 abundance of the detected genes, max |Delta| over internal adjacent pairs with
    both genes detected, and ``n_flags`` = flagged OP pairs with both genes in the TU."""
    g = genes.drop_duplicates("locus_tag").set_index("locus_tag")
    partial = g["partial"].astype("boolean").fillna(False).astype(bool).to_dict()
    detected = g["detected"].astype("boolean").fillna(False).astype(bool).to_dict()
    abund = pd.to_numeric(g["abundance_log2"], errors="coerce").to_dict()
    delta = {
        frozenset((str(a), str(b))): abs(float(d))
        for a, b, d in zip(pairs["gene_a"], pairs["gene_b"], pairs["delta_log2"], strict=True)
        if pd.notna(d)
    }
    op = flagged.loc[flagged["uniop_class"] == "OP"]
    flag_pairs = list(zip(op["gene_a"].astype(str), op["gene_b"].astype(str), strict=True))
    recs = []
    for t in tus.itertuples(index=False):
        tags = _gene_list(t.genes)
        eligible = [x for x in tags if not partial.get(x, False)]
        det = [x for x in tags if detected.get(x, False)]
        vals = [abund[x] for x in det if x in abund and not math.isnan(abund[x])]
        ds = [
            delta[frozenset((a, b))]
            for a, b in itertools.pairwise(tags)
            if detected.get(a, False) and detected.get(b, False) and frozenset((a, b)) in delta
        ]
        members = set(tags)
        recs.append(
            {
                "tu_id": t.tu_id,
                "contig": t.contig,
                "start": t.start,
                "end": t.end,
                "strand": t.strand,
                "genes": tags,
                "n_genes": len(tags),
                "detected_frac": (
                    sum(x in det for x in eligible) / len(eligible) if eligible else math.nan
                ),
                "median_abundance_log2": float(np.median(vals)) if vals else math.nan,
                "max_abs_delta_log2": max(ds) if ds else math.nan,
                "n_flags": sum(a in members and b in members for a, b in flag_pairs),
            }
        )
    if not recs:
        return pd.DataFrame({c: pd.Series(dtype=object) for c in OPERON_COLUMNS})
    out = pd.DataFrame(recs, columns=list(OPERON_COLUMNS))
    for c in _INT:
        out[c] = out[c].astype("int64")
    rank = {c: k for k, c in enumerate(dict.fromkeys(tus["contig"]))}
    out = out.assign(_rank=out["contig"].map(rank))
    out = out.sort_values(["_rank", "start"], kind="stable").drop(columns="_rank")
    return out.reset_index(drop=True)
