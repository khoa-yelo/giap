"""Independent checks for giapsim's oriC, ori-ter gradient and null scenarios (issue M5-04).

Written from docs/specs/M5-04-sim-oric-gradient.md (GIAP_SPEC.md §7.3.5(e), §11.3), without reading
the simulator: the replichore truth, the gradient shapes, circular interval overlap, the chr1
intergenic gaps, Moran's I at lag 1 and a hand-built ``Layout`` for the spec's coordinate examples.
Coordinates are 0-based half-open; ``end`` may pass the replicon length (unwrapped, §6.5).
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

ORIC_LEN = 232
#: C15 columns (the ``parse_bakta`` output), in order.
FEATURE_COLUMNS: tuple[str, ...] = (
    "feature_id", "type", "contig", "start", "end", "strand", "gene", "product", "amino_acid",
    "anticodon",
)  # fmt: skip
GENE_COLUMNS: tuple[str, ...] = (
    "locus_tag", "gene", "product", "feature_type", "contig", "start", "end", "strand",
    "partial", "crosses_origin", "hypothetical", "aa", "nt", "start_codon",
)  # fmt: skip


# --------------------------------------------------------------------------- replichores


def replichore(m: float, o: float, length: float) -> tuple[str, float]:
    """(arm, position) of a point ``m`` for an oriC midpoint ``o`` on a circle of ``length``:
    clockwise distance d = (m - o) mod L; ``right`` with d / (L/2) when d <= L/2, else ``left``
    with (L - d) / (L/2)."""
    d = (m - o) % length
    half = length / 2.0
    if d <= half:
        return "right", d / half
    return "left", (length - d) / half


def midpoint(start: float, end: float, length: float) -> float:
    """Midpoint of an (unwrapped) interval, reduced mod ``length``."""
    return ((start + end) / 2.0) % length


def gradient(p: float, oriter_log2: float, shape: str) -> float:
    """TU-level shift at replichore position ``p``: cosine (a/2) cos(pi p), linear a (0.5 - p)."""
    if shape == "cosine":
        return oriter_log2 / 2.0 * math.cos(math.pi * p)
    if shape == "linear":
        return oriter_log2 * (0.5 - p)
    raise ValueError(shape)


# --------------------------------------------------------------------------- intervals


def overlaps_circular(a: tuple[int, int], b: tuple[int, int], length: int) -> bool:
    """Do half-open intervals ``a`` and ``b`` (either may be unwrapped past ``length``) overlap
    on the circle?"""
    for k in (-1, 0, 1):
        s, e = b[0] + k * length, b[1] + k * length
        if a[0] < e and s < a[1]:
            return True
    return False


def chr1_gaps(genes: pd.DataFrame, length: int) -> list[tuple[int, int]]:
    """Intergenic gaps of chr1 as (left end, right start) in circular genome order, unwrapped so
    that right > left; the wrap gap (last gene to first gene) included. Overlapping neighbours
    give right <= left (zero or negative width)."""
    g = genes[genes["contig"] == "chr1"]
    iv = sorted((int(s), int(e)) for s, e in zip(g["start"], g["end"], strict=True))
    out = []
    for (_, e0), (s1, _) in itertools.pairwise(iv):
        out.append((e0, s1))
    # wrap gap: the furthest end (mod L, after the last start) to the first start (+ L)
    last_end = max(e for _, e in iv)
    first_start = iv[0][0]
    if last_end > length:  # a gene crosses the origin
        out.append((last_end - length, first_start))
    else:
        out.append((last_end, first_start + length))
    return out


def circ_dist(a: float, b: float, length: float) -> float:
    d = abs(a - b) % length
    return min(d, length - d)


# --------------------------------------------------------------------------- statistics


def morans_i_lag1(x: np.ndarray) -> float:
    """Moran's I with lag-1 neighbours along a chain (w_{i,i+1} = w_{i+1,i} = 1)."""
    z = np.asarray(x, dtype=float) - float(np.mean(x))
    n = len(z)
    return float(n * np.sum(z[1:] * z[:-1]) / ((n - 1) * np.sum(z * z)))


def pearson(x: Sequence[float], y: Sequence[float]) -> float:
    return float(np.corrcoef(np.asarray(x, float), np.asarray(y, float))[0, 1])


# --------------------------------------------------------------------------- hand-built layouts


def hand_layout(
    length: int,
    genes: Sequence[tuple[int, int]],
    oric: tuple[int, int] | None,
    *,
    tus: Sequence[Sequence[int]] | None = None,
    plasmid: tuple[int, Sequence[tuple[int, int]]] | None = None,
) -> Any:
    """A ``giapsim.layout.Layout`` by hand: chr1 of ``length`` with ``genes`` (start, end) on ``+``
    (``end > length`` marks an origin-crossing gene), TUs ``tus`` (gene index lists, default
    singletons), an optional plasmid ``p1`` (length, genes) and an optional oriC (start, end)."""
    from giapsim.layout import Layout

    rows: list[tuple[Any, ...]] = []
    contigs = [("chr1", length, genes)]
    if plasmid is not None:
        contigs.append(("p1", plasmid[0], plasmid[1]))
    k = 0
    for cid, _, gs in contigs:
        for s, e in gs:
            k += 1
            aa = "M" + "AKLEGRSTVD" * 6
            rows.append((f"GSIM_{5 * k:05d}", None, f"protein {k}", "cds", cid, int(s), int(e),
                         "+", False, e > length if cid == "chr1" else False, False, aa,
                         "ATG" + "GCT" * 60 + "TAA", "ATG"))  # fmt: skip
    gdf = pd.DataFrame(rows, columns=list(GENE_COLUMNS)).astype({"start": "int64", "end": "int64"})
    gdf["amr"] = False
    tags = gdf["locus_tag"].tolist()
    n_chr = len(genes)
    groups = [list(t) for t in tus] if tus is not None else [[i] for i in range(n_chr)]
    groups += [[i] for i in range(n_chr, len(tags))]
    tu_rows = []
    for j, grp in enumerate(groups):
        sub = gdf.iloc[grp]
        tu_rows.append((f"TU_{j + 1}", str(sub["contig"].iloc[0]), int(sub["start"].min()),
                        int(sub["end"].max()), "+", [tags[i] for i in grp]))  # fmt: skip
    tus_df = pd.DataFrame(tu_rows, columns=["tu_id", "contig", "start", "end", "strand", "genes"])
    reps = pd.DataFrame(
        [(cid, n, "circular", "chromosome" if cid == "chr1" else "plasmid")
         for cid, n, _ in contigs],
        columns=["contig", "length", "topology", "kind"],
    ).astype({"length": "int64"})  # fmt: skip
    regions = pd.DataFrame(
        {c: pd.Series([], dtype=object) for c in
         ("region_id", "region_type", "contig", "start", "end", "source_tool", "score")}
    )  # fmt: skip
    region_genes = pd.DataFrame(
        {c: pd.Series([], dtype=object) for c in
         ("region_id", "locus_tag", "overlap_frac", "subset", "core_biosynthetic")}
    )  # fmt: skip
    pairs = pd.DataFrame({"gene_a": [], "gene_b": [], "uniop_p": []})
    tu_of = {tags[i]: f"TU_{j + 1}" for j, grp in enumerate(groups) for i in grp}
    truth = pd.DataFrame({"locus_tag": tags, "tu_id": [tu_of[t] for t in tags],
                          "region_ids": [[] for _ in tags],
                          "region_effect": np.zeros(len(tags))})  # fmt: skip
    xrefs = pd.DataFrame({"locus_tag": [], "db": [], "id": []})
    feats = pd.DataFrame({c: pd.Series([], dtype=object) for c in FEATURE_COLUMNS})
    if oric is not None:
        feats = pd.DataFrame(
            [("chr1_oric1", "oric", "chr1", int(oric[0]), int(oric[1]), ".", None, None, None,
              None)],
            columns=list(FEATURE_COLUMNS),
        ).astype({"start": "int64", "end": "int64"})  # fmt: skip
    return Layout(
        replicons=reps, genes=gdf, tus=tus_df, operon_pairs=pairs, regions=regions,
        region_genes=region_genes, truth=truth, xrefs=xrefs, features=feats,
    )  # fmt: skip


def oric_of(layout: Any) -> tuple[int, int, float] | None:
    """(start, end, midpoint mod L) of the layout's chr1 oriC row, None without one."""
    f = layout.features
    f = f[(f["type"] == "oric") & (f["contig"] == "chr1")]
    if len(f) == 0:
        return None
    s, e = int(f["start"].iloc[0]), int(f["end"].iloc[0])
    length = int(layout.replicons.set_index("contig").loc["chr1", "length"])
    return s, e, midpoint(s, e, length)
