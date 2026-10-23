"""The M9-26 20-Mb stress fixture: a synthetic ``OUT/`` (``RunView`` frames only, not a pipeline
run) with 25,000 genes on one circular contig (docs/specs/M9-26-scale-and-determinism.md work
item 2, example 1; GIAP_SPEC.md §10.1 "Plots aggregate into bins when a track would draw more
than 20,000 elements").

Gene coordinates come from the giapsim sequence-free layout (``giapsim.layout.simulate_layout``,
no regions, 25,000 chromosome genes at a 220-aa median: one ≈ 20.5-Mb ``chr1``). Detection and
abundance are drawn here from ``numpy.random.default_rng(seed)``: ``DETECTED_FRAC`` of the genes
are ``detected`` with a measured (not imputed) ``abundance_log2 ~ N(20, 2)``, so V2 has more than
20,000 points on its one panel; the rest are ``not_detected`` (no abundance). ``genes.tsv`` is
written 1-based closed under the ``genes`` schema (``helpers.viz_out_tree``), ``preflight.json``
names the replicon. Nothing here touches ``$GIAP_DB_DIR`` or the network.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from helpers.viz_landscape import write_preflight
from helpers.viz_out_tree import gene_row, make_out, write_schema_table

N_GENES = 25_000
#: share of genes detected with a measured abundance (> 20,000 V2 points)
DETECTED_FRAC = 0.9
MAX_ELEMENTS = 20_000
SEED = 2026


@dataclass(frozen=True)
class StressRun:
    out: Path
    contig: str
    length: int
    n_genes: int
    n_plus: int
    n_minus: int
    n_points: int  # detected, measured, not imputed, not partial


def stress_out(tmp: Path, seed: int = SEED, n_genes: int = N_GENES) -> StressRun:
    """Write the stress ``OUT/`` under ``tmp`` and return its counts."""
    from giapsim.layout import LayoutConfig, simulate_layout

    lay = simulate_layout(
        LayoutConfig(n_genes_chromosome=n_genes, gene_len_median_aa=220.0, regions=()),
        np.random.default_rng(seed),
    )
    rep = lay.replicons.iloc[0]
    contig, length = str(rep["contig"]), int(rep["length"])
    genes = lay.genes
    rng = np.random.default_rng(seed + 1)
    det = rng.random(len(genes)) < DETECTED_FRAC
    log2 = rng.normal(20.0, 2.0, size=len(genes))
    rows = []
    for i, g in enumerate(genes.itertuples(index=False)):
        rows.append(
            gene_row(
                locus_tag=str(g.locus_tag),
                contig=contig,
                start=int(g.start) + 1,  # layout frames are 0-based half-open
                end=int(g.end),
                strand=str(g.strand),
                feature_type="cds",
                product=None if g.product is None else str(g.product),
                detection_status="detected" if det[i] else "not_detected",
                detected=bool(det[i]),
                partial=False,
                crosses_origin=False,
                replicon_type="chromosome",
                abundance_log2=float(log2[i]) if det[i] else None,
                abundance_raw=float(2.0 ** log2[i]) if det[i] else None,
                abundance_imputed=False,
                abundance_type="ibaq" if det[i] else "none",
                cog_category="S",
            )
        )
    out = make_out(tmp)
    write_schema_table(out / "tables" / "genes.tsv", "genes", rows)
    write_preflight(out, [(contig, length, "circular", "chromosome")])
    strands = genes["strand"].astype(str)
    return StressRun(
        out=out,
        contig=contig,
        length=length,
        n_genes=len(genes),
        n_plus=int((strands == "+").sum()),
        n_minus=int((strands == "-").sum()),
        n_points=int(det.sum()),
    )
