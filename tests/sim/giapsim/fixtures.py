"""M4 work directories for the integration, golden and property tests (issue M4-20).

``build_m4_fixture`` writes the M4-06 contract tables C1-C14 of one synthetic genome into a work
directory, in both layouts the M4 steps read: flat (``giap.io.work_tables.default_paths``, used by
``run_context_stats`` and ``run_sensitivity``) and with the three M3 TSVs also under ``mapping/``
(``giap.steps.gene_table.work_inputs``). The result is a pure function of ``(size, seed, evidence,
layout)``, so the whole M4 chain can be run on it twice and compared byte for byte (§11.6 TE7).

Sizes (spec "Interface"):

* ``sim50k``  -- the M0 50-kb genome ``simulate(SimConfig(), seed)`` as an M4-07 layout
  (``layout_from_simulation``, ~40-50 genes depending on the seed) plus synthetic regions: one
  prophage block and one scattered AMR set, both sized as a fraction of the *actual* gene count
  (M4-20 lead decision 5), so every M4 output has rows. This is the §11.4 golden fixture.
* ``sim1000`` / ``sim4500`` -- sequence-free ``simulate_layout`` genomes of 1,000 / 4,500
  chromosomal genes with the default region set; ``layout={...}`` passes ``LayoutConfig`` fields
  (``linear_contigs``, ``plasmids``, ...) through (lead decision 4).

The 50-kb genome is small: with ~29 genes carrying peptide evidence the picked gene q-value cannot
fall below 1/29 = 0.0345, so a run at the default ``stats.alpha.identification`` = 0.01 detects
nothing and the golden tables are recorded at 0.05 (lead decision 1). That is a property of the
fixture, not a defect of the statistics.

Randomness comes only from ``seed`` (three independent ``SeedSequence`` children for the layout,
the expression truth and the evidence). Nothing here imports ``giap``.
"""

from __future__ import annotations

import shutil
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from giapsim.contract_tables import write_contract_tables
from giapsim.evidence import EvidenceConfig, simulate_evidence, write_evidence_tables
from giapsim.expression import ExpressionConfig, simulate_expression
from giapsim.layout import (
    Layout,
    LayoutConfig,
    RegionSpec,
    layout_from_simulation,
    plant_regions,
    simulate_layout,
)

__all__ = ["M3_TSVS", "SIZES", "build_m4_fixture", "fixture_layout"]

#: fixture sizes of the M4-20 spec
SIZES: tuple[str, ...] = ("sim50k", "sim1000", "sim4500")
#: chromosomal genes of the sequence-free sizes (``sim50k`` gets its genes from the M0 simulator)
_N_GENES: Mapping[str, int] = {"sim1000": 1000, "sim4500": 4500}
#: the M3 tables ``gene_table`` reads from ``mapping/`` and the other steps from ``work_dir``
M3_TSVS: tuple[str, ...] = ("peptide_hits.tsv", "protein_hits.tsv", "unmapped.tsv")
#: ``sim50k`` planted regions as a share of the gene count, with floors (lead decision 5)
_PROPHAGE_SHARE, _PROPHAGE_MIN = 0.15, 3
_AMR_SHARE, _AMR_MIN = 0.06, 2


def _streams(seed: int) -> tuple[np.random.Generator, ...]:
    """Three independent generators from ``seed`` (layout, expression, evidence)."""
    return tuple(np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(3))


def _planted(n_genes: int) -> tuple[RegionSpec, ...]:
    """One prophage block and one scattered AMR set sized for ``n_genes`` (decision 5)."""
    return (
        RegionSpec("prophage", max(_PROPHAGE_MIN, round(_PROPHAGE_SHARE * n_genes))),
        RegionSpec("amr", max(_AMR_MIN, round(_AMR_SHARE * n_genes)), scattered=True),
    )


def _sim50k_layout(seed: int, rng: np.random.Generator) -> Layout:
    """The M0 50-kb genome as an M4-07 layout with the planted prophage and AMR regions."""
    from giapsim.simulate import SimConfig, simulate

    with tempfile.TemporaryDirectory(prefix="giapsim-m4-") as tmp:
        layout = layout_from_simulation(simulate(SimConfig(), seed, Path(tmp)))
    return plant_regions(layout, _planted(len(layout.genes)), rng)


def fixture_layout(
    seed: int, size: str = "sim1000", layout: Mapping[str, Any] | None = None
) -> Layout:
    """The layout of the fixture ``(size, seed, layout)`` (``layout``: ``LayoutConfig`` fields)."""
    if size not in SIZES:
        raise ValueError(f"unknown fixture size {size!r}; expected one of {', '.join(SIZES)}")
    rng = _streams(seed)[0]
    if size == "sim50k":
        if layout:
            raise ValueError("layout= is for the sequence-free sizes; sim50k comes from the M0 "
                             "simulator (use sim1000 for linear contigs or plasmids)")  # fmt: skip
        return _sim50k_layout(seed, rng)
    cfg = LayoutConfig(n_genes_chromosome=_N_GENES[size], **dict(layout or {}))
    return simulate_layout(cfg, rng)


def _write_id_map(path: Path, locus_tags: Sequence[str]) -> None:
    """C14 ``id_map.tsv``: the simulator's annotation IDs are the locus tags themselves."""
    lines = ["source_id\tlocus_tag", *(f"{t}\t{t}" for t in locus_tags)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_m4_fixture(
    work_dir: Path,
    *,
    seed: int,
    size: str = "sim1000",
    evidence: EvidenceConfig | None = None,
    layout: Mapping[str, Any] | None = None,
) -> None:
    """Write the M4-06 contract tables of the fixture ``(size, seed, evidence, layout)``.

    ``work_dir`` is created if needed and receives C1-C8 (``replicons``, ``genes``, ``xrefs``,
    ``annotation_meta.json``, ``tus``, ``layers``, ``regions``, ``region_genes``), C15
    ``features``, C9 ``peptides.norm.parquet`` (C10 ``abundance.norm.parquet`` only with
    ``EvidenceConfig.abundance_table``), C11-C13 under ``mapping/`` and as flat copies, and C14
    ``id_map.tsv``.
    """
    _, rng_expr, rng_ev = _streams(seed)
    lay = fixture_layout(seed, size, layout)
    expr = simulate_expression(lay, ExpressionConfig(), rng_expr)
    ev = simulate_evidence(lay, expr, evidence or EvidenceConfig(), rng_ev)
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    write_contract_tables(lay, work_dir)
    write_evidence_tables(ev, work_dir)
    for name in M3_TSVS:
        shutil.copyfile(work_dir / "mapping" / name, work_dir / name)
    _write_id_map(work_dir / "id_map.tsv", lay.genes["locus_tag"].astype(str).tolist())
