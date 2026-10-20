"""Step ``detectability`` (GIAP_SPEC.md §6.2 "gene table → p_detect (cross-fitted)", §7.3.3;
docs/specs/M4-15-gene-table-assembly.md): adds the M4-13 cross-fitted ``p_detect`` to
``genes_pre`` and writes the ``detectability`` rows of ``detection_summary``.

Eligible genes are the non-partial ones (§7.3); their ``detected`` value is the M4-12 call. The
model draws only from the step's own stream ``giap.rng.for_step("detectability", seed)``. With a
constant model every eligible gene gets the flag ``detectability_constant``.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]

from giap.config import GiapConfig
from giap.genetable.assemble import GENES_COLUMNS, gene_flags
from giap.genetable.covariates import COVARIATE_COLUMNS
from giap.genetable.detectability import COLUMNS as DIAG_COLUMNS
from giap.genetable.detectability import DetectabilityResult, fit_detectability
from giap.io.genes_tsv import read_genes_parquet, write_genes_parquet
from giap.io.parquet import write_parquet_deterministic
from giap.steps.registry import StepContext, run_step, step

__all__ = ["add_p_detect", "run_detectability"]

_DIAG_TYPES = (pa.string(),) * 4 + (pa.float64(),) * 3


def add_p_detect(pre: pd.DataFrame, res: DetectabilityResult) -> pd.DataFrame:
    """``pre`` with ``p_detect`` in its §9.3 place and ``detectability_constant`` when constant."""
    out = pre.copy()
    tags = out["locus_tag"].astype(str)
    p = res.p_detect.set_axis(res.p_detect.index.astype(str))
    out["p_detect"] = p.reindex(tags).to_numpy(dtype="float64")
    if res.constant:
        out["flags"] = [
            gene_flags(
                {"genes_pre": list(f), "detectability": [] if part else ["detectability_constant"]}
            )
            for f, part in zip(out["flags"], out["partial"].astype(bool), strict=True)
        ]
    return out[[c for c in GENES_COLUMNS if c in out.columns]]


@step(
    "detectability",
    inputs=["genes_pre"],
    outputs=["genes_final", "detectability_summary"],
)
def detectability(ctx: StepContext) -> None:
    """``genes_pre`` → gene table with ``p_detect`` and the diagnostics rows."""
    pre = read_genes_parquet(Path(str(ctx.inputs["genes_pre"])))
    idx = pd.Index(pre["locus_tag"].astype(str))
    detected = pd.Series(pre["detected"].astype(bool).to_numpy(), index=idx, dtype="boolean")
    detected[pre["partial"].astype(bool).to_numpy()] = pd.NA
    cov = pre[["locus_tag", *COVARIATE_COLUMNS]].astype(dict.fromkeys(COVARIATE_COLUMNS, "float64"))
    tus = pd.Series(pre["tu_id"].to_numpy(), index=idx)
    res = fit_detectability(detected, cov, tus, ctx.cfg.detectability, ctx.rng)
    write_genes_parquet(ctx.outputs["genes_final"], add_p_detect(pre, res))
    diag = res.diagnostics
    rows = [
        tuple(
            (None if pd.isna(v) else str(v)) if i < 4 else (None if pd.isna(v) else float(v))
            for i, v in enumerate(r)
        )
        for r in diag[DIAG_COLUMNS].itertuples(index=False, name=None)
    ]
    write_parquet_deterministic(
        ctx.outputs["detectability_summary"], list(zip(DIAG_COLUMNS, _DIAG_TYPES, strict=True)),
        rows,
    )  # fmt: skip
    ctx.log.info("detectability: constant=%s, %d diagnostics rows", res.constant, len(rows))


def run_detectability(work_dir: Path, cfg: GiapConfig) -> Path:
    """Run ``detectability`` on ``work_dir/genetable/genes_pre.parquet``; returns
    ``work_dir/genetable/genes.parquet`` (diagnostics in ``detection_summary_parts/``)."""
    gt = work_dir / "genetable"
    out = gt / "genes.parquet"
    outputs = {
        "genes_final": out,
        "detectability_summary": gt / "detection_summary_parts" / "detectability.parquet",
    }
    run_step("detectability", cfg, {"genes_pre": gt / "genes_pre.parquet"}, outputs)
    return out
