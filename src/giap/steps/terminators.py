"""Terminator steps (GIAP_SPEC.md §6.2, §7.1.5; M2-15): ``transterm_input``, the external
``terminatornet`` / ``transtermhp`` / ``bactermfinder`` (opt-in) and ``terminator_layer``.

``terminators.parquet`` (0-based half-open) holds the primary tool's calls (+ BacTermFinder when
enabled), never the cross-check tool's. ``terminators_crosscheck.tsv`` is header-only when no
cross-check runs (``--terminators transtermhp`` or ``cross_check_transtermhp: false``). With
``annotation.precomputed.terminators`` all these steps are replaced by the layer.
"""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]

from giap.annotate.terminators import (
    Terminator,
    crosscheck,
    parse_bactermfinder,
    parse_terminatornet,
    parse_transterm_tt,
    write_transterm_coords,
)
from giap.io.parquet import write_parquet_deterministic
from giap.io.tsv import write_tsv
from giap.model.coords import Interval, Replicon, interval_sort_key
from giap.model.genes import Gene
from giap.steps.registry import InputPath, StepContext, external, step

_STR, _INT, _F64 = pa.string(), pa.int64(), pa.float64()
_LAYER_COLS = [
    ("contig", _STR), ("start", _INT), ("end", _INT), ("strand", _STR), ("confidence", _F64),
    ("source", _STR), ("hairpin_dg", _F64),
]  # fmt: skip
_XC_COLS = ["contig", "n_terminatornet", "n_transtermhp", "n_matched", "share_matched"]
_PRE = "terminators"


def _path(p: InputPath) -> Path:
    if isinstance(p, list):
        raise ValueError(f"expected one file, got {len(p)}")
    return p


def _replicons(path: Path) -> dict[str, Replicon]:
    rows = pq.read_table(path).to_pylist()
    return {r["contig"]: Replicon(r["contig"], r["length"], r["topology"], r["kind"]) for r in rows}


def _genes(path: Path) -> list[Gene]:
    cols = ["locus_tag", "feature_type", "contig", "start", "end", "strand", "aa", "start_codon"]
    return [
        Gene(
            locus_tag=r["locus_tag"],
            feature_type=r["feature_type"],
            interval=Interval(r["contig"], r["start"], r["end"], r["strand"]),
            aa=r["aa"] or "",
            nt=None,
            start_codon=r["start_codon"],
        )
        for r in pq.read_table(path, columns=cols).to_pylist()
    ]


external("terminatornet", inputs=["genome"], outputs=["terminatornet"], precomputed=_PRE)


@step("transterm_input", inputs=["genes"], outputs=["genes_coords"], precomputed=_PRE)
def transterm_input(ctx: StepContext) -> None:
    """genes → TransTermHP ``genes.coords`` (CDS, sORFs and pseudogenes)."""
    write_transterm_coords(ctx.outputs["genes_coords"], _genes(_path(ctx.inputs["genes"])))


external(
    "transtermhp", inputs=["genome", "genes_coords"], outputs=["transtermhp"], precomputed=_PRE
)
external(
    "bactermfinder",
    inputs=["genome"],
    outputs=["bactermfinder"],
    module="opt_in",
    precomputed=_PRE,
)


def _bactermfinder(p: Path, reps: dict[str, Replicon]) -> list[Terminator]:
    """A directory without a CSV → no calls; exactly one CSV (the ``*_mean.csv`` table) →
    parse it; several CSVs → ValueError naming them (no arbitrary pick)."""
    if p.is_dir():
        files = sorted(f for f in p.iterdir() if f.is_file() and f.suffix == ".csv")
        if not files:
            return []
        if len(files) > 1:
            names = ", ".join(f.name for f in files)
            raise ValueError(f"several BacTermFinder tables in {p}: {names}")
        p = files[0]
    return parse_bactermfinder(p, reps)


@step(
    "terminator_layer",
    inputs=["replicons"],
    optional_inputs=["terminatornet", "transtermhp", "bactermfinder"],
    outputs=["terminators", "terminators_crosscheck"],
    precomputed=_PRE,
)
def terminator_layer(ctx: StepContext) -> None:
    """Primary-tool calls (``annotation.terminators.tool``) + optional BacTermFinder, and the
    TerminatorNet/TransTermHP cross-check table."""
    cfg = ctx.cfg.annotation.terminators
    reps = _replicons(_path(ctx.inputs["replicons"]))
    primary_key = "transtermhp" if cfg.tool == "transtermhp" else "terminatornet"
    if primary_key not in ctx.inputs:
        raise ValueError(f"terminator_layer: tool {cfg.tool} selected but {primary_key} missing")
    tn = tt = None
    if "terminatornet" in ctx.inputs:
        tn = parse_terminatornet(_path(ctx.inputs["terminatornet"]), reps)
    if "transtermhp" in ctx.inputs:
        tt = parse_transterm_tt(_path(ctx.inputs["transtermhp"]), reps)
    primary = tt if cfg.tool == "transtermhp" else tn
    assert primary is not None
    layer = list(primary)
    if cfg.bactermfinder and "bactermfinder" in ctx.inputs:
        layer += _bactermfinder(_path(ctx.inputs["bactermfinder"]), reps)
    rank = {c: k for k, c in enumerate(reps)}
    layer.sort(key=lambda t: (interval_sort_key(t.interval, rank, t.source), t.confidence))
    rows = [
        (t.interval.contig, t.interval.start, t.interval.end, t.interval.strand, t.confidence,
         t.source, t.hairpin_dg)
        for t in layer
    ]  # fmt: skip
    write_parquet_deterministic(ctx.outputs["terminators"], _LAYER_COLS, rows)
    xc_rows: list[tuple[object, ...]] = []
    if cfg.tool == "terminatornet" and cfg.cross_check_transtermhp and tn is not None:
        if tt is None:
            ctx.log.warning("terminator_layer: no TransTermHP output; cross-check table empty")
        else:
            xc_rows = [
                (r.contig, r.n_terminatornet, r.n_transtermhp, r.n_matched, r.share_matched)
                for r in crosscheck(tn, tt, contig_order=list(reps))
            ]
    write_tsv(ctx.outputs["terminators_crosscheck"], _XC_COLS, xc_rows)
