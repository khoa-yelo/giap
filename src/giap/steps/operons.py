"""Operon steps (GIAP_SPEC.md §6.2, §7.1.3, §7.4; M2-12, M6-08): ``uniop_input``, ``uniop``
(ext), ``operon_pairs`` and ``operons`` (Module 4).

``uniop_input`` writes the Prodigal-style GFF and its ID map; the external ``uniop`` step runs
UniOP (never imported, ADR-0008). ``operon_pairs`` merges ``uniop.pred`` (row numbers resolved
through ``uniop_cds``) with the distance model (parameters from ``calibration.yaml``), builds TUs
and their upstream FASTAs and records the pair source in ``operon_meta.json``. With
``annotation.precomputed.operon_pairs`` the layer is read directly and its ``source`` kept.
Work-directory tables are 0-based half-open; ``tus.parquet`` ``genes`` is a native list column.

``operons`` chains M6-04 ``pair_frame``, M6-05 ``discordance`` (pass 1: no detection-limit model),
M6-06 ``missed_operons`` and M6-07 ``merge_tests`` / ``flag_pairs`` / ``operons_table`` and writes
``tables/pairs.tsv``, ``operons.tsv`` and ``operon_flags.tsv``. Each absent evidence layer
(``promoter_windows`` → ``promoters``, ``terminators``, ``tx_rates``) gives NA evidence and the
WARNING run flag ``layer_missing:<layer>``; a precomputed layer arrives as its M2-02 TSV under the
same key. The §6.7 event ``operons.missed_operons`` records π₀ and the decoy bin counts (M6-08
decision 6). ``peptides``, ``peptide_hits``, ``peptide_runs`` and ``detection_summary`` are for
M6-02 (pass 2) and are not read yet. Without abundance every test is ``none``.
"""

from __future__ import annotations

import json
from typing import Any

import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]

from giap.annotate.layers import assign_tus, genes_of, replicons_of, terminators_of, tus_of
from giap.annotate.operons import DISTANCE_FLAG, build_pairs, build_tus, write_tu_regions
from giap.annotate.precomputed import load_precomputed
from giap.annotate.uniop_adapter import PairProb, parse_uniop_pred, write_prodigal_gff
from giap.config import GiapConfig
from giap.data import calibration
from giap.io.fasta import read_genome_fasta
from giap.io.genes_tsv import read_genes_parquet
from giap.io.operon_tables import write_operon_flags_tsv, write_operons_tsv, write_pairs_tsv
from giap.io.parquet import write_parquet_deterministic
from giap.io.tsv import read_tsv, write_tsv
from giap.model.coords import Replicon
from giap.operons.discordance import NO_DL_MODEL, discordance
from giap.operons.flags import LAYER_MISSING, flag_pairs, merge_tests
from giap.operons.missed import missed_operons
from giap.operons.pairs import RNA_TYPES, Evidence, pair_frame
from giap.operons.tu_table import operons_table
from giap.steps.common import one_path as _path
from giap.steps.common import read_genes as _genes
from giap.steps.common import read_replicons as _replicons
from giap.steps.registry import StepContext, external, step

_STR, _INT, _F64 = pa.string(), pa.int64(), pa.float64()
_PAIR_COLS = [("gene_a", _STR), ("gene_b", _STR), ("uniop_p", _F64), ("source", _STR)]
_TU_COLS = [
    ("tu_id", _STR), ("contig", _STR), ("start", _INT), ("end", _INT), ("strand", _STR),
    ("genes", pa.list_(_STR)), ("upstream_start", _INT), ("upstream_end", _INT),
]  # fmt: skip


@step(
    "uniop_input",
    inputs=["genes", "replicons"],
    outputs=["uniop_cds", "uniop_map"],
    precomputed="operon_pairs",
)
def uniop_input(ctx: StepContext) -> None:
    """genes → Prodigal-style GFF (``uniop_cds``) and ``prodigal_id → locus_tag`` map."""
    reps = _replicons(_path(ctx.inputs["replicons"]))
    genes = _genes(_path(ctx.inputs["genes"]))
    id_map = write_prodigal_gff(ctx.outputs["uniop_cds"], genes, list(reps.values()))
    write_tsv(ctx.outputs["uniop_map"], ["prodigal_id", "locus_tag"], id_map)


external("uniop", inputs=["uniop_cds"], outputs=["uniop_pred"], precomputed="operon_pairs")


def _uniop_pairs(ctx: StepContext, reps: dict[str, Replicon]) -> list[PairProb] | None:
    layer = ctx.cfg.annotation.precomputed.operon_pairs
    if layer is not None:
        df = load_precomputed("operon_pairs", layer, reps)
        return [
            PairProb(str(a), str(b), float(p), src)
            for a, b, p, src in df[["gene_a", "gene_b", "uniop_p", "source"]].itertuples(
                index=False
            )
        ]
    pred = ctx.inputs.get("uniop_pred")
    if pred is None:
        return None
    _, rows = read_tsv(_path(ctx.inputs["uniop_map"]))
    id_map = {str(r[0]): str(r[1]) for r in rows}
    return parse_uniop_pred(_path(pred), _path(ctx.inputs["uniop_cds"]), id_map)


@step(
    "operon_pairs",
    inputs=["genes", "uniop_map", "uniop_cds", "replicons", "genome"],
    optional_inputs=["uniop_pred"],
    outputs=["pairs_raw", "tus", "tu_upstream", "operon_meta"],
)
def operon_pairs(ctx: StepContext) -> None:
    """UniOP / distance-model pairs, TUs (threshold ``annotation.uniop.pair_threshold``)."""
    reps = _replicons(_path(ctx.inputs["replicons"]))
    genes = _genes(_path(ctx.inputs["genes"]))
    model = calibration.load()["operon_distance"]
    uniop = _uniop_pairs(ctx, reps)
    pairs, flags = build_pairs(genes, reps, uniop, float(model["d0"]), float(model["s"]))
    if ctx.cfg.annotation.precomputed.operon_pairs is not None:
        meta = {"operon_source": "precomputed", "reason": "annotation.precomputed.operon_pairs"}
    elif DISTANCE_FLAG in flags:
        why = "no uniop.pred" if "uniop_pred" not in ctx.inputs else "uniop.pred has no pair"
        meta = {"operon_source": "distance_model", "reason": f"UniOP failed ({why})"}
        ctx.log.warning(
            "operon_pairs: %s; distance model used, run flag %s",
            why,
            DISTANCE_FLAG,
            extra={"event": "operon_pairs.run_flag", "data": {"flag": DISTANCE_FLAG}},
        )
    else:
        meta = {"operon_source": "uniop", "reason": "UniOP probabilities (pin 3ed14dfc)"}
    ctx.outputs["operon_meta"].write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    ann = ctx.cfg.annotation
    up_bp = (ann.tu_upstream_bp[0], ann.tu_upstream_bp[1])
    tus = build_tus(genes, pairs, ann.uniop.pair_threshold, reps, up_bp)
    out = ctx.outputs
    pair_rows = [(p.gene_a, p.gene_b, p.uniop_p, p.source) for p in pairs]
    write_parquet_deterministic(out["pairs_raw"], _PAIR_COLS, pair_rows)
    tu_rows = [
        (t.tu_id, t.contig, t.interval.start, t.interval.end, t.strand, list(t.genes),
         t.upstream.start, t.upstream.end)
        for t in tus
    ]  # fmt: skip
    write_parquet_deterministic(out["tus"], _TU_COLS, tu_rows)
    genome = {r.id: r.sequence for r in read_genome_fasta(_path(ctx.inputs["genome"]))}
    write_tu_regions(out["tu_upstream"], tus, genome)


#: optional evidence key -> layer name of ``layer_missing:<layer>`` (as ``context_layers``)
_LAYERS = {"promoter_windows": "promoters", "terminators": "terminators", "tx_rates": "tx_rates"}
_MISSED_EVENT = "operons.missed_operons"


def run_operons(
    gene_table: pd.DataFrame,
    genes: pd.DataFrame,
    replicons: pd.DataFrame,
    tus: pd.DataFrame,
    pairs_raw: pd.DataFrame,
    features: pd.DataFrame,
    layers: dict[str, pd.DataFrame],
    cfg: GiapConfig,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """The ``pairs``, ``operons`` and ``operon_flags`` frames and the decoy ``bin_counts``
    (lead decision 4's chain). ``layers`` maps each present evidence key (``promoter_windows``,
    ``terminators``, ``tx_rates``, in the producer shapes) to its frame. Pure."""
    terms = layers.get("terminators")
    ev = Evidence(
        windows=layers.get("promoter_windows"),
        terminators=None if terms is None else tuple(terminators_of(terms)),
        tx=layers.get("tx_rates"),
        rna=features[features["type"].isin(RNA_TYPES)].reset_index(drop=True),
    )
    pairs = pair_frame(gene_table, genes_of(genes), replicons_of(replicons), pairs_raw, ev, cfg)
    d = discordance(pairs, gene_table, cfg, None)
    m = missed_operons(pairs, cfg, genes=gene_table)
    missing = [name for key, name in _LAYERS.items() if key not in layers]
    merged = merge_tests(pairs, d.rows, m.rows, gene_table, missing_layers=missing)
    flagged = flag_pairs(merged, cfg)
    return merged, operons_table(tus, merged, gene_table, flagged), flagged, m.bin_counts


def _layer(ctx: StepContext, key: str, reps: dict[str, Replicon], tus: pd.DataFrame) -> Any:
    """The evidence layer under ``key`` (producer Parquet, or the M2-02 TSV when precomputed)."""
    pre = "promoters" if key == "promoter_windows" else key
    path = _path(ctx.inputs[key])
    if getattr(ctx.cfg.annotation.precomputed, pre) is None:
        return pd.read_parquet(path)
    df = load_precomputed(pre, path, reps)
    return assign_tus(df, tus_of(tus), reps) if key == "tx_rates" else df


def _plain(v: Any) -> Any:
    return v.item() if hasattr(v, "item") else v


@step(
    "operons",
    inputs=["gene_table", "genes", "replicons", "tus", "pairs_raw", "features"],
    optional_inputs=["promoter_windows", "terminators", "tx_rates", "peptides", "peptide_hits",
                     "peptide_runs", "detection_summary"],
    outputs=["pairs", "operons", "operon_flags"],
)  # fmt: skip
def operons_step(ctx: StepContext) -> None:
    """Pair table, discordance and missed-operon tests, flags and the TU table (Module 4)."""
    ins = {k: _path(v) for k, v in ctx.inputs.items()}
    rd = pd.read_parquet
    replicons, tus = rd(ins["replicons"]), rd(ins["tus"])
    reps = replicons_of(replicons)
    layers = {k: _layer(ctx, k, reps, tus) for k in _LAYERS if k in ins}
    for key, name in _LAYERS.items():
        if key not in layers:
            ctx.log.warning("%s%s: layer not available; its pair evidence is NA",
                            LAYER_MISSING, name, extra={"event": "operons.layer_missing", "data":
                                {"layer": name, "flag": f"{LAYER_MISSING}{name}"}})  # fmt: skip
    ctx.log.info("%s: pass 1 has no detection-limit model; one-detected pairs are untested",
                 NO_DL_MODEL,
                 extra={"event": "operons.run_flag", "data": {"flag": NO_DL_MODEL}})  # fmt: skip
    gt = read_genes_parquet(ins["gene_table"])
    pairs, ops, flagged, bins = run_operons(
        gt, rd(ins["genes"]), replicons, tus, rd(ins["pairs_raw"]), rd(ins["features"]), layers,
        ctx.cfg,
    )  # fmt: skip
    data = {
        "pi0": float(ctx.cfg.operons.pi0),
        "bin_counts": [{k: _plain(v) for k, v in r.items()} for r in bins.to_dict("records")],
    }
    ctx.log.info("operons: decoy bin counts recorded (pi0 = %s)", data["pi0"],
                 extra={"event": _MISSED_EVENT, "data": data})  # fmt: skip
    write_pairs_tsv(ctx.outputs["pairs"], pairs)
    write_operons_tsv(ctx.outputs["operons"], ops)
    write_operon_flags_tsv(ctx.outputs["operon_flags"], flagged)
