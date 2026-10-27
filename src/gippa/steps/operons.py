"""Operon steps (GIPPA_SPEC.md §6.2, §7.1.3, §7.4; M2-12, M6-08): ``uniop_input``, ``uniop``
(ext), ``operon_pairs`` and ``operons`` (Module 4).

``uniop_input`` writes the Prodigal-style GFF and its ID map; the external ``uniop`` step runs
UniOP (never imported, ADR-0008). ``operon_pairs`` merges ``uniop.pred`` (row numbers resolved
through ``uniop_cds``) with the distance model (parameters from ``calibration.yaml``), builds TUs
and their upstream FASTAs and records the pair source in ``operon_meta.json``. With
``annotation.precomputed.operon_pairs`` the layer is read directly and its ``source`` kept.
Work-directory tables are 0-based half-open; ``tus.parquet`` ``genes`` is a native list column.

``operons`` chains M6-04 ``pair_frame``, M6-05 ``discordance``, M6-06 ``missed_operons`` and
M6-07 ``merge_tests`` / ``flag_pairs`` / ``operons_table`` and writes
``tables/pairs.tsv``, ``operons.tsv`` and ``operon_flags.tsv``; with ``operon_meta`` saying
``distance_model`` every ``both_detected_t`` row of ``pairs.tsv`` carries
``tu_source_distance_model`` (ADR-0061). Each absent evidence layer
(``promoter_windows`` → ``promoters``, ``terminators``, ``tx_rates``) gives NA evidence and the
WARNING run flag ``layer_missing:<layer>``; a precomputed layer arrives as its M2-02 TSV under the
same key. The §6.7 event ``operons.missed_operons`` records π₀ and the decoy bin counts (M6-08
decision 6). With ``peptide_hits`` (C11) the step fits M6-02's detection-limit model
(``operons.detection_limit`` RNG stream) and, when it is usable, scores the one-detected OP pairs
(exactly one gene with >= ``k_min`` observable peptides seen; M6-02 lead decision 9) with
``P(miss)`` from the out-of-fold model of the seen partner's TU (M6-05H; no q, outside the BH
family). Without peptide hits or a usable model those pairs stay ``none`` with the run flag
``no_detection_limit_model`` (and the event stays ``detected_a xor detected_b``). A model in
force that needs replicate covariates the step does not join yet leaves them ``none`` with the
pair flag ``pmiss_covariates_missing`` and a WARNING. ``peptides``, ``peptide_runs`` (the
replicate model) and ``detection_summary`` are not read yet. Without abundance every test is
``none``. The §6.7 event ``operons.discordance_fit`` records the within-operon null in force
(M6-13: the log-linear scale, or the constant fallback and why).
"""

from __future__ import annotations

import json
import logging
import math
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]

from gippa.annotate.layers import assign_tus, genes_of, replicons_of, terminators_of, tus_of
from gippa.annotate.operons import DISTANCE_FLAG, build_pairs, build_tus, write_tu_regions
from gippa.annotate.precomputed import load_precomputed
from gippa.annotate.uniop_adapter import PairProb, parse_uniop_pred, write_prodigal_gff
from gippa.config import GippaConfig
from gippa.data import calibration
from gippa.io.fasta import read_genome_fasta
from gippa.io.genes_tsv import read_genes_parquet
from gippa.io.operon_tables import write_operon_flags_tsv, write_operons_tsv, write_pairs_tsv
from gippa.io.parquet import write_parquet_deterministic
from gippa.io.tsv import read_tsv, write_tsv
from gippa.io.work_tables import _read_tsv_table
from gippa.model.coords import Replicon
from gippa.operons import detection_limit as dl
from gippa.operons.discordance import FALLBACK_WARN, NO_DL_MODEL, SIGMA_MODEL, PMissFn, discordance
from gippa.operons.flags import LAYER_MISSING, add_tu_source_flag, flag_pairs, merge_tests
from gippa.operons.missed import missed_operons
from gippa.operons.pairs import RNA_TYPES, Evidence, pair_frame
from gippa.operons.tu_table import operons_table
from gippa.rng import for_step
from gippa.stats import tdist
from gippa.stats.tdist import ScaledStudentT
from gippa.steps.common import one_path as _path
from gippa.steps.common import operon_source
from gippa.steps.common import read_genes as _genes
from gippa.steps.common import read_replicons as _replicons
from gippa.steps.registry import StepContext, external, step

#: ``operon_meta.json`` reasons of the distance model (M11-14d; ``UniOP failed (...)`` when UniOP
#: ran and produced no pair)
REASON_NOT_INSTALLED = "UniOP not installed"
REASON_DISABLED = "UniOP disabled (annotation.uniop.enabled=false)"
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
    uni = ctx.cfg.annotation.uniop
    precomputed = ctx.cfg.annotation.precomputed.operon_pairs is not None
    off = not precomputed and uni.enabled is False  # a stray uniop.pred is ignored
    uniop = None if off else _uniop_pairs(ctx, reps)
    pairs, flags = build_pairs(genes, reps, uniop, float(model["d0"]), float(model["s"]))
    if precomputed:
        meta = {"operon_source": "precomputed", "reason": "annotation.precomputed.operon_pairs"}
    elif DISTANCE_FLAG in flags:
        if off:
            reason = REASON_NOT_INSTALLED if uni.resolved_from == "auto" else REASON_DISABLED
        else:
            why = "no uniop.pred" if "uniop_pred" not in ctx.inputs else "uniop.pred has no pair"
            reason = f"UniOP failed ({why})"
        meta = {"operon_source": "distance_model", "reason": reason}
        ctx.log.warning(
            "operon_pairs: %s; distance model used, run flag %s",
            reason,
            DISTANCE_FLAG,
            extra={"event": "operon_pairs.run_flag", "data": {"flag": DISTANCE_FLAG}},
        )
        ctx.log.warning(
            "%suniop: operon pairs and TUs come from the distance model",
            LAYER_MISSING,
            extra={"event": "operon_pairs.layer_missing",
                   "data": {"layer": "uniop", "flag": f"{LAYER_MISSING}uniop"}},
        )  # fmt: skip
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
LAYERS = {"promoter_windows": "promoters", "terminators": "terminators", "tx_rates": "tx_rates"}
_MISSED_EVENT = "operons.missed_operons"
DL_STEP = "operons.detection_limit"
#: pair flag: the model in force needs replicate covariates the step does not join (M6-05H, S2)
COV_MISSING = "pmiss_covariates_missing"


def _with_aa(gene_table: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    aa = genes.drop_duplicates("locus_tag").set_index("locus_tag")["aa"]
    return gene_table.assign(aa=gene_table["locus_tag"].map(aa))


def detection_limit(
    gene_table: pd.DataFrame, genes: pd.DataFrame, peptide_hits: pd.DataFrame, cfg: GippaConfig
) -> tuple[dl.DetectionLimitFit, PMissFn | None]:
    """M6-02's fit on the gene table (plus the C2 ``aa``) and C11, and the ``p_miss`` callable of
    :func:`discordance` when the model in force is usable and needs no covariate the step does
    not join (else ``None``). Deterministic: the RNG is
    ``for_step("operons.detection_limit", cfg.seed)``."""
    g = _with_aa(gene_table, genes)
    fit = dl.fit_detection_limit(g, peptide_hits, None, cfg, for_step(DL_STEP, cfg.seed))
    ok = fit.model.usable and not dl.needs_covariates(fit.model)
    return fit, (dl.p_miss_for_genes(fit.model, g, cfg) if ok else None)


def run_operons(
    gene_table: pd.DataFrame,
    genes: pd.DataFrame,
    replicons: pd.DataFrame,
    tus: pd.DataFrame,
    pairs_raw: pd.DataFrame,
    features: pd.DataFrame,
    layers: dict[str, pd.DataFrame],
    cfg: GippaConfig,
    p_miss: PMissFn | None = None,
    *,
    n_seen: Mapping[str, int] | None = None,
    k_min: int = 1,
    no_model_flag: str = NO_DL_MODEL,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, ScaledStudentT]:
    """The ``pairs``, ``operons`` and ``operon_flags`` frames, the decoy ``bin_counts`` and the
    discordance null ``Discordance.scaled`` (lead decision 4's chain; M6-13). ``layers`` maps
    each present evidence key (``promoter_windows``, ``terminators``, ``tx_rates``, in the
    producer shapes) to its frame; ``p_miss`` is the
    detection-limit callable (:func:`detection_limit`), ``None`` without one; ``n_seen``,
    ``k_min`` and ``no_model_flag`` go to :func:`discordance`. Pure."""
    terms = layers.get("terminators")
    ev = Evidence(
        windows=layers.get("promoter_windows"),
        terminators=None if terms is None else tuple(terminators_of(terms)),
        tx=layers.get("tx_rates"),
        rna=features[features["type"].isin(RNA_TYPES)].reset_index(drop=True),
    )
    pairs = pair_frame(gene_table, genes_of(genes), replicons_of(replicons), pairs_raw, ev, cfg)
    d = discordance(pairs, gene_table, cfg, p_miss, n_seen=n_seen, k_min=k_min,
                    no_model_flag=no_model_flag)  # fmt: skip
    m = missed_operons(pairs, cfg, genes=gene_table)
    missing = [name for key, name in LAYERS.items() if key not in layers]
    merged = merge_tests(pairs, d.rows, m.rows, gene_table, missing_layers=missing)
    flagged = flag_pairs(merged, cfg)
    return merged, operons_table(tus, merged, gene_table, flagged), flagged, m.bin_counts, d.scaled


def read_layer(
    cfg: GippaConfig, key: str, path: Path, reps: dict[str, Replicon], tus: pd.DataFrame
) -> Any:
    """The evidence layer under ``key`` (producer Parquet, or the M2-02 TSV when precomputed)."""
    pre = "promoters" if key == "promoter_windows" else key
    if getattr(cfg.annotation.precomputed, pre) is None:
        return pd.read_parquet(path)
    df = load_precomputed(pre, path, reps)
    return assign_tus(df, tus_of(tus), reps) if key == "tx_rates" else df


def _plain(v: Any) -> Any:
    """numpy scalars → Python; NaN and ±∞ → None so the event line is strict JSON (M6-08)."""
    v = v.item() if hasattr(v, "item") else v
    return None if isinstance(v, float) and not math.isfinite(v) else v


def discordance_fit_data(s: ScaledStudentT) -> dict[str, Any]:
    """The ``operons.discordance_fit`` payload (M6-13; strict JSON)."""
    return {k: _plain(v) for k, v in {
        "sigma_model": SIGMA_MODEL, "sigma_model_used": s.model, "reason": s.reason, "n_fit": s.n,
        "loc_log2": s.loc, "scale0_log2": s.scale0, "gamma_per_sd": s.gamma,
        "gamma_per_log2": s.gamma / s.sd if s.model != "none" and s.sd > 0 else math.nan,
        "df": s.df, "df_at_bound": s.model != "none" and s.df in tdist.DF_BOUNDS,
        "gamma_at_bound": s.model == "loglinear" and s.gamma in tdist.GAMMA_BOUNDS,
        "abundance_center_log2": s.center, "abundance_sd_log2": s.sd, "z_min": s.z_min,
        "z_max": s.z_max, "lr_stat": s.lr_stat, "lr_p": s.lr_p,
    }.items()}  # fmt: skip


@step(
    "operons",
    inputs=["gene_table", "genes", "replicons", "tus", "pairs_raw", "features"],
    optional_inputs=["promoter_windows", "terminators", "tx_rates", "peptides", "peptide_hits",
                     "peptide_runs", "detection_summary", "operon_meta"],
    outputs=["pairs", "operons", "operon_flags"],
)  # fmt: skip
def operons_step(ctx: StepContext) -> None:
    """Pair table, discordance and missed-operon tests, flags and the TU table (Module 4)."""
    ins = {k: _path(v) for k, v in ctx.inputs.items()}
    rd = pd.read_parquet
    replicons, tus = rd(ins["replicons"]), rd(ins["tus"])
    reps = replicons_of(replicons)
    layers = {k: read_layer(ctx.cfg, k, ins[k], reps, tus) for k in LAYERS if k in ins}
    for key, name in LAYERS.items():
        if key not in layers:
            ctx.log.warning("%s%s: layer not available; its pair evidence is NA",
                            LAYER_MISSING, name, extra={"event": "operons.layer_missing", "data":
                                {"layer": name, "flag": f"{LAYER_MISSING}{name}"}})  # fmt: skip
    gt, genes = read_genes_parquet(ins["gene_table"]), rd(ins["genes"])
    p_miss: PMissFn | None = None
    seen_kw: dict[str, Any] = {}
    if "peptide_hits" in ins:
        hits = _read_tsv_table("peptide_hits", ins["peptide_hits"])
        fit, p_miss = detection_limit(gt, genes, hits, ctx.cfg)
        m = fit.model
        dl_data = {"kind": m.kind, "usable": m.usable, "n_fit": m.n_fit,
                "auc": m.auc if math.isfinite(m.auc) else None, "reason": m.reason}  # fmt: skip
        ctx.log.info("operons: detection-limit model %s (usable=%s, n_fit=%d, auc=%.4g%s)",
                     m.kind, m.usable, m.n_fit, m.auc, f", {m.reason}" if m.reason else "",
                     extra={"event": DL_STEP, "data": dl_data})  # fmt: skip
        if m.usable:  # lead decision: the k_min event only with a usable model
            seen_kw = {"n_seen": dl.observable_peptides_seen(_with_aa(gt, genes), hits, ctx.cfg),
                       "k_min": m.k_min}  # fmt: skip
            if p_miss is None:  # needs covariates (S2)
                seen_kw["no_model_flag"] = COV_MISSING
                covs = sorted(k for k in m.coef if k not in ("intercept", "abundance_log2"))
                ctx.log.warning("%s: the %s detection-limit model needs covariates (%s) that are "
                                "not joined yet; one-detected pairs are untested", COV_MISSING,
                                m.kind, ", ".join(covs),
                                extra={"event": "operons.pmiss_covariates_missing", "data":
                                       {"flag": COV_MISSING, "covariates": covs}})  # fmt: skip
    if p_miss is None and "no_model_flag" not in seen_kw:
        why = "no peptide_hits" if "peptide_hits" not in ins else "model not usable"
        ctx.log.info("%s: %s; one-detected pairs are untested", NO_DL_MODEL, why,
                     extra={"event": "operons.run_flag",
                            "data": {"flag": NO_DL_MODEL}})  # fmt: skip
    pairs, ops, flagged, bins, scaled = run_operons(
        gt, genes, replicons, tus, rd(ins["pairs_raw"]), rd(ins["features"]), layers, ctx.cfg,
        p_miss, **seen_kw,
    )  # fmt: skip
    dfit = discordance_fit_data(scaled)
    ctx.log.log(logging.WARNING if scaled.reason in FALLBACK_WARN else logging.INFO,
                "operons: within-operon discordance null",
                extra={"event": "operons.discordance_fit", "data": dfit})  # fmt: skip
    data = {
        "pi0": float(ctx.cfg.operons.pi0),
        "bin_counts": [{k: _plain(v) for k, v in r.items()} for r in bins.to_dict("records")],
    }
    ctx.log.info("operons: decoy bin counts recorded (pi0 = %s)", data["pi0"],
                 extra={"event": _MISSED_EVENT, "data": data})  # fmt: skip
    # ADR-0061 decision 1: both_detected_t rows are conservative under distance-model TUs
    pairs = add_tu_source_flag(pairs, "pairs", operon_source(ctx.inputs.get("operon_meta")))
    write_pairs_tsv(ctx.outputs["pairs"], pairs)
    write_operons_tsv(ctx.outputs["operons"], ops)
    write_operon_flags_tsv(ctx.outputs["operon_flags"], flagged)
