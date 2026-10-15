"""Mapping steps (GIAP_SPEC.md §6.1, §6.2, §7.2; M3-17).

Six-frame index ─► exact peptides ─► tolerant (+ decoys) ─► placement FDR ─► loci; proteins
``> 40 aa`` ─► miniprot ─► tblastn (leftovers). The pure functions of M3-02 … M3-16 do the work;
these steps only read and write the work files. Keys: roots ``genome``, ``peptides``,
``proteins_long``, ``reference_peptides``; ``genes``/``features``/``replicons`` from
``parse_bakta``; intermediates ``sixframe`` (``.npz``), ``hits_exact``, ``hits_tolerant``,
``placement``, ``mapping_fdr`` (Parquet), ``miniprot_gff``, ``proteins_leftover``,
``tblastn_tsv``; ``loci`` writes the seven §7.2.9 files into the directory ``mapping``.

Queries (Q2, Q3): non-decoy rows with ``q_value <= input.q_value_max`` or no q-value and length in
``mapping.exact.[min_len, max_len]``; a row source ``auto`` (or missing) is ``input.peptide_source``
or, when that is ``auto`` too, ``database`` if the table has any q-value or PEP, else ``denovo``.
"""

from __future__ import annotations

import hashlib
from collections.abc import Collection, Iterable, Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq  # type: ignore[import-untyped]

from giap.config import GiapConfig
from giap.io.fasta import read_genome_fasta, write_fasta
from giap.io.peptides import read_peptide_table, validate_peptides
from giap.io.protein_fasta import read_protein_fasta
from giap.mapping import denovo_precision as dnp
from giap.mapping.chance import residue_freqs
from giap.mapping.contaminants import match_contaminant_peptides, unmapped_table
from giap.mapping.decoys import DecoyType, build_decoys
from giap.mapping.exact import match_exact
from giap.mapping.frameshift import match_frameshift
from giap.mapping.hits import hits_frame
from giap.mapping.loci import annotate_peptide_hits, annotate_protein_hits
from giap.mapping.miniprot import ProteinAlignment, parse_miniprot_gff
from giap.mapping.outputs import write_mapping_outputs
from giap.mapping.placement_fdr import best_costs, placement
from giap.mapping.protein_assign import PROTEIN_HITS_COLUMNS, assign_alignments
from giap.mapping.sixframe import Frame, SixFrameIndex, build_sixframe, load_index, save_index
from giap.mapping.tblastn import merge_hsps, parse_tblastn
from giap.mapping.tolerant import match_tolerant
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene, GeneModel, NcFeature
from giap.steps.registry import InputPath, StepContext, external, step

_SOURCES = ("database", "denovo")


def _path(p: InputPath) -> Path:
    if isinstance(p, list):
        raise ValueError(f"expected one file, got {len(p)}")
    return p


def _opt(ctx: StepContext, key: str) -> Path | None:
    p = ctx.inputs.get(key)
    return None if p is None else _path(p)


# ---------------------------------------------------------------- query selection (pure)


def resolve_sources(peptides: pd.DataFrame, cfg: GiapConfig) -> pd.Series:
    """Per-row source, never ``auto`` (Q3 default applied table-wide)."""
    default = cfg.input.peptide_source
    if default == "auto":
        scored = any(
            c in peptides.columns and pd.to_numeric(peptides[c], errors="coerce").notna().any()
            for c in ("q_value", "pep")
        )
        default = "database" if scored else "denovo"
    if "source" not in peptides.columns:
        return pd.Series([default] * len(peptides), index=peptides.index, dtype=object)
    src = peptides["source"].astype(object)
    return src.where(src.isin(_SOURCES), default).astype(object)


def select_queries(peptides: pd.DataFrame, cfg: GiapConfig) -> pd.DataFrame:
    """Query peptides (one row each, sorted by ``peptide``) with a resolved ``source``.

    Raises:
        ValueError: when the table has no ``peptide`` column.
    """
    if "peptide" not in peptides.columns:
        raise ValueError("peptide table has no 'peptide' column")
    lo, hi = cfg.mapping.exact.min_len, cfg.mapping.exact.max_len
    pep = peptides["peptide"].astype(str)
    keep = pep.str.len().between(lo, hi).to_numpy(dtype=bool).copy()
    if "is_decoy" in peptides.columns:
        keep &= ~peptides["is_decoy"].fillna(False).astype(bool).to_numpy()
    if "q_value" in peptides.columns:
        q = pd.to_numeric(peptides["q_value"], errors="coerce").to_numpy(dtype=float)
        keep &= ~(q > cfg.input.q_value_max)
    out = peptides.loc[keep].copy()
    out["peptide"] = pep[keep]
    out["source"] = resolve_sources(peptides, cfg)[keep]
    out = out.drop_duplicates("peptide").sort_values("peptide", kind="mergesort")
    return out.reset_index(drop=True)


def decoy_types(queries: pd.DataFrame, cfg: GiapConfig) -> dict[str, DecoyType]:
    """Peptide -> decoy type of its source (``mapping.decoys``)."""
    d = cfg.mapping.decoys
    kind: dict[str, DecoyType] = {"database": d.database_input, "denovo": d.denovo_input}
    return {
        str(p): kind[str(s)] for p, s in zip(queries["peptide"], queries["source"], strict=True)
    }


def length_label(n: int, classes: Sequence[Sequence[int]]) -> str | None:
    """Label of the first ``mapping.length_classes`` class holding ``n`` residues."""
    for b in classes:
        if int(b[0]) <= n <= int(b[1]):
            return dnp.class_label(b)
    return None


SAME_DB_MIN = 0.95  # M3-10 decision 4


def _id_list(value: object) -> list[str]:
    """Engine protein IDs of one row (contract C9 list/array, or a ``;``-joined string)."""
    if value is None or isinstance(value, float):
        return []
    if isinstance(value, str):
        return [x for x in value.split(";") if x]
    if isinstance(value, Iterable):
        return [str(x) for x in value if x is not None and str(x)]
    return [str(value)]


def resolve_search_db(peptides: pd.DataFrame, database_ids: Collection[str]) -> str:
    """``search_db: auto`` for the table (M3-10 decision 4): ``same`` when at least 95 % of the
    target peptides with engine protein IDs name an ID in ``database_ids``, else ``other``."""
    df = peptides
    if "is_decoy" in df.columns:
        df = df[~df["is_decoy"].fillna(False).astype(bool)]
    ids = [_id_list(v) for v in df["protein_ids"]] if "protein_ids" in df.columns else []
    ids = [i for i in ids if i]
    if not ids:
        return "other"
    hit = sum(any(x in database_ids for x in i) for i in ids)
    return "same" if hit >= SAME_DB_MIN * len(ids) - 1e-9 else "other"


# ---------------------------------------------------------------- work-file readers


def _queries(ctx: StepContext) -> pd.DataFrame:
    return select_queries(pd.read_parquet(_path(ctx.inputs["peptides"])), ctx.cfg)


def _resolved_queries(ctx: StepContext, gm: GeneModel) -> pd.DataFrame:
    """Queries with ``search_db: auto`` (or missing) resolved table-wide and logged."""
    table = pd.read_parquet(_path(ctx.inputs["peptides"]))
    queries = select_queries(table, ctx.cfg)
    sdb = queries["search_db"] if "search_db" in queries.columns else None
    auto = (
        np.ones(len(queries), dtype=bool)
        if sdb is None
        else ~sdb.astype(object).isin(["same", "other"]).to_numpy(dtype=bool)
    )
    if auto.any():
        value = resolve_search_db(table, {g.locus_tag for g in gm.genes})
        queries["search_db"] = value if sdb is None else sdb.astype(object)
        queries.loc[auto, "search_db"] = value
        ctx.log.info("%s: search_db auto resolved to %s for %d queries",
                     ctx.log.name.rsplit(".", 1)[-1], value, int(auto.sum()))  # fmt: skip
    return queries


def _genome(path: Path) -> dict[str, str]:
    return {r.id: r.sequence for r in read_genome_fasta(path)}


def _replicons(ctx: StepContext) -> dict[str, Replicon]:
    rows = pq.read_table(_path(ctx.inputs["replicons"])).to_pylist()
    return {r["contig"]: Replicon(r["contig"], r["length"], r["topology"], r["kind"]) for r in rows}


def _gene_model(ctx: StepContext) -> GeneModel:
    """GeneModel from the ``parse_bakta`` tables (0-based half-open)."""
    genes = tuple(
        Gene(
            r["locus_tag"], r["feature_type"],
            Interval(r["contig"], r["start"], r["end"], r["strand"]), r["aa"] or "", r["nt"],
            r["start_codon"], r["gene"], r["product"] or "", bool(r["partial"]),
            bool(r["hypothetical"]),
        )
        for r in pq.read_table(_path(ctx.inputs["genes"])).to_pylist()
    )  # fmt: skip
    feats = tuple(
        NcFeature(
            r["feature_id"], r["type"], Interval(r["contig"], r["start"], r["end"], r["strand"]),
            r["gene"], r["product"] or "", r["amino_acid"], r["anticodon"],
        )
        for r in pq.read_table(_path(ctx.inputs["features"])).to_pylist()
    )  # fmt: skip
    return GeneModel(replicons=_replicons(ctx), genes=genes, features=feats)


def _hits(ctx: StepContext, *keys: str) -> pd.DataFrame:
    parts = [pd.read_parquet(_path(ctx.inputs[k])) for k in keys]
    parts = [p for p in parts if len(p)]
    return pd.concat(parts, ignore_index=True) if parts else hits_frame()


# ---------------------------------------------------------------- peptide steps


@step("sixframe_index", inputs=["genome", "replicons", "peptides"], outputs=["sixframe"])
def sixframe_index(ctx: StepContext) -> None:
    """Six-frame translation plus the decoy frames the queries' sources need (§7.2.1)."""
    reps = _replicons(ctx)
    index = build_sixframe(_genome(_path(ctx.inputs["genome"])), reps)
    queries = _queries(ctx)
    kinds = sorted(set(decoy_types(queries, ctx.cfg).values()))
    if not kinds:
        kinds = [ctx.cfg.mapping.decoys.database_input]
    window = ctx.cfg.mapping.decoys.shuffle_window_aa
    decoys = {k: build_decoys(index, k, window, ctx.rng) for k in kinds}
    save_index(ctx.outputs["sixframe"], index, decoys)
    ctx.log.info(
        "sixframe_index: %d frames on %d contigs; decoys %s for %d queries",
        len(index.frames), len(reps), ",".join(kinds), len(queries),
    )  # fmt: skip


def _log_hits(ctx: StepContext, n_queries: int, hits: pd.DataFrame) -> None:
    by = hits.groupby("search")["peptide"].nunique() if len(hits) else pd.Series(dtype=int)
    ctx.log.info(
        "%s: %d queries; %d target rows (%d peptides), %d decoy rows (%d peptides)",
        ctx.log.name.rsplit(".", 1)[-1], n_queries,
        int((hits["search"] == "target").sum()), int(by.get("target", 0)),
        int((hits["search"] == "decoy").sum()), int(by.get("decoy", 0)),
    )  # fmt: skip


@step("map_exact", inputs=["peptides", "sixframe"], outputs=["hits_exact"])
def map_exact(ctx: StepContext) -> None:
    """Route 1: exact target and decoy hits of every query (§7.2.2)."""
    queries = _queries(ctx)
    index, decoys = load_index(_path(ctx.inputs["sixframe"]))
    ex = ctx.cfg.mapping.exact
    hits = match_exact(
        queries["peptide"].tolist(), index, decoys, decoy_types(queries, ctx.cfg),
        min_len=ex.min_len, max_len=ex.max_len,
    )  # fmt: skip
    hits.to_parquet(ctx.outputs["hits_exact"], index=False)
    _log_hits(ctx, len(queries), hits)


@step("map_tolerant", inputs=["peptides", "sixframe", "hits_exact"], outputs=["hits_tolerant"])
def map_tolerant(ctx: StepContext) -> None:
    """Route 2: tolerant target search for exact leftovers, decoys for every query (§7.2.4)."""
    queries = _queries(ctx)
    tcfg = ctx.cfg.mapping.tolerant
    peps = [p for p in queries["peptide"].tolist() if len(p) >= tcfg.min_len]
    exact = pd.read_parquet(_path(ctx.inputs["hits_exact"]))
    left = set(peps) - set(exact.loc[exact["search"] == "target", "peptide"])
    parts: list[pd.DataFrame] = []
    if tcfg.enabled and peps:
        index, decoys = load_index(_path(ctx.inputs["sixframe"]))
        decoy_of = decoy_types(queries, ctx.cfg)
        parts.append(match_tolerant(peps, index, decoys, decoy_of, left, tcfg))
        if tcfg.frameshift.enabled:
            parts.append(match_frameshift(peps, index, decoys, decoy_of, left, tcfg))
    parts = [p for p in parts if len(p)]
    hits = pd.concat(parts, ignore_index=True) if parts else hits_frame()
    hits.to_parquet(ctx.outputs["hits_tolerant"], index=False)
    _log_hits(ctx, len(peps), hits)


def _reference_exact(
    path: Path, cfg: GiapConfig, index: SixFrameIndex, decoys: Mapping[DecoyType, tuple[Frame, ...]]
) -> pd.DataFrame:
    """Reference peptides (generic table) with ``target_exact`` (exact target hit, I = L)."""
    raw, _ = read_peptide_table(path)
    ref, _ = validate_peptides(raw, cfg.input)
    ref = ref[["peptide", "q_value", "is_decoy"]].copy()
    ref["peptide"] = ref["peptide"].astype(str)
    ref["is_decoy"] = ref["is_decoy"].fillna(False).astype(bool)
    ref["q_value"] = pd.to_numeric(ref["q_value"], errors="coerce").astype(float)
    peps = sorted(set(ref["peptide"]))
    kind = next(iter(decoys))
    ex = cfg.mapping.exact
    hx = match_exact(
        peps, index, decoys, dict.fromkeys(peps, kind), min_len=ex.min_len, max_len=ex.max_len
    )
    found = set(hx.loc[hx["search"] == "target", "peptide"])
    ref["target_exact"] = ref["peptide"].isin(found)
    return ref


def _denovo(
    ctx: StepContext,
    dn: pd.DataFrame,
    exact: pd.DataFrame,
    index: SixFrameIndex,
    decoys: Mapping[DecoyType, tuple[Frame, ...]],
) -> tuple[np.ndarray, list[dict[str, object]]]:
    """π̂_wrong per de novo query (ADR-0010) and the ``denovo_precision`` rows (§7.2.5)."""
    classes = ctx.cfg.mapping.length_classes
    tx = set(exact.loc[exact["search"] == "target", "peptide"])
    dx = set(exact.loc[exact["search"] == "decoy", "peptide"])
    labels = [length_label(len(p), classes) for p in dn["peptide"]]
    score = (
        pd.to_numeric(dn["score"], errors="coerce").to_numpy(dtype=float)
        if "score" in dn.columns
        else np.full(len(dn), np.nan)
    )
    df = pd.DataFrame(
        {
            "peptide": dn["peptide"].to_numpy(),
            "score": score,
            "length_class": labels,
            "target_exact": dn["peptide"].isin(tx).to_numpy(),
            "decoy_exact": dn["peptide"].isin(dx).to_numpy(),
        }
    )
    has_cls = np.array([lab is not None for lab in labels], dtype=bool)
    df = df[has_cls].reset_index(drop=True)
    ref_path = _opt(ctx, "reference_peptides")
    ref = None if ref_path is None else _reference_exact(ref_path, ctx.cfg, index, decoys)
    c, assumed = dnp.c_by_class(ref, classes)
    if assumed:
        ctx.log.warning(
            "placement: no usable reference peptides; c_l = %.2f assumed (run flag c_assumed)",
            dnp.C_DEFAULT,
        )
    bins = dnp.denovo_precision(df, c, assumed, ctx.rng)
    pi = np.ones(len(dn))
    pi[has_cls] = dnp.pi_wrong_denovo(bins, dnp.score_bins(df["score"].to_numpy(dtype=float)))
    dtype = ctx.cfg.mapping.decoys.denovo_input
    rows: list[dict[str, object]] = [
        {
            "section": "denovo_precision", "decoy_type": dtype, "score_bin": b.score_bin,
            "precision": b.precision, "ci_low": b.ci_low, "ci_high": b.ci_high,
            "c_value": b.c_value, "flags": ";".join(b.flags),
        }
        for b in bins
    ]  # fmt: skip
    return pi, rows


@step(
    "placement",
    inputs=[
        "peptides", "sixframe", "hits_exact", "hits_tolerant", "genes", "features", "replicons",
    ],
    optional_inputs=["reference_peptides"],
    outputs=["placement", "mapping_fdr"],
)  # fmt: skip
def placement_step(ctx: StepContext) -> None:
    """Best costs, de novo precision, placement q-values and ``mapping_fdr`` rows (§7.2.5)."""
    queries = _resolved_queries(ctx, _gene_model(ctx))
    index, decoys = load_index(_path(ctx.inputs["sixframe"]))
    exact = pd.read_parquet(_path(ctx.inputs["hits_exact"]))
    queries["pi_denovo"] = np.nan
    is_dn = (queries["source"] == "denovo").to_numpy()
    dn_rows: list[dict[str, object]] = []
    if is_dn.any():
        pi, dn_rows = _denovo(ctx, queries[is_dn], exact, index, decoys)
        queries.loc[is_dn, "pi_denovo"] = pi
    costs = best_costs(_hits(ctx, "hits_exact", "hits_tolerant"), queries, residue_freqs(index))
    per, fdr = placement(costs, ctx.cfg)
    if dn_rows:
        fdr = pd.concat([fdr, pd.DataFrame(dn_rows)], ignore_index=True)
    per.to_parquet(ctx.outputs["placement"], index=False)
    fdr.to_parquet(ctx.outputs["mapping_fdr"], index=False)
    ctx.log.info(
        "placement: %d queries (%d de novo), %d with a target hit, %d accepted",
        len(queries), int(is_dn.sum()), len(per), int(per["accepted"].sum()),
    )  # fmt: skip


# ---------------------------------------------------------------- protein route


external("miniprot", inputs=["genome", "proteins_long"], outputs=["miniprot_gff"])


def _miniprot(ctx: StepContext, gm: GeneModel, lens: Mapping[str, int]) -> list[ProteinAlignment]:
    gff = _opt(ctx, "miniprot_gff")
    return [] if gff is None else parse_miniprot_gff(gff, gm.replicons, lens)


@step(
    "protein_leftovers",
    inputs=["proteins_long", "miniprot_gff", "genes", "features", "replicons"],
    outputs=["proteins_leftover"],
)
def protein_leftovers(ctx: StepContext) -> None:
    """Proteins without an assigned miniprot alignment, the tblastn queries (§7.2, §7.2.3)."""
    recs = read_protein_fasta(_path(ctx.inputs["proteins_long"]))
    gm = _gene_model(ctx)
    alns = _miniprot(ctx, gm, {r.id: len(r.sequence) for r in recs})
    rows, _ = assign_alignments(alns, gm, ctx.cfg.mapping.proteins)
    assigned = set(rows.loc[rows["gene"].notna(), "protein_id"])
    left = [(r.id, r.sequence) for r in recs if r.id not in assigned]
    write_fasta(ctx.outputs["proteins_leftover"], left)
    ctx.log.info("protein_leftovers: %d proteins, %d assigned by miniprot, %d to tblastn",
                 len(recs), len(assigned), len(left))  # fmt: skip


external("tblastn", inputs=["genome", "proteins_leftover"], outputs=["tblastn_tsv"])


def _proteins(ctx: StepContext, gm: GeneModel) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(protein_hits rows, protein state for the unmapped bucket)."""
    state_cols = ["protein_id", "sequence", "assigned", "has_alignment"]
    path = _opt(ctx, "proteins_long")
    if path is None:
        return pd.DataFrame(columns=list(PROTEIN_HITS_COLUMNS)), pd.DataFrame(columns=state_cols)
    recs = read_protein_fasta(path)
    alns = _miniprot(ctx, gm, {r.id: len(r.sequence) for r in recs})
    rows, _ = assign_alignments(alns, gm, ctx.cfg.mapping.proteins)
    assigned = set(rows.loc[rows["gene"].notna(), "protein_id"])
    aligned = {a.protein_id for a in alns}
    tsv = _opt(ctx, "tblastn_tsv")
    if tsv is not None:
        hsps = [h for h in parse_tblastn(tsv) if h.qseqid not in assigned]
        chains = merge_hsps(hsps, gm.replicons)
        more, _ = assign_alignments(chains, gm, ctx.cfg.mapping.proteins)
        assigned |= set(more.loc[more["gene"].notna(), "protein_id"])
        aligned |= {a.protein_id for a in chains}
        rows = pd.concat([p for p in (rows, more) if len(p)] or [rows], ignore_index=True)
    state = pd.DataFrame(
        [(r.id, r.sequence, r.id in assigned, r.id in aligned) for r in recs], columns=state_cols
    )
    return annotate_protein_hits(rows, gm), state


# ---------------------------------------------------------------- loci and outputs


def _peptide_rows(
    hits: pd.DataFrame, per: pd.DataFrame, queries: pd.DataFrame, cfg: GiapConfig
) -> pd.DataFrame:
    """Target hits at the best cost of every accepted peptide, with input and placement columns."""
    acc = per[per["accepted"].astype(bool)].set_index("peptide")
    tgt = hits[(hits["search"] == "target") & hits["peptide"].isin(acc.index)]
    best = tgt["peptide"].map(acc["m_target"])
    rows = tgt[tgt["n_diff"].to_numpy() == best.to_numpy()].copy().reset_index(drop=True)
    q = queries.set_index("peptide")
    p = rows["peptide"]
    rows["peptide_len"] = p.str.len()
    rows["source"] = p.map(acc["source"])
    rows["placement_q"] = p.map(acc["placement_q"])
    rows["flags"] = p.map(acc["flags"])
    for col, src in (("best_q_input", "q_value"), ("pep_input", "pep"), ("score_input", "score")):
        rows[col] = p.map(q[src]) if src in q.columns else np.nan
    rows["length_class"] = [length_label(len(s), cfg.mapping.length_classes) for s in p]
    return rows


@step(
    "loci",
    inputs=[
        "placement", "mapping_fdr", "hits_exact", "hits_tolerant", "peptides", "genome",
        "genes", "features", "replicons",
    ],
    optional_inputs=["proteins_long", "miniprot_gff", "tblastn_tsv", "contaminants"],
    outputs=["mapping"],
)  # fmt: skip
def loci(ctx: StepContext) -> None:
    """Locus classes and genes, protein assignment, unmapped bucket and the seven §7.2.9 files."""
    cfg = ctx.cfg
    gm = _gene_model(ctx)
    queries = _resolved_queries(ctx, gm)
    per = pd.read_parquet(_path(ctx.inputs["placement"]))
    gpath = _path(ctx.inputs["genome"])
    genome = _genome(gpath)
    rows = _peptide_rows(_hits(ctx, "hits_exact", "hits_tolerant"), per, queries, cfg)
    pep_hits = annotate_peptide_hits(
        rows, gm, genome, queries, cfg.genome_relation,
        decoy_prefixes=cfg.input.decoy_prefixes,
        contaminant_prefixes=cfg.input.contaminant_prefixes,
    )  # fmt: skip
    prot_hits, state = _proteins(ctx, gm)
    lib = _opt(ctx, "contaminants")
    placed = set(per.loc[per["accepted"].astype(bool), "peptide"])
    left = [p for p in queries["peptide"] if p not in placed]
    pep_matches = match_contaminant_peptides(left, read_protein_fasta(lib)) if lib else {}
    unmapped = unmapped_table(
        queries, per, state, pep_matches, {}, cfg.input.contaminant_prefixes,
        tolerant_min_len=cfg.mapping.tolerant.min_len, q_value_max=cfg.input.q_value_max,
    )  # fmt: skip
    fdr = pd.read_parquet(_path(ctx.inputs["mapping_fdr"]))
    out = ctx.outputs["mapping"]
    out.mkdir(parents=True, exist_ok=True)
    sha = hashlib.sha256(gpath.read_bytes()).hexdigest()
    write_mapping_outputs(out, pep_hits, prot_hits, unmapped, fdr, gm.replicons, sha)
    n_pep = pep_hits["peptide"].nunique() if len(pep_hits) else 0
    ctx.log.info(
        "loci: %d peptide rows (%d peptides), %d protein rows, %d unmapped",
        len(pep_hits), n_pep, len(prot_hits), len(unmapped),
    )  # fmt: skip
