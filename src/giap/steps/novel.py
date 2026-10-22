"""Novel-ORF steps and the second round (GIAP_SPEC.md §6.2, §7.6.2; M7-08, ADR-0024):
``novel_candidates``, ``promotech_novel`` (ext), ``novel_promoter_layer``, ``ostir_novel`` (ext),
``novel_rbs_layer`` and ``novel_orfs``.

``novel_candidates`` runs the M7-06 library on the loci tables and writes the ORFs
(``novel_candidates.parquet``, one row per ORF x peptide with the peptide's mapped interval),
the second-round Promotech windows
(``novel_windows/windows_NNN.fa`` + ``novel_windows_map.parquet``) and ``ostir_input_novel.csv``.
The two tools are never imported (GPL-3.0). ``novel_rbs_layer`` keys OSTIR features by
``(contig, strand, start)``, ``start`` the 0-based first base of the start codon in
transcription direction. ``novel_orfs`` joins the second-round scores by position and writes the
four M7-07 files after the M7-09 evidence (gene-table RBS threshold, terminator layer).
Under the local executor the precomputed ``promoters`` and ``rbs_novel`` TSVs
replace the tool and layer steps and arrive here instead of the Parquet files.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]

from giap.annotate.layers import terminators_of
from giap.annotate.ostir import parse_ostir_starts, read_ostir_input, write_ostir_input
from giap.annotate.precomputed import load_precomputed
from giap.annotate.promotech_cmd import parse_promotech
from giap.annotate.promoters import gene_windows, promoter_layer, write_chunks
from giap.io.fasta import read_genome_fasta
from giap.io.novel_tables import (
    write_annotation_patch,
    write_novel_faa,
    write_novel_gff3,
    write_novel_orfs,
    write_pragma_gff3,
)
from giap.io.parquet import write_parquet_deterministic
from giap.io.protein_fasta import read_protein_fasta
from giap.io.work_tables import FILE_NAMES, _read_tsv_table
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene, GeneModel
from giap.novel.candidates import (
    Orf,
    candidate_rows,
    minimal_call,
    ostir_rows,
    peptide_loci,
    promoter_windows,
    reconstruct,
    window_chunks,
)
from giap.novel.evidence import apply_evidence
from giap.novel.outputs import orf_table
from giap.steps.common import one_path as _path
from giap.steps.common import read_genes as _genes
from giap.steps.common import read_replicons as _replicons
from giap.steps.mapping import _gene_model
from giap.steps.promoters import _LAYER_COLS, _prediction_files, _rows
from giap.steps.rbs import _RBS_COLS
from giap.steps.registry import StepContext, external, step

_STR, _INT, _F64, _BOOL = pa.string(), pa.int64(), pa.float64(), pa.bool_()
_CAND_COLS = [
    ("orf_id", _STR), ("peptide", _STR), ("class", _STR), ("contig", _STR), ("strand", _STR),
    ("frame", _INT), ("stop", _INT), ("starts", pa.list_(_INT)), ("start_firsts", pa.list_(_INT)),
    ("start_longest", _INT), ("end_5p_open", _BOOL), ("end_3p_open", _BOOL), ("best_q", _F64),
    ("overlapping_gene", _STR), ("first_peptide", _INT),
    ("peptide_start", _INT), ("peptide_end", _INT),  # M7-09 d5: the peptide's mapped interval
]  # fmt: skip
_MAP_COLS = [
    ("header", _STR), ("contig", _STR), ("strand", _STR), ("end3p", _INT), ("owner", _STR),
    ("k", _INT), ("orf_id", _STR),
]  # fmt: skip
_RBS_NOVEL_COLS = [("contig", _STR), ("strand", _STR), ("start", _INT), *_RBS_COLS[1:]]
_PARQUET_MAGIC = b"PAR1"


def first_base(strand: str, lowest: int, n: int) -> int:
    """First base, in transcription direction, of the codon whose lowest coordinate is
    ``lowest`` (on ``-`` its highest coordinate, wrapped on a replicon of length ``n``)."""
    return lowest if strand == "+" else (lowest + 2) % n


# ---------------------------------------------------------------- novel_candidates


def _candidate_records(
    orfs: Sequence[Orf],
    reps: Mapping[str, Replicon],
    loci: Mapping[tuple[str, str], Interval],
) -> list[tuple[Any, ...]]:
    out = []
    for o in orfs:
        firsts = [first_base(o.strand, x, reps[o.contig].length) for x in o.starts]
        for p in o.peptides:
            iv = loci.get((o.contig, p))
            out.append((o.orf_id, p, o.cls, o.contig, o.strand, o.frame, o.stop, list(o.starts),
                        firsts, o.start_longest, o.end_5p_open, o.end_3p_open, o.best_q,
                        o.overlapping_gene, o.first_peptide, None if iv is None else iv.start,
                        None if iv is None else iv.end))  # fmt: skip
    return out


def _opt_int(v: Any) -> int | None:
    return None if pd.isna(v) else int(v)


def read_candidates(path: Path) -> tuple[list[Orf], dict[str, list[int]]]:
    """The ORFs of ``novel_candidates.parquet`` (file order) and the first base, in transcription
    direction, of each ORF's starts (5' to 3')."""
    df = pd.read_parquet(path)
    orfs: list[Orf] = []
    firsts: dict[str, list[int]] = {}
    for oid, g in df.groupby("orf_id", sort=False):
        r = g.iloc[0]
        orfs.append(Orf(
            orf_id=str(oid), cls=r["class"], contig=str(r["contig"]), strand=r["strand"],
            frame=int(r["frame"]), stop=int(r["stop"]), starts=tuple(int(x) for x in r["starts"]),
            start_longest=_opt_int(r["start_longest"]), end_5p_open=bool(r["end_5p_open"]),
            end_3p_open=bool(r["end_3p_open"]), peptides=tuple(sorted(map(str, g["peptide"]))),
            best_q=None if pd.isna(r["best_q"]) else float(r["best_q"]),
            overlapping_gene=None if pd.isna(r["overlapping_gene"]) else str(r["overlapping_gene"]),
            first_peptide=_opt_int(r["first_peptide"]),
        ))  # fmt: skip
        firsts[str(oid)] = [int(x) for x in r["start_firsts"]]
    return orfs, firsts


def read_peptide_intervals(path: Path) -> dict[str, list[Interval]]:
    """Per ORF, the mapped intervals of its peptides (``novel_candidates.parquet`` columns
    ``peptide_start``/``peptide_end``; empty for a file without them)."""
    df = pd.read_parquet(path)
    if not {"peptide_start", "peptide_end"} <= set(df.columns):
        return {}
    out: dict[str, list[Interval]] = {}
    cols = (df[c] for c in ("orf_id", "contig", "strand", "peptide_start", "peptide_end"))
    for oid, c, s, a, b in zip(*cols, strict=True):
        if pd.notna(a) and pd.notna(b):
            out.setdefault(str(oid), []).append(Interval(str(c), int(a), int(b), s))
    return out


@step(
    "novel_candidates",
    inputs=["mapping", "peptides", "genes", "features", "replicons", "genome", "annotation_meta"],
    optional_inputs=["contaminants"],
    outputs=["novel_candidates", "novel_windows", "novel_windows_map", "ostir_input_novel"],
)
def novel_candidates(ctx: StepContext) -> None:
    """Loci tables → candidate ORFs, second-round Promotech windows and OSTIR rows (M7-06)."""
    ins, cfg = ctx.inputs, ctx.cfg
    hits_path = _path(ins["mapping"]) / FILE_NAMES["peptide_hits"]
    hits = _read_tsv_table("peptide_hits", hits_path)
    peptides = pd.read_parquet(_path(ins["peptides"]))
    lib = ins.get("contaminants")
    contaminants = None if lib is None else read_protein_fasta(_path(lib))
    gm = _gene_model(ctx)
    genome = {r.id: r.sequence for r in read_genome_fasta(_path(ins["genome"]))}
    rows = candidate_rows(hits, peptides, contaminants)
    orfs = reconstruct(rows, gm, genome, cfg)
    write_parquet_deterministic(
        ctx.outputs["novel_candidates"],
        _CAND_COLS,
        _candidate_records(orfs, gm.replicons, peptide_loci(rows)),
    )
    chunks, wmap = window_chunks(promoter_windows(orfs, gm, genome, cfg))
    write_chunks(ctx.outputs["novel_windows"], chunks)
    write_parquet_deterministic(ctx.outputs["novel_windows_map"], _MAP_COLS, _rows(wmap))
    anti_sd = str(json.loads(_path(ins["annotation_meta"]).read_text())["anti_sd"])
    o_rows = ostir_rows(orfs, gm, genome, anti_sd, cfg)
    write_ostir_input(ctx.outputs["ostir_input_novel"], o_rows)
    ctx.log.info("novel_candidates: %d ORFs, %d windows, %d OSTIR rows",
                 len(orfs), sum(map(len, chunks)), len(o_rows))  # fmt: skip


# ---------------------------------------------------------------- second-round layers

external("promotech_novel", inputs=["novel_windows"], outputs=["promotech_novel"],
         precomputed="promoters")  # fmt: skip


@step(
    "novel_promoter_layer",
    inputs=["novel_windows_map", "promotech_novel"],
    outputs=["novel_promoter_windows"],
    precomputed="promoters",
)
def novel_promoter_layer(ctx: StepContext) -> None:
    """M2-13 ``promoter_layer`` on the second-round window map."""
    wmap = pd.read_parquet(_path(ctx.inputs["novel_windows_map"]))
    scores: dict[str, float] = {}
    for path in _prediction_files(ctx.inputs["promotech_novel"]):
        scores.update(parse_promotech(path))
    layer = promoter_layer(wmap, scores)
    write_parquet_deterministic(ctx.outputs["novel_promoter_windows"], _LAYER_COLS, _rows(layer))


external("ostir_novel", inputs=["ostir_input_novel"], outputs=["ostir_novel"],
         precomputed="rbs_novel")  # fmt: skip


def _window_starts(seq: str, start: int, end: int, codons: frozenset[str]) -> list[int]:
    """1-based in-frame start-codon positions of an OSTIR window in ``[start, end]``."""
    return [p for p in range(start, end + 1, 3) if seq[p - 1 : p + 2].upper() in codons]


@step(
    "novel_rbs_layer",
    inputs=["ostir_novel", "ostir_input_novel", "novel_candidates"],
    outputs=["rbs_novel"],
    precomputed="rbs_novel",
)
def novel_rbs_layer(ctx: StepContext) -> None:
    """``ostir_novel.csv`` → ``rbs_novel.parquet``: one row per start OSTIR scored."""
    codons = frozenset(c.upper() for c in ctx.cfg.annotation.start_codons)
    orfs, firsts = read_candidates(_path(ctx.inputs["novel_candidates"]))
    by_id = {o.orf_id: o for o in orfs}
    rows = read_ostir_input(_path(ctx.inputs["ostir_input_novel"]))
    req: dict[str, dict[int, int]] = {}  # id -> window position -> genome first base
    for r in rows:
        ps = _window_starts(r.seq, r.start, r.end, codons)
        if len(ps) != len(firsts[r.id]):
            raise ValueError(f"{r.id}: {len(ps)} window starts, ORF has {len(firsts[r.id])}")
        req[r.id] = dict(zip(ps, firsts[r.id], strict=True))
    feats = parse_ostir_starts(_path(ctx.inputs["ostir_novel"]),
                               {k: set(v) for k, v in req.items()})  # fmt: skip
    out = []
    for r in rows:
        o, got = by_id[r.id], feats.get(r.id, {})
        for p, pos in req[r.id].items():
            if (f := got.get(p)) is not None:
                out.append((o.contig, o.strand, pos, f.tir_ostir, f.sd_dg, f.mrna_dg, f.spacing_dg,
                            f.standby_dg, f.spacer_bp, r.seq[p - 1 : p + 2].upper()))  # fmt: skip
    write_parquet_deterministic(ctx.outputs["rbs_novel"], _RBS_NOVEL_COLS, out)


# ---------------------------------------------------------------- novel_orfs


def _is_parquet(path: Path) -> bool:
    with path.open("rb") as fh:
        return fh.read(4) == _PARQUET_MAGIC


def _layer(
    ctx: StepContext, key: str, pre: str, reps: Mapping[str, Replicon]
) -> pd.DataFrame | None:
    """The Parquet a layer step wrote, or the precomputed M2-02 TSV (``load_precomputed``)."""
    given = ctx.inputs.get(key)
    if given is None:
        return None
    path = _path(given)
    if _is_parquet(path):
        return pd.read_parquet(path)
    df = load_precomputed(pre, path, reps)
    if pre == "promoters":  # back to the M2-13 1-based end3p of the window map and Parquet
        df = df.assign(end3p=df["end3p"] + 1)
    return df


def _promoter_scores(frames: Sequence[pd.DataFrame | None]) -> dict[tuple[str, str, int], float]:
    out: dict[tuple[str, str, int], float] = {}
    for df in frames:
        if df is None:
            continue
        for c, s, e, v in zip(df["contig"], df["strand"], df["end3p"], df["score"], strict=True):
            key = (str(c), str(s), int(e))
            if pd.notna(v) and (key not in out or float(v) > out[key]):
                out[key] = float(v)
    return out


def second_round_scores(
    orfs: Sequence[Orf],
    rbs: Mapping[tuple[str, str, int], float] | None,
    promoters: Mapping[tuple[str, str, int], float] | None,
    gm: GeneModel,
    genome: Mapping[str, str],
    k_range: tuple[int, int],
    width: int,
) -> pd.DataFrame:
    """Per ORF: ``start_best_ostir`` (lowest codon coordinate; highest ``tir_ostir``, ties →
    the more upstream start), ``ostir_tir`` and the ``promoter_score`` at the used start (best,
    else longest): the maximum over its windows' ``(contig, strand, end3p)``."""
    recs = []
    for o in orfs:
        n = gm.replicons[o.contig].length
        best, tir = None, float("nan")
        for x in o.starts if rbs is not None else ():
            v = rbs.get((o.contig, o.strand, first_base(o.strand, x, n)))  # type: ignore[union-attr]
            if v is not None and (best is None or v > tir):
                best, tir = x, v
        used = best if best is not None else o.start_longest
        score = float("nan")
        if used is not None and promoters is not None:
            anchor = Gene(
                o.orf_id, "cds", Interval(o.contig, used, used + 3, o.strand), "", None, None
            )
            rep = gm.replicons[o.contig]
            wins = gene_windows(anchor, rep, genome[o.contig], k_range, width)
            vals = [promoters[k] for w in wins if (k := (w.contig, w.strand, w.end3p)) in promoters]
            score = max(vals) if vals else score
        recs.append((o.orf_id, best, tir, score))
    df = pd.DataFrame(recs, columns=["orf_id", "start_best_ostir", "ostir_tir", "promoter_score"])
    return df.astype({"start_best_ostir": "object"})


@step(
    "novel_orfs",
    inputs=["novel_candidates", "gene_table", "genes", "replicons", "genome"],
    optional_inputs=["novel_promoter_windows", "rbs_novel", "promoter_windows", "terminators"],
    outputs=["novel_orfs", "novel_orfs_gff3", "novel_orfs_faa", "annotation_patch"],
)
def novel_orfs(ctx: StepContext) -> None:
    """Candidates + second-round scores + evidence (M7-09) → ``novel_orfs.tsv``/``.gff3``/
    ``.faa`` and ``annotation_patch.gff3`` (class FDR: M7-11)."""
    ins, cfg = ctx.inputs, ctx.cfg
    reps = _replicons(_path(ins["replicons"]))
    gm = GeneModel(replicons=reps, genes=tuple(_genes(_path(ins["genes"]))), features=())
    genome = {r.id: r.sequence for r in read_genome_fasta(_path(ins["genome"]))}
    orfs, _ = read_candidates(_path(ins["novel_candidates"]))
    missing = [k for k in ("rbs_novel", "terminators") if k not in ins]
    if "novel_promoter_windows" not in ins and "promoter_windows" not in ins:
        missing.append("promoters")
    for key in sorted(missing):
        ctx.log.warning(
            "layer_missing:%s: no %s layer; its novel_orfs columns are NA", key, key,
            extra={"event": "novel_orfs.layer_missing",
                   "data": {"layer": key, "flag": f"layer_missing:{key}"}},
        )  # fmt: skip
    rbs_df = _layer(ctx, "rbs_novel", "rbs_novel", reps)
    rbs = None
    if rbs_df is not None:
        cols = (rbs_df[c] for c in ("contig", "strand", "start", "tir_ostir"))
        rbs = {(str(c), str(s), int(p)): float(v) for c, s, p, v in zip(*cols, strict=True)
               if pd.notna(v)}  # fmt: skip
    frames = [_layer(ctx, k, "promoters", reps) for k in ("novel_promoter_windows",
                                                          "promoter_windows")]  # fmt: skip
    prom = None if all(f is None for f in frames) else _promoter_scores(frames)
    pcfg = cfg.annotation.promoters
    k_range = (pcfg.window_3p_range[0], pcfg.window_3p_range[1])
    sec = second_round_scores(orfs, rbs, prom, gm, genome, k_range, pcfg.window_nt)
    scored = minimal_call(orfs, cfg).merge(sec, on="orf_id", how="left")
    table = orf_table(orfs, scored, genome=genome, replicons=reps)
    terms_df = _layer(ctx, "terminators", "terminators", reps)
    table = apply_evidence(
        table, orfs, gene_table=pd.read_parquet(_path(ins["gene_table"])),
        terminators=None if terms_df is None else terminators_of(terms_df), gm=gm,
        genome=genome, cfg=cfg,
        peptide_intervals=read_peptide_intervals(_path(ins["novel_candidates"])),
    )  # fmt: skip
    write_novel_orfs(ctx.outputs["novel_orfs"], table)
    write_novel_faa(ctx.outputs["novel_orfs_faa"], table, genome, reps)
    if len(table):
        write_novel_gff3(ctx.outputs["novel_orfs_gff3"], table, reps)
        write_annotation_patch(ctx.outputs["annotation_patch"], table, gm)
    else:  # rule 3: no candidate -> pragma-only GFF3s
        write_pragma_gff3(ctx.outputs["novel_orfs_gff3"], reps)
        write_pragma_gff3(ctx.outputs["annotation_patch"], reps)
    ctx.log.info("novel_orfs: %d ORFs, %d called", len(table), int(table["called"].sum()))
