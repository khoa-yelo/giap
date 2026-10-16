"""Oracles for the M2-20 ``context_layers`` tests (test author).

Two kinds:

* **Independent** (no ``giap.annotate.layers``): TU position / size, ``uniop_p_prev`` /
  ``uniop_p_next`` in transcription order straight from ``pairs_raw`` (spec column definitions),
  ``region_ids`` from ``region_genes``, the gene ``flags`` set.
* **Merged producers** (the functions M2-12/13/14/15/17/19 already ship): the per-gene promoter
  / tx / terminator features and the ``annotation/*`` writers. M2-20 only wires them, so its
  files must equal what these writers produce from the same producer outputs.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pandas as pd

from giap.annotate.operons import Tu
from giap.annotate.terminators import Terminator
from giap.annotate.uniop_adapter import PairProb
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene

# ----------------------------------------------------------------------------- frame -> objects


def replicons(df: pd.DataFrame) -> dict[str, Replicon]:
    return {
        str(r.contig): Replicon(str(r.contig), int(r.length), r.topology, r.kind)
        for r in df.itertuples(index=False)
    }


def genes(df: pd.DataFrame) -> list[Gene]:
    return [
        Gene(
            locus_tag=str(r.locus_tag),
            feature_type=r.feature_type,
            interval=Interval(str(r.contig), int(r.start), int(r.end), r.strand),
            aa=r.aa or "",
            nt=r.nt,
            start_codon=r.start_codon,
            gene=r.gene,
            product=r.product or "",
            partial=bool(r.partial),
        )
        for r in df.itertuples(index=False)
    ]


def tus(df: pd.DataFrame) -> list[Tu]:
    return [
        Tu(
            str(r.tu_id), str(r.contig),
            Interval(str(r.contig), int(r.start), int(r.end), r.strand), r.strand, tuple(r.genes),
            Interval(str(r.contig), int(r.upstream_start), int(r.upstream_end), r.strand),
        )
        for r in df.itertuples(index=False)
    ]  # fmt: skip


def pairs(df: pd.DataFrame) -> list[PairProb]:
    return [
        PairProb(str(a), str(b), float(p), s)
        for a, b, p, s in df[["gene_a", "gene_b", "uniop_p", "source"]].itertuples(index=False)
    ]


def terminators(df: pd.DataFrame) -> list[Terminator]:
    out = []
    for r in df.itertuples(index=False):
        dg = getattr(r, "hairpin_dg", None)
        out.append(
            Terminator(
                Interval(str(r.contig), int(r.start), int(r.end), r.strand),
                float(r.confidence), r.source,
                None if dg is None or (isinstance(dg, float) and math.isnan(dg)) else float(dg),
            )
        )  # fmt: skip
    return out


def regions(df: pd.DataFrame) -> list[Any]:
    from giap.annotate.regions import Region

    out = []
    for r in df.itertuples(index=False):
        has = not (pd.isna(r.start) or pd.isna(r.end))
        iv = Interval(str(r.contig), int(r.start), int(r.end), ".") if has else None
        score = None if r.score is None or pd.isna(r.score) else float(r.score)
        out.append(
            Region(str(r.region_id), r.region_type, str(r.contig), iv, r.source_tool, score,
                   tuple(r.genes if r.genes is not None else ()), r.notes or "")
        )  # fmt: skip
    return out


# ----------------------------------------------------------------------------- independent


def tu_positions(tus_df: pd.DataFrame) -> dict[str, tuple[str, int, int]]:
    """locus_tag -> (tu_id, 1-based position in transcription order, TU size)."""
    out = {}
    for tu_id, gs in zip(tus_df["tu_id"], tus_df["genes"], strict=True):
        gs = list(gs)
        for i, g in enumerate(gs, start=1):
            out[str(g)] = (str(tu_id), i, len(gs))
    return out


def prev_next(genes_df: pd.DataFrame, pairs_df: pd.DataFrame) -> dict[str, tuple[float, float]]:
    """locus_tag -> (uniop_p_prev, uniop_p_next) in transcription order (NaN = NA).

    Pairs are (gene_a, gene_b) in genome order: on ``+`` the previous gene is ``gene_a`` of the
    pair whose ``gene_b`` is the gene; on ``-`` it is ``gene_b`` of the pair whose ``gene_a`` is
    the gene."""
    as_a = {str(a): float(p) for a, p in zip(pairs_df["gene_a"], pairs_df["uniop_p"], strict=True)}
    as_b = {str(b): float(p) for b, p in zip(pairs_df["gene_b"], pairs_df["uniop_p"], strict=True)}
    nan = float("nan")
    out = {}
    for tag, strand in zip(genes_df["locus_tag"], genes_df["strand"], strict=True):
        tag = str(tag)
        if strand == "+":
            out[tag] = (as_b.get(tag, nan), as_a.get(tag, nan))
        else:
            out[tag] = (as_a.get(tag, nan), as_b.get(tag, nan))
    return out


def region_ids(region_genes_df: pd.DataFrame | None) -> dict[str, list[str]]:
    out: dict[str, set[str]] = {}
    if region_genes_df is not None:
        for rid, tag in zip(
            region_genes_df["region_id"], region_genes_df["locus_tag"], strict=True
        ):
            out.setdefault(str(tag), set()).add(str(rid))
    return {k: sorted(v) for k, v in out.items()}


def gene_flags(genes_df: pd.DataFrame, codon_df: pd.DataFrame | None) -> dict[str, set[str]]:
    cai = set()
    if codon_df is not None:
        cai = {
            str(t)
            for t, c in zip(codon_df["locus_tag"], codon_df["cai_reference"], strict=True)
            if c
        }
    out = {}
    for r in genes_df.itertuples(index=False):
        f = set()
        if r.partial:
            f.add("partial")
        if r.crosses_origin:
            f.add("crosses_origin")
        if r.feature_type == "pseudogene":
            f.add("pseudogene")
        if str(r.locus_tag) in cai:
            f.add("cai_reference")
        out[str(r.locus_tag)] = f
    return out


def gc_local_dev(
    genes_df: pd.DataFrame,
    replicons_df: pd.DataFrame,
    seqs: Mapping[str, str],
    window: int = 10_000,
) -> list[float]:
    """Brute-force ``gc_local_dev`` (GIAP_SPEC §7.8.1 / §9.3: 10-kb window GC - genome GC), one
    value per ``genes_df`` row (NaN = NA). Pinned by the test author (lead decision 4):

    * GC fraction = (G + C) / (A + C + G + T), case-insensitive; every other letter (N, IUPAC
      ambiguity codes) is excluded from numerator and denominator. Values are fractions, not %.
    * genome GC over the whole genome (every sequence of the FASTA), not per replicon.
    * window ``[m - window // 2, m - window // 2 + window)`` around the gene midpoint
      ``m = (start + end) // 2`` (0-based half-open; origin-crossing genes are stored unwrapped,
      so ``m`` is taken modulo the length on circular replicons). Circular replicons wrap; a
      circular replicon not longer than the window counts each base once (the whole replicon);
      linear replicons truncate the window at their ends.
    * NA when the window has no A/C/G/T base. Strand plays no role.
    """
    acgt, gc = set("ACGT"), set("GC")

    def frac(bases: list[str]) -> float:
        d = sum(b in acgt for b in bases)
        return float("nan") if d == 0 else sum(b in gc for b in bases) / d

    genome_gc = frac([b for s in seqs.values() for b in s.upper()])
    reps = {str(r.contig): (int(r.length), r.topology) for r in replicons_df.itertuples()}
    half = window // 2
    out = []
    for r in genes_df.itertuples(index=False):
        seq = seqs[str(r.contig)].upper()
        length, topo = reps[str(r.contig)]
        m = (int(r.start) + int(r.end)) // 2
        if topo == "circular":
            m %= length
            pos = (
                list(range(length))
                if length <= window
                else [p % length for p in range(m - half, m - half + window)]
            )
        else:
            pos = list(range(max(0, m - half), min(length, m - half + window)))
        out.append(frac([seq[p] for p in pos]) - genome_gc)
    return out


# ----------------------------------------------------------------------------- merged producers


def promoter_features(frames: Mapping[str, pd.DataFrame | None], yes: float = 0.5) -> pd.DataFrame:
    from giap.annotate.promoters import gene_promoter_features

    reps = replicons(frames["replicons"])  # type: ignore[arg-type]
    return gene_promoter_features(
        genes(frames["genes"]), tus(frames["tus"]), frames["promoter_windows"], yes,  # type: ignore[arg-type]
        replicons=reps,
    ).set_index("locus_tag")  # fmt: skip


def tx_features(frames: Mapping[str, pd.DataFrame | None]) -> pd.DataFrame:
    from giap.annotate.promoter_calculator import tu_tx_features

    reps = replicons(frames["replicons"])  # type: ignore[arg-type]
    return tu_tx_features(tus(frames["tus"]), frames["tx_rates"], reps).set_index("tu_id")  # type: ignore[arg-type]


def terminator_downstream_with_partner(
    frames: Mapping[str, pd.DataFrame | None],
) -> dict[str, float]:
    """``terminator_conf_downstream`` for genes whose next gene in transcription order is a
    same-strand adjacent gene (a pair partner): then "the next same-strand gene" is unambiguous.
    NaN = NA."""
    from giap.annotate.terminators import terminator_conf_downstream

    gs = {g.locus_tag: g for g in genes(frames["genes"])}  # type: ignore[arg-type]
    reps = replicons(frames["replicons"])  # type: ignore[arg-type]
    terms = terminators(frames["terminators"])  # type: ignore[arg-type]
    out = {}
    for p in pairs(frames["pairs_raw"]):  # type: ignore[arg-type]
        a, b = gs[p.gene_a], gs[p.gene_b]
        up, down = (a, b) if a.interval.strand == "+" else (b, a)
        v = terminator_conf_downstream(up, down, terms, reps[up.interval.contig])
        out[up.locus_tag] = float("nan") if v is None else float(v)
    return out


def expected_annotation_files(
    chain_inputs: Mapping[str, Path], out: Path, yes: float = 0.5
) -> dict[str, Path]:
    """The §9.3 TSVs and the three BEDs written by the merged writers from the producer files."""
    from giap.annotate.operons import write_operon_pairs_tsv, write_tus_tsv
    from giap.annotate.promoters import promoters_bed
    from giap.annotate.regions import regions_bed, write_regions_tsv
    from giap.annotate.terminators import terminators_bed
    from giap.io.bed import write_bed

    out.mkdir(parents=True, exist_ok=True)
    reps = replicons(pd.read_parquet(chain_inputs["replicons"]))
    write_operon_pairs_tsv(
        out / "operon_pairs.tsv", pairs(pd.read_parquet(chain_inputs["pairs_raw"]))
    )
    write_tus_tsv(out / "tus.tsv", tus(pd.read_parquet(chain_inputs["tus"])))
    regs = regions(pd.read_parquet(chain_inputs["regions"]))
    write_regions_tsv(out / "regions.tsv", regs)
    write_bed(out / "regions.bed", regions_bed(regs, reps), reps)
    pw = pd.read_parquet(chain_inputs["promoter_windows"])
    write_bed(out / "promoters.bed", promoters_bed(pw, reps, yes), reps)
    terms = terminators(pd.read_parquet(chain_inputs["terminators"]))
    write_bed(out / "terminators.bed", terminators_bed(terms, reps), reps)
    return {p.name: p for p in out.iterdir()}


def promoter_record_scores(
    windows: pd.DataFrame, records: Any, reps: Mapping[str, Replicon], yes: float = 0.5
) -> dict[str, float]:
    """BED name -> the tool's own score of a merged ``promoters.bed`` record: the maximum
    Promotech score of the windows (>= ``yes``, 40 nt, M2-13 window placement) inside it."""
    from giap.annotate.promoters import WINDOW_NT

    out = {}
    for r in records:
        iv = r.interval
        rep = reps[iv.contig]
        sub = windows[(windows["contig"] == iv.contig) & (windows["strand"] == iv.strand)]
        best = -math.inf
        for p, x in zip(sub["end3p"], sub["score"], strict=True):
            if x < yes:
                continue
            lo = int(p) - WINDOW_NT if iv.strand == "+" else int(p) - 1
            if rep.topology == "circular":
                lo %= rep.length
            shifts = (0, rep.length) if rep.topology == "circular" else (0,)
            if any(iv.start <= lo + s and lo + s + WINDOW_NT <= iv.end for s in shifts):
                best = max(best, float(x))
        assert best > -math.inf, r.name
        out[r.name] = best
    return out


def isna(x: Any) -> bool:
    return x is None or (isinstance(x, float) and math.isnan(x)) or x is pd.NA


def same(a: Any, b: Any) -> bool:
    """Equal values, NA == NA, floats compared with a 1e-12 relative tolerance."""
    if isna(a) or isna(b):
        return isna(a) and isna(b)
    if isinstance(a, float) or isinstance(b, float):
        return math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-12)
    return bool(a == b)
