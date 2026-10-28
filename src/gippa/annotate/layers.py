"""Per-gene context layers and the ``annotation/*`` files (GIPPA_SPEC.md §7.1.8, §9.1-§9.3; M2-20).

``gene_layers`` joins the parsed layers of the annotation track into one row per gene of
``genes.parquet`` (``LAYER_COLUMNS``, the C6 ``layers.parquet`` contract). Inputs are the producer
frames (0-based half-open); a missing optional layer is ``None`` and gives NA columns.
``uniop_p_prev`` / ``uniop_p_next`` are the pair probabilities with the previous / next gene in
transcription order (genome order on ``+``, reversed on ``-``). ``gc_local_dev`` is the GC
fraction of the 10-kb window centred on the gene midpoint minus the whole-genome GC (A/C/G/T only,
case-insensitive; circular replicons wrap, linear ones truncate; lead decision 7).
``merged_gff3_features`` builds the source-tagged ``annotation/genes.gff3`` model;
``write_annotation`` writes every §7.1.8 file with the producers' own writers.
"""

from __future__ import annotations

import shutil
from collections.abc import Mapping, Sequence
from itertools import pairwise
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from gippa.annotate.operons import Tu, write_operon_pairs_tsv, write_tus_tsv
from gippa.annotate.promoter_calculator import tu_tx_features
from gippa.annotate.promoters import K_RANGE, gene_promoter_features
from gippa.annotate.regions import Region, regions_bed, write_regions_tsv
from gippa.annotate.terminators import Terminator, terminator_conf_downstream, terminators_bed
from gippa.annotate.uniop_adapter import PairProb
from gippa.config import GippaConfig
from gippa.io.bed import BedRecord, write_bed
from gippa.io.gff3 import GffFeature, write_gff3
from gippa.io.tsv import format_value
from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene

LAYER_COLUMNS: tuple[str, ...] = (
    "locus_tag", "tu_id", "tu_position", "tu_size", "uniop_p_prev", "uniop_p_next", "region_ids",
    "gc_local_dev", "promoter_score", "promoter_score_tu", "n_promoter_windows_tu",
    "tx_rate_max_tu", "tx_rate_sum_tu", "tss_distance_bp", "terminator_conf_downstream",
    "tir_ostir", "sd_dg", "mrna_dg", "spacing_dg", "standby_dg", "spacer_bp", "start_codon",
    "cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "rare_codons_2_10", "flags",
)  # fmt: skip
INT_COLUMNS = ("tu_position", "tu_size", "n_promoter_windows_tu", "tss_distance_bp",
               "spacer_bp", "rare_codons_2_10")  # fmt: skip
LIST_COLUMNS = ("region_ids", "flags")
STR_COLUMNS = ("locus_tag", "tu_id", "start_codon")
GC_WINDOW_BP = 10_000
_PROM = ("promoter_score", "promoter_score_tu", "n_promoter_windows_tu")
_RBS = ("tir_ostir", "sd_dg", "mrna_dg", "spacing_dg", "standby_dg", "spacer_bp")
_CODON = ("cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "rare_codons_2_10")
_GFF_REGION_TYPE = {
    "prophage": "prophage", "plasmid": "plasmid", "bgc": "biosynthetic_gene_cluster",
    "crispr": "CRISPR",
}  # fmt: skip


# ----------------------------------------------------------------------------- frame -> objects


def replicons_of(df: pd.DataFrame) -> dict[str, Replicon]:
    return {
        str(r.contig): Replicon(str(r.contig), int(r.length), r.topology, r.kind)
        for r in df.itertuples(index=False)
    }


def _str(v: Any) -> str | None:
    return v if isinstance(v, str) and v else None


def genes_of(df: pd.DataFrame) -> list[Gene]:
    return [
        Gene(
            locus_tag=str(r.locus_tag), feature_type=r.feature_type,
            interval=Interval(str(r.contig), int(r.start), int(r.end), r.strand),
            aa="", nt=None, start_codon=_str(getattr(r, "start_codon", None)),
            gene=_str(getattr(r, "gene", None)), product=_str(getattr(r, "product", None)) or "",
            partial=bool(getattr(r, "partial", False)),
        )
        for r in df.itertuples(index=False)
    ]  # fmt: skip


def tus_of(df: pd.DataFrame) -> list[Tu]:
    return [
        Tu(
            str(r.tu_id), str(r.contig),
            Interval(str(r.contig), int(r.start), int(r.end), r.strand), r.strand,
            tuple(str(g) for g in r.genes),
            Interval(str(r.contig), int(r.upstream_start), int(r.upstream_end), r.strand),
        )
        for r in df.itertuples(index=False)
    ]  # fmt: skip


def pairs_of(df: pd.DataFrame) -> list[PairProb]:
    cols = df[["gene_a", "gene_b", "uniop_p", "source"]].itertuples(index=False)
    return [PairProb(str(a), str(b), float(p), s) for a, b, p, s in cols]


def terminators_of(df: pd.DataFrame) -> list[Terminator]:
    out = []
    for r in df.itertuples(index=False):
        dg = getattr(r, "hairpin_dg", None)
        out.append(
            Terminator(
                Interval(str(r.contig), int(r.start), int(r.end), r.strand), float(r.confidence),
                r.source, None if dg is None or pd.isna(dg) else float(dg),
            )
        )  # fmt: skip
    return out


def regions_of(df: pd.DataFrame) -> list[Region]:
    out = []
    for r in df.itertuples(index=False):
        has = not (pd.isna(r.start) or pd.isna(r.end))
        iv = Interval(str(r.contig), int(r.start), int(r.end), ".") if has else None
        score = None if r.score is None or pd.isna(r.score) else float(r.score)
        genes = tuple(r.genes if r.genes is not None else ())
        notes = r.notes if isinstance(r.notes, str) else ""
        out.append(Region(str(r.region_id), r.region_type, str(r.contig), iv, r.source_tool,
                          score, genes, notes))  # fmt: skip
    return out


def assign_tus(
    tx: pd.DataFrame, tus: Sequence[Tu], replicons: Mapping[str, Replicon]
) -> pd.DataFrame:
    """Precomputed ``tx_rates`` rows (no ``tu_id``) → one row per TU whose upstream region holds
    the TSS on the same strand (the ``tx_layer`` shape)."""
    rows = []
    for tu in tus:
        up, rep = tu.upstream, replicons[tu.contig]
        n = up.end - up.start
        sel = tx[(tx["contig"] == tu.contig) & (tx["strand"] == tu.strand)]
        off = sel["tss"] - up.start
        if rep.topology == "circular":
            off = off % rep.length
        hit = (off >= 0) & (off < n)
        for t, rate in zip(sel["tss"][hit], sel["tx_rate"][hit], strict=True):
            rows.append((tu.contig, tu.strand, int(t), float(rate), tu.tu_id))
    return pd.DataFrame(rows, columns=["contig", "strand", "tss", "tx_rate", "tu_id"])


# ----------------------------------------------------------------------------- gc_local_dev


def _counts(seq: str) -> tuple[np.ndarray, np.ndarray]:
    """Prefix sums of A/C/G/T and G/C bases (case-insensitive)."""
    b = np.frombuffer(seq.upper().encode("ascii", "replace"), dtype=np.uint8)
    gc = (b == ord("G")) | (b == ord("C"))
    acgt = gc | (b == ord("A")) | (b == ord("T"))
    z = np.zeros(1, dtype=np.int64)
    return np.concatenate([z, np.cumsum(acgt)]), np.concatenate([z, np.cumsum(gc)])


def gc_local_dev(
    genes: pd.DataFrame,
    replicons: pd.DataFrame,
    genome: Mapping[str, str],
    window: int = GC_WINDOW_BP,
) -> list[float]:
    """Window GC minus whole-genome GC per gene row (NaN when the window has no A/C/G/T)."""
    sums = {c: _counts(s) for c, s in genome.items()}
    tot = sum(int(a[-1]) for a, _ in sums.values())
    genome_gc = sum(int(g[-1]) for _, g in sums.values()) / tot if tot else float("nan")
    reps = replicons_of(replicons)
    half = window // 2
    out = []
    for c, s, e in zip(genes["contig"], genes["start"], genes["end"], strict=True):
        acgt, gc = sums[str(c)]
        rep = reps[str(c)]
        n, m = rep.length, (int(s) + int(e)) // 2
        if rep.topology == "circular" and n <= window:
            spans = [(0, n)]
        elif rep.topology == "circular":
            lo = (m - half) % n
            hi = lo + window
            spans = [(lo, min(hi, n))] + ([(0, hi - n)] if hi > n else [])
        else:
            spans = [(max(0, m - half), min(n, m - half + window))]
        d = sum(int(acgt[b] - acgt[a]) for a, b in spans if b > a)
        g = sum(int(gc[b] - gc[a]) for a, b in spans if b > a)
        out.append(g / d - genome_gc if d else float("nan"))
    return out


# ----------------------------------------------------------------------------- per-gene layers


def _next_same_strand(
    genes: Sequence[Gene], replicons: Mapping[str, Replicon]
) -> dict[str, Gene | None]:
    """Next gene on the same contig and strand in transcription order (wrapping on circular)."""
    groups: dict[tuple[str, str], list[Gene]] = {}
    for g in sorted(genes, key=lambda g: (g.interval.start, g.interval.end, g.locus_tag)):
        groups.setdefault((g.interval.contig, g.interval.strand), []).append(g)
    out: dict[str, Gene | None] = {}
    for (contig, strand), gs in groups.items():
        ring = replicons[contig].topology == "circular" and len(gs) > 1
        step = 1 if strand == "+" else -1
        for k, g in enumerate(gs):
            j = k + step
            out[g.locus_tag] = gs[j % len(gs)] if ring or 0 <= j < len(gs) else None
    return out


def _by_tag(df: pd.DataFrame | None, cols: Sequence[str]) -> dict[str, tuple[Any, ...]]:
    if df is None:
        return {}
    return {str(r[0]): tuple(r[1:]) for r in df[["locus_tag", *cols]].itertuples(index=False)}


def _na(v: Any) -> Any:
    return None if v is None or (not isinstance(v, str) and pd.isna(v)) else v


def gene_layers(
    genes: pd.DataFrame,
    tus: pd.DataFrame,
    pairs: pd.DataFrame,
    promoters: pd.DataFrame | None,
    tx: pd.DataFrame | None,
    terms: pd.DataFrame | None,
    rbs: pd.DataFrame | None,
    codon: pd.DataFrame | None,
    region_genes: pd.DataFrame | None,
    replicons: pd.DataFrame,
    *,
    genome: Mapping[str, str] | None = None,
    yes: float = 0.5,
    k_range: tuple[int, int] = K_RANGE,
) -> pd.DataFrame:
    """One ``LAYER_COLUMNS`` row per gene of ``genes`` (same order)."""
    reps = replicons_of(replicons)
    gs = genes_of(genes)
    tu_objs = tus_of(tus)
    pos = {t: (tu.tu_id, i, len(tu.genes)) for tu in tu_objs for i, t in enumerate(tu.genes, 1)}
    as_a = dict(zip(pairs["gene_a"].astype(str), pairs["uniop_p"].astype(float), strict=True))
    as_b = dict(zip(pairs["gene_b"].astype(str), pairs["uniop_p"].astype(float), strict=True))
    nan = float("nan")
    cols: dict[str, list[Any]] = {c: [] for c in LAYER_COLUMNS}
    rid: dict[str, set[str]] = {}
    if region_genes is not None:
        for r, t in zip(region_genes["region_id"], region_genes["locus_tag"], strict=True):
            rid.setdefault(str(t), set()).add(str(r))
    prom = {} if promoters is None else _by_tag(
        gene_promoter_features(gs, tu_objs, promoters, yes, k_range, replicons=reps), _PROM
    )  # fmt: skip
    txf = {} if tx is None else {
        str(r[0]): tuple(r[1:])
        for r in tu_tx_features(tu_objs, tx, reps).itertuples(index=False)
    }  # fmt: skip
    nxt = _next_same_strand(gs, reps) if terms is not None else {}
    term_objs = terminators_of(terms) if terms is not None else []
    rb, cd = _by_tag(rbs, _RBS), _by_tag(codon, _CODON)
    cai_ref = set()
    if codon is not None and "cai_reference" in codon:
        cai_ref = set(codon.loc[codon["cai_reference"].fillna(False).astype(bool), "locus_tag"])
    gc = gc_local_dev(genes, replicons, genome) if genome is not None else [nan] * len(gs)
    crosses = genes["crosses_origin"] if "crosses_origin" in genes else [False] * len(gs)
    for g, gcd, cross in zip(gs, gc, crosses, strict=True):
        t = g.locus_tag
        tu_id, k, size = pos[t]
        up, down = (as_b.get(t, nan), as_a.get(t, nan))
        prev_p, next_p = (up, down) if g.interval.strand == "+" else (down, up)
        tc = None
        if terms is not None:
            tc = terminator_conf_downstream(g, nxt[t], term_objs, reps[g.interval.contig])
        tx_row = txf.get(tu_id, (None, None, None))
        flags = [f for f, on in (("partial", g.partial), ("crosses_origin", bool(cross)),
                                 ("pseudogene", g.feature_type == "pseudogene"),
                                 ("cai_reference", t in cai_ref)) if on]  # fmt: skip
        vals = (t, tu_id, k, size, prev_p, next_p, sorted(rid.get(t, ())), gcd,
                *prom.get(t, (None,) * 3), *tx_row, tc, *rb.get(t, (None,) * 6), g.start_codon,
                *cd.get(t, (None,) * 8), flags)  # fmt: skip
        for c, v in zip(LAYER_COLUMNS, vals, strict=True):
            cols[c].append(v if c in LIST_COLUMNS else _na(v))
    df = pd.DataFrame(cols, columns=list(LAYER_COLUMNS))
    for c in LAYER_COLUMNS:
        if c in INT_COLUMNS:
            df[c] = df[c].astype("Int64")
        elif c not in STR_COLUMNS and c not in LIST_COLUMNS:
            df[c] = df[c].astype("float64")
    return df


def layer_flags(cfg: GippaConfig, present: set[str]) -> tuple[str, ...]:
    """``layer_missing:<layer>`` run flags for the ``context_layers`` input keys ``present``."""
    missing = {
        k for k in ("tx_rates", "terminators", "rbs", "codon", "regions") if k not in present
    }
    if not present & {"promoter_windows", "promoters"}:
        missing.add("promoters")
    ann = cfg.annotation
    if not ann.promoter_calculator.enabled and ann.precomputed.tx_rates is None:
        missing.add("tx_rates")
    if "mobile" not in cfg.modules:
        missing.add("genomad")
    if "bgc" not in cfg.modules:
        missing.add("antismash")
    return tuple(sorted(f"layer_missing:{k}" for k in missing))


# ----------------------------------------------------------------------------- genes.gff3


def tu_sources(tus: Sequence[Tu], pairs: Sequence[PairProb]) -> dict[str, str]:
    """GFF3 source per TU: ``distance_model`` when a linking pair (or, for a singleton, every
    pair of the run) comes from the distance model, else ``UniOP``. Pairs are keyed in genome
    order, so a ``-`` TU (genes in transcription order) is walked in reverse."""
    src = {(p.gene_a, p.gene_b): p.source for p in pairs}
    all_dm = bool(pairs) and all(p.source == "distance_model" for p in pairs)
    out = {}
    for tu in tus:
        gs = list(tu.genes) if tu.strand == "+" else list(reversed(tu.genes))
        links = [src.get(ab) for ab in pairwise(gs)]
        dm = "distance_model" in links if links else all_dm
        out[tu.tu_id] = "distance_model" if dm else "UniOP"
    return out


def terminator_records(
    terms: Sequence[Terminator], replicons: Mapping[str, Replicon]
) -> list[tuple[BedRecord, Terminator]]:
    """``terminators.bed`` records paired with their calls (the BED order)."""
    rank = {c: k for k, c in enumerate(replicons)}
    ordered = sorted(terms, key=lambda t: (rank[t.interval.contig], t.interval.start,
                                           t.interval.end, t.interval.strand, t.source,
                                           t.confidence))  # fmt: skip
    return list(zip(terminators_bed(ordered, replicons), ordered, strict=True))


def merged_gff3_features(
    genes: Sequence[Gene],
    tus: Sequence[Tu],
    promoters: Sequence[tuple[BedRecord, float]],
    terms: Sequence[tuple[BedRecord, Terminator]],
    regions: Sequence[Region],
    *,
    tu_source: Mapping[str, str] | None = None,
) -> list[GffFeature]:
    """Genes (``Bakta``), TUs (``operon``), promoters (BED record, tool score), terminators and
    regions as GFF3 features. Bakta's IDs are kept; GIPPA's carry their kind (``operon:``,
    ``promoter:``, ``terminator:``, ``region:``) so IDs stay unique per file (§9.5)."""
    out = []
    for g in genes:
        attrs: dict[str, str | list[str]] = {"ID": g.locus_tag, "locus_tag": g.locus_tag}
        if g.gene:
            attrs["gene"] = g.gene
        if g.product:
            attrs["product"] = g.product
        attrs["feature_type"] = g.feature_type
        cds = g.feature_type != "pseudogene"
        out.append(GffFeature(g.interval, "Bakta", "CDS" if cds else "pseudogene",
                              phase=0 if cds else None, attributes=attrs))  # fmt: skip
    for tu in tus:
        up = f"{tu.upstream.start + 1}-{tu.upstream.end}"
        src = (tu_source or {}).get(tu.tu_id, "UniOP")
        out.append(
            GffFeature(
                tu.interval,
                src,
                "operon",
                attributes={"ID": f"operon:{tu.tu_id}", "genes": list(tu.genes), "upstream": up},
            )
        )
    for b, x in promoters:
        out.append(GffFeature(b.interval, "Promotech", "promoter", score=x,
                              attributes={"ID": f"promoter:{b.name}"}))  # fmt: skip
    for b, t in terms:
        out.append(
            GffFeature(
                b.interval,
                t.source,
                "terminator",
                attributes={"ID": f"terminator:{b.name}", "confidence": format_value(t.confidence)},
            )
        )
    for r in regions:
        if r.interval is None or r.region_type not in _GFF_REGION_TYPE:
            continue
        attrs = {"ID": f"region:{r.region_id}", **({"genes": list(r.genes)} if r.genes else {})}
        out.append(GffFeature(r.interval, r.source_tool, _GFF_REGION_TYPE[r.region_type],
                              attributes=attrs))  # fmt: skip
    return out


# ----------------------------------------------------------------------------- files


def copy_native(tool_dirs: Mapping[str, Path], dest: Path) -> None:
    """Byte copies of native tool outputs to ``dest/<tool>/`` (no gzip, no rename)."""
    for tool, src in sorted(tool_dirs.items()):
        if src.is_dir():
            shutil.copytree(src, dest / tool, dirs_exist_ok=True)
        else:
            (dest / tool).mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dest / tool / src.name)


def write_annotation(
    outdir: Path,
    *,
    replicons: Mapping[str, Replicon],
    features: Sequence[GffFeature],
    pairs: Sequence[PairProb],
    tus: Sequence[Tu],
    regions: Sequence[Region] | None = None,
    promoters_bed: Sequence[BedRecord] | None = None,
    term_bed: Sequence[BedRecord] | None = None,
    copies: Mapping[str, Path] | None = None,
    tool_dirs: Mapping[str, Path] | None = None,
) -> None:
    """Write every §7.1.8 file of ``annotation/`` (``copies``: file name → byte-copied source)."""
    outdir.mkdir(parents=True, exist_ok=True)
    write_gff3(outdir / "genes.gff3", features, list(replicons.values()))
    write_operon_pairs_tsv(outdir / "operon_pairs.tsv", pairs)
    write_tus_tsv(outdir / "tus.tsv", tus)
    if regions is not None:
        write_regions_tsv(outdir / "regions.tsv", regions)
        write_bed(outdir / "regions.bed", regions_bed(regions, replicons), replicons)
    if promoters_bed is not None:
        write_bed(outdir / "promoters.bed", promoters_bed, replicons)
    if term_bed is not None:
        write_bed(outdir / "terminators.bed", term_bed, replicons)
    for name, src in sorted((copies or {}).items()):
        shutil.copyfile(src, outdir / name)
    if tool_dirs:
        copy_native(tool_dirs, outdir / "tools")
