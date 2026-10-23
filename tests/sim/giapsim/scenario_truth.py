"""Truth tables of the M7-12 scenarios: variant sites, indel sites and novel-ORF denominators.

docs/specs/M7-12-sim-variant-novel-scenarios.md rules 1-3 (lead decisions 2, 7 and 8). Written
only with ``SimConfig.m7_truth`` (and RBSs planted only with ``SimConfig.plant_novel_rbs``), so
every pre-M7-12 config keeps its outputs byte for byte.

* ``truth/alt_sites.tsv``: one row per non-synonymous SNV of ``truth/strain_differences.tsv``
  whose residues are both sense (nonsense SNVs are left out), on the **reference** codon (the
  genome GIAP receives): ``codon_start`` is the lowest 1-based coordinate of the codon, ``strand``
  and ``gene`` the lifted gene's, ``ref_aa`` / ``alt_aa`` the two letters of ``aa_change``. Rows
  in (contig order, codon_start, strand) order.
* ``truth/indel_sites.tsv``: the in-frame indels inside genes, columns copied from
  ``strain_differences.tsv`` (decision 2: kept out of ``alt_sites``).
* ``truth/novel_orfs.tsv``: the seven M3-18 columns, then ``stop`` (0-based lowest base of the
  stop codon), ``true_start`` (1-based first base of the start codon in transcription direction),
  ``n_true_peptides`` and ``rbs_planted``; the N-terminal extensions follow the novel ORFs as
  class ``nterm_extension`` rows (true start to stop, true protein).

``n_true_peptides`` (rule 2, decision 7) counts distinct I/L-folded correct database peptides of
``truth/peptides.tsv`` (``is_wrong = false``, not decoys) covering the object: for a site, the
gene's peptides whose sample interval holds the SNV; for an intergenic/out-of-frame/antisense
ORF, the peptides placed in it; for an N-terminal extension, the gene's peptides with at least
one residue upstream of the annotated start. Coordinates of novel ORFs are on ``genome.fasta``
(the M7 strain scenario plants none).
"""

from __future__ import annotations

import csv
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from giapsim.codons import codons_of, revcomp, translate_cds
from giapsim.features import Feat
from giapsim.signals import RBS, clean_rbs, put

__all__ = [
    "ALT_COLUMNS",
    "INDEL_COLUMNS",
    "NOVEL_COLUMNS",
    "PLANTED_TIR",
    "plant_novel_rbs",
    "strong_tirs",
    "write_m7_truth",
]

ALT_COLUMNS = ["contig", "codon_start", "strand", "gene", "ref_aa", "alt_aa", "n_true_peptides"]
INDEL_COLUMNS = ["contig", "pos_sample", "pos_reference", "ref", "alt", "type", "gene", "aa_change"]
NOVEL_COLUMNS = ["orf_id", "class", "contig", "start", "end", "strand", "aa",
                 "stop", "true_start", "n_true_peptides", "rbs_planted"]  # fmt: skip
#: ``(mu, sigma)`` of the log TIR of planted novel starts in ``rbs_novel`` (decision 8): stronger
#: than first-round TU RBSs (lognormal(8, 0.8)) so ≥ 90 % of planted starts clear the detected
#: genes' 10 % quantile by design.
PLANTED_TIR = (9.0, 0.5)
_SNV = re.compile(r"^([A-Z*])(\d+)([A-Z*])$")

Key = tuple[str, str, int]  # contig, strand, 0-based first base of the start codon


def _tsv(path: Path, header: Sequence[str], rows: Sequence[Sequence[Any]]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def _nterm_feats(feats: Sequence[Feat]) -> list[Feat]:
    return [f for f in feats if f.kind == "cds" and f.true_start and not f.pseudo_reason]


def plant_novel_rbs(
    rng: np.random.Generator,
    seqs: dict[str, bytearray],
    novel: Sequence[Sequence[Any]],
    feats: Sequence[Feat],
) -> set[Key]:
    """Plant ``AGGAGG`` (spacer 5-10, second cores broken as in signals.py) before every planted
    intergenic ORF and every N-terminal true start; returns their start keys."""
    spans = [(str(r[2]), int(r[3]) - 1, int(r[4]), str(r[5])) for r in novel
             if r[1] == "intergenic"]  # fmt: skip
    spans += [(f.contig, *f.span(), f.strand) for f in _nterm_feats(feats)]
    keys = set()
    for contig, s0, e0, strand in spans:
        sp = int(rng.integers(5, 11))
        put(seqs[contig], s0, e0, strand, -sp - len(RBS), RBS)
        clean_rbs(seqs[contig], s0, e0, strand, -sp - len(RBS), sp)
        keys.add((contig, strand, s0 if strand == "+" else e0 - 1))
    return keys


def strong_tirs(rng: np.random.Generator, rows: list[list[Any]], keys: set[Key]) -> None:
    """Redraw ``tir_ostir`` of the ``rbs_novel`` rows at ``keys`` from ``PLANTED_TIR``."""
    at = {(r[0], r[1], int(r[2]) - 1): i for i, r in enumerate(rows)}
    for k in sorted(keys):
        if k in at:
            rows[at[k]][3] = round(float(rng.lognormal(*PLANTED_TIR)), 4)


def _correct(out: Path) -> dict[str, list[tuple[str, int, int, str]]]:
    """Correct database peptides by locus: (contig, start, end, I/L-folded peptide)."""
    path = out / "truth" / "peptides.tsv"
    by: dict[str, list[tuple[str, int, int, str]]] = {}
    if not path.exists():
        return by
    with path.open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["origin"] == "decoy" or r["is_wrong"] == "true":
                continue
            by.setdefault(r["locus_tag"], []).append(
                (r["contig"], int(r["start"]), int(r["end"]), r["peptide"].replace("I", "L"))
            )
    return by


def _covers(s1: int, e1: int, p0: int, n: int) -> bool:
    return any(s1 - 1 <= q < e1 for q in (p0, p0 + n))


def _codon_start(g: Feat, p0: int, n: int) -> int:
    """1-based lowest coordinate of ``g``'s codon holding 0-based ``p0`` (contig length ``n``)."""
    if p0 < g.start:
        p0 += n
    if g.strand == "+":
        c0 = g.start + 3 * ((p0 - g.start) // 3)
    else:
        c0 = g.end - 3 * ((g.end - 1 - p0) // 3 + 1)
    return c0 % n + 1


def alt_site_rows(
    strain_rows: Sequence[Sequence[Any]],
    lifted: Sequence[Feat],
    ref_len: Mapping[str, int],
    sample_len: Mapping[str, int],
    peps: Mapping[str, list[tuple[str, int, int, str]]],
) -> list[list[Any]]:
    by_locus = {g.locus: g for g in lifted}
    order = {c: k for k, c in enumerate(ref_len)}
    out = []
    for cid, pos_s, pos_r, _ref, _alt, typ, gene, change in strain_rows:
        m = _SNV.match(str(change))
        if typ != "snv" or not gene or m is None or "*" in (m[1], m[3]):
            continue
        cs = _codon_start(by_locus[gene], int(pos_r) - 1, ref_len[cid])
        hits = {p for c, s1, e1, p in peps.get(gene, ())
                if c == cid and _covers(s1, e1, int(pos_s) - 1, sample_len[cid])}  # fmt: skip
        out.append([cid, cs, by_locus[gene].strand, gene, m[1], m[3], len(hits)])
    out.sort(key=lambda r: (order[r[0]], r[1], r[2]))
    return out


def novel_rows(
    novel: Sequence[Sequence[Any]],
    feats: Sequence[Feat],
    text: Mapping[str, str],
    planted: set[Key],
    peps: Mapping[str, list[tuple[str, int, int, str]]],
) -> list[list[Any]]:
    rows = []
    for orf_id, cls, cid, s1, e0, strand, aa in novel:
        ts = int(s1) if strand == "+" else int(e0)
        stop = int(e0) - 3 if strand == "+" else int(s1) - 1
        n = len({p for *_, p in peps.get(str(orf_id), ())})
        rows.append([orf_id, cls, cid, s1, e0, strand, aa, stop, ts, n, (cid, strand, ts - 1)])
    for f in _nterm_feats(feats):
        s0, e0 = f.span()
        nt = text[f.contig][s0:e0]
        nt = nt if f.strand == "+" else revcomp(nt)
        plus, a1, b1 = f.strand == "+", f.start + 1, f.end  # annotated CDS, 1-based closed
        n = len({p for c, s1, e1, p in peps.get(f.locus, ())
                 if c == f.contig and (s1 < a1 if plus else e1 > b1)})  # fmt: skip
        stop = e0 - 3 if f.strand == "+" else s0
        assert f.true_start is not None
        rows.append([f"NOVEL_{len(rows) + 1:04d}", "nterm_extension", f.contig, s0 + 1, e0,
                     f.strand, translate_cds(codons_of(nt)), stop, f.true_start, n,
                     (f.contig, f.strand, f.true_start - 1)])  # fmt: skip
    for r in rows:
        r[-1] = "true" if r[-1] in planted else "false"
    return rows


def write_m7_truth(
    out: Path,
    novel: Sequence[Sequence[Any]],
    feats: Sequence[Feat],
    text: Mapping[str, str],
    planted: set[Key],
    strain: tuple[Sequence[Sequence[Any]], Sequence[Feat], Mapping[str, str]] | None,
) -> dict[str, Path]:
    """Write the M7-12 truth tables (after proteomics); ``strain`` = (strain_differences rows,
    lifted features, reference sequences) or None. Returns the truth paths by key."""
    peps = _correct(out)
    paths = {"alt_sites": out / "truth" / "alt_sites.tsv", "novel_orfs": out / "truth" /
             "novel_orfs.tsv"}  # fmt: skip
    alt: list[list[Any]] = []
    if strain is not None:
        rows, lifted, ref = strain
        sample_len = {c: len(s) for c, s in text.items()}
        alt = alt_site_rows(rows, lifted, {c: len(s) for c, s in ref.items()}, sample_len, peps)
        paths["indel_sites"] = out / "truth" / "indel_sites.tsv"
        indels = [r for r in rows if r[5] in ("ins", "del") and r[6]
                  and (len(str(r[3])) - len(str(r[4]))) % 3 == 0]  # fmt: skip
        _tsv(paths["indel_sites"], INDEL_COLUMNS, indels)
    _tsv(paths["alt_sites"], ALT_COLUMNS, alt)
    _tsv(paths["novel_orfs"], NOVEL_COLUMNS, novel_rows(novel, feats, text, planted, peps))
    return paths
