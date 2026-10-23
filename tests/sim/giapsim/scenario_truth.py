"""Truth tables of the M7-12 scenarios: variant sites, indel sites and novel-ORF denominators.

docs/specs/M7-12-sim-variant-novel-scenarios.md rules 1-3 (lead decisions 2, 7 and 8). Written
only with ``SimConfig.m7_truth`` (and RBSs planted only with ``SimConfig.plant_novel_rbs``), so
every pre-M7-12 config keeps its outputs byte for byte.

* ``truth/alt_sites.tsv``: one row per codon holding an SNV of ``truth/strain_differences.tsv``
  inside a gene, taken in the **sample** gene's reading frame and lifted to the reference
  (decision 11): ``codon_start`` = 1-based lowest coordinate of ``cs_sample + pos_reference -
  pos_sample``, ``ref_aa`` / ``alt_aa`` the code-11 translations of that codon on the reference /
  sample; rows synonymous or with a stop on either side are left out. ``frameshift_upstream`` is
  true when a frameshifting indel of the same gene lies upstream (transcription direction) of the
  SNV, so the reference frame there differs from the annotated one. Rows in (contig order,
  codon_start, strand) order, one per codon.
* ``truth/indel_sites.tsv``: every indel inside a gene, columns copied from
  ``strain_differences.tsv`` plus ``frameshift`` (length change not a multiple of 3) and
  ``junction`` (``pos_reference + 1``: the first reference base after the VCF anchor), so M7-13
  can stratify or mask; decision 2 keeps them out of ``alt_sites``.
* ``truth/novel_orfs.tsv``: the seven M3-18 columns, then ``stop`` (1-based lowest base of the
  stop codon, decision 16), ``true_start`` (1-based first base of the start codon in
  transcription direction), ``n_true_peptides`` and ``rbs_planted``; the N-terminal extensions
  follow the novel ORFs as class ``nterm_extension`` rows (true start to stop, true protein).
* ``truth/novel_orfs.json`` (lead decision 26), only when novel ORFs or N-terminal extensions are
  requested: ``{"requested": {class: int}, "realised": {class: int}, "shortfall": {class: int}}``
  over the four classes of ``NOVEL_CLASSES``; ``realised`` is the class's row count in
  ``truth/novel_orfs.tsv`` and ``shortfall = requested - realised`` (intergenic planting stops at
  the number of usable gaps, decision 15).

``n_true_peptides`` (rule 2, decision 7) counts distinct I/L-folded correct database peptides of
``truth/peptides.tsv`` (``is_wrong = false``, not decoys) covering the object: for a site, the
gene's peptides whose sample interval holds the SNV; for an intergenic/out-of-frame/antisense
ORF, the peptides placed in it; for an N-terminal extension, the gene's peptides with at least
one residue upstream of the annotated start. Coordinates of novel ORFs are on ``genome.fasta``
(the M7 strain scenario plants none).
"""

from __future__ import annotations

import csv
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from giapsim.codons import CODE_11, codons_of, revcomp, translate_cds
from giapsim.features import STARTS, Feat
from giapsim.signals import RBS, clean_rbs, put

__all__ = [
    "ALT_COLUMNS",
    "INDEL_COLUMNS",
    "NOVEL_CLASSES",
    "NOVEL_COLUMNS",
    "novel_report",
    "plant_novel_rbs",
    "strong_tirs",
    "write_m7_truth",
]

ALT_COLUMNS = ["contig", "codon_start", "strand", "gene", "ref_aa", "alt_aa", "n_true_peptides",
               "frameshift_upstream"]  # fmt: skip
INDEL_COLUMNS = ["contig", "pos_sample", "pos_reference", "ref", "alt", "type", "gene", "aa_change",
                 "frameshift", "junction"]  # fmt: skip
NOVEL_COLUMNS = ["orf_id", "class", "contig", "start", "end", "strand", "aa",
                 "stop", "true_start", "n_true_peptides", "rbs_planted"]  # fmt: skip
NOVEL_CLASSES = ("intergenic", "cds_outframe", "antisense", "nterm_extension")

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


def strong_tirs(
    rng: np.random.Generator,
    rows: list[list[Any]],
    keys: set[Key],
    tir: tuple[float, float],
) -> None:
    """Redraw ``tir_ostir`` of the ``rbs_novel`` rows at ``keys`` from lognormal(``tir``)
    (``SimConfig.planted_tir``, decision 14). One draw per key whatever the arm, so changing it
    changes only these values."""
    at = {(r[0], r[1], int(r[2]) - 1): i for i, r in enumerate(rows)}
    for k in sorted(keys):
        if k in at:
            rows[at[k]][3] = round(float(rng.lognormal(*tir)), 4)


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


def _sample_codon(g: Feat, p0: int, n: int) -> int:
    """0-based lowest coordinate (unwrapped, may exceed ``n``) of the codon of sample gene ``g``
    holding 0-based ``p0`` (contig length ``n``), in ``g``'s annotated frame."""
    if p0 < g.start:
        p0 += n
    if g.strand == "+":
        return g.start + 3 * ((p0 - g.start) // 3)
    return g.end - 3 * ((g.end - 1 - p0) // 3 + 1)


def _codon(seq: str, c0: int, strand: str) -> str:
    n = len(seq)
    nt = "".join(seq[(c0 + i) % n] for i in range(3))
    return revcomp(nt) if strand == "-" else nt


def _residue(codon: str, first: bool) -> str:
    """Code-11 residue; a valid start codon at the annotated start reads M, as in
    ``strain_differences.tsv`` ``aa_change``."""
    return "M" if first and codon in STARTS else CODE_11[codon]


def _is_indel(r: Sequence[Any]) -> bool:
    return r[5] in ("ins", "del")


def _shift(r: Sequence[Any]) -> int:
    return len(str(r[4])) - len(str(r[3]))


def alt_site_rows(
    strain_rows: Sequence[Sequence[Any]],
    genes: Sequence[Feat],
    sample: Mapping[str, str],
    ref: Mapping[str, str],
    peps: Mapping[str, list[tuple[str, int, int, str]]],
) -> list[list[Any]]:
    """``truth/alt_sites.tsv`` rows (decision 11): each SNV codon in the sample gene's frame
    (``genes`` = sample features), lifted to the reference by ``pos_reference - pos_sample``."""
    by_locus = {g.locus: g for g in genes}
    order = {c: k for k, c in enumerate(ref)}
    fs_at: dict[str, list[int]] = {}  # gene -> 0-based sample positions of frameshifting indels
    for r in strain_rows:
        if _is_indel(r) and r[6] and _shift(r) % 3:
            fs_at.setdefault(str(r[6]), []).append(int(r[1]) - 1)
    out: dict[tuple[str, int, str, str], list[Any]] = {}
    for cid, pos_s, pos_r, _ref, _alt, typ, gene, _change in strain_rows:
        if typ != "snv" or not gene:
            continue
        g, ns, nr = by_locus[str(gene)], len(sample[cid]), len(ref[cid])
        ps, pr = int(pos_s) - 1, int(pos_r) - 1
        cs = _sample_codon(g, ps, ns)
        rc = (cs + (pr - ps)) % nr
        first = cs == (g.start if g.strand == "+" else g.end - 3)  # initiator reads M (strain)
        ra = _residue(_codon(ref[cid], rc, g.strand), first)
        aa = _residue(_codon(sample[cid], cs, g.strand), first)
        if ra == aa or "*" in (ra, aa):
            continue
        key = (str(cid), rc + 1, g.strand, str(gene))
        if key in out:
            continue  # one row per codon
        plus = g.strand == "+"
        fs = any((q < ps) if plus else (q > ps) for q in fs_at.get(str(gene), ()))
        hits = {p for c, s1, e1, p in peps.get(str(gene), ())
                if c == cid and _covers(s1, e1, ps, ns)}  # fmt: skip
        out[key] = [cid, rc + 1, g.strand, gene, ra, aa, len(hits), "true" if fs else "false"]
    return sorted(out.values(), key=lambda r: (order[r[0]], r[1], r[2]))


def indel_site_rows(strain_rows: Sequence[Sequence[Any]]) -> list[list[Any]]:
    """``truth/indel_sites.tsv`` rows: every indel in a gene, ``frameshift`` and ``junction``."""
    return [[*r, "true" if _shift(r) % 3 else "false", int(r[2]) + 1]
            for r in strain_rows if _is_indel(r) and r[6]]  # fmt: skip


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
        stop = int(e0) - 2 if strand == "+" else int(s1)  # 1-based (decision 16)
        n = len({p for *_, p in peps.get(str(orf_id), ())})
        rows.append([orf_id, cls, cid, s1, e0, strand, aa, stop, ts, n, (cid, strand, ts - 1)])
    for f in _nterm_feats(feats):
        s0, e0 = f.span()
        nt = text[f.contig][s0:e0]
        nt = nt if f.strand == "+" else revcomp(nt)
        plus, a1, b1 = f.strand == "+", f.start + 1, f.end  # annotated CDS, 1-based closed
        n = len({p for c, s1, e1, p in peps.get(f.locus, ())
                 if c == f.contig and (s1 < a1 if plus else e1 > b1)})  # fmt: skip
        stop = e0 - 2 if f.strand == "+" else s0 + 1
        assert f.true_start is not None
        rows.append([f"NOVEL_{len(rows) + 1:04d}", "nterm_extension", f.contig, s0 + 1, e0,
                     f.strand, translate_cds(codons_of(nt)), stop, f.true_start, n,
                     (f.contig, f.strand, f.true_start - 1)])  # fmt: skip
    for r in rows:
        r[-1] = "true" if r[-1] in planted else "false"
    return rows


def novel_report(
    requested: Mapping[str, int], rows: Sequence[Sequence[Any]]
) -> dict[str, dict[str, int]]:
    """Requested, realised and shortfall counts per class of ``NOVEL_CLASSES`` (decision 26);
    ``rows`` are ``truth/novel_orfs.tsv`` rows (class in column 1)."""
    req = {c: int(requested.get(c, 0)) for c in NOVEL_CLASSES}
    got = {c: sum(1 for r in rows if r[1] == c) for c in NOVEL_CLASSES}
    return {"requested": req, "realised": got,
            "shortfall": {c: req[c] - got[c] for c in NOVEL_CLASSES}}  # fmt: skip


def write_m7_truth(
    out: Path,
    novel: Sequence[Sequence[Any]],
    feats: Sequence[Feat],
    text: Mapping[str, str],
    planted: set[Key],
    strain: tuple[Sequence[Sequence[Any]], Mapping[str, str]] | None,
    requested: Mapping[str, int] | None = None,
) -> dict[str, Path]:
    """Write the M7-12 truth tables (after proteomics); ``strain`` = (strain_differences rows,
    reference sequences) or None; ``requested`` = per-class novel counts asked for (the
    ``truth/novel_orfs.json`` report is written when any is positive). Returns the truth paths
    by key."""
    peps = _correct(out)
    paths = {"alt_sites": out / "truth" / "alt_sites.tsv", "novel_orfs": out / "truth" /
             "novel_orfs.tsv"}  # fmt: skip
    alt: list[list[Any]] = []
    if strain is not None:
        rows, ref = strain
        alt = alt_site_rows(rows, feats, text, ref, peps)
        paths["indel_sites"] = out / "truth" / "indel_sites.tsv"
        _tsv(paths["indel_sites"], INDEL_COLUMNS, indel_site_rows(rows))
    _tsv(paths["alt_sites"], ALT_COLUMNS, alt)
    nrows = novel_rows(novel, feats, text, planted, peps)
    _tsv(paths["novel_orfs"], NOVEL_COLUMNS, nrows)
    if requested is not None and any(v > 0 for v in requested.values()):
        paths["novel_report"] = out / "truth" / "novel_orfs.json"
        rep = novel_report(requested, nrows)
        paths["novel_report"].write_text(json.dumps(rep, indent=2) + "\n", encoding="utf-8")
    return paths
