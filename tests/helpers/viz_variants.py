"""Small synthetic ``OUT/`` trees for the V12 ``variants`` tests (M9-18).

Every builder takes 1-based closed coordinates, as the published TSVs carry them; ``RunView``
converts ``peptide_hits.start`` / ``genes.tsv`` ``start`` to 0-based on read. ``aa_pos`` and
``window_start_aa`` are 1-based protein positions. Tables go through
``helpers.viz_out_tree.write_schema_table`` (schema header; absent cells ``NA``). Nothing here
touches ``$GIAP_DB_DIR`` or the network.

Fixture genome (``chr1`` 10,000 bp circular, ``chr2`` 5,000 bp linear):

- ``GSIM_00105`` ``+`` [1001, 1903]: 300 aa; residues 210-219 are ``LVAEQTRGLK`` (Q at 214).
  Residue ``p`` is the codon at ``1001 + 3 (p - 1)``.
- ``GSIM_00200`` ``-`` [5001, 5600]: 199 aa; residue ``p`` is the codon whose lowest base is
  ``5601 - 3 p`` (the protein is in transcription order, high to low coordinates).
- ``GSIM_00900`` ``+`` [9901, 10299] (crosses the origin, end > 10,000): 132 aa; residue ``p``
  is the unwrapped codon ``9901 + 3 (p - 1)``.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from helpers.viz_locus import hit, preflight
from helpers.viz_out_tree import gene_row, make_out, write_schema_table

AAS = "ACDEFGHIKLMNPQRSTVWY"
L_CHR1 = 10_000
REPLICONS = (
    {"id": "chr1", "kind": "chromosome", "length": L_CHR1, "topology": "circular"},
    {"id": "chr2", "kind": "unknown", "length": 5_000, "topology": "linear"},
)
#: the §9.4 / M7-05 flag order
FLAG_ORDER = (
    "mass_equivalent", "snv", "multi_nt", "strain_difference", "database_derived",
    "start_codon", "stop_codon", "vcf_unrepresentable",
)  # fmt: skip


def protein(n: int, seed: int) -> str:
    """A deterministic random protein of ``n`` residues (no stop)."""
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list(AAS), size=n))


def _with(seq: str, pos1: int, residues: str) -> str:
    """``seq`` with ``residues`` written from 1-based position ``pos1``."""
    i = pos1 - 1
    return seq[:i] + residues + seq[i + len(residues) :]


P105 = _with(protein(300, 105), 210, "LVAEQTRGLK")
P200 = _with(protein(199, 200), 50, "D")
P900 = _with(protein(132, 900), 35, "S")
PROTEINS: Mapping[str, str] = {"GSIM_00105": P105, "GSIM_00200": P200, "GSIM_00900": P900}


def plus_codon(gene_start: int, pos: int) -> int:
    """1-based lowest base of residue ``pos`` of a ``+`` gene starting at ``gene_start``."""
    return gene_start + 3 * (pos - 1)


def minus_codon(gene_end: int, pos: int) -> int:
    """1-based lowest base of residue ``pos`` of a ``-`` gene ending at ``gene_end``."""
    return gene_end - 3 * pos + 1


def plus_hit(peptide: str, gene_start: int, first_res: int, contig: str = "chr1") -> dict:
    """A ``+`` placement whose first residue is gene residue ``first_res``."""
    start = plus_codon(gene_start, first_res)
    return hit(peptide, contig, start, start + 3 * len(peptide) - 1, "+")


def minus_hit(peptide: str, gene_end: int, first_res: int, contig: str = "chr1") -> dict:
    """A ``-`` placement whose first residue (highest codon) is gene residue ``first_res``."""
    last_res = first_res + len(peptide) - 1
    return hit(peptide, contig, minus_codon(gene_end, last_res), gene_end - 3 * (first_res - 1),
               "-")  # fmt: skip


def variant(
    vid: str,
    *,
    gene: str = "GSIM_00105",
    aa_pos: int = 214,
    ref_aa: str = "Q",
    alt_aa: str = "K",
    alt_type: str = "substitution",
    call: str = "ALT",
    confidence: str = "high",
    n_alt: int = 3,
    n_ref: int = 0,
    class_error_rate: str = "0.004",
    flags: str | None = "snv",
    contig: str = "chr1",
    strand: str = "+",
    codon_start: int | None = None,
) -> dict[str, object]:
    """One ``variants.tsv`` row; ``codon_start`` defaults to the fixture gene's codon."""
    if codon_start is None:
        codon_start = {
            "GSIM_00105": plus_codon(1001, aa_pos),
            "GSIM_00200": minus_codon(5600, aa_pos),
            "GSIM_00900": (plus_codon(9901, aa_pos) - 1) % L_CHR1 + 1,
        }.get(gene, 1)
    return {
        "variant_id": vid,
        "contig": contig,
        "codon_start": codon_start,
        "strand": strand,
        "gene": gene,
        "aa_pos": aa_pos,
        "ref_codon": "CAA",
        "ref_aa": ref_aa,
        "alt_aa": alt_aa,
        "alt_type": alt_type,
        "call": call,
        "confidence": confidence,
        "n_ref": n_ref,
        "n_alt": n_alt,
        "nt_changes": 1 if alt_type == "substitution" else None,
        "alt_codon_vcf": None,
        "class_error_rate": class_error_rate,
        "class_error_n": 40,
        "source": "database",
        "flags": flags,
    }


def fixture_genes() -> list[dict[str, object]]:
    """``genes.tsv`` rows of the three fixture genes (1-based closed)."""
    return [
        gene_row(locus_tag="GSIM_00105", contig="chr1", start=1001, end=1903, strand="+",
                 length_aa=300),
        gene_row(locus_tag="GSIM_00200", contig="chr1", start=5001, end=5600, strand="-",
                 length_aa=199),
        gene_row(locus_tag="GSIM_00900", contig="chr1", start=9901, end=10299, strand="+",
                 length_aa=132, crosses_origin=True),
    ]  # fmt: skip


def genes_parquet_frame(proteins: Mapping[str, str] = PROTEINS) -> pd.DataFrame:
    """A C2 ``genes.parquet``-shaped frame (0-based half-open) of the fixture genes."""
    coords = {
        "GSIM_00105": ("chr1", 1000, 1903, "+"),
        "GSIM_00200": ("chr1", 5000, 5600, "-"),
        "GSIM_00900": ("chr1", 9900, 10299, "+"),
    }
    rows = []
    for tag, aa in proteins.items():
        contig, start, end, strand = coords.get(tag, ("chr2", 0, 3 * len(aa) + 3, "+"))
        rows.append((tag, "CDS", contig, start, end, strand, aa))
    cols = ["locus_tag", "feature_type", "contig", "start", "end", "strand", "aa"]
    return pd.DataFrame(rows, columns=cols).astype({"start": "int64", "end": "int64"})


def window(aa: str, pos: int, flank: int = 10) -> tuple[int, str]:
    """Oracle: ``(window_start_aa, ref_window)``, residues ``pos ± flank`` clipped, 1-based."""
    lo = max(1, pos - flank)
    hi = min(len(aa), pos + flank)
    return lo, aa[lo - 1 : hi]


def windows_frame(
    sites: Sequence[Mapping[str, object]], proteins: Mapping[str, str] = PROTEINS
) -> pd.DataFrame:
    """Oracle ``variant_windows`` rows for ``sites`` whose gene is in ``proteins``."""
    rows = []
    for s in sites:
        aa = proteins.get(str(s["gene"]))
        if aa is None:
            continue
        lo, seq = window(aa, int(str(s["aa_pos"])))
        rows.append((s["variant_id"], s["gene"], int(str(s["aa_pos"])), lo, seq))
    cols = ["variant_id", "gene", "aa_pos", "window_start_aa", "ref_window"]
    return pd.DataFrame(rows, columns=cols)


def build_out(
    tmp: Path,
    variants: Sequence[Mapping[str, object]] | None,
    hits: Sequence[Mapping[str, object]] | None = None,
    windows: pd.DataFrame | None = None,
    genes: Sequence[Mapping[str, object]] | None = None,
) -> Path:
    """An ``OUT/`` with ``preflight.json`` and ``genes.tsv``; ``variants=None`` leaves
    ``variants.tsv`` out; ``windows`` is written to ``plots/support/variant_windows.tsv``."""
    from giap.io.tsv import write_tsv

    out = make_out(tmp)
    pf = out / "preflight" / "preflight.json"
    pf.parent.mkdir(parents=True, exist_ok=True)
    pf.write_text(json.dumps(preflight(REPLICONS), sort_keys=True), encoding="utf-8")
    write_schema_table(
        out / "tables" / "genes.tsv", "genes", fixture_genes() if genes is None else genes
    )
    if variants is not None:
        write_schema_table(out / "tables" / "variants.tsv", "variants", variants)
    if hits is not None:
        write_schema_table(out / "mapping" / "peptide_hits.tsv", "peptide_hits", hits)
    if windows is not None:
        sup = out / "plots" / "support"
        sup.mkdir(parents=True, exist_ok=True)
        write_tsv(
            sup / "variant_windows.tsv",
            [str(c) for c in windows.columns],
            windows.itertuples(index=False, name=None),
        )
    return out


# ---------------------------------------------------------------- the main fixture

#: VAR0001 site peptides (GSIM_00105 aa 214 Q->K)
PEP1 = "LVAEKTRGLK"  # residues 210-219, K at index 4
PEP2 = "AEKTR"  # residues 212-216, K at index 2
PEP_LONG = _with(P105[199:239], 15, "K")  # residues 200-239, K at index 14
#: VAR0002 site peptide (GSIM_00200 aa 50 D->N, - strand): residues 45-54, N at index 5
PEP_MINUS = _with(P200[44:54], 6, "N")
#: VAR0003 frameshift peptides (GSIM_00105 aa 250)
PEP_FS = "QWERTYHKLMNP"
#: VAR0004 origin-crossing peptide (GSIM_00900 aa 35 S->T): residues 31-40, T at index 4
PEP_ORIGIN = _with(P900[30:40], 5, "T")


def main_variants() -> list[dict[str, object]]:
    """VAR0001-VAR0004 ALT (panel (b) order 1, 2, 3, 4) and VAR0005 mixed."""
    return [
        variant("VAR0001"),
        variant("VAR0002", gene="GSIM_00200", aa_pos=50, ref_aa="D", alt_aa="N", strand="-",
                n_alt=2, n_ref=1, class_error_rate="0.01", flags="snv;start_codon"),
        variant("VAR0003", aa_pos=250, ref_aa=P105[249], alt_aa="fs+1", alt_type="frameshift",
                confidence="candidate", n_alt=2, class_error_rate="unestimated",
                flags="database_derived"),
        variant("VAR0004", gene="GSIM_00900", aa_pos=35, ref_aa="S", alt_aa="T",
                confidence="candidate", n_alt=1, class_error_rate="0.02", flags=None),
        variant("VAR0005", aa_pos=120, ref_aa=P105[119], alt_aa="W", call="mixed",
                confidence="candidate", n_alt=1, n_ref=2, flags=None),
    ]  # fmt: skip


def main_hits() -> list[dict[str, object]]:
    """Covering and non-covering placements for the main fixture."""
    origin = plus_hit(PEP_ORIGIN, 9901, 31)  # unwrapped 9991-10020: crosses the origin
    origin["crosses_origin"] = True
    dup = plus_hit(PEP1, 1001, 210)
    dup["uniqueness"] = "multi_identical"
    return [
        plus_hit(PEP1, 1001, 210),
        dup,  # the same peptide twice: one distinct peptide
        plus_hit(PEP2, 1001, 212),
        plus_hit(PEP_LONG, 1001, 200),
        # not covering VAR0001: upstream, ends just before the codon, other strand, other contig
        plus_hit(P105[99:110], 1001, 100),
        plus_hit(P105[199:213], 1001, 200),
        {**plus_hit("LVAEQTRGLR", 1001, 210), "strand": "-"},
        plus_hit("LVAEKTRGLM", 1001, 210, contig="chr2"),
        minus_hit(PEP_MINUS, 5600, 45),
        {**minus_hit("LVAEKTRGLW", 5600, 45), "strand": "+"},  # wrong strand for VAR0002
        plus_hit(PEP_FS, 1001, 245),
        origin,
    ]


def covering() -> dict[str, set[str]]:
    """Oracle: the distinct covering peptides of each ALT site of the main fixture."""
    return {
        "VAR0001": {PEP1, PEP2, PEP_LONG},
        "VAR0002": {PEP_MINUS},
        "VAR0003": {PEP_FS},
        "VAR0004": {PEP_ORIGIN},
    }
