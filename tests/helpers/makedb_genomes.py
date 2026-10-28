"""Hand-built genomes and Bakta-like JSON for the `gippa makedb` tests (issue M3-23).

Background is ``WALL`` (a stop in every frame of both strands in every 22 nt), so background
ORFs stay < 20 aa. Planted stretches are built from 15-nt units whose frame 0 has no stop while
the other five frames each have a stop *inside* every unit (``UNITS``), so a planted stretch
only opens the intended frame. ``UNITS2`` also leave the minus frame aligned with frame 0 open
(for an antisense ORF inside a plus-strand CDS).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from makedb_oracle import revcomp, translate  # tests/helpers is put on sys.path by the tests

WALL = "TAAATAAATAATTATTTATTTA"
UNITS = (
    "AACTATTTAACTAGG",
    "AATTATCTAGTTAAC",
    "ACTAGCTATCTAACT",
    "ACTAGCTTAATCTAC",
    "ACTAGCTTAATCTAT",
    "ATTAACTATTTAAAC",
    "CCTAGCCTAAATCAC",
    "CTAACTCATGCTAAG",
)
UNITS2 = (
    "AATTACGTGATTAAC",
    "ACTAATCATATAGCT",
    "ACTGTAACTAATTAC",
    "AGAGTAGCTAATCAT",
)


def stretch(n_codons: int, k: int = 0, units: tuple[str, ...] = UNITS) -> str:
    """``n_codons`` stop-free codons (frame 0) built from ``units`` starting at unit ``k``."""
    out = ""
    i = k
    while len(out) < 3 * n_codons:
        out += units[i % len(units)]
        i += 1
    return out[: 3 * n_codons]


@dataclass
class Contig:
    id: str
    length: int
    topology: str
    kind: str = "chromosome"
    seq: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if not self.seq:
            self.seq = list((WALL * (self.length // len(WALL) + 1))[: self.length])

    def put(self, start: int, s: str) -> None:
        for k, c in enumerate(s):
            if self.topology == "linear" and start + k >= self.length:
                raise ValueError("write beyond linear contig end")
            self.seq[(start + k) % self.length] = c

    def put_minus(self, start: int, s_minus: str) -> None:
        """Write ``s_minus`` (read 5'->3' on the minus strand) onto ``[start, start+len)``."""
        self.put(start, revcomp(s_minus))

    @property
    def text(self) -> str:
        return "".join(self.seq)


def feature(
    kind: str,
    contig: str,
    start0: int,
    end0: int,
    strand: str,
    *,
    locus: str | None = None,
    aa: str | None = None,
    product: str | None = None,
    length: int | None = None,
    **extra: Any,
) -> dict[str, Any]:
    """Bakta JSON feature (1-based closed; origin-crossing written as ``start > stop``)."""
    stop = end0 if length is None or end0 <= length else end0 - length
    f: dict[str, Any] = {
        "type": kind,
        "sequence": contig,
        "start": start0 + 1,
        "stop": stop,
        "strand": strand,
    }
    if locus is not None:
        f["locus"] = locus
        f["id"] = locus
    if aa is not None:
        f["aa"] = aa
    if product is not None:
        f["product"] = product
    f.update(extra)
    return f


def write_bakta_json(path: Path, contigs: list[Contig], features: list[dict[str, Any]]) -> Path:
    data = {
        "version": {"bakta": "1.12.2", "db": {"version": "6.0", "type": "full"}},
        "sequences": [
            {
                "id": c.id,
                "description": "",
                "nt": c.text,
                "length": c.length,
                "complete": True,
                "type": c.kind,
                "topology": c.topology,
            }
            for c in contigs
        ],
        "features": features,
    }
    path.write_text(json.dumps(data, indent=1) + "\n")
    return path


def write_fasta(path: Path, contigs: list[Contig]) -> Path:
    lines = []
    for c in contigs:
        lines.append(f">{c.id}")
        t = c.text
        lines.extend(t[i : i + 80] for i in range(0, len(t), 80))
    path.write_text("\n".join(lines) + "\n")
    return path


# ---------------------------------------------------------------- the example genome

L_CHR = 3300  # 150 x 22: the WALL background is periodic across the origin
L_LIN = 660
G1, G2, G3, G4, G5 = "GSIM_00005", "GSIM_00010", "GSIM_00015", "GSIM_00020", "GSIM_00025"
S1, P1 = "GSIM_00030", "GSIM_00035"
G1_PRODUCT = "DNA gyrase subunit A"


@dataclass
class ExampleGenome:
    contigs: list[Contig]
    features: list[dict[str, Any]]
    aa: dict[str, str]  # locus tag -> annotated protein (no terminal '*')
    expect: dict[str, str]  # name -> expected header or sequence (see build_example)

    @property
    def genome(self) -> dict[str, str]:
        return {c.id: c.text for c in self.contigs}


def _gene_nt(n_body: int, k: int, start: str = "ATG") -> str:
    return start + stretch(n_body, k) + "TAA"


def build_example() -> ExampleGenome:
    """chr1 (circular, 3,300 nt) and lin1 (linear, 660 nt) with the spec's examples 4-7, 9.

    0-based layout on chr1 (plus coordinates):

    * G2 ``+ [240, 480)``; inside it the minus-strand ORF ``[300, 390)`` (30 aa, antisense, ex. 5)
      bounded by ``TTA`` (minus TAA) at 297 and 390.
    * G1 ``+ [600, 900)``, product ``DNA gyrase subunit A``; in-frame stop at 561 and a 12-codon
      stop-free upstream part ``[564, 600)`` whose first codon is GTG (ex. 7; patch start).
    * G4 ``+ [1000, 1150)``; its JSON ``aa`` carries a terminal ``*``.
    * intergenic ORF ``+ [1200, 1290)`` (30 aa, first codon ATG; ex. 4 and patch ex. 9).
    * 19-aa ORF ``+ [1500, 1557)`` (ex. 6) and 20-aa ORF ``+ [1700, 1760)`` (boundary).
    * sORF S1 ``+ [1900, 1960)`` (19 aa) inside the 20-aa stop-to-stop frame ``[1897, 1957)``.
    * G3 ``- [2000, 2300)`` with a 10-codon upstream part ``[2300, 2330)`` on the minus strand.
    * tRNA T1 ``+ [2500, 2576)`` under the 30-aa ORF ``+ [2490, 2580)`` (rna_overlap: omitted).
    * pseudogene P1 ``+ [2700, 2900)`` with a 30-aa ORF ``+ [2731, 2821)`` in another frame.
    * G5 ``- [3000, 3150)``: same protein as G4 (duplicate sequences are all written).
    * origin-crossing intergenic ORF ``+ [3240, 3390)`` (50 aa; longer than the matcher's wrap).

    lin1 (linear): 19-aa ORF truncated at the contig start ``[0, 57)``; 25-aa ORF truncated at
    the contig end ``[585, 660)``; a 60-codon stretch ``[203, 383)`` with one ``N`` codon.
    """
    chr1 = Contig("chr1", L_CHR, "circular", "chromosome")
    lin1 = Contig("lin1", L_LIN, "linear", "contig")
    feats: list[dict[str, Any]] = []
    aa: dict[str, str] = {}
    expect: dict[str, str] = {}

    def protein(nt: str) -> str:
        return "M" + translate(nt[3:-3])

    # G2 with the antisense ORF (example 5)
    g2 = "ATG" + stretch(18, 1) + "TTA" + stretch(30, 0, UNITS2) + "TTA" + stretch(28, 2) + "TAA"
    assert len(g2) == 240
    chr1.put(234, "TAAGCA")
    chr1.put(240, g2)
    aa[G2] = protein(g2)
    feats.append(feature("cds", "chr1", 240, 480, "+", locus=G2, aa=aa[G2]))  # no product
    expect["antisense"] = "GIPPAORF_antisense_chr1_301_390_m"
    expect["antisense_aa"] = translate(revcomp(chr1.text[300:390]))

    # G1 with a 12-aa in-frame upstream extension (example 7)
    up = "GTG" + stretch(11, 3)
    g1 = _gene_nt(98, 4)
    chr1.put(561, "TAA" + up + g1)
    aa[G1] = protein(g1)
    feats.append(feature("cds", "chr1", 600, 900, "+", locus=G1, aa=aa[G1], product=G1_PRODUCT))
    expect["nterm_plus_up"] = translate(up)
    expect["nterm_plus_aa"] = translate(up) + aa[G1][:30]
    expect["patch_nterm_aa"] = "M" + translate(chr1.text[567:897])

    # G4 (+) and G5 (-) encode the same protein; G4's JSON aa ends with '*'
    g4 = _gene_nt(48, 5)
    chr1.put(994, "TAAGCA" + g4)
    aa[G4] = protein(g4)
    feats.append(feature("cds", "chr1", 1000, 1150, "+", locus=G4, aa=aa[G4] + "*", product="P"))
    chr1.put_minus(3000, "TAAGCT" + g4)  # minus string: stop, 1 codon, gene
    # the minus string occupies [3000, 3156): gene [3000, 3150), upstream [3150, 3156)
    aa[G5] = aa[G4]
    feats.append(feature("cds", "chr1", 3000, 3150, "-", locus=G5, aa=aa[G5], product="P"))

    # intergenic ORFs (examples 4 and 6, the 20-aa boundary)
    ex4 = "ATG" + stretch(29, 6)
    chr1.put(1197, "TAA" + ex4 + "TAA")
    expect["intergenic"] = "GIPPAORF_intergenic_chr1_1201_1290_p"
    expect["intergenic_aa"] = translate(ex4)
    chr1.put(1497, "TAA" + stretch(19, 7) + "TAA")
    expect["orf19"] = "GIPPAORF_intergenic_chr1_1501_1557_p"
    chr1.put(1697, "TAA" + stretch(20, 1) + "TAA")
    expect["orf20"] = "GIPPAORF_intergenic_chr1_1701_1760_p"

    # sORF S1 in frame with a 20-aa stop-to-stop segment
    s1 = _gene_nt(18, 2)
    chr1.put(1894, "TAAGCA" + s1)
    aa[S1] = protein(s1)
    feats.append(feature("sorf", "chr1", 1900, 1960, "+", locus=S1, aa=aa[S1], product="sORF"))
    expect["sorf_frame"] = "_chr1_1898_1957_p"

    # G3 on the minus strand with a 10-aa upstream extension
    up3 = stretch(10, 5)
    g3 = _gene_nt(98, 6)
    chr1.put_minus(2000, "TAA" + up3 + g3)  # occupies [2000, 2333)
    aa[G3] = protein(g3)
    feats.append(feature("cds", "chr1", 2000, 2300, "-", locus=G3, aa=aa[G3], product="Q"))
    expect["nterm_minus_aa"] = translate(up3) + aa[G3][:30]

    # tRNA under an ORF -> rna_overlap (omitted)
    chr1.put(2487, "TAA" + stretch(30, 3) + "TAA")
    feats.append(
        feature("tRNA", "chr1", 2500, 2576, "+", id="T1", product="tRNA-Ala", amino_acid="Ala")
    )
    expect["rna"] = "_chr1_2491_2580_p"

    # pseudogene with an ORF in another frame on the same strand -> pseudogene (omitted)
    chr1.put(2700, "ATG")
    chr1.put(2728, "TAA" + stretch(30, 4) + "TAA")
    feats.append(feature("cds", "chr1", 2700, 2900, "+", locus=P1, aa="M" + "A" * 65, pseudo=True))
    expect["pseudo"] = "_chr1_2732_2821_p"

    # origin-crossing ORF: [3240, 3390) = 50 aa, stops at 3237 and 90 (unwrapped 3390)
    orig = stretch(50, 2)
    chr1.put(3237, "TAA" + orig + "TAA")
    expect["origin"] = "GIPPAORF_intergenic_chr1_3241_3390_p"
    expect["origin_aa"] = translate(orig)

    # lin1 (linear)
    lin1.put(0, stretch(19, 0) + "TAA")
    expect["lin_start19"] = "_lin1_1_57_p"
    with_n = stretch(60, 1)
    with_n = with_n[:90] + "GNA" + with_n[93:]
    lin1.put(200, "TAA" + with_n + "TAA")
    lin1.put(582, "TAA" + stretch(25, 3))
    expect["lin_end"] = "GIPPAORF_intergenic_lin1_586_660_p"
    expect["lin_end_aa"] = translate(stretch(25, 3))

    # JSON in non-coordinate order (gene-model order is by coordinate)
    feats.sort(key=lambda f: str(f.get("locus", f.get("id"))), reverse=True)
    return ExampleGenome([chr1, lin1], feats, aa, expect)


def write_patch_gff3(path: Path, genome: ExampleGenome) -> Path:
    """Patch with a novel ORF (example 9) and a G1 N-terminal extension (new GTG start at 564)."""
    lines = [
        "##gff-version 3",
        f"##sequence-region chr1 1 {L_CHR}",
        f"##sequence-region lin1 1 {L_LIN}",
        f"chr1\tgippa\tCDS\t565\t900\t.\t+\t0\tID={G1};locus_tag={G1};Note=nterm_extension",
        "chr1\tgippa\tCDS\t1201\t1293\t.\t+\t0\tID=NOVEL_chr1_1201_p;product=novel ORF",
    ]
    path.write_text("\n".join(lines) + "\n")
    return path
