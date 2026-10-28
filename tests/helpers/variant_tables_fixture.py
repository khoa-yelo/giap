"""Fixture genome, builders and oracles for the M7-03 tests (docs/specs/M7-03-*.md).

Written by the test author from the spec text. The oracle (:func:`minimal_codon`) does not import
``gippa``; inputs reach ``vcf_records`` and the writers only through the public M7-02 functions
``pileup`` and ``sites`` (whatever ``pileup`` returns is passed through untouched), so these tests
do not depend on how the pileup carries its per-codon ALT detail (M7-03 lead decision 2).

Genome (forward strand):

* ``chr``, linear, 2,400 nt.
  - ``NISS_00012`` ``+`` [868, 1168): codon ``i`` at ``868 + 3i``; i=44 (1000) ``GCT`` A (aa_pos
    45, Appendix B.7), i=50 (1018) ``AAT`` N, i=60/61 (1048/1051) ``TCT GAT`` SD, i=70/71
    (1078/1081) ``GGT ACT`` GT, i=80 (1108) ``CTT`` L.
  - ``NISS_00020`` ``-`` [1300, 1600): codon ``i`` at ``1600 - 3(i+1)``; i=10 (1567) gene strand
    ``GCT`` (forward ``AGC``), i=20/21 (1537/1534) gene strand ``TCT GAT`` SD.
  - no gene at [180, 260).
* ``c2``, circular, 1,000 nt, no genes; codon 996 ``TCT`` and the origin-crossing codon 999
  ``GAT`` (bases 999, 0, 1).
"""

from __future__ import annotations

import importlib
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from types import ModuleType
from typing import Any

import pandas as pd

from gippa.config.model import GippaConfig
from gippa.model import Replicon
from gippa.model.coords import Interval
from gippa.model.genes import Gene, GeneModel
from helpers.align_oracle import place
from helpers.fs_genomes import background, codons_for
from helpers.pileup_oracle import make_hit

# --------------------------------------------------------------------------- table 11 oracle

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE: dict[str, str] = {
    a + b + c: _AAS[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
CODONS64: tuple[str, ...] = tuple(sorted(CODE))
RESIDUES21: tuple[str, ...] = tuple(sorted(set(CODE.values())))  # 20 + '*'
_RC = str.maketrans("ACGT", "TGCA")


def rc(s: str) -> str:
    return s.translate(_RC)[::-1]


def translate(nt: str) -> str:
    return "".join(CODE[nt[i : i + 3]] for i in range(0, len(nt) - 2, 3))


def encodes(codon: str, aa: str) -> bool:
    return CODE[codon] in ("IL" if aa == "J" else aa)


def minimal_codon(ref: str, aa: str) -> tuple[str, int]:
    """Brute force G14 + lead decision 1: fewest differences, then lexicographic."""
    best = min((sum(x != y for x, y in zip(c, ref, strict=True)), c) for c in CODONS64
               if encodes(c, aa))  # fmt: skip
    return best[1], best[0]


# --------------------------------------------------------------------------- genome

L_CHR, L_C2 = 2400, 1000
SAMPLE = "pilot_A3"
SHA = "0123456789abcdef" * 4
VERSION = "0.0.0-test"
FILLER = "AGSTEKRDQPFYHWMV"


@dataclass(frozen=True)
class GSpec:
    tag: str
    start: int
    end: int
    strand: str

    def codon(self, i: int) -> int:
        return self.start + 3 * i if self.strand == "+" else self.end - 3 * (i + 1)

    @property
    def n(self) -> int:
        return (self.end - self.start) // 3


GP = GSpec("NISS_00012", 868, 1168, "+")
GM_ = GSpec("NISS_00020", 1300, 1600, "-")

_GP_CODONS = {44: "GCT", 50: "AAT", 60: "TCT", 61: "GAT", 70: "GGT", 71: "ACT", 80: "CTT"}
_GM_CODONS = {10: "GCT", 20: "TCT", 21: "GAT"}


def _gene_nt(g: GSpec, fixed: Mapping[int, str], seed: int) -> str:
    prot = "".join(FILLER[(seed + 7 * i) % len(FILLER)] for i in range(g.n))
    cod = codons_for("M" + prot[1:])
    for i, c in fixed.items():
        cod[i] = c
    return "".join(cod)


def _build() -> tuple[str, str]:
    seq = background(L_CHR, seed=303)
    seq = place(seq, GP.start, _gene_nt(GP, _GP_CODONS, 1), "+")
    seq = place(seq, GM_.start, _gene_nt(GM_, _GM_CODONS, 4), "-")
    return seq, place(background(L_C2, seed=304), 996, "TCTGAT", "+")


SEQ_CHR, SEQ_C2 = _build()
GENOME: dict[str, str] = {"chr": SEQ_CHR, "c2": SEQ_C2}
REPS: dict[str, Replicon] = {
    "chr": Replicon("chr", L_CHR, "linear"),
    "c2": Replicon("c2", L_C2, "circular"),
}


def gene_codon_nt(g: GSpec, i: int) -> str:
    c = g.codon(i)
    nt = SEQ_CHR[c : c + 3]
    return rc(nt) if g.strand == "-" else nt


def _gene(g: GSpec) -> Gene:
    nt = "".join(gene_codon_nt(g, i) for i in range(g.n))
    aa = "M" + translate(nt)[1:]
    iv = Interval("chr", g.start, g.end, g.strand)  # type: ignore[arg-type]
    return Gene(g.tag, "cds", iv, aa, None, "ATG")


GMODEL = GeneModel(REPS, (_gene(GP), _gene(GM_)), ())
CFG = GippaConfig.model_validate({})

assert gene_codon_nt(GP, 44) == "GCT" and GP.codon(44) == 1000
assert SEQ_CHR[1567:1570] == "AGC" and gene_codon_nt(GM_, 10) == "GCT"
assert GM_.codon(20) == 1537 and GM_.codon(21) == 1534
assert SEQ_C2[996:999] == "TCT" and SEQ_C2[999] + SEQ_C2[:2] == "GAT"

# --------------------------------------------------------------------------- modules


def pl() -> ModuleType:
    return importlib.import_module("gippa.variants.pileup")


def al() -> ModuleType:
    return importlib.import_module("gippa.variants.align")


def vr() -> ModuleType:
    return importlib.import_module("gippa.variants.vcf_records")


def vt() -> ModuleType:
    return importlib.import_module("gippa.io.variant_tables")


def vcf_io() -> ModuleType:
    return importlib.import_module("gippa.io.vcf")


# --------------------------------------------------------------------------- hit builders

Op = tuple[str, tuple[int, ...], str, str]


def aa_at(contig: str, codon: int, strand: str) -> str:
    seq, n = GENOME[contig], REPS[contig].length
    nt = (seq + seq)[codon % n : codon % n + 3]
    return translate(rc(nt) if strand == "-" else nt)


def match_op(contig: str, codon: int, strand: str) -> Op:
    r = aa_at(contig, codon, strand)
    return ("match", (codon,), r, r)


def gene_peptide(g: GSpec, first: int, last: int, at: Mapping[int, Op]) -> list[Op]:
    """Ops over gene residues ``[first, last]``; residue ``i`` in ``at`` uses that op (which may
    consume the next residue's codon too: a two-codon op at ``i`` skips ``i + 1``)."""
    ops: list[Op] = []
    i = first
    while i <= last:
        if i in at:
            op = at[i]
            ops.append(op)
            i += len(op[1])
        else:
            ops.append(match_op("chr", g.codon(i), g.strand))
            i += 1
    return ops


def sub(g: GSpec, i: int, alt: str) -> Op:
    c = g.codon(i)
    return ("substitution", (c,), aa_at("chr", c, g.strand), alt)


def two(g: GSpec, i: int, kind: str, alt: str) -> Op:
    c1, c2 = g.codon(i), g.codon(i + 1)
    ref = aa_at("chr", c1, g.strand) + aa_at("chr", c2, g.strand)
    return (kind, (c1, c2), ref, alt)


def peptides_at(g: GSpec, i: int, op: Op, n: int, shift: int = 0) -> list[Any]:
    """``n`` distinct peptides carrying ``op`` at residue ``i`` (different left flanks)."""
    span = len(op[1])
    return [
        make_hit(al(), "chr", g.strand, gene_peptide(g, i - 1 - k - shift, i + span, {i: op}),
                 length=L_CHR)
        for k in range(n)
    ]  # fmt: skip


def ref_peptides(g: GSpec, i: int, n: int) -> list[Any]:
    """``n`` distinct REF peptides covering residue ``i`` (different right flanks)."""
    return [
        make_hit(al(), "chr", g.strand, gene_peptide(g, i - 1, i + 2 + k, {}), length=L_CHR)
        for k in range(n)
    ]


def raw_hit(contig: str, strand: str, ops: Sequence[Op], frameshift: Any = None) -> Any:
    return make_hit(al(), contig, strand, ops, frameshift=frameshift,
                    length=REPS[contig].length)  # fmt: skip


def support_of(aligns: Iterable[Any], source: str = "denovo") -> dict[str, Any]:
    m = pl()
    return {a.peptide: m.Support(a.peptide, source, "same", 0.001, "unique") for a in aligns}


# --------------------------------------------------------------------------- pipeline


def run(
    aligns: Sequence[Any],
    gm: GeneModel = GMODEL,
    genome: Mapping[str, str] = GENOME,
) -> tuple[Any, pd.DataFrame]:
    """``(pile, sites)`` from the public M7-02 functions; ``pile`` is passed through as is."""
    m = pl()
    support = support_of(aligns)
    pile = m.pileup(list(aligns), gm, genome, CFG, support=support)
    return pile, m.sites(pile, CFG)


def records(
    aligns: Sequence[Any],
    extra: Mapping[str, Sequence[object]] | None = None,
    gm: GeneModel = GMODEL,
    genome: Mapping[str, str] = GENOME,
    reps: Mapping[str, Replicon] = REPS,
) -> tuple[list[Any], pd.DataFrame, Any, pd.DataFrame]:
    """``(vcf records, variants table, pile, sites)``; ``extra`` adds later columns to sites."""
    pile, sites = run(aligns, gm, genome)
    for col, vals in (extra or {}).items():
        sites[col] = list(vals)
    recs, table = vr().vcf_records(sites, pile, genome, reps)
    return sorted(recs, key=lambda r: (r.contig, r.pos0, len(r.ref))), table, pile, sites


def write(path: Any, recs: Sequence[Any], reps: Mapping[str, Replicon] = REPS,
          genome: Mapping[str, str] = GENOME, sample: str = SAMPLE) -> None:  # fmt: skip
    v = vr()
    vcf_io().write_vcf(
        path, recs, list(reps.values()), gippa_version=VERSION, genome_sha256=SHA, sample=sample,
        info_fields=v.INFO_FIELDS, format_fields=v.FORMAT_FIELDS, sequences=genome,
    )  # fmt: skip


def body_lines(path: Any) -> list[list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    return [ln.split("\t") for ln in lines if not ln.startswith("#")]


def info_dict(text: str) -> dict[str, str]:
    if text == ".":
        return {}
    out = {}
    for item in text.split(";"):
        k, _, v = item.partition("=")
        out[k] = v
    return out


def is_na(v: object) -> bool:
    if v is None or v is pd.NA:
        return True
    try:
        return bool(pd.isna(v))  # type: ignore[call-overload]
    except (TypeError, ValueError):
        return False


def flag_set(v: object) -> set[str]:
    """A ``flags`` cell as a set: tuple/list, ``;``-joined string, or empty/NA."""
    if isinstance(v, list | tuple | set | frozenset):
        return {str(x) for x in v}
    if is_na(v) or v == "":
        return set()
    return set(str(v).split(";"))


def site_row(table: pd.DataFrame, contig: str, codon: int, strand: str,
             alt_type: str | None = None) -> list[dict[str, Any]]:  # fmt: skip
    sel = table[(table["contig"] == contig) & (table["codon_start"] == codon)
                & (table["strand"] == strand)]  # fmt: skip
    if alt_type is not None:
        sel = sel[sel["alt_type"] == alt_type]
    return [dict(r) for _, r in sel.iterrows()]


VARIANTS_COLUMNS: tuple[str, ...] = (
    "variant_id", "contig", "codon_start", "strand", "gene", "aa_pos", "ref_codon", "ref_aa",
    "alt_aa", "alt_type", "call", "confidence", "n_ref", "n_alt", "nt_changes", "alt_codon_vcf",
    "class_error_rate", "class_error_n", "source", "flags",
)  # fmt: skip
PILEUP_COLUMNS: tuple[str, ...] = (
    "contig", "codon_start", "strand", "gene", "aa_pos", "ref_aa", "n_ref", "n_il",
    "alt_counts", "call",
)  # fmt: skip
VARIANT_FLAGS: tuple[str, ...] = (
    "mass_equivalent", "snv", "multi_nt", "strain_difference", "database_derived",
    "start_codon", "stop_codon", "vcf_unrepresentable",
)  # fmt: skip
