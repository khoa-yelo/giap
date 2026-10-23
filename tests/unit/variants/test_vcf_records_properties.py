"""M7-03 properties (docs/specs/M7-03-variant-tables-and-vcf.md, acceptance criterion 3).

For random genomes and random ALT sites (substitutions, ``J``, swaps, 2<->2 blocks; both
strands; several alleles per codon): every record's REF matches the genome (``write_vcf`` is
called with ``sequences=genome``), and applying each ALT to the genome and translating on the
gene's strand gives that allele's ``AA_ALT`` (``J`` = I or L) with ``NT_CHANGES`` nucleotide
differences that are minimal (G14). Inputs go through the public M7-02 ``pileup``/``sites``.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any

from helpers.variant_tables_fixture import (
    CODE,
    CODONS64,
    SAMPLE,
    SHA,
    VERSION,
    al,
    encodes,
    make_hit,
    pl,
    rc,
    support_of,
    translate,
    vcf_io,
    vr,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.config.model import GiapConfig
from giap.model import Replicon
from giap.model.genes import GeneModel

L = 120
CFG = GiapConfig.model_validate({})
SENSE = sorted(c for c in CODONS64 if CODE[c] != "*")
PEPTIDE_AAS = "ACDEFGHIKLMNPQRSTVWY"


def _aa(seq: str, codon: int, strand: str) -> str:
    nt = seq[codon : codon + 3]
    return translate(rc(nt) if strand == "-" else nt)


def _fold(peptide: str) -> str:
    return peptide.replace("I", "L")


@st.composite
def scenario(draw: st.DrawFn) -> dict[str, Any]:
    seq = "".join(draw(st.lists(st.sampled_from(SENSE), min_size=L // 3, max_size=L // 3)))
    strand = draw(st.sampled_from("+-"))
    c = 3 * draw(st.integers(4, L // 3 - 6))  # codon (lowest coordinate) of the site
    kind = draw(st.sampled_from(["substitution", "swap", "block"]))
    step = 3 if strand == "+" else -3
    codons = (c,) if kind == "substitution" else (c, c + step)
    ref = "".join(_aa(seq, x, strand) for x in codons)
    if kind == "swap" and ref[0] == ref[1]:
        kind, codons, ref = "substitution", (c,), ref[0]
    n_alleles = draw(st.integers(1, 3))
    alleles: list[tuple[str, int]] = []
    for _ in range(n_alleles):
        if kind == "swap":
            alt = ref[::-1]
        else:
            alt = "".join(draw(st.sampled_from(PEPTIDE_AAS)) for _ in ref)
            if any(a == r or (a in "IL" and r in "IL") for a, r in zip(alt, ref, strict=True)):
                continue
        # alleles equal after I/L folding are one allele (block AI / AL, substitution I / L)
        if _fold(alt) in {_fold(a) for a, _ in alleles}:
            continue
        alleles.append((alt, draw(st.integers(1, 3))))
    if not alleles:
        alleles = [("W" if ref[0] != "W" else "C", 1)] if kind == "substitution" else []
        if not alleles:
            kind, codons, ref = "substitution", (c,), ref[0]
            alleles = [("W" if ref != "W" else "C", 1)]
    return {"seq": seq, "strand": strand, "codons": codons, "kind": kind, "ref": ref,
            "alleles": alleles}  # fmt: skip


def _hits(sc: dict[str, Any]) -> list[Any]:
    seq, strand, codons = sc["seq"], sc["strand"], sc["codons"]
    step = 3 if strand == "+" else -3
    first = codons[0]
    out = []
    for alt, n in sc["alleles"]:
        for k in range(n):
            ops = []
            for j in range(k + 1, 0, -1):  # left flank, k + 1 codons
                x = first - j * step
                r = _aa(seq, x, strand)
                ops.append(("match", (x,), r, r))
            ops.append((sc["kind"], tuple(codons), sc["ref"], alt))
            nxt = codons[-1] + step
            r = _aa(seq, nxt, strand)
            ops.append(("match", (nxt,), r, r))
            out.append(make_hit(al(), "r", strand, ops, length=L))
    return out


@given(scenario())
@settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
def test_alt_translates_to_aa_alt_and_ref_matches_genome(sc: dict[str, Any]) -> None:
    rep = Replicon("r", L, "linear")
    genome = {"r": sc["seq"]}
    gm = GeneModel({"r": rep}, (), ())
    aligns = _hits(sc)
    m = pl()
    pile = m.pileup(aligns, gm, genome, CFG, support=support_of(aligns))
    sites = m.sites(pile, CFG)
    v = vr()
    recs, table = v.vcf_records(sites, pile, genome, {"r": rep})
    assert len(table) == len(sites)
    assert len(recs) == 1
    (r,) = recs
    with tempfile.TemporaryDirectory() as d:
        vcf_io().write_vcf(
            Path(d) / "v.vcf", recs, [rep], giap_version=VERSION, genome_sha256=SHA,
            sample=SAMPLE, info_fields=v.INFO_FIELDS, format_fields=v.FORMAT_FIELDS,
            sequences=genome,
        )  # fmt: skip
    width = 3 * len(sc["codons"])
    assert r.pos0 == min(sc["codons"]) and len(r.ref) == width
    seq = sc["seq"]
    assert seq[r.pos0 : r.pos0 + width] == r.ref
    aa_alt = r.info["AA_ALT"]
    nt_changes = r.info["NT_CHANGES"]
    aa_alt = [aa_alt] if isinstance(aa_alt, str) else list(aa_alt)
    nt_changes = [nt_changes] if isinstance(nt_changes, int) else list(nt_changes)
    assert len(aa_alt) == len(r.alts) == len(nt_changes) == len(sc["alleles"])
    strand = sc["strand"]
    gene_ref = rc(r.ref) if strand == "-" else r.ref
    assert translate(gene_ref) == r.info["AA_REF"]
    for alt_nt, want, n in zip(r.alts, aa_alt, nt_changes, strict=True):
        assert len(alt_nt) == width
        mutated = seq[: r.pos0] + alt_nt + seq[r.pos0 + width :]
        gene_alt = mutated[r.pos0 : r.pos0 + width]
        gene_alt = rc(gene_alt) if strand == "-" else gene_alt
        for i, aa in enumerate(want):
            assert encodes(gene_alt[3 * i : 3 * i + 3], aa), (gene_alt, want)
        assert n == sum(a != b for a, b in zip(gene_alt, gene_ref, strict=True))
        # G14: no codon set encoding the ALT residues is closer to REF (per codon)
        for i, aa in enumerate(want):
            rc3 = gene_ref[3 * i : 3 * i + 3]
            best = min(sum(x != y for x, y in zip(c, rc3, strict=True)) for c in CODONS64
                       if encodes(c, aa))  # fmt: skip
            got = sum(x != y for x, y in zip(gene_alt[3 * i : 3 * i + 3], rc3, strict=True))
            assert got == best
    n_alt = sorted((n for _, n in sc["alleles"]), reverse=True)
    assert list(r.ad[1:]) == n_alt and r.ad[0] == 0 and r.dp == sum(n_alt)
