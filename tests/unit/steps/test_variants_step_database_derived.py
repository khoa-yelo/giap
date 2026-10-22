"""M7-05 review B1 (PR #175): ``database_derived`` through the wired step, and N1.

docs/specs/M7-05-variant-flags.md rule 4: ``database_derived`` when at least one ALT-supporting
peptide is ``source = database`` with ``search_db = other`` (under ``genome_relation = related``);
degenerate row "frameshift event": rules 3-4 apply to frameshifts too; rule 5 and the pinned
``start_codon`` reading ("``end - 3`` on -, wrapped on circular replicons").

The giapsim step fixture is all ``database``/``same``, so these tests drive
``giap.steps.variants.run_variants`` on a hand-built genome instead:

* linear ``c1`` (400 nt), CDS ``G1`` + [30, 330) built from fixed codons (``PROT``);
* four ALT sites, in ``sites`` order: substitution A -> S at codon 39 (database/other),
  substitution at codon 96 (denovo), a +1 frameshift whose last prefix codon is 168
  (database/other), substitution at codon 210 (database/same). The supports alternate, so
  reversing or emptying the per-row supporting peptides, or dropping ``search_db``, changes at
  least one row's flags;
* N1: circular ``c2`` (600 nt) with a minus-strand CDS ``G2`` [450, 606) crossing the origin and
  an annotated GTG start; its first codon in transcription direction is [603, 606), key 3 after
  wrapping. An N-terminal peptide reading M there gets ``start_codon``.
"""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest

from giap.config import GiapConfig
from giap.mapping.sixframe import translate11_text
from giap.model import Replicon
from giap.model.coords import Interval
from giap.model.genes import Gene, GeneModel
from giap.steps.variants import run_variants

CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
    "K": "AAA", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT",
    "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}  # fmt: skip
_COMP = str.maketrans("ACGT", "TGCA")

# ---- c1: linear, G1 + [30, 330)
L1, G1_START = 400, 30
PROT = "M" + ("KDAEGTVNQRSWYFHPCK" * 6)[:99]
SEQ1 = "A" * G1_START + "".join(CODON[a] for a in PROT) + "A" * (L1 - G1_START - 3 * len(PROT))

# ---- c2: circular, G2 - [450, 606): gene-strand nt start with GTG; [600, 606) wraps to [0, 6)
L2, G2_START, G2_END = 600, 450, 606
PROT2 = "V" + ("KDAEGTSNQRHWYFTPCE" * 3)[:51]  # GTG start, so the genome reads V there
_G2_NT = "GTG" + "".join(CODON[a] for a in PROT2[1:])
_G2_PLUS = _G2_NT.translate(_COMP)[::-1]  # + strand bases of [450, 606)
SEQ2 = _G2_PLUS[150:] + "C" * (G2_START - 6) + _G2_PLUS[:150]


def _codon(k: int) -> int:
    """Genome coordinate of residue ``k`` of G1."""
    return G1_START + 3 * k


def _sub(lo: int, hi: int, k: int, alt: str) -> tuple[str, int, int]:
    """Peptide over G1 residues [lo, hi) with residue ``k`` replaced by ``alt``, and its locus."""
    pep = PROT[lo:k] + alt + PROT[k + 1 : hi]
    return pep, _codon(lo), _codon(hi)


# the four ALT sites of c1, in sites order
S1 = _sub(1, 11, 3, "S")  # A -> S at codon 39
S2 = _sub(18, 28, 22, "R")  # codon 96
S3 = _sub(55, 65, 60, "W")  # codon 210
FS_PRE = (_codon(40), _codon(47))  # [150, 171): G1 residues 40-46
FS_SUF = (FS_PRE[1] + 1, FS_PRE[1] + 19)  # [172, 190): one nt skipped (+1)
FS_PEP = PROT[40:47] + translate11_text(SEQ1[FS_SUF[0] : FS_SUF[1]])
FS_CODON = FS_PRE[1] - 3  # last prefix codon, 168

# N1: the N-terminal 8 residues of G2 lie on [582, 606), across the origin of c2
N1_PEP = "M" + PROT2[1:8]


def _hit(pep: str, contig: str, start: int, end: int, strand: str, source: str,
         dt: str) -> dict[str, Any]:  # fmt: skip
    return {
        "peptide": pep, "peptide_len": len(pep), "source": source, "contig": contig,
        "start": start, "end": end, "strand": strand, "n_diff": 1, "diff_types": dt,
        "placement_q": 0.0, "uniqueness": "unique",
    }  # fmt: skip


def _gm() -> GeneModel:
    reps = {"c1": Replicon("c1", L1, "linear"), "c2": Replicon("c2", L2, "circular")}
    g1 = Gene("G1", "cds", Interval("c1", G1_START, G1_START + 3 * len(PROT), "+"), PROT, None,
              "ATG")  # fmt: skip
    g2 = Gene("G2", "cds", Interval("c2", G2_START, G2_END, "-"), "M" + PROT2[1:], None, "GTG")
    return GeneModel(reps, (g1, g2), ())


def _cfg() -> GiapConfig:
    return GiapConfig.model_validate({"genome_relation": "related", "sample_name": "S1"})


SUPPORT: dict[str, tuple[str, str]] = {  # peptide -> (source, search_db)
    S1[0]: ("database", "other"),
    S2[0]: ("denovo", "same"),
    FS_PEP: ("database", "other"),
    S3[0]: ("database", "same"),
    N1_PEP: ("denovo", "same"),
}


def _run(
    extra: tuple[tuple[tuple[str, int, int], str, str], ...] = (), drop: tuple[str, ...] = ()
) -> pd.DataFrame:
    """``variants.tsv`` frame of ``run_variants`` on both contigs. ``extra`` adds G1 peptides
    as ``((peptide, start, end), source, search_db)``; ``drop`` removes peptides."""
    hits = [
        _hit(S1[0], "c1", S1[1], S1[2], "+", SUPPORT[S1[0]][0], "substitution"),
        _hit(S2[0], "c1", S2[1], S2[2], "+", SUPPORT[S2[0]][0], "substitution"),
        _hit(FS_PEP, "c1", FS_PRE[0], FS_SUF[1], "+", SUPPORT[FS_PEP][0], "frameshift"),
        _hit(S3[0], "c1", S3[1], S3[2], "+", SUPPORT[S3[0]][0], "substitution"),
        _hit(N1_PEP, "c2", G2_END - 24, G2_END, "-", SUPPORT[N1_PEP][0], "substitution"),
    ]
    hits = [h for h in hits if h["peptide"] not in drop]
    sdb = {p: s for p, (_, s) in SUPPORT.items() if p not in drop}
    for (pep, start, end), source, search_db in extra:
        hits.append(_hit(pep, "c1", start, end, "+", source, "substitution"))
        sdb[pep] = search_db
    tol = pd.DataFrame([{
        "peptide": FS_PEP, "search": "target", "contig": "c1", "start": FS_PRE[0],
        "end": FS_SUF[1], "strand": "+", "blocks": f"{FS_PRE[0]}-{FS_PRE[1]};"
        f"{FS_SUF[0]}-{FS_SUF[1]}",
    }])  # fmt: skip
    peptides = pd.DataFrame({"peptide": list(sdb), "search_db": list(sdb.values())})
    genome = {"c1": SEQ1, "c2": SEQ2}
    _, table, _ = run_variants(pd.DataFrame(hits), tol, peptides, _gm(), genome, _cfg())
    return table


def _flags(cell: object) -> list[str]:
    if isinstance(cell, str):
        return [f for f in cell.split(";") if f]
    return [str(f) for f in cell] if isinstance(cell, (tuple, list)) else []


@pytest.fixture(scope="module")
def table() -> pd.DataFrame:
    return _run()


def _site(table: pd.DataFrame, contig: str, codon: int, alt_type: str) -> dict[str, Any]:
    sel = table[(table["contig"] == contig) & (table["codon_start"] == codon)
                & (table["alt_type"] == alt_type)]  # fmt: skip
    assert len(sel) == 1, f"expected one {alt_type} site at {contig}:{codon}, got {len(sel)}"
    return sel.iloc[0].to_dict()


def test_fixture_has_the_four_alt_sites_in_alternating_support_order(table: pd.DataFrame) -> None:
    c1 = table[table["contig"] == "c1"]
    got = list(zip(c1["codon_start"], c1["alt_type"], c1["alt_aa"], strict=True))
    assert got == [(39, "substitution", "S"), (96, "substitution", "R"),
                   (FS_CODON, "frameshift", "fs+1"), (210, "substitution", "W")]  # fmt: skip
    assert (table["call"] == "ALT").all()


def test_substitution_from_a_database_other_peptide_is_database_derived(
    table: pd.DataFrame,
) -> None:
    row = _site(table, "c1", 39, "substitution")
    assert row["ref_aa"] == "A" and row["alt_aa"] == "S"
    assert _flags(row["flags"]) == ["mass_equivalent", "snv", "strain_difference",
                                    "database_derived"]  # fmt: skip


def test_substitution_from_a_denovo_peptide_is_not_database_derived(table: pd.DataFrame) -> None:
    flags = _flags(_site(table, "c1", 96, "substitution")["flags"])
    assert "strain_difference" in flags
    assert "database_derived" not in flags


def test_substitution_from_a_database_same_peptide_is_not_database_derived(
    table: pd.DataFrame,
) -> None:
    flags = _flags(_site(table, "c1", 210, "substitution")["flags"])
    assert "strain_difference" in flags
    assert "database_derived" not in flags


def test_frameshift_from_a_database_other_peptide_is_database_derived(
    table: pd.DataFrame,
) -> None:
    row = _site(table, "c1", FS_CODON, "frameshift")
    assert _flags(row["flags"]) == ["strain_difference", "database_derived"]


def test_each_row_carries_the_flags_of_its_own_supporting_peptides(table: pd.DataFrame) -> None:
    c1 = table[table["contig"] == "c1"]
    got = ["database_derived" in _flags(f) for f in c1["flags"]]
    assert got == [True, False, True, False]


def test_mixed_sources_at_one_alt_are_database_derived() -> None:
    # degenerate row "mixed sources at one ALT": a shorter de novo peptide with S1's A -> S
    shorter = _sub(2, 11, 3, "S")
    t = _run(extra=((shorter, "denovo", "same"),))
    row = _site(t, "c1", 39, "substitution")
    assert row["n_alt"] == 2
    assert "database_derived" in _flags(row["flags"])
    # the de novo peptide alone (S1 absent from the support) gives no database_derived
    t2 = _run(extra=((shorter, "denovo", "same"),), drop=(S1[0],))
    row2 = _site(t2, "c1", 39, "substitution")
    assert row2["n_alt"] == 1
    assert "database_derived" not in _flags(row2["flags"])


def test_n1_minus_strand_gtg_start_across_the_origin_gets_start_codon(
    table: pd.DataFrame,
) -> None:
    row = _site(table, "c2", 3, "substitution")  # [603, 606) wrapped
    assert (row["strand"], row["ref_aa"], row["alt_aa"]) == ("-", "V", "M")
    assert row["ref_codon"] == "GTG"
    flags = _flags(row["flags"])
    assert "start_codon" in flags and "snv" in flags
    assert "database_derived" not in flags
