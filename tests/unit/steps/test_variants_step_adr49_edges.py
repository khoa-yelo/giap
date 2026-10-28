"""ADR-0049 A6a through ``gippa.steps.variants.run_variants_qc`` at the edges the review of PR #269
names (items 2b and 2e; docs/adr/ADR-0049-m7-variant-precision.md).

* 2b: a substitution at the first codon of a ``+`` gene that starts at coordinate 0 of a linear
  contig, and at the first (transcription order) codon of a ``-`` gene that ends at the contig's
  end. There is no residue before either site: the step must not crash (the window would run off
  the contig) and must not flag the site ``mass_equivalent``.
* 2e: one codon with two ALT rows. ``G>Q`` comes from peptides that also read as the left 2->1
  block ``AG -> Q``; ``G>W`` from peptides that do not. Only the ``G>Q`` row is flagged: each row
  is read against the alignments of its own ALT peptides.

Hand-built linear genome, 0-based half-open hits, fixed codons per residue.
"""

from __future__ import annotations

from typing import Any

import pandas as pd

from gippa.config import GippaConfig
from gippa.model import Replicon
from gippa.model.coords import Interval
from gippa.model.genes import Gene, GeneModel
from gippa.steps.variants import run_variants_qc

ME = "mass_equivalent"
CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
    "I": "ATT", "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA",
    "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}  # fmt: skip
_COMP = str.maketrans("ACGT", "TGCA")

# e1: + gene from coordinate 0; e2: - gene ending at the contig end. Site: M (codon 0 of the
# protein) read as Q by two peptides.
EDGE_PROT = "MSVATEPIKRSCS"
EDGE_PEPS = ("QSVATEPIK", "QSVATEPIKR")
SEQ_E1 = "".join(CODON[a] for a in EDGE_PROT) + "A" * 30
E2_START = 30
SEQ_E2 = "A" * E2_START + "".join(CODON[a] for a in EDGE_PROT).translate(_COMP)[::-1]
E2_END = len(SEQ_E2)

# c3: + gene at [30, ...); site G (protein index 7) read as Q (block AG -> Q) and as W (no block)
C3_START = 30
C3_PROT = "MSCSDWAGSVATEPIKRSCS"
C3_SITE = 7
C3_END = C3_START + 3 * len(C3_PROT)
C3_PEPS = {"Q": ("QSVATEPIK", "QSVATEPIKR"), "W": ("WSVATEPIK", "WSVATEPIKR")}
SEQ_C3 = "A" * C3_START + "".join(CODON[a] for a in C3_PROT) + "A" * 30


def _hit(pep: str, contig: str, start: int, end: int, strand: str) -> dict[str, Any]:
    return {
        "peptide": pep, "peptide_len": len(pep), "source": "denovo", "contig": contig,
        "start": start, "end": end, "strand": strand, "n_diff": 1,
        "diff_types": "substitution", "placement_q": 0.0, "uniqueness": "unique",
    }  # fmt: skip


def _gm() -> GeneModel:
    n = 3 * len(EDGE_PROT)
    reps = {
        "e1": Replicon("e1", len(SEQ_E1), "linear"),
        "e2": Replicon("e2", len(SEQ_E2), "linear"),
        "c3": Replicon("c3", len(SEQ_C3), "linear"),
    }
    genes = (
        Gene("E1", "cds", Interval("e1", 0, n, "+"), EDGE_PROT, None, "ATG"),
        Gene("E2", "cds", Interval("e2", E2_START, E2_END, "-"), EDGE_PROT, None, "ATG"),
        Gene("G3", "cds", Interval("c3", C3_START, C3_END, "+"), C3_PROT, None, "ATG"),
    )
    return GeneModel(reps, genes, ())


def _hits() -> list[dict[str, Any]]:
    out = []
    for pep in EDGE_PEPS:
        out.append(_hit(pep, "e1", 0, 3 * len(pep), "+"))
        out.append(_hit(pep, "e2", E2_END - 3 * len(pep), E2_END, "-"))
    for peps in C3_PEPS.values():
        for pep in peps:
            lo = C3_START + 3 * C3_SITE
            out.append(_hit(pep, "c3", lo, lo + 3 * len(pep), "+"))
    return out


def _table() -> pd.DataFrame:
    hits = _hits()
    peptides = pd.DataFrame({"peptide": sorted({h["peptide"] for h in hits}),
                             "search_db": "same"})  # fmt: skip
    cfg = GippaConfig.model_validate({"genome_relation": "related", "sample_name": "S1"})
    genome = {"e1": SEQ_E1, "e2": SEQ_E2, "c3": SEQ_C3}
    _, table, _, _ = run_variants_qc(pd.DataFrame(hits), None, peptides, _gm(), genome, cfg)
    return table


def _flags(cell: object) -> list[str]:
    if isinstance(cell, str):
        return [f for f in cell.split(";") if f]
    return [str(f) for f in cell] if isinstance(cell, (tuple, list)) else []


def _rows(table: pd.DataFrame, contig: str, codon: int) -> list[dict[str, Any]]:
    sel = table[(table["contig"] == contig) & (table["codon_start"] == codon)
                & (table["alt_type"] == "substitution")]  # fmt: skip
    return [r for r in sel.to_dict("records")]


def test_fixture_peptides_read_one_substitution() -> None:
    for pep in EDGE_PEPS:
        assert pep[1:] == EDGE_PROT[1 : len(pep)]
    for alt, peps in C3_PEPS.items():
        for pep in peps:
            assert pep == alt + C3_PROT[C3_SITE + 1 : C3_SITE + len(pep)]
    assert C3_PROT[C3_SITE - 1 : C3_SITE + 1] == "AG"


def test_plus_site_at_the_first_codon_of_a_linear_contig_is_called_and_not_flagged() -> None:
    rows = _rows(_table(), "e1", 0)
    assert len(rows) == 1, rows
    r = rows[0]
    assert (r["ref_aa"], r["alt_aa"], r["strand"], r["call"]) == ("M", "Q", "+", "ALT"), r
    assert ME not in _flags(r["flags"]), r


def test_minus_site_at_the_last_codon_of_a_linear_contig_is_called_and_not_flagged() -> None:
    rows = _rows(_table(), "e2", E2_END - 3)
    assert len(rows) == 1, rows
    r = rows[0]
    assert (r["ref_aa"], r["alt_aa"], r["strand"], r["call"]) == ("M", "Q", "-", "ALT"), r
    assert ME not in _flags(r["flags"]), r


def test_two_alt_rows_at_one_codon_are_flagged_from_their_own_peptides() -> None:
    rows = _rows(_table(), "c3", C3_START + 3 * C3_SITE)
    by_alt = {r["alt_aa"]: r for r in rows}
    assert sorted(by_alt) == ["Q", "W"], rows
    assert all(r["ref_aa"] == "G" for r in rows), rows
    assert ME in _flags(by_alt["Q"]["flags"]), by_alt["Q"]
    assert ME not in _flags(by_alt["W"]["flags"]), by_alt["W"]
