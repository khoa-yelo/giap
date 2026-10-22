"""M7-07 review (PR #165) B1: a no-start ORF's length comes from its first peptide's mapped
position, so a tolerant-route peptide (not verbatim in the genome translation) neither shortens
the ORF nor stops the step."""

from __future__ import annotations

import pytest
from helpers.novel_oracle import Case, hit_row, translate
from helpers.novel_outputs_fixture import orfs_of

from giap.model.coords import Replicon
from giap.novel.outputs import gff_features, orf_table

N = 300
_BODY = "GCT" * 10 + "GCTATTGAAAAACGTGATTCTCCGCATTGG" + "TTCAACCAGGACTACAAGCGCGTTGAGCAC"
_BODY += "TTCAACCAG" * 6
_SEQ = "TAA" + _BODY + "TGA"
_SEQ += "GCT" * ((N - len(_SEQ)) // 3)
_REPS = {"c": Replicon("c", N, "circular", "chromosome")}
_STOP_END = len("TAA" + _BODY) + 3  # exclusive end of the stop codon
_TOL = "W" + translate(_SEQ[33:63])[1:]  # one substitution against the genome translation


def _hits(with_exact: bool) -> list[dict[str, object]]:
    hits = [hit_row(_TOL, "c", 33, 63, "+", "intergenic", N, route="tolerant")]
    if with_exact:
        hits.append(hit_row(translate(_SEQ[63:93]), "c", 63, 93, "+", "intergenic", N))
    return hits


@pytest.mark.parametrize("with_exact", [True, False], ids=["tolerant+exact", "tolerant-only"])
def test_no_start_length_uses_the_tolerant_first_peptide(with_exact: bool) -> None:
    case = Case({"c": _SEQ}, _REPS, _hits(with_exact))
    (orf,) = orfs_of(case)
    assert not orf.starts and orf.first_peptide == 33
    table = orf_table([orf], None, genome=case.genome, replicons=_REPS)
    assert int(table["length_aa_longest"].iloc[0]) == (_STOP_END - 3 - 33) // 3
    (feat,) = gff_features(table)
    assert (feat.interval.start, feat.interval.end) == (33, _STOP_END)
