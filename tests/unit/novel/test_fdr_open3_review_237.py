"""M7-11 picked pairing of open-3' ORFs (code review #237 B1; lead decision 13).

An ORF whose downstream walk hits a linear contig end has ``end_3p_open`` and ``Orf.stop`` = the
lowest coordinate of its last whole codon (M7-06). makedb writes the same reading frame as the
segment ``[s, e)`` ending at that codon's end (no closing stop). The picked pairing key must match
the two: on ``+`` the ORF key is ``o.stop + 3`` = ``e``; on ``-`` it is ``o.stop - 3`` = ``s - 3``.

Both sides are built by the real code: makedb ``six_frame_orfs`` for the entries and M7-06
``candidate_rows`` / ``reconstruct`` for the ORF, on a linear 300-nt contig with no stop codon
downstream of the peptides (the M7-06 ``open3`` example). The review's repro, ``TAA + GCT*40``,
gives makedb ``[3, 123)`` and ``Orf.stop = 120``; here ``+`` gives ``[201, 300)`` and
``Orf.stop = 297``, ``-`` its mirror ``[0, 99)`` and ``Orf.stop = 0``.
"""

from __future__ import annotations

import math

import pytest
from helpers.novel_fdr_fixture import pep_row, peptides_frame
from helpers.novel_oracle import Case, Seq, hit_row, mirror_case, translate

from giap.config.model import GiapConfig
from giap.makedb.build import six_frame_orfs
from giap.mapping.sixframe import build_sixframe
from giap.model.coords import Replicon
from giap.novel.candidates import candidate_rows, reconstruct
from giap.novel.fdr import picked_class_fdr

N = 300
LIN = {"lin": Replicon("lin", N, "linear", "plasmid")}
PEP_CODONS = ["GCT", "ATT", "GAA", "AAA", "CGT", "GAT", "TCT", "CCG", "CAT", "TGG"]
TAIL = ["TTC", "AAC", "CAG", "GAC", "TAC", "AAG"]


def open3_case() -> Case:
    """Linear 300 nt: stop TAA [198, 201), ATG [210, 213), peptides [240, 270), [270, 288); no
    stop downstream before the contig end (last whole codon [297, 300))."""
    s = Seq(N).put(198, "TAA").fill(201, N).put(210, "ATG")
    seq = s.put(240, "".join(PEP_CODONS)).put(270, "".join(TAIL)).text()
    hits = [
        hit_row(translate(seq[240:270]), "lin", 240, 270, "+", "intergenic", N),
        hit_row(translate(seq[270:288]), "lin", 270, 288, "+", "intergenic", N),
    ]
    return Case({"lin": seq}, dict(LIN), hits)


def _world(strand: str):  # type: ignore[no-untyped-def]
    case = open3_case() if strand == "+" else mirror_case(open3_case())
    rows = candidate_rows(case.hits_df(), case.peptides_df(), None)
    (orf,) = reconstruct(rows, case.gm, case.genome, GiapConfig(seed=1))
    index = build_sixframe(dict(case.genome), case.gm.replicons)
    entries = six_frame_orfs(index, case.gm, dict(case.genome), 10)
    return case, orf, entries


@pytest.mark.parametrize("strand", ["+", "-"])
def test_open3_fixture_has_the_reviewed_geometry(strand: str) -> None:
    """Guards the setup: the ORF is 3'-open and makedb writes the same frame up to the end."""
    _, orf, entries = _world(strand)
    assert orf.end_3p_open is True
    assert orf.strand == strand
    assert orf.stop == (N - 3 if strand == "+" else 0)
    seg = [(e.start, e.end) for e in entries if e.strand == strand and e.cls == "intergenic"]
    assert ((201, N) if strand == "+" else (0, N - 201)) in seg


def _pair_and_compete(strand: str):  # type: ignore[no-untyped-def]
    _, orf, entries = _world(strand)
    want = (201, N) if strand == "+" else (0, N - 201)
    (entry,) = [e for e in entries if e.strand == strand and (e.start, e.end) == want]
    did = entry.db_id
    # Target peptides of the ORF beat the entry's own decoy (input q 1e-5 < 1e-3).
    rows = [pep_row(p, [did], 1 + k, is_decoy=False) for k, p in enumerate(orf.peptides)]
    rows.append(pep_row("DECYPEPTDE", ["rev_" + did], 100, is_decoy=True))
    id_map = None  # entries read from the tagged IDs (decision: id_map=None)
    res = picked_class_fdr([orf], peptides_frame(rows), id_map)
    return orf, entry, res


@pytest.mark.parametrize("strand", ["+", "-"])
def test_open3_orf_pairs_with_its_makedb_entry(strand: str) -> None:
    """B1: the open-3' ORF is paired (no ``no_makedb_pair``) and wins its class (q = 1/1)."""
    orf, entry, res = _pair_and_compete(strand)
    assert res.reason.get(orf.orf_id) is None, (
        f"{strand}: Orf.stop={orf.stop}, makedb [{entry.start}, {entry.end}) -> "
        f"reason {res.reason.get(orf.orf_id)!r}"
    )
    assert not math.isnan(res.q[orf.orf_id])
    assert res.q[orf.orf_id] == pytest.approx(1.0)  # (0 decoy wins / f + 1) / 1 target win
    st = res.classes[entry.cls]
    assert (st.n_target_wins, st.n_decoy_wins, st.estimated) == (1, 0, True)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_open3_orf_pairs_with_replicon_lengths_given(strand: str) -> None:
    """B1 with ``replicon_lengths``: a linear contig's open end pairs the same way."""
    _, orf, entries = _world(strand)
    want = (201, N) if strand == "+" else (0, N - 201)
    (entry,) = [e for e in entries if e.strand == strand and (e.start, e.end) == want]
    did = entry.db_id
    rows = [pep_row(p, [did], 1 + k, is_decoy=False) for k, p in enumerate(orf.peptides)]
    rows.append(pep_row("DECYPEPTDE", ["rev_" + did], 100, is_decoy=True))
    res = picked_class_fdr([orf], peptides_frame(rows), None, replicon_lengths={"lin": N})
    assert res.reason.get(orf.orf_id) is None
    assert res.q[orf.orf_id] == pytest.approx(1.0)
