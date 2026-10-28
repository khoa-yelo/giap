"""M3-15 lead decisions after code review (PR #63), items 1, 2 and 4.

1. Only searched peptides enter the unmapped bucket: engine decoys (``is_decoy``) and peptides with
   ``q_value > q_value_max`` are excluded inside ``unmapped_table`` (keyword ``q_value_max``,
   default ``input.q_value_max``). ``q_value == q_value_max`` is kept (inclusive, as in the config
   comment "peptides with q above this are not evidence"). Test-author choice: a missing q-value
   (NaN; engines that report only a score) is not "above" the cut-off, so the peptide is kept.
2. ``qseqid`` is resolved with ``resolve_protein_id`` like ``sseqid``: a raw ``sp|ACC|NAME`` query
   header matches ``protein_state.protein_id == ACC``.
4. One test per mutant that survived the reviewer's mutation run.
"""

from __future__ import annotations

import importlib
import inspect
import math
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
from contaminant_oracle import blast_line

from gippa.io.protein_fasta import LetterCounts, ProteinRecord
from gippa.mapping.placement_fdr import PER_PEPTIDE_COLUMNS

PREFIXES = ["Cont_", "CON__", "CONT_", "cRAP"]


def _c() -> ModuleType:
    return importlib.import_module("gippa.mapping.contaminants")


def _peptides(rows: list[tuple[str, list[str], bool, float]]) -> pd.DataFrame:
    """(peptide, protein_ids, is_decoy, q_value)."""
    return pd.DataFrame(
        {
            "peptide": [r[0] for r in rows],
            "protein_ids": pd.Series([list(r[1]) for r in rows], dtype=object),
            "is_decoy": pd.Series([r[2] for r in rows], dtype=bool),
            "is_contaminant": pd.Series([False] * len(rows), dtype=bool),
            "q_value": pd.Series([r[3] for r in rows], dtype=float),
        }
    )


def _placements(rows: list[tuple[str, bool]]) -> pd.DataFrame:
    """(peptide, accepted) -> PER_PEPTIDE_COLUMNS frame."""
    data: dict[str, list[Any]] = {
        "peptide": [p for p, _ in rows],
        "source": ["database"] * len(rows),
        "decoy_type": ["reverse"] * len(rows),
        "length": [len(p) for p, _ in rows],
        "m_target": [0] * len(rows),
        "e_target": [1e-6] * len(rows),
        "placement_q": [0.001 if a else 0.2 for _, a in rows],
        "accepted": [a for _, a in rows],
        "flags": [""] * len(rows),
    }
    df = pd.DataFrame({c: data[c] for c in PER_PEPTIDE_COLUMNS}, columns=list(PER_PEPTIDE_COLUMNS))
    df["accepted"] = df["accepted"].astype(bool)
    return df


def _proteins(rows: list[tuple[str, str, bool, bool]]) -> pd.DataFrame:
    """(protein_id, sequence, assigned, has_alignment)."""
    return pd.DataFrame(
        {
            "protein_id": [r[0] for r in rows],
            "sequence": [r[1] for r in rows],
            "assigned": pd.Series([r[2] for r in rows], dtype=bool),
            "has_alignment": pd.Series([r[3] for r in rows], dtype=bool),
        }
    )


def _table(
    peptides: pd.DataFrame | None = None,
    placements: pd.DataFrame | None = None,
    proteins: pd.DataFrame | None = None,
    pep_matches: dict[str, str] | None = None,
    prot_matches: dict[str, str] | None = None,
    prefixes: list[str] | None = None,
    **kw: Any,
) -> list[Any]:
    return list(
        _c().unmapped_table(
            peptides if peptides is not None else _peptides([]),
            placements if placements is not None else _placements([]),
            proteins if proteins is not None else _proteins([]),
            pep_matches or {},
            prot_matches or {},
            PREFIXES if prefixes is None else prefixes,
            **kw,
        )
    )


def _ids(rows: list[Any], kind: str = "peptide") -> list[str]:
    return [r.sequence_id for r in rows if r.kind == kind]


def _rec(rid: str, seq: str) -> ProteinRecord:
    return ProteinRecord(rid, (rid,), rid, seq, 0, LetterCounts(), LetterCounts())


def _write_blast(tmp_path: Path, lines: list[str]) -> Path:
    p = tmp_path / "blastp.tsv"
    p.write_text("".join(line + "\n" for line in lines))
    return p


# ------------------------------------------------------- decision 1: decoys and q > q_value_max


def test_q_value_max_keyword_defaults_to_config() -> None:
    from gippa.config.model import InputConfig

    param = inspect.signature(_c().unmapped_table).parameters["q_value_max"]
    assert param.kind is inspect.Parameter.KEYWORD_ONLY
    assert param.default == InputConfig().q_value_max == 0.01


def test_reviewer_decoy_rplsvtavk_not_in_bucket() -> None:
    # the reviewer's case: an engine decoy with no target hit used to land in the bucket as no_hit
    peps = _peptides([("RPLSVTAVK", ["rev_P12345"], True, 0.001)])
    assert _table(peps) == []


def test_decoys_excluded_even_when_contaminant_like() -> None:
    peps = _peptides(
        [
            ("SLNNQFASFIDK", ["CON__P04264"], True, 0.0),  # decoy, library match, prefixed IDs
            ("QWERTYHNMPAK", [], True, 0.0),  # decoy with a rejected target placement
        ]
    )
    plc = _placements([("QWERTYHNMPAK", False)])
    rows = _table(peps, plc, pep_matches={"SLNNQFASFIDK": "P04264"})
    assert rows == []


def test_q_above_max_excluded_at_max_kept() -> None:
    peps = _peptides(
        [
            ("AAAAAAAAAAAK", [], False, 0.01),  # exactly q_value_max: kept
            ("CCCCCCCCCCCK", [], False, 0.0100001),  # above: excluded
            ("DDDDDDDDDDDK", [], False, 0.5),  # far above: excluded
            ("EEEEEEEEEEEK", [], False, 0.0),
        ]
    )
    assert _ids(_table(peps)) == ["AAAAAAAAAAAK", "EEEEEEEEEEEK"]


def test_q_value_max_keyword_used() -> None:
    peps = _peptides(
        [
            ("AAAAAAAAAAAK", [], False, 0.05),
            ("CCCCCCCCCCCK", [], False, 0.0500001),
            ("DDDDDDDDDDDK", [], False, 0.01),
        ]
    )
    assert _ids(_table(peps, q_value_max=0.05)) == ["AAAAAAAAAAAK", "DDDDDDDDDDDK"]
    assert _ids(_table(peps, q_value_max=0.001)) == []


def test_q_above_max_excluded_for_every_reason_and_label() -> None:
    peps = _peptides(
        [
            ("VATVSLPR", [], False, 0.2),  # library match, too_short
            ("FFFFFFFFFFFK", ["CON__P02768"], False, 0.2),  # contaminant_id
            ("GGGGGGGGGGGK", [], False, 0.2),  # rejected placement
        ]
    )
    plc = _placements([("GGGGGGGGGGGK", False)])
    assert _table(peps, plc, pep_matches={"VATVSLPR": "P00761"}) == []


def test_missing_q_value_kept() -> None:
    peps = _peptides([("AAAAAAAAAAAK", [], False, math.nan)])
    (row,) = _table(peps)
    assert (row.sequence_id, row.label, row.reason) == ("AAAAAAAAAAAK", "unmapped", "no_hit")


def test_q_and_decoy_filter_does_not_touch_proteins() -> None:
    prots = _proteins([("P1", "M" + "A" * 60, False, False)])
    peps = _peptides([("RPLSVTAVK", [], True, 0.9)])
    rows = _table(peps, proteins=prots, q_value_max=0.0)
    assert [(r.kind, r.sequence_id) for r in rows] == [("protein", "P1")]


@settings(max_examples=150, deadline=None)
@given(
    data=st.lists(
        st.tuples(
            st.booleans(),  # is_decoy
            st.sampled_from([0.0, 0.005, 0.01, 0.0100001, 0.02, 0.5, math.nan]),
            st.sampled_from(["none", "rejected", "accepted"]),
        ),
        max_size=12,
    ),
    qmax=st.sampled_from([0.0, 0.01, 0.02]),
)
def test_property_bucket_peptides_are_searched_and_unplaced(
    data: list[tuple[bool, float, str]], qmax: float
) -> None:
    peps = [f"{'ACDEFGHKLMNPQRSTVWY'[i]}" * 11 + "K" for i in range(len(data))]
    frame = _peptides([(p, [], d, q) for p, (d, q, _) in zip(peps, data, strict=True)])
    plc = _placements(
        [(p, s == "accepted") for p, (_, _, s) in zip(peps, data, strict=True) if s != "none"]
    )
    expected = sorted(
        p
        for p, (decoy, q, state) in zip(peps, data, strict=True)
        if not decoy and not (q > qmax) and state != "accepted"
    )
    assert _ids(_table(frame, plc, q_value_max=qmax)) == expected


# ------------------------------------------------------------ decision 2: qseqid resolution


def test_blastp_raw_uniprot_qseqid_resolved(tmp_path: Path) -> None:
    p = _write_blast(
        tmp_path,
        [
            blast_line(qseqid="sp|Q9XYZ1|PROT_ECOLI", sseqid="sp|P02769|ALBU_BOVIN"),
            blast_line(qseqid="tr|A0A000|A0A000_9BACT", sseqid="P00761"),
        ],
    )
    assert _c().match_contaminant_proteins(p) == {"Q9XYZ1": "P02769", "A0A000": "P00761"}


def test_blastp_raw_header_query_matches_protein_state(tmp_path: Path) -> None:
    c = _c()
    p = _write_blast(
        tmp_path, [blast_line(qseqid="sp|Q9XYZ1|PROT_ECOLI", sseqid="sp|P02769|ALBU_BOVIN")]
    )
    prot_matches = c.match_contaminant_proteins(p)
    prots = _proteins([("Q9XYZ1", "M" + "A" * 60, False, True)])
    (row,) = _table(proteins=prots, prot_matches=prot_matches)
    assert (row.sequence_id, row.label, row.reason, row.contaminant_match) == (
        "Q9XYZ1",
        "contaminant",
        "below_identity",
        "P02769",
    )


def test_blastp_raw_and_plain_qseqid_same_protein_tie_rule(tmp_path: Path) -> None:
    # the raw header and the bare accession are the same query: best bitscore wins across both
    p = _write_blast(
        tmp_path,
        [
            blast_line(qseqid="sp|Q9XYZ1|PROT_ECOLI", sseqid="LIB_B", bitscore="300"),
            blast_line(qseqid="Q9XYZ1", sseqid="LIB_A", bitscore="200"),
        ],
    )
    assert _c().match_contaminant_proteins(p) == {"Q9XYZ1": "LIB_B"}


# ---------------------------------------------------------- decision 4: surviving mutants


def test_repeated_protein_rows_dedup() -> None:
    # mutant: `if bool(assigned) or pid in seen` -> `if bool(assigned)`
    seq = "M" + "A" * 60
    prots = _proteins([("P1", seq, False, False), ("P1", seq, False, False)])
    rows = _table(proteins=prots)
    assert _ids(rows, "protein") == ["P1"]


def test_protein_ids_merged_across_repeated_peptide_rows() -> None:
    # mutant: `cur.extend(...)` -> `cur[:] = _as_ids(ids)` (last row wins)
    pep = "AAAAAAAAAAAK"
    peps = _peptides([(pep, ["GSIM_1"], False, 0.0), (pep, ["CON__P02768"], False, 0.0)])
    (row,) = _table(peps)
    assert (row.label, row.reason, row.contaminant_match) == ("unmapped", "no_hit", None)
    # reversed order: the first ID of the merged list is the prefixed one, still not all prefixed
    peps = _peptides([(pep, ["CON__P02768"], False, 0.0), (pep, ["GSIM_1"], False, 0.0)])
    (row,) = _table(peps)
    assert (row.label, row.reason) == ("unmapped", "no_hit")
    # all prefixed across rows: first ID of the merged list is the match
    peps = _peptides([(pep, ["CON__B"], False, 0.0), (pep, ["Cont_A", "CON__B"], False, 0.0)])
    (row,) = _table(peps)
    assert (row.label, row.reason, row.contaminant_match) == (
        "contaminant",
        "contaminant_id",
        "CON__B",
    )


@pytest.mark.parametrize("prefixes", [[""], ["", "CON__"], ["CON__", ""]])
def test_empty_string_prefix_rejected(prefixes: list[str]) -> None:
    # mutant: `tuple(p for p in contaminant_prefixes if p)` -> `tuple(contaminant_prefixes)`
    peps = _peptides([("AAAAAAAAAAAK", ["GSIM_1"], False, 0.0)])
    prots = _proteins([("P1", "M" + "A" * 60, False, False)])
    rows = _table(peps, proteins=prots, prefixes=prefixes)
    assert {r.label for r in rows} == {"unmapped"}
    assert len(rows) == 2


@pytest.mark.parametrize("qlen", ["0", "-5"])
def test_blastp_qlen_not_positive_raises(tmp_path: Path, qlen: str) -> None:
    # mutant: `if qlen <= 0: raise` -> `if False:` (qlen 0 would be ZeroDivisionError,
    # negative qlen a silently negative coverage)
    line = blast_line(qseqid="q1", sseqid="s1").split("\t")
    line[12] = qlen
    p = _write_blast(tmp_path, ["\t".join(line)])
    with pytest.raises(ValueError, match="line 1"):
        _c().match_contaminant_proteins(p)


def test_lowercase_peptide_matched_case_insensitive() -> None:
    # mutant: `seq.upper().replace(...)` -> `seq.replace(...)`
    lib = [_rec("LIB1", "MKWVTFISLLLLFSSAYSRGVFRR")]
    got = _c().match_contaminant_peptides(["wvtfisllllf", "gvfrr"], lib)
    assert got == {"wvtfisllllf": "LIB1", "gvfrr": "LIB1"}


@pytest.mark.parametrize(
    "rows",
    [
        [True, False],
        [False, True],
        [False, True, False],
    ],
)
def test_any_accepted_placement_row_keeps_peptide_out(rows: list[bool]) -> None:
    # mutant: `state.get(pep, False) or bool(acc)` -> `bool(acc)` (last row wins)
    pep = "AAAAAAAAAAAK"
    peps = _peptides([(pep, [], False, 0.0)])
    plc = _placements([(pep, a) for a in rows])
    assert _table(peps, plc) == []


def test_all_rejected_placement_rows_give_rejected_reason() -> None:
    pep = "AAAAAAAAAAAK"
    peps = _peptides([(pep, [], False, 0.0)])
    plc = _placements([(pep, False), (pep, False)])
    (row,) = _table(peps, plc)
    assert row.reason == "rejected_placement_fdr"
