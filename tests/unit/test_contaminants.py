"""M3-15: unmapped bucket and contaminants (``giap.mapping.contaminants``), examples 1-10.

Input contracts used here (the spec names the frames but not their columns; test-author choices,
reported to the lead):

- ``peptides``: the aggregated query peptides (M1-06/M1-12 columns ``peptide``, ``protein_ids`` =
  native list of strings, empty when missing (M3-15 lead decision 1), ``is_decoy``,
  ``is_contaminant``); every row is a peptide that went through mapping.
- ``placements``: M3-10 ``placement`` per-peptide frame (``PER_PEPTIDE_COLUMNS``); a peptide is in
  it iff it has a target hit; ``accepted`` is the placement decision.
- ``protein_state``: one row per mapped protein (> 40 aa): ``protein_id``, ``sequence``,
  ``assigned`` (bool, assigned by miniprot or tblastn), ``has_alignment`` (bool, any alignment from
  either route). Assigned proteins never enter the bucket.
- ``unmapped_table(..., tolerant_min_len=10)``: keyword (default 10 = ``mapping.tolerant.min_len``)
  so ``too_short`` follows the configured tolerant minimum (M3-15 Q2 / M3-index decision).
- blastp tabular output uses the 15 tblastn columns (``giap.mapping.tblastn.OUTFMT``); identity =
  ``pident / 100``, coverage = ``(qend - qstart + 1) / qlen`` of one HSP; both thresholds inclusive;
  several passing hits: highest bitscore, identity, coverage, then lowest ID (lead decision 2).
"""

from __future__ import annotations

import dataclasses
import importlib
import random
import shlex
import sys
from pathlib import Path
from types import ModuleType
from typing import Any, get_args

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
from contaminant_oracle import blast_line

from giap.io.protein_fasta import ProteinRecord, read_protein_fasta
from giap.mapping.placement_fdr import PER_PEPTIDE_COLUMNS

MINI = Path(__file__).resolve().parents[1] / "data" / "contaminants_mini.faa"
TRYPSIN, KRT1, KRT10, ALBUMIN = "P00761", "P04264", "P13645", "P02769"
PREFIXES = ["Cont_", "CON__", "CONT_", "cRAP"]
UNMAPPED_FIELDS = ("sequence_id", "sequence", "kind", "label", "reason", "contaminant_match")


def _c() -> ModuleType:
    return importlib.import_module("giap.mapping.contaminants")


@pytest.fixture(scope="module")
def library() -> list[ProteinRecord]:
    return read_protein_fasta(MINI)


def _peptides(rows: list[tuple[str, str | None]]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "peptide": [p for p, _ in rows],
            "protein_ids": pd.Series(
                [ids.split(";") if ids else [] for _, ids in rows], dtype=object
            ),
            "is_decoy": [False] * len(rows),
            "is_contaminant": [
                ids is not None and all(i.startswith(tuple(PREFIXES)) for i in ids.split(";"))
                for _, ids in rows
            ],
        }
    )


def _placements(rows: list[tuple[str, float, bool]]) -> pd.DataFrame:
    """(peptide, placement_q, accepted) -> a PER_PEPTIDE_COLUMNS frame (database, m = 0)."""
    data: dict[str, list[Any]] = {
        "peptide": [p for p, _, _ in rows],
        "source": ["database"] * len(rows),
        "decoy_type": ["reverse"] * len(rows),
        "length": [len(p) for p, _, _ in rows],
        "m_target": [0] * len(rows),
        "e_target": [1e-6] * len(rows),
        "placement_q": [q for _, q, _ in rows],
        "accepted": [a for _, _, a in rows],
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


def _empty_placements() -> pd.DataFrame:
    return _placements([])


def _empty_proteins() -> pd.DataFrame:
    return _proteins([])


def _table(
    peptides: pd.DataFrame | None = None,
    placements: pd.DataFrame | None = None,
    proteins: pd.DataFrame | None = None,
    pep_matches: dict[str, str] | None = None,
    prot_matches: dict[str, str] | None = None,
    **kw: Any,
) -> list[Any]:
    return list(
        _c().unmapped_table(
            peptides if peptides is not None else _peptides([]),
            placements if placements is not None else _empty_placements(),
            proteins if proteins is not None else _empty_proteins(),
            pep_matches or {},
            prot_matches or {},
            PREFIXES,
            **kw,
        )
    )


def _by_id(rows: list[Any]) -> dict[str, Any]:
    out = {r.sequence_id: r for r in rows}
    assert len(out) == len(rows), "sequence_id repeated"
    return out


def _write_blast(tmp_path: Path, lines: list[str]) -> Path:
    p = tmp_path / "blastp.tsv"
    p.write_text("".join(line + "\n" for line in lines))
    return p


# ------------------------------------------------------------------------------- interface


def test_unmapped_row_fields_and_frozen() -> None:
    c = _c()
    assert dataclasses.is_dataclass(c.UnmappedRow)
    assert tuple(f.name for f in dataclasses.fields(c.UnmappedRow)) == UNMAPPED_FIELDS
    row = c.UnmappedRow("PEPTIDEK", "PEPTIDEK", "peptide", "unmapped", "no_hit", None)
    with pytest.raises(dataclasses.FrozenInstanceError):
        row.label = "contaminant"  # type: ignore[misc]
    assert not hasattr(row, "__dict__")  # slots=True


def test_reason_vocabulary() -> None:
    assert set(get_args(_c().Reason)) == {
        "no_hit",
        "rejected_placement_fdr",
        "below_identity",
        "too_short",
        "contaminant_id",
    }


# ------------------------------------------------------------------------- peptide_reason


@pytest.mark.parametrize(
    ("has_hit", "accepted", "length", "expected"),
    [
        (True, True, 12, None),  # placed and accepted: never in the bucket
        (True, True, 8, None),
        (True, False, 12, "rejected_placement_fdr"),  # example 4
        (True, False, 8, "rejected_placement_fdr"),  # short exact hit, rejected
        (False, False, 8, "too_short"),  # example 6
        (False, False, 9, "too_short"),
        (False, False, 10, "no_hit"),  # = tolerant minimum: was searched tolerantly
        (False, False, 14, "no_hit"),  # example 7
    ],
)
def test_peptide_reason(has_hit: bool, accepted: bool, length: int, expected: str | None) -> None:
    assert _c().peptide_reason(has_hit, accepted, length, 10) == expected


def test_peptide_reason_follows_tolerant_min_len() -> None:
    c = _c()
    assert c.peptide_reason(False, False, 11, 12) == "too_short"
    assert c.peptide_reason(False, False, 12, 12) == "no_hit"
    assert c.peptide_reason(False, False, 7, 7) == "no_hit"


# ---------------------------------------------------------------- match_contaminant_peptides


def test_fixture_library_ids(library: list[ProteinRecord]) -> None:
    assert sorted(r.id for r in library) == sorted([TRYPSIN, KRT1, KRT10, ALBUMIN])


def test_example1_trypsin_peptide_matched(library: list[ProteinRecord]) -> None:
    assert _c().match_contaminant_peptides(["VATVSLPR"], library) == {"VATVSLPR": TRYPSIN}


def test_example2_keratin_peptide_matched(library: list[ProteinRecord]) -> None:
    assert _c().match_contaminant_peptides(["SLNNQFASFIDK"], library) == {"SLNNQFASFIDK": KRT1}


def test_example3_i_equals_l_query_side(library: list[ProteinRecord]) -> None:
    # library has SLNNQ... ; the peptide is written with I
    got = _c().match_contaminant_peptides(["SINNQFASFIDK", "SLNNQFASFLDK"], library)
    assert got == {"SINNQFASFIDK": KRT1, "SLNNQFASFLDK": KRT1}


def test_i_equals_l_library_side(tmp_path: Path) -> None:
    fa = tmp_path / "lib.faa"
    fa.write_text(">CONX1 x\nMKAAIIGGRK\n")
    lib = read_protein_fasta(fa)
    assert _c().match_contaminant_peptides(["AALLGGR", "AALIGGR"], lib) == {
        "AALLGGR": "CONX1",
        "AALIGGR": "CONX1",
    }


def test_unmatched_peptides_absent(library: list[ProteinRecord]) -> None:
    got = _c().match_contaminant_peptides(["PEPTIDEKAAW", "VATVSLPR"], library)
    assert got == {"VATVSLPR": TRYPSIN}


def test_several_entries_lowest_id(tmp_path: Path) -> None:
    fa = tmp_path / "lib.faa"
    fa.write_text(">Z9 z\nMMMSHAREDPEPKMMM\n>B2 b\nWWSHAREDPEPKWW\n>Q5 q\nSHAREDPEPK\n")
    lib = read_protein_fasta(fa)
    c = _c()
    assert c.match_contaminant_peptides(["SHAREDPEPK"], lib) == {"SHAREDPEPK": "B2"}
    # independent of library order
    assert c.match_contaminant_peptides(["SHAREDPEPK"], list(reversed(lib))) == {"SHAREDPEPK": "B2"}


def test_library_absent_no_matches() -> None:
    assert _c().match_contaminant_peptides(["VATVSLPR", "SLNNQFASFIDK"], []) == {}


def test_no_peptides() -> None:
    assert _c().match_contaminant_peptides([], read_protein_fasta(MINI)) == {}


# ------------------------------------------------------------------------------ blastp_cmd


def test_blastp_cmd(tmp_path: Path) -> None:
    from giap.mapping.tblastn import OUTFMT

    q, lib = tmp_path / "leftovers.faa", tmp_path / "contaminants.fasta"
    cmd = _c().blastp_cmd(q, lib, 4)
    assert all(isinstance(x, str) for x in cmd)
    assert Path(cmd[0]).name == "blastp"
    assert "-query" in cmd and cmd[cmd.index("-query") + 1] == str(q)
    assert str(lib) in cmd or any(str(lib) in x for x in cmd)
    assert "-outfmt" in cmd and cmd[cmd.index("-outfmt") + 1] == OUTFMT
    assert "4" in cmd
    assert shlex.join(cmd)  # a plain argument vector


# ---------------------------------------------------------------- match_contaminant_proteins


def test_example9_albumin_92_identity_60_coverage(tmp_path: Path) -> None:
    p = _write_blast(
        tmp_path,
        [blast_line(qseqid="PROT1", sseqid=ALBUMIN, pident=92.0, qstart=1, qend=60, qlen=100)],
    )
    assert _c().match_contaminant_proteins(p) == {"PROT1": ALBUMIN}


def test_example10_low_coverage_not_contaminant(tmp_path: Path) -> None:
    p = _write_blast(
        tmp_path,
        [blast_line(qseqid="PROT1", sseqid=ALBUMIN, pident=95.0, qstart=11, qend=50, qlen=100)],
    )
    assert _c().match_contaminant_proteins(p) == {}


@pytest.mark.parametrize(
    ("pident", "qstart", "qend", "ok"),
    [
        (90.0, 1, 50, True),  # both thresholds inclusive
        (89.9, 1, 100, False),
        (100.0, 1, 49, False),
        (100.0, 51, 100, True),
        (70.0, 1, 100, False),  # example 8's identity
    ],
)
def test_protein_thresholds(
    tmp_path: Path, pident: float, qstart: int, qend: int, ok: bool
) -> None:
    p = _write_blast(
        tmp_path,
        [blast_line(qseqid="P", sseqid=ALBUMIN, pident=pident, qstart=qstart, qend=qend, qlen=100)],
    )
    assert _c().match_contaminant_proteins(p) == ({"P": ALBUMIN} if ok else {})


def test_protein_custom_thresholds(tmp_path: Path) -> None:
    p = _write_blast(
        tmp_path, [blast_line(qseqid="P", sseqid="S", pident=80.0, qstart=1, qend=30, qlen=100)]
    )
    c = _c()
    assert c.match_contaminant_proteins(p) == {}
    assert c.match_contaminant_proteins(p, min_identity=0.8, min_coverage=0.3) == {"P": "S"}


def test_protein_passing_subject_chosen_over_failing(tmp_path: Path) -> None:
    p = _write_blast(
        tmp_path,
        [
            blast_line(qseqid="P", sseqid="A0", pident=99.0, qstart=1, qend=20, qlen=100),
            blast_line(qseqid="P", sseqid="Z9", pident=95.0, qstart=1, qend=90, qlen=100),
        ],
    )
    assert _c().match_contaminant_proteins(p) == {"P": "Z9"}


@pytest.mark.parametrize(
    ("a", "b", "winner"),
    [
        # (sseqid, bitscore, pident, qend) of two passing HSPs of one query (qlen 100)
        (("A1", 300, 91.0, 60), ("Z9", 400, 91.0, 60), "Z9"),  # highest bitscore
        (("A1", 300, 91.0, 60), ("Z9", 300, 97.0, 60), "Z9"),  # then highest identity
        (("A1", 300, 95.0, 60), ("Z9", 300, 95.0, 90), "Z9"),  # then highest coverage
        (("Z9", 300, 95.0, 80), ("A1", 300, 95.0, 80), "A1"),  # then lowest subject ID
    ],
)
def test_protein_tie_rule(
    tmp_path: Path,
    a: tuple[str, int, float, int],
    b: tuple[str, int, float, int],
    winner: str,
) -> None:
    lines = [
        blast_line(qseqid="P", sseqid=s, bitscore=bs, pident=pi, qstart=1, qend=qe, qlen=100)
        for s, bs, pi, qe in (a, b)
    ]
    c = _c()
    assert c.match_contaminant_proteins(_write_blast(tmp_path, lines)) == {"P": winner}
    assert c.match_contaminant_proteins(_write_blast(tmp_path, lines[::-1])) == {"P": winner}


def test_protein_tie_rule_failing_hsp_never_wins(tmp_path: Path) -> None:
    lines = [
        blast_line(qseqid="P", sseqid="A1", bitscore=900, pident=89.0, qend=100, qlen=100),
        blast_line(qseqid="P", sseqid="B2", bitscore=900, pident=99.0, qend=40, qlen=100),
        blast_line(qseqid="P", sseqid="Z9", bitscore=150, pident=90.0, qend=50, qlen=100),
    ]
    assert _c().match_contaminant_proteins(_write_blast(tmp_path, lines)) == {"P": "Z9"}


def test_protein_matches_independent_of_row_order(tmp_path: Path) -> None:
    lines = [
        blast_line(qseqid=f"P{i % 4}", sseqid=f"S{i}", pident=90.0 + i % 10, qlen=100, qend=50 + i)
        for i in range(20)
    ]
    c = _c()
    ref = c.match_contaminant_proteins(_write_blast(tmp_path, lines))
    assert set(ref) == {"P0", "P1", "P2", "P3"}
    for seed in range(5):
        random.Random(seed).shuffle(lines)
        assert c.match_contaminant_proteins(_write_blast(tmp_path, lines)) == ref


def test_protein_empty_output(tmp_path: Path) -> None:
    assert _c().match_contaminant_proteins(_write_blast(tmp_path, [])) == {}


def test_protein_bad_column_count_raises(tmp_path: Path) -> None:
    bad = "\t".join(blast_line(qseqid="P", sseqid="S").split("\t")[:12])
    with pytest.raises(ValueError):
        _c().match_contaminant_proteins(_write_blast(tmp_path, [bad]))


# ------------------------------------------------------------------- unmapped_table: peptides


def test_example1_row(library: list[ProteinRecord]) -> None:
    pep = "VATVSLPR"
    rows = _table(
        _peptides([(pep, None)]),
        pep_matches=_c().match_contaminant_peptides([pep], library),
    )
    assert rows == [_c().UnmappedRow(pep, pep, "peptide", "contaminant", "too_short", TRYPSIN)]


def test_example2_row() -> None:
    pep = "SLNNQFASFIDK"
    (row,) = _table(_peptides([(pep, None)]), pep_matches={pep: KRT1})
    assert (row.label, row.reason, row.contaminant_match) == ("contaminant", "no_hit", KRT1)


def test_example3_row(library: list[ProteinRecord]) -> None:
    pep = "SINNQFASFIDK"
    (row,) = _table(
        _peptides([(pep, None)]), pep_matches=_c().match_contaminant_peptides([pep], library)
    )
    assert (row.sequence_id, row.label, row.reason, row.contaminant_match) == (
        pep,
        "contaminant",
        "no_hit",
        KRT1,
    )


def test_example4_rejected_placement() -> None:
    pep = "AGHTWQEPLNVK"
    (row,) = _table(_peptides([(pep, "GSIM_00010")]), _placements([(pep, 0.03, False)]))
    assert row == _c().UnmappedRow(pep, pep, "peptide", "unmapped", "rejected_placement_fdr", None)


def test_example5_contaminant_id() -> None:
    pep = "AWHMNPQGTSEK"
    (row,) = _table(_peptides([(pep, "CON__P02768")]))
    assert row == _c().UnmappedRow(
        pep, pep, "peptide", "contaminant", "contaminant_id", "CON__P02768"
    )


def test_example6_too_short() -> None:
    (row,) = _table(_peptides([("WHMNPQGK", None)]))
    assert row == _c().UnmappedRow("WHMNPQGK", "WHMNPQGK", "peptide", "unmapped", "too_short", None)


def test_example7_no_hit() -> None:
    pep = "WHMNPQGKTWHMNR"
    assert len(pep) == 14
    (row,) = _table(_peptides([(pep, None)]))
    assert row == _c().UnmappedRow(pep, pep, "peptide", "unmapped", "no_hit", None)


def test_decoy_only_peptides() -> None:
    # a decoy hit but no target hit: absent from placements -> no_hit / too_short (M3-10 ex. 9)
    long_, short = "QWERTYHNMPAK", "QWERTYHK"
    rows = _by_id(_table(_peptides([(long_, None), (short, None)])))
    assert rows[long_].reason == "no_hit"
    assert rows[short].reason == "too_short"
    assert {r.label for r in rows.values()} == {"unmapped"}


def test_too_short_follows_tolerant_min_len_kwarg() -> None:
    pep = "WHMNPQGKTWH"  # 11 aa
    (row,) = _table(_peptides([(pep, None)]), tolerant_min_len=12)
    assert row.reason == "too_short"
    (row,) = _table(_peptides([(pep, None)]), tolerant_min_len=10)
    assert row.reason == "no_hit"


def test_accepted_peptides_never_in_bucket() -> None:
    # Q1: accepted placements stay out, even with contaminant protein_ids or a library match
    a, b = "VATVSLPR", "AWHMNPQGTSEK"
    rows = _table(
        _peptides([(a, None), (b, "CON__P02768")]),
        _placements([(a, 0.0, True), (b, 0.001, True)]),
        pep_matches={a: TRYPSIN},
    )
    assert rows == []


def test_library_match_beats_prefix_and_keeps_reason() -> None:
    pep = "SLNNQFASFIDK"
    (row,) = _table(
        _peptides([(pep, "CON__P04264")]),
        _placements([(pep, 0.2, False)]),
        pep_matches={pep: KRT1},
    )
    assert (row.label, row.reason, row.contaminant_match) == (
        "contaminant",
        "rejected_placement_fdr",
        KRT1,
    )


def test_rejected_with_contaminant_prefix() -> None:
    pep = "AGHTWQEPLNVK"
    (row,) = _table(_peptides([(pep, "CON__P00761")]), _placements([(pep, 0.5, False)]))
    assert (row.label, row.reason, row.contaminant_match) == (
        "contaminant",
        "contaminant_id",
        "CON__P00761",
    )


def test_prefix_needs_every_protein_id() -> None:
    pep = "AWHMNPQGTSEK"
    (row,) = _table(_peptides([(pep, "CON__P02768;GSIM_00010")]))
    assert (row.label, row.reason, row.contaminant_match) == ("unmapped", "no_hit", None)


def test_prefix_first_id_and_other_prefixes() -> None:
    a, b = "AWHMNPQGTSEK", "QWERTYHNMPAK"
    rows = _by_id(_table(_peptides([(a, "CON__P00761;CON__P02768"), (b, "Cont_ALBU;cRAP001")])))
    assert (rows[a].reason, rows[a].contaminant_match) == ("contaminant_id", "CON__P00761")
    assert (rows[b].reason, rows[b].contaminant_match) == ("contaminant_id", "Cont_ALBU")


def test_prefixes_come_from_argument() -> None:
    pep = "AWHMNPQGTSEK"
    peps = _peptides([(pep, "CON__P02768")])
    c = _c()
    (row,) = c.unmapped_table(peps, _empty_placements(), _empty_proteins(), {}, {}, ["XYZ_"])
    assert (row.label, row.reason, row.contaminant_match) == ("unmapped", "no_hit", None)


def test_library_absent_everything_unmapped_with_reasons() -> None:
    peps = _peptides([("VATVSLPR", None), ("SLNNQFASFIDK", None), ("AGHTWQEPLNVK", None)])
    prots = _proteins(
        [("PROT1", "M" + "A" * 80, False, False), ("PROT2", "M" + "W" * 80, False, True)]
    )
    rows = _by_id(_table(peps, _placements([("AGHTWQEPLNVK", 0.03, False)]), prots))
    assert {k: (r.label, r.reason, r.contaminant_match) for k, r in rows.items()} == {
        "VATVSLPR": ("unmapped", "too_short", None),
        "SLNNQFASFIDK": ("unmapped", "no_hit", None),
        "AGHTWQEPLNVK": ("unmapped", "rejected_placement_fdr", None),
        "PROT1": ("unmapped", "no_hit", None),
        "PROT2": ("unmapped", "below_identity", None),
    }


# ------------------------------------------------------------------- unmapped_table: proteins

SEQ_A = "M" + "KAVLDEGTWR" * 7
SEQ_B = "M" + "PQSTNHYFRE" * 6


def test_example8_below_identity() -> None:
    (row,) = _table(proteins=_proteins([("PROT8", SEQ_A, False, True)]))
    assert row == _c().UnmappedRow("PROT8", SEQ_A, "protein", "unmapped", "below_identity", None)


def test_protein_no_alignment_no_hit() -> None:
    (row,) = _table(proteins=_proteins([("PROT7", SEQ_A, False, False)]))
    assert (row.kind, row.label, row.reason) == ("protein", "unmapped", "no_hit")


@pytest.mark.parametrize(("aligned", "reason"), [(False, "no_hit"), (True, "below_identity")])
def test_example9_albumin_contaminant_keeps_reason(aligned: bool, reason: str) -> None:
    (row,) = _table(
        proteins=_proteins([("PROT9", SEQ_B, False, aligned)]), prot_matches={"PROT9": ALBUMIN}
    )
    assert row == _c().UnmappedRow("PROT9", SEQ_B, "protein", "contaminant", reason, ALBUMIN)


def test_assigned_proteins_not_in_bucket() -> None:
    rows = _table(
        proteins=_proteins([("PROTA", SEQ_A, True, True), ("PROTB", SEQ_B, False, True)]),
        prot_matches={"PROTA": ALBUMIN},
    )
    assert [r.sequence_id for r in rows] == ["PROTB"]


def test_protein_contaminant_prefix_id() -> None:
    (row,) = _table(proteins=_proteins([("CON__P02768", SEQ_B, False, False)]))
    assert (row.label, row.reason, row.contaminant_match) == (
        "contaminant",
        "contaminant_id",
        "CON__P02768",
    )


# ----------------------------------------------------------------- ordering and determinism


def _mixed_inputs() -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, str], dict[str, str]
]:
    peps = _peptides(
        [
            ("VATVSLPR", None),
            ("SLNNQFASFIDK", None),
            ("AGHTWQEPLNVK", "GSIM_00010"),
            ("AWHMNPQGTSEK", "CON__P02768"),
            ("WHMNPQGK", None),
            ("WHMNPQGKTWHMNR", None),
            ("ACCEPTEDPEPK", "GSIM_00020"),
        ]
    )
    plc = _placements([("AGHTWQEPLNVK", 0.03, False), ("ACCEPTEDPEPK", 0.0, True)])
    prots = _proteins(
        [
            ("PROT9", SEQ_B, False, False),
            ("PROT8", SEQ_A, False, True),
            ("GSIM_00030", SEQ_A, True, True),
        ]
    )
    return peps, plc, prots, {"VATVSLPR": TRYPSIN, "SLNNQFASFIDK": KRT1}, {"PROT9": ALBUMIN}


def test_sorted_by_kind_then_sequence_id() -> None:
    rows = _table(*_mixed_inputs())
    keys = [(r.kind, r.sequence_id) for r in rows]
    assert keys == sorted(keys)
    assert [r.kind for r in rows] == ["peptide"] * 6 + ["protein"] * 2
    assert "ACCEPTEDPEPK" not in {r.sequence_id for r in rows}


def test_row_value_types() -> None:
    c = _c()
    reasons = set(get_args(c.Reason))
    for r in _table(*_mixed_inputs()):
        assert isinstance(r, c.UnmappedRow)
        assert type(r.sequence_id) is str and type(r.sequence) is str
        assert r.kind in ("peptide", "protein")
        assert r.label in ("contaminant", "unmapped")
        assert r.reason in reasons
        assert (r.contaminant_match is None) == (r.label == "unmapped")
        if r.contaminant_match is not None:
            assert type(r.contaminant_match) is str
        if r.kind == "peptide":
            assert r.sequence_id == r.sequence


def test_order_independent_and_deterministic() -> None:
    peps, plc, prots, pm, prm = _mixed_inputs()
    ref = _table(peps, plc, prots, pm, prm)
    assert _table(peps, plc, prots, pm, prm) == ref
    for seed in range(5):
        assert (
            _table(
                peps.sample(frac=1, random_state=seed).reset_index(drop=True),
                plc.sample(frac=1, random_state=seed).reset_index(drop=True),
                prots.sample(frac=1, random_state=seed).reset_index(drop=True),
                dict(reversed(list(pm.items()))),
                prm,
            )
            == ref
        )


def test_inputs_not_mutated() -> None:
    peps, plc, prots, pm, prm = _mixed_inputs()
    copies = (peps.copy(), plc.copy(), prots.copy(), dict(pm), dict(prm))
    _table(peps, plc, prots, pm, prm)
    pd.testing.assert_frame_equal(peps, copies[0])
    pd.testing.assert_frame_equal(plc, copies[1])
    pd.testing.assert_frame_equal(prots, copies[2])
    assert (pm, prm) == copies[3:]


def test_empty_inputs() -> None:
    assert _table() == []
