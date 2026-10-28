"""The ``denovo_only_unsupported`` flag (docs/specs/M7-17-denovo-support-flag.md; ADR-0050
option 3; GIPPA_SPEC.md §7.6.3, §9.4): the rule, ``db_supported`` and ``apply_evidence``.

Acceptance criterion 1 (unit), criterion 4 (schema) and the ``FLAG_ORDER`` / writer half of
"The flag": examples 1-5, 7-11 and the decoy, contaminant, absent-peptide, I/L and ``auto``
edge rows.

Readings pinned by the test author:

* ``denovo_only_unsupported(cls, has_denovo, db_support, rbs)`` is
  ``has_denovo and not db_support and not (rbs is True and cls not in RBS_CLASSES)``; ``rbs``
  NA (``None``) is no support (rule 4).
* ``db_supported(orfs, peptides, sources)`` returns a pair for **every** ORF passed (an ORF with
  no matched peptide maps to ``(False, False)``).
* ``apply_evidence(..., denovo_support=None)`` equals the call without the argument; an ORF
  absent from a given mapping is not tested here (the step passes every ORF).
* "Nothing else changes": every column but ``flags`` is equal to the run without
  ``denovo_support``, and ``flags`` differs only by the new token.
"""

from __future__ import annotations

import itertools
import json
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.novel_evidence_fixture import (
    GT15,
    GT30,
    ID1,
    ID6,
    cfg,
    ev,
    ex1_case,
    flags_of,
    is_na,
    nterm2,
)
from helpers.novel_fdr_fixture import make_orf, pep_row, peptides_frame
from helpers.novel_oracle import Case
from helpers.novel_outputs_fixture import C1, ex1_hits, ex1_seq, orfs_of, scored, table
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.config.model import GippaConfig
from gippa.novel.candidates import NOVEL_CLASSES

FLAG = "denovo_only_unsupported"
RBS_CLASSES = frozenset({"intergenic", "nterm_extension"})
SCHEMA = Path(__file__).resolve().parents[3] / "src/gippa/schemas/tables/novel_orfs.schema.json"
OLD_FLAGS = (
    "junction_spanned", "homopolymer_junction", "class_fdr_unestimated", "rbs_supported",
    "no_start_codon", "partial", "ostir_missing", "rbs_threshold_unestimated", "junction_unknown",
    "promoter_supported", "terminator_supported", "readthrough_candidate",
)  # fmt: skip


def cc() -> Any:
    import gippa.novel.class_call as m

    return m


def _sources(peptides: pd.DataFrame, config: GippaConfig | None = None) -> pd.Series:
    from gippa.steps.mapping import resolve_sources

    return resolve_sources(peptides, config or cfg())


def _support(orfs: Sequence[Any], peptides: pd.DataFrame,
             config: GippaConfig | None = None) -> dict[str, tuple[bool, bool]]:  # fmt: skip
    out = cc().db_supported(orfs, peptides, _sources(peptides, config))
    return {k: (bool(a), bool(b)) for k, (a, b) in out.items()}


# ---------------------------------------------------------------- the rule (truth table)


def _expected(cls: str, has_dn: bool, db: bool, rbs: bool | None) -> bool:
    return has_dn and not db and not (rbs is True and cls not in RBS_CLASSES)


@pytest.mark.parametrize(
    ("cls", "has_dn", "db", "rbs"),
    list(itertools.product(NOVEL_CLASSES, (True, False), (True, False), (True, False, None))),
)
def test_rule_truth_table(cls: str, has_dn: bool, db: bool, rbs: bool | None) -> None:
    got = ev().denovo_only_unsupported(cls, has_dn, db, rbs)
    assert isinstance(got, bool)
    assert got is _expected(cls, has_dn, db, rbs)


@pytest.mark.parametrize("cls", sorted(RBS_CLASSES))
def test_rbs_classes_with_rbs_support_are_flagged(cls: str) -> None:
    """Q1 / lead decision 1: RBS is the class rule there, not support beyond the call."""
    assert ev().denovo_only_unsupported(cls, True, False, True) is True


@pytest.mark.parametrize("cls", ["cds_outframe", "antisense", "pseudogene", "sorf",
                                 "cterm_extension"])  # fmt: skip
def test_rbs_support_exempts_outside_the_rbs_classes(cls: str) -> None:
    assert ev().denovo_only_unsupported(cls, True, False, True) is False
    assert ev().denovo_only_unsupported(cls, True, False, False) is True
    assert ev().denovo_only_unsupported(cls, True, False, None) is True


def test_flag_is_last_in_flag_order() -> None:
    from gippa.novel.outputs import FLAG_ORDER

    assert tuple(FLAG_ORDER) == (*OLD_FLAGS, FLAG)


# ---------------------------------------------------------------- db_supported


def _orf(oid: str, peptides: Sequence[str], cls: str = "cds_outframe", stop: int = 1200) -> Any:
    return make_orf(oid, cls, "c1", "+", stop, peptides)


def _dn(p: str, k: int = 1, **kw: Any) -> dict[str, Any]:
    return pep_row(p, ["x"], k, is_decoy=False, source="denovo", q=math.nan, **kw)


def _db(p: str, k: int = 1, **kw: Any) -> dict[str, Any]:
    return pep_row(p, ["b0001"], k, is_decoy=False, source="database", **kw)


def test_example_5_mixed_db_and_dn_peptides_is_db_support() -> None:
    o = _orf("A", ["AEPLK", "GGTVSR"], cls="antisense")
    peps = peptides_frame([_db("AEPLK"), _dn("GGTVSR", 2)])
    assert _support([o], peps) == {"A": (True, True)}


def test_example_8_source_less_table_resolves_to_denovo() -> None:
    """``input.peptide_source: auto``, no q-value or PEP -> every row ``denovo``."""
    o = _orf("S", ["MKTAYK", "GGSWNK"], cls="sorf")
    rows = [pep_row(p, ["x"], k, is_decoy=False, source="auto", q=math.nan)
            for k, p in enumerate(o.peptides, 1)]  # fmt: skip
    peps = peptides_frame(rows)
    peps["q_value"] = math.nan
    assert _support([o], peps) == {"S": (True, False)}


def test_auto_source_with_a_q_value_resolves_to_database() -> None:
    o = _orf("S", ["MKTAYK", "GGSWNK"], cls="sorf")
    rows = [pep_row(p, ["x"], k, is_decoy=False, source="auto")
            for k, p in enumerate(o.peptides, 1)]  # fmt: skip
    assert _support([o], peptides_frame(rows)) == {"S": (False, True)}


def test_auto_source_follows_input_peptide_source() -> None:
    o = _orf("S", ["MKTAYK"], cls="sorf")
    peps = peptides_frame([pep_row("MKTAYK", ["x"], 1, is_decoy=False, source="auto")])
    c = GippaConfig(seed=1, input={"peptide_source": "denovo"})
    assert _support([o], peps, c) == {"S": (True, False)}


def test_example_9_mixed_source_peptide_aggregated_database() -> None:
    o = _orf("M", ["ELVLLK", "GGTVSR"])
    peps = peptides_frame([_db("ELVLLK"), _dn("GGTVSR", 2)])
    assert _support([o], peps) == {"M": (True, True)}


def test_example_10_il_folded_match_is_db_support() -> None:
    o = _orf("F", ["LLGK"])
    assert o.peptides == ("LLGK",)
    assert _support([o], peptides_frame([_db("ILGK")])) == {"F": (False, True)}


def test_folded_key_with_several_rows_takes_the_union() -> None:
    o = _orf("U", ["ELVLLK"])
    peps = peptides_frame([_dn("ELVILK"), _db("ELVLLK", 2)])
    assert _support([o], peps) == {"U": (True, True)}


def test_all_denovo_peptides() -> None:
    o = _orf("D", ["AAAGK", "CCCGK", "DDDGK"])
    peps = peptides_frame([_dn(p, k) for k, p in enumerate(o.peptides, 1)])
    assert _support([o], peps) == {"D": (True, False)}


def test_decoy_rows_are_ignored() -> None:
    o = _orf("D", ["AAAGK"])
    dec = pep_row("AAAGK", ["rev_x"], 3, is_decoy=True, source="database")
    assert _support([o], peptides_frame([_dn("AAAGK"), dec])) == {"D": (True, False)}
    assert _support([o], peptides_frame([dec])) == {"D": (False, False)}


def test_contaminant_rows_are_ignored() -> None:
    o = _orf("D", ["AAAGK"])
    con = {**_db("AAAGK", 3), "is_contaminant": True}
    assert _support([o], peptides_frame([_dn("AAAGK"), con])) == {"D": (True, False)}
    assert _support([o], peptides_frame([con])) == {"D": (False, False)}


def test_orf_peptide_absent_from_c9() -> None:
    a = _orf("A", ["WWWGK"])
    b = _orf("B", ["AAAGK", "WWWGK"], stop=2400)
    out = _support([a, b], peptides_frame([_dn("AAAGK")]))
    assert out == {"A": (False, False), "B": (True, False)}


def test_only_the_orfs_own_peptides_count() -> None:
    """Example 11's db half: the extended gene's in-frame peptides are not ORF peptides."""
    o = _orf("N", ["AAAGK", "CCCGK"], cls="nterm_extension")
    rows = [_dn("AAAGK"), _dn("CCCGK", 2), *(_db(f"GENE{c}K", 3 + i) for i, c in
                                              enumerate("ACDEFGHKM"))]  # fmt: skip
    assert _support([o], peptides_frame(rows)) == {"N": (True, False)}


def test_every_orf_gets_a_pair() -> None:
    orfs = [_orf(f"O{i}", [f"PEP{c}K"], stop=300 * (i + 1)) for i, c in enumerate("ACDEF")]
    out = _support(orfs, peptides_frame([_dn("PEPAK"), _db("PEPCK")]))
    assert set(out) == {o.orf_id for o in orfs}


# Hypothesis: db_supported equals a direct oracle and does not depend on C9 row order.

_RES = st.sampled_from(list("AGILK"))
_PEP = st.lists(_RES, min_size=3, max_size=4).map(lambda x: "".join(x) + "R")
_ROW = st.tuples(_PEP, st.sampled_from(["database", "denovo"]), st.booleans(), st.booleans())


def _fold(p: str) -> str:
    return p.upper().replace("I", "L")


def _oracle(orfs: Sequence[Any], rows: Sequence[tuple[str, str, bool, bool]]) -> dict:
    src: dict[str, set[str]] = {}
    for p, s, dec, con in rows:
        if not dec and not con:
            src.setdefault(_fold(p), set()).add(s)
    out = {}
    for o in orfs:
        s = set().union(*(src.get(p, set()) for p in o.peptides))
        out[o.orf_id] = ("denovo" in s, "database" in s)
    return out


def _frame(rows: Sequence[tuple[str, str, bool, bool]]) -> pd.DataFrame:
    recs = []
    for k, (p, s, dec, con) in enumerate(rows, 1):
        r = pep_row(p, ["rev_x" if dec else "x"], k, is_decoy=dec, source=s)
        recs.append({**r, "is_contaminant": con})
    return peptides_frame(recs)


@settings(max_examples=60, deadline=None)
@given(
    rows=st.lists(_ROW, min_size=1, max_size=12),
    orf_peps=st.lists(st.lists(_PEP, min_size=1, max_size=4), min_size=1, max_size=5),
    data=st.data(),
)
def test_db_supported_matches_oracle_and_ignores_row_order(
    rows: list[tuple[str, str, bool, bool]], orf_peps: list[list[str]], data: Any
) -> None:
    orfs = [_orf(f"O{i}", ps, stop=300 * (i + 1)) for i, ps in enumerate(orf_peps)]
    want = _oracle(orfs, rows)
    assert _support(orfs, _frame(rows)) == want
    perm = data.draw(st.permutations(rows))
    assert _support(orfs, _frame(perm)) == want


# ---------------------------------------------------------------- apply_evidence


Support = Mapping[str, tuple[bool, bool]] | str | None


def _evidence(case: Case, orfs: list[Any], *, best: Mapping[str, int | None],
              tir: Mapping[str, float], gt: pd.DataFrame = GT30,
              support: Support = "absent") -> pd.DataFrame:  # fmt: skip
    cols: dict[str, dict[str, Any]] = {"ostir_tir": dict(tir)}
    t = table(case, orfs, scored(orfs, best=dict(best), **cols))
    kw: dict[str, Any] = {} if support == "absent" else {"denovo_support": support}
    return ev().apply_evidence(t, orfs, gene_table=gt, terminators=[], gm=case.gm,
                               genome=case.genome, cfg=cfg(), peptide_intervals=None,
                               **kw)  # fmt: skip


def _pair(case: Case, *, best: Mapping[str, int | None], tir: Mapping[str, float],
          support: Mapping[str, tuple[bool, bool]],
          gt: pd.DataFrame = GT30) -> tuple[dict[str, Any], dict[str, Any]]:  # fmt: skip
    """(row without ``denovo_support``, row with it) of the single ORF of ``case``."""
    orfs = orfs_of(case)
    assert len(orfs) == 1, orfs
    (a,) = _evidence(case, orfs, best=best, tir=tir, gt=gt).to_dict("records")
    (b,) = _evidence(case, orfs, best=best, tir=tir, gt=gt, support=support).to_dict("records")
    return a, b


def _same_but_flag(a: Mapping[str, Any], b: Mapping[str, Any]) -> None:
    """Every column but ``flags`` equal; ``flags`` equal up to the new token."""
    for k in a:
        if k == "flags":
            continue
        va, vb = a[k], b[k]
        assert (is_na(va) and is_na(vb)) or va == vb, (k, va, vb)
    assert [f for f in flags_of(b) if f != FLAG] == flags_of(a)


def _ordered(row: Mapping[str, Any]) -> None:
    from gippa.novel.outputs import FLAG_ORDER

    f = flags_of(row)
    assert f == sorted(f, key=list(FLAG_ORDER).index), f


DN_ONLY = (True, False)


def test_example_1_intergenic_with_rbs_is_flagged() -> None:
    a, b = _pair(ex1_case(), best={ID1: 60}, tir={ID1: 500.0}, support={ID1: DN_ONLY})
    assert b["rbs_supported"] == True  # noqa: E712
    f = flags_of(b)
    assert "rbs_supported" in f and f[-1] == FLAG, f
    _same_but_flag(a, b)
    _ordered(b)


def test_example_2_cds_outframe_with_rbs_is_exempt() -> None:
    case = ex1_case("cds_outframe")
    a, b = _pair(case, best={ID1: 60}, tir={ID1: 500.0}, support={ID1: DN_ONLY})
    assert b["rbs_supported"] == True and b["called"] == True  # noqa: E712
    assert FLAG not in flags_of(b)
    _same_but_flag(a, b)


def test_example_3_cds_outframe_below_rbs_threshold_is_flagged() -> None:
    case = ex1_case("cds_outframe")
    a, b = _pair(case, best={ID1: 60}, tir={ID1: 120.0}, support={ID1: DN_ONLY})
    assert b["rbs_supported"] == False and b["called"] == True  # noqa: E712
    assert flags_of(b)[-1] == FLAG
    _same_but_flag(a, b)


def test_example_4_no_start_codon_rbs_na_is_flagged() -> None:
    seq = ex1_seq(atg=False, gtg=False).text()
    case = Case({"c1": seq}, dict(C1), ex1_hits(seq, cls="cds_outframe"))
    (orf,) = orfs_of(case)
    a, b = _pair(case, best={orf.orf_id: None}, tir={}, support={orf.orf_id: DN_ONLY})
    assert is_na(b["rbs_supported"])
    f = flags_of(b)
    assert "no_start_codon" in f and f[-1] == FLAG, f
    assert b["called"] == True  # noqa: E712
    _same_but_flag(a, b)
    _ordered(b)


def test_threshold_unestimated_is_no_exemption() -> None:
    case = ex1_case("cds_outframe")
    a, b = _pair(case, best={ID1: 60}, tir={ID1: 500.0}, gt=GT15, support={ID1: DN_ONLY})
    assert is_na(b["rbs_supported"])
    assert flags_of(b)[-1] == FLAG
    _same_but_flag(a, b)


def test_ostir_missing_is_no_exemption() -> None:
    case = ex1_case("cds_outframe")
    a, b = _pair(case, best={ID1: 60}, tir={}, support={ID1: DN_ONLY})
    assert is_na(b["rbs_supported"])
    assert flags_of(b)[-1] == FLAG
    _same_but_flag(a, b)


def test_example_7_uncalled_row_is_flagged() -> None:
    """One peptide: the cds_outframe class rule fails, so ``called`` is false; flag set anyway."""
    case = ex1_case("cds_outframe", n_hits=1)
    a, b = _pair(case, best={ID1: 60}, tir={ID1: 50.0}, support={ID1: DN_ONLY})
    assert b["called"] == False  # noqa: E712
    assert flags_of(b)[-1] == FLAG
    _same_but_flag(a, b)


def test_example_11_nterm_extension_with_rbs_is_flagged() -> None:
    a, b = _pair(nterm2(), best={ID6: 249}, tir={ID6: 700.0}, support={ID6: DN_ONLY})
    assert b["rbs_supported"] == True and b["called"] == True  # noqa: E712
    assert flags_of(b)[-1] == FLAG
    _same_but_flag(a, b)


@pytest.mark.parametrize("cls", ["intergenic", "cds_outframe"])
@pytest.mark.parametrize("pair", [(True, True), (False, True), (False, False)])
@pytest.mark.parametrize("tir", [500.0, 120.0])
def test_database_support_or_no_denovo_never_flags(
    cls: str, pair: tuple[bool, bool], tir: float
) -> None:
    a, b = _pair(ex1_case(cls), best={ID1: 60}, tir={ID1: tir}, support={ID1: pair})
    assert FLAG not in flags_of(b)
    _same_but_flag(a, b)


def test_none_support_is_the_library_default() -> None:
    case = ex1_case("cds_outframe")
    orfs = orfs_of(case)
    a = _evidence(case, orfs, best={ID1: 60}, tir={ID1: 120.0})
    b = _evidence(case, orfs, best={ID1: 60}, tir={ID1: 120.0}, support=None)
    pd.testing.assert_frame_equal(a, b)
    assert FLAG not in flags_of(b.iloc[0])


def test_flag_does_not_change_rank_or_called_on_a_multi_orf_table() -> None:
    from helpers.novel_outputs_fixture import two_contigs

    case = two_contigs()
    orfs = orfs_of(case)
    assert len(orfs) >= 2
    best = {o.orf_id: o.starts[-1] for o in orfs}
    tir = {o.orf_id: 500.0 for o in orfs}
    sup = {o.orf_id: (DN_ONLY if k % 2 == 0 else (True, True)) for k, o in enumerate(orfs)}
    a = _evidence(case, orfs, best=best, tir=tir)
    b = _evidence(case, orfs, best=best, tir=tir, support=sup)
    assert list(a["orf_id"]) == list(b["orf_id"])
    by_id = {o.orf_id: o for o in orfs}  # the table is in output order, not orfs_of order
    assert set(by_id) == set(b["orf_id"])
    for ra, rb in zip(a.to_dict("records"), b.to_dict("records"), strict=True):
        o = by_id[rb["orf_id"]]
        _same_but_flag(ra, rb)
        assert (FLAG in flags_of(rb)) is ev().denovo_only_unsupported(
            o.cls, *sup[o.orf_id], None if is_na(rb["rbs_supported"]) else
            bool(rb["rbs_supported"]))  # fmt: skip


# ---------------------------------------------------------------- writer and schema


def _flags_field() -> dict[str, Any]:
    d = json.loads(SCHEMA.read_text(encoding="utf-8"))
    return next(f for f in d["fields"] if f["name"] == "flags")


def test_schema_pattern_lists_the_flag_last_in_both_alternations() -> None:
    pat = _flags_field()["constraints"]["pattern"]
    first = pat.split("(?:", 2)[1].split(")", 1)[0].split("|")
    second = pat.rsplit("(?:", 1)[1].split(")", 1)[0].split("|")
    assert tuple(first) == (*OLD_FLAGS, FLAG)
    assert tuple(second) == (*OLD_FLAGS, FLAG)


def test_schema_description_names_the_flag() -> None:
    desc = _flags_field()["description"]
    assert f"`{FLAG}`" in desc and "ADR-0050" in desc, desc
    assert "de novo" in desc and "database" in desc and "RBS" in desc, desc


def test_writer_output_carries_the_flag_and_validates(tmp_path: Path) -> None:
    from helpers import novel_steps_fixture as nf

    from gippa.io.novel_tables import read_novel_orfs, write_novel_orfs

    case = ex1_case("cds_outframe")
    orfs = orfs_of(case)
    t = _evidence(case, orfs, best={ID1: 60}, tir={ID1: 120.0}, support={ID1: DN_ONLY})
    p = tmp_path / "novel_orfs.tsv"
    write_novel_orfs(p, t)
    rows = nf.tsv_dicts(p)
    assert str(rows[0]["flags"]).split(";")[-1] == FLAG
    assert flags_of(read_novel_orfs(p).iloc[0])[-1] == FLAG
    nf.assert_valid_novel_orfs(p)


def test_schema_rejects_a_misspelt_token(tmp_path: Path) -> None:
    frictionless = pytest.importorskip("frictionless")
    from gippa.io.novel_tables import write_novel_orfs

    case = ex1_case("cds_outframe")
    orfs = orfs_of(case)
    t = _evidence(case, orfs, best={ID1: 60}, tir={ID1: 120.0}, support={ID1: DN_ONLY})
    t["flags"] = t["flags"].str.replace(FLAG, "denovo_only_unsuported", regex=False)
    p = tmp_path / "novel_orfs.tsv"
    write_novel_orfs(p, t)
    d = {k: v for k, v in json.loads(SCHEMA.read_text(encoding="utf-8")).items()
         if k != "coordinates"}  # fmt: skip
    res = frictionless.Resource(path=p.name, basepath=str(p.parent), format="tsv",
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    assert not res.validate().valid
