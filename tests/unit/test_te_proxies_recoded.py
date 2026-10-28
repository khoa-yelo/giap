"""TE1's ``sec_site`` label and TE2's note and ``alt_stop_codon`` count (M7-16c; spec rules 24-25,
example 31, parked defaults Q25-Q27; test plan C10).

* ``te1_result``: a miss whose peptide has at least one ``peptide_hits`` row and every row
  carries ``sec_site`` is an exception of ``kind = "sec_site"`` (``locus_tag``, ``peptide``,
  ``hits`` as today); every other miss stays ``kind = "miss"`` (no row: ``hits = "unplaced"``).
  ``counts["sec_site"]`` counts those misses before the cap. ``misses``, ``exact_inframe``, the
  share and ``passed`` are unchanged. Without a ``flags`` column the rows are unflagged (Q25).
  ``per_peptide`` (default None, Q26): each ``miss`` / ``sec_site`` exception whose peptide is in
  it also gets ``placement_q``, ``m_target`` and ``accepted``.
* ``te2_result``: ``counts["alt_stop_codon"]`` = ALT rows flagged ``stop_codon`` at any
  confidence (Q27); ``TE2_NO_STOP_NOTE`` is rule 25's string; the pass rule is unchanged.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from gippa.model.coords import Interval

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
te: Any = importlib.import_module("validation.proxies.te_proxies")

GENES = {"G1": Interval("chr", 100, 400, "+")}
CIRC = {"chr": 3000}
SOURCES = {"pS": {"G1"}, "pT": {"G1"}, "pU": {"G1"}}
NOTE = (
    "no high stop_codon site: annotated stops and recoded Sec/Pyl codons are consumable "
    "(M7-16, ADR-0058 items 1-6), but no stop_codon ALT site reached high confidence (on Nissle "
    "the two fdnG Sec sites are candidate, one unique peptide each; fdhF has no tryptic peptide "
    "over its Sec codon)"
)


def _hit(pep: str, route: str, s0: int, e0: int, flags: str | None = None) -> dict[str, Any]:
    row = {
        "peptide": pep,
        "route": route,
        "locus_class": "cds_inframe",
        "gene": "G1",
        "contig": "chr",
        "start": str(s0 + 1),
        "end": str(e0),
        "strand": "+",
    }
    if flags is not None:
        row["flags"] = flags
    return row


def example_31(with_flags: bool = True) -> pd.DataFrame:
    rows = [
        _hit("pS", "tolerant", 103, 148, "sec_site" if with_flags else None),
        _hit("pS", "tolerant", 203, 248, "multi;sec_site" if with_flags else None),
        _hit("pT", "tolerant", 106, 151, "" if with_flags else None),
    ]
    return pd.DataFrame(rows)


def _per(rows: list[tuple[str, int, float, bool]]) -> pd.DataFrame:
    from gippa.mapping.placement_fdr import PER_PEPTIDE_COLUMNS

    return pd.DataFrame(
        [
            {
                "peptide": p,
                "source": "database",
                "decoy_type": "reverse",
                "length": 28,
                "m_target": m,
                "e_target": 1e-6,
                "placement_q": q,
                "accepted": acc,
                "flags": "",
            }
            for p, m, q, acc in rows
        ],
        columns=list(PER_PEPTIDE_COLUMNS),
    )


# ---------------------------------------------------------------- TE1, example 31


def test_example_31_counts_and_labels() -> None:
    r = te.te1_result(SOURCES, {}, example_31(), GENES, CIRC)
    assert r.counts["misses"] == 3 and r.counts["sec_site"] == 1
    assert r.counts["exact_inframe"] == 0
    by = {x["peptide"]: x for x in r.exceptions}
    assert by["pS"]["kind"] == "sec_site" and by["pS"]["locus_tag"] == "G1"
    assert by["pS"]["hits"] == "tolerant:cds_inframe:G1;tolerant:cds_inframe:G1"
    assert by["pT"]["kind"] == "miss" and by["pT"]["hits"] == "tolerant:cds_inframe:G1"
    assert by["pU"] == {"kind": "miss", "locus_tag": "G1", "peptide": "pU", "hits": "unplaced"}
    assert [x["peptide"] for x in r.exceptions] == ["pS", "pT", "pU"]


def test_example_31_share_and_passed_as_without_flags() -> None:
    a = te.te1_result(SOURCES, {}, example_31(True), GENES, CIRC)
    b = te.te1_result(SOURCES, {}, example_31(False), GENES, CIRC)
    assert a.shares == b.shares and a.passed == b.passed
    assert {k: v for k, v in a.counts.items() if k != "sec_site"} == {
        k: v for k, v in b.counts.items() if k != "sec_site"
    }


def test_without_a_flags_column_everything_stays_miss_and_sec_site_is_zero() -> None:
    """Q25: callers whose frame has no ``flags`` column keep today's result plus ``sec_site`` 0."""
    r = te.te1_result(SOURCES, {}, example_31(False), GENES, CIRC)
    assert r.counts == {
        "peptides": 3,
        "non_atg_nterm": 0,
        "exact_inframe": 0,
        "misses": 3,
        "exact_inframe_masked": 0,
        "sec_site": 0,
    }
    assert {x["kind"] for x in r.exceptions} == {"miss"}


def test_a_peptide_with_one_unflagged_row_is_a_miss() -> None:
    hits = pd.DataFrame(
        [_hit("pS", "tolerant", 103, 148, "sec_site"), _hit("pS", "tolerant", 203, 248, "multi")]
    )
    r = te.te1_result({"pS": {"G1"}}, {}, hits, GENES, CIRC)
    assert r.counts["sec_site"] == 0 and r.exceptions[0]["kind"] == "miss"


def test_pyl_site_and_readthrough_are_not_sec_site() -> None:
    hits = pd.DataFrame(
        [
            _hit("pS", "tolerant", 103, 148, "pyl_site"),
            _hit("pT", "tolerant", 106, 151, "readthrough"),
        ]
    )
    r = te.te1_result({"pS": {"G1"}, "pT": {"G1"}}, {}, hits, GENES, CIRC)
    assert r.counts["sec_site"] == 0
    assert {x["kind"] for x in r.exceptions} == {"miss"}


def test_an_exact_inframe_peptide_with_sec_site_rows_is_not_a_miss() -> None:
    hits = pd.DataFrame([_hit("pS", "exact", 103, 148, "sec_site")])
    r = te.te1_result({"pS": {"G1"}}, {}, hits, GENES, CIRC)
    assert r.counts["exact_inframe"] == 1 and r.counts["sec_site"] == 0
    assert r.exceptions == []


def test_sec_site_count_is_taken_before_the_exception_cap() -> None:
    sources = {f"p{i:03d}": {"G1"} for i in range(80)}
    hits = pd.DataFrame([_hit(p, "tolerant", 103, 148, "sec_site") for p in sources])
    r = te.te1_result(sources, {}, hits, GENES, CIRC)
    assert r.counts["misses"] == 80 and r.counts["sec_site"] == 80
    assert len(r.exceptions) == 50 and {x["kind"] for x in r.exceptions} == {"sec_site"}


def test_sec_site_list_flags_are_read() -> None:
    """The contract reader gives ``flags`` as lists; both forms are accepted."""
    hits = example_31()
    hits["flags"] = [["sec_site"], ["multi", "sec_site"], []]
    r = te.te1_result(SOURCES, {}, hits, GENES, CIRC)
    assert r.counts["sec_site"] == 1


# ---------------------------------------------------------------- TE1, per_peptide (Q26)


def test_per_peptide_q_values_on_miss_and_sec_site_rows() -> None:
    per = _per([("pS", 1, 0.0012, True), ("pT", 1, 0.004, True)])
    r = te.te1_result(SOURCES, {}, example_31(), GENES, CIRC, per_peptide=per)
    by = {x["peptide"]: x for x in r.exceptions}
    assert (by["pS"]["kind"], by["pS"]["placement_q"], by["pS"]["m_target"]) == (
        "sec_site",
        pytest.approx(0.0012),
        1,
    )
    assert bool(by["pS"]["accepted"]) is True
    assert by["pT"]["placement_q"] == pytest.approx(0.004) and int(by["pT"]["m_target"]) == 1
    assert by["pU"].get("placement_q") is None  # not in the placement table
    assert by["pU"]["hits"] == "unplaced"


def test_per_peptide_rejected_peptide_records_accepted_false() -> None:
    """Spec "Rejection is a finding": a rejected Sec peptide has no ``peptide_hits`` row, stays
    ``miss`` / ``unplaced`` and records ``accepted = False``."""
    per = _per([("pU", 1, 0.2, False)])
    r = te.te1_result(SOURCES, {}, example_31(), GENES, CIRC, per_peptide=per)
    u = {x["peptide"]: x for x in r.exceptions}["pU"]
    assert (u["kind"], u["hits"]) == ("miss", "unplaced")
    assert not bool(u["accepted"]) and u["placement_q"] == pytest.approx(0.2)


def test_per_peptide_none_keeps_the_rows_as_today() -> None:
    r = te.te1_result(SOURCES, {}, example_31(), GENES, CIRC, per_peptide=None)
    for x in r.exceptions:
        assert not {"placement_q", "m_target", "accepted"} & set(x)


def test_per_peptide_does_not_touch_non_atg_rows() -> None:
    per = _per([("pS", 1, 0.001, True), ("pX", 0, 0.0, True)])
    sources = {**SOURCES, "pX": {"G1"}}
    r = te.te1_result(sources, {"pX": "G1"}, example_31(), GENES, CIRC, per_peptide=per)
    x = [e for e in r.exceptions if e["kind"] == "non_atg_nterm"]
    assert x == [{"kind": "non_atg_nterm", "locus_tag": "G1", "peptide": "pX"}]


def test_te1_proxy_passes_the_placement_table() -> None:
    """``te1_proxy`` reads the run's per-peptide placement table and passes it (rule 24)."""
    import inspect

    src = inspect.getsource(te.te1_proxy)
    assert "per_peptide" in src
    assert "per_peptide" in inspect.signature(te.te1_result).parameters
    assert inspect.signature(te.te1_result).parameters["per_peptide"].default is None


# ---------------------------------------------------------------- TE2 (rule 25)


def _row(i: int, conf: str, flags: tuple[str, ...], ref: str = "N", alt: str = "D") -> dict:  # type: ignore[type-arg]
    from gippa.variants.vcf_records import VARIANTS_COLUMNS

    row: dict[str, Any] = dict.fromkeys(VARIANTS_COLUMNS, pd.NA)
    row.update(
        variant_id=f"VAR_{i}",
        contig="chr1",
        codon_start=3 * i,
        strand="+",
        gene=f"G{i}",
        aa_pos=10 + i,
        ref_codon="TGA",
        ref_aa=ref,
        alt_aa=alt,
        alt_type="substitution",
        call="ALT",
        confidence=conf,
        n_ref=0,
        n_alt=2 if conf == "high" else 1,
        source="database",
        flags=tuple(flags),
    )
    return row


def _sites(rows: list[dict]) -> pd.DataFrame:  # type: ignore[type-arg]
    from gippa.variants.vcf_records import VARIANTS_COLUMNS

    return pd.DataFrame(rows, columns=list(VARIANTS_COLUMNS))


def test_te2_note_is_rule_25_text() -> None:
    assert te.TE2_NO_STOP_NOTE == NOTE
    r = te.te2_result(_sites([_row(1, "high", ("snv", "start_codon"), "V", "M")]))
    assert r.notes == [NOTE]
    assert r.to_json()["notes"] == [NOTE]


def test_te2_alt_stop_codon_counts_every_confidence() -> None:
    rows = [
        _row(1, "candidate", ("stop_codon", "vcf_unrepresentable"), "*", "U"),
        _row(2, "candidate", ("stop_codon", "vcf_unrepresentable"), "*", "U"),
        _row(3, "high", ("snv", "start_codon"), "V", "M"),
        _row(4, "candidate", ("snv",)),
    ]
    r = te.te2_result(_sites(rows))
    assert r.counts == {
        "alt_sites": 4,
        "high": 1,
        "high_start_codon": 1,
        "high_stop_codon": 0,
        "violations": 0,
        "alt_stop_codon": 2,
    }
    assert r.passed is True
    assert r.notes == [NOTE]


def test_te2_alt_stop_codon_with_a_high_stop_site_and_no_note() -> None:
    rows = [_row(1, "high", ("snv", "stop_codon"), "*", "W"), _row(2, "candidate", ("snv",))]
    r = te.te2_result(_sites(rows))
    assert (r.counts["alt_stop_codon"], r.counts["high_stop_codon"]) == (1, 1)
    assert r.notes == [] and r.passed is True


def test_te2_empty_sites() -> None:
    r = te.te2_result(_sites([]))
    assert r.counts["alt_stop_codon"] == 0 and r.passed is True


def test_module_docstring_no_longer_says_unmapped() -> None:
    doc = " ".join((te.__doc__ or "").split())
    assert "stay unmapped; ADR-0058" not in doc
