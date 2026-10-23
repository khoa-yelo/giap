"""M7-10 class-specific false-ALT rates from paired de novo and database results
(docs/specs/M7-10-class-error-rates.md rules 1-4, examples 1-6, every degenerate row; lead
decisions 2026-10-22; GIAP_SPEC.md §7.5.4, §8.1 row "ALT call is real", §7.5.5 ``CLASS_ERR``).

Interface pinned here (test-author reading of the spec; see the report to the lead):

* ``giap.variants.error_rates`` is pure. ``ClassRate(alt_type, mass_equivalent, n, errors,
  rate)`` is a frozen, slotted dataclass in that field order; ``rate = errors / n`` and ``None``
  when ``n == 0``.
* ``paired_observations(denovo, reference, alt_support, q_max=0.01)``:

  - ``denovo``: the de novo rows, columns ``run``, ``spectrum_id``, ``peptide`` (stripped);
  - ``reference``: the ``--reference-peptides`` rows, columns ``run``, ``spectrum_id``,
    ``peptide``, ``q_value`` (float) and ``is_decoy`` (bool);
  - ``alt_support``: one row per (supporting peptide, ALT site), columns ``variant_id``,
    ``peptide``, ``alt_type`` and ``mass_equivalent`` (bool);
  - extra columns (``source`` here, as C17 carries it) are ignored.

  It returns one row per (de novo row, class) with at least the columns ``run``, ``spectrum_id``,
  ``alt_type``, ``mass_equivalent`` and ``error`` (bool). A de novo row is an observation when
  its peptide supports at least one ALT and its ``(run, spectrum_id)`` equals that of a target
  reference row with ``q_value <= q_max`` (rule 1; equality at the threshold pairs). Its classes
  are the distinct ``(alt_type, mass_equivalent)`` of the ALTs it supports: a row supporting two
  ALTs of one class is one observation of that class (rule 2 counts rows, not sites). ``error``
  is "the two peptides differ after I -> L folding" (rule 3). A de novo row without a
  ``spectrum_id`` never pairs.
* ``class_rates(obs)`` returns one ``ClassRate`` per class present in ``obs`` (order free); any
  extra entry it returns has ``n == 0`` and ``rate is None``.
* ``annotate_sites(sites, rates)`` returns a copy of ``sites`` (input untouched, same rows in
  the same order, every other column unchanged) with ``class_error_rate`` (the float rate, or the
  string ``unestimated``) and ``class_error_n`` (integer, NA when unestimated). A site's class is
  ``(alt_type, "mass_equivalent" in flags)`` with ``flags`` the M7-05 tuple column; frameshift
  rows are class ``(frameshift, False)`` (rule 4). A class absent from ``rates`` or with ``n == 0``
  is ``unestimated``. This is the ``sites`` column convention ``vcf_records`` already reads
  (M7-03), so example 6 goes through the merged ``vcf_records`` unchanged.
"""

from __future__ import annotations

import copy
import dataclasses
import importlib
import math
from collections.abc import Iterable, Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pytest
from helpers.variant_tables_fixture import GP, body_lines, info_dict, peptides_at, records, sub
from helpers.variant_tables_fixture import write as write_vcf_file

UNEST = "unestimated"
SUB_F = ("substitution", False)
SUB_T = ("substitution", True)
BLOCK_T = ("block", True)
SWAP_F = ("swap", False)
FS_F = ("frameshift", False)


def er() -> ModuleType:
    return importlib.import_module("giap.variants.error_rates")


# --------------------------------------------------------------------------- builders


def denovo(rows: Iterable[tuple[str, str | None, str]]) -> pd.DataFrame:
    """De novo rows ``(run, spectrum_id, peptide)``."""
    rows = list(rows)
    return pd.DataFrame(
        {
            "run": pd.Series([r[0] for r in rows], dtype=object),
            "spectrum_id": pd.Series([r[1] for r in rows], dtype=object),
            "peptide": pd.Series([r[2] for r in rows], dtype=object),
            "source": pd.Series(["denovo"] * len(rows), dtype=object),
        }
    )


def reference(rows: Iterable[tuple[Any, ...]]) -> pd.DataFrame:
    """Reference rows ``(run, spectrum_id, peptide[, q_value[, is_decoy]])``; q 0.001, target."""
    full = [(*r, 0.001, False)[:5] if len(r) == 3 else (*r, False)[:5] for r in rows]
    return pd.DataFrame(
        {
            "run": pd.Series([r[0] for r in full], dtype=object),
            "spectrum_id": pd.Series([r[1] for r in full], dtype=object),
            "peptide": pd.Series([r[2] for r in full], dtype=object),
            "source": pd.Series(["database"] * len(full), dtype=object),
            "q_value": pd.Series([float(r[3]) for r in full], dtype="float64"),
            "is_decoy": pd.Series([bool(r[4]) for r in full], dtype=bool),
        }
    )


def support(rows: Iterable[tuple[str, str, str, bool]]) -> pd.DataFrame:
    """ALT support ``(variant_id, peptide, alt_type, mass_equivalent)``."""
    rows = list(rows)
    return pd.DataFrame(
        {
            "variant_id": pd.Series([r[0] for r in rows], dtype=object),
            "peptide": pd.Series([r[1] for r in rows], dtype=object),
            "alt_type": pd.Series([r[2] for r in rows], dtype=object),
            "mass_equivalent": pd.Series([bool(r[3]) for r in rows], dtype=bool),
        }
    )


def sites(rows: Sequence[tuple[str, str, tuple[str, ...]]]) -> pd.DataFrame:
    """A ``sites``-like frame: ``(variant_id, alt_type, flags)`` plus plausible M7-02 columns."""
    alt = {"substitution": "S", "block": "GG", "swap": "TA", "frameshift": "fs+1"}
    return pd.DataFrame(
        {
            "variant_id": [r[0] for r in rows],
            "contig": ["chr1"] * len(rows),
            "codon_start": [30 * (i + 1) for i in range(len(rows))],
            "strand": ["+"] * len(rows),
            "gene": ["G1"] * len(rows),
            "aa_pos": pd.array([10 * (i + 1) for i in range(len(rows))], dtype="Int64"),
            "ref_codon": ["GCT"] * len(rows),
            "ref_aa": ["A"] * len(rows),
            "alt_aa": [alt[r[1]] for r in rows],
            "alt_type": [r[1] for r in rows],
            "call": ["ALT" if i % 2 == 0 else "mixed" for i in range(len(rows))],
            "confidence": ["high"] * len(rows),
            "n_ref": [i % 2 for i in range(len(rows))],
            "n_alt": [2] * len(rows),
            "source": ["denovo"] * len(rows),
            "flags": [r[2] for r in rows],
        }
    )


def rate(alt_type: str, me: bool, n: int, errors: int) -> Any:
    return er().ClassRate(alt_type, me, n, errors, errors / n if n else None)


def by_class(rates: Iterable[Any]) -> dict[tuple[str, bool], tuple[int, int, float | None]]:
    """``{(alt_type, mass_equivalent): (n, errors, rate)}`` of the entries with ``n > 0``."""
    out: dict[tuple[str, bool], tuple[int, int, float | None]] = {}
    for r in rates:
        if r.n == 0:
            assert r.rate is None, r
            continue
        key = (str(r.alt_type), bool(r.mass_equivalent))
        assert key not in out, f"class {key} returned twice"
        out[key] = (int(r.n), int(r.errors), r.rate)
    return out


def rates_of(dn: pd.DataFrame, ref: pd.DataFrame, sup: pd.DataFrame, **kw: Any) -> dict:
    m = er()
    return by_class(m.class_rates(m.paired_observations(dn, ref, sup, **kw)))


def is_na(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


def cell(df: pd.DataFrame, vid: str, col: str) -> Any:
    rows = df[df["variant_id"] == vid]
    assert len(rows) == 1, f"{vid}: {len(rows)} rows"
    return rows.iloc[0][col]


def assert_rate(df: pd.DataFrame, vid: str, want: float, n: int) -> None:
    got = cell(df, vid, "class_error_rate")
    assert not isinstance(got, str), f"{vid}: {got!r}, expected {want}"
    assert float(got) == pytest.approx(want, abs=1e-12), vid
    assert int(cell(df, vid, "class_error_n")) == n, vid


def assert_unestimated(df: pd.DataFrame, vid: str) -> None:
    assert cell(df, vid, "class_error_rate") == UNEST, vid
    assert is_na(cell(df, vid, "class_error_n")), vid


# --------------------------------------------------------------------------- example 1 data

# 10 spectra of run r1 whose de novo peptides support substitution sites VAR_1 / VAR_2 (not
# mass_equivalent); the reference agrees on 8 (two of them only up to I/L) and differs on 2.
EX1_DN = [("r1", f"s{i}", p) for i, p in enumerate(
    ["AAASAAK", "AAASAAK", "CCCSCCK", "CCCSCCK", "DDDSDDK", "DDDSDDK", "EEESEEK", "EEESEEK",
     "FFLSFFK", "FFLSFFK"], start=1)]  # fmt: skip
EX1_REF = [
    ("r1", "s1", "AAASAAK"), ("r1", "s2", "AAAAAAK"), ("r1", "s3", "CCCSCCK"),
    ("r1", "s4", "CCCSCCK"), ("r1", "s5", "DDDSDDK"), ("r1", "s6", "DDDADDK"),
    ("r1", "s7", "EEESEEK"), ("r1", "s8", "EEESEEK"), ("r1", "s9", "FFISFFK"),
    ("r1", "s10", "FFLSFFK"),
]  # fmt: skip
EX1_SUP = [
    ("VAR_1", "AAASAAK", "substitution", False), ("VAR_1", "CCCSCCK", "substitution", False),
    ("VAR_2", "DDDSDDK", "substitution", False), ("VAR_2", "EEESEEK", "substitution", False),
    ("VAR_2", "FFLSFFK", "substitution", False),
]  # fmt: skip


def ex1() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    return denovo(EX1_DN), reference(EX1_REF), support(EX1_SUP)


# --------------------------------------------------------------------------- the interface


def test_class_rate_is_a_frozen_slotted_dataclass() -> None:
    m = er()
    assert dataclasses.is_dataclass(m.ClassRate)
    names = [f.name for f in dataclasses.fields(m.ClassRate)]
    assert names == ["alt_type", "mass_equivalent", "n", "errors", "rate"]
    r = m.ClassRate("substitution", False, 10, 2, 0.2)
    with pytest.raises(dataclasses.FrozenInstanceError):
        r.n = 3  # type: ignore[misc]
    assert not hasattr(r, "__dict__"), "ClassRate should use slots"


def test_paired_observations_default_q_max_is_0_01() -> None:
    import inspect

    sig = inspect.signature(er().paired_observations)
    assert list(sig.parameters)[:3] == ["denovo", "reference", "alt_support"]
    assert sig.parameters["q_max"].default == 0.01


def test_observation_columns() -> None:
    obs = er().paired_observations(*ex1())
    assert {"run", "spectrum_id", "alt_type", "mass_equivalent", "error"} <= set(obs.columns)
    assert len(obs) == 10
    assert obs["error"].map(bool).sum() == 2


# --------------------------------------------------------------------------- examples


def test_example1_ten_substitutions_two_disagree() -> None:
    got = rates_of(*ex1())
    assert set(got) == {SUB_F}
    n, errors, r = got[SUB_F]
    assert (n, errors) == (10, 2)
    assert r == pytest.approx(0.2)


def test_example1_rate_on_every_row_of_the_class() -> None:
    m = er()
    rates = m.class_rates(m.paired_observations(*ex1()))
    st = sites([("VAR_1", "substitution", ("snv",)), ("VAR_2", "substitution", ("multi_nt",))])
    out = m.annotate_sites(st, rates)
    for vid in ("VAR_1", "VAR_2"):
        assert_rate(out, vid, 0.2, 10)


def test_example2_four_blocks_all_disagree() -> None:
    dn = denovo([("r1", f"b{i}", f"GG{c}PEPK") for i, c in enumerate("ACDE", start=1)])
    ref = reference([("r1", f"b{i}", f"N{c}PEPK") for i, c in enumerate("ACDE", start=1)])
    sup = support([("VAR_7", f"GG{c}PEPK", "block", True) for c in "ACDE"])
    m = er()
    rates = m.class_rates(m.paired_observations(dn, ref, sup))
    assert by_class(rates) == {BLOCK_T: (4, 4, pytest.approx(1.0))}
    out = m.annotate_sites(sites([("VAR_7", "block", ("mass_equivalent",))]), rates)
    assert_rate(out, "VAR_7", 1.0, 4)


def test_example3_il_difference_is_not_an_error() -> None:
    dn = denovo([("r1", "s1", "PEPTIDEK")])
    ref = reference([("r1", "s1", "PEPTLDEK")])
    sup = support([("VAR_1", "PEPTIDEK", "substitution", False)])
    obs = er().paired_observations(dn, ref, sup)
    assert len(obs) == 1
    assert not bool(obs["error"].iloc[0])
    assert rates_of(dn, ref, sup) == {SUB_F: (1, 0, pytest.approx(0.0))}


def test_example4_class_without_paired_observation_is_unestimated() -> None:
    m = er()
    rates = m.class_rates(m.paired_observations(*ex1()))
    st = sites([("VAR_1", "substitution", ()), ("VAR_3", "swap", ()),
                ("VAR_4", "substitution", ("mass_equivalent", "snv"))])  # fmt: skip
    out = m.annotate_sites(st, rates)
    assert_rate(out, "VAR_1", 0.2, 10)
    assert_unestimated(out, "VAR_3")  # (swap, False): no observation
    assert_unestimated(out, "VAR_4")  # (substitution, True): another class than VAR_1's


def test_example5_reference_q_0_02_is_ignored() -> None:
    dn = denovo([("r1", "s1", "AAASAAK"), ("r1", "s2", "AAASAAK")])
    sup = support([("VAR_1", "AAASAAK", "substitution", False)])
    m = er()
    ref = reference([("r1", "s1", "AAAAAAK", 0.02), ("r1", "s2", "AAAAAAK", 0.02)])
    assert len(m.paired_observations(dn, ref, sup)) == 0
    rates = m.class_rates(m.paired_observations(dn, ref, sup))
    assert by_class(rates) == {}
    out = m.annotate_sites(sites([("VAR_1", "substitution", ())]), rates)
    assert_unestimated(out, "VAR_1")
    # one confident pair left: the class is estimated from it alone
    ref2 = reference([("r1", "s1", "AAAAAAK", 0.02), ("r1", "s2", "AAAAAAK", 0.004)])
    assert rates_of(dn, ref2, sup) == {SUB_F: (1, 1, pytest.approx(1.0))}


def test_q_exactly_at_the_threshold_pairs() -> None:
    dn = denovo([("r1", "s1", "AAASAAK")])
    sup = support([("VAR_1", "AAASAAK", "substitution", False)])
    assert rates_of(dn, reference([("r1", "s1", "AAASAAK", 0.01)]), sup) == {
        SUB_F: (1, 0, pytest.approx(0.0))
    }


def test_q_max_is_honoured() -> None:
    dn = denovo([("r1", "s1", "AAASAAK")])
    sup = support([("VAR_1", "AAASAAK", "substitution", False)])
    ref = reference([("r1", "s1", "AAAAAAK", 0.03)])
    assert rates_of(dn, ref, sup) == {}
    assert rates_of(dn, ref, sup, q_max=0.05) == {SUB_F: (1, 1, pytest.approx(1.0))}


def test_example6_vcf_class_err_is_0_2(tmp_path: Path) -> None:
    """Example 6 through the merged M7-03 ``vcf_records`` and ``write_vcf``."""
    m = er()
    rates = m.class_rates(m.paired_observations(*ex1()))
    recs0, _table0, pile, st = records(peptides_at(GP, 44, sub(GP, 44, "W"), 2))
    assert len(st) == 1 and st["alt_type"].iloc[0] == "substitution"
    st["flags"] = [()]  # (substitution, False): example 1's class
    out = m.annotate_sites(st, rates)
    assert float(out["class_error_rate"].iloc[0]) == pytest.approx(0.2)
    from helpers.variant_tables_fixture import GENOME, REPS, vr

    recs, table = vr().vcf_records(out, pile, GENOME, REPS)
    assert table["class_error_rate"].tolist() == ["0.2"]
    assert table["class_error_n"].tolist() == [10]
    write_vcf_file(tmp_path / "v.vcf", recs)
    info = info_dict(body_lines(tmp_path / "v.vcf")[0][7])
    assert info["CLASS_ERR"] == "0.2"
    assert len(recs0) == 1


# --------------------------------------------------------------------------- degenerate rows


def test_no_reference_peptides_every_alt_row_unestimated() -> None:
    m = er()
    st = sites([("VAR_1", "substitution", ()), ("VAR_2", "block", ("mass_equivalent",)),
                ("VAR_3", "frameshift", ())])  # fmt: skip
    out = m.annotate_sites(st, [])
    for vid in ("VAR_1", "VAR_2", "VAR_3"):
        assert_unestimated(out, vid)
    dn, _ref, sup = ex1()
    empty_ref = reference([])
    assert len(m.paired_observations(dn, empty_ref, sup)) == 0
    assert by_class(m.class_rates(m.paired_observations(dn, empty_ref, sup))) == {}


def test_reference_rows_without_spectrum_id_never_pair() -> None:
    dn, _ref, sup = ex1()
    ref = reference([("r1", None, p) for _r, _s, p in EX1_REF])
    assert rates_of(dn, ref, sup) == {}


def test_denovo_rows_without_spectrum_id_never_pair() -> None:
    dn = denovo([("r1", None, "AAASAAK"), ("r1", "s2", "AAASAAK")])
    ref = reference([("r1", "s2", "AAAAAAK")])
    sup = support([("VAR_1", "AAASAAK", "substitution", False)])
    assert rates_of(dn, ref, sup) == {SUB_F: (1, 1, pytest.approx(1.0))}


def test_database_only_input_is_unestimated() -> None:
    """No de novo rows: §7.5.4 measures de novo errors only."""
    _dn, ref, sup = ex1()
    m = er()
    obs = m.paired_observations(denovo([]), ref, sup)
    assert len(obs) == 0
    rates = m.class_rates(obs)
    assert by_class(rates) == {}
    out = m.annotate_sites(sites([("VAR_1", "substitution", ())]), rates)
    assert_unestimated(out, "VAR_1")


def test_class_with_pairs_but_no_error_is_rate_0_with_n() -> None:
    dn = denovo([("r1", f"s{i}", "AAASAAK") for i in range(3)])
    ref = reference([("r1", f"s{i}", "AAASAAK") for i in range(3)])
    sup = support([("VAR_1", "AAASAAK", "substitution", False)])
    m = er()
    rates = m.class_rates(m.paired_observations(dn, ref, sup))
    assert by_class(rates) == {SUB_F: (3, 0, pytest.approx(0.0))}
    out = m.annotate_sites(sites([("VAR_1", "substitution", ())]), rates)
    got = cell(out, "VAR_1", "class_error_rate")
    assert got != UNEST and float(got) == 0.0
    assert int(cell(out, "VAR_1", "class_error_n")) == 3


def test_decoy_reference_row_is_not_a_pair() -> None:
    dn = denovo([("r1", "s1", "AAASAAK")])
    sup = support([("VAR_1", "AAASAAK", "substitution", False)])
    assert rates_of(dn, reference([("r1", "s1", "KAASAAA", 0.001, True)]), sup) == {}


def test_pairing_needs_the_same_run() -> None:
    dn = denovo([("r1", "s1", "AAASAAK")])
    sup = support([("VAR_1", "AAASAAK", "substitution", False)])
    assert rates_of(dn, reference([("r2", "s1", "AAAAAAK")]), sup) == {}
    assert rates_of(dn, reference([("r1", "s1", "AAAAAAK")]), sup) == {
        SUB_F: (1, 1, pytest.approx(1.0))
    }


def test_peptide_supporting_no_alt_is_not_an_observation() -> None:
    dn = denovo([("r1", "s1", "AAASAAK"), ("r1", "s2", "CCCCCCK")])
    ref = reference([("r1", "s1", "AAASAAK"), ("r1", "s2", "CCCACCK")])
    sup = support([("VAR_1", "AAASAAK", "substitution", False)])
    assert rates_of(dn, ref, sup) == {SUB_F: (1, 0, pytest.approx(0.0))}


# --------------------------------------------------------------------------- rule 2 classes


def test_peptide_supporting_two_classes_counts_in_both() -> None:
    dn = denovo([("r1", "s1", "GGASAAK")])
    ref = reference([("r1", "s1", "NAAAAAK")])
    sup = support([("VAR_1", "GGASAAK", "block", True),
                   ("VAR_2", "GGASAAK", "substitution", False)])  # fmt: skip
    assert rates_of(dn, ref, sup) == {
        BLOCK_T: (1, 1, pytest.approx(1.0)),
        SUB_F: (1, 1, pytest.approx(1.0)),
    }


def test_two_alts_of_one_class_are_one_observation() -> None:
    """Rule 2 counts de novo rows: two substitution ALTs on one peptide are one observation."""
    dn = denovo([("r1", "s1", "AWASAAK")])
    ref = reference([("r1", "s1", "AAAAAAK")])
    sup = support([("VAR_1", "AWASAAK", "substitution", False),
                   ("VAR_2", "AWASAAK", "substitution", False)])  # fmt: skip
    assert rates_of(dn, ref, sup) == {SUB_F: (1, 1, pytest.approx(1.0))}


def test_one_peptide_on_two_spectra_is_two_observations() -> None:
    dn = denovo([("r1", "s1", "AAASAAK"), ("r2", "s1", "AAASAAK")])
    ref = reference([("r1", "s1", "AAASAAK"), ("r2", "s1", "AAAAAAK")])
    sup = support([("VAR_1", "AAASAAK", "substitution", False)])
    assert rates_of(dn, ref, sup) == {SUB_F: (2, 1, pytest.approx(0.5))}


def test_mass_equivalent_splits_the_substitution_class() -> None:
    dn = denovo([("r1", "s1", "AAASAAK"), ("r1", "s2", "CCCWCCK")])
    ref = reference([("r1", "s1", "AAAAAAK"), ("r1", "s2", "CCCWCCK")])
    sup = support([("VAR_1", "AAASAAK", "substitution", True),
                   ("VAR_2", "CCCWCCK", "substitution", False)])  # fmt: skip
    m = er()
    rates = m.class_rates(m.paired_observations(dn, ref, sup))
    assert by_class(rates) == {SUB_T: (1, 1, pytest.approx(1.0)), SUB_F: (1, 0, 0.0)}
    st = sites([("VAR_1", "substitution", ("mass_equivalent", "snv")),
                ("VAR_2", "substitution", ("snv",))])  # fmt: skip
    out = m.annotate_sites(st, rates)
    assert_rate(out, "VAR_1", 1.0, 1)
    assert_rate(out, "VAR_2", 0.0, 1)


def test_frameshift_rows_are_their_own_class() -> None:
    dn = denovo([("r1", "s1", "AAAKPPPR"), ("r1", "s2", "AAASAAK")])
    ref = reference([("r1", "s1", "AAAKPPLR"), ("r1", "s2", "AAASAAK")])
    sup = support([("VAR_5", "AAAKPPPR", "frameshift", False),
                   ("VAR_1", "AAASAAK", "substitution", False)])  # fmt: skip
    m = er()
    rates = m.class_rates(m.paired_observations(dn, ref, sup))
    assert by_class(rates) == {FS_F: (1, 1, pytest.approx(1.0)), SUB_F: (1, 0, 0.0)}
    st = sites([("VAR_1", "substitution", ()), ("VAR_5", "frameshift", ())])
    out = m.annotate_sites(st, rates)
    assert_rate(out, "VAR_5", 1.0, 1)
    assert_rate(out, "VAR_1", 0.0, 1)


def test_class_rates_on_hand_built_observations() -> None:
    obs = pd.DataFrame(
        {
            "run": ["r1"] * 5,
            "spectrum_id": ["s1", "s2", "s3", "s4", "s5"],
            "alt_type": ["swap", "swap", "swap", "block", "block"],
            "mass_equivalent": [False, False, False, True, True],
            "error": [True, False, False, True, True],
        }
    )
    got = by_class(er().class_rates(obs))
    assert got[SWAP_F][:2] == (3, 1) and got[SWAP_F][2] == pytest.approx(1 / 3)
    assert got[BLOCK_T] == (2, 2, pytest.approx(1.0))
    assert set(got) == {SWAP_F, BLOCK_T}


def test_class_rates_of_no_observation() -> None:
    obs = er().paired_observations(denovo([]), reference([]), support([]))
    assert by_class(er().class_rates(obs)) == {}


# --------------------------------------------------------------------------- annotate_sites


def test_annotate_sites_keeps_rows_columns_and_input() -> None:
    m = er()
    st = sites([("VAR_1", "substitution", ()), ("VAR_2", "swap", ()),
                ("VAR_3", "block", ("mass_equivalent",)), ("VAR_4", "frameshift", ())])  # fmt: skip
    before = copy.deepcopy(st)
    rates = [rate("substitution", False, 10, 2), rate("block", True, 4, 4),
             rate("swap", False, 0, 0)]  # fmt: skip
    out = m.annotate_sites(st, rates)
    pd.testing.assert_frame_equal(st, before)  # pure: the input is untouched
    assert out["variant_id"].tolist() == ["VAR_1", "VAR_2", "VAR_3", "VAR_4"]
    for col in st.columns:
        assert out[col].tolist() == st[col].tolist(), col
    assert {"class_error_rate", "class_error_n"} <= set(out.columns)
    assert_rate(out, "VAR_1", 0.2, 10)
    assert_unestimated(out, "VAR_2")  # n == 0
    assert_rate(out, "VAR_3", 1.0, 4)
    assert_unestimated(out, "VAR_4")  # no frameshift rate
    n = out["class_error_n"]
    assert pd.api.types.is_integer_dtype(n.dtype), n.dtype


def test_annotate_sites_on_an_empty_frame() -> None:
    out = er().annotate_sites(sites([]), [rate("substitution", False, 10, 2)])
    assert len(out) == 0
    assert {"class_error_rate", "class_error_n"} <= set(out.columns)


def test_rates_written_by_the_merged_writer(tmp_path: Path) -> None:
    """``annotate_sites`` output -> M7-03 ``vcf_records`` -> ``write_variants``: a rate is written
    ``%.6g``, an unestimated row ``unestimated`` with ``class_error_n`` ``NA``."""
    from helpers.variant_tables_fixture import GENOME, REPS, vr

    from giap.io.variant_tables import write_variants

    m = er()
    aligns = peptides_at(GP, 44, sub(GP, 44, "W"), 2) + peptides_at(
        GP, 50, ("block", (GP.codon(50),), "N", "GG"), 1
    )
    _recs, _t, pile, st = records(aligns)
    st["flags"] = [("mass_equivalent",) if t == "block" else () for t in st["alt_type"]]
    out = m.annotate_sites(st, [rate("substitution", False, 3, 1)])
    _recs2, table = vr().vcf_records(out, pile, GENOME, REPS)
    write_variants(tmp_path / "variants.tsv", table)
    got = pd.read_csv(tmp_path / "variants.tsv", sep="\t", dtype=str, keep_default_na=False)
    by_type = dict(zip(got["alt_type"], zip(got["class_error_rate"], got["class_error_n"],
                                            strict=True), strict=True))  # fmt: skip
    assert by_type["substitution"] == ("0.333333", "3")
    assert by_type["block"] == (UNEST, "NA")
