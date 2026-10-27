"""Residue masses, the generated mass-equivalence list and the difference-operation tables
(issue M3-05; GIPPA_SPEC.md Appendix B.5, §7.2.4, §7.5.3).

Contract choices made by the test author where the issue spec is silent (see the PR body):

- ``MassEquivalence`` sides may come in either order; residue compositions are compared as
  multisets (``"DA"`` == ``"AD"``), modified forms are written ``"<residue>+<Unimod name>"``.
- ``BlockTables`` exposes its numba arrays as ``numpy.uint8`` attributes and a decoded view
  ``pairs() -> Iterable[tuple[str, str, str]]`` of ``(left, right, kind)`` residue strings
  (M3-02 letters, ``I`` -> ``L``), kind in {``block``, ``near_block``}.
"""

from __future__ import annotations

import dataclasses
import os
import subprocess
import sys
import typing
from collections.abc import Mapping
from importlib.resources import files

import numpy as np
import pytest

from gippa.mapping import diffops, masses
from gippa.mapping.diffops import BLOCKS, DIFF_TYPES, NEAR_BLOCKS, block_tables
from gippa.mapping.masses import MOD_DELTA, RESIDUE_MASS, MassEquivalence, mass_equivalents

# Appendix B.5, verbatim.
B5_RESIDUES: dict[str, float] = {
    "G": 57.02146,
    "A": 71.03711,
    "S": 87.03203,
    "P": 97.05276,
    "V": 99.06841,
    "T": 101.04768,
    "C": 103.00919,
    "L": 113.08406,
    "I": 113.08406,
    "N": 114.04293,
    "D": 115.02694,
    "Q": 128.05858,
    "K": 128.09496,
    "E": 129.04259,
    "M": 131.04049,
    "H": 137.05891,
    "F": 147.06841,
    "R": 156.10111,
    "Y": 163.06333,
    "W": 186.07931,
    "U": 150.95364,
    "O": 237.14773,
}
B5_MODS: dict[str, float] = {
    "Oxidation": 15.994915,
    "Dioxidation": 31.989829,
    "Deamidated": 0.984016,
    "Methyl": 14.015650,
    "Dimethyl": 28.031300,
    "Trimethyl": 42.046950,
    "Acetyl": 42.010565,
    "Carbamidomethyl": 57.021464,
    "Carbamyl": 43.005814,
    "Formyl": 27.994915,
    "Phospho": 79.966331,
    "Amidated": -0.984016,
    "Dehydrated": -18.010565,
    "Ammonia-loss": -17.026549,
    "Cation:Na": 21.981943,
}
STANDARD = "ACDEFGHKLMNPQRSTVWY"  # 20 standard residues with I merged into L


def _canon(side: str) -> str:
    """Multiset form of a residue composition; modified forms are kept as written."""
    return side if "+" in side else "".join(sorted(side))


def _mass(side: str) -> float:
    if "+" in side:
        res, mod = side.split("+", 1)
        return RESIDUE_MASS[res] + MOD_DELTA[mod]
    return sum(RESIDUE_MASS[c] for c in side)


def _index(entries: list[MassEquivalence]) -> dict[frozenset[str], MassEquivalence]:
    out: dict[frozenset[str], MassEquivalence] = {}
    for e in entries:
        key = frozenset({_canon(e.left), _canon(e.right)})
        assert key not in out, f"duplicate entry for {sorted(key)}"
        out[key] = e
    return out


def _find(entries: list[MassEquivalence], a: str, b: str) -> MassEquivalence | None:
    return _index(entries).get(frozenset({_canon(a), _canon(b)}))


@pytest.fixture(scope="module")
def default_list() -> list[MassEquivalence]:
    return mass_equivalents()


# --------------------------------------------------------------------------- mass tables


def test_residue_masses_match_b5_exactly() -> None:
    assert isinstance(RESIDUE_MASS, Mapping)
    assert dict(RESIDUE_MASS) == B5_RESIDUES


def test_i_and_l_share_one_mass() -> None:
    assert RESIDUE_MASS["I"] == RESIDUE_MASS["L"] == 113.08406


def test_modification_deltas_match_b5_exactly() -> None:
    assert isinstance(MOD_DELTA, Mapping)
    assert dict(MOD_DELTA) == B5_MODS


@pytest.mark.parametrize(
    ("name", "source", "values"),
    [
        ("residue_masses.tsv", "B.5", B5_RESIDUES),
        ("modifications.tsv", "Unimod", B5_MODS),
    ],
)
def test_data_tsv_has_source_header_and_b5_values(
    name: str, source: str, values: dict[str, float]
) -> None:
    text = files("gippa.data").joinpath(name).read_text(encoding="utf-8")
    lines = text.splitlines()
    comments = [ln for ln in lines if ln.startswith("#")]
    assert comments, f"{name} must start with a header comment naming its source"
    assert lines[0].startswith("#")
    header = "\n".join(comments)
    assert "B.5" in header
    assert source in header
    rows = [ln.split("\t") for ln in lines if ln and not ln.startswith("#")]
    by_key = {r[0]: r for r in rows}
    for key, value in values.items():
        assert key in by_key, f"{name}: no row for {key}"
        assert any(_num(f) == value for f in by_key[key][1:]), key


def _num(text: str) -> float | None:
    try:
        return float(text.replace(chr(0x2212), "-"))  # tolerate a typographic minus
    except ValueError:
        return None


# --------------------------------------------------------------------------- examples 1-7


def test_example1_n_minus_gg() -> None:
    assert RESIDUE_MASS["N"] - 2 * RESIDUE_MASS["G"] == pytest.approx(0.00001, abs=1e-9)


@pytest.mark.parametrize(
    ("a", "b", "delta", "kind"),
    [
        # Appendix B.5 "Checks the generated list MUST pass" (Δ 0 means < 0.001 here)
        ("N", "GG", 0.00001, "block"),
        ("Q", "AG", 0.00001, "block"),
        ("GT", "AS", 0.0, "block"),
        ("AV", "GL", 0.0, "block"),
        ("AN", "GQ", 0.0, "block"),
        ("W", "AD", 0.01526, "near_block"),
        ("W", "EG", 0.01526, "near_block"),
        ("R", "GV", 0.01124, "near_block"),
        ("K", "Q", 0.03638, "substitution"),
        ("A+Oxidation", "S", 0.000005, "modification"),
        ("N+Deamidated", "D", 0.000006, "modification"),
        ("G+Methyl", "A", 0.0, "modification"),
        ("M+Oxidation", "F", 0.033005, "modification"),
    ],
)
def test_b5_checks_in_default_list(
    default_list: list[MassEquivalence], a: str, b: str, delta: float, kind: str
) -> None:
    e = _find(default_list, a, b)
    assert e is not None, f"{a} vs {b} missing from mass_equivalents()"
    assert e.kind == kind
    assert e.delta_da == pytest.approx(delta, abs=2e-6)
    assert e.delta_da == pytest.approx(abs(_mass(a) - _mass(b)), abs=1e-9)


def test_example2_reported_values_round_as_in_b5(default_list: list[MassEquivalence]) -> None:
    for a, b, reported in [("W", "AD", 0.0153), ("W", "EG", 0.0153), ("R", "GV", 0.0112)]:
        e = _find(default_list, a, b)
        assert e is not None
        assert round(e.delta_da, 4) == reported
    kq = _find(default_list, "K", "Q")
    assert kq is not None and round(kq.delta_da, 4) == 0.0364
    mf = _find(default_list, "M+Oxidation", "F")
    assert mf is not None and abs(mf.delta_da - 0.0329) <= 2e-4


def test_example6_tolerance_cuts_m_oxidation_vs_f() -> None:
    assert _find(mass_equivalents(0.05), "M+Oxidation", "F") is not None
    small = mass_equivalents(tol_da=0.03)
    assert _find(small, "M+Oxidation", "F") is None
    assert _find(small, "K", "Q") is None  # Δ 0.0364
    assert _find(small, "W", "AD") is not None  # Δ 0.0153
    assert _find(small, "N", "GG") is not None


@pytest.mark.parametrize(
    ("a", "b"),
    [
        # §7.5.3 expected hits of the generated list
        ("N+Deamidated", "D"),
        ("Q+Deamidated", "E"),
        ("A+Oxidation", "S"),
        ("F+Oxidation", "Y"),
        ("G+Methyl", "A"),
        ("S+Methyl", "T"),
        ("V+Methyl", "L"),
        ("D+Methyl", "E"),
        ("N+Methyl", "Q"),
    ],
)
def test_section_7_5_3_modification_hits(
    default_list: list[MassEquivalence], a: str, b: str
) -> None:
    e = _find(default_list, a, b)
    assert e is not None, f"{a} vs {b}"
    assert e.kind == "modification"


# --------------------------------------------------------------------------- list invariants


def test_i_l_merged_u_o_excluded(default_list: list[MassEquivalence]) -> None:
    for e in default_list:
        for side in (e.left, e.right):
            residues = side.split("+", 1)[0]
            assert set(residues) <= set(STANDARD), (e, "only standard residues, I written as L")
    assert _find(default_list, "I", "L") is None


def test_no_trivial_entries(default_list: list[MassEquivalence]) -> None:
    for e in default_list:
        assert _canon(e.left) != _canon(e.right), e


def test_entries_within_tolerance_and_delta_recomputed(
    default_list: list[MassEquivalence],
) -> None:
    assert default_list
    for e in default_list:
        assert 0.0 <= e.delta_da <= 0.05
        assert e.delta_da == pytest.approx(abs(_mass(e.left) - _mass(e.right)), abs=1e-9)


def test_kind_matches_shape(default_list: list[MassEquivalence]) -> None:
    for e in default_list:
        if "+" in e.left or "+" in e.right:
            assert e.kind == "modification", e
            assert ("+" in e.left) != ("+" in e.right), e
            plain = e.right if "+" in e.left else e.left
            assert len(plain) == 1
        elif len(e.left) == len(e.right) == 1:
            assert e.kind == "substitution", e
        else:
            assert {len(e.left), len(e.right)} in ({1, 2}, {2}), e
            assert e.kind in ("block", "near_block"), e


def test_sorted_by_kind_delta_left_right(default_list: list[MassEquivalence]) -> None:
    keys = [(e.kind, e.delta_da, e.left, e.right) for e in default_list]
    assert keys == sorted(keys)


def test_deterministic_within_process() -> None:
    assert mass_equivalents() == mass_equivalents()
    assert mass_equivalents(0.02) == mass_equivalents(0.02)


def test_deterministic_across_hash_seeds() -> None:
    code = (
        "from gippa.mapping.masses import mass_equivalents as m;"
        "print(repr([(e.left, e.right, e.delta_da, e.kind) for e in m()]))"
    )
    outs = []
    for seed in ("0", "12345"):
        env = {**os.environ, "PYTHONHASHSEED": seed}
        r = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, env=env, check=True
        )
        outs.append(r.stdout)
    assert outs[0] == outs[1]
    assert outs[0].strip() not in ("", "[]")


def test_mass_equivalence_is_frozen_slotted_dataclass() -> None:
    assert dataclasses.is_dataclass(MassEquivalence)
    names = [f.name for f in dataclasses.fields(MassEquivalence)]
    assert names == ["left", "right", "delta_da", "kind"]
    e = MassEquivalence(left="K", right="Q", delta_da=0.03638, kind="substitution")
    with pytest.raises(dataclasses.FrozenInstanceError):
        e.left = "Q"  # type: ignore[misc]
    assert not hasattr(e, "__dict__")


# --------------------------------------------------------------------------- diffops


def test_diff_types_order() -> None:
    assert DIFF_TYPES == ("substitution", "swap", "block", "near_block", "frameshift")
    assert set(typing.get_args(diffops.DiffType)) == set(DIFF_TYPES)


def test_blocks_fixed_by_spec() -> None:
    assert BLOCKS == (("N", "GG"), ("Q", "AG"), ("GT", "AS"), ("AV", "GL"), ("AN", "GQ"))
    assert NEAR_BLOCKS == (("W", "AD"), ("W", "EG"), ("R", "GV"))
    for a, b in BLOCKS + NEAR_BLOCKS:
        assert a == "".join(sorted(a)) and b == "".join(sorted(b))


def test_blocks_consistent_with_generated_list(default_list: list[MassEquivalence]) -> None:
    for a, b in BLOCKS:
        e = _find(default_list, a, b)
        assert e is not None, (a, b)
        assert e.delta_da < 0.001
        assert e.kind == "block"
    for a, b in NEAR_BLOCKS:
        e = _find(default_list, a, b)
        assert e is not None, (a, b)
        assert e.delta_da <= 0.05
        assert e.kind == "near_block"


def _orders(comp: str) -> set[str]:
    return {comp, comp[::-1]}


def _expected_table_pairs() -> set[tuple[frozenset[str], str]]:
    out: set[tuple[frozenset[str], str]] = set()
    for pairs, kind in ((BLOCKS, "block"), (NEAR_BLOCKS, "near_block")):
        for a, b in pairs:
            for x in _orders(a):
                for y in _orders(b):
                    out.add((frozenset({x, y}), kind))
    return out


@pytest.fixture(scope="module")
def table_pairs() -> set[tuple[frozenset[str], str]]:
    tables = block_tables()
    got = {(frozenset({left, right}), kind) for left, right, kind in tables.pairs()}
    return got


def test_block_tables_expand_both_orders(table_pairs: set[tuple[frozenset[str], str]]) -> None:
    expected = _expected_table_pairs()
    assert len(expected) == 15 + 6
    assert table_pairs == expected


@pytest.mark.parametrize(
    ("a", "b", "kind"),
    [
        ("N", "GG", "block"),
        ("Q", "AG", "block"),
        ("Q", "GA", "block"),
        ("GT", "AS", "block"),
        ("TG", "SA", "block"),
        ("GT", "SA", "block"),
        ("AV", "GL", "block"),
        ("VA", "LG", "block"),
        ("AN", "QG", "block"),
        ("W", "DA", "near_block"),
    ],
)
def test_example8_block_tables_contents(
    table_pairs: set[tuple[frozenset[str], str]], a: str, b: str, kind: str
) -> None:
    assert (frozenset({a, b}), kind) in table_pairs


def test_example8_block_tables_exclusions(table_pairs: set[tuple[frozenset[str], str]]) -> None:
    flat = {s for pair, _ in table_pairs for s in pair}
    assert "I" not in "".join(flat)  # I -> L applied
    assert all(frozenset({"K", "Q"}) != pair for pair, _ in table_pairs)
    assert all(len(pair) == 2 for pair, _ in table_pairs)


def test_block_tables_masses_within_tolerance(
    table_pairs: set[tuple[frozenset[str], str]],
) -> None:
    for pair, kind in table_pairs:
        a, b = sorted(pair)
        d = abs(_mass(a) - _mass(b))
        assert d < (0.001 if kind == "block" else 0.05 + 1e-12), (a, b, kind)


def test_block_tables_arrays_are_uint8() -> None:
    tables = block_tables()
    arrays = [v for v in _attr_values(tables) if isinstance(v, np.ndarray)]
    assert arrays, "BlockTables must carry numpy arrays for the numba verifier"
    for arr in arrays:
        assert arr.dtype == np.uint8
        assert arr.flags.c_contiguous


def test_block_tables_deterministic() -> None:
    a, b = block_tables(), block_tables()
    assert sorted(a.pairs()) == sorted(b.pairs())
    for x, y in zip(_attr_values(a), _attr_values(b), strict=True):
        if isinstance(x, np.ndarray):
            np.testing.assert_array_equal(x, y)


def _attr_values(obj: object) -> list[object]:
    if dataclasses.is_dataclass(obj):
        return [getattr(obj, f.name) for f in dataclasses.fields(obj)]
    if hasattr(obj, "_fields"):
        return list(obj)  # type: ignore[call-overload]
    return list(vars(obj).values())


def test_modules_exist_with_public_names() -> None:
    for name in ("RESIDUE_MASS", "MOD_DELTA", "MassEquivalence", "mass_equivalents"):
        assert hasattr(masses, name)
    for name in ("DiffType", "DIFF_TYPES", "BLOCKS", "NEAR_BLOCKS", "block_tables"):
        assert hasattr(diffops, name)
