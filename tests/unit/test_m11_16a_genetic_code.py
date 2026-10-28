"""M11-16a ``gippa.model.genetic_code``: AC 5 and the item 1 errors
(docs/specs/M11-16-translation-tables-4-25.md, Definition item 1).

Interface (spec item 1): ``SUPPORTED_TABLES``, ``codon_table``, ``stop_codons``,
``minus_stop_codons``, ``recoded_codons`` and ``table_name``; every function raises
``ValueError("unsupported translation table {t}; GIPPA supports 4, 11 and 25")`` for any other
table; results are cached; the module is pure and does not import Biopython.
"""

from __future__ import annotations

import ast
import itertools
import re
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from Bio.Data import CodonTable
from hypothesis import given
from hypothesis import strategies as st

ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "src" / "gippa" / "model" / "genetic_code.py"
TABLES = (4, 11, 25)
BAD = (1, 2, 12, 26, 0, -1, 3, 33)
TCAG = {
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    25: "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}
CODONS = ["".join(c) for c in itertools.product("TCAG", repeat=3)]
_COMP = str.maketrans("ACGT", "TGCA")


def gc() -> Any:
    import gippa.model.genetic_code as mod

    return mod


def _revcomp(codon: str) -> str:
    return codon.translate(_COMP)[::-1]


def _biopython(table: int) -> dict[str, str]:
    t = CodonTable.unambiguous_dna_by_id[table]
    return {**t.forward_table, **dict.fromkeys(t.stop_codons, "*")}


def test_supported_tables() -> None:
    assert gc().SUPPORTED_TABLES == (4, 11, 25)
    assert isinstance(gc().SUPPORTED_TABLES, tuple)


@pytest.mark.parametrize("table", TABLES)
def test_ac5_codon_table_equals_biopython(table: int) -> None:
    got = dict(gc().codon_table(table))
    assert got == _biopython(table)
    assert len(got) == 64
    assert all(re.fullmatch("[ACGT]{3}", k) for k in got)
    assert all(len(v) == 1 for v in got.values())


@pytest.mark.parametrize("table", TABLES)
def test_codon_table_is_the_ncbi_tcag_string(table: int) -> None:
    table_map = gc().codon_table(table)
    assert "".join(table_map[c] for c in CODONS) == TCAG[table]


def test_only_tga_differs_between_the_tables() -> None:
    t11, t4, t25 = (gc().codon_table(t) for t in (11, 4, 25))
    diff = {c for c in CODONS if not (t11[c] == t4[c] == t25[c])}
    assert diff == {"TGA"}
    assert (t11["TGA"], t4["TGA"], t25["TGA"]) == ("*", "W", "G")


@pytest.mark.parametrize(
    ("table", "stops", "minus"),
    [
        (11, {"TAA", "TAG", "TGA"}, {"TTA", "CTA", "TCA"}),
        (4, {"TAA", "TAG"}, {"TTA", "CTA"}),
        (25, {"TAA", "TAG"}, {"TTA", "CTA"}),
    ],
)
def test_stop_and_minus_stop_codons(table: int, stops: set[str], minus: set[str]) -> None:
    s = gc().stop_codons(table)
    m = gc().minus_stop_codons(table)
    assert isinstance(s, frozenset) and isinstance(m, frozenset)
    assert s == stops
    assert m == minus
    assert m == frozenset(_revcomp(c) for c in s)
    assert s == frozenset(c for c, aa in gc().codon_table(table).items() if aa == "*")


@pytest.mark.parametrize(
    ("table", "want"),
    [(11, {"U": "TGA", "O": "TAG"}), (4, {"O": "TAG"}), (25, {"O": "TAG"})],
)
def test_recoded_codons(table: int, want: dict[str, str]) -> None:
    got = dict(gc().recoded_codons(table))
    assert got == want
    assert set(got.values()) <= gc().stop_codons(table)


@pytest.mark.parametrize(
    ("table", "name"),
    [
        (11, "Bacterial, Archaeal and Plant Plastid"),
        (4, "Mycoplasma and Spiroplasma"),
        (25, "Candidate Division SR1 and Gracilibacteria"),
    ],
)
def test_table_name(table: int, name: str) -> None:
    assert gc().table_name(table) == name


FUNCS = ["codon_table", "stop_codons", "minus_stop_codons", "recoded_codons", "table_name"]


@pytest.mark.parametrize("func", FUNCS)
@pytest.mark.parametrize("table", BAD)
def test_item1_unsupported_table_raises(func: str, table: int) -> None:
    fn: Callable[[int], object] = getattr(gc(), func)
    want = f"unsupported translation table {table}; GIPPA supports 4, 11 and 25"
    with pytest.raises(ValueError, match=f"^{re.escape(want)}$"):
        fn(table)


@given(table=st.integers(-1000, 1000).filter(lambda t: t not in TABLES))
def test_item1_every_other_integer_raises(table: int) -> None:
    for func in FUNCS:
        with pytest.raises(ValueError, match="unsupported translation table"):
            getattr(gc(), func)(table)


@pytest.mark.parametrize("func", FUNCS)
def test_results_are_cached(func: str) -> None:
    fn = getattr(gc(), func)
    for table in TABLES:
        assert fn(table) is fn(table)


@pytest.mark.parametrize("func", FUNCS)
def test_functions_are_exported(func: str) -> None:
    assert func in gc().__all__
    assert "SUPPORTED_TABLES" in gc().__all__


def test_module_is_pure_and_does_not_import_biopython() -> None:
    """Spec item 1: built from the TCAG strings, not imported from Biopython at run time;
    no I/O."""
    tree = ast.parse(MODULE.read_text(encoding="utf-8"))
    roots: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            roots |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            roots.add(node.module.split(".")[0])
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            assert node.func.id != "open", "the table module must not do I/O"
    assert "Bio" not in roots
    assert roots <= {"__future__", "functools", "collections", "typing", "types", "itertools"}
