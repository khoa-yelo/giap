"""Coordinates, Interval and Replicon (issue M1-01; GIAP_SPEC.md §2.1 decision 16, §6.5).

Internal coordinates are 0-based half-open; TSV/GFF3/VCF use 1-based closed. On a circular
replicon an origin-crossing interval is stored unwrapped (``end > length``); linear contigs never
wrap. ``crosses_origin`` is computed from the replicon, not stored (lead default, question 1).
"""

from __future__ import annotations

import ast
import dataclasses
from pathlib import Path

import pytest

import giap.model
from giap.model import (
    Interval,
    Replicon,
    check_interval,
    from_one_based_closed,
    interval_sort_key,
    split_at_origin,
    to_one_based_closed,
    wrap_position,
)

L = Replicon("c1", 1000, "linear")
C = Replicon("c1", 1000, "circular")


# --- Examples 1-10 (one test each) -------------------------------------------------------------


def test_example_01_to_one_based_closed() -> None:
    assert to_one_based_closed(Interval("c1", 99, 200, "+")) == (100, 200)


def test_example_02_from_one_based_closed_linear() -> None:
    assert from_one_based_closed("c1", 100, 200, "+", L) == Interval("c1", 99, 200, "+")


def test_example_03_minus_strand_round_trip_and_length() -> None:
    iv = Interval("c1", 10, 19, "-")
    assert to_one_based_closed(iv) == (11, 19)
    assert iv.strand == "-"
    assert iv.length_bp == 9


def test_example_04_unwrapped_origin_crossing_input() -> None:
    iv = from_one_based_closed("c1", 951, 1050, "-", C)
    assert iv == Interval("c1", 950, 1050, "-")
    assert iv.crosses_origin(C) is True


def test_example_05_wrapped_input_is_unwrapped_on_circular() -> None:
    assert from_one_based_closed("c1", 951, 50, "+", C) == Interval("c1", 950, 1050, "+")


def test_example_06_split_at_origin() -> None:
    assert split_at_origin(Interval("c1", 950, 1050, "-"), C) == (
        Interval("c1", 950, 1000, "-"),
        Interval("c1", 0, 50, "-"),
    )


def test_example_07_linear_contigs_never_wrap() -> None:
    with pytest.raises(ValueError):
        check_interval(Interval("c1", 950, 1050, "+"), L)


def test_example_08_wrap_position() -> None:
    assert wrap_position(1003, C) == 3
    with pytest.raises(ValueError):
        wrap_position(1003, L)


def test_example_09_empty_interval_rejected() -> None:
    with pytest.raises(ValueError):
        Interval("c1", 5, 5)


def test_example_10_sort_by_fasta_order_then_start_end_strand() -> None:
    ivs = [Interval("c2", 0, 9, "+"), Interval("c1", 5, 9, "-"), Interval("c1", 5, 9, "+")]
    order = {"c2": 0, "c1": 1}
    got = sorted(ivs, key=lambda iv: interval_sort_key(iv, order))
    assert got == [Interval("c2", 0, 9, "+"), Interval("c1", 5, 9, "+"), Interval("c1", 5, 9, "-")]


# --- Interface details ------------------------------------------------------------------------


def test_defaults() -> None:
    assert Interval("c1", 0, 1).strand == "."
    assert Replicon("p1", 5000, "circular").kind == "unknown"


def test_dataclasses_are_frozen() -> None:
    iv = Interval("c1", 0, 10)
    with pytest.raises(dataclasses.FrozenInstanceError):
        iv.start = 1  # type: ignore[misc]
    rep = Replicon("c1", 10, "linear", "chromosome")
    with pytest.raises(dataclasses.FrozenInstanceError):
        rep.length = 5  # type: ignore[misc]


def test_crosses_origin_is_computed_not_stored() -> None:
    names = {f.name for f in dataclasses.fields(Interval)}
    assert names == {"contig", "start", "end", "strand"}


@pytest.mark.parametrize("length", [0, -1])
def test_replicon_length_must_be_positive(length: int) -> None:
    with pytest.raises(ValueError):
        Replicon("c1", length, "linear")


def test_replicon_length_one_is_valid() -> None:
    assert Replicon("c1", 1, "circular").length == 1


def test_coords_module_matches_package_reexport() -> None:
    from giap.model import coords

    for name in (
        "Interval",
        "Replicon",
        "check_interval",
        "split_at_origin",
        "to_one_based_closed",
        "from_one_based_closed",
        "wrap_position",
        "interval_sort_key",
    ):
        assert getattr(giap.model, name) is getattr(coords, name)


# --- Edge cases -------------------------------------------------------------------------------


def test_whole_circular_replicon_from_zero_does_not_cross() -> None:
    iv = Interval("c1", 0, 1000)
    check_interval(iv, C)
    assert iv.crosses_origin(C) is False
    assert split_at_origin(iv, C) == (iv,)


def test_end_equal_length_does_not_cross_but_length_plus_one_does() -> None:
    assert Interval("c1", 900, 1000).crosses_origin(C) is False
    assert Interval("c1", 900, 1001).crosses_origin(C) is True


def test_whole_replicon_unwrapped_with_positive_start_is_valid() -> None:
    iv = Interval("c1", 400, 1400, "+")
    check_interval(iv, C)
    assert iv.crosses_origin(C) is True
    assert split_at_origin(iv, C) == (Interval("c1", 400, 1000, "+"), Interval("c1", 0, 400, "+"))


def test_longer_than_circular_replicon_rejected() -> None:
    with pytest.raises(ValueError):
        check_interval(Interval("c1", 400, 1401), C)


def test_negative_start_rejected() -> None:
    with pytest.raises(ValueError):
        Interval("c1", -1, 10)


def test_end_before_start_rejected() -> None:
    with pytest.raises(ValueError):
        Interval("c1", 10, 5)


def test_contig_mismatch_rejected() -> None:
    with pytest.raises(ValueError):
        check_interval(Interval("c2", 0, 10), C)


@pytest.mark.parametrize("rep", [L, C], ids=["linear", "circular"])
def test_start_must_be_inside_replicon(rep: Replicon) -> None:
    with pytest.raises(ValueError):
        check_interval(Interval("c1", 1000, 1001), rep)


def test_linear_end_equal_length_is_valid() -> None:
    check_interval(Interval("c1", 0, 1000), L)


def test_split_on_linear_returns_input() -> None:
    iv = Interval("c1", 10, 1000, "-")
    check_interval(iv, L)
    assert split_at_origin(iv, L) == (iv,)


def test_split_order_is_same_for_plus_strand() -> None:
    assert split_at_origin(Interval("c1", 950, 1050, "+"), C) == (
        Interval("c1", 950, 1000, "+"),
        Interval("c1", 0, 50, "+"),
    )


def test_wrapped_input_on_linear_rejected() -> None:
    with pytest.raises(ValueError):
        from_one_based_closed("c1", 951, 50, "+", L)


def test_from_one_based_closed_single_base() -> None:
    assert from_one_based_closed("c1", 1, 1, ".", L) == Interval("c1", 0, 1, ".")


def test_wrap_position_inside_range_is_identity() -> None:
    assert wrap_position(0, C) == 0
    assert wrap_position(999, C) == 999
    assert wrap_position(999, L) == 999
    assert wrap_position(1000, C) == 0


def test_sort_strand_rank_plus_minus_dot() -> None:
    ivs = [Interval("c1", 0, 9, "."), Interval("c1", 0, 9, "-"), Interval("c1", 0, 9, "+")]
    got = sorted(ivs, key=lambda iv: interval_sort_key(iv, {"c1": 0}))
    assert [iv.strand for iv in got] == ["+", "-", "."]


def test_sort_by_start_then_end_then_id() -> None:
    order = {"c1": 0}
    keys = [
        interval_sort_key(Interval("c1", 5, 20, "+"), order, "a"),
        interval_sort_key(Interval("c1", 5, 10, "+"), order, "b"),
        interval_sort_key(Interval("c1", 5, 10, "+"), order, "a"),
        interval_sort_key(Interval("c1", 1, 30, "-"), order, "z"),
    ]
    assert sorted(range(4), key=lambda i: keys[i]) == [3, 2, 1, 0]


# --- Acceptance: no file I/O in giap.model ------------------------------------------------------

_IO_MODULES = {"io", "os", "pathlib", "shutil", "gzip", "bz2", "lzma", "tempfile"}


def test_model_package_does_no_file_io() -> None:
    model_dir = Path(giap.model.__file__).parent
    for path in sorted(model_dir.glob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                assert not {a.name.split(".")[0] for a in node.names} & _IO_MODULES, path
            elif isinstance(node, ast.ImportFrom) and node.module:
                assert node.module.split(".")[0] not in _IO_MODULES, path
            elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                assert node.func.id != "open", path
