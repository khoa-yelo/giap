"""Per-step random streams (GIAP_SPEC.md §6.6)."""

from __future__ import annotations

import zlib

import numpy as np
import pytest

from giap.rng import for_step


def test_returns_numpy_generator() -> None:
    assert isinstance(for_step("peptide_mapping", 1), np.random.Generator)


def test_same_name_and_seed_give_identical_draws() -> None:
    a = for_step("stats", 20261014).random(100)
    b = for_step("stats", 20261014).random(100)
    np.testing.assert_array_equal(a, b)


def test_different_names_give_different_streams() -> None:
    a = for_step("stats", 7).random(50)
    b = for_step("bootstrap", 7).random(50)
    assert not np.array_equal(a, b)


def test_different_seeds_give_different_streams() -> None:
    a = for_step("stats", 7).random(50)
    b = for_step("stats", 8).random(50)
    assert not np.array_equal(a, b)


@pytest.mark.parametrize("name", ["stats", "permutation_test", "", "étape"])
@pytest.mark.parametrize("seed", [0, 1, 20261014, 2**40])
def test_matches_spec_formula_exactly(name: str, seed: int) -> None:
    def expected() -> np.random.Generator:
        return np.random.default_rng(np.random.SeedSequence([seed, zlib.crc32(name.encode())]))

    np.testing.assert_array_equal(for_step(name, seed).random(20), expected().random(20))
    np.testing.assert_array_equal(
        for_step(name, seed).integers(0, 2**62, 20), expected().integers(0, 2**62, 20)
    )


def test_adding_a_step_does_not_change_another_steps_draws() -> None:
    before = for_step("stats", 3).random(10)
    for_step("new_step", 3).random(1000)
    after = for_step("stats", 3).random(10)
    np.testing.assert_array_equal(before, after)
