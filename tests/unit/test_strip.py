"""Modification stripping for converters (issue M1-07; GIPPA_SPEC.md §5.4)."""

from __future__ import annotations

import re

import pytest
from hypothesis import given
from hypothesis import strategies as st

from gippa.io.converters.strip import strip_modifications

UPPER = re.compile(r"^[A-Z]*$")


@pytest.mark.parametrize(
    ("seq", "expected"),
    [
        # spec examples 1-4
        ("_PEPM(UniMod:35)IDEK_", "PEPMIDEK"),
        ("n[42]AC[+57.021]DEK", "ACDEK"),
        ("PEPS[UNIMOD:21]TIDE+15.995K", "PEPSTIDEK"),
        ("PEPmIDEK", "PEPIDEK"),
        # §5.4 list, one syntax each
        ("PEPM[+15.995]IDEK", "PEPMIDEK"),
        ("PEPM[UNIMOD:35]IDEK", "PEPMIDEK"),
        ("PEPM+15.995IDEK", "PEPMIDEK"),
        ("PEPM-17.03IDEK", "PEPMIDEK"),
        ("PEPM+16IDEK", "PEPMIDEK"),
        ("n[42]PEPTIDEK", "PEPTIDEK"),
        ("PEPTIDEKc[-0.984]", "PEPTIDEK"),
        ("K.PEPTIDEK.R", "PEPTIDEK"),  # flanks removed (M1-07 review)
        ("PEPTIDEK", "PEPTIDEK"),
        ("", ""),
        ("pepTIDEk", "TIDE"),
        ("AC(UniMod:4)M(UniMod:35)K", "ACMK"),
    ],
)
def test_strip_examples(seq: str, expected: str) -> None:
    assert strip_modifications(seq) == expected


@given(st.text(max_size=60))
def test_output_is_uppercase_letters_only(seq: str) -> None:
    assert UPPER.match(strip_modifications(seq))


@given(st.text(alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ", max_size=60))
def test_plain_uppercase_is_identity(seq: str) -> None:
    assert strip_modifications(seq) == seq


@given(st.text(max_size=60))
def test_idempotent(seq: str) -> None:
    once = strip_modifications(seq)
    assert strip_modifications(once) == once
