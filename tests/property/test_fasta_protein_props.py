"""Hypothesis properties for the protein FASTA reader (issue M1-03; GIPPA_SPEC.md §5.3)."""

from __future__ import annotations

import tempfile
from pathlib import Path

from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.io.protein_fasta import (
    STANDARD_AA,
    read_protein_fasta,
    resolve_protein_id,
    route_by_length,
)

FAST = settings(max_examples=150, deadline=None)
STD_SEQS = st.text(alphabet=sorted(STANDARD_AA), min_size=0, max_size=300)
ACC = st.text(alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_", min_size=1, max_size=12)
PART = st.one_of(
    ACC,
    st.builds(lambda db, a, n: f"{db}|{a}|{n}", st.sampled_from(["sp", "tr"]), ACC, ACC),
)
TOKENS = st.lists(PART, min_size=1, max_size=5).map(";".join)
BOUNDS = st.tuples(st.integers(1, 60), st.integers(0, 60)).map(lambda t: (t[0], t[0] + t[1]))


@FAST
@given(seq=STD_SEQS)
def test_standard_sequence_has_no_counts(seq: str) -> None:
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "p.faa"
        p.write_text(f">a\n{seq.lower()}\n")
        (r,) = read_protein_fasta(p)
    assert r.sequence == seq
    assert sum(r.nonstandard.values()) == 0
    assert dict(r.other) == {}
    assert r.internal_stops == 0


@FAST
@given(token=TOKENS)
def test_resolve_protein_id_is_idempotent(token: str) -> None:
    group_id, members = resolve_protein_id(token)
    assert group_id == ";".join(members)
    # members are de-duplicated after resolution, first one kept (M1-03 review decision)
    assert len(members) == len(set(members))
    assert len(members) <= token.count(";") + 1
    assert all("|" not in m and m for m in members)
    assert resolve_protein_id(group_id) == (group_id, members)


@FAST
@given(length=st.integers(0, 200), bounds=BOUNDS)
def test_route_by_length_partitions_lengths(length: int, bounds: tuple[int, int]) -> None:
    lo, hi = bounds
    route = route_by_length(length, peptide_len=bounds)
    expected = "dropped" if length < lo else "peptide" if length <= hi else "protein"
    assert route == expected
