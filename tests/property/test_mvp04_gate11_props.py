"""Invariants of the fixed MVP-04 gate 11 scorer (contract: tests/unit/test_mvp04_gate11_scorer.py).

* I/L folding: relabelling the ALT residue of true I/L sites among ``I``, ``L`` and ``J`` (on
  either side) never changes precision or recall.
* initiator exclusion: adding initiator rows (``start_codon``-flagged, or ALT M over a
  GTG/TTG/CTG novel-ORF start) never changes precision or recall and raises
  ``n_initiator_m_excluded`` by exactly their number.
* recall denominator: raising ``mapping.tolerant.min_len`` never adds a recall site.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

from helpers.mvp04_gate11_data import FLAGS, N_SITES, base_case, gene
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

SETTINGS = settings(max_examples=25, deadline=None,
                    suppress_health_check=[HealthCheck.function_scoped_fixture])  # fmt: skip
IL = st.sampled_from("ILJ")


@SETTINGS
@given(st.lists(st.tuples(IL, IL), min_size=N_SITES, max_size=N_SITES),
       st.integers(0, N_SITES))  # fmt: skip
def test_il_relabelling_never_changes_the_scores(labels: list[tuple[str, str]], n_il: int) -> None:
    with tempfile.TemporaryDirectory() as d:
        c = base_case(Path(d))
        ref = c.gate()
        for k, (t, a) in enumerate(labels[:n_il]):
            c.alt[k][5], c.variants[k][8] = t, a
        res = c.gate()
    assert res["precision"] == ref["precision"] == 1.0
    assert res["recall"] == ref["recall"] == 1.0


@SETTINGS
@given(st.integers(0, 6), st.integers(0, 6), st.sampled_from(["GTG", "TTG", "CTG"]),
       st.integers(0, 5))  # fmt: skip
def test_initiator_rows_only_move_the_initiator_count(
    n_flag: int, n_novel: int, start: str, n_false: int
) -> None:
    with tempfile.TemporaryDirectory() as d:
        c = base_case(Path(d))
        for k in range(n_false):
            c.var(300_000 + 3 * k, "W")
        ref = c.gate()
        for k in range(n_flag):
            c.var(400_000 + 300 * k, "M", gene=gene(500 + k), ref_codon="GTG", ref="V",
                  flags=FLAGS + ";start_codon")  # fmt: skip
        for k in range(n_novel):
            plus = k % 2 == 0
            c.var(100_001 if plus else 200_298, "M", strand="+" if plus else "-",
                  ref_codon=start, ref="L" if start != "GTG" else "V")  # fmt: skip
        res = c.gate()
    assert res["precision"] == ref["precision"]
    assert res["recall"] == ref["recall"]
    assert res["n_initiator_m_excluded"] == ref["n_initiator_m_excluded"] + n_flag + n_novel
    assert res["value"] == ref["value"] + n_flag + n_novel


@SETTINGS
@given(st.lists(st.integers(7, 14), min_size=1, max_size=8), st.integers(6, 14),
       st.integers(0, 4))  # fmt: skip
def test_raising_min_len_never_adds_a_recall_site(lens: list[int], lo: int, step: int) -> None:
    with tempfile.TemporaryDirectory() as d:
        c = base_case(Path(d))
        for j, n in enumerate(lens):
            k = N_SITES + j
            c.alt.append([c.alt[0][0], 1000 + 600 * k, "+", gene(k), "A", "V", 2, "false"])
            c.true_peps(gene(k), 1000 + 600 * k, length=n)
        c.min_len = lo
        a = c.gate()["n_recall_sites"]
        c.min_len = lo + step
        b = c.gate()["n_recall_sites"]
    assert b <= a
    assert a == N_SITES * (lo <= 12) + sum(n >= lo for n in lens)
