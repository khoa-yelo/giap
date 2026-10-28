"""The generic GO-elim walker (issue M8-01; GIPPA_SPEC.md §7.7.3 "GO redundancy (elim, Alexa
2006)", §7.7.5 "Unit: elim on a toy DAG").

Interface pinned here (docs/specs/M8-01-stats-wallenius-elim.md, module ``gippa.function.elim``):
``elim(term_genes, parents, test, threshold=0.01) -> {term: (p_classic, p_elim)}``.

- terms are visited deepest first (level = longest path to a root; ties by term ID);
- a term is tested on its genes minus those marked removed for it; a term whose elim p is
  strictly below ``threshold`` (M8-index: strict, as topGO) marks its current genes removed for
  every ancestor (transitive closure of ``parents``);
- ``p_classic`` is the test on the full gene set; an empty remaining set gives elim p NaN and the
  test is never called on it; cycles raise ``ValueError`` naming the cycle; the test is never
  called twice with the same set for one term.

Expected values below are worked by hand. The toy test is a lookup table keyed by the exact gene
set, so a call with any unexpected set fails loudly (``KeyError``).
"""

from __future__ import annotations

import copy
import inspect
import math
from collections import Counter
from collections.abc import Mapping

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.function.elim import elim


def fs(*x: int) -> frozenset[int]:
    return frozenset(x)


class Table:
    """Toy test: p from a table keyed by gene set; records every call."""

    def __init__(self, table: Mapping[frozenset[int], float]) -> None:
        self.table = dict(table)
        self.calls: list[frozenset[int]] = []

    def __call__(self, genes: frozenset[int]) -> float:
        assert isinstance(genes, frozenset)
        self.calls.append(genes)
        return self.table[genes]


# Example 8 DAG: R <- A <- {A1, A2}; R <- B.
GENES8: dict[str, frozenset[int]] = {
    "A1": fs(1, 2, 3),
    "A2": fs(4),
    "A": fs(1, 2, 3, 4, 5),
    "B": fs(6),
    "R": fs(1, 2, 3, 4, 5, 6),
}
PARENTS8: dict[str, frozenset[str]] = {
    "A1": frozenset({"A"}),
    "A2": frozenset({"A"}),
    "A": frozenset({"R"}),
    "B": frozenset({"R"}),
    "R": frozenset(),
}


def table8(p_a1: float = 0.001, p_a_elim: float = 0.04) -> Table:
    return Table(
        {
            fs(1, 2, 3): p_a1,
            fs(4): 0.3,
            fs(1, 2, 3, 4, 5): 0.2,
            fs(4, 5): p_a_elim,
            fs(6): 0.6,
            fs(1, 2, 3, 4, 5, 6): 0.1,
            fs(4, 5, 6): 0.07,
        }
    )


def test_signature_and_default_threshold() -> None:
    sig = inspect.signature(elim)
    assert list(sig.parameters) == ["term_genes", "parents", "test", "threshold"]
    assert sig.parameters["threshold"].default == 0.01


def test_example8_removal_propagates_to_all_ancestors() -> None:
    t = table8()
    out = elim(GENES8, PARENTS8, t)
    assert set(out) == set(GENES8)
    assert out["A1"] == (0.001, 0.001)
    assert out["A2"] == (0.3, 0.3)
    assert out["A"] == (0.2, 0.04)  # classic on {1..5}, elim on {4, 5}
    assert out["B"] == (0.6, 0.6)
    assert out["R"] == (0.1, 0.07)  # elim on {4, 5, 6}: A's p = 0.04 is not < 0.01
    assert all(isinstance(v, float) for pair in out.values() for v in pair)


def test_example8_each_set_tested_once() -> None:
    t = table8()
    elim(GENES8, PARENTS8, t)
    counts = Counter(t.calls)
    assert set(counts) == set(t.table)
    assert all(c == 1 for c in counts.values()), counts


def test_example8_chain_of_significant_terms() -> None:
    # A significant on {4, 5} as well: R loses {1..5} and is tested on {6}.
    t = table8(p_a_elim=0.004)
    out = elim(GENES8, PARENTS8, t)
    assert out["A"] == (0.2, 0.004)
    assert out["R"] == (0.1, 0.6)
    assert fs(4, 5, 6) not in t.calls


def test_example9_threshold_boundary_is_not_significant() -> None:
    t = table8(p_a1=0.01)
    out = elim(GENES8, PARENTS8, t)
    assert out["A1"] == (0.01, 0.01)
    assert out["A"] == (0.2, 0.2)  # tested on its full set
    assert out["R"] == (0.1, 0.1)
    assert fs(4, 5) not in t.calls


def test_custom_threshold() -> None:
    t = table8(p_a1=0.01)
    out = elim(GENES8, PARENTS8, t, threshold=0.05)
    assert out["A1"] == (0.01, 0.01)
    assert out["A"] == (0.2, 0.04)  # 0.04 < 0.05: A is now significant too
    assert out["R"] == (0.1, 0.6)


def test_significance_uses_elim_p_not_classic_p() -> None:
    # R <- A <- A1. A's classic p is tiny but its elim p (on {2}) is not, so A removes nothing.
    genes = {"A1": fs(1), "A": fs(1, 2), "R": fs(1, 2, 3)}
    parents = {"A1": frozenset({"A"}), "A": frozenset({"R"}), "R": frozenset()}
    t = Table({fs(1): 0.001, fs(1, 2): 0.001, fs(2): 0.5, fs(2, 3): 0.25, fs(1, 2, 3): 0.1})
    out = elim(genes, parents, t)
    assert out["A"] == (0.001, 0.5)
    assert out["R"] == (0.1, 0.25)


def test_level_is_longest_path_to_root() -> None:
    # D's parents are A and R, so D's level is 2 (D -> A -> R) and it is visited before A.
    # With shortest-path levels D and A would tie and A ("A" < "D") would be tested on {4, 5}.
    genes = {"D": fs(5), "A": fs(4, 5), "R": fs(4, 5, 6)}
    parents = {"D": frozenset({"A", "R"}), "A": frozenset({"R"}), "R": frozenset()}
    t = Table({fs(5): 0.001, fs(4, 5): 0.2, fs(4): 0.3, fs(4, 6): 0.4, fs(4, 5, 6): 0.1})
    out = elim(genes, parents, t)
    assert out["D"] == (0.001, 0.001)
    assert out["A"] == (0.2, 0.3)
    assert out["R"] == (0.1, 0.4)


def test_empty_after_removal_is_nan_and_not_tested() -> None:
    # R <- A <- A1 with A's genes equal to A1's: A has nothing left after A1 is significant.
    genes = {"A1": fs(1, 2), "A": fs(1, 2), "R": fs(1, 2, 3)}
    parents = {"A1": frozenset({"A"}), "A": frozenset({"R"}), "R": frozenset()}
    t = Table({fs(1, 2): 0.001, fs(3): 0.3, fs(1, 2, 3): 0.2})
    out = elim(genes, parents, t)
    assert out["A1"] == (0.001, 0.001)
    p_classic, p_elim = out["A"]
    assert p_classic == 0.001
    assert math.isnan(p_elim)
    assert out["R"] == (0.2, 0.3)  # only A1's genes are removed; the NaN term adds nothing
    assert frozenset() not in t.calls


def test_two_roots_sharing_a_child() -> None:
    genes = {"C": fs(1), "R1": fs(1, 2), "R2": fs(1, 3)}
    parents = {"C": frozenset({"R1", "R2"}), "R1": frozenset(), "R2": frozenset()}
    t = Table({fs(1): 0.001, fs(1, 2): 0.2, fs(2): 0.5, fs(1, 3): 0.3, fs(3): 0.7})
    out = elim(genes, parents, t)
    assert out["C"] == (0.001, 0.001)
    assert out["R1"] == (0.2, 0.5)
    assert out["R2"] == (0.3, 0.7)


def test_removal_does_not_reach_siblings_or_descendants() -> None:
    # A1 significant: its sibling A2 and the unrelated B keep their full sets.
    genes = {"A1": fs(1, 2), "A2": fs(2, 3), "A": fs(1, 2, 3), "B": fs(1)}
    parents = {
        "A1": frozenset({"A"}),
        "A2": frozenset({"A"}),
        "A": frozenset(),
        "B": frozenset(),
    }
    t = Table({fs(1, 2): 0.001, fs(2, 3): 0.5, fs(1, 2, 3): 0.2, fs(3): 0.6, fs(1): 0.4})
    out = elim(genes, parents, t)
    assert out["A2"] == (0.5, 0.5)
    assert out["B"] == (0.4, 0.4)
    assert out["A"] == (0.2, 0.6)


@pytest.mark.parametrize(
    "parents",
    [
        {"GO:0000001": frozenset({"GO:0000002"}), "GO:0000002": frozenset({"GO:0000001"})},
        {
            "GO:0000001": frozenset({"GO:0000002"}),
            "GO:0000002": frozenset({"GO:0000003"}),
            "GO:0000003": frozenset({"GO:0000001"}),
        },
    ],
)
def test_example10_cycle_raises_naming_it(parents: dict[str, frozenset[str]]) -> None:
    genes = {t: fs(1) for t in parents}
    with pytest.raises(ValueError) as info:
        elim(genes, parents, Table({fs(1): 0.5}))
    for term in parents:
        assert term in str(info.value)


def test_self_loop_raises() -> None:
    with pytest.raises(ValueError, match="GO:0000009"):
        elim(
            {"GO:0000009": fs(1)},
            {"GO:0000009": frozenset({"GO:0000009"})},
            Table({fs(1): 0.5}),
        )


def test_cycle_above_an_acyclic_part_raises() -> None:
    genes = {"L": fs(1), "X": fs(1, 2), "Y": fs(1, 3)}
    parents = {"L": frozenset({"X"}), "X": frozenset({"Y"}), "Y": frozenset({"X"})}
    with pytest.raises(ValueError, match="X"):
        elim(genes, parents, Table({}))


def test_inputs_not_mutated_and_deterministic() -> None:
    genes = copy.deepcopy(GENES8)
    parents = copy.deepcopy(PARENTS8)
    out1 = elim(genes, parents, table8())
    out2 = elim(genes, parents, table8())
    assert out1 == out2
    assert genes == GENES8 and parents == PARENTS8
    # dict insertion order must not change the result
    rev_genes = dict(reversed(list(GENES8.items())))
    rev_parents = dict(reversed(list(PARENTS8.items())))
    assert elim(rev_genes, rev_parents, table8()) == out1


# ----------------------------------------------------------------------------- properties


def _hash_p(genes: frozenset[int]) -> float:
    """Deterministic p in [0, 1) from the gene set."""
    h = 0
    for g in sorted(genes):
        h = (h * 1_000_003 + g + 1) % 2_147_483_647
    return (h % 10_007) / 10_007


@st.composite
def _dag(draw: st.DrawFn) -> tuple[dict[str, frozenset[int]], dict[str, frozenset[str]]]:
    m = draw(st.integers(1, 8))
    names = [f"T{i:02d}" for i in range(m)]
    genes: dict[str, frozenset[int]] = {}
    parents: dict[str, frozenset[str]] = {}
    for i, name in enumerate(names):
        g = draw(st.frozensets(st.integers(0, 11), min_size=1, max_size=8))
        ps = draw(st.frozensets(st.sampled_from(names[:i]), max_size=3)) if i else frozenset()
        genes[name] = g
        parents[name] = ps
    # true-path rule: a parent annotates every gene of its children
    for name in reversed(names):
        for p in parents[name]:
            genes[p] = genes[p] | genes[name]
    return genes, parents


def _leaves(parents: Mapping[str, frozenset[str]]) -> set[str]:
    has_child = {p for ps in parents.values() for p in ps}
    return set(parents) - has_child


@settings(max_examples=60, deadline=None)
@given(_dag())
def test_property_threshold_zero_gives_classic(
    dag: tuple[dict[str, frozenset[int]], dict[str, frozenset[str]]],
) -> None:
    genes, parents = dag
    out = elim(genes, parents, _hash_p, threshold=0.0)
    assert set(out) == set(genes)
    for term, (p_classic, p_elim) in out.items():
        assert p_classic == _hash_p(genes[term])
        assert p_elim == p_classic


@settings(max_examples=60, deadline=None)
@given(_dag(), st.floats(0.0, 1.0))
def test_property_elim_sets_are_subsets_and_leaves_unchanged(
    dag: tuple[dict[str, frozenset[int]], dict[str, frozenset[str]]], threshold: float
) -> None:
    genes, parents = dag
    calls: list[frozenset[int]] = []

    def test(g: frozenset[int]) -> float:
        calls.append(g)
        return _hash_p(g)

    out = elim(genes, parents, test, threshold=threshold)
    assert frozenset() not in calls
    assert all(any(c <= g for g in genes.values()) for c in calls)
    for term in _leaves(parents):
        assert out[term][1] == out[term][0]
    for p_classic, p_elim in out.values():
        assert 0.0 <= p_classic <= 1.0
        assert math.isnan(p_elim) or 0.0 <= p_elim <= 1.0
