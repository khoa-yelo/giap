"""M8-05 code review (PR #171): ``parent_group`` with multiple parents and non-significant
ancestors, ``redundant`` decided by the BH-substituted q (not the re-test p alone), and the
nullable boolean dtype of ``redundant``."""

# ruff: noqa: F811
from __future__ import annotations

from collections.abc import Callable

from helpers.gsea_data import gene_sets
from test_gsea import FakeEngine, _family, _go_frame, _row, _run, fake  # noqa: F401


def test_parent_group_prefers_a_significant_ancestor_over_a_shallower_self(
    fake: Callable[..., FakeEngine],
) -> None:
    """Z <- Y <- A <- T and T <- X (top level): depth(T) = 1 < depth(A) = 2. With A and T
    significant, T groups under A, not under itself (B1)."""
    genes, tags = _go_frame()
    sets = {"GO:T": tags[:11], "GO:A": tags[:20], "GO:Y": tags[:30], "GO:Z": tags[:40],
            "GO:X": tags[:11] + tags[-39:]}  # fmt: skip
    parents = {"GO:T": {"GO:A", "GO:X"}, "GO:A": {"GO:Y"}, "GO:Y": {"GO:Z"}, "GO:Z": set(),
               "GO:X": set()}  # fmt: skip
    fake({11: 0.001, 20: 0.001, 30: 0.9, 40: 0.9, 50: 0.9, 9: 0.001})
    out = _run(genes, [gene_sets("go_bp", sets, parents=parents)])
    assert _row(out, "GO:T")["parent_group"] == "GO:A"
    assert _row(out, "GO:A")["parent_group"] == "GO:A"


def test_parent_group_skips_non_significant_ancestors(fake: Callable[..., FakeEngine]) -> None:
    """G <- M <- P <- Q with G and P significant, M and Q not: G and P group under P."""
    genes, tags = _go_frame()
    sets = {"GO:G": tags[:11], "GO:M": tags[:24], "GO:P": tags[:30], "GO:Q": tags[:45]}
    parents = {"GO:G": {"GO:M"}, "GO:M": {"GO:P"}, "GO:P": {"GO:Q"}, "GO:Q": set()}
    fake({11: 0.001, 24: 0.9, 30: 0.001, 45: 0.9, 19: 0.001})
    out = _run(genes, [gene_sets("go_bp", sets, parents=parents)])
    assert _row(out, "GO:G")["parent_group"] == "GO:P"
    assert _row(out, "GO:P")["parent_group"] == "GO:P"


def test_redundant_uses_the_substituted_bh_q_not_the_retest_p(
    fake: Callable[..., FakeEngine],
) -> None:
    """Re-test p 0.04 < alpha, but BH of [0.001, 0.04, 0.5, 0.6, 0.7] gives q 0.1 > alpha:
    redundant (decisions 3 and 6)."""
    genes, tags = _go_frame()
    sets, parents = _family(tags, 15, 3, [40, 41, 42])
    fake({15: 0.001, 18: 0.002, 3: 0.04, 40: 0.5, 41: 0.6, 42: 0.7})
    out = _run(genes, [gene_sets("go_bp", sets, parents=parents)])
    assert _row(out, "GO:P")["redundant"] == True  # noqa: E712
    assert str(out["redundant"].dtype) == "boolean"
