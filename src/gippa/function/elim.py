"""The GO-elim walker (Alexa et al. 2006; GIPPA_SPEC.md §7.7.3 "GO redundancy", §7.7.5).

Generic over any DAG and any test. Terms are visited deepest first (level = longest path to a
root; ties by term ID). Each term is tested on its genes minus those removed for it; a term whose
elim p is strictly below ``threshold`` (as topGO's ``< cutOff``) removes its tested genes from
every ancestor. Cycles are detected before any test runs.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping

TestFn = Callable[[frozenset[int]], float]


def _levels(parents: Mapping[str, frozenset[str]], terms: list[str]) -> dict[str, int]:
    """Longest path to a root per term; ``ValueError`` naming the terms of any cycle."""
    level: dict[str, int] = {}
    on_path: list[str] = []
    for start in terms:
        if start in level:
            continue
        stack: list[tuple[str, list[str]]] = [(start, sorted(parents.get(start, ())))]
        on_path.append(start)
        while stack:
            term, todo = stack[-1]
            if todo:
                nxt = todo.pop()
                if nxt in on_path:
                    cycle = [*on_path[on_path.index(nxt) :], nxt]
                    raise ValueError(f"cycle in the term DAG: {' -> '.join(cycle)}")
                if nxt not in level:
                    stack.append((nxt, sorted(parents.get(nxt, ()))))
                    on_path.append(nxt)
                continue
            level[term] = 1 + max((level[q] for q in parents.get(term, ())), default=-1)
            stack.pop()
            on_path.pop()
    return level


def _ancestors(parents: Mapping[str, frozenset[str]], order: list[str]) -> dict[str, set[str]]:
    """Transitive closure of ``parents``; ``order`` lists every term after its parents."""
    anc: dict[str, set[str]] = {}
    for term in order:
        acc: set[str] = set()
        for q in parents.get(term, ()):
            acc.add(q)
            acc |= anc[q]
        anc[term] = acc
    return anc


def elim(
    term_genes: Mapping[str, frozenset[int]],
    parents: Mapping[str, frozenset[str]],
    test: TestFn,
    threshold: float = 0.01,
) -> dict[str, tuple[float, float]]:
    """{term: (p_classic, p_elim)}; elim p is NaN when nothing is left after removal."""
    unknown = sorted({q for ps in parents.values() for q in ps} - set(term_genes))
    unknown += sorted(set(parents) - set(term_genes))
    if unknown:
        raise ValueError(f"parents name terms without a gene set: {unknown}")
    terms = sorted(term_genes)
    level = _levels(parents, terms)
    anc = _ancestors(parents, sorted(terms, key=lambda t: (level[t], t)))
    removed: dict[str, set[int]] = {t: set() for t in terms}
    out: dict[str, tuple[float, float]] = {}
    for term in sorted(terms, key=lambda t: (-level[t], t)):
        genes = frozenset(term_genes[term])
        p_classic = float(test(genes)) if genes else float("nan")
        left = genes - removed[term]
        if not left:
            out[term] = (p_classic, float("nan"))
            continue
        p_elim = p_classic if left == genes else float(test(left))
        out[term] = (p_classic, p_elim)
        if p_elim < threshold:
            for a in anc[term]:
                removed[a] |= left
    return out
