"""M11-16c: the M11-14w fixture generator gains code 25
(docs/specs/M11-16-translation-tables-4-25.md item 8: "The M11-14w fixture generator
(``n_genes = 300``, 400 codons, ...) gains code 25 (TGA -> Gly, rate ``r``). The code-4 path
keeps its seeds and values"; the code-4/11 values are pinned in ``test_m11_16c_identity.py``).

Interface assumed: ``helpers.genetic_code_fixture.synthetic_genome(seed, gc, 25, tga_rate=r)``
builds the code-25 fixture with the code-4 layout: ``spacer`` bp, then a gene of ``n_codons``
codons (ATG, sense codons with TGA at rate ``r``, TAA), repeated ``n_genes`` times, plus a final
spacer; half of the genes reverse-complemented. Spec example 9 ("the synthetic code-25 fixture
at 11 -> warn") needs its long-ORF ratio below 0.55, checked here with the brute-force reference.
"""

from __future__ import annotations

import math

import pytest
from helpers.genetic_code_fixture import S4, S11, reference_cover, revcomp, synthetic_genome

N_GENES, N_CODONS, SPACER = 300, 400, 60


def _genes(seq: str) -> list[str]:
    """The genes of a fixture in transcription orientation."""
    assert len(seq) == N_GENES * (SPACER + 3 * N_CODONS) + SPACER
    out = []
    for i in range(N_GENES):
        a = SPACER + i * (SPACER + 3 * N_CODONS)
        g = seq[a : a + 3 * N_CODONS]
        out.append(g if g.startswith("ATG") and g.endswith("TAA") else revcomp(g))
    return out


def _interior(genes: list[str]) -> list[str]:
    return [g[i : i + 3] for g in genes for i in range(3, 3 * N_CODONS - 3, 3)]


@pytest.mark.parametrize(("seed", "gc", "rate"), [(7, 0.28, 0.01), (3, 0.5, 0.01), (5, 0.4, 0.05)])
def test_code25_fixture_has_tga_at_rate(seed: int, gc: float, rate: float) -> None:
    seq = synthetic_genome(seed, gc, 25, tga_rate=rate)
    genes = _genes(seq)
    minus = sum(not seq[SPACER + i * (SPACER + 3 * N_CODONS) :].startswith("ATG")
                for i in range(N_GENES))  # fmt: skip
    assert minus == N_GENES // 2
    for g in genes:
        assert g.startswith("ATG") and g.endswith("TAA")
    cod = _interior(genes)
    assert not {"TAA", "TAG"} & set(cod)  # no stop of table 25 inside a gene
    k, n = cod.count("TGA"), len(cod)
    se = math.sqrt(rate * (1 - rate) / n)
    assert abs(k / n - rate) <= 3 * se, (k, n, rate)


def test_code25_fixture_long_orf_ratio_below_cut() -> None:
    """At table 11 the code-25 fixture looks like a TGA-sense genome (ratio < 0.55)."""
    seq = synthetic_genome(7, 0.28, 25)
    bp11, _ = reference_cover(seq, S11, 300)
    bp4, n4 = reference_cover(seq, S4, 300)
    assert n4 >= 20
    assert 20 * bp11 < 11 * bp4, (bp11, bp4)


def test_code25_fixture_is_seeded_and_not_code11() -> None:
    assert synthetic_genome(11, 0.4, 25) == synthetic_genome(11, 0.4, 25)
    assert synthetic_genome(11, 0.4, 25) != synthetic_genome(12, 0.4, 25)
    assert synthetic_genome(11, 0.4, 25) != synthetic_genome(11, 0.4, 11)
