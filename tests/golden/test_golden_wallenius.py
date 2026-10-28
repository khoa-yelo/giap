"""M8-13 §1.1: GIPPA's Wallenius ORA against goseq and the exact reference, on the committed
golden data (golden layer: every PR, CI, no R). ADR-0044 (accepted 2026-10-23) amends the gate.

- Universe: all 3,000 genes of ``inputs/ora_genes.tsv.gz`` (goseq ran with
  ``use_genes_without_cat = TRUE``), ``p_detect`` as the PWF verbatim.
- Counts equal goseq's ``numInCat`` / ``numDEInCat`` and ω equals goseq's weight identity
  (M8-01 example 3) to relative 1e-12, on every one of the 150 categories (none skipped).
- Criterion 1 (ADR-0044 decision 1): ``|p - p_exact| ≤ wallenius_rel_tol · p_exact`` on all 300
  tails against ``tests/data/golden/wallenius/exact_tails.tsv.gz``; that file is independently
  checked against a test-side urn recursion from the definition (``helpers.golden_r_gate``).
- goseq as second reference (ADR-0044 decision 2):
  ``|p - p_goseq| ≤ max(wallenius_rel_tol · p_goseq, wallenius_goseq_abs_tol)``.
- Example 3: ``detection_enrichment``'s ``p_classic`` is ``wallenius_ora`` exactly.
- Example 10 (mutation, ω inverted): the exact rule catches all 150 categories, the goseq rule
  all but the 5 deep-tail ones.
- A5 vacuity floor (ADR-0044 decision 4) and the §5 ``wallenius`` report section.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.golden_r_gate import (
    GOLDEN,
    TAILS,
    decade,
    golden_r_tol,
    num,
    read_exact_tails,
    read_gz,
    update_report,
    urn_tails,
)

N_CATEGORIES = 150
DEEP_TAIL = {"CAT_004", "CAT_010", "CAT_011", "CAT_012", "CAT_013"}  # ADR-0044 decision 2


@dataclass(frozen=True)
class OraWorld:
    genes: list[str]
    detected: np.ndarray
    p_detect: np.ndarray
    cats: dict[str, np.ndarray]  # category -> in-category mask, file order
    goseq: dict[str, dict[str, str]]
    gippa: dict[str, Any]  # category -> OraResult


@pytest.fixture(scope="module")
def world() -> OraWorld:
    from gippa.stats.wallenius import wallenius_ora

    rows = read_gz(GOLDEN / "inputs" / "ora_genes.tsv.gz")
    genes = [r["gene"] for r in rows]
    det = np.array([int(r["detected"]) for r in rows], dtype=bool)
    p = np.array([float(r["p_detect"]) for r in rows], dtype=np.float64)
    index = {g: i for i, g in enumerate(genes)}
    cats: dict[str, np.ndarray] = {}
    for r in read_gz(GOLDEN / "inputs" / "ora_categories.tsv.gz"):
        mask = cats.setdefault(r["category"], np.zeros(len(genes), dtype=bool))
        mask[index[r["gene"]]] = True
    goseq = {r["category"]: r for r in read_gz(GOLDEN / "outputs" / "goseq.tsv.gz")}
    gippa = {c: wallenius_ora(det, p, m) for c, m in cats.items()}
    return OraWorld(genes, det, p, cats, goseq, gippa)


@pytest.fixture(scope="module")
def urn(world: OraWorld) -> dict[tuple[str, str], float]:
    out = {}
    for c, o in world.gippa.items():
        over, under = urn_tails(o.k, o.M, o.n, o.N, o.odds)
        out[(c, "over")], out[(c, "under")] = over, under
    return out


def _p(o: Any, tail: str) -> float:
    return float(o.p_over if tail == "over" else o.p_under)


def _goseq_p(row: dict[str, str], tail: str) -> float:
    return num(row[f"{tail}_represented_pvalue"])


def _rel(p: float, ref: float) -> float:
    return abs(p - ref) / ref if ref > 0 else (0.0 if p == ref else math.inf)


# ---------------------------------------------------------------- universe, counts, odds


def test_universe_is_every_golden_gene(world: OraWorld) -> None:
    assert len(world.cats) == N_CATEGORIES == len(world.goseq)
    assert list(world.cats) == list(world.goseq)
    for o in world.gippa.values():
        assert o.M == len(world.genes) == 3000
        assert int(world.detected.sum()) == o.N


def test_example_1_counts_equal_goseq_on_all_150(world: OraWorld) -> None:
    bad = [
        c
        for c, o in world.gippa.items()
        if (o.n, o.k) != (int(world.goseq[c]["numInCat"]), int(world.goseq[c]["numDEInCat"]))
    ]
    assert bad == []


def test_odds_equal_the_goseq_weight_identity(world: OraWorld) -> None:
    """avg_weight·(N_genes - n)/(Σpwf - n·avg_weight), avg_weight = mean pwf in the category."""
    total = math.fsum(world.p_detect.tolist())
    for c, mask in world.cats.items():
        n = int(mask.sum())
        avg = math.fsum(world.p_detect[mask].tolist()) / n
        identity = avg * (mask.size - n) / (total - n * avg)
        assert world.gippa[c].odds == pytest.approx(identity, rel=1e-12, abs=0.0), c


# ---------------------------------------------------------------- criterion 1: the exact reference


def test_example_2_tails_match_the_exact_reference(world: OraWorld) -> None:
    """ADR-0044 decision 1: all 300 tails within wallenius_rel_tol of the committed exact p."""
    tol = golden_r_tol()["wallenius_rel_tol"]
    exact = read_exact_tails()
    assert set(exact) == {(c, t) for c in world.cats for t in TAILS}
    bad = {
        (c, t): (_p(o, t), exact[(c, t)])
        for c, o in world.gippa.items()
        for t in TAILS
        if not abs(_p(o, t) - exact[(c, t)]) <= tol * exact[(c, t)]
    }
    assert bad == {}


def test_exact_reference_agrees_with_the_urn_definition(
    urn: dict[tuple[str, str], float],
) -> None:
    """The committed reference is checked against an independent test-side computation."""
    exact = read_exact_tails()
    bad = {key: (v, urn[key]) for key, v in exact.items() if not _rel(v, urn[key]) <= 1e-12}
    assert bad == {}


def test_tails_match_the_urn_definition(world: OraWorld, urn: dict[tuple[str, str], float]) -> None:
    """GIPPA against the definition at wallenius_rel_tol (needs no committed reference)."""
    tol = golden_r_tol()["wallenius_rel_tol"]
    bad = {
        (c, t): (_p(o, t), urn[(c, t)])
        for c, o in world.gippa.items()
        for t in TAILS
        if not abs(_p(o, t) - urn[(c, t)]) <= tol * urn[(c, t)]
    }
    assert bad == {}


# ---------------------------------------------------------------- goseq, second reference


def _goseq_rule(p: float, p_goseq: float, tol: dict[str, float]) -> bool:
    bound = max(tol["wallenius_rel_tol"] * p_goseq, tol["wallenius_goseq_abs_tol"])
    return abs(p - p_goseq) <= bound


def _goseq_tol() -> dict[str, float]:
    tol = golden_r_tol()
    assert "wallenius_goseq_abs_tol" in tol, (
        "tolerances.yaml golden_r lacks wallenius_goseq_abs_tol (1.0e-5, ADR-0044 decision 2; "
        "added with Gate-Change-Approved: ADR-0044)"
    )
    return tol


def test_tails_within_goseq_documented_accuracy(world: OraWorld) -> None:
    tol = _goseq_tol()
    bad = {
        (c, t): (_p(o, t), _goseq_p(world.goseq[c], t))
        for c, o in world.gippa.items()
        for t in TAILS
        if not _goseq_rule(_p(o, t), _goseq_p(world.goseq[c], t), tol)
    }
    assert bad == {}


def test_example_10_inverted_odds_is_caught(
    world: OraWorld, urn: dict[tuple[str, str], float]
) -> None:
    """Mutation: ``category_odds`` returning out/in. The definition rule catches every category;
    the goseq rule every category but the deep-tail five (ADR-0044: 145 of 150)."""
    from gippa.stats.wallenius import wallenius_tails

    tol = _goseq_tol()
    by_urn, by_goseq = set(), set()
    for c, o in world.gippa.items():
        mutant = wallenius_tails(o.k, o.M, o.n, o.N, 1.0 / o.odds)
        for t, p in zip(TAILS, mutant, strict=True):
            if not abs(p - urn[(c, t)]) <= tol["wallenius_rel_tol"] * urn[(c, t)]:
                by_urn.add(c)
            if not _goseq_rule(p, _goseq_p(world.goseq[c], t), tol):
                by_goseq.add(c)
    assert by_urn == set(world.cats)
    assert set(world.cats) - by_goseq <= DEEP_TAIL


# ---------------------------------------------------------------- example 3: p_classic wiring


def test_example_3_p_classic_is_wallenius_ora(world: OraWorld) -> None:
    """One gene per TU, all eligible, no folds: ``p_classic`` == ``wallenius_ora`` (float ==)."""
    from gippa.function.collections import GeneSets
    from gippa.function.ora import detection_enrichment

    genes = pd.DataFrame(
        {
            "locus_tag": world.genes,
            "tu_id": ["TU_" + g for g in world.genes],
            "p_detect": world.p_detect,
            "detected": world.detected,
            "partial": False,
        }
    )
    sets = {c: frozenset(np.asarray(world.genes)[m].tolist()) for c, m in world.cats.items()}
    coll = GeneSets("kegg", sets, {}, {}, 1.0, 0)
    df = detection_enrichment(genes, [coll], (1, 3000), 0.05, folds=None)
    got = {(r.term, r.direction): r.p_classic for r in df.itertuples()}
    assert set(got) == {(c, t) for c in world.cats for t in TAILS}
    bad = {key: (v, _p(world.gippa[key[0]], key[1])) for key, v in got.items()
           if v != _p(world.gippa[key[0]], key[1])}  # fmt: skip
    assert bad == {}


# ---------------------------------------------------------------- A5 vacuity floor and §5 report


def test_vacuity_floor_on_the_exact_tails() -> None:
    exact = read_exact_tails().values()
    assert sum(0.0 < p < 1e-20 for p in exact) >= 2
    assert sum(1e-20 <= p < 1e-8 for p in exact) >= 1


def test_report_wallenius_section(world: OraWorld, urn: dict[tuple[str, str], float]) -> None:
    exact = read_exact_tails()
    tol = golden_r_tol()["wallenius_rel_tol"]
    refs = {
        "goseq": {(c, t): _goseq_p(world.goseq[c], t) for c in world.cats for t in TAILS},
        "exact": exact,
        "urn": urn,
    }
    payload: dict[str, Any] = {"categories": len(world.cats)}
    for name, ref in refs.items():
        part: dict[str, Any] = {}
        for t in TAILS:
            rel = {
                c: _rel(_p(world.gippa[c], t), ref[(c, t)]) for c in world.cats if ref[(c, t)] > 0
            }  # p_ref = 0 rows are listed separately
            worst = max(rel, key=rel.__getitem__)
            part[t] = {"max_rel": rel[worst], "argmax": worst,
                       "n_over_tol": sum(v > tol for v in rel.values())}  # fmt: skip
        part["p_ref_zero"] = sorted(f"{c}:{t}" for (c, t), v in ref.items() if v == 0.0)
        part["decades"] = dict(Counter(decade(v) for v in ref.values()))
        payload[name] = part
    path = update_report("wallenius", payload)
    assert path.is_file()
    assert payload["exact"]["over"]["n_over_tol"] == payload["exact"]["under"]["n_over_tol"] == 0
