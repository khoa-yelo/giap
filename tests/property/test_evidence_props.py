"""Hypothesis properties of M4-11 accepted evidence (docs/specs/M4-11-accepted-evidence.md).

Acceptance criterion: relabelling peptides leaves every count unchanged; a duplicate C11 row for the
same peptide and locus changes nothing; ``n_accepted >= n_unique``. Also: agreement with the
independent oracle (tests/helpers/evidence_oracle.py) on random small scenarios, including genes
crossing the origin of a circular replicon, and monotonicity in ``input.q_value_max``.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import evidence_oracle as eo

SETTINGS = settings(max_examples=80, deadline=None, suppress_health_check=[HealthCheck.too_slow])
L = 3000
COUNT_COLS = (
    "n_peptides_unique",
    "n_peptides_multi_identical",
    "n_peptides_shared",
    "n_peptides_tolerant",
    "n_accepted",
)
CLASSES = (
    "cds_inframe",
    "sorf",
    "nterm_extension",
    "cterm_extension",
    "pseudogene",
    "cds_outframe",
    "antisense",
)
ASSIGN = {"cds_inframe", "sorf", "nterm_extension", "cterm_extension", "pseudogene"}


def _ev() -> Any:
    from giap.genetable import evidence

    return evidence


def _cfg(q: float, relation: str) -> Any:
    from giap.config import GiapConfig

    return GiapConfig(genome_relation=relation, input={"q_value_max": q})


@st.composite
def scenarios(draw: st.DrawFn) -> dict[str, Any]:
    n_genes = draw(st.integers(1, 4))
    genes = []
    for i in range(n_genes):
        start = 100 + 600 * i
        end = start + 300
        crosses = False
        if i == n_genes - 1 and draw(st.booleans()):
            start, end, crosses = L - 150, L + 150, True  # origin-crossing gene
        genes.append(
            eo.gene(
                f"GSIM_{5 * (i + 1):05d}",
                start,
                end,
                partial=draw(st.booleans()),
                crosses_origin=crosses,
            )
        )
    n_peps = draw(st.integers(0, 7))
    hits: list[dict[str, Any]] = []
    peps: list[dict[str, Any]] = []
    for k in range(n_peps):
        name = "PEP" + eo.letters(k) + "K"
        on = draw(st.lists(st.integers(0, n_genes - 1), min_size=1, max_size=2, unique=True))
        uniq = "unique" if len(on) == 1 else draw(st.sampled_from(["multi_identical", "shared"]))
        q = draw(st.sampled_from([None, 0.0005, 0.004, 0.01, 0.02, 0.2]))
        peps.append(eo.pep(name, q))
        for gi in on:
            g = genes[gi]
            cls = draw(st.sampled_from(CLASSES))
            n_aa = draw(st.integers(7, 12))
            off = 3 * draw(st.integers(0, 95))
            s = int(g["start"]) + off
            e = s + 3 * n_aa
            if s >= L:
                s, e = s - L, e - L
            route = draw(st.sampled_from(["exact", "exact", "tolerant"]))
            hits.append(
                eo.hit(
                    name,
                    str(g["locus_tag"]) if cls in ASSIGN or uniq == "shared" else None,
                    s,
                    e,
                    locus_class=cls if uniq != "shared" else "cds_inframe",
                    uniqueness=uniq,
                    route=route,
                    placement_q=draw(st.sampled_from([0.0, 0.004, 0.01, 0.03, None])),
                    flags=draw(st.sampled_from([[], ["same_db"]])),
                    crosses_origin=e > L,
                )
            )
    return {
        "genes": genes,
        "hits": hits,
        "peps": peps,
        "q": draw(st.sampled_from([0.01, 0.05])),
        "relation": draw(st.sampled_from(["related", "same"])),
        "min_peptides": draw(st.integers(1, 3)),
    }


def _run(sc: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    ev = _ev()
    hits = eo.hits_frame(sc["hits"])
    acc = ev.accepted_peptides(hits, eo.peptides_frame(sc["peps"]), _cfg(sc["q"], sc["relation"]))
    out = ev.gene_evidence(
        eo.genes_frame(sc["genes"]),
        hits,
        acc,
        eo.replicons_frame((("chr1", L, "circular"),)),
        min_peptides=sc["min_peptides"],
    )
    return acc, out


def _rows(out: pd.DataFrame, sc: dict[str, Any]) -> dict[str, dict[str, Any]]:
    tags = [str(g["locus_tag"]) for g in sc["genes"]]
    assert len(out) == len(tags)
    if "locus_tag" in out.columns:
        return {str(r["locus_tag"]): r for r in out.to_dict("records")}
    return dict(zip(tags, out.to_dict("records"), strict=True))


def _flags(cell: Any) -> frozenset[str]:
    if cell is None or (not isinstance(cell, (list, tuple, np.ndarray)) and pd.isna(cell)):
        return frozenset()
    return frozenset(str(f) for f in cell)


def _summary(out: pd.DataFrame, sc: dict[str, Any]) -> dict[str, tuple[Any, ...]]:
    res = {}
    for tag, r in _rows(out, sc).items():
        bq = r["best_peptide_q"]
        bq = None if bq is None or pd.isna(bq) else round(float(bq), 12)
        res[tag] = (
            *(int(r[c]) for c in COUNT_COLS),
            round(float(r["coverage_frac"]), 12),
            bq,
            bool(r["detected_2pep"]),
            _flags(r["evidence_flags"]),
        )
    return res


@SETTINGS
@given(scenarios())
def test_matches_oracle(sc: dict[str, Any]) -> None:
    _, out = _run(sc)
    rules = eo.Rules(sc["q"], 0.01, sc["relation"])
    exp = eo.gene_evidence(
        eo.genes_frame(sc["genes"]),
        eo.hits_frame(sc["hits"]),
        eo.peptides_frame(sc["peps"]),
        eo.replicons_frame((("chr1", L, "circular"),)),
        rules,
        sc["min_peptides"],
    )
    for tag, r in _rows(out, sc).items():
        e = exp[tag]
        for c in COUNT_COLS:
            assert int(r[c]) == e[c], (tag, c)
        assert math.isclose(float(r["coverage_frac"]), e["coverage_frac"], abs_tol=1e-12), tag
        assert bool(r["detected_2pep"]) == e["detected_2pep"], tag
        assert _flags(r["evidence_flags"]) == e["evidence_flags"], tag
        bq = r["best_peptide_q"]
        if math.isnan(e["best_peptide_q"]):
            assert bq is None or pd.isna(bq)
        else:
            assert math.isclose(float(bq), e["best_peptide_q"]), tag


@SETTINGS
@given(scenarios(), st.randoms(use_true_random=False))
def test_relabelling_peptides_changes_nothing(sc: dict[str, Any], rnd: Any) -> None:
    names = sorted({p["peptide"] for p in sc["peps"]})
    new = ["REL" + eo.letters(i) + "AK" for i in range(len(names))]
    rnd.shuffle(new)
    mapping = dict(zip(names, new, strict=True))
    sc2 = {
        **sc,
        "hits": [{**h, "peptide": mapping[h["peptide"]]} for h in sc["hits"]],
        "peps": [{**p, "peptide": mapping[p["peptide"]]} for p in sc["peps"]],
    }
    rnd.shuffle(sc2["hits"])
    rnd.shuffle(sc2["peps"])
    assert _summary(_run(sc)[1], sc) == _summary(_run(sc2)[1], sc2)


@SETTINGS
@given(scenarios(), st.data())
def test_duplicate_c11_row_changes_nothing(sc: dict[str, Any], data: st.DataObject) -> None:
    if not sc["hits"]:
        return
    i = data.draw(st.integers(0, len(sc["hits"]) - 1))
    sc2 = {**sc, "hits": [*sc["hits"], dict(sc["hits"][i])]}
    acc1, out1 = _run(sc)
    acc2, out2 = _run(sc2)
    assert _summary(out1, sc) == _summary(out2, sc2)
    key = ["gene", "peptide"]
    assert len(acc1) == len(acc2)
    assert sorted(map(tuple, acc1[key].astype(str).values.tolist())) == sorted(
        map(tuple, acc2[key].astype(str).values.tolist())
    )


@SETTINGS
@given(scenarios())
def test_count_invariants(sc: dict[str, Any]) -> None:
    acc, out = _run(sc)
    for tag, r in _rows(out, sc).items():
        n_acc, n_uni, n_multi = (
            int(r["n_accepted"]),
            int(r["n_peptides_unique"]),
            int(r["n_peptides_multi_identical"]),
        )
        assert n_acc >= n_uni
        assert n_acc == n_uni + n_multi  # uniqueness is constant per peptide here
        assert int(r["n_peptides_tolerant"]) <= n_acc
        cov = float(r["coverage_frac"])
        assert 0.0 <= cov <= 1.0
        if n_acc == 0:
            assert cov == 0.0
            assert r["best_peptide_q"] is None or pd.isna(r["best_peptide_q"])
            assert not bool(r["detected_2pep"])
        if bool(r["detected_2pep"]):
            assert n_acc >= sc["min_peptides"]
        assert n_acc == int((acc["gene"].astype(str) == tag).sum())
    # accepted pairs are unique and only unique/multi_identical
    assert not acc.duplicated(["gene", "peptide"]).any()
    assert set(acc["uniqueness"]) <= {"unique", "multi_identical"}
    elig = acc["abundance_eligible"].astype(bool)
    assert not (elig & (acc["uniqueness"] != "unique")).any()
    assert not (elig & (acc["locus_class"] == "nterm_extension")).any()


@SETTINGS
@given(scenarios())
def test_tighter_q_never_accepts_more(sc: dict[str, Any]) -> None:
    loose = _summary(_run({**sc, "q": 0.05})[1], sc)
    tight = _summary(_run({**sc, "q": 0.01})[1], sc)
    for tag in loose:
        assert tight[tag][4] <= loose[tag][4]  # n_accepted
        assert tight[tag][5] <= loose[tag][5] + 1e-12  # coverage
