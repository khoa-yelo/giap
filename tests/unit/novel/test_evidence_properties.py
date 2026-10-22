"""Properties of the M7-09 evidence and ranking (docs/specs/M7-09-novel-evidence-and-ranking.md
acceptance criteria 2 and 3; rules 2, 9, 10):

* ``called`` (and ``rbs_supported``) never depend on ``promoter_score`` or the terminator layer:
  permuting the promoter scores across ORFs and redrawing the terminator calls changes neither;
* ``rank`` is a strict total order on the called ORFs (exactly ``1..k``, following ``rank_key``),
  NA on the others, and does not depend on the input row order (table rows and ``orfs`` list
  shuffled together);
* ``rbs_threshold`` is the ``numpy.quantile`` (linear) of the eligible detected scores, ``None``
  below 20 of them, and does not depend on the gene-table row order.

The genome has one example-1 ORF per contig (random peptide codons, so no peptide maps twice);
hypothesis draws each ORF's class, peptide count, chosen start, TIR, q and promoter score.
"""

from __future__ import annotations

import random
from typing import Any

import numpy as np
import pandas as pd
from helpers.novel_evidence_fixture import cfg, ev, gene_table, is_na, term
from helpers.novel_oracle import STOPS, Case, il, translate
from helpers.novel_outputs_fixture import ex1_hits, ex1_seq, orfs_of, scored, table
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.model.coords import Replicon

K = 6
CONTIGS = [f"k{i}" for i in range(K)]
_POOL = sorted(a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS)
SETTINGS = settings(
    max_examples=40, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
CLASSES = ("intergenic", "cds_outframe", "antisense", "sorf")


def _codons() -> dict[str, list[str]]:
    rnd = random.Random(909)
    out: dict[str, list[str]] = {}
    seen: set[str] = set()
    for c in CONTIGS:
        while True:
            cod = [rnd.choice(_POOL) for _ in range(13)]
            seq = ex1_seq(codons=cod).text()
            peps = {il(translate(seq[72:102])), il(translate(seq[90:111]))}
            if not peps & seen:
                seen |= peps
                out[c] = cod
                break
    return out


CODONS = _codons()

orf_spec = st.fixed_dictionaries({
    "cls": st.sampled_from(CLASSES),
    "n_hits": st.sampled_from([1, 2, 2]),
    "best": st.sampled_from([None, 42, 60]),
    "tir": st.one_of(st.none(), st.sampled_from([100.0, 390.0, 500.0, 2000.0])),
    "q": st.sampled_from([1e-5, 1e-3, 0.01]),
    "promoter": st.one_of(st.none(), st.sampled_from([0.1, 0.5, 0.9])),
    "term": st.one_of(st.none(), st.sampled_from([0.2, 0.5, 0.8])),
})  # fmt: skip
specs_st = st.lists(orf_spec, min_size=K, max_size=K)


def _case(specs: list[dict[str, Any]]) -> Case:
    genome, hits, q = {}, [], {}
    for c, s in zip(CONTIGS, specs, strict=True):
        seq = ex1_seq(codons=CODONS[c]).text()
        genome[c] = seq
        hs = ex1_hits(seq, contig=c, cls=s["cls"])[: s["n_hits"]]
        hits += hs
        q.update({h["peptide"]: s["q"] for h in hs})
    reps = {c: Replicon(c, 1000, "circular") for c in CONTIGS}
    return Case(genome, reps, hits, q=q)


def _run(
    case: Case,
    specs: list[dict[str, Any]],
    *,
    promoter: list[float | None] | None = None,
    terms: list[float | None] | None = None,
    perm: list[int] | None = None,
) -> pd.DataFrame:
    orfs = orfs_of(case)
    assert [o.contig for o in orfs] == CONTIGS
    ids = [o.orf_id for o in orfs]
    prom = promoter if promoter is not None else [s["promoter"] for s in specs]
    tconf = terms if terms is not None else [s["term"] for s in specs]
    sc = scored(
        orfs,
        best={i: s["best"] for i, s in zip(ids, specs, strict=True)},
        ostir_tir={i: (s["tir"] if s["best"] is not None else None)
                   for i, s in zip(ids, specs, strict=True)},
        promoter_score=dict(zip(ids, prom, strict=True)),
    )  # fmt: skip
    t = table(case, orfs, sc)
    calls = [term(c, 199, 230, "+", v) for c, v in zip(CONTIGS, tconf, strict=True) if v]
    if perm is not None:
        t = t.iloc[perm].reset_index(drop=True)
        orfs = [orfs[k] for k in perm]
    return ev().apply_evidence(
        t, orfs, gene_table=gene_table([100.0 * k for k in range(1, 31)]), terminators=calls,
        gm=case.gm, genome=case.genome, cfg=cfg(),
    )  # fmt: skip


def _by_id(t: pd.DataFrame, col: str) -> dict[str, Any]:
    return {i: (None if is_na(v) else v) for i, v in zip(t["orf_id"], t[col], strict=True)}


@SETTINGS
@given(specs=specs_st, data=st.data())
def test_called_never_depends_on_promoter_or_terminator(specs: list[dict], data: Any) -> None:
    case = _case(specs)
    base = _run(case, specs)
    prom = data.draw(st.permutations([s["promoter"] for s in specs]))
    terms = data.draw(st.lists(st.one_of(st.none(), st.floats(0, 1)), min_size=K, max_size=K))
    other = _run(case, specs, promoter=list(prom), terms=terms)
    assert _by_id(base, "called") == _by_id(other, "called")
    assert _by_id(base, "rbs_supported") == _by_id(other, "rbs_supported")
    none = _run(case, specs, promoter=[None] * K, terms=[None] * K)
    assert _by_id(base, "called") == _by_id(none, "called")


@SETTINGS
@given(specs=specs_st, data=st.data())
def test_rank_is_a_strict_total_order_on_called_and_row_order_free(
    specs: list[dict], data: Any
) -> None:
    case = _case(specs)
    t = _run(case, specs)
    recs = t.to_dict("records")
    called = [r for r in recs if r["called"]]
    assert sorted(int(r["rank"]) for r in called) == list(range(1, len(called) + 1))
    assert all(is_na(r["rank"]) for r in recs if not r["called"])
    want = [r["orf_id"] for r in sorted(called, key=ev().rank_key)]
    assert [r["orf_id"] for r in sorted(called, key=lambda r: int(r["rank"]))] == want
    perm = data.draw(st.permutations(list(range(K))))
    shuffled = _run(case, specs, perm=list(perm))
    assert list(shuffled["orf_id"]) == [recs[k]["orf_id"] for k in perm]
    for col in ("rank", "called", "flags", "rbs_supported", "terminator_conf"):
        assert _by_id(shuffled, col) == _by_id(t, col), col


@SETTINGS
@given(
    tirs=st.lists(st.floats(0, 1e5, allow_nan=False), min_size=0, max_size=60),
    noise=st.lists(
        st.tuples(
            st.one_of(st.none(), st.floats(0, 1e5)),
            st.sampled_from([True, False, None]),
            st.booleans(),
        ),
        max_size=20,
    ),
    q=st.floats(0, 1),
    seed=st.integers(0, 2**16),
)
def test_rbs_threshold_property(
    tirs: list[float], noise: list[tuple[Any, Any, bool]], q: float, seed: int
) -> None:
    gt = gene_table(tirs, extra=noise)
    elig = list(tirs) + [t for t, d, p in noise if t is not None and d is True and not p]
    got = ev().rbs_threshold(gt, q)
    if len(elig) < 20:
        assert got is None
    else:
        assert got is not None
        assert np.isclose(got, float(np.quantile(np.array(elig), q)), rtol=1e-12, atol=1e-9)
    shuffled = gt.sample(frac=1.0, random_state=seed).reset_index(drop=True)
    again = ev().rbs_threshold(shuffled, q)
    assert (got is None and again is None) or np.isclose(got, again)  # type: ignore[arg-type]
