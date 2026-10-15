"""giapsim evidence: tests from the lead decisions after code review of PR #80 (issue M4-09).

docs/specs/M4-09-sim-evidence.md, "Lead decisions after code review", items 4-5:

4. Kill the surviving mutants: origin wrap in ``_span`` (``origin_crossing=True``, seed 1,
   GSIM_05000 on the minus strand, every hit against ``span_coords`` mod length and
   ``crosses_origin``); an sORF fixture (sorf -> ``locus_class = sorf``, not ``cds_inframe``);
   mean/SD of the C10 value minus the true log2 abundance; the ``multi_identical`` partner
   offset; ``_digest`` against the oracle's ``n_obs`` (the false-target / decoy weights).
5. ``EvidenceConfig`` rejects ``mu1 <= mu0`` (PEP/q must be monotone in the score).

Expectations come from tests/helpers/sim_evidence_oracle.py (written from the spec).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import sim_evidence_oracle as eo


def _layout(seed: int, **kw: Any) -> Any:
    from giapsim.layout import LayoutConfig, simulate_layout

    return simulate_layout(LayoutConfig(**kw), np.random.default_rng(seed))


def _express(lay: Any, seed: int) -> Any:
    from giapsim.expression import ExpressionConfig, simulate_expression

    return simulate_expression(lay, ExpressionConfig(), np.random.default_rng(seed))


def _force_detected(ex: Any, tags: set[str] | None = None) -> Any:
    """Expression with ``detected_true`` set for ``tags`` (all genes if None); input untouched."""
    from giapsim.expression import Expression

    g = ex.genes
    det = g["detected_true"].astype(bool) | (
        True if tags is None else g["locus_tag"].astype(str).isin(tags)
    )
    return Expression(g.assign(detected_true=det), ex.islands)


def _evidence(lay: Any, ex: Any, seed: int, **kw: Any) -> Any:
    from giapsim.evidence import EvidenceConfig, simulate_evidence

    return simulate_evidence(lay, ex, EvidenceConfig(**kw), np.random.default_rng(seed))


def _wrap(s: int, e: int, length: int) -> tuple[int, int, bool]:
    """Oracle: shift a span starting past the origin back by one replicon length; a span crosses
    the origin when it starts before and ends strictly after it."""
    shift = length * (s // length)
    return s - shift, e - shift, s - shift < length < e - shift


# ----------------------------------------------------------------------------- item 4: origin wrap


ORIGIN_GENE = "GSIM_05000"


@pytest.fixture(scope="module")
def origin_case() -> tuple[Any, Any]:
    lay = _layout(1, origin_crossing=True)
    ex = _force_detected(_express(lay, 1), {ORIGIN_GENE})
    return lay, _evidence(lay, ex, 1)


def test_origin_gene_is_the_minus_strand_crossing_gene(origin_case: tuple[Any, Any]) -> None:
    lay, _ = origin_case
    g = eo.gene_rows(lay)[ORIGIN_GENE]
    length = int(
        dict(zip(lay.replicons["contig"], lay.replicons["length"], strict=True))[g["contig"]]
    )
    assert g["strand"] == "-"
    assert int(g["start"]) < length < int(g["end"])


def test_origin_hits_wrap_mod_length_and_flag_crossing(origin_case: tuple[Any, Any]) -> None:
    """Every hit on the origin-crossing gene: span_coords reduced mod the replicon length, and
    crosses_origin exactly when the span straddles the origin."""
    lay, ev = origin_case
    g = eo.gene_rows(lay)[ORIGIN_GENE]
    length = int(
        dict(zip(lay.replicons["contig"], lay.replicons["length"], strict=True))[g["contig"]]
    )
    h = ev.peptide_hits
    rows = h[h["gene"].astype(str) == ORIGIN_GENE].to_dict("records")
    assert len(rows) >= 3
    kinds = set()
    for r in rows:
        p = str(r["peptide"])
        offs = eo.offsets(str(g["aa"]), p)
        assert offs, r
        want = {_wrap(*eo.span_coords(int(g["start"]), int(g["end"]), "-", a, a + len(p)), length)
                for a in offs}  # fmt: skip
        got = (int(r["start"]), int(r["end"]), bool(r["crosses_origin"]))
        assert got in want, (got, want)
        assert 0 <= got[0] < length
        kinds.add(
            "cross" if got[2] else ("after" if got[1] <= int(g["end"]) - length else "before")
        )
    # the fixture exercises spans past the origin, across it and before it
    assert kinds == {"cross", "after", "before"}, kinds


@pytest.mark.parametrize(
    ("start", "end", "strand", "a", "b", "length"),
    [
        (985, 1045, "-", 10, 15, 1000),  # [1000, 1015): starts exactly at the origin -> [0, 15)
        (985, 1045, "-", 15, 20, 1000),  # [985, 1000): ends exactly at the origin, not crossing
        (985, 1045, "-", 14, 19, 1000),  # [988, 1003): crosses
        (985, 1045, "+", 0, 3, 1000),  # [985, 994): before the origin
        (985, 1045, "+", 1, 5, 1000),  # [988, 1000): ends exactly at the origin, not crossing
        (985, 1045, "+", 4, 7, 1000),  # [997, 1006): crosses
        (985, 1045, "+", 6, 9, 1000),  # [1003, 1012): wholly past -> [3, 12)
    ],
)
def test_span_boundaries(start: int, end: int, strand: str, a: int, b: int, length: int) -> None:
    """``_span`` at the origin boundaries (an end exactly at the origin does not cross)."""
    from giapsim.evidence import _span

    g = {"start": start, "end": end, "strand": strand}
    want = _wrap(*eo.span_coords(start, end, strand, a, b), length)
    assert _span(g, a, b, length) == want


# ----------------------------------------------------------------------------- item 4: sORF class


def test_sorf_peptides_have_locus_class_sorf() -> None:
    lay = _layout(7, n_genes_chromosome=300, frac_sorf=0.6)
    ev = _evidence(lay, _force_detected(_express(lay, 7)), 7)
    genes = eo.gene_rows(lay)
    h = ev.peptide_hits
    ft = [str(genes[str(t)]["feature_type"]) for t in h["gene"]]
    assert ft.count("sorf") >= 5
    for f, c in zip(ft, h["locus_class"], strict=True):
        assert c == eo.LOCUS_CLASS[f], (f, c)
    assert eo.LOCUS_CLASS["sorf"] == "sorf"


# ----------------------------------------------------------------------------- item 4: C10 noise


def test_c10_value_minus_true_abundance_centred_and_small() -> None:
    """C10 value = true log2 abundance + noise: mean difference 0 (within 4 SE), SD small."""
    lay = _layout(1)
    ex = _express(lay, 1)
    ev = _evidence(lay, ex, 1, abundance_table="ibaq")
    ab = ev.abundance
    assert ab is not None and len(ab) >= 100
    true = dict(
        zip(ex.genes["locus_tag"].astype(str), ex.genes["abundance_log2_true"], strict=True)
    )
    d = ab["value"].to_numpy(float) - np.array([true[str(p)] for p in ab["protein_id"]])
    sd = float(d.std(ddof=1))
    assert 0.05 < sd < 0.5, sd
    assert abs(float(d.mean())) <= 4 * sd / math.sqrt(len(d)), (d.mean(), sd, len(d))


# --------------------------------------------------------------------------- item 4: partner offset


def _offset(g: Any, s: int, e: int) -> tuple[int, int]:
    """Residue span [a, b) of a gene from its 0-based genomic span (no origin crossing)."""
    if g["strand"] == "+":
        a3, b3 = s - int(g["start"]), e - int(g["start"])
    else:
        a3, b3 = int(g["end"]) - e, int(g["end"]) - s
    assert a3 % 3 == 0 and b3 % 3 == 0
    return a3 // 3, b3 // 3


def test_multi_identical_partner_at_same_residue_offset() -> None:
    """Both rows of a multi_identical peptide sit at the same residue offset of their gene, the
    offset where the peptide occurs in its true gene (the copies are identical)."""
    lay = _layout(1)
    ev = _evidence(lay, _express(lay, 1), 1, multi_identical_pairs=2)
    genes = eo.gene_rows(lay)
    gene_true = dict(zip(ev.truth["peptide"], ev.truth["gene_true"], strict=True))
    h = ev.peptide_hits
    mi = h[h["uniqueness"] == "multi_identical"]
    assert mi["peptide"].nunique() >= 4
    n_nonzero = 0
    for p, grp in mi.groupby("peptide"):
        spans = {}
        for r in grp.to_dict("records"):
            g = genes[str(r["gene"])]
            a, b = _offset(g, int(r["start"]), int(r["end"]))
            assert b - a == len(str(p)) and a >= 0 and b <= len(str(g["aa"]))
            spans[str(r["gene"])] = a
        assert len(spans) == 2
        tg = str(gene_true[p])
        assert tg in spans
        assert spans[tg] in eo.offsets(str(genes[tg]["aa"]), str(p))
        assert len(set(spans.values())) == 1, (p, spans)
        n_nonzero += spans[tg] > 0
    assert n_nonzero >= 1


# ----------------------------------------------------------------------------- item 4: weights


@pytest.fixture(scope="module")
def sorf_layout() -> Any:
    return _layout(7, n_genes_chromosome=300, frac_sorf=0.6)


def test_digest_matches_oracle(sorf_layout: Any) -> None:
    """``_digest``: the 7-30 aa, <= 1 missed Trypsin/P peptides at their first offset, and the
    weight n_observable (distinct 0-missed peptides) of the oracle."""
    from giapsim.evidence import _digest

    for aa in [*sorf_layout.genes["aa"], None, "", "K" * 3, "AAAAAAK" + "PAAAAAK"]:
        cands, w = _digest(aa)
        assert w == eo.n_obs(aa), aa
        assert {p for p, _ in cands} == eo.peptides_of(aa), aa
        assert len(cands) == len({p for p, _ in cands})
        for p, a in cands:
            assert a == eo.offsets(str(aa), p)[0]


def test_null_rows_never_on_genes_without_observable_peptides(sorf_layout: Any) -> None:
    """False targets and decoys are drawn proportional to n_observable: a gene whose only 7-30 aa
    peptides need a missed cleavage (n_observable = 0) never receives one."""
    lay = sorf_layout
    aa = dict(zip(lay.genes["locus_tag"].astype(str), lay.genes["aa"], strict=True))
    zero = {t for t, a in aa.items() if eo.n_obs(a) == 0 and eo.digestible(a)}
    assert len(zero) >= 5
    ev = _evidence(lay, _express(lay, 7), 7, pi_wrong=0.5, n_decoys=3000)
    t = ev.truth
    dec = {str(g) for g in t["decoy_gene"] if isinstance(g, str)}
    bad = eo.false_set(t)
    h = ev.peptide_hits
    fal = {str(g) for p, g in zip(h["peptide"], h["gene"], strict=True) if p in bad}
    assert len(dec) >= 100 and len(fal) >= 100
    assert not dec & zero, sorted(dec & zero)
    assert not fal & zero, sorted(fal & zero)


# ----------------------------------------------------------------------------- item 5


@pytest.mark.parametrize("mixture", [(3.0, 0.0), (1.0, 1.0), (0.0, -0.5)])
def test_config_rejects_mu1_not_above_mu0(mixture: tuple[float, float]) -> None:
    from giapsim.evidence import EvidenceConfig

    with pytest.raises(ValueError):
        EvidenceConfig(mixture=mixture)


def test_config_accepts_mu1_above_mu0() -> None:
    from giapsim.evidence import EvidenceConfig

    assert EvidenceConfig(mixture=(-1.0, 0.5)).mixture == (-1.0, 0.5)
