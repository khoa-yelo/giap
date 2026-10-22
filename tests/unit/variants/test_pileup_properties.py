"""M7-02 properties (docs/specs/M7-02-pileup-and-calls.md, acceptance criteria 2-4).

* Random hand-built alignments on G1/G2 of ``helpers.pileup_oracle`` (I/L-free stretches, so
  folding is the identity there) with random sources, ``search_db``, placement q, uniqueness and
  genome relation: the pileup equals the rule 1-7 oracle and the sites equal the rows derived
  from it (rule 7 confidence, rule 11 columns).
* Permuting the input rows leaves both frames byte-identical.
* Adding a REF-only peptide never turns ``REF`` into anything else; adding an ALT peptide never
  lowers a confidence.
* Same-strain: peptides cut from the genome's own genes (both strands, across the origin, I/L
  flipped at random) and aligned by M7-01 ``align_hit``: every covered codon is ``REF``, or
  ``no_call`` where the genome residue is I/L and a covering peptide reads the other letter;
  no site rows.

Not statistical (M7-02 "counting rules"): no calibration test is owed here; the F1 gate is M7-13.
"""

from __future__ import annotations

import importlib
from collections.abc import Mapping, Sequence
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
from helpers.align_oracle import place
from helpers.fs_genomes import background, codons_for
from helpers.pileup_oracle import (
    G1,
    G2,
    SEQ,
    GeneSpec,
    L,
    Sup,
    alt_counts_text,
    expected,
    gene_hit,
)
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

from giap.config.model import GiapConfig
from giap.model import Replicon
from giap.model.coords import Interval
from giap.model.genes import Gene, GeneModel

REP = Replicon("c1", L, "circular")
SETTINGS = settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
ALT_LETTERS = "DEKQWY"
WINDOWS = {"G1": (15, 60), "G2": (10, 60)}
GENES = {"G1": G1, "G2": G2}


def pl() -> ModuleType:
    return importlib.import_module("giap.variants.pileup")


def al() -> ModuleType:
    return importlib.import_module("giap.variants.align")


def _gene(spec: GeneSpec) -> Gene:
    iv = Interval("c1", spec.start, spec.end, spec.strand)  # type: ignore[arg-type]
    return Gene(spec.tag, "cds", iv, "M" + spec.protein[1:], None, "ATG")


GM = GeneModel({"c1": REP}, (_gene(G1), _gene(G2)), ())


def run(
    aligns: Sequence[Any],
    support: Mapping[str, Sup],
    relation: str = "related",
    high: int = 2,
    gm: GeneModel = GM,
    genome: Mapping[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    m = pl()
    cfg = GiapConfig.model_validate(
        {"genome_relation": relation, "variants": {"high_min_peptides": high}}
    )
    s = {p: m.Support(p, v.source, v.search_db, v.q, v.uniqueness) for p, v in support.items()}
    el = m.eligible(list(aligns), s, cfg)
    pile = m.pileup(el, gm, genome or {"c1": SEQ}, cfg, support=s)
    return pile, m.sites(pile, cfg)


# --------------------------------------------------------------------------- strategies


@st.composite
def hit_specs(draw: st.DrawFn) -> tuple[str, int, int, dict[int, str]]:
    tag = draw(st.sampled_from(sorted(GENES)))
    lo, hi = WINDOWS[tag]
    a = draw(st.integers(lo, hi - 4))
    b = draw(st.integers(a + 4, min(a + 10, hi)))
    g = GENES[tag]
    positions = draw(st.lists(st.integers(a, b - 1), max_size=2, unique=True))
    reads = {}
    for i in positions:
        letters = [x for x in ALT_LETTERS if x != g.protein[i]]
        reads[i] = draw(st.sampled_from(letters))
    return tag, a, b, reads


sup_st = st.builds(
    Sup,
    source=st.sampled_from(["denovo", "database"]),
    search_db=st.sampled_from(["same", "other"]),
    q=st.sampled_from([0.0005, 0.01, 0.02]),
    uniqueness=st.sampled_from(["unique", "unique", "multi_identical", "shared"]),
)


@st.composite
def scenarios(draw: st.DrawFn) -> tuple[list[Any], dict[str, Sup], str, int]:
    a = al()
    specs = draw(st.lists(hit_specs(), min_size=1, max_size=12))
    hits = [gene_hit(a, GENES[t], lo, hi, r) for t, lo, hi, r in specs]
    support = {p: draw(sup_st) for p in sorted({h.peptide for h in hits})}
    relation = draw(st.sampled_from(["related", "same"]))
    high = draw(st.sampled_from([2, 3]))
    return hits, support, relation, high


def oracle_or_reject(
    hits: Sequence[Any], support: Mapping[str, Sup], relation: str
) -> dict[tuple[str, str, int], dict[str, Any]]:
    try:
        return expected(hits, support, relation)
    except AssertionError:
        assume(False)
        raise


# --------------------------------------------------------------------------- oracle agreement


def _alt_type(key: str) -> tuple[str, str]:
    return ("substitution", key) if ">" not in key else ("block", key.split(">")[1])


@SETTINGS
@given(scenarios())
def test_pileup_and_sites_match_the_oracle(sc: tuple[list[Any], dict[str, Sup], str, int]) -> None:
    hits, support, relation, high = sc
    want = oracle_or_reject(hits, support, relation)
    pile, sites = run(hits, support, relation, high)
    got = {(r["contig"], r["strand"], int(r["codon_start"])): r for _, r in pile.iterrows()}
    assert set(got) == set(want)
    for key, w in want.items():
        r = got[key]
        assert (int(r["n_ref"]), int(r["n_il"]), r["call"]) == (w["n_ref"], w["n_il"], w["call"])
        text = r["alt_counts"]
        text = "" if not isinstance(text, str) else text
        assert text == alt_counts_text(w["alts"])
    want_sites = sorted(
        (key[2], key[1], alt, n, w["n_ref"], w["call"], "high" if n >= high else "candidate")
        for key, w in want.items()
        if w["call"] in ("ALT", "mixed")
        for alt, n in w["alts"].items()
    )
    got_sites = sorted(
        (int(s["codon_start"]), s["strand"], s["alt_aa"], int(s["n_alt"]), int(s["n_ref"]),
         s["call"], s["confidence"])
        for _, s in sites.iterrows()
    )  # fmt: skip
    assert got_sites == want_sites
    assert sites["variant_id"].tolist() == [f"VAR_{i}" for i in range(1, len(sites) + 1)]
    assert set(sites["alt_type"]) <= {"substitution"}


# --------------------------------------------------------------------------- invariants


@SETTINGS
@given(scenarios(), st.randoms(use_true_random=False))
def test_permuting_rows_is_byte_identical(
    sc: tuple[list[Any], dict[str, Sup], str, int], rnd: Any
) -> None:
    hits, support, relation, high = sc
    pile1, sites1 = run(hits, support, relation, high)
    shuffled = list(hits)
    rnd.shuffle(shuffled)
    support_rev = dict(reversed(list(support.items())))
    pile2, sites2 = run(shuffled, support_rev, relation, high)
    assert pile1.to_csv(index=False) == pile2.to_csv(index=False)
    assert sites1.to_csv(index=False) == sites2.to_csv(index=False)
    pd.testing.assert_frame_equal(pile1, pile2)
    pd.testing.assert_frame_equal(sites1, sites2)


@SETTINGS
@given(scenarios(), hit_specs())
def test_adding_a_ref_only_peptide_keeps_ref(
    sc: tuple[list[Any], dict[str, Sup], str, int],
    extra: tuple[str, int, int, dict[int, str]],
) -> None:
    hits, support, relation, high = sc
    tag, a, b, _ = extra
    new = gene_hit(al(), GENES[tag], a, b)
    support2 = {**support, new.peptide: support.get(new.peptide, Sup())}
    oracle_or_reject([*hits, new], support2, relation)
    pile1, _ = run(hits, support2, relation, high)
    pile2, _ = run([*hits, new], support2, relation, high)
    before = {
        (r["contig"], r["strand"], int(r["codon_start"])): r["call"] for _, r in pile1.iterrows()
    }
    after = {
        (r["contig"], r["strand"], int(r["codon_start"])): r["call"] for _, r in pile2.iterrows()
    }
    for key, call in before.items():
        if call == "REF":
            assert after[key] == "REF"


@SETTINGS
@given(scenarios(), hit_specs())
def test_adding_an_alt_peptide_never_lowers_a_confidence(
    sc: tuple[list[Any], dict[str, Sup], str, int],
    extra: tuple[str, int, int, dict[int, str]],
) -> None:
    hits, support, relation, high = sc
    tag, a, b, reads = extra
    assume(reads)
    new = gene_hit(al(), GENES[tag], a, b, reads)
    support2 = {**support, new.peptide: support.get(new.peptide, Sup())}
    oracle_or_reject([*hits, new], support2, relation)
    _, s1 = run(hits, support2, relation, high)
    _, s2 = run([*hits, new], support2, relation, high)
    rank = {"candidate": 0, "high": 1}

    def conf(df: pd.DataFrame) -> dict[tuple[Any, ...], int]:
        return {
            (s["contig"], int(s["codon_start"]), s["strand"], s["alt_type"], s["alt_aa"]): rank[
                s["confidence"]
            ]
            for _, s in df.iterrows()
        }

    c1, c2 = conf(s1), conf(s2)
    for key, v in c1.items():
        assert key in c2 and c2[key] >= v


# --------------------------------------------------------------------------- same strain


@st.composite
def own_genome(draw: st.DrawFn) -> tuple[str, GeneModel, list[Any]]:
    """Three genes in disjoint 300-nt windows (rotated, so one may cross the origin) and
    peptides cut from their own translations, some I/L flipped."""
    seed = draw(st.integers(0, 2**31 - 1))
    rng = np.random.default_rng(seed)
    n = 900
    rot = draw(st.integers(0, n - 1))
    seq = background(n, seed=seed % 1000)
    rep = Replicon("c1", n, "circular")
    letters = "ACDEFGHIKLMNPQRSTVWYLI"
    genes, hits = [], []
    a_mod = al()
    for k in range(3):
        n_cod = draw(st.integers(20, 90))
        strand = draw(st.sampled_from(["+", "-"]))
        gs = (300 * k + rot) % n
        ge = gs + 3 * n_cod
        protein = "".join(letters[int(i)] for i in rng.integers(len(letters), size=n_cod))
        seq = place(seq, gs, "".join(codons_for(protein, rng)), strand)
        genes.append(
            Gene(f"T{k}", "cds", Interval("c1", gs, ge, strand), protein, None, "ATG")  # type: ignore[arg-type]
        )
    for g in genes:
        gs, ge, strand = g.interval.start, g.interval.end, g.interval.strand
        protein = g.aa
        for _ in range(draw(st.integers(1, 5))):
            a = draw(st.integers(0, len(protein) - 5))
            b = draw(st.integers(a + 5, min(a + 15, len(protein))))
            pep = list(protein[a:b])
            for j, ch in enumerate(pep):
                if ch in "IL" and draw(st.booleans()):
                    pep[j] = "L" if ch == "I" else "I"
            s, e = (gs + 3 * a, gs + 3 * b) if strand == "+" else (ge - 3 * b, ge - 3 * a)
            if s >= n:
                s, e = s - n, e - n
            hits.append(a_mod.align_hit("".join(pep), "c1", s, e, strand, 0, None, seq, rep))
    genes.sort(key=lambda x: (x.interval.start, x.locus_tag))
    return seq, GeneModel({"c1": rep}, tuple(genes), ()), hits


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(own_genome())
def test_same_strain_peptides_give_only_ref_or_il_no_call(
    case: tuple[str, GeneModel, list[Any]],
) -> None:
    seq, gm, hits = case
    support = {h.peptide: Sup() for h in hits}
    pile, sites = run(hits, support, "same", gm=gm, genome={"c1": seq})
    assert len(sites) == 0
    assert set(pile["call"]) <= {"REF", "no_call"}
    reads: dict[tuple[str, int], set[tuple[str, str]]] = {}
    for h in hits:
        for op in h.ops:
            reads.setdefault((h.strand, op.codons[0] % len(seq)), set()).add((op.ref, op.alt))
    for _, r in pile.iterrows():
        assert empty_text(r["alt_counts"])
        if r["call"] == "no_call":
            obs = reads[(r["strand"], int(r["codon_start"]))]
            assert r["ref_aa"] in ("I", "L")
            assert any(alt != ref for ref, alt in obs)


def empty_text(v: object) -> bool:
    return not isinstance(v, str) or v == ""
