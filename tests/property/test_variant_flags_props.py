"""M7-05 properties (docs/specs/M7-05-variant-flags.md rule 7 and acceptance criterion 3): flags
never change ``call``, ``confidence``, ``n_ref`` or ``n_alt``.

* Site level: ``site_flags`` never mutates its ``sites`` row, is deterministic, returns known
  flags once each in the §9.4 order, and follows rules 2-4 for every alt type.
* Chain level: the M7-04 ``run_variants`` (which now fills ``flags``) and the same chain without
  flags (``sites`` -> ``vcf_records`` with no ``flags`` column) give the same rows with the same
  ``variant_id``, ``call``, ``confidence``, ``n_ref`` and ``n_alt``; the unflagged
  ``vcf_unrepresentable`` flags are kept; under ``related`` every row carries
  ``strain_difference``. Peptide hits are hand-built on the M7-02 test genome
  (``helpers.pileup_oracle``): exact windows of G1 (+) and G2 (-) with at most one substitution.
"""

from __future__ import annotations

import copy
import importlib
from types import ModuleType
from typing import Any

import pandas as pd
from helpers.pileup_oracle import G1, G2, SEQ, GeneSpec, L, codon_aa, codon_nt
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.config.model import GiapConfig
from giap.model import Replicon
from giap.model.coords import Interval
from giap.model.genes import Gene, GeneModel

ORDER = (
    "mass_equivalent", "snv", "multi_nt", "strain_difference", "database_derived",
    "start_codon", "stop_codon",
)  # fmt: skip
VCF_FLAG = "vcf_unrepresentable"
REP = Replicon("c1", L, "circular")
GENOME = {"c1": SEQ}
ALTS = "ACDEFGHIKLMNPQRSTVWY"
SOURCES = st.tuples(st.sampled_from(["denovo", "database"]), st.sampled_from(["same", "other"]))


def fl() -> ModuleType:
    return importlib.import_module("giap.variants.flags")


def _gene(g: GeneSpec, sc: str) -> Gene:
    return Gene(g.tag, "cds", Interval("c1", g.start, g.end, g.strand), "M" + g.protein[1:],  # type: ignore[arg-type]
                None, sc)  # fmt: skip


GM = GeneModel({"c1": REP}, (_gene(G1, "GTG"), _gene(G2, "ATG")), ())


def _cfg(relation: str, tol: float = 0.05, k: int = 2) -> GiapConfig:
    return GiapConfig.model_validate({
        "genome_relation": relation,
        "variants": {"mass_equivalent_tol_da": tol, "high_min_peptides": k},
    })  # fmt: skip


# --------------------------------------------------------------------------- site level


@st.composite
def sites(draw: st.DrawFn) -> dict[str, Any]:
    kind = draw(st.sampled_from(["substitution", "swap", "block", "frameshift"]))
    g = draw(st.sampled_from([G1, G2]))
    i = draw(st.integers(0, 90))
    codon = g.codon(i)
    width = 2 if kind == "swap" or (kind == "block" and draw(st.booleans())) else 1
    ref_codon = "".join(codon_nt(SEQ, g.codon(i + j), g.strand) for j in range(width))
    if kind == "frameshift":
        alt = draw(st.sampled_from(["fs+1", "fs-1"]))
    elif kind == "swap":
        r = [codon_aa(SEQ, g.codon(i + j), g.strand) for j in range(2)]
        alt = r[1] + r[0]
    elif kind == "block":
        alt = draw(st.text(ALTS, min_size=3 - width, max_size=3 - width))
    else:
        alt = draw(st.sampled_from(ALTS.replace("I", "").replace("L", "") + "J"))
    n_ref = draw(st.integers(0, 3))
    n_alt = draw(st.integers(1, 4))
    return {
        "variant_id": "VAR_1", "contig": "c1", "codon_start": codon, "strand": g.strand,
        "gene": g.tag, "aa_pos": i + 1, "ref_codon": ref_codon,
        "ref_aa": codon_aa(SEQ, codon, g.strand), "alt_aa": alt, "alt_type": kind,
        "call": "mixed" if n_ref else "ALT", "confidence": "high" if n_alt >= 2 else "candidate",
        "n_ref": n_ref, "n_alt": n_alt, "source": "denovo",
    }  # fmt: skip


@settings(max_examples=300, deadline=None)
@given(
    s=sites(),
    relation=st.sampled_from(["related", "same"]),
    tol=st.sampled_from([0.0, 0.02, 0.05, 0.1]),
    sources=st.frozensets(SOURCES, max_size=4),
)
def test_site_flags_never_touch_the_call(
    s: dict[str, Any], relation: str, tol: float, sources: frozenset[tuple[str, str]]
) -> None:
    before = copy.deepcopy(s)
    cfg = _cfg(relation, tol)
    got = fl().site_flags(s, GM, cfg, sources)
    assert s == before, "site_flags mutated its site"
    assert fl().site_flags(s, GM, cfg, sources) == got
    assert isinstance(got, list)
    assert set(got) <= set(ORDER), got
    assert got == sorted(set(got), key=ORDER.index), got
    assert ("strain_difference" in got) == (relation == "related")
    if relation == "related":
        assert ("database_derived" in got) == (("database", "other") in sources)
    nt = [f for f in got if f in ("snv", "multi_nt")]
    assert len(nt) == (1 if s["alt_type"] == "substitution" else 0), got
    if s["alt_type"] == "frameshift":
        assert "mass_equivalent" not in got
    if s["alt_type"] == "swap":  # ADR-0049 decision 2: same composition, always flagged
        assert "mass_equivalent" in got
    if s["alt_type"] != "substitution":
        assert "start_codon" not in got
    if s["ref_aa"] != "*":
        assert "stop_codon" not in got


# --------------------------------------------------------------------------- chain level


@st.composite
def hits(draw: st.DrawFn) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: dict[str, dict[str, Any]] = {}
    srcs: dict[str, tuple[str, str]] = {}
    for _ in range(draw(st.integers(1, 6))):
        g = draw(st.sampled_from([G1, G2]))
        n = draw(st.integers(7, 12))
        a = draw(st.integers(0, 99 - n))
        ref = [codon_aa(SEQ, g.codon(i), g.strand) for i in range(a, a + n)]
        pep = list(ref)
        sub = draw(st.none() | st.integers(0, n - 1))
        n_diff = 0
        if sub is not None:
            pep[sub] = draw(st.sampled_from(ALTS))
            if pep[sub] != ref[sub] and not {pep[sub], ref[sub]} <= {"I", "L"}:
                n_diff = 1
        peptide = "".join(pep)
        lo, hi = sorted((g.codon(a), g.codon(a + n - 1)))
        source, sdb = draw(SOURCES)
        rows.setdefault(peptide, {
            "peptide": peptide, "contig": "c1", "start": lo, "end": hi + 3, "strand": g.strand,
            "n_diff": n_diff, "diff_types": "substitution" if n_diff else pd.NA,
            "source": source, "placement_q": draw(st.sampled_from([0.001, 0.005, 0.02])),
            "uniqueness": "unique",
        })  # fmt: skip
        srcs.setdefault(peptide, (source, sdb))
    ph = pd.DataFrame(list(rows.values()))
    peps = pd.DataFrame({"peptide": list(srcs), "search_db": [s[1] for s in srcs.values()],
                         "protein_ids": "X_1"})  # fmt: skip
    return ph, peps


def _unflagged(ph: pd.DataFrame, peps: pd.DataFrame, cfg: GiapConfig) -> pd.DataFrame:
    from giap.steps.variants import support_of
    from giap.variants.align import align_hits
    from giap.variants.pileup import eligible, pileup, sites
    from giap.variants.vcf_records import vcf_records

    support = support_of(ph, peps, {g.locus_tag for g in GM.genes})
    aligns = align_hits(ph, None, GENOME, GM.replicons)
    pile = pileup(eligible(aligns, support, cfg), GM, GENOME, cfg, support=support)
    _, table = vcf_records(sites(pile, cfg), pile, GENOME, GM.replicons)
    return table


CORE = ["variant_id", "contig", "codon_start", "strand", "alt_aa", "alt_type", "call",
        "confidence", "n_ref", "n_alt"]  # fmt: skip


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    data=hits(),
    relation=st.sampled_from(["related", "same"]),
    k=st.integers(1, 3),
    tol=st.sampled_from([0.0, 0.05]),
)
def test_flags_never_change_call_confidence_or_counts(
    data: tuple[pd.DataFrame, pd.DataFrame], relation: str, k: int, tol: float
) -> None:
    from giap.steps.variants import run_variants

    ph, peps = data
    cfg = _cfg(relation, tol, k)
    base = _unflagged(ph, peps, cfg)
    _, table, _ = run_variants(ph, None, peps, GM, GENOME, cfg)
    assert len(table) == len(base)
    pd.testing.assert_frame_equal(
        table[CORE].reset_index(drop=True), base[CORE].reset_index(drop=True)
    )
    for got, old in zip(table["flags"], base["flags"], strict=True):
        got, old = tuple(got), tuple(old)
        assert set(old) <= set(got)
        assert (VCF_FLAG in got) == (VCF_FLAG in old)
        m705 = [f for f in got if f != VCF_FLAG]
        assert m705 == sorted(set(m705), key=ORDER.index), got
        if relation == "related":
            assert "strain_difference" in got, got
        else:
            assert "strain_difference" not in got, got
