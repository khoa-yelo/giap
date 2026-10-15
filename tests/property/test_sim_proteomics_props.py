"""Properties of giapsim's proteomics module (issue M3-19; GIAP_SPEC.md §11.2-§11.3).

* ``digest`` equals the independent Trypsin/P oracle (tests/helpers/peptide_oracle.py) for any
  protein, length range and number of missed cleavages;
* for random seeds and feature mixes, the simulated truth agrees with the genome: every correct
  peptide is the translation of its stated span inside its stated true protein, wrong peptides
  occur in no protein of the genome (up to I/L), decoys are reversed tryptic peptides of the
  named protein, PEP/q follow the score mixture, q is monotone, and the table passes M1-06.

Conventions as in tests/unit/test_sim_proteomics.py.
"""

from __future__ import annotations

import itertools
import tempfile
from pathlib import Path
from typing import Any

import pytest
from helpers.peptide_oracle import (
    Run,
    empty,
    il,
    maps_to,
    mixture_fdr,
    mixture_pep,
    truthy,
    trypsin_p,
    trypsin_p_peptides,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

import giapsim
from giap.config.model import InputConfig
from giap.io.peptides import read_peptide_table, validate_peptides

AA = "ACDEFGHIKLMNPQRSTVWY"


@settings(max_examples=300, deadline=None)
@given(
    protein=st.text(alphabet=st.sampled_from(AA + "KRKRP"), min_size=0, max_size=120),
    lo=st.integers(1, 12),
    width=st.integers(0, 30),
    missed=st.integers(0, 2),
)
def test_digest_equals_oracle(protein: str, lo: int, width: int, missed: int) -> None:
    from giapsim.proteomics import digest

    got = sorted(digest(protein, digest_len=(lo, lo + width), missed_cleavages=missed))
    want = sorted(protein[a:b] for a, b in trypsin_p(protein, (lo, lo + width), missed))
    assert sorted(set(got)) == sorted(set(want))
    for p in got:
        assert p in protein and lo <= len(p) <= lo + width


@settings(max_examples=100, deadline=None)
@given(protein=st.text(alphabet=st.sampled_from(AA), min_size=1, max_size=200))
def test_digest_default_rule(protein: str) -> None:
    from giapsim.proteomics import digest

    for p in digest(protein, digest_len=(7, 30), missed_cleavages=1):
        assert 7 <= len(p) <= 30
        assert sum(c in "KR" for c in p[:-1]) <= 1
        assert p[-1] in "KR" or protein.endswith(p)


@st.composite
def sim_cases(draw: Any) -> dict[str, Any]:
    from giapsim.proteomics import ProteomicsConfig

    prot = ProteomicsConfig(
        pi_wrong=draw(st.sampled_from([0.0, 0.01, 0.05, 0.2])),
        n_runs=draw(st.integers(1, 3)),
        detect_logit=draw(st.sampled_from([(-2.0, 0.08, 0.9), (50.0, 0.0, 0.0)])),
    )
    kw: dict[str, Any] = {
        "genome_length": draw(st.integers(40, 70)) * 1_000,
        "n_sorfs": draw(st.integers(0, 3)),
        "n_pseudogenes": draw(st.integers(0, 2)),
        "n_nterm_extensions": draw(st.integers(0, 2)),
        "origin_crossing_gene": draw(st.booleans()),
        "n_plasmids": draw(st.integers(0, 1)),
        "plasmid_length": (5_000, 8_000),
    }
    if draw(st.booleans()):
        kw["n_novel_orfs"] = (("intergenic", 1), ("antisense", 1))
    return {"seed": draw(st.integers(0, 2**31 - 1)), "prot": prot, "kw": kw}


@settings(max_examples=12, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(case=sim_cases())
def test_truth_agrees_with_genome(case: dict[str, Any]) -> None:
    with tempfile.TemporaryDirectory() as d:
        cfg = giapsim.SimConfig(proteomics=case["prot"], **case["kw"])
        giapsim.simulate(cfg, case["seed"], Path(d))
        run = Run(Path(d))
        pi = case["prot"].pi_wrong
        mix = tuple(case["prot"].score_mixture)
        hay = run.proteome_and_annotation_il
        for r in run.truth:
            if r["origin"] == "gene":
                p = run.proteins.get(r["locus_tag"])
                assert p is not None and r["locus_tag"] not in run.pseudo_loci, r
                assert r["peptide"] in p.aa, r
                assert maps_to(
                    r["peptide"], run.genome, r["contig"], int(r["start"]), int(r["end"]),
                    r["strand"],
                ), r  # fmt: skip
                assert int(r["n_loci_genome"]) >= 1
            elif r["origin"] in ("shuffled", "other_strain"):
                assert truthy(r["is_wrong"]) and il(r["peptide"]) not in hay, r
        if pi == 0.0:
            assert not [r for r in run.truth if truthy(r["is_wrong"])]
        for r in run.decoy_rows():
            assert all(empty(r.get(c)) for c in ("score", "pep", "q_value")), r
            locs = [i[4:] for i in r["protein_ids"].split(";")]
            assert all(i.startswith("rev_") for i in r["protein_ids"].split(";"))
            prots = [run.proteins[x].aa if x in run.proteins else run.json_aa[x] for x in locs]
            assert any(r["peptide"][::-1] in trypsin_p_peptides(a) for a in prots), r
        targets = sorted(run.target_rows(), key=lambda r: -float(r["score"]))
        qs = [float(r["q_value"]) for r in targets]
        assert all(a <= b + 1e-12 for a, b in itertools.pairwise(qs))
        for r in targets:
            s = float(r["score"])
            assert float(r["pep"]) == pytest.approx(mixture_pep(s, pi, mix), rel=1e-3, abs=1e-9)
            assert float(r["q_value"]) == pytest.approx(mixture_fdr(s, pi, mix), rel=1e-3, abs=1e-9)
        runs = {r["run"] for r in run.peptides}
        assert len(runs) == case["prot"].n_runs or not run.peptides
        df, issues = read_peptide_table(run.out / "peptides.tsv")
        _, issues2 = validate_peptides(df, InputConfig())
        assert not [i for i in [*issues, *issues2] if i.severity == "error"]
