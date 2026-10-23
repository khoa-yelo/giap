"""Properties of the M7-12 giapsim scenarios over random seeds (GIAP_SPEC.md §11.3; issue M7-12).

For any seed: the variant truth agrees with ``strain_differences.tsv`` and its rule-2 counts with
a brute-force recount; every planted novel ORF reads start..stop on the receiving genome with its
``stop``/``true_start`` columns consistent; isobaric de novo errors stay within the four isobaric
kinds. Conventions as in tests/unit/test_sim_m7_scenarios.py; the oracle is
tests/helpers/m7_sim_oracle.py. Few examples (each draw is a 100-kb simulation).
"""

from __future__ import annotations

import tempfile
from pathlib import Path

from helpers import m7_sim_oracle as O
from hypothesis import HealthCheck, example, given, settings
from hypothesis import strategies as st

import giapsim

seeds = st.integers(min_value=0, max_value=2**31 - 1)
_SETTINGS = settings(max_examples=3, deadline=None,
                     suppress_health_check=[HealthCheck.too_slow])  # fmt: skip


def _sim(name: str, seed: int, d: Path) -> Path:
    import giapsim.scenarios as sc

    giapsim.simulate(getattr(sc, name)(seed), seed, d)
    return d


@_SETTINGS
@given(seed=seeds)
@example(seed=29297)  # non-ATG initiator SNV (review #225 R2)
def test_alt_sites_consistent_for_any_seed(seed: int) -> None:
    with tempfile.TemporaryDirectory() as td:
        out = _sim("variant_f1", seed, Path(td))
        sample = O.read_fasta(out / "genome.fasta")
        ref_genes = O.bakta_features(out / "reference_bakta_like.json")
        ref_seq = O.read_fasta(out / "reference.fasta")
        snv = {}
        for s in O.nonsyn_snvs(out):
            cs = O.ref_codon_start(ref_genes[s["gene"]], int(s["pos_reference"]),
                                   len(ref_seq[s["contig"]]))  # fmt: skip
            snv[(s["contig"], cs, s["gene"])] = s
        correct = O.correct_peptide_rows(out)
        for r in O.rows(out / "truth" / "alt_sites.tsv"):
            if O.truthy(r["frameshift_upstream"]):
                continue  # sample-frame lift, checked by the round-2 tests (decision 11)
            s = snv[(r["contig"], int(r["codon_start"]), r["gene"])]
            assert (r["ref_aa"], r["alt_aa"]) == (s["ref_aa"], s["alt_aa"])
            n = O.recount_site_peptides(out, s["gene"], s["contig"], int(s["pos_sample"]),
                                        len(sample[s["contig"]]), correct)  # fmt: skip
            assert int(r["n_true_peptides"]) == n


@_SETTINGS
@given(seed=seeds)
def test_novel_truth_consistent_for_any_seed(seed: int) -> None:
    with tempfile.TemporaryDirectory() as td:
        out = _sim("novel_orf_power", seed, Path(td))
        genome = O.receiving_genome(out)
        correct = O.correct_peptide_rows(out)
        for r in O.rows(out / "truth" / "novel_orfs.tsv"):
            s1, e1, st_ = int(r["start"]), int(r["end"]), r["strand"]
            nt = O.strand_slice(genome[r["contig"]], s1 - 1, e1, st_)
            assert nt[:3] in O.STARTS and nt[-3:] in O.STOPS
            assert O.translate_orf(nt) == r["aa"]
            assert int(r["stop"]) == ((e1 - 2) if st_ == "+" else s1)  # 1-based (decision 16)
            assert int(r["true_start"]) == (s1 if st_ == "+" else e1)
            assert int(r["n_true_peptides"]) == O.recount_orf_peptides(out, r, correct)


@_SETTINGS
@given(seed=seeds)
def test_isobaric_errors_only_for_any_seed(seed: int) -> None:
    with tempfile.TemporaryDirectory() as td:
        out = _sim("isobaric_attack", seed, Path(td))
        kinds = {t["error_type"] for t in O.rows(out / "truth" / "denovo.tsv")} - {"none"}
        assert kinds <= O.ISOBARIC_TYPES
        assert O.rows(out / "truth" / "alt_sites.tsv") == []
