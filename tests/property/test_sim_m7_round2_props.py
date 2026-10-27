"""Properties of M7-12 round 2 (lead decisions 11, 12, 15 and 16; issue M7-12, PR #225) over
random seeds at 100 kb (GIPPA_SPEC.md §11.3).

For any seed: ``variant_f1``'s ``alt_sites.tsv`` equals the sample-frame lift of every
non-synonymous SNV (ambiguous codons near indels / at first codons excluded) and its
``frameshift_upstream`` agrees with the gene's upstream indels; ``novel_orf_power`` simulates
without raising (decision 15), writes ``stop`` 1-based (decision 16) and its six-frame wrong
peptides (decision 12) translate from their truth locus, out of frame of every coding span and
absent from every true protein. Interface as in tests/unit/test_sim_m7_round2.py; the oracle is
tests/helpers/m7_round2_oracle.py. Few examples (each draw is a 100-kb simulation).
"""

from __future__ import annotations

import tempfile
from pathlib import Path

from helpers import m7_round2_oracle as R
from helpers import m7_sim_oracle as O
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

import gippasim

seeds = st.integers(min_value=0, max_value=2**31 - 1)
_SETTINGS = settings(max_examples=3, deadline=None,
                     suppress_health_check=[HealthCheck.too_slow])  # fmt: skip


def _sim(name: str, seed: int, d: Path) -> Path:
    import gippasim.scenarios as sc

    gippasim.simulate(getattr(sc, name)(seed), seed, d)
    return d


@_SETTINGS
@given(seed=seeds)
def test_alt_sites_sample_frame_for_any_seed(seed: int) -> None:
    with tempfile.TemporaryDirectory() as td:
        out = _sim("variant_f1", seed, Path(td))
        ft = R.sample_frame_alt_sites(out)
        got = {}
        for r in O.rows(out / "truth" / "alt_sites.tsv"):
            site = (r["contig"], int(r["codon_start"]), r["gene"])
            if site not in ft.ambiguous:
                got[(r["contig"], int(r["codon_start"]), r["strand"], r["gene"], r["ref_aa"],
                     r["alt_aa"])] = O.truthy(r["frameshift_upstream"])  # fmt: skip
        want = {k for k in ft.want if (k[0], k[1], k[3]) not in ft.ambiguous}
        assert set(got) == want
        for k, flag in got.items():
            any_fs, net = ft.shift[(k[0], k[1], k[3])]
            assert flag or not net, k
            assert any_fs or not flag, k


@_SETTINGS
@given(seed=seeds)
def test_novel_orf_power_never_raises_and_truth_holds(seed: int) -> None:
    with tempfile.TemporaryDirectory() as td:
        out = _sim("novel_orf_power", seed, Path(td))
        genome = O.receiving_genome(out)
        for r in O.rows(out / "truth" / "novel_orfs.tsv"):
            s1, e1, st_, stop = int(r["start"]), int(r["end"]), r["strand"], int(r["stop"])
            assert stop == (e1 - 2 if st_ == "+" else s1)
            assert O.strand_slice(genome[r["contig"]], stop - 1, stop + 2, st_) in O.STOPS
        coding = R.coding_spans(out)
        prots = [f["aa"] for f in O.bakta_features(out / "bakta_like.json").values()
                 if f.get("aa")]  # fmt: skip
        hay = "|".join(O.il(p) for p in prots + [r["aa"] for r in
                       O.rows(out / "truth" / "novel_orfs.tsv")])  # fmt: skip
        wrong = [r for r in O.rows(out / "truth" / "peptides.tsv")
                 if O.truthy(r["is_wrong"]) and r["origin"] != "decoy"]  # fmt: skip
        assert wrong
        for r in wrong:
            assert r["origin"] == "sixframe_noncoding"
            s1, e1 = int(r["start"]), int(r["end"])
            nt = O.strand_slice(genome[r["contig"]], s1 - 1, e1, r["strand"])
            assert R.translate_nt(nt) == r["peptide"]
            assert O.il(r["peptide"]) not in hay
            span = R.Span(r["contig"], s1, e1, r["strand"])
            assert not any(R.in_frame_with(span, g) for g in coding)
