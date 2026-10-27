"""Properties of M7-12 round 3 (lead decisions 21 and 22; PR #225) over random seeds at 100 kb
(GIPPA_SPEC.md §7.6.6, §7.9.2, §11.3).

These hold for any seed of ``novel_orf_power``:

* every six-frame decoy row carries ``rev_`` + the ``db_id`` of the makedb entry
  (``gippa makedb --six-frame-orfs``, ``min_orf_aa = 7``) whose protein holds the reversed
  peptide, with the entry's class in the ID;
* every six-frame wrong target lies in exactly one makedb entry;
* decoy score, q and PEP are finite, and q and PEP sit in [0, 1];
* every other decoy is ``rev_<locus_tag>`` of an annotated protein;
* the four decision 22 territories each hold at least n // 4 - 2 of the n wrong peptides.

The interface is the one in tests/unit/test_sim_m7_round3.py, and the oracle is
tests/helpers/m7_round3_oracle.py. Each draw is a 100-kb simulation plus a makedb build, so few
examples are run.
"""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

from helpers import m7_round2_oracle as R
from helpers import m7_round3_oracle as Q
from helpers import m7_sim_oracle as O
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

import gippasim

seeds = st.integers(min_value=0, max_value=2**31 - 1)
_SETTINGS = settings(max_examples=3, deadline=None,
                     suppress_health_check=[HealthCheck.too_slow])  # fmt: skip


@_SETTINGS
@given(seed=seeds)
def test_sixframe_decoys_and_territories_for_any_seed(seed: int) -> None:
    import gippasim.scenarios as sc

    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "sim"
        gippasim.simulate(sc.novel_orf_power(seed), seed, out)
        entries = Q.makedb_entries(out, Path(td) / "db")
        index = Q.EntryIndex(entries, Q.lengths(out))
        t, d6, dp = Q.split_rows(out)
        assert d6 and t
        for r in d6:
            for pid in Q.ids(r["protein_ids"]):
                m = Q.DECOY_ID.match(pid)
                assert m is not None, r
                e = entries.get(pid[len("rev_") :])
                assert e is not None and e.cls == m.group(1), r
                assert r["peptide"][::-1] in e.aa, (r, e)
            vals = [float(r[c]) for c in ("score", "q_value", "pep")]
            assert all(math.isfinite(v) for v in vals), r
            assert 0.0 <= vals[1] <= 1.0 and 0.0 <= vals[2] <= 1.0, r
        prot = {k: f["aa"] for k, f in O.bakta_features(out / "bakta_like.json").items()}
        for r in dp:
            for pid in Q.ids(r["protein_ids"]):
                assert pid.startswith("rev_") and pid[4:] in prot, r
        loci = Q.wrong_loci(out)
        for s in loci.values():
            assert len(index.entry_of(s.contig, s.s1, s.e1, s.strand)) == 1, s
        genome = O.receiving_genome(out)
        cds, genes = R.cds_spans(out), Q.true_starts(out)
        by = {c: 0 for c in Q.CLASSES4}
        for s in loci.values():
            by[Q.wrong_class4(s, cds, genes, genome)] += 1
        n = sum(by.values())
        for c in Q.CLASSES4:
            assert by[c] >= n // 4 - 2, by
