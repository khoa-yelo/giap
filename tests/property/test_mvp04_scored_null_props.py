"""Invariants of ``decoys = "scored_null"`` for any seed (MVP-04 part 2; spec step 6, rule 6).

On the test-time 200-kb excerpt (tests/helpers/mvp04_data.py), for random seeds:

* each null PSM (wrong peptide, run) is written exactly once, as its target row or as the decoy
  row of its entry; every decoy row is such a null PSM (no ``match_wrong`` extras);
* a decoy row is the reversed null peptide, ``rev_<entry>`` (an expressible recorded locus tag or
  a ``GIAPORF_…`` six-frame entry), with finite q / PEP / score and intensity NA;
* no decoy equals a target peptide (I = L); ``decoy_collisions`` is a non-negative int bounded by
  the null PSMs whose reversal collides;
* the q / PEP of every null row are the mixture's at its score and ``pi_wrong_realised``;
* lead decision 9: every protein-entry null is ``undetected_product``, a Trypsin/P product of the
  sample protein of its attributed locus, which has no true peptide, and its truth row maps it
  to that locus on the genome.

Interface: as tests/unit/test_mvp04_scored_null.py.
"""

from __future__ import annotations

import functools
import importlib
import math
import tempfile
from dataclasses import replace
from pathlib import Path

from helpers import mvp04_data as M
from helpers import mvp04_decoys as D
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giapsim.expression import IslandSpec

ISLANDS = ((2, 3.0), (18, -3.0), (34, 3.0), (50, -3.0), (66, 3.0), (82, -3.0))


@functools.cache
def _excerpt() -> M.Excerpt:
    return M.make_excerpt(Path(tempfile.mkdtemp(prefix="mvp04_p2_props_")))


@settings(
    max_examples=5, deadline=None, suppress_health_check=[HealthCheck.too_slow], derandomize=True
)
@given(seed=st.integers(min_value=0, max_value=2**31 - 1))
def test_scored_null_invariants_hold_for_any_seed(seed: int) -> None:
    ex = _excerpt()
    base = importlib.import_module("giapsim.scenarios").k12_evidence()
    expr = replace(
        base.expression, islands=tuple(IslandSpec(M.CONTIG, f, 8, e) for f, e in ISLANDS)
    )
    cfg = replace(base, expression=expr, n_novel_per_class=10,
                  proteomics=replace(base.proteomics, decoys="scored_null"))  # fmt: skip
    out = Path(tempfile.mkdtemp(prefix=f"mvp04_p2_props_{seed}_"))
    importlib.import_module("giapsim.real_genome").simulate_on_genome(
        cfg, seed, out, genome_fasta=ex.fasta, bakta_json=ex.bakta_json, tu_set=None,
        check_hashes=False,
    )  # fmt: skip
    v = D.view(out)
    assert v.provenance["decoys"] == "scored_null"
    bad = D.accounting_violations(v)
    assert not bad, bad[:10]

    recorded = D.annotation(ex.bakta_json)
    targets = {D.il(r["peptide"]) for r in v.target_rows} | {D.il(p) for p in v.correct}
    for n in v.nulls.values():
        for r in n.decoys.values():
            assert r["peptide"] == n.peptide[::-1], r
            assert D.il(r["peptide"]) not in targets, r
            assert r["intensity"] == "NA", r
            assert all(D.finite(r[c]) for c in ("q_value", "pep", "score")), r
            (pid,) = D.ids(r["protein_ids"])
            assert pid.startswith("rev_"), r
            entry = pid[len("rev_") :]
            if n.kind == "sixframe":
                assert D.ORF_ID.match(entry), r
            else:
                assert entry in recorded and D.expressible(recorded[entry]), r

    colliders = [n for n in v.nulls.values() if D.il(n.peptide[::-1]) in targets]
    k = v.provenance["decoy_collisions"]
    assert isinstance(k, int) and 0 <= k <= len(colliders) * len(v.runs)

    pi = float(v.meta["pi_wrong_realised"])
    mix = tuple(float(x) for x in cfg.proteomics.score_mixture)
    for _, _, r, _ in D.null_psms(v):
        s = float(r["score"])
        assert math.isclose(float(r["q_value"]), D.mixture_q(s, pi, mix), rel_tol=1e-6,
                            abs_tol=1e-12), r  # fmt: skip
        assert math.isclose(float(r["pep"]), D.mixture_pep(s, pi, mix), rel_tol=1e-6,
                            abs_tol=1e-12), r  # fmt: skip

    sample = D.read_seq(out / "truth" / "sample_genome.fna")
    proteins = D.locus_proteins(out, recorded, sample)
    nb = D.protein_null_violations(v, recorded, proteins, sample, ex.seq)
    assert not any(nb.values()), {k: x[:3] for k, x in nb.items() if x}
