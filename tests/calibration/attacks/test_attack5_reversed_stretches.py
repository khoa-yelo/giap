"""Attack 5 (GIAP_SPEC.md §8.5): de novo errors that reverse a stretch of sequence -> placement FDR
with shuffled decoys (a reversed decoy would inflate D). Issue M3-21; docs/methods/placement-fdr.md
"Why shuffle for de novo input"; M3-21 lead decisions 1-5. "A failing attack blocks the merge."

Generator: 100-kb giapsim genome (1 Mb full), M3-19 table (pi_wrong 0.05, ``search_db = other``),
M3-20 de novo reads of the same spectra with ``error_weights = {reversed_stretch: 1.0}`` (every
erroneous read carries one reversed 3-6-residue stretch), reference peptides. De novo input ->
shuffle decoys (``mapping.decoys.denovo_input`` default). The false-placement definition, the
per-list and pooled gates and the unique guard are those of ``_placement_suite.run_fdr_suite``;
a de novo read's true origin is its spectrum's true peptide whatever the error (decision 2).
Exercised (measured at 100 kb, 2026-10-20): shuffle m = 1 (about 30-40 per replicate) and m = 2
(about 500) -- reversed stretches are mostly cost 2 at the origin.

Documentation, not a gate: the first ``DOC_N`` replicates are mapped again with reversed decoys
for the de novo input; decoy-hit counts D per stratum under both decoys (and the FDP under the
reversed decoy) go to ``placement_m3_21_diagnostics.jsonl``.
"""

from __future__ import annotations

from pathlib import Path

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _placement_suite import run_fdr_suite, write_diagnostics
from helpers.placement_sims import run_chain, score_chain

pytestmark = pytest.mark.calibration

LEVELS = C.load_tolerances().levels
DOC_N = 10
BASE = "placement.attack5"


def _kw(full: bool) -> dict[str, object]:
    from giapsim.denovo import DenovoConfig
    from giapsim.proteomics import ProteomicsConfig

    return {
        "length": 1_000_000 if full else 100_000,
        "proteomics": ProteomicsConfig(pi_wrong=0.05, search_db="other"),
        "denovo": DenovoConfig(error_weights={"reversed_stretch": 1.0}),
    }


def _decoys_by_stratum(avail: dict[str, Path]) -> list[int]:
    f = pd.read_parquet(avail["mapping_fdr"])
    f = f[(f["section"] == "placement") & (f["length_class"] == "all")]
    f = f[f["decoy_type"] == "shuffle"].set_index("n_diff")  # de novo rows (label by source)
    return [int(f.loc[m, "decoys"]) for m in (0, 1, 2)]


def _run(full: bool, tmp: Path) -> None:
    kw = _kw(full)
    d_shuffle: dict[int, list[int]] = {}

    def replicate(seed: int, d: Path):  # type: ignore[no-untyped-def]
        ch = run_chain("denovo", seed, d, **kw)  # type: ignore[arg-type]
        d_shuffle[seed] = _decoys_by_stratum(ch.avail)
        return score_chain(ch, LEVELS, d)

    memo = run_fdr_suite(BASE, replicate, tmp, [("shuffle", 1), ("shuffle", 2)], full=full)
    doc = []
    for i in sorted(memo)[:DOC_N]:
        seed = int(C.replicate_rng(BASE, i).integers(0, 2**31 - 1))
        ch = run_chain(
            "denovo", seed, tmp / f"rev{i}", overrides=["mapping.decoys.denovo_input=reverse"],
            **kw,  # type: ignore[arg-type]
        )  # fmt: skip
        d_rev = _decoys_by_stratum(ch.avail)
        rev = score_chain(ch, LEVELS, tmp / f"rev{i}")
        doc.append((d_shuffle[seed], d_rev, rev.fdp["pooled"]))
    if doc:
        write_diagnostics(
            {
                "suite": f"{BASE}.reverse_vs_shuffle_doc",
                "n": len(doc),
                "mean_D_shuffle_by_m": np.mean([d[0] for d in doc], axis=0).tolist(),
                "mean_D_reverse_by_m": np.mean([d[1] for d in doc], axis=0).tolist(),
                "pooled_fdp_reverse": np.mean([d[2] for d in doc], axis=0).tolist(),
                "pooled_fdp_shuffle": np.mean(
                    [memo[i].fdp["pooled"] for i in sorted(memo)[:DOC_N]], axis=0
                ).tolist(),
            }
        )


def test_attack5_reversed_stretches_shuffle_decoys(tmp_path: Path) -> None:
    _run(False, tmp_path)


@pytest.mark.calibration_full
def test_attack5_reversed_stretches_shuffle_decoys_full(tmp_path: Path) -> None:
    _run(True, tmp_path)
