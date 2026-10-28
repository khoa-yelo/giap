"""Realised placement FDR on gippasim genomes, exact route (issue M3-10; GIPPA_SPEC.md §7.2.5, §8.4
"Placement FDR" row, fast tier: 100 kb).

Per replicate i (R = ``replicates_fast.placement``; genome seed drawn from ``replicate_rng`` of
stream ``placement.genome_exact``, so ``_criteria.calibrate`` gates it with ADR-0015 two-block
confirmation and sharding): a 100-kb gippasim genome (M3-18 simulator); correct peptides are the
tryptic peptides (7-20 aa) of its
annotated genes; wrong peptides are i.i.d. 7- and 8-mers drawn from the genome's six-frame residue
frequencies (``residue_freqs``), minus any that equal a correct peptide. All are database input
(reverse decoys, ``search_db = other``) and every PEP equals the known wrong share. Since ADR-0019
the PEP is not read: π̂ is the per-length-class exact-share estimate, and the wrong 7-8-mers
concentrate the no-origin share in classes 7 and 8 (the length-correlated case of ADR-0019,
which a pooled π̂ cannot handle). The full chain
runs: ``build_sixframe`` -> ``build_decoys`` -> ``match_exact`` ->
``best_costs`` -> ``placement``. Criterion (§8.4 "FDR"): mean FDP over R <= q + k SE(FDP) at every
level in ``levels``. The genome-level suite with tolerant hits and de novo errors is M3-21.
"""

from __future__ import annotations

import re
from pathlib import Path

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels

from gippa.config import GippaConfig
from gippa.io.fasta import read_genome_fasta
from gippa.mapping.chance import residue_freqs
from gippa.mapping.decoys import build_decoys
from gippa.mapping.exact import match_exact
from gippa.mapping.placement_fdr import best_costs, placement
from gippa.mapping.sixframe import ALPHABET, build_sixframe
from gippa.model import Replicon
from gippa.model.genes import gene_model_from_bakta_json
from gippa.rng import for_step
from gippasim.simulate import SimConfig, simulate

pytestmark = pytest.mark.calibration

SUITE = "placement"

GENOME_BP = 100_000
N_WRONG = 20_000
RES = ALPHABET[:19]


def _tryptic(aa: str) -> list[str]:
    return [p for p in re.findall(r".*?[KR](?!P)|.+$", aa) if 7 <= len(p) <= 20]


def _replicate(seed: int, outdir: Path) -> tuple[pd.DataFrame, set[str], int]:
    res = simulate(SimConfig(genome_length=GENOME_BP), seed=seed, outdir=outdir)
    (rec,) = read_genome_fasta(res.genome_fasta)
    rep = Replicon(rec.id, len(rec.sequence), "circular")
    model = gene_model_from_bakta_json(res.bakta_json, {rec.id: rep})
    idx = build_sixframe({rec.id: rec.sequence}, {rec.id: rep})
    freqs = residue_freqs(idx)
    correct = {p for g in model.genes for p in _tryptic(g.aa) if set(p) <= set(RES)}
    rng = for_step("calibration_placement_genome", seed)
    lengths = rng.integers(7, 9, N_WRONG)
    codes = rng.choice(19, size=(N_WRONG, 8), p=freqs.f)
    wrong = {"".join(RES[c] for c in row[:L]) for row, L in zip(codes, lengths, strict=True)}
    wrong -= correct
    queries = sorted(correct | wrong)
    decoys = {"reverse": build_decoys(idx, "reverse", 50, for_step("sixframe_index", seed))}
    hits = match_exact(queries, idx, decoys, dict.fromkeys(queries, "reverse"))
    pi = len(wrong) / len(queries)
    peptides = pd.DataFrame(
        {
            "peptide": queries,
            "q_value": 0.001,
            "pep": pi,
            "is_decoy": False,
            "source": "database",
            "search_db": "other",
            "pi_denovo": np.nan,
        }
    )
    costs = best_costs(hits, peptides, freqs)
    per, _ = placement(costs, GippaConfig())
    n_wrong_target = int(hits[(hits["search"] == "target")]["peptide"].isin(wrong).sum())
    return per, wrong, n_wrong_target


def test_realised_placement_fdr_on_gippasim(tmp_path: Path) -> None:
    levels = C.load_tolerances().levels
    side: list[tuple[int, int]] = []  # (wrong target hits, accepted at loosest q), replicates run

    def fdp(rng: np.random.Generator) -> np.ndarray:
        seed = int(rng.integers(0, 2**31 - 1))
        per, wrong, n_wt = _replicate(seed, tmp_path / f"rep{seed}")
        q = per["placement_q"].astype(float).to_numpy()
        is_wrong = per["peptide"].isin(wrong).to_numpy()
        out = np.zeros(len(levels))
        for j, a in enumerate(levels):
            called = ~np.isnan(q) & (q <= a)
            n = int(called.sum())
            out[j] = (is_wrong & called).sum() / n if n else 0.0
        side.append((n_wt, int((~np.isnan(q) & (q <= levels[-1])).sum())))
        return out

    r = C.replicates(SUITE, full=False)
    stream = f"{SUITE}.genome_exact"
    calibrate_levels(SUITE, "fdr", levels, fdp, r, stream=stream).assert_ok()
    if side:
        n_wt, n_acc = np.mean(side, axis=0)
        assert n_wt > 1  # the null is exercised: wrong peptides do hit the target by chance
        assert n_acc > 100
