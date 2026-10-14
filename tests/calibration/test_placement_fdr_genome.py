"""Realised placement FDR on giapsim genomes, exact route (issue M3-10; GIAP_SPEC.md §7.2.5, §8.4
"Placement FDR" row, fast tier: 100 kb).

Per replicate i (R = ``replicates_fast.placement``, seed ``seeds.calibration_base`` + i): a 100-kb
giapsim genome (M3-18 simulator); correct peptides are the tryptic peptides (7-20 aa) of its
annotated genes; wrong peptides are i.i.d. 7- and 8-mers drawn from the genome's six-frame residue
frequencies (``residue_freqs``), minus any that equal a correct peptide. All are database input
(reverse decoys, ``search_db = other``) and every PEP equals the known wrong share, so pi_wrong is
exact. The full chain runs: ``build_sixframe`` -> ``build_decoys`` -> ``match_exact`` ->
``best_costs`` -> ``placement``. Criterion (§8.4 "FDR"): mean FDP over R <= q + k SE(FDP) at every
level in ``levels``. The genome-level suite with tolerant hits and de novo errors is M3-21.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from giap.config import GiapConfig
from giap.io.fasta import read_genome_fasta
from giap.mapping.chance import residue_freqs
from giap.mapping.decoys import build_decoys
from giap.mapping.exact import match_exact
from giap.mapping.placement_fdr import best_costs, placement
from giap.mapping.sixframe import ALPHABET, build_sixframe
from giap.model import Replicon
from giap.model.genes import gene_model_from_bakta_json
from giap.rng import for_step
from giapsim.simulate import SimConfig, simulate

pytestmark = pytest.mark.calibration

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
    per, _ = placement(costs, GiapConfig())
    n_wrong_target = int(hits[(hits["search"] == "target")]["peptide"].isin(wrong).sum())
    return per, wrong, n_wrong_target


def test_realised_placement_fdr_on_giapsim(tolerances: dict[str, object], tmp_path: Path) -> None:
    levels = [float(x) for x in tolerances["levels"]]  # type: ignore[union-attr]
    reps, seeds = tolerances["replicates_fast"], tolerances["seeds"]
    assert isinstance(reps, dict) and isinstance(seeds, dict)
    r, base = int(reps["placement"]), int(seeds["calibration_base"])
    k = float(tolerances["se_multiplier"])  # type: ignore[arg-type]
    fdp = np.zeros((r, len(levels)))
    n_acc = np.zeros(r)
    n_false_hits = 0
    for i in range(r):
        per, wrong, n_wt = _replicate(base + i, tmp_path / f"rep{i}")
        n_false_hits += n_wt
        q = per["placement_q"].astype(float).to_numpy()
        is_wrong = per["peptide"].isin(wrong).to_numpy()
        for j, a in enumerate(levels):
            called = ~np.isnan(q) & (q <= a)
            n = int(called.sum())
            fdp[i, j] = (is_wrong & called).sum() / n if n else 0.0
        n_acc[i] = int((~np.isnan(q) & (q <= levels[-1])).sum())
    # the null is exercised: wrong peptides do hit the target by chance
    assert n_false_hits > r
    assert n_acc.mean() > 100
    mean = fdp.mean(axis=0)
    se = fdp.std(axis=0, ddof=1) / np.sqrt(r)
    for j, a in enumerate(levels):
        assert mean[j] <= a + k * se[j], (
            f"q={a}: mean FDP {mean[j]:.4f} > {a} + {k}*{se[j]:.4f} over R={r}"
        )
