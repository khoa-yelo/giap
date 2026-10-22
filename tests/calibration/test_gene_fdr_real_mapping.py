"""One gene-FDR sanity replicate through the real M3 mapping (docs/specs/M4-21-m2-m3-conformance.md
work item 5; GIAP_SPEC.md §7.3.2, §8.4 suite "Gene FDR (picked, MAYU)").

Sanity, not gated: one replicate, nothing appended to the suite log (``_criteria.log_result``), no
two-block confirmation. Suite ``gene_fdr`` (SUITE_REGISTRY file entry).

Generator: giapsim with a 500-kb genome (453 genes) and M3-19 proteomics at
``pi_wrong`` = 0.05 (the suite's larger pi) and a lower detection intercept (-5 instead of -2) so
that about a quarter of the genes have no correct peptide; seed ``suite_seed(STREAM, 0)``. Then
``giap run --executor local`` maps the peptides for real (``sixframe_index`` ... ``loci``) and
``run_gene_table`` reruns detection on that work directory at ``stats.alpha.identification`` = q
for each q in ``tolerances.yaml`` ``levels``.

Truth: a gene is present when the sample holds at least one correct peptide of it
(``truth/peptides.tsv``: ``origin = gene``, ``is_wrong = false``); every giapsim protein is
expressed, so "absent" means "not identified by any correct peptide". FDP = detected absent genes /
detected genes.

Bound (one replicate): FDP <= q + ``se_multiplier`` * sqrt(q (1 - q) / n_detected), the binomial
SE of one replicate's FDP at the nominal rate -- the fast-suite FDR criterion with R = 1.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import _criteria as C
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m4_chain as mc
import m4_real as mr

STREAM = "gene_fdr.real_mapping"
GENOME_BP = 500_000


@pytest.mark.calibration
def test_gene_fdr_one_replicate_through_real_mapping(tmp_path: Path) -> None:
    from giap.steps.gene_table import run_gene_table
    from giapsim.proteomics import ProteomicsConfig
    from giapsim.simulate import SimConfig, simulate

    tol = C.load_tolerances()
    seed = C.suite_seed(STREAM, 0) % 2**31
    proteomics = ProteomicsConfig(pi_wrong=0.05, detect_logit=(-5.0, 0.08, 0.9))
    sim = simulate(SimConfig(genome_length=GENOME_BP, proteomics=proteomics,
                             write_precomputed=True), seed, tmp_path / "sim")  # fmt: skip
    db = mc.test_db(tmp_path / "db")
    work = mr.giap_run(sim, tmp_path / "OUT", db, seed)

    truth = pd.read_csv(sim.truth["peptides"], sep="\t", dtype=str, keep_default_na=False)
    correct = truth[(truth["origin"] == "gene") & (truth["is_wrong"].str.lower() == "false")]
    present = set(correct["locus_tag"])
    assert (truth["is_wrong"].str.lower() == "true").sum() > 0, "no wrong peptide planted"

    failures = []
    for q in tol.levels:
        cfg = mc.config(seed, stats={"alpha": {"identification": q}})
        with mc.db_env(db):
            pre = pd.read_parquet(run_gene_table(work, cfg))
        absent = ~pre["locus_tag"].astype(str).isin(present)
        assert absent.sum() >= 50, f"only {int(absent.sum())} absent genes: FDP would be vacuous"
        det = pre["detected"].astype(bool)
        assert set(pre.loc[det, "detection_tier"]) <= {"decoy_picked"}, "decoys were emitted"
        n_det = int(det.sum())
        assert n_det > 0, f"nothing detected at q = {q}"
        fdp = int((det & absent).sum()) / n_det
        bound = q + tol.se_multiplier * math.sqrt(q * (1 - q) / n_det)
        if fdp > bound:
            failures.append((q, n_det, round(fdp, 4), round(bound, 4)))
    assert not failures, f"FDP above the one-replicate bound (q, detected, FDP, bound): {failures}"
