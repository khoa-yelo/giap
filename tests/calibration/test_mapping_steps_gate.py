"""Mapping chain through the registered steps on giapsim: recall and placement FDR (issue M3-17;
GIAP_SPEC.md §7.2.9 "Gate M3", §8.4 "Placement FDR" row; definitions of docs/specs/M3-21).

M3-21 owns the full genome-level suites; these two tests check that the *wiring* (query
selection, routes, decoy types, π̂, placement, loci; helpers.mapping_steps.drive with
``giap.steps.registry.run_step``) keeps the guarantees the pure functions have:

* **Recall** (claim §13 Gate M3, threshold ``mapping_recall_min``): 100-kb genome, M3-19 database
  errors (``pi_wrong = 0.01``, ``search_db = other``), seed ``seeds.calibration_base``. Among
  correct peptides with ``n_loci_genome = 1``, not an alternative-start N-terminus, from an
  annotated CDS/sORF and with q ≤ ``input.q_value_max``, the share with a ``peptide_hits.tsv`` row
  (accepted, 1-based TSV coordinates as the truth table) at the true locus with ``locus_class``
  in {``cds_inframe``, ``sorf``}. Failures print the first 20 missed peptides.
* **Placement FDR** (claim §8.4: mean FDP over R ≤ q + 2 SE at every level of ``levels``): R =
  ``replicates_fast.placement`` 50-kb genomes (seed from ``replicate_rng`` of stream
  ``placement.steps_db``; ADR-0015 confirmation via ``calibrate``), M3-19 shuffled wrong peptides
  with ``pi_wrong = 0.05``, ``search_db = other``. FDP at level a = share of truly wrong peptides
  (truth ``is_wrong``) among the peptides of the ``placement`` work table with
  ``placement_q ≤ a``. 50 kb (not 100 kb) keeps the per-PR budget; M3-21 runs 100 kb / 1 Mb.
"""

from __future__ import annotations

from pathlib import Path

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels
from helpers.mapping_steps import (
    drive,
    fresh_steps,
    key_like,
    module2_files,
    read_tsv_df,
    roots_of,
    simulate_with_proteomics,
)

from giap.config import resolve_config

pytestmark = pytest.mark.calibration

SUITE = "placement"


def _chain(tmp: Path, seed: int, length: int, pi_wrong: float) -> tuple[dict, dict, object]:
    from giapsim.proteomics import ProteomicsConfig

    run = simulate_with_proteomics(
        tmp / "in",
        seed,
        length=length,
        proteomics=ProteomicsConfig(pi_wrong=pi_wrong, search_db="other"),
    )
    cfg = resolve_config(run.config_yaml, [f"seed={seed}"])
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        avail = drive(reg, cfg, tmp / "work", roots_of(run))
        keys = {"placement": key_like(reg.get("placement").outputs, "placement")}
    return avail, keys, run


def test_recall_of_planted_peptides_through_the_steps(
    tmp_path: Path, tolerances: dict[str, object]
) -> None:
    seed = int(tolerances["seeds"]["calibration_base"])  # type: ignore[index]
    avail, _, run = _chain(tmp_path, seed, 100_000, 0.01)
    files = module2_files(tmp_path / "work", avail)
    hits = read_tsv_df(files["peptide_hits.tsv"])
    truth = read_tsv_df(run.res.truth["peptides"])  # type: ignore[attr-defined]
    peps = pd.read_parquet(run.peptides)
    qmax = resolve_config(run.config_yaml).input.q_value_max
    confident = set(peps.loc[~peps["is_decoy"] & (peps["q_value"] <= qmax), "peptide"])
    feats = read_tsv_df(run.res.genes_truth)  # type: ignore[attr-defined]
    coding = set(feats["locus_tag"]) if "locus_tag" in feats.columns else None
    t = truth[
        (truth["is_wrong"] == "false")
        & (truth["n_loci_genome"] == "1")
        & (truth["alt_start_nterm"] == "false")
        & truth["peptide"].isin(confident)
    ]
    if coding is not None:
        t = t[t["locus_tag"].isin(coding)]
    ok_cls = hits["locus_class"].isin(["cds_inframe", "sorf"])
    found = {
        (r.peptide, r.contig, r.start, r.end, r.strand)
        for r in hits[ok_cls].itertuples(index=False)
    }
    planted = [(r.peptide, r.contig, r.start, r.end, r.strand) for r in t.itertuples(index=False)]
    assert len(planted) > 200
    missed = [p for p in planted if p not in found]
    recall = 1 - len(missed) / len(planted)
    need = float(tolerances["mapping_recall_min"])  # type: ignore[arg-type]
    assert recall >= need, f"recall {recall:.5f} < {need}; first missed: {missed[:20]}"


def test_placement_fdr_through_the_steps(tmp_path: Path) -> None:
    levels = C.load_tolerances().levels
    accepted: list[int] = []

    def fdp(rng: np.random.Generator) -> np.ndarray:
        seed = int(rng.integers(0, 2**31 - 1))
        avail, keys, run = _chain(tmp_path / f"rep{seed}", seed, 50_000, 0.05)
        per = pd.read_parquet(avail[keys["placement"]]).drop_duplicates("peptide")
        truth = read_tsv_df(run.res.truth["peptides"])  # type: ignore[attr-defined]
        is_w = (truth["is_wrong"] == "true") & (truth["origin"] != "decoy")
        wrong = set(truth.loc[is_w, "peptide"])
        q = per["placement_q"].astype(float).to_numpy()
        is_wrong = per["peptide"].isin(wrong).to_numpy()
        out = np.zeros(len(levels))
        for j, a in enumerate(levels):
            called = ~np.isnan(q) & (q <= a)
            n = int(called.sum())
            out[j] = (is_wrong & called).sum() / n if n else 0.0
        accepted.append(int((~np.isnan(q) & (q <= levels[-1])).sum()))
        return out

    r = C.replicates(SUITE, full=False)
    calibrate_levels(SUITE, "fdr", levels, fdp, r, stream=f"{SUITE}.steps_db").assert_ok()
    if accepted:
        assert np.mean(accepted) > 100  # the criterion is not vacuous
