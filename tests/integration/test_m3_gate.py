"""Gate M3 recall (GIAP_SPEC.md §13 Gate M3 "≥ 99.9 % recall of planted peptides"; issue M3-21
"recall gate | M3-19, pi_wrong = 0.01, 100 kb, 5 seeds | recall ≥ 0.999 (mapping_recall_min)").

Definition (docs/specs/M3-21 "Recall"): among correct, planted peptides from annotated CDS/sORFs
with exactly one genomic locus (``n_loci_genome = 1``) and input q ≤ ``input.q_value_max``, the
share with an accepted placement (a ``peptide_hits.tsv`` row) at the true locus (same contig,
1-based start/end and strand as the truth table) with ``locus_class`` in {``cds_inframe``,
``sorf``}. Alternative-start N-terminal peptides (``alt_start_nterm``; the six-frame translation
reads the start codon as V/L, not M) are excluded as in M3-17. Seeds ``calibration_base + k``,
k = 0..4; the recall is pooled over the five genomes and each genome must also reach the bar.
Failures list the first 20 missed peptides with their truth loci.

The golden half of the gate (the four mapping TSVs of the 50-kb fixture) is
tests/golden/test_mapping_golden.py (M3-17).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from helpers.mapping_steps import module2_files, read_tsv_df
from helpers.placement_sims import run_chain

N_SEEDS = 5


def _recall(tmp: Path, seed: int) -> tuple[int, list[tuple[str, ...]]]:
    from giapsim.proteomics import ProteomicsConfig

    prot = ProteomicsConfig(pi_wrong=0.01, search_db="other")
    ch = run_chain("database", seed, tmp, proteomics=prot)
    files = module2_files(ch.avail["mapping"].parent, ch.avail)
    hits = read_tsv_df(files["peptide_hits.tsv"])
    truth = read_tsv_df(ch.res.truth["peptides"])
    peps = pd.read_parquet(ch.avail["peptides"])
    qmax = ch.cfg.input.q_value_max
    confident = set(peps.loc[~peps["is_decoy"].astype(bool) & (peps["q_value"] <= qmax),
                             "peptide"])  # fmt: skip
    feats = read_tsv_df(ch.res.truth["features"])
    coding = set(feats.loc[feats["type"].isin(["cds", "sorf"]) & (feats["pseudo_reason"] == ""),
                           "feature_id"])  # fmt: skip
    t = truth[
        (truth["is_wrong"] == "false")
        & (truth["origin"] == "gene")
        & (truth["n_loci_genome"] == "1")
        & (truth["alt_start_nterm"] == "false")
        & truth["peptide"].isin(confident)
        & truth["locus_tag"].isin(coding)
    ]
    ok = hits["locus_class"].isin(["cds_inframe", "sorf"])
    found = {
        (r.peptide, r.contig, str(r.start), str(r.end), r.strand)
        for r in hits[ok].itertuples(index=False)
    }
    planted = [
        (r.peptide, r.contig, str(r.start), str(r.end), r.strand, r.locus_tag)
        for r in t.itertuples(index=False)
    ]
    missed = [p for p in planted if p[:5] not in found]
    return len(planted), missed


def test_m3_gate_recall_of_planted_peptides(tmp_path: Path, tolerances: dict[str, object]) -> None:
    base = int(tolerances["seeds"]["calibration_base"])  # type: ignore[index]
    need = float(tolerances["mapping_recall_min"])  # type: ignore[arg-type]
    total, missed_all, per_seed = 0, [], []
    for k in range(N_SEEDS):
        n, missed = _recall(tmp_path / f"s{k}", base + k)
        assert n > 200, f"seed {base + k}: only {n} planted peptides"
        total += n
        missed_all += missed
        per_seed.append((base + k, n, 1 - len(missed) / n))
    recall = 1 - len(missed_all) / total
    low = [s for s in per_seed if s[2] < need]
    assert recall >= need and not low, (
        f"recall {recall:.5f} (per seed {per_seed}) < {need}; first missed (peptide, contig,"
        f" start, end, strand, locus_tag): {missed_all[:20]}"
    )
