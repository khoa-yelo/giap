# Detection tiers and gene-level FDR (M4-12)

Status: statistician re-review 2026-10-15 (PR #92, after f297600): PEP-only path approved;
**changes requested** on the contaminant-decoy rule (lead decision 11), see Sign-off. Code: `src/giap/genetable/detection.py`.
Spec: GIAP_SPEC.md §2.1 decision 4, §7.3, §7.3.2, §9.3; ADR-0009; docs/specs/M4-12-detection-tiers.md.

## Tier choice (once per run, from the input)
1. `decoy_picked`: at least one C9 decoy row with `protein_ids` passes q <= `input.q_value_max`
   (NA q passes), and the paired share f = paired / (paired + unpaired) is at least
   max(`mapping.proteins.min_identity`, `mapping.proteins.min_coverage`) (0.80 by default).
   Decoys above `q_value_max` are neither paired nor unpaired.
2. `mayu`: decoys exist but f is below the cut-off (warning `decoy_pairing_low`), or an accepted
   target carries a q-value or a PEP.
3. `two_peptide`: C11 has rows but no accepted q or PEP.
4. `protein_list`: no peptides at all.

## Protein to gene map
Priority: (1) the ID is a C2 `locus_tag`; (2) the `gene` of its C12 rank-1 hits (one gene only);
(3) peptide vote: among accepted unique target peptides whose C9 `protein_ids` contain the ID, the
gene with the most distinct peptides (a tie leaves the ID unmapped). Protein groups `A;B` are split;
contaminant-prefixed IDs are never mapped.

## Decoy pairing
The longest matching `input.decoy_prefixes` entry is stripped from every member of the decoy's
`protein_ids`. A decoy with any stripped member carrying an `input.contaminant_prefixes` prefix is
a contaminant decoy and is neither paired nor unpaired (lead decisions 6, 11). Otherwise the decoy
pairs when every member maps, and all map to the same gene.

## PEP-only input (lead decision 5)
When every accepted q-value is missing and PEPs are given, peptide q is the running mean of the
sorted distinct accepted PEPs (q_i = mean of the PEPs <= PEP_i), and MAYU uses those q for the
thresholds and sum PEP for F(t). Mixed input keeps the given q (decision 11).

Asymmetry (lead question 1). Decoys pass the same q <= `q_value_max` filter as targets, but targets
must also be placed on a gene and decoys cannot be. False targets that do not place are lost while
their decoys stay, so D is too large and the picked estimate errs conservative. In giapsim every
false target places (`false_source = random_gene`), which is the least favourable case for this
asymmetry.

What 1/f assumes (ADR-0009). Decoys are lost **at random with respect to score** while their target
counterparts stay. When loss is symmetric (a reference protein absent from the strain loses its
false targets and its decoy alike, as in cross-strain searches), 1/f over-corrects and the estimate
is conservative. Loss that depends on score is not corrected. Stress run below: losing the
*best-scoring* 10 % / 20 % of passing decoys gives FDP 1.5x / 2x q at q = 0.01, and the 80 % guard
does not catch it. Nothing in the input should make pairing depend on score. If an engine ever
does that, a KS check of paired against unpaired decoy scores would be the diagnostic (backlog).

## Picked gene-level FDR (Savitski 2015, ADR-0009)
Per eligible (non-partial) gene: best accepted target and best paired decoy by the key
(q ascending, PEP ascending, score descending; NA last), dense-ranked over the union (no -log10 0).
The better one is kept; a tie goes to the decoy. With f the paired share,
FDR(t) = (D(t)/f + 1) / T(t), made monotone. Decoy winners and genes without accepted
peptides get q NA (`not_detected`). `detection_summary.tsv` reports f (`decoy_paired_share`) and 1/f
(`decoy_win_weight`).

## MAYU tier
`giap.stats.mayu.mayu_gene_fdr` (see mayu-gene-fdr.md) on eligible genes, with sizes
`n_observable_peptides`, P(t) = distinct accepted peptides, F = sum of PEP when every accepted peptide has
one. Otherwise F = t |P(t)|, with warning `pep_partial` when only some have one. `fdr_2pep` is reported
under its own name and is not the gene FDR. Warning `unspecific_search`: more than 50 % of distinct
accepted peptides (placements ending at the gene's C-terminus excepted) do not end in K/R.

## Unestimated tiers
`two_peptide`: >= `detection.min_peptides_sequences_only` distinct accepted peptides -> `detected`;
exactly one -> `single_peptide_unestimated` (not detected in the primary analysis). `protein_list`:
a C12 rank-1 hit with a gene -> `detected`. Both carry no gene q.
Contaminant evidence counts toward detection in every tier (lead decision 2).

## Calibration evidence
### Implementer run (fast suite, R = 50, block 1; mean FDP +- SE; Slurm job 1849087)
| claim | pi_wrong | q = 0.01 | q = 0.05 | q = 0.10 |
|---|---|---|---|---|
| picked | 0.01 | 0.0094 +- 0.0005 | 0.0488 +- 0.0009 | 0.0522 +- 0.0010 |
| picked | 0.05 | 0.0081 +- 0.0004 | 0.0468 +- 0.0011 | 0.0515 +- 0.0009 |
| MAYU (PEP and no PEP) | 0.01 | 0.0104 +- 0.0004 | 0.0490 +- 0.0010 | 0.0541 +- 0.0010 |
| MAYU (PEP and no PEP) | 0.05 | 0.0091 +- 0.0004 | 0.0482 +- 0.0010 | 0.0538 +- 0.0009 |
| picked, 10 % decoy loss | 0.05 | 0.0094 +- 0.0005 | 0.0485 +- 0.0013 | 0.0526 +- 0.0010 |
| picked, 20 % decoy loss | 0.05 | 0.0091 +- 0.0006 | 0.0458 +- 0.0011 | 0.0510 +- 0.0010 |
| 30 % decoy loss (any tier) | 0.05 | 0.0100 +- 0.0004 | 0.0499 +- 0.0008 | 0.0552 +- 0.0008 |

All claims pass (mean <= q + k SE). Notes for review: the MAYU PEP and no-PEP variants are identical
because giapsim's q-values are running means of PEP, so sum PEP = t |P(t)| at every threshold. The FDP
saturates near 0.05 at q = 0.10 because accepted peptides already pass q <= 0.01. Recorded, not gated:
`two_peptide` FDP 0.003 (pi 0.01) and 0.053 (pi 0.05).

### Statistician replication (Slurm job 1849140; R = 1,000 at 2,000 genes, R = 300 at 4,500 genes)
The pipeline is the suite's own (`detection_tiers_calsim`: M4-07 -> M4-08 -> M4-09 -> M4-10/11 ->
`detect`). Seeds are `SeedSequence([20261015, 1000 pi, n_genes, rep])`, independent of the harness
streams. The script is local/scratch/review-m4-12/exp.py. The table gives mean FDP +- SE at
q = 0.01, and in brackets the ratio of mean FDP to the mean estimated FDR of the last gene detected.

| claim | 2,000 genes, pi 0.01 | 2,000, pi 0.05 | 4,500, pi 0.01 | 4,500, pi 0.05 |
|---|---|---|---|---|
| picked | 0.00979 +- 0.00013 | 0.00926 +- 0.00013 | 0.00970 +- 0.00016 | 0.00949 +- 0.00015 |
| MAYU (PEP = no PEP) | 0.00996 +- 0.00009 [1.08] | 0.00979 +- 0.00009 [1.04] | 0.01004 +- 0.00011 [1.04] | 0.00979 +- 0.00011 [1.01] |
| picked, q from TDA D/T | 0.01013 +- 0.00013 | 0.00960 +- 0.00013 | | |
| MAYU PEP, q from TDA | 0.00963 +- 0.00009 | 0.00947 +- 0.00009 | | |
| MAYU no PEP, q from TDA | 0.01010 +- 0.00011 | 0.00993 +- 0.00011 | | |
| picked, 10 % / 20 % random loss | 0.00980 / 0.00985 | 0.00932 / 0.00947 | | |
| same, without 1/f (pre-ADR mutant) | 0.01073 / 0.01216 | 0.01030 / 0.01168 | | |
| picked, best-scoring 10 % / 20 % lost | 0.0148 / 0.0197 | 0.0145 / 0.0193 | | |

At q = 0.05 every row except the mutant lies within 0.0475-0.0500. The pre-ADR mutant at 20 % loss reaches 0.0526 /
0.0522 (z ~ +11-15). At q = 0.10 every row equals its q = 0.05 FDP (saturation).

Reading:
1. **MAYU at q = 0.01 is at level, with no margin.** The implementer's 0.0104 +- 0.0004 was a +1 SE
   block. With 3x-6x the replicates the mean is 0.0100 (pi 0.01) and 0.0098 (pi 0.05). The estimate
   itself runs 4-8 % low (ratio > 1, as in mayu-gene-fdr.md: noisy a_g, bin averaging). The step
   between gene q values absorbs this, so the realised mean lands on q. This is noise around an
   unbiased-at-q procedure, not material anti-conservatism. Nightly false-alarm cost (normal
   approximation, two-block rule): MAYU pi 0.01 at q = 0.01 and R = 300 fails one block with
   probability ~0.05 and both with ~0.0024. Every other gene_fdr claim, the attack included, is
   below 1e-3. The suite adds ~0.003 to the ~0.04 nightly budget of ADR-0015.
2. **PEP against no PEP.** In giapsim q is the running mean of PEP, so sum PEP equals t |P(t)| and
   the two gated claims are one claim. The M4-09 simulator is the one already in use, so it does not
   separate them on its own. Once q comes from target-decoy competition (D/T over the emitted
   decoys), the two differ (0.0096 against 0.0101 at pi 0.01) and both stay at level. An
   `EvidenceConfig` mode with an engine-like q is needed before the two can be counted as
   independent evidence (backlog).
3. **Decoy-loss attack (ADR-0009).** With 1/f, random loss of 10 % or 20 % keeps the level. Without
   1/f, 20 % loss inflates FDP by 1.2-1.25x. The attack catches that mutant with probability ~0.4-0.5
   per level in the fast suite (both blocks must fail) and ~1.0 in the full suite. At 10 % loss its
   power is low: ~0.004 fast and ~0.04 full at q = 0.01. With floor(0.3 n) decoys lost, f = 0.70 <
   0.80 holds in every replicate, so the "30 %" claim always runs MAYU. It tests the guard, not 1/f.
4. **Levels 0.05 and 0.10 are weak and vacuous.** Peptides are capped at q <= 0.01, so the loosest
   gene FDR reachable is ~0.053 and the q = 0.10 claims sit ~200 SE below the bound. Giving them
   teeth needs a suite variant with `input.q_value_max` = 0.10 (backlog).
5. Null structure: false targets and decoys land on genes in proportion to `n_observable` (§8.4),
   which is MAYU's own model. Departures from it are covered by the M4-02 attacks (mayu-gene-fdr.md).
   TU blocks and autocorrelation do not enter gene detection. One test family per run (all eligible
   genes), as in §8.1. Effect size: the gene q is reported next to M4-11's `n_accepted` and
   `coverage_frac`. The FDR curve rows carry no CI (`ci_low`/`ci_high` NA). They are estimates of an
   expectation, as in mayu-gene-fdr.md.

### Statistician re-review after f297600 (Slurm job 1850323 + 1850540; R = 1,000 at 2,000 genes, R = 300 at 4,500)
Script local/scratch/review-92/exp.py; seeds `SeedSequence([20261015, 92, 1000 pi, n_genes, rep])`.
Mean FDP +- SE.

PEP-only input (decoys dropped, C9 `q_value` and C11 `best_q_input` NA, PEP kept; tier `mayu` in
every replicate). Every target peptide passes the q filter, so q = 0.05 / 0.10 are reachable:

| n genes, pi | q = 0.01 | q = 0.05 | q = 0.10 |
|---|---|---|---|
| 2,000, 0.01 | 0.00973 +- 0.00009 | 0.04886 +- 0.00018 | 0.05514 +- 0.00018 |
| 2,000, 0.05 | 0.00965 +- 0.00009 | 0.04926 +- 0.00021 | 0.09831 +- 0.00028 |
| 4,500, 0.01 | 0.00983 +- 0.00010 | 0.04940 +- 0.00022 | 0.05538 +- 0.00020 |
| 4,500, 0.05 | 0.00995 +- 0.00011 | 0.04945 +- 0.00025 | 0.09871 +- 0.00035 |

At q <= 0.05 it equals the `mayu_pep` claim (giapsim q is the running mean of PEP, so the derived
q reproduces it on the accepted prefix). At q = 0.10 and pi 0.05 it is the first gene-FDR claim of
M4-12 that is informative at 0.10 (FDP 0.098, at level). Approved.

Contaminant-shared groups (attack; picked input, a share X of the genes has `CON__<g>` /
`rev_CON__<g>` appended to target / decoy `protein_ids`; targets still place by sequence):

| n genes, pi | X = 0.2, q = 0.01 | X = 0.4, q = 0.01 | X = 0.4, q = 0.05 |
|---|---|---|---|
| 2,000, 0.01 | 0.01184 +- 0.00015 | 0.01597 +- 0.00019 | 0.05302 +- 0.00018 |
| 2,000, 0.05 | 0.01176 +- 0.00015 | 0.01575 +- 0.00019 | 0.05316 +- 0.00020 |
| 4,500, 0.01 | 0.01199 +- 0.00017 | 0.01621 +- 0.00021 | 0.05348 +- 0.00020 |
| 4,500, 0.05 | 0.01232 +- 0.00019 | 0.01657 +- 0.00024 | 0.05334 +- 0.00024 |

The tier is `decoy_picked` in every replicate: the mixed decoys leave both counts, so f stays ~1
while D loses a share X of the gene decoys and FDP ~ q / (1 - X). The rule rewards the input with
the higher tier exactly when its decoy side is depleted, and nothing bounds X (an organism whose
proteins overlap the contaminant FASTA, or a contaminant FASTA that carries organism proteins).
Gated as tests/calibration/attacks/test_detection_tiers_attack_inputs.py (X = 0.3): fast run
(job 1850411) is a confirmed fail at q = 0.01 (0.0143 / 0.0138 against bound 0.0115) and at
q = 0.05; the `pep_only` claim of the same file passes at all three levels.

## Sign-off
Re-review (2026-10-15, after f297600 and lead decision 11): **changes requested.**
1. Decision 5 (PEP-only q) is implemented and calibrated (table above); approved.
2. Decision 11 as implemented ("ANY stripped ID is a contaminant -> out of both counts") is
   anti-conservative for decoys of mixed groups `rev_P;rev_CON__X`: their targets `P;CON__X`
   still place on P's gene by sequence, so dropping only the decoy is unaccounted decoy loss.
   Required: a decoy is a contaminant decoy only when ALL its stripped members carry a contaminant
   prefix; for mixed decoys drop the contaminant members and pair on the rest (or, at least,
   count them as unpaired so 1/f and the 80 % guard see them). The attack file above must pass.
3. Decisions 7-9 verified (longest prefix, `Detection.protein_gene_map`, mutant tests).

Statistician (M4-12 review, 2026-10-15): the gene-level FDR claims (picked with ADR-0009's 1/f, and
MAYU) are **approved** as calibrated on the evidence above. Sign-off of the module is **conditional**
on two code fixes, listed in the PR #92 review: (a) PEP-only input reaches the `mayu` tier but
detects no gene, because the gene thresholds come from q alone; (b) decoys of contaminant proteins
count as unpaired and lower f.
