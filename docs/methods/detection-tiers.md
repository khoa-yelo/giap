# Detection tiers and gene-level FDR (M4-12)

Status: implementer draft, awaiting statistician sign-off. Code: `src/giap/genetable/detection.py`.
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
The first matching `input.decoy_prefixes` entry is stripped from every member of the decoy's
`protein_ids`. The decoy pairs when every member maps, and all map to the same gene.

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

## Calibration evidence (fast suite, R = 50, block 1; mean FDP +- SE; Slurm job 1849087)
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
