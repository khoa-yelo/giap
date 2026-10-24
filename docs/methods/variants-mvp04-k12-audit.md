# Variants on the MVP-04 K-12 run: flag audit (statistician, 2026-10-24)

Run `local/runs/mvp04/20261023_a1` (read-only), inputs `local/mvp04/20261023` (giapsim
`real_genome`, seed 20261023). Aggregates only.

## Inputs are not "K-12 against K-12"
`provenance.json`: giapsim derives a *sample* strain from MG1655 with 2,332 planted SNVs
(1,488 non-synonymous gene codons in `truth/alt_sites.tsv`) and 1,200 planted novel ORFs. The
peptides are a simulated database search against the sample's proteome: all 52,278 rows have
`source = database`, `search_db = other`; the config sets `genome_relation: related`.

## Flags: expected, not a bug
`flags_for` (src/giap/variants/flags.py:189-191) implements §7.5.3 literally: `strain_difference`
= any ALT under `related`; `database_derived` = an ALT supported by a `(database, other)` peptide.
Both are therefore true for every ALT of this run, and MVP-04 gate 11 requires exactly that
(`n_rows_missing_flags = 0`). Flags are not uniform: 128/428 `mass_equivalent`, 35 `start_codon`,
0 `multi_nt` (the sample carries SNVs only). With a single-source input both flags are constant
and carry no information; they discriminate only in mixed inputs. On a real same-strain K-12
search (`search_db = same`, `genome_relation = same`) neither flag would appear.

## The 428 ALT/mixed rows (70 high, 358 candidate)
| class | high | candidate |
|---|---|---|
| truth codon (sample differs from genome) | 50 | 198 |
| real sample SNV inside a planted novel ORF (no gene; not in `alt_sites.tsv`) | 5 | 7 |
| ALT M at an annotated GTG/TTG start (`aa_pos` 1, `start_codon`) | 6 | 29 |
| ALT M at a planted novel-ORF non-ATG start (no gene, no flag) | 9 | 118 |
| genome = sample at the codon: ALT peptide from one IS3 / IS421 copy piled up at its identical copies (`mixed`, n_alt 1) | 0 | 6 |

The only ALTs not explained by the sample are the 6 repeat-copy candidates (6/428 = 0.014,
Wilson 95 % [0.006, 0.030]); no `high` row is false. I/L: the 30 `J` rows are the folded I/L
ALT token (M7-02 rule 5, M7-18), not a residue error; 29 of 30 are truth codons.

## Gate 11 (MVP-04) scoring
`scripts/dev/mvp04_summary.py` `gate11` now gives precision 40/70 = 0.571 [0.455, 0.681] and
recall 40/60 = 0.667 [0.541, 0.773]: fail. The causes are in the scorer and the truth table:
1. ALT `J` vs truth `I`/`L` is compared literally (10 `high` rows). Folded: 50/70.
2. Initiator M at annotated or novel non-ATG starts counted as false (15 `high`); giapsim itself
   drops initiator codons from truth. Excluding them: 50/55 = 0.909 [0.804, 0.961].
3. Truth covers annotated genes only; 5 `high` rows are real SNVs in novel ORFs. Adding them:
   55/55 = 1.0 [0.935, 1.0].
4. Recall: `n_true_peptides` counts true peptides shorter than `mapping.tolerant.min_len` (10),
   which can never support an ALT. With folding and length >= 10: 50/53 = 0.943 [0.846, 0.981].

## Display
V12 panel (b) `_header` (src/giap/viz/variants.py:244) formats `gene` and `aa_pos` without an NA
case, so no-gene sites print "None aa None". Spec M9-18 decision 11 covers only the peptide line.
