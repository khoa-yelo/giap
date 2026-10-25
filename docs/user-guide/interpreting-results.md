# Interpreting results

Every claim GIAP makes carries one of four labels:

- a q-value;
- an estimated false discovery rate (FDR);
- a family-wise error rate (FWER);
- the label `unestimated`.

Nothing is reported as trusted without one of these. This page explains what each means in GIAP,
and what a detection does and does not tell you. The statistical methods behind each estimate
are documented one note per method under [Methods](../methods/index.md).

## q-values, FDR, FWER and `unestimated`

**FDR.** The false discovery rate is the expected share of false claims among the claims accepted
at a threshold. GIAP estimates it in two ways:

- from decoys, which are matched false hypotheses such as reversed or shuffled genome
  translations, or convergent gene pairs for operons;
- from p-values, with Benjamini–Hochberg (BH) adjustment within a pre-declared family of tests.

**q-value.** The q-value of a claim is the smallest FDR at which that claim would be accepted.
Accepting every row with `q_value ≤ 0.05` gives a list with an estimated FDR of at most 5 %. A
q-value is a property of the list, not a probability that one row is wrong.

**FWER.** The family-wise error rate is the probability of at least one false claim in the whole
family. GIAP uses FWER for abundance islands: it compares each window with the maximum statistic
over circular block-bootstrap resamples of the replicon. An island at FWER 0.05 is a stronger statement than a q-value of 0.05.

**`unestimated`.** GIAP labels a claim `unestimated` when the input cannot support an error
estimate, for example:

- gene detection from peptide sequences alone;
- a test whose null has too few units.

The claim is still reported, so you can use it, but no error rate stands behind it. Flags such as
`ci_unestimated`, `class_fdr_unestimated` or `not_estimable` say which part is missing (see the
[flag reference](../reference/flags.md)).

**Families and defaults.** BH applies within the families fixed in code. Here are the main claims,
with the default threshold of each (`stats.alpha.*` in the
[configuration reference](../reference/configuration.md#stats)):

| Claim | Error control | Default |
|---|---|---|
| Peptide placed at a locus | placement q-value, from genome decoys, per mismatch stratum | 0.01 |
| Gene detected | gene q-value (tiers `decoy_picked`, `mayu`) or `unestimated` | 0.01 |
| Novel ORF | FDR per class | 0.01 |
| ALT call | class error rate, else `unestimated` | 0.01 |
| Region effect | BH within region type | 0.05 |
| Autocorrelation | BH across lags | 0.05 |
| Abundance island | FWER | 0.05 |
| Operon discordance, missed operon | BH; decoy FDR | 0.1 |
| Category detection, GSEA | BH per ontology and direction; per collection | 0.05 |
| Determinant model | permutation p-value | 0.05 |

Every p-value and q-value is reported, not only the significant ones. Every test row also carries
an effect size with a 95 % confidence interval (`effect`, `ci_low`, `ci_high`), the method and
the null model.

## Detection tiers

How GIAP decides that a gene is detected depends on what the input contains. The tier is chosen
once per run and recorded in the manifest and in `genes.tsv` (`detection_tier`). It is the first
of the four below that the input allows:

| Tier | When | Detected means | Gene error |
|---|---|---|---|
| `decoy_picked` | decoy rows with `protein_ids` exist, and enough of them pair with a gene | gene q ≤ 0.01 by picked target–decoy competition over genes | gene q-value |
| `mayu` | otherwise, the accepted peptides carry q-values or PEPs | gene q ≤ 0.01 by a MAYU-style estimate of false genes | gene q-value |
| `two_peptide` | otherwise, peptide sequences only | at least 2 distinct accepted peptides | `unestimated` |
| `protein_list` | proteins only, no peptides | a mapped input protein | `unestimated` (inherited from the search engine) |

Under `two_peptide`, a gene with only one peptide gets `detection_status =
single_peptide_unestimated`. It counts as not detected in the primary analysis and as detected in
the sensitivity rerun. Partial genes at linear contig ends are `partial_excluded`.
`detected_2pep` (at least 2 distinct accepted peptides) is reported in every tier.

To get the strongest tier, keep the search engine's decoy rows (see
[Inputs](inputs.md#peptide-table-peptidestsv)). Details:
[detection tiers](../methods/detection-tiers.md) and
[MAYU gene FDR](../methods/mayu-gene-fdr.md).

## Absence of detection is not absence of expression

GIAP says "detected", never "expressed". A gene can go undetected for reasons unrelated to
whether it is made. The protein may be:

- short, with few observable peptides;
- membrane-bound, or hard to ionise;
- below the instrument's detection limit;
- absent from the search database.

To account for this, GIAP fits a detectability model. It predicts each gene's probability of
detection (`p_detect`) from protein covariates, with cross-fitting. Detection tests then compare
observed detections with the sum of `p_detect`, never with raw counts. A region "detected less
than expected" is detected less than its proteins' detectability predicts. That is a statement
about detection, not about transcription or translation.

When detectability cannot be modelled, every gene carries `detectability_constant`. This happens,
for example, when no gene is detected. Detection-dependent tests are then not estimable.

## Placement FDR is not search-engine FDR

Two error rates apply to a peptide, and they answer different questions.

- **The search engine's FDR** (the input `q_value`, `pep`) asks: is this spectrum matched to the
  right peptide sequence? GIAP inherits it and does not re-estimate it.
- **GIAP's placement FDR** (`placement_q` in `peptide_hits.tsv`) asks: given the peptide
  sequence, is it placed on the right locus of this genome? A short peptide can match a genome
  window by chance, and more so when mismatches are allowed. GIAP estimates this rate by
  searching the same peptides against decoy translations of the genome (reversed frames, or
  shuffled frames for de novo input). It works within strata of 0, 1 and 2 differences.

A peptide with a search q-value of 0.001 can still have a placement q-value of 1, for example when
it matches the genome only with two substitutions. A peptide counts as evidence for a gene only when all of
these hold:

- it passes the search threshold (`q_value ≤ input.q_value_max`, when a q-value exists);
- its placement is accepted (`placement_q ≤ stats.alpha.placement`, or `same_db`);
- its locus class assigns it to the gene;
- its `uniqueness` is `unique` or `multi_identical`;
- under `--genome-relation same`, it is an exact hit. Peptides rejected by placement FDR are in `mapping/unmapped.tsv` with
`reason = rejected_placement_fdr`. `mapping/mapping_fdr.tsv` gives, per length class and number
of differences, the targets, decoys, π̂_wrong and stratum FDR. Details:
[placement FDR](../methods/placement-fdr.md).

**Same-strain input.** Peptides searched against this genome's own proteins (`search_db = same`)
skip the decoy step for exact hits: `placement_q = 0` and the flag `same_db`. Tolerant hits still
go through it.

## Kinds of placement

Each row of `mapping/peptide_hits.tsv` is one peptide at one accepted locus. The columns `route`,
`n_diff`, `diff_types` and `flags` say how the peptide was placed.

| Kind | How to recognise it | What it means |
|---|---|---|
| Exact | `route = exact`, `n_diff = 0`, empty `diff_types` | the peptide occurs verbatim in a six-frame translation (I and L treated as equal) |
| Tolerant | `route = tolerant`, `n_diff` 1 or 2, `diff_types` from `substitution`, `swap`, `block`, `near_block` | the peptide matches with up to 1 difference at 10–11 aa, or up to 2 at 12–40 aa; peptides shorter than 10 aa are never placed this way |
| Frameshift | `route = tolerant`, `diff_types` contains `frameshift` | a ±1-nt frameshift: the prefix lies in one frame and the suffix in the adjacent frame, with at least 4 aa on each side |

The difference types count one each:

- `substitution`: one residue differs;
- `swap`: two adjacent residues swapped;
- `block`: a mass-equivalent exchange, such as N for GG or Q for AG;
- `near_block`: a near-equal mass exchange, such as W for AD or EG;
- `frameshift`: the frame shift itself.

Tolerant placements are judged at π̂_wrong = 1, which is conservative. They carry the flag
`strain_difference` when `--genome-relation related`, and genes with such evidence carry
`tolerant_evidence`.

In `mapping/mappings.gff3`, a frameshift placement has two `match_part` children, one per frame.
It also has the attribute `frameshift_reread` (the base read twice) or `frameshift_skipped` (the
base skipped), giving the 1-based position of that base. In `mappings.bed`, the two parts are
BED12 blocks.

**Flags on placements** (`flags` column; full list in the
[flag reference](../reference/flags.md#peptide_hits)):

- `readthrough`: the tolerant alignment substituted a residue for the annotated stop codon at the
  end of a CDS or sORF (cost 1, recorded as a `substitution`). This is evidence of stop-codon
  read-through, or a C-terminal extension (locus class `cterm_extension`). `mappings.gff3` gives
  the codon's 1-based position as `readthrough_codon`. Only annotated terminal stops can be read
  this way; every other stop codon stays unmatchable. Read-through rows are part of the ordinary
  placement FDR of their stratum. That FDR does not control the read-through subset on its own.
- `sec_site` and `pyl_site` mark peptides that consume a recoded TGA (selenocysteine, U) or TAG
  (pyrrolysine, O) codon. They are specified in M7-16 (part 16c) but are not in the 1.0
  `peptide_hits` schema. A 1.0 run does not write them.
- `same_db`: an exact hit from same-strain input, accepted without the decoy step.
- `low_complexity`: the peptide has at most 3 distinct residues. Its chance-hit score is less
  reliable, but its q-value is unchanged
  ([ADR-0011](../adr/ADR-0011-low-complexity-placement-flag.md)).
- `origin_spanning`: the placement crosses the origin of a circular replicon.
- `multi`, `multi_identical`: the peptide has accepted loci in more than one locus class, or in
  identical gene copies. `uniqueness` (`unique`, `multi_identical`, `shared`) says whether the
  peptide supports one gene.
- `engine_protein_mismatch`: with same-strain input, the gene GIAP chose is not among the
  engine's protein IDs.
- `nterm_extension_evidence`: the peptide lies upstream of the annotated start, in frame.

## Empty results

**A table with only a header** means the step ran and found nothing at its threshold. For
example, a header-only `tables/novel_orfs.tsv` means no novel ORF or start-site candidate passed
its class FDR. It is a result, not a failure. Also check:

- that the input allowed the test;
- the flags on the rows and the run, such as `class_fdr_unestimated`.

Likewise, an empty `islands.tsv` or `operon_flags.tsv` means no window or pair passed. If no gene
is detected, the detection-dependent tests are `not_estimable`.

**A file that is missing** usually means its module was switched off. `run_manifest.json` lists
switched-off modules under `modules_off`, and they write no files. A missing tool layer (for
example, geNOMAD without its database) appears as the run flag `layer_missing:<layer>` in
`logs/events.jsonl`, and the tests that need it are skipped or marked. The plots of a switched-off
module or missing layer are skipped with a reason (`module_off:<module>`,
`layer_missing:<layer>`).

## Wording

GIAP's own text follows one rule: "detected", never "expressed". The report, the plot titles and
the column names all use it. Every number in the report that comes from a statistic carries its
error estimate or `unestimated`. Please keep the same wording when you write about your results.
