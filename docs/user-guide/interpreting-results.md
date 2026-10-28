# Interpreting results

Every claim GIPPA makes carries one of four labels:

- a q-value;
- an estimated false discovery rate (FDR);
- a family-wise error rate (FWER);
- the label `unestimated`.

That is the reporting rule. It does not mean every method is calibrated in every situation. Some
methods ship in 1.0 with a known calibration failure. Their rows carry a flag, and their p- and
q-values are not discovery claims. They are listed in
[Methods shipped flagged uncalibrated](#methods-shipped-flagged-uncalibrated) below. Read that
section before you report any function, determinant, island, operon or variant result.

This page explains what each label means in GIPPA, and what a detection does and does not tell
you. The statistical methods behind each estimate are documented one note per method under
[Methods](../methods/index.md).

## q-values, FDR, FWER and `unestimated`

**FDR.** The false discovery rate is the expected share of false claims among the claims accepted
at a threshold. GIPPA estimates it in two ways:

- from decoys, which are matched false hypotheses:
  - reversed or shuffled genome translations, for placement;
  - convergent and → ← → sandwich pairs, for missed operons;
- from p-values, with Benjamini–Hochberg (BH) adjustment within a pre-declared family of tests.

**q-value.** The q-value of a claim is the smallest FDR at which that claim would be accepted.
Accepting every row with `q_value ≤ 0.05` gives a list with an estimated FDR of at most 5 %. A
q-value is a property of the list, not a probability that one row is wrong.

BH and decoy q-values control the FDR **within one family only**. Nothing is controlled across
families. If you scan several tables, expect roughly α × (number of discoveries) false ones in
each family ([multiple testing](../methods/multiple-testing.md)).

**FWER.** The family-wise error rate is the probability of at least one false claim in the whole
family. GIPPA uses FWER for abundance islands. One family covers the whole run: each window is
compared with the maximum statistic over the windows of **every** scanned replicon. The null comes
from a circular TU-block bootstrap of the detrended TU scores. The FWER claim is conditional on
a dependence regime (see [islands](#islands-the-fwer-holds-in-a-stated-regime) below). It is not
an unconditional guarantee.

**`unestimated`.** GIPPA labels a claim `unestimated` when the input cannot support an error
estimate, for example:

- gene detection from peptide sequences alone;
- a test whose null has too few units.

The claim is still reported, so you can inspect it. Do not report it as error-controlled. Flags
such as `ci_unestimated`, `class_fdr_unestimated`, `not_estimable`, `ora_uncalibrated` and
`null_uncalibrated` say which part is missing or uncalibrated (see the
[flag reference](../reference/flags.md)).

**Families and defaults.** BH applies within the families fixed in code. Here are the main claims,
with the default threshold of each (`stats.alpha.*` in the
[configuration reference](../reference/configuration.md#stats)):

| Claim | Error control | Default |
|---|---|---|
| Peptide placed at a locus | placement q-value from genome decoys: one list per decoy type × stratum (0, 1, 2 differences), length classes pooled | 0.01 |
| Gene detected | gene q-value (tiers `decoy_picked`, `mayu`) or `unestimated` | 0.01 |
| Novel ORF | FDR per class | 0.01 |
| ALT call | the 0.01 is the placement-q filter on supporting peptides (`stats.alpha.variants`); the class error rate is reported, never thresholded, else `unestimated` (see below) | 0.01 |
| Region effect | BH within each (region type, statistic) family | 0.05 |
| Autocorrelation | BH across lags | 0.05 |
| Abundance island | FWER, run-wide over replicons; conditional on the dependence regime | 0.05 |
| Operon discordance | BH over the both-detected pairs (`both_detected_t`); pooled level under smooth, monotone σ(abundance) only | 0.1 |
| Missed operon | decoy FDR | 0.1 |
| Category detection (ORA) | BH per ontology and direction; **every row `ora_uncalibrated`: not a discovery claim** | 0.05 |
| GSEA | BH per collection; conditional on detection, with no p_detect adjustment | 0.05 |
| Determinant model | permutation p-value; **every null p `null_uncalibrated`: not a discovery claim** | 0.05 |

Every p-value and q-value is reported, not only the significant ones. Every test row has the
columns `effect`, `ci_low`, `ci_high`, `method` and `null_model`. Where no CI is estimated,
`ci_low` and `ci_high` are NA and the row carries `ci_unestimated`. This is the case for:

- every GSEA row (NES has no CI, [ADR-0027](../adr/ADR-0027-gsea-nes-without-ci.md));
- every ORA row (`effect` = k − E[K | N], with no CI);
- the autocorrelation lag rows;
- the operon discordance rows;
- the P(miss) rows, whose `effect` is also NA.

The island CI is a different case. It is descriptive and conditional on selection: the window was
chosen as the run-wide maximum. The CI therefore under-covers, essentially never contains zero,
and its coverage has not been measured. The island's `p_value` is the only test.

**ALT calls.** The class error rate estimates the share of erroneous reads among ALT-supporting
de novo reads whose spectrum a database search identifies confidently. It is not the false-ALT
rate among all ALT reads. It is biased upwards, up to a degenerate value of 1, when the reference
search could not contain the ALT peptide. It is biased downwards by spectrum-quality selection.
Read a rate of 1 from a reference-proteome search as "not measured"
([variants](../methods/variants.md)).

## Methods shipped flagged uncalibrated

These results are reported, but their p- and q-values are **not error-controlled**. Do not report
them as discoveries.

- **Detection enrichment (ORA): `ora_uncalibrated` on every row.** The detection-probability
  recalibration is anti-conservative for categories selected on TU size or detection probability
  ([ADR-0035](../adr/ADR-0035-ora-tu-detection-recalibration.md),
  [ADR-0037](../adr/ADR-0037-ora-flexible-tu-recalibration.md);
  [function](../methods/function.md)). A skipped fit (`tu_recal_skipped`) is a second reason for
  the flag.
- **Determinant model null: `null_uncalibrated` on every null p.** The permutation p is valid only
  when TUs are exchangeable. It fails under autocorrelation along the genome, positional covariates
  and detection selection ([ADR-0034](../adr/ADR-0034-determinants-null-uncalibrated.md)).
- **Operon discordance outside smooth noise.** The discordance level holds under smooth, monotone
  σ(abundance). With a step change in noise at some abundance (the D-1b-step attack), it is
  exceeded, and no row is flagged
  ([ADR-0055 Amendment 1](../adr/ADR-0055-operon-discordance-abundance-noise.md)).
- **Islands at strong, short-resolved autocorrelation.** The FWER breaches its level at TU-level
  ρ = 0.9 when the bootstrap block is short. No row is flagged (see the next section).
- **Variant-call precision (gate M7).** The simulated precision gate G1 and the de novo error
  attacks A6 are strict expected failures. Gate M7 is not recorded as passed. Mass-equivalent and
  swap ALTs stay `high` with the flag `mass_equivalent`
  ([ADR-0049](../adr/ADR-0049-m7-variant-precision.md)).

## Islands: the FWER holds in a stated regime

The island FWER is claimed for these cases
([islands › Validity domain](../methods/islands.md#validity-domain-adr-0021-amendment-1),
[ADR-0021](../adr/ADR-0021-island-oriter-block-length.md) Amendment 1):

- TU-level autocorrelation ρ ≤ 0.5;
- ρ = 0.9 on a chromosome of about 4,500 genes, where the block length is long enough.

It is breached at ρ = 0.9 when the block length ℓ′ is 7 or less (0.095 to 0.208 at α = 0.05).
It does not apply at all when `context.islands.block_len` is fixed. No per-row flag marks a
breach, because no statistic reliably detects one. The report instead prints, per replicon:

- the estimated dependence length ℓ̂ (0 means "not estimated", never "no autocorrelation");
- the block length `block_len`;
- the lag-1 Moran's I of the detrended TU scores.

Use these to judge where your data sit.

## Operons: discordance and P(miss)

- **Discordance** (`test = both_detected_t` in `pairs.tsv`) is one BH family over the pairs with
  both genes detected. It holds level only under smooth, monotone σ(abundance), and it has little
  power at realistic effect sizes ([operons](../methods/operons.md)).
- **P(miss)** (`test = one_detected_pmiss`) is the model's probability that this partner would be
  missed if it were co-expressed. It is not a significance level. These rows have `q_value` NA,
  sit outside the BH family, and never make a discordance call
  ([ADR-0056](../adr/ADR-0056-pmiss-outside-the-bh-family.md)). P(miss) is anti-conservative when
  detectability varies between proteins in ways the model does not see. Thresholding it at 0.05
  flagged 0.198 of one-detected rows in a simulated null with such variation, where about 0.05 was
  expected.

## GSEA is conditional on detection

Abundance GSEA ranks detected genes only. Its results are **conditional on detection, with no
p_detect adjustment** ([ADR-0026](../adr/ADR-0026-gsea-conditional-on-detection.md)): when a set's
genes are harder to detect, its detected members are selected for high abundance, and the set can
look enriched without any true difference. Islands and operon discordance likewise carry no
p_detect adjustment.

- "No significant set" never means "no enrichment".
- `detection_confounded` and `genomic_cluster` are screens, not guarantees: an unflagged row is
  not certified ([function GSEA › Limits](../methods/function-gsea.md#limits)).
- In a collection that has flagged rows, no FDR claim is made for the unflagged rows. Flagged rows
  enter the same BH family, and their effect on the other rows has not been re-measured for the
  R engine.

## Detection tiers

How GIPPA decides that a gene is detected depends on what the input contains. The tier is chosen
once per run and recorded in the manifest and in `genes.tsv` (`detection_tier`). It is the first
of the four below that the input allows:

| Tier | When | Detected means | Gene error |
|---|---|---|---|
| `decoy_picked` | decoy rows with `protein_ids` exist, at least 80 % of them pair with a gene, and the decoys share a ranking column (q-value, PEP or score) with the targets | gene q ≤ 0.01 by picked target–decoy competition over genes | gene q-value |
| `mayu` | otherwise, the accepted peptides carry q-values or PEPs; also chosen when decoys exist but pair too rarely (warning `decoy_pairing_low`) | gene q ≤ 0.01 by a MAYU-style estimate of false genes | gene q-value |
| `two_peptide` | otherwise, peptide sequences only | at least 2 distinct accepted peptides | `unestimated` |
| `protein_list` | proteins only, no peptides | a mapped input protein | `unestimated` (inherited from the search engine) |

Decoys without a ranking column cannot be picked. Preflight warns about them, and the tier falls
back.

Under `two_peptide`, a gene with only one peptide gets `detection_status =
single_peptide_unestimated`. It counts as not detected in the primary analysis and as detected in
the sensitivity rerun. Partial genes at linear contig ends are `partial_excluded`.
`detected_2pep` (at least 2 distinct accepted peptides) is reported in every tier.

To get the strongest tier, keep the search engine's decoy rows (see
[Inputs](inputs.md#peptide-table-peptidestsv)). Details:
[detection tiers](../methods/detection-tiers.md) and
[MAYU gene FDR](../methods/mayu-gene-fdr.md).

## Absence of detection is not absence of expression

GIPPA says "detected", never "expressed". A gene can go undetected for reasons unrelated to
whether it is made. The protein may be:

- short, with few observable peptides;
- membrane-bound, or hard to ionise;
- below the instrument's detection limit;
- absent from the search database.

To account for this, GIPPA fits a detectability model. It predicts each gene's probability of
detection (`p_detect`) from protein covariates, with cross-fitting. Region O − E tests and ORA
compare observed detections with the sum of `p_detect`, never with raw counts. A region "detected
less than expected" is detected less than its proteins' detectability predicts. That is a
statement about detection, not about transcription or translation. GSEA, islands and operon
discordance work on abundances and have no p_detect adjustment (see above).

**`detectability_constant`.** The model is replaced by a constant in any of these cases:

- fewer than 50 detected eligible genes;
- fewer than 50 undetected eligible genes;
- fewer TUs than cross-fitting folds.

`p_detect` is then the overall detection rate for every gene, and every gene carries
`detectability_constant`. The detection tests still run, but **without any adjustment for
detectability**, and ORA uses ω = 1. Treat their results as unadjusted comparisons. "No gene
detected" is only an extreme case of this, where the detection-dependent tests are
`not_estimable`.

## Placement FDR is not search-engine FDR

Two error rates apply to a peptide, and they answer different questions.

- **The search engine's FDR** (the input `q_value`, `pep`) asks: is this spectrum matched to the
  right peptide sequence? GIPPA inherits it and does not re-estimate it. Engine precursor q-values
  may themselves under-control FDR, and they do not control gene-level error (GIPPA_SPEC.md §8.7).
  This is why gene detection uses the gene-level tiers. GIPPA's pilot validation measures the
  realised false discovery proportion by entrapment.
- **GIPPA's placement FDR** (`placement_q` in `peptide_hits.tsv`) asks: given the peptide
  sequence, is it placed on the right locus of this genome? A short peptide can match a genome
  window by chance, and more so when mismatches are allowed. GIPPA estimates this rate by
  searching the same peptides against decoy translations of the genome (reversed frames, or
  shuffled frames for de novo input). There is one list per decoy type and per stratum of 0, 1 or
  2 differences, with length classes pooled.

A peptide with a search q-value of 0.001 can still have a placement q-value of 1, for example when
it matches the genome only with two substitutions. A peptide counts as evidence for a gene only
when all of these hold:

- it passes the search threshold (`q_value ≤ input.q_value_max`, when a q-value exists);
- its placement is accepted (`placement_q ≤ stats.alpha.placement`, or `same_db`);
- its locus class assigns it to the gene;
- its `uniqueness` is `unique` or `multi_identical`;
- under `--genome-relation same`, it is an exact hit.

Peptides rejected by placement FDR are in `mapping/unmapped.tsv` with
`reason = rejected_placement_fdr`.

`mapping/mapping_fdr.tsv` lists the targets, decoys, π̂_wrong and FDR per length class and number
of differences. Its per-length-class rows are descriptive. They are **not** the q-values: the
q-values come from one pooled list per decoy type and stratum. Details:
[placement FDR › Strata and ranking](../methods/placement-fdr.md#strata-and-ranking-m3-10).

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
  the codon's 1-based position as `readthrough_codon`. Only annotated terminal stops and the
  annotated recoded codons below can be read this way; every other stop codon stays unmatchable.
  Read-through rows are part of the ordinary placement FDR of their stratum. That FDR does not
  control the read-through subset on its own.
- `sec_site` and `pyl_site`: the peptide consumed an annotated recoded codon inside a gene,
  selenocysteine (`sec_site`) or pyrrolysine (`pyl_site`). Only codons the annotation itself
  recodes count: a `U` in the annotated protein over a TGA codon, or an `O` over a TAG codon, in
  a CDS or sORF (not the start codon, and never in a pseudogene). A `U` or `O` over any other
  codon is skipped. Every other internal TGA or TAG stays unmatchable.
  - The peptide consumes the codon when its alignment holds the whole codon, in frame, with one
    residue standing for it (cost 1, recorded as a `substitution`). Frameshift placements never
    consume one.
  - The flag follows the codon, not the residue the search engine reported: TGA gives `sec_site`
    and TAG gives `pyl_site`, so a peptide that reads `W` over a selenocysteine TGA still gets
    `sec_site`. Such a row gets `sec_site` or `pyl_site` instead of `readthrough`. A peptide that
    consumes both a terminal stop and a recoded codon gets both flags.
  - `mappings.gff3` gives the codon's 1-based position as `recoded_codon` (comma-separated if
    there are several), in the same way as `readthrough_codon`.
  - Like read-through rows, these rows are part of the ordinary placement FDR of their stratum,
    and that FDR does not control the `sec_site` or `pyl_site` subset on its own. The recoded
    codons are mirrored in every decoy, so they do not bias the list FDR.
- `same_db`: an exact hit from same-strain input, accepted without the decoy step.
- `low_complexity`: the peptide has at most 3 distinct residues, with I and L counted as one.
  For such peptides, the chance-hit score E can underestimate chance hits by 1.9× to 15×. The
  decoys keep the q-value valid, and the flag changes nothing in it
  ([ADR-0011](../adr/ADR-0011-low-complexity-placement-flag.md)).
- `origin_spanning`: the placement crosses the origin of a circular replicon.
- `multi`, `multi_identical`: the peptide has accepted loci in more than one locus class, or in
  identical gene copies. `uniqueness` (`unique`, `multi_identical`, `shared`) says whether the
  peptide supports one gene.
- `engine_protein_mismatch`: with same-strain input, the gene GIPPA chose is not among the
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

An empty result is not evidence of absence. Several tests have little or no power by
construction:

- **Low islands.** At normal detection rates, the low direction is a structural zero. Even a fully
  undetected window of 20 TUs cannot reach the critical value. No low island is therefore not
  evidence that no silenced region exists, unless fewer than about a quarter of genes are
  undetected ([islands](../methods/islands.md#the-low-direction-is-a-structural-zero-at-normal-detection-rates)).
- **Missed operons.** With the shipped decoys (convergent and sandwich pairs), this test made no
  call in the simulation gate or in any attack. It is near-powerless by construction, so an empty
  list says nothing about missed operons ([operons](../methods/operons.md#power-what-this-family-can-and-cannot-find)).
- **Operon discordance.** Power is near zero at realistic effect sizes, so few or no flagged pairs
  is the expected outcome even when discordant pairs exist.
- **GSEA.** "No significant set" never means "no enrichment" (see
  [GSEA is conditional on detection](#gsea-is-conditional-on-detection)).

**A file that is missing** usually means its module was switched off. `run_manifest.json` lists
switched-off modules under `modules_off`, and they write no files. A missing tool layer (for
example, geNOMAD without its database) appears as the run flag `layer_missing:<layer>` in
`logs/events.jsonl`, and the tests that need it are skipped or marked. The plots of a switched-off
module or missing layer are skipped with a reason (`module_off:<module>`,
`layer_missing:<layer>`).

## Wording

GIPPA's own text follows one rule: "detected", never "expressed". The report, the plot titles and
the column names all use it. Every number in the report that comes from a statistic carries its
error estimate or `unestimated`. Please keep the same wording when you write about your results.
