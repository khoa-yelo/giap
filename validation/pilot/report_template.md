# GIAP pilot report, pass 1 (arms A1 and A2)

Rendered by `validation/pilot/metrics.py render` from `metrics.json` (docs/specs/M10-05-pilot-metrics.md).
Every value below is followed by the file and key it comes from; none is typed by hand.
This template holds no numbers. Only aggregates appear here: no sequence, no gene list, no per-gene value.

## 1. Status

| §15.6 bullet | Status |
|---|---|
| GIAP runs end to end on A1, A2 and A3 | {{checklist.0.status}} |
| ≥ 99.5 % of A2 peptides map exactly (TE1) | {{checklist.1.status}} |
| A1 and A2 detections and abundances are compared | {{checklist.2.status}} |
| Measured runtimes replace the estimates | {{checklist.3.status}} |
| TE2–TE5 are reported; failures explained | {{checklist.4.status}} |

Spec SHA: {{spec}}.

## 2. Data and arms

PXD078459 wild-type runs (M10-01 manifest), DIA-NN 2.3.1 Academia (ADR-0051). A1 searched UniProt
K-12 (UP000000625) plus contaminants; A2 searched the Nissle strain database (`strain_db.faa`) plus
contaminants. Job IDs of M10-02/03/04 are in PROGRESS.md and the `SACCT` table.

## 3. Metrics

### Metric 1: identifications at 1 % FDR (union over runs)

Only the precursor count is at a controlled 1 % FDR (global precursor q-value). Peptide counts are
not FDR-controlled exactly. Protein groups are groups with ≥ 1 accepted precursor (no protein-level FDR control).
Per-run counts are "accepted at global q ≤ 0.01 among the rows DIA-NN reported for that run", not a per-run FDR.

| | A1 | A2 |
|---|---|---|
| precursors | {{m1.A1.union.precursors.all}} | {{m1.A2.union.precursors.all}} |
| peptides (all) | {{m1.A1.union.peptides.all}} | {{m1.A2.union.peptides.all}} |
| peptides (non-contaminant) | {{m1.A1.union.peptides.non_contaminant}} | {{m1.A2.union.peptides.non_contaminant}} |
| peptides (contaminant-only) | {{m1.A1.union.peptides.contaminant_only}} | {{m1.A2.union.peptides.contaminant_only}} |
| protein groups | {{m1.A1.union.protein_groups.all}} | {{m1.A2.union.protein_groups.all}} |
| protein groups with ≥ 2 accepted peptides | {{m1.A1.union.protein_groups_min2.all}} | {{m1.A2.union.protein_groups_min2.all}} |

### Metric 2: mapping shares (non-contaminant query peptides)

| | A1 | A2 |
|---|---|---|
| query peptides (U) | {{m2.A1.n_universe}} | {{m2.A2.n_universe}} |
| contaminant-only (C) | {{m2.A1.n_contaminant_only}} | {{m2.A2.n_contaminant_only}} |
| denominator (D) | {{m2.A1.n_denominator}} | {{m2.A2.n_denominator}} |
| exact | {{m2.A1.shares.exact}} | {{m2.A2.shares.exact}} |
| tolerant only | {{m2.A1.shares.tolerant}} | {{m2.A2.shares.tolerant}} |
| unmapped | {{m2.A1.shares.unmapped}} | {{m2.A2.shares.unmapped}} |

Genome-inconsistent peptides in A1 (`detection_summary.tsv`): rate {{m2.A1.genome_inconsistent.value}}, count {{m2.A1.genome_inconsistent.count}}; equal to tolerant + unmapped: {{m2.A1.genome_inconsistent.consistent}}.

### Metric 3: coverage of the strain CDS

Shares with nominal 95 % Wilson intervals, descriptive only. A1 and A2 are paired on the same genes,
so the arms are not compared by interval overlap; hypothesis 1 is the paired comparison.
"Strain-specific" is defined by peptide identity (no shared in-silico tryptic peptide with the K-12
digest, I/L collapsed), a proxy for "no K-12 counterpart"; it does not mean "undetectable by A1",
which can still detect such genes through tolerant placements.

| subset | A1 | A2 |
|---|---|---|
| all | {{m3.A1.subsets.all.share}} | {{m3.A2.subsets.all.share}} |
| chromosome | {{m3.A1.subsets.chromosome.share}} | {{m3.A2.subsets.chromosome.share}} |
| plasmid | {{m3.A1.subsets.plasmid.share}} | {{m3.A2.subsets.plasmid.share}} |
| prophage | {{m3.A1.subsets.prophage.share}} | {{m3.A2.subsets.prophage.share}} |
| BGC | {{m3.A1.subsets.bgc.share}} | {{m3.A2.subsets.bgc.share}} |
| strain-specific | {{m3.A1.subsets.strain_specific.share}} | {{m3.A2.subsets.strain_specific.share}} |
| K-12 shared | {{m3.A1.subsets.k12_shared.share}} | {{m3.A2.subsets.k12_shared.share}} |

### Metric 4: quantitative agreement (shared genes, from `giap compare`)

Spearman rho {{m4.spearman_rho}} [{{m4.ci_low}}, {{m4.ci_high}}] on {{m4.n_shared}} shared genes (TU-block bootstrap); Bland–Altman bias {{m4.ba_bias_log2}} log2; flags: {{m4.flags}}.

### Metric 5: stability of conclusions (Jaccard, from `giap compare`)

| family | Jaccard | note |
|---|---|---|
| islands | {{m5.islands.jaccard}} | {{m5.islands.note}} |
| operon flags | {{m5.operon_flags.jaccard}} | {{m5.operon_flags.note}} |
| enriched categories | {{m5.enriched_categories.jaccard}} | {{m5.enriched_categories.note}} |

### Metric 7: runtime and peak memory (pilot measurements on a shared cluster, not the benchmark)

| | A1 | A2 |
|---|---|---|
| wall time | {{m7.A1.wall_time_s}} | {{m7.A2.wall_time_s}} |
| sum of step wall times | {{m7.A1.sum_task_wall_s}} | {{m7.A2.sum_task_wall_s}} |
| max step peak RSS (one task, not the arm's total memory) | {{m7.A1.max_step_peak_rss_bytes}} | {{m7.A2.max_step_peak_rss_bytes}} |
| complete (no cached step; else totals are lower bounds) | {{m7.A1.runtime_complete}} | {{m7.A2.runtime_complete}} |

<!-- if m8.peptide -->
### Metric 8: A1's free entrapment (from `entrapment.py`, M10-07)

The entrapment is free: it is the K-12-only part of A1's own database, the K-12 peptides absent from
Nissle's translation after the I/L and mass-equivalent rules. Estimator: Wen et al. 2025 (combined
method). r = K-12-only / shared in-silico peptides = {{m8.universe.r}} ({{m8.universe.k12_only}} K-12-only, {{m8.universe.shared}} shared); at precursor level
r_prec = {{m8.universe.r_prec}} (from the configured charge and m/z bounds, not the library).
The assessment is made only at precursor level, the unit DIA-NN's global q-value controls; the peptide level is §15.5's number and is not assessed.

| | peptide (§15.5) | precursor (assessed) |
|---|---|---|
| N_T (shared) | {{m8.peptide.N_T}} | {{m8.precursor.N_T}} |
| N_E (K-12-only) | {{m8.peptide.N_E}} | {{m8.precursor.N_E}} |
| excluded as contaminant | {{m8.peptide.excluded_contaminant}} | {{m8.precursor.excluded_contaminant}} |
| not in the in-silico universe | {{m8.peptide.not_in_universe}} | {{m8.precursor.not_in_universe}} |
| `fdp_combined` (estimated upper bound in expectation) | {{m8.peptide.fdp_combined}} | {{m8.precursor.fdp_combined}} |
| `fdp_lower` (realised lower bound, conditional on the caveats) | {{m8.peptide.fdp_lower}} | {{m8.precursor.fdp_lower}} |
| with 95 % confidence, under equal chance, the FDP is at most | {{m8.peptide.ci.ucb95}} | {{m8.precursor.ci.ucb95}} |
| assessment | {{m8.peptide.assessment}} | {{m8.precursor.assessment}} |

Caveats:

- {{m8.caveats.0}}
- {{m8.caveats.1}}
- {{m8.caveats.2}}
- {{m8.caveats.3}}

<!-- end m8.peptide -->
## 4. TE1 on real data (A2)

Every A2 target identification, true or false, is a substring of an annotated protein, so TE1 checks
that database, genome and mapping agree; it says nothing about identification error.
The expected misses are structural (N-terminal products of non-ATG starts, Sec sites); anything else is a mapping defect.

Share exact in frame at an own gene: {{te1.share_exact_inframe}} of {{te1.n_peptides}} peptides, nominal 95 % Wilson interval [{{te1.ci_low}}, {{te1.ci_high}}] (misses cluster by gene, so its coverage is optimistic; the interval is not a test).
Pass (point estimate at or above the threshold): {{te1.passed}}; lower end below the threshold: {{te1.ci_low_below_threshold}}.

Misses: non-ATG N-terminal {{te1.misses.non_atg_nterm}}, exact at another gene {{te1.misses.exact_other_gene}}, tolerant only {{te1.misses.tolerant_only}}; masked exact in frame {{te1.exact_inframe_masked}}.

Secondary, never the verdict: excluding non-ATG products {{te1.share_excluding_non_atg}}; without peptides whose accepted rows are all at the cut-off {{te1.share_without_q_at_cutoff}}.
A1 diagnostic (exact in frame at any gene, no pass flag): {{te1_a1_diagnostic.share}}.

## 5. Hypotheses

Test family for this issue: {H1}, size 1, no BH.

**H1** (A2 detects more strain CDS than A1, mostly plasmid, prophage and BGC genes): verdict {{h1.verdict}}; failed clauses: {{h1.failed_clauses}}.
A verdict other than "supported" means not supported by the pre-registered rule.
A2-only share of discordant genes {{h1.prop_a2_only}}; coverage difference {{h1.coverage_diff}}.
TU-level test: {{h1.b_tu}} against {{h1.c_tu}} units ({{h1.mixed_tu}} mixed, left out), p_tu {{h1.p_tu}}, unit effect {{h1.unit_effect}}.
Gene-level p_gene {{h1.p_gene}} (assumes independent genes; decides nothing).
Share of A2-only genes in plasmid, prophage or BGC {{h1.mostly}} (base rate {{h1.mostly_base_rate}}).
The arms are two database searches of the same spectra, so discordance is largely fixed by sequence
divergence rather than sampled, and dependence at the scale of a region remains; p is not a replication probability.

**H3** (context conclusions agree between A1 and A2): verdict {{h3.verdict}}; rho criterion {{h3.criteria.rho.status}}; islands {{h3.criteria.islands.status}}; operon flags {{h3.criteria.operon_flags.status}}; enriched categories {{h3.criteria.enriched_categories.status}}.

<!-- if m8.peptide -->
**H2** (A1's genome-inconsistent FDP is at least its entrapment FDP): deferred (it needs the full validation, FDRBench); shown side by side only: A1 genome-inconsistent share {{m2.A1.genome_inconsistent.value}}, free-entrapment `fdp_combined` {{m8.peptide.fdp_combined}}.
**H4**: deferred (it needs A3 and M10-10).
<!-- end m8.peptide -->
<!-- ifnot m8.peptide -->
**H2** and **H4**: deferred (they need A3, M10-07/M10-10 and the full validation).
<!-- end m8.peptide -->

## 6. Pending in pass 2

Metric 6 (M10-10),
<!-- ifnot m8.peptide -->
metric 8 (M10-07),
<!-- end m8.peptide -->
TE2–TE5 and TE7 (M10-06), A3 (M10-08–11), benchmark (M10-12), report-freeze re-runs (M10-14).

## 7. Deviations from the spec

- DIA-NN 2.3.1 Academia instead of 2.7.0 (ADR-0051).
- Inputs and metric files are not committed (Khoa's rule of 2026-10-23).
- Protein groups are groups with at least one accepted precursor, with no protein-level FDR control (metric 1).
- Strain-specific genes are defined by peptide identity, not orthology (metric 3).
- Region tests are not compared in pass 1, and the Jaccards are not restricted to genes detected in
  both arms (H3; a conservative stand-in, biased against H3).

## 8. Notes

<!-- notes -->
