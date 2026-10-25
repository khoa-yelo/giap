# GIAP pilot report (arms A1, A2 and A3)

Rendered by `validation/pilot/metrics.py render --final` from `metrics.json` (docs/specs/M10-05-pilot-metrics.md, docs/specs/M10-15-pilot-gate.md).
Every value below is followed by the file and key it comes from; none is typed by hand.
This template holds no numbers. Only aggregates appear here: no sequence, no gene list, no per-gene value.

## 1. Status

| §15.6 bullet | Status | Evidence |
|---|---|---|
| GIAP runs end to end on A1, A2 and A3 | {{checklist.0.status}} | {{checklist.0.evidence}} |
| ≥ 99.5 % of A2 peptides map exactly (TE1) | {{checklist.1.status}} | {{checklist.1.evidence}} |
| A1 and A2 detections and abundances are compared | {{checklist.2.status}} | {{checklist.2.evidence}} |
| Measured runtimes replace the estimates | {{checklist.3.status}} | {{checklist.3.evidence}} |
| TE2–TE5 are reported; failures explained | {{checklist.4.status}} | {{checklist.4.evidence}} |

A failing bullet is explained in the notes below; no threshold was moved to make a bullet pass (R3).
GIAP commit: {{giap_commit}}. Spec SHA: {{spec}}.
Submitted for Khoa's review (§14.8 item 4).

## 2. Data and arms

PXD078459 wild-type runs (M10-01 manifest), DIA-NN 2.3.1 Academia (ADR-0051). A1 searched UniProt
K-12 (UP000000625) plus contaminants; A2 searched the Nissle strain database (`strain_db.faa`) plus
contaminants. A3 is de novo: Cascadia 0.0.7 on the same runs or, when M10-08's verdict rules it
out, Casanovo 5.2.1 on PXD053727 / UTI89 (fallback, a different strain).
A3 branch of the runbook: {{evidence.runbook.a3_branch}}; M10-08 verdict {{evidence.a3_feasibility.verdict}} (reason: {{evidence.a3_feasibility.reason}}).
Job IDs are in `reproduce.md` (Evidence table) and in the status table above.

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
<!-- m3 rows -->

Subsets: all of G (cds, sorf), by replicon (`plasmid:<contig>` per plasmid; `unknown`: neither
chromosome nor plasmid), by region (prophage, bgc), by K-12 class (`no_peptides`: no in-silico
tryptic peptide); `pseudogene` is outside G.

### Metric 4: quantitative agreement (shared genes, from `giap compare`)

Spearman rho {{m4.spearman_rho}} [{{m4.ci_low}}, {{m4.ci_high}}] on {{m4.n_shared}} shared genes (TU-block bootstrap); Bland–Altman bias {{m4.ba_bias_log2}} log2; flags: {{m4.flags}}.

### Metric 5: stability of conclusions (Jaccard, from `giap compare`)

| family | Jaccard | note |
|---|---|---|
| islands | {{m5.islands.jaccard}} | {{m5.islands.note}} |
| operon flags | {{m5.operon_flags.jaccard}} | {{m5.operon_flags.note}} |
| enriched categories | {{m5.enriched_categories.jaccard}} | {{m5.enriched_categories.note}} |

### Metric 6: novel ORFs and variants, A2/A3 overlap (M10-10)

Descriptive, unestimated: counts and Jaccards carry no error estimate.

| family | Jaccard A2/A3 |
|---|---|
| novel ORFs (called) | {{m6.overlap_a2_a3.novel_orfs.called.jaccard}} |
| ALT sites (high confidence in both arms) | {{m6.overlap_a2_a3.alt_sites.high.jaccard}} |

<!-- if m8.universe -->
### Metric 8: A1's free entrapment (from `entrapment.py`, M10-07)

The entrapment is free: it is the K-12-only part of A1's own database, the K-12 peptides absent from
Nissle's translation after the I/L and mass-equivalent rules. Estimator: Wen et al. (combined
method; reference in the methods note). r = K-12-only / shared in-silico peptides = {{m8.universe.r}} ({{m8.universe.k12_only}} K-12-only, {{m8.universe.shared}} shared); at precursor level
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

The 95 % bound assumes {{m8.precursor.ci.assumes}}.
Warning codes of `entrapment.py` (details in `entrapment.json`; a digest mismatch is a lead finding): {{m8.warning_codes}}.

<!-- end m8.universe -->
<!-- ifnot m8.universe -->
### Metric 8: free entrapment FDP on A1 (M10-07)

| level | combined FDP | lower bound | assessment |
|---|---|---|---|
| peptide (primary) | {{m8.peptide.fdp_combined}} | {{m8.peptide.fdp_lower}} | {{m8.peptide.assessment}} |
| precursor | {{m8.precursor.fdp_combined}} | {{m8.precursor.fdp_lower}} | {{m8.precursor.assessment}} |

<!-- end m8.universe -->
## 4. TE1 on real data (A2)

Every A2 target identification, true or false, is a substring of an annotated protein, so TE1 checks
that database, genome and mapping agree; it says nothing about identification error.
The expected misses are structural (N-terminal products of non-ATG starts, Sec sites); anything else is a mapping defect.

Share exact in frame at an own gene: {{te1.share_exact_inframe}} of {{te1.n_peptides}} peptides, nominal 95 % Wilson interval [{{te1.ci_low}}, {{te1.ci_high}}] (misses cluster by gene, so its coverage is optimistic; the interval is not a test).
Pass (point estimate at or above the threshold): {{te1.passed}}; lower end below the threshold: {{te1.ci_low_below_threshold}}.

Misses ({{te1.n_misses}}, first matching cause): non-ATG N-terminal {{te1.misses.non_atg_nterm}}, exact at another gene {{te1.misses.exact_other_gene}}, exact in another class {{te1.misses.exact_other_class}}, tolerant only {{te1.misses.tolerant_only}}, unplaced {{te1.misses.unplaced}}; masked exact in frame {{te1.exact_inframe_masked}}.
Non-ATG N-terminal products in the denominator: {{te1.non_atg_nterm.peptides}}. Peptides with no locus tag of `genes.tsv` among their sources: {{te1.unknown_source}}; with only NA q-values: {{te1.target_q_na}}.

Secondary, never the verdict: excluding non-ATG products {{te1.share_excluding_non_atg}}; without peptides whose accepted rows are all at the cut-off {{te1.share_without_q_at_cutoff}}.
Removed {{te1.share_without_q_at_cutoff.removed}} peptides (A2 q_crossings {{te1.share_without_q_at_cutoff.q_crossings}}); a pass that becomes a fail under this superset: {{te1.share_without_q_at_cutoff.pass_becomes_fail}}.
A1 diagnostic (exact in frame at any gene, no pass flag): {{te1_a1_diagnostic.share}}.

## 5. True examples (TE2–TE5) and determinism (TE7)

From `te_real.json` and `te7.json` (docs/specs/M10-06-te-real-and-te7.md); every definition was fixed
before the real data were read (R3). Shares carry nominal 95 % Wilson intervals, descriptive only.
Every TE that did not pass has a `## Diagnosis: TE<n>` section in the notes below, naming its cause class.

| TE | Status | Headline |
|---|---|---|
| TE2: no unexplained `high` ALT call in A2 | {{te_real.te2.status}} | violations {{te_real.te2.counts.violations}} |
| TE3: ribosomal proteins detected (≥ 90 %) and in the top quartile | {{te_real.te3.status}} | detected {{te_real.te3.detection.share}}; median percentile rank {{te_real.te3.abundance.percentile_rank}} |
| TE4: operon stoichiometry (iBAQ, molar scale) | {{te_real.te4.status}} | log2 RplL/RplJ {{te_real.te4.pairs.rplL/rplJ.log2_ratio}} in [1, 3]; log2 AtpA/AtpD {{te_real.te4.pairs.atpA/atpD.log2_ratio}} in [−1, 1] |
| TE5: trypsin and keratin peptides in the contaminant bucket (A1 and A2) | {{te_real.te5.status}} | in the bucket: A1 {{te_real.te5.A1.share_bucket}}, A2 {{te_real.te5.A2.share_bucket}} |
| TE7: a second A2 run gives byte-identical tables | {{te7.giap.status}} | gate files compared {{te7.giap.n_gate}} |

TE7 is decided by GIAP's own tables only; the DIA-NN layer (a second search of the same spectra) is an
upstream finding, never the verdict: {{te7.diann.status}}. Raw data to GIAP tables reproduce byte for byte: {{te7.end_to_end}}.
A `known_lapack` difference is still a fail (M4-20 decision 9; writer rounding is the decision 10 backlog item).
TE7 is not a §15.6 bullet: a TE7 fail is shown here and in section 9 with its diagnosis and listed as an open item.

## 6. Hypotheses

Test family for this issue: {H1}, size 1, no BH.

**H1** (A2 detects more strain CDS than A1, mostly plasmid, prophage and BGC genes): verdict {{h1.verdict}}; failed clauses: {{h1.failed_clauses}}.
`not_supported` means not supported by the pre-registered rule; `not_evaluable` means the rule could not be applied (no discordant TU).
A2-only genes b {{h1.b}}, A1-only genes c {{h1.c}}; A2-only share of discordant genes {{h1.prop_a2_only}}; coverage difference {{h1.coverage_diff}}.
These are nominal 95 % Wilson intervals; they describe the effects and are not the test (a
two-sided 95 % interval and the one-sided 0.05 test can disagree near the boundary).
TU-level test: {{h1.b_tu}} against {{h1.c_tu}} units ({{h1.mixed_tu}} mixed, left out), p_tu {{h1.p_tu}}, unit effect {{h1.unit_effect}}.
Gene-level p_gene {{h1.p_gene}} (assumes independent genes; decides nothing).
Share of A2-only genes in plasmid, prophage or BGC {{h1.mostly}} (base rate {{h1.mostly_base_rate}}); strain-specific share of A2-only genes {{h1.strain_specific_share}}.
The arms are two database searches of the same spectra, so discordance is largely fixed by sequence
divergence rather than sampled, and dependence at the scale of a region remains; p is not a replication probability.

**H3** (context conclusions agree between A1 and A2): verdict {{h3.verdict}}; rho criterion {{h3.criteria.rho.status}}; islands {{h3.criteria.islands.status}}; operon flags {{h3.criteria.operon_flags.status}}; enriched categories {{h3.criteria.enriched_categories.status}}.
The numbers are shown whatever the verdict: rho {{h3.criteria.rho.rho}} [{{h3.criteria.rho.ci_low}}, {{h3.criteria.rho.ci_high}}].

| family | n_a | n_b | n_shared | union | Jaccard | status |
|---|---|---|---|---|---|---|
| islands | {{h3.criteria.islands.n_a}} | {{h3.criteria.islands.n_b}} | {{h3.criteria.islands.n_shared}} | {{h3.criteria.islands.union}} | {{h3.criteria.islands.jaccard}} | {{h3.criteria.islands.status}} |
| operon flags | {{h3.criteria.operon_flags.n_a}} | {{h3.criteria.operon_flags.n_b}} | {{h3.criteria.operon_flags.n_shared}} | {{h3.criteria.operon_flags.union}} | {{h3.criteria.operon_flags.jaccard}} | {{h3.criteria.operon_flags.status}} |
| enriched categories | {{h3.criteria.enriched_categories.n_a}} | {{h3.criteria.enriched_categories.n_b}} | {{h3.criteria.enriched_categories.n_shared}} | {{h3.criteria.enriched_categories.union}} | {{h3.criteria.enriched_categories.jaccard}} | {{h3.criteria.enriched_categories.status}} |

**H2** (the strain-specific database lowers genome-inconsistent identifications): deferred (full validation, FDRBench); shown side by side only: A1 genome-inconsistent share {{m2.A1.genome_inconsistent.value}} and the free entrapment estimate (peptide level) {{m8.peptide.fdp_combined}}.

**H4** (genome-anchored de novo precision against the precision measured on A2): verdict {{h4.verdict}}.
Only a lower bound on precision exists (sequence-level confirmation by A2), so the data can count against H4 but never confirm it; it is `unestimated` unless A3 is the Cascadia branch.

## 7. Not evaluated

One line per metric or hypothesis whose value is null with a reason: the pilot data do not allow it (§15.7 "reported where the pilot data allow").

<!-- not evaluated -->

## 8. Runtimes

### Metric 7: runtime and peak memory (pilot measurements on a shared cluster, not the benchmark)

| | A1 | A2 |
|---|---|---|
| wall time | {{m7.A1.wall_time_s}} | {{m7.A2.wall_time_s}} |
| sum of step wall times | {{m7.A1.sum_task_wall_s}} | {{m7.A2.sum_task_wall_s}} |
| max step peak RSS (one task, not the arm's total memory) | {{m7.A1.max_step_peak_rss_bytes}} | {{m7.A2.max_step_peak_rss_bytes}} |
| complete (no cached step; else totals are lower bounds) | {{m7.A1.runtime_complete}} | {{m7.A2.runtime_complete}} |

### Benchmark (M10-12, `benchmarks/baseline.json`, wheat nodes)

Medians of the standard run against the estimates of §13.2 (35 min, range 16–55, and 25 min); a median at an estimate is reported, not judged.

| cell | median wall time | commit | CPU |
|---|---|---|---|
| K-12, standard | {{evidence.benchmark.cells.k12_standard.median_s}} | {{evidence.benchmark.cells.k12_standard.commit}} | {{evidence.benchmark.cells.k12_standard.cpu}} |
| K-12, core | {{evidence.benchmark.cells.k12_core.median_s}} | {{evidence.benchmark.cells.k12_core.commit}} | {{evidence.benchmark.cells.k12_core.cpu}} |
| simulated, standard | {{evidence.benchmark.cells.sim5m_standard.median_s}} | {{evidence.benchmark.cells.sim5m_standard.commit}} | {{evidence.benchmark.cells.sim5m_standard.cpu}} |
| simulated, core | {{evidence.benchmark.cells.sim5m_core.median_s}} | {{evidence.benchmark.cells.sim5m_core.commit}} | {{evidence.benchmark.cells.sim5m_core.cpu}} |

Benchmark definition (ADR-M10-D): {{evidence.benchmark.adr}}, status {{evidence.benchmark.adr_status}}.

## 9. Reproducibility

Freeze commit of the runbook: {{evidence.runbook.freeze_commit}} (runbook status {{evidence.runbook.status}}); every arm, the compare and the metrics ran at it (`validation/pilot/reproduce.md`, Evidence table).
The runbook's check writes `verify-<commit>.json` beside the freeze runs; it is part of the review packet.
TE7 at the freeze commit: {{te7.giap.status}}.
The benchmark was not re-run at the freeze commit (M10-12 L12); its commit is in section 8.
Row-level outputs stay on FarmShare and are removed only by the runbook's guarded clean-up step.

## 10. Deviations from the spec

- DIA-NN 2.3.1 Academia instead of 2.7.0 (ADR-0051).
- Row-level inputs and outputs are not committed (Khoa's rule of 2026-10-23); aggregates and hashes are (2026-10-24).
- Protein groups are groups with at least one accepted precursor, with no protein-level FDR control (metric 1).
- Strain-specific genes are defined by peptide identity, not orthology (metric 3).
- Region tests are not compared, and the Jaccards are not restricted to genes detected in
  both arms (H3; a conservative stand-in, biased against H3).
- H2 is deferred to the full validation; metric 8 stands beside it without a verdict.
- The A3 branch when it is the fallback or absent, every deferral made at the freeze, and any
  accepted ADR that changes a pilot setting are stated in the notes.

## Notes

<!-- notes -->
