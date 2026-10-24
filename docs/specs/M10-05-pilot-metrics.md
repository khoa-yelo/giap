# M10-05: pilot metrics for A1/A2 (§15.5 metrics 1–5 and 7), TE1 on real data, `REPORT.md` draft

Milestone M10 · **pass 1** (MVP closer of the metrics side) · depends on: M10-04 (run directories
`out_a1`, `out_a2` with exit 0, the compare directory of `giap compare out_a1 out_a2 --labels
A1,A2`), M10-03 (the converted tables and `<a>.convert.json` under `$GIAP_ROOT/pilot/inputs/`),
M10-02 (`summary.json` per arm) · estimate ≈ 360 code lines + ≈ 450 test lines · **`stats`**
(light: Wilson intervals, one exact binomial test; see "Statistician review") · roles:
validation (`validation/pilot/metrics.py`, `validation/pilot/report_template.md`), test-author
(`tests/unit/test_m10_05_*.py`, `tests/integration/test_m10_05_*.py`), statistician (review) ·
Khoa: one question (Q-K1, what `REPORT.md` and the metric files may contain and whether they are
committed); the default does not wait for him.

## Spec sections implemented

- §15.5 metrics 1–5 and 7, for the arms A1 and A2 only. Metric 6 (novel ORFs, A2/A3 overlap) is
  M10-10's; metric 8 (free entrapment) is M10-07's.
- §11.6 TE1 on real data: "≥ 99.5 % of A2 peptides map exactly, in frame, to annotated CDS".
  Index question L5 (lead decision 1) fixes the denominator.
- §15.6 pass criteria: the checklist in the draft report, with the status pass 1 can give each
  bullet.
- §15.7 hypotheses 1 and 3, where A1/A2 data allow (hypotheses 2 and 4 wait for M10-07/M10-10
  and the full validation).
- Index decision 7: "REPORT.md is written from tables, never from screenshots or numbers typed by
  hand; each number cites the file it comes from."
- Khoa's 2026-10-23 rule (M10-02 lead decisions 4 and 6, M10-03 Q-K1): nothing derived from the
  DIA-NN reports is committed (report samples, column lists, values, converted tables). This spec
  reads it as covering every metric file and the rendered report until Khoa answers Q-K1.
- AGENTS.md conventions: pure functions for statistics; file I/O only at the edges; R3 (no
  threshold tuned on the pilot data), R4 (state what was not estimated), R11.

## What exists (main 3b2f1ffa; M10-02 `m10/diann-searches`, M10-03 `m10/diann-convert-tests`)

| Piece | State | Used for |
|---|---|---|
| `validation/proxies/te_proxies.py` | in-silico TE1/TE2/TE8. Pure: `te1_result(sources, excluded, peptide_hits, genes, circular)`, `non_atg_nterm(genes)`, `coding_genes`, `inframe_inside`, `ProxyResult`. Constant `TE1_MIN = 0.999` (the proxy's bar) | TE1 on real data **reuses `te1_result`** and `non_atg_nterm`; no copy of the in-frame logic |
| `giap.steps.mapping.select_queries(peptides, cfg)` | the mapping step's query set: non-decoy rows with `q_value ≤ input.q_value_max` (0.01) or NA, length in `mapping.exact.[min_len, max_len]` (7–40), one row per distinct `peptide` | the peptide universe of metric 2 and TE1, so the metrics count exactly what GIAP queried |
| `giap.context.summary.qc_rows` | writes `detection_summary.tsv` row `section=qc, metric=genome_inconsistent_rate` (share of distinct target, non-contaminant `search_db = other` peptides without an exact placement; `n` = that denominator) | metric 2 "genome-inconsistent peptides in A1" is **read**, not recomputed |
| `giap compare` (M9-22) | `compare_abundance.tsv` (ρ on shared genes: detected and not imputed in both; TU-block bootstrap CI; Bland–Altman), `compare_stability.tsv` (Jaccard per family `islands`, `operon_flags`, `enriched_categories`, `novel_orfs`, `alt_sites`; flags `both_empty`, `missing_in:<label>`), `compare_detection.tsv` (gene × run matrix) | metrics 4 and 5 are **read**; hypothesis 1 uses the matrix |
| `run_manifest.json` | `started_at`, `finished_at`, `exit_code`, `config` (resolved), `steps[]` with `name`, `status` (`completed`/`cached`/`failed`), `n_tasks`, `n_retries`, `wall_time_s`, `cpu_time_s`, `peak_rss_bytes` (from `logs/trace.tsv`, last attempt per task; absent metrics left out, never 0) | metric 7 |
| `genes.tsv` | `locus_tag`, `feature_type` {cds, sorf, pseudogene}, `contig`, `start`/`end` (1-based closed), `strand`, `crosses_origin`, `replicon_type` {chromosome, plasmid, unknown}, `region_ids` (`<type>_<n>`, §9.2), `detection_status` {detected, not_detected, single_peptide_unestimated, partial_excluded} | metric 3 |
| `peptide_hits.tsv`, `unmapped.tsv` | accepted placements (`route` {exact, tolerant}, `locus_class`, `gene`, coordinates 1-based closed); unmapped and contaminant sequences (`label`, `kind`, `reason`) | metric 2, TE1 |
| `giap.stats.cutoffs.wilson_ci(k, n, level)` | Wilson score interval, exact 0/1 endpoints, NaN for n = 0 | every share's CI |
| M10-03 `<a>.convert.json` | counts of the converted table incl. `distinct_peptides_target_q01`, `q_crossings` | cross-check of metric 1; `q_crossings` reported next to TE1 |
| M10-02 `summary.json` | `distinct_stripped_q01` on the **raw** report q-values | cross-check of metric 1 |
| M10-04 | spec not yet written (2026-10-23). This spec assumes `$GIAP_ROOT/pilot/runs/out_a1`, `…/out_a2`, `…/cmp_a1_a2` (question L1); `metrics.py` takes every path as an argument, so only the runbook depends on it | inputs |

## Inputs (all on FarmShare; none committed)

| Name | Path (default) | Metric |
|---|---|---|
| `OUT_A1`, `OUT_A2` | `$GIAP_ROOT/pilot/runs/out_a{1,2}/` | 2, 3, 7, TE1 |
| `CMP` | `$GIAP_ROOT/pilot/runs/cmp_a1_a2/` | 4, 5, H1, H3 |
| `PEP_A1`, `PEP_A2` | `$GIAP_ROOT/pilot/inputs/a{1,2}.tsv.gz` (M10-03) | 1, 2, TE1 |
| `CONV_A1`, `CONV_A2` | `$GIAP_ROOT/pilot/inputs/a{1,2}.convert.json` | 1 (cross-check), TE1 (`q_crossings`) |
| `SUM_A1`, `SUM_A2` | M10-02 `summary.json` of each arm (path as M10-02 lead decision 4 leaves it) | 1 (cross-check) |
| `BAKTA_JSON` | `tests/data/recorded/bakta/1.12.2-full/nissle/…json.gz` (committed) | TE1 (`non_atg_nterm`, gene intervals), strain-specific set |
| `GENOME` | `tests/data/genomes/nissle1917.fna.gz` (committed, circular-tagged); its sha256 must equal the runs' genome hash | TE1 (`read_coding_genes`, circularity) |
| `K12_FASTA` | `$GIAP_DATA_DIR/proteomes/UP000000625.fasta` | strain-specific set (metric 3) |
| `SACCT` (optional) | a TSV written by M10-04 from `sacct -P -o JobID,JobName,State,Elapsed,MaxRSS,NodeList,AllocCPUS` for each arm's head job | metric 7 (arm level) |

`metrics.py` checks before computing: both run manifests have `exit_code = 0`; the genome sha256
of both manifests is equal (as `giap compare` already demands); each `PEP_*` sha256 equals the
input sha256 the run's manifest recorded for `--peptides` (the metrics describe the table the run
used). Any failure exits 2 naming the file; nothing is written.

## Metric definitions

Notation: for an arm, **P** is its converted peptide table; a row is a *target* row when
`is_decoy` is false; *accepted* means `q_value ≤ 0.01` (inclusive; the value as written by the
converter, so the same set GIAP's preflight and mapping see). A `protein_ids` value is
*contaminant-only* when every `;` item starts with one of the run's `input.contaminant_prefixes`
(read from `run_manifest.json` `config`; default `Cont_`, `CON__`, `CONT_`, `cRAP`), the preflight
rule. Shares are given with a 95 % Wilson interval (`wilson_ci`); a share whose denominator is 0 is
`null` with `reason`, never 0.

### Metric 1: identifications at 1 % FDR per arm

From **P** only (accepted target rows):

- `precursors` = distinct (`modified_peptide`, `charge`);
- `peptides` = distinct `peptide` (stripped sequence);
- `protein_groups` = distinct `protein_ids` values after sorting the `;` items of each value (so
  `B;A` and `A;B` are one group).

Each count is given three ways: `all`, `contaminant_only`, `non_contaminant` (= all −
contaminant-only), and both over the union of runs and per run (`run` value). The FDR is DIA-NN's
global precursor q-value (`Global.Q.Value`, M10-03 F3); protein groups are **the groups of
accepted precursors, not groups at 1 % protein-group FDR**: the converter does not carry
`PG.Q.Value` (question L2). The label in the report says so.

Cross-checks (recorded, not failed): `peptides.all` must equal `CONV.counts.distinct_peptides_target_q01`
(else exit 1: the metric and the converter disagree on the same file); the difference to
`SUM.distinct_stripped_q01` is reported as `raw_vs_written_q01` and must equal the number of
distinct peptides among `CONV.counts.q_crossings` rows (rounded across the cut-off by `%.6g`).

### Metric 2: mapping shares; genome-inconsistent peptides in A1

Universe **U** = `select_queries(P, cfg)` (the run's cfg). **C** = the peptides of **U** whose
accepted target rows are all contaminant-only. **D** = **U** \ **C**. For each peptide of **D**,
from `peptide_hits.tsv`:

- `exact`: at least one row with `route = exact`;
- `tolerant`: no exact row, at least one `route = tolerant` row;
- `unmapped`: neither.

The three shares of |**D**| sum to 1. Reported per arm: |**U**|, |**C**|, |**D**|, the three
counts and shares with CIs, and the `unmapped` count split by `unmapped.tsv` `reason` for those
peptides (`no_hit`, `rejected_placement_fdr`, `below_identity`, `too_short`; a peptide with
several reasons counts under the first in that order; one absent from `unmapped.tsv` counts as
`not_listed`, expected 0). |**C**| is cross-checked against the distinct `unmapped.tsv` rows with
`label = contaminant, kind = peptide` (difference reported).

**Genome-inconsistent (A1):** `value` and `n` of the A1 `detection_summary.tsv` row
`section = qc, metric = genome_inconsistent_rate`, and the count `round(value × n)`. Under the
same universe this equals `tolerant + unmapped` of A1's **D**; the metric reports both and
`consistent = (difference = 0)` in counts. A mismatch is not a failure; the report explains it
(it would mean `qc_rows` and `select_queries` see different peptides). For A2 the row is NA
(`search_db = same`), which is expected.

### Metric 3: coverage of the strain CDS

Universe **G** = `genes.tsv` rows with `feature_type` in {cds, sorf} (the same locus tags in both
arms; `giap compare` already fails otherwise). Pseudogenes are reported as their own row, not in
**G**. *Detected* = `detection_status = detected`; the other statuses are counted per status, all
inside the denominator. Per arm and per subset: genes, detected, share and CI.

| Subset | Definition (per gene) |
|---|---|
| `all` | **G** |
| `chromosome`, `plasmid` | `replicon_type`; `plasmid` also split per `contig` (pMUT1, pMUT2 by their contig names) |
| `prophage` | some `region_ids` item starts with `prophage_` |
| `bgc` | some `region_ids` item starts with `bgc_` |
| `strain_specific` | the gene's protein has ≥ 1 in-silico tryptic peptide, and none of them occurs in the K-12 digest (below) |
| `k12_shared` | ≥ 1 in-silico peptide occurs in the K-12 digest |
| `no_peptides` | the protein yields no in-silico peptide (reported, in neither of the two above) |

**Strain-specific, defined by peptide identity, not orthology** (no homology tool is allowed, §16
"No MMseqs2"). Digest: `giap.mapping.digest.digest` with its defaults (trypsin_p, ≤ 1 missed
cleavage, 7–30 aa; the TE proxies' settings), I and L collapsed to one letter on both sides. K-12
side: every protein of `K12_FASTA`. Nissle side: the Bakta protein of each gene of **G**
(`coding_genes` of the parsed Bakta JSON, joined by locus tag). This is the set A1 cannot detect by
an exact match **by construction**, which is what hypothesis 1 asks about; the report states the
definition and that it is a proxy (question L3). The set and its size are a pure function of two
public inputs; its size, not its members, goes into `metrics.json`.

### Metric 4: quantitative agreement

The `compare_abundance.tsv` row with {`run_a`, `run_b`} = {A1, A2}: `n_shared`, `n_units`,
`spearman_rho`, `ci_low`, `ci_high`, `method`, `n_boot`, `ba_bias_log2`, `ba_loa_low`,
`ba_loa_high`, `flags`, copied unchanged. "Shared genes" is compare's definition (detected and not
imputed in both). `metrics.py` does not recompute ρ (one implementation; M9-22 owns it and its
review). Missing row or NA ρ: copied with its flags (`too_few_shared`, `missing_in:<label>`).

### Metric 5: stability of conclusions

The `compare_stability.tsv` rows for the A1/A2 pair and the families `islands`, `operon_flags`,
`enriched_categories`: `n_a`, `n_b`, `n_shared`, `jaccard`, `flags`, copied unchanged. The
`novel_orfs` and `alt_sites` rows are copied too, under `other_families` (metric 6 is M10-10's).
`enriched_categories` needs the function module (index L10, on); if it is off in a run, the row
is NA with `missing_in`, and the report says so.

### Metric 7: runtime and peak memory per arm and per step

Per arm, from `run_manifest.json`:

- per step: `name`, `status`, `n_tasks`, `n_retries`, `wall_time_s`, `cpu_time_s`,
  `peak_rss_bytes` (absent stays absent, written `null`);
- `wall_time_s` of the arm = `finished_at − started_at`;
- `sum_task_wall_s` = sum of the steps' `wall_time_s`; `max_step_peak_rss_bytes` = max of the
  steps' `peak_rss_bytes` (tasks run concurrently on separate allocations, so the max, not a sum,
  is the peak of one task; the report says it is not the arm's total memory);
- `cached_steps`: names of steps with `status = cached`. If any, the arm is
  `runtime_complete: false` and its totals are marked as lower bounds;
- from `SACCT` when given: the head job's `Elapsed`, `MaxRSS`, `NodeList`, `AllocCPUS`, `State`.
- `hosts`: the distinct `hostname` values of `logs/trace.tsv` (node class context for M10-12).

These are pilot measurements on a shared cluster, not the benchmark (M10-12 replaces the §13.2
estimates; §15.6 bullet 4 stays `pending` here).

### TE1 on real data (A2)

**Definition** (index L5, lead decision 1, unchanged): the denominator **D₂** is metric 2's **D**
for A2: distinct target peptides with `q_value ≤ 0.01` (`Global.Q.Value` as converted), decoys
out, contaminant-only peptides out, within the mapping length range. N-terminal products of strict
non-ATG genes stay **in** the denominator; they are counted (`non_atg_nterm`) and reported.

**Computation**, reusing the proxy:

1. `sources`: peptide → the set of `protein_ids` items of its accepted target rows that are not
   contaminant-prefixed (for A2 these are Bakta locus tags; mixed contaminant/strain peptides keep
   their strain items). Restricted to **D₂**.
2. `genes`, `circular`: from the proxy's own reader `read_coding_genes(BAKTA_JSON, GENOME)`
   (locus tag → the gene's 0-based half-open `Interval`; contig → length when circular, else
   None), exactly as `te1_proxy` builds them. No new coordinate code. `metrics.py` checks that
   these locus tags equal the cds/sorf/pseudogene locus tags of `out_a2/tables/genes.tsv` (exit 2
   otherwise: the run used another annotation).
3. `te1_result(sources, excluded={}, peptide_hits, genes, circular)` from
   `validation.proxies.te_proxies` (excluded is empty because of L5). Its `counts.exact_inframe`
   is the numerator; its `passed` flag (the proxy's 0.999) is **ignored**.
4. `non_atg_nterm(coding_genes(parse_bakta_json(BAKTA_JSON)))`: how many peptides of **D₂** are
   such products, and how many of those are misses.
5. Misses split by cause (first match): `non_atg_nterm`; `exact_other_gene` (an exact in-frame
   CDS/sORF row exists, but at a gene not in the peptide's `sources`); `exact_other_class` (exact
   rows only of other `locus_class`es, counted per class); `tolerant_only`; `unplaced` (per
   `unmapped.tsv` reason as in metric 2); `unknown_source` (no `sources` item is a locus tag of
   `genes.tsv`; expected 0).

**Judgement:** `share_exact_inframe = exact_inframe / |D₂|`; **pass iff share ≥ 0.995** (§11.6,
point estimate, as written; R3: not tuned). The Wilson 95 % CI is reported beside it, and the
report states when the CI's lower end is below 0.995 although the point estimate passes. Also
reported: `q_crossings` of A2 (M10-03) and the share recomputed with those peptides removed, so
a pass that depends on `%.6g` rounding at the cut-off is visible (question L4). `|D₂| = 0` gives
share `null`, `passed: false`, `reason: empty_denominator`.

TE1 is also computed for A1 as a **diagnostic only** (`te1_a1_diagnostic`, no pass flag): A1's
`sources` are UniProt accessions, not locus tags, so the "own gene" test does not apply; the
diagnostic reports the share of A1's **D** with an exact in-frame CDS/sORF row at any gene.

### Hypotheses 1 and 3 (reported, not gated)

**H1** "A2 detects more strain CDS than A1, mostly plasmid, prophage and BGC genes". From the
`compare_detection.tsv` matrix on **G**: b = genes detected in A2 only, c = in A1 only. Test:
exact one-sided binomial, p = P(X ≥ b | X ~ Bin(b + c, ½)) (`scipy.stats.binomtest`,
`alternative="greater"`); b + c = 0 gives p `null`. "Mostly": the share of the b genes in
plasmid ∪ prophage ∪ BGC, with its Wilson CI, and the same share among all of **G** for context;
also the share of the b genes that are strain-specific. **Verdict** `supported` iff b > c,
p ≤ 0.05 and the share > 0.5; `not_supported` otherwise; `not_evaluable` when b + c = 0. The
thresholds are fixed here, before M10-04's runs report numbers (question K2).

**H3** "context conclusions agree between A1 and A2 on shared genes": the metric 5 Jaccards and the
metric 4 ρ with CI, **without a verdict** (the spec fixes no agreement threshold; question K2).

## Outputs

All under `METRICS_DIR` (default `$GIAP_ROOT/pilot/metrics/`), never under the checkout:

| File | Content | Row-level data from DIA-NN? |
|---|---|---|
| `metrics.json` | `{"schema": "giap-pilot-metrics/1", "inputs": {name: {path, sha256}}, "spec": "<git SHA of this spec at compute time>", "m1": …, "m2": …, "m3": …, "m4": …, "m5": …, "m7": …, "te1": …, "te1_a1_diagnostic": …, "h1": …, "h3": …, "checklist": […], "warnings": […]}` | no: counts, shares, CIs, statistics, step timings |
| `te1.json` | the `te1` object alone, with `share_exact_inframe` (the field the index's evidence table names) | no |
| `te1_exceptions.tsv` | one row per miss: `kind`, `locus_tag`, `peptide`, `hits` (the proxy's exception rows, **without** the 50-row cap) | **yes** (peptide sequences): never committed under any Q-K1 answer |
| `REPORT.md` | rendered from `validation/pilot/report_template.md` and `metrics.json` | no (by the template rules below) |

JSON is written with sorted keys, two-space indent, floats via `repr` of Python floats, a final
newline, and **no timestamp**, so the same inputs give the same bytes (TE7 spirit). Lists are
sorted (steps in manifest order, which is deterministic; subsets in table order above).

## `REPORT.md` draft

### What it may contain (default; Q-K1 for Khoa)

**Default: aggregates only.** Counts, shares and their CIs, test statistics and p-values, ρ,
Jaccard, runtimes and memory, the §15.6/§15.7 statuses, file paths, job IDs, the GIAP commit, and
the sha256 of inputs that are already public or committed (genome, Bakta JSON, `strain_db.faa.gz`,
UP000000625). **Not:** any peptide or precursor sequence, protein or locus-tag list, per-gene
value, per-row table, column list of a DIA-NN report, or a sha256 of a DIA-NN-derived file (those
stay in `metrics.json` `inputs` on FarmShare). Named examples ("rplL detected in A2") are not in
the draft; M10-15 may add a few if Khoa allows (Q-K1 b).

**Default: nothing rendered is committed** until Khoa answers. `REPORT.md` is written to
`$GIAP_ROOT/pilot/metrics/REPORT.md`; `validation/pilot/REPORT.md` and `validation/pilot/metrics/`
do not exist in the repository (a guard test enforces it). This is the reversible choice: if Khoa
says the aggregates may be committed, one later commit copies `REPORT.md`, `metrics.json` and
`te1.json` and drops the guard line; a commit made too early cannot be taken back. The index's
evidence table (`validation/pilot/metrics/te1.json`) changes to `$GIAP_ROOT/pilot/metrics/te1.json`
under the default.

### Template (`validation/pilot/report_template.md`, committed; no numbers)

Placeholders `{{m1.A2.union.peptides.non_contaminant}}` (dotted paths into `metrics.json`); every
rendered number is followed by its source, e.g. `0.996 [0.990, 0.998] (metrics.json › te1)`. A
placeholder whose path is missing renders `**missing: <path>**` and the render exits 1 (no silent
gaps). Floats render with 3 significant digits for shares and ρ, `p` with 2 significant digits,
times as `h:mm:ss`, memory in GiB with 1 decimal. A hand-written notes file
`$GIAP_ROOT/pilot/metrics/notes.md` (prose by the validation agent: explanations of failures,
caveats) is included verbatim under "Notes" when present; the render rejects it if it contains a
digit sequence of 4 or more characters that is not in `metrics.json` (numbers come from tables,
index decision 7) or a string matching `[ACDEFGHIKLMNPQRSTVWY]{7,}` (a peptide).

Sections:

1. **Status**: the §15.6 checklist table (below), the GIAP commit, the spec SHA.
2. **Data and arms**: PXD078459 runs (from M10-01's committed manifest), DIA-NN 2.3.1 Academia
   (ADR-0051), the databases, job IDs of M10-02/03/04 (from `SACCT` and PROGRESS).
3. **Metric 1–5, 7** tables, one per metric, A1 and A2 side by side.
4. **TE1** with its miss breakdown and the `q_crossings` sensitivity line.
5. **Hypotheses**: H1 verdict and numbers; H3 numbers; H2 and H4 `deferred` with the reason.
6. **Pending in pass 2**: metric 6 (M10-10), metric 8 (M10-07), TE2–TE5 and TE7 (M10-06), A3
   (M10-08–11), benchmark (M10-12), report-freeze re-runs (M10-14).
7. **Deviations from the spec**: DIA-NN version (ADR-0051), inputs not committed (Khoa's rule),
   protein groups at precursor FDR (metric 1), strain-specific by peptide identity (metric 3).
8. **Notes** (optional, from `notes.md`).

### §15.6 checklist in pass 1

| Bullet | Status rule in `metrics.py` | Expected after M10-05 |
|---|---|---|
| GIAP runs end to end on A1, A2 and A3 | `pass` only when all three manifests exist with exit 0; with A1 and A2 at exit 0 and no A3: `partial (A1, A2; A3 pending M10-10)` | partial |
| ≥ 99.5 % of A2 peptides map exactly (TE1) | `pass` / `fail` from `te1.passed` | pass or fail |
| A1 and A2 detections and abundances are compared | `pass` when `compare_detection.tsv`, `compare_abundance.tsv` and `compare_stability.tsv` exist and validate against their schemas and the A1/A2 abundance row exists (ρ may be NA with a flag, which is reported) | pass |
| Measured runtimes replace the estimates | always `pending (M10-12)`; metric 7 is shown as pilot measurements | pending |
| TE2–TE5 are reported; failures explained | always `pending (M10-06)` | pending |

## Interface (`validation/pilot/metrics.py`)

Pure functions (DataFrames and dicts in, dicts out; no file access): `identifications(P, prefixes)`,
`mapping_shares(P, peptide_hits, unmapped, cfg)`, `genome_inconsistent(detection_summary)`,
`coverage(genes, strain_specific)`, `strain_specific_set(nissle_proteins, k12_proteins)`,
`abundance_agreement(compare_abundance, labels)`, `stability(compare_stability, labels)`,
`runtime(manifest, trace_rows, sacct_rows)`, `te1_real(P, peptide_hits, unmapped, genes,
circular, non_atg, cfg)` (genes and circular from `read_coding_genes`), `hypothesis1(compare_detection, coverage_subsets)`,
`checklist(metrics, run_states)`, `render(template, metrics, notes) -> str`. A thin CLI does the
reading and writing:

```text
python validation/pilot/metrics.py compute --a1 OUT_A1 --a2 OUT_A2 --cmp CMP \
    --peptides-a1 PEP_A1 --peptides-a2 PEP_A2 --convert-a1 CONV_A1 --convert-a2 CONV_A2 \
    [--summary-a1 SUM_A1 --summary-a2 SUM_A2] --bakta-json BAKTA_JSON --genome GENOME \
    --k12-fasta K12_FASTA \
    [--sacct SACCT] [--a3 OUT_A3] --out METRICS_DIR
python validation/pilot/metrics.py render --metrics METRICS_DIR/metrics.json \
    [--notes METRICS_DIR/notes.md] --out METRICS_DIR/REPORT.md
```

Exit codes: 0 written; 1 an internal inconsistency (metric 1 vs `convert.json`, a render
placeholder missing, a rejected notes file); 2 a precondition (missing file, run exit ≠ 0, genome
or peptide-table hash mismatch, `--out` inside the checkout). Run with the shared venv Python and
`PYTHONPATH=$REPO/src:$REPO` (for `validation.proxies.te_proxies`); no `uv run`. It reads only
GIAP's own outputs and committed or public inputs, no downloaded archive. Expected cost: one CPU,
< 2 min (the largest step is digesting ~9,000 proteins), so it runs inside an allocation, not on
a login node; if it exceeds 10 min it moves to a `dev` job (question L5).

## Statistician review

The PR carries the `stats` label. The statistician reviews, before the test step (methods first,
§14.1), the sections "Metric definitions" (CIs, the strain-specific proxy, the metric 2 and TE1
universes) and "Hypotheses 1 and 3" (the one-sided exact binomial test, the "mostly" rule, the
decision to give H3 no verdict); and in the PR: that `wilson_ci` and `binomtest` are called as
specified, that no threshold depends on the data (R3), and that the report wording does not read
an unsupported hypothesis as refuted or a CI as a test. ρ, its bootstrap CI and the Jaccards are
M9-22's (already reviewed) and are copied, not re-reviewed. No methods note is added: this section
and the template wording are the note (question L6).

## Examples

Hand-made values (not from the real reports).

| # | Input | Expected |
|---|---|---|
| 1 | **Metric 1.** A2 table: WT_1 `PEPTIDEK`/`PEPTIDEK`/z2/q 0.001/target/`JMJJOG_00010`; WT_1 `PEPTIDEK`/z3/q 0.004/target/`JMJJOG_00010`; WT_2 `PEPTIDEK`/z2/q 0.002/target/`JMJJOG_00010`; WT_2 `AAAAGGGK`/z2/q 0.02/target/`JMJJOG_00020`; WT_2 `KEDITPEP`/z2/q 0.005/**decoy**; WT_2 `LLLGGGR`/z2/q **0.01**/target/`Cont_P00761` | union: precursors all 3 (contaminant-only 1, non-contaminant 2); peptides all 2 (1, 1); protein groups all 2 (1, 1). Per run: WT_1 precursors 2, peptides 1; WT_2 precursors 2 (`PEPTIDEK`/2, `LLLGGGR`/2), peptides 2. `AAAAGGGK` (q > 0.01) and the decoy are not counted; q = 0.01 is |
| 2 | **Metric 1 groups.** two accepted rows with `protein_ids` `JMJJOG_00031;JMJJOG_00030` and `JMJJOG_00030;JMJJOG_00031` | `protein_groups.all` = 1 |
| 3 | **Metric 2.** A1 **D** of 10 peptides: 7 with an exact row, 2 with tolerant rows only, 1 in `unmapped.tsv` as `no_hit`; one more peptide with `protein_ids = Cont_P02769` | \|U\| 11, \|C\| 1, \|D\| 10; exact 0.7 [0.397, 0.892], tolerant 0.2, unmapped 0.1 (no_hit 1). `detection_summary` row `genome_inconsistent_rate` value 0.3, n 10 → count 3 = tolerant + unmapped, `consistent: true` |
| 4 | **Metric 3.** **G** = 100 genes (90 chromosome, 10 plasmid); A2 detects 60 (8 plasmid), A1 50 (2 plasmid) | A2 all 0.60 [0.502, 0.691]; A1 all 0.50 [0.404, 0.596]; plasmid A2 0.8 [0.490, 0.943], A1 0.2 [0.057, 0.510] |
| 5 | **Strain-specific.** Nissle gene X digests to {`AAAGGGLK`, `PEPTLDEK`}; K-12 digest contains `PEPTIDEK`; gene Y digests to {`QQQWWWR`} absent from K-12; gene Z is 5 aa long (no peptide of 7–30 aa) | X `k12_shared` (I/L collapsed, `PEPTLDEK` ≡ `PEPTIDEK`); Y `strain_specific`; Z `no_peptides` |
| 6 | **H1.** b = 12 A2-only, c = 2 A1-only; 7 of the 12 in plasmid ∪ prophage ∪ BGC | p = 106/16384 = 0.00647; share 0.583 [0.320, 0.807]; verdict `supported` |
| 7 | **Metrics 4, 5.** `compare_abundance` A1,A2: n_shared 45, ρ 0.83 [0.76, 0.88]; `compare_stability`: islands jaccard 0.5; operon_flags NA `both_empty`; enriched_categories NA `missing_in:A1` | values copied unchanged with their flags; report: islands 0.5; operon flags "none in either arm"; enriched categories "function module off in A1" |
| 8 | **Metric 7.** manifest `started_at` 10:00:00, `finished_at` 10:42:30; steps `mapping` completed 120 s, 2 GiB; `bakta` cached; `gene_table` completed 30 s, `peak_rss_bytes` absent | arm wall 2550 s; `sum_task_wall_s` 150; `max_step_peak_rss_bytes` 2147483648; `cached_steps` ["bakta"], `runtime_complete: false`; `gene_table` peak `null` |
| 9 | **TE1 pass.** \|D₂\| = 1,000: 996 exact in frame at an own gene (3 of them masked by `rna_overlap`), 2 `non_atg_nterm` misses, 1 exact only as `cds_outframe`, 1 `no_hit`; `q_crossings` 0 | share 0.996 [0.990, 0.998], `passed: true`, note "CI lower end below 0.995"; misses `non_atg_nterm` 2, `exact_other_class` {cds_outframe: 1}, `unplaced` {no_hit: 1}; `exact_inframe_masked` 3 |
| 10 | **TE1 fail and contaminants.** \|D₂\| = 1,000 with 994 exact in frame; a peptide `Cont_P00722;JMJJOG_00350` is in **D₂** with sources {`JMJJOG_00350`}; a `Cont_P00761`-only peptide is not | share 0.994 [0.987, 0.997], `passed: false`; checklist TE1 `fail`; the fail is listed for the notes |

## Edge cases

- **Empty input.** No accepted target peptide in an arm: metric 1 counts 0; metric 2 shares
  `null` (`reason: empty_denominator`); TE1 `passed: false` with that reason; checklist TE1
  `fail (no data)`. A run directory missing a table that the metric needs (for example no
  `compare_stability.tsv`): that metric is `null` with `reason: missing:<file>`, the others are
  computed, exit 0, and the checklist bullet that needs it is `fail`.
- **Circular and linear contigs.** The Nissle chromosome and plasmids are circular-tagged. A hit or
  gene across the origin: `te1_result` → `inframe_inside` already tries shifts ±L (unit-tested in
  M3-22); `circular` comes from `read_coding_genes` (the genome's topology tags). A linear contig
  gets `None` and no shift. Coverage and metric 2 read no coordinates. Test fixture: a gene with `crosses_origin =
  true` whose masked exact hit is reported unwrapped.
- **Minus strand.** The in-frame test for `-` genes uses the gene end (`(fe − b) % 3`), in the
  reused function. `peptide_hits.tsv` is 1-based closed; `te1_result` already converts
  (`int(s1) - 1`). `metrics.py` adds no coordinate conversion; a fixture with a minus-strand gene
  and a masked exact hit checks that the reused path is wired with the right inputs.
- **Ties and ordering.** Distinct-counts have no ties. Miss causes use a fixed first-match order.
  The unmapped-reason split uses a fixed order. H1's test is exact (no rank ties). ρ ties are
  compare's (average ranks). `te1_exceptions.tsv` is sorted by (`kind`, `locus_tag`, `peptide`).
- **q at the cut-off.** `q = 0.01` is accepted. Values above 0.01 rounded to `0.01` by the
  converter (`q_crossings`) are accepted because GIAP accepts them; TE1 reports the share without
  them as well.
- **Target rows with NA q-value.** Not expected (M10-03 F3). `select_queries` would keep them; the
  metrics count them in a warning `target_q_na` and, for metric 1, do not count them as accepted.
- **Peptide in several runs or charges.** One peptide in metric 2 and TE1; several precursors in
  metric 1.
- **I/L.** DIA-NN's `Stripped.Sequence` is literal; mapping is I/L-tolerant by its own rules. Only
  the strain-specific digest collapses I/L explicitly.
- **Selenocysteine and other recoded genes.** Peptides over a Sec site cannot map exactly in frame
  (TE2 proxy note, M7-16): they show up as misses (`exact_other_class` or `unplaced`) and the notes
  name the cause from the aggregate; no special case.
- **A1 `protein_ids` are UniProt accessions.** Metric 2 and the A1 TE1 diagnostic do not use them
  as genes; only TE1 on A2 checks own genes.
- **Run resumed (`-resume`).** Cached steps carry no fresh timing: `runtime_complete: false`.
- **Local executor (no `trace.tsv`).** `peak_rss_bytes` absent for every step, `hosts` empty;
  metric 7 still writes.
- **A3 present (`--a3`).** Only the checklist's first bullet reads it; no A3 metric is computed here.
- **`--out` inside the repository.** Exit 2 (the outputs hold DIA-NN-derived data).

## Acceptance criteria

On CI (no real data; every input is built in `tmp_path`):

1. Each pure function reproduces its example above exactly (examples 1–10), with CIs equal to
   `wilson_ci` to 1e-12 and the H1 p-value to 1e-12.
2. Metric 2's universe equals `select_queries` on the same table and cfg (a fixture with a decoy,
   a q = 0.02 row, a 41-aa peptide and a q-NA target shows each filtered or kept as
   `select_queries` does), and its three shares sum to 1.
3. TE1 calls `validation.proxies.te_proxies.te1_result` (asserted with a spy), passes
   `excluded={}`, and ignores the proxy's `passed`: a fixture with share 0.997 passes (proxy would
   fail at 0.999); 0.994 fails.
4. Coordinates: the minus-strand and origin-crossing fixtures (edge cases) give the masked
   exact hit as in frame; the same hit shifted by one base is a miss (`exact_other_class`).
5. `compute` on a hand-made pair of run directories plus a compare directory (minimal
   `genes.tsv`, `peptide_hits.tsv`, `unmapped.tsv`, `detection_summary.tsv`, `run_manifest.json`,
   `logs/trace.tsv`, `preflight/preflight.json`, the three compare tables, valid against their
   schemas) writes `metrics.json`, `te1.json`, `te1_exceptions.tsv`; a second call gives
   byte-identical files.
6. Preconditions: run exit ≠ 0, genome hashes differ, peptide-table hash ≠ manifest, a missing
   input, `--out` under the repository → exit 2 and nothing written. Metric 1 disagreeing with
   `convert.json` → exit 1.
7. `render`: every placeholder of the committed template resolves on the fixture's
   `metrics.json`; a removed key gives exit 1 naming the path; every number in the rendered file
   is followed by its source; `notes.md` with a 7-letter amino-acid string or an unknown 4-digit
   number is rejected (exit 1).
8. Checklist statuses follow the table above for: A1+A2 only (`partial`), A1+A2+A3 at exit 0
   (`pass`), TE1 fail, compare table missing (`fail`); bullets 4 and 5 are always `pending`.
9. Guard: `git ls-files` lists no `validation/pilot/REPORT.md`, nothing under
   `validation/pilot/metrics/`, and no `te1_exceptions.tsv`, `metrics.json` or `te1.json` anywhere
   (extends M10-03's guard).
10. Integration (marked `slow`; Slurm via `test_branch.sbatch`, not the session): two
    `giap run --executor local` runs on giapsim SIM100K (`tests/helpers/m4_real.py`): an "A2" run
    with the simulator's same-strain peptides plus two planted peptides (one from an out-of-frame
    translation, one non-ATG N-terminal product of a planted GTG-start gene), and an "A1" run of
    the same table with `search_db = other` and a few residues substituted so that tolerant and
    unmapped peptides exist; then `giap compare`; then `compute` and `render`. TE1 counts equal the
    planted truth (the two planted misses, everything else exact in frame); metric 2's A1
    tolerant + unmapped equals `genome_inconsistent_rate × n`; metric 3 subsets equal giapsim's
    planted prophage/BGC/plasmid gene sets; metric 7 has one entry per manifest step.

On FarmShare (evidence; PR description and PROGRESS, not CI):

11. `compute` and `render` on the real M10-04 outputs exit 0; `metrics.json` `warnings` are listed
    in the PR; `te1.json` `share_exact_inframe` and `passed` are reported in PROGRESS › Measured
    numbers (as an aggregate, per Q-K1); the rendered `REPORT.md` sits at
    `$GIAP_ROOT/pilot/metrics/REPORT.md`. A TE1 fail does not fail the issue: it is explained in
    `notes.md` (§15.6, "investigated and explained") and becomes a lead finding.

Link to §15.6: this issue gives evidence for bullet 2 (TE1) and the computed part of bullet 3, and
the checklist for all five.

## Test plan

Owner: the test-author, tests first (R2), new files only:

- `tests/unit/test_m10_05_metrics.py`: criteria 1–4 (pure functions, hand-made frames).
- `tests/unit/test_m10_05_cli.py`: criteria 5–8 (`tests/helpers/` gets a small builder for a
  synthetic run directory under a new name; schemas validated with the existing validator).
- `tests/unit/test_m10_05_guard.py`: criterion 9.
- `tests/integration/test_m10_05_sim_pilot.py`: criterion 10 (`slow`).

Local runs are targeted only (venv Python, `PYTHONPATH=$PWD/src:$PWD/tests/sim:$PWD`); the
integration test runs on Slurm.

## Estimate

| Piece | Code lines |
|---|---|
| metrics 1, 2 (+ genome-inconsistent) | 70 |
| metric 3 + strain-specific digest | 55 |
| metrics 4, 5 (readers of compare tables) | 25 |
| metric 7 | 35 |
| TE1 real (+ A1 diagnostic, miss causes) | 60 |
| H1, checklist | 35 |
| render + notes checks | 40 |
| CLI, preconditions, writers | 40 |
| **Total** | **≈ 360** (index: 360) + template (prose) |

Tests ≈ 450 lines. One PR.

## Questions

### For the lead (reversible defaults in use)

| # | Question | Default |
|---|---|---|
| L1 | M10-04 is not specced yet: run and compare directory paths? | `$GIAP_ROOT/pilot/runs/out_a1`, `out_a2`, `cmp_a1_a2`; M10-04 adopts them or this runbook follows M10-04 (arguments make it a one-line change) |
| L2 | Metric 1 "protein groups at 1 % FDR": the converter drops `PG.Q.Value`. Add it to the converted table (a converter change in M10-03), read the raw report (deleted at M10-14), or report groups of accepted precursors? | groups of accepted precursors, labelled so in the report |
| L3 | "Strain-specific genes": peptide identity against the K-12 digest (no shared in-silico tryptic peptide, I/L collapsed) as the proxy for "no K-12 counterpart"? | yes; the definition is printed in the report |
| L4 | TE1 when `q_crossings > 0`: the verdict uses the share as GIAP saw it (crossings included); the share without them is shown. Should the verdict use the stricter share? | verdict on the share as GIAP saw it; both shown |
| L5 | `metrics.py compute` in the session allocation (1 CPU, < 2 min), or always a `dev` job? | session allocation; `dev` job if > 10 min |
| L6 | Is this spec's statistics section enough as the methods note, or does the statistician want `docs/methods/pilot_metrics.md`? | the spec section; a note only if the statistician asks |
| L7 | The index evidence table names `validation/pilot/metrics/te1.json` (committed). Under Q-K1's default it is `$GIAP_ROOT/pilot/metrics/te1.json`. Update the index row now? | yes, in the implementation PR |

### For Khoa (next batch; the defaults do not block)

| # | Question | Default (reversible) |
|---|---|---|
| Q-K1 | (a) May `REPORT.md`, `metrics.json` and `te1.json` be committed when they hold aggregates only (counts, shares, CIs, ρ, Jaccard, p-values, runtimes), with no sequence, no ID list, no per-gene value and no hash of a DIA-NN-derived file? (b) May the final report name a few genes as examples (for example "the ribosomal proteins")? | (a) not committed until you answer; kept under `$GIAP_ROOT/pilot/metrics/`, review packet in `local/review/` (M10-15); content already restricted to aggregates. (b) no names |
| Q-K2 | §15.7 calls the hypotheses "fixed before running", but gives no decision rule. This spec fixes H1 (`supported` iff b > c, one-sided exact binomial p ≤ 0.05, and > 50 % of A2-only genes in plasmid ∪ prophage ∪ BGC) and gives H3 no verdict. Do you want a threshold for H3 (for example every Jaccard ≥ 0.5 and ρ's CI lower end ≥ 0.7) fixed before the M10-04 numbers are seen? | H1 as stated; H3 descriptive only. A rule set after the numbers are seen would not be pre-specified, so it is reported as post hoc |
