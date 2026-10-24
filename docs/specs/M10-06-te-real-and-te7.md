# M10-06: real-data TE2–TE5 (`validation/pilot/te_real.py`) and TE7 (a second A2 run, byte comparison)

Milestone M10 · **pass 2** · depends on: M10-05 (`metrics.py`, `metrics.json`, the report template
and its notes rules), M10-04 (the accepted A2 attempt under `$GIAP_ROOT/runs/pilot/pass1/`,
`ATTEMPTS.tsv`), M10-02 and M10-03 (`diann_pilot.sbatch`, `diann_convert.sbatch`, both honouring
`OUT_DIR`) · estimate ≈ 310 code lines (index: 280) + ≈ 380 test lines · **not `stats`** (no test
statistic, no p-value; one reused Wilson interval, descriptive) · roles: validation
(`validation/pilot/te_real.py`, the template section, `notes.md` diagnoses), test-author
(`tests/unit/test_m10_06_*.py`, `tests/integration/test_m10_06_*.py`), lead (the three Slurm jobs,
R14 budget, PROGRESS) · Khoa: nothing new; Q-K1 of M10-05 (what may be committed) still governs.

## Spec sections implemented

- §11.6 TE2 "no `high` ALT calls from A2 except at sites flagged `start_codon` or `stop_codon`"
  (pilot).
- §11.6 TE3 "ribosomal proteins ≥ 90 % detected in A2, with median abundance in the top quartile"
  (pilot).
- §11.6 TE4 "stoichiometry within operons: log2(RplL/RplJ) in [1, 3] (about 4×, from the rplKAJL
  operon); log2(AtpA/AtpD) in [−1, 1] (α3β3)" (pilot, iBAQ).
- §11.6 TE5 "trypsin and keratin peptides land in the contaminant bucket" (pilot).
- §11.6 TE7 "two runs with the same seed give byte-identical tables" and §11.7 "Any
  non-deterministic output fails CI", applied to a full real-data A2 run.
- §15.6 bullet 5: "TE2–TE5 are reported. A failure is investigated and explained in the report
  before Khoa's review." The index's evidence row: "`te_real.json` has one entry per TE with
  `pass`/`fail` and, for each fail, a `REPORT.md` section with the diagnosis".
- §7.3.4 (iBAQ; "with a relative scale … stoichiometry checks are skipped"), §7.1.2
  (ribosomal-protein rule), M3-15 (the unmapped/contaminant bucket rules), M7-14 (`te2_result`,
  `id_mixups`), M4-20 decision 9 (LAPACK results are not bit-portable across CPU types).
- R3 (no threshold or definition is tuned on the pilot data: every definition below is fixed
  before `te_real.py` runs on real data), R4 (state what was not evaluated), R9, R11, R14, and
  Khoa's 2026-10-23 rule (nothing derived from the DIA-NN reports is committed).

## What exists (main 1f567901; M10-02 to M10-05 on their branches)

| Piece | State | Used for |
|---|---|---|
| `validation/proxies/te_proxies.py` | `te2_result(sites) -> ProxyResult` (counts `alt_sites`, `high`, `high_start_codon`, `high_stop_codon`, `violations`; exemption set `TE2_EXEMPT`; up to 50 example violations). No TE3/TE4 proxy; the TE5 proxy is a test (`tests/integration/test_te5_proxy.py`, M3-15) | TE2 **reuses `te2_result`**; no copy of its rule |
| `giap.annotate.bakta_layers` | `RIBOSOMAL_PRODUCT = ^(30S|50S) ribosomal protein`, `RIBOSOMAL_GENE = ^(rps|rpl|rpm)[A-Z]`, `is_ribosomal_protein(g)` (§7.1.2) | TE3's gene set (see the gap below) |
| `giap.steps.mapping.select_queries` | the queried peptide set (M10-05 metric 2's **U**) | TE5 universe |
| `genes.tsv` | `locus_tag`, `gene`, `product`, `feature_type`, `detection_status`, `abundance_log2`, `abundance_imputed`, `abundance_type`, `abundance_scale`, `top3_log2`, `n_peptides_unique`, `n_observable_peptides`, `region_ids`, `replicon_type` | TE2 breakdowns, TE3, TE4 |
| `variants.tsv` | `variant_id`, `gene`, `aa_pos`, `ref_aa`, `alt_aa`, `alt_type`, `call`, `confidence` {high, candidate}, `n_alt`, `class_error_rate`, `flags` | TE2 |
| `unmapped.tsv` | `sequence_id`, `kind`, `label` {contaminant, unmapped}, `reason`, `contaminant_match` | TE5 |
| `detection_summary.tsv` | QC row `id_mixups` (ALT sites of any confidence with same-database support and neither flag; M7-14) | TE2 cross-check |
| `$GIAP_DB_DIR/contaminants/contaminants.fasta` | 381 records, `>sp|Cont_<ACC>|<NAME> <description> OS=…`; the file both DIA-NN arms searched (M10-03) | TE5 families |
| `diann_pilot.sbatch A2` (M10-02) | writes everything to `${OUT_DIR:-$GIAP_ROOT/pilot}/A2/`, including `outputs.sha256` (one line per output file) and `provenance.json` (node, image sha256, cfg sha256). The first A2 search: Slurm 1879956, rye-02 (`CPU_SKU:9754`), 16 CPUs, 14 min 21 s | the second search (TE7, DIA-NN layer) |
| `diann_convert.sbatch A2` (M10-03) | reads `${OUT_DIR}/A2/A2_report.parquet`, writes `${OUT_DIR}/inputs/a2.tsv.gz` + `a2.convert.json`; conversion is byte-deterministic (its own determinism step) | converts the second search |
| `giap_run.sbatch` | head job; runs the main checkout; `run_manifest.json` has `git_commit`, `git_dirty`, `resume`, `inputs[]` sha256, `seeds`, `config`, `databases`, `steps[]`; `logs/trace.tsv` has `hostname` per task | the second GIAP run |
| `tests/helpers/byte_identity.py` | `tree_bytes(root)`, `differing(a, b)` | the TE7 tree comparison reuses the same idea (in `te_real.py`, which may not import `tests/`) |
| `scripts/measure_runtimes.py` | host → node-class lookup from the trace `hostname` | TE7 diagnosis (which CPU class ran each step) |
| M10-05 `metrics.py` | `compute`, `render`, the §15.6 checklist (bullet 5 "always `pending (M10-06)`"), `notes.md` rules (no unknown 4-digit number, no 7-letter amino-acid string) | the hook below |

### Gap found while writing this spec: §7.1.2's ribosomal regex catches modification enzymes

On the committed Nissle annotation (`tests/data/recorded/bakta/1.12.2-full/nissle/`, public data,
read 2026-10-23), `is_ribosomal_protein` selects **59** CDS. Six of them are enzymes that act on
ribosomal proteins, not ribosomal proteins: `prmA` "50S ribosomal protein L11 methyltransferase",
`prmB`, `rimL` "50S ribosomal protein L7/L12-serine acetyltransferase", `ycaO`, `rimK`, `rimO`.
They are low-abundance and would pull TE3's detected share and median down for a reason that is
not about the pipeline. TE3 therefore uses a stricter rule (below, question L1), which selects
**53** CDS on that annotation (including the zinc-limitation paralog `rpmE2` "type B 50S ribosomal
protein L31"). The same annotation names no `rplQ`, `rplT` or `rpsH` gene and names two genes
`rpsD`; TE3 takes the annotation as it is and the report states the set size. The §7.1.2 rule
itself is not changed (it is GIAP's CAI reference, out of scope); a backlog item records it.

## Inputs (all on FarmShare; none committed)

| Name | Path (default) | Used by |
|---|---|---|
| `OUT_A2` | the accepted pass-1 A2 attempt, resolved through `$GIAP_ROOT/runs/pilot/pass1/ATTEMPTS.tsv` (M10-05 lead decision 1) | TE2–TE5, TE7 reference |
| `OUT_A1` | the accepted pass-1 A1 attempt (same resolution) | TE5 (A1 half) |
| `PEP_A1`, `PEP_A2` | `$GIAP_ROOT/pilot/inputs/a{1,2}.tsv.gz` | TE5 (accepted rows and `protein_ids`) |
| `CONTAM` | `$GIAP_DB_DIR/contaminants/contaminants.fasta` | TE5 families |
| `OUT_REP` | `$GIAP_ROOT/runs/pilot/te7/out_a2_rep1` (and `out_a2_rep2` when needed, see TE7) | TE7 |
| `DIANN_A`, `DIANN_B` | `$GIAP_ROOT/pilot/A2/` and `$GIAP_ROOT/pilot/te7/A2/` (search directories) | TE7, DIA-NN layer |
| `CONV_B` | `$GIAP_ROOT/pilot/te7/inputs/a2.tsv.gz` and `a2.convert.json` | TE7, DIA-NN layer |

Preconditions of `te_real.py compute` (exit 2, nothing written): each run's manifest has
`exit_code = 0`; the `--peptides` sha256 in each manifest equals the sha256 of `PEP_*`; the
contaminant FASTA's sha256 equals the one the manifests record for the contaminants database
(when the manifest records it; otherwise a `warnings` entry `contaminants_unpinned`).

## Definitions

Notation as in M10-05: *accepted* = `q_value ≤ 0.01`; *target* = `is_decoy` false; **G** = the
`genes.tsv` rows with `feature_type` in {cds, sorf}; *measured* = `abundance_log2` not NA and
`abundance_imputed` false. Every TE has `status` ∈ {`pass`, `fail`, `not_evaluable`} and, when not
`pass`, a `reason` list. Shares carry a 95 % Wilson interval (`giap.stats.cutoffs.wilson_ci`),
labelled nominal and descriptive (as TE1: a census share, judged on the point value; the CI is
never the test).

### TE2: no unexplained `high` ALT call in A2

1. `te2_result(OUT_A2/tables/variants.tsv)` from `validation.proxies.te_proxies`, unchanged:
   counts, `passed`, the first 50 violations. The proxy's `notes` are dropped (they describe the
   in-silico proxy, not the pilot).
2. **status** `pass` iff `violations = 0`; `fail` otherwise. A missing `variants.tsv` (variants
   module off or failed) is `not_evaluable`, `reason: missing:tables/variants.tsv`. A file with
   zero rows is `pass` with `vacuous: true` (reported as such).
3. **Diagnosis aggregates** for the violations (all of them, not only 50): counts by `alt_type`, by
   `call`, by the sorted set of other flags present (for example `mass_equivalent`), by `n_alt`
   bin (1, 2, ≥ 3), by replicon type and by region type of the gene (`region_ids` prefixes
   `prophage_`, `bgc_`, …; join on `variants.gene` = `genes.locus_tag`), and the count at
   `aa_pos = 1` without `start_codon` (a start-codon rule miss would show here).
4. **Cross-check:** `violations ≤ id_mixups` (the QC row counts every confidence); otherwise a
   `warnings` entry `te2_id_mixups_inconsistent`. `id_mixups` NA: entry `id_mixups_na`.

### TE3: ribosomal proteins detected and abundant in A2

**Set R** (question L1): the genes of **G** whose `product` matches

```text
^(type [AB] )?(30S|50S) ribosomal protein [LS][0-9]+(/L[0-9]+)?( [A-Z][^ ]*)?$
```

(the protein name, optionally one capitalised trailing token such as `RpsL`; anything after it,
such as "methyltransferase", excludes the gene), **or** whose `gene` matches §7.1.2's
`RIBOSOMAL_GENE` (`^(rps|rpl|rpm)[A-Z]`, reused, not copied). Pseudogenes are out.

- **Clause 1 (detection):** `share = |R ∩ detected| / |R|`, detected = `detection_status =
  detected` (as M10-05 metric 3). Holds iff `share ≥ 0.90`.
- **Clause 2 (abundance):** **Q** = the measured genes of **G** (R's measured genes included).
  `q3 = numpy.quantile(Q.abundance_log2, 0.75)` (default linear method), `m =
  numpy.median(abundance_log2 of R ∩ measured)`. Holds iff `m ≥ q3`. Reported with `m`, `q3`,
  `|R ∩ measured|`, `|Q|` and the percentile rank of `m` in **Q** (share of **Q** with value ≤ `m`).
  Rank-based, so any `abundance_type` other than `none` is evaluable (the type is reported).
- **status:** `pass` iff both clauses hold; `fail` iff an evaluable clause fails;
  `not_evaluable` when `|R| = 0` (`reason: no_ribosomal_genes`), or clause 1 holds and clause 2
  is not evaluable (`|R ∩ measured| = 0` or `abundance_type = none`).
- **Pre-registered secondary numbers** (never the verdict): clause 1 and 2 on `R` without the
  `type B` paralogs (expressed under zinc limitation); clause 1 and 2 on the §7.1.2 set
  (`is_ribosomal_protein`), so the effect of the rule choice is visible.
- **Per-gene list** for the investigation (`te_real_exceptions.tsv`, kind `te3_not_detected`):
  locus tag, gene, product, `detection_status`, `n_observable_peptides`.

### TE4: stoichiometry within two operons (A2)

Pairs (fixed by §11.6): `(rplL, rplJ, [1, 3])` and `(atpA, atpD, [−1, 1])`, closed intervals.

- **Lookup:** the **G** row whose `gene` equals the name exactly (case-sensitive). No row:
  `not_evaluable`, `gene_missing:<name>`; more than one: `gene_ambiguous:<name>`. (On the committed
  Nissle annotation each of the four names occurs once.)
- **Requirements:** both genes measured; `abundance_scale = molar` (§7.3.4: stoichiometry checks
  are skipped on a relative scale; `reason: relative_scale`). The `abundance_type` is reported;
  `ibaq` is the expected default.
- **Value:** `log2_ratio = abundance_log2(a) − abundance_log2(b)` (float64, no rounding before the
  comparison). A pair holds iff `low ≤ log2_ratio ≤ high`.
- **status:** `pass` iff both pairs are evaluable and hold; `fail` iff an evaluable pair fails;
  `not_evaluable` otherwise (each pair's status reported).
- **Pre-registered secondary:** the same ratio on `top3_log2` (when both are present), and the
  `n_peptides_unique` and `n_observable_peptides` of the four genes, for the diagnosis.

### TE5: trypsin and keratin peptides land in the contaminant bucket (A2 and A1)

**Families**, from `CONTAM` headers (description = the header text after the ID token, up to
` OS=`):

| Family | Rule on the description | In the verdict |
|---|---|---|
| `trypsin` | `^(Trypsin|Anionic trypsin|Cationic trypsin)$` (porcine `Cont_P00761`, bovine `Cont_Q29463`; trypsin inhibitors and chymotrypsinogens are out) | yes |
| `keratin` | `^Keratin, type (I|II) ` (86 records) | yes |
| `keratin_associated` | `^Keratin-associated protein` (91 records) | no: reported (question L5) |

**Set S** (per arm): the peptides of **U** = `select_queries(P, cfg)` with at least one accepted
target row whose `protein_ids` contains an ID of a verdict family. For each peptide of S, in this
first-match order:

1. `genome_exact`: an accepted `route = exact` row in `peptide_hits.tsv` (the sequence is in the
   genome; M3-15 keeps accepted placements out of the bucket). **Exempt**, counted.
2. `bucket`: an `unmapped.tsv` row with `kind = peptide`, `sequence_id` = the peptide and `label =
   contaminant`. **Pass.** Sub-count `match_other_entry` when `contaminant_match` is not an ID of
   the peptide's own family (legal: M3-15 takes the lowest library ID that contains it).
3. Violations: `placed_tolerant` (an accepted tolerant placement and no exact one: the peptide did
   not land in the bucket), `labelled_unmapped` (in `unmapped.tsv` with `label = unmapped`),
   `absent` (neither placed nor listed; expected 0).

Per arm: `|S|`, per family, the counts above, `share_bucket = bucket / (|S| − genome_exact)` with
its Wilson CI. Arm status: `pass` iff `|S| − genome_exact > 0` and no violation; `fail` iff a
violation; `not_evaluable` iff `|S| − genome_exact = 0` (`reason: no_trypsin_keratin_peptides`).
**TE5 status:** `pass` iff both arms pass; `fail` iff either fails; `not_evaluable` otherwise.
The `keratin_associated` family gets the same counts under `secondary`.

### TE7: a second A2 run, byte for byte

Two layers. **Only the GIAP layer decides TE7**; the DIA-NN layer is an upstream reproducibility
finding, reported beside it and never folded into the verdict.

#### GIAP layer (the verdict)

**The pair.** The reference is the accepted pass-1 A2 attempt when it is *eligible*: its manifest
has `resume: false`, no step `cached`, `git_dirty: false`, and its `git_commit` is
*code-equivalent* to the replicate's commit, i.e. `git diff --quiet REF REP -- src/ workflow/
pyproject.toml uv.lock scripts/slurm/giap_run.sbatch` succeeds. Otherwise two fresh replicates
(`out_a2_rep1`, `out_a2_rep2`) on one commit form the pair (question L8). Each replicate is a
**fresh** run: a new outdir (so a new Nextflow work directory), no `--resume`, the same command as
M10-04's A2 (same `--peptides` path `$GIAP_ROOT/pilot/inputs/a2.tsv.gz`, genome, annotation,
`--db-dir`, `--threads 14 --memory 54.GB`, default seed).

**Pair preconditions** (exit 2, nothing written, the message names the field): both manifests
`exit_code = 0`; every `inputs[]` sha256 equal; `seeds` equal; `databases` names, versions and
checksums equal; resolved `config` equal after removing the keys that hold the outdir; eligibility
as above.

**File classes** (paths relative to OUT, matched in this order):

| Class | Paths | Treatment |
|---|---|---|
| `excluded` | `work/**`, `logs/**`, `.nextflow*`, `run_manifest.json`, `report/**`, `preflight/**` | not compared (timestamps, absolute paths, host names); listed by count |
| `external` | `annotation/tools/**` (native outputs of Bakta, UniOP, Promotech, geNOMAD, antiSMASH, …) | compared and reported; **not** in the verdict (tools outside GIAP, R5) |
| `figures` | `plots/**` except `plots/data/**` | compared and reported; not in the verdict (TE7 names tables) |
| `gate` | every other regular file: `tables/**`, `mapping/**`, `annotation/*` (top level), `plots/data/**`, `schemas/datapackage.json` | **byte comparison decides TE7** |

**Verdict:** `pass` iff the `gate` file sets are equal and every gate file has the same sha256 in
both runs; `fail` otherwise. There is no tolerance (R3; §11.7).

**Diagnosis of each differing gate file** (first match; descriptive only, never a tolerance):

1. `only_in:<run>`;
2. `gzip_header_only`: `.gz` files whose decompressed bytes are equal;
3. `row_order_only`: the header is equal and the sorted lists of data lines are equal;
4. `numeric_only`: same header and row count, and every differing cell parses as a float in both
   runs. Reported: the column names, the number of differing cells, and the maximum absolute and
   relative difference. Sub-class `known_lapack` when every differing column is in the committed
   list of LAPACK-dependent columns (`validation/pilot/te7_known_columns.tsv`, holding the
   `regions_stats.tsv` float columns of M4-20 decision 9 and the join-count columns of M5-11
   decisions 44, 46 and 47) **and** the producing step ran on different CPU classes in the two runs;
5. `content`: anything else (the first differing line number only, no content).

The producing step of a file is the registered step that declares it as an output (the step
registry); its CPU class in each run comes from the trace `hostname` of that step's tasks through
`measure_runtimes.py`'s host-to-class lookup. When no step declares the file, `producer: unknown`
and the diagnosis lists the CPU classes per step of both runs instead. A `known_lapack`-only result is still `fail`; the report says it is the known
non-portability (M4-20 decision 9) and links the M4-20 decision 10 backlog item (writer rounding).
Whether it should count as a pass is question L9.

#### DIA-NN layer (reported, never the verdict)

A second A2 search with the M10-02 chain, unchanged except for the output directory and the node
class: `OUT_DIR=$GIAP_ROOT/pilot/te7 sbatch -J giap-diann-A2-te7 --constraint=CPU_SKU:9754
scripts/slurm/diann_pilot.sbatch A2` (both steps, library prediction included, 16 CPUs, the class
of the first search: question L10), then its check, then `OUT_DIR=$GIAP_ROOT/pilot/te7 sbatch -J
giap-diann-convert-A2-te7 scripts/slurm/diann_convert.sbatch A2`. Compared, in order:

| Level | Compared | Result |
|---|---|---|
| L0 | every file name of the two `outputs.sha256` lists, equal or not (names and booleans only; logs and manifests are expected to differ by time stamps and paths) | `files_equal`, `files_differ` |
| L1 | sha256 of `A2_report.parquet` | equal → `identical_report` |
| L2 | sha256 of the converted `a2.tsv.gz` (the converter reads only its mapped columns and is byte-deterministic) | equal → `identical_converted` |
| L3 | when L2 differs: aggregate difference of the two converted tables, keyed by (`run`, `modified_peptide`, `charge`, `is_decoy`): rows only in A, only in B, in both; Jaccard of the accepted target precursors and of the accepted target peptides; rows in both whose `q_value` differs, whose acceptance flips, whose `intensity` differs, and the maximum absolute log2 intensity ratio | `nondeterministic` with these aggregates |

`diann.status` is the first of `identical_report`, `identical_converted`, `nondeterministic`.
**End-to-end claim:** when `diann.status` ≠ `nondeterministic` and the GIAP layer passes, the
report states that raw data → GIAP tables reproduce byte for byte (the replicate's converted table
has the same bytes as the input the GIAP pair used, so the GIAP pair covers it). When DIA-NN is
nondeterministic, the GIAP pair still runs on the original `a2.tsv.gz` (TE7 needs identical
inputs), the report shows the L3 aggregates as an upstream finding, and no GIAP run on the second
table is made in this issue (question L11). Same node class, same image, same cfg, same thread
count: a difference is DIA-NN's own (threading), not the hardware's; the report says which.

## Investigation of failures (§15.6 "investigated and explained")

For each TE2–TE5 with `status` ≠ `pass`, and for TE7 when it fails, the validation agent writes a
section `## Diagnosis: TE<n>` in `$GIAP_ROOT/pilot/metrics/notes.md`, built from `te_real.json`,
`te7.json` and `te_real_exceptions.tsv` (aggregates only in the text). Each section names one
**cause class**, the evidence, and the action:

| Cause class | Meaning | Required action |
|---|---|---|
| `pipeline_bug` | a GIAP defect (mapping, variant flags, bucket rule, gene table, determinism) | a backlog entry or issue linked in the section; the fix in its own PR; the TE re-evaluated on a re-run (M10-04 resume rules), both results kept |
| `annotation` | the Bakta model (a name, a start, a missing gene) explains it | stated; no change to the annotation |
| `search` | the DIA-NN identification or settings explain it (for example no N-terminal acetylation searched) | stated; a settings change is an ADR-0051 amendment for Khoa, not done here |
| `biology` | the sample explains it (growth phase, a paralog not expressed) | stated with the evidence |
| `definition` | the TE's definition does not fit the data | stated as a deviation for Khoa; **no** new definition or threshold (R3); only the pre-registered secondary numbers above may be cited |

`render` (M10-05) gains one rule: a TE with `status` ≠ `pass` and no `## Diagnosis: TE<n>`
section in `notes.md` is exit 1, naming the TE. The section must name one cause class from the
table (exit 1 otherwise).

## Outputs

All under `METRICS_DIR` (default `$GIAP_ROOT/pilot/metrics/`), never under the checkout (exit 2
for `--out` inside the repository, as M10-05):

| File | Content | Row-level DIA-NN-derived data? |
|---|---|---|
| `te_real.json` | `{"schema": "giap-pilot-te-real/1", "inputs": {name: {path, sha256}}, "spec": <git SHA of this spec>, "te2": …, "te3": …, "te4": …, "te5": {"A1": …, "A2": …, "status": …}, "warnings": […]}` | no: counts, shares, CIs, four log2 values and ratios |
| `te7.json` | `{"schema": "giap-pilot-te7/1", "pair": {run: {outdir, git_commit, job}}, "giap": {"status", "n_gate", "n_external", "n_figures", "n_excluded", "differing": [{path, class, sub_class, columns, n_cells, max_abs, max_rel, producer, cpu_classes}], "external_differing": […], "figures_differing": […]}, "diann": {"status", "l0": {...}, "l3": {...} or null, "nodes": […]}, "end_to_end": bool}` | no values; GIAP table **paths** and column names, and DIA-NN file names |
| `te_real_exceptions.tsv` | `te`, `kind`, `id` (locus tag or peptide), `detail`; sorted by (`te`, `kind`, `id`) | **yes**: never committed |
| `te7_files.tsv` | `path`, `class`, `sha256_a`, `sha256_b`, `equal` | hashes of derived files: never committed |

JSON as M10-05: sorted keys, two-space indent, final newline, no timestamp; the same inputs give
the same bytes.

**`metrics.py` hook** (≈ 15 lines in M10-05's file): `compute` takes optional `--te-real FILE`
and `--te7 FILE` and copies them under `te_real` and `te7`. Checklist bullet 5 becomes: `pending
(M10-06)` without `--te-real`; `pass` when TE2–TE5 all pass; `reported (k not passed, explained)`
otherwise, with the explanation enforced by `render`. TE7 is shown in its own report section
(it is not a §15.6 bullet). The template (`validation/pilot/report_template.md`) gets a section
"True examples (TE2–TE5) and determinism (TE7)" with one table row per TE (status, the headline
number, its source) and drops TE2–TE5 and TE7 from "Pending in pass 2".

**Committed:** `te_real.py`, `te7_known_columns.tsv` (column names of GIAP's own tables, from the
two cited decisions), the template section, the `metrics.py` hook, tests, and a README section
"TE2–TE5 and TE7 (M10-06)" with job IDs, states, nodes, elapsed times and exit codes only.
PROGRESS gets the job IDs and, as aggregates per Q-K1's default, only the five statuses.

## Interface (`validation/pilot/te_real.py`)

Pure functions: `te2(variants, genes, detection_summary)`, `ribosomal_set(genes)`, `te3(genes)`,
`te4(genes)`, `contaminant_families(headers)`, `te5_arm(P, peptide_hits, unmapped, families, cfg)`,
`classify_tree(files_a, files_b, known_columns, producers) -> dict` (file maps of path → bytes
reader), `diann_layer(sums_a, sums_b, conv_a, conv_b) -> dict` (frames for L3). Thin CLI:

```text
python validation/pilot/te_real.py compute --a2 OUT_A2 --a1 OUT_A1 \
    --peptides-a1 PEP_A1 --peptides-a2 PEP_A2 --contaminants CONTAM --out METRICS_DIR
python validation/pilot/te_real.py te7 --ref OUT_A --rep OUT_B --repo GIAP_REPO \
    [--diann-a DIANN_A --diann-b DIANN_B --peptides-a PEP_A2 --peptides-b CONV_B] --out METRICS_DIR
```

Exit codes: 0 written, whatever the TE statuses (a fail is a finding, not a crash); 1 an internal
inconsistency (a cross-check that must hold by construction, such as `|S|` per family exceeding
`|S|`); 2 a precondition. Run with the shared venv Python and `PYTHONPATH=$REPO/src:$REPO`; no
`uv run`. Cost: `compute` < 1 min on 1 CPU; `te7` hashes two run trees (a few GB with the
`external` class) and, when L2 differs, reads two converted tables: a `dev` job (1 CPU, 8 GB,
30 min limit), not the session.

## Runbook (lead; R14 and order)

1. After M10-05 merged and the accepted pass-1 A2 attempt exists: `te_real.py compute` (session
   or `dev`), then `metrics.py compute --te-real …` and `render`.
2. The second DIA-NN search (16 CPUs, ~15 min) and its conversion (4 CPUs, minutes), chained with
   `--dependency=afterok`. Submit only when the sum S of running and pending `giap-*` CPUs is
   ≤ 112.
3. The GIAP replicate: `sbatch -J giap-run-pilot-a2-te7 scripts/slurm/giap_run.sbatch …
   --outdir $GIAP_ROOT/runs/pilot/te7/out_a2_rep1` with M10-04's A2 arguments, when S ≤ 14 (one
   arm needs up to 114 CPUs) and no other pilot arm runs (index decision 3). It does not wait for
   step 2: its input is the original table. A second replicate, if eligibility fails, follows the
   same rule, one at a time.
4. `te_real.py te7` in a `dev` job; then `metrics.py compute --te7 …` and `render`; diagnoses
   into `notes.md` for every non-pass; render again.
5. README section and PROGRESS (job IDs at submission, statuses at the end).

## Examples

Hand-made values (not from the real reports).

| # | Input | Expected |
|---|---|---|
| 1 | **TE2 pass.** `variants.tsv`: `high` M at a GTG start, flags `start_codon;database_derived`; `high` U→* readthrough, flags `stop_codon`; `candidate` N→D, flags `database_derived`. `id_mixups` = 1 | counts `alt_sites` 3, `high` 2, `high_start_codon` 1, `high_stop_codon` 1, `violations` 0; `status: pass`; cross-check 0 ≤ 1, no warning |
| 2 | **TE2 fail.** One extra `high` N→D at `aa_pos` 57 of a minus-strand plasmid gene in `prophage_2`, flags `database_derived;mass_equivalent`, `n_alt` 2; `id_mixups` = 2 | `violations` 1, `status: fail`; breakdown `alt_type` {substitution: 1}, flags {`database_derived;mass_equivalent`: 1}, `n_alt` {2: 1}, replicon {plasmid: 1}, region {prophage: 1}, `aa_pos_1_unflagged` 0; one `te2_violation` row in `te_real_exceptions.tsv`; checklist bullet 5 `reported (1 not passed, explained)` only once `notes.md` has `## Diagnosis: TE2` |
| 3 | **TE3 pass, then clause 2 fails.** **G** has 20 measured genes with `abundance_log2` 1, 2, …, 20; R = 10 genes, 9 detected; the 9 detected r-proteins are measured at 12 … 20 | clause 1: 9/10 = 0.90 [0.596, 0.982], holds (≥ 0.90 inclusive); `q3` = 15.25, `m` = 16.0, percentile rank 0.80: holds; `status: pass`. With the r-proteins at 4 … 12 instead (others at 1–3 and 13–20): `m` = 8.0 < 15.25, rank 0.40, `status: fail`, `reason: [abundance]` |
| 4 | **TE3 detection boundary and secondaries.** \|R\| = 52, 46 detected; the 6 missed include `rpmE2` (type B) | share 0.885 [0.770, 0.946]: clause 1 fails, `status: fail`. Secondary without type B: 46/51 = 0.902 [0.790, 0.957] (reported, verdict unchanged). The §7.1.2-set secondary is computed on the genes `is_ribosomal_protein` selects in the same fixture |
| 5 | **TE4.** `rplL` 22.0, `rplJ` 20.1 → 1.9 (float64 1.8999999999999986) in [1, 3]; `atpA` 19.0, `atpD` 19.6 → −0.6 in [−1, 1]; scale molar, type ibaq | both pairs hold, `status: pass`; Top3 secondaries reported. Variants: `rplL` 23.0, `rplJ` 20.0 → 3.0 holds; 23.25 → 3.25 fails (`status: fail`); `rplJ` imputed → that pair `not_evaluable`, `status: not_evaluable`; `abundance_scale = relative` → both `relative_scale`, `not_evaluable` |
| 6 | **TE5 (A2).** Accepted target peptides of **U**: `VATVSLPR` (`Cont_P00761`) in `unmapped.tsv` as `contaminant`, `no_hit`, match `Cont_P00761`; `SLNNQFASFIDK` (`Cont_P04264`) in the bucket with match `Cont_P13645`; `LSSPATLNSR` (`Cont_P00761;JMJJOG_01234`) with an accepted exact hit; a `Cont_P02769` (albumin) peptide; `LASYLDK` (`Cont_P13645`) with an accepted tolerant placement only | S = 4 (albumin is not a family); `genome_exact` 1 (exempt); `bucket` 2 (`match_other_entry` 1); `placed_tolerant` 1; `share_bucket` 2/3; arm `fail`; TE5 `fail`. Without the tolerant peptide: 2/2 = 1.0 [0.342, 1.0], arm `pass` |
| 7 | **TE7 pass.** Reference and replicate at the same commit, both fresh; `run_manifest.json`, `logs/trace.tsv` and `report/report.html` differ; every `tables/`, `mapping/`, `annotation/*`, `plots/data/` file and `schemas/datapackage.json` has the same sha256; one `annotation/tools/promotech/…` file differs | `giap.status: pass`; `external_differing` lists the Promotech file (not in the verdict); excluded counted |
| 8 | **TE7 fail, classified.** `tables/genes.tsv` differs in 3 `abundance_log2` cells by ≤ 2e-15; `tables/pileup.tsv.gz` differs, decompressed equal; `tables/islands.tsv` has the same lines in another order; `tables/regions_stats.tsv` differs only in its M4-20 float columns, the producing step ran on `CPU_SKU:6148` vs `CPU_SKU:9754`; `mapping/extra.tsv` only in the replicate | `giap.status: fail`; classes: genes `numeric_only` (columns [abundance_log2], 3 cells, max_abs 2e-15), pileup `gzip_header_only`, islands `row_order_only`, regions_stats `numeric_only`/`known_lapack`, extra `only_in:rep1`. The report lists five diagnoses; each non-LAPACK class is a `pipeline_bug` candidate |
| 9 | **DIA-NN identical after conversion.** `A2_report.parquet` sha256 differs (the file embeds the output path); `a2.tsv.gz` sha256 equal; GIAP layer passes | `diann.status: identical_converted`; `l3: null`; `end_to_end: true`; the report says raw data → GIAP tables reproduce byte for byte, and why the report file itself differs (L0 names) |
| 10 | **DIA-NN nondeterministic.** Converted tables differ: 12 rows only in A, 9 only in B, 3 acceptance flips, max \|log2 intensity ratio\| 0.004; GIAP pair (original table twice) identical | `diann.status: nondeterministic` with those aggregates and the precursor and peptide Jaccards; `giap.status: pass`; `end_to_end: false`; TE7 is reported as passed; DIA-NN's threading nondeterminism is a separate finding in the report (same node class, so not hardware) |

## Edge cases

- **Empty input.** No `high` call: TE2 `pass`, `vacuous: true` when `variants.tsv` has no row at
  all. R empty: TE3 `not_evaluable`. No family peptide: TE5 arm `not_evaluable`. All of these need
  a diagnosis section (they are not `pass`), except a vacuous TE2 pass, which the report labels.
- **Circular and linear contigs.** TE2–TE5 read no coordinates: TE2 judges flags (the variant
  step already handled origin-crossing codons), TE3/TE4 join on locus tags and names, TE5 on
  sequences. A fixture puts a TE2 violation on a gene with `crosses_origin = true` on a circular
  plasmid and another on a linear contig: both counted once, replicon breakdown correct. TE7
  compares bytes, whatever the topology.
- **Minus strand.** Same: no strand-dependent logic is added. Example 2's violation sits on a
  minus-strand gene; `aa_pos` is taken as written (protein coordinates, strand-independent).
- **Ties and boundaries.** TE3: `share = 0.90` holds; `m = q3` holds; an even number of measured
  r-proteins takes NumPy's mean of the two middle values; duplicate abundance values need no rule
  (quantile and median are order statistics). TE4: values exactly on 1, 3, −1 or 1 hold (closed).
  TE5: a peptide in both families counts once in S and once in each family row; first-match order
  decides its class. TE7: files are compared in sorted path order; diagnosis classes are
  first-match. Exceptions are sorted by (`te`, `kind`, `id`).
- **Duplicate names.** Two `rpsD` rows in the Nissle model are both in R (two genes). A TE4 name
  found twice is `gene_ambiguous`, never resolved by position.
- **A1 in TE5.** A1's `protein_ids` are UniProt accessions plus the same `Cont_` contaminant IDs;
  the family rule is the same. A K-12 peptide that happens to be keratin-identical cannot occur
  (S needs a family ID).
- **TE7 reference resumed or on another commit.** Ineligible: two replicates on one commit. A
  replicate that fails (exit ≠ 0) is an M10-04-style attempt; R9 applies (three attempts, then
  Blocked).
- **TE7 with a module switched off** in one run only: the pair precondition on `config` fails
  (exit 2). With different task hosts but no differing file: `pass`; hosts are recorded.
- **`.gz` gate files.** A differing gzip header is a fail (GIAP writes `mtime = 0`, so it is a
  determinism bug), classified `gzip_header_only`.
- **DIA-NN replicate fails or lands on another class.** A failed search is retried once (R9 counts
  it), then the DIA-NN layer is `not_run` with the reason, and TE7 is still decided by the GIAP
  pair. A different node (the constraint cannot be met) is recorded in `diann.nodes`, and the
  report says the hardware may contribute.
- **Contaminant FASTA changed since the searches.** Its sha256 differs from the manifest's: exit 2
  (TE5 would judge families against another library).

## Acceptance criteria

On CI (no real data; every input is synthetic or hand-made in `tmp_path`; giapsim for the
integration test):

1. `te2` calls `validation.proxies.te_proxies.te2_result` (spy) and reproduces examples 1–2,
   including every breakdown and the `id_mixups` cross-check warning when violations exceed it;
   a missing `variants.tsv` gives `not_evaluable`; an empty one `pass` with `vacuous: true`.
2. `ribosomal_set` on a fixture with the products of the gap section (L11 methyltransferase,
   L7/L12-serine acetyltransferase, S6--L-glutamate ligase, `S12 RpsL`, `type B 50S ribosomal
   protein L31`, a gene named `rpmH` with a non-standard product, a pseudogene `rplA`) selects
   exactly the r-proteins and none of the enzymes or the pseudogene; on the committed Nissle
   annotation it selects 53 genes, and `is_ribosomal_protein` 59 (both asserted; public input).
3. `te3` reproduces examples 3–4 (shares and Wilson CIs to 1e-12, `q3`, `m`, ranks) and the
   boundary cases (`share = 0.90`, `m = q3` hold); the secondaries never change `status`.
4. `te4` reproduces example 5 and its variants, including `gene_missing`, `gene_ambiguous`,
   `relative_scale` and the closed bounds.
5. `contaminant_families` on the real header strings of the committed `tests/data/contaminants_mini.faa`
   plus hand-written headers (trypsin inhibitor, chymotrypsinogen, keratin-associated protein)
   assigns exactly the families above; `te5_arm` reproduces example 6, and a peptide outside
   `select_queries` (41 aa, or a decoy) is not in S.
6. `classify_tree` reproduces examples 7–8 class by class on two hand-made trees; `known_lapack`
   needs both conditions (columns listed and CPU classes different); a `known_lapack`-only pair is
   still `fail`.
7. `diann_layer` reproduces examples 9–10 on hand-made `outputs.sha256` lists and converted
   tables; the L3 aggregates match an independent count in the test.
8. Preconditions: input sha256, seeds, databases or config differing; a reference with `resume:
   true` or a cached step; non-code-equivalent commits (a temporary git repository with a change
   under `src/`); `--out` in the repository → exit 2, nothing written. A change only under `docs/`
   keeps the pair eligible.
9. `metrics.py` hook: bullet 5 is `pending (M10-06)` without `--te-real`, `pass` with four passes,
   and `reported (k not passed, explained)` otherwise; `render` exits 1 when a non-pass TE has no
   `## Diagnosis: TE<n>` section or the section names no cause class.
10. Determinism of the script: `compute` and `te7` twice on the same inputs give byte-identical
    `te_real.json`, `te7.json` and TSVs.
11. Guard (extends M10-05's): `git ls-files` lists no `te_real.json`, `te7.json`,
    `te_real_exceptions.tsv` or `te7_files.tsv`.
12. Integration (`slow`, Slurm via `test_branch.sbatch`): two `giap run --executor local` runs of
    giapsim SIM100K with the same seed in two directories, with planted trypsin and keratin
    peptides from `contaminants_mini.faa` added to the peptide table; `te7` gives `pass`; TE5 gives
    `pass` with S equal to the planted family peptides not in the genome (oracle
    `tests/helpers/contaminant_oracle.py`); TE2 runs on the run's `variants.tsv`. A copy of one run
    with one `genes.tsv` cell changed gives `fail` with `numeric_only` on that column.

On FarmShare (evidence; PR description, README section and PROGRESS, not CI):

13. The second DIA-NN search and its conversion exit 0, and the GIAP replicate passes M10-04's
    rule set A; `te_real.py compute` and `te7` exit 0; `metrics.py compute --te-real --te7` and
    `render` exit 0, which, by criterion 9, means every non-pass TE has its diagnosis.
14. The PR lists the five statuses (TE2–TE5, TE7), `diann.status` and the cause class of every
    diagnosis; numbers stay under `$GIAP_ROOT/pilot/metrics/` (Q-K1). A failing TE does not fail
    the issue; an unexplained one does. Every `pipeline_bug` diagnosis has a linked backlog entry
    or issue.

Link to §15.6: bullet 5 gets its evidence here; TE7 is §11.6/§11.7 evidence on real data.

## Test plan

Owner: the test-author, tests first (R2), new files only:

- `tests/unit/test_m10_06_te_real.py`: criteria 1–5 (pure functions, hand-made frames; the
  committed Nissle annotation for criterion 2).
- `tests/unit/test_m10_06_te7.py`: criteria 6–8, 10 (hand-made trees, a temporary git repository).
- `tests/unit/test_m10_06_hook.py`: criteria 9, 11 (needs M10-05's `metrics.py`).
- `tests/integration/test_m10_06_sim_te7.py`: criterion 12 (`slow`, Slurm).

Local runs are targeted only (venv Python, `PYTHONPATH=$PWD/src:$PWD/tests/sim:$PWD`, `-n 2` at
most); the integration test runs on Slurm.

## Estimate

| Piece | Code lines |
|---|---|
| TE2 (reuse + breakdowns + cross-check) | 30 |
| TE3 (set, two clauses, secondaries) | 45 |
| TE4 | 30 |
| TE5 (families, classes, two arms) | 50 |
| TE7 GIAP layer (pair checks, classes, diagnosis) | 65 |
| TE7 DIA-NN layer (L0–L3) | 35 |
| CLI, preconditions, writers | 40 |
| `metrics.py` hook and the render rule | 15 |
| **Total** | **≈ 310** (index: 280) + `te7_known_columns.tsv`, template section (prose) |

Tests ≈ 380 lines. One PR. If review pushes it over ~400 lines, TE7 splits off as M10-06b (its
own CLI subcommand already).

## Questions for the lead (reversible defaults in use)

| # | Question | Default |
|---|---|---|
| L1 | TE3's "ribosomal proteins": the §7.1.2 rule picks up six modification enzymes on Nissle (gap section). Use the stricter product rule (53 genes) for the verdict, with the §7.1.2 set as a secondary? | yes; backlog item for §7.1.2's CAI reference (not changed here) |
| L2 | TE3 "median abundance in the top quartile": median `abundance_log2` of the measured r-proteins ≥ the 75th percentile (NumPy linear) of all measured strain CDS/sORF, r-proteins included? | yes |
| L3 | TE3: keep the zinc-limitation paralog (`rpmE2`, "type B") in R, with a secondary without it? | yes |
| L4 | TE4 on `abundance_log2` (iBAQ by default), requiring a molar scale and both genes measured (not imputed); Top3 as a secondary only? | yes |
| L5 | TE5 families: trypsin = description exactly "Trypsin"/"Anionic trypsin"/"Cationic trypsin"; keratin = "Keratin, type I/II …"; keratin-associated proteins reported but not judged? | yes |
| L6 | TE5 on both arms (the spec names no arm) or A2 only? | both; TE5 passes only when both arms pass |
| L7 | TE5: does a family peptide with an accepted tolerant placement violate "land in the contaminant bucket"? (M3-15 keeps any accepted placement out of the bucket; an exact genome match is exempt.) | yes, a violation (`placed_tolerant`) |
| L8 | TE7 pair: compare against the pass-1 A2 attempt when it is fresh and code-equivalent (no diff under `src/`, `workflow/`, `pyproject.toml`, `uv.lock`, `giap_run.sbatch`), else two new replicates on one commit? | yes |
| L9 | Cross-node LAPACK differences (M4-20 decision 9) can make a cluster TE7 fail on `regions_stats.tsv` and the join-count columns, because tasks land on different CPU classes. Count a `known_lapack`-only result as a pass, pin both runs' tasks to one CPU class (a profile change, out of scope here), or report it as a fail with that diagnosis? | report `fail` with the `known_lapack` diagnosis and the M4-20 decision 10 backlog link (TE7 says byte-identical; R3); pinning stays a backlog option |
| L10 | DIA-NN replicate: rerun both steps (library prediction too) in a new `OUT_DIR`, 16 CPUs, `--constraint=CPU_SKU:9754` (the class of Slurm 1879956), so a difference is DIA-NN's own and not the hardware's? | yes |
| L11 | If DIA-NN is nondeterministic, run GIAP a third time on the second converted table (a sensitivity analysis of A2's conclusions to search noise)? | no in this issue: the L3 aggregates are reported; a backlog item if the differences are not trivial |
| L12 | TE7 gate set: `tables/`, `mapping/`, top-level `annotation/`, `plots/data/`, `schemas/datapackage.json`; external tool outputs and figures reported but not judged; manifest, logs, report and preflight excluded? | yes |
| L13 | Enforce "investigated and explained" in `render` (exit 1 when a non-pass TE has no diagnosis section with a cause class)? | yes |
| L14 | Size: 310 lines is above the index's 280. One PR, or TE7 split now? | one PR; split only if review pushes past ~400 |

No new question for Khoa: what may be committed stays M10-05's Q-K1 (default: nothing derived from
the pilot runs, only job IDs, states and the five statuses in PROGRESS). A `search` diagnosis that
asks for different DIA-NN settings would be an ADR-0051 amendment for him.
