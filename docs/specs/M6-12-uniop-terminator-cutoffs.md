# M6-12: UniOP AUROC and cut-off, the distance-model refit, the TerminatorNet cut-off, ADR-0022

Milestone M6 · track V (validation) · **ready** (M6-11 merged as #141, 2026-10-21) · depends on:
M6-11 (truth sets), merged M2-12 (UniOP adapter, distance model), M2-15 (TerminatorNet parser,
terminator windows) · **statistical** (label `stats`) · four PRs (see "Split"): about 140 code lines
in `src/giap/stats/`, about 230 in `validation/`, about 25 in the protected-file PR · produces
ADR-0022 and, after Khoa accepts it, a protected-file change.

Release of record: **RegulonDB 14.5.0** (release 2026-01-28; Khoa, 2026-10-21, PROGRESS.md ›
Questions item 1). GIAP_SPEC.md §15.3 and §7.4.6 say "v12"; the spec is read-only, and 14.5.0 is
the release that was fetched and is cited.

## Spec sections implemented

- §15.4 steps 3–5, verbatim: UniOP — "Record the AUROC. Choose the smallest threshold with
  precision ≥ 0.90 (if none reaches 0.90, the threshold with maximum F1). Refit the
  distance-model parameters d₀ and s." TerminatorNet — "Precision = the share of calls in
  intergenic regions downstream of genes that lie within ±50 bp of a RegulonDB terminator on the
  same strand. Choose the smallest confidence with precision ≥ 0.90 (else maximum F1). RegulonDB
  is incomplete, so this precision is a lower bound; say so in the ADR." Step 5: "Record the
  cut-offs in `src/giap/data/calibration.yaml`, with RegulonDB version, date, AUROC, precision and
  recall. The ADR goes to Khoa (protected file)."
- §7.4.6 "Truth check": the UniOP-alone AUROC is the third clause of the M6 gate (§7.4.7, §13 row
  M6). The UniOP-plus-co-abundance comparison needs K-12 proteomics and stays deferred.
- §7.1.3: the distance model p(operon | d) = 1 / (1 + exp((d − d₀)/s)), defaults d₀ = 50, s = 20,
  "refitted on RegulonDB K-12 pairs during calibration (§15.4) and stored in
  `src/giap/data/calibration.yaml`".
- §7.1.5 and §7.4.5: "a TerminatorNet terminator at or above the calibrated confidence cut-off".

## Inputs, as they exist now

Every input is already on disk. **No tool needs to run again** (see "Where compute runs").

| Input | Path | Identity |
|---|---|---|
| Pair truth | `$GIAP_DATA_DIR/regulondb/14.5.0/truth/pair_truth.tsv` | columns `gene_a, gene_b, label, tu_ids`; `label ∈ {positive, negative, ambiguous}`; `tu_ids` `;`-joined; one row per same-strand adjacent pair of `giap.annotate.operons.adjacent_pairs` on contig `U00096.3`, wrap pair included |
| Terminator truth | `$GIAP_DATA_DIR/regulondb/14.5.0/truth/terminator_truth.tsv` | columns `contig, start, end, strand, source_id`; **0-based half-open**; 375 rows (177 `+`, 198 `-`) |
| Truth summary | `$GIAP_DATA_DIR/regulondb/14.5.0/truth/truth_summary.json` | `regulondb_version` 14.5.0, `release_date` 2026-01-28, and the SHA-256s of `TUSet` (c96c8a3b…), `TerminatorSet` (9d916ffe…) and `GeneProductSet` (3c2c6be3…) |
| Gene model | `tests/data/recorded/bakta/1.12.2-full/k12/k12_mg1655.json.gz` + `tests/data/genomes/k12_mg1655.fna.gz` | the **circular-tagged** Bakta 1.12.2 full-DB run (Slurm 1849096); 4,373 genes; replicon `U00096.3`, circular, 4,641,652 bp. M6-11 built the truth from this same model. |
| UniOP run | `$GIAP_ROOT/recordings/uniop/k12_mg1655/{uniop.pred, uniop_cds.gff}` | pin 3ed14dfc (ADR-0008), Slurm 1849513_8; sha256 `uniop.pred` 7d852a93…, `uniop_cds.gff` ea587089…; 4,372 input genes |
| TerminatorNet run | `tests/data/recorded/terminatornet/721803095459/k12/terms.tsv.gz` (committed, verbatim `terminatornet.tsv`) | commit 72180309 (ADR-0008), `-s -n 17`, wheat-04, wall 3,603 s (job 1847778); genome `genome.norm.fa` sha256 918f35c7…, the same file Bakta's preflight wrote; 8,294 calls (4,166 `+`, 4,128 `-`) |

M6-11's counts (truth_summary.json): 3,067 pairs, 1,730 positive, 1,009 negative, 328 ambiguous
(10.7 %), so **2,739 labelled pairs**. Gene mapping share is 0.926.

### The UniOP recording's `uniop_map.tsv` must not be used

The recorded UniOP input was built from the **non-circular** Bakta run
(`$GIAP_ROOT/recordings/bakta_full/k12_mg1655`). In that run `KOIAIB_00001` is the CDS at
190..255. The circular-tagged run that the truth uses calls one more CDS, at 3..98, and that CDS
takes `KOIAIB_00001`. **Every later locus tag shifts by one.** So the recorded
`uniop_map.tsv` (`U00096.3_1 → KOIAIB_00001`) labels each UniOP gene with its neighbour's tag.

Joining on that map makes 1,873 of the 2,739 labelled pairs match **by key, with the wrong
probability**. The other 866 do not match at all. This was checked in-session.

**Rule:** the driver rebuilds the map **by coordinates**. Each `uniop_cds.gff` row (1-based
closed `start..end`, strand) is matched to the gene model's gene whose
`(interval.start + 1, interval.end, strand)` is equal. Every GFF row must match exactly one gene;
otherwise the driver stops with `ValueError` naming the first row that fails.

Checked in-session: 4,372 of 4,372 rows match. The one model gene without a UniOP row is
`KOIAIB_00001` (3..98). With this map, all 2,739 labelled pairs have a UniOP probability. Only two
ambiguous pairs lack one: `KOIAIB_00001–00002` and the wrap pair `KOIAIB_04716–00001`.

The driver reads `uniop.pred` with `giap.annotate.uniop_adapter.parse_uniop_pred(pred, gff,
id_map)`, using the rebuilt `id_map`.

## How each number is computed

Notation: L = the labelled pairs (label ≠ ambiguous) that have a UniOP p; y = 1 for positive, 0 for
negative; P = Σy, N = |L| − P. A pair is **predicted operon at threshold t iff p ≥ t**. This is the
comparison the pipeline uses everywhere: `classify` in `operons/pairs.py`, `build_tus` in
`annotate/operons.py`, and the flag code.

### 1. UniOP AUROC with its CI

- **Point estimate.** The Mann–Whitney AUROC: AUROC = (1/(P·N)) Σ_{i∈pos} Σ_{j∈neg} [p_i > p_j] +
  ½[p_i = p_j]. Ties count one half. UniOP's output is heavily tied (153 distinct values over
  2,739 labelled pairs; 395 pairs at exactly 0.010000), so the tie rule matters. Compute it
  through ranks (average ranks), not a pairwise loop.
- **CI: a 95 % percentile bootstrap that resamples whole directons.** A directon is a maximal run
  of consecutive genes on one strand in genome order. On a circular replicon the last and first
  runs join when they share a strand. Every same-strand adjacent pair belongs to exactly one
  directon. Pairs inside a directon are correlated, because operons span several pairs, so
  resampling single pairs would make the CI too narrow.
  - Use `giap.stats.bootstrap.tu_block_bootstrap(df, tu_col="directon_id", stat_fn, B=2000,
    rng, level=0.95)`, without strata.
  - A replicate with only one class gives NaN and is dropped. Fewer than B/2 valid replicates give
    a CI of (NaN, NaN), which `tu_block_bootstrap` already does.
  - `rng = giap.rng.for_step("m6_12_cutoffs", seed)`, with `seed` = the config default
    (20261014).
- **Reported:** the estimate and the CI to 3 decimals in docs and in full precision in JSON, with
  P, N and the number of directons.

### 2. The cut-off rule, the same for UniOP and TerminatorNet

- **Candidates.** The distinct observed scores, sorted in descending order. For each candidate t:
  TP(t) = #{y = 1, score ≥ t}, FP(t) = #{y = 0, score ≥ t}, precision = TP/(TP+FP),
  recall = TP/P, and F1 = 2TP/(2TP+FP+FN).
- **Rule.** t* = the **smallest** candidate t with precision(t) ≥ 0.90. If no candidate reaches
  0.90, t* = the candidate with maximum F1; among equal F1 values, take the smallest t.
- **Exact arithmetic.** "Precision ≥ 0.90" is tested in integers as 10·TP ≥ 9·(TP+FP), and F1
  ties are compared as fractions, so a precision of exactly 0.9 counts.
- **Literal reading.** This is the §15.4 text taken literally. Precision is not monotone in t,
  so the rule can choose a t whose neighbour just above it is below 0.90 (example 4). The
  alternative "stable" reading is question 1.
- **Recorded value.** t* is the observed score itself, written as Python's shortest round-trip
  `repr` of the float. It is **never rounded up**: TerminatorNet has 2,304 calls at exactly
  0.8812014060630281, so rounding to 0.882 would silently drop all of them through the `≥`
  comparison.
- **Reported at t*:** precision, recall, F1, TP, FP and the rule that fired (`precision_0.90` or
  `max_f1`).
  - UniOP precision and recall get 95 % directon-bootstrap CIs, with t* **held fixed** (not
    re-selected in each replicate).
  - TerminatorNet precision and recall get Wilson 95 % intervals.
  - The ADR notes that the precision at t* is optimistic, because t* was chosen to reach 0.90.
- **The full curve** (every candidate's TP, FP, precision, recall and F1) is written beside the
  JSON (see "Outputs").

### 3. The d₀ / s refit

- **Data.** The same L, with d = `giap.annotate.operons.intergenic_distance(a, b, replicon)`.
  This is `b.start − a.end` in genome order, so it is negative for overlaps, and it is taken
  across the origin on a circular replicon. Every labelled pair is used, whatever its UniOP p
  is.
- **Model.** y ~ Bernoulli(σ(β₀ + β₁d)), fitted by unpenalised maximum likelihood with
  Newton–Raphson (IRLS). Start at the current constants (β₀ = 50/20, β₁ = −1/20); stop when
  max|Δβ| < 1e-10, with at most 100 iterations. Then **d₀ = −β₀/β₁** and **s = −1/β₁**, because
  σ(β₀ + β₁d) = 1/(1 + exp((d − d₀)/s)).
- **Failure.** Any of the following is a refit failure: no convergence, β₁ ≥ 0 (operon odds not
  falling with distance), a non-finite β, or quasi-separation (|β₁| > 1, i.e. s < 1 bp). A failure
  keeps d₀ = 50 and s = 20 with `calibrated: false`, and the ADR says so.
- **CI.** The same directon bootstrap (same B and RNG stream, drawn after the AUROC replicates)
  refits each replicate. A failed replicate gives NaN and is dropped; the CI is the 95 %
  percentile CI of d₀ and of s.
- **Distance-model AUROC.** Reported beside UniOP's, because it measures the fallback's quality.
  The model is monotone decreasing in d, so its AUROC is the Mann–Whitney AUROC of **−d**. It does
  not depend on the fitted values, and it gets the same bootstrap CI.
- **Also reported:** the log-loss of the refitted model and of the placeholder (50, 20) on L, so
  the ADR can show the gain.

### 4. The TerminatorNet cut-off

- **Calls.** `giap.annotate.terminators.parse_terminatornet(terms.tsv.gz, replicons)` gives
  0-based half-open intervals on the forward strand with the strand and `confidence` (= `pred`).
  Only contig `U00096.3` is used.
- **Downstream windows: what "in intergenic regions downstream of genes" means.**
  - Every gene g of the model (CDS, sORF, pseudogene, as in `genes.parquet`) has one downstream
    window. It runs from g's 3′ end, moved 50 bp back into g (the §7.1.5 slack), up to the
    nearest end of g's **genome-order neighbour on either strand** in the direction of
    transcription.
  - This is `giap.annotate.terminators._window(g, neighbour, replicon, 50)`.
    - For a `+` gene the neighbour is the next gene in genome order, and the window is
      `[g.end − 50, neighbour.start)`.
    - For a `-` gene the neighbour is the previous gene, and the window is
      `[neighbour.end, g.start + 50)`.
  - On a circular replicon, the neighbour wraps past the origin and the window is unrolled.
  - On a linear contig, a gene with no neighbour gets 300 bp past its stop
    (`DOWNSTREAM_NO_NEXT_BP`), clipped to the contig.
  - A window with `lo ≥ hi` is empty.
- **Eligible calls.** A call on strand s is eligible iff its whole interval lies inside a window
  of some gene on strand s. "Inside" means `lo ≤ start + k` and `end + k ≤ hi` for a shift k in
  {0, ±length} on a circular replicon. This is `_max_in`'s containment rule.
- **Match.** An eligible call matches when some RegulonDB terminator on the same contig and
  strand is within **±50 bp**: the interval gap max(0, r.start − c.end, c.start − r.end) ≤ 50.
  Shifts of ±length are tried on a circular replicon. An opposite-strand terminator never matches.
- **Precision and the cut-off.** Labels are y = 1 for matched and y = 0 for unmatched eligible
  calls. Precision(t) = matched eligible calls with confidence ≥ t divided by eligible calls with
  confidence ≥ t, which is §15.4's "share of calls". Several calls near one terminator each count.
  t* comes from rule 2, over the distinct confidences of the eligible calls.
  - **The ADR states that this precision is a lower bound**: an unmatched call may be a real
    terminator RegulonDB lacks.
  - TerminatorNet AUROC (matched vs unmatched eligible calls) is reported as a diagnostic only.
    It does not appear in calibration.yaml.
- **Recall.** Eligible truth = RegulonDB terminators that lie inside a same-strand window, by the
  same containment rule. Recall(t) = eligible truth matched by at least one eligible call with
  confidence ≥ t, divided by eligible truth. The recall against all 375 terminators is reported
  beside it.
  - RegulonDB includes Rho-dependent terminators, which TerminatorNet does not model (§7.4.5), so
    recall is low by construction. Recall split by RegulonDB `type` is reported. The type comes
    from `TerminatorSet` through `validation/regulondb_truth.verified_release` and
    `parse_terminator_set`, because `terminator_truth.tsv` does not carry it.

## Where compute runs, and what it costs

| Work | Where | Cost |
|---|---|---|
| UniOP on K-12 | **no new run**: the recording (Slurm 1849513_8) is reused, with the map rebuilt by coordinates | 0. A fresh run on the circular-tagged Bakta model (question 2) would be one array task of `scripts/slurm/nf_process_real.sbatch` (task 8, UNIOP; 7 s of tool wall time; the job reserves 16 CPUs / 64 GB for under 10 min, about 3 CPU-h reserved). That run first needs an input override, because the script reads `recordings/bakta_full/`, the non-circular run. |
| TerminatorNet on K-12 | **no new run**: the committed recording (job 1847778) is reused. Its input is the genome alone, the genome sha256 equals Bakta's preflight genome, and the commit equals the ADR-0008 pin. | 0. A re-run is needed only if the pin or the genome changes: 60 min wall on 17 CPUs, about 17 CPU-h, via `record_tools.sbatch terminatornet`, job name `giap-*`. Count it against the R14 128-CPU budget first. |
| The calibration itself (`validation/calibrate_cutoffs.py`, B = 2,000) | session allocation, 1 CPU | estimated under 2 min: 2,000 × (AUROC on 2,739 pairs + one IRLS fit) and 8,294 × 375 interval checks. If a real run exceeds 10 min, it moves to `sbatch` (spec §3.8). |
| Fast suite after PR-C | Slurm, `scripts/slurm/test_branch.sbatch` | about 13 min (Slurm 1864288 took 13:10) |

## Outputs

`validation/calibrate_cutoffs.py --regulondb-dir $GIAP_DATA_DIR/regulondb/14.5.0 --bakta-json …
--genome … --uniop-dir $GIAP_ROOT/recordings/uniop/k12_mg1655 --terminatornet
tests/data/recorded/terminatornet/721803095459/k12/terms.tsv.gz --out-dir
$GIAP_DATA_DIR/regulondb/14.5.0/calibration` writes the files below. None of them is committed
(R6/R11). The ADR quotes the numbers and the SHA-256 of `cutoffs.json`.

- `cutoffs.json`:
  - `regulondb_version`, `release_date`, the three truth SHA-256s, `date` (UTC of the run), `seed`,
    `B`, `cluster: "directon"`, `n_directons`;
  - `uniop`: `{pin, pred_sha256, gff_sha256, n_labelled, n_positive, n_negative, n_without_p,
    auroc, auroc_ci, threshold, rule, precision, precision_ci, recall, recall_ci, f1, tp, fp}`;
  - `distance_model`: `{d0, d0_ci, s, s_ci, converged, auroc, auroc_ci, logloss,
    logloss_placeholder}`;
  - `terminatornet`: `{commit, terms_sha256, n_calls, n_eligible, n_truth, n_truth_eligible,
    threshold, rule, precision, precision_ci, precision_is_lower_bound: true, recall, recall_ci,
    recall_all, recall_by_type, auroc_diagnostic, tp, fp}`.
- `uniop_curve.tsv` and `terminatornet_curve.tsv`, with columns `threshold, tp, fp, precision,
  recall, f1`, in descending threshold order.
- `terminatornet_calls.tsv`: every eligible call, as `contig, start, end, strand, confidence,
  gene, matched_source_id` (empty when unmatched). Reviewers can audit matches from it.

All JSON floats are finite or `null`, written with `giap.manifest.finite`. Re-running with the same
inputs and seed gives a byte-identical `cutoffs.json` apart from `date`.

## What changes in protected files, and only after ADR-0022 is accepted

The run reads its thresholds from **two** places:
- `annotation.uniop.pair_threshold` and `annotation.terminators.terminatornet_yes` in
  `src/giap/config/defaults.yaml` (and the matching defaults in `src/giap/config/model.py`,
  enforced by `test_model_defaults_equal_defaults_yaml`);
- `operon_distance` in `src/giap/data/calibration.yaml`, read by `steps/operons.py`.

`calibration.yaml`'s `terminatornet_yes.value` is not read by any code today. The manifest
copies `{calibrated, provenance}` for each `calibration.yaml` key and emits the run flag
`uncalibrated:<key>` while `calibrated` is false (MVP-03c). Note that ADR-0021 Amendment 1
dropped a different flag, the per-row `uncalibrated` flag.

PR-C therefore changes:

```yaml
# src/giap/data/calibration.yaml (PR-C)
operon_distance:
  d0: <fit>            # or 50 with calibrated: false if the refit fails
  s: <fit>
  provenance: "RegulonDB 14.5.0 K-12 refit, ADR-0022 (M6-12)"
  calibrated: true
  regulondb_version: "14.5.0"
  date: "<run date>"
  auroc: <distance-model AUROC>
  n_pairs: 2739
uniop_pair_threshold:  # new key, mirrored in defaults.yaml annotation.uniop.pair_threshold
  value: <t*>
  provenance: "RegulonDB 14.5.0 K-12, ADR-0022 (M6-12)"
  calibrated: true
  regulondb_version: "14.5.0"
  date: "<run date>"
  rule: precision_0.90 | max_f1
  auroc: <UniOP AUROC>
  precision: <…>
  recall: <…>
  n_pairs: 2739
terminatornet_yes:     # mirrored in defaults.yaml annotation.terminators.terminatornet_yes
  value: <t*>
  provenance: "RegulonDB 14.5.0 K-12, ADR-0022 (M6-12); precision is a lower bound"
  calibrated: true
  regulondb_version: "14.5.0"
  date: "<run date>"
  rule: precision_0.90 | max_f1
  precision: <…>       # lower bound
  recall: <…>
  n_calls: <eligible calls>
```

PR-C also changes `defaults.yaml` (two values) and `config/model.py` (two defaults). Nothing else
changes there: no code changes and no schema changes. The manifest schema already accepts any
key with `{calibrated, provenance}`.

## Split (one PR per row; counted lines = ADR-0007 code under `src/`, `workflow/`, `scripts/`)

| PR | Branch | Content | Roles | Counted lines | Waits for |
|---|---|---|---|---|---|
| A | `m6/cutoffs-stats` | (1) Statistician: a methods section "K-12 cut-off calibration (M6-12)" in `docs/methods/operons.md` (rules 1–4 above, the directon bootstrap, the lower-bound caveat), reviewed first (`stats`). (2) Test-author: unit tests on synthetic data (examples 1–5, degenerate rows marked A). (3) Implementer: `src/giap/stats/cutoffs.py`, pure functions: `auroc(scores, y)`, `threshold_curve(scores, y)`, `choose_cutoff(curve, precision=0.90)`, `fit_distance_model(d, y) -> DistanceFit`. | statistician → TA → impl | ≈ 140 | nothing |
| B | `m6/cutoffs-k12` (stacked on A, or after A merges) | (1) TA: tests for the driver on toy genomes (examples 6–9, the remaining degenerate rows), all offline. (2) Impl: `validation/calibrate_cutoffs.py` (coordinate re-map, directons, windows and matching, the bootstrap, writers, CLI). (3) The real run in the session; numbers go into PROGRESS.md and the methods note. (4) **`docs/adr/ADR-0022-uniop-terminatornet-cutoffs.md`, status `proposed`, Needs Khoa: yes**, with every number. (5) `changelog.d/` fragment. | TA → impl → lead | 0 (`validation/` is not counted; ≈ 230 actual) | A |
| C0 | `m6/calibration-tests-agnostic` | Lead only (role guard: existing tests change): rewrite the tests pinned to the placeholder values so they read from `calibration.load()`. These are `tests/unit/test_operons.py` ≈ l. 530–548 and the `{"operon_distance", "terminatornet_yes"}` key sets and `uncalibrated:*` flag assertions in `tests/unit/manifest/test_content_cli_{local,nextflow}.py`. Add one test: for every `calibration.yaml` entry with `calibrated: true` and a `defaults.yaml` mirror, the two values are equal. It passes vacuously today. | lead | 0 (tests) | nothing; may land in parallel with A |
| — | (on main) | Khoa sets ADR-0022 `Status: accepted-by-khoa`, as a tiny docs PR like #249 | Khoa / lead | 0 | B |
| C | `m6/calibration-adr0022` | `calibration.yaml` (3 entries), `defaults.yaml` (2 values), `config/model.py` (2 defaults). Every commit carries `Gate-Change-Approved: ADR-0022`. The fast suite runs on Slurm. | lead or impl | ≈ 25 | ADR-0022 accepted on main, C0 merged |

Everything except PR-C lands without waiting for Khoa. If the TerminatorNet half degenerates (see
below), PR-C carries only the UniOP and distance entries, and the TerminatorNet entry stays
`calibrated: false`.

## Examples (input → expected output)

| # | Input | Expected |
|---|---|---|
| 1 | AUROC with ties: positives p = {0.9, 0.8, 0.5}, negatives p = {0.5, 0.2} | AUROC = 5.5/6 = **0.9167**. The tie 0.5 = 0.5 counts ½. |
| 2 | Rule, a threshold reaches 0.90: (p, y) = (0.99,1) (0.95,1) (0.90,1) (0.90,0) (0.80,1) (0.70,1) (0.60,0) (0.50,1) (0.40,0) (0.10,0) | Curve: precision is 1.0 at 0.99, 1.0 at 0.95, 0.75 at 0.90, 0.80 at 0.80, 0.833 at 0.70 … **t\* = 0.95**, `rule = precision_0.90`, precision 1.000, recall 0.333 (TP 2, FP 0). The two pairs tied at 0.90 enter together. AUROC = 0.8125. |
| 3 | Rule, no threshold reaches 0.90: (0.9,0) (0.8,1) (0.7,1) (0.6,0) (0.5,1) | Maximum precision is 0.667. F1 is 0, 0.4, 0.667, 0.571 and 0.75, so **t\* = 0.5**, `rule = max_f1`, precision 0.600, recall 1.000. The ADR states that 0.90 was not reached. |
| 4 | Literal rule with non-monotone precision: eight positives at 0.99, 0.98, …, 0.92, then (0.91,0), (0.90,1), (0.50,0) | Precision is 1.0 down to 0.92, then 8/9 = 0.889 at 0.91, then **exactly 9/10 at 0.90**, then 0.818 at 0.50. The literal rule gives **t\* = 0.90** (precision 0.900 by the integer test 10·9 ≥ 9·10; recall 1.000). The "stable" reading of question 1 would give 0.92. |
| 5 | Distance refit on expected counts: d ∈ {−20, −10, …, 200}, 1,000 pairs at each d, round(1000·σ) positives from d₀ = 40, s = 15 | converged; \|d₀ − 40\| < 0.2 and \|s − 15\| < 0.1 (the in-session check gave 39.991 and 14.962). With all y = 1 below d = 60 and all y = 0 above (separation), the fit fails: `converged: false`, and the yaml keeps 50 / 20, `calibrated: false`. |
| 6 | TerminatorNet, linear 10 kb contig. Genes g1 + [100,1000), g2 + [1500,2500), g3 − [3000,4000), g4 − [5000,6000). Windows: g1 [950,1500), g2 [2450,3000), g3 [2500,3050), g4 [4000,5050). Calls: c1 + [1100,1140) 0.95; c2 + [2460,2490) 0.90; c3 − [2600,2640) 0.85; c4 + [3500,3540) 0.99; c5 − [4100,4130) 0.60; c6 + [1460,1520) 0.99. Truth: T1 + [1120,1150), T2 + [2560,2580), T3 − [2680,2700), T4 + [4100,4130). | Eligible: c1, c2, c3, c5. c4 is in no `+` window, and c6 crosses g2's start. Matches: c1 to T1 (overlap); c2 none (gap to T2 = 70 > 50); c3 to T3 (gap 40, minus strand); c5 none (T4 is on the other strand). Precision is 1.0 at 0.95, 0.5 at 0.90, 0.667 at 0.85, 0.5 at 0.60, so **t\* = 0.95**, precision 1.000 (lower bound), recall 1/3. Eligible truth is T1, T2, T3 (T4 lies in no `+` window); `recall_all` = 1/4. |
| 7 | Circular 10 kb replicon: gA + [200,800), gB + [9500,9900); a call + [20,60) at 0.8; truth + [80,100) | gB's neighbour wraps to gA, so gB's window is [9450, 10200) unrolled. The call shifted by +10,000 is [10020,10060), inside the window, so it is **eligible**. The gap to the truth is 20, so it **matches**. The pair gB→gA has d = 200 − 9900 + 10000 = **300** in the refit. |
| 8 | UniOP join on a toy recording whose map is shifted by one (GFF row `c_1` at 190..255 +, map `c_1 → T1`, but the model has T1 at [2,98) and T2 at [189,255)) | The coordinate map gives `c_1 → T2`; `uniop_map.tsv` is ignored. A GFF row with no exact coordinate match, or one matching two genes, raises `ValueError` naming that row. On the real K-12 inputs: 4,372 / 4,372 matched, `n_without_p` = 0 among labelled pairs. |
| 9 | The real K-12 run (PR-B) | `cutoffs.json` with `n_labelled = 2739`, `n_positive = 1730`, `n_negative = 1009`; the UniOP AUROC with CI to 3 decimals; t\* and its rule; d₀ and s with CIs; the TerminatorNet t\* with `precision_is_lower_bound: true`. Running again with the same seed gives the same file apart from `date`. |
| 10 | A `giap run` after PR-C | `manifest.calibration` has three keys, each `calibrated: true` with the ADR-0022 provenance. `flags` has no `uncalibrated:*` entry. TUs use the new `pair_threshold`; terminator "yes" uses the new `terminatornet_yes`. |

## Edge and degenerate cases

| Case | Behaviour | PR |
|---|---|---|
| Labelled pair without a UniOP p | excluded from L and counted in `n_without_p`. Never filled from the distance model: the gate needs the UniOP-alone AUROC. If more than 1 % of labelled pairs lack a p, the run escalates to the lead instead of writing numbers. | B |
| Fewer than 500 labelled pairs, or P = 0 or N = 0 | `ValueError`; no calibration is claimed (M6-11 example 5 guard). On K-12 there are 2,739. | A, B |
| All scores equal | AUROC = 0.5; one candidate, so t\* = that value by max F1 unless precision ≥ 0.90 | A |
| Ties at a candidate | candidates are distinct values, so all tied items enter together (examples 2 and 4) | A |
| Precision is exactly 0.90 | counts as reached (integer test) | A |
| Equal maximum F1 at several t | the smallest t | A |
| Non-monotone precision | literal rule (example 4); see question 1 | A |
| A bootstrap replicate with one class, or a failed refit | NaN, dropped; fewer than B/2 valid gives a (NaN, NaN) CI, reported as `null` | A, B |
| Distance refit: no convergence, β₁ ≥ 0, non-finite β, or s < 1 bp | `converged: false`; the yaml keeps 50 / 20 with `calibrated: false`; the ADR says so | A, C |
| Minus strand | pairs: d is computed in genome order and is strand-independent. Terminators: windows mirror (`[neighbour.end, g.start + 50)`), and matching requires the same strand (example 6). | B |
| Circular replicon | the wrap pair's d crosses the origin; windows and matching try ±length shifts (example 7). On K-12 the wrap pair is ambiguous, so it is excluded, but the code path must handle a labelled wrap pair. | B |
| Linear contig, gene at a contig end | a window of 300 bp past the stop, clipped to the contig | B |
| Overlapping or touching genes | the window is only the 50-bp tail (`lo ≥ hi` gives no window) | B |
| A call overlapping a window edge | not eligible (containment, example 6 c6) | B |
| Several calls match one truth terminator | each counts for precision (§15.4 "share of calls"); the terminator counts once for recall | B |
| No eligible TerminatorNet calls, or zero matches at every t | the TerminatorNet half is not calibrated: `terminatornet_yes` stays 0.5, `calibrated: false`, and the ADR says why. The UniOP half still lands. | B, C |
| Empty `uniop.pred` (title and header only) | `parse_uniop_pred` gives `[]`, so n_without_p = all; escalate | B |
| TerminatorNet `NO TERMINATORS FOUND` lines | skipped by the parser | B |
| Recording identity mismatch | `ValueError` and no output. This covers: TerminatorNet commit ≠ ADR-0008 pin; genome sha256 ≠ `918f35c7…`; a UniOP sha256 ≠ the PROVENANCE value; truth SHA-256s ≠ `truth_summary.json`. | B |
| Contigs other than `U00096.3` | ignored, with a count logged (K-12 has none) | B |
| A later RegulonDB release | run again into `regulondb/<release>/calibration/`; new ADR numbers; the provenance fields make the change visible | — |

## Acceptance criteria (each checkable by a test or a file check)

- [ ] `giap.stats.cutoffs` reproduces examples 1–5 exactly. AUROC agrees with
  `sklearn.metrics.roc_auc_score`, used as a **test-only** oracle if it is in the dev extras, to
  1e-12 on 100 random tied inputs (Hypothesis).
- [ ] `choose_cutoff` uses integer comparisons: a property test over random (TP, FP) curves never
  disagrees with `Fraction` arithmetic.
- [ ] `fit_distance_model` recovers (40, 15) within the tolerances of example 5 and reports the
  separation case as `converged = False`.
- [ ] The driver reproduces examples 6–8 on toy inputs written by the tests. No test reads
  `$GIAP_DATA_DIR` or the network.
- [ ] The driver raises on every identity mismatch in the degenerate table, and on a shifted map
  whose coordinates do not match.
- [ ] Two runs with the same seed give `cutoffs.json` files that are equal apart from `date`
  (test on toy inputs).
- [ ] The real run's `cutoffs.json` has `n_labelled = 2739`, `n_without_p = 0` and a finite AUROC
  CI. The UniOP AUROC (3 decimals, with CI) is recorded in PROGRESS.md and in
  `docs/methods/operons.md`. This is the gate's third clause.
- [ ] ADR-0022 exists, with status `proposed` and Needs Khoa: yes. It contains every number in
  `cutoffs.json`, the sha256 of `cutoffs.json`, both rules that fired, and the sentence "RegulonDB
  is incomplete, so this precision is a lower bound". It also lists the downstream changes:
  `pair_threshold`, `terminatornet_yes`, d₀ / s, and the cleared `uncalibrated:*` run flags.
- [ ] After PR-C: `calibration.load()` returns the three entries with `calibrated: true` and the
  fields listed above. The C0 agreement test passes with non-vacuous content.
  `test_model_defaults_equal_defaults_yaml` passes. The fast suite passes on Slurm.
  `check_protected.py` passes with `Gate-Change-Approved: ADR-0022` against main, where ADR-0022
  is `accepted-by-khoa`.
- [ ] After PR-C, a `giap run` on SIM100K (example 10) has no `uncalibrated:*` flag. M6-09's fast
  calibration suites are unchanged: the simulations plant their own probabilities and pass explicit
  thresholds. Any test that read the default threshold implicitly is listed in PR-C.

## Hardening owed (not in M6-12)

- "Adding co-abundance does not lower K-12 operon AUROC" is deferred to the full validation by
  §13.
- How stable t\* is across bootstrap replicates (the share of replicates choosing the same
  candidate) would make a useful ADR footnote. It is left out to keep PR-B small.
- The committed M2-21 UniOP excerpt `tests/data/recorded/uniop/3ed14dfc/k12/uniop_map.tsv`
  carries the non-circular run's locus tags while `tests/data/recorded/bakta/1.12.2-full/k12/` is
  the circular-tagged run. No test was found joining the two by path, but any future join would be
  off by one. This goes to `docs/backlog.md` for the lead.
- The §8.6 `full` sensitivity variants stay at UniOP thresholds 0.3 and 0.7, as the spec text
  fixes them (§8.6), even if t\* moves far from 0.5. See question 7.

## Questions for the lead (each with the default used if unanswered)

1. **The precision rule's reading.** Literal "smallest t with precision ≥ 0.90" can pick a t
   whose next candidate up is below 0.90 (example 4). The alternative is "smallest t such that
   every t′ ≥ t also has precision ≥ 0.90". **Default: literal** (the §15.4 text). The ADR shows
   both values when they differ.
2. **UniOP input.** The default is to re-map the existing recording by coordinates (exact: 4,372 /
   4,372; no compute). The alternative is a fresh UniOP run from the circular-tagged Bakta model,
   which needs an input override in `nf_process_real.sbatch` and costs about 3 CPU-h reserved.
   **Default: re-map.** The only difference is gene `KOIAIB_00001`, and both of its pairs are
   ambiguous.
3. **Bootstrap unit.** Directon-cluster bootstrap or iid pairs? **Default: directon** (correlated
   pairs), B = 2,000, percentile 95 %. The ADR adds the iid CI as one line for comparison.
4. **"Within ±50 bp".** Interval gap ≤ 50 bp, or endpoint (3′-end) distance? **Default: interval
   gap** (overlap counts as 0).
5. **"Intergenic regions downstream of genes".** The default window is from the stop − 50 bp to the
   nearest genome-order neighbour on either strand, over genes in the gene model (CDS / sORF /
   pseudogene), with containment. Should Bakta RNA features (tRNA / rRNA operons, which RegulonDB
   terminates) also define windows? **Default: gene model only**, consistent with
   `terminator_between`.
6. **Terminator truth types and the recall denominator.** **Default: all RegulonDB types** count
   as matches (literal). Recall uses window-eligible truth, with `recall_all` and the per-type
   split reported. The alternative is to restrict matches to `rho-independent`.
7. **Sensitivity variants.** If t\* is far from 0.5 (for example ≥ 0.7), the §8.6 variants 0.3 and
   0.7 both fall below the default. **Default: leave them as the spec fixes them**, and file a
   backlog note for Khoa if t\* ≥ 0.7.
8. **Where `pair_threshold` is recorded.** The default is a new `uniop_pair_threshold` entry in
   `calibration.yaml` for provenance, with `defaults.yaml` staying the value the code reads; the C0
   test enforces that the two agree. The alternative is to make the code read `calibration.yaml`,
   which is a code change in operons, flags and sensitivity. **Default: the mirror plus the
   agreement test.**
9. **PR-C0 ownership.** Existing tests pinned to placeholders must change. Under the role guard
   that is lead work. **Default: the lead does C0 in parallel with PR-A.**
10. **Partial calibration.** If the TerminatorNet half or the refit degenerates, does PR-C land
    with only the halves that worked? **Default: yes.** Each entry has its own `calibrated` flag,
    and the gate's clauses are separable.

## Lead decisions (2026-10-23)
1. **Q1, precision rule:** the stable reading. The cut-off is the smallest t such that every
   candidate t′ ≥ t has precision ≥ 0.90, so a low-threshold spike in a non-monotone precision
   curve cannot set the cut-off. The ADR also shows the literal smallest t when the two differ.
   If no t qualifies, the max-F1 fallback applies as written.
2. **Q2–Q8, Q10:** the proposed defaults are accepted. That means:
   - Rebuild the UniOP map from GFF coordinates, and stop on any unmatched row.
   - Bootstrap over directons with B = 2,000, with the single-pair CI as a comparison line.
   - "Within ±50 bp" is the interval gap.
   - Windows come from the gene model only.
   - Match all terminator types, and report recall by type.
   - Keep the §8.6 variants; file a backlog note if the cut-off is ≥ 0.7.
   - Mirror the threshold in `calibration.yaml`, with an agreement test.
   - PR C lands with whichever halves worked.
3. **Q9:** PR C0 is lead work and runs in parallel with A.
4. **Hardening owed:** the shifted M2-21 UniOP excerpt goes to the backlog.
5. **Order:** this spec merges; then A (the statistician's methods section, then the
   test-author, then the implementer), C0 by the lead, then B. ADR-0022 then goes to Khoa as
   proposed, and C follows Khoa's acceptance.
6. **Statistician design review (46124d74; methods section "K-12 cut-off calibration (M6-12)"):
   B1 and B2 accepted; the methods section is now the authority wherever it is more specific than
   this spec.**
   - **B1, refit:** before fitting, run an exact 1-D separation check. Start from (logit ȳ, 0).
     Halve the Newton step until the log-likelihood does not fall (at most 30 halvings). A
     `LinAlgError` is a failure, not an exception. |β₁| > 1 is a backstop. The test that
     separates the two methods is d ∈ {−20, …, 400}, with 100 pairs per d drawn from d₀ = 100,
     s = 50, and must give 100.14 ± 0.05 and 50.06 ± 0.05. Preliminary K-12 fit: d₀ ≈ 102.7,
     s ≈ 47.2.
   - **B2, TerminatorNet precision:** add a placement null with R = 200. Each call moves to a
     random same-strand window and keeps its offset from the window's 3′ side. Report
     `coincidence_rate` c and `precision_lower_bound` = (obs − c)/(1 − c). The §15.4 sentence in
     the ADR is qualified accordingly.
   - **Example 4** under decision 1 expects 0.92. A top tie group below 0.90 sends the rule to
     max-F1, and that case has its own test.
   - **Accepted:**
     - R1: `precision_cv` (10-fold directon CV) and the bootstrap quantiles of t*.
     - R2: a window-cluster CI for TerminatorNet precision.
     - R4: the paired UniOP − distance AUROC difference with its CI.
     - R5: the ADR wording.
   - **R3:** TerminatorNet's F1 curve is flat (best precision 0.229, so the rule falls back to
     max-F1 and t* ≈ 1 − 4e-12). The rule stays as specified. Whether to use a one-SE rule instead
     changes §15.4, so it goes to Khoa in the ADR-0022 question. Pending Khoa; default: max-F1
     as specified.
   - **Attack nulls,** test-author, under `tests/calibration/attacks/`:
     - a directon-clustered null with AUROC 0.5 (CI coverage ≥ 93 % at R = 500);
     - the chance-match null;
     - a top-group spike that sends the rule to max-F1;
     - the hard refit start, plus separated and quasi-separated cases, none of which raises.
   - **For PR B:** `parse_terminatornet` reads plain text, so the driver decompresses
     `terms.tsv.gz` first.
