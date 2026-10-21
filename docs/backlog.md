# Backlog

Issues may live on GitHub; this file mirrors the milestone plan (GIAP_SPEC.md §13).

- [ ] M0 scaffold (in progress)
- [ ] M1 inputs and preflight
- [ ] M2 annotation wrappers, full Bakta DB, recorded fixtures
- [ ] M3 mapping and decoy FDR; pilot kit
- [ ] M4 gene table, detectability, region tests
- [ ] M5 autocorrelation and islands
- [ ] M6 operons; RegulonDB calibration (source currently unreachable)
- [ ] M7 variants and novel ORFs (declare GIAP's VCF FILTER values with descriptions; M1-05 writes generic ones)
- [ ] M8 function and determinants
- [ ] M9 visualisation and report
- [ ] M10 pilot validation
- [ ] M11 release 1.0

- M3: `mapping.tolerant.swaps` and `mapping.tolerant.mass_equivalents` config keys are not consumed by M3-07; decide which issue implements them or remove them (from M3-07 review, PR #48).
- M6: the simulator plants stop-codon readthrough (C-terminal extensions, ADR-0013) for cterm_extension calibration (from M3-19 test-author note).
- M2-09: antiSMASH databases with `url: "image"` (databases bundled in the image) — setup edge case from M2-08 review.
- M2-09 (from M2-08 review, low): delete the geNOMAD tarball after unpacking (§3.9); reject non-200 responses before hashing; cap redirects; test `_provenance` requires prev status ok; drop redundant tar `filter` weakening check.
- M2-09 follow-up: GO release directories can announce a different `data-version` date in the OBO header (2026-06-19 → 2026-06-15, 2026-08-05 → 2026-07-26), so `verify_go` rejects them and the pin stays at 2026-05-19; verify the header against a pinned `data_version` instead of the directory date.
- M4: migrate existing calibration tests to `_criteria.calibrate()` so ADR-0015 two-block confirmation and sharding apply to all suites (M4-05 decision 7; test-author).
- M4 (calibration migration follow-up): suites marked only `calibration` (stats_*, placement lists/genome) have no `calibration_full` variant, so the nightly full array skips them; adding one needs a new R — statistician decision.
- M4 (statistician, #69): gate the sim_proteomics q-value FDR check and the one-sided de novo coverage checks through `calibrate`.
- M4 (statistician, #69 minor): the 1.5q teeth test at q = 0.01 has power < 1 at fast R (seed-fragile) — raise R or use 2q at q = 0.01; teeth tests rerun in every shard task — run them once.
- Declare threadpoolctl (used by giap.genetable.detectability) as a direct dependency at the next uv lock update (M4-13).
- M3-17 follow-up: wire the blastp contaminant check for unassigned proteins into the mapping steps (only peptide library matching is wired).
- M3-17 follow-up: example 5 (nextflow executor end to end) is xfail(strict=False) until M2-20 (context_layers) merges; remove the xfail then (decision 16).
- M4-20: simulator `multi_identical_pairs` should give paired genes identical `aa` in the layout (M4-09 decision 2).
- M4-14: record the declared abundance type (e.g. maxlfq) in the input metadata; C10 has no scale column (M4-09 decision 3).
- M4-09: add a benchmarks/ entry for simulate_evidence on 4,500 genes (budget 3 s; unit test keeps a 30 s guard).
- M4-20 (with M4-09 decision 2): for shared_frac, the simulator's second hit sequence is not in that gene's aa (M4-09 review).
- M4-14: re-run the abundance Spearman sanity claim on M4-09's evidence simulator (logistic detectability, missed cleavages) once #80 is merged.
- M3-17/M4: resolve makedb protein IDs (GIAPORF_, NOVEL_, *_patch) through db_id_map.tsv in annotate_protein_hits instead of flagging cross_strain_protein (M3-23 review).
- M2: the full-DB Bakta recordings (bakta_full/*) were made without a replicon table, so every replicon is 'linear'; re-record with preflight's replicons.tsv (--replicons) so circular wraps are exercised (M3-23 review).
- M2: one shared IUPAC-aware reverse complement (ostir.py _COMP, bakta.py, operons.py only complement ACGTN; preflight allows IUPAC) (M2-16 review).
- M2: move the step helpers (_genes, _path, _replicons) from giap.steps.operons to a public shared module before more step modules import them (M2-16 review).
- M2-15: test the manifest `uncalibrated` flag and step-level `cross_check_transtermhp: false` once manifest and optional-output support exist.
- M2-15: parse_bactermfinder's input format is provisional (no recording); record a real BacTermFinder run (opt-in tool, env not built) and pin the parser to it.
- M4-20: replace M4-16's test stand-ins (p_detect_true, fixed COG table, tier/detectability columns) with the real M4-12/M4-13 outputs once merged.
- M4-14 (statistician): re-measure abundance.ibaq_denominator on M4-09 data before switching (fails 10–40 % at p_obs = 0.3); restate as a pooled mean Δρ if needed; add the p_obs = 0.3 attack.
- M4 notes: median imputation pulls detected-but-unquantified genes toward the centre; the min/exclude sensitivity reruns bracket it (statistician, #82).
- M7: compute the M4-16 QC metric id_mixups (§7.5.1: same-database ALT sites that are not start/stop sites); M4-16 reports NA + requires_variants until then.
- M4-16 (statistician #87 sign-off, non-blocking): (a) also flag ci_few_units by the number of TUs carrying a detection (non-coverage 0.07–0.12 for 1–5-gene TUs, 0.36 for 10-gene TUs at E≈10; job 1850360); (b) genome_inconsistent_rate only under --genome-relation same (§7.5.1), not related; (c) `hypothetical` parsed with astype(bool): the string "false" reads as True — parse booleans explicitly.
- Benchmarks: nothing reads benchmarks/baseline.json yet — add a check that compares benchmark results with it (M4-13 review).
- At the next uv lock update: scikit-learn>=1.8 (l1_ratio API) and threadpoolctl as direct dependencies (M4-13).
- M4-15/M4-17/M8 attack-2 suites must gate on M4-13 p_detect; even the correct model slightly over-predicts detection for the shortest genes (decile z ≈ −0.66) (statistician, #83).
- m3/mapping-nf: write the search_db auto-resolution event from the Nextflow `giap _step` path too (only --executor local writes it, M3-17 decision 14).
- M4-12 (statistician #92): giapsim mode with decoy-derived q (D/T) to separate MAYU with/without PEP; calibration run with q_value_max = 0.1; document 'best-scoring decoys lost' as a limitation of the 80 % pairing guard.
- M4-14: use Detection.protein_gene_map from M4-12 instead of re-deriving the map (statistician #92).
- M2-10: BacTermFinder recording (env not built) and record_fixture.py examples 1–3 deferred.
- M2-10 (#93 review): re-record the antiSMASH fixtures on the circular K-12 GBFF (current ones used the linear GBFF; ~18 s each) or note it in their PROVENANCE.
- M2-10 (#93 review, R6): scripts/fetch_test_genomes.py must verify checksums by default (NCBI Datasets md5sum.txt or mandatory --compare); tag_circular should check records instead of tagging all.
- TE6: pin Bakta's K-12 oriC call count at two (true origin + dnaA DnaA-box cluster); the replichore analysis (§1447) needs a rule for choosing the true origin.
- Spec note: Nissle GCA_043228065.1 has three plasmids (p15A-ILA-dapA, likely engineered; pMUT1; pMUT2) — the spec's edge case names two.
- Bakta fixture PROVENANCE: add attribution (Bakta DB CC BY 4.0, doi 10.5281/zenodo.14916843). Later recording waves need a size budget (wave A uses ~48 % of 50 MB).
- nf_process_real (job 1849282 K-12, 1849365 Nissle): PROMOTECH and UNIOP report 'failed' on K-12/Nissle only because no recorded inputs exist there. Needs (a) record_tools_inputs.py builders: Promotech targeted windows via the merged M2-13 promoter_windows step, UniOP CDS GFF from the Bakta JSON (format of tests/data/recorded/uniop/*/sim100k/uniop_cds.gff); (b) recordings of both tools on K-12/Nissle under $GIAP_ROOT/recordings (record_tools.sbatch: drop the sim100k-only promotech skip, add uniop); (c) the array's PROMOTECH/UNIOP branches reading them.
- Calibration log (tests/calibration/…/regions.jsonl, append-only): each record needs a run ID + commit so stale records (e.g. pre-pooling crispr hl_shift 'fail') are distinguishable (statistician #89).
- M4-17 nightly: the A-LIN full attack (test_alin_linear_replicon_hl_at_level) takes ~57 min per block on one core — shard it in the nightly/weekly calibration job (statistician #89 sign-off).
- M2-20 follow-ups (implementer report, real run 1851270): annotate.nf passes no `tool_dirs`, so OUT/annotation/tools/ is empty (§7.1.8 wants native tool outputs untouched); `cai`/`fop` all NA on sim100k layers — check whether expected (no reference set on the simulator) or a wiring gap; consider splitting schemas + gc_local_dev if review asks (658 gated lines).
- M2-20 (#104 review): annotate.nf passes no tool_dirs → OUT/annotation/tools/ empty (§7.1.8 native outputs untouched); nothing produces id_map, so --annotation GFF3/GBFF runs lack annotation/id_map.tsv. Wire both.
- M4-18 (lead decision 4): the dosage-adjusted plasmid row treats replicon depth as known (CI shifted by log2(plasmid / reference depth)). Deferred: propagate depth uncertainty (e.g. bootstrap or delta-method term from per-window coverage); the manifest should record a missing `--replicon-depth` (the flag itself is wired into `context_stats` since MVP-01).
- M4-17/M4-18 (lead decision 17; statistician #110 finding 7): a shared `ci_few_units` flag for region-test and plasmid TU-block CIs computed from fewer than 10 TUs (docs/methods/bootstrap.md: percentile CIs from < 10 TUs under-cover), as M4-16 already does for detection_summary.
- M4-15 decision 10: "file I/O only through io/" is enforced for genetable/assemble.py only. Step-level convention open: merged steps (mapping, context, operons, codon, ...) read Parquet directly with pd.read_parquet; decide whether steps must go through giap.io readers and add an AST check like test_assemble_module_is_pure if so.
- §9.5 per-field `unit` in the table schema resources (src/giap/schemas/tables/*.schema.json): none of them has it yet (M4-15 decision 16; project-wide).
- Promotech sim100k recording ($GIAP_ROOT/recordings/promotech/sim100k) is a genome-wide -g run with every row labelled chr1,plasmid_1; re-record in GIAP's -s window mode (M2-21 test-author).
- Nightly: test the newest open PR branches too (M2-21 decision 5, deferred). containers.yml / GHCR images deferred (M2-21 decision 2).
- giap run --annotation <bakta.json> does not stop BAKTA from running under Nextflow (M2-21 implementer); workaround --set annotation.precomputed.bakta_json=…. Make --annotation map to the precomputed path.
- Benchmarks: gate only steps with a baseline ≥ 1 s (M2-21 decision 10); implement the threshold in the benchmark check.
- M4-19 decision 7: the run manifest records `output.sensitivity = off`; decision 8: no rows for a module that is off. (Decision 16, registering and wiring `run_sensitivity`, is done: MVP-01.)
- M4-18 (#110 review): context_stats reads genes.parquet with pd.read_parquet; switch to M4-15's read_genes_parquet (output identical today) with the M5 context_stats work.
- Flaky timing asserts on loaded nodes: tests/unit/genetable/test_evidence.py::test_m409_default_fixture_detected_true_have_evidence (1.0 s; seen 1.01–1.13 s) and test_sim_layout::test_speed_1000_genes (0.2 s; seen 0.23 s). Move to benchmarks/ or loosen.
- Slurm test runs on barley nodes are ~8× slower (subprocess timeouts in test_cli_step_load, test_exact_review, test_criteria_review): raise those timeouts or exclude barley in scripts/slurm/test_branch.sbatch.
- M3-22 (#115 review, non-blocking): methods text: all 4 TE1 misses are selenocysteine (fdnG) peptides; split the weekly TE8 fetch and proxy into separate steps and send fetch stderr to te8.log; retry IncompleteRead and validate X-Total-Results in validation/fetch_uniprot_k12.py; dead branch at te_proxies.py:139-140; pooled_decoy_fdr docstring: valid only for search_db = other.
- M5-01 (#114 review, non-blocking test gaps): tied/unsorted 2-point pos accepted; keyword-only arguments untested; min_points < 3 unvalidated (LoessFit.k reports 1); loess_span_points does not validate n/frac; NaN in Blom censored_mask counts as censored.
- Calibration full tier exceeds 2 h in one job: shard it across array tasks.
- M4-19 (#117 statistician): sampled-shift region tests share one family stream across regions, so CRN breaks above `exhaustive_shift_max_genes`; give `test_region` (M4-17) a per-region shift stream. Sensitivity reruns compute 1,000-replicate bootstrap CIs that flagged sets never use; shrink or skip them (primary recomputed with the same count).
- PRIORITY (M5-05 decision 31, now load-bearing): shard the full calibration tier across nightly array
  tasks. The nightly full tier now carries three ~65-minute autocorrelation attack arms (a, b, c; 1,000
  replicates at 4,500 genes each, and a block-2 confirmation doubles an arm) on top of the main full
  suite and the M3 placement full checks (~38 min). M5-06's island attacks land next.
- local/dev-helpers/regen_registry.py writes only workflow/tests/data/annotate/step_registry.json; fold in
  the mapping twin (M3-25 regenerated it by hand with the same filter the test uses).
- M3-25 (#124 review, non-blocking): `_stratum0_pooled` takes n/k/r from the class table but mean PEP from
  all of P (unreachable with the default config); tests/helpers/placement_sims.py `prefilter=` silently
  discards `inject_short=`; the per-class π̂ gate has zero tolerance for a class seen in one replicate.
- Two homes for the K-12 reference proteome: scripts/slurm/fetch_resources.sbatch writes
  $GIAP_DATA_DIR/proteomes/UP000000625.fasta (used by the M3-24 pilot cfgs), while
  validation/fetch_uniprot_k12.py and scripts/slurm/weekly.sbatch use $GIAP_DATA_DIR/uniprot/. Pick one.
- MVP-01 pass 2 (found while wiring): with a precomputed `regions` layer ANNOTATE has no C7/C8
  Parquet, so CONTEXT (context_stats, sensitivity) skips itself silently — the run should say so,
  and ideally the precomputed layer should reach M4 as the contract tables.
- MVP-01 pass 2: nothing produces `id_map` under Nextflow (C14 is an optional input of
  `gene_table` and `context_stats`; the local executor stages `annotation/id_map.tsv`), and the
  `tables/` outputs are published but not validated inside the pipeline (contract assertions).
- MVP-01 pass 2: the registered `sensitivity` step writes a header-only `tables/sensitivity.tsv`
  when `output.sensitivity = off` (a declared step output must exist), while `run_sensitivity`
  leaves no file under the local executor. Decide which one the §9.4 manifest should describe.
- MVP-01 real run (jobs 1861957, 1861972, sim100k): every step runs and 1,646 peptide hits reach
  the gene table, but no gene is detected: the `decoy_picked` tier has no decoy win, so every
  eligible gene gets the picked-FDR floor q = 1/85 = 0.0118 > `stats.alpha.identification` = 0.01.
  An 87-gene genome cannot reach the identification alpha at all — the statistician should say
  whether the MVP evidence run needs a bigger genome, a different tier or a flag.
- MVP-01 pass 2: `detect` (M4-12) builds tier and q-threshold curve rows, but `gene_table` has no
  diagnostics output key, so those rows never reach `tables/detection_summary.tsv` (only the
  basic / qc / detectability sections are there).
  *M4-20 (decision 3) closed the Python half*: `gene_table` now writes the output `tier_summary`
  (`genetable/detection_summary_parts/tier.parquet`) and `context_stats` reads it as an optional
  input, so `run_context_stats` emits the `tier` section. The Nextflow half is open: passing it to
  CONTEXT needs an eleventh `take:` channel, and `workflow/tests/context.nf.test` binds exactly ten
  (`input[0..9]`), so the nf-test has to change with it — test-author work, not an implementer edit.
- M4-20: the §11.4 golden `regions_stats.tsv` cannot be byte-identical across machines. The
  `hl_shift` family (`hl_shift_log2`, `effect`, `ci_low`, `ci_high`) differs in the last 1-2 ulp per
  CPU because `giap.context.regions.residualise` (ADR-0018 decision 1) calls `np.linalg.lstsq`
  (LAPACK SVD); the three CI jobs (ubuntu 3.11, ubuntu 3.13, macOS 3.12) disagree with each other, so
  no recording machine satisfies the matrix. The other three golden tables and every other column of
  this one are byte-identical. Fix options: compare the golden's float columns numerically
  (test-author), round what the `regions_stats` writer emits (output-format decision, helps users
  too), or make the OLS BLAS-free (moves every calibration number). Evidence and bit patterns in
  PR #133.
- M4-20 (decision 2): `run_detectability` now writes `tables/genes.tsv` (`write_genes_tsv` had no
  caller at all), so the M4 chain and the local executor produce the §9.3 master table. It is *not*
  a declared step output yet, so `giap run` does not publish it: adding the key `genes_tsv` to
  `context_stats` (the step that owns `tables/`) breaks the hand-written `outputs=` dicts and the
  partial gene tables of `tests/integration/test_context_stats_step.py`, and adding it to
  `detectability` breaks `main_genetable_context.nf.test` ("tables/ published without C7/C8").
  Settle it with the publishDir item below; §6.2 gives `genes.tsv` to `finalize`, which takes the
  call over when it arrives.
- MVP-01 pass 2 (#129 review): `stepIn`'s `need` defaults to `[]`, so a step with no `need` emits
  an empty task — which is exactly the MVP-01 failure — as soon as one upstream channel is empty.
  `map_exact`, `map_tolerant`, `loci` and `parse_bakta` are each one empty upstream away from it.
  `need` should default to the step's required keys (the registry already lists them).
- MVP-01 pass 2 (#129 review): `OUT/genetable/genes_pre.parquet` is an intermediate published into
  the user's output tree, and §9.1 has no `genetable/` directory at all (it specifies
  `tables/genes.tsv`, which `giap run` still does not publish — see the M4-20 item above). Settle
  publishDir discipline — which keys are user-facing, which stay in work/ — with the M6 writers.
- MVP-01 pass 2 (#129 review): `cog/cog-20.def.tab` is a literal in `main.nf` (`cogTable()`), while
  Python reads the path from `db_sources.yaml` and the db manifest records it per item.
  `tests/unit/test_nf_main_installed_items.py` now pins the two literals against `SOURCES`, in the
  style of `test_nf_annotate_registry.py`; pass 2 should drop the duplication instead, by having
  `installed()` read the manifest entry's own `path`.
- `np.linalg.lstsq` in `giap.context.regions.residualise` is not bit-portable across CPUs (LAPACK SVD),
  so `regions_stats.tsv` float columns differ by 1–2 ulp between machines and between CI legs (M4-20
  decision 9, probe 1863133). Writer-level rounding is the real fix — see the §9.5 per-field precision
  item — and belongs with the M6 writers.
- CONTEXT cannot receive `tier_summary`: it would need an 11th `take:` channel and
  workflow/tests/context.nf.test binds exactly input[0..9] (M4-20 decision 3). Test-author work.
- tests/helpers/islands_data.py (m5/islands) and tests/helpers/oriter_data.py (m5/oriter) both carry a
  `tu_frame` and a giapsim-scenario driver; whichever merges second should deduplicate. Note the ori–ter
  driver deliberately keeps the simulator's real coordinates, while the island one re-lays genes through
  `Contig`, which would destroy the ori–ter geometry — so the shared helper must preserve both modes.
- Post-merge code review of #131 (M5-06 islands), pass-2 items (the origin-crossing TU bug is being
  fixed in its own PR): (1) the high-before-low tie-break sort at `islands.py:585` is untested;
  (2) `detected_frac` is not pinned (add a `missing = exclude` case with a detected NaN-abundance gene
  inside an island); (3) `_eligible` duplicates the `replicon_index` rule in `regions.py:83-86`;
  (4) gene order of a wrapping island's `genes` list (sorted by start, not window order) — decide and
  document; (5) `TuSeries.genes` docstring says genome order, the spec says transcription order;
  (6) the `min_tus` guard is untested for `min_tus ≥ 3`; (7) the B = 9 hand-computed p check is not
  literally present; (8) the PR was 3,062 lines against the ~400 guideline.
- M6-01 hardening (statistician, #140): `_KDE_MAX_GRID` coarsens silently; a bin with targets but no
  decoys gets weight 0 in f0 (fix before lfdr is ever thresholded); the lfdr can come from a different
  decoy set than `which` reports; three proposed attacks are in docs/methods/operons.md for M6-09.
