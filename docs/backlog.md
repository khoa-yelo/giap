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
- `tests/unit/test_setup_resources.py::test_contaminant_library_is_not_vendored` walks the whole
  checkout and counts `.nf-test/.../work/**/contaminants.fasta` (a symlink to the real
  $GIAP_DB_DIR library staged by mapping_nf_main_review.nf.test), so it fails after a local
  nf-test run in the same worktree. CI never sees it (separate jobs); add `.nf-test` to its skip
  set.
