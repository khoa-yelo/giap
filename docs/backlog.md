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
- M7-01: co-optimal alignments with equal cost, diff-type multiset and difference positions but different op shapes are broken by DP move order (1:1, 2:2, 1:2, 2:1), not pinned by a test; the frameshift property oracle sorts kinds alphabetically, not by DIFF_TYPES precedence (no divergence found in 200k random junction cases).
- M7-01: align_hits loops rows in Python (~110 µs/hit); batch the DP in numba if Module 5 ever aligns far more than 10^5 hits.
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
- M5-09 (lead decision 3): `context_stats` writes `tables/genes.tsv`, so a run that skips it (no
  C7/C8, as the "precomputed regions layer" nf-test does) produces no genes.tsv at all. Resolved
  when §6.2's `finalize` step arrives and takes ownership of the file; `gene_context.parquet` and
  the gene table are already on disk for it to join.
- M5-09: `giap.steps.registry` gained `optional_outputs` (a subset of `outputs` a caller may leave
  unwired) so the M5 keys could be added to `context_stats` without rewiring M4-18's registry-path
  callers. It has no dedicated unit test yet — the test-author should cover the decorator check,
  `_check_keys` and the "wired together or not at all" rule in `_m5_context`.
- `np.linalg.lstsq` in `giap.context.regions.residualise` is not bit-portable across CPUs (LAPACK SVD),
  so `regions_stats.tsv` float columns differ by 1–2 ulp between machines and between CI legs (M4-20
  decision 9, probe 1863133). Writer-level rounding is the real fix — see the §9.5 per-field precision
  item — and belongs with the M6 writers.
- tests/helpers/islands_data.py (m5/islands) and tests/helpers/oriter_data.py (m5/oriter) both carry a
  `tu_frame` and a giapsim-scenario driver; whichever merges second should deduplicate. Note the ori–ter
  driver deliberately keeps the simulator's real coordinates, while the island one re-lays genes through
  `Contig`, which would destroy the ori–ter geometry — so the shared helper must preserve both modes.
- M4-21 pass 2 (code review of #145): F2 unit tests for `precomputed_members` (constant overlap_frac,
  skipped unknown tag, core_biosynthetic mutants survive); F3 no real-path test of prophage subset
  classes beyond `cargo`; F4 log or `notes`-tag that precomputed regions have product-only subsets and
  `core_biosynthetic = False`; F5 unit-test `topo_order` with `members`; F6 `_precomputed` ignores
  `context.region_types` (the tool path filters it); F7 the sensitivity goldens (sim50k, sim100k_real)
  have zero primary/variant calls; F8 stale flat files win over keyed ones in `contract_paths` /
  `work_inputs`. `PUBLISHED_FILES` repeats the Nextflow names (pinned by a test) — one shared source.
- MVP-02 pass 2 (code review of #146): `-resume` reuses a stale region_layer result when only the
  precomputed TSV's content changes (stage it as a declared input, decision 2, which also replaces the
  absolute-path fix in `write_params`); example 4 ran under farmshare_local, the farmshare,apptainer
  conclusion (GIAP_STEP runs natively) is from reading the config, not from a run.
- RegulonDB TLS: the server omits its intermediate, so any plain client fails verification. GIAP completes
  the chain itself (fetch_regulondb.py pins GlobalSign GCC R46 OV TLS CA 2025, sha256
  d160e2de4e56cb10b66c1cb0adcb79cf93c78dcd9db30c2018220262c04063f9, valid to 2029-06-23, verified against
  the already-trusted GlobalSign Root R46). Re-check when the intermediate expires, or drop the workaround
  if RegulonDB ever publishes a complete chain.
- M6-11 pass 2 (non-blocking items from the code review of #141; lead decision 7):
  - fetch_regulondb: when the endpoint reports no release, fall back to the date + hash (ADR-0017)
    instead of the directory name `unknown`;
  - record `(date, sha256)` per retrieval, not only the dates, so a changed download is visible;
  - the unreachable-host error names the PROGRESS.md › Blocked entry, as the spec's degenerate row says;
  - clean up `.part` files on failure and fsync before the rename in `write_atomically`;
  - bound the intermediate download (read at most a few kB);
  - type-check the GraphQL payload fields (content is `str`, metadata `str | None`);
  - catch `UnicodeDecodeError` in the retry loop;
  - record the genome-version equivalence NC_000913.3 = U00096.3 (RegulonDB vs the recorded Bakta contig);
  - the tests live in `tests/unit/test_m6_11_regulondb_*.py`, not the spec's `validation/tests/`
    (pytest's testpaths is `tests`): amend the spec or move them;
  - the "no RegulonDB bytes committed" check matches file names only; add a content or size check.
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
- M6-01 pass 2 (code review of #140, non-blocking): (N1) replace `_COLLAPSE_FRAC` with the exact
  scale-free criterion: the Student-t likelihood is unbounded when the largest tie count k satisfies
  k > df_lo * (n - k); (N2) take the KDE grid span of `_reflected_density` from the sample only, not
  from the evaluation points; `_bins` silently truncates float bins; `bin_index` accepts NaN edges;
  the +inf distance in the `bin_index` docstring.
- M6-04 polish (non-behavioural): `promoter_between_tracks` exposes the private `_Track` alias in its signature (make a public `PromoterTrack` type); pin `pair_frame` dtypes (`detected_*` bool, `flags` list) through a pairs.tsv schema when M6-08 writes the table; consider one shared index object for the M2 `*_between` helpers so M6-07 reuses it.
- M8-07 dedupe: `giap.determinants.features._eligible` is a local copy of the §7.3 regions-table rule
  (`regions_table._eligible`); switch it to `giap.genetable.eligible` once M8-03 (#159) merges
  (lead decision 2 of docs/specs/M8-07-determinants-features.md).
- M8-07 review of #162 (approve; lead: `_lookup` deduplicates before indexing):
  - **M8-08 brief:** on sim100k, 14 encoded columns are constant (`start_codon_other`, `replicon_type_unknown`, 12 `nterm_residue_*`). The elastic net's standardisation must handle zero variance and NaN cells (impute or drop), and must drop any column that is all-NaN after encoding;
  - an `nterm_residue` outside the 20 levels (X, U, `*`, lowercase) gives all-zero one-hot rather than NaN; a lowercase `aa` skips the Met rule;
  - test gaps: the C2 `ribosomal_protein` column winning, and each fallback arm alone; one-hot NA → NaN; §6.5 tie-breaks (same start with a different end; strand or locus_tag ties);
  - methods: say that `layer_missing:<layer>` is inferred from which features are absent from the returned map.

- M6-02 (test-author): `tests/unit/genetable/test_contract.py` pins `TABLE_KEYS` (no `peptide_runs`)
  and `m4_contract_fixture.build_frames` has no C16 frame. `test_tables_keys` and
  `test_missing_column[peptide_runs-*]` fail against the C16 contract that
  `test_contract_peptide_runs.py` requires. Add C16 the way M5-04 added C15.
- M6-02 (lead to confirm): `fit_detection_limit` puts the peptide model in force when it is usable,
  else the replicate model when it is usable, else the unusable peptide model. The spec does not say
  what to do when the peptide model exists but is unusable while runs exist.
- M6-02: the replicate model takes R as the number of distinct C16 runs for every gene. The method
  doc's "R is the gene's own number of runs" needs a per-gene run-eligibility input, which C16 does
  not carry.
- M6-02: `detectability.covariates` absent from the `genes` frame are silently left out of the
  replicate model; M6-05H should pass the M4-10 covariate table joined in.
- M6-02: the method doc asks the detection summary to count detected genes whose evidence lies
  only outside the observable set (not training genes); not yet emitted.
- M6-02: `k_min = 2` whenever any gene carries `detection_tier == "two_peptide"` (the tier is
  run-wide in M4-12).
- M6-02 review of #157 (lead fixes: a non-finite abundance gives NaN in both models; a k_min = 2 recovery test; `_newton` requires a vanishing Newton step and a negative definite H at convergence):
  - the replicate model needs `cov=` at scoring, which the spec signature lacks. A NaN covariate scores NaN, but training imputes a median the model does not store; store it;
  - a stray run label inflates R (the M1-06 `""` run plus 2 real runs gives R = 3, usable). Reject or flag it in C16 validation;
  - a constant length or GRAVY gives a singular H, reported as `separated_or_not_converged`. Drop zero-variance features as the replicate model does;
  - NaN `tu_id` values merge into one fold group; the replicate folds depend on how much rng the peptide fit consumed; I/L variants inside one gene are deduplicated raw but matched folded; non-standard letters create unmatchable observable peptides;
  - test gaps (surviving mutants): I/L matching; the decision 8 in-force rule; the log of `length_aa` at scoring; the `one_class` guard; the MAX_ITER and line-search-failure paths; the convergence tolerance; the out-of-fold AUC.
- M6-02 statistician review (56d5a72), non-blocking: a stray `""` run inflates R (C16 validation); store the covariate median in the model; under `two_peptide` the replicate model counts a run with ≥ 1 peptide as seen; detected genes with imputed abundances train on fill values; the DL1 docstring still says the naive oracle is "recorded, not gated".

- M8-03 review of #159 (approve; lead: ADR-0023 adds map01310/map01320, methods note):
  - the GO closure does not filter by namespace (go-basic has no cross-namespace edges): filter in `up()` or assert in `read_go`;
  - `if kegg_names` tests a Path's truthiness; a missing names file gives all-NA names with no flag;
  - `build_collections` turns an integer index array into a wrong bool mask without error, and a wrong-length mask raises a bare IndexError; check shape and dtype;
  - `_ids` treats a plain string cell ("GO:1;GO:2") as empty;
  - test gaps: roots in `parents`; KEGG and COG `n_unknown_ids`; slim names that differ from the GO names.

- M8-13 (from M8-02): the random-ES kernel walks each prefix in O(K²) instead of fgsea's block
  convex hull (0.61 s for the 2,500 x 1,500 x 1,000 benchmark, so no need yet); its statEps
  weighting of zero stats is ported but unpinned, and R's long-double sums may decide near-ties
  differently — check against the golden env.
- M8-02: `fgsea_simple` raises an even `sample_size` by one (fgseaMultilevel does); GIAP's config is
  fixed at 101, so this only matters if the key becomes configurable.
- M8-02 review of #156 (approve; for M8-13 or pass 2):
  - test gaps (survived mutants; the code is correct): stable tie order with more than 16 tied values; sampling without replacement (an exact-distribution test with K ≥ 2 and small N); the modeFraction < 10 boundary; "NES mean = 0 gives NaN"; the asymmetric ES = 0 branches for nMoreExtreme and modeFraction; a `_needs_multilevel` table near its boundary; statEps with zero stats; even-sampleSize rounding;
  - M8-13 parity: the O(K²) prefix walk equals exact ES. Where fgsea's C++ hull differs, the cause is its exact integer tie arithmetic (port ±x against C++ 0, e.g. n = 9, S = [5, 7, 2, 3, 6, 8]) and a possible hull underestimate at small k, not only long double;
  - the single-set path: R's `calcGseaStat` takes no statEps when only one set survives, while the port always uses the statEps kernel (they differ only when there are zero stats);
  - `prepare` should reject empty-string labels, as `preparePathways` does;
  - put the full MIT text in the module headers and in `tests/helpers/fgsea_oracle.py`, not "(full text in NOTICE.md)";
  - the `benchmarks/baseline.json` `commit` field is free text, not a hash.
- M7-01 review of #154 (non-blocking):
  - untested: two same-direction 1↔2 blocks (|n−m| = 2, align.py:184); leftmost tie-break across frameshift junctions; the 2→1 block start (`pi = i - 2`); the move-3 bounds guard;
  - `_check` could also compare the recomputed diff-type multiset with the row's `diff_types`;
  - frameshift work-table lookup: filter to `frameshift` rows, and check blocks against the row's start/end;
  - `test_align_sim100k.py`: assert that frameshift rows occur;
  - test_align_props.py:266: sort kinds by RANK, not by name (unobservable today).
- M6-06 review of #148, owed to M6-08/M6-09:
  - ADR draft for Khoa: exclude pairs with an imputed gene from targets and decoys (statistician: a co-imputed NO surplus gives ~198 false calls per replicate on a 2,000-pair genome; exclusion restores the level);
  - attacks M-1 (co-imputation surplus; fails until the ADR), M-2 (class-dependent abundance/detection mix), M-3 (M6-01 attacks A and B);
  - validate `operons.pi0` in `OperonsConfig` so the degenerate paths reject it too;
  - `genes=None` means molar and all eligible: the step must always pass `genes` (or make it required);
  - drop the internal `which` column before `pairs.tsv` (M6-07 decision 4);
  - mutant test gaps: `partial` on `gene_a` only; `N_NO` counting only targets; `null_model` name with a dropped set; deterministic decoys-before-targets order; decoys' `partial` symmetry.
- M6-06 (pass 2): `bin_counts` and the π₀ actually used must reach the §9.4 run report and the manifest when the operons step is wired (M6-07/M6-08); the function returns them but nothing writes them yet.
- M6-06 (pass 2): the missed-operon benchmark baseline was measured in the iron-15 session allocation only (0.23 s); re-measure on wheat with the other M6 entries.
- M6-06 (pass 2): the lfdr of targets in a bin with no decoys is too small (methods note, "Known asymmetry"); set it NA or 1 there before lfdr is ever thresholded.
- M6-05 review of #147 (code reviewer, statistician), owed to M6-05H/M6-09:
  - skip the t fit when the scale is relative (today it fits and may warn "too_few_op_pairs" on relative rows);
  - log the 1e-6 scale floor in the non-convergence warning, not the raw MAD;
  - document that one-detected rows get no `no_detection_limit_model` on the too-few/relative paths;
  - test gaps (surviving mutants): undetected gene with a finite imputed value; `unestimated` on the too-few and relative paths; one side partial only; pair flags (`rna_between`) kept; flags of partial rows;
  - statistician F3/D-1: abundance-dependent Δ noise breaks the lowest-quartile level (2.4–3× nominal) — method ADR or narrowed claim, Khoa decides; D-2 (TU-chained Δ) and D-3 (power record) as attacks;
  - F5: keeping P(miss) out of the BH family in M6-05H needs an ADR; F6: flag imputed-but-detected tested pairs `abundance_imputed`.
- M6-05 pass 1 simplifications (implementer): `discordance()` accepts `cfg` and `p_miss` but uses
  neither — a NaN abundance is untested under every `abundance.missing` mode, not only `exclude`, and
  the one-detected branch is always `unestimated (no_detection_limit_model)`; `PMissModel` is a
  placeholder Protocol until M6-02 merges (replace with its real type in M6-05H).
- M8-03 review of #159 (approve; lead: ADR-0023 adds map01310/map01320, methods note):
  - the GO closure does not filter by namespace (go-basic has no cross-namespace edges): filter in `up()` or assert in `read_go`;
  - `if kegg_names` tests a Path's truthiness; a missing names file gives all-NA names with no flag;
  - `build_collections` turns an integer index array into a wrong bool mask without error, and a wrong-length mask raises a bare IndexError; check shape and dtype;
  - `_ids` treats a plain string cell ("GO:1;GO:2") as empty;
  - test gaps: roots in `parents`; KEGG and COG `n_unknown_ids`; slim names that differ from the GO names.


- M6-06 review of #148, owed to M6-08/M6-09:
  - ADR draft for Khoa: exclude pairs with an imputed gene from targets and decoys (statistician: a co-imputed NO surplus gives ~198 false calls per replicate on a 2,000-pair genome; exclusion restores the level);
  - attacks M-1 (co-imputation surplus; fails until the ADR), M-2 (class-dependent abundance/detection mix), M-3 (M6-01 attacks A and B);
  - validate `operons.pi0` in `OperonsConfig` so the degenerate paths reject it too;
  - `genes=None` means molar and all eligible: the step must always pass `genes` (or make it required);
  - drop the internal `which` column before `pairs.tsv` (M6-07 decision 4);
  - mutant test gaps: `partial` on `gene_a` only; `N_NO` counting only targets; `null_model` name with a dropped set; deterministic decoys-before-targets order; decoys' `partial` symmetry.
- M6-07 pass 1 simplifications (implementer):
  - `flags.py` imports the private `_bool`, `_flags`, `_per_gene` helpers from `missed.py` (also duplicated in `discordance.py`); move them to one shared `giap.operons._util` module;
  - `operons_table` restricts `n_flags` to flagged **OP** pairs with both genes in the TU, so `sum(n_flags) <= flagged OP pairs` holds by construction even if a C5 TU ever contained an NO pair; confirm with the lead or drop the restriction;
  - `detected_frac` for an all-partial TU is NA (0/0), as the spec lists it owed to pass 2;
  - `merge_tests` keeps the first row when a `pair_id` appears in both M6-05 and M6-06 rows (cannot happen: OP vs NO), with no warning;
  - `evidence_note` wording is placeholder prose (spec "Hardening owed"); `layer_missing:*` is not added to SW rows (pass 2);
  - `flag_pairs` builds notes with `iterrows` (fine at thousands of flagged pairs; vectorise if a benchmark says so).
- M6-06 (pass 2): `bin_counts` and the π₀ actually used must reach the §9.4 run report and the manifest when the operons step is wired (M6-07/M6-08); the function returns them but nothing writes them yet.
- M6-06 (pass 2): the missed-operon benchmark baseline was measured in the iron-15 session allocation only (0.23 s); re-measure on wheat with the other M6 entries.
- M6-06 (pass 2): the lfdr of targets in a bin with no decoys is too small (methods note, "Known asymmetry"); set it NA or 1 there before lfdr is ever thresholded.
- M6-07 review of #152 (non-blocking):
  - tests to add: the terminator cut-off on its own (promotech 0.5 / terminatornet 0.9 / score 0.8); a pair's own flags kept by `merge_tests`; the both-detected guard in `max_abs_delta_log2`; a gene missing from the gene table counted as undetected; the fully linked circular ring (wrap OP pair inside one TU);
  - validate `missing_layers=` ∈ {promoters, terminators}; guard `merge_tests` against empty inputs that have no columns;
  - move `_bool`/`_flags`/`_per_gene` to `giap.operons._util`;
  - M6-10: restricting `n_flags` to OP pairs matters once pairs are reclassified at UniOP 0.3/0.7 while the TUs stay at the default.

## M8-05 review nits (PR #171)
- Surviving mutants (test gaps):
  - re-tests on a fresh rng;
  - unestimated rows set to p = 1 in the substituted BH vector;
  - boundaries: `k <= min_size`, `bh(sub) >= alpha`, `q < alpha`;
  - the `parent_group` ID tie-break.
- A missing `tu_id` crashes with a TypeError. Fail closed with a ValueError, as context does.
- A NaN `tu_position` makes gene order row-order dependent.
- An infinite `abundance_log2` passes `notna()`.
- Duplicate `locus_tag` rows are counted twice.
- A missing collection name becomes the literal "NA".
- M8-15 should adopt the statistician's GO-family FDR arm (local/stat/m805code/probe_code.py `go_null`/`go_alt`), and cite the R = 4,000 level runs (Slurm 1867122).
