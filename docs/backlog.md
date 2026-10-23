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
- M2-15: test the manifest `uncalibrated` flag and step-level `cross_check_transtermhp: false` once manifest and optional-output support exist.
- M2-15: parse_bactermfinder's input format is provisional (no recording); record a real BacTermFinder run (opt-in tool, env not built) and pin the parser to it.
- M4-20: replace M4-16's test stand-ins (p_detect_true, fixed COG table, tier/detectability columns) with the real M4-12/M4-13 outputs once merged.
- M4-14 (statistician): re-measure abundance.ibaq_denominator on M4-09 data before switching (fails 10–40 % at p_obs = 0.3); restate as a pooled mean Δρ if needed; add the p_obs = 0.3 attack.
- M4 notes: median imputation pulls detected-but-unquantified genes toward the centre; the min/exclude sensitivity reruns bracket it (statistician, #82).
- M7: compute the M4-16 QC metric id_mixups (§7.5.1: same-database ALT sites that are not start/stop sites); M4-16 reports NA + requires_variants until then.
- M7-01: co-optimal alignments with equal cost, diff-type multiset and difference positions but different op shapes are broken by DP move order (1:1, 2:2, 1:2, 2:1), not pinned by a test; the frameshift property oracle sorts kinds alphabetically, not by DIFF_TYPES precedence (no divergence found in 200k random junction cases).
- M7-01: align_hits loops rows in Python (~110 µs/hit); batch the DP in numba if Module 5 ever aligns far more than 10^5 hits.
- M7-02: frameshift site `alt_aa` is the placeholder `fs+1`/`fs-1` and `ref_codon` the last prefix codon; settle both in M7-03 with the VCF writer (lead decision 3). The M7-03 writer must drop the internal `_ref_codon`/`_alts` pileup columns; pandas `attrs` (frameshifts, contig order) do not survive `pd.concat`.
- M7-02: a folded peptide with several categories at one codon (only possible with conflicting alignments of one sequence) is counted as its smallest ALT key, not pinned by a test (the oracle rejects such cases).
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
  *Khoa (2026-10-22): flag it.* The `tier` row `q_floor` = 1/T and the warning `q_floor_above_alpha`;
  q and alpha unchanged (fix/q-floor-flag).
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
- M7-06 review of #161 (lead fixes: an antisense `overlapping_gene` from the opposite strand; `cds_outframe` takes the smallest same-strand tag; I/L twins share the contaminant ban and the multi-locus count):
  - M7-09 must attribute shared Promotech windows by `(contig, strand, end3p)` against each ORF's starts, not by `owner` (first-window-wins);
  - the M7-15 step must pass the in-memory 0-based C11 frame, not the published 1-based `peptide_hits.tsv`;
  - frameshift (tolerant) hits walk downstream in the prefix frame; exclude them or record it;
  - `nterm_extension` assumes the next canonical stop is the gene's own stop (selenoproteins, recoding, partial genes);
  - N bases are not stops (undocumented); the whole-turn error names no frame number;
  - test gaps: the `#` library separator; library-side I/L folding; a non-11 translation table; `frame` across an origin with L % 3 ≠ 0; the window dedupe and sort key; `window_chunks`.


## M7-07 review nits (PR #165), for pass 2 (M7-15)
- **B2:** the − strand no-start 3'-open heuristic in `_span` (`outputs.py`). Add an `open_3p` flag or column after an ADR on the §9.4 vocabulary. Lead decision 4.
- The translation table is hard-coded to 11 in `outputs.py`, while `reconstruct` uses `cfg.annotation.translation_table`. Assert it, or thread the table through.
- Flags and `class_fdr` strings copied from `scored` are not checked against the schema vocabulary at write time.
- `start_best_ostir` means the codon's lowest coordinate in `scored` but the transcription-direction first base in the table. Document this for M7-08.
- Make `_strand_view`, `_walk_up`, `_stops` and `_il` a public helper in candidates.
- The sort key uses the transcription-direction `end` rather than the §6.5 highest coordinate, which only matters for strand ties.
- Rows in `scored` with no matching ORF are silently ignored.
- The `nterm_extension` patch does not check that the gene's strand matches the ORF's.
- Surviving mutants (test gaps) — pin each with an explicit example:
  - the `start_best_ostir`-not-a-start ValueError;
  - `ostir_missing` on no-start ORFs;
  - `partial` on 5'-open ORFs;
  - start and stop codons straddling the origin;
  - `called` when `scored` is None;
  - the strand tie-break;
  - the patch wrap across the origin;
  - flags carried from `scored`;
  - `n_peptides` taken from `scored`.

## M7-08 review (PR #168, approved): test gaps and nits
- Surviving mutants (correct code, untested):
  - the TSV promoter `+1` (use a single scored window, and assert NA one position away);
  - a minus-strand start codon crossing the origin (`first_base` wrap);
  - a multi-start ORF whose best-OSTIR and longest starts carry different promoter scores;
  - a NaN `tir_ostir` row;
  - the `novel_rbs_layer` count-mismatch guard.
- The empty second round under Nextflow logs a spurious `layer_missing:rbs_novel` WARNING.
- PROVENANCE.md cites the unreachable WIP commit 336cc1b. Cite 4760059 instead.
- `annotation.ostir.downstream_nt` has no lower bound. Below 3 it cuts the last start codon, in the first round too.
- Paths only a real run exercises (the header-only OSTIR CSV, an empty `novel_windows/`) are covered only by the owed example 7.
- The simulator's `rbs_novel` rows leave out codons crossing the origin. Note this in the M7-12 scenario.
- M8-01 (pass 2): `wallenius_tails` evaluates scipy's Wallenius pmf over the whole support on every
  call (both tails are direct sums). M8-04 calls it once per category with the same (M, N) and many
  (n, ω), so a per-(M, n, N, ω) cache is useless; if the 3,000-category budget tightens, sum only
  the shorter tail plus a pmf evaluation of the other near k, or move to the log-space pmf the
  far-tail xfails need anyway (M8-13 ADR). Measured 2.30 s / 5 s on iron-15 (EPYC 7742).
- M8-01: `elim` raises on a parent named in `parents` but absent from `term_genes` and gives
  (NaN, NaN) without calling the test for a term with an empty gene set; neither is pinned by a
  test yet (test-author). `wallenius_tails` raises on NaN/negative odds and on n or N outside
  [0, M]; also unpinned.
- M8-01 review of #153 (non-blocking):
  - pin the untested fail-closed contracts: `wallenius_tails` NaN/inf/negative odds and n/N outside [0, M]; `elim` unknown parent/key; an empty gene set gives (NaN, NaN) without a test call;
  - `category_odds`: test the finite-p filter and the shape check; range-check p_detect ∈ [0, 1];
  - test ω = 0 with lo > 0, and N < n with k in (N, n];
  - non-integer k raises a TypeError deep in the slice; validate it;
  - M8-04: catch `ValueError` from extreme odds and report the category `unestimated` (`odds_out_of_range`); M8-13: exact pmf for extreme odds.
- M8-07 dedupe: `giap.determinants.features._eligible` is a local copy of the §7.3 regions-table rule
  (`regions_table._eligible`); switch it to `giap.genetable.eligible` once M8-03 (#159) merges
  (lead decision 2 of docs/specs/M8-07-determinants-features.md).
- M8-07 review of #162 (approve; lead: `_lookup` deduplicates before indexing):
  - **M8-08 brief:** on sim100k, 14 encoded columns are constant (`start_codon_other`, `replicon_type_unknown`, 12 `nterm_residue_*`). The elastic net's standardisation must handle zero variance and NaN cells (impute or drop), and must drop any column that is all-NaN after encoding;
  - an `nterm_residue` outside the 20 levels (X, U, `*`, lowercase) gives all-zero one-hot rather than NaN; a lowercase `aa` skips the Met rule;
  - test gaps: the C2 `ribosomal_protein` column winning, and each fallback arm alone; one-hot NA → NaN; §6.5 tie-breaks (same start with a different end; strand or locus_tag ties);
  - methods: say that `layer_missing:<layer>` is inferred from which features are absent from the returned map.


## M8-08 (determinants models), owed
- The acceptance benchmark, 2,500 rows × 35 columns in under 120 s, recorded in `benchmarks/baseline.json`.
- The `docs/methods/determinants.md` section "Models".
- The calibration tests need the protected `tolerances.yaml` key `determinants_delta_r2_max_abs_bias`, which is owed to Khoa (decision 6).
- The `too_few_rows` rows carry `n` = NA. Pass the model-row count to `model_rows_table`.

## M8-08 review nits (PR #166)
- Make `tu=` a required keyword in `cross_validate`, or assert in the M8-09 and M8-10 tests that it is passed.
- `cross_validate` has no input validation: the lengths of y, folds and tu are never checked against X.
- Document that `X.index` must be `locus_tag`, because a RangeIndex breaks canonical order.
- A fold with no usable columns appends no booster, so M8-09's repeat × fold indexing shifts.
- `cv_r2_sd` with R = 1 is NaN and carries no flag.
- Untested parameters: LightGBM `deterministic`, `force_row_wise` and the seed derivation; the fallback that predicts the training mean when no columns are usable.
- Constant features give AUROC 0.5 with CI [0.5, 0.5]. Consider a `constant_feature` flag.
- An NA `tu_id` in `univariate` becomes one shared pseudo-TU. Raise instead.
- `univariate` bootstraps even when the point estimate is NaN.
- The methods note must cover:
  - the inner-CV preprocessing;
  - the ΔR² caveat (the repeats share data);
  - univariate CIs are unadjusted and make no discovery claim.
- §7.8.2's separate report of ribosomal genes: confirm which issue owns it (M8-09 adds rows flagged `cai_reference`).

## M8-09 review nits (#173, one-round MVP review)
- Tests at rule edges: `>=` top-decile cut, SD ddof 1, `>` at 2 SD, ribosomal genes imputed/undetected with finite abundance.
- Untested error paths: OOF size vs model rows, booster count ≠ R·K, `len(X)`.
- `null_p` should check `X`/`y`/`tu`/`folds_rep1` lengths.
- `b_count` derived from p; `n_perm` is B not the non-NaN count — compute both directly.
- A NaN observed R² gives p NaN unflagged; add `unestimated`.
- Pin the LightGBM refit seed with a bagging-on test.
- Decide and test `top_decile` under `too_few_rows`.
- `_canonical` recomputed per permutation and model.
- Owed (acceptance): B = 1,000 benchmark in baseline.json; methods "Null and attribution" section.
- M8-10 carry-over (statistician #173): `cv_r2_mean` (the null p's effect) gets `ci_unestimated` (ADR-0031 pattern); `n_perm` should count non-NaN null values as `perm_p` does; optional exact test against the observed repeat-1 R².
- Hardening (ADR-0034): spatially contiguous CV folds for observed and null, a positional-feature decision (`replichore_position`), a detection-selection model; `null_uncalibrated` comes off only when every arm of test_attack_determinants_null.py passes.
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
- M7-02 review of #158 (approve; lead fix: `sites` raises when the pileup's attrs are missing):
  - **decide with M7-05:** I/L splits the swap and block ALT keys (`IK>KI` against `IK>KL` give two `candidate` sites where one `high` ALT is expected). Fold I/L, or write `J`, in swap and block tokens. This hurts attack 6 / M7-13 sensitivity. It also produces the multi-category-per-folded-peptide case, so the earlier backlog note calling that "only conflicting alignments" is wrong;
  - a NaN `placement_q` is silently eligible (`NaN > alpha` is False). Write an explicit rule (reject NaN unless the peptide is exact/`same_db`, as M4-11 does) and add a test;
  - M7-03: replace the `_ref_codon`/`_alts` columns and the `attrs` channel with a frozen `Pileup(table, alts, frameshifts, contig_order)` or an explicit `frameshifts=` argument;
  - `pileup()` does not apply `eligible`; the M7-04 step must call it first, or the docstring must say so loudly;
  - bare KeyErrors on a contig missing from the genome or an unknown strand; an empty `source` when support is missing; the −L gene shift is dead code;
  - under `same`/`other`, a stripped frameshift peptide still adds REF at its off-frame suffix codons (M7-05);
  - pin the str/object string dtypes in the M7-03 schema;
  - test gaps (surviving mutants): `il` ops kept under same/other; frameshift next codon on − strand (731/734 without 728); folded frameshift `n_alt`; `alt_aa` tie-break order; Int64 dtypes on header-only frames; a block peptide spanning the frameshift codons counting in `n_ref`.
- M7-03 review of #160 (lead fix: the datapackage resource profiles, plus a package-validation test):
  - replace the pileup's `_alts` column and the `attrs` channel with an explicit structure (a frozen `Pileup`); the M7-02 review test must change first;
  - `read_variants` is a byte inverse, not a frame inverse: `class_error_rate` is stored as a `%.6g` string and comes back as a float, and text dtypes differ. Keep the float in the table and format only in the writers;
  - §9.1 puts the descriptor at `OUT/schemas/datapackage.json`, so the `tables/...` paths resolve under `schemas/`. M7-04/finalize must settle where the descriptor and schemas are copied;
  - a bare `KeyError` when `sites` and `pile` mismatch; a `;`-joined string `flags` cell is exploded into characters; `read_variants` validates nothing;
  - weak schema patterns: `class_error_rate`, and `alt_aa` for `fs±1`;
  - for − two-codon events, TSV `codon_start` (the first codon in peptide order) differs from VCF POS (the lower codon); note it in the methods. A codon with both a substitution and a swap gets one GT per length-class record; document it;
  - test gaps: a − two-codon event whose second codon wraps past 0; duplicate or unsorted input flags; `read_variants` flag splitting.

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

## M7-04 review nits (PR #164)
- `OUT/schemas/datapackage.json` names `tables/...`, which Frictionless resolves under `OUT/schemas/`, not `OUT`. Lead decision 5 already owes this; the test checks only `OUT / path`.
- The table-schema copies into `OUT/schemas/tables/` (`runner/local.py`) are untested.
- `sample_name`: a peptide file named `auto.tsv` resolves to `auto`, and `write_vcf` then fails with exit 4. An explicit empty `sample_name` gives an empty VCF sample column. Guard both.
- Replicon order in variants output is untested because the fixture has one contig. Add a two-replicon case.
- `placement_q` is untested at the step level. Add a case with q above `stats.alpha.variants`.
- The stub VCF `#CHROM` line hard-codes `sample` (`variants.nf:37`).
- A single quote in a stub header would break the `printf` quoting. This predates M7-04.
- Spec drift:
  - the nf-test path is `workflow/tests/variants.nf.test`;
  - the process alias is `VARIANTS_STEP`.
- `read_gene_model` features are untested by the loci tests (mutant M22).
- Still owed:
  - the manifest `no_peptides` entry;
  - the real sim100k `giap_run.sbatch` run;
  - the step benchmark.

## M7-05 review nits (#175, one-round MVP review)
- Test a minus-strand CDS crossing the origin for `start_codon` (the `first %= rep.length` mutant survived).
- M7-16: decide `stop_codon` and `mass_equivalent` for blocks and swaps whose REF contains `*`.
- The §9.4 flag order is defined twice (`flags.FLAG_ORDER` and `vcf_records._FLAG_ORDER`); derive one from the other.
- `site_flags` rebuilds `start_codons` per call; point bulk callers to `flags_for` in its docstring.
- The FLAGS INFO description in vcf_records.py and variants.nf has no sync test (only the IDs are compared).
- MVP real run on 813 genes (Slurm 1865149, validation): operon tests now run (515 t, 102 decoy), but:
  - detection is saturated: 804/813 genes are detected, `p_detect` is the constant 0.989 fallback, and picked FDR q = 1/804 for all. The giapsim real-run scenario needs decoys with null scores (the earlier item) and lower peptide coverage;
  - Promoter Calculator runs one Nextflow task per TU (288 tasks for 0.9 Mb); batch it;
  - `--annotation bakta.json` now feeds `annotation.precomputed.bakta_json` (mvp/annotation-json). The previous runs re-ran Bakta, so `engine_protein_mismatch_rate` was 0.9997; re-run after the merge to confirm it falls.


- M6-06 review of #148, owed to M6-08/M6-09:
  - ADR draft for Khoa: exclude pairs with an imputed gene from targets and decoys (statistician: a co-imputed NO surplus gives ~198 false calls per replicate on a 2,000-pair genome; exclusion restores the level);
  - attacks M-1 (co-imputation surplus; fails until the ADR), M-2 (class-dependent abundance/detection mix), M-3 (M6-01 attacks A and B);
  - validate `operons.pi0` in `OperonsConfig` so the degenerate paths reject it too;
  - `genes=None` means molar and all eligible: the operons step always passes `genes` (M6-08); make it required in pass 2;
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
- M6-06 (pass 2): `bin_counts` and π₀ now reach `logs/events.jsonl` as the §6.7 event `operons.missed_operons` (M6-08); the §9.4 run report and `run_manifest.json` still have no writer that reads them. The first bin's `lower_bp` is written as `-Infinity` (Python `json` default), which is not strict JSON.
- M6-06 (pass 2): `bin_counts` and the π₀ actually used must reach the §9.4 run report and the manifest when the operons step is wired (M6-07/M6-08); the function returns them but nothing writes them yet.
- M6-06 (pass 2): the missed-operon benchmark baseline was measured in the iron-15 session allocation only (0.23 s); re-measure on wheat with the other M6 entries.
- M6-06 (pass 2): the lfdr of targets in a bin with no decoys is too small (methods note, "Known asymmetry"); set it NA or 1 there before lfdr is ever thresholded.
- M6-07 review of #152 (non-blocking):
  - tests to add: the terminator cut-off on its own (promotech 0.5 / terminatornet 0.9 / score 0.8); a pair's own flags kept by `merge_tests`; the both-detected guard in `max_abs_delta_log2`; a gene missing from the gene table counted as undetected; the fully linked circular ring (wrap OP pair inside one TU);
  - validate `missing_layers=` ∈ {promoters, terminators}; guard `merge_tests` against empty inputs that have no columns;
  - move `_bool`/`_flags`/`_per_gene` to `giap.operons._util`;
  - M6-10: restricting `n_flags` to OP pairs matters once pairs are reclassified at UniOP 0.3/0.7 while the TUs stay at the default.
- M6-08 pass 1 (implementer):
  - merge time: add `pairs`, `operons`, `operon_flags` → `tables/*.tsv` to `giap.runner.local.PUBLISHED_FILES` (#145; `test_local_publish_names.py` and the slow M6-08 local tests need them);
  - `run_operons` returns a 4-tuple (pairs, operons, operon_flags, bin_counts), not the spec's 3-tuple, so the event needs no second missed_operons call;
  - `schemas/datapackage.json` does not exist in the repo; the three resources are per-table files like the others (spec acceptance item);
  - `benchmarks/baseline.json` runtime of the step on a 4,500-gene genome is not recorded yet (spec acceptance item);
  - the publishDir rule (only `tables/*.tsv` in a subworkflow's out spec) is written in `operons.nf`; `genetable/genes_pre.parquet` from MVP-01 still publishes;
  - `peptides`, `peptide_hits`, `peptide_runs`, `detection_summary` are declared but not staged by OPERONS nor read (M6-02, pass 2);
  - the `giap.steps.common` importers alias `one_path as _path` etc.; rename the call sites when the modules are next touched.
- M6-08 review of #155 (approve; pass 2):
  - **must-do before the manifest writer:** open bin bounds `±Infinity` in the `operons.missed_operons` event are not strict JSON; write them as `null` in `_plain`;
  - local executor stages `peptides` and `detection_summary` into `operons` (and orders it after `context_stats`); OPERONS under Nextflow stages neither; reconcile before M6-02 reads them;
  - missing-layer WARNINGs have no §6.7 `operons.layer_missing` event;
  - `assign_tus` on precomputed `tx_rates` in `_layer` is dead work (`pair_frame` never reads `tu_id`);
  - `_cell` stringifies `None` items inside lists, bypassing `join_list`'s missing-item guard;
  - nf-test lives at `workflow/tests/operons.nf.test`, not `workflow/tests/subworkflows/`;
  - test gaps (surviving mutants): `missed_operons` without `genes=` (needs a relative-scale/partial-gene fixture); precomputed branch of `_layer` (promoters → `promoter_windows`) has no unit or nf-test; `no_detection_limit_model` log level unpinned.
- MVP real run (Slurm 1864792) is degenerate by input, not by code (validation diagnosis):
  - 82 candidate genes, so the specified picked FDR floor (0 + 1)/82 = 0.0122 > α = 0.01 (§7.3.2, GIAP_SPEC.md:1331, worked example at :3253), and no gene can be detected;
  - then `p_detect` = 0, the detectability constant fallback; `detected_2pep` and `abundance_log2` are independent of the call;
  - giapsim (tooling): the real-run scenario needs several hundred genes, and decoy rows are written with NA q/pep/score so a decoy can never win; draw them from a null;
  - M4: warn when the number of eligible genes is < 1/α; a new flag (e.g. `picked_floor_above_alpha`) is new vocabulary and needs Khoa/ADR;
  - M0 preflight `evidence.py:286`: the ranking-column check passes on target q-values; check the decoy rows' values.

## M8-04a review nits (#172, one-round MVP review)
- N2: the `lo > hi` branch is indistinguishable from the k-range check; pin it by error message.
- N3: Newton vs pure bisection in `_solve_theta` is untested (speed only).
- N4: `_solve_theta` returns at the iteration cap without checking the 1e-9·total tolerance.
- N6: `tilt` is typed `total: float`; the spec signature says `int`.
- N7: the benchmark docstring quotes 0.4 ms/call; recorded 0.66–0.76 ms (within budget).
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

## M8-05b review nits (#174, one-round MVP review)
- The detection |z| rule's "ES NaN" branch never fires on real rows (fgsea_simple returns finite ES for unestimated rows): key it on `unestimated` (p NaN) or reword spec 9b / methods to "ES = 0" (statistician asked).
- `cluster_deltas` rebuilds contig map, positions and P_w per collection; build once per `gsea_table` call.
- A NaN `contig` becomes a linear pseudo-replicon "nan"; guard or comment.
- M8-06 checklist: register `genomic_cluster` and `detection_confounded` in the `enrichment_gsea` schema (§9.4).
- `test_gsea_table_flags_need_no_rng` uses the fake engine; add a real-engine p-value comparison to pin the RNG stream.

## M5-10 reviews (#184; code approved, statistician BLOCK on the note, B1–B3 fixed in the PR), nits
- Code: `islands_family` copies `cfg` without re-validation (unlike `with_overrides`); re-validate or comment why it is safe.
- Tests: an island wrapping the origin of a circular replicon; pin the replicon order of draws (two replicons near threshold); a direct `islands_family` test with `stats.permutations != 199` and an automatic block length (mutants: topology flip, length + 1, reversed replicon order survived; $GIAP_ROOT/review/pr184-mut/).
- `islands_frame` computes 199-resample CIs the flagged set never uses; skip them if island counts grow.
- No golden has a non-empty islands set; add one with a planted island once M5 goldens exist.
- Stats: a Monte Carlo baseline row (the primary rerun on an independent stream) as the reference for variant rows (§8.6 / lead decision).
- Stats: one random substream per replicon in M5-06, so a T change elsewhere leaves a replicon's draws alone (changes islands.tsv; ADR-0021 note).
- Stats: `alpha.context` < 1/200 makes every islands set empty and every Jaccard NA; warn or document.
## M5-07 reviews (#188; code CHANGES B1–B2 and statistician LIFT, fixed in the PR), nits
- **Owed before the M5 gate is ticked:** run the `islands` array task of `scripts/slurm/calibration_full.sbatch` on main, then fill "Main run: Slurm ⟨ID⟩" in docs/methods/islands.md and PROGRESS.md.
- `test_planted_island_avoids_regions` checks nothing (empty region list); delete it now that `test_power_layout_has_no_regions` pins `regions == ()`, and amend the spec edge case "island_power guarantees no overlap" (false: overlap is left to the caller).
- Arm (a)'s docstring claims ℓ′ is "capped at 20 (asserted in every replicate)"; the code asserts ≤ 20. The attack arms should log ℓ̂ / ℓ′ through `_log_block_lengths`.
- The recorder should write per cell: hit count, Wilson CI, q, median 95th percentile of M_b, `low_direction_ceiling(q, n_tus)`, and the null FWER (lead decision 30).
- Replace the deepcopy replay in `_replicate` by a helper that returns the truth with the contigs (the replay is now asserted against `scenario_contigs`).
- The spec names `test_attack_islands.py`; the arms live in `test_attack_islands_scenarios.py`. The recorded ρ 0.9 cells overlap `test_recorded_autocorrelation_rho09_sizes`; arm (b) is the existing gradient arm at 4,500 genes. Move `replicate` / `gate_pair` into tests/helpers/islands_data.py (N6).
- `test_low_direction_is_a_structural_zero` (50 s locally, near the 60 s budget) should assert ceiling < M_b 95th percentile per replicate.
- `islands_recorded.jsonl` lines lack `stream` / `n_genes`.
- Optional full-tier teeth test: a fixed block of 3 at ρ 0.9 (0.765 at R = 200 in the probe).
- The effect-CI coverage the islands note promised for M5-07 was never measured.
- Recorder CLI: NaN → None mapping untested; `--replicates 0` and bad `--cells` need argparse validation; `--run-id` should default to `$SLURM_JOB_ID` and record a dirty commit; append (not prepend) to sys.path or only under `__main__`; `test_arms_are_the_spec_scenarios` inherits `calibration_full` though it takes milliseconds.
- Power-grid sbatch wrapper (36 cells, xargs -P 36) belongs in scripts/slurm/.
## M8-06a review nits (#181; one round. B1 and N1 fixed in the PR)
- **N2** `n_unknown_ids` is never non-zero in any fixture (mutant M09). Add a gene with an unknown GO ID or KO.
- **N3** `ora_size` and `gsea_size` are tested only at their defaults (M11, M17). Add a run with a non-default size.
- **N4** The function events are logged before the writers run. If a writer fails, `events.jsonl` still holds events for a step that exited 4.
- **N5** benchmarks/baseline.json records the commit as "0b8ee09 plus uncommitted". Re-record it against a real commit.
- **N6** `_natural_key` gives `TU_1` and `TU_01` the same key. Add a string tiebreak.
- Mutation evidence: Slurm 1867940, $GIAP_ROOT/review/pr181-mut/.

## M7-09 review nits (#180; one round, blockers were test-only)
- **N1 Homopolymer window edges are untested** (mutants m26, m27): add runs at [591, 596) and [607, 612). Leftward run extension (m64) and the linear contig end (m65) are untested too.
- **N2 Untested equalities:** `terminator_downstream` at conf == `terminatornet_yes` (m42), and BacTermFinder skipping (m17).
- **N3 `FLAG_ORDER`:**
  - the test checks the helper's copy rather than `outputs.FLAG_ORDER` itself;
  - no test writes `partial` together with `ostir_missing` (m55);
  - unknown-flags-last is untested (m53).
- **N4** `sorf` rows sent through the pseudogene junction logic go undetected (m61).
- **N5** `read_peptide_intervals` on a pre-M7-09 parquet returns `{}` silently; it should log a warning.
- **N6** The `orf_id` rank tie-break is lexicographic (`_1000_` < `_999_`). Document it.
- **N7** `has_start` in the intergenic rule is equivalent in the pipeline (m08); the test can't detect its removal.
- **N8** For a frameshift pseudogene (length not a multiple of 3), the junction scan skips the last whole codon. Document it.
- Mutation evidence: Slurm 1867907, $GIAP_ROOT/review/pr180-mut/.
## M8-04b reviews (#179: code review CHANGES, statistician BLOCK, fixed in decision 13), pass 2 and nits
- **Pass 2: lift `ora_uncalibrated`.**
  - The statistician's flexible cross-fitted recalibration passes every cell of all 8 worlds at R = 3000 ($GIAP_ROOT/stat/m804b/probe_flex.py). Its terms are 4 quantile hinges in logit q_raw, a one-hot TU size, a logit × log n interaction, and ridge 1e-4.
  - It needs an ADR-0035 amendment (Khoa).
  - The flag can be lifted only on a pass at R = 1000 in both blocks of 4b–4f, and the gate needs rate logging.
- **N1** Comp still underflows to 0 for a very large or heavily clipped TU (`np.exp(log_comp)`, `expit(-eta)` once eta > 745). Carry log_comp into the tilt.
- **N2** Duplicate locus_tag rows are accepted: the gene is counted twice in n_genes.
- **N3** Clip semantics: `p_detect_clipped` goes only on rows whose full set holds the clipped TU, and the clip also moves p̂ in (0, 1e-12). State which rows count as "affected".
- **N4** ora.py imports gsea's private `_closure`, `_depth` and `_natural_key`, and `_closure` duplicates `elim._ancestors`. Move them to a shared module.
- **N5** Untested boundaries:
  - `MIN_CLASS` at 1 and 200, and `<` vs `<=`;
  - the inclusive size upper bound;
  - `q_value <= alpha`;
  - the parent_group ID tie-break;
  - the strict `<` in the elim replay;
  - `CI_MIN_UNITS` at 9 and `<=`;
  - `_PAIR_ATOL`;
  - comp computed as 1 − q on the raw and fitted paths;
  - `folds={}` skips the fit with `tu_recal_skipped`, as documented; the statistician's re-check corrected the earlier claim that it raises.
- **N6** `_classic`: a NaN `category_odds` gives `p_classic` NA with no flag, and the bare `except ValueError` can hide real errors.
- **N7** The end-to-end `detection_enrichment` benchmark is still owed; the current one covers only the primitive.
- **N8** test_attack_ora_tusize.py:96 recovers the replicate index from `seed_seq.entropy`, which is fragile.
- **N9** `_collection_rows` rebuilds the gene maps for every collection, and GO `_closure` is built twice.
- **Lifting the flag:** the commit that lifts `ora_uncalibrated` must add it explicitly on the skip path; test_ora_review_179.py pins this.
- **N10** The attack arms' `ora_uncalibrated` fallback cannot trigger on its own; the flag is a hand-edited constant.
- **Statistician:** category_summary's bootstrap undercovers small categories (10 TUs: 0.886; 15 TUs: 0.90). Raise `CI_MIN_UNITS` to about 30, or use BCa or a studentized bootstrap.
- **Statistician:** cross-fitting does not remove leakage when a category is a large share of one TU-size group. In the tested worlds this was only conservative.
- **Flaky?** `tests/unit/function/test_gsea_flags.py::test_detection_flag_with_the_real_engine_d1_d2` failed only in the 16-worker fast suite (Slurm 1867909 and 1867927, es = −0.5 for top=True). It passes alone and in its file. Check main and look for state leaking between tests.

## M8-06b review nits (#187; one round, approved)
- Groovy `installed()` uses `path.exists()`, Python `is_file()`: a directory at a resource path is staged by Nextflow only.
- `installed()` re-parses `db_manifest.json` per call (five more callers); parse once.
- The local executor stages function resources (up to 32 MB) even when the step cannot run (no peptides); Nextflow does not.
- `genetable.nf` header comment line is ~150 characters.
- No test pins that `tu_folds.parquet` is published under `--modules core` by both executors.
- Example 1 builds GO and KEGG xrefs but asserts only `cog` rows; assert kegg / go_slim too.
- Mutation evidence: $GIAP_ROOT/review/pr187-mut/ (5/5 killed).
## M8-12 review nits (#192; one round, B1–B3 fixed in the PR)
- `check_golden_provenance.py` could check that the lock sha agrees across PROVENANCE, VERSIONS and the lock; add a fixture where the base branch moved.
- Guard bypasses: the Rscript guard is case-sensitive and misses `R CMD`, `["R", "--vanilla"]` and concatenation; `check_imports` misses `importlib.import_module("rpy2")` / `__import__`.
- golden-r.yml: pin actions by SHA, move `contents: write` to the job, add the R / package version-pin check the sbatch does.
- Add an input set with exactly 10 names matched after absent names drop (sharper `minSize`).
- Precision test: per-field `repr(float(s))` round trip.
- VERSIONS.txt carries host-specific absolute paths (BLAS); the FarmShare and workflow routes always differ there.
- Open with Khoa: R6 (bioconda post-link downloads from bioconductor.org) and the PR-size gate counting `*.explicit.txt` (ci/gates.yaml is protected).
## M9 design follow-ups (lead, 2026-10-22)
- MVP-03: build the §6.2 / §6.7 `manifest` and `finalize` steps (run_manifest.json; KEGG retrieval date, provenance), before M9-16.
- M4: should `detection_summary` carry a region-**type** O − E row (§10.2 V7 lists region types; no table has them)?
- UpSetPlot 0.9.0 fails on pandas 3.0.6 (`Invalid RGBA argument: nan`); V16 uses the custom matplotlib fallback (§10.2).
- M8-17 (ADR-0038 item 1): BioContainers image for the fgsea env under apptainer profiles — deferred; conda-only for now (lead decision 9).

## M9-01 (style and plot framework), owed
- `RunView.table` leaves list columns (`genes`, `flags`, …) as `;`-joined strings (lead decision 2, 2026-10-22). Split them on `;` (ADR-0006) once the schemas mark which string columns are lists; until then each plot splits what it draws.

## M9-01 review nits (#191; one round, B1–B4 fixed in the PR)
- `test_figures_are_closed_after_writing` and the Agg-backend test are vacuous (Figure.number is always None; Agg passes without `matplotlib.use`).
- Origin-crossing BED features come back as `_part1` / `_part2`, not rejoined to end > length.
- `tus.tsv` / `regions.tsv` are read untyped beyond the coordinate columns; a header-only file gives `object` columns.
- The `work/resolved_config.yaml` fallback applies only in `from_outdir`.
- `COG_SUPER` grouping: cite the NCBI source of the super-category grouping (fun-20 group order).
- `normalise_ids("")` returns `[]`: decide whether an empty `--plots` is an error.
- `PATHS` lacks `mappings.bed`, `genes.gff3`, `operon_pairs.tsv`, `regions.bed`; add `codon_start` (pileup, variants) to `ONE_BASED_STARTS` when those keys exist.
- Example 10 runs only on FarmShare; the colour cycle leaves out black and grey (spec wording).
- `read_tsv` now treats an empty cell as missing in genes.tsv too (matches `giap.io.genes_tsv`).
- Mutation evidence: $GIAP_ROOT/review/pr191-mut/ (22/25 killed, 1 equivalent).

## M9-03 review nits (#210; one round, oriC rounding and the 2^20 gzip boundary pinned in the PR)
- Untested guards in `giap.viz.support`: `window_bp < 1`, the missing-contig `ValueError`, the empty-sequence skip.
- `datapackage.json` lists `plots/support/gc_skew.tsv` uncompressed; a `.tsv.gz` write (not reachable at 10 kb windows for bacteria) would leave the descriptor pointing at a missing file — write the descriptor path from the file actually written.
- `plots` step: with `replicons` absent but a `plots_support` input passed, plots render from the incoming directory, not this run's output.
- No test checks that `write_datapackage` copies `schemas/plot_support/*.json` into OUT.
- ADR-0042: when the M7-10/M7-13 (`calsuite_variants`) and M9-22 (`calsuite_compare`) calibration tests land, extend scripts/slurm/calibration_full.sbatch to `--array=0-14` (task 13 = variants, 14 = compare); until then those two keys have no tests and the array skips them.
- M9-04 V17 review nits (#197): when every gene is partial, (e) draws 0 % bars and should show an empty-state note; integer tick locators on the count axes; partial genes are left out of (d) although the spec restricts only (e); the "unmapped" layer of (c) includes contaminant and too_short peptides; the length-bound lines belong at 6.5/40.5 (spec + test); an empty `replicons` draws an empty axis with a warning; the TSV message rows drop the drawn text; preflight `genome.contiguity` has an empty `message` (a preflight bug); "abundance not available" is not in the spec's titles; "1 preflight warnings" should be singular; the imputed legend shows when nothing is imputed; 8 surviving mutants (exactly 10 warnings, rule fallback, (b) floor, circular markers, constant-abundance padding, (e) inverted y, absent-genes title).
- ADR-0041 follow-up (after #196 merges): move the four tokens from RUN_FLAGS_PENDING to the accepted run vocabulary and update the vocabulary test (lead).

## M9-02 review nits (#205; one round, B1–B5 and N1 local half fixed in the PR)
- N1 Nextflow half: REPORT gets no `events`, so `layer_missing` skips fire under the local executor and `giap report` but not in a Nextflow `giap run`. Fix before the first Nextflow-run plot that declares layers.
- Config digest reads "not available" in `giap run` but appears after `giap report`, so `report.html` changes on a rebuild of an untouched run.
- Executor-equality test compares only `plots/data` (empty until plots register): also compare `report.html` bytes, and the run's report with a rebuilt one.
- Peptide "of N input" counts `kind=protein` rows from `unmapped`.
- Stale plot files from an earlier run are shown in place of a skip reason that holds now; with `formats=[png]` the report says "plot file absent" although the PNG exists.
- PLOTS has no resource label (1 CPU / 2 GB) while REPORT_HTML gets `process_low`.
- `ALL_PUBLISHED_KEYS` declares optional inputs never read (`peptide_hits`, `tus_tsv`); SVG id separator `-` is not provably unique.
- Surviving mutants: `<=` boundary in the card counts, the two `nunique` peptide counts. Move `cards()` to `giap/report/cards.py` ahead of M9-16.
- Card/table read failures are notes in the HTML only (build.py is pure); consider an event from the step.
## M8-17a review nits (#204; one round, B1 R_LIBS_SITE fixed in the PR)
- 17b: `read_results` must require rows to match the in-range sets exactly (a dropped row now vanishes silently); `read_versions` should check the `seed`/`jobs` lines and the data.table version.
- Set `TMPDIR` to the call directory in `child_env` so R's `Rtmp*` is not left behind after a timeout kill; give the post-kill `proc.communicate()` a timeout.
- `write_call` accepts IDs with tab, newline or `;` (corrupts the R inputs): reject them.
- `scripts/build_envs.sh`: a failed lock check on an existing env leaves the `.giap-built` marker; the marker hashes only `fgsea.yml`, so a lock-only change never rebuilds.
- `giap_fgsea.R`: an empty input file exits 1, not 3.
- `FgseaEngineError` imitates a `StepError` and is then wrapped in a second one by `run_step`.
- Surviving mutants: duplicate leading-edge check, size-range check, trailing-newline check, `R_HOME` filter, `set(edge) <= ranked` check.
## M9-11 review nits (#208; one round, B1 overlap fixed, B2–B5 pinned by test_v11_review.py)
- Title medians come from pairs.tsv, so the data TSV cannot reproduce them (index decision 5): widen decision 4's exception to cover medians, or add per-pair (a) rows; (a) rows carry a bin index but no edges.
- Flag legend handles take the first point's class colour; use a neutral face (`other` can be OP or NO).
- Surviving minor mutants: outline colour, flagged points drawn twice, flagged markers not class-coloured, `set_xlim` removed, degenerate (c) range padding, legend `ncol`.
- (a)/(b) x-label "|Δ| partner log2 abundance" reads awkwardly; the spec's title template gives "1 pairs".
## M9-07 review nits (#203; one round, B1 out-of-range p and NITs 2–3 fixed in the PR)
- V18 fit title and corner text ignore `recalibration_separated = 1` (prints "slope 41"): append "; recalibration separated".
- "1 test families shown": singular for k = 1.
- Carry the slope/intercept Wald CIs (`ci_low`/`ci_high`) into the TSV and corner text (pass 2).
- Both real fixtures are constant fits; add a non-constant-fit golden when M4-13's fixture allows (the synthetic render's bin n sum 875 ≠ n_eligible 100).
- Open in the spec: where `n_outside_0_1` sorts (now after the 20 bin rows); a family with every p outside [0, 1] now draws 20 empty bins plus the outside count.
## M9-06 review nits (#202; one round, B1 fixed in the PR)
- V7 title: when no row has an estimated CI, say "O − E CI unestimated in every group" (not the vacuous "matches … with a CI"); basic-only title uses `%.1f` (→ "-0.0"), switch to `%.2g` and prefix the block ("COG J:").
- Dead E = 0 branches (`detection.py` o/e helpers); `_pct` duplicates `base.title_pct` (widen it to floats); share the `flags` parsing.
- Region rows with NA O − E show only "CI unestimated": add "O − E NA" to the label.
- `ci_few_units` rows are drawn as ordinary CIs and can headline the title; `small_region` rows can headline too — lead reading needed in hardening.
- Figure height grows 0.17 in per row uncapped (26 COG + ~100 regions > 20 in): cap or hand to M9-26.
- V8 "none" grey equals the "poorly" super-category grey (two identical legend patches): hatch "none" or merge the entries.
- V8 strips lack the median tick the violins have; title wording at `detection.py:328/331`.
- M9-05 V6: a run with targets under both `reverse` and `shuffle` decoy types (database + de novo placements) — split or overlay the strata (lead decision 3, 2026-10-22).
- M9-05 review nits (#198):
  - **Labels:** N4, the "rejected by FDR" label is pushed beside the placed bar with no leader line, and the (b) "0" label sits on the α line.
  - **Test gaps:**
    - N5b, no case with targets under both decoy types;
    - N6, no test with a contaminant label and `reason=no_hit`;
    - N7, no title test with FDR = α, and none with every stratum estimated and above α;
    - N8, the reviewer's alluvial geometry check (scratchpad geom.py) should become a test.
  - **Wording:** N9, the sky-blue "rejected" ribbons are hard to tell from the placed ones; N10, "strata 0" should be singular, and the (c) label should say "95 % CI".
  - **Size:** about 25 lines of simplifications listed in the review.

## M9-12 review nits (#219; one round, decisions 11–15)
- V13 start fallback (NA `start_longest`): the span starts at the first placed peptide anywhere on the
  contig; a multi-placement peptide can stretch it, and origin crossing is ignored there.
- One failing `render_locus` call aborts all of V13 (panel (a) included); V3 behaves the same. Render
  each view in isolation and record a per-view skip.
## M9-14 review nits (#221; one round)
- No test pins escaping of data values in twin JSON independent of plotly's own `<` escaping (a
  `json.dumps`-built test twin would kill the `_script_json`-removed mutant); the probe was inert.
- A string-column `NA` in a data TSV reaches a twin as NaN and is written as `"nan"`: M9-15 twins
  should read data TSVs with `keep_default_na=False` or map NA explicitly.
- JS vs Python string order (code point vs UTF-16 unit) and lowercasing (final sigma) differ in rare
  cases; document in the methods/report notes. The 5 `tables.js` mutants survive (decision 1).
- `giap run --strict` does not map `report.size_exceeded` to exit 3 (only preflight warnings reach
  it today); wire step-level strict warnings once a second one exists.
- `remove_uids=True` in `figure_json` is unpinned; an instance with a twin but no SVG still gets
  "plot file absent".
### M9-08 (V01 circos) after the review round
- No unit test renders V01 with empty promoter/terminator/TU tables (crash fixed in 8faffdb was
  caught only by integration tests).
- Inner oriC/ter labels of small adjacent replicons crowd each other near the centre (five-replicon
  review render); consider dropping them below a sector-size threshold.
- Counted lines 629 (> 400 warn): split the label-placement code into its own module in hardening.
## M8-17b review nits (#214; one round each, code review and statistical review)
- M8-15 (statistical re-review, before the gate): re-run the `test_attack_gsea_flags.py` arms on the R
  engine, flag spillover included (unflagged rows in a collection that also has flagged rows; only
  measured under the port, Slurm 1867330), and re-measure the N4 detection-screen and
  `genomic_cluster` φ = 0.7 numbers on R. The methods note makes no claim for these until then.
- Weekly `gsea` array task: 4 CPU / 2 h now also runs the `requires_r` replicates (Slurm 1869053 alone
  took 47 min on 8 CPUs) and needs `giap-fgsea` built; `tests/calibration/_criteria.py` does not
  exclude `requires_r`, so the harness's `gsea` suite errors where the env is absent.
- Untested in `giap.function.gsea`: dropping `lo <= n <= hi` for whole-ranking sets; the
  `len(tus) < min_size` and `bh(sub)[i] > alpha` boundaries; `SIMPLE_ENGINE` re-seeding per job.
- `giap.io.fgsea_r`: an explicit `GIAP_FGSEA_PREFIX` skips the `.giap-built` check (an unbuilt env
  says "is not a directory"); untested with no `R` in `bin/`; `engine.provenance` keeping its value
  after a failed call is unpinned.
- The `function.gsea_engine` event could log the seed R ran with.
## M9-09 review nits (#218; one round, decisions 21–26)
- V2 trend memory: up to 2.2 GB peak RSS at 20,000 points (the PLOTS process gets 2 GB on the first
  attempt); fit the trend on ≤ 5,000 bin medians.
- V2 rug-only panels show meaningless y ticks (±0.04); hide them.
- V2: oriC at position 1 is hidden by the left spine; the landmark −1 is untested (decision 3).
- V10 data TSV keeps the producer's values for `not_estimable` rows (flags alongside); the plot
  blanks them (decision 21). Revisit if a reader takes the zeros as estimates.
- Statistician note (TU ordinals in the islands detrend vs bp positions in V2) is in the PR body;
  move it into docs/methods with the M9 methods index.
## M5-11 (statistician, hardening)
- Island power under giapsim's default `covariates_only` detection: no island at +3/+5/+10 SD on 1,000 genes (ℓ̂ = 1, ℓ' = 3; ~45 % censored regardless of abundance, rank scores saturate). The M5-11 golden uses sim4500 with `detection="abundance"`; decide whether the power table should report the covariates-only regime.

## M5-11 review nits (#217; one round)
- The `giap-bench` sbatch script used for the `context_m5` benchmark (Slurm 1869098) is not in
  scripts/slurm/; commit it beside the other benchmark launchers.
- M5 gate evidence: full `islands` records exist (Slurm 1868687, R = 500, all pass); the full
  `autocorrelation` record does not — submit `calibration_full.sbatch --array=4` on main. Calibration
  records carry no commit field, so "on main" cannot be checked from them.
- The island golden's join-count rows depend on the lbfgs detectability fit (decision 44); a
  BLAS-stable `p_detect` would let them return to byte equality.

### M9-18 (V12) untested decisions
- Tests for M9-18 decisions 6–9: NA window for a gene missing from genes.parquet, in-frame cover
  rule, "(N peptides)" line, header-only `variant_windows.tsv`.

### M9-18 review nits (PR #227)
- N1 peptide ending exactly on the codon's last base; N2 ref/pep column alignment (incl. a window
  clipped at the protein start); N3 `_best` untested; N4 duplicate `locus_tag` rows untested.
- N5 `stop_codon` flag shares vermillion with the bold ALT residue.
- N6 title edge cases: only mixed sites ("0 … 0 … unestimated"); min == max prints "0.004–0.004".
- N7 panel (a) y-label "ALT sites" also covers mixed bars.
- N8 "(N peptides)" can differ from `n_alt` (q > 0.01 placements, IL-only, other ALT residues); a
  block before the site misplaces it.
- N9 misleading comment about sites past the protein end (support.py).
### M9-16b follow-ups
- Reword "run_manifest.json not written by this GIAP version" once the tail re-render lands (only
  shown when the manifest write failed or in the work-dir copy); test-author step first.
- Schedule MVP-03c (tools[] / databases[] in the manifest) before the M9-25 approval packet.
