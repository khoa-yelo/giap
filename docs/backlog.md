# Backlog

Issues may live on GitHub; this file mirrors the milestone plan (GIPPA_SPEC.md §13).

- [ ] M0 scaffold (in progress)
- [ ] M1 inputs and preflight
- [ ] M2 annotation wrappers, full Bakta DB, recorded fixtures
- [ ] M3 mapping and decoy FDR; pilot kit
- [ ] M4 gene table, detectability, region tests
- [ ] M5 autocorrelation and islands
- [ ] M6 operons; RegulonDB calibration (source currently unreachable)
- [ ] M7 variants and novel ORFs (declare GIPPA's VCF FILTER values with descriptions; M1-05 writes generic ones)
- [ ] M8 function and determinants
- [ ] M9 visualisation and report
- [ ] M10 pilot validation
- [ ] M11 release 1.0

- M3: `mapping.tolerant.swaps` and `mapping.tolerant.mass_equivalents` config keys are not consumed by M3-07; decide which issue implements them or remove them (from M3-07 review, PR #48).
- M7 (moved from M6 by M6-03 lead decision 2; M7-12 or M7-16): the simulator plants stop-codon readthrough (C-terminal extensions, ADR-0013) for cterm_extension calibration (from M3-19 test-author note).
- M2-09: antiSMASH databases with `url: "image"` (databases bundled in the image) — setup edge case from M2-08 review.
- M2-09 (from M2-08 review, low): delete the geNOMAD tarball after unpacking (§3.9); reject non-200 responses before hashing; cap redirects; test `_provenance` requires prev status ok; drop redundant tar `filter` weakening check.
- M2-09 follow-up: GO release directories can announce a different `data-version` date in the OBO header (2026-06-19 → 2026-06-15, 2026-08-05 → 2026-07-26), so `verify_go` rejects them and the pin stays at 2026-05-19; verify the header against a pinned `data_version` instead of the directory date.
- M4: migrate existing calibration tests to `_criteria.calibrate()` so ADR-0015 two-block confirmation and sharding apply to all suites (M4-05 decision 7; test-author).
- M4 (calibration migration follow-up): suites marked only `calibration` (stats_*, placement lists/genome) have no `calibration_full` variant, so the nightly full array skips them; adding one needs a new R — statistician decision.
- M4 (statistician, #69): gate the sim_proteomics q-value FDR check and the one-sided de novo coverage checks through `calibrate`.
- M4 (statistician, #69 minor): the 1.5q teeth test at q = 0.01 has power < 1 at fast R (seed-fragile) — raise R or use 2q at q = 0.01; teeth tests rerun in every shard task — run them once.
- Declare threadpoolctl (used by gippa.genetable.detectability) as a direct dependency at the next uv lock update (M4-13).
- M3-17 follow-up: wire the blastp contaminant check for unassigned proteins into the mapping steps (only peptide library matching is wired).
- M3-17 follow-up: example 5 (nextflow executor end to end) is xfail(strict=False) until M2-20 (context_layers) merges; remove the xfail then (decision 16).
- M4-20: simulator `multi_identical_pairs` should give paired genes identical `aa` in the layout (M4-09 decision 2).
- M4-14: record the declared abundance type (e.g. maxlfq) in the input metadata; C10 has no scale column (M4-09 decision 3).
- M4-09: add a benchmarks/ entry for simulate_evidence on 4,500 genes (budget 3 s; unit test keeps a 30 s guard).
- M4-20 (with M4-09 decision 2): for shared_frac, the simulator's second hit sequence is not in that gene's aa (M4-09 review).
- M4-14: re-run the abundance Spearman sanity claim on M4-09's evidence simulator (logistic detectability, missed cleavages) once #80 is merged.
- M3-17/M4: resolve makedb protein IDs (GIPPAORF_, NOVEL_, *_patch) through db_id_map.tsv in annotate_protein_hits instead of flagging cross_strain_protein (M3-23 review).
- M2: the full-DB Bakta recordings (bakta_full/*) were made without a replicon table, so every replicon is 'linear'; re-record with preflight's replicons.tsv (--replicons) so circular wraps are exercised (M3-23 review).
- M2: one shared IUPAC-aware reverse complement (ostir.py _COMP, bakta.py, operons.py only complement ACGTN; preflight allows IUPAC) (M2-16 review).
- M2-15: test the manifest `uncalibrated` flag and step-level `cross_check_transtermhp: false` once manifest and optional-output support exist. **The `uncalibrated` flag part is closed by MVP-03c (#236); the cross_check part is open.**
- M2-15: parse_bactermfinder's input format is provisional (no recording); record a real BacTermFinder run (opt-in tool, env not built) and pin the parser to it.
- M4-20: replace M4-16's test stand-ins (p_detect_true, fixed COG table, tier/detectability columns) with the real M4-12/M4-13 outputs once merged.
- M4-14 (statistician): re-measure abundance.ibaq_denominator on M4-09 data before switching (fails 10–40 % at p_obs = 0.3); restate as a pooled mean Δρ if needed; add the p_obs = 0.3 attack.
- M4 notes: median imputation pulls detected-but-unquantified genes toward the centre; the min/exclude sensitivity reruns bracket it (statistician, #82).
- M7-01: co-optimal alignments with equal cost, diff-type multiset and difference positions but different op shapes are broken by DP move order (1:1, 2:2, 1:2, 2:1), not pinned by a test; the frameshift property oracle sorts kinds alphabetically, not by DIFF_TYPES precedence (no divergence found in 200k random junction cases).
- M7-01: align_hits loops rows in Python (~110 µs/hit); batch the DP in numba if Module 5 ever aligns far more than 10^5 hits.
- M7-02: frameshift site `alt_aa` is the placeholder `fs+1`/`fs-1` and `ref_codon` the last prefix codon; settle both in M7-03 with the VCF writer (lead decision 3). The M7-03 writer must drop the internal `_ref_codon`/`_alts` pileup columns; pandas `attrs` (frameshifts, contig order) do not survive `pd.concat`.
- M7-02: a folded peptide with several categories at one codon (only possible with conflicting alignments of one sequence) is counted as its smallest ALT key, not pinned by a test (the oracle rejects such cases).
- M4-16 (statistician #87 sign-off, non-blocking): (a) also flag ci_few_units by the number of TUs carrying a detection (non-coverage 0.07–0.12 for 1–5-gene TUs, 0.36 for 10-gene TUs at E≈10; job 1850360); (b) genome_inconsistent_rate only under --genome-relation same (§7.5.1), not related; (c) `hypothetical` parsed with astype(bool): the string "false" reads as True — parse booleans explicitly.
- Benchmarks: nothing reads benchmarks/baseline.json yet — add a check that compares benchmark results with it (M4-13 review).
- At the next uv lock update: scikit-learn>=1.8 (l1_ratio API) and threadpoolctl as direct dependencies (M4-13).
- M4-15/M4-17/M8 attack-2 suites must gate on M4-13 p_detect; even the correct model slightly over-predicts detection for the shortest genes (decile z ≈ −0.66) (statistician, #83).
- m3/mapping-nf: write the search_db auto-resolution event from the Nextflow `gippa _step` path too (only --executor local writes it, M3-17 decision 14).
- M4-12 (statistician #92): gippasim mode with decoy-derived q (D/T) to separate MAYU with/without PEP; calibration run with q_value_max = 0.1; document 'best-scoring decoys lost' as a limitation of the 80 % pairing guard.
- M4-14: use Detection.protein_gene_map from M4-12 instead of re-deriving the map (statistician #92).
- M2-10: BacTermFinder recording (env not built) and record_fixture.py examples 1–3 deferred.
- M2-10 (#93 review): re-record the antiSMASH fixtures on the circular K-12 GBFF (current ones used the linear GBFF; ~18 s each) or note it in their PROVENANCE.
- M2-10 (#93 review, R6): scripts/fetch_test_genomes.py must verify checksums by default (NCBI Datasets md5sum.txt or mandatory --compare); tag_circular should check records instead of tagging all.
- TE6: pin Bakta's K-12 oriC call count at two (true origin + dnaA DnaA-box cluster); the replichore analysis (§1447) needs a rule for choosing the true origin.
- Spec note: Nissle GCA_043228065.1 has three plasmids (p15A-ILA-dapA, likely engineered; pMUT1; pMUT2) — the spec's edge case names two.
- Bakta fixture PROVENANCE: add attribution (Bakta DB CC BY 4.0, doi 10.5281/zenodo.14916843). Later recording waves need a size budget (wave A uses ~48 % of 50 MB).
- nf_process_real (job 1849282 K-12, 1849365 Nissle): PROMOTECH and UNIOP report 'failed' on K-12/Nissle only because no recorded inputs exist there. Needs (a) record_tools_inputs.py builders: Promotech targeted windows via the merged M2-13 promoter_windows step, UniOP CDS GFF from the Bakta JSON (format of tests/data/recorded/uniop/*/sim100k/uniop_cds.gff); (b) recordings of both tools on K-12/Nissle under $GIPPA_ROOT/recordings (record_tools.sbatch: drop the sim100k-only promotech skip, add uniop); (c) the array's PROMOTECH/UNIOP branches reading them.
- Calibration log (tests/calibration/…/regions.jsonl, append-only): each record needs a run ID + commit so stale records (e.g. pre-pooling crispr hl_shift 'fail') are distinguishable (statistician #89).
- M4-17 nightly: the A-LIN full attack (test_alin_linear_replicon_hl_at_level) takes ~57 min per block on one core — shard it in the nightly/weekly calibration job (statistician #89 sign-off).
- M2-20 follow-ups (implementer report, real run 1851270): annotate.nf passes no `tool_dirs`, so OUT/annotation/tools/ is empty (§7.1.8 wants native tool outputs untouched); `cai`/`fop` all NA on sim100k layers — check whether expected (no reference set on the simulator) or a wiring gap; consider splitting schemas + gc_local_dev if review asks (658 gated lines).
- M2-20 (#104 review): annotate.nf passes no tool_dirs → OUT/annotation/tools/ empty (§7.1.8 native outputs untouched); nothing produces id_map, so --annotation GFF3/GBFF runs lack annotation/id_map.tsv. Wire both.
- M4-18 (lead decision 4): the dosage-adjusted plasmid row treats replicon depth as known (CI shifted by log2(plasmid / reference depth)). Deferred: propagate depth uncertainty (e.g. bootstrap or delta-method term from per-window coverage); the manifest should record a missing `--replicon-depth` (the flag itself is wired into `context_stats` since MVP-01).
- M4-17/M4-18 (lead decision 17; statistician #110 finding 7): a shared `ci_few_units` flag for region-test and plasmid TU-block CIs computed from fewer than 10 TUs (docs/methods/bootstrap.md: percentile CIs from < 10 TUs under-cover), as M4-16 already does for detection_summary.
- M4-15 decision 10: "file I/O only through io/" is enforced for genetable/assemble.py only. Step-level convention open: merged steps (mapping, context, operons, codon, ...) read Parquet directly with pd.read_parquet; decide whether steps must go through gippa.io readers and add an AST check like test_assemble_module_is_pure if so.
- §9.5 per-field `unit` in the table schema resources (src/gippa/schemas/tables/*.schema.json): none of them has it yet (M4-15 decision 16; project-wide).
- Promotech sim100k recording ($GIPPA_ROOT/recordings/promotech/sim100k) is a genome-wide -g run with every row labelled chr1,plasmid_1; re-record in GIPPA's -s window mode (M2-21 test-author).
- Nightly: test the newest open PR branches too (M2-21 decision 5, deferred). containers.yml / GHCR images deferred (M2-21 decision 2).
- Benchmarks: gate only steps with a baseline ≥ 1 s (M2-21 decision 10); implement the threshold in the benchmark check.
- M4-19 decision 7: the run manifest records `output.sensitivity = off`; decision 8: no rows for a module that is off. (Decision 16, registering and wiring `run_sensitivity`, is done: MVP-01.) **Closed by MVP-03c (#236): `sensitivity.level` and `n_rows: 0`.**
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
  $GIPPA_DATA_DIR/proteomes/UP000000625.fasta (used by the M3-24 pilot cfgs), while
  validation/fetch_uniprot_k12.py and scripts/slurm/weekly.sbatch use $GIPPA_DATA_DIR/uniprot/. Pick one.
- MVP-01 pass 2 (found while wiring): with a precomputed `regions` layer ANNOTATE has no C7/C8
  Parquet, so CONTEXT (context_stats, sensitivity) skips itself silently — the run should say so,
  and ideally the precomputed layer should reach M4 as the contract tables.
- MVP-01 pass 2: nothing produces `id_map` under Nextflow (C14 is an optional input of
  `gene_table` and `context_stats`; the local executor stages `annotation/id_map.tsv`), and the
  `tables/` outputs are published but not validated inside the pipeline (contract assertions).
- MVP-01 pass 2: the registered `sensitivity` step writes a header-only `tables/sensitivity.tsv` **Closed by MVP-03c (#236).**
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
  CPU because `gippa.context.regions.residualise` (ADR-0018 decision 1) calls `np.linalg.lstsq`
  (LAPACK SVD); the three CI jobs (ubuntu 3.11, ubuntu 3.13, macOS 3.12) disagree with each other, so
  no recording machine satisfies the matrix. The other three golden tables and every other column of
  this one are byte-identical. Fix options: compare the golden's float columns numerically
  (test-author), round what the `regions_stats` writer emits (output-format decision, helps users
  too), or make the OLS BLAS-free (moves every calibration number). Evidence and bit patterns in
  PR #133.
- M4-20 (decision 2): `run_detectability` now writes `tables/genes.tsv` (`write_genes_tsv` had no
  caller at all), so the M4 chain and the local executor produce the §9.3 master table. It is *not*
  a declared step output yet, so `gippa run` does not publish it: adding the key `genes_tsv` to
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
  `tables/genes.tsv`, which `gippa run` still does not publish — see the M4-20 item above). Settle
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
- M5-09: `gippa.steps.registry` gained `optional_outputs` (a subset of `outputs` a caller may leave
  unwired) so the M5 keys could be added to `context_stats` without rewiring M4-18's registry-path
  callers. It has no dedicated unit test yet — the test-author should cover the decorator check,
  `_check_keys` and the "wired together or not at all" rule in `_m5_context`.
- `np.linalg.lstsq` in `gippa.context.regions.residualise` is not bit-portable across CPUs (LAPACK SVD),
  so `regions_stats.tsv` float columns differ by 1–2 ulp between machines and between CI legs (M4-20
  decision 9, probe 1863133). Writer-level rounding is the real fix — see the §9.5 per-field precision
  item — and belongs with the M6 writers.
- tests/helpers/islands_data.py (m5/islands) and tests/helpers/oriter_data.py (m5/oriter) both carry a
  `tu_frame` and a gippasim-scenario driver; whichever merges second should deduplicate. Note the ori–ter
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
  conclusion (GIPPA_STEP runs natively) is from reading the config, not from a run.
- RegulonDB TLS: the server omits its intermediate, so any plain client fails verification. GIPPA completes
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
- M8-07 dedupe: `gippa.determinants.features._eligible` is a local copy of the §7.3 regions-table rule
  (`regions_table._eligible`); switch it to `gippa.genetable.eligible` once M8-03 (#159) merges
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

## M8-10 review nits (#178, one-round MVP review: APPROVE)
- **N1** Log the `--publication` runtime at run time; today only the CLI help states it (lead decision 3).
- **N2** `_threads` is capped by sched_getaffinity, not by `task.cpus`. Under the local executor or a container it can oversubscribe; determinism is not affected.
- **N3** farmshare.config: a `process_high` retry asks for 64 GB over 16 CPUs, which is 4096 MB/CPU and over FarmShare's 4000 MB limit. Fix it in farmshare.config; base.config's `cpus = 16` is pinned by a test.
- **N4** The step uses private helpers from steps/codon.py and steps/operons.py, plus its own C5 reader. Move them into a public contract reader in io/.
- **N5** The `_cell` `""`→NA branch is redundant or untested (mutant m07 survived).
- **N6** Nothing checks that CAI reference rows exist (m04 and m12 survived).
- **N7** Only the bootstrap stream is pinned. Pin the CV stream name and the null seed too (m10, m11).
- **N8** Add a step-level test for NaN `n_perm` and `unestimated` on null_p (m02). The `features_dropped` test checks a subset rather than equality (m06).
- **N9** main_determinants.nf.test should assert where `tus` and `replicons` come from.
- **N10** One main.nf header line is over 100 characters.
- **N11** Still owed: a real gippa_run.sbatch run on sim100k, and the 4,500-gene runtime in benchmarks/baseline.json.
- Mutation evidence: Slurm 1867779; scratch is in $GIPPA_ROOT/scratch/review178_mut/.
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
- M8-02: `fgsea_simple` raises an even `sample_size` by one (fgseaMultilevel does); GIPPA's config is
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
- M6-06 (pass 2): `bin_counts` and the π₀ actually used must reach the §9.4 run report and the manifest when the operons step is wired (M6-07/M6-08); the function returns them but nothing writes them yet. **Closed by MVP-03c (#236).**
- M6-06 (pass 2): the missed-operon benchmark baseline was measured in the iron-15 session allocation only (0.23 s); re-measure on wheat with the other M6 entries.
- M6-06 (pass 2): the lfdr of targets in a bin with no decoys is too small (methods note, "Known asymmetry"); set it NA or 1 there before lfdr is ever thresholded.
- M6-05 review of #147 (code reviewer, statistician), owed to M6-05H/M6-09:
  - skip the t fit when the scale is relative (today it fits and may warn "too_few_op_pairs" on relative rows);
  - log the 1e-6 scale floor in the non-convergence warning, not the raw MAD;
  - document that one-detected rows get no `no_detection_limit_model` on the too-few/relative paths;
  - test gaps (surviving mutants): (undetected gene with a finite imputed value: done in M6-05H, #311); `unestimated` on the too-few and relative paths; one side partial only; pair flags (`rna_between`) kept; flags of partial rows;
  - statistician F3/D-1: abundance-dependent Δ noise breaks the lowest-quartile level (2.4–3× nominal) — method ADR or narrowed claim, Khoa decides; D-2 (TU-chained Δ) and D-3 (power record) as attacks;
  - F5: keeping P(miss) out of the BH family in M6-05H needs an ADR; F6: flag imputed-but-detected tested pairs `abundance_imputed`.
- M6-05 pass 1 simplifications (implementer): `discordance()` accepts `cfg` and `p_miss` but uses
  neither — a NaN abundance is untested under every `abundance.missing` mode, not only `exclude`, and
  the one-detected branch is always `unestimated (no_detection_limit_model)`; `PMissModel` is a
  placeholder Protocol until M6-02 merges (replace with its real type in M6-05H). **Closed by
  M6-05H (#311):** the Protocol stays as the adaptor, renamed `PMissFn` (`p_miss(u, a_v, ratios, v)`);
  `detection_limit.p_miss_for_genes` binds M6-02's `PMissModel` to it.
- M6-05H (#311): F5 ADR owed: keeping P(miss) out of the BH family (M6-index lead decision 1) departs
  from §7.4.3's "both kinds together"; write it with M6-09 attack C, which measures the mixture.
- M6-05H (#311), deferred: the `operons` step does not read `peptide_runs` (C16) yet, so the
  replicate model is never fitted there, and it does not join the M4-10 covariates. A replicate model
  with covariate coefficients is guarded (`pmiss_covariates_missing` + WARNING, review S2); wire
  `peptide_runs` and pass `cov=` from the covariate table together.

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
  - the real sim100k `gippa_run.sbatch` run;
  - the step benchmark.

## M7-05 review nits (#175, one-round MVP review)
- Test a minus-strand CDS crossing the origin for `start_codon` (the `first %= rep.length` mutant survived).
- M7-16: decide `stop_codon` and `mass_equivalent` for blocks and swaps whose REF contains `*`.
- The §9.4 flag order is defined twice (`flags.FLAG_ORDER` and `vcf_records._FLAG_ORDER`); derive one from the other.
- `site_flags` rebuilds `start_codons` per call; point bulk callers to `flags_for` in its docstring.
- The FLAGS INFO description in vcf_records.py and variants.nf has no sync test (only the IDs are compared).
- MVP real run on 813 genes (Slurm 1865149, validation): operon tests now run (515 t, 102 decoy), but:
  - detection is saturated: 804/813 genes are detected, `p_detect` is the constant 0.989 fallback, and picked FDR q = 1/804 for all. The gippasim real-run scenario needs decoys with null scores (the earlier item) and lower peptide coverage;
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
  - `flags.py` imports the private `_bool`, `_flags`, `_per_gene` helpers from `missed.py` (also duplicated in `discordance.py`); move them to one shared `gippa.operons._util` module;
  - `operons_table` restricts `n_flags` to flagged **OP** pairs with both genes in the TU, so `sum(n_flags) <= flagged OP pairs` holds by construction even if a C5 TU ever contained an NO pair; confirm with the lead or drop the restriction;
  - `detected_frac` for an all-partial TU is NA (0/0), as the spec lists it owed to pass 2;
  - `merge_tests` keeps the first row when a `pair_id` appears in both M6-05 and M6-06 rows (cannot happen: OP vs NO), with no warning;
  - `evidence_note` wording is placeholder prose (spec "Hardening owed"); `layer_missing:*` is not added to SW rows (pass 2);
  - `flag_pairs` builds notes with `iterrows` (fine at thousands of flagged pairs; vectorise if a benchmark says so).
- M6-06 (pass 2): `bin_counts` and π₀ now reach `logs/events.jsonl` as the §6.7 event `operons.missed_operons` (M6-08); the §9.4 run report and `run_manifest.json` still have no writer that reads them. The first bin's `lower_bp` is written as `-Infinity` (Python `json` default), which is not strict JSON. **Closed by MVP-03c (#236).**
- M6-06 (pass 2): `bin_counts` and the π₀ actually used must reach the §9.4 run report and the manifest when the operons step is wired (M6-07/M6-08); the function returns them but nothing writes them yet. **Closed by MVP-03c (#236).**
- M6-06 (pass 2): the missed-operon benchmark baseline was measured in the iron-15 session allocation only (0.23 s); re-measure on wheat with the other M6 entries.
- M6-06 (pass 2): the lfdr of targets in a bin with no decoys is too small (methods note, "Known asymmetry"); set it NA or 1 there before lfdr is ever thresholded.
- M6-07 review of #152 (non-blocking):
  - tests to add: the terminator cut-off on its own (promotech 0.5 / terminatornet 0.9 / score 0.8); a pair's own flags kept by `merge_tests`; the both-detected guard in `max_abs_delta_log2`; a gene missing from the gene table counted as undetected; the fully linked circular ring (wrap OP pair inside one TU);
  - validate `missing_layers=` ∈ {promoters, terminators}; guard `merge_tests` against empty inputs that have no columns;
  - move `_bool`/`_flags`/`_per_gene` to `gippa.operons._util`;
  - M6-10: restricting `n_flags` to OP pairs matters once pairs are reclassified at UniOP 0.3/0.7 while the TUs stay at the default.
- M6-08 pass 1 (implementer):
  - merge time: add `pairs`, `operons`, `operon_flags` → `tables/*.tsv` to `gippa.runner.local.PUBLISHED_FILES` (#145; `test_local_publish_names.py` and the slow M6-08 local tests need them);
  - `run_operons` returns a 4-tuple (pairs, operons, operon_flags, bin_counts), not the spec's 3-tuple, so the event needs no second missed_operons call;
  - `schemas/datapackage.json` does not exist in the repo; the three resources are per-table files like the others (spec acceptance item);
  - `benchmarks/baseline.json` runtime of the step on a 4,500-gene genome is not recorded yet (spec acceptance item);
  - the publishDir rule (only `tables/*.tsv` in a subworkflow's out spec) is written in `operons.nf`; `genetable/genes_pre.parquet` from MVP-01 still publishes;
  - `peptides`, `peptide_hits`, `peptide_runs`, `detection_summary` are declared but not staged by OPERONS nor read (M6-02, pass 2);
  - the `gippa.steps.common` importers alias `one_path as _path` etc.; rename the call sites when the modules are next touched.
- M6-08 review of #155 (approve; pass 2):
  - **must-do before the manifest writer:** open bin bounds `±Infinity` in the `operons.missed_operons` event are not strict JSON; write them as `null` in `_plain`; **Closed by MVP-03c (#236).**
  - local executor stages `peptides` and `detection_summary` into `operons` (and orders it after `context_stats`); OPERONS under Nextflow stages neither; reconcile before M6-02 reads them;
  - missing-layer WARNINGs have no §6.7 `operons.layer_missing` event;
  - `assign_tus` on precomputed `tx_rates` in `_layer` is dead work (`pair_frame` never reads `tu_id`);
  - `_cell` stringifies `None` items inside lists, bypassing `join_list`'s missing-item guard;
  - nf-test lives at `workflow/tests/operons.nf.test`, not `workflow/tests/subworkflows/`;
  - test gaps (surviving mutants): `missed_operons` without `genes=` (needs a relative-scale/partial-gene fixture); precomputed branch of `_layer` (promoters → `promoter_windows`) has no unit or nf-test; `no_detection_limit_model` log level unpinned.
- MVP real run (Slurm 1864792) is degenerate by input, not by code (validation diagnosis):
  - 82 candidate genes, so the specified picked FDR floor (0 + 1)/82 = 0.0122 > α = 0.01 (§7.3.2, GIPPA_SPEC.md:1331, worked example at :3253), and no gene can be detected;
  - then `p_detect` = 0, the detectability constant fallback; `detected_2pep` and `abundance_log2` are independent of the call;
  - gippasim (tooling): the real-run scenario needs several hundred genes, and decoy rows are written with NA q/pep/score so a decoy can never win; draw them from a null;
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
- Tests: an island wrapping the origin of a circular replicon; pin the replicon order of draws (two replicons near threshold); a direct `islands_family` test with `stats.permutations != 199` and an automatic block length (mutants: topology flip, length + 1, reversed replicon order survived; $GIPPA_ROOT/review/pr184-mut/).
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
- Mutation evidence: Slurm 1867940, $GIPPA_ROOT/review/pr181-mut/.

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
- Mutation evidence: Slurm 1867907, $GIPPA_ROOT/review/pr180-mut/.
## M8-04b reviews (#179: code review CHANGES, statistician BLOCK, fixed in decision 13), pass 2 and nits
- **Pass 2: lift `ora_uncalibrated`.**
  - The statistician's flexible cross-fitted recalibration passes every cell of all 8 worlds at R = 3000 ($GIPPA_ROOT/stat/m804b/probe_flex.py). Its terms are 4 quantile hinges in logit q_raw, a one-hot TU size, a logit × log n interaction, and ridge 1e-4.
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
- Mutation evidence: $GIPPA_ROOT/review/pr187-mut/ (5/5 killed).
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
- `read_tsv` now treats an empty cell as missing in genes.tsv too (matches `gippa.io.genes_tsv`).
- Mutation evidence: $GIPPA_ROOT/review/pr191-mut/ (22/25 killed, 1 equivalent).

## M9-03 review nits (#210; one round, oriC rounding and the 2^20 gzip boundary pinned in the PR)
- Untested guards in `gippa.viz.support`: `window_bp < 1`, the missing-contig `ValueError`, the empty-sequence skip.
- `datapackage.json` lists `plots/support/gc_skew.tsv` uncompressed; a `.tsv.gz` write (not reachable at 10 kb windows for bacteria) would leave the descriptor pointing at a missing file — write the descriptor path from the file actually written.
- `plots` step: with `replicons` absent but a `plots_support` input passed, plots render from the incoming directory, not this run's output.
- No test checks that `write_datapackage` copies `schemas/plot_support/*.json` into OUT.
- ADR-0042: when the M7-10/M7-13 (`calsuite_variants`) and M9-22 (`calsuite_compare`) calibration tests land, extend scripts/slurm/calibration_full.sbatch to `--array=0-14` (task 13 = variants, 14 = compare); until then those two keys have no tests and the array skips them.
- M9-04 V17 review nits (#197): when every gene is partial, (e) draws 0 % bars and should show an empty-state note; integer tick locators on the count axes; partial genes are left out of (d) although the spec restricts only (e); the "unmapped" layer of (c) includes contaminant and too_short peptides; the length-bound lines belong at 6.5/40.5 (spec + test); an empty `replicons` draws an empty axis with a warning; the TSV message rows drop the drawn text; preflight `genome.contiguity` has an empty `message` (a preflight bug); "abundance not available" is not in the spec's titles; "1 preflight warnings" should be singular; the imputed legend shows when nothing is imputed; 8 surviving mutants (exactly 10 warnings, rule fallback, (b) floor, circular markers, constant-abundance padding, (e) inverted y, absent-genes title).

## M9-02 review nits (#205; one round, B1–B5 and N1 local half fixed in the PR)
- N1 Nextflow half: REPORT gets no `events`, so `layer_missing` skips fire under the local executor and `gippa report` but not in a Nextflow `gippa run`. Fix before the first Nextflow-run plot that declares layers.
- Config digest reads "not available" in `gippa run` but appears after `gippa report`, so `report.html` changes on a rebuild of an untouched run.
- Executor-equality test compares only `plots/data` (empty until plots register): also compare `report.html` bytes, and the run's report with a rebuilt one.
- Peptide "of N input" counts `kind=protein` rows from `unmapped`.
- Stale plot files from an earlier run are shown in place of a skip reason that holds now; with `formats=[png]` the report says "plot file absent" although the PNG exists.
- PLOTS has no resource label (1 CPU / 2 GB) while REPORT_HTML gets `process_low`.
- `ALL_PUBLISHED_KEYS` declares optional inputs never read (`peptide_hits`, `tus_tsv`); SVG id separator `-` is not provably unique.
- Surviving mutants: `<=` boundary in the card counts, the two `nunique` peptide counts. Move `cards()` to `gippa/report/cards.py` ahead of M9-16.
- Card/table read failures are notes in the HTML only (build.py is pure); consider an event from the step.
## M8-17a review nits (#204; one round, B1 R_LIBS_SITE fixed in the PR)
- 17b: `read_results` must require rows to match the in-range sets exactly (a dropped row now vanishes silently); `read_versions` should check the `seed`/`jobs` lines and the data.table version.
- Set `TMPDIR` to the call directory in `child_env` so R's `Rtmp*` is not left behind after a timeout kill; give the post-kill `proc.communicate()` a timeout.
- `write_call` accepts IDs with tab, newline or `;` (corrupts the R inputs): reject them.
- `scripts/build_envs.sh`: a failed lock check on an existing env leaves the `.gippa-built` marker; the marker hashes only `fgsea.yml`, so a lock-only change never rebuilds.
- `gippa_fgsea.R`: an empty input file exits 1, not 3.
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

### M7-12 statistician NITs (PR #225) → M7-13
- F1 matching folds I/L/J; precision over all truth sites (n_true_peptides counts database
  peptides only); use the method's own peptide rule as the N-terminal-extension denominator.
- Isobaric attack: a "systematic" error option keyed by protein site; run with
  `--genome-relation same`. Document `variant_f1` sample peptides labelled `search_db = same`.
- `same_strain_insilico`: all 11 non-ATG N-terminal peptides unmapped; the start-codon exception
  is never exercised.
- Three attack nulls proposed in the comment for M7-13 (tests/calibration/attacks/).
- GIPPA: genomes below ~200 genes hit a 1/n gene-level q floor above 0.01 and detect 0 genes;
  decide the small-genome behaviour (warn, or a method that can reach q ≤ 0.01).
- (M7-12 round 2) per-class expected false novel candidates in the truth; an I/L `ms_observable` column in alt_sites.

### M7-12 statistician round 2 (PR #225)
- **Possible cause of the sim100k 0-genes floor:** the older (non-six-frame) gippasim decoy rows
  carry no score/q/pep, so gene-level `decoy_picked` q is 1/T. Check this before the M9-13 evidence
  run and the q-floor question to Khoa.
- M7-13: the de novo arm of the class-FDR gate gets 0–2 false candidates per class per replicate at
  300k; run it at ≥ 1M.
- M7-13 G1: a policy is needed for `mass_equivalent` calls. On variant_f1 seed 7, high-confidence
  substitution precision is 0.82, and 5 of 6 FPs are A>S / N>D modification errors.
- Attack-null proposals in comment 6358601933.
- M7-12 N1 (statistician round 3): six-frame target q/PEP come from the wrong share before the coin (about 2× conservative). PSM-level calibration on six-frame scenarios needs the post-coin share.
- M7-12 N3: G4 for `nterm_extension` must go through M7-11's locus-tag pairing (planted peptides carry the gene's locus tag).
- M7-13: the report says a class with fewer than about 100 candidates cannot be called at α 0.01 (class q floor).
### M6-09 (from the M6-03 statistician review, PR #222)
- **Method ADR needed (Khoa):** the SW decoy floor (≈2 SW decoys vs 131 targets in the 50–150 bp
  bin) keeps min q ≥ 0.88 under the merged decoy-FDR method, so the missed-operon FDR gate passes
  vacuously. Options: one phantom decoy in total, a CV-only or reweighted-maximum rule, coarser
  bins; otherwise document the gate as "vacuous pass: no calls". Probes:
  statistician scripts probe.py / floor.py / teeth.py (lift into tests/calibration/attacks/).
- Record one arm at the defaults (`p_same_strand` 0.6, `DEFAULT_BETA`); the power table reports
  the realised planted count (≈15 % of both-detected NO targets, not 10 %).
- (M6-09, statistician round 2) example 1 is infeasible at 0.1 SD under every arm (merged engine: 0
  calls everywhere; CV-only 1-bin: 2.5 calls); add a realistic 10 % detection-attack block arm;
  report realised planted counts; lift stat603/probe2.py and bias.py into tests/calibration/attacks/.

### M6-03 code-review nits (PR #222)
- Decision 9 (odd `n_planted`) untested: `discordance_power(1800, 41, 2.0)` gives 21/20.
- Untested edges: SW triple across the origin; convergent boundary at a second replicon's first gene;
  2-gene circular replicon; the inclusive 300-bp upper gap bound; background SD from both-detected pairs only.
- The circular wrap boundary never gets the short convergent gap, though layout.py and the methods
  note say "every" boundary: fix the wording or the generator.
- Methods note: the R = 1,000 null OP tail z = −2.05 needs a "within chance" remark.
- `_tsv_cell` / `read_operon_truth` `"int"` dtype branch unused.
- mypy over tests/sim (MYPYPATH=tests/sim) reports 5 errors from tests/sim/gippasim/layout.py:663.
### M8 Wallenius (statistician, PR #220)
- **ADR-0044 (Needs Khoa):** M8-13 criterion 1 reference = exact Wallenius (urn recursion /
  exact rational integral, committed TSV + pure-Python generator); goseq second reference at
  `golden_r.wallenius_goseq_abs_tol: 1.0e-5`, replacing `wallenius_abs_floor`.
- **ADR-0045 (lead recommends):** GIPPA far tails (< 1e-6) by the exact urn recursion; scipy loses
  relative accuracy below ~1e-9 and returns 0 below ~1e-14 (40-urn spot check, stat813/scipy_far.py).
  Hardening pass; does not block M8-13.
- G3b NES across seeds: reported as per-set z, alert at max |z| > 4.5; p rule alert at ≥ 3
  exceedances; alerts open a `stats` issue.
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
- `gippa run --strict` does not map `report.size_exceeded` to exit 3 (only preflight warnings reach
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
  took 47 min on 8 CPUs) and needs `gippa-fgsea` built; `tests/calibration/_criteria.py` does not
  exclude `requires_r`, so the harness's `gsea` suite errors where the env is absent.
- Untested in `gippa.function.gsea`: dropping `lo <= n <= hi` for whole-ranking sets; the
  `len(tus) < min_size` and `bh(sub)[i] > alpha` boundaries; `SIMPLE_ENGINE` re-seeding per job.
- `gippa.io.fgsea_r`: an explicit `GIPPA_FGSEA_PREFIX` skips the `.gippa-built` check (an unbuilt env
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
- Island power under gippasim's default `covariates_only` detection: no island at +3/+5/+10 SD on 1,000 genes (ℓ̂ = 1, ℓ' = 3; ~45 % censored regardless of abundance, rank scores saturate). The M5-11 golden uses sim4500 with `detection="abundance"`; decide whether the power table should report the covariates-only regime.

## M5-11 review nits (#217; one round)
- The `gippa-bench` sbatch script used for the `context_m5` benchmark (Slurm 1869098) is not in
  scripts/slurm/; commit it beside the other benchmark launchers.
- M5 gate evidence: full `islands` records exist (Slurm 1868687, R = 500, all pass); the full
  `autocorrelation` record does not — submit `calibration_full.sbatch --array=4` on main. Calibration
  records carry no commit field, so "on main" cannot be checked from them.
- The island golden's join-count rows depend on the lbfgs detectability fit (decision 44); a
  BLAS-stable `p_detect` would let them return to byte equality.

### M9-16b follow-ups
- Reword "run_manifest.json not written by this GIPPA version" once the tail re-render lands (only
  shown when the manifest write failed or in the work-dir copy); test-author step first.
- Schedule MVP-03c (tools[] / databases[] in the manifest) before the M9-25 approval packet.
### M9 pass 2 (from M9-16 tests)
- §7.8.2's separate ribosomal-gene report table: no M9-16 behaviour; assign to a pass-2 issue.
- M9-16 decisions 5 and 7 (ALT card with mixed `unestimated` rows; `provenance_digest` line) are
  implemented but untested: add tests in the hardening pass.

### M9-16 review nits (PR #228)
- Untested (surviving mutants): `call == "ALT"` filter (add a `mixed` row), "(N unestimated)" tail,
  `provenance_digest` line, `git_dirty` → "(uncommitted changes)", sensitivity "both zero" when one
  is zero, last vs first `operons.missed_operons` event, whitespace-only estimate in
  `unlabelled_numbers`, KEGG `retrieved` vs `version` (fixture gives both the same date).
- `_diagnostics` crashes the report via `render_l_hat` on an unknown/missing `block_rule`: note it.
- Function section loses the KEGG caveat when `category_summary` is absent but GSEA is present.
- `genomic_cluster` flag matched as a substring; split on `;`.
- Wording scanner: attribute text (`title`, `alt`), words split by inline tags, unclosed `<q>`.
- Sync script: missing end marker, second block, unsupported Markdown, `--docs` to a missing dir
  without `--check` deletes every snippet.
- `content.snippets()` reads package data from report/ (like templates): documented exception.
### M7-11 follow-ups
- M7-11b: novel_orfs step wiring of `gippa.novel.fdr` (decisions 6–9), step-level tests.
- Origin-crossing makedb entries without replicon lengths are unestimated (decision 8); pass lengths.
- M7-11 statistician NIT 1: c is computed over hit peptides only; decide whether to widen it (pending lead/statistician, not blocking).
- M7-11 vs M4-24: when only some targets lack q, M7-11 falls back to score while M4-24 keeps q and ranks them last; align the two (statistician #237 round 2, safe).
- Review #237 nits (M7-11):
  - 1: a class with no seen entries reports `paired_share_low` even at f = 1, where it should be `no_target_score`.
  - 2: hoist the loop-invariant `seen` (fdr.py:294).
  - 3: without lengths, entries whose stop codon lies past the origin are unmatched; name them in the methods note.
  - 4: `min_identity` / `min_coverage` have no lower bound, so 0/0 gives f = 0 and `run_detection` raises; this needs a config bound (protected, ADR) or a guard.
  - 5: the methods-note status line says decisions 1–11; it should say 1–13.
  - 6: `_decoy_tags` could return the parsed class.
  - 7: `denovo_class_q` uses an iloc loop and does not cast strand to str.
### M9-24 follow-ups
- **Diagnosed:** V08 (and V02, V03, V07) show no measured abundance on sim100k because 0 of 87
  genes are detected: picked-protein q has floor 1/n_win = 1/85 > alpha 0.01
  (genetable/detection.py:344-349; #186 added the `q_floor_above_alpha` warning, q unchanged).
  Abundance input is intact. The sim100k end-to-end fixture therefore never exercises any
  detected-gene branch. Options (Khoa): a larger evidence genome (> 100 eligible genes), a
  sim100k alpha override, or a two-peptide fallback when the floor exceeds alpha. V08's title
  could cite `q_floor_above_alpha` from run flags.
- The checklist reads tick labels through matplotlib's private `Axis._update_ticks()`; replace when
  matplotlib offers a public way.

### M9-24 review nits (PR #230)
- N1 any axhline/axvline is exempt from the palette check in any colour (a red axvline passes);
  `style.RC` defines no greys, so restrict the exemption to grey/black reference lines and test it.
- N2 COG + detection check: palettes share #0072B2/#E69F00 → false fail and false pass cases.
- N3 units allow-list too permissive (generic suffixes "rate", "density", "fraction"; any
  parenthetical counts as a unit).
- N4 overlap: clipped texts outside the view false-fail; text off the figure edge undetected;
  rotated text uses an axis-aligned box (false fails at 30–45°).
- N5 surviving mutants: `_update_ticks` vs `get_major_ticks`, axis labels in the text list,
  collections by offsets only, `ax.texts` in the count, each half of `is_empty_state`, `pid_`
  prefix in `_tsv`, `score()` with nothing gated.
- N6 `check_instance` rebinds `inst.figure.canvas`; `_tsv` calls `Path.is_file()` from viz/.
- N7 the snapshot sweep skips until M9-25 fixtures; sim100k exercises only no-abundance branches.
- N8 pin matplotlib in pyproject.toml while the checklist uses `Axis._update_ticks()`.
### M9-15 follow-ups
- V11 twin bars at bin index; write `bin_lo`/`bin_hi` into the V11 data TSV so the twin matches
  the static edges.
- Tests for M9-15 decision 5 (`meta`: α lines, kb/Mb choice, origin-crossing islands).

### M9-15 code-review nits (PR #234, comment 6358615453)
- Decision 5 (`meta`) untested: no test fails when `figure_json` never passes meta. 12 of 15 manual
  mutants survive (V2 1-bp offsets hidden by tolerance, origin split, meta dispatch, `columns=`,
  `html.unescape`, V7 `lo == hi`, V10/V11 alpha, V8 n = 50 boundary, FASTA order, V18 bar offset).
- `inspect` dispatch counts keyword-only parameters and `*args`, so `(data, title, *, extra=None)`
  raises TypeError; decide at registration instead.
- `columns=`: a twin whose columns drift silently disappears; log a warning.
- V11 (a) x-axis still says "|Δ| …" though the bars sit at the bin index.
- V7: rows with the same label collapse onto one category, and some CI bars are hidden that the static draws.
- V5 `arrangement="fixed"` without node x/y loses the M9-05 order.
- V8 has no median mark or legend; the V10 no-meta fallback differs for flagged rows; V18's constant state has no message.
- PR size about 2,300 lines in total (711 counted).
### M7-10 follow-ups
- Switch the class-error calibration to the M7-12 `variant_f1` scenario once #225 merges (decision 7).
- M7-13: class error rate reference-proteome arm (moved from M7-10 #239 B2); also consider a run flag when every class rate is 1 (statistician #239).
- M7-10: `paired_observations` keeps the last reference row per spectrum (`keep="last"`), so raw callers with conflicting rows depend on row order (pipeline safe: preflight drops conflicts). Also a run flag when every class rate is 1, and a `low_n` flag or interval columns (statistician #239).
- Review #239 nits (M7-10):
  - 1: `_peptide_spectra` swallows reference `ValueError`s and validate issues, so a malformed reference file is silently `unestimated`.
  - 2: edge-case tests are missing for a frameshift with `mass_equivalent`, a raw NaN reference q, NA `is_decoy`, and NA `spectrum_id` on both sides.
  - 3: the methods note "Observations" bullet should say distinct `(run, spectrum_id, peptide)`.
  - 4: `site_class` ignores `;`-joined string flags cells; document this or split them.
  - 5: the `# fmt: skip` on the cli.py `files` dict is unneeded.
### M9-13 review nits (PR #224)
- N1 checks.py test gaps (surviving mutants): `<script defer src=…>`, `url(http…)` outside style,
  uppercase `DATA:` URI, query-only `href="?q=1"`.
- N2 external loads not detected: `srcset`, `<object data>`, `<video poster>`, SVG `<feImage href>`,
  `<base href>`, meta refresh — extend before M9-14/M9-25 rely on checks.py.
- N3 `relative_links` treats root-absolute paths as relative (docstring); tab/newline URL cleanup.
- N4 examples 2–3 checked on the local report only, not the Nextflow-published one.
- N5 untested: Nextflow run with a module off gives `module_off`; main.nf passes `events` to PLOTS.
- N6 main.nf repeats the `function` condition for FUNCTION's events; fold into the existing `if`.
- N7 test_branch.sbatch warm-up ≥ 25 s should warn loudly (30 s probe in tests/conftest.py).
### M9-22 hardening
- ρ-CI calibration on full gippasim run pairs (now a simulated TU world at gippasim sizes).
- Test the section order (matrix, then intersections) and the `none` row's place among ties.

### M9-22 statistician nits (PR #232)
- TUs from the pair's first run only; a union of both runs' TU maps is safer.
- compare.md: Bland–Altman limits are as uninterpretable as the bias under `scale_differs`;
  optionally a TU-block CI for the BA bias.
- `enriched_categories` deflates Jaccard when only one run has GSEA; GSEA `nes == 0` maps to
  "down"; the operon_flags key ignores the `flag` column.
- Unit tests for `ci_unestimated` and NA-TU genes as their own unit; recorded heavy-tailed and
  15-TU calibration arms; attack nulls A1–A5 (comment) in tests/calibration/attacks/.

### M9-22 code-review nits (PR #232, comment 6358611226)
- Untested (surviving mutants):
  - the 9-unit boundary of decision 11;
  - decision 7 (TUs from the wrong run);
  - genes without a TU as their own units;
  - strand in the ALT key;
  - the `end`/`strand` gene-order tie-breaks;
  - same-size, different locus-tag sets;
  - `abundance_type == "none"`, its sort position and its zero-count row.
- `test_function_off_in_one_run_is_missing_in` deletes the tables, so the `resolved_config.yaml`
  module check is never exercised; a run with only one of its two enrichment tables is untested.
- A malformed manifest/preflight JSON or an unreadable stability table gives a traceback, not an
  exit code (`gippa report` exits 4).
- A manifest without a genome sha256 should fall back to preflight; otherwise say so in ADR-0046.
- The `n_units` schema description says NA when there is no CI, but the code writes 1 alongside an
  NA CI. Fix the description.
- One shared stream across pairs: adding a run shifts later pairs' CIs. Spawn per pair.
- GSEA NES 0/NaN counted as "down".
- `gene_order` duplicates the §6.5 key.
- No manifest for `CMP/`.
- compare.md still says "Statistician sign-off: pending".
### MVP-03c follow-ups
- `gippa _step manifest` passes no db manifest, so `databases` is `[]` on the Nextflow step path (decision 20).

### MVP-03c code-review nits (PR #236, comment 6358790241)
- Decision 17 `mapping = null` for stub events is untested (the test accepts either).
- `_tools`: the "first non-null in trace order" merge can give one row both a digest and a lock hash; the digest winner depends on completion order.
- A tag-only container still gets a lock hash; `*.pip.txt` locks are not hashed.
- A calibration key without a `calibrated` field shows `false` but adds no `uncalibrated:` flag.
- `databases`:
  - a non-string provisional version fails the schema;
  - a top-level `schema_version` becomes a row;
  - a non-dict `items` entry raises;
  - `retrieved: ""` should be null.
- No unit test writes a 0-byte output, so the 0-byte rule is covered only by test_mapping_nf_run_events.
- The plots step's `oric_position` events are labelled `step: context_stats`.
- `_plain` maps ±∞ but not NaN.
- Schema: `oric` items lack `additionalProperties: false`; `bin_counts` items accept any numeric key.
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
- Reword "run_manifest.json not written by this GIPPA version" once the tail re-render lands (only
  shown when the manifest write failed or in the work-dir copy); test-author step first.
- Schedule MVP-03c (tools[] / databases[] in the manifest) before the M9-25 approval packet.
- M8-15: a synthetic go_slim OBO for the sixth collection (M8-14 gap 4).

### M4-24 follow-up (gippasim scored null decoys)
- gippasim M3-19 proteomics writes decoy rows with q, PEP and score all NA, so the picked tier never
  sees decoy competition (D = 0 above every target; after M4-24 such runs fall back to `mayu`). Give
  decoys null scores: each null PSM lands on a target or a decoy with probability ½ (decoy count
  Binomial(2 · n_wrong_rows, ½)), with score, PEP and q drawn from the null component like the wrong
  targets, and add null-only target genes (genes whose only evidence is null PSMs), so calibration
  exercises decoy wins on genes that also have a target. Evidence: Slurm 1869554 and
  local/scratch/qfloor/ (patched proteomics.py: `decoy_score_na` false, decoy wins still only on
  target-less genes because no gene is null-only). Then restore `decoy_picked` on the sim100k_real
  golden and in tests/calibration/test_gene_fdr_real_mapping.py (M4-24 question 3).
- M4-24 decision 7: preflight `candidate_tier` should use the shared key columns, so it does not
  preview `decoy_picked` for all-unscored decoys.
- M9-23: switch existing plot titles to `gippa.viz.base.plural` as their tests are next touched.
- M9-23: hide V16's `Skip` (`input_absent:compare`) from ordinary run reports; add a test for the runs-section detection tier per run.
- Review #225 nits (M7-12):
  - N2: the methods-note FDP table is at 200 ORFs per class; rerun it at the gate count of 300.
  - N3: `OrfEntries` does not restate `NTERM_MIN_UPSTREAM_AA` and assumes `min_orf_aa` = 7; assert `digest_len[0]` ≥ 7 in six-frame mode.
  - N4: drop the throwaway `Builder(10)` in the mixed-decoy test.
  - N5: validate that `n_per_class` ≥ 0.
  - N6: add a 1M gate-count case to the benchmark and refresh the `novel_orf_power` baseline.
  - N7: `_protein_decoys` checks the target set without I/L folding.
  - N8: the spec names `tests/sim/test_m7_scenarios.py`; the file is `tests/unit/test_sim_m7_scenarios.py`.
- MVP-04: origin-crossing novel ORFs excluded from the planted set on K-12 (spec question 7).
### M7-12 statistician NITs (PR #225) → M7-13
- F1 matching folds I/L/J; precision over all truth sites (n_true_peptides counts database
  peptides only); use the method's own peptide rule as the N-terminal-extension denominator.
- Isobaric attack: a "systematic" error option keyed by protein site; run with
  `--genome-relation same`. Document `variant_f1` sample peptides labelled `search_db = same`.
- `same_strain_insilico`: all 11 non-ATG N-terminal peptides unmapped; the start-codon exception
  is never exercised.
- Three attack nulls proposed in the comment for M7-13 (tests/calibration/attacks/).
- GIPPA: genomes below ~200 genes hit a 1/n gene-level q floor above 0.01 and detect 0 genes;
  decide the small-genome behaviour (warn, or a method that can reach q ≤ 0.01).
- (M7-12 round 2) per-class expected false novel candidates in the truth; an I/L `ms_observable` column in alt_sites.

### M7-12 statistician round 2 (PR #225)
- **Possible cause of the sim100k 0-genes floor:** the older (non-six-frame) gippasim decoy rows
  carry no score/q/pep, so gene-level `decoy_picked` q is 1/T. Check this before the M9-13 evidence
  run and the q-floor question to Khoa.
- M7-13: the de novo arm of the class-FDR gate gets 0–2 false candidates per class per replicate at
  300k; run it at ≥ 1M.
- M7-13 G1: a policy is needed for `mass_equivalent` calls. On variant_f1 seed 7, high-confidence
  substitution precision is 0.82, and 5 of 6 FPs are A>S / N>D modification errors.
- Attack-null proposals in comment 6358601933.
- M7-12 N1 (statistician round 3): six-frame target q/PEP come from the wrong share before the coin (about 2× conservative). PSM-level calibration on six-frame scenarios needs the post-coin share.
- M7-12 N3: G4 for `nterm_extension` must go through M7-11's locus-tag pairing (planted peptides carry the gene's locus tag).
- M7-13: the report says a class with fewer than about 100 candidates cannot be called at α 0.01 (class q floor).
- Review #225 nits (M7-12):
  - N2: the methods-note FDP table is at 200 ORFs per class; rerun it at the gate count of 300.
  - N3: `OrfEntries` does not restate `NTERM_MIN_UPSTREAM_AA` and assumes `min_orf_aa` = 7; assert `digest_len[0]` ≥ 7 in six-frame mode.
  - N4: drop the throwaway `Builder(10)` in the mixed-decoy test.
  - N5: validate that `n_per_class` ≥ 0.
  - N6: add a 1M gate-count case to the benchmark and refresh the `novel_orf_power` baseline.
  - N7: `_protein_decoys` checks the target set without I/L folding.
  - N8: the spec names `tests/sim/test_m7_scenarios.py`; the file is `tests/unit/test_sim_m7_scenarios.py`.
- MVP-04: origin-crossing novel ORFs excluded from the planted set on K-12 (spec question 7).

### M6-10 follow-ups
- Runtime: build the pair frame once per threshold and keep reruns under 2× a single run; not pinned (timing is flaky on the shared node).
- Review #243 nits:
  - N1: `gippa.context.sensitivity` imports `gippa.steps.operons`, the first context→steps import. Move `run_operons` to a pure module and `read_layer` to io/.
  - N2: the "relative" log has no `extra` event, and its wording fits only this family.
  - N3: no unit or step test passes precomputed TSV layers (mutant M9).
  - N4: hardening still owed: `--sensitivity full` in the nightly e2e, and the M9 text item on threshold sensitivity.
  - N5: `test_unpatched_variants_unchanged` has an odd `assert [c for c in ...]`.
- M7-11b: no run input carries makedb's `db_id_map.tsv`, so `nterm_extension` pairing comes from the run's gene model and the 80 % pairing guard cannot fire from the pipeline. A real input needs a new `defaults.yaml` key (protected, ADR, Khoa).
- M8-15 (from M8-14 statistician, #251): build ORA sets through M8-03 (world.sets keeps 10–500 TUs, ORA tests 5–500); keep the confounded set's GO ancestors out of level checks and null pools; compare the truth's dominant block to M8-09 by signed argmax on top-decile genes; wire planted features into gene-table columns `promoter_score`, `tir_ostir`, `cai` (simulator CAI uses the bundled codon usage); record a position-decay 0.5 arm on the R engine; gippasim `resid` should read `ecfg.tu_sd` (hard-coded 1).
- M7-14 review nits (#248): a step test mixing same-db and other support across sites (alt_peps order, variants.py:220); te2_result `notes` empty when high_stop_codon > 0; merge the duplicate `_flag_set`/`_flag_names` helpers (drop the unreachable `;` branch); `variants_qc` in `PUBLISHED_FILES`/`run_view.PATHS` though never published; test_id_mixups_sim100k returns early on missing inputs (make it skip); an nf-test with a non-empty `variants_qc` into CONTEXT.
- M7 (from MVP-04 statistician, #247): nterm-extension decoys pair with a gene only after a target row is accepted (M7-12 decision 28), so an all-decoy null gene vanishes instead of counting as a decoy win (anti-conservative, ≈ +0.8 false genes at K-12). Pair through makedb's entry→gene relation.
- M4 (from M4-24 statistician, #238): ADR-0047 q-filter symmetry: switch the decoy q filter on only when *every* C9 target carries q (with 50 % NA-q targets FDP at q 0.05 is bimodal, top decile ≈ 0.135; Slurm 1871630). Needs an ADR-0047 amendment (Khoa).
- MVP-04 review nits (#247): the scored_null collision rule is tested only vacuously (0 collisions everywhere): force a palindromic null; `gate7` writes NaN (invalid JSON) when nothing is detected; `_publish` merges into an existing outdir (stale files survive); `_read_fasta` KeyErrors on sequence before the first header.
- M4 (from MVP-04 part 2, #258): gene_table logs the detection warnings `decoy_pairing_low`, `decoy_unscored`, `q_floor_above_alpha` as plain text with no event or run flag, so run_manifest.json cannot carry them; log them with an event (and decide flag tokens, ADR-0041).
- M9-23: compare.html has no size warning (the 25 MB check covers report.html only); the V16 TSV/twin carry every shared gene per pair. Add the check or bin the twin past N points.
- M9-23: V16's skip on ordinary runs happens inside the plot function, so `skip_reason()` is None and report.html's skipped notes omit it; move it to a registered skip reason.
- M9-23 review nits (#254): PNG-only output says "plot file absent" on compare.html; binned V16 cells have per-cell colour scales and no count legend; V16 upper-cell scatters have no axis values; a locus tag literally `NA` is misread on the twin's TSV round trip; drop `NA_CELL` (= GREY) and `_mappable`.
- (M7-11b, statistician #256 note 2) Floor warning for picked classes quotes `n_targets`; quote `n_target_wins` (the floor's denominator). Implementer follow-up.
- (M7-11b, statistician #256 note 1) Fallback calls (`no_makedb_pair`, `no_target_score`, `decoys_unscored`, untagged) have no error control; M7-13 reports the per-class fallback count.
- (MVP-04 gate 13, statistician) Multiple strong islands mask each other under the stationary block-bootstrap maxT null (islands inflate ℓ̂, ℓ′ ≥ w copies them into M_b). Candidates: ℓ̂ with top windows masked, or a step-down scan. Needs an ADR-0021 amendment + calibration.
- (M7-11b review N1) `class_call.py:133-137` re-derives `fdr_class` from private `fdr` helpers; have `PickedClassFdr` return the entry per ORF.
- (M7-11b review N2/N3) `reason` typed `str | None` (spec: `Reason | None`); redundant guard at `class_call.py:145`.
- (M7-11b review N4, next to Q11) `best_costs` keeps the first row's `source` for duplicate peptides while `orf_regimes` uses the union; a db+de novo peptide without makedb decoys can end `no_target_hit` (conservative).
- (M7-11b review N5, lead) close backlog line 777 (M7-11b), re-point 688, add the Q10 V13 floor-caveat line.
- (M7-11b TA, decision 9) `no_target_score` is unreachable from the novel step (orf_regimes and picked_class_fdr read the same target rows); test forces it via monkeypatch. Either drop the reason or document it as library-only.
- (M8-14 decision 6) packed set with < 10 scoring TUs could trip the finite-p assert in `_packed_tu_rep` (~1e-5/set) at large full-tier R; consider raising PACKED_SINGLES.
- (M8-14 decision 7, for M8-15) gippasim pseudogenes without `aa` record n_obs=0 / GRAVY NaN while GIPPA translates them from `nt`; M8-15's p_detect refit sees different covariates for them.
- (MVP-04 #258 review) boundary tests for summary gates 4, 6, 8, 11, 13 (p) and a gate-13 no-high-island case; sbatch seed from `K12_EVIDENCE_SEED` as `MVP04_SEED`, tested as passed through; gate 15 counts the `empty_state` SVG as content (N2); N3–N9 tidying in `scripts/dev/mvp04_summary.py`.
- (M2-14b item 7) `PC_BATCH = 100` (regions per PROMOTER_CALCULATOR task) is a constant in annotate.nf; make it a config key (ADR, defaults.yaml) if a genome needs a different batch size.

### M7-17 follow-ups (PR #291)
- (Q6) The flagged share of true ORFs is fixed by construction in the two existing arms (0 picked, ≈ 1 de novo). Add a mixed database + de novo input arm (M7-13c or MVP-04) to measure it.
- (Q8) `orf_regimes` treats an unresolved `auto` C9 source as non-`denovo`, while mapping and the flag resolve it table-wide (`resolve_sources`); a source-less de novo table gives `unestimated` ORFs instead of `denovo_subgroup`. M7-11 follow-up: resolve sources the same way (this changes regimes and q, so it needs its own issue).
- (Decision 4) `test_denovo_support_flag_gate.py` reruns its own G4-k3 arm (one extra simulation per replicate); reuse the G4-k3 runs.
- (#291 review nits)
  - report card `sources` misses `novel_orfs:flags` (build.py:304);
  - `_FLAG` is defined after its first use in steps/novel.py;
  - `rbs is True` in evidence.py:60 should be `rbs is not None and bool(rbs)`;
  - the Hypothesis ORF-side alphabet includes `I`;
  - the gate test should count `side_flag != tsv_flag` as a violation.
- (M7-17 gate, statistician 7f8e4ed0) `inject_false_orfs` should skip segments at truth keys, as `null_segments` does: about 0.5 % of G4-k3 injections land on a planted ORF's stop and add no false ORF. Makes G4-k3 slightly stricter and changes its stream (rerun); needs Khoa's batch. Do not exclude every CDS-overlapping segment (replicate 130 is correct behaviour).

### 2026-10-24 review nits (deferred)
- (#297 review nit 2) Size: about 1,280 test lines against the ~400 guide; next time, ship the Wallenius half (`test_golden_wallenius.py` plus most of `tests/helpers/golden_r_gate.py`) as its own PR (process note, no code change).
- (#297 review nit 3) Empty leading edges: `GoldenReplayEngine` maps an empty `leadingEdge` cell to `()`, but `_leading_edges_in_fgsea_order`, `test_es_oracle_leading_edges_match_golden` and `_g3b` (`tests/golden/test_golden_gsea_wiring.py:225`, `:299`, `:585`) call `.split(";")` and get `[""]`; add a shared `_edge(cell)` helper and exercise the oracle's `top == -bottom → ES 0` branch.
- (#297 review nit 4) `es_oracle_pass` (`tests/golden/test_golden_gsea_wiring.py:605`) is `not any("oracle" in b ...)`, so leading-edge failures against the oracle also clear it; rename it `oracle_pass` or match only the ES message.
- (#297 review nit 5) `decade()` (`tests/helpers/golden_r_gate.py:100`) has no direct test (an off-by-one mutant survived); add a three-case parametrised test (0.5, 0.05, 0).
- (#297 review nit 6) `g1_pass`/`g2_pass` re-run the G1/G2 test functions through `_passes` (`tests/golden/test_golden_gsea_wiring.py:382-395`), and `g2_pass` leaves out `test_g2_written_table_validates` and the NA-flag test; record outcomes in a `pytest_runtest_logreport` hook instead.
- (#297 review nit 7) Dead field: `GoldenReplayEngine.calls` (`tests/helpers/golden_gsea_fixture.py:122`) is appended to and never read.
- (#297 review nit 8) The Slurm test `test_example9_uncentred_scores_fail_g3b` (`tests/golden/test_golden_gsea_wiring.py:659`) repeats the CI test `test_example9_uncentred_scores_fail_the_g3b_floor` at about 2 min of live R in the weekly task; keep it for its end-to-end floor assertion, but weigh the runtime.
- (#299 review nit 3) `tests/unit/determinants/test_determinants_sensitivity_agreement_none.py`: the guard message at line 152 ("fixture: no explained gene") and the docstring at line 137 ("seed 3 has a top-decile gene without a positive block") are stale under the |SHAP| rule; say in the docstring that "None in both runs is excluded" is covered only by the hand-built tests.
- (#299 review nit 5) `determinants_family` looks up `sums.loc[t]` once per row (`src/gippa/determinants/sensitivity.py:105`); vectorise it with `np.abs(sums.to_numpy())` and a row-wise argmax, and assert or document that `agreement_block` expects `len(sums) == len(BLOCKS)`.
- (#300 review nit 14) Trap timing in `scripts/slurm/benchmark_standard.sbatch`: bash defers the TERM trap (`:66`) until the foreground `srun` or `sleep 60` (`:169`) exits, which can outlast KillWait; use `sleep 60 & wait $!`, `srun … & wait $!` with the trap forwarding to the srun PID, or `#SBATCH --signal=B:USR1@600`, and stop a TERM during the sacct retries from marking an already finished run `timeout` (the `t0`-before-`cur` part is already fixed at `:190`).
- (#300 review nit 19) `snap()` quoting: `sed "s|^|$d/|"` breaks on `|` or `&` in a cache path; the merged `snap()` (`scripts/slurm/benchmark_standard.sbatch:149`) already uses a `printf '%s/%s\n'` loop, so check and close.
- (#300 review nit 21) `summary.md` shows timeouts as "timeout (exit 124)" (`scripts/dev/benchmark_summary.py:242`), while the spec says they are reported as "> 5 h".
- (#300 review nit 22) Before the 36-hour benchmark array, submit one short smoke cell (`BENCH_ONLY=sim5m/core`, a small `BENCH_RUN_TIMEOUT`, a scratch `BENCH_ID`) to confirm that an `srun --mem=64000M` step starts next to the batch step under `CR_CORE_MEMORY` and that step IDs start at 0.
- (#298 review nit 3) The GBFF sibling lookup (`workflow/subworkflows/local/annotate.nf:79-81`) only finds lower-case `.gbff`/`.gbff.gz`, so `x.JSON.GZ` next to `x.GBFF` skips antiSMASH with `layer_missing` on a case-sensitive filesystem; say so in the warning or docs, or add `.GBFF`/`.GBFF.GZ` candidates.
- (#298 review nit 4) The `--executor local` test (`tests/unit/test_cli_run_annotation_json_gz.py:197`) stops at `_local_roots`; a `stage_roots` + `parse_bakta` round trip on a gzipped root would close the loop (optional; gzip is detected by magic bytes).
- (M10-06 review N4) `_diagnose` (`validation/pilot/te7.py:122-123` on m10/spec-06) returns `content` with `first_differing_line: None` when the bytes differ but `splitlines()` is equal (CRLF vs LF, trailing newline); give a line number or a distinct reason.
- (M10-06 review N5) `_numeric` (`validation/pilot/te7.py:58` on m10/spec-06) counts a NaN-vs-number cell in `n_cells` but adds nothing to `max_abs`/`max_rel`; report an `n_nan_mismatch` count or document it.
- (M10-06 review N8) `_tree` (`validation/pilot/te7.py:263` on m10/spec-06) walks all of `work/` (millions of Nextflow task files) only to count it as excluded; prune the excluded top-level dirs, or note the cost.
- (M10-06 review N9) `_diagnose` decompresses and `splitlines()` a whole differing gate file (`validation/pilot/te7.py:102-109` on m10/spec-06), so a multi-GB pileup could exhaust the 8 GB `dev` job; cap the size or stream.
- (M10-06 review N10) `metrics.py compute --te-real/--te7` (`validation/pilot/metrics.py:1155-1157` on m10/spec-06) does not check `schema` (`gippa-pilot-te-real/1` / `gippa-pilot-te7/1`), so a swapped file passes silently.
- (M10-05) The unit fixture `tests/helpers/pilot_metrics_m10.py` writes `peptide_hits.tsv`/`unmapped.tsv` under `tables/` (lines 371-378 on m10/spec-05), but GIPPA writes them under `mapping/`, and `validation/pilot/metrics.py:932` reads both; move the fixture to `mapping/` and drop the `tables/` read.
- (Slurm-only failures, jobs 1879953, 1880533, 1881053) `tests/integration/test_nf_review.py::test_farmshare_local_user_limits_vs_allocation` (`:243`) and `::test_farmshare_user_limits_below_profile_win` (`:233`) exit 137 on gippa-test nodes but pass in CI; find out why.
- (M10-08 note) Cascadia 0.0.7 fetches Unimod over plain HTTP twice per run (pyteomics at the end, psims `GenericModification(...).mass` at the start), and the offline launcher `scripts/pilot/cascadia_offline.py` (#301) now blocks both; upstream (Noble-Lab/cascadia issues #6 and #12) may want to know. Note only: no one has been contacted.
- (M10-09 review F8) `float()` in the cascadia score parser also accepts `9_0e-2` and padded values (` 0.9 `); tighten it to a decimal regex if it ever matters.
- (M10-09 review F9) one `.ssl` file carrying two different `file` values becomes two runs with no note (W2 says one value per file); consider a warning.
- (M10-09 TA) the validator's `non_letters` message says "(dropped)", but the row is written with `peptide` missing; align the wording.
- (M10-10 spec finding 1) `mzml_pilot.sbatch` asks `--mem=32G` and `denovo_pilot.sbatch` asks `--mem=64G` with 8 CPUs; at MaxMemPerCPU=4000 Slurm raises them to 9 and 17 CPUs. Use `--mem=32000M` and size denovo from the smoke MaxRSS (M10-10 overrides both on the command line).
- (M10-10 spec finding 3) the plus-strand `orf_id` depends on where the peptides fall, so compare's `novel_orfs` Jaccard can count one ORF as two across arms; key on the stop codon (as M10-10 metric 6 does).
- (M10-10 statistician S6) add a per-bin n column to `mapping_fdr.tsv` (M3 owner), so the de novo bins can be cross-checked directly.

### 2026-10-24 afternoon (lead)
- (M6-13 statistician) a variance-floor σ model σ² = σ_f² + σ₁²·e^(2γz), gated on D-1b-floor, D-1c and the null quartile check cut on a. It lowers the D-1c top-quartile excess from 1.7× to about 1.4×, at a 1.03–1.10× level cost in the null a-quartiles.
- (M6-13 Q1) P(miss) from the scaled fit's rescaled residuals (ADR-0055 Consequences). Measure with attack C under D-1 and D-1c noise.
- (#310 review N-4) no test checks that gippasim's own pair table agrees with `build_pairs` (classes, distances, wrap-around pairs).
- (M10-07 review) N1, the methods-note sentence; N3, the `outside_bounds` count. `metrics.py`'s `_outside_repo` stays strict, because main's M10-05 tests pin exit 2 for a git-ignored path; revisit whether it should accept git-ignored paths, as entrapment.py does.
- (M10-08b review N10) timsconvert's lock pulls in `pyBaf2Sql`, which bundles Bruker's baf2sql library. ADR-0054 covers only the TDF-SDK. Confirm the licence note, or record that baf2sql is unused for .d/TDF input.
- (Khoa, 2026-10-24) running GIPPA on macOS is untested. The Python package passes on macOS CI; the full pipeline would need `-profile docker`, which spec §3 keeps but doesn't test. Candidate for M11: a Docker-profile test run, and documenting the platforms.

### 2026-10-24 evening (lead)
- (M9-25, Khoa's approval) Two checklist failures were left: V16 uses #555555, which is not in the §10.1 palette, and on K-12 V15 panel (c)'s title overlaps the CV R² line. Fix both and re-baseline.
- (M7-16a statistician) Relabel or retire the linear local Nissle Bakta run `local/recordings/bakta_full/nissle1917`; the committed recording `tests/data/recorded/bakta/1.12.2-full/nissle` is circular.
- (M7-16a statistician, #331) At 1 Mb, `related_other_strain` accepts 0–5 true read-through rows per replicate. Confirm the gippasim strain model puts those differences at the stop, and that 16b's e2e test turns them into `stop_codon` ALT sites or `cterm_extension` evidence.
- (#331 review nits) (a) `Search(marked)` and `load_stop_cells`: bounds-check cells read from disk, because an out-of-frame index lands on the stop gap between frames; (b) test the `save_index` missing-decoy-cells error; (c) the `to_interval` round-trip check in `frame_cells` is unreachable; (d) `match_frameshift`'s reverse re-mapping could call `mirror_cells`, and `decoy_marks` repeats `save_index`'s key check.
- (MVP-04 gate 16) The K-12 run took 14.8 CPU-hours against the expected 3–8; profile it.
- (build_envs) `scripts/slurm/build_envs.sh` removes the old env before the new one solves (job 1888853 lost the cascadia env and its locks). Build into a temporary prefix and swap, and set `CONDA_OVERRIDE_CUDA` for CUDA pytorch envs on CPU nodes.
- (R14 memory) Three sbatch scripts ask 16G at 4 CPUs, and their tests pin it; at MaxMemPerCPU = 4000 MB, Slurm keeps them at 4 CPUs only by luck of rounding. Ask `--mem` ≤ 4000M × CPUs.
- (FarmShare) Khoa's `~/.bashrc` exports `RUN`, and sbatch passes it into every job. Scripts must never read `$RUN` (denovo_pilot uses `A3_RUN`); add a guard test for other scripts.
- (A3) Record the `--reference-peptides` database used for A3 (Cascadia gives de novo sequences only).
- (M11-01 check_dist nits) The remaining reviewer nits on `scripts/release/check_dist.py` (see #330's review).
- (#332 snapshot cache) Parallel byte-identity across xdist workers, and when the cache is invalidated.
- (Cascadia CPU) Batch cost varies about 6× across a run; the equal-block array (STEP=predict) is kept only for probes. The work queue (STEP=work) is the supported route.
- (env locks) `workflow/envs/locks` has locks for 6 of 15 envs; decide whether every env gets one.
- (M7-16 statistician, forwarded; `local/stat/m7-16-review.md`, Slurm job 1889319) on the real Nissle genome the exact expected regular chance mass gives a shuffle decoy/target ratio of 0.98–0.99 under the i.i.d. law but 0.93–0.96 under the tryptic law at L = 10–12, and 0.84–0.88 at L = 16, m = 2 (reverse is 1.03–1.07, conservative). The shuffle under-counts tryptic chance hits, which may make de novo placement anti-conservative on a real genome. gippasim does not show it (regular D/T 0.97–1.08). This is multiplicity-weighted mass, not distinct words, at one genome and three seeds: check it properly in the M3-21 / TE8 work before acting on it.
- 16b-3 (M7-16 rule 12): a read-through sORF merges into a `sorf`-class novel ORF (CLASS_PRIORITY), so it never gets a `cterm_extension` row or `readthrough_candidate`. Decide whether rule 12 should keep a read-through flag on the merged sORF. Spec, then test-author → implementer.

### 2026-10-25 (lead)
- (M7-16b-3, rule 14; deferred from #335 by a lead decision) `SimConfig.n_readthrough` with `truth/readthrough.tsv`, and its sim100k e2e test. The feature loop started on 2026-10-25; drop this line when it merges.
- (R14, 2026-10-25 01:52) `scripts/slurm/nightly_launch.sh` submits its 65 CPUs without looking at the other `gippa-*` jobs; that night the nightly, Cascadia (112) and a 16b test job (8) came to 170 CPUs and the lead cancelled nightly tasks 0 and 2, the benchmark resubmit 1891690 and the summary 1891598. Sum the running and pending `gippa-*` CPUs first and defer (or shrink) when the total would pass 128.
- (Pilot commit mismatch, 2026-10-25) `pilot_pass1`-style watchers ran A2 at 0667e366 and A1 at 6f969ee7, because `gippa_run.sbatch` reads the code live from the main checkout and #335 merged between the arms. Watchers should run every arm from a worktree pinned at one commit (as `pilot_a2_rerun.sh` does) instead of the main checkout.
- (Merge queue, 2026-10-25) `queue.T.sh` fast-forwards the main checkout on every merge, so a live `gippa_run` can pick up new code mid-run (harmless for #334, which touched no `src/` or workflow modules). Have it skip the fast-forward while a `gippa-run-*` job reads the main checkout, or run pipeline jobs from pinned worktrees only.
- (#309 follow-up) `scripts/dev/progress_bars.py` misses running M10 branches (it showed 0 running for M10 while branches were open).
- (#335 review) In the e2e test's `consumed()`, add the `+n` shift only on circular contigs.
- (#335 review) `docs/methods/variants.md`: say that a strain stop-loss under `search_db = other` gives ordinary pairs, and that paired read-through raises the pooled class rate.
- (#335 review) A `makedb` patch that ends on the old stop leaves a trailing `X`; handle or document it.
- (#335 review) Add a step-level test of the `loci` circular map.
- (#335 review) The digest pin skips missing files; make a missing file fail, or say why it is skipped.
- (#335 review) Report the per-arm count of REF `*` sites and `cterm_extension` calls as a gate diagnostic (statistician: about 0–0.3 per replicate expected in the G1 arm from natural stop-loss, 0 in the same-strain and novel arms).
- (#335 review) The observer's `drop_duplicates` runs across decoy types; filter to one decoy type first.
- (#334 review nit 5) The offline guard (`scripts/pilot/cascadia_offline.py`) does not block DNS lookups (`socket.getaddrinfo`, `_socket`); block them too.
- (#334 review) Claim Cascadia work-queue chunks by grid cell rather than by chunk name.
- (#334 review) If `local/handoff/run/cascadia_gens.sh` is ever reused, its missing count must become `miss=len(S.pending_chunks(w,20))`.
- #346 review S4/S5 (M7-16c): pin the gippasim selenoprotein exclusion rules against the surviving mutants (C11, `test_gippasim_selenoproteins.py`, incl. M33 the N-extension RBS window); make C12's vacuous acceptance at 100 kb visible (warn or non-strict xfail when `n_acc == 0`). Test-author, small PR.

### 2026-10-25 afternoon (lead)
- **Off-FarmShare install (release blocker, Khoa question):** a `pip install gippa` outside FarmShare cannot pass preflight, because the tool envs and images come from `scripts/slurm/build_envs.sbatch` in a checkout. Default for 1.0: document FarmShare from a checkout as the supported path. Proposal: a `gippa setup --tools` bootstrap that builds the micromamba envs from the bundled `workflow/envs/*.yml` and pulls the pinned images from the bundled `containers.config`. **Closed by M11-14 (spec #367; 14a–14u2, 14d2; docs 14i): `gippa setup --images --profile docker|apptainer` pulls the pinned images and builds UniOP locally, and preflight checks them without a checkout ([Install](get-started/install.md)). Left: the first GHCR publish (14j, Khoa) and AC 16 (wheel-only Slurm job).**
- **V03 locus order mixes scales** (`viz/locus.py::_order`): the pool is sorted by `error`, which is a region q, an island FWER p or an operon-flag q depending on the kind (decisions 3, 5, 8). Sort within kind, or rank by kind first; needs a spec decision (statistician), then test-author → implementer. Spec: [docs/specs/V03-order-within-kind.md](specs/V03-order-within-kind.md) (`stats`).
- **Quickstart demo peptide table:** ship a small peptide table (< 100 KB) under `docs/tutorials/data/` so the quickstart can run `gippa run` with `--peptides`. Prefer gippasim-generated peptides (no sourcing question) over in-silico tryptic K-12 peptides.
- **14 `slow` calibration tests that are not `calibration_full`** run in no scheduled job (M11-03b Q2). Schedule them in the weekly job or mark them `calibration_full`.
- **`test_m10_09_cascadia_contract::test_c4_scan[sample]`** reads `local/pilot/A3_cascadia` outside the repo and fails with "scan not strictly increasing" on the current local data. Skip it when the data is not the pinned one, or pin the data by digest.
- **`test_sim_function::test_example7_under_one_second_for_4000_genes`** times out (about 1.12 s) on busy Slurm nodes. Make it a benchmark or loosen it only with the statistician's and lead's sign-off (R3 does not cover timing, but keep it honest).

### M11-14 spec review nits (#367)
- #370 review N4 (14u2): preflight hashes the UniOP `.sif` (≈ 1 GB) on every run; cache the hash by `(size, mtime)` in `local_build`, or record that the cost is accepted.
- #370 review N5 (14u2): `scripts/build_envs.sh` builds the FarmShare prebuilt `gippa-uniop` with an unchecked `git clone`, while `setup --uniop` uses the sha256-checked archive; move `build_envs.sh` to the same verified archive, or record why the two paths may differ.
- Code review N3: the published images redistribute GPL-3 tools (OSTIR, PromoTech, Promoter Calculator, BacTermFinder) and GPL-2 R (fgsea image). Add one line to the spec and NOTICE on how corresponding source is offered (the `org.opencontainers.image.source` label, the pinned commit and archive URL, NOTICE links), or cite M11-06 AC 4 and extend it to fgsea/R. Docs role, with 14c or 14i.

### 2026-10-25 evening (lead)
- #350 review backlog (M7-16b-3): spec lines on `n_readthrough = 10` (around the 16b e2e items) and `tests/unit/test_sim_readthrough.py` still cite 10, not 40; de novo errors give two `multi_nt` candidate ALT rows at planted stops (assert every REF `*` site is a planted stop?); add floors (`assert called` in the patch round trip; planted stops ⊆ index terminal sites); give the consuming-rows test its own floor constant; fix the `expressed_proteins` docstring on how the interval grows; pin `pick_readthrough`'s feature order; drop stale `# type: ignore` on `n_readthrough`.
- #349 review backlog (M11-03b): record part names are not tied to their task; squeue states and two launchers racing; a space in the log path breaks the JUnit option; a small temp-dir leak; `pyproject.toml`, `uv.lock` and `tests/conftest.py` are missing from `CALIB_PATHS` and the sbatch `dirty` check; the launcher does not pass `--sacct-json`, so OOM/TIMEOUT states only show when the aggregator is run by hand.
- Nightly launcher has no R14 budget check (Khoa question in `khoa_q.md`; the full-tier launcher avoids the nightly window instead).

### 2026-10-27 (lead)
- **#445 review nits (calibration plan; S1–S3 fixed in the PR):**
  - N1: `probe_log_lower_bounds` across several probe logs keeps the last file's bound, not the largest (`calibration_lanes.py` main, `lower.update`); under `--max` a `:lower-bound` row is never replaced by a smaller real measurement; interrupted JUnit lanes give no lower bounds.
  - N2: a lower bound assumes the test ran until the log's last stamp; a silently dead worker overstates it. Add a docstring caveat.
  - N3: mutant M8 (lower bounds override measured values) survives; add a two-file CLI test (finished in log A, unfinished in log B).
  - N4: mutant M9 (`>` to `>=` in `--max`) survives; a tie relabels the source. Harmless.
  - N5: `test_durations_rows_are_the_max_of_old_and_new` checks one row; compare against `origin/main`'s `calibration_durations.tsv` and assert no row went down.
  - N6: the per-model determinants rows (`1901050:measured`, 23579 / 78976 s) are hand-entered and cover block 1 only, so the ⅔ test passes determinants on block 1 only (the PR's known gap).
- **Flaky test watch:** `tests/calibration` `test_interrupt_rescans_for_pytests_started_after_usr1` failed three times on CI under load (group start vs USR1 rescan); #448 tolerates SIGINT during the late pytest's start-up. If it recurs, reopen.
- **V18 macOS SVG determinism (possible R7 bug):** #446 skips the FarmShare-only tests on macOS, but `V18_stat_diagnostics.svg` was not byte-identical across two renders on macOS (PYTHONHASHSEED 3 vs 4; Linux passes). Render V18 on macOS several times and diff.
- **Setup-test network guard:** setup tests made real apptainer pulls (about 6 GB each) and filled /tmp (`test_setup_run_never_starts_a_real_quota`, fixed in #441). Add a shared fake `apptainer` that refuses `pull` for all setup tests (drafting in branch `test/setup-network-guard`).
- **Calibration ⅔ check ignores the merge phase:** `test_lpt_estimate_within_two_thirds_of_time` weighs only the lanes, but the merge runs inside the same job time (dry run 1900347: novel_orfs merge killed at 0.95 h after 6 shards × 4.7 h). Add the merge's duration to the estimate.
- **NFS TMPDIR partial-dir tests missing:** no test covers the `.partial` staging paths (`setupdb/run.py`, `uniop_build.py`) when TMPDIR is on NFS.
- **Preflight GBFF warning:** warn when the bgc module is on and no `<stem>.gbff[.gz]` sits next to the precomputed Bakta JSON (antiSMASH is otherwise skipped with `layer_missing`); name the expected path.
- **TerminatorNet sharding:** TerminatorNet runs super-linear on large high-GC genomes; a spec is in progress in branch `spec/terminatornet-shard`.

### Deferred review items, 2026-10-21 to 2026-10-23 (lead's parked list)
Merged from the lead's local list on 2026-10-25. Skipped as already done or already in this file: the bash-4 check (#270), Promoter Calculator batching (#282), the M9-17 Nextflow wiring (#283), the M7-18 X/U crash (#284), and M7-17 Q6, Q8 and the G4-k3 reuse (under "M7-17 follow-ups (PR #291)").
- (M7-11b decision 10) gippasim/makedb ID mismatches left: annotated peptides in the first 30 codons do not also list the gene's nterm entry; 27 non-ATG-start peptides beginning with M are absent from makedb's literal translation; 18 rows name entries shorter than min_orf_aa=20.
- M8-15 attack 2: keep `member_shift` −1.5 as the separation arm, add an overlap arm at −0.75 (unweighted 0.95 at 0.05, R=150; confirmed fail at R≥100); run both with the default `DetectabilityConfig()` on covariates recomputed by `covariate_table`, not `n_obs_true` (statistician #251).
- ADR (M8-04): the fitted ORA treats q̂ as known, so it is conservative and loses power when membership is predictable from covariates; proposal: score test of the category indicator in the M4-13 GLM (residual-projection variance) (statistician #251).
- gippasim `_confound`: K/R moves change `tm_segments` in 66 % of members (+1.56) and pI by ≤0.22; document or constrain the moves to keep `tm_segments` (statistician #251).
- M8-15: the oracle arm through the ORA engine is slightly conservative at α 0.10 (0.025–0.067); check M8-04's TU recalibration/TU unit on this world (statistician #251).
- M8-14 code-review nits (comment 6367524340 on #251).
- M6-12 `precision_cv`: add a reference test rebuilding the folds (`permutation(n_clusters) % k`) and comparing pooled held-out precision exactly (row-fold and in-sample mutants survive) (statistician #263).
- M6-12 `fit_distance_model`: after 30 failed halvings the 2⁻³⁰ step is judged on |h·δ|; judge convergence on |δ|/gradient or fail with `line_search` (statistician #263, unreachable in practice).
- M6-12 attacks: strengthen the directon-clustered null so the iid CI under-covers; add the mixed-length/both-strand/circular coincidence-null case (statistician #263).
- Nightly task 13 (M7 variants gate) is red while ADR-0049 is open: add a merge-side known-failures list naming ADR-0049 so the nightly summary separates expected from new failures (lead, #255).
- M6-12 code-review nits (#263 comment 6367759779): public alias for `_window` (+ docstring "same-strand"), contig checks in `window_of`/`_containing`/`coincidence_rate`, step-halving floor comment, `threshold_boot_q` comments, document the coincidence stream order, nested-gene windows, input validation (`wilson_ci` k>n, `distance_log_loss` s=0, `precision_cv` k=1, empty `choose_cutoff`, all-equal d reported as separation).
- pyproject `pythonpath = ["tests/sim"]` shadows an explicit PYTHONPATH gippasim; cross-branch measuring needs `-o pythonpath=…`; document in AGENTS.md or test_branch.sbatch (M7-13b implementer).
- M6-12 B code-review nits (#264 comment 6369206328): atomic outputs; guarded `sys.path` (or `python -m`); `main` catches only ValueError; hash the truth TSVs and refuse unknown labels; `y` from `matched_source_id`; `--R` flag and R ≥ 1 check; separate streams so a TerminatorNet re-record doesn't move the UniOP `precision_cv`; driver size cuts (move `_hits`/`_called`/`truth_recall` to gippa.stats.cutoffs, `coordinate_id_map` to uniop_adapter, `directon_ids` to operons, `open_text` for gz, `to_csv` writers, split `run()`).
- M8-15 follow-ups: centring arm (ADR-B) and annotated-universe arm (spec, not yet written).
- M6-12 B (#264 stat comment 6369349361): ADR-0022 quotes `precision_lower_bound_window` as the 97.5 % one-sided bound; the point bound exceeds true precision in ~9.5 % of worlds, so it is never quoted as a bound.
- ADR (M8-09): LightGBM at the §5.2 defaults overfits at ~800 model rows (CV R² −0.002 at true 0.2, 1,500 genes; OLS 0.183); options: inner-fold early stopping, 100 trees, rate 0.01 + min child 50 + 7 leaves (statistician M8-15 P2, local/stat/m815-review/).
- M8-04 projection-variance ADR: also covers the oracle arm's shortfall at −1.5 from the TU recalibration (0.02 vs 0.05 without it, α 0.10).
- M8-15 recovery arm: route the planted features through the M8-07 columns (`promoter_score`, `tir_ostir`, `cai`) as determinants.md says, or amend the note.
- Re-file the item "M8-15 should adopt the statistician's GO-family FDR arm" (M8 section above) under hardening; §8.4 is level-only.
- GSEA null pooled band too narrow for 100 dependent sets per replicate (0.0396 vs 0.0438 at α 0.05); recorded only.
- M8-15 code-review nits (#266 comment 6370059803): forking after LightGBM's OpenMP can hang (use spawn or prefetch before any model fit); duplicated helpers; per-collection records take their collection list from replicate 0 only; document that the fast tiers reuse the full tier's first worlds; the baseline was timed at 16 CPUs, not the spec's 4.
- M8 gate: run the M8-15 full tier (calibration_full: ora/gsea/determinants) on Slurm and record it in PROGRESS before marking gate M8 met.
- Calibration teeth (house `teeth_log`) rerun both full blocks in every shard (~250 s per tag in M8-15); consider skipping teeth under --shard and running them once in --merge (#266 test-author).
- M8-16: no nf-test for SENSITIVITY_RERUNS (sensitivity.nf) after it left CONTEXT; add one. Also `docs/methods/sensitivity.md` "Enriched categories"/"Determinants" sections (spec acceptance) and a `requires_r` test of the real fgsea engine in the family (#265 implementer).
- Calibration records carry no commit field; add one so gate evidence can prove the commit (M7-13c ledger gap 6). Remove the `@GATE` switch from test_m7_gate.py::ledger once #255 and M7-13b merge (gap 4).
- test_m7_gate.py: re-derive merged records with C.merge_shards to catch a stale merge (#268 nit 4); docstring pointer that strand/circular-origin claims live in O3 (nit 7).
- M7 nightly: shard the novel_orfs block 2 (serial G4-k3 block 2 ≈ 3.2 h vs --time 04:00; M7-13b decision 12). #267 follow-up `test_attack_oracle_picked_class_fdp_dense` "no calls at 0.01, intergenic": resolved in PR #449 (stale oracle key in the test).
- gippasim scored_null mode: M7-11b makedb-ID mapping is off under the default wrong_source (correct novel-ORF peptides name planted IDs; protein nulls upstream of the start name the locus tag). Kept #247 behaviour so K-12 scored-null outputs stay unchanged; decide whether the guarantee should hold in scored mode (#267 merge, lead default: keep).
- ADR-0049 #269 nits: `site_flags` never applies the left-block rule (document or pass through); `op.alt` never compared with the site alt_aa; mention `identifiable` in the M7-05 spec; the `lo >= 0` guard in left_block_read is provably dead (remove or keep as a guard).
- M9-17/M9-19 (V14 implementer 7410ec5b): `gippa.steps.plots.plots_step` swaps `inputs["plots_support"]` for its own `out/support` whenever `write_plot_support` writes (replicons or variants present), so `gsea_running` published by the function step may be invisible to V14 in a real run. Check when M9-17 publishing lands; add a step test with replicons + function support.
- M7-13b #271 review nits: shard merge is xfail-unaware (a sharded task 11 is red nightly on the six ADR-0050 cells; note beside ADR-0049 "known consequence"); injection helpers lack unit tests and `_reverse_stretch` leaves palindromic stretches unchanged (~5% at length 3), blurring ADR-0050's exact vs reversed split; `fdp_of` negative-count and bool branches untested.
- M9-26 gap 4: V2 rug is unbinned; a genome with >20,000 undetected genes breaks §10.1. Default: bin by status share as V1.
- M9-20 nits: N2 "spread over 25 CV models" comes from plotting cfg not the run (use recorded fit count); N5 "jitter by rank" comment wording; M8 elastic-net out-of-fold predictions so V15 (d) can plot the primary model.
- M9-19 nits: N9 `_draw_set` xlim from x.max(); N10 split `_draw_b` label layout; add "es" to checklist UNITLESS ("running ES (enrichment score)" passes only via parentheses).
- M9-21 nit: V18 captions print p 0.996 as "1" (.2g); use .3g above 0.1.
- M9-17 nits: N4 gsea_table rank labels via DataFrame.attrs (return explicitly); N6 pred_contrib computed twice (derive shap_blocks from shap_features); N7 redundant method!=descriptive filter; drop PD grid points <= 0 for LOG_FEATURES.
- ADR-0050 option 1 (accepted for later): picked ORF-level target-vs-mirrored-decoy competition to remove the reversed-stretch k-chances excess at q 0.01 (FDP 0.0135–0.0140 vs ~0.011–0.0125).
- M9-26 nits: byte_identity.py:91 puts repo src ahead of PYTHONPATH; V1 binning decided on total gene count, not per track (§10.1); needs_nextflow requires nf-test unnecessarily for the resume test.
- M9-19 V14 after the review round: the pinned upper-left legend covers the cog J point and its leader in the main case (needs more y headroom or a 2-column legend; check at M9-25); q = 0 and a q > alpha row can share the cap height (only the ^ marker separates them); the twin alpha line is fixed at 0.05 (could read the alpha from the title).
- M7-18 follow-ups: B/Z/J as ambiguous matches (Module 2 first); preflight count of B/J/O/U/X/Z in peptide tables; Sec/Pyl readthrough calls (U over TGA, O over TAG) belong to the M7-16 ADR (§7.5.1 vs M7-15 decision 9).
- M9-20 V15 after the review round: in the big-genome case 1,025 flagged points cover the hexbin and the "best" legend sits over bins (smaller or translucent flagged markers in binned mode; check at M9-25); the twin can't show the null-p and univariate footnotes (not in its data TSV); the colourbar test only searches figure.axes.
- M9-17 follow-ups: write_datapackage copies only the schemas resources name, so null_regions_grid.json is missing from OUT/schemas although datapackage.json mentions it; only the local publish drops the stale .tsv/.tsv.gz sibling (the Nextflow side belongs to m9/support-nf-wiring).
- chore: farmshare.config queueSize 8 × 16 CPUs plus the head job can reach 130 CPUs (over R14's 128); cap it (M10 L1).
- M9-21 V18: the constant-null bar (width 1/30) is a near-invisible hairline on wide x ranges; give it a minimum drawn width (check at M9-25).
- M7-13b #271 follow-ups: gippasim names origin-crossing entries b > n (the exemption checks a > b; fails closed); G4-k3 reach pooled over classes (split per class); review nits 4–7, incl. merge_shards reporting MERGE_FAIL for the six ADR-0050 cells once task 11 is sharded.
- chore: needs_diann marker and probe_env don't know the verified DIA-NN copy (M10-02 Q8).
- M10-01: test the 1/K0 range (needs TDF mobility calibration) and the diaPASEF-variant notes (Q3).
- M7-15 follow-ups: static guard for the step→config-section map (registry.py:184); GIPPA code version in the GIPPA_STEP task hash; preflight keep-unchanged hunk + all-cached assertion; record cache 'deep' hashing cost.
- M10-02: test the free-space check (rule 4) and the mv -n race between two install jobs (need fake df/mv hooks).
- M9-17nf: unit test for plots.copy_staged_support; determinants.nf emits no events into the PLOTS events pair.

#### M7-18 review nits (#284, 2026-10-23)
- Test the folded-peptide dedupe (`AKXREFGHILWD` + `AKXREFGHLLWD` → 1 site, 1 peptide, both originals in examples), pileup.py:245.
- First skipped letter wins via setdefault (pileup.py:245): order-dependent by_residue/by_rule for overlapping hits; use min() or pin first-wins.
- Pin rare/ambiguous for mixed ops (`KU` rare, `UX` ambiguous), pileup.py:150.
- Assert the event message wording (decision 4) in test_step_emits_one_info_event_when_a_site_is_skipped.
- Acceptance criteria: say output byte equality rests on the goldens, not the property.

#### M2-14b review nits (#282, 2026-10-23)
- test_pc_batch_resume_stub.py docstring (and spec lines 21, 39): the guarantee is "unchanged inputs fully cached, batch membership stable"; a real region change reruns every batch (OPERON_PAIRS new work paths; TU add/remove shifts boundaries).
- promoter_calculator_batch.nf.test example 6: put the failing region in the middle ([tu0001, tu0004 FAIL, tu0003]) and assert 2 `Sequence length:` lines (tu0003 never ran).
- PC_BATCH copied in annotate.nf, workflow/tests/real/promoter_calculator.nf and test_pc_batching_real.py DRIVER: make the static test assert they match.
- Record the K-12 PROMOTER_CALCULATOR wall time in PROGRESS.md after the first real run (MVP-04 / M6-12C); expected ~4–5 min per 100-region batch.

#### Fast-suite test-side speedups (perf agent, 2026-10-23; after #286)
1. gippasim `function_sim.draw` ~1,000 CPU-s/run (test_sim_function_calibration examples 1-3, attacks): cache worlds per replicate where arms share a stream (check streams align), or vectorise layout/_confound (gippasim goldens as identity check).
2. test_gene_fdr_one_replicate_through_real_mapping (265 s after #286): reads only gene_table output; restrict modules / lower determinants.permutations via --set → ~20-40 s.
3. test_mapping_steps_review: module-scoped shared sim50 run for read-only items (11, d13, d14), ~90 CPU-s.
4. Local-run integration tests (m4_real_path, gene_table_local, mapping_local, novel steps, id_mixups): session-scoped shared sim100k run dir, ~150 CPU-s.
5. test_nf_review: share one nextflow config render between the two farmshare tests, ~45 s.

#### Flaky property test (Slurm 1879780, 2026-10-23)
- tests/property/test_m5_fixture_props.py::test_oric_is_one_chr1_feature_near_its_fraction: Hypothesis seed=759 puts the oriC midpoint 0.0236 from the target fraction (tolerance 0.02). Simulator fixture issue (gippasim), not an M7-17 change: tighten the fixture's placement or derive the tolerance from the gene spacing; pin the failing example with @example.

#### #286 review nits (perf, 2026-10-23)
- N1: models.py:207-219 relies on sklearn private _path_residuals/_pre_fit; pin sklearn in pyproject or assert Gram contiguity/dtype; comment the verified version (1.9.1).
- N2: _rho uses np.corrcoef(vstack) where scipy uses corrcoef(column_stack, rowvar=False): use scipy's exact call (BLAS layout on macOS).
- N3: _ranks/_auroc assume clean inputs: convert y to bool, assert no NaN.
- N4: bootstrap.py:56 stat_fn type for columns= path (Callable[..., float] or overloads).
- N5: _fit_predict threads param always 1; null.column() used once.
- N6: ADR-0053: bit identity with main shown empirically only; old 16-thread LightGBM CV was not guaranteed thread-count-invariant on large inputs.
- N7: benchmark null budget 45 s > 33.74 s before: set ~20 s.
- N8: null_p holds all B permuted y copies: docstring line.

#### M9-16b review (#287, 2026-10-23), non-blocking
- Test temp-file cleanup: monkeypatch os.replace to raise OSError; old page bytes kept, no .report.html.tmp left (plots.py:190-191).
- `gippa report` size warning untested now that it goes through rerender_report (plots.py:192): caplog test for one report.size_exceeded WARNING with data {"bytes": size}.
- Behaviour 2/7 negatives: spy tests that exit 3 and exit 5 never call rerender_report; Nextflow test globs work/*/*/report/report.html for the in-graph "not written" sentence.
- cli.py:380 reads the manifest outside the try (a read failure is exit 5, not the tail's exit 4); comment the mid-function `plots as P` import; gippa report builds RunView twice; docstring: an oversized page logs report.size_exceeded twice (graph + tail).
- backlog: _orf_span no-start fallback can give start < 1 near the replicon start (+ strand); handle wrap (pending)

#### M8-12b statistician (#289), non-blocking
- PROVENANCE prints Σf |Σ−1| as 0.000e+00 (float underflow; the true value is about 2^-1228): print log2 of the bound.
- Optional: the exact integer fixed-point urn (scratchpad fixrec.py, 640-bit, about 15 s) as a CI cross-check on all 300 tails.
- Spec §1.3 nit: SET_026/027 keep their structure, but their members changed.

#### M8-12b code review (#289), deferred
- S1: the §3.4 checks never fail in a test (run(full=True) with a perturbed R1; _check; CLI exits for unknown category or gene); 9 mutants survived.
- S2: pin 2-3 R1-regime tail strings (or the dtype), since a float64 reference survived.
- S4: golden_wallenius.sbatch should require GOLDEN_WALLENIUS_REASON (with :?) plus a Supersedes line, like golden_r.
- Nits: the precision loop overshoots P (one-shot estimate); duplicate R1 work in exact_tails; _read_tsv has no header check and no duplicate (category, gene) check; _numpy_odds rereads the files; dead k-outside-support exit; untested r2_pmf d_num == 0 and k == lo branches; make_golden_inputs held.add collision untested; golden_r -e check is per file; golden-r.yml has no ORA guard or Supersedes line; vacuity printed as a dict repr; split run() (about 95 lines).

#### M10-01 code review (#290), deferred
- Tests for: --proto-redir on GET/HEAD; exactly one analysis.tdf (2 .d folders in one archive); files outside <run>.d dropped; the unpacked skip needs tdf_bin; an oversized partial discarded; a genotype annotation modification ("wild type ΔfliC"); WT token left boundary (SWT_1.d.zip); acquisition_mode source_folder sidecar.
- Selector: WT_ROS also matches the WT token. Consider rejecting treatment tokens (ROS) explicitly rather than relying on MLA order.
- Nits: --resolve can't exit 0 on the real record (http:// links; https_route only upgrades ftp://); the property test's _expected copies classify; a missing tdf drops the dia_window_rows/pasef_rows keys; IsolationWidth/2 on an int column; immutable=1; an unknown checksum algorithm gives a traceback, not exit 2; run_label .d is case-sensitive; the override URL should use download.iprox.cn.
- CLI: --annotation *.json.gz silently runs Bakta live (cli._bakta_json checks only .json). Accept .json.gz (decompress, or have the step read gz); test-author → implementer, small PR (M10-04 Q3).

#### M11 wave B reviews (2026-10-25), deferred
- **#338 M11-06:**
  - S1 (Khoa): a digest PR opened with GITHUB_TOKEN triggers no CI and needs "Allow Actions to create PRs". Use an App/PAT, or document a manual close/reopen.
  - S2: a dispatch from a non-main ref opens a PR against main.
  - S3: the `git-<sha8>` tag comes from the env comment, not the Dockerfile ARG.
  - S4: the process→tool map is duplicated (containers.yml and update_digests.py).
  - S5: pinning only PROMOTECH_NOVEL drops promotech from the conda-only comment.
  - S6 (Khoa): the licences of conda packages bundled in public images (PyQt5 GPL-3.0, ViennaRNA).
  - Nits: a double digest is accepted; writes happen outside the try; a re-run after the branch was pushed fails; `set -f` for the tools input; two back-to-back pushes give two conflicting PRs.
- **#339 M11-07:**
  - S3: tests for the hidden-option filter, the empty `multiple` default, and `_status` multi-line text.
  - S4: the hand-written Methods/Decisions nav conflicts between parallel agents. Generate the nav, or list only the index pages.
  - Nits: exclude adr/TEMPLATE.md; flags.md shows a raw regex; the `_unset` sentinel is detected by type name; CRLF passes `--check`; no site_url/repo_url; envvar not shown.
- **#341 M11-12:**
  - B1 follow-up: the no-changelog label takes effect only on the next push. Read live labels with `gh pr view --json labels`, or add labeled/unlabeled triggers.
  - S2: a rename inside changelog.d/ passes as a new fragment; the `C` status branch is dead.
  - S3 tests: `..` vs `...` with main moved on; renamed or edited fragments; pattern anchors; the `": "` separator and `refactor!:`; exact label match; a merge-only range.
  - Nits: a bad BASE gives a traceback; non-ASCII paths are quoted; a module moved out of src/ needs no fragment; `open` exists only on macOS; README line length; the fragment slug wording; link regex vs parentheses.
- **#342 M11-04:**
  - S1 (Khoa): the actor guard can't tell agents from Khoa (the FarmShare gh login is khoa-yelo). Use a separate account or a token without actions:write; also check `github.triggering_actor`.
  - S2: readiness fatal for pypi is untested.
  - S3: skip-existing, packages-dir and `--draft` are untested.
  - S4: a test against `${{ inputs.* }}` in run scripts; `$GITHUB_SHA` at :192.
  - S5: notes.py tests for the preamble, "No significant changes." and rc with `--draft`.
  - Nits: `::error::` echoes the unvalidated version; stderr vs stdout; `persist-credentials: false`; notes.py fences and `# ` headings; a multi-digit rc patch; example 5 needs a built section.
- **#337 M11-03:**
  - Item 11: no job writes `calibration_full` / `calibration_attacks` rows to nightly.tsv, so readiness R3/R4 stay pending. Add those nightly rows before M11-13.
  - A `--check-runs` file without `head_sha`.
  - An evidence path recorded as a basename only.
- **M11-04 / M11-02:** an M11 milestone fragment will break test_m11_02's "exactly 12 milestone fragments"; update the count when M11's fragment lands.
- **Earlier parked:** datapackage.json lacks `genes` and others; config.schema.json fields have no descriptions (gen_reference warns).
- **Rebrand (Khoa, 2026-10-25):** after development, restructure and rename to GIPPA (repo, package, code, simplification). Do it before the first real publish (M11-13). Sequence: an ADR (overrides the spec name), one scripted rename PR verified byte-identical apart from names (PR-size exception), then simplification PRs gated by goldens and calibration. `gippa` was free on PyPI, TestPyPI, bioconda, conda-forge and GitHub on 2026-10-25.
- **V03 order review nits (#360, one round; statistician BLOCK fixed in the PR, code review APPROVE):**
  - `V03_locus.tsv`'s `kind` column cannot tell the two operon-flag families apart; append an `estimate` column later.
  - The spec's Edge cases table lacks the row on re-rendering at a higher `alpha.context` (the lead decisions call for it).
  - Q4: show the V03 and V13 locus plots in the report in index-rank order.
  - `FAMILY_ORDER.index` raises `ValueError` for an `operon_flag` Locus built without a family (only `_flag_loci` builds them today).
  - `_order` looks up the family position per locus and `_index` recomputes the tier per row (negligible).
  - The user guide could say "BH q-value within each region type".
- **#359 review (demo peptides; B1 and S1, S3, S4 fixed in the PR):**
  - Bug (src, own PR): `run_manifest.json` `detection_tier` copies preflight's provisional tier (`src/gippa/manifest.py:374`), so it says `decoy_picked` while genes.tsv and detection_summary.tsv say `mayu` (demo run 1894327: unscored decoys, rule 4(c)). The `decoy_unscored` warning is also missing from the manifest warnings.
  - S2: no CI test runs `scripts/docs/make_demo_peptides.py --check`; a later gippasim change could silently break the committed table. Test-author.
  - Option: have the generator emit scored decoys so the demo shows the `decoy_picked` tier.
  - Nits: the `cp ../tests/...` step assumes the repo root; download URLs point at unpinned `main`; the gz checksum can vary with the zlib build (`--check` compares uncompressed bytes); README keeps `--profile docker` while the quickstart says FarmShare/apptainer; `MAX_BYTES = 100_000` vs "100 KB".
  - Demo run warning: one work dir had an unreadable `versions.yml` (ScannerError), skipped by `gippa.io.manifest_sources`.
- **M10-04 (#296) pass-1 findings (validation, 2026-10-25):**
  - antiSMASH did not run in either arm: the precomputed Bakta JSON has no GBFF next to it, so the step is skipped while the bgc module is on (A9 passes because `modules_off` is empty). Write or locate the GBFF for pass 2, or make the skip a `layers_missing` warning in the report.
  - A2 lacks `rbs_novel`; its novel-ORF promoter and RBS layers are `not_run`.
  - `pilot_compare.sbatch` does not record the commit it ran at; record `git -C $GIPPA_REPO rev-parse HEAD` in compare.json.
  - Superseded outdirs renamed after the run keep `out_a2` in their manifests' `command_line`.
- **Review follow-ups (lead, 2026-10-25):** preflight `_has_ranking` reports `decoy_picked` when the decoy ranking columns are all NA; gippasim `scored_null` yields 0 protein-entry nulls; M11-08 spec line 85 is stale; the 14 slow calibration tests need a Slurm home.
- **#364 M10-10a review (APPROVE; S1/S2 tests added in the PR):**
  - S3: `pilot_compare.sbatch` asks `--mem=8G` (8192M on 2 CPUs allocates 3 CPUs at MaxMemPerCPU 4000); the M10-10 Jobs table says 8000M. M10-04's `test_m10_04_compare_job.py:201` pins 8G; change both together.
  - `cascadia_convert.py` docstring says `src/` is put on `sys.path` (it is, only through the `diann_convert` import); `a3.convert.json` `command` is rebuilt, not the real argv; the 2-label `compare.json` has three extra keys; `A3_OUT_DIR` also redirects the casanovo branch; `$A3` has no outside-the-checkout guard; `--arm A3` A5 requires `reference_peptides` only when passed; a stale `a3.FAILED` stays when the gate or a precondition exits; seven untested boundary cases (review comment).
- **#363 (manifest detection_tier):** the registry path `gippa _step manifest` (`manifest_step`) has no `genes` input, so it still records preflight's tier; `gippa run` (`run_tail`) is fixed. Add an optional input if that path is ever used for real runs.
- **#363 review nits:** `read_detection_tier` does not warn when genes.tsv and detection_summary.tsv disagree; `test_manifest_cli_local.py:223` and `test_manifest_cli_nextflow.py:162` still assert manifest tier == preflight tier (pass only because the sim's tiers agree; point them at genes.tsv); `test_manifest_detection_tier.py` mixes import styles; `manifest_step`'s docstring should say it records the provisional tier.
- **#366 M10-10b statistical review (LIFT; one round, nits deferred):**
  - H4 level depends on bin size: `a2_false_confirm` false `not_supported` 0.00 at 300–1,000 peptides/bin, 0.02 at 2,000, 0.27 at 5,000 (scratch seeds, statistician). Add an H4 / methods-table caveat (valid up to ~2,000 peptides per bin at 1 % A2 false confirmations), a reported (not gated) `a2_false_confirm` n=5,000 attack scenario, and Khoa's "correct `L_b` for A2's FDR" amendment before H4 is reused on a larger dataset.
  - `validation/pilot/report_template.md:226` should say `not_supported` needs two consecutive reportable bins.
  - No unit test pins the strict `<` at `genome.ci_high == L_b_adj.ci_low` (equal reads `consistent`).
  - `tests/calibration/attacks/test_attack_h4_lower_bound.py` docstring still describes the one-bin rule (statistician-only docstring edit).
- **#374/M11-14p review nits (one round):** update_digests follows only line-leading `includeConfig` (inline profile includes are not; the docstring should say so); with docker `-u`, `HOME=/` is not writable in the container and the tools' cache writes are untested (the probe runs `bakta --version` only); a `db_dir` with spaces or `:` breaks the bind argument (apptainer too); `docker_db.nf.test` without `GIPPA_SMOKE_DB_DIR` leaves a temp dir per load; the actions warn about the Node 20 deprecation (as the existing job).

## Genetic code 4 support (post-1.0)

Khoa, 2026-10-27 (Q-K30): GIPPA 1.0 supports genetic code 11 only; preflight refuses any other
`annotation.translation_table` (#413) and warns when a default-11 genome looks like code 4
(M11-14w, `genome.genetic_code`). Code 4 (Mycoplasma, Spiroplasma, Mesoplasma, Ureaplasma:
TGA = Trp) is planned for a later release. Places that assume table 11 today:
- Bakta: pass `--translation-table` (Bakta accepts 4, 11 and 25) from `annotation.translation_table`.
- Six-frame index stop set: `src/gippa/mapping/sixframe.py` (`translate11`, `STOP`, the
  consumable-stop cells); the index format version must change with the table.
- `src/gippa/annotate/codon_usage.py` (raises for any table but 11; codon-bias is called with
  `genetic_code=11`) and `src/gippa/steps/codon.py`.
- Novel candidates: `src/gippa/novel/candidates.py` (`_stops(table)` is already table-aware; check
  its callers and `outputs.py`, which hard-codes 11, see the M7 nits above).
- `src/gippa/mapping/loci.py`: `_PLUS_STOPS` is hard-coded to table 11; `_RECODED_FLAG` (TGA =
  `sec_site`) means something else under code 4.
- `src/gippa/genetable/covariates.py` (reads the table; check every use).
- Preflight: allow 4 in `_translation_table_check` (`src/gippa/preflight/run.py`), and let the
  M11-14w warning suggest `annotation.translation_table: 4`.
- Possible post-annotation corroboration for the warning: the median CDS length, TGA-terminated
  CDSs followed in frame by a same-strand CDS, and W enrichment among unmapped peptides.
- Evidence (M. hyorhinis SK76, GCF_000313635.1, PR #411; docs/validation/datasets/README.md on
  `validation/more-datasets`): at table 11 Bakta called 1,474 features against about 750 in
  RefSeq (genes cut at TGA); the median protein was 90 aa against about 287 elsewhere; all 628
  unmapped organism peptides contain W; 246 `strain_difference` and 354
  `nterm_extension_evidence` flags were artefacts. The M11-14w long-ORF ratio is 0.227 for this
  genome, against 0.82–0.88 for three code-11 genomes.

- **#420 review nits (nightly task-2 shards; N4 noted in the nightly.sbatch header):**
  - N5: tests for shard 2 recording first (shard 1 writes the row; its gate runs in the shard-1 worktree), the "coverage not measured (k of N shard data files)" path, and two concurrent record + merge calls giving exactly one complete merge.
  - N6: after the first sharded night, compare the shards' `--durations`/junit times and record them in PROGRESS; balance by duration or raise `NIGHTLY_CAL_SHARDS` (≤ 4, R14) if a shard passes ~4 h. The nightly now takes 81 CPUs at 01:30 without checking other `gippa-*` jobs (calibration_full overlap).

- **`gippa setup --from <dir>`: offline install from pre-downloaded archives (docs PR "download the
  databases by hand", 2026-10-27).** Today setup adopts a copy already unpacked at its path when it
  verifies (`install` → `verify`), so most by-hand installs work (docs/get-started/manual-databases.md).
  Gaps: (1) KEGG cannot be adopted: unpinned, `verify` trusts only a checksum recorded by setup's own
  download, so a hand-placed `kegg/*.tsv` is `absent`; (2) Bakta and geNOMAD are checked by their
  version file only, after the archive the user verified by hand is gone; (3) `--check` never records,
  so users must know to rerun plain `setup`. Proposal: `--from DIR` takes the archives and files as
  published (`db-light.tar.xz`, `genomad_db_v1.9.tar.gz`, `RF-HOT.zip`, the GO/COG/contaminant files,
  the UniOP tarball), checks each against the `db_sources.yaml` pin before unpacking with the existing
  safe `_unpack`, and records `url` as the official URL plus `message: "from <DIR>"`; for KEGG it records
  the sha256 of the given bytes and the file mtime date as version (ADR-0017 style), flagged as
  user-supplied. No network. Needs a spec item and an ADR for the KEGG provenance rule.

## From-spectra guide gaps (docs/get-started/from-spectra.md, M11-15 item 15)

- **No DDA converter.** `gippa convert --from` reads `diann` (DIA), `cascadia` (DIA de novo) and
  `generic` only; DDA users must write the generic table by hand. Candidates: MaxQuant
  `evidence.txt`, FragPipe/MSFragger `psm.tsv`, Sage `results.sage.tsv` (each with its decoy
  rows, q-values and intensities mapped to the generic columns; contract tests from recorded
  outputs, as for DIA-NN).
- **`--from casanovo` and `--from instanovo` exit 5** (accepted choices, no converter). The pilot's
  documented A3 fallback is Casanovo 5.2.1 (`denovo_pilot.sbatch casanovo`, mzTab), so a Casanovo
  mzTab converter is the first one to add. Record Casanovo's licence in NOTICE.md with it.
- **The Bakta JSON is not published.** `docs/user-guide/outputs.md` lists `annotation/tools/`
  (`bakta/`, …), but `workflow/modules/local/bakta/main.nf` has no publishDir and runs have no
  `annotation/tools/` (e.g. the pilot's pass-3 `out_a2`). `gippa makedb` needs `--annotation`
  (exit 2 without it), so the from-spectra routes tell users to run Bakta by hand. Publish
  `annotation/tools/bakta/<sample>.json` and `.gbff` (then `makedb --annotation OUT/...`), or fix
  outputs.md.
- **Six-frame route never run on real data.** The pilot databases had no six-frame ORFs, so the
  guide has no run time or size for DIA-NN against `makedb --six-frame-orfs`, and the picked
  class FDR has not been checked on DIA-NN's decoy rows (do they carry the `GIPPAORF_…` IDs so
  that `novel_orfs.decoys_not_makedb` stays quiet?). One Nissle pass-3 run would settle both.
- **Not recorded:** the wall time of Cascadia per run (the pilot's CPU-route `.ssl` came from
  work-queue arrays), the output size of a `gippa demo` run, and `gippa setup` wall time for the
  default set. The quick start gives only the measured parts.
