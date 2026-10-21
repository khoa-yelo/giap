# GIAP progress

> **Resuming after a timeout:** start an allocation, run `bash scripts/resume_status.sh`, then use the
> resume prompt (GIAP_SPEC.md §0.3). Every worktree's next feature-loop step is printed; uncommitted
> files mean an agent stopped mid-step — inspect them, then re-dispatch that step (test-author →
> implementer → code-reviewer). Slurm jobs survive the session; their logs are in `local/logs/`.
_Last updated: 2026-10-20 16:50 PDT · session https://claude.ai/code/session_01FctC7kZEMou2AYwYy1KhJa_

## Milestone progress
<!-- progress:start -->
🟩 done (merged) · 🟨 running (open branch) · ⬜ pending; designed issues: 84/91 done (92 %), 2 running · generated 2026-10-21 10:35 by `scripts/dev/progress_bars.py`

| Milestone | Progress | Done | Running | Pending |
|---|---|---|---|---|
| M0 Scaffold | 🟩 | 1 | 0 | 0 |
| M1 Inputs and preflight | 🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩 100 % | 12 | 0 | 0 |
| M2 Annotation wrappers | 🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩 100 % | 21 | 0 | 0 |
| M3 Mapping and decoy FDR | 🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩 100 % | 25 | 0 | 0 |
| M4 Gene table and region tests | 🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟨⬜ 90 % | 19 | 1 | 1 |
| M5 Autocorrelation and islands | 🟩🟩🟩🟩🟩🟩🟨⬜⬜⬜⬜ 55 % | 6 | 1 | 4 |
| M6 Operons | _not designed yet_ | - | - | - |
| M7 Variants and novel ORFs | _not designed yet_ | - | - | - |
| M8 Function and determinants | _not designed yet_ | - | - | - |
| M9 Visualisation and report | _not designed yet_ | - | - | - |
| M10 Pilot validation | _not designed yet_ | - | - | - |
| M11 Release 1.0 | _not designed yet_ | - | - | - |
<!-- progress:end -->

## Current milestone
M0–M1: done. M2: done (21/21; M2-21 closed it, #107).
M3: 22/24 — M3-24 needs DIA-NN from Khoa; ADR-0019 implementation is the last code item.
M4: 19/21 after #117 — M4-20, M4-21 not started.
M5 (autocorrelation and islands): designed (#108, 11 issues); M5-01 merged, M5-02/03/04 in review.
M6–M11: not designed yet.

## Done (with evidence)
- [M0] `gh` authenticated as khoa-yelo (SSH protocol); repository is public (Actions minutes unlimited).
- [M0] Repository connected to git@github.com:khoa-yelo/giap.git; branch `m0/scaffold` pushed.
- [M0] Bootstrap ran (`scripts/bootstrap_farmshare.sh`): uv 0.13.1, Nextflow 26.10.0, nf-test 0.9.5,
  gh 2.102.0, OpenJDK 21, all under `local/` (ADR-0002).
- [M0] `uv.lock` resolved (131 packages); venv in `local/venv`.
- [M0] Probe in the session allocation (iron-01) and as Slurm job 1845260 (`giap-probe`, wheat-08,
  COMPLETED in 31 s): batch path proven.
- [M0] Appendix A files, ADR-0001 (imports §2.1/§2.4), CI guards (`scripts/ci/check_*.py`), `ci.yml`,
  capability-marker plugin with `--require-capabilities` (tests/unit/test_capabilities.py, 2 passed).

## In progress
Session 2026-10-20 (lead + up to 6 subagents in parallel, Khoa's instruction). Merged today: #107 M2-21
(nightly array, weekly, measured runtimes), #108 M5 design, #109 one-subprocess step-load test (52.6 s →
2.7 s), #110 M4-18 regions stats, #111 M4-15 gene table, #112 M5 decisions + ADR-0020/0021, #113 M3-21
calibration suites, #114 M5-01 Blom/LOESS, #115 M3-22 digest + TE1/TE8.

Open PRs: #116 M5-02 (autocorr engine; code + statistician approved, CI), #117 M4-19 sensitivity (both
approved, CI), #118 M5-03 scan (code approved; statistician asked for the islands methods note —
implementer finishing), #119 M5-04 simulator oriC/gradient/C15 (code review running), #120 chore: four
wall-clock asserts moved to benchmarks/, #126 M5-05 autocorrelation table (implementation + benchmark
baseline; awaiting code and statistician review).

M5-05 (#126, implementer step done): fast suite 8668 passed (Slurm 1861597), calibration fast tier
including all four attack arms 23 passed (1861598), **full tier R = 1,000 at 4,500 genes 12 passed in
1:21** (1861606). Lead decision 34 is resolved: `tu_null_4500.morans_i.tu_block.lag1` is 0.113 against a
bound of 0.119 at α 0.10 (0.055 / 0.011 at 0.05 / 0.01), all on block 1, reproduced alone in 1861757 —
the fast-tier `unconfirmed` (0.19 then 0.13 against 0.16) is the 1-SE noise decision 34 expected, no
bound touched. Benchmark baseline `autocorrelation_rows_b1000_n4500_20_lags` = 2.07 s of a 30 s budget
(wheat-07, 4 CPUs, Slurm 1861601).

M5-08 (#132, implementer step done): `giap.context.oriter` — oriC resolution (circular mean inside
0.05·L, `multiple_oric` event on both branches), replichore geometry, the two genes.tsv columns and the
three descriptive OLS slopes with a circular block-bootstrap CI; block length clamped per fit with that
fit's own T (lead decision 1 of 2026-10-21), ℓ̂ once per replicon from the pooled residuals. Unit +
property tests 42 passed, module coverage 97 %; fast suite **8872 passed, 29 skipped** (Slurm 1863080,
6:37); the recorded fast coverage arm (example 7, 2 log2, ρ 0.5, R = 50) **coverage 0.96, estimable 1.00,
mean CI width 1.11** (Slurm 1863081, 27 s) — the test-author's 0.96 reproduced. Plus
docs/methods/oriter.md. `features` stays a plain parameter; the step plumbing is M5-09's.
M5-06 (**PR #131**, implementer step done, 2026-10-21): `giap.context.islands` on the merged M5-01/02/03
engines — one **genome-wide** M_b over the scanned replicons (decisions 21/46), ℓ' = min(20, 2ℓ̂+1, T),
the best window as the reported span (decision 42), the `flags` column with `ci_few_units` only, and
ADR-0021 Amendment 1 (no `uncalibrated`; `IslandScan` carries ℓ̂, ℓ' and the lag-1 Moran's I per replicon
for the §9.4 run report). Unit 34 + props 5 + sim 6 (slow) passed; islands calibration fast tier
4 passed (ρ 0 → 0.000, ρ 0.5 → 0.020, chromosome + plasmid → 0.020 at α 0.05, bound 0.112 at R = 50 —
equal to the test-author's pre-implementation oracle), teeth confirm the spec's ℓ̂ at **0.72**; the
non-nightly attack arms 6 passed (the decision-50 ceiling 3.35–3.42 against M_b q95 4.49–4.92, best low
p 1.000 at 10/20/40 silenced TUs). CI all green: 8760 passed, coverage 96.25 % (islands.py 95 %).
Benchmark baseline `islands_4500_genes_b1000` = 1.79 s of a 20 s budget (4,500 genes, 2,482 TUs,
B = 1,000 + 1,000). **No Slurm run**: the shared account sat at its 512-CPU QOS cap all day
(`QOSMaxCpuPerUserLimit`, 32 × 16-CPU non-GIAP `n1170_tools` array), so jobs 1863043/1863044 never
started and were cancelled; CI covers the fast suite, the fast calibration tier and the coverage gate.
Open for review: the M5-06 section of `docs/methods/islands.md` (bootstrap null, genome-wide family,
Amendment 1's validity domain, the low-direction ceiling, effect CI, degenerate cases) is not written
yet — `docs/methods/` is the statistician's.

M6-11 (implementer step done, 2026-10-21): RegulonDB is **fetched and the K-12 truth sets are built**.
`validation/fetch_regulondb.py` completes the TLS chain itself — it downloads the GlobalSign
intermediate RegulonDB omits, refuses it unless its SHA-256 equals the pin
`d160e2de…63f9` (GlobalSign GCC R46 OV TLS CA 2025, valid to 2029-06-23), checks with
`openssl verify` that it is signed by `GlobalSign Root R46` (already in the trust store), and only then
adds it to the CA bundle; verification is never disabled (a test greps the source for `verify=False`,
`CERT_NONE` and friends). Real fetch, twice, identical hashes: **RegulonDB 14.5.0** (release 2026-01-28,
EcoCyc 29.5, genome NC_000913.3) into `$GIAP_DATA_DIR/regulondb/14.5.0/` — `TUSet.txt` 446,341 B
(sha256 c96c8a3b…), `TerminatorSet.txt` 60,282 B (9d916ffe…), `GeneProductSet.txt` 1,290,719 B
(3c2c6be3…), plus `provenance.json`. Nothing committed (R6/R11); no Slurm job (the whole download takes
about two seconds, as the spec predicted).

**Spec example 5, the evidence M6-12 consumes** (after the #141 review fixes, lead decisions 1–8;
`validation/regulondb_truth.py` on the recorded `tests/data/recorded/bakta/1.12.2-full/k12/` JSON,
contig U00096.3, ±30 bp, every file verified against `provenance.json` first): **3,067 same-strand
adjacent pairs** (the pipeline's `adjacent_pairs`, **wrap pair included**: yjtD → KOIAIB_00001, labelled
ambiguous because KOIAIB_00001 is a Bakta-only sORF at 3–98) — **1,730 positive, 1,009 negative, 2,739
labelled** (≥ 2,000 as required), 328 ambiguous and excluded (10.7 %); **375 terminators**. Gene mapping
**4,392 of 4,745 RegulonDB genes = 0.926 overall** (4,217 onto Bakta CDS + 175 RNA genes onto Bakta
ncRNA/tRNA/rRNA features); **CDS share 0.923** (4,217 / 4,570); both above the 0.8 escalation floor;
353 unmapped. dnaK–dnaJ, ftsA–ftsZ and pdhR–aceE are now positive (decision 3). Two RegulonDB names
occur twice and are kept per gene id (decision 5). The summary carries release 14.5.0 and the three
SHA-256s. Unit tests 91 passed, nothing touches the network. Release of record: 14.5.0, not v12
(decision 8; PROGRESS › Questions for Khoa).

ADRs accepted today by Khoa: ADR-0019 (stratum-0 π̂ per length class, Clopper–Pearson upper bound, no PEP),
ADR-0020 (autocorrelation: re-detrend per permutation, residual join count, CI NA) **+ Amendment 1**
(Pearson residual (d−p)/√(p(1−p)); the raw residual rejected at 0.074–0.105 when low p_detect clusters),
ADR-0021 (block length ℓ' = min(20, 2ℓ̂+1, T)) **+ Amendment 1** (the validity domain is a dependence
regime, not a replicon size; `uncalibrated` dropped — the §9.4 run report carries ℓ̂, ℓ' and the lag-1
Moran's I per replicon instead).

Infrastructure: scrontab reinstalled from scripts/slurm/nightly.scrontab (01:30 launcher → array 0-3 +
summary with afterany); first real-cluster dry run all green (1860190–1860196). Hold watchdog running.
CI ≈ 20 min.

M3-25 (ADR-0019 implementation) is **PR #124** (implementer done, CI running). Slurm: fast suite 1861427
(6 failures, all node-load subprocess timeouts that pass locally, plus the golden `mapping_fdr.tsv`
regenerated on purpose — the 50-kb fixture is all `same_db`, so π̂ = c = 1 on empty lists); fast-tier
placement calibration 1861428 (41 passed); the two flipped `calibration_full` estimate checks 1861430
(both pass: related_dbstrain Ŵ0 1.551 vs V0 0.79, ΣV0 = 79 ≥ 30; db_stressed Ŵ0 4.306 vs V0 2.75,
ΣV0 = 275 ≥ 100). Round 2 (decisions 11–14 + the statistician's extras): config guard on
`mapping.length_classes`, named ValueErrors in `stratum0_pi_classes`, per-class `ci_high` from the
Clopper–Pearson `upper`, and the methods notes (sign-offs recorded, post-ADR numbers, deflated-PEP
caveat, "read π̂_ℓ from `stratum0_classes`"). Fast suite **1861618: 8786 passed, 0 failed**; every CI
check green (guards included, after the lead's EXEMPT_TESTS line). **Left for a lead/test-author
commit (R2 — the implementer may not touch existing or new test files):** decision 11's deletion of
`pi_wrong_database` (six unit tests in tests/unit/test_placement_fdr.py still import and pin it; the
function is dead in `src/` and its docstring says so), decision 12's two edge-case tests (written and
mutation-checked, patch at
`/tmp/claude-430805/-home-classes-bios-270-khoa-giap/b9f44921-ae38-43a6-9b00-e1601c824b49/scratchpad/decision12-edge-tests.patch`),
and the statistician's `write_diagnostics` extra (sketch beside it as
`statistician-diagnostics.patch`).

Worktrees in flight: adr0019 (ADR-0019 implementation, PR #124), m5-scan (#118), m5-sim-oric (#119),
m5-moran (#116), sensitivity (#117), timing-bench (#120), digest-te/gene-table/regions-stats/calibration-gate
(merged, removable).

Follow-ups queued (docs/backlog.md has the full list):
- Per-process real-tool Slurm job array for workflow modules (inputs from recordings) — after #65 (Khoa).
- calibration_full variants for migrated suites (statistician); GO data-version pin (M2-09).
- M3-21: c_ref − c_true, chance rate of wrong peptides, PEP–length correlation (ADR-0010).
- M4-12/M4-17 forwarded items; M6: simulate readthrough (cterm_extension).
Lessons (keep):
- A merge loop must not merge right after pushing a main merge (stale green snapshot); merge_queue now waits
  until the PR head contains origin/main and CI re-ran.
- Launch failures came from SLURM_GET_USER_ENV/SLURM_EXPORT_ENV inherited from the interactive allocation
  (#102), not from the nodes.
- Lead decisions on stats rules need the statistician's attack check: decision 11 ("any") was gameable.
- Hypothesis found a separation case (no MLE) the spec never covered: specs must say what happens when an
  estimator does not exist.
- A PR that conflicts with main gets NO CI run; waitmerge now stops on it (#90). Merge main, then retry.
- Fixture genomes need `[topology=circular]` (spec §~770); plain NCBI FASTA → linear replicons. Nissle
  GCA_043228065.1 has 3 plasmids (p15A, pMUT1, pMUT2).
- Generated golden outputs (tests/data/golden/) are exempt from the role check, like recorded fixtures.
- Shared test files across tool-layer PRs: scratchpad merge_main_nf.sh resolves BUILDERS + step_registry.
- Real `giap run` (SIM100K, farmshare_local,apptainer, real $GIAP_DB_DIR): preflight ok, BAKTA (13.5 min),
  GENOMAD, TERMINATORNET ok, `--resume` caches; stops at the first unregistered py step (codon_usage, M2-17).
- Tests built on invented fixtures hid contract bugs (M4-13, M4-14): fixtures must come from the real
  producers' outputs. Test module basenames must be unique (scripts/ci/check_test_names.py in CI).
- Slurm "user env retrieval failed / requeued held" happens repeatedly on batch nodes: release held giap-*
  jobs (`scontrol release <id>`); a 2-min watchdog does this during sessions.
- Shared files across tool-layer PRs (tests/unit/test_nf_cmd_sync.py BUILDERS,
  workflow/tests/data/annotate/step_registry.json): merge main, keep both builders, regenerate the registry.
- Never rebase + force-push: merge origin/main into the branch and push plainly (force-push is blocked;
  Khoa, 2026-10-15). Recover a locally rebased branch by cherry-picking onto origin's branch.
- Slurm work dirs live under $GIAP_ROOT (scratch/…), never node-local /tmp or the session scratchpad.
- Rebuilding a shared micromamba env from one branch breaks `test_envs_built` on all other branches until
  that branch merges (promotech/PyQt5, #65).
- Every tool wrapper gets a real-tool run before code review; stub fixtures hid ~10 bugs in M2-07.
- Keep ≤ 128 CPUs: Slurm rounds each job up by 2 CPUs; chain extra jobs with --dependency.
- Slurm jobs run from the main checkout only; never remove a worktree a job may be using.
- Subagents never install packages; set PYTHONPATH per worktree; unset PYTHONHOME after sourcing the env.

## Next actions (ordered)
1. Merge #116, #117, then #118, #119, #120 (waitmerge; merge origin/main into each branch first —
   #116 and #120 both touch benchmarks/baseline.json).
2. ADR-0019 implementation (branch m3/adr0019-pi-hat): class estimator, detection_summary QC rows,
   `input_prefiltered` warning, and flip the two strict xfails (1-Mb related_dbstrain estimate check and
   test_database_stressed_stratum0_estimate_covers_false_placements; check the ΣV0 floors still hold,
   ≈ 285 expected). Then M3 is done except M3-24.
3. M5-05 (autocorrelation tables) and M5-06 (islands) test-authors after #116/#118 merge; then M5-07,
   M5-08 (needs M5-02 + M5-04), M5-09.
4. M4-20, M4-21.
5. Design M6 (operons).

## Blocked
| Since | What | Host or tool | Needed from Khoa |
|---|---|---|---|
| 2026-10-14 | ~~RegulonDB v12 download~~ **RESOLVED 2026-10-21** | Not an outage: RegulonDB serves a cert (GlobalSign GCC R46 OV TLS CA 2025, issued 2026-10-14) **without the intermediate**, so verification fails. The intermediate from GlobalSign's own CA endpoint verifies against `GlobalSign Root R46`, already in FarmShare's trust store; with it, the site and `/datasets` return HTTP 200. Khoa approved completing the chain (2026-10-21). | nothing — M6-11/M6-12 unblocked |

## Answered by Khoa
- 2026-10-15: force-push blocked → merge main into feature branches instead of rebasing (option 3).
- 2026-10-15: light mode (no new agents) 01:30–03:30 for the usage reset; per-process Slurm testing for M2-07 wanted.
- 2026-10-14: ADR-0013 cterm_extension, ADR-0015 two-block confirmation, ADR-0016 size gate excludes workflow/tests (all accepted-by-khoa).
- 2026-10-14: novel ORFs are checked for a nearby promoter and terminator — evidence columns, flags and a
  ranking key, not a calling requirement (ADR-0012, accepted-by-khoa; implemented in M6).
- 2026-10-14: ADR-0011 (`low_complexity` placement flag) accepted as proposed.
- 2026-10-14: M3-13 tblastn HSP chaining limits (max gap 300 nt, query overlap ≤ 10 aa) are fixed
  constants documented in the methods note — no config keys.
- 2026-10-14: replace the slow Bakta full-DB download with a resumable one (done, #33).
- 2026-10-14: §7.3.2 picked FDR rescales decoy wins by the paired share (ADR-0009, accepted-by-khoa).
- 2026-10-14: subagent commit blocked by the permission system → Khoa approved the lead committing the M2-02 test fix.
- 2026-10-14: ADR-0005 item 4 (new preflight warning thresholds) rejected: decoy share and free RAM/disk are Info only.
- 2026-10-14: ADR-0004 (M0 size exemption; hook realpath fix) and ADR-0007 (size guard counts code only) accepted; "merge the stack".

## Questions for Khoa
| # | Question | Options | Default if no answer | Asked |
|---|---|---|---|---|
| 1 | GIAP_SPEC.md §15.3/§15.4 name "RegulonDB v12", but v12.0.0 is the 2023-10-23 release and the live one is **14.5.0** (2026-01-28). The fetcher records the release the endpoint reports. | (a) use the current release (14.5.0) and cite it in ADR-0022; (b) pin v12 and ask RegulonDB for the archived files | (a) — ADR-0017: record the version and the hash of what we actually fetched | 2026-10-21 (M6-11) |

## Environment notes
- capabilities: allocation cpus=8, mem=31.2 GB, time left at 16:32 UTC=4:18; local/ used=1.75 GB;
  apptainer=1.4.4, micromamba=2.3.3, nextflow=26.10.0, nf-test=0.9.5, gh=2.102.0, sbatch=ok,
  gpu partition=visible
- databases present: Bakta light+full v6.0, geNOMAD, Promotech RF-HOT, antiSMASH (4.1 GB); genomes K-12,
  Nissle; UP000000625 proteome
- sources: all §3.2 sources reachable except regulondb.ccg.unam.mx

## Measured numbers
| Step | Genome | Wall | Peak RSS | Node class (Slurm job ID) |
|---|---|---|---|---|
| bootstrap (micromamba env giap-wf) | — | ~3 min | — | iron-01 (session) |
| uv sync --all-extras | — | 1 min 48 s | — | iron-01 (session) |
| Bakta full DB | Nissle 1917 | 13.5 min | — | wheat (1846846) |
| TerminatorNet (-n 17) | K-12 MG1655 | 60 min | — | wheat-04 (1847778) |
| giap setup GO/KEGG/COG/contaminants | — | — | go-basic 32.2 MB, KEGG links 2.3 MB | session (review #67) |
| Bakta 1.12.2 full DB 6.0, circular-tagged, `--replicons` (M2-10) | K-12 MG1655 | 9.9 min (592 s) | 20.6 GB (job, with SIM100K) | barley-01, 16 CPU (1849096) |
| Bakta 1.12.2 full DB 6.0, circular-tagged (M2-10) | SIM100K | 5.3 min (317 s) | (same job) | barley-01, 16 CPU (1849096) |
| Bakta 1.12.2 full DB 6.0, circular-tagged (M2-10) | Nissle 1917 | 9.9 min (595 s) | 20.5 GB | barley-01, 16 CPU (1849097) |
| Bakta 1.12.2 light DB 6.0 (M2-10) | SIM100K | 34 s | 5.9 GB | barley-01, 16 CPU (1849098) |
