# M10-04: pass-1 pilot runs, `giap run` A1 and A2 and `giap compare` (the MVP closer)

Milestone M10 · **pass 1** (MVP) · depends on: M10-03 (`$GIAP_ROOT/pilot/inputs/a{1,2}.tsv.gz`
written by `diann_convert.sbatch` with `convert.json` `failures: []`); MVP-04's acceptance run
(soft: see "Order of work") · estimate ≈ 85 code lines (index: 80) + ≈ 200 test lines · not `stats` ·
roles: validation (`validation/pilot/run_evidence.py`, README section), implementer
(`scripts/slurm/pilot_compare.sbatch`), test-author (CI tests on hand-made run directories), lead
(submits the jobs, R14 budget, PROGRESS) · Khoa: none blocking (one parked question, Q-K1).

## Spec sections implemented

- §15.2 step 8, the A1/A2 half: "Run `giap run` for A1, A2 and A3 as batch jobs; then
  `giap compare out_a1 out_a2 out_a3 --labels A1,A2,A3`". A3 is M10-10; this issue compares A1 and A2.
  `REPORT.md` is M10-05.
- §15.2 "GIAP runs" block (genome, annotation, peptides, default `--genome-relation related`).
- §15.6 bullets 1 (A1 and A2 part: "GIAP runs end to end") and 3 ("A1 and A2 detections and
  abundances are compared (`giap compare`)"). Evidence table in the index.
- §15.5 metric 7: runtime and peak memory per arm and per step (`OUT/logs/trace.tsv`, measured with
  `scripts/measure_runtimes.py`). Recorded here; M10-05 tabulates it.
- §7.9.3 `giap compare`; ADR-0046 (locus-tag join, genome identity, exit 1 on a different genome);
  ADR-0042; M9-22 (`compare_*.tsv`) and M9-23 (V16, `report/compare.html`).
- §9.1 output tree; §5.1 exit codes (0 success, 1 invalid input, 2 environment, 4 step failed,
  5 internal).
- §3.8 (batch jobs, job IDs in PROGRESS, polling); R9 (three attempts); R11 and Khoa's 2026-10-23
  rule (nothing derived from the DIA-NN reports is committed); R14 (CPU cap).
- Index lead decisions 1 (L1: one arm at a time; L10: the function module on) and M10-03 lead
  decision 1 Q4 (no `--abundance`) and Q5 (`search_db` already set by the converter).

## What exists (main 3b2f1ffa)

| Piece | State |
|---|---|
| `scripts/slurm/giap_run.sbatch` | head job `giap-run`, `normal`, 2 CPU, 8 GB, 2 days; `cd $GIAP_REPO` (default: the **main checkout**) and `uv run --no-sync giap run --profile farmshare,apptainer "$@"`. Its arguments are passed verbatim, so `$SLURM_JOB_ID` cannot appear in `--outdir` |
| `giap run` | `--threads` (default 16), `--memory` (default `64.GB`), `--modules` (default `core,mobile,bgc,function`, `defaults.yaml`), `--seed` (default 20261014), `--resume`, `--db-dir`. Exit codes per §5.1; `run_manifest.json` has `status`, `exit_code`, `failed_step`, `git_commit`, `git_dirty`, `inputs[]` (sha256), `annotation`, `databases`, `detection_tier`, `seeds`, `steps`, `modules_off`, `warnings`. Nextflow's work dir is `OUT/work`; the trace is `OUT/logs/trace.tsv`; the report is `OUT/report/report.html` (re-rendered after the manifest, M9-16b) |
| `giap compare RUN… --outdir CMP --labels …` | exit 1 (`CompareInputError`) on fewer than two runs, a label count mismatch, a run with no genome sha256, **different genome digests**, no `tables/genes.tsv`, or different locus-tag sets. Writes `CMP/tables/compare_{detection,abundance,stability}.tsv`, `CMP/plots/V16_compare.*`, `CMP/report/compare.html`. Bootstrap through `giap.rng.for_step("compare", seed)` |
| `giap.io.output_checks.check_tsv(key, path, stem)` | checks a TSV against `giap/schemas/tables/<stem>.schema.json` (columns, types, `required`, `enum`, `minimum`, `pattern`); raises `ValueError` |
| `workflow/conf/farmshare.config` | `queueSize = 8`; per task `cpus = max(label cpus, ceil(memory / 4 GB))`, clamped by `resourceLimits.cpus = min(16, threads)` and `memory = min(64 GB, --memory)`. Its comment says FarmShare enforces `MaxMemPerCPU = 4000 MB` |
| Genome | `tests/data/genomes/nissle1917.fna.gz` (GCA_043228065.1, sha256 `fe748016…23ee`): `CP171242.1` chromosome 5,051,730 bp and **three** plasmids `CP171243.1` (4,523 bp), `CP171244.1` (3,173 bp), `CP171245.1` (5,514 bp), all `[topology=circular]` |
| Annotation | `tests/data/recorded/bakta/1.12.2-full/nissle/nissle1917.json.gz` (sha256 `a012798…b8a7`; Bakta 1.12.2, full DB 6.0; 4,982 `JMJJOG_` locus tags). `validation/pilot/db_a2/strain_db.faa.gz` (4,692 proteins) was built from this JSON, so A2's protein IDs are its locus tags |
| Peptides | `$GIAP_ROOT/pilot/inputs/a1.tsv.gz`, `a2.tsv.gz` (M10-03; not yet produced). `search_db` is `other` (A1) and `same` (A2) in the table |
| Databases | `$GIAP_DB_DIR` holds `contaminants/`, `cog/`, `go/`, `kegg/`, `genomad_db/`, `promotech/`, `antismash/`, `bakta/` |
| Real-genome run times | K-12 real-tool run (M7-15): Slurm 1873779 COMPLETED in **5 h 59 min** (head), after 1871633 FAILED at 3 h 34 min; 2,647 tasks. MVP-04 expects 3–8 h. The index's "35–60 min per arm" (§13.2, 16 threads, standard run) is not what a full real-tool run measures today |

### Gap found while writing this spec: `--annotation *.json.gz`

`cli._bakta_json` sets `annotation.precomputed.bakta_json` only when the `--annotation` suffix is
`.json`. The recorded JSON is `.json.gz`, so passing it as in §15.2 would leave the precomputed key
unset. The Nextflow path would then treat it as a user annotation and **run Bakta live** (and the
`--executor local` path exits 1). The run would still finish, but on a fresh Bakta annotation, not
the recorded one, and the pilot inputs would differ from the README. **Default (Q3):** this issue
decompresses the JSON once into `$GIAP_ROOT/pilot/inputs/nissle1917.bakta.json` (mode 0444, sha256
recorded) and passes that `.json`. A two-line CLI fix (`.json.gz` accepted) goes to the backlog as
its own PR. The evidence check (rule A4 below) fails any arm whose trace contains a `BAKTA` task.

## Deliverables

| # | File | Role | What |
|---|---|---|---|
| 1 | `validation/pilot/run_evidence.py` | validation | pure functions `arm_evidence(outdir, sacct_rows, expect) -> dict` and `compare_evidence(cmp, cmp2, labels) -> dict`; thin CLI `arm` / `compare`; writes JSON, prints one markdown row per subject (≈ 60 lines) |
| 2 | `scripts/slurm/pilot_compare.sbatch` | implementer | runs `giap compare` twice (into `CMP` and `$TMPDIR/cmp2`) and `run_evidence.py compare` (≈ 25 lines) |
| 3 | `validation/pilot/README.md` › "Pass-1 runs (M10-04)" | validation | commands, job IDs, exit codes, `sacct` State/Elapsed/MaxRSS/NodeList, the main commit, input sha256s, `report.html` and `compare.html` sha256s and byte sizes |
| 4 | `PROGRESS.md` | lead | job IDs on submission; Measured numbers rows (per arm head wall, the top steps from `measure_runtimes.py`); Done entry with the evidence paths |
| 5 | tests (below) | test-author | CI only, on hand-made run directories |

No change to `src/giap`. The arm runs use the existing `giap_run.sbatch`, unchanged.

### What is committed (and what is not)

| Item | Where | Committed |
|---|---|---|
| run directories `out_a1`, `out_a2` (tables, plots, `report.html`, work dirs) | `$GIAP_ROOT/runs/pilot/pass1/` | **no**: every table is derived from the DIA-NN peptides (Khoa's rule), and the work dirs are tens of GB |
| comparison directory `cmp_a1_a2` | same | **no** (derived) |
| evidence JSON `evidence/{a1,a2,compare}.json` | same | **no**: it holds paths, hashes, exit codes and resource numbers only, but it lives with the runs |
| decompressed annotation `nissle1917.bakta.json` | `$GIAP_ROOT/pilot/inputs/` | no (the `.json.gz` is already committed) |
| README section (deliverable 3), `run_evidence.py`, the sbatch, tests | repository | **yes**. Only provenance: IDs, exit codes, states, runtimes, memory, sha256s and byte sizes. **No biological number** (no counts of detected genes, no ρ, no Jaccard): those are M10-05's, and whether they may be committed is Q-K1 |

## Layout under `$GIAP_ROOT/runs/pilot/`

```text
$GIAP_ROOT/runs/pilot/
  pass1/                      this issue (MVP evidence; M10-14's report-freeze runs go to freeze/)
    out_a2/                   first attempt of A2;   out_a2.r2/, out_a2.r3/ for fresh re-runs
    out_a1/                   first attempt of A1;   out_a1.r2/ …
    cmp_a1_a2/                giap compare out of the two accepted attempts
    evidence/a1.json, a2.json, compare.json
    ATTEMPTS.tsv              one row per submission: arm, attempt, outdir, job, resume, git_commit, state, exit_code
```

Names follow §15.2 (`out_a1`, `out_a2`) so the commands read as the spec's. The job ID is not in
the name: `giap_run.sbatch` passes its arguments verbatim, and a `--resume` resubmission must reuse
the same directory. `ATTEMPTS.tsv` maps each directory to its jobs. `.r<n>` suffixes appear only for
fresh attempts (Failure handling). `pass1/` keeps these runs apart from M10-06's second A2 run (TE7)
and M10-14's freeze runs.

## The runs

### Before submitting (lead, every time)

1. **Inputs.** `a1.tsv.gz` and `a2.tsv.gz` exist, and their sha256 equals the one in
   `a{1,2}.convert.json` (M10-03). The decompressed annotation exists with its recorded sha256.
2. **Code.** The main checkout is fast-forwarded to `origin/main`, includes M10-03, and is clean
   (`git status --porcelain` empty). `giap_run.sbatch` runs that checkout. Its short sha is written
   to `ATTEMPTS.tsv`, and `run_manifest.json` `git_commit` must later equal it.
3. **R14 budget.** S = the CPUs of all running and pending `giap-*` jobs of every agent
   (`squeue -h -u $USER -n` over `giap-*` names, `%C`). One arm at the settings below needs at most
   8 × 14 + 2 = **114** CPUs. Submit only if S + 114 ≤ 128, that is S ≤ 14. Otherwise wait and
   re-check (polls ≤ 10 min apart). Do not lower `--threads` to fit (Q1).
4. **Order.** A2 first, then A1 after A2's head job ends (any state). Never both at once (index
   decision 3: 2 × 114 > 128). A2 first because it is the best case and feeds TE1 (M10-05): a
   pipeline bug shows on it as well as on A1, and its result is the one M10-05 needs first.

### Commands

```bash
P=$GIAP_ROOT/pilot/inputs; R=$GIAP_ROOT/runs/pilot/pass1
G=$GIAP_REPO/tests/data/genomes/nissle1917.fna.gz
gzip -dc $GIAP_REPO/tests/data/recorded/bakta/1.12.2-full/nissle/nissle1917.json.gz > $P/nissle1917.bakta.json
chmod 0444 $P/nissle1917.bakta.json; sha256sum $P/nissle1917.bakta.json      # once; recorded in the README

sbatch -J giap-run-pilot-a2 scripts/slurm/giap_run.sbatch \
  --genome $G --annotation $P/nissle1917.bakta.json --peptides $P/a2.tsv.gz \
  --db-dir $GIAP_DB_DIR --threads 14 --memory 54.GB --outdir $R/out_a2
# after A2's head job has ended:
sbatch -J giap-run-pilot-a1 scripts/slurm/giap_run.sbatch \
  --genome $G --annotation $P/nissle1917.bakta.json --peptides $P/a1.tsv.gz \
  --db-dir $GIAP_DB_DIR --threads 14 --memory 54.GB --outdir $R/out_a1
# after both arms pass rule set A:
sbatch -J giap-pilot-compare scripts/slurm/pilot_compare.sbatch $R/out_a1 $R/out_a2 $R/cmp_a1_a2
```

Settled choices:

- **`--threads 14 --memory 54.GB`** (Q1). `--threads 14` alone does not bound the run at 114
  CPUs: a task whose memory is clamped to 64 GB at 14 CPUs asks 4,681 MB per CPU, above FarmShare's
  `MaxMemPerCPU = 4000 MB`, and Slurm then raises its CPU count to 16 or more. 54 GB is 55,296 MB,
  or 3,950 MB per CPU at 14 CPUs. Every Nissle task measured so far fits: Bakta full is not run,
  DETERMINANTS asks 48 GB, and the K-12 tasks peaked at 7.5 GB. The evidence (rule A6) checks the
  largest `AllocCPUS` of any task job.
- **Modules: the default `core,mobile,bgc,function`** (index L10: function on). No `--modules`
  flag. The run uses `$GIAP_DB_DIR`'s GO, KEGG and COG; the KEGG retrieval date comes from the
  manifest `databases` (ADR-0017). No KEGG-derived content is committed (none is: rule above).
- **Seed: the default 20261014**, the same for both arms (recorded from `seeds`). No `--seed`
  flag.
- **No `--abundance`** (M10-03 Q4): gene abundance comes from peptide intensities.
  **No `--reference-peptides`** (that is A3's).
- **`--genome-relation related`** (the default; §15.2) for both arms. A1 searched K-12, so many
  of its peptides are unmapped or tolerant; that is a result, not an error.
- **No `--strict`**: warnings are recorded, not fatal (definition of done).

### The compare job (`scripts/slurm/pilot_compare.sbatch OUT_A1 OUT_A2 CMP`)

| `-J` | Partition / QOS | CPUs | Mem | Time limit (expected) |
|---|---|---|---|---|
| `giap-pilot-compare` | normal / normal | 2 | 8 GB | 01:00:00 (minutes: 4,982 genes, 1,000 bootstrap draws) |

Header as `scripts/slurm/_header.sbatch`, `set -euo pipefail`, the environment sourced, `set -f`.
In order:

1. Exactly three arguments; `OUT_A1` and `OUT_A2` are directories, and `CMP` does not exist
   (a re-run uses a new name or the lead removes the old one: no merging into stale output). Else
   exit 2.
2. `giap compare OUT_A1 OUT_A2 --outdir CMP --labels A1,A2`, then the same into
   `$TMPDIR/cmp2`. giap runs as `"$GIAP_PYTHON" -c 'from giap.cli import cli; cli()'`, with
   `PYTHONPATH=$GIAP_REPO/src` and `GIAP_PYTHON` defaulting to the shared venv's python. No `uv run`.
   A non-zero exit of the first call is the job's exit code (1 for `CompareInputError`).
3. `"$GIAP_PYTHON" validation/pilot/run_evidence.py compare CMP --second $TMPDIR/cmp2 --labels
   A1,A2 --out <R>/evidence/compare.json`. Its exit code is the job's.

## `run_evidence.py`

### `arm OUT --arm A1|A2 --job JOBID --expect-commit SHA --out FILE`

The CLI calls `sacct -P -n -j JOBID -o JobID,State,ExitCode,Elapsed,MaxRSS,NodeList` for the head
job. It also queries the task jobs: the `native_id` column of `OUT/logs/trace.tsv`, with
`-o JobID,AllocCPUS,MaxRSS,State`. Then it calls the pure `arm_evidence`. The JSON holds:

- `arm`, `outdir`, `job` (head: id, state, exit code, elapsed, MaxRSS, node);
- `manifest`: `status`, `exit_code`, `failed_step`, `failed_log`, `git_commit`, `git_dirty`,
  `profile`, `resume`, `started_at`, `finished_at`, `seeds`, `modules_off`, `detection_tier`,
  `annotation` (version, DB version, DB type), `databases` (name, version, date), `warnings`;
- `inputs`: role → sha256 from the manifest;
- `report`: path, sha256 and bytes of `report/report.html`;
- `tables`: the §9.1 `tables/*.tsv` present, each with sha256 and bytes (no content);
- `tasks`: count, failed count, retried count, `max_alloc_cpus`, peak `MaxRSS`;
- `failures`: list of rule strings (below). The exit code is 0 when `failures` is empty, else 1.

**Rule set A (an arm passes when every rule holds):**

| # | Rule | Failure string |
|---|---|---|
| A1 | `run_manifest.json` exists and validates against `run_manifest.schema.json` | `no_manifest` / `manifest_invalid:<path>` |
| A2 | `status = completed`, `exit_code = 0`, the head job `COMPLETED` with `ExitCode 0:0`, and no step `failed` in `steps` | `run_failed:<failed_step or exit code>` |
| A3 | `git_commit` equals `--expect-commit`, and `git_dirty` is false | `commit_mismatch:<sha>` / `dirty_checkout` |
| A4 | the annotation is the recorded one: manifest `annotation.bakta_version = 1.12.2`, `db_type = full`, and no `BAKTA` task in the trace | `live_bakta` |
| A5 | the input sha256 for roles genome and peptides equal the expected files (`--expect-input role=sha`, given for genome, annotation and peptides) | `input_mismatch:<role>` |
| A6 | `max_alloc_cpus ≤ 14` | `cpu_cap:<n>` |
| A7 | `report/report.html` exists and is non-empty | `no_report` |
| A8 | `tables/genes.tsv` exists and passes `check_tsv("genes", …, "genes")` | `genes_table:<message>` |
| A9 | `modules_off` is empty (all four modules ran) | `module_off:<name>` |

A step that ran and produced an empty result is not a failure: for example no novel ORF at class
q ≤ 0.01, or no enriched category. Warnings are listed, never failed.

### `compare CMP --second CMP2 --labels A1,A2 --out FILE`

**Rule set C:**

| # | Rule | Failure string |
|---|---|---|
| C1 | `CMP/tables/compare_detection.tsv`, `compare_abundance.tsv` and `compare_stability.tsv` exist and each passes `check_tsv` against its schema | `table:<name>:<message>` |
| C2 | `compare_abundance` has exactly one row, `run_a = A1`, `run_b = A2` | `abundance_rows:<n>` |
| C3 | `compare_detection`'s `matrix` section has exactly 2 × G rows, where G is the row count of `out_a1/tables/genes.tsv`; each locus tag appears once per run | `matrix_rows:<n>` |
| C4 | `compare_stability` has one row per family in its schema enum, for the pair (A1, A2). A family with no item in either run is a row with `n_a = n_b = 0`, or the flag the schema defines; it is not a missing row | `stability_family:<name>` |
| C5 | `CMP/report/compare.html` and `CMP/plots/V16_compare.svg` exist and are non-empty | `no_compare_report` |
| C6 | the three tables in `CMP` and `CMP2` are byte-identical (same seed: §6.6, TE7 in small) | `nondeterministic:<name>` |

The JSON holds the tables' sha256 and bytes, `compare.html` sha256 and bytes, and `failures`. Like
the arm JSON, it holds no value from the tables (no ρ). The stability-family rule C4 is checked
against M9-22's actual behaviour before it is fixed in a test (Q6).

## What "done" means (MVP integration pass)

M10-04 is done when all of these hold:

1. Both arms have **one accepted attempt** that passes rule set A, after at most three attempts
   each (R9).
2. `pilot_compare.sbatch` exits 0 on the two accepted attempts (rule set C).
3. The README section, `ATTEMPTS.tsv`, PROGRESS (job IDs, Measured numbers, Done) and the
   evidence JSON paths are written. The PR description lists the evidence `failures: []` for all
   three subjects.
4. Every pipeline bug the runs exposed has a backlog entry or its own issue, **fixed in a separate
   PR**. This PR stays the evidence script, the sbatch and the README. A fix that needed a re-run
   is reflected in the accepted attempt's `git_commit`.

What done does **not** require, by design (MVP first, then hardening):

- no biological threshold: TE1, ρ, detection counts and Jaccards are M10-05's;
- no check of `report.html` content beyond existence: the report content PRs (M9-16, M9-17,
  M9-19/20, M9-26) are open, and the index says `REPORT.md` is written from tables;
- no runtime gate: the runtimes are recorded. A run longer than 8 h is a backlog item for M10-12/13,
  not a failure;
- no byte reproducibility of the arms (M10-06's TE7) and no frozen methods (M10-14). M6-12 C, M7
  and M8 hardening are open, so pass-1 outputs may carry `uncalibrated` flags. That is expected,
  and the README says so.

## Failure handling and resubmission

| Symptom (head job) | Class | Action |
|---|---|---|
| `TIMEOUT` (2 days), `NODE_FAIL`, `PREEMPTED`, or the head `CANCELLED` by the system | infrastructure | resubmit the **same command plus `--resume`** into the same outdir. Not counted toward R9's three attempts; recorded in `ATTEMPTS.tsv` |
| exit 2 (environment: nextflow, a container or a database missing) | environment | fix the environment (`build_envs.sbatch`, `giap setup`), resubmit with `--resume`; counts as an attempt |
| exit 1 at preflight (e.g. no target evidence at q ≤ 0.01, or IDs not in the gene model) | input | **not a GIAP bug to patch here**: a finding for M10-02/M10-03. The arm waits for corrected inputs and then re-runs fresh (`.r<n>`) |
| exit 4 naming a step (manifest `failed_step`, `failed_log`) | pipeline | read `failed_log`. A memory or time failure after Nextflow's retries means resubmitting with `--resume` (and a lead decision if the 54 GB cap is the cause, Q1). A code bug means a separate fix PR, the main checkout fast-forwarded, then the resume rule below. Counts as an attempt |
| exit 5 (`OUT/logs/internal_error.log`) | bug | as exit 4 with a code bug |
| exit 0 but a rule of set A fails | evidence | A3 or A5 (wrong commit or input): fresh re-run. A4 (`live_bakta`): fix the command, fresh re-run. A6 (`cpu_cap`): the run is valid evidence but breached R14. Record it, open a backlog item for the profile (index L1), and do not resubmit for that reason |

**Resume after a code fix (Q4).** `--resume` reuses cached tasks whose inputs and task hash did not
change. M7-15's per-step cache keys cover the config slice, not the Python code, so a cached
upstream task may be from older code. Rule: resume only when `git diff --name-only <old>..<new>`
touches only the failed step's module and modules downstream of it in the step registry. Otherwise
start a fresh `.r<n>` attempt. The rule applied is written into `ATTEMPTS.tsv` (`resume` column:
`fresh`, `resume-infra` or `resume-downstream`). The manifest of a resumed run records
`resume: true` and the new commit; rule A3 compares it with the expected (new) commit.

**R9.** Three counted attempts of one arm failing on the same cause: stop and file a Blocked row
(or a backlog entry, if it is a known bug with an owner). Then continue with M10-08/09 work.

**Compare failure.** Exit 1 with "different genomes" or "locus-tag sets differ" means the arms were
not run on the same genome and annotation, which the commands rule out. It is a bug finding, not
something to work around (no `--labels` or path changes). A rule C6 failure (non-determinism) is a
`giap compare` bug for M9-22's owner. In both cases the arms stay accepted and only the compare job
re-runs after the fix.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | S = 0 (no `giap-*` job), inputs and commit checked; A2 command above | `sbatch` prints a job ID J2; `ATTEMPTS.tsv` row `A2 1 out_a2 J2 fresh <sha> PENDING NA`; J2 goes into PROGRESS before any polling |
| 2 | S = 32 (M10-08's env build plus another agent's 16-CPU job) | 32 + 114 > 128: no submission; re-check later; nothing is lowered to fit |
| 3 | J2 ends `COMPLETED 0:0`; manifest `completed`, commit = expected, no `BAKTA` task, largest task AllocCPUS 14, `report.html` 2.1 MB, `modules_off: []` | `run_evidence.py arm … --arm A2` exits 0; `a2.json` `failures: []`; one markdown row printed for the README |
| 4 | J2 `TIMEOUT` after 2 days with 1,900 of ~2,650 tasks done | resubmit the same command plus `--resume` into `out_a2`; `ATTEMPTS.tsv` row 2 `resume-infra`; Nextflow skips the completed tasks; not an R9 attempt |
| 5 | A1 head exits 4; manifest `failed_step = OPERONS_STEP`; `failed_log` shows a `KeyError` on a plasmid TU | a backlog item or issue with the log path; fix PR merged; main fast-forwarded; the diff touches only `src/giap/operons/`, downstream of the failed step, so A1 resumes (`resume-downstream`); attempt 2 |
| 6 | A1 run with `--annotation …/nissle1917.json.gz` by mistake | exit 0 after a live Bakta task; rule A4 gives `live_bakta`, exit 1 from `run_evidence.py`; A1 re-runs fresh as `out_a1.r2` with the `.json` path |
| 7 | Both arms accepted; `pilot_compare.sbatch out_a1 out_a2 cmp_a1_a2` | `giap compare` exit 0; three tables, `plots/V16_compare.{svg,png,…}`, `report/compare.html`; `compare_abundance` one row `A1 A2 …`; matrix 2 × 4,982 rows (if `genes.tsv` lists 4,982 genes); the second call is byte-identical; `compare.json` `failures: []`; job exit 0 |
| 8 | `out_a1` accidentally produced from a genome FASTA with different headers but the same sequences | `giap compare` exit 1, "runs A1 and A2 are on different genomes …" (the digest covers the file); the job exits 1; A1 re-runs on the committed genome |
| 9 | `pilot_compare.sbatch` re-submitted with `cmp_a1_a2` already present | exit 2, "CMP exists"; nothing written |
| 10 | A2 preflight exits 1 "no target peptide at q ≤ 0.01" (a converter regression) | input class: finding for M10-03 with the preflight message; no GIAP change here; after the corrected tables, A2 re-runs fresh as `out_a2.r2` |

## Edge cases

- **Empty input.** A peptide table with no target row at q ≤ 0.01 cannot have passed M10-03's
  summary rules. If it reaches `giap run`, preflight exits 1 (no evidence). This is the input class
  above, never a GIAP patch. An arm with targets but zero detected genes would pass rule set A (the
  pipeline ran) and be reported by M10-05. `giap compare` then gives a matrix of `not_detected` and
  a ρ row with NA and its flag. C1 still holds if the schema allows NA, so the test-author pins this
  case on a fixture.
- **A1's cross-strain peptides.** Many A1 peptides are K-12-only: unmapped or tolerant, and some
  tolerant hits become `strain_difference` variants. These are expected outputs, not failures.
  Genome-inconsistent counts are metric 2 (M10-05).
- **Circular replicons and the origin.** All four Nissle replicons are circular-tagged. Genes and
  novel ORFs that span the origin are handled by the pipeline. `giap compare` joins on locus tags,
  so no coordinate is read here. The compare matrix covers every replicon, including the three
  small plasmids.
- **Plasmid count.** §15.1 names "pMUT1 and pMUT2"; the assembly GCA_043228065.1 has three
  plasmid records (4.5, 3.2 and 5.5 kb). The runs use the assembly as committed. The README notes
  the difference, and M10-05 reports plasmid genes per replicon, not per name.
- **Minus strand.** Not touched by this issue: mapping and annotation are GIAP's. A1 and A2 use the
  same gene model, so a strand error would show in both and in M10-05's TE1, not in compare.
- **Ties.** `giap compare` handles tied abundances (Spearman with average ranks, as M9-22). Two
  attempts of one arm both passing: the **latest** passing attempt is the accepted one, and the
  earlier is kept but marked `superseded` in `ATTEMPTS.tsv`. Labels cannot tie: `A1,A2` are fixed.
- **Detection tier differs between arms** (e.g. A2 `decoy_picked`, A1 `mayu`): legal. Compare
  shows each run's tier in `compare.html`. Recorded in the evidence (`detection_tier`), not failed.
- **Duplicate peptides across runs.** The tables hold three raw runs (`run` = `WT_1..3`).
  Aggregation is preflight's (§5.3); no `--replicate` handling is needed here (ADR-0042 applies to
  compare's replicate counts, and both arms carry the same three runs).
- **Head job killed before the tail** (no `run_manifest.json`): rule A1 `no_manifest`. Resume
  writes the manifest on completion.
- **`report.html` over its size budget.** `giap run` does not fail on it, so it is recorded (bytes),
  not gated. `giap report --strict` is not run here.
- **Session ends while jobs run.** Job IDs are in PROGRESS and `ATTEMPTS.tsv` before polling. The
  evidence is computed after the fact from the outdir and `sacct`, so nothing depends on the
  session.
- **Other agents' jobs start after submission** and push the sum above 128. R14 is checked at
  submission. A breach caused by later submissions belongs to them. This arm's own bound is 114
  CPUs.
- **Disk.** Each arm's `OUT/work` may reach tens of GB. Before each submission, `du -sh` of
  `runs/pilot/` goes into PROGRESS (§3.9). M10-14 runs `nextflow clean` on the pass-1 work dirs.

## Acceptance criteria

On CI (every commit; no Slurm, no real data, no `$GIAP_ROOT`):

1. `arm_evidence` on a hand-made run directory has `failures: []` and the listed fields. The
   directory holds a valid minimal `run_manifest.json`, a `trace.tsv` with two tasks and
   `native_id`s, `report/report.html`, a `tables/genes.tsv` that passes `check_tsv`, plus
   injected `sacct` rows.
2. Each rule A1–A9 fails on a fixture built to break it, with the exact failure string, and the CLI
   exits 1. Cases: a missing manifest; `status: failed` with `failed_step`; another commit; `git_dirty`
   true; a `BAKTA` trace row; a wrong peptides sha256; a task with AllocCPUS 16; no `report.html`; a
   `genes.tsv` missing a required column; `modules_off: ["function"]`.
3. A run with empty results (a `genes.tsv` with all genes `not_detected`, no novel ORF table rows)
   still passes rule set A.
4. `compare_evidence` on a hand-made `CMP` passes, and each of C1–C6 fails on its fixture: a column
   missing; two abundance rows; a matrix row short; a family missing; no `compare.html`; one byte
   changed in `CMP2`'s detection table. The fixtures are better made by running the real `giap
   compare` on two tiny synthetic run directories, as M9-22's tests do (no hand-copied real values).
5. `pilot_compare.sbatch` in a sandbox (`GIAP_PYTHON` and `PYTHONPATH` set, the two tiny run
   directories): exit 0 and `compare.json`; exit 2 on a wrong argument count, a missing run
   directory or an existing `CMP`; exit 1 when the two runs have different genome digests. Nothing
   written under the sandbox repo.
6. House rules for the new sbatch: a `giap-*` name; `--partition=normal`, `--qos=normal`,
   `--cpus-per-task=2`, `--mem=8G`, `--time=01:00:00`; `--output` under `local/logs/`; the
   environment sourced; `set -euo pipefail`; no `uv run`, `uv sync` or `pip`.
7. Guard (extends M10-03's): no tracked file under `validation/pilot/` matches `out_a*`, `cmp_*`,
   `evidence/*.json`, `ATTEMPTS.tsv`, `report.html` or `compare.html`. The README section contains no
   number outside the provenance fields: a test reads the section and allows only job IDs, sha256s,
   byte sizes, durations, memory, exit codes, dates, CPU counts and the seed (Q-K1).

On FarmShare (evidence; PR description, PROGRESS, README; not CI):

8. Both arms' accepted attempts: `run_evidence.py arm` exit 0. Head job IDs, `sacct` State
   `COMPLETED`, ExitCode `0:0`, Elapsed, MaxRSS and NodeList are in PROGRESS › Measured numbers, and
   the top five steps by wall time come from `measure_runtimes.py`.
9. `pilot_compare.sbatch` exit 0; `compare.json` `failures: []`.
10. The README section lists, per arm: the command, every attempt's job ID and class, the accepted
    outdir, `git_commit`, input sha256s (genome, annotation `.json`, peptides), and the sha256 and
    bytes of `report.html`. For compare: the job ID and `compare.html` sha256.
11. `ATTEMPTS.tsv` has one row per submission; no arm has more than three counted attempts.

Link to §15.6: criteria 8–9 are the evidence for bullet 1 (A1 and A2 parts) and bullet 3. M10-05
reads `out_a1`, `out_a2` and `cmp_a1_a2` from the accepted paths in `evidence/*.json`.

## Test plan

Owner: the test-author, tests first (R2), new files only:

- `tests/unit/test_m10_04_run_evidence.py`: criteria 1–4. Fixtures built in `tmp_path`. The
  manifest comes from a minimal valid document against `run_manifest.schema.json`; `sacct` rows are
  passed as data, never from a `sacct` call. The compare fixtures come from running `giap compare`
  on two tiny run directories made with the existing test helpers (as `tests/unit/compare/` and
  `tests/unit/test_cli_compare.py` do).
- `tests/unit/test_m10_04_compare_job.py`: criteria 5 and 6, in a sandbox like
  `helpers.pilot_kit.pilot_sandbox`, with the real `giap compare`.
- `tests/unit/test_m10_04_guard.py`: criterion 7.

Local runs: targeted only, with the venv python and `PYTHONPATH=$PWD/src:$PWD/tests/sim`; never
`uv run` in a worktree. The full suite goes to Slurm (`test_branch.sbatch`).

## Order of work

1. Spec (this), then the test-author's tests and the implementer's script and sbatch. These are
   independent of the data and can merge before M10-03's jobs run.
2. **MVP-04 acceptance (soft dependency).** If its K-12 run is pending or running when M10-04 is
   ready, the pilot does not run alongside it: its 114 CPUs plus MVP-04's do not fit under 128. The
   pilot waits, and gains MVP-04's lessons (fixes merged to main). If MVP-04's acceptance has not
   been submitted when M10-03's tables exist, the pilot does **not** wait for it. The pilot run is
   also an end-to-end real-genome run, and its bugs are equally useful (Q7).
3. A2, then A1, then compare; evidence; README; PROGRESS.

## Estimate

| Piece | Code lines |
|---|---|
| `validation/pilot/run_evidence.py` (`arm_evidence`, `compare_evidence`, CLI) | 60 |
| `scripts/slurm/pilot_compare.sbatch` | 25 |
| README section | prose |
| **Total** | **≈ 85** (index: 80) |

Tests ≈ 200 lines. One PR. Wall time is dominated by the runs, about 3–8 h per arm (M7-15: 6 h on
K-12), so ≈ 6–16 h for both in sequence plus queue waits. This is not the index's 35–60 min
(§13.2). M10-12 measures the standard run.

## Questions

### For Khoa (parked in `local/handoff/khoa_q.md`; the defaults do not block)

| # | Question | Default (reversible) |
|---|---|---|
| Q-K1 | Your rule keeps everything derived from the DIA-NN reports out of the repository. Do **aggregate numbers** count as derived: genes detected per arm, ρ between A1 and A2, Jaccard of islands? M10-04 commits none; M10-05's `REPORT.md` and `validation/pilot/metrics/` (index decision 5, §15.5) need them | M10-04: provenance only (IDs, exit codes, hashes, sizes, runtimes). M10-05 asks again with a concrete list; until then, aggregates are treated as allowed **in `REPORT.md` only**, with no row-level table committed |

### For the lead

| # | Question | Default (reversible) |
|---|---|---|
| Q1 | The MVP-04 default was `--threads 14`. Add `--memory 54.GB`, so that no task exceeds 14 CPUs under FarmShare's `MaxMemPerCPU = 4000 MB` (one arm ≤ 114 CPUs)? And if a task fails at the 54 GB cap? | yes, `--threads 14 --memory 54.GB`. A failure at the cap gets a lead decision: re-run that arm with `--memory 64.GB` only while no other `giap-*` job runs (≤ 130 CPUs worst case), recorded as an R14 exception; the profile fix (index L1) moves up the backlog |
| Q2 | Order and concurrency: A2 first, then A1, never both at once? (Two arms at `--threads 7 --memory 27.GB` would fit 116 CPUs together, but DETERMINANTS asks 48 GB) | sequential, A2 first |
| Q3 | `--annotation *.json.gz` silently runs Bakta live (`cli._bakta_json` accepts only `.json`). Decompress for the pilot and fix the CLI in a separate small PR, or fix it first? | decompress now (`$GIAP_ROOT/pilot/inputs/nissle1917.bakta.json`); backlog chore "accept `.json.gz` as a Bakta JSON in `giap run`" with a test; the pilot does not wait for it. Rule A4 guards either way |
| Q4 | Resume after a code fix: allowed only when the diff touches the failed step and downstream steps, else a fresh attempt? | yes, as written; the class goes in `ATTEMPTS.tsv` |
| Q5 | The two arms repeat the same annotation-layer work (TerminatorNet, Promoter Calculator, geNOMAD, antiSMASH: most of the hours). Feed A1 the layers of A2 as precomputed layers (MVP-02) to halve the wall time? | no for pass 1: two independent runs are the cleaner evidence, and the freeze runs (M10-14) repeat them anyway. Backlog item for M10-12/13 |
| Q6 | C4: how does `compare_stability` write a family with no item in either run (a row with `n_a = n_b = 0` and `jaccard` NA plus a flag, or no row)? | the test-author reads M9-22's behaviour and pins it. If it writes no row, C4 becomes "every family present in either run has a row" |
| Q7 | MVP-04's acceptance run is a soft dependency: wait for it, or run the pilot first if M10-03 finishes first? | do not wait; never run both at once (R14) |
| Q8 | "Done" has no biological or report-content gate (those are M10-05 and the M9 PRs), and pipeline bugs are fixed in separate PRs. Agreed? | yes |
| Q9 | Outdirs under `$GIAP_ROOT/runs/pilot/pass1/` with §15.2's names (`out_a1`, `out_a2`, `cmp_a1_a2`), `.r<n>` for fresh re-runs and job IDs in `ATTEMPTS.tsv`, rather than a job ID in the name? | yes |
