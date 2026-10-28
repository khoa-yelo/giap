# M2-15b: shard the TerminatorNet scan (split → chunk → merge)

Milestone M2 (follow-up; Khoa approved 2026-10-27, performance pass) · depends on: M2-15, ADR-0062 ·
estimate ≈ 150 code lines (script + 3 processes + subworkflow) and ≈ 200 test lines

## Why
`ANNOTATE:TERMINATORNET` sits on the critical path of every real run. A 5-Mb genome takes 50 min at
16 CPUs (ADR-0062), and on slow nodes it takes 107–175 min at 8 CPUs. msmeg failed with exit 140
after 12 h (PROGRESS.md, #411). Three things make it slow:
- everything runs in one task, capped at 16 CPUs (`resourceLimits`);
- the + strand and the − strand run one after the other;
- `predict_proba` runs serially in the parent, one row at a time, after the two pool phases.

A 12-h timeout throws away the whole scan, and a 16-CPU request waits longer in the queue than a
small one.

## Upstream facts this design relies on (pinned commit 72180309, read 2026-10-27)
Source: `$CONDA_PREFIX/share/terminatornet/{TerminatorNet.py,Features.py}` in `gippa-terminatornet`.

| Fact | Where | Consequence |
|---|---|---|
| `HAIRPIN_LENGTH, TAIL_LENGTH = 52, 8`; `WINDOW = 60`, `STEP = 1`; windows `j in range(0, max(L-WINDOW+1, 1))` | Features.py 313–314, 372 | `nwin(L) = max(L-59, 1)` windows per strand sequence |
| A window reads only its own neighbourhood of the sequence passed in. `scan_oneprocess_1(j)` and `scan_oneprocess_2` read the module globals `Features.SEQUENCE` and `Features.GC_PERCENT` | Features.py | a chunk that sets `SEQUENCE` to the **whole** strand sequence and maps only its own *j* computes exactly the values the unsharded run computes for those *j* |
| `GC = 0.5` is forced after `clean_sequences` | TerminatorNet.py `command_line_arguments` | no genome-wide statistic crosses a chunk boundary |
| `clean_sequences` upper-cases, maps U→T, **drops non-ACGT characters** (when ≥ 75 % are ACGT) and blanks sequences shorter than 4 nt; `CONTAINS_LONG_SEQUENCE` = some cleaned sequence is ≥ 501 nt | TerminatorNet.py | coordinates are in the *cleaned* sequence, so every task must clean exactly as upstream does |
| `read_in_sequences` drops records with no sequence lines | TerminatorNet.py | an empty record produces no output line |
| `-n ≥ 2` with a long sequence: `scan_multiprocess` (pool phase 1 = TT score ≥ 30, pool phase 2 = features, then serial per-row `scaler.transform` + `predict_proba`, keep `pred ≥ 0.5`, then **stable sort by stop**). With no long sequence, or with `-n 1`: `scan()`, which keeps window order and does not sort | Features.py, `run_TerminatorNet` | the order of the calls entering `merge_overlapping` depends on the mode |
| `merge_overlapping` is a greedy sweep in list order. A later call replaces the best one when `pred >= best_pred`. With `-s` the − strand is mapped back with `seq_length` and its block is sorted by forward stop | TerminatorNet.py | merging must run **once per (contig, strand)**, over all of its chunks, in the upstream order |
| Output: every contig's + block in FASTA order, then every contig's − block. `TERM_k` restarts in each block. An empty block is written as `name\tNO TERMINATORS FOUND` | `run_TerminatorNet`, `output_terms` | the merge writes blocks in this order |

Probe (2026-10-27, 4-record FASTA, `-s -n 2`): the empty record wrote nothing. The 2-nt record wrote
`NO TERMINATORS FOUND` on both strands. A record with Ns and lower-case letters was reported in
N-stripped coordinates.

## Reference output (lead decision 2)
The reference is the upstream output at **`-n ≥ 2`**:
- if any cleaned sequence is ≥ 501 nt, the calls are stop-sorted before merging;
- otherwise they stay in window order.

That is what the deployed module produces (`-n ${task.cpus}`, 16 on FarmShare). The sharded step
reproduces it at any CPU count, including 1. Unsharded `-n 1` can differ from it, and is not a
reference.

## What changes
1. **Script.** New file `workflow/bin/gippa_terminatornet_shard.py` (on the PATH of Nextflow tasks).
   It runs in the `gippa-terminatornet` env and is not part of `src/gippa`. It imports `TerminatorNet`
   and `Features` from `$TERMINATORNET_HOME` (MIT; fetched, not vendored, R5) only inside the
   subcommands. The planning and ordering functions are pure and import nothing from upstream, so the
   fast suite can test them by loading the file by path. The constant `CORE = 250_000` windows is not a
   config key; it is a parameter of the pure functions so tests can use small values.
2. **SPLIT (`plan`).** Reads `genome.fa` with upstream `read_in_sequences` and `clean_sequences`. For
   each record *i* (FASTA order) and each strand (+ = cleaned sequence S, − = upstream
   `reverse_complement(S)`, which has the same length L), core *k* = windows
   `[k·CORE, min((k+1)·CORE, nwin(L)))`, for k = 0, 1, ….

   `plan.tsv` has a header and one row per chunk:
   `chunk_id seq_index seq_name strand seq_length core_start core_end mode`.
   - `chunk_id` = `c{seq_index:05d}{p|m}{k:04d}`, e.g. `c00000p0000`;
   - `mode` = `long` if upstream `CONTAINS_LONG_SEQUENCE`, else `short` (the same on every row);
   - rows are ordered by (seq_index, strand + before −, k).

   If no record has any sequence, the step exits 1 with upstream's "could not read in genomic
   sequences" message, as the unsharded run does.
3. **CHUNK (`chunk <chunk_id>`), whole-strand scanning (lead decision 3).** Reads and cleans the genome
   again exactly as in `plan` (it is deterministic and takes seconds), then takes the plan row's strand
   sequence (S, or the reverse complement of S).
   - It sets `Features.SEQUENCE` to that **whole** strand sequence and `Features.GC_PERCENT = 0.5`.
   - It maps upstream `scan_oneprocess_1` over the **owned global window indices** `range(core_start,
     core_end)` with `multiprocessing.Pool(task.cpus)`.
   - It keeps `j ≥ 0`, maps `scan_oneprocess_2`, and runs the same serial per-row
     `scaler.transform(np.array(x).reshape(1, -1))` and `model.predict_proba(x)[0, 1]`, keeping
     `pred ≥ 0.5`.

   There is no padding, no slicing and no coordinate shift. Each window sees exactly the sequence it
   sees in the unsharded run (including the contig-end clamping and the shortened hairpin when
   L < 60), so the per-window values are identical by construction.

   Ownership is decided by the window start *j*. The call tuples do not carry *j*, so the filter has
   to be the index list itself, and every window is evaluated by exactly one chunk.

   Rows are **not** batched: batching could change the last bits of `pred`, and the TSV prints `repr`
   (backlog, item 10).

   Output: `calls_<chunk_id>.pkl` (pickle protocol 4) holding `chunk_id`, `core_start`, `core_end`
   and the call tuples in window order, with upstream types unchanged (`np.float64` `pred`). The model
   is loaded from `model.pkl` in the cwd (symlinked as in the unsharded module). The chunk sets no BLAS
   thread variable, so its environment matches the unsharded module.
4. **MERGE (`merge`).** For each (seq_index, strand) in upstream output order (every + block in FASTA
   order, then every − block):
   - concatenate its chunks' calls in k order;
   - if `mode = long`, stable-sort them by stop (itemgetter(1)), as `scan_multiprocess` does;
   - set `TerminatorNet.BOTH_STRANDS = True` and call upstream `merge_overlapping(calls, L, strand)`
     (the + strand uses upstream's default arguments);
   - write the block with upstream `output_terms` to `terminatornet.tsv`.

   Within a chunk the calls are in window order, and chunks are concatenated in core order, so the
   concatenation is in global window order. A stable sort of it therefore equals upstream's sort.
   Before writing anything the merge checks that every plan `chunk_id` has exactly one calls file
   whose `core_start`/`core_end` match the plan. A missing, extra or mismatched file exits 1. The
   merge emits `versions.yml` with the same `terminatornet: <GIPPA_COMMIT>` entry as before.
5. **Params (lead decisions 1 and 4).** Two new Nextflow params in `workflow/nextflow.config`. They are
   not defaults.yaml keys.
   - `terminatornet_sharded = true`. With `false`, annotate.nf runs the old unsharded
     `modules/local/terminatornet` module, unchanged, and that module's resources come back with it.
     The old module also stays as the acceptance reference and as `workflow/tests/real/terminatornet.nf`.
   - `terminatornet_chunk_cpus = 4`. This is `TERMINATORNET_CHUNK`'s `cpus`, clamped by
     `resourceLimits`. It can be changed per run (`--terminatornet_chunk_cpus N`), and msmeg is
     measured at 4.
6. **Nextflow.** New modules `workflow/modules/local/terminatornet_{plan,chunk,merge}/main.nf` (conda
   `envs/terminatornet.yml`, the same container) and a subworkflow
   `workflow/subworkflows/local/terminatornet_sharded.nf`:

   | Process | Input | Output | Resources (base.config `withName`) |
   |---|---|---|---|
   | `TERMINATORNET_PLAN` | `tuple val(meta_id), path(genome, stageAs: 'genome.fa')` | `path 'plan.tsv'` | `process_single` (1 CPU, 2 GB) |
   | `TERMINATORNET_CHUNK` | `tuple val(chunk_id), path(genome, stageAs: 'genome.fa'), path(plan)`; tag `chunk_id` | `path "calls_${chunk_id}.pkl"` | `cpus = { params.terminatornet_chunk_cpus }`; `memory = { 4.GB * task.attempt }`; FarmShare: `time = { 2.h * task.attempt }` |
   | `TERMINATORNET_MERGE` | `tuple val(meta_id), path(genome, stageAs: 'genome.fa'), path(plan)`, `path(calls)` (collected) | `path 'terminatornet.tsv', emit: tsv`; `path 'versions.yml'` | `process_single` |

   - The chunk channel is `plan.splitCsv(header: true, sep: '\t').map { it.chunk_id }`, combined with
     the genome and plan.
   - The merge gets `CHUNK.out.calls.collect()`. It reads the files in plan order, so the channel
     order does not matter.
   - annotate.nf chooses `TERMINATORNET_SHARDED(gin)` or `TERMINATORNET(gin)` from
     `params.terminatornet_sharded` and keeps `kv('terminatornet', … .tsv)`, so `TERMINATOR_LAYER` and
     its inputs are unchanged.
   - The container, conda and prebuilt-env selectors are widened to `TERMINATORNET.*`.
   - FarmShare's `max(cpus, ceil(memory / 4000 MB))` rule is repeated for the CHUNK selector. At 4 CPUs
     and 4 GB a chunk asks for 4 CPUs; a retry at 8 GB still asks for 4.
   - The run-review process→step mapping maps all three processes to `terminatornet`.
   - Stubs: PLAN writes `params.stub.terminatornet_chunks ?: 1` rows (the key is added to
     `params.stub`), CHUNK touches its pickle, and MERGE writes `genome\tNO TERMINATORS FOUND` as the
     old stub does.
7. **Retry and resume.** The 137/140/143 policy applies per chunk, so a timeout re-runs one chunk
   (about 5 min at 4 CPUs) instead of the whole scan. Chunk hashes depend only on the genome, the plan,
   `chunk_id` and the CPU count, so `-resume` keeps finished chunks.
8. **Topology.** TerminatorNet treats every sequence as linear. Circular replicons are not wrapped:
   there are no windows across the origin, exactly as in the unsharded run.
9. **R14 impact.**
   - Chunk count = Σ over contigs of 2·⌈nwin/CORE⌉: about 40 for 5 Mb, about 56 for msmeg (≈ 7 Mb),
     and 4 for sim100k (chromosome + plasmid, two strands each).
   - Each chunk takes about 76 CPU-min/Mb-strand × 0.25 Mb ≈ 19 CPU-min, so about 5 min at 4 CPUs
     (from ADR-0062's 50.3 min × 15.1 CPUs ÷ 10 Mb-strands).
   - With FarmShare `queueSize = 8`, the step uses at most 8 × 4 = **32 CPUs** at a time, up from 16.
     The run's bound of 8 × 16 = 128 CPUs is unchanged.
   - The chunks can fill all 8 slots, so other ready tasks wait behind them. Lowering
     `terminatornet_chunk_cpus` does not free slots; only fewer chunks would.
   - Submitting 40–56 jobs at `submitRateLimit = 10/1min` takes 4–6 min, overlapped with execution.
   - Memory falls from 16 GB in one task to 4 GB per chunk.
10. **Out of scope (backlog).** Batched `predict_proba`, a config key for CORE, scanning across the
    origin of circular replicons.

## ADR
No new ADR (lead decision 4). The implementer adds an **amendment to ADR-0062**. It records:
- that `TERMINATORNET` is replaced by `TERMINATORNET_PLAN`/`_CHUNK`/`_MERGE` behind
  `params.terminatornet_sharded = true`, with the unsharded module as the fallback;
- the chunk resources (`params.terminatornet_chunk_cpus = 4`, 4 GB × attempt, 2 h × attempt on
  FarmShare) and the 32-CPU concurrency cap under `queueSize = 8`;
- the correction that "output does not depend on `-n`" holds only for `-n ≥ 2` (`-n 1` skips the stop
  sort), and that the `-n ≥ 2` output is the reference;
- the msmeg measurement once it exists.

## Examples (CORE = 250,000 unless stated)
| # | Input | Expected |
|---|---|---|
| 1 | one linear contig, L = 1,000,000, any mode | nwin = 999,941. + cores `[0,250000) [250000,500000) [500000,750000) [750000,999941)`. The same 4 cores on − in RC coordinates. 8 plan rows `c00000p0000…c00000m0003` |
| 2 | plasmid shorter than one chunk, L = 120,000 | one core per strand, `[0,119941)`; 2 rows |
| 3 | contig shorter than a window, L = 40 | nwin = 1: core `[0,1)`. The single upstream window uses the shortened hairpin `s[0:32]` |
| 4 | empty record (`>e` with no sequence) between two contigs | no plan row and no output line for `e`. `seq_index` counts only the records upstream keeps |
| 5 | 2-nt record `>t AC` | cleaned to `''`: core `[0,1)`. Output `t\tNO TERMINATORS FOUND` in the + section and again in the − section |
| 6 | boundary: + strand, CORE = 1,000, a hairpin from window j = 995 at `[1003, 1030]` (U-tract to 1,038) | evaluated only by chunk 0 (owns j < 1,000), reading the whole sequence. Chunk 1 does not scan j = 995. The call appears once, exactly as in the single-chunk run |
| 7 | tie across a chunk edge: call A (window in core 0; start 980, stop 1,010, pred 0.93) and call B (core 1; start 1,000, stop 1,030, pred 0.93) | one merged call, B (`pred >= best_pred`: the later one wins), the same as the unsharded run. With B.pred = 0.92 the result is A |
| 8 | lower-case letters and Ns: `>c ACGTNNNNacgt…` | plan coordinates are in the cleaned (N-stripped, upper-case) sequence. Output equals the unsharded output (probe values above) |
| 9 | the same chunk run at `task.cpus` 1, 4 and 8 | byte-identical `calls_*.pkl` |
| 10 | genome whose records are all < 501 nt (`mode = short`) | no stop sort, calls in window order. Output equals the unsharded `-n 2` run |

## Tests (test-author)
- **Fast (pure, file loaded by path).**
  - The plan's cores **partition `[0, nwin)` exactly** for every (contig, strand): sorted by k, the
    first starts at 0, each starts where the previous one ends, the last ends at nwin, and none is
    empty. This holds over lengths 0–3,000 and CORE ∈ {1, 7, 61, 1000}, plus examples 1–5.
  - Chunk ids are unique and in sort order.
  - Merge order (+ blocks, then − blocks) and the mode-dependent stop sort, checked with fake call
    tuples (example 7 and its 0.92 variant).
  - Merge validation: a missing, extra or mismatched pickle exits 1.
- **`needs_conda` (Slurm).** Examples 3, 5, 6, 7, 8, 9 and 10 with the real model on small synthetic
  FASTAs: sharded output, compared with `cmp`, against `TerminatorNet.py -s -n 2`. Then sim100k at
  CORE = 10,000: byte-identical to `tests/data/recorded/terminators/terminatornet/sim100k/terminatornet.tsv`.
- **nf-test (stub).**
  - 3 stub chunks give 3 CHUNK tasks, and MERGE emits `terminatornet.tsv`.
  - `--terminatornet_sharded false` runs the single `TERMINATORNET` task.
  - `--terminatornet_chunk_cpus 2` gives CHUNK 2 CPUs.

  The lead updates the existing tests that name `TERMINATORNET`: `annotate.nf.test`,
  `main_genetable_context.nf.test`, `lib/annotate_stub.groovy`, `test_tool_resources_pins.py`,
  `test_nf_conda_profile.py`, `test_m11_06_*`, `test_m11_14_images_json.py` and the runner trace mapping.

## Acceptance
- Examples 1–10 hold, and the fast suite, nf-test, ruff and mypy pass.
- `terminatornet.tsv` from the sharded pipeline (`-profile farmshare,apptainer`, CORE = 250,000,
  4-CPU chunks) is byte-identical (`cmp`) to these unsharded references:

  | Sample | Reference | sha256 |
  |---|---|---|
  | sim100k | `$GIPPA_ROOT/recordings/terminatornet/sim100k/terminatornet.tsv` (= the committed fixture) | `d0898b10…` |
  | Nissle 1917 | `$GIPPA_ROOT/recordings/terminatornet/nissle1917/terminatornet.tsv` | `c319d52f…` |
  | K-12 MG1655 | `$GIPPA_ROOT/recordings/terminatornet/k12_mg1655/terminatornet.tsv` | `67307825…` |
  | mhyo, gsulf | the fresh unsharded runs (prerequisite below) | recorded in the PR |

  Each recording's `genome.fa` is the input. Checked on 2026-10-27, the recordings are at the current
  pin:
  - `versions.txt` gives commit `72180309…`, the same as the env's `GIPPA_COMMIT`, recorded with
    `-n 17`;
  - the explicit lock sha256 `d437c792…` equals the current `workflow/envs/locks/terminatornet.explicit.txt`;
  - the recorded pip lock (`ab22b7d9…`) differs from the current one only by packages that #383
    removed and TerminatorNet does not import (biopython, plotnine, …). numpy, scipy and
    scikit-learn are unchanged.
- **Validation prerequisite: fresh unsharded mhyo and gsulf reference runs** (the old module,
  `-n ≥ 2`, the same pin), copied to `$GIPPA_ROOT/recordings/terminatornet/{mhyo,gsulf}/` together
  with their `genome.fa`. There is no pinned reference under `$GIPPA_ROOT/runs` or `validation`. The
  #411 outputs exist only as Nextflow work files:
  - `local/more_datasets/mhyo/out/work/e9/564ba3…` (`-n 3`, sha256 `17fd2323…`);
  - `local/more_datasets/gsulf/out/work/dc/af3e1d…` (`-n 6`, sha256 `cfa108c9…`).

  Both are at pin 72180309, but their staged `genome.fa` is gone. They serve as a cross-check only.
- On **msmeg** (timed only), PROGRESS.md records the following at `terminatornet_chunk_cpus = 4`, next
  to the failed unsharded attempt (exit 140 at 12 h):
  - the wall time of the sharded step (first PLAN submit to MERGE end);
  - the chunk count and the max chunk wall time;
  - max RSS and peak concurrent CPUs (from the trace).
- No run exceeds the R14 caps (≤ 128 CPUs, ≤ 20 jobs per user across GIPPA jobs).
- The implementing PR adds `changelog.d/<PR>.perf.md` (or the matching type).

## Lead decisions (2026-10-27, on #451)
1. CPUs per chunk come from `params.terminatornet_chunk_cpus`, default 4 (cap of 32 CPUs under
   `queueSize = 8`). msmeg is measured at 4, and the value can be changed.
2. The reference is the `-n ≥ 2` (stop-sorted) output; this goes in the ADR-0062 amendment.
3. Chunks scan their own windows on the whole cleaned strand sequence: no overlap and no coordinate
   shift. The partition test replaces the overlap property test.
4. Sharding is always on, with the fallback `params.terminatornet_sharded = true` (false = the old
   module). No new ADR; an ADR-0062 amendment instead.
5. CORE = 250,000 windows.
6. References: `$GIPPA_ROOT/recordings/terminatornet/{sim100k,nissle1917,k12_mg1655}` (verified at
   the pin); fresh unsharded mhyo and gsulf runs as a prerequisite; msmeg timed only.

## Open questions for the lead
None blocking. One non-blocking note: whether 56 four-CPU msmeg chunks holding all 8 queue slots
delays other ready tasks enough to matter shows up in the msmeg trace. If it does, raise CORE rather
than lowering the CPUs.
