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
2-CPU one.

## Upstream facts this design relies on (pinned commit 72180309, read 2026-10-27)
Source: `$CONDA_PREFIX/share/terminatornet/{TerminatorNet.py,Features.py}` in `gippa-terminatornet`.

| Fact | Where | Consequence |
|---|---|---|
| `HAIRPIN_LENGTH, TAIL_LENGTH = 52, 8`; `WINDOW = 60`, `STEP = 1`; windows `j in range(0, max(L-WINDOW+1, 1))` | Features.py 313–314, 372 | `nwin(L) = max(L-59, 1)` windows per strand sequence |
| Window *j* reads `s[j : j+52]` (hairpin, `min(j+52, L-8)`), `tail_TT` reads `s[j+52 : j+67]` (`TT_TAIL_LENGTH = 15`; positions past the end count as "not T"), U-tract `s[j+he+1 : j+he+9]` with `he ≤ 51`, A-tract `s[max(j+hs-8, 0) : j+hs]` with `hs ≥ 0` | Features.py 160–166, 320–324 | window *j* depends only on `s[j-8 : j+67)`, i.e. **8 nt to the left and 67 nt to the right** (75 nt in total) |
| `GC = 0.5` is forced after `clean_sequences` | TerminatorNet.py `command_line_arguments` | no genome-wide statistic crosses a chunk boundary |
| `clean_sequences` upper-cases, maps U→T, **drops non-ACGT characters** (when ≥ 75 % are ACGT) and blanks sequences shorter than 4 nt; `CONTAINS_LONG_SEQUENCE` = some cleaned sequence is ≥ 501 nt | TerminatorNet.py | coordinates are in the *cleaned* sequence, so planning must clean exactly as upstream does |
| `read_in_sequences` drops records with no sequence lines | TerminatorNet.py | an empty record produces no output line |
| `-n ≥ 2` with a long sequence: `scan_multiprocess` (pool phase 1 = TT score ≥ 30, pool phase 2 = features, then serial per-row `scaler.transform` + `predict_proba`, keep `pred ≥ 0.5`, then **stable sort by stop**). With no long sequence, or with `-n 1`: `scan()`, which keeps window order and does not sort | Features.py, `run_TerminatorNet` | the order of the calls entering `merge_overlapping` depends on the mode |
| `merge_overlapping` is a greedy sweep in list order. A later call replaces the best one when `pred >= best_pred`. With `-s` the − strand is mapped back with `seq_length` and its block is sorted by forward stop | TerminatorNet.py | merging must run **once per (contig, strand)**, over all of its chunks, in the upstream order |
| Output: every contig's + block in FASTA order, then every contig's − block. `TERM_k` restarts in each block. An empty block is written as `name\tNO TERMINATORS FOUND` | `run_TerminatorNet`, `output_terms` | the merge writes blocks in this order |

Probe (2026-10-27, 4-record FASTA, `-s -n 2`): the empty record wrote nothing. The 2-nt record wrote
`NO TERMINATORS FOUND` on both strands. A record with Ns and lower-case letters was reported in
N-stripped coordinates.

## What changes
1. **Script.** New file `workflow/bin/gippa_terminatornet_shard.py` (on the PATH of Nextflow tasks).
   It runs in the `gippa-terminatornet` env and is not part of `src/gippa`. It imports `TerminatorNet`
   and `Features` from `$TERMINATORNET_HOME` (MIT; fetched, not vendored, R5) only inside the
   subcommands. The planning and ordering functions are pure and import nothing from upstream, so the
   fast suite can test them by loading the file by path. Constants: `CORE = 250_000` windows and
   `PAD = 200` nt. Neither is a config key (defaults.yaml is protected); both are parameters of the pure
   functions so tests can use small values.
2. **SPLIT (`plan`).** Reads `genome.fa` with upstream `read_in_sequences` and `clean_sequences`. For
   each record *i* (FASTA order) and each strand (+ = cleaned sequence S, − = upstream
   `reverse_complement(S)`, which has the same length L):
   - core *k* = windows `[k·CORE, min((k+1)·CORE, nwin(L)))`, for k = 0, 1, …;
   - chunk slice = `[max(0, core_start − PAD), min(L, core_end + PAD))` in that strand's own coordinates.

   `plan.tsv` has a header and one row per chunk:
   `chunk_id seq_index seq_name strand seq_length core_start core_end chunk_start chunk_end mode`.
   - `chunk_id` = `c{seq_index:05d}{p|m}{k:04d}`, e.g. `c00000p0000`;
   - `mode` = `long` if upstream `CONTAINS_LONG_SEQUENCE`, else `short` (the same on every row);
   - rows are ordered by (seq_index, strand + before −, k).

   If no record has any sequence, the step exits 1 with upstream's "could not read in genomic
   sequences" message, as the unsharded run does.
3. **Overlap size.** Window *j* reads `[j−8, j+67)`. Owned windows satisfy `j < core_end`, and the slice
   reaches `core_end + PAD`. So `PAD ≥ 67` is enough on the right and `PAD ≥ 8` on the left. That also
   keeps `min(j+52, len−8) = j+52` for every owned *j* away from the contig end. `PAD = 200` gives
   about 3× margin. It costs nothing, because padding windows are never scanned (item 4). At contig
   edges the slice is clamped to `[0, L)`, which reproduces upstream's own clamping (A-tract `max(…, 0)`,
   tail past the end = "not T", the shortened hairpin when L < 60). At start-up the script asserts
   `PAD ≥ max(TAIL_LENGTH, HAIRPIN_LENGTH + TT_TAIL_LENGTH)` against the imported upstream constants,
   so a re-pin that widens the reach fails loudly instead of drifting.
4. **CHUNK (`chunk <chunk_id>`).** Reads and cleans the genome again exactly as in `plan` (it is
   deterministic and takes seconds), then takes the plan row and its strand sequence.
   - It sets `Features.SEQUENCE` to the slice and `Features.GC_PERCENT = 0.5`.
   - It maps upstream `scan_oneprocess_1` over the **owned window indices only**,
     `j_local = j − chunk_start` for `j` in the core, with `multiprocessing.Pool(task.cpus)`.
   - It keeps `j ≥ 0` and maps `scan_oneprocess_2`.
   - It then runs the same serial per-row `scaler.transform(np.array(x).reshape(1, -1))` and
     `model.predict_proba(x)[0, 1]`, keeping `pred ≥ 0.5`. Rows are **not** batched: batching could
     change the last bits of `pred`, and the TSV prints `repr` (backlog item, item 10).
   - It shifts `start` and `stop` by `+chunk_start`.

   Ownership is decided by the window start *j*, not by the hairpin position: the call tuples do not
   carry *j*, so the filter has to be the index list itself. Every window is then evaluated by exactly
   one chunk. Output: `calls_<chunk_id>.pkl` (pickle protocol 4) holding `chunk_id`, `core_start`,
   `core_end` and the call tuples in window order, with upstream types unchanged (`np.float64` `pred`).
   The model is loaded from `model.pkl` in the cwd (symlinked as in the unsharded module). The chunk
   sets no BLAS thread variable, so its environment matches the unsharded module.
5. **MERGE (`merge`).** For each (seq_index, strand) in upstream output order (every + block in FASTA
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
6. **Nextflow.** New modules `workflow/modules/local/terminatornet_{plan,chunk,merge}/main.nf` (conda
   `envs/terminatornet.yml`, the same container) and a subworkflow
   `workflow/subworkflows/local/terminatornet_sharded.nf`:

   | Process | Input | Output | Resources (base.config `withName`) |
   |---|---|---|---|
   | `TERMINATORNET_PLAN` | `tuple val(meta_id), path(genome, stageAs: 'genome.fa')` | `path 'plan.tsv'` | `process_single` (1 CPU, 2 GB) |
   | `TERMINATORNET_CHUNK` | `tuple val(chunk_id), path(genome, stageAs: 'genome.fa'), path(plan)`; tag `chunk_id` | `path "calls_${chunk_id}.pkl"` | `cpus = 2; memory = { 4.GB * task.attempt }`; FarmShare: `time = { 2.h * task.attempt }` |
   | `TERMINATORNET_MERGE` | `tuple val(meta_id), path(genome, stageAs: 'genome.fa'), path(plan)`, `path(calls)` (collected) | `path 'terminatornet.tsv', emit: tsv`; `path 'versions.yml'` | `process_single` |

   - The chunk channel is `plan.splitCsv(header: true, sep: '\t').map { it.chunk_id }`, combined with
     the genome and plan.
   - The merge gets `CHUNK.out.calls.collect()`. It reads the files in plan order, so the channel
     order does not matter.
   - annotate.nf swaps `TERMINATORNET(gin)` for `TERMINATORNET_SHARDED(gin)` and keeps
     `kv('terminatornet', … .tsv)`, so `TERMINATOR_LAYER` and its inputs are unchanged.
   - The 'TERMINATORNET' `withName` line in base.config (16 CPUs) is removed.
   - The container, conda and prebuilt-env selectors are widened to `TERMINATORNET.*`.
   - The run-review process→step mapping maps all three processes to `terminatornet`.
   - The old `modules/local/terminatornet` module stays, unchanged, as the unsharded reference for
     `workflow/tests/real/terminatornet.nf` and the acceptance comparison.
   - Stubs: PLAN writes `params.stub?.terminatornet_chunks ?: 1` rows, CHUNK touches its pickle, and
     MERGE writes `genome\tNO TERMINATORS FOUND` as the old stub does.
7. **Retry and resume.** The 137/140/143 policy applies per chunk, so a timeout re-runs one ~10-min
   chunk instead of the whole scan. Chunk hashes depend only on the genome, the plan and `chunk_id`, so
   `-resume` keeps finished chunks.
8. **Topology.** TerminatorNet treats every sequence as linear. Circular replicons are not wrapped:
   there are no windows across the origin, exactly as in the unsharded run.
9. **R14 impact.**
   - Chunk count = Σ over contigs of 2·⌈nwin/CORE⌉: about 40 for 5 Mb, about 56 for msmeg (≈ 7 Mb),
     2 for sim100k's plasmid plus 2 for its chromosome.
   - Each chunk takes about 76 CPU-min/Mb-strand × 0.25 Mb ≈ 19 CPU-min, so about 10 min at 2 CPUs
     (from ADR-0062's 50.3 min × 15.1 CPUs ÷ 10 Mb-strands).
   - With FarmShare `queueSize = 8`, the step uses at most 8 × 2 = 16 CPUs at a time, down from
     today's one task of up to 16. The run's bound of 8 × 16 = 128 CPUs is unchanged.
   - Submitting 40–56 jobs at `submitRateLimit = 10/1min` takes 4–6 min, overlapped with execution.
   - Memory falls from 16 GB in one task to 4 GB per chunk (2 CPUs under the 4000 MB/CPU rule).
10. **Out of scope (backlog).** Batched `predict_proba`, a config key for CORE, scanning across the
    origin of circular replicons, and the `-n 1` ordering difference (Q2).

## Examples (CORE = 250,000, PAD = 200 unless stated)
| # | Input | Expected |
|---|---|---|
| 1 | one linear contig, L = 1,000,000, any mode | nwin = 999,941. + cores `[0,250000) [250000,500000) [500000,750000) [750000,999941)`. Slices `[0,250200) [249800,500200) [499800,750200) [749800,1000000)`. The same 4 cores on − in RC coordinates. 8 plan rows `c00000p0000…c00000m0003` |
| 2 | plasmid shorter than one chunk, L = 120,000 | one core per strand `[0,119941)`, slice `[0,120000)`; 2 rows |
| 3 | contig shorter than a window, L = 40 | nwin = 1: core `[0,1)`, slice `[0,40)`. The single upstream window uses the shortened hairpin `s[0:32]` |
| 4 | empty record (`>e` with no sequence) between two contigs | no plan row and no output line for `e`. `seq_index` counts only the records upstream keeps |
| 5 | 2-nt record `>t AC` | cleaned to `''`: core `[0,1)`, slice `[0,0)`. Output `t\tNO TERMINATORS FOUND` in the + section and again in the − section |
| 6 | boundary: + strand, CORE = 1,000, a hairpin from window j = 995 at `[1003, 1030]` (U-tract to 1,038) | evaluated only by chunk 0 (its slice reaches 1,200 ≥ j+67 = 1,062). Chunk 1's slice `[800, 2200)` contains it, but chunk 1 does not scan j = 995. The call appears once, exactly as in the single-chunk run |
| 7 | tie across a chunk edge: call A (window in core 0; start 980, stop 1,010, pred 0.93) and call B (core 1; start 1,000, stop 1,030, pred 0.93) | one merged call, B (`pred >= best_pred`: the later one wins), the same as the unsharded run. With B.pred = 0.92 the result is A |
| 8 | lower-case letters and Ns: `>c ACGTNNNNacgt…` | plan coordinates are in the cleaned (N-stripped, upper-case) sequence. Output equals the unsharded output (probe values above) |
| 9 | the same chunk run at `task.cpus` 1, 2 and 8 | byte-identical `calls_*.pkl` |
| 10 | genome whose records are all < 501 nt (`mode = short`) | no stop sort, calls in window order. Output equals the unsharded `-n 2` run |

## Tests (test-author)
- **Fast (pure, file loaded by path).**
  - Plan rows for examples 1–5: cores partition `[0, nwin)` with no gap or overlap, and slices are
    clamped to `[0, L)`.
  - Chunk ids are unique and in sort order.
  - With a fake deterministic per-window scorer (for example, a hash of `s[j−8:j+67]` gives a call
    with start/stop/pred), the sharded result equals the single-chunk result. This is a property test
    over random sequences of length 0–3,000, CORE ∈ {1, 7, 61, 1000}, PAD ∈ {67, 200}, both strands
    and both modes. Make it fail with PAD = 66 (right reach) and with a 7-nt left reach.
  - Merge order (+ blocks, then − blocks) and validation (a missing, extra or mismatched pickle exits 1).
- **`needs_conda` (Slurm).** Examples 3, 5, 6, 7, 8 and 10 with the real model on small synthetic
  FASTAs: sharded output, compared with `cmp`, against `TerminatorNet.py -s -n 2`. Then sim100k at
  CORE = 10,000: byte-identical to `tests/data/recorded/terminators/terminatornet/sim100k/terminatornet.tsv`.
- **nf-test (stub).** 3 stub chunks give 3 CHUNK tasks, and MERGE emits `terminatornet.tsv`. The lead
  updates the existing tests that name `TERMINATORNET`: `annotate.nf.test`,
  `main_genetable_context.nf.test`, `lib/annotate_stub.groovy`, `test_tool_resources_pins.py`,
  `test_nf_conda_profile.py`, `test_m11_06_*`, `test_m11_14_images_json.py` and the runner trace mapping.

## Acceptance
- Examples 1–10 hold, and the fast suite, nf-test, ruff and mypy pass.
- `terminatornet.tsv` from the sharded pipeline (`-profile farmshare,apptainer`, real CORE) is
  byte-identical (`cmp`) to the unsharded module's output at the same pin for **sim100k, Nissle 1917,
  mhyo and gsulf**. The references are `$GIPPA_ROOT/recordings/terminatornet/<sample>/terminatornet.tsv`
  where they exist, otherwise a fresh run of `modules/local/terminatornet`. The sha256 values go in
  the PR.
- On **msmeg**, PROGRESS.md records the wall time of the sharded step (first PLAN submit to MERGE end),
  the chunk count, the max chunk wall time, max RSS and peak concurrent CPUs (from the trace), next to
  the failed unsharded attempt (exit 140 at 12 h).
- No run exceeds the R14 caps (≤ 128 CPUs, ≤ 20 jobs per user across GIPPA jobs).
- The PR adds `changelog.d/<PR>.perf.md` (or the matching type).

## Open questions for the lead
1. **Chunk CPUs vs `queueSize = 8`.** With 2-CPU chunks the step never uses more than 16 CPUs at once,
   the same as today. The gain is then limited to:
   - the parallel `predict_proba` (today ~6 % serial: 1,508 % of 1,600 %);
   - shorter queue waits for small jobs (ADR-0062: ≤ 0.2 min vs 6 min);
   - slow-node isolation and per-chunk retry, which fixes msmeg's failure mode.

   A real cut in wall time needs `cpus × 8 > 16`: 4-CPU chunks give ≤ 32 CPUs, 8-CPU chunks ≤ 64, both
   within R14. The chunks also take all 8 slots from other ready tasks. Keep 2 (the brief) and decide
   after the msmeg measurement, or start at 4?
2. **The `-n 1` ordering.** Unsharded `-n 1` uses `scan()` (no stop sort) even for long sequences, so
   its output can differ from `-n ≥ 2`. ADR-0062's "output does not depend on `-n`" holds only for
   `-n ≥ 2`. The sharded step reproduces `-n ≥ 2` at any `--threads`, including 1. Is that the
   intended reference? An ADR-0062 note may be needed.
3. **Slice vs whole sequence.** A chunk could set `Features.SEQUENCE` to the whole strand sequence and
   scan only its owned global *j*. That needs no PAD or coordinate shift and gives the same output by
   construction, at the same read cost (every chunk re-reads the FASTA anyway). This spec keeps the
   approved slice design. Should it switch?
4. **ADR.** Do the new process names, the removed `withName 'TERMINATORNET'` line and the always-on
   sharding (no config switch, since defaults.yaml is protected) need an ADR or an ADR-0062 amendment,
   or is this spec enough? Is an escape hatch back to the unsharded module wanted?
5. **CORE = 250,000.** This sets the chunk to about 10 min at 2 CPUs, the tail to ≤ 1 chunk and 40–56
   Slurm jobs per run. Is that job count acceptable on FarmShare, or should CORE be 500,000 (half the
   jobs, about a 20-min tail)?
6. **The mhyo/gsulf/msmeg references.** Are the #411 unsharded `terminatornet.tsv` files kept under
   `$GIPPA_ROOT` (and at pin 72180309)? If not, the acceptance runs need fresh unsharded jobs (about
   1–2 h at 16 CPUs each for gsulf; msmeg needs > 12 h and is only timed, not compared).
