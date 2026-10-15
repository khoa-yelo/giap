# M2-07: Nextflow processes for the external annotation tools and the `ANNOTATE` subworkflow

Milestone M2 · depends on: M2-01 (env files, container pins), M2-06 (skeleton) · estimate ≈ 400 code
lines (Groovy; ≈ 30 per module × 10 + subworkflow ≈ 100) · capabilities: `needs_nextflow` (stub mode only)

## Spec sections implemented

- §6.3: "One process per external tool, and one per Python step or step group." Each declares `tag`,
  a `label`, `conda` + `container`, `versions.yml`, and "a `stub:` block that writes minimal valid
  outputs, used by nf-test in stub mode". "Memory-heavy tools (Bakta, geNOMAD, antiSMASH, genome-wide
  Promotech) declare memory". "**Optional modules.** `mobile`, `bgc` and `function` are switched with
  `params.modules`. Downstream steps accept a missing optional layer and record `layer_missing` in the manifest."
- §6.1 data flow (annotation track): "Bakta (or --annotation JSON) ─► gene model ─┬─► UniOP → pairs → TUs
  ├─► Promotech (targeted windows) ├─► Promoter Calculator (TU upstream) ├─► TerminatorNet (+ TransTermHP;
  BacTermFinder opt-in) ├─► OSTIR (anti-SD from 16S) ; codon-bias └─► antiSMASH [bgc]"; "geNOMAD [mobile]
  (parallel with Bakta)".
- §7.1.1 commands (verbatim per tool; see table below); §7.1.4 "chunk them into files of 100,000
  windows and run the chunks in parallel"; G11 Promoter Calculator "once per TU upstream region, in parallel".
- §6.4 / §5.2 `annotation.precomputed`: a given layer skips its tool.
- §17: Apptainer fixes "bind a writable directory, set `--containall` / `HOME` … record the fix in the module's `meta.yml`".

## Public interface

`workflow/modules/local/<tool>/{main.nf,meta.yml}` for `bakta, uniop, promotech, promoter_calculator,
terminatornet, transtermhp, bactermfinder, ostir, genomad, antismash`; `workflow/subworkflows/local/annotate.nf`;
`workflow/tests/annotate.nf.test`; `tests/unit/test_nf_cmd_sync.py`.

| Process | Label | Script (essentials, from §7.1.1) | Stub outputs |
|---|---|---|---|
| `BAKTA` | `process_high` | `bakta --db $db --output bakta --prefix $prefix --threads $task.cpus --keep-contig-headers (--complete \| --replicons $replicons) [--regions $regions] $genome` | `bakta/<P>.json` copied from `assets/stub/bakta.json` (minimal valid JSON), empty `.gff3 .gbff .faa .tsv` |
| `UNIOP` | `process_low` | `python $UNIOP_HOME/src/UniOP -f $gff --gff_source prodigal_gff -t out` | `uniop.pred` header only |
| `PROMOTECH` (per chunk) | `process_medium` | `python $PROMOTECH_HOME/promotech.py -s -m RF-HOT -f $chunk -o out/` | `promotech_NNN/` with header-only native file |
| `PROMOTER_CALCULATOR` (per TU) | `process_single` | `promoter-calculator -i $region -o ${tu_id}.csv -j 1` | header-only CSV |
| `TERMINATORNET` | `process_medium` | `python $TERMINATORNET_HOME/TerminatorNet.py -f genome.fa -s -n $task.cpus -o terms.tsv` (genome staged as `genome.fa`) | `NO TERMINATORS FOUND` line |
| `TRANSTERMHP` | `process_single` | `transterm -p $EXPTERM genome.fa genes.coords > transtermhp.tt` | empty `.tt` |
| `BACTERMFINDER` (opt-in) | `process_medium` | per its README at the pinned commit | empty dir |
| `OSTIR` | `process_medium` | `ostir -i ostir_input.csv -o ostir.csv -j $task.cpus` | header-only CSV |
| `GENOMAD` [mobile] | `process_high` | `genomad end-to-end --cleanup --splits 8 --threads $task.cpus genome.fa out $db/genomad_db` | empty summary TSVs with headers |
| `ANTISMASH` [bgc] | `process_high` | `antismash --genefinding-tool none --minimal --cpus $task.cpus --output-dir out $gbff` | JSON `{"records": []}` |

Every process: `tag "${meta_id}"`; `conda "${projectDir}/envs/<tool>.yml"`; container from
`conf/containers.config` (`withName`); emits `versions.yml` (`"${task.process}": {<tool>: <version>}`),
the version read from the tool where it reports one, else the pinned value. Memory-heavy processes set
`memory` explicitly (Bakta 48.GB, geNOMAD 48.GB, antiSMASH 32.GB).

`ANNOTATE` (subworkflow) wires the index's step contract: `GIAP_STEP` calls for py steps;
`params.annotation.precomputed.<key>` non-null → `Channel.fromPath(<file>)` replaces the tool and its
parse step; `params.modules` without `mobile`/`bgc` → `GENOMAD`/`ANTISMASH` not run and
`region_layer` gets no optional input; `bactermfinder` only when `params.annotation.terminators.bactermfinder`.
Promotech chunks and Promoter Calculator regions scatter with `flatten()` and gather sorted by name
(`collect(sort: true)`) so downstream inputs are order-stable.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `nf-test … annotate.nf.test -stub`, all modules on, no precomputed | every process above (except `BACTERMFINDER`) runs once; `GIAP_STEP` runs for each py step of the contract table |
| 2 | stub with `precomputed.bakta_json = sim/bakta_like.json` | `BAKTA` not run; `parse_bakta` receives that file |
| 3 | stub with `modules = [core]` | no `GENOMAD`, no `ANTISMASH`; `region_layer` still runs |
| 4 | stub with 3 window chunks and 5 TU regions | `PROMOTECH` × 3, `PROMOTER_CALCULATOR` × 5; `promoter_layer` input list sorted `promotech_000 … 002` |
| 5 | `terminators.bactermfinder = true` | `BACTERMFINDER` × 1 |
| 6 | `test_nf_cmd_sync.py`: `terminatornet_cmd.terminatornet_command(...)` flags vs. `modules/local/terminatornet/main.nf` script | every flag of the builder (`-f`, `-s`, `-n`, `-o`) appears in the script block; `-s` present |
| 7 | grep every `main.nf` | has `stub:`, `versions.yml`, `label 'process_*'`, `conda`, `tag` |

## Edge cases

- TerminatorNet's parser matches flags by substring (§7.1.1): staged names have no spaces, only
  documented flags are passed; covered by example 6.
- A failing optional-module process fails the run (exit 4) unless the module is off; no silent skip.
- Bakta with `--regions` only when `--annotation` is a GFF3/GBFF (M2-11 output present).
- `meta.yml` records licence (R5) and any Apptainer fix (`--containall`, writable `HOME`).

## Acceptance criteria

- [ ] `workflow/tests/annotate.nf.test`: examples 1–5 in stub mode (CI and session).
- [ ] `tests/unit/test_nf_cmd_sync.py` (example 6) for every tool whose builder exists; tools whose
  builder lands later (M2-12…M2-18) add their pair to the test in their own PR.
- [ ] `tests/unit/test_nf_modules_lint.py` (example 7) and: every `withName` in `containers.config`
  matches a process; every process name maps to a `TOOLS` key (M2-01).

## Questions for the lead

1. Command builders in Python vs. scripts in Groovy duplicate the command. **Proposed default:** keep
   both (tool images do not contain GIAP); example 6 keeps them in sync.
2. Stub Bakta JSON: **Proposed default:** `workflow/assets/stub/bakta.json` = a 3-gene JSON from the
   simulator (≤ 5 kB), so downstream stubs run real parsing in stub mode.

## Lead decisions (implementation)
1. Example 4's order check cannot read `ln -s` staging lines: Nextflow 26 stages inputs in descending name
   order regardless of channel order. The `GIAP_STEP` stub writes its `--in` entries in order to
   `.giap_inputs` (one `key=path` per line); the nf-test reads that (`H.stepInputs`).
2. ADR-0007 size: 619 code lines → two stacked PRs. PR A: the 10 tool modules (workflow/modules/local/<tool>/
   incl. meta.yml) + stub asset, with the module lint and cmd-sync tests. PR B: annotate subworkflow,
   main.nf, giap_step, configs, `giap run` wiring, with the remaining M2-07 tests (stacked on A).

## Lead decisions after code review (PR #65, changes requested — real-tool run failed)
Fixed in this PR (#65):
1. ANNOTATE keys follow the merged M2-12 step registry (`steps/operons.py`): `uniop_input` gets `replicons`;
   `operon_pairs` gets `uniop_cds`, `replicons`, `genome` and emits `operon_meta`; `tu_upstream` is a
   directory; the precomputed branch passes only declared keys. A test checks every ANNOTATE `py step` call's
   inputs/outputs against the step registry (all registered steps, not just M2-12).
2. Bakta auto `--complete`: replicons.tsv topology column holds `c`/`l` (preflight); the check compares to
   `c`, and the test fixture uses the real preflight format.
3. Containers without `ps` (bakta, antiSMASH): Nextflow's task metrics need `ps`. Fix without new images:
   the apptainer profile binds the host's `ps` (procps) into those two processes only if a test run proves
   it works; otherwise disable per-task metrics for them (`process.withName` … or trace fields that avoid
   `ps`) — implementer verifies with a real run on Slurm. `meta.yml` notes updated.
4. `farmshare,apptainer` turns conda on for the conda-only tools (uniop, promotech, promoter_calculator,
   terminatornet, ostir), using the prebuilt envs under $GIAP_ROOT (micromamba, `conda.useMicromamba`),
   and containers for the rest. Config test.
5. Promotech: the env (workflow/envs/promotech.yml) gains the missing PyQt5 dependency (or the
   headless equivalent that makes `import PyQt5` succeed), lock updated, env rebuilt via
   `sbatch scripts/slurm/build_envs.sbatch promotech`. Stub writes `sequences_predictions.csv`
   (tab, `CHROM SEQ PRED`) like the real tool.
6. Promoter Calculator: when the tool writes no CSV (no promoter), the process writes the real 27-column
   header; the stub uses the same header.
7. Stub headers equal the real outputs (OSTIR incl. `RBS_distance_bp` and order; geNOMAD virus_genes.tsv 20
   columns) — a test compares stub headers with the recordings under $GIAP_ROOT/recordings (skip if absent)
   and with a committed copy of each real header.
8. Precomputed bakta_json without a GBFF: antiSMASH is skipped and `layer_missing` recorded.
9. Bakta `--regions` file is a staged path input (part of the task hash).
10. Acceptance: a real-tool ANNOTATE run on SIM100K via Slurm (`scripts/slurm/` script, work dir under
    $GIAP_ROOT/scratch — never node-local /tmp) completes every enabled tool; Bakta output equals the
    recording; result noted in the PR.
Follow-up PR (stacked, branch m2/run-wiring): `giap run` writes `inputs.replicons` (runs preflight) and
`params.db_dir` so main.nf takes the ANNOTATE branch (review item 1), and `giap _step` calls
`load_step_modules()` (item 5, pre-existing). Tests for both.

## Lead decisions after the re-review (PR #65, changes requested)
11. Bakta always gets `--replicons <replicons.tsv>` when preflight wrote it (its `c`/`l` topology and the
    replicon types come from that table, so the chromosome is never re-typed by length); `--complete` stays
    only for the all-circular case per §7.1.1. `bakta_cmd.py` and `bakta/main.nf` agree (cmd-sync test).
    The real-run comparison includes `sequences[]` types/topology (follow-up harness PR).
12. BACTERMFINDER (conda-only, opt-in) gets a prebuilt-env entry like the other five conda-only tools
    (six in total); conda is never built on the fly under farmshare or farmshare,apptainer.
13. The `ps` shim bind is read-only (`:ro`). Decision 3 is amended: the shim is a `/proc`-based script
    (`workflow/assets/ps/ps`) printing `pid ppid` for Nextflow's process tree, bound for BAKTA and ANTISMASH.

### Lead decisions for the run-wiring follow-up (branch m2/run-wiring)
1. `giap run` runs preflight in full, including the evidence rule (`inputs.evidence_present`): a run with
   neither peptides nor proteins exits 1 before Nextflow starts. Older genome-only `giap run` tests are
   updated to pass evidence (test-author).
2. The test-author's interface choices are accepted: absolute `inputs.replicons` in preflight's 5-column
   c/l format matching the normalised genome's IDs (the run passes the normalised genome); `params.db_dir`
   = --db-dir, else $GIAP_DB_DIR; a found step without required `--in` fails exit 5 "required input".
3. After the code review of #78 (real `giap run`, Slurm 1848591):
   a. `--resume` must cache: preflight rewrites `work/preflight/*` only when the content changed (unchanged
      files keep their mtime); a stub test runs `giap run` twice with `--resume` and asserts cached > 0.
   b. The DB directory is visible in containers: the apptainer profile binds `params.db_dir` read-only
      (`apptainer.runOptions`/`containerOptions` `-B <db_dir>:<db_dir>:ro`) when set. Config test.
   c. `--profile` accepts comma-joined known profiles (e.g. `farmshare_local,apptainer`,
      `farmshare,apptainer`); unknown names still exit 1.
   d. Order: preflight runs before the `nextflow` presence check (a bad genome exits 1 with a report even
      without Nextflow).
   e. `db_dir` is passed to preflight/`read_manifest` explicitly; no `os.environ` mutation.
   f. Tests: `--strict` warning exits 3 before launch; missing manifest exits 2 with no launch.
4. After the re-review of #78 (real `giap run`, Slurm 1848630: exit 2 at `environment.tools`):
   a. `environment.tools` checks binaries/environments from the M2-01 tool registry, not from the DB
      manifest: a tool is present when its pinned container image exists under `$GIAP_ROOT/containers`
      (apptainer profiles) or, for conda-only tools, its prebuilt micromamba env exists with a build
      marker matching `workflow/envs/<tool>.yml`; under the `test`/stub profiles tools are not required.
   b. Databases and models come from the `giap setup` manifest, with item names mapped to requirement
      names: `cog2020_def`+`cog2020_fun` → `cog`, `go_basic`+`goslim_generic` → `go`, `kegg_*` → `kegg`,
      `promotech_models` → `promotech` (model), `contaminants` → `contaminants` (all parts of a group
      must be ok).
   c. Tests use a manifest in the real `giap setup` format (copy the structure of a real run, versions
      faked) and fake container/env directories; never the real $GIAP_DB_DIR or $GIAP_ROOT.
   d. Acceptance: a real `giap run` on SIM100K via Slurm passes preflight and reaches ANNOTATE.
