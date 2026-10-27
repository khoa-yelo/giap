# M2-06: Nextflow skeleton: configuration, profiles, labels, `GIPPA_STEP`, launcher, stub test

Milestone M2 · depends on: M2-02 · estimate ≈ 330 code lines (Groovy + Python; tests not counted) ·
capabilities: `needs_nextflow` (nf-test stub mode; available in the session allocation and in CI)

## Spec sections implemented

- §6.3 Launch: "`gippa run` validates inputs and config in Python, writes `OUT/work/params.json` and
  `resolved_config.yaml`, then launches: `nextflow run <pkg>/workflow/main.nf -params-file
  params.json -profile <profile> -work-dir OUT/work -with-trace OUT/logs/trace.tsv -with-timeline
  OUT/logs/timeline.html -with-report OUT/logs/nf_report.html [-resume]`. Nextflow never parses user
  parameters on its own. Do not use Nextflow plugins in version 1."
- §6.3 Process conventions: `tag`; labels "`process_single` (1 CPU, 2 GB), `process_low` (2 CPU,
  6 GB), `process_medium` (8 CPU, 16 GB) or `process_high` (16 CPU, 48 GB)"; `conda` and
  `container`; "nf-core style `versions.yml`"; "a `stub:` block that writes minimal valid outputs".
- §6.3 Resources: "`process.resourceLimits = [cpus: params.threads, memory: params.memory]`".
- §6.3 Profiles: `apptainer` ("`apptainer.enabled = true`, `apptainer.autoMounts = true`, images
  `docker://…@sha256:…` cached in `$NXF_APPTAINER_CACHEDIR`"), `conda` ("`conda.useMicromamba =
  true`"), `docker`; `farmshare` ("`executor = 'slurm'`, `queue = 'normal'`, `clusterOptions =
  '--qos=normal'`, `executor.queueSize = 20`, `executor.submitRateLimit = '10/1min'`, job names
  `gippa-<process>`, `process.scratch = '/tmp'` (node-local), `resourceLimits = [cpus: 64, memory:
  240.GB, time: 48.h]`"); `farmshare_local` ("local executor inside the session allocation, capped at
  `SLURM_CPUS_ON_NODE` CPUs and 30 GB"); `test`; `standard` ("local executor, all cores").
- §6.3 Error strategy: "Retry once on exit codes 137 or 143 with memory × 2; otherwise terminate,
  and `gippa run` exits 4 naming the failed step and its log." Publishing: "`publishDir` copies final
  files into `OUT/` (mode `copy`)".
- §3.8: "`scripts/slurm/gippa_run.sbatch`: 2 CPU, 8 GB, up to 2 days" with `-profile farmshare,apptainer`.

## Public interface

Files: `workflow/main.nf` (entry workflow; calls `ANNOTATE` from M2-07 when present, else only
`GIPPA_STEP` smoke), `workflow/nextflow.config` (params defaults: `threads`, `memory`, `outdir`,
`modules`, `annotation`; profiles; `includeConfig`), `workflow/conf/{base,test,farmshare,
farmshare_local,standard}.config`, `workflow/modules/local/gippa_step/main.nf` + `meta.yml`,
`workflow/tests/main.nf.test` (stub), `scripts/slurm/gippa_run.sbatch`.

`GIPPA_STEP` process: input `tuple val(step), path(inputs, stageAs: …)`, `val(in_keys)`,
`val(out_spec)`; script `gippa _step ${step} --config ${resolved} --in k=path … --out k=path …`;
stub writes an empty file (or empty Parquet via `gippa _step --stub`, see Q2) for each output.

`src/gippa/runner/nextflow.py`
```python
def write_params(outdir: Path, cfg: GippaConfig, inputs: Mapping[str, Path]) -> tuple[Path, Path]  # params.json, resolved_config.yaml
def launch_command(workflow: Path, params: Path, profile: str, outdir: Path, resume: bool) -> list[str]
def failed_step(trace_tsv: Path) -> tuple[str, str] | None        # (step, workdir) of the first FAILED/ABORTED row
def run_nextflow(cmd: Sequence[str], outdir: Path) -> ExitCode    # 0, or 4 with message naming step + log
```

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `launch_command(W, OUT/"work/params.json", "farmshare,apptainer", OUT, resume=False)` | `["nextflow","run",str(W),"-params-file",str(OUT/"work/params.json"),"-profile","farmshare,apptainer","-work-dir",str(OUT/"work"),"-with-trace",str(OUT/"logs/trace.tsv"),"-with-timeline",str(OUT/"logs/timeline.html"),"-with-report",str(OUT/"logs/nf_report.html")]` |
| 2 | same with `resume=True` | the same list + `"-resume"` |
| 3 | `write_params` twice with the same config | byte-identical `params.json` (sorted keys, no timestamps, relative paths only where the spec allows) |
| 4 | trace with rows `…\tCOMPLETED\t…ANNOTATE:PARSE_BAKTA…` and `…\tFAILED\t…ANNOTATE:TERMINATORNET (chr1)…` | `failed_step` → `("terminatornet", "<workdir hash path>")`; `run_nextflow` returns 4 and prints `step terminatornet failed; log: <workdir>/.command.log` |
| 5 | `nextflow config workflow -profile farmshare` | `process.executor = 'slurm'`, `process.queue = 'normal'`, `process.clusterOptions = '--qos=normal'`, `executor.queueSize = 20`, `executor.submitRateLimit = '10/1min'`, `process.scratch = '/tmp'`, `process.resourceLimits = [cpus:64, memory:240.GB, time:48.h]` |
| 6 | `SLURM_CPUS_ON_NODE=8 nextflow config workflow -profile farmshare_local` | `executor = 'local'`, `resourceLimits = [cpus: 8, memory: 30.GB]` |
| 7 | process with `label 'process_medium'` under `base.config` | `cpus = 8`, `memory = 16.GB`; on exit 137 → retried once with 32.GB; exit 1 → terminate |
| 8 | `nf-test test workflow/tests/main.nf.test --profile test` (stub) | passes; `GIPPA_STEP` stub outputs exist; `versions.yml` emitted |

## Edge cases

- Job names: `executor.jobName = { "gippa-" + task.process.tokenize(':')[-1].toLowerCase() }` (R14:
  every GIPPA job name starts with `gippa-`); test via `nextflow config` output.
- `params.threads`/`params.memory` come only from `params.json`; no `params.*` has a CLI-derived
  default that differs from `defaults.yaml` (a test compares them).
- No `plugins {}` block anywhere in `workflow/` (grep test).
- `farmshare_local` without `SLURM_CPUS_ON_NODE` falls back to 4 CPUs and warns.

## Acceptance criteria

- [ ] `tests/unit/test_runner_nextflow.py`: examples 1–4 (pure Python, no Nextflow).
- [ ] `tests/integration/test_nf_config.py` (`needs_nextflow`): examples 5–7 by parsing
  `nextflow config -flat` output.
- [ ] `workflow/tests/main.nf.test` (example 8) runs in `ci.yml` (stub mode) and in the session.
- [ ] `scripts/slurm/gippa_run.sbatch` follows the §3.8 header with 2 CPU, 8 GB, `--time=2-00:00:00`.

## Questions for the lead

1. `docker` profile is required by §6.3 but cannot be tested on FarmShare. **Proposed default:**
   define it (3 lines) and test only `nextflow config` parsing.
2. Stubs for Parquet outputs: an empty file is not valid Parquet. **Proposed default:** `gippa _step
   NAME --stub` writes schema-valid empty outputs (zero-row Parquet with the declared columns); it
   needs GIPPA in the stub environment, which the `GIPPA_STEP` process always has.

## Lead decisions (implementation, accepted)
1. `resolved_config.yaml` is written to `OUT/` (next to the results), not `OUT/work`.
2. `params.memory` defaults to `"64 GB"` (the CLI default; defaults.yaml has no memory key).
3. The `GIPPA_STEP` stub only `touch`es declared outputs (CI runners have no `gippa`); the `smoke` step
   behind main.nf exists only under `-stub-run` until real registry steps are wired (M2-07).
4. The farmshare_local warning prints via an Elvis-chain `println` because Nextflow 26's config parser
   rejects `if` in config files.

## Lead decisions after code review (PR #37, changes requested)
1. `nextflow.config` sets `trace.fields` explicitly (incl. `task_id,hash,name,tag,status,exit,workdir`);
   the runner reads rows with `.get` and falls back to building the work dir from `hash` under the
   run's work dir. A test uses a trace written with Nextflow's DEFAULT field set (no `workdir`), and an
   integration test drives a real failing stub run through `run_nextflow` (exit 4, step named).
2. `failed_step` = the last FAILED row whose process name has no later COMPLETED row (a retried-then-
   succeeded attempt never counts). Test: step A fails once then succeeds, step B fails → B.
3. R14 under `farmshare`: worst-case concurrent CPUs ≤ 128. `executor.queueSize = 8` and per-task
   `resourceLimits` cpus ≤ 16, memory ≤ 64 GB; a task's cpus are `max(label cpus, ceil(memory / 4 GB))`
   (FarmShare MaxMemPerCPU 4000 MB) and stay within the cap after a memory ×2 retry. A unit test parses
   `nextflow config -profile farmshare` and asserts queueSize × max cpus ≤ 128. Larger-memory steps
   (bigmem) get their own label with `maxForks 1` when they appear (M2-07+).
4. `farmshare_local`: `executor { cpus = SLURM_CPUS_PER_TASK (default 4); memory = SLURM_MEM_PER_NODE
   − 2 GB (default 30 GB) }` so parallel tasks share the allocation, not the host.
5. User `--threads`/`--memory` and profile limits combine as the minimum of the two (no silent override);
   documented in the runner docstring and the spec (no ADR: defaults.yaml unchanged).
6. Tests kill the surviving mutants: input paths absolute in params.json (Nextflow's cwd is OUT),
   missing-trace guard, runner messages on stderr.
7. `apptainer.cacheDir` defaults to `$GIPPA_ROOT/containers` when `NXF_APPTAINER_CACHEDIR` is unset.
   Wiring `gippa run` to the runner is M2-07 (note in PROGRESS); `gippa_run.sbatch` stays.
