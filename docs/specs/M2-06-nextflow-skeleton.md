# M2-06: Nextflow skeleton: configuration, profiles, labels, `GIAP_STEP`, launcher, stub test

Milestone M2 · depends on: M2-02 · estimate ≈ 330 code lines (Groovy + Python; tests not counted) ·
capabilities: `needs_nextflow` (nf-test stub mode; available in the session allocation and in CI)

## Spec sections implemented

- §6.3 Launch: "`giap run` validates inputs and config in Python, writes `OUT/work/params.json` and
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
  `giap-<process>`, `process.scratch = '/tmp'` (node-local), `resourceLimits = [cpus: 64, memory:
  240.GB, time: 48.h]`"); `farmshare_local` ("local executor inside the session allocation, capped at
  `SLURM_CPUS_ON_NODE` CPUs and 30 GB"); `test`; `standard` ("local executor, all cores").
- §6.3 Error strategy: "Retry once on exit codes 137 or 143 with memory × 2; otherwise terminate,
  and `giap run` exits 4 naming the failed step and its log." Publishing: "`publishDir` copies final
  files into `OUT/` (mode `copy`)".
- §3.8: "`scripts/slurm/giap_run.sbatch`: 2 CPU, 8 GB, up to 2 days" with `-profile farmshare,apptainer`.

## Public interface

Files: `workflow/main.nf` (entry workflow; calls `ANNOTATE` from M2-07 when present, else only
`GIAP_STEP` smoke), `workflow/nextflow.config` (params defaults: `threads`, `memory`, `outdir`,
`modules`, `annotation`; profiles; `includeConfig`), `workflow/conf/{base,test,farmshare,
farmshare_local,standard}.config`, `workflow/modules/local/giap_step/main.nf` + `meta.yml`,
`workflow/tests/main.nf.test` (stub), `scripts/slurm/giap_run.sbatch`.

`GIAP_STEP` process: input `tuple val(step), path(inputs, stageAs: …)`, `val(in_keys)`,
`val(out_spec)`; script `giap _step ${step} --config ${resolved} --in k=path … --out k=path …`;
stub writes an empty file (or empty Parquet via `giap _step --stub`, see Q2) for each output.

`src/giap/runner/nextflow.py`
```python
def write_params(outdir: Path, cfg: GiapConfig, inputs: Mapping[str, Path]) -> tuple[Path, Path]  # params.json, resolved_config.yaml
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
| 8 | `nf-test test workflow/tests/main.nf.test --profile test` (stub) | passes; `GIAP_STEP` stub outputs exist; `versions.yml` emitted |

## Edge cases

- Job names: `executor.jobName = { "giap-" + task.process.tokenize(':')[-1].toLowerCase() }` (R14:
  every GIAP job name starts with `giap-`); test via `nextflow config` output.
- `params.threads`/`params.memory` come only from `params.json`; no `params.*` has a CLI-derived
  default that differs from `defaults.yaml` (a test compares them).
- No `plugins {}` block anywhere in `workflow/` (grep test).
- `farmshare_local` without `SLURM_CPUS_ON_NODE` falls back to 4 CPUs and warns.

## Acceptance criteria

- [ ] `tests/unit/test_runner_nextflow.py`: examples 1–4 (pure Python, no Nextflow).
- [ ] `tests/integration/test_nf_config.py` (`needs_nextflow`): examples 5–7 by parsing
  `nextflow config -flat` output.
- [ ] `workflow/tests/main.nf.test` (example 8) runs in `ci.yml` (stub mode) and in the session.
- [ ] `scripts/slurm/giap_run.sbatch` follows the §3.8 header with 2 CPU, 8 GB, `--time=2-00:00:00`.

## Questions for the lead

1. `docker` profile is required by §6.3 but cannot be tested on FarmShare. **Proposed default:**
   define it (3 lines) and test only `nextflow config` parsing.
2. Stubs for Parquet outputs: an empty file is not valid Parquet. **Proposed default:** `giap _step
   NAME --stub` writes schema-valid empty outputs (zero-row Parquet with the declared columns); it
   needs GIAP in the stub environment, which the `GIAP_STEP` process always has.

## Lead decisions (implementation, accepted)
1. `resolved_config.yaml` is written to `OUT/` (next to the results), not `OUT/work`.
2. `params.memory` defaults to `"64 GB"` (the CLI default; defaults.yaml has no memory key).
3. The `GIAP_STEP` stub only `touch`es declared outputs (CI runners have no `giap`); the `smoke` step
   behind main.nf exists only under `-stub-run` until real registry steps are wired (M2-07).
4. The farmshare_local warning prints via an Elvis-chain `println` because Nextflow 26's config parser
   rejects `if` in config files.
