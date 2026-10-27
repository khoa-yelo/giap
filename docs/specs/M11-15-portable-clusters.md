# M11-15: portable install and clusters (auto profile, demo, zero-flag setup, generic Slurm, Singularity, step image, `--nf-config`, site configs)

Milestone M11 · **pass 2, after the GIPPA rename** · depends on: M11-01 (wheel
bundles `workflow/`), M11-06 and M11-14 (GHCR images, `images.json`, `image_cache()`, `setup
--images`, 14p profile parity, 14u1/14u2 local UniOP), M11-09 (`docker-smoke.yml`), M11-09b (bare
refs, ADR-0059), M11-10 (Bioconda recipe: referenced, not duplicated), ADR-0060 (image tags, no
retagging), ADR-0064 (rename) · estimate ≈ 2,570 code lines + ≈ 2,680 test lines in 13 code/test
PRs, one CI PR (Khoa) and one docs PR, plus a validation run on a second cluster (table at the
end) · not `stats`, except item 8 AC 6 (byte identity of step outputs between host and image) ·
roles: implementer (`src/gippa/cli.py`, `src/gippa/runner/`, `src/gippa/preflight/`,
`src/gippa/setupdb/`, `src/gippa/data/`, `workflow/`, `scripts/release/`, `scripts/dev/`,
`scripts/docs/`), test-author (`tests/unit/test_m11_15*_*.py`,
`tests/integration/test_m11_15*_*.py`; also updates the existing tests that pin today's defaults,
item 6), docs (`docs/get-started/`, `docs/user-guide/platforms.md`, new
`docs/user-guide/add-your-cluster.md`) · **Khoa: yes**: the CI changes (15n), the first
step-image publish (with 14j), the second cluster, and the questions at the end.

**Naming.** Names follow ADR-0064 (command, package, `GIPPA_*` variables, `gippa-*` job names);
the spec is `GIPPA_SPEC.md` (renamed by the rename tool, Q-K32). `<pkg>` means the PyPI
project name read from `pyproject.toml`, `<owner>` the GHCR owner (`khoa-yelo`, M11-14 Q-K14).

## Goal and headline acceptance

Khoa (2026-10-27): users `pip install` and run the full pipeline on clusters other than Stanford
FarmShare, with the **"lowest barrier to entry, simple installation"**. Default target until Khoa
says otherwise: a laptop/workstation with Docker or Apptainer, or a Slurm cluster with Apptainer or
SingularityCE whose compute nodes may have no internet and share a filesystem with the launching
host.

**Headline acceptance (HA).** A new user, with no root, no `--profile` and no hand-written
config, runs:

```bash
conda install -c bioconda gippa        # or: pip install gippa
gippa demo                             # minutes; prints the report path
gippa setup                            # databases + images for the detected runtime
gippa run --genome g.fna --peptides p.tsv --outdir out
```

and gets `out/report/report.html` with `run_manifest.json` `status: completed`, on (a) Linux with
Docker, (b) Linux with Apptainer, (c) a Slurm cluster with Apptainer or SingularityCE (offline
compute nodes), (d) FarmShare (detected; the `farmshare` profile is used). macOS with Docker
Desktop is (a) under emulation. The pip path additionally needs Java ≥ 17 (item 5); the conda path
brings Java and Nextflow. HA is checked by: CI (a) in 15n, the FarmShare nightly for (b)/(d), and
the 15v validation run for (c).

## Spec sections implemented

- §3 Containers row ("The `docker` profile is kept for users elsewhere"), G19; §3.9 (sizes);
  §5.1 (`gippa setup`, `gippa run` options and defaults); §6.3 Launch ("no plugins"), Processes
  ("Python processes call `gippa _step`"), Profiles; §6.4 (the local executor, used by `demo`);
  §7.0 Environment rows (`environment.tools`); §12 Releases ("containers pinned by digest"),
  `containers.yml`; §14.6; §17 ("Nextflow on NFS … `process.scratch = '/tmp'`").
- R1 (launch, packaging and install only; no analysis change), R5/R11 (licences; UniOP never
  redistributed), R6 (official sources only: the Nextflow launcher), R12 (no personal emails in
  site files), R14 (FarmShare caps stay in `farmshare`).
- `docs/get-started/install.md` › "A Slurm (or other) executor off FarmShare" (replaced) and the
  five-step install (replaced by the four commands).
- Khoa's decisions 2026-10-27: licence stays MIT (Q-K34, closed); "turn everything on and just
  document the license stuff, this mostly targets academic use" (Q-K35, closed; item 6).

## What exists today (main 5f5f6306, 2026-10-27; written with the ADR-0064 names)

| Piece | State |
|---|---|
| Default profile | `gippa run --profile` default `farmshare,apptainer`; `gippa preflight --profile` and `gippa setup --profile` default none. Accepted names: the hard-coded `_PROFILES` (`src/gippa/cli.py:211`); `test`/`stub` rejected. |
| Profiles | `standard`, `test`, `stub`, `farmshare`, `farmshare_local`, `apptainer`, `conda`, `docker` in `workflow/nextflow.config` |
| Slurm | only `conf/farmshare.config`: `queue = 'normal'`, `--qos=normal --exclude=…`, `queueSize = 8`, MaxMemPerCPU 4000 MB CPU rule, caps 16 CPUs / 64 GB / 48 h, `scratch = '/tmp'`, job names `gippa-<process>`, `prebuilt_envs.config` |
| Singularity | no profile; `_uniop_argv` already falls back to a `singularity` binary |
| User Nextflow config | none (`-c` only for `DIR/local_tools.config`); install.md: "`gippa run` has no `-c` option" |
| `gippa setup` | `--db-dir` **required**; defaults (`SetupOptions`, `src/gippa/setupdb/items.py:84`): `bakta light`, `genomad` on, `promotech_models` on; `antismash`, `go`, `kegg`, `cog`, `contaminants`, `uniop` off; `--images` off; images need an explicit `--profile apptainer\|docker` |
| Nextflow/Java | prerequisites (M11-14 Q-K17); `gippa run` exits 2 when `nextflow` is not on `PATH`; Java not checked |
| Demo data | `tests/data/genomes/sim100k.fna.gz` (34,034 B), `docs/tutorials/data/sim100k_peptides.tsv.gz` (68,421 B, synthetic, `scripts/docs/make_demo_peptides.py --check`); neither in the wheel. No `gippa demo`. No precomputed layers for SIM100K committed (gippasim, under `tests/sim/`, writes them in tests) |
| `GIPPA_STEP` | no `container`: every Python step (≈ 25 aliases) runs `gippa _step` with the host's `gippa`. `FUNCTION_STEP`/`SENSITIVITY` start R fgsea (`gippa.io.fgsea_r`) |
| Step inputs outside the task dir | `gippa _step --config ${params.config}` reads `OUT/resolved_config.yaml` by absolute path; precomputed layers are absolute paths in it |
| Images | 7 upstream pins; `images.json` `{}` until 14j. M11-06 Q-L8: "no `gippa` core image" |
| Config tests | unit tests parse configs as text (`tests/unit/test_m11_14_profile_parity.py` reader); `nextflow config -flat` only under `@pytest.mark.needs_nextflow` (`tests/integration/test_nf_config.py`) and in ci.yml's `package` job |

## Definition

### 1. Profile auto-detection (15a)

When `--profile` is absent, `gippa run`, `gippa preflight` and `gippa setup` choose one with a pure
function `detect_profile(probe) -> Detected(profile, reason, hints)`. `probe` is an injected
record of `PATH` lookups, `docker info` success, `SLURM_JOB_ID`, the Slurm cluster name and the
`setup.runtime` that `db_manifest.json` records. The default of `gippa run --profile` becomes
`auto` (shown in `--help`); `gippa preflight` keeps "no tool check" only for `--profile none`.

**Two independent choices, then combined.** The executor comes from the launching host only; the
container runtime may come from what `gippa setup` recorded (the images were pulled for it).
Setup on a transfer node says nothing about the executor of the launching host, so setup records
**only the runtime** (`setup.runtime`, schema minor), never an executor.

*Executor, first match wins:*

| # | Condition | Executor | Reason printed |
|---|---|---|---|
| D1 | Slurm cluster name is `farmshare` (pinned constant; `SLURM_CLUSTER_NAME`, else `scontrol show config` `ClusterName`) | `farmshare` (runtime forced to `apptainer`) | `FarmShare detected` |
| D2 | a site profile's header `// detect: cluster_name=<x>` (item 10) equals the cluster name | `<site>` (runtime from its `// runtime:` header) | `site <site> detected` |
| D3 | `sbatch` on `PATH` and the runtime below is R0–R2 (Apptainer/Singularity) or R4 (conda) | `slurm` | `sbatch on PATH` |
| D4 | otherwise | local (no executor profile) | — |

*Runtime (skipped for D1/D2). R0 first; then, on a Slurm host, R1, R2, R4 (docker cannot
run Slurm tasks, edge case (a)); on a host without Slurm, R3, R1, R2, R4:*

| # | Condition | Runtime |
|---|---|---|
| R0 | `setup.runtime` recorded in the DB directory and still usable here (binary on `PATH`; for docker `docker info` exits 0) | that runtime (`as recorded by gippa setup in <db_dir>`) |
| R1 | `apptainer` on `PATH` | `apptainer` |
| R2 | `singularity` on `PATH` | `singularity` |
| R3 | `docker` on `PATH` and `docker info` exits 0 within 10 s | `docker` |
| R4 | `conda`, `mamba` or `micromamba` on `PATH` | `conda` |
| R5 | none | `gippa run`/`preflight`: exit 2, `no container runtime or conda found: install Apptainer (no root needed on most Linux systems) or Docker; see <docs URL>`; `gippa setup`: see edge case (c) |

- **Combination.** D3 needs R0–R2 (or R4 → `slurm,conda`); with only docker (R3) the executor is
  local and a hint is printed (edge case (a)). Docker before Apptainer on a host without Slurm is
  Q-L27. R0 never selects docker for a Slurm executor.
- **Printed** to stderr before preflight: `gippa run: profile slurm,apptainer (auto: sbatch on
  PATH; apptainer as recorded by gippa setup; pass --profile to choose)`, then any hint.
  Recorded in `run_manifest.json` (`launch.profile`, `launch.profile_source: auto|flag`,
  `launch.profile_reason`).
- `--profile` given → used as is (validated as in item 3); no detection.
- **Within 15a**, R0 sees only `--db-dir`/`GIPPA_DB_DIR` (the item 6 pointer file arrives with
  15c1), and D2 is tested with a fixture site directory (the first real site lands with 15j).

**Examples**

| # | Probe | Expected |
|---|---|---|
| 1.1 | laptop: `docker` on PATH, `docker info` 0 | `docker` |
| 1.2 | workstation: `apptainer` only | `apptainer` |
| 1.3 | cluster: `sbatch`, `apptainer` | `slurm,apptainer` |
| 1.4 | SingularityCE cluster: `sbatch`, `singularity` | `slurm,singularity` |
| 1.5 | FarmShare: cluster name `farmshare` | `farmshare,apptainer`, same as today's default |
| 1.6 | `docker` and `apptainer`, no Slurm, daemon reachable | `docker`; daemon not reachable (user not in the `docker` group) → `apptainer` |
| 1.7 | `sbatch`, no container runtime, `micromamba` | `slurm,conda` |
| 1.8 | nothing on PATH, `gippa run` | exit 2 with the R5 message |
| 1.9 | `setup.runtime: docker`, docker reachable, host also has apptainer, no Slurm | `docker` (R0) |
| 1.10 | `setup.runtime: apptainer`, apptainer gone, docker reachable | R0 skipped with warning `setup used apptainer, not available here`; `docker` |
| 1.11 | FarmShare, `setup.runtime: apptainer` (setup ran on a login node) | `farmshare,apptainer` (D1), never local apptainer |
| 1.12 | `sbatch` and docker only | `docker` (local executor) + hint (a) |
| 1.13 | `sbatch`, `apptainer`, `SLURM_JOB_ID` unset | `slurm,apptainer` + hint (b) |
| 1.14 | `gippa setup` on a transfer node with no runtime | DBs downloaded; images skipped with warning (c); exit 0 if the DBs succeed |

**Edge cases**
- (a) **Slurm and Docker, no Apptainer/Singularity:** the local executor is used and a hint is
  printed: `sbatch on PATH but only docker found: tasks run on this host; install Apptainer or
  pass --profile`.
- (b) **D3 on a login node** (`SLURM_JOB_ID` unset): hint `the Nextflow head process runs on this
  host for the whole run; submit gippa run as a job (see Install › Run on a cluster)`.
- (c) **`gippa setup` with no runtime (R5)** (the item 9 setup host may be a transfer node): it
  does not exit 2. It downloads the DBs, skips the images with warning `images: no container
  runtime found here; rerun gippa setup --images --profile apptainer|singularity|docker on a host
  that has one (same --db-dir)`, records no `setup.runtime`, and exits by the DB results.
- (d) **Test isolation.** About 48 unit/integration files call `preflight`, many without
  `--profile`; under the new default they would probe the real host (GitHub runners have Docker;
  FarmShare matches D1). 15a's test PR adds an **autouse fixture** in `tests/conftest.py` that
  pins the probe to "nothing found", sets `XDG_CONFIG_HOME` and `XDG_DATA_HOME` to `tmp_path`
  subdirectories and unsets `GIPPA_DB_DIR`, so neither the host nor the item 6 pointer file leaks
  into tests. Tests that need a probe override it explicitly.
- **Inside an allocation** (`SLURM_JOB_ID` set): D3; Nextflow submits from the allocation (the
  head-job pattern). Where compute nodes may not submit, the first `sbatch` fails (exit 4, log
  named) and the docs say `--profile apptainer` (local executor in the allocation); a generic
  `slurm_local` is Q-L23.
- **Slurm without containers** (`slurm,conda`): the conda limits of item 12 apply.
- `docker info` hanging on a broken daemon, `scontrol` absent or slow: 10 s timeouts, treated as
  not found.
- Rootless Docker/Podman `docker` shims count as docker if `docker info` exits 0. Podman maps UIDs
  differently under `-u $(id -u):$(id -g)`; platforms.md lists it as untested.

**Acceptance criteria**
1. `detect_profile` table test: every executor and runtime rule and examples 1.1–1.14, with an
   injected probe (no real `PATH`, no subprocess).
2. CLI tests (CliRunner, fake `nextflow` recording argv): no `--profile` → the detected profile
   on the command line, the printed line and hints; `--profile docker` → no detection.
3. Manifest test for `launch.profile_source`/`profile_reason`; setup test for `setup.runtime`
   (and its absence in 1.14).
4. The autouse fixture of edge case (d) exists, and one test asserts that a bare `preflight`
   call in the suite sees the pinned probe.
5. FarmShare: a static test on the pinned cluster-name constant `farmshare`, plus a one-off
   check (recorded in the PR) that `gippa_run.sbatch` with its `--profile` removed resolves to
   `farmshare,apptainer`.

### 2. `gippa demo` (15b)

- **ADR needed (with 15b).** GIPPA_SPEC §5.1 documents the local executor "under Developer
  guide, not in the user guide" and §6.4 calls it developer-only. `gippa demo` exposes it to
  users, so 15b ships a short ADR (or a section of the Q-K39 ADR) recording that departure:
  the demo runs it on bundled layers only; `gippa run --executor local` stays developer-only.
- **Command.** `gippa demo [--outdir DIR] [--full] [--open]`, default `--outdir ./gippa-demo`.
  Without `--full`: runs the bundled SIM100K demo through the **local executor** (§6.4) with
  bundled precomputed layers: no Nextflow, no containers, no databases, no network. It prints
  `gippa demo: done in 48 s; report: /abs/gippa-demo/report/report.html`. `--open` opens it with
  `webbrowser.open` (no error when no browser).
- **`--full`** (after `gippa setup`): runs the same inputs through `gippa run` with the detected
  profile and the real tools, i.e. a check of the whole install (preflight, images, DBs, Nextflow).
- **Modules.** `core,mobile,bgc` from the bundled layers (gippasim writes Bakta-like JSON, operon
  pairs, promoters, Tx rates, terminators, RBS, the novel-ORF RBS layer `rbs_novel` (ADR-0024)
  and regions — every `PrecomputedConfig` field the demo modules read); `function` is left out without the
  GO/KEGG/COG resources and is included when `--full` finds them installed. The report shows a
  banner "Synthetic demo data (SIM100K, gippasim seed 20261014)"; the manifest records
  `demo: true`.
- **Data ships in the wheel** under `gippa/data/demo/`: `sim100k.fna.gz` (34 KB),
  `sim100k_peptides.tsv.gz` (68 KB), the precomputed layers (gzip, target ≤ 600 KB total), and
  `PROVENANCE.md`. **Justification:** ≈ 0.7 MB is small against the wheel and its dependencies
  (NumPy alone is ≈ 15 MB); fetching would need a hosting URL (an R6-acceptable one is a GitHub
  release asset or Zenodo), the network on first use, and one more failure mode in the very first
  command a user runs. A test caps the demo directory at 1 MiB.
- **Generation.** `scripts/docs/make_demo_layers.py [--check]` writes the layers from gippasim
  (SIM100K, the seed above), as `make_demo_peptides.py` does for the peptides; the committed
  files are byte-checked by `--check`. The genome stays the one SIM100K genome shared with the
  tests and the tutorial (`tests/data/genomes/sim100k.fna.gz`; a test asserts the wheel copy is
  byte-identical), so the demo inherits SIM100K's layout as it is.

**Examples**

| # | Input | Expected |
|---|---|---|
| 2.1 | `gippa demo` on a host with no Nextflow, Docker or DBs | exit 0 in ≤ 2 min on 2 CPUs; `gippa-demo/report/report.html` exists; path printed |
| 2.2 | `gippa demo` again, same outdir | refuses `gippa-demo exists; use --outdir or remove it` (exit 1), no overwrite |
| 2.3 | `gippa demo --full` before `gippa setup` | preflight exit 2 naming the missing DBs and `run: gippa setup` |
| 2.4 | `gippa demo --full` after setup, auto profile `docker` | `gippa run` with `--profile docker`; `status: completed` |
| 2.5 | `gippa demo --open` with no display | report path printed, exit 0 |

**Edge cases**: the SIM100K genome is two circular replicons (100 kb chromosome, 8 kb plasmid)
with genes on both strands, so the demo exercises minus-strand coordinates. **It has no gene
across the origin** (`origin_crossing_gene=False`), so the demo does not cover circular wrap;
that stays with the existing origin-crossing tests. A read-only cwd (`--outdir` needed; the error says
so); the demo under `--strict` (no warnings expected; a warning is a bug).

**Acceptance criteria**
1. Wheel content test (`scripts/release/check_dist.py`): `gippa/data/demo/` present, ≤ 1 MiB.
2. `make_demo_layers.py --check` passes in CI.
3. An integration test runs `gippa demo --outdir tmp` (local executor) and asserts the report and
   the manifest `demo: true`, and that the tables are byte-identical to a second, independent
   `gippa run --executor local` on the same bundled layers (the existing §6.4 check pattern). No
   committed sha256 listing, so table changes and lock bumps do not break it.
4. CLI tests for 2.2, 2.3 (fake preflight failure) and 2.5.

### 3. A generic `slurm` profile and executor flags (15a profile; 15e flags)

**Executor values are params written by `gippa run`, read by the profiles.** `gippa run` writes
five new keys into `OUT/work/params.json`: `slurm_queue`, `slurm_account`, `slurm_options`,
`slurm_max_jobs`, `slurm_max_time` (defaults `null` in `nextflow.config`'s `params {}`). The
profiles read them, as the container profiles already read `params.db_dir`. No `-c` file is
generated, so a profile's own values (FarmShare's QoS, excludes and R14 caps) are combined, never
replaced. Params enter no task hash unless a process script reads them, and none does.

**Profile** `workflow/conf/slurm.config`, profile `slurm`:

```groovy
process {
    executor       = 'slurm'
    queue          = params.slurm_queue ?: null
    clusterOptions = [params.slurm_options, params.slurm_account ? "--account=${params.slurm_account}" : null].findAll().join(' ')
    resourceLimits = [cpus: params.threads as int, memory: params.memory as nextflow.util.MemoryUnit,
                      time: (params.slurm_max_time ?: '48h') as nextflow.util.Duration]
}
executor {
    queueSize       = (params.slurm_max_jobs ?: 20) as int
    submitRateLimit = '10/1min'
    jobName         = { "gippa-" + task.process.tokenize(':')[-1].toLowerCase() }
}
```

- **Nothing FarmShare-specific:** no QoS, `--exclude`, MaxMemPerCPU rule, `scratch` or
  `prebuilt_envs.config`. The cluster's default partition is used unless `--queue` is given.
  `base.config`'s retry policy applies. (The implementer checks on Nextflow 26.10.0 that
  `queue = null` means "no `-p`"; else the `queue` line is written only via a ternary on a
  non-null value.)
- **`farmshare.config` reads the same params and keeps its own values:**
  `queue = params.slurm_queue ?: 'normal'`; `clusterOptions` = `'--qos=normal
  --exclude=barley-01,barley-03,rye-02'` followed by `slurm_options` and `--account=…` when set;
  `resourceLimits` keeps `Math.min(16, …)` and `[64.GB, …].min()` and takes
  `time: [48.h, max_time].min()`; `queueSize = Math.min(8, params.slurm_max_jobs ?: 8)` (R14:
  never more than 8 × 16 CPUs).

**Flags (not `--set`).** `--set` edits the GIPPA config (`resolved_config.yaml`, hashed per step
slice); executor settings must never enter it.

| Flag | Env fallback (read by `gippa`) | Param | Default |
|---|---|---|---|
| `--queue NAME` | `GIPPA_SLURM_QUEUE` | `slurm_queue` | profile's (none for `slurm`) |
| `--account NAME` | `GIPPA_SLURM_ACCOUNT` | `slurm_account` | unset |
| `--cluster-options STR` | `GIPPA_SLURM_OPTIONS` | `slurm_options` (appended to the profile's options) | unset |
| `--max-jobs N` | none | `slurm_max_jobs` | 20 (`slurm`), 8 (`farmshare`, also the cap) |
| `--max-time DUR` | none | `slurm_max_time` | 48 h (also FarmShare's cap) |
| `--threads`, `--memory` (exist) | none | `threads`, `memory`, as today | as today |

- The values are recorded in `run_manifest.json` (`launch.executor_options`; flag/env values
  only, never the rest of the environment).
- **Validation** (usage errors, Click exit 2, before preflight). The values end up in a Groovy
  config expression and in the `#SBATCH` header, so they are **rejected, not escaped**, when they
  contain a quote (`'` or `"`), `$`, a backslash, a backtick, `#`, a newline or any other control
  character: `invalid --cluster-options: quotes, $, \, `, # and control characters are not
  allowed`. Further: `--queue`/`--account` `[A-Za-z0-9_.,+-]{1,64}` (the 64-character cap is
  deliberate; a site with longer names uses `--cluster-options` or a site file);
  `--max-jobs` 1–1000; `--max-time` `^\d+(\.\d+)?\s*(s|m|h|d)$`; any of them without a
  Slurm-executor profile (`slurm`, `farmshare`, a Slurm site) → `--queue needs a Slurm profile
  (slurm, farmshare or a site profile)`.
- **Preflight:** with a Slurm-executor profile, `sbatch`, `squeue`, `scancel` on the launching
  host's `PATH`, else error `slurm executor: sbatch not found on PATH`. Queue and account are not
  queried (a wrong value fails at the first `sbatch`: exit 4 with the log).
- **Profile registry and combination rules (15a).** The accepted names are parsed from the bundled
  `nextflow.config` `profiles {}` (built-ins + sites), replacing `_PROFILES`; `test`/`stub` stay
  rejected. At most one executor profile (`standard`, `farmshare`, `farmshare_local`, `slurm`,
  sites): `--profile slurm,farmshare` → `more than one executor profile: slurm, farmshare`. At
  most one container runtime (`apptainer`, `singularity`, `docker`): today `container_runtime()`
  silently prefers apptainer; this becomes an error (Q-L16).

**Examples**

| # | Input | Expected |
|---|---|---|
| 3.1 | `--profile slurm,apptainer`, no flags | params `slurm_*` all null; `nextflow config -flat`: `process.executor = 'slurm'`, no queue, empty `clusterOptions`, `executor.queueSize = 20` |
| 3.2 | `… --queue batch --account lab42 --max-jobs 50` | params `slurm_queue = "batch"`, `slurm_account = "lab42"`, `slurm_max_jobs = 50`; merged: `process.queue = 'batch'`, `process.clusterOptions = '--account=lab42'`, `executor.queueSize = 50` |
| 3.3 | `GIPPA_SLURM_ACCOUNT=lab42` and `--account other` | `slurm_account = "other"` |
| 3.4 | `--cluster-options '--qos=long' --account lab42` | merged `clusterOptions = '--qos=long --account=lab42'` |
| 3.5 | `--queue batch --profile apptainer` | usage error, exit 2 |
| 3.6 | `--cluster-options $'--qos=x\n#SBATCH --mem=1T'` | usage error naming the control character |
| 3.7 | `slurm,apptainer`, no `sbatch` | preflight exit 2 |
| 3.8 | `farmshare,apptainer --queue bigmem` | merged `queue = 'bigmem'`; `clusterOptions` still `'--qos=normal --exclude=barley-01,barley-03,rye-02'` |
| 3.9 | `farmshare,apptainer --account lab42 --max-time 72h --max-jobs 50` | merged `clusterOptions = '--qos=normal --exclude=barley-01,barley-03,rye-02 --account=lab42'`; `resourceLimits` cpus `min(16, threads)`, memory `min(64 GB, memory)`, time `48h` (FarmShare's cap wins); `queueSize = 8` |
| 3.10 | hostile `--cluster-options "x' + \"\${System.exit(0)}\" + '"` | usage error, exit 2; no params file written, Nextflow not started |

**Edge cases**: `--max-time` below a label's `time` (clamped; 140 retries stop at the cap, as
today); `--threads` larger than a node (sbatch rejects; exit 4 names the step); `--account`
given both in `--cluster-options` and as a flag (both written, Slurm takes the last;
documented); `--resume` with other executor flags (task hashes unchanged, cache reused); a user
`--nf-config` that sets `process.clusterOptions` replaces the merged value (it is a `-c` file;
documented in item 4).

**Acceptance criteria**
1. Static test on `slurm.config` (the parity-test reader): `executor = 'slurm'`; no `qos`,
   `exclude`, `4000.MB`, `scratch`, `prebuilt_envs`; reads exactly the five `slurm_*` params;
   `jobName` starts `"gippa-"`; `submitRateLimit` present.
2. Static test on `farmshare.config`: the QoS/exclude literal, `Math.min(16`, `[64.GB` and
   `Math.min(8` are still there, and the five params are read.
3. `write_params` tests: the five keys present (null when unset), flag beats env (3.3).
4. Validation tests: 3.5, 3.6, 3.10 and one value per rejected character class.
5. CLI tests for 3.1, 3.7 and both combination errors.
6. `needs_nextflow`: `nextflow config workflow -flat -profile <p> -params-file <file>` shows the
   merged values of 3.1, 3.2, 3.4, 3.8 and 3.9.

### 4. `gippa run --nf-config FILE` (repeatable) (15e)

- Each file is passed as `-c <abs path>`, in the order given, after every built-in `-c` file.
  `gippa preflight` takes the same option to lint them.
- **Precedence, lowest first (install.md states it verbatim):** `$NXF_HOME/config`
  (`~/.nextflow/config`) < `workflow/nextflow.config` + the selected profiles (which read the
  `slurm_*` params of item 3) < `OUT/nextflow.config` (the launch directory: `gippa run` runs
  Nextflow in OUT) < `local_tools.config` < `--nf-config` files in the order given (last wins).
  The `needs_nextflow` test confirms this order on Nextflow 26.10.0 with one key set at every
  level.
  **Params are not overridable this way:** `-params-file` beats config `params`; use flags or
  `--set`.
- **Preflight rule `environment.nf_config`** (new rule id) for each file: missing / not a regular
  file / unreadable / not UTF-8 → error (exit 2), e.g. `--nf-config /x/site.config: not found`;
  a top-level `params` block → warning `… sets params; gippa run writes the params (use gippa flags
  or --set)`; `container =` or `conda =` inside a `withName` → warning `… overrides the container
  of <PROCESS>; preflight's image check does not cover it`; a `profiles` block → warning
  `… defines profiles; --profile accepts only the bundled and site profiles` (Q-L17). An existing
  `OUT/nextflow.config` → warning naming it. The scan reuses `_strip_comments`/`_blocks` of
  `src/gippa/preflight/environment.py` and never executes Groovy; syntax errors are Nextflow's
  (exit 4 with the log).
- `run_manifest.json` `launch.nf_configs = [{path, sha256}]` in order.

**Examples**: 4.1 `--nf-config site.config` → argv `… -c /abs/site.config -work-dir …`;
4.2 two files → two `-c` in order; 4.3 a missing file → exit 2 before Nextflow; 4.4
`params { threads = 4 }` → warning (exit 0; `--strict` exit 3), `params.json` keeps the CLI value;
4.5 `withName: 'BAKTA' { container = 'x' }` → warning naming BAKTA; 4.6 a relative path →
resolved against the shell's cwd.

**Edge cases**: the same file twice (passed twice); a file that `includeConfig`s another (only the
given file is scanned; documented); a 0-byte file (allowed, no warning); a directory (`not a
regular file`).

**Acceptance criteria**: scan unit tests on 4.3–4.5 and the edge cases (`tmp_path`; the
unreadable case skipped as root); CLI argv tests for 4.1, 4.2, 4.6 and the `-c` order
(`local_tools.config` < `--nf-config…`); manifest test; a `needs_nextflow` test that sets one key
(e.g. `process.queue`) in each of the five levels and asserts the winner as each top level is
removed; docs test greps the precedence list.

### 5. Nextflow and Java for pip users (15d)

- **conda path:** the Bioconda recipe's `run:` has `nextflow` and `openjdk >=17` (M11-10 Q-L10
  default; M11-10 owns the recipe). Nothing else here.
- **pip path:** `gippa setup` (zero-flag included) checks `nextflow` on `PATH`; if missing it
  downloads the **official self-contained launcher** `nextflow-26.10.0-dist` from the
  `nextflow-io/nextflow` GitHub release `v26.10.0` (R6) and installs it as
  **`<db_dir>/bin/nextflow`** (mode 0755), inside the DB directory, so it is an ordinary
  `db_sources.yaml` item (`nextflow`, `path: bin/nextflow`; `_load_sources` keeps rejecting paths
  outside the DB directory) and its `db_manifest.json` record describes a file next to it. The
  `-dist` file needs no further download at first run (the plain launcher fetches jars, which
  breaks offline nodes). `--no-nextflow` skips it.
- **Pin.** The sha256 in `db_sources.yaml` is computed once by the implementer from the official
  asset **and** cross-checked against the `digest` field that the GitHub releases API reports for
  that asset; the PR records both values and the URL.
- **Allowlist.** The asset is served from a GitHub release-asset host (the implementer verifies
  which: `release-assets.githubusercontent.com` or `objects.githubusercontent.com`, and pins only
  that one). Adding it to `REDIRECT_HOSTS["github.com"]` would also widen the UniOP archive's
  redirect allowlist, so the allowlist becomes **per item**: a `db_sources.yaml` row may name
  `redirect_hosts: [...]`, and `REDIRECT_HOSTS` stays the default for rows that do not. The
  UniOP row keeps exactly `codeload.github.com`; the `nextflow` row gets only the asset host.
- **Lookup in `gippa run`/`demo --full`:** `nextflow` on `PATH`, else `<db_dir>/bin/nextflow`;
  the manifest records which. Version advice unchanged (≥ 25.04.8; pinned 26.10.0).
- **Java:** never downloaded. `java -version` (or `$JAVA_HOME/bin/java`) is run (10 s timeout) and
  the major version parsed. In **`gippa run` and `gippa demo --full`** a missing or old Java is an
  error (exit 2): `Java 17 or newer is needed for Nextflow (found 11 at /usr/bin/java); install
  it, e.g.: micromamba create -n gippa-java -c conda-forge openjdk=21` (or `not found`). In
  **`gippa setup`** it is a warning only (the DB downloads need no Java, and transfer nodes often
  have none). On macOS, the `/usr/bin/java` stub without a JDK exits 1 (and may open a dialog):
  any nonzero exit or unparseable output counts as "not found". This supersedes M11-14 Q-K17's
  "prerequisites" for Nextflow only (Q-K42).
- `<data dir>` = `$XDG_DATA_HOME/gippa`, else `~/.local/share/gippa` (item 6's default DB
  directory lives under it).

**Examples**: 5.1 no `nextflow`, Java 21 → setup downloads, checks, installs
`<db_dir>/bin/nextflow`, records it; a later `gippa run` uses it; 5.2 checksum mismatch → item
`failed`, `pulled <a>, pinned <b>`, exit 2, nothing installed; 5.3 `gippa run` with
`openjdk version "11.0.2"` → exit 2 with the message; 5.4 `"1.8.0_292"` → major 8 → same error;
5.5 `nextflow` already on PATH → no download, the record names the found version; 5.6 `gippa
setup` with no Java → warning, DBs downloaded, exit by the DB results; 5.7 macOS stub (exit 1,
no version line) → "not found".

**Acceptance criteria**: Java version parser tests (`"21.0.1"`, `"17"`, `"1.8.0_292"`,
`"11.0.2"`, the macOS stub's exit 1, garbage → not found); setup test with a fake fetcher for
5.1/5.2/5.6 (no network); per-item redirect test (the UniOP row still refuses the asset host;
the `nextflow` row refuses `codeload.github.com`); `_load_sources` still rejects a path outside
the DB directory; run lookup-order test.

### 6. Zero-flag `gippa setup` (15c1: DB directory and plan; 15c2: defaults and licence record)

Khoa (Q-K35, 2026-10-27): "turn everything on and just document the license stuff, this mostly
targets academic use".

- **`--db-dir` becomes optional (15c1).** Default: `$GIPPA_DB_DIR`, else `<data dir>/db` (item
  5). Setup writes the chosen directory to `$XDG_CONFIG_HOME/gippa/db_dir` (else
  `~/.config/gippa/db_dir`), one line. Every reader of the DB directory — `gippa run`, `gippa
  preflight`, `gippa demo --full` **and `gippa makedb`** (which reads it for the contaminants,
  `cli.py` `makedb`) — takes `--db-dir`, else `$GIPPA_DB_DIR`, else that file, else
  `<data dir>/db`, and prints `gippa run: DB directory /x (from gippa setup)` when the file was
  used. One helper (`resolve_db_dir(flag, environ, config_home)`) serves all four.
- **New defaults (15c2)** (`SetupOptions` and the Click options; each gets a `--no-*` form):

  | Item | Today | After 15c2 |
  |---|---|---|
  | Bakta light | on | on |
  | Bakta full (84 GB) | `--bakta full` | unchanged (opt-in) |
  | geNOMAD DB | on | on |
  | antiSMASH DBs | off (`--antismash`) | **on** (`--no-antismash`) |
  | Promotech models | on | on |
  | GO, KEGG, COG, contaminants | off (`--go`, `--kegg`, `--cog`, `--contaminants`) | **on** (`--no-go`, `--no-kegg`, `--no-cog`, `--no-contaminants`) |
  | images | off (`--images`, needs `--profile`) | **on, for the detected runtime** (item 1; recorded as `setup.runtime`); `--no-images`; `--images all` unchanged; no runtime → skipped with a warning (item 1 edge case (c)) |
  | Nextflow launcher | none | on when missing (item 5); `--no-nextflow` |
  | UniOP | built under `--images` after consent | unchanged: the terms prompt / `--accept-uniop-terms`; without consent the run uses the distance model (`annotation.uniop.enabled: auto`, already the default) |
  | BacTermFinder | opt-in at run time (`--add-bactermfinder`) | unchanged, pending Khoa (Q-K43) |

- **`--minimal`** (15c2; small disks): Bakta light, Promotech models, GO, COG, KEGG, contaminants, the
  images the `core,function` modules need; no geNOMAD, no antiSMASH. Recorded as
  `setup.set: minimal`. A later `gippa run` without `--modules` on that DB directory uses
  `core,function` and prints `modules core,function (DB directory set up with --minimal; mobile
  and bgc need gippa setup without --minimal)`, recorded in the manifest (Q-L25).
- **Before any download (15c1)** setup prints the plan: one line per item with its size (from
  `db_sources.yaml` `size_gib`; the implementer fills the "not recorded" sizes of
  setup-databases.md by measuring once on FarmShare, and adds per-image sizes to
  `images.json`/`containers.config` comments), the **temporary space** a SIF build needs (the
  largest image, ≈ 9 GB for antiSMASH, in `APPTAINER_TMPDIR`, default `/tmp`), the total, the
  free space of the DB directory's, the image cache's and the temp directory's filesystems, and
  these two lines (15c2):
  - `geNOMAD (module mobile) is for non-commercial use only; for commercial use: gippa setup
    --no-genomad, gippa run --modules core,bgc,function`;
  - `KEGG is free for academic use only; for other use: gippa setup --no-kegg`.
  Free space < total × 1.1 → exit 2 before downloading, `needs 31.4 GiB, 12.0 GiB free on
  /home/u: pass --db-dir (or set GIPPA_DB_DIR) on a larger disk, or use --minimal`. On a terminal
  it asks `Proceed? [Y/n]` (default yes; `--yes` or no TTY → proceed).
- **Quotas (15c1).** `shutil.disk_usage` reports the whole filesystem (often TBs on NFS, Lustre,
  GPFS), not the user's quota, so it does not catch `$HOME` quotas. Therefore:
  - setup tries the user quota where it is cheap — `quota -w` and, on Lustre, `lfs quota -u $USER
    <mount>` — with a 10 s timeout, best effort; a parsed limit lowers the "free" figure;
  - when the DB directory is the **default** one (no `--db-dir`, no `GIPPA_DB_DIR`) and it is
    under `$HOME` on a `shared` filesystem (item 9 `fs_kind`), the plan prints `the default DB
    directory is in your home directory, which usually has a small quota on clusters; pass
    --db-dir on a project or scratch filesystem` and asks for confirmation (no TTY: proceed with
    the warning);
  - an `EDQUOT` or `ENOSPC` during a download or pull marks that item `failed` with `disk quota
    exceeded or disk full at <path>: pass --db-dir …`, removes the partial file (downloads go to
    `<name>.part` and are renamed only when complete) and goes on with the next item; exit 2.
- **Apptainer's own cache (15c1).** `apptainer pull` (`setupdb/images.py`) also stores the layer
  blobs in `~/.apptainer/cache`, roughly doubling the image bytes under `$HOME`. Setup pulls with
  `--disable-cache` unless `APPTAINER_CACHEDIR` is set by the user (then it is left alone);
  `singularity pull` likewise (`SINGULARITY_CACHEDIR`).
  **Every default is one line** in `SetupOptions` + its Click option, so flipping geNOMAD or KEGG
  back to opt-in is a one-line change and a test update.
- **Record of use (15c2).** `run_manifest.json` gains `licence_notes` (e.g. `["geNOMAD:
  non-commercial", "KEGG: academic use"]`) for the items the run actually used; the report's
  Methods footer prints the same line. (UniOP's "no licence" note from 14u2 joins the same list.)
- **Existing tests and docs that pin today's defaults** (the test-author updates them in 15c2; the
  list is from `grep` and is to be confirmed): `tests/unit/test_cli_setup.py` (`NOTHING_ELSE`),
  `test_setupdb.py`, `test_setupdb_review.py`, `test_setup_resources.py`,
  `test_setup_resources_review.py`, `test_m11_14f2_setup_images.py`,
  `test_m11_14u2_uniop_terms.py`, `test_m11_14u1_uniop_fetch.py`, `test_m11_14_install_gaps.py`,
  `test_cli_edges.py`, `test_cli_exit_codes.py`, `test_cli.py`, `test_m10_07_cli.py`,
  `test_m11_07_reference_sources.py`, `test_m11_08_docs.py`, `test_m10_12_sbatch.py`;
  `tests/integration/test_setup_resources_live.py`; and for the run default profile (item 1):
  `test_cli_run_nextflow.py`, `test_cli_run_review.py`, `test_runner_nextflow.py`,
  `test_m11_14f1_preflight_images.py`, `test_m11_14d_*.py`. Docs: `README.md` quick start,
  `docs/get-started/{install,setup-databases,upgrade}.md`. `scripts/slurm/gippa_setup.sbatch`
  keeps its explicit flags (a static test asserts the resulting set is unchanged on FarmShare).

**Examples**

| # | Input | Expected |
|---|---|---|
| 6.1 | `gippa setup` on a laptop with Docker, 200 GB free | plan printed with both licence lines; downloads Bakta light, geNOMAD, antiSMASH, Promotech, GO, KEGG, COG, contaminants; pulls images for `docker`; asks for UniOP consent; writes `~/.config/gippa/db_dir` |
| 6.2 | `gippa setup --no-genomad --no-kegg` | neither downloaded; their licence lines not printed |
| 6.3 | `gippa setup --minimal` | no geNOMAD, no antiSMASH; `setup.set: minimal`; a later `gippa run` uses `core,function` with the notice |
| 6.4 | 12 GiB free, 31 GiB needed | exit 2 before any download, message as above |
| 6.5 | `gippa setup` then `gippa run … ` with no `--db-dir` | run uses the recorded directory and says so |
| 6.6 | a run that used geNOMAD and KEGG | manifest `licence_notes` lists both; report footer shows them |
| 6.7 | `gippa setup --check` (no other flag) | verifies the default set; downloads nothing |
| 6.8 | quota hit while downloading antiSMASH (fake fetcher raises `OSError(EDQUOT)`) | `antismash` item `failed` with the `--db-dir` advice; no `*.part` left; the next items still run; exit 2 |
| 6.9 | `--profile apptainer`, `APPTAINER_CACHEDIR` unset | argv `apptainer pull --disable-cache <cache>/<name>.img docker://<ref>`; with `APPTAINER_CACHEDIR=/x` set: no `--disable-cache` |
| 6.10 | default DB dir under `$HOME` on `nfs4`, no TTY | the home-quota warning in the plan; proceeds |
| 6.11 | `gippa makedb` with no `--db-dir` after setup | uses the recorded DB directory, says so |

**Edge cases**: `--db-dir` given and different from the recorded one (the new one is recorded;
the old is not touched); `GIPPA_DB_DIR` set and the file pointing elsewhere (the variable wins,
as for every reader); a DB directory on a filesystem whose free space cannot be read (warning,
proceed); `--images all` with `--minimal` (images all, DBs minimal: allowed); an empty or
unreadable `db_dir` pointer file (ignored with a warning).

**Acceptance criteria**: 15c1: `resolve_db_dir` tests for the four readers; plan/size/temp/
total/free-space tests with an injected `disk_usage` and quota runner (`quota`/`lfs` output
fixtures, timeout); examples 6.8–6.11; the `.part` rule. 15c2: `SetupOptions()` default test (the
table); the two licence lines exactly; `--minimal` set and the
run-time module notice; pointer-file read/write tests (`tmp_path` as `XDG_CONFIG_HOME`);
manifest `licence_notes` and the report footer (snapshot text); `gippa_setup.sbatch` set unchanged.

### 7. A `singularity` profile (15f)

SingularityCE ships no `apptainer` binary, so `-profile apptainer` fails there.

```groovy
singularity {
    singularity.enabled    = true
    singularity.autoMounts = true
    singularity.cacheDir   = System.getenv('NXF_SINGULARITY_CACHEDIR') ?: System.getenv('GIPPA_ROOT')?.concat('/containers')
    singularity.runOptions = params.db_dir ? "-B ${params.db_dir}:${params.db_dir}:ro" : ''
    includeConfig 'conf/container_common.config'
    process { withName: 'BAKTA|ANTISMASH' { containerOptions = "-B ${projectDir}/assets/ps/ps:/usr/bin/ps:ro" } }
}
```

- Parity test (14p) extended to `docker`, `apptainer`, `singularity`.
- `image_cache(environ, db_dir, runtime)`: `singularity` → `NXF_SINGULARITY_CACHEDIR`, else
  `$GIPPA_ROOT/containers`, else `<db_dir>/containers`; apptainer unchanged. File names are the
  same `simpleName` (ADR-0059: `ApptainerCache` extends `SingularityCache`). `launch_env` sets
  `NXF_SINGULARITY_CACHEDIR` under `singularity` as it sets `NXF_APPTAINER_CACHEDIR` today.
- `setup --images --profile singularity`: `singularity pull <cache>/<name>.img docker://<ref>`;
  `setup --uniop --profile singularity`: `singularity build --fakeroot` (14u2 retry and failure
  rules), a `singularity` block in `local_tools.config`; `uniop_runtime()` maps it like apptainer.
- Preflight: `container_runtime()` returns `"singularity"`; the apptainer branches of
  `missing_tool_reasons` apply with that cache; hint `run: gippa setup --images --profile
  singularity`.
- **Runtime binary check (both):** `--profile apptainer` with only `singularity` on PATH → error
  `apptainer not found on PATH; this host has singularity: use --profile singularity` (auto
  detection never makes this choice); no binary → error. With a Slurm executor both are warnings
  (some clusters load the runtime only on compute nodes; Q-L18).

**Examples**: 7.1 `NXF_SINGULARITY_CACHEDIR=/c` → cache `/c`; 7.2 only `NXF_APPTAINER_CACHEDIR`
set, runtime singularity → not read (GIPPA_ROOT or `<db>/containers`); 7.3 miniprot `.img` absent →
`image quay.io/biocontainers/miniprot@sha256:… not in the singularity cache /c; run: gippa setup
--images --profile singularity`; 7.4 the fake-`singularity` pull argv; 7.5
`--profile singularity,apptainer` → usage error.

**Acceptance criteria**: parity test over three profiles (dropping the singularity DB bind fails
it); `image_cache` tests (existing apptainer tests unchanged); `missing_tool_reasons` and
runtime-binary tests with monkeypatched `shutil.which`; setup argv tests; `needs_nextflow`:
`nextflow config -flat -profile test,singularity` shows `singularity.enabled = true`.

### 8. The step image: `GIPPA_STEP` in a container (15g, 15h)

**8a. Images (15g).**
- `workflow/containers/step/Dockerfile`: `FROM python:3.12-slim@sha256:<pinned>`; `apt-get
  install --no-install-recommends procps ca-certificates` (Nextflow's task metrics need `ps`; the
  slim base has none — the same reason as the BAKTA/ANTISMASH shim; procps is GPL-2 system
  software, aggregated like the Debian base); then `pip install --require-hashes -r
  requirements.lock` (from `uv export --frozen --no-dev --no-emit-project`, hashes kept,
  committed, checked against `uv.lock` by a test); then, as a separate invocation, `pip install
  --no-deps dist/<pkg>-<version>-py3-none-any.whl`. `HOME=/tmp`; runs as any UID. A marker file
  `/opt/gippa/STEP_IMAGE` (see the R steps). Labels
  `org.opencontainers.image.{version,source,revision,licenses=MIT}`, `org.gippa.wheel.sha256`.
- **Contents rule (MIT):** only the package and its permissive Python dependencies on a Debian
  slim base with `procps` and `ca-certificates` (system software, aggregated). No bioinformatics
  tool, no R, no GPL/AGPL/non-commercial/unlicensed Python package.
  `scripts/release/check_image_licences.py` runs inside the built image and fails on a Python
  distribution licence not in `workflow/containers/step/licences.allow`.
- **R steps.** `FUNCTION_STEP` and `SENSITIVITY` start R fgsea and a container cannot start a
  sibling container. They use `<pkg>-step-r`: the fgsea image's locked R env (R GPL-2, fgsea MIT;
  M11-14 14c) plus the same wheel. R is started by `subprocess`, never imported (R5). Q-K37 asks
  whether this aggregate is acceptable; the fallback is "these two steps run on the host" (8c,
  per-step `host`).
- **The FarmShare fgsea prefix must not leak in.** `prebuilt_envs.config` (included by
  `farmshare` and `farmshare_local`) sets `env { GIPPA_FGSEA_PREFIX = <host path> }`, and
  Nextflow exports the `env` scope into every task, also inside a container. In image mode that
  host path does not exist. So both step images carry the marker file `/opt/gippa/STEP_IMAGE`,
  whose content is the image's own fgsea prefix (`/opt/conda` in step-r, empty in step). When the
  file exists, `gippa.io.fgsea_r.resolve_launcher` uses that prefix and
  **ignores `GIPPA_FGSEA_PREFIX`** (empty content → `FgseaEngineUnavailable`, which cannot happen
  because the R selector routes those two steps to step-r). A file, not an environment variable,
  because the `env` scope and Apptainer's environment pass-through can override variables but not
  files. On the host (no marker) the resolution order of M11-14 14c is unchanged.
- **Names and tags (amends ADR-0060 decision 1, in the Q-K39 ADR).** ADR-0060 names images
  `ghcr.io/<owner>/<pkg>-<tool>` with the single tag `<version>-<key>`. The step images follow the
  name rule — `ghcr.io/<owner>/<pkg>-step` and `ghcr.io/<owner>/<pkg>-step-r`, also for dev
  builds — but their tag is the package version alone, `<version>` = `gippa.__version__`: the
  image is a function of the release wheel and the committed lock, and a PyPI version is
  immutable, so the version is the content key. The wheel cannot hold the digest of an image
  that contains it, so it names the tag. ADR-0060's no-retag rule applies unchanged.
- **Pinned by digest at run time.** `release.yml` (15n) builds both images from the release wheel
  before the PyPI upload and publishes their digests as a release asset
  `step-images-<version>.json` (an official R6 source). `gippa setup` pulls by tag, checks the
  label `version == __version__`, verifies the pulled digest against that asset when it is
  reachable (else records it on first use with the warning `step image digest not verified:
  release asset unreachable`), and records `image:step`, `image:step-r` in `db_manifest.json`.
  `gippa run` then passes `params.step_image = <ref>:<version>@<recorded digest>`, so every run
  is pinned by digest (§12); the departure is narrowed to "the wheel names a tag".
- **Dev builds:** `scripts/dev/build_step_image.sh [--runtime docker|apptainer|singularity]`
  builds the wheel and both images from the checkout, tags `<pkg>-step:dev-<sha12 of the wheel>`
  and `<pkg>-step-r:dev-<sha12>`, and prints the refs for `--step-image`.

**8b. Wiring (15h).**
- `GIPPA_STEP` gains label `gippa_step` (keeps `process_single`). `container_common.config` gains
  `withLabel: 'gippa_step' { container = { params.step_image ?: null } }` and
  `withName: 'FUNCTION_STEP|SENSITIVITY' { container = { params.step_image_r ?: params.step_image ?: null } }`.
  The implementer verifies on Nextflow 26.10.0 that a `null` container runs the task on the host,
  else selects one of two config files from `gippa run`; the tests assert behaviour.
- `params.step_image`, `params.step_image_r` written only in image mode.
- **Binds are files, not `OUT`.** The task directories are `OUT/work/<xx>/<hash>` and must stay
  writable; binding `OUT` read-only would cover them (and under Apptainer the result would depend
  on bind order). So `params.step_binds` lists **files**: `params.config`
  (`OUT/resolved_config.yaml`) and each resolved `annotation.precomputed.*` path, bound read-only
  at the same path for the `gippa_step` label (`-B f:f:ro` / `-v f:f:ro`). The DB directory is
  already bound read-only by the profiles; `OUT/work` is mounted read-write by Nextflow itself;
  `publishDir` runs outside the container.
- `container_common.config`'s `beforeScript` also exports `NUMBA_CACHE_DIR=$PWD/tmp/numba`: numba's
  `cache=True` (`mapping/qgram.py`, `frameshift.py`) cannot write to the image's read-only
  site-packages and would fall back to `~/.cache` in the bound `$HOME`. `docker.envWhitelist`
  gains `NUMBA_CACHE_DIR,GIPPA_DB_DIR,GIPPA_CONTAINER_RUNTIME`.
- Manifest `launch.step_runtime` (`host`|`image`) and `launch.step_image` (ref + digest).

**8c. Which runtime: `--step-runtime auto|image|host` (default `auto`), `--step-image REF`.**
Pure `resolve_step_runtime(profile, installed, version, flag, override)`:

| Flag | Container profile | Package | Result |
|---|---|---|---|
| `host` | any | any | host |
| `image` | none | any | usage error `--step-runtime image needs a container profile` |
| `image` | yes | final release, or `--step-image` | image |
| `image` | yes | checkout / dev version, no `--step-image` | preflight error `no published step image for <version>; build one with scripts/dev/build_step_image.sh and pass --step-image` |
| `auto` | none | any | host |
| `auto` | yes | checkout (`workflow_main()` outside site-packages) or `.dev`/`+local` version | host (FarmShare development unchanged) |
| `auto` | yes | installed final release | image |
| any but `host` | yes | `--step-image REF` | image with REF |

- Preflight: image mode checks the step images like pinned tools (`step image
  ghcr.io/…/gippa-step:1.0.0 not in the apptainer cache …; run: gippa setup`); host mode with a
  non-FarmShare Slurm executor → warning `Python steps run with the host gippa (<sys.executable>);
  every compute node must see this path and the same version`. A `--step-image` whose version
  label ≠ the host's → error `step image version 1.0.1 ≠ host 1.0.0`.
- `setup` (images on) includes the step images for the installed version; skipped in a checkout
  with a note.

**Examples**: 8.1 wheel 1.0.0, `slurm,apptainer`, image cached → `params.step_image =
ghcr.io/khoa-yelo/gippa-step:1.0.0@sha256:<recorded>`, `step_runtime: image`; 8.2 checkout, `farmshare,apptainer` → host,
no warning; 8.3 wheel, image not pulled → exit 2 with `run: gippa setup`; 8.4 checkout,
`--step-runtime image` → exit 2 (`no published step image for 0.0.1 …`); 8.5 checkout,
`--step-image gippa-step:dev-abc123 --profile docker` → image; 8.6 `--profile conda --step-runtime
image` → usage error; 8.7 wheel, `--step-runtime host`, `slurm,apptainer` → host + warning;
8.8 a stub run with a local step image → every `GIPPA_STEP` trace row's `container` equals it;
8.9 `--annotation /data/x/bakta.json` → the file `/data/x/bakta.json` (and
`OUT/resolved_config.yaml`) in `step_binds`, never `/data/x` or `OUT`; 8.10 FarmShare, installed
wheel, `farmshare,apptainer` (auto → image mode): `FUNCTION_STEP` runs in step-r and uses
`/opt/conda` although the task environment has FarmShare's `GIPPA_FGSEA_PREFIX`; 8.11 a run in
host mode, then `-resume` in image mode → every Python step reruns (the `container` directive is
part of the task hash), tool tasks stay cached.

**Edge cases**: host and image dependency versions differ (tables may differ if NumPy, numba or
LightGBM numerics change; AC 6 checks identity, a difference is a finding for the statistician,
never a tolerance change, R3); `PYTHONNOUSERSITE=1` already set; Apptainer binds `$HOME` (the
image must not read it); docker `-u` with an arbitrary UID; a precomputed path that is a symlink
(the resolved target file is bound); switching between host and image mode changes every
`GIPPA_STEP` task hash, so `-resume` reruns the Python steps (documented in install.md).

**Acceptance criteria**
1. Static: `GIPPA_STEP` labels; the `gippa_step` selector reached by all three container profiles;
   the R selector names exactly `FUNCTION_STEP|SENSITIVITY`, both `GIPPA_STEP` aliases.
2. `resolve_step_runtime` table test.
3. `write_params`: `step_image*`, `step_binds` iff image mode; sorted, unique, absolute.
4. Dockerfile static test: base by digest; `procps` installed; `pip install --require-hashes -r
   requirements.lock`, then the wheel with `--no-deps`; the marker file; labels; no conda, R or
   tool installs (step-r: exactly the fgsea lock added).
4b. `resolve_launcher` test: with a marker file (in `tmp_path`, path injected) and
   `GIPPA_FGSEA_PREFIX=/host/path` set, the marker's prefix is used; without it, M11-14's order.
4c. (`needs_nextflow`, docker) a stub run under the locally built step image writes its outputs
   in the task directory (OUT/work writable) and reads `OUT/resolved_config.yaml` through the
   file bind.
5. `requirements.lock` agrees with `uv.lock` (re-export with `--no-emit-project`, hashes kept).
6. (CI 15n, advisory) one real step (`parse_bakta` on the fixture `bakta.json`) under
   `-profile test,docker` in the locally built image is byte-identical to the host run.
7. Licence check script tested on a fake site-packages (one GPL dist fails and is named).

### 9. Offline compute nodes (15i)

| Where | What | Network |
|---|---|---|
| login / transfer node | `pip install` / `conda install`; `gippa setup` (DBs, images, step images, Nextflow, UniOP build) | yes |
| launching host (head sbatch job or a login node allowed long processes) | `gippa preflight`, `gippa run --offline` | no |
| compute nodes | tasks | no |

**On a filesystem every node sees, same path:** the DB directory, the image cache, `OUT` (incl.
`OUT/work`; Nextflow's Slurm executor requires it), the Python environment holding `gippa`
(site-packages holds `workflow/` = `projectDir`; the ps shim is bound from it into
BAKTA/ANTISMASH), every input and precomputed file. `$NXF_HOME` only on the launching host.

- **`--offline`** (implied by an exported `NXF_OFFLINE=true`): `launch_env` sets
  `NXF_OFFLINE=true`; preflight turns "will be pulled at run time" into errors: docker tool images
  and (with 15h) the step images not in the local daemon (`docker image inspect` fails) → `image <ref> not pulled and --offline: run
  gippa setup where the network is`; plain `conda` processes Nextflow would create → `--offline:
  conda would create <tool>'s env at run time`.
- **Shared-FS warning** (Slurm-executor profiles only): pure `fs_kind(path, mountinfo_text)` →
  `shared` (`nfs`, `nfs4`, `lustre`, `gpfs`, `beegfs`, `cephfs`, `panfs`, `wekafs`,
  `fuse.glusterfs`), `local` (`ext2/3/4`, `xfs`, `btrfs`, `tmpfs`, `overlay`, `zfs`) or
  `unknown` (no warning); longest mount-point prefix wins; `\040` escapes decoded. `local` for
  OUT, DB dir, image cache or `workflow_main().parent` → warning `environment.shared_fs: <role>
  <path> is on a node-local filesystem (<fstype>); compute nodes will not see it`. Never an error.
- A static test asserts no module `script:` block contains `curl`, `wget`, `pip install`,
  `conda install`, `git clone` or `https://`.

**Examples**: 9.1 `slurm,apptainer --offline`, all cached → `NXF_OFFLINE=true`, exit 0;
9.2 `docker --offline`, miniprot not pulled → exit 2; 9.3 OUT on tmpfs under Slurm → warning;
9.4 OUT on nfs4 → none; 9.5 local executor, OUT on tmpfs → none; 9.6 exported `NXF_OFFLINE=true`
→ treated as `--offline`.

**Edge cases**: autofs not yet mounted (`unknown`); compute nodes that do have internet
(`--offline` optional); the head job's time limit (documented: it must outlive the run;
`-resume` after a kill).

**Acceptance criteria**: `fs_kind` table test on a mountinfo fixture (nested mounts, escaped
space); preflight tests 9.2, 9.3, 9.5 and the conda row; `launch_env` tests 9.1, 9.6; the
no-network static test.

### 10. Site configs, nf-core/configs style (15j; FarmShare move 15m)

- `workflow/conf/sites/<site>.config`, one per cluster, listed in `nextflow.config` as
  `<site> { includeConfig 'conf/sites/<site>.config' }` (Nextflow cannot glob); selectable as
  `--profile <site>,<runtime>` or auto-detected (item 1 D2). `_template.config` (leading `_`) is
  not a profile.
- **Header** (parsed by the test): `// site: <name>` = file stem, `[a-z][a-z0-9_]{1,31}`;
  `// scheduler: slurm`; `// runtime: apptainer|singularity|docker`; `// maintainer:
  @<github-handle>` (no email, R12); `// detect: cluster_name=<x>` (optional);
  `// tested: <YYYY-MM-DD> <pkg version> validation/platforms/<site>/<file>.json` or
  `// tested: no`.
- **May set:** `process.executor`, `process.resourceLimits` (required), `queue`,
  `clusterOptions`, `scratch`, `executor.queueSize`, `submitRateLimit`, `jobName` (starts
  `"gippa-"`), resource settings in `withLabel`/`withName`. **Must not set:** `container`,
  `conda`, `params`, `errorStrategy`, `maxRetries`, `*.cacheDir` (the env variables keep setup,
  preflight and Nextflow in agreement), `*.enabled`, `includeConfig` of anything but
  `'../prebuilt_envs.config'` (the FarmShare sites, below) or another file under `conf/sites/`,
  `System.getenv` of anything but `SLURM_*`, `USER`, `HOME`, `GIPPA_ROOT`, `TMPDIR`.
- **FarmShare is the first site.** 15m `git mv`s `conf/farmshare.config` and
  `conf/farmshare_local.config` to `sites/` (names unchanged; header added with `detect:` = the D1
  constant; their `includeConfig 'prebuilt_envs.config'` becomes `'../prebuilt_envs.config'`,
  since includes resolve relative to the including file) and updates the 9 tests that cite the old paths (`test_antismash_resources_pin.py`,
  `test_farmshare_mem_per_cpu.py`, `test_tool_resources_pins.py`, `test_m2_21_nightly_launch.py`,
  `test_nf_skeleton_files.py`, `test_m8_17b_wiring.py`, `test_retry_exit_140.py`,
  `test_nf_modules_lint.py`, `tests/integration/test_m7_15_retry_policy.py`) (Q-L19). Until then the validator treats them as sites in place.
- **Tests:** static (header, allowed/forbidden keys, 1:1 files ↔ profile entries, names) and
  `needs_nextflow` + the CI `package` job: `nextflow config <wf> -flat -profile <site>,<runtime>`
  exits 0 and shows `process.executor` for every site.
- **Docs:** `docs/user-guide/add-your-cluster.md`: (1) auto-detection and `--profile
  slurm,<runtime>` with flags; (2) `--nf-config` for stable settings; (3) contributing
  `sites/<site>.config` from the template with a validation run (item 11): checklist, header,
  the tests that run.

**Examples**: 10.1 complete `sites/acme.config` + entry → tests pass, `--profile acme,singularity`
accepted; 10.2 no profile entry → `site acme has no profile in nextflow.config`; 10.3
`singularity.cacheDir = '/x'` → `forbidden in site acme: cacheDir`; 10.4 `// maintainer:
someone@uni.edu` → fails; 10.5 `--profile acme,slurm` → more than one executor profile; 10.6
`detect: cluster_name=acme` and the cluster name `acme` → auto profile `acme,<runtime>`.

**Acceptance criteria**: the static tests over the template and a fixture directory of bad files
(10.2–10.4) each producing its message; the CLI help lists `slurm`, `singularity` and every site.

### 11. A second tested platform and CI (15v; CI in 15n)

**Validation run** (evidence for a `tested` row):
1. Fresh env from the **release-candidate wheel** (TestPyPI or the `release.yml` dry-run
   artefact), not a checkout.
2. The HA commands on that cluster: `gippa demo`; `gippa setup` on the login/transfer node (exit
   code and `db_manifest.json` recorded).
3. Compute nodes offline (a task's `curl -sI https://ghcr.io` fails; recorded); `gippa preflight`
   with the auto profile → 0, and the printed profile.
4. Stub run: `nextflow run <installed>/workflow/main.nf -stub -profile test,<auto profile>`.
5. `gippa run --genome <sim100k> --peptides <demo peptides> --outdir out --offline` (real tools):
   `status: completed`.
6. Tables compared with the FarmShare run of the same wheel and inputs (TE7): byte-identical, or
   each differing file listed and explained.
7. `validation/platforms/<site>/<date>.json` (site, scheduler + version, runtime + version,
   Nextflow, Java, package version, wheel sha256, commands, exit codes, job IDs, wall times,
   offline evidence, TE7 result), validated by `scripts/release/check_platform_evidence.py`
   against a JSON Schema; the platforms-page test requires each `tested` cell to name a passing
   evidence file.

**Acceptance criteria**: the checker's tests (one valid, three invalid files); platforms test
extended; the run itself is a Khoa/lead step recorded in PROGRESS with job IDs.

### 12. Conda and Bioconda (no new code; M11-10 owns the recipe)

| Tool | Pin | Package on Bioconda/conda-forge | Under `--profile conda` |
|---|---|---|---|
| bakta | 1.12.2 | yes (`bakta`) | Nextflow builds the env from `workflow/envs/bakta.yml` |
| antismash | 8.0.4 | yes | same |
| genomad | 1.12.0 | yes | same |
| miniprot | 0.18 | yes | same |
| blast (TBLASTN) | 2.16.0 | yes | same |
| transtermhp | 2.09 | yes | same |
| fgsea | 1.36.2 | yes (`bioconductor-fgsea`, R from conda-forge) | prebuilt env (started by the step, not a process) |
| ostir, promoter_calculator | 1.1.3, 1.2.4 | no: conda base + `pip` from PyPI | prebuilt env required by preflight (as today) |
| promotech, terminatornet, bactermfinder | git commits | no: GitHub sources | prebuilt env, or the GHCR image after 14j |
| uniop | git commit | no; never redistributed | local build (14u1/14u2) |

- **For M11-10 (Q-L20):** `run:` = Python deps + `nextflow` + `openjdk >=17` (its Q-L10
  default, which item 5 relies on). Not the tools: their pins conflict in one env, and AGPL
  antiSMASH / non-commercial geNOMAD would be installed for users who opt out.
- The HA conda path therefore still runs the tools from **containers** (auto runtimes R0–R3); a
  host with conda and no container runtime (runtime R4: `conda` or `slurm,conda`) passes preflight only for the six Bioconda
  tools and lists the others (Q-L24). Offline + conda needs envs created beforehand (item 9).

### 13. Licence

GIPPA stays MIT (Khoa, Q-K34, 2026-10-27); a GPL switch is not needed. The package never
imports GPL/AGPL/non-commercial/unlicensed tools (`scripts/ci/check_imports.py`, R5), there is no
all-in-one image of the GPL tools, and the step image holds only the package and permissive
Python dependencies (item 8a). Licence notices for geNOMAD and KEGG are printed and recorded
(item 6).

### 14. Later (not built): in-process bindings for permissive tools

Candidates must be MIT/BSD (or public domain) and have an official Python binding. From
`src/gippa/annotate/tools.py`: miniprot (MIT; no official binding known), BLAST+ (public domain; no
official binding), TerminatorNet (MIT; design rule 4 keeps it a separate process), fgsea (MIT; R),
codon_bias (MIT; already in-process). **None qualifies today**; a backlog note only (Q-L21).

### 15. Docs and the one-page quick start (15l)

- **`docs/get-started/quickstart.md`** (new, first in the nav; README's quick start points to
  it): the four HA commands, one short paragraph each (what it does, how long, disk), "what was
  detected" (the printed profile line and how to override it), the licence lines, and three
  "if this fails" boxes (no container runtime: install Apptainer/Docker; Java < 17; not enough
  disk: `--db-dir`/`--minimal`). ≤ 1 page (a docs test caps it at 120 lines).
  **Amended (Khoa, 2026-10-27):** the quick start is a guided path, one section per step, each
  with the exact command, what it prints, how long it takes and the disk it uses:
  1. install; 2. set up (download the databases and images; link the by-hand page
  `manual-databases.md`); 3. `gippa demo` (bundled test data, the quick "does it run" check,
  what a good report looks like); 4. example runs on the shipped example inputs (genome plus
  peptides; genome plus proteins), with a short tour of the inputs and of the outputs (gene
  table, regions, report); 5. "your own data": which input you have decides the route.
  The 120-line cap moves to the four-command block at its top ("the short version"); the page
  as a whole has no cap.
- **`docs/get-started/from-spectra.md`** (new, "From raw spectra to a GIPPA run"; linked from
  quick start step 5): for users with mass spectra but no protein or peptide table yet. One
  section per route, each with *when to use it*, the commands (database from `gippa makedb`,
  the search, `gippa convert`, `gippa run`) and its caveats:
  - DIA-NN against the called genes (Bakta CDS + contaminants): the default for annotated
    genomes; fastest, best-calibrated FDR; misses unannotated ORFs.
  - DIA-NN against six-frame ORFs (`makedb --six-frame-orfs`): novel-ORF discovery; larger
    search space, so stricter FDR and longer runs; feeds the novel-ORF module.
  - De novo sequencing (the tool and converter GIPPA supports): no database or poorly annotated
    or divergent strains; lower sensitivity; peptides placed by GIPPA's own tolerant mapper.
  - DDA engines GIPPA's `convert` already reads, listed from the code.
  A summary table (route × input × typical use × run time × what GIPPA modules it feeds), the
  pilot (PXD078459) as the worked example, and the licence of each search tool (DIA-NN's own
  terms; GIPPA never bundles or runs it unless the user supplies it).
- `install.md`: the four commands replace the five steps; "Run on a cluster" replaces "A Slurm
  (or other) executor off FarmShare" (item 9 table, flags, `--nf-config` precedence, `--offline`,
  shared FS, step image, a generic head-job sbatch example); remove "`gippa run` has no `-c`
  option"; the profile table gains `slurm`, `singularity`, sites and auto-detection.
- `setup-databases.md`: new defaults, sizes filled in, `--minimal`, licence lines, the default DB
  directory.
- `platforms.md`: rows "Linux + Slurm + Apptainer/SingularityCE (generic)" (`untested` until 15v)
  and the second cluster; the 15n CI jobs; the step image.
- `add-your-cluster.md` (item 10); `release.md` step 5 gains the step images (with 14j); CLI
  reference regenerated; changelog fragments per PR.

## CI changes for Khoa (15n, `.github/` only)

1. `release.yml`: build the step and step-r images from the release wheel, run the licence check
   inside, push `ghcr.io/<owner>/<pkg>-step:<version>` and `…/<pkg>-step-r:<version>`, publish
   `step-images-<version>.json` (both digests) as a release asset, in the `ghcr`
   environment, verify public, then upload to PyPI. Dry runs build and check, never push.
2. `docker-smoke.yml` job `stub-docker`: build the step image from the PR wheel (local tag), run
   the whole workflow `-stub -profile test,docker` with `--modules core,mobile,function`
   (antiSMASH's ~9 GB image does not fit a hosted runner) and `params.step_image` = the local tag.
   Advisory, `timeout-minutes: 30`, same path filter and weekly schedule.
3. `docker-smoke.yml` job `step-image-real`: item 8 AC 6.
4. `docker-smoke.yml` job `ha-laptop` (HA (a), advisory, weekly): fresh venv, `pip install` the
   PR wheel, `gippa demo` (needs no Docker), then `gippa preflight` on the demo inputs with no
   `--profile` (asserts `docker` is detected and the missing-DB message names `gippa setup`).
   A full `gippa setup` (≈ 30 GB) does not fit a hosted runner; HA (a) end to end is checked by
   hand on a workstation and recorded with the 15v evidence.
5. `ci.yml` `package` job: `gippa demo` from the installed wheel (local executor, no network), and
   `nextflow config "$wf" -profile <p>` for `slurm,apptainer`, `slurm,singularity` and every site.

## Edge cases (cross-cutting)

- **Circular/linear, minus strand, empty input, ties:** no coordinates, statistics or tables
  change. The demo genome (two circular replicons, genes on both strands; **no gene across the
  origin**, so circular wrap stays covered by the existing origin-crossing tests, not by the demo)
  and the byte-identity checks (item 2 AC 3, item 8 AC 6, item 11 step 6) show the launch and
  install changes did not alter them. "Empty input" here: no executor flags, an empty `--nf-config`, no
  `--db-dir`, no `--profile` (each specified above). "Ties" in auto-detection: several runtimes at
  once resolve by the fixed rule order (item 1), never by `PATH` order.
- Nextflow's profile application order is not relied on: a static test asserts that the key sets
  of the executor profiles and the runtime profiles are disjoint **except for a named allow-list**:
  the shared `container_common.config` block, and the keys that `prebuilt_envs.config` gives
  `farmshare`/`farmshare_local` (`conda.enabled`, `conda.useMicromamba`, the conda-only tools'
  `withName … conda`, `env.GIPPA_FGSEA_PREFIX`), which the `conda` runtime profile also sets to
  the same values. Each allow-listed key is set identically wherever it appears, so order cannot
  change the result; the step images ignore `env.GIPPA_FGSEA_PREFIX` through their marker file
  (item 8a).

## Estimate and PR split (each ≤ ~400 code lines; tests separate, ADR-0007)

| PR | Scope | Code | Tests | Depends on | Label |
|---|---|---|---|---|---|
| 15a | profile registry parsed from `nextflow.config`; combination rules; `slurm` profile (reading the `slurm_*` params, null defaults); `detect_profile` (executor + runtime) + default `auto` in run/preflight/setup; hints; printed choice; manifest; `setup.runtime`; the autouse test fixture (item 1 edge case (d)) | 320 | 340 | none | |
| 15b | `gippa demo` (local executor, `--full`, `--open`); demo data in the wheel; `make_demo_layers.py` (incl. `rbs_novel`); `check_dist` rule; the short ADR for the local executor | 250 | 220 | none | |
| 15c1 | `resolve_db_dir` + pointer file (run, preflight, demo, makedb); setup plan: sizes, temp space, total, free space, best-effort quota, home warning + prompt; `EDQUOT`/`ENOSPC` + `.part`; `--disable-cache` pulls | 230 | 260 | 15a | |
| 15c2 | flipped defaults + `--no-*`; images for the detected runtime, skip without runtime; `--minimal` + run notice; licence lines; `licence_notes` + report footer; update the pinned tests listed in item 6 | 170 | 200 (+ edits) | 15c1 | |
| 15d | Nextflow `-dist` into `<db_dir>/bin` + pin (cross-checked with the API digest) + per-item redirect allowlist; run lookup; Java check (setup warns, run errors) | 160 | 160 | none | R6 allowlist change |
| 15e | executor flags → `slurm_*` params, rejection rules, `farmshare.config` merge; `--nf-config` + `environment.nf_config`; `-c` order; precedence test; `sbatch` check; manifest | 300 | 320 | 15a | |
| 15f | `singularity` profile; `image_cache(…, runtime)`; `launch_env`; setup `--images/--uniop` for singularity; R(P); runtime-binary check; parity × 3 | 200 | 220 | 15a | |
| 15g | step Dockerfiles (procps, hashed lock, marker file), `requirements.lock`, `licences.allow`, `check_image_licences.py`, `build_step_image.sh`; the Q-K39 ADR | 150 | 140 | none | |
| 15h | `gippa_step` label + selectors, `resolve_step_runtime`, `--step-runtime/--step-image`, params (`<tag>@<digest>`), file binds, `NUMBA_CACHE_DIR`, envWhitelist, marker handling in `fgsea_r`, setup (release-asset digest check) + preflight for step images, manifest | 320 | 320 | 15f, 15g | `stats` review of item 8 AC 6 only |
| 15i | `--offline`, offline preflight rows (tool and step images), `fs_kind` + `environment.shared_fs`, no-network test | 150 | 180 | 15f, 15h | |
| 15j | `sites/` + template, site validator, `detect:` header wiring | 60 | 150 | 15a | |
| 15n | CI changes 1–5 | 190 | 40 | 15b, 15g, 15h | CI change, **Khoa approval** |
| 15l | docs incl. quick start; changelog fragments | 0 (≈ 400 prose) | 30 | 15a–15j | |
| 15m | FarmShare configs → `sites/` (`'../prebuilt_envs.config'`), 9 tests, docs | 10 | 40 (edits) | 15j | mechanical |
| 15v | evidence checker + schema; the validation run; platforms row | 60 | 60 | all, Q-K36 | **Khoa** (cluster access) |

Total ≈ 2,570 code (≈ 2,380 without 15n), ≈ 2,680 test lines, ≈ 400 prose.

**Order (largest barrier reduction per line first):** W1: **15a** (auto profile), **15b**
(demo), 15d, 15g. W2: 15c1 then 15c2 (zero-flag setup), 15e, 15f, 15j. W3: 15h, then 15i.
W4: 15n (Khoa), 15l.
Then the first GHCR publish with the step images (14j + 15n item 1), 15m, and 15v on
the RC wheel. 15c1, 15e, 15f and 15j all build on 15a's profile registry: merge 15a first.

## Questions

Numbering continues after Q-K35 and Q-L15 (renumber if taken elsewhere). Each has a reversible
default the PRs use.

### Closed by Khoa (2026-10-27)

- Q-K34 licence: MIT stays; GPL not needed.
- Q-K35 setup defaults: "turn everything on and just document the license stuff" (item 6).

### For Khoa

| # | Question | Default in use | PR |
|---|---|---|---|
| Q-K36 | Which second cluster for 15v, and whose account? | Stanford Sherlock (Slurm + Apptainer, offline compute nodes) if Khoa has access; else any Slurm + Apptainer/SingularityCE cluster a collaborator offers; the row stays `untested` until then | 15v |
| Q-K37 | `<pkg>-step-r` aggregates R (GPL-2) with the MIT package in one image (R started as a separate program). Acceptable? | yes, as the published fgsea image already ships R; if no, `FUNCTION_STEP`/`SENSITIVITY` run on the host in image mode | 15g, 15h |
| Q-K38 | Schedulers beyond Slurm (PBS/Torque, SGE, LSF, HTCondor, cloud batch)? | Slurm only in 1.x; others via `--nf-config` (Nextflow supports them), documented `untested`, no profile, no detection | 15a, 15j |
| Q-K39 | One ADR for: the step images reverse M11-06 Q-L8 ("no `gippa` core image (R1)"); they amend ADR-0060 decision 1 (name `<pkg>-step`/`<pkg>-step-r`, tag `<version>` instead of `<version>-<key>`); the wheel names a tag while runs pin the setup-recorded digest. Accept? | yes, one ADR with 15g | 15g |
| Q-K40 | Actions minutes for `stub-docker` (~30 min), `step-image-real`, `ha-laptop` (advisory; path-filtered PRs + weekly) | yes | 15n |
| Q-K41 | Step images pushed by `release.yml` (built from the exact release wheel, before PyPI) rather than `containers.yml` | `release.yml` | 15n |
| Q-K42 | `gippa setup` downloads the official Nextflow `-dist` launcher for pip users (reverses M11-14 Q-K17's "prerequisites" for Nextflow; Java stays a prerequisite) | yes (lead direction from Khoa's "simple installation") | 15d |
| Q-K43 | BacTermFinder: keep opt-in (`--add-bactermfinder`) or turn on with everything else? | opt-in (pending Khoa, per the lead) | 15c2 |

### For the lead

| # | Question | Default in use | PR |
|---|---|---|---|
| Q-L16 | Two container runtimes in `--profile` (silently apptainer today) | usage error | 15a |
| Q-L17 | `--profile` names defined inside a `--nf-config` file | not accepted | 15e |
| Q-L18 | Runtime binary missing on the launching host with a Slurm executor | warning (error with local executors) | 15f |
| Q-L19 | Move FarmShare configs under `sites/` | yes, one mechanical PR (15m) | 15m |
| Q-L20 | M11-10 recipe run deps | Python deps + `nextflow` + `openjdk`; no tools | M11-10 |
| Q-L21 | In-process bindings | backlog note only | none |
| Q-L22 | Executor settings: flags vs a `--set` section | flags + env fallbacks → `slurm_*` params read by the profiles, never in `resolved_config.yaml` or task hashes | 15e |
| Q-L23 | A generic `slurm_local` (local executor capped to the allocation, like `farmshare_local`) | not in 15a; `--profile apptainer` inside an allocation is documented | later |
| Q-L24 | Host with conda but no container runtime (runtime R4: `conda` or `slurm,conda`): add `gippa setup --conda-envs` to prebuild the conda-only envs? | no; preflight lists them and the message suggests installing Apptainer | later |
| Q-L25 | `--minimal` DB directory + no `--modules`: auto `core,function` with a notice, or an error | auto with notice, recorded in the manifest | 15c2 |
| Q-L26 | `gippa demo` uses the developer-only local executor (GIPPA_SPEC §5.1, §6.4) for a user command | yes, recorded in a short ADR written with 15b (item 2); `gippa run --executor local` stays developer-only | 15b |
| Q-L27 | Docker preferred over Apptainer on a host without Slurm when both are usable | yes (lead's order); R0 (`setup.runtime`) keeps setup and run consistent | 15a |

### Gaps found while writing

- **The step image cannot be pinned by digest in the wheel** (the image contains the wheel): the
  wheel names `<pkg>-step:<version>`; setup verifies the digest against the release asset and
  records it, and runs pass `<tag>@<digest>`, so runs stay pinned (§12); the Q-K39 ADR records
  the departure and the ADR-0060 amendment.
- **Steps read files outside their task directory** (`OUT/resolved_config.yaml`, precomputed
  layers); it works on the host and under Apptainer only when the paths happen to be bound.
  `step_binds` binds exactly those files read-only (never `OUT`, which holds the writable task
  directories); staging the config as an input would change every task hash (rejected).
- **FarmShare's `env.GIPPA_FGSEA_PREFIX` reaches containers** (Nextflow exports the `env` scope
  into every task); the step images' marker file makes `resolve_launcher` ignore it.
- **`FUNCTION_STEP`/`SENSITIVITY` cannot reach fgsea from inside a container** (no nested
  runtime): the step-r image (Q-K37).
- **Host vs image dependency versions** may change step bytes; CI checks one step, 15v checks a
  full run.
- **`projectDir` on compute nodes:** the ps shim is bound from the installed package, so the
  Python environment must be on the shared filesystem even in image mode. Copying the shim into
  `OUT/work` would remove that need (follow-up).
- **Home-directory quotas:** the zero-flag default DB directory is under `$HOME`, which on most
  clusters is too small for ≈ 30 GB of DBs and images, and `shutil.disk_usage` sees the
  filesystem, not the quota. Item 6 adds a best-effort quota query, a warning and confirmation
  for a default DB directory under `$HOME` on a shared filesystem, `EDQUOT`/`ENOSPC` handling
  without partial files, and `--disable-cache` pulls (Apptainer's blob cache would double the
  bytes under `$HOME`). A cluster-scratch default is not guessable.
- **DB and image sizes** "not recorded" in setup-databases.md (antiSMASH, Promotech, GO, KEGG,
  COG, contaminants; per-image sizes) must be measured before the plan line can state a total.
- `install.md` says "`gippa run` has no `-c` option" while 14u1 already passes `-c
  local_tools.config`; 15e/15l make the docs match.
- §6.3's `farmshare` bullet (`queueSize = 20`, 64 CPUs / 240 GB) differs from `farmshare.config`
  (8, 16 CPUs / 64 GB, R14); not changed here. The generic `slurm` profile uses §6.3's
  `queueSize = 20` and `submitRateLimit`.
