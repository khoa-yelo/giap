# M11-15: portable install and clusters (auto profile, demo, zero-flag setup, generic Slurm, Singularity, step image, `--nf-config`, site configs)

Milestone M11 · **pass 2, after the GIPPA rename** · depends on: M11-01 (wheel
bundles `workflow/`), M11-06 and M11-14 (GHCR images, `images.json`, `image_cache()`, `setup
--images`, 14p profile parity, 14u1/14u2 local UniOP), M11-09 (`docker-smoke.yml`), M11-09b (bare
refs, ADR-0059), M11-10 (Bioconda recipe: referenced, not duplicated), ADR-0060 (image tags, no
retagging), ADR-0064 (rename) · estimate ≈ 2,380 code lines + ≈ 2,440 test lines in 12 code/test
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
`GIAP_SPEC.md` keeps its name (the rename script's keep list). `<pkg>` means the PyPI project
name read from `pyproject.toml`, `<owner>` the GHCR
owner (`khoa-yelo`, M11-14 Q-K14).

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
| Default profile | `gippa run --profile` default `farmshare,apptainer`; `gippa preflight --profile` and `gippa setup --profile` default none. Accepted names: the hard-coded `_PROFILES` (`src/gippa/cli.py:211`); `test`/`stub` rejected. `scripts/slurm/annotate_real.sbatch` runs `gippa run` without `--profile` |
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
function `detect_profile(probe) -> Detected(profile, reason)`; `probe` is an injected record of
`PATH` lookups, `docker info` success, Slurm env variables, the Slurm cluster name and the
`db_manifest.json` setup record. The default of `gippa run --profile` becomes `auto` (shown in
`--help`); `gippa preflight` keeps "no tool check" only for `--profile none`.

**Rules, first match wins:**

| # | Condition | Profile | Reason printed |
|---|---|---|---|
| D0 | the DB directory's `db_manifest.json` records `setup.profile` P (written by `gippa setup`) and P's runtime and executor are still available | P | `as recorded by gippa setup in <db_dir>` |
| D1 | Slurm cluster name equals FarmShare's (pinned constant; read from `SLURM_CLUSTER_NAME`, else `scontrol show config` `ClusterName`) | `farmshare,apptainer` | `FarmShare detected` |
| D2 | a site profile's header `// detect: cluster_name=<x>` (item 10) equals the cluster name | `<site>,<its runtime>` | `site <site> detected` |
| D3 | `sbatch` on `PATH` and `apptainer` on `PATH` | `slurm,apptainer` | `sbatch and apptainer on PATH` |
| D4 | `sbatch` and `singularity` (no `apptainer`) | `slurm,singularity` | … |
| D5 | `docker` on `PATH` and `docker info` exits 0 (daemon reachable) | `docker` | `docker daemon reachable` |
| D6 | `apptainer` on `PATH` | `apptainer` | … |
| D7 | `singularity` on `PATH` | `singularity` | … |
| D8 | `sbatch` on `PATH`, no container runtime | `slurm,conda` | `sbatch on PATH; no container runtime` |
| D9 | `conda`, `mamba` or `micromamba` on `PATH` | `conda` | `no container runtime; conda found` |
| D10 | nothing | none → exit 2 | `no container runtime or conda found: install Apptainer (no root needed on most Linux systems) or Docker; see <docs URL>` |

- **Printed** to stderr before preflight: `gippa run: profile slurm,apptainer (auto: sbatch and
  apptainer on PATH; pass --profile to choose)`. Recorded in `run_manifest.json`
  (`launch.profile`, `launch.profile_source: auto|flag`, `launch.profile_reason`).
- `gippa setup` records the profile it used in `db_manifest.json` (`setup.profile`, schema minor),
  so D0 makes `setup` and `run` agree (the images were pulled for that runtime).
- `--profile` given → used as is (validated as in item 3); no detection.

**Examples**

| # | Probe | Expected |
|---|---|---|
| 1.1 | laptop: `docker` on PATH, `docker info` 0 | `docker` (D5) |
| 1.2 | workstation: `apptainer` only | `apptainer` (D6) |
| 1.3 | cluster login node: `sbatch`, `apptainer` | `slurm,apptainer` (D3) |
| 1.4 | SingularityCE cluster: `sbatch`, `singularity` | `slurm,singularity` (D4) |
| 1.5 | FarmShare: cluster name matches | `farmshare,apptainer` (D1), same as today's default |
| 1.6 | `docker` and `apptainer` both, no Slurm, daemon reachable | `docker` (D5); daemon not reachable (user not in the `docker` group) → `apptainer` (D6) |
| 1.7 | `sbatch`, no runtime, `micromamba` | `slurm,conda` (D8) |
| 1.8 | nothing on PATH | exit 2 with the D10 message |
| 1.9 | `db_manifest.json` says `setup.profile: docker`, docker still reachable, host also has apptainer | `docker` (D0) |
| 1.10 | `db_manifest.json` says `apptainer`, but apptainer is gone | D0 skipped; next rule; a warning `setup used apptainer, not available here` |

**Edge cases**
- **Inside an allocation** (`SLURM_JOB_ID` set): still D3/D4: Nextflow submits tasks from the
  allocation (the head-job pattern). On clusters that forbid submission from compute nodes, the
  first `sbatch` fails (exit 4, log named) and the docs say `--profile apptainer` (local executor
  within the allocation). A generic `slurm_local` (like `farmshare_local`) is Q-L23.
- **Slurm without containers** (D8): the conda path's limits apply (item 12): preflight lists the
  conda-only tools that need prebuilt envs.
- `docker info` can hang on a broken daemon: 10 s timeout → treated as not reachable.
- `scontrol` absent or slow: 10 s timeout → no cluster name; D1/D2 skipped.
- Rootless Docker/Podman's `docker` shim: counts as docker if `docker info` exits 0.

**Acceptance criteria**
1. `detect_profile` table test: every row D0–D10 and examples 1.1–1.10, with an injected probe
   (no real `PATH`, no subprocess).
2. CLI tests (CliRunner, fake `nextflow` recording argv, monkeypatched probe): no `--profile` →
   the detected profile on the `nextflow` command line and the printed line; `--profile docker`
   → no detection.
3. Manifest test for `launch.profile_source`/`profile_reason`; setup test for `setup.profile`.
4. The FarmShare nightly (unchanged scripts pass `--profile` explicitly) and
   `annotate_real.sbatch` (no `--profile`) resolve to `farmshare,apptainer` on FarmShare
   (static test on the pinned cluster-name constant + a one-off check recorded in the PR).

### 2. `gippa demo` (15b)

- **Command.** `gippa demo [--outdir DIR] [--full] [--open]`, default `--outdir ./gippa-demo`.
  Without `--full`: runs the bundled SIM100K demo through the **local executor** (§6.4) with
  bundled precomputed layers: no Nextflow, no containers, no databases, no network. It prints
  `gippa demo: done in 48 s; report: /abs/gippa-demo/report/report.html`. `--open` opens it with
  `webbrowser.open` (no error when no browser).
- **`--full`** (after `gippa setup`): runs the same inputs through `gippa run` with the detected
  profile and the real tools, i.e. a check of the whole install (preflight, images, DBs, Nextflow).
- **Modules.** `core,mobile,bgc` from the bundled layers (gippasim writes Bakta-like JSON, operon
  pairs, promoters, Tx rates, terminators, RBS and regions); `function` is left out without the
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
  files are byte-checked by `--check`. The genome and peptide files are copied from their current
  places (one source of truth: the tests and the tutorial read `gippa/data/demo/` after this PR, or
  a test asserts the copies are byte-identical, the implementer's choice).

**Examples**

| # | Input | Expected |
|---|---|---|
| 2.1 | `gippa demo` on a host with no Nextflow, Docker or DBs | exit 0 in ≤ 2 min on 2 CPUs; `gippa-demo/report/report.html` exists; path printed |
| 2.2 | `gippa demo` again, same outdir | refuses `gippa-demo exists; use --outdir or remove it` (exit 1), no overwrite |
| 2.3 | `gippa demo --full` before `gippa setup` | preflight exit 2 naming the missing DBs and `run: gippa setup` |
| 2.4 | `gippa demo --full` after setup, auto profile `docker` | `gippa run` with `--profile docker`; `status: completed` |
| 2.5 | `gippa demo --open` with no display | report path printed, exit 0 |

**Edge cases**: the SIM100K genome is two circular replicons (100 kb chromosome, 8 kb plasmid)
with genes on both strands and one spanning the chromosome origin, so the demo exercises
circular wrap and minus-strand coordinates; a read-only cwd (`--outdir` needed; the error says
so); the demo under `--strict` (no warnings expected; a warning is a bug).

**Acceptance criteria**
1. Wheel content test (`scripts/release/check_dist.py`): `gippa/data/demo/` present, ≤ 1 MiB.
2. `make_demo_layers.py --check` passes in CI.
3. An integration test runs `gippa demo --outdir tmp` (local executor) and asserts the report, the
   manifest `demo: true`, and the tables byte-identical to a committed expected listing
   (sha256 per file; TE7-style).
4. CLI tests for 2.2, 2.3 (fake preflight failure) and 2.5.

### 3. A generic `slurm` profile and executor flags (15a profile; 15e flags)

**Profile** `workflow/conf/slurm.config`, profile `slurm`:

```groovy
process {
    executor       = 'slurm'
    resourceLimits = [cpus: params.threads as int, memory: params.memory as nextflow.util.MemoryUnit, time: 48.h]
}
executor {
    queueSize       = 20
    submitRateLimit = '10/1min'
    jobName         = { "gippa-" + task.process.tokenize(':')[-1].toLowerCase() }
}
```

- **Nothing FarmShare-specific:** no `queue`, `clusterOptions`, QoS, `--exclude`, MaxMemPerCPU
  rule, `scratch` or `prebuilt_envs.config`; those stay in `farmshare` (unchanged). The cluster's
  default partition is used unless `--queue` is given. `base.config`'s retry policy applies.

**Flags (not `--set`).** `--set` edits the GIPPA config (`resolved_config.yaml`, hashed per step
slice); executor settings must never enter it or a task hash.

| Flag | Env fallback (read by `gippa`) | Nextflow setting | Default |
|---|---|---|---|
| `--queue NAME` | `GIPPA_SLURM_QUEUE` | `process.queue` | unset |
| `--account NAME` | `GIPPA_SLURM_ACCOUNT` | `process.clusterOptions` gains `--account=NAME` | unset |
| `--cluster-options STR` | `GIPPA_SLURM_OPTIONS` | `process.clusterOptions` (account appended) | unset |
| `--max-jobs N` | none | `executor.queueSize` | 20 |
| `--max-time DUR` | none | `process.resourceLimits.time` (the whole map rewritten, cpus/memory kept) | 48 h |
| `--threads`, `--memory` (exist) | none | per-task caps via `params`, as today | as today |

- `gippa run` renders them with a pure `render_executor_config(opts) -> str | None` into
  `OUT/work/executor.config`, passed as `-c` after `local_tools.config` and before `--nf-config`
  files; no file without options. Values recorded in `run_manifest.json`
  (`launch.executor_options`; flag values only, never the environment).
- **Validation** (usage errors, Click exit 2, before preflight): `--queue`/`--account`
  `[A-Za-z0-9_.,+-]{1,64}`; `--cluster-options` without control characters or `#` (Nextflow
  writes it into the `#SBATCH` header); `--max-jobs` 1–1000; `--max-time`
  `^\d+(\.\d+)?\s*(s|m|h|d)$`; any of them without a Slurm-executor profile (`slurm`,
  `farmshare`, a Slurm site) → `--queue needs a Slurm profile (slurm, farmshare or a site
  profile)`.
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
| 3.1 | `--profile slurm,apptainer`, no flags | no `executor.config` |
| 3.2 | `… --queue batch --account lab42 --max-jobs 50` | `process.queue = 'batch'`, `process.clusterOptions = '--account=lab42'`, `executor.queueSize = 50`; `-c OUT/work/executor.config` after `local_tools.config` |
| 3.3 | `GIPPA_SLURM_ACCOUNT=lab42` and `--account other` | `--account=other` |
| 3.4 | `--cluster-options '--qos=long' --account lab42` | `'--qos=long --account=lab42'` |
| 3.5 | `--queue batch --profile apptainer` | usage error, exit 2 |
| 3.6 | `--cluster-options $'--qos=x\n#SBATCH --mem=1T'` | usage error naming the control character |
| 3.7 | `slurm,apptainer`, no `sbatch` | preflight exit 2 |
| 3.8 | `farmshare,apptainer --queue bigmem` | allowed; QoS/exclude unchanged (only an explicit `--cluster-options` replaces them; documented) |

**Edge cases**: `--max-time` below a label's `time` (clamped; 140 retries stop at the cap, as
today); `--threads` larger than a node (sbatch rejects; exit 4 names the step); `--account` given
both in `--cluster-options` and as a flag (both written, Slurm takes the last; documented);
`--resume` with other executor flags (task hashes unchanged, cache reused).

**Acceptance criteria**
1. Static test on `slurm.config` (the parity-test reader): `executor = 'slurm'`; none of `queue`,
   `clusterOptions`, `scratch`, `qos`, `exclude`, `4000.MB`, `prebuilt_envs`; `jobName` starts
   `"gippa-"`; `queueSize`, `submitRateLimit` present.
2. `render_executor_config` byte tests for 3.2–3.4 and `None` for no options.
3. `launch_command` `-c` order: `local_tools.config` < `executor.config` < `--nf-config…`.
4. CLI tests for 3.1, 3.5, 3.6, 3.7 and both combination errors.
5. `needs_nextflow`: `nextflow config workflow -flat -profile slurm,apptainer -c <3.2 file>`
   shows `process.executor = 'slurm'`, `process.queue = 'batch'`, `executor.queueSize = 50`.

### 4. `gippa run --nf-config FILE` (repeatable) (15e)

- Each file is passed as `-c <abs path>`, in the order given, after every built-in `-c` file.
  `gippa preflight` takes the same option to lint them.
- **Precedence (install.md states it verbatim):** `workflow/nextflow.config` + selected profiles <
  `$NXF_HOME/config` and `OUT/nextflow.config` (Nextflow's own lookup; `gippa run` runs Nextflow in
  OUT) < `local_tools.config` < `executor.config` < `--nf-config` files in order (last wins).
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
unreadable case skipped as root); CLI argv tests for 4.1, 4.2, 4.6; manifest test; docs test greps
the precedence list.

### 5. Nextflow and Java for pip users (15d)

- **conda path:** the Bioconda recipe's `run:` has `nextflow` and `openjdk >=17` (M11-10 Q-L10
  default; M11-10 owns the recipe). Nothing else here.
- **pip path:** `gippa setup` (zero-flag included) checks `nextflow` on `PATH`; if missing it
  downloads the **official self-contained launcher** `nextflow-26.10.0-dist` from the
  `nextflow-io/nextflow` GitHub release `v26.10.0` (R6), checks its sha256 against a pin in
  `src/gippa/data/db_sources.yaml` (item `nextflow`, computed once by the implementer from the
  official asset and recorded with the URL in the PR), installs it as
  `<data dir>/bin/nextflow` (mode 0755) and records `tool:nextflow` in `db_manifest.json`. The
  `-dist` file needs no further download at first run (the plain launcher would fetch jars, which
  breaks offline nodes). `REDIRECT_HOSTS` gains `github.com → objects.githubusercontent.com,
  release-assets.githubusercontent.com` (the release-asset hosts; implementer verifies which one
  GitHub uses and pins only that). `--no-nextflow` skips it.
- **Lookup in `gippa run`/`demo --full`:** `nextflow` on `PATH`, else `<data dir>/bin/nextflow`;
  the manifest records which. Version advice unchanged (≥ 25.04.8; pinned 26.10.0).
- **Java:** never downloaded. `gippa setup` and `gippa run` run `java -version` (or
  `$JAVA_HOME/bin/java`), parse the major version, and fail with exit 2: `Java 17 or newer is
  needed for Nextflow (found 11 at /usr/bin/java); install it, e.g.: micromamba create -n gippa-java
  -c conda-forge openjdk=21` (or `not found`). This supersedes M11-14 Q-K17's "prerequisites" for
  Nextflow only (Q-K42).
- `<data dir>` = `$XDG_DATA_HOME/gippa`, else `~/.local/share/gippa` (item 6 uses the same).

**Examples**: 5.1 no `nextflow`, Java 21 → setup downloads, checks, installs, records; a later
`gippa run` uses it; 5.2 checksum mismatch → item `failed`, `pulled <a>, pinned <b>`, exit 2,
nothing installed; 5.3 `java -version` prints `openjdk version "11.0.2"` → exit 2 with the
message; 5.4 Java prints `"1.8.0_292"` → major 8 → same error; 5.5 `nextflow` already on PATH →
no download, `tool:nextflow` records the found version.

**Acceptance criteria**: Java version parser tests (`"21.0.1"`, `"17"`, `"1.8.0_292"`,
`"11.0.2"`, garbage → not found); setup test with a fake fetcher for 5.1/5.2 (no network);
`REDIRECT_HOSTS` test; run lookup-order test.

### 6. Zero-flag `gippa setup` (15c)

Khoa (Q-K35, 2026-10-27): "turn everything on and just document the license stuff, this mostly
targets academic use".

- **`--db-dir` becomes optional.** Default: `$GIPPA_DB_DIR`, else `<data dir>/db` (item 5).
  Setup writes the chosen directory to `$XDG_CONFIG_HOME/gippa/db_dir` (else
  `~/.config/gippa/db_dir`), one line. `gippa run`, `gippa preflight` and `gippa demo --full` take
  `--db-dir`, else `$GIPPA_DB_DIR`, else that file, else `<data dir>/db`, and print `gippa run:
  DB directory /x (from gippa setup)` when the file was used.
- **New defaults** (`SetupOptions` and the Click options; each gets a `--no-*` form):

  | Item | Today | After 15c |
  |---|---|---|
  | Bakta light | on | on |
  | Bakta full (84 GB) | `--bakta full` | unchanged (opt-in) |
  | geNOMAD DB | on | on |
  | antiSMASH DBs | off (`--antismash`) | **on** (`--no-antismash`) |
  | Promotech models | on | on |
  | GO, KEGG, COG, contaminants | off (`--go`, `--kegg`, `--cog`, `--contaminants`) | **on** (`--no-go`, `--no-kegg`, `--no-cog`, `--no-contaminants`) |
  | images | off (`--images`, needs `--profile`) | **on, for the detected profile** (item 1); `--no-images`; `--images all` unchanged |
  | Nextflow launcher | none | on when missing (item 5); `--no-nextflow` |
  | UniOP | built under `--images` after consent | unchanged: the terms prompt / `--accept-uniop-terms`; without consent the run uses the distance model (`annotation.uniop.enabled: auto`, already the default) |
  | BacTermFinder | opt-in at run time (`--add-bactermfinder`) | unchanged, pending Khoa (Q-K43) |

- **`--minimal`** (small disks): Bakta light, Promotech models, GO, COG, KEGG, contaminants, the
  images the `core,function` modules need; no geNOMAD, no antiSMASH. Recorded as
  `setup.set: minimal`. A later `gippa run` without `--modules` on that DB directory uses
  `core,function` and prints `modules core,function (DB directory set up with --minimal; mobile
  and bgc need gippa setup without --minimal)`, recorded in the manifest (Q-L25).
- **Before any download** setup prints the plan: one line per item with its size (from
  `db_sources.yaml` `size_gib`; the implementer fills the "not recorded" sizes of
  setup-databases.md by measuring once on FarmShare, and adds per-image sizes to
  `images.json`/`containers.config` comments), the total, the free space of the DB directory's
  and the image cache's filesystems, and these two lines:
  - `geNOMAD (module mobile) is for non-commercial use only; for commercial use: gippa setup
    --no-genomad, gippa run --modules core,bgc,function`;
  - `KEGG is free for academic use only; for other use: gippa setup --no-kegg`.
  Free space < total × 1.1 → exit 2 before downloading, `needs 31.4 GiB, 12.0 GiB free on
  /home/u: pass --db-dir (or set GIPPA_DB_DIR) on a larger disk, or use --minimal`. On a terminal
  it asks `Proceed? [Y/n]` (default yes; `--yes` or no TTY → proceed).
  **Every default is one line** in `SetupOptions` + its Click option, so flipping geNOMAD or KEGG
  back to opt-in is a one-line change and a test update.
- **Record of use.** `run_manifest.json` gains `licence_notes` (e.g. `["geNOMAD:
  non-commercial", "KEGG: academic use"]`) for the items the run actually used; the report's
  Methods footer prints the same line. (UniOP's "no licence" note from 14u2 joins the same list.)
- **Existing tests and docs that pin today's defaults** (the test-author updates them in 15c; the
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

**Edge cases**: `--db-dir` given and different from the recorded one (the new one is recorded;
the old is not touched); `GIPPA_DB_DIR` set and the file pointing elsewhere (the variable wins,
as for every reader); a DB directory on a filesystem whose free space cannot be read (warning,
proceed); `--images all` with `--minimal` (images all, DBs minimal: allowed); an empty or
unreadable `db_dir` pointer file (ignored with a warning).

**Acceptance criteria**: `SetupOptions()` default test (the table); plan/size/total/free-space
tests with an injected `disk_usage`; the two licence lines exactly; `--minimal` set and the
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
- `workflow/containers/step/Dockerfile`: `FROM python:3.12-slim@sha256:<pinned>`; installs the
  locked runtime dependencies (`requirements.lock` from `uv export --frozen --no-dev
  --no-hashes`, committed, checked against `uv.lock` by a test), then `pip install --no-deps
  dist/<pkg>-<version>-py3-none-any.whl`; `HOME=/tmp`; runs as any UID. Labels
  `org.opencontainers.image.{version,source,revision,licenses=MIT}`, `org.gippa.wheel.sha256`.
- **Contents rule (MIT):** only the package and its permissive Python dependencies on a Debian
  slim base (system software, aggregated). No bioinformatics tool, no R, no GPL/AGPL/
  non-commercial/unlicensed Python package. `scripts/release/check_image_licences.py` runs inside
  the built image and fails on a distribution licence not in
  `workflow/containers/step/licences.allow`.
- **R steps.** `FUNCTION_STEP` and `SENSITIVITY` start R fgsea and a container cannot start a
  sibling container. They use `<pkg>-step-r`: the fgsea image's locked R env (R GPL-2, fgsea MIT;
  M11-14 14c) plus the same wheel, with `GIPPA_FGSEA_PREFIX=/opt/conda`, so `resolve_launcher`
  takes the prefix and never the nested-runtime branch. R is started by `subprocess`, never
  imported (R5). Q-K37 asks whether this aggregate is acceptable; the fallback is "these two
  steps run on the host" (8c, per-step `host`).
- **Tags:** `ghcr.io/<owner>/<pkg>:<version>` and `…/<pkg>-step-r:<version>`, `<version>` =
  `gippa.__version__`, built from **the release wheel** by `release.yml` before the PyPI upload
  (15n). The wheel cannot hold the digest of an image that contains it, so it names the tag;
  `setup` pulls it, checks the label `version == __version__`, and records the digest
  (`image:step`, `image:step-r`) in `db_manifest.json`; preflight checks the cached file against
  that record. ADR-0060's no-retag rule covers these tags (Gaps).
- **Dev builds:** `scripts/dev/build_step_image.sh [--runtime docker|apptainer|singularity]`
  builds the wheel and image from the checkout, tag `<pkg>-step:dev-<sha12 of the wheel>`, and
  prints the ref for `--step-image`.

**8b. Wiring (15h).**
- `GIPPA_STEP` gains label `gippa_step` (keeps `process_single`). `container_common.config` gains
  `withLabel: 'gippa_step' { container = { params.step_image ?: null } }` and
  `withName: 'FUNCTION_STEP|SENSITIVITY' { container = { params.step_image_r ?: params.step_image ?: null } }`.
  The implementer verifies on Nextflow 26.10.0 that a `null` container runs the task on the host,
  else selects one of two config files from `gippa run`; the tests assert behaviour.
- `params.step_image`, `params.step_image_r` written only in image mode.
- **Binds:** `params.step_binds` = sorted unique absolute directories: `OUT`, `db_dir`, the
  resolved parent of each `annotation.precomputed.*` path; bound read-only for the `gippa_step`
  label in all three container profiles (`publishDir` runs outside the container).
  `docker.envWhitelist` gains `GIPPA_DB_DIR,GIPPA_CONTAINER_RUNTIME`.
- Manifest `launch.step_runtime` (`host`|`image`) and `launch.step_image` (ref + recorded digest).

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
  ghcr.io/…/gippa:1.0.0 not in the apptainer cache …; run: gippa setup`); host mode with a
  non-FarmShare Slurm executor → warning `Python steps run with the host gippa (<sys.executable>);
  every compute node must see this path and the same version`. A `--step-image` whose version
  label ≠ the host's → error `step image version 1.0.1 ≠ host 1.0.0`.
- `setup` (images on) includes the step images for the installed version; skipped in a checkout
  with a note.

**Examples**: 8.1 wheel 1.0.0, `slurm,apptainer`, image cached → `params.step_image =
ghcr.io/khoa-yelo/gippa:1.0.0`, `step_runtime: image`; 8.2 checkout, `farmshare,apptainer` → host,
no warning; 8.3 wheel, image not pulled → exit 2 with `run: gippa setup`; 8.4 checkout,
`--step-runtime image` → exit 2 (`no published step image for 0.0.1 …`); 8.5 checkout,
`--step-image gippa-step:dev-abc123 --profile docker` → image; 8.6 `--profile conda --step-runtime
image` → usage error; 8.7 wheel, `--step-runtime host`, `slurm,apptainer` → host + warning;
8.8 a stub run with a local step image → every `GIPPA_STEP` trace row's `container` equals it;
8.9 `--annotation /data/x/bakta.json` → `/data/x` in `step_binds`.

**Edge cases**: host and image dependency versions differ (tables may differ if NumPy, numba or
LightGBM numerics change; AC 6 checks identity, a difference is a finding for the statistician,
never a tolerance change, R3); `PYTHONNOUSERSITE=1` already set; Apptainer binds `$HOME` (the
image must not read it); docker `-u` with an arbitrary UID; a precomputed path that is a symlink
(the resolved parent is bound).

**Acceptance criteria**
1. Static: `GIPPA_STEP` labels; the `gippa_step` selector reached by all three container profiles;
   the R selector names exactly `FUNCTION_STEP|SENSITIVITY`, both `GIPPA_STEP` aliases.
2. `resolve_step_runtime` table test.
3. `write_params`: `step_image*`, `step_binds` iff image mode; sorted, unique, absolute.
4. Dockerfile static test (base by digest; `requirements.lock` then `--no-deps` wheel; labels; no
   conda, R or tool installs).
5. `requirements.lock` agrees with `uv.lock`.
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
  `NXF_OFFLINE=true`; preflight turns "will be pulled at run time" into errors: docker images not
  in the local daemon (`docker image inspect` fails) → `image <ref> not pulled and --offline: run
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
  preflight and Nextflow in agreement), `*.enabled`, `includeConfig` outside `conf/`,
  `System.getenv` of anything but `SLURM_*`, `USER`, `HOME`, `GIPPA_ROOT`, `TMPDIR`.
- **FarmShare is the first site.** 15m `git mv`s `conf/farmshare.config` and
  `conf/farmshare_local.config` to `sites/` (names unchanged; header added with `detect:` = the D1
  constant) and updates the 9 tests that cite the old paths (`test_antismash_resources_pin.py`,
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
- The HA conda path therefore still runs the tools from **containers** (auto profile D3–D7); a
  host with conda and no container runtime (D8/D9) passes preflight only for the six Bioconda
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
   inside, push `ghcr.io/<owner>/<pkg>:<version>` and `…-step-r:<version>` in the `ghcr`
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
  change. The demo genome (two circular replicons, both strands, a gene across the origin) and
  the byte-identity checks (item 2 AC 3, item 8 AC 6, item 11 step 6) show the launch and install
  changes did not alter them. "Empty input" here: no executor flags, an empty `--nf-config`, no
  `--db-dir`, no `--profile` (each specified above). "Ties" in auto-detection: several runtimes at
  once resolve by the fixed rule order (item 1), never by `PATH` order.
- Nextflow's profile application order is not relied on: no two combinable profiles set the same
  key (a static test asserts the key sets of executor and runtime profiles are disjoint, except
  the shared `container_common.config` block).

## Estimate and PR split (each ≤ ~400 code lines; tests separate, ADR-0007)

| PR | Scope | Code | Tests | Depends on | Label |
|---|---|---|---|---|---|
| 15a | profile registry parsed from `nextflow.config`; combination rules; `slurm` profile; `detect_profile` + default `auto` in run/preflight/setup; printed choice; manifest; `setup.profile` | 300 | 320 | none | |
| 15b | `gippa demo` (local executor, `--full`, `--open`); demo data in the wheel; `make_demo_layers.py`; `check_dist` rule | 250 | 220 | none | |
| 15c | zero-flag setup: optional `--db-dir` + pointer file, flipped defaults, `--no-*`, `--minimal` + run notice, plan/sizes/free space/prompt, licence lines, `licence_notes` + report footer; update the pinned tests listed in item 6 | 280 | 300 (+ edits) | 15a (detected profile for images) | |
| 15d | Nextflow `-dist` fetch + pin + `REDIRECT_HOSTS`; run lookup; Java check | 150 | 150 | none | R6 allowlist change |
| 15e | executor flags + `render_executor_config`; `--nf-config` + `environment.nf_config`; `-c` order; `sbatch` check; manifest | 300 | 300 | 15a | |
| 15f | `singularity` profile; `image_cache(…, runtime)`; `launch_env`; setup `--images/--uniop` for singularity; R(P); runtime-binary check; parity × 3 | 200 | 220 | 15a | |
| 15g | step Dockerfiles, `requirements.lock`, `licences.allow`, `check_image_licences.py`, `build_step_image.sh` | 140 | 130 | none | |
| 15h | `gippa_step` label + selectors, `resolve_step_runtime`, `--step-runtime/--step-image`, params, binds, envWhitelist, setup + preflight for step images, manifest | 300 | 300 | 15f, 15g | `stats` review of item 8 AC 6 only |
| 15i | `--offline`, offline preflight rows, `fs_kind` + `environment.shared_fs`, no-network test | 150 | 180 | 15f | |
| 15j | `sites/` + template, site validator, `detect:` header wiring | 60 | 150 | 15a | |
| 15n | CI changes 1–5 | 180 | 40 | 15b, 15g, 15h | CI change, **Khoa approval** |
| 15l | docs incl. quick start; changelog fragments | 0 (≈ 400 prose) | 30 | 15a–15j | |
| 15m | FarmShare configs → `sites/`, 9 tests, docs | 10 | 40 (edits) | 15j | mechanical |
| 15v | evidence checker + schema; the validation run; platforms row | 60 | 60 | all, Q-K36 | **Khoa** (cluster access) |

Total ≈ 2,380 code (≈ 2,200 without 15n), ≈ 2,440 test lines, ≈ 400 prose.

**Order (largest barrier reduction per line first):** W1: **15a** (auto profile), **15b**
(demo), 15d, 15g. W2: 15c (zero-flag setup), 15e, 15f, 15j. W3: 15h, 15i. W4: 15n (Khoa), 15l.
Then the first GHCR publish with the step images (14j + 15n item 1), 15m, and 15v on
the RC wheel. 15c, 15e, 15f and 15j all build on 15a's profile registry: merge 15a first.

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
| Q-K39 | A package image on GHCR reverses M11-06 Q-L8 ("no `gippa` core image (R1)"); and it is pinned by tag, not digest, in the wheel. Accept by ADR? | yes, both recorded in one ADR | 15g |
| Q-K40 | Actions minutes for `stub-docker` (~30 min), `step-image-real`, `ha-laptop` (advisory; path-filtered PRs + weekly) | yes | 15n |
| Q-K41 | Step images pushed by `release.yml` (built from the exact release wheel, before PyPI) rather than `containers.yml` | `release.yml` | 15n |
| Q-K42 | `gippa setup` downloads the official Nextflow `-dist` launcher for pip users (reverses M11-14 Q-K17's "prerequisites" for Nextflow; Java stays a prerequisite) | yes (lead direction from Khoa's "simple installation") | 15d |
| Q-K43 | BacTermFinder: keep opt-in (`--add-bactermfinder`) or turn on with everything else? | opt-in (pending Khoa, per the lead) | 15c |

### For the lead

| # | Question | Default in use | PR |
|---|---|---|---|
| Q-L16 | Two container runtimes in `--profile` (silently apptainer today) | usage error | 15a |
| Q-L17 | `--profile` names defined inside a `--nf-config` file | not accepted | 15e |
| Q-L18 | Runtime binary missing on the launching host with a Slurm executor | warning (error with local executors) | 15f |
| Q-L19 | Move FarmShare configs under `sites/` | yes, one mechanical PR (15m) | 15m |
| Q-L20 | M11-10 recipe run deps | Python deps + `nextflow` + `openjdk`; no tools | M11-10 |
| Q-L21 | In-process bindings | backlog note only | none |
| Q-L22 | Executor settings: flags vs a `--set` section | flags + env fallbacks → `executor.config`, never in `resolved_config.yaml` or task hashes | 15e |
| Q-L23 | A generic `slurm_local` (local executor capped to the allocation, like `farmshare_local`) | not in 15a; `--profile apptainer` inside an allocation is documented | later |
| Q-L24 | Host with conda but no container runtime (D8/D9): add `gippa setup --conda-envs` to prebuild the conda-only envs? | no; preflight lists them and the message suggests installing Apptainer | later |
| Q-L25 | `--minimal` DB directory + no `--modules`: auto `core,function` with a notice, or an error | auto with notice, recorded in the manifest | 15c |
| Q-L26 | `gippa demo` uses the developer-only local executor (§6.4) for a user command | yes: no Nextflow/containers/DBs needed for the first command; `--full` runs Nextflow | 15b |
| Q-L27 | Docker preferred over Apptainer when both are usable (D5 before D6) | yes (lead's order); D0 keeps setup and run consistent | 15a |

### Gaps found while writing

- **The step image cannot be pinned by digest in the wheel** (the image contains the wheel):
  version tag + no-retag policy + digest recorded by setup and checked by preflight; departs from
  §3/§14.6 "by digest" for this image (Q-K39 ADR).
- **Steps read files outside their task directory** (`OUT/resolved_config.yaml`, precomputed
  layers); it works on the host and under Apptainer only when the paths happen to be bound.
  `step_binds` fixes it; staging the config as an input would change every task hash (rejected).
- **`FUNCTION_STEP`/`SENSITIVITY` cannot reach fgsea from inside a container** (no nested
  runtime): the step-r image (Q-K37).
- **Host vs image dependency versions** may change step bytes; CI checks one step, 15v checks a
  full run.
- **`projectDir` on compute nodes:** the ps shim is bound from the installed package, so the
  Python environment must be on the shared filesystem even in image mode. Copying the shim into
  `OUT/work` would remove that need (follow-up).
- **Home-directory quotas:** the zero-flag default DB directory is under `$HOME`, which on most
  clusters is too small for ≈ 30 GB of DBs and images; the free-space check stops early and the
  message names `--db-dir`/`GIPPA_DB_DIR`/`--minimal`. A cluster-scratch default is not
  guessable.
- **DB and image sizes** "not recorded" in setup-databases.md (antiSMASH, Promotech, GO, KEGG,
  COG, contaminants; per-image sizes) must be measured before the plan line can state a total.
- `install.md` says "`gippa run` has no `-c` option" while 14u1 already passes `-c
  local_tools.config`; 15e/15l make the docs match.
- §6.3's `farmshare` bullet (`queueSize = 20`, 64 CPUs / 240 GB) differs from `farmshare.config`
  (8, 16 CPUs / 64 GB, R14); not changed here. The generic `slurm` profile uses §6.3's
  `queueSize = 20` and `submitRateLimit`.
