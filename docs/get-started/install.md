# Install

GIPPA has two parts:

- the Python package, which provides the `gippa` command and bundles the Nextflow workflow, its
  container pins (`workflow/conf/containers.config`, `workflow/conf/images.json`) and the database
  pins (`src/gippa/data/db_sources.yaml`);
- the tools it runs as separate programs (Bakta, miniprot, geNOMAD, antiSMASH and others). Each
  tool runs from a container image pinned by digest, or, on FarmShare, from a prebuilt micromamba
  environment.

`gippa preflight`, `gippa convert`, `gippa compare` and `gippa report` need only the Python package.
`gippa setup` needs no Nextflow (it installs the pinned launcher when none is on `PATH`). `gippa
run` needs Nextflow, Java 17 or newer, a container runtime, the tool images and the databases.

## Tested platforms

The full pipeline is tested on Stanford FarmShare. Every other row is as stated:

--8<-- "docs/user-guide/platforms.md:platforms"

The cell words (tested, smoke-tested, untested) and the CI jobs are explained on
[Platforms](../user-guide/platforms.md).

!!! warning "Before the first image publish"
    Four core tools, Promotech, Promoter Calculator, TerminatorNet and OSTIR, have Dockerfiles
    under `workflow/containers/` but no published image yet: the first GHCR publish is a manual
    release step ([Releasing › step 5](../release.md#5-ghcr-publishing-and-package-visibility-khoa)).
    Until it is done, `gippa run` under `--profile docker` or `--profile apptainer` alone stops in
    preflight with `conda-only tool under a container profile without conda: <tool>` (exit code 2)
    for each of them. Two opt-in tools are missing the same way: BacTermFinder (with
    `--add-bactermfinder`), and fgsea (module `function` with `--abundance`; `images.json` pins
    no image for it yet). On FarmShare, use the [checkout path](#farmshare-from-a-checkout). Once a
    release pins those images, the five steps below need nothing else.

## Prerequisites

| Requirement | Version | Where it is pinned |
|---|---|---|
| Python | Python ≥ 3.11 | `requires-python = ">=3.11"` in `pyproject.toml` |
| Nextflow | 26.10.0 is the pinned, tested version; at least 25.04.8 (advice: GIPPA does not check the version) | `scripts/bootstrap_farmshare.sh`, [ADR-0002](../adr/ADR-0002-farmshare-env-pins.md) |
| Java | 17 or newer, which Nextflow needs; the bootstrap installs OpenJDK 21 | `scripts/bootstrap_farmshare.sh` |
| Docker or Apptainer | the `docker` CLI with access to a daemon, or `apptainer` | see [Platforms](../user-guide/platforms.md) |

**Nextflow.** When `nextflow` is not on `PATH`, `gippa setup` downloads the official
self-contained launcher `nextflow-26.10.0-dist` from the `nextflow-io/nextflow` GitHub release
`v26.10.0`, checks its pinned sha256 and installs it as `DIR/bin/nextflow` (recorded in
`DIR/db_manifest.json` as item `nextflow`). The `-dist` launcher downloads nothing more at its
first run, so it also works on nodes without internet access. `--no-nextflow` skips it. When a
`nextflow` is on `PATH`, setup names its version and path and downloads nothing. `gippa run` uses
`nextflow` on `PATH`, else `DIR/bin/nextflow`, and exits with code 2 when it finds neither.

**Java** is never downloaded. `gippa run` runs `$JAVA_HOME/bin/java -version` (else `java` on
`PATH`) and exits with code 2 before Nextflow starts when Java is missing or older than 17;
`gippa setup` only warns, since the database downloads need no Java. Install Java from
conda-forge (never the `defaults` channel), for example with micromamba:

```bash
micromamba create -n gippa-java -c conda-forge openjdk=21
micromamba activate gippa-java
java -version
```

To install Nextflow yourself instead, use the same environment:
`micromamba create -n gippa-nf --override-channels -c conda-forge -c bioconda nextflow=26.10.0 openjdk=21`.

GIPPA 1.0 supports genetic code 11 only. Preflight refuses any other
`annotation.translation_table`, and warns (`genome.genetic_code`) when a genome looks like it uses
code 4 ([Inputs › Genome FASTA](../user-guide/inputs.md#genome-fasta)).

## The five steps

```bash
# 1. the Python package
pip install gippa
gippa --version

# 2. Nextflow, Java and Docker or Apptainer on PATH (Prerequisites above)

# 3. pinned images, databases and the local UniOP build
export NXF_APPTAINER_CACHEDIR=/big/disk/gippa-images     # optional, apptainer only (below)
gippa setup --db-dir /path/to/db --profile apptainer   # or: --profile docker

# 4. check the inputs and the installed tools
GIPPA_DB_DIR=/path/to/db gippa preflight --genome genome.fna --peptides peptides.tsv \
    --profile apptainer                                  # or: --profile docker

# 5. run
gippa run --genome genome.fna --peptides peptides.tsv --db-dir /path/to/db \
    --outdir out --profile apptainer                     # or: --profile docker
```

`uv tool install gippa` works in place of `pip install gippa`. The extra `gippa[docs]` installs
MkDocs Material for building this site. **Bioconda:** the recipe is pending (M11-10); until it
is merged, use PyPI.

**Step 3.** `gippa setup` pulls the images, then downloads the databases. Every item is on by
default; each has an opt-out such as `--no-antismash` (the flags are on
[Set up the databases](setup-databases.md)). Before downloading, setup prints the licence terms
of geNOMAD (non-commercial use only) and KEGG (academic use only).

- With `--profile`, setup pulls the pinned images that the configuration's modules need;
  `--images all` pulls every pinned image and `--no-images` none. Without `--profile` the images
  are skipped with a warning. `--config` and `--set` select the configuration, as for
  `gippa run`.
- `--profile apptainer` runs `apptainer pull <cache>/<name>.img docker://<ref>`, where `<name>` is
  the file name Nextflow looks for ([ADR-0059](../adr/ADR-0059-bare-image-refs.md)).
  `--profile docker` runs `docker pull <ref>` and checks that the pulled digest equals the pin.
  A comma-joined profile is accepted for the images, but the UniOP step needs one runtime:
  `--profile docker,conda` exits 2 with `ambiguous UniOP runtime for docker,conda: use one of
  conda | docker | apptainer`.
- Each image is recorded as an `image:<tool>` item in `DIR/db_manifest.json`. A failed pull or a
  digest mismatch (`pulled <d1>, pinned <d2>`) is exit code 2; the other images are still pulled.
- **The Apptainer cache** is one directory: `NXF_APPTAINER_CACHEDIR`, else
  `$GIPPA_ROOT/containers`, else `DIR/containers` (`setup --cache-dir` overrides it for that
  `setup` command). Every reader uses that same directory:
    - setup's image pull and UniOP build;
    - setup's Bakta and antiSMASH downloaders, which run inside the pulled images;
    - preflight under `--profile` (`gippa preflight --profile` and `gippa run`);
    - Nextflow: `gippa run` exports `NXF_APPTAINER_CACHEDIR` to the cache it checked when the
      variable is not set (`apptainer.cacheDir`).

  `DIR/containers` therefore works with no variable set. Export `NXF_APPTAINER_CACHEDIR` to put
  the images on another disk, or to share one cache with Nextflow runs outside `gippa`; keep it
  set for setup, preflight and run. A `--cache-dir` other than that directory is seen only by
  that `setup` command: preflight then reports the images `not in the apptainer cache`.
- On a machine with Apptainer and no Docker, the Bakta (and, with `--antismash`, antiSMASH)
  downloaders run inside the tool image: setup pulls the images first, so one `setup --images`
  run installs those databases too.
- Compute nodes without network: run step 3 where the network is. Preflight checks the cache,
  not the registry.

**Step 4.** `gippa preflight` validates the inputs in seconds and reads the installed versions from
`$GIPPA_DB_DIR/db_manifest.json`. With `--profile`, it also checks that every tool is installed
under that profile, exactly as `gippa run --profile` does before Nextflow starts (next section):
the same profile names, the DB directory from `$GIPPA_DB_DIR` and the prebuilt envs under
`$GIPPA_ROOT`. A missing tool is exit code 2; for missing images the message says `run: gippa setup
--images --profile apptainer`. Without `--profile`, no tool is checked.

**Step 5.** The default `--profile` is `farmshare,apptainer`, so off FarmShare always pass it.
`--memory` defaults to `64.GB`; set it to what the machine has. See [Run](../user-guide/run.md).

### UniOP

UniOP, the operon predictor, has no licence file, so GIPPA never redistributes it: there is no
UniOP image. `setup --images` (or `setup --uniop` alone) fetches the authors' archive at the pinned
commit from GitHub, checks its sha256, and builds it on your machine. Nothing is pushed anywhere.

- **Terms.** Setup prints this notice whenever it builds or checks a local build (not under a
  `farmshare*` profile, which uses the prebuilt env):
  "UniOP (github.com/hongsua/UniOP @3ed14dfc) is fetched from its authors and built on this
  machine; it has no licence file, so gippa does not distribute it and whether your use is
  permitted is yours to decide." On a terminal it asks `Build UniOP for your own use? [y/N]`; the
  default is no. For scripts and batch jobs, pass `--accept-uniop-terms`. Without either, nothing
  is fetched: `setup --images` warns `image:uniop: terms not accepted; rerun with
  --accept-uniop-terms` and exits 0; `setup --uniop` exits 2.
- **Runtime.** micromamba when a profile enables conda, else Apptainer, else Docker. Apptainer
  runs `apptainer build --fakeroot` (plain `apptainer build` as root), and retries once as
  `apptainer build --fakeroot --ignore-fakeroot-command`. The build is named in
  `DIR/local_tools.config`, which `gippa run` passes to Nextflow.
- **When the build fails** (no network, no Docker daemon access, no fakeroot), setup records the
  reason and exits 2 (once the terms are accepted, under `--images` and `--uniop` alike). The ways
  out:
    - run `setup` again on a machine that can build, with the same DB directory and cache on a
      filesystem the run host sees (preflight checks the run host: the Docker image ID or the
      `.sif` sha256 must equal the recorded build);
    - switch profile, as the Apptainer message says: `use --profile docker or conda, or set
      annotation.uniop.enabled=auto`;
    - or run without UniOP. `annotation.uniop.enabled` is `auto` by default: with UniOP missing,
      preflight warns `uniop not installed (unlicensed, not distributed); operons use the distance
      model; run: gippa setup --uniop --profile <profile> --db-dir <db_dir> --accept-uniop-terms`
      and the run goes on. `--set annotation.uniop.enabled=true` makes a missing UniOP an error
      (exit 2) instead.
- **What the distance model loses.** Without UniOP, the operon pairs and the TUs come from the
  intergenic-distance model of [ADR-0022](../adr/ADR-0022-uniop-terminatornet-cutoffs.md)
  (linked when d ≤ 102 bp). On the RegulonDB K-12 pairs it has precision 0.828 [0.805, 0.847]
  and recall 0.916 [0.901, 0.931], against UniOP's precision 0.90 at its cut-off. The run carries
  the flag `operon_source_distance_model` and the manifest `operon_source: distance_model`. Per
  [ADR-0061](../adr/ADR-0061-distance-model-fallback-scope.md), Moran's I rows under the TU-block
  null and the discordance `both_detected_t` rows carry `tu_source_distance_model`: Moran's I is
  anti-conservative at short lags (rejection 0.34 at α = 0.05 in simulation), and discordance is
  conservative (fewer pairs are found). Under `--publication`, `auto` with UniOP missing is an
  error; set `annotation.uniop.enabled: false` to publish with the distance model.
- `gippa compare` warns `compare.operon_source` when two runs used different operon sources.

## Execution profiles

`gippa run --profile` takes one or more of these Nextflow profiles, joined by commas: `farmshare`,
`farmshare_local`, `slurm`, `apptainer`, `conda`, `docker`, `standard`. The default is
`farmshare,apptainer`.

| Profile | What it does | Status |
|---|---|---|
| `apptainer` | runs every pinned tool from its image in the Apptainer cache; binds the DB directory read-only | tested (on FarmShare) |
| `docker` | the same images through Docker: DB directory bound read-only, tasks run as your user (`-u $(id -u):$(id -g)`) | smoke-tested (`docker-smoke`); see [Platforms](../user-guide/platforms.md) |
| `conda` | lets Nextflow create each tool's environment from `workflow/envs/*.yml`; preflight still requires the prebuilt environments of the conda-only tools | not tested on its own |
| `farmshare` | Slurm executor on Stanford FarmShare, plus the prebuilt tool environments under `$GIPPA_ROOT` | tested |
| `farmshare_local` | the same environments, but all tasks run inside the current Slurm allocation | tested |
| `slurm` | generic Slurm executor: the cluster's default partition, no site-specific options; set the queue, account, extra `sbatch` options, job cap and time limit with the `params.slurm_*` options (`slurm_queue`, `slurm_account`, `slurm_options`, `slurm_max_jobs`, `slurm_max_time`); combine it with one container runtime, e.g. `slurm,apptainer` | not tested |
| `standard` | Nextflow's local executor, with no container or Conda settings | not tested |

**A Slurm (or other) executor off FarmShare.** The analysis steps (`GIPPA_STEP`) run `gippa` on the
host, not in a container. Every node that runs tasks therefore needs the same Python environment
with the same `gippa` version on `PATH`, for example one venv on a shared filesystem. `gippa run`
has no `-c` option for your own Nextflow config; Nextflow also reads `$NXF_HOME/config` (default
`~/.nextflow/config`), where an executor can be set. The DB directory and the Apptainer cache
(`NXF_APPTAINER_CACHEDIR`) must also be on a filesystem that every node sees.

## What `gippa run` needs before preflight passes

In run mode, and in `gippa preflight --profile`, preflight checks that every tool the run will use
is installed under `--profile` (`missing_tools` in `src/gippa/preflight/environment.py`). A
missing tool is an environment error: `gippa run` exits with code 2 before Nextflow starts, and
the message names the tool and the reason.

| Tool kind | `docker` | `apptainer` (no conda profile) | `farmshare*`, or with `conda` |
|---|---|---|---|
| pinned in `containers.config` | present (Nextflow pulls by digest) | present iff its `.img` is in the cache | with `apptainer` in the list: as in that column; else not checked |
| no image pin (conda-only) | missing: `conda-only tool under a container profile without conda` | same | needs the prebuilt env `$GIPPA_ROOT/micromamba/envs/gippa-<tool>/` |
| UniOP | the local build in `local_tools.config` (image ID checked) | the local `.sif` (sha256 checked) | `farmshare*`: the prebuilt env; plain `conda`: the micromamba build |
| fgsea (`function` with `--abundance`) | prebuilt env, else its `images.json` image | prebuilt env, else that image in the cache | prebuilt env |
| in-process tools | never missing | never missing | never missing |

A prebuilt env counts only with its marker file `.gippa-built`, holding the SHA-256 of
`workflow/envs/<tool>.yml`. Preflight also reads `DIR/db_manifest.json` for the databases
([Set up the databases](setup-databases.md)). The antiSMASH image alone is about 9 GB
([Set up the databases › Layout and sizes](setup-databases.md#layout-and-sizes)).

## What `gippa setup` needs

`gippa setup` fetches most items itself over HTTPS (the Nextflow launcher too, when none is on
`PATH`) and needs no Nextflow, no Java and no checkout. Two
downloaders run inside a tool environment or image: `bakta_db download` (which also runs
`amrfinder_update`) and antiSMASH's database downloader (`--antismash`). For each,
`ToolRunner` (`src/gippa/setupdb/run.py`) takes the first of:

1. a micromamba environment `gippa-bakta` or `gippa-antismash` under `$MAMBA_ROOT_PREFIX/envs/`
   (default `$GIPPA_ROOT/micromamba`);
2. the pinned image in the Apptainer cache (as in step 3), run with `apptainer exec`;
3. the pinned image through `docker run -u <uid>:<gid>`, when `docker` is on `PATH`.

When none is found, the item stays `absent` with the message `no gippa-<tool> environment, no
pinned <tool> image … run gippa setup --images --profile apptainer|docker`, and setup exits 2.

## FarmShare from a checkout

On Stanford FarmShare, GIPPA runs from a repository checkout with prebuilt environments
(`AGENTS.md › Setup`). Inside a Slurm allocation, never on a login node:

```bash
source scripts/farmshare_env.sh            # sets GIPPA_ROOT, NXF_APPTAINER_CACHEDIR=$GIPPA_ROOT/containers
bash scripts/bootstrap_farmshare.sh        # once: uv, Nextflow, nf-test, gh in user space
uv sync --frozen --all-extras
sbatch scripts/slurm/build_envs.sbatch <tool> [...]   # micromamba envs and Apptainer pulls
sbatch scripts/slurm/gippa_setup.sbatch     # the databases
sbatch scripts/slurm/gippa_run.sbatch <args>           # -profile farmshare,apptainer
```

`build_envs.sbatch` builds `$GIPPA_ROOT/micromamba/envs/gippa-<tool>/` with its `.gippa-built`
marker, including the conda-only tools and UniOP, and pulls every pinned image into
`$GIPPA_ROOT/containers`. Under `farmshare,apptainer`, the conda-only tools use their prebuilt
environments and the image-pinned tools use their images: `workflow/conf/container_common.config`,
shared by the `apptainer` and `docker` profiles, sets `conda = ''` for every image-pinned process,
so an environment never shadows an image.

## Licences

GIPPA is MIT-licensed; the tools keep their own licences ([Licences](../user-guide/licences.md)).
geNOMAD (module `mobile`) is for non-commercial use only: for commercial use, run with
`--modules core,bgc,function` and install with `--no-genomad`. KEGG (`--kegg`) is free for
academic use only. The GPL and AGPL tools (Bakta, Promotech, Promoter Calculator, OSTIR,
TransTermHP, BacTermFinder; antiSMASH) run as separate programs in their own images or
environments; GIPPA never imports or vendors their code.

## Check the installation

```bash
gippa --version
gippa setup --db-dir /path/to/db --check          # verify what is installed, download nothing
```

The [CLI reference](../reference/cli.md) lists every command and option. To move to a newer
release, see [Upgrade](upgrade.md).
