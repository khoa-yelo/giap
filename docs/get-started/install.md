# Install

GIAP has two parts:

- the Python package, which provides the `giap` command;
- the tools it runs as separate programs (Bakta, miniprot, geNOMAD, antiSMASH and others). Each
  tool runs from a pinned container image or a prebuilt micromamba environment.

`giap preflight`, `giap convert`, `giap compare` and `giap report` need only the Python package.
`giap makedb` also needs the contaminant library from `giap setup`, unless you pass
`--no-contaminants`. `giap setup` needs a Bakta and, with `--antismash`, an antiSMASH tool
environment or image (see below), but not Nextflow. `giap run` needs Nextflow, every tool
environment and image, and the databases.

!!! warning "Supported platform for `giap run` in 1.0"
    The full pipeline is supported on Stanford FarmShare, from a checkout of the repository. Its
    tool environments are built by the repository's scripts. Outside FarmShare, no supported way
    to build the conda-only tool environments exists yet; this is tracked as an open issue. See
    [What `giap run` needs](#what-giap-run-needs-before-preflight-passes) and
    [Platforms](../user-guide/platforms.md).

## Requirements

| Requirement | Version | Where it is pinned |
|---|---|---|
| Python | Python ≥ 3.11 | `requires-python = ">=3.11"` in `pyproject.toml` |
| Nextflow | 26.10.0 (the pinned version) or newer | `scripts/bootstrap_farmshare.sh`, [ADR-0002](../adr/ADR-0002-farmshare-env-pins.md) |
| Java | 17 or newer, which Nextflow needs; the bootstrap installs OpenJDK 21 | `scripts/bootstrap_farmshare.sh` |
| A container runtime or Conda | Apptainer (tested), Conda/micromamba or Docker | see [Platforms](../user-guide/platforms.md) |

## The Python package

From PyPI, with pip or uv:

```bash
pip install giap            # into the active environment
uv tool install giap        # or: an isolated tool environment with uv
giap --version
```

The PyPI package is published with each release. The extra `giap[docs]` installs MkDocs Material
for building this site.

**Bioconda.** The recipe is pending (M11-10). Once it is merged, the install command will be
`micromamba install -c conda-forge -c bioconda giap`. Until then, use PyPI.

## Nextflow and Java

Install Nextflow and a Java runtime from conda-forge and bioconda (never the `defaults` channel):

```bash
micromamba create -n giap-nf --override-channels -c conda-forge -c bioconda nextflow=26.10.0 openjdk=21
micromamba activate giap-nf
nextflow -version
```

`giap run` exits with code 2 when `nextflow` is not on `PATH`.

## Execution profiles

`giap run --profile` takes one or more of these Nextflow profiles, joined by commas: `farmshare`,
`farmshare_local`, `apptainer`, `conda`, `docker`, `standard`. The default is
`farmshare,apptainer`.

| Profile | What it does | Status |
|---|---|---|
| `apptainer` | runs the container-pinned tools (Bakta, TransTermHP, geNOMAD, antiSMASH, miniprot, BLAST+) from images pulled by digest | tested (on FarmShare) |
| `conda` | lets Nextflow create each tool's environment from `workflow/envs/*.yml`; preflight still requires the prebuilt environments of the conda-only tools (below) | not tested on its own; see [Platforms](../user-guide/platforms.md) |
| `docker` | the same images as `apptainer`, through Docker | not tested; see [Platforms](../user-guide/platforms.md) |
| `farmshare` | Slurm executor on Stanford FarmShare, plus the prebuilt tool environments under `$GIAP_ROOT` | tested |
| `farmshare_local` | the same environments, but all tasks run inside the current Slurm allocation | tested |
| `standard` | Nextflow's local executor, with no container or Conda settings | not tested |

Some tools have no container image yet: UniOP, Promotech, Promoter Calculator, TerminatorNet,
OSTIR, BacTermFinder (opt-in) and the R fgsea runner (used with `function` and `--abundance`).
These are the conda-only tools. [Platforms](../user-guide/platforms.md) lists which profile
combinations are tested.

## What `giap run` needs before preflight passes

In run mode, preflight checks that every tool the run will use is installed. A missing tool is an
environment error: `giap run` exits with code 2 before Nextflow starts. The checks are in
`src/giap/preflight/environment.py` (`missing_tools`).

1. **`$GIAP_ROOT`.** Preflight looks for tool environments and images under `$GIAP_ROOT`. When it
   is not set, every tool that needs an environment or image is reported missing.
2. **The conda-only tool environments, under every profile.** Each conda-only tool the run uses
   needs `$GIAP_ROOT/micromamba/envs/giap-<tool>/`. That directory must contain the marker file
   `.giap-built`, holding the SHA-256 of `workflow/envs/<tool>.yml`. `scripts/build_envs.sh`
   builds the environments and writes the markers. On FarmShare, run it as
   `sbatch scripts/slurm/build_envs.sbatch` from a repository checkout. It sources
   `scripts/farmshare_env.sh` and updates the lock files in `workflow/envs/locks/`. It is written
   for FarmShare: a supported way to build these environments elsewhere does not exist yet.
3. **The pinned images, under the `apptainer` profile.** Each container-pinned tool (Bakta,
   TransTermHP, geNOMAD, antiSMASH, miniprot, BLAST+) needs its image file in
   `$GIAP_ROOT/containers/`. The file is named after the reference in
   `workflow/conf/containers.config`, with `/` and `:` replaced by `-` and `.img` appended (a
   `docker://` prefix, if present, is dropped first). Preflight reads only this directory.
   `scripts/slurm/build_envs.sbatch` pulls every pinned image into `$NXF_APPTAINER_CACHEDIR`,
   which `scripts/farmshare_env.sh` sets to `$GIAP_ROOT/containers`. Under the `docker` and `conda` profiles, preflight does not check
   the pinned images: Nextflow resolves them itself.
4. **The databases**, from `giap setup` (next section).

The antiSMASH image alone is about 9 GB (GIAP_SPEC.md §3.9).

## What `giap setup` needs

`giap setup` fetches most items itself over HTTPS. It does not need Nextflow. Two downloaders run
inside a tool environment:

- `bakta_db download`, which also runs `amrfinder_update` (for `--bakta light` or `--bakta full`);
- antiSMASH's database downloader (for `--antismash`).

For each of these, `giap setup` looks for one of two things (`src/giap/setupdb/run.py`,
`ToolRunner`):

- a micromamba environment `giap-bakta` or `giap-antismash` under `$MAMBA_ROOT_PREFIX/envs/`
  (default `$GIAP_ROOT/micromamba`);
- the pinned image, already pulled into `$NXF_APPTAINER_CACHEDIR` (default
  `$GIAP_ROOT/containers`), which it runs with `apptainer exec`.

To find the image reference, it reads `workflow/conf/containers.config` from `$GIAP_REPO`. From
an installed wheel, the default location is wrong (a known issue). Set `GIAP_REPO` to a checkout
of the repository at the same version. When neither the environment nor the image is found,
the item stays absent, its record in `db_manifest.json` names the missing tool, and `giap setup`
exits with code 2. On FarmShare, the environments and images come
from `sbatch scripts/slurm/build_envs.sbatch bakta antismash`.

## Databases

The databases and models are downloaded once with `giap setup`. See
[Set up the databases](setup-databases.md).

## Check the installation

```bash
giap --version
giap setup --db-dir /path/to/db --check          # verify what is installed, download nothing
giap preflight --genome genome.fna --proteins proteins.faa
```

The [CLI reference](../reference/cli.md) lists every command and option.
