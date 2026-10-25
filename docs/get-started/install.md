# Install

GIAP has two parts:

- the Python package, which provides the `giap` command;
- the tools it runs as separate programs (Bakta, miniprot, geNOMAD, antiSMASH and others). Nextflow
  starts each tool from a pinned container or Conda environment.

`giap preflight`, `giap convert`, `giap makedb`, `giap compare` and `giap report` need only the
Python package. `giap run` and `giap setup` also need Nextflow and the tool environments.

## Requirements

| Requirement | Version | Where it is pinned |
|---|---|---|
| Python | 3.11 or newer (Python ≥ 3.11) | `requires-python = ">=3.11"` in `pyproject.toml` |
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
micromamba create -n giap-nf -c conda-forge -c bioconda nextflow=26.10.0 openjdk=21
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
| `conda` | builds each tool's environment from `workflow/envs/*.yml` with micromamba | not tested on its own; see [Platforms](../user-guide/platforms.md) |
| `docker` | the same images as `apptainer`, through Docker | not tested; see [Platforms](../user-guide/platforms.md) |
| `farmshare` | Slurm executor on Stanford FarmShare, plus the prebuilt tool environments under `$GIAP_ROOT` | tested |
| `farmshare_local` | the same environments, but all tasks run inside the current Slurm allocation | tested |
| `standard` | Nextflow's local executor, with no container or Conda settings | not tested |

Some tools have no container image yet. These are UniOP, Promotech, Promoter Calculator,
TerminatorNet, OSTIR, BacTermFinder and the R fgsea runner. Outside FarmShare they need Conda next
to the container profile. [Platforms](../user-guide/platforms.md) lists which combinations are
tested.

Apptainer caches images in `$NXF_APPTAINER_CACHEDIR`, or in `$GIAP_ROOT/containers` when that is
not set. The antiSMASH image alone is about 9 GB (GIAP_SPEC.md §3.9).

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
