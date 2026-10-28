# Platforms

Where GIPPA is tested and which execution profiles are supported. Each cell is one of:

- **tested**: run by a CI job (GitHub Actions, `.github/workflows/`) or a FarmShare Slurm suite
  (`scripts/slurm/`) named in the cell;
- **smoke-tested**: a reduced run on that platform (a few real tools, not the whole pipeline);
- **untested**: expected to work in principle, but nothing in GIPPA's test suites runs it.

<!-- --8<-- [start:platforms] -->
| Platform | Python package (`gippa preflight`, `convert`, `report`, `compare`) | Full pipeline |
|---|---|---|
| Linux x86-64 + Apptainer (FarmShare) | tested (`test`, `package` on ubuntu; `nightly`) | tested (`-profile farmshare_local,apptainer`; `nightly`, `weekly`) |
| Linux x86-64 + conda | tested (`test`, `test-compat` on ubuntu) | tested for the conda-only tools only (prebuilt micromamba envs; `nightly`) |
| Linux x86-64 + Docker | tested (`test`, `package` on ubuntu) | smoke-tested (`docker-smoke`: MINIPROT and TBLASTN under `-profile test,docker`; `params.db_dir` readable in the BAKTA container) |
| macOS (arm64 / x86-64) | tested (`test-compat`, macOS, Python 3.12) | untested; Docker Desktop runs the amd64 images under emulation |
| Windows | untested; use WSL2 and follow the Linux rows | untested |
<!-- --8<-- [end:platforms] -->

## What the cells mean

- **CI jobs.** `test` runs the unit, property, golden, contract and snapshot tests on ubuntu with
  Python 3.13 on every pull request and push to main; `test-compat` repeats them with Python 3.11
  on ubuntu and Python 3.12 on macOS after each merge to main; `package` builds the wheel,
  installs it in a clean environment and runs `gippa preflight` from it.
- **`docker-smoke`** (`.github/workflows/docker-smoke.yml`) runs the real MINIPROT and TBLASTN
  processes in their pinned images (by digest, `workflow/conf/containers.config`) on the
  100-kb synthetic test genome and checks that the outputs have the format GIPPA's parsers read.
  It runs on pull requests that change `workflow/` and once a week on main. It is advisory: a
  registry outage or rate limit does not block a merge.
- **FarmShare suites.** `nightly` (`scripts/slurm/nightly.sbatch`) and `weekly`
  (`scripts/slurm/weekly.sbatch`) run the real pipeline with `gippa run` inside a Slurm allocation;
  user runs on FarmShare go through `scripts/slurm/gippa_run.sbatch`.

## Profiles and tools outside the containers

- The **pinned-container tools** (Bakta, TransTermHP, geNomad, antiSMASH, miniprot, BLAST+) run
  from the images in `workflow/conf/containers.config` under `-profile apptainer` or
  `-profile docker`.
- The **conda-only tools** have no published image yet: UniOP, Promotech, Promoter Calculator,
  TerminatorNet, BacTermFinder and OSTIR keep their per-process conda environments
  (`workflow/envs/`). Their images are built from `workflow/containers/<tool>/` and pinned by a
  release once published (UniOP never: `gippa setup --images` builds it locally, see
  [Install › UniOP](../get-started/install.md#uniop)). Under `-profile docker` or plain
  `-profile apptainer` such a tool counts as missing, because it would run on the host without its
  environment. fgsea (R, the function step) runs from its prebuilt environment, else from its
  `workflow/conf/images.json` image once one is pinned.
- **DIA-NN** is not shipped with GIPPA: you run it yourself and pass its report to
  `gippa convert --from diann`.

## Docker (Linux and macOS)

- **Docker 20.10 or later.** Under `-profile docker` every task runs with
  `docker run --pull=never` (the first option of `docker.runOptions` in
  `workflow/nextflow.config`), so runs never pull images: only `gippa setup` does. Pull the
  images first with `gippa setup --images --profile docker`. Preflight (`gippa preflight
  --profile docker` and `gippa run`) checks each pinned image the selected modules need with
  `docker image inspect <ref>` and stops with exit code 2 when one is not local, when the Docker
  daemon is not reachable, or when the `docker` client is older than 20.10 (its `docker run`
  has no `--pull` option).
- **podman.** With the `podman-docker` package, `docker` runs podman 3.0 or later, which has
  `--pull=never` and `image inspect`. This is best effort and untested: no CI job or FarmShare
  suite runs podman.
- **Your own `docker.runOptions`.** A `docker.runOptions` set in your own Nextflow config
  (`$NXF_HOME/config` or a `nextflow.config` in the launch directory) replaces GIPPA's string.
  Keep `--pull=never` as its first option, or Docker pulls a missing image in the middle of the
  run.

## macOS and Windows

- The Python package installs and is tested on macOS. The full pipeline on macOS needs Docker
  Desktop and `-profile docker`; the images are amd64-only, so on Apple silicon (arm64) Docker
  Desktop runs them under emulation, which is slower. The antiSMASH image and the Bakta full
  database are large; check disk space and Docker Desktop's disk and memory limits first.
- On Windows, use WSL2 and treat it as Linux; nothing in GIPPA's suites runs on Windows or WSL2.

## Image references and registries

`workflow/conf/containers.config` pins each image as a bare `<image>@sha256:<digest>` reference,
with no `docker://` prefix: Docker takes it as is, and Nextflow adds `docker://` for Apptainer.
The Bakta and antiSMASH images use Docker Hub short names (`oschwengers/bakta`,
`antismash/standalone`). If your own configuration sets `docker.registry` or `apptainer.registry`
(nf-core configs often set `quay.io`), Nextflow prefixes those two names with that registry, where
they do not exist. Leave both settings unset, or override the two processes' `container` with
`docker.io/…` references.
