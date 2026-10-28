# Set up the databases

`gippa setup` downloads the databases and models into one directory and records each item in
`DIR/db_manifest.json`. Run it once per installation. Then pass the directory to `gippa run` with
`--db-dir DIR`, or set `GIPPA_DB_DIR=DIR`.

## The command for the default modules

`gippa run` runs all four modules by default: `core`, `mobile` (geNOMAD), `bgc` (antiSMASH) and
`function` (COG, GO, KEGG). `gippa setup` installs everything they need by default:

```bash
gippa setup --db-dir /path/to/db --profile apptainer   # or --profile docker
```

That is the Bakta light database, the geNOMAD database, the antiSMASH databases, the Promotech
models, GO, KEGG, COG, the contaminant library and, with `--profile`, the container images. Each
item has an opt-out: `--no-genomad`, `--no-antismash`, `--no-go`, `--no-kegg`, `--no-cog`,
`--no-contaminants`, `--no-images`.

`--minimal` is for small disks: it leaves out geNOMAD and antiSMASH and pulls only the images of
the `core` and `function` modules. It is recorded as `setup.set: minimal` in `db_manifest.json`;
a later `gippa run` on that directory without `--modules` runs `core,function` and says so.

Before any download, setup prints the licence terms of the planned items:

- `geNOMAD (module mobile) is for non-commercial use only; for commercial use: gippa setup
  --no-genomad, gippa run --modules core,bgc,function`;
- `KEGG is free for academic use only; for other use: gippa setup --no-kegg`.

Each run records the terms of what it used in `run_manifest.json` (`licence_notes`), and the
report's Methods section shows them. See [Licences](../user-guide/licences.md).

To use the Bakta full database, add `--bakta full`. Runs annotated with the light database carry the
run flag `bakta_db_light`. With `gippa run --annotation bakta.json`, Bakta is not run. Preflight
still checks that a Bakta database is installed, though, so do not use `--bakta skip` for a run.

Before the first `gippa setup`, read [Install › What `gippa setup` needs](install.md#what-gippa-setup-needs):
the Bakta and antiSMASH downloaders run in a tool environment or image that `gippa setup` does not
install.

On a shared cluster, run `gippa setup` as a batch job: `--bakta full` alone downloads about 84 GB.

## Flags

Every flag is listed in the [CLI reference](../reference/cli.md#gippa-setup).

| Flag | Default | Installs | Needed by |
|---|---|---|---|
| `--db-dir` | required | the target directory | every item |
| `--bakta` | `light` | `light`, `full` or `skip`: Bakta DB v6.0 (Zenodo record 14916843; `bakta_db download` also runs `amrfinder_update`) | annotation, unless `--annotation bakta.json` is given |
| `--genomad` / `--no-genomad` | on | geNOMAD DB v1.9 (Zenodo record 14886553) | module `mobile` |
| `--antismash` / `--no-antismash` | on | antiSMASH 8.0.4 databases (skipped when the image already contains them) | module `bgc` |
| `--promotech-models` / `--no-promotech-models` | on | Promotech RF-HOT model | promoter predictions (core) |
| `--go` / `--no-go` | on | `go-basic.obo` and `goslim_generic.obo` | module `function` |
| `--kegg` / `--no-kegg` | on | KEGG KO-to-pathway links and pathway names | module `function` |
| `--cog` / `--no-cog` | on | COG 2020 category tables (`cog-20.def.tab`, `fun-20.tab`) | `cog_category` in `genes.tsv`; module `function` |
| `--contaminants` / `--no-contaminants` | on | the universal protein contaminant FASTA | contaminant labels in mapping; `gippa makedb` |
| `--minimal` | off | Bakta light, Promotech, GO, KEGG, COG, contaminants; no geNOMAD, no antiSMASH | modules `core`, `function` |
| `--images [all]` / `--no-images` | on with `--profile` | the pinned container images the modules need (`all`: every pinned image) | `gippa run` with a container profile |
| `--check` | off | nothing: only verifies what is installed | — |
| `--nextflow` / `--no-nextflow` | on | the pinned Nextflow 26.10.0 `-dist` launcher (GitHub release `v26.10.0`), only when no `nextflow` is on `PATH`; never under `--uniop` or `--check` | `gippa run` without its own Nextflow |

The shared options `--config`, `--set`, `--log-level` and `--version` work as for every command.

## Layout and sizes

The paths are relative to `--db-dir`, from `src/gippa/data/db_sources.yaml`. Approximate sizes on
disk are from GIPPA_SPEC.md §3.9. Where §3.9 gives no size, the table says "not recorded". In that
case, measure the size after your first download.

| Item | Path | Size on disk (approx.) |
|---|---|---|
| Bakta light DB v6.0 | `bakta/db-light` | 3.9 GB unpacked (the 1.3 GB archive can be deleted) |
| Bakta full DB v6.0 | `bakta/db` | about 84 GB unpacked |
| geNOMAD DB v1.9 | `genomad_db` | about 0.8 GB; `gippa setup` asks for 4 GiB of free space before downloading |
| antiSMASH databases | `antismash` | not recorded |
| Promotech RF-HOT | `promotech` | not recorded |
| GO | `go/go-basic.obo`, `go/goslim_generic.obo` | not recorded |
| KEGG | `kegg/ko_pathway.tsv`, `kegg/pathway_names.tsv` | not recorded |
| COG 2020 | `cog/cog-20.def.tab`, `cog/fun-20.tab` | not recorded |
| Contaminants | `contaminants/contaminants.fasta` | not recorded |
| Nextflow launcher 26.10.0 | `bin/nextflow` | 47 MiB (49,084,616 bytes) |

The Apptainer images take a separate 15–25 GB in the image cache (§3.9).

## `db_manifest.json`

`gippa setup` writes `DIR/db_manifest.json`, next to the items. Each item gets a record with its
name, version, URL, checksum, size and retrieval date. The schema is
`src/gippa/schemas/db_manifest.schema.json`. `gippa run` copies the database versions and checksums
into the run's `run_manifest.json`. Preflight reads the manifest to check that the database each
module needs is installed.

KEGG has no release version. For KEGG, the manifest records the retrieval date and the SHA-256 of
the downloaded bytes ([ADR-0017](../adr/ADR-0017-unversioned-sources.md)). The contaminant
library is pinned: commit `bddf6da9`, with a fixed SHA-256 in `db_sources.yaml`.

## Failures

`gippa setup` exits with code 2 in two cases: a source cannot be reached, or a checksum does not
match. The message names the host, and the item stays absent. GIPPA never falls back to another
source. Run the command again later; `--check` shows which items are present.

## What `gippa setup` does not download

- **RegulonDB.** Users never download RegulonDB. GIPPA's developers use RegulonDB v12 only to
  calibrate cut-offs on E. coli K-12, and GIPPA does not redistribute it.
- **Search-engine software and spectra.** GIPPA never searches spectra. You run DIA-NN or a de novo
  tool yourself, under your own licence (see [Licences](../user-guide/licences.md)).
- **Genomes.** You provide the genome FASTA.

## Licences of the downloaded data

geNOMAD is for non-commercial use only. KEGG is free for academic use; other users need a licence
from KEGG. Read [Licences](../user-guide/licences.md) before you install either.
