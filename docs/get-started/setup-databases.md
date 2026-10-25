# Set up the databases

`giap setup` downloads the databases and models into one directory and records each item in
`DIR/db_manifest.json`. Run it once per installation. Then pass the directory to `giap run` with
`--db-dir DIR`, or set `GIAP_DB_DIR=DIR`.

## The command for the default modules

`giap run` runs all four modules by default: `core`, `mobile` (geNOMAD), `bgc` (antiSMASH) and
`function` (COG, GO, KEGG). They need this setup:

```bash
giap setup --db-dir /path/to/db --antismash --go --kegg --cog --contaminants
```

Without the extra flags, `giap setup` installs the Bakta light database, the geNOMAD database and
the Promotech models. antiSMASH, GO, KEGG, COG and the contaminant library are opt-in.

To use the Bakta full database, add `--bakta full`. Runs annotated with the light database carry the
run flag `bakta_db_light`. Bakta is skipped when you give `giap run --annotation bakta.json`, so
`--bakta skip` is enough in that case.

On a shared cluster, run `giap setup` as a batch job: `--bakta full` alone downloads about 84 GB.

## Flags

Every flag is listed in the [CLI reference](../reference/cli.md#giap-setup).

| Flag | Default | Installs | Needed by |
|---|---|---|---|
| `--db-dir` | required | the target directory | every item |
| `--bakta` | `light` | `light`, `full` or `skip`: Bakta DB v6.0 (Zenodo record 14916843; `bakta_db download` also runs `amrfinder_update`) | annotation, unless `--annotation bakta.json` is given |
| `--genomad` / `--no-genomad` | on | geNOMAD DB v1.9 (Zenodo record 14886553) | module `mobile` |
| `--antismash` / `--no-antismash` | off | antiSMASH 8.0.4 databases (skipped when the image already contains them) | module `bgc` |
| `--promotech-models` / `--no-promotech-models` | on | Promotech RF-HOT model | promoter predictions (core) |
| `--go` | off | `go-basic.obo` and `goslim_generic.obo` | module `function` |
| `--kegg` | off | KEGG KO-to-pathway links and pathway names | module `function` |
| `--cog` | off | COG 2020 category tables (`cog-20.def.tab`, `fun-20.tab`) | `cog_category` in `genes.tsv`; module `function` |
| `--contaminants` | off | the universal protein contaminant FASTA | contaminant labels in mapping; `giap makedb` |
| `--check` | off | nothing: only verifies what is installed | — |

The shared options `--config`, `--set`, `--log-level` and `--version` work as for every command.

## Layout and sizes

The paths are relative to `--db-dir`, from `src/giap/data/db_sources.yaml`. Approximate sizes on
disk are from GIAP_SPEC.md §3.9. Where §3.9 gives no size, the table says "not recorded". In that
case, measure the size after your first download.

| Item | Path | Size on disk (approx.) |
|---|---|---|
| Bakta light DB v6.0 | `bakta/db-light` | 3.9 GB unpacked (the 1.3 GB archive can be deleted) |
| Bakta full DB v6.0 | `bakta/db` | about 84 GB unpacked |
| geNOMAD DB v1.9 | `genomad_db` | about 0.8 GB |
| antiSMASH databases | `antismash` | not recorded |
| Promotech RF-HOT | `promotech` | not recorded |
| GO | `go/go-basic.obo`, `go/goslim_generic.obo` | not recorded |
| KEGG | `kegg/ko_pathway.tsv`, `kegg/pathway_names.tsv` | not recorded |
| COG 2020 | `cog/cog-20.def.tab`, `cog/fun-20.tab` | not recorded |
| Contaminants | `contaminants/contaminants.fasta` | not recorded |

The Apptainer images take a separate 15–25 GB in the image cache (§3.9).

## `db_manifest.json`

`giap setup` writes `DIR/db_manifest.json`, next to the items. Each item gets a record with its
name, version, URL, checksum, size and retrieval date. The schema is
`src/giap/schemas/db_manifest.schema.json`. `giap run` copies the database versions and checksums
into the run's `run_manifest.json`. Preflight reads the manifest to check that the database each
module needs is installed.

KEGG and the contaminant library have no release version. For these, the manifest records the
retrieval date and the SHA-256 of the downloaded bytes
([ADR-0017](../adr/ADR-0017-unversioned-sources.md)).

## Failures

`giap setup` exits with code 2 in two cases: a source cannot be reached, or a checksum does not
match. The message names the host, and the item stays absent. GIAP never falls back to another
source. Run the command again later; `--check` shows which items are present.

## What `giap setup` does not download

- **RegulonDB.** Users never download RegulonDB. GIAP's developers use RegulonDB v12 only to
  calibrate cut-offs on E. coli K-12, and GIAP does not redistribute it.
- **Search-engine software and spectra.** GIAP never searches spectra. You run DIA-NN or a de novo
  tool yourself, under your own licence (see [Licences](../user-guide/licences.md)).
- **Genomes.** You provide the genome FASTA.

## Licences of the downloaded data

geNOMAD is for non-commercial use only. KEGG is free for academic use; other users need a licence
from KEGG. Read [Licences](../user-guide/licences.md) before you install either.
