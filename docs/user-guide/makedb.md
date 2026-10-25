# Building a search database

`giap makedb` writes a protein FASTA of the strain for your own search engine. It is optional and
a convenience only. GIAP never searches spectra (GIAP_SPEC.md §16), and no FDR is built on this
database: decoys and FDR are the search engine's job. Every option is listed in the
[CLI reference](../reference/cli.md#giap-makedb).

```bash
giap makedb --annotation bakta.json --db-dir /path/to/db --outdir strain_db
giap makedb --annotation bakta.json --six-frame-orfs --min-orf-aa 20 --outdir strain_db_orfs
```

| Option | Meaning |
|---|---|
| `--annotation` | the Bakta JSON of the strain; required in 1.0 (see below) |
| `--genome` | genome FASTA; default: the sequences inside the JSON |
| `--outdir` | the database directory to write (required) |
| `--db-dir` | the `giap setup` directory with the contaminant library; default `$GIAP_DB_DIR` |
| `--six-frame-orfs` | add class-tagged six-frame ORFs |
| `--min-orf-aa` | minimum ORF length in aa (default 20; at least 1) |
| `--patch` | an `annotation_patch.gff3` from a previous `giap run`, whose patched proteins are added |
| `--no-contaminants` | do not append the contaminant proteins |

## What it writes

`--outdir` receives three files:

- **`strain_db.faa`**, the FASTA, wrapped at 60 columns, with no `*` and no decoys. It holds:
  - the Bakta CDS and sORF proteins, as `>{locus_tag} {product}`;
  - with `--patch`, the patched proteins;
  - with `--six-frame-orfs`, the six-frame ORFs of at least `--min-orf-aa` residues that are not in
    frame with an annotated CDS, as `>GIAPORF_{class}_{contig}_{start}_{end}_{p|m}`, with
    1-based closed coordinates. The classes are `intergenic`, `cds_outframe`, `antisense` and
    `nterm_extension`;
  - the contaminants, as `>CONT_{id}`, unless `--no-contaminants` is given.
- **`db_id_map.tsv`**: one row per FASTA record. Coordinates are 1-based closed, and they are
  unwrapped across the origin (`end > length`), as in the tables.
- **`makedb_manifest.json`**: the inputs, options and counts.

## Requirements and errors

- `--annotation` is required in 1.0. Without it, `giap makedb` exits with code 2: running Bakta
  from `makedb` is not available. Run `giap run` first, or Bakta yourself.
- Without `--no-contaminants`, the contaminant library must be installed
  (`giap setup --contaminants`). Otherwise the exit code is 2.
- An unreadable or inconsistent annotation or patch exits with code 1.

## Searching with it

Search the FASTA with your engine, letting it generate its own decoys. With DIA-NN, use
`diann --report-decoys` so that the decoy rows reach GIAP. Then convert the results with
[`giap convert`](convert.md). Pass `--search-db same` when the database was built from the
genome you run GIAP on. Exact hits then skip the placement decoy step (`same_db`; see
[Interpreting results](interpreting-results.md#placement-fdr-is-not-search-engine-fdr)).
