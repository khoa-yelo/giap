# Quickstart

This page runs GIAP once, from start to finish, on a small synthetic demo: the SIM100K genome
(a 100 kb circular chromosome and an 8 kb plasmid) and a peptide table simulated for it. Both were
made by giapsim, GIAP's test simulator, so the run needs no download from NCBI or UniProt and no
search engine. You need an installation ([Install](install.md)) and the databases
([Set up the databases](setup-databases.md)).

!!! warning "Where this runs"
    `giap run` passes preflight only when the tool environments and images are in place under
    `$GIAP_ROOT` (see [Install › What `giap run` needs](install.md#what-giap-run-needs-before-preflight-passes)).
    In 1.0 this is supported on FarmShare from a repository checkout. Steps 1 and 2 work anywhere
    the Python package is installed.

!!! note "The demo inputs are synthetic"
    `docs/tutorials/data/sim100k_peptides.tsv.gz` (68 KB, 1,685 rows of one run, `run1`) is a
    generic peptide table ([Inputs › Peptide table](../user-guide/inputs.md#peptide-table-peptidestsv))
    with q-values, PEPs, scores, intensities and a few decoy rows. giapsim generated every row from
    `tests/data/genomes/sim100k.fna.gz`: it contains no real mass spectrometry data. How it was made,
    and the one command that remakes it, are in its
    [provenance note](../tutorials/data/PROVENANCE.md).

    **The results describe a simulated genome.** Use the run as a smoke test of the installation
    and to learn the output layout, not to interpret biology. No run of exactly this command has
    been recorded yet.

## 1. Get the inputs

From a repository checkout, both files are already there:

```bash
mkdir -p quickstart && cd quickstart
cp ../tests/data/genomes/sim100k.fna.gz ../docs/tutorials/data/sim100k_peptides.tsv.gz .
```

Without a checkout, download them from the repository:

```bash
mkdir -p quickstart && cd quickstart
curl -fLO https://raw.githubusercontent.com/khoa-yelo/giap/main/tests/data/genomes/sim100k.fna.gz
curl -fLO https://raw.githubusercontent.com/khoa-yelo/giap/main/docs/tutorials/data/sim100k_peptides.tsv.gz
```

The genome headers carry ` [topology=circular]`, so both replicons are read as circular. GIAP
reads gzipped FASTA and peptide tables directly.

## 2. Check the inputs

```bash
giap preflight --genome sim100k.fna.gz --peptides sim100k_peptides.tsv.gz --outdir preflight_check
```

Preflight takes seconds. It prints a table of checks (`ok`, `info`, `warn`, `error`) and writes
`preflight_check/preflight/preflight.json` and `preflight.md`. Exit code 0 means you can go on.
Any other code is explained in [Running GIAP › Exit codes](../user-guide/run.md#exit-codes); for
an input error (1), the table names the rule that failed.

On the demo, preflight exits 0 with one warning, `genome.contiguity` (a 100 kb chromosome is far
smaller than a bacterial genome), and info rows for the topology and the small plasmid.
`peptides.content` counts 1,685 rows with none dropped. The evidence tier is reported as
`decoy_picked`, provisionally: decoy pairing is checked after annotation, and the tier the run
uses is recorded in `run_manifest.json` and in the `detection_tier` column of `genes.tsv` (see
[Interpreting results](../user-guide/interpreting-results.md#detection-tiers)).

## 3. Run

```bash
giap run --genome sim100k.fna.gz \
  --peptides sim100k_peptides.tsv.gz \
  --db-dir /path/to/db \
  --threads 16 --memory 64.GB \
  --outdir out
```

`--profile` defaults to `farmshare,apptainer`. Off FarmShare, choose the profiles that fit your
machine: see [Install › Execution profiles](install.md#execution-profiles) and
[Platforms](../user-guide/platforms.md). If a run stops, rerun the same command with `--resume`.
Measured runtimes: see [Runtimes](../user-guide/runtimes.md).

### Your own data

To run on a real genome, pass its FASTA and your search results. Convert DIA-NN or Cascadia
output with [`giap convert`](../user-guide/convert.md), then pass the table with `--peptides`.
Without a peptide table, a protein list works too: for example the E. coli K-12 MG1655 assembly
GCA_000005845.2 from NCBI with the UniProt reference proteome UP000000625:

```bash
curl -fLO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz
curl -fL -o UP000000625.fasta \
  "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000000625%29"
giap run --genome GCA_000005845.2_ASM584v2_genomic.fna.gz \
  --proteins UP000000625.fasta \
  --set annotation.bakta.complete=true \
  --db-dir /path/to/db --outdir out_k12
```

A protein list gives the `protein_list` detection tier, whose gene detection is `unestimated`.
Giving a whole reference proteome as the list marks every mapped gene `detected`, so the region,
island and enrichment contrasts have nothing to compare against.
NCBI headers carry no topology tag; `--set annotation.bakta.complete=true` marks the replicons of
this complete assembly circular (see [Inputs](../user-guide/inputs.md#genome-fasta)).

## 4. The output tree

```text
out/
├── run_manifest.json      tool and database versions, parameters, input SHA-256, seeds, timings
├── resolved_config.yaml   the configuration the run used
├── preflight/             preflight.json, preflight.md
├── annotation/            genes.gff3, regions.tsv, regions.bed, promoters.bed, terminators.bed,
│                          operon_pairs.tsv, tus.tsv, tools/ (native tool outputs)
├── mapping/               peptide_hits.tsv, protein_hits.tsv, unmapped.tsv, mapping_fdr.tsv,
│                          mappings.bed, mappings.probed, mappings.gff3
├── tables/                genes.tsv (the master table), detection_summary.tsv, …
├── plots/                 Vnn_<name>.svg|png|pdf, data/Vnn_<name>.tsv, support/
├── report/report.html
├── schemas/datapackage.json
└── logs/                  events.jsonl, trace.tsv, timeline.html, nf_report.html, <step>.log
```

[Outputs](../user-guide/outputs.md) describes every directory. Each table is documented in the
[table reference](../reference/tables.md).

## 5. Open the report

`report/report.html` is one self-contained file. Plotly.js, the CSS and the scripts are inlined,
so the page loads nothing from the network. Open it in any browser:

```bash
xdg-open out/report/report.html      # Linux desktop
open out/report/report.html          # macOS
```

On a cluster, copy the file to your computer first, for example with `scp` or `rsync`. The report
links to the full TSV tables by relative path, so copy the whole `out/` directory if you want those
links to work.

## Next

- [Inputs](../user-guide/inputs.md): peptide tables, abundance and existing annotations.
- [Running GIAP](../user-guide/run.md): modules, seeds, `--publication` and exit codes.
- [Interpreting results](../user-guide/interpreting-results.md): what the q-values mean.
- [Tutorial: E. coli K-12](../tutorials/ecoli-k12.md).
