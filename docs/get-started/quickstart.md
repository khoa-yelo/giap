# Quickstart

This page runs GIAP once, from start to finish, on public inputs: the E. coli K-12 MG1655 genome
from NCBI and the K-12 reference proteome from UniProt. You need an installation
([Install](install.md)) and the databases ([Set up the databases](setup-databases.md)).

!!! note "Why a protein list"
    `giap run` needs evidence: `--peptides`, `--proteins` or both. A genome alone fails preflight
    with exit code 1 (`inputs.evidence_present`). GIAP ships no demo peptide table. This quickstart
    therefore uses the UniProt proteome as a protein list. That gives the `protein_list` detection
    tier, whose gene detection is `unestimated` (see
    [Interpreting results](../user-guide/interpreting-results.md#detection-tiers)). To use your own
    search results, convert them with [`giap convert`](../user-guide/convert.md) and pass
    `--peptides`.

## 1. Download the inputs

```bash
mkdir -p quickstart && cd quickstart
# genome: NCBI assembly GCA_000005845.2 (one circular chromosome, U00096.3)
curl -fLO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz
# proteins: UniProt reference proteome UP000000625 (E. coli K-12 MG1655)
curl -fL -o UP000000625.fasta \
  "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000000625%29"
```

NCBI headers carry no topology tag, so preflight reports the chromosome as linear and adds an info
row (`genome.topology`). This assembly is complete. You can mark its replicons circular in either
of two ways:

- add ` [topology=circular]` to the header;
- pass `--set annotation.bakta.complete=true`, as below.

See [Inputs](../user-guide/inputs.md#genome-fasta).

## 2. Check the inputs

```bash
giap preflight --genome GCA_000005845.2_ASM584v2_genomic.fna.gz \
  --proteins UP000000625.fasta --set annotation.bakta.complete=true --outdir preflight_check
```

Preflight takes seconds. It prints a table of checks (`ok`, `info`, `warn`, `error`) and writes
`preflight_check/preflight/preflight.json` and `preflight.md`. Exit code 0 means you can go on.
Exit code 1 means an input error: the table names the rule that failed.

## 3. Run

```bash
giap run --genome GCA_000005845.2_ASM584v2_genomic.fna.gz \
  --proteins UP000000625.fasta \
  --set annotation.bakta.complete=true \
  --db-dir /path/to/db \
  --threads 16 --memory 64.GB \
  --outdir out
```

`--profile` defaults to `farmshare,apptainer`. Off FarmShare, choose the profiles that fit your
machine: see [Install › Execution profiles](install.md#execution-profiles) and
[Platforms](../user-guide/platforms.md). If a run stops, rerun the same command with `--resume`.
Measured runtimes: see [Runtimes](../user-guide/runtimes.md).

## 4. The output tree

```text
out/
├── run_manifest.json      tool and database versions, parameters, input SHA-256, seeds, timings
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
