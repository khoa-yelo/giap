# Outputs

Every run writes the same tree. Every table has a published schema, and every plot comes with the
data behind it. Each table's columns, types, units and constraints are in the
[table reference](../reference/tables.md). The controlled flag values are in the
[flag reference](../reference/flags.md).

## The output tree

```text
OUT/
├── run_manifest.json          versions, database versions and checksums, parameters,
│                              input SHA-256, seeds, per-step time and peak memory
├── resolved_config.yaml       the configuration the run used
├── preflight/                 preflight.json, preflight.md
├── annotation/
│   ├── genes.gff3             merged gene model with source-tagged layers
│   ├── regions.tsv, regions.bed
│   ├── promoters.bed, terminators.bed, terminators_crosscheck.tsv
│   ├── operon_pairs.tsv, tus.tsv, id_map.tsv (only with --annotation GFF3/GBFF)
│   └── tools/                 native outputs of each tool (bakta/, uniop/, genomad/, …)
├── mapping/                   peptide_hits.tsv, protein_hits.tsv, unmapped.tsv, mapping_fdr.tsv,
│                              mappings.bed, mappings.probed, mappings.gff3
├── tables/                    genes.tsv (master), detection_summary.tsv, regions_stats.tsv,
│                              autocorrelation.tsv, islands.tsv, pairs.tsv, operons.tsv,
│                              operon_flags.tsv, pileup.tsv.gz, variants.tsv, variants.vcf,
│                              novel_orfs.tsv, novel_orfs.gff3, novel_orfs.faa,
│                              annotation_patch.gff3, category_summary.tsv,
│                              enrichment_detection.tsv, enrichment_gsea.tsv,
│                              determinants_features.tsv, determinants_model.tsv,
│                              determinants_attribution.tsv, sensitivity.tsv
├── plots/                     Vnn_<name>.svg|png|pdf, data/Vnn_<name>.tsv, support/*.tsv
├── report/                    report.html (self-contained)
├── schemas/datapackage.json   Frictionless descriptor of every published table
└── logs/                      events.jsonl, trace.tsv, timeline.html, nf_report.html, <step>.log
```

Some files sit outside this list. `genetable/`, `context/` and `novel/` hold Parquet
intermediates that later steps and `giap report` read. `work/` is Nextflow's work directory. These
are not part of the published table set, and their layout may change between versions.

Optional modules that are switched off write no files. `run_manifest.json` lists them under
`modules_off` (`mobile`, `bgc`, `function`). The `function` module writes `category_summary.tsv`,
`enrichment_detection.tsv` and `enrichment_gsea.tsv`.

## Main tables

| Table | One row per | Reference |
|---|---|---|
| `tables/genes.tsv` | annotated CDS, sORF or pseudogene: identity, context, evidence, detection, abundance, function, flags | [genes](../reference/tables.md#genes) |
| `tables/detection_summary.tsv` | detection by tier, threshold, feature class and replicon | [detection_summary](../reference/tables.md#detection_summary) |
| `mapping/peptide_hits.tsv` | peptide × accepted locus | [peptide_hits](../reference/tables.md#peptide_hits) |
| `mapping/protein_hits.tsv` | protein alignment (miniprot, tblastn) | [protein_hits](../reference/tables.md#protein_hits) |
| `mapping/unmapped.tsv` | unmapped or contaminant peptide or protein, with the reason | [unmapped](../reference/tables.md#unmapped) |
| `mapping/mapping_fdr.tsv` | placement FDR stratum, or de novo precision bin | [mapping_fdr](../reference/tables.md#mapping_fdr) |
| `tables/regions_stats.tsv` | region or region-type test | [regions_stats](../reference/tables.md#regions_stats) |
| `tables/autocorrelation.tsv` | lag of the autocorrelation tests | [autocorrelation](../reference/tables.md#autocorrelation) |
| `tables/islands.tsv` | abundance island | [islands](../reference/tables.md#islands) |
| `tables/pairs.tsv`, `operons.tsv`, `operon_flags.tsv` | adjacent gene pair; TU; flagged pair | [pairs](../reference/tables.md#pairs), [operons](../reference/tables.md#operons), [operon_flags](../reference/tables.md#operon_flags) |
| `tables/variants.tsv`, `variants.vcf`, `pileup.tsv.gz` | REF/ALT call at a covered codon; pileup | [variants](../reference/tables.md#variants), [pileup](../reference/tables.md#pileup) |
| `tables/novel_orfs.tsv` (+ `.gff3`, `.faa`, `annotation_patch.gff3`) | novel ORF or start-site candidate | [novel_orfs](../reference/tables.md#novel_orfs) |
| `tables/category_summary.tsv`, `enrichment_detection.tsv`, `enrichment_gsea.tsv` | function category; enrichment test | [category_summary](../reference/tables.md#category_summary), [enrichment_detection](../reference/tables.md#enrichment_detection), [enrichment_gsea](../reference/tables.md#enrichment_gsea) |
| `tables/determinants_*.tsv` | gene features; model; attribution | [determinants_features](../reference/tables.md#determinants_features), [determinants_model](../reference/tables.md#determinants_model), [determinants_attribution](../reference/tables.md#determinants_attribution) |
| `tables/sensitivity.tsv` | sensitivity rerun of a key result | [sensitivity](../reference/tables.md#sensitivity) |

## Table conventions

| Aspect | Rule |
|---|---|
| Format | TSV, UTF-8, one header row, `NA` for missing values, one entity per row |
| Sort order | contig (in input FASTA order), then start, end, strand and ID |
| Column names | snake_case, with units in the name: `length_aa`, `distance_bp`, `abundance_log2`, `mw_da` |
| IDs | Bakta locus tags for genes; `NOVEL_<contig>_<start>_<p\|m>` novel ORFs; `<type>_<n>` regions; `TU_<n>`; `ISL_<n>` islands |
| Booleans | `true` / `false` |
| Flags | a `flags` column of `;`-separated values from the controlled [vocabulary](../reference/flags.md) |
| Test results | every test row has `stat`, `effect`, `ci_low`, `ci_high`, `p_value`, `q_value`, `method`, `null_model`, `n_perm` |
| Parquet | with `--parquet`, an optional Parquet twin of large tables |
| Versions | `schema_version` in the manifest; breaking column changes only in major versions |

`schemas/datapackage.json` describes every published table: each column's type, unit and
constraints, and the table's coordinate convention (`1-based-closed` or `0-based-half-open`).

## Coordinates

| Files | Convention |
|---|---|
| TSV tables, GFF3 (`genes.gff3`, `mappings.gff3`, `novel_orfs.gff3`, `annotation_patch.gff3`) | **1-based closed**: `start` and `end` are the first and last base |
| VCF (`variants.vcf`) | 1-based, positions within `[1, length]` |
| BED and proBed (`regions.bed`, `promoters.bed`, `terminators.bed`, `mappings.bed`, `mappings.probed`) | **0-based half-open**: `start` is the 1-based start minus 1, and `end` is exclusive (equal to the 1-based last base) |
| `plots/data/*.tsv` | 1-based closed, like the tables |

**Origin-crossing features.** On a circular replicon, a feature can run over the origin:

- TSV and GFF3 keep one record with `end > length`, an "unwrapped" coordinate. GFF3 marks the
  replicon with `Is_circular=true` on its `region` line. `genes.tsv` sets
  `crosses_origin = true` and the flag `crosses_origin`.
- BED and proBed write two records, named `<name>_part1` (from the start to the end of the
  replicon) and `<name>_part2` (from base 0 onward).

Linear contigs never wrap.

**Example.** The values below are illustrative, not from a real run. Take a circular replicon
`chr1` of 10,000 bp with two genes:

- `G_0001` on the minus strand, covering bases 1,201–2,100;
- `G_0002` on the plus strand, crossing the origin, covering bases 9,801–10,000 and then 1–300.

`tables/genes.tsv` (1-based closed; selected columns):

```text
locus_tag  contig  start  end    strand  crosses_origin
G_0001     chr1    1201   2100   -       false
G_0002     chr1    9801   10300  +       true
```

`annotation/genes.gff3` (1-based closed; the origin-crossing gene keeps `end > length`):

```text
##sequence-region chr1 1 10000
chr1  giap  region  1     10000  .  .  .  ID=chr1;Is_circular=true
chr1  Bakta  CDS     1201  2100   .  -  0  ID=G_0001
chr1  Bakta  CDS     9801  10300  .  +  0  ID=G_0002
```

A BED file of the same genes (0-based half-open; the origin-crossing gene is split in two):

```text
chr1  0     300    G_0002_part2  0  +
chr1  1200  2100   G_0001        0  -
chr1  9800  10000  G_0002_part1  0  +
```

Peptide placements follow the same rules. In `mapping/peptide_hits.tsv`, a placement across the
origin has `end` greater than the replicon length, `crosses_origin = true` and the flag
`origin_spanning`.

## Plots and report

`plots/` holds each plot as `Vnn_<name>.svg`, `.png` (300 dpi) and `.pdf`, as set by
`report.formats` and `report.dpi`. The exact data behind each plot is in `plots/data/Vnn_<name>.tsv`.
`report/report.html` is a self-contained page. See [Plots and report](plots-and-report.md).

## Manifest

`run_manifest.json` validates against `src/giap/schemas/run_manifest.schema.json`. It records:

- the status (`completed` or `failed`) and the exit code;
- the command line;
- tool versions, and database versions with their checksums;
- the resolved parameters and the SHA-256 of every input;
- the seeds and the detection tier;
- per-step time and peak memory;
- skipped steps and `modules_off`.

`giap compare` uses the genome SHA-256 recorded here to check that runs share a genome.
