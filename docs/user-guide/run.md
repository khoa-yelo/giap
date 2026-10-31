# Running GIPPA

`gippa run` checks the inputs (preflight), then runs the analysis through Nextflow. It writes the
output tree described in [Outputs](outputs.md). This page covers the common options. The
[CLI reference](../reference/cli.md#gippa-run) lists every option with its type and default.

```bash
gippa run --genome genome.fna --peptides peptides.tsv --db-dir /path/to/db --outdir out
```

## Check the inputs first

```bash
gippa preflight --genome genome.fna --peptides peptides.tsv --abundance abundance.tsv --outdir out
```

`gippa preflight` runs the same input checks as the first step of `gippa run`, in seconds, and runs
nothing else. It writes `out/preflight/preflight.json` and `preflight.md`. `gippa run` repeats the
checks and stops before Nextflow starts when one fails. In run mode, preflight also checks that
the tools and databases of the chosen modules and profile are installed.

## Inputs

| Option | Input |
|---|---|
| `--genome` | genome FASTA (required) |
| `--peptides` | generic peptide table |
| `--proteins` | protein FASTA |
| `--abundance` | protein abundance table |
| `--annotation` | Bakta JSON, GFF3 or GenBank |
| `--reference-peptides` | same-strain database peptides for the same spectra |
| `--replicon-depth` | `contig`, `mean_depth` table |
| `--outdir` | output directory (required) |
| `--db-dir` | database directory of `gippa setup`; default `$GIPPA_DB_DIR` |

The formats are described in [Inputs](inputs.md).

## Modules

`--modules` takes a comma-separated list chosen from `core`, `mobile`, `bgc` and `function`. All
four run by default.

| Module | What it adds | Needs from `gippa setup` |
|---|---|---|
| `core` | annotation and context layers, mapping, the gene table, detection, context statistics, operons, variants, novel ORFs and determinants; always on, even when it is not listed | Bakta DB (checked even with `--annotation bakta.json`) and Promotech models; contaminants recommended (they label contaminant peptides) |
| `mobile` | prophage and plasmid regions from geNOMAD | `--genomad` |
| `bgc` | biosynthetic gene clusters from antiSMASH | `--antismash` |
| `function` | COG, GO and KEGG categories, detection enrichment, abundance GSEA | `--go --kegg --cog` |

```bash
gippa run --genome genome.fna --peptides peptides.tsv --modules core,function --outdir out
```

A module that is switched off writes no files, and the manifest lists it under `modules_off`.
Later steps run without it. For example, the region tests then have no prophage or BGC regions.
[Pipeline processes](#pipeline-processes) shows which Nextflow processes each module adds.

## Pipeline processes

`gippa run` hands the analysis to Nextflow (`workflow/main.nf`). The diagram below shows every
Nextflow process of that workflow, grouped by subworkflow (`ANNOTATE`, `MAPPING`, `GENETABLE`,
`VARIANTS`, `CONTEXT`, `OPERONS`, `DETERMINANTS`, `SENSITIVITY_RERUNS`, `FUNCTION`, `NOVEL` and
`REPORT`). Each card is one process: its name (in the Nextflow trace and `logs/trace.tsv` it
carries its subworkflow, as in `ANNOTATE:BAKTA`) and its main output files. The edges are the
channels between processes.

![GIPPA Nextflow processes: ANNOTATE runs Bakta and the promoter, terminator, RBS, operon and region tools, each followed by a layer step, and CONTEXT_LAYERS merges the layers; MAPPING builds the six-frame index, maps peptides exactly and tolerantly, places them with decoys and assigns loci, with optional miniprot and tblastn for proteins; GENETABLE builds the gene table and detectability; VARIANTS, CONTEXT, OPERONS, DETERMINANTS, SENSITIVITY_RERUNS, FUNCTION and NOVEL read it; REPORT draws the plots and report.html.](../assets/gippa-processes-light.svg#only-light)
![GIPPA Nextflow processes (dark theme).](../assets/gippa-processes-dark.svg#only-dark)

### How to read it

| Mark | Meaning |
|---|---|
| filled card | an external tool's process, run in its container or environment (Bakta, UniOP, Promotech, …) |
| outlined card | a `gippa _step` process: a Python step run through the `GIPPA_STEP` module under its own name |
| dashed card or input | optional: an optional input (`--proteins`: miniprot, tblastn; `--annotation` with a GFF3 or GBFF, whose CDS Bakta re-calls), a module that can be switched off (all four run by default; `mobile`: geNOMAD, `bgc`: antiSMASH, `function`) or a switch that is off by default (`annotation.terminators.bactermfinder`) |
| **G** badge | the process reads the normalised genome, `genome.norm.fa` |
| **A** badge | the process reads Bakta's gene model (`genes`, `features`, `replicons.parquet`) or the tables `ANNOTATE` emits (`xrefs`, `tus`, `layers`, `regions`) |
| **P** badge | the process reads the peptide table, `peptides.norm.parquet` |

The badges replace edges for these three inputs, which most processes read. The diagram does not
draw the database paths (`--db-dir`), the per-step event logs, or every second reader of a file:
for example, `sixframe.npz` also goes to `MAP_TOLERANT`, `PLACEMENT` and `LOCI`, the hits tables
to `VARIANTS_STEP` and `NOVEL_ORFS`, `PLACEMENT`'s `stratum0_classes.parquet` and
`GENE_TABLE`'s `tier.parquet` to `CONTEXT_STATS`, and `proteins.long.faa` to
`PROTEIN_LEFTOVERS` and `LOCI`.

Some solid cards can be skipped by a configuration switch that is on by default:
`annotation.uniop.enabled` (`UNIOP`; `auto` runs it when UniOP is installed),
`annotation.promoter_calculator.enabled` (`PROMOTER_CALCULATOR`, `TX_LAYER`),
`annotation.terminators.tool` (`TERMINATORNET_SHARDED` for `terminatornet`) and
`annotation.terminators.cross_check_transtermhp` (`TRANSTERM_INPUT`, `TRANSTERMHP`, which also run
when the tool is `transtermhp`).

Some processes are collapsed or replaced:

- `TERMINATORNET_SHARDED` is a subworkflow of three processes, `TERMINATORNET_PLAN`,
  `TERMINATORNET_CHUNK` (one task per chunk) and `TERMINATORNET_MERGE`. With
  `params.terminatornet_sharded = false`, a single `TERMINATORNET` process runs instead.
- A precomputed annotation layer (`annotation.precomputed.*`) replaces its tool and parse
  processes. For example, a precomputed `bakta_json` skips `BAKTA`.
- Without a peptide table, only `ANNOTATE` and `REPORT` run.

### Where the graph comes from

The process list is the one Nextflow itself reports, with every module switched on:

```bash
nextflow run workflow/main.nf -profile test -stub -preview -params-file params.json -with-dag dag.mmd
```

`-preview` builds the graph without running a task. The edges follow the channel wiring in
`workflow/subworkflows/local/*.nf`. The figure is drawn by `scripts/docs/process_figure.py`;
`python scripts/docs/process_figure.py --check` fails when the committed images are out of date.

## Analysis options

| Option | Values | Config key |
|---|---|---|
| `--genome-relation` | `related` (default) or `same`: is the proteomics from this exact strain? | `genome_relation` |
| `--abundance-type` | `ibaq`, `top3`, `maxlfq`, `intensity`, `unknown` | `abundance.type` |
| `--abundance-scale` | `auto`, `log2`, `linear` | `abundance.scale` |
| `--abundance-missing` | `median`, `exclude`, `min`, `zero`, `mean` | `abundance.missing` |
| `--promoters` | `targeted` (upstream windows; default) or `genome` (genome-wide) | `annotation.promoters.mode` |
| `--terminators` | `terminatornet` (default) or `transtermhp` | `annotation.terminators.tool` |
| `--add-bactermfinder` | also run BacTermFinder (opt-in) | `annotation.terminators.bactermfinder` |
| `--sensitivity` | `basic` (default), `full`, `off` | `output.sensitivity` |
| `--parquet` | also write Parquet twins of large tables | `output.parquet` |

Every option above is a shortcut for a configuration key. Any key can be set with
`--config FILE` or `--set key=value` (repeatable; the value is YAML). The order of precedence is: defaults,
then the config file, then the command-line flags, then `--set`. An unknown key exits with code 1.
The resolved values are written to `OUT/resolved_config.yaml` and to the manifest. All keys are listed
in the [configuration reference](../reference/configuration.md).

## Resources and execution

| Option | Default | Meaning |
|---|---|---|
| `--threads` | 16 (config `threads`) | threads; the upper limit per task |
| `--memory` | `64.GB` | memory for the workflow |
| `--profile` | `farmshare,apptainer` | comma-joined Nextflow profiles (see [Install](../get-started/install.md#execution-profiles)) |
| `--resume` | off | resume a stopped run in the same `--outdir` |

On a local executor (no `farmshare` or `slurm` profile), a `threads` you did not set is capped to
the CPUs available, with one line on stderr; a value from `--threads`, the config file or `--set`
is used as given.

**Resuming.** `--resume` passes `nextflow -resume` to Nextflow. Tasks whose inputs did not change
are taken from Nextflow's cache in `OUT/work`. Give the same inputs and options as the run you
resume. Nextflow's work directories can reach tens of GB per run (GIPPA_SPEC.md §3.9). Delete
`OUT/work` once the run is complete and you no longer need to resume it.

The `--executor` option is a developer feature. It is described in the
[Developer guide](../developer.md), not here.

## Seeds and reproducibility

Every random draw (permutations, bootstraps, decoy shuffles, cross-fitting folds) comes from one
master seed. `--seed` sets it, and the default is the config value `seed`. Each step derives its
own stream from the master seed and the step name. The same inputs, options and seed therefore
give byte-identical tables. The seed is recorded in `run_manifest.json`.

```bash
gippa run --genome genome.fna --peptides peptides.tsv --seed 42 --outdir out_seed42
```

## Publication mode

`--publication` raises the resamples to 10,000 permutations and 10,000 bootstrap draws
(`stats.permutations`, `stats.bootstrap`, and the determinants null `determinants.permutations`).
It also records the setting in the manifest. The default is 1,000 of each. Permutation p-values are
(1 + b)/(1 + B), where b counts the permuted statistics at least as extreme as the observed one,
out of B permutations. With B = 1,000, the smallest attainable p-value is therefore about 0.001.
Publication mode lowers that floor. The run takes longer; measured runtimes are on
[Runtimes](runtimes.md).

## Strict mode

`--strict` turns preflight warnings into a failure: the run stops before Nextflow starts, with
exit code 3. Without `--strict`, warnings are reported in `preflight/preflight.md` and the run
continues. `gippa preflight --strict` behaves the same way. `gippa report --strict` exits with
code 3 when the report exceeds its size budget.

## Exit codes

Every `gippa` command uses the same exit codes:

| Exit code | Meaning |
|---|---|
| 0 | success |
| 1 | invalid input: preflight errors, a bad option or an unknown configuration key |
| 2 | environment error: a tool, database or model is missing (for example `nextflow` not on `PATH`, or a `gippa setup` source that cannot be reached) |
| 3 | warnings were raised and `--strict` was given |
| 4 | a pipeline step failed |
| 5 | internal error (a bug); the message asks you to report it |

When a step fails, the manifest names it and points to its log: `failed_step` and `failed_log`,
from `OUT/logs/trace.tsv`. An internal error writes `OUT/logs/internal_error.log`. Please report
it at <https://github.com/khoa-yelo/gippa/issues> with the command line and that log.

## Logs

`OUT/logs/` holds the following:

- `events.jsonl`: structured warnings and events, one per line;
- `trace.tsv`, `timeline.html` and `nf_report.html`: Nextflow's task trace and reports;
- `<step>.log` for a step that failed.

`--log-level` (`DEBUG`, `INFO`, `WARNING`, `ERROR`) sets how much is printed to the terminal.
