# Comparing runs

`gippa compare` compares two or more runs on the same genome. A typical use is the same sample
searched in different ways: a database search, a de novo search, or two search engines. Every
option is listed in the [CLI reference](../reference/cli.md#gippa-compare).

```bash
gippa compare out_db out_denovo out_hybrid --outdir cmp --labels db,denovo,hybrid
```

| Option | Meaning |
|---|---|
| `RUNS` | two or more `gippa run` output directories |
| `--outdir` | the comparison directory to write (required) |
| `--labels` | comma-separated run labels, one per run, in the order of `RUNS`; default: the directory names (repeats get `_2`, `_3`, …) |

## Checks (exit code 1)

`gippa compare` stops with exit code 1 when any of these fails:

- fewer than two runs are given;
- `--labels` has the wrong number of labels, or labels that are repeated or empty;
- a label contains `&` or `;`, or is `none` or `other`;
- a run has no genome SHA-256 in `run_manifest.json` or `preflight/preflight.json`;
- the runs' genome SHA-256 values differ;
- a run has no `tables/genes.tsv`;
- the runs' locus-tag sets differ.

The last check means the runs must share one gene model. Use the same `--annotation` (or the same
Bakta database) for every run you want to compare.

## Outputs

```text
CMP/
├── tables/compare_detection.tsv    gene x run detection matrix, plus the exact intersections
├── tables/compare_abundance.tsv    pairwise abundance concordance
├── tables/compare_stability.tsv    pairwise Jaccard overlap of result sets
├── plots/V16_compare.svg|png|pdf   plus plots/data/V16_compare.tsv
└── report/compare.html             a self-contained comparison page
```

The columns are listed in the table reference:
[compare_detection](../reference/tables.md#compare_detection),
[compare_abundance](../reference/tables.md#compare_abundance) and
[compare_stability](../reference/tables.md#compare_stability).

- **`compare_detection.tsv`** has two kinds of rows. `matrix` rows give a gene, a run and whether
  the gene was detected. `intersection` rows give the number of genes detected in exactly a given
  set of runs, as in an UpSet plot.
- **`compare_abundance.tsv`** has one row per pair of runs. It uses the genes detected with a
  measured (not imputed) abundance in both. It gives:
  - Spearman ρ, with a 95 % TU-block bootstrap CI (`stats.bootstrap` resamples);
  - the Bland–Altman bias and limits of agreement on log2 abundance.

  The flag `too_few_shared` marks pairs with too few shared genes, and `scale_differs` marks runs
  whose abundance scales differ.
- **`compare_stability.tsv`** gives, per pair of runs and per result family, the Jaccard overlap
  of the flagged sets. The families are islands, operon flags, enriched categories, novel ORFs and
  ALT sites. `both_empty` and `missing_in:<run>` mark sets that are empty, or tables that are
  absent.

V16 shows an UpSet plot of the detected genes across runs, and the pairwise abundance concordance.
`compare.html` lists each run's directory, genome SHA-256 and detection tier, together with V16
and the tables.

Detection in each run keeps that run's own error control: its tier and gene q-values. A gene that
is detected in one run and not another is not, by itself, a significant difference between them.
Read the detection tier of each run before you compare counts (see
[Interpreting results](interpreting-results.md#detection-tiers)).
