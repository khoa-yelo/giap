# Plots and report

GIAP 1.0 builds 16 plots: V1–V3, V5–V8 and V10–V18. Each plot is written to `plots/` as
`Vnn_<name>.svg`, `.png` and `.pdf`. The exact data drawn is in `plots/data/Vnn_<name>.tsv`, with
coordinates 1-based closed. Plot titles state the finding, and they are generated from the tables.
The palettes are colourblind-safe, and a category keeps one colour across all plots.

## The plots

| ID | Name | Question | Data table |
|---|---|---|---|
| V01 | `circos` | Where in the genome is the proteome? | `plots/data/V01_circos.tsv` |
| V02 | `abundance_landscape` | Is abundance clustered along the chromosome? | `plots/data/V02_abundance_landscape.tsv` |
| V03 | `locus` | What happens at this prophage, BGC, island or flagged operon? | `plots/data/V03_locus.tsv` (index) and `plots/data/V03_locus_<id>.tsv` per locus |
| V05 | `mapping_sankey` | Where did my sequences go (route and outcome)? | `plots/data/V05_mapping_sankey.tsv` |
| V06 | `mapping_fdr` | Can I trust tolerant hits and de novo peptides? | `plots/data/V06_mapping_fdr.tsv` |
| V07 | `obs_vs_exp` | Which regions and categories are detected less than detectability predicts? | `plots/data/V07_obs_vs_exp.tsv` |
| V08 | `cog_abundance` | Which functions dominate abundance? | `plots/data/V08_cog_abundance.tsv` |
| V10 | `correlogram` | Is there structure along the chromosome beyond operons? | `plots/data/V10_correlogram.tsv` |
| V11 | `operon_pairs` | Do operon partners look alike, and which other neighbours do too? | `plots/data/V11_operon_pairs.tsv` |
| V12 | `variants` | Where do peptides support ALT residues? | `plots/data/V12_variants.tsv` |
| V13 | `novel_orfs` | Which new genes or upstream extensions are supported? | `plots/data/V13_novel_orfs.tsv` |
| V14 | `enrichment` | Which functions are high or low (GSEA, detection enrichment)? | `plots/data/V14_enrichment.tsv` |
| V15 | `determinants` | What explains abundance (promoter, RBS, codon usage)? | `plots/data/V15_determinants.tsv` |
| V16 | `compare` | How do runs or search strategies differ? | `CMP/plots/data/V16_compare.tsv` |
| V17 | `qc` | Are the inputs sane? | `plots/data/V17_qc.tsv` |
| V18 | `stat_diagnostics` | Are the tests calibrated on this dataset? | `plots/data/V18_stat_diagnostics.tsv` |

Notes on particular plots:

- **V04 and V09 are deferred.** V04 was a locus browser and V09 a mass-fraction treemap. Neither
  is built in 1.0.
- **V03** draws at most `report.max_locus_plots` loci (50 by default), ordered by q-value. Novel
  ORFs are drawn by V13.
- **V14** belongs to the `function` module, so it is skipped when that module is off.
- **V16** is written only by `giap compare` (see [Comparing runs](compare.md)). Every other plot
  is written by `giap run` and `giap report`.

A plot whose input is missing is skipped with a reason. The reasons are `input_absent:<table>`,
`module_off:<module>` and `layer_missing:<layer>`. When there is nothing to draw, the plot shows an
empty state.

Some plots need more data than the plot TSV holds. These extra tables are in `plots/support/`:
`replicon_landmarks.tsv`, `gc_skew.tsv`, `variant_windows.tsv`, `gsea_running.tsv`,
`determinants_shap.tsv`, `determinants_pd.tsv`, `null_regions.tsv`, `null_islands.tsv` and
`null_determinants.tsv`. They are described in the [table reference](../reference/tables.md).

Formats and resolution come from `report.formats` (default `svg`, `png`, `pdf`) and `report.dpi`
(default 300).

## The report

`report/report.html` is one self-contained HTML file:

- Plotly.js, the CSS and the table-sorting script are inlined;
- it loads nothing from the network;
- large tables are linked by relative path rather than embedded.

The report has these sections:

1. Summary
2. Inputs and quality control (V17)
3. Peptide mapping (V5, V6)
4. Gene detection (V7, V18)
5. Genome context (V1, V2, V10, V3; islands, region tests, autocorrelation)
6. Operons (V11)
7. Variants (V12)
8. Novel ORFs (V13)
9. Function (V8, V14)
10. Determinants of detection (V15)
11. Sensitivity
12. Methods
13. Provenance (a digest of the manifest)

Table sections show the top rows, sortable and filterable, with a link to the full TSV. These
plots have an interactive Plotly version in the report: V2, V5, V6, V7, V8, V10, V11, V14, V15,
V17 and V18. The others are static SVG.

## Rebuilding plots and the report

`giap report` rebuilds `plots/` and `report/` from the tables of an existing run, without
recomputing any statistics. Every option is listed in the
[CLI reference](../reference/cli.md#giap-report).

```bash
giap report out                       # all plots and the report
giap report out --plots V02,V07       # only these plots, then the report
giap report out --no-interactive      # a static report (no Plotly)
giap report out --strict              # exit 3 if the report exceeds its size budget
```

By default, `giap report` reads the run's own `resolved_config.yaml`, so thresholds, formats and
dpi match the run. `--config` or `--set` override them, for example
`--set report.formats=[svg]`. Running `giap report` twice gives identical `plots/data/*.tsv`.

The report's size budget is 25 MB. Beyond it, `--strict` exits with code 3.
