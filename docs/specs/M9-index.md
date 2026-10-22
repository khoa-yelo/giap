# M9 "Visualisation and report": issue index

Spec: GIAP_SPEC.md §13 row M9. Deliverables: "16 plots, style module, report, `giap report`,
`giap compare` + V16". **Gate:** "**Khoa approves** the plot snapshots and report layout" (§14.8
checkpoint 3). After approval, a style change is a reviewed change (§10.1).

Sources read: §10 (whole: §10.1 conventions and the colour table, §10.2 the 16-plot catalogue,
§10.3 the report), §2.1 decision 13 (interactive self-contained report), §2.3, §1 rule 6
("detected", never "expressed"), §4 repository layout (`viz/`, `report/`, `compare/`,
`subworkflows/local/report.nf`), §5.1 (`giap report`, `giap compare`, exit codes), §5.2 `report:`
keys, `stats.alpha.*`, `stats.bootstrap`, §6.2 steps `plots`, `report`, `manifest`, `finalize`,
§6.3/§6.4 (executors, optional modules, `layer_missing`), §6.5 (coordinates, origin-crossing
features), §6.6/§6.7 (determinism, events, manifest), §7.3.4–§7.3.5, §7.4.5, §7.5.3, §7.6.5,
§7.7.5, §7.8.2 caveats, §7.9.3 (`giap compare`), §7.9.4 (`giap report`), §8.2 (R4 wording), §8.6,
§9.1–§9.5 (every table the plots read), §11.2 (snapshot layer, < 2 min), §11.4 (fixtures, R11
sizes), §12 (`ci.yml` runs the snapshot layer), §14.4 (visualisation loop and the style
checklist), §14.8, §16. ADRs read for their effect on plots and report text: ADR-0006 (value
encoding), ADR-0007 (size counts code), ADR-0012 (V3 draws promoters and terminators), ADR-0017
(the report states the KEGG retrieval date), ADR-0018 (region O − E), ADR-0020 (V10 envelopes,
`ci_unestimated`), ADR-0021 Amendment 1 (the run report carries ℓ̂, ℓ′, Î₁ per replicon; no
`uncalibrated` flag), ADR-0025 (V14 must not show `genomic_cluster` rows as calibrated), ADR-0026
(GSEA is "conditional on detection"), ADR-0027 (NES without CI), ADR-0029/0035 (ORA is TU-level,
`method = cpb_tu_ora`), ADR-0032 (univariate CIs uncalibrated), ADR-0033 (elastic-net ΔR² is the
primary attribution in V15), ADR-0034 (`null_uncalibrated`).

One issue = one PR; size counts **code** only (ADR-0007). "Stats" = the plot or table computes a
statistic: label `stats`, the statistician reviews it (§14.1). "Stats (light)" = a display
computation the statistician checks for agreement with the producing test (a smoother, a re-derived
running sum, a partial-dependence curve), not a new inference.

## How M9 is sequenced (Khoa's standing direction: MVP first)

- **Pass 1 (integration), 01–13.** The style module and the plot framework (01), the `plots` and
  `report` steps with a minimal static report and `giap report` (02), the small "plot support"
  files that carry run inputs which are not in the §9.1 tree (03), then every plot whose producing
  table is **on main today** (04–12), and an integration closer (13) that renders every available
  plot end-to-end from a real sim100k `giap run` under both executors.
- **Pass 2 (hardening and the gate), 14–26.** The interactive report (Plotly), the full report
  content, plot-support outputs from the analysis steps (model and null artefacts that tables do not
  hold), the plots whose tables land later (V12, V14, V15, V18 null panels), `giap compare` and V16,
  the automated style checklist and accessibility, the snapshot baselines for **Khoa's approval**,
  and scale and determinism hardening.

Nothing is dropped; it is scheduled. A plot whose input is missing is skipped with a recorded reason
(decision 4); it never fails the run.

## How this design lands (one branch, docs only)

| Commit | Files |
|---|---|
| 1/2 | this index, M9-01 … M9-13 (pass 1) |
| 2/2 | M9-14 … M9-26 (pass 2) |

## Issues

| Order | Issue | Pass | Scope | Depends on | Code lines | Stats |
|---|---|---|---|---|---|---|
| 1 | [M9-01](M9-01-style-and-plot-framework.md) | 1 | `giap.viz.style` (§10.1 palettes and widths), plot registry, `RunView` table loader, writers (svg/png/pdf + data TSV), empty states, binning, title helpers | none | 360 | no |
| 2 | [M9-02](M9-02-plots-report-steps-and-cli.md) | 1 | `plots` and `report` steps, minimal static `report.html` (13 sections), `giap report`, `REPORT` subworkflow, local publish | 01 | 380 | no |
| 3 | [M9-03](M9-03-plot-support-run-inputs.md) | 1 | `plots/support/`: replicon landmarks (oriC, terminus), GC skew, from run inputs that are not in `OUT/` (ADR-M9-A) | 02 | 160 | no |
| 4 | [M9-04](M9-04-v17-qc.md) | 1 | V17 `qc` | 01 | 220 | no |
| 5 | [M9-05](M9-05-v05-v06-mapping.md) | 1 | V5 `mapping_sankey` (custom alluvial), V6 `mapping_fdr` | 01 | 340 | no |
| 6 | [M9-06](M9-06-v07-v08-detection.md) | 1 | V7 `obs_vs_exp` forest, V8 `cog_abundance` | 01 | 320 | no |
| 7 | [M9-07](M9-07-v18-diagnostics-pass1.md) | 1 | V18 `stat_diagnostics` panels (a) reliability and (b) p-value histograms | 01 | 220 | no |
| 8 | [M9-08](M9-08-v01-circos.md) | 1 | V1 `circos` (pyCirclize) | 01, 03 | 380 | no |
| 9 | [M9-09](M9-09-v02-v10-landscape-correlogram.md) | 1 | V2 `abundance_landscape`, V10 `correlogram` | 01, 03 | 340 | stats (light) |
| 10 | [M9-10](M9-10-v03-locus-plots.md) | 1 | V3 `locus_<id>` renderer (pyGenomeViz), locus selection and ordering | 01 | 340 | no |
| 11 | [M9-11](M9-11-v11-operon-pairs.md) | 1 | V11 `operon_pairs` | 01 | 260 | no |
| 12 | [M9-12](M9-12-v13-novel-orfs.md) | 1 | V13 `novel_orfs` (a) counts, (b) top-20 locus views | 10 | 160 | no |
| 13 | [M9-13](M9-13-pass1-integration.md) | 1 | **pass-1 closer:** real sim100k `giap run` with plots and report, both executors; `giap report` idempotence; link and offline checks | 02–12 | 100 | no |
| 14 | [M9-14](M9-14-interactive-report.md) | 2 | Plotly.js inlined once, figure embedding, sortable/filterable tables, `--no-interactive`, the 25 MB budget | 02 | 320 | no |
| 15 | [M9-15](M9-15-plotly-twins.md) | 2 | interactive twins of V2, V5 (Sankey), V6, V7, V8, V10, V11, V17, V18 | 14, 04–09, 11 | 380 | no |
| 16 | [M9-16](M9-16-report-content.md) | 2 | all summary cards, methods and caveats snippets, sensitivity section, provenance, ADR-required caveats, wording checks | 02 | 360 | no |
| 17 | [M9-17](M9-17-plot-support-analysis.md) | 2 | plot support from analysis steps: GSEA running-score data, per-feature SHAP, partial dependence, null draws (regions, islands, determinants) | 03; #178, #181 merged | 340 | stats (light) |
| 18 | [M9-18](M9-18-v12-variants.md) | 2 | V12 `variants` (+ `variant_windows.tsv` support) | 03, 10; #164, #175 merged | 320 | no |
| 19 | [M9-19](M9-19-v14-enrichment.md) | 2 | V14 `enrichment` | 17; #181 and its `FUNCTION` wiring merged | 340 | stats (light) |
| 20 | [M9-20](M9-20-v15-determinants.md) | 2 | V15 `determinants` | 17; #178 merged | 340 | stats (light) |
| 21 | [M9-21](M9-21-v18-null-panels.md) | 2 | V18 panel (c): null distributions against the observed statistic | 07, 17 | 160 | no |
| 22 | [M9-22](M9-22-compare-tables.md) | 2 | `giap compare`: genome check, `compare_detection`, `compare_abundance` (ρ + TU-bootstrap CI, Bland–Altman), `compare_stability` (Jaccard), schemas | none | 380 | **yes** |
| 23 | [M9-23](M9-23-v16-compare-plot.md) | 2 | V16 `compare` (custom matplotlib UpSet, concordance), `report/compare.html` | 22, 14 | 320 | no |
| 24 | [M9-24](M9-24-style-checklist.md) | 2 | the §14.4 style checklist as code: units, overlaps, font size, palette, deuteranopia, element count, title-states-finding | 01 | 300 | no |
| 25 | [M9-25](M9-25-snapshot-baselines.md) | 2 | **gate:** pytest-mpl snapshot suite on a frozen fixture, report-structure snapshot, CI `--mpl` job, the approval packet for Khoa | 13, 24 (and each plot) | 150 | no |
| 26 | [M9-26](M9-26-scale-and-determinism.md) | 2 | 5-Mb render, K-12 real-data report, byte-identical images, `-resume`, benchmarks | 13 | 140 | no |

Total ≈ 7,400 code lines (≈ 3,600 pass 1, ≈ 3,800 pass 2), plus an estimated 6,000 test lines. That
is more than the spec's 1.5 weeks at one PR at a time; it fits only if the plot issues run in
parallel (04–12 are independent once 01 is merged) and pass 2's plot issues start as soon as their
producing PRs merge.

## Dependency DAG

```text
Pass 1   01 ─► 02 ─► 03 ─┬─► 08 (V1)
          │               └─► 09 (V2, V10)
          ├─► 04 (V17)  05 (V5, V6)  06 (V7, V8)  07 (V18 a,b)  11 (V11)
          └─► 10 (V3) ─► 12 (V13)
         02 … 12 ─► 13 (real sim100k run: every available plot in report.html)
Pass 2   02 ─► 14 ─► 15 (Plotly twins)      02 ─► 16 (report content)
         03 + #178 + #181 ─► 17 ─┬─► 19 (V14)   ◄─ #181 + FUNCTION wiring
                                 ├─► 20 (V15)   ◄─ #178
                                 └─► 21 (V18 c) ◄─ 07
         03 + 10 + #164/#175 ─► 18 (V12)
         22 ─► 23 (V16, compare.html) ◄─ 14
         01 ─► 24 (style checklist) ─► 25 (snapshot baselines; Khoa approves) ◄─ 13 and every plot
         13 ─► 26
```

## Waves

| Wave | Issues | Notes |
|---|---|---|
| A (now) | 01, 22 | 22 needs no plot code; it reads `OUT/` trees and the stats engine |
| B | 02, 04, 05, 06, 07, 10, 11 | six plot PRs in parallel once 01 merges; each is testable on in-memory frames |
| C | 03, 08, 09, 12, 24 | 08 and 09 draw the landmark and GC-skew tracks only when 03's files exist |
| D | 13 | pass-1 evidence: one real `sbatch scripts/slurm/giap_run.sbatch` on sim100k with `plots/` and `report/report.html` |
| E (hardening) | 14, 16, 17, 18, 23; then 15, 19, 20, 21 | 17–20 wait for #164/#175/#178/#181; 18–20 are coded against the branch schemas listed below |
| F (gate) | 25, 26 | 25 renders the approval packet; Khoa approves; the baselines are locked |

## Already built (checked on main, cf41ca4)

| Piece | State |
|---|---|
| `src/giap/viz/__init__.py`, `src/giap/report/__init__.py`, `src/giap/compare/__init__.py` | docstring-only packages |
| `giap report OUT [--no-interactive] [--plots …]`, `giap compare RUNS … --outdir --labels` | Click commands exist with the §5.1 options; both call `_not_implemented(…, "M9")` (exit 5) |
| `report:` config keys (`interactive`, `max_locus_plots`, `formats`, `dpi`) | in `defaults.yaml` and the pydantic model; **no change needed** (protected) |
| `plots` / `report` / `manifest` / `finalize` steps | **not registered**; `workflow/subworkflows/local/report.nf` does not exist |
| `tests/snapshot/` | does not exist; `ci.yml` runs the tests tree **without** `--mpl`, so pytest-mpl tests would run as smoke tests only |
| `giap.determinants.CAVEATS` | on the #178 branch (M8-10): the §7.8.2 caveats as constant text for the report |
| `run_manifest.json` | schema exists (`schemas/run_manifest.schema.json`), **no writer**: no step writes it yet (decision 26) |

Dependencies: matplotlib 3.11.2 (PSF-style), pyCirclize 1.10.1 (MIT), pyGenomeViz 1.7.0 (MIT),
UpSetPlot 0.9.0 (BSD-3), plotly 7.1.0 (MIT; bundles plotly.js 4.1.1, MIT, 4.8 MB), Jinja2 3.1.6
(BSD), pytest-mpl 0.19.0 (BSD-2) are **all already declared** in `pyproject.toml` and installed in
the shared venv (checked in the session). None is GPL; `scripts/ci/check_imports.py` forbids none
of them. **M9 needs no new dependency.** Not available and not proposed: kaleido (Plotly static
export; matplotlib draws every static file), colorspacious/daltonlens (the deuteranopia check uses
the published Machado 2009 matrix in numpy, M9-24), markdown (installed only transitively through
the docs extra, so not usable at run time; M9-16 uses plain HTML snippets).

**UpSetPlot is broken on the locked pandas.** In the session, `upsetplot.UpSet(…).plot()` on
pandas 3.0.6 raised `ValueError: Invalid RGBA argument: nan` (a Copy-on-Write `fillna(inplace=True)`
in `upsetplot/plotting.py:798`). §10.2 names the fallback ("custom matplotlib"); M9-23 uses it.
Removing the unused dependency is a later lock change with its own ADR (decision 29).

## What M9 consumes (producing table → plot)

Status: **main** = written by a `giap run` on main today; **PR** = the producing step is in an open
PR (schema read from that branch with `git show origin/<branch>:<path>`); **none** = no table holds
it (plot-support files, decision 3).

| Plot | Tables read (§9.1 path) | Producer (milestone, status) | Missing today |
|---|---|---|---|
| V1 circos | `preflight/preflight.json` (replicons), `tables/genes.tsv`, `mapping/peptide_hits.tsv`, `annotation/{tus.tsv, promoters.bed, terminators.bed, regions.tsv}`, `tables/islands.tsv`; `tables/variants.tsv`, `tables/novel_orfs.tsv` | M1-12, M4-15/M5-09, M3-16, M2-20, M5-06 (**main**); M7-04 (**PR #164/#175**); M7-08 (**main**), M7-09 (**PR #180**, columns already in the schema) | GC skew (needs the genome), oriC/ter: M9-03 support |
| V2 abundance_landscape | `tables/genes.tsv`, `tables/islands.tsv`, `preflight/preflight.json` | M4-15, M5-06 (**main**) | oriC/ter lines: M9-03 support |
| V3 locus | `tables/genes.tsv`, `mapping/peptide_hits.tsv`, `annotation/{promoters.bed, terminators.bed, regions.tsv}`, `tables/{regions_stats, islands, operon_flags, pairs}.tsv` | M2-20, M3-16, M4-17/18, M5-06, M6-07 (**main**) | — |
| V5 mapping_sankey | `mapping/{peptide_hits, protein_hits, unmapped}.tsv` | M3-16 (**main**) | — |
| V6 mapping_fdr | `mapping/mapping_fdr.tsv` | M3-10/M3-11/M3-16 (**main**) | — |
| V7 obs_vs_exp | `tables/detection_summary.tsv` (section `basic`), `tables/regions_stats.tsv` | M4-16, M4-17/18 (**main**) | region-**type** rows: no table (Q in M9-06) |
| V8 cog_abundance | `tables/genes.tsv` (`cog_category`, `abundance_log2`, `abundance_imputed`, `detected`) | M4-15 (**main**) | — |
| V10 correlogram | `tables/autocorrelation.tsv` | M5-05/M5-09 (**main**) | — |
| V11 operon_pairs | `tables/pairs.tsv`, `tables/operon_flags.tsv` | M6-04/M6-07/M6-08 (**main**) | — |
| V12 variants | `tables/variants.tsv` (+ `peptide_hits.tsv`) | M7-03 writer (**main**), M7-04/M7-05 step (**PR #164, #175**) | reference windows: M9-18 support |
| V13 novel_orfs | `tables/novel_orfs.tsv` | M7-07/M7-08 (**main**); M7-09 evidence columns (**PR #180**) | — |
| V14 enrichment | `tables/enrichment_gsea.tsv`, `tables/enrichment_detection.tsv` | M8-06a (**PR #181**; `FUNCTION` Nextflow wiring is a later M8-06b) | running-score data: M9-17 |
| V15 determinants | `tables/determinants_model.tsv`, `tables/determinants_attribution.tsv` | M8-10 (**PR #178**) | per-feature SHAP, partial dependence: M9-17 |
| V16 compare | `CMP/tables/compare_*.tsv` | M9-22 | — |
| V17 qc | `preflight/preflight.json`, `tables/genes.tsv`, `mapping/{peptide_hits, unmapped}.tsv` | M1-12, M4-15, M3-16 (**main**) | — |
| V18 stat_diagnostics | `tables/detection_summary.tsv` (section `detectability`), every table with `p_value` | M4-13/M4-16 and all test tables (**main** / PR) | null draws: M9-17 |
| report | all of the above, `tables/sensitivity.tsv` (M4-19, **main**), `context/island_diagnostics.parquet` (M5-09, **main**; ADR-0021), `logs/` events, `run_manifest.json` (**none**) | | manifest writer (decision 26) |

## Decisions taken in these specs (for the lead to confirm)

1. **Plot IDs and names.** Files use the two-digit ID (`V01_circos`); `--plots` accepts `V1` or
   `V01` and normalises; an unknown ID is exit 1 (M9-02 Q1).
2. **Two steps, one input map.** `plots` and `report` are §6.2 steps whose inputs are the published
   files, keyed by the work keys that `giap.runner.local.PUBLISHED` / `PUBLISHED_FILES` already map
   to `OUT/` paths. `giap report OUT` runs the same two step functions with inputs resolved from
   `OUT/` (M9-01 `RunView`, M9-02).
3. **ADR-M9-A "plot support files".** A plot input that no §9.1 file holds (GC skew, oriC/ter,
   variant reference windows, GSEA running scores, per-feature SHAP, partial dependence, null draws)
   is written once to `OUT/plots/support/<name>.tsv[.gz]` by the step that can compute it, listed
   in `schemas/datapackage.json` with `role: plot_support`, and read by `giap report`. No plot refits
   a model or reruns a null (§7.9.4) (M9-03, M9-17).
4. **Missing input, three states.** (a) module off or the required file absent → no plot files; an
   event `plots.skipped {plot, reason}` with `reason` ∈ {`module_off:<m>`, `layer_missing:<l>`,
   `input_absent:<file>`}; the report section says why. (b) the input exists but has no rows → an
   **empty-state** figure ("No islands in this run") and a header-only data TSV. (c) a plot raising
   an exception fails the `plots` step (exit 4); no silent fallback (M9-01 Q2).
5. **Titles** are produced by a per-plot function from the plot's own data TSV; they state numbers
   and a claim only where a table row carries its error estimate; otherwise they say `unestimated`
   (§10.1, R4).
6. **Byte-identical outputs.** Data TSVs use the §9.2 writer (ADR-0006, `%.6g`) and are identical
   across executors and reruns (§7.9.4). Images are made deterministic too: `svg.hashsalt` fixed,
   SVG/PDF dates removed, PNG metadata stripped (M9-01).
7. **Fonts and sizes.** DejaVu Sans (ships with matplotlib, no system dependency); 7 pt base, 8 pt
   titles; widths 3.5 in or 7.2 in (§10.1); backend Agg (M9-01 Q3).
8. **COG super-categories** use the four NCBI groups with Okabe–Ito hues: information storage and
   processing `#0072B2`, cellular processes and signalling `#009E73`, metabolism `#E69F00`, poorly
   characterised `#999999` (the Okabe–Ito grey); letters are labelled. The 26 letter descriptions
   are bundled in `giap/data/cog_categories.tsv` (§4 lists "COG categories" as bundled data) (M9-01
   Q4).
9. **Detection-status colours:** `detected` `#0072B2`, `not_detected` and `partial_excluded`
   `#BBBBBB` (partial hatched), `single_peptide_unestimated` `#E69F00` (§10.1 "unestimated").
10. **Island hatching:** high `///`, low `\\\` (§10.1 says "with hatching", no pattern).
11. **Data-TSV coordinates are 1-based closed** like every TSV (§9.2); origin-crossing features keep
    `end > length` (§6.5).
12. **The 20,000-element rule is per track** (§10.1 "when a track would draw more than 20,000
    elements"); the §14.4 checklist item is evaluated per axes/track. A track over budget is binned
    into equal-width bins (M9-01 Q5, M9-24).
13. **Abundance plots draw measured values only** (`abundance_imputed = false`); imputed genes are
    counted in the data TSV and the title, never drawn as points (V2, V3, V8).
14. **V5 routes and outcomes** are derived from the three mapping tables as tabled in M9-05; an
    unmapped sequence goes through route `unmapped`.
15. **V7 uses O − E** (the tables carry its CI; O/E has none) and does not read
    `category_summary.tsv` (one COG estimate per plot) (M9-06 Q1).
16. **V11 densities are density-normalised histograms** with shared bins, not KDEs (no bandwidth
    choice) (M9-11 Q1).
17. **V2's trend is `giap.stats.smooth` LOESS** at `context.islands.detrend_span_frac`, circular on
    circular replicons; the statistician confirms it is the displayed analogue of the islands
    detrend (M9-09 Q1).
18. **Locus pool and order (V3):** prophage and BGC regions (`regions_stats`, `stat = o_minus_e`,
    `subset = all`), islands, and the TU of each flagged operon pair, ordered by their own error
    estimate (q, or FWER for islands), NA last, ties by type then ID; at most
    `report.max_locus_plots`. Novel ORFs are V13's (decision 19) (M9-10 Q1).
19. **V13(b)** writes `V13_novel_orfs_locus_<orf_id>.*` for the top 20 `called` ORFs by `rank`,
    with the V3 renderer (M9-12 Q1).
20. **Locus window** = feature span ± 2 kb, clipped on linear contigs, unwrapped on circular ones
    (M9-10 Q2).
21. **Report embedding:** static figures inline as `<svg>`; tables show the top 20 rows with a
    relative link to the full TSV; the 13 sections have fixed IDs (M9-02).
22. **Snapshot comparisons run on one CI leg** (ubuntu-latest, Python 3.12, `--mpl`), on a frozen
    OUT-tree fixture under `tests/data/snapshot/`; the other legs run the snapshot tests without
    `--mpl` (font rasterisation differs across OS) (M9-25 Q1).
23. **Plotly.js** is inlined once from `plotly.offline.get_plotlyjs()`; no MathJax, no topojson
    fetches; `--no-interactive` / `report.interactive: false` omits it (M9-14).
24. **The table-sorting script** is ≈ 40 lines of hand-written vanilla JS in the template; no
    third-party JS (M9-14).
25. **Wording:** a test fails if a title, a report string or a template contains "express"
    (§1 rule 6) (M9-16).
26. **Provenance without a manifest.** No step writes `run_manifest.json` yet. The report shows the
    manifest digest when the file exists; otherwise the preflight input SHA-256s, the GIAP version
    and the SHA-256 of `resolved_config.yaml`, marked "manifest not written". `giap compare` checks
    genome identity through `preflight.json` until the manifest exists (M9-16 Q1, M9-22 Q1).
27. **`giap compare` tables** all live under `CMP/tables/`; `compare_detection.tsv` is long format
    with `section` ∈ {`matrix`, `intersection`} (dynamic per-run columns cannot have a fixed schema)
    (ADR-M9-B, M9-22 Q2).
28. **`giap compare` joins genes by `locus_tag`** and requires the same locus-tag set in every run
    (exit 1 otherwise, with counts) (M9-22 Q3).
29. **UpSetPlot is not used** (broken on pandas 3; above); V16 draws a custom UpSet. The dependency
    stays in `pyproject.toml` until the next lock update (no lock change in M9).
30. **Colour-blind check:** the Machado 2009 deuteranopia matrix (severity 1) applied to the
    rendered PNG; the minimum pairwise colour distance between the categories of each plot is
    **recorded**, not gated, until Khoa sets a threshold (M9-24 Q1).
31. **Minimum font size 6 pt** at the final figure size (legible at single-column width) (M9-24 Q2).
32. **RNG streams:** `compare.bootstrap` only (M9-22); every plot is deterministic.

## ADRs M9 owes (numbers: next free at merge time)

- **ADR-M9-A "Plot support files under `plots/support/`"** (03, 17, 18): adds a directory to the
  §9.1 tree; non-breaking; lead-level (decision 3).
- **ADR-M9-B "`giap compare` output layout"** (22): fills §7.9.3's gaps (paths, long format,
  locus-tag join, genome check through `preflight.json` until the manifest exists); lead-level.
- Conditional: if Khoa wants a colour-distance gate (decision 30) or a different font minimum
  (decision 31), the checklist thresholds get an ADR (they become reviewed style rules after the
  gate, §10.1).

## Khoa's checkpoint (§14.8 item 3)

M9-25 produces `local/review/m9-approval/` (gitignored): every plot rendered on the frozen fixture
and on the K-12 run (M9-26), its deuteranopia simulation, the checklist score, and `report.html`.
Khoa approves the baselines and the layout; the lead then merges M9-25 and the snapshots are locked.
Iterations follow the §14.4 loop (checklist first, then Khoa's preference); after 3 rejected rounds
on one plot, R9 applies.

## Backlog items folded in (docs/backlog.md)

| Backlog line | Where it lands |
|---|---|
| "M9 visualisation and report" | this index; the checkbox closes when M9-25 is approved and merged |
| ADR-0021: "the §9.4 run report carries ℓ̂, ℓ′ and Î₁ per replicon" | M9-16 (context section, from `context/island_diagnostics.parquet`) |
| M6-06: "`bin_counts` and π₀ must reach the §9.4 run report" | M9-16 (operons section, from the `operons.missed_operons` event); the manifest half stays with the manifest step |
| §7.8.2 "Caveats (printed in the report)", ADR-0017 KEGG retrieval date, §7.4 σ70/intrinsic-only note | M9-16 |
| "§7.8.2's separate report of ribosomal genes" | M9-20 (V15 marks `cai_reference` genes separately) and M9-16 (a table) |

## Risks

- **Producing PRs still open.** V12, V14, V15 are coded against the branch schemas
  (`git show origin/m8/function-step:…`, `origin/m8/determinants-step:…`, `origin/m7/variant-flags:…`).
  A schema change in review means a follow-up; the `RunView` reader validates columns by name so a
  mismatch fails loudly.
- **sim100k is degenerate for many plots.** The 2026-10-22 validation run
  (`local/runs/val-mvp-2026-10-22`) has 87 genes, no islands, no operon flags, no protein hits, and
  every `basic` row with `n_detected = 0`. Pass 1's acceptance is "every plot renders or records its
  skip reason", not "every plot is interesting"; the snapshot fixture (M9-25) is a richer frozen
  simulation, and the K-12 run (M9-26) is the real-data check.
- **Snapshot fragility across platforms.** RMS ≤ 2 is tight for text rendered by different FreeType
  builds; decision 22 confines the comparison to one leg.
- **Report size.** Plotly.js is 4.8 MB; 50 locus SVGs and up to 20 V13 views at ≈ 100 kB each add
  ≈ 7 MB; the 25 MB budget holds on K-12 if tables are linked, not embedded (M9-14 measures it).
- **No manifest writer** (decision 26): provenance and the KEGG retrieval date are incomplete until
  someone builds the §6.7 `manifest` step.

## Feasibility on FarmShare

Every plot is matplotlib on Agg, single-threaded, seconds per plot on sim100k (pyCirclize rendered a
two-sector, three-track SVG in 0.5 s in the session). The `plots` step runs at `process_low`. A
5-Mb genome needs binning only for the per-gene tracks (≈ 4,500 genes, far under 20,000). No Slurm
job beyond the pass-1 `giap_run.sbatch` and the M9-26 K-12 run.

## Questions for the lead (summary; details per issue)

1. ADR-M9-A (support files) — accept the directory and the "no refit in `giap report`" rule?
   **Proposed default:** yes.
2. Who builds the §6.2 `manifest` and `finalize` steps? They are not in any milestone's
   deliverables, and M9's provenance section and `giap compare` want the manifest. **Proposed
   default:** the lead schedules them as an MVP issue before M9-16; M9 reads the manifest when
   present (decision 26).
3. V7 "region types": §10.2 lists them, but no table holds a region-**type** O − E (§7.3.5(a) has
   replicon, feature class and COG; `regions_stats` has individual regions). **Proposed default:**
   V7 omits region-type rows; a backlog line asks M4 owners whether `detection_summary` should gain
   them (computing them in the plot would recompute a statistic).
4. Pass placement of V12, V14, V15: the lead's brief put plots whose tables land later in pass 2.
   **Proposed default:** keep them in pass 2 but let each start the day its producing PR merges.

## Questions that need Khoa (not the lead)

1. The M9 gate itself: approval of the snapshot baselines and the report layout (M9-25 packet).
2. `giap compare`'s ρ CI is a new statistical claim; gating its coverage needs a
   `replicates_fast/full.compare` entry in the protected `tests/calibration/tolerances.yaml`
   (an ADR accepted by Khoa). Until then M9-22 records coverage without gating (M9-22 Q2).
3. Whether a colour-distance threshold for the deuteranopia check (decision 30) or a different
   minimum font size (decision 31) should become a gated style rule.

