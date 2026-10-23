# M8 "Function and determinants": issue index

Spec: GIAP_SPEC.md §13 row M8. Deliverables: "Modules 7 and 8; fgsea port; Wallenius ORA; R golden
data". **Gate:** "GSEA and ORA match fgsea and goseq on golden data; null p-values uniform", i.e.
the §7.7.5 golden tolerances (`tolerances.yaml` `golden_r`: Wallenius relative difference ≤ 10⁻⁶,
fgsea ES ≤ 10⁻⁹, NES ≤ 0.02, |log₁₀ p − log₁₀ p_R| ≤ 3 × combined `log2err` / log₂ 10) and the
§8.4 suites `ora` (R = 100 / 1,000), `gsea` (50 / 300) and `determinants` (20 / 100) at level.
Nothing in M8 is blocked; the R environment for the golden data has to be **built** (M8-12).

Sources: §7.7 (whole Module 7: §7.7.1 annotation and resources, §7.7.2 `category_summary.tsv`,
§7.7.3 Wallenius ORA and GO elim, §7.7.4 TU-collapsed fgsea multilevel and redundancy, §7.7.5
outputs and tests), §7.8 (whole Module 8: §7.8.1 features, §7.8.2 the six analysis steps and the
caveats, §7.8.3 outputs and tests), §2.1 decisions 10 and 14, §2.4 G15, §3.2 sources
(`purl.obolibrary.org`, `rest.kegg.jp`), §3.7 ("R golden references … conda R env, or a short batch
job"), §5.2 `function:` and `determinants:` keys, `stats.alpha.enrichment` / `.determinants`,
`stats.bootstrap`, §6.2 steps `function` and `determinants` (and `sensitivity`), §6.3/§6.4, §6.5,
§6.6 (LightGBM determinism, numba reduction order), §7.3 (eligible genes, `p_detect`, imputation),
§7.9.1 (`giap setup --go --kegg --cog`), §8.1 rows "Category detection", "GSEA", "Determinant
model", §8.2, §8.3 (`wallenius_ora`, `fgsea_multilevel`, `tu_block_bootstrap`), §8.4 suites `ora`,
`gsea`, `determinants`, §8.5 attacks 2 and 3, §8.6 families "enriched categories" and
"determinants", §9.1–§9.5 (the six tables, flags `abundance_not_explained_by_sequence`,
`cai_reference`, `bakta_db_light`, `layer_missing:<layer>`), §10.2 V13–V15/V18 (M9 consumers),
§11.2 golden layer, §11.3 ("Determinants. Planted additive effects …"), §11.4 (full-DB fixtures),
§11.5 (golden references from R), §12 (`golden-r.yml`), §16, §17 (fgsea fallback). ADRs read:
ADR-0006 (value encoding), ADR-0007 (size counts code), ADR-0014/0017 (KEGG is unversioned:
recorded, not pinned), ADR-0015 (two-block confirmation), ADR-0018 (detectability confounding —
the ORA attack is the same failure mode), ADR-0019, ADR-0021 Amendment 1 (no `uncalibrated` flag).

One issue = one PR; size counts **code** only (ADR-0007; `tests/sim/` lines are tooling and listed
separately). "Stats" = statistical change: label `stats`, the statistician reviews the methods note
and the simulations before the code (§14.1).

## How M8 is sequenced (Khoa, 2026-10-20/21)

"Write the spec as usual, but prioritize end to end integration first, then work on the rigour
after." So this index carries the **whole** of M8 — every table, test, golden file, calibration and
ADR the spec owes — and splits it into two passes, as M6 and MVP-01 do:

- **Pass 1 (integration):** 01–10, in two parallel tracks that each end in a wiring issue.
  - *Function track* 01, 02, 03 → 04, 05 → **06** (`function` step wired). GSEA runs on the fgsea
    *simple* stage, labelled `method = fgsea_simple`, with its honest 1/(nPermSimple + 1) floor.
  - *Determinants track* 07 → 08 → 09 → **10** (`determinants` step wired, permutation null on).
  - End state: `giap run` on sim100k writes all six M8 tables, both executors, byte-identical.
- **Pass 2 (hardening and the gate):** 11 (fgsea multilevel; GSEA switches to
  `method = fgsea_multilevel`), 12 (the R environment and the golden data), 13 (**gate, half 1**:
  golden equivalence), 14 (simulator scenarios), 15 (**gate, half 2**: the three null suites and
  attacks 2 and 3), 16 (sensitivity families), and every per-issue "Hardening owed" section.

Nothing is dropped; it is scheduled. Every issue that produces a table states how it reaches
`giap run`: 06 for the three Module 7 tables, 10 for the three Module 8 tables, 16 for the
sensitivity rows.

## How this design lands (one branch, three commits, docs only)

| Commit | Files |
|---|---|
| 1/3 | this index, M8-01, M8-02 (the two engine issues) |
| 2/3 | M8-03 … M8-10 (the rest of pass 1) |
| 3/3 | M8-11 … M8-16 (pass 2: multilevel, golden data, simulator, calibration, sensitivity) |

## Issues

| Order | Issue | Track | Pass | Scope | Depends on | Code lines | Stats |
|---|---|---|---|---|---|---|---|
| 1 | [M8-01](M8-01-stats-wallenius-elim.md) | S | 1 | `stats/`: Wallenius ORA with p_detect odds (direct tail sums), generic GO-elim walker | none | 230 | yes |
| 2 | [M8-02](M8-02-stats-fgsea-simple.md) | S | 1 | `stats/`: fgsea port I — ES, leading edge, the simple permutation stage, NES, simple `log2err` (numba) | none | 320 | yes |
| 3 | [M8-03](M8-03-function-collections.md) | R | 1 | resources → collections: GO (propagated, BP/MF/CC), GO slim, KEGG pathways, COG letters, coverage | none | 300 | no |
| 4 | [M8-04](M8-04-category-summary-ora.md) | R | 1 | `category_summary` (O − E, TU-bootstrap CI) and `enrichment_detection` (ORA, BH families, elim, `parent_group`) | 01, 03 | 300 | yes |
| 5 | [M8-05](M8-05-gsea-tu-collapsed.md) | R | 1 | `enrichment_gsea`: TU collapse, set sizes in TUs, BH per collection, leading-edge genes, `redundant` | 02, 03 | 260 | yes |
| 6 | [M8-06](M8-06-function-step-and-run.md) | T | 1 | the `function` step, three writers + schemas, `detectability`'s `tu_folds` output (ADR-0035), `function.nf`, `main.nf` (first `params.modules` switch), five setup-manifest roots, a real `giap run`; split 06a (step + tables) / 06b (both executors) | 03–05 (04b = #179) | 420 (≈ 125 JSON) | no |
| 7 | [M8-07](M8-07-determinants-features.md) | R | 1 | `determinants_features.tsv`: the five blocks, the Met rule, `dist_to_tu_start_bp`, region and context columns | none | 280 | no |
| 8 | [M8-08](M8-08-determinants-models.md) | R | 1 | univariate ρ / AUROC with TU-bootstrap CIs; LightGBM + elastic net under repeated `GroupKFold(5)` by TU; drop-one-block ΔR² | 07 | 340 | yes |
| 9 | [M8-09](M8-09-determinants-attribution-null.md) | R | 1 | SHAP block sums (`pred_contrib`), dominant block, residual flag; TU-size-preserving permutation null | 08 | 260 | yes |
| 10 | [M8-10](M8-10-determinants-step-and-run.md) | T | 1 | the `determinants` step, three writers + schemas, `determinants.nf`, `main.nf`, a real `giap run` | 07–09 | 300 | no |
| 11 | [M8-11](M8-11-stats-fgsea-multilevel.md) | S | 2 | `stats/`: fgsea port II — multilevel splitting (`sampleSize = 101`, `eps = 0`), multilevel `log2err` | 02 | 330 | yes |
| 12 | [M8-12](M8-12-golden-r-env.md) | G | 2 | R env (conda-forge + bioconda pins), frozen golden inputs, `make_golden_r.R`, `golden_r.sbatch`, `golden-r.yml`, committed outputs | 02 (inputs only) | 220 | no |
| 13 | [M8-13](M8-13-golden-equivalence.md) | C | 2 | **gate, half 1**: Wallenius vs goseq, ES / NES / p vs fgsea, on the committed golden data | 01, 12, 12b, 17 (#214) | 20 (+380 tests) | yes |
| 14 | [M8-14](M8-14-sim-function-determinants.md) | F | 2 | giapsim: synthetic ontology + xrefs, ORA null and confounding attack, GSEA null and operon-packed attack, planted determinant effects, determinants null | none | 0 (+330 sim) | yes (light) |
| 15 | [M8-15](M8-15-calibration-gate.md) | C | 2 | **gate, half 2**: `ora`, `gsea`, `determinants` at level; attacks 2 and 3; planted-R² recovery recorded | 04, 05, 09, 11, 14 | 60 (+tests) | yes |
| 16 | [M8-16](M8-16-sensitivity-function-determinants.md) | R | 2 | §8.6 families "enriched categories" and "determinants" (CV R², dominant-block `agreement`) | 06, 10 | 160 | yes (light) |

Total ≈ 3,750 code lines in `src/`, `workflow/` and `scripts/` (+ ≈ 300 simulator lines,
≈ 4,500 test lines); ≈ 2,900 of the code lines are pass 1. The spec's 1.5 weeks fits only if the
two pass-1 tracks run in parallel and 12 (an R environment build on Slurm) is started early, since
it has no code dependency beyond the frozen inputs.

## Dependency DAG

```text
Pass 1  function:      01 ─┐
                       03 ─┼─► 04 ─┐
                       02 ─┴─► 05 ─┴─► 06 ─► (giap run writes the three Module 7 tables)
        determinants:  07 ─► 08 ─► 09 ─► 10 ─► (giap run writes the three Module 8 tables)
Pass 2  02 ─► 11 ─► 05 switches to multilevel (one call site)
        12 (env build + golden data; Slurm) ─► 13 (gate half 1)  ◄─ 01, 11
        14 ─► 15 (gate half 2)  ◄─ 04, 05, 09, 11
        06, 10 ─► 16
```

01, 02, 03 and 07 have no dependency and can start at once (four independent PRs); 12's environment
build can be submitted in wave A as well, because building and freezing inputs needs no GIAP code
beyond giapsim.

## Waves

| Wave | Issues | Notes |
|---|---|---|
| A (now) | 01, 02, 03, 07; start M8-12 with its env yml + `golden_r.sbatch` | four independent PRs; the env build is a Slurm job (≤ 4 CPUs) whose result the rest of M8-12 waits for |
| B | 04, 05, 08 | 05 needs only 02's simple stage |
| C | 06, 09 | 06 is the function integration milestone (a real `sbatch scripts/slurm/giap_run.sbatch`) |
| D | 10 | the determinants integration milestone; same evidence rule |
| E (hardening) | 11, 12, 14 in parallel; then 13, 15, 16 | 15 submits `calibration_full.sbatch` tasks `ora`, `gsea`, `determinants` (≤ 3 jobs, ≤ 48 CPUs; sum running `giap-*` CPUs first, R14) |

## What M8 consumes (checked against merged main, c4d8841)

| ID | Artefact | Producer (status) | What M8 relies on |
|---|---|---|---|
| G | `genetable/genes.parquet` | M4-15 + M4-13 `detectability` (merged; wired by MVP-01) | `locus_tag`, `contig`, `start`, `end`, `strand`, `partial`, `feature_type`, `tu_id`, `tu_position`, `tu_size`, `region_ids`, `detected`, `detection_status`, `p_detect`, `abundance_log2`, `abundance_imputed`, `abundance_scale`, the layer and protein feature columns, `cog_id`, `cog_category`, `go_terms`, `kegg_ko`, `ec`, `annotation_source`, `flags` (incl. `cai_reference`) |
| C2 | `genes.parquet` | M2-04 (merged) | `aa` (the N-terminal rule), `ribosomal_protein` (CAI reference set), intervals |
| C5 | `tus.parquet` | M2-12 (merged) | TU membership and transcription order (`dist_to_tu_start_bp`, TU collapse) |
| C6 | `layers.parquet` | M2-20 (merged) | `n_promoter_windows_tu` (= `n_promoters_tu_to_gene`), `spacing_dg`/`standby_dg` (not used), the layer columns the gene table already carries |
| C7/C8 | `regions.parquet`, `region_genes.parquet` | M2-19 (merged) | `in_prophage`, `in_bgc` |
| X | `context/gene_context.parquet` | M5-09 (merged, #138) | `replichore_position` (NA without an oriC call) |
| R | `$GIAP_DB_DIR` items `go_basic`, `goslim_generic` (GO 2026-05-19, sha256-pinned), `kegg_ko_pathway`, `kegg_pathway_names` (retrieved 2026-10-15, ADR-0017), `cog2020_def`, `cog2020_fun` | M2-09 `giap setup` (merged; **all installed on FarmShare**, `db_manifest.json` status ok) | read through the setup manifest, as `cog_table()` already does (M4-15 decision 5) |
| E | `giap.stats` | merged | `bh` (`multitest`), `tu_block_bootstrap`, `percentile_ci`, `perm_p`, `poisson_binomial_sf` (fallback only) |
| T | calibration harness | M4-05 (merged) | `calibrate`, `check_level`, `check_level_pooled`; `replicates_fast/full.{ora, gsea, determinants}` and `golden_r` already exist in `tolerances.yaml` (protected, **unchanged**) |
| K | config | merged `FunctionConfig` (`go_slim`, `ora_size`, `gsea_size`, `gsea.{sample_size, eps, n_perm_simple}`), `DeterminantsConfig` (`outer_folds`, `repeats`, `permutations`, `lightgbm`, `elastic_net`), `stats.alpha.{enrichment, determinants}`, `stats.bootstrap`, `modules` (incl. `function`) | **no `defaults.yaml` change** (it is protected) |
| D | dependencies | `pyproject.toml` (merged) | goatools (BSD-2), lightgbm (MIT), scikit-learn, scipy ≥ 1.18 (`nchypergeom_wallenius`), numba — all already runtime dependencies; no lock change |
| W | `workflow/subworkflows/local/{annotate,genetable,context}.nf`, `GIAP_STEP`, `stepIn`/`kv`/`pick`, `installed()` in `main.nf`; `giap.runner.local` roots (`cli.py` stages `cog2020_def`) | MVP-01 (merged) | 06 and 10 copy the CONTEXT pattern; 06 adds five setup-manifest roots the way `cogTable()` does |

M6-01 (engine, landing) is **not** needed by M8. M6-08's `operons.nf` is not on main yet; M8's
subworkflows do not depend on it, but 06 and 10 must merge `origin/main` before editing `main.nf`
if M6-08 lands first (one conflict hunk, expected).

## Licences, and how the R golden data is produced

| Package | Version (bioconda / conda-forge, 2026-10-21) | Licence (official source) | Role in GIAP |
|---|---|---|---|
| fgsea | 1.36.2 (Bioconductor release is 1.38.0; bioconda lags) | **MIT + file LICENCE** (bioconductor.org/packages/fgsea; LICENCE: "2016-2019, Alexey Sergushichev") | **ported** to Python + numba from its source (G15). MIT permits this; the port carries fgsea's copyright and permission notice in the module header and in NOTICE.md. Never a runtime dependency |
| goseq | 1.62.0 (Bioconductor release 1.64.0) | **LGPL (≥ 2)** | **not ported, not read.** GIAP's Wallenius ORA is written from the spec formula and the goseq paper (Young 2010) on top of scipy's `nchypergeom_wallenius` (BSD). goseq runs only in the golden environment |
| BiasedUrn (goseq's Wallenius) | 2.0.12 | **GPL-3** | never imported or copied; golden environment only |
| r-base | 4.5.x (the r45 builds of both Bioconductor packages pin `r-base >=4.5,<4.6`) | GPL-2+ | golden environment only |

So the licence design is: **one MIT port (fgsea) with attribution, one clean-room implementation
(Wallenius ORA), and an R environment that exists only to generate numbers.** Numbers produced by
GPL/LGPL programs are not covered by those licences, so committing the golden outputs is fine. No
`src/` file may `import` or shell out to R; `scripts/ci/check_imports.py` gains `rpy2` in its
FORBIDDEN set (M8-12) so an accidental bridge fails CI.

**Is the R environment available on FarmShare?** No. `local/micromamba/envs/` holds `giap-ostir`,
`giap-promotech`, `giap-promoter_calculator`, `giap-terminatornet`, `giap-uniop` and `giap-wf`
only. The Lmod module `r/4.5.2` exists, but it has no Bioconductor packages, and installing them
from bioconductor.org or CRAN would use sources outside §3.2 (R6). §11.5 names the route:
conda-forge and bioconda. M8-12 therefore adds `scripts/golden/golden_r.yml` (env `giap-golden-r`)
and builds it in a Slurm job; the env is an estimated 1–2 GB (GO.db and geneLenDataBase are goseq
dependencies), well inside the §3.9 budget. Exact procedure, pins and provenance: M8-12.

## KEGG

`giap setup --kegg` (M2-09, merged) fetches `https://rest.kegg.jp/link/pathway/ko` and
`/list/pathway` with no pinned checksum (ADR-0017: retrieval date + sha256 of the bytes recorded in
`db_manifest.json`; both are installed on FarmShare, retrieved 2026-10-15, 49,071 KO→`path:map`
links plus the same number of `path:ko` duplicates, 587 pathway names). KEGG's terms (rest.kegg.jp,
checked 2026-10-21): "KEGG API at rest.kegg.jp is made available only for academic use by academic
users"; non-academic use needs a licence. Consequences for M8, all already consistent with R11 and
NOTICE.md:

- **No KEGG file or KEGG-derived table is ever committed** — not in `tests/data/`, not in the
  golden data. Every KEGG test uses a synthetic two-pathway table written by the test itself.
- The golden inputs (M8-12) use **synthetic gene sets** only, so they carry no KEGG, GO or COG
  content at all; the Wallenius/fgsea arithmetic does not care what a set is called.
- Run outputs contain KEGG pathway names for the user's own run; the report and the user guide
  state the academic-use restriction and the retrieval date (ADR-0017), and the manifest lists it.
- M8 needs no new KEGG endpoint. It uses only `path:map…` IDs (the `path:ko…` rows are exact
  duplicates) so that IDs join `list/pathway` (decision 6).

## What happens when an estimator does not exist (summary; details per issue)

| Estimator | Does not exist when | Output |
|---|---|---|
| a collection | its resource is not installed (`giap setup --go/--kegg/--cog` not run) | no rows for it; run flag `layer_missing:go` / `layer_missing:kegg` / `layer_missing:cog` |
| Wallenius ω | the category is the whole universe, or Σ p_detect outside it is 0 | the category is not tested; row `flags` `unestimated` (index decision 12) |
| ORA p | category size (full set) outside `function.ora_size` | not tested, no row; an elim-reduced set of an in-range term is still tested (M8-04 Q3) |
| any GSEA row | no abundance input, or < 2 ranked TUs | `enrichment_gsea.tsv` header-only; run flag `layer_missing:abundance` or `too_few_ranked_tus` |
| GSEA p | fgsea's rule "modeFraction < 10" (unbalanced ES signs) | `p_value`, `q_value`, `nes`, `log2err` NA, `flags` `unestimated` |
| GSEA `log2err` | multilevel with `isCpGeHalf = FALSE` | `log2err` NA, p kept, `flags` `log2err_unestimated` |
| determinants models | < 100 model rows or < 25 TUs (decision 18) | no CV / null / block rows; univariate rows still written; model table row `flags` `too_few_rows` |
| a feature | its layer is missing (no Promotech, no OSTIR …) | the column NA; the feature dropped from both models and listed in the manifest; a block with no feature left is not reported |
| dominant block | no block with a positive SHAP sum | `dominant_block` NA |

No silent fallbacks (the M4-17 rule): every NA carries a flag, and no NA p enters a BH family.

## Calibration suites owned by M8 (replicates from `tolerances.yaml`, unchanged)

| Suite (key) | Issue | Generator (M8-14 scenario) | Criterion | R fast / full |
|---|---|---|---|---|
| `ora` (**gate**) | M8-15 | `ora_null` — detection ~ Bernoulli(p_detect) only, categories independent of detection | level, pooled over categories and α ∈ {0.01, 0.05, 0.10} | 100 / 1,000 |
| attack 2 (gated) | M8-15 | `ora_attack_confounded` — a membrane-like category with low p_detect | level **with** weights; **fails without** them (teeth) | as `ora` |
| `gsea` (**gate**) | M8-15 | `gsea_null` — abundance independent of sets | level, pooled | 50 / 300 |
| attack 3 (gated) | M8-15 | `gsea_attack_operon_packed` — sets packed into single operons | level after TU collapse; gene-level ranking shown anti-conservative (teeth) | as `gsea` |
| `determinants` (**gate**) | M8-15 | `determinants_null` — response independent of features, TU structure kept | level (discrete permutation p, §8.4) | 20 / 100 |
| recovery (recorded) | M8-15 | `determinants_planted` — known block effects and true R² | CV R² within its CI of the truth; dominant block recovered; recorded | 20 |

## Methods notes (statistician sign-off, §14.1, §14.9)

New: `docs/methods/function.md` (the universe, ω and why it is goseq's weighting applied to
detectability, the direct tail sums, elim and its raw-p threshold, BH families, TU collapse and why
gene-label permutation is too liberal, the centring of TU scores, redundancy, coverage bias),
`docs/methods/fgsea-port.md` (what was ported from which fgsea version and file, the RNG difference
from R and why the golden tolerances are Monte-Carlo-aware), `docs/methods/determinants.md`
(features and their sources, the grouped CV, the null and why one CV repeat per permutation is
conservative, SHAP block sums, the caveats of §7.8.2 verbatim). Updates: `multiple-testing.md` (the
ORA and GSEA families), `calibration.md` (three suites, two attacks), `sensitivity.md` (two
families).

## ADRs M8 owes (numbers: next free at merge time; M6 reserved ADR-0022)

- **ADR-A "M8 tables carry a `flags` column"** (06, 10; decision 12): §9.3 gives
  `enrichment_detection.tsv`, `enrichment_gsea.tsv` and `determinants_model.tsv` no `flags`
  column, but R4 needs `unestimated` somewhere. Adding a column is non-breaking (§9.2 semver).
  Lead-level (fills a gap; no default changes).
- **ADR-B "GSEA TU scores are median-centred"** (05; decision 8): changes the fgsea weights, so it
  changes a statistic; statistician evidence from M8-15 (level and power, centred vs raw).
- **ADR-C "KEGG overview maps excluded"** (03; decision 6): only if the lead accepts the default.
- Conditional: if M8-13 misses a golden tolerance after 3 iterations (R9), the §7.7.4 fallback ADR
  ("run R fgsea in its own environment as an external tool") goes to Khoa; if M8-15's attack 2
  fails, an ADR on the ORA null (Poisson-binomial on Σ d given individual p_detect, the engine's
  `poisson_binomial_sf`) — never a tolerance change (R3).

## Backlog items folded in (docs/backlog.md, not duplicated)

| Backlog line | Where it lands |
|---|---|
| "M8 function and determinants" | this index; the checkbox closes when 13 and 15 are green |
| "M4-15/M4-17/M8 attack-2 suites must gate on M4-13 p_detect …" | M8-15: attack 2 uses the **fitted** cross-fitted p_detect from M4-13, never the simulator's truth; the truth arm is recorded beside it |
| "M4 notes: median imputation pulls detected-but-unquantified genes toward the centre …" | M8-05 decision 9 (imputed values tie in GSEA) and M8-16 (the `min` / `exclude` reruns bracket it) |
| "At the next uv lock update: scikit-learn>=1.8 …" | M8-08 uses `ElasticNetCV(l1_ratio=…)`; it must work on the **locked** version, and the lock is not touched by M8 |

## Risks

- **scipy's `sf` underflows.** Verified in the session on scipy 1.18.1:
  `nchypergeom_wallenius.sf(k − 1, …)` returns **0.0** for a tail of 3.2 × 10⁻¹⁴ (the tail sum is
  computed as 1 − cdf). Taken literally, §7.7.3 would give p = 0 for every strong enrichment and
  fail the 10⁻⁶ relative tolerance. M8-01 computes both tails as direct pmf sums (decision 1); the
  semantics are exactly the spec's.
- **Wallenius with a mean odds is an approximation.** It is exact only when p_detect is constant
  inside and outside C. A category with heterogeneous p_detect (attack 2) is the risk; M8-15
  measures it with the fitted p_detect and the ADR fallback is named above.
- **The golden comparison is stochastic on one side.** fgsea's NES and p depend on R's RNG; the
  port cannot reproduce R's stream. The tolerances (NES 0.02, 3 × combined `log2err`) were chosen
  for that, but a set with a large `log2err` on both sides makes the check weak. M8-13 reports how
  many sets fall in each p decade, so a vacuous pass is visible.
- **Determinants cost.** B = 1,000 permutations × (5 LightGBM + 5 elastic-net fits) on ~2,500
  rows is ~10⁴ fits per model — an estimated 30–40 CPU-minutes per run at the §5.2 parameters
  (to be measured in M8-09). The step needs
  `process_high` and parallel permutations with spawned seeds (M8-09); the calibration fast tier
  cannot afford the defaults (M8-15 decision 23).
- **Small integration fixture.** sim100k has ~100 genes, so determinants models are below the
  100-row minimum and most GSEA sets are below 10 TUs: pass 1's acceptance is "the tables exist
  with the right schema and the right `unestimated` states", not "results are non-empty". A
  meaningful real-data run (K-12 full-DB Bakta recording + in-silico peptides) is hardening.
- **Coverage leans toward well-studied genes** (§7.7.1): K-12's full-DB fixture has GO on 90 %,
  COG on 57 % and KEGG on 38 % of CDS (counted in the session). The report must print coverage on
  every row, and the universe question (decision 2) matters most for KEGG.

## Feasibility on FarmShare

Module 7 needs < 2 GB and seconds to minutes (go-basic.obo is 32 MB; GO propagation over ~4,500
genes is a set union; fgsea simple on ~2,500 TUs × ~1,500 sets is numba-bound). Module 8's null
dominates (above). The R env build is one ≤ 4-CPU job of ≤ 1 h; generating golden data is seconds.
Full calibration: `ora` (1,000 replicates, analytical p) and `gsea` (300) are cheap;
`determinants` (100 replicates × B = 999 at the §5.2 parameters, 1,500 genes; the fast tier runs
reduced, decision 23) is the only suite that needs sharding inside the R14 caps.

## Decisions taken in these specs (for the lead to confirm)

1. **Wallenius tails by direct pmf sums** (`Σ_{j≥k} pmf(j)` and `Σ_{j≤k} pmf(j)`), not
   `.sf(k − 1)` / `.cdf(k)` (M8-01 Q1). Same quantity; no cancellation.
2. **The ORA universe is the spec's: all eligible genes**, annotated or not, per collection; the
   golden goseq call therefore uses `use_genes_without_cat = TRUE` (M8-04 Q1).
3. **GO propagation uses `is_a` and `part_of`** (the GO.db / topGO convention), via goatools
   (M8-03 Q1).
4. **elim applies to `direction = over` only**; `under` rows report `p_value = p_classic`
   (topGO's elim is defined for over-representation) (M8-04 Q2).
5. **GO slim is its own collection `go_slim`**, tested with classic p (no elim) in ORA and GSEA
   (M8-03 Q3).
6. **KEGG uses `map` pathway IDs only and excludes the global/overview maps `map011xx` and
   `map012xx`** (M8-03 Q2; ADR-C).
7. **Collection values** in the enrichment tables: `cog`, `kegg`, `go_bp`, `go_mf`, `go_cc`,
   `go_slim` (the six GSEA collections of §7.7.4, reused for ORA); `category_summary.tsv` keeps
   §9.3's enum {cog, kegg, go_slim} (M8-04 Q4).
8. **TU scores are median-centred before fgsea** (stats are otherwise all positive under
   `scoreType = "std"`, where fgsea warns) (M8-05 Q1; ADR-B).
9. **Imputed abundances stay in the TU means** (the §8.6 `min`/`exclude` reruns imply they are in
   the primary analysis); ties are ordered by TU ID (M8-05 Q2).
10. **Redundancy re-test:** the parent's re-test p replaces its p in the collection's BH vector and
    the resulting q is compared with α (M8-05 Q3).
11. **`parent_group`** = the most general significant ancestor in the same namespace (the term
    itself when it has none); NA for non-significant rows and non-GO collections (M8-04 Q5).
12. **A `flags` column** on `enrichment_detection.tsv`, `enrichment_gsea.tsv` and
    `determinants_model.tsv`, with new values `unestimated`, `size_filtered`,
    `log2err_unestimated`, `too_few_rows`, `null_degenerate`, `ci_few_units` (reused from M4-16),
    added through the schema (ADR-A).
13. **Pass-1 GSEA ships the fgsea simple stage** (`method = fgsea_simple`, p floor 1/1001); M8-11
    switches the call to multilevel (M8-05 Q4).
14. **Elastic-net missing values:** median imputation fitted inside each training fold; an all-NA
    feature is dropped from both models (M8-08 Q1).
15. **CV R² per repeat = 1 − SSE/SST of the pooled out-of-fold predictions**; mean and SD over the
    5 repeats (M8-08 Q2).
16. **The null reuses the repeat-1 fold assignment** for every permutation (M8-09 Q1).
17. **The residual flag is two-sided:** |residual| > 2 SD (M8-09 Q2).
18. **Models need ≥ 100 rows and ≥ 25 TUs** (≥ 5 TUs per outer fold) (M8-08 Q3).
19. **`nterm_residue` is a feature only.** §16 excludes a "Met removal" *catalogue* for mapping and
    annotation; §7.8.1 asks for this one rule as a covariate, and nothing else uses it (M8-07 Q1).
20. **`abundance_not_explained_by_sequence`** lives in `determinants_attribution.tsv` now;
    `tables/genes.tsv` gains it when the §6.2 `finalize` step exists (`genes.tsv` is written by
    `context_stats` today, before `determinants` runs) (M8-10 Q1).
21. **The golden R env is `scripts/golden/golden_r.yml`**, built by `scripts/slurm/golden_r.sbatch`,
    not a `workflow/envs/` tool (it is not a pipeline tool and must not join `build_envs.sh all`)
    (M8-12 Q1).
22. **`golden-r.yml` (manual GitHub workflow) is created** as §12 lists it; the FarmShare job is
    the primary route (M8-12 Q2).
23. **Determinants calibration fast tier runs at reduced model size** (400 genes, B = 99,
    `n_estimators` 50); the full tier uses the §5.2 parameters. `tolerances.yaml` R is unchanged
    (M8-15 Q1).
24. **"Null p-values uniform" is read as §8.4's level criterion** for all three suites (their p are
    discrete: hypergeometric-type, permutation, multilevel-from-permutation); KS is recorded, not
    gated (M8-15 Q2).
25. **Sensitivity:** the "enriched categories" family reruns GSEA under `min` / `exclude` and ORA
    under the `two_peptide` variant only (M8-16 Q1).
26. **RNG streams:** `function.gsea` (fgsea), `function.bootstrap` (category CIs),
    `determinants.cv`, `determinants.null`, `determinants.bootstrap`, each from
    `giap.rng.for_step`; the ORA and elim are deterministic.

## Questions that need Khoa (not the lead)

1. If the fgsea port misses a golden tolerance after 3 iterations (R9), the §7.7.4 fallback ADR
   (R fgsea as an external tool in its own env) is a dependency change and goes to Khoa.
2. ADR-B (centring) changes a statistic; the lead may accept it, but if the lead prefers the
   literal spec (uncentred scores), Khoa should know the fgsea warning will appear on every run.

## Lead decisions (2026-10-21)
1. The 26 decisions above are **accepted** as proposed, together with each issue's smaller proposed defaults.
2. **ADR-B (median-centring the TU scores before fgsea)** is accepted by the lead. M8-05 writes it; Khoa is informed in the PR. The literal uncentred spec would raise the fgsea warning on every run.
3. **Licences:**
   - fgsea is MIT: port it with attribution in each file header and in NOTICE.md.
   - goseq (LGPL) and BiasedUrn (GPL-3) are never ported, read or imported. R runs only to produce golden numbers, in its own env (`giap-golden-r`, M8-12).
   - The guard against `rpy2` and `Rscript` lands with M8-12.
4. **KEGG:** no KEGG file or KEGG-derived table is ever committed (KEGG academic terms; ADR-0017). Tests and golden data use synthetic sets.
5. **The scipy Wallenius tail bug:** M8-01 sums the pmf directly. The `sf(k−1) == 0` case at 3.2e-14 is pinned as a regression test.
6. **fgsea version:** port from the **1.36.2** sources, the bioconda version that the golden env pins. The implementer diffs master's details (integer rounding, sign-conditional p, reversed negative leading edge) against 1.36.2 before relying on them.
7. **For Khoa:** the §7.7.4 fallback ADR (R fgsea as an external tool) is raised only if M8-13 misses a tolerance three times (R9).
