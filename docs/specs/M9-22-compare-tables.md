# M9-22: `giap compare`: genome check and the three comparison tables

Milestone M9 · **pass 2** · depends on: none (reads `OUT/` trees; uses `giap.stats`) · estimate
≈ 380 code lines (+ ≈ 380 test lines) · **stats** (Spearman ρ with a TU-block bootstrap CI,
Bland–Altman limits, Jaccard indices) · writes **ADR-M9-B**

## Spec sections implemented

- §7.9.3: "Inputs. Two or more run directories on the same genome, checked through the genome
  SHA-256 in each manifest; a mismatch is exit code 1. Writes to `CMP/`: `tables/compare_detection.tsv`:
  a gene × run detection matrix plus the UpSet intersections; `compare_abundance.tsv`: pairwise
  Spearman ρ on shared detected genes with TU-bootstrap CIs, and Bland–Altman bias and limits;
  `compare_stability.tsv`: Jaccard overlap of island gene sets, operon-flag pairs, enriched
  categories (q ≤ 0.05), novel ORFs and ALT sites; `plots/V16_compare.*`; `report/compare.html`."
  (V16 and the HTML are M9-23.)
- §5.1 `giap compare OUT_A OUT_B [OUT_C ...] --outdir CMP [--labels A1,A2,A3]`, exit codes; §8.2
  (CIs, R4), §8.3 (`tu_block_bootstrap`, `percentile_ci`), §5.2 `stats.bootstrap`,
  `stats.alpha.enrichment`, `stats.alpha.context`; §6.6 (`giap.rng.for_step("compare", seed)`);
  §9.2 conventions; §9.5 schemas; §15 (M10 runs A1/A2/A3 through this command).

## ADR-M9-B (written by this issue): `giap compare` output layout

1. All three tables live under `CMP/tables/` (§7.9.3 places only the first there explicitly).
2. `compare_detection.tsv` is **long**: `section` enum{matrix, intersection}. Per-run columns would
   make the schema depend on the run count.
3. Genes are joined by `locus_tag`; every run must have the same locus-tag set (exit 1 otherwise).
4. Genome identity: the manifest's genome SHA-256 when `run_manifest.json` exists, else
   `preflight/preflight.json` `inputs[role = genome].sha256` (index decision 26).

## Public interface

`src/giap/compare/` (pure): `detection.py`, `abundance.py`, `stability.py`; `src/giap/io/compare_io.py`
(read runs, write tables); schemas `src/giap/schemas/tables/compare_{detection,abundance,stability}.schema.json`;
`cli.compare` replaces `_not_implemented`. Methods note `docs/methods/compare.md`.

```python
def detection_matrix(runs: Mapping[str, pd.DataFrame]) -> pd.DataFrame
def intersections(matrix: pd.DataFrame, labels: Sequence[str]) -> pd.DataFrame
def abundance_pairs(runs: Mapping[str, pd.DataFrame], tus: pd.Series, B: int,
                    rng: np.random.Generator) -> pd.DataFrame
def stability(runs: Mapping[str, RunView], labels: Sequence[str], cfg: GiapConfig) -> pd.DataFrame
```

## Tables

**`compare_detection.tsv`**: `section`, `locus_tag`, `run`, `detected` (bool), `detection_status`,
`intersection` (runs joined by `&` in label order), `n_genes`. Matrix rows: one per gene × run
(`intersection`, `n_genes` NA). Intersection rows: one per non-empty exact intersection of the
detected sets (a gene detected in exactly A1 and A3 belongs to `A1&A3` only), `locus_tag`, `run`,
`detected` NA; plus `intersection = none` for genes detected in no run. Sorted by section, then
genes in §6.5 order, then run in label order; intersections by `n_genes` descending then name.

**`compare_abundance.tsv`**: `run_a`, `run_b`, `n_shared`, `spearman_rho`, `ci_low`, `ci_high`,
`ba_bias_log2`, `ba_loa_low`, `ba_loa_high`, `method` (`spearman_tu_bootstrap`), `n_boot`,
`flags`. One row per unordered pair (label order). Shared genes = detected in both with
`abundance_imputed = false` in both. ρ = Spearman (average ranks for ties); CI = percentile CI of
B = `stats.bootstrap` TU-block resamples (the TUs of the first run of the pair, resampled with
replacement; genes carried with their TU); Bland–Altman on log2 abundance: bias = mean(a − b),
limits = bias ± 1.96 × SD(a − b) (ddof = 1).

**`compare_stability.tsv`**: `family` enum{islands, operon_flags, enriched_categories, novel_orfs,
alt_sites}, `run_a`, `run_b`, `n_a`, `n_b`, `n_shared`, `jaccard`, `flags`.

| Family | Set per run |
|---|---|
| islands | genes in islands with `p_value ≤ stats.alpha.context` |
| operon_flags | (`gene_a`, `gene_b`) of `operon_flags` rows with `q_value ≤ stats.alpha.operon_flags` |
| enriched_categories | (`table`, `collection`, `term`, `direction`) with `q_value ≤ stats.alpha.enrichment` in `enrichment_detection` and `enrichment_gsea` (direction = sign of NES) |
| novel_orfs | `orf_id` with `called = true` |
| alt_sites | (`contig`, `codon_start`, `strand`, `alt_aa`) of `variants.tsv` rows with `call = ALT` |

Jaccard = n_shared / |A ∪ B|; both empty → NA, flag `both_empty`; a table absent in a run → NA,
flag `missing_in:<label>`.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| one run given | exit 1 (§7.9.3 "two or more") |
| genome SHA-256s differ | exit 1 naming both labels and hashes |
| neither manifest nor preflight hash | exit 1 ("cannot verify the genome") |
| locus-tag sets differ | exit 1 with counts in each run and the symmetric difference size |
| `--labels` count ≠ runs, or duplicate labels | exit 1 |
| no `--labels` | directory basenames; duplicates suffixed `_2`, `_3` |
| `n_shared` < 3 | ρ, CI, BA NA, flag `too_few_shared` |
| `n_shared` TUs < 2 | CI NA, flag `ci_few_units` (the M4-16 flag) |
| abundance types or scales differ between runs (`abundance_type`, `abundance_scale`) | BA columns written, flag `scale_differs` (bias is not interpretable) |
| a run without abundance | that pair's abundance row NA, flag `missing_in:<label>` |
| `function` module off in one run | `enriched_categories` NA, `missing_in` |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | A1 detects {g1, g2, g3}, A2 {g2, g3, g4}, A3 {g3} on genes g1–g5 | intersections: `A1&A2&A3` → g3 (1); `A1` → 1; `A1&A2` → g2 (1); `A2` → 1; `none` → g5 (1) |
| 2 | shared measured log2: A = 1, 2, 3, 4, 5; B = 2, 1, 4, 3, 5 | ρ = 0.8; bias 0; SD 1; limits −1.96, 1.96 |
| 3 | island genes A {g1, g2, g3}, B {g2, g3, g4} | jaccard 0.5, n_shared 2 |
| 4 | both runs no novel ORF called | jaccard NA, flag `both_empty` |
| 5 | genome hashes differ | exit 1 |
| 6 | same seed, same runs, twice | byte-identical tables (bootstrap from `for_step("compare", seed)`) |
| 7 | A1 `maxlfq`, A2 `ibaq` | flag `scale_differs` |

## Acceptance criteria

- [ ] `tests/unit/compare/test_detection.py`, `test_abundance.py`, `test_stability.py`: examples
  1–4, 7; ties in Spearman; ddof.
- [ ] `tests/unit/test_cli_compare.py`: examples 5, 6 and every exit-1 row.
- [ ] Calibration (statistician): `tests/calibration/test_compare_rho_ci.py` — coverage of the
  TU-block bootstrap CI for ρ on giapsim pairs with a known common abundance and TU-correlated
  noise; R from a new `tolerances.yaml` entry **proposed** to Khoa (protected file; Q2), recorded
  until then.
- [ ] Schemas validate; ADR-M9-B committed; `docs/methods/compare.md` signed off.

## Questions for the lead

1. Genome check through `preflight.json` until the manifest exists (ADR-M9-B item 4). **Proposed
   default:** yes.
2. The ρ CI is a new statistical claim needing calibration; its replicate count belongs in the
   protected `tolerances.yaml`. **Proposed default:** a Khoa-approved ADR adds
   `replicates_fast/full.compare` (e.g. 100 / 1,000 like `ora`); until then the test records
   coverage and does not gate (the ADR-0032 pattern).
3. "q ≤ 0.05" in §7.9.3 vs the config key: **Proposed default:** `stats.alpha.enrichment` (0.05 by
   default), so a user's alpha applies.
4. Exact (UpSet) intersections vs inclusive ones: **Proposed default:** exact, as UpSet plots draw.

## Lead decisions (2026-10-22)

1. The genome check reads `run_manifest.json` (written since MVP-03a, #196) when present and falls
   back to `preflight.json` otherwise; both runs must agree on the genome digest.
2. The ρ CI's coverage is recorded, not gated, until a Khoa-approved ADR adds
   `replicates_fast/full.compare` to `tolerances.yaml` (the ADR-0032 pattern; raised with Khoa).
3. The threshold is the config key `stats.alpha.enrichment` (0.05 by default).
4. Intersections are exact (UpSet semantics).
5. Feature loop: test-author → implementer → statistician (stats label) → code review.
