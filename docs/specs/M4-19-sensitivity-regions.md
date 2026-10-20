# M4-19: Sensitivity framework and the `regions` family (`sensitivity.tsv`)

Milestone M4 · track R · depends on: M4-12, M4-13, M4-14, M4-18 · estimate ≈ 250 code lines (+ ≈ 250
test lines) · **can start now** on fixtures once dependencies merge · **not on the M4 gate path**
(see Q1)

## Spec sections implemented

- §8.6 quoted: "`basic` (default) reruns the context statistics, operons, GSEA and the determinants CV
  R² … under `abundance.missing ∈ {min, exclude}`. Under the `two_peptide` tier it also counts
  single-peptide genes as detected." "`full` adds UniOP thresholds 0.3 and 0.7, and Top3 instead of
  iBAQ." "Permutations. Sensitivity reruns use B = 199, since they only compare flagged sets."
  "Columns: `family`, `variant`, `n_primary`, `n_variant`, `n_shared`, `jaccard`; plus `agreement` for
  the dominant determinant block." "Families: islands (gene sets), regions (q ≤ 0.05 IDs), operon flags
  (pair IDs), enriched categories (q ≤ 0.05), determinants."
- §2.4 G9: "`tables/sensitivity.tsv`, with Jaccard stability under `min` and `exclude`".
- §7.3.2: single-peptide genes "as detected in a sensitivity rerun".
- §5.2 `output.sensitivity: basic | full | off`, `abundance.sensitivity: [min, exclude]`.

## Public interface

`src/giap/context/sensitivity.py`
```python
@dataclass(frozen=True, slots=True)
class Variant:
    name: str                                  # "missing=min", "missing=exclude", "single_peptide_detected", "top3"
    overrides: Mapping[str, object]            # config keys changed
    single_peptide_detected: bool = False

Family = Callable[[pd.DataFrame, GiapConfig, np.random.Generator], frozenset[str]]  # flagged IDs
FAMILIES: dict[str, Family]                    # M4 registers "regions"; M5–M8 add theirs

def variants(cfg: GiapConfig, tier: str) -> list[Variant]: ...
def jaccard(a: frozenset[str], b: frozenset[str]) -> float | None: ...
def sensitivity_rows(primary: Mapping[str, frozenset[str]], reruns: Mapping[tuple[str, str], frozenset[str]]) -> pd.DataFrame: ...
def regions_family(genes: pd.DataFrame, cfg: GiapConfig, rng) -> frozenset[str]: ...   # "<region_id>:<stat>" with q ≤ alpha.context
```
`src/giap/steps/sensitivity.py`: `run_sensitivity(work_dir, cfg)`: for each variant, rebuild the
abundance (M4-14 with `policy`), and for `single_peptide_detected` also the detection column and a
refitted `p_detect` (M4-13), then rerun each registered family with `stats.permutations = 199`;
writes `tables/sensitivity.tsv` (`io/sensitivity_tsv.py`, schema `src/giap/schemas/tables/sensitivity.json`).

- `variants`: `off` → none; `basic` → `missing=<p>` for each p in `abundance.sensitivity` (skipping
  the primary policy), plus `single_peptide_detected` when tier is `two_peptide`; `full` adds `top3`
  (when the source is peptide intensities) — UniOP thresholds are M6's family and are not run here.
- Primary set for comparison is recomputed with B = 199 too, so differences reflect the choice,
  not Monte Carlo error at different B (see Q2).
- `jaccard(∅, ∅)` = None → written `NA`; otherwise |A ∩ B| / |A ∪ B|.
- RNG: `giap.rng.for_step("sensitivity", seed)`; each (family, variant) draws from a child stream
  keyed by its name.

## Examples

| # | Input | Expected row |
|---|---|---|
| 1 | primary {prophage_1:o_minus_e, bgc_1:hl_shift}, `missing=min` {prophage_1:o_minus_e} | `regions`, `missing=min`, n_primary 2, n_variant 1, n_shared 1, jaccard 0.5 |
| 2 | both empty | jaccard `NA`, counts 0 |
| 3 | tier `mayu`, `output.sensitivity = basic`, primary policy median | variants `missing=min`, `missing=exclude` |
| 4 | tier `two_peptide` | adds `single_peptide_detected` |
| 5 | primary policy `min` | variants `missing=exclude` only (no rerun of the primary) |
| 6 | `output.sensitivity = off` | no file written; manifest records it |
| 7 | variant `missing=exclude` on the fixture | O − E rows identical to primary (detection unchanged); only HL rows may differ |

## Edge cases

- `exclude` with every detected gene observed → identical to primary (jaccard 1 or NA).
- `top3` variant without peptide intensities → skipped with a log line.
- Families registered by later milestones but whose module is off → no rows.

## Acceptance criteria

- [ ] `tests/unit/context/test_sensitivity.py`: examples 1–7.
- [ ] Property: `jaccard` symmetric, in [0, 1], 1 iff sets equal and non-empty.
- [ ] Integration: `run_sensitivity` on the fixture writes a schema-valid `sensitivity.tsv` with the
  `regions` family rows; determinism with the same seed.

## Methods note (statistician sign-off)

`docs/methods/sensitivity.md`: variants, B = 199, primary recomputed at the same B, Jaccard of flagged
sets, the single-peptide rerun refits detectability.

## Questions for the lead

1. Build the sensitivity framework in M4 (only the `regions` family) or defer to M5? **Proposed
   default:** M4, so M5–M8 only register families; it is off the gate path and can slip.
2. Compare variants with a primary recomputed at B = 199 (as above) or with the reported primary
   (B = 1,000)? **Proposed default:** recomputed at 199; otherwise permutation noise reads as instability.

## Lead decisions (2026-10-20, after the statistician's attack check)
1. Q1 is already decided by M4-index decision 7 (framework built in M4, regions family only, off the
   gate path).
2. Q2 accepted, with common random numbers: one child stream keyed by family (not by family × variant),
   including the M4-13 folds for the `single_peptide_detected` refit. With
   `exhaustive_shift_max_genes = 20000` every bacterial replicon uses all shifts, so the recomputed
   primary equals the reported one — test that. The methods note must say that above 20,000 genes the
   p floor is 0.005 and `n_primary` may differ from the count in `regions_stats.tsv`.
3. Family IDs are `<region_id>:<contig>:<stat>` (M4-17 `split_contigs` regions get one row per contig).
4. Flagged set = rows with `null_model = circular_shift` and q ≤ `alpha.context`; rows with p NA are
   never flagged.

## Lead decisions after the test-author's report (2026-10-20)
5. The family signature stays `(genes, cfg, rng)`; `run_sensitivity` binds C1 replicons, C7 regions and
   C8 region_genes with `functools.partial` (as M5-10 decision 38), so the regions family uses the same
   inputs and topology as `regions_stats`.
6. The single-peptide refit reuses the primary's folds (`for_step("detectability", seed)`); all other
   draws of a family come from its one family stream (CRN). The primary is never refit.
7. Example 6: no sensitivity.tsv is written; the step logs that the family is off. Recording it in the run
   manifest is done by `finalize` later (backlog).
8. A family whose module is disabled produces no rows (backlog test once module switches exist).
9. `top3` is not rerun when `abundance.method` is already `top3`; it is skipped with a "top3" log line when
   there are no peptide intensities.
10. Writer `giap.io.sensitivity_tsv`, output `tables/sensitivity.tsv`, schema
    `src/giap/schemas/tables/sensitivity.schema.json` (repo convention); the test may stop accepting
    `sensitivity.json` later.
11. The implementer starts after M4-15 (#111) merges, since the integration tests need its steps.

## Lead decisions after the code review of #117 (2026-10-20)
12. The step loads only C1 (genes), C7 (regions, region_genes) and C8 (replicons) and raises on their
    error issues (`check_issues`); a missing required table is an error, never empty rows.
13. No fallback: `regions_family` requires `regions`, `region_genes` and `replicons` (keyword-only, no
    default); the test helper `call_regions_family` passes them. Topology always comes from C8.
14. Tests (test-author) kill M09 (`top3` must read `top3_log2`), M27 (decision-5 binding: a fixture whose
    C7/C8 differ from what `region_ids` imply, including a linear replicon), M24 (refitting the unchanged
    primary reproduces its `p_detect` exactly), plus a missing-file test for each of C1/C7/C8.
15. The `top3` variant table sets `abundance_type` = top3 and carries the top3 raw values and flags.
16. Backlog: M4-19 step registration and Nextflow wiring; decision 7 (manifest records `off`); decision 8
    (no rows when a module is off).
17. Decision 12 table names use the M4-06 numbering: C1 replicons, C2 genes, C7 regions, C8 region_genes
    (the step loads the gene table plus regions, region_genes and replicons). Decision 15 reading
    confirmed: under `top3`, `abundance_raw` holds the Top3 raw value, so log2(abundance_raw) = top3_log2
    on observed genes.
18. (statistician review of #117) Methods-note fixes only: common random numbers hold only with
    exhaustive shifts (≤ `exhaustive_shift_max_genes`, every bacterial replicon); with sampled shifts a
    change in one region shifts later regions' draws (measured: O−E p moved by up to 0.075). Single-peptide
    genes keep their observed one-peptide abundance. Add an interpretation paragraph. Backlog: per-region
    shift stream in M4-17 `test_region`; skip or shrink bootstrap CIs in sensitivity reruns.
