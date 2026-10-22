# M8-05: Abundance GSEA — TU-collapsed ranking, set sizes in TUs, BH per collection, `redundant`

Milestone M8 · track R · **pass 1** · depends on: M8-02 (fgsea simple), M8-03 (collections) ·
estimate ≈ 260 code lines (+ ≈ 350 test lines) · **statistical** (label `stats`)

## Spec sections implemented

- §7.7.4 "Ranking unit. Genes in one TU share a transcript, so permuting gene labels would be too
  liberal. The ranked units are TUs: each TU's score is the mean log2 abundance of its detected
  genes, and a TU belongs to a set when any of its detected genes does. Set sizes are counted in
  TUs and limited to `function.gsea_size` ([10, 500])."
- §7.7.4 "Reported per set: ES, NES, p, `log2err`, BH q within each collection (COG, KEGG, GO BP,
  GO MF, GO CC, GO slim), and the leading edge (genes)."
- §7.7.4 "Redundancy. A significant GO parent is flagged `redundant` when re-testing it without the
  genes of its significant descendants gives q > α."
- §8.1 row "GSEA | gene set | TU-collapsed scores, fgsea multilevel | BH per collection | 0.05";
  §8.2 "Resample TUs, not genes"; §8.5 attack 3; §9.3 `enrichment_gsea.tsv`.

## Public interface

`src/giap/function/gsea.py` (pure)
```python
@dataclass(frozen=True, slots=True)
class TuRanking:
    tu_ids: tuple[str, ...]; scores: F64            # centred (decision 8)
    tu_genes: Mapping[str, tuple[str, ...]]         # detected member genes per TU

def tu_ranking(genes: pd.DataFrame) -> TuRanking: ...
def tu_sets(rank: TuRanking, coll: GeneSets) -> dict[str, frozenset[str]]: ...
def gsea_table(genes: pd.DataFrame, colls: Sequence[GeneSets], size: tuple[int, int],
               cfg: GseaConfig, alpha: float, rng: np.random.Generator) -> pd.DataFrame: ...
```

1. **TU scores.** Detected eligible genes with a non-NA `abundance_log2` (imputed included —
   decision 9) are grouped by `tu_id`; score = their mean. TUs without such a gene are not ranked.
   Scores are **median-centred** (score − median over ranked TUs; decision 8, ADR-B) so fgsea's
   `std` mode sees both signs; ranks are unchanged by centring, only the weights are.
2. **Order and ties.** The input to `prepare` is ordered by TU ID (natural sort of `TU_<n>`), so
   tied scores (frequent under median imputation) are broken by TU ID, identically in R (M8-12
   writes the same order).
3. **Set membership in TUs.** A TU is in set S when any of its detected genes is in S (the gene-
   level `GeneSets`). Size = number of member TUs; sets outside `function.gsea_size` are dropped.
4. **fgsea.** Per collection, `prepare` + `fgsea_simple` (M8-02) with `n_perm =
   cfg.n_perm_simple`, `sample_size = cfg.sample_size`; pass 1 reports the simple p and
   `method = fgsea_simple` (decision 13); M8-11 adds the multilevel call for `needs_multilevel`
   sets and `method = fgsea_multilevel`. `eps = cfg.eps` (0) is passed through to M8-11.
5. **Columns.** `size_tus`, `es`, `nes`, `log2err`; `stat = es`, `effect = nes`, `ci_low/ci_high`
   NA (fgsea has no CI; the methods note says so and `log2err` is the precision measure);
   `p_value`; `q_value` = BH within the collection (§7.7.4 — one family per collection, both ES
   signs together, as fgsea's `padj`); `null_model = tu_set_permutation`; `n_perm =
   cfg.n_perm_simple`.
6. **Leading edge (genes).** The leading-edge TUs from M8-02 expanded to the detected genes of
   those TUs that are members of S, in leading-edge TU order, then gene order within the TU;
   `;`-joined in the writer.
7. **Redundancy** (GO namespaces): for each set with q ≤ α that has ≥ 1 significant descendant,
   remove the genes of all its significant descendants, rebuild its TU membership (a TU stays when
   another of its detected genes is still in the reduced set), re-run fgsea on that one set with
   the same `n_perm` (stream `function.gsea`, after the main runs), substitute its p into the
   collection's BH vector and recompute **its** q (decision 10). `redundant = true` when that q
   > α. Every other row: `redundant = false`; non-GO collections: `false`.
8. **`parent_group`**, **`coverage`**: as M8-04 (most general significant ancestor; the
   collection's coverage).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no abundance input, or `abundance_scale` absent | no rows; run flag `layer_missing:abundance` (M8-06 writes the header-only file) |
| < 2 ranked TUs | no rows; run flag `too_few_ranked_tus` |
| fgsea `modeFraction < 10` | p, q, NES, `log2err` NA; `flags = unestimated` (index decision 12) |
| NES undefined (same-sign mean 0) | as above |
| ES = 0 tie | `leading_edge` empty (`NA`), row otherwise normal |
| the redundancy re-test drops the set below 10 TUs | re-tested anyway (the filter applies to the original set); recorded in `flags` as `size_filtered` |
| `relative` abundance scale | GSEA still runs (ranks are scale-free; §7.3.4 restricts only Module 4's ratio tests) |
| sim100k | almost no set reaches 10 TUs: the table is header-only or nearly — a **passing** pass-1 result |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | TU_1 {g1 18.0, g2 16.0 (detected)}, TU_2 {g3 15.0}, TU_3 {g4 undetected} | TU scores TU_1 17.0, TU_2 15.0; TU_3 not ranked; median 16.0 → centred **+1.0, −1.0** |
| 2 | set S = {g2} | TU_1 ∈ S (g2 is detected); size 1 TU |
| 3 | set S = {g4} | no TU (g4 undetected): size 0, dropped |
| 4 | two TUs tied at the median after imputation, IDs TU_10 and TU_9 | order TU_9 before TU_10 (natural sort) |
| 5 | a 600-TU set | dropped (`gsea_size` [10, 500]) |
| 6 | GO parent P (q 0.01) with significant child C; P's TUs are exactly C's TUs plus 3 others | P re-tested on the 3 others (size 3), q recomputed with P's new p → if > 0.05, `redundant = true` |
| 7 | attack fixture: 30 genes of one operon in one set, gene-level vs TU-level ranking | TU-level set size 1 TU → dropped; gene-level size 30 → tested; documents why TUs (M8-15 measures the level) |
| 8 | same inputs and seed twice; 1 vs 8 threads | identical tables |

## Acceptance criteria

- [ ] `tests/unit/function/test_gsea.py` covers examples 1–8 and every degenerate row.
- [ ] Property: permuting gene rows within the gene table changes nothing; permuting TU labels of
  the scores changes nothing; adding an undetected gene to a set changes nothing.
- [ ] Property: every leading-edge gene is a detected member of the set.
- [ ] `docs/methods/function.md` §"GSEA": TU collapse, centring and why (ADR-B), the family, the
  redundancy test, what `log2err` means.
- [ ] ADR-B written (proposed) with the statistician's evidence plan (M8-15 compares centred and
  uncentred level and power).

## Calibration or attack work owed

- `gsea` suite and attack 3 "gene sets packed into one operon" (M8-15): stays at level after TU
  collapse; the gene-label ranking arm must be shown anti-conservative (teeth).

## Hardening owed (pass 2)

- Switch to multilevel for `needs_multilevel` sets once M8-11 merges (one call site, ≤ 20 lines,
  in M8-11's PR).
- Golden: the TU ranking written by M8-12 must be byte-identical to this function's output on the
  golden gene frame, so the R and Python runs see the same vector.

## Questions for the lead

1. Median-centre the TU scores (proposed; ADR-B) or rank the raw mean log2 abundances (all
   positive; fgsea warns "maybe you should switch to scoreType = pos")?
2. Imputed abundances in the TU mean: included (proposed; §8.6's `min`/`exclude` reruns imply it) or
   excluded as in Module 8's model rows?
3. Redundancy re-test: substitute the re-test p into the collection's BH vector and recompute that
   q (proposed), or BH the re-tests as their own family?
4. Pass 1 ships `method = fgsea_simple` with the 1/1001 floor (proposed), or GSEA waits for M8-11?

## Lead decisions (2026-10-21)

1. **Q1:** median-centre the TU scores, as proposed. With raw abundances every score is positive and the test would be one-sided under `scoreType = std`.
   - This is a methods choice and goes to the statistician review.
   - ADR-B is owed to Khoa. Until then, results carry no calibrated claim beyond the M8-02 p-values.
2. **Q2:** imputed abundances count in the TU mean, as proposed. The §8.6 `min`/`exclude` sensitivity reruns are where this is varied.
3. **Q3:** a re-test p replaces its value in the collection's BH vector, and that q is recomputed. Accepted as proposed.
4. **Q4:** pass 1 ships `method = fgsea_simple` with the 1/1001 floor (MVP first), and M8-11 replaces it.
   - The NES CI is not estimated, so rows carry `ci_unestimated`.
   - Ties are ordered with the seeded random order of M8-02 decision 7.
5. **Stacking:** the branch m8/gsea is m8/fgsea (#156) with m8/collections (#159) merged in. The lead merges main once those PRs land.
6. **The test author's choices (commit 69eb309) are accepted as the contract.**
   - **Ties:** a seeded order from `rng` is used when ties enter fgsea. `TuRanking.tu_ids` stays in natural TU-ID order.
   - **Run flags:** `layer_missing:abundance` and `too_few_ranked_tus` go in `table.attrs["run_flags"]`, the same channel as M7-03 decision 4. An empty table keeps every column.
   - **No abundance:** this means a missing or null `abundance_scale` with all-NA `abundance_log2`. Fewer than 2 ranked TUs gives `too_few_ranked_tus`.
   - **`tu_genes`:** the scoring genes, in `tu_position` order.
   - **Cells:** `leading_edge` and `flags` are sequences, every row carries `ci_unestimated`, and `flags` is the last column.
   - **Decision 3, made precise:**
     - the reported `q_value` is the main-run BH q;
     - the substituted BH, with the parent's re-test p in place of its main-run p, decides only `redundant`;
     - each re-test is one fgsea call with the same `n_perm`, made after all main runs. A re-test below `size[0]` still runs and flags `size_filtered`.
   - **Engine:** `fgsea_simple(prep, n_perm, sample_size, rng)` is called by name, so tests can patch it.
   - **`parent_group`:** the most general significant ancestor (by depth, then ID), or the term's own ID. It is NA for non-significant rows and non-GO collections.
