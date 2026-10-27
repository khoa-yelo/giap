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

`src/gippa/function/gsea.py` (pure)
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
2. **Order and ties.** `TuRanking.tu_ids` is in natural TU-ID order (`TU_<n>`). Ties enter `prepare` in one seeded random order per run, drawn from `rng` (amended by lead decisions 4, 6 and 7, after M8-02 decision 7). M8-12 writes the post-shuffle order so R can reproduce it.
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
- [ ] `docs/methods/function-gsea.md`: TU collapse, centring and why (ADR-B), the family, the
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
7. **After the statistician's pre-code review (BLOCK on the test contract):**
   - **Rule 2 is amended.** Ties enter fgsea in a seeded random order (decisions 4 and 6), never in TU-ID order.
   - **Redundancy re-tests:**
     - they call `prepare` with min_size 1;
     - a re-test with 0 TUs gives `redundant = true` and the flag `size_filtered`;
     - a NaN re-test p gives `redundant` = NA.
   - **Attacks** (tests/calibration/attacks/test_attack_gsea_structure.py):
     - positional autocorrelation with scattered sets is gated;
     - positional autocorrelation with clustered sets (one run of 15 TUs) is a strict xfail in `calibration_full`;
     - detectability confounded with set membership is a strict xfail;
     - the default median-imputation ties arm is gated.
   - **Restricted claim**, which the methods note must state:
     - p-values are valid under TU exchangeability;
     - there is no level claim for genomically clustered sets under positional autocorrelation;
     - results are conditional on detection, with no `p_detect` adjustment;
     - NES has no CI (a §8.2 deviation);
     - "significant and non-redundant" is not an FDR-controlled list;
     - "no significant set" never means "no enrichment" (the floor is about 1/500).
   - **ADRs owed to Khoa:**
     - the autocorrelation null (rotation null, detrending, a `genomic_cluster` flag, or accept);
     - detectability selection (censored Blom scores, a `p_detect`-weighted null, or a flag);
     - the §8.2 NES CI;
     - ADR-B (centring).
8. **Khoa's decisions (2026-10-21, ADR-0025 to ADR-0028):**
   - **Clustered sets:** a `genomic_cluster` flag plus a restricted claim in pass 1. The circular-shift null comes in the hardening pass.
   - **Detection:** a restricted claim, "conditional on detection", plus a flag.
   - **NES:** accepted without a CI.
   - **Centring:** median centring is accepted (ADR-B).
9. **The `genomic_cluster` and `detection_confounded` flags (ADR-0025, ADR-0026). Accepted from the statistician's proposal.**
   - The flags are built as a stacked follow-up, **M8-05b** (branch m8/gsea-flags), so that M8-05 can be reviewed on its own.
   - The proposal's limits are accepted for pass 1:
     - the cluster flag is calibrated only for φ ≤ 0.5;
     - the detection flag is a screen. In one mixed world it is marginal at α = 0.01 (0.0147 against 0.0142).
   - The hardening ADRs (the rotation null, and censored or p_detect-weighted ranks) remove both limits.
   - The full proposal follows; evidence is in Slurm jobs 1866662 and 1866683.

### 9. Row flags `genomic_cluster` and `detection_confounded` (ADR-0025, ADR-0026)

Both flags are computed for every main-run row, which means every tested set (`size_tus` in
`function.gsea_size`). The redundancy re-tests compute no flags. Both are pure functions of the
inputs below and use no RNG. They are appended to `flags` after the existing values, in the order
`ci_unestimated`, `unestimated`, `size_filtered`, `genomic_cluster`, `detection_confounded`. Both
are new §9.4 values (scope "Function", or "Run and table rows" if the lead prefers) and go
through the schema.

#### 9a. `genomic_cluster`: the set's AR(1) variance inflation along the genome

**Meaning.** δ estimates how much a set's member TUs inflate the variance of its mean score,
relative to the TU-label permutation null, if TU scores were AR(1) with φ₀ = 0.5 along the ranked
TU order. Rows with δ ≥ τ carry no level claim under positional autocorrelation (ADR-0025).

**TU order.** Take the ranked TUs (`TuRanking.tu_ids`) and split them by replicon (`contig` of
their scoring genes; all the genes of a TU share one contig). Within a replicon r, keep the
natural TU-ID order and index the TUs i = 0 … n_r − 1. Let n = Σ_r n_r be the number of ranked TUs.
- *Do TU IDs follow genome order?* Yes. `gippa.annotate.operons.build_tus` numbers `TU_<n>` by
  contig (FASTA order), then start, then strand (+ before −). Natural sort of the IDs is therefore
  genome order within each contig. A TU that crosses the origin starts at its first gene, near the
  end of the ring, so it gets the last number on that replicon. It is still adjacent to TU 0
  through the circular distance, so the definition below is correct for it.
- Unranked TUs (undetected, or without a scoring gene) are skipped. Two members with only unranked
  TUs between them count as adjacent (d = 1). This is conservative, because it flags more.

**Distance.** For members a and b on the same replicon, d = |i_a − i_b|. On a circular replicon
(C1 `topology == "circular"`), d = min(d, n_r − d). Members on different replicons are never paired.

**Statistic.** Set S has k member TUs (its `size_tus`). Let w(d) = 0.5^d for 1 ≤ d ≤ 10, and 0
otherwise.
- C = (2/k) · Σ_{pairs a<b in S, same replicon} w(d_ab)
- P_w = Σ_r Σ_{unordered TU pairs i<j of replicon r} w(d_ij) (all ranked TUs, the same distance
  rule). In closed form:
  - circular: Σ_{d=1}^{min(10, ⌊n_r/2⌋)} c_d · 0.5^d, where c_d = n_r when 2d < n_r and n_r/2
    when 2d = n_r;
  - linear: Σ_{d=1}^{min(10, n_r−1)} (n_r − d) · 0.5^d.
- E[C] = 2 (k − 1) · P_w / (n (n − 1)), the mean of C when k members are placed at random among
  the n ranked TUs.
- **δ = (C − E[C]) / (1 − k/n)**, and δ = 0 when k < 2 or k ≥ n.
- **Flag: δ ≥ τ = 0.15.**

*Why this form.* Write V(φ) = (1 + φ)/(1 − φ). Under AR(1) with φ = φ₀, the variance of Σ_S(x − x̄)
is k(1 + C) − (k²/n)·V(φ₀) (with C taken to all lags). The permutation variance is about
k(1 − k/n), so their ratio is 1 + (C − E[C])/(1 − k/n) = 1 + δ. δ is centred at 0 for randomly
placed sets of any size, so large scattered sets (k up to 500) are not flagged by chance. It also
counts gaps: every-other-TU sets have δ ≈ 0.4, which a "runs of ≥ m adjacent TUs" rule misses
(below).

**Inputs.** `TuRanking` (tu_ids, tu_genes), the gene frame's `contig` column, the C1 replicons
(`contig`, `topology`), and the per-set TU membership from `tu_sets`. `gsea_table` gains one
argument, `replicons` (a `contig → topology` mapping or the C1 frame). Any topology other than
`"circular"` is linear, as in `gippa.context.autocorrelation`. Cost is O(Σ_S k·10) after sorting:
a member only needs its next ≤ 10 members on the same replicon. The circular seam needs care.

#### 9b. `detection_confounded`: the detection difference that can produce the observed direction

**Meaning.** Only detected genes are ranked. If a set's genes are detected less often than the
background for reasons other than abundance, its detected members over-represent high abundance,
which pushes ES up. A higher detection rate pushes ES down. The flag marks rows whose detection
rate differs from the universe *in the direction that would produce the observed ES sign*. For
those rows, the result may be a detection artefact (ADR-0026: "conditional on detection").

**Universe and counts (genes, not TUs).**
- Eligible genes are the genes the ranking could use if they were detected: non-partial, with
  finite `p_detect` when that column exists. This is `_scoring` without the detected/abundance
  conditions, and it is also the ORA universe of §7.7.3.
- M = number of eligible genes; N = number of scoring genes (detected, non-NA `abundance_log2`);
- n = number of eligible genes in the gene-level set S (from `GeneSets`, undetected members
  included); x = number of scoring genes in S.

**Statistic.** f = N/M, v = n f (1 − f) (M − n)/(M − 1), z = (x − n f)/√v, and z = 0 when v = 0
(the hypergeometric z of the detected count).

**Flag** (z₀ = 1):
- ES finite and ≠ 0 (`unestimated` rows included): flag when z · ES < 0 and |z| ≥ 1.
- ES = 0: flag when |z| ≥ 1. A NaN ES, which `fgsea_simple` does not return, is treated as 0
  (a defensive guard; statistician #174 Q3).

**Inputs.** The gene frame (`detected`, `abundance_log2`, `partial`, `p_detect`), the gene-level
sets, and the row's ES. No new configuration is needed. If z₀ should be configurable, it would
be `function.gsea_detection_z`, default 1.0, a protected statistical constant.

## Evidence (probe, R = 2,000 to 16,000 independent single-set replicates per arm, two seeds pooled)

The cells give the rejection rate at α 0.01 / 0.05 / 0.10. An asterisk marks a rate above the §8.4
bound α + 2√(α(1−α)/R) for that cell's row count. The bounds for the mixed arms are: R ≈ 2,200 gives
0.0142 / 0.0593 / 0.1128; R ≈ 4,300 gives 0.0130 / 0.0566 / 0.1092.

### `genomic_cluster` (δ ≥ 0.15)

The world has 1,500 TUs with E. coli-like operon lengths, within-TU ρ 0.8, and between-TU AR(1) φ.
"mix" arms draw k from 10 to 40 and a random shape for each replicate: scattered; j adjacent pairs;
r runs of 3 to 5; a spaced block plus scattered TUs; one run plus scattered TUs; or one run.

| arm | flag rate | flagged rows | unflagged rows |
|---|---|---|---|
| φ 0.5, one run of 15 | 1.000 | 0.080* / 0.194* / 0.284* | none |
| φ 0.5, five runs of 3 | 1.000 | 0.038* / 0.127* / 0.203* | none |
| φ 0.5, seven adjacent pairs | 1.000 | 0.025* / 0.090* / 0.155* | none |
| φ 0.5, three adjacent pairs | 1.000 | 0.020* / 0.067* / 0.122* | none |
| φ 0.5, every other TU (15) | 1.000 | 0.025* / 0.106* / 0.174* | none |
| φ 0.2, one run of 15 | 1.000 | 0.025* / 0.091* / 0.150* | none |
| φ 0.5, 15 scattered | 0.001 | 4 rows | 0.0110 / 0.0513 / 0.1034 (R 3,996) |
| φ 0.5, mix | 0.723 | 0.045* / 0.126* / 0.197* | **0.0104 / 0.0519 / 0.0961** (R 4,435) |
| φ 0.2, mix | 0.731 | 0.017* / 0.074* / 0.133* | **0.0088 / 0.0492 / 0.0980** (R 2,154) |
| φ 0 (exchangeable), mix | 0.731 | 0.011 / 0.054 / 0.103 | 0.0102 / 0.0451 / 0.1013 (R 2,152) |
| φ 0.7, mix (limit) | 0.734 | 0.092* / 0.192* / 0.268* | 0.0165* / 0.0574 / 0.1031 (R 2,125) |

- **Comparison rule** "≥ 50 % of members in runs of ≥ 3 adjacent TUs": in the φ 0.5 mix its
  unflagged rows reject at 0.0166* / 0.0700* / 0.1248*. It never flags seven pairs (0.090*) or
  every-other-TU sets (0.106*). Rejected.
- **The threshold.** By δ bin, in the φ 0.5 mix at α 0.05:
  - [0.10, 0.15): 0.059 (tune) and 0.059 (confirm);
  - [0.15, 0.20): 0.055 and 0.081;
  - [0.20, 0.30): 0.068 and 0.086.

  τ = 0.2 would leave three-pair sets (δ ≈ 0.18–0.20, 0.067*) unflagged, so τ = 0.15.
- **Flag rate under the null (random placement).** On the E. coli K-12 ring (2,816 RegulonDB
  14.5 operons as TUs), the flag rate is ≤ 0.001 for k ∈ {10, 15, 30, 60, 100, 200, 300, 500}.
  With 60 % of TUs ranked (n = 1,710) it is ≤ 0.002. The sd of δ is 0.02–0.03 at every k.
- **Realistic sets.** These are E. coli gene-name families and product-keyword sets with
  ≥ 10 TUs, standing in for COG/GO sets. With all TUs ranked, 10 of 39 are flagged:
  - phage (δ 1.43), flagellar (0.53), lys (0.43), fimbrial (0.42), met (0.25), rpl (0.24),
    rps (0.21), cys (0.18), arg (0.17), ins/IS elements (0.16).
  - None of the broad scattered sets is flagged: membrane 0.03, transcriptional 0.01, kinase
    −0.01, ABC transporter 0.08, tRNA 0.09, dehydrogenase 0.08, pur −0.01.

  With 60 % of TUs ranked, 2 of 17 are flagged (phage, ins).

### `detection_confounded` (z · ES < 0 and |z| ≥ 1)

The world has 2,000 TUs (E. coli-like, ρ 0.8). Detection follows logistic(β(a − s·[in set] + 0.5)),
where a is the abundance offset, and only detected genes are ranked. `p_detect` is a logistic fit
on a covariate −s·[in set] + σ·N(0, 1):
- σ = 0: the covariates explain the confounding;
- σ = 1: they explain part of it;
- σ = ∞: the covariates are uninformative, so the confounding is hidden.

"mix" arms draw s ~ U(0, 1.5), β ∈ {1, 2} and a set of 15 to 120 TUs for each replicate.

| world | flag rate | flagged rows | unflagged rows |
|---|---|---|---|
| s 1.5, β 2 (σ 0 / 1 / ∞) | 0.999 / 0.999 / 0.998 | 0.52* / 0.79* / 0.88* | ≤ 6 rows, 0 rejections |
| s 1.0, β 1, σ 0 | 0.818 | 0.041* / 0.168* / 0.274* | 0.0014 / 0.0178 / 0.0438 (R 730) |
| s 1.0, β 1, σ 1 | 0.813 | 0.051* / 0.176* / 0.267* | 0.0013 / 0.0134 / 0.0321 (R 748) |
| s 1.0, β 1, σ ∞ | 0.824 | 0.049* / 0.166* / 0.266* | 0.0028 / 0.0100 / 0.0213 (R 703) |
| null s 0, β 2, σ 0 | 0.160 | 0.000 / 0.028 / 0.064 | 0.0107 / 0.0545 / 0.1078 (R 3,359) |
| null s 0, β 2, σ ∞ | 0.160 | 0.003 / 0.022 / 0.050 | 0.0089 / 0.0533 / 0.0971 (R 3,359) |
| null s 0, β 1, k 15–120 | 0.133 | 0.006 / 0.030 / 0.069 | 0.0123 / 0.0544 / 0.1080 (R 5,184) |
| mix, σ 0 | 0.716 | 0.169* / 0.328* / 0.430* | 0.0147* / 0.0545 / 0.1077 (R 2,182) |
| mix, σ 1 | 0.715 | 0.170* / 0.329* / 0.430* | 0.0132 / 0.0587 / 0.1069 (R 2,199) |
| mix, σ ∞ | 0.723 | 0.168* / 0.332* / 0.436* | 0.0106 / 0.0544 / 0.1025 (R 4,263) |
| true effect +0.3, no confounding | 0.018 | (power arm) | |
| true effect −0.3, no confounding | 0.073 | (power arm) | |

**Rejected alternatives** (same records):
- **|z| ≥ 2 regardless of ES sign.**
  - The mix unflagged rows reject at 0.024* / 0.083* / 0.147*.
  - Rows with |z| < 1 still reject at about 0.08, because the detected-fraction z is about as noisy
    as the bias it predicts.
  - The rule flags 45–49 % of true ±0.3 effects.
- **|z_pd| ≥ 2 on mean `p_detect`.**
  - It works when the covariates carry the confounding: the σ 1 mix gives 0.055 unflagged.
  - It misses hidden confounding completely: in the σ ∞ mix the unflagged rows reject at 0.25*.
- **z₀ = 0.75 or 0.5.** The unconfounded null fails, for example 0.0131* / 0.0557 / 0.1106*. The
  sign condition selects against rows whose detection and ES move together by chance, and those
  rows are enriched for rejections. z₀ = 1 is the smallest value that keeps the unconfounded
  null at level.

## Test cases (inputs → flag) for the test author

### `genomic_cluster` (δ to 6 decimals; independent brute force and probe implementation agree)

Positions are indices in the ranked order of one replicon unless stated.

| # | ranked TUs per replicon (topology) | member positions | δ | flag |
|---|---|---|---|---|
| G1 | 100 (circular) | 0–14 (a run of 15) | 1.706342 | yes |
| G2 | 100 (circular) | 0, 11, 22, …, 88 (k 9) | −0.177427 | no |
| G3 | 100 (circular) | 97, 98, 99, 0, 1, 2 (across the seam) | 1.322169 | yes |
| G4 | 100 (**linear**) | 97, 98, 99, 0, 1, 2 | 0.781309 | yes; the 99–0 pair does not count |
| G5 | A 60 + B 40 (both circular) | A59, A20, B0, B20 | −0.063070 | no; A59 and B0 are never adjacent |
| G6 | 100 (circular) | 59, 60, 20, 80 (the G5 positions on one ring) | 0.197347 | yes |
| G7 | 100 (circular) | 0, 2, 4, …, 28 (every other TU, k 15) | 0.381740 | yes |
| G8 | 1,500 (circular) | five runs of 3 at 0, 300, 600, 900, 1200 | 0.822901 | yes |
| G9 | 1,500 (circular) | pairs at 0, 100, 200 + 9 singletons at 400, 500, …, 1200 | 0.183171 | yes |
| G10 | 1,500 (circular) | pairs at 0, 100 + 11 singletons at 300, 400, …, 1300 | 0.115831 | no |
| G11 | 1,500 (circular) | 0, 100, …, 1400 (k 15) | −0.018849 | no |
| G12 | 2 (circular) + 1 (circular) | both TUs of the 2-ring | 1.000000 | yes (c_1 = n_r/2 = 1, P_w 0.5) |
| G13 | any | k < 2, or k = n | 0 | no |

- **G14 (order and ranking).** The TU IDs are TU_1 … TU_12 and only the odd ones are ranked;
  ranked positions 0–5 are TU_1, TU_3, …, TU_11. Members {TU_1, TU_3, TU_5} form a run (d = 1).
  This gives the same δ as positions {0, 1, 2} in a 6-TU ranking.
- **G15 (natural sort).** TU_9 and TU_10 are adjacent; TU_1 and TU_10 are not.
- **Property tests:**
  - permuting gene rows, or relabelling the TU scores, leaves δ unchanged;
  - δ is invariant to rotating all positions on a circular replicon;
  - for a random placement, the mean of C over many draws equals E[C] (Monte Carlo).

### `detection_confounded` (gene counts; z to 4 decimals)

| # | M | N | n | x | ES | z | flag |
|---|---|---|---|---|---|---|---|
| D1 | 1000 | 500 | 40 | 12 | +0.4 | −2.5807 | yes |
| D2 | 1000 | 500 | 40 | 12 | −0.4 | −2.5807 | no (the bias works against the observed sign) |
| D3 | 1000 | 500 | 40 | 18 | +0.4 | −0.6452 | no (\|z\| < 1) |
| D4 | 1000 | 500 | 40 | 27 | −0.3 | +2.2581 | yes |
| D5 | 1000 | 500 | 40 | 27 | +0.3 | +2.2581 | no |
| D6 | 1000 | 600 | 30 | 15 | +0.5 | −1.1346 | yes |
| D7 | 1000 | 600 | 30 | 16 | +0.5 | −0.7564 | no (boundary pair with D6) |
| D8 | 1000 | 500 | 40 | 12 | NaN ES (defensive; not produced by `fgsea_simple`) | −2.5807 | yes |
| D9 | 1000 | 1000 | 40 | 40 | +0.4 | 0 (v = 0) | no |
| D10 | n = M (the set is the universe) | | | | any | 0 | no |

- **D11.** A partial gene in S counts in neither n nor M.
- **D12.** An undetected member gene counts in n but not in x. A detected member with NA
  `abundance_log2` also counts in n but not in x, since it is not a scoring gene.
- **D13.** Adding an undetected gene that is in no set changes M and therefore f, z and possibly
  the flag. This is the intended "background" behaviour.
- **D14.** Redundancy re-test rows do not exist as rows, so they carry no flags. The flag of a
  re-tested parent comes from its main run.

## Attack follow-ups (statistician, after implementation; tests/calibration/attacks/test_attack_gsea_structure.py)

- **Clustered sets.** Keep (d) "one run of 15" as a strict xfail on *all* rows, as ADR-0025
  requires. Add a gated arm that checks **unflagged rows hold level** in the φ 0.5 mix, plus a
  gated teeth check that a run of 15 and five runs of 3 are flagged in ≥ 99 % of replicates.
- **Detection confounding.** Keep (e) as a strict xfail on all rows. Add a gated arm that
  checks unflagged rows hold level in the σ ∞ mix, plus a gated teeth check that the flag
  covers ≥ 99 % of s 1.5 rows and ≥ 75 % of s 1.0 rows.
- **Calibration of the flags.** Add a gated arm for the unconfounded null, so that unflagged
  rows hold level and the flag rate is ≤ 0.25.

## Limits (for the methods note)

1. **`genomic_cluster` is calibrated for φ ≤ 0.5.** At φ 0.7 the unflagged rows are
   anti-conservative at α 0.01 (0.0165 vs 0.0143). The kernel fixes φ₀ = 0.5, because the step
   does not estimate φ.
   - A data-adaptive variant would use the lag-1 autocorrelation of ranked TU scores in place of
     φ₀. It is left to the hardening rotation null (ADR-0025 §2), which removes the need for the
     flag.
2. **The flag uses ranked-TU order, not bp distance.** Long gaps of unranked TUs count as
   adjacency. This is conservative. Abundance correlation that decays with distance in bp rather
   than in TU count is not modelled.
3. **`detection_confounded` is a screen, not a guarantee.**
   - In mixed-confounding worlds (s ~ U(0, 1.5), β ∈ {1, 2}) unflagged rows reject at about
     1.3α / 1.15α / 1.07α at α = 0.01 / 0.05 / 0.10 (pooled R ≈ 12,600: 0.0128 / 0.0567 / 0.1073),
     above the §8.4 bound; unscreened rows in the same worlds reject at about 0.25 at α = 0.05
     (Slurm 1867332). The committed arms gate unflagged rows at ≤ 1.5α.
   - Rows with |z| < 1 and the "wrong" sign keep a small bias that no detected-fraction threshold
     can screen out. The detected fraction is as noisy as the bias it predicts.
   - The fix is the censored-rank / p_detect-weighted null left to a later ADR.
4. **The detection flag reacts to *observed* detection,** which mixes abundance and detectability.
   It is sign-aware, so real effects whose detection difference has the same sign are rarely
   flagged (1.8 % and 7.3 % for ±0.3). A real effect with an opposite-sign detection difference
   is flagged, because the two cannot be told apart.
5. **Null flag rate and overdispersion.** Without confounding, 13–16 % of rows are flagged, about
   ½·P(|Z| ≥ 1) = 0.16; flagged null rows are conservative (0.02–0.03 at α 0.05), so unflagged
   null rows reject at about 1.1α (0.0551 / 0.1088 at α 0.05 / 0.10). The z uses the
   gene-level hypergeometric variance and ignores the within-TU correlation of detection. In
   genomes with long, tightly co-detected operons this raises the flag rate; it costs no level.
6. **The BH family includes flagged rows.** q values of unflagged rows come from a family that
   includes flagged rows; with 10 flagged sets among 40 the FWER of unflagged null rows rose from
   0.004–0.02 to 0.022–0.057 (Slurm 1867330; still ≤ α under the 1/500 floor; M8-11 re-checks).
   Both flag claims are per-row level, not FDR among unflagged rows.
7. **Unflagged ≠ calibrated beyond these worlds.** Both flags are tuned on the probe's AR(1) and
   logistic-detection worlds. Real sets can combine both structures; a row can carry both flags.
8. **Realistic-set check.** RegulonDB operons stand in for GIPPA TUs, and gene-name and keyword
   families stand in for COG/GO sets, because no COG/GO gene mapping for K-12 is on disk. Flag
   rates on real collections should be recorded by M8-15 on the first real E. coli run.
10. **After the code and statistician reviews (PR #171):**
    - **Code review B1, fixed by the lead.** `parent_group` took the term itself as a candidate, so a term with a shallower path through another parent grouped under itself even when it had a significant ancestor. It now uses the term only when no ancestor is significant. Depth is the shortest path.
    - **New tests** (tests/unit/function/test_gsea_review.py):
      - a term with two parents;
      - non-significant ancestors;
      - `redundant` decided by the substituted BH q, not the re-test p alone;
      - the `boolean` dtype.
    - **Statistician B1 and code review D1, fixed.** The methods note is docs/methods/function-gsea.md, and it states the restricted claim.
    - **Simple-stage floor (about 1/500):** this is a temporary departure from §7.7.4 and §8.1. M8-11 is a release gate for `enrichment_gsea.tsv`.
    - **Merge order:** M8-06 must not merge before M8-05b (the flags).
    - The remaining nits go to docs/backlog.md.

