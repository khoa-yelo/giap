# ADR-0056: P(miss) is reported without a q and stays out of the discordance BH family
- Status: proposed
- Date: 2026-10-24
- Spec sections: §7.4.3 "One detected, one not" and "Multiple testing", §8.1 row "Operon
  discordance", §8.5, §9.3 `pairs.tsv`
- Needs Khoa: yes. It departs from §7.4.3's "BH across all tested OP pairs (both kinds together)"
  and from the §8.1 row, which names the detection-limit model as part of the discordance family.

## Context
§7.4.3 tests two kinds of OP pair. A pair with both genes detected gets a two-sided Student-t
p-value against the within-operon Δ distribution. A pair with exactly one gene detected gets
p = P(miss): the probability, under M6-02's detection-limit model, that the undetected partner u
would be missed, given the detected partner's abundance a_v and the empirical within-operon ratios.
The spec then puts both kinds into one BH family at q ≤ 0.1.

P(miss) is a model probability of an event that has already happened. It is computed only on the
pairs where u *was* missed, so conditional on entering the family it is not super-uniform under any
null. M6-index lead decision 1 and M6-02 lead decision 4 therefore reported P(miss) with no q and
kept it out of the BH family, pending M6-09 attack C and an ADR. M6-05H (PR #311) implements that:
a `one_detected_pmiss` row has `p_value` = P(miss), `method = detection_limit`,
`null_model = p_detect_model` and `q_value` NA, and the BH family holds only the `both_detected_t`
rows.

M6-09 recorded attack C with the *oracle* P(miss) (the generating model). This ADR adds attack C
through the step, with the M6-02 model *fitted*
(`tests/calibration/attacks/test_attack_operons_pmiss.py`).

## Decision
1. P(miss) stays in `pairs.tsv` `p_value` for `test = one_detected_pmiss` rows, as §7.4.3 and §9.3
   name it, with `method = detection_limit` and `null_model = p_detect_model`. It gets **no q** and
   **never enters the BH family**; `operon_flags.tsv` makes no discordance call from it.
2. The §8.1 discordance family is every `both_detected_t` row of the run, one BH family at
   `stats.alpha.operon_flags`.
3. The report and the methods note describe P(miss) as "the model's probability that this partner
   would be missed if it were co-expressed", never as a significance level, and say that it is
   anti-conservative when detectability varies between proteins in ways the model does not see.
4. A scored `one_detected_pmiss` row carries `ci_unestimated`, as the `both_detected_t` rows do.
   Its `effect` and CI stay NA (M6-05 item 4), and §8's "effect sizes with 95 % CIs accompany every
   p or q" is then met by the flag rather than silently.
5. Pooling is revisited only with a method that makes the one-detected arm a valid p-value under
   protein-level detectability heterogeneity. For example, a model with a per-protein random
   effect, or a replicate-model P(miss) where runs identify that effect. Any such method must pass
   the `gene_effect` arms of attack C with the pooled family's any-call share within the §8.4
   bound.

## Evidence (attack C through the step; fitted M6-02 model)
Setup: `operon_null(1800)` genomes (TU blocks and the `OPERON_T` within-TU walk). Detection is
redrawn at peptide level, with logit = −11.5 + 0.5 a − 0.03 L + 0.8 GRAVY + o_g, and a gene is
detected at ≥ 1 seen peptide (k_min = 1). C11 holds the seen peptides, and P(miss) comes from
`giap.steps.operons.detection_limit` → `discordance`. Arms:
- `null`: o_g = 0;
- `gene_effect`: o_g ~ N(0, 1), a per-protein detectability the model does not see;
- `gene_effect_05`: o_g ~ N(0, 0.5²);
- `block`: o_g = −1.5 on the ADR-0018 prophage block of `operon_attack_detection`.

Full tier: Slurm 1884524, R = 300 genomes per arm, 2,331 s on rye-01. It ran at commit 65d9fedd
on branch stat/m6-05h-attack-c: #311 round 1 (fb1410cf, u scored with **u's** TU model) merged
into #310 (bbd6d4a4). The fast tier (R = 20) was rerun on #311's review round (4551ab61: v's TU,
the k_min event) and reproduces the round-1 fast numbers to every printed digit.
- On these simulations OP pairs never cross a TU boundary (pure UniOP classes), so the v-TU fix
  cannot move them.
- At k_min = 1, "fewer than k_min observable peptides seen" is the simulated non-detection, so the
  event fix cannot move them either.
- On real data, OP pairs that cross a C5 TU boundary can move.

Rates are at α = 0.01 / 0.05 / 0.10. The step's P(miss) matches an independent recomputation
to ≤ 6·10⁻¹⁶ on every row.

| arm | mean P(miss) vs miss rate | recalibration (intercept, slope) | share of one-detected rows with P(miss) ≤ α | rate over at-risk pairs (u missed and P(miss) ≤ α) | pooled per-pair rate (t + P(miss)) | any BH call at q ≤ 0.1: t family / pooled |
|---|---|---|---|---|---|---|
| null | 0.197 vs 0.196 | (−0.03, 0.97) | 0.0037 / 0.028 / 0.068 | 0.0007 / 0.0055 / 0.0134 | 0.0082 / 0.0430 / 0.0895 | 0.077 / 0.063 |
| gene_effect_05 | 0.172 vs 0.209 | (0.18, 0.83) | 0.0114 / 0.061 / 0.125 | 0.0024 / 0.0128 / 0.0261 | 0.0104 / **0.0544** / **0.1091** | 0.050 / 0.093 |
| gene_effect | **0.119 vs 0.236** | **(0.40, 0.56)** | **0.067 / 0.198 / 0.317** | 0.0159 / 0.0467 / 0.0747 | **0.0318 / 0.107 / 0.183** | 0.067 / **0.717** |
| block | 0.216 vs 0.252 | (0.24, 0.92) | 0.0061 / 0.037 / 0.084 | 0.0015 / 0.0093 / 0.0211 | 0.0090 / 0.0431 / 0.0912 | 0.110 / **0.140** |

- Pooled sizes: 65,738 to 79,281 one-detected rows and 260,456 to 347,647 at-risk orientations
  per arm.
- In every arm the oracle (the true model with its offset, Δ from the generating t) is calibrated:
  recalibration slope 0.96 to 0.99, intercept within ±0.03.
- The step's ratio sample (both-detected Δ only) costs nothing: with the true model it gives a
  mean P(miss) equal to the oracle's to 3 decimals. The departures are all in the fitted model.
- The any-call bound for R = 300 is 0.1 + 2 SE + 1/R ≈ 0.138. Flags per genome at P(miss) ≤ 0.05
  are 6.4 (null), 14.8 (SD 0.5), 52.3 (SD 1) and 8.1 (block).

Reading:
- **Correct model: P(miss) is a calibrated probability.** As a discrete p-value of the event "u
  missed", it is conservative over the at-risk pairs. The share of *flagged* one-detected pairs is
  below α, and pooling would not break the family on this arm.
- **Protein-level detectability breaks it.** With o_g ~ N(0, 1) the fitted P(miss) is roughly
  half the true miss rate. The conditional flag share is several times α, the per-pair rate of
  the pooled mixture exceeds α, and the pooled BH family makes a false call in most genomes; the
  shipped t family does not.
  - The mechanism is the conditional (Huggins) likelihood. It assumes that detectability is a
    function of the three features alone. Genes that were detected have high o_g, so the model
    over-predicts detection for the genes that were not.
  - Real proteomes have such heterogeneity: ionisation, digestion, membrane location and peptide
    sharing all vary beyond length and GRAVY.
- Spatially clustered detectability (`block`) biases P(miss) downward more mildly (intercept
  +0.24). It already takes the pooled family's any-call share to 0.140, just over the bound.
- With SD 0.5 the pooled family stays inside its any-call bound, but the per-pair mixture rate is
  already above α at 0.05 and 0.10.

The oracle arm of M6-09 (generating model, gene-level detection) stays as recorded there:
0 / 0.0127 / 0.0832 on `operon_null` and 0.0176 / 0.0649 / 0.1189 on the detection attack.

## Options considered
- **(a) Pool, as §7.4.3 says.** This is safe only when M6-02 is correctly specified. Under a
  per-protein detectability SD of 1 logit, the pooled family's FWER/FDR under the complete null is
  far above 0.1. Rejected.
- **(b) A separate BH family over the P(miss) rows.** Conditional on entering the family, P(miss)
  is not super-uniform, so BH over the one-detected rows alone is not valid. The valid discrete
  family would have to count every at-risk pair, including those where u was detected (p = 1). It
  would still inherit (a)'s misspecification failure. Rejected.
- **(c) No q, outside the family, labelled a probability (this ADR).** No claim is made that can
  fail, and the number is still reported for ranking and inspection.

## Consequences
- The §8.1 discordance family is narrower than the spec's. The M6 discordance claim is about
  both-detected pairs only (and is itself pending ADR-0055).
- One-detected pairs are never flagged as discordant. A user who thresholds `p_value` on
  `one_detected_pmiss` rows gets the false-flag shares in the table above. The report must say so.
- `docs/methods/operons.md` "One detected (M6-05H)" and "Family" cite this ADR in place of
  "needs an ADR (pending M6-09 attack C)". The two *Superseded* notes in docs/specs/M6-02 and M6-05
  point here.
- When the #311 review round lands (u scored with v's TU model; the one-detected event = fewer than
  k_min observable peptides seen), attack C through the step is rerun and the table updated. On
  real data, pairs that cross a TU boundary can move.
