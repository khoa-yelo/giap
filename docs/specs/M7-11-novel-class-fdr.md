# M7-11: Novel-ORF error control — picked class FDR with decoy rows, subgroup FDR for de novo

Milestone M7 · track S (statistics) · **pass 2** · depends on: M7-06, M7-09 (the wired
`novel_orfs` step), merged `giap.stats.tda.picked_fdr` and `monotone_q`, M3-10 placement FDR
(`stratum_qvalues`, ADR-0010), M3-23 makedb IDs · estimate ≈ 300 code lines (+ ≈ 500 test lines)
· **stats** (statistician reviews the methods note and simulations before the code)

## Spec sections implemented

- §7.6.4: "Novel peptides have a lower prior than annotated ones, so the global search FDR
  understates their error. **With decoy rows in the input** (for example after a search against the
  `makedb` database): a picked target–decoy FDR at ORF level, within each class. Decoys pair
  through the class-tagged `makedb` IDs (§7.9.2). A call needs class q ≤ 0.01. **De novo input:** a
  subgroup placement FDR. Each decoy hit is assigned the locus class of its mirrored target position
  (reversed or shuffled windows keep the frame), then FDR_class = π̂_wrong (D_class + 1) / T_class.
  Calls need q ≤ 0.01. **Neither:** a call needs ≥ 2 distinct peptides and carries
  `class_fdr = unestimated`."
- §8.1 row "Novel ORF | ORF | decoy rows in the input (picked, per class); genome decoy for de novo
  | subgroup FDR | 0.01"; §7.9.2 (`>GIAPORF_{class}_{contig}_{start1}_{end1}_{p|m}`; "No decoys.
  Search engines generate their own; with DIA-NN, use `--report-decoys`"); §8.2 ("decoys mimic
  targets"); `stats.alpha.novel_orfs` (0.01); ADR-0009 (picked FDR rescaled by the paired share);
  ADR-0010 (per-peptide π̂ weights).

## Public interface

`src/giap/novel/fdr.py` (pure)
```python
Regime = Literal["picked_decoy", "denovo_subgroup", "unestimated"]

def regime(peptides: pd.DataFrame, orfs: Sequence[Orf]) -> Regime: ...
def picked_class_q(orfs: Sequence[Orf], peptides: pd.DataFrame, id_map: pd.DataFrame | None
                   ) -> dict[str, float]: ...              # orf_id -> class q
def denovo_class_q(orfs: Sequence[Orf], target_hits: pd.DataFrame, decoy_hits: pd.DataFrame,
                   classifier: LocusClassifier, pi: Mapping[str, float]) -> dict[str, float]: ...
```

## Rules

1. **Regime per run.** `picked_decoy` when C9 holds decoy rows whose protein IDs (prefix stripped,
   `input.decoy_prefixes`) are makedb class-tagged IDs; `denovo_subgroup` when the candidate
   peptides are `source = denovo`; otherwise `unestimated` (the pass-1 behaviour). A run with both
   database and de novo candidates applies the regime per ORF from its peptides' source
   (decision 1).
2. **Picked, per class.** Target ORF ↔ makedb entry: match the reconstructed ORF to the makedb
   six-frame ORF with the same `(contig, strand, stop)` (the makedb ORF is stop-to-stop, so the
   stop identifies it; `nterm_extension` entries by locus tag). Score of a target ORF = the best
   (minimum) input q, then maximum score, over its peptides; score of its decoy = the same over
   decoy rows whose stripped ID is `<decoy_prefix>GIAPORF_…` of that entry. Within each class,
   `picked_fdr` (ties go to the decoy), with ADR-0009's rescaling by the paired share of the
   class's decoys. Classes without makedb tags (`pseudogene`, `sorf` by locus tag only,
   `cterm_extension`) fall back to `unestimated` (decision 2).
3. **De novo subgroup.** Decoy hits come from the mapping work tables (`hits_exact`,
   `hits_tolerant`, `search = decoy`) at each peptide's best decoy cost; each is assigned the locus
   class of its mirrored target interval with the M3-14 `LocusClassifier` (decoy rows already carry
   the mirrored target interval, M3-04). Within each novel class, a peptide-level list ordered by E
   (the M3-10 ranking score) with the ADR-0010 weighted form `FDR_class(e) = (W_class(e) + c) /
   T_class(e)`, monotone q. ORF class q = the **minimum** peptide class q among its peptides
   (decision 3).
4. **Call.** `called` = class rule (M7-09) **and** class q ≤ `stats.alpha.novel_orfs`; under
   `unestimated`, the ≥ 2-peptide rule (M7-09 rule 8). `class_fdr` = the class q (float) or
   `unestimated`; the flag `class_fdr_unestimated` only in the last case.
5. **Reporting.** `detection_summary`-style rows are not added; the per-class T, D, paired share
   and regime go to `logs/novel_orfs.log` and `events.jsonl` (`novel_orfs.class_fdr`), and to the
   manifest.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| decoy rows present but none class-tagged (search not against `makedb`) | `unestimated` + warning `decoys_not_makedb` |
| a class with targets and zero decoys | q = (0 + 1)/T per the +1 rule (picked: as `picked_fdr`); small classes cannot pass 0.01 (B.6 logic) |
| paired share of a class's decoys < 80 % | that class falls back to `unestimated` (ADR-0009's guard, by analogy with MAYU fallback) |
| de novo with no decoy hit in a class | q = c / T (ADR-0010) |
| ORF whose peptides mix database and de novo | the database regime when it is `picked_decoy`, else de novo (decision 1) |
| π̂ unavailable for a de novo peptide | 1 (ADR-0010) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | intergenic class: 150 target ORFs outscore their decoys before the first decoy win | class q of those 150 = 1/150 = **0.0067** → callable at 0.01 (B.6 analogue) |
| 2 | intergenic class: 50 targets, no decoy wins | q = 1/50 = **0.02** → none called |
| 3 | `cds_outframe` ORF matched to `GIAPORF_cds_outframe_c1_101_400_p`, its decoy `rev_GIAPORF_cds_outframe_c1_101_400_p` scores higher | the target loses the pick (counts as a decoy win) |
| 4 | de novo: a decoy hit whose mirrored target interval is intergenic | counts in D_intergenic only |
| 5 | de novo: ORF with peptide class q 0.004 and 0.03 | ORF class q **0.004** |
| 6 | search against Bakta proteins with `--report-decoys` (`rev_<locus_tag>` decoys only) | `unestimated` + `decoys_not_makedb` |

## Acceptance criteria

- [ ] `tests/unit/novel/test_fdr.py`: examples 1–6 and every degenerate row.
- [ ] **Calibration suite `novel_orfs`** (§8.4 "Novel-ORF class FDR | planted ORFs + false peptides
  + de novo errors | FDR | 50 | 300"): in M7-13, per regime, mean FDP ≤ q + 2 SE at q ∈ {0.01, 0.05,
  0.10}, ADR-0015 two-block confirmation.
- [ ] `docs/methods/novel-orfs.md` §"Class FDR" with the three regimes, the ORF ↔ makedb matching,
  decision 3 and its risk; statistician sign-off.

## Calibration and attacks owed

- The `novel_orfs` suite (M7-13), both regimes.
- Attack (new, M7-13): de novo errors that reverse a stretch (§8.5 attack 5) concentrated in
  intergenic space — the class FDR must stay at level with shuffled decoys.
- Attack (new): planted false peptides drawn preferentially from GC-poor intergenic regions,
  where the shuffled decoy composition differs most from the targets.

## Questions for the lead

1. Mixed-source runs (database and de novo candidates). **Proposed default:** per-ORF regime by its
   peptides' source, database first; the methods note states that the two classes of q-values are
   not pooled into one list.
2. Classes without makedb tags (`pseudogene`, `sorf`, `cterm_extension`). **Proposed default:**
   `unestimated` with the ≥ 2-peptide rule. Pairing sORFs by `rev_<locus_tag>` is possible but they
   are annotated genes, already covered by the gene-level FDR of M4-12.
3. ORF-level q from peptide-level subgroup q (de novo): the minimum (proposed) is the natural
   "best evidence" statistic, but it is anti-conservative at ORF level when ORFs have many
   peptides. **Proposed default:** minimum, *and* the statistician measures ORF-level FDP in M7-13;
   if it fails, the ADR options are a picked ORF-level competition on the genome decoy or a
   Bonferroni-on-peptides-per-ORF correction.
4. The 80 % paired-share guard per class is borrowed from §7.3.2 by analogy. **Proposed default:**
   adopt it, recorded in the methods note; change only by ADR.

## Lead decisions (2026-10-22)

1. **Mixed-source runs:** the proposed default. Each ORF's regime follows its peptides' source,
   database first; database and de novo q-values are never pooled into one list (methods note).
2. **Untagged classes** (`pseudogene`, `sorf`, `cterm_extension`): `unestimated` with the
   ≥ 2-peptide rule; no `rev_<locus_tag>` pairing in v1.
3. **ORF-level q (de novo):** the minimum of its peptides' subgroup q, and the statistician measures
   ORF-level FDP in M7-13. If that fails, an ADR chooses between a picked ORF-level competition on
   the genome decoy and a per-ORF Bonferroni over peptides. The methods note names the
   anti-conservative direction now.
4. **80 % paired-share guard per class:** adopted from §7.3.2 by analogy and recorded in the
   methods note; changes only by ADR.
5. **Test-author resolutions (11a6ca30) accepted:** frame layouts, NaN for unestimated or pick-losing
   ORFs, paired = stripped ID in `id_map`, `id_map=None` parses the tagged IDs, c = largest searched
   π̂ (missing = 1), q before score with ties to the decoy.
6. **Split (MVP first):** this issue ships the pure `giap.novel.fdr` module that the tests pin.
   Wiring it into the `novel_orfs` step (the `decoys_not_makedb` warning, the rule 4 call, the
   `class_fdr` and flag columns, rule 5 logging and manifest) is **M7-11b**, with step-level tests
   written once the implementer has fixed the step inputs.
7. **A tagged-class ORF with no makedb entry** is `unestimated`, with a `no_makedb_pair` reason.
8. **Origin-crossing makedb entries** are matched modulo the replicon length when lengths are
   given. Without lengths the ORF is `unestimated` (backlog item); `picked_class_q` keeps its
   signature.
9. **Class disagreement:** when the reconstructed ORF class differs from the entry's tag, the entry
   tag wins, because pairing and the class list follow the makedb entry. The ORF is flagged
   `class_mismatch` (M7-11b column).
10. **Implementer choices (49811ae3), accepted:**
    - Decision 2 beats decision 9: an untagged own class stays `untagged_class`.
    - Locus-tag pairing is used only for `nterm_extension` without a (contig, strand, stop) match.
    - The paired share counts distinct stripped decoy IDs.
    - The picked score uses only non-de-novo target peptides.
    - A de novo target's class is its ORF's class.
    - Decoys are counted once per (peptide, class) at the minimum E.
    - `regime()` of mixed ORFs returns the first of picked / de novo / unestimated.
    - `stats.tda.picked_fdr_rescaled` now holds the ADR-0009 rule, and M4-12 `_picked_q` uses it
      (same arithmetic). A paired share of 0 now raises instead of giving NaN.
    - M7-11b step tests pin decisions 7–9.
11. **Statistician block on #237 (comment 6359036495).**
    - B1, picked: a class is ranked on (q, score) only when every target and every decoy row has a
      finite q. Otherwise it is ranked on score alone, if all rows have a score. Otherwise the class
      is `unestimated` with reason `decoys_unscored`, mirroring the M4-24 shared-key rule.
    - B1, de novo: a NaN E raises `ValueError`.
    - B2: `docs/methods/novel-orfs.md` §"Class FDR" records the decision-3 anti-conservatism with
      the statistician's numbers. M7-13 G4 measures it with false ORFs of ≥ 3 peptides, per class.
    - Implementer round 2. The lead removes the two Block-1 strict xfails after the fix; the
      decision-3 xfail stays until M7-13.
12. **Statistician round 2 (4704a61e), B3:** in the all-q-finite branch, score breaks ties only
    when every decoy row in the class has a score. Otherwise the class ranks on q alone, and ties
    go to the decoy (decision 5). This matches M4-24's shared-key rule. A target without a score
    may still rank last, which is conservative.
    - The strict xfail `test_picked_unscored_decoys_do_not_lose_q_ties` is removed by the lead
      after the fix.
    - Not blocking: when only some targets lack q, M7-11 ranks on score alone while M4-24 keeps q.
      This is safe and goes to the backlog for alignment.
