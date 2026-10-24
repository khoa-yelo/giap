# Novel ORFs (M7)

Status: §"Class FDR" written for the statistician's round-2 review of PR #237 (M7-11). Code:
`src/giap/novel/fdr.py` (pure); the `novel_orfs` step wiring (M7-11b, PR #256) is
`src/giap/novel/class_call.py` and §"Wiring in the `novel_orfs` step". Spec: GIAP_SPEC.md §7.6.4,
§7.9.2, §8.1, §8.2; ADR-0009, ADR-0010; docs/specs/M7-11-novel-class-fdr.md, lead decisions 1–14;
docs/specs/M7-11b-class-fdr-wiring.md, lead decisions 1–7.

## Class FDR

Novel peptides have a lower prior than annotated ones, so the global search FDR understates their
error (§7.6.4). GIAP therefore estimates the error of a novel-ORF call **within its class**
(`intergenic`, `cds_outframe`, `antisense`, `nterm_extension`), with `stats.alpha.novel_orfs`
(0.01) as the call threshold.

### Regimes (rule 1, decision 1)
The regime is chosen per ORF from the source of its peptides, database first.

- `picked_decoy`: the run has decoy rows whose stripped protein IDs (longest
  `input.decoy_prefixes` match removed) are class-tagged `makedb` IDs
  `GIAPORF_{class}_{contig}_{start1}_{end1}_{p|m}`, and the ORF has at least one database
  (non-`denovo`) peptide.
- `denovo_subgroup`: otherwise, when the ORF has a de novo peptide.
- `unestimated`: everything else, for example decoys from a Bakta-protein search
  (`rev_<locus_tag>`). The pass-1 rule (≥ 2 distinct peptides) applies, and
  `class_fdr = unestimated`.

Database and de novo q-values are never pooled into one list. Their decoys answer different
questions (a search-engine competition and a genome placement), so a pooled list would have no
single decoy model.

### Picked class FDR (rule 2; ADR-0009)
**Matching.** Each target ORF is matched to the makedb six-frame ORF with the same
`(contig, strand, stop)`. A makedb ORF is stop to stop, so the codon position just past its 3′
end identifies it: the entry's 0-based end on the plus strand, and 3 nt below the entry's start
on the minus strand. For an ORF closed by a stop codon, that position is where `Orf.stop` (the
stop's lowest coordinate) lies. A 3′-open ORF (linear contig end or an `X` codon) has no closing
stop: `Orf.stop` is its last whole codon, and makedb ends the segment at that codon's end. Its key
is therefore `Orf.stop + 3` on the plus strand and `Orf.stop − 3` on the minus strand
(decision 13). When replicon lengths are given, both keys are wrapped modulo the length, so
origin-crossing entries match (decisions 8 and 14). Without lengths, an ORF whose entry crosses
the origin is unestimated. An `nterm_extension` ORF with no stop
match is paired by locus tag. Its makedb entry covers the upstream part plus the first 30 codons,
so its stop lies inside the gene.

The pairing rules have these consequences:
- When the ORF's class and the entry's tag differ, the entry tag wins, because pairing and the class
  list follow the entry. The ORF is flagged `class_mismatch` (decision 9).
- An ORF whose own class is untagged (`pseudogene`, `sorf`, `cterm_extension`) is unestimated
  (`untagged_class`, decisions 2 and 10).
- A tagged ORF with no entry is unestimated (`no_makedb_pair`, decision 7).

**Scores.** The score of an entry's target side is its matched ORF's best peptide. Only
non-de-novo peptides count, because makedb decoys mimic database-search targets only; adding
de novo target evidence with no decoy counterpart would be anti-conservative (decision 10). The
score of its decoy side is the best decoy row whose stripped ID is the entry. Within each entry
class, `picked_fdr` keeps the better side per entry, and ties go to the decoy. Entries hit only by
their decoy count as decoy wins.

**Unscored rows (decision 11, as M4-24).** C9 allows a missing `q_value` and a missing `score`.
Mapping a missing q to +∞ would rank an unscored decoy below every target and silently give D = 0.
In the statistician's all-null attack, 400 of 400 null ORFs then reached q ≤ 0.05. GIAP therefore
ranks a class only on a key that every target and decoy row of that class carries:
- (q ascending, then score descending) when every q is finite;
- otherwise score alone, when every score is finite;
- otherwise the class is unestimated with reason `decoys_unscored`.

**Paired share and rescaling.** f is the share of the class's distinct stripped decoy IDs that are
makedb entries. The decoy's class is the tag in its ID. Unpaired decoys never enter the list, and
each decoy win counts 1/f:

  FDR(t) = (D(t)/f + 1) / T(t), made monotone (ADR-0009, via `stats.tda.picked_fdr_rescaled`).

The unit is distinct IDs because decoy wins are counted per entry: the share of lost entries is
what the correction needs, while a row count would weight entries by their number of peptides. An
entry missing from `id_map` also loses its target (`no_makedb_pair`). That loss is symmetric, so
1/f is conservative there. 1/f matters when a decoy is lost but its target stays.

**80 % guard (decision 4).** A class with f < 0.80 falls back to unestimated
(`paired_share_low`), by analogy with the §7.3.2 tier guard; changing the guard needs an ADR. Below
the guard the cost is power (the ≥ 2-peptide rule), not validity. In the decoy-loss attack (10 % and
20 % of decoy rows re-pointed to unmapped tagged IDs), FDP stays at level.

A class with targets and no decoys gets q = 1/T, so it needs ≥ 100 target wins before any ORF
reaches 0.01 (the B.6 logic).

### De novo subgroup FDR (rule 3; ADR-0010)
Decoy hits come from the mapping work tables (`search = decoy`). Each peptide is taken at its best
decoy cost (minimum `n_diff`). Each hit takes the M3-14 `LocusClassifier` class of its mirrored
target interval, and is counted once per (peptide, class) at its minimum E (decision 10). A target
peptide enters the list of its ORF's class with its E. Per class, peptides are ordered by E:

  FDR_class(e) = (W_class(e) + c) / T_class(e), monotone q (`stratum_qvalues`),

where W sums the decoys' own π̂ (a missing π̂ counts as 1), and c is the largest π̂ over the
searched peptides, those with a target or decoy hit (decision 5), with a missing π̂ again
counting as 1. A class with no decoy hit gets q = c/T. A hit with no E (`e_target` or `e_decoy`
NaN) raises `ValueError` naming the peptide (decision 11). Dropping such a decoy would shrink W
and be anti-conservative: with 30 % of decoy E missing, FDP was 0.0133 / 0.069 / 0.104 against
bounds of 0.0115 / 0.053 / 0.103.

### Why the classes are never pooled
The false share differs strongly between classes. `intergenic` ORFs are mostly false, while
`cds_outframe` calls next to expressed genes are mostly true. A pooled list lends the clean
class's target wins to the noisy class. In the statistician's attack (intergenic 3 % true next to
cds_outframe 80 % true, FDP measured inside intergenic), the per-class estimate is at level. A
pooled list gives FDP **0.068 / 0.29 / 0.46** at q = 0.01 / 0.05 / 0.10. The de novo attack
(80 % of chance hits intergenic, correct peptides 5 % intergenic) is also at level per class.

### ORF q = minimum peptide q (decision 3): anti-conservative
For de novo ORFs, the ORF class q is the minimum of its peptides' subgroup q. That is the
best-evidence statistic, but it is **anti-conservative at ORF level**: a false ORF with k wrong
peptides gets k chances to draw a small q. In the review's generator, ORF FDP at
q = 0.01 / 0.05 / 0.10 was:

| Wrong peptides per false ORF | ORF FDP | q + 2 SE |
|---|---|---|
| 1–5 | 0.0185 / 0.078 / 0.113 | 0.0117 / 0.053 / 0.103 |
| 1–10 | 0.030 / 0.097 / 0.131 | — |

The 1–5 case failed in both confirmation blocks. That generator is a worst case: in real data,
chance hits of wrong de novo peptides rarely fall in one ORF. The cases that matter are systematic
de novo errors that cluster, such as reversed stretches (§8.5 attack 5) in intergenic space.
**M7-13 G4 must measure ORF-level FDP per class with false ORFs carrying ≥ 3 peptides**, including
attack 5. If it fails, an ADR chooses between a picked ORF-level competition on the genome decoy
and a per-ORF Bonferroni over peptides (decision 3). A strict xfail in
`tests/calibration/attacks/test_attack_novel_class_fdr.py` flags any change to decision 3.

### Wiring in the `novel_orfs` step (M7-11b)
Code: `src/giap/novel/class_call.py`, `giap.novel.evidence.class_rules` and `_call`,
`src/giap/steps/novel.py`. Spec: docs/specs/M7-11b-class-fdr-wiring.md, lead decisions 1–7.

**Inputs.** The step reads the C9 `peptides` (picked regime), and for the de novo regime it reads
`placement_costs`, `hits_exact` and `hits_tolerant`. All are optional. Without `peptides`, every
ORF is `unestimated` (`no_decoys`), which is the pass-1 output byte for byte. One regime is
assigned per ORF, with no cascade: a picked-regime ORF whose picked result is unestimated never
falls back to its de novo peptides (Q6). Each (regime, class) list is one §8.1 "Novel ORF" family,
and α is `stats.alpha.novel_orfs` everywhere.

**Call rule (rule 6).**
- Numeric `class_fdr` = q: `called = class rule ∧ q ≤ α`. q is compared as the unformatted float.
  The TSV shows `%.6g`, and q is never parsed back from the cell. There is no peptide minimum
  (§7.6.4), but the `cds_outframe` and `antisense` class rules still require 2 peptides.
- `class_fdr = unestimated`: `called = class rule ∧ n_peptides ≥ min_peptides_without_decoys`
  (§7.6.4 "Neither"), with the flag `class_fdr_unestimated`.
- `lost_pick` (Q1): the class is estimated, but the ORF's own decoy beat it or tied it. The cell is
  **1** and the ORF is never called.
  - Why it must not fall back to the ≥ 2-peptide rule: losing the pick is the competition's
    evidence that the target looks like chance. Falling back would call exactly the ORFs that the
    competition rejected.
  - In the statistician's attack (cds_outframe ORFs of 2–3 peptides, 25 % true, exchangeable
    false targets and decoys; `tests/calibration/attacks/test_attack_class_call_wiring.py`), the
    fallback reading gives FDP **0.98 / 0.31 / 0.35** of `called` at α = 0.01 / 0.05 / 0.10. The
    wired rule gives 0.0015 / 0.043 / 0.096, against bounds of 0.012 / 0.058 / 0.112 (R = 50).
  - Why 1 and not "the q at the decoy's rank": the target never enters the picked list, so it has
    no q. 1 is the conservative bound, and it keeps the cell numeric.
- The FDR claim covers the set {q ≤ α}. The class rule then keeps a subset of it. This preserves
  control as long as the rule does not prefer false ORFs. The rules (start codon, length, RBS
  support, ≥ 2 peptides) favour real ORFs. M7-13 G4 measures the called set itself,
  `class_rule_met ∧ q ≤ a`, from the side table.
- Fallback calls in a makedb run carry no error control. These are the `no_makedb_pair`,
  `no_target_score` and untagged-class ORFs. They are labelled `unestimated` (§8: never silently
  trusted), and G4's prediction set ("class q ≤ q") excludes them. M7-13 should therefore also
  report how many calls in each class came through the fallback.

**Class q floor (rule 8).**
- No picked q can fall below 1/T_wins, where T_wins is the class's entries won by their target.
  At any threshold, FDR = (D/f + 1)/T ≥ 1/T ≥ 1/T_wins, and the monotone q is a minimum of these
  values.
- No de novo q can fall below c/T, where T is the class's target peptides, because W ≥ 0.
- When the floor is **strictly** above α, nothing in the class can be called. The run then logs
  the WARNING `novel_orfs.class_q_floor_above_alpha`: **a class with fewer than about 1/α
  candidates cannot be called at α**. q and α are unchanged (Khoa, 2026-10-22).
  - T_wins = 100 at α = 0.01 gives a floor equal to α. That raises no warning, because the best
    ORF can still be called.
  - A class with no target win (picked) or no target (de novo), and an unestimated class, has no
    floor and raises no warning.
- The warning goes to the manifest `warnings` list (`--strict` exits 3), not to a run-flag token
  (Q4). It mirrors the gene-level `q_floor_above_alpha`.
- The picked floor counts entries won. The warning text quotes `n_targets`, the ORFs entering the
  list. When several ORFs share one entry, or many targets lose the pick, the two differ: read
  `q_floor` in the event data, not the target count.

**Paired share in v1 (rule 2, Q5).**
- No run input carries makedb's `db_id_map.tsv`. `makedb_id_map` therefore builds the map from
  every class-tagged ID in `protein_ids` (targets, and decoys with their prefix stripped). The
  `nterm_extension` locus tags come from `nterm_entries` on the run's own gene model.
- So every tagged decoy is an entry, the paired share is f = 1, and the 80 % guard cannot fire
  from the pipeline (backlog).
- This is **not anti-conservative**. The guard and the 1/f rescaling exist because an unpaired
  decoy is dropped from the list while targets stay, so D is undercounted. In v1 no tagged decoy is
  dropped: a decoy whose entry has no target ORF competes alone, and so counts as a decoy win.
  - In the re-pointed-decoys attack (20 % and 50 % of decoy rows moved to tagged IDs with no
    target, through `makedb_id_map`), FDP is 0.0081 / 0.044 / 0.089 at 20 % and
    0.0062 / 0.034 / 0.071 at 50 %, against bounds of about 0.012 / 0.053 / 0.104. FDP falls as
    the loss grows.
- Pairing can still be lost on the **target** side.
  - A `nterm_extension` ORF pairs through the locus tag of the run's own `nterm_entries`, and the
    other classes pair through the `(contig, strand, stop)` key. If the run's annotation or genome
    differs from the one makedb used, the ORF is `no_makedb_pair`.
  - Its decoy stays in the list, so the rest of the class becomes more conservative. The ORF
    itself, however, drops to the uncontrolled fallback above.
  - A real `db_id_map.tsv` would not change this, because f counts observed decoys.
- No share detects decoy rows that are missing from the input altogether, for example an engine
  that exports only some of its decoys. On the calibration generator, halving the tagged decoy
  rows gives FDP 0.017 / 0.096 / 0.20. Input decoys are assumed complete, as for the gene-level
  picked tier.

**De novo inputs (rule 3, Q11).**
- `denovo_inputs` keeps only `source = denovo` rows of `placement_costs`:
  - targets: ORF peptides with a finite `e_target`;
  - decoys: those peptides' decoy hit rows at their best decoy cost, carrying `e_decoy`;
  - weights: π̂ from `pi_denovo`;
  - c: the largest π̂ over the de novo peptides with a target or decoy hit.
- Database peptides' decoy hits and their π̂ (NaN, which would count as 1) never enter a de novo
  list. This is valid, because the list estimates the chance hits of wrong de novo peptides.
- The target E is the peptide's best target E over the genome, as M7-11 assumed. The step does not
  restrict it to the ORF's interval. `candidate_rows` keeps only peptides with a single novel locus
  and no hit in an excluding class, so the two usually coincide. They are not guaranteed to.
- Decoys of de novo peptides that are not in a de novo ORF still count (for example, banned,
  multi-locus or picked-regime ORF peptides). That is conservative.

**Calibration of the wired path.** `tests/calibration/test_class_call_calibration.py` runs
`class_fdr`, with its own `makedb_id_map` and `denovo_inputs`, plus `class_fdr_cell`.
- In the statistician's mutation probe (PR #256), these mutations fail the criterion:
  - de novo decoys restricted to ORF peptides (0.0175 / 0.084 / 0.116 against
    0.012 / 0.053 / 0.103);
  - decoy weights halved;
  - c = 0 (at α = 0.01);
  - half of the picked decoy rows dropped.
- Two mutations pass:
  - dropping the picked +1, whose effect is too small with T in the hundreds;
  - pooling the classes, because the generator gives every class the same prior. The M7-11
    `low_prior_class` attack covers that case.
- The suite measures numeric cells only. `test_attack_class_call_wiring.py` adds the `called`
  column.

**Diagnostics.** `novel/novel_class_fdr.parquet` has one row per `novel_orfs.tsv` row, in the same
order.
- Columns: `orf_id`, `class`, `fdr_class`, `regime`, `q` (float, NaN unless estimated and won),
  `class_fdr` (the TSV cell), `reason`, `class_mismatch`, `class_rule_met` and `called`.
- G4 calls at level a with `class_rule_met ∧ q ≤ a`. A `lost_pick` row has `q` NaN, so it is never
  in the set.
- One `novel_orfs.class_fdr` INFO event per (regime, class) list carries the counts, f, c, the
  floor and `n_called`. The manifest `novel_orfs.classes` block holds the same data.

Statistician (PR #256, M7-11b review, 2026-10-23): wiring, call rule and floor rule reviewed.
Verdict recorded on the PR.

### Exact systematic de novo errors (ADR-0050)

The class FDR treats wrong peptides as decoy hits treat them: a wrong read lands on the genome by chance, so a reversed or shuffled genome
decoy is as likely to receive it. A de novo read that is wrong but **exactly equals** the
translation of a non-coding six-frame segment breaks that assumption. It places with zero
cost, like a true peptide, and no decoy can match it except by chance, so no decoy-based
estimate counts it. When several such reads fall in one segment, the ORF also gets several
chances at a small q (decision 3). In M7-13 G4-k3 (20 injected false ORFs per class, 3–5 wrong
reads each) the ORF-level FDP was 0.075 for `cds_outframe` and 0.072 for `antisense` at both
q = 0.01 and 0.05, and 0.017 for `intergenic` at 0.01; the exact reads alone account for all of
it. GIAP therefore does not claim class-FDR control against exact systematic de novo errors.
q values and calls are unchanged. Calls that this limitation can affect carry
`denovo_only_unsupported`: every peptide is de novo, none is a database PSM, and there is no RBS
support beyond what the class rule already requires. RBS support at the 0.10 quantile is weak
evidence: in the simulator a non-coding start passes it with probability 0.060 for
`cds_outframe` and 0.075 for `antisense` (`intergenic` is never exempt; Slurm 1879351, fast tier,
R = 10), so the flag is a pointer to the evidence gap, not a second error estimate. Users who need such calls
should confirm them by a database search against the `giap makedb` six-frame database, or by
independent evidence.

### Limitations
- **c over hit peptides only.** ADR-0010 reads "the peptides … searched"; c here is the largest π̂
  over peptides with a target or decoy hit (decision 5). This is anti-conservative only when the
  highest-π̂ peptides have no hit at all: negligible at N in the thousands, visible in B.6-sized
  lists.
- **Class asymmetry.** A de novo target takes its ORF's class (the highest-priority member class),
  while its decoy analogue takes its own interval's class. A wrong peptide that joins an ORF of
  another class moves T toward that class, and no D moves with it. The effect is small; M7-13 G4
  reports FDP per class, not only per regime.
- **Target E.** A target peptide's E is its best E over all its target hits, not its hit inside the
  ORF. The M7-11b step does not restrict it, and documents this under "De novo inputs" above.
- `regime()` of a mixed set of ORFs returns the first of picked / de novo / unestimated held. The
  M7-11b step uses `orf_regimes` per ORF. The regime of each ORF is in the side table, and the
  per-class events hold the list counts. No separate regime-count event is logged.
