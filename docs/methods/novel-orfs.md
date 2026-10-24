# Novel ORFs (M7)

Status: §"Class FDR" written for the statistician's round-2 review of PR #237 (M7-11). Code:
`src/giap/novel/fdr.py` (pure; the `novel_orfs` step wiring is M7-11b). Spec: GIAP_SPEC.md §7.6.4,
§7.9.2, §8.1, §8.2; ADR-0009, ADR-0010; docs/specs/M7-11-novel-class-fdr.md, lead decisions 1–14.

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
  ORF. The two coincide only if `target_hits` is restricted to the ORF's interval, which the step
  (M7-11b, from M7-06) must guarantee or document.
- `regime()` of a mixed set of ORFs returns the first of picked / de novo / unestimated held. M7-11b
  logs the per-ORF regime counts and the share of decoys that were scored.
