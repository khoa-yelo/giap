# M6-03: giapsim — operon nulls, planted co-abundance, planted discordance and the attacks

Milestone M6 · track F (simulator; tooling written by the implementer, tested by the test-author) ·
**pass 2** · depends on: nothing (merged giapsim M4-07/08/09, M5-04) · estimate ≈ 40 code lines in
`src/` (contract additions) + ≈ 260 lines in `tests/sim/` (+ ≈ 250 test lines) · statistical
(light): the statistician signs off the generators M6-09 gates on

## Spec sections implemented

- §11.3 "Planted features" and the abundance model: TU sizes, intra-operon gaps −20 to +40 bp,
  inter-TU gaps 50–300 bp, "strand runs with set convergent and divergent frequencies" (so CV and
  SW decoys exist in quantity), σ70 motifs upstream of TU starts and terminators after TU ends.
- §8.4 generator column for the suites `discordance` ("OP pairs from the t distribution") and
  `missed_operons` ("NO pairs independent; planted co-abundant pairs for power").
- §8.5 attack 1 (autocorrelated abundance) applied to Module 4.
- §7.4.7 "Calibration" items, which name exactly the three criteria M6-09 checks.

## What the simulator must add

`tests/sim/giapsim/scenarios.py`
```python
def operon_null(n_genes: int, *, rho: float = 0.0) -> Pair: ...
def missed_operon_power(n_genes: int, n_planted: int, effect_sd: float) -> Pair: ...
def discordance_power(n_genes: int, n_planted: int, effect_log2: float) -> Pair: ...
def operon_attack_ar1(n_genes: int, rho: float) -> Pair: ...
def operon_attack_detection(n_genes: int, short_genes_in: str) -> Pair: ...
```

1. **`operon_null`.** Within-TU gene abundances drawn so that Δ of OP pairs follows a Student-t
   (the `gene_sd` / `stoich_sd` machinery already produces a symmetric within-TU spread; this
   scenario fixes the parameters and records the truth), **no** internal promoters or terminators
   planted between operon members, NO pairs independent of each other. Truth: the generating t
   parameters and the list of OP, NO, CV and SW pairs.
2. **`missed_operon_power`.** `n_planted` NO pairs (adjacent, same strand, UniOP-like p < 0.5) are
   given nearly equal abundance (|Δ| ~ `effect_sd` × the background SD); the rest are independent.
   Truth: `truth/missed_operons.tsv` with the planted pair IDs.
3. **`discordance_power`.** `n_planted` OP pairs get a Δ of `effect_log2`, half with an internal
   promoter planted in the gap and half with an internal terminator, so the §7.4.5 flag rows 1 and
   2 have a truth set. Truth: pair ID, direction, which element was planted.
4. **`operon_attack_ar1`.** AR(1) abundance along the chromosome at ρ (0.5, 0.9) with **no** operon
   effect at all: neighbours look alike because of position, not transcription. This is the attack
   that decides whether the distance-reweighted decoys are adequate.
5. **`operon_attack_detection`.** Short, hydrophobic genes concentrated in one block (the ADR-0018
   lesson) with small intergenic gaps, so detection and distance are confounded: the decoy
   reweighting bins on distance, and the detection filter removes targets and decoys at different
   rates.
6. **UniOP-like probabilities.** The merged `giapsim.layers._operon_pairs` already writes
   `precomputed/operon_pairs.tsv`; these scenarios must guarantee a configurable share of pairs
   above and below 0.5 and record the *true* operon membership separately, so "UniOP is wrong"
   can be simulated (a planted co-abundant NO pair is exactly a UniOP false negative).
7. **Peptide-level truth** for M6-02: the merged evidence simulator (M4-09) already emits observable
   peptides and per-run observations; these scenarios set a detection regime in which a known set of
   genes is undetected while its partner is detected.
8. **Contract entry.** The new truth files are declared in `giapsim.contract_tables`, as M5-04
   declared C15.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| `n_planted` > available NO pairs | `ValueError` naming both counts |
| a scenario on a genome too small for 20 OP pairs | allowed, but the docstring says the suite needs ≥ 1,000 genes |
| ρ = 0 | `operon_attack_ar1` degenerates to `operon_null` and the test asserts that |
| no convergent pair generated (strand-run frequency 0) | `ValueError`: the decoys would not exist |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `operon_null(1800)` | ≥ 300 OP pairs, ≥ 200 NO pairs, ≥ 100 CV pairs, ≥ 30 SW pairs with both genes detected |
| 2 | `operon_null(1800)` Δ of OP pairs | `fit_student_t` recovers the generating parameters within tolerance; a KS test against the generating t passes at p > 0.001 |
| 3 | `missed_operon_power(1800, 50, 0.1)` | the 50 planted pairs have the smallest \|Δ\| in ≥ 80 % of seeds |
| 4 | `discordance_power(1800, 40, 2.0)` | 20 pairs with a planted promoter and Δ > 0, 20 with a planted terminator and Δ < 0, all in the truth file |
| 5 | `operon_attack_ar1(1800, 0.9)` | lag-1 Moran's I of TU abundance ≥ 0.5; **no** planted operon effect in the truth |
| 6 | `operon_attack_detection(1800, "prophage")` | detected fraction inside the block is ≤ half that outside, and the block's median intergenic gap is below the genome median |
| 7 | the same seed twice | identical outputs (the standing giapsim property) |

## Acceptance criteria

- [ ] `tests/unit/sim/test_operon_scenarios.py` (test-author): examples 1–7 and the degenerate
  rows; the truth files round-trip through `contract_tables`.
- [ ] A calibration-style sanity test (`tests/calibration/test_sim_operon_calibration.py`, R = 20
  fast): the generators do what they claim — the null has no excess of small |Δ| among NO pairs
  relative to decoys, and the power scenarios do.
- [ ] Runtime: a 1,800-gene scenario in < 5 s (the giapsim budget), recorded in `benchmarks/`.
- [ ] Documented in `docs/methods/operons.md` §"Simulations", including what each attack is built
  to break.

## Questions for the lead

1. Should `operon_null` plant *no* internal promoters/terminators at all (proposed), or plant them
   at the genome-wide false-positive rate the other layers use? **Proposed default:** none, so the
   level check is clean; the attack scenarios carry the realistic rates.
2. The backlog line "M6: the simulator plants stop-codon readthrough (C-terminal extensions,
   ADR-0013) for `cterm_extension` calibration" sits under the M6 heading but belongs to the
   novel-ORF / variant work. **Proposed default:** move it to the M7 backlog line; M6-03 does not
   implement it. Needs the lead's word so it is not silently dropped.

## Lead decisions (2026-10-22)

1. `operon_null` plants **no** internal promoters or terminators between operon members, so the
   level check is clean; the attack scenarios carry the realistic false-positive rates.
2. The stop-codon read-through simulator line moves to the M7 backlog (M7-12 / M7-16); M6-03 does not
   implement it. Done in this commit.
3. Feature loop: test-author → implementer (giapsim is tooling the implementer writes) → the
   statistician signs off the generators M6-09 gates on → code review.

## Lead decisions on the tests (2026-10-22, after 1e35b08)

4. **Pure classes in the null.** In `operon_null`, true non-operon same-strand pairs draw their
   UniOP-like p below 0.5 and true operon pairs at or above it, so the OP class M6-09 checks is
   exactly the true operons. "UniOP is wrong" is simulated by `missed_operon_power` and the attacks.
5. Example 3 is accepted as the test reads it. "Background SD" is the SD of Δ among unplanted
   true-NO pairs with both genes detected.
6. At least 60 % of planted pairs have both genes detected (examples 3 and 4): accepted.
7. `operon_attack_ar1(n, ρ)` returns the same configs as `operon_null(n, rho=ρ)` for every ρ; it is
   the named attack entry point (the tests assert ρ = 0; equality for all ρ is the implementation).
8. **Generating t.** The within-TU increment is Δ = loc + scale · t_df with loc −0.3, scale 0.5 and
   df 5, the parameters of the methods-note simulation (Slurm 1864589), recorded in
   `truth/operon_t.tsv`. The implementer adds this increment; the merged normal stoich/noise stays
   for the other scenarios.
9. An odd `n_planted` in `discordance_power` plants ceil(n/2) promoters and floor(n/2) terminators.

## Lead decisions on the statistician's block (2026-10-22, PR #222 comment 6356254106)

The discordance generators (`operon_null` for discordance, `discordance_power`) are approved as
built. The missed-operon generators change as follows, so M6-09's FDR gate can fail.

10. **Detection attack with teeth.** In `operon_attack_detection`, detection depends on abundance
    and is drawn **after** the within-TU walk; `operon_t` no longer requires `covariates_only`
    detection. The block's shift is mild enough that at least 25 % of the block's adjacent pairs
    have both genes detected (test at the default seed range).
11. **Short-gap decoys.** A new opt-in keyword, default off so existing scenarios are byte-identical,
    draws convergent TU-boundary gaps from a range reaching below 50 bp (default when on:
    uniform 10–150 bp). `operon_attack_detection` and `missed_operon_power` turn it on.
12. **Attack sanity test.** At R = 20, `operon_attack_detection` yields a larger fraction of small
    |Δ| (below the null's 25th percentile) among NO targets than among decoys in the same distance
    bin, one-sided z ≥ 2 pooled over replicates.
13. **Planted pairs on the attacks.** `operon_attack_ar1` and `operon_attack_detection` take
    `n_planted` and `effect_sd` (default 0, no `missed_operons.tsv` rows); when non-zero they plant
    as `missed_operon_power` does and write `truth/missed_operons.tsv`.
14. **Null sanity test.** Add a direct config assertion: the null has ρ = 0 and no planted pairs.
15. **Evidence commit.** The R = 1,000 numbers are re-run on the implementation commit and the log
    header names it.
16. **Not in M6-03.** Under the merged method the SW decoy floor keeps the minimum q ≥ 0.88 in the
    50–150 bp bin even with 10–15; the remedy is a method ADR for M6-09 (backlog, raised with Khoa).
    The statistician's NITs (choice 1 default-arm record, choice 4 realised planted count) go to
    M6-09's spec.

## Lead decisions on the round-2 tests (2026-10-22, after 2ea7253)

17. **Example 6 relaxed.** Decision 10 (≥ 25 % both-detected block pairs) cannot hold with block
    detection ≤ 0.5 × outside; example 6 now requires block detection < 0.9 × outside (lead edit of
    the existing test). The attack's teeth are checked by decision 12, not by the detection ratio.
18. **Golden digests** pinned to the locked numpy/pandas are accepted as the byte-identical guard;
    a dependency bump that changes float formatting re-pins them in the same PR.
19. **Decision 12 cut-off** is one global 25th percentile of |Δ| among the null's both-detected
    true NO pairs (as the tests read it), not per bin.

## Lead decisions on round 2 (2026-10-22, after d5ae874)

20. **CPU-independent golden digests.** The expression floats differ in the last ulp between
    nodes (barley-04 vs iron-01, Slurm 1869367), at fd3f3a3 as at d5ae874. The digests now hash
    floats at 9 significant digits; both nodes give the re-pinned values (Slurm 1869378).
21. **Decision 12 gates in the full tier only.** At R = 20 the expected CMH z is ≈ 1.1 (R = 1,000:
    z = 7.78, Slurm 1869345), so the fast tier prints z and checks that a sub-50-bp bin holds
    targets and decoys; z ≥ 2 is asserted at R = 1,000.
22. The implementer's deviations (block = a third of the chromosome, shift −6 on the logit, slope
    ≤ 2) are accepted for the statistician's re-check.
