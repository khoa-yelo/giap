# Placement FDR: genome decoys, chance hits and per-stratum q-values

Implements the genome side of GIAP_SPEC.md §7.2.1 (decoy translations), §7.2.5 (placement FDR, G3),
§8.1 rows "Peptide placed at a locus", §8.2 "decoys mimic targets", §8.5 attack 5 and Appendix
B.2/B.6. Issues: M3-03 (decoys, this section), M3-09 (chance hits), M3-10 (strata and q-values),
M3-21 (genome-level calibration). Engine for the q-values: `docs/methods/target-decoy.md`.

Sections: Genome decoys (M3-03, **signed off below**); Chance hits and the ranking score (M3-09,
**signed off below**); π̂_wrong in the stratum FDR (M3-10, **rule signed off below**); Strata and
ranking, π̂_wrong sources, Same-strain input, Calibration (M3-10, **signed below; implementation
approved at 031378c after two same-strain fixes**, see "Implementation review (M3-10)"). Genome-level
calibration with tolerant hits and de novo errors, and attack 5, are M3-21.

## Genome decoys (M3-03)

### What the decoy has to estimate

A placement is wrong when a wrong peptide (a misidentified spectrum or a de novo error) happens to
occur in the six-frame translation. For one wrong peptide the number of such chance hits depends
only on the search space: the set of *distinct* stop- and X-free L-mers of the six frames
(each peptide counts once at its best cost, §7.2.5 step 1, so multiplicity of a word does not
matter), and the peptide law. The decoy search space must give a wrong peptide the same hit
probability as the target space, while giving a correct peptide (which is in the target by
construction) only that same chance probability. Then D (peptides with a decoy hit at cost m)
estimates the number of wrong peptides with a chance target hit at cost m, and π̂_wrong (D + 1)/T
is the stratum FDR (M3-10).

Formally, for a peptide law Q and a space with distinct words W, P(hit) = Σ_{w ∈ W} Q(w) (exact
search; for ≤ m differences replace W by its m-neighbourhood). The decoy is valid for Q when
Σ_{W_D} Q = Σ_{W_T} Q.

### The two constructions

Both act on the I→L-folded uint8 frames of `build_sixframe`, frame by frame, and return frames of
the same contig, strand, offset and length (`build_decoys`).

| | reverse (database input) | shuffle (de novo input) |
|---|---|---|
| Operation | frame string reversed | inside each stop-free segment, non-X residues cut into consecutive windows of 50 from the segment start (last window shorter); each window permuted uniformly (Fisher–Yates), seeded `for_step("sixframe_index", seed)` |
| `*` and `X` | mirrored (position n−1−i) | fixed in place |
| Frame composition | identical | identical |
| Per-window composition | identical (mirrored windows) | identical (exactly, per 50-aa window) |
| Stop-free segment lengths and X positions | same multiset, mirrored | identical |
| Number of L-mer positions without `*`/`X` | identical | identical |
| Local order | kept, reversed (dipeptide xy → yx) | destroyed below 50 aa; composition gradient above 50 aa kept |
| Repeats / low-complexity runs | kept (a repeat stays a repeat) | homopolymer runs kept; periodic repeats and multi-copy elements destroyed |
| Mirror map (§7.6.4, `mirror_aa`) | [n − end, n − start) | identity |

**Same search-space size.** The decoy frames have the same frames, the same length, the same `*`
and X geometry (identical or mirrored) and therefore exactly the same number of candidate L-mer
positions for every L (verified for L = 7, 8 in every replicate of the attack suite). The I→L fold
is applied before decoying, so neither space can match I/L differently. The circular extension
(first 120 nt appended, §7.2.1) is decoyed with its frame: for reverse the duplicated words stay
duplicated, so distinct-word counts agree exactly; for shuffle the extension is shuffled
independently of the origin region, which adds at most 40 distinct decoy words per frame per
circular replicon (≈ 2 × 10⁻⁵ of a 5-Mb space, conservative, negligible).

### Why reverse for database input

- **Exact under any i.i.d. peptide law.** If Q(w) = Π Q₁(wⱼ) then Q(rev w) = Q(w), and W_D = rev(W_T)
  word for word, so Σ_{W_D} Q = Σ_{W_T} Q *exactly* — for every genome, including repeat-rich and
  low-complexity genomes, and for every difference budget (the m-neighbourhood of rev w is the
  reverse of the m-neighbourhood of w for substitutions; adjacent swaps and the mass-equivalence
  operations are closed under reversal because the operation tables list both orders). Positional
  laws that are not reversal-symmetric (tryptic C-terminal K/R) still agree in expectation when the
  genome's L-mers are reversal-exchangeable; the attack suite confirms it on giapsim genomes.
- **It keeps everything structural.** Composition per window, repeats, low-complexity runs and the
  stop geometry are kept, so the decoy reproduces the de-duplication that repeats impose on the
  target (see the low-complexity results below: shuffle over-counts by 40–70 % there, reverse is
  exact).
- **Deterministic and seed-free**, and its own inverse (`mirror_aa` is an involution).
- It is also the standard choice of proteomic search engines, so database peptides are not
  systematically "close to" reversed sequence (input decoy rows are never mapped, M3-17).

### Why shuffle for de novo input (§8.5 attack 5)

De novo errors frequently reverse a short stretch of a correct peptide (b/y-ion ambiguity). A
peptide with a reversed internal stretch, or a fully reversed one, matches the *reversed* frame near
the true locus far more often than chance, so a reversed decoy would inflate D with correct-locus
peptides and make the de novo FDR uselessly conservative and the de novo precision estimator
(§7.2.5, r_bℓ) biased. A window shuffle breaks all order below 50 aa, so such peptides hit the
decoy only at chance rate, while the 50-aa windows keep local composition (codon bias, hydrophobic
stretches, AT-rich non-coding frames), which is what per-position match probability depends on.

**Window size.** 50 aa is a compromise: windows must be much longer than the longest tolerant
peptide window that matters for chance hits (≤ 40 aa) so that a peptide-length stretch is a
near-random draw from its window, and short enough to keep composition gradients (genes vs.
intergenic, GC skew, transmembrane segments). Windows restart at each stop and skip X, so a window
never mixes two stop-free segments. With window 1 the shuffle is the identity (unit-tested);
with an unbounded window (whole segment) the attack results were indistinguishable at 100 kb
(mean D/T expected-hit ratio 1.003 vs. 1.000 at L = 7; exploratory, 20 genomes).

### What each decoy is a null for

| Matching | reverse | shuffle |
|---|---|---|
| Exact (Route 1) | exact null for any i.i.d. (and any reversal-symmetric) wrong-peptide law, any genome | null for i.i.d. wrong peptides on genomes without heavy repeats; conservative with repeats |
| Tolerant, ≤ m differences (Route 3) | same, because the five cost-1 operations are closed under reversal | same as exact, with the m-neighbourhood |
| Frameshift splits (M3-08) | decoy splits by mirrored adjacency (M3-08) | identity map keeps frame adjacency |

Neither decoy is a null for wrong peptides that share *sequence order* with the genome's own
proteins (see "Limitations").

### Evidence

**Implementation review (M3-03, `src/giap/mapping/decoys.py`).**

- The shuffle swaps slot l of a window with a uniform slot in [0, l] in increasing l (Durstenfeld's
  forward variant of Fisher–Yates), which is uniform over the l! orders by induction; swaps never
  leave the window (`j = s − local[s] + draw`), so per-window multisets, `*` and X positions are
  preserved exactly (property-tested, `tests/property/test_decoys_props.py`).
- Windows restart at `*` but not at X (X skipped), as the issue specifies.
- Determinism: one generator for all frames in index order, one draw per shuffled residue
  (`rng.integers(0, local + 1)`); reverse draws nothing. Seed-identical runs are byte-identical.
  Note: the draws are not `rng.permutation` as the issue text says; this is statistically
  equivalent, but the stream is fixed by this implementation and must not change silently (golden
  tables, §6.6).

**Attack suite (`tests/calibration/attacks/test_attack_genome_decoys.py`).** Fast suite: R = 50
giapsim genomes of 100 kb (seeds 20261014 + i); full suite R = 300. For each genome the *exact*
expected hit probability of an i.i.d. peptide drawn from the genome's six-frame residue
frequencies is computed on the distinct L-mers of the target and of each decoy (no peptide sampling
noise), and three sampled peptide sets are searched.

| Null (fast, R = 50) | reverse D/T | shuffle D/T |
|---|---|---|
| exact expected hits, L = 7 (mean ± SE; range) | 1 (exact) | 1.003 ± 0.002 (0.976–1.022) |
| exact expected hits, L = 8 | 1 (exact) | 1.003 ± 0.002 (0.964–1.037) |
| sampled i.i.d. 7-mers, 10⁶ peptides: T = 1,332 | D = 1,306 (p = 0.63) | D = 1,341 (p = 0.88) |
| sampled tryptic 7-mers (C-term K/R): T = 1,239 | D = 1,261 (p = 0.67) | D = 1,213 (p = 0.61) |
| sampled 8-mers ≤ 1 substitution, 10⁵ peptides: T = 953 | D = 919 (p = 0.45) | D = 944 (p = 0.85) |

p = two-sided binomial test of D vs T given D + T. Per-genome scatter of the shuffle ratio
(SD ≈ 1.2–1.6 % at 100 kb; one 1-Mb genome gave 1.003 / 0.998) is the decoy's own Monte Carlo error, unbiased in the
mean; it shrinks with genome size.

**Low-complexity genome** (same genomes plus 24 simple-sequence repeats of 360 nt — (CAG)n, (GCC)n,
(AT)n, (GA)n, A-runs, hexamer repeats — and ten copies of one 1.5-kb element, ≈ 24 % repetitive;
R = 10 fast):

| | reverse D/T | shuffle D/T |
|---|---|---|
| exact expected hits, L = 7 / L = 8 | 1 / 1 (exact) | 1.40 / 1.69 |
| sampled i.i.d. 7-mers | 221 / 227 | 296 / 227 |
| sampled 8-mers ≤ 1 substitution | 148 / 158 | 221 / 158 |

Repeats shrink the set of distinct target words by ≈ 23 %; the shuffle restores them, so the
shuffled decoy over-counts chance hits. This is conservative (FDR over-estimated, power lost for de
novo input on repeat-rich genomes), never anti-conservative; the gate asserts exactly that.

**Shuffle exactness.** Joint uniformity over two 4-aa windows separated by an X (576 cells,
23,040 draws), the window after a `*`, independence of windows, and position × residue uniformity
in a full 50-aa window: all χ² p > 0.001. The test-author's unit test (one 3-residue window, 600
seeds) is adequate as a smoke test (it catches the classic "swap with any slot" bias already at
n = 3) but tests neither independence between windows nor the X/`*` geometry; the attack suite
covers those.

### Limitations (documented, not CI gates)

Exploratory runs, 1-Mb synthetic genomes, 4 × 10⁵ sampled 7-mers each:

| Violation | D/T reverse | D/T shuffle | Consequence |
|---|---|---|---|
| Wrong peptides share the genome's order-1 (dipeptide) structure: genes and peptides from one asymmetric Markov chain, strength 0.1 / 0.2 / 0.4 | 0.96 / 0.88 / 0.48 | 0.99 / 0.96 / 0.62 | Both decoys under-count chance hits of "protein-like" wrong peptides. Real proteome dipeptide asymmetry is weak (≈ the 0.1 row); database wrong IDs from a related proteome (`search_db = other`) are the realistic case. M3-21 must measure this with real protein sequences (TE8 proxy, UniProt K-12 peptides vs. the Nissle genome). |
| Periodic wrong peptides (dipeptide repeats XYXY…) on the low-complexity genome | 1.00 | 0.88 | Shuffle is anti-conservative for periodic wrong peptides hitting periodic target repeats (period-3 repeats go the other way, 4.2). Recommendation for M3-10: flag `low_complexity` placements (e.g. ≤ 3 distinct residues in the peptide) and report them separately rather than relying on the shuffle decoy. |
| De novo reversed stretches with reverse decoys (attack 5) | inflated D | — | Why de novo uses shuffle; quantified in M3-21. |

### How decoy hits enter the q-values (M3-10)

Completed in "Strata and ranking" and "π̂_wrong in the stratum FDR" below.

### Sign-off (Genome decoys)

Statistician (M3-03 review, 2026-10-14): **approved** for `reverse_frame`, `shuffle_frame`,
`build_decoys`, `mirror_aa` and the `sixframe_index` work file as implemented on `m3/decoys`.
Forwarded: the Markov-structure caveat to M3-21 (real-protein calibration), the low-complexity
recommendation to M3-10, and the repeat over-count of the shuffle (power, not validity).

## Chance hits and the ranking score (M3-09)

### The model

A wrong peptide of length L is modelled as L residues drawn independently from the genome's
six-frame composition f, and every searchable target window as L residues drawn independently from
the same f, independently of the peptide. Then one position matches with probability
s = Σ_a f_a², a window is within ≤ m substitutions with probability
Σ_{j≤m} C(L, j)(1 − s)^j s^(L−j) (j mismatching positions out of L, each mismatching with
probability 1 − s), and the expected number of chance windows is

E(L, m) = N Σ_{j≤m} C(L, j)(1 − s)^j s^(L−j)  (§7.2.5, B.2; `expected_hits_uniform`).

The per-peptide ranking score E_i (`expected_hits`) conditions on the peptide's residues: position i
mismatches a random target residue with probability 1 − f_{a_i}, the mismatch count is
Poisson-binomial, and E_i = N · P(≤ m mismatches). Averaged over i.i.d.-f peptides, E_i equals
E(L, m) exactly (positions are independent and E_f[1 − f_a] = 1 − s); the calibration test
confirms it by simulation at L = 8.

**Composition f** (`residue_freqs`): counts of the 19 residue codes over the target frames, I and L
merged (the frames are I→L folded before anything else, so neither target nor decoy can match I/L
differently); `*` and X are excluded from f. Query letters outside the 19 (U, O, B, Z, J, X) get
f = 0: they never match, so they always count as a mismatch.

**N** (`n_positions`) is the number of six-frame codon positions whose first nucleotide lies in the
replicon: for a circular replicon one per nucleotide and strand (≈ 2 × genome length), the 120-nt
circular extension excluded, but every position near the origin still starts a full window that
reads through the extension (origin-crossing windows are searchable and counted once); for a linear
contig the last ≤ L − 1 positions of each frame have no full window but are counted (negligible).
As the spec defines it, N counts **all** positions, including `*` and X and windows that contain
them. Those windows can never match, so the number of searchable windows is W_L < N (next section).

**Substitutions only.** E counts substitution neighbours. The tolerant matcher also accepts
adjacent swaps and mass-equivalent blocks at cost 1, so the true cost-≤ m neighbourhood is larger
than E assumes (lead default, M3-09 Q1). This does not matter for validity: E is a **ranking
score**, not an FDR estimate. Target and decoy E are computed by the same function from the
peptide sequence and the same genome f (decoys keep the target composition exactly, see Genome
decoys), so E is exchangeable between a wrong peptide's target and decoy hits; the decoy count D
(which uses all operations) calibrates the FDR. An error in E costs power (a worse ordering), never
validity. E also counts hit *positions* with multiplicity, whereas a peptide counts once in the
strata: P(≥ 1 hit) ≤ E, ≈ E when E ≪ 1, which is the range that is accepted.

### Reproduction of B.2

`expected_hits_uniform(L, m, 0.069, 1e7)` against exact rational arithmetic: relative error
≤ 1.3 × 10⁻¹⁵ in all 21 cells; every cell rounds to the B.2 table (e.g. (7, 0) 0.0744635,
(8, 1) 0.559742, (9, 2) 2.36692, (10, 1) 0.00332505, (12, 2) 0.00141835, (14, 1) 1.05295 × 10⁻⁷;
L = 6: 1.08, 88.4, 3.04 × 10³). Monte Carlo (`tests/calibration/test_chance_hits.py`, random
19-letter target of 10⁶ residues with Σf² = 0.069, scaled to N = 10⁷, ≤ 1 substitution):

| L | simulated (± SE) | E(L, 1) | ratio | B.2 ≤ 1 |
|---|---|---|---|---|
| 8 | 0.561 ± 0.006 | 0.5597 | 1.002 | 0.56 |
| 9 | 0.0440 ± 0.0007 | 0.0434 | 1.015 | 0.043 |
| 10 | 0.00339 ± 0.00009 | 0.00333 | 1.020 | 0.0033 |

All within the 20 % gate with SE ≤ 2.7 %. The counting kernel is checked against brute force. This
reproduces the *formula*; it says nothing about real six-frame translations, which the attack below
measures.

### When the independence assumption fails

Attack `tests/calibration/attacks/test_attack_chance_hits.py` (statistician). For a peptide drawn
i.i.d. from the genome's own f, the exact expected number of hit windows is
Σ_{searchable windows t} P(≤ m mismatches | t), mismatch probability 1 − f(t_j) at position j
(closed form for m ≤ 2; no peptide sampling, cross-checked against sampled exact hits). It is
compared with E(L, m) computed by `residue_freqs` + `expected_hits_uniform` on genomes built through
`build_sixframe` (circular, 100 kb, R = 8 per family, seeds from `calibration_base`). Ratios:
stop = W_L / N; cond = observed / (W_L × tail), the error given the stop correction;
total = observed / E = stop × cond. Mean over replicates, m = 0 / 1 / 2:

| Genome | L | stop | cond (m = 0 / 1 / 2) | total (m = 0 / 1 / 2) |
|---|---|---|---|---|
| i.i.d. nucleotides, GC 50 % | 8 | 0.681 | 0.999 / 0.999 / 0.999 | 0.681 / 0.681 / 0.681 |
| | 12 | 0.562 | 0.988 / 0.991 / 0.994 | 0.555 / 0.557 / 0.559 |
| i.i.d. nucleotides, GC 75 % (s ≈ 0.10) | 8 | 0.895 | 0.993 / 0.995 / 0.997 | 0.889 / 0.891 / 0.893 |
| | 12 | 0.847 | 0.989 / 0.992 / 0.994 | 0.838 / 0.840 / 0.842 |
| giapsim (genes, codon structure, AT-rich intergenic) | 7 | 0.748 | 1.019 / 1.012 / 1.006 | 0.762 / 0.757 / 0.752 |
| | 8 | 0.719 | 1.027 / 1.018 / 1.011 | 0.739 / 0.732 / 0.727 |
| | 10 | 0.666 | 1.047 / 1.036 / 1.025 | 0.697 / 0.690 / 0.683 |
| | 12 | 0.618 | 1.068 / 1.055 / 1.043 | 0.660 / 0.652 / 0.644 |
| giapsim + 24 simple-sequence repeats of 360 nt (≈ 9 %) | 7 | 0.766 | 1.82 / 1.42 / 1.21 | 1.39 / 1.09 / 0.93 |
| | 8 | 0.740 | 2.54 / 1.80 / 1.41 | 1.88 / 1.33 / 1.04 |
| | 10 | 0.691 | 6.7 / 3.9 / 2.5 | 4.6 / 2.7 / 1.7 |
| | 12 | 0.647 | 22.6 / 11.7 / 6.4 | 14.6 / 7.6 / 4.1 |

(Replicate SE of cond on the i.i.d. genomes ≈ 0.3–0.6 % at L = 8–12; the 1 % deviations at L = 12
are within 2 SE, confirmed with R = 30 at GC 75 %: 0.998 ± 0.002 at L = 8, 0.999 ± 0.005 at L = 12.)

Three effects, with their directions:

1. **Stops and X** (factor W_L / N ≈ (1 − p_stop)^L, here 0.56–0.91; p_stop ≈ 3/64 for GC 50 %): N counts positions that
   cannot start a match. E **over**-states the absolute chance count — conservative — by a factor 1.1–1.8,
   more for longer L and AT-rich genomes (more stop codons in non-coding frames). With i.i.d.
   codons this is the only error (cond = 1 at any composition, including a strongly skewed GC 75 %
   genome: the model handles composition skew through f).
2. **Heterogeneous composition and dipeptide structure.** Coding, antisense and intergenic frames
   have different residue laws g_k; for an i.i.d.-f peptide E[Π f(t_j)] = Σ_k w_k (Σ_a f_a g_k(a))^L
   ≥ s^L by Jensen, and codon/dipeptide structure inside genes adds to it. E **under**-states the
   conditional count by 2–7 % (growing with L, shrinking with m). On realistic genomes the stop
   factor dominates, so E stays conservative overall (total 0.64–0.76).
3. **Low complexity.** Simple-sequence repeats put long runs of a few abundant residues into the
   target (poly-A, poly-P, (AR)n, …); windows there have Π f(t_j) ≫ s^L. E **under**-states the
   absolute count, strongly and increasingly with L (×1.9 at L = 8, ×15 at L = 12, exact matches).
   By construction the excess comes from peptides that are themselves rich in the repeat
   residues (an exact or 1-substitution match to a repeat window needs a near-repeat peptide):
   those hit far more often than their E_i says, while typical peptides are barely affected. This is the same structure that made the shuffled decoy
   anti-conservative for periodic wrong peptides (Genome decoys › Limitations, D/T 0.88) and that
   the reverse decoy preserves exactly.

**Consequence.** Since E is only a ranking score, none of these errors affects FDR validity as long
as the decoy is a valid null (reverse: exact for any i.i.d. law on any genome; shuffle: see Genome
decoys). They affect (a) power, through the ordering: the stop factor varies with L, so pooling
lengths within a stratum (G3) ranks long peptides slightly too pessimistically relative to short
ones (≈ 3–5 % per residue against a factor ≈ 1/s ≈ 14 per residue in E itself: immaterial);
low-complexity peptides are ranked too optimistically, which is the case for the `low_complexity`
placement flag recommended to M3-10; (b) any *absolute* reading of E (a reported "expected chance
hits" column, or the §7.2.5 table read as a property of a real genome): E over-states real
genomes by 1.3–1.8× and under-states low-complexity peptides by up to an order of magnitude. If E is
ever reported as an absolute expectation, use W_L instead of N and flag low-complexity peptides.

### Numerical accuracy of the Poisson-binomial

`giap.stats.poisson_binomial` (shared with M4-02) computes the pmf by the standard DP
pmf_{i+1}[j] = pmf_i[j](1 − p_i) + pmf_i[j − 1] p_i, O(n²) time, O(n) memory. Every step is a
convex combination of non-negative numbers, so there is no cancellation: each entry has relative
error ≤ ≈ 2n·u (u = 1.1 × 10⁻¹⁶) whatever its size, and tiny tails stay exact to relative precision
rather than to absolute 10⁻¹⁶. The cdf and sf are summed directly from the pmf (never 1 − the
other), as §8.3 requires. Checks (review, 2026-10-14):

- brute-force enumeration of all 2ⁿ outcomes, n = 0–12, 260 random p vectors including
  p^20-skewed ones: max relative pmf error 2.5 × 10⁻¹⁵; cdf and sf for every k ∈ [−1, n + 1] within
  10⁻¹⁴ absolute;
- exact rational arithmetic, n = 60, p ∈ [0, 0.1]: max relative error 6.0 × 10⁻¹⁶ over the whole
  pmf, including the extreme tail P(X = 60) = 4.75 × 10⁻⁸³;
- n = 5,000 (M4-02 scale): Σ pmf − 1 = 4 × 10⁻¹⁶, 36 ms;
- E at L = 40, m = 0 stays positive (7 × 10⁻⁴⁵ for a uniform-f poly-A query at N = 10⁷).

Limits (outside the peptide range, documented only): entries below ≈ 10⁻³⁰⁸ underflow to 0 (n ≳ 230
at p = 0.05 for P(X = n)); `expected_hits_uniform` uses exact integer binomial coefficients and
raises `OverflowError` for L ≳ 1,000 with m ≈ L/2. The vectorised in-place update reads its
right-hand side before assigning, so there is no aliasing.

### Review findings (M3-09, `src/giap/mapping/chance.py`, `src/giap/stats/poisson_binomial.py`)

- B.2: all 21 cells exact to 1.3 × 10⁻¹⁵; examples 1–8 of the issue pass; Poisson-binomial correct
  against brute force and exact arithmetic (above). Edge cases: m ∉ [0, L] and p ∉ [0, 1] or NaN
  raise; empty p gives pmf [1]; k < 0 gives cdf 0 and sf 1; k > n gives cdf 1 and sf 0.
- **Origin off-by-two (minor, requested fix).** `_core_residues` keeps a minus-strand codon when its
  first (5′-on-minus, highest-coordinate) nucleotide is < replicon length. The two minus-strand
  codons whose 5′ nucleotide is replicon position 0 or 1 lie wholly in the extension in the frame
  coordinates and are dropped, although they are unique origin-crossing codons, so a circular
  replicon gets N = 2·length − 2 instead of 2·length (checked at length 9, 10, 11, 300 and 100 kb).
  Fix: keep a minus-strand codon when its lowest coordinate, `first_nt − 2`, is < length. Effect on
  E and f: relative 10⁻⁵, none on rankings or q-values (target and decoy share N). The attack
  accepts either value until fixed.
- N includes `*`/X positions as the spec defines it; documented above with its direction
  (conservative ×0.56–0.91), not a defect.

### Sign-off (Chance hits)

Statistician (M3-09 review, 2026-10-14): **approved** for `expected_hits_uniform`,
`expected_hits`, `residue_freqs` and `giap.stats.poisson_binomial` as implemented on
`m3/chance-hits`, with the origin off-by-two to be fixed (non-blocking: no effect on any statistical
output). Forwarded to M3-10: E is a ranking score only (never report it as an absolute expectation
without the W_L correction); low-complexity peptides are ranked too optimistically by E, which
reinforces the `low_complexity` placement flag.

## π̂_wrong in the stratum FDR (M3-10)

Scope: how a π̂_wrong that varies between peptides enters FDR_m(e) (§7.2.5 step 3). In practice
this is de novo input in stratum 0, where π̂ is per score bin (step 4, M3-11). Database input in
stratum 0 (one π̂ for all) and every stratum m ≥ 1 (π̂ = 1) have a constant π̂, and the rule
below reduces exactly to the spec formula there. This section **replaces** lead decision M3-10 Q2
("mean π̂ of the targets with E ≤ e"), which is anti-conservative.

### Model and assumptions

Fix a stratum m and a decoy type. Every mapped query peptide j of that decoy type is searched
against target and decoy (lead decision M3-10 Q1). Let π_j = P(j wrong | its score bin) and
π̂_j its estimate.

- **A1 (decoy hit ⇔ chance target hit of the same peptide).** For a wrong peptide j, (target hit at
  cost m, its E) and (decoy hit at cost m, its E) have the same law: p_j(e) = P(hit with E ≤ e)
  for both (Genome decoys; E is computed by the same function for both, M3-09).
- **A2 (decoy hits do not depend on correctness).** The decoy hit of j has the same law p_j(e)
  whether j is correct or wrong: a correct peptide's decoy hits are chance hits too.
- **A3 (π̂ is a peptide-level label).** π̂_j depends on the peptide's score bin only, not on its
  hits; it is the same for the peptide's target and decoy hit.

p_j(e) may vary freely between peptides (length, composition) and may be correlated with π_j.

### Derivation

The number of false targets with E ≤ e is V(e) = Σ_j 1{j wrong} 1{chance target hit of j, E ≤ e},
so by A1 E[V(e)] = Σ_j π_j p_j(e). Candidate estimators of V(e) and their expectations:

| Estimator of V(e) | Expectation | Unbiased when |
|---|---|---|
| Q2 as written: π̄_T(e) · D(e), π̄_T(e) = mean π̂ of targets with E ≤ e | π̄_T(e) Σ_j p_j(e) | π̄_T(e) = Σ π_j p_j / Σ p_j; fails in general |
| (a) π̄_S · D(e), π̄_S = mean π̂ over all searched peptides | π̄_S Σ_j p_j(e) | Cov_j(π_j, p_j(e)) = 0 |
| (b) W(e) = Σ_{decoy hits j with E ≤ e} π̂_j | Σ_j π_j p_j(e) (A2, A3) | always (given A1–A3 and π̂ = π) |

Why Q2 fails: a correct peptide always has its target hit, a wrong one only by chance, so the
targets with small E are mostly correct peptides from the good (low-π) bins and π̄_T(e) is far
below the decoy-hit-weighted π̄ that D(e) needs; D counts chance hits of *all* searched peptides,
most of which come from the bad bins. The bias is largest at the top of the list, i.e. at the
operating point α = 0.01. Why (a) is not enough: chance-hit rates vary by a factor ≈ 14 per
residue of length (M3-09), and low-score de novo bins plausibly hold shorter peptides, so
Cov(π, p) > 0 is the realistic case and (a) is then anti-conservative by Σπp / (π̄ Σp).

The "+1" (§7.2.5, B.6) guards the top of the list against D = 0 by pre-counting one decoy. With
weights, one more decoy hit can weigh up to max π̂, so the pre-count must be c = max π̂. Using the
mean π̂ instead (c = π̄_S) left a small excess on B.6-sized lists with very unequal bins (mean FDP
0.0113 ± 0.0004 at q = 0.01, N = 150, π = 0.02 / 0.9, R = 3,000).

### Signed rule

Within stratum m and decoy type, one list ordered by E (unchanged), at every target E = e:

  FDR_m(e) = (W_m(e) + c) / T_m(e),  W_m(e) = Σ over peptides j with a decoy hit at cost m and
  E_D,j ≤ e of π̂_j,  c = max π̂ over all searched peptides of that decoy type in the stratum's
  input,

and q_i = min over e ≥ E_i of FDR_m(e), clipped to 1; ties as before (a decoy tied with a target
counts). Definitions:

- π̂_j is the peptide's own π̂ by the step-4 route, whether it contributes a target, a decoy or
  both: database m = 0 the one `pi_wrong_database` value; de novo m = 0 its `pi_denovo` (1 if
  missing); m ≥ 1 always 1. A decoy-only peptide carries its own π̂.
- "Searched peptides of that decoy type" = every peptide of that source in the costs table that
  enters the decoy search, i.e. all rows of the source except `same_db` exact peptides in
  stratum 0 (they leave both stratum-0 lists). If the set is empty, c = 1. Any de novo peptide
  without π̂ makes c = 1.
- Constant π̂ = π: W = π D, c = π, FDR = π (D + 1)/T, identical to the spec and to
  `tda_qvalues(−E_T, −E_D, π, plus_one=True)`.
- Monotonicity: q is non-decreasing in E within the stratum; adding a decoy hit, raising any π̂_j
  of a decoy, or raising c never lowers a q.
- `mapping_fdr` rows (per decoy type × length class × m, no threshold): `decoys` = D (count, as
  before), `fdr_stratum` = (W + c)/T with W summed over the class's decoy hits and the stratum's
  c (NA when T = 0), `pi_wrong` = (W + c)/(D + 1), the effective π̂, so `fdr_stratum` =
  `pi_wrong` (D + 1)/T still holds and `pi_wrong` equals the constant π̂ for database input and
  m ≥ 1.

Suggested interface (implementer's choice, the semantics above bind):
`stratum_qvalues(e_target, e_decoy, pi_decoy, pi_plus)`, `pi_decoy` one weight per decoy,
`pi_plus` = c (scalar in (0, 1]).

### Evidence

Statistics-only simulation (generator of `tests/calibration/test_placement_fdr_lists.py`: every
peptide draws a chance decoy hit with probability p, wrong peptides a chance target hit from the
same law, correct peptides a low-E target hit), de novo stratum 0, N = 4,000, R = 300, mean FDP ±
SE at q = 0.01 / 0.05 / 0.10 (\* = above q + 2 SE); n = mean accepted at q = 0.01.

| Scenario | Q2 (mean over targets) | (a) mean over searched | **(b) signed** | n at 0.01 (b) |
|---|---|---|---|---|
| 10 bins π 0.02–0.9, p = 0.25 (calibration test) | 0.0143\* / 0.069\* / 0.130\* | 0.0098 / 0.050 / 0.100 | 0.0095 / 0.050 / 0.100 | 1,063 |
| 2 bins π 0.02 / 0.9, p = 0.25 | 0.0354\* / 0.119\* / 0.175\* | 0.0097 / 0.050 / 0.100 | 0.0095 / 0.050 / 0.100 | 1,069 |
| constant π = 0.3 | 0.0100 / 0.050 / 0.097 | identical | identical | 1,932 |
| 10 bins, p = 0.05 + 0.4 π | 0.0182\* / 0.086\* / 0.159\* | 0.0126\* / 0.065\* / 0.130\* | 0.0094 / 0.050 / 0.100 | 949 |
| 2 bins, p = 0.05 + 0.4 π | 0.0540\* / 0.168\* / 0.253\* | 0.0163\* / 0.085\* / 0.172\* | 0.0093 / 0.049 / 0.100 | 796 |
| 10 bins, π̂ = 1.5 π | 0.0099 / 0.049 / 0.094 | 0.0072 / 0.037 / 0.073 | 0.0069 / 0.036 / 0.073 | 878 |
| 2 bins, π̂ = 1.5 π (capped at 1) | 0.0312\* / 0.108\* / 0.171\* | 0.0086 / 0.045 / 0.089 | 0.0085 / 0.044 / 0.089 | 1,005 |
| 10 bins, π̂ = 0.5 π | 0.0286\* / 0.130\* / 0.176\* | 0.0199\* / 0.100\* / 0.176\* | 0.0196\* / 0.100\* / 0.176\* | 1,497 |
| 2 bins, π̂ = 0.5 π | 0.0615\* / 0.175\* / 0.175\* | 0.0198\* / 0.100\* / 0.175\* | 0.0195\* / 0.100\* / 0.175\* | 1,496 |

Small lists (R = 3,000), signed rule: 2 bins N = 150: 0.0000 / 0.044 / 0.095 (n = 0 at 0.01: with
c ≈ 0.9 no list of ≤ 90 targets can reach 0.01); N = 300: 0.0054 / 0.048 / 0.097 (n = 33);
10 bins N = 300: 0.0025 / 0.046 / 0.096 (n = 17). With c = mean π̂ instead: 0.0113\* / 0.0099 /
0.0089 at q = 0.01 for these three, n = 29 / 70 / 65.

**Misestimation.** Under-estimating π by 50 % doubles the FDP for every rule (E[W] = ½ E[V]):
validity needs π̂ ≥ π, which is why step 4 uses 1 − the *lower* 95 % bound of the de novo
precision. Over-estimating by 50 % is conservative for (b) (FDP ≈ q/1.5); Q2 is anti-conservative
even then in the two-bin case.

**Attack** (`tests/calibration/attacks/test_attack_placement_pi.py`, §8.5; N = 2,000 unless noted;
R = 50 fast / 300 full): two_bins, two_bins_corr, ten_bins_corr, two_bins_small (N = 300),
ten_bins_over (π̂ = 1.5 π), constant must stay at level through `placement`; a power check asserts
that Q2 and (a) exceed q + 2 SE at q = 0.01 and 0.05 on two_bins_corr (full suite: Q2 0.053 /
0.168, (a) 0.016 / 0.085). The signed rule, run on the same generators through a reference
implementation, gave (full suite, q = 0.01 / 0.05 / 0.10): two_bins 0.0093 / 0.050 / 0.101,
two_bins_corr 0.0088 / 0.049 / 0.100, ten_bins_corr 0.0086 / 0.050 / 0.100, two_bins_small
0.0055 / 0.049 / 0.100, ten_bins_over 0.0064 / 0.036 / 0.073, constant 0.0100 / 0.050 / 0.096.

### Alternative considered: separate lists per π̂ value (score bin)

Running stratum 0 separately per score bin (the Q2 alternative) applies the spec formula with a
constant π̂ in each list. It is also valid in every scenario above (it needs A1–A3 only within a
bin) and has more power at α = 0.01 on large inputs with unequal bins (N = 4,000: 1,336 vs 1,063
accepted with 10 bins, 1,969 vs 1,069 with 2 bins), less at q ≥ 0.05 and on small inputs (10 bins
N = 300 at 0.01: 55 vs 17 for (b), 21 vs 0 at N = 150). It splits a stratum into several lists,
so q is no longer monotone in E within a stratum (§7.2.5 step 3 and the M3-10 acceptance
criterion), and `mapping_fdr` would need a bin dimension. It is therefore a spec change for the
lead (ADR), not the default; if adopted, the attack above applies unchanged.

### Forwarded

- **Database input, mean PEP.** The spec's single π̂ = mean PEP is estimator (a): unbiased only if
  PEP is uncorrelated with the chance-hit rate. Short peptides have both higher PEP and many more
  chance hits, so mean PEP can under-state the false-placement share (10-bin correlated scenario:
  0.0126 at q = 0.01; database PEPs are small, so the absolute effect is bounded by the PEPs). The
  same weighting (Σ PEP_j over decoy hits + max PEP) removes the bias but departs from the spec
  text: lead decision / ADR; M3-21 should measure Cov(PEP, chance-hit rate) on A2.
- **π̂ from the same decoys.** The de novo π̂ (M3-11) uses decoy-match rates r_bℓ of the same
  shuffle; a bin with many chance decoy hits gets a lower precision and a higher π̂, so W and π̂
  are positively dependent: conservative.

### Tests that must change (test-author)

1. `tests/unit/test_placement_fdr.py::test_per_target_pi_uses_mean_over_targets_up_to_e` → weights
   on decoys: targets E = 1, 2, 3; one decoy E = 2.5 with π̂ 0.5; c = 0.6:
   FDR = 0.6/1, 0.6/2, (0.5 + 0.6)/3 → q = [0.3, 0.3, 0.366667].
2. `…::test_denovo_stratum0_uses_per_peptide_pi` → de novo m = 0, length 7: targets (E, π̂) =
   (1, 0.5), (2, 0.6), (3, 0.3); decoy-only peptide E = 2.5, π̂ 0.5; c = max = 0.6:
   q = [0.3, 0.3, 0.366667]; row (`7`, shuffle): targets 3, decoys 1, `fdr_stratum` = 1.1/3 =
   0.366667, `pi_wrong` = 1.1/2 = 0.55. (With the old π̂ 1.0 for the second target: c = 1,
   q = [0.5, 0.5, 0.5], `fdr_stratum` 0.5, `pi_wrong` 0.75.)
3. `…::test_constant_pi_equals_m4_tda_qvalues` and the module docstring (FDR(e) line): call with
   decoy weights `np.full(ed.size, pi)` and c = pi; equality with `tda_qvalues` unchanged.
   Every other unit test (examples 1–10, ties, clipping, π̂ = 1 when missing) is unchanged: their
   π̂ is constant within the list.
4. `tests/property/test_placement_fdr_props.py`: the `stratum` strategy draws one π̂ per *decoy*
   and a scalar c ∈ (0, 1] (≥ max decoy π̂ to mirror `placement`); `test_constant_pi_equals_…`
   as in 3; `test_input_order_invariant` permutes decoy weights with the decoys; add: raising any
   decoy π̂ or c never lowers a q. The other properties (q ∈ [0, 1], monotone in E, adding a decoy
   never lowers q, scale invariance, all `placement` properties) stay.
5. `tests/calibration/test_placement_fdr_lists.py::test_denovo_per_bin_pi_controls_fdr`: generator
   unchanged (expected to pass, 0.0095 / 0.050 / 0.100 at R = 300); docstring: drop "most
   favourable case for … mean pi over targets with E <= e" and cite this section.
6. Spec text `docs/specs/M3-10-placement-fdr.md` (bullets "π̂ per target" and `mapping_fdr` rows,
   Q2) and lead decision 2 in `docs/specs/M3-index.md`: lead updates them to this rule.

### Sign-off (π̂_wrong in the stratum FDR)

Statistician (M3-10 pre-implementation decision, 2026-10-14): **signed** for the rule above
(decoy-weighted W_m(e) + c, c = max π̂ over the searched peptides of the decoy type). Lead
decision M3-10 Q2 (mean π̂ over targets with E ≤ e) is **rejected** as anti-conservative and
blocks merge if implemented. The remaining M3-10 sections (strata and ranking, π̂_wrong sources,
same-strain input, calibration) are signed with the implementation review.

## Strata and ranking (M3-10)

- **Search.** Every mapped query peptide (non-decoy, q ≤ `input.q_value_max` or q missing) is
  searched against the decoy of its own source (database → reverse, de novo → shuffle; both built
  for mixed input) exactly as against the target (exact for all, tolerant for all ≥ 10 aa; lead
  decision M3-10 Q1). Each peptide counts once at its best target cost m_T and, independently,
  once at its best decoy cost m_D (separate search, not competition: a peptide may be a target in
  stratum 0 and a decoy in stratum 1). Decoy hits are de-duplicated on their mirrored target
  interval, never placed and never written to `peptide_hits.tsv`.
- **Lists.** One list per decoy type × m ∈ {0, 1, 2}: T = peptides with m_T = m, D = peptides with
  m_D = m, both ranked by E ascending (M3-09; target and decoy E from the same function and the
  same genome f). Length classes are **pooled** inside a list (G3): E already carries the length
  (≈ ×14 per residue), so pooling only orders long against short peptides by their chance rate;
  per-class lists would have small T and a +c that dominates (B.6). The per length class × m rows
  of `mapping_fdr.tsv` are descriptive: T and D of the class, `fdr_stratum` = (W_class + c)/T_class
  with the stratum's c, no threshold; they are not the q-values.
- **FDR and q.** FDR_m(e) = (W_m(e) + c)/T_m(e) at every target E (rule and derivation in
  "π̂_wrong in the stratum FDR"); q_i = min over e ≥ E_i; ties: all entries with E ≤ e count, so a
  decoy tied with a target counts against it and tied targets get one q; clipped to [0, 1] by
  M4-01's `monotone_q`, so with constant π̂ the q-values equal `tda_qvalues(−E_T, −E_D, π̂)`
  (unit- and property-tested). Accept at q ≤ `stats.alpha.placement` (0.01, inclusive).
- **Families (§8.1).** Each list is one FDR family: "placed (exact)" = stratum 0, "placed with ≥ 1
  difference" = strata 1 and 2 with π̂ = 1, separately per decoy type. Each accepted set has FDR
  ≤ α; their union has marginal FDR E[V]/E[R] ≤ α exactly and FDR ≤ α up to the weighting of the
  ratios (each list pays its own +c, which makes the union slightly conservative). A peptide accepted in stratum m is accepted at
  all its loci of cost m: multi-locus ambiguity is uniqueness (M3-14), not FDR.
- **Effect size.** The per-peptide effect measure next to `placement_q` is the ranking score E
  (expected chance windows, conservative ×1.3–1.8 on real genomes, see M3-09) and m_T; the list
  quantities T, D, W, c are in `mapping_fdr.tsv`. A target-decoy q-value has no CI of its own; its
  uncertainty is the +c pre-count (B.6).

**B.6.** Stratum 1, π̂ = 1 (c = 1): 150 targets before the first decoy give FDR = (0 + 1)/150 =
0.0067, all 150 accepted at 0.01; 50 targets and no decoys give 1/50 = 0.02, none accepted. With
per-peptide de novo π̂ the floor is c/T with c = max π̂, so a de novo stratum-0 list needs
T ≥ c/α targets before the first decoy (90 at c = 0.9): small de novo inputs are conservative by
design (attack `two_bins_small`: FDP 0.0063 at 0.01).

## π̂_wrong sources (M3-10)

| Input, stratum | π̂_j | Population |
|---|---|---|
| database, m = 0 | **superseded by ADR-0019** (see "Cross-strain input and run collapse"): the exact-match estimate of the peptide's length class. Until then: one value, mean PEP if any PEP; else max q ≤ `input.q_value_max`; else 1 | the database peptides **of the stratum-0 lists**, i.e. all database rows except `same_db` exact peptides (ADR-0019 also drops contaminants; see "Same-strain input") |
| de novo, m = 0 | `pi_denovo` of the peptide's score bin (1 − lower 95 % precision bound, M3-11); missing → 1 | per peptide |
| any, m ≥ 1 | 1 | — |

Mean PEP is estimator (a) of the derivation: unbiased when PEP is calibrated and uncorrelated with
the chance-hit rate (forwarded to M3-21). Max q estimates the FDR of the whole input list at its
own threshold (≤ `q_value_max`), an estimate of the wrong share of the peptides that enter mapping
(test `max_q_route`: FDP 0.0020 for a true wrong share 0.008, conservative). π̂ = 0 (all PEPs exactly 0) is allowed by the formula and accepts
every target; engines report PEP 0 only through underflow, so a mean over thousands of peptides
is never 0 in practice. c = max π̂ over the same population (1 if empty or any de novo π̂
missing).

## Same-strain input (M3-10)

§7.2.5 step 5: database peptides from a search against this genome's own proteins
(`search_db = same`) with an exact target hit get `placement_q = 0` and flag `same_db`, and leave
both stratum-0 lists. Rationale: the engine's own target-decoy already controls the identification
error against exactly these sequences, and an exact hit of an own-protein peptide is its source
locus; the remaining question is uniqueness (M3-14). Their decoy hits at m ≥ 1 still count in D
(conservative: a correct peptide's decoy hit is a chance hit, A2). Tolerant hits of `same` peptides
go through steps 1–4 like any other.

Two conditions are needed for the rest of the list to stay valid (attack
`tests/calibration/attacks/test_attack_placement_same_db.py`):

1. **Database input only.** `search_db` is filled from `input.search_db` for every row, de novo
   rows included (`io/peptides.py`). A de novo peptide is not from a protein search; the shortcut
   must require `source = database`, otherwise every exact de novo hit gets q = 0.
2. **π̂ on the listed population.** The `same_db` exact peptides have left the stratum-0 lists, so
   they must also leave the mean PEP (and the max-q fallback). They are selected for being correct
   (an exact hit in their own genome) and carry small PEPs, so averaging them in lowers π̂ for the
   `other` peptides that remain in the list. This changes nothing for single-`search_db` inputs
   (all `other`: nothing removed; all `same`: the stratum-0 database target list is empty), only
   for mixed inputs such as A1 (K-12 search) and A2 (Nissle's own proteins) in one run, which the
   per-row `search_db` column allows.

| Attack (R = 50 fast / 300 full; FDP over list-accepted peptides) | `placement` at 9ec53d1 (before fix) | corrected rule (reference; = `placement` at 031378c) |
|---|---|---|
| `denovo_labelled_same`: de novo, π 0.02 / 0.9, every row `search_db = same` | 0.175 / 0.175 / 0.175 (full 0.176) | 0.0097 / 0.049 / 0.097 (full 0.0093 / 0.049 / 0.099) |
| `mixed_same_other`: 2,000 `same` exact (PEP 0.001) + 2,000 `other` (π = PEP = 0.3) | 0.0196 / 0.097 / 0.097 (full 0.0199 / 0.096 / 0.096) | 0.0098 / 0.050 / 0.097 (full 0.0100 / 0.050 / 0.096) |

(mean FDP at q = 0.01 / 0.05 / 0.10; SE ≤ 0.0016 fast, ≤ 0.0007 full.) Required fixes:
`same_exact = is_db & (search_db == "same") & (m_target == 0)` and
`pi_wrong_database(df[is_db & ~same_exact], …)`.

An unresolved `search_db = auto` reaching `placement` is treated as `other`: valid (conservative,
the peptide is FDR-controlled), but it loses the `same_db` shortcut; `auto` must be resolved to
`same`/`other` before mapping (M3-17).

## Cross-strain input and run collapse (ADR-0019, M3-25)

**Status:** signed off by the statistician (M3-25 review of #124, 2026-10-20; ADR-0019 decisions
1–5). Decision and evidence:
`docs/adr/ADR-0019-stratum0-pi-hat-exact-share.md` (accepted by Khoa 2026-10-20; Slurm 1860123,
1860138, 1860180). This section supersedes the mean-PEP row of "π̂_wrong sources" for database
input in stratum 0 and resolves the item forwarded from ADR-0010 (π̂ correlated with the
chance-hit rate).

**Why mean PEP fails.** Two defects, both measured in M3-21:

1. *Cross-strain input* (`search_db = other`). Some **correct** identifications report the
   database strain's sequence across a strain difference. They have no exact true locus, so an
   exact chance hit elsewhere is a false stratum-0 placement, while their PEP measures
   identification error only and is small. Mean PEP counts only the identification errors among
   the peptides without an exact origin and under-estimates the false placements by about
   (π_wrong + δ)/π_wrong, δ = the spanning share of the correct identifications.
2. *Run collapse.* §5.3 collapses a peptide's runs with `pep = min`. The mean of min-collapsed
   PEPs is not a peptide-level wrong share: with 3 runs M3-21 measured π̂ ≈ 0.0026 against a true
   wrong share ≈ 0.024. No PSM-level aggregation (min, mean, max, product) is calibrated in
   general, so placement stops reading PEP instead of changing §5.3's aggregation, which ranking,
   detection and evidence still use.

**The rule.** Population P: all database rows except `same_db` exact peptides (unchanged from
M3-10) and except contaminants (M3-25 decision 7 — they are not genome-derived, so they carry no
information about the chance-hit rate). Per length class ℓ of `mapping.length_classes`:

- n_ℓ = peptides of P in ℓ; k_ℓ = those with no exact target hit (m_T ≠ 0 or no hit at all);
  r_ℓ = **share of n_ℓ** with an exact decoy hit;
- U_ℓ = one-sided 95 % Clopper–Pearson upper bound of k_ℓ/n_ℓ (1 when k_ℓ = n_ℓ), mirroring
  step 4's "1 − lower 95 % bound";
- π̂_ℓ = 1 if r_ℓ = 1, else min(1, U_ℓ/(1 − r_ℓ)).

This is the database analogue of the de novo precision estimator with c_l = 1: a peptide with an
exact origin always matches, one without matches at the decoy rate. Every peptide of P carries
its class's π̂ into ADR-0010's FDR_0(e) = (W_0(e) + c)/T_0(e), with c = max π̂_ℓ over the classes
with n_ℓ > 0. A peptide outside every class weighs 1 and does **not** enter c; an empty class has
no weight. Edge cases: r_ℓ = 0 gives π̂ = U_ℓ, and m < r (fewer exact target than exact decoy
matches) clips to 1. m ≥ 1 is unchanged (π̂ = 1). Why per class and not pooled: the no-origin
share correlates with the chance-hit rate, which changes by ≈ ×14 per residue of length; an
oracle that knows the true **pooled** no-origin share still fails the short-peptide and padding
scenarios (2.19×, 2.02×).

**Evidence** (failing cells of 18 list-level scenarios, δ × π × q; ADR-0019 Evidence):
mean PEP fails 6–14 cells per scenario, a pooled exact share 0–13, the per-class rule **0**
everywhere, including PEP deflated 10× (the run collapse) and the length-correlated inputs.
Power at δ = 0 relative to mean PEP is 0.984–1.0 at N = 4,000 and 0.85–0.97 at N = 300 (π̂ = 1
instead gives 0.43–0.84 and 0.25). At genome scale the estimated/true false count moves from
0.03–0.46 under mean PEP to 0.96–2.4 under this rule; acceptance at q = 0.01 is unchanged, small
or binding lists lose up to ~15 % power.

**Limitation and the run warning.** If the upstream pipeline kept only genome-matching peptides,
k is ~0 by construction and no genome-based rule can be calibrated — mean PEP included (FDP 0.061
at q = 0.01, δ = 0.2). `placement` therefore raises the run warning `input_prefiltered` when k
falls below the 0.001 quantile of Binomial(n, mean PEP × (1 − r)) over P pooled; this is the one
place the row-level mean PEP is still read, and it is read directly, not through the old
`pi_wrong_database` (whose max-q and 1.0 fallbacks would misfire here). The detector inherits
PEP's weakness in one direction: a **deflated** PEP (the §5.3 min over runs, defect 2 above, or
an engine that under-reports) lowers the expected count n × mean PEP × (1 − r) towards the
realised k, so the test loses power and a genuinely pre-filtered input can pass unwarned. It
never fires spuriously for that reason — too small an expectation only hides the warning — and
the π̂ it guards does not read PEP at all. An inflated PEP is the conservative direction.

**Reporting.** `mapping_fdr.tsv` columns are unchanged, with `pi_wrong` = (W + c)/(D + 1) now
reflecting the class π̂. Read the per-class π̂ from the class table, **not** from `mapping_fdr`:
a class with no decoy hit has W = 0 and D = 0, so its `pi_wrong` row is c, the largest class
estimate, not that class's π̂_ℓ. The class table is the `placement` step's `stratum0_classes`
output (`mapping/stratum0_classes.parquet`, columns `length_class, n, k, r, upper, pi_wrong`),
an optional input of `context_stats`, which turns it
into `detection_summary.tsv` QC rows `stratum0_unmapped_share` = k/n (the §7.5.1
genome-inconsistent rate of this population), one for the run and one per length class, each
per-class row carrying that class's `upper` as `ci_high`. It sits
beside `genome_inconsistent_rate`, which measures another population (C9 `other`, non-contaminant,
over accepted exact placements) and is unchanged.

## Calibration (M3-10)

Criterion (§8.4 "FDR"): mean FDP over R replicates ≤ q + 2 SE at q ∈ {0.01, 0.05, 0.10}; R =
`replicates_fast.placement` (50) per PR, `replicates_full.placement` (300) nightly; seeds from
`tolerances.yaml`.

**List level** (`tests/calibration/test_placement_fdr_lists.py`, through `placement`): correct
peptides get a low-E target hit, every peptide a chance decoy hit with probability p and wrong
peptides an independent chance target hit from the same E law (A1–A2 hold by construction).
Review run, 9ec53d1, R = 50, mean FDP ± SE at 0.01 / 0.05 / 0.10:

| Scenario | FDP | n accepted at 0.01 |
|---|---|---|
| known π = 0.3, m = 0 | 0.0095 ± 0.0003 / 0.049 ± 0.001 / 0.097 ± 0.001 | 1,898 |
| known π = 0.3, m = 1 (π̂ = 1, conservative) | 0.0026 / 0.015 / 0.030 | 1,001 |
| calibrated PEP ~ Beta(0.6, 1.4) | 0.0103 ± 0.0004 / 0.050 ± 0.001 / 0.096 ± 0.001 | 1,959 |
| max-q route (true share 0.008) | 0.0020 / 0.0020 / 0.0020 | 3,977 |
| de novo, 10 bins | 0.0086 ± 0.0005 / 0.049 ± 0.001 / 0.099 ± 0.001 | 1,037 |

**Attacks** (`test_attack_placement_pi.py`, through `placement`, R = 50): two_bins 0.0080 / 0.050 /
0.101; two_bins_corr 0.0085 / 0.049 / 0.099; ten_bins_corr 0.0089 / 0.050 / 0.099;
two_bins_small 0.0063 / 0.045 / 0.095; ten_bins_over 0.0066 / 0.036 / 0.072; constant 0.0093 /
0.048 / 0.095 (SE ≤ 0.0017 at 0.01). The rejected rules still fail (power check passes). Same-strain
attack: see the table above (failed at 9ec53d1, passes at 031378c).

**Genome level, exact route** (`tests/calibration/test_placement_fdr_genome.py`): 100-kb giapsim
genomes, tryptic peptides of annotated genes as correct, 20,000 i.i.d. 7–8-mers from the genome's
own f as wrong, database input, reverse decoys, PEP = the known wrong share. R = 50: 0.0102 ± 0.0007
/ 0.0236 ± 0.0009 / 0.0236 ± 0.0009, 592 accepted at 0.01; 719 chance target hits of wrong
peptides in total (the null is exercised). At 0.05 and 0.10 every target is accepted, so the FDP
equals the wrong share among targets (0.024 < q). At 0.01 the mean is within 1 SE of q, which is
expected: with constant π̂ and an exact reverse decoy the estimator is unbiased up to the +1.
Tolerant strata, de novo errors with shuffle decoys, attack 5 (reversed stretches), 1-Mb genomes
and real-protein wrong peptides (TE8 proxy) are M3-21.

## Implementation review (M3-10, `src/giap/mapping/placement_fdr.py`, 9ec53d1)

Checked against ADR-0010, the final section of `docs/specs/M3-10-placement-fdr.md` and this note:

- `stratum_qvalues`: T(e) and W(e) by `searchsorted(side="right")` (ties counted against the
  target, tied targets share one q), `monotone_q` from M4-01 (running minimum from the loose end,
  clip to [0, 1]); NaN targets → NaN, NaN decoys ignored; weights and c validated in [0, 1].
  Correct. Reduction to `tda_qvalues` with constant π̂: unit- and property-tested.
- W and c: per decoy type and stratum; W sums each decoy hit's own π̂ (de novo `pi_denovo`,
  missing → 1, clipped to [0, 1]; database the one value; m ≥ 1 → 1); c = max π̂ over **all**
  peptides of the source in the costs table, including peptides without any hit, minus `same_db`
  exact peptides at m = 0; 1 for m ≥ 1 or an empty population. Correct.
- `mapping_fdr` rows: decoy type × m × length class plus `all`, `fdr_stratum` = (W + c)/T (NaN at
  T = 0), `pi_wrong` = (W + c)/(D + 1), `n_accepted` per class; peptides outside every length class
  count only in `all`. Correct.
- Flags joined with `;`: matches §9 ("`;`-separated controlled vocabulary"). `low_complexity` (≤ 3
  distinct residues, I = L) is the flag this note recommended (M3-03, M3-09), but it is not in the
  §9.4 placement vocabulary, which may be extended "only through the schema": the lead must add it
  to the schema/ADR (M3-16) before outputs are written. Not a statistical blocker.
- **Blocking at 9ec53d1, fixed in 031378c:** the two same-strain conditions above. The fix is exactly
  `same_exact = is_db & (search_db == "same") & (m_target == 0)` and
  `pi_wrong_database(df[is_db & ~same_exact], …)`. Re-run at 1dafac0 (R = 50 / 300, q = 0.01 / 0.05 /
  0.10): `denovo_labelled_same` 0.0097 / 0.049 / 0.097 (full 0.0093 / 0.049 / 0.099),
  `mixed_same_other` 0.0098 / 0.050 / 0.097 (full 0.0099 / 0.050 / 0.096). The list, `placement_pi`
  and genome numbers above are unchanged (single-`search_db` inputs are unaffected). The
  `low_complexity` flag was accepted by ADR-0011.

### Sign-off (Strata and ranking, π̂_wrong sources, Same-strain input, Calibration)

Statistician (M3-10 implementation review, 2026-10-14): the four sections above are **signed** as
the method. The implementation is **approved** at 031378c (re-verified at 1dafac0): the `same_db` shortcut
requires `source = database`, the database π̂ is computed over the stratum-0 population, and
`test_attack_placement_same_db.py` passes (at 9ec53d1 the FDP at q = 0.01 was 0.175 and 0.0196).
Forwarded to M3-21: Cov(PEP, chance-hit rate) on A2; the list-level calibration of tolerant strata
with real decoy hits; attack 5.

## Genome-level calibration (M3-21)

Suites: `tests/calibration/test_placement_fdr_sim.py`, `test_denovo_precision_sim.py`,
`attacks/test_attack5_reversed_stretches.py`; truth oracle `tests/helpers/placement_truth.py`
(origins per row type, overlap rule of lead decision 2, unique guard, cost-≤ m_T+1 mutant).

**What "exercised" means.** A list (decoy type × m) is exercised in a replicate by the number of
peptides in it that are false placements (no accepted locus overlaps a true origin), at any q.
A gate on a list with no false peptides passes trivially and is a **smoke check**: it shows the
chain runs and accepts, not that the FDR is controlled.

**At 100 kb the placement lists are not exercised enough to test the bound.** In the fast tier
(R = 50, job 1860994) the database-shuffled (π 0.01, 0.05), related other-strain and related
db-strain suites had 0 false peptides in every list. The de novo error-model suite and attack 5 had a mean of
0.14 and 0.12 per replicate in shuffle m = 0 and ≤ 0.02 elsewhere. These gates are therefore smoke checks
(`claimed = ()`). The genome-level evidence for stratum 0 comes from two places only:

- the **stressed stratum-0 arm** (`test_database_stressed_stratum0_placement_fdr`, full tier at
  1 Mb): short wrong database peptides (7/8-mers of protein sequences that are absent from the
  genome, as in b2sim `--mode short`) are injected test-locally with PEP set to the true wrong
  share. This arm claims reverse m = 0. With δ = 0 (exact PEP) its FDP gate
  is expected to pass. At 100 kb (K = 2,000, R = 50, job 1860994): 1.6 false peptides per
  replicate, mean FDP 0.0009 at every q. At 1 Mb (K = 400, R = 300, job 1860995): 2.85 false per
  replicate, mean FDP 0.00016 at every q, about 17,960 accepted. The fast arm also asserts the
  estimate check mean(Ŵ0) ≥ mean(V0) − 2 SE: Ŵ0 2.44 vs V0 1.6, because the injected share
  dominates the pooled π̂ at 100 kb. At 1 Mb the estimate check is a separate test,
  `test_database_stressed_stratum0_estimate_covers_false_placements` (`calibration_full`,
  R = 100, ΣV0 ≥ 100). It was a strict xfail "ADR-0019 implementation pending" until M3-25
  implemented the ADR: mean PEP gave Ŵ0 0.57 vs V0 2.85. The cause is the length-correlated no-origin share. 75 % of the injected
  peptides are 7-mers, and the 7-mer class has a no-origin share of 0.22 against 0.031 pooled,
  so a pooled π̂ times all decoy hits under-counts. ADR-0019's per-class π̂ corrects it; the
  statistician's reference rule gives Ŵ0 3.65 vs V0 3.50 over 6 replicates (review of #113,
  round 2). **Post-ADR, through `placement`** (M3-25, job 1861430, R = 100): Ŵ0 **4.306** vs
  V0 **2.75**, SE 0.191, ΣV0 = 275 (floor 100) — the check passes, with a ratio of 1.57 against
  0.20 under mean PEP. The FDP stays far below q under mean PEP only because the stratum is large;
- the **related db-strain estimate check**
  (`test_related_dbstrain_stratum0_estimate_covers_false_placements`, full tier, 1 Mb, 5 SNV/kb,
  R = 100). It requires ΣV0 ≥ 30 false stratum-0 placements and mean W0 ≥ mean V0 − 2 SE (SE of V0 − W0 over lists). Observed
  (job 1860995, after the B1 origin fix): ΣV0 = 79, mean V0 0.79, mean W0 0.150, SE 0.099. It was
  a strict xfail until ADR-0019 (exact-match π̂ per length class) was implemented in M3-25. Correct
  identifications that span a strain difference have no exact origin but small PEP, so mean PEP
  under-counts their exact chance hits (ADR-0019 Context, defect 1). **Post-ADR, through
  `placement`** (M3-25, job 1861430, R = 100): Ŵ0 **1.551** vs V0 **0.79**, SE 0.107,
  ΣV0 = 79 (floor 30) — the check passes, with a ratio of 1.96 against 0.19 under mean PEP.
  Both post-ADR ratios (1.57 and 1.96) sit inside ADR-0019's expected 0.96–2.4 band; the excess
  over 1 is the Clopper–Pearson bound plus c, and c is set by the most sparsely populated class
  (a class with k = n forces π̂_ℓ = 1), so most of the conservative margin can come from one
  small class — read `stratum0_classes` to see which.

**FDR evidence at m ≥ 1 rests on the list-level suites** (`test_placement_fdr_lists.py`,
`attacks/test_attack_placement_pi.py`), not on these genome suites.

**Attack 5.** Its FDP gate is a smoke check as well. Its content is the decoy-count documentation:
de novo reversed stretches searched with reverse decoys inflate the m = 0 decoy count D
(mean 9.1 reversed vs 0.5 shuffled per replicate). That is why de novo input uses shuffle decoys.
