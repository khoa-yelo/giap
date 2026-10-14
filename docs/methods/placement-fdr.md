# Placement FDR: genome decoys, chance hits and per-stratum q-values

Implements the genome side of GIAP_SPEC.md §7.2.1 (decoy translations), §7.2.5 (placement FDR, G3),
§8.1 rows "Peptide placed at a locus", §8.2 "decoys mimic targets", §8.5 attack 5 and Appendix
B.2/B.6. Issues: M3-03 (decoys, this section), M3-09 (chance hits), M3-10 (strata and q-values),
M3-21 (genome-level calibration). Engine for the q-values: `docs/methods/target-decoy.md`.

Sections: Genome decoys (M3-03, **signed off below**); Chance hits and the ranking score (M3-09,
**signed off below**); Strata and ranking (M3-10, *to be completed*); π̂_wrong (M3-10, *to be
completed*); Same-strain input (M3-10, *to be completed*); Calibration incl. attack 5 (M3-10/M3-21,
*to be completed*).

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

### How decoy hits enter the q-values (M3-10, *to be completed*)

- Every mapped query peptide is searched against the decoy of its own source (database → reverse,
  de novo → shuffle; both built for mixed input), exactly as against the target (exact for all,
  tolerant for all ≥ 10 aa; lead decision M3-10 Q1). Each peptide counts once at its best decoy
  cost m_D, independently of its target cost m_T (separate search, not competition).
- Decoy hits are de-duplicated on their mirrored target interval (`mirror_aa`, modulo the
  replicon), are never placed and never written to `peptide_hits.tsv`.
- Stratum m: T = peptides with m_T = m, D = peptides with m_D = m, ranked by E (M3-09; target and
  decoy E computed by the same function from the same genome frequencies, since decoys have the
  target's composition). FDR_m(e) = π̂_wrong (D_m(e) + 1)/T_m(e), q by `tda_qvalues(−E_T, −E_D, π̂)`;
  ties counted against the target. Lists are per `decoy_type`.
- *To be completed with M3-10:* π̂_wrong (mean PEP, max q, de
  novo lower bound; per-bin π̂ as the mean over targets with E ≤ e, lead decision M3-10 Q2);
  `same_db`; the B.6 example; the list-level calibration (exchangeable target/decoy E, FDR ≤ q +
  2 SE) and the genome-level calibration and attack 5 (M3-21).

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
