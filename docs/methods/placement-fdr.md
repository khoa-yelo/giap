# Placement FDR: genome decoys, chance hits and per-stratum q-values

Implements the genome side of GIAP_SPEC.md §7.2.1 (decoy translations), §7.2.5 (placement FDR, G3),
§8.1 rows "Peptide placed at a locus", §8.2 "decoys mimic targets", §8.5 attack 5 and Appendix
B.2/B.6. Issues: M3-03 (decoys, this section), M3-09 (chance hits), M3-10 (strata and q-values),
M3-21 (genome-level calibration). Engine for the q-values: `docs/methods/target-decoy.md`.

Sections: Genome decoys (M3-03, **signed off below**); Chance hits and the ranking score (M3-09,
*to be completed*); Strata and ranking (M3-10, *to be completed*); π̂_wrong (M3-10, *to be
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
- *To be completed with M3-09/M3-10:* the E-ranking justification; π̂_wrong (mean PEP, max q, de
  novo lower bound; per-bin π̂ as the mean over targets with E ≤ e, lead decision M3-10 Q2);
  `same_db`; the B.6 example; the list-level calibration (exchangeable target/decoy E, FDR ≤ q +
  2 SE) and the genome-level calibration and attack 5 (M3-21).

### Sign-off (Genome decoys)

Statistician (M3-03 review, 2026-10-14): **approved** for `reverse_frame`, `shuffle_frame`,
`build_decoys`, `mirror_aa` and the `sixframe_index` work file as implemented on `m3/decoys`.
Forwarded: the Markov-structure caveat to M3-21 (real-protein calibration), the low-complexity
recommendation to M3-10, and the repeat over-count of the shuffle (power, not validity).
