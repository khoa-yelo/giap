# M2-17: Codon usage: CAI, tAI, ENC, Fop, GC3, codon-pair bias, `at_2_10`, `rare_codons_2_10`

Milestone M2 · depends on: M2-04 (genes, ribosomal set, tRNA counts) · estimate ≈ 220 code lines
(tests not counted) · capabilities: none (codon-bias 0.5.0 is MIT and in the main env)

## Spec sections implemented

- §7.1.1: "Codon usage | codon-bias 0.5.0 (MIT, PyPI) | main GIAP env, in-process | library calls |
  CAI, tAI, ENC, Fop, GC3, codon-pair bias".
- §7.1.6: "**codon-bias.** CAI (reference = ribosomal-protein genes); tAI (tRNA gene copy numbers);
  ENC; Fop (optimal codons from the reference set); GC3; codon-pair bias (reference = all CDS).
  **Computed directly:** AT content of codons 2–10; `rare_codons_2_10`: codons whose share within
  their synonymous family is < 10 % genome-wide."
- §6.2 `codon_usage` | py | "genes, `trna_counts.tsv` → `codon.parquet` (uses codon-bias in-process; MIT)".
- §5.2 `annotation.translation_table: 11`, `codon: {cai_reference: ribosomal_proteins}`. §9.3
  `genes.tsv`: `cai`, `tai`, `enc`, `fop`, `gc3`, `cpb`, `at_2_10`, `rare_codons_2_10`; §9.4 gene flag `cai_reference`.

## Public interface

`src/giap/annotate/codon_usage.py`; step `codon_usage`.
```python
@dataclass(frozen=True, slots=True)
class CodonFeatures:
    cai: float | None; tai: float | None; enc: float | None; fop: float | None
    gc3: float | None; cpb: float | None; at_2_10: float | None; rare_codons_2_10: int | None

def scored_sequence(gene: Gene) -> str | None
    # nt without the start and the stop codon (codons 2 … n−1); None if pseudogene, partial,
    # length % 3 != 0, or < 10 codons
def synonymous_shares(seqs: Iterable[str], table: int = 11) -> dict[str, float]   # codon → share within its amino-acid family
def rare_codons(shares: Mapping[str, float], cutoff: float = 0.10) -> frozenset[str]
def at_2_10(nt: str) -> float | None        # AT fraction of codons 2–10 (27 nt) of the full CDS nt
def rare_2_10(nt: str, rare: frozenset[str]) -> int | None
def gc3(seq: str) -> float
def codon_usage(genes: Sequence[Gene], trna_counts: Mapping[str, int], table: int = 11) -> tuple[dict[str, CodonFeatures], tuple[str, ...]]
    # (locus_tag → features, run flags)
```
codon-bias calls (all with `genetic_code=11`; the library default is 1):
`CodonAdaptationIndex(ref_seq=<ribosomal scored seqs>, genetic_code=11)`;
`TrnaAdaptationIndex(tGCN=DataFrame(anti_codon, GCN), prokaryote=True, genetic_code=11)`;
`EffectiveNumberOfCodons(genetic_code=11)`; `FrequencyOfOptimalCodons(ref_seq=<ribosomal>, genetic_code=11)`;
`CodonPairBias(ref_seq=<all eligible CDS>, genetic_code=11)`. `codon.parquet` columns: `locus_tag`
+ the 8 fields; genes with `scored_sequence is None` get all-NA except `at_2_10`/`rare_codons_2_10`
when ≥ 10 codons exist.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `at_2_10("ATG" + "AAA"*9 + "TAA")` | `1.0` |
| 2 | `at_2_10("ATG" + "GCC"*9 + "TAA")` | `0.0` |
| 3 | CDS with 8 codons | `at_2_10 → None`, `rare_2_10 → None` |
| 4 | `gc3("GCCGCGGCA")` | `2/3 = 0.666667` |
| 5 | genome-wide Leu codon counts CTG 50, TTA 13, TTG 13, CTT 11, CTC 10, CTA 3 (total 100) | shares CTA 0.03, CTC 0.10 → `rare_codons` contains `CTA`, not `CTC` (< 10 % strict) |
| 6 | codons 2–10 `CTA CTG CTA AAA …` with `rare = {CTA}` | `rare_2_10 = 2` |
| 7 | genes whose ribosomal set is empty | `cai = fop = None` for every gene; flags `("cai_reference_missing",)` (Q2) |
| 8 | `trna_counts = {}` | `tai = None` for all genes |
| 9 | ribosomal gene `rplB` | carries gene flag `cai_reference` (exported in `codon.parquet` as `cai_reference=true`) |
| 10 | two runs on the same genes | identical `codon.parquet` bytes |

## Edge cases

- Single-codon families (Met ATG, Trp TGG) have share 1.0 and are never rare; stops are excluded.
- Alternative starts (GTG/TTG) are excluded from scores because codon 1 is dropped.
- Codons with `N` or other IUPAC letters are skipped by `synonymous_shares`; `at_2_10` is A+T over
  A+C+G+T bases only (ambiguous bases excluded from numerator and denominator).
- sORFs with ≥ 10 codons are scored; pseudogenes never.

## Acceptance criteria

- [ ] `tests/unit/test_codon_usage.py`: examples 1–10.
- [ ] `tests/golden/test_codon_k12.py`: on the recorded K-12 full JSON, CAI of `rplB`, `rpsA`,
  `tufA` lies in the top decile, and median CAI of ribosomal genes > median of all genes
  (`pending_sample` until M2-10).
- [ ] A test asserts `genetic_code=11` reaches every codon-bias constructor (monkeypatched spy).

## Questions for the lead

1. Score codons 2…n−1 (drop start and stop) for every index. **Proposed default:** yes, consistently.
2. No ribosomal-protein genes found (fragmented draft). **Proposed default:** CAI/Fop `NA` and a run
   flag `cai_reference_missing` (new run flag in the §9.4 vocabulary, added via the schema in M2-20).

## Lead decisions (after the test-author's pass)
1. Q1/Q2: the proposed defaults are accepted.
2. The test-author's resolutions are accepted: "< 10 codons" counted on the full CDS; CAI reference = the
   scored ribosomal genes (`cai_reference`); genome-wide shares and the CPB reference over all scored
   genes, unseen codons in a present family have share 0.0; `codon.parquet` columns `locus_tag`, the 8
   fields, `cai_reference`; NaN written as null.
3. The codon-bias library must pass scripts/ci/check_imports.py (licence) and already be in uv.lock; if
   not, stop and tell the lead (no new dependencies).
