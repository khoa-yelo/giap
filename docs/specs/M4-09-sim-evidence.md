# M4-09: giapsim evidence: peptides, identification errors, decoys and placed hits

Milestone M4 · track F (simulator) · depends on: M4-08 · estimate ≈ 380 lines in `tests/sim/giapsim/`
(tooling) · **can start now**

Produces the evidence half of the M4-06 contract (C9 `peptides.norm.parquet`, C10
`abundance.norm.parquet`, and the M3-format final tables C11 `mapping/peptide_hits.tsv`, C12
`mapping/protein_hits.tsv`, C13 `mapping/unmapped.tsv`, 1-based) directly from the
expression truth, as if M1 preflight and M3 mapping had run perfectly (every target peptide placed at
its true locus). Mapping errors are M3's business and calibrated there; this generator isolates
identification error at the gene level.

## Spec sections implemented

- §11.3 Peptides: "tryptic digest (7–30 aa, ≤ 1 missed cleavage); peptide detectability logistic in
  peptide features and log abundance; log-normal intensities." "Database-search errors: false
  identifications at rate π_wrong, drawn from a shuffled proteome or from a related 'other strain'
  proteome. Scores come from a two-component mixture, from which q-values and PEPs are computed
  consistently. Decoy rows (reversed sequences) carry null scores."
- §8.4: "Gene FDR (picked, MAYU) | false peptides on random genes, size-proportional | FDR".
- §8.7: "false precursors pile up on otherwise undetected genes".
- §7.2.6 uniqueness classes and §7.3.2 decoy pairing ("Strip the decoy prefix from a decoy row's
  protein ID and map the resulting target protein to its gene").

## Interface (`tests/sim/giapsim/evidence.py`)

```python
@dataclass(frozen=True)
class EvidenceConfig:
    pi_wrong: float = 0.01                 # share of false target identifications among accepted-looking rows
    false_source: Literal["random_gene", "unplaceable"] = "random_gene"
    n_decoys: Literal["matched"] | int = "matched"       # matched = as many decoys as false targets (TDA)
    mixture: tuple[float, float] = (0.0, 3.0)            # N(mu0,1) wrong/decoy, N(mu1,1) correct scores
    emit: frozenset[str] = frozenset({"q_value", "pep", "is_decoy", "intensity", "protein_ids"})
    protein_ids: Literal["locus_tag", "foreign"] = "locus_tag"   # foreign = "UP_<n>" accessions (cross-strain)
    multi_identical_pairs: int = 0          # duplicate gene pairs (IS copies) sharing identical sequences
    shared_frac: float = 0.0               # share of peptides that also occur in a second, different gene
    tolerant_frac: float = 0.0             # share of true peptides placed via route=tolerant (n_diff 1)
    proteins_only: bool = False            # protein_list tier: emit C12 rows, no peptides
    abundance_table: Literal["none", "ibaq", "maxlfq"] = "none"

@dataclass(frozen=True)
class Evidence:
    peptides: pd.DataFrame       # C9 rows (targets + decoys)
    peptide_hits: pd.DataFrame   # C11 rows (targets only, accepted placements)
    protein_hits: pd.DataFrame | None
    abundance: pd.DataFrame | None
    truth: pd.DataFrame          # peptide, is_false, gene_true, decoy_gene (paired gene)
    gene_truth: pd.DataFrame     # locus_tag, present (≥ 1 correct peptide), n_true_peptides

def simulate_evidence(layout: Layout, expr: Expression, cfg: EvidenceConfig,
                      rng: np.random.Generator) -> Evidence: ...
def write_evidence_tables(ev: Evidence, work_dir: Path) -> None: ...
```
Rules:
- True peptides: digest each `detected_true` gene (trypsin/P, 7–30 aa, ≤ 1 missed cleavage); each
  observed with probability `sigmoid(a + b·abundance_log2_true + peptide noise)`; at least one peptide
  for every `detected_true` gene. Intensity ~ LogNormal(log(2^abundance) + N(0, 0.5²)).
- False targets: n_false = round(π_wrong/(1 − π_wrong) · n_true). With `random_gene`, each lands on gene
  g with probability ∝ n_observable_g (all genes, detected or not) and takes one of g's tryptic peptides
  not already observed, so it maps; with `unplaceable`, a random tryptic-like sequence that maps nowhere
  (does not enter C11). Scores: correct ~ N(mu1, 1), false and decoy ~ N(mu0, 1).
- Decoys: reversed tryptic peptides (C-terminal K/R kept) of genes drawn ∝ n_observable;
  `protein_ids = "rev_" + <that gene's protein ID>`; decoy rows never enter C11.
- PEP = π0·f0(s) / (π0·f0(s) + (1 − π0)·f1(s)) with π0 = π_wrong; q(s) = mean PEP over targets with
  score ≥ s (consistent Bayesian FDR). Fields not in `emit` are written as NA (to drive tier choice).
- C11 rows: `locus_class = cds_inframe` (`sorf` for sORFs, `pseudogene` for pseudogenes), `route`
  exact (or tolerant with `n_diff = 1`, flag `strain_difference`), `placement_q = 0`, `uniqueness`
  `unique`, or `multi_identical`/`shared` as configured (one row per gene hit).
- `abundance_table = "ibaq"`: C10 rows per protein ID with log2 iBAQ from true abundance + noise;
  `"maxlfq"`: same values, declared type maxlfq (relative).

## Examples

| # | Call | Expected |
|---|---|---|
| 1 | default config, default layout and expression, seed 1 | every `detected_true` gene has ≥ 1 correct peptide in C11; C9 decoy rows count = false target count (matched) |
| 2 | `pi_wrong = 0.05` | false share among target C9 rows within 0.05 ± 0.01 |
| 3 | example 1, false peptides' genes | share landing on genes with `detected_true = False` ≈ share of n_observable held by those genes (± 0.05) |
| 4 | decoy row `protein_ids = "rev_GSIM_00015"` | stripping `rev_` gives a locus tag present in C2 |
| 5 | `emit = {"q_value"}` | C9 `pep` all NA, `is_decoy` false and no decoy rows; drives the `mayu` tier |
| 6 | `emit = set()` (sequences only) | C9 has only peptide strings and counts; drives `two_peptide` |
| 7 | `proteins_only = True` | C9 and C11 empty; C12 has one row per `detected_true` gene (`class identical`, `rank 1`); drives `protein_list` |
| 8 | `multi_identical_pairs = 2` | 2 gene pairs with identical `aa`; their peptides appear twice in C11 with `uniqueness = multi_identical` |
| 9 | q/PEP consistency | over targets, mean PEP of rows with q ≤ 0.01 is ≤ 0.01; realised false share among those within 0.01 + 3·SE |
| 10 | same seed twice | byte-identical tables |

## Edge cases

- A gene with no tryptic peptide of 7–30 aa (e.g. a 25-aa sORF) cannot be detected by peptides;
  `detected_true` for it is forced to false and recorded in `gene_truth`.
- `pi_wrong = 0` → no false targets, no decoys when `matched`.
- `n_decoys` integer overrides matched (to test unpaired-decoy fallbacks, combine with
  `protein_ids = "foreign"` for unmappable decoys).

## Acceptance criteria

- [ ] `tests/sim/test_evidence.py`: examples 1–10.
- [ ] Output passes `giap.genetable.contract.validate_table` for C9–C12 (cross-check, once M4-06 merges).
- [ ] No import of `giap` in `giapsim`.
- [ ] 4,500 genes, default config: < 3 s.

## Methods note

Described in `docs/methods/detection-tiers.md` (M4-12) as the gene-FDR null; statistician reviews.

## Questions for the lead

1. Is it acceptable that the M4 gene-FDR calibration bypasses real mapping (perfect placement)?
   **Proposed default:** yes; M3 calibrates placement FDR separately, and M4-21 reruns one replicate
   end-to-end through real M3 once it lands.
2. Default π_wrong for calibration: {0.01, 0.05}? **Proposed default:** fast suite 0.01 and 0.05 (half
   the replicates each); nightly adds 0.10.
3. M3-19 (`giapsim/proteomics.py`) defines the same score mixture, PEP and q model on the full
   simulator. **Proposed default:** one shared helper `giapsim/scores.py` (mixture → PEP, q); whichever
   of M3-19 / M4-09 merges first creates it and the other imports it, so both FDR calibrations use
   one error model.
