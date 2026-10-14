# M4-12: Detection tiers: picked gene FDR, MAYU, two-peptide, protein list

Milestone M4 · track G · depends on: M4-01, M4-02, M4-11 (calibration needs M4-09) · estimate ≈ 380
code lines (+ ≈ 400 test lines) · **can start now** (fixtures) · gate item: "gene-level FDR calibrated
(picked, MAYU)"

## Spec sections implemented

- §2.1 decision 4: "Three tiers, chosen from what the input carries: decoy rows → picked gene-level
  target–decoy FDR; q-values or PEPs only → MAYU-style gene FDR; sequences only → genes with ≥ 2
  distinct peptides, single-peptide genes `unestimated`. The ≥ 2-peptide tier is always reported too."
- §7.3.2 table (quoted): `decoy_picked` "decoy rows with `protein_ids` exist and ≥ 80 % of them pair
  with a gene"; `mayu` "otherwise, q-values or PEPs exist"; `two_peptide` "otherwise, peptide sequences
  only … single-peptide genes get `detection_status = single_peptide_unestimated`"; `protein_list`
  "proteins only, no peptides | detected = a mapped input protein | `unestimated`".
  "The tier is chosen once per run from the input and recorded in the manifest and in `genes.tsv`."
  "For downstream statistics, single-peptide genes under `two_peptide` count as not detected in the
  primary analysis and as detected in a sensitivity rerun."
- §7.3.2 picked steps 1–4: "s_g is the best accepted target peptide of gene g, by −log₁₀ q (ties broken
  by lower PEP, then higher `score`)." "Strip the decoy prefix … map the resulting target protein to
  its gene … Unpaired decoys are counted and reported; if more than 20 % of decoys are unpaired, fall
  back to the `mayu` tier with a warning." "A gene is detected at gene q ≤ `stats.alpha.identification`."
- §7.3.2 MAYU (M4-02 computes it): "Detected at q ≤ 0.01." "≥ 2-peptide variant … reported in
  `detection_summary.tsv` only."
- §7.3 "Eligible genes for detection statistics are the non-partial protein-coding rows."
- §9.3 `detection_status` enum{detected, not_detected, single_peptide_unestimated, partial_excluded};
  `detection_summary.tsv` section `tier` "the tier and its FDR curve points".
- §8.4 suite "Gene FDR (picked, MAYU) | false peptides on random genes, size-proportional | FDR | 50 | 300".

## Public interface

`src/giap/genetable/detection.py` (pure)
```python
Tier = Literal["decoy_picked", "mayu", "two_peptide", "protein_list"]

def protein_gene_map(genes, peptides, peptide_hits, protein_hits) -> pd.DataFrame: ...
    # protein_id → gene, method ∈ {locus_tag, protein_hit, peptide_vote}
def pair_decoys(peptides: pd.DataFrame, pmap: pd.DataFrame, decoy_prefixes: Sequence[str]) -> pd.DataFrame: ...
    # one row per decoy peptide: peptide, q_value, pep, score, gene (NA if unpaired)
def choose_tier(peptides, peptide_hits, protein_hits, pairing: pd.DataFrame) -> tuple[Tier, list[str]]: ...
    # (tier, warnings)
@dataclass(frozen=True, slots=True)
class Detection:
    tier: Tier
    table: pd.DataFrame        # locus_tag, detected, detection_status, gene_q_value, detection_tier
    summary_rows: pd.DataFrame # detection_summary section "tier" rows (M4-16 writes them)
    warnings: tuple[str, ...]
def detect(genes, evidence, accepted, peptides, peptide_hits, protein_hits, covariates,
           cfg: GiapConfig) -> Detection: ...
```
- **Protein → gene map** (priority): (1) the protein ID is a C2 `locus_tag`; (2) C12 `gene` of its
  rank-1 hit; (3) peptide vote: among accepted unique target peptides whose C9 `protein_ids` contain
  the protein, the gene with most distinct peptides (a tie → unmapped).
- **Decoy pairing:** a C9 row with `is_decoy` and q ≤ `q_value_max` (when q exists); strip the first
  matching prefix of `input.decoy_prefixes` from each of its `protein_ids`; paired when the stripped
  IDs map to exactly one gene, else unpaired.
- **Score key** (targets and decoys alike): the lexicographic tuple (q ascending, PEP ascending, score
  descending), NA last; converted to one float by dense rank over the union of target-best and
  decoy-best keys (higher = better) before `picked_fdr` — equivalent to "−log₁₀ q, ties by PEP then
  score" without infinities at q = 0.
- **Tier choice:** `decoy_picked` if ≥ 1 decoy row has `protein_ids` and paired/(paired + unpaired)
  ≥ 0.8; if decoys exist but < 0.8 pair → `mayu` with warning `decoy_pairing_low`; else `mayu` if any
  accepted target has q or PEP; else `two_peptide` if C11 has rows; else `protein_list`.
- **Status:** non-eligible (partial) → `partial_excluded`, `detected` false, q NA. Then per tier:
  `decoy_picked`/`mayu`: `detected` = gene q ≤ `stats.alpha.identification`; a gene whose decoy wins
  → `not_detected`, q NA; no accepted peptide → `not_detected`, q NA. `two_peptide`: `detected_2pep` →
  `detected`; exactly 1 accepted → `single_peptide_unestimated` (`detected` false); 0 → `not_detected`.
  `protein_list`: a C12 rank-1 hit with `gene` → `detected`.
- MAYU inputs: eligible genes; sizes = `n_observable_peptides` (M4-10); P(t) = distinct accepted
  peptides; F from PEP if every accepted peptide has PEP, else t·|P(t)|.
- `summary_rows` (section `tier`, §9.3 long format): `group_type = "tier"`, `group = <tier>`, metrics
  `n_detected`, `n_detected_2pep`, `n_decoys_paired`, `n_decoys_unpaired`; and the FDR curve as
  `group_type = "q_threshold"`, `group ∈ {0.001, 0.005, 0.01, 0.02, 0.05, 0.1}`, metrics `n_detected`,
  `fdr` (and `fdr_2pep` for MAYU).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | M4-01 example 8 encoded as genes G1…G6 with best target / decoy keys | G1, G2 q 0.5; G4 0.667; G6 0.75; G3, G5 `not_detected` q NA |
| 2 | 150 genes with best targets above every decoy, 1 decoy below all (B.6 shape) | all 150 q = 1/150 → `detected` at 0.01 |
| 3 | decoys: 10 rows, 9 with `rev_<locus tag>`, 1 with `rev_UNKNOWN` | paired 9/10 = 0.9 ≥ 0.8 → `decoy_picked`; `n_decoys_unpaired = 1` |
| 4 | decoys 10 rows, 7 pair | `mayu`, warning `decoy_pairing_low` |
| 5 | decoy `protein_ids = "rev_P1;rev_P2"`, P1 → G, P2 → H | unpaired (two genes) |
| 6 | no decoys, q-values present; B.3 configuration (5,000 eligible genes, uniform sizes) | tier `mayu`; curve row at t covering all 3,000 genes reports `fdr ≈ 0.0067` |
| 7 | sequences only; G with 2 accepted, H with 1, K with 0 | tier `two_peptide`; G `detected`; H `single_peptide_unestimated`, `detected` false; K `not_detected`; all q NA |
| 8 | proteins only: C12 rank-1 hit on G, none on H | tier `protein_list`; G `detected`, H `not_detected` |
| 9 | partial gene with 5 accepted peptides | `partial_excluded`, `detected` false |
| 10 | target best key (q 0.001, PEP 0.01) vs paired decoy (q 0.001, PEP 0.01, same score) | tie → decoy wins → `not_detected` |

## Edge cases

- Zero eligible genes → tier still chosen; empty curve; warning `no_eligible_genes`.
- q = 0 for many targets: handled by ranks (no −log₁₀ 0).
- Contaminant-prefixed protein IDs (`input.contaminant_prefixes`) are never mapped to genes.
- Protein IDs that are protein groups `A;B` (M1-03 keeps the group ID): each member is looked up.
- MAYU with PEP present for some peptides only → t·|P(t)| for the run (warning `pep_partial`).

## Acceptance criteria

- [ ] `tests/unit/genetable/test_detection.py`: examples 1–10.
- [ ] Property: adding a decoy never turns a `not_detected` gene into `detected`; `detected` implies
  `n_accepted ≥ 1`; tier is a function of the input columns only (row permutation invariant).
- [ ] **Calibration** `tests/calibration/test_gene_fdr.py` (suite `gene_fdr`, R = 50 fast / 300 full,
  seeds from the harness): M4-07 layout (2,000 genes fast, 4,500 full) → M4-08 → M4-09 with
  π_wrong ∈ {0.01, 0.05} → M4-10/11 → `detect`. FDP = (detected genes with `present = false`) /
  max(1, detected). Criterion `check_fdr` at q ∈ {0.01, 0.05, 0.10} (detection rerun at each q),
  separately for `decoy_picked` (decoys emitted) and `mayu` (decoys dropped; with PEP and without PEP).
- [ ] Recorded (not gated): empirical FDP of the `two_peptide` tier on the same replicates, written to
  the calibration log (it is `unestimated` by design).
- [ ] Fast suite for this file < 4 min on 4 CPUs.

## Methods note (statistician sign-off)

`docs/methods/detection-tiers.md`: tier choice, protein→gene map, decoy pairing and the 80 % rule,
picked competition at gene level (Savitski 2015), MAYU use (cross-reference `mayu-gene-fdr.md`),
two-peptide and protein-list as `unestimated`, and the calibration evidence (FDP tables per tier).

## Questions for the lead

1. Decoy rows filtered by the same engine threshold (q ≤ `q_value_max`) as targets before picking?
   **Proposed default:** yes, symmetric with target acceptance (targets also need placement, which
   decoys cannot have; the statistician checks the calibration covers this asymmetry).
2. The 80 % pairing threshold has no config key. **Proposed default:** a module constant now; a
   `detection.min_decoy_pairing` key later only via ADR (defaults.yaml is protected).
3. In peptide tiers, do mapped long proteins (C12) also count as detection? **Proposed default:** no;
   they only feed the protein→gene map and protein-level abundance. §7.3.2 uses proteins only in
   `protein_list`.
4. The peptide-vote step of the protein→gene map is not in the spec (needed when engine IDs are
   UniProt accessions and no `proteins.faa` is given, as in pilot arm A1). **Proposed default:** keep it.
