# M3-10: Placement FDR (G3): strata, π̂_wrong, q-values, `mapping_fdr` rows

Milestone M3 · depends on: M3-09 (E_i); consumes the hit schema of M3-04 (synthetic hit tables
suffice for this issue's tests) and the de novo π̂ of M3-11 through an argument · creates
`src/giap/mapping/placement_fdr.py` · pure Python · estimate ≈ 330 code lines · **statistical**:
`docs/methods/placement-fdr.md` (sections "Strata and ranking", "π̂_wrong", "Same-strain input",
"Calibration"), signed off by the statistician (§14.1); the statistician also owns attack 5
(§8.5), run in M3-21

## Spec sections implemented

- §7.2.5 Procedure 1–6, quoted in full:
  1. "Each peptide counts once, at its best cost, in the target search (m_T) and in the decoy
     search (m_D). The strata are m ∈ {0, 1, 2}."
  2. Ranking score E_i (M3-09).
  3. "Within stratum m, take the target list (peptides with m_T = m) and the decoy list (peptides
     with m_D = m), both ordered by E ascending. At each threshold e:
     FDR_m(e) = π̂_wrong (D_m(e) + 1) / T_m(e). The placement q-value of peptide i is
     q_i = min over e ≥ E_i of FDR_m(e). Accept when q ≤ `stats.alpha.placement` (0.01)."
  4. "π̂_wrong: m = 0, database input: the mean PEP of the input peptides if PEP is present;
     otherwise the maximum input q-value (≤ `input.q_value_max`); otherwise 1. m = 0, de novo
     input: 1 − the lower 95 % bound of the de novo precision for the peptide's score bin; 1 if
     that cannot be estimated. m ≥ 1: always 1 (conservative)."
  5. "**Same-strain input.** Peptides from a search against this genome's own proteins
     (`search_db = same`) skip step 3 for exact hits (`placement_q = 0`, flag `same_db`). Their
     only ambiguity is multi-locus placement. Tolerant hits still go through steps 1–4."
  6. "`mapping/mapping_fdr.tsv` reports one row per length class × m: T, D, π̂_wrong, the stratum
     FDR π̂(D + 1)/T and n_accepted. Per-m totals are reported too."
- G3 (strata by mismatch count, ranking by E pools length classes; per length × m table still
  reported). Decision 13 (α = 0.01). §8.1 rows "Peptide placed at a locus (exact)" and "… with ≥ 1
  difference". Appendix B.6.
- §9.3 `mapping_fdr.tsv`: `length_class`, `n_diff`, `targets`, `decoys`, `pi_wrong`,
  `fdr_stratum`, `n_accepted`, `decoy_type` enum{reverse, shuffle}; "Per-m total rows carry
  `length_class = all`."

## Public interface

```python
@dataclass(frozen=True, slots=True)
class PeptideCosts:          # one row per mapped query peptide (built from the hit tables)
    peptide: str; length: int; source: Literal["database", "denovo"]; search_db: Literal["same", "other", "auto"]
    m_target: int | None; m_decoy: int | None; e_target: float | None; e_decoy: float | None
    pep: float | None; q_value: float | None; pi_denovo: float | None   # from M3-11 (de novo only)

def best_costs(hits: pd.DataFrame, peptides: pd.DataFrame, freqs: ResidueFreqs) -> pd.DataFrame: ...
def pi_wrong_database(peptides: pd.DataFrame, q_value_max: float) -> float: ...
def stratum_qvalues(e_target: np.ndarray, e_decoy: np.ndarray, pi_target: np.ndarray) -> np.ndarray: ...
def placement(costs: pd.DataFrame, cfg: GiapConfig) -> tuple[pd.DataFrame, pd.DataFrame]: ...
#  → (per-peptide: peptide, m_target, placement_q, accepted, flags) , mapping_fdr rows (section = "placement")
```

- `stratum_qvalues`: `T(e) = #{targets with E ≤ e}`, `D(e) = #{decoys with E ≤ e}`; the FDR is
  evaluated at every target E; `q_i = min_{e ≥ E_i} FDR(e)`, then clipped to ≤ 1. Ties in E
  between targets and decoys: all entries with E ≤ e count (a decoy tied with a target is
  counted at that threshold). With constant π̂ this equals M4-01's
  `tda_qvalues(-e_target, -e_decoy, pi, plus_one=True)` (§8.3); reuse it if merged, and keep a
  test asserting the equality.
- **Stratum lists are per `decoy_type`** (database peptides vs reverse, de novo vs shuffle).
- **π̂ per target.** Database, m = 0: one value for all (`pi_wrong_database`: mean PEP over the
  mapped evidence peptides with PEP; else max q ≤ `q_value_max`; else 1). De novo, m = 0:
  `pi_denovo` of the peptide's score bin. m ≥ 1: 1. With per-peptide π̂ (de novo), FDR(e) uses the
  mean π̂ over targets with E ≤ e (reduces to the spec formula when π̂ is constant; see Q2).
- **same_db.** Peptides with `search_db = same` and `m_target = 0` get `placement_q = 0`, flag
  `same_db`, and are left out of the stratum-0 target **and** decoy lists.
- **mapping_fdr rows**: per `decoy_type` × length class (`mapping.length_classes`) × m ∈ {0,1,2}:
  `targets` = T, `decoys` = D (whole class, no threshold), `pi_wrong` (mean π̂ over the class's
  targets), `fdr_stratum` = π̂(D + 1)/T (NA when T = 0), `n_accepted` (q ≤ α); plus one
  `length_class = all` row per m and decoy type. Rows are sorted by decoy type, m, length class.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | B.6: stratum m = 1, π̂ = 1, 150 targets with E below the first decoy | FDR = (0 + 1)/150 = `0.00667`; those 150 have q = `0.00667`, accepted at 0.01 |
| 2 | B.6: 50 targets, no decoys, π̂ = 1 | FDR = 1/50 = `0.02`; none accepted |
| 3 | stratum 1, π̂ = 1, targets E = 0.001, 0.002, 0.003, decoy E = 0.0025 | FDR at the target thresholds `1, 0.5, 0.667`; q = `0.5, 0.5, 0.667` |
| 4 | same lists in stratum 0 with π̂ = 0.02 | q = `0.01, 0.01, 0.01333`; the first two accepted (q ≤ 0.01 inclusive) |
| 5 | database input, PEP of the mapped peptides `[0.001, 0.003, 0.008]` | π̂ = `0.004` |
| 6 | database input, no PEP, q-values `[0.001, 0.0098, 0.03]`, `q_value_max = 0.01` | π̂ = `0.0098` |
| 7 | database input, neither PEP nor q | π̂ = `1` |
| 8 | `search_db = same`, exact target hit | `placement_q = 0`, flags `same_db`, not counted in stratum-0 T or D; its 1-difference twin with `search_db = same` goes through steps 1–4 |
| 9 | peptide with `m_target = None` (no target hit) and `m_decoy = 1` | contributes to D of stratum 1 only; no `placement_q` (goes to `unmapped`, reason `no_hit`) |
| 10 | one 8-aa and one 12-aa target in stratum 1 | one `length_class` row each (`8`, `12-14`) and one `all` row with `targets = 2` |

## Edge cases

- A stratum with T = 0: no q-values; the `mapping_fdr` row has `targets = 0`, `fdr_stratum = NA`.
- A peptide accepted in stratum m is accepted at all its target loci at cost m (multi-locus
  ambiguity is uniqueness, M3-14, not FDR).
- Input decoy rows (`is_decoy`) are never mapped (M3-17); they do not enter these lists.
- Order independence: shuffling input rows gives identical outputs.

## Acceptance criteria

- [ ] `tests/unit/test_placement_fdr.py`: examples 1–10.
- [ ] Property: q-values are monotone non-decreasing in E within a stratum; adding a decoy never
  lowers any q; q ∈ [0, 1].
- [ ] Calibration, statistics-only (`tests/calibration/test_placement_fdr_lists.py`, marker
  `calibration`, R = 50 fast / 300 full from `tolerances.yaml`): correct peptides get target E
  from a low-E distribution and no decoy hit; wrong peptides get exchangeable target/decoy E
  (independent draws from one distribution), with known π_wrong; mean FDP over R replicates
  ≤ q + 2 SE(FDP) at q ∈ {0.01, 0.05, 0.10} (§8.4 "FDR"). The genome-level calibration with
  simulated peptides is M3-21.
- [ ] Methods note `docs/methods/placement-fdr.md` exists (sections listed above, with the B.6
  example and the calibration design) and the statistician has signed it off before merge.

## Questions for the lead

1. Which peptides are searched against the decoys? **Proposed default:** every mapped query
   peptide (exact decoy search for all, tolerant decoy search for all ≥ 10 aa), so D estimates the
   chance hits of the whole input; the target tolerant search runs only for exact-target
   leftovers. With π̂ = 1 for m ≥ 1 this is conservative.
2. De novo π̂ varies per score bin: how does it enter FDR_0(e)? **Proposed default:** the mean π̂
   of targets with E ≤ e (equals the spec formula when π̂ is constant). Alternative: run stratum 0
   separately per score bin (smaller lists, more conservative with the +1).
3. "Mean PEP of the input peptides": over which peptides? **Proposed default:** the aggregated,
   non-decoy peptides that enter mapping (q ≤ `q_value_max` or q missing).
4. `length_class` labels in TSV. **Proposed default:** `7`, `8`, `9`, `10-11`, `12-14`, `15-19`,
   `20-40` (single value when lo = hi), plus `all`.

## Lead decision after statistician review (pre-implementation): ADR-0010 supersedes Q2
Per-peptide π̂: FDR_m(e) = (W_m(e) + c) / T_m(e), W = Σ π̂_j over decoy hits with E ≤ e (own π̂;
missing → 1; m ≥ 1 → 1), c = max π̂ over searched peptides of that decoy type. Interface
`stratum_qvalues(e_target, e_decoy, pi_decoy, pi_plus)`. `mapping_fdr`: `fdr_stratum` = (W + c)/T,
`pi_wrong` = (W + c)/(D + 1). Where this section and earlier bullets (the "mean π̂ over targets" text,
`pi_wrong` = mean over the class's targets) disagree, this section wins. Test changes as listed by the
statistician: `test_per_target_pi_uses_mean_over_targets_up_to_e` → targets E = 1, 2, 3, one decoy E =
2.5 with weight 0.5, c = 0.6 → q = [0.3, 0.3, 0.366667]; `test_denovo_stratum0_uses_per_peptide_pi` →
second target π̂ 0.6: q = [0.3, 0.3, 0.366667], fdr_stratum 0.366667, pi_wrong 0.55;
`test_constant_pi_equals_m4_tda_qvalues` → decoy weights `np.full(ed.size, pi)`, c = pi; property tests
take per-decoy π̂ and scalar c ≥ max weight, permute weights with decoys, and add "raising any weight
or c never lowers a q".

## Lead decisions after statistician review (PR #41, BLOCK)
1. The same-strain shortcut applies to database rows only: `same_exact = is_db & (search_db == "same") &
   (m_t == 0)`. De novo rows never get it, whatever `search_db` the input reader filled in.
2. `pi_wrong_database` is computed over database rows that stay in the stratum-0 lists:
   `df[is_db & ~same_exact]` (same peptides as c). Single-`search_db` inputs are unchanged.
   The statistician's attack tests/calibration/attacks/test_attack_placement_same_db.py must pass.
3. `low_complexity` is added to the §9.4 placement flag vocabulary by ADR-0011 (accepted by Khoa 2026-10-14; it changes
   an output vocabulary).
4. Unresolved `search_db = auto` is treated as `other` (conservative); M3-17 resolves `auto` before M3-10.
