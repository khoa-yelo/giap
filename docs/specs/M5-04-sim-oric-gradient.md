# M5-04: Simulator: oriC feature, ori–ter gradient, `features` contract table (C15), null scenarios

Milestone M5 · track F (simulator, tooling) · depends on: none (extends merged M4-07 layout, M4-08
expression) · estimate ≈ 50 code lines in `src/` (contract entry) + ≈ 200 simulator lines in
`tests/sim/` (+ ≈ 200 test lines) · **not statistical** (tooling; the statistician reads the scenario
definitions) · **can start now** · no tools, no Slurm

## Spec sections implemented

- §11.3 "Non-coding features: … an oriC-like region next to a marker gene"; "TU-level log abundance = μ +
  replicon effect + region effects + planted islands (± effect over 5–20 TUs) + an AR(1) background
  (ρ tunable) + noise".
- §7.3.5(e) "oriC comes from Bakta; the terminus is taken as oriC + L/2. Replichore position = circular
  distance from oriC / (L/2), in [0, 1]."
- §8.4 rows "Autocorrelation | gene-level noise (label null); TU-structured noise without inter-TU
  structure (TU null)" and "Islands | AR(1) backgrounds without islands"; §8.5 attack 1 "Autocorrelated
  abundance with no islands → islands, autocorrelation".
- §6.2 `parse_bakta` "→ `features.parquet` (rRNA, tRNA, CRISPR, oriC, gaps)" — merged in M2-04
  (`giap.steps.annotation._FEATURE_COLS`) but not part of the M4 input contract (M4-06 C1–C14).

## Public interface

**Contract (src, implementer):** `giap.genetable.contract.TABLES["features"]` = `feature_id` str, `type`
str (enum of `giap.model.genes.NcType`), `contig` str, `start` int, `end` int, `strand` str
(`+`, `-`, `.`), `gene` str?, `product` str?, `amino_acid` str?, `anticodon` str? — the columns
`parse_bakta` writes. It is **optional** in `validate_inputs` (absent file = no oriC → M5-08 `no_oric`);
rules: `contig` known, 0 ≤ start < end ≤ 2·length. Called **C15** in the M5 contract.

**Simulator (tests/sim/giapsim, implementer):**
```python
# layout.py
LayoutConfig.oric_frac: float | None = None   # oriC midpoint on chr1 as a fraction of its length; None = no oriC
Layout.features: pd.DataFrame                 # C15 columns; one "oric" row (232 bp, strand ".") or empty
# expression.py
ExpressionConfig.oriter_log2: float = 0.0     # ori minus ter difference of TU-level log2 abundance
ExpressionConfig.oriter_shape: Literal["cosine", "linear"] = "cosine"
# genes frame gains replichore_true {right, left, NA} and replichore_position_true (float, NA)
# contract_tables.py: write_contract_tables also writes features.parquet (empty table when no oriC)
# scenarios.py (new)
def label_null(n_genes: int) -> tuple[LayoutConfig, ExpressionConfig]: ...
def tu_null(n_genes: int, detection: str = "abundance") -> tuple[LayoutConfig, ExpressionConfig]: ...
def ar1_background(n_genes: int, rho: float) -> tuple[LayoutConfig, ExpressionConfig]: ...
def island_power(n_genes: int, rho: float, n_tus: int, effect_sd: float) -> tuple[LayoutConfig, ExpressionConfig]: ...
def oriter_gradient(n_genes: int, rho: float, oriter_log2: float) -> tuple[LayoutConfig, ExpressionConfig]: ...
```
- **oriC placement:** the 232-bp feature is centred in the intergenic gap of chr1 whose midpoint is
  nearest to `oric_frac × length` (never overlapping a gene); `feature_id` `<chr1 prefix>_oric1`.
- **Replichore truth:** for each gene, midpoint m (mod length); clockwise distance d = (m − o) mod L from
  the oriC midpoint o; `right` with position d/(L/2) when d ≤ L/2, else `left` with (L − d)/(L/2).
  Plasmids and layouts without oriC: NA.
- **Gradient:** TU-level shift g(p) at the TU midpoint's replichore position p: cosine
  (`oriter_log2`/2)·cos(π p); linear `oriter_log2`·(0.5 − p). Added to `tu_log2` before islands and the
  AR(1); requires an oriC (`ValueError` otherwise). `oriter_log2 = 0` leaves every existing output
  byte-identical (no extra random draws).
- **Scenarios** (one place for the calibration and power generators; all use `circular` chr1,
  `detection = "covariates_only"` unless stated):
  - `label_null`: `tu_sd = 0`, `position_decay_log2 = 0`, `stoich_sd` and `gene_sd` iid → no structure
    at all (both nulls hold);
  - `tu_null`: `ar1_rho = 0`, `tu_sd = 2`, decay 0.2 → within-TU correlation only (TU-block null holds,
    gene-label null false); `detection = "abundance"` so detection inherits the TU structure;
  - `ar1_background(ρ)`: default expression at the given ρ, no islands, no region effects;
  - `island_power`: `ar1_background` plus one `IslandSpec` at a random-but-seeded TU start away from
    regions; effect sign gives the direction;
  - `oriter_gradient`: `ar1_background` with `oric_frac = 0.8` and the cosine gradient.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `LayoutConfig(oric_frac=None)` (default) | `features` empty with C15 columns; every pre-M5 output byte-identical to before |
| 2 | `oric_frac = 0.5` on a 1,000-gene chr1 | one `oric` row, length 232, no gene overlap, midpoint within one inter-TU gap of 0.5·L |
| 3 | oriC midpoint o = 900 on L = 1,000; gene midpoints 100, 400, 500, 650 | (right, 0.4), (right, 1.0), (left, 0.8), (left, 0.5) |
| 4 | gene crossing the origin, start 990, end 1,010 (unwrapped), o = 0 | midpoint 0 → (right, 0.0) |
| 5 | `oriter_log2 = 2`, cosine, `tu_sd = 0`, `gene_sd = stoich_sd = 0`, decay 0 | TU at p = 0 has μ + 1, at p = 1 μ − 1, at p = 0.5 μ |
| 6 | `oriter_log2 = 1` without oriC | `ValueError` |
| 7 | `write_contract_tables` | `features.parquet` present; `validate_inputs` reports no issue; without the file, no issue either (optional) |
| 8 | `tu_null(1000)` | lag-1 correlation of gene `abundance_log2_true` within TUs > 0.5; between the last gene of a TU and the first of the next ≈ 0 (\|r\| < 0.1 over 20 seeds pooled) |
| 9 | `label_null(1000)` | Moran's I lag 1 of `abundance_log2_true` ≈ 0 (\|I\| < 0.08 for 95 % of 100 seeds) |

## Edge cases

- `oric_frac` outside [0, 1) → `ValueError`; `oric_frac` with no chr1 → `ValueError`.
- The `origin_crossing` layout option and an oriC near 0 together: the oriC is still placed in a gap,
  the crossing gene keeps `end > length`.
- `layout_from_simulation(sim)` (M0 sequence simulator, which already writes an oriC-like block into
  `bakta_like.json`) fills `Layout.features` from the JSON features of type `oriC` (and `crispr`).

## Acceptance criteria

- [ ] `tests/sim/test_layout_oric.py`, `tests/sim/test_expression_oriter.py`,
  `tests/sim/test_scenarios.py` (test-author): examples 1–9.
- [ ] `tests/unit/genetable/test_contract.py` gains C15 cases (column list equals giapsim's literal copy;
  optional absence; out-of-range rows rejected).
- [ ] Determinism: same seed → byte-identical `features.parquet` and expression frame.

## Questions for the lead

1. C15 is an addition to the M4 input contract read by M5 only. **Proposed default:** add it to
   `TABLES` as optional (no ADR: it documents an existing M2 output, nothing changes for M4).
2. Gradient applied at TU level (shared by the TU's genes) rather than per gene? **Proposed default:**
   TU level, consistent with the region and island effects of M4-08.

## Lead decisions (2026-10-20, after the statistician's attack check)

11. Q1: C15 `features` is optional in the contract, no ADR. Bakta writes oriC strand `?` (K-12,
    Nissle); the contract test checks it arrives as `.` through `parse_bakta`.
12. Q2: the ori–ter gradient applies at TU level.
13. Add scenarios: `label_null` + gradient (`tu_sd = 0`, cosine 2 log2); `tu_null` with
    `short_genes_in = "prophage"` and `detection = covariates_only`; `ar1_background(1000, 0.9)`;
    `draft_ar1` (500-gene chromosome + 133 linear 30-gene contigs); `island_power(..., detection =
    "abundance")`.
