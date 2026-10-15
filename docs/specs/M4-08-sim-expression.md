# M4-08: giapsim expression and detection truth

Milestone M4 · track F (simulator) · depends on: M4-07 · estimate ≈ 300 lines in `tests/sim/giapsim/`
(tooling) · **can start now**

## Spec sections implemented

- §11.3 Expression: "TU-level log abundance = μ + replicon effect (plasmid copy number) + region
  effects (for example, prophage −3) + planted islands (± effect over 5–20 TUs) + an AR(1) background
  (ρ tunable) + noise." "Gene level = TU level + position-in-TU decay + a stoichiometry factor +
  noise." "peptide detectability logistic in peptide features and log abundance".
- §8.4 generators: "Detectability | logistic truth on physicochemical covariates"; "Region tests | AR(1)
  abundance, detection from p_detect only"; "Islands | AR(1) backgrounds without islands".
- §8.5 attack nulls 2 ("Detectability confounded with region membership (for example, short genes
  concentrated in prophages)") and 4 ("High-copy plasmids").
- §7.3.6 power: "regions with detection shifts of 10–50 %".

## Interface (`tests/sim/giapsim/expression.py`)

```python
@dataclass(frozen=True)
class ExpressionConfig:
    mu_log2: float = 20.0
    tu_sd: float = 2.0; gene_sd: float = 0.5; stoich_sd: float = 0.7
    ar1_rho: float = 0.5                       # along TU order within each replicon (circular-stationary start)
    position_decay_log2: float = 0.2           # per position in the TU
    plasmid_copy_log2: Mapping[str, float] = {}        # e.g. {"p1": log2(10)}
    region_effect_log2: Mapping[str, float] = {}       # region_type → abundance shift, e.g. {"prophage": -3}
    region_detect_shift: Mapping[str, float] = {}      # region_type → additive shift on the detection logit
    islands: tuple[IslandSpec, ...] = ()               # (replicon, first TU index, n_tus, effect in SD)
    detection: Literal["covariates_only", "abundance"] = "covariates_only"
    beta: Mapping[str, float] = DEFAULT_BETA   # intercept, log_length, log1p_n_obs, gravy, log2_abundance

@dataclass(frozen=True)
class Expression:
    genes: pd.DataFrame   # locus_tag, tu_log2, abundance_log2_true, p_detect_true, detected_true,
                          # n_obs_true, gravy_true, in_island, region_effect_log2
    islands: pd.DataFrame # truth for M5

def simulate_expression(layout: Layout, cfg: ExpressionConfig, rng: np.random.Generator) -> Expression: ...
def simple_covariates(aa: str) -> tuple[int, float]: ...   # (n_obs tryptic 7–30, GRAVY), own minimal code
```
- `p_detect_true = sigmoid(β0 + β1·log(length_aa) + β2·log1p(n_obs) + β3·GRAVY [+ β4·abundance_log2
  when detection = "abundance"] + region_detect_shift)`. `detected_true ~ Bernoulli(p_detect_true)`.
  Covariates are computed by `simple_covariates`, independent of `giap` (M4-10 recomputes them).
- Partial genes keep values (they are excluded later by `giap`, not by the simulator).
- **Full-simulator hook:** M3-19 gives the nucleotide simulator a plain log-normal
  `gene_log_abundance` and says M4 replaces it "behind the same function". This issue provides
  `expression_for_simulation(sim_genes, tus, cfg, rng)` returning the same per-gene log abundance, and
  swaps it in if M3-19 has merged (else M4-20 does).
- AR(1): x_0 ~ N(0, σ²/(1−ρ²)), x_t = ρ x_{t−1} + ε_t over TUs in replicon order.

## Examples

| # | Call | Expected |
|---|---|---|
| 1 | default config on the default 1,000-gene layout, seed 1 | `p_detect_true` ∈ (0, 1); mean detected share within 0.55 ± 0.15 with `DEFAULT_BETA` |
| 2 | `ar1_rho = 0.9` | lag-1 autocorrelation of TU-level `tu_log2` (before noise) within 0.9 ± 0.05 on a 4,500-gene layout |
| 3 | `region_effect_log2 = {"prophage": -3}` | mean `abundance_log2_true` of prophage genes − mean elsewhere within −3 ± 0.5 |
| 4 | `plasmid_copy_log2 = {"p1": log2(10)}` | plasmid genes shifted by +3.32 ± 0.5 |
| 5 | `region_detect_shift = {"bgc": logit(0.7) − logit(0.5)}` | for genes with p ≈ 0.5 the BGC detection rate rises toward 0.7 (power input, 20 % shift) |
| 6 | `detection = "covariates_only"` | `detected_true` independent of `abundance_log2_true` given covariates (logistic fit coefficient on abundance ≈ 0, |z| < 3) |
| 7 | `simple_covariates("MAGLLKAAGPEPTIDERPAAAGGSTLKAAAAAAAK")` | n_obs = 3 (`AAGPEPTIDER`, `PAAAGGSTLK`, `AAAAAAAK`; `MAGLLK` too short), trypsin/P rule |
| 8 | same layout, config, seed twice | identical frames |

## Edge cases

- Empty `islands` → `islands` truth frame has zero rows and the columns.
- `ar1_rho` must be in [0, 1); 0 → iid TU effects.
- A replicon with one TU: AR(1) reduces to its stationary draw.

## Acceptance criteria

- [ ] `tests/sim/test_expression.py`: examples 1–8 (seeded; tolerances as stated).
- [ ] No import of `giap` (grep test shared with M4-07).
- [ ] 4,500 genes in < 1 s.

## Methods note

The generators are described in each consumer's methods note (`detectability.md`, `region-tests.md`);
the statistician reviews this issue because the nulls it produces gate M4.

## Questions for the lead

1. Default β for detection (`DEFAULT_BETA`): log length +1.2, log1p n_obs +0.8, GRAVY −0.6, and
   the intercept tuned once so the default layout gives ≈ 55 % detected (≈ −8.9 for 280-aa median
   genes with ≈ 13 observable peptides), an E. coli-like single-run share? **Proposed default:** yes;
   the statistician may tune before the calibration tests are written.

## Lead decisions (after the test-author's pass)
1. The test-author's resolutions are accepted: `tu_sd` = AR(1) innovation σ; `tu_log2` = μ + replicon +
   region + island + AR(1), one value per TU, no extra TU noise; planted shifts/islands are deterministic
   (paired same-seed runs differ by exactly the shift); island first-TU index 0-based in replicon TU order,
   effect = `effect × tu_sd` at ρ = 0; `n_obs` per M4-10 (distinct tryptic peptides 7–30 aa); natural logs;
   ρ ∉ [0, 1) or unknown detection mode → ValueError; `expression_for_simulation(genes, tus, cfg, rng)`.
2. Examples 3 and 4 (±0.5): one seed's scatter (SD ≈ 0.8) exceeds ±0.5, so the spec's tolerance is tested
   on the mean over 30 seeds, and the exact shift on paired same-seed runs.
3. The calibration file runs in suite `regions`; its detection checks in `detectability`.
