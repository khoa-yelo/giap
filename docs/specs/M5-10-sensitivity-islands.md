# M5-10: Sensitivity: the `islands` family (Jaccard of island gene sets under `min` and `exclude`)

Milestone M5 · track R · depends on: M5-06, M5-09; **M4-19** (sensitivity framework, not merged yet) ·
estimate ≈ 100 code lines (+ ≈ 150 test lines) · **statistical (light)**: label `stats`, statistician
reviews the flagged-set definition · **off the gate path** (like M4-19) · no tools, no Slurm

## Spec sections implemented

- §8.6 quoted: "`basic` (default) reruns the context statistics … under `abundance.missing ∈ {min,
  exclude}`." "Permutations. Sensitivity reruns use B = 199, since they only compare flagged sets."
  "Families: islands (gene sets), regions (q ≤ 0.05 IDs), …"
- §2.4 G9 (`tables/sensitivity.tsv`, Jaccard stability).

## Public interface

`src/giap/context/sensitivity.py` (M4-19) gains
```python
def islands_family(genes: pd.DataFrame, cfg: GiapConfig, rng: np.random.Generator,
                   *, tus: pd.DataFrame, replicons: pd.DataFrame) -> frozenset[str]: ...
FAMILIES["islands"] = ...   # bound to the context inputs by run_sensitivity (Q1)
```
- Flagged set = { `"<direction>:<locus_tag>"` for every gene of every island of M5-06 } (an island is
  reported only when its FWER p ≤ `stats.alpha.context`). Direction is part of the ID so a gene moving from
  a high to a low island counts as a change.
- Reruns use `stats.bootstrap = 199` and `stats.permutations = 199` (block length) for the primary and
  each variant (M4-19 rule: the primary is recomputed at the same B).
- RNG: ~~the child stream of M4-19 keyed `islands/<variant>`~~ — replaced by lead decision 1: one family stream `for_step("sensitivity/islands", seed)`.
- Rows: `family = islands`, `variant`, `n_primary`, `n_variant`, `n_shared`, `jaccard` (NA when both sets
  are empty).

## Examples

| # | Input | Expected row |
|---|---|---|
| 1 | primary {high:G1, high:G2, high:G3}, `missing=min` {high:G2, high:G3, high:G4} | n_primary 3, n_variant 3, n_shared 2, jaccard 0.5 |
| 2 | no island in either | jaccard NA, counts 0 |
| 3 | primary {low:G7}, variant {high:G7} | n_shared 0, jaccard 0.0 |
| 4 | `output.sensitivity = off` | no `islands` rows |
| 5 | `missing = exclude` | detected genes without abundance drop out of the scores (M5-06 step 1); rows still produced |

## Edge cases

- `exclude` when every detected gene has an abundance → variant equals primary (jaccard 1 or NA).
- The flagged set uses the best-window gene list (M5-06 Q3); merging differences show up as changed sets.

## Acceptance criteria

- [ ] `tests/unit/context/test_sensitivity_islands.py`: examples 1–5.
- [ ] Integration: `run_sensitivity` on the M5 fixture writes `islands` rows into a schema-valid
  `sensitivity.tsv`; deterministic for a fixed seed.

## Methods note

`docs/methods/sensitivity.md` §"Islands": the flagged-set definition, B = 199 for both primary and
variants.

## Questions for the lead

1. M4-19's `Family` signature takes only the gene table, but islands also need the TUs and replicon
   topology. **Proposed default:** `run_sensitivity` binds the extra context inputs with
   `functools.partial` when it calls a family that declares them (keyword-only parameters); no change for
   `regions`.
2. Should autocorrelation also be a family? §8.6 lists only islands for M5. **Proposed default:** no (R1).

## Lead decisions (2026-10-20, after the statistician's attack check)

38. Q1 (`functools.partial`) and Q2 (no autocorrelation family): accept.

## Lead decisions (2026-10-22, after the test author)
1. The test author's resolutions are accepted:
   - **one stream per family:** `for_step("sensitivity/islands", seed)`, as M4-19 decision 2 and the merged step use; this replaces the spec's `islands/<variant>`;
   - **B = 199** for both `stats.permutations` and `stats.bootstrap`;
   - **every island:** the flagged set uses every gene of every `islands_frame` row;
   - **`tus` is required;**
   - **M5-06 is reached** through `giap.context.islands`.
2. **Existing pins.** The lead filtered three existing sensitivity tests to `family == "regions"` and widened the FAMILIES pin. The implementer:
   - adds `tus` to the registry pin and to context.nf's `stepIn('sensitivity', ...)`;
   - regenerates the golden `sensitivity.tsv` files with GIAP_UPDATE_GOLDEN=1, noting it in the PR.
