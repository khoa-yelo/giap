# M4-21: Conformance with the real M2/M3 outputs; re-record goldens through the real path

Milestone M4 · track T · depends on: M4-20, **and M2 (`parse_bakta`, `operon_pairs`/TUs,
`context_layers`/regions) and M3 (`loci`) merged** · estimate ≈ 150 code lines (adapter fixes, if any)
(+ ≈ 150 test lines) · **cannot start before M2 and M3 merge**

## Spec sections implemented

- §6.1: "The only join between the mapping track and the annotation track is the Bakta locus tag."
- §6.4: "Both executors call the same step functions, and an integration test checks that they
  produce byte-identical tables on the synthetic genome." (M4 part: the local executor path through
  M2 → M3 → M4 on the simulated genome with `annotation.precomputed`.)
- §11.4 golden fixture; §11.6 TE7.

## Work

1. Run M2's `parse_bakta`, TU and region steps and M3's mapping steps on the 100-kb synthetic genome
   (precomputed layers, local executor) and validate every file against `giap.genetable.contract`
   (M4-06): unmark `tests/contract/test_m4_input_contract.py::test_real_outputs` (`pending_sample`).
2. Fix mismatches on whichever side deviates from the agreed contract (contract changes go through
   one PR reviewed by the M2/M3 owners).
3. Run M4 end to end on those real outputs; frictionless-validate all M4 tables.
4. Add `tests/golden/sim50k_real/` produced through the real path; keep `sim50k/` (fixture path) for
   M4-only regressions.
5. One gene-FDR calibration replicate through real M3 mapping (sanity, not gated): FDP within the
   fast-suite bound.

## Examples

| # | Scenario | Expected |
|---|---|---|
| 1 | M2 `genes.parquet` for the 100-kb genome | `validate_table("genes", …)` returns no errors |
| 2 | M2 TUs | every gene in exactly one TU (M4-06 Q1 honoured) |
| 3 | M3 `peptide_hits.parquet` from in-silico tryptic peptides of the simulated proteome | assigning classes have `gene`; M4-11 accepts ≥ 99.9 % of planted peptides |
| 4 | M4 on real outputs, run twice | byte-identical tables (TE7) |

## Acceptance criteria

- [ ] Contract test green on real outputs; `pending_sample` marker removed.
- [ ] `tests/golden/test_m4_golden.py` covers both `sim50k` and `sim50k_real`.
- [ ] `PROGRESS.md` records the conformance run and any contract changes (with PR links).

## Methods note

None.

## Questions for the lead

1. If M2's region membership or TU rules differ from the simulator's (e.g. TU singletons), which side
   adapts? **Proposed default:** the contract (M4-06) is the arbiter; the producer adapts unless the
   lead decides the contract was wrong.

## Lead decisions (2026-10-21)

1. Q1: the M4-06 contract (`giap.genetable.contract`) is the arbiter. The producer adapts unless the
   lead decides the contract was wrong; such a case is reported to the lead, never fixed by changing
   the contract in this issue.

## Lead decisions after the test-author step (2026-10-21)

2. **M2 adapts (C7/C8).** With `annotation.precomputed.regions` set, `region_layer` still writes the
   C7 `regions` and C8 `region_genes` Parquet, built from the precomputed TSV, so `context_stats` is
   planned on precomputed runs. The contract is unchanged.
3. **M4 consumers read the keyed layout.** `run_context_stats` and `run_sensitivity` resolve inputs
   through `work_inputs`, as `run_gene_table` already does, so the real `OUT/work` layout works.
4. **`context_stats` must be planned under `giap run --executor local`.** Wire its `peptide_hits` /
   `protein_hits` / `unmapped` inputs to the producer `loci` actually emits (the `mapping` outputs).
   This is integration work and outranks everything else in this issue.
5. The real-path golden case is **`sim100k_real`** (the 50-kb genome plants no regions); keep `sim50k`
   for M4-only regressions. Record it with `GIAP_UPDATE_GOLDEN=1` once 2–4 pass.
