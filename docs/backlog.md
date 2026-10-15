# Backlog

Issues may live on GitHub; this file mirrors the milestone plan (GIAP_SPEC.md §13).

- [ ] M0 scaffold (in progress)
- [ ] M1 inputs and preflight
- [ ] M2 annotation wrappers, full Bakta DB, recorded fixtures
- [ ] M3 mapping and decoy FDR; pilot kit
- [ ] M4 gene table, detectability, region tests
- [ ] M5 autocorrelation and islands
- [ ] M6 operons; RegulonDB calibration (source currently unreachable)
- [ ] M7 variants and novel ORFs (declare GIAP's VCF FILTER values with descriptions; M1-05 writes generic ones)
- [ ] M8 function and determinants
- [ ] M9 visualisation and report
- [ ] M10 pilot validation
- [ ] M11 release 1.0

- M3: `mapping.tolerant.swaps` and `mapping.tolerant.mass_equivalents` config keys are not consumed by M3-07; decide which issue implements them or remove them (from M3-07 review, PR #48).
- M6: the simulator plants stop-codon readthrough (C-terminal extensions, ADR-0013) for cterm_extension calibration (from M3-19 test-author note).
- M2-09: antiSMASH databases with `url: "image"` (databases bundled in the image) — setup edge case from M2-08 review.
- M2-09 (from M2-08 review, low): delete the geNOMAD tarball after unpacking (§3.9); reject non-200 responses before hashing; cap redirects; test `_provenance` requires prev status ok; drop redundant tar `filter` weakening check.
- M2-09 follow-up: GO release directories can announce a different `data-version` date in the OBO header (2026-06-19 → 2026-06-15, 2026-08-05 → 2026-07-26), so `verify_go` rejects them and the pin stays at 2026-05-19; verify the header against a pinned `data_version` instead of the directory date.
- M4: migrate existing calibration tests to `_criteria.calibrate()` so ADR-0015 two-block confirmation and sharding apply to all suites (M4-05 decision 7; test-author).
- M4 (calibration migration follow-up): suites marked only `calibration` (stats_*, placement lists/genome) have no `calibration_full` variant, so the nightly full array skips them; adding one needs a new R — statistician decision.
- M4 (statistician, #69): gate the sim_proteomics q-value FDR check and the one-sided de novo coverage checks through `calibrate`.
- M4 (statistician, #69 minor): the 1.5q teeth test at q = 0.01 has power < 1 at fast R (seed-fragile) — raise R or use 2q at q = 0.01; teeth tests rerun in every shard task — run them once.
- M4-20: simulator `multi_identical_pairs` should give paired genes identical `aa` in the layout (M4-09 decision 2).
- M4-14: record the declared abundance type (e.g. maxlfq) in the input metadata; C10 has no scale column (M4-09 decision 3).
- M4-09: add a benchmarks/ entry for simulate_evidence on 4,500 genes (budget 3 s; unit test keeps a 30 s guard).
- M4-20 (with M4-09 decision 2): for shared_frac, the simulator's second hit sequence is not in that gene's aa (M4-09 review).
