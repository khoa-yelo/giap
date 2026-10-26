# OSTIR full recordings re-recorded without the user site (2026-10-26)

Agent-Role: validation. Follow-up to #390 (`PYTHONNOUSERSITE=1` for tool runs).

## What changed

The full OSTIR 1.1.3 recordings for K-12 MG1655 and Nissle 1917
(`$GIAP_ROOT/recordings/ostir/{k12_mg1655,nissle1917}/`, recorded 2026-10-15) were folded with
ViennaRNA 2.7.1 imported from `~/.local`. They were not folded with the locked env's 2.6.4,
although `run.log` printed 2.6.4. Both samples were re-recorded with the unchanged recorder
`scripts/slurm/record_tools.sbatch ostir`, using `PYTHONNOUSERSITE=1` and the `giap-ostir` env
(lock hashes in `versions.txt` unchanged; `RNA.__version__` 2.6.4 from the env's site-packages).

- Slurm job **1896079** (`giap-record-ostir`, 8 CPUs, rye-02, 19 min, COMPLETED, 0 failures).
  The job wrote to a staging directory (`RECORD_OUT=$GIAP_ROOT/staging/ostir-rerecord-2026-10-26`)
  because the recorder `rm -rf`s an incomplete target. Changed files were then moved into the
  recording directories.
- Every replaced file is kept next to its replacement as `<file>.usersite-2026-10-26.bak`
  (`ostir.csv`, `command.txt`, `runtime.txt`, `SHA256SUMS`). `ostir_input.csv`, `versions.txt`,
  `run.log` and `.giap-complete` are byte-identical and were left in place. The new `SHA256SUMS`
  verifies (`sha256sum -c`).
- `command.txt` and `runtime.txt` differ only in `-j 8` (was `-j 9`) and in wall time and node.

| sample | input rows | output rows | changed rows | old `ostir.csv` sha256 | new `ostir.csv` sha256 |
|---|---|---|---|---|---|
| k12_mg1655 | 4,372 | 4,371 | 3 | `67893813…374c` | `6c307fe8…7b10` |
| nissle1917 | 4,689 | 4,688 | 3 | `265b17d5…fc8b` | `d634a3d7…c68e` |

The inputs are byte-identical to the old ones: sha256 `897f7382…e3f` (K-12) and `96ed3f49…40b5` (Nissle).

## Evidence: the 6 changed rows

All 6 rows are ATG starts at `start_position` 51. Only the folding-dependent columns change
(`dG_mRNA`, `dG_start_codon` and `start_codon` are equal). Old = user-site ViennaRNA 2.7.1,
new = env ViennaRNA 2.6.4.

| gene | product | expression old → new | RBS_distance_bp | dG_total | dG_rRNA:mRNA | dG_spacing | dG_standby |
|---|---|---|---|---|---|---|---|
| KOIAIB_00226 | yafE | 0.4262 → 0.1148 | 4 → 10 | 20.3289 → 23.6087 | -1.181 → -0.881 | 0.0039 → 2.9837 | 0.0 → 0.0 |
| KOIAIB_01380 | sapC | 3.3507 → 194.2349 | 2 → 5 | 15.1741 → 5.025 | -9.581 → -3.281 | 16.4491 → 0.0 | 0.0 → 0.0 |
| KOIAIB_02406 | rhmR | 0.0386 → 0.0291 | 0 → -1 | 26.3339 → 27.0346 | -11.081 → -10.981 | 17.2089 → 17.2096 | -4.3 → -4.9 |
| JMJJOG_00751 | FCD domain protein | 43.3505 → 2833.4105 | 2 → 5 | 8.7741 → -1.675 | -11.781 → -5.781 | 16.4491 → 0.0 | 0.0 → 0.0 |
| JMJJOG_01590 | rhmR | 0.0386 → 0.0291 | 0 → -1 | 26.3339 → 27.0346 | -11.081 → -10.981 | 17.2089 → 17.2096 | -4.3 → -4.9 |
| JMJJOG_02638 | sapC | 3.3507 → 194.2349 | 2 → 5 | 15.1741 → 5.025 | -9.581 → -3.281 | 16.4491 → 0.0 | 0.0 → 0.0 |

`rhmR` and `sapC` have the same 5′ window in both strains, which is why they change the same way.
Three of the six rows change TIR (`expression`) by 58× or more, so the error is not a small
numerical drift.

Cross-checks:
- The new `ostir.csv` files are byte-identical to the lead investigation's isolated rerun
  (session scratchpad `ostir_full/{k12_mg1655,nissle1917}`). The investigation's user-site rerun
  (`*_user`) is byte-identical to the old recordings, so the user site explains the whole
  difference.
- The committed excerpts in `tests/data/recorded/ostir/1.1.3/{k12_mg1655,nissle1917}/ostir.csv`
  contain none of the 6 genes, and every excerpt row is still present verbatim in the new full files.
- `tests/contract/test_ostir_contract.py`: 14 passed against the new full recordings. These are
  the full-recording parse and input-writer tests, plus the committed-excerpt tests.

## Consumers

| consumer | uses these recordings? | affected? | action |
|---|---|---|---|
| `tests/contract/test_ostir_contract.py` (full-run tests via `$GIAP_ROOT`) | yes | no: parse/input checks only | reran: 14 passed |
| `tests/data/recorded/ostir/1.1.3/{k12_mg1655,nissle1917}` excerpts | yes (copied rows) | no: byte-identical | none |
| `tests/data/recorded/ostir/PROVENANCE.md` | yes (full-file sha256s) | yes: the output hashes were stale | updated in this PR |
| `tests/data/recorded/ostir/1.1.3/sim100k`, `sim100k_novel` | separate recordings | no: byte-identical on rerun (investigation) | none |
| `local/runs/pilot/pass1/{out_a1,out_a2,out_a2_0667e366,out_a2_1549a930}`, `pass2/out_a3`, `pass3/out_a2` (Nissle) | no: OSTIR is run by the pipeline | **yes**: all three Nissle genes carry the old user-site values (pass3 as JNOIEE_00751/01591/02641) | rerun needed for `rbs` / `determinants` (`tir_ostir`) outputs; full pipeline run, not started |
| `local/runs/mvp04/20261023_{a1,a2}` (K-12) | no: pipeline OSTIR | **yes**: all three K-12 genes carry the old values (as KOIAIB_00227/01381/02407) | rerun needed if the TIR/determinants outputs are used again; not started |
| `validation/pilot/README.md`, `report_template.md`, `metrics/` | no | no TIR/RBS metric reported (proteomics metrics only); determinants outputs are not quoted | none |
| `docs/validation/mvp-run-2026-10-22.md`, `local/runs/val-mvp-*`, `m806b-sim100k`, `demo-quickstart` | no: simulated genomes | none of the 6 windows occur in their inputs; any other usersite effect on sim genes is unverified | none (later runs pick up #390) |
| calibration (`tests/calibration/test_sim_m7_*`) | no: simulated, pipeline OSTIR | not from these recordings | none |

The pipeline runs were matched to the recordings by 5′-window sequence in their published
`ostir_input.csv`. Their OSTIR tasks predate #390 (merged 2026-10-26 01:15 PDT). pilot pass3's
OSTIR task was launched at 01:03 and its `versions.yml` has no `viennarna:` line. Only 3 of about
4,700 genes per genome change, so the pilot's gene-level determinant results should move very
little. The exact effect is unmeasured until those runs are repeated under #390.

No code, thresholds or holdout data were touched. This report contains no row-level genome or
pilot data beyond the 6 rows above.
