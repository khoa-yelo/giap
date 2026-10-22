# M9-26: Scale and determinism: a 5-Mb render, a K-12 report, byte-identical images, `-resume`

Milestone M9 · **pass 2** · depends on: M9-13 · estimate ≈ 140 code lines (+ ≈ 200 test lines) ·
not statistical

## Spec sections implemented

- §10.1 "Large genomes. Plots aggregate into bins when a track would draw more than 20,000
  elements." §14.4 (speed loop after M3/M8/M10; the visualisation loop data "rendered from the
  simulated and E. coli runs"), §11.2 benchmarks, §12 (benchmarks against `baseline.json` on one
  node class), R7/TE7, §6.3 `-resume`, §10.3 size target < 25 MB.

## Work

1. **5-Mb render.** `benchmarks/test_plots_bench.py`: the `plots` and `report` steps on a 5-Mb
   giapsim run (≈ 4,500 genes, ≈ 50,000 placements): wall time and peak RSS recorded in
   `benchmarks/baseline.json`; every track ≤ 20,000 elements (M9-24 item).
2. **A 20-Mb stress fixture** (synthetic `RunView` frames only, not a pipeline run): 25,000 genes on
   one contig → the gene tracks of V1 and V2's points bin; a test asserts the binned rows.
3. **K-12 report.** One `sbatch scripts/slurm/giap_run.sbatch` on the K-12 full-DB Bakta recording
   (`tests/data/recorded/bakta/1.12.2-full/k12`) with M3-22 in-silico peptides and a simulated
   abundance column (the M8-06 hardening run): `report.html` size recorded (< 25 MB), the plots go
   into the M9-25 approval packet. Job ID in PROGRESS.md.
4. **Byte-identical images.** A test renders every plot twice on the snapshot fixture and compares
   `.svg`, `.pdf`, `.png` bytes (M9-01 item 7); plus `report.html` byte identity across two
   `giap report` calls.
5. **`-resume`.** An nf-test (stub) shows `PLOTS` and `REPORT_HTML` are cached on `-resume` when no
   table changed, and rerun when one table changed.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| render time > 10 min on K-12 | profile; the speed loop (§14.4) applies; record before/after in an ADR if `baseline.json` changes |
| PNG bytes differ across matplotlib patch versions | expected; identity is asserted within one environment only |
| report > 25 MB on K-12 | warning event; reduce embedded rows or locus count via `report.max_locus_plots` (config, no default change) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | 25,000 genes on 20 Mb | V1 gene tracks `binned = true`, ≤ 20,000 rows each |
| 2 | two renders of V07 | identical bytes in all formats |
| 3 | `-resume` with no change | `PLOTS` cached |
| 4 | K-12 run | `report.html` < 25 MB; size in PROGRESS.md |

## Acceptance criteria

- [ ] Benchmarks recorded; examples 1–3 as tests; example 4 as Slurm evidence.

## Questions for the lead

1. Should the K-12 run here double as M8-06's "meaningful real-data run"? **Proposed default:**
   yes, one job serves both.
