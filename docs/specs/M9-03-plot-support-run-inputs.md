# M9-03: Plot support files from run inputs (`plots/support/`): landmarks and GC skew

Milestone M9 · **pass 1** · depends on: M9-02 · estimate ≈ 160 code lines (+ ≈ 200 test lines) ·
not statistical · writes **ADR-M9-A**

## The gap this closes

V1 draws GC skew and oriC/terminus marks, and V2 draws oriC and terminus lines. None of these is in
a §9.1 file: GC skew needs the genome sequence, and the oriC call is in the intermediate
`features.parquet`. `giap report OUT` must rebuild the plots from `OUT/` alone (§7.9.4). So the
`plots` step, which receives the genome, features and replicons as inputs during `giap run`, writes
two small support tables once. `giap report` reads them back and never needs the genome.

## Spec sections implemented

- §10.2 V1 ("GC skew"), V2 ("oriC and terminus lines").
- §7.3.5(e) ("oriC comes from Bakta; the terminus is taken as oriC + L/2"), as implemented by
  M5-08 (`giap.context.oriter.oric_position`, `ORIC_SPREAD_FRAC`, the `no_oric` rule).
- §7.9.4 (rebuild without recomputing statistics), §9.5 (descriptor), §6.5 (coordinates), §6.6.

## ADR-M9-A (written by this issue)

"Plot support files under `OUT/plots/support/`". A support file holds data that a plot needs, that
no §9.1 file holds, and that only a step with the run's inputs (or a fitted model, M9-17) can
compute. Rules: TSV per §9.2 (`.tsv.gz` with gzip `mtime=0` when larger than 1 MB, so bytes are
deterministic); one Frictionless resource each in `schemas/datapackage.json` with the custom
property `role: plot_support`; written by the producing step during `giap run`; read, never
written, by `giap report`; not a result table (no §9.2 test columns, no claims). Non-breaking
addition to §9.1. Status: lead-level.

## Public interface

`src/giap/viz/support.py` (pure)
```python
def replicon_landmarks(replicons: pd.DataFrame, features: pd.DataFrame | None) -> pd.DataFrame
def gc_skew(genome: Mapping[str, str], replicons: pd.DataFrame, window_bp: int) -> pd.DataFrame
```
`src/giap/io/plot_support.py`: `write_support(name, df, plots_dir)`; schemas
`src/giap/schemas/plot_support/{replicon_landmarks,gc_skew}.json`. The `plots` step (M9-02) calls
both before rendering when `genome`, `features` and `replicons` are among its inputs.

## Tables

**`plots/support/replicon_landmarks.tsv`**: `contig`, `length_bp`, `topology` enum{circular,
linear}, `oric_position` (1-based; NA), `ter_position` (1-based; NA), `status` enum{called,
no_oric, linear}. One row per replicon in FASTA order.

- `oric_position` = M5-08's `oric_position` on that contig (the same function, so V1/V2 and the
  ori–ter rows can never disagree), written 1-based (`+ 1`).
- `ter_position` = ((oriC₀ + L/2) mod L) + 1, with oriC₀ the 0-based position and integer
  division `L // 2`.
- Linear contigs and contigs without an oriC call: both NA, `status` `linear` or `no_oric`.

**`plots/support/gc_skew.tsv`**: `contig`, `start`, `end` (1-based closed window), `gc_skew`
((G − C)/(G + C), NA when G + C = 0), `cumulative_gc_skew` (running sum along the contig from
position 1). Windows of `window_bp` tiled from position 1; the last window is shorter; windows do
not wrap.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| `giap report` on `OUT/` without `plots/support/` | V1 omits the GC-skew track and V2/V1 the landmark marks; each notes "not available" |
| several oriC calls within `ORIC_SPREAD_FRAC` · L | their circular mean (M5-08 rule) |
| several oriC calls farther apart | `no_oric` |
| contig shorter than `window_bp` | one window covering it |
| a window of only N | `gc_skew` NA; cumulative carries the previous value |
| lowercase sequence | counted case-insensitively |
| `genome` input absent (a hand-built step call) | no `gc_skew.tsv`; landmarks still written |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | contig `c1`, L = 100,000, circular, one oriC feature 0-based [99,900, 100,100) (crosses origin; midpoint 0) | `oric_position` 1, `ter_position` 50,001, `status` called |
| 2 | contig `c2`, L = 8,000, circular, no oriC feature | `oric_position` NA, `ter_position` NA, `status` no_oric |
| 3 | contig `lin`, linear | both NA, `status` linear |
| 4 | sequence `GGGC` + 9,996 × `A`, `window_bp` 10,000 | one window 1–10,000, `gc_skew` (3 − 1)/4 = 0.5 |
| 5 | 25,000 bp with windows 10,000 | three windows: 1–10,000, 10,001–20,000, 20,001–25,000 |
| 6 | windows with skew 0.5, −0.25, NA | cumulative 0.5, 0.25, 0.25 |
| 7 | `local/runs/val-mvp-2026-10-22/OUT` replicons (chr1 100,000; plasmid_1 8,000) with sim100k features | two rows, in that order |
| 8 | the same inputs twice | byte-identical files |

## Acceptance criteria

- [ ] `tests/unit/viz/test_support.py`: examples 1–6; `oric_position` agrees with
  `giap.context.oriter` on every giapsim seed in a small property test.
- [ ] `tests/unit/io/test_plot_support.py`: example 8; gzip for > 1 MB with `mtime=0`; schema
  validation with frictionless.
- [ ] `tests/unit/steps/test_plots_step.py` (extends M9-02): support written when `genome` is
  supplied, absent otherwise; `giap report` never writes `plots/support/`.
- [ ] ADR-M9-A committed (`docs/adr/ADR-00NN-plot-support-files.md`).

## Questions for the lead

1. GC-skew window: §10.2 gives none. **Proposed default:** 10 kb, the V1 peptide-density bin, as a
   module constant (not a config key; adding a key would touch the protected `defaults.yaml`).
2. Should `giap report` recompute GC skew when `--config` names a genome? **Proposed default:** no;
   §5.1 gives `giap report` no genome argument.

## Lead decisions (2026-10-22)
1. GC-skew window 10 kb (the V1 peptide-density bin), a module constant, not a config key.
2. `giap report` never recomputes GC skew; with no support file the GC-skew track takes its empty state.
