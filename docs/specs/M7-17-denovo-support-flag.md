# M7-17: Flag unsupported de novo-only novel-ORF calls (ADR-0050 option 3)

Milestone M7 · track S (statistics) / P (pipeline) · **pass 2** · depends on: M7-11b (#267, the
`novel_orfs` step wiring of the class FDR and the `novel_class_fdr` side table; must merge first),
M7-09 (`apply_evidence`, `rbs_supported`), M3-17 (`steps.mapping.resolve_sources`); the
calibration part also needs M7-13b (the G4-k3 harness `tests/helpers/m7_gate_novel_records.py`,
branch `m7/gate-b`) · blocks: M7-13c (the gate record lists the ADR-0050 cells with this flag as
the shipped mitigation) · estimate ≈ 85 counted code lines (+ ≈ 280 test lines, not gated) ·
**stats** label (the statistician reviews the rule and the methods-note paragraph)

## Spec sections implemented

- ADR-0050 option 3, accepted by Khoa on 2026-10-23: "Keep decision 3 and declare a documented
  limitation: exact systematic de novo errors are outside what the class FDR controls. … De novo-only
  calls with no database or RBS/TIR support are flagged." (The ADR file on `m7/gate-b` still reads
  "proposed"; the lead updates it.)
- §9.4 "Flag vocabulary (controlled; extend only through the schema)", Novel ORFs row: one new
  flag, added through `schemas/tables/novel_orfs.schema.json`, as M7-06/07/09 did.
- §7.6.3 "**RBS support** means the OSTIR initiation score at the chosen start is ≥ the
  `novel_orfs.rbs_min_quantile` (default 0.10) quantile of scores at annotated starts of detected
  genes." This is the only RBS/TIR evidence GIAP computes for novel ORFs: OSTIR's `expression` is
  the TIR (§7.1.6), so "RBS support" and "TIR support" are one test, `rbs_supported`.
- §7.6.4 and M7-11 decision 3 (de novo ORF q = minimum peptide subgroup q): **unchanged**.
- §7.6.6 "Calibration: class FDR on simulations with planted novel ORFs and de novo errors";
  M7-13 rule 8 (G4-k3).

## What exists today (origin/m7/class-fdr-wiring-b, 8d3d3ae8)

- `orf_regimes` (`novel/fdr.py`) reads each ORF's peptides' aggregated C9 `source` (I/L folded,
  non-decoy rows); `picked_decoy` needs makedb decoys and a non-`denovo` source; else
  `denovo_subgroup` with any `denovo` source; else `unestimated`. A source of `auto` counts as
  non-`denovo` there.
- `steps.mapping.resolve_sources(peptides, cfg)` resolves `auto` table-wide: `input.peptide_source`,
  else `database` when the table has any q-value or PEP, else `denovo`. Mapping uses this.
- `apply_evidence` sets `rbs_supported = ostir_tir ≥ rbs_threshold(gene_table,
  cfg.novel_orfs.rbs_min_quantile)` (NA when either is NA; threshold `None` below `RBS_MIN_GENES`
  = 20 detected genes with a TIR) and the flag `rbs_supported` when true. The class rule requires
  `rbs_supported` for `intergenic` and `nterm_extension` (`RBS_CLASSES`).
- `FLAG_ORDER` ends with `readthrough_candidate`; the schema `flags` pattern lists the same 12.
- The step writes `novel_class_fdr.parquet` (`_SIDE_COLS`, 10 columns).

## The flag

**Name:** `denovo_only_unsupported` (§9.4 style: lower snake case, states the evidence, as
`class_fdr_unestimated` and `rbs_supported` do).

**Rule** (per `novel_orfs` row, every row, called or not):

```
denovo_only_unsupported  :=  has_denovo  and  not db_support  and  not rbs_exempt
```

1. **Peptide source.** `src = resolve_sources(peptides, cfg)` on the aggregated C9 `peptides`
   frame the step already reads (M7-11b step 1). Rows with `is_decoy` or `is_contaminant` true
   are ignored. Per ORF, its `Orf.peptides` are matched to C9 rows by the I/L-folded peptide
   (`_il`, as `orf_regimes` does); a folded key with several C9 rows takes the union of their
   sources.
2. **`db_support`** = at least one of the ORF's peptides has resolved source `database`.
   - This is "≥ 1 database PSM": the C9 aggregation already gives `database` priority over
     `denovo` for a peptide seen by both (`io.peptides._PRIORITY`, `mixed_source` warning), and
     every ORF peptide already passed `input.q_value_max` (default 0.01) as a mapping query.
     No new q threshold.
   - It counts whatever the database was (makedb six-frame, Bakta proteins, `search_db` same or
     other) and whatever the placement (exact or tolerant).
   - Only the ORF's own peptides count. For `nterm_extension` the extended gene's in-frame
     peptides are `cds_inframe` and never ORF peptides (§7.6.1), so they do not count.
3. **`has_denovo`** = at least one of the ORF's peptides has resolved source `denovo`. Together with
   rule 2: every ORF peptide that is found in C9 is de novo.
4. **`rbs_exempt`** = `rbs_supported is True` **and** the class is not in `RBS_CLASSES`
   (`intergenic`, `nterm_extension`). For those two classes RBS support is part of the class rule,
   so every called ORF there has it and it is not support beyond the call; for them the flag
   reduces to `has_denovo and not db_support` (Q1).
   - `rbs_supported` NA (no OSTIR layer, no start codon, `ostir_missing`, or
     `rbs_threshold_unestimated`) is **no** support.
   - Thresholds: `novel_orfs.rbs_min_quantile` (0.10) and `RBS_MIN_GENES` (20), both existing. No
     new threshold or config key.
5. **Not support:** `promoter_supported`, `terminator_supported`, `junction_spanned`, the number
   of peptides, and the class q itself (ADR-0050 names only database and RBS/TIR support; Q9).
6. **Regime-independent.** The rule reads sources and `rbs_supported`, not `regime`.
   - `picked_decoy` ORFs always have `db_support`, so never carry the flag.
   - `denovo_subgroup` ORFs can have `db_support` when the run has no makedb decoys (a database
     peptide then does not make the ORF picked); they are then not flagged.
   - `unestimated` ORFs (no decoys, ≥ 2-peptide rule) of de novo peptides are flagged too: exact
     systematic errors defeat that rule as well (Q3).
7. **Without the `peptides` input** (M7-11b step 1: pass-1 behaviour) the flag is never set and
   the side-table cell is null: the sources are unknown.

**Nothing else changes.** `class_fdr`, `q`, `called`, `rank` and the class rules are computed
exactly as in M7-11b; the flag is applied after them and read by nothing in GIAP's call path.
So G3 power, every G4 number and the ADR-0050 strict xfails (G4-k3 cds_outframe and antisense at
0.01 and 0.05, intergenic and GC at 0.01) are unaffected; a fix to the method still shows as
XPASS.

## Public interface

### `src/giap/novel/class_call.py` (pure; additions)

```python
def db_supported(orfs: Sequence[Orf], peptides: pd.DataFrame,
                 sources: pd.Series) -> dict[str, tuple[bool, bool]]: ...
    # orf_id -> (has_denovo, db_support); sources = resolve_sources(peptides, cfg),
    # index-aligned with peptides; decoy and contaminant rows ignored (rules 1-3)
```

### `src/giap/novel/evidence.py`

- `apply_evidence(..., denovo_support: Mapping[str, tuple[bool, bool]] | None = None)`. When it
  is `None` (library callers, all existing tests) nothing changes. Otherwise, after
  `rbs_supported` is set, the row gets `denovo_only_unsupported` by rule 4.
- `denovo_only_unsupported(cls: str, has_denovo: bool, db_support: bool,
  rbs: bool | None) -> bool` — the rule, one line, exported for the unit tests.

### `src/giap/novel/outputs.py` and the schema

- `FLAG_ORDER` gains `denovo_only_unsupported` **last** (after `readthrough_candidate`). If the
  pending `class_mismatch` (M7-11b Q3) lands, it goes after `class_fdr_unestimated` as planned;
  the two do not conflict.
- `novel_orfs.schema.json` `flags` pattern gains the token in both alternations; the `flags`
  description gains "; `denovo_only_unsupported`: every peptide is de novo, none is a database
  PSM, and there is no RBS support beyond the class rule (ADR-0050)".
- `orf_table` keeps a `denovo_only_unsupported` it receives (it only rewrites
  `class_fdr_unestimated`), so no change there.

### `novel_orfs` step (`steps/novel.py`)

1. `_class_fdr` already reads `peptides`; it also returns `db_supported(orfs, peptides,
   resolve_sources(peptides, cfg))` (or `None` without `peptides`).
2. `apply_evidence(..., denovo_support=…)`.
3. `_SIDE_COLS` gains `("denovo_only_unsupported", _BOOL)` last: true/false, null without
   `peptides`. The side table stays in `novel_orfs.tsv` row order.
4. One INFO event `novel_orfs.denovo_support` per run with `peptides` (none without):
   `data = {n_orfs, n_denovo_only, n_flagged, n_called, n_flagged_called}`, where `n_denovo_only`
   counts `has_denovo and not db_support`. Message: `"novel_orfs: <n_flagged_called> of
   <n_called> calls are de novo-only without database or RBS support (ADR-0050: exact de novo
   errors are outside the class FDR)"`. Not a WARNING and not a run-flag token (Q4).

### Report (`report/build.py`, `_novel_fdr` card)

When some called row carries the flag, the novel-ORF card's estimate gains the suffix
`"; <k> of <m> calls de novo-only without database or RBS support (exact de novo errors are not
covered by the class FDR)"`. The same suffix follows `unestimated`. With no flagged call the card
is byte-identical to today (Q5).

### Not changed

`novel_orfs.gff3` (attributes are `rank`, `n_peptides` only), `novel_orfs.faa`,
`annotation_patch.gff3`, the manifest, `events` other than the new one, Nextflow wiring (no new
input or output).

## Examples

α = 0.01, RBS threshold `rbs_min` = 400 (20+ detected genes), I/L folded peptides, every ORF
meets its class rule unless stated. "db" = resolved source `database`, "dn" = `denovo`. Flags are
listed in `FLAG_ORDER`; other flags of the row are elided as "…".

| # | Input | Expected |
|---|---|---|
| 1 | `intergenic`, 2 dn peptides, `denovo_subgroup` q 0.004, `ostir_tir` 500 | `rbs_supported` true; flags `rbs_supported;…;denovo_only_unsupported`; `called` true, q and rank as M7-11b (RBS is the class rule, not extra support) |
| 2 | `cds_outframe`, 3 dn, q 0.002, `ostir_tir` 500 | `rbs_supported` true, **no** flag (`rbs_exempt`); called true |
| 3 | as 2 with `ostir_tir` 120 | `rbs_supported` false; flag set; called true |
| 4 | as 2 with no start codon (`ostir_tir` NA) | `rbs_supported` NA; flags `…;no_start_codon;…;denovo_only_unsupported`; called true |
| 5 | `antisense`, peptides `AEPLK` (db) and `GGTVSR` (dn); run without makedb decoys → `denovo_subgroup` | `db_support` true → no flag |
| 6 | makedb run, `intergenic` ORF with 2 db peptides, `picked_decoy` q 0.006 | no flag; side cell false |
| 7 | `cds_outframe`, 3 dn, q 0.03 (> α), `ostir_tir` 50 | flag set although `called` false (G4 reads q ≤ 0.05) |
| 8 | source-less table (`input.peptide_source: auto`, no q-value or PEP column → resolved `denovo`); `sorf`, 2 peptides, no decoys → `unestimated` | flag set (regime-independent); called true by the ≥ 2-peptide rule; `class_fdr_unestimated;…;denovo_only_unsupported` |
| 9 | peptide `ELVLLK` has C9 rows from database and de novo (`mixed_source`; aggregated `database`); it is one of the 2 peptides of a `cds_outframe` ORF; the other is dn | no flag |
| 10 | ORF peptide `LLGK`, C9 row `ILGK` source db | I/L folded match → `db_support`, no flag |
| 11 | `nterm_extension` of `b0002`, 2 dn upstream peptides, `ostir_tir` 700, gene `b0002` has 9 db in-frame peptides | flag set: the gene's peptides are not ORF peptides and RBS is the class rule |
| 12 | step without `peptides` (the M7-09 fixture) | no flag; side column null; no `novel_orfs.denovo_support` event; `novel_orfs.tsv`, GFF3s, FAA byte-identical to M7-11b |
| 13 | run of examples 1, 2, 3, 6, 7 | event `{n_orfs: 5, n_denovo_only: 4, n_flagged: 3, n_called: 4, n_flagged_called: 2}`; card `"class FDR ≤ 0.006; 2 of 4 calls de novo-only without database or RBS support (exact de novo errors are not covered by the class FDR)"` |

## Edge cases

| Situation | Behaviour |
|---|---|
| ORF with mixed db and dn peptides | `db_support` → never flagged, in any regime (examples 5, 9) |
| absent `rbs_novel` (OSTIR) layer | `ostir_tir`, `rbs_supported` NA → no exemption → every de novo-only ORF flagged (the RBS classes are not called anyway); the existing `layer_missing:rbs_novel` WARNING is unchanged |
| fewer than 20 detected genes with a TIR | threshold `None`, `rbs_supported` NA, `rbs_threshold_unestimated` as today → no exemption |
| absent `peptides` input | flag never set, side cell null, no event (rule 7, example 12) |
| annotated-gene search database (Bakta proteins, no makedb) | database peptides rarely reach novel loci; when one does it is `db_support`. ORFs are `unestimated` (`decoys_not_makedb`) or `denovo_subgroup`; the flag follows the sources only |
| `nterm_extension` (annotated-gene class) | only the upstream peptides count (example 11); RBS is never an exemption |
| `pseudogene`, `sorf` | RBS can exempt (not in `RBS_CLASSES`) |
| `cterm_extension` | flag computed as for others; never called before M7-16 |
| ORF peptide absent from C9 (should not happen) | contributes no source; an ORF with no matched peptide has `has_denovo` false → no flag |
| decoy or contaminant C9 row with the same folded peptide | ignored (decoys never support; contaminant peptides never reach a novel ORF, §7.6.1) |
| circular replicon, origin-crossing ORF | flag independent of coordinates; same rule |
| linear contig, open-3′ ORF | same rule; `partial` and the flag can co-occur |
| minus strand | same rule (no coordinate is read) |
| no candidate ORF | no row; event `{n_orfs: 0, …: 0}` when `peptides` present; card unchanged |
| ties in `rank_key` | unchanged: the flag is not a rank key and does not reorder rows |
| `class_mismatch` ORF | the flag reads the ORF's own peptides and its TSV `class` (the class rule's class) |
| resume duplicates the event | one event per step run; resume reruns the step and writes it once |

## Line estimate (counted: `src`, `workflow`, `scripts`, added + deleted)

| Part | Lines |
|---|---|
| `novel/class_call.py` `db_supported` | 22 |
| `novel/evidence.py` rule + `apply_evidence` parameter | 12 |
| `novel/outputs.py` `FLAG_ORDER`; `novel_orfs.schema.json` pattern and description | 4 |
| `steps/novel.py` sources, side column, event | 25 |
| `report/build.py` card suffix | 12 |
| `docs/methods/novel-orfs.md` (not counted) | — |
| **Total** | **≈ 75–85** |

Tests ≈ 280, not gated: unit ≈ 140, step ≈ 70, report ≈ 20, calibration ≈ 50.

## Test plan and acceptance criteria

1. **Unit** `tests/unit/novel/test_denovo_support_flag.py`:
   - `denovo_only_unsupported` truth table: 7 classes × `has_denovo` × `db_support` ×
     `rbs ∈ {True, False, None}`; in particular `intergenic`/`nterm_extension` with `rbs` true
     are flagged and `cds_outframe` with `rbs` true is not.
   - `db_supported`: examples 5, 8, 9, 10 and the decoy, contaminant, absent-peptide and I/L edge
     rows; `auto` resolution both ways (q-value present → database; none → denovo).
   - `apply_evidence` with `denovo_support`: examples 1–4, 7, 11; `called`, `rank`, `class_fdr` and
     every other column equal to the run without `denovo_support` (the "nothing else changes"
     check, on the same frames); flags ordered by `FLAG_ORDER`.
2. **Step** `tests/unit/steps/test_novel_denovo_support_step.py`: examples 6, 12, 13 through
   `novel_orfs`; side-table column name, dtype, position and null without `peptides`; the event's
   exact `data` and message; no WARNING is added.
3. **Report**: the card of example 13 and the unchanged card with no flagged call.
4. **Schema**: a `novel_orfs.tsv` carrying the flag validates with `frictionless validate`; a
   misspelt token does not.
5. **Existing tests pass unchanged**, in particular `test_outputs.py`, `test_novel_evidence.py`,
   `test_class_call.py`, `tests/unit/steps/test_novel_*`, the report snapshots and
   `test_novel_local_vs_nextflow.py`; the ADR-0050 strict xfails stay strict xfails (they are not
   touched and must not XPASS).
6. **Calibration** (`tests/calibration/test_novel_orfs_gate.py` gains a G4-k3 companion, or a new
   `test_denovo_support_flag_gate.py` reusing `inject_false_orfs` and `k3_record`; stream
   `novel_orfs.g4k3.flag`; fast R = `GATE_R["fast"]`, full R = 300; it reuses the G4-k3 runs, so
   no extra simulation in the full tier):
   - *Asserted, every replicate:* for each six-frame class, every injected false ORF (an
     `inject_false_orfs` key, all of which carry ≥ 2 exact products) whose side row has
     `q ≤ 0.05` and `class_rule_met` carries `denovo_only_unsupported`, **except** a
     `cds_outframe` or `antisense` row with `rbs_supported` true (rule 4). For `intergenic` this
     is every such row. A violation is an `AssertionError`; zero injected rows reaching q ≤ 0.05
     in a class over the block is a `RuntimeError` (vacuous), as G4-k3's reach check.
   - *Recorded* (diagnostics jsonl and PROGRESS.md), per class: the number of such injected rows,
     the share flagged, and the **escape share** (rows exempted by chance RBS support).
   - *Recorded*, per regime (`picked_decoy` from the G3/G4 database arm, `denovo_subgroup` from
     the G4 de novo arm) and class: the share of **true** ORFs (truth matches, rule 5 of M7-13)
     with `q ≤ 0.01` and the class rule met that carry the flag. Expected by construction: 0 in
     the picked arm, 1 for `intergenic`/`nterm_extension` in the de novo arm, and 1 − (chance
     RBS rate) for `cds_outframe`/`antisense` there (giapsim plants RBS only at intergenic and
     N-terminal true starts).
7. **Methods note** `docs/methods/novel-orfs.md` §"Class FDR" gains §"Exact systematic de novo
   errors (ADR-0050)" with the paragraph below and the recorded shares; the statistician signs off.
8. **Backlog / gate:** M7-13c's gate record names the flag next to the ADR-0050 cells; a backlog
   line records Q6 and Q8.
9. Lint, format, mypy clean; `check_imports.py` clean.

## Methods-note paragraph (to be added verbatim, numbers filled from criterion 6)

> **Exact systematic de novo errors (ADR-0050).** The class FDR treats wrong peptides as decoy
> hits treat them: a wrong read lands on the genome by chance, so a reversed or shuffled genome
> decoy is as likely to receive it. A de novo read that is wrong but **exactly equals** the
> translation of a non-coding six-frame segment breaks that assumption. It places with zero
> cost, like a true peptide, and no decoy can match it except by chance, so no decoy-based
> estimate counts it. When several such reads fall in one segment, the ORF also gets several
> chances at a small q (decision 3). In M7-13 G4-k3 (20 injected false ORFs per class, 3–5 wrong
> reads each) the ORF-level FDP was 0.075 for `cds_outframe` and 0.072 for `antisense` at both
> q = 0.01 and 0.05, and 0.017 for `intergenic` at 0.01; the exact reads alone account for all of
> it. GIAP therefore does not claim class-FDR control against exact systematic de novo errors.
> q values and calls are unchanged. Calls that this limitation can affect carry
> `denovo_only_unsupported`: every peptide is de novo, none is a database PSM, and there is no RBS
> support beyond what the class rule already requires. RBS support at the 0.10 quantile is weak
> evidence: in the simulator a non-coding start passes it with probability <escape share>, so the
> flag is a pointer to the evidence gap, not a second error estimate. Users who need such calls
> should confirm them by a database search against the `giap makedb` six-frame database, or by
> independent evidence.

## Questions for the lead (each with the reversible default in use)

1. **RBS exemption vs the requested calibration assertion (the main gap).** ADR-0050 says "no
   database or RBS/TIR support". Read literally, two things go wrong:
   - `intergenic` and `nterm_extension` calls always have RBS support (class rule), so they
     would never be flagged, yet `intergenic` is one of the failing G4-k3 cells;
   - in giapsim a non-coding start's TIR is lognormal(6, 0.8) against a threshold near the 10 %
     quantile of gene TIRs (mostly lognormal(8, 0.8)); with several candidate starts and the best
     one chosen, an injected false `cds_outframe`/`antisense` ORF passes RBS by chance often
     (unmeasured; plausibly 10–40 %). Those ORFs escape the flag, so "every exact-arm ORF called
     at q ≤ 0.05 carries the flag" **cannot hold** for those classes under any RBS exemption.
   **Default:** rule 4 (RBS exempts only outside `RBS_CLASSES`); the calibration asserts the rule
   on real pipeline output (all of `intergenic`; the others except chance-RBS rows) and records
   the escape share. **Alternative (A):** drop the RBS exemption entirely (flag = de novo-only
   without database support). Then the assertion holds in every class, the flag is simpler, but
   it departs from Khoa's wording, so it would go to Khoa. One line changes either way.
2. **Flag on uncalled rows.** **Default:** yes, every row, so G4-style readers at q ≤ 0.05 (> α)
   see it; the event and card count called rows separately.
3. **`unestimated`-regime de novo ORFs.** **Default:** flagged (regime-independent, rule 6).
   Alternative: only `denovo_subgroup`, which is narrower than "de novo-only calls".
4. **WARNING or token.** **Default:** INFO event only. A WARNING would fire on nearly every de novo
   run and make `--strict` exit 3; a run-flag token beyond §9.4 needs Khoa (ADR-0041).
5. **Report text.** **Default:** the card suffix above, shown only when a called row is flagged.
   If M9 renders a V13 caveat panel later, it reads the same count.
6. **Flagged share of true ORFs in a mixed-input arm.** In the two existing arms the share is
   fixed by construction (0 picked, ≈ 1 de novo). An informative number needs database + de novo
   input in one run. **Default:** record the two arms only; backlog a mixed arm for M7-13c or
   MVP-04.
7. **Flag name.** **Default:** `denovo_only_unsupported`. Alternatives: `denovo_only`,
   `unsupported_denovo`.
8. **`auto` source inconsistency** (found while speccing): `orf_regimes` treats an unresolved
   `auto` C9 source as non-`denovo`, while mapping resolves it table-wide. A source-less de novo
   table thus gives `unestimated` ORFs instead of `denovo_subgroup`. **Default:** the flag uses
   `resolve_sources` (example 8); `orf_regimes` is left alone here (changing it changes regimes and
   q, which this issue must not do) and a backlog line asks M7-11 follow-up to resolve sources
   the same way.
9. **Other evidence as support** (promoter, terminator, homolog). **Default:** no; ADR-0050 names
   database and RBS/TIR only, and promoter/terminator are labels that never gate (ADR-0012).
10. **ADR status.** ADR-0050 still says "proposed" and its interim default says the flag "is not
    added until Khoa decides". **Default:** the lead marks it accepted (option 3) before this
    spec's test-author starts.

## Lead decisions (spec-writer 2e7f3fd4)
1. **Q1:** the default stands. The RBS exemption applies, as Khoa's option 3 wording says, but not for intergenic or nterm_extension. The
   calibration test exempts injected cds_outframe/antisense rows that pass RBS by chance, and records
   that escape share per class. If the share exceeds 0.2 in either class, the lead takes "drop the
   RBS exemption" to Khoa with the numbers.
2. Q2–Q9: the defaults are accepted. The mixed-arm and `orf_regimes` auto-source fixes go to the backlog.
3. Q10: the lead marks ADR-0050 accepted (option 3 with the flag, option 1 later) before the test-author starts.
