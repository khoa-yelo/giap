# M7-11b: Novel-ORF class FDR wiring: `class_fdr`, the class-q call rule and the q-floor warning

Milestone M7 · track S (statistics) / P (pipeline) · **pass 2** · depends on: M7-11 (#237, the
pure `giap.novel.fdr`; must merge first), M7-09 (the wired `novel_orfs` step and `apply_evidence`),
M3-10 (`placement`, π̂ per de novo peptide), M3-14 (`LocusClassifier`), M3-23 (makedb IDs) ·
blocks: M7-13b (G3, G4 and G4-k3 in docs/specs/M7-13-calibration-gate.md) · estimate ≈ 330 counted
code lines (+ ≈ 450 test lines, not gated) · **stats** label (statistician reviews the call rule
and the floor rule, not the arithmetic, which M7-11 already owns)

## Spec sections implemented

- §7.6.4, all three branches: "A call needs class q ≤ 0.01" (picked, per class); "Calls need
  q ≤ 0.01" (de novo subgroup); "**Neither:** a call needs ≥ 2 distinct peptides and carries
  `class_fdr = unestimated`."
- §9.3 `tables/novel_orfs.tsv`: "`class_fdr` (float or `unestimated`), `called`, `rank`".
- §9.4 Novel-ORF flag `class_fdr_unestimated`.
- §8.1 row "Novel ORF | ORF | … | 0.01" through `stats.alpha.novel_orfs` (never a literal 0.01).
- M7-11 rules 4 and 5 and lead decisions 6–9 and 14. Decision 6 assigns the wiring to this issue.
  Decision 10 says "M7-11b step tests pin decisions 7–9".
- M7-12 decision 26 and M7-13 Q9: the picked class q cannot fall below 1/T, so a class with
  fewer than about 1/α candidates cannot be called. The run states this as a warning. This mirrors
  the gene-level `q_floor_above_alpha` warning (Khoa, 2026-10-22: "flag it"; q and α unchanged).

## What exists today (main at 9da78562 + #237)

- `apply_evidence` writes `class_fdr = "unestimated"` for every ORF and calls an ORF when it has
  ≥ `min_peptides_without_decoys` peptides and meets the class rule (`evidence.py`, M7-09 rule 8).
- `orf_table` already adds `class_fdr_unestimated` exactly when the cell is `unestimated`, and it
  formats a float cell with `%.6g` (`outputs._class_fdr`).
- `giap.novel.fdr` (#237) provides `orf_regimes`, `picked_class_fdr` (q, reason, `class_mismatch`,
  `ClassStats`) and `denovo_class_q`, but nothing calls them.
- The `novel_orfs` step takes no `peptides`, `features`, hit tables or π̂. `placement` computes
  π̂ per de novo peptide and E at the best target and decoy cost, but persists neither.

## Public interface

### `src/giap/novel/class_call.py` (new, pure)

```python
Reason = Literal["no_decoys", "decoys_not_makedb", "untagged_class", "no_makedb_pair",
                 "paired_share_low", "decoys_unscored", "no_target_score", "lost_pick",
                 "no_target_hit", "denovo_inputs_missing"]

@dataclass(frozen=True)
class OrfClassFdr:
    orf_id: str
    regime: Regime                 # giap.novel.fdr.Regime
    fdr_class: str                 # the list it competed in: the entry tag (picked), else ORF class
    q: float                       # NaN unless estimated and won (picked) / has a target hit (de novo)
    reason: Reason | None          # None iff q is finite
    class_mismatch: bool           # decision 9

@dataclass(frozen=True)
class ClassSummary:
    regime: Literal["picked_decoy", "denovo_subgroup"]
    cls: str
    estimated: bool
    reason: str | None             # ClassStats.reason (picked), else None
    n_targets: int                 # picked: ORFs entering the list (won + lost_pick); de novo: target peptides
    n_target_wins: int | None      # picked only
    n_decoys: int                  # picked: decoy wins; de novo: decoy peptides of the class
    paired_share: float | None     # picked only
    c: float | None                # de novo only (ADR-0010)
    q_floor: float | None          # picked 1/n_target_wins, de novo c/n_targets; None when 0 targets
    floor_above_alpha: bool        # q_floor is not None and q_floor > alpha

@dataclass(frozen=True)
class DenovoInputs:
    target_hits: pd.DataFrame      # peptide, e_target (de novo peptides with a finite e_target)
    decoy_hits: pd.DataFrame       # HIT_COLUMNS rows, search = decoy, de novo peptides, + e_decoy
    pi: Mapping[str, float]        # de novo peptide -> pi_denovo (NaN dropped: missing = 1)
    classifier: Classifier         # giap.mapping.loci.LocusClassifier

@dataclass(frozen=True)
class ClassFdrResult:
    orfs: dict[str, OrfClassFdr]           # every input ORF
    classes: tuple[ClassSummary, ...]      # sorted by (regime, cls)
    n_decoy_rows: int                      # is_decoy rows in peptides
    decoys_not_makedb: bool                # n_decoy_rows > 0 and no class-tagged decoy ID

def makedb_id_map(peptides: pd.DataFrame, gm: GeneModel, genome: Mapping[str, str],
                  decoy_prefixes: Sequence[str]) -> pd.DataFrame: ...
def denovo_inputs(costs: pd.DataFrame, hits: pd.DataFrame, classifier: Classifier) -> DenovoInputs: ...
def class_fdr(orfs: Sequence[Orf], peptides: pd.DataFrame, *, gm: GeneModel,
              genome: Mapping[str, str], cfg: GiapConfig,
              denovo: DenovoInputs | None) -> ClassFdrResult: ...
def class_fdr_cell(r: OrfClassFdr) -> float | str: ...   # q | 1.0 (lost_pick) | "unestimated"
```

### Changes to existing modules

- `giap.novel.fdr`: add `denovo_class_fdr(...) -> DenovoClassFdr(q, classes)`, where
  `classes[cls] = DenovoClassStats(n_targets, n_decoys, w_decoys, c)`. It is computed from the
  same lists as today. `denovo_class_q` becomes `denovo_class_fdr(...).q`, so its results do not
  change. `picked_class_fdr` is not changed.
- `giap.makedb.build`: add a public `nterm_entries(gm, genome)`, an alias of `_nterm_entries`.
- `giap.novel.evidence`: factor out `class_rules(table, orfs, *, gene_table, gm, genome, cfg,
  peptide_intervals) -> dict[str, bool]` (the M7-09 rule-4 `ok`), which `apply_evidence` uses.
  `apply_evidence` stops overwriting `class_fdr`: it keeps the cell it receives. It then applies
  the call rule below.
- `giap.novel.outputs.FLAG_ORDER` and `schemas/tables/novel_orfs.schema.json` (flag pattern) gain
  `class_mismatch`, placed after `class_fdr_unestimated` (Q3).
- `steps.mapping.placement` gains the output `placement_costs` (`placement_costs.parquet`). It is
  the `best_costs` frame that the step already builds, with `pi_denovo` filled, written as is.

### `novel_orfs` step (registry)

```
inputs          = [novel_candidates, gene_table, genes, replicons, genome]          (unchanged)
optional_inputs = [novel_promoter_windows, rbs_novel, promoter_windows, terminators,  (unchanged)
                   peptides, features, hits_exact, hits_tolerant, placement_costs]   (new)
outputs         = [novel_orfs, novel_orfs_gff3, novel_orfs_faa, annotation_patch,     (unchanged)
                   novel_class_fdr]                                                  (new)
```

The new inputs are optional, so the existing M7-09 step tests (which pass none of them) keep the
pass-1 behaviour byte for byte (acceptance criterion 3).

Step body, after `second_round_scores`:

1. When `peptides` is absent, every ORF gets `regime = unestimated` and `reason = no_decoys`,
   which is today's output. Otherwise `gm` is built with `features`, as `novel_candidates` does
   (`_gene_model(ctx)`).
2. `denovo` is built only when some ORF's regime is `denovo_subgroup`. It needs `hits_exact`,
   `hits_tolerant` (concatenated, as `placement._hits` does) and `placement_costs`. If any of the
   three is missing, the de novo ORFs get `reason = denovo_inputs_missing`, and the step logs
   one WARNING `novel_orfs.layer_missing` with flag `layer_missing:placement_costs` (or the
   missing key). The flag token is already registered.
3. `res = class_fdr(orfs, peptides, gm=gm, genome=genome, cfg=cfg, denovo=denovo)`.
4. In `scored` (from `minimal_call`), set `class_fdr` to `class_fdr_cell(res.orfs[id])` and add
   `class_mismatch` to `flags` where it is set. `orf_table` derives `class_fdr_unestimated`.
5. `class_rules(...)`, then `apply_evidence(...)`, then the writers, as today. Write
   `novel_class_fdr` (below).
6. Events and warnings (below).

## Rules

1. **Regime per ORF** = `orf_regimes(peptides, orfs, cfg.input.decoy_prefixes)` (M7-11 rule 1,
   decision 1). There is no cascade: an ORF in the picked regime whose picked result is
   unestimated does not fall back to its de novo peptides (Q6).
2. **Picked.** `picked_class_fdr(picked_orfs, peptides, makedb_id_map(...),
   decoy_prefixes=cfg.input.decoy_prefixes, replicon_lengths={id: rep.length for every
   replicon})`. Lengths are passed for linear replicons too (decision 14).
   - `makedb_id_map` returns `ID_MAP_COLUMNS` rows with `kind = orf`, one per distinct
     class-tagged ID in `protein_ids` (targets and stripped decoys). `class`, `contig`, `start`,
     `end` and `strand` come from `parse_orf_id`, converted to 1-based closed.
   - `locus_tag` is set only for `nterm_extension` IDs. It is the locus tag of the
     `nterm_entries(gm, genome)` record with that `db_id`, so it reproduces makedb on the run's
     own annotation. With no such record it is NA, and the ORF becomes `no_makedb_pair`
     (decision 7).
   - Consequence: every tagged decoy is in the map, so the paired share is 1 in v1. The 80 % guard
     is reachable only with a real `db_id_map.tsv` (Q5).
3. **De novo.** `denovo_class_fdr(denovo_orfs, d.target_hits, d.decoy_hits, d.classifier, d.pi)`.
   - `denovo_inputs` restricts both hit frames and `pi` to `source = denovo` peptides, so the
     decoys of database peptides never enter a de novo list (Q11).
   - It attaches `e_decoy` (and `e_target`) from `placement_costs` by peptide.
     `denovo_class_fdr` keeps only each peptide's best-cost decoy rows, so the E at `m_decoy` is
     the right one.
   - A NaN E on a kept row raises `ValueError` (M7-11 decision 11). The step does not catch it.
4. **Unestimated regime.** The reason is `decoys_not_makedb` when the run has decoy rows,
   otherwise `no_decoys`.
5. **Cell value** (`class_fdr_cell`):
   - a finite q gives that float (written `%.6g`);
   - `lost_pick` gives **1.0**: the class is estimated, but the ORF lost to its decoy (Q1);
   - every other reason gives `"unestimated"`.
6. **Call** (replaces M7-09 rule 8 in `apply_evidence`). Let `rule` be `class_rules[orf_id]`
   and α be `cfg.stats.alpha.novel_orfs`.
   - Numeric `class_fdr` = q: `called = rule and q ≤ α`, compared as floats before formatting.
     There is no peptide minimum (§7.6.4). The class rule of `cds_outframe` and `antisense`
     still asks for 2 peptides.
   - `class_fdr = "unestimated"`: `called = rule and n_peptides ≥
     novel_orfs.min_peptides_without_decoys`.
   - `rank` is unchanged (it covers called ORFs only, M7-09 rule 10).
7. **Class mismatch** (decision 9).
   - The ORF competes in its entry's class (`fdr_class`).
   - Its TSV `class` stays the reconstructed class, and its class rule uses that class (Q9).
   - It carries the flag `class_mismatch`.
8. **Class q floor.**
   - Picked: `q_floor = 1 / n_target_wins`. This is the smallest picked q: with no decoy win
     above the last target, FDR = (0/f + 1)/T.
   - De novo: `q_floor = c / n_targets`, the ADR-0010 floor with W = 0.
   - `floor_above_alpha = q_floor > α`, a float comparison, so T = 100 at α = 0.01 does not warn.
   - A class with 0 targets (or 0 target wins) has `q_floor = None` and no warning. Unestimated
     classes have no floor.
9. **Events** (`ctx.log`, with `extra={"event": …, "data": …}`, which feed `events.jsonl` and
   `logs/novel_orfs.log`):
   - `novel_orfs.class_fdr`, INFO, one per `ClassSummary`. Its `data` holds every
     `ClassSummary` field, plus `alpha` and `n_called`. `cls` is written as `class`. NaN
     values are written as null.
   - `novel_orfs.class_q_floor_above_alpha`, WARNING, one per summary with
     `floor_above_alpha`. Data: `{regime, class, n_targets, q_floor, alpha}`. Message:
     `"novel_orfs: class <class> (<regime>) has <T> targets; its q floor <q_floor> > alpha
     <alpha>, so none can be called (a class with fewer than about 1/alpha = <1/α> candidates
     cannot be called at alpha = <α>)"`.
   - `novel_orfs.decoys_not_makedb`, WARNING, once, when `res.decoys_not_makedb`. Data:
     `{n_decoy_rows}`. It fires whether or not the run has any ORF.
   - The two warnings reach the manifest `warnings` list (`novel_orfs:<event>`). They are not
     run-flag tokens (Q4), and they make `--strict` exit 3 ("warnings were raised and `--strict` was given").
10. **Manifest.** `_content` gains `"novel_orfs": {"classes": [...]}`, the de-duplicated
    `novel_orfs.class_fdr` event data (`_event_data`), sorted by `(regime, class)`. Runs without
    the step get `"novel_orfs": null`. `run_manifest.schema.json` gains the block (Q7).
11. **`novel_class_fdr` table** (work output, published at `novel/novel_class_fdr.parquet`, Q2).
    - One row per ORF, in `novel_orfs.tsv` row order.
    - Columns: `orf_id` str, `class` str, `fdr_class` str, `regime` str, `q` float64 (NaN),
      `class_fdr` str (the TSV cell), `reason` str|null, `class_mismatch` bool, `class_rule_met`
      bool, `called` bool.
    - Written with `write_parquet_deterministic`. M7-13 G4 reads `q` and `class_rule_met` from
      it, because at q ∈ {0.05, 0.10} `called` cannot tell whether the class rule was met.

## Nextflow wiring

- `mapping.nf`:
  - PLACEMENT's output list gains `'placement_costs=placement_costs.parquet'`.
  - `emit:` gains `hits_exact = pick(MAP_EXACT.out.outputs, 'hits_exact.parquet')` and
    `placement_costs = pick(P, 'placement_costs.parquet')`. `hits_tolerant` is already emitted.
- `novel.nf`:
  - `take:` gains `hits` (`[hits_exact, hits_tolerant, placement_costs]`, a value channel).
  - The `novel_orfs` stepIn gains `kv('peptides', peptides)`, `kv('features', …)` (as
    `novel_candidates` gets it), and `kv('hits_exact', …)`, `kv('hits_tolerant', …)`,
    `kv('placement_costs', …)`. None of them goes into the required-keys list.
  - The output list gains `'novel_class_fdr=novel/novel_class_fdr.parquet'`, and `emit:` gains
    `novel_class_fdr`.
- `main.nf`: pass `MAPPING.out.hits_exact`, `.hits_tolerant` and `.placement_costs` to NOVEL, and
  publish `kv('novel_class_fdr', NOVEL.out.novel_class_fdr)`.
- `runner/local.py` `PUBLISHED_FILES` gains `"novel_class_fdr": "novel/novel_class_fdr.parquet"`.
  `test_local_publish_names.py` pins it against `novel.nf`.
- Stub runs: the stub writes the Parquet with its header schema, as the other Parquet stubs do.

## Examples

ORFs and entries are given 1-based in IDs. Picked scores use the #237 fixture convention (rank r:
`q_value = r·1e-5`, `score = 1000 − r`). α = 0.01 unless stated, and every ORF meets its class
rule unless stated. "Event" means `novel_orfs.class_fdr` data.

| # | Input | Expected |
|---|---|---|
| 1 | makedb run; 150 `intergenic` ORFs, each paired with its own `GIAPORF_intergenic_c1_…_p` entry, ranks 1–150; 30 decoy-only entries `rev_GIAPORF_intergenic_…`, ranks 151–180 | all 150: `class_fdr` **`0.00666667`**, `called` true, no `class_fdr_unestimated`; event `{regime: picked_decoy, class: intergenic, estimated: true, n_targets: 150, n_target_wins: 150, n_decoys: 30, paired_share: 1.0, q_floor: 0.006666…, floor_above_alpha: false, n_called: 150}`; no WARNING |
| 2 | as 1 with 50 targets and no decoys | all 50: **`0.02`**, `called` false; one WARNING `novel_orfs.class_q_floor_above_alpha` `{regime: picked_decoy, class: intergenic, n_targets: 50, q_floor: 0.02, alpha: 0.01}` |
| 3 | as 2 with exactly 100 targets | all 100: **`0.01`**, `called` true; no WARNING (0.01 > 0.01 is false) |
| 4 | `cds_outframe` ORF, + strand, stop codon at 0-based 400 (`GIAPORF_cds_outframe_c1_101_400_p`), 3 peptides, best target rank 2; decoy row `rev_GIAPORF_cds_outframe_c1_101_400_p` at rank 1 (or tied: same q and score) | `class_fdr` **`1`**, `called` false, `reason` `lost_pick`, flags without `class_fdr_unestimated` (the ≥ 2-peptide rule is **not** applied) |
| 5 | search against Bakta proteins, decoys `rev_b0001`, `rev_b0002` only (7 decoy rows); one `cds_outframe` ORF with 2 peptides and one `intergenic` ORF with 1 peptide | both `unestimated` + `class_fdr_unestimated`, `reason` `decoys_not_makedb`; cds_outframe `called` true, intergenic `called` false; one WARNING `novel_orfs.decoys_not_makedb` `{n_decoy_rows: 7}`; no `novel_orfs.class_fdr` event |
| 6 | database input with no decoy row (the M7-09 step fixture) | `novel_orfs.tsv`, `.gff3`, `.faa` and patch **byte-identical** to main; every side-table row `regime` `unestimated`, `reason` `no_decoys`; no new WARNING |
| 7 | makedb run; a `sorf` ORF with 2 peptides and a `pseudogene` ORF with 1 peptide | both `unestimated`, `reason` `untagged_class`; sorf `called` true, pseudogene `called` false |
| 8 | reconstructed `intergenic` ORF, c1 +, stop codon at 0-based 400, matched to entry `GIAPORF_cds_outframe_c1_101_400_p` | side table `class` `intergenic`, `fdr_class` `cds_outframe`; TSV `class` `intergenic`, flags contain `class_mismatch`; its q comes from the `cds_outframe` list; its class rule is the intergenic one |
| 9 | de novo only; 120 `intergenic` ORFs, one target peptide each; no decoy hit classifies as intergenic; `pi_denovo` = 0.5 for every de novo peptide | c = 0.5; all 120: **`0.00416667`**, `called` true; event `{regime: denovo_subgroup, class: intergenic, n_targets: 120, n_decoys: 0, c: 0.5, q_floor: 0.004166…}` |
| 10 | de novo only; 40 `antisense` target peptides (20 ORFs × 2), `pi_denovo` missing for all (= 1) | c = 1, q_floor 0.025 > 0.01: WARNING `class_q_floor_above_alpha` `{regime: denovo_subgroup, class: antisense, n_targets: 40, q_floor: 0.025, alpha: 0.01}`; no antisense ORF called |
| 11 | makedb run; gene `b0002` (c1, +) whose `nterm_entries` record is `GIAPORF_nterm_extension_c1_<s1>_<e1>_p`; an `nterm_extension` ORF with `overlapping_gene = b0002` and peptides on that entry; a second run whose gene model lacks `b0002` | run 1: paired by locus tag, numeric `class_fdr`; run 2: `unestimated`, `reason` `no_makedb_pair` |
| 12 | mixed run (makedb decoys and de novo rows); ORF A has database peptides, ORF B only de novo, both `intergenic` | A `regime` `picked_decoy`, B `denovo_subgroup`; two `novel_orfs.class_fdr` events for `intergenic`, one per regime; B's de novo list holds only de novo ORFs' peptides (lists not pooled, decision 1) |
| 13 | as 1 with `stats.alpha.novel_orfs = 0.005` | all 150 `0.00666667`, `called` false; WARNING with `alpha: 0.005` |

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no candidate ORF | empty TSV with header, empty side table with schema; pragma-only GFF3s (as today); no class event; `decoys_not_makedb` still fires when applicable |
| `peptides` input absent | as example 6 (all `no_decoys`); no de novo inputs are read |
| de novo ORFs but `placement_costs` (or a hit table) absent | those ORFs `unestimated`, `denovo_inputs_missing`; WARNING `novel_orfs.layer_missing` with `layer_missing:<key>` |
| de novo ORF whose peptides have no finite `e_target` | `unestimated`, `no_target_hit` |
| NaN `e_decoy` on a best-cost decoy row | `ValueError` (decision 11); the step fails |
| class whose every target loses the pick | all `1`, none called; `n_target_wins` 0 → `q_floor` null, no floor WARNING |
| class with targets and zero decoys | picked q = 1/T, de novo q = c/T (examples 1–3, 9) |
| target and decoy tie on (q, score) | decoy wins (decision 5) → `lost_pick` → `1` |
| decoy rows with NaN q in a class | class ranked on score; with NaN score too → `decoys_unscored` → `unestimated` (decision 11) |
| aggregated peptide that is target in one row and decoy in another (`is_decoy` false after `all`) | a target; its `rev_` IDs are ignored (`_decoy_tags` reads `is_decoy` rows only) |
| two ORFs paired to one entry | both get the entry's q (one list member) |
| origin-crossing ORF on a circular replicon | key wrapped with the replicon length (decision 8), paired normally |
| open-3′ ORF on a linear contig, + and − | key `stop ± 3` (decision 13), wrapped (decision 14) |
| minus-strand ORF | key = lowest 0-based coordinate of the stop codon, `GIAPORF_…_m` ID: same rule as + |
| `cterm_extension` | untagged → `unestimated`; its class rule is false, so never called |
| α changed in config | every comparison uses `stats.alpha.novel_orfs` (example 13) |
| `--strict` with a floor or `decoys_not_makedb` warning | exit 3 (warnings were raised) |
| resume duplicates a class event | `_event_data` de-duplicates it; manifest unchanged |

## Line estimate (counted: `src`, `workflow`, `scripts`, added + deleted)

| Part | Lines |
|---|---|
| `novel/class_call.py` (dataclasses, `makedb_id_map`, `denovo_inputs`, `class_fdr`, floor, cell) | 130 |
| `novel/fdr.py` `denovo_class_fdr` + stats, wrapper | 45 |
| `novel/evidence.py` `class_rules` extraction + call rule | 20 |
| `novel/outputs.py` flag, schema pattern, `makedb/build.py` alias | 6 |
| `steps/novel.py` inputs, de novo loading, classifier, events, side table | 70 |
| `steps/mapping.py` `placement_costs` | 8 |
| `manifest.py` + `run_manifest.schema.json` | 28 |
| Nextflow (`mapping.nf`, `novel.nf`, `main.nf`) + `runner/local.py` | 22 |
| **Total** | **≈ 330** |

Tests (≈ 450, not gated): unit `class_call`, step tests, manifest block, nf-test stub, and the
local-vs-Nextflow parity. One PR fits the 400 warning. If review grows it past 400, split off
**M7-11b-2** (the manifest block and the `novel_class_fdr` publication, ≈ 45 lines). M7-13b needs
only the TSV and the work-dir side table.

## Acceptance criteria

1. `tests/unit/novel/test_class_call.py`: examples 1–5 and 7–13 at the `class_fdr` +
   `class_rules` + `apply_evidence` level, and every degenerate row except those marked as step
   rows. It checks the exact TSV cell strings and the `called` values.
2. `tests/unit/steps/test_novel_class_fdr_step.py`: decisions 7–9 through the step (pinned by
   decision 10), examples 5, 6, 10 and 12, the three events and two WARNINGs with exact `data`,
   the missing-input degenerate rows, and the side-table columns, dtypes and row order.
3. Every existing test passes unchanged, in particular `tests/unit/novel/test_outputs.py`,
   `test_candidates.py`, `test_novel_evidence.py`, `tests/unit/steps/test_novel_*`, and
   `tests/integration/test_novel_local_vs_nextflow.py`. Example 6 is the byte-identity check.
4. `novel_orfs.tsv` with `class_mismatch` and numeric `class_fdr` validates with
   `frictionless validate` against the updated schema. The manifest with the `novel_orfs` block
   validates against `run_manifest.schema.json`.
5. `workflow/tests/novel.nf.test` and `main_novel.nf.test` (stub) see `novel_class_fdr` published.
   The local-vs-Nextflow test is extended to a giapsim `novel_orf_power` 100k run with makedb
   decoys: the four tables and the side table are byte-identical between the executors.
6. Smoke check (integration, `slow`), one `novel_orf_power(seed, "100k")` run, database input:
   - `intergenic`, `cds_outframe`, `antisense` and `nterm_extension` have no `unestimated` ORF
     other than `no_makedb_pair` caused by a gene-model miss. Zero is expected, since M7-13
     treats it as a bug.
   - At least one class event has `estimated: true`.
   - For every summary, `floor_above_alpha` equals the floor computed from its own
     `n_target_wins`.
7. `docs/methods/novel-orfs.md` §"Class FDR" gains:
   - "Call rule" (rule 6, Q1's `lost_pick` = 1 and why it must not fall back);
   - "Class q floor" (rule 8, with the sentence "a class with fewer than about 1/α candidates
     cannot be called at α");
   - the v1 paired share = 1 consequence (rule 2).
   The statistician signs off.
8. Backlog: the M7-11b line, line 688 (the GIAP-warning half; the gate-record half stays with
   M7-13c), and M7-11 nit 1 are closed or re-pointed. Q5 and Q10 go to the backlog.
9. Lint, format, mypy clean; ruff passes on the new module; `check_imports.py` clean.

## Questions for the lead

1. **`lost_pick` cell value.** The ORF's class is estimated, but its decoy beat it, so picked FDR
   assigns it no q. **Proposed default:** `class_fdr = 1` and not called. Writing `unestimated`
   would send it to the ≥ 2-peptide rule, so ORFs beaten by their decoys would be called and
   G4's FDP would inflate. An alternative is the class q at the decoy's rank, but it is
   undefined for a target that never enters the list.
2. **Per-ORF diagnostics (`regime`, `reason`, `fdr_class`, `class_rule_met`) as `novel_orfs.tsv`
   columns — pending Khoa** (adds §9.3 columns; ADR-0012 and ADR-0036 needed Khoa).
   **Reversible default:** the TSV keeps exactly the §9.3 columns. The diagnostics live in the
   work/published side table `novel/novel_class_fdr.parquet`, which M7-13b reads. If Khoa adds
   the columns, the side table becomes redundant and is dropped in one PR.
3. **`class_mismatch` flag** (decision 9 says "column"). **Proposed default:** a §9.4 Novel-ORF
   flag added through the schema, like the M7-06/07/09 additions (`no_start_codon`, `partial`, …),
   plus the side-table bool. If the lead holds that a flag beyond §9.4 needs Khoa, treat it as
   pending Khoa: the side-table bool alone is the reversible fallback.
4. **Floor warning as a run-flag token** (`RUN_FLAGS`; ADR-0041 shows tokens beyond §9.4 need
   Khoa). **Proposed default:** WARNING event and manifest `warnings` only, no token. This
   mirrors the gene-level `q_floor_above_alpha`, which is also a warning and not a token. If a
   token is wanted, it is pending Khoa.
5. **Real `db_id_map.tsv` input.** v1 derives `nterm_extension` locus tags from the run's own gene
   model and treats every tagged decoy as paired (paired share 1). A user-supplied map would make
   the 80 % guard live, but it needs a new config key (`defaults.yaml` is protected; adding a key
   needs an ADR). **Proposed default:** no input in v1; a backlog line names the guard as
   unreachable from the pipeline.
6. **No picked → de novo cascade** for a picked-regime ORF whose picked result is unestimated
   (for example `no_makedb_pair`) but which also has de novo peptides. **Proposed default:** no
   cascade. That keeps one regime per ORF and keeps the de novo lists free of database-first
   ORFs. Revisit if M7-13's mixed arm shows losses.
7. **Manifest block** (M7-11 rule 5 says "and to the manifest"). **Proposed default:** include
   it (≈ 28 lines). If the PR passes 400, defer it with the publication to M7-11b-2. Manifest
   goldens that now gain `"novel_orfs": null` are updated by the lead, not the implementer
   (role guard).
8. **M7-13 runs and `--strict`.** Floor warnings are expected at 1M for intergenic and for small
   de novo classes. **Proposed default:** M7-13 never passes `--strict`, and the gate harness
   ignores the warning (it computes vacuity itself, M7-13 rule 6).
9. **TSV `class` of a mismatched ORF.** **Proposed default:** keep the reconstructed class (GIAP's
   classification, used by the class rule, ranking and patch). The entry class appears only as
   `fdr_class`. M7-13 rule 5 then counts it as FP + FN, as decided there.
10. **Report caveat** (MVP-04 gate 15 mentions a q-floor caveat in the report). **Proposed
    default:** out of scope. A backlog line asks M9 to render `class_q_floor_above_alpha`
    warnings on the V13 card.
11. **De novo decoy population in mixed runs.** **Proposed default:** only de novo peptides'
    decoy hits and π̂ enter the de novo lists, and c is the maximum π̂ over de novo peptides.
    This leaves M7-11 NIT 1 (c over hit peptides only) where it is, and the backlog line stays
    open.

## Lead decisions (2026-10-23)
1. **Q1, Q4, Q6–Q11:** the proposed defaults are accepted as written.
2. **Q2, new `novel_orfs.tsv` columns: pending Khoa.** The reversible default is in use: the TSV
   keeps exactly its §9.3 columns, and `regime`, `reason`, `fdr_class` and `class_rule_met` live in
   `novel/novel_class_fdr.parquet`. Listed in PROGRESS.md › Questions for Khoa.
3. **Q3, `class_mismatch` as a novel-ORF flag: pending Khoa.** It is a flag beyond §9.4. For now
   it is only the side-table bool, and no flag is added. If Khoa accepts the flag, a one-line follow-up
   adds it.
4. **Q5:** no `db_id_map.tsv` input in v1. A backlog line records that the 80 % pairing guard cannot
   fire from the pipeline.
5. **Order:** test-author on a branch stacked on m7/novel-class-fdr (#237, in the merge queue),
   then the implementer, then the statistician and code review. The PR is labelled stats. If the
   PR goes past 400 counted lines, the manifest block and side-table publication split out as
   M7-11b-2.
