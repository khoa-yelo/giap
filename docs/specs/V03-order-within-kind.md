# V03-order: rank V3 loci within their own kind, not across scales

Backlog item: "V03 locus order mixes scales" (docs/backlog.md, 2026-10-25 afternoon) · depends
on: M9-10 (merged) · estimate ≈ 40 code lines (+ ≈ 200 test lines) · **`stats`** (it changes which
loci a reader sees first and which loci the cap renders) · roles: statistician (review of the rule
before the test step), test-author (new tests), lead (edits to existing tests, snapshot baseline),
implementer (`src/giap/viz/locus.py`, user guide, changelog)

## Problem

`src/giap/viz/locus.py::_order` sorts the whole V3 pool by `Locus.error` ascending, NA last, then
kind (`KIND_ORDER`), then `locus_id` (M9-10 lead decisions 3, 5 and 8). `error` is:

| Kind | `error` | Family (GIAP_SPEC.md §8.1) | Level |
|---|---|---|---|
| prophage, bgc | `regions_stats.q_value` (`o_minus_e`, `subset = all`) | region effect, BH within type | `stats.alpha.context` (0.05) |
| island | `islands.p_value` (`q_value` is NA, §9.3 `islands.tsv`) | islands, FWER over windows (max statistic) | `stats.alpha.context` (0.05) |
| operon_flag | `operon_flags.q_value` (pool keeps only `q ≤ alpha`) | operon discordance / missed operon, BH or decoy FDR | `stats.alpha.operon_flags` (0.10) |

A BH q-value, a family-wise p-value and a q-value from a family with a different level are not on
one scale. A q of 0.004 and an FWER p of 0.004 do not mean the same evidence. The two families
also have different sizes, and the flag family is controlled at twice the level. So "ISL_1 before
prophage_1 because 0.005 < 0.0197" has no meaning. Under the cap (`report.max_locus_plots`, 50) the
mixed sort can also hide a whole kind: 60 flags with q < 0.03 push every prophage and BGC out of
the rendered set (example 4).

The user guide already says this (`docs/user-guide/plots-and-report.md`: "so one ordering mixes
two scales").

## Spec sections

- GIAP_SPEC.md §10.2 V3: "At most `report.max_locus_plots` (50), ordered by q". Islands have no q
  (`q_value = NA`, §9.3 `islands.tsv`), so a literal reading would put every island in the NA tail.
  §10.2 does not say how to order across kinds. This spec fills that gap. It does not override
  §10.2: within each kind the order is still by that kind's own q, or FWER p for islands.
- §8.1 claims and their error control (the family levels above); `stats.alpha.context`, `stats.alpha.operon_flags`
  (§5.2).
- M9-10 lead decisions 3, 5, 8 (this spec replaces 3 and refines 8; 5 is unchanged).

## Proposed rule (for the statistician)

Each locus gets the sort key `(tier, kind_rank, kind_index, locus_id)`:

1. **tier**:
   - 0 when `error ≤ alpha(kind)`. The comparison is a plain float `≤` with no tolerance. `alpha`
     is `stats.alpha.context` for prophage, bgc and island, and `stats.alpha.operon_flags` for
     operon_flag. Both are read from the `cfg` passed to `select_loci`.
   - 1 when `error > alpha(kind)`.
   - 2 when `error` is NA (`None`).
2. **kind_rank**: the 1-based position of the locus among the loci of the **same kind and the
   same tier**, sorted by `(error, locus_id)`. In tier 2 it is 0 for every locus, so tier 2 sorts
   by kind, then ID. This is decision 8 unchanged.
3. **kind_index**: `KIND_ORDER.index(kind)`, so prophage < bgc < island < operon_flag.
4. **locus_id**: plain string order, as today.

In words:
- Significant loci come first, each at its own family's level.
- Within a tier the kinds take turns: the best of each kind, then the second best of each kind,
  and so on.
- A kind with no loci left in the tier drops out of the rotation.
- No step ever compares an error from one family with an error from another. Tier membership
  compares an error only with its own family's level.

Unchanged:
- Pool membership: the same rows, filters and merging (operon flags merge by TU, and a merged
  locus keeps its best flag's q, decision 6).
- `select_loci` returns the whole ordered pool (decision 5).
- The first `report.max_locus_plots` loci are rendered.
- The index TSV columns (decision 4), labels and titles.
- Novel ORFs stay out of the pool (decision 1).

Changed:
- The order of `select_loci`, and with it the index TSV `rank`, which loci the cap renders, and
  the order of the rendered instances in the plot result.
- The index figure text names the estimate: `1. prophage_1 [prophage] q 0.0197`,
  `2. ISL_1 [island] FWER p 0.005`, `NA` when missing. Its content is not pinned (decision 4), so
  this is free. It tells the reader that the column is not one scale.
- The docstrings of `select_loci` and `Locus.error`, and the V03 bullet of
  `docs/user-guide/plots-and-report.md`.

Why not the alternatives:
- **Kind first** (all prophages, then BGCs, ...): simple, but under the cap the last kind
  (operon flags) can be starved. Non-significant regions (regions are not filtered by alpha,
  M9-10) would also come before every significant flag.
- **Plain interleave without tiers**: a region at q 0.9 would be rendered before the second
  island at FWER p 0.001.
- **Converting to one scale** (for example, BH over all p-values): this would be a new statistic
  on the plot layer. It also cannot be done for islands, whose per-window p-values are not in the
  table.

## Examples

`alpha.context = 0.05`, `alpha.operon_flags = 0.10` unless stated otherwise. "Now" is the order of
the current `_order`.

| # | Pool (kind, error) | Expected order | Now |
|---|---|---|---|
| 1 | prophage_1 q 0.003; ISL_1 p 0.01; PAIR_7 q 0.02; prophage_2 q NA (M9-10 ex. 1) | prophage_1, ISL_1, PAIR_7, prophage_2 | same |
| 2 | prophage_1 q 0.04; bgc_1 q 0.001; ISL_1 p 0.0005; ISL_2 p 0.002; PAIR_1 q 0.08 | prophage_1, bgc_1, ISL_1, PAIR_1, ISL_2 | ISL_1, bgc_1, ISL_2, prophage_1, PAIR_1 |
| 3 | prophage_1 q 0.06 (not significant at 0.05); PAIR_1 q 0.09 (significant at 0.10) | PAIR_1, prophage_1 | prophage_1, PAIR_1 |
| 4 | PAIR_001 … PAIR_060 with q = i / 10 000 (each its own TU); prophage_1 q 0.03; bgc_1 q 0.04; cap 50 | prophage_1, bgc_1, PAIR_001 … PAIR_060; rendered: prophage_1, bgc_1, PAIR_001 … PAIR_048; 62 index rows, the last 12 `rendered = false` | the 60 flags, then prophage_1, bgc_1; rendered PAIR_001 … PAIR_050 only |
| 5 | ISL_2 p 0.005; ISL_1 p 0.005 (tie); prophage_1 q 0.01 | prophage_1, ISL_1, ISL_2 | ISL_1, ISL_2, prophage_1 |
| 6 | prophage_2 q NA; bgc_1 q 0.9; ISL_2 p NA | bgc_1, prophage_2, ISL_2 | same |
| 7 | bgc_1, prophage_2, prophage_1, ISL_1, PAIR_1, all at 0.05 (ties at the level) | prophage_1, bgc_1, ISL_1, PAIR_1, prophage_2 | prophage_1, prophage_2, bgc_1, ISL_1, PAIR_1 |
| 8 | prophage_1 q 0.01; prophage_2 q 0.02; prophage_3 q 0.5; prophage_4 q 0.7; ISL_1 p 0.001; bgc_1 q 0.3 | prophage_1, ISL_1, prophage_2, prophage_3, bgc_1, prophage_4 (tier 1 restarts the rotation: prophage_3 and bgc_1 are both rank 1 there) | ISL_1, prophage_1, prophage_2, bgc_1, prophage_3, prophage_4 |
| 9 | snapshot fixture `sim_plots`: ISL_1 p 0.005; prophage_1 q 0.0197; bgc_1 q 0.0551; PAIR_176 q 0.881 (dropped, > 0.10) | prophage_1, ISL_1, bgc_1 | ISL_1, prophage_1, bgc_1 |
| 10 | `stats.alpha.context = 0.01`: prophage_1 q 0.03; ISL_1 p 0.005; ISL_2 p 0.02 | ISL_1, prophage_1, ISL_2 | ISL_1, ISL_2, prophage_1 |

## Edge cases

| Situation | Behaviour |
|---|---|
| empty pool (validation run, header-only tables) | `[]`; the empty state "no loci to show" is unchanged |
| a kind with no loci, or none in a tier | it drops out of that tier's rotation; ranks are not reserved for it (example 3: no prophage at rank 1 of tier 0, PAIR_1 is first) |
| only one kind | the order is `(tier, error, locus_id)` = today's order for that kind (M9-10 example 2, 60 prophages, is unchanged) |
| `error` NA, NaN or ±inf | `_num` already maps NaN and inf to `None`. All of these go to tier 2, which sorts by kind, then ID |
| `error == alpha(kind)` exactly | tier 0 (`≤`, the same comparison as the flag filter and `compare/stability.py`) |
| ties in `error` within a kind and tier | `locus_id` order (example 5) |
| equal `kind_rank` across kinds | kind order (examples 2, 7) |
| island with `p > alpha.context` (for example, `islands.tsv` produced at a looser level) | tier 1, like any other non-significant locus |
| operon flags | always tier 0, because the pool already drops `q > alpha.operon_flags`. A merged TU locus ranks by its best flag's q |
| `locus_id` digit order | plain string order (`prophage_10` < `prophage_2`), as today; not in scope |
| circular or linear contig, origin-crossing span, minus strand | no effect: the key uses only kind, error and ID, never coordinates or strand |
| the same `locus_id` in two kinds | kind order breaks the tie; the key stays total |
| input row order (`regions_stats`, `islands`, `operon_flags` shuffled) | no effect: the key is total, so the order is deterministic |
| `report.max_locus_plots = 0` | none rendered; the index lists all in the new order |
| `max_locus_plots` ≥ pool size | all rendered; only the order changes |

## Acceptance criteria

- [ ] `select_loci` returns the pool in the order of the "Proposed rule" key. Examples 1–10 hold
  in a new `tests/unit/viz/test_locus_order_within_kind.py`.
- [ ] **No cross-family comparison (metamorphic property).** Take one kind and apply any strictly
  increasing map to its errors that keeps each error on the same side of that kind's alpha. The
  order of the whole pool does not change. Applying the same kind of map to the errors of every
  kind at once also leaves the order unchanged.
- [ ] **Within-kind order (property).** For each kind, the subsequence of the result is sorted by
  `(error is None, error, locus_id)`.
- [ ] **Rotation (property).** In every prefix of the tier-0 block, the counts of any two kinds
  that still have tier-0 loci left after the prefix differ by at most 1. The same holds for
  tier 1.
- [ ] **Tiers.** Every tier-0 locus comes before every tier-1 locus, and every tier-1 locus comes
  before every NA locus. The alpha values come from the `cfg` passed in (example 10 with a
  `stats.alpha.context=0.01` override).
- [ ] **Determinism.** Shuffling the input rows of the three tables gives the same order.
- [ ] Pool membership, labels, spans, the index columns, the cap and the empty state are
  unchanged. Every existing test in `test_locus_select.py` except the two listed below passes
  without edits.
- [ ] The index TSV `rank` and `rendered` columns follow the new order (example 4: 62 rows, ranks
  1–62, `rendered` true for the first 50).
- [ ] The index figure text names the estimate per row (`q` for prophage, bgc and operon_flag;
  `FWER p` for island).
- [ ] The `docs/user-guide/plots-and-report.md` V03 bullet describes the new rule, and the
  sentence "one ordering mixes two scales" is removed. A `changelog.d/<PR>.bugfix.md` fragment is
  added.
- [ ] `ruff`, `ruff format --check` and `mypy src/giap` are clean.

## Test plan

### Existing tests the lead edits (test-authors and implementers do not)

1. `tests/unit/viz/test_locus_select.py::test_ties_break_by_kind_then_id`: the expected order
   becomes `["prophage_1", "bgc_1", "ISL_1", "PAIR_1", "prophage_2"]` (example 7). Rename it to
   `test_ties_break_by_kind_rank_then_kind_then_id`.
2. `tests/unit/viz/test_locus_select.py::test_property_order_is_error_na_last_then_kind_then_id`:
   its oracle is the old key. Either replace the oracle with the new key, with the alphas from
   `resolve_config()`, or delete the test in favour of the test-author's properties below. Lead's
   choice; the second avoids two oracles.
3. The module docstring of `test_locus_select.py` cites decision 3 implicitly ("order"). Add a
   pointer to this spec. Optional.
4. The snapshot baseline `tests/snapshot/baseline/sim_plots-V03_locus.png` (the index figure):
   - The first two lines swap (example 9), and the estimate name is added to each line, so the
     baseline changes.
   - The lead regenerates it on Slurm with `--mpl-generate-path` and Khoa approves it (§10.1:
     after M9 a baseline change is a reviewed change).
   - No test code changes.

Unchanged, checked against the new rule:
- `test_example1_order_across_kinds`, `test_na_errors_sort_last_then_by_kind`,
  `test_example2_sixty_candidates_fifty_rendered` (q = i/1000, so prophage_1 … prophage_50 are
  tier 0 because `50/1000 ≤ 0.05` in floating point, and prophage_51 … prophage_60 are tier 1),
  the cap tests and `test_example1_index_tsv_is_one_based_with_kind_and_na_error`.
- `test_locus_na_coords.py`, `test_locus_review.py` (it compares sets), `test_locus_render.py`
  and `test_real_number_format.py:252` (it compares sets).
- `tests/snapshot/test_plots.py`: the first `LOCUS_KEEP = 2` instances are still {prophage_1,
  ISL_1}, now in that order. The baselines `sim_plots-V03_locus_ISL_1.png` and
  `sim_plots-V03_locus_prophage_1.png` are unchanged, and the `sim_empty` baselines are unchanged.
- `tests/snapshot/report_structure.json`: unchanged. `io/plot_files.read_plot_files` lists the
  instance files by sorted name, so the report shows `V03_locus`, `V03_locus_ISL_1`,
  `V03_locus_prophage_1` in that order whatever the rank (see question 4).

### New tests (test-author), in `tests/unit/viz/test_locus_order_within_kind.py`

- Examples 2, 3, 4, 5, 7, 8 and 10 through `select_loci`, built with `tests/helpers/viz_locus`
  as in `test_locus_select.py`.
- Example 4 also through the plot function: 50 rendered instances in the expected order, and the
  index `rank` and `rendered` columns.
- Example 9 on the `sim_plots` snapshot fixture (`tests/helpers/snapshot_fixture`), `select_loci`
  IDs only.
- Hypothesis properties: the metamorphic property, the within-kind order, the rotation, the tier
  order and input-row shuffling (acceptance criteria). Pools mix all four kinds, with errors on
  both sides of each alpha, exactly at alpha, and NA.
- The index figure text: each line contains `q` or `FWER p` according to the kind.
- One test that the alphas come from `cfg` (example 10), and one with
  `stats.alpha.operon_flags=0.05`: a flag at q 0.06 is dropped, which is pool membership and
  unchanged, and a flag at q 0.04 stays tier 0.

The tests are fast unit tests with no Slurm, apart from the lead's baseline regeneration.

## Questions for the lead and statistician

1. **The rule.** Tiered round-robin (significant at the family's own level, then non-significant,
   then NA; kinds rotate within a tier). **Default:** as proposed. It is reversible: only `_order`
   and its tests change. Alternatives: kind-first blocks, or a plain rotation without tiers (see
   "Why not the alternatives").
2. **Tier levels.** Regions and islands use `stats.alpha.context`, and flags use
   `stats.alpha.operon_flags`. **Default:** these two keys. They are the same thresholds that
   `compare/stability.py` and the islands step use. Should islands instead use a hard 0.05
   (§7.3.5 "windows with p_w ≤ 0.05")? Today the islands step reads `alpha.context`, so the two
   agree.
3. **Index TSV columns.** Add an `estimate` column (`q` or `fwer_p`) and/or a `kind_rank` column
   to `V03_locus.tsv`? **Default:** no, because decision 4 pins the columns exactly. `kind` already
   tells a reader which estimate a row uses, and the figure text names it. Adding columns later
   only means appending them.
4. **Report order (found while writing this spec).** The HTML report lists V03 instances in
   file-name order (`read_plot_files` sorts the stems), not in rank order. So the order chosen here
   only decides which loci are rendered (the cap) and the index order. It does not decide what a
   report reader sees first. **Default:** out of scope here; add a backlog item to show V03 (and
   V13) locus instances in index-rank order in the report. That change would alter
   `report_structure.json`.
5. **Decision bookkeeping.** **Default:** on acceptance, the lead adds a "Lead decisions" line to
   `docs/specs/M9-10-v03-locus-plots.md` saying that decision 3 is replaced and decision 8 refined
   by this spec. No ADR, because §10.2 says nothing about the order across kinds and this rule
   keeps "ordered by q" within each kind. Write an ADR instead if the statistician reads §10.2 as
   one global sort by q.
6. **Numeric ID order** (`prophage_10` before `prophage_2`). **Default:** unchanged string order.
   A natural sort is a separate, cosmetic change.

## Statistician review (2026-10-25)

**Verdict: CHANGES.** The tiered round-robin is the right kind of rule: every comparison in the key
is either an error against its own family's level (tier) or an error against errors of the same
family (rank), so no number from one family is ever compared with a number from another. All ten
examples reproduce under the key as written (checked with an independent script). But the
rotation unit is the *kind*, and one kind, `operon_flag`, holds two §8.1 families. Items 1–3 must
be in the spec before the test step. Items 4–5 are smaller. No re-review is needed if they are
adopted as written below. This is a display rule that adds no new p, q or estimate, so §8.2
(effect sizes with CIs) and calibration gates are unaffected, and attack nulls (§8.5) do not apply.

1. **(blocking) Rotate by §8.1 family, not by kind.** `operon_flags.tsv` mixes two families:
   `uniop_class = OP` rows carry the discordance BH q (§8.1 "Operon discordance", one BH family per
   run, `docs/methods/operons.md` "Family"), and `uniop_class = NO` rows carry the missed-operon
   decoy-FDR q (§8.1 "Missed operon"). `flag_for` emits flags only for OP and NO. Ranking all flags
   by one q list is the same cross-family comparison this spec removes. It is just hidden inside
   one kind. Required:
   - The rotation unit is the **family**. There are five: prophage, bgc, island,
     `operon_discordance` (OP) and `missed_operon` (NO). The key is
     `(tier, family_rank, family_index, locus_id)`, and `family_index` follows the §8.1 row order:
     prophage < bgc < island < operon_discordance < missed_operon. `kind` and every index
     column stay as they are (decision 4). The family is internal, read from `uniop_class`.
   - Operon-flag merging (decision 6) groups by `(family, contig, span)`. Today a merged locus
     takes "its best flag's q" over both families, which is again a cross-family minimum. In
     practice OP pairs take a TU span and NO pairs span two TUs, so cross-family merges should
     not happen. The key makes that guaranteed instead of likely.
   - The edge case "only one kind" becomes "only one family". A pool made only of flags from both
     families now interleaves them; it is no longer today's order.
   - Add example 11: prophage_1 q 0.01; PAIR_1 (OP) q 0.08; PAIR_2 (NO) q 0.01; PAIR_3 (OP)
     q 0.02 → prophage_1, PAIR_3, PAIR_2, PAIR_1 (now: prophage_1, PAIR_2, PAIR_3, PAIR_1).
   - The within-kind property, the metamorphic property and the rotation property all read
     "family" where they say "kind".
2. **(blocking) The rotation property is too weak.** "Counts of any two kinds with loci left differ
   by at most 1" accepts A1 B1 B2 A2 B3 for A = 2 and B = 3 loci. That is not the rule
   (A1 B1 A2 B2 B3), because the property does not fix the order within a round. Replace it with
   the exact round structure. Split each tier block into maximal runs where `family_rank` is
   constant. Run r must contain exactly one locus from each family with at least r loci in that
   tier, in `family_index` order. Compute `family_rank` in the test from the within-family sorted
   subsequence, not from the implementation. Keep the balance property as a readable second
   check if wanted.
3. **(blocking) Reader-facing text must not imply one ranking.** After this change the TSV `rank`
   and the "1., 2., …" in the index figure are a *display order*. A reader who sees "3. PAIR_1 q
   0.09" above "4. prophage_1 q 0.06" (example 3) will read it as an evidence ranking unless told
   otherwise. Required:
   - The index figure prints a header line for each non-empty tier ("significant at each
     family's own level", "not significant", "no estimate"). It also prints one caption line:
     "Order alternates between test families; values in different kinds are on different scales
     and are not compared." The content is not pinned (decision 4), so this is free.
   - The estimate names are per family: `q` for prophage, bgc and operon_discordance,
     `decoy-FDR q` for missed_operon, `FWER p` for island. Item 1 makes the two flag families
     differ, and the acceptance criterion "q for operon_flag" changes to match.
   - The user-guide V03 bullet, the `select_loci` docstring and the `Locus.error` docstring say
     that `rank` is a display order, not a ranking of evidence. They also say that `error` holds
     each family's own estimate, and that §10.2's "ordered by q" applies within each family.
     Nothing may say "ranked by significance" or "most significant first" across kinds.
4. **The metamorphic property needs exact wording.** The map g applied to one family must be
   strictly increasing and must keep the predicate `e ≤ alpha` for every error. So a value at
   alpha may map to alpha or below, never above, and NA stays NA. Add two siblings, both
   order-invariant:
   - (a) change one family's alpha to any value that no locus of that family crosses;
   - (b) a cross-family stress: within a tier, send every error of one family below, then above,
     every error of all other families (for example g(e) = e·10⁻⁶ in tier 0). The test then
     detects any surviving cross-family comparison directly.

   The generators must produce flags of both OP and NO classes, ties within a family, values
   exactly at alpha (and the float neighbours either side via `math.nextafter`) and the same
   `locus_id` in two kinds. Examples 1–10 and the remaining properties (the tiers, input-row
   shuffling, and the order within a family by `(error is None, error, locus_id)`) are the right
   invariants.
5. **(non-blocking) Tier is relative to the `cfg` passed, but the pools are censored at run time.**
   The islands step writes only windows with p ≤ the run's `alpha.context`, and `_flag_loci`
   filters at the `cfg`'s `alpha.operon_flags`. If `select_loci` gets a lower `alpha.context`
   than the run used, islands move to tier 1, which is correct. With a higher one, islands
   between the two levels are absent while regions there are tier 0. This is a known asymmetry
   of re-rendering, not a defect of the rule. Add one edge-case row; no code change.

**Answers.**

- **Q1 (the rule).** Adopt it, with item 1. Kind-first blocks fail the cap (they starve the last
  family). A rotation without tiers interleaves non-significant regions with significant
  islands. A common scale (BH over pooled p, or any p → q conversion) would be a new statistic on
  the plot layer. It would also be invalid: the island p is already family-wise and the flag q
  values come from a decoy procedure. Equal turns per family are a display policy, not an
  inference. That is acceptable as long as item 3 tells the reader so. One side effect: under the
  cap, a family with few loci gets the same share of the first rounds as a family with many.
  This is intended (example 4).
- **Q2 (tier levels).** Use `stats.alpha.context` for islands, not a hard 0.05. Three reasons:
  - §5.2 calls `stats.alpha` "the single source of every significance level", and the §8.1
    "Default" column gives defaults, not constants.
  - The islands step (`docs/methods/islands.md`, "a window is reported when p_w ≤
    `stats.alpha.context`") and `compare/stability.py` use the same key. A hard 0.05 would put
    reported islands in tier 1 whenever a user sets `alpha.context` above 0.05.
  - The "0.05" in §7.3.5(d) step 5 (the Moran's I block-length rule) is an internal tuning
    constant, not the claim level. It is unrelated.

  Flags at `stats.alpha.operon_flags` are right for both flag families. Both §8.1 rows default to
  0.1, and the pool filter uses that key.
- **Q5 (bookkeeping).** No ADR. §10.2's "ordered by q" cannot be a single global sort. Islands
  have no q, so a literal reading sends every island to the NA tail. §8.2's pre-declared families
  also make q values from different families incomparable, so a global sort would have no
  meaning. Reading §10.2 as "ordered by each family's own estimate, within the family" is the only
  coherent reading, and this rule satisfies it. A lead-decisions line in
  `docs/specs/M9-10-v03-locus-plots.md` is enough. It should say that decision 3 is replaced and
  decision 8 refined, and name the family (item 1) as the rotation unit.
