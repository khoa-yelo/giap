# M9-23: V16 `compare` and `report/compare.html`

Milestone M9 · **pass 2** · depends on: M9-22, M9-14 · estimate ≈ 320 code lines (+ ≈ 260 test
lines) · not statistical

## Spec sections implemented

- §10.2 V16: "UpSet of detected genes across runs; pairwise abundance concordance matrix
  (Spearman ρ, scatter)"; library "UpSetPlot (fallback: custom matplotlib)"; "Plotly
  (concordance)". "V16 is written only by `giap compare`."
- §7.9.3 (`plots/V16_compare.*`, `report/compare.html`), §10.1, §10.3 (self-contained HTML rules).

## Why custom matplotlib

UpSetPlot 0.9.0 fails on the locked pandas 3.0.6 (`ValueError: Invalid RGBA argument: nan` from a
Copy-on-Write `fillna(inplace=True)`, reproduced in the session). §10.2's fallback applies; no
`import upsetplot` anywhere (index decision 29).

## Panels (double width)

- **(a) UpSet:** from `compare_detection` intersection rows: bars of `n_genes` per exact
  intersection (sorted by `n_genes` descending, then name; `none` excluded and printed as a
  footnote), a dot matrix below (filled `#0072B2` dot = run in the intersection, joined by a line;
  empty dots `#BBBBBB`), set-size bars at the left (detected genes per run). More than 30
  intersections → the 30 largest drawn, the rest summed as "other" in the TSV.
- **(b) Concordance matrix:** runs × runs; lower triangle cells coloured by ρ (cividis, 0 to 1,
  negative values clipped with a marker), text "ρ `<r>` [`<lo>`, `<hi>`]"; diagonal labels; upper
  triangle: small scatters of shared measured log2 (a vs b) with the identity line, binned per
  M9-01 when > 20,000 points.

Data TSV `V16_compare.tsv`: `panel` enum{a, b}; (a) `intersection`, `n_genes`, `runs` (list); (b)
`run_a`, `run_b`, `spearman_rho`, `ci_low`, `ci_high`, `n_shared`, and scatter rows with
`locus_tag`, `x`, `y`. Sorted by panel and the drawing order.

Title: "`<n>` genes detected in all `<k>` runs; median pairwise ρ `<m>`" (counts and ρ, each ρ with
its CI in the figure).

## `report/compare.html`

Same template machinery (M9-02, M9-14): sections `runs` (labels, directories, genome SHA-256,
detection tier per run), `detection` (V16 a, the intersection table), `abundance` (V16 b, Plotly
twin of the concordance scatter, `compare_abundance` table with BA columns), `stability`
(`compare_stability` table), `provenance`. Self-contained; same offline and link checks.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| two runs | a 2 × 2 matrix; UpSet with ≤ 3 intersections |
| ρ NA (`too_few_shared`) | grey cell "ρ NA (too few shared)" |
| all intersections empty (no detections anywhere) | empty state "No gene detected in any run" |
| 12 runs (2¹² intersections possible) | the 30 largest drawn |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | M9-22 example 1 | (a) bars 1, 1, 1, 1 (`none` in footnote); set sizes A1 3, A2 3, A3 1; title "1 gene detected in all 3 runs; …" |
| 2 | ρ(A1, A2) 0.8 [0.6, 0.93] | cell text "ρ 0.8 [0.6, 0.93]" |
| 3 | 40 intersections | 30 drawn; TSV has an `other` row with the remainder's sum |
| 4 | `compare.html` | `external_loads` = []; relative links resolve |

## Acceptance criteria

- [ ] `tests/unit/viz/test_v16_compare.py`: examples 1–3; no `upsetplot` import in `src/`
  (an AST test).
- [ ] `tests/unit/report/test_compare_html.py`: example 4; section ids.
- [ ] M9-25 snapshots V16 on a three-run fixture.

## Questions for the lead

1. Wording "genes" for n = 1 ("1 genes") — **Proposed default:** a tiny pluralisation helper in
   `giap.viz.base` used by every title.

## Lead decisions (2026-10-23)
1. **Question 1 accepted:** a small pluralisation helper in `giap.viz.base` (for example
   `plural(n, "gene")` → "1 gene" / "2 genes"). Every title uses it. Existing titles switch to it
   only when their tests are next touched; that is a backlog item, not this PR.
2. #232 (M9-22 `giap compare`) has merged, so the test-author starts from main.
