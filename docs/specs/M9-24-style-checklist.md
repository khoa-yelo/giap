# M9-24: The §14.4 style checklist as code, and accessibility

Milestone M9 · **pass 2** · depends on: M9-01 (and each plot it checks) · estimate ≈ 300 code lines
(+ ≈ 260 test lines) · not statistical

## Spec sections implemented

- §14.4 visualisation loop ("Metric: style-checklist score, then Khoa's preference … Stop rule: Khoa
  approves and the snapshot is locked") and the checklist: "the title states the finding; axes carry
  units; palette as in §10.1; no overlapping labels; legible at single-column width; data TSV
  written; colourblind check (simulated deuteranopia); ≤ 20,000 drawn elements".
- §10.1 (palettes, widths, 20,000 per track), §1 rule 6.

## Public interface

`src/giap/viz/checklist.py`
```python
@dataclass(frozen=True)
class CheckResult:
    plot: str; item: str; passed: bool | None; detail: str   # None = recorded, not gated
def check_instance(inst: Instance, pid: str, files: Sequence[Path]) -> list[CheckResult]
def score(results: Sequence[CheckResult]) -> float          # passed / gated items
```
`scripts/dev/style_checklist.py OUT` prints a table per plot and writes
`local/review/<name>/checklist.tsv` (not committed).

## The eight items

| Item | Check | Gated |
|---|---|---|
| title states the finding | title non-empty, contains a digit or one of the empty-state phrases, and no "express" | yes |
| axes carry units | every visible axis label with a quantity has a unit in parentheses or is in a declared unitless list (`log2 abundance`, `q`, `lfdr`, `ρ`, `NES`, counts) | yes |
| palette as §10.1 | every categorical face/edge colour drawn ∈ §10.1 table ∪ Okabe–Ito ∪ {black, white, `#BBBBBB`}; continuous mappables use viridis or cividis; COG and detection palettes never in one figure (M9-01 Q4) | yes |
| no overlapping labels | pairwise intersection of the rendered bounding boxes of all visible text artists (tick labels, annotations, titles) on the Agg renderer is empty | yes |
| legible at single-column width | every text ≥ 6 pt at the saved figure size (index decision 31); figure width ∈ {3.5, 7.2} in | yes |
| data TSV written | `data/Vnn_<name>[_suffix].tsv` exists, and has rows unless the instance is an empty state | yes |
| colour-blind check | the PNG converted with the Machado et al. 2009 deuteranopia matrix (severity 1.0, applied in linear RGB); the minimum pairwise CIE76 ΔE between the plot's categorical colours, before and after, is reported; the simulated PNG is written next to the original under `local/review/` | **recorded** (index decision 30) |
| ≤ 20,000 drawn elements | per axes: patches + lines + collection offsets + texts ≤ 20,000 | yes |

## Degenerate cases

| Situation | Behaviour |
|---|---|
| empty-state instance | units, palette and element items `passed = None` ("not applicable") |
| text rotated 90° | its rotated bounding box is used |
| SVG-only formats | colour-blind check renders an in-memory PNG at 100 dpi |
| V1 polar axes (pyCirclize) | element count per track axes; overlap check on sector labels |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | a figure with tick labels "1", "10" overlapping | `no overlapping labels` fails with both strings in `detail` |
| 2 | an axis labelled "length" | units item fails ("length" needs a unit) |
| 3 | a bar in `#FF0000` | palette item fails naming `#FF0000` |
| 4 | a 5 pt annotation | legibility fails |
| 5 | the §10.1 detection colours under deuteranopia | ΔE values recorded (no pass/fail) |
| 6 | 25,000 points on one axes | element item fails |

## Acceptance criteria

- [ ] `tests/unit/viz/test_checklist.py`: examples 1–6 on crafted figures; the Machado matrix
  constants equal the published values (Machado, Oliveira & Fernandes 2009, IEEE TVCG 15(6), the
  deuteranomaly matrix at severity 1.0) and are cited in the module docstring.
- [ ] `tests/snapshot/test_checklist_all.py`: every registered plot on the snapshot fixture passes
  every gated item (this is what M9-25's approval packet reports as the score).

## Questions for the lead

1. A gate on the deuteranopia ΔE needs a threshold; the spec gives none. **Proposed default:**
   record only; Khoa decides at M9-25 (index decision 30).
2. Minimum 6 pt (index decision 31). **Proposed default:** 6 pt; a stricter value can follow Khoa's
   review.
3. "Title states the finding" is checked only mechanically (a number or a declared phrase).
   **Proposed default:** mechanical check plus Khoa's review; no language model judges titles.
