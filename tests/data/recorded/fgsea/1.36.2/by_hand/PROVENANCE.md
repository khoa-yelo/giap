# by_hand: recorded-by-hand (not R output)

- Written by the test author on 2026-10-22 without R (the session cannot run R).
- Format: docs/specs/M8-17-r-fgsea-engine.md §3 (`results.tsv` header `job_id pathway size ES NES
  pval log2err leadingEdge`; numbers as R's `format(x, digits = 17)`; `NA`; empty leading edge).
- Content: 12 TUs (one tie at 1.8), jobs `main:cog` and `main:kegg` with sizes [3, 5]; `K00001`
  is an unbalanced-signs row (`pval`, `NES`, `log2err` NA), `K00002` has `log2err` NA with a
  finite tiny `pval`, `K00003` has ES 0 and an empty leading edge.
- The statistics are invented and are **not** what fgsea 1.36.2 would return.
