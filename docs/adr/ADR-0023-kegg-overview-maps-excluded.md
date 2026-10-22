# ADR-0023: KEGG global and overview maps are excluded from the KEGG collection
- Status: accepted
- Date: 2026-10-21
- Spec sections: §7.7.1, §7.7.4
- Needs Khoa: no. This fills a gap in §7.7.4, which says "KEGG pathways" without saying which, and it changes no default, threshold or interface. It is the M8 index's ADR-C, lead-level once the lead accepted the M8-03 default (lead decision 1).

## Context
KEGG's "Global and overview maps" are unions of many metabolic pathways. Some examples:
- map01100, Metabolic pathways
- map01110, Biosynthesis of secondary metabolites
- map01200, Carbon metabolism
- map01310, Nitrogen cycle
- map01320, Sulfur cycle

They overlap nearly every metabolic set. Where they survive the 500-gene size cap, they dominate the BH family with redundant, uninterpretable hits. In the KEGG link file, every `path:map` row is also duplicated as a `path:ko` row.

## Decision
The KEGG collection keeps only `path:map` links. It drops these maps:
- `map011xx`;
- `map012xx`;
- map01310;
- map01320.

These are KEGG's "Global and overview maps" group as of the installed release. A KO whose only links are to these maps belongs to no KEGG set. The exclusion is recorded here and in docs/methods/function-collections.md. Pass 1 emits no per-run count, so spec M8-03's "excluded and counted" becomes "excluded and documented".

## Evidence
The code review of #159 found that map01310 and map01320 (with 39 and 53 KO links) passed the original `map011xx`/`map012xx` rule. tests/unit/function/test_collections_review.py pins the full list, along with neighbouring maps that are kept.

## Consequences
- A KEGG release that adds a new overview map must update the pattern (`giap.function.collections._OVERVIEW`) and this ADR.
- A future `kegg_overview_excluded: N` run flag would be an additive change.
