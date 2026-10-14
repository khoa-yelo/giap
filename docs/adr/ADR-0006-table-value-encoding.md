# ADR-0006: Encoding of special values in GIAP tables, and strict BED writing
- Status: accepted
- Date: 2026-10-14
- Spec sections: §6.6, §9.2, §9.5
- Needs Khoa: no (fills a gap; no default, threshold or interface changes)

## Context
§6.6 fixes `%.6g` for floats but does not say how infinities, NaN or list cells are written.
Tables are validated with Frictionless (§9.5), whose default number type reads `INF`, `-INF`
and `NaN`, not Python's lowercase `inf`. The M1-04 review also found that malformed BED12 blocks
and thick ranges were silently clamped.

## Decision
1. Floats: finite values `%.6g`; `+inf` → `INF`, `-inf` → `-INF`; NaN and missing → `NA`.
2. List cells are joined with `;`; an item containing `;` or a missing item is an error.
   proBed's `modifications` field uses `,` as its format requires.
3. BED/proBed writers validate instead of clamping: BED12 blocks must start at 0, be sorted and
   non-overlapping, and end at the record length; the thick range must lie inside the interval;
   names may not contain whitespace; scores are integers 0–1000 (bool rejected). Asking for BED6
   output while passing thick or colour fields is an error.

## Evidence (before → after: metrics, tests, simulations)
Review of M1-04 (manual mutations survived; clamping reproduced on 100–200 → 110–160).

## Consequences
Readers map `INF`/`-INF` back to ±inf and `NA` to missing.
