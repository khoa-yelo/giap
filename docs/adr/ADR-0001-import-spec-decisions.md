# ADR-0001: Import the settled decisions of GIPPA_SPEC.md §2.1 and §2.4
- Status: accepted-by-khoa
- Date: 2026-10-14
- Spec sections: §2.1, §2.4
- Needs Khoa: no (the decisions were accepted when the spec was written)

## Context
GIPPA_SPEC.md v1.1 (FarmShare edition) records 20 design decisions (§2.1) and the implementation
choices G1–G19 (§2.4), all accepted by Khoa before development started.

## Decision
All of §2.1 (decisions 1–20) and §2.4 (G1–G19) are in force as written. Later ADRs may supersede
individual items; each such ADR names the item it replaces.

## Evidence (before → after: metrics, tests, simulations)
None needed: this ADR imports decisions rather than changing them.

## Consequences
Any change to a default, threshold, interface or dependency fixed by §2 needs a new ADR that
Khoa accepts (R10, §14.8).
