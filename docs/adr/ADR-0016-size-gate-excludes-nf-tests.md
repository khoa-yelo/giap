# ADR-0016: The PR-size gate does not count workflow/tests/ (amends ADR-0007)
- Status: accepted-by-khoa (2026-10-14, in session: "Exclude workflow/tests/")
- Date: 2026-10-14
- Spec sections: R-size rule (§14.6), ADR-0007, ci/gates.yaml (values unchanged)
- Needs Khoa: yes (changes what an accepted gate counts) — answered

## Context
ADR-0007 gates PR size on code only (src/, workflow/, scripts/). nf-test specs and their data must live
under workflow/tests/, so they were counted as code: M2-07 is 762 code lines but 1,228 by the gate. The
role guard already treats workflow/tests/ as tests (PR #40).

## Decision
The size gate counts src/, workflow/ (excluding workflow/tests/) and scripts/. workflow/tests/ is reported
with the other non-gated lines. pr_lines_warn / pr_lines_fail in ci/gates.yaml are unchanged.

## Consequences
One CI step changes (.github/workflows/ci.yml, "PR size"). Code under workflow/ (modules, subworkflows,
configs, main.nf) still counts.
