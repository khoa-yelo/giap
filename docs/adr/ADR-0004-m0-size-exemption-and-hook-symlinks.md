# ADR-0004: One-off PR-size exemption for M0, and symlink-safe path protection
- Status: proposed
- Date: 2026-10-14
- Spec sections: §0.2 R3, §14.6, Appendix A.4
- Needs Khoa: yes (a gate exception and a change to a protected file)

## Context
1. **PR size.** The guard in `ci.yml` fails a PR above 800 changed lines (`ci/gates.yaml`,
   `pr_lines_fail`), excluding lock files and fixtures. The M0 scaffold PR (#1) changes about
   4,800 lines plus the 3,262-line `GIAP_SPEC.md`: the spec itself, the Appendix A templates, CI,
   the probe, the full CLI option surface of §5.1, the config model of §5.2 and the simulator.
   Splitting it into 6+ PRs would add review overhead without adding review value, because most
   of it is verbatim spec content.
2. **Hook bypass.** `.claude/hooks/protect_paths.py` compares paths with `os.path.abspath`, which
   does not resolve symlinks. On FarmShare `/farmshare/home/classes` is a symlink to
   `/home/classes`, so an Edit addressed through the `/farmshare/...` path to `ci/gates.yaml` is not
   blocked (reproduced by the code reviewer). The hook is protected (`.claude/hooks/*`), so the
   agent may not change it on its own.

## Decision (proposed)
1. PR #1 (M0 scaffold) may merge with the PR-size guard failing. No gate value changes; later
   PRs follow the 400/800-line rule.
2. Change two lines in `.claude/hooks/protect_paths.py`:
   ```python
   root = os.path.realpath(os.environ.get("CLAUDE_PROJECT_DIR", os.getcwd()))
   ...
   rel = os.path.relpath(os.path.realpath(path), root)
   ```
   Merged with the trailer `Gate-Change-Approved: ADR-0004` once Khoa sets this ADR to
   `accepted-by-khoa`.

## Evidence (before → after: metrics, tests, simulations)
- PR #1 guard log: `changed lines: 8058 (warn 400, fail 800)` (includes the spec).
- Hook today (rc 2 = blocked): project root and file path spelled the same way → 2; spelled
  differently (`/home/classes/...` vs `/farmshare/home/classes/...`, either direction) → 0.
  The session's `$PWD` is the `/farmshare/...` spelling, so the mismatch is easy to hit.
  With `realpath` on both, all four combinations give 2.

## Consequences
CI remains the backstop for protected files in any case (`scripts/ci/check_protected.py`).
