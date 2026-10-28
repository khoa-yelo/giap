# M11-13: the M11 gate — release candidate, readiness on the final commit, Khoa's approval

Milestone M11 · **pass 2, the gate** · depends on: M11-01…M11-12 (M11-09–M11-11 soft), the M1–M10
gates recorded in PROGRESS (Q-K12) · estimate ≈ 30 code lines (the version line and towncrier
output) + ≈ 60 test lines, as three PRs (13a release prep, 13b readiness evidence, 13c post-release)
· not `stats` · roles: docs (13a changelog build, known-limitations fragment), validation (13b
evidence runs and `readiness.py compute`), lead (review packet, PROGRESS gate block, Khoa
hand-off, transcribing Khoa's answer), implementer (13c CITATION DOI, README badge) · **Khoa: yes,
§14.8 item 5.** No agent, the lead included, marks this gate done; only Khoa's dispatch publishes.

## Spec sections implemented

- §13 row M11 gate: "the acceptance criteria below; **Khoa approves** the release"; §13.2.
- §3.10 "never push tags … Tags and releases come from `release.yml` (workflow_dispatch) or from
  Khoa"; §12 `release.yml`; §12 Releases "a Zenodo DOI each time"; §11 "a real E. coli run gates
  every release" and e2e "weekly and before release".
- §14.8 item 5; R3 (no threshold moves to make a row pass; waivers are Khoa's ADRs); R8.

## Steps (the runbook of M11-05, section B–D, applied to 1.0.0)

| Step | PR / action | Who | Output |
|---|---|---|---|
| 1 | **13a** `chore(release): 1.0.0`: `__version__ = "1.0.0"`; `+known-limitations.misc.md` updated with the open non-gate issues (Q-K12); `towncrier build --version 1.0.0 --date <merge day> --yes` (fragments removed, `CHANGELOG.md` gains `## GIPPA 1.0.0 (<date>)`); `CITATION.cff` `version: 1.0.0`, `date-released` | docs + lead | merged commit **R** |
| 2 | evidence on R: `calibration_full.sbatch`, nightly calibration subset + attacks, `weekly.sbatch` (K-12, Nissle, TE8), `test_branch.sbatch` full suite; R14 CPU sum before each sbatch | validation | job IDs in PROGRESS; rows in `nightly.tsv`/`weekly.tsv` for R |
| 3 | **13b** `docs(release): readiness 1.0.0`: `readiness.py compute --version 1.0.0 --freeze-run <A2 freeze dir> --out validation/release/1.0.0/` | validation | `readiness.json`, `READINESS.md`; paperwork-only diff |
| 4 | `gh workflow run release.yml -f version=1.0.0 -f dry_run=true` on main | lead | green run URL; artefact sha256s |
| 5 | review packet `local/review/m11-release/r<k>/`: READINESS.md, the dry-run summary, `notes.md`, artefact sha256s, the open-issue list; PROGRESS › Questions for Khoa row "Approve release 1.0.0"; the 13b PR labelled `needs-khoa` | lead | Khoa notified in chat |
| 6 | approval in writing; if a row needs a waiver, Khoa accepts the waiver ADR first and step 3 reruns | **Khoa** | PROGRESS transcription |
| 7 | rehearsal `dry_run=false target=testpypi`, then `target=pypi` | **Khoa** | PyPI release, tag `v1.0.0`, GitHub release, Zenodo DOI |
| 8 | PyPI install smoke in a fresh venv inside the allocation | validation | PROGRESS line |
| 9 | **13c** `docs(release): DOI and badges for 1.0.0`: `CITATION.cff › identifiers` (concept + version DOI), README badges; Bioconda recipe rendered (M11-10) and handed to Khoa | implementer + docs | merged; Khoa opens the Bioconda PR |
| 10 | PROGRESS › Done `### Gate M11 (Release 1.0.0, <date>)`: links to PyPI, release, DOI, readiness, Khoa's words; `progress_bars.py --write` | lead | gate record |

## Examples

| # | Situation | Expected |
|---|---|---|
| 1 | a docs-only PR (`docs/user-guide/run.md`) merges after R, before the dispatch | readiness still names R; `readiness.py verify` passes by M11-03's paperwork allow-list; no rerun |
| 2 | a `fix:` PR to `src/` merges after R | readiness becomes stale (`fail (stale: src/…)`); steps 1–3 repeat on the new commit (no new version number unless 13a must change; the `CHANGELOG.md` section is edited in the fix PR) |
| 3 | R7 fails (k12 standard 58 min) | no dispatch; the lead proposes either a speed-loop issue or a waiver ADR for Khoa (Q-K6) |
| 4 | Khoa answers "approved" | step 7 is Khoa's; agents only watch the run and do 8–10 |
| 5 | PyPI publish succeeded, `github-release` failed | Khoa re-dispatches the same inputs (M11-04 edge case); no new tag by any agent |

## Edge cases

- **Main moves during the gate** (other agents merging): the lead pauses merges touching `src/`,
  `workflow/`, `tests/calibration/` between step 1 and step 7 (a PROGRESS note "release freeze");
  docs and paperwork PRs may merge.
- **Readiness `pending` rows that can never be filled** (a TE whose data are unavailable): those are
  `deferred` with the reason by M11-03's rules, not `pending`.
- **Tag exists** (someone created `v1.0.0` by hand): `release.yml` refuses; ask Khoa.
- **Zenodo did not archive** (integration off): the DOI PR waits; PROGRESS › Questions for Khoa.
- Genomic edge cases: not applicable.

## Acceptance criteria (the gate)

1. `validation/release/1.0.0/readiness.json` on main has `overall: pass` (waived rows cite
   accepted-by-khoa ADRs) and `readiness.py verify` passes at the commit Khoa dispatched.
2. A green `dry_run=true` run on that commit precedes Khoa's non-dry dispatch.
3. Khoa's written approval is transcribed in PROGRESS before the dispatch.
4. After Khoa's dispatch: PyPI has `gippa 1.0.0` (wheel + sdist with the dry run's sha256s), tag
   `v1.0.0` points at the dispatched commit, the GitHub release has the notes, and Zenodo has a DOI;
   `CITATION.cff` carries the DOI (13c).
5. No agent commit, command log or PR in M11 creates or pushes a tag (the lead checks
   `git log --all --format='%H %s'` and the Actions history for tag events by non-Khoa actors).

## Estimate

Code ≈ 30; tests ≈ 60 (a test that the 13a state is consistent: version, CFF, changelog heading,
empty `changelog.d/`). Evidence runs: one full calibration job array, one weekly job, one full-suite
job; wall time is queue-bound.

## Questions for the lead (reversible defaults in use)

- Q-K12 (Khoa) open non-gate issues listed as known limitations. Q-K4 version. Q-K5/Q-K6 criteria.
- Release freeze between steps 1 and 7. Default: yes, for `src/`, `workflow/`, `tests/calibration/`.
