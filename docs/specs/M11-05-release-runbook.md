# M11-05: the release runbook (`docs/release.md`)

Milestone M11 · **pass 1 (MVP), last issue** · depends on: M11-04 · estimate 0 code lines +
≈ 40 test lines + ≈ 200 prose lines · not `stats` · roles: docs (`docs/release.md`), test-author
(`tests/unit/test_m11_05_runbook.py`) · Khoa: performs the steps marked **Khoa**.

## Spec sections implemented

- §3.10 "Never push tags"; "Tags and releases come from `release.yml` (workflow_dispatch) or from
  Khoa". §12 `release.yml` row and "Bioconda recipe PR as a manual follow-up". §14.8 items 5–7
  ("Each release"; "Licence and account steps only Khoa can take"; "making the repository public").
- R10 (the agent decides everything else), R12 (no tokens in the runbook).

## Definition

`docs/release.md`, linked from the site's Developer guide (M11-07). Numbered steps; each says
**who** (agent role / lead / Khoa), the exact command or web page, and the evidence to record in
PROGRESS. No step tells an agent to tag, upload, or dispatch with `dry_run=false`.

### A. Once, before the first release (Khoa)

1. PyPI → Account → Publishing → "Add a new pending publisher": project `gippa`, owner
   `khoa-yelo`, repository `gippa`, workflow `release.yml`, environment `pypi`. Same on TestPyPI with
   environment `testpypi` (Q-K1).
2. GitHub → Settings → Environments: create `pypi` and `testpypi`, required reviewer `khoa-yelo`,
   deployment branch `main` only (Q-K2).
3. Zenodo → log in with GitHub → GitHub tab → switch on `khoa-yelo/gippa` (Q-K3).
4. GitHub → Settings → Pages → Source "GitHub Actions" (Q-K9; needed by M11-07's deploy job).
5. GHCR: after the first `containers.yml` run, each `gippa-<tool>` package → Package settings →
   visibility Public (Q-K8; never `gippa-uniop`, which is never pushed).

### B. Each release: agent and lead

6. Release-prep PR (`m11/release-<v>`, M11-13): `__version__ = "<v>"`, `towncrier build --version
   <v> --date <YYYY-MM-DD> --yes`, `CITATION.cff` `version` + `date-released`, known-limitations
   fragment updated before the build. Merged by the lead after CI.
7. Evidence on the merged commit R (Slurm, R14 caps): `calibration_full.sbatch`, the nightly
   calibration subset, `weekly.sbatch` for TE8/e2e, `test_branch.sbatch` full suite; job IDs in
   PROGRESS.
8. `readiness.py compute --version <v> --out validation/release/<v>/` on FarmShare; commit the two
   files in a paperwork-only PR (`docs(release): readiness <v>`); merge.
9. Lead: `gh workflow run release.yml -f version=<v> -f dry_run=true`; record the run URL.
10. Lead: the review packet (`local/review/m11-release/`: READINESS.md, the dry-run summary, the
    notes, artefact sha256s), PR or issue labelled `needs-khoa`, message to Khoa.

### C. Each release: Khoa

11. Approve (§14.8 item 5) in writing; the lead copies the words into PROGRESS.
12. Rehearsal (first release; optional later): Actions → release → Run workflow with
    `version=<v> dry_run=false target=testpypi`; approve the `testpypi` environment; check
    `pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ gippa==<v>`.
13. Real: `version=<v> dry_run=false target=pypi`; approve the `pypi` environment. The workflow
    uploads, tags `v<v>` and creates the GitHub release; Zenodo mints the DOI.

### D. After the release (agent, then Khoa)

14. Agent: PyPI install smoke in a fresh venv on FarmShare (inside the allocation):
    `gippa --version`, `gippa version --json`, `gippa preflight` on the tiny genome.
15. Agent: post-release PR: `CITATION.cff › identifiers` (Zenodo concept DOI and version DOI from
    the Zenodo record page), README badge, `__version__` stays `<v>` until the next release.
16. Agent: render the Bioconda recipe (M11-10) with the PyPI sdist sha256 and hand it to Khoa;
    **Khoa** opens the `bioconda/bioconda-recipes` PR from their fork (needs their account).
17. Agent: PROGRESS › Done "Release <v>" with the PyPI URL, tag, release URL, DOI, Bioconda PR.

### E. If something fails

- PyPI upload done, release step failed → re-dispatch the same inputs (M11-04 edge case).
- Wrong artefact on PyPI → never delete; yank (Khoa, PyPI web UI) and release `<v+patch>`.
- Readiness `fail` → no publish; fix through the normal loop or a waiver ADR (Khoa).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | a reader following B–C for 1.0.0 | every command is copy-pasteable with `<v>` = `1.0.0`; no command contains `git tag` or `git push --tags` |
| 2 | test scans `docs/release.md` | every step header names exactly one actor in **bold** (`agent`, `lead`, `Khoa`) |

## Edge cases

- **Step 12 uses the real version**: a TestPyPI upload of `<v>` occupies `<v>` on TestPyPI only; PyPI is a
  separate index, so the real upload is unaffected.
- **Pages or Zenodo not enabled** (Khoa skipped A3/A4): the release still works; the DOI or the site
  is missing and the post-release PR waits; recorded under PROGRESS › Questions for Khoa.
- Genomic edge cases: not applicable.

## Acceptance criteria

1. `docs/release.md` exists with sections A–E and steps 1–17.
2. A test asserts: no `git tag`, `git push --tags`, `--follow-tags` or `dry_run=false` in any step
   whose actor is not Khoa; every step names an actor; the workflow, environment and project names
   match `release.yml` (parsed) — `release.yml`, `pypi`, `testpypi`, `gippa`.
3. AGENTS.md gains one line (lead edit, not in this PR's role): "Agents never dispatch release.yml
   with dry_run=false and never push tags (docs/release.md)."

## Estimate

0 code; ≈ 200 prose; tests ≈ 40.
