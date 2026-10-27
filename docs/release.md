# Releasing

How a GIAP release is prepared, approved and published. Agents prepare everything and run dry
runs. Khoa is the only one who approves, publishes and sets up accounts.

!!! danger "Two rules for every agent"
    Agents never push tags, and never dispatch `release.yml` with `dry_run=false`.
    Tags and GitHub releases come only from the `github-release` job of `release.yml`, run by
    Khoa's dispatch (GIAP_SPEC.md §3.10). An agent that finds a step it cannot do without breaking
    one of these rules stops and asks Khoa through PROGRESS › Questions for Khoa.

Each numbered step names one actor in bold:

- **agent**: a subagent in its role (docs, validation, implementer).
- **lead**: the main session.
- **Khoa**: the account owner.

Each step also gives the exact command or web page, and the evidence to record in `PROGRESS.md`.
`<v>` is the version (PEP 440, matching `^\d+\.\d+\.\d+(rc\d+)?$`, e.g. `1.0.0` or `1.0.0rc1`).
`<base>` is `<v>` without its `rc` suffix (`1.0.0` for both `1.0.0rc1` and `1.0.0`). No
step needs a token in a file or a command line (R12): trusted publishing uses OIDC, and the
workflow's `GITHUB_TOKEN` is the only credential it uses.

The workflow inputs, from `.github/workflows/release.yml`:

| input | type | default | meaning |
|---|---|---|---|
| `version` | string, required | none | must equal `giap.__version__` on the dispatched commit |
| `dry_run` | boolean | `true` | verify, build, smoke and notes only: no upload, no tag |
| `target` | choice `testpypi`, `pypi` | `testpypi` | the index for an upload when `dry_run` is false |

The workflow refuses to run in these cases:

- it was dispatched on any ref other than `refs/heads/main`;
- the version does not match the pattern above;
- the version input differs from `giap.__version__`;
- tag `v<v>` already exists.

The `publish` job needs three things: `dry_run` false, `github.actor == 'khoa-yelo'`, and approval
of the `target` environment. The `github-release` job also needs `target == 'pypi'`. A TestPyPI
upload therefore never creates a tag.

## A. Once, before the first release

### 1. PyPI and TestPyPI pending publishers (**Khoa**)

- PyPI: <https://pypi.org/manage/account/publishing/> → "Add a new pending publisher" → GitHub.
  Enter PyPI project name `giap`, owner `khoa-yelo`, repository `giap`, workflow name
  `release.yml`, environment name `pypi`.
- TestPyPI: <https://test.pypi.org/manage/account/publishing/> → the same form. Enter PyPI
  project name `giap`, owner `khoa-yelo`, repository `giap`, workflow name `release.yml`,
  environment name `testpypi`.
- Evidence: PROGRESS line "Q-K1 done: pending publishers on PyPI and TestPyPI (date)".

### 2. GitHub environments `pypi` and `testpypi` (**Khoa**)

- GitHub → `khoa-yelo/giap` → Settings → Environments → New environment. Create `pypi`, then
  `testpypi`. For each one:
    - Deployment protection rules → Required reviewers → `khoa-yelo`.
    - Leave "Prevent self-review" off: Khoa both dispatches and approves.
    - Deployment branches and tags → Selected branches → `main`.
- Add no environment secrets. Trusted publishing needs none.
- Evidence: PROGRESS line "Q-K2 done: environments pypi/testpypi, reviewer khoa-yelo, main only".

### 3. Zenodo GitHub integration (**Khoa**)

- <https://zenodo.org> → log in with GitHub → account menu → GitHub → switch on `khoa-yelo/giap`.
  Zenodo lists only public repositories, so this needs §14.8 item 7 (making the repository public)
  done first. Zenodo reads its metadata from `CITATION.cff`. There is no `.zenodo.json`.
- From then on, every published GitHub release (step 13) is archived and gets a DOI.
- If this step is skipped, the release still works. The DOI is missing, the post-release PR (step
  15) waits, and the lead records it under PROGRESS › Questions for Khoa.
- Evidence: PROGRESS line "Q-K3 done: Zenodo integration on".

### 4. GitHub Pages (**Khoa**)

- GitHub → Settings → Pages → Build and deployment → Source "GitHub Actions" (Q-K9).
- GitHub → Settings → Secrets and variables → Actions → Variables → New repository variable
  `PAGES_ENABLED` = `true`. The `deploy` job of `docs.yml` runs only on `main` and only when
  `vars.PAGES_ENABLED == 'true'`. Until then, `docs.yml` only builds the site
  (`mkdocs build --strict`) and uploads the Pages artefact.
- Pages on a private repository needs a paid GitHub plan, so on a free plan this step also waits
  for §14.8 item 7 (making the repository public).
- GitHub creates the `github-pages` environment itself. Leave it as GitHub sets it up.
- Evidence: PROGRESS line "Q-K9 done", and the URL of the first green `deploy` job.

### 5. GHCR package visibility (**Khoa**)

- Do this after all five packages exist: `ghcr.io/khoa-yelo/giap-<tool>` for `promotech`,
  `promoter_calculator`, `terminatornet`, `bactermfinder` and `ostir`. A `containers.yml` run
  started by a push builds only the tools whose `workflow/containers/<tool>/` or
  `workflow/envs/<tool>.yml` changed. To build all five, the lead dispatches
  `gh workflow run containers.yml -f tools=all`. Packages are private on first push.
- For each of the five: GitHub → profile → Packages → `giap-<tool>` → Package settings → Danger
  Zone → Change visibility → Public (Q-K8).
- There is never a `giap-uniop` package. UniOP is unlicensed and `containers.yml` refuses to build
  it (R5, R11).
- Until a package is public, Apptainer pulls it on FarmShare with 401 or 403, the M11-06 smoke
  reports `private package`, and the digest PR stays a draft.
- Evidence: PROGRESS line "Q-K8 done: five GHCR packages public".

## B. Each release: agent and lead

### 6. Release-prep PR `m11/release-<v>` (**agent**)

- Role: docs. The lead merges the PR after CI. Changes on the branch:
    - `src/giap/__init__.py`: `__version__ = "<v>"`.
    - `changelog.d/+known-limitations.misc.md`: update it before the build. It lists the open
      non-gate issues (Q-K12). It must also carry the licence limit "geNOMAD is for non-commercial
      use only", which the fragment does not have today.
    - Run `towncrier build --version <base> --date <YYYY-MM-DD> --yes`. This turns the fragments
      into `## GIAP <base> (<YYYY-MM-DD>)` in `CHANGELOG.md` and removes them.
    - `CITATION.cff`: set `version: <v>` and `date-released: <YYYY-MM-DD>`.
- **The rc path**: the towncrier `--version` is always `<base>`, never `<v>` with `rc`.
    - `notes.py 1.0.0rc1` reads the `## GIAP 1.0.0 (` section and adds "Release candidate
      1.0.0rc1 of GIAP 1.0.0." in front of it. A `## GIAP 1.0.0rc1 (…)` heading is never read, and
      the build job would exit 1.
    - So one `## GIAP 1.0.0 (…)` section serves `1.0.0rc1` and the final `1.0.0`.
    - The final `1.0.0` prep PR after an rc changes only `__version__` and `CITATION.cff`
      (`version`, `date-released`). It builds no changelog: the fragments are already consumed,
      and the section stays as the rc PR wrote it.
    - That PR changes `src/`, so it is a new R. Steps 7–8 run again and overwrite
      `validation/release/1.0.0/`.
- `scripts/ci/check_changelog.py` exempts this PR only through its commit subjects: every
  non-merge commit has to be of type `chore`, e.g. `chore(release): <v>`. If the check fails
  anyway, add the `no-changelog` label and push again. Labels are read from the push event, and a
  re-run reuses the old labels.
- The merged commit is **R**. **Release freeze** from R until step 13, recorded in PROGRESS:
    - Merge only PRs whose every changed path is paperwork in `readiness.py`'s sense:
      `PROGRESS.md`, `README.md`, `CITATION.cff`, `CHANGELOG.md`, or a path under
      `validation/release/`, `docs/` (but not `docs/methods/`) or `changelog.d/`.
    - Any other path, such as `src/`, `workflow/`, `tests/`, `scripts/`, `.github/`, `ci/`,
      `benchmarks/`, `pyproject.toml` or `uv.lock`, makes `readiness.py verify` fail with
      `fail (stale: <path>)`.
    - `release.yml` builds main's HEAD, so such a merge makes the PyPI publish fail at verify.
- Evidence: PR URL and the SHA of R.

### 7. Evidence runs on R (**agent**)

- Role: validation, on FarmShare. Fast-forward the main checkout to R first. Source
  `scripts/farmshare_env.sh`. Before each `sbatch`, add up the running and pending CPUs of the
  `giap-*` jobs (R14: at most 128 CPUs, 4 GPUs and 20 jobs).
    - `bash scripts/slurm/calibration_full_launch.sh origin/main`: pins a worktree at R, checks the
      worktree's `calibration_plan.tsv` and the R14 budget (a launch queues 16 jobs and Σ = 98
      CPUs; it refuses when its worst-case end, the longest suite time plus 20 min, reaches the
      01:30 nightly), then submits one `calibration_suite.sbatch` job per suite and an `afterany`
      aggregator that writes the `calibration_full` and `calibration_attacks` rows to
      `$GIAP_ROOT/logs/nightly.tsv`. It prints `jobs=<j0>,…,<j14> rows=<id> log=<dir> cpus=<Σ>`.
    - `sbatch --array=2 scripts/slurm/nightly.sbatch`: the fast suite with the
      `calibration_nightly` tier, on `origin/main`. Its row goes to `$GIAP_ROOT/logs/nightly.tsv`,
      labelled `calibration_nightly`. A single task also posts a one-task table to the "Nightly
      (FarmShare)" issue. Add `--export=ALL,NIGHTLY_POST=dry` to print the table instead.
    - `sbatch scripts/slurm/weekly.sbatch`: K-12 and Nissle end to end with TE6, TE8 and TE2, on
      `origin/main`. Its rows go to `$GIAP_ROOT/logs/weekly.tsv`.
    - `sbatch scripts/slurm/test_branch.sbatch "$PWD"`: the full fast suite on 16 CPUs. The last log
      line is `PYTEST_EXIT=<code>`.
- Poll with `squeue -j <id>` and `sacct -j <id>`.
- Evidence: every job ID and its exit status in PROGRESS.

### 8. Readiness on R, paperwork-only PR (**agent**)

- Role: validation, on FarmShare.
- Wait until every CI check run of R on `main` has finished. A check that is still running has a
  `null` conclusion, and then H6 is `pending`. `verify` does not recompute H6 later, because H6's
  source is `farmshare`.
- Work in a worktree on a new branch at R, never by committing in the main checkout:

    ```bash
    git worktree add -b m11/readiness-<v> local/worktrees/readiness-<v> <R>
    ```

- In that worktree, run:

    ```bash
    python scripts/release/readiness.py compute --version <v> --logs "$GIAP_ROOT/logs" \
        --freeze-run <A2 freeze run dir> --out validation/release/<base>/
    ```

- This writes `readiness.json` and `READINESS.md`. Exit codes: 0 `pass`, 1 `fail`, 3 `pending`.
- Do not use `--offline`. Row H6 reads the commit's check runs with `gh api`. Without network,
  save that response and pass it with `--check-runs <file>`.
- The directory is `<base>`, also for an rc. `release.yml` reads
  `validation/release/${VERSION%%rc*}/readiness.json`.
- Open a PR `docs(release): readiness <v>` that changes only these two files, and merge it. It is
  paperwork, so `readiness.py verify` still accepts R as the readiness commit.
- Check on the merged `main`:

    ```bash
    python scripts/release/readiness.py verify validation/release/<base>/readiness.json
    ```

    This prints `commit rule: pass` and `overall: pass`, with exit 0.
- Evidence: PR URL, `overall`, the rows that are not `pass`.

### 9. Dry run of `release.yml` (**lead**)

```bash
gh workflow run release.yml -f version=<v> -f dry_run=true
gh run list --workflow release.yml --event workflow_dispatch --user khoa-yelo --limit 1 \
    --json databaseId,url,headSha,status,conclusion
gh run download <run id> -n release -D local/review/m11-release/r<k>/
```

- Check that the listed run's `headSha` is main's HEAD, so it is the run just started.
- A dry run renders the notes with `notes.py <v> --draft` and publishes nothing.
- The step summary should say `dry run: nothing published`. It also lists the readiness overall
  and the sha256 of each artefact.
- Evidence: the run URL and the head SHA in PROGRESS.

### 10. Review packet and hand-off (**lead**)

- `local/review/m11-release/r<k>/` holds:
    - `READINESS.md`;
    - the dry-run step summary;
    - `notes.md` and `dist/`, from step 9;
    - `sha256sum dist/*`;
    - the list of open non-gate issues.
- Add a PROGRESS › Questions for Khoa row "Approve release `<v>`". The step-8 PR is already
  merged, so open an issue `Release <v>: approval` that links the packet, label it `needs-khoa`,
  and message Khoa.
- Evidence: the packet path and the PROGRESS row.

## C. Each release: Khoa

### 11. Written approval (**Khoa**)

- Khoa replies "approved" (§14.8 item 5). The lead copies Khoa's words into PROGRESS before any
  dispatch.
- A row that needs a waiver gets an ADR with the header `Status: accepted-by-khoa` and the line
  `Waives readiness <id>` first. Then step 8 runs again.

### 12. Rehearsal on TestPyPI (**Khoa**)

- Required for the first release, optional later.
- Start the run from Actions → release → Run workflow (branch `main`) with `version=<v>`,
  `dry_run=false`, `target=testpypi`. From Khoa's own shell, the same run is
  `gh workflow run release.yml -f version=<v> -f dry_run=false -f target=testpypi`.
- On the run page: Review deployments → `testpypi` → Approve and deploy.
- Check the install in a fresh venv:

    ```bash
    pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ giap==<v>
    ```

- This occupies `<v>` on TestPyPI only. PyPI is a separate index, so step 13 is unaffected.
- No tag and no GitHub release are created.

### 13. The release (**Khoa**)

- Start the run from Actions → release → Run workflow (branch `main`) with `version=<v>`,
  `dry_run=false`, `target=pypi`, then approve the `pypi` environment.
- Verify fails if readiness is missing, does not verify, or is not `pass` (waived rows count as
  `pass` only through accepted ADRs).
- The workflow then:
    - uploads the wheel and sdist;
    - tags `v<v>` on the dispatched commit;
    - creates the GitHub release "GIAP `<v>`" with `notes.md` (a pre-release for an `rc`).
- Zenodo mints the DOI.
- Evidence: the run URL, which the lead records.

## D. After the release

### 14. PyPI install smoke on FarmShare (**agent**)

- Role: validation, inside the Slurm allocation, never on a login node. Use a fresh venv outside
  the checkout, with `PYTHONPATH` unset, and run from the main checkout's root (for the two input
  files):

    ```bash
    smoke=$(mktemp -d /tmp/giap-smoke.XXXXXX)
    uv venv "$smoke/venv" --python 3.12
    VIRTUAL_ENV="$smoke/venv" uv pip install giap==<v>
    "$smoke/venv/bin/giap" --version
    "$smoke/venv/bin/giap" version --json
    "$smoke/venv/bin/giap" preflight --genome tests/data/genomes/sim100k.fna.gz \
        --proteins tests/data/recorded/miniprot/0.18/sim100k/proteins.faa
    ```

- The expected results:
    - `giap --version` prints `giap, version <v>`;
    - the JSON has `"giap": "<v>"`;
    - preflight exits 0.
- Evidence: a PROGRESS line with the three outputs.

### 15. Post-release PR (**agent**)

- Role: implementer and docs. Open `docs(release): DOI and badges for <v>`.
- Add `CITATION.cff › identifiers`: the Zenodo concept DOI and the version DOI, both from the
  Zenodo record page.
- Add the README DOI badge.
- `__version__` stays `<v>` until the next release-prep PR.
- If Zenodo has no record (step 3 skipped), the PR waits and the lead notes it under PROGRESS ›
  Questions for Khoa.

### 16. Bioconda recipe (**agent**)

- Render the recipe with M11-10's `scripts/release/bioconda_recipe.py`, using the PyPI sdist sha256
  (`https://pypi.org/pypi/giap/<v>/json`) and the DOI. This script is not merged yet; take its
  options from M11-10 once it is. Its spec writes `<out>/giap/meta.yaml`, so `--out recipes`
  gives `recipes/giap/meta.yaml`.
- Hand `recipes/giap/meta.yaml` to Khoa.
- Khoa then forks `bioconda/bioconda-recipes`, adds the recipe and opens the PR from their own
  account. Agents never open PRs in other organisations' repositories. The bioconda CI does the
  lint and build checks there.

### 17. PROGRESS record (**lead**)

- Add PROGRESS › Done "Release `<v>`" with:
    - the PyPI URL, the tag `v<v>`, the release URL and the DOI;
    - the Bioconda PR, the readiness file and Khoa's approval.
- End the release freeze.
- Run `git fetch`, then `python scripts/dev/progress_bars.py --write`, then commit and push.

## E. If something fails

- **PyPI upload done, `github-release` failed**:
    - **First choice:** Khoa opens the original run and uses "Re-run failed jobs". This reuses the
      same commit and the `release` artefact that run built. Only `github-release` and `summary`
      run again, and `github.actor` stays the original dispatcher.
    - **Only if that is impossible:** Khoa dispatches again with the same inputs.
        - A new dispatch rebuilds from main's HEAD at that moment. If main has moved past the
          paperwork set, readiness is stale and verify fails.
        - The tag does not exist yet, so the tag check passes.
        - `publish` reads `https://pypi.org/pypi/giap/<v>/json`. It skips the upload when every
          file already on the index has the same sha256 as the new build.
        - It stops with `<file> is on the index with a different sha256` when the rebuild is not
          byte-identical. Then ask Khoa: the options are a patch release, or Khoa creating the
          GitHub release by hand.
- **Wrong artefact on PyPI**: never delete it. Khoa yanks it in the PyPI web UI, and the fix ships
  as the next patch version.
- **Readiness `fail` or `pending`**: nothing is published. Fix it through the normal feature loop,
  or through a waiver ADR that Khoa accepts. Then repeat steps 6–8 on the new commit.
- **Tag `v<v>` already exists** (someone made it by hand): `release.yml` refuses to run. Ask Khoa.
- **Any non-paperwork PR merged after R** (see the release freeze in step 6): `readiness.py verify`
  fails with `fail (stale: <path>)`. Repeat steps 6–8 on the new commit.

## Script reference

These options are read from the argparse definitions in `scripts/release/`.

`readiness.py compute`:

| option | required | meaning |
|---|---|---|
| `--version V` | yes | version recorded in the JSON |
| `--repo DIR` | no, default `.` | checkout read at HEAD (the readiness commit) |
| `--logs DIR` | no | directory with `nightly.tsv` and `weekly.tsv` (`$GIAP_ROOT/logs`) |
| `--freeze-run DIR` | no | A2 freeze run; R2 falls back to `mapping/mapping_fdr.tsv` there |
| `--out DIR` | yes | writes `readiness.json` and `READINESS.md` |
| `--offline` | no | touches no network; H6 becomes `pending` |
| `--check-runs FILE` | no | a saved `gh api …/commits/<sha>/check-runs` response for H6 |

`readiness.py verify READINESS_JSON [--repo DIR]` checks two things:

- the commit rule: the readiness commit is HEAD, or an ancestor that differs from HEAD only in
  paperwork (`PROGRESS.md`, `README.md`, `CITATION.cff`, `CHANGELOG.md`, `validation/release/`,
  `docs/` except `docs/methods/`, and `changelog.d/`);
- the `repo` rows, which it recomputes on HEAD.

It exits 0 on `pass` and 1 otherwise.

Other scripts:

- `notes.py VERSION [--draft]`: run from the repository root. It reads the
  `## GIAP <version> (...)` section of `CHANGELOG.md` (the base version's section for an `rc`) and
  writes `notes.md`. `--draft` falls back to `towncrier build --draft` and heads the notes `DRAFT`.
  It exits 1 when the section is missing (without `--draft`) or empty.
- `check_dist.py DIST_DIR [--version V]`: checks the wheel and sdist. It prints `ok <artefact>` or
  `FAIL <artefact>: <rule>: <path>`. Exit codes: 0 clean, 1 violations, 2 when DIST_DIR has no
  wheel or no sdist.
- `update_digests.py --digests FILE --date YYYY-MM-DD --source-commit SHA [--config FILE]
  [--nextflow-config FILE]`: used by `containers.yml`. It pins `ghcr.io/khoa-yelo/giap-<tool>@sha256:…`
  in `workflow/conf/containers.config`. With `--nextflow-config`, it also adds the pinned processes
  to the one `conda = ''` selector, found by following that file's `includeConfig` lines (it lives
  in `workflow/conf/container_common.config`, shared by the `apptainer` and `docker` profiles).
- `check_changelog.py BASE` (in `scripts/ci/`): a PR that changes `src/` or `workflow/` (except
  `workflow/tests/`) must add `changelog.d/<PR>.<type>.md` or `changelog.d/+<slug>.<type>.md`. The
  types are `feature`, `bugfix`, `doc`, `removal` and `misc`. A PR is exempt with the label
  `no-changelog` (read from the push event, so push again after adding it), or when every
  non-merge commit subject is of type `test`, `ci`, `chore` or `refactor`.

## Output tables

`readiness.json` has the top-level keys `version`, `commit`, `computed`, `overall` and `rows`. It
also has `notes`, a list. `overall` is `fail` if any row fails, else `pending` if any row is
pending, else `pass`.

Each row has these fields:

| field | content |
|---|---|
| `id` | R1–R10 (§13.2 criteria), H1–H6 (housekeeping: gates, version and tag, changelog, digests, licences, CI) |
| `criterion` | the spec's words |
| `status` | `pass`, `fail`, `pending` or `waived` |
| `value`, `threshold` | the reading and the rule |
| `evidence` | files, log rows, ADRs |
| `waiver` | `ADR-NNNN` or null |
| `source` | `repo`, `farmshare` or `both`; `verify` recomputes the `repo` rows |

`READINESS.md` renders the same rows with the columns
`id | criterion | status | value | threshold | evidence | source`. It is headed by the commit, the
computed time and the overall status, and ends with Notes. Runs of seven or more amino-acid letters
are redacted (M10-05 guard).

The `release.yml` step summary has two tables:

- `field | value`: `version`, `sha`, `readiness overall`, and `published`. `published` reads
  `dry run: nothing published`, or `publish skipped: actor <x> is not khoa-yelo`, or the publish
  and release results.
- `artefact | sha256`, one row per file in `dist/`.

## Licences of the tools GIAP runs

GIAP is MIT-licensed. This section follows `NOTICE.md`, which is the authoritative list, and
readiness row H5 checks the registry against it. `scripts/ci/check_imports.py` fails CI if GIAP
imports a GPL, AGPL or unlicensed tool.

| Tool | Licence (as in `NOTICE.md`) | How GIAP uses it |
|---|---|---|
| Bakta | GPL-3.0 | separate process |
| Promotech | GPL-3.0 | separate process |
| Promoter Calculator | GPL-3.0 (or later) | separate process |
| OSTIR | GPL-3.0 | separate process |
| TransTermHP | GPL-3.0 (or later) | separate process |
| BacTermFinder | GPL-3.0 | opt-in; separate process |
| antiSMASH | AGPL-3.0 | separate process |
| geNOMAD | non-commercial | separate process; **non-commercial use only** |
| UniOP | unlicensed (no licence file) | fetched at a pinned commit when environments are built; never vendored |
| TerminatorNet, miniprot | MIT | separate process |
| BLAST+ | public-domain (NCBI) | separate process |
| fgsea | MIT | run by `workflow/bin/giap_fgsea.R` in R 4.5 (GPL-2) as a separate process; R and fgsea are never imported |
| codon-bias | MIT | imported in-process (PyPI) |
| DIA-NN | proprietary (academic licence) | used only by users' own searches; binaries never redistributed |

What this means for a release:

- **GPL and AGPL tools run as separate programs.** GIAP never imports or vendors their code (R5).
  The wheel and sdist carry none of their code.
- **fgsea is also ported.** Its MIT code is ported to Python in `src/giap/stats/fgsea.py`, which
  runs in-process. `NOTICE.md` carries the MIT licence text, and no fgsea file is vendored.
- **Four GHCR images contain GPL-3.0 tools**, fetched from upstream unmodified:
    - Promotech;
    - Promoter Calculator (GPL-3.0-or-later);
    - OSTIR;
    - BacTermFinder.

    `NOTICE.md` gives their source links at the pinned versions and commits. Each image's
    Dockerfile is in `workflow/containers/<tool>/`. The fifth image, TerminatorNet, is MIT.
- **UniOP is never published.** There is no GHCR image for it (R5, R11).
- **geNOMAD: non-commercial use only.** Step 6 adds this to the known-limitations fragment, so the
  release notes say it.
- **Data resources:**
    - KEGG: academic use; `giap setup --kegg` downloads it on the user's machine;
    - RegulonDB: not redistributed;
    - Gene Ontology: CC BY 4.0;
    - NCBI COG 2020: NCBI FTP.
- **The contaminant library** (Frankenfield et al. 2022) has no licence file. `giap setup
  --contaminants` downloads it, and it is never vendored (R11).
- **plotly.js** (MIT) is inlined once in `report.html` when interactive figures are shown.

## Open questions for Khoa

- **The actor guard cannot tell agents from Khoa.** The FarmShare `gh` login is `khoa-yelo`, so a
  dispatch from an agent session also has `github.actor == 'khoa-yelo'`. The `publish` job's only
  remaining barrier is then the environment approval, and the same login could approve it through
  the API. Two fixes are recommended:
    - a separate agent account with write access, for the FarmShare `gh` login: the actor guard then
      works, and dry runs still dispatch;
    - or a fine-grained token without `actions:write`: agents can then not dispatch at all, and
      Khoa runs step 9.
- **PRs opened with `GITHUB_TOKEN` do not trigger CI.** `containers.yml` opens the digest PR
  `ci/digests-<sha7>` as a draft with the workflow token, so no CI starts when it opens.
    - The normal path: `ci.yml` also runs on `ready_for_review`. Once the FarmShare Apptainer smoke
      passes, `gh pr ready <PR>` from the `khoa-yelo` login starts CI.
    - Fallbacks: a plain `git commit --allow-empty` pushed to that branch from FarmShare, or Khoa
      closing and reopening the PR.
    - A GitHub App token, or the agent account above, would remove the gap.
- **R3 and R4 writers (M11-03b, done).** `nightly.sbatch` task 2 writes `calibration_nightly`;
  `calibration_full_launch.sh` (step 7) writes `calibration_full` and `calibration_attacks`. A
  collection error, a missing suite or an attack row with no attack tests fails closed. The first
  real full-tier run should confirm that task 13's attack test is never all-skipped.
- **Who may dispatch dry runs.** The default is that the lead may: a dry run uploads nothing and
  creates no tag (M11-04). `AGENTS.md › Releases` already states the two rules above.
