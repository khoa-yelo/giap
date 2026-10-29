# M11 "Release 1.0": issue index

Spec: GIPPA_SPEC.md §13 row M11. Deliverables: "packages, containers by digest, MkDocs site,
changelog, citation file, preprint outline". **Gate:** "the acceptance criteria below [§13.2];
**Khoa approves** the release" (§14.8 checkpoint 5).

Sources read: §13 row M11 and §13.2 (all seven bullets), §12 (`release.yml`: "verify gates,
changelog and the version; build sdist and wheel; publish to PyPI (trusted publishing); create the
tag and GitHub release (the Zenodo DOI comes from the GitHub integration). Bioconda recipe PR as a
manual follow-up"; `containers.yml`: "build per-tool images, push to `ghcr.io/<owner>/gippa-<tool>`,
open a PR updating the digests"; the **Releases** paragraph: "semantic versioning, changelog
fragments (towncrier), PyPI and Bioconda packages, containers pinned by digest, and a Zenodo DOI
each time"), §3.10 ("never push tags"; "Tags and releases come from `release.yml`
(workflow_dispatch) or from Khoa"), §4 (layout: `CHANGELOG.md`, `CITATION.cff`, `NOTICE.md`,
`docs/` = "MkDocs: user guide, methods/, adr/, specs/", extras `docs`), §9.2 ("semantic
versioning; breaking column changes only in major versions"), §11 ("a real E. coli run gates every
release"; e2e "weekly and before release"), §14.1 (docs role), §14.2 ("Every PR … includes a
changelog fragment"), §14.6 (pins by digest), §14.8 items 5–7, §16 (nothing here is deferred),
R1, R3, R5 ("geNOMAD is non-commercial: document that in the user guide"), R6, R11, R12, R13,
Appendix A.5 (docs agent: "Document every CLI option and output table from the code and schemas …
licence notes"). Backlog: "(Khoa, 2026-10-24) running GIPPA on macOS is untested … Candidate for
M11: a Docker-profile test run, and documenting the platforms." Prior specs: M2-01 (images, ADR-0008),
M10-12 (`docs/user-guide/runtimes.md`, question L9), M10-15 (`metrics.py checklist()`).

One issue = one PR unless noted; size counts **code** only (ADR-0007: `src/`, `workflow/`,
`scripts/`; workflows under `.github/` are counted here too, to be safe). Prose (user guide,
runbook, outline, fragments) is not gated but is listed.

## How M11 is sequenced (Khoa's standing direction: MVP first)

- **Pass 1 (MVP: an installable, tag-ready release path), 01–05.** A wheel and sdist that install
  into a clean environment and find their own Nextflow workflow (01); a 1.0 changelog and a complete
  citation file (02); a release-readiness checker that turns §13.2 into a pass/fail table with
  evidence (03); `release.yml` that verifies, builds and, only when Khoa dispatches it with
  `dry_run=false`, publishes and tags (04); the runbook of the steps only Khoa can take (05).
  **Pass 1 is done when a `dry_run=true` dispatch of `release.yml` on main builds both artefacts,
  installs the wheel, renders the release notes and prints the readiness table** (the table may
  still show `pending` rows: the evidence comes from M6–M10 and from M11-13).
- **Pass 2 (polish and the gate), 06–13.** GHCR images by digest for the conda-only tools (06),
  the MkDocs site with generated references (07), the user guide and tutorials (08), platforms and a
  Docker-profile smoke (09), the Bioconda recipe (10), the preprint outline (11), the
  changelog-fragment check and the README (12), and the gate: release candidate, readiness on the
  final commit, Khoa's approval and Khoa's dispatch (13).

## Already built (checked on main a6bf2e96, 2026-10-24)

| Piece | State |
|---|---|
| `pyproject.toml` | hatchling; `name = "gippa"`, **`version = "0.0.1"` written twice** (pyproject and `src/gippa/__init__.py`; also `CITATION.cff`); extras `dev`, `docs` (`mkdocs-material`); `gippa` console script; `[tool.towncrier]` with `directory = "changelog.d"`; wheel `packages = ["src/gippa"]` |
| Workflow in the wheel | **missing.** `gippa.runner.nextflow.workflow_main()` looks for `<checkout>/workflow/main.nf`, then `gippa/workflow/main.nf` inside the package, but nothing copies `workflow/` into the wheel. An installed `gippa run` finds no workflow |
| `CHANGELOG.md` | header and the towncrier marker only; **no release section** |
| `changelog.d/` | 2 fragments (`m0-scaffold.feature.md`, `m6-11-regulondb.feature.md`) for 317 merged PRs; §14.2's "every PR includes a changelog fragment" is not enforced |
| `CITATION.cff` | `cff-version 1.2.0`, title, one author (Khoa Hoang), licence, repository, `version: 0.0.1`; no `type`, `date-released`, `abstract`, `keywords`, `identifiers` |
| `NOTICE.md`, `LICENSE` (MIT) | present; NOTICE lists every tool, the plotly.js inline, the fgsea port |
| `README.md` | says "Status: under development (milestone M0, scaffold)" |
| `.github/workflows/` | `ci.yml`, `golden-r.yml` only. **No `release.yml`, `containers.yml` or docs workflow** |
| Containers | `workflow/conf/containers.config`: 7 processes by digest (Bakta, TransTermHP, geNOMAD, antiSMASH, miniprot, BLAST, TBLASTN). Six Dockerfiles under `workflow/containers/` (uniop, promotech, promoter_calculator, terminatornet, bactermfinder, ostir) are **never built**; those processes stay conda-only. The UniOP Dockerfile says "built for local use only, never pushed publicly" (unlicensed, R5/R11). fgsea is conda-only by M8-17 decision 9 |
| MkDocs | **no `mkdocs.yml`, no user guide.** `docs/` holds `adr/` (57), `methods/` (29 notes, 4 with `$…$` math), `specs/`, `backlog.md`. `mkdocs` 1.6.1 and `mkdocs-material` 9.7.7 are in the shared venv. `docs/user-guide/runtimes.md` is promised by M10-12 (L9) |
| Generated references | `scripts/write_config_schema.py --check` (config JSON Schema) exists; `src/gippa/schemas/` has the Frictionless `datapackage.json`, `tables/`, manifest schemas. No CLI or table reference page |
| Release evidence | `validation/pilot/metrics.py checklist()` (§15.6 rows), `validation/pilot/metrics/metrics.json` (committed by M10-15b), `benchmarks/baseline.json`, `$GIPPA_ROOT/logs/nightly.tsv` and `weekly.tsv` (date, job, commit, exit, summary). No release checker |
| Accounts | repository `khoa-yelo/gippa` is **public** (the spec assumed private; Pages and Zenodo both need public). PyPI `gippa` and Bioconda `gippa`: **not taken** (404 on 2026-10-24). No Pages site. No tags, no releases. All runtime deps exist on conda-forge or bioconda (`codon-bias` on bioconda only) |

## Issues

| Order | Issue | Pass | Scope | Depends on | Code lines | Test lines | Khoa |
|---|---|---|---|---|---|---|---|
| 1 | [M11-01](M11-01-package-build.md) | 1 | single-source version; `workflow/` bundled as `gippa/workflow`; sdist/wheel content rules; `scripts/release/check_dist.py`; CI `package` job (build, check, clean-venv install smoke) | none | 190 | 230 | no |
| 2 | [M11-02](M11-02-changelog-citation.md) | 1 | towncrier rendering rules; retroactive milestone fragments M0–M10; `CITATION.cff` complete and schema-checked; version agreement test | none | 30 | 160 | Q-K7 (authors) |
| 3 | [M11-03](M11-03-release-readiness.md) | 1 | `scripts/release/readiness.py`: §13.2 bullets + release hygiene → `validation/release/<v>/readiness.json` + `READINESS.md` | 01, 02 (soft) | 340 | 320 | Q-K5, Q-K6 (criteria readings) |
| 3b | [M11-03b](M11-03b-nightly-readiness-labels.md) | 1 | `calibration_full` / `calibration_attacks` rows: task records + JUnit in `calibration_full.sbatch`, `calibration_full_launch.sh`, `scripts/ci/calibration_rows.py` | 03 | 300 | 300 | no |
| 3c | [M11-03c](M11-03c-calibration-full-sharding.md) | 2 | full tier fits its time limit: one job per suite with per-core lanes (replicate shards + unsharded test groups), shard-0-only tests auto-routed to group lanes, committed plan, block 2 on lanes (stats), launcher CPU fit; 5 PRs | 03b; blocks 13 | 640 | 900 | no |
| 4 | [M11-04](M11-04-release-workflow.md) | 1 | `.github/workflows/release.yml` (workflow_dispatch, `dry_run` default true, actor and environment guards, trusted publishing, tag + GitHub release by the workflow) | 01, 03 | 190 | 170 | Q-K1–Q-K4 (account set-up) |
| 5 | [M11-05](M11-05-release-runbook.md) | 1 | `docs/release.md`: the agent's prep steps and Khoa's account and dispatch steps, in order; post-release follow-ups | 04 | 0 (≈ 200 prose) | 40 | runs it |
| 6 | [M11-06](M11-06-ghcr-containers.md) | 2 | `containers.yml`: build and push 5 GHCR images (not UniOP), digest PR; `scripts/release/update_digests.py`; Apptainer smoke on Slurm | 01 | 210 | 160 | Q-K8 (package visibility) |
| 7 | [M11-07](M11-07-mkdocs-site.md) | 2 | `mkdocs.yml`, nav, `scripts/docs/gen_reference.py` (CLI, config, tables; `--check`), `docs.yml` (strict build; Pages deploy) | 01 | 260 | 180 | Q-K9 (Pages) |
| 8 | [M11-08](M11-08-user-guide.md) | 2 | user guide and tutorials (install, setup, quickstart, inputs, outputs, interpretation, licences, runtimes) | 07 | 0 (≈ 1,200 prose) | 60 | no |
| 9 | [M11-09](M11-09-platforms-docker-smoke.md) | 2 | platforms page; Actions `docker-smoke` job (two pinned images under `-profile docker`) | 07, 06 (soft) | 90 | 60 | Q-K10 (Actions minutes) |
| 9b | [M11-09b](M11-09b-bare-image-refs.md) | 2 | drop `docker://` from `containers.config` pins (bare refs; Apptainer cache names unchanged); writers, readers, H4 guard; ADR-0059 | 06, 03; blocks 09 | 40 | 150 | no |
| 10 | [M11-10](M11-10-bioconda-recipe.md) | 2 | `packaging/bioconda/meta.yaml` rendered from pyproject; dependency agreement test; Khoa's recipe PR | 01, 04 | 130 | 110 | opens the bioconda PR |
| 11 | [M11-11](M11-11-preprint-outline.md) | 2 | `docs/preprint/outline.md`: sections, figure plan, claim → evidence map, limitations | 03, 08 (soft) | 0 (≈ 250 prose) | 50 | Q-K11 (venue, authors, AI statement) |
| 12 | [M11-12](M11-12-changelog-check-readme.md) | 2 | `scripts/ci/check_changelog.py` + CI step; README for 1.0 | 02 | 120 | 130 | no |
| 13 | [M11-13](M11-13-release-gate.md) | 2, gate | release-prep PR (1.0.0, built changelog, CITATION date), evidence jobs on the RC commit, committed readiness, review packet; Khoa approves and dispatches; post-release PRs (DOI, Bioconda, PyPI smoke) | all | 30 | 60 | **yes** (§14.8 item 5) |
| 14 | [M11-14](M11-14-container-install.md) | 2, after the rebrand | container-complete pip install: content-addressed immutable GHCR tags (`<version>-<key>`, `images.json`, authenticated existence check), fgsea image, UniOP built locally from the authors' archive by `setup` (sha256-checked, terms accepted by the user; 14u1/14u2), `enabled: true|false|auto` for failed builds (`auto` default only after the statistician's attack and operating point, 14s1/14s2/14d2), `setup --images`, preflight rule for pinned images, docker/apptainer profile parity, version bump policy, rc and `reusable(C, R)` evidence reuse, compare version warnings, install + upgrade pages (14a–14j; 14b/14g/14h/14p are CI changes); ADR-0060 | 01, 04, 06, 09, 09b, 03; rename PR before 14j | 1,590 | 2,220 | **yes** (first GHCR publish, repo settings, CI changes, Q-K13–Q-K19) |
| 15 | [M11-15](M11-15-portable-clusters.md) | 2, after the rename | lowest-barrier install and clusters: profile auto-detection, `gippa demo`, zero-flag `setup` (all DBs on, licence notices), Nextflow launcher fetch + Java check, generic `slurm` profile + executor flags, `--nf-config`, `singularity` profile, step image (`GIPPA_STEP` in a container), offline nodes, site configs, second tested cluster (15a–15v, 15c split in 15c1/15c2; 15n is a CI change) | 01, 10, 14, 09, 09b | 2,570 | 2,680 | **yes** (CI 15n, step-image publish, second cluster, Q-K36–Q-K43) |
| 16 | [M11-16](M11-16-translation-tables-4-25.md) | 2 | genetic codes 4 and 25 alongside 11 (`annotation.translation_table`; supersedes Q-K30): table module, preflight (allowed set, tool-option conflicts, `genome.genetic_code` both directions), Bakta/miniprot/tblastn plumbing, six-frame index format 3 keyed by the table, stop sets and recoded sites, codon usage, novel, variants, makedb, manifest/report, gippasim `genetic_code`, calibration arms, M. hyorhinis code-4 validation (16a–16g); ADR-0067 | 14w, #413 | 1,070 | 1,750 | **yes** (D1: protected `defaults.yaml` tblastn default) |
| 19 | [M11-19](M11-19-nextflow-docker-pull-never.md) | 2 | N6b: only `gippa setup` pulls. `--pull=never` in the docker profile's `docker.runOptions`; preflight requires every pinned `containers.config` ref the selected modules use to be local under docker (`docker image inspect`, reasons naming `gippa setup --images --profile docker`), Docker >= 20.10 (podman >= 3.0, best effort); shared with the #465 fgsea and local-UniOP docker arms; docs + changelog (19a) | 14 (14f1, 14f2, 14p), #465 | 80 | 380 | no (Q1-Q12, lead) |
| 20 | [M11-20](M11-20-initiator-methionine.md) | 2 | initiator M at annotated non-ATG starts: start sites and cells (mirrored into decoys, as M7-16), exact-route credit for an M-initial query at a cell, six-frame index format 4, flag `initiator_m`, variants alignment/pileup, TE1/TE8 proxies; ADR-0068 (`stats`) | 16, #481 | 190 | 420 | **yes** (Q-K1: ADR-0068 vs §2.2 improvement 6) |

Total ≈ 1,590 code lines (≈ 750 in pass 1), ≈ 1,730 test lines, ≈ 1,950 prose lines, before M11-14
(added 2026-10-25: ≈ 1,590 code, ≈ 2,220 test, ≈ 380 prose; its questions Q-K13–Q-K19 and Q-L11–Q-L15
are in the table below and in the issue file).

## Dependency DAG

```text
Pass 1   01 ─┬─► 03 ─► 04 ─► 05
         02 ─┘
Pass 2   01 ─► 06 ─┐
         01 ─► 07 ─► 08 ─► 11
               07 ─► 09
         04 ─► 10
         02 ─► 12
                   └─────────► 13 (all of the above; Khoa approves, Khoa dispatches)
```

## Waves (≤ the lead's agent cap)

| Wave | Issues | Notes |
|---|---|---|
| A | 01, 02 | disjoint files (pyproject packaging vs. changelog.d/CITATION); 02's pyproject edit is the `[tool.towncrier]` table only, merge 01 first |
| B | 03, 06, 07 | 06 and 07 share no files with 03 |
| C | 04, 08, 12 | 04 reads 03's JSON format |
| D | 05, 09, 10, 11 | 05 closes the MVP pass |
| E | 13 | after M6–M10 gates are recorded; ends at Khoa |

## What M11 needs from other milestones

| Evidence (§13.2 bullet) | Comes from | Status 2026-10-24 |
|---|---|---|
| A2 TE1 ≥ 99.5 %, mapping FDR ≤ 1 % | M10-15b `metrics.json` (`te1`), the A2 freeze run's `mapping/mapping_fdr.tsv` | M10 freeze pending |
| calibration suites and attack nulls | `calibration_full.sbatch` + `-m calibration` on the RC commit (M11-13) | gates M6-09 (D-1), M7, M8 open items |
| TE assertions | M10-06 (`te_real`, `te7`), weekly TE8 | M10-06 merged; TE2–TE5 strictness is Q-K5 |
| byte-identical tables | TE7 (integration and pilot) | in place |
| ≈ 35 min (16–55) on 16 threads / 64 GB | `benchmarks/baseline.json › standard_run.wheat.k12_standard` (M10-12) | M10-12 evidence PR pending |
| pilot passes | M10-15 gate block in PROGRESS with Khoa's review | pending |

## Questions for the lead

Numbered Q-K (needs Khoa: an account step or a §14.8 checkpoint) and Q-L (the lead decides). Each
has a reversible default the specs already use. Details are in the issue files.

| # | Question | Default in use | Issue |
|---|---|---|---|
| Q-K1 | PyPI: register a **pending trusted publisher** for project `gippa` (owner `khoa-yelo`, repo `gippa`, workflow `release.yml`, environment `pypi`), and the same on TestPyPI? | yes, both; the first non-dry dispatch targets TestPyPI | 04, 05 |
| Q-K2 | Who may dispatch a publishing run? | `github.actor == 'khoa-yelo'` **and** GitHub environments `pypi`/`testpypi` with Khoa as required reviewer | 04 |
| Q-K3 | Zenodo: enable the GitHub integration for `khoa-yelo/gippa` before the first release; metadata from `CITATION.cff` (no `.zenodo.json`)? | yes | 05, 13 |
| Q-K4 | Version string | `1.0.0` (PEP 440 = semver here); no `rc` on PyPI — the rehearsal uploads `1.0.0` itself to TestPyPI (a separate index) before the PyPI dispatch | 04, 13 |
| Q-K5 | §13.2 says every available true-example assertion "passes"; §15.6 lets TE2–TE5 be "reported; failures explained". Which reading gates the release? | strict: a TE that is `reported (not passed)` makes readiness `fail`; Khoa may waive by an accepted ADR, which readiness then cites | 03 |
| Q-K6 | Is the runtime bullet ("about 35 min (range 16–55)") a pass/fail criterion? | yes: median `k12_standard` in [16, 55] min is `pass`; outside is `fail` unless waived by ADR | 03 |
| Q-K7 | `CITATION.cff` authors, ORCID, affiliation; acknowledge AI assistance? | Khoa Hoang only, no ORCID until given, no email (R12), no affiliation; no AI statement in CFF (the preprint decides) | 02 |
| Q-K8 | GHCR packages created by `containers.yml` are private on first push; make the 5 images public? | yes, Khoa flips each package to public once; UniOP is never pushed | 06 |
| Q-K9 | Host the docs on GitHub Pages (Settings › Pages › Source: GitHub Actions)? | yes; until enabled the deploy job is skipped and CI only builds | 07 |
| Q-K10 | Spend ≈ 5 Actions minutes per PR on a Docker-profile smoke (Khoa's backlog item)? | yes, only on PRs touching `workflow/` and nightly on main | 09 |
| Q-K11 | Preprint: venue, author list, AI-use statement, which pilot figures may appear | bioRxiv (Bioinformatics); Khoa sole author pending; an AI-use statement section left for Khoa; aggregates only (Khoa's data rule) | 11 |
| Q-K12 | Release 1.0 with M6–M9 issues still open (M6-09 D-1, M7-16, M8-11, M9-25)? | readiness requires every milestone gate recorded in PROGRESS; non-gate open issues are listed under "Known limitations" in the 1.0 notes for Khoa to accept | 03, 13 |
| Q-K13 | UniOP for pip users: licence request, local build, distance-model fallback; behaviour under `--publication` | **answered by Khoa 2026-10-25 (install mechanism):** `setup` fetches UniOP from the authors (sha256-checked) and builds it locally after the user accepts the terms, so preflight passes with UniOP (14u1/14u2); `enabled: auto` only for failed builds, the default only after 14s1/14s2 (14d2); `--publication` unchanged; optional: Khoa asks the authors for a licence (then published, local build a fallback) | 14 |
| Q-K14 | GHCR owner after the rename | `ghcr.io/khoa-yelo/gippa-<tool>` | 14 |
| Q-K15 | Release candidates on TestPyPI only, or also PyPI pre-releases? | TestPyPI only | 14 |
| Q-K16 | GHCR retention | never prune; if Khoa prunes: unreleased digests only, ≥ 180 days, never the latest per tool | 14 |
| Q-K17 | Nextflow + Java as prerequisites, or a `setup --nextflow` download? | prerequisites | 14 |
| Q-K18 | Approve every `ghcr` build, or only the first? | every build until Khoa removes the reviewer | 14 |
| Q-K19 | fgsea image vs. M8-17 decision 9 (conda-only) | image | 14 |
| Q-L1 | sdist contents | `src/`, `workflow/` (without `workflow/tests`), README, LICENSE, NOTICE, CHANGELOG, CITATION, pyproject; **no** tests, docs, validation, scripts, benchmarks | 01 |
| Q-L2 | Installed-package smoke cannot use `gippa run --profile test,stub` (`_PROFILES` rejects `test`/`stub`) | the smoke parses the installed workflow with `nextflow config -profile test,stub` and `-profile farmshare,apptainer`, and asserts `workflow_main()` resolves inside site-packages; no CLI change | 01 |
| Q-L3 | Retroactive changelog granularity | one fragment per milestone M0–M10 (plus the two existing), written by the docs role from PROGRESS and the indexes | 02 |
| Q-L4 | Readiness is computed on FarmShare (it reads `$GIPPA_ROOT/logs`) and committed; `release.yml` cannot see FarmShare | the committed `readiness.json` names commit R; `release.yml` passes only if HEAD = R or R..HEAD touches only release-paperwork paths | 03, 04 |
| Q-L5 | Pin third-party Actions by commit SHA in the new workflows | yes in `release.yml`, `containers.yml`, `docs.yml` (tag in a comment); `ci.yml` unchanged | 04, 06, 07 |
| Q-L6 | Put ADRs and specs on the site? | ADRs yes ("Decisions"), specs no (`exclude_docs`), backlog no | 07 |
| Q-L7 | Math in methods notes (4 notes use `$…$`) | `pymdownx.arithmatex` (generic) + MathJax 3 from jsDelivr on the site only | 07 |
| Q-L8 | A `gippa` core image on GHCR? | no: the spec asks for per-tool images only (R1) | 06 |
| Q-L9 | Changelog-fragment check severity | fail for PRs that change `src/` or `workflow/` without a fragment; label `no-changelog` exempts | 12 |
| Q-L10 | Bioconda `run:` dependencies beyond pyproject | add `nextflow` (bioconda) and `openjdk >=17`; Apptainer/conda left to the user | 10 |
| Q-L11 | Image cache without `$GIPPA_ROOT` | `NXF_APPTAINER_CACHEDIR` → `$GIPPA_ROOT/containers` → `<db_dir>/containers` | 14 |
| Q-L12 | BacTermFinder lock files | resolve on FarmShare in 14b | 14 |
| Q-L13 | Separate `output_schema_version` from the config's `schema_version` | yes, `1.0` from 1.0.0 | 14 |
| Q-L14 | Keep the M11-06 `:<sha7>` tags | no (`images.json` `built_from` + OCI label) | 14 |
| Q-L15 | `NON_STATS_PATHS` and the never list | as in M11-14 item 4 (statistician-reviewed on #367); changes need `stats` review | 14 |
