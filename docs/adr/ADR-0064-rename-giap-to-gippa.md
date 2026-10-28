# ADR-0064: Rename GIAP to GIPPA
- Status: accepted-by-khoa (2026-10-27, in session, Q-K31: accept as proposed; rename all output names now; one-off PR-size exception granted; Khoa-only steps handed over when the rename PR is ready)
- Date: 2026-10-27
- Spec sections: title and every section that names GIAP (§0.2 R3, R14; §3; §5.1; §6.3; §7.9.2;
  §9.4; §12; §14.6); M11-13; M11-14 §6 ("Sequencing with the rebrand"); ADR-0004, ADR-0007;
  ADR-0060 §5
- Needs Khoa: yes (overrides the spec's name; changes the package, CLI, env-var and output names;
  needs a PR-size exception, `.github/` edits and GitHub, PyPI and FarmShare actions that only
  Khoa can take)

## Context
On 2026-10-25 Khoa decided to rebrand GIAP to GIPPA before the first real publish (M11-13). The
rebrand covers the repository, the Python package, the CLI and the code
(docs/backlog.md › Rebrand). Planned order: this ADR, then one scripted rename PR whose outputs
are byte-identical apart from names, then simplification PRs gated by goldens and calibration.
GIAP_SPEC.md is read-only and names GIAP throughout, so the rename needs an ADR that overrides it.
On 2026-10-27, read-only checks of `https://pypi.org/pypi/gippa/json` and
`https://test.pypi.org/pypi/gippa/json` both returned HTTP 404, so the name `gippa` is free on
both. Nothing was registered.

The repository has 22,181 case-insensitive occurrences of `giap` in 2,243 of its 2,920 tracked
files. The measured inventory and the list of places where a blind `sed` would do damage are in
`docs/specs/rename-gippa.md`. That page is outside the site, so read it
in the repository.

## Decision
1. **Name.** The product is GIPPA. The distribution, import package and CLI are `gippa`. The
   source tree is `src/gippa/`, and `giap.cli:main` becomes `gippa.cli:main`. Classes with Giap in
   the name are renamed (`GiapConfig` → `GippaConfig`). The simulator is renamed too
   (`tests/sim/giapsim` → `tests/sim/gippasim`, `import gippasim`).
2. **One mechanical PR, made by a committed script.** `scripts/dev/rename_gippa.py` takes an
   explicit keep list, maps case-preserving substrings (`giap`→`gippa`, `GIAP`→`GIPPA`,
   `Giap`→`Gippa`) and renames paths with `git mv`. It changes no logic.
   - It uses plain substrings, not `\b` word boundaries. `\b` misses `"\tgiap"`, `"\ngiap"` and
     `b"giap"`, and these occur 20+ times in GFF writers and tests.
   - It never edits gzip or other binary files, the keep-list paths or the `GIAP_SPEC` token.
   - `--check` lists every `giap` left after the rename by keep category, and fails on any
     uncategorised hit.
   - The PR contains only the script's output, the script, the comparator (item 8) and the
     regenerated files (config schema, reference pages, `uv.lock`). Simplification comes later,
     in separate PRs.
3. **What stays as written.**
   - GIAP_SPEC.md, including its file name. The 1,171 `GIAP_SPEC` references stay valid, and
     this ADR maps the name.
   - The historical ADRs (ADR-0001 to ADR-0063), docs/specs/, the PROGRESS.md history and
     docs/backlog.md history.
   - Recorded tool data and every PROVENANCE file under tests/data/. These are recordings, and
     several are pinned by sha256.
   - validation/ reports and pilot artefacts from past runs.
   - Merged commit messages and the tag-less history.

   Living docs are renamed: README, docs/user-guide, docs/methods, tutorials, mkdocs.yml,
   CITATION.cff, NOTICE.md, AGENTS.md, CLAUDE.md and .claude/. Pending changelog.d fragments are
   renamed too, because they render into the first CHANGELOG.
4. **Environment variables: a clean break in code, an alias only in the FarmShare shell.**
   - All 63 `GIAP_*` names become `GIPPA_*`, except `GIAP_SPEC`.
   - `src/` and `workflow/` read only `GIPPA_*`. They have no fallback, because no release
     exists and nobody outside this project sets the old names. A fallback added now would also
     have to be removed later, and that removal is another src change that invalidates evidence.
   - `scripts/farmshare_env.sh` exports each `GIPPA_*` and also, for a transition window, the old
     `GIAP_*` name with the same value. The old exports serve jobs already in flight, old
     worktrees and the installed scrontab.
   - The old exports are removed after the first publish. Their removal is a scripts/ change, so
     it rides with the next change that already needs new evidence.
   - `$GIAP_ROOT` and `$GIAP_DB_DIR` keep their on-disk locations. Only the variable names change.
5. **FarmShare and Slurm.**
   - Job names become `gippa-*`: 44 distinct `giap-*` names in scripts/slurm and the launchers.
   - R14 accounting (scripts/resume_status.sh, session_start.sh, the CPU-budget check and
     AGENTS.md) matches `(giap|gippa)-*` until no `giap-*` job is left. After that it matches
     `gippa-*` only.
   - The repository checkout path `/home/classes/bios/270/khoa/giap` and the contents of
     `$GIAP_ROOT` do not move.
   - The 15 micromamba envs `giap-<tool>` (including `giap-wf`, which is on `PATH`) become
     `gippa-<tool>`. The marker `.giap-built` becomes `.gippa-built`, and the file `GIAP_COMMIT`
     inside each env and image becomes `GIPPA_COMMIT`. All envs are rebuilt with
     `build_envs.sbatch`. Measured builds took 0:28 to 6:28 min each.
   - The old envs are deleted only after the new evidence passes.
   - Khoa reinstalls the scrontab, because its job name and log name change.
6. **Containers.**
   - Image names become `ghcr.io/khoa-yelo/gippa-<tool>`, as M11-14 §6 and Q-K14 already
     assume.
   - Dockerfile labels (`org.gippa.image.key`, `image.title`), `/opt/gippa/` and
     `workflow/bin/gippa_fgsea.R` change with them.
   - Every image key changes, because the key hashes those files and their paths. That is
     expected: `images.json` is `{}` and nothing has been published.
   - Locally built Apptainer images (fgsea, UniOP via `setup`) are rebuilt.
7. **Output names change now, before 1.0.0.** No `v*` tag exists, so the M11-14 §4 schema-diff
   rule does not apply yet. The changes are these:

   | Output | Before | After |
   |---|---|---|
   | `run_manifest.json`, `preflight.json`, `makedb_manifest.json`, VCF header | `giap_version` | `gippa_version` |
   | VCF `##source` | `##source=giap-<v>` | `##source=gippa-<v>` |
   | GFF3 column 2 | `giap` / `GIAP` | `gippa` / `GIPPA` (each writer keeps its case) |
   | six-frame ORF IDs and decoys | `GIAPORF_…` / `rev_GIAPORF_…` | `GIPPAORF_…` / `rev_GIPPAORF_…` |
   | datapackage name | `giap-run` | `gippa-run` |
   | report and compare HTML | title and text "GIAP" | "GIPPA" |
   | `version --json` | `{"giap": …}` | `{"gippa": …}` |
   | manifest tool row | `giap` | `gippa` |
   | schema `$id` | `github.com/khoa-yelo/giap/…` | `github.com/khoa-yelo/gippa/…` |
   | Nextflow processes | `GIAP_STEP`, `GIAP_MAPPING` | `GIPPA_STEP`, `GIPPA_MAPPING` |

   - Table column names and config keys contain no `giap`, so neither changes. The config
     schema changes only its title (`GiapConfig`).
   - Readers do not accept the old names. Pre-rename runs are not inputs to anything after the
     rename, because all evidence is regenerated (item 9).
   - `compare` on a pre-rename run treats the version as not recorded.
   - The pilot's A2 database was built with `six_frame_orfs: false`, so no pilot search depends
     on `GIAPORF_`.
   - `output_schema_version` stays `1.0` at the first release.
   - **One exception: the SVG id salt stays `"giap"`.** It is `HASHSALT` in `io/plot_files.py`
     and `svg.hashsalt` in `viz/style.py`. It is an invisible constant, and changing it would
     change every SVG id, so SVGs could not be shown to be byte-identical. A code comment
     records why it stays.
8. **Verification: outputs byte-identical apart from names.**
   - **Commits.** B is the main tip the rename branches from. A is the head of the rename PR.
   - **Inputs.** One input set is made once with B and used for both runs:
     - giapsim sim100k at the nightly seed;
     - the two snapshot fixtures `sim_empty` and `sim_plots`;
     - `makedb --six-frame-orfs` on sim100k, which exercises the ORF IDs;
     - a `compare` of two runs.

     Each is run at B (`giap …`) and at A (`gippa …`). Both runs use
     `-profile farmshare,apptainer`, the same seed and the same thread count, on the same node
     class, through `giap_run.sbatch`.
   - **Comparator** (`scripts/dev/rename_compare.py`, committed in the PR):
     - Pairs files by relative path, with `N` applied to the path. The two file sets must be
       equal.
     - `.gz` files are decompressed, and Parquet is read to canonical TSV.
     - Text is compared after `N` maps both sides' `GIPPA|Gippa|gippa|GIAP|Giap|giap` to one
       token per case (`@N@`, `@n@`, `@Nn@`).
     - **Sequence guard.** FASTA sequence lines and sequence or peptide columns are compared
       exactly, without `N`, because GIAP is a valid tetrapeptide. It occurs 33–38 times in
       each recorded K-12 and Nissle Bakta proteome.
     - PNG and PDF must be byte-equal.
     - The comparator ignores only `git_commit`, `git_dirty`, the provenance digest and the
       volatile fields that the R6 byte-identity check already ignores. Every sha256 recorded in
       a manifest is recomputed on A's file and must match.
     - **Pass** means zero differences. Any difference, such as a sort order moved by the longer
       ORF prefix, is listed in the PR, and the PR is fixed or the difference is justified.
   - **Suites at A.** On Slurm: the every-commit suite, the snapshot suite (`--mpl`, baselines
     unchanged), the goldens, the integration suite and `nf-test`. Then `rename_gippa.py --check`.
   - **Re-pins.** Sim-output sha256 pins can change only where the simulator writes its own name,
     as in `source=giapsim`. Each such re-pin is listed in the PR with the comparator line
     showing that its content is equal under `N`.
9. **Release evidence is regenerated on the renamed commit.** The rename touches src/, tests/sim/,
   workflow/, ci/ and `scripts/release/readiness.py`. Readiness is self-protected: SELF_PATH
   changes, and its NON_STATS_PATHS and never-list entries `src/giap/…` become `src/gippa/…`,
   which needs `stats` review (Q-L15). So `reusable(C, R)` and `code_equal` reuse no pre-rename
   row.

   | Row | Evidence | Re-run | Recent wall time |
   |---|---|---|---|
   | R3 | `calibration_full` + `calibration_nightly` rows | `calibration_full_launch.sh`: 15 suites × 4 CPU (60 CPU) | ≈ 6 h (array 1894441: 13 min to 4 h per task, 44.6 h serial test time; limit now 6 h, #422) |
   | R4 | attack rows of the same runs | same jobs | (included) |
   | R5 (TE8) | `weekly.tsv` row `k12_mg1655` (code-equal on RUN_PATHS) | `weekly.sbatch` | ≈ 1.5–2 h (giap-weekly-annot 1:17–2:08 h) |
   | R6 | nightly `integration` row (code-equal on INTEG_PATHS) | `nightly.sbatch` (4 × 16 CPU) | ≈ 3 h sharded (#420; tasks timed out at 6 h before sharding) |
   | R1, R2, R9 | `validation/pilot/metrics/metrics.json`, `REPORT.md` (file-read, not commit-gated) | pilot pass 4 (A1 + A2 + compare) | ≈ 3–4 h (head jobs 1:37–1:52 h each, in parallel; tool jobs re-run because the process names change the Nextflow task hashes) |
   | R7, R8 | `benchmarks/baseline.json` `standard_run` (file-read) | nightly task 1 benchmarks; the standard run comes from the pilot or weekly run | inside the runs above |
   | H6 | CI on the release commit | automatic | — |

   - R1, R2 and R9 are not commit-gated by readiness, but the pilot is re-run anyway, so the
     release evidence carries the release's name and commit.
   - R14: calibration (60 CPU) and nightly (64 CPU) together use 124 of the 128 CPUs. Pilot and
     weekly therefore run after one of them.
   - Total ≈ 6–7 h with everything overlapped, or ≈ 12 h sequential, plus the env rebuilds
     (< 1 h).
10. **PR-size exception.** The rename PR changes about 15,600 lines in about 1,810 files, of which
    about 2,000 lines are in src/ and workflow/. That is far above the ~400-line rule (ADR-0007).
    By the ADR-0004 precedent, Khoa grants a one-off exemption. Review covers the script, the
    keep list, the `--check` output, the comparator report and `git diff --stat`, not every line.
11. **Rollout order.**
    1. Khoa accepts this ADR.
    2. Freeze: merge or park every open PR, and let running agents and GIAP Slurm jobs finish.
    3. Run the B-side runs of item 8.
    4. The rename PR. Khoa makes or approves its `.github/` part (the list is in the plan); agents
       never edit .github/.
    5. Run the A side, the comparator and the suites. Merge with Khoa's approval: it is a CI
       change.
    6. Khoa renames the GitHub repository `khoa-yelo/giap` → `khoa-yelo/gippa` (git and web URLs
       redirect; GitHub Pages does not), updates `origin` in the FarmShare checkout and adds
       PyPI and TestPyPI trusted publishers for `gippa` → `khoa-yelo/gippa`.
    7. On FarmShare:
       - update `farmshare_env.sh` (in the PR);
       - fast-forward the main checkout and `uv sync` there, which makes the editable install
         `gippa`;
       - `build_envs.sbatch` for all envs;
       - rebuild the local images;
       - Khoa reinstalls the scrontab.
    8. Regenerate the evidence (item 9), run the M11-13 readiness check, then M11-14 14j (the
       first GHCR publish) and the PyPI publish.
    9. Simplification PRs, gated by goldens and calibration.
12. **Rollback.**
    - Before the first publish, revert the single squash commit in a PR and rename the
      repository back.
    - The old `giap-*` envs and the transition aliases are still in place until the evidence
      passes, so FarmShare keeps working.
    - After the first publish there is no rollback of names, because PyPI and GHCR names are
      permanent. Fix forward.

## Evidence (before → after: metrics, tests, simulations)
- Measured inventory, comparator design and hazard list:
  `docs/specs/rename-gippa.md`.
- The "after" evidence is item 8's comparator report and suite results. They are attached to the
  rename PR, not to this ADR.

## Consequences
- The release is one renamed-commit evidence cycle later, ≈ 6–12 h of Slurm wall time plus
  review.
- The published package, CLI, images and outputs never carry the name GIAP. Only the read-only
  spec, the SVG salt and the historical records keep it.
- Open branches made before the rename must merge main and re-run `rename_gippa.py` on their own
  diff before they can merge.
- Nextflow `-resume` caches from before the rename do not hit.
- Readers of old docs follow GIAP_SPEC.md, and this ADR maps the old name to the new one.

## Amendment 2026-10-27 (Q-K32): every document is rewritten; new expansion
Khoa decided (Q-K32, in session: "Re write everything") that after the rename the documentation
does not mention GIAP at all. The decision text above is kept as accepted. This amendment changes:

- item 2: "never edits ... the `GIAP_SPEC` token" no longer holds; the token is renamed;
- item 3: the spec and the historical records are rewritten (below);
- item 4: "except `GIAP_SPEC`" no longer holds; `GIAP_SPEC` is renamed with the spec;
- the second Consequences bullet ("Only the read-only spec, the SVG salt and the historical
  records keep it"): only the SVG salt, the data keep paths and this rename record keep GIAP;
- the last Consequences bullet ("Readers of old docs follow GIAP_SPEC.md"): it is false after the
  rename, because the file moves to `GIPPA_SPEC.md` and no document names GIAP. Readers of old
  docs (merged commits, old checkouts) follow this ADR, which maps the old name to the new one.

The rest of the decision stands.

- **Rewritten too:** GIAP_SPEC.md, which moves to `GIPPA_SPEC.md` (its content and the
  `GIAP_SPEC` tokens are renamed with it, so the token is no longer protected), the historical
  ADRs ADR-0001 to ADR-0063, docs/adr/index.md, docs/specs/, PROGRESS.md, docs/backlog.md,
  docs/validation/ and validation/reports/. validation/release/ holds only `.gitkeep`; a file
  added there is rewritten unless it is sha256-pinned.
- **Kept as written (data only):** PROVENANCE files, recorded tool data, the pilot artefacts
  (REPORT.md, metrics, db_a2, inputs, the sha256-pinned DIA-NN cfgs), binary and sequence files,
  the SVG salt, `.giap-complete`, `.giap_setup.lock`, the checkout path, the rename tooling and
  `.github/` (Khoa). The rename record also keeps the old name: this ADR, its index line and the
  inventory `docs/specs/rename-gippa.md`.
- **New expansion:** GIPPA stands for "genome integrated prokaryotic proteomics analysis". Every
  variant of the old expansion ("genome-integrated analysis of proteomics data", with or without
  "for bacterial genomes") is replaced by it, with the first letter capitalised where the
  original was.
- `rename_gippa.py --check` also fails on any `giap` left in README.md, CITATION.cff, mkdocs.yml,
  AGENTS.md, CLAUDE.md, GIPPA_SPEC.md, PROGRESS.md, docs/ and every other tracked `*.md` / `*.rst`
  (NOTICE.md, CHANGELOG.md, changelog.d/, .claude/agents/, validation/), outside the rename
  record, the data keep paths and the kept data tokens (and the AGENTS.md R14 note of item 5,
  which names the old job prefix).
- Moving GIAP_SPEC.md and editing ci/protected.txt are protected-file changes; the rename PR
  carries `Gate-Change-Approved: ADR-0064`.
