# M11-14: container-complete install with versioned artifacts

Milestone M11 · **pass 2, after the rebrand** · depends on: M11-01 (wheel bundles `workflow/`),
M11-04 (`release.yml`), M11-06 (`containers.yml`, `update_digests.py`), M11-09 (docker smoke),
M11-09b (bare refs, ADR-0059), M11-03 (readiness H4), ADR-0022 (operon distance model); the
GIAP → GIPPA rename PR (backlog "Rebrand") before any publish · estimate ≈ 1,590 code lines +
≈ 2,220 test lines in 14 code/test PRs and one docs PR (table at the end), plus one Khoa action (14j)
· **`stats` (14d, 14d2, 14s1, 14s2, 14g)**: 14d changes which model builds the TUs, which are the
exchangeability blocks of most §8 nulls; 14g decides when calibration evidence is reused ·
roles: implementer (`scripts/release/`, `src/giap/preflight/`, `src/giap/setupdb/`,
`src/giap/io/fgsea_r.py`, `src/giap/compare/`, `src/giap/report/`, `workflow/`), test-author
(`tests/unit/test_m11_14_*.py`), **statistician** (owns
`tests/calibration/attacks/test_attack_operon_distance_fallback.py`, the ADR-0022 addendum and the
`NON_STATS_PATHS` list), docs (`docs/get-started/install.md`, `docs/get-started/upgrade.md`,
`docs/release.md`) · versioning policy:
[ADR-0060](../adr/ADR-0060-release-versioning-and-image-keys.md) (proposed) · **Khoa: yes**: the
first GHCR publish and the repository settings (14j), the CI changes in 14b, 14g and 14p, and
Q-K13–Q-K19.

Khoa's direction (2026-10-25): "`pip install` + preflight passes" with every tool in a container,
and "iterate on the final version fast. Make sure that whatever is built can be iterated to have
1.1, 1.2, 2.0, and so on."

## Spec sections implemented

- §3 Containers row ("referenced by digest"), §6.3 container field
  `ghcr.io/<owner>/<pkg>-<tool>@sha256:…`, §6.x profiles (`docker`, `apptainer`, "for users
  elsewhere"), G19.
- §7.0 Environment rows (`environment.tools`: "binaries, environments and DB versions of enabled
  modules found → else error (exit 2)").
- §8 nulls whose blocks are TUs (§8.5 attacks); §7.1.3/§7.4 operons (ADR-0022 fallback); R3, R4.
- §9.2 Versioning row: "`schema_version` in the manifest; semantic versioning; breaking column
  changes only in major versions"; §12 Releases: "semantic versioning … containers pinned by
  digest"; §12 `containers.yml`; §14.6 "`containers.config` by digest"; §17 "prefer `ghcr.io`
  images built by `containers.yml`".
- §3.10 and AGENTS.md: agents never push tags and never publish; only Khoa publishes.
- R5/R11 (UniOP unlicensed: never redistributed), R6 (official sources only), R14 (cluster caps).
- Backlog: "Off-FarmShare install (release blocker, Khoa question)" (this issue answers it) and
  "Rebrand" (images are published only under the new name).

## The target 1.0 user path

```text
1. pip install gippa                       # Python package; workflow/, containers.config, images.json inside
2. Docker or Apptainer, Nextflow ≥ 25.04.8 + Java ≥ 17 on PATH   (prerequisites, Q-K17)
3. gippa setup --db-dir DB --images --profile docker|apptainer [--accept-uniop-terms]
                                           # databases, pre-pulled images, local UniOP build (item 1b)
4. gippa preflight ... --profile docker|apptainer                  # passes, with UniOP
5. gippa run ... --profile docker|apptainer
```

No repository checkout, no `build_envs.sbatch`, no `$GIAP_ROOT/micromamba`. Nextflow pulls every
tool image by digest (docker), or runs it from the cache that step 3 filled (apptainer). Below,
`<pkg>` is the PyPI project name read from `pyproject.toml` (`giap` today, `gippa` after the
rename); `<owner>` is the GitHub owner (`khoa-yelo`, Q-K14).

**UniOP (Khoa, 2026-10-25: "UniOP is published research, we should be able to allow a mechanism
for installing it").** UniOP has no licence file, so GIAP never redistributes it (no published
image, no vendoring; R5/R11). Step 3 instead fetches the pinned commit from the authors' GitHub on
the user's machine, checks its sha256 and builds it locally (item 1b, 14u1/14u2), after the user accepts
the terms. Step 4 then passes **with UniOP** under the default `annotation.uniop.enabled: true`.
The distance-model fallback (`enabled: auto`, item 1 (c)) is kept only for installs where that
build fails (offline, no build rights), and it becomes the default only after 14s1/14s2 (14d2).

## Definition

### 1. Image coverage (every non-in-process tool)

| Tool | Image today | After M11-14 | PR |
|---|---|---|---|
| bakta, transtermhp, genomad, antismash, miniprot, blast (BLAST, TBLASTN) | upstream, by digest | unchanged | none |
| promotech, promoter_calculator, terminatornet, bactermfinder, ostir | Dockerfile, never built (M11-06) | `ghcr.io/<owner>/<pkg>-<tool>` by digest, content-addressed tag (item 2) | 14a, 14b, 14j |
| fgsea | conda-only, started by the `function` step (M8-17 decision 9, ADR-0038) | `ghcr.io/<owner>/<pkg>-fgsea` with `giap_fgsea.R` baked in; run by `giap.io.fgsea_r` when no env or prefix is found | 14c |
| uniop | conda-only; unlicensed | **never published.** Built locally from the authors' archive by `<pkg> setup --uniop`, which `setup --images` runs by default (14u1, 14u2); `annotation.uniop.enabled: true \| false \| auto` (14d) for the failure path; `auto` becomes the default only after 14s1 + 14s2 (14d2); Khoa may ask the authors for a licence | 14u1, 14u2, 14d, 14d2, 14s1, 14s2 |
| codon_bias | in-process | unchanged | none |
| cascadia, timsconvert (`PILOT_TOOLS`), DIA-NN | pilot only | out of scope: not on the user path | none |

**fgsea (14c).**
- The image is built from `workflow/containers/fgsea/Dockerfile`, which installs
  `workflow/envs/locks/fgsea.explicit.txt` and copies `workflow/bin/giap_fgsea.R` to
  `/opt/<pkg>/giap_fgsea.R`. So `workflow/bin/giap_fgsea.R` is a key input (item 2): a script
  change rebuilds the image.
- Its ref lives **only** in `workflow/conf/images.json` (`"processes": []`). No
  `withName: 'FGSEA'` line is written, because no such Nextflow process exists and the selector
  would log "no process matching config selector" on every run. H4 and `update_digests.py` accept
  an `images.json` entry with no process for `fgsea` only.
- Resolution order in `resolve_launcher`:
  1. `GIAP_FGSEA_PREFIX`. An explicit prefix **never falls through**: a bad prefix still raises
     as today, so the FarmShare profiles (`prebuilt_envs.config` sets it) behave as before.
  2. The prebuilt `giap-fgsea` env (marker equal to the env file's sha256).
  3. The image, run through the runtime that `giap run` passes in `GIAP_CONTAINER_RUNTIME`
     (`docker` | `apptainer`; unset → `FgseaEngineUnavailable`). The cache directory comes from
     the same helper as preflight (item 3; 14a).
- The image argv (AC 7 pins it):
  - docker: `docker run --rm --cidfile <call_dir>.cid -u <uid>:<gid> -v <call_dir>:<call_dir> -w <call_dir> -e LC_ALL=C <ref> Rscript --vanilla /opt/<pkg>/giap_fgsea.R <the nine options of command()>`;
  - apptainer: `apptainer exec --cleanenv -B <call_dir> --pwd <call_dir> <cache>/<simpleName>.img Rscript --vanilla /opt/<pkg>/giap_fgsea.R <options>`.

  Lead amendment (2026-10-26, #386 review B1/S1): `--vanilla` stays on the image path, because
  apptainer binds `$HOME`, so the user's `~/.Rprofile`, `~/.Renviron` and `~/R` library would
  otherwise load. A killed `docker run` client does not stop its container (dockerd owns it, and R
  as PID 1 ignores SIGTERM), so docker writes its container ID to `<call_dir>.cid`, a sibling of
  the bound directory that survives its removal. On timeout or any exception the engine runs
  `docker kill <id>` and then `docker rm -f <id>`; failures there are logged and do not hide the
  original error. The `.cid` file is removed afterwards. Apptainer stays in the process group and
  needs no extra step.

  The options name files relative to the call directory (as `command()` does today), so only that
  directory is bound. The shebang front end is not read under the image path; the image runs
  `Rscript` explicitly.
- No GSEA Nextflow process is added (M8-17 decision 9 stays: the step starts R itself). R is GPL-2
  and fgsea MIT; the image label and NOTICE say so.

**UniOP (14u1, 14u2, 14d, 14d2, 14s1, 14s2).** Three options, not exclusive:

| Option | What happens | Default |
|---|---|---|
| (a) licence from the authors | publish `<pkg>-uniop` like the others; the local build (b) becomes a fallback | optional: Khoa may ask (Q-K13) |
| (b) local build (14u1, 14u2) | item 1b below: fetched from the authors, built on the user's machine | **the install path** (`setup --images` includes it) |
| (c) degrade (14d) | `enabled: auto` resolves to `true` when UniOP is installed for the profile, else `false`; operons use the ADR-0022 distance model | **failures only** (offline, no build rights); available in 14d, the default only after 14d2 |

*(c) in detail (14d, `stats`).*
- **Config.** `annotation.uniop.enabled: true | false | auto` (new; config model, `defaults.yaml`
  and the JSON Schema). 14d ships the default `true` (today's behaviour). 14d2 changes the default
  to `auto` in a one-line PR. 14d2 may merge only after 14s1 and 14s2 have merged and the
  statistician has signed off in the PR.
- **Preflight resolution** (in `src/giap/preflight/`; the resolved value goes into
  `resolved_config.yaml`):

  | `enabled` | UniOP installed | `--publication` | Result |
  |---|---|---|---|
  | `true` | yes | any | run UniOP |
  | `true` | no | any | exit 2, `environment.tools: uniop` missing |
  | `auto` | yes | any | resolves to `true` |
  | `auto` | no | no | resolves to `false`; warning `environment.tools: uniop not installed (unlicensed, not distributed); operons use the distance model; run: <pkg> setup --uniop …` (item 3) |
  | `auto` | no | yes | exit 2, `environment.tools: uniop not installed; set annotation.uniop.enabled: false to publish with the distance-model operons (see Methods)` |
  | `false` | any | any | distance model, flag and caveat below (accepted under `--publication`) |

- **Run.** `giap run` passes `params.run_uniop` (from the resolved value), and `UNIOP` is skipped
  when it is false.
- **`operon_pairs`** takes its reason from the resolved config, not from the missing input. The
  reasons:
  - `UniOP not installed` when `auto` resolved to false;
  - `UniOP disabled (annotation.uniop.enabled=false)` for an explicit false;
  - `UniOP failed (…)` (today's text) only when UniOP ran and produced no pair.

  In each of these cases the step emits the run flag `operon_source_distance_model` (`DISTANCE_FLAG`,
  already in `RUN_FLAGS`) **and** `layer_missing:uniop`, so the manifest's `layers_missing` lists
  `uniop`.
- **Manifest** (schema minor) gains `operon_source` (`uniop` | `distance_model` | `precomputed`)
  and `operon_source_reason`, copied from `operon_meta.json`.
- **Report.** When the source is `distance_model`, the summary card and the Methods section say:
  "Operons/TUs: intergenic-distance model (<reason>); precision X [CI], recall Y [CI] on
  RegulonDB K-12 at the pair threshold; TU-based tests use these units". X and Y are read from
  `calibration.yaml: operon_distance` (14s2). Today `report/` has no operon-source text.
- **Compare** (14h) warns `compare.operon_source: run A uniop, run B distance_model` when the
  sources differ. `giap compare` bootstraps over the first run's TUs, so mixing sources changes
  the blocks.

*Statistical evidence that gates 14d2 (statistician-owned; R3: no bound, α, R or tolerance
moves).*
- **14s1: attack `tests/calibration/attacks/test_attack_operon_distance_fallback.py`** (§8.5).
  The fast tier runs under `-m calibration`; the R of the existing gates goes under
  `calibration_full`.
  - **Arm A, class re-derivation.** Run giapsim `operon_null`, `operon_attack_ar1` and
    `operon_attack_detection` with `uniop_p` replaced by
    `distance_model_p(intergenic_bp, d₀, s)` from `calibration.yaml` and `source =
    distance_model`. Then run the M6-09 discordance level and missed-operon decoy-FDR checks at
    their existing α, bounds and R.
  - **Arm B, K-12-like gaps.** Draw each same-strand pair's intergenic distance from the
    class-conditional ADR-0022 bin table (operon and non-operon counts per bin), so the distance
    model splits and merges true TUs at K-12 rates. Truth is generated from the true TUs, and the
    distance-model TUs are the analysis blocks. The arm feeds them into:
    - the TU-collapsed GSEA and the TU-counted ORA at level (the arms of
      `test_attack_gsea_structure` (a)–(c) and `test_attack_ora_tusize` 4b);
    - the determinants TU-block bootstrap CI level
      (`test_attack_determinants_tu_null::test_univariate_ci_level`);
    - the TU-blocked autocorrelation null.

    The same arm also runs once with UniOP-like errors (the Beta classes), so the two sources are
    compared on the same data.
  - **Outcome rule.** Every arm that passes at the existing bounds clears the default. If an arm
    fails, the default may still ship only if the affected rows carry a flag
    (`tu_source_distance_model`) with an R4 caveat and an ADR scoping them as uncalibrated under
    the fallback. The bound never moves.
- **14s2: operating point.**
  - Evaluate precision and recall of `p_dist ≥ pair_threshold` (0.5, Khoa's decision, unchanged)
    on RegulonDB 14.5.0 K-12, with the directon-bootstrap CIs (the M6-12 driver's
    `precision_ci`/`recall_ci`, t fixed).
  - Write the results to `calibration.yaml: operon_distance` (`precision_at_pair_threshold`,
    `recall_at_pair_threshold`, each with `ci_low`/`ci_high`) and to an ADR-0022 addendum. No
    new threshold and no re-selection.
  - Add one arm to `weekly.sbatch`: K-12 standard with `annotation.uniop.enabled=false`. Record
    the TU agreement with RegulonDB next to the UniOP arm in `weekly.tsv`. The arm is reported,
    not gated, and adds no bound.

**1b. UniOP local build (14u1, 14u2).**
- **Command.** `<pkg> setup --uniop --db-dir DB --profile P [--cache-dir DIR]
  [--accept-uniop-terms]`. It replaces the 14e draft's `--tools uniop --runtime`.
- **Runtime from the profile list (one rule for setup, preflight and `giap run`).** `P` may be a
  comma list. The resolved runtime `R(P)` is the first rule that matches, in this order:
  1. **micromamba**, if any profile in the list enables conda (`conda`, `farmshare`,
     `farmshare_local`);
  2. **apptainer**, if `apptainer` is in the list;
  3. **docker**, if `docker` is in the list;
  4. otherwise none: `uniop` is missing (`no container or conda profile`).

  So the micromamba rule wins: `farmshare,apptainer` (the `giap_run.sbatch` profile) resolves to
  micromamba, which is also the only FarmShare path that needs no fakeroot. Lists that combine
  `conda` with `docker` or `apptainer` (and no `farmshare*`) are not supported for UniOP. Preflight
  reports `uniop` missing with `ambiguous UniOP runtime for <P>: use one of conda | docker | apptainer`.
- **Included by default.** `setup --images --profile P` runs the same step whenever `uniop` is in
  `required_tools(cfg)` for the resolved config (true for the defaults; false with
  `annotation.uniop.enabled: false` or precomputed operon pairs). `--no-uniop` skips it.
- **Fetch (one download, one check, three runtimes; 14u1).**
  - The archive is a new `DbItem` row of `src/giap/data/db_sources.yaml`, laid out exactly like the
    others:
    - `name: uniop_source`;
    - `version: git:3ed14dfc2a5d69b993d99a26e50c1f8e9723c087`;
    - `url: https://github.com/hongsua/UniOP/archive/3ed14dfc2a5d69b993d99a26e50c1f8e9723c087.tar.gz`;
    - `host: github.com`;
    - `checksum: sha256:<hex>`;
    - `path: uniop/UniOP-3ed14dfc2a5d69b993d99a26e50c1f8e9723c087.tar.gz` (in `db_dir`, so it passes
      `_safe_name`);
    - `tool: null`;
    - `check: sha256`;
    - `size_gib` left at the default.

    The source is official (R6), and its commit equals `annotate/tools.py` (ADR-0008). `plan()`
    includes the row only under `--uniop`, or under `--images` when UniOP is in the set, never in
    plain `setup`.
  - **Redirect.** GitHub answers `302 → codeload.github.com`. 14u1 adds
    `"github.com": frozenset({"codeload.github.com"})` to `HttpsFetcher.REDIRECT_HOSTS` (an
    R6-reviewed allowlist change). Codeload streams without `Content-Length`. The fetcher needs
    none, and no `Content-Length` handling is added.
  - The existing `check_file` verifies the sha256 before anything is built. A mismatch deletes the
    file and exits 2: `uniop: archive sha256 <got> != pinned <want>`.
  - No builder fetches anything itself. `setup` assembles a staging context
    `<cache>/uniop/context/` from the bundled files, at the repo-relative paths the Dockerfile
    `COPY`s (`workflow/envs/uniop.yml`, its locks), plus the Dockerfile, `uniop.def` and the
    verified archive as `uniop-src.tar.gz`.
  - **Permissions.** `<cache>/uniop/` and `<db_dir>/uniop/` are created `0700`; the archive and
    the `.sif` are `0600`. A lab-wide cache therefore never shares one user's UniOP copy, and the
    install page says so.
- **Image key (14a, 14u1).**
  - `key(t)` lives in `src/giap/io/image_key.py`, so an installed wheel can compute it from the
    bundled `workflow/`. 14a puts it there. `scripts/release/image_key.py` is the CLI over it, and
    readiness H4 calls the same function.
  - The file mode comes from the exec bit (`100755` if any x bit is set, else `100644`), which
    matches git's two modes without needing git.
  - The uniop key gains one extra line, `uniop-src.tar.gz\t100644\t<pinned sha256>`, so a pin
    bump changes the tag. The AC 1 `COPY`-source test accepts `uniop-src.tar.gz` through that
    line.
  - Preflight compares `local_build.key` with `key("uniop")` computed from the installed package.
- **Terms notice (always printed, one line; 14u2):** `UniOP (github.com/hongsua/UniOP @3ed14dfc)
  is fetched from its authors and built on this machine; it has no licence file, so <pkg> does not
  distribute it and whether your use is permitted is yours to decide.`
  - **Default: not accepted unless the user says so.** On a TTY, setup asks
    `Build UniOP for your own use? [y/N]`. Without a TTY and without `--accept-uniop-terms`, UniOP
    is not built. Its item is recorded `absent` with the message
    `terms not accepted; rerun with --accept-uniop-terms`.
  - **Exit codes.** Under explicit `--uniop`, a refusal or a missing acceptance exits 2. Under the
    implicit `setup --images`, the other items still install and setup exits 0 with a warning, so
    a scripted install never hangs and never fails on UniOP alone. Preflight then names the
    command (item 3). Once the terms are accepted, a failed fetch, checksum or build exits 2 in
    both forms, after the other items are installed.
  - **Why this default.** With no licence, GIAP has no permission to pass on and cannot accept on
    the user's behalf. A silent yes in a non-interactive install would make GIAP the one deciding.
    A default of no costs a scripted install one flag, and pip users following the install page
    copy that flag knowingly.
- **Build per runtime.** Nothing is ever pushed.
  - **docker (14u1):** `docker build -f <context>/Dockerfile -t localhost/<pkg>-uniop:<version>-<key> <context>`.
    The docker CLI only. Podman is out of scope for 1.0. On FarmShare (podman, no docker) the
    routes are the prebuilt env or apptainer.
  - **apptainer (14u2):** `workflow/containers/uniop/uniop.def` is new and is a key input.
    - It mirrors the Dockerfile: `Bootstrap: docker` from the same digest-pinned micromamba base,
      `%files` for `uniop.yml` and `uniop-src.tar.gz`, and `%post` doing the same install and
      extraction into `/opt/conda/share/uniop` plus `GIAP_COMMIT`, then `chmod -R a+rX` on that
      directory.
    - `%post` has **no `chown`**. In the fallback below only one uid is mapped, and the Dockerfile's
      `--chown=$MAMBA_USER` has no `.def` equivalent.
    - **Build and fallback.** As root: `apptainer build <sif> <def>`. Otherwise:
      1. `apptainer build --fakeroot <sif> <def>`;
      2. if that exits non-zero, retry once with `apptainer build --fakeroot --ignore-fakeroot-command <sif> <def>`.

      The `<sif>` is `<cache>/localhost-<pkg>-uniop-<version>-<key>.sif`. Step 2 is the mode that
      works on FarmShare's SUID apptainer 1.4.4: plain `--fakeroot` falls back to libfakeroot there
      and fails with `GLIBC_2.38 not found` (#370 review S1).
    - **"No fakeroot" means both attempts failed.** There is no `/etc/subuid` pre-check, because it
      says "no" on FarmShare even though step 2 succeeds. The item is then recorded `absent`
      with `apptainer build failed (fakeroot unavailable or build error): <last stderr line>; use --profile docker or conda, or set annotation.uniop.enabled=auto`.
  - **micromamba (14u2), plain `conda` profile:** creates `<cache-dir>/envs/giap-uniop` from the
    lock, then extracts the archive into `<prefix>/share/uniop` with `GIAP_COMMIT`, which is where
    `UNIOP_HOME` defaults in `modules/local/uniop/main.nf`.
  - **micromamba, `farmshare*` lists (the prebuilt case):** `prebuilt_envs.config` already pins
    `UNIOP { conda = $GIAP_ROOT/micromamba/envs/giap-uniop }`. `setup --uniop` builds nothing
    there. It checks that this prefix exists with `share/uniop`, and records it with
    `local_build.runtime: prebuilt`. Otherwise it exits 2 with
    `sbatch scripts/slurm/build_envs.sbatch uniop`.
- **`local_tools.config` (14u1 writes it; each builder adds its block).**
  - The file is one `profiles { … }` with one block per runtime that was built, named for the
    profile:

    ```groovy
    profiles {
        docker    { process { withName: 'UNIOP' { container = 'localhost/<pkg>-uniop:<version>-<key>'; conda = '' } } }
        apptainer { process { withName: 'UNIOP' { container = '<absolute .sif path>'; conda = '' } } }
        conda     { process { withName: 'UNIOP' { conda = '<prefix>' } } }
    }
    ```
  - It **never sets `conda.enabled`**, because conda is a global switch. Conda stays on only where
    the `conda` and `farmshare*` profiles turn it on.
  - **When `giap run` adds it.** `giap run` adds `-c <db_dir>/local_tools.config` only when the
    file exists and `R(P)` is docker, apptainer, or micromamba under plain `conda`. It **never**
    adds the file under a `farmshare*` list. There `prebuilt_envs.config` owns UNIOP, so a local
    `.sif` and the host env are never combined.
  - A block applies only under its own profile. A docker build followed by `-profile apptainer`
    therefore finds no apptainer block, and preflight reports `uniop` missing; Apptainer never
    tries to pull `docker://localhost/…`.
- **Dockerfile (14u1).** `ADD https://…/archive/${GIT_SHA}.tar.gz` becomes
  `COPY uniop-src.tar.gz`, so the build reads only the verified archive. 14b's `ADD --checksum`
  change therefore skips uniop. The header says "built only locally by `<pkg> setup --uniop` from
  the authors' archive; never pushed (unlicensed)". `update_digests.py` keeps refusing `UNIOP`,
  and `containers.yml` never builds it.
- **CI never builds real UniOP.** CI, the docker smoke (M11-09) and `weekly.sbatch` never run the
  real `setup --uniop`, and never upload a built image, `.sif` or archive as an artifact. The
  weekly UniOP arm keeps using the FarmShare prebuilt env.
- **What is recorded (14u2).**
  - The `db_manifest.json` item `image:uniop` has `version: git:3ed14dfc…`, the archive `url`,
    `checksum: sha256:<archive>`, `path` (the `.sif`, the env prefix, or `""` for docker) and
    `status`.
  - A new optional object `local_build` holds `{"runtime": "docker" | "apptainer" | "micromamba"
    | "prebuilt", "ref", "digest", "key", "terms": "flag" | "prompt" | ""}`. `terms` is `""` for
    `prebuilt`, which setup did not fetch.
  - The run manifest's `uniop` entry has `version`, `local_build: true`, `runtime`,
    `source_sha256` and `container_digest` (the same schema minor as 14d). `container_digest`
    means something different per runtime: for docker, `sha256:<image ID>` (a local config digest,
    not a registry manifest digest); for apptainer, `sha256:<.sif file bytes>`; for micromamba and
    prebuilt, it is absent.
  - Compare (14h) therefore does not warn on the uniop `container_digest`. It warns on
    `source_sha256`, the key or the runtime, so a rebuild of the same key on another host is not
    reported as a tool change.

### 2. Content-addressed, immutable image tags (14a, 14b; ADR-0060 §1)

- **Image key.**
  - For tool `t`, the inputs are every file under `workflow/containers/t/`, plus
    `workflow/envs/t.yml`, plus `workflow/envs/locks/t.explicit.txt` and `t.pip.txt` when
    present. fgsea adds `workflow/bin/giap_fgsea.R`.
  - `key(t)` is the first 12 hex of sha256 over the UTF-8 text made of one line
    `<repo-relative POSIX path>\t<git file mode: 100644 or 100755>\t<sha256 hex of the file bytes>\n`
    per input.
  - The lines are sorted by the Unicode code points of the path (Python `sorted` on `str`). The
    file mode is included because `COPY` preserves the exec bit.
  - A static test parses each Dockerfile's `COPY` sources and asserts each one is an input, so
    the key cannot miss a file the build reads.
- **Tag.** `<version>-<key>`. `<version>` is the registry version, with `git:<40 hex>` written
  `git-<first 8 hex>`, and every character outside `[A-Za-z0-9_.-]` replaced by `_`. Examples:
  `ostir:1.1.3-3f9c0a1b2d4e`, `promotech:git-56251ad9-0b7e19c4aa21`. A tag must fully match
  `[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}`.
- **Locked, checked builds (14b).**
  - Each Dockerfile installs from the explicit lock (and pip lock) when one exists, as
    `scripts/build_envs.sh` does.
  - Each `ADD https://github.com/.../archive/<sha>.tar.gz` gains `--checksum=sha256:<hex>`, so an
    upstream archive cannot change under a key.
  - BacTermFinder has no lock today. 14b adds `locks/bactermfinder.{explicit,pip}.txt`, resolved
    on FarmShare by the existing lock script, or states why it cannot (Q-L12).
  - A static test asserts that each Dockerfile's `GIT_SHA` and each env-file version pin equal
    the registry entry in `annotate/tools.py`.
- **Immutable tags (14b).**
  - `containers.yml` pushes only the `<version>-<key>` tag. It drops M11-06's `:<sha7>` tag and
    the bare `:<tool version>` tag, which is mutable.
  - **Existence check.** The plan job logs in to GHCR with `GITHUB_TOKEN` (`packages: read`) and
    runs `docker buildx imagetools inspect ghcr.io/<owner>/<pkg>-t:<tag>`. The result is
    classified:
    - exit 0 means **exists**: the tool is not rebuilt, and its digest goes into the digest PR;
    - an error whose message matches `manifest unknown|not found|NAME_UNKNOWN` (registry 404)
      means **absent**;
    - **any other error** (401/403, network, rate limit) **fails the job**. An anonymous failure
      on a still-private package can therefore never be read as "absent".
  - **Re-check before push.** The build job repeats the same check, with the same
    classification, immediately before `docker/build-push-action`. It pushes only if the tag is
    still absent, else it fails `tag exists: <ref>`.
  - **Race.** The concurrency group becomes the global `containers` (not per ref), with
    `cancel-in-progress: false`, so no two runs of the workflow overlap. Pushes from outside the
    workflow are forbidden (ADR-0060 §2), so no window remains.
  - The anonymous check lives only in `release.yml` verify (14g), where it means "not public".
- **Only changed tools rebuild.** A GIAP release in which `ostir`'s inputs did not change keeps the
  same `ostir` tag and digest as the release before.
- **Lock file.**
  - `workflow/conf/images.json` is new and bundled in the wheel. Per GHCR tool it lists
    `{"key", "tag", "ref" (bare digest ref), "processes", "built_from" (commit), "licence"}`.
  - `update_digests.py` writes it in the same PR as `containers.config`.
  - `scripts/release/image_key.py --check` recomputes every key. It fails with
    `stale image: <tool> (inputs changed since <key>)` when a key input changed without a new pin.
- **No deletion, no retagging (ADR-0060 §2).** Agents never delete or retag GHCR packages or
  versions. Every digest named in a released `containers.config` or `images.json` (any `v*` tag) is
  kept forever. Digests that no release names may be pruned only by Khoa, not before 180 days after
  their build, and never the latest one per tool. Public GHCR storage is free, so the default is
  "prune nothing" (Q-K16).

### 3. Preflight and setup (14f1, 14f2)

`missing_tools(tools, profile, giap_root)` becomes `missing_tools(tools, profile, image_cache,
env_root, local_tools)` (14f1). The image-cache helper `image_cache(environ, db_dir)` is added in
14a, so 14c and 14f1 share it.

Conda is enabled only by the `conda` profile and by `prebuilt_envs.config` (the `farmshare*`
profiles). Under `docker` alone **and under plain `apptainer`**, a process that has only a
`conda` directive runs on the host without its environment. The table follows from that:

| Tool kind | `docker` | `apptainer` (no conda profile) | `conda`, or any `farmshare*` profile |
|---|---|---|---|
| pinned image (upstream or GHCR) | present (Nextflow pulls by digest) | present iff `<image_cache>/<simpleName>.img` exists | as today |
| pinned image **and** a prebuilt env | the image is used; the env is ignored | same | as today (env) |
| conda-only (after M11-14: only a micromamba-built local UniOP) | **missing**: `conda-only tool under a container profile without conda: <tool>` | same | needs the prebuilt env (as today) |
| uniop with a local build (`local_tools.config` + `local_build`) | present iff the local tag exists and its image ID equals `local_build.digest` (`docker image inspect`, a fake in tests) | present iff the `.sif` exists and its sha256 equals `local_build.digest` | plain `conda`: the `conda` block's prefix with `share/uniop/GIAP_COMMIT` equal to the pin; `farmshare*`: the prebuilt env (as today; `local_tools.config` not used) |
| uniop, not built, stale or wrong runtime (no item, `absent`, `local_build.key` ≠ `key("uniop")` of the installed package, no block for `R(P)`, or an ambiguous list) | missing; the message names the command (below) | same | same |
| uniop with `enabled: auto` / `false` | the item 1 table | same | same |
| fgsea | prefix → env → image (item 1) | same | as today (prefix) |
| in-process | never missing | never missing | never missing |

- `image_cache` is `NXF_APPTAINER_CACHEDIR`, else `$GIAP_ROOT/containers`, else
  `<db_dir>/containers`. A pip user with no `$GIAP_ROOT` therefore has a defined place (Q-L11).
  `giap run` passes the same directory to Nextflow (`apptainer.cacheDir`).
- The "no prebuilt env needed" rule applies exactly to tools that have a `withName` pin in the
  bundled `containers.config` (plus fgsea through `images.json`). The FarmShare behaviour is
  unchanged.
- **`<pkg> setup --images --profile apptainer|docker [--cache-dir DIR]` (14f2).**
  - Which images: those `required_tools(cfg, abundance=True)` names for the config that `setup`
    resolves from its existing `--config`/`--set` options (default: the default modules). With
    `--images all`, every pinned image. fgsea is included when `function` is in `modules`.
  - Apptainer: `apptainer pull <cache>/<simpleName>.img docker://<ref>`, the file name Nextflow
    looks for (ADR-0059). Docker: `docker pull <ref>`.
  - It verifies that the pulled digest equals the pin, and records
    `{"name": "image:<tool>", "version": "<tag or upstream tag>", "checksum": "<digest>",
    "status": "ok"}` items in `db_manifest.json`.
  - This and `setup --uniop` (which `--images` includes, item 1b) are the only commands that
    touch a registry or the network for images.
- **UniOP message.** When `uniop` is missing under `enabled: true`, preflight exits 2 with
  `environment.tools: uniop not installed; run: <pkg> setup --uniop --profile <profile> --db-dir <db_dir> --accept-uniop-terms`,
  filled in with the run's own profile and DB directory. A stale build says `uniop: local build
  <key_old> does not match this release (<key_new>)` and the same command. The reason from
  `local_build`/the item message is appended when there is one (for example `terms not accepted`,
  `no fakeroot`). Under `enabled: auto` the same text is the warning (item 1 (c)). The
  `--publication` rule is unchanged.
- **`ToolRunner` (14f2)** runs the Bakta and antiSMASH downloaders. It gains a `docker run -u
  <uid>:<gid> -v …` path, and drops its requirement for `GIAP_ROOT` or a checkout
  (`setupdb/run.py:108-113`). It looks up the micromamba env, then the apptainer cache from
  `image_cache`, then docker, so `<pkg> setup` itself needs no checkout and no micromamba.
- **Database pins per release.** `src/giap/data/db_sources.yaml` ships in the wheel and pins
  every database per release: version, URL, checksum. `<pkg> setup --check` reports an installed
  item whose version or checksum differs from this release's pin as `mismatch`, with the message
  `installed <v_old>, this release pins <v_new>`. This uses today's status set; no new status is
  added. An unchanged pin across releases means the DB directory is reused as is.

### 3b. Profile parity: docker and plain apptainer (14p)

The `docker` profile today is only `docker.enabled = true`. It lacks four things that the
`apptainer` profile has or that docker needs, so the core tools cannot run under it. 14p adds them
to `workflow/nextflow.config`:

| Need | `apptainer` today | `docker` after 14p |
|---|---|---|
| DB visible in the task (BAKTA takes `val(db)`) | `runOptions = "-B ${params.db_dir}:${params.db_dir}:ro"` | `docker.runOptions` gains `-v ${params.db_dir}:${params.db_dir}:ro` when `params.db_dir` is set |
| `ps` for the BAKTA/ANTISMASH images (M2-07 decision 13) | `containerOptions = "-B …/assets/ps/ps:/usr/bin/ps:ro"` | `containerOptions = "-v ${projectDir}/assets/ps/ps:/usr/bin/ps:ro"` for `BAKTA\|ANTISMASH` |
| per-task `TMPDIR` | `beforeScript = 'mkdir -p tmp; export TMPDIR=$PWD/tmp'` | the same `beforeScript` |
| files owned by the user | (Apptainer runs as the user) | `docker.runOptions` gains `-u $(id -u):$(id -g)` (M11-09 left it conditional; now always) |
| pinned processes never use conda | `withName: '<pinned set>' { conda = '' }` | the same selector, so `-profile docker,conda` cannot shadow an image |

- The pinned-process selector moves to one shared `process` block, which the `apptainer` and
  `docker` profiles both include, so the two cannot drift. `update_digests.py --nextflow-config`
  keeps editing that one selector.
- **Docker smoke (14p, a CI change).** `docker-smoke.yml` gains one job that runs `BAKTA` on the
  100-kb test genome with a Bakta light DB fixture, or checks that `params.db_dir` is visible and
  readable inside the BAKTA container (`test -r <db_dir>/version.json`) when the full light DB is
  too large for a runner. It is advisory, like the existing jobs.

### 4. Versioning, release candidates and compatibility (14g, 14h; ADR-0060 §3–§6)

**What bumps what** (ADR-0060 §3). The deciding question is whether a user's script or comparison
breaks.

| Change | Bump |
|---|---|
| a column, table, manifest key, flag token or config key is removed, renamed, retyped or changes meaning; a CLI option removed; `output_schema_version` major; a DB pin whose old DB no longer works | **major** (2.0) |
| a tool version or image key change; a DB pin bump; a new column, table, flag, option or module; a default changed (14d2 is one); any change that alters a default run's output bytes on the golden/TE7 runs; a calibration rerun that moves a gate number | **minor** (1.1) |
| a bug fix, crash fix, docs, performance with byte-identical goldens and TE7 | **patch** (1.0.1) |

- **Output schema.** 14h adds `output_schema_version` (tables + manifest), `MAJOR.MINOR`,
  starting at `1.0` with the 1.0.0 release. A new column bumps the minor; a breaking change bumps
  the major. From 1.0.0 on, its major equals the package major. Before 1.0.0 no rule applies
  (no `v*` tag exists). The config's `schema_version` keeps its meaning, the config format
  (Q-L13).
- **Schema diff test (14h).**
  - It compares the column dictionaries in `src/giap/schemas/tables/` and the manifest schema with
    the copies at the newest `v*` tag reachable from HEAD, read with `git show v<prev>:…` (no
    network).
  - It fails when a column was removed or renamed without a major bump of `__version__`, or a
    column was added without a minor bump of `output_schema_version`.
  - **No `v*` tag in the repository** (before the first release): the test skips, with the reason
    `no previous release`.
  - **Shallow clone or no tags fetched** (`git rev-parse --is-shallow-repository` prints `true`): the test **fails** with
    `fetch tags: actions/checkout with fetch-depth: 0` instead of skipping. 14h sets
    `fetch-depth: 0` and `fetch-tags: true` on the CI job that runs it (a CI change, so Khoa
    approves it).

**Release candidates** (ADR-0060 §4). `1.1.0rc1`, `1.1.0rc2`, … are cut on main by the same
release-prep PR as M11-13 step 1. Agents dispatch only `dry_run=true`. Khoa dispatches
`dry_run=false target=testpypi` for an rc; rc versions are uploaded to TestPyPI only, never to
PyPI (Q-K15). The final `1.1.0` follows on the commit of the last rc plus the version line.

**Evidence reuse: one predicate for every reusable gate (14g).** Evidence produced at commit C is
reusable at R iff `reusable(C, R)` holds:

1. Let `D` = `git diff --name-only --no-renames C R`. With `--no-renames`, a move lists both the
   old and the new path, and a deletion lists the deleted path; every listed path counts as
   changed.
2. Every path in `D` lies under an entry of `NON_STATS_PATHS`, matched on **whole path
   components**: an entry `a/b/` matches `a/b/x` but not `a/bx`; a file entry matches only itself.
3. **Version-line exception.** `src/giap/__init__.py` is allowed in `D` iff its diff (`git diff
   -U0 C R -- src/giap/__init__.py`) is exactly one hunk with one `-` and one `+` line, both
   matching `^__version__ = "[^"]+"$`. The version enters no statistic: `rng.for_step(step,
   seed)` does not hash it.
4. **Self-protection.** If `D` contains `scripts/release/readiness.py`, the predicate is false,
   whatever else changed. That file holds `NON_STATS_PATHS` and `CALIB_PATHS`, so a widening
   commit can never certify itself. Changes to the list need `stats` review.

`NON_STATS_PATHS` (a constant in `readiness.py`, 14g):

```text
docs/  changelog.d/  README.md  CHANGELOG.md  mkdocs.yml  PROGRESS.md  AGENTS.md  CLAUDE.md
.claude/  .github/  scripts/dev/  scripts/docs/  src/giap/report/  src/giap/viz/
```

The **never list** is also a constant. A test asserts it is disjoint from `NON_STATS_PATHS`, so a
later edit cannot slip these in:

```text
src/giap/preflight/  src/giap/runner/  src/giap/setupdb/  src/giap/cli.py  src/giap/data/
src/giap/config/  src/giap/rng.py  src/giap/io/  workflow/  pyproject.toml  uv.lock
tests/sim/  tests/calibration/  ci/  scripts/release/readiness.py
```

Why the first four are excluded:

- `preflight/` resolves the abundance scale, the topology, the evidence tier and (with 14d) the
  operon source;
- `runner/` writes `params.json` and `resolved_config.yaml`;
- `setupdb/` holds the DB pins;
- `cli.py` maps `--seed` and `--publication` and runs compare's statistics inline.

`src/giap/report/` and `src/giap/viz/` only render existing values; their regressions belong to
the snapshot and unit tests.

| Gate (M11-03 rows) | patch | minor | major |
|---|---|---|---|
| unit/lint/types, package smoke (CI) | rerun | rerun | rerun |
| `image_key.py --check`, H4, static schema diff test | rerun | rerun | rerun |
| calibration (nightly subset + attacks; `calibration_full`) | reuse iff `reusable(C, R)`, else rerun (subset; full if `stats`-labelled code, `calibration.yaml` or a tool feeding a statistic changed) | same rule | full rerun |
| weekly e2e (K-12, Nissle, TE8) | rerun one (K-12 standard) | rerun | rerun |
| A2 pilot gate (TE1, mapping FDR) | reuse iff `reusable(C, R)`, else rerun | same rule | rerun |
| benchmark `k12_standard` window | reuse iff `reusable(C, R)`, else rerun | same rule | rerun |
| Khoa approval (§14.8 item 5) | yes | yes | yes |

Readiness prints `reused from <C7> (reusable: diff ⊆ NON_STATS_PATHS)` in each reused row.

**Compatibility check (14h).** The manifest already records `giap_version`, each tool's `version`
and `container_digest` (sha256 only, `run_manifest.schema.json`) and the databases. 14d adds
`operon_source`. `<pkg> compare` adds `compare_versions` to `CMP/report/compare.html` and to the
event log:

| Runs differ in | Effect |
|---|---|
| package patch only | nothing |
| package minor | warning `compare.cross_version: runs from 1.0.2 and 1.1.0` |
| any tool `version` or `container_digest` (for a local UniOP build: `source_sha256`, key or runtime instead, item 1b) | warning `compare.tool_versions: ostir 1.1.3 sha256:3f9c0a1b2d4e vs 1.1.4 sha256:77a0c1d2e3f4` (version + the first 12 hex of the digest; one line per tool). The manifest has no tag; `images.json` maps the digest to its tag |
| database version | warning `compare.db_versions: bakta_db 6.0 vs 6.1` |
| `operon_source` | warning `compare.operon_source: …` |
| `output_schema_version` major | error, exit 1: `compare: output schema 1.x and 2.x differ; rerun one with the same major` (override `--allow-cross-major`, which keeps the warning) |

### 5. The release flow (14b, 14g, 14j; ADR-0060 §5)

```text
feature PR edits workflow/containers/<t>/ or envs/<t>.yml (or locks, or giap_fgsea.R)
   └─ PR CI: image_key.py --plan prints "will rebuild: <t> (<old key> → <new key>)"; no push
merge to main
   └─ containers.yml (push trigger; global concurrency group): plan, authenticated check of each key's tag
        ├─ vars.GHCR_PUBLISH != 'true'  → stop after plan ("publishing off"); nothing is pushed
        └─ build (environment `ghcr`, Khoa approves; Q-K18) → re-check → push <version>-<key>
             └─ digest PR ci/digests-<sha7>: containers.config + images.json (+ the shared selector)
                  └─ FarmShare Apptainer smoke (M11-06 item 5) → job IDs in the PR → merge (waitmerge)
release-prep PR (version, changelog, CITATION) → evidence by the table above → readiness
   └─ release.yml dry_run=true (lead): verify adds `image_key.py --check` and an anonymous
      `imagetools inspect` of every GHCR pin (fails `not public: <ref>`)
        └─ Khoa: approval → testpypi (rc) → pypi (final) → tag v<v> by the workflow
```

- **Agents never:**
  - push tags;
  - dispatch `release.yml` with `dry_run=false`;
  - create or set `GHCR_PUBLISH`;
  - create or approve the `ghcr` environment;
  - change package visibility;
  - delete packages.
- **CI changes need Khoa's approval.** 14b (`containers.yml`), 14g (`release.yml` verify), 14p
  (`docker-smoke.yml`) and 14h (one checkout option in `ci.yml`) edit `.github/workflows`. Each of
  these PRs carries the label `needs-khoa` and is merged only after Khoa approves it.
- **14j: Khoa's checklist.** It starts after the rename is merged and 14b is merged.
  1. Create the GitHub environment `ghcr`: required reviewer `khoa-yelo`, deployment branch `main`
     only.
  2. Create the repository variable `GHCR_PUBLISH` = `true`.
  3. Dispatch `containers.yml` with `tools=all`, and approve the `ghcr` deployment.
  4. Make each package public (Q-K8).

  The agents then merge the digest PR and run AC 13. Until step 2, `containers.yml` runs its plan
  job only.
- The wheel ships `workflow/conf/containers.config` and `images.json` (M11-01 bundles
  `workflow/`), so each PyPI release pins exact digests forever.

### 6. Sequencing with the rebrand

- The code PRs land in either order relative to the rename PR: they are code only, and nothing is
  published. The rename script also renames the image prefix (`giap-` → `gippa-`) in
  `containers.yml`, `update_digests.py`, the Dockerfile labels and the `IMAGE_NS`-derived names.
- `containers.yml` derives the prefix from `pyproject.toml` `[project] name`. A static test asserts
  that every GHCR ref in `containers.config` and `images.json` starts with
  `ghcr.io/<owner>/<project name>-`. So nothing can be published as `giap-*` after the rename, and
  no `giap-*` pin can survive it.
- 14j (the first publish) is gated on "rename merged" and on Khoa's steps 1–2. No image is ever
  published under `giap-`.

### 7. Docs (14i)

- **`docs/get-started/install.md`** is rewritten for the five-step path above:
  - a prerequisites table (Python, Nextflow, Java, Docker or Apptainer);
  - `pip install`, `setup --images`, preflight, run;
  - FarmShare from a checkout as a second section;
  - the "Supported platform" warning replaced by the tested-platforms table (M11-09);
  - UniOP: `setup --images` builds it from the authors' archive (14u1, 14u2); the terms notice and
    `--accept-uniop-terms` for scripted installs; what to do when the build fails (offline, no
    fakeroot): rerun elsewhere, or `--set annotation.uniop.enabled=auto` (the default after 14d2)
    with what the distance model loses;
  - on a Slurm executor, every node needs the same Python environment, because `GIAP_STEP` runs
    on the host (Q-L8: no core image).
- **`docs/get-started/upgrade.md`** is new, placed after install in the nav. It covers:
  - `pip install -U <pkg>==1.y`;
  - which images change: the `images.json` diff, shown by `<pkg> setup --images --check`;
  - `setup --check` for DB pins;
  - comparing runs across versions (the compare warnings);
  - the bump table from ADR-0060;
  - rollback: `pip install <pkg>==1.x`, whose own pins still resolve (item 2).
- **`docs/release.md`:**
  - step 5 is rewritten for 14j's checklist;
  - a new "Cutting 1.x.y" section holds the gate table and the reuse predicate;
  - the backlog entry "Off-FarmShare install" is closed.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `image_key.py ostir` on a tree where `workflow/containers/ostir/Dockerfile`, `workflow/envs/ostir.yml`, `locks/ostir.explicit.txt`, `locks/ostir.pip.txt` exist | 12 hex key over 4 code-point-sorted `path\tmode\tsha256` lines; tag `1.1.3-<key>` |
| 2 | the same tree with one byte changed in `ostir.yml` | different key; `--check` against the old `images.json` exits 1 `stale image: ostir (inputs changed since <old key>)`; other tools unaffected |
| 2b | the same tree with only the exec bit of the Dockerfile flipped | different key |
| 3 | a comment added to `workflow/envs/bakta.yml` (upstream-image tool) | no GHCR key changes; `--plan` prints `will rebuild: none` |
| 4 | `promotech` registry version `git:56251ad9b883…` | tag `git-56251ad9-<key>` |
| 5 | plan where `inspect` of `ghcr.io/khoa-yelo/gippa-ostir:1.1.3-<key>` exits 0 (fake script) | `ostir` not built; its digest goes into `digests.json` |
| 5b | `inspect` fails with `unauthorized` | the plan job fails; nothing is built |
| 5c | `inspect` fails with `manifest unknown` at plan time, and exits 0 at the pre-push re-check | build job fails `tag exists: <ref>`; no push |
| 6 | `missing_tools([ostir], "docker", …)` with an `OSTIR` pin and no env | `[]` |
| 7 | `missing_tools([ostir], "apptainer", cache=fake)` with the pin and no `.img` in the fake cache | `["ostir"]`; message names `<pkg> setup --images --profile apptainer` |
| 7b | a tool with only a prebuilt env, profile `apptainer` (no conda profile) | missing: `conda-only tool under a container profile without conda` |
| 8 | preflight, `enabled: auto`, profile `docker`, no UniOP, no `--publication` | exit 0; warning `uniop not installed …`; `params.run_uniop = false`; after the run `layers_missing` contains `uniop`, the flag `operon_source_distance_model` is set, manifest `operon_source: distance_model`, `operon_source_reason: UniOP not installed` |
| 8b | the same with `--publication` | exit 2, `uniop not installed; set annotation.uniop.enabled: false to publish …` |
| 9 | `enabled: true`, no UniOP | exit 2, `environment.tools: uniop` missing |
| 9b | `enabled: false` with `--publication` | exit 0; reason `UniOP disabled (annotation.uniop.enabled=false)`; flag and report caveat present |
| 10 | `compare` of a 1.0.0 run and a 1.1.0 run whose OSTIR digests differ | exit 0; warnings `compare.cross_version` and `compare.tool_versions: ostir 1.1.3 sha256:<12> vs …` |
| 11 | `compare` of `output_schema_version` 1.3 and 2.0 runs | exit 1 unless `--allow-cross-major` |
| 12 | readiness for 1.0.1 where `D` = {`src/giap/viz/locus.py`, `docs/x.md`, `src/giap/__init__.py` (version line only)} | calibration, pilot and benchmark rows `pass (reused from <C7> …)` |
| 13 | the same with `src/giap/cli.py` (or `src/giap/stats/fdr.py`) in `D` | those rows `pending` until rerun on R |
| 13b | `D` contains `scripts/release/readiness.py` | not reusable, even if only `NON_STATS_PATHS` changed |
| 13c | a move `src/giap/stats/x.py` → `docs/x.py` | `D` lists both paths (`--no-renames`); not reusable |
| 13d | `D` = {`src/giap/report_stats.py`} | not reusable (component match: `src/giap/report/` does not match) |
| 14 | `update_digests.py` with `{"UNIOP": "ghcr.io/…/gippa-uniop@sha256:…"}` | exit 1 (M11-06 example 5 still holds) |
| 15 | `fgsea_r.command` under `GIAP_CONTAINER_RUNTIME=docker`, no prefix, no env | the exact docker argv of item 1, `--cidfile <call_dir>.cid`, `-u <uid>:<gid>`, `Rscript --vanilla /opt/<pkg>/giap_fgsea.R` |
| 15b | `GIAP_FGSEA_PREFIX` set to a directory without R, image available | `FgseaEngineUnavailable` (no fall-through) |
| 16 | `setup --images --profile docker --accept-uniop-terms`, fake fetcher serving the fixture archive (sha256 equal to the pin) behind a fake `302` from `github.com` to `codeload.github.com`, fake `docker` | one fetch of the pinned URL, the redirect followed; archive at `<db_dir>/uniop/UniOP-<commit>.tar.gz` (mode `0600`); `docker build … -t localhost/<pkg>-uniop:git-3ed14dfc-<key> <context>`; `local_tools.config` has only the `docker` block, with `conda = ''`, and no `conda.enabled`; `image:uniop` `ok` with `local_build.runtime: docker`, `digest` = the fake image ID, `terms: flag`; no `push` |
| 16a | a redirect to any other host | `refused …` (the allowlist holds only `codeload.github.com` for `github.com`) |
| 16b | the same, profile `apptainer`, not root, fake `apptainer` whose `--fakeroot` build exits 1 with `GLIBC_2.38 not found` and whose `--fakeroot --ignore-fakeroot-command` build succeeds | exactly two build calls, in that order; the `apptainer` block names the absolute `.sif` path, with `conda = ''`; `digest` = the `.sif` sha256 |
| 16c | the same, profile `conda` | env at `<cache-dir>/envs/giap-uniop`; `<prefix>/share/uniop/GIAP_COMMIT` = the commit; the `conda` block `conda = '<prefix>'`; no `conda.enabled` |
| 16d | `setup --uniop --profile farmshare,apptainer --accept-uniop-terms`, a fake `$GIAP_ROOT` with `micromamba/envs/giap-uniop/share/uniop` | `R = micromamba`; no fetch, no build, no `local_tools.config` change; `local_build.runtime: prebuilt`, `ref` = that prefix. Without the prefix: exit 2 naming `sbatch scripts/slurm/build_envs.sbatch uniop` |
| 16e | `R(P)` for `docker`, `apptainer`, `conda`, `farmshare`, `farmshare_local,apptainer`, `apptainer,docker`, `conda,docker`, `standard` | docker, apptainer, micromamba, micromamba, micromamba, apptainer, ambiguous (missing), none (missing) |
| 17 | the fixture archive with one byte changed | exit 2 `uniop: archive sha256 <got> != pinned <want>`; the file is deleted; no build command runs; `local_tools.config` unchanged |
| 18 | `setup --images --profile docker`, no TTY, no `--accept-uniop-terms` | the notice is printed; no fetch and no build; `image:uniop` `absent`, message `terms not accepted; rerun with --accept-uniop-terms`; other images pulled; exit 0 with a warning |
| 18b | `setup --uniop --profile docker`, no TTY, no flag | exit 2 with the same message |
| 18c | TTY, the user answers `n` (or presses Enter) | as 18 (`--images`) or 18b (`--uniop`); answering `y` builds, with `terms: prompt` |
| 19 | preflight, default `enabled: true`, profile `docker`, after example 16 | exit 0; `uniop` present; `giap run`'s command line has `-c <db_dir>/local_tools.config` |
| 19a | `giap run --profile farmshare,apptainer` with a `local_tools.config` present | no `-c <db_dir>/local_tools.config` on the Nextflow command line |
| 19b | preflight, `enabled: true`, `--db-dir /data/db`, profile `apptainer`, nothing built | exit 2 `environment.tools: uniop not installed; run: <pkg> setup --uniop --profile apptainer --db-dir /data/db --accept-uniop-terms` |
| 19c | after example 16, the local tag removed, its image ID differing from `local_build.digest`, or `local_build.key` differing from `key("uniop")` | `uniop` missing, the same command (stale: `local build <old> does not match this release (<new>)`) |
| 19d | after example 16 (docker only), preflight under `-profile apptainer` | `uniop` missing, the apptainer command; nothing tries `docker://localhost/…` |
| 20 | `setup --images --profile apptainer --accept-uniop-terms`, both fake apptainer builds exit 1 | the other images are installed first; `image:uniop` `absent` with `apptainer build failed (fakeroot unavailable or build error): <last stderr line>; use --profile docker or conda, or set annotation.uniop.enabled=auto`; setup exits 2; preflight under `enabled: auto` then gives example 8 |
| 21 | stub run (integration) after example 16, profile `test,docker,stub` | manifest `operon_source: uniop`; the `uniop` entry has `local_build: true`, `runtime: docker`, `source_sha256`, `container_digest` = `sha256:<fake image ID>` |

## Edge cases

- **Two processes, one image** (`PROMOTECH`/`PROMOTECH_NOVEL`, `OSTIR`/`OSTIR_NOVEL`,
  `BLAST`/`TBLASTN`): one key, one digest, both lines; `images.json` lists both processes.
- **An image with no process** (fgsea): only in `images.json`; H4 accepts it for `fgsea` only.
- **Key unchanged, upstream source changed.** This cannot happen:
  - for the four git tools, the archive URL names the pinned commit and `ADD --checksum` pins its
    bytes (14b);
  - for pip and conda packages, the locks are in the key;
  - the base image is pinned by digest inside the Dockerfile.
- **Rebuilding the same key gives a different digest** (non-reproducible layers, timestamps): the
  tag already exists, so no rebuild happens; the first digest wins forever.
- **A tool fails to build** (Promotech's Python 3.6 env, §17): that tool keeps its previous pin,
  or stays conda-only on FarmShare. The digest PR says so, `--check` keeps failing until it builds,
  and readiness H4 is `fail` (no release with a stale pin).
- **Private package** (Q-K8 not done): release.yml verify fails `not public: <ref>`, and the
  FarmShare smoke reports `private package`. The authenticated plan check is unaffected.
- **Offline compute nodes** (HPC): `setup --images` runs where the network is, and preflight checks
  the cache, not the registry.
- **Docker on arm64 (Apple Silicon)**: images are `linux/amd64` only, so Docker runs them emulated
  and slowly; the platforms page says so (M11-06 question, unchanged).
- **A user with an old `$GIAP_ROOT/micromamba` env and a new image**: under the container profiles
  the image wins (item 3 table); under `farmshare*` the env wins, as today.
- **Empty input** (no enabled tool needs an image, e.g. every layer precomputed):
  `setup --images` pulls nothing and records nothing, and preflight passes.
- **Precomputed operon pairs** (`annotation.precomputed.operon_pairs`): UniOP is not required and
  `enabled` is ignored. `operon_source` is `precomputed` (as today). `setup --images` does not
  build UniOP for that config.
- **UniOP: offline or no build rights** (no network, no docker daemon access, no fakeroot): the
  item is `absent` with the reason, and setup exits 2 (item 1b). The documented way out is to
  build on a machine that can (docker: `docker save`/`load`; apptainer: copy the `.sif`, then
  rerun `setup --uniop` there to record it), or `enabled: auto`/`false` (item 1 (c)).
- **UniOP: the build host is not the run host** (docker on a laptop, run on a server; an
  apptainer `.sif` on a path compute nodes cannot see): preflight checks the run host, so the
  image-ID or `.sif` check fails there and names the command.
- **UniOP: package upgrade.** A new release whose uniop key differs (Dockerfile, `.def`, env,
  locks or archive pin changed) reports the old build as stale (example 19c). A release with the
  same key reuses the build.
- **UniOP: GitHub regenerates the archive.** The bytes of a commit archive can change while the
  commit stays the same. The sha256 check then fails (example 17), and nothing is built from
  unverified bytes. The fix is a reviewed pin update (R6 official source), not a bypass flag.
- **UniOP: podman.** Out of scope for 1.0: only the docker CLI is supported. On FarmShare
  (podman, no docker) the routes are the prebuilt env (`farmshare*`) or apptainer.
- **UniOP: shared caches.** A lab-wide `NXF_APPTAINER_CACHEDIR` never exposes one user's archive
  or `.sif`, because of the `0700`/`0600` modes (item 1b).
- **UniOP: rerun after success.** `setup --uniop` with an `ok` item whose key and digest still
  match does nothing (idempotent); `--check` verifies without building.
- **GHCR outage or rate limit**: the plan job fails (B3 classification). R6: retry once later, then
  PROGRESS › Blocked.
- Genomic edge cases: circular replicons matter only through the operon fallback, where
  `distance_model_p` uses the wrap-around intergenic distance as today (preflight topology,
  unchanged). Minus-strand pairs are same-strand pairs like any other. Ties at `p_dist = 0.5`
  count as linked (`≥`, as in ADR-0022). Linear/circular, strand and tie behaviour are otherwise
  not applicable to this issue.

## Acceptance criteria

All tests are offline. No test pulls from or queries a registry: `imagetools inspect`,
`docker pull/build/image inspect`, `apptainer pull/build` are replaced by fakes on `PATH` or by
injected callables. Preflight tests use temporary fake `GIAP_ROOT`/cache/DB directories.

1. `image_key.py`: examples 1–4 and 2b on a temporary tree. Keys are stable across runs and
   independent of mtime and directory-listing order. A Dockerfile `COPY` of a file outside the
   inputs fails the static test, and so does a Dockerfile `GIT_SHA` or env pin that differs from
   `annotate/tools.py`.
2. `containers.yml` static test:
   - only the `<version>-<key>` tag is pushed;
   - the plan job has `packages: read` and logs in;
   - the classification follows examples 5, 5b and 5c (a fake `docker` script driven by the
     workflow's own shell steps, extracted and run by the test);
   - the build job re-checks before the push;
   - the concurrency group is the global `containers`;
   - `build` needs `vars.GHCR_PUBLISH == 'true'` and the `ghcr` environment;
   - the prefix is derived from `pyproject.toml`;
   - no `uniop`; every third-party action is SHA-pinned; nothing pushes to main;
   - every `ADD` of a URL has `--checksum`.
3. `update_digests.py` writes `images.json` and `containers.config` consistently, including the
   process-less fgsea entry. It is idempotent, and M11-06 examples 1–5 still pass.
   `image_key.py --check` passes on the repo after the digest PR.
4. A static test: every GHCR ref starts with `ghcr.io/<owner>/<project name>-`, and every
   non-in-process tool in `TOOLS` other than `uniop` has a pin (`withName`, or `images.json` for
   fgsea). Until 14j merges, the test lists the expected-missing set explicitly, and that set may
   only shrink.
5. Preflight: one test per cell of the item 3 table on fake roots; examples 6, 7 and 7b; the item 1
   resolution table, including examples 8, 8b, 9 and 9b.
6. `setup --images`, with fake `apptainer`/`docker` binaries:
   - the file names equal Nextflow's `simpleName` (ADR-0059);
   - a digest mismatch fails `pulled <d1>, pinned <d2>`;
   - the manifest items are written;
   - the image set follows `--config`/`--set` and `--images all`;
   - `setup --check` reports a DB pin change as `mismatch` with both versions;
   - `ToolRunner` works with no `GIAP_ROOT` and no checkout, through the docker path.
7. fgsea:
   - the `resolve_launcher` order (prefix → env → image), with an explicit prefix never falling
     through (examples 15 and 15b);
   - the docker and apptainer argv are exact (item 1);
   - `workflow/bin/giap_fgsea.R` is a key input of `fgsea`;
   - the GSEA golden (`golden-r.yml`) is unchanged when run through the image on FarmShare (one
     Slurm job, ID in the PR).
8. UniOP degrade (14d):
   - the stub-run integration test with `run_uniop=false` produces `layers_missing: [uniop]`, the
     `operon_source_distance_model` flag, manifest `operon_source`/`operon_source_reason`, and a
     distance-model `operon_pairs.tsv`;
   - each of the three reasons appears in its case;
   - the report fixture shows the operon-source text with the numbers from `calibration.yaml`.
9. UniOP local build (14u1, 14u2), with an injected fake fetcher and fake
   `docker`/`apptainer`/`micromamba` on `PATH`. The fixture archive is a few-byte tarball made by
   the test, never UniOP code (R5).
   - **No network, enforced.** `tests/conftest.py` has no network guard today, so 14u1 adds a
     non-autouse `no_network` fixture to `tests/unit/conftest.py`. It monkeypatches
     `socket.socket.connect`, `socket.create_connection` and `socket.getaddrinfo` to raise. Every
     14u test module sets `pytestmark = pytest.mark.usefixtures("no_network")`, and one test shows
     that a real `HttpsFetcher` fetch raises under it.
   - **14u1:**
     - examples 16, 16a, 16e, 17 and 19a;
     - the `uniop_source` row loads through `_load_sources`, and `plan()` leaves it out of plain
       `setup`;
     - the Dockerfile has no `ADD` of a URL;
     - the archive and the `uniop/` directories have modes `0600` and `0700`;
     - `key("uniop")` computed from an installed wheel equals the repo's key (package smoke);
       changing `uniop.def` or the archive pin changes the key;
     - the generated `local_tools.config` never contains `conda.enabled`, every container line
       carries `conda = ''`, and it parses under `nextflow config -profile docker` (a CI job that
       already has Nextflow).
   - **14u2:**
     - examples 16b, 16c, 16d, 18, 18b, 18c, 19, 19b, 19c, 19d and 20;
     - `uniop.def` has no `chown` and no `curl`/`wget`/URL;
     - the notice is printed on every path, and the prompt default is no;
     - the `local_build` object and the preflight rows (item 3) agree;
     - no `push` is ever called; `update_digests.py` still refuses `UNIOP` (example 14).
   - a micromamba-only build under `-profile docker` is missing (no `docker` block), with the docker command.
9a. **UniOP integration (14u2):** example 21, as a stub-run test in `tests/integration` with a
    fake `docker`.
9b. **One manual Slurm job on FarmShare** (ID in the 14u2 PR, not a CI test). It runs the real
    `setup --uniop --profile apptainer --accept-uniop-terms`, which takes the
    `--ignore-fakeroot-command` fallback, then `--profile farmshare,apptainer` (prebuilt), and
    `preflight` for each. It records the archive sha256 that 14u1's pin carries. The 14u1
    implementer computes that value first, in a short job.
10. **14s1/14s2 (statistician):**
    - arms A and B run at the existing bounds; the outcome is recorded in the PR;
    - `calibration.yaml: operon_distance` carries the precision and recall at 0.5 with CIs;
    - the ADR-0022 addendum is merged;
    - the weekly arm writes its `weekly.tsv` row.

    14d2 (the default flip) has a test asserting that `defaults.yaml` says `auto` only if those
    keys exist in `calibration.yaml`.
11. Profile parity (14p): a static test parses `nextflow.config` and asserts, for `docker` and
    `apptainer`:
    - the DB bind, the ps shim on `BAKTA|ANTISMASH`, the `TMPDIR` `beforeScript` and the
      pinned-process conda-off selector (shared block);
    - for docker, `-u $(id -u):$(id -g)`.

    An `-stub-run` integration test under `-profile test,docker,stub` with a logging fake `docker`
    sees `-v <db_dir>:<db_dir>:ro` and `-u` on the BAKTA call. The docker-smoke job's DB-visibility
    check is green once (run URL in the PR).
12. Compare: examples 10 and 11 on two fake run directories whose manifests differ only in the
    listed fields, including the `operon_source` warning.
13. Readiness (14g):
    - examples 12, 13, 13b, 13c and 13d;
    - the version-line exception accepts exactly one `__version__` hunk and rejects any second
      changed line in `src/giap/__init__.py`;
    - `NON_STATS_PATHS` and the never list are disjoint, and every entry of both exists;
    - H4 fails on a stale key.
14. Schema diff test (14h): a removed column with an unchanged major fails, and an added column
    with an unchanged schema minor fails (fixture trees). No `v*` tag skips with
    `no previous release`; a shallow fixture repo fails with `fetch tags`.
15. Docs: `mkdocs build --strict -d <scratch>` passes with `upgrade.md` in the nav, and the
    install page's commands match the CLI options (the generated CLI reference, M11-07).
16. After 14j (Khoa), all in one Slurm job, job ID in PROGRESS: on a FarmShare node with only the
    wheel installed in a fresh venv (no checkout on `PYTHONPATH`), `setup --images --profile
    apptainer`, `preflight` and a 100-kb `run` pass.

## Estimate and PR split (each ≤ ~400 lines of code; tests counted separately, ADR-0007)

| PR | Scope | Code | Tests | Depends on | Label |
|---|---|---|---|---|---|
| 14a | `key()` in `src/giap/io/image_key.py` (mode from the exec bit, so a wheel can compute it) and the `scripts/release/image_key.py` CLI over it (tag, `--plan`, `--check`), `images.json` (process-less entries), `update_digests.py` writes it, `image_cache()` helper | 190 | 220 | M11-06 | |
| 14b | `containers.yml`: key tags only, authenticated check + re-check, global concurrency, `GHCR_PUBLISH` + `ghcr` gate, prefix from pyproject; Dockerfiles from locks with `ADD --checksum` (not uniop: 14u1 replaces its `ADD` with `COPY`); BacTermFinder locks; pin-equality test | 160 | 150 | 14a | CI change, **Khoa approval** |
| 14p | profile parity: docker DB bind, ps shim, `TMPDIR`, `-u`, shared pinned-process block; docker-smoke DB-visibility job | 70 | 110 | none | CI change, **Khoa approval** |
| 14f1 | `missing_tools` new signature and the item 3 table (conda-on-host under docker and plain apptainer; local tools) | 110 | 160 | 14a | |
| 14f2 | `setup --images` (config-driven set, digest check, manifest items); `ToolRunner` docker path, no `GIAP_ROOT`/checkout requirement | 170 | 180 | 14a, 14f1 | |
| 14d | `annotation.uniop.enabled` (default `true`), preflight resolution incl. `--publication`, `params.run_uniop`, workflow skip, reasons, `layer_missing:uniop`, manifest `operon_source*`, report text | 170 | 200 | 14f1 (same function: merge 14f1 first) | **`stats`** |
| 14s1 | `tests/calibration/attacks/test_attack_operon_distance_fallback.py`, arms A and B | 0 | 260 | none | **`stats`**, statistician |
| 14s2 | precision/recall at 0.5 with CIs in `calibration.yaml`; ADR-0022 addendum; weekly arm `uniop.enabled=false` | 40 | 40 | none | **`stats`**, statistician |
| 14d2 | default `auto` (`defaults.yaml`, model); guard test | 10 | 20 | 14d, 14s1, 14s2, statistician sign-off | **`stats`** |
| 14u1 | `setup --uniop` and its inclusion in `setup --images` (`--no-uniop`); `R(P)`; `uniop_source` row; `REDIRECT_HOSTS` github → codeload; verified fetch; staging context and permissions; Dockerfile `COPY` + header and `uniop.def`; the archive key line (`key()` in `src/giap/io/image_key.py` from 14a); docker builder; `local_tools.config` profile blocks; `giap run -c` rule; `no_network` fixture (replaces the 14e draft) | 150 | 180 | 14a, 14f2 | R6 allowlist change |
| 14u2 | apptainer builder (fakeroot fallback) and micromamba/prebuilt; terms notice, prompt and `--accept-uniop-terms`; `local_build` item and the run-manifest `uniop` entry; preflight rows and message; example 21 integration test; the AC 9b Slurm job | 130 | 170 | 14u1, 14d | |
| 14c | fgsea Dockerfile with the script baked in; `fgsea_r` image path; `GIAP_CONTAINER_RUNTIME` from `giap run`; NOTICE line | 130 | 150 | 14a, 14f1 (cache helper, runtime) | |
| 14g | readiness: `reusable(C, R)` (components, `--no-renames`, version-line exception, self-protection), `NON_STATS_PATHS` + never list, H4 key check and process-less fgsea; `release.yml` verify: `--check` + anonymous inspect | 140 | 170 | 14a | **`stats`**, CI change, **Khoa approval** |
| 14h | compare version, digest and operon-source warnings; `output_schema_version`; schema diff test; `fetch-depth: 0` in CI | 120 | 190 | 14d (operon_source) | CI change, **Khoa approval** |
| 14i | docs: install, upgrade, release.md, backlog entry closed; changelog fragments | 0 (≈ 380 prose) | 20 | all above | |
| 14j | **Khoa**: `ghcr` environment, `GHCR_PUBLISH`, dispatch `tools=all`, approve, make public; then the digest PR and AC 16 | 0 | 0 | rename, 14b | Khoa |

Total ≈ 1,590 code, ≈ 2,220 test lines, ≈ 380 prose.

**Waves:**
- W1: 14a, 14p, 14f1, 14s1, 14s2.
- W2: 14b, 14f2, 14d, 14g, 14c.
- W3: 14u1 then 14u2, and 14h.
- W4: 14d2, after the statistician's sign-off, and 14i.
- Then the rename, then 14j.

14d merges after 14f1, because both edit `missing_tools`. 14c needs the cache helper (14a) and the
runtime plumbing (14f1).

## Questions

Numbering continues M11-index (Q-K12, Q-L10); the rows are copied there. Each has a reversible
default that the PRs use.

| # | Question | Default in use | PR |
|---|---|---|---|
| Q-K13 | UniOP for pip users: (a) ask the authors for a licence, (b) local build, (c) degrade to the distance model when absent? And under `--publication`? | **answered by Khoa 2026-10-25 (install mechanism):** (b) is the install path: pip install → `setup` (incl. UniOP, fetched from the authors with a sha256 check, built locally after the user accepts the terms) → preflight passes with UniOP (14u1, 14u2). (c) `enabled: auto` is kept only for failed builds, still the default only after 14s1/14s2 and the statistician (14d2). `--publication` unchanged: `auto` must resolve to `true` or exit 2; an explicit `enabled: false` is accepted, with the flag and the caveat. (a) optional: Khoa may ask the authors; if granted, `<pkg>-uniop` is published like the others and the local build becomes a fallback | 14u1, 14u2, 14d, 14d2 |
| Q-K14 | GHCR owner after the rename: `ghcr.io/khoa-yelo/gippa-<tool>` or an organisation? | `khoa-yelo` (same as the repo owner today) | 14b, 14j |
| Q-K15 | Release candidates: TestPyPI only, or also PyPI pre-releases (`pip install --pre`)? | TestPyPI only; PyPI gets finals only (M11 Q-K4 for 1.0, extended to every release) | 14g, 14i |
| Q-K16 | GHCR retention: prune unreleased digests after 180 days, or never prune? | never prune (public storage is free); the 180-day rule applies only if Khoa chooses to prune | ADR-0060 |
| Q-K17 | Nextflow and Java are not pip-installable. Accept them as prerequisites (step 2), or add `<pkg> setup --nextflow` that downloads the pinned launcher from the official release (still needs Java)? | prerequisites, documented with a micromamba one-liner; no download command in 1.0 | 14i |
| Q-K18 | Every `containers.yml` build needs Khoa's approval of the `ghcr` environment, or only the first one (Khoa then removes the reviewer)? | every build, until Khoa removes the reviewer (it is a publish) | 14b, 14j |
| Q-K19 | fgsea in an image (14c), or keep M8-17 decision 9 (conda-only; GSEA `layers_missing` for pip users)? | image (14c); without it, `--abundance` runs break "preflight passes" | 14c |
| Q-L11 | Image cache without `$GIAP_ROOT` | `NXF_APPTAINER_CACHEDIR` → `$GIAP_ROOT/containers` → `<db_dir>/containers` | 14a, 14f1 |
| Q-L12 | BacTermFinder has no lock files: resolve one on FarmShare? | yes (14b), with the build log in the PR | 14b |
| Q-L13 | `schema_version` today is the config's field copied into the manifest. Split out `output_schema_version`? | yes, starting at `1.0` with 1.0.0; the config's `schema_version` keeps its meaning (config format) | 14h |
| Q-L14 | Keep the M11-06 `:<sha7>` tags for traceability? | no: `images.json` `built_from` and the OCI `revision` label record the commit | 14b |
| Q-L15 | `NON_STATS_PATHS` contents and the never list | as in item 4 after the statistician's review of #367; changes need `stats` review | 14g |

### Gaps found while writing (for the lead)

- **Conda-only processes run on the host.** Under `-profile docker` alone **and under plain
  `-profile apptainer`**, a conda-only process has its `conda` directive ignored: conda is enabled
  only by the `conda` profile and by `prebuilt_envs.config`. Such a process runs on the host
  without its tool, while preflight reports it present when the env exists. Item 3 and 14f1 make
  it a `missing` with a reason.
- **The `docker` profile cannot run Bakta today.** It has no DB bind, no ps shim, no `TMPDIR` and
  no `-u`, and `docker-smoke.yml` exercises only MINIPROT and TBLASTN. Item 3b and 14p fix this.
- **Image builds are not reproducible.** The M11-06 Dockerfiles solve `workflow/envs/<tool>.yml`
  at build time instead of installing the locks, and their `ADD`s carry no checksum. 14b fixes
  both before the first publish.
- **The Python environment must match on every node.** `GIAP_STEP` runs on the host with the
  host's `giap`, so on a Slurm executor every node needs the same environment. The install page
  says so.
- **UniOP archive sha256 not yet known** (14u1). The pin needs the sha256 of the official commit
  archive, computed once on FarmShare by the 14u1 implementer (the Slurm job of AC 9b). GitHub
  does not strictly guarantee byte-stable archives; if one changes, the check fails closed
  (edge cases). The alternative, a `git clone` at the commit with a tree-hash check, needs `git`
  on the user's machine, so it is not the default.
- **Non-interactive default for the terms** (14u2): no, unless `--accept-uniop-terms` or a TTY
  "y". The lead may prefer yes-with-notice; that is a one-line change to the default, but it
  means GIAP accepts on the user's behalf, which is why it is not the default here.
- **The distance-model fallback has not been tested as a default.** No calibration test or attack
  runs the TU-block nulls with `source = distance_model`. giapsim's gap layout never makes the
  model split a true TU, and `pair_threshold = 0.5` was calibrated on UniOP scores. 14s1 and 14s2
  close this before 14d2.
