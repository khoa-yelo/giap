# M2-01: Tool registry, environment files, container pins, `build_envs`

Milestone M2 · depends on: M0 · estimate ≈ 380 code lines (tests not counted, ADR-0007) ·
capabilities: unit tests none; env builds `needs_slurm`, `needs_conda`, `needs_net`; image pulls `needs_apptainer`

## Spec sections implemented

- §7.1.1 tool table (pinned): Bakta 1.12.2 "conda `bakta=1.12.2` or container `oschwengers/bakta`
  (pin by digest; Apptainer on FarmShare)"; UniOP "conda env with Python 3.11, pandas, numpy,
  scikit-learn; source cloned at a pinned commit when the env or image is built"; Promotech "frozen
  env: Python 3.6.10, TensorFlow 2.1.0, scikit-learn 0.23.1"; Promoter Calculator 1.2.4 "(PyPI
  `promotercalculator`, GPL) own env"; TerminatorNet "own env; scikit-learn pinned to the version that
  pickled `model.pkl` (G12)"; TransTermHP "conda `transtermhp=2.09`"; BacTermFinder "pinned commit, own
  env"; OSTIR "own env with ViennaRNA (conda-forge/bioconda) + `pip install ostir==1.1.3`"; codon-bias
  0.5.0 "main GIAP env, in-process"; geNOMAD "conda `genomad=1.12.0`"; antiSMASH "container
  `antismash/standalone` (pin by digest …) or conda `antismash=8.0.4`"; miniprot 0.18, BLAST+ pinned.
- §4: "`workflow/envs/<tool>.yml` conda env per tool (conda-forge + bioconda, exact pins)";
  "`workflow/containers/<tool>/` Dockerfile per tool (built in CI, pushed to GHCR, referenced by digest)".
- §6.3: each process declares "`conda` (`workflow/envs/<tool>.yml`) and `container`
  (`ghcr.io/<owner>/giap-<tool>@sha256:…`)"; §3.1 "Images are referenced by digest and run as
  `docker://…@sha256:…` through Apptainer".
- §3.2: conda only from conda-forge/bioconda, "never the `defaults` channel"; §3.9 "run `micromamba
  clean -a` after builds"; §12 `containers.yml` "open a PR updating the digests in
  `workflow/conf/containers.config`"; R5 licences; §2.1 decision 18 (UniOP never vendored).
- AGENTS.md: `sbatch scripts/slurm/build_envs.sbatch <tool> [...]` "(micromamba envs
  (workflow/envs/*.yml) and Apptainer pulls)".

## Public interface

`src/giap/annotate/tools.py`
```python
Module = Literal["core", "mobile", "bgc", "function", "opt_in", "mapping"]

@dataclass(frozen=True, slots=True)
class ToolSpec:
    name: str                  # registry key = env suffix = workflow/modules/local/<name>
    version: str               # "1.12.2", or "git:<40-hex>" for git-pinned tools
    licence: str               # "GPL-3.0", "AGPL-3.0", "MIT", "unlicensed", "non-commercial", ...
    env_file: str              # "workflow/envs/<name>.yml"
    env_name: str              # "giap-<name>"
    version_cmd: tuple[str, ...]
    module: Module
    source_repo: str | None = None   # git URL for git-pinned tools
    in_process: bool = False         # codon-bias only

TOOLS: Mapping[str, ToolSpec]       # bakta uniop promotech promoter_calculator terminatornet transtermhp
                                    # bactermfinder ostir genomad antismash miniprot blast codon_bias
def required_tools(cfg: GiapConfig) -> list[ToolSpec]      # for T-PRE-23 (M1-11) and giap setup

@dataclass(frozen=True, slots=True)
class EnvFile: name: str; channels: tuple[str, ...]; conda: tuple[str, ...]; pip: tuple[str, ...]
def read_env_file(path: Path) -> EnvFile
def env_problems(env: EnvFile) -> list[str]                # [] when pinned and channel-clean
def read_containers_config(path: Path) -> dict[str, str]   # process name → container string
```
Files: `workflow/envs/{uniop,promotech,promoter_calculator,terminatornet,transtermhp,bactermfinder,
ostir,genomad,antismash,miniprot,blast}.yml` (bakta.yml exists); `workflow/conf/containers.config`;
`workflow/containers/{uniop,promotech,promoter_calculator,terminatornet,bactermfinder,ostir}/Dockerfile`;
`scripts/build_envs.sh <tool>...` (idempotent: create if absent, clone git sources at the pinned commit
into `$MAMBA_ROOT_PREFIX/envs/giap-<tool>/share/<tool>/`, run `version_cmd`, write
`workflow/envs/locks/<tool>.explicit.txt` via `micromamba env export --explicit`, `micromamba clean
-a -y`); `scripts/slurm/build_envs.sbatch` (§3.8 header, `--cpus-per-task=4 --mem=16G --time=08:00:00`,
also `apptainer pull` of every pinned image into `$NXF_APPTAINER_CACHEDIR`); `scripts/resolve_digest.py
<image:tag>` (registry v2 API, prints `docker://image@sha256:…`); rows in `NOTICE.md`.

`env_problems` rules: channels exactly `[conda-forge, bioconda]` (order free, `defaults` or
`nodefaults`-less extra channels rejected); every conda dependency `name=version` or
`name=version=build`; every pip dependency `name==version`; `name` = `giap-<file stem>`.

## Examples

| # | Input | Expected output |
|---|---|---|
| 1 | `read_env_file("workflow/envs/bakta.yml")` | `EnvFile("giap-bakta", ("conda-forge","bioconda"), ("bakta=1.12.2",), ())`; `env_problems` → `[]` |
| 2 | env with `- pandas` | `["pandas: conda dependency without exact version"]` |
| 3 | env with `channels: [defaults, bioconda]` | contains `"channel 'defaults' is not allowed"` |
| 4 | env with `pip: [ostir>=1.1]` | `["ostir>=1.1: pip dependency must be pinned with =="]` |
| 5 | `required_tools(default config)` | names `bakta, uniop, promotech, promoter_calculator, terminatornet, transtermhp, ostir, genomad, antismash` (bgc and mobile on by default; `bactermfinder` off; `codon_bias` in-process; miniprot/blast listed under `mapping`, returned too) |
| 6 | `--set annotation.precomputed.bakta_json=x.json --set modules=[core]` | no `bakta`, `genomad`, `antismash` |
| 7 | `--set annotation.terminators.cross_check_transtermhp=false` | no `transtermhp` |
| 8 | `containers.config` line `withName: 'BAKTA' { container = 'docker://oschwengers/bakta@sha256:<64hex>' }` | `read_containers_config` → `{"BAKTA": "docker://oschwengers/bakta@sha256:<64hex>"}` |
| 9 | `bash scripts/build_envs.sh bakta` twice | second run creates nothing, prints `giap-bakta exists`, exit 0 |

## Edge cases

- A container string with a tag but no digest (`…:v1.12.2`) is rejected by the test (digests only).
- Tools whose GHCR image does not exist yet have **no** `withName` entry; the list of such tools is
  explicit in `containers.config` comments and in a test constant (see Q3).
- `version` for a git tool must be a full 40-hex SHA, never a branch or tag.
- `build_envs.sh` with an unknown tool name exits 2 and lists valid names.
- TerminatorNet's scikit-learn pin is found by unpickling `model.pkl` once with a recent sklearn
  and reading the "trained with version X" warning; the pin and how it was found go in the env
  file comment.
- Promotech fallback (§17): conda `python=3.6` + pip wheels; if that fails, the Dockerfile path.

## Acceptance criteria

- [ ] `tests/unit/test_tool_envs.py`: every `workflow/envs/*.yml` passes `env_problems`; every
  `TOOLS` entry has an env file and vice versa; examples 1–8.
- [ ] Every git-pinned tool's `version` is `git:<40-hex>`; every container value matches
  `^docker://[^@\s]+@sha256:[0-9a-f]{64}$`.
- [ ] `tests/integration/test_envs_built.py` (`needs_conda`): for each env present under
  `$MAMBA_ROOT_PREFIX/envs`, `micromamba run -n giap-<tool> <version_cmd>` prints the pinned version.
- [ ] `build_envs.sbatch` ran on FarmShare for all 12 envs; job ID, wall time and env sizes recorded
  in `PROGRESS.md` (budget §3.9: 6–9 GB); lock files committed.
- [ ] `NOTICE.md` lists every tool with licence; geNOMAD marked non-commercial (R5).
- [ ] `check_imports.py` still passes; no tool source copied into the repository (R5, R11).

## Questions for the lead

1. Which commits pin UniOP, Promotech, TerminatorNet, BacTermFinder? **Proposed default:** HEAD of
   each default branch on the build date (`git ls-remote`), recorded in an ADR with the date.
2. Lock files: the spec asks the manifest for a "conda lock hash". **Proposed default:** commit
   `workflow/envs/locks/<tool>.explicit.txt` (small text) and hash it into the manifest.
3. Six tools have no upstream image (UniOP, Promotech, Promoter Calculator, TerminatorNet,
   BacTermFinder, OSTIR) until `containers.yml` pushes GHCR images. **Proposed default:** those
   processes use their conda env even under `-profile apptainer` (per-process `conda` fallback),
   and `containers.yml` is added in M2-21; nightly runs `farmshare,apptainer` for the rest.
4. `check_imports.py` does not list `terminatornet`/`bactermfinder` (unlicensed). **Proposed
   default:** add both names (one-line change, not a protected file).

## Lead decisions after code review
1. Lock files also capture pip packages: `micromamba run -n <env> pip freeze` → `workflow/envs/locks/<tool>.pip.txt`;
   the manifest's lock hash covers both files.
2. A successful build writes `<env>/.giap-built` holding the env yml's sha256; a missing or different marker
   triggers a clean rebuild (half-built envs and changed ymls are repaired).
3. The version command runs as its own step; failure or a mismatch with the registry version (or commit SHA
   for git tools) fails the build for that tool, and no lock file is written.
4. `all` builds conda envs only for tools without a pinned container, and skips opt-in tools (bactermfinder);
   `build_envs.sh <tool>` still builds any env explicitly (conda profile on demand). §3.9 disk budget.
5. The Apptainer cache file name matches Nextflow's own naming for the pinned Nextflow version (verify in its
   source/docs); if uncertain, let Nextflow pull and skip the manual pull.
6. `tests/integration/test_envs_built.py` (marker `needs_conda`) checks each built env's version command and
   lock files.
7. The `farmshare_env.sh` comment says worktree jobs share the main `local/` only when `GIAP_ROOT` is exported.
