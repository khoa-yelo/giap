# ADR-0008: Annotation tool pins, image digests and licences (M2-01)
- Status: accepted
- Date: 2026-10-14
- Spec sections: §3.1, §3.2, §4, §6.3, §7.1.1, §17, R5; M2-index lead decisions 1 and 8
- Needs Khoa: no (lead decision 1 accepted pinning source-only tools to the default-branch commit)

## Context
§7.1.1 fixes release versions for most tools but says only "pinned commit" for the four
source-only tools, and "pin by digest" for images. M2-01 needs exact values for the registry
(`src/gippa/annotate/tools.py`), the env files (`workflow/envs/*.yml`), `workflow/conf/containers.config`
and the Dockerfiles. All values below were resolved on 2026-10-14 with read-only queries against the
§3.2/§18 sources (`git ls-remote` and the GitHub API, `micromamba search` on conda-forge + bioconda,
the PyPI JSON API, and registry v2 manifests via `scripts/resolve_digest.py`).

## Decision

**Source-only tools: default-branch HEAD on 2026-10-14.** Licences read from each repository's
LICENSE at the pinned commit.

| Tool | Repository | Branch | Commit | Commit date | Licence |
|---|---|---|---|---|---|
| UniOP | https://github.com/hongsua/UniOP | main | `3ed14dfc2a5d69b993d99a26e50c1f8e9723c087` | 2026-09-02 | none (no LICENSE file in the tree) → `unlicensed` |
| Promotech | https://github.com/BioinformaticsLabAtMUN/PromoTech | master | `56251ad9b883ef831b4753fc623d5ec970fe65e0` | 2023-06-19 | GPL-3.0 |
| TerminatorNet | https://github.com/btjaden/TerminatorNet | main | `72180309545988dba63b8bd5db5d0b09bf4dba70` | 2025-06-29 | MIT ("Copyright (c) 2025 btjaden") |
| BacTermFinder | https://github.com/BioinformaticsLabAtMUN/BacTermFinder (cited by PMC11890068) | main | `024870447f6288ebc9d29b314e957a6dfdfeb5f6` | 2026-05-06 | GPL-3.0 |

TerminatorNet and BacTermFinder stay in `scripts/ci/check_imports.py` `FORBIDDEN`: they run only as
separate processes (design rule 4), whatever their licence.

**Released tools.**

| Tool | Version | Conda / PyPI pin | Licence (source) |
|---|---|---|---|
| Bakta | 1.12.2 | `bakta=1.12.2` (bioconda, `pyhdfd78af_0`) | GPL-3.0 (GitHub) |
| Promoter Calculator | 1.2.4 | pip `promotercalculator==1.2.4` | GPL-3.0-or-later (wheel classifier) |
| TransTermHP | 2.09 | `transtermhp=2.09` (bioconda, `h9948957_2`) | GPL-3.0-or-later (bioconda recipe) |
| OSTIR | 1.1.3 | pip `ostir==1.1.3` + `viennarna=2.6.4` | GPL-3.0 (GitHub) |
| codon-bias | 0.5.0 | main GIPPA env (in-process) | MIT (PyPI, GitHub) |
| geNOMAD | 1.12.0 | `genomad=1.12.0` (bioconda, `pyhdfd78af_0`) | non-commercial (§7.1.1; GitHub reports NOASSERTION) |
| antiSMASH | 8.0.4 | `antismash=8.0.4` (bioconda) | AGPL-3.0 (GitHub) |
| miniprot | 0.18 | `miniprot=0.18` (bioconda, `h577a1d6_0`) | MIT (GitHub) |
| BLAST+ | 2.16.0 | `blast=2.16.0` (bioconda) | public domain (NCBI-PD, bioconda) |

BLAST+ 2.16.0 rather than the newer 2.17.0 so the env matches FarmShare's `blast-plus/2.16.0` module.

**Env-file choices** (full lists in each `workflow/envs/<tool>.yml`):
- UniOP: Python 3.11.5 (README), pandas 2.1.1, numpy 1.26.0, scikit-learn 1.3.0.
- Promotech: Python 3.6.10, scikit-learn 0.23.1 and the other packages from conda-forge/bioconda at
  upstream `promotech_env.yml` versions; TensorFlow 2.1.0 is not on conda-forge, so it is a pip
  cp36 wheel with `protobuf==3.19.6` and `progressbar2==3.37.1` (§17 fallback "conda python=3.6 +
  pip wheels"; the Dockerfile is the next fallback).
- TerminatorNet: scikit-learn **1.1.2**, the version that pickled `model.pkl` (G12). Read by
  disassembling `model.pkl` at the pinned commit with `pickletools.dis` (no unpickling): the
  `MLPClassifier` state carries `_sklearn_version = '1.1.2'`. Python 3.10.14, numpy 1.23.5, scipy 1.9.3.
- BacTermFinder: Python 3.9.19 (README), the `requirements.txt` versions of the packages its scripts
  import (numpy 1.23.0, pandas 1.5.0, biopython 1.81, h5py 3.7.0, joblib 1.2.0, tqdm 4.64.1; pip
  tensorflow 2.8.0, keras 2.8.0, protobuf 3.20.0).
- Promoter Calculator: Python 3.10.14 with biopython 1.84, numpy 1.26.4, scipy 1.14.1,
  scikit-learn 1.5.2 (the wheel declares them without versions and ships no pickled model).
- OSTIR: Python 3.10.14, numpy 1.26.4, viennarna 2.6.4.

**Image digests** (`workflow/conf/containers.config`; `docker://…@sha256:…`, multi-arch index digest):

| Process | Tag resolved | Digest |
|---|---|---|
| BAKTA | `oschwengers/bakta:v1.12.2` | `sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18` |
| TRANSTERMHP | `quay.io/biocontainers/transtermhp:2.09--h9948957_2` | `sha256:fc6fe0404be9340a883ec7a7b90415dfc3826c2bca7333a0314335c7bfee2d03` |
| GENOMAD | `quay.io/biocontainers/genomad:1.12.0--pyhdfd78af_0` | `sha256:d54652ff868af72604ea0e54c541c5c92168c7ac821a422de05dcbad65509e5c` |
| ANTISMASH | `antismash/standalone:8.0.4` | `sha256:d3f390e687cfe02f1387be40256dddbf678e91bc56b96c36da07b56ec8bda29d` |
| MINIPROT | `quay.io/biocontainers/miniprot:0.18--h577a1d6_0` | `sha256:2eb53fea53743b56b1c10ac4b202469b92db959ffbb8bafd83156fdbde80e22f` |
| BLAST | `quay.io/biocontainers/blast:2.16.0--h66d330f_5` | `sha256:9dfc69f990c0aeb936276ee591ed32919a79f46dfa34060055e05a050a17959c` |
| (Dockerfile base) | `mambaorg/micromamba:2.3.3-ubuntu24.04` | `sha256:f90ebbd10a6dabe1ab87497babd33085d6cf743da6176cf6bf119b9431cdf251` |

UniOP, Promotech, Promoter Calculator, TerminatorNet, BacTermFinder and OSTIR have no upstream image:
they keep their conda env under `-profile apptainer` until `containers.yml` (M2-21) publishes
`ghcr.io/<owner>/gippa-<tool>` from `workflow/containers/<tool>/Dockerfile`. The Dockerfiles fetch
the pinned commit's GitHub archive at build time (never vendored, R5). The UniOP image has no
licence to redistribute and must not be pushed to a public registry without the authors' consent.

## Evidence (before → after: metrics, tests, simulations)
`tests/unit/annotate/` (M2-01) checks every pin, digest format, licence and env file. The envs are
not built yet: `sbatch scripts/slurm/build_envs.sbatch all` will confirm that each env solves, and
will write `workflow/envs/locks/<tool>.explicit.txt`.

## Consequences
Upgrading any pin, digest or commit needs a new ADR, a re-run of `build_envs.sbatch` and the
contract-test diffs (§14.6). If an env does not solve (most likely Promotech or BacTermFinder), the
next fallback is the Dockerfile path (§17), recorded in a follow-up ADR.
