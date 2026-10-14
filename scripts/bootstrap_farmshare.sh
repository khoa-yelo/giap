#!/bin/bash
# GIAP one-time bootstrap on FarmShare (spec §3.4). User space only; idempotent.
set -euo pipefail
source "$(dirname "$0")/farmshare_env.sh"
mkdir -p "$GIAP_ROOT"/{bin,db,data,tools,logs,cache,containers,work,worktrees,ci} "$MAMBA_ROOT_PREFIX"
grep -qx '/local/' "$GIAP_REPO/.gitignore" 2>/dev/null || echo '/local/' >> "$GIAP_REPO/.gitignore"
# uv: the official standalone installer, into local/bin
command -v uv >/dev/null || curl -LsSf https://astral.sh/uv/install.sh \
  | env UV_NO_MODIFY_PATH=1 UV_INSTALL_DIR="$GIAP_ROOT/bin" sh
# Nextflow, nf-test, Java and gh from conda-forge and bioconda only, never the 'defaults' channel
micromamba env list | grep -q 'giap-wf' || \
  micromamba create -y -n giap-wf -c conda-forge -c bioconda --override-channels \
    "nextflow=${NEXTFLOW_VERSION:-26.10.0}" nf-test=0.9.5 openjdk=21 gh
micromamba clean -a -y >/dev/null
nextflow -version && nf-test version && gh --version && uv --version
