# GIPPA environment on FarmShare. Source it; never execute it.
unset PYTHONHOME PYTHONSTARTUP     # a user-profile PYTHONHOME can stop the venv Python from starting
export GIPPA_REPO="${GIPPA_REPO:-${GIAP_REPO:-/home/classes/bios/270/khoa/giap}}"
export GIPPA_ROOT="${GIPPA_ROOT:-${GIAP_ROOT:-$GIPPA_REPO/local}}"   # gitignored; worktree jobs share the main local/ only when GIPPA_ROOT is exported
export GIPPA_DB_DIR="$GIPPA_ROOT/db"
export GIPPA_DATA_DIR="$GIPPA_ROOT/data"
export GIPPA_TOOLS_DIR="$GIPPA_ROOT/tools"            # licence-gated binaries (DIA-NN), installed by Khoa
export GIPPA_LOG_DIR="$GIPPA_ROOT/logs"
export MAMBA_ROOT_PREFIX="$GIPPA_ROOT/micromamba"    # GIPPA's own root, separate from other projects
export CONDA_PKGS_DIRS="$MAMBA_ROOT_PREFIX/pkgs"
export UV_CACHE_DIR="$GIPPA_ROOT/cache/uv"
export UV_PYTHON_INSTALL_DIR="$GIPPA_ROOT/uv-python"
export UV_PROJECT_ENVIRONMENT="$GIPPA_ROOT/venv"
export UV_LINK_MODE=copy
export PIP_CACHE_DIR="$GIPPA_ROOT/cache/pip"
export APPTAINER_CACHEDIR="$GIPPA_ROOT/cache/apptainer"
export NXF_HOME="$GIPPA_ROOT/nextflow"
export NXF_APPTAINER_CACHEDIR="$GIPPA_ROOT/containers"
export NXF_CONDA_CACHEDIR="$MAMBA_ROOT_PREFIX/nxf-envs"
export NXF_WORK="$GIPPA_ROOT/work"
export NXF_OPTS="-Xms512m -Xmx4g"
export PYTHONNOUSERSITE=1   # py3.10 tool envs would otherwise import ~/.local packages (ViennaRNA, biopython)
export NUMBA_CACHE_DIR="/tmp/gippa-numba-${USER:-$(id -un)}"     # node-local; NFS caches are shared across nodes
export PATH="$GIPPA_ROOT/bin:$MAMBA_ROOT_PREFIX/envs/gippa-wf/bin:$MAMBA_ROOT_PREFIX/envs/giap-wf/bin:$HOME/.local/bin:$PATH"
module load apptainer/1.4.4 >/dev/null 2>&1 || true
# The interactive allocation exports SLURM_GET_USER_ENV=1 and SLURM_EXPORT_ENV=NONE; every sbatch from the
# session inherited them, so batch nodes had to rebuild the login environment at launch and ~1 in 5 jobs
# failed with "user env retrieval failed" (2026-10-15, wheat-04/barley-02). Jobs get this script's env instead.
unset SLURM_GET_USER_ENV SLURM_EXPORT_ENV
# ADR-0064 decision 4: the old names, same values, for jobs in flight, old worktrees
# and the installed scrontab; removed after the first publish.
export GIAP_REPO="$GIPPA_REPO"
export GIAP_ROOT="$GIPPA_ROOT"
export GIAP_DB_DIR="$GIPPA_DB_DIR"
export GIAP_DATA_DIR="$GIPPA_DATA_DIR"
export GIAP_TOOLS_DIR="$GIPPA_TOOLS_DIR"
export GIAP_LOG_DIR="$GIPPA_LOG_DIR"
