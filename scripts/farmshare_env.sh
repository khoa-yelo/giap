# GIAP environment on FarmShare. Source it; never execute it.
unset PYTHONHOME PYTHONSTARTUP     # a user-profile PYTHONHOME can stop the venv Python from starting
export GIAP_REPO="${GIAP_REPO:-/home/classes/bios/270/khoa/giap}"
export GIAP_ROOT="$GIAP_REPO/local"                   # gitignored; all data and environments live here
export GIAP_DB_DIR="$GIAP_ROOT/db"
export GIAP_DATA_DIR="$GIAP_ROOT/data"
export GIAP_TOOLS_DIR="$GIAP_ROOT/tools"            # licence-gated binaries (DIA-NN), installed by Khoa
export GIAP_LOG_DIR="$GIAP_ROOT/logs"
export MAMBA_ROOT_PREFIX="$GIAP_ROOT/micromamba"    # GIAP's own root, separate from other projects
export CONDA_PKGS_DIRS="$MAMBA_ROOT_PREFIX/pkgs"
export UV_CACHE_DIR="$GIAP_ROOT/cache/uv"
export UV_PYTHON_INSTALL_DIR="$GIAP_ROOT/uv-python"
export UV_PROJECT_ENVIRONMENT="$GIAP_ROOT/venv"
export UV_LINK_MODE=copy
export PIP_CACHE_DIR="$GIAP_ROOT/cache/pip"
export APPTAINER_CACHEDIR="$GIAP_ROOT/cache/apptainer"
export NXF_HOME="$GIAP_ROOT/nextflow"
export NXF_APPTAINER_CACHEDIR="$GIAP_ROOT/containers"
export NXF_CONDA_CACHEDIR="$MAMBA_ROOT_PREFIX/nxf-envs"
export NXF_WORK="$GIAP_ROOT/work"
export NXF_OPTS="-Xms512m -Xmx4g"
export NUMBA_CACHE_DIR="/tmp/giap-numba-${USER}"     # node-local; NFS caches are shared across nodes
export PATH="$GIAP_ROOT/bin:$MAMBA_ROOT_PREFIX/envs/giap-wf/bin:$HOME/.local/bin:$PATH"
module load apptainer/1.4.4 >/dev/null 2>&1 || true
