#!/bin/bash
# Build the pinned micromamba env of each named tool (GIAP_SPEC.md §3.2, §3.9, §7.1.1; ADR-0008).
#   bash scripts/build_envs.sh <tool> [...]      (all: every workflow/envs/*.yml)
# Idempotent: an existing env giap-<tool> is left alone ("giap-<tool> exists"). For git-pinned
# tools the source is cloned at the pinned commit into $MAMBA_ROOT_PREFIX/envs/giap-<tool>/share/<tool>/.
# Then the tool's version command runs, the explicit lock is written to
# workflow/envs/locks/<tool>.explicit.txt, and the package cache is cleaned.
# Run it through `sbatch scripts/slurm/build_envs.sbatch`, never on a login node.
set -euo pipefail

GIAP_REPO="${GIAP_REPO:-/home/classes/bios/270/khoa/giap}"
ENV_DIR="$GIAP_REPO/workflow/envs"
mapfile -t VALID < <(find "$ENV_DIR" -maxdepth 1 -name '*.yml' -printf '%f\n' | sed 's/\.yml$//' | sort)

usage() { echo "usage: build_envs.sh <tool>... | all   (valid: ${VALID[*]})" >&2; }
[ "$#" -gt 0 ] || { usage; exit 2; }
if [ "$1" = "all" ]; then set -- "${VALID[@]}"; fi
for tool in "$@"; do
  if [[ ! " ${VALID[*]} " == *" $tool "* ]]; then
    echo "unknown tool: $tool" >&2
    usage
    exit 2
  fi
done

# shellcheck source=/dev/null
source "$GIAP_REPO/scripts/farmshare_env.sh"
PY="${GIAP_PYTHON:-$UV_PROJECT_ENVIRONMENT/bin/python}"
LOCKS="$GIAP_REPO/workflow/envs/locks"
mkdir -p "$LOCKS"

# Registry fields (single source of truth: src/giap/annotate/tools.py), one per line.
registry() {
  PYTHONPATH="$GIAP_REPO/src" "$PY" - "$1" <<'PYEOF'
import shlex, sys
from giap.annotate.tools import TOOLS
t = TOOLS[sys.argv[1]]
print(t.env_name)
print(t.version)
print(t.source_repo or "")
print(shlex.join(t.version_cmd))
PYEOF
}

for tool in "$@"; do
  mapfile -t reg < <(registry "$tool")
  env_name="${reg[0]}" version="${reg[1]}" repo="${reg[2]}" version_cmd="${reg[3]}"
  prefix="$MAMBA_ROOT_PREFIX/envs/$env_name"
  if micromamba env list | awk '{print $1}' | grep -qx "$env_name"; then
    echo "$env_name exists"
  else
    echo "creating $env_name from workflow/envs/$tool.yml"
    micromamba create -y -n "$env_name" -c conda-forge -c bioconda --override-channels \
      -f "$ENV_DIR/$tool.yml"
  fi
  if [ -n "$repo" ]; then
    sha="${version#git:}"
    src="$prefix/share/$tool"
    if [ "$(cat "$src/GIAP_COMMIT" 2>/dev/null || true)" != "$sha" ]; then
      rm -rf "$src"
      git clone --quiet "$repo" "$src"
      git -C "$src" checkout --quiet "$sha"
      echo "$sha" > "$src/GIAP_COMMIT"
    fi
  fi
  echo "$env_name version: $(eval "micromamba run -n $env_name $version_cmd" 2>&1 | head -n 3)"
  micromamba env export -n "$env_name" --explicit > "$LOCKS/$tool.explicit.txt"
done
micromamba clean -a -y > /dev/null
