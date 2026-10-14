#!/bin/bash
# Build the pinned micromamba env of each named tool (GIAP_SPEC.md §3.2, §3.9, §7.1.1; ADR-0008).
#   bash scripts/build_envs.sh <tool> [...]
#   bash scripts/build_envs.sh all   (conda envs only for tools without a pinned container in
#                                     workflow/conf/containers.config, opt-in tools skipped; §3.9)
# A named tool is always built (conda profile on demand).
# Per tool (M2-01 review decisions 1-3):
#   - $prefix/.giap-built holds the sha256 of workflow/envs/<tool>.yml after a successful build;
#     a missing or different marker removes the env and rebuilds it (half-built envs, changed ymls).
#   - git-pinned tools: source cloned at the pinned commit into $prefix/share/<tool>/.
#   - the version command must succeed and print the registry version (commit SHA for git tools);
#     otherwise that tool fails and no lock is written.
#   - locks: workflow/envs/locks/<tool>.explicit.txt (conda) and <tool>.pip.txt (pip freeze).
# The script continues with the next tool after a failure and exits 1 if any tool failed.
# Run it through `sbatch scripts/slurm/build_envs.sbatch`, never on a login node.
set -euo pipefail

GIAP_REPO="${GIAP_REPO:-/home/classes/bios/270/khoa/giap}"
ENV_DIR="$GIAP_REPO/workflow/envs"
mapfile -t VALID < <(find "$ENV_DIR" -maxdepth 1 -name '*.yml' -printf '%f\n' | sed 's/\.yml$//' | sort)

usage() { echo "usage: build_envs.sh <tool>... | all   (valid: ${VALID[*]})" >&2; }
[ "$#" -gt 0 ] || { usage; exit 2; }
for tool in "$@"; do
  if [ "$tool" != "all" ] && [[ ! " ${VALID[*]} " == *" $tool "* ]]; then
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

# `all`: conda tools that have no pinned container and are not opt-in (decision 4).
if [ "$1" = "all" ]; then
  mapfile -t ALL < <(PYTHONPATH="$GIAP_REPO/src" "$PY" - "$GIAP_REPO" <<'PYEOF'
import sys
from pathlib import Path
from giap.annotate.tools import TOOLS, read_containers_config
pinned = {k.lower() for k in read_containers_config(Path(sys.argv[1]) / "workflow/conf/containers.config")}
for name, t in sorted(TOOLS.items()):
    if not t.in_process and t.module != "opt_in" and name not in pinned:
        print(name)
PYEOF
)
  set -- "${ALL[@]}"
  echo "all: ${*}"
fi

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

build_one() {
  local tool="$1" reg env_name version repo version_cmd prefix marker want sha src out
  mapfile -t reg < <(registry "$tool")
  env_name="${reg[0]}" version="${reg[1]}" repo="${reg[2]}" version_cmd="${reg[3]}"
  prefix="$MAMBA_ROOT_PREFIX/envs/$env_name"
  marker="$prefix/.giap-built"
  want="$(sha256sum "$ENV_DIR/$tool.yml" | awk '{print $1}')"
  rm -f "$LOCKS/$tool.explicit.txt" "$LOCKS/$tool.pip.txt"

  if [ -d "$prefix" ] && [ "$(cat "$marker" 2>/dev/null || true)" = "$want" ]; then
    echo "$env_name exists"
  else
    if [ -e "$prefix" ]; then
      echo "rebuilding $env_name (build marker missing or stale)"
      micromamba env remove -y -n "$env_name" > /dev/null || true
      rm -rf "$prefix"
    fi
    echo "creating $env_name from workflow/envs/$tool.yml"
    micromamba create -y -n "$env_name" -c conda-forge -c bioconda --override-channels \
      -f "$ENV_DIR/$tool.yml" || return 1
  fi
  if [ -n "$repo" ]; then
    sha="${version#git:}"
    src="$prefix/share/$tool"
    if [ "$(cat "$src/GIAP_COMMIT" 2>/dev/null || true)" != "$sha" ]; then
      rm -rf "$src"
      git clone --quiet "$repo" "$src" || return 1
      git -C "$src" checkout --quiet "$sha" || return 1
      echo "$sha" > "$src/GIAP_COMMIT"
    fi
  fi

  # Version step (decision 3).
  if ! out="$(eval "micromamba run -n $env_name $version_cmd" 2>&1)"; then
    echo "$env_name: version command failed: $(head -n 3 <<< "$out")" >&2
    return 1
  fi
  echo "$env_name version: $(head -n 3 <<< "$out")"
  if [[ "$out" != *"${version#git:}"* ]]; then
    echo "$env_name: version mismatch, expected ${version#git:}" >&2
    return 1
  fi

  # Locks (decision 1), then the marker (decision 2).
  micromamba env export -n "$env_name" --explicit > "$LOCKS/$tool.explicit.txt.tmp" || return 1
  micromamba run -n "$env_name" pip freeze > "$LOCKS/$tool.pip.txt.tmp" || return 1
  mv "$LOCKS/$tool.explicit.txt.tmp" "$LOCKS/$tool.explicit.txt"
  mv "$LOCKS/$tool.pip.txt.tmp" "$LOCKS/$tool.pip.txt"
  echo "$want" > "$marker"
}

failed=()
for tool in "$@"; do
  if ! build_one "$tool"; then
    rm -f "$LOCKS/$tool".*.tmp "$LOCKS/$tool.explicit.txt" "$LOCKS/$tool.pip.txt"
    failed+=("$tool")
    echo "FAILED: $tool" >&2
  fi
done
micromamba clean -a -y > /dev/null
if [ "${#failed[@]}" -gt 0 ]; then
  echo "build failed for: ${failed[*]}" >&2
  exit 1
fi
