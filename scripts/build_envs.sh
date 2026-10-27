#!/bin/bash
# Build the pinned micromamba env of each named tool (GIPPA_SPEC.md §3.2, §3.9, §7.1.1; ADR-0008).
#   bash scripts/build_envs.sh <tool> [...]
#   bash scripts/build_envs.sh all   (conda envs only for tools without a pinned container in
#                                     workflow/conf/containers.config, opt-in tools skipped; §3.9)
# A named tool is always built (conda profile on demand).
# Pilot-only envs (gippa.annotate.tools.PILOT_TOOLS, M10-08) build by name with conda-forge only.
# Per tool (M2-01 review decisions 1-3):
#   - $prefix/.gippa-built holds the sha256 of workflow/envs/<tool>.yml after a successful build;
#     a missing or different marker removes the env and rebuilds it (half-built envs, changed ymls).
#   - git-pinned tools: source cloned at the pinned commit into $prefix/share/<tool>/.
#   - the version command must succeed and print the registry version (commit SHA for git tools);
#     otherwise that tool fails and no lock is written.
#   - locks: workflow/envs/locks/<tool>.explicit.txt (conda) and <tool>.pip.txt (pip freeze).
#   - frozen tools (a `# frozen-lock` line in the yml; M8-17a: fgsea) are installed from their
#     committed explicit lock instead of a solve; the export must list exactly the lock's
#     packages, the lock is never rewritten or removed, and there is no pip lock.
# The script continues with the next tool after a failure and exits 1 if any tool failed.
# Run it through `sbatch scripts/slurm/build_envs.sbatch`, never on a login node.
set -euo pipefail

GIPPA_REPO="${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}"
ENV_DIR="$GIPPA_REPO/workflow/envs"
# Portable to bash 3.2 (macOS CI): no mapfile, associative arrays or GNU-only find/sha256sum.
VALID=""
for f in "$ENV_DIR"/*.yml; do
  [ -e "$f" ] || continue
  b="${f##*/}"
  VALID="$VALID ${b%.yml}"
done
VALID="${VALID# }"

usage() { echo "usage: build_envs.sh <tool>... | all   (valid: $VALID)" >&2; }
[ "$#" -gt 0 ] || { usage; exit 2; }
for tool in "$@"; do
  if [ "$tool" != "all" ] && [[ ! " $VALID " == *" $tool "* ]]; then
    echo "unknown tool: $tool" >&2
    usage
    exit 2
  fi
done

# shellcheck source=/dev/null
source "$GIPPA_REPO/scripts/farmshare_env.sh"
# The env's Python must never see ~/.local: a user-site package would leak into the version step
# and into `pip freeze`, so the pip lock would pin it (M11-14b; e.g. ViennaRNA 2.7.1 over conda's).
export PYTHONNOUSERSITE=1
PY="${GIPPA_PYTHON:-$UV_PROJECT_ENVIRONMENT/bin/python}"
LOCKS="$GIPPA_REPO/workflow/envs/locks"
mkdir -p "$LOCKS"

# `all`: conda tools that have no pinned container and are not opt-in (decision 4).
if [ "$1" = "all" ]; then
  all_tools() {
    PYTHONPATH="$GIPPA_REPO/src" "$PY" - "$GIPPA_REPO" <<'PYEOF'
import sys
from pathlib import Path
from gippa.annotate.tools import TOOLS, read_containers_config
pinned = {k.lower() for k in read_containers_config(Path(sys.argv[1]) / "workflow/conf/containers.config")}
for name, t in sorted(TOOLS.items()):
    if not t.in_process and t.module != "opt_in" and name not in pinned:
        print(name)
PYEOF
  }
  ALL="$(all_tools)"
  # shellcheck disable=SC2086  # tool names are plain words
  set -- $ALL
  echo "all: $*"
fi

# Registry fields (single source of truth: src/gippa/annotate/tools.py), one per line.
registry() {
  PYTHONPATH="$GIPPA_REPO/src" "$PY" - "$1" "$ENV_DIR" <<'PYEOF'
import hashlib, shlex, sys
from pathlib import Path
from gippa.annotate.tools import PILOT_TOOLS, TOOLS
from gippa.io.fgsea_r import FGSEA_TOOL
frozen = {FGSEA_TOOL.name: FGSEA_TOOL}  # joins TOOLS in M8-17b
t = {**frozen, **PILOT_TOOLS, **TOOLS}[sys.argv[1]]  # PILOT_TOOLS: M10-08, never in `all`
print(t.env_name)
print(t.version)
print(t.source_repo or "")
print(shlex.join(t.version_cmd))
print(hashlib.sha256((Path(sys.argv[2]) / f"{sys.argv[1]}.yml").read_bytes()).hexdigest())
print("-c conda-forge" if t.module == "pilot" else "-c conda-forge -c bioconda")  # M10-08
PYEOF
}

# The committed explicit lock is the source of the env (never rewritten or removed here).
frozen() { grep -q '^# frozen-lock' "$ENV_DIR/$1.yml"; }

build_one() {
  local tool="$1" reg env_name version repo version_cmd prefix marker want sha src out mode chans
  reg="$(registry "$tool")" || return 1
  { IFS= read -r env_name; IFS= read -r version; IFS= read -r repo
    IFS= read -r version_cmd; IFS= read -r want; IFS= read -r chans; } <<< "$reg"
  mode=solve; if frozen "$tool"; then mode=frozen; fi
  prefix="$MAMBA_ROOT_PREFIX/envs/$env_name"
  marker="$prefix/.gippa-built"
  if [ "$mode" = "frozen" ]; then
    [ -s "$LOCKS/$tool.explicit.txt" ] || { echo "$tool: no committed lock" >&2; return 1; }
  else
    rm -f "$LOCKS/$tool.explicit.txt" "$LOCKS/$tool.pip.txt"
  fi

  if [ -d "$prefix" ] && [ "$(cat "$marker" 2>/dev/null || true)" = "$want" ]; then
    echo "$env_name exists"
  else
    if [ -e "$prefix" ]; then
      echo "rebuilding $env_name (build marker missing or stale)"
      micromamba env remove -y -n "$env_name" > /dev/null || true
      rm -rf "$prefix"
    fi
    if [ "$mode" = "frozen" ]; then
      echo "creating $env_name from workflow/envs/locks/$tool.explicit.txt"
      micromamba create -y -n "$env_name" --override-channels \
        -f "$LOCKS/$tool.explicit.txt" || return 1
    else
      echo "creating $env_name from workflow/envs/$tool.yml"
      # shellcheck disable=SC2086  # $chans: "-c conda-forge [-c bioconda]"
      micromamba create -y -n "$env_name" $chans --override-channels \
        -f "$ENV_DIR/$tool.yml" || return 1
    fi
  fi
  if [ -n "$repo" ]; then
    sha="${version#git:}"
    src="$prefix/share/$tool"
    if [ "$(cat "$src/GIPPA_COMMIT" 2>/dev/null || true)" != "$sha" ]; then
      rm -rf "$src"
      git clone --quiet "$repo" "$src" || return 1
      git -C "$src" checkout --quiet "$sha" || return 1
      echo "$sha" > "$src/GIPPA_COMMIT"
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

  if [ "$mode" = "frozen" ]; then  # the env must be exactly the committed lock
    micromamba env export -n "$env_name" --explicit > "$LOCKS/$tool.explicit.txt.tmp" || return 1
    if ! diff <(grep '^https://' "$LOCKS/$tool.explicit.txt" | sort) \
              <(grep '^https://' "$LOCKS/$tool.explicit.txt.tmp" | sort); then
      echo "$env_name: installed packages differ from locks/$tool.explicit.txt" >&2
      return 1
    fi
    rm -f "$LOCKS/$tool.explicit.txt.tmp"
    echo "$want" > "$marker"
    return 0
  fi

  # Locks (decision 1), then the marker (decision 2).
  micromamba env export -n "$env_name" --explicit > "$LOCKS/$tool.explicit.txt.tmp" || return 1
  PYTHONNOUSERSITE=1 micromamba run -n "$env_name" pip freeze > "$LOCKS/$tool.pip.txt.tmp" \
    || return 1
  mv "$LOCKS/$tool.explicit.txt.tmp" "$LOCKS/$tool.explicit.txt"
  mv "$LOCKS/$tool.pip.txt.tmp" "$LOCKS/$tool.pip.txt"
  echo "$want" > "$marker"
}

failed=
for tool in "$@"; do
  if ! build_one "$tool"; then
    rm -f "$LOCKS/$tool".*.tmp
    if ! frozen "$tool"; then
      rm -f "$LOCKS/$tool.explicit.txt" "$LOCKS/$tool.pip.txt"
    fi
    failed="$failed $tool"
    echo "FAILED: $tool" >&2
  fi
done
micromamba clean -a -y > /dev/null
if [ -n "$failed" ]; then
  echo "build failed for:$failed" >&2
  exit 1
fi
