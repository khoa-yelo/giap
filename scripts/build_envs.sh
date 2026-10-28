#!/bin/bash
# Build the pinned micromamba env of each named tool (GIPPA_SPEC.md §3.2, §3.9, §7.1.1; ADR-0008).
#   bash scripts/build_envs.sh <tool> [...]
#   bash scripts/build_envs.sh all   (conda envs only for tools without a pinned container in
#                                     workflow/conf/containers.config, opt-in tools skipped; §3.9)
# A named tool is always built (conda profile on demand).
# Pilot-only envs (gippa.annotate.tools.PILOT_TOOLS, M10-08) build by name with conda-forge only.
# Per tool (M2-01 review decisions 1-3):
#   - $prefix/.gippa-built holds the sha256 of workflow/envs/<tool>.yml after a successful build;
#     a missing or different marker rebuilds the env (half-built envs, changed ymls).
#   - a rebuild never deletes the old env first (backlog, job 1888853): the old prefix is renamed
#     to $prefix.gippa-old (same directory, so the same filesystem), the new env is created at
#     the real prefix and verified, and only then is the old one removed.  If the create, the
#     git checkout, the version check or the lock export fails, the new prefix is removed and the
#     old one renamed back, unchanged.  Conda envs are not relocatable (the prefix is baked into
#     shebangs, entry points and .pc/activate files), so an env is never used at any path other
#     than the one it was built at: the new env is never built at a temp path and moved in, and
#     the old env only ever runs at its own path again.  A run killed mid-build (time limit)
#     restores the old env on SIGTERM, or failing that at the start of the next run.
#   - CUDA builds on CPU nodes: if the yml asks for a CUDA build string (e.g. `pytorch=2.0.*=cuda*`)
#     and `nvidia-smi -L` lists no GPU, the solve runs with CONDA_OVERRIDE_CUDA (default 11.8; an
#     exported CONDA_OVERRIDE_CUDA wins), so the __cuda virtual package exists (job 1901404).
#   - git-pinned tools: source cloned at the pinned commit into $prefix/share/<tool>/.
#   - the version command must succeed and print the registry version (commit SHA for git tools);
#     otherwise that tool fails and no lock is written.
#   - locks: workflow/envs/locks/<tool>.explicit.txt (conda) and <tool>.pip.txt (pip freeze),
#     written to .tmp files and renamed over the old locks only after the build verified; a
#     failed build leaves the old locks as they were (they still describe the restored env).
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

# A dependency with a CUDA build string: `name=ver=*cuda*`, `name ver *cuda*` or `name[build=*cuda*]`.
wants_cuda() {
  grep -Eq '^[[:space:]]*-[[:space:]]*[A-Za-z0-9_.-]+(=[^=[:space:]#]*=[^[:space:]#]*cuda|[[:space:]]+[^[:space:]#]+[[:space:]]+[^[:space:]#]*cuda|[^#]*\[[^]#]*build[^]#]*cuda)' \
    "$ENV_DIR/$1.yml"
}
has_gpu() { nvidia-smi -L 2> /dev/null | grep -q '^GPU '; }
# CONDA_OVERRIDE_CUDA for this tool's create, or empty (never export it empty: an empty value
# hides a real GPU from the solver).
cuda_override() {
  if wants_cuda "$1" && ! has_gpu; then echo "${CONDA_OVERRIDE_CUDA:-11.8}"; fi
}
mm_create() {
  local cuda="$1"; shift
  if [ -n "$cuda" ]; then
    echo "no GPU on this node: CONDA_OVERRIDE_CUDA=$cuda for the CUDA build"
    CONDA_OVERRIDE_CUDA="$cuda" micromamba create "$@"
  else
    micromamba create "$@"
  fi
}

# The env being rebuilt (its old copy at "$CUR_PREFIX.gippa-old"), for restore_old.
CUR_PREFIX=
# Failed or interrupted rebuild: drop the new prefix, rename the old env back (unchanged).
restore_old() {
  local prefix="$1"
  [ -n "$prefix" ] && [ -d "$prefix.gippa-old" ] || return 0
  rm -rf "$prefix"
  mv "$prefix.gippa-old" "$prefix"
  echo "restored the previous ${prefix##*/}" >&2
}
# Leftover of a run killed mid-rebuild: keep the new env only if it finished (marker written).
recover_old() {
  local prefix="$1" want="$2"
  [ -d "$prefix.gippa-old" ] || return 0
  if [ -d "$prefix" ] && [ "$(cat "$prefix/.gippa-built" 2>/dev/null || true)" = "$want" ]; then
    rm -rf "$prefix.gippa-old"
  else
    restore_old "$prefix"
  fi
}
trap 'restore_old "$CUR_PREFIX"; exit 143' TERM INT

build_one() {
  local tool="$1" reg env_name version repo version_cmd prefix marker want sha src out mode chans
  local cuda
  reg="$(registry "$tool")" || return 1
  { IFS= read -r env_name; IFS= read -r version; IFS= read -r repo
    IFS= read -r version_cmd; IFS= read -r want; IFS= read -r chans; } <<< "$reg"
  mode=solve; if frozen "$tool"; then mode=frozen; fi
  prefix="$MAMBA_ROOT_PREFIX/envs/$env_name"
  marker="$prefix/.gippa-built"
  if [ "$mode" = "frozen" ]; then
    [ -s "$LOCKS/$tool.explicit.txt" ] || { echo "$tool: no committed lock" >&2; return 1; }
  fi
  recover_old "$prefix" "$want"
  CUR_PREFIX="$prefix"

  if [ -d "$prefix" ] && [ "$(cat "$marker" 2>/dev/null || true)" = "$want" ]; then
    echo "$env_name exists"
  else
    if [ -e "$prefix" ]; then
      echo "rebuilding $env_name (build marker missing or stale); the old env waits at" \
        "$prefix.gippa-old until the new one verifies"
      mv "$prefix" "$prefix.gippa-old" || return 1
    fi
    cuda="$(cuda_override "$tool")"
    if [ "$mode" = "frozen" ]; then
      echo "creating $env_name from workflow/envs/locks/$tool.explicit.txt"
      mm_create "$cuda" -y -n "$env_name" --override-channels \
        -f "$LOCKS/$tool.explicit.txt" || return 1
    else
      echo "creating $env_name from workflow/envs/$tool.yml"
      # shellcheck disable=SC2086  # $chans: "-c conda-forge [-c bioconda]"
      mm_create "$cuda" -y -n "$env_name" $chans --override-channels \
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
    rm -rf "$prefix.gippa-old"
    return 0
  fi

  # Locks (decision 1), then the marker (decision 2).
  micromamba env export -n "$env_name" --explicit > "$LOCKS/$tool.explicit.txt.tmp" || return 1
  PYTHONNOUSERSITE=1 micromamba run -n "$env_name" pip freeze > "$LOCKS/$tool.pip.txt.tmp" \
    || return 1
  mv "$LOCKS/$tool.explicit.txt.tmp" "$LOCKS/$tool.explicit.txt"
  mv "$LOCKS/$tool.pip.txt.tmp" "$LOCKS/$tool.pip.txt"
  echo "$want" > "$marker"
  rm -rf "$prefix.gippa-old"
}

failed=
for tool in "$@"; do
  CUR_PREFIX=
  if ! build_one "$tool"; then
    rm -f "$LOCKS/$tool".*.tmp
    restore_old "$CUR_PREFIX"  # the old env and its locks stay as they were
    failed="$failed $tool"
    echo "FAILED: $tool" >&2
  fi
  CUR_PREFIX=
done
micromamba clean -a -y > /dev/null
if [ -n "$failed" ]; then
  echo "build failed for:$failed" >&2
  exit 1
fi
