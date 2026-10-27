#!/usr/bin/env bash
# Dev build of the step images (M11-15 item 8a, slice 15g; the Q-K39 ADR, ADR-0065).
#   scripts/dev/build_step_image.sh [--runtime docker|apptainer|singularity]
# Builds the wheel and both images from this checkout (a temporary context of the locks and the
# wheel only), tags <pkg>-step:dev-<sha12 of the wheel>
# and <pkg>-step-r:dev-<sha12>, and prints the refs for --step-image / --step-image-r.
# docker builds the images; apptainer/singularity build with docker first (a Dockerfile needs a
# Docker-compatible builder) and convert to SIF files under dist/ (docker-daemon://).
set -euo pipefail

usage() { echo "usage: $0 [--runtime docker|apptainer|singularity]" >&2; exit 2; }

runtime=docker
while [ $# -gt 0 ]; do
  case "$1" in
    --runtime) [ $# -ge 2 ] || usage; runtime="$2"; shift 2 ;;
    --runtime=*) runtime="${1#--runtime=}"; shift ;;
    -h|--help) echo "usage: $0 [--runtime docker|apptainer|singularity]"; exit 0 ;;
    *) echo "build_step_image: unknown argument: $1" >&2; usage ;;
  esac
done
case "$runtime" in
  docker|apptainer|singularity) ;;
  *) echo "build_step_image: unknown runtime '$runtime' (docker, apptainer or singularity)" >&2
     exit 2 ;;
esac

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$root"
pkg="$(sed -n 's/^name = "\(.*\)"$/\1/p' pyproject.toml | head -n 1)"
for bin in uv docker "$runtime"; do
  command -v "$bin" >/dev/null || { echo "build_step_image: $bin not found" >&2; exit 1; }
done

out="$(mktemp -d "${TMPDIR:-/tmp}/step-wheel.XXXXXX")"
trap 'rm -rf "$out"' EXIT
uv build --wheel --out-dir "$out/wheel" >&2
wheel_path="$(ls "$out"/wheel/*.whl)"
wheel="$(basename "$wheel_path")"
mkdir -p dist
cp "$wheel_path" "dist/$wheel"
sha="$(sha256sum "dist/$wheel" | cut -d' ' -f1)"
sha12="$(echo "$sha" | cut -c1-12)"
version="$(echo "$wheel" | cut -d- -f2)"
revision="$(git rev-parse HEAD 2>/dev/null || echo unknown)"
# A minimal build context (never the checkout: local/ holds the venv and the worktrees), with the
# repository-relative paths the Dockerfiles bind-mount.
ctx="$out/context"
for f in workflow/containers/step/requirements.lock workflow/containers/step/build-requirements.lock \
         workflow/envs/locks/fgsea.explicit.txt "dist/$wheel"; do
  mkdir -p "$ctx/$(dirname "$f")"
  cp "$f" "$ctx/$f"
done

build() {  # build <Dockerfile> <ref>; prints the ref to pass to gippa run
  local dockerfile="$1" ref="$2" sif
  docker build -f "$dockerfile" \
    --build-arg WHEEL="$wheel" --build-arg WHEEL_SHA256="$sha" \
    --build-arg VERSION="$version" --build-arg REVISION="$revision" -t "$ref" "$ctx" >&2
  if [ "$runtime" = docker ]; then
    echo "$ref"
  else
    sif="$root/dist/$(echo "$ref" | tr ':' '-').sif"
    "$runtime" build --force "$sif" "docker-daemon://$ref" >&2
    echo "$sif"
  fi
}
step="$(build workflow/containers/step/Dockerfile "$pkg-step:dev-$sha12")"
step_r="$(build workflow/containers/step-r/Dockerfile "$pkg-step-r:dev-$sha12")"
echo "--step-image $step"
echo "--step-image-r $step_r"
