#!/bin/bash
#SBATCH --job-name=giap-fetch-bakta-db
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Builds the Bakta env and downloads the light and full DB v6.x into $GIAP_DB_DIR/bakta
# (spec §3.9, §7.9.1, G17). Long-latency prefetch for M2; `giap setup` later verifies and records it.
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
cd "$GIAP_REPO"
micromamba env list | grep -q 'giap-bakta' || \
  micromamba create -y -n giap-bakta -c conda-forge -c bioconda --override-channels -f workflow/envs/bakta.yml
micromamba clean -a -y > /dev/null
mkdir -p "$GIAP_DB_DIR/bakta"
run() { micromamba run -n giap-bakta "$@"; }
run bakta --version
if [ ! -d "$GIAP_DB_DIR/bakta/db-light" ]; then
  run bakta_db download --output "$GIAP_DB_DIR/bakta" --type light
fi
if [ ! -d "$GIAP_DB_DIR/bakta/db" ]; then
  run bakta_db download --output "$GIAP_DB_DIR/bakta" --type full
fi
ls -la "$GIAP_DB_DIR/bakta"; du -sh "$GIAP_DB_DIR/bakta"/*
cat "$GIAP_DB_DIR/bakta"/db*/version.json 2>/dev/null || true
