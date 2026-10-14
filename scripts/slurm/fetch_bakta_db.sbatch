#!/bin/bash
#SBATCH --job-name=giap-fetch-bakta-db
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Downloads the Bakta light and full DB v6.0 into $GIAP_DB_DIR/bakta (spec §3.9, §7.9.1, G17).
# Long-latency prefetch for M2; `giap setup` later verifies and records it.
# Resumable: aria2c keeps a .aria2 control file next to each tarball, so resubmitting after a
# timeout continues where it stopped. Source and checksums: the official Zenodo record (R6).
# Install uses the pinned Bakta container (workflow/conf/containers.config) — `bakta_db install`.
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
cd "$GIAP_REPO"

RECORD=14916843                                   # Bakta DB v6.0, 2025-02-24, doi:10.5281/zenodo.14916843
BASE="https://zenodo.org/records/$RECORD/files"
declare -A MD5=(
  [db-light.tar.xz]=4a6e059ded39e9c5537ef4137d2f5648
  [db.tar.xz]=4c1115e40abfa2b464ae5dd988bdd88e
)
declare -A DIR=([db-light.tar.xz]=db-light [db.tar.xz]=db)
OUT="$GIAP_DB_DIR/bakta"
IMG=$(ls "$GIAP_ROOT"/containers/oschwengers-bakta@sha256-*.img | head -1)
mkdir -p "$OUT"

for f in db-light.tar.xz db.tar.xz; do
  if [ -d "$OUT/${DIR[$f]}" ] && [ -f "$OUT/${DIR[$f]}/version.json" ]; then
    echo "$f: ${DIR[$f]} already installed"; continue
  fi
  # A partial file without an aria2 control file (e.g. from `bakta_db download`) cannot be resumed safely.
  if [ -f "$OUT/$f" ] && [ ! -f "$OUT/$f.aria2" ]; then
    if ! echo "${MD5[$f]}  $OUT/$f" | md5sum -c --status; then
      echo "$f: removing partial download without a resume record"; rm -f "$OUT/$f"
    fi
  fi
  if [ ! -f "$OUT/$f" ] || [ -f "$OUT/$f.aria2" ]; then
    echo "$f: downloading (resumable) from $BASE/$f"
    aria2c --continue=true --max-connection-per-server=8 --split=8 --min-split-size=64M \
      --max-tries=0 --retry-wait=30 --timeout=120 --file-allocation=none \
      --auto-file-renaming=false --allow-overwrite=false --summary-interval=600 \
      --console-log-level=warn --dir="$OUT" --out="$f" "$BASE/$f?download=1"
  fi
  echo "$f: verifying md5"
  echo "${MD5[$f]}  $OUT/$f" | md5sum -c -
  echo "$f: installing with $(basename "$IMG")"
  apptainer exec --bind "$OUT" "$IMG" bakta_db install --db-file "$OUT/$f" --output "$OUT"
  [ -f "$OUT/${DIR[$f]}/version.json" ] && rm -f "$OUT/$f"
done
ls -la "$OUT"; du -sh "$OUT"/*
cat "$OUT"/db*/version.json 2>/dev/null | grep -E '"(type|major|minor|date)"' || true
