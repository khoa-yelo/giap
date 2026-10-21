#!/bin/bash
#SBATCH --job-name=giap-fetch-pxd078459
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Download the three wild-type runs of PXD078459 (GIAP_SPEC.md §15.1, §15.2 step 4; M3-24).
#
#   sbatch scripts/slurm/fetch_pxd078459.sbatch
#
# About 5-10 GB of timsTOF (diaPASEF) data go to $GIAP_DATA_DIR/pxd078459/, never into the
# repository (R11). Every downloaded file is checksummed into SHA256SUMS next to it.
#
# URLs (M3-24 decision 11): https only, never plain ftp:// and never with TLS verification off
# (R6). The exact file names of the deposition could not be verified when this script was written
# (M3 has no network), so they are UNVERIFIED defaults: before the M10 run, confirm them against
# the ProteomeXchange record and override PXD078459_BASE_URL / PXD078459_FILES if they differ.
# The run also has to confirm the acquisition mode of the first run (§15.1: diaPASEF is expected).
#
# Environment: PXD078459_BASE_URL, PXD078459_FILES (space separated), GIAP_DATA_DIR.
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"

ACCESSION=PXD078459
RECORD="https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=$ACCESSION"
BASE_URL="${PXD078459_BASE_URL:-https://download.iprox.cn/dataset/$ACCESSION}"
FILES="${PXD078459_FILES:-WT_1.d.zip WT_2.d.zip WT_3.d.zip}"
DEST="$GIAP_DATA_DIR/pxd078459"

echo "== $ACCESSION -> $DEST"
echo "== record: $RECORD"
echo "== UNVERIFIED file names (decision 11): confirm them on the record above before trusting this run"
case "$BASE_URL" in
  https://* | ftps://*) ;;
  *) echo "refusing a non-TLS source: $BASE_URL (R6)" >&2; exit 2 ;;
esac

sha256_of() {  # GNU sha256sum or the BSD shasum of macOS CI; both print "<hash>  <name>"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum -- "$@"
  else
    shasum -a 256 -- "$@"
  fi
}

mkdir -p "$DEST"
NAMES=()
NAME_N=0
for name in $FILES; do
  url="$BASE_URL/$name"
  if [ -s "$DEST/$name" ]; then
    echo "already there, skipping: $DEST/$name"
  else
    echo "downloading $url"
    curl --fail --location --silent --show-error --remove-on-error --retry 5 --retry-delay 10 --max-time 86400 \
      --output "$DEST/$name" "$url"
  fi
  NAMES+=("$name")
  NAME_N=$((NAME_N + 1))
done
[ "$NAME_N" -gt 0 ] || { echo "PXD078459_FILES is empty: nothing to download" >&2; exit 2; }

(cd "$DEST" || exit 1; sha256_of "${NAMES[@]}" > SHA256SUMS)
echo "== SHA256SUMS"
cat "$DEST/SHA256SUMS"

# The runs are needed as *.d directories by scripts/slurm/diann_pilot.sbatch; unpack the archives
# when they really are archives, and say so plainly when they are not.
for archive in "$DEST"/*.zip; do
  [ -e "$archive" ] || continue
  if unzip -tq "$archive" >/dev/null 2>&1; then
    unzip -oq "$archive" -d "$DEST"
  else
    echo "warning: $archive is not a readable zip; unpack the runs by hand before the DIA-NN job" >&2
  fi
done

echo "== $DEST"
du -sh "$DEST" 2>/dev/null || true
ls -1 "$DEST"
echo "== next: confirm the acquisition mode of the first run (§15.2 step 4), then sbatch scripts/slurm/diann_pilot.sbatch A1"
