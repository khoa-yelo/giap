#!/bin/bash
#SBATCH --job-name=giap-fetch-pxd078459
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# The three wild-type runs of PXD078459 (GIAP_SPEC.md §15.1, §15.2 step 4; M3-24, M10-01).
#
#   sbatch --qos=dev --cpus-per-task=1 --mem=2G --time=00:30:00 \
#     scripts/slurm/fetch_pxd078459.sbatch --resolve   # 1. the record -> a draft manifest
#   sbatch scripts/slurm/fetch_pxd078459.sbatch        # 2. download what the manifest lists
#
# 1. --resolve fetches the ProteomeXchange record into $GIAP_DATA_DIR/pxd078459/record/, records
#    its sha256 and UTC retrieval date, selects the wild-type runs (validation/pilot/px_record.py)
#    and writes record/pxd078459_files.draft.tsv. It downloads no raw data. The agent reviews the
#    draft and commits it as validation/pilot/inputs/pxd078459_files.tsv.
# 2. The download reads that committed manifest (PXD078459_FILES with PXD078459_BASE_URL, or
#    PXD078459_BASE_URL alone, override it), checks the free space, downloads with curl (resuming
#    partial files), verifies the record checksums (one retry), writes SHA256SUMS and the filled
#    manifest next to the data, and unpacks every run to $GIAP_DATA_DIR/pxd078459/<run>.d/.
#
# R6 / M3-24 decision 11: https or ftps only, never TLS verification off, official hosts only.
# R11: raw MS data never go into the repository.
# Exit codes: 0 done; 1 a download, verification or unpack failed; 2 a precondition.
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"

ACCESSION=PXD078459
RECORD_URL="https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=$ACCESSION&outputMode=XML"
DEST="$GIAP_DATA_DIR/pxd078459"
PILOT="$GIAP_REPO/validation/pilot"
MANIFEST="$PILOT/inputs/pxd078459_files.tsv"

sha256_of() {  # GNU sha256sum or the BSD shasum of macOS CI; both print "<hash>  <name>"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum -- "$@"
  else
    shasum -a 256 -- "$@"
  fi
}

mkdir -p "$DEST"
echo "== $ACCESSION -> $DEST"

if [ "${1:-}" = "--resolve" ]; then
  REC="$DEST/record"
  mkdir -p "$REC"
  rm -f "$REC/$ACCESSION.xml.part" "$REC/selection.part" "$REC/pxd078459_files.draft.tsv"
  echo "== record: $RECORD_URL"
  if ! curl --fail --location --silent --show-error --proto '=https' --proto-redir '=https' \
      --retry 1 --retry-delay 60 --max-time 600 --output "$REC/$ACCESSION.xml.part" "$RECORD_URL"; then
    rm -f "$REC/$ACCESSION.xml.part"
    echo "the ProteomeXchange record is unreachable: add a Blocked row (R6), retry once later" >&2
    exit 1
  fi
  mv "$REC/$ACCESSION.xml.part" "$REC/$ACCESSION.xml"
  record_sha=$(sha256_of "$REC/$ACCESSION.xml" | cut -d' ' -f1)
  printf 'url\t%s\nretrieved_utc\t%s\nsha256\t%s\n' "$RECORD_URL" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    "$record_sha" > "$REC/$ACCESSION.xml.info"
  cat "$REC/$ACCESSION.xml.info"
  rc=0
  python3 -I "$PILOT/px_record.py" "$REC/$ACCESSION.xml" --out "$REC/selection.part" || rc=$?
  [ "$rc" -eq 0 ] || exit "$rc"
  python3 -I "$PILOT/pxd_fetch.py" https "$REC/selection.part" \
    --out "$REC/pxd078459_files.draft.tsv" || rc=$?
  rm -f "$REC/selection.part"
  echo "== draft manifest: $REC/pxd078459_files.draft.tsv (review it, then commit it as $MANIFEST)"
  exit "$rc"
fi

rc=0
python3 -I "$PILOT/pxd_fetch.py" download --dest "$DEST" --manifest "$MANIFEST" \
  --files "${PXD078459_FILES:-}" --base-url "${PXD078459_BASE_URL:-}" || rc=$?
if [ -f "$DEST/SHA256SUMS" ]; then
  echo "== SHA256SUMS"
  cat "$DEST/SHA256SUMS"
fi
echo "== $DEST"
du -sh "$DEST" 2>/dev/null || true
ls -1 "$DEST"
if [ "$rc" -eq 0 ]; then
  echo "== next: python3 -I validation/pilot/acquisition_mode.py $DEST/WT_{1,2,3}.d --out validation/pilot/inputs/acquisition.json (dev QOS)"
fi
exit "$rc"
