#!/bin/bash
#SBATCH --job-name=gippa-diann-image
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Verify Khoa's DIA-NN 2.3.1 Academia image and install the read-only copy every DIA-NN job runs
# (docs/specs/M10-02-diann-searches.md "Image verification"; index lead decision 2; ADR-0051).
#
#   sbatch scripts/slurm/diann_image.sbatch
#
# DIANN_SIF_SHARED is the source (world-writable, so never run); DIANN_SIF is the copy. The pin is
# read from validation/pilot/diann_image.sha256 only. The copy, not the source, is hashed; an
# existing copy is never replaced (a tampered one is evidence). The banner is checked on the
# temporary copy, so a copy is installed only after every check passed (PR #293 review S3). Exit 0
# installed or already installed, 2 any precondition or mismatch; nothing is retried automatically.
set -euo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
module load apptainer/1.4.4 >/dev/null 2>&1 || true

PIN_FILE="$GIPPA_REPO/validation/pilot/diann_image.sha256"
SHARED="${DIANN_SIF_SHARED:-/farmshare/home/classes/bios/270/jordan/diann_2.3.1.sif}"
SIF="${DIANN_SIF:-$GIPPA_ROOT/images/diann_2.3.1.sif}"
IMAGES="$(dirname "$SIF")"
PY="${GIPPA_PYTHON:-$UV_PROJECT_ENVIRONMENT/bin/python}"
BANNER_WANT="DIA-NN 2.3.1 Academia"

die() { echo "diann_image: $*" >&2; exit 2; }
sha() { sha256sum "$1" | cut -d' ' -f1; }
read_only() { local m; m=$(stat -c %a "$1"); [ $((0$m & 0222)) -eq 0 ]; }

# 1. the pin: one line, "<64 hex>  diann_2.3.1.sif"
[ -f "$PIN_FILE" ] || die "missing pin file $PIN_FILE (validation/pilot/diann_image.sha256)"
if [ "$(grep -c '' "$PIN_FILE")" -ne 1 ] || ! grep -Eqx '[0-9a-f]{64}  diann_2\.3\.1\.sif' "$PIN_FILE"; then
  die "malformed pin file $PIN_FILE: want one line '<sha256>  diann_2.3.1.sif'"
fi
PIN="$(cut -d' ' -f1 "$PIN_FILE")"

# 2. an existing copy: accepted only when it matches the pin and is read-only; never replaced
if [ -e "$SIF" ]; then
  if [ -f "$SIF" ] && [ "$(sha "$SIF")" = "$PIN" ] && read_only "$SIF"; then
    echo "already installed: $SIF"
    exit 0
  fi
  die "the installed copy $SIF does not match the pin $PIN (or is writable); inspect it and delete it by hand"
fi

# 3. the source
if [ ! -f "$SHARED" ] || [ ! -r "$SHARED" ]; then
  die "the shared image $SHARED is not a readable file; ask Khoa for the authentic 2.3.1 Academia image; never fetch DIA-NN elsewhere (R6)"
fi

# 4. space for the copy
mkdir -p -m 0755 "$IMAGES"
need_kb=$(( $(stat -c %s "$SHARED") * 2 / 1024 ))
free_kb=$(df -Pk "$IMAGES" | awk 'NR == 2 { print $4 }')
[ "$free_kb" -ge "$need_kb" ] || die "only ${free_kb} KB free under $IMAGES; the copy needs ${need_kb} KB (2 x the image)"

# 5. copy, then hash the copy (the source may change during the copy)
TMP="$IMAGES/.$(basename "$SIF").tmp.${SLURM_JOB_ID:-$$}"
trap 'rm -f "${TMP:?}"' EXIT
cp --preserve=timestamps "$SHARED" "$TMP"
COPY_SHA="$(sha "$TMP")"
if [ "$COPY_SHA" != "$PIN" ]; then
  rm -f "${TMP:?}"
  die "the shared image no longer matches the pin: $SHARED copied with sha256 $COPY_SHA, the pin ($PIN_FILE) is $PIN; ask Khoa for the authentic 2.3.1 Academia image; never fetch DIA-NN elsewhere (R6)"
fi

# 6. read-only, then the banner check on the temporary copy: a failed banner installs nothing
chmod 0444 "$TMP"
BANNER="$(apptainer exec "$TMP" /diann-2.3.1/diann-linux --version 2>&1 || true)"
# the binary prints a blank line before its banner: compare the first non-empty line
BANNER1="$(printf '%s\n' "$BANNER" | grep -m 1 . || true)"
case "$BANNER1" in
  *"$BANNER_WANT"*) ;;
  *) rm -f "${TMP:?}"; die "the copy's banner is not '$BANNER_WANT': $BANNER1; nothing installed" ;;
esac
LABELS="$(apptainer inspect --json "$TMP" 2>/dev/null || true)"
APPTAINER_VERSION="$(apptainer --version 2>/dev/null || true)"

# 7. an atomic rename that never overwrites; re-check whatever is there now
mv -n "$TMP" "$SIF" || true
rm -f "${TMP:?}"
if [ "$(sha "$SIF")" != "$PIN" ] || ! read_only "$SIF"; then
  die "the installed copy $SIF does not match the pin $PIN (or is writable); inspect it and delete it by hand"
fi

# 8. the install record
"$PY" -I - "$SHARED" "$SIF" "$PIN" "$PIN_FILE" "$COPY_SHA" "$(hostname)" "$LABELS" \
  "$APPTAINER_VERSION" "$BANNER" "${SLURM_JOB_ID:-}" >"$SIF.json" <<'PY'
import grp, json, os, pwd, stat, sys, time

src, copy, pin, pin_file, sha, host, labels, version, banner, job = sys.argv[1:]

def record(path, **extra):
    st = os.stat(path)
    return {"path": path, "size": st.st_size, "mode": format(stat.S_IMODE(st.st_mode), "o"),
            **extra}

source = record(src, sha256=sha)
st = os.stat(src)
source.update(owner=pwd.getpwuid(st.st_uid).pw_name, group=grp.getgrgid(st.st_gid).gr_name,
              mtime=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(st.st_mtime)))
try:
    parsed = json.loads(labels)["data"]["attributes"]["labels"]
except (ValueError, KeyError, TypeError):
    parsed = labels
print(json.dumps({
    "source": source, "copy": record(copy, sha256=sha), "pin": pin, "pin_file": pin_file,
    "installed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "job_id": job,
    "hostname": host, "labels": parsed, "apptainer_version": version,
    "banner": [line for line in banner.splitlines() if line.strip()][:2],
}, indent=2))
PY
echo "installed: $SIF (sha256 $PIN, mode 0444); record $SIF.json"
