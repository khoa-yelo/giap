#!/bin/bash
#SBATCH --job-name=gippa-fetch-a3-assets
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# The external A3 assets (M10-08 deliverable 6; ADR-0054 = ADR-M10-B). The only A3 script that
# reaches the network.
#
#   A3_WEIGHTS_HOST_OK=1 sbatch scripts/slurm/fetch_a3_assets.sbatch --resolve   # 1. print sha256
#   A3_WEIGHTS_HOST_OK=1 sbatch scripts/slurm/fetch_a3_assets.sbatch             # 2. verified copy
#
# - Cascadia weights: HTTPS from the Drive file ID in validation/pilot/a3_assets.tsv to
#   $GIPPA_DB_DIR/cascadia/cascadia_0.0.7.ckpt (0444). Drive is not a §3.2 source: only once
#   ADR-M10-B is accepted (A3_WEIGHTS_HOST_OK=1 exactly), otherwise exit 2 before any request.
#   The body must be a zip (Lightning checkpoint, PK\x03\x04), never an HTML page, in both modes.
#   curl: 2 attempts x 3000 s + 60 s delay = 6060 s < the 2 h --time.
# - Unimod table: psims' vendored unimod_tables.xml.gz, gunzipped from the built gippa-cascadia
#   env to $GIPPA_DB_DIR/cascadia/unimod_tables.xml (0444). No network, so it runs without the
#   weights gate.
# - --resolve prints each asset's sha256 and bytes (trust on first use: no checksum is
#   published) for the lead to commit into a3_assets.tsv; it installs nothing. Without --resolve
#   every file must match the committed sha256 (R6), else exit 1 and no file at the final path.
# - No image is pulled: route B1 (ProteoWizard, five more vendor EULAs) is not accepted.
# R6: https only, TLS verification never off. R11: nothing is written into the checkout.
# Exit codes: 0 done; 1 a download or verification failed; 2 a precondition.
set -euo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"

ASSETS="$GIPPA_REPO/validation/pilot/a3_assets.tsv"
DEST="$GIPPA_DB_DIR/cascadia"
CASCADIA_ENV="$MAMBA_ROOT_PREFIX/envs/gippa-cascadia"
RESOLVE=0
case "${1:-}" in
  "") ;;
  --resolve) RESOLVE=1 ;;
  *) echo "usage: fetch_a3_assets.sbatch [--resolve]" >&2; exit 2 ;;
esac
[ -f "$ASSETS" ] || { echo "no $ASSETS" >&2; exit 2; }

field() {  # field ASSET COLUMN: one cell of a3_assets.tsv
  awk -F'\t' -v a="$1" -v c="$2" \
    'NR == 1 { for (i = 1; i <= NF; i++) if ($i == c) k = i; next } $1 == a { print $k; exit }' \
    "$ASSETS"
}
sha256_of() {  # GNU sha256sum or the BSD shasum of macOS CI
  if command -v sha256sum >/dev/null 2>&1; then sha256sum -- "$1"; else shasum -a 256 -- "$1"; fi \
    | cut -d' ' -f1
}

# install_asset ASSET PART FINAL [zip]: verify PART against the committed sha256 (and, with zip,
# the zip magic), then move it into place.
install_asset() {
  local asset="$1" part="$2" final="$3" want got bytes
  if [ "${4:-}" = zip ] && [ "$(head -c 4 "$part" | od -An -tx1 | tr -d ' \n')" != 504b0304 ]; then
    rm -f "$part"
    echo "$asset: Drive returned a page, not the file (no zip magic): nothing installed" >&2
    return 1
  fi
  want="$(field "$asset" sha256)"
  got="$(sha256_of "$part")"
  bytes="$(wc -c < "$part" | tr -d ' ')"
  echo "$asset sha256 $got bytes $bytes"
  if [ "$RESOLVE" -eq 1 ]; then
    rm -f "$part"
    return 0
  fi
  if [ "$got" != "$want" ]; then
    rm -f "$part"
    echo "$asset: sha256 $got does not match a3_assets.tsv ($want): nothing installed" >&2
    return 1
  fi
  chmod 0444 "$part"
  mv -f "$part" "$final"
  echo "$asset -> $final (0444, retrieved $(date -u +%Y-%m-%d))"
}

needs_sha() {  # without --resolve, an asset with no committed sha256 is never fetched (R6)
  local want
  want="$(field "$1" sha256)"
  if [ "$RESOLVE" -eq 0 ] && ! [[ "$want" =~ ^[0-9a-f]{64}$ ]]; then
    echo "$1: no sha256 in a3_assets.tsv: run with --resolve and commit it first" >&2
    return 1
  fi
}

mkdir -p "$DEST"
rc=0
u_part="" w_part=""
cleanup() { local p; for p in "$u_part" "$w_part"; do [ -z "$p" ] || rm -f -- "$p"; done; }
trap cleanup EXIT
newpart() { mktemp "$DEST/.$1.part.XXXXXX"; }  # unique per job (N8)
fail() { [ "$rc" -ne 0 ] || rc="$1"; }  # the first failure decides the exit code

# Unimod table (offline, from the env build).
if needs_sha unimod_table; then
  gz="$(find "$CASCADIA_ENV" -path '*/psims/controlled_vocabulary/vendor/unimod_tables.xml.gz' \
    2>/dev/null | head -n 1 || true)"
  if [ -z "$gz" ]; then
    echo "no psims unimod_tables.xml.gz under $CASCADIA_ENV: build_envs.sbatch cascadia" >&2
    exit 2
  fi
  echo "== unimod_table from $gz"
  u_part="$(newpart unimod_tables.xml)"; part="$u_part"
  if gzip -dc "$gz" > "$part"; then
    install_asset unimod_table "$part" "$DEST/unimod_tables.xml" || fail 1
  else
    rm -f "$part"
    fail 1
  fi
else
  fail 2
fi

# Cascadia weights (HTTPS, Drive).
if [ "${A3_WEIGHTS_HOST_OK:-}" != 1 ]; then
  echo "weights host not in §3.2: ADR-M10-B (ADR-0054) must be accepted first;" \
    "then set A3_WEIGHTS_HOST_OK=1" >&2
  fail 2
  exit "$rc"
fi
url="$(field cascadia_weights url)"
case "$url" in https://*) ;; *) echo "cascadia_weights: not an https URL: $url" >&2; exit 2 ;; esac
if needs_sha cascadia_weights; then
  w_part="$(newpart cascadia_0.0.7.ckpt)"; part="$w_part"
  echo "== cascadia_weights from $url"
  if info="$(curl --fail --location --silent --show-error --proto '=https' --proto-redir '=https' \
      --retry 1 --retry-delay 60 --max-time 3000 --write-out '%{content_type} %{size_download}' \
      --output "$part" "$url")"; then
    echo "cascadia_weights: content type and bytes: $info"
    case "$info" in
      text/html*) rm -f "$part"
        echo "cascadia_weights: Drive returned a page, not the file: nothing installed" >&2
        fail 1 ;;
      *) install_asset cascadia_weights "$part" "$DEST/cascadia_0.0.7.ckpt" zip || fail 1 ;;
    esac
  else
    rm -f "$part"
    echo "cascadia_weights: download failed: retry once later, then a Blocked row (R6)" >&2
    fail 1
  fi
else
  fail 2
fi
ls -l "$DEST"
exit "$rc"
