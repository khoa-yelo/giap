#!/bin/bash
#SBATCH --job-name=giap-fetch-resources
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Prefetch of the §3.2 resources that later steps need (spec §3.2, §3.9; R6): geNOMAD DB v1.9, the
# Promotech RF-HOT model, the K-12 MG1655 and Nissle 1917 genomes, the K-12 UniProt proteome and the
# antiSMASH databases. Official sources only; every file is checksummed and logged in
# $GIAP_DB_DIR/prefetch_manifest.tsv (source, file, checksum, how it was verified). `giap setup`
# (M2-08) later verifies and records them in db_manifest.json. Re-running skips finished items and
# resumes partial downloads (aria2c control files). A failing item is logged and skipped (R6).
set -uo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
cd "$GIAP_REPO"
DB="$GIAP_DB_DIR"; DATA="${GIAP_DATA_DIR:-$GIAP_ROOT/data}"
MAN="$DB/prefetch_manifest.tsv"
mkdir -p "$DB" "$DATA"
[ -f "$MAN" ] || printf 'item\tsource\tfile\tchecksum\tverified\tdate\n' > "$MAN"
log() { printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" "$(date -u +%FT%TZ)" >> "$MAN"; }
fetch() {  # url dir name
  aria2c --continue=true --max-connection-per-server=4 --split=4 --min-split-size=32M --max-tries=0 \
    --retry-wait=30 --timeout=120 --file-allocation=none --auto-file-renaming=false \
    --summary-interval=600 --console-log-level=warn --dir="$2" --out="$3" "$1"
}
md5ok() { echo "$1  $2" | md5sum -c --status; }
fail=0

# 1. geNOMAD DB v1.9 (Zenodo 14886553; md5 from the record)
if [ ! -d "$DB/genomad_db" ]; then  # §7.1.1: DB/genomad_db
  u="https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz?download=1"
  if fetch "$u" "$DB" genomad_db_v1.9.tar.gz && md5ok 67244b528bb8bed464d1ca147136d33e "$DB/genomad_db_v1.9.tar.gz"; then
    tar -xzf "$DB/genomad_db_v1.9.tar.gz" -C "$DB" && rm -f "$DB/genomad_db_v1.9.tar.gz"
    log genomad_db "$u" genomad_db_v1.9.tar.gz md5:67244b528bb8bed464d1ca147136d33e "zenodo md5"
  else echo "BLOCKED genomad_db"; fail=1; fi
fi

# 2. Promotech RF-HOT model (www.cs.mun.ca; no published checksum: sha256 recorded on first download)
if [ ! -d "$DB/promotech/RF-HOT" ] && [ ! -f "$DB/promotech/RF-HOT.zip.sha256" ]; then
  u="https://www.cs.mun.ca/~lourdes/public/PromoTech_models/RF-HOT.zip"
  mkdir -p "$DB/promotech"
  if fetch "$u" "$DB/promotech" RF-HOT.zip && unzip -tq "$DB/promotech/RF-HOT.zip"; then
    s=$(sha256sum "$DB/promotech/RF-HOT.zip" | cut -d' ' -f1); echo "$s" > "$DB/promotech/RF-HOT.zip.sha256"
    unzip -oq "$DB/promotech/RF-HOT.zip" -d "$DB/promotech"
    log promotech_rf_hot "$u" RF-HOT.zip "sha256:$s" "zip integrity; first-download sha256"
  else echo "BLOCKED promotech_rf_hot"; fail=1; fi
fi

# 3. Reference genomes (NCBI FTP; md5 from the assembly's md5checksums.txt)
for spec in "k12_mg1655 GCA_000005845.2_ASM584v2 7e69874199f23fd21b060dc0b2b72321" \
            "nissle1917 GCA_043228065.1_ASM4322806v1 ff062688cf999ade769f4a27c2cb0949"; do
  set -- $spec; name=$1; asm=$2; md5=$3; acc=${asm%%_ASM*}
  p=$(echo "$acc" | sed -E 's/GCA_([0-9]{3})([0-9]{3})([0-9]{3}).*/\1\/\2\/\3/')
  f="${asm}_genomic.fna.gz"; u="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/$p/$asm/$f"
  mkdir -p "$DATA/genomes/$name"
  if [ -f "$DATA/genomes/$name/$f" ] && md5ok "$md5" "$DATA/genomes/$name/$f"; then continue; fi
  if fetch "$u" "$DATA/genomes/$name" "$f" && md5ok "$md5" "$DATA/genomes/$name/$f"; then
    log "genome_$name" "$u" "$f" "md5:$md5" "ncbi md5checksums.txt"
  else echo "BLOCKED genome_$name"; fail=1; fi
done

# 4. K-12 proteome UP000000625 (UniProt REST; release recorded, sha256 on download)
if [ ! -f "$DATA/proteomes/UP000000625.fasta" ]; then
  u='https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000000625%29'
  mkdir -p "$DATA/proteomes"
  rel=$(curl -s -D - -o "$DATA/proteomes/UP000000625.fasta.part" --max-time 1800 "$u" | tr -d '\r' | awk -F': ' 'tolower($1)=="x-uniprot-release"{print $2}')
  n=$(grep -c '^>' "$DATA/proteomes/UP000000625.fasta.part" 2>/dev/null || echo 0)
  if [ -n "$rel" ] && [ "$n" -gt 4000 ]; then
    mv "$DATA/proteomes/UP000000625.fasta.part" "$DATA/proteomes/UP000000625.fasta"
    s=$(sha256sum "$DATA/proteomes/UP000000625.fasta" | cut -d' ' -f1)
    log proteome_UP000000625 "$u" UP000000625.fasta "sha256:$s" "UniProt release $rel; $n entries"
  else echo "BLOCKED proteome_UP000000625 (release=$rel entries=$n)"; fail=1; fi
fi

# 5. antiSMASH databases (dl.secondarymetabolites.org) via the pinned image's own downloader
IMG=$(ls "$GIAP_ROOT"/containers/antismash-standalone@sha256-*.img 2>/dev/null | head -1)
if [ -n "$IMG" ] && [ ! -f "$DB/antismash/.giap-complete" ]; then
  mkdir -p "$DB/antismash"
  if apptainer exec --bind "$DB/antismash" "$IMG" download-antismash-databases --database-dir "$DB/antismash" \
     && apptainer exec --bind "$DB/antismash" "$IMG" antismash --databases "$DB/antismash" --check-prereqs; then
    touch "$DB/antismash/.giap-complete"
    log antismash_dbs "dl.secondarymetabolites.org (download-databases)" "$DB/antismash" "-" "antismash --check-prereqs"
  else echo "BLOCKED antismash_dbs"; fail=1; fi
fi

echo "== manifest"; column -t -s $'\t' "$MAN" | cut -c1-200
du -sh "$DB"/* "$DATA"/* 2>/dev/null
exit $fail
