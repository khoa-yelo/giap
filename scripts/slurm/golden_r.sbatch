#!/bin/bash
#SBATCH --job-name=giap-golden-r
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# R golden references (M8-12; GIAP_SPEC.md §3.8, §11.5). Usage: sbatch golden_r.sbatch [checkout]
#   1. builds env giap-golden-r from scripts/golden/golden_r.yml when its .giap-built marker
#      (sha256 of the yml, the build_envs.sh pattern) is missing or stale; conda-forge + bioconda only;
#   2. checks packageVersion() of each package equals its pin, exports golden_r.explicit.txt;
#   3. writes the inputs and runs make_golden_r.R twice, checking byte-identical outputs;
#   4. writes tests/data/golden/r/{inputs,outputs,VERSIONS.txt,PROVENANCE.md} in <checkout>.
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
WT="$(cd "${1:-$GIAP_REPO}" && pwd)"
export TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
cd "$WT"
PY="$UV_PROJECT_ENVIRONMENT/bin/python"
ENV=giap-golden-r
YML="$WT/scripts/golden/golden_r.yml"
LOCK="$WT/scripts/golden/golden_r.explicit.txt"
DEST="$WT/tests/data/golden/r"
prefix="$MAMBA_ROOT_PREFIX/envs/$ENV"
marker="$prefix/.giap-built"
want="$(sha256sum "$YML" | cut -d' ' -f1)"
: "${SLURM_JOB_ID:?run through sbatch: PROVENANCE.md records the job}"
echo "checkout=$WT head=$(git rev-parse --short HEAD) node=$(hostname) job=${SLURM_JOB_ID:-none}"

if [ -d "$prefix" ] && [ "$(cat "$marker" 2>/dev/null || true)" = "$want" ]; then
  echo "$ENV exists"
else
  [ -e "$prefix" ] && { echo "rebuilding $ENV (marker missing or stale)"; rm -rf "$prefix"; }
  # One retry if conda.anaconda.org is unreachable (R6: no other source).
  micromamba create -y -n "$ENV" -c conda-forge -c bioconda --override-channels -f "$YML" ||
    { sleep 60; micromamba create -y -n "$ENV" -c conda-forge -c bioconda --override-channels -f "$YML"; }
fi
micromamba run -n "$ENV" Rscript -e '
pins <- c(fgsea = "1.36.2", goseq = "1.62.0", BiasedUrn = "2.0.12")
for (p in names(pins)) if (as.character(packageVersion(p)) != pins[[p]])
  stop(p, " ", as.character(packageVersion(p)), " != pin ", pins[[p]])
if (paste(R.version$major, R.version$minor, sep = ".") != "4.5.3") stop("R is not 4.5.3")
cat("versions match the pins\n")'
micromamba env export -n "$ENV" --explicit > "$LOCK.tmp"
mv "$LOCK.tmp" "$LOCK"
echo "$want" > "$marker"
lock_sha="$(sha256sum "$LOCK" | cut -d' ' -f1)"

seed="$("$PY" -c 'import sys, yaml; print(yaml.safe_load(open(sys.argv[1]))["seeds"]["GOLDEN_SEED"])' \
  tests/calibration/tolerances.yaml)"
"$PY" scripts/golden/make_golden_inputs.py --out "$TMPDIR/in1" --seed "$seed"
"$PY" scripts/golden/make_golden_inputs.py --out "$TMPDIR/in2" --seed "$seed"
diff -r "$TMPDIR/in1" "$TMPDIR/in2"
for run in 1 2; do
  micromamba run -n "$ENV" Rscript --vanilla scripts/golden/make_golden_r.R \
    "$TMPDIR/in1" "$TMPDIR/out$run" "$TMPDIR/VERSIONS$run.txt" "$LOCK" "$seed"
done
diff -r "$TMPDIR/out1" "$TMPDIR/out2"
echo "inputs and outputs byte-identical across two runs"
size="$(du -cb "$TMPDIR"/out1/* | tail -1 | cut -f1)"
[ "$size" -le $((5 * 1024 * 1024)) ] || { echo "outputs exceed 5 MB (R11)" >&2; exit 1; }

rm -rf "$DEST/inputs" "$DEST/outputs"
mkdir -p "$DEST"
cp -r "$TMPDIR/in1" "$DEST/inputs"
cp -r "$TMPDIR/out1" "$DEST/outputs"
cp "$TMPDIR/VERSIONS1.txt" "$DEST/VERSIONS.txt"
dirty=""
[ -z "$(git status --porcelain -- scripts/golden scripts/slurm/golden_r.sbatch)" ] ||
  dirty=" (scripts modified after this commit)"
cat > "$DEST/PROVENANCE.md" <<MD
# R golden references: fgsea and goseq (M8-12; GIAP_SPEC.md §11.4, §11.5)

- Command: \`sbatch scripts/slurm/golden_r.sbatch\` (runs \`scripts/golden/make_golden_inputs.py
  --seed $seed\` and \`Rscript --vanilla scripts/golden/make_golden_r.R\` in env \`$ENV\`).
- Slurm job: ${SLURM_JOB_ID:?run through sbatch} on $(hostname), $(date -u +%Y-%m-%d) ($(date -u +%H:%M) UTC).
- Env lock: \`scripts/golden/golden_r.explicit.txt\`, sha256 \`$lock_sha\`
  (R $(micromamba run -n "$ENV" Rscript -e 'cat(paste(R.version$major, R.version$minor, sep="."))'),
  fgsea 1.36.2, goseq 1.62.0, BiasedUrn 2.0.12; conda-forge and bioconda only).
- Scripts: git commit \`$(git rev-parse HEAD)\`$dirty.
- Seed: GOLDEN_SEED = $seed (\`tests/calibration/tolerances.yaml\`); \`set.seed(GOLDEN_SEED)\`
  before \`fgseaMultilevel\` and before \`goseq\`; \`nproc = 1\`.
- Reproducibility: the inputs and the R outputs were each generated twice in this job and were
  byte-identical (M8-12 examples 2 and 3).
- Inputs are synthetic (random IDs and memberships): no KEGG, GO or COG content. Only numbers
  are committed; no R package source.
- Regeneration only through this job or \`.github/workflows/golden-r.yml\`, with a reason (a
  package upgrade needs an ADR); CI fails a change to \`tests/data/golden/r/\` that does not
  update this file.
MD
du -ab "$DEST"
echo "golden R data written to $DEST"
