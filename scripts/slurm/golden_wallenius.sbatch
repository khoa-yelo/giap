#!/bin/bash
#SBATCH --job-name=giap-golden-wallenius
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# The exact Wallenius reference (M8-12b §3.6, §4.2; ADR-0044 decision 1).
# Usage: sbatch golden_wallenius.sbatch [checkout]
#   runs scripts/golden/make_wallenius_exact.py twice (--workers 4) outside the checkout, checks
#   the two files are byte-identical, then writes tests/data/golden/wallenius/{exact_tails.tsv.gz,
#   PROVENANCE.md}. GOLDEN_WALLENIUS_REASON (optional) is recorded as the reason.
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
WT="$(cd "${1:-$GIAP_REPO}" && pwd)"
: "${SLURM_JOB_ID:?run through sbatch: PROVENANCE.md records the job}"
export TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
cd "$WT"
PY="$UV_PROJECT_ENVIRONMENT/bin/python"
GEN=scripts/golden/make_wallenius_exact.py
INPUTS="$WT/tests/data/golden/r/inputs"
DEST="$WT/tests/data/golden/wallenius"
reason="${GOLDEN_WALLENIUS_REASON:-M8-12b: the exact Wallenius reference (ADR-0044 decision 1)}"
echo "checkout=$WT head=$(git rev-parse --short HEAD) node=$(hostname) job=$SLURM_JOB_ID"

for run in 1 2; do
  "$PY" "$GEN" --inputs "$INPUTS" --out "$TMPDIR/exact$run.tsv.gz" --workers 4 \
    --log "$TMPDIR/log$run.json"
done
cmp "$TMPDIR/exact1.tsv.gz" "$TMPDIR/exact2.tsv.gz"
echo "two runs byte-identical"

checks="$("$PY" - "$TMPDIR/log1.json" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
c, env = d.get("checks", {}), d.get("environment", {})
print(f"- Environment: Python {env.get('python')}, numpy {env.get('numpy')}, "
      f"np.longdouble {env.get('longdouble')}, machine {env.get('machine')}.")
print("- Checks (M8-12b §3.4):")
print(f"  - R1 80-bit vs float64, all tails: max rel {c.get('r1_longdouble_vs_float64_max_rel')} (limit 1e-13);")
print(f"  - R1 80-bit vs R2 on the R2 tails ({', '.join(d.get('r2_tails', []))}): "
      f"max rel {c.get('r1_vs_r2_max_rel')} (limit 1e-13);")
spot = c.get("spot_checks", {})
print(f"  - R2 spot checks on {len(spot)} tails >= 1e-3 (limit 1e-13): "
      + "; ".join(f"{k} rel {v}" for k, v in spot.items()) + ";")
sumf = c.get("sum_f", {})
print(f"  - Σf = 1 by R2 on {len(sumf)} categories (limit 1e-50): "
      + "; ".join(f"{k} |Σ - 1| {v}" for k, v in sumf.items()) + ";")
print(f"  - vacuity floor (ADR-0044 decision 4): {c.get('vacuity')};")
print(f"  - ω max ulp difference against the numpy mean ratio: {d.get('omega_max_ulp_vs_numpy')}.")
PY
)"
sha_genes="$(sha256sum "$INPUTS/ora_genes.tsv.gz" | cut -d' ' -f1)"
sha_cats="$(sha256sum "$INPUTS/ora_categories.tsv.gz" | cut -d' ' -f1)"
dirty=""
[ -z "$(git status --porcelain -- "$GEN" scripts/slurm/golden_wallenius.sbatch)" ] ||
  dirty=" (scripts modified after this commit)"

mkdir -p "$DEST"
cp "$TMPDIR/exact1.tsv.gz" "$DEST/exact_tails.tsv.gz"
cat > "$DEST/PROVENANCE.md" <<MD
# Exact Wallenius reference (M8-12b; ADR-0044 decision 1; GIAP_SPEC.md §7.7.5, §11.4)

- Command: \`sbatch scripts/slurm/golden_wallenius.sbatch\` (runs \`$GEN --workers 4\` twice).
- Slurm job: $SLURM_JOB_ID on $(hostname), $(date -u +%Y-%m-%d) ($(date -u +%H:%M) UTC).
- Reason: $reason
- Scripts: git commit \`$(git rev-parse HEAD)\`$dirty.
$checks
- Inputs: \`tests/data/golden/r/inputs/ora_genes.tsv.gz\` sha256 \`$sha_genes\`;
  \`tests/data/golden/r/inputs/ora_categories.tsv.gz\` sha256 \`$sha_cats\`.
- Method: R1 (urn recursion in 80-bit long double) for tails >= 1e-3; R2 (the exact rational
  integral, fixed point with >= 200 surviving bits) for tails < 1e-3; the larger tail is
  1 - small + pmf(k); ω is the exact rational mean ratio, correctly rounded once to float64.
  The reference is ADR-0044 decision 1. No goseq, BiasedUrn or scipy code is used.
- Reproducibility: two runs in this job were byte-identical.
- Regeneration only through this job and only with a reason; a change to the ORA inputs requires
  it. CI fails a change to \`tests/data/golden/wallenius/\` that does not update this file.
MD
du -ab "$DEST"
echo "exact Wallenius reference written to $DEST"
