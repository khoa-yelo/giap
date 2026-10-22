#!/bin/bash
#SBATCH --job-name=giap-record-ostir-novel
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# M7-08 rule 9: one real OSTIR 1.1.3 run on a multi-start window of three sim100k novel
# candidates, for tests/contract/test_ostir_range.py. Input from scripts/dev/record_ostir_novel_input.py
# (novel_candidates' ostir_input_novel.csv on the M7-08 sim100k); the command is the OSTIR module's
# (workflow/modules/local/ostir/main.nf) in the official micromamba env giap-ostir
# (workflow/envs/ostir.yml + lock), as record_tools.sbatch runs it.
# Output: $GIAP_ROOT/recordings/ostir/sim100k_novel/ with ostir_input.csv, ostir.csv, command.txt,
# versions.txt, run.log, runtime.txt, SHA256SUMS. Usage: sbatch record_ostir_novel.sbatch <worktree>
set -euo pipefail
WT="${1:?worktree dir}"
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
export TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
export NUMBA_CACHE_DIR="$TMPDIR/numba"
N="${SLURM_CPUS_PER_TASK:-4}"
D="$GIAP_ROOT/recordings/ostir/sim100k_novel"
rm -rf "$D"; mkdir -p "$D"; cd "$D"
echo "worktree=$WT head=$(git -C "$WT" rev-parse --short HEAD) node=$(hostname)"
PYTHONPATH="$WT/src:$WT/tests/sim:$WT/tests" "$UV_PROJECT_ENVIRONMENT/bin/python" \
  "$WT/scripts/dev/record_ostir_novel_input.py" "$D"
c=(ostir -i ostir_input.csv -o ostir.csv -j "$N")
{ echo "# ostir_input.csv: scripts/dev/record_ostir_novel_input.py (giap $(git -C "$WT" rev-parse --short HEAD))"
  printf '%q ' "${c[@]}"; echo; } > command.txt
t0=$(date +%s)
micromamba run -n giap-ostir "${c[@]}" 2>&1 | tee run.log
echo "wall_s=$(( $(date +%s) - t0 )) cpus=$N node=$(hostname)" > runtime.txt
{ echo "ostir $(micromamba run -n giap-ostir python -c 'import importlib.metadata as m; print(m.version("ostir"))')"
  (cd "$WT/workflow/envs/locks" && sha256sum ostir.explicit.txt ostir.pip.txt); } > versions.txt
find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
cat versions.txt runtime.txt
