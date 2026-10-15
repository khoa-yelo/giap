#!/bin/bash
#SBATCH --job-name=giap-annotate-real
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Real-tool ANNOTATE run on SIM100K (M2-07 review decision 10): main.nf with inputs.replicons and
# db_dir (the `giap run` wiring is the follow-up PR), -profile farmshare_local,apptainer, work dir
# under $GIAP_ROOT/scratch (never node-local /tmp). Bakta is compared with the recording
# $GIAP_ROOT/recordings/bakta_full/sim100k.
# Until `giap _step` loads the step modules and the M2-13..M2-18 steps exist, a `giap` wrapper
# first on PATH loads them and stands in for unregistered py steps: tool inputs from
# scripts/dev/record_tools_inputs.py on the recorded Bakta JSON, empty files for the layer steps.
# Usage: sbatch [--export=ALL,WT=<worktree>] annotate_real.sbatch
set -uo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
module load apptainer 2>/dev/null || true
WT="${WT:-$GIAP_REPO/local/worktrees/nf-tools}"
export PYTHONPATH="$WT/src:$WT/tests/sim" PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH" NXF_ANSI_LOG=false
REC="$GIAP_ROOT/recordings/bakta_full/sim100k"
PRE="$REC/preflight/work/preflight"
W="$GIAP_ROOT/scratch/annotate-real/${SLURM_JOB_ID:-manual}"
mkdir -p "$W/bin" && cd "$W" || exit 1
export TMPDIR="$W/tmp" GIAP_HARNESS_BAKTA_JSON="$REC/bakta/sim100k.json" GIAP_HARNESS_REPO="$WT"
mkdir -p "$TMPDIR"
echo "worktree=$WT ($(git -C "$WT" rev-parse --short HEAD)) work=$W"

cat > "$W/bin/giap" <<EOF
#!$UV_PROJECT_ENVIRONMENT/bin/python
import sys
sys.path.insert(0, "$W/bin")
from harness_step import main
main()
EOF
chmod +x "$W/bin/giap"
cat > "$W/bin/harness_step.py" <<'EOF'
import os, sys
from pathlib import Path

def _kv(args, flag):
    return [a.split("=", 1) for f, a in zip(args, args[1:]) if f == flag]

def placeholder(name, args):
    sys.path.insert(0, os.environ["GIAP_HARNESS_REPO"] + "/scripts/dev")
    import record_tools_inputs as rti
    js = Path(os.environ["GIAP_HARNESS_BAKTA_JSON"])
    for key, path in _kv(args, "--out"):
        out = Path(path.replace("*", "000"))
        out.parent.mkdir(parents=True, exist_ok=True)
        if key == "genes_coords":
            rti.write_coords(js, out)
        elif key == "ostir_input":
            rti.write_ostir(js, out)
        elif key == "windows":   # one 40-nt window upstream of every gene start
            seqs, genes = rti.load(js)
            with out.open("w") as fh:
                for g in genes:
                    nt = seqs[rti.contig(g)]["nt"]
                    if g["strand"] == "+":
                        w = nt[max(g["start"] - 41, 0): g["start"] - 1]
                    else:
                        w = rti.revcomp(nt[g["stop"]: g["stop"] + 40])
                    if len(w) == 40:
                        fh.write(f">{g['locus']}\n{w.upper()}\n")
        else:
            out.touch()
    print(f"harness: placeholder for unregistered step {name}", file=sys.stderr)

def main():
    args = sys.argv[1:]
    if args[:1] == ["_step"]:
        from giap.steps import load_step_modules
        from giap.steps.registry import all_steps
        load_step_modules()
        if args[1] not in {s.name for s in all_steps()}:
            placeholder(args[1], args[2:])
            return
    from giap.cli import main as cli
    cli()
EOF
export PATH="$W/bin:$PATH"

python - "$W" "$PRE" <<'PY'
import json, os, sys
from pathlib import Path
from giap.config import resolve_config
from giap.runner.nextflow import write_params
w, pre = Path(sys.argv[1]), Path(sys.argv[2])
params, _ = write_params(w / "out", resolve_config(), {"genome": pre / "genome.norm.fa", "replicons": pre / "replicons.tsv"})
p = json.loads(params.read_text()); p["db_dir"] = os.environ["GIAP_DB_DIR"]
params.write_text(json.dumps(p, indent=2, sort_keys=True) + "\n")
PY

SECONDS=0
nextflow run "$WT/workflow/main.nf" -params-file "$W/out/work/params.json" \
  -profile farmshare_local,apptainer -work-dir "$W/work" -with-trace "$W/trace.tsv"
echo "nextflow exit=$? wall=${SECONDS}s"
cut -f4,6,7,10 "$W/trace.tsv"

python - "$W/trace.tsv" "$REC/bakta/sim100k.json" <<'PY'
import csv, json, sys
from pathlib import Path
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
bakta = [r for r in rows if r["name"].startswith("ANNOTATE:BAKTA") and r["status"] == "COMPLETED"]
if not bakta:
    sys.exit("BAKTA did not complete")
got = json.loads(Path(bakta[0]["workdir"], "bakta", "genome.json").read_text())
want = json.loads(Path(sys.argv[2]).read_text())
drop = {"id"}
f = lambda d: sorted(json.dumps({k: v for k, v in x.items() if k not in drop}, sort_keys=True) for x in d["features"])
g, r = f(got), f(want)
print(f"bakta features: run {len(g)}, recording {len(r)}, equal={g == r}")
print(f"bakta sequences equal: {[s['nt'] for s in got['sequences']] == [s['nt'] for s in want['sequences']]}")
diff = sorted(set(g) ^ set(r))
for d in diff[:10]:
    print("  diff:", d[:300])
PY
