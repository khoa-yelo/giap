#!/bin/bash
#SBATCH --job-name=giap-pilot-cleanup
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# S15, the pilot clean-up (docs/specs/M10-14-pilot-reproducibility.md, "The clean-up"; §3.8, §3.9).
#
#   sbatch scripts/slurm/pilot_cleanup.sbatch            # DRY_RUN=1 (default): plan.tsv only
#   CLEAN=1 sbatch scripts/slurm/pilot_cleanup.sbatch    # executes the plan
#
# Plan step 1: `nextflow clean -f -k` in every pilot/bench run directory with a .nextflow/, then
# its work/ is deleted. Step 4: the mzML directory. Step 2 (the DIA-NN reports and libraries) and
# the de novo outputs (*.ssl, *.mztab) are listed as "kept pending K2" and never deleted (lead
# decisions 2026-10-24); steps 3, 5 and 6 are kept. Everything goes to $GIAP_ROOT/pilot/cleanup/: plan.tsv
# and the append-only ledger DELETED.tsv (a row before each delete, a second row after it, so a
# re-run after a kill continues; a target already gone is `absent`). F is read from the runbook
# header ($RUNBOOK, default $GIAP_REPO/validation/pilot/reproduce.md). Exit 0 done or dry run;
# 1 a deletion or a `nextflow clean` failed; 2 a precondition (nothing deleted).
set -euo pipefail
set -f
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/giap.XXXXXX")"
export TMPDIR
trap 'rm -rf "$TMPDIR"' EXIT

die2() { echo "pilot_cleanup: $*" >&2; exit 2; }

GIAP_PYTHON="${GIAP_PYTHON:-$UV_PROJECT_ENVIRONMENT/bin/python}"
RUNBOOK="${RUNBOOK:-$GIAP_REPO/validation/pilot/reproduce.md}"
DRY_RUN="${DRY_RUN:-1}"
[ "${CLEAN:-0}" = 1 ] && DRY_RUN=0
JOB="${SLURM_JOB_ID:-$$}"
ROOT="$(realpath -e -- "$GIAP_ROOT")" || die2 "GIAP_ROOT=$GIAP_ROOT does not exist"
OUT="$ROOT/pilot/cleanup"
[ -f "$RUNBOOK" ] || die2 "no runbook at $RUNBOOK"

# ---- the freeze commit and the protected paths (Inputs rows, Outputs rows with kept = yes)
"$GIAP_PYTHON" -I "$GIAP_REPO/validation/pilot/reproduce_check.py" protected "$RUNBOOK" \
  --repo "$GIAP_REPO" > "$TMPDIR/protected" || die2 "cannot read $RUNBOOK"
F="$(head -n 1 "$TMPDIR/protected")"
case "$F" in
  *[!0-9a-f]*|"") die2 "the runbook names no freeze commit (freeze_commit=$F)" ;;
esac
[ "${#F}" -eq 40 ] || die2 "the runbook names no freeze commit (freeze_commit=$F)"
F7="$(printf '%s' "$F" | cut -c1-7)"

# ---- no pilot job pending or running (this job excluded)
squeue -u "${USER:-$(id -un)}" -h -o '%i %j' > "$TMPDIR/jobs" || die2 "squeue failed"
while read -r id name; do
  [ -n "$id" ] && [ "$id" != "$JOB" ] || continue
  case "$name" in
    giap-run-*pilot*|giap-run-freeze-*|giap-diann-*|giap-entrapment*|giap-denovo*|giap-pilot-*|giap-bench*)
      die2 "pilot job $id ($name) is pending or running" ;;
  esac
done < "$TMPDIR/jobs"

# ---- the evidence preconditions and the plan (plan.tsv: path, kind, step, action)
mkdir -p "$OUT"
# Any error in this step (an unreadable or malformed JSON included) is a precondition: exit 2.
"$GIAP_PYTHON" -I - "$ROOT" "$F" "${GIAP_DATA_DIR:-$ROOT/data}" "$TMPDIR/protected" \
  "$OUT/plan.tsv" <<'PY' || exit 2
import glob, json, os, sys
root, F, data, protected_file, plan = sys.argv[1:6]
f7 = F[:7]
def die(msg):
    print(f"pilot_cleanup: {msg}", file=sys.stderr); sys.exit(2)
def canon(p):  # ancestors resolved, the last component kept (a link stays a link)
    return os.path.join(os.path.realpath(os.path.dirname(p)), os.path.basename(p))
R = f"{root}/runs/pilot/freeze/{f7}"
ev = sorted(glob.glob(f"{R}/evidence/*.json"))
if not ev:
    die(f"no freeze evidence under {R}/evidence/")
for p in ev:
    if json.load(open(p)).get("failures") != []:
        die(f"{p} has failures")
if not os.path.isfile(f"{root}/pilot/metrics/freeze-{f7}/metrics.json"):
    die(f"no {root}/pilot/metrics/freeze-{f7}/metrics.json")
vj = f"{root}/pilot/reproduce/verify-{f7}.json"
doc = json.load(open(vj))
if doc.get("freeze_commit") != F or doc.get("failures") != [] or doc.get("static_ok") is not True:
    die(f"{vj}: not a passing S-check at {F}")
rec = doc.get("reconvert") or {}
if not rec or any(v not in ("identical", "changed", "absent") for v in rec.values()):
    die(f"{vj}: re-conversion results are not recorded")
rows = []
runs = [d for b in ("pass1", "te7", f"freeze/{f7}") for d in glob.glob(f"{root}/runs/pilot/{b}/*")]
runs += glob.glob(f"{root}/bench/*/r*/*/out")
for d in sorted(runs):
    if os.path.isdir(f"{d}/.nextflow") and not os.path.islink(d):
        rows.append((d, "run", "1", "nextflow_clean"))
        w = f"{d}/work"
        rows.append((w, "dir", "1", "skipped_symlink" if os.path.islink(w) else "delete"))
k2 = ("*.speclib", "*.parquet", "*matrix*", "quant")
for d in (f"{root}/pilot/A1", f"{root}/pilot/A2", f"{root}/pilot/te7/A2"):
    for pat in k2:
        for p in sorted(glob.glob(f"{d}/{pat}")):
            rows.append((p, "dir" if os.path.isdir(p) else "file", "2", "kept_pending_K2"))
for pat in ("A3_cascadia/*/*.ssl", "A3_casanovo/*/*.mztab"):
    rows += [(p, "file", "4", "kept_pending_K2") for p in sorted(glob.glob(f"{root}/pilot/{pat}"))]
mz = canon(f"{data}/pxd078459/mzml")
if os.path.lexists(mz):
    kind = "dir" if os.path.isdir(mz) else "file"
    rows.append((mz, kind, "4", "skipped_symlink" if os.path.islink(mz) else "delete"))
lines = open(protected_file).read().splitlines()[1:]
keep = [os.path.realpath(k) for k in lines if k.strip()]
for p, _k, _s, action in rows:
    if "\t" in p or "\n" in p or not os.path.isabs(p):
        die(f"refusing a plan path that is relative or holds a tab or newline: {p!r}")
    if action not in ("delete", "nextflow_clean"):
        continue
    rp = os.path.realpath(p)
    if not (canon(p).startswith(root + "/") and rp.startswith(root + "/")):
        die(f"target {p} (resolving to {rp}) is not inside {root}")
    for k in keep if action == "delete" else []:
        if k == rp or k.startswith(rp + "/") or rp.startswith(k + "/"):
            die(f"target {p} is (or holds) a runbook Inputs/Outputs row: {k}")
with open(plan, "w") as fh:
    fh.write("path\tkind\tstep\taction\n")
    fh.writelines("\t".join(r) + "\n" for r in rows)
PY

# ---- every target resolves inside $GIAP_ROOT
while IFS="$(printf '\t')" read -r path kind step action; do
  [ "$action" = delete ] || [ "$action" = nextflow_clean ] || continue
  case "$path" in "$ROOT"/?*) ;; *) die2 "$path is not inside $ROOT" ;; esac
  [ -e "$path" ] || continue
  real="$(realpath -e -- "$path")" || die2 "cannot resolve $path"
  case "$real" in "$ROOT"/?*) ;; *) die2 "$path resolves to $real, not inside $ROOT" ;; esac
done < <(tail -n +2 "$OUT/plan.tsv")

echo "== plan ($OUT/plan.tsv)"
tail -n +2 "$OUT/plan.tsv" | awk -F'\t' '{print "step " $3 "\t" $4 "\t" $1}'
if [ "$DRY_RUN" != 0 ]; then
  echo "pilot_cleanup: dry run, nothing deleted (CLEAN=1 executes the plan)"
  exit 0
fi

# ---- execute
LEDGER="$OUT/DELETED.tsv"
[ -s "$LEDGER" ] || printf 'path\tkind\tbytes\tsha256\taction\tutc\tjob\n' > "$LEDGER"
row() { printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$JOB" >> "$LEDGER"; }
echo "== du -sb before: $(du -sb "$ROOT" | cut -f1)"
rc=0
while IFS="$(printf '\t')" read -r path kind step action; do
  case "$action" in
    nextflow_clean)
      row "$path" "$kind" NA NA nextflow_clean
      if (cd "$path" && nextflow clean -f -k); then
        row "$path" "$kind" NA NA nextflow_clean_done
      else
        row "$path" "$kind" NA NA nextflow_clean_failed
        echo "pilot_cleanup: nextflow clean failed in $path" >&2
        rc=1
      fi
      continue ;;
    skipped_symlink) row "$path" "$kind" NA NA skipped_symlink; continue ;;
    delete) ;;
    *) continue ;;
  esac
  if [ -L "$path" ]; then row "$path" "$kind" NA NA skipped_symlink; continue; fi
  if [ ! -e "$path" ]; then row "$path" "$kind" 0 NA absent; continue; fi
  bytes="$(du -sb -- "$path" | cut -f1)"
  sha=NA
  if [ -f "$path" ] && [ "$bytes" -le 10737418240 ]; then
    sha="$(sha256sum -- "$path" | awk '{print $1}')"
  fi
  row "$path" "$kind" "$bytes" "$sha" deleting
  if rm -rf -- "$path"; then
    row "$path" "$kind" "$bytes" "$sha" deleted
  else
    row "$path" "$kind" "$bytes" "$sha" failed
    rc=1
  fi
done < <(tail -n +2 "$OUT/plan.tsv")
echo "== du -sb after: $(du -sb "$ROOT" | cut -f1)"
exit "$rc"
