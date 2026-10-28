#!/bin/bash
#SBATCH --job-name=gippa-cal-suite
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --signal=B:USR1@600
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# One full-tier calibration suite as one Slurm job, run as lanes (issue M11-03c; spec
# docs/specs/M11-03c-calibration-full-sharding.md, "calibration_suite.sbatch"). Submitted by
# calibration_full_launch.sh with --job-name=gippa-cal-<suite>, --cpus-per-task, --mem and --time
# from scripts/slurm/calibration_plan.tsv; no --array.
# Environment: GIPPA_SUITE (a replicates_full key), GIPPA_WORKTREE (the checkout to run; its plan,
# durations and scripts/ci/calibration_lanes.py are used), GIPPA_CALIBRATION_LOG_DIR.
# Test seams: GIPPA_CALIBRATION_CPUSET (CPU list replacing sched_getaffinity),
# GIPPA_CALIBRATION_SYSFS (replaces /sys/devices/system/cpu), GIPPA_CALIBRATION_GRACE (seconds
# between SIGINT to the pytests and the record on USR1/TERM; default 120).
# Steps: (1) plan row and whole cores (cpus_per_lane / 2 per lane); (2) collection with
# --calibration-lane-report, lanes from `calibration_lanes.py lanes`, one --collect-only per lane
# and verify_partition; (3) every lane in the background under `taskset -c <its cores>` with
# OMP/OPENBLAS/MKL/NUMBA_NUM_THREADS and SLURM_CPUS_PER_TASK = cpus_per_lane and PYTEST_ADDOPTS =
# its deselect list + --junitxml=<log>/task-<i>.<shard<k>|group<g>>.junit.xml (shard lanes run
# _criteria.py --shard k/N; group lanes log to <log>/lanes/<suite>-group<g>); (4) after every
# lane, when shards > 0, --merge with the shard deselect list on the shard lanes' cores,
# GIPPA_CALIBRATION_CONFIRM_LANES=<shards> and GIPPA_CALIBRATION_CONFIRM_JUNIT=<log>/task-<i>.merge;
# then counts, confirm_parts.json, the exit fold and task-<i>.json (written once).
# (5) Every step runs as `( trap - INT; exec ... ) & wait`, so USR1 (10 min before the limit) and
# TERM are handled at once: SIGINT to the pytest descendants of this shell only (found through
# /proc parent links), wait for them up to the grace, record exit 143, SIGKILL what is left.
# Python runs with the venv interpreter only (no uv between the shell and pytest). Lane output:
# <log>/lanes/<suite>-<role>.out.
set -uo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)" || exit 2
export TMPDIR
STATE="$TMPDIR/state"
JOBPY="$TMPDIR/suite_job.py"
mkdir -p "$STATE"
SUITE="${GIPPA_SUITE:-}"
CHECKOUT="${GIPPA_WORKTREE:-$GIPPA_REPO}"
LOG="${GIPPA_CALIBRATION_LOG_DIR:-$GIPPA_ROOT/logs/calibration/${SLURM_JOB_ID:-local}}"
export GIPPA_CALIBRATION_LOG_DIR="$LOG"
PY="${UV_PROJECT_ENVIRONMENT:-}/bin/python"
GRACE="${GIPPA_CALIBRATION_GRACE:-120}"
TASK=""; recorded=""; in_signal=""; step_pids=()
unset GIPPA_CALIBRATION_CONFIRM GIPPA_CALIBRATION_CONFIRM_LANES GIPPA_CALIBRATION_CONFIRM_JUNIT \
  GIPPA_CALIBRATION_CONFIRM_SLICE GIPPA_CALIBRATION_REPLICATES

cat > "$JOBPY" <<'PY'
"""Helpers of calibration_suite.sbatch (M11-03c PR 3); argv[1] is the command."""
import dataclasses, datetime, importlib.util, json, os, shlex, subprocess, sys
import xml.etree.ElementTree as ET
from pathlib import Path


def lanes_mod(checkout):
    path = Path(checkout) / "scripts" / "ci" / "calibration_lanes.py"
    spec = importlib.util.spec_from_file_location("gippa_calibration_lanes", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def q(v):
    return shlex.quote(str(v))


def cpulist(text):
    out = set()
    for part in text.strip().split(","):
        part = part.strip()
        if part:
            a, _, b = part.partition("-")
            out.update(range(int(a), int(b or a) + 1))
    return out


def jload(path, default=None):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return default


def plan(state, checkout, suite, log):
    """meta.json, plan.json and plan.sh (TASK, SHARDS, NGROUPS, CPL) of the suite's row."""
    m = lanes_mod(checkout)
    try:
        rows = {r.suite: r for r in m.read_plan(Path(checkout) / "scripts/slurm/calibration_plan.tsv",
                                                m._suites())}
    except (m.PlanError, OSError, KeyError, TypeError) as exc:
        print(f"calibration_suite: {exc}", file=sys.stderr)
        return 2
    if suite not in rows:
        print(f"calibration_suite: {suite!r} is not a suite of the plan", file=sys.stderr)
        return 2
    row = rows[suite]
    meta = {"suite": suite, "task": row.task, "job": os.environ.get("SLURM_JOB_ID", ""),
            "checkout": checkout, "log": log}
    Path(state, "meta.json").write_text(json.dumps(meta))
    Path(state, "plan.json").write_text(json.dumps(dataclasses.asdict(row)))
    Path(state, "plan.sh").write_text(f"TASK={row.task}; SHARDS={row.shards}; NGROUPS={row.groups}; "
                                      f"CPL={row.cpus_per_lane}\n")
    return 0


def cores(state, cpuset, sysfs, nlanes, cpl):
    """L_CPUS[pos]: ``cpl / 2`` whole cores (both siblings inside the job's set) per lane."""
    have = cpulist(cpuset) if cpuset else set(os.sched_getaffinity(0))
    found, seen = [], set()
    for cpu in sorted(have):
        try:
            sib = cpulist(Path(sysfs, f"cpu{cpu}", "topology", "thread_siblings_list").read_text())
        except (OSError, ValueError):
            continue
        if sib <= have and min(sib) not in seen:
            seen.add(min(sib))
            found.append(sorted(sib))
    per = int(cpl) // 2
    need = int(nlanes) * per
    if len(found) < need:
        print(f"calibration_suite: {nlanes} lanes need {need} whole cores, the job's CPUs "
              f"{sorted(have)} hold {len(found)}", file=sys.stderr)
        return 2
    lines = []
    for pos in range(int(nlanes)):
        cpus = sorted(c for core in found[pos * per:(pos + 1) * per] for c in core)
        lines.append(f"L_CPUS[{pos}]={q(','.join(map(str, cpus)))}\n")
    Path(state, "cores.sh").write_text("".join(lines))
    return 0


def shlanes(state):
    """lanes.sh from lanes.json: kind, index, n, part ('' = empty lane), deselect, junit opts."""
    meta = jload(Path(state, "meta.json"))
    report = jload(Path(state, "report.json"))
    if not isinstance(report, dict) or not report:
        print("calibration_suite: the collection is empty", file=sys.stderr)
        return 2
    lanes = jload(Path(state, "lanes.json"))
    if not isinstance(lanes, list) or not lanes:
        print("calibration_suite: no lanes", file=sys.stderr)
        return 2
    out = [f"NLANES={len(lanes)}\n"]
    for pos, ln in enumerate(lanes):
        desel = " ".join(f"--deselect {q(d)}" for d in ln["deselect"])
        junit = f"--junitxml={q(os.path.join(meta['log'], ln['part']))}" if ln["part"] else ""
        out.append(f"L_KIND[{pos}]={q(ln['kind'])}; L_IDX[{pos}]={q(ln['index'])}; "
                   f"L_N[{pos}]={q(ln['n'])}; L_PART[{pos}]={q(ln['part'] or '')}; "
                   f"L_DESEL[{pos}]={q(desel)}; L_JUNIT[{pos}]={q(junit)}\n")
    Path(state, "lanes.sh").write_text("".join(out))
    return 0


def _row_lanes(state):
    meta = jload(Path(state, "meta.json"))
    m = lanes_mod(meta["checkout"])
    row = m.PlanRow(**jload(Path(state, "plan.json")))
    raw = jload(Path(state, "lanes.json"), [])
    lanes = [m.Lane(d["kind"], d["index"], d["n"], tuple(d["deselect"]), d["part"],
                    tuple(d["expected"])) for d in raw]
    return m, meta, row, lanes


def _checks(state):
    return list(jload(Path(state, "checks.json"), []))


def _add_checks(state, problems):
    Path(state, "checks.json").write_text(json.dumps(_checks(state) + list(problems)))


def collected_ids(path):
    """Nodeids of ``pytest --collect-only -q`` output: the lines before the first blank one."""
    ids = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            break
        if "::" in line:
            ids.append(line.strip())
    return ids


def partition(state):
    m, meta, row, lanes = _row_lanes(state)
    report = jload(Path(state, "report.json"), {})
    collected = {"all": sorted(report)}
    problems = []
    for pos, ln in enumerate(lanes):
        if ln.part is None:
            continue
        rc = Path(state, f"collect.{pos}.rc").read_text().strip()
        if rc != "0":
            problems.append(f"{row.suite}: --collect-only of {ln.part} exit {rc}")
        collected[ln.part] = collected_ids(Path(state, f"collect.{pos}.out"))
    problems += m.verify_partition(row, lanes, collected)
    if problems:
        _add_checks(state, problems)
        print("calibration_suite: partition problems:\n  " + "\n  ".join(problems),
              file=sys.stderr)
        return 2
    return 0


def _junit_count(path):
    """Testcases, without the nameless childless one of an interrupted lane (no phantom pass); a
    nameless one with a failure, error or skip child still counts."""
    try:
        return sum(1 for c in ET.parse(path).getroot().iter("testcase")
                   if c.get("name") or any(k.tag in ("failure", "error", "skipped") for k in c))
    except (OSError, ET.ParseError):
        return None


def _rc(path):
    try:
        return int(Path(path).read_text().strip())
    except (OSError, ValueError):
        return None


def _confirm_parts(log, suite):
    """(names, problem) of <log>/<suite>.confirm_parts.json."""
    path = Path(log, f"{suite}.confirm_parts.json")
    if not path.exists():
        return [], f"{suite}: --merge wrote no {path.name}"
    doc = jload(path)
    if not isinstance(doc, list) or not all(isinstance(p, str) and p and "/" not in p
                                            and not p.startswith(".") for p in doc):
        return [], f"{suite}: {path.name} is not a list of plain JUnit names"
    return list(doc), None


def record(state, forced):
    """Write task-<i>.json once (temporary name, renamed); print the exit it records."""
    m, meta, row, lanes = _row_lanes(state)
    log, suite = meta["log"], meta["suite"]
    launched = Path(state, "launched").exists()
    merged = Path(state, "merge.launched").exists()
    checks = _checks(state)
    merge_rc = _rc(Path(state, "exit.merge")) if merged else None
    if merged and merge_rc is None:  # started, but its exit was lost
        merge_rc = 2
    merge_parts = []
    if merged:
        merge_parts, problem = _confirm_parts(log, suite)
        if problem and forced == "-" and merge_rc in (0, 1):
            checks.append(problem)
    entries, exits, counts = [], [], {}
    for pos, ln in enumerate(lanes):
        rc = _rc(Path(state, f"exit.{pos}")) if ln.part else None
        n = _junit_count(Path(log, ln.part)) if ln.part and launched else None
        if n is not None:
            counts[ln.part] = n
        if ln.part and launched:
            exits.append(2 if rc is None else rc)
        entries.append({"kind": ln.kind, "index": ln.index, "exit": rc,
                        "tests": len(ln.expected), "junit_tests": n})
    if forced == "-":
        checks += m.verify_counts(lanes, counts)
        rc = m.fold_exit(exits, merge_rc)
        if rc in (0, 1) and checks:
            rc = 2
    else:
        rc = int(forced)
    commit, dirty = "", True
    git = ["git", "-C", meta["checkout"]]
    try:
        commit = subprocess.run([*git, "rev-parse", "HEAD"], capture_output=True, text=True,
                                check=True).stdout.strip()
        st = subprocess.run([*git, "status", "--porcelain", "--", "src", "tests/calibration",
                             "tests/sim", "tests/helpers", "ci"], capture_output=True, text=True,
                            check=True).stdout
        dirty = bool(st.strip())
    except (OSError, subprocess.CalledProcessError):
        pass
    rec = {"task": meta["task"], "suite": suite, "job": meta["job"], "commit": commit or None,
           "dirty": dirty, "exit": rc, "args": [], "sharded": row.shards > 0,
           "parts": m.record_parts(lanes, merge_parts) if launched else [],
           "date": datetime.date.today().isoformat(), "lanes": entries,
           "plan": dataclasses.asdict(row), "checks": checks, "merge_exit": merge_rc}
    out = Path(log, f"task-{meta['task']}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(f"{out}.tmp.{os.getpid()}")
    tmp.write_text(json.dumps(rec) + "\n")
    os.replace(tmp, out)
    print(rc)
    return 0


def _is_pytest(pid):
    try:
        argv = Path(f"/proc/{pid}/cmdline").read_bytes().decode(errors="replace").split("\0")
    except OSError:
        return False
    names = {os.path.basename(a) for a in argv[:2]}
    return bool(names & {"pytest", "py.test"}) or any(
        a == "-m" and b == "pytest" for a, b in zip(argv, argv[1:]))


def tree(root, which):
    """Descendants of ``root`` (``pytest`` ones only, or ``all``), this walker excluded."""
    parent = {}
    for d in os.listdir("/proc"):
        if d.isdigit():
            try:
                stat = Path(f"/proc/{d}/stat").read_text()
            except OSError:
                continue
            parent[int(d)] = int(stat.rsplit(")", 1)[1].split()[1])
    root = int(root)
    mine, p = set(), os.getpid()
    while p in parent and p != root and p > 1:
        mine.add(p)
        p = parent[p]
    kids = {}
    for c, pp in parent.items():
        kids.setdefault(pp, []).append(c)
    out, todo = [], list(kids.get(root, []))
    while todo:
        c = todo.pop()
        if c in mine:
            continue
        todo.extend(kids.get(c, []))
        if which == "all" or _is_pytest(c):
            out.append(c)
    print(" ".join(map(str, sorted(out))))
    return 0


CMDS = {"plan": plan, "cores": cores, "shlanes": shlanes, "partition": partition,
        "record": record, "tree": tree}
sys.exit(CMDS[sys.argv[1]](*sys.argv[2:]))
PY

write_record() {  # write_record <exit|->: task-<i>.json once ("-": fold); prints the exit
  [ -z "$recorded" ] && [ -n "$TASK" ] || return 0
  if "$PY" "$JOBPY" record "$STATE" "$1"; then recorded=1
  else echo "calibration_suite: could not write task-$TASK.json" >&2; fi
}
now_us() { local t=${EPOCHREALTIME//[^0-9]/}; echo $(( 10#$t )); }  # bash 5 wall clock, µs
interrupt() {  # SIGINT to this shell's pytest descendants; wait up to the grace; SIGKILL the rest
  local end p live sent=""
  end=$(( $(now_us) + GRACE * 1000000 ))  # microseconds: whole seconds could cut the grace by 1 s
  while :; do
    for p in $("$PY" -I "$JOBPY" tree "$$" pytest); do
      case " $sent " in *" $p "*) ;; *) kill -INT "$p" 2>/dev/null; sent="$sent $p" ;; esac
    done
    live=""
    for p in ${step_pids[@]+"${step_pids[@]}"} $sent; do
      if kill -0 "$p" 2>/dev/null; then live=1; fi
    done
    [ -n "$live" ] && [ "$(now_us)" -lt "$end" ] || break
    sleep 0.2
  done
  for p in $("$PY" -I "$JOBPY" tree "$$" all); do kill -KILL "$p" 2>/dev/null; done
}
on_signal() {  # on_signal USR1|TERM; TERM records first: Slurm may SIGKILL before the grace ends
  [ "${BASHPID:-$$}" = "$$" ] || exit 143
  [ -z "$in_signal" ] || return 0
  in_signal=1
  if [ "$1" = TERM ]; then write_record 143 >/dev/null; fi
  echo "calibration_suite: signal received, interrupting the lanes (grace ${GRACE} s)" >&2
  interrupt
  write_record 143 >/dev/null
  exit 143
}
on_exit() {
  local rc=$?
  [ "${BASHPID:-$$}" = "$$" ] || exit "$rc"
  write_record "$rc" >/dev/null
  rm -rf "$TMPDIR"
  exit "$rc"
}
trap on_exit EXIT
trap 'on_signal USR1' USR1
trap 'on_signal TERM' TERM

await() {  # await <pid>: its exit status, waiting again when a trapped signal ended `wait`
  local rc
  while :; do
    wait "$1"; rc=$?
    if [ "$rc" -gt 128 ] && kill -0 "$1" 2>/dev/null; then continue; fi
    return "$rc"
  done
}
thread_env() {  # in a lane's subshell: thread pools = the lane's CPUs (spec S1)
  export OMP_NUM_THREADS="$CPL" OPENBLAS_NUM_THREADS="$CPL" MKL_NUM_THREADS="$CPL" \
    NUMBA_NUM_THREADS="$CPL" SLURM_CPUS_PER_TASK="$CPL"
}
collect() {  # collect <out> <rc file> <PYTEST_ADDOPTS> [pytest args]: one --collect-only, in the background
  local out="$1" rcf="$2" opts="$3" rc pid; shift 3
  ( trap - INT; export PYTEST_ADDOPTS="$opts"
    exec "$PY" -m pytest tests/calibration -m calibration_full -k "calsuite_$SUITE" --collect-only -q "$@" >"$out" ) &
  pid=$!; step_pids+=("$pid")
  await "$pid"; rc=$?
  printf '%s\n' "$rc" > "$rcf"
  return "$rc"
}

if [ -z "$SUITE" ] || [ ! -d "$CHECKOUT" ] || [ ! -x "$PY" ]; then
  echo "calibration_suite: need GIPPA_SUITE, GIPPA_WORKTREE and \$UV_PROJECT_ENVIRONMENT/bin/python" >&2
  exit 2
fi
cd "$CHECKOUT" || exit 2
# the worktree's gippa, not the shared venv's editable install (_criteria.py checks gippa.__file__)
export PYTHONPATH="$CHECKOUT/src:$CHECKOUT/tests/sim:$CHECKOUT/tests"
mkdir -p "$LOG/lanes" || exit 2

# 1. plan row and whole cores
"$PY" "$JOBPY" plan "$STATE" "$CHECKOUT" "$SUITE" "$LOG" || exit 2
. "$STATE/plan.sh"
"$PY" "$JOBPY" cores "$STATE" "${GIPPA_CALIBRATION_CPUSET:-}" \
  "${GIPPA_CALIBRATION_SYSFS:-/sys/devices/system/cpu}" "$((SHARDS + NGROUPS))" "$CPL" || exit 2
. "$STATE/cores.sh"

# 2. collection, lanes, per-lane --collect-only, partition check
collect "$STATE/collect.all.out" "$STATE/collect.all.rc" "" \
  "--calibration-lane-report=$STATE/report.json" || { echo "calibration_suite: collection failed" >&2; exit 2; }
durations=()
if [ -f scripts/slurm/calibration_durations.tsv ]; then
  durations=(--durations scripts/slurm/calibration_durations.tsv)
fi
"$PY" scripts/ci/calibration_lanes.py lanes --suite "$SUITE" --report "$STATE/report.json" \
  --plan scripts/slurm/calibration_plan.tsv ${durations[@]+"${durations[@]}"} > "$STATE/lanes.json" || exit 2
"$PY" "$JOBPY" shlanes "$STATE" || exit 2
. "$STATE/lanes.sh"
pos=0
while [ "$pos" -lt "$NLANES" ]; do
  if [ -n "${L_PART[$pos]}" ]; then
    collect "$STATE/collect.$pos.out" "$STATE/collect.$pos.rc" "${L_DESEL[$pos]}"
  fi
  pos=$((pos + 1))
done
"$PY" "$JOBPY" partition "$STATE" || exit 2

# 3. every lane in the background on its own cores
: > "$STATE/launched"
lane_pids=(); lane_pos=(); shard_cpus=""; shard_desel=""
pos=0
while [ "$pos" -lt "$NLANES" ]; do
  if [ -n "${L_PART[$pos]}" ]; then
    role="${L_KIND[$pos]}${L_IDX[$pos]}"
    lane_log="$LOG"; extra=()
    if [ "${L_KIND[$pos]}" = shard ]; then
      extra=(--shard "${L_IDX[$pos]}/${L_N[$pos]}")
    else
      lane_log="$LOG/lanes/$SUITE-$role"; mkdir -p "$lane_log"
    fi
    ( trap - INT; thread_env
      export PYTEST_ADDOPTS="${L_DESEL[$pos]:+${L_DESEL[$pos]} }${L_JUNIT[$pos]}"
      export GIPPA_CALIBRATION_LOG_DIR="$lane_log"
      exec taskset -c "${L_CPUS[$pos]}" "$PY" tests/calibration/_criteria.py --suite "$SUITE" \
        ${extra[@]+"${extra[@]}"} >"$LOG/lanes/$SUITE-$role.out" 2>&1 ) &
    lane_pids+=("$!"); lane_pos+=("$pos"); step_pids+=("$!")
    echo "calibration_suite: $SUITE $role on CPUs ${L_CPUS[$pos]} (pid $!)"
  fi
  if [ "${L_KIND[$pos]}" = shard ]; then
    shard_cpus="${shard_cpus:+$shard_cpus,}${L_CPUS[$pos]}"; shard_desel="${L_DESEL[$pos]}"
  fi
  pos=$((pos + 1))
done

# 4. wait for every lane, then the merge (shards > 0), then counts, fold and record
j=0
while [ "$j" -lt "${#lane_pids[@]}" ]; do
  await "${lane_pids[$j]}"; rc=$?
  printf '%s\n' "$rc" > "$STATE/exit.${lane_pos[$j]}"
  echo "calibration_suite: lane ${L_KIND[${lane_pos[$j]}]}${L_IDX[${lane_pos[$j]}]} exit $rc"
  j=$((j + 1))
done
if [ "$SHARDS" -gt 0 ]; then
  rm -f "$LOG/$SUITE.confirm_parts.json"
  : > "$STATE/merge.launched"
  ( trap - INT; thread_env
    export PYTEST_ADDOPTS="$shard_desel"
    export GIPPA_CALIBRATION_CONFIRM_LANES="$SHARDS" GIPPA_CALIBRATION_CONFIRM_JUNIT="$LOG/task-$TASK.merge"
    exec taskset -c "$shard_cpus" "$PY" tests/calibration/_criteria.py --suite "$SUITE" --merge \
      >"$LOG/lanes/$SUITE-merge.out" 2>&1 ) &
  mpid=$!; step_pids+=("$mpid")
  await "$mpid"; rc=$?
  printf '%s\n' "$rc" > "$STATE/exit.merge"
  echo "calibration_suite: merge exit $rc"
fi
final=$("$PY" "$JOBPY" record "$STATE" -) && recorded=1 || final=""
case "$final" in
  ''|*[!0-9]*) echo "calibration_suite: could not write task-$TASK.json" >&2; final=2 ;;
esac
echo "calibration_suite: $SUITE (task $TASK) exit $final"
exit "$final"
