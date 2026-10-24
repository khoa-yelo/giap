"""Check the pilot runbook ``validation/pilot/reproduce.md`` (M10-14, stdlib only).

``static RUNBOOK [--repo DIR]`` checks the runbook against a checkout: scripts, job classes,
committed input hashes and their records, stray hashes, the A3 branch and, for a final runbook,
the Evidence table (docs/specs/M10-14-pilot-reproducibility.md, criterion 1). ``verify RUNBOOK
--out J [--repo DIR]`` hashes the Inputs rows that exist on FarmShare (``$GIAP_ROOT``) and checks
the records it can read. Exit 0 no failure, 1 failures, 2 usage. The only warning is
``duplicate_hash:<id>,<id>``. ``protected RUNBOOK`` prints the freeze commit, then the expanded
paths of every Inputs row and every Outputs row with ``kept = yes`` (for the clean-up job).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import stat
import sys
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

HEX64 = re.compile(r"(?<![0-9a-fA-F])[0-9a-f]{64}(?![0-9a-fA-F])")
HEADER = re.compile(r"<!-- reproduce: (.*) -->")
STEP = re.compile(r"^### (S\d+|S-[a-z]+) (.*?) \(([^)]*)\)\s*$", re.M)
SECTION = re.compile(r"^#{2,3} ", re.M)
FENCE = re.compile(r"^```bash\n(.*?)^```", re.M | re.S)
A3 = re.compile(r"^### A3-(cascadia|fallback) \(status: (\w+)\)", re.M)
PATH_TOKEN = re.compile(
    r"(?:\$\{?GIAP_REPO\}?/)?((?:scripts|validation|workflow|tests/data)/[^\s\"'`;|)(<>]+)"
)
SBATCH_LINE = re.compile(r"^\s*(\w+=\S*\s+)*sbatch\b")
SRUN_LINE = re.compile(r"^\s*(\w+=\S*\s+)*srun\b")
FORBIDDEN = re.compile(r"\buv\s+(run|sync)\b|\bpip3?\s|\bpython3? -m pip\b")
JOB_FIELDS = (
    *("script", "name", "partition", "qos", "cpus", "gpus", "mem", "time"),
    *("constraint", "array"),
)
# Job field -> (#SBATCH long option, short option)
CLASS_KEYS = {
    "partition": ("partition", "p"),
    "qos": ("qos", None),
    "cpus": ("cpus-per-task", "c"),
    "mem": ("mem", None),
    "time": ("time", "t"),
    "constraint": ("constraint", "C"),
    "array": ("array", "a"),
}
ROOTS = ("$GIAP_ROOT", "$GIAP_DATA_DIR", "$GIAP_DB_DIR")
NOT_HASHED = (".tdf_bin",)  # gigabytes; M10-02 decision, spec L11
READ_ONLY = (".sif",)  # the DIA-NN copy must be mode 0444
FINAL_COMMIT_STEPS = {f"S{n}" for n in range(8, 13)}


@dataclass
class Step:
    id: str
    title: str
    issue: str
    job: dict[str, str] | None
    commands: list[str] = field(default_factory=list)
    reads: list[str] = field(default_factory=list)
    writes: list[str] = field(default_factory=list)
    job_lines: list[str] = field(default_factory=list)


@dataclass
class Runbook:
    status: str
    freeze_commit: str
    a3_branch: str
    inputs: list[dict[str, str]]
    evidence: list[dict[str, str]]
    outputs: list[dict[str, str]]
    a3: dict[str, str]
    steps: list[Step]
    text: str = ""


def _section(text: str, heading: str) -> str:
    m = re.search(rf"^## {re.escape(heading)}\s*$(.*?)(?=^## |\Z)", text, re.M | re.S)
    return m.group(1) if m else ""


def _table(body: str) -> list[dict[str, str]]:
    lines = [ln.strip() for ln in body.splitlines() if ln.strip().startswith("|")]
    if len(lines) < 2:
        return []

    def cells(ln: str) -> list[str]:
        return [c.strip() for c in ln.strip().strip("|").split("|")]

    cols = cells(lines[0])
    return [dict(zip(cols, cells(ln), strict=False)) for ln in lines[2:]]


def _list(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip() not in ("", "-", "none")]


def _parse_step(m: re.Match[str], body: str) -> Step:
    job_lines = re.findall(r"^Job:\s*(.*)$", body, re.M)
    job: dict[str, str] | None = None
    if job_lines and not job_lines[0].startswith("none"):
        job = dict(kv.split("=", 1) for kv in job_lines[0].split() if "=" in kv)
    reads = re.search(r"^Reads:\s*(.*)$", body, re.M)
    writes = re.search(r"^Writes:\s*(.*)$", body, re.M)
    return Step(
        id=m.group(1),
        title=m.group(2),
        issue=m.group(3),
        job=job,
        commands=FENCE.findall(body),
        reads=_list(reads.group(1)) if reads else [],
        writes=_list(writes.group(1)) if writes else [],
        job_lines=job_lines,
    )


def parse(text: str) -> Runbook:
    first = text.splitlines()[0] if text else ""
    hm = HEADER.fullmatch(first.strip())
    header = dict(kv.split("=", 1) for kv in hm.group(1).split() if "=" in kv) if hm else {}
    steps = []
    for m in STEP.finditer(text):
        nxt = SECTION.search(text, m.end())
        steps.append(_parse_step(m, text[m.end() : nxt.start() if nxt else len(text)]))
    return Runbook(
        status=header.get("status", ""),
        freeze_commit=header.get("freeze_commit", ""),
        a3_branch=header.get("a3_branch", ""),
        inputs=_table(_section(text, "Inputs")),
        evidence=_table(_section(text, "Evidence")),
        outputs=_table(_section(text, "Outputs")),
        a3=dict(A3.findall(text)),
        steps=steps,
        text=text,
    )


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sbatch_header(path: Path) -> dict[str, str]:
    pat = re.compile(r"^#SBATCH\s+--([\w-]+)(?:[= ](\S+))?", re.M)
    return {m.group(1): m.group(2) or "" for m in pat.finditer(path.read_text())}


def _logical(block: str) -> list[str]:
    return block.replace("\\\n", " ").splitlines()


def _option(cmds: str, long: str, short: str | None) -> str | None:
    pat = rf"(?<!\S)--{long}[= ]\"?([^\s\"]+)"
    if short:
        pat += rf"|(?<!\S)-{short}\s+\"?([^\s\"]+)"
    m = re.search(pat, cmds)
    return next(g for g in m.groups() if g is not None) if m else None


def _gpus(hdr: dict[str, str]) -> str:
    if "gpus" in hdr:
        return hdr["gpus"].rsplit(":", 1)[-1]
    gres = hdr.get("gres", "")
    return gres.rsplit(":", 1)[-1] if gres.startswith("gpu") else "0"


def _check_job(step: Step, repo: Path) -> list[str]:
    job = step.job
    if job is None:
        return []
    out = [f"job_field_missing:{step.id}:{f}" for f in JOB_FIELDS if f not in job]
    name = job.get("name", "")
    if not name.startswith("giap-"):
        out.append(f"job_name:{step.id}:{name}")
    script = job.get("script", "-")
    if script == "-":  # a dev job without a script (srun): nothing to compare
        return out
    if not (repo / script).is_file():
        return [*out, f"missing_script:{script} ({step.id})"]
    hdr = sbatch_header(repo / script)
    lines = [ln for c in step.commands for ln in _logical(c)]
    cmds = " ".join(ln for ln in lines if SBATCH_LINE.match(ln) or SRUN_LINE.match(ln))
    expected = {f: _option(cmds, *CLASS_KEYS[f]) or hdr.get(CLASS_KEYS[f][0]) for f in CLASS_KEYS}
    expected["constraint"] = expected["constraint"] or "-"
    expected["array"] = expected["array"] or "-"
    expected["gpus"] = _option(cmds, "gpus", "G") or _gpus(hdr)
    for f, want in expected.items():
        if f in job and job[f] != want:
            out.append(f"job_class_mismatch:{step.id}:{f}:{job[f]}!={want}")
    jnames = re.findall(r"(?:-J|--job-name)[= ]\s*\"?([^\s\"]+)", cmds)
    if name and (jnames or hdr.get("job-name")) and name not in (jnames or [hdr["job-name"]]):
        out.append(f"job_class_mismatch:{step.id}:name:{name}!={'|'.join(jnames)}")
    return out


def _check_blocks(rb: Runbook, repo: Path) -> list[str]:
    out: list[str] = []
    owners: list[tuple[str, str]] = []
    for s in rb.steps:
        owners += [(s.id, c) for c in s.commands] + [(s.id, j) for j in s.job_lines]
    owned = {c for _, c in owners}
    owners += [("-", c) for c in FENCE.findall(rb.text) if c not in owned]
    for sid, chunk in owners:
        for tok in PATH_TOKEN.findall(chunk):
            tok = tok.rstrip(".,:")
            if not any(ch in tok for ch in "*?[{$") and not (repo / tok).exists():
                out.append(f"missing_path:{tok} ({sid})")
        if chunk.startswith(("script=", "none")):  # a Job line: no commands
            continue
        for ln in _logical(chunk):
            if FORBIDDEN.search(ln):
                out.append(f"forbidden_command:{sid}:{ln.strip()}")
            if not SBATCH_LINE.match(ln):
                continue
            scripts = re.findall(r"(scripts/slurm/[\w.-]+\.sbatch)", ln)
            if not scripts:
                out.append(f"sbatch_without_script:{sid}:{ln.strip()}")
            for sc in scripts:
                if not (repo / sc).is_file():
                    out.append(f"missing_script:{sc} ({sid})")
            for jn in re.findall(r"(?:-J|--job-name)[= ]\s*\"?([^\s\"]+)", ln):
                if not jn.startswith("giap-"):
                    out.append(f"job_name:{sid}:{jn}")
    return out


def _check_rows(rows: list[dict[str, str]], repo: Path, kind: str) -> list[str]:
    out: list[str] = []
    seen: dict[str, str] = {}
    by_hash: dict[str, tuple[str, str]] = {}
    for r in rows:
        rid, path, sha = r.get("id", ""), r.get("path", ""), r.get("sha256", "")
        if not re.fullmatch(r"[a-z0-9_]+", rid):
            out.append(f"bad_id:{rid}")
        if rid in seen:
            out.append(f"duplicate_id:{rid}")
        seen[rid] = path
        if not re.fullmatch(r"[0-9a-f]{64}", sha):
            out.append(f"bad_sha256:{rid}")
        if sha in by_hash and by_hash[sha][1] != path:
            out.append(f"duplicate_hash:{by_hash[sha][0]},{rid}")
        by_hash.setdefault(sha, (rid, path))
        committed = not path.startswith(ROOTS)
        if committed and (path.startswith("/") or not (repo / path).is_file()):
            out.append(f"bad_path:{rid}:{path}")
            committed = False
        elif committed and sha256_file(repo / path) != sha:
            out.append(f"hash_mismatch:{rid}")
        recs = [x.strip() for x in r.get("recorded_in", "").split(";") if x.strip()]
        if kind == "outputs":
            recs = [x for x in recs if x != "-"]
            kept = r.get("kept", "")
            if not re.fullmatch(r"yes|superseded|deleted:\d{4}-\d{2}-\d{2}", kept):
                out.append(f"bad_kept:{rid}:{kept}")
        elif not recs:
            out.append(f"no_record:{rid}")
        if kind == "inputs" and committed and recs and all(x == path for x in recs):
            out.append(f"self_record:{rid}")
        for rec in recs:
            if rec.startswith("$"):
                continue  # FarmShare record: only verify reads it
            p = repo / rec
            if not p.is_file():
                out.append(f"missing_record:{rid}:{rec}")
            elif sha not in p.read_text(errors="replace"):
                out.append(f"not_recorded:{rid}:{rec}")
    return out


def _check_a3(rb: Runbook, repo: Path) -> list[str]:
    out = [f"a3_missing:{b}" for b in ("cascadia", "fallback") if b not in rb.a3]
    out += [f"a3_status:{b}:{s}" for b, s in rb.a3.items() if s not in ("active", "inactive")]
    active = sorted(b for b, s in rb.a3.items() if s == "active")
    want = [] if rb.a3_branch == "pending" else [rb.a3_branch]
    if active != want:
        out.append(f"a3_header_mismatch:header={rb.a3_branch},active={','.join(active) or '-'}")
    if rb.a3_branch == "cascadia":  # the fallback may also follow an E1-E3 outcome
        p = repo / "validation" / "pilot" / "metrics" / "a3_feasibility.json"
        verdict = json.loads(p.read_text()).get("verdict") if p.is_file() else "none"
        if verdict != "pass":
            out.append(f"a3_branch_mismatch:verdict={verdict},active=cascadia")
    return out


def _check_final(rb: Runbook) -> list[str]:
    if rb.status != "final":
        return []
    out = []
    if not re.fullmatch(r"[0-9a-f]{40}", rb.freeze_commit):
        out.append(f"final_freeze_commit:{rb.freeze_commit}")
    for r in rb.evidence:
        sid = r.get("step", "")
        for col, val in r.items():
            if val.strip() in ("", "pending", "TBD"):
                out.append(f"evidence_incomplete:{sid}:{col}")
        late = sid in FINAL_COMMIT_STEPS and r.get("state") != "superseded"
        if late and r.get("commit") != rb.freeze_commit:
            out.append(f"evidence_commit:{sid}:{r.get('commit')}")
    return out


def static_check(runbook: Runbook, repo: Path) -> list[str]:
    """Failure strings for ``runbook`` against the checkout ``repo`` (``duplicate_hash`` warns)."""
    rb = runbook
    out: list[str] = []
    if rb.status not in ("draft", "final"):
        out.append(f"bad_header:status={rb.status}")
    if not (rb.freeze_commit == "pending" or re.fullmatch(r"[0-9a-f]{40}", rb.freeze_commit)):
        out.append(f"bad_header:freeze_commit={rb.freeze_commit}")
    if rb.a3_branch not in ("cascadia", "fallback", "pending"):
        out.append(f"bad_header:a3_branch={rb.a3_branch}")
    if not rb.steps:
        out.append("no_steps")
    if not rb.inputs:
        out.append("no_inputs")
    ids = [s.id for s in rb.steps]
    out += [f"duplicate_step:{i}" for i in sorted({i for i in ids if ids.count(i) > 1})]
    for s in rb.steps:
        if len(s.job_lines) != 1:
            out.append(f"job_lines:{s.id}:{len(s.job_lines)}")
        out += _check_job(s, repo)
    out += _check_blocks(rb, repo)
    out += _check_rows(rb.inputs, repo, "inputs")
    out += _check_rows(rb.outputs, repo, "outputs")
    known = {r.get("sha256") for r in rb.inputs + rb.outputs}
    out += [f"stray_hash:{h}" for h in sorted(set(HEX64.findall(rb.text)) - known)]
    out += _check_a3(rb, repo)
    out += _check_final(rb)
    return list(dict.fromkeys(out))


def _resolve(path: str, giap_root: Path, env: Mapping[str, str], repo: Path) -> Path:
    bases = {
        "$GIAP_ROOT": str(giap_root),
        "$GIAP_DATA_DIR": env.get("GIAP_DATA_DIR") or str(giap_root / "data"),
        "$GIAP_DB_DIR": env.get("GIAP_DB_DIR") or str(giap_root / "db"),
    }
    for var, base in bases.items():
        if path == var or path.startswith(var + "/"):
            return Path(base + path[len(var) :])
    return repo / path


def verify(
    runbook: Runbook, giap_root: Path, repo: Path, env: Mapping[str, str] | None = None
) -> dict[str, Any]:
    """Hash every Inputs row that exists and check each record it can read (S-check step 2)."""
    env = {} if env is None else env
    failures: list[str] = []
    rows: dict[str, dict[str, Any]] = {}
    for r in runbook.inputs:
        rid, sha = r["id"], r["sha256"]
        p = _resolve(r["path"], giap_root, env, repo)
        row: dict[str, Any] = {"path": str(p)}
        rows[rid] = row
        if r["path"].endswith(NOT_HASHED):
            row["status"] = "not_checkable"
            continue
        if not p.is_file():
            row["status"] = "absent"
            continue
        got = sha256_file(p)
        row["sha256"] = got
        row["status"] = "ok" if got == sha else "mismatch"
        if got != sha:
            failures.append(f"mismatch:{rid}")
        if r["path"].endswith(READ_ONLY):
            mode = stat.S_IMODE(p.stat().st_mode)
            row["mode"] = f"{mode:o}"
            if mode != 0o444:
                row["status"] = "mismatch"
                failures.append(f"mode:{rid}:{mode:o}")
        for rec in (x.strip() for x in r.get("recorded_in", "").split(";") if x.strip()):
            rp = _resolve(rec, giap_root, env, repo)
            if not rp.is_file():
                failures.append(f"record_absent:{rid}:{rec}")
                row["status"] = "mismatch"
            elif sha not in rp.read_text(errors="replace"):
                failures.append(f"not_recorded:{rid}:{rec}")
                row["status"] = "mismatch"
    return {"freeze_commit": runbook.freeze_commit, "failures": failures, "inputs": rows}


def protected(runbook: Runbook, giap_root: Path, env: Mapping[str, str], repo: Path) -> list[str]:
    rows = runbook.inputs + [r for r in runbook.outputs if r.get("kept") == "yes"]
    return [str(_resolve(r["path"], giap_root, env, repo)).rstrip("/") for r in rows]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    default_repo = Path(__file__).resolve().parents[2]
    for name in ("static", "verify", "protected"):
        p = sub.add_parser(name)
        p.add_argument("runbook", type=Path)
        p.add_argument("--repo", type=Path, default=default_repo)
        if name == "verify":
            p.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    rb = parse(a.runbook.read_text())
    if a.cmd == "static":
        issues = static_check(rb, a.repo)
        for s in issues:
            print(("warning: " if s.startswith("duplicate_hash:") else "") + s)
        return 1 if any(not s.startswith("duplicate_hash:") for s in issues) else 0
    root = os.environ.get("GIAP_ROOT", "")
    if not root:
        print("reproduce_check: GIAP_ROOT is not set", file=sys.stderr)
        return 2
    if a.cmd == "protected":
        print(rb.freeze_commit)
        print("\n".join(protected(rb, Path(root), os.environ, a.repo)))
        return 0
    res = verify(rb, Path(root), a.repo, os.environ)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=2, sort_keys=True) + "\n")
    for f in res["failures"]:
        print(f, file=sys.stderr)
    return 1 if res["failures"] else 0


if __name__ == "__main__":
    sys.exit(main())
