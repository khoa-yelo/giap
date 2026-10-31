#!/usr/bin/env python3
"""Release readiness table (issue M11-03; GIPPA_SPEC.md §13.2, §12 ``release.yml``, §11).

``compute`` reads the working tree of ``--repo`` (commit = HEAD), the FarmShare logs and an
optional A2 freeze run, and writes ``readiness.json`` plus ``READINESS.md`` (``render`` of the
JSON). ``verify`` re-checks the commit rule and recomputes every ``repo`` row on HEAD. Thresholds
are the spec's own literals below (R3); outputs hold aggregates only (R11). Exit codes: 0 pass,
1 fail, 3 pending (``compute``); 0/1 (``verify``). With ``--offline`` nothing touches the network.
"""

from __future__ import annotations

import argparse
import datetime as dt
import functools
import json
import math
import re
import subprocess
import sys
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

TE1_MIN = 0.995  # §13.2 "≥ 99.5 % of A2 peptides map exactly" (inclusive)
MAPPING_FDR_MAX = 0.01  # §13.2 "the decoy-estimated mapping FDR is ≤ 1 % at defaults"
RUNTIME_MIN_S = 960  # §13.2 "about 35 min (range 16-55 …)": 16 min (Q-K6, inclusive)
RUNTIME_MAX_S = 3300  # 55 min

METRICS = "validation/pilot/metrics/metrics.json"
REPORT = "validation/pilot/REPORT.md"
BASELINE = "benchmarks/baseline.json"
# a change here invalidates a calibration row (tests/sim, tests/helpers: M11-03b Q6)
# (R3/R4 reuse rows by ``reusable`` below since M11-14; the constant stays for its Q6 test)
CALIB_PATHS = ("src", "tests/calibration", "tests/sim", "tests/helpers", "ci")
RUN_PATHS = ("src", "workflow", "ci", "validation/proxies")  # ... and a weekly TE8 row
INTEG_PATHS = ("src", "workflow", "ci", "tests/integration", "tests/contract")  # ... integration
# Registry processes without an image pin (containers.config header; M8-17 decision 9; in-process)
NOT_CONTAINERISED = frozenset(
    {"uniop", "promotech", "promoter_calculator", "terminatornet", "bactermfinder", "ostir"}
    | {"fgsea", "codon_bias"}
)
PAPERWORK_FILES = frozenset({"PROGRESS.md", "README.md", "CITATION.cff", "CHANGELOG.md"})
PAPERWORK_TREES = ("validation/release/", "docs/", "changelog.d/")
# Evidence reuse (M11-14 item 4; ADR-0060 §4): calibration evidence from commit C stands at R iff
# every path of ``git diff --name-only --no-renames C R`` lies under NON_STATS_PATHS (whole path
# components; a directory entry ends with "/"). Changes to either list need `stats` review (Q-L15).
NON_STATS_PATHS: tuple[str, ...] = (
    "docs/", "changelog.d/", "README.md", "CHANGELOG.md", "CITATION.cff", "mkdocs.yml",
    "PROGRESS.md", "AGENTS.md", "CLAUDE.md", ".claude/", ".github/",
    # single files: scripts/dev/island_power.py feeds a calibration test (never list)
    "scripts/dev/approval_packet.py", "scripts/dev/progress_bars.py", "scripts/dev/replant.sh",
    "scripts/dev/waitmerge.sh", "scripts/docs/", "src/gippa/report/", "src/gippa/viz/",
    "validation/release/",
)  # fmt: skip
# Never reusable; disjoint from NON_STATS_PATHS (a unit test asserts it).
NEVER_PATHS: tuple[str, ...] = (
    "src/gippa/preflight/", "src/gippa/runner/", "src/gippa/setupdb/", "src/gippa/cli.py",
    "src/gippa/data/", "src/gippa/config/", "src/gippa/rng.py", "src/gippa/io/", "workflow/",
    "pyproject.toml", "uv.lock", "tests/sim/", "tests/calibration/", "ci/",
    "scripts/release/readiness.py", "scripts/dev/island_power.py",
)  # fmt: skip
SELF_PATH = "scripts/release/readiness.py"  # self-protection: a widening commit never certifies
VERSION_FILE = "src/gippa/__init__.py"
_VERSION_LINE = re.compile(r'__version__ = "[^"]+"')
REUSED_NOTE = "(reusable: diff ⊆ NON_STATS_PATHS)"
IMAGES = "workflow/conf/images.json"
UNPUBLISHED = frozenset({"uniop"})  # never published (R5); scripts/release/image_key.py
STEP_IMAGES = frozenset({"step", "step-r"})  # tagged <version>, not keyed (ADR-0065, Q-K39)
PROCESS_LESS = frozenset({"fgsea"})  # pinned in images.json only, no withName line
# A bare digest reference (ADR-0059); group 1 catches the docker:// prefix `docker run` rejects.
_DIGEST = re.compile(r"(docker://)?[a-z0-9][a-z0-9._/-]*@sha256:[0-9a-f]{64}")
_WITH_NAME = re.compile(
    r"withName\s*:\s*['\"]([^'\"]+)['\"]\s*\{[^}]*?\bcontainer\s*=\s*['\"]([^'\"]+)['\"]"
)
_CONTAINER = re.compile(r"\bcontainer\s*=\s*['\"]([^'\"]+)['\"]")
_PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")  # the M10-05 guard rule
# Words with a peptide-like run that are names, not peptides: file names anywhere, the process
# names of containers.config in the H4 row only; any other such run is redacted (R11).
ALLOWED_WORDS = frozenset({"CHANGELOG", "LICENSE", "READINESS"})
H4_WORDS = ALLOWED_WORDS | {"TRANSTERMHP", "ANTISMASH"}
REDACTED = "[redacted]"

CRITERIA = {
    "R1": "≥ 99.5 % of A2 peptides map exactly",
    "R2": "the decoy-estimated mapping FDR is ≤ 1 % at defaults",
    "R3": "Every statistical method passes its calibration suite (§8.4)",
    "R4": "every attack null (§8.5) stays at level",
    "R5": "Every true-example assertion that the available data allow (§11.6) passes",
    "R6": "The same seed gives byte-identical tables",
    "R7": "A standard run on 16 threads and 64 GB takes about 35 min (range 16-55 …)",
    "R8": "Core modules alone take about 25 min",
    "R9": "The pilot passes (§15.6)",
    "R10": "co-abundance … moves to the full validation",
    "H1": "milestone gates M1-M10 recorded in PROGRESS.md",
    "H2": "version agreement; release tag not yet created",
    "H3": "changelog built for the version",
    "H4": "containers pinned by digest",
    "H5": "licences listed; GIPPA is MIT",
    "H6": "CI check runs on the commit",
}
SOURCES = {"R2": "farmshare", "R3": "farmshare", "R4": "farmshare", "H6": "farmshare"}
SOURCES |= {"R5": "both", "R6": "both"}

Row = dict[str, Any]


def _row(rid: str, status: str, value: str, threshold: str = "", ev: Iterable[str] = ()) -> Row:
    return {"id": rid, "criterion": CRITERIA[rid], "status": status, "value": value,
            "threshold": threshold, "evidence": list(ev), "waiver": None,
            "source": SOURCES.get(rid, "repo")}  # fmt: skip


def overall(statuses: Iterable[str]) -> str:
    """``fail`` if any row fails, else ``pending`` if any is pending, else ``pass``."""
    s = set(statuses)
    return "fail" if "fail" in s else "pending" if "pending" in s else "pass"


def paperwork_path(path: str) -> bool:
    """Release paperwork that may change after ``commit`` (commit rule of ``verify``)."""
    if path in PAPERWORK_FILES:
        return True
    return path.startswith(PAPERWORK_TREES) and not path.startswith("docs/methods/")


def non_stats_path(path: str) -> bool:
    """``path`` lies under an entry of ``NON_STATS_PATHS`` on whole path components (pure)."""
    return bool(path) and any(
        path.startswith(e) if e.endswith("/") else path == e for e in NON_STATS_PATHS
    )


# ------------------------------------------------------------------ git and logs


def _git(repo: Path, *args: str, check: bool = True) -> str | None:
    p = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True)
    if p.returncode != 0:
        if check:
            raise RuntimeError(f"git {' '.join(args)}: {p.stderr.strip()}")
        return None
    return p.stdout.strip()


@functools.cache
def code_equal(repo: Path, commit: str, paths: tuple[str, ...]) -> bool:
    """``commit`` is HEAD or an ancestor with ``git diff --quiet commit HEAD -- paths`` (X..R)."""
    c = _git(repo, "rev-parse", "--verify", "--quiet", f"{commit}^{{commit}}", check=False) or ""
    diff = ("diff", "--quiet", "--no-renames", c, "HEAD", "--", *paths)
    tests = (("merge-base", "--is-ancestor", c, "HEAD"), diff)
    return bool(c) and all(_git(repo, *a, check=False) is not None for a in tests)


def _resolve(repo: Path, rev: str) -> str | None:
    return _git(repo, "rev-parse", "--verify", "--quiet", f"{rev}^{{commit}}", check=False) or None


@functools.cache
def _changed(repo: Path, c: str, r: str) -> tuple[str, ...] | None:
    """``git diff --name-only --no-renames c r`` (full SHAs); None if git fails."""
    out = _git(repo, "diff", "-z", "--no-ext-diff", "--name-only", "--no-renames", c, r,
               check=False)  # fmt: skip
    return None if out is None else tuple(p for p in out.split("\0") if p)


@functools.cache
def _version_line_only(repo: Path, c: str, r: str) -> bool:
    """The ``VERSION_FILE`` diff is one hunk: one ``-`` and one ``+`` ``__version__`` line."""
    out = _git(repo, "diff", "--no-ext-diff", "--no-color", "--no-textconv", "-U0", c, r, "--",
               VERSION_FILE, check=False)  # fmt: skip
    if out is None:
        return False
    lines = out.splitlines()
    start = next((i for i, ln in enumerate(lines) if ln.startswith("@@")), None)
    if start is None:
        return False
    body = lines[start:]
    if sum(ln.startswith("@@") for ln in body) != 1 or len(body) != 3:
        return False
    old, new = body[1], body[2]
    return (old[:1], new[:1]) == ("-", "+") and all(
        _VERSION_LINE.fullmatch(x[1:]) for x in (old, new)
    )


def reusable(repo: Path, c: str, r: str) -> bool:
    """Evidence from commit ``c`` stands at ``r`` (M11-14 item 4): every path of ``git diff
    --name-only --no-renames c r`` is a ``non_stats_path``, except a one-hunk ``__version__``
    change of ``src/gippa/__init__.py``; never if ``scripts/release/readiness.py`` changed or
    either commit is unknown. Reads ``NON_STATS_PATHS`` at call time."""
    cc, rr = _resolve(repo, c), _resolve(repo, r)
    if cc is None or rr is None:
        return False
    changed = _changed(repo, cc, rr)
    if changed is None or SELF_PATH in changed:
        return False
    return all(
        _version_line_only(repo, cc, rr) if p == VERSION_FILE else non_stats_path(p)
        for p in changed
    )


def reusable_ancestor(repo: Path, commit: str) -> bool:
    """``commit`` is HEAD or an ancestor of HEAD with ``reusable(commit, HEAD)`` (R3, R4)."""
    c = _resolve(repo, commit)
    if c is None or _git(repo, "merge-base", "--is-ancestor", c, "HEAD", check=False) is None:
        return False
    return reusable(repo, c, "HEAD")


@dataclass(frozen=True)
class LogRow:
    date: str
    job: str
    commit: str
    exit: int
    summary: str

    @property
    def label(self) -> str:
        return self.summary.split(":", 1)[0].strip()

    def ok(self) -> bool:
        bad = re.findall(r"(\d+) (?:failed|errors?)\b", self.summary)
        return self.exit == 0 and all(int(n) == 0 for n in bad)

    def cite(self, log: str) -> str:
        return f"{log} {self.date} {self.job} {self.commit} exit {self.exit}"


def parse_log(text: str) -> tuple[list[LogRow], int]:
    """Rows of a nightly/weekly TSV and the number of malformed lines (skipped, never a pass)."""
    rows, bad = [], 0
    for line in text.splitlines():
        f = line.split("\t", 4)
        if not line.strip():
            continue
        if len(f) < 5 or not re.fullmatch(r"-?\d+", f[3].strip()) or len(f[2].strip()) < 7:
            bad += 1
            continue
        rows.append(LogRow(f[0].strip(), f[1].strip(), f[2].strip(), int(f[3]), f[4].strip()))
    return rows, bad


def newest(
    rows: list[LogRow], commit_ok: Callable[[str], bool], pred: Callable[[LogRow], bool]
) -> LogRow | None:
    """Newest row (date, then job ID, then file order) matching ``pred`` with ``commit_ok``."""
    hits = [(r.date, [int(n) for n in re.findall(r"\d+", r.job)], i, r)
            for i, r in enumerate(rows) if pred(r) and commit_ok(r.commit)]  # fmt: skip
    return max(hits, key=lambda h: h[:3])[3] if hits else None


# ------------------------------------------------------------------ context


@dataclass
class Ctx:
    repo: Path
    version: str
    logs: Path | None = None
    freeze: Path | None = None
    offline: bool = True
    check_runs: Path | None = None
    notes: list[str] = field(default_factory=list)
    _cache: dict[str, Any] = field(default_factory=dict)

    def text(self, rel: str) -> str | None:
        p = self.repo / rel
        return p.read_text(encoding="utf-8") if p.is_file() else None

    def json(self, rel: str) -> Any:
        t = self.text(rel)
        try:
            return None if t is None else json.loads(t)
        except ValueError:
            return None

    def log(self, name: str) -> list[LogRow] | None:
        if name not in self._cache:
            p = None if self.logs is None else self.logs / name
            rows, bad = parse_log(p.read_text(encoding="utf-8")) if p and p.is_file() else ([], 0)
            self._cache[name] = rows if p and p.is_file() else None
            if bad:
                self.notes.append(f"{name}: {bad} malformed lines skipped (never read as a pass)")
        return self._cache[name]  # type: ignore[no-any-return]

    def commits(self, paths: tuple[str, ...]) -> Callable[[str], bool]:
        return lambda commit: code_equal(self.repo, commit, paths)

    def reusable_from(self, commit: str) -> bool:
        return reusable_ancestor(self.repo, commit)


def _missing(rid: str, what: str, thr: str = "") -> Row:
    return _row(rid, "pending", f"missing: {what}", thr)


# ------------------------------------------------------------------ R rows


def r1(c: Ctx) -> Row:
    thr, m = "≥ 99.5 % and te1.passed", c.json(METRICS)
    te1 = (m or {}).get("te1") or {}
    if "share_exact_inframe" not in te1:
        return _missing("R1", METRICS, thr)
    share = te1["share_exact_inframe"]
    if not isinstance(share, int | float) or isinstance(share, bool) or not math.isfinite(share):
        why = te1.get("reason") or "not a finite number"
        return _row("R1", "fail", f"share_exact_inframe {share} ({why})", thr, [METRICS])
    val = f"{share * 100:.2f} %"
    if share < TE1_MIN:
        return _row("R1", "fail", f"{val} < 99.5 %", thr, [METRICS])
    ok = te1.get("passed") is True
    return _row("R1", "pass" if ok else "fail", val + ("" if ok else "; te1.passed is not true"),
                thr, [METRICS])  # fmt: skip


def r2(c: Ctx) -> Row:
    # G1: metrics.json > mapping_fdr (an aggregate the validation role adds) is the reading. The
    # M3-10 mapping_fdr.tsv holds per-stratum placement rows (decoy_type x n_diff x length_class)
    # with no single accepted-set FDR at the default threshold, so it never yields a pass here.
    thr, m = "≤ 1 %", c.json(METRICS) or {}
    v = m.get("mapping_fdr")
    v = v.get("value") if isinstance(v, dict) else v
    if isinstance(v, int | float) and not isinstance(v, bool):
        st = "pass" if v <= MAPPING_FDR_MAX else "fail"
        return _row("R2", st, f"{v * 100:.2f} %", thr, [f"{METRICS} > mapping_fdr"])
    if c.freeze is None:
        return _row("R2", "pending", f"no mapping_fdr in {METRICS} and no --freeze-run", thr)
    tsv, rel = c.freeze / "mapping" / "mapping_fdr.tsv", "<freeze-run>/mapping/mapping_fdr.tsv"
    if not tsv.is_file():
        return _missing("R2", rel, thr)
    return _row("R2", "pending", "mapping_fdr.tsv has per-stratum rows only (gap G1)", thr, [rel])


def _labelled(label: str) -> Callable[[LogRow], bool]:
    return lambda r: r.label == label


def _calib(
    c: Ctx, rid: str, labels: tuple[str, ...], none: str, paths: tuple[str, ...] | None = None
) -> Row:
    """Newest nightly row per label for R, or for an ancestor C: code-equal on ``paths`` when
    given (R6), else ``reusable(C, R)`` (R3, R4; M11-14 item 4)."""
    thr = ("exit 0 and 0 failed for R or a code-equal ancestor" if paths is not None
           else "exit 0 and 0 failed for R or an ancestor C with reusable(C, R)")  # fmt: skip
    rows = c.log("nightly.tsv")
    if rows is None:
        return _missing(rid, "logs/nightly.tsv", thr)
    commits = c.commits(paths) if paths is not None else c.reusable_from
    hits = [newest(rows, commits, _labelled(lb)) for lb in labels]
    if any(h is None for h in hits):
        miss = [lb for lb, h in zip(labels, hits, strict=True) if h is None]
        return _row(rid, "pending", f"{none} ({', '.join(miss)})", thr)
    found = [h for h in hits if h is not None]
    bad = [h.label for h in found if not h.ok()]
    val = f"failed: {', '.join(bad)}" if bad else "; ".join(h.summary for h in found)
    if paths is None:
        head = _resolve(c.repo, "HEAD")
        olds = dict.fromkeys(x for h in found if (x := _resolve(c.repo, h.commit)) != head)
        val += "".join(f"; reused from {x[:7]} {REUSED_NOTE}" for x in olds if x)
    return _row(rid, "fail" if bad else "pass", val, thr, [h.cite("nightly.tsv") for h in found])


def r3(c: Ctx) -> Row:
    return _calib(c, "R3", ("calibration_full", "calibration_nightly"), "no calibration run for R")


def r4(c: Ctx) -> Row:
    return _calib(c, "R4", ("calibration_attacks",), "no attack row")


def _te(te: Any, status: Callable[[dict[str, Any]], Any], doc: dict[str, Any]) -> str:
    """One TE verdict: pass, pending (``doc`` without schema, or ``pending …``), deferred, fail."""
    if not isinstance(te, dict) or "schema" not in doc:
        return "pending (placeholder, no schema)"  # metrics.te_placeholders() before M10-06 data
    s = str(status(te) if te else "missing")
    reason = ", ".join(map(str, te.get("reason") or [])) or "no data"
    if s == "pass" or s.startswith("pending"):
        return s
    return {"missing": "pending", "not_evaluable": f"deferred ({reason})"}.get(s, f"fail ({s})")


def _te7(m: dict[str, Any]) -> str:
    te7 = m.get("te7") or {}
    return _te(te7, lambda t: (t.get("gippa") or {}).get("status", "missing"), te7)


def r5(c: Ctx) -> Row:
    m, parts, ev = c.json(METRICS), {}, [METRICS]
    if m is None:
        return _missing("R5", METRICS)
    parts["TE1"] = r1(c)["status"]
    real = m.get("te_real") or {}
    for k in ("te2", "te3", "te4", "te5"):
        parts[k.upper()] = _te(real.get(k) or {}, lambda t: t.get("status", "missing"), real)
    parts["TE7"] = _te7(m)
    wk = c.log("weekly.tsv")
    hit = None if wk is None else newest(wk, c.commits(RUN_PATHS), _labelled("k12_mg1655"))
    if hit is None:
        parts["TE8"] = "pending (no weekly k12_mg1655 row for R)"
    else:  # the newest k12 row decides, even when it died before the TE8 step
        ev.append(hit.cite("weekly.tsv"))
        ok = hit.ok() and "TE8" in hit.summary
        parts["TE8"] = "pass" if ok else f"fail (exit {hit.exit}; {hit.summary})"
    st = overall(s.split(" ", 1)[0] for s in parts.values())
    return _row("R5", st, "; ".join(f"{k} {s}" for k, s in parts.items()),
                "each TE pass; no data → deferred; reported → fail (Q-K5)", ev)  # fmt: skip


def r6(c: Ctx) -> Row:
    m = c.json(METRICS)
    if m is None:
        return _missing("R6", METRICS)
    calib = _calib(c, "R6", ("integration",), "no integration run for R", INTEG_PATHS)
    te7 = _te7(m)
    st = overall([calib["status"], te7.split(" ", 1)[0]])
    val = f"TE7 {te7}; integration: {calib['value']}"
    return _row("R6", st, val, "te7.gippa.status pass and the integration run",
                [METRICS, *calib["evidence"]])  # fmt: skip


def _hmm(s: float) -> str:
    mins = round(s / 60)
    return f"{mins // 60}:{mins % 60:02d}"


def _baseline_cell(c: Ctx, cell: str) -> tuple[str, dict[str, Any]] | None:
    run = (c.json(BASELINE) or {}).get("standard_run") or {}  # the spec's class (wheat) first
    for cls, cells in sorted(run.items(), key=lambda kv: kv[0] != "wheat"):
        if isinstance(cells, dict) and "median_s" in (cells.get(cell) or {}):
            return cls, cells[cell]
    return None


def r7(c: Ctx) -> Row:
    thr, hit = (
        f"{_hmm(RUNTIME_MIN_S)}-{_hmm(RUNTIME_MAX_S)} (inclusive)",
        _baseline_cell(c, "k12_standard"),
    )
    if hit is None:
        return _missing("R7", f"{BASELINE} > standard_run.*.k12_standard", thr)
    cls, cell = hit
    s, ev = float(cell["median_s"]), [f"{BASELINE} > standard_run.{cls}.k12_standard"]
    where = f"node class {cls} ({cell.get('threads', '?')} threads)"
    out = (f" > {_hmm(RUNTIME_MAX_S)}" if s > RUNTIME_MAX_S
           else f" < {_hmm(RUNTIME_MIN_S)}" if s < RUNTIME_MIN_S else "")  # fmt: skip
    return _row("R7", "fail" if out else "pass", f"{_hmm(s)}{out} ({s:.0f} s); {where}", thr, ev)


def r8(c: Ctx) -> Row:
    hit = _baseline_cell(c, "k12_core")
    if hit is None:
        return _missing("R8", f"{BASELINE} > standard_run.*.k12_core", "reported only")
    s = float(hit[1]["median_s"])
    return _row("R8", "pass", f"{_hmm(s)} ({s:.0f} s); node class {hit[0]}", "reported only",
                [f"{BASELINE} > standard_run.{hit[0]}.k12_core"])  # fmt: skip


def _gate(progress: str, k: int) -> str | None:
    m = re.search(rf"(?ms)^### Gate M{k}(?![\w-])(.*?)(?=^#{{1,3}} |\Z)", progress)
    return None if m is None else m.group(1)


def r9(c: Ctx) -> Row:
    thr, prog, rep = (
        "Gate M10 with a Khoa review line; no fail/pending in REPORT.md",
        c.text("PROGRESS.md"),
        c.text(REPORT),
    )
    if rep is None:
        return _missing("R9", REPORT, thr)
    block = _gate(prog or "", 10)
    if block is None or not re.search(r"(?im)^\s*[-*]?\s*Khoa review\b", block):
        return _row("R9", "pending", "no Khoa review line under ### Gate M10", thr, ["PROGRESS.md"])
    lines, sts, col = rep.splitlines(), [], None
    tab = [[x.strip().lower() for x in ln.strip().strip("|").split("|")]
           if ln.lstrip().startswith("|") else None for ln in lines]  # fmt: skip
    for i, cells in enumerate(tab):  # the first table with a "Status…" header cell, and only it
        if (
            col := next((j for j, x in enumerate(cells or []) if x.startswith("status")), None)
        ) is not None:
            for body in tab[i + 2 :]:
                if body is None:
                    break
                sts += body[col : col + 1]
            break
    if col is None or not sts:
        why = "no status table" if col is None else "empty status table"
        return _row("R9", "pending", f"{why} in {REPORT}", thr, [REPORT])
    bad = [s for s in sts if s.startswith(("fail", "pending"))]
    st = overall("fail" if s.startswith("fail") else "pending" for s in bad)
    return _row(
        "R9", st, f"{len(sts)} bullets; not passed: {len(bad)}", thr, ["PROGRESS.md", REPORT]
    )


def r10(c: Ctx) -> Row:
    return _row("R10", "deferred", "spec: co-abundance moves to the full validation", "none")


# ------------------------------------------------------------------ H rows


def h1(c: Ctx) -> Row:
    prog = c.text("PROGRESS.md")
    if prog is None:
        return _missing("H1", "PROGRESS.md", "### Gate M1 … M10")
    miss = [f"M{k}" for k in range(1, 11) if _gate(prog, k) is None]
    val = f"missing: {', '.join(miss)}" if miss else "10 gate blocks"
    return _row("H1", "fail" if miss else "pass", val, "### Gate M1 … M10", ["PROGRESS.md"])


def h2(c: Ctx) -> Row:
    init = re.search(
        r"(?m)^__version__\s*=\s*['\"]([^'\"]+)", c.text("src/gippa/__init__.py") or ""
    )
    cff = re.search(r"(?m)^version:\s*['\"]?([^'\"\s]+)", c.text("CITATION.cff") or "")
    got = {"gippa.__version__": init and init.group(1), "CITATION.cff": cff and cff.group(1)}
    bad = [f"{k} = {v}" for k, v in got.items() if v != c.version]
    if _git(c.repo, "tag", "-l", f"v{c.version}"):
        bad.append(f"tag v{c.version} exists")
    val = "; ".join(bad) if bad else f"{c.version} everywhere; no tag v{c.version}"
    return _row("H2", "fail" if bad else "pass", val, f"= {c.version}",
                ["src/gippa/__init__.py", "CITATION.cff", "git tag -l"])  # fmt: skip


def h3(c: Ctx) -> Row:
    head = f"## GIPPA {c.version} ("
    d = c.repo / "changelog.d"
    frags = (
        sorted(p.name for p in d.rglob("*") if p.is_file() and p.name != ".gitkeep")
        if d.is_dir()
        else []
    )
    bad = ([] if head in (c.text("CHANGELOG.md") or "") else [f"no '{head}' section"]) + (
        [f"unbuilt fragments: {', '.join(frags)}"] if frags else []
    )
    val = "; ".join(bad) or "built"
    return _row("H3", "fail" if bad else "pass", val, "section present, no fragments",
                ["CHANGELOG.md", "changelog.d/"])  # fmt: skip


def registry(text: str) -> dict[str, bool]:
    """Tool names of ``tools.py``'s ``TOOLS`` registry -> in-process flag (parsed, not imported)."""
    seg = text.split("TOOLS: Mapping", 1)[-1].split("PILOT_TOOLS", 1)[0]
    names = dict.fromkeys(re.findall(r'_tool\(\s*"([\w-]+)"', seg), False)
    for blk in re.findall(r"ToolSpec\((.*?)\n\s*\)", seg, re.S):
        if m := re.search(r'name="([\w-]+)"', blk):
            names[m.group(1)] = "in_process=True" in blk
    return names


_H4_EXPECTED = "<image>@sha256:<64 hex>, no docker:// prefix"


def _h4_reason(name: str, img: str) -> str | None:
    """Why a ``container =`` value fails H4 (None if it is a bare digest reference)."""
    m = _DIGEST.fullmatch(img)
    if m is None:
        return f"tag, not digest: {name}"
    return f"docker:// prefix: {name}" if m.group(1) else None


def h4(c: Ctx) -> Row:
    cfg_rel, tools_rel = "workflow/conf/containers.config", "src/gippa/annotate/tools.py"
    cfg, reg = c.text(cfg_rel), registry(c.text(tools_rel) or "")
    if cfg is None or not reg:
        return _missing("H4", cfg_rel if cfg is None else f"{tools_rel} registry")
    text = "\n".join(ln for ln in cfg.splitlines() if not ln.lstrip().startswith("//"))
    pins = dict(_WITH_NAME.findall(text))
    bad = [r for p, img in pins.items() if (r := _h4_reason(p, img))]
    bad += [
        r
        for img in dict.fromkeys(_CONTAINER.findall(text))
        if img not in pins.values() and (r := _h4_reason(img, img))
    ]  # every `container =` (§12)
    bad += [f"not pinned, not conda-only: {n}" for n, inproc in reg.items()
            if not inproc and n not in NOT_CONTAINERISED and n.upper() not in pins]  # fmt: skip
    bad += _h4_images(c.repo)
    unpinned = sum(n.upper() not in pins for n in reg)
    val = "; ".join(bad) or f"{len(pins)} digests; {unpinned} conda-only/in-process"
    return _row("H4", "fail" if bad else "pass", val, _H4_EXPECTED, [cfg_rel, tools_rel, IMAGES])


def _image_key(tool: str, root: Path) -> str:
    """``gippa.io.image_key.key`` of the checkout next to this script (ValueError: no inputs)."""
    src = str(Path(__file__).resolve().parents[2] / "src")
    if src not in sys.path:
        sys.path.insert(0, src)  # the checkout's image_key, not an installed gippa's
    from gippa.io.image_key import key

    return key(tool, root)


def _h4_images(root: Path) -> list[str]:
    """H4 reasons from ``images.json`` (M11-14 item 2): every GHCR tool (a directory under
    ``workflow/containers/`` but ``uniop`` and the step images) has a current entry; every entry
    has image inputs and processes (fgsea excepted). Never raises."""
    base = root / "workflow" / "containers"
    skip = UNPUBLISHED | STEP_IMAGES
    tools = sorted(p.name for p in base.iterdir() if p.is_dir() and p.name not in skip) \
        if base.is_dir() else []  # fmt: skip
    path = root / IMAGES
    if not path.is_file():
        return [f"not in {IMAGES} (absent): {t}" for t in tools]
    try:
        images = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        return [f"{IMAGES} unreadable ({type(e).__name__})"]
    if not isinstance(images, dict):
        return [f"{IMAGES} is not an object"]
    bad = [f"not in {IMAGES}: {t}" for t in tools if t not in images]
    for tool, entry in sorted(images.items()):
        if not isinstance(entry, dict) or not isinstance(entry.get("key"), str):
            bad.append(f"invalid {IMAGES} entry: {tool}")
            continue
        try:
            k = _image_key(str(tool), root)
        except (OSError, ValueError):
            bad.append(f"no image inputs: {tool}")
            continue
        if k != entry["key"]:
            bad.append(f"stale image: {tool} (inputs changed since {entry['key']})")
        procs = entry.get("processes")
        if not isinstance(procs, list) or (not procs and tool not in PROCESS_LESS):
            bad.append(f"no processes in {IMAGES}: {tool}")
    return bad


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", s.lower())


def h5(c: Ctx) -> Row:
    notice, lic = c.text("NOTICE.md") or "", (c.text("LICENSE") or "").strip()
    listed = {_norm(x) for x in re.findall(r"(?m)^\|\s*([^|]+?)\s*\|", notice)}
    bad = [f"not in NOTICE.md: {n}" for n in registry(c.text("src/gippa/annotate/tools.py") or "")
           if _norm(n) not in listed]  # fmt: skip
    if not lic.startswith("MIT License"):
        bad.append("LICENSE is not MIT")
    return _row("H5", "fail" if bad else "pass", "; ".join(bad) or "all tools listed; MIT",
                "registry ⊆ NOTICE.md table; MIT", ["NOTICE.md", "LICENSE"])  # fmt: skip


H6_PASS = frozenset({"success", "skipped"})  # skipped: PR-only jobs (guards, summary) on main


def h6(c: Ctx) -> Row:
    # "all success" read as success or skipped (PR-only jobs never run on a main commit); any other
    # conclusion fails; a null conclusion (still running) is pending; at least one success.
    thr, sha = "every check run success or skipped", _git(c.repo, "rev-parse", "HEAD") or ""
    if c.check_runs is not None:
        doc = json.loads(c.check_runs.read_text(encoding="utf-8"))
        got = doc.get("check_runs", [])
        runs = [(str(r.get("name")), r.get("conclusion") or "none") for r in got]
        ev = [c.check_runs.name]
        if other := sorted({str(r.get("head_sha")) for r in got} - {sha}):
            return _row("H6", "fail", f"check runs for {', '.join(other)}, not {sha}", thr, ev)
    elif c.offline:
        return _row("H6", "pending", "offline: GitHub check runs not read", thr)
    else:
        q = '.check_runs[] | [.name, .conclusion // "none"] | @tsv'
        out = _git(c.repo, "rev-parse", "--show-toplevel")  # gh resolves {owner}/{repo} here
        api = f"repos/{{owner}}/{{repo}}/commits/{sha}/check-runs"
        try:
            p = subprocess.run(["gh", "api", "--paginate", api, "--jq", q],
                               capture_output=True, text=True, cwd=out)  # fmt: skip
        except FileNotFoundError:
            return _row("H6", "pending", "gh not on PATH", thr)
        if p.returncode != 0:
            return _row("H6", "pending", f"gh api failed: {p.stderr.strip()[:200]}", thr)
        runs = [(n, s) for n, _, s in (ln.partition("\t") for ln in p.stdout.splitlines()) if s]
        ev = [f"gh api …/commits/{sha}/check-runs"]
    bad = [f"{n}: {s}" for n, s in runs if s not in H6_PASS and s != "none"]
    running = [f"{n}: running" for n, s in runs if s == "none"]
    if bad or running or not any(s == "success" for _, s in runs):
        st = "fail" if bad else "pending"
        return _row("H6", st, "; ".join(bad + running) or "no successful check run", thr, ev)
    skip = [n for n, s in runs if s == "skipped"]
    val = f"{len(runs) - len(skip)} check runs success"
    return _row("H6", "pass", val + (f"; skipped: {', '.join(skip)}" if skip else ""), thr, ev)


ROWS: tuple[Callable[[Ctx], Row], ...] = (
    r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, h1, h2, h3, h4, h5, h6,
)  # fmt: skip


def apply_waivers(repo: Path, rows: list[Row]) -> list[Row]:
    """An accepted ADR "Waives readiness <id>" turns that row's fail/pending into ``waived``."""
    adrs = sorted((repo / "docs" / "adr").glob("ADR-*.md")) if (repo / "docs/adr").is_dir() else []
    for r in rows:
        for p in adrs:
            text = p.read_text(encoding="utf-8")
            if r["status"] not in ("fail", "pending") or not re.search(
                rf"Waives readiness {r['id']}\b", text
            ):
                continue
            num = re.match(r"ADR-\d{4}", p.name)
            status = re.search(r"(?im)^\W*Status:\s*(\S+)", text)  # the header line only
            if status and status.group(1).lower() == "accepted-by-khoa" and num:
                r["evidence"].append(
                    f"docs/adr/{p.name} (accepted-by-khoa; waived from {r['status']})"
                )
                r.update(status="waived", waiver=num.group())
            else:
                r["evidence"].append(f"docs/adr/{p.name}: not accepted")
    return rows


# ------------------------------------------------------------------ render, compute, verify


def redact(text: str, allow: frozenset[str] = ALLOWED_WORDS) -> str:
    """Replace every peptide-like run (M10-05 rule) outside an ``allow`` word, verbatim
    otherwise: no case change, so the guard still sees a leak that slipped through (R11)."""

    def word(m: re.Match[str]) -> str:
        return m.group() if m.group() in allow else _PEPTIDE_LIKE.sub(REDACTED, m.group())

    return re.sub(r"[A-Za-z]{7,}", word, text)


def _redact_doc(doc: Any, allow: frozenset[str] = ALLOWED_WORDS) -> Any:
    if isinstance(doc, str):
        return redact(doc, allow)
    if isinstance(doc, list):
        return [_redact_doc(x, allow) for x in doc]
    if isinstance(doc, dict):
        allow = H4_WORDS if doc.get("id") == "H4" else allow
        return {k: _redact_doc(v, allow) for k, v in doc.items()}
    return doc


def _cell(s: object) -> str:
    return str(s).replace("|", "\\|").replace("\n", " ")


def render(doc: dict[str, Any]) -> str:
    """``READINESS.md``: a pure function of ``readiness.json``; text verbatim but ``redact``-ed,
    line by line (process names only in the H4 row); no row-level data is ever read (R11)."""
    out = [f"# GIPPA {doc['version']} release check", "", f"- Commit: `{doc['commit']}`",
           f"- Computed: {doc['computed']}", f"- Overall: **{doc['overall']}**", "",
           "| id | criterion | status | value | threshold | evidence | source |",
           "|---|---|---|---|---|---|---|"]  # fmt: skip
    out = [redact(ln) for ln in out]
    for r in doc["rows"]:
        ev = [*r["evidence"], *([f"waiver {r['waiver']}"] if r.get("waiver") else [])]
        cells = (r["id"], r["criterion"], r["status"], r["value"], r["threshold"],
                 "<br>".join(ev) or "none", r["source"])  # fmt: skip
        line = "| " + " | ".join(_cell(x) for x in cells) + " |"
        out.append(redact(line, H4_WORDS if r["id"] == "H4" else ALLOWED_WORDS))
    out += ["", "## Notes", "", *[redact(f"- {n}") for n in doc.get("notes") or []],
            f"- Peptide guard (M10-05): runs of 7+ amino-acid letters are shown as {REDACTED}, "
            "except names on the fixed allow-list.", ""]  # fmt: skip
    return "\n".join(out)


def compute(ctx: Ctx) -> dict[str, Any]:
    rows = apply_waivers(ctx.repo, [f(ctx) for f in ROWS])
    doc = {"version": ctx.version, "commit": _git(ctx.repo, "rev-parse", "HEAD"),
           "computed": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
           "overall": overall(r["status"] for r in rows), "rows": rows,
           "notes": ctx.notes}  # fmt: skip
    clean = _redact_doc(doc)
    if clean != doc:
        clean["notes"].append(f"peptide-like text {REDACTED} in readiness.json (M10-05 guard)")
    return clean  # type: ignore[no-any-return]


def commit_rule(repo: Path, commit: str) -> str:
    if (
        not re.fullmatch(r"[0-9a-f]{40}", commit)
        or _git(repo, "cat-file", "-e", f"{commit}^{{commit}}", check=False) is None
    ):
        return f"fail (unknown commit {commit})"
    if _git(repo, "rev-parse", "HEAD") == commit:
        return "pass"
    if _git(repo, "merge-base", "--is-ancestor", commit, "HEAD", check=False) is None:
        return f"fail (commit {commit} is not an ancestor of HEAD)"
    stale = [
        p
        for p in (
            _git(repo, "diff", "--name-only", "--no-renames", commit, "HEAD") or ""
        ).splitlines()
        if not paperwork_path(p)
    ]
    return f"fail (stale: {stale[0]})" if stale else "pass"


def verify(path: Path, repo: Path) -> int:
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
        rule, ctx = commit_rule(repo, str(doc.get("commit", ""))), Ctx(repo, str(doc["version"]))
    except (OSError, ValueError, KeyError, AttributeError) as e:
        print(f"overall: fail (unreadable {path.name}: {e!r})")
        return 1
    fresh = {
        r["id"]: r
        for r in apply_waivers(
            repo, [f(ctx) for f in ROWS if SOURCES.get(f.__name__.upper(), "repo") == "repo"]
        )
    }
    rows = [fresh.get(r["id"], r) for r in doc["rows"]]
    print(f"commit rule: {rule}")
    for r in rows:
        print(redact(f"{r['id']}\t{r['status']}\t{r['value']}"))
    ids_ok = [r["id"] for r in rows] == [f.__name__.upper() for f in ROWS]
    statuses = [r["status"] for r in rows] + ["pass" if rule == "pass" and ids_ok else "fail"]
    print(f"overall: {overall(statuses)}")
    return 0 if overall(statuses) == "pass" else 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="readiness.py", description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compute")
    c.add_argument("--version", required=True)
    c.add_argument("--repo", type=Path, default=Path("."))
    c.add_argument("--logs", type=Path)
    c.add_argument("--freeze-run", type=Path)
    c.add_argument("--out", type=Path, required=True)
    c.add_argument("--offline", action="store_true")
    c.add_argument("--check-runs", type=Path, help="saved gh api …/check-runs response (H6)")
    v = sub.add_parser("verify")
    v.add_argument("readiness", type=Path)
    v.add_argument("--repo", type=Path, default=Path("."))
    a = ap.parse_args(argv)
    if a.cmd == "verify":
        return verify(a.readiness, a.repo)
    doc = compute(Ctx(a.repo, a.version, a.logs, a.freeze_run, a.offline, a.check_runs))
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "readiness.json").write_text(
        json.dumps(doc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (a.out / "READINESS.md").write_text(render(doc), encoding="utf-8")
    print(f"overall: {doc['overall']} ({a.out})")
    return {"pass": 0, "fail": 1}.get(doc["overall"], 3)


if __name__ == "__main__":
    sys.exit(main())
