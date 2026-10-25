#!/usr/bin/env python3
"""Release readiness table (issue M11-03; GIAP_SPEC.md §13.2, §12 ``release.yml``, §11).

``compute`` reads the working tree of ``--repo`` (commit = HEAD), the FarmShare logs and an
optional A2 freeze run, and writes ``readiness.json`` plus ``READINESS.md`` (``render`` of the
JSON). ``verify`` re-checks the commit rule and recomputes every ``repo`` row on HEAD. Thresholds
are the spec's own literals below (R3); outputs hold aggregates only (R11). Exit codes: 0 pass,
1 fail, 3 pending (``compute``); 0/1 (``verify``). With ``--offline`` nothing touches the network.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
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
CALIB_PATHS = ("src", "tests/calibration", "ci")  # a change here invalidates a calibration row
RUN_PATHS = ("src", "workflow", "ci", "validation/proxies")  # ... and a weekly TE8 row
# Registry processes without an image pin (containers.config header; M8-17 decision 9; in-process)
NOT_CONTAINERISED = frozenset(
    {"uniop", "promotech", "promoter_calculator", "terminatornet", "bactermfinder", "ostir"}
    | {"fgsea", "codon_bias"}
)
PAPERWORK_FILES = frozenset({"PROGRESS.md", "README.md", "CITATION.cff", "CHANGELOG.md"})
PAPERWORK_TREES = ("validation/release/", "docs/", "changelog.d/")
_DIGEST = re.compile(r"docker://[^@\s]+@sha256:[0-9a-f]{64}")
_WITH_NAME = re.compile(
    r"withName\s*:\s*['\"]([^'\"]+)['\"]\s*\{[^}]*?\bcontainer\s*=\s*['\"]([^'\"]+)['\"]"
)
_PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")  # the M10-05 guard rule

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
    "H5": "licences listed; GIAP is MIT",
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


# ------------------------------------------------------------------ git and logs


def _git(repo: Path, *args: str, check: bool = True) -> str | None:
    p = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True)
    if p.returncode != 0:
        if check:
            raise RuntimeError(f"git {' '.join(args)}: {p.stderr.strip()}")
        return None
    return p.stdout.strip()


def eligible(repo: Path, paths: tuple[str, ...]) -> list[str]:
    """HEAD and its ancestors with no change under ``paths`` since them (newest first)."""
    last = _git(repo, "rev-list", "-n1", "HEAD", "--", *paths)
    if not last:
        return (_git(repo, "rev-list", "HEAD") or "").split()
    return [*(_git(repo, "rev-list", f"{last}..HEAD") or "").split(), last]


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


def newest(rows: list[LogRow], commits: list[str], pred: Callable[[LogRow], bool]) -> LogRow | None:
    """Newest row (date, then job ID, then file order) for one of ``commits`` matching ``pred``."""
    hits = [(r.date, [int(n) for n in re.findall(r"\d+", r.job)], i, r)
            for i, r in enumerate(rows)
            if pred(r) and any(c.startswith(r.commit) for c in commits)]  # fmt: skip
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

    def commits(self, paths: tuple[str, ...]) -> list[str]:
        return self._cache.setdefault(str(paths), eligible(self.repo, paths))  # type: ignore[no-any-return]


def _missing(rid: str, what: str, thr: str = "") -> Row:
    return _row(rid, "pending", f"missing: {what}", thr)


# ------------------------------------------------------------------ R rows


def r1(c: Ctx) -> Row:
    thr, m = "≥ 99.5 % and te1.passed", c.json(METRICS)
    te1 = (m or {}).get("te1") or {}
    if "share_exact_inframe" not in te1:
        return _missing("R1", METRICS, thr)
    share = float(te1["share_exact_inframe"])
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
    if isinstance(v, int | float):
        st = "pass" if v <= MAPPING_FDR_MAX else "fail"
        return _row("R2", st, f"{v * 100:.2f} %", thr, [f"{METRICS} > mapping_fdr"])
    if c.freeze is None:
        return _row("R2", "pending", f"no mapping_fdr in {METRICS} and no --freeze-run", thr)
    tsv = c.freeze / "mapping" / "mapping_fdr.tsv"
    if not tsv.is_file():
        return _missing("R2", str(tsv), thr)
    return _row(
        "R2", "pending", "mapping_fdr.tsv has per-stratum rows only (gap G1)", thr, [str(tsv)]
    )


def _labelled(label: str) -> Callable[[LogRow], bool]:
    return lambda r: r.label == label


def _calib(c: Ctx, rid: str, labels: tuple[str, ...], none: str) -> Row:
    thr, rows = "exit 0 and 0 failed for R or a code-equal ancestor", c.log("nightly.tsv")
    if rows is None:
        return _missing(rid, "logs/nightly.tsv", thr)
    commits = c.commits(CALIB_PATHS)
    hits = [newest(rows, commits, _labelled(lb)) for lb in labels]
    if any(h is None for h in hits):
        miss = [lb for lb, h in zip(labels, hits, strict=True) if h is None]
        return _row(rid, "pending", f"{none} ({', '.join(miss)})", thr)
    found = [h for h in hits if h is not None]
    bad = [h.label for h in found if not h.ok()]
    val = f"failed: {', '.join(bad)}" if bad else "; ".join(h.summary for h in found)
    return _row(rid, "fail" if bad else "pass", val, thr, [h.cite("nightly.tsv") for h in found])


def r3(c: Ctx) -> Row:
    return _calib(c, "R3", ("calibration_full", "calibration_nightly"), "no calibration run for R")


def r4(c: Ctx) -> Row:
    return _calib(c, "R4", ("calibration_attacks",), "no attack row")


def r5(c: Ctx) -> Row:
    m, parts, ev = c.json(METRICS), {}, [METRICS]
    if m is None:
        return _missing("R5", METRICS)
    parts["TE1"] = r1(c)["status"]
    for k in ("te2", "te3", "te4", "te5"):
        te = (m.get("te_real") or {}).get(k) or {}
        s = str(te.get("status", "missing"))
        reason = ", ".join(map(str, te.get("reason") or [])) or "no data"
        known = {"pass": "pass", "missing": "pending", "not_evaluable": f"deferred ({reason})"}
        parts[k.upper()] = known.get(s, f"fail ({s})")  # reported (not passed) is fail (Q-K5)
    parts["TE7"] = "pass" if (m.get("te7") or {}).get("passed") is True else "fail"
    wk = c.log("weekly.tsv")
    hit = None if wk is None else newest(wk, c.commits(RUN_PATHS), lambda r: "TE8" in r.summary)
    if hit is None:
        parts["TE8"] = "pending (no weekly TE8 row for R)"
    else:
        parts["TE8"] = "pass" if hit.exit == 0 else f"fail (weekly exit {hit.exit})"
        ev.append(hit.cite("weekly.tsv"))
    st = overall(s.split(" ", 1)[0] for s in parts.values())
    return _row("R5", st, "; ".join(f"{k} {s}" for k, s in parts.items()),
                "each TE pass; no data → deferred; reported → fail (Q-K5)", ev)  # fmt: skip


def r6(c: Ctx) -> Row:
    m = c.json(METRICS)
    if m is None:
        return _missing("R6", METRICS)
    calib = _calib(c, "R6", ("integration",), "no integration run for R")
    te7 = (m.get("te7") or {}).get("passed") is True
    st = overall([calib["status"], "pass" if te7 else "fail"])
    val = f"te7.passed {str(te7).lower()}; integration: {calib['value']}"
    return _row("R6", st, val, "te7.passed and the integration run", [METRICS, *calib["evidence"]])


def _hmm(s: float) -> str:
    mins = round(s / 60)
    return f"{mins // 60}:{mins % 60:02d}"


def _baseline_cell(c: Ctx, cell: str) -> tuple[str, dict[str, Any]] | None:
    for cls, cells in ((c.json(BASELINE) or {}).get("standard_run") or {}).items():
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
    if s > RUNTIME_MAX_S:
        return _row(
            "R7", "fail", f"{_hmm(s)} > {_hmm(RUNTIME_MAX_S)} ({s:.0f} s); {where}", thr, ev
        )
    if s < RUNTIME_MIN_S:
        return _row(
            "R7", "fail", f"{_hmm(s)} < {_hmm(RUNTIME_MIN_S)} ({s:.0f} s); {where}", thr, ev
        )
    return _row("R7", "pass", f"{_hmm(s)} ({s:.0f} s); {where}", thr, ev)


def r8(c: Ctx) -> Row:
    hit = _baseline_cell(c, "k12_core")
    if hit is None:
        return _missing("R8", f"{BASELINE} > standard_run.*.k12_core", "reported only")
    s = float(hit[1]["median_s"])
    return _row("R8", "pass", f"{_hmm(s)} ({s:.0f} s); node class {hit[0]}", "reported only",
                [f"{BASELINE} > standard_run.{hit[0]}.k12_core"])  # fmt: skip


def _gate(progress: str, k: int) -> str | None:
    m = re.search(rf"(?ms)^### Gate M{k}\b(.*?)(?=^#{{1,3}} |\Z)", progress)
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
    lines = [ln for ln in rep.splitlines() if ln.lstrip().startswith("|")]
    head = next((i for i, ln in enumerate(lines) if "status" in ln.lower()), None)
    if head is None:
        return _row("R9", "pending", f"no status table in {REPORT}", thr, [REPORT])
    col = [x.strip().lower() for x in lines[head].strip("| ").split("|")].index("status")
    cells = [ln.strip("| ").split("|") for ln in lines[head + 2 :]]
    sts = [x[col].strip().lower() for x in cells if len(x) > col]
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
    init = re.search(r"(?m)^__version__\s*=\s*['\"]([^'\"]+)", c.text("src/giap/__init__.py") or "")
    cff = re.search(r"(?m)^version:\s*['\"]?([^'\"\s]+)", c.text("CITATION.cff") or "")
    got = {"giap.__version__": init and init.group(1), "CITATION.cff": cff and cff.group(1)}
    bad = [f"{k} = {v}" for k, v in got.items() if v != c.version]
    if _git(c.repo, "tag", "-l", f"v{c.version}"):
        bad.append(f"tag v{c.version} exists")
    val = "; ".join(bad) if bad else f"{c.version} everywhere; no tag v{c.version}"
    return _row("H2", "fail" if bad else "pass", val, f"= {c.version}",
                ["src/giap/__init__.py", "CITATION.cff", "git tag -l"])  # fmt: skip


def h3(c: Ctx) -> Row:
    head = f"## GIAP {c.version} ("
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


def h4(c: Ctx) -> Row:
    cfg_rel, tools_rel = "workflow/conf/containers.config", "src/giap/annotate/tools.py"
    cfg, reg = c.text(cfg_rel), registry(c.text(tools_rel) or "")
    if cfg is None or not reg:
        return _missing("H4", cfg_rel if cfg is None else f"{tools_rel} registry")
    text = "\n".join(ln for ln in cfg.splitlines() if not ln.lstrip().startswith("//"))
    pins = dict(_WITH_NAME.findall(text))
    bad = [f"tag, not digest: {p}" for p, img in pins.items() if not _DIGEST.fullmatch(img)]
    bad += [f"not pinned, not conda-only: {n}" for n, inproc in reg.items()
            if not inproc and n not in NOT_CONTAINERISED and n.upper() not in pins]  # fmt: skip
    val = "; ".join(bad) or f"{len(pins)} digests; {len(reg) - len(pins)} conda-only/in-process"
    return _row(
        "H4", "fail" if bad else "pass", val, "docker://…@sha256:<64 hex>", [cfg_rel, tools_rel]
    )


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", s.lower())


def h5(c: Ctx) -> Row:
    notice, lic = c.text("NOTICE.md") or "", (c.text("LICENSE") or "").strip()
    listed = {_norm(x) for x in re.findall(r"(?m)^\|\s*([^|]+?)\s*\|", notice)}
    bad = [f"not in NOTICE.md: {n}" for n in registry(c.text("src/giap/annotate/tools.py") or "")
           if _norm(n) not in listed]  # fmt: skip
    if not lic.startswith("MIT License"):
        bad.append("LICENSE is not MIT")
    return _row("H5", "fail" if bad else "pass", "; ".join(bad) or "all tools listed; MIT",
                "registry ⊆ NOTICE.md table; MIT", ["NOTICE.md", "LICENSE"])  # fmt: skip


def h6(c: Ctx) -> Row:
    thr, sha = "every check run success", _git(c.repo, "rev-parse", "HEAD") or ""
    if c.check_runs is not None:
        doc = json.loads(c.check_runs.read_text(encoding="utf-8"))
        runs = [(r["name"], r.get("conclusion")) for r in doc.get("check_runs", [])]
        ev = [str(c.check_runs)]
    elif c.offline:
        return _row("H6", "pending", "offline: GitHub check runs not read", thr)
    else:
        q = '.check_runs[] | [.name, .conclusion // "none"] | @tsv'
        out = _git(c.repo, "rev-parse", "--show-toplevel")  # gh resolves {owner}/{repo} here
        api = f"repos/{{owner}}/{{repo}}/commits/{sha}/check-runs"
        p = subprocess.run(["gh", "api", "--paginate", api, "--jq", q],
                           capture_output=True, text=True, cwd=out)  # fmt: skip
        if p.returncode != 0:
            return _row("H6", "pending", f"gh api failed: {p.stderr.strip()[:200]}", thr)
        runs = [tuple(ln.split("\t", 1)) for ln in p.stdout.splitlines() if "\t" in ln]  # type: ignore[misc]
        ev = [f"gh api …/commits/{sha}/check-runs"]
    if not runs:
        return _row("H6", "pending", "no check runs", thr, ev)
    bad = [f"{n}: {s}" for n, s in runs if s != "success"]
    val = "; ".join(bad) or f"{len(runs)} check runs success"
    return _row("H6", "fail" if bad else "pass", val, thr, ev)


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
            if re.search(r"(?im)^\W*Status:\s*accepted-by-khoa\b", text) and num:
                r.update(status="waived", waiver=num.group())
                r["evidence"].append(f"docs/adr/{p.name} (accepted-by-khoa)")
            else:
                r["evidence"].append(f"docs/adr/{p.name}: not accepted")
    return rows


# ------------------------------------------------------------------ render, compute, verify


def _cell(s: object) -> str:
    return str(s).replace("|", "\\|").replace("\n", " ")


def render(doc: dict[str, Any]) -> str:
    """``READINESS.md``: a pure function of ``readiness.json``."""
    out = [f"# GIAP {doc['version']} release check", "", f"- Commit: `{doc['commit']}`",
           f"- Computed: {doc['computed']}", f"- Overall: **{doc['overall']}**", "",
           "| id | criterion | status | value | threshold | evidence | source |",
           "|---|---|---|---|---|---|---|"]  # fmt: skip
    for r in doc["rows"]:
        ev = [*r["evidence"], *([f"waiver {r['waiver']}"] if r.get("waiver") else [])]
        cells = (r["id"], r["criterion"], r["status"], r["value"], r["threshold"],
                 "<br>".join(ev) or "none", r["source"])  # fmt: skip
        out.append("| " + " | ".join(_cell(x) for x in cells) + " |")
    out += ["", "## Notes", "", *[f"- {n}" for n in doc.get("notes") or ["none"]], ""]
    # File and tool names such as the changelog or licence file are lowercased here so that no
    # uppercase run can pass for a peptide; no row-level data is ever read (R11).
    return _PEPTIDE_LIKE.sub(lambda m: m.group().lower(), "\n".join(out))


def compute(ctx: Ctx) -> dict[str, Any]:
    rows = apply_waivers(ctx.repo, [f(ctx) for f in ROWS])
    return {"version": ctx.version, "commit": _git(ctx.repo, "rev-parse", "HEAD"),
            "computed": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
            "overall": overall(r["status"] for r in rows), "rows": rows,
            "notes": ctx.notes}  # fmt: skip


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
        for p in (_git(repo, "diff", "--name-only", commit, "HEAD") or "").splitlines()
        if not paperwork_path(p)
    ]
    return f"fail (stale: {stale[0]})" if stale else "pass"


def verify(path: Path, repo: Path) -> int:
    doc = json.loads(path.read_text(encoding="utf-8"))
    rule = commit_rule(repo, str(doc.get("commit", "")))
    ctx = Ctx(repo, str(doc["version"]))
    fresh = {
        r["id"]: r
        for r in apply_waivers(
            repo, [f(ctx) for f in ROWS if SOURCES.get(f.__name__.upper(), "repo") == "repo"]
        )
    }
    rows = [fresh.get(r["id"], r) for r in doc["rows"]]
    print(f"commit rule: {rule}")
    for r in rows:
        print(f"{r['id']}\t{r['status']}\t{r['value']}")
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
