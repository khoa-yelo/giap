"""Builders for the M10-15a tests (the pilot gate tooling; docs/specs/M10-15-pilot-gate.md).

Everything is synthetic and hand-made in ``tmp``: the runbook (``reproduce.md``, M10-14's
format: the header line, ``## Evidence`` and the ``## A3`` subsections), ``baseline.json``
(``standard_run``), an ADR file, ``a3_feasibility.json``, a small final-mode template, a minimal
``metrics.json`` and the FarmShare side of a review packet (per-arm ``report.html``,
``compare.html``, ``verify-<F7>.json``, the benchmark ``summary.md``). Nothing comes from
``$GIPPA_ROOT``, ``$GIPPA_DB_DIR`` or a DIA-NN report; no peptide sequence is written.

Interface pinned here (the spec gives the command lines; the rest is fixed by these tests; the
"Lead decisions after the 15a code review" of the spec are binding):

* ``metrics.py compute`` takes ``--baseline``, ``--benchmark-adr`` and ``--a3-feasibility`` but
  **not** ``--runbook`` (B1: argparse exit 2). ``metrics.json: evidence`` holds ``benchmark``,
  ``a3_feasibility`` and ``exit_codes`` (``{"A1": int|null, "A2": int|null, "A3": int|null}``,
  the arms' ``run_manifest.json: exit_code``, null for an arm not given) and never the runbook
  (no ``runbook``/``steps`` content, no ``inputs.runbook``). Without ``--a3-feasibility`` the
  block is ``{"verdict": null, "reason": "pending (M10-08 verdict)"}`` (S7). Checklist row 0 of
  ``compute`` keeps M10-05's plain strings (``pass`` / ``partial (A1, A2; A3 pending M10-10)`` /
  ``fail``); the evidence strings of ``compute`` hold no job ID.
* ``metrics.py render --final --metrics M --runbook R [--notes N] --out O [--template T]``:
  ``--runbook`` is required in final mode (absent: exit 2). Render recomputes bullet 1 and the job
  IDs of the evidence strings (rows 1-3: S8/S9/S11, S14, S12) from the runbook's Evidence rows,
  ignoring ``state == superseded`` rows (B2), and fills ``evidence.runbook``/``evidence.steps``
  for the template. Bullet 1 (S2): an arm exit ``n != 0`` gives ``fail (<arm> exit <n>)``; the
  label comes only from ``a3_branch``: ``pending``/unknown gives ``partial (...)``; no A3 run
  gives ``fail (A3 not run: <feasibility reason>)``; ``cascadia`` with feasibility verdict
  ``fail`` gives ``fail (A3 feasibility fail)``; otherwise ``pass (A3 = Cascadia)`` or
  ``pass (A3 = Casanovo fallback on UTI89, different strain)``.
  Exit 1 for: a runbook not ``status=final``; a ``pending``/``partial`` checklist status; a
  ``pending`` or absent TE status (``te_real.te<n>.status``, ``te7.gippa.status``), whatever the
  ``schema`` key (``final_problems``, its message names the path); an empty ``evidence`` string
  (message names ``evidence``); an absent placeholder path; a ``**missing:`` string (the message
  names what is missing); a non-pass TE without a ``## Diagnosis: TE<n>`` section naming one
  cause class. Exit 2 for a ``metrics.json`` without ``checklist`` or a final render without
  ``--runbook``. Nothing written on exit 1 or 2. A placeholder under a metric whose ``reason`` is
  a finding (``different_genome``, ``different_strain``, ``empty_discoveries``, ``a3_not_run``)
  renders ``not evaluable (<reason>)``; in draft mode a null under a ``pending`` reason renders
  that reason (N4), never a bare ``NA``.
* checklist rows are dicts ``{"bullet", "status", "evidence"}``.
* A benchmark cell is **complete** when ``standard_run.wheat.<cell>`` exists and its ``median_s``
  is not null (M10-12 writes entries for complete cells only and gives them no ``status`` key);
  the fixtures mark an incomplete one with ``median_s: null`` and ``status: "incomplete"``.
  ADR-M10-D ``rejected``/``superseded`` gives ``fail (ADR-M10-D <status>)`` (N1).
* ``review_packet.py --metrics-dir D --runbook R --round K --out OUT --repo REPO --a1 DIR --a2 DIR
  [--a3 DIR] --cmp DIR --verify FILE --bench-summary FILE --adr FILE [--adr FILE ...]``. The
  committed files (``validation/pilot/REPORT.md``, ``validation/pilot/metrics/metrics.json``,
  ``validation/pilot/report_notes.md``, ``validation/pilot/metrics/a3_feasibility.json``,
  ``validation/pilot/report_template.md``, ``PROGRESS.md``) are read under ``--repo``; the
  per-arm HTML is ``<DIR>/report/report.html``, the compare HTML ``<CMP>/report/compare.html``.
  The round is written to ``OUT/r<K>/``. Check 1 re-renders with ``--runbook R`` and compares
  bytes (N7), and its message carries the refusal reason (S5); the freeze directory's
  ``metrics.json`` must equal the committed one byte for byte. A3's HTML is required when the
  re-rendered bullet 1 reports an A3 run (``pass (A3 = ...)`` or ``fail (A3 exit <n>)``).
  ``OUT`` must lie outside ``--repo`` or be git-ignored by it (S3, exit 2). An existing round is
  refused with a message starting ``6:`` (exit 1); an error while copying exits 2, with no
  traceback and no ``OUT/.r<K>.tmp`` left (S4). The script sets ``sys.dont_write_bytecode``
  before importing ``validation.pilot`` (N6).
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from helpers.pilot_metrics_m10 import GIPPA_COMMIT, REPO

REVIEW_PACKET_PY = REPO / "validation" / "pilot" / "review_packet.py"
F = GIPPA_COMMIT
F7 = F[:7]
OTHER_COMMIT = "fedcba9876543210fedcba9876543210fedcba98"
#: the source separator of the report (``metrics.SOURCE_SEP``)
SEP = " \u203a "
CELLS = ("k12_standard", "k12_core", "sim5m_standard", "sim5m_core")
CAUSE_CLASSES = ("pipeline_bug", "annotation", "search", "biology", "definition")

#: Evidence rows of a final runbook: step -> (label, job)
EVIDENCE_JOBS: dict[str, tuple[str, str]] = {
    "S8": ("A2", "9100008"),
    "S9": ("A1", "9100009"),
    "S10": ("A2-rep", "9100010"),
    "S11": ("A3", "9100011"),
    "S12": ("compare", "9100012"),
    "S13": ("evidence", "9100013"),
    "S14": ("metrics", "9100014"),
    "S-bench": ("bench", "9100099"),
}

BULLETS = (
    "GIPPA runs end to end on A1, A2 and A3",
    ">= 99.5 % of A2 peptides map exactly (TE1)",
    "A1 and A2 detections and abundances are compared",
    "Measured runtimes replace the estimates",
    "TE2-TE5 are reported; failures explained",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_text(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def write_json(path: Path, doc: Any) -> Path:
    return write_text(path, json.dumps(doc, sort_keys=True, indent=2) + "\n")


# ------------------------------------------------------------------ committed inputs


def runbook_text(
    *,
    status: str = "final",
    freeze_commit: str = F,
    a3_branch: str = "cascadia",
    steps: Sequence[str] | None = None,
    commits: Mapping[str, str] | None = None,
    superseded: Sequence[tuple[str, str, str]] = (),
    superseded_first: bool = False,
) -> str:
    """A runbook in M10-14's format with an Evidence row per step of ``steps`` (default: every
    step of ``EVIDENCE_JOBS``, without S11 when ``a3_branch`` is not ``cascadia``/``fallback``).
    ``commits`` overrides the ``commit`` cell of a step (default ``freeze_commit``).
    ``superseded``: extra ``(step, job, commit)`` rows with ``state`` ``superseded`` and exit 1
    (M10-14 §4), written after the accepted rows (before them with ``superseded_first``)."""
    if steps is None:
        steps = [s for s in EVIDENCE_JOBS if s != "S11" or a3_branch in ("cascadia", "fallback")]
    commits = dict(commits or {})
    head = (
        f"<!-- reproduce: status={status} freeze_commit={freeze_commit} a3_branch={a3_branch} -->"
    )
    rows = [
        "| step | label | job | state | exit | elapsed | maxrss | node | commit |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for s in steps:
        label, job = EVIDENCE_JOBS[s]
        c = commits.get(s, freeze_commit)
        rows.append(f"| {s} | {label} | {job} | COMPLETED | 0 | 01:02:03 | 4G | wheat-05 | {c} |")
    old = [
        f"| {s} | {EVIDENCE_JOBS[s][0]} | {job} | superseded | 1 | 00:10:00 | 1G | wheat-01 | {c} |"
        for s, job, c in superseded
    ]
    rows = rows[:2] + old + rows[2:] if superseded_first else rows + old
    active = {"cascadia": ("active", "inactive"), "fallback": ("inactive", "active")}.get(
        a3_branch, ("inactive", "inactive")
    )
    return "\n".join(
        [
            head,
            "# Pilot reproduction runbook (fixture)",
            "",
            "## Evidence",
            "",
            *rows,
            "",
            "## A3",
            "",
            f"### A3-cascadia (status: {active[0]})",
            "",
            f"### A3-fallback (status: {active[1]})",
            "",
        ]
    )


def baseline_entry(median_s: float | None = 3912.0) -> dict[str, Any]:
    """A ``standard_run.wheat.<cell>`` entry in M10-12's shape (``median_s`` null: incomplete)."""
    entry: dict[str, Any] = {
        "median_s": median_s,
        "runs_s": [3700.0, 3912.0, 4106.0],
        "spread": 0.1038,
        "noisy": False,
        "nodes": ["wheat-05"] * 3,
        "cpu": "Intel Xeon Gold 6148 (CPU_SKU:6148)",
        "threads": 16,
        "mem_mb": 64000,
        "task_memory": "56.GB",
        "profile": "farmshare",
        "modules": "core,mobile,bgc,function",
        "bakta": "full",
        "steps": {},
        "estimate_s": 2100,
        "ratio_to_estimate": 1.863,
        "within_estimate_range": False,
        "slurm_jobs": ["9100099_0"],
        "commit": OTHER_COMMIT,
        "date": "2026-10-25",
        "rounds": 3,
        "inputs_sha256": "a" * 64,
        "adr": "ADR-0099",
    }
    if median_s is None:
        entry.update(runs_s=[3700.0, 3912.0], status="incomplete")
    return entry


def write_baseline(path: Path, cells: Mapping[str, str | None] | None = None) -> Path:
    """``baseline.json`` with an unrelated micro-benchmark key and ``standard_run.wheat``:
    ``cells`` maps a cell to ``complete``/``incomplete``/``None`` (absent); default all
    complete. ``cells == {}`` writes no ``standard_run`` key at all."""
    doc: dict[str, Any] = {"annotation_steps": {"wheat": {"median_s": 12.5}}}
    states = {c: "complete" for c in CELLS} if cells is None else dict(cells)
    if states:
        doc["standard_run"] = {
            "wheat": {
                c: baseline_entry(None if st == "incomplete" else 3912.0)
                for c, st in states.items()
                if st is not None
            }
        }
    return write_json(path, doc)


def write_adr(
    path: Path, status: str, title: str = "Benchmark definition and first baseline"
) -> Path:
    return write_text(
        path,
        f"# {path.stem}: {title}\n- Status: {status}\n- Date: 2026-10-25\n"
        "- Needs Khoa: yes\n\n## Context\nFixture.\n",
    )


def write_feasibility(path: Path, verdict: str, reason: str | None = None) -> Path:
    return write_json(
        path, {"schema": "gippa-a3-feasibility/1", "verdict": verdict, "reason": reason}
    )


# ------------------------------------------------------------------ final-mode render fixtures

FIXTURE_TEMPLATE = """# Pilot report (fixture template)

## 1. Status

| bullet | status | evidence |
|---|---|---|
| runs | {{checklist.0.status}} | {{checklist.0.evidence}} |
| TE1 | {{checklist.1.status}} | {{checklist.1.evidence}} |
| compare | {{checklist.2.status}} | {{checklist.2.evidence}} |
| runtimes | {{checklist.3.status}} | {{checklist.3.evidence}} |
| TE2-TE5 | {{checklist.4.status}} | {{checklist.4.evidence}} |

GIPPA commit {{gippa_commit}}; runbook status {{evidence.runbook.status}}.
Submitted for Khoa's review (§14.8 item 4).

## 3. Metrics

Agreement: {{m4.spearman_rho}}.

## 5. TE2-TE5 and TE7

| TE | status |
|---|---|
| TE2 | {{te_real.te2.status}} |
| TE3 | {{te_real.te3.status}} |
| TE4 | {{te_real.te4.status}} |
| TE5 | {{te_real.te5.status}} |
| TE7 | {{te7.gippa.status}} |

## 7. Not evaluated

<!-- not evaluated -->

## 10. Notes

<!-- notes -->
"""


#: checklist row 0 as ``compute`` writes it (M10-05's plain strings; render --final relabels it)
DRAFT_BULLET1 = {0: "pass", None: "partial (A1, A2; A3 pending M10-10)"}


def checklist_rows(
    statuses: Sequence[str] | None = None, a3_exit: int | None = 0
) -> list[dict[str, str]]:
    """The checklist as ``compute`` writes it: no job ID in any evidence string (B1)."""
    row0 = DRAFT_BULLET1.get(a3_exit, "fail")
    statuses = list(statuses or [row0, "pass", "pass", "pass", "pass"])
    ev = [
        f"metrics.json{SEP}evidence.exit_codes; reproduce.md{SEP}Evidence",
        f"metrics.json{SEP}te1.share_exact_inframe, te1.passed",
        f"metrics.json{SEP}m4, m5",
        f"baseline.json{SEP}standard_run.wheat.k12_standard.median_s; ADR-0099 accepted",
        f"metrics.json{SEP}te_real.te2.status",
    ]
    return [
        {"bullet": b, "status": s, "evidence": e}
        for b, s, e in zip(BULLETS, statuses, ev, strict=True)
    ]


def gate_block(
    a3_exit: int | None = 0,
    verdict: str | None = "pass",
    reason: str | None = None,
    a1_exit: int | None = 0,
) -> dict[str, Any]:
    """``metrics.json: evidence`` as ``compute`` writes it (no runbook)."""
    feas: dict[str, Any] = {"verdict": verdict, "reason": reason}
    if verdict is None and reason is None:
        feas["reason"] = "pending (M10-08 verdict)"
    return {
        "exit_codes": {"A1": a1_exit, "A2": 0, "A3": a3_exit},
        "a3_feasibility": feas,
        "benchmark": {
            "cells": {
                c: {"median_s": 3912.0, "commit": OTHER_COMMIT, "cpu": "6148"} for c in CELLS
            },
            "incomplete": [],
            "adr": "ADR-0099-benchmark-definition",
            "adr_status": "accepted",
        },
    }


def final_metrics(
    *,
    statuses: Sequence[str] | None = None,
    te: Mapping[str, str] | None = None,
    te7: str = "pass",
    rho: float = 0.83,
    gippa_commit: str = F,
    a3_exit: int | None = 0,
    verdict: str | None = "pass",
    reason: str | None = None,
) -> dict[str, Any]:
    """A minimal ``metrics.json`` for the fixture template: the checklist, ``evidence`` (exit
    codes, M10-08's verdict, the benchmark), ``te_real`` and ``te7`` in M10-06's shape
    (statuses only), one m4 number."""
    te = {"te2": "pass", "te3": "pass", "te4": "pass", "te5": "pass", **(te or {})}
    return {
        "schema": "gippa-pilot-metrics/1",
        "spec": "0" * 40,
        "gippa_commit": gippa_commit,
        "inputs": {},
        "warnings": [],
        "evidence": gate_block(a3_exit, verdict, reason),
        "checklist": checklist_rows(statuses, a3_exit),
        "m4": {"spearman_rho": rho, "ci_low": 0.76, "ci_high": 0.88},
        "te_real": {"schema": "gippa-pilot-te-real/1", **{k: {"status": v} for k, v in te.items()}},
        "te7": {"schema": "gippa-pilot-te7/1", "gippa": {"status": te7}, "diann": {"status": "pass"}},
    }


def diagnosis(te: str, cause: str = "biology") -> str:
    return (
        f"## Diagnosis: {te}\n\nCause class: {cause}\n\nThe gene is not expressed in this "
        "growth phase; stated with the evidence of the detection table.\n"
    )


NOTES_OK = "The pilot ran at the freeze commit.\n"


def render_final(
    metrics: Path,
    out: Path,
    *,
    notes: Path | None = None,
    template: Path | None = None,
    runbook: Path | None = None,
    final: bool = True,
) -> subprocess.CompletedProcess[str]:
    """``metrics.py render [--final] [--runbook R] …``, called in-process through
    ``metrics.main`` (the CLI's own entry point; a subprocess per call costs ~2 s of imports).
    An argparse error is exit 2."""
    import contextlib
    import io

    from helpers.pilot_metrics_m10 import load_metrics

    args = ["render", *(["--final"] if final else []), "--metrics", str(metrics), "--out", str(out)]
    if notes is not None:
        args += ["--notes", str(notes)]
    if template is not None:
        args += ["--template", str(template)]
    if runbook is not None:
        args += ["--runbook", str(runbook)]
    so, se = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(so), contextlib.redirect_stderr(se):
        try:
            code = load_metrics().main(args)
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else 1
    return subprocess.CompletedProcess(args, int(code or 0), so.getvalue(), se.getvalue())


SOURCE = re.compile(r"\s*\([^()]*\u203a[^()]*\)$")


def status_row(report: str, label: str) -> list[str]:
    """The cells of the status-table row whose first cell starts with ``label``; a trailing
    source ``(file SEP key)`` that render appends to a value is dropped."""
    for ln in report.splitlines():
        cells = [SOURCE.sub("", c.strip()) for c in ln.strip().strip("|").split("|")]
        if ln.startswith("|") and cells and cells[0].startswith(label):
            return cells
    raise AssertionError(f"no status row {label!r} in the report")


# ------------------------------------------------------------------ compute inputs (M10-06/07/10)


def write_te_inputs(d: Path) -> tuple[Path, Path]:
    """Hand-made ``te_real.json`` and ``te7.json`` (every TE passes) with every path of the
    final template (M10-06's shape; no row-level content)."""
    te_real = {
        "schema": "gippa-pilot-te-real/1",
        "te2": {"status": "pass", "counts": {"violations": 0}},
        "te3": {
            "status": "pass",
            "detection": {"share": 0.9},
            "abundance": {"percentile_rank": 0.8},
        },
        "te4": {
            "status": "pass",
            "pairs": {"rplL/rplJ": {"log2_ratio": 0.5}, "atpA/atpD": {"log2_ratio": 0.2}},
        },
        "te5": {"status": "pass", "A1": {"share_bucket": "low"}, "A2": {"share_bucket": "low"}},
    }
    te7 = {
        "schema": "gippa-pilot-te7/1",
        "gippa": {"status": "pass", "n_gate": 12},
        "diann": {"status": "pass"},
        "end_to_end": True,
    }
    return write_json(d / "te_real.json", te_real), write_json(d / "te7.json", te7)


def with_later_hooks(doc: dict[str, Any], m6_reason: str | None) -> dict[str, Any]:
    """``doc`` with the keys M10-10 (``m6``, ``h4``) and M10-07 (``m8``) will add: ``m6``/``h4``
    a reason-only finding when ``m6_reason`` is given (``a3_not_run``, ``different_strain``),
    else numbers (kept as computed when compute had the A3 inputs, M10-10b); ``m8`` with
    numbers."""
    level = {"fdp_combined": 0.012, "fdp_lower": 0.004, "assessment": "within_target"}
    out = dict(doc)
    computed = "reason" not in (doc.get("m6") or {"reason": None})
    if m6_reason is None and not computed:
        out["m6"] = {
            "overlap_a2_a3": {
                "novel_orfs": {"called": {"jaccard": 0.4}},
                "alt_sites": {"high": {"jaccard": 0.3}},
            }
        }
        out["h4"] = {"verdict": "holds"}
    elif m6_reason is not None:
        out["m6"] = {"reason": m6_reason}
        out["h4"] = {"verdict": "unestimated", "reason": m6_reason}
    out["m8"] = {"peptide": dict(level), "precursor": dict(level)}
    return out


# ------------------------------------------------------------------ the review packet sandbox


def run_packet(
    *args: str, root: Path, cwd: Path, script: Path | None = None, pythonpath: str | None = None
) -> subprocess.CompletedProcess[str]:
    """``python validation/pilot/review_packet.py ARGS`` with ``$GIPPA_ROOT``/``$GIPPA_DB_DIR``
    pointing into the sandbox. Fails the calling test when the script is absent. ``script`` and
    ``pythonpath``: run a copy of the script (the N6 test); the environment then has no
    ``PYTHONDONTWRITEBYTECODE``."""
    import pytest

    script = script or REVIEW_PACKET_PY
    if not script.is_file():
        pytest.fail(f"{script} does not exist")
    env = dict(os.environ)
    env["PYTHONPATH"] = pythonpath or os.pathsep.join(
        [str(REPO / "src"), str(REPO), str(REPO / "tests" / "sim")]
    )
    if pythonpath is not None:
        env.pop("PYTHONDONTWRITEBYTECODE", None)
        env.pop("PYTHONPYCACHEPREFIX", None)
    env["GIPPA_ROOT"] = str(root)
    env["GIPPA_DB_DIR"] = str(root / "db")
    env["GIPPA_DATA_DIR"] = str(root / "data")
    return subprocess.run(
        [sys.executable, str(script), *args],
        capture_output=True,
        text=True,
        env=env,
        cwd=str(cwd),
        check=False,
    )


@dataclass
class Sandbox:
    base: Path
    repo: Path
    root: Path
    out: Path
    metrics_dir: Path
    runbook: Path
    a1: Path
    a2: Path
    a3: Path | None
    cmp: Path
    verify: Path
    bench_summary: Path
    adrs: list[Path] = field(default_factory=list)

    @property
    def report(self) -> Path:
        return self.repo / "validation" / "pilot" / "REPORT.md"

    @property
    def metrics(self) -> Path:
        return self.repo / "validation" / "pilot" / "metrics" / "metrics.json"

    @property
    def notes(self) -> Path:
        return self.repo / "validation" / "pilot" / "report_notes.md"

    @property
    def template(self) -> Path:
        return self.repo / "validation" / "pilot" / "report_template.md"

    def args(self, round_k: int = 1) -> list[str]:
        out = [
            "--metrics-dir",
            str(self.metrics_dir),
            "--runbook",
            str(self.runbook),
            "--round",
            str(round_k),
            "--out",
            str(self.out),
            "--repo",
            str(self.repo),
            "--a1",
            str(self.a1),
            "--a2",
            str(self.a2),
            "--cmp",
            str(self.cmp),
            "--verify",
            str(self.verify),
            "--bench-summary",
            str(self.bench_summary),
        ]
        if self.a3 is not None:
            out += ["--a3", str(self.a3)]
        for adr in self.adrs:
            out += ["--adr", str(adr)]
        return out

    def run(self, round_k: int = 1, *extra: str) -> subprocess.CompletedProcess[str]:
        return run_packet(*self.args(round_k), *extra, root=self.root, cwd=self.base)

    def write_metrics(self, doc: Mapping[str, Any]) -> None:
        """Writes ``doc`` as the committed and the freeze ``metrics.json`` (same bytes)."""
        write_json(self.metrics, doc)
        (self.metrics_dir / "metrics.json").write_bytes(self.metrics.read_bytes())

    def rerender(self) -> None:
        """Re-render the committed report from the sandbox's committed metrics, notes and
        runbook."""
        tmp = self.base / "rendered.md"
        res = render_final(
            self.metrics, tmp, notes=self.notes, template=self.template, runbook=self.runbook
        )
        assert res.returncode == 0, f"render --final failed: {res.returncode}\n{res.stderr}"
        self.report.write_bytes(tmp.read_bytes())
        (self.metrics_dir / "REPORT.md").write_bytes(tmp.read_bytes())
        tmp.unlink()


def build_sandbox(
    base: Path,
    *,
    a3: bool = True,
    a3_exit: int = 0,
    statuses: Sequence[str] | None = None,
    runbook: str | None = None,
    notes: str = NOTES_OK,
    git: bool = False,
) -> Sandbox:
    """A checkout-like ``repo/`` with the committed files, a ``root/`` standing in for
    ``$GIPPA_ROOT`` and the packet output ``out/`` (not created). The committed ``REPORT.md`` is
    rendered by ``metrics.py render --final`` (so the builder needs the final mode).

    ``a3=False``: example 6 (M10-08 ``fail``/``route_unavailable``, runbook ``fallback`` without
    S11, no A3 run); else A3 ran with exit ``a3_exit`` on the Cascadia branch. ``git``: ``repo/``
    is a git checkout whose ``.gitignore`` ignores ``/local/``."""
    repo, root = base / "repo", base / "root"
    pilot = repo / "validation" / "pilot"
    exit3 = a3_exit if a3 else None
    write_text(pilot / "report_template.md", FIXTURE_TEMPLATE)
    write_json(
        pilot / "metrics" / "metrics.json",
        final_metrics(
            statuses=statuses,
            a3_exit=exit3,
            verdict="pass" if a3 else "fail",
            reason=None if a3 else "route_unavailable",
        ),
    )
    write_feasibility(
        pilot / "metrics" / "a3_feasibility.json",
        "pass" if a3 else "fail",
        None if a3 else "route_unavailable",
    )
    write_text(pilot / "report_notes.md", notes)
    rb = write_text(
        pilot / "reproduce.md",
        runbook
        if runbook is not None
        else runbook_text(
            a3_branch="cascadia" if a3 else "fallback",
            steps=None if a3 else [s for s in EVIDENCE_JOBS if s != "S11"],
        ),
    )
    adrs = [
        write_adr(
            repo / "docs" / "adr" / "ADR-0051-diann-version.md", "proposed", "DIA-NN 2.3.1 Academia"
        ),
        write_adr(repo / "docs" / "adr" / "ADR-0099-benchmark-definition.md", "accepted"),
    ]
    write_text(repo / "PROGRESS.md", "# PROGRESS\n\n## Questions for Khoa\n\n- none open\n")
    if git:
        write_text(repo / ".gitignore", "/local/\n__pycache__/\n")
        subprocess.run(["git", "init", "-q", str(repo)], check=True, capture_output=True)
    runs = root / "runs" / "pilot" / "freeze" / F7
    arms = {}
    for name in ("out_a1", "out_a2", *(["out_a3"] if a3 else [])):
        arms[name] = runs / name
        write_text(runs / name / "report" / "report.html", f"<html><body>{name}</body></html>\n")
    cmp = runs / "cmp"
    write_text(cmp / "report" / "compare.html", "<html><body>compare</body></html>\n")
    verify = write_json(
        root / "pilot" / "reproduce" / f"verify-{F7}.json", {"failures": [], "warnings": []}
    )
    bench = write_text(root / "bench" / "b1" / "summary.md", "| cell | median |\n|---|---|\n")
    metrics_dir = root / "pilot" / "metrics" / f"freeze-{F7}"
    metrics_dir.mkdir(parents=True)
    (metrics_dir / "metrics.json").write_bytes((pilot / "metrics" / "metrics.json").read_bytes())
    sb = Sandbox(
        base=base,
        repo=repo,
        root=root,
        out=base / "out" / "m10-pilot",
        metrics_dir=metrics_dir,
        runbook=rb,
        a1=arms["out_a1"],
        a2=arms["out_a2"],
        a3=arms.get("out_a3"),
        cmp=cmp,
        verify=verify,
        bench_summary=bench,
        adrs=adrs,
    )
    sb.rerender()
    return sb


def tree(root: Path, exclude: Path | None = None) -> dict[str, str]:
    """``relative path -> sha256`` of every file under ``root`` (outside ``exclude``)."""
    out = {}
    for p in sorted(root.rglob("*")):
        if p.is_file() and (exclude is None or exclude not in p.parents):
            out[str(p.relative_to(root))] = sha256(p)
    return out


def packet_files(sb: Sandbox, *, a3: bool = True) -> set[str]:
    """The packet table of the spec ("The review packet"), relative to the round directory."""
    files = {
        "index.md",
        "REPORT.md",
        "metrics.json",
        "report_notes.md",
        "runs/A1/report.html",
        "runs/A2/report.html",
        "runs/compare.html",
        "reproduce.md",
        f"verify-{F7}.json",
        "benchmark/summary.md",
        "a3_feasibility.json",
        "MANIFEST.sha256",
    }
    if a3:
        files.add("runs/A3/report.html")
    files |= {f"adrs/{p.name}" for p in sb.adrs}
    return files
