"""Summary of one DIA-NN conversion for the pilot (M10-03, docs/specs/M10-03-diann-convert.md).

``summarize`` is pure: report columns, the converted tables read back, and the converter's issues
in; counts, issues and rule failures out. The ``summary`` CLI reads the files, adds paths, hashes,
the giap version and the job, and writes ``<a>.convert.json`` (paths, hashes and counts only, no
row; Q-K1). Exit 1 when a summary rule fails.

    python -I validation/pilot/diann_convert.py summary --arm A1 --report R --peptides P \\
        --abundance B --log L --out J

``python -I`` ignores ``PYTHONPATH``, so this checkout's ``src/`` is put on ``sys.path`` here.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import socket
import subprocess
import sys
from pathlib import Path
from typing import Any

_REPO = Path(__file__).resolve().parents[2]
_SRC = _REPO / "src"
if not _SRC.is_dir() and os.environ.get("GIAP_REPO"):
    _SRC = Path(os.environ["GIAP_REPO"]) / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import pandas as pd  # noqa: E402

from giap.config.model import InputConfig  # noqa: E402
from giap.io.issues import Issue  # noqa: E402
from giap.io.tsv import format_value  # noqa: E402

SEARCH_DB = {"A1": "other", "A2": "same"}
Q_CUT = 0.01
_ISSUE = re.compile(r"^giap convert: (error|warn|info): (\S+?): (.*)$")
_ROWS = re.compile(r"\((\d+) rows\)")


def _decoy(col: pd.Series) -> pd.Series:
    return col.fillna(0).astype(int).astype(bool)


def summarize(
    report_cols: pd.DataFrame,
    peptides: pd.DataFrame,
    abundance: pd.DataFrame,
    issues: list[Issue],
) -> dict[str, Any]:
    """Counts over the converted table, the issue table and the failures of summary rules 1-4."""
    decoy = peptides["is_decoy"].fillna(False).astype(bool)
    q = pd.to_numeric(peptides["q_value"], errors="coerce")
    prefixes = tuple(InputConfig().contaminant_prefixes)
    n_cont = (
        peptides["protein_ids"]
        .fillna("")
        .astype(str)
        .map(lambda s: [i.startswith(prefixes) for i in s.split(";") if i])
    )
    rep_q = pd.to_numeric(report_cols["Global.Q.Value"], errors="coerce").astype("float64")
    written = rep_q.dropna().map(lambda v: float(format_value(float(v))))
    inten = pd.to_numeric(peptides["intensity"], errors="coerce")
    runs = peptides["run"].astype(str).value_counts()
    counts: dict[str, Any] = {
        "rows": len(peptides),
        "target_rows": int((~decoy).sum()),
        "decoy_rows": int(decoy.sum()),
        "decoy_rows_q_na": int((decoy & q.isna()).sum()),
        "contaminant_only_rows": int(n_cont.map(lambda f: bool(f) and all(f)).sum()),
        "mixed_contaminant_rows": int(n_cont.map(lambda f: any(f) and not all(f)).sum()),
        "rows_per_run": {k: int(runs[k]) for k in sorted(runs.index)},
        "distinct_peptides_target_q01": int(
            peptides.loc[~decoy & (q <= Q_CUT), "peptide"].nunique()
        ),
        "q_crossings": int(((rep_q.dropna() > Q_CUT) & (written <= Q_CUT)).sum()),
        "intensity_na_rows": int(inten.isna().sum()),
        "abundance_groups": len(abundance),
    }
    table: dict[str, dict[str, Any]] = {}
    for i in issues:
        entry = table.setdefault(i.code, {"severity": i.severity, "rows": 0})
        entry["rows"] += i.n_rows
        if i.severity == "error":
            entry["severity"] = "error"
    failures = []
    if counts["rows"] != len(report_cols):
        failures.append(f"rule 1: {counts['rows']} table rows, {len(report_cols)} report rows")
    rep_decoys = int(_decoy(report_cols["Decoy"]).sum()) if "Decoy" in report_cols else 0
    if counts["decoy_rows"] != rep_decoys or rep_decoys == 0:
        failures.append(f"rule 2: {counts['decoy_rows']} decoy rows, report has {rep_decoys}")
    errors = sorted(c for c, e in table.items() if e["severity"] == "error")
    if errors:
        failures.append(f"rule 3: error issues {errors}")
    rep_runs = sorted(set(report_cols["Run"].astype(str)))
    if rep_runs != sorted(counts["rows_per_run"]):
        failures.append(f"rule 4: runs {sorted(counts['rows_per_run'])} != report {rep_runs}")
    return {"counts": counts, "issues": table, "failures": failures}


def parse_log(text: str) -> list[Issue]:
    """The issue lines ``giap convert`` printed to stderr."""
    out = []
    for line in text.splitlines():
        m = _ISSUE.match(line)
        if m:
            sev, code, msg = m.groups()
            n = _ROWS.search(msg)
            out.append(Issue(code, sev, msg, int(n.group(1)) if n else 0))  # type: ignore[arg-type]
    return out


def _final(path: Path) -> Path:
    """``.<name>.tmp.<job>[.gz]`` → ``<name>``: the name the job moves the file to."""
    n = path.name
    return path.with_name(n[1:].split(".tmp.")[0]) if n.startswith(".") and ".tmp." in n else path


def _file(path: Path, rows: int) -> dict[str, Any]:
    data = path.read_bytes()
    return {
        "path": str(_final(path)),
        "sha256": hashlib.sha256(data).hexdigest(),
        "bytes": len(data),
        "rows": rows,
    }


def _giap() -> dict[str, Any]:
    import platform

    import giap

    repo = Path(os.environ.get("GIAP_REPO") or _REPO)

    def git(*a: str) -> str:
        r = subprocess.run(["git", "-C", str(repo), *a], capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else ""

    status = git("status", "--porcelain", "--untracked-files=no")
    return {
        "version": {"giap": giap.__version__, "python": platform.python_version()},
        "commit": git("rev-parse", "HEAD"),
        "dirty": bool(status),
    }


def _now() -> str:
    return dt.datetime.now(dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def main(argv: list[str] | None = None) -> int:
    import pyarrow.parquet as pq

    from giap.io.peptides import read_peptide_table

    ap = argparse.ArgumentParser(prog="diann_convert.py")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("summary")
    s.add_argument("--arm", required=True, choices=sorted(SEARCH_DB))
    for flag in ("--report", "--peptides", "--abundance", "--log", "--out"):
        s.add_argument(flag, required=True, type=Path)
    a = ap.parse_args(argv)

    schema = pq.read_schema(a.report).names
    cols = [c for c in ("Run", "Decoy", "Global.Q.Value") if c in schema]
    report = pq.read_table(a.report, columns=cols).to_pandas()
    peptides, _ = read_peptide_table(a.peptides)
    abundance = pd.read_csv(a.abundance, sep="\t")
    log = a.log.read_text()
    out = summarize(report, peptides, abundance, parse_log(log))
    command = ["giap", "convert", "--from", "diann", str(a.report), "--out"]
    command += [str(_final(a.peptides)), "--abundance-out", str(_final(a.abundance))]
    command += ["--search-db", SEARCH_DB[a.arm]]
    doc = {
        "arm": a.arm,
        "search_db": SEARCH_DB[a.arm],
        "report": _file(a.report, len(report)),
        "outputs": {
            "peptides": _file(a.peptides, len(peptides)),
            "abundance": _file(a.abundance, len(abundance)),
            "log": _file(a.log, len(log.splitlines())),
        },
        "giap": _giap(),
        "command": command,
        "job": {
            "id": os.environ.get("SLURM_JOB_ID", ""),
            "node": os.environ.get("SLURMD_NODENAME") or socket.gethostname(),
            "start": os.environ.get("GIAP_JOB_START", ""),
            "end": _now(),
        },
        **out,
    }
    a.out.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n")
    for f in out["failures"]:
        print(f"diann_convert summary: {f}", file=sys.stderr)
    return 1 if out["failures"] else 0


if __name__ == "__main__":
    sys.exit(main())
