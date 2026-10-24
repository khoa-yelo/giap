"""The acquisition mode of every timsTOF run, read from ``<run>.d/analysis.tdf`` (M10-01).

Stdlib only (``sqlite3``), and the database is opened read-only. ``MsMsType`` in ``Frames``:
0 = MS1, 8 = ddaPASEF MS/MS, 9 = diaPASEF MS/MS (the open-source TDF readers' documentation).

    python -I validation/pilot/acquisition_mode.py RUN.d [RUN.d ...] --out acquisition.json

Exit 0 only when every run is diaPASEF and complete; 3 for a ddaPASEF dataset (the DIA-NN
settings must change, §17; stop before M10-02); 1 otherwise. The JSON is written in every case.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any
from urllib.parse import quote

METADATA_KEYS = (
    "InstrumentName",
    "InstrumentFamily",
    "AcquisitionSoftwareVersion",
    "AcquisitionDateTime",
)


def classify(counts: dict[int, int], dia_rows: int, pasef_rows: int) -> str:
    """The mode of one run from its frame-type counts and the PASEF table sizes."""
    present = {t for t, n in counts.items() if n}
    if 8 in present and 9 in present:
        return "mixed"
    if 9 in present and dia_rows:
        return "diaPASEF"
    if 8 in present and pasef_rows:
        return "ddaPASEF"
    if present == {0}:
        return "ms1_only"
    return "unknown"


def _count(con: sqlite3.Connection, tables: set[str], table: str) -> int:
    if table not in tables:
        return 0
    return int(con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0])


def read_tdf(tdf: Path) -> dict[str, Any]:
    """The acquisition facts of one ``analysis.tdf`` (never written to)."""
    con = sqlite3.connect(f"file:{quote(str(tdf.resolve()))}?mode=ro", uri=True)
    try:
        tables = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        counts: dict[int, int] = {}
        gradient = None
        if "Frames" in tables:
            q = "SELECT MsMsType, COUNT(*) FROM Frames GROUP BY MsMsType"
            counts = {int(t): int(n) for t, n in con.execute(q)}
            gradient = con.execute("SELECT MAX(Time) FROM Frames").fetchone()[0]
        dia_rows = _count(con, tables, "DiaFrameMsMsWindows")
        pasef_rows = _count(con, tables, "PasefFrameMsMsInfo")
        groups, mz_range = 0, None
        if dia_rows:
            groups, lo, hi = con.execute(
                "SELECT COUNT(DISTINCT WindowGroup), MIN(IsolationMz - IsolationWidth / 2),"
                " MAX(IsolationMz + IsolationWidth / 2) FROM DiaFrameMsMsWindows"
            ).fetchone()
            mz_range = [lo, hi]
        metadata: dict[str, str] = {}
        if "GlobalMetadata" in tables:
            rows = con.execute("SELECT Key, Value FROM GlobalMetadata").fetchall()
            metadata = {k: v for k, v in rows if k in METADATA_KEYS}
    finally:
        con.close()
    return {
        "mode": classify(counts, dia_rows, pasef_rows),
        "msms_type_counts": {str(t): n for t, n in sorted(counts.items()) if n},
        "dia_window_rows": dia_rows,
        "pasef_rows": pasef_rows,
        "window_groups": int(groups),
        "mz_range": mz_range,
        "gradient_length_s": gradient,
        "metadata": metadata,
    }


def inspect_run(run_dir: Path) -> dict[str, Any]:
    run = run_dir.name.removesuffix(".d")
    sidecar = run_dir.parent / f"{run_dir.name}.source_folder"
    source = sidecar.read_text().strip() if sidecar.is_file() else run_dir.name
    tdf = run_dir / "analysis.tdf"
    out: dict[str, Any] = {"run": run, "source_folder": source, "notes": []}
    if not tdf.is_file():
        out.update(mode="unknown", msms_type_counts={}, window_groups=0, mz_range=None,
                   gradient_length_s=None, metadata={})  # fmt: skip
        out["notes"].append("not a Bruker TDF run: no analysis.tdf")
    else:
        try:
            out.update(read_tdf(tdf))
        except sqlite3.DatabaseError as err:
            out.update(mode="unknown", msms_type_counts={}, window_groups=0, mz_range=None,
                       gradient_length_s=None, metadata={})  # fmt: skip
            out["notes"].append(f"analysis.tdf is not readable: {err}")
    out["incomplete_run"] = not (run_dir / "analysis.tdf_bin").is_file()
    if out["incomplete_run"]:
        out["notes"].append("analysis.tdf_bin missing: DIA-NN needs the binary")
    return out


def dataset_mode(runs: Sequence[dict[str, Any]]) -> str:
    modes = {r["mode"] for r in runs}
    return modes.pop() if len(modes) == 1 else "inconsistent"


def exit_code(doc: dict[str, Any]) -> int:
    if doc["dataset_mode"] == "ddaPASEF":
        return 3
    if doc["dataset_mode"] == "diaPASEF" and not any(r["incomplete_run"] for r in doc["runs"]):
        return 0
    return 1


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("runs", nargs="+", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    runs = [inspect_run(d) for d in args.runs]
    doc = {"dataset_mode": dataset_mode(runs), "runs": runs}
    args.out.write_text(json.dumps(doc, indent=2) + "\n")
    for r in runs:
        notes = "; ".join(r["notes"])
        print(f"{r['run']}\t{r['mode']}\t{r['msms_type_counts']}\t{notes}")
    print(f"dataset_mode: {doc['dataset_mode']}")
    code = exit_code(doc)
    if code == 3:
        print("ddaPASEF: settings must change: §17; stop before M10-02", file=sys.stderr)
    elif code == 1:
        print("the runs are not all complete diaPASEF runs: stop before M10-02", file=sys.stderr)
    return code


if __name__ == "__main__":
    sys.exit(main())
