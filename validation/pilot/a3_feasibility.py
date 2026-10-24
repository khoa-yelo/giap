"""A3 feasibility: is Cascadia 0.0.7 usable on the diaPASEF pilot runs? (M10-08).

    python validation/pilot/a3_feasibility.py mzml-check MZML \\
        --acquisition validation/pilot/inputs/acquisition.json --run WT_1 --out J

``mzml-check`` (08a) tests a converted mzML against the run's ``acquisition.json`` facts with the
rules M1-M8 of docs/specs/M10-08-a3-feasibility.md. It streams the file with the stdlib
``iterparse`` and never decodes a binary array (it reads ``defaultArrayLength`` only). ``J`` is
canonical JSON (sorted keys, 2-space indent, trailing newline) with ``rules``, ``pass``,
``failed`` and ``empty_ms2_share``. ``evaluate`` and ``early-fail`` come with 08b.

Exit codes: 0 pass; 1 a failed check or an unreadable mzML (a route failure, never a Cascadia
verdict); 2 a precondition (missing file, unknown run, unknown option). Never imports Cascadia.
"""

from __future__ import annotations

import argparse
import json
import sys
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

RULES = tuple(f"M{i}" for i in range(1, 9))
MS_LEVEL, CENTROID, PROFILE, SCAN_START = "MS:1000511", "MS:1000127", "MS:1000128", "MS:1000016"
TARGET, LOWER, UPPER = "MS:1000827", "MS:1000828", "MS:1000829"
# scan start time unit -> seconds, by unitAccession or else unitName; no unit is never seconds.
RT_SCALE = {"UO:0000031": 60.0, "minute": 60.0, "UO:0000010": 1.0, "second": 1.0}


@dataclass
class Summary:
    """What the rules need from one mzML (order-independent counts, sets and maxima)."""

    n_ms1: int = 0
    n_ms2: int = 0
    n_ms2_empty: int = 0
    n_bad_precursor: int = 0  # MS2 spectra without exactly one complete, positive window
    n_not_centroid: int = 0
    max_width: float = 0.0
    max_rt_s: float | None = None
    n_rt_no_unit: int = 0  # scan start times without a known unit (M7 fails)
    windows: set[tuple[float, float, float]] = field(default_factory=set)


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _params(elem: ET.Element) -> dict[str, ET.Element]:
    return {p.get("accession", ""): p for p in elem if _local(p.tag) == "cvParam"}


def _float(p: ET.Element | None) -> float | None:
    try:
        return float(p.get("value", "")) if p is not None else None
    except ValueError:
        return None


def _add_spectrum(s: Summary, spec: ET.Element) -> None:
    own = _params(spec)
    if CENTROID not in own or PROFILE in own:
        s.n_not_centroid += 1
    for scan in (e for e in spec.iter() if _local(e.tag) == "scan"):
        start = _params(scan).get(SCAN_START)
        rt = _float(start)
        if start is None or rt is None:
            continue
        scale = RT_SCALE.get(start.get("unitAccession") or start.get("unitName") or "")
        if scale is None:
            s.n_rt_no_unit += 1
            continue
        s.max_rt_s = rt * scale if s.max_rt_s is None else max(s.max_rt_s, rt * scale)
    level = own.get(MS_LEVEL)
    if level is None or level.get("value") == "1":
        s.n_ms1 += level is not None
        return
    s.n_ms2 += 1
    s.n_ms2_empty += int(spec.get("defaultArrayLength", "0") or 0) < 1  # ValueError: unreadable
    precursors = [e for e in spec.iter() if _local(e.tag) == "precursor"]
    good = 0
    for prec in precursors:
        win = next((e for e in prec if _local(e.tag) == "isolationWindow"), None)
        p = _params(win) if win is not None else {}
        t, lo, up = (_float(p.get(a)) for a in (TARGET, LOWER, UPPER))
        if t is None or lo is None or up is None:
            continue
        s.windows.add((round(t, 2), round(lo, 2), round(up, 2)))
        s.max_width = max(s.max_width, lo + up)
        good += lo > 0 and up > 0
    s.n_bad_precursor += not (len(precursors) == 1 and good == 1)


def summarise(path: Path) -> Summary:
    """Stream ``path``; raises ``ET.ParseError`` on malformed XML, ``ValueError`` on a
    non-integer ``defaultArrayLength``."""
    s = Summary()
    for _, elem in ET.iterparse(path, events=("end",)):
        if _local(elem.tag) == "spectrum":
            _add_spectrum(s, elem)
            elem.clear()
    return s


def check(s: Summary, acq: dict[str, Any]) -> dict[str, Any]:
    """The rules M1-M8 for one run's ``acquisition.json`` entry (pure)."""
    counts = acq["msms_type_counts"]
    n0, n9 = int(counts.get("0", 0)), int(counts.get("9", 0))
    rows, groups = int(acq["dia_window_rows"]), int(acq["window_groups"])
    lo_mz, hi_mz = (float(x) for x in acq["mz_range"])
    rules = {
        "M1": s.n_ms1 == n0,
        "M2": s.n_ms2 > 0 and n9 <= s.n_ms2 <= n9 * rows,
        "M3": s.n_bad_precursor == 0,
        "M4": s.max_width <= (hi_mz - lo_mz) / 2,
        "M5": groups <= len(s.windows) <= rows,
        "M6": s.n_not_centroid == 0,
        "M7": s.n_rt_no_unit == 0
        and s.max_rt_s is not None
        and s.max_rt_s >= 0.9 * float(acq["gradient_length_s"]),
        "M8": s.n_ms2_empty == 0,
    }
    failed = sorted(r for r, ok in rules.items() if not ok)
    return {
        "rules": rules,
        "pass": not failed,
        "failed": failed,
        "empty_ms2_share": s.n_ms2_empty / s.n_ms2 if s.n_ms2 else 0.0,
        "n_ms1": s.n_ms1,
        "n_ms2": s.n_ms2,
        "n_windows": len(s.windows),
        "max_isolation_width": s.max_width,
        "max_scan_start_s": s.max_rt_s,
    }


def _mzml_check(args: argparse.Namespace) -> int:
    for p in (args.mzml, args.acquisition):
        if not p.is_file():
            print(f"mzml-check: no such file: {p}", file=sys.stderr)
            return 2
    runs = {r["run"]: r for r in json.loads(args.acquisition.read_text())["runs"]}
    if args.run not in runs:
        print(f"mzml-check: run {args.run!r} not in {args.acquisition}", file=sys.stderr)
        return 2
    try:
        summary = summarise(args.mzml)
    except (ET.ParseError, ValueError) as exc:
        print(f"mzml-check: unreadable mzML {args.mzml}: {exc}", file=sys.stderr)
        return 1
    if summary.n_rt_no_unit:
        print(
            f"mzml-check: {summary.n_rt_no_unit} scan start times without a unit", file=sys.stderr
        )
    j = {"run": args.run, **check(summary, runs[args.run])}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(j, indent=2, sort_keys=True) + "\n")
    print(f"mzml-check {args.run}: " + " ".join(f"{r}={j['rules'][r]}" for r in RULES))
    if j["failed"]:
        print(f"mzml-check {args.run}: failed {', '.join(j['failed'])}", file=sys.stderr)
        return 1
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    mc = sub.add_parser("mzml-check", help="structural check of a converted mzML (M1-M8)")
    mc.add_argument("mzml", type=Path)
    mc.add_argument("--acquisition", type=Path, required=True)
    mc.add_argument("--run", required=True)
    mc.add_argument("--out", type=Path, required=True)
    try:
        args = ap.parse_args(argv)
    except SystemExit as exc:
        return int(exc.code or 0)
    return _mzml_check(args)


if __name__ == "__main__":
    sys.exit(main())
