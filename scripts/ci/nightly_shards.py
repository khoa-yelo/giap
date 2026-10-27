#!/usr/bin/env python3
"""Fold the shards of nightly task 2 into its one ``calibration_nightly`` row.

Usage: ``nightly_shards.py record --dir D --shard K --exit RC --summary TEXT`` (one shard's result
as ``D/shard-K.json``) and ``nightly_shards.py merge --dir D --of N [--label calibration_nightly]``
(prints ``<exit>\\t<label>: <summary>`` and exits 0 once all N shards recorded, else exits 3 and
prints nothing). ``scripts/slurm/nightly.sbatch`` runs both under one flock, so exactly one shard
(the last to finish) sees the complete set and writes the night's row.

The merged summary keeps the pytest summary-line format that ``scripts/release/readiness.py``
reads: counts summed per kind (``N failed, M passed, …``), the longest shard's time; a shard
without a pytest summary line (killed, failed early) adds its own text. The exit is the largest
shard exit. Stdlib only.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

#: pytest's summary order (terminal.py ``KNOWN_TYPES``), singular and plural both summed
KINDS = ("failed", "passed", "skipped", "deselected", "xfailed", "xpassed", "warnings", "errors")
_COUNT = re.compile(r"(\d+) (failed|passed|skipped|deselected|xfailed|xpassed|warnings?|errors?)\b")
_TIME = re.compile(r"\bin (\d+(?:\.\d+)?)s\b")


def record(d: Path, shard: int, rc: int, summary: str) -> Path:
    """Write ``d/shard-<shard>.json`` atomically (a rerun of the same shard replaces it)."""
    d.mkdir(parents=True, exist_ok=True)
    out = d / f"shard-{shard}.json"
    tmp = out.with_suffix(f".tmp{os.getpid()}")
    tmp.write_text(json.dumps({"shard": shard, "exit": rc, "summary": summary}))
    tmp.replace(out)
    return out


def _body(summary: str) -> str:
    """``label: text`` → ``text`` (the shard rows carry the task label)."""
    head, sep, rest = summary.partition(":")
    return rest.strip() if sep and " " not in head.strip() else summary.strip()


def _hms(secs: float) -> str:
    s = int(secs)
    return f"{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}"


def fold(results: list[dict[str, object]], of: int) -> tuple[int, str]:
    """``(exit, summary)`` of the shards: summed pytest counts, longest time, other shard notes."""
    counts = dict.fromkeys(KINDS, 0)
    secs: list[float] = []
    notes: list[str] = []
    for r in sorted(results, key=lambda r: int(str(r["shard"]))):
        body = _body(str(r["summary"]))
        found = _COUNT.findall(body)
        t = _TIME.search(body)
        if found and t:
            for n, kind in found:
                counts[kind if kind.endswith("s") or kind in KINDS else kind + "s"] += int(n)
            secs.append(float(t[1]))
        else:
            notes.append(f"shard {r['shard']}/{of}: {body}")
    parts = [f"{counts[k]} {k}" for k in KINDS if counts[k]]
    text = ", ".join(parts)
    if secs:
        text += f" in {max(secs):.2f}s ({_hms(max(secs))}, longest of {len(secs)} shards)"
    text = "; ".join(x for x in [text, *notes] if x)
    rc = max(int(str(r["exit"])) for r in results)
    return rc, text


def merge(d: Path, of: int) -> tuple[int, str] | None:
    """The folded result once shards 1..``of`` all recorded, else None."""
    files = [d / f"shard-{k}.json" for k in range(1, of + 1)]
    if not all(f.is_file() for f in files):
        return None
    return fold([json.loads(f.read_text()) for f in files], of)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    rec = sub.add_parser("record")
    rec.add_argument("--dir", type=Path, required=True)
    rec.add_argument("--shard", type=int, required=True)
    rec.add_argument("--exit", type=int, required=True)
    rec.add_argument("--summary", required=True)
    mer = sub.add_parser("merge")
    mer.add_argument("--dir", type=Path, required=True)
    mer.add_argument("--of", type=int, required=True)
    mer.add_argument("--label", default="calibration_nightly")
    a = ap.parse_args(argv)
    if a.cmd == "record":
        record(a.dir, a.shard, a.exit, a.summary)
        return 0
    got = merge(a.dir, a.of)
    if got is None:
        return 3
    print(f"{got[0]}\t{a.label}: {got[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
