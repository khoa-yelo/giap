"""Entry gate and conversion summary of the A3 run (M10-10a, docs/specs/M10-10-a3-run.md).

``gate`` checks entry rule 1 on M10-08's committed ``a3_feasibility.json``: exit 0 for a
well-formed ``pass`` whose weights sha256 equals the ``cascadia_weights`` row of
``a3_assets.tsv``; 3 for a well-formed ``fail`` (M10-10 is skipped); 2 for a missing or malformed
file or a weights mismatch. ``summarize`` is pure: the ``.ssl`` data-row counts, the converted
table, the converter's issues and the acquisition run labels in; counts, issues and the failures
of summary rules 1-4 out. The ``summary`` CLI adds paths, hashes, the gippa version and the job and
writes ``a3.convert.json`` (paths, hashes and counts only, no row). Exit 1 when a rule fails.

    python -I validation/pilot/cascadia_convert.py gate [--feasibility F] [--assets T]
    python -I validation/pilot/cascadia_convert.py summary --ssl-sums S --peptides P --log L \\
        --acquisition A --out J

``python -I`` ignores ``PYTHONPATH``, so this checkout (and its ``src/``) is put on ``sys.path``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import socket
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

import pandas as pd  # noqa: E402
from validation.pilot.diann_convert import _file, _final, _gippa, _now, parse_log  # noqa: E402

from gippa.io.converters.cascadia import CASCADIA_DEFAULT_CUTOFF  # noqa: E402
from gippa.io.issues import Issue  # noqa: E402

PILOT = _REPO / "validation" / "pilot"
FEASIBILITY = PILOT / "metrics" / "a3_feasibility.json"
ASSETS = PILOT / "a3_assets.tsv"
GATE_RUN = "WT_1"
FAIL_REASONS = frozenset(
    {
        "too_few_peptides",
        "low_genome_match",
        "both",
        "route_unavailable",
        "structurally_incompatible",
        "tool_incompatible",
    }
)
FORBIDDEN_COLUMNS = ("q_value", "pep", "is_decoy", "intensity")


class GateError(Exception):
    """Entry rule 1 cannot be checked (exit 2)."""


def weights_sha(assets: Path) -> str:
    """The sha256 of the ``cascadia_weights`` row of ``a3_assets.tsv`` (GateError if none)."""
    try:
        with assets.open(encoding="utf-8", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
    except OSError as exc:
        raise GateError(f"cannot read {assets}: {exc}") from exc
    shas = [(r.get("sha256") or "").strip() for r in rows if r.get("asset") == "cascadia_weights"]
    if len(shas) != 1 or not shas[0]:
        raise GateError(f"{assets}: no single cascadia_weights row with a sha256")
    return shas[0]


def gate_verdict(doc: object, assets_sha: str) -> str | None:
    """``None`` for a pass, the fail reason for a fail; GateError when malformed or mismatched."""
    if not isinstance(doc, dict):
        raise GateError("not a JSON object")
    sv = doc.get("schema_version")
    if type(sv) is not int or sv != 1:
        raise GateError(f"schema_version {sv!r}, expected 1")
    if doc.get("run") != GATE_RUN:
        raise GateError(f"run {doc.get('run')!r}, expected {GATE_RUN}")
    if "weights_sha256" not in doc:
        raise GateError("no weights_sha256")
    verdict, reason = doc.get("verdict"), doc.get("reason")
    if verdict == "fail":
        if reason not in FAIL_REASONS:
            raise GateError(f"verdict fail with an unknown reason {reason!r}")
        return str(reason)
    if verdict != "pass":
        raise GateError(f"verdict {verdict!r}, expected pass or fail")
    if reason is not None:
        raise GateError(f"verdict pass with reason {reason!r}")
    if doc["weights_sha256"] != assets_sha:
        raise GateError(f"weights_sha256 {doc['weights_sha256']!r} != a3_assets.tsv {assets_sha}")
    return None


def gate(feasibility: Path, assets: Path) -> int:
    try:
        try:
            doc = json.loads(feasibility.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise GateError(f"cannot read {feasibility}: {exc}") from exc
        reason = gate_verdict(doc, weights_sha(assets))
    except GateError as exc:
        print(f"cascadia_convert gate: M10-08 is not done ({exc}); submit nothing", file=sys.stderr)
        return 2
    if reason is not None:
        print(f"M10-08 verdict fail ({reason}): M10-10 is skipped; open M10-11a", file=sys.stderr)
        return 3
    print("cascadia_convert gate: M10-08 verdict pass", file=sys.stderr)
    return 0


def _issue_table(issues: Sequence[Issue]) -> dict[str, dict[str, Any]]:
    table: dict[str, dict[str, Any]] = {}
    for i in issues:
        entry = table.setdefault(i.code, {"severity": i.severity, "rows": 0})
        entry["rows"] += i.n_rows
        if i.severity == "error":
            entry["severity"] = "error"
    return table


def summarize(
    ssl_rows: Mapping[str, int],
    peptides: pd.DataFrame,
    issues: Sequence[Issue],
    acquisition_runs: Sequence[str],
) -> dict[str, Any]:
    """Counts over the converted table, the issue table and the failures of summary rules 1-4."""
    pep = peptides["peptide"].astype("string")
    mod = peptides["modified_peptide"].astype("string")
    score = pd.to_numeric(peptides["score"], errors="coerce").astype("float64")
    charge = pd.to_numeric(peptides["charge"], errors="coerce").astype("Int64")
    runs = peptides["run"].astype("string").fillna("NA").value_counts()
    charges = charge.dropna().astype(int).value_counts()
    lengths = pep.dropna().drop_duplicates().str.len()
    counts: dict[str, Any] = {
        "rows": len(peptides),
        "rows_per_run": {str(k): int(runs[k]) for k in sorted(runs.index)},
        "distinct_peptides": int(pep.nunique(dropna=True)),
        "distinct_peptides_7_30": int(((lengths >= 7) & (lengths <= 30)).sum()),
        "modified_rows": int((mod.fillna("") != pep.fillna("")).sum()),
        "charge_counts": {str(k): int(charges[k]) for k in sorted(charges.index)},
        "score_min": float(score.min()) if score.notna().any() else None,
        "score_median": float(score.median()) if score.notna().any() else None,
        "score_max": float(score.max()) if score.notna().any() else None,
        "score_at_or_below_default_cutoff": int((score <= CASCADIA_DEFAULT_CUTOFF).sum()),
    }
    table = _issue_table(issues)
    failures = []
    n_ssl = sum(ssl_rows.values())
    if counts["rows"] != n_ssl:
        failures.append(f"rule 1: {counts['rows']} table rows, {n_ssl} .ssl data rows")
    if sorted(counts["rows_per_run"]) != sorted(acquisition_runs):
        failures.append(
            f"rule 2: runs {sorted(counts['rows_per_run'])} != acquisition {list(acquisition_runs)}"
        )
    errors = sorted(c for c, e in table.items() if e["severity"] == "error")
    if errors:
        failures.append(f"rule 3: error issues {errors}")
    sources = sorted(set(peptides["source"].astype("string").fillna("NA")) - {"denovo"})
    extra = [c for c in FORBIDDEN_COLUMNS if c in peptides.columns]
    if sources or extra:
        failures.append(f"rule 4: sources other than denovo {sources}, columns {extra}")
    return {"counts": counts, "issues": table, "failures": failures}


def ssl_entries(sums: Path) -> tuple[dict[str, dict[str, Any]], list[str]]:
    """``{run: {path, sha256, bytes, rows}}`` from a ``sha256sum`` file; mismatches as failures."""
    entries: dict[str, dict[str, Any]] = {}
    failures = []
    for line in sums.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        want, name = line.split(maxsplit=1)
        path = Path(name.removeprefix("*"))
        data = path.read_bytes()
        got = hashlib.sha256(data).hexdigest()
        if got != want:
            failures.append(f"ssl: sha256 of {path} is {got}, {sums.name} lists {want}")
        rows = sum(1 for ln in data.decode("utf-8").splitlines()[1:] if ln.strip())
        run = path.name.removesuffix(".ssl")
        entries[run] = {"path": str(path), "sha256": got, "bytes": len(data), "rows": rows}
    return entries, failures


def acquisition_runs(path: Path) -> list[str]:
    """The run labels of ``acquisition.json``, in file order."""
    return [str(r["run"]) for r in json.loads(path.read_text(encoding="utf-8"))["runs"]]


def summary(a: argparse.Namespace) -> int:
    from gippa.io.peptides import read_peptide_table

    ssl, ssl_failures = ssl_entries(a.ssl_sums)
    peptides, _ = read_peptide_table(a.peptides)
    log = a.log.read_text(encoding="utf-8")
    out = summarize(
        {r: e["rows"] for r, e in ssl.items()},
        peptides,
        parse_log(log),
        acquisition_runs(a.acquisition),
    )
    out["failures"] = ssl_failures + out["failures"]
    command = ["gippa", "convert", "--from", "cascadia", *(e["path"] for e in ssl.values())]
    command += ["--out", str(_final(a.peptides)), "--source", "denovo"]
    doc = {
        "ssl": ssl,
        "outputs": {
            "peptides": _file(a.peptides, len(peptides)),
            "log": _file(a.log, len(log.splitlines())),
        },
        "gippa": _gippa(),
        "command": command,
        "job": {
            "id": os.environ.get("SLURM_JOB_ID", ""),
            "node": os.environ.get("SLURMD_NODENAME") or socket.gethostname(),
            "start": os.environ.get("GIPPA_JOB_START", ""),
            "end": _now(),
        },
        **out,
    }
    a.out.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for f in out["failures"]:
        print(f"cascadia_convert summary: {f}", file=sys.stderr)
    return 1 if out["failures"] else 0


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="cascadia_convert.py")
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gate", help="entry rule 1 (M10-08 verdict)")
    g.add_argument("--feasibility", type=Path, default=FEASIBILITY)
    g.add_argument("--assets", type=Path, default=ASSETS)
    s = sub.add_parser("summary", help="a3.convert.json")
    for flag in ("--ssl-sums", "--peptides", "--log", "--acquisition", "--out"):
        s.add_argument(flag, required=True, type=Path)
    a = ap.parse_args(list(argv) if argv is not None else None)
    if a.cmd == "gate":
        return gate(a.feasibility, a.assets)
    return summary(a)


if __name__ == "__main__":
    sys.exit(main())
