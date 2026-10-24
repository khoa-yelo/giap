"""A3 feasibility: is Cascadia 0.0.7 usable on the diaPASEF pilot runs? (M10-08).

    python validation/pilot/a3_feasibility.py mzml-check MZML \\
        --acquisition validation/pilot/inputs/acquisition.json --run WT_1 --out J
    python validation/pilot/a3_feasibility.py evaluate RUN.ssl --run WT_1 --out J \\
        [--genome FNA --bakta-json JSON --contaminants FASTA --cutoff X] [provenance options]
    python validation/pilot/a3_feasibility.py early-fail --reason R --evidence TEXT --out J

``mzml-check`` (08a) tests a converted mzML against the run's ``acquisition.json`` facts with the
rules M1-M8 of docs/specs/M10-08-a3-feasibility.md. It streams the file with the stdlib
``iterparse`` and never decodes a binary array (it reads ``defaultArrayLength`` only). ``J`` is
canonical JSON (sorted keys, 2-space indent, trailing newline) with ``rules``, ``pass``,
``failed`` and ``empty_ms2_share``.

``evaluate`` (08b) applies the failure rule L6 to a Cascadia ``.ssl``: distinct stripped, I/L-folded
7-30-mers above the score cut-off (strict ``>``, default 0.8), the top 1,000 matched exactly
against the six-frame translation with ``giap.mapping`` (contaminant hits count as non-matches).
``early-fail`` records E1-E3. Both write the verdict JSON (aggregates only, no sequence).

Exit codes: 0 pass; 3 fail (L6 or early-fail); 1 an error (failed check, unreadable mzML or
``.ssl``, missing column, bad score); 2 a precondition (missing or empty file, unknown run,
unknown option or reason). Never imports Cascadia.
"""

import argparse
import csv
import datetime
import hashlib
import json
import math
import os
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from giap.annotate.bakta import parse_bakta_json
from giap.config import resolve_config
from giap.io.fasta import read_genome_fasta
from giap.io.protein_fasta import read_protein_fasta
from giap.mapping.contaminants import match_contaminant_peptides
from giap.mapping.decoys import build_decoys
from giap.mapping.exact import match_exact
from giap.mapping.sixframe import build_sixframe
from giap.rng import for_step
from giap.stats.cutoffs import wilson_ci

REPO = Path(__file__).resolve().parents[2]
NISSLE_FNA = REPO / "tests/data/genomes/nissle1917.fna.gz"
NISSLE_JSON = REPO / "tests/data/recorded/bakta/1.12.2-full/nissle/nissle1917.json.gz"
RULE = {"min_peptides": 500, "min_share": 0.20, "top_n": 1000}
EARLY = ("route_unavailable", "structurally_incompatible", "tool_incompatible")
STANDARD = frozenset("ACDEFGHKLMNPQRSTVWY")  # I is folded into L
_GROUPS, _MASSES = re.compile(r"\[[^\]]*\]|\([^)]*\)"), re.compile(r"[+-]\d+(?:\.\d+)?|-")

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


# ------------------------------------------------------------------ evaluate (08b)


class InputError(ValueError):
    """An unreadable ``.ssl`` (exit 1)."""


def fold(seq: str) -> str:
    """Spec step 2: drop ``[…]``/``(…)``, signed masses and ``-``; uppercase; keep A-Z; I -> L."""
    s = _MASSES.sub("", _GROUPS.sub("", seq)).upper()
    return re.sub(r"[^A-Z]", "", s).replace("I", "L")


@dataclass
class Tally:
    """Steps 1-7: row counts and the distinct passing peptides, ranked."""

    n_rows: int = 0
    n_len_7_30: int = 0
    n_below_cutoff: int = 0
    best: dict[str, tuple[float, str]] = field(default_factory=dict)  # folded -> (score, written)

    def ranked(self) -> list[str]:
        return sorted(self.best, key=lambda p: (-self.best[p][0], p))


def tally(rows: Sequence[tuple[str, float]], cutoff: float) -> Tally:
    """Pure: count, length-filter, cut off (strict ``>``) and collapse ``(sequence, score)``."""
    t = Tally(n_rows=len(rows))
    for written, score in rows:
        pep = fold(written)
        if not 7 <= len(pep) <= 30:
            continue
        t.n_len_7_30 += 1
        if score <= cutoff:
            t.n_below_cutoff += 1
        elif pep not in t.best or score > t.best[pep][0]:
            t.best[pep] = (score, t.best[pep][1] if pep in t.best else written)
    return t


def read_ssl(path: Path) -> list[tuple[str, float]]:
    """``(sequence, score)`` per row; columns by name; :class:`InputError` names the line."""
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not {"sequence", "score"} <= set(reader.fieldnames or ()):
            raise InputError(f"{path}: the header needs 'sequence' and 'score'")
        rows = []
        for r in reader:
            try:
                score = float(r["score"] or "")
            except ValueError:
                score = math.nan
            if not math.isfinite(score):
                raise InputError(f"{path} line {reader.line_num}: score {r['score']!r}")
            rows.append((r["sequence"] or "", score))
    return rows


def genome_hits(
    top: list[str], seqs: dict[str, str], replicons: Any
) -> tuple[set[str], set[str], Any]:
    """Steps 9-12 with ``giap.mapping``: target and decoy peptides, and the target hit table."""
    cfg = resolve_config()
    index = build_sixframe(seqs, replicons)
    decoys = build_decoys(index, "shuffle", cfg.mapping.decoys.shuffle_window_aa,
                          for_step("a3_feasibility", cfg.seed))  # fmt: skip
    hits = match_exact(top, index, {"shuffle": decoys}, dict.fromkeys(top, "shuffle"),
                       min_len=7, max_len=30)  # fmt: skip
    tgt = hits[hits["search"] == "target"]
    return set(tgt["peptide"]), set(hits.loc[hits["search"] == "decoy", "peptide"]), tgt


def verdict(n_pass: int, n_top: int, n_match: int) -> dict[str, Any]:
    """Pure: rule L6 (both rules say "fewer than"); ``low_genome_match`` is null at ``n_top`` 0."""
    few = n_pass < RULE["min_peptides"]
    low = None if n_top == 0 else 100 * n_match < 20 * n_top  # min_share 0.20, exact integers
    reason = {(False, False): None, (True, False): "too_few_peptides", (False, True):
              "low_genome_match", (True, True): "both"}[(few, bool(low))]  # fmt: skip
    return {"fails": {"too_few_peptides": few, "low_genome_match": low}, "reason": reason,
            "verdict": "pass" if reason is None else "fail"}  # fmt: skip


def _record(args: argparse.Namespace, **fields: Any) -> dict[str, Any]:
    """The verdict JSON skeleton (every field) with the job's provenance options."""
    check = json.loads(args.mzml_check.read_text())["rules"] if args.mzml_check else None
    nulls = ("genome_sha256", "cutoff", "cutoff_source", "n_rows", "n_below_cutoff", "n_len_7_30",
             "n_pass", "n_top", "n_match", "n_nonstandard", "share", "share_ci95", "decoy_share",
             "contaminant_share", "matches_by_replicon", "length_hist", "evidence")  # fmt: skip
    return {**dict.fromkeys(nulls), "schema_version": 1, "run": args.run,
            "cascadia_version": "0.0.7", "weights_sha256": args.weights_sha256,
            "route": args.route, "route_version": args.route_version,
            "mzml_sha256": args.mzml_sha256, "mzml_check": check, "rule": RULE,
            "job_ids": args.job_id, "date": args.date, **fields}  # fmt: skip


def _write(out: Path, j: dict[str, Any]) -> int:
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(j, indent=2, sort_keys=True) + "\n")
    print(f"a3_feasibility {j['run']}: verdict {j['verdict']} ({j['reason']})")
    return 0 if j["verdict"] == "pass" else 3


def _evaluate(args: argparse.Namespace) -> int:
    db = os.environ.get("GIAP_DB_DIR")
    cont = args.contaminants or (Path(db) / "contaminants/contaminants.fasta" if db else None)
    for p in (args.ssl, args.genome, args.bakta_json, cont):
        if p is None or not p.is_file() or p.stat().st_size == 0:
            print(f"evaluate: missing or empty: {p} (no Cascadia result for {args.run}?)",
                  file=sys.stderr)  # fmt: skip
            return 2
    try:
        t = tally(read_ssl(args.ssl), args.cutoff)
    except (InputError, csv.Error, UnicodeDecodeError) as exc:
        print(f"evaluate: {exc}", file=sys.stderr)
        return 1
    ranked = t.ranked()
    top = ranked[: RULE["top_n"]]
    seqs = {r.id: r.sequence for r in read_genome_fasta(args.genome)}  # step 8
    replicons = parse_bakta_json(args.bakta_json, genome=seqs).model.replicons
    target, decoy, tgt = genome_hits(top, seqs, replicons) if top else (set(), set(), None)
    contaminant = set(match_contaminant_peptides(top, read_protein_fasta(cont)))
    matched = target - contaminant  # conservative: a contaminant is never a genome match
    n_top, n_match = len(top), len(matched)
    by_rep = dict.fromkeys(replicons, 0)
    if tgt is not None:
        for contig, g in tgt[tgt["peptide"].isin(matched)].groupby("contig"):
            by_rep[str(contig)] = int(g["peptide"].nunique())
    if args.detail:
        args.detail.parent.mkdir(parents=True, exist_ok=True)
        head = "rank\tsequence\tfolded\tscore\tmatch\tcontaminant\tdecoy\n"
        args.detail.write_text(head + "".join(
            f"{i}\t{t.best[p][1]}\t{p}\t{t.best[p][0]}\t{p in matched:d}\t{p in contaminant:d}"
            f"\t{p in decoy:d}\n" for i, p in enumerate(top, 1)))  # fmt: skip
    share = n_match / n_top if n_top else None
    j = _record(
        args, **verdict(len(ranked), n_top, n_match),
        genome_sha256=hashlib.sha256(args.genome.read_bytes()).hexdigest(),
        cutoff=args.cutoff, cutoff_source="override" if args.cutoff_given else "cli_default",
        n_rows=t.n_rows, n_below_cutoff=t.n_below_cutoff, n_len_7_30=t.n_len_7_30,
        n_pass=len(ranked), n_top=n_top, n_match=n_match,
        n_nonstandard=sum(1 for p in ranked if set(p) - STANDARD), share=share,
        share_ci95=list(wilson_ci(n_match, n_top)) if n_top else None,
        decoy_share=len(decoy) / n_top if n_top else None,
        contaminant_share=len(contaminant) / n_top if n_top else None,
        matches_by_replicon=by_rep,
        length_hist={str(n): c for n, c in sorted(Counter(map(len, ranked)).items())},
    )  # fmt: skip
    return _write(args.out, j)


def _early_fail(args: argparse.Namespace) -> int:
    j = _record(args, fails={"too_few_peptides": None, "low_genome_match": None},
                verdict="fail", reason=args.reason, evidence=args.evidence)  # fmt: skip
    return _write(args.out, j)


def _provenance(p: argparse.ArgumentParser, run: str | None = None) -> None:
    """Filled in by the validation run from the job's environment and logs (spec, 08b)."""
    p.add_argument("--run", required=run is None, default=run)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--weights-sha256")
    p.add_argument("--route", choices=("B2", "B1"))
    p.add_argument("--route-version")
    p.add_argument("--mzml-sha256")
    p.add_argument("--mzml-check", type=Path, help="mzml-check's J (its rules are copied)")
    p.add_argument("--job-id", action="append", default=[])
    p.add_argument("--date", default=datetime.datetime.now(datetime.UTC).date().isoformat())


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    mc = sub.add_parser("mzml-check", help="structural check of a converted mzML (M1-M8)")
    mc.add_argument("mzml", type=Path)
    mc.add_argument("--acquisition", type=Path, required=True)
    mc.add_argument("--run", required=True)
    mc.add_argument("--out", type=Path, required=True)
    ev = sub.add_parser("evaluate", help="rule L6 on a Cascadia .ssl (verdict JSON)")
    ev.add_argument("ssl", type=Path)
    ev.add_argument("--genome", type=Path, default=NISSLE_FNA)
    ev.add_argument("--bakta-json", type=Path, default=NISSLE_JSON)
    ev.add_argument("--contaminants", type=Path, help="default $GIAP_DB_DIR/contaminants/...")
    ev.add_argument("--cutoff", type=float, help="Cascadia's -t (default its CLI default 0.8)")
    ev.add_argument("--detail", type=Path, help="per-peptide top list (under $GIAP_ROOT only)")
    _provenance(ev)
    ef = sub.add_parser("early-fail", help="record E1-E3 (verdict fail, no counts)")
    ef.add_argument("--reason", choices=EARLY, required=True)
    ef.add_argument("--evidence", required=True)
    _provenance(ef, run="WT_1")  # L6 is decided on run 1
    try:
        args = ap.parse_args(argv)
    except SystemExit as exc:
        return int(exc.code or 0)
    if args.cmd == "evaluate":
        args.cutoff_given = args.cutoff is not None
        args.cutoff = 0.8 if args.cutoff is None else args.cutoff
        return _evaluate(args)
    return _early_fail(args) if args.cmd == "early-fail" else _mzml_check(args)


if __name__ == "__main__":
    sys.exit(main())
