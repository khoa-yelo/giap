r"""TU pair agreement of ``giap run`` outputs with RegulonDB (M11-14 14s2, weekly fallback arm).

The truth is ``pair_truth.tsv`` of ``validation/regulondb_truth.py`` (GIAP_SPEC.md §15.4 step 2:
adjacent same-strand K-12 pairs labelled positive / negative / ambiguous on the recorded full-DB
Bakta gene model), the same pairs and labels 14s2 used for the operating point. A pair is called
co-transcribed when both genes lie in one TU of the run's ``tables/operons.tsv``. Over the
labelled (non-ambiguous) pairs: precision = TP / (TP + FP), recall = TP / (TP + FN).

Prints one line for the weekly row, e.g.::

    TU pairs vs RegulonDB 14.5.0: fallback P=0.828 R=0.916 (TP 1584 FP 328 FN 146); uniop P=...

A missing truth or output is reported as "not available" (exit 0: reported, never gated).

Usage::

    python scripts/tu_agreement.py --regulondb-dir $GIAP_DATA_DIR/regulondb/14.5.0 \
        --arm fallback=<outdir> --arm "uniop=<outdir>"
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def truth_labels(rows: Iterable[Mapping[str, str]]) -> dict[frozenset[str], bool]:
    """Labelled pairs of ``pair_truth.tsv``: ``{gene_a, gene_b} -> positive``; ambiguous dropped."""
    out: dict[frozenset[str], bool] = {}
    conflicting: set[frozenset[str]] = set()
    for r in rows:
        if r.get("label") in ("positive", "negative"):
            pair, positive = frozenset((r["gene_a"], r["gene_b"])), r["label"] == "positive"
            if out.get(pair, positive) != positive:
                conflicting.add(pair)  # labelled both ways: treated as ambiguous
            out[pair] = positive
    return {p: v for p, v in out.items() if p not in conflicting}


def tu_of(rows: Iterable[Mapping[str, str]]) -> dict[str, str]:
    """``locus_tag -> tu_id`` from ``operons.tsv`` (``genes`` joined with ``;`` or ``,``)."""
    out: dict[str, str] = {}
    for r in rows:
        if not r.get("tu_id"):
            raise ValueError("operons.tsv row without a tu_id")
        for g in (r.get("genes") or "").replace(",", ";").split(";"):
            if g.strip() and g.strip() != "NA":
                out[g.strip()] = r["tu_id"]
    return out


def agreement(truth: Mapping[frozenset[str], bool], tus: Mapping[str, str]) -> dict[str, float]:
    """TP, FP, FN, precision and recall of "both genes in one called TU" against ``truth``."""
    tp = fp = fn = 0
    for pair, positive in truth.items():
        a, b = sorted(pair)
        called = a in tus and tus[a] == tus.get(b)
        tp += called and positive
        fp += called and not positive
        fn += positive and not called
    prec = tp / (tp + fp) if tp + fp else float("nan")
    rec = tp / (tp + fn) if tp + fn else float("nan")
    return {"tp": tp, "fp": fp, "fn": fn, "precision": prec, "recall": rec}


def _arm(name: str, outdir: Path, truth: Mapping[frozenset[str], bool]) -> str:
    ops = outdir / "tables" / "operons.tsv"
    if not ops.is_file():
        return f"{name}: not available (no {ops})"
    try:
        a = agreement(truth, tu_of(_rows(ops)))
    except (OSError, ValueError, KeyError) as exc:
        return f"{name}: not available ({exc})"
    return (f"{name} P={a['precision']:.3f} R={a['recall']:.3f} "
            f"(TP {a['tp']} FP {a['fp']} FN {a['fn']})")  # fmt: skip


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--regulondb-dir", type=Path, required=True)
    ap.add_argument("--arm", action="append", default=[], help="NAME=OUTDIR (repeatable)")
    a = ap.parse_args(argv)
    if bad := [s for s in a.arm if "=" not in s]:
        ap.error(f"--arm needs NAME=OUTDIR: {bad[0]}")
    arms = [s.split("=", 1) for s in a.arm]
    path = a.regulondb_dir / "truth" / "pair_truth.tsv"
    if not path.is_file():
        names = "; ".join(f"{n[0]}: not available" for n in arms)
        print(f"TU pairs vs regulondb: not available (no {path}); {names}")
        return 0
    try:
        version = json.loads((a.regulondb_dir / "provenance.json").read_text())["regulondb_version"]
    except (OSError, ValueError, KeyError, TypeError):
        version = a.regulondb_dir.name
    truth = truth_labels(_rows(path))
    parts = [_arm(n[0], Path(n[-1]), truth) for n in arms]
    print(f"TU pairs vs RegulonDB {version} ({len(truth)} labelled): " + "; ".join(parts))
    return 0


if __name__ == "__main__":
    sys.exit(main())
