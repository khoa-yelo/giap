#!/usr/bin/env python
"""Sharded TerminatorNet scan: plan -> chunk -> merge (docs/specs/M2-15b-terminatornet-sharding.md).

Runs in the ``gippa-terminatornet`` env (Python 3.10), not part of ``src/gippa``. Upstream
``TerminatorNet`` and ``Features`` (MIT; fetched, not vendored, R5) are imported from
``$TERMINATORNET_HOME`` inside the subcommands only; the planning and ordering functions are pure.

* ``plan [--core N]``: ``genome.fa`` -> ``plan.tsv``, one row per chunk core of windows;
* ``chunk <chunk_id> --cpus N``: scans the core's windows on the whole cleaned strand sequence with
  upstream's own pool functions -> ``calls_<chunk_id>.pkl`` (protocol 4);
* ``merge``: every block (+ blocks in FASTA order, then - blocks) through upstream
  ``merge_overlapping`` / ``output_terms`` -> ``terminatornet.tsv``, byte-identical to
  ``TerminatorNet.py -s -n N`` (N >= 2).
"""

from __future__ import annotations

import argparse
import csv
import multiprocessing
import os
import pickle
import sys
from collections.abc import Sequence
from operator import itemgetter
from pathlib import Path
from types import ModuleType
from typing import Any

CORE = 250_000  # windows per chunk core (spec lead decision 5; not a config key)
WINDOW = 60  # upstream Features.WINDOW = HAIRPIN_LENGTH + TAIL_LENGTH, STEP = 1
PLAN_COLUMNS = (
    "chunk_id",
    "seq_index",
    "seq_name",
    "strand",
    "seq_length",
    "core_start",
    "core_end",
    "mode",
)
_INT_COLUMNS = ("seq_index", "seq_length", "core_start", "core_end")
GENOME, PLAN, MODEL, OUT = "genome.fa", "plan.tsv", "model.pkl", "terminatornet.tsv"

PlanRow = dict[str, Any]
Call = tuple[Any, ...]


# ---------------------------------------------------------------- pure functions


def nwin(length: int) -> int:
    """Upstream's window count per strand sequence: ``range(0, max(L - WINDOW + 1, 1))``."""
    return max(length - WINDOW + 1, 1)


def chunk_id(seq_index: int, strand: str, k: int) -> str:
    """``c{seq_index:05d}{p|m}{k:04d}``."""
    return f"c{seq_index:05d}{'p' if strand == '+' else 'm'}{k:04d}"


def plan_rows(records: Sequence[tuple[str, int]], mode: str, core: int = CORE) -> list[PlanRow]:
    """One row per core ``[k*core, min((k+1)*core, nwin))`` of each (record, strand), ordered by
    (seq_index, + before -, k); ``records`` = ``[(name, cleaned_length), ...]`` as upstream keeps
    them."""
    if core < 1:
        raise ValueError(f"core must be >= 1, got {core}")
    rows: list[PlanRow] = []
    for i, (name, length) in enumerate(records):
        n = nwin(length)
        for strand in "+-":
            for k, start in enumerate(range(0, n, core)):
                rows.append(
                    {
                        "chunk_id": chunk_id(i, strand, k),
                        "seq_index": i,
                        "seq_name": name,
                        "strand": strand,
                        "seq_length": length,
                        "core_start": start,
                        "core_end": min(start + core, n),
                        "mode": mode,
                    }
                )
    return rows


def output_blocks(rows: Sequence[PlanRow]) -> list[list[PlanRow]]:
    """Rows grouped per (seq_index, strand) in upstream output order (every + block in FASTA
    order, then every - block), each group in core (k) order."""
    groups: dict[tuple[int, int], list[PlanRow]] = {}
    for r in rows:
        groups.setdefault(("+-".index(r["strand"]), int(r["seq_index"])), []).append(r)
    return [sorted(groups[key], key=lambda r: int(r["core_start"])) for key in sorted(groups)]


def order_calls(calls_by_chunk: Sequence[Sequence[Call]], mode: str) -> list[Call]:
    """The chunks' calls concatenated in the given (core) order; in ``long`` mode stable-sorted by
    stop, as upstream ``scan_multiprocess`` does (``short`` = upstream ``scan``: window order)."""
    calls = [c for chunk in calls_by_chunk for c in chunk]
    if mode == "long":
        calls.sort(key=itemgetter(1))
    return calls


# ---------------------------------------------------------------- upstream and files


def _fail(msg: str) -> None:
    sys.stderr.write(f"\nError - {msg}\n\n")
    sys.exit(1)


def _upstream() -> tuple[ModuleType, ModuleType]:
    home = os.environ.get("TERMINATORNET_HOME") or os.path.join(
        os.environ.get("CONDA_PREFIX", ""), "share", "terminatornet"
    )
    sys.path.insert(0, home)
    import Features  # type: ignore[import-not-found]
    import TerminatorNet  # type: ignore[import-not-found]

    return TerminatorNet, Features


def _genome(tn: ModuleType) -> tuple[list[str], list[str], bool]:
    """Cleaned sequences, names and CONTAINS_LONG_SEQUENCE, exactly as upstream reads them."""
    seqs, names = tn.read_in_sequences(GENOME)
    if len(seqs) == 0:
        _fail(f"could not read in genomic sequences from file {GENOME}")
    seqs, _gc, long_seq = tn.clean_sequences(seqs)
    return seqs, names, bool(long_seq)


def _read_plan() -> list[PlanRow]:
    with open(PLAN, newline="") as fh:
        rows: list[PlanRow] = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        for c in _INT_COLUMNS:
            r[c] = int(r[c])
    return rows


# ---------------------------------------------------------------- subcommands


def cmd_plan(core: int) -> None:
    tn, _ = _upstream()
    seqs, names, long_seq = _genome(tn)
    rows = plan_rows(
        list(zip(names, map(len, seqs), strict=True)), "long" if long_seq else "short", core
    )
    with open(PLAN, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(PLAN_COLUMNS)
        w.writerows([r[c] for c in PLAN_COLUMNS] for r in rows)


def cmd_chunk(cid: str, cpus: int) -> None:
    import numpy as np

    tn, features = _upstream()
    row = next((r for r in _read_plan() if r["chunk_id"] == cid), None)
    if row is None:
        _fail(f"chunk {cid} is not in {PLAN}")
        return
    seqs, _names, _long = _genome(tn)
    seq = seqs[row["seq_index"]]
    # Upstream scan_multiprocess: module globals read by the pool functions (inherited by fork).
    features.SEQUENCE = seq if row["strand"] == "+" else features.reverse_complement(seq)
    features.GC_PERCENT = 0.5
    with open(MODEL, "rb") as fh:
        model, scaler = pickle.load(fh)
    ctx = multiprocessing.get_context("fork")
    with ctx.Pool(cpus) as pool:
        candidates = pool.map(
            features.scan_oneprocess_1, list(range(row["core_start"], row["core_end"]))
        )
    candidates2 = [(hrpn, j, tt) for hrpn, j, tt in candidates if j >= 0]
    with ctx.Pool(cpus) as pool:
        candidates3 = pool.map(features.scan_oneprocess_2, candidates2)
    calls: list[Call] = []
    for feats, hairpin, u_tract, a_tract, _tt, j in candidates3:
        x = scaler.transform(np.array([feats[f] for f in features.FEATURES]).reshape(1, -1))
        pred = model.predict_proba(x)[0, 1]  # row by row, as upstream (no batching)
        if pred >= 0.5:
            call = (
                j + feats["hairpin_start"],
                j + feats["hairpin_end"],
                a_tract,
                hairpin,
                u_tract,
                feats["structure"],
                feats["GH"],
                pred,
            )
            # Pool batching decides which equal strings share an object, and pickle's memo is
            # identity-based: intern them so the bytes do not depend on the CPU count.
            calls.append(tuple(sys.intern(v) if type(v) is str else v for v in call))
    payload = {
        "chunk_id": cid,
        "core_start": row["core_start"],
        "core_end": row["core_end"],
        "calls": calls,
    }
    with open(f"calls_{cid}.pkl", "wb") as fh:
        pickle.dump(payload, fh, protocol=4)


def cmd_merge() -> None:
    rows = _read_plan()
    plan_ids = {r["chunk_id"]: r for r in rows}
    files = {p.name[len("calls_") : -len(".pkl")]: p for p in Path.cwd().glob("calls_*.pkl")}
    missing, extra = sorted(plan_ids.keys() - files.keys()), sorted(files.keys() - plan_ids.keys())
    if missing or extra:
        _fail(f"calls pickles do not match {PLAN}: missing {missing}, extra {extra}")
    payloads: dict[str, dict[str, Any]] = {}
    for cid, r in plan_ids.items():
        with files[cid].open("rb") as fh:
            p = pickle.load(fh)
        got = (p.get("chunk_id"), p.get("core_start"), p.get("core_end"))
        if got != (cid, r["core_start"], r["core_end"]):
            _fail(f"calls_{cid}.pkl holds {got}, plan has {(cid, r['core_start'], r['core_end'])}")
        payloads[cid] = p
    tn, _ = _upstream()
    seqs, names, _long = _genome(tn)
    tn.BOTH_STRANDS = True
    with open(OUT, "w") as out:
        tn.OUT = out
        for block in output_blocks(rows):
            first = block[0]
            calls = order_calls([payloads[r["chunk_id"]]["calls"] for r in block], first["mode"])
            i = first["seq_index"]
            if first["strand"] == "+":
                terms = tn.merge_overlapping(calls)
            else:
                terms = tn.merge_overlapping(calls, len(seqs[i]), "-")
            tn.output_terms(terms, names[i])


def main(argv: Sequence[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("plan", help=f"{GENOME} -> {PLAN}")
    p.add_argument("--core", type=int, default=CORE, help="windows per chunk core")
    c = sub.add_parser("chunk", help="scan one chunk -> calls_<chunk_id>.pkl")
    c.add_argument("chunk_id")
    c.add_argument("--cpus", type=int, default=1)
    sub.add_parser("merge", help=f"calls_*.pkl -> {OUT}")
    args = ap.parse_args(argv)
    if args.cmd == "plan":
        cmd_plan(args.core)
    elif args.cmd == "chunk":
        cmd_chunk(args.chunk_id, args.cpus)
    else:
        cmd_merge()


if __name__ == "__main__":
    main()
