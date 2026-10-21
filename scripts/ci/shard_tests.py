#!/usr/bin/env python3
"""Split a pytest selection into stable shards by test id (CI-01 lead decision 1).

Usage: ``pytest --collect-only -q -m calibration_full | shard_tests.py --shard K --of N``
(or ``shard_tests.py --shard K --of N collected.txt``); the test ids of shard ``K`` are printed
one per line, in collection order. ``--count`` prints ``shard<TAB>tests`` for every shard instead.

``scripts/slurm/calibration_full.sbatch`` shards *by test id*, not by file: a single
``calibration_full`` test is the expensive unit (one attack arm runs for about an hour), so a file
is far too coarse. The assignment is ``sha1(test id) % shards``, which is stable across processes,
machines and runs (unlike :func:`hash`, which is salted per process), so every array task computes
the same partition from the same collection without a shared plan file. Renaming or adding a test
reshuffles only its own shard membership, never the whole partition's correctness: each shard runs
whatever lands in it, and together the shards are exactly the selection.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
from collections.abc import Iterable, Sequence
from pathlib import Path

NODE_ID = re.compile(r"^(?P<id>[^\s<>=]+\.py::\S.*?)\s*$")


def node_ids(text: str) -> list[str]:
    """The test ids in ``pytest --collect-only -q`` output, in collection order, deduplicated.

    Everything that is not a ``<path>.py::<test>`` line (blank lines, the ``N/M tests collected``
    footer, warnings) is ignored.
    """
    out: list[str] = []
    seen: set[str] = set()
    for raw in text.splitlines():
        match = NODE_ID.match(raw.strip())
        if match is None:
            continue
        node = match.group("id")
        if node not in seen:
            seen.add(node)
            out.append(node)
    return out


def shard_of(node_id: str, shards: int) -> int:
    """Shard of one test id: ``sha1(id) % shards``, the same everywhere and always."""
    if shards < 1:
        raise ValueError(f"shards must be >= 1, not {shards}")
    digest = hashlib.sha1(node_id.encode("utf-8"), usedforsecurity=False).digest()
    return int.from_bytes(digest, "big") % shards


def shard(node_ids_: Iterable[str], index: int, shards: int) -> list[str]:
    """The test ids of shard ``index`` of ``shards``, in the order they were collected."""
    if not 0 <= index < shards:
        raise ValueError(f"shard {index} outside 0..{shards - 1}")
    return [n for n in node_ids_ if shard_of(n, shards) == index]


def counts(node_ids_: Iterable[str], shards: int) -> list[int]:
    """How many tests each shard gets."""
    out = [0] * shards
    for node in node_ids_:
        out[shard_of(node, shards)] += 1
    return out


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("collected", nargs="?", type=Path, help="pytest --collect-only -q output")
    ap.add_argument("--shard", type=int, default=0, help="this shard's index (0-based)")
    ap.add_argument("--of", type=int, required=True, dest="shards", help="number of shards")
    ap.add_argument("--count", action="store_true", help="print the size of every shard instead")
    args = ap.parse_args(argv)
    text = args.collected.read_text() if args.collected else sys.stdin.read()
    ids = node_ids(text)
    if not ids:
        print("shard_tests: no test ids in the collection", file=sys.stderr)
        return 1
    try:
        if args.count:
            for i, n in enumerate(counts(ids, args.shards)):
                print(f"{i}\t{n}")
        else:
            for node in shard(ids, args.shard, args.shards):
                print(node)
    except ValueError as exc:
        print(f"shard_tests: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
