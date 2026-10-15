"""Promoter Calculator test fixtures (M2-14; test author, independent of ``src/``).

* ``HEADER``: the 27-column CSV header of Promoter Calculator 1.2.4 (recorded; also in
  ``tests/data/recorded/headers/promoter_calculator.csv.txt``).
* ``write_pc_csv``: a hand-written tool CSV (CRLF like the tool, or LF like the nf stub).
* ``raw_tss``: the value the *tool* writes in its ``TSS`` column for a 1-based region position.
  Promoter Calculator 1.2.4 writes, for ``+`` hits, the 0-based index of the +1 base
  (``ITR = seq[TSS:TSS+20]``) and, for ``-`` hits, ``len - TSS_rev`` = the 1-based forward
  coordinate of the +1 base (``wrapper.py``; checked on all 244,836 recorded rows). ``PcHit.tss`` is
  1-based (spec), so the parser adds 1 to ``+`` rows only.
* ``recorder_tus``: the TU regions of ``scripts/dev/record_tools_inputs.py regions`` (distance
  model, -300..+50, linear clipping) rebuilt from a Bakta JSON, as ``Tu`` objects; reproduces every
  recorded ``regions/<TU>.fa`` of sim100k and K-12 byte for byte.
* ``read_tus_tsv``: ``tests/data/recorded/promoter_calculator/1.2.4/<sample>/tus.tsv`` → ``Tu``.
"""

from __future__ import annotations

import collections
import csv
import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from giap.annotate.operons import Tu
from giap.model.coords import Interval, Replicon

HEADER = (
    "promoter_sequence,TSS,UP,hex35,spacer,hex10,disc,ITR,dG_total,dG_10,dG_35,dG_disc,dG_ITR,"
    "dG_ext10,dG_spacer,dG_bind,dG_UP,Tx_rate,UP_position,hex35_position,spacer_position,"
    "hex10_position,disc_position,fraction,strand,drop,length"
)
COLUMNS = HEADER.split(",")
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
RECORDED = Path(__file__).resolve().parents[1] / "data" / "recorded" / "promoter_calculator"


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def raw_tss(tss1: int, strand: str) -> int:
    """Tool ``TSS`` column for a hit whose +1 base is the 1-based region position ``tss1``."""
    return tss1 - 1 if strand == "+" else tss1


def _row(tss_raw: int, strand: str, tx_rate: float) -> list[str]:
    vals = {c: "0" for c in COLUMNS}
    vals.update(
        promoter_sequence="A" * 83,
        TSS=str(tss_raw),
        UP="A" * 24,
        hex35="TTGACA",
        spacer="A" * 17,
        hex10="TATAAT",
        disc="AAAAAAA",
        ITR="A" * 20,
        Tx_rate=repr(float(tx_rate)),
        UP_position="[1, 25]",
        hex35_position="[25, 31]",
        spacer_position="[31, 48]",
        hex10_position="[48, 54]",
        disc_position="[54, 61]",
        fraction="1",
        strand=strand,
        drop="False",
        length="83",
    )
    return [vals[c] for c in COLUMNS]


def write_pc_csv(
    path: Path, hits: Sequence[tuple[int, str, float]], *, crlf: bool = True, raw: bool = False
) -> Path:
    """Write a tool-format CSV; ``hits`` = (1-based region TSS, strand, Tx_rate) unless ``raw``
    (then the first item is written verbatim to the ``TSS`` column)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\r\n" if crlf else "\n")
        w.writerow(COLUMNS)
        for t, s, x in hits:
            w.writerow(_row(t if raw else raw_tss(t, s), s, x))
    return path


def header_only(path: Path, *, crlf: bool = True) -> Path:
    """The file ``PROMOTER_CALCULATOR`` writes when the tool found no promoter (M2-07)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes((HEADER + ("\r\n" if crlf else "\n")).encode())
    return path


def read_dicts(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def region_seq(fa: Path) -> str:
    return "".join(ln.strip() for ln in fa.read_text().splitlines() if not ln.startswith(">"))


def itr_at(region: str, tss1: int, strand: str) -> str:
    """20 nt starting at the +1 base ``tss1`` (1-based, forward region coordinates)."""
    if strand == "+":
        return region[tss1 - 1 : tss1 - 1 + 20]
    n = len(region)
    return revcomp(region)[n - tss1 : n - tss1 + 20]


def genome_itr(seq: str, tss0: int, strand: str) -> str:
    """20 nt transcribed from the 0-based genome TSS ``tss0`` on ``strand`` (circular read)."""
    ring = seq + seq
    if strand == "+":
        return ring[tss0 : tss0 + 20]
    lo = tss0 - 19
    return revcomp(ring[lo + len(seq) : tss0 + 1 + len(seq)] if lo < 0 else ring[lo : tss0 + 1])


def make_tu(tu_id: str, contig: str, strand: str, tu: tuple[int, int], up: tuple[int, int]) -> Tu:
    return Tu(
        tu_id=tu_id,
        contig=contig,
        interval=Interval(contig, tu[0], tu[1], strand),  # type: ignore[arg-type]
        strand=strand,
        genes=(f"{tu_id}_g1",),
        upstream=Interval(contig, up[0], up[1], strand),  # type: ignore[arg-type]
    )


def read_tus_tsv(path: Path) -> tuple[list[Tu], dict[str, Replicon], dict[str, int]]:
    """``tus.tsv`` → TUs, replicons and the 0-based first base of each TU's first start codon."""
    tus, reps, sc = [], {}, {}
    with path.open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            reps[r["contig"]] = Replicon(r["contig"], int(r["contig_length"]), r["topology"])  # type: ignore[arg-type]
            tus.append(
                make_tu(
                    r["tu_id"], r["contig"], r["strand"],
                    (int(r["tu_start"]), int(r["tu_end"])),
                    (int(r["region_start"]), int(r["region_end"])),
                )
            )  # fmt: skip
            sc[r["tu_id"]] = int(r["start_codon"])
    return tus, reps, sc


def recorder_tus(
    bakta_json: Path, up: int = 300, down: int = 50, d_max: int = 50
) -> tuple[list[Tu], dict[str, Replicon], dict[str, str]]:
    """TUs and regions exactly as ``record_tools_inputs.py regions`` built them; also returns the
    contig sequences (upper case)."""
    data: dict[str, Any] = json.loads(bakta_json.read_text())
    seqs = {s["id"]: s for s in data["sequences"]}
    genes = [f for f in data["features"] if f["type"] in {"cds", "sorf"}]
    genes.sort(key=lambda f: (f.get("sequence") or f.get("contig"), f["start"], f["stop"]))
    by: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for g in genes:
        by[str(g.get("sequence") or g.get("contig"))].append(g)
    tus: list[Tu] = []
    reps: dict[str, Replicon] = {}
    nts: dict[str, str] = {}
    n = 0
    for cid, gs in by.items():
        nt = str(seqs[cid]["nt"]).upper()
        size = len(nt)
        nts[cid] = nt
        reps[cid] = Replicon(cid, size, "linear")  # the recorder clipped linearly
        runs: list[list[dict[str, Any]]] = []
        for g in gs:
            prev = runs[-1][-1] if runs else None
            if (
                prev is not None
                and prev["strand"] == g["strand"]
                and (g["start"] - prev["stop"] - 1 <= d_max)
            ):
                runs[-1].append(g)
            else:
                runs.append([g])
        for run in runs:
            n += 1
            tid = f"TU_{n:05d}"
            if run[0]["strand"] == "+":
                s0 = run[0]["start"] - 1
                ends = [o["stop"] for o in gs if o["stop"] <= s0]
                a, b = max(s0 - up, max(ends, default=0), 0), min(s0 + down, size)
                tus.append(make_tu(tid, cid, "+", (s0, run[-1]["stop"]), (a, b)))
            else:
                e0 = run[-1]["stop"]
                starts = [o["start"] - 1 for o in gs if o["start"] - 1 >= e0]
                b = min(e0 + up, min(starts, default=size), size)
                tus.append(
                    make_tu(tid, cid, "-", (run[0]["start"] - 1, e0), (max(e0 - down, 0), b))
                )
    return tus, reps, nts
