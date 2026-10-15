"""Regenerate the Bakta JSON excerpts used by the M2-19 region tests (run by hand, not by pytest).

Usage: ``python tests/data/recorded/regions/regions_make_excerpts.py`` with ``$GIAP_ROOT`` set.
Selection (see PROVENANCE.md): every CDS/sORF overlapping a geNOMAD provirus or an antiSMASH
area, the nearest gene on each side outside them, every gene within 3 kb of a CRISPR array, every
gene on a contig shorter than 10 kb, every AMRFinderPlus gene; every CRISPR feature. Fields kept:
``type, sequence, start, stop, strand, locus, gene, product, expert`` (``type``/``source`` only),
``pseudo``; CRISPR features keep ``id`` instead of ``locus``. Sequences keep ``id, type,
topology, length`` (no ``nt``: the tests pass the replicons).
"""

from __future__ import annotations

import json
import os
from pathlib import Path

from giap.annotate.antismash import parse_antismash_json
from giap.annotate.genomad import parse_genomad
from giap.model.coords import Replicon

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ["GIAP_ROOT"]) / "recordings"
SAMPLES = ["k12_mg1655", "nissle1917"]
GENE_TYPES = {"cds", "sorf"}


def _iv(f: dict) -> tuple[str, int, int]:
    return f["sequence"], f["start"] - 1, f["stop"]


def _trim(f: dict) -> dict:
    keep = {k: f[k] for k in ("type", "sequence", "start", "stop", "strand") if k in f}
    if f["type"] in GENE_TYPES:
        keep["locus"] = f["locus"]
    else:
        keep["id"] = f.get("id") or f.get("locus")
    for k in ("gene", "product", "pseudo"):
        if f.get(k) is not None:
            keep[k] = f[k]
    if f.get("expert"):
        keep["expert"] = [
            {k: e[k] for k in ("type", "source") if e.get(k) is not None} for e in f["expert"]
        ]
    return keep


def main() -> None:
    for s in SAMPLES:
        data = json.loads((ROOT / "bakta_full" / s / "bakta" / f"{s}.json").read_text())
        seqs = data["sequences"]
        reps = {q["id"]: Replicon(q["id"], q["length"], q.get("topology", "linear")) for q in seqs}
        gm = parse_genomad(ROOT / "genomad" / s / "genomad", "genome", reps)
        bgc = parse_antismash_json(ROOT / "antismash" / s / "antismash" / "bakta.json", reps)
        spans = [
            (p.interval.contig, p.interval.start, p.interval.end)
            for p in gm.proviruses
            if p.interval is not None
        ]
        spans += [(b.contig, b.interval.start, b.interval.end) for b in bgc]
        crisprs = [_iv(f) for f in data["features"] if f["type"] == "crispr"]
        genes = [f for f in data["features"] if f["type"] in GENE_TYPES]
        keep: set[str] = set()
        for c, a, b in spans:
            on = sorted((g for g in genes if g["sequence"] == c), key=lambda g: g["start"])
            inside = [i for i, g in enumerate(on) if _iv(g)[1] < b and _iv(g)[2] > a]
            if inside:
                lo, hi = max(inside[0] - 1, 0), min(inside[-1] + 1, len(on) - 1)
                keep.update(g["locus"] for g in on[lo : hi + 1])
        for c, a, b in crisprs:
            keep.update(
                g["locus"]
                for g in genes
                if g["sequence"] == c and _iv(g)[1] < b + 3000 and _iv(g)[2] > a - 3000
            )
        small = {q["id"] for q in seqs if q["length"] < 10_000}
        keep.update(g["locus"] for g in genes if g["sequence"] in small)
        keep.update(
            g["locus"]
            for g in genes
            if any(e.get("type") == "amrfinder" for e in g.get("expert") or [])
        )
        feats = [
            _trim(f)
            for f in data["features"]
            if f["type"] == "crispr" or (f["type"] in GENE_TYPES and f["locus"] in keep)
        ]
        out = {
            "version": data.get("version"),
            "sequences": [{k: q[k] for k in ("id", "type", "topology", "length")} for q in seqs],
            "features": feats,
        }
        (HERE / f"{s}.bakta_excerpt.json").write_text(json.dumps(out, indent=0) + "\n")
        print(s, len(feats))


if __name__ == "__main__":
    main()
