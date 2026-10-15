"""Tool inputs for scripts/slurm/record_tools.sbatch, built from a recorded Bakta JSON.

Stand-ins for the GIAP steps that do not exist yet (``transterm_input`` M2-15, ``ostir_input``
M2-16, TU upstream regions M2-14), following the spec definitions so the recorded tool outputs have
the real formats. Stdlib only.

    record_tools_inputs.py coords  BAKTA_JSON OUT.coords            # TransTermHP genes.coords
    record_tools_inputs.py ostir   BAKTA_JSON OUT.csv               # OSTIR input (spec §7.1.6)
    record_tools_inputs.py regions BAKTA_JSON OUTDIR                # one FASTA per TU (spec §7.1.3)

TUs use the distance-model fallback with its default parameters (spec §7.1.3: p = 1 / (1 +
exp((d - 50) / 20)) >= 0.5, i.e. same-strand neighbours with intergenic distance d <= 50 bp), since
no UniOP recording is wired in here.
"""

from __future__ import annotations

import collections
import json
import sys
from pathlib import Path
from typing import Any

GENE_TYPES = {"cds", "sorf"}
ANTI_SD_DEFAULT = "ACCTCCTTA"
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def load(path: Path) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    data = json.loads(path.read_text())
    seqs = {s["id"]: s for s in data["sequences"]}
    genes = [f for f in data["features"] if f["type"] in GENE_TYPES]
    genes.sort(key=lambda f: (f.get("sequence") or f.get("contig"), f["start"], f["stop"]))
    return seqs, genes


def contig(f: dict[str, Any]) -> str:
    return str(f.get("sequence") or f.get("contig"))


def anti_sd(path: Path) -> str:
    """Most common 3'-terminal 9-mer of the 16S rRNAs (spec §7.1.2); E. coli default otherwise."""
    feats = json.loads(path.read_text())["features"]
    rrna16 = [f for f in feats if f["type"] == "rRNA" and "16S" in f.get("product", "")]
    tails = [f["nt"][-9:].upper() for f in rrna16]
    if not tails:
        return ANTI_SD_DEFAULT
    best = collections.Counter(tails).most_common(1)[0][0]
    return best if "CCTCC" in best else ANTI_SD_DEFAULT


def window(seq: str, circular: bool, a: int, b: int) -> tuple[str, int]:
    """genome[a, b), wrapped (circular) or clipped (linear); also the bases clipped at a."""
    n = len(seq)
    if circular:
        return "".join(seq[i % n] for i in range(a, b)), 0
    lo, hi = max(a, 0), min(b, n)
    return seq[lo:hi], lo - a


def write_coords(js: Path, out: Path) -> None:
    _, genes = load(js)
    with out.open("w") as fh:
        for g in genes:
            s, e = (g["start"], g["stop"]) if g["strand"] == "+" else (g["stop"], g["start"])
            fh.write(f"{g['locus']} {s} {e} {contig(g)}\n")


def write_ostir(js: Path, out: Path) -> None:
    seqs, genes = load(js)
    asd = anti_sd(js)
    with out.open("w") as fh:
        fh.write("id,seq,start,end,anti-Shine-Dalgarno\n")
        for g in genes:
            rec = seqs[contig(g)]
            circ = rec["topology"] == "circular"
            if g["strand"] == "+":
                s0 = g["start"] - 1
                seq, clip = window(rec["nt"], circ, s0 - 50, s0 + 50)
            else:
                e0 = g["stop"]
                fwd, _ = window(rec["nt"], circ, e0 - 50, e0 + 50)
                clip = 0 if circ else max(e0 + 50 - len(rec["nt"]), 0)
                seq = revcomp(fwd)
            pos = 51 - clip
            fh.write(f"{g['locus']},{seq.upper()},{pos},{pos},{asd}\n")


def write_regions(js: Path, outdir: Path, up: int = 300, down: int = 50, d_max: int = 50) -> None:
    seqs, genes = load(js)
    outdir.mkdir(parents=True, exist_ok=True)
    by_contig: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for g in genes:
        by_contig[contig(g)].append(g)
    n_tu = 0
    for cid, gs in by_contig.items():
        nt = seqs[cid]["nt"]
        tus: list[list[dict[str, Any]]] = []
        for g in gs:
            prev = tus[-1][-1] if tus else None
            linked = prev is not None and prev["strand"] == g["strand"]
            if linked and prev is not None and g["start"] - prev["stop"] - 1 <= d_max:
                tus[-1].append(g)
            else:
                tus.append([g])
        for tu in tus:
            n_tu += 1
            tu_id = f"TU_{n_tu:05d}"
            if tu[0]["strand"] == "+":
                s0 = tu[0]["start"] - 1
                ends = [o["stop"] for o in gs if o["stop"] <= s0]
                a = max(s0 - up, max(ends, default=0), 0)
                region = nt[a : min(s0 + down, len(nt))]
            else:
                e0 = tu[-1]["stop"]
                starts = [o["start"] - 1 for o in gs if o["start"] - 1 >= e0]
                b = min(e0 + up, min(starts, default=len(nt)), len(nt))
                region = revcomp(nt[max(e0 - down, 0) : b])
            (outdir / f"{tu_id}.fa").write_text(f">{tu_id}\n{region.upper()}\n")


def main(argv: list[str]) -> None:
    cmd, js, out = argv[1], Path(argv[2]), Path(argv[3])
    {"coords": write_coords, "ostir": write_ostir, "regions": write_regions}[cmd](js, out)


if __name__ == "__main__":
    main(sys.argv)
