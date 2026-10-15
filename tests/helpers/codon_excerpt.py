"""K-12 codon-usage excerpt (M2-17): build it from the recorded Bakta JSON, turn it back into a
small self-contained Bakta JSON plus genome FASTA for the ``parse_bakta`` → ``codon_usage`` steps.

The excerpt (``tests/data/codon_k12_excerpt.json.gz``, provenance next to it) keeps, from the
recorded K-12 MG1655 full-DB Bakta 1.12.2 JSON: every ribosomal-protein gene (so the CAI/Fop
reference set equals the whole genome's), every tRNA anticodon (so the tAI weights equal the
whole genome's), the EF-Tu genes, every pseudogene, every sORF shorter than 46 nt, a few named
genes with alternative starts, and every 50th remaining CDS. ``to_bakta`` places the genes on one
synthetic linear contig (plus strand, random 60-nt spacers from a fixed seed) so the step can run
without the 4.6 Mb genome; gene nucleotide sequences are the recorded ones.

Regenerate (only when the recording changes)::

    REC=$GIAP_ROOT/recordings/bakta_full/k12_mg1655/bakta/k12_mg1655.json
    python -P tests/helpers/codon_excerpt.py $REC tests/data/codon_k12_excerpt.json.gz
"""

from __future__ import annotations

import gzip
import json
import random
import sys
from pathlib import Path
from typing import Any

EXCERPT = Path(__file__).resolve().parents[1] / "data" / "codon_k12_excerpt.json.gz"
NAMED = ("rpsT", "nhaA", "carA", "lacZ", "thrA")
CONTIG = "k12x"


def build_excerpt(bakta_json: Path) -> dict[str, Any]:
    """Recorded Bakta JSON → excerpt dict (genes in genome order, tRNA anticodons)."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
    from giap.annotate.bakta import parse_bakta_json
    from giap.annotate.bakta_layers import is_ribosomal_protein

    ann = parse_bakta_json(bakta_json)
    keep = []
    others = 0
    for g in ann.model.genes:
        special = (
            is_ribosomal_protein(g)
            or "elongation factor Tu" in g.product
            or "EF-Tu" in g.product
            or g.feature_type == "pseudogene"
            or (g.feature_type == "sorf" and len(g.nt or "") < 46)
            or g.gene in NAMED
        )
        if not special:
            others += 1
            if others % 50:
                continue
        keep.append(
            {
                "locus": g.locus_tag,
                "type": g.feature_type,
                "gene": g.gene,
                "product": g.product,
                "partial": g.partial,
                "start_type": g.start_codon,
                "nt": g.nt,
            }
        )
    trnas = [
        {"anti_codon": f.anticodon, "amino_acid": f.amino_acid}
        for f in ann.model.features
        if f.type == "trna" and f.anticodon
    ]
    return {
        "source": (
            f"{bakta_json.name} (Bakta {ann.bakta_version}, DB {ann.db_version} {ann.db_type})"
        ),
        "genes": keep,
        "trnas": trnas,
    }


def load_excerpt(path: Path = EXCERPT) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        data: dict[str, Any] = json.load(fh)
    return data


def to_bakta(excerpt: dict[str, Any], seed: int = 17) -> tuple[dict[str, Any], str]:
    """Excerpt → (Bakta-like JSON dict, genome FASTA text) on one synthetic linear contig."""
    rnd = random.Random(seed)

    def spacer(n: int = 60) -> str:
        return "".join(rnd.choice("ACGT") for _ in range(n))

    parts: list[str] = [spacer()]
    pos = len(parts[0])
    feats: list[dict[str, Any]] = []
    for g in excerpt["genes"]:
        nt = g["nt"]
        f: dict[str, Any] = {
            "type": "sorf" if g["type"] == "sorf" else "cds",
            "sequence": CONTIG,
            "start": pos + 1,
            "stop": pos + len(nt),
            "strand": "+",
            "locus": g["locus"],
            "gene": g["gene"],
            "product": g["product"],
            "aa": "M",
            "nt": nt,
            "start_type": g["start_type"],
        }
        if g["type"] == "pseudogene":
            f["pseudo"] = True
        if g["partial"]:
            f["truncated"] = "5'"
        feats.append(f)
        parts += [nt, spacer()]
        pos += len(nt) + 60
    for k, t in enumerate(excerpt["trnas"], 1):
        feats.append(
            {
                "type": "tRNA",
                "sequence": CONTIG,
                "start": pos + 1,
                "stop": pos + 76,
                "strand": "+",
                "id": f"k12x_t{k}",
                "product": f"tRNA-{t['amino_acid']}",
                "amino_acid": t["amino_acid"],
                "anti_codon": t["anti_codon"].lower(),
            }
        )
        parts += [spacer(76), spacer()]
        pos += 136
    genome = "".join(parts)
    assert len(genome) == pos
    data = {
        "sequences": [
            {"id": CONTIG, "type": "chromosome", "topology": "linear", "length": len(genome)}
        ],
        "features": feats,
        "version": {"bakta": "1.12.2", "db": {"version": "6.0", "type": "full"}},
    }
    return data, f">{CONTIG}\n{genome}\n"


def write_inputs(directory: Path, excerpt: dict[str, Any] | None = None) -> tuple[Path, Path]:
    """Write ``bakta.json`` and ``genome.norm.fa`` for the excerpt into ``directory``."""
    data, fasta = to_bakta(excerpt or load_excerpt())
    directory.mkdir(parents=True, exist_ok=True)
    js, fa = directory / "bakta.json", directory / "genome.norm.fa"
    js.write_text(json.dumps(data), encoding="utf-8")
    fa.write_text(fasta, encoding="utf-8")
    return js, fa


if __name__ == "__main__":
    src, dst = Path(sys.argv[1]), Path(sys.argv[2])
    ex = build_excerpt(src)
    with dst.open("wb") as raw, gzip.GzipFile("", "wb", fileobj=raw, mtime=0) as fh:
        fh.write((json.dumps(ex, indent=0, sort_keys=True) + "\n").encode())
    print(f"{dst}: {len(ex['genes'])} genes, {len(ex['trnas'])} tRNAs")
