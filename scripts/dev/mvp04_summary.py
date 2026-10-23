"""MVP-04 evidence-run summary (docs/specs/MVP-04-ecoli-evidence-run.md "The run", acceptance
gate 7; lead decisions 11 and 15).

Reads only published GIAP outputs (``tables/genes.tsv``, ``mapping/peptide_hits.tsv``) and the
giapsim inputs (``peptides.tsv``, ``truth/peptides.tsv``) and writes
``<outdir>/evidence_summary.json``. So far it carries gate 7 under ``gate7``:

* a gene is **true** when it is the ``locus_tag`` of a truth row with ``is_wrong = false``
  (decision 11; a protein-entry null's truth row carries its locus tag but ``is_wrong = true``);
* ``realised_gene_fdp`` = detected genes that are not true / detected genes;
* ``false_detected_by_channel`` counts the false-detected genes with at least one
  ``peptide_hits`` row (``gene`` = the gene) of a null peptide of the channel; a gene may count
  in several (decision 15): ``protein_null`` (truth ``origin = undetected_product``),
  ``nterm_extension_sixframe_null`` (``origin = sixframe_noncoding`` whose target rows in the
  input carry a ``GIAPORF_nterm_extension_*`` entry), ``tolerant_sixframe_null`` (any other
  six-frame null whose hit on the gene has ``route = tolerant``), and ``other`` (none of these).

    python scripts/dev/mvp04_summary.py <outdir> <inputs dir>
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

PROTEIN_NULL_ORIGIN = "undetected_product"  # MVP-04 lead decision 9
SIXFRAME_NULL_ORIGIN = "sixframe_noncoding"
NTERM_PREFIX = "GIAPORF_nterm_extension_"
CHANNELS = ("protein_null", "nterm_extension_sixframe_null", "tolerant_sixframe_null", "other")
SUMMARY_NAME = "evidence_summary.json"


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _ids(cell: str | None) -> list[str]:
    return [x for x in (cell or "").split(";") if x and x != "NA"]


def _true(cell: str | None) -> bool:
    return (cell or "").strip().lower() == "true"


def gate7_from_rows(
    genes: Iterable[Mapping[str, str]],
    hits: Iterable[Mapping[str, str]],
    peptides: Iterable[Mapping[str, str]],
    truth: Iterable[Mapping[str, str]],
) -> dict[str, Any]:
    """Gate 7 from the four tables' rows (module doc)."""
    detected = {g["locus_tag"] for g in genes if _true(g.get("detected"))}
    true_genes: set[str] = set()
    origin: dict[str, str] = {}
    for t in truth:
        if _true(t.get("is_wrong")):
            if t.get("origin") in (PROTEIN_NULL_ORIGIN, SIXFRAME_NULL_ORIGIN):
                origin[t["peptide"]] = str(t["origin"])
        else:
            true_genes.update(_ids(t.get("locus_tag")))
    nterm = {
        p["peptide"]
        for p in peptides
        if not _true(p.get("is_decoy"))
        and any(e.startswith(NTERM_PREFIX) for e in _ids(p.get("protein_ids")))
    }
    false_genes = detected - true_genes
    by: dict[str, set[str]] = {k: set() for k in CHANNELS}
    for h in hits:
        gene, pep = h.get("gene", ""), h["peptide"]
        if gene not in false_genes or pep not in origin:
            continue
        if origin[pep] == PROTEIN_NULL_ORIGIN:
            by["protein_null"].add(gene)
        elif pep in nterm:
            by["nterm_extension_sixframe_null"].add(gene)
        elif h.get("route") == "tolerant":
            by["tolerant_sixframe_null"].add(gene)
    by["other"] = false_genes - set().union(*(by[k] for k in CHANNELS[:-1]))
    n_det, n_false = len(detected), len(false_genes)
    return {
        "n_detected": n_det,
        "n_false_detected": n_false,
        "realised_gene_fdp": n_false / n_det if n_det else float("nan"),
        "false_detected_by_channel": {k: len(by[k]) for k in CHANNELS},
    }


def gate7(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Gate 7 of the run in ``outdir`` against the giapsim inputs in ``inputs``."""
    outdir, inputs = Path(outdir), Path(inputs)
    return gate7_from_rows(
        _rows(outdir / "tables" / "genes.tsv"),
        _rows(outdir / "mapping" / "peptide_hits.tsv"),
        _rows(inputs / "peptides.tsv"),
        _rows(inputs / "truth" / "peptides.tsv"),
    )


def summarise(outdir: Path, inputs: Path) -> dict[str, Any]:
    """The evidence summary (gate 7 so far)."""
    return {"gate7": gate7(outdir, inputs)}


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("outdir", type=Path)
    ap.add_argument("inputs", type=Path)
    a = ap.parse_args(argv)
    path = a.outdir / SUMMARY_NAME
    path.write_text(json.dumps(summarise(a.outdir, a.inputs), indent=2, sort_keys=True) + "\n")
    print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
