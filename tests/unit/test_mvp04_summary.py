"""Gate 7 in ``scripts/dev/mvp04_summary.py`` (docs/specs/MVP-04-ecoli-evidence-run.md acceptance
7, lead decision 11; statistician review of PR #247, B3(a) and item 9).

Decision 11: a "true peptide" is a ``truth/peptides.tsv`` row with ``is_wrong = false``. A
protein-entry null's truth row carries its locus tag but has ``is_wrong = true``, so it never
makes a gene true (a summary that counted any truth row would score every null gene as true,
FDP = 0, and gate 7 would be vacuous). The summary also reports the false-detected genes by
channel.

Assumed interface (test-author reading; the implementer may add keys, not rename these):

``gate7(outdir: Path, inputs: Path) -> dict`` reads only published outputs
(``tables/genes.tsv``: ``locus_tag``, ``detected``; ``mapping/peptide_hits.tsv``: ``peptide``,
``route``, ``gene``) and the gippasim inputs (``peptides.tsv``: ``peptide``, ``protein_ids``,
``is_decoy``; ``truth/peptides.tsv``: ``peptide``, ``is_wrong``, ``origin``, ``locus_tag``) and
returns

* ``n_detected``: genes with ``detected = true``;
* ``n_false_detected``: detected genes that are no ``locus_tag`` of a truth row with
  ``is_wrong = false``;
* ``realised_gene_fdp`` = ``n_false_detected / n_detected``;
* ``false_detected_by_channel``: false-detected genes with at least one ``peptide_hits`` row
  (``gene`` = the gene) of a null peptide of that channel; a gene may count in several:

  - ``protein_null``: truth ``origin = undetected_product`` (decision 9);
  - ``nterm_extension_sixframe_null``: truth ``origin = sixframe_noncoding`` whose target rows in
    the input ``peptides.tsv`` carry a ``GIPPAORF_nterm_extension_*`` entry (any route);
  - ``tolerant_sixframe_null``: any other six-frame null whose hit on the gene has
    ``route = tolerant`` (the one-mismatch route under ``--genome-relation related``);
  - ``other``: false-detected genes with none of the three.

``evidence_summary.json`` carries this dict under ``gate7`` (not tested here: the full summary
needs every gate's outputs).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "dev" / "mvp04_summary.py"
NT = "GIPPAORF_nterm_extension_U00096.3_4901_5300_p"
AS = "GIPPAORF_antisense_U00096.3_9001_9150_m"
CHANNELS = ("protein_null", "nterm_extension_sixframe_null", "tolerant_sixframe_null", "other")


def _summary() -> ModuleType:
    spec = importlib.util.spec_from_file_location("mvp04_summary", SCRIPT)
    assert spec is not None and spec.loader is not None, SCRIPT
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("mvp04_summary", mod)
    spec.loader.exec_module(mod)
    return mod


def _tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join("\t".join(r) for r in [header, *rows]) + "\n")


# gene: (detected, [(peptide, route)] hits)
GENES: dict[str, tuple[bool, list[tuple[str, str]]]] = {
    "G1": (True, [("TRUEPEPK", "exact")]),  # true
    "G2": (True, [("PROTNULLK", "exact")]),  # protein null only: false (decision 11)
    "G3": (True, [("NTERMNULLK", "exact")]),  # nterm six-frame null
    "G4": (True, [("ANTINULLK", "tolerant")]),  # tolerant six-frame null
    "G5": (True, [("PROTNULLTWOK", "exact"), ("ANTINULLTWOK", "tolerant")]),  # both channels
    "G6": (True, [("TRUEPEPK", "exact")]),  # false, evidence from G1's true peptide: other
    "G7": (False, [("PROTNULLSEVK", "exact")]),  # not detected
    "G8": (True, [("TRUEEIGHTK", "exact"), ("PROTNULLEIGHTK", "exact")]),  # true + a null
    "G9": (True, [("ANTIEXACTK", "exact")]),  # six-frame null, exact, not nterm: other
}
TRUTH = [  # peptide, is_wrong, origin, locus_tag
    ("TRUEPEPK", "false", "cds", "G1"),
    ("TRUEEIGHTK", "false", "cds", "G8"),
    ("PROTNULLK", "true", "undetected_product", "G2"),
    ("PROTNULLTWOK", "true", "undetected_product", "G5"),
    ("PROTNULLSEVK", "true", "undetected_product", "G7"),
    ("PROTNULLEIGHTK", "true", "undetected_product", "G8"),
    ("NTERMNULLK", "true", "sixframe_noncoding", "NA"),
    ("ANTINULLK", "true", "sixframe_noncoding", "NA"),
    ("ANTINULLTWOK", "true", "sixframe_noncoding", "NA"),
    ("ANTIEXACTK", "true", "sixframe_noncoding", "NA"),
]
ENTRY = {"TRUEPEPK": "G1", "TRUEEIGHTK": "G8", "PROTNULLK": "G2", "PROTNULLTWOK": "G5",
         "PROTNULLSEVK": "G7", "PROTNULLEIGHTK": "G8", "NTERMNULLK": NT, "ANTINULLK": AS,
         "ANTINULLTWOK": AS, "ANTIEXACTK": AS}  # fmt: skip


@pytest.fixture
def dirs(tmp_path: Path) -> tuple[Path, Path]:
    out, inp = tmp_path / "run", tmp_path / "inputs"
    _tsv(out / "tables" / "genes.tsv", ["locus_tag", "detected", "detection_status"],
         [[g, "true" if d else "false", "detected" if d else "not_detected"]
          for g, (d, _) in GENES.items()])  # fmt: skip
    hits = [[p, route, "1" if route == "tolerant" else "0", "cds_inframe", "unique", g]
            for g, (_, hs) in GENES.items() for p, route in hs]  # fmt: skip
    _tsv(out / "mapping" / "peptide_hits.tsv",
         ["peptide", "route", "n_diff", "locus_class", "uniqueness", "gene"], hits)  # fmt: skip
    _tsv(inp / "truth" / "peptides.tsv",
         ["peptide", "is_wrong", "origin", "locus_tag", "contig", "start", "end", "strand",
          "n_loci_genome", "alt_start_nterm"],
         [[p, w, o, t, "NA", "NA", "NA", "NA", "1", "false"] for p, w, o, t in TRUTH])  # fmt: skip
    peps = [[p, e, "0.001", "0.01", "5.0", "false", "1000", "run1"] for p, e in ENTRY.items()]
    peps.append(["KLLUNMRETN", "rev_" + NT, "0.002", "0.02", "4.0", "true", "NA", "run2"])
    _tsv(inp / "peptides.tsv",
         ["peptide", "protein_ids", "q_value", "pep", "score", "is_decoy", "intensity", "run"],
         peps)  # fmt: skip
    return out, inp


def _gate7(dirs: tuple[Path, Path]) -> dict[str, Any]:
    res: dict[str, Any] = _summary().gate7(*dirs)
    return res


def test_summary_script_exists() -> None:
    assert SCRIPT.is_file(), SCRIPT


def test_true_gene_needs_a_truth_row_with_is_wrong_false(dirs: tuple[Path, Path]) -> None:
    """G2 has a truth row with its locus tag (a protein null) and is still false."""
    res = _gate7(dirs)
    assert res["n_detected"] == 8
    assert res["n_false_detected"] == 6  # G2, G3, G4, G5, G6, G9
    assert res["realised_gene_fdp"] == pytest.approx(6 / 8)


def test_false_detected_genes_by_channel(dirs: tuple[Path, Path]) -> None:
    by = _gate7(dirs)["false_detected_by_channel"]
    assert set(by) >= set(CHANNELS), by
    assert {k: by[k] for k in CHANNELS} == {
        "protein_null": 2,  # G2, G5
        "nterm_extension_sixframe_null": 1,  # G3
        "tolerant_sixframe_null": 2,  # G4, G5
        "other": 2,  # G6 (another gene's true peptide), G9 (exact non-nterm six-frame null)
    }


def test_a_gene_with_a_true_peptide_and_a_null_is_true(dirs: tuple[Path, Path]) -> None:
    """G8: a true peptide makes the gene true whatever null evidence it also has; undetected G7
    counts nowhere."""
    out, _inp = dirs
    genes = out / "tables" / "genes.tsv"
    drop = ("G2", "G3", "G4", "G5", "G6", "G9")
    keep = [ln for ln in genes.read_text().splitlines() if not ln.startswith(drop)]
    genes.write_text("\n".join(keep) + "\n")
    res = _gate7(dirs)
    assert (res["n_detected"], res["n_false_detected"]) == (2, 0)
    assert res["realised_gene_fdp"] == 0.0
    assert all(res["false_detected_by_channel"][k] == 0 for k in CHANNELS)
