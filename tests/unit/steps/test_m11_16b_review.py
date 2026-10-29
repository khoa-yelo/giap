"""#476 review (M11-16b): three tests for surviving mutants.

- M07 (``table=_table(ctx)`` dropped from ``annotate_peptide_hits`` in the ``loci`` step): a
  peptide upstream of a CDS across an in-frame TGA, through the real ``loci`` step, on both
  strands. At 4 it is ``nterm_extension``; at 11 the walk stops at the TGA (``intergenic``).
- M12: a format-3 index whose meta has ``translation_table: true`` is an unsupported format, not
  a table mismatch (``bool`` is not a table).
- M14: ``build_sixframe`` on an empty genome still refuses an unsupported table.

The N-terminal worlds are ``_plus_world`` / ``_minus_world`` of ``test_m11_16b_stops`` with the
hit window rewritten to a unique peptide, so the exact route places it at one locus only.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from helpers.readthrough_toys import CODON, rc
from test_m11_16b_mapping_steps import Chain
from test_m11_16b_stops import _minus_world, _plus_world

from gippa.mapping import sixframe as sf
from gippa.mapping.decoys import build_decoys
from gippa.model import Replicon
from gippa.rng import for_step
from gippa.steps import load_step_modules

PEP = "KDEFHPQRSV"  # 10 aa, 30 nt
PEP_NT = "".join(CODON[a] for a in PEP)


@pytest.fixture(autouse=True)
def _steps() -> None:
    load_step_modules()


def _with_hit(world: tuple[dict[str, str], Any], strand: str) -> tuple[dict[str, str], Any]:
    genome, gm = world
    nt = genome["n"]
    # the hit window: [0, 30) + or [60, 90) -
    nt = PEP_NT + nt[30:] if strand == "+" else nt[:60] + rc(PEP_NT)
    assert len(nt) == 90 and nt[45:48] in ("TGA", "TCA")
    return {"n": nt}, gm


@pytest.mark.parametrize(
    ("strand", "table", "want"),
    [
        ("+", 4, "nterm_extension"),
        ("+", 11, "intergenic"),
        ("-", 4, "nterm_extension"),
        ("-", 11, "intergenic"),
    ],
)
def test_m07_loci_step_passes_the_table_to_the_walk(
    tmp_path: Path, strand: str, table: int, want: str
) -> None:
    world = _plus_world("TGA") if strand == "+" else _minus_world("TCA")
    genome, gm = _with_hit(world, strand)
    ch = Chain(tmp_path, genome, gm, [PEP], table).index()
    hx = ch.exact()
    rows = hx[hx["peptide"] == PEP]
    span = (0, 30) if strand == "+" else (60, 90)
    got = sorted(zip(rows["start"], rows["end"], rows["strand"], strict=True))
    assert got == [(*span, strand)]
    (row,) = ch.loci().query("peptide == @PEP").to_dict("records")
    assert row["locus_class"] == want
    assert row["gene"] == ("G" if want == "nterm_extension" else "NA")
    assert ("nterm_extension_evidence" in row["flags"].split(";")) is (table == 4)


def test_m12_format_4_with_a_bool_table_is_unsupported(tmp_path: Path) -> None:
    nt = "ATGGCCTGAGCCTAA"
    index = sf.build_sixframe({"c": nt}, {"c": Replicon("c", 15, "linear")}, table=11)
    path = tmp_path / "six.npz"
    rev = build_decoys(index, "reverse", 50, for_step("t", 1))
    sf.save_index(path, index, {"reverse": rev})
    with np.load(path) as npz:
        arrays = {k: npz[k] for k in npz.files}
    meta = json.loads(str(arrays["meta"]))
    assert meta["format"] == 4  # M11-20
    meta["translation_table"] = True
    arrays["meta"] = np.array(json.dumps(meta))
    bad = tmp_path / "bad.npz"
    with bad.open("wb") as fh:
        np.savez_compressed(fh, **arrays)  # type: ignore[arg-type]
    for fn in (sf.load_index, sf.load_stop_cells, sf.load_start_cells):
        for table in (1, 11, 4):
            with pytest.raises(ValueError, match="unsupported sixframe_index format") as ei:
                fn(bad, table=table)
            assert "built for translation table" not in str(ei.value)


def test_m14_empty_genome_still_refuses_an_unsupported_table() -> None:
    with pytest.raises(ValueError, match="unsupported translation table 12"):
        sf.build_sixframe({}, {}, table=12)
