"""The ``sec_site`` / ``pyl_site`` flags and the ``recoded_codon`` GFF3 attribute (M7-16c; spec
rule 19, rule 21 (f), rule 26, examples 20, 27 and 30 (mapping part); test plan C8).

``annotate_peptide_hits(..., consumed=...)``: a row consuming a site of kind ``recoded`` gets
``sec_site`` when the site's codon is TGA and ``pyl_site`` when it is TAG, **instead of**
``readthrough``; a row consuming a terminal and a recoded site gets both flags. The flag depends on
the consumed site, never on the query letter (a ``W`` over the recoded TGA still gets
``sec_site``). ``mappings.gff3`` carries ``recoded_codon=<p1>[,<p2>]`` (1-based ``codon_start +
1``, ascending) after ``flags``, as ``readthrough_codon``; the work column never reaches
``peptide_hits.tsv``.

Toy world: the 16b contig ``c1`` with ``G1`` annotated as one CDS ``[0, 93)`` whose protein has
``U`` at residue 18 (recoded TGA ``[54, 57)``, terminal TAA ``[90, 93)``).
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.readthrough_toys import back_translate, gene, model
from test_mapping_outputs_readthrough import annotate, attrs, feature_of, hit, tsv, write

from giap.mapping.outputs import PEPTIDE_HITS_COLUMNS
from giap.mapping.readthrough import StopSite
from giap.model import Replicon

N = 93
SEC = StopSite("c1", "+", 54, ("G1",), "TGA", "recoded")
PYL = StopSite("c1", "+", 54, ("G1",), "TAG", "recoded")
TERM = StopSite("c1", "+", 90, ("G1",), "TAA", "terminal")
EX20 = "WQLNEFKUVPAGHTL"  # example 20's shape: U over the recoded TGA, aa [11, 26)
EX30 = "WQLNEFKWVPAGHTL"  # example 30's shape: W over the same TGA


def recoded_world(codon: str = "TGA", letter: str = "U") -> Any:
    nt = back_translate("MSDKPLTEAGRWQLNEFK*VPAGHTLDMRK*", {18: codon, 30: "TAA"})
    reps = {"c1": Replicon("c1", N, "linear")}
    aa = "MSDKPLTEAGRWQLNEFK" + letter + "VPAGHTLDMRK"
    return {"c1": nt}, reps, model(reps, [replace(gene("G1", "c1", 0, N, "+"), aa=aa)])


def flags_of(ann: pd.DataFrame, peptide: str) -> set[str]:
    v = ann.set_index("peptide").loc[peptide, "flags"]
    return set(str(v).split(";")) - {""}


# ---------------------------------------------------------------- the flags


@pytest.mark.parametrize(
    ("relation", "want"), [("same", "sec_site"), ("related", "sec_site;strain_difference")]
)
def test_example_20_sec_site_and_no_readthrough(relation: str, want: str) -> None:
    ann = annotate([hit(EX20, "c1", 33, 78)], [(SEC,)], relation, world=recoded_world())
    (r,) = ann.to_dict("records")
    assert r["flags"] == want
    assert r["locus_class"] == "cds_inframe" and r["gene"] == "G1"


def test_pyl_site_for_a_recoded_tag() -> None:
    world = recoded_world("TAG", "O")
    ann = annotate([hit("WQLNEFKOVPAGHTL", "c1", 33, 78)], [(PYL,)], world=world)
    assert ann["flags"].iloc[0] == "pyl_site"


def test_terminal_and_recoded_site_in_one_row_carry_both_flags() -> None:
    pep = "AGRWQLNEFKUVPAGHTLDMRKQ"[:20]
    row = hit(pep, "c1", 33, 93, diff_types="substitution;substitution")
    ann = annotate([row], [(SEC, TERM)], world=recoded_world())
    assert flags_of(ann, pep) == {"readthrough", "sec_site"}
    assert ann["flags"].iloc[0] == "readthrough;sec_site"


def test_example_30_w_over_the_recoded_tga_still_gets_sec_site() -> None:
    ann = annotate([hit(EX30, "c1", 33, 78)], [(SEC,)], world=recoded_world())
    assert ann["flags"].iloc[0] == "sec_site"


def test_flags_through_consumed_sites_of_the_recoded_site_set() -> None:
    """The ``loci`` path: ``consumed_sites`` over ``annotated_stops(..., recoded=True)``."""
    from giap.mapping.readthrough import annotated_stops, consumed_sites

    genome, reps, gm = recoded_world()
    sites, _ = annotated_stops(gm, genome, recoded=True)
    assert sites == (SEC, TERM)
    rows = [
        hit(EX20, "c1", 33, 78),
        hit(EX30, "c1", 33, 78),
        hit("MSDKPLTEAG", "c1", 0, 30, route="exact"),
    ]
    df = pd.DataFrame(rows)
    got = consumed_sites(df, sites, {"c1": N}, {"c1": False})
    assert got == [(SEC,), (SEC,), ()]
    ann = annotate(rows, got, world=(genome, reps, gm))
    assert flags_of(ann, EX20) == {"sec_site"} and flags_of(ann, EX30) == {"sec_site"}
    assert flags_of(ann, "MSDKPLTEAG") == set()


def test_example_27_merged_site_row_gets_sec_site_not_readthrough() -> None:
    """Rule 21 (f) / Q22: G5's terminal TGA is G1's recoded codon; the merged site is
    ``recoded``, so a row over ``[15, 18)`` gets ``sec_site`` and never ``readthrough``."""
    from test_readthrough_recoded import ex27

    from giap.mapping.readthrough import annotated_stops, consumed_sites

    genome, gm = ex27()
    sites, _ = annotated_stops(gm, genome, recoded=True)
    pep = "AMAAAUAAAA"
    rows = [hit(pep, "c", 0, 30)]
    got = consumed_sites(pd.DataFrame(rows), sites, {"c": 60}, {"c": False})
    assert [[(s.codon_start, s.kind) for s in t] for t in got] == [[(15, "recoded")]]
    ann = annotate(rows, got, world=(genome, dict(gm.replicons), gm))
    assert ann["flags"].iloc[0] == "sec_site"


# ---------------------------------------------------------------- the writer


def test_example_20_gff3_recoded_codon_after_flags(tmp_path: Path) -> None:
    world = recoded_world()
    ann = annotate([hit(EX20, "c1", 33, 78)], [(SEC,)], world=world)
    out = write(tmp_path / "m", ann, world[1])
    kv = attrs(feature_of(out, EX20)[8])
    keys = [k for k, _ in kv]
    assert dict(kv)["flags"] == "sec_site"
    assert dict(kv)["recoded_codon"] == "55"
    assert keys.index("recoded_codon") == keys.index("flags") + 1
    assert "readthrough_codon" not in dict(kv)


def test_both_attributes_on_a_two_site_row(tmp_path: Path) -> None:
    world = recoded_world()
    pep = "AGRWQLNEFKUVPAGHTLDM"
    row = hit(pep, "c1", 33, 93, diff_types="substitution;substitution")
    plain = write(tmp_path / "a", annotate([row], None, world=world), world[1])
    out = write(tmp_path / "b", annotate([row], [(SEC, TERM)], world=world), world[1])
    kv = dict(attrs(feature_of(out, pep)[8]))
    assert kv["readthrough_codon"] == "91" and kv["recoded_codon"] == "55"
    assert set(kv["flags"].split(",")) == {"readthrough", "sec_site"}
    ka = [k for k, _ in attrs(feature_of(plain, pep)[8])]
    kb = [k for k, _ in attrs(feature_of(out, pep)[8])]
    assert kb[: len(ka)] == ka  # the existing order is unchanged
    assert set(kb[len(ka) :]) == {"readthrough_codon", "recoded_codon"}


def test_two_recoded_cells_list_both_positions_ascending(tmp_path: Path) -> None:
    a = StopSite("c1", "+", 15, ("A",), "TGA", "recoded")
    b = StopSite("c1", "+", 30, ("A",), "TGA", "recoded")
    pep = "AAAAAUAAAAUAAAA"
    ann = annotate([hit(pep, "c1", 0, 45, diff_types="substitution;substitution")], [(a, b)])
    assert ann["flags"].iloc[0] == "sec_site"
    out = write(tmp_path / "m", ann)
    assert dict(attrs(feature_of(out, pep)[8]))["recoded_codon"] == "16,31"


def test_tsv_columns_unchanged_and_recoded_work_column_not_written(tmp_path: Path) -> None:
    world = recoded_world()
    ann = annotate([hit(EX20, "c1", 33, 78)], [(SEC,)], world=world)
    out = write(tmp_path / "m", ann, world[1])
    header = (out / "peptide_hits.tsv").read_text().splitlines()[0].split("\t")
    assert tuple(header) == PEPTIDE_HITS_COLUMNS
    assert "recoded_codon" not in (out / "peptide_hits.tsv").read_text()
    assert tsv(out).set_index("peptide").loc[EX20, "flags"] == "sec_site"


def test_bed_probed_and_unmapped_unchanged_by_the_recoded_flag(tmp_path: Path) -> None:
    world = recoded_world()
    rows = [hit(EX20, "c1", 33, 78)]
    a = write(tmp_path / "a", annotate(rows, None, world=world), world[1])
    b = write(tmp_path / "b", annotate(rows, [(SEC,)], world=world), world[1])
    for name in ("mappings.bed", "mappings.probed", "unmapped.tsv", "protein_hits.tsv"):
        assert (a / name).read_bytes() == (b / name).read_bytes(), name


def test_frameshift_rows_never_carry_sec_site(tmp_path: Path) -> None:
    world = recoded_world()
    fs = hit("MSDKPLTEAGRWQLNE", "c1", 0, 49, diff_types="frameshift", blocks="0-24;25-49")
    ann = annotate([hit(EX20, "c1", 33, 78), fs], [(SEC,), ()], world=world)
    assert "sec_site" not in flags_of(ann, "MSDKPLTEAGRWQLNE")
    out = write(tmp_path / "m", ann, world[1])
    for r in (ln.split("\t") for ln in (out / "mappings.gff3").read_text().splitlines()):
        if len(r) > 8 and r[2] == "match_part":
            assert "recoded_codon=" not in r[8]
