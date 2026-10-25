"""The ``readthrough`` flag and the ``readthrough_codon`` GFF3 attribute (M7-16b rule 16; spec
examples 11, 12, 15 and the edge cases "Multi-locus peptide", "Two cells in one row").

``annotate_peptide_hits(..., consumed=...)`` takes ``consumed_sites``' output (one tuple of
``StopSite`` per input row, input order). A row with at least one consumed site of kind
``terminal`` gets the flag ``readthrough`` (flags stay ``;``-joined alphabetically). The consumed
codon positions travel in a work column ``readthrough_codon`` that ``write_mapping_outputs`` never
writes to ``peptide_hits.tsv``; ``mappings.gff3`` carries ``readthrough_codon=<p1>[,<p2>]`` (the
1-based ``codon_start + 1``, ascending) right after ``flags``, on flagged peptide features only,
never on ``match_part`` children. BED, proBed and ``unmapped.tsv`` are unchanged.

The work column's own encoding is the implementer's choice; these tests read only what the
writer makes of it. Toy contig ``c1`` / gene ``G1`` of the spec (stop TGA ``[54, 57)``).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.readthrough_toys import gene, model, toy_c1

import giap.mapping.loci as _loci_mod
from giap.mapping.outputs import (
    MAPPING_FDR_COLUMNS,
    PEPTIDE_HITS_COLUMNS,
    PROTEIN_HITS_COLUMNS,
    write_mapping_outputs,
)
from giap.mapping.readthrough import StopSite
from giap.model import Replicon

G1_SITE = StopSite("c1", "+", 54, ("G1",), "TGA", "terminal")
SHA = "cd" * 32
EX1 = "WQLNEFKWVPAGHTL"


def hit(
    peptide: str,
    contig: str,
    start: int,
    end: int,
    strand: str = "+",
    *,
    route: str = "tolerant",
    diff_types: str | None = "substitution",
    blocks: str | None = None,
    flags: str | None = None,
) -> dict[str, Any]:
    """One ``_peptide_rows`` row: hit columns plus the input and placement columns."""
    n_diff = 0 if route == "exact" else len(str(diff_types).split(";"))
    return {
        "peptide": peptide,
        "query": peptide,
        "search": "target",
        "decoy_type": None,
        "contig": contig,
        "frame_id": 0,
        "aa_start": 0,
        "aa_end": (end - start) // 3,
        "start": start,
        "end": end,
        "strand": strand,
        "frame": 1,
        "route": route,
        "n_diff": n_diff,
        "diff_types": None if route == "exact" else diff_types,
        "blocks": blocks,
        "peptide_len": len(peptide),
        "source": "denovo",
        "placement_q": 0.001,
        "flags": flags,
        "best_q_input": float("nan"),
        "pep_input": float("nan"),
        "score_input": 0.9,
        "length_class": "15-19",
    }


def _peptides(peps: Sequence[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {"peptide": list(peps), "search_db": "other", "protein_ids": [[] for _ in peps]}
    )


def annotate(
    rows: Sequence[dict[str, Any]],
    consumed: Sequence[tuple[StopSite, ...]] | None,
    relation: str = "same",
    world: Any = None,
) -> pd.DataFrame:
    genome, _reps, gm = world or toy_c1()
    df = pd.DataFrame(list(rows))
    kw = {} if consumed is None else {"consumed": list(consumed)}
    return _loci_mod.annotate_peptide_hits(
        df, gm, genome, _peptides(df["peptide"].unique()), relation, **kw
    )  # type: ignore[arg-type]


def write(out: Path, ann: pd.DataFrame, reps: Mapping[str, Replicon] | None = None) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    reps = reps or toy_c1()[1]
    prot = pd.DataFrame(columns=list(PROTEIN_HITS_COLUMNS))
    fdr = pd.DataFrame(columns=list(MAPPING_FDR_COLUMNS))
    write_mapping_outputs(out, ann, prot, [], fdr, reps, SHA)
    return out


def gff_rows(path: Path) -> list[list[str]]:
    return [ln.split("\t") for ln in path.read_text().splitlines() if ln and ln[0] != "#"]


def attrs(col9: str) -> list[tuple[str, str]]:
    return [tuple(kv.split("=", 1)) for kv in col9.split(";")]  # type: ignore[misc]


def feature_of(out: Path, peptide: str) -> list[str]:
    rows = [
        r
        for r in gff_rows(out / "mappings.gff3")
        if r[2] == "nucleotide_to_protein_match" and dict(attrs(r[8]))["Name"] == peptide
    ]
    assert len(rows) == 1, rows
    return rows[0]


def tsv(out: Path) -> pd.DataFrame:
    return pd.read_csv(out / "peptide_hits.tsv", sep="\t", dtype=str, keep_default_na=False)


ROWS = [hit(EX1, "c1", 33, 78), hit("MSDKPLTEAG", "c1", 0, 30, route="exact")]


# ---------------------------------------------------------------- the flag


@pytest.mark.parametrize(
    ("relation", "want"), [("same", "readthrough"), ("related", "readthrough;strain_difference")]
)
def test_example_11_flag_and_alphabetical_order(relation: str, want: str) -> None:
    ann = annotate(ROWS, [(G1_SITE,), ()], relation)
    by = ann.set_index("peptide")
    assert by.loc[EX1, "flags"] == want
    assert "readthrough" not in str(by.loc["MSDKPLTEAG", "flags"] or "")
    assert by.loc[EX1, "locus_class"] == "cterm_extension"


def test_example_12_window_ending_on_the_stop_is_flagged_and_cds_inframe() -> None:
    ann = annotate([hit("GRWQLNEFKW", "c1", 27, 57)], [(G1_SITE,)])
    (r,) = ann.to_dict("records")
    assert r["flags"] == "readthrough" and r["locus_class"] == "cds_inframe"


def test_no_consumed_argument_flags_nothing() -> None:
    ann = annotate(ROWS, None)
    assert not any("readthrough" in str(f or "") for f in ann["flags"])


def test_consumed_aligns_with_input_rows_not_the_sorted_output() -> None:
    """``annotate_peptide_hits`` sorts by peptide (``MSDKPLTEAG`` before ``WQLN...``);
    ``consumed`` follows the *input* order (``WQLN...`` first)."""
    ann = annotate(ROWS, [(G1_SITE,), ()])
    assert list(ann["peptide"]) == ["MSDKPLTEAG", EX1]
    assert [("readthrough" in str(f or "")) for f in ann["flags"]] == [False, True]


def test_recoded_sites_never_give_readthrough() -> None:
    """Only ``terminal`` sites give ``readthrough`` (rule 16; 16c maps ``recoded`` to
    ``sec_site``/``pyl_site`` instead, rule 19)."""
    rec = StopSite("c1", "+", 54, ("G1",), "TGA", "recoded")
    ann = annotate([ROWS[0]], [(rec,)])
    assert "readthrough" not in str(ann["flags"].iloc[0] or "")


def test_input_flags_are_merged_with_readthrough() -> None:
    ann = annotate([hit(EX1, "c1", 33, 78, flags="low_complexity")], [(G1_SITE,)])
    assert ann["flags"].iloc[0] == "low_complexity;readthrough"


def test_work_column_set_exactly_on_flagged_rows() -> None:
    ann = annotate(ROWS, [(G1_SITE,), ()])
    assert "readthrough_codon" in ann.columns
    by = ann.set_index("peptide")["readthrough_codon"]
    assert not _missing(by.loc[EX1])
    assert _missing(by.loc["MSDKPLTEAG"])


def _missing(v: object) -> bool:
    if v is None or v is pd.NA:
        return True
    if isinstance(v, float) and v != v:
        return True
    return isinstance(v, str | tuple | list) and len(v) == 0


# ---------------------------------------------------------------- the writer


def test_example_11_gff3_attribute_after_flags(tmp_path: Path) -> None:
    out = write(tmp_path / "m", annotate(ROWS, [(G1_SITE,), ()]))
    kv = attrs(feature_of(out, EX1)[8])
    keys = [k for k, _ in kv]
    assert dict(kv)["flags"] == "readthrough"
    assert dict(kv)["readthrough_codon"] == "55"
    assert keys.index("readthrough_codon") == keys.index("flags") + 1
    plain = dict(attrs(feature_of(out, "MSDKPLTEAG")[8]))
    assert "readthrough_codon" not in plain


def test_existing_attribute_order_is_unchanged(tmp_path: Path) -> None:
    rows = [hit(EX1, "c1", 33, 78, flags="low_complexity"), ROWS[1]]
    a = write(tmp_path / "a", annotate(rows, None))
    b = write(tmp_path / "b", annotate(rows, [(G1_SITE,), ()]))
    ka = [k for k, _ in attrs(feature_of(a, EX1)[8])]
    kb = [k for k, _ in attrs(feature_of(b, EX1)[8])]
    assert kb == [*ka, "readthrough_codon"]


def test_tsv_columns_unchanged_and_work_column_never_written(tmp_path: Path) -> None:
    out = write(tmp_path / "m", annotate(ROWS, [(G1_SITE,), ()]))
    header = (out / "peptide_hits.tsv").read_text().splitlines()[0].split("\t")
    assert tuple(header) == PEPTIDE_HITS_COLUMNS
    assert "readthrough_codon" not in PEPTIDE_HITS_COLUMNS
    assert "readthrough_codon" not in (out / "peptide_hits.tsv").read_text()
    t = tsv(out).set_index("peptide")
    assert t.loc[EX1, "flags"] == "readthrough"
    assert (t.loc[EX1, "start"], t.loc[EX1, "end"]) == ("34", "78")


def test_bed_probed_and_unmapped_unchanged(tmp_path: Path) -> None:
    a = write(tmp_path / "a", annotate(ROWS, None))
    b = write(tmp_path / "b", annotate(ROWS, [(G1_SITE,), ()]))
    for name in (
        "mappings.bed",
        "mappings.probed",
        "unmapped.tsv",
        "protein_hits.tsv",
        "mapping_fdr.tsv",
    ):
        assert (a / name).read_bytes() == (b / name).read_bytes(), name
    ta, tb = tsv(a), tsv(b)
    assert list(ta.columns) == list(tb.columns)
    others = [c for c in ta.columns if c != "flags"]
    assert ta[others].equals(tb[others])


def test_example_12_attribute_on_a_window_ending_on_the_stop(tmp_path: Path) -> None:
    out = write(tmp_path / "m", annotate([hit("GRWQLNEFKW", "c1", 27, 57)], [(G1_SITE,)]))
    f = feature_of(out, "GRWQLNEFKW")
    assert (f[3], f[4]) == ("28", "57")
    assert dict(attrs(f[8]))["readthrough_codon"] == "55"


def _circular_world() -> Any:
    """Circular contig ``o`` (300 nt) with a ``+`` gene ``O1`` whose TAA stop is ``[0, 3)``
    (the gene runs ``[240, 303)`` unwrapped)."""
    reps = {"o": Replicon("o", 300, "circular")}
    seq = ("GCT" * 100)[:300]
    seq = "TAA" + seq[3:]
    return {"o": seq}, reps, model(reps, [gene("O1", "o", 240, 303, "+")])


def test_example_15_origin_copy_attribute_and_origin_spanning(tmp_path: Path) -> None:
    """Unwrapped ``[282, 312) +`` over the site at ``codon_start`` 0: ``readthrough_codon=1``
    (wrapped ``codon_start + 1``), plus ``origin_spanning`` as today."""
    world = _circular_world()
    site = StopSite("o", "+", 0, ("O1",), "TAA", "terminal")
    ann = annotate([hit("AAAAAAWAAA", "o", 282, 312)], [(site,)], world=world)
    assert {"origin_spanning", "readthrough"} <= set(ann["flags"].iloc[0].split(";"))
    out = write(tmp_path / "m", ann, world[1])
    kv = dict(attrs(feature_of(out, "AAAAAAWAAA")[8]))
    assert kv["readthrough_codon"] == "1"
    assert {"origin_spanning", "readthrough"} <= set(kv["flags"].split(","))


def test_two_cells_one_row_list_both_positions_ascending(tmp_path: Path) -> None:
    """Edge case "Two cells in one row": one flag, both 1-based positions, ``,``-separated."""
    a = StopSite("c1", "+", 15, ("A",), "TAA", "terminal")
    b = StopSite("c1", "+", 30, ("B",), "TGA", "terminal")
    ann = annotate(
        [hit("AAAAAWAAAAWAAAA", "c1", 0, 45, diff_types="substitution;substitution")], [(a, b)]
    )
    assert ann["flags"].iloc[0] == "readthrough"
    out = write(tmp_path / "m", ann)
    assert dict(attrs(feature_of(out, "AAAAAWAAAAWAAAA")[8]))["readthrough_codon"] == "16,31"


def test_multi_locus_peptide_flag_is_per_row(tmp_path: Path) -> None:
    """Edge case "Multi-locus peptide": only the row that consumes a cell is flagged."""
    rows = [hit(EX1, "c1", 33, 78), hit(EX1, "c1", 0, 45, diff_types="swap")]
    ann = annotate(rows, [(G1_SITE,), ()])
    flagged = ann[ann["flags"].fillna("").str.contains("readthrough")]
    assert list(zip(flagged["start"], flagged["end"], strict=True)) == [(33, 78)]
    out = write(tmp_path / "m", ann)
    feats = [r for r in gff_rows(out / "mappings.gff3") if r[2] == "nucleotide_to_protein_match"]
    with_attr = [(r[3], r[4]) for r in feats if "readthrough_codon=" in r[8]]
    assert with_attr == [("34", "78")]


def test_frameshift_rows_and_match_parts_never_carry_it(tmp_path: Path) -> None:
    fs = hit("MSDKPLTEAGRWQLNE", "c1", 0, 49, diff_types="frameshift", blocks="0-24;25-49")
    ann = annotate([*ROWS, fs], [(G1_SITE,), (), ()])
    out = write(tmp_path / "m", ann)
    for r in gff_rows(out / "mappings.gff3"):
        if r[2] == "match_part" or dict(attrs(r[8])).get("Name") != EX1:
            assert "readthrough_codon=" not in r[8], r
