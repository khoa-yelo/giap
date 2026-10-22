"""M7-05 the wired ``variants`` step fills ``flags`` (docs/specs/M7-05-variant-flags.md rule 7,
acceptance criteria 3-4; GIAP_SPEC.md §7.5.3, §7.5.5, §9.4, B.7).

On the M7-04 giapsim chain (sample-strain peptides mapped to the reference genome, default
``genome_relation: related``):

* every ``variants.tsv`` row has a non-empty ``flags`` cell, ``;``-joined, each flag once, in the
  §9.4 order (``mass_equivalent, snv, multi_nt, strain_difference, database_derived,
  start_codon, stop_codon``, then M7-03's ``vcf_unrepresentable``);
* each cell equals an independent oracle of rules 1-6 built from the row itself (``ref_codon``,
  ``ref_aa``, ``alt_aa``, ``alt_type``) and the gene model. The fixture's peptides are all
  ``database`` with ``search_db: same`` (giapsim), so ``database_derived`` is absent; if a
  peptide with another ``search_db`` ever appears, rows with database support skip the check;
* the VCF ``FLAGS`` INFO field is the ``|``-joined union of the flags of the record's alleles
  (order not pinned here; M7-03 sorts it);
* flags change nothing else: the step's rows equal the M7-04 chain without flags
  (``test_variants_step.oracle``) in every column but ``flags``.
"""

from __future__ import annotations

from collections.abc import Mapping
from itertools import product
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from test_variants_step import (
    Run,
    _cfg,
    _chain,
    _gene_model,
    _replicons,
    _vcf_lines,
    oracle,
    run_variants_step,
    step_inputs,
)

ORDER = (
    "mass_equivalent", "snv", "multi_nt", "strain_difference", "database_derived",
    "start_codon", "stop_codon", "vcf_unrepresentable",
)  # fmt: skip

B5: dict[str, float] = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919, "L": 113.08406, "I": 113.08406, "N": 114.04293,
    "D": 115.02694, "Q": 128.05858, "K": 128.09496, "E": 129.04259, "M": 131.04049,
    "H": 137.05891, "F": 147.06841, "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}  # fmt: skip
MODS = (
    15.994915, 31.989829, 0.984016, 14.015650, 28.031300, 42.046950, 42.010565, 57.021464,
    43.005814, 27.994915, 79.966331, -0.984016, -18.010565, -17.026549, 21.981943,
)  # fmt: skip
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {
    a + b + c: _AAS[16 * i + 4 * j + k]
    for (i, a), (j, b), (k, c) in product(enumerate("TCAG"), repeat=3)
}


def _translate(nt: str) -> str:
    return "".join(CODE.get(nt[i : i + 3], "X") for i in range(0, len(nt) - 2, 3))


def _mass_eq(ref: str, alt: str, tol: float = 0.05) -> bool:
    if any(c not in B5 and c != "J" for c in ref + alt):
        return False
    mr = sum(B5[c] for c in ref.replace("J", "L"))
    ma = sum(B5[c] for c in alt.replace("J", "L"))
    return min(abs(ma - mr - d) for d in (0.0, *MODS)) <= tol + 1e-9


def _snv(codon: str, aa: str) -> bool:
    want = "IL" if aa == "J" else aa
    return any(sum(x != y for x, y in zip(codon, c, strict=True)) <= 1
               for c, a in CODE.items() if a in want)  # fmt: skip


def _starts(gm: Any) -> set[tuple[str, str, int]]:
    out = set()
    for g in gm.genes:
        if g.feature_type not in ("cds", "sorf"):
            continue
        iv, rep = g.interval, gm.replicons[g.interval.contig]
        first = iv.start if iv.strand == "+" else iv.end - 3
        out.add((iv.contig, iv.strand, first % rep.length))
    return out


def expected_flags(
    r: Mapping[str, Any], starts: set[tuple[str, str, int]], db_other: bool
) -> list[str] | None:
    """Rules 1-6 for one ``variants.tsv`` row (1-based ``codon_start``) under ``related``;
    ``None`` when ``database_derived`` cannot be decided from the row (database peptides while
    some peptide of the input is not ``search_db: same``)."""
    kind, alt, ref_nt = r["alt_type"], r["alt_aa"], r["ref_codon"]
    want: set[str] = {"strain_difference"}
    if kind == "substitution":
        if _mass_eq(r["ref_aa"], alt):
            want.add("mass_equivalent")
        want.add("snv" if _snv(ref_nt, alt) else "multi_nt")
        if alt == "M" and (r["contig"], r["strand"], int(r["codon_start"]) - 1) in starts:
            want.add("start_codon")
        if r["ref_aa"] == "*":
            want.add("stop_codon")
    elif kind == "block" and _mass_eq(_translate(ref_nt), alt):
        want.add("mass_equivalent")
    if "vcf_unrepresentable" in r["flags"].split(";"):
        want.add("vcf_unrepresentable")
    if db_other and "database" in r["source"].split(";"):
        return None
    return [f for f in ORDER if f in want]


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    return _chain(Path(tmp_path_factory.mktemp("variants_flags_chain")), strain=True)


@pytest.fixture(scope="module")
def full(chain: dict[str, Path], tmp_path_factory: pytest.TempPathFactory) -> Run:
    out = Path(tmp_path_factory.mktemp("variants_flags"))
    return run_variants_step(step_inputs(chain), out, _cfg())


@pytest.fixture(scope="module")
def table(full: Run) -> pd.DataFrame:
    df = full.frame("variants")
    assert len(df) >= 10, "the fixture is meant to give variant sites"
    return df


def test_every_row_has_flags(table: pd.DataFrame) -> None:
    empty = table[table["flags"].isin(["", "NA"])]
    assert empty.empty, f"{len(empty)} of {len(table)} rows without flags:\n{empty.head()}"


def test_flags_are_semicolon_joined_in_section_9_4_order(table: pd.DataFrame) -> None:
    for cell in table["flags"]:
        items = cell.split(";")
        assert all(items), cell
        assert set(items) <= set(ORDER), cell
        assert items == sorted(set(items), key=ORDER.index), cell


def test_flags_follow_rules_1_to_6(chain: dict[str, Path], table: pd.DataFrame) -> None:
    gm = _gene_model(chain["genes"], _replicons(chain["replicons"]))
    starts = _starts(gm)
    db_other = bool((pd.read_parquet(chain["peptides"])["search_db"] != "same").any())
    checked = 0
    for r in table.to_dict("records"):
        want = expected_flags(r, starts, db_other)
        if want is None:
            continue
        assert r["flags"] == ";".join(want), r
        checked += 1
    assert checked >= 10, "too few rows to check"


def test_fixture_exercises_snv_and_strain_difference(table: pd.DataFrame) -> None:
    flags = [set(c.split(";")) for c in table["flags"]]
    assert all("strain_difference" in f for f in flags), "related: every ALT is a strain difference"
    sub = table["alt_type"] == "substitution"
    assert any("snv" in f for f, s in zip(flags, sub, strict=True) if s)


def test_snv_matches_nt_changes_on_substitutions(table: pd.DataFrame) -> None:
    """On substitutions ``snv`` is exactly ``nt_changes <= 1`` (M7-03's minimal-change codon)."""
    for r in table[table["alt_type"] == "substitution"].to_dict("records"):
        items = r["flags"].split(";")
        assert ("snv" in items) == (int(r["nt_changes"]) <= 1), r
        assert ("multi_nt" in items) == (int(r["nt_changes"]) > 1), r


def _info(text: str) -> dict[str, str]:
    return dict(x.split("=", 1) for x in text.split(";") if "=" in x)


def test_vcf_flags_are_the_pipe_joined_union_of_the_alleles(full: Run, table: pd.DataFrame) -> None:
    _, body = _vcf_lines(full.paths["variants_vcf"])
    assert body
    with_flags = 0
    for line in body:
        chrom, pos, _id, ref, alts = line.split("\t")[:5]
        info = _info(line.split("\t")[7])
        rows = []
        for r in table.to_dict("records"):
            if r["contig"] != chrom or r["alt_codon_vcf"] not in alts.split(","):
                continue
            width = len(r["ref_codon"])
            low = int(r["codon_start"]) - 1
            if r["strand"] == "-" and width == 6:
                low -= 3
            if low + 1 == int(pos) and width == len(ref):
                rows.append(r)
        assert rows, line
        union = {f for r in rows for f in r["flags"].split(";")}
        assert "FLAGS" in info, line
        got = info["FLAGS"].split("|")
        assert len(got) == len(set(got)), line
        assert set(got) == union, (line, union)
        with_flags += 1
    assert with_flags == len(body)


def test_flags_change_nothing_else(full: Run, tmp_path_factory: pytest.TempPathFactory) -> None:
    """Rule 7: the step equals the M7-04 chain without flags in every other column."""
    base = oracle(full.inputs, _cfg(), Path(tmp_path_factory.mktemp("variants_flags_oracle")))
    old = pd.read_csv(base["paths"]["variants"], sep="\t", dtype=str, keep_default_na=False)
    new = full.frame("variants")
    others = [c for c in new.columns if c != "flags"]
    pd.testing.assert_frame_equal(new[others], old[others])
    for a, b in zip(new["flags"], old["flags"], strict=True):
        assert set(x for x in b.split(";") if x not in ("", "NA")) <= set(a.split(";"))
