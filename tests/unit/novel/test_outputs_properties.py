"""Properties of the M7-07 novel-ORF outputs (docs/specs/M7-07-novel-orf-outputs.md rules 2, 4-7;
lead decision 1): example 1 on either strand, rotated to any origin position (so ORFs cross the
origin), with any OSTIR choice, gives transcription-direction TSV coordinates, a GFF3/patch
interval from the chosen start to the stop, and proteins that do not depend on strand or origin.
Row order does not depend on the order of the inputs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from helpers.novel_outputs_fixture import (
    L,
    ex1,
    feature_lines,
    mirror_case,
    orfs_of,
    outputs,
    protein_from,
    proteins,
    rotate_case,
    scored,
    table,
    tsv_rows,
    two_contigs,
    write_all,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

PLUS = {"longest": 42, "gtg": 60, "stop": 120}
MINUS = {"longest": 955, "gtg": 937, "stop": 877}
N_AA = {"longest": 26, "gtg": 20}
SETTINGS = settings(
    max_examples=40, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)


def _expected_proteins() -> dict[str, str]:
    seq = ex1().genome["c1"]
    return {"longest": protein_from(seq, 42, 120, "+"), "gtg": protein_from(seq, 60, 120, "+")}


@SETTINGS
@given(
    minus=st.booleans(),
    shift=st.integers(0, L - 1),
    choice=st.sampled_from([None, "longest", "gtg"]),
)
def test_coordinates_lengths_and_proteins(
    tmp_path_factory: object, minus: bool, shift: int, choice: str | None
) -> None:
    base = mirror_case(ex1()) if minus else ex1()
    case = rotate_case(base, shift)
    pos = {k: (v + shift) % L for k, v in (MINUS if minus else PLUS).items()}
    (orf,) = orfs_of(case)
    assert orf.start_longest == pos["longest"] and orf.stop == pos["stop"]
    best = None if choice is None else pos[choice]
    t = table(case, [orf], scored([orf], best={orf.orf_id: best}, called={orf.orf_id: True}))
    out = tmp_path_factory.mktemp("p")  # type: ignore[attr-defined]
    paths = write_all(Path(out), case, t)
    (row,) = tsv_rows(paths["novel_orfs.tsv"])

    def first_base(x: int) -> int:  # 1-based first base of the codon at x, transcription order
        return (x + 2) % L + 1 if minus else x + 1

    def last_base_stop(x: int) -> int:
        return x + 1 if minus else (x + 2) % L + 1

    assert row["start_longest"] == str(first_base(pos["longest"]))
    assert row["end"] == str(last_base_stop(pos["stop"]))
    assert row["start_best_ostir"] == (None if best is None else str(first_base(best)))
    used = choice or "longest"
    assert (row["length_aa_longest"], row["length_aa_best"]) == ("26", str(N_AA[used]))
    assert row["start_codon_best"] == ("GTG" if used == "gtg" else "ATG")
    assert ("ostir_missing" in (row["flags"] or "").split(";")) == (best is None)

    span = 3 * (N_AA[used] + 1)
    lo = pos["stop"] if minus else pos[used]
    want_cols = [str(lo + 1), str(lo + span), ".", "-" if minus else "+", "0"]
    (g,) = feature_lines(paths["novel_orfs.gff3"])
    (p,) = feature_lines(paths["annotation_patch.gff3"])
    assert g[3:8] == want_cols and p[3:8] == want_cols

    prot = _expected_proteins()
    got = proteins(case, t)
    if choice in (None,):
        assert got == {f"{orf.orf_id}|longest": prot["longest"]}
    elif choice == "longest":
        assert got == {f"{orf.orf_id}|longest,best_ostir": prot["longest"]}
    else:
        assert got == {f"{orf.orf_id}|longest": prot["longest"],
                       f"{orf.orf_id}|best_ostir": prot["gtg"]}  # fmt: skip


@settings(max_examples=30, deadline=None)
@given(perm=st.permutations(range(3)), sperm=st.permutations(range(3)))
def test_table_does_not_depend_on_input_order(perm: list[int], sperm: list[int]) -> None:
    case = two_contigs()
    orfs = orfs_of(case)
    sc = scored(orfs, best={"NOVEL_c1_43_p": 60})
    ref = table(case, orfs, sc)
    got = table(case, [orfs[i] for i in perm], sc.iloc[list(sperm)].reset_index(drop=True))
    pd.testing.assert_frame_equal(got.reset_index(drop=True), ref.reset_index(drop=True))
    assert outputs().gff_features(got) == outputs().gff_features(ref)
