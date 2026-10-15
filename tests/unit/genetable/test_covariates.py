"""M4-10 protein covariates: examples 1-10, edge cases, covariate_table (docs/specs/M4-10-...).

Spec: GIAP_SPEC.md §7.3.3 (covariates, pseudogene rule), §7.3.4 (observable peptides), §5.2
(`abundance.enzyme`, `abundance.observable_len`), §9.3 (`genes.tsv` protein columns), §2.4 G6.
Hand-computed values are given next to each example; Biopython 1.88 is the reference for GRAVY,
pI and MW (spec: "equal to Biopython 1.88 to 1e-9").
"""

from __future__ import annotations

import dataclasses
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import covariates_oracle as oracle

COV_COLUMNS = ("length_aa", "n_observable_peptides", "gravy", "pi", "mw_da", "tm_segments")
EX1 = "MAGLLKAAGPEPTIDERPAAAGGSTLKAAAAAAAK"
EX3 = "MAGLLKAAGPEPTIDERPAAAGGSTLK" + "W" * 31 + "KAAAAAAAK"


def _cov():
    from giap.genetable import covariates as cov

    return cov


# ---------------------------------------------------------------- examples


def test_example1_digest_trypsin_p() -> None:
    assert _cov().digest(EX1, "trypsin_p") == ["MAGLLK", "AAGPEPTIDER", "PAAAGGSTLK", "AAAAAAAK"]


def test_example2_digest_trypsin_no_cleavage_before_p() -> None:
    assert _cov().digest(EX1, "trypsin") == ["MAGLLK", "AAGPEPTIDERPAAAGGSTLK", "AAAAAAAK"]


def test_example3_observable_length_window() -> None:
    cov = _cov()
    # trypsin_p: AAGPEPTIDER (11), PAAAGGSTLK (10), AAAAAAAK (8); MAGLLK (6) and W*31+K (32) out.
    assert cov.n_observable_peptides(EX3, "trypsin_p") == 3
    # trypsin: AAGPEPTIDERPAAAGGSTLK (21), AAAAAAAK (8).
    assert cov.n_observable_peptides(EX3, "trypsin") == 2
    assert cov.n_observable_peptides(EX3, "trypsin_p", (7, 30)) == 3


def test_example4_distinct_peptides_counted_once() -> None:
    assert _cov().n_observable_peptides("AAAAAAAKAAAAAAAKGGGGGGGR", "trypsin_p") == 2


def test_example5_physicochemical_matches_biopython() -> None:
    c = _cov().covariates("MKTAYIAKQR", "trypsin_p", (7, 30))
    # GRAVY by hand: (1.9-3.9-0.7+1.8-1.3+4.5+1.8-3.9-3.5-4.5)/10 = -7.8/10.
    assert c.gravy == pytest.approx(-0.78, abs=1e-9)
    assert c.pi == pytest.approx(10.288976860046386, abs=1e-9)
    assert c.mw_da == pytest.approx(1209.4613, abs=1e-9)
    assert round(c.pi, 3) == 10.289 and round(c.mw_da, 2) == 1209.46
    assert c.length_aa == 10
    # MK | TAYIAK | QR: no fragment of 7-30 aa.
    assert c.n_observable_peptides == 0
    assert c.tm_segments == 0


def test_example6_single_tm_segment() -> None:
    assert _cov().tm_segments("L" * 25) == 1


def test_example7_two_separated_segments() -> None:
    assert _cov().tm_segments("L" * 20 + "K" * 20 + "L" * 20) == 2


def test_example8_short_gap_merges_runs() -> None:
    # A window with 3 K and 16 L: (16*3.8 - 3*3.9)/19 = 2.584 >= 1.6, so every window passes.
    assert _cov().tm_segments("L" * 20 + "K" * 3 + "L" * 20) == 1


def test_example9_short_and_capped() -> None:
    cov = _cov()
    assert cov.tm_segments("L" * 18) == 0
    assert cov.tm_segments(("L" * 20 + "K" * 20) * 12) == 10
    assert cov.tm_segments(("L" * 20 + "K" * 20) * 12, cap=1000) == 12


@pytest.mark.parametrize("missing", [None, pd.NA, np.nan])
def test_example10_pseudogene_translation(missing: object) -> None:
    nt = "GTG" + "GCT" * 5 + "TAA" + "GCT" * 3 + "TAA"
    assert _cov().protein_sequence(missing, nt, "pseudogene") == "MAAAAA"  # type: ignore[arg-type]


# ---------------------------------------------------------------- digest details


def test_digest_basic_boundaries() -> None:
    cov = _cov()
    assert cov.digest("AAK", "trypsin_p") == ["AAK"]  # no empty trailing fragment
    assert cov.digest("KR", "trypsin_p") == ["K", "R"]
    assert cov.digest("KP", "trypsin_p") == ["K", "P"]
    assert cov.digest("KP", "trypsin") == ["KP"]
    assert cov.digest("AKPARPGK", "trypsin") == ["AKPARPGK"]
    assert cov.digest("AKPARPGK", "trypsin_p") == ["AK", "PAR", "PGK"]
    assert cov.digest("MSTV", "trypsin") == ["MSTV"]


def test_digest_missed_cleavages_adds_concatenations() -> None:
    cov = _cov()
    one = cov.digest(EX1, "trypsin_p", missed=1)
    assert sorted(one) == sorted(
        [
            *["MAGLLK", "AAGPEPTIDER", "PAAAGGSTLK", "AAAAAAAK"],
            *["MAGLLKAAGPEPTIDER", "AAGPEPTIDERPAAAGGSTLK", "PAAAGGSTLKAAAAAAAK"],
        ]
    )
    two = cov.digest(EX1, "trypsin_p", missed=2)
    assert set(two) == set(one) | {
        "MAGLLKAAGPEPTIDERPAAAGGSTLK",
        "AAGPEPTIDERPAAAGGSTLKAAAAAAAK",
    }
    assert len(two) == 4 + 3 + 2
    assert cov.digest(EX1, "trypsin_p", missed=0) == cov.digest(EX1, "trypsin_p")


def test_observable_length_bounds_inclusive() -> None:
    cov = _cov()
    seq = "A" * 6 + "K" + "C" * 29 + "K" + "D" * 30 + "K" + "E" * 5 + "K"  # 7, 30, 31, 6
    assert cov.n_observable_peptides(seq, "trypsin_p") == 2
    assert cov.n_observable_peptides(seq, "trypsin_p", (8, 31)) == 2  # 30 and 31
    assert cov.n_observable_peptides(seq, "trypsin_p", (6, 31)) == 4


def test_observable_i_and_l_kept_distinct() -> None:
    assert _cov().n_observable_peptides("AAAAAAIK" + "AAAAAALK", "trypsin_p") == 2


@pytest.mark.parametrize("bad", ["chymotrypsin", "Trypsin", "lys-c", ""])
def test_unknown_enzyme_rejected(bad: str) -> None:
    cov = _cov()
    with pytest.raises(ValueError):
        cov.digest(EX1, bad)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        cov.n_observable_peptides(EX1, bad)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        cov.covariates(EX1, bad, (7, 30))  # type: ignore[arg-type]


# ---------------------------------------------------------------- tm_segments details


def test_tm_exact_window_length() -> None:
    cov = _cov()
    assert cov.tm_segments("L" * 19) == 1
    assert cov.tm_segments("") == 0


def test_tm_unknown_letters_score_zero() -> None:
    cov = _cov()
    assert cov.tm_segments("X" * 19) == 0
    assert cov.tm_segments("A" * 19) == 1  # 1.8 >= 1.6
    # 17 A + 2 X: 30.6/19 = 1.611 >= 1.6; 16 A + 3 X: 28.8/19 = 1.516 < 1.6 (unknowns count as 0,
    # they are not skipped: skipping would give 1.8 for both).
    assert cov.tm_segments("A" * 17 + "XX") == 1
    assert cov.tm_segments("A" * 16 + "XXX") == 0
    assert cov.tm_segments("A" * 16 + "UBZ") == 0


def test_tm_parameters_respected() -> None:
    cov = _cov()
    assert cov.tm_segments("L" * 10, window=10) == 1
    assert cov.tm_segments("L" * 10, window=11) == 0
    assert cov.tm_segments("A" * 19, threshold=2.0) == 0
    assert cov.tm_segments(("L" * 20 + "K" * 20) * 4, cap=3) == 3


def test_tm_merging_boundary_gap() -> None:
    cov = _cov()
    # A window with k K and 19-k L passes iff 72.2 - 7.7 k >= 30.4, i.e. k <= 5. With a 6-K gap the
    # windows holding all 6 K have mean (49.4 - 23.4)/19 = 1.368 < 1.6, so the runs split; with a
    # 5-K gap every window passes ((53.2 - 19.5)/19 = 1.774).
    assert cov.tm_segments("L" * 20 + "K" * 6 + "L" * 20) == 2
    assert cov.tm_segments("L" * 20 + "K" * 5 + "L" * 20) == 1
    for gap in range(0, 25):
        s = "L" * 20 + "K" * gap + "L" * 20
        assert cov.tm_segments(s) == oracle.tm_segments(s), gap


# ---------------------------------------------------------------- protein_sequence


def test_protein_sequence_prefers_aa() -> None:
    assert _cov().protein_sequence("MKT", "ATGAAAACC", "cds") == "MKT"
    assert _cov().protein_sequence("MKTW", "ATGAAAACCTGG", "pseudogene") == "MKTW"


def test_protein_sequence_frameshift_whole_codons_only() -> None:
    cov = _cov()
    assert cov.protein_sequence(None, "ATG" + "GCT" * 2 + "GC", "pseudogene") == "MAA"
    assert cov.protein_sequence(None, "TTG" + "AAA" + "T", "pseudogene") == "MK"


def test_protein_sequence_table11_stops_and_start() -> None:
    cov = _cov()
    assert cov.protein_sequence(None, "ATGTGGTGAGCT", "pseudogene") == "MW"  # TGA stops
    assert cov.protein_sequence(None, "ATGTGGTAGGCT", "pseudogene") == "MW"  # TAG stops
    assert cov.protein_sequence(None, "TTGCTG", "pseudogene") == "ML"  # TTG start read as M
    assert cov.protein_sequence(None, "ATGGCTGCT", "pseudogene") == "MAA"  # no stop: all codons
    assert cov.protein_sequence(None, "ATG", "pseudogene") == "M"


# ---------------------------------------------------------------- covariates edge cases


def test_covariates_fields_and_frozen() -> None:
    cov = _cov()
    c = cov.covariates("MKTAYIAKQR", "trypsin_p", (7, 30))
    assert isinstance(c, cov.Covariates)
    assert tuple(f.name for f in dataclasses.fields(c)) == COV_COLUMNS
    with pytest.raises(dataclasses.FrozenInstanceError):
        c.gravy = 0.0  # type: ignore[misc]


def test_covariates_trailing_and_leading_stop_stripped() -> None:
    cov = _cov()
    ref = cov.covariates("MKTAYIAKQR", "trypsin_p", (7, 30))
    assert cov.covariates("MKTAYIAKQR*", "trypsin_p", (7, 30)) == ref
    assert cov.covariates("*MKTAYIAKQR", "trypsin_p", (7, 30)) == ref


@pytest.mark.parametrize("extra", ["X", "B", "Z", "J", "U", "O"])
def test_covariates_nonstandard_letters_removed_but_counted(extra: str) -> None:
    cov = _cov()
    seq = "MKTAY" + extra + "IAKQR"
    c = cov.covariates(seq, "trypsin_p", (7, 30))
    assert c.length_aa == 11
    assert c.gravy == pytest.approx(-0.78, abs=1e-9)
    assert c.pi == pytest.approx(10.288976860046386, abs=1e-9)
    assert c.mw_da == pytest.approx(1209.4613, abs=1e-9)


def test_covariates_selenoprotein() -> None:
    cov = _cov()
    seq = "MSTVLAKGUCDEFGHIKR*"
    c = cov.covariates(seq, "trypsin_p", (7, 30))
    g, p, m = oracle.biopython_props(seq)
    assert c.length_aa == 18
    assert (c.gravy, c.pi, c.mw_da) == pytest.approx((g, p, m), abs=1e-9)


def test_covariates_very_short_protein() -> None:
    cov = _cov()
    c = cov.covariates("M", "trypsin_p", (7, 30))
    _g, p, m = oracle.biopython_props("M")
    assert c.length_aa == 1 and c.n_observable_peptides == 0 and c.tm_segments == 0
    assert c.gravy == pytest.approx(1.9, abs=1e-12)
    assert (c.pi, c.mw_da) == pytest.approx((p, m), abs=1e-9)


def test_covariates_empty_protein_nan() -> None:
    c = _cov().covariates("", "trypsin_p", (7, 30))
    assert c.length_aa == 0
    assert math.isnan(c.gravy) and math.isnan(c.pi) and math.isnan(c.mw_da)
    # "other covariates NaN": counts may be 0 or NaN, never a positive number.
    for v in (c.n_observable_peptides, c.tm_segments):
        assert v == 0 or math.isnan(v)


@pytest.mark.parametrize("seq", ["XXXX", "*", "X*"])
def test_covariates_nothing_left_after_cleaning(seq: str) -> None:
    c = _cov().covariates(seq, "trypsin_p", (7, 30))
    assert math.isnan(c.gravy) and math.isnan(c.pi) and math.isnan(c.mw_da)


@pytest.mark.parametrize("seq", ["mktayiakqr", "MKTAYiAKQR"])
def test_covariates_lower_case_rejected(seq: str) -> None:
    with pytest.raises(ValueError):
        _cov().covariates(seq, "trypsin_p", (7, 30))


def test_covariates_observable_uses_given_range_and_enzyme() -> None:
    cov = _cov()
    assert cov.covariates(EX3, "trypsin_p", (7, 30)).n_observable_peptides == 3
    assert cov.covariates(EX3, "trypsin", (7, 30)).n_observable_peptides == 2
    assert cov.covariates(EX3, "trypsin_p", (9, 30)).n_observable_peptides == 2
    assert cov.covariates(EX3, "trypsin_p", (7, 32)).n_observable_peptides == 4


def test_covariates_deterministic() -> None:
    cov = _cov()
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ"
    assert cov.covariates(seq, "trypsin_p", (7, 30)) == cov.covariates(seq, "trypsin_p", (7, 30))


# ---------------------------------------------------------------- covariate_table

C2_COLUMNS = (
    "locus_tag", "gene", "product", "feature_type", "contig", "start", "end", "strand",
    "partial", "crosses_origin", "hypothetical", "aa", "nt", "start_codon",
)  # fmt: skip


def _back_translate(aa: str) -> str:
    return "".join(oracle.BACK[c] for c in aa) + "TAA"


def build_genes(n: int = 1000, seed: int = 7) -> pd.DataFrame:
    """C2-shaped table on one contig: log-normal lengths (median 280 aa), every 50th gene a
    pseudogene with ``aa = NA`` and an internal stop, every 40th a sORF (< 30 aa)."""
    rng = np.random.default_rng(seed)
    letters = np.array(list(oracle.STANDARD_AA))
    rows = []
    pos = 100
    for i in range(n):
        if i % 40 == 39:
            length = int(rng.integers(10, 30))
            ftype = "sorf"
        else:
            length = max(31, round(280 * math.exp(0.45 * rng.standard_normal())))
            ftype = "cds"
        aa = "M" + "".join(rng.choice(letters, size=length - 1))
        nt = _back_translate(aa)
        if i % 50 == 25:
            ftype = "pseudogene"
            cut = 3 * (length // 2)
            nt = nt[:cut] + "TAA" + nt[cut:]
            aa_val: object = pd.NA
        else:
            aa_val = aa
        start, end = pos, pos + len(nt)
        pos = end + int(rng.integers(50, 300))
        rows.append(
            {
                "locus_tag": f"GSIM_{5 * (i + 1):05d}",
                "gene": pd.NA,
                "product": "hypothetical protein",
                "feature_type": ftype,
                "contig": "chr1",
                "start": start,
                "end": end,
                "strand": "+" if (i // 3) % 2 == 0 else "-",
                "partial": False,
                "crosses_origin": False,
                "hypothetical": True,
                "aa": aa_val,
                "nt": nt,
                "start_codon": nt[:3],
            }
        )
    return pd.DataFrame(rows, columns=list(C2_COLUMNS))


def run_table(genes: pd.DataFrame, cfg: object | None = None) -> pd.DataFrame:
    from giap.config import GiapConfig
    from giap.genetable.covariates import covariate_table

    return covariate_table(genes, cfg if cfg is not None else GiapConfig())


def _expected_row(aa: object, nt: str, enzyme: str, length: tuple[int, int]) -> dict[str, float]:
    seq = oracle.translate_pseudogene(nt) if aa is None or pd.isna(aa) else str(aa)
    g, p, m = oracle.biopython_props(seq)
    return {
        "length_aa": len(seq.strip("*")),
        "n_observable_peptides": oracle.n_observable(seq, enzyme, length),
        "gravy": g,
        "pi": p,
        "mw_da": m,
        "tm_segments": oracle.tm_segments(seq),
    }


def test_covariate_table_columns_rows_order_and_values() -> None:
    genes = build_genes(120)
    out = run_table(genes)
    assert list(out.columns) == ["locus_tag", *COV_COLUMNS]
    assert list(out["locus_tag"]) == list(genes["locus_tag"])
    for (_, g), (_, r) in zip(genes.iterrows(), out.iterrows(), strict=True):
        exp = _expected_row(g["aa"], g["nt"], "trypsin_p", (7, 30))
        for k in ("length_aa", "n_observable_peptides", "tm_segments"):
            assert int(r[k]) == exp[k], (g["locus_tag"], k)
        for k in ("gravy", "pi", "mw_da"):
            assert float(r[k]) == pytest.approx(exp[k], abs=1e-9), (g["locus_tag"], k)


def test_covariate_table_pseudogene_rule() -> None:
    genes = build_genes(60)
    pseudo = genes[genes["aa"].isna()]
    assert len(pseudo) >= 1
    out = run_table(genes).set_index("locus_tag")
    for _, g in pseudo.iterrows():
        seq = oracle.translate_pseudogene(g["nt"])
        assert int(out.loc[g["locus_tag"], "length_aa"]) == len(seq)
        assert seq.startswith("M") and "*" not in seq


def test_covariate_table_example10_row() -> None:
    genes = build_genes(3)
    nt = "GTG" + "GCT" * 5 + "TAA" + "GCT" * 3 + "TAA"
    genes.loc[1, ["feature_type", "aa", "nt", "start_codon"]] = ["pseudogene", pd.NA, nt, "GTG"]
    out = run_table(genes)
    row = out.iloc[1]
    g, p, m = oracle.biopython_props("MAAAAA")
    assert int(row["length_aa"]) == 6
    assert (float(row["gravy"]), float(row["pi"]), float(row["mw_da"])) == pytest.approx(
        (g, p, m), abs=1e-9
    )


def test_covariate_table_uses_config_enzyme_and_length() -> None:
    from giap.config import GiapConfig

    genes = build_genes(3)
    genes.loc[0, ["aa", "nt"]] = [EX3, _back_translate(EX3)]
    default = run_table(genes)
    assert int(default.loc[0, "n_observable_peptides"]) == 3
    trypsin = run_table(genes, GiapConfig.model_validate({"abundance": {"enzyme": "trypsin"}}))
    assert int(trypsin.loc[0, "n_observable_peptides"]) == 2
    narrow = run_table(genes, GiapConfig.model_validate({"abundance": {"observable_len": [9, 30]}}))
    assert int(narrow.loc[0, "n_observable_peptides"]) == 2


def test_covariate_table_empty_protein_kept() -> None:
    genes = build_genes(3)
    genes.loc[2, ["feature_type", "aa"]] = ["pseudogene", ""]
    out = run_table(genes)
    assert len(out) == 3
    assert int(out.loc[2, "length_aa"]) == 0
    assert out.loc[2, ["gravy", "pi", "mw_da"]].isna().all()


def test_covariate_table_does_not_mutate_input_and_is_deterministic() -> None:
    genes = build_genes(80)
    before = genes.copy()
    a = run_table(genes)
    b = run_table(genes)
    pd.testing.assert_frame_equal(genes, before)
    pd.testing.assert_frame_equal(a, b)


def time_thousand(genes: pd.DataFrame) -> float:
    t0 = time.perf_counter()
    out = run_table(genes)
    secs = time.perf_counter() - t0
    assert len(out) == len(genes)
    return secs


def test_covariate_table_thousand_genes_loose_guard() -> None:
    # Budget is < 2 s (benchmarks/test_covariates.py); this per-commit guard is 10x looser.
    genes = build_genes(1000)
    out = run_table(genes)
    assert len(out) == 1000 and out["locus_tag"].is_unique
    assert time_thousand(genes) < 20.0
