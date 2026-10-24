# ruff: noqa: RUF001, RUF002, RUF003  (± and – are the plots' own characters)
"""K-12 plot survey, batch 1 item 5: numbers in titles and labels.

V02, V03, V10, V13 and V15 printed raw floats ("Moran's I 0.492986 (q = 0.00285429)"); the report
index shows the same titles (it reads each SVG's ``<title>``) and the V3 locus index. V15 printed a
small, non-zero SD as "± 0.00" and a narrow interval as "0.02–0.02".

Interface assumed (test-author choice, reported to the lead): ``giap.viz.base.sig(value,
digits=2) -> str`` is the shared formatter (V07's "q = 0.23" style): ``value`` rounded to
``digits`` significant digits, no exponent for ``1e-3 <= |value| < 1e4``, ``"NA"`` for ``None`` /
NaN (the ``base.title_pct`` NA convention). Titles may use 2 or 3 significant digits, so the plot
tests check that every decimal number has at most 3 significant digits and stays within rounding
of the input.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd
import pytest
from helpers import viz_determinants as vd
from helpers import viz_locus as vl
from helpers.viz_landscape import correlogram_rows, gene, island, landscape_out, plot_spec
from hypothesis import given
from hypothesis import strategies as st
from matplotlib.text import Text

from giap.config import GiapConfig, resolve_config
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.base import Instance, Skip

MINUS = "−"
_DEC = re.compile(r"[-−]?\d*\.\d+(?:e[-+]?\d+)?")


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


def _parse(s: str) -> float:
    return float(s.strip().replace(MINUS, "-").replace(",", ""))


def _sig_digits(token: str) -> int:
    mant = token.replace(MINUS, "").replace("-", "").split("e")[0]
    digits = mant.replace(".", "").lstrip("0")
    if "." not in mant:
        digits = digits.rstrip("0")
    return len(digits)


def _too_precise(text: str) -> list[str]:
    """Decimal numbers in ``text`` with more than 3 significant digits."""
    return [t for t in _DEC.findall(text) if _sig_digits(t) > 3]


def _sig() -> object:
    fn = getattr(base, "sig", None)
    assert callable(fn), "giap.viz.base.sig (the shared number formatter) is missing"
    return fn


# --------------------------------------------------------------------- the shared helper


@pytest.mark.parametrize(
    ("value", "want"),
    [(0.00285429, "0.0029"), (0.492986, "0.49"), (0.23, "0.23"), (0.0412, "0.041"),
     (0.00123, "0.0012"), (0.5, "0.5")],
)  # fmt: skip
def test_sig_two_significant_digits(value: float, want: str) -> None:
    assert _sig()(value) == want  # type: ignore[operator]


def test_sig_three_digits_when_asked() -> None:
    assert _sig()(0.00285429, 3) == "0.00285"  # type: ignore[operator]
    assert _sig()(0.492986, 3) == "0.493"  # type: ignore[operator]


@pytest.mark.parametrize("value", [None, math.nan, pd.NA])
def test_sig_na_is_the_house_na(value: object) -> None:
    assert _sig()(value) == "NA"  # type: ignore[operator]


def test_sig_negative_keeps_its_sign() -> None:
    assert _parse(_sig()(-0.0312)) == pytest.approx(-0.031)  # type: ignore[operator]


@pytest.mark.parametrize("value", [0.001, 0.0123, 1.0, 12.34, 99.7, 999.9, 1234.5, 9876.0])
def test_sig_has_no_exponent_in_the_title_range(value: float) -> None:
    out = _sig()(value)  # type: ignore[operator]
    assert "e" not in out.lower(), out
    assert _parse(out) == pytest.approx(value, rel=0.05)


@given(st.floats(min_value=1e-4, max_value=9.9, allow_nan=False), st.sampled_from([2, 3]))
def test_property_sig_rounds_to_the_digits(value: float, digits: int) -> None:
    out = _sig()(value, digits)  # type: ignore[operator]
    got = _parse(out)
    step = 10 ** (math.floor(math.log10(value)) - digits + 1)
    assert abs(got - value) <= 0.5 * step * (1 + 1e-9) + 1e-15, (value, out)
    assert _sig_digits(out) <= digits, (value, out)


# --------------------------------------------------------------------- V02


def test_v02_title_rounds_morans_i_and_q(tmp_path: Path) -> None:
    cfg1 = resolve_config(overrides=["context.islands.detrend_span_frac=1.0"])
    rows = correlogram_rows(
        overrides={
            ("morans_i", "gene_label", 1): (0.492986, -0.002, -0.05, 0.05, 0.00285429),
            ("morans_i", "tu_block", 1): (0.492986, -0.004, -0.05, 0.12, 0.00285429),
        }
    )
    genes = [gene(f"g{i}", m, y) for i, (m, y) in
             enumerate(zip([10_000, 30_000, 50_000, 70_000, 90_000], [1.0, 2.0, 3.0, 2.0, 1.0],
                           strict=True))]  # fmt: skip
    out = landscape_out(tmp_path, genes=genes,
                        islands=[island("ISL_1", "high", 20_001, 40_000, 0.01)],
                        autocorr=rows)  # fmt: skip
    result = plot_spec("V02").fn(RunView.from_outdir(out), cfg1)
    assert not isinstance(result, Skip)
    title = result[0].title
    assert "Moran's I" in title, title
    assert _too_precise(title) == [], title
    assert ("0.49" in title and "0.0029" in title) or ("0.493" in title and "0.00285" in title)


# --------------------------------------------------------------------- V10


def test_v10_title_rounds_morans_i_and_q(tmp_path: Path, cfg: GiapConfig) -> None:
    rows = correlogram_rows(
        "chrA",
        overrides={
            ("morans_i", "gene_label", 1): (0.492986, -0.002, -0.01, 0.01, 0.001),
            ("morans_i", "tu_block", 1): (0.492986, -0.004, -0.05, 0.12, 0.00285429),
        },
    )
    out = landscape_out(tmp_path, replicons=(("chrA", 100_000, "circular", "chromosome"),),
                        genes=[gene("g1", 10_000, 1.0)], autocorr=rows)  # fmt: skip
    result = plot_spec("V10").fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip)
    title = result[0].title
    assert title.startswith("Lag-1 Moran's I "), title
    assert _too_precise(title) == [], title
    assert "above" in title


# --------------------------------------------------------------------- V03 and its index


def _v03(tmp: Path, cfg: GiapConfig) -> list[Instance]:
    import giap.viz.locus  # noqa: F401  (registers V03)

    out = vl.build_out(
        tmp,
        genes=[vl.gene("g1", "chr1", 10_101, 11_000, log2=20.0),
               vl.gene("g2", "chr1", 40_101, 41_000, log2=21.0)],
        regions_stats=[vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 12_000,
                                      0.00285429)],
        islands=[vl.island("ISL_1", "chr1", 40_001, 42_000, 0.0123456, mean=21.53872)],
    )  # fmt: skip
    result = base.registered()["V03"].fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    return list(result)


def test_v03_locus_titles_round_q_mean_and_p(tmp_path: Path, cfg: GiapConfig) -> None:
    insts = _v03(tmp_path, cfg)
    loci = [i for i in insts if i.suffix is not None]
    assert {i.suffix for i in loci} == {"prophage_1", "ISL_1"}
    for inst in loci:
        assert _too_precise(inst.title) == [], inst.title
        sup = inst.figure.get_suptitle()
        assert _too_precise(sup) == [], sup


def test_v03_index_rounds_the_error_column_text(tmp_path: Path, cfg: GiapConfig) -> None:
    (index,) = [i for i in _v03(tmp_path, cfg) if i.suffix is None]
    texts = [t.get_text() for t in index.figure.findobj(Text) if t.get_text().strip()]
    assert any("prophage_1" in t for t in texts), texts
    bad = [t for t in texts if _too_precise(t)]
    assert bad == [], bad


# --------------------------------------------------------------------- V13 locus label


def test_v13_locus_label_rounds_best_q(tmp_path: Path, cfg: GiapConfig) -> None:
    import giap.viz.novel_orfs  # noqa: F401  (registers V13)

    row = {
        "orf_id": "NOVEL_chr1_5001_p", "class": "intergenic", "contig": "chr1", "strand": "+",
        "frame": 1, "start_longest": 5_001, "start_best_ostir": 5_001, "end": 5_300,
        "length_aa_longest": 99, "length_aa_best": 99, "start_codon_best": "ATG",
        "n_peptides": 2, "peptides": "PEPTWDEK", "best_peptide_q": 0.000123456,
        "class_fdr": "0.00412345", "called": True, "rank": 1,
    }  # fmt: skip
    out = vl.build_out(tmp_path, novel_orfs=[row])
    result = base.registered()["V13"].fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip)
    views = [i for i in result if i.suffix and i.suffix.startswith("locus_")]
    assert len(views) == 1
    assert _too_precise(views[0].title) == [], views[0].title


# --------------------------------------------------------------------- V15


def _v15_title(tmp: Path, cfg: GiapConfig, **kw: object) -> str:
    out = vd.build_out(tmp, vd.model_table(**kw))  # type: ignore[arg-type]
    result = plot_spec("V15").fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    return result[0].title


def test_v15_small_nonzero_sd_is_not_printed_as_zero(tmp_path: Path, cfg: GiapConfig) -> None:
    """K-12: "Elastic net CV R² 0.04 ± 0.00"."""
    t = _v15_title(tmp_path, cfg, en_cv=(0.0412, 0.00123, 0.02))
    m = re.search(r"CV R² (\S+) ± (\S+?)[;,( ]", t)
    assert m, t
    assert _parse(m.group(2)) == pytest.approx(0.00123, rel=0.05), t
    assert _parse(m.group(1)) == pytest.approx(0.0412, rel=0.05), t


def test_v15_narrow_interval_ends_stay_distinct(tmp_path: Path, cfg: GiapConfig) -> None:
    """K-12: "protein 0.02 (repeat t-interval 0.02–0.02, uncalibrated)"."""
    en = {"transcription": (0.005, -0.01, 0.02), "protein": (0.0203, 0.0151, 0.0249)}
    t = _v15_title(tmp_path, cfg, en=en)
    m = re.search(r"protein (\S+) \(repeat t-interval (\S+)–(\S+),", t)
    assert m, t
    v, lo, hi = (_parse(g) for g in m.groups())
    assert lo < hi, t
    assert (v, lo, hi) == pytest.approx((0.0203, 0.0151, 0.0249), rel=0.05), t
    assert _too_precise(t) == [], t
