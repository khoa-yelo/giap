"""PR #264 review (comment 6369206328): the ADR-0022 outputs, the stream order after stream 8, the
``terms.tsv.gz`` pin and the review's cheap gaps (``validation/calibrate_cutoffs.py``, M6-12 B).

Lead decisions pinned here, on top of ``test_m6_12_calibrate_cutoffs.py``:

- **Stream order after stream 8.** The same ``rng = for_step(STEP, seed)`` continues with
  (9) ``threshold_boot_q(p, y, directon, rng=rng, B=B)``; (10) the coincidence null, **one** call
  ``coincidence_rate(eligible calls, confidence, windows, truth, replicon, t*, rng=rng, R=R,
  return_rates=True)`` (eligible calls sorted by (start, end, strand), ``truth`` all on-contig
  truth intervals); (11) ``precision_cv(p, y, directon, rng=rng)`` (only with >= 10 directons);
  (12) the window-cluster CI ``fixed_cutoff_ci(confidence, matched, gene, t*, rng=rng, B=B)``,
  clustered by the call's **gene** (window), not by row. Draws 10 and 12 happen only when the
  TerminatorNet half is calibrated.
- ``coincidence_rate`` = mean(rates); ``coincidence_rate_mc_se`` = sd(rates, ddof=1) / sqrt(R).
  This replaces the 20-call batch means, so any R >= 2 gives a finite SE.
- ``precision_lower_bound`` = lb(precision, c); ``precision_lower_bound_wilson`` =
  lb(Wilson low of tp / (tp + fp), c); ``precision_lower_bound_window`` = lb(window CI low, c);
  lb(lo, c) = max(0, (lo - c) / (1 - c)) (``giap.stats.cutoffs.precision_lower_bound``).
- ``precision_ci_window`` is ``[null, null]`` when uncalibrated, like every other CI key.
- **B2:** ``Pins`` gains ``terms_sha256``, the SHA-256 of the ``terms.tsv.gz`` bytes. ``K12_PINS``
  holds the committed recording's value, which the committed ``PROVENANCE.md`` also lists. A
  different terms file next to a correct ``PROVENANCE.md`` is a ``ValueError``, and nothing is
  written. ``tests/helpers/k12_cutoffs_world.py`` puts ``terms_sha256`` in ``World.pins``.
- ``main`` returns 1 with a one-line ``calibrate_cutoffs:`` message on stderr, without raising,
  when an input file is missing (``FileNotFoundError``).
"""

from __future__ import annotations

import gzip
import importlib
import logging
import math
import random
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.k12_cutoffs_world import (
    CONTIG,
    EX6_CALLS,
    EX6_GENES,
    EX6_LABELS,
    EX6_PROBS,
    EX6_TRUTH,
    LIN10K,
    BigSpec,
    Call,
    G,
    Term,
    World,
    big_spec,
    big_world,
    build_world,
    example6,
    gene_obj,
    refresh_pins,
    sha256,
    write_uniop,
)

from giap.annotate.operons import intergenic_distance
from giap.model.coords import Interval, Replicon
from giap.rng import for_step
from giap.stats.cutoffs import (
    auroc_ci,
    auroc_diff_ci,
    choose_cutoff,
    coincidence_rate,
    distance_fit_ci,
    downstream_windows,
    fit_distance_model,
    fixed_cutoff_ci,
    precision_cv,
    precision_lower_bound,
    terminator_match,
    threshold_boot_q,
    threshold_curve,
    wilson_ci,
    window_of,
)

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
cc: Any = importlib.import_module("validation.calibrate_cutoffs")

SEED = 20261014
B = 30
R = 40  # a multiple of 20 (the reviewer's ask), and not 20 itself
K12_TERMS = REPO / "tests/data/recorded/terminatornet/721803095459/k12/terms.tsv.gz"
K12_TERMS_SHA = "9be61c75c0d7835f7d8c8f3c3948b95adfa34445003ec570ff9bbb2ca1e3a2ff"


@pytest.fixture(autouse=True)
def _no_real_data(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    missing = tmp_path / "no-such-root"
    for var in ("GIAP_ROOT", "GIAP_DB_DIR", "GIAP_DATA_DIR"):
        monkeypatch.setenv(var, str(missing))


def run(w: World, out: Path, **kw: Any) -> dict[str, Any]:
    kw.setdefault("pins", cc.Pins(**w.pins))
    kw.setdefault("B", B)
    kw.setdefault("R", R)
    result: dict[str, Any] = cc.run(*w.run_args(out), **kw)
    return result


def ci(b: Any) -> list[float]:
    return [b.ci_low, b.ci_high]


def approx_ci(got: Any, want: list[float]) -> None:
    assert isinstance(got, list) and len(got) == 2
    for g, w in zip(got, want, strict=True):
        if math.isnan(w):
            assert g is None
        else:
            assert g == pytest.approx(w, rel=1e-9, abs=1e-12)


def approx(got: Any, want: float) -> None:
    if math.isnan(want):
        assert got is None
    else:
        assert got == pytest.approx(want, rel=1e-9, abs=1e-12)


def iv(start: int, end: int, strand: str, contig: str = CONTIG) -> Interval:
    return Interval(contig, start, end, strand)  # type: ignore[arg-type]


def no_outputs(out: Path) -> bool:
    return not (out / "cutoffs.json").exists()


# ============================================ a big world where genes hold several eligible calls


def two_call_spec(seed: int = 7) -> BigSpec:
    """``big_spec`` with a second call in every other window that holds one, and confidences
    set so that the TerminatorNet t* keeps some false calls.

    - In ``big_spec`` each window holds at most one eligible call, so gene clusters equal row
      clusters and mutant M12 (the window CI by row) survives. Each second call sits 44-62 bp
      past the stop. Every such window has >= 70 bp of room, so the call stays eligible. It
      overlaps the window's truth terminator when there is one.
    - Matched eligible calls get confidence in [0.6, 1). One unmatched call in 25 gets
      [0.6, 0.9), and the rest get [0, 0.55). The 0.90 rule then ends near 0.6 with tp ~ 100
      and fp > 0, so neither the Wilson CI nor the window CI is the degenerate [1, 1].
    """
    spec = big_spec(seed)
    rnd = random.Random(seed + 1)
    by_end = {(g.end, g.strand): g for g in spec.genes}
    by_start = {(g.start, g.strand): g for g in spec.genes}
    extra: list[Call] = []
    k = 0
    for c in spec.calls:
        g = by_end.get((c.start - 8, "+")) if c.strand == "+" else by_start.get((c.end + 8, "-"))
        if g is None:
            continue
        k += 1
        if k % 2:
            continue
        if c.strand == "+":
            extra.append(Call(g.end + 44, g.end + 62, "+", 0.0))
        else:
            extra.append(Call(g.start - 62, g.start - 44, "-", 0.0))
    rep = Replicon(CONTIG, spec.length, "circular")
    windows = downstream_windows([gene_obj(g) for g in spec.genes], rep)
    truth = [iv(t.start, t.end, t.strand) for t in spec.truth]
    calls, n_un = [], 0
    for c in sorted([*spec.calls, *extra], key=lambda c: (c.start, c.end, c.strand)):
        x = iv(c.start, c.end, c.strand)
        u = rnd.random()
        if window_of(x, windows, rep) is None:
            conf = u
        elif terminator_match(x, truth, rep) is not None:
            conf = 0.6 + 0.4 * u
        else:
            n_un += 1
            conf = 0.6 + 0.3 * u if n_un % 25 == 0 else 0.55 * u
        calls.append(Call(c.start, c.end, c.strand, round(conf, 4)))
    return BigSpec(
        spec.genes, spec.labels, spec.probs, calls, spec.truth, spec.directon_of, spec.length
    )


def two_call_world(root: Path) -> tuple[World, BigSpec]:
    spec = two_call_spec()
    rep = Replicon(CONTIG, spec.length, "circular")
    w = build_world(root, spec.genes, [rep], spec.labels, spec.probs, spec.calls, spec.truth)
    return w, spec


def labelled(w: World, spec: BigSpec, *, with_p_only: bool = True) -> pd.DataFrame:
    by = {g.locus_tag: g for g in w.genes}
    rep = w.replicons[CONTIG]
    rows = [
        (
            int(lab == "positive"),
            spec.probs.get((a, b_), math.nan),
            float(intergenic_distance(by[a], by[b_], rep)),
            spec.directon_of[a],
        )
        for a, b_, lab in w.pairs
        if lab != "ambiguous" and ((a, b_) in spec.probs or not with_p_only)
    ]
    return pd.DataFrame(rows, columns=["y", "p", "d", "dir"])


def oracle(
    w: World,
    spec: BigSpec,
    *,
    seed: int = SEED,
    b: int = B,
    r: int = R,
    rows: frozenset[str] = frozenset(),
) -> dict[str, Any]:
    """Every post-stream-8 output, recomputed from part A in the accepted order.

    ``rows`` names the draws to cluster by row instead (``"q"``, ``"cv"`` or ``"window"``). It
    is used only for teeth, to show the world tells the two apart. The window CI is the whole
    ``fixed_cutoff_ci`` call, so its unused recall block is drawn too. ``rng_after`` is the
    generator after the last draw.
    """
    L = labelled(w, spec)
    p, y, d, cl = L.p.to_numpy(), L.y.to_numpy(), L.d.to_numpy(), L["dir"].to_numpy()
    rid = np.arange(len(L))
    t_u = choose_cutoff(threshold_curve(p, y)).threshold
    rng = for_step("m6_12_cutoffs", seed)
    for c in (cl, rid):  # streams 1-8, only to advance rng
        auroc_ci(p, y, c, rng=rng, B=b)
        fixed_cutoff_ci(p, y, c, t_u, rng=rng, B=b)
        auroc_ci(-d, y, c, rng=rng, B=b)
        auroc_diff_ci(p, -d, y, c, rng=rng, B=b)
        distance_fit_ci(d, y, c, rng=rng, B=b)
    q_cl = rid if "q" in rows else cl
    out: dict[str, Any] = {"q": list(threshold_boot_q(p, y, q_cl, rng=rng, B=b))}

    rep = w.replicons[CONTIG]
    windows = downstream_windows(w.genes, rep)
    truth = [iv(t.start, t.end, t.strand) for t in spec.truth]
    calls = sorted(spec.calls, key=lambda c: (c.start, c.end, c.strand))
    el = [(iv(c.start, c.end, c.strand), c.conf) for c in calls]
    el = [(x, q) for x, q in el if window_of(x, windows, rep) is not None]
    ivs, conf = [x for x, _ in el], np.array([q for _, q in el])
    gene = np.array([window_of(x, windows, rep) for x in ivs])
    ym = np.array([int(terminator_match(x, truth, rep) is not None) for x in ivs])
    cut = choose_cutoff(threshold_curve(conf, ym))
    rates = coincidence_rate(
        ivs, conf, windows, truth, rep, cut.threshold, rng=rng, R=r, return_rates=True
    )
    c_hat = float(np.mean(rates))
    out |= {
        "rates": rates,
        "c": c_hat,
        "se": float(np.std(rates, ddof=1) / math.sqrt(r)),
        "cv": precision_cv(p, y, rid if "cv" in rows else cl, rng=rng),
    }
    unit = np.arange(len(ivs)) if "window" in rows else gene
    win, win_recall = fixed_cutoff_ci(conf, ym, unit, cut.threshold, rng=rng, B=b)
    wilson_lo, wilson_hi = wilson_ci(cut.tp, cut.tp + cut.fp)
    out |= {
        "cut": cut,
        "ivs": ivs,
        "gene": gene,
        "ym": ym,
        "conf": conf,
        "window_ci": ci(win),
        "window_recall_ci": ci(win_recall),
        "lb": precision_lower_bound(cut.precision, c_hat),
        "lb_wilson": precision_lower_bound(wilson_lo, c_hat),
        "lb_window": precision_lower_bound(win.ci_low, c_hat),
        "wilson_hi": wilson_hi,
        "window_hi": win.ci_high,
        "rng_after": rng,
    }
    return out


@pytest.fixture(scope="module")
def two(tmp_path_factory: pytest.TempPathFactory) -> tuple[World, BigSpec, dict[str, Any], Any]:
    root = tmp_path_factory.mktemp("two")
    w, spec = two_call_world(root / "w")
    doc = cc.run(*w.run_args(root / "out"), pins=cc.Pins(**w.pins), seed=SEED, B=B, R=R)
    return w, spec, doc, oracle(w, spec)


def test_two_call_world_discriminates(two: Any) -> None:
    """Preconditions that make the oracle tests below kill the reviewers' mutants."""
    _w, _spec, doc, o = two
    tn = doc["terminatornet"]
    assert tn["calibrated"] is True
    assert o["cut"].fp > 0 and o["cut"].tp >= 30
    counts = pd.Series(o["gene"]).value_counts()
    assert (counts >= 2).sum() >= 10  # gene clusters != row clusters (M12)
    assert 0.0 < o["c"] < 1.0 and o["se"] > 0.0
    # lower vs upper limits differ after the correction (M10 / M11)
    assert precision_lower_bound(o["wilson_hi"], o["c"]) != pytest.approx(o["lb_wilson"])
    assert precision_lower_bound(o["window_hi"], o["c"]) != pytest.approx(o["lb_window"])
    assert o["lb_wilson"] > 0.0 and o["lb_window"] > 0.0  # not clamped to 0
    assert o["window_recall_ci"] != pytest.approx(o["window_ci"])  # recall half != precision
    # c over every eligible call (threshold 0) is another number: fewer calls are kept at t*
    # (so the exact comparison in the stream-10 test catches it), and c halved is far off
    n_kept = int((o["conf"] >= o["cut"].threshold).sum())
    assert 0 < n_kept < len(o["conf"])
    assert o["c"] / 2 < o["c"] - 5 * o["se"]


@pytest.mark.parametrize("which", ["q", "cv", "window"])
def test_row_clusters_would_differ(two: Any, which: str) -> None:
    """Teeth: clustering by row at the same rng state changes ``threshold_boot_q``,
    ``precision_cv`` or the window CI, so the oracle tests detect the wrong cluster unit."""
    w, spec, _doc, o = two
    alt = oracle(w, spec, rows=frozenset({which}))
    key = {"q": "q", "cv": "cv", "window": "window_ci"}[which]
    assert alt[key] != pytest.approx(o[key])


def test_stream9_threshold_boot_q(two: Any) -> None:
    """Stream 9 continues the shared rng and clusters by directon (M15; rows)."""
    _w, _spec, doc, o = two
    q = doc["uniop"]["threshold_boot_q"]
    assert len(q) == 3
    for g, want in zip(q, o["q"], strict=True):
        approx(g, want)


def test_stream10_coincidence_one_call_mean_and_sd_over_root_r(two: Any) -> None:
    """c = mean of the R per-placement rates of one call over the calls >= t*; MC SE =
    sd(ddof=1)/sqrt(R). The 20-call batch means (M8 / M9), c halved, c over every eligible
    call and a wrong SE denominator all give other numbers."""
    _w, _spec, doc, o = two
    tn = doc["terminatornet"]
    assert doc["R"] == R
    approx(tn["coincidence_rate"], o["c"])
    approx(tn["coincidence_rate_mc_se"], o["se"])


def test_stream11_precision_cv(two: Any) -> None:
    """``precision_cv`` folds directons, after the one coincidence call (M16; rows)."""
    _w, _spec, doc, o = two
    approx(doc["uniop"]["precision_cv"], o["cv"])


def test_stream12_window_ci_by_gene(two: Any) -> None:
    """The window CI is the *precision* half of ``fixed_cutoff_ci``, resampling genes (windows),
    not rows (M12)."""
    _w, _spec, doc, o = two
    approx_ci(doc["terminatornet"]["precision_ci_window"], o["window_ci"])


def test_draws_end_after_the_window_recall_block(
    two: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The driver's one generator ends exactly where the oracle's does. This pins the unused
    recall block of the window ``fixed_cutoff_ci`` and every other draw count, even where no
    written number shows it."""
    w, _spec, _doc, o = two
    seen: list[np.random.Generator] = []

    def spy(step: str, seed: int) -> np.random.Generator:
        g = for_step(step, seed)
        seen.append(g)
        return g

    monkeypatch.setattr(cc, "for_step", spy)
    run(w, tmp_path / "out", seed=SEED)
    assert len(seen) == 1
    assert seen[0].bit_generator.state == o["rng_after"].bit_generator.state


def test_lower_bounds_against_the_oracle(two: Any) -> None:
    _w, _spec, doc, o = two
    tn = doc["terminatornet"]
    approx(tn["precision_lower_bound"], o["lb"])
    approx(tn["precision_lower_bound_wilson"], o["lb_wilson"])
    approx(tn["precision_lower_bound_window"], o["lb_window"])


def test_lower_bounds_use_the_low_limits_of_the_written_cis(two: Any) -> None:
    """Self-consistency of the written JSON, independent of the draws: each bound is
    (lo - c)/(1 - c) of its own CI's *low* limit (M10: Wilson high, M11: window high)."""
    _w, _spec, doc, _o = two
    tn = doc["terminatornet"]
    c = tn["coincidence_rate"]
    assert tn["precision_lower_bound"] == pytest.approx(precision_lower_bound(tn["precision"], c))
    assert tn["precision_lower_bound_wilson"] == pytest.approx(
        precision_lower_bound(tn["precision_ci"][0], c)
    )
    assert tn["precision_lower_bound_window"] == pytest.approx(
        precision_lower_bound(tn["precision_ci_window"][0], c)
    )
    assert tn["precision_lower_bound_wilson"] <= tn["precision_lower_bound"]


def test_mc_se_finite_when_r_is_not_a_multiple_of_20(tmp_path: Path) -> None:
    """One call at R = 7 has an SE (sd/sqrt(7)); the batch design gave one batch and null."""
    w, spec = two_call_world(tmp_path / "w")
    doc = run(w, tmp_path / "out", R=7)
    o = oracle(w, spec, r=7)
    approx(doc["terminatornet"]["coincidence_rate_mc_se"], o["se"])
    assert doc["terminatornet"]["coincidence_rate_mc_se"] is not None


# ================================================================== schema of the new keys


NEW_TN_KEYS = {
    "coincidence_rate",
    "coincidence_rate_mc_se",
    "precision_ci_window",
    "precision_lower_bound",
    "precision_lower_bound_wilson",
    "precision_lower_bound_window",
}


def test_schema_new_terminatornet_keys(two: Any) -> None:
    _w, _spec, doc, _o = two
    tn = doc["terminatornet"]
    assert set(tn) >= NEW_TN_KEYS
    assert isinstance(tn["precision_ci_window"], list) and len(tn["precision_ci_window"]) == 2
    assert all(isinstance(x, float) for x in tn["precision_ci_window"])
    for k in NEW_TN_KEYS - {"precision_ci_window"}:
        assert isinstance(tn[k], float), k
    assert set(doc["pins"]) >= {"terms_sha256"}
    assert {"R", "seed", "B"} <= set(doc)


def test_uncalibrated_window_ci_is_a_null_pair(tmp_path: Path) -> None:
    """Uncalibrated, ``precision_ci_window`` is [null, null] like every other CI key (nit 8)."""
    calls = [EX6_CALLS[3]]  # c4: in no + window
    w = build_world(tmp_path / "w", EX6_GENES, [LIN10K], EX6_LABELS, EX6_PROBS, calls, EX6_TRUTH)
    tn = run(w, tmp_path / "out", min_labelled=2)["terminatornet"]
    assert tn["calibrated"] is False
    assert set(tn) >= NEW_TN_KEYS
    assert tn["precision_ci_window"] == [None, None]
    for k in NEW_TN_KEYS - {"precision_ci_window"}:
        assert tn[k] is None, k


# ================================================================== B2: terms.tsv.gz pinned


def test_k12_pins_hold_the_terms_sha256() -> None:
    assert cc.K12_PINS.terms_sha256 == K12_TERMS_SHA


def test_terms_pin_is_the_committed_file_and_its_provenance() -> None:
    assert sha256(K12_TERMS) == K12_TERMS_SHA
    prov = (K12_TERMS.parent / "PROVENANCE.md").read_text(encoding="utf-8")
    listed = re.findall(r"^([0-9a-f]{64})  terms\.tsv\.gz$", prov, flags=re.M)
    assert listed == [cc.K12_PINS.terms_sha256]


def test_world_pins_carry_the_terms_sha256(tmp_path: Path) -> None:
    w = example6(tmp_path / "w")
    assert w.pins["terms_sha256"] == sha256(w.terminatornet)
    assert cc.Pins(**w.pins).terms_sha256 == sha256(w.terminatornet)


def test_other_terms_file_with_the_right_provenance_is_refused(tmp_path: Path) -> None:
    """Swap in different calls; PROVENANCE.md still names the pinned commit."""
    w = example6(tmp_path / "w")
    pins = cc.Pins(**w.pins)
    prov = (w.terminatornet.parent / "PROVENANCE.md").read_text()
    with gzip.open(w.terminatornet, "rt") as fh:
        lines = fh.read().splitlines()
    with gzip.open(w.terminatornet, "wt") as fh:
        fh.write("\n".join(lines[:-2] + lines[-1:]) + "\n")  # drop one call
    assert (w.terminatornet.parent / "PROVENANCE.md").read_text() == prov
    assert cc.terminatornet_commit(w.terminatornet) == pins.terminatornet_commit
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="terms"):
        cc.run(*w.run_args(out), pins=pins, B=B, R=R, min_labelled=2)
    assert not out.exists() or not any(out.iterdir())


def test_terms_sha_written_is_the_pin(tmp_path: Path) -> None:
    w = example6(tmp_path / "w")
    doc = run(w, tmp_path / "out", min_labelled=2)
    assert doc["terminatornet"]["terms_sha256"] == w.pins["terms_sha256"]
    assert doc["pins"]["terms_sha256"] == w.pins["terms_sha256"]


# ================================================================== contig filters (M1, M19)


def _add_truth_rows(w: World, rows: list[tuple[str, int, int, str, str]]) -> None:
    path = w.regulondb_dir / "truth" / "terminator_truth.tsv"
    text = path.read_text(encoding="utf-8")
    text += "".join(f"{c}\t{s}\t{e}\t{st}\t{i}\n" for c, s, e, st, i in rows)
    path.write_text(text, encoding="utf-8")


def test_off_contig_truth_rows_are_ignored(tmp_path: Path) -> None:
    """M1. A pX row on top of g2's + window, near c2 (gap 0) and at c5's strand-matched spot.
    Unfiltered, it becomes a U00096.3 interval: n_truth goes from 4 to 6, and c2 / c5 match."""
    w = example6(tmp_path / "w")
    _add_truth_rows(w, [("pX", 2460, 2490, "+", "TX1"), ("pX", 4100, 4130, "-", "TX2")])
    out = tmp_path / "out"
    tn = run(w, out, min_labelled=2)["terminatornet"]
    assert (tn["n_truth"], tn["n_truth_eligible"]) == (4, 3)
    assert (tn["threshold"], tn["tp"], tn["fp"]) == (0.95, 1, 0)
    assert tn["recall_all"] == pytest.approx(1 / 4)
    calls = pd.read_csv(out / "terminatornet_calls.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert list(calls["matched_source_id"]) == ["T1", "", "T3", ""]


def test_off_contig_genes_do_not_change_windows(tmp_path: Path) -> None:
    """M19. A gene on pX at the coordinates between c1 and g2 would cut g1's window at 1300 if
    the contig filter were missing anywhere (windows, eligibility, directons)."""
    plasmid = Replicon("pX", 10_000, "linear")
    genes = [*EX6_GENES, G("px1", 1200, 1300, "+", "pX"), G("px2", 2550, 2700, "-", "pX")]
    calls = [*EX6_CALLS, Call(1350, 1380, "+", 0.97)]  # in g1's window only if px1 is ignored
    truth = [*EX6_TRUTH, Term("T5", 1390, 1400, "+")]
    w = build_world(tmp_path / "w", genes, [LIN10K, plasmid], EX6_LABELS, EX6_PROBS, calls, truth)
    doc = run(w, tmp_path / "out", min_labelled=2)
    tn = doc["terminatornet"]
    assert (tn["n_calls"], tn["n_eligible"], tn["n_truth"], tn["n_truth_eligible"]) == (7, 5, 5, 4)
    assert (tn["tp"], tn["fp"]) == (2, 0)
    assert doc["n_directons"] == 2


# ================================================================== recall_all vs recall (M23)


def test_recall_all_counts_a_hit_on_an_ineligible_terminator(tmp_path: Path) -> None:
    """Call cX [1440, 1470) + (0.97) is in g1's window. Truth TX [1490, 1510) + crosses g2's start,
    so it is ineligible, but it lies 20 bp from cX. At t* = 0.95, the called calls are cX and c1.
    Eligible truth is T1, T2, T3, of which T1 is hit: recall = 1/3. All truth is 5 rows, of which
    T1 and TX are hit: recall_all = 2/5. The eligible hit count would give 1/5."""
    calls = [*EX6_CALLS, Call(1440, 1470, "+", 0.97)]
    truth = [*EX6_TRUTH, Term("TX", 1490, 1510, "+")]
    w = build_world(tmp_path / "w", EX6_GENES, [LIN10K], EX6_LABELS, EX6_PROBS, calls, truth)
    tn = run(w, tmp_path / "out", min_labelled=2)["terminatornet"]
    assert (tn["n_truth"], tn["n_truth_eligible"]) == (5, 3)
    assert (tn["threshold"], tn["tp"], tn["fp"]) == (0.95, 2, 0)
    assert tn["recall"] == pytest.approx(1 / 3)
    assert tn["recall_all"] == pytest.approx(2 / 5)
    approx_ci(tn["recall_all_ci"], list(wilson_ci(2, 5)))


# ================================================================== the 1 % boundary (M6)


def _boundary_world(root: Path, n_drop: int) -> tuple[World, int]:
    """``big_spec`` with labels trimmed to 500 labelled pairs, then ``n_drop`` UniOP p removed."""
    spec = big_spec()
    labels = dict(spec.labels)
    labelled = [k for k, v in labels.items() if v != "ambiguous" and k in spec.probs]
    for k in labelled[500:]:
        labels[k] = "ambiguous"
    rep = Replicon(CONTIG, spec.length, "circular")
    w = build_world(root, spec.genes, [rep], labels, spec.probs, spec.calls, spec.truth)
    n_lab = sum(1 for p in w.pairs if p[2] != "ambiguous")
    gone = set(labelled[10 : 10 + n_drop])
    probs = {k: v for k, v in spec.probs.items() if k not in gone}
    for f in ("uniop.pred", "uniop_cds.gff", "uniop_map.tsv"):
        (w.uniop_dir / f).unlink()
    write_uniop(w.uniop_dir, w.genes, w.replicons, probs)
    refresh_pins(w)
    return w, n_lab


def test_exactly_one_percent_without_p_passes(tmp_path: Path) -> None:
    """500 labelled pairs, 5 without p (exactly 1 %): excluded and counted, not an escalation."""
    w, n_lab = _boundary_world(tmp_path / "w", 5)
    assert n_lab == 500
    u = run(w, tmp_path / "out", min_labelled=400)["uniop"]
    assert (u["n_without_p"], u["n_labelled"]) == (5, 495)


def test_one_percent_plus_one_without_p_escalates(tmp_path: Path) -> None:
    w, n_lab = _boundary_world(tmp_path / "w", 6)
    assert n_lab == 500
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="escalate"):
        run(w, out, min_labelled=400)
    assert no_outputs(out)


# ================================================================== guards write nothing


@pytest.mark.parametrize("label", ["negative", "positive"])
def test_one_class_raises_and_writes_nothing(tmp_path: Path, label: str) -> None:
    labels = {k: label for k in EX6_LABELS}
    w = build_world(tmp_path / "w", EX6_GENES, [LIN10K], labels, EX6_PROBS, EX6_CALLS, EX6_TRUTH)
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run(w, out, min_labelled=2)
    assert not out.exists()


# ================================================================== main: missing inputs


@pytest.mark.parametrize("flag", ["--uniop-dir", "--bakta-json", "--genome", "--terminatornet"])
def test_main_missing_input_exits_1_without_a_traceback(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    flag: str,
) -> None:
    w, _ = big_world(tmp_path / "w")  # >= 500 labelled pairs: no guard before the read
    monkeypatch.setattr(cc, "K12_PINS", cc.Pins(**w.pins))
    out = tmp_path / "out"
    args = w.cli_args(out)
    i = args.index(flag)
    target = Path(args[i + 1])
    missing = target.parent / ("missing_" + target.name)
    args[i + 1] = str(missing)
    if flag == "--terminatornet":  # PROVENANCE.md is still beside the missing file
        assert (missing.parent / "PROVENANCE.md").is_file()
    assert cc.main([*args, "--B", "20"]) == 1  # returns, does not raise
    err = capsys.readouterr().err
    assert err.startswith("calibrate_cutoffs: ")
    assert "Traceback" not in err
    assert no_outputs(out)


# ================================================================== statistician asks (#264)


def test_window_bound_present_and_non_null_when_calibrated(two: Any) -> None:
    """The ADR quotes ``precision_lower_bound_window``. When calibrated, it is a finite float
    in [0, 1] and lies at or below the point bound."""
    _w, _spec, doc, _o = two
    tn = doc["terminatornet"]
    assert "precision_lower_bound_window" in tn
    lb = tn["precision_lower_bound_window"]
    assert isinstance(lb, float) and 0.0 <= lb <= 1.0
    assert lb <= tn["precision_lower_bound"]


def test_refit_population_is_l_when_a_pair_lacks_p(tmp_path: Path) -> None:
    """d0 and s are fitted on L (labelled pairs with a UniOP p). One labelled pair without p is
    left out, and that changes the fit (teeth)."""
    spec = big_spec()
    labelled_keys = [k for k, v in spec.labels.items() if v != "ambiguous" and k in spec.probs]
    gone = labelled_keys[5]
    rep = Replicon(CONTIG, spec.length, "circular")
    w = build_world(
        tmp_path / "w", spec.genes, [rep], spec.labels, spec.probs, spec.calls, spec.truth
    )
    probs = {k: v for k, v in spec.probs.items() if k != gone}
    for f in ("uniop.pred", "uniop_cds.gff", "uniop_map.tsv"):
        (w.uniop_dir / f).unlink()
    write_uniop(w.uniop_dir, w.genes, w.replicons, probs)
    refresh_pins(w)
    doc = run(w, tmp_path / "out")
    assert doc["uniop"]["n_without_p"] == 1
    trimmed = BigSpec(
        spec.genes, spec.labels, probs, spec.calls, spec.truth, spec.directon_of, spec.length
    )
    L = labelled(w, trimmed)
    every = labelled(w, trimmed, with_p_only=False)
    assert len(every) == len(L) + 1
    fit, fit_all = fit_distance_model(L.d, L.y), fit_distance_model(every.d, every.y)
    dm = doc["distance_model"]
    assert dm["converged"] is True
    assert (dm["d0"], dm["s"]) == pytest.approx((fit.d0, fit.s), rel=1e-9)
    assert (fit.d0, fit.s) != pytest.approx((fit_all.d0, fit_all.s), rel=1e-9)


def test_off_contig_truth_rows_are_dropped_with_a_logged_count(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """The spec asks for a logged count. Two pX truth rows give one INFO record from the
    ``calibrate_cutoffs`` logger that names terminator truth and the count 2."""
    w = example6(tmp_path / "w")
    _add_truth_rows(w, [("pX", 2460, 2490, "+", "TX1"), ("pX", 4100, 4130, "-", "TX2")])
    with caplog.at_level(logging.INFO, logger="calibrate_cutoffs"):
        run(w, tmp_path / "out", min_labelled=2)
    msgs = [
        r.getMessage()
        for r in caplog.records
        if r.name == "calibrate_cutoffs" and r.levelno == logging.INFO
    ]
    hits = [m for m in msgs if "truth" in m.lower() and re.search(r"\b2\b", m)]
    assert hits, msgs
