r"""K-12 cut-off calibration: UniOP AUROC and cut-off, the d0 / s refit, the TerminatorNet cut-off.

Issue M6-12 part B (GIAP_SPEC.md §15.4 steps 3-5, §7.1.3, §7.4.6); spec
``docs/specs/M6-12-uniop-terminator-cutoffs.md``; methods ``docs/methods/operons.md`` "K-12 cut-off
calibration (M6-12)". This driver does the I/O, the identity checks, the coordinate re-map,
directons, eligibility and matching; every statistic is a pure function of ``giap.stats.cutoffs``.
Randomness comes from one ``giap.rng.for_step(STEP, seed)`` generator, drawn in the methods'
fixed stream order. Nothing written here is committed (R6 / R11).

Usage::

    python validation/calibrate_cutoffs.py \
        --regulondb-dir $GIAP_DATA_DIR/regulondb/14.5.0 \
        --bakta-json tests/data/recorded/bakta/1.12.2-full/k12/k12_mg1655.json.gz \
        --genome tests/data/genomes/k12_mg1655.fna.gz \
        --uniop-dir $GIAP_ROOT/recordings/uniop/k12_mg1655 \
        --terminatornet tests/data/recorded/terminatornet/721803095459/k12/terms.tsv.gz \
        --out-dir $GIAP_DATA_DIR/regulondb/14.5.0/calibration
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import re
import sys
import tempfile
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from validation.regulondb_truth import (
    DATASETS,
    parse_terminator_set,
    verified_release,
)

from giap.annotate.bakta import parse_bakta_json
from giap.annotate.operons import adjacent_pairs, intergenic_distance
from giap.annotate.terminators import parse_terminatornet
from giap.annotate.uniop_adapter import PairProb, genome_order, parse_uniop_pred
from giap.io.fasta import read_genome_fasta, write_fasta
from giap.manifest import finite
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene
from giap.rng import for_step
from giap.stats import cutoffs as sc

STEP = sc.STEP
UNIOP_PIN = "3ed14dfc2a5d69b993d99a26e50c1f8e9723c087"  # ADR-0008
TOL_BP = 50
MAX_WITHOUT_P = 0.01
_COMMIT = re.compile(r"tool commit `([0-9a-f]{40})`")
log = logging.getLogger("calibrate_cutoffs")


@dataclass(frozen=True, slots=True)
class Pins:
    """The expected recording identities (spec "Inputs, as they exist now")."""

    terminatornet_commit: str
    genome_sha256: str
    uniop_pred_sha256: str
    uniop_gff_sha256: str
    terms_sha256: str


K12_PINS = Pins(
    terminatornet_commit="72180309545988dba63b8bd5db5d0b09bf4dba70",
    genome_sha256="918f35c7bfc8a3c8d824ed281a4185903025a0361242219b6b0e0ea6ecaa14f2",
    uniop_pred_sha256="7d852a93a436a9bcef2cf4cf38d608bcd5d59054d98bd1b1ddc21196690b5ee0",
    uniop_gff_sha256="ea5870894e26c649b4e7520c65cf57e1138e1786688db7663a7b3e3eea062906",
    terms_sha256="9be61c75c0d7835f7d8c8f3c3948b95adfa34445003ec570ff9bbb2ca1e3a2ff",
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalized_genome_sha256(path: Path) -> str:
    """SHA-256 of the genome as preflight writes ``genome.norm.fa``."""
    with tempfile.TemporaryDirectory() as tmp:
        norm = Path(tmp) / "genome.norm.fa"
        write_fasta(norm, ((r.id, r.sequence) for r in read_genome_fasta(path)))
        return _sha(norm)


def terminatornet_commit(terms_path: Path) -> str:
    """The ``tool commit `<sha>``` recorded in the ``PROVENANCE.md`` beside ``terms_path``."""
    prov = terms_path.parent / "PROVENANCE.md"
    if not prov.is_file():
        raise ValueError(f"{prov}: no PROVENANCE.md beside the TerminatorNet recording")
    m = _COMMIT.search(prov.read_text(encoding="utf-8"))
    if m is None:
        raise ValueError(f"{prov}: PROVENANCE.md records no tool commit")
    return m.group(1)


def coordinate_id_map(gff: Path, genes: Sequence[Gene]) -> dict[str, str]:
    """GFF ID -> locus tag by exact coordinates; an unmatched or ambiguous row is a ValueError."""
    by: dict[tuple[int, int, str, str], list[str]] = {}
    for g in genes:
        iv = g.interval
        by.setdefault((iv.start + 1, iv.end, iv.strand, iv.contig), []).append(g.locus_tag)
    out: dict[str, str] = {}
    for ln in gff.read_text(encoding="utf-8").splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        f = ln.split("\t")
        gid = f[8].split(";")[0].removeprefix("ID=")
        hits = by.get((int(f[3]), int(f[4]), f[6], f[0]), [])
        if len(hits) != 1:
            raise ValueError(f"UniOP GFF row {gid}: {len(hits)} genes at {f[0]}:{f[3]}..{f[4]}")
        out[gid] = hits[0]
    return out


def directon_ids(genes: Sequence[Gene], replicon: Replicon) -> dict[str, int]:
    """Locus tag -> directon (maximal same-strand run in genome order; circular wrap joined)."""
    gs = genome_order([g for g in genes if g.interval.contig == replicon.id], [replicon.id])
    gs_ = gs[replicon.id]
    run, out = -1, {}
    for k, g in enumerate(gs_):
        if k == 0 or g.interval.strand != gs_[k - 1].interval.strand:
            run += 1
        out[g.locus_tag] = run
    ring = replicon.topology == "circular" and run > 0
    if ring and gs_[0].interval.strand == gs_[-1].interval.strand:
        out = {t: 0 if r == run else r for t, r in out.items()}
    return out


def pair_table(
    genes: Sequence[Gene], replicon: Replicon, truth: pd.DataFrame, probs: Sequence[PairProb]
) -> pd.DataFrame:
    """Labelled truth rows in file order, with ``y``, ``uniop_p``, ``d`` and ``directon_id``."""
    on = [g for g in genes if g.interval.contig == replicon.id]
    by = {g.locus_tag: g for g in on}
    adjacent = {(a.locus_tag, b.locus_tag) for a, b in adjacent_pairs(on, {replicon.id: replicon})}
    p = {(x.gene_a, x.gene_b): x.uniop_p for x in probs}
    directon = directon_ids(on, replicon)
    rows = []
    for a, b, label in zip(truth["gene_a"], truth["gene_b"], truth["label"], strict=True):
        if a not in by or b not in by or (a, b) not in adjacent:
            raise ValueError(f"truth pair {a}-{b} is not an adjacent same-strand pair of the model")
        if label in ("positive", "negative"):
            d = intergenic_distance(by[a], by[b], replicon)
            y = int(label == "positive")
            rows.append((a, b, label, y, p.get((a, b), math.nan), d, directon[a]))
    cols = ["gene_a", "gene_b", "label", "y", "uniop_p", "d", "directon_id"]
    return pd.DataFrame(rows, columns=cols)


def _hits(called: Sequence[Interval], truth: Sequence[Interval], rep: Replicon, tol: int) -> int:
    """Truth intervals within ``tol`` bp (interval gap, same contig and strand) of any call."""
    return sum(
        any(
            c.contig == t.contig and c.strand == t.strand and sc.interval_gap(c, t, rep) <= tol
            for c in called
        )
        for t in truth
    )


def truth_recall(
    calls: Sequence[Interval],
    confidence: Sequence[float],
    truth: Sequence[Interval],
    windows: Sequence[sc.Window],
    replicon: Replicon,
    threshold: float,
    *,
    tol_bp: int = TOL_BP,
) -> tuple[int, int]:
    """(k, n): n eligible truth terminators, k of them near an eligible call >= threshold."""
    called = _called(calls, confidence, windows, replicon, threshold)
    eligible = [t for t in truth if sc.window_of(t, windows, replicon) is not None]
    return _hits(called, eligible, replicon, tol_bp), len(eligible)


def _called(
    calls: Sequence[Interval],
    conf: Sequence[float],
    windows: Sequence[sc.Window],
    rep: Replicon,
    t: float,
) -> list[Interval]:
    return [
        c
        for c, p in zip(calls, conf, strict=True)
        if p >= t and sc.window_of(c, windows, rep) is not None
    ]


# ------------------------------------------------------------------------------- inputs


def _check(what: str, got: str, want: str) -> None:
    if got != want:
        raise ValueError(f"{what}: {got} != expected {want}")


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _verify(
    regulondb_dir: Path, genome: Path, uniop_dir: Path, terminatornet: Path, pins: Pins
) -> tuple[dict[str, str], dict[str, Any]]:
    """Every identity check of the degenerate table; returns the release tables and summary."""
    tables, _meta = verified_release(regulondb_dir)
    summary = json.loads((regulondb_dir / "truth" / "truth_summary.json").read_text())
    for ds in DATASETS:
        got = hashlib.sha256(tables[ds].encode("utf-8")).hexdigest()
        _check(f"{ds} sha256 vs truth_summary.json", got, str(summary["sha256"].get(ds)))
    _check("TerminatorNet commit", terminatornet_commit(terminatornet), pins.terminatornet_commit)
    _check("terms.tsv.gz sha256", _sha(terminatornet), pins.terms_sha256)
    _check("genome.norm.fa sha256", normalized_genome_sha256(genome), pins.genome_sha256)
    _check("uniop.pred sha256", _sha(uniop_dir / "uniop.pred"), pins.uniop_pred_sha256)
    _check("uniop_cds.gff sha256", _sha(uniop_dir / "uniop_cds.gff"), pins.uniop_gff_sha256)
    return tables, summary


# ------------------------------------------------------------------------------- statistics


def _ci(b: Any) -> list[float]:
    return [b.ci_low, b.ci_high]


def _uniop_streams(L: pd.DataFrame, t: float, rng: np.random.Generator, B: int) -> dict[str, Any]:
    """Streams 1-8 (methods "Bootstrap streams"): directon clusters, then row_id clusters."""
    p, y, d = L["uniop_p"].to_numpy(), L["y"].to_numpy(), L["d"].to_numpy(dtype=float)
    out: dict[str, Any] = {}
    for suffix, cl in (("", L["directon_id"].to_numpy()), ("_iid", np.arange(len(L)))):
        out["auroc_ci" + suffix] = _ci(sc.auroc_ci(p, y, cl, rng=rng, B=B))
        prec, rec = sc.fixed_cutoff_ci(p, y, cl, t, rng=rng, B=B)
        out["precision_ci" + suffix], out["recall_ci" + suffix] = _ci(prec), _ci(rec)
        out["d_auroc_ci" + suffix] = _ci(sc.auroc_ci(-d, y, cl, rng=rng, B=B))
        out["d_auroc_diff_ci" + suffix] = _ci(sc.auroc_diff_ci(p, -d, y, cl, rng=rng, B=B))
        d0, s = sc.distance_fit_ci(d, y, cl, rng=rng, B=B)
        out["d0_ci" + suffix], out["s_ci" + suffix] = _ci(d0), _ci(s)
    return out


def _cut_fields(cut: sc.Cutoff) -> dict[str, Any]:
    keys = ("threshold", "literal_threshold", "rule", "precision", "recall", "f1", "tp", "fp")
    return {k: getattr(cut, k) for k in keys}


def _coincidence(
    elig: Sequence[Interval],
    conf: Sequence[float],
    ctx: _TnCtx,
    t: float,
    rng: np.random.Generator,
    R: int,
) -> tuple[float, float]:
    """(c, Monte Carlo SE of c): one call; the mean and sd(ddof=1)/sqrt(R) of the R rates."""
    rates = np.asarray(
        sc.coincidence_rate(
            elig, conf, ctx.windows, ctx.truth, ctx.rep, t, rng=rng, R=R, return_rates=True
        )
    )
    se = float(np.std(rates, ddof=1) / math.sqrt(R)) if R > 1 else math.nan
    return float(np.mean(rates)), se


@dataclass(frozen=True)
class _TnCtx:
    rep: Replicon
    windows: list[sc.Window]
    truth: list[Interval]
    truth_ids: list[str]
    truth_types: list[str]


def _tn_point(
    calls: Sequence[tuple[Interval, float]], ctx: _TnCtx
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame | None, sc.Cutoff | None]:
    """Eligibility, matching, the cut-off and truth recall (no randomness)."""
    rows = []
    for c, p in sorted(calls, key=lambda x: (x[0].start, x[0].end, x[0].strand)):
        gene = sc.window_of(c, ctx.windows, ctx.rep)
        if gene is not None:
            j = sc.terminator_match(c, ctx.truth, ctx.rep, TOL_BP)
            mid = "" if j is None else ctx.truth_ids[j]
            rows.append((c.contig, c.start, c.end, c.strand, p, gene, mid, c))
    cols = ["contig", "start", "end", "strand", "confidence", "gene", "matched_source_id", "iv"]
    el = pd.DataFrame(rows, columns=cols)
    y = (el["matched_source_id"] != "").astype(int).to_numpy()
    n_truth = len(ctx.truth)
    elig_truth = [k for k, t in enumerate(ctx.truth) if sc.window_of(t, ctx.windows, ctx.rep)]
    doc: dict[str, Any] = {
        "n_calls": len(calls),
        "n_eligible": len(el),
        "n_truth": n_truth,
        "n_truth_eligible": len(elig_truth),
        "calibrated": bool(y.sum() > 0),
    }
    if not doc["calibrated"]:
        return doc, el, None, None
    conf = el["confidence"].to_numpy(dtype=float)
    curve = sc.threshold_curve(conf, y)
    cut = sc.choose_cutoff(curve)
    called = [c for c, p in zip(el["iv"], conf, strict=True) if p >= cut.threshold]
    k_el = _hits(called, [ctx.truth[k] for k in elig_truth], ctx.rep, TOL_BP)
    k_all = _hits(called, ctx.truth, ctx.rep, TOL_BP)
    n_el = len(elig_truth)
    by_type: dict[str, dict[str, Any]] = {}
    for kind in sorted({ctx.truth_types[k] for k in elig_truth}):
        ts = [ctx.truth[k] for k in elig_truth if ctx.truth_types[k] == kind]
        kk = _hits(called, ts, ctx.rep, TOL_BP)
        ci = list(sc.wilson_ci(kk, len(ts)))
        by_type[kind] = {"k": kk, "n": len(ts), "recall": kk / len(ts), "recall_ci": ci}
    doc |= _cut_fields(cut) | {
        "precision_ci": list(sc.wilson_ci(cut.tp, cut.tp + cut.fp)),
        "precision_is_lower_bound": True,
        "recall": k_el / n_el if n_el else math.nan,
        "recall_ci": list(sc.wilson_ci(k_el, n_el)),
        "recall_all": k_all / n_truth,
        "recall_all_ci": list(sc.wilson_ci(k_all, n_truth)),
        "recall_by_type": by_type,
        "auroc_diagnostic": sc.auroc(conf, y),
    }
    return doc, el, curve, cut


_TN_NULL = (
    "threshold", "literal_threshold", "rule", "precision", "precision_ci", "coincidence_rate",
    "coincidence_rate_mc_se", "precision_lower_bound", "precision_lower_bound_wilson",
    "precision_lower_bound_window", "precision_ci_window", "recall", "recall_ci", "recall_all",
    "recall_all_ci", "recall_by_type", "auroc_diagnostic", "tp", "fp",
)  # fmt: skip


def _curve_frame(curve: pd.DataFrame | None) -> pd.DataFrame:
    if curve is None:
        return pd.DataFrame(columns=["threshold", "tp", "fp", "precision", "recall", "f1"])
    P = int(curve["tp"].iloc[-1])
    tp, fp = curve["tp"], curve["fp"]
    rec = tp / P if P else math.nan
    return curve.assign(precision=tp / (tp + fp), recall=rec, f1=2 * tp / (tp + fp + P))


# ------------------------------------------------------------------------------- the run


def run(
    regulondb_dir: Path,
    bakta_json: Path,
    genome: Path,
    uniop_dir: Path,
    terminatornet: Path,
    out_dir: Path,
    *,
    seed: int = sc.DEFAULT_SEED,
    B: int = sc.DEFAULT_B,
    R: int = 200,
    pins: Pins = K12_PINS,
    contig: str = "U00096.3",
    min_labelled: int = sc.MIN_LABELLED,
) -> dict[str, Any]:
    """Verify, compute and write the four outputs; returns the ``cutoffs.json`` content."""
    tables, summary = _verify(regulondb_dir, genome, uniop_dir, terminatornet, pins)
    seqs = {r.id: r.sequence for r in read_genome_fasta(genome)}
    model = parse_bakta_json(bakta_json, genome=seqs).model
    rep = model.replicons[contig]
    genes = [g for g in model.genes if g.interval.contig == contig]
    gff, pred = uniop_dir / "uniop_cds.gff", uniop_dir / "uniop.pred"
    probs = parse_uniop_pred(pred, gff, coordinate_id_map(gff, model.genes))
    tab = pair_table(genes, rep, _read_tsv(regulondb_dir / "truth" / "pair_truth.tsv"), probs)
    has_p = tab["uniop_p"].notna()
    n_without_p = int((~has_p).sum())
    if n_without_p > MAX_WITHOUT_P * len(tab):
        raise ValueError(
            f"{n_without_p} of {len(tab)} labelled pairs lack a UniOP p (> 1 %): escalate"
        )
    L = tab[has_p].reset_index(drop=True)
    n_pos, n_neg = sc.require_labelled(L["y"], min_labelled)
    p, y, d = L["uniop_p"].to_numpy(), L["y"].to_numpy(), L["d"].to_numpy(dtype=float)
    u_curve = sc.threshold_curve(p, y)
    u_cut = sc.choose_cutoff(u_curve)
    fit = sc.fit_distance_model(d, y)

    tdf = _read_tsv(regulondb_dir / "truth" / "terminator_truth.tsv")
    if (off := int((tdf["contig"] != contig).sum())) > 0:
        log.info("dropped %d terminator truth rows on other contigs", off)
    tdf = tdf[tdf["contig"] == contig]
    cols = (tdf["start"], tdf["end"], tdf["strand"])
    rdb_terms = parse_terminator_set(tables["TerminatorSet"], contig)
    kinds = {t.terminator_id: t.terminator_type for t in rdb_terms}
    ctx = _TnCtx(
        rep,
        sc.downstream_windows(genes, rep),
        [Interval(contig, int(a), int(b), s) for a, b, s in zip(*cols, strict=True)],
        list(tdf["source_id"]),
        [kinds.get(i) or "unknown" for i in tdf["source_id"]],
    )
    terms = parse_terminatornet(terminatornet, model.replicons)
    all_calls = [(t.interval, t.confidence) for t in terms]
    calls = [(c, q) for c, q in all_calls if c.contig == contig]
    if len(calls) < len(all_calls):
        log.info("ignored %d TerminatorNet calls on other contigs", len(all_calls) - len(calls))
    tn, el, t_curve, t_cut = _tn_point(calls, ctx)

    rng = for_step(STEP, seed)
    boot = _uniop_streams(L, u_cut.threshold, rng, B)
    cl = L["directon_id"].to_numpy()
    q = sc.threshold_boot_q(p, y, cl, rng=rng, B=B)
    nulls = {k: [math.nan] * 2 if k.endswith(("_ci", "_ci_window")) else None for k in _TN_NULL}
    tn = nulls | tn
    if t_cut is not None:
        conf = el["confidence"].to_numpy(dtype=float)
        c, se = _coincidence(list(el["iv"]), list(conf), ctx, t_cut.threshold, rng, R)
        tn |= {"coincidence_rate": c, "coincidence_rate_mc_se": se}
    n_dir = int(L["directon_id"].nunique())
    cv = sc.precision_cv(p, y, cl, rng=rng) if n_dir >= 10 else math.nan
    if t_cut is not None:
        ym = (el["matched_source_id"] != "").astype(int).to_numpy()
        win, _ = sc.fixed_cutoff_ci(conf, ym, el["gene"].to_numpy(), t_cut.threshold, rng=rng, B=B)
        lb = sc.precision_lower_bound
        c = tn["coincidence_rate"]
        tn |= {
            "precision_ci_window": _ci(win),
            "precision_lower_bound": lb(t_cut.precision, c),
            "precision_lower_bound_wilson": lb(tn["precision_ci"][0], c),
            "precision_lower_bound_window": lb(win.ci_low, c),
        }

    uniop = {
        "pin": UNIOP_PIN,
        "pred_sha256": pins.uniop_pred_sha256,
        "gff_sha256": pins.uniop_gff_sha256,
        "n_labelled": len(L),
        "n_positive": n_pos,
        "n_negative": n_neg,
        "n_without_p": n_without_p,
        "auroc": sc.auroc(p, y),
        **_cut_fields(u_cut),
        "precision_cv": cv,
        "threshold_boot_q": list(q),
        **{k: v for k, v in boot.items() if not k.startswith(("d_", "d0", "s_"))},
    }
    distance = {
        "d0": fit.d0,
        "s": fit.s,
        "converged": fit.converged,
        "reason": fit.reason,
        "auroc": sc.auroc(-d, y),
        "auroc_diff": sc.auroc(p, y) - sc.auroc(-d, y),
        "logloss": sc.distance_log_loss(d, y, fit.d0, fit.s),
        "logloss_placeholder": sc.distance_log_loss(d, y, sc.D0_DEFAULT, sc.S_DEFAULT),
        **{k.removeprefix("d_"): v for k, v in boot.items() if k.startswith(("d_", "d0", "s_"))},
    }
    tn = {"commit": pins.terminatornet_commit, "terms_sha256": _sha(terminatornet), **tn}
    tn["precision_is_lower_bound"] = True
    doc = {
        "regulondb_version": summary.get("regulondb_version"),
        "release_date": summary.get("release_date"),
        "truth_sha256": summary["sha256"],
        "date": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "seed": seed,
        "B": B,
        "R": R,
        "cluster": "directon",
        "n_directons": n_dir,
        "pins": asdict(pins),
        "uniop": uniop,
        "distance_model": distance,
        "terminatornet": tn,
    }
    text = json.dumps(finite(doc), indent=2) + "\n"
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv: dict[str, Any] = {"sep": "\t", "index": False, "lineterminator": "\n"}
    _curve_frame(u_curve).to_csv(out_dir / "uniop_curve.tsv", **tsv)
    _curve_frame(t_curve).to_csv(out_dir / "terminatornet_curve.tsv", **tsv)
    el.drop(columns="iv").to_csv(out_dir / "terminatornet_calls.tsv", **tsv)
    (out_dir / "cutoffs.json").write_text(text, encoding="utf-8")
    result: dict[str, Any] = json.loads(text)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    for flag in ("--regulondb-dir", "--bakta-json", "--genome", "--uniop-dir", "--terminatornet"):
        ap.add_argument(flag, type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=sc.DEFAULT_SEED)
    ap.add_argument("--B", type=int, default=sc.DEFAULT_B)
    a = ap.parse_args(argv)
    paths = (a.regulondb_dir, a.bakta_json, a.genome, a.uniop_dir, a.terminatornet, a.out_dir)
    try:
        run(*paths, seed=a.seed, B=a.B, pins=K12_PINS)
    except (ValueError, FileNotFoundError) as exc:
        print(f"calibrate_cutoffs: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
