"""Inputs and an oracle for the M5-08 origin-to-terminus gradient (docs/specs/M5-08-oriter.md).

Nothing here reads ``gippa.context.oriter``: the oracle is written from the spec text (§7.3.5(e),
lead decisions 31-35, ADR-0021 decision 1 with Amendment 1) on top of the *merged* engines
``gippa.stats.autocorr.choose_block_len`` / ``block_len_prime`` (M5-02),
``gippa.stats.bootstrap.circular_block_bootstrap`` / ``percentile_ci`` (M4-04) and
``gippa.context.regions.replicon_index`` (M4-17), each tested by its own issue. The replichore
geometry is taken from the independent M5-04 helper ``oric_oracle.replichore``.

Inputs. Hand-built gene tables reuse the M4-17 builder (``region_tests_data.gene_table``, §9.3
``genes.tsv`` columns); :func:`tu_frame` derives the C5 ``tus`` frame (0-based half-open) and
:func:`features` the C15 ``features`` frame (``type = oric`` rows). :func:`recorded_features`
parses the committed Bakta 1.12.2 JSONs (``tests/data/recorded/bakta/1.12.2-full/``) with
``parse_bakta_json`` -- no tool run.

Interface assumed of ``gippa.context.oriter`` (spec "Public interface"):

``oric_position(features, contig, length) -> float | None``
``replichore(mid, oric, length) -> (labels, positions)``
``gene_replichores(genes, replicons, features) -> DataFrame[locus_tag, replichore,
    replichore_position]``
``oriter_rows(genes, tus, replicons, features, cfg, rng) -> DataFrame[AUTOCORR_COLUMNS]``

and the module constants ``STATISTIC``, ``METHOD``, ``NULL_MODEL``, ``ARMS``, mirrored below.

Draw order assumed of ``oriter_rows`` (spec "RNG", lead decision 35): replicons in C1 order;
per analysed replicon the block-length permutations first (only when
``context.islands.block_len`` is ``"auto"`` and the pooled fit exists), then the bootstrap index
replicates of the pooled, right and left fits, skipping a fit that does not exist.
"""

from __future__ import annotations

import gzip
import math
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import oric_oracle as oo
import region_tests_data as rd
from region_tests_data import Contig, flag_set, is_na  # noqa: F401  (re-exported)

REPO = Path(__file__).resolve().parents[2]
RECORDED = REPO / "tests" / "data" / "recorded" / "bakta" / "1.12.2-full"
GENOMES = REPO / "tests" / "data" / "genomes"

#: §6.6 random stream of the ori-ter part
STEP = "context_stats.oriter"
STATISTIC = "oriter_slope"
METHOD = "ols_replichore"
NULL_MODEL = "bootstrap"
#: fit order = draw order (spec "RNG")
ARMS = ("pooled", "right", "left")
DESCRIPTIVE = "descriptive"
NO_ORIC = "no_oric"
NOT_ESTIMABLE = "not_estimable"
MULTIPLE_ORIC_EVENT = "context_stats.oriter.multiple_oric"
#: lead decision 31: several oriC calls are averaged only inside this share of the replicon
ORIC_SPREAD_FRAC = 0.05
#: the §9.3 autocorrelation columns plus the §9.2 test columns (M5-05 ``AUTOCORR_COLUMNS``)
COLUMNS = (
    "contig", "statistic", "lag", "null_model", "value", "null_mean", "null_lo", "null_hi",
    "stat", "effect", "ci_low", "ci_high", "p_value", "q_value", "method", "n_perm", "flags",
)  # fmt: skip
#: C15 ``features`` columns, in contract order
FEATURE_COLUMNS = oo.FEATURE_COLUMNS
MIN_TUS = 3  # degenerate table: a fit with fewer than 3 TUs with values does not exist


def config(
    *,
    bootstrap: int = 199,
    permutations: int = 99,
    block_len: int | str = "auto",
    seed: int = 7,
) -> Any:
    """A ``GippaConfig`` with the §5.2 knobs the ori-ter part reads."""
    from gippa.config.model import GippaConfig

    return GippaConfig.model_validate(
        {
            "seed": seed,
            "stats": {"bootstrap": int(bootstrap), "permutations": int(permutations)},
            "context": {"islands": {"block_len": block_len}},
        }
    )


# --------------------------------------------------------------------------- input frames


def features(rows: Sequence[tuple[str, str, int, int]] = ()) -> pd.DataFrame:
    """A C15 ``features`` frame from ``(type, contig, start, end)`` rows (0-based half-open)."""
    out = pd.DataFrame(
        [
            (f"F{i + 1:04d}", t, c, int(s), int(e), ".", None, None, None, None)
            for i, (t, c, s, e) in enumerate(rows)
        ],
        columns=list(FEATURE_COLUMNS),
    )
    if not len(out):
        out = pd.DataFrame({c: pd.Series([], dtype=object) for c in FEATURE_COLUMNS})
    return out.astype({"start": "int64", "end": "int64"}) if len(out) else out


def oric(contig: str, *spans: tuple[int, int]) -> pd.DataFrame:
    """A C15 frame with one ``oric`` row per span on ``contig``."""
    return features([("oric", contig, s, e) for s, e in spans])


def oric_at(contig: str, *midpoints: float, width: int = oo.ORIC_LEN) -> pd.DataFrame:
    """A C15 frame with one ``oric`` row centred on each midpoint (width ``width``, even)."""
    half = width // 2
    return oric(contig, *((int(m) - half, int(m) + half) for m in midpoints))


def tu_frame(genes: pd.DataFrame) -> pd.DataFrame:
    """C5 ``tus`` from a gene table: 0-based half-open span of each (contig, tu_id)."""
    rows = []
    for (contig, tu), g in genes.groupby(["contig", "tu_id"], sort=False):
        rows.append(
            (
                str(tu),
                str(contig),
                int(g["start"].min()),
                int(g["end"].max()),
                str(g["strand"].iloc[0]),
                [str(t) for t in g.sort_values("start")["locus_tag"]],
            )
        )
    out = pd.DataFrame(rows, columns=["tu_id", "contig", "start", "end", "strand", "genes"])
    out = out.sort_values(["contig", "start"], kind="mergesort").reset_index(drop=True)
    return out.astype({"start": "int64", "end": "int64"})


def frames(contigs: Sequence[Contig]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """(genes, tus, replicons) for ``contigs`` (M4-17 ``gene_table`` + :func:`tu_frame`)."""
    genes, replicons = rd.gene_table(list(contigs))
    return genes, tu_frame(genes), replicons


def series_contig(
    name: str,
    values: Sequence[float],
    *,
    circular: bool = True,
    kind: str | None = None,
    detected: Sequence[int] | None = None,
    p: float | Sequence[float] = 0.5,
    tu_size: int = 1,
) -> Contig:
    """A replicon whose detected genes carry ``values`` as ``abundance_log2`` (TU runs)."""
    n = len(values)
    det = [1] * n if detected is None else [int(x) for x in detected]
    probs = [float(p)] * n if isinstance(p, (int, float)) else [float(x) for x in p]
    return Contig(
        name=name,
        d=det,
        p=probs,
        a=[float(v) for v in values],
        circular=circular,
        kind=kind,
        tu=[f"TU_{name}_{i // tu_size}" for i in range(n)],
    )


@dataclass(frozen=True)
class Placed:
    """A replicon with genes at chosen midpoints (the examples need exact positions).

    ``mids`` are gene midpoints in 0-based coordinates; gene i spans
    [mid - width/2, mid + width/2) and may pass ``length`` (unwrapped, §6.5). ``values`` are the
    log2 abundances of the detected genes (default 10 + i), ``tu`` the TU label per gene
    (default one TU per gene).
    """

    name: str
    length: int
    mids: Sequence[float]
    values: Sequence[float] | None = None
    detected: Sequence[int] | None = None
    width: int = 100
    tu: Sequence[str] | None = None
    topology: str = "circular"
    kind: str = "chromosome"


def placed(specs: Sequence[Placed]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """(genes, tus, replicons) for replicons whose genes sit at given midpoints."""
    rows: list[dict[str, Any]] = []
    reps = []
    for spec in specs:
        n = len(spec.mids)
        det = [1] * n if spec.detected is None else [int(v) for v in spec.detected]
        vals = (
            [10.0 + i for i in range(n)] if spec.values is None else [float(v) for v in spec.values]
        )
        tus = (
            [f"TU_{spec.name}_{i}" for i in range(n)]
            if spec.tu is None
            else [str(t) for t in spec.tu]
        )
        half = spec.width / 2.0
        for i, mid in enumerate(spec.mids):
            start, end = int(mid - half), int(mid + half)
            if start < 0:  # an origin-crossing gene is written unwrapped (§6.5)
                start, end = start + spec.length, end + spec.length
            rows.append(
                {
                    "locus_tag": f"{spec.name}_g{i:05d}",
                    "feature_type": "cds",
                    "contig": spec.name,
                    "start": start,
                    "end": end,
                    "strand": "+",
                    "partial": False,
                    "crosses_origin": end > spec.length,
                    "tu_id": tus[i],
                    "region_ids": [],
                    "detected": bool(det[i]),
                    "detection_status": "detected" if det[i] else "not_detected",
                    "p_detect": 0.5,
                    "abundance_log2": vals[i] if det[i] else math.nan,
                    "length_aa": 99,
                    "n_observable_peptides": 5,
                    "gravy": -0.2,
                }
            )
        reps.append((spec.name, int(spec.length), spec.topology, spec.kind))
    genes = pd.DataFrame(rows, columns=list(rd.GENE_COLUMNS)).astype(
        {
            "start": "int64",
            "end": "int64",
            "partial": bool,
            "detected": bool,
            "length_aa": "int64",
            "n_observable_peptides": "int64",
            "gravy": "float64",
        }
    )
    replicons = pd.DataFrame(reps, columns=["contig", "length", "topology", "kind"]).astype(
        {"length": "int64"}
    )
    return genes, tu_frame(genes), replicons


def right_arm_mids(oric: float, length: int, positions: Sequence[float]) -> list[float]:
    """Gene midpoints whose right-replichore positions are exactly ``positions``."""
    return [(oric + p * length / 2.0) % length for p in positions]


def replicon_length(replicons: pd.DataFrame, contig: str) -> int:
    """C1 length of ``contig``."""
    row = replicons[replicons["contig"].astype(str) == contig]
    return int(row["length"].iloc[0])


# --------------------------------------------------------------------------- recorded Bakta


def _fasta(path: Path) -> dict[str, str]:
    raw = path.read_bytes()
    text = (gzip.decompress(raw) if path.suffix == ".gz" else raw).decode()
    out: dict[str, list[str]] = {}
    cur = ""
    for line in text.splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out[cur] = []
        elif line:
            out[cur].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


@cache
def recorded(name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(C15 ``features``, C1 ``replicons``) of a recorded Bakta 1.12.2 full-DB JSON.

    ``name`` is ``k12`` | ``nissle`` | ``sim100k``. The JSON carries no nucleotides (M2 lead
    decision 4), so the committed genome FASTA supplies the replicon lengths.
    """
    stem = {"k12": "k12_mg1655", "nissle": "nissle1917", "sim100k": "sim100k"}[name]
    from gippa.annotate.bakta import parse_bakta_json

    ann = parse_bakta_json(
        RECORDED / name / f"{stem}.json.gz", genome=_fasta(GENOMES / f"{stem}.fna.gz")
    )
    feats = features(
        [(f.type, f.interval.contig, f.interval.start, f.interval.end) for f in ann.model.features]
    )
    feats["feature_id"] = [f.id for f in ann.model.features]
    reps = pd.DataFrame(
        [(c, r.length, r.topology, r.kind) for c, r in ann.model.replicons.items()],
        columns=["contig", "length", "topology", "kind"],
    ).astype({"length": "int64"})
    return feats, reps


# --------------------------------------------------------------------------- oriC oracle


def circular_mean(midpoints: Sequence[float], length: float) -> float:
    """Circular mean of ``midpoints`` on a circle of ``length`` (angles 2 pi m / L)."""
    ang = 2.0 * math.pi * np.asarray(midpoints, dtype=float) / length
    mean = math.atan2(float(np.sin(ang).mean()), float(np.cos(ang).mean()))
    return float((mean / (2.0 * math.pi) * length) % length)


def oracle_oric(feats: pd.DataFrame | None, contig: str, length: int) -> float | None:
    """The oriC midpoint used for ``contig``, or None (spec "oriC", lead decision 31).

    One ``oric`` row on the contig gives its midpoint; several give their circular mean, but
    only when every pair lies within ``ORIC_SPREAD_FRAC`` * L circular distance -- otherwise the
    analysis is skipped (``no_oric``).
    """
    if feats is None or not len(feats):
        return None
    f = feats[(feats["type"].astype(str) == "oric") & (feats["contig"].astype(str) == contig)]
    mids = [oo.midpoint(float(s), float(e), float(length)) for s, e in zip(f["start"], f["end"],
                                                                          strict=True)]  # fmt: skip
    if not mids:
        return None
    if len(mids) == 1:
        return float(mids[0])
    spread = max(oo.circ_dist(a, b, float(length)) for a in mids for b in mids)
    if spread > ORIC_SPREAD_FRAC * length:
        return None
    return circular_mean(mids, float(length))


def analysed(replicons: pd.DataFrame) -> list[str]:
    """Replicons the ori-ter analysis covers: circular chromosomes, in C1 order (Q3)."""
    sel = (replicons["kind"].astype(str) == "chromosome") & (
        replicons["topology"].astype(str) == "circular"
    )
    return [str(c) for c in replicons.loc[sel, "contig"]]


# --------------------------------------------------------------------------- TU oracle


@dataclass(frozen=True)
class OracleTus:
    """TUs of one replicon with a value: ids in genome order, x, y and the arm label."""

    contig: str
    length: int
    oric: float | None
    tu_id: tuple[str, ...]
    x: np.ndarray
    y: np.ndarray
    arm: tuple[str, ...]


def tu_values(genes: pd.DataFrame, contig: str) -> dict[str, float]:
    """y_t per TU: mean ``abundance_log2`` of the TU's detected eligible genes with a value."""
    from gippa.context.regions import replicon_index

    ix = replicon_index(genes, contig)
    out: dict[str, float] = {}
    if not len(ix):
        return out
    a = pd.to_numeric(ix["a"], errors="coerce").to_numpy(float)
    for tu, sub in pd.DataFrame({"tu": ix["tu_id"].astype(str), "a": a}).groupby("tu", sort=False):
        vals = sub["a"].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        if vals.size:
            out[str(tu)] = float(vals.mean())
    return out


def oracle_tus(
    genes: pd.DataFrame,
    tus: pd.DataFrame,
    replicons: pd.DataFrame,
    feats: pd.DataFrame | None,
    contig: str,
) -> OracleTus:
    """The TU table of one replicon: genome order, x = replichore position of the TU midpoint."""
    length = replicon_length(replicons, contig)
    o = oracle_oric(feats, contig, length)
    vals = tu_values(genes, contig)
    t = tus[tus["contig"].astype(str) == contig].sort_values(["start", "tu_id"], kind="mergesort")
    ids, xs, ys, arms = [], [], [], []
    for tu, s, e in zip(t["tu_id"].astype(str), t["start"], t["end"], strict=True):
        if tu not in vals or o is None:
            continue
        mid = oo.midpoint(float(s), float(e), float(length))
        arm, pos = oo.replichore(mid, o, float(length))
        ids.append(tu)
        xs.append(pos)
        ys.append(vals[tu])
        arms.append(arm)
    return OracleTus(
        contig=contig,
        length=length,
        oric=o,
        tu_id=tuple(ids),
        x=np.asarray(xs, dtype=float),
        y=np.asarray(ys, dtype=float),
        arm=tuple(arms),
    )


def arm_order(ot: OracleTus, arm: str) -> np.ndarray:
    """Row positions of one fit's sequence: pooled in genome order, an arm oriC -> terminus."""
    if arm == "pooled":
        return np.arange(ot.x.size, dtype=int)
    sel = np.flatnonzero(np.asarray(ot.arm, dtype=object) == arm)
    return sel[np.argsort(ot.x[sel], kind="mergesort")]


def ols_slope(x: np.ndarray, y: np.ndarray) -> float:
    """OLS slope of y on x; NaN when fewer than 3 points or x has no variance."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < MIN_TUS:
        return math.nan
    sxx = float(((x - x.mean()) ** 2).sum())
    if sxx <= 0.0:
        return math.nan
    return float(((x - x.mean()) * (y - y.mean())).sum() / sxx)


def ols_resid(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Residuals of the OLS fit (NaN-free; the intercept makes them centred)."""
    b = ols_slope(x, y)
    if math.isnan(b):
        return np.zeros(0)
    a = float(y.mean() - b * x.mean())
    return np.asarray(y - (a + b * x), dtype=float)


# --------------------------------------------------------------------------- row oracle


def _na_row(contig: str, arm: str, flags: Sequence[str]) -> dict[str, Any]:
    return {
        "contig": contig, "statistic": STATISTIC, "lag": math.nan, "null_model": NULL_MODEL,
        "value": math.nan, "null_mean": math.nan, "null_lo": math.nan, "null_hi": math.nan,
        "stat": f"{STATISTIC}_{arm}", "effect": math.nan, "ci_low": math.nan,
        "ci_high": math.nan, "p_value": math.nan, "q_value": math.nan, "method": METHOD,
        "n_perm": 0, "flags": tuple(sorted(flags)),
    }  # fmt: skip


def block_length(ell: int | str, l_hat: int, t: int) -> int:
    """ADR-0021 decision 1: an integer config wins (clamped to T), else l' = min(20, 2 l + 1, T)."""
    from gippa.stats.autocorr import block_len_prime

    if isinstance(ell, int):
        return int(min(ell, t))
    return int(block_len_prime(l_hat, t))


def oracle_replicon(
    ot: OracleTus, cfg: Any, rng: np.random.Generator
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Rows of one circular chromosome plus a diagnostics dict (``l_hat``, ``block_len``)."""
    from gippa.stats.autocorr import choose_block_len
    from gippa.stats.bootstrap import circular_block_bootstrap, percentile_ci

    diag: dict[str, Any] = {"l_hat": None, "block_len": {}}
    if ot.oric is None:
        return [_na_row(ot.contig, "pooled", (DESCRIPTIVE, NO_ORIC))], diag

    b = int(cfg.stats.bootstrap)
    conf = cfg.context.islands.block_len
    pooled = arm_order(ot, "pooled")
    pooled_slope = ols_slope(ot.x[pooled], ot.y[pooled])
    l_hat = 0
    if not math.isnan(pooled_slope) and not isinstance(conf, int):
        l_hat = int(
            choose_block_len(
                ols_resid(ot.x[pooled], ot.y[pooled]),
                True,
                int(cfg.stats.permutations),
                rng,
                alpha=0.05,
                max_len=20,
            )
        )
        diag["l_hat"] = l_hat

    rows = []
    for arm in ARMS:
        idx = arm_order(ot, arm)
        x, y = ot.x[idx], ot.y[idx]
        slope = ols_slope(x, y)
        if math.isnan(slope):
            rows.append(_na_row(ot.contig, arm, (DESCRIPTIVE, NOT_ESTIMABLE)))
            continue
        ell = block_length(conf, l_hat, int(x.size))
        diag["block_len"][arm] = ell
        draws = circular_block_bootstrap(np.arange(x.size), ell, b, rng).astype(int)
        slopes = np.array([ols_slope(x[d], y[d]) for d in draws], dtype=float)
        valid = slopes[np.isfinite(slopes)]
        flags = [DESCRIPTIVE]
        if valid.size < b / 2:
            lo = hi = math.nan
            flags.append(NOT_ESTIMABLE)
        else:
            lo, hi = percentile_ci(valid, 0.95)
        rows.append(
            {
                "contig": ot.contig,
                "statistic": STATISTIC,
                "lag": math.nan,
                "null_model": NULL_MODEL,
                "value": slope,
                "null_mean": math.nan,
                "null_lo": math.nan,
                "null_hi": math.nan,
                "stat": f"{STATISTIC}_{arm}",
                "effect": slope,
                "ci_low": lo,
                "ci_high": hi,
                "p_value": math.nan,
                "q_value": math.nan,
                "method": METHOD,
                "n_perm": int(valid.size),
                "flags": tuple(sorted(flags)),
            }
        )
    return rows, diag


def oracle_rows(
    genes: pd.DataFrame,
    tus: pd.DataFrame,
    replicons: pd.DataFrame,
    feats: pd.DataFrame | None,
    cfg: Any,
    seed: int | None = None,
) -> pd.DataFrame:
    """The expected ori-ter rows of the whole genome (replicons in C1 order)."""
    from gippa.rng import for_step

    rng = for_step(STEP, cfg.seed if seed is None else seed)
    out: list[dict[str, Any]] = []
    for contig in analysed(replicons):
        rows, _ = oracle_replicon(oracle_tus(genes, tus, replicons, feats, contig), cfg, rng)
        out.extend(rows)
    return pd.DataFrame(out, columns=list(COLUMNS))


def oracle_gene_columns(
    genes: pd.DataFrame, replicons: pd.DataFrame, feats: pd.DataFrame | None
) -> pd.DataFrame:
    """Expected ``gene_replichores``: one row per gene of the gene table, NA off the analysis."""
    oric_of = {c: oracle_oric(feats, c, replicon_length(replicons, c)) for c in analysed(replicons)}
    length_of = dict(zip(replicons["contig"].astype(str), replicons["length"], strict=True))
    arms, pos = [], []
    for contig, s, e in zip(genes["contig"].astype(str), genes["start"], genes["end"],
                            strict=True):  # fmt: skip
        o = oric_of.get(contig)
        if o is None:
            arms.append(None)
            pos.append(math.nan)
            continue
        length = float(length_of[contig])
        arm, p = oo.replichore(oo.midpoint(float(s), float(e), length), o, length)
        arms.append(arm)
        pos.append(p)
    return pd.DataFrame(
        {
            "locus_tag": genes["locus_tag"].astype(str).to_numpy(),
            "replichore": arms,
            "replichore_position": pos,
        }
    )


# --------------------------------------------------------------------------- driving the step


def rows(
    genes: pd.DataFrame,
    tus: pd.DataFrame,
    replicons: pd.DataFrame,
    feats: pd.DataFrame | None,
    cfg: Any,
    seed: int | None = None,
) -> pd.DataFrame:
    """``gippa.context.oriter.oriter_rows`` with the spec's stream."""
    from gippa.context.oriter import oriter_rows
    from gippa.rng import for_step

    rng = for_step(STEP, cfg.seed if seed is None else seed)
    return oriter_rows(genes, tus, replicons, feats, cfg, rng)  # type: ignore[no-any-return]


def gene_columns(
    genes: pd.DataFrame, replicons: pd.DataFrame, feats: pd.DataFrame | None
) -> pd.DataFrame:
    """``gippa.context.oriter.gene_replichores``."""
    from gippa.context.oriter import gene_replichores

    return gene_replichores(genes, replicons, feats)  # type: ignore[no-any-return]


def row_of(frame: pd.DataFrame, arm: str, contig: str | None = None) -> Mapping[str, Any]:
    """The single row of one fit (``stat = oriter_slope_<arm>``) of ``contig``."""
    sel = frame[frame["stat"].astype(str) == f"{STATISTIC}_{arm}"]
    if contig is not None:
        sel = sel[sel["contig"].astype(str) == contig]
    if len(sel) != 1:
        raise AssertionError(f"expected one {arm} row for {contig}, got {len(sel)}")
    return dict(sel.iloc[0])


def scenario_frames(
    rng: np.random.Generator, pair: Any
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """(genes, tus, replicons, features) drawn from a gippasim ``scenarios`` pair.

    Unlike the island helper this keeps the simulator's real coordinates (the ori-ter geometry
    is the whole point): the gene table carries ``lay.genes`` intervals, the C5 frame is
    ``lay.tus`` and the C15 frame is ``lay.features`` (with the oriC of ``oric_frac``). Evidence
    columns come from the M4-08 truth (``detected_true``, ``p_detect_true``,
    ``abundance_log2_true``), as ``region_tests_data.tables_from`` does with ``p_source =
    "true"``.
    """
    import gippasim.expression as ex
    import gippasim.layout as la

    lcfg, ecfg = pair
    lay = la.simulate_layout(lcfg, rng)
    e = ex.simulate_expression(lay, ecfg, rng)
    tu_of = {
        str(g): str(t) for t, gs in zip(lay.tus["tu_id"], lay.tus["genes"], strict=True) for g in gs
    }
    g = lay.genes
    eg = e.genes.set_index("locus_tag").loc[g["locus_tag"].astype(str)]
    partial = g["partial"].astype(bool).to_numpy()
    det = eg["detected_true"].astype(bool).to_numpy() & ~partial
    a = pd.to_numeric(eg["abundance_log2_true"], errors="coerce").to_numpy(float)
    genes = pd.DataFrame(
        {
            "locus_tag": g["locus_tag"].astype(str).to_numpy(),
            "feature_type": "cds",
            "contig": g["contig"].astype(str).to_numpy(),
            "start": g["start"].to_numpy(),
            "end": g["end"].to_numpy(),
            "strand": g["strand"].astype(str).to_numpy(),
            "partial": partial,
            "crosses_origin": g["crosses_origin"].astype(bool).to_numpy(),
            "tu_id": [tu_of[t] for t in g["locus_tag"].astype(str)],
            "region_ids": [[] for _ in range(len(g))],
            "detected": det,
            "detection_status": np.where(
                partial, "partial_excluded", np.where(det, "detected", "not_detected")
            ),
            "p_detect": np.where(
                partial, np.nan, pd.to_numeric(eg["p_detect_true"], errors="coerce").to_numpy()
            ),
            "abundance_log2": np.where(det, a, np.nan),
            "length_aa": 99,
            "n_observable_peptides": 5,
            "gravy": -0.2,
        },
        columns=list(rd.GENE_COLUMNS),
    )
    return genes, lay.tus, lay.replicons, lay.features


def planted_slope(ot: OracleTus, amplitude: float) -> float:
    """OLS slope of the planted cosine curve (a/2) cos(pi x) over the TUs of the pooled fit.

    This is the estimand of the coverage check (spec example 7): the simulator shifts every chr1
    TU by ``(a / 2) cos(pi p)`` at its replichore position p (gippasim ``_gradient``), so the
    descriptive slope the analysis targets is the OLS slope of that deterministic curve on the
    TUs that enter the fit.
    """
    x = ot.x
    return ols_slope(x, amplitude / 2.0 * np.cos(math.pi * x))


def pooled_ci(frame: pd.DataFrame, contig: str = "chr1") -> tuple[float, float, float]:
    """(value, ci_low, ci_high) of the pooled row of ``contig``."""
    row = row_of(frame, "pooled", contig)
    return float(row["value"]), float(row["ci_low"]), float(row["ci_high"])


def gradient_contig(
    rng: np.random.Generator,
    name: str,
    n_genes: int,
    *,
    amplitude: float = 0.0,
    rho: float = 0.0,
    detected_frac: float = 1.0,
    oric_frac: float = 0.0,
    noise: float = 0.5,
) -> tuple[Contig, float]:
    """A hand-built circular chromosome with a cosine ori-ter gradient of ``amplitude`` log2.

    One gene per TU, abundance = 20 + (a/2) cos(pi p) + AR(1) noise, with p the replichore
    position of the gene midpoint for an oriC at ``oric_frac`` of the replicon. Returns the
    contig and the oriC midpoint (the gene layout of ``region_tests_data.gene_table`` is
    deterministic: gene i spans [400 i + 100, 400 i + 400) and L = 400 n + 200).
    """
    length = 400 * n_genes + 200
    o = float(int(oric_frac * length))
    mids = np.array([400 * i + 250 for i in range(n_genes)], dtype=float)
    pos = np.array([oo.replichore(m, o, float(length))[1] for m in mids])
    x = np.zeros(n_genes)
    e = rng.standard_normal(n_genes)
    if n_genes:
        x[0] = e[0] / math.sqrt(1.0 - rho * rho) if rho else e[0]
        for t in range(1, n_genes):
            x[t] = rho * x[t - 1] + e[t]
    a = 20.0 + amplitude / 2.0 * np.cos(math.pi * pos) + noise * x
    d = rng.random(n_genes) < detected_frac
    return (
        Contig(
            name=name,
            d=[int(v) for v in d],
            p=[float(detected_frac)] * n_genes,
            a=[float(v) for v in a],
            circular=True,
            kind="chromosome",
            tu=[f"TU_{name}_{i}" for i in range(n_genes)],
        ),
        o,
    )
