#!/usr/bin/env python3
"""Write the frozen synthetic inputs of the R golden references (M8-12 item 1; GIAP_SPEC.md §11.5).

    python scripts/golden/make_golden_inputs.py [--out DIR] [--seed INT]

``--out`` defaults to ``tests/data/golden/r/inputs``; ``--seed`` to ``GOLDEN_SEED`` from
``tests/calibration/tolerances.yaml``. Four gzip TSVs (mtime 0, no file name, so the bytes depend
only on the seed): ``gsea_stats``, ``gsea_sets``, ``ora_genes``, ``ora_categories``. numpy, the
standard library and PyYAML only, never ``giap``. The IDs and memberships are random: no KEGG,
GO or COG content. Numbers are written as fixed 6-decimal strings so R and Python parse them
to the same doubles.
"""

from __future__ import annotations

import argparse
import gzip
import io
import sys
from pathlib import Path

import numpy as np
import yaml

REPO = Path(__file__).resolve().parents[2]
TOLERANCES = REPO / "tests" / "calibration" / "tolerances.yaml"
DEFAULT_OUT = REPO / "tests" / "data" / "golden" / "r" / "inputs"
SCALE = 10**6  # numbers are stored as integer micro-units, written with 6 decimals

N_TU, N_SETS, N_PLANTED, N_TINY = 2000, 120, 20, 5
TIE_COUNT = 50  # two values x 50 TUs = 5 % exact ties
N_ZEROS = 6
N_GENES, N_CATS = 3000, 150
FRAC_UNCATEGORISED = 0.30

Table = list[tuple[str, ...]]


def golden_seed() -> int:
    data = yaml.safe_load(TOLERANCES.read_text())
    return int(data["seeds"]["GOLDEN_SEED"])


def micro_str(m: int) -> str:
    sign = "-" if m < 0 else ""
    a = abs(int(m))
    return f"{sign}{a // SCALE}.{a % SCALE:06d}"


def _make_unique(micro: np.ndarray, frozen: np.ndarray, taken: set[int]) -> np.ndarray:
    """Nudge non-frozen entries by 1e-6 until they are distinct, non-zero and not in ``taken``."""
    out = micro.copy()
    seen = set(taken)
    for i in np.argsort(out, kind="stable"):
        if frozen[i]:
            continue
        m = int(out[i])
        while m == 0 or m in seen:
            m += 1
        seen.add(m)
        out[i] = m
    return out


def gsea_inputs(rng: np.random.Generator) -> tuple[Table, Table]:
    score = rng.standard_normal(N_TU)
    score -= np.median(score)
    micro = np.rint(score * SCALE).astype(np.int64)
    frozen = np.zeros(N_TU, dtype=bool)
    special = rng.permutation(N_TU)[: 2 * TIE_COUNT + N_ZEROS]
    tie_values = (int(0.612345 * SCALE), int(-0.734567 * SCALE))
    micro[special[:TIE_COUNT]] = tie_values[0]
    micro[special[TIE_COUNT : 2 * TIE_COUNT]] = tie_values[1]
    micro[special[2 * TIE_COUNT :]] = 0
    frozen[special] = True
    micro = _make_unique(micro, frozen, {0, *tie_values})
    values = micro / SCALE
    if abs(float(np.median(values))) > 0.05:
        raise SystemExit("scores are not centred")
    ids = [f"TU_{i:04d}" for i in range(1, N_TU + 1)]
    stats: Table = [(tu, micro_str(m)) for tu, m in zip(ids, micro, strict=True)]

    sets: list[list[str]] = []
    # 20 planted sets: z = mean / sd * sqrt(n) from ~2.5 (p ~1e-2) to ~14 (p far below 1e-30),
    # alternating top and bottom, by exponential tilting of the sampling weights.
    sizes = np.rint(np.geomspace(20, 200, N_PLANTED)).astype(int)
    z_target = np.geomspace(2.5, 14.0, N_PLANTED)
    sd = float(values.std())
    for j in range(N_PLANTED):
        beta = z_target[j] / np.sqrt(sizes[j]) / sd
        sign = 1.0 if j % 2 == 0 else -1.0
        w = np.exp(sign * beta * values)
        for _ in range(1000):  # redraw until the realised z is within 10 % of the target
            idx = rng.choice(N_TU, size=sizes[j], replace=False, p=w / w.sum())
            z = abs(values[idx].mean()) / sd * np.sqrt(sizes[j])
            if abs(z - z_target[j]) <= 0.1 * z_target[j]:
                break
        else:
            raise SystemExit(f"planted set {j + 1}: z {z:.2f} misses {z_target[j]:.2f}")
        sets.append([ids[i] for i in sorted(idx)])
    # 5 tiny-but-legal sets (size 10).
    for _ in range(N_TINY):
        sets.append([ids[i] for i in sorted(rng.choice(N_TU, size=10, replace=False))])
    # 2 sets with names absent from the scores: one stays legal (25 of 30 match), one falls
    # below minSize after matching (8 of 12 match) and must be dropped by both sides.
    for n, n_absent, tag in ((30, 5, "A"), (12, 4, "B")):
        members = [ids[i] for i in sorted(rng.choice(N_TU, size=n - n_absent, replace=False))]
        members += [f"TU_ABSENT_{tag}{k}" for k in range(1, n_absent + 1)]
        sets.append(members)
    # The rest: null sets with log-uniform sizes 10-500 (one of exactly 500).
    n_null = N_SETS - len(sets)
    null_sizes = np.rint(np.exp(rng.uniform(np.log(10), np.log(500), n_null))).astype(int)
    null_sizes[0] = 500
    for n in null_sizes:
        sets.append([ids[i] for i in sorted(rng.choice(N_TU, size=int(n), replace=False))])
    long: Table = [(f"SET_{s:03d}", m) for s, members in enumerate(sets, 1) for m in members]
    return stats, long


def ora_inputs(rng: np.random.Generator) -> tuple[Table, Table]:
    p = rng.beta(2.0, 2.0, N_GENES)
    p_micro = np.clip(np.rint(p * SCALE).astype(np.int64), 1, SCALE - 1)
    edge = rng.permutation(N_GENES)[:6]
    p_micro[edge[:3]] = 0
    p_micro[edge[3:]] = SCALE
    pd = p_micro / SCALE
    detected = (rng.random(N_GENES) < pd).astype(int)
    genes = [f"GENE_{i:04d}" for i in range(1, N_GENES + 1)]
    genes_table: Table = [
        (g, micro_str(m), str(d)) for g, m, d in zip(genes, p_micro, detected, strict=True)
    ]

    order = rng.permutation(N_GENES)
    n_out = round(FRAC_UNCATEGORISED * N_GENES)
    pool = np.sort(order[n_out:])  # genes that regular categories draw from
    det_pool = pool[detected[pool] == 1]
    mid_det = det_pool[(pd[det_pool] > 0.3) & (pd[det_pool] < 0.7)]
    mid = pool[(pd[pool] > 0.3) & (pd[pool] < 0.7)]

    def pick(source: np.ndarray, n: int) -> list[int]:
        return sorted(int(i) for i in rng.choice(source, size=n, replace=False))

    def mixed(n: int, frac_det: float) -> list[int]:
        k = round(n * frac_det)
        mid_undet = np.setdiff1d(mid, mid_det)
        return sorted(pick(mid_det, k) + pick(mid_undet, n - k))

    cats: list[list[int]] = []
    # All categorised genes minus one (a mid-p_detect gene, so neither weight is 0 nor infinite);
    # the 30 % uncategorised genes stay outside every category (PR #192 B1, lead decision).
    drop = int(rng.choice(mid))
    cats.append([])  # filled below, once every other category is drawn
    # omega >> 1 and omega << 1.
    for n, lo in ((40, 0.85), (120, 0.80), (300, 0.70)):
        cats.append(pick(pool[pd[pool] > lo], n))
    for n, hi in ((40, 0.15), (120, 0.20), (300, 0.30)):
        cats.append(pick(pool[(pd[pool] < hi) & (pd[pool] > 0)], n))
    # k = n and k = 0 (small, and k = n strong: 150 all detected, p ~1e-45).
    cats.append(pick(mid_det, 8))
    cats.append(pick(np.setdiff1d(mid, mid_det), 8))
    cats.append(pick(mid_det, 150))
    # Strong planted enrichments and one depletion, at omega ~ 1.
    cats.append(mixed(200, 0.85))
    cats.append(mixed(100, 0.80))
    cats.append(mixed(200, 0.20))
    # The p_detect = 0 and = 1 genes together with 20 others.
    cats.append(sorted({int(i) for i in edge} | set(pick(mid, 20))))
    # The rest: random categories with log-uniform sizes 5-500.
    n_rand = N_CATS - len(cats)
    sizes = np.rint(np.exp(rng.uniform(np.log(5), np.log(500), n_rand))).astype(int)
    for n in sizes:
        cats.append(pick(pool, int(n)))
    categorised = set(pool.tolist()).union(*map(set, cats[1:]))
    cats[0] = sorted(int(i) for i in categorised - {drop})
    long: Table = [(f"CAT_{c:03d}", genes[g]) for c, members in enumerate(cats, 1) for g in members]
    return genes_table, long


def write_tsv(path: Path, header: tuple[str, ...], rows: Table) -> None:
    text = "".join("\t".join(r) + "\n" for r in [header, *rows])
    buf = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buf, compresslevel=9, mtime=0) as gz:
        gz.write(text.encode("utf-8"))
    path.write_bytes(buf.getvalue())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--seed", type=int, default=None)
    args = ap.parse_args(argv)
    seed = golden_seed() if args.seed is None else args.seed
    rng_gsea, rng_ora = (np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(2))
    stats, sets = gsea_inputs(rng_gsea)
    genes, cats = ora_inputs(rng_ora)
    args.out.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out / "gsea_stats.tsv.gz", ("tu_id", "score"), stats)
    write_tsv(args.out / "gsea_sets.tsv.gz", ("set_id", "member"), sets)
    write_tsv(args.out / "ora_genes.tsv.gz", ("gene", "p_detect", "detected"), genes)
    write_tsv(args.out / "ora_categories.tsv.gz", ("category", "gene"), cats)
    return 0


if __name__ == "__main__":
    sys.exit(main())
