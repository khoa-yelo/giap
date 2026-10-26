r"""Operating point of the distance-model operon fallback on RegulonDB K-12 (M11-14 14s2).

Spec ``docs/specs/M11-14-container-install.md`` 14s2; ADR-0022 addendum. Precision and recall of
``p_dist >= pair_threshold`` with p_dist = 1 / (1 + exp((d - d0) / s)), d0 and s read from the
bundled ``calibration.yaml`` and ``pair_threshold`` from ``defaults.yaml`` (0.5, Khoa 2026-10-24),
on the same labelled pair set L as the M6-12 driver (``validation/calibrate_cutoffs.py``, identity
checks included), with directon-bootstrap percentile CIs (t fixed, never re-selected; the M6-12
``fixed_cutoff_ci``). UniOP at the same threshold and the paired differences UniOP - distance
(same directon replicates) are reported next to it. No p or q value is produced (no §8.1 family).
Nothing written here is committed (R6 / R11): the output names counts and estimates only.

Usage::

    python validation/operon_fallback_operating_point.py \
        --regulondb-dir $GIAP_DATA_DIR/regulondb/14.5.0 \
        --bakta-json tests/data/recorded/bakta/1.12.2-full/k12/k12_mg1655.json.gz \
        --genome tests/data/genomes/k12_mg1655.fna.gz \
        --uniop-dir $GIAP_ROOT/recordings/uniop/k12_mg1655 \
        --terminatornet tests/data/recorded/terminatornet/721803095459/k12/terms.tsv.gz \
        --out $GIAP_ROOT/stat/m11-14s2/operating_point.json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from validation.calibrate_cutoffs import K12_PINS, _read_tsv, _verify, coordinate_id_map, pair_table

from giap.annotate.bakta import parse_bakta_json
from giap.annotate.operons import distance_model_p
from giap.annotate.uniop_adapter import parse_uniop_pred
from giap.config import GiapConfig
from giap.data.calibration import load as load_calibration
from giap.io.fasta import read_genome_fasta
from giap.manifest import finite
from giap.rng import for_step
from giap.stats import cutoffs as sc
from giap.stats.bootstrap import tu_block_bootstrap

STEP = "m11_14s2_operon_fallback"


def labelled_pairs(
    regulondb_dir: Path, bakta_json: Path, genome: Path, uniop_dir: Path, terminatornet: Path,
    contig: str = "U00096.3",
) -> pd.DataFrame:  # fmt: skip
    """L of the M6-12 driver: labelled adjacent same-strand pairs with ``y``, ``uniop_p``, ``d``
    and ``directon_id`` (identities verified against the ADR-0022 pins)."""
    _verify(regulondb_dir, genome, uniop_dir, terminatornet, K12_PINS)
    seqs = {r.id: r.sequence for r in read_genome_fasta(genome)}
    model = parse_bakta_json(bakta_json, genome=seqs).model
    rep = model.replicons[contig]
    genes = [g for g in model.genes if g.interval.contig == contig]
    gff, pred = uniop_dir / "uniop_cds.gff", uniop_dir / "uniop.pred"
    probs = parse_uniop_pred(pred, gff, coordinate_id_map(gff, model.genes))
    tab = pair_table(genes, rep, _read_tsv(regulondb_dir / "truth" / "pair_truth.tsv"), probs)
    return tab[tab["uniop_p"].notna()].reset_index(drop=True)


def _point(score: np.ndarray, y: np.ndarray, t: float) -> dict[str, Any]:
    called = score >= t
    tp, fp = int((called & (y == 1)).sum()), int((called & (y == 0)).sum())
    p = int(y.sum())
    return {"precision": tp / (tp + fp) if tp + fp else math.nan, "recall": tp / p,
            "f1": 2 * tp / (tp + fp + p), "tp": tp, "fp": fp, "n_called": tp + fp}  # fmt: skip


def operating_point(
    L: pd.DataFrame, d0: float, s: float, t: float, *, seed: int = sc.DEFAULT_SEED,
    B: int = sc.DEFAULT_B,
) -> dict[str, Any]:  # fmt: skip
    """Point estimates and directon / iid bootstrap CIs. Stream order: per cluster scheme
    (directon, then iid), distance precision and recall (``fixed_cutoff_ci``), UniOP precision and
    recall, then the paired precision and recall differences."""
    y = L["y"].to_numpy(dtype=np.int64)
    d = L["d"].to_numpy(dtype=float)
    pd_ = np.array([distance_model_p(int(v), d0, s) for v in d])
    pu = L["uniop_p"].to_numpy(dtype=float)
    rng = for_step(STEP, seed)
    out: dict[str, Any] = {
        "distance_model": _point(pd_, y, t),
        "uniop": _point(pu, y, t),
        "d_max_called_bp": int(d[pd_ >= t].max()),
    }
    for suffix, cl in (("", L["directon_id"].to_numpy()), ("_iid", np.arange(len(L)))):
        for name, score in (("distance_model", pd_), ("uniop", pu)):
            prec, rec = sc.fixed_cutoff_ci(score, y, cl, t, rng=rng, B=B)
            out[name]["precision_ci" + suffix] = [prec.ci_low, prec.ci_high]
            out[name]["recall_ci" + suffix] = [rec.ci_low, rec.ci_high]
        df = pd.DataFrame({"a": pu, "b": pd_, "y": y, "cluster": cl})
        for stat, fn in (("precision", sc.precision_at), ("recall", sc.recall_at)):
            ci = tu_block_bootstrap(
                df, "cluster", lambda x, f=fn: f(x["a"], x["y"], t) - f(x["b"], x["y"], t),
                B, rng, 0.95,
            )  # fmt: skip
            out.setdefault("uniop_minus_distance", {})[stat] = (
                out["uniop"][stat] - out["distance_model"][stat]
            )
            out["uniop_minus_distance"][f"{stat}_ci{suffix}"] = [ci.ci_low, ci.ci_high]
    out |= {"n_labelled": len(L), "n_positive": int(y.sum()), "n_negative": int((y == 0).sum()),
            "n_directons": int(L["directon_id"].nunique()), "pair_threshold": t, "d0": d0, "s": s,
            "seed": seed, "B": B, "step": STEP}  # fmt: skip
    return out


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    for flag in ("--regulondb-dir", "--bakta-json", "--genome", "--uniop-dir", "--terminatornet"):
        ap.add_argument(flag, type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=sc.DEFAULT_SEED)
    ap.add_argument("--B", type=int, default=sc.DEFAULT_B)
    a = ap.parse_args(argv)
    L = labelled_pairs(a.regulondb_dir, a.bakta_json, a.genome, a.uniop_dir, a.terminatornet)
    model = load_calibration()["operon_distance"]
    t = float(GiapConfig().annotation.uniop.pair_threshold)
    res = operating_point(L, float(model["d0"]), float(model["s"]), t, seed=a.seed, B=a.B)
    res["date"] = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(finite(res), indent=2) + "\n", encoding="utf-8")
    print(json.dumps(finite(res), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
