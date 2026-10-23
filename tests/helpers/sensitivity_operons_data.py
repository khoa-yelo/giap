"""Inputs and oracles for the M6-10 ``operon_flags`` sensitivity family
(docs/specs/M6-10-sensitivity-operon-flags.md, lead decisions 2026-10-22; GIAP_SPEC.md §8.6, G9;
M4-19 lead decisions 2 and 5; test author).

Interface assumed (spec Behaviour 1-3 and lead decision 1, test-author reading):

* ``giap.context.sensitivity.operon_flags_family(genes, cfg, rng, *, pairs_raw, tus, replicons,
  gene_models, features, promoter_windows=None, terminators=None, tx_rates=None)
  -> frozenset[str]`` registered as ``FAMILIES["operon_flags"]``. ``genes`` is the M4-15 gene
  table (the M4-19 positional argument), so the C2 gene model, which the step's contract calls
  ``genes``, is bound under ``gene_models`` (the name clash is the spec's; ``FAMILY_KW`` maps
  step key -> family keyword in one place). An evidence layer left at ``None`` is a missing layer
  (``layer_missing:<layer>``), as an absent optional input of the ``operons`` step.
* The flagged set is the ``pair_id`` column of the ``operon_flags`` frame that the merged M6-08
  chain (``giap.steps.operons.run_operons``) gives on the same tables and config.
* ``Variant`` gains ``config_patch`` (dotted keys, falsy by default); ``variants(cfg, tier)``
  under ``full`` appends ``pair_threshold=<t>`` for ``t`` in (0.3, 0.7) not equal to the primary
  ``annotation.uniop.pair_threshold``, each with ``config_patch =
  {"annotation.uniop.pair_threshold": t}`` and no ``overrides``.

Tables: ``build_chain`` (``helpers.layers_chain``, the merged producers on giapsim SIM100K) with
``pairs_raw`` probabilities respread over (0, 1) by ``spread_pairs_raw`` (the recorded UniOP run
puts every pair near 0.73, so no threshold in [0.3, 0.7] would move a class). ``gene_table``
builds a §9.3 gene table over the C2 genes whose abundances are constant on each same-strand run
(OP / NO pairs co-abundant, CV / DV / SW decoys not), with planted jumps inside runs.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from giap.config import GiapConfig

FAMILY = "operon_flags"
#: M4-19 decision 2 (common random numbers): one stream per family
STREAM = "sensitivity/operon_flags"
SENS_B = 199
#: spec Behaviour 3 (``full``)
UNIOP_THRESHOLDS = (0.3, 0.7)
PATCH_KEY = "annotation.uniop.pair_threshold"
#: step / chain key -> family keyword (C2 ``genes`` clashes with the positional gene table)
FAMILY_KW: Mapping[str, str] = {
    "pairs_raw": "pairs_raw", "tus": "tus", "replicons": "replicons", "genes": "gene_models",
    "features": "features", "promoter_windows": "promoter_windows",
    "terminators": "terminators", "tx_rates": "tx_rates",
}  # fmt: skip
LAYER_KEYS = ("promoter_windows", "terminators", "tx_rates")
LAYER_NAME = {"promoter_windows": "promoters", "terminators": "terminators", "tx_rates": "tx_rates"}
#: every gene whose index is a multiple of this is undetected
UNDETECTED_EVERY = 9
#: a gene whose position in its strand run is a multiple of this jumps by ``JUMP`` log2
JUMP_EVERY, JUMP = 4, 3.0


def variant_name(t: float) -> str:
    return f"pair_threshold={t:g}"


def config(**sections: Mapping[str, Any]) -> GiapConfig:
    """``GiapConfig(seed=1)`` with ``abundance.missing = min`` (sensitivity ``[min, exclude]``
    then reruns only ``missing=exclude``) and the given sections merged in."""
    raw: dict[str, Any] = {"seed": 1, "abundance": {"missing": "min"}}
    for k, v in sections.items():
        raw[k] = {**raw.get(k, {}), **dict(v)}
    return GiapConfig.model_validate(raw)


def with_threshold(cfg: GiapConfig, t: float) -> GiapConfig:
    ann = cfg.annotation
    uniop = ann.uniop.model_copy(update={"pair_threshold": t})
    return cfg.model_copy(update={"annotation": ann.model_copy(update={"uniop": uniop})})


# --------------------------------------------------------------------------- tables


def spread_pairs_raw(pairs_raw: pd.DataFrame) -> pd.DataFrame:
    """``pairs_raw`` with ``uniop_p`` respread deterministically over (0, 1)."""
    out = pairs_raw.copy()
    n = len(out)
    out["uniop_p"] = ((np.arange(n) * 37) % n + 0.5) / n
    return out


def _strand_runs(c2: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Per C2 row: index of its same-strand run (genome order per contig) and its position in it."""
    order = c2.sort_values(["contig", "start"], kind="stable").index
    run = np.zeros(len(c2), dtype=np.int64)
    pos = np.zeros(len(c2), dtype=np.int64)
    r, p, prev = -1, 0, None
    for i in order:
        key = (c2.at[i, "contig"], c2.at[i, "strand"])
        if key != prev:
            r, p, prev = r + 1, 0, key
        run[c2.index.get_loc(i)] = r
        pos[c2.index.get_loc(i)] = p
        p += 1
    return run, pos


def gene_table(
    c2: pd.DataFrame,
    *,
    abundance: bool = True,
    imputed: Sequence[int] = (),
    single_peptide: Sequence[int] = (),
    tier: str = "mayu",
    scale: str = "absolute",
) -> pd.DataFrame:
    """A §9.3 gene table (every ``GENES_COLUMNS`` column) over the C2 genes ``c2``.

    Genes ``i % UNDETECTED_EVERY == 0`` are not detected. ``imputed`` (detected) genes have no
    observed abundance and carry the primary ``min`` fill (``abundance_imputed``, flag
    ``imputed_abundance``); ``single_peptide`` genes are ``single_peptide_unestimated`` (not
    detected) with an observed abundance, which the ``single_peptide_detected`` variant uses.
    """
    from giap.genetable.assemble import GENES_COLUMNS, GENES_TYPES

    n = len(c2)
    default = {"list": lambda: [], "bool": lambda: False, "int": lambda: None,
               "float": lambda: math.nan, "str": lambda: None}  # fmt: skip
    df = pd.DataFrame({c: [default[GENES_TYPES[c]]() for _ in range(n)] for c in GENES_COLUMNS})
    for c in ("locus_tag", "gene", "product", "feature_type", "contig", "start", "end", "strand",
              "partial", "crosses_origin", "hypothetical"):  # fmt: skip
        df[c] = c2[c].to_numpy()
    idx = np.arange(n)
    det = idx % UNDETECTED_EVERY != 0
    sp = np.zeros(n, dtype=bool)
    sp[list(single_peptide)] = True
    det &= ~sp
    imp = np.zeros(n, dtype=bool)
    imp[list(imputed)] = True
    assert not (imp & ~det).any(), "imputed genes must be detected"
    status = np.where(det, "detected", "not_detected").astype(object)
    status[sp] = "single_peptide_unestimated"
    df["detected"] = det
    df["detection_status"] = status
    df["detection_tier"] = tier
    df["p_detect"] = np.linspace(0.2, 0.9, n)
    df["tu_id"] = [f"TU_{i}" for i in idx]
    df["abundance_imputed"] = False
    df["flags"] = [[] for _ in range(n)]
    if abundance:
        run, pos = _strand_runs(c2.reset_index(drop=True))
        level = 10.0 + 2.0 * (run % 5) + 0.01 * np.sin(idx) + JUMP * (pos % JUMP_EVERY == 3)
        observed = (det & ~imp) | sp
        log2 = np.where(observed, level, np.nan)
        fill = float(np.min(log2[det & ~imp])) if (det & ~imp).any() else math.nan
        log2 = np.where(imp, fill, log2)
        df["abundance_log2"] = log2
        df["abundance_raw"] = np.exp2(log2)
        df["top3_log2"] = np.where(observed, level + 0.25 * np.cos(idx), np.nan)
        df["abundance_type"] = "ibaq"
        df["abundance_scale"] = scale
        df["imputation_policy"] = "min"
        df["abundance_imputed"] = imp
        df["flags"] = [["imputed_abundance"] if m else [] for m in imp]
    return df


def roundtrip(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    """``df`` written with the M4-15 writer and read back (the dtypes the step sees)."""
    from giap.io.genes_tsv import read_genes_parquet, write_genes_parquet

    path.parent.mkdir(parents=True, exist_ok=True)
    write_genes_parquet(path, df)
    return read_genes_parquet(path)


def tables(chain: Any, *, pairs_raw: pd.DataFrame | None = None,
           drop: Sequence[str] = ()) -> dict[str, pd.DataFrame]:  # fmt: skip
    """Step key -> frame for every table the family takes (``drop`` leaves layers out)."""
    out = {k: pd.read_parquet(chain.inputs[k]) for k in FAMILY_KW if k not in drop}
    out["pairs_raw"] = spread_pairs_raw(out["pairs_raw"]) if pairs_raw is None else pairs_raw
    return out


def call_family(genes: pd.DataFrame, cfg: GiapConfig, rng: np.random.Generator,
                t: Mapping[str, pd.DataFrame]) -> frozenset[str]:  # fmt: skip
    from giap.context import sensitivity as S

    return S.operon_flags_family(genes, cfg, rng, **{FAMILY_KW[k]: v for k, v in t.items()})  # type: ignore[no-any-return]


# --------------------------------------------------------------------------- oracle


def oracle_frames(
    genes: pd.DataFrame, cfg: GiapConfig, t: Mapping[str, pd.DataFrame]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(merged pair frame, ``operon_flags`` frame) of the merged M6-08 chain on ``t``."""
    from giap.steps.operons import run_operons

    layers = {k: t[k] for k in LAYER_KEYS if k in t}
    merged, _, flagged, _ = run_operons(genes, t["genes"], t["replicons"], t["tus"],
                                        t["pairs_raw"], t["features"], layers, cfg)  # fmt: skip
    return merged, flagged


def oracle_flagged(
    genes: pd.DataFrame, cfg: GiapConfig, t: Mapping[str, pd.DataFrame]
) -> frozenset[str]:
    return frozenset(str(p) for p in oracle_frames(genes, cfg, t)[1]["pair_id"])


# --------------------------------------------------------------------------- step


def write_inputs(chain: Any, root: Path, gt: pd.DataFrame, *,
                 drop: Sequence[str] = ()) -> dict[str, Path]:  # fmt: skip
    """Explicit ``sensitivity`` step inputs: the gene table, C1 / C5 / C7 / C8 and the pair
    tables (the spread ``pairs_raw`` written next to the gene table)."""
    from giap.io.genes_tsv import write_genes_parquet

    root.mkdir(parents=True, exist_ok=True)
    gt_path = root / "genes_table.parquet"
    write_genes_parquet(gt_path, gt)
    pr = root / "pairs_raw.parquet"
    spread_pairs_raw(pd.read_parquet(chain.inputs["pairs_raw"])).to_parquet(pr, index=False)
    keys = ("replicons", "regions", "region_genes", "tus", "genes", "features", *LAYER_KEYS)
    ins = {k: Path(chain.inputs[k]) for k in keys if k not in drop}
    ins["pairs_raw"] = pr
    ins["gene_table"] = gt_path
    return ins
