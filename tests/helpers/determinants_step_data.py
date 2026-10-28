"""Inputs and checks for the M8-10 tests (docs/specs/M8-10-determinants-step-and-run.md).

Contract choices the test author made where the spec leaves them open (each one is reported to
the lead):

* **Interface.** ``gippa.steps.determinants.determinants_step`` is registered as ``determinants``
  (kind ``py``, module ``core``) with the spec's input / optional-input / output keys; the tests
  call it through ``gippa.steps.registry.run_step`` (what ``gippa _step`` and the local executor
  call). ``gippa.determinants.CAVEATS`` is the three §7.8.2 caveats, verbatim, in order.
* **Writers.** ``gippa.io.determinants_tables`` has ``write_determinants_features``,
  ``write_determinants_model`` and ``write_determinants_attribution`` (``(path, frame)``, like
  ``write_sensitivity_tsv``); the schemas are ``src/gippa/schemas/tables/
  determinants_{features,model,attribution}.schema.json`` (the ``<table>.schema.json``
  convention, M4-15 lead decision 16).
* **Columns.** ``determinants_features.tsv``: ``gippa.determinants.features.FEATURE_COLUMNS``;
  ``determinants_model.tsv``: ``gippa.determinants.models.TABLE_COLUMNS`` (§9.3 plus ``flags``,
  index decision 12); ``determinants_attribution.tsv``:
  ``gippa.determinants.attribution.ATTRIBUTION_COLUMNS`` (§9.3).
* **B.** The null's B is ``determinants.permutations`` (the §5.2 key); ``--publication`` sets it
  to 10,000 (lead decision 3), as it already sets ``stats.permutations``. The univariate
  bootstrap B is ``stats.bootstrap``.
* **Streams.** The univariate rows equal ``univariate(features, cfg.stats.bootstrap,
  for_step("determinants.bootstrap", cfg.seed))`` exactly (the stream is not shared with the CV
  or the null), where ``features`` is ``feature_frame`` on the step's inputs.
* **Order.** Model rows sorted by model (lightgbm, elastic_net, univariate), then metric
  (cv_r2_mean, cv_r2_sd, null_p, delta_r2, spearman_rho, auroc), then ``feature_or_block`` as a
  plain string sort (the spec says "then feature_or_block" without naming another order).
* **Flags.** ``flags`` is ``;``-joined; a row with no flag may be written ``NA`` or empty.
  ``null_p``: ``null_uncalibrated`` always, ``unestimated`` when the observed R² is NaN, CI NA,
  ``n_perm`` = the non-NaN null values; ``cv_r2_mean``: ``ci_unestimated`` and NA CI (decision 4).
  Fewer than 100 model rows: exactly one row per model, ``cv_r2_mean`` / ``all`` / value NA,
  flags ``too_few_rows`` (and ``ci_unestimated``), no ``cv_r2_sd`` / ``null_p`` / ``delta_r2``.
* **No abundance** means every ``abundance_log2`` of the gene table is NA (the M4-14 ``none``
  table). The run flag is a log record whose message contains ``layer_missing:abundance`` (the
  ``context_stats`` ``layer_missing:features`` pattern).
* **Events.** ``determinants.caveats`` with ``data["caveats"] == list(CAVEATS)`` and
  ``determinants.features_dropped`` with ``data["features"]`` a list holding at least every
  feature whose ``determinants_features.tsv`` column is NA in every row.
"""

from __future__ import annotations

import importlib
import json
import math
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import determinants_features_data as fd
import m4_chain as mc

STEP = "determinants"
MODULE = "gippa.steps.determinants"
REQUIRED = ("gene_table", "genes", "tus", "layers", "replicons")
OPTIONAL = ("regions", "region_genes", "gene_context")
OUTPUTS = ("determinants_features", "determinants_model", "determinants_attribution")
#: §9.1 names (spec step 7: ``tables/determinants_*.tsv`` only)
PUBLISHED = {k: f"tables/{k}.tsv" for k in OUTPUTS}
SCHEMAS = {k: k for k in OUTPUTS}

MODELS = ("lightgbm", "elastic_net")
MODEL_ORDER = ("lightgbm", "elastic_net", "univariate")
METRIC_ORDER = ("cv_r2_mean", "cv_r2_sd", "null_p", "delta_r2", "spearman_rho", "auroc")
MODEL_COLUMNS = (
    "model", "feature_or_block", "metric", "value", "ci_low", "ci_high", "n", "n_perm", "flags",
)  # fmt: skip
BLOCKS = ("transcription", "translation", "codon", "protein", "context")
ATTRIBUTION_COLUMNS = (
    "locus_tag", "abundance_log2", "predicted_log2", "residual",
    *(f"shap_{b}" for b in BLOCKS), "dominant_block", "top_decile", "flags",
)  # fmt: skip
MODEL_FLAGS = frozenset(
    {"too_few_rows", "unestimated", "null_degenerate", "null_uncalibrated", "ci_unestimated"}
)
ATTRIBUTION_FLAGS = frozenset(
    {"too_few_rows", "unestimated", "abundance_not_explained_by_sequence", "cai_reference"}
)
#: §7.8.2 "Caveats (printed in the report)", verbatim
CAVEATS = (
    "In a synthetic GFP library, mRNA folding near the RBS explained over half the variance; "
    "native genes add regulation and turnover, so expect a lower R².",
    "All features are predictions, and the data come from one condition.",
    "CAI, tAI and GC3 are collinear, so interpret blocks rather than single features.",
)
CAVEATS_EVENT = "determinants.caveats"
DROPPED_EVENT = "determinants.features_dropped"
RUN_FLAG = "layer_missing:abundance"


def mod() -> ModuleType:
    """``gippa.steps.determinants`` (first, so a missing module fails before any fixture work)."""
    return importlib.import_module(MODULE)


# ---------------------------------------------------------------------------- fixtures


def build_work(root: Path, seed: int, size: str = "sim1000") -> Path:
    """The M4-07/M4-08 fixture (``build_m4_fixture``) through the M4 chain: the contract tables,
    ``genetable/genes.parquet`` and M5-09's ``context/gene_context.parquet``."""
    work = mc.build(root / "work", seed=seed, size=size)
    db = mc.test_db(root / "db")
    mc.run_m4(work, mc.config(seed, stats={"permutations": 99, "bootstrap": 100}), db)
    return work


def inputs(work: Path, *, gene_table: Path | None = None, optional: bool = True) -> dict[str, Path]:
    """Step inputs of a ``build_work`` directory (``gene_table`` replaces the M4 gene table)."""
    out = {
        "gene_table": gene_table or work / "genetable" / "genes.parquet",
        "genes": work / "genes.parquet",
        "tus": work / "tus.parquet",
        "layers": work / "layers.parquet",
        "replicons": work / "replicons.parquet",
    }
    if optional:
        out |= {
            "regions": work / "regions.parquet",
            "region_genes": work / "region_genes.parquet",
            "gene_context": work / "context" / "gene_context.parquet",
        }
    return out


def outputs(out: Path) -> dict[str, Path]:
    return {k: out / PUBLISHED[k] for k in OUTPUTS}


def config(seed: int, *, B: int = 49, threads: int = 1, boot: int = 49, **det: Any) -> Any:
    """A small-but-complete determinants configuration (unit budget): B permutations, two CV
    repeats of three folds, 30 boosting rounds; ``stats.permutations`` is left at its default so
    a step reading the wrong key shows up as ``n_perm != B``."""
    d: dict[str, Any] = {
        "permutations": B,
        "repeats": 2,
        "outer_folds": 3,
        "lightgbm": {"n_estimators": 30},
        "elastic_net": {"inner_folds": 2},
    }
    d.update(det)
    return mc.config(seed, threads=threads, stats={"bootstrap": boot}, determinants=d)


def run(cfg: Any, ins: Mapping[str, Path], out: Path) -> dict[str, Path]:
    """``run_step("determinants", ...)``; returns the output paths."""
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    outs = outputs(out)
    reg.run_step(STEP, cfg, dict(ins), outs)
    return outs


def gene_table_variant(work: Path, dest: Path, *, keep_abundance: int | None) -> Path:
    """A copy of the M4 gene table with ``abundance_log2`` NA everywhere (``keep_abundance=None``:
    the no-abundance run of spec step 5) or kept only on the first ``keep_abundance`` detected,
    non-imputed genes (a below-the-minimum run, like sim100k)."""
    g = pd.read_parquet(work / "genetable" / "genes.parquet")
    det = g["detected"].astype("boolean").fillna(False).to_numpy(bool)
    imp = g["abundance_imputed"].astype("boolean").fillna(False).to_numpy(bool)
    keep = pd.Series(False, index=g.index)
    if keep_abundance is not None:
        idx = g.index[det & ~imp & g["abundance_log2"].notna().to_numpy()][:keep_abundance]
        keep.loc[idx] = True
    g.loc[~keep, "abundance_log2"] = math.nan
    if "abundance_raw" in g.columns:
        g.loc[~keep, "abundance_raw"] = math.nan
    g.loc[~keep, "abundance_imputed"] = False
    if keep_abundance is None and "abundance_type" in g.columns:
        g["abundance_type"] = "none"
    dest.parent.mkdir(parents=True, exist_ok=True)
    g.to_parquet(dest, index=False)
    return dest


def expected_features(ins: Mapping[str, Path]) -> pd.DataFrame:
    """``feature_frame`` on the step's input files (the M8-07 function the step must call)."""
    from gippa.determinants.features import feature_frame

    regions = None
    if "regions" in ins and "region_genes" in ins:
        regions = fd.regions_arg(pd.read_parquet(ins["regions"]),
                                 pd.read_parquet(ins["region_genes"]))  # fmt: skip
    ctx = pd.read_parquet(ins["gene_context"]) if "gene_context" in ins else None
    return feature_frame(
        pd.read_parquet(ins["gene_table"]),
        fd.genes_from_frame(pd.read_parquet(ins["genes"])),
        fd.tus_from_frame(pd.read_parquet(ins["tus"])),
        pd.read_parquet(ins["layers"]),
        regions,
        ctx,
        fd.replicons_from_frame(pd.read_parquet(ins["replicons"])),
    )


# ---------------------------------------------------------------------------- reading


def read(path: Path) -> pd.DataFrame:
    """A written TSV as strings (``NA`` verbatim)."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def header(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").splitlines()[0].split("\t")


def flags(cell: str) -> set[str]:
    return set() if cell in ("", "NA") else set(cell.split(";"))


def num(cell: str) -> float:
    return math.nan if cell in ("", "NA") else float(cell)


def rows(tbl: pd.DataFrame, **eq: str) -> pd.DataFrame:
    m = pd.Series(True, index=tbl.index)
    for k, v in eq.items():
        m &= tbl[k] == v
    return tbl[m]


def sort_key(tbl: pd.DataFrame) -> list[tuple[int, int, str]]:
    return [
        (MODEL_ORDER.index(m), METRIC_ORDER.index(k), f)
        for m, k, f in zip(tbl["model"], tbl["metric"], tbl["feature_or_block"], strict=True)
    ]


def fmt(v: Any) -> str:
    """One cell as the M1-04 writer formats it."""
    from gippa.io.tsv import format_value

    if isinstance(v, float) and math.isnan(v):
        return "NA"
    if v is None or (not isinstance(v, str) and pd.isna(v)):
        return "NA"
    return format_value(v)


def schema_dir() -> Path:
    """``schemas/tables`` of the imported ``gippa`` package."""
    import gippa

    return Path(gippa.__file__).resolve().parent / "schemas" / "tables"


def descriptor(stem: str) -> dict[str, Any]:
    """The Frictionless resource ``<stem>.schema.json`` (``coordinates`` removed)."""
    path = schema_dir() / f"{stem}.schema.json"
    assert path.is_file(), f"no {path.name} in {path.parent}"
    d: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return {k: v for k, v in d.items() if k != "coordinates"}


def validate(path: Path, stem: str) -> None:
    """§9.5: the written table validates with ``frictionless`` against its resource."""
    import frictionless

    schema = frictionless.Schema.from_descriptor(descriptor(stem))
    resource = frictionless.Resource(path=path.name, basepath=str(path.parent), schema=schema)
    report = resource.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def field_names(stem: str) -> list[str]:
    return [f["name"] for f in descriptor(stem)["fields"]]


def enum_of(stem: str, field: str) -> Sequence[str] | None:
    for f in descriptor(stem)["fields"]:
        if f["name"] == field:
            enum: Sequence[str] | None = f.get("constraints", {}).get("enum")
            return enum
    raise AssertionError(f"{stem}: no field {field}")
