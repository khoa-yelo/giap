"""M6-05H at step level: the ``operons`` step fits M6-02's detection-limit model when peptide
evidence is given, and ``tables/pairs.tsv`` then carries ``one_detected_pmiss`` rows
(docs/specs/M6-05-discordance.md item 4; M6-02 acceptance "Integration: with M6-05H,
``tables/pairs.tsv`` ... contains ``one_detected_pmiss`` rows"; M6-08 degenerate row "no peptides
and one run"; GIPPA_SPEC.md §7.4.3, §9.3 ``pairs.tsv`` ``test`` enum).

Interface pinned here (test-author reading):

* ``peptide_hits`` (an existing optional input of the step) is the C11 TSV
  ``mapping/peptide_hits.tsv`` (1-based, the M4-06 contract columns), read the way
  ``gippa.io.work_tables`` reads it. M6-02 needs only its ``gene`` and ``peptide`` columns.
* The step fits ``gippa.operons.detection_limit.fit_detection_limit`` on the gene table (with
  C2 ``aa``), C11 and C16 (absent here), with ``gippa.rng.for_step("operons.detection_limit",
  cfg.seed)`` (M6-02 item 7). Each one-detected OP pair is scored with the out-of-fold model of
  the **detected partner v's** TU (``model_for_tu(model, tu(v))``; M6-02 lead decision 9, review
  B1) and the undetected gene u's ``observable_peptides``.
* The one-detected event is M6-02 lead decision 9's "< k_min observable peptides seen" (review
  S1): u has fewer than k_min of its observable peptides assigned in C11, v at least k_min, and
  the pair has no both-detected t test. ``_seen_counts`` below is the test's own oracle of that
  count.
* Without peptide evidence, or with a model that is not usable, one-detected rows stay ``none``
  with ``no_detection_limit_model``.
* A model in force that needs replicate covariates the step does not join (the replicate model
  with ``detectability.covariates`` coefficients, possible once ``peptide_runs`` is given) must
  not turn every P(miss) into NaN silently (review S2): its one-detected rows are ``none`` with
  the pair flag ``pmiss_covariates_missing`` and the step logs a WARNING naming that flag.
* P(miss) has no q and is outside the BH family (M6-index lead decision 1); see
  ``PMISS_IN_BH_FAMILY`` in ``tests/unit/operons/test_discordance_one_detected.py``.

Inputs: ``helpers.layers_chain.build_chain`` (gippasim SIM100K). Three genes in four are detected,
with ``abundance_log2 = 12 + sin(i)``. ``tu_id`` comes from C5. The peptide hits are drawn
(fixed numpy seed, test data only) from a logistic in abundance and length, so the peptide
model is usable (about 1,000 training rows).
"""

from __future__ import annotations

import json
import logging
import math
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.layers_chain import Chain, build_chain
from scipy.special import expit

from gippa.config import GippaConfig
from gippa.genetable.contract import TABLES
from gippa.io.tsv import write_tsv
from gippa.stats.multitest import bh

SCHEMA_DIR = Path(__file__).resolve().parents[3] / "src" / "gippa" / "schemas" / "tables"
OUTPUTS = {"pairs": "pairs.tsv", "operons": "operons.tsv", "operon_flags": "operon_flags.tsv"}
PMISS, TTEST = "one_detected_pmiss", "both_detected_t"
PMISS_IN_BH_FAMILY = False  # M6-index lead decision 1 (see the unit test module)
DL_STEP = "operons.detection_limit"


def _cfg() -> GippaConfig:
    return GippaConfig(seed=1)


# ----------------------------------------------------------------------------- fixtures


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(Path(tmp_path_factory.mktemp("operons_pmiss_chain")))


def _tu_of(chain: Chain) -> dict[str, str]:
    return {str(t): str(r.tu_id) for r in chain.frame("tus").itertuples() for t in r.genes}


def write_gene_table(chain: Chain, path: Path, *, split_tus: bool = False) -> Path:
    """A §9.3 gene table over the C2 genes with abundance and C5 ``tu_id``. ``split_tus`` moves
    every undetected gene (index % 4 == 3) into a TU of its own, so its OP pairs cross a TU
    boundary (UniOp pairs and C5 TUs are independent sources) and u's TU trained nothing."""
    from gippa.genetable.assemble import GENES_COLUMNS, GENES_TYPES
    from gippa.io.genes_tsv import write_genes_parquet

    c2 = chain.frame("genes")
    n = len(c2)
    default = {"list": lambda: [], "bool": lambda: False, "int": lambda: None,
               "float": lambda: math.nan, "str": lambda: None}  # fmt: skip
    df = pd.DataFrame({c: [default[GENES_TYPES[c]]() for _ in range(n)] for c in GENES_COLUMNS})
    for c in ("locus_tag", "gene", "product", "feature_type", "contig", "start", "end", "strand",
              "partial", "crosses_origin", "hypothetical"):  # fmt: skip
        df[c] = c2[c].to_numpy()
    idx = np.arange(n)
    det = idx % 4 != 3
    df["detected"] = det
    df["detection_status"] = np.where(det, "detected", "not_detected")
    df["p_detect"] = np.linspace(0.2, 0.9, n)
    df["abundance_log2"] = np.where(det, 12.0 + np.sin(idx), np.nan)
    df["abundance_type"] = "ibaq"
    df["abundance_scale"] = "absolute"
    df["abundance_imputed"] = False
    tu = _tu_of(chain)
    df["tu_id"] = [tu.get(str(t)) for t in df["locus_tag"]]
    if split_tus:
        solo = zip(df.loc[~det, "tu_id"], idx[~det], strict=True)
        df.loc[~det, "tu_id"] = [f"{t}_solo_{i}" for t, i in solo]
    path.parent.mkdir(parents=True, exist_ok=True)
    write_genes_parquet(path, df)
    return path


def write_peptide_hits(
    chain: Chain, path: Path, *, max_rows: int | None = None, extra: frozenset[str] = frozenset()
) -> Path:
    """C11 rows: each observable peptide of a detected gene is seen with
    ``expit(-0.5 + 1.2 (a - 12) - 0.03 (L - 15))``; every observable peptide of a gene in
    ``extra`` (undetected genes) is seen too."""
    from gippa.operons.detection_limit import _observable_seqs

    cfg = _cfg()
    c2 = chain.frame("genes")
    rng = np.random.default_rng(5)
    cols = [c.name for c in TABLES["peptide_hits"]]
    rows: list[list[Any]] = []
    for i, r in enumerate(c2.itertuples(index=False)):
        if not isinstance(r.aa, str) or (i % 4 == 3 and r.locus_tag not in extra):
            continue
        a = 12.0 + math.sin(i)
        for s in _observable_seqs(r.aa, cfg):
            if i % 4 == 3:
                pass  # an ``extra`` gene: every observable peptide seen
            elif rng.random() >= expit(-0.5 + 1.2 * (a - 12.0) - 0.03 * (len(s) - 15)):
                continue
            cell = {
                "peptide": s, "peptide_len": len(s), "source": None, "best_q_input": 0.001,
                "pep_input": None, "score_input": None, "contig": r.contig,
                "start": int(r.start) + 1, "end": int(r.end), "strand": r.strand, "frame": None,
                "crosses_origin": False, "route": "exact", "n_diff": 0, "diff_types": None,
                "placement_q": None, "length_class": None, "locus_class": "cds_inframe",
                "uniqueness": "unique", "gene": r.locus_tag, "novel_candidate_id": None,
                "flags": None,
            }  # fmt: skip
            rows.append([cell[c] for c in cols])
    if max_rows is not None:
        rows = rows[:max_rows]
    path.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(path, cols, rows)
    return path


def step_inputs(chain: Chain, gene_table: Path, extra: Mapping[str, Path]) -> dict[str, Path]:
    keys = ("genes", "replicons", "tus", "pairs_raw", "features", "promoter_windows",
            "terminators", "tx_rates")  # fmt: skip
    ins = {k: chain.inputs[k] for k in keys}
    ins["gene_table"] = gene_table
    ins.update(extra)
    return ins


class _Recorder(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


class Run:
    def __init__(self, out: Path, records: list[logging.LogRecord], inputs: dict[str, Path]):
        self.out, self.records, self.inputs = out, records, inputs
        self.paths = {k: out / "tables" / v for k, v in OUTPUTS.items()}

    def frame(self, key: str = "pairs") -> pd.DataFrame:
        return pd.read_csv(self.paths[key], sep="\t", dtype=str, keep_default_na=False)


def run_step(inputs: Mapping[str, Path], out: Path, cfg: GippaConfig) -> Run:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    outputs = {k: out / "tables" / v for k, v in OUTPUTS.items()}
    rec = _Recorder()
    log = logging.getLogger("gippa.steps.operons")
    old = log.level
    log.addHandler(rec)
    log.setLevel(logging.DEBUG)
    try:
        reg.run_step("operons", cfg, dict(inputs), outputs)
    finally:
        log.removeHandler(rec)
        log.setLevel(old)
    return Run(out, rec.records, dict(inputs))


@pytest.fixture(scope="module")
def with_peptides(chain: Chain, tmp_path_factory: pytest.TempPathFactory) -> Run:
    root = tmp_path_factory.mktemp("operons_pmiss")
    gt = write_gene_table(chain, root / "genetable" / "genes.parquet")
    hits = write_peptide_hits(chain, root / "mapping" / "peptide_hits.tsv")
    return run_step(step_inputs(chain, gt, {"peptide_hits": hits}), root / "out", _cfg())


@pytest.fixture(scope="module")
def without_peptides(chain: Chain, tmp_path_factory: pytest.TempPathFactory) -> Run:
    root = tmp_path_factory.mktemp("operons_nopep")
    gt = write_gene_table(chain, root / "genetable" / "genes.parquet")
    return run_step(step_inputs(chain, gt, {}), root / "out", _cfg())


def _flags(s: str) -> set[str]:
    return set() if s == "NA" else set(s.split(";"))


def _f(s: str) -> float:
    return math.nan if s == "NA" else float(s)


def _one_detected_op(df: pd.DataFrame) -> pd.Series:
    """OP rows with exactly one detected gene, neither partial (the rows that carried
    ``no_detection_limit_model`` in pass 1)."""
    op = df["class"] == "OP"
    return (
        op
        & (df["detected_a"] != df["detected_b"])
        & df["flags"].map(lambda s: "too_few_op_pairs" not in _flags(s))
    )


def _seen_counts(chain: Chain, hits: Path) -> dict[str, int]:
    """The test's oracle of M6-02's count: per gene, its distinct observable peptides with a C11
    assignment to it (I and L the same)."""
    from gippa.operons.detection_limit import _observable_seqs

    h = pd.read_csv(hits, sep="\t", dtype=str, keep_default_na=False)
    seen = {(g, q.replace("I", "L")) for g, q in zip(h["gene"], h["peptide"], strict=True)}
    out: dict[str, int] = {}
    for r in chain.frame("genes").itertuples(index=False):
        seqs = _observable_seqs(r.aa, _cfg()) if isinstance(r.aa, str) else []
        out[str(r.locus_tag)] = sum((str(r.locus_tag), q.replace("I", "L")) in seen for q in seqs)
    return out


def _seen_one_detected(df: pd.DataFrame, seen: Mapping[str, int], k_min: int = 1) -> pd.Series:
    """OP rows of the S1 event: exactly one gene with >= k_min observable peptides seen, no t
    test, not ``too_few_op_pairs``."""
    sa = df["gene_a"].map(lambda g: seen.get(g, 0) >= k_min)
    sb = df["gene_b"].map(lambda g: seen.get(g, 0) >= k_min)
    ok = df["flags"].map(lambda s: "too_few_op_pairs" not in _flags(s))
    return (df["class"] == "OP") & (sa != sb) & (df["test"] != TTEST) & ok


# ----------------------------------------------------------------------------- the live branch


def test_model_fixture_is_usable(chain: Chain, with_peptides: Run) -> None:
    """Guard on the test data: M6-02 fits a usable peptide model on these inputs."""
    from gippa.io.work_tables import _read_tsv_table
    from gippa.operons import detection_limit as dl
    from gippa.rng import for_step

    fit = dl.fit_detection_limit(
        _dl_genes(chain, with_peptides), _read_tsv_table("peptide_hits",
        with_peptides.inputs["peptide_hits"]), None, _cfg(), for_step(DL_STEP, _cfg().seed),
    )  # fmt: skip
    assert fit.model.usable, fit.model.reason


def test_pairs_tsv_has_one_detected_pmiss_rows(with_peptides: Run) -> None:
    df = with_peptides.frame()
    assert (df["test"] == PMISS).sum() > 0, df["test"].value_counts().to_dict()


def test_every_one_detected_op_pair_is_scored_or_flagged(with_peptides: Run, chain: Chain) -> None:
    """Every row of the S1 event is scored, or ``none`` + ``no_detection_limit_model`` (u has no
    observable peptide), and no row outside the event is scored."""
    df = with_peptides.frame()
    one = _seen_one_detected(df, _seen_counts(chain, with_peptides.inputs["peptide_hits"]))
    assert one.sum() > 0
    assert not (df.loc[~one, "test"] == PMISS).any(), df.loc[~one & (df["test"] == PMISS),
                                                            "pair_id"].tolist()  # fmt: skip
    sub = df[one]
    ok = (sub["test"] == PMISS) | (
        (sub["test"] == "none")
        & sub["flags"].map(lambda s: "no_detection_limit_model" in _flags(s))
    )
    assert ok.all(), sub.loc[~ok, ["pair_id", "test", "flags"]].to_string()
    assert (sub["test"] == PMISS).mean() > 0.5  # most undetected partners have peptides


def test_pmiss_rows_follow_item4(with_peptides: Run) -> None:
    df = with_peptides.frame()
    pm = df[df["test"] == PMISS]
    assert len(pm) > 0
    p = pm["p_value"].map(_f).to_numpy()
    assert np.isfinite(p).all() and (p >= 0).all() and (p <= 1).all()
    assert (pm["method"] == "detection_limit").all()
    assert (pm["null_model"] == "p_detect_model").all()
    assert (pm["effect"] == "NA").all()
    for r in pm.itertuples(index=False):
        a_v = r.abundance_log2_a if r.detected_a == "true" else r.abundance_log2_b
        assert _f(r.stat) == pytest.approx(_f(a_v), rel=1e-5), r.pair_id
        assert "no_detection_limit_model" not in _flags(r.flags), r.pair_id
    if not PMISS_IN_BH_FAMILY:
        assert (pm["q_value"] == "NA").all()


def test_bh_family_from_the_step(with_peptides: Run) -> None:
    df = with_peptides.frame()
    fam = df["test"].isin({TTEST, PMISS} if PMISS_IN_BH_FAMILY else {TTEST}) & (df["class"] == "OP")
    assert fam.sum() > 0
    p = df.loc[fam, "p_value"].map(_f).to_numpy()
    q = df.loc[fam, "q_value"].map(_f).to_numpy()
    np.testing.assert_allclose(q, bh(p), rtol=1e-4, atol=1e-6)  # %.6g in the TSV


def test_pairs_tsv_validates_with_the_pmiss_enum(with_peptides: Run) -> None:
    """§9.3: ``one_detected_pmiss`` is a value of the ``test`` enum and the table validates."""
    frictionless = pytest.importorskip("frictionless")
    d = json.loads((SCHEMA_DIR / "pairs.schema.json").read_text(encoding="utf-8"))
    test_field = next(f for f in d["fields"] if f["name"] == "test")
    assert PMISS in test_field["constraints"]["enum"]
    d = {k: v for k, v in d.items() if k != "coordinates"}
    path = with_peptides.paths["pairs"]
    res = frictionless.Resource(path=path.name, basepath=str(path.parent),
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])
    assert (with_peptides.frame()["test"] == PMISS).any()


def _dl_genes(chain: Chain, run: Run) -> pd.DataFrame:
    from gippa.io.genes_tsv import read_genes_parquet

    gt = read_genes_parquet(run.inputs["gene_table"])
    aa = chain.frame("genes").set_index("locus_tag")["aa"]
    out = gt[["locus_tag", "tu_id", "detected", "abundance_log2", "detection_tier"]].copy()
    out["aa"] = out["locus_tag"].map(aa)
    return out


def test_pmiss_values_equal_the_m6_02_oracle(chain: Chain, with_peptides: Run) -> None:
    """Each P(miss) is ``p_miss(model_for_tu(model, tu(v)), a_v, r, observable_peptides(u))``
    with r = ±Δ of the fit pairs (the both-detected OP rows here; no gene is imputed) and v the
    detected partner (M6-02 lead decision 9, review B1)."""
    from gippa.io.work_tables import _read_tsv_table
    from gippa.operons import detection_limit as dl
    from gippa.rng import for_step

    cfg = _cfg()
    genes = _dl_genes(chain, with_peptides)
    hits = _read_tsv_table("peptide_hits", with_peptides.inputs["peptide_hits"])
    model = dl.fit_detection_limit(genes, hits, None, cfg, for_step(DL_STEP, cfg.seed)).model
    df = with_peptides.frame()
    delta = df.loc[(df["class"] == "OP") & (df["test"] == TTEST), "delta_log2"].map(_f).to_numpy()
    tu, aa = genes.set_index("locus_tag")["tu_id"], genes.set_index("locus_tag")["aa"]
    pm = df[df["test"] == PMISS]
    assert len(pm) > 0
    for r in pm.itertuples(index=False):
        want, _ = _oracle(r, model, delta, tu, aa, cfg)
        assert _f(r.p_value) == pytest.approx(want, abs=1e-4), r.pair_id


def _oracle(
    r: Any, model: Any, delta: np.ndarray, tu: pd.Series, aa: pd.Series, cfg: GippaConfig
) -> tuple[float, float]:
    """(P(miss) with v's TU model, P(miss) with u's TU model) for a pairs.tsv row."""
    from gippa.operons import detection_limit as dl

    down = r.detected_a == "true"  # gene_b undetected: u downstream, r = Δ
    u, v = (r.gene_b, r.gene_a) if down else (r.gene_a, r.gene_b)
    a_v = _f(r.abundance_log2_a if down else r.abundance_log2_b)
    ratios = delta if down else -delta
    peps = dl.observable_peptides(str(aa[u]), cfg)
    return (
        dl.p_miss(dl.model_for_tu(model, str(tu[v])), a_v, ratios, peps),
        dl.p_miss(dl.model_for_tu(model, str(tu[u])), a_v, ratios, peps),
    )


def test_pmiss_uses_the_partners_tu_when_the_pair_crosses_tus(chain: Chain, tmp_path: Path) -> None:
    """B1 at step level: every undetected gene sits in a TU of its own, so each one-detected pair
    crosses a TU boundary and u's TU trained nothing (``model_for_tu`` falls back to the full
    model). The P(miss) are those of v's out-of-fold model, and the two readings differ on at
    least one row (so the test tells them apart)."""
    from gippa.io.work_tables import _read_tsv_table
    from gippa.operons import detection_limit as dl
    from gippa.rng import for_step

    cfg = _cfg()
    gt = write_gene_table(chain, tmp_path / "genetable" / "genes.parquet", split_tus=True)
    hits = write_peptide_hits(chain, tmp_path / "mapping" / "peptide_hits.tsv")
    run = run_step(step_inputs(chain, gt, {"peptide_hits": hits}), tmp_path / "out", cfg)
    genes = _dl_genes(chain, run)
    model = dl.fit_detection_limit(genes, _read_tsv_table("peptide_hits", hits), None, cfg,
                                   for_step(DL_STEP, cfg.seed)).model  # fmt: skip
    assert model.usable and model.fold_coef
    df = run.frame()
    delta = df.loc[(df["class"] == "OP") & (df["test"] == TTEST), "delta_log2"].map(_f).to_numpy()
    tu, aa = genes.set_index("locus_tag")["tu_id"], genes.set_index("locus_tag")["aa"]
    pm = df[df["test"] == PMISS]
    assert len(pm) > 0
    differ = 0
    for r in pm.itertuples(index=False):
        want, u_reading = _oracle(r, model, delta, tu, aa, cfg)
        assert _f(r.p_value) == pytest.approx(want, abs=1e-4), r.pair_id
        differ += abs(want - u_reading) > 1e-3
    assert differ > 0


def test_same_inputs_give_byte_identical_tables(with_peptides: Run, tmp_path: Path) -> None:
    again = run_step(with_peptides.inputs, tmp_path / "again", _cfg())
    for key in OUTPUTS:
        assert again.paths[key].read_bytes() == with_peptides.paths[key].read_bytes(), key


def test_usable_model_drops_the_no_model_run_flag(with_peptides: Run) -> None:
    flags = [getattr(r, "data", {}).get("flag") for r in with_peptides.records
             if getattr(r, "event", None) == "operons.run_flag"]  # fmt: skip
    assert "no_detection_limit_model" not in flags


# ----------------------------------------------------------------------------- unestimated paths


def test_no_peptides_one_detected_rows_are_unestimated(without_peptides: Run) -> None:
    """M6-08 degenerate row / §7.4.3: one run and no peptides → ``none`` +
    ``no_detection_limit_model`` (regression guard: pass-1 behaviour must survive)."""
    df = without_peptides.frame()
    sub = df[_one_detected_op(df)]
    assert len(sub) > 0
    assert (sub["test"] == "none").all()
    assert sub["flags"].map(lambda s: "no_detection_limit_model" in _flags(s)).all()
    assert (df["test"] != PMISS).all()


def test_unusable_model_one_detected_rows_are_unestimated(chain: Chain, tmp_path: Path) -> None:
    """Too few peptide hits (< 50 training rows) → model not usable → no P(miss) anywhere."""
    gt = write_gene_table(chain, tmp_path / "genetable" / "genes.parquet")
    hits = write_peptide_hits(chain, tmp_path / "mapping" / "peptide_hits.tsv", max_rows=5)
    run = run_step(step_inputs(chain, gt, {"peptide_hits": hits}), tmp_path / "out", _cfg())
    df = run.frame()
    sub = df[_one_detected_op(df)]
    assert len(sub) > 0
    assert (df["test"] != PMISS).all()
    assert sub["flags"].map(lambda s: "no_detection_limit_model" in _flags(s)).all()
    assert any("no_detection_limit_model" in r.getMessage() for r in run.records)


# ----------------------------------------------------------------------------- review of #311


def test_undetected_gene_with_peptides_seen_is_not_one_detected(
    chain: Chain, with_peptides: Run, tmp_path: Path
) -> None:
    """S1 at step level (M6-02 lead decision 9): give C11 hits to a few genes that are not
    ``detected`` in the gene table and whose pairs were scored. With >= k_min (1) of their
    observable peptides seen they were not missed, so those pairs are no longer
    ``one_detected_pmiss``; the other scored pairs are unchanged."""
    base = with_peptides.frame()
    pm = base[base["test"] == PMISS]
    u = [b if a == "true" else g for a, b, g in zip(pm["detected_a"], pm["gene_b"], pm["gene_a"],
                                                  strict=True)]  # fmt: skip
    chosen = frozenset(sorted(set(u))[:3])
    assert chosen
    gt = write_gene_table(chain, tmp_path / "genetable" / "genes.parquet")
    hits = write_peptide_hits(chain, tmp_path / "mapping" / "peptide_hits.tsv", extra=chosen)
    seen = _seen_counts(chain, hits)
    assert all(seen[g] >= 1 for g in chosen)
    df = run_step(step_inputs(chain, gt, {"peptide_hits": hits}), tmp_path / "out", _cfg()).frame()
    touched = df["gene_a"].isin(chosen) | df["gene_b"].isin(chosen)
    by_id = df.set_index("pair_id")
    for pid in base.loc[(base["test"] == PMISS) & touched.to_numpy(), "pair_id"]:
        assert by_id.loc[pid, "test"] != PMISS, pid
    still = base.loc[(base["test"] == PMISS) & ~touched.to_numpy(), "pair_id"]
    assert len(still) > 0
    assert (by_id.loc[still, "test"] == PMISS).all()


COV_FLAG = "pmiss_covariates_missing"


def _replicate_fit(coef: Mapping[str, float]) -> Any:
    from gippa.genetable.detectability import COLUMNS
    from gippa.operons import detection_limit as dl

    model = dl.PMissModel(kind="replicate", coef=dict(coef), k_min=1, n_fit=600, auc=0.75,
                          usable=True)  # fmt: skip
    summary = pd.DataFrame(columns=COLUMNS).astype(
        {"value": "float64", "ci_low": "float64", "ci_high": "float64"}
    )
    return dl.DetectionLimitFit(model, None, model, summary)


def write_peptide_runs(hits: Path, path: Path) -> Path:
    """C16 ``peptide_runs.parquet``: each C11 peptide in 1 to 3 of three runs (test data)."""
    h = pd.read_csv(hits, sep="\t", dtype=str, keep_default_na=False)
    peps = sorted(set(h["peptide"]))
    rows = [(q, f"run{j}", 1000.0 * (j + 1), 1) for i, q in enumerate(peps)
            for j in range(1 + i % 3)]  # fmt: skip
    df = pd.DataFrame(rows, columns=["peptide", "run", "intensity", "n_rows"])
    df["n_rows"] = df["n_rows"].astype("int64")
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)
    return path


def _replicate_run(chain: Chain, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
                   coef: Mapping[str, float]) -> tuple[Run, dict[str, int]]:  # fmt: skip
    from gippa.operons import detection_limit as dl

    monkeypatch.setattr(dl, "fit_detection_limit", lambda *a, **k: _replicate_fit(coef))
    gt = write_gene_table(chain, tmp_path / "genetable" / "genes.parquet")
    hits = write_peptide_hits(chain, tmp_path / "mapping" / "peptide_hits.tsv")
    runs = write_peptide_runs(hits, tmp_path / "work" / "peptide_runs.parquet")
    ins = step_inputs(chain, gt, {"peptide_hits": hits, "peptide_runs": runs})
    return run_step(ins, tmp_path / "out", _cfg()), _seen_counts(chain, hits)


def test_replicate_model_needing_covariates_is_flagged_not_silently_nan(
    chain: Chain, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """S2: the model in force is the replicate model with a ``detectability.covariates``
    coefficient (``length_aa``), which the step does not join. The one-detected rows are
    ``none`` with ``pmiss_covariates_missing`` (not a silent NaN → ``no_detection_limit_model``),
    no P(miss) is written, the step logs a WARNING naming the flag, and pairs.tsv validates."""
    run, seen = _replicate_run(
        chain, tmp_path, monkeypatch, {"intercept": -9.0, "abundance_log2": 0.8, "length_aa": 0.1}
    )
    df = run.frame()
    one = _seen_one_detected(df, seen)
    assert one.sum() > 0
    assert (df["test"] != PMISS).all()
    sub = df[one]
    assert (sub["test"] == "none").all()
    assert sub["p_value"].eq("NA").all() and sub["q_value"].eq("NA").all()
    assert sub["flags"].map(lambda s: COV_FLAG in _flags(s)).all(), sub[["pair_id", "flags"]]
    t = df[df["test"] == TTEST]
    assert len(t) > 0 and not t["flags"].map(lambda s: COV_FLAG in _flags(s)).any()
    warned = [r for r in run.records if r.levelno >= logging.WARNING and COV_FLAG in r.getMessage()]
    assert warned, [r.getMessage() for r in run.records]
    frictionless = pytest.importorskip("frictionless")
    d = json.loads((SCHEMA_DIR / "pairs.schema.json").read_text(encoding="utf-8"))
    d = {k: v for k, v in d.items() if k != "coordinates"}
    pth = run.paths["pairs"]
    res = frictionless.Resource(path=pth.name, basepath=str(pth.parent),
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_replicate_model_without_covariates_is_scored(
    chain: Chain, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """S2 control: a replicate model on abundance alone needs no covariate, so it scores the
    one-detected rows and nothing is flagged ``pmiss_covariates_missing``."""
    run, _ = _replicate_run(chain, tmp_path, monkeypatch,
                            {"intercept": -9.0, "abundance_log2": 0.8})  # fmt: skip
    df = run.frame()
    assert (df["test"] == PMISS).sum() > 0
    assert not df["flags"].map(lambda s: COV_FLAG in _flags(s)).any()
    assert not any(COV_FLAG in r.getMessage() for r in run.records)
