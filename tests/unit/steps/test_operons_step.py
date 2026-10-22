"""M6-08 the ``operons`` step on the sim100k producer chain
(docs/specs/M6-08-operons-step-and-run.md;
examples 1, 3, 4 and 6 at step level, degenerate cases, lead decisions 3, 4 and 6; GIAP_SPEC.md
§6.2 step row, §6.5, §6.7, §9.3).

Interface pinned here (test-author reading of the spec):

* ``@step("operons", inputs={gene_table, genes, replicons, tus, pairs_raw, features},
  optional_inputs={promoter_windows, terminators, tx_rates, peptides, peptide_hits,
  peptide_runs, detection_summary}, outputs={pairs, operons, operon_flags})``. The module is
  core, the kind is py, and no precomputed layer replaces the step.
* The step chains M6-04 ``pair_frame``, M6-05 ``discordance`` (``p_miss=None`` in pass 1), M6-06
  ``missed_operons`` (always ``genes=``), and M6-07 ``merge_tests``, ``flag_pairs`` and
  ``operons_table`` (lead decision 4). It writes the results with ``giap.io.operon_tables``.
  The oracle below calls exactly that chain on the step's inputs, so the step's tables must equal
  the oracle's byte for byte.
* ``missing_layers`` names the absent optional layers the way ``context_layers`` does:
  ``promoter_windows`` gives ``promoters``, ``terminators`` gives ``terminators`` and
  ``tx_rates`` gives ``tx_rates``. Each absent layer is also logged as the run flag
  ``layer_missing:<layer>``.
* Pass 1 logs ``no_detection_limit_model`` (spec item 1).
* Lead decision 6 (π₀ and the decoy ``bin_counts`` in the manifest). No ``run_manifest.json``
  writer exists yet. §6.7 events feed it, so the step records a §6.7 event
  ``operons.missed_operons`` through ``ctx.log`` (``extra={"event": ..., "data": ...}``) whose
  ``data`` holds ``pi0`` (``cfg.operons.pi0``) and ``bin_counts`` (``MissedOperons.bin_counts``
  as a list of per-bin records). ``data`` must be plain JSON.

Inputs come from ``helpers.layers_chain.build_chain`` (the merged producers on giapsim SIM100K,
with recorded UniOP, TerminatorNet and TransTermHP outputs). The gene table is written with
``giap.io.genes_tsv.write_genes_parquet``. Three genes in four are detected, with finite abundances
when the run has abundance. With abundance, >= 20 OP pairs are usable and ``operon_flags.tsv``
has rows. Without abundance, every test is ``none`` and ``operon_flags.tsv`` is header-only
(example 3, lead decision 3).
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

from giap.config import GiapConfig
from giap.operons.flags import FLAG_COLUMNS
from giap.operons.pairs import PAIR_COLUMNS, RNA_TYPES
from giap.operons.tu_table import OPERON_COLUMNS

SCHEMA_DIR = Path(__file__).resolve().parents[3] / "src" / "giap" / "schemas" / "tables"
REQUIRED = {"gene_table", "genes", "replicons", "tus", "pairs_raw", "features"}
OPTIONAL = {"promoter_windows", "terminators", "tx_rates", "peptides", "peptide_hits",
            "peptide_runs", "detection_summary"}  # fmt: skip
OUTPUTS = {"pairs": "pairs.tsv", "operons": "operons.tsv", "operon_flags": "operon_flags.tsv"}
LAYER_NAME = {"promoter_windows": "promoters", "terminators": "terminators", "tx_rates": "tx_rates"}
EVENT = "operons.missed_operons"


def _cfg() -> GiapConfig:
    return GiapConfig(seed=1)


# ----------------------------------------------------------------------------- fixtures


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(Path(tmp_path_factory.mktemp("operons_chain")))


def write_gene_table(chain: Chain, path: Path, *, abundance: bool) -> Path:
    """A §9.3 gene table over the C2 genes (every ``GENES_COLUMNS`` column)."""
    from giap.genetable.assemble import GENES_COLUMNS, GENES_TYPES
    from giap.io.genes_tsv import write_genes_parquet

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
    if abundance:
        df["abundance_log2"] = np.where(det, 12.0 + np.sin(idx), np.nan)
        df["abundance_type"] = "ibaq"
        df["abundance_scale"] = "absolute"
    df["abundance_imputed"] = False
    path.parent.mkdir(parents=True, exist_ok=True)
    write_genes_parquet(path, df)
    return path


def step_inputs(chain: Chain, gene_table: Path, *, drop: tuple[str, ...] = ()) -> dict[str, Path]:
    keys = ("genes", "replicons", "tus", "pairs_raw", "features", "promoter_windows",
            "terminators", "tx_rates")  # fmt: skip
    ins = {k: chain.inputs[k] for k in keys if k not in drop}
    ins["gene_table"] = gene_table
    return ins


class Run:
    def __init__(self, out: Path, records: list[logging.LogRecord], inputs: dict[str, Path]):
        self.out = out
        self.records = records
        self.inputs = inputs
        self.paths = {k: out / "tables" / v for k, v in OUTPUTS.items()}

    def frame(self, key: str) -> pd.DataFrame:
        return pd.read_csv(self.paths[key], sep="\t", dtype=str, keep_default_na=False)

    def messages(self) -> list[str]:
        return [r.getMessage() for r in self.records]


def run_operons_step(inputs: Mapping[str, Path], out: Path, cfg: GiapConfig) -> Run:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    outputs = {k: out / "tables" / v for k, v in OUTPUTS.items()}
    rec = _Recorder()
    log = logging.getLogger("giap.steps.operons")
    old = log.level
    log.addHandler(rec)
    log.setLevel(logging.DEBUG)
    try:
        reg.run_step("operons", cfg, dict(inputs), outputs)
    finally:
        log.removeHandler(rec)
        log.setLevel(old)
    return Run(out, rec.records, dict(inputs))


class _Recorder(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


def oracle(inputs: Mapping[str, Path], cfg: GiapConfig, out: Path) -> dict[str, Any]:
    """Lead decision 4's chain on the step's own inputs, written with the M6-08 writers."""
    from giap.annotate.layers import genes_of, replicons_of, terminators_of
    from giap.io.genes_tsv import read_genes_parquet
    from giap.io.operon_tables import write_operon_flags_tsv, write_operons_tsv, write_pairs_tsv
    from giap.operons.discordance import discordance
    from giap.operons.flags import flag_pairs, merge_tests
    from giap.operons.missed import missed_operons
    from giap.operons.pairs import Evidence, pair_frame
    from giap.operons.tu_table import operons_table

    rd = pd.read_parquet
    gt = read_genes_parquet(inputs["gene_table"])
    feats = rd(inputs["features"])
    ev = Evidence(
        windows=rd(inputs["promoter_windows"]) if "promoter_windows" in inputs else None,
        terminators=(tuple(terminators_of(rd(inputs["terminators"])))
                     if "terminators" in inputs else None),
        tx=rd(inputs["tx_rates"]) if "tx_rates" in inputs else None,
        rna=feats[feats["type"].isin(RNA_TYPES)].reset_index(drop=True),
    )  # fmt: skip
    pairs = pair_frame(gt, genes_of(rd(inputs["genes"])), replicons_of(rd(inputs["replicons"])),
                       rd(inputs["pairs_raw"]), ev, cfg)  # fmt: skip
    d = discordance(pairs, gt, cfg, None)
    m = missed_operons(pairs, cfg, genes=gt)
    missing = [LAYER_NAME[k] for k in LAYER_NAME if k not in inputs]
    merged = merge_tests(pairs, d.rows, m.rows, gt, missing_layers=missing)
    flagged = flag_pairs(merged, cfg)
    ops = operons_table(rd(inputs["tus"]), merged, gt, flagged)
    out.mkdir(parents=True, exist_ok=True)
    paths = {k: out / v for k, v in OUTPUTS.items()}
    write_pairs_tsv(paths["pairs"], merged)
    write_operons_tsv(paths["operons"], ops)
    write_operon_flags_tsv(paths["operon_flags"], flagged)
    return {"paths": paths, "missed": m, "merged": merged, "flagged": flagged}


def _validate(path: Path, stem: str) -> Any:
    frictionless = pytest.importorskip("frictionless")
    d = json.loads((SCHEMA_DIR / f"{stem}.schema.json").read_text(encoding="utf-8"))
    d = {k: v for k, v in d.items() if k != "coordinates"}
    res = frictionless.Resource(path=path.name, basepath=str(path.parent),
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    return res.validate()


@pytest.fixture(scope="module")
def full(chain: Chain, tmp_path_factory: pytest.TempPathFactory) -> Run:
    """Example 1 at step level: every layer, with abundance."""
    root = tmp_path_factory.mktemp("operons_full")
    gt = write_gene_table(chain, root / "genetable" / "genes.parquet", abundance=True)
    return run_operons_step(step_inputs(chain, gt), root / "out", _cfg())


@pytest.fixture(scope="module")
def no_abundance(chain: Chain, tmp_path_factory: pytest.TempPathFactory) -> Run:
    """Example 3: no abundance."""
    root = tmp_path_factory.mktemp("operons_noab")
    gt = write_gene_table(chain, root / "genetable" / "genes.parquet", abundance=False)
    return run_operons_step(step_inputs(chain, gt), root / "out", _cfg())


@pytest.fixture(scope="module")
def no_terminators(chain: Chain, tmp_path_factory: pytest.TempPathFactory) -> Run:
    """Example 4: no terminator layer."""
    root = tmp_path_factory.mktemp("operons_noterm")
    gt = write_gene_table(chain, root / "genetable" / "genes.parquet", abundance=True)
    ins = step_inputs(chain, gt, drop=("terminators",))
    return run_operons_step(ins, root / "out", _cfg())


# ----------------------------------------------------------------------------- registry


def test_operons_step_is_registered_with_the_spec_keys() -> None:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    spec = reg.get("operons")
    assert spec.kind == "py"
    assert spec.module == "core"  # lead decision 2: no --modules flag of its own
    assert spec.precomputed is None
    assert set(spec.inputs) == REQUIRED
    assert set(spec.optional_inputs) == OPTIONAL
    assert set(spec.outputs) == set(OUTPUTS)


def test_operons_runs_after_its_producers_in_the_local_plan() -> None:
    """The local executor orders ``operons`` after ``detectability`` (``gene_table``) and
    ``operon_pairs`` (``pairs_raw``, ``tus``)."""
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    names = [s.name for s in reg.all_steps()]
    order = reg.topo_order(names)
    for producer in ("detectability", "operon_pairs", "parse_bakta", "terminator_layer",
                     "promoter_layer", "tx_layer"):  # fmt: skip
        assert order.index(producer) < order.index("operons"), producer


# ----------------------------------------------------------------------------- example 1


def test_example_1_three_tables_written(full: Run) -> None:
    for key, path in full.paths.items():
        assert path.is_file(), f"{key} not written"
    assert full.frame("pairs").columns.tolist() == list(PAIR_COLUMNS)
    assert full.frame("operons").columns.tolist() == list(OPERON_COLUMNS)
    assert full.frame("operon_flags").columns.tolist() == list(FLAG_COLUMNS)


def test_example_1_tables_equal_the_decision_4_chain(full: Run, tmp_path: Path) -> None:
    want = oracle(full.inputs, _cfg(), tmp_path / "oracle")
    for key, path in full.paths.items():
        assert path.read_bytes() == want["paths"][key].read_bytes(), key
    assert len(want["flagged"]) > 0, "fixture meant to give flagged pairs"


def test_pairs_tsv_has_no_which_column(full: Run, no_abundance: Run) -> None:
    for run in (full, no_abundance):
        header = run.paths["pairs"].read_text(encoding="utf-8").split("\n", 1)[0].split("\t")
        assert "which" not in header


def test_pairs_one_row_per_adjacent_pair_plus_sandwiches(chain: Chain, full: Run) -> None:
    from giap.annotate.layers import genes_of, replicons_of
    from giap.operons.pairs import adjacent_chain, sandwich_pairs

    genes = genes_of(chain.frame("genes"))
    reps = replicons_of(chain.frame("replicons"))
    n_adj, n_sw = len(adjacent_chain(genes, reps)), len(sandwich_pairs(genes, reps))
    df = full.frame("pairs")
    assert len(df) == n_adj + n_sw
    assert int((df["class"] == "SW").sum()) == n_sw
    assert df["pair_id"].is_unique
    nums = [int(p.rsplit("_", 1)[1]) for p in df["pair_id"]]
    assert nums == list(range(1, len(df) + 1))


def test_operons_coordinates_are_one_based_closed(chain: Chain, full: Run) -> None:
    """§6.5: the TSV start is C5 ``start + 1``; ``end`` is unchanged."""
    tus = chain.frame("tus").set_index("tu_id")
    df = full.frame("operons")
    assert set(df["tu_id"]) == set(tus.index)
    for r in df.itertuples(index=False):
        assert int(r.start) == int(tus.at[r.tu_id, "start"]) + 1, r.tu_id
        assert int(r.end) == int(tus.at[r.tu_id, "end"]), r.tu_id
        assert r.genes.split(";") == list(tus.at[r.tu_id, "genes"]), r.tu_id


def test_example_6_three_outputs_validate(full: Run, no_abundance: Run) -> None:
    for run in (full, no_abundance):
        for key in OUTPUTS:
            report = _validate(run.paths[key], key)
            assert report.valid, (key, report.flatten(["rowNumber", "fieldName", "type", "note"]))


def test_same_inputs_give_byte_identical_tables(chain: Chain, full: Run, tmp_path: Path) -> None:
    again = run_operons_step(full.inputs, tmp_path / "again", _cfg())
    for key in OUTPUTS:
        assert again.paths[key].read_bytes() == full.paths[key].read_bytes(), key


def test_pass_1_logs_no_detection_limit_model(full: Run) -> None:
    assert any("no_detection_limit_model" in m for m in full.messages()), full.messages()


# ----------------------------------------------------------------------------- decision 6


def _event(run: Run) -> dict[str, Any]:
    found = [r for r in run.records if getattr(r, "event", None) == EVENT]
    assert len(found) == 1, f"expected one {EVENT} event, got {len(found)}"
    data = getattr(found[0], "data", None)
    assert isinstance(data, dict), data
    return data


def test_decision_6_event_records_pi0_and_bin_counts(full: Run, tmp_path: Path) -> None:
    data = _event(full)
    assert data["pi0"] == _cfg().operons.pi0
    want = oracle(full.inputs, _cfg(), tmp_path / "oracle")["missed"].bin_counts
    got = pd.DataFrame(data["bin_counts"])
    assert got.columns.tolist() == want.columns.tolist()
    assert len(got) == len(want) == len(_cfg().operons.distance_bins_bp) - 1
    for col in want.columns:
        assert [float(x) for x in got[col]] == [float(x) for x in want[col]], col
    json.dumps(data)  # §6.7 events are plain JSON (no numpy scalars)


def test_decision_6_pi0_follows_the_config(chain: Chain, full: Run, tmp_path: Path) -> None:
    cfg = GiapConfig.model_validate({"seed": 1, "operons": {"pi0": 0.8}})
    run = run_operons_step(full.inputs, tmp_path / "pi0", cfg)
    assert _event(run)["pi0"] == pytest.approx(0.8)


# ----------------------------------------------------------------------------- example 3


def test_example_3_no_abundance_every_test_none(no_abundance: Run) -> None:
    df = no_abundance.frame("pairs")
    assert len(df) > 0
    assert set(df["test"]) == {"none"}
    for col in ("stat", "effect", "ci_low", "ci_high", "p_value", "q_value", "lfdr", "n_perm"):
        assert set(df[col]) == {"NA"}, col


def test_example_3_operon_flags_header_only(no_abundance: Run) -> None:
    text = no_abundance.paths["operon_flags"].read_text(encoding="utf-8")
    assert text == "\t".join(FLAG_COLUMNS) + "\n"


def test_example_3_too_few_pairs_is_the_normal_outcome(no_abundance: Run) -> None:
    """Degenerate case "sim100k": ``too_few_op_pairs`` / ``too_few_no_pairs``, and the step
    passes."""
    df = no_abundance.frame("pairs")
    flags = df["flags"].map(lambda s: set() if s == "NA" else set(s.split(";")))
    op, no = df["class"] == "OP", df["class"] == "NO"
    assert flags[op].map(lambda f: "too_few_op_pairs" in f).any()
    assert flags[no].map(lambda f: "too_few_no_pairs" in f).all()


def test_example_3_equals_the_decision_4_chain(no_abundance: Run, tmp_path: Path) -> None:
    want = oracle(no_abundance.inputs, _cfg(), tmp_path / "oracle")
    for key, path in no_abundance.paths.items():
        assert path.read_bytes() == want["paths"][key].read_bytes(), key


# ----------------------------------------------------------------------------- example 4


def test_example_4_terminator_between_na_everywhere(no_terminators: Run) -> None:
    df = no_terminators.frame("pairs")
    assert set(df["terminator_between"]) == {"NA"}


def test_example_4_layer_missing_on_op_and_no_rows_only(no_terminators: Run) -> None:
    df = no_terminators.frame("pairs")
    has = df["flags"].map(lambda s: "layer_missing:terminators" in s.split(";"))
    op_no = df["class"].isin(["OP", "NO"])
    assert has[op_no].all()
    assert not has[~op_no].any()


def test_example_4_run_flag_logged(no_terminators: Run) -> None:
    hits = [r for r in no_terminators.records
            if "layer_missing:terminators" in r.getMessage()]  # fmt: skip
    assert hits, no_terminators.messages()
    assert all(r.levelno >= logging.WARNING for r in hits)


def test_example_4_equals_the_decision_4_chain(no_terminators: Run, tmp_path: Path) -> None:
    want = oracle(no_terminators.inputs, _cfg(), tmp_path / "oracle")
    for key, path in no_terminators.paths.items():
        assert path.read_bytes() == want["paths"][key].read_bytes(), key


def test_all_three_layers_missing_names_each(chain: Chain, tmp_path: Path) -> None:
    gt = write_gene_table(chain, tmp_path / "genetable" / "genes.parquet", abundance=True)
    ins = step_inputs(chain, gt, drop=tuple(LAYER_NAME))
    run = run_operons_step(ins, tmp_path / "out", _cfg())
    df = run.frame("pairs")
    for col in ("promoter_between", "tx_rate_between", "terminator_between"):
        assert set(df[col]) == {"NA"}, col
    want = {f"layer_missing:{v}" for v in LAYER_NAME.values()}
    op_no = df.loc[df["class"].isin(["OP", "NO"]), "flags"]
    assert op_no.map(lambda s: want <= set(s.split(";"))).all()
    logged = " ".join(run.messages())
    for flag in want:
        assert flag in logged, flag
    assert run.paths["operon_flags"].is_file()
