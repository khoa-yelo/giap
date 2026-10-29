"""Mapping steps, query selection and routing (issue M3-17; GIPPA_SPEC.md §6.1, §6.2, §7.2,
§7.3.1; docs/specs/M3-17-mapping-steps.md examples 1 and 6, plus step I/O and determinism).

The registered steps run through ``gippa.steps.registry.run_step`` (M2-02) on a 50-kb gippasim
genome with M3-19 proteomics after ``gippa preflight`` (M1-12 work files); external tools never
run: miniprot / tblastn outputs are injected from ``tests/data/recorded`` (helpers.mapping_steps).
Each step's output is compared with the merged pure functions (M3-02 … M3-16) on the same input.

Ambiguities resolved here (spec silent; see helpers.mapping_steps for the key conventions):

* ``select_queries`` returns a frame with one row per query and a ``peptide`` column holding the
  input text; ``source`` in the result is resolved (never ``auto``): Q3 default, *table-level* —
  ``database`` when any row of the table has a ``q_value`` or ``pep``, else ``denovo``.
* Q2 (lead decision 4, provisional): q exactly at ``input.q_value_max`` is mapped (``≤``).
* "Run flag ``c_assumed`` recorded" (example 6) = the ``denovo_precision`` rows of
  ``mapping_fdr.tsv`` carry ``c_assumed`` in ``flags`` and a mapping step logs a WARNING that
  names ``c_assumed``.
* "Every step logs counts": each mapping py step emits at least one INFO record with a number on
  its ``gippa.steps.<name>`` logger (``run_step``'s logger; ``logs/<step>.log`` is the CLI's sink).
* The shuffled decoy translations in the ``sixframe`` work file are
  ``build_decoys(index, "shuffle", mapping.decoys.shuffle_window_aa,
  for_step("sixframe_index", seed))`` (the step's own stream, §6.6).
"""

from __future__ import annotations

import logging
import math
import re
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.mapping_steps import (
    MAPPING_EXT_STEPS,
    MAPPING_PY_STEPS,
    MAPPING_TSVS,
    RECORDED_MINIPROT,
    RECORDED_TBLASTN,
    SimRun,
    drive,
    fresh_steps,
    key_like,
    mapping_module,
    module2_files,
    proteins_key,
    read_tsv_df,
    roots_of,
    simulate_with_proteomics,
)
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.annotate.bakta import parse_bakta_json
from gippa.config import GippaConfig, resolve_config
from gippa.genetable.contract import validate_table
from gippa.io.work_tables import _read_tsv_table
from gippa.mapping.chance import residue_freqs
from gippa.mapping.decoys import build_decoys
from gippa.mapping.exact import match_exact
from gippa.mapping.frameshift import match_frameshift
from gippa.mapping.placement_fdr import best_costs, placement
from gippa.mapping.sixframe import build_sixframe, load_index, load_start_cells
from gippa.mapping.tolerant import match_tolerant
from gippa.rng import for_step

HIT_KEY = ["peptide", "search", "decoy_type", "contig", "start", "end", "strand", "n_diff"]
SEED = 11
AA = "ACDEFGHKLMNPQRSTVWY"


@pytest.fixture
def steps() -> Iterator[ModuleType]:
    with pytest.MonkeyPatch.context() as mp:
        yield fresh_steps(mp)


# ---------------------------------------------------------------- peptide tables


def norm_table(rows: list[dict[str, Any]]) -> pd.DataFrame:
    """A ``peptides.norm.parquet``-shaped frame (M1-12 / contract C9 columns)."""
    base = {
        "q_value": np.nan,
        "pep": np.nan,
        "score": np.nan,
        "is_decoy": False,
        "is_contaminant": False,
        "protein_ids": ["P1"],
        "modified_peptides": [],
        "n_rows": 1,
        "n_runs": 1,
        "intensity": np.nan,
        "source": "database",
        "search_db": "other",
    }
    df = pd.DataFrame([{**base, **r} for r in rows])
    for c in ("q_value", "pep", "score", "intensity"):
        df[c] = df[c].astype("float64")
    df["is_decoy"] = df["is_decoy"].astype(bool)
    df["source"] = df["source"].astype("string")
    df["search_db"] = df["search_db"].astype("string")
    return df


def expected_queries(df: pd.DataFrame, cfg: GippaConfig) -> set[str]:
    """Oracle of the spec's query rule (Behaviour, Q2)."""
    lo, hi = cfg.mapping.exact.min_len, cfg.mapping.exact.max_len
    qmax = cfg.input.q_value_max
    out = set()
    for r in df.itertuples(index=False):
        q = r.q_value
        if bool(r.is_decoy) or not lo <= len(r.peptide) <= hi:
            continue
        if q is None or (isinstance(q, float) and math.isnan(q)) or q <= qmax:
            out.add(str(r.peptide))
    return out


def select(df: pd.DataFrame, cfg: GippaConfig | None = None) -> pd.DataFrame:
    out = mapping_module().select_queries(df, cfg or GippaConfig())
    assert isinstance(out, pd.DataFrame)
    assert "peptide" in out.columns
    return out


# ---------------------------------------------------------------- declarations (§6.2)


def test_mapping_steps_are_registered(steps: ModuleType) -> None:
    for name in MAPPING_PY_STEPS:
        s = steps.get(name)
        assert (s.kind, s.module) == ("py", "core"), name
        assert s.func is not None
    for name in MAPPING_EXT_STEPS:
        s = steps.get(name)
        assert (s.kind, s.module, s.func) == ("ext", "core", None), name
        assert "genome" in s.inputs, name


def test_step_keys_follow_the_section_6_2_files(steps: ModuleType) -> None:
    six = key_like(steps.get("sixframe_index").outputs, "sixframe")
    assert "genome" in steps.get("sixframe_index").inputs
    ex = steps.get("map_exact")
    assert {"peptides", six} <= set(ex.inputs)
    hx = key_like(ex.outputs, "hits_exact")
    tol = steps.get("map_tolerant")
    # leftovers need the exact hits; the index is searched again
    assert {"peptides", six, hx} <= set(tol.inputs)
    ht = key_like(tol.outputs, "hits_tolerant")
    pl = steps.get("placement")
    assert {hx, ht} <= set(pl.inputs) | set(pl.optional_inputs)
    assert "peptides" in set(pl.inputs) | set(pl.optional_inputs)  # input q / PEP
    key_like(pl.outputs, "placement")
    key_like(pl.outputs, "mapping_fdr")
    key_like(pl.outputs, "stratum0_classes")  # ADR-0019 class table (M3-25 decision 3)
    lo = steps.get("loci")
    assert key_like(pl.outputs, "placement") in lo.inputs
    assert "genes" in set(lo.inputs) | set(lo.optional_inputs)


def test_miniprot_reads_the_long_proteins(steps: ModuleType) -> None:
    key = proteins_key(steps)
    assert Path(key).name.split(".")[0] in {"proteins_long", "proteins"}, key


def test_reference_peptides_are_an_optional_input(steps: ModuleType) -> None:
    """``--reference-peptides`` feeds the de novo precision (c_l); never required (§7.2.5)."""
    holders = [n for n in MAPPING_PY_STEPS if "reference_peptides" in steps.get(n).optional_inputs]
    assert holders, "no mapping step declares reference_peptides as an optional input"
    assert all("reference_peptides" not in steps.get(n).inputs for n in MAPPING_PY_STEPS)


def test_topological_order_follows_the_data_flow(steps: ModuleType) -> None:
    """§6.1: six-frame index → exact → tolerant (+ decoys) → placement FDR → loci."""
    order = steps.topo_order(MAPPING_PY_STEPS)
    pos = {n: i for i, n in enumerate(order)}
    chain = list(MAPPING_PY_STEPS)
    assert [pos[n] for n in chain] == sorted(pos[n] for n in chain), order


# ---------------------------------------------------------------- example 1, Q2, Q3


def test_example_1_only_confident_target_peptides_of_valid_length() -> None:
    df = norm_table(
        [
            {"peptide": "PEPTIDEK", "q_value": 0.001, "pep": 0.001},
            {"peptide": "KEDITPEP", "q_value": 0.001, "is_decoy": True, "protein_ids": ["rev_P1"]},
            {"peptide": "AAAK", "q_value": 0.001},
            {"peptide": "LONGPEPTIDEK", "q_value": 0.2},
        ]
    )
    out = select(df)
    assert out["peptide"].tolist() == ["PEPTIDEK"]


def test_q2_boundary_and_missing_q() -> None:
    cfg = GippaConfig()
    qmax = cfg.input.q_value_max
    df = norm_table(
        [
            {"peptide": "AAAAAAAK", "q_value": qmax},  # = cut-off: mapped
            {"peptide": "CCCCCCCK", "q_value": np.nextafter(qmax, 1.0)},  # just above: not
            {"peptide": "DDDDDDDK", "q_value": np.nan},  # no q: mapped
        ]
    )
    assert set(select(df, cfg)["peptide"]) == {"AAAAAAAK", "DDDDDDDK"}


def test_q2_uses_the_configured_cut_off() -> None:
    cfg = resolve_config(None, ["input.q_value_max=0.25"])
    df = norm_table([{"peptide": "LONGPEPTIDEK", "q_value": 0.2}])
    assert select(df, cfg)["peptide"].tolist() == ["LONGPEPTIDEK"]


def test_length_window_comes_from_mapping_exact() -> None:
    df = norm_table(
        [
            {"peptide": "A" * 6 + "K", "q_value": 0.001},  # 7: min_len
            {"peptide": "A" * 39 + "K", "q_value": 0.001},  # 40: max_len
            {"peptide": "A" * 40 + "K", "q_value": 0.001},  # 41
        ]
    )
    assert {len(p) for p in select(df)["peptide"]} == {7, 40}
    cfg = resolve_config(None, ["mapping.exact.min_len=8", "mapping.exact.max_len=30"])
    assert select(df, cfg).empty


def test_q3_auto_source_is_database_when_the_table_has_q_or_pep() -> None:
    df = norm_table(
        [
            {"peptide": "PEPTIDEK", "q_value": 0.001, "source": "auto"},
            {"peptide": "PEPTLDER", "pep": 0.01, "source": "auto"},
        ]
    )
    out = select(df)
    assert "source" in out.columns
    assert set(out["source"]) == {"database"}


def test_q3_auto_source_is_denovo_without_q_and_pep() -> None:
    df = norm_table(
        [
            {"peptide": "PEPTIDEK", "score": 0.9, "source": "auto"},
            {"peptide": "PEPTLDER", "score": 0.5, "source": "auto"},
        ]
    )
    out = select(df)
    assert set(out["source"]) == {"denovo"}


def test_explicit_sources_are_kept() -> None:
    df = norm_table(
        [
            {"peptide": "PEPTIDEK", "q_value": 0.001, "source": "database"},
            {"peptide": "PEPTLDER", "score": 0.5, "source": "denovo"},
        ]
    )
    out = select(df).set_index("peptide")
    assert out.loc["PEPTIDEK", "source"] == "database"
    assert out.loc["PEPTLDER", "source"] == "denovo"


_row = st.fixed_dictionaries(
    {
        "peptide": st.text(alphabet=AA, min_size=1, max_size=45),
        "q_value": st.one_of(st.just(np.nan), st.floats(0.0, 1.0)),
        "is_decoy": st.booleans(),
    }
)


@settings(max_examples=60, deadline=None)
@given(st.lists(_row, min_size=0, max_size=25, unique_by=lambda r: r["peptide"]))
def test_property_selection_is_exactly_the_query_rule(rows: list[dict[str, Any]]) -> None:
    cfg = GippaConfig()
    df = norm_table(rows) if rows else norm_table([{"peptide": "X"}]).iloc[0:0]
    out = select(df, cfg)
    assert set(out["peptide"]) == expected_queries(df, cfg)
    assert out["peptide"].is_unique


# ---------------------------------------------------------------- chain on a 50-kb simulation


class Chain:
    def __init__(
        self,
        run: SimRun,
        cfg: GippaConfig,
        avail: dict[str, Path],
        reg: ModuleType,
        records: list[logging.LogRecord],
    ):
        self.records = records
        self.run = run
        self.cfg = cfg
        self.avail = avail
        self.work = next(iter(avail.values())).parent
        self.keys = {
            "sixframe": key_like(reg.get("sixframe_index").outputs, "sixframe"),
            "hits_exact": key_like(reg.get("map_exact").outputs, "hits_exact"),
            "hits_tolerant": key_like(reg.get("map_tolerant").outputs, "hits_tolerant"),
            "placement": key_like(reg.get("placement").outputs, "placement"),
        }
        self.files = module2_files(self.work, avail)

    def path(self, stem: str) -> Path:
        return self.avail[self.keys[stem]]


def _with_near_misses(run: SimRun, out: Path) -> Path:
    """Add one-substitution variants of 10 correct ≥ 12-aa peptides (tolerant route)."""
    df = pd.read_parquet(run.peptides)
    keep = df[~df["is_decoy"] & (df["peptide"].str.len() >= 12)].head(10)
    rows = []
    for r in keep.to_dict("records"):
        p = str(r["peptide"])
        i = len(p) // 2
        sub = "W" if p[i] != "W" else "Y"
        rows.append({**r, "peptide": p[:i] + sub + p[i + 1 :]})
    new = pd.concat([df, pd.DataFrame(rows)], ignore_index=True)
    path = out / "peptides.near.parquet"
    new.to_parquet(path, index=False)
    return path


def _run_chain(tmp: Path, seed: int = SEED, denovo: bool = False) -> Chain:
    from gippasim.proteomics import ProteomicsConfig

    run = simulate_with_proteomics(
        tmp / "in", seed, proteomics=ProteomicsConfig(pi_wrong=0.05, search_db="other")
    )
    peps = _with_near_misses(run, tmp)
    if denovo:
        df = pd.read_parquet(peps)
        df["source"] = pd.Series(["denovo"] * len(df), dtype="string")
        df["q_value"] = np.nan
        df["pep"] = np.nan
        peps = tmp / "peptides.denovo.parquet"
        df.to_parquet(peps, index=False)
    run = SimRun(run.res, run.genome, peps, run.proteins_long, run.config_yaml)
    cfg = resolve_config(run.config_yaml, [f"seed={seed}"])
    records: list[logging.LogRecord] = []
    handler = _ListHandler(records)
    log = logging.getLogger("gippa.steps")
    old = log.level
    log.addHandler(handler)
    log.setLevel(logging.INFO)
    try:
        with pytest.MonkeyPatch.context() as mp:
            reg = fresh_steps(mp)
            avail = drive(reg, cfg, tmp / "work", roots_of(run))
            return Chain(run, cfg, avail, reg, records)
    finally:
        log.removeHandler(handler)
        log.setLevel(old)


class _ListHandler(logging.Handler):
    def __init__(self, sink: list[logging.LogRecord]) -> None:
        super().__init__(logging.INFO)
        self.sink = sink

    def emit(self, record: logging.LogRecord) -> None:
        self.sink.append(record)


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return _run_chain(tmp_path_factory.mktemp("chain"))


@pytest.fixture(scope="module")
def oracle(chain: Chain) -> dict[str, Any]:
    """The same chain through the pure M3 functions (M11-20: ``match_exact`` with the index's
    initiator start cells, as ``map_exact`` passes them)."""
    cfg = chain.cfg
    idx, decoys = load_index(chain.path("sixframe"), table=11)
    peps = pd.read_parquet(chain.run.peptides)
    queries = sorted(expected_queries(peps, cfg))
    decoy_of = dict.fromkeys(queries, "reverse")
    hx = match_exact(
        queries,
        idx,
        decoys,
        decoy_of,
        min_len=cfg.mapping.exact.min_len,
        max_len=cfg.mapping.exact.max_len,
        starts=load_start_cells(chain.path("sixframe"), table=11),
    )
    with_target = set(hx.loc[hx["search"] == "target", "peptide"])
    tol_q = [q for q in queries if len(q) >= cfg.mapping.tolerant.min_len]
    left = {q for q in tol_q if q not in with_target}
    tcfg = cfg.mapping.tolerant
    ht = pd.concat(
        [
            match_tolerant(tol_q, idx, decoys, decoy_of, left, tcfg),
            match_frameshift(tol_q, idx, decoys, decoy_of, left, tcfg),
        ],
        ignore_index=True,
    )
    qdf = peps[peps["peptide"].isin(queries)].copy()
    qdf["pi_denovo"] = np.nan
    costs = best_costs(pd.concat([hx, ht], ignore_index=True), qdf, residue_freqs(idx))
    per, fdr = placement(costs, cfg)
    return {
        "index": idx,
        "decoys": decoys,
        "queries": queries,
        "hx": hx,
        "ht": ht,
        "left": left,
        "per": per,
        "fdr": fdr,
    }


def _keyset(df: pd.DataFrame) -> set[tuple[Any, ...]]:
    d = df[HIT_KEY].copy()
    d["decoy_type"] = d["decoy_type"].where(d["search"] == "decoy", "")
    return {tuple(r) for r in d.astype(str).itertuples(index=False)}


def test_sixframe_work_file_is_the_index_and_its_decoys(chain: Chain) -> None:
    idx, decoys = load_index(chain.path("sixframe"), table=11)
    model = parse_bakta_json(chain.run.res.bakta_json).model
    reps = dict(model.replicons)
    assert {c: (r.length, r.topology) for c, r in idx.replicons.items()} == {
        c: (r.length, r.topology) for c, r in reps.items()
    }
    want = build_sixframe(_fasta(chain.run.genome), reps)
    assert len(idx.frames) == len(want.frames)
    for a, b in zip(idx.frames, want.frames, strict=True):
        assert (a.contig, a.strand, a.offset) == (b.contig, b.strand, b.offset)
        assert np.array_equal(a.aa, b.aa)
    assert "reverse" in decoys  # database input -> reverse decoys (mapping.decoys)
    rev = build_decoys(want, "reverse", 50, for_step("sixframe_index", chain.cfg.seed))
    for a, b in zip(decoys["reverse"], rev, strict=True):
        assert np.array_equal(a.aa, b.aa)


def _fasta(path: Path) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    cur = ""
    for ln in path.read_text().splitlines():
        if ln.startswith(">"):
            cur = ln[1:].split()[0]
            out[cur] = []
        elif ln:
            out[cur].append(ln.strip())
    return {k: "".join(v) for k, v in out.items()}


def test_map_exact_equals_match_exact_on_the_queries(chain: Chain, oracle: dict) -> None:
    got = pd.read_parquet(chain.path("hits_exact"))
    assert set(got["route"]) <= {"exact"}
    assert _keyset(got) == _keyset(oracle["hx"])
    # §7.2.1 decoys are counted, never placed; input decoy rows and q > cut-off never searched
    assert set(got["peptide"]) <= set(oracle["queries"])


def test_map_tolerant_searches_targets_only_for_exact_leftovers(chain: Chain, oracle: dict) -> None:
    got = pd.read_parquet(chain.path("hits_tolerant"))
    tgt = got[got["search"] == "target"]
    assert set(tgt["peptide"]) <= oracle["left"], "§7.2: only failures move on"
    want = oracle["ht"]
    want_t = want[want["search"] == "target"]
    assert set(tgt["peptide"]) == set(want_t["peptide"])
    for pep, m in want_t.groupby("peptide")["n_diff"].min().items():
        assert tgt.loc[tgt["peptide"] == pep, "n_diff"].min() == m, pep
    # the planted one-substitution variants take the tolerant route
    assert len(tgt) > 0


def test_map_tolerant_searches_decoys_for_every_query(chain: Chain, oracle: dict) -> None:
    """M3-10 Q1 (lead decision 1): tolerant decoy search for every query ≥ min_len."""
    got = pd.read_parquet(chain.path("hits_tolerant"))
    dec = got[got["search"] == "decoy"]
    want = oracle["ht"]
    want_d = want[want["search"] == "decoy"]
    assert set(dec["peptide"]) == set(want_d["peptide"])
    for pep, m in want_d.groupby("peptide")["n_diff"].min().items():
        assert dec.loc[dec["peptide"] == pep, "n_diff"].min() == m, pep


def test_placement_equals_the_pure_functions(chain: Chain, oracle: dict) -> None:
    got = pd.read_parquet(chain.path("placement"))
    assert {"peptide", "placement_q"} <= set(got.columns)
    want = oracle["per"].set_index("peptide")["placement_q"].astype(float)
    have = got.drop_duplicates("peptide").set_index("peptide")["placement_q"].astype(float)
    assert set(have.index) == set(want.index)
    np.testing.assert_allclose(have.loc[want.index].to_numpy(), want.to_numpy(), rtol=1e-12)


def test_mapping_fdr_has_the_placement_rows(chain: Chain, oracle: dict) -> None:
    fdr = read_tsv_df(chain.files["mapping_fdr.tsv"])
    pl = fdr[fdr["section"] == "placement"]
    assert set(pl["decoy_type"]) >= {"reverse"}
    got = pl[(pl["decoy_type"] == "reverse") & (pl["length_class"] == "all")]
    want = oracle["fdr"]
    want = want[(want["decoy_type"] == "reverse") & (want["length_class"] == "all")]
    assert sorted(zip(got["n_diff"].astype(int), got["targets"].astype(int), strict=True)) == (
        sorted(zip(want["n_diff"], want["targets"], strict=True))
    )


def test_seven_module_2_outputs_and_contract(chain: Chain, oracle: dict) -> None:
    files = chain.files
    hits = _read_tsv_table("peptide_hits", files["peptide_hits.tsv"])
    errors = [i for i in validate_table("peptide_hits", hits) if i.severity == "error"]
    assert not errors, errors
    prots = _read_tsv_table("protein_hits", files["protein_hits.tsv"])
    assert len(prots) == 0  # peptides only
    un = _read_tsv_table("unmapped", files["unmapped.tsv"])
    assert not [i for i in validate_table("unmapped", un) if i.severity == "error"]
    gm = parse_bakta_json(chain.run.res.bakta_json).model
    tags = {g.locus_tag for g in gm.genes}
    assert set(hits["gene"].dropna()) <= tags
    for f in ("mappings.bed", "mappings.probed", "mappings.gff3"):
        assert files[f].is_file()


def test_only_queries_reach_the_hit_and_unmapped_tables(chain: Chain, oracle: dict) -> None:
    """§7.3.1 / Q2: decoy rows and q > cut-off are never mapped, never in the unmapped bucket."""
    peps = pd.read_parquet(chain.run.peptides)
    hits = read_tsv_df(chain.files["peptide_hits.tsv"])
    un = read_tsv_df(chain.files["unmapped.tsv"])
    un_p = set(un.loc[un["kind"] == "peptide", "sequence"])
    queries = set(oracle["queries"])
    decoys = set(peps.loc[peps["is_decoy"], "peptide"])
    assert not set(hits["peptide"]) & decoys
    assert set(hits["peptide"]) <= queries
    assert un_p <= queries
    assert not set(hits["peptide"]) & un_p
    assert queries <= set(hits["peptide"]) | un_p
    q = pd.to_numeric(hits["best_q_input"], errors="coerce")
    assert (q.isna() | (q <= chain.cfg.input.q_value_max)).all()


def test_every_mapping_step_logs_counts(chain: Chain) -> None:
    for name in MAPPING_PY_STEPS:
        recs = [r for r in chain.records if r.name == f"gippa.steps.{name}"]
        assert any(r.levelno >= logging.INFO and re.search(r"\d", r.getMessage()) for r in recs), (
            f"{name}: no count logged"
        )


def test_rerun_is_byte_identical(chain: Chain, tmp_path: Path) -> None:
    """§6.6: same inputs and seed → identical mapping TSVs and work tables."""
    again = _run_chain(tmp_path)
    for f in MAPPING_TSVS:
        assert again.files[f].read_bytes() == chain.files[f].read_bytes(), f
    for stem in ("hits_exact", "hits_tolerant", "placement"):
        a = pd.read_parquet(again.path(stem))
        b = pd.read_parquet(chain.path(stem))
        pd.testing.assert_frame_equal(a, b)


def test_tables_carry_no_absolute_paths(chain: Chain) -> None:
    for f in MAPPING_TSVS:
        text = chain.files[f].read_text()
        assert str(chain.work) not in text and str(Path.home()) not in text, f


# ---------------------------------------------------------------- example 6 (de novo, no ref.)


@pytest.fixture(scope="module")
def denovo_chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return _run_chain(tmp_path_factory.mktemp("denovo"), denovo=True)


def test_example_6_denovo_without_reference_assumes_c(denovo_chain: Chain) -> None:
    fdr = read_tsv_df(denovo_chain.files["mapping_fdr.tsv"])
    dn = fdr[fdr["section"] == "denovo_precision"]
    assert len(dn) > 0, "no de novo precision rows in mapping_fdr.tsv"
    assert dn["flags"].fillna("").str.split(";").map(lambda f: "c_assumed" in f).all()
    pl = fdr[fdr["section"] == "placement"]
    assert "shuffle" in set(pl["decoy_type"])  # mapping.decoys.denovo_input


def test_example_6_c_assumed_is_logged(denovo_chain: Chain) -> None:
    assert any(
        r.levelno >= logging.WARNING and "c_assumed" in r.getMessage()
        for r in denovo_chain.records
        if r.name.startswith("gippa.steps.")
    )


def test_denovo_decoys_are_the_step_stream_shuffle(denovo_chain: Chain) -> None:
    idx, decoys = load_index(denovo_chain.path("sixframe"), table=11)
    assert "shuffle" in decoys
    cfg = denovo_chain.cfg
    want = build_decoys(
        idx, "shuffle", cfg.mapping.decoys.shuffle_window_aa, for_step("sixframe_index", cfg.seed)
    )
    for a, b in zip(decoys["shuffle"], want, strict=True):
        assert np.array_equal(a.aa, b.aa)


def test_denovo_rerun_is_byte_identical(denovo_chain: Chain, tmp_path: Path) -> None:
    again = _run_chain(tmp_path, denovo=True)
    for f in MAPPING_TSVS:
        assert again.files[f].read_bytes() == denovo_chain.files[f].read_bytes(), f


# ---------------------------------------------------------------- protein route (injection)


def _gff_without(src: Path, drop: set[str], dst: Path) -> Path:
    """The recorded miniprot GFF minus every alignment block (##PAF + features) of ``drop``."""
    out: list[str] = []
    keep = True
    for ln in src.read_text().splitlines():
        if ln.startswith("##PAF"):
            keep = ln.split("\t")[1] not in drop
        elif ln.startswith("##"):
            keep = True
        if keep:
            out.append(ln)
    dst.write_text("\n".join(out) + "\n")
    return dst


def _tsv_only(src: Path, ids: set[str], dst: Path) -> Path:
    rows = [ln for ln in src.read_text().splitlines() if ln.split("\t")[0] in ids]
    dst.write_text("\n".join(rows) + "\n")
    return dst


@pytest.fixture(scope="module")
def protein_chain(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    import dataclasses
    import hashlib

    from helpers.mapping_steps import preflight

    from gippasim.proteomics import ProteomicsConfig
    from gippasim.simulate import SIM100K, SIM100K_SEED, simulate

    tmp = tmp_path_factory.mktemp("proteins")
    cfg_sim = dataclasses.replace(SIM100K, proteomics=ProteomicsConfig())
    res = simulate(cfg_sim, SIM100K_SEED, tmp / "sim")
    sha = hashlib.sha256(res.genome_fasta.read_bytes()).hexdigest()
    assert sha.startswith("dd10216e"), "SIM100K genome differs from the recorded fixtures"
    run = preflight(res, tmp / "pf")
    ids = [ln[1:].split()[0] for ln in (RECORDED_MINIPROT / "proteins.faa").read_text().split("\n")
           if ln.startswith(">")]  # fmt: skip
    drop = set(ids[::2])
    gff = _gff_without(RECORDED_MINIPROT / "miniprot.gff", drop, tmp / "miniprot.part.gff")
    tbl = _tsv_only(RECORDED_TBLASTN / "tblastn.tsv", drop, tmp / "tblastn.part.tsv")
    cfg = resolve_config(run.config_yaml, ["seed=3"])
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        roots = roots_of(run) | {proteins_key(reg): RECORDED_MINIPROT / "proteins.faa"}
        (mp_key,) = reg.get("miniprot").outputs
        (tb_key,) = reg.get("tblastn").outputs
        tb_in = [k for k in reg.get("tblastn").inputs if k != "genome"]
        avail = drive(
            reg,
            cfg,
            tmp / "work",
            roots,
            inject={"miniprot": {mp_key: gff}, "tblastn": {tb_key: tbl}},
        )
    return {
        "ids": ids,
        "drop": drop,
        "files": module2_files(tmp / "work", avail),
        "leftover": [avail[k] for k in tb_in if k in avail],
    }


def test_proteins_take_the_cheapest_route(protein_chain: dict[str, Any]) -> None:
    """§7.2: miniprot first; only proteins without an assignment go to tblastn."""
    hits = read_tsv_df(protein_chain["files"]["protein_hits.tsv"])
    drop = protein_chain["drop"]
    routes = hits.groupby("protein_id")["route"].agg(set)
    assert set(routes.index) == set(protein_chain["ids"])
    for pid, r in routes.items():
        assert r == ({"tblastn"} if pid in drop else {"miniprot"}), (pid, r)
    first = hits[hits["rank"].astype(int) == 1].set_index("protein_id")
    assert (first["gene"] == first.index).all()  # IDs are the locus tags (PROVENANCE.md)
    assert set(first.loc[sorted(drop), "class"]) == {"divergent"}
    assert set(first.loc[sorted(set(protein_chain["ids"]) - drop), "class"]) == {"identical"}


def test_tblastn_gets_exactly_the_miniprot_leftovers(protein_chain: dict[str, Any]) -> None:
    (left,) = protein_chain["leftover"]
    ids = {ln[1:].split()[0] for ln in left.read_text().splitlines() if ln.startswith(">")}
    assert ids == protein_chain["drop"]


def test_assigned_proteins_are_not_unmapped(protein_chain: dict[str, Any]) -> None:
    un = read_tsv_df(protein_chain["files"]["unmapped.tsv"])
    assert un[un["kind"] == "protein"].empty


# ---------------------------------------------------------------- error exits


def test_peptides_without_a_peptide_column_fail_the_step(
    steps: ModuleType, chain: Chain, tmp_path: Path
) -> None:
    bad = tmp_path / "bad.parquet"
    pd.DataFrame({"sequence": ["PEPTIDEK"], "q_value": [0.001]}).to_parquet(bad)
    spec = steps.get("map_exact")
    ins = {k: chain.avail[k] for k in spec.inputs} | {"peptides": bad}
    outs = {k: tmp_path / k for k in spec.outputs}
    with pytest.raises(steps.StepError):
        steps.run_step("map_exact", chain.cfg, ins, outs)
    assert not any(p.exists() for p in outs.values())


def test_gippa_step_cli_exits_4_on_a_corrupt_index(
    steps: ModuleType, chain: Chain, tmp_path: Path
) -> None:
    """``steps`` = ``fresh_steps``: the real steps whatever registry an earlier test left behind
    (M3-17 review decision 12)."""
    from click.testing import CliRunner

    from gippa.cli import cli

    spec = steps.get("map_exact")
    six = key_like(spec.inputs, "sixframe")
    bad = tmp_path / "sixframe.npz"
    bad.write_bytes(b"not an npz")
    args = ["_step", "map_exact"]
    for k in spec.inputs:
        args += ["--in", f"{k}={bad if k == six else chain.avail[k]}"]
    for k in spec.outputs:
        args += ["--out", f"{k}={tmp_path / 'out' / k}"]
    r = CliRunner().invoke(cli, args)
    assert r.exit_code == 4, r.output
