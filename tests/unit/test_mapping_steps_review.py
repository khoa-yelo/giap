"""M3-17 lead decisions after code review (PR #77), items 6-12 of
docs/specs/M3-17-mapping-steps.md ("Lead decisions after code review"); item 10 (``_engine_ids``)
is pinned in tests/unit/test_loci_engine_ids.py.

Conventions pinned here where the decisions are silent:

* Item 6. ``giap.steps.mapping.resolve_search_db(peptides, database_ids) -> "same" | "other"`` is
  the table-level resolution of ``search_db: auto``: ``same`` when at least 95 % of the target
  (``is_decoy`` false) peptides have their ``protein_ids`` (a list, contract C9) among
  ``database_ids`` (the gene model's locus tags), else ``other``. Engine decoy rows never count.
  (The tests give every peptide one ID, so "any ID" and "all IDs" agree.) Rows with an explicit
  ``same`` / ``other`` keep it. "Recorded in the run metadata" = a mapping step logs a record on
  its ``giap.steps.<name>`` logger naming ``search_db`` and the resolved value (the CLI's sink
  writes it to ``logs/<step>.log`` and ``events.jsonl``).
* Item 7. ``giap run --executor local`` (documented codes, §5.1 / §6.3 / §6.4): no evidence → exit 1
  with the preflight message ``no evidence``; proteins only (miniprot cannot run locally) → exit 4
  naming ``miniprot`` (example 3); a GFF3 ``--annotation`` without
  ``annotation.precomputed.bakta_json`` → exit 1 naming ``bakta_json`` (§6.4 "a missing layer is an
  error naming the key"). Never exit 0 without the Module 2 outputs.
* Item 8. A leftover ``OUT/work/<key>`` from an earlier run (``proteins_long``, ``contaminants``)
  never enters a rerun: the rerun's mapping TSVs equal the first run's.
* Item 9. The contaminant library is ``$GIAP_DB_DIR/contaminants/contaminants.fasta`` with a
  ``contaminants`` manifest entry (M3-15); a confident query peptide that lands nowhere and lies in
  a library entry is labelled ``contaminant`` in ``unmapped.tsv`` with that entry in
  ``contaminant_match``.
* Item 11. The seven §7.2.9 files are in ``OUT/mapping`` (§9.1).
* Item 12. One test per mutant the reviewer saw survive (local/scratch/review-m3-17/mutants.py:
  M05 tolerant decoys per query, M06 frameshift wiring, M07 best-cost filter, M08 de novo π̂,
  M14 tolerant ``min_len`` boundary, M15 reference peptides, M16 leftover assignment, M17 tblastn
  leftover filter, M18 ``protein_ids`` separator, M23 genome SHA-256). M05, M06 and M14 use a
  hand-built genome; M07 calls ``giap.steps.mapping._peptide_rows`` (rows, placement, queries,
  config) directly.
* Re-review decisions 13-15. 13: probe E (``--peptides``, no ``--annotation``, no
  ``annotation.precomputed.bakta_json``) → exit 1 naming ``bakta_json``; missing optional roots
  still exit 0 without naming them. 14: ``OUT/logs/events.jsonl`` holds an event whose name or
  ``data`` keys contain ``search_db`` and whose ``data`` values include the resolved value. 15:
  N10 (``loci`` uses the resolved queries) and M10 (contaminant matching of every non-accepted
  query, incl. ``rejected_placement_fdr``) of local/scratch/review-m3-17/r2/mutants2.py.

Every preflight here uses a temporary ``$GIAP_DB_DIR`` with a test manifest
(``helpers.mapping_steps.write_test_db``).
"""

from __future__ import annotations

import dataclasses
import hashlib
import logging
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.mapping_steps import (
    CONTAMINANTS_MINI,
    MAPPING_TSVS,
    MODULE2_FILES,
    RECORDED_MINIPROT,
    RECORDED_TBLASTN,
    drive,
    fresh_steps,
    mapping_module,
    module2_files,
    preflight,
    proteins_key,
    read_tsv_df,
    roots_of,
    write_test_db,
)
from test_mapping_steps import norm_table  # type: ignore[import-not-found]

from giap.config import GiapConfig, resolve_config
from giap.mapping.chance import residue_freqs
from giap.mapping.decoys import build_decoys
from giap.mapping.hits import hits_frame
from giap.mapping.placement_fdr import best_costs, placement
from giap.mapping.sixframe import build_sixframe, load_index, save_index
from giap.model.coords import Replicon
from giap.rng import for_step

SEED = 21


@pytest.fixture
def steps() -> Iterator[ModuleType]:
    with pytest.MonkeyPatch.context() as mp:
        yield fresh_steps(mp)


class _ListHandler(logging.Handler):
    def __init__(self, sink: list[logging.LogRecord]) -> None:
        super().__init__(logging.INFO)
        self.sink = sink

    def emit(self, record: logging.LogRecord) -> None:
        self.sink.append(record)


def _drive_logged(
    cfg: GiapConfig, work: Path, roots: dict[str, Path], **kw: Any
) -> tuple[dict[str, Path], list[logging.LogRecord]]:
    records: list[logging.LogRecord] = []
    handler = _ListHandler(records)
    log = logging.getLogger("giap.steps")
    old = log.level
    log.addHandler(handler)
    log.setLevel(logging.INFO)
    try:
        with pytest.MonkeyPatch.context() as mp:
            reg = fresh_steps(mp)
            return drive(reg, cfg, work, roots, **kw), records
    finally:
        log.removeHandler(handler)
        log.setLevel(old)


def _flags(x: object) -> set[str]:
    if x is None or (not isinstance(x, str) and pd.isna(x)):
        return set()
    return {f for f in str(x).split(";") if f}


# ================================================================ item 6: search_db auto (pure)


def _db_table(n_target: int, n_in_db: int, n_decoy: int = 0) -> tuple[pd.DataFrame, set[str]]:
    """``n_target`` target peptides, the first ``n_in_db`` with a gene-model ID, the rest with a
    foreign ID; ``n_decoy`` engine decoys (``rev_`` IDs, never in the model)."""
    rng = np.random.default_rng(n_target * 1000 + n_in_db)
    aa = list("ACDEFGHKMNPQRSTVWY")
    seen: set[str] = set()
    rows: list[dict[str, Any]] = []
    while len(seen) < n_target + n_decoy:
        p = "".join(rng.choice(aa, 11)) + "K"
        if p in seen:
            continue
        seen.add(p)
        i = len(seen) - 1
        if i < n_target:
            pid = f"GSIM_{i:05d}" if i < n_in_db else f"OTHER_{i:05d}"
            rows.append({"peptide": p, "q_value": 0.001, "protein_ids": [pid], "search_db": "auto"})
        else:
            rows.append({"peptide": p, "q_value": 0.001, "is_decoy": True,
                         "protein_ids": [f"rev_GSIM_{i:05d}"], "search_db": "auto"})  # fmt: skip
    model = {f"GSIM_{i:05d}" for i in range(n_target + n_decoy)}
    return norm_table(rows), model


@pytest.mark.parametrize(
    ("n_target", "n_in_db", "want"),
    [
        (100, 100, "same"),
        (100, 95, "same"),  # exactly 95 %
        (100, 94, "other"),  # just below
        (20, 19, "same"),
        (20, 18, "other"),
        (100, 0, "other"),
    ],
)
def test_item6_search_db_auto_95_percent_rule(n_target: int, n_in_db: int, want: str) -> None:
    df, model = _db_table(n_target, n_in_db)
    got = mapping_module().resolve_search_db(df, model)
    assert got == want


def test_item6_engine_decoys_do_not_count() -> None:
    """95/100 targets in the model plus 10 decoys (``rev_`` IDs): still ``same``."""
    df, model = _db_table(100, 95, n_decoy=10)
    assert mapping_module().resolve_search_db(df, model) == "same"


def test_item6_protein_ids_may_be_arrays() -> None:
    """Parquet gives ``protein_ids`` back as numpy arrays (C9 list column)."""
    df, model = _db_table(20, 19)
    df["protein_ids"] = [np.asarray(v, dtype=object) for v in df["protein_ids"]]
    assert mapping_module().resolve_search_db(df, model) == "same"


# ================================================================ chains on a same-strain sim


def _sim(tmp: Path, seed: int, **prot: Any) -> Any:
    from giapsim.proteomics import ProteomicsConfig
    from giapsim.simulate import SimConfig, simulate

    cfg = SimConfig(proteomics=ProteomicsConfig(**prot), write_precomputed=True)
    return simulate(cfg, seed, tmp / "sim")


def _with_parquet(run: Any, df: pd.DataFrame, path: Path) -> Any:
    df.to_parquet(path, index=False)
    return dataclasses.replace(run, peptides=path)


EXPLICIT_SAME = 3


@pytest.fixture(scope="module")
def same_chain(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """50-kb sim searched against its own proteins; ``search_db`` set to ``auto`` on every row
    except three target peptides that say ``same`` explicitly and carry an extra foreign
    ``protein_ids`` entry (``XOTHER``)."""
    tmp = tmp_path_factory.mktemp("same")
    run = preflight(_sim(tmp, SEED, search_db="same"), tmp / "pf")
    df = pd.read_parquet(run.peptides)
    df["search_db"] = pd.Series(["auto"] * len(df), dtype="string")
    tgt = df.index[~df["is_decoy"].astype(bool) & (df["peptide"].str.len() >= 8)][:EXPLICIT_SAME]
    df.loc[tgt, "search_db"] = "same"
    df["protein_ids"] = [
        [*list(v), "XOTHER"] if i in set(tgt) else list(v)
        for i, v in zip(df.index, df["protein_ids"], strict=True)
    ]
    run = _with_parquet(run, df, tmp / "peptides.auto.parquet")
    cfg = resolve_config(run.config_yaml, [f"seed={SEED}"])
    avail, records = _drive_logged(cfg, tmp / "work", roots_of(run))
    return {
        "run": run,
        "cfg": cfg,
        "avail": avail,
        "records": records,
        "files": module2_files(tmp / "work", avail),
        "explicit": set(df.loc[tgt, "peptide"]),
        "ids": dict(zip(df["peptide"], df["protein_ids"], strict=True)),
    }


def _auto_exact(chain: dict[str, Any]) -> pd.DataFrame:
    per = pd.read_parquet(chain["avail"]["placement"])
    per = per[(per["source"] == "database") & (per["m_target"].astype(int) == 0)]
    return per[~per["peptide"].isin(chain["explicit"])]


def test_item6_auto_resolves_to_same_on_a_same_strain_search(same_chain: dict[str, Any]) -> None:
    """≥ 95 % of the IDs are locus tags: exact database hits take the same-strain shortcut."""
    per = _auto_exact(same_chain)
    assert len(per) > 100
    assert per["flags"].map(lambda f: "same_db" in _flags(f)).all()
    assert (per["placement_q"].astype(float) == 0).all()


def test_item6_resolution_is_recorded(same_chain: dict[str, Any]) -> None:
    recs = [
        r.getMessage()
        for r in same_chain["records"]
        if r.name.startswith("giap.steps.") and r.levelno >= logging.INFO
    ]
    assert any("search_db" in m and "same" in m for m in recs), recs


@pytest.fixture(scope="module")
def other_chain(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """The same-strain sim with 10 % of the target peptides' IDs replaced: ``auto`` → ``other``."""
    tmp = tmp_path_factory.mktemp("other")
    run = preflight(_sim(tmp, SEED, search_db="same"), tmp / "pf")
    df = pd.read_parquet(run.peptides)
    df["search_db"] = pd.Series(["auto"] * len(df), dtype="string")
    tgt = df.index[~df["is_decoy"].astype(bool)]
    foreign = set(tgt[::10])
    df["protein_ids"] = [
        [f"FOREIGN_{i}"] if i in foreign else list(v)
        for i, v in zip(df.index, df["protein_ids"], strict=True)
    ]
    run = _with_parquet(run, df, tmp / "peptides.other.parquet")
    cfg = resolve_config(run.config_yaml, [f"seed={SEED}"])
    avail, records = _drive_logged(cfg, tmp / "work", roots_of(run), targets=("placement",))
    return {"avail": avail, "records": records, "explicit": set()}


def test_item6_auto_resolves_to_other_below_95_percent(other_chain: dict[str, Any]) -> None:
    per = _auto_exact(other_chain)
    assert len(per) > 100
    assert not per["flags"].map(lambda f: "same_db" in _flags(f)).any()
    msgs = [r.getMessage() for r in other_chain["records"] if r.name.startswith("giap.steps.")]
    assert any("search_db" in m and "other" in m for m in msgs), msgs


# ---------------------------------------------------------------- M18 / item 10, M23


def test_m18_engine_ids_list_with_the_hit_gene_is_no_mismatch(same_chain: dict[str, Any]) -> None:
    """A ``same`` peptide whose ``protein_ids`` list holds the hit gene (plus another ID) is not
    ``engine_protein_mismatch``; the list must never be joined with a wrong separator."""
    hits = read_tsv_df(same_chain["files"]["peptide_hits.tsv"])
    rows = hits[hits["peptide"].isin(same_chain["explicit"]) & hits["gene"].notna()]
    rows = rows[[g in set(same_chain["ids"][p]) for p, g in zip(rows["peptide"], rows["gene"],
                                                                strict=True)]]  # fmt: skip
    assert len(rows) > 0, "no explicit same-strain peptide hit its own gene"
    assert not rows["flags"].map(lambda f: "engine_protein_mismatch" in _flags(f)).any()


def test_item10_ids_join_workaround_is_removed() -> None:
    assert not hasattr(mapping_module(), "_ids")


def test_m23_probed_records_carry_the_genome_sha256(same_chain: dict[str, Any]) -> None:
    """proBed ``genomeReferenceVersion`` = ``sha256:<hex>`` of ``genome.norm.fa``, every record."""
    from giap.io.probed import PROBED_EXTRA

    sha = hashlib.sha256(same_chain["run"].genome.read_bytes()).hexdigest()
    col = 12 + PROBED_EXTRA.index("genomeReferenceVersion")
    recs = [
        ln.split("\t") for ln in same_chain["files"]["mappings.probed"].read_text().splitlines()
        if ln and not ln.startswith(("#", "track", "browser"))
    ]  # fmt: skip
    assert len(recs) > 0
    assert {r[col] for r in recs} == {f"sha256:{sha}"}


# ================================================================ M05, M06, M14: built genome

CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
    "I": "ATT", "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCG", "Q": "CAG",
    "R": "CGT", "S": "TCT", "T": "ACC", "V": "GTT", "W": "TGG", "Y": "TAT",
}  # fmt: skip
P_EXACT = "MPYTDWHKFERS"  # exact target hit; 1-substitution reverse decoy hit (M05)
P_EXACT_DECOY = P_EXACT[:6] + "C" + P_EXACT[7:]  # placed reversed in the target frame
P_TEN = "YWMCHPFDTR"  # 10 aa = tolerant.min_len, one substitution from the genome (M14)
P_TEN_GENOME = P_TEN[:4] + "E" + P_TEN[5:]
FS_A, FS_B = "FWYCMH", "PDTRWE"  # A | +1 nt | B: a frameshift-split peptide (M06)
P_FS = FS_A + FS_B


def _codons(aa: str) -> str:
    return "".join(CODON[a] for a in aa)


@pytest.fixture(scope="module")
def built(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """A 3-kb linear contig with the planted segments in frame +0, separated by random spacers;
    map_exact and map_tolerant run through ``run_step`` on its index (reverse decoys)."""
    tmp = tmp_path_factory.mktemp("built")
    rng = np.random.default_rng(5)

    def spacer(n: int = 600) -> str:
        return "".join(rng.choice(list("ACGT"), n))

    parts = [
        spacer(), _codons(P_EXACT), spacer(), _codons(P_EXACT_DECOY[::-1]), spacer(),
        _codons(P_TEN_GENOME), spacer(), _codons(FS_A) + "G" + _codons(FS_B), spacer(),
    ]  # fmt: skip
    seq = "".join(parts)
    seq += "A" * (-len(seq) % 3)
    reps = {"chr1": Replicon("chr1", len(seq), "linear", "chromosome")}
    index = build_sixframe({"chr1": seq}, reps)
    cfg = GiapConfig(seed=3)
    decoys = {"reverse": build_decoys(index, "reverse", 50, for_step("sixframe_index", 3))}
    six = tmp / "sixframe"
    save_index(six, index, decoys)
    peps = tmp / "peptides"
    rows = [{"peptide": p, "q_value": 0.001, "pep": 0.001} for p in (P_EXACT, P_TEN, P_FS)]
    norm_table(rows).to_parquet(peps, index=False)
    avail = {"peptides": peps, "sixframe": six}
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        for name in ("map_exact", "map_tolerant"):
            spec = reg.get(name)
            outs = {k: tmp / k for k in spec.outputs}
            reg.run_step(name, cfg, {k: avail[k] for k in spec.inputs}, outs)
            avail.update(outs)
    return {
        "exact": pd.read_parquet(avail["hits_exact"]),
        "tol": pd.read_parquet(avail["hits_tolerant"]),
    }


def _is_fs(df: pd.DataFrame) -> pd.Series:
    return df["diff_types"].fillna("").str.split(";").map(lambda t: "frameshift" in t)


def test_built_genome_sanity(built: dict[str, Any]) -> None:
    ex = built["exact"]
    assert set(ex.loc[ex["search"] == "target", "peptide"]) == {P_EXACT}


def test_m05_tolerant_decoys_are_searched_for_queries_with_an_exact_hit(
    built: dict[str, Any],
) -> None:
    """M3-10 Q1: the tolerant decoy search covers every query ≥ min_len, also exact-hit ones."""
    tol = built["tol"]
    d = tol[(tol["peptide"] == P_EXACT) & (tol["search"] == "decoy") & ~_is_fs(tol)]
    assert len(d) >= 1, tol
    assert set(d["decoy_type"]) == {"reverse"}
    assert d["n_diff"].min() == 1
    # its target side is not searched again (§7.2: only failures move on)
    assert tol[(tol["peptide"] == P_EXACT) & (tol["search"] == "target")].empty


def test_m06_frameshift_split_is_wired(built: dict[str, Any]) -> None:
    tol = built["tol"]
    t = tol[(tol["peptide"] == P_FS) & (tol["search"] == "target") & _is_fs(tol)]
    assert len(t) >= 1, tol
    assert t["blocks"].notna().all()


def test_m14_tolerant_min_len_is_inclusive(built: dict[str, Any]) -> None:
    tol = built["tol"]
    t = tol[(tol["peptide"] == P_TEN) & (tol["search"] == "target") & ~_is_fs(tol)]
    assert len(t) == 1, tol
    assert int(t["n_diff"].iloc[0]) == 1


# ================================================================ M07: best-cost rows only


def test_m07_peptide_rows_keep_only_the_best_cost_hits() -> None:
    pep = "PEPTLDEKAAWK"

    def hit(start: int, n_diff: int, search: str = "target") -> dict[str, Any]:
        return {
            "peptide": pep, "query": pep, "search": search,
            "decoy_type": "reverse" if search == "decoy" else None, "contig": "chr1",
            "frame_id": 0, "aa_start": start // 3, "aa_end": start // 3 + 12, "start": start,
            "end": start + 36, "strand": "+", "frame": 1, "route": "tolerant", "n_diff": n_diff,
            "diff_types": "substitution" * (n_diff > 0), "blocks": None,
        }  # fmt: skip

    hits = hits_frame([hit(300, 1), hit(900, 2), hit(1500, 1), hit(2100, 1, "decoy")])
    per = pd.DataFrame(
        {
            "peptide": [pep], "source": ["database"], "decoy_type": ["reverse"], "length": [12],
            "m_target": [1], "e_target": [0.01], "placement_q": [0.001], "accepted": [True],
            "flags": [""],
        }
    )  # fmt: skip
    queries = norm_table([{"peptide": pep, "q_value": 0.002, "pep": 0.003}])
    rows = mapping_module()._peptide_rows(hits, per, queries, GiapConfig())
    assert sorted(rows["start"].astype(int)) == [300, 1500]
    assert set(rows["n_diff"].astype(int)) == {1}
    assert set(rows["search"]) == {"target"}


# ================================================================ M08, M15: scored de novo input


@pytest.fixture(scope="module")
def denovo_scored(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """giapsim de novo reads (with scores) as the peptide input; placement with and without the
    simulator's ``reference_peptides.tsv``."""
    from giapsim.denovo import DenovoConfig
    from giapsim.proteomics import ProteomicsConfig
    from giapsim.simulate import SimConfig, simulate

    tmp = tmp_path_factory.mktemp("denovo_scored")
    sim = simulate(
        SimConfig(proteomics=ProteomicsConfig(), denovo=DenovoConfig(), write_precomputed=True),
        SEED,
        tmp / "sim",
    )
    res = dataclasses.replace(sim, peptides=sim.outdir / "peptides_denovo.tsv")
    run = preflight(res, tmp / "pf")
    cfg = resolve_config(run.config_yaml, [f"seed={SEED}"])
    roots = roots_of(run)
    no_ref, _ = _drive_logged(cfg, tmp / "noref", roots, targets=("placement",))
    ref = roots | {"reference_peptides": sim.outdir / "reference_peptides.tsv"}
    with_ref, _ = _drive_logged(cfg, tmp / "ref", ref, targets=("placement",))
    return {"cfg": cfg, "run": run, "no_ref": no_ref, "ref": with_ref}


def _denovo_rows(avail: dict[str, Path]) -> pd.DataFrame:
    fdr = pd.read_parquet(avail["mapping_fdr"])
    return fdr[fdr["section"] == "denovo_precision"]


def test_m15_reference_peptides_replace_the_assumed_c(denovo_scored: dict[str, Any]) -> None:
    without = _denovo_rows(denovo_scored["no_ref"])
    assert len(without) > 0
    assert without["flags"].map(lambda f: "c_assumed" in _flags(f)).all()
    with_ref = _denovo_rows(denovo_scored["ref"])
    assert len(with_ref) > 0
    assert not with_ref["flags"].map(lambda f: "c_assumed" in _flags(f)).any()


def test_m08_denovo_pi_lowers_the_placement_q(denovo_scored: dict[str, Any]) -> None:
    """π̂_wrong from the de novo precision enters placement: compared with π = 1 (the pure
    functions with ``pi_denovo`` missing) no de novo q rises and some fall (§7.2.5, ADR-0010)."""
    cfg = denovo_scored["cfg"]
    avail = denovo_scored["no_ref"]
    index, _ = load_index(avail["sixframe"])
    queries = mapping_module().select_queries(pd.read_parquet(avail["peptides"]), cfg)
    queries["pi_denovo"] = np.nan
    hits = pd.concat(
        [pd.read_parquet(avail["hits_exact"]), pd.read_parquet(avail["hits_tolerant"])],
        ignore_index=True,
    )
    flat, _ = placement(best_costs(hits, queries, residue_freqs(index)), cfg)
    got = pd.read_parquet(avail["placement"])
    a = got[got["source"] == "denovo"].set_index("peptide")["placement_q"].astype(float)
    b = flat[flat["source"] == "denovo"].set_index("peptide")["placement_q"].astype(float)
    assert set(a.index) == set(b.index) and len(a) > 0
    b = b.loc[a.index]
    ok = a.notna() & b.notna()
    assert (a[ok] <= b[ok] + 1e-12).all()
    assert (a[ok] < b[ok] - 1e-9).any(), "π̂ never lowered a de novo placement q"


# ================================================================ M16, M17: protein route


@pytest.fixture(scope="module")
def protein_full(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """SIM100K with the recorded miniprot GFF and the *full* recorded tblastn table injected
    (tblastn hits exist for proteins miniprot already assigned)."""
    from giapsim.proteomics import ProteomicsConfig
    from giapsim.simulate import SIM100K, SIM100K_SEED, simulate

    tmp = tmp_path_factory.mktemp("protein_full")
    cfg_sim = dataclasses.replace(SIM100K, proteomics=ProteomicsConfig())
    res = simulate(cfg_sim, SIM100K_SEED, tmp / "sim")
    sha = hashlib.sha256(res.genome_fasta.read_bytes()).hexdigest()
    assert sha.startswith("dd10216e"), "SIM100K genome differs from the recorded fixtures"
    run = preflight(res, tmp / "pf")
    cfg = resolve_config(run.config_yaml, ["seed=3"])
    faa = RECORDED_MINIPROT / "proteins.faa"
    ids = [ln[1:].split()[0] for ln in faa.read_text().splitlines() if ln.startswith(">")]
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        roots = roots_of(run) | {proteins_key(reg): faa}
        (mp_key,) = reg.get("miniprot").outputs
        (tb_key,) = reg.get("tblastn").outputs
        inject = {
            "miniprot": {mp_key: RECORDED_MINIPROT / "miniprot.gff"},
            "tblastn": {tb_key: RECORDED_TBLASTN / "tblastn.tsv"},
        }
        avail = drive(reg, cfg, tmp / "work", roots, inject=inject)
    return {"cfg": cfg, "ids": ids, "avail": avail, "tmp": tmp,
            "files": module2_files(tmp / "work", avail)}  # fmt: skip


def test_m17_tblastn_hits_of_miniprot_assigned_proteins_are_ignored(
    protein_full: dict[str, Any],
) -> None:
    tsv = (RECORDED_TBLASTN / "tblastn.tsv").read_text().splitlines()
    in_tsv = {ln.split("\t")[0] for ln in tsv if ln and not ln.startswith("#")}
    assert in_tsv & set(protein_full["ids"]), "recorded tblastn has no assigned protein"
    hits = read_tsv_df(protein_full["files"]["protein_hits.tsv"])
    routes = hits.groupby("protein_id")["route"].agg(set)
    assert set(routes.index) == set(protein_full["ids"])
    assert all(r == {"miniprot"} for r in routes), routes[routes.map(lambda r: r != {"miniprot"})]


def test_m16_aligned_but_unassigned_proteins_go_to_tblastn(
    protein_full: dict[str, Any],
) -> None:
    """A protein whose miniprot alignment passes identity and coverage but lands on no gene
    (its gene removed from the model) is a leftover (§7.2: only failures move on)."""
    avail = protein_full["avail"]
    victim = protein_full["ids"][1]
    genes = pd.read_parquet(avail["genes"])
    assert victim in set(genes["locus_tag"])
    tmp = protein_full["tmp"] / "m16"
    tmp.mkdir()
    genes[genes["locus_tag"] != victim].to_parquet(tmp / "genes", index=False)
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        spec = reg.get("protein_leftovers")
        ins = {k: avail[k] for k in (*spec.inputs, *spec.optional_inputs) if k in avail}
        ins["genes"] = tmp / "genes"
        (out_key,) = spec.outputs
        reg.run_step("protein_leftovers", protein_full["cfg"], ins, {out_key: tmp / out_key})
    text = (tmp / out_key).read_text()
    left = {ln[1:].split()[0] for ln in text.splitlines() if ln.startswith(">")}
    assert left == {victim}


# ================================================================ items 7, 8, 9, 11: giap run


@pytest.fixture(scope="module")
def sim50(tmp_path_factory: pytest.TempPathFactory) -> Any:
    from giapsim.proteomics import ProteomicsConfig
    from giapsim.simulate import SimConfig, simulate

    out = tmp_path_factory.mktemp("sim50_review")
    return simulate(
        SimConfig(proteomics=ProteomicsConfig(), write_precomputed=True), 7, out / "sim"
    )


def _giap_run(sim: Any, out: Path, db: Path, *extra: str, config: Path | None = None) -> Any:
    from click.testing import CliRunner

    from giap.cli import cli

    args = ["run", "--genome", str(sim.genome_fasta), "--config", str(config or sim.config_yaml)]
    args += ["--outdir", str(out), "--modules", "core", "--seed", "7", "--executor", "local"]
    return CliRunner().invoke(cli, [*args, *extra], env={"GIAP_DB_DIR": str(db)})


def _no_module2(out: Path) -> bool:
    return not list(out.rglob("peptide_hits.tsv"))


def test_item7_no_evidence_exits_1(sim50: Any, tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    r = _giap_run(sim50, out, write_test_db(tmp_path / "db"), "--annotation", str(sim50.bakta_json))
    assert r.exit_code == 1, r.output
    assert "no evidence" in r.output
    assert _no_module2(out)


def test_item7_proteins_only_exits_4_naming_miniprot(sim50: Any, tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    r = _giap_run(
        sim50, out, write_test_db(tmp_path / "db"),
        "--proteins", str(sim50.proteins_faa), "--annotation", str(sim50.bakta_json),
    )  # fmt: skip
    assert r.exit_code == 4, r.output
    assert "miniprot" in r.output
    assert _no_module2(out)


def test_item7_gff3_annotation_without_bakta_json_layer_exits_1(sim50: Any, tmp_path: Path) -> None:
    assert sim50.annotation_gff3 is not None
    lines = sim50.config_yaml.read_text().splitlines()
    cfg = tmp_path / "no_bakta_layer.yaml"
    cfg.write_text("\n".join(ln for ln in lines if "bakta_json" not in ln) + "\n")
    out = tmp_path / "OUT"
    r = _giap_run(
        sim50, out, write_test_db(tmp_path / "db"),
        "--peptides", str(sim50.peptides), "--annotation", str(sim50.annotation_gff3), config=cfg,
    )  # fmt: skip
    assert r.exit_code == 1, r.output
    assert "bakta_json" in r.output
    assert _no_module2(out)


def _peptide_run(sim: Any, out: Path, db: Path, peptides: Path | None = None) -> Any:
    return _giap_run(
        sim, out, db, "--peptides", str(peptides or sim.peptides),
        "--annotation", str(sim.bakta_json),
    )  # fmt: skip


def test_item11_local_outputs_are_published_to_out_mapping(sim50: Any, tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    r = _peptide_run(sim50, out, write_test_db(tmp_path / "db"))
    assert r.exit_code == 0, r.output
    for f in MODULE2_FILES:
        assert (out / "mapping" / f).is_file(), f"OUT/mapping/{f} missing"


def test_item8_stale_work_files_are_ignored(sim50: Any, tmp_path: Path) -> None:
    db = write_test_db(tmp_path / "db")
    out = tmp_path / "OUT"
    r = _peptide_run(sim50, out, db)
    assert r.exit_code == 0, r.output
    first = {f: p.read_bytes() for f, p in module2_files(out).items() if f in MAPPING_TSVS}
    peps = pd.read_parquet(next(out.rglob("peptides.norm.parquet")))
    some = [p for p in peps["peptide"] if 8 <= len(p) <= 30][:5]
    work = out / "work"
    work.mkdir(exist_ok=True)
    (work / "proteins_long").write_text(">STALE1\n" + "MKT" * 30 + "\n")
    (work / "contaminants").write_text(">STALE_CONT\n" + "GG".join(some) + "\n")
    (work / "reference_peptides").write_text("peptide\tq_value\nSTALEPEPTIDEK\t0.001\n")
    r = _peptide_run(sim50, out, db)
    assert r.exit_code == 0, r.output
    again = module2_files(out)
    for f, data in first.items():
        assert again[f].read_bytes() == data, f
    un = read_tsv_df(again["unmapped.tsv"])
    assert not un["sequence_id"].astype(str).str.startswith("STALE").any()


TRYPSIN_PEPTIDES = ("LSSPATLNSR", "VATVSLPR")  # porcine trypsin, tests/data/contaminants_mini.faa


def test_item9_contaminant_library_is_staged_and_matched(sim50: Any, tmp_path: Path) -> None:
    lib = CONTAMINANTS_MINI.read_text().replace("\n", "")
    assert all(p in lib for p in TRYPSIN_PEPTIDES)
    src = pd.read_csv(sim50.peptides, sep="\t", dtype=str, keep_default_na=False)
    extra = []
    for i, p in enumerate(TRYPSIN_PEPTIDES):
        row = dict.fromkeys(src.columns, "NA")
        row |= {"peptide": p, "protein_ids": f"XPROT{i}", "q_value": "0.0001", "pep": "0.0001",
                "score": "5", "is_decoy": "false", "intensity": "10", "run": "run1",
                "charge": "2", "spectrum_id": f"run1.cont{i}", "source": "database",
                "search_db": "other"}  # fmt: skip
        extra.append(row)
    peps = tmp_path / "peptides_cont.tsv"
    pd.concat([src, pd.DataFrame(extra)], ignore_index=True).to_csv(peps, sep="\t", index=False)
    db = write_test_db(tmp_path / "db", contaminants=CONTAMINANTS_MINI)
    out = tmp_path / "OUT"
    r = _peptide_run(sim50, out, db, peptides=peps)
    assert r.exit_code == 0, r.output
    un = read_tsv_df(module2_files(out)["unmapped.tsv"]).set_index("sequence")
    for p in TRYPSIN_PEPTIDES:
        assert p in un.index, f"{p} not in unmapped.tsv"
        assert un.loc[p, "label"] == "contaminant", un.loc[p].to_dict()
        assert "P00761" in str(un.loc[p, "contaminant_match"]), un.loc[p].to_dict()


# ================================================================ decisions 13-15 (re-review)
# docs/specs/M3-17-mapping-steps.md, "Lead decisions after the re-review (PR #77)". Probes and
# mutants: local/scratch/review-m3-17/r2/ (probe E; N10, M10).


def _no_bakta_layer(sim: Any, tmp: Path) -> Path:
    cfg = tmp / "no_bakta_layer.yaml"
    lines = sim.config_yaml.read_text().splitlines()
    cfg.write_text("\n".join(ln for ln in lines if "bakta_json" not in ln) + "\n")
    return cfg


def test_d13_peptides_without_bakta_json_exits_1_naming_the_key(sim50: Any, tmp_path: Path) -> None:
    """Probe E: ``--peptides`` with neither ``--annotation`` nor
    ``annotation.precomputed.bakta_json``: the mapping targets are pruned only because the
    REQUIRED root ``bakta_json`` is missing → exit 1 naming it, never a silent exit 0."""
    out = tmp_path / "OUT"
    r = _giap_run(
        sim50, out, write_test_db(tmp_path / "db"), "--peptides", str(sim50.peptides),
        config=_no_bakta_layer(sim50, tmp_path),
    )  # fmt: skip
    assert r.exit_code == 1, r.output
    assert "bakta_json" in r.output, r.output
    assert _no_module2(out)


def test_d13_missing_optional_roots_are_still_pruned_silently(sim50: Any, tmp_path: Path) -> None:
    """No ``--proteins`` (``proteins_long``), no contaminant library in ``$GIAP_DB_DIR``
    (``contaminants``), no reference peptides: exit 0 with the Module 2 files, nothing named."""
    out = tmp_path / "OUT"
    r = _peptide_run(sim50, out, write_test_db(tmp_path / "db"))
    assert r.exit_code == 0, r.output
    for f in MODULE2_FILES:
        assert (out / "mapping" / f).is_file(), f"OUT/mapping/{f} missing"
    for key in ("proteins_long", "contaminants", "reference_peptides"):
        assert key not in r.output, r.output


def _search_db_events(events: Path) -> list[dict[str, Any]]:
    import json

    found = []
    for ln in events.read_text(encoding="utf-8").splitlines():
        if not ln.strip():
            continue
        ev = json.loads(ln)
        data = ev.get("data") if isinstance(ev.get("data"), dict) else {}
        if "search_db" in str(ev.get("event", "")) or "search_db" in data:
            found.append(ev)
    return found


def test_d14_search_db_auto_resolution_is_an_event(sim50: Any, tmp_path: Path) -> None:
    """``search_db`` absent (= ``auto``, C9) on a same-strain search (probe A):
    ``OUT/logs/events.jsonl`` holds an event naming ``search_db`` (event name or data key) whose
    data carries the resolved value ``same``."""
    src = pd.read_csv(sim50.peptides, sep="\t", dtype=str, keep_default_na=False)
    src = src.drop(columns=["search_db"])
    peps = tmp_path / "peptides_auto.tsv"
    src.to_csv(peps, sep="\t", index=False)
    out = tmp_path / "OUT"
    r = _peptide_run(sim50, out, write_test_db(tmp_path / "db"), peptides=peps)
    assert r.exit_code == 0, r.output
    events = out / "logs" / "events.jsonl"
    assert events.is_file()
    evs = _search_db_events(events)
    assert evs, events.read_text()[-2000:]
    assert any("same" in {str(v) for v in (e.get("data") or {}).values()} for e in evs), evs


N_MISMATCH = 40  # every 40th target peptide (2.5 %): the table still resolves to ``same``


@pytest.fixture(scope="module")
def auto_chain(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """Same-strain sim, ``search_db: auto`` on every row; every 40th target peptide's
    ``protein_ids`` is swapped for another gene-model locus tag (still in the database, so the
    table resolves to ``same``) → its hits on its true gene are ``engine_protein_mismatch``."""
    tmp = tmp_path_factory.mktemp("auto_chain")
    run = preflight(_sim(tmp, SEED, search_db="same"), tmp / "pf")
    df = pd.read_parquet(run.peptides)
    df["search_db"] = pd.Series(["auto"] * len(df), dtype="string")
    tgt = list(df.index[~df["is_decoy"].astype(bool)])
    ids = {i: [str(x) for x in df.at[i, "protein_ids"]] for i in df.index}
    swapped: dict[str, set[str]] = {}
    for k in range(0, len(tgt), N_MISMATCH):
        i = tgt[k]
        other = next(
            (ids[j] for j in tgt[k + 1 :] + tgt[:k] if ids[j] and not set(ids[j]) & set(ids[i])),
            None,
        )
        if other is not None:
            ids[i] = list(other)
            swapped[str(df.at[i, "peptide"])] = set(other)
    df["protein_ids"] = [ids[i] for i in df.index]
    run = _with_parquet(run, df, tmp / "peptides.auto.parquet")
    cfg = resolve_config(run.config_yaml, [f"seed={SEED}"])
    avail, _ = _drive_logged(cfg, tmp / "work", roots_of(run))
    return {"cfg": cfg, "avail": avail, "tmp": tmp, "swapped": swapped,
            "files": module2_files(tmp / "work", avail)}  # fmt: skip


def test_d15_n10_loci_flags_mismatch_on_auto_resolved_same(auto_chain: dict[str, Any]) -> None:
    """N10: ``loci`` must use the resolved queries, so ``auto`` → ``same`` reaches M3-14's
    ``engine_protein_mismatch`` check (a hit gene not among the engine's protein IDs)."""
    hits = read_tsv_df(auto_chain["files"]["peptide_hits.tsv"])
    sw = auto_chain["swapped"]
    rows = hits[hits["peptide"].isin(sw.keys()) & hits["gene"].notna()]
    rows = rows[[g not in sw[p] for p, g in zip(rows["peptide"], rows["gene"], strict=True)]]
    assert len(rows) > 0, "no swapped peptide hit a gene"
    flagged = rows["flags"].map(lambda f: "engine_protein_mismatch" in _flags(f))
    assert flagged.all(), rows.loc[~flagged, ["peptide", "gene", "flags"]]
    rest = hits[~hits["peptide"].isin(sw.keys())]
    assert len(rest) > 0


def test_d15_m10_rejected_placements_are_matched_against_contaminants(
    auto_chain: dict[str, Any],
) -> None:
    """M10: every query not ACCEPTED (incl. a target hit rejected by placement FDR) is matched
    against the contaminant library; such a peptide is labelled ``contaminant`` with reason
    ``rejected_placement_fdr`` and the library entry in ``contaminant_match``."""
    avail = auto_chain["avail"]
    per = pd.read_parquet(avail["placement"])
    state = per.groupby("peptide")["accepted"].agg(lambda s: bool(s.astype(bool).any()))
    rejected = [p for p, a in state.items() if not a and len(str(p)) >= 8][:3]
    assert rejected, "the chain has no target hit rejected by placement FDR"
    tmp = auto_chain["tmp"] / "m10"
    tmp.mkdir()
    lib = tmp / "contaminants.fasta"
    lib.write_text("".join(f">CONT_{i:02d} test\nGGG{p}GGG\n" for i, p in enumerate(rejected)))
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        spec = reg.get("loci")
        ins = {k: avail[k] for k in (*spec.inputs, *spec.optional_inputs) if k in avail}
        ins["contaminants"] = lib
        (out_key,) = spec.outputs
        reg.run_step("loci", auto_chain["cfg"], ins, {out_key: tmp / out_key})
    un = read_tsv_df(module2_files(tmp / out_key)["unmapped.tsv"]).set_index("sequence")
    for i, p in enumerate(rejected):
        assert p in un.index, f"{p} not in unmapped.tsv"
        row = un.loc[p]
        assert row["label"] == "contaminant", row.to_dict()
        assert row["reason"] == "rejected_placement_fdr", row.to_dict()
        assert str(row["contaminant_match"]) == f"CONT_{i:02d}", row.to_dict()
