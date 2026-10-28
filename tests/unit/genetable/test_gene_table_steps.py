"""M4-15 steps ``gene_table`` and ``detectability`` (docs/specs/M4-15-gene-table-assembly.md,
lead decisions 1-3; GIPPA_SPEC.md §6.2, §6.6, §9.4; M2-02 registry).

Pinned by the test author (see the hand-back report):
- Registry keys follow the existing producers: ``gene_table`` reads ``genes`` / ``layers`` /
  ``mapping`` (the ``loci`` output directory holding ``peptide_hits.tsv`` and friends) and may read
  only keys of ``KNOWN_INPUTS``; it writes the key ``genes_pre``; the COG 2020 table arrives as the
  optional key ``cog2020_def`` (the setup item name, as ``contaminants`` does for ``loci``) — the
  raw NCBI ``cog-20.def.tab``. ``detectability`` reads ``genes_pre``.
- Without the COG table ``cog_category`` is NA and a warning carrying ``layer_missing:cog`` is
  logged (the M2-19/M2-20 convention for run flags).
- ``run_gene_table(work_dir, cfg)`` / ``run_detectability(work_dir, cfg)`` find the contract files
  in ``work_dir`` in the real ``gippa run`` layout (M3 TSVs under ``mapping/``, C14 under
  ``annotation/``) and in the flat ``default_paths`` layout (lead decision 8), and return
  ``work_dir/genetable/genes_pre.parquet`` / ``work_dir/genetable/genes.parquet``.
- A contract ``error`` issue makes the step fail (``StepError`` through ``run_step``; the spec's
  "exit code 1" is reported as an inconsistency, the runner maps ``StepError`` to 4).
"""

from __future__ import annotations

import logging
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gene_table_fixture as gt
import m4_contract_fixture as mf
from gene_table_fixture import as_list, is_na, tag

KNOWN_INPUTS = {
    "genes", "replicons", "xrefs", "annotation_meta", "tus", "layers", "regions", "region_genes",
    "peptides", "abundance", "mapping", "id_map", "preflight", "cog2020_def",
}  # fmt: skip


def _reg() -> Any:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    return reg


def _frames(dup: bool = False) -> dict[str, Any]:
    fr = gt.frames_with()
    fr["genes"].loc[5, "partial"] = True
    fr["layers"].at[5, "flags"] = ["partial"]
    if dup:
        fr["genes"] = pd.concat([fr["genes"], fr["genes"].iloc[[3]]], ignore_index=True)
    return fr


def write_work(root: Path, frames: dict[str, Any], layout: str = "both") -> dict[str, Path]:
    """Contract files under ``root``; ``layout`` ``flat`` (``default_paths``), ``run`` (M3 TSVs
    under ``mapping/``, C14 under ``annotation/``) or ``both``. Returns the step input keys."""
    from gippa.io.work_tables import default_paths

    flat = default_paths(root) | {"preflight": root / "preflight.json"}
    mf.write_work(frames, flat)
    moved = {"peptide_hits.tsv": "mapping", "protein_hits.tsv": "mapping",
             "unmapped.tsv": "mapping", "id_map.tsv": "annotation"}  # fmt: skip
    for name, sub in moved.items():
        if layout != "flat":
            (root / sub).mkdir(exist_ok=True)
            shutil.copyfile(root / name, root / sub / name)
        if layout == "run":
            (root / name).unlink()
    keys = {k: p for k, p in flat.items() if k not in ("peptide_hits", "protein_hits", "unmapped")}
    return keys | {"mapping": root / "mapping", "id_map": root / "annotation" / "id_map.tsv"}


def run_gene_table_step(
    root: Path, frames: dict[str, Any], cog: Path | None = None, seed: int = 1
) -> Path:
    from gippa.config import GippaConfig

    reg = _reg()
    spec = reg.get("gene_table")
    keys = write_work(root / "work", frames, layout="run")
    if cog is not None:
        keys["cog2020_def"] = cog
    allowed = set(spec.inputs) | set(spec.optional_inputs)
    inputs = {k: p for k, p in keys.items() if k in allowed}
    out = root / "out"
    outputs = {k: out / k for k in spec.outputs}
    outputs["genes_pre"] = out / "genetable" / "genes_pre.parquet"
    reg.run_step("gene_table", GippaConfig(seed=seed), inputs, outputs)
    return outputs["genes_pre"]


# ----------------------------------------------------------------------------- registry


def test_steps_registered() -> None:
    reg = _reg()
    gtab, det = reg.get("gene_table"), reg.get("detectability")
    for s in (gtab, det):
        assert (s.kind, s.module, s.precomputed) == ("py", "core", None)
        assert s.func is not None
    assert {"genes", "layers", "mapping"} <= set(gtab.inputs)
    assert set(gtab.inputs) | set(gtab.optional_inputs) <= KNOWN_INPUTS
    assert "cog2020_def" in gtab.optional_inputs
    assert "genes_pre" in gtab.outputs
    assert "genes_pre" in det.inputs
    assert reg.topo_order(["detectability", "gene_table", "loci", "context_layers"]) == [
        "context_layers", "loci", "gene_table", "detectability",
    ]  # fmt: skip


# ----------------------------------------------------------------------------- gene_table


def test_gene_table_step_with_cog_table(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    cog = gt.write_cog_def(tmp_path / "db" / "cog" / "cog-20.def.tab")
    with caplog.at_level(logging.WARNING):
        path = run_gene_table_step(tmp_path, _frames(), cog=cog)
    assert not [r for r in caplog.records if "layer_missing:cog" in r.getMessage()]
    df = pd.read_parquet(path)
    assert len(df) == 20 and df["locus_tag"].tolist() == [tag(i) for i in range(20)]
    assert set(gt.GENES_COLUMNS) - {"p_detect"} <= set(df.columns)
    assert "aa" not in df.columns and "nt" not in df.columns
    by = df.set_index("locus_tag")
    assert (by.loc[tag(0), "cog_id"], by.loc[tag(0), "cog_category"]) == ("COG0000", "J")
    assert (by.loc[tag(2), "cog_id"], by.loc[tag(2), "cog_category"]) == ("COG0002", "KT")
    assert is_na(by.loc[tag(4), "cog_category"])  # COG0004 is not in the table
    assert is_na(by.loc[tag(1), "cog_id"]) and is_na(by.loc[tag(1), "cog_category"])
    assert int(by.loc[tag(0), "start"]) == 100  # internal 0-based in the parquet
    assert by.loc[tag(5), "detection_status"] == "partial_excluded"
    if "p_detect" in df.columns:
        assert df["p_detect"].isna().all()  # genes_pre: everything but p_detect
    assert not any("detectability_constant" in as_list(f) for f in df["flags"])


def test_gene_table_step_without_cog_table(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Lead decision 3: no COG 2020 table -> ``cog_category`` NA and run flag
    ``layer_missing:cog``."""
    with caplog.at_level(logging.WARNING):
        path = run_gene_table_step(tmp_path, _frames())
    msgs = [r.getMessage() for r in caplog.records]
    assert any("layer_missing:cog" in m for m in msgs), msgs
    df = pd.read_parquet(path)
    assert df["cog_category"].isna().all()
    assert df.set_index("locus_tag").loc[tag(0), "cog_id"] == "COG0000"


def test_gene_table_step_fails_on_contract_error(tmp_path: Path) -> None:
    reg = _reg()
    with pytest.raises(reg.StepError) as info:
        run_gene_table_step(tmp_path, _frames(dup=True))
    assert "contract" in str(info.value) or "locus_tag" in str(info.value)
    assert not (tmp_path / "out" / "genetable" / "genes_pre.parquet").exists()


# ----------------------------------------------------------------------------- run_* functions


def _assert_nulls_kept(pre: Path, out: Path) -> None:
    """Lead decision 11: string NAs of genes_pre are still null in genes.parquet."""
    import pyarrow.parquet as pq

    tp, to = pq.read_table(pre), pq.read_table(out)
    n = to.num_rows
    for c in ("island_id", "replichore"):
        assert to.column(c).null_count == n, c
    for c in ("cog_id", "cog_category", "gene", "product"):
        values = to.column(c).to_pylist()
        assert "nan" not in [str(v).lower() for v in values if v is not None], c
        assert to.column(c).null_count == tp.column(c).null_count > 0, c


def _setup_db(root: Path, *, listed: bool, raw: bytes | None = None) -> Path:
    """A tmp ``$GIPPA_DB_DIR``: fixture ``cog/cog-20.def.tab`` and a ``gippa setup`` manifest that
    lists ``cog2020_def`` (``listed``) or only another item."""
    import json

    db = root / "db"
    cog = gt.write_cog_def(db / "cog" / "cog-20.def.tab")
    if raw is not None:
        cog.write_bytes(raw)
    item = {"name": "cog2020_def" if listed else "kegg_ko_list", "version": "2020", "url": "",
            "checksum": "", "size_bytes": cog.stat().st_size, "retrieved": "",
            "path": "cog/cog-20.def.tab" if listed else "kegg/ko_list.tsv", "status": "ok",
            "message": ""}  # fmt: skip
    doc = {"schema_version": "1.0", "items": [item]}
    (db / "db_manifest.json").write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    return db


@pytest.mark.parametrize("cp1252", [False, True])
def test_cog_table_found_through_setup_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture,
    cp1252: bool,
) -> None:  # fmt: skip
    """Lead decisions 12/13: ``$GIPPA_DB_DIR`` manifest lists ``cog2020_def`` -> categories set
    (also when the file has a cp1252 byte, as the real NCBI file does)."""
    from gippa.config import GippaConfig
    from gippa.steps.gene_table import run_gene_table

    raw = None
    if cp1252:
        raw = b"".join("\t".join(r).encode("utf-8") + b"\n" for r in gt.COG_DEF_ROWS)
        raw = raw.replace(b"DnaA\t", b"DnaA\x92\t", 1)
    monkeypatch.setenv("GIPPA_DB_DIR", str(_setup_db(tmp_path, listed=True, raw=raw)))
    work = tmp_path / "work"
    work.mkdir()
    write_work(work, _frames(), layout="run")
    with caplog.at_level(logging.WARNING):
        pre = run_gene_table(work, GippaConfig(seed=3))
    assert not [r for r in caplog.records if "layer_missing:cog" in r.getMessage()]
    by = pd.read_parquet(pre).set_index("locus_tag")
    assert (by.loc[tag(0), "cog_category"], by.loc[tag(2), "cog_category"]) == ("J", "KT")


def test_cog_file_without_manifest_entry_is_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """Lead decision 12: the file alone (setup never recorded it) -> NA + ``layer_missing:cog``."""
    from gippa.config import GippaConfig
    from gippa.steps.gene_table import run_gene_table

    monkeypatch.setenv("GIPPA_DB_DIR", str(_setup_db(tmp_path, listed=False)))
    work = tmp_path / "work"
    work.mkdir()
    write_work(work, _frames(), layout="run")
    with caplog.at_level(logging.WARNING):
        pre = run_gene_table(work, GippaConfig(seed=3))
    assert any("layer_missing:cog" in r.getMessage() for r in caplog.records)
    assert pd.read_parquet(pre)["cog_category"].isna().all()


def test_m4_paths_preflight_directory(tmp_path: Path) -> None:
    """Lead decision 20(a): a ``preflight`` input that is a directory (the local executor's
    ``work/preflight/``) resolves to the ``preflight.json`` inside it."""
    from gippa.steps.gene_table import m4_paths

    pre = tmp_path / "preflight"
    pre.mkdir()
    (pre / "preflight.json").write_text("{}\n", encoding="utf-8")
    paths = m4_paths({"preflight": pre, "mapping": tmp_path / "mapping"})
    assert paths["preflight"] == pre / "preflight.json"
    assert m4_paths({"preflight": pre / "preflight.json"})["preflight"] == pre / "preflight.json"


@pytest.mark.parametrize("as_dir", [True, False])
def test_abundance_scale_read_from_preflight(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, as_dir: bool
) -> None:
    """Lead decision 20(a): with no peptide intensities the C10 route needs the preflight scale
    (``abundance.scale = auto``); it is read from ``preflight.json`` given as a file or as the
    directory holding it, and no ``contract.preflight.abundance_scale`` issue is logged."""
    import json

    from gippa.config import GippaConfig

    fr = _frames()
    fr["peptides"]["intensity"] = float("nan")  # force the protein (C10) route
    reg = _reg()
    spec = reg.get("gene_table")
    keys = write_work(tmp_path / "work", fr, layout="run")
    pre_dir = tmp_path / "work" / "preflight_dir"
    pre_dir.mkdir()
    (pre_dir / "preflight.json").write_text(
        json.dumps({"routing": {"abundance_scale": "log2"}}) + "\n", encoding="utf-8"
    )
    keys["preflight"] = pre_dir if as_dir else pre_dir / "preflight.json"
    allowed = set(spec.inputs) | set(spec.optional_inputs)
    inputs = {k: p for k, p in keys.items() if k in allowed}
    assert "preflight" in inputs and "abundance" in inputs
    out = tmp_path / "out" / "genes_pre.parquet"
    outputs = {k: tmp_path / "out" / k for k in spec.outputs} | {"genes_pre": out}
    with caplog.at_level(logging.WARNING):
        reg.run_step("gene_table", GippaConfig(seed=1), inputs, outputs)
    assert not [r for r in caplog.records if "contract.preflight" in r.getMessage()]
    by = pd.read_parquet(out).set_index("locus_tag")
    assert set(by["abundance_type"]) != {"none"}
    assert float(by.loc[tag(0), "abundance_raw"]) == 20.0  # C10 value, log2 scale kept


@pytest.fixture
def no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Point ``$GIPPA_DB_DIR`` at an empty directory: tests never read the real databases."""
    empty = tmp_path / "empty_db"
    empty.mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(empty))
    return empty


@pytest.mark.parametrize("layout", ["run", "flat"])
def test_run_functions_write_genetable_files(
    tmp_path: Path, no_db: Path, caplog: pytest.LogCaptureFixture, layout: str
) -> None:
    from gippa.config import GippaConfig
    from gippa.steps.detectability import run_detectability
    from gippa.steps.gene_table import run_gene_table

    work = tmp_path / "work"
    work.mkdir()
    write_work(work, _frames(), layout=layout)
    cfg = GippaConfig(seed=3)
    with caplog.at_level(logging.WARNING):
        pre = run_gene_table(work, cfg)
    assert pre == work / "genetable" / "genes_pre.parquet" and pre.is_file()
    pre_df = pd.read_parquet(pre)
    assert len(pre_df) == 20
    assert pre_df.set_index("locus_tag").loc[tag(3), "source_id"] == "orig_3"  # C14 found
    assert any("layer_missing:cog" in r.getMessage() for r in caplog.records)

    out = run_detectability(work, cfg)
    assert out == work / "genetable" / "genes.parquet" and out.is_file()
    _assert_nulls_kept(pre, out)
    df = pd.read_parquet(out)
    assert set(gt.GENES_COLUMNS) <= set(df.columns)
    by = df.set_index("locus_tag")
    eligible = ~by["partial"].astype(bool)
    # 20 genes < 50 per class: M4-13 constant p_detect = detected share of eligible genes
    rate = float(by.loc[eligible, "detected"].astype(bool).mean())
    assert np.allclose(by.loc[eligible, "p_detect"].astype(float), rate)
    assert pd.isna(by.loc[tag(5), "p_detect"])
    for t in by.index:
        flags = as_list(by.loc[t, "flags"])
        assert ("detectability_constant" in flags) is bool(eligible[t]), (t, flags)
    assert as_list(by.loc[tag(5), "flags"]) == ["partial"]

    parts = work / "genetable" / "detection_summary_parts" / "detectability.parquet"
    diag = pd.read_parquet(parts)
    assert len(diag) > 0 and set(diag["section"]) == {"detectability"}
