"""M4-06 lead decisions after code review (PR #51), items 1-3 of
docs/specs/M4-06-input-contract.md "Lead decisions after code review".

1. C10 follows merged M1-12: ``abundance.norm.parquet`` has ``protein_id``, ``raw_id`` (optional)
   and ``value``; the scale comes from ``preflight.json`` (``routing.abundance_scale``) into
   ``M4Inputs.abundance_scale``. Pinned against the real output of ``giap preflight``.
2. ``cterm_extension`` (ADR-0013) is a locus class that assigns its gene.
3. The four mutants that survived review: layers rows for unknown genes; an AMR region with
   exactly one NA coordinate; ``1.5`` in an object-dtype int column; TSV bools other than
   ``true``/``false``.

Loader key for the report: ``paths["preflight"]`` (the M1-12 ``OUT/preflight/preflight.json``).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner
from m4_contract_fixture import _hit, _hits_frame, build_frames, make_inputs, tag, write_work

from giap.io.issues import Issue

LOW = ["--set", "preflight.min_n50_warn=1", "--set", "preflight.min_largest_contig_warn=1"]


def errors(issues: list[Issue]) -> list[Issue]:
    return [i for i in issues if i.severity == "error"]


def text(i: Issue) -> str:
    return i.message + " " + " ".join(i.examples)


def with_code(issues: list[Issue], code: str) -> list[Issue]:
    return [i for i in issues if i.code == code]


def _contract() -> Any:
    from giap.genetable import contract

    return contract


# --------------------------------------------------------------------------- 1. C10 vs M1-12


@pytest.fixture
def preflight_out(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Run merged M1-12 ``giap preflight`` on a small genome with peptides and abundance."""
    from giap.cli import cli

    for d in ("db_empty", "repo_empty"):
        (tmp_path / d).mkdir()
    monkeypatch.setenv("GIAP_DB_DIR", str(tmp_path / "db_empty"))
    monkeypatch.setenv("GIAP_REPO", str(tmp_path / "repo_empty"))
    rng = np.random.default_rng(7)
    seq = "".join(np.array(list("ACGT"))[rng.integers(0, 4, 120_000)])
    genome = tmp_path / "g.fa"
    genome.write_text(f">c circular=true\n{seq}\n", encoding="utf-8")
    peps = ["peptide\tprotein_ids\trun\tintensity\tq_value\tis_decoy"]
    for k, p in enumerate(("PEPTIDEK", "AAGLIEVQAPILSR", "LLDEGKTVR", "VQAPILSRGDGTK")):
        ids = f"P{k}" + (f";P{k + 1}" if k == 0 else "")
        peps.append(f"{p}\t{ids}\tr1\t1000\t0.001\tfalse")
        peps.append(f"{p}\t{ids}\tr2\t{'' if k == 2 else 500}\t0.002\tfalse")
    peps.append("KAAADITPEP\tREV_P1\tr1\t10\t0.5\ttrue")
    (tmp_path / "pep.tsv").write_text("\n".join(peps) + "\n", encoding="utf-8")
    ab = ["protein_id\tvalue"] + [f"P{k}\t{k % 7 + 1}.5" if k % 5 else f"P{k}\t" for k in range(12)]
    (tmp_path / "abund.tsv").write_text("\n".join(ab) + "\n", encoding="utf-8")
    out = tmp_path / "out"
    args = ["preflight", "--genome", str(genome), "--peptides", str(tmp_path / "pep.tsv")]
    args += ["--abundance", str(tmp_path / "abund.tsv"), "--outdir", str(out), *LOW]
    r = CliRunner().invoke(cli, args)
    assert (out / "work" / "preflight" / "abundance.norm.parquet").is_file(), (
        r.exit_code,
        r.output,
        r.exception,
    )
    return out


def test_c10_columns_follow_m1_12() -> None:
    cols = {c.name: c for c in _contract().TABLES["abundance"]}
    assert {n for n, c in cols.items() if not c.optional} == {"protein_id", "value"}
    assert "scale" not in cols or cols["scale"].optional
    assert "raw_id" in cols and cols["raw_id"].optional and cols["raw_id"].dtype == "str"


@pytest.mark.parametrize("with_raw_id", [True, False])
def test_c10_raw_id_optional_and_no_scale_column(with_raw_id: bool) -> None:
    df = pd.DataFrame(
        {
            "protein_id": pd.array(["P1", "P2"], dtype=pd.StringDtype()),
            "raw_id": pd.array(["p1", "p2"], dtype=pd.StringDtype()),
            "value": [1.5, np.nan],
        }
    )
    if not with_raw_id:
        df = df.drop(columns="raw_id")
    assert errors(_contract().validate_table("abundance", df)) == []


def test_real_preflight_abundance_passes_contract(preflight_out: Path) -> None:
    real = pd.read_parquet(preflight_out / "work" / "preflight" / "abundance.norm.parquet")
    assert "raw_id" in real.columns and "scale" not in real.columns  # M1-12 as merged
    assert errors(_contract().validate_table("abundance", real)) == []


def test_loader_takes_scale_from_preflight_json(preflight_out: Path, tmp_path: Path) -> None:
    from giap.io.work_tables import M4Inputs, default_paths, load_m4_inputs

    report = json.loads((preflight_out / "preflight" / "preflight.json").read_text("utf-8"))
    scale = report["routing"]["abundance_scale"]
    assert scale in ("log2", "linear")
    paths: dict[str, Path] = dict(default_paths(tmp_path / "work"))
    write_work(build_frames(20), paths, skip=("abundance",))
    paths["abundance"] = preflight_out / "work" / "preflight" / "abundance.norm.parquet"
    paths["preflight"] = preflight_out / "preflight" / "preflight.json"
    inputs, issues = load_m4_inputs(paths)
    assert errors(issues) == [], [text(i) for i in errors(issues)]
    assert "abundance_scale" in {f for f in M4Inputs.__dataclass_fields__}
    assert inputs.abundance is not None and len(inputs.abundance) == 12
    assert inputs.abundance_scale == scale


@pytest.mark.parametrize("scale", ["log2", "linear"])
def test_loader_scale_from_fixture_preflight_json(tmp_path: Path, scale: str) -> None:
    from giap.io.work_tables import default_paths, load_m4_inputs

    f = build_frames(20)
    f["abundance_scale"] = scale
    paths: dict[str, Path] = dict(default_paths(tmp_path / "work"))
    paths["preflight"] = tmp_path / "preflight" / "preflight.json"
    write_work(f, paths)
    inputs, issues = load_m4_inputs(paths)
    assert errors(issues) == []
    assert inputs.abundance_scale == scale


# --------------------------------------------------------------------------- 2. cterm_extension


def test_cterm_extension_is_an_assigning_locus_class() -> None:
    c = _contract()
    assert "cterm_extension" in c.LOCUS_CLASSES
    assert "cterm_extension" in c.ASSIGNING_CLASSES
    col = next(x for x in c.TABLES["peptide_hits"] if x.name == "locus_class")
    assert col.enum is not None and "cterm_extension" in col.enum


def _cterm_frames(gene: str | None) -> dict[str, Any]:
    f = build_frames(20)
    g = f["genes"]
    i = 9
    start = int(g.loc[i, "end"]) - 6 if g.loc[i, "strand"] == "+" else int(g.loc[i, "start"]) - 18
    row = _hit("CTERMEXTK", start, str(g.loc[i, "strand"]), "cterm_extension", "unique", gene)
    f["peptide_hits"] = pd.concat([f["peptide_hits"], _hits_frame([row])], ignore_index=True)
    return f


def _all_issues(f: dict[str, Any]) -> list[Issue]:
    c = _contract()
    return [
        *c.validate_table("peptide_hits", f["peptide_hits"]),
        *c.validate_inputs(make_inputs(f)),
    ]


def test_gene_assigned_cterm_extension_row_is_valid() -> None:
    issues = _all_issues(_cterm_frames(tag(9)))
    assert errors(issues) == [], [text(i) for i in errors(issues)]


def test_cterm_extension_without_gene_is_assignment_error() -> None:
    hit = with_code(_all_issues(_cterm_frames(None)), "contract.peptide_hits.assignment")
    assert hit and all(i.severity == "error" for i in hit)


# --------------------------------------------------------------------------- 3. surviving mutants


def test_layers_row_for_unknown_gene_is_error() -> None:
    f = build_frames(20)
    lay = f["layers"]
    extra = lay.iloc[[0]].copy()
    extra["locus_tag"] = "GSIM_99999"
    f["layers"] = pd.concat([lay, extra], ignore_index=True)
    hit = with_code(_contract().validate_inputs(make_inputs(f)), "contract.layers.gene_coverage")
    assert hit and all(i.severity == "error" for i in hit)
    assert any("GSIM_99999" in text(i) for i in hit)


@pytest.mark.parametrize("which", ["start", "end"])
def test_amr_region_with_one_na_coordinate_is_error(which: str) -> None:
    f = build_frames(20)
    r = f["regions"]
    i = r.index[r["region_type"] == "amr"][0]
    assert pd.isna(r.loc[i, "start"]) and pd.isna(r.loc[i, "end"])
    r.loc[i, which] = 1000
    hit = with_code(_contract().validate_table("regions", r), "contract.regions.amr_coordinates")
    assert hit and hit[0].severity == "error" and hit[0].n_rows == 1


def test_fractional_float_in_object_int_column_is_dtype_error() -> None:
    g = build_frames(20)["genes"]
    g["start"] = g["start"].astype(object)
    g.loc[2, "start"] = 1.5
    hit = with_code(_contract().validate_table("genes", g), "contract.genes.dtype")
    assert hit and any("start" in text(i) for i in hit)
    assert any("3" in i.examples for i in hit)  # 1-based row of position 2


def test_integral_float_in_object_int_column_is_accepted() -> None:
    g = build_frames(20)["genes"]
    g["start"] = g["start"].astype(object)
    g.loc[2, "start"] = float(g.loc[2, "start"])
    assert with_code(_contract().validate_table("genes", g), "contract.genes.dtype") == []


def _rewrite_tsv_cell(path: Path, column: str, row: int, value: str) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    data = [k for k, ln in enumerate(lines) if k > 0 and not ln.startswith("#")]
    cells = lines[data[row]].split("\t")
    cells[header.index(column)] = value
    lines[data[row]] = "\t".join(cells)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.mark.parametrize("bad", ["TRUE", "1", "True", "0", "FALSE"])
def test_tsv_bool_accepts_only_lowercase_true_false(tmp_path: Path, bad: str) -> None:
    from giap.io.work_tables import default_paths, load_m4_inputs

    paths = default_paths(tmp_path / "work")
    write_work(build_frames(20), paths)
    assert "crosses_origin" in paths["peptide_hits"].read_text(encoding="utf-8")
    _rewrite_tsv_cell(paths["peptide_hits"], "crosses_origin", 0, bad)
    _, issues = load_m4_inputs(paths)
    hit = with_code(issues, "contract.peptide_hits.dtype")
    assert hit and hit[0].severity == "error"
    assert any("crosses_origin" in text(i) for i in hit)


def test_tsv_bool_lowercase_is_valid(tmp_path: Path) -> None:
    from giap.io.work_tables import default_paths, load_m4_inputs

    paths = default_paths(tmp_path / "work")
    write_work(build_frames(20), paths)
    _rewrite_tsv_cell(paths["peptide_hits"], "crosses_origin", 0, "false")
    inputs, issues = load_m4_inputs(paths)
    assert with_code(issues, "contract.peptide_hits.dtype") == []
    assert not bool(inputs.peptide_hits["crosses_origin"].iloc[0])
