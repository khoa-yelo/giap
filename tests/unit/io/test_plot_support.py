"""``gippa.io.plot_support.write_support`` and the support schemas (M9-03 acceptance: example 8;
gzip above 1 MB with ``mtime=0``; schema validation with frictionless; ADR-M9-A rules).

Interface assumed beyond the spec text (M9-03 test author):

- ``write_support(name, df, plots_dir)`` writes ``plots_dir/support/<name>.tsv`` (or
  ``<name>.tsv.gz`` when the TSV is larger than 1 MB), creating ``support/``; the frame is
  already in its written form (1-based, see tests/unit/viz/test_support.py), §9.2 cells
  (``NA``, ``%.6g`` floats, integers as integers) through ``gippa.io.tsv``.
- Exactly one of ``<name>.tsv`` / ``<name>.tsv.gz`` exists after a write (a stale twin would be
  read first by ``RunView.support``, which tries ``.tsv`` before ``.tsv.gz``).
- Schemas live at ``src/gippa/schemas/plot_support/<name>.json``; ``schemas/datapackage.json``
  has one resource per support table named ``<name>`` with ``role: plot_support`` and a ``path``
  under ``plots/support/``, whose ``schema`` resolves (relative to ``src/gippa/schemas``) to that
  file. The 1 MB threshold is tested well clear of both readings (10^6 and 2^20 bytes).
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

import gippa.schemas

SCHEMAS = Path(gippa.schemas.__file__).resolve().parent
NAMES = ("replicon_landmarks", "gc_skew")


def _write(name: str, df: pd.DataFrame, plots_dir: Path) -> Any:
    from gippa.io.plot_support import write_support

    return write_support(name, df, plots_dir)


def landmarks() -> pd.DataFrame:
    from gippa.viz.support import replicon_landmarks

    reps = pd.DataFrame(
        [
            ("chr1", 100_000, "circular", "chromosome"),
            ("plasmid_1", 8_000, "circular", "plasmid"),
            ("lin", 3_000, "linear", "unknown"),
        ],
        columns=["contig", "length", "topology", "kind"],
    )
    feats = pd.DataFrame(
        [("F1", "oric", "chr1", 12_445, 12_677, ".", None, "oriC", None, None)],
        columns=["feature_id", "type", "contig", "start", "end", "strand", "gene", "product",
                 "amino_acid", "anticodon"],
    )  # fmt: skip
    return replicon_landmarks(reps, feats)


def skew(n_bp: int = 25_000, window: int = 10_000, seed: int = 7) -> pd.DataFrame:
    from gippa.viz.support import gc_skew

    rng = np.random.default_rng(seed)
    seq = "".join(rng.choice(list("ACGT"), size=n_bp))
    if n_bp > 2 * window:
        seq = seq[:window] + "N" * window + seq[2 * window :]  # window 2 is all N (NA)
    reps = pd.DataFrame(
        [("chr1", len(seq), "circular", "chromosome")],
        columns=["contig", "length", "topology", "kind"],
    )
    return gc_skew({"chr1": seq}, reps, window)


def big_skew_frame(rows: int = 60_000) -> pd.DataFrame:
    """A gc_skew-shaped frame of well over 2 MB of TSV."""
    rng = np.random.default_rng(3)
    start = np.arange(rows, dtype=np.int64) * 10 + 1
    vals = rng.uniform(-1, 1, size=rows)
    return pd.DataFrame(
        {
            "contig": "chromosome_with_a_long_name",
            "start": start,
            "end": start + 9,
            "gc_skew": vals,
            "cumulative_gc_skew": np.cumsum(vals),
        }
    )


def _support_files(plots_dir: Path) -> list[str]:
    base = plots_dir / "support"
    return sorted(p.name for p in base.iterdir()) if base.is_dir() else []


# ---------------------------------------------------------------- writing


def test_landmarks_written_under_support_as_plain_tsv(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    _write("replicon_landmarks", landmarks(), plots)
    assert _support_files(plots) == ["replicon_landmarks.tsv"]
    lines = (plots / "support" / "replicon_landmarks.tsv").read_text(encoding="utf-8").split("\n")
    assert lines[0] == "contig\tlength_bp\ttopology\toric_position\tter_position\tstatus"
    assert lines[1] == "chr1\t100000\tcircular\t12562\t62562\tcalled"
    assert lines[2] == "plasmid_1\t8000\tcircular\tNA\tNA\tno_oric"
    assert lines[3] == "lin\t3000\tlinear\tNA\tNA\tlinear"
    assert lines[4:] == [""]  # one trailing newline, \n endings


def test_gc_skew_na_and_integer_cells(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    _write("gc_skew", skew(), plots)
    text = (plots / "support" / "gc_skew.tsv").read_text(encoding="utf-8")
    rows = [line.split("\t") for line in text.splitlines()]
    assert rows[0] == ["contig", "start", "end", "gc_skew", "cumulative_gc_skew"]
    assert [r[1:3] for r in rows[1:]] == [["1", "10000"], ["10001", "20000"], ["20001", "25000"]]
    assert rows[1][3] != "NA" and rows[2][3] == "NA"  # window 2 is all N
    assert rows[2][4] == rows[1][4]  # cumulative carries over the NA window


def test_example_8_same_inputs_twice_give_byte_identical_files(tmp_path: Path) -> None:
    for d in ("a", "b"):
        _write("replicon_landmarks", landmarks(), tmp_path / d)
        _write("gc_skew", skew(), tmp_path / d)
    for name in ("replicon_landmarks.tsv", "gc_skew.tsv"):
        a = (tmp_path / "a" / "support" / name).read_bytes()
        assert a == (tmp_path / "b" / "support" / name).read_bytes()


def test_large_table_is_gzipped_with_mtime_zero_and_deterministic(tmp_path: Path) -> None:
    df = big_skew_frame()
    _write("gc_skew", df, tmp_path / "a")
    _write("gc_skew", df, tmp_path / "b")
    assert _support_files(tmp_path / "a") == ["gc_skew.tsv.gz"]
    raw = (tmp_path / "a" / "support" / "gc_skew.tsv.gz").read_bytes()
    assert raw[:2] == b"\x1f\x8b"
    assert raw[4:8] == b"\x00\x00\x00\x00"  # gzip MTIME field
    assert raw == (tmp_path / "b" / "support" / "gc_skew.tsv.gz").read_bytes()
    text = gzip.decompress(raw).decode("utf-8")
    assert len(text.encode("utf-8")) > 2**21
    lines = text.splitlines()
    assert lines[0] == "contig\tstart\tend\tgc_skew\tcumulative_gc_skew"
    assert len(lines) == len(df) + 1


def test_small_table_is_not_gzipped(tmp_path: Path) -> None:
    _write("gc_skew", big_skew_frame(rows=2_000), tmp_path)  # well under 10^6 bytes
    assert _support_files(tmp_path) == ["gc_skew.tsv"]


def test_rewrite_leaves_no_stale_twin(tmp_path: Path) -> None:
    _write("gc_skew", big_skew_frame(rows=2_000), tmp_path)
    _write("gc_skew", big_skew_frame(), tmp_path)
    assert _support_files(tmp_path) == ["gc_skew.tsv.gz"]
    _write("gc_skew", big_skew_frame(rows=2_000), tmp_path)
    assert _support_files(tmp_path) == ["gc_skew.tsv"]


def test_run_view_reads_back_what_was_written(tmp_path: Path) -> None:
    from gippa.io.run_view import RunView

    out = tmp_path / "OUT"
    (out / "tables").mkdir(parents=True)
    lm, gs = landmarks(), skew()
    _write("replicon_landmarks", lm, out / "plots")
    _write("gc_skew", gs, out / "plots")
    view = RunView.from_outdir(out)
    got_lm = view.support("replicon_landmarks")
    got_gs = view.support("gc_skew")
    assert got_lm is not None and got_gs is not None
    assert list(got_lm["contig"]) == ["chr1", "plasmid_1", "lin"]
    assert list(got_lm["status"]) == ["called", "no_oric", "linear"]
    assert [int(x) for x in got_gs["start"]] == [int(x) for x in gs["start"]]
    assert np.allclose(
        got_gs["gc_skew"].astype(float), gs["gc_skew"].astype(float), equal_nan=True, rtol=1e-5
    )


def test_big_table_read_back_through_run_view(tmp_path: Path) -> None:
    from gippa.io.run_view import RunView

    out = tmp_path / "OUT"
    (out / "tables").mkdir(parents=True)
    _write("gc_skew", big_skew_frame(), out / "plots")
    got = RunView.from_outdir(out).support("gc_skew")
    assert got is not None and len(got) == 60_000


# ---------------------------------------------------------------- schemas and descriptor


def _schema(name: str) -> dict[str, Any]:
    path = SCHEMAS / "plot_support" / f"{name}.json"
    assert path.is_file(), path
    data: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return data


@pytest.mark.parametrize("name", NAMES)
def test_schema_fields_follow_the_spec_columns(name: str) -> None:
    cols = {
        "replicon_landmarks": [
            "contig", "length_bp", "topology", "oric_position", "ter_position", "status",
        ],
        "gc_skew": ["contig", "start", "end", "gc_skew", "cumulative_gc_skew"],
    }[name]  # fmt: skip
    schema = _schema(name)
    assert [f["name"] for f in schema["fields"]] == cols
    assert "NA" in schema.get("missingValues", [])
    # not a result table: no §9.2 test columns
    assert not {"p_value", "q_value", "stat", "effect", "null_model"} & set(cols)


def test_landmark_schema_enums() -> None:
    fields = {f["name"]: f for f in _schema("replicon_landmarks")["fields"]}
    assert set(fields["topology"]["constraints"]["enum"]) == {"circular", "linear"}
    assert set(fields["status"]["constraints"]["enum"]) == {"called", "no_oric", "linear"}
    for col in ("length_bp", "oric_position", "ter_position"):
        assert fields[col]["type"] == "integer"


def test_gc_skew_schema_types() -> None:
    fields = {f["name"]: f for f in _schema("gc_skew")["fields"]}
    assert fields["start"]["type"] == fields["end"]["type"] == "integer"
    assert fields["gc_skew"]["type"] == fields["cumulative_gc_skew"]["type"] == "number"


@pytest.mark.parametrize("name", NAMES)
def test_datapackage_lists_each_support_table_with_role_plot_support(name: str) -> None:
    pkg = json.loads((SCHEMAS / "datapackage.json").read_text(encoding="utf-8"))
    res = [r for r in pkg["resources"] if r["name"] == name]
    assert len(res) == 1, [r["name"] for r in pkg["resources"]]
    r = res[0]
    assert r["role"] == "plot_support"
    assert r["path"].startswith(f"plots/support/{name}.tsv")
    assert r.get("profile") == "tabular-data-resource"
    assert (SCHEMAS / r["schema"]).resolve() == (SCHEMAS / "plot_support" / f"{name}.json")


def _validate(path: Path, name: str) -> None:
    frictionless = pytest.importorskip("frictionless")
    resource = frictionless.Resource(
        path=path.name,
        basepath=str(path.parent),
        format="tsv",
        schema=frictionless.Schema.from_descriptor(_schema(name)),
    )
    report = resource.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_written_landmarks_validate_against_their_schema(tmp_path: Path) -> None:
    _write("replicon_landmarks", landmarks(), tmp_path)
    _validate(tmp_path / "support" / "replicon_landmarks.tsv", "replicon_landmarks")


def test_written_gc_skew_validates_against_its_schema(tmp_path: Path) -> None:
    _write("gc_skew", skew(), tmp_path)
    _validate(tmp_path / "support" / "gc_skew.tsv", "gc_skew")


def test_gzipped_gc_skew_validates_against_its_schema(tmp_path: Path) -> None:
    _write("gc_skew", big_skew_frame(rows=70_000), tmp_path)
    _validate(tmp_path / "support" / "gc_skew.tsv.gz", "gc_skew")


def test_landmark_schema_rejects_an_unknown_status(tmp_path: Path) -> None:
    frictionless = pytest.importorskip("frictionless")
    bad = tmp_path / "bad.tsv"
    bad.write_text(
        "contig\tlength_bp\ttopology\toric_position\tter_position\tstatus\n"
        "chr1\t100\tcircular\tNA\tNA\tmaybe\n",
        encoding="utf-8",
    )
    resource = frictionless.Resource(
        path=bad.name, basepath=str(tmp_path), format="tsv",
        schema=frictionless.Schema.from_descriptor(_schema("replicon_landmarks")),
    )  # fmt: skip
    assert not resource.validate().valid


@pytest.mark.parametrize(("size", "gzipped"), [(2**20, False), (2**20 + 1, True)])
def test_gzip_threshold_is_two_to_the_twenty_bytes(
    tmp_path: Path, size: int, gzipped: bool
) -> None:
    # one column "x": header "x\n" plus one cell and its newline (lead decision 6)
    df = pd.DataFrame({"x": ["a" * (size - 3)]})
    path = Path(_write("gc_skew", df, tmp_path))
    assert path.name.endswith(".tsv.gz") is gzipped
    raw = gzip.decompress(path.read_bytes()) if gzipped else path.read_bytes()
    assert len(raw) == size
