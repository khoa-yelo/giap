"""``giap compare OUT_A OUT_B [OUT_C ...] --outdir CMP [--labels ...]`` (M9-22 examples 5 and 6,
every exit-1 row of "Degenerate cases", ADR-M9-B items 1-4, lead decision 1; GIAP_SPEC.md §5.1,
§6.6, §7.9.3, §9.5).

Interface assumed beyond the spec text: tables are written to ``CMP/tables/compare_{detection,
abundance,stability}.tsv`` (ADR-M9-B item 1); error messages go to the combined CliRunner output;
no table is written when the command exits 1; the bootstrap draws from
``giap.rng.for_step("compare", cfg.seed)`` with ``B = stats.bootstrap`` and the first run's
``tu_id`` as TUs. Runs are small ``OUT/`` trees from ``tests/helpers/compare_data.py``. Never
touches the real ``$GIAP_DB_DIR`` or the network.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner, Result
from helpers.compare_data import (
    ABUNDANCE_COLUMNS,
    DETECTION_COLUMNS,
    SHA_A,
    SHA_B,
    STABILITY_COLUMNS,
    gene,
    genes_frame,
    tus_of,
    write_run,
)

from giap.exitcodes import ExitCode

SCHEMAS = Path(__file__).resolve().parents[2] / "src" / "giap" / "schemas" / "tables"
TABLES = ("compare_detection", "compare_abundance", "compare_stability")
FAST = ["--set", "stats.bootstrap=200"]
N = 12


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIAP_DB_DIR", str(tmp_path / "db_empty"))
    monkeypatch.setenv("GIAP_ROOT", str(tmp_path / "db_empty"))


def _genes(shift: float, detected_upto: int = N, names: list[str] | None = None) -> list[dict]:
    names = names or [f"g{i:02d}" for i in range(N)]
    vals = np.sin(np.arange(len(names)) * 1.7) * 3 + np.arange(len(names)) * 0.4
    pairs = enumerate(zip(names, vals, strict=True))
    return [gene(n, start=1000 * i, detected=i < detected_upto, log2=float(v + shift * (i % 3)),
                 tu=f"TU_{i // 2}") for i, (n, v) in pairs]  # fmt: skip


def _two(tmp: Path, **kw: object) -> tuple[Path, Path]:
    a = write_run(tmp / "runA", _genes(0.0), **kw)  # type: ignore[arg-type]
    b = write_run(tmp / "runB", _genes(0.5, detected_upto=10), **kw)  # type: ignore[arg-type]
    return a, b


def _compare(*args: str) -> Result:
    from giap.cli import cli

    res = CliRunner().invoke(cli, ["compare", *args])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


def _read(cmp: Path, name: str) -> pd.DataFrame:
    return pd.read_csv(cmp / "tables" / f"{name}.tsv", sep="\t", na_values=["NA"],
                       keep_default_na=False)  # fmt: skip


def _no_tables(cmp: Path) -> bool:
    return not any((cmp / "tables" / f"{t}.tsv").exists() for t in TABLES)


# ----------------------------------------------------------------- happy path


def test_two_runs_write_three_tables(tmp_path: Path) -> None:
    a, b = _two(tmp_path)
    cmp = tmp_path / "CMP"
    res = _compare(str(a), str(b), "--outdir", str(cmp), "--labels", "A1,A2", *FAST)
    assert res.exit_code == ExitCode.SUCCESS, res.output
    det, ab, st = (_read(cmp, t) for t in TABLES)
    assert list(det.columns) == list(DETECTION_COLUMNS)
    assert list(ab.columns) == list(ABUNDANCE_COLUMNS)
    assert list(st.columns) == list(STABILITY_COLUMNS)
    assert set(det["run"].dropna()) == {"A1", "A2"}
    assert (det["section"] == "matrix").sum() == 2 * N
    assert len(ab) == 1 and (ab["run_a"].iloc[0], ab["run_b"].iloc[0]) == ("A1", "A2")
    assert int(ab["n_shared"].iloc[0]) == 10
    assert int(ab["n_boot"].iloc[0]) == 200
    assert len(st) == 5


def test_detection_sections_are_contiguous_blocks(tmp_path: Path) -> None:
    """Sorted by section first (enum or alphabetical block order is not pinned)."""
    a, b = _two(tmp_path)
    cmp = tmp_path / "CMP"
    assert _compare(str(a), str(b), "--outdir", str(cmp), *FAST).exit_code == 0
    det = _read(cmp, "compare_detection")
    sections = list(det["section"])
    assert sections in (sorted(sections), sorted(sections, reverse=True))
    inter = det[det["section"] == "intersection"]
    counts = dict(zip(inter["intersection"], inter["n_genes"], strict=True))
    if counts.get("none") == 0:
        del counts["none"]
    assert counts == {"runA&runB": 10, "runA": 2}


@pytest.mark.parametrize("name", TABLES)
def test_schema_exists_and_tables_validate(tmp_path: Path, name: str) -> None:
    frictionless = pytest.importorskip("frictionless")
    path = SCHEMAS / f"{name}.schema.json"
    assert path.is_file(), path
    descriptor = {k: v for k, v in json.loads(path.read_text()).items() if k != "coordinates"}
    a, b = _two(tmp_path)
    c = write_run(tmp_path / "runC", _genes(1.0, detected_upto=6))
    cmp = tmp_path / "CMP"
    assert _compare(str(a), str(b), str(c), "--outdir", str(cmp), *FAST).exit_code == 0
    schema = frictionless.Schema.from_descriptor(descriptor)
    resource = frictionless.Resource(
        path=f"{name}.tsv", basepath=str(cmp / "tables"), schema=schema
    )
    report = resource.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_schema_enums() -> None:
    def field(stem: str, name: str) -> dict:
        doc = json.loads((SCHEMAS / f"{stem}.schema.json").read_text())
        (f,) = [f for f in doc["fields"] if f["name"] == name]
        return f  # type: ignore[no-any-return]

    assert set(field("compare_detection", "section")["constraints"]["enum"]) == {
        "matrix", "intersection",
    }  # fmt: skip
    assert set(field("compare_stability", "family")["constraints"]["enum"]) == {
        "islands", "operon_flags", "enriched_categories", "novel_orfs", "alt_sites",
    }  # fmt: skip
    assert field("compare_abundance", "method")["constraints"]["enum"] == ["spearman_tu_bootstrap"]
    for stem in TABLES:
        doc = json.loads((SCHEMAS / f"{stem}.schema.json").read_text())
        assert "NA" in doc["missingValues"]


def test_three_runs_pairs_and_families(tmp_path: Path) -> None:
    a, b = _two(tmp_path)
    c = write_run(tmp_path / "runC", _genes(1.0, detected_upto=6))
    cmp = tmp_path / "CMP"
    res = _compare(str(a), str(b), str(c), "--outdir", str(cmp), "--labels", "A1,A2,A3", *FAST)
    assert res.exit_code == 0, res.output
    ab = _read(cmp, "compare_abundance")
    assert list(zip(ab["run_a"], ab["run_b"], strict=True)) == [("A1", "A2"), ("A1", "A3"),
                                                               ("A2", "A3")]  # fmt: skip
    assert len(_read(cmp, "compare_stability")) == 15


def test_default_labels_are_basenames_with_suffixes(tmp_path: Path) -> None:
    a = write_run(tmp_path / "x" / "out", _genes(0.0))
    b = write_run(tmp_path / "y" / "out", _genes(0.5))
    c = write_run(tmp_path / "z" / "out", _genes(0.2))
    cmp = tmp_path / "CMP"
    assert _compare(str(a), str(b), str(c), "--outdir", str(cmp), *FAST).exit_code == 0
    runs = list(dict.fromkeys(_read(cmp, "compare_detection")["run"].dropna()))
    assert runs == ["out", "out_2", "out_3"]


# ----------------------------------------------------------- genome identity


def test_example5_genome_hashes_differ_exit_1(tmp_path: Path) -> None:
    a = write_run(tmp_path / "runA", _genes(0.0), preflight_sha=SHA_A)
    b = write_run(tmp_path / "runB", _genes(0.5), preflight_sha=SHA_B)
    cmp = tmp_path / "CMP"
    res = _compare(str(a), str(b), "--outdir", str(cmp), "--labels", "A1,A2", *FAST)
    assert res.exit_code == ExitCode.INVALID_INPUT
    for token in ("A1", "A2", SHA_A, SHA_B):
        assert token in res.output, token
    assert _no_tables(cmp)


def test_manifest_hash_mismatch_exit_1(tmp_path: Path) -> None:
    a = write_run(tmp_path / "runA", _genes(0.0), preflight_sha=None, manifest_sha=SHA_A)
    b = write_run(tmp_path / "runB", _genes(0.5), preflight_sha=None, manifest_sha=SHA_B)
    res = _compare(str(a), str(b), "--outdir", str(tmp_path / "CMP"), *FAST)
    assert res.exit_code == ExitCode.INVALID_INPUT


def test_manifest_wins_over_preflight(tmp_path: Path) -> None:
    """Decision 1: ``run_manifest.json`` when present, ``preflight.json`` otherwise."""
    a = write_run(tmp_path / "runA", _genes(0.0), preflight_sha=SHA_A, manifest_sha=SHA_A)
    b = write_run(tmp_path / "runB", _genes(0.5), preflight_sha=SHA_B, manifest_sha=SHA_A)
    res = _compare(str(a), str(b), "--outdir", str(tmp_path / "CMP"), *FAST)
    assert res.exit_code == 0, res.output


def test_manifest_in_one_run_preflight_in_the_other(tmp_path: Path) -> None:
    a = write_run(tmp_path / "runA", _genes(0.0), preflight_sha=None, manifest_sha=SHA_A)
    b = write_run(tmp_path / "runB", _genes(0.5), preflight_sha=SHA_A)
    assert _compare(str(a), str(b), "--outdir", str(tmp_path / "C1"), *FAST).exit_code == 0
    c = write_run(tmp_path / "runC", _genes(0.5), preflight_sha=SHA_B)
    res = _compare(str(a), str(c), "--outdir", str(tmp_path / "C2"), *FAST)
    assert res.exit_code == ExitCode.INVALID_INPUT


def test_no_genome_hash_cannot_verify(tmp_path: Path) -> None:
    a = write_run(tmp_path / "runA", _genes(0.0), preflight_sha=None)  # preflight w/o genome row
    b = write_run(tmp_path / "runB", _genes(0.5))
    cmp = tmp_path / "CMP"
    res = _compare(str(a), str(b), "--outdir", str(cmp), *FAST)
    assert res.exit_code == ExitCode.INVALID_INPUT
    assert "cannot verify the genome" in res.output
    assert _no_tables(cmp)


def test_no_preflight_and_no_manifest_cannot_verify(tmp_path: Path) -> None:
    a, b = _two(tmp_path)
    (a / "preflight" / "preflight.json").unlink()
    res = _compare(str(a), str(b), "--outdir", str(tmp_path / "CMP"), *FAST)
    assert res.exit_code == ExitCode.INVALID_INPUT
    assert "cannot verify the genome" in res.output


# ------------------------------------------------------------ other exit 1


def test_one_run_exit_1(tmp_path: Path) -> None:
    a, _ = _two(tmp_path)
    cmp = tmp_path / "CMP"
    res = _compare(str(a), "--outdir", str(cmp), *FAST)
    assert res.exit_code == ExitCode.INVALID_INPUT
    assert _no_tables(cmp)


def test_locus_tag_sets_differ_exit_1(tmp_path: Path) -> None:
    a = write_run(tmp_path / "runA", _genes(0.0))
    names = [f"g{i:02d}" for i in range(N - 1)] + ["gXX"]
    vals = [*_genes(0.5, names=names), gene("gYY", start=99_000, log2=1.0)]
    b = write_run(tmp_path / "runB", vals)
    cmp = tmp_path / "CMP"
    res = _compare(str(a), str(b), "--outdir", str(cmp), *FAST)
    assert res.exit_code == ExitCode.INVALID_INPUT
    # counts in each run (12, 13) and the symmetric difference size (g11, gXX, gYY -> 3)
    for token in ("12", "13", "3"):
        assert token in res.output, token
    assert _no_tables(cmp)


@pytest.mark.parametrize("labels", ["A1", "A1,A2,A3", "A1,A1"])
def test_bad_labels_exit_1(tmp_path: Path, labels: str) -> None:
    a, b = _two(tmp_path)
    cmp = tmp_path / "CMP"
    res = _compare(str(a), str(b), "--outdir", str(cmp), "--labels", labels, *FAST)
    assert res.exit_code == ExitCode.INVALID_INPUT
    assert _no_tables(cmp)


def test_not_implemented_is_gone(tmp_path: Path) -> None:
    a, b = _two(tmp_path)
    res = _compare(str(a), str(b), "--outdir", str(tmp_path / "CMP"), *FAST)
    assert res.exit_code != ExitCode.INTERNAL, res.output


# ------------------------------------------------------------- determinism


def test_example6_same_seed_byte_identical(tmp_path: Path) -> None:
    a, b = _two(tmp_path)
    for name in ("C1", "C2"):
        res = _compare(str(a), str(b), "--outdir", str(tmp_path / name), "--labels", "A1,A2",
                       *FAST, "--set", "seed=99")  # fmt: skip
        assert res.exit_code == 0, res.output
    for t in TABLES:
        one = (tmp_path / "C1" / "tables" / f"{t}.tsv").read_bytes()
        two = (tmp_path / "C2" / "tables" / f"{t}.tsv").read_bytes()
        assert one == two, t


def test_bootstrap_uses_for_step_compare_and_stats_bootstrap(tmp_path: Path) -> None:
    from giap.compare.abundance import abundance_pairs
    from giap.rng import for_step

    a, b = _two(tmp_path)
    cmp = tmp_path / "CMP"
    res = _compare(str(a), str(b), "--outdir", str(cmp), "--labels", "A1,A2",
                   "--set", "stats.bootstrap=300", "--set", "seed=4242")  # fmt: skip
    assert res.exit_code == 0, res.output
    got = _read(cmp, "compare_abundance").iloc[0]
    fa, fb = genes_frame(_genes(0.0)), genes_frame(_genes(0.5, detected_upto=10))
    want = abundance_pairs({"A1": fa, "A2": fb}, tus_of(fa), 300, for_step("compare", 4242))
    want_row = want.iloc[0]
    assert int(got["n_boot"]) == 300
    for col in ("spearman_rho", "ci_low", "ci_high", "ba_bias_log2"):
        assert got[col] == pytest.approx(float(want_row[col]), rel=1e-5, abs=1e-6), col


def test_different_seed_changes_the_ci(tmp_path: Path) -> None:
    a = write_run(tmp_path / "runA", [gene(f"g{i:03d}", start=100 * i, log2=float(v),
                                           tu=f"TU_{i // 2}")
                                      for i, v in enumerate(np.linspace(0, 8, 60))])  # fmt: skip
    noise = np.random.default_rng(1).normal(scale=2.0, size=60)
    b = write_run(tmp_path / "runB", [gene(f"g{i:03d}", start=100 * i, log2=float(v + e),
                                           tu=f"TU_{i // 2}")
                                      for i, (v, e) in enumerate(zip(np.linspace(0, 8, 60), noise,
                                                                     strict=True))])  # fmt: skip
    rows = []
    for seed in ("1", "2"):
        cmp = tmp_path / f"C{seed}"
        assert _compare(str(a), str(b), "--outdir", str(cmp), *FAST, "--set",
                        f"seed={seed}").exit_code == 0  # fmt: skip
        rows.append(_read(cmp, "compare_abundance").iloc[0])
    assert rows[0]["spearman_rho"] == rows[1]["spearman_rho"]
    assert (rows[0]["ci_low"], rows[0]["ci_high"]) != (rows[1]["ci_low"], rows[1]["ci_high"])
