# ruff: noqa: F811  (pytest fixtures imported from test_cli_run_wiring)
"""``giap run --annotation bakta.json.gz`` is the precomputed Bakta layer, like ``bakta.json``
(GIAP_SPEC.md §5.3 "Bakta JSON: used as is, and Bakta is not run"; backlog fix
fix/annotation-json-gz, lead decisions 1 and 3).

Before the fix ``cli._bakta_json`` accepted only the ``.json`` suffix, so a gzipped Bakta JSON
left ``annotation.precomputed.bakta_json`` null and ``workflow/subworkflows/local/annotate.nf``
ran live Bakta (its ``else`` branch). Lead decisions:

1. ``--annotation x.json.gz`` is accepted exactly like ``x.json`` (any letter case): the
   precomputed path is used and Bakta does not run. The CLI may pass the ``.gz`` itself or a
   copy decompressed once into the run's work area (under ``--outdir``); either way the file
   named by ``annotation.precomputed.bakta_json`` must parse to the same annotation.
2. Withdrawn by the lead: GFF3 and GenBank annotations stay accepted as GIAP_SPEC.md §5.1 says
   (their ``--regions`` wiring is a separate backlog item); only ``.json.gz`` changes here.
3. ``run_manifest.json`` records the sha256 of the file as given (the compressed bytes).

A fake ``nextflow`` records its argv; nothing real is launched and no network or real database
is touched (``$GIAP_ROOT`` and ``$GIAP_DB_DIR`` are fakes from ``test_cli_run_wiring``).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import tempfile
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from unit.test_cli_run_wiring import (  # noqa: F401  (fixtures)
    _fake_giap_root,
    _run,
    db,
    fake_nf,
    genome,
    proteins,
)

from giap.annotate.bakta import parse_bakta_json
from giap.exitcodes import ExitCode

FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "recorded"
    / "makedb"
    / "sim100k_plasmid_1.bakta.json"
)


def _gzip_fixture(dest: Path) -> Path:
    """The committed Bakta JSON fixture, gzipped at test time to ``dest``."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    with FIXTURE.open("rb") as src, gzip.open(dest, "wb") as out:
        shutil.copyfileobj(src, out)
    return dest


def _params(out: Path) -> dict[str, object]:
    path = out / "work" / "params.json"
    assert path.is_file(), "giap run wrote no params.json"
    data = json.loads(path.read_text())
    assert isinstance(data, dict)
    return data


def _precomputed_bakta(out: Path) -> str | None:
    ann = _params(out)["annotation"]
    assert isinstance(ann, dict)
    value = ann["precomputed"]["bakta_json"]
    assert value is None or isinstance(value, str)
    return value


def _same_annotation(path: Path) -> None:
    """``path`` (plain or gzipped) parses to the fixture's annotation."""
    got, want = parse_bakta_json(path), parse_bakta_json(FIXTURE)
    assert got.model == want.model
    assert dict(got.sequences) == dict(want.sequences)
    assert (got.bakta_version, got.db_version, got.db_type) == (
        want.bakta_version,
        want.db_version,
        want.db_type,
    )


def _assert_precomputed_from(gz: Path, out: Path, bakta: str | None) -> None:
    """Decision 1: the precomputed layer is the given ``.gz`` or its copy under ``out``."""
    assert bakta is not None, (
        f"--annotation {gz.name} left annotation.precomputed.bakta_json null: "
        "annotate.nf would run live Bakta"
    )
    path = Path(bakta)
    assert path.is_absolute(), bakta
    assert path.is_file(), f"precomputed bakta_json {bakta} does not exist"
    if path != gz.absolute():
        assert path.is_relative_to(out.absolute()), (
            f"decompressed copy {bakta} is not in the run's work area under {out}"
        )
    _same_annotation(path)


# ---------------------------------------------------------------- the loader


def test_loader_reads_gzipped_fixture_like_plain(tmp_path: Path) -> None:
    gz = _gzip_fixture(tmp_path / "sim100k_plasmid_1.bakta.json.gz")
    assert gz.read_bytes()[:2] == b"\x1f\x8b"
    _same_annotation(gz)


# ---------------------------------------------------------------- nextflow executor params


@pytest.mark.parametrize("name", ["bakta.json.gz", "BAKTA.JSON.GZ", "sample.v2.Json.Gz"])
def test_annotation_json_gz_becomes_precomputed_bakta_json(
    name: str, proteins: Path, genome: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    gz = _gzip_fixture(tmp_path / "in" / name)
    out = tmp_path / "OUT"
    code, text = _run(
        ["--genome", str(genome), "--annotation", str(gz), "--outdir", str(out)], proteins
    )
    assert code == ExitCode.SUCCESS, text
    _assert_precomputed_from(gz, out, _precomputed_bakta(out))


def test_annotation_json_gz_never_yields_live_bakta_params(
    proteins: Path, genome: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    """annotate.nf runs BAKTA exactly when ``precomputed.bakta_json`` is falsy; a ``.gz`` must
    never leave it so, nor stage a Bakta ``--regions`` file."""
    gz = _gzip_fixture(tmp_path / "in" / "bakta.json.gz")
    out = tmp_path / "OUT"
    code, text = _run(
        ["--genome", str(genome), "--annotation", str(gz), "--outdir", str(out)], proteins
    )
    assert code == ExitCode.SUCCESS, text
    assert fake_nf.is_file(), "nextflow was not launched"
    params = _params(out)
    assert params["annotation"]["precomputed"]["bakta_json"], "live Bakta params"  # type: ignore[index]
    inputs = params["inputs"]
    assert isinstance(inputs, dict)
    assert "regions_gff3" not in inputs, "a Bakta --regions file was staged for a Bakta JSON"


def test_set_still_wins_over_a_json_gz_annotation(
    proteins: Path, genome: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    gz = _gzip_fixture(tmp_path / "in" / "bakta.json.gz")
    other = tmp_path / "in" / "other.json"
    shutil.copyfile(FIXTURE, other)
    out = tmp_path / "OUT"
    argv = ["--genome", str(genome), "--annotation", str(gz), "--outdir", str(out),
            "--set", f"annotation.precomputed.bakta_json={other}"]  # fmt: skip
    code, text = _run(argv, proteins)
    assert code == ExitCode.SUCCESS, text
    assert _precomputed_bakta(out) == str(other.absolute())


@settings(
    max_examples=6,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
@given(
    stem=st.from_regex(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,12}", fullmatch=True),
    suffix=st.sampled_from([".json", ".json.gz"]).flatmap(
        lambda s: st.tuples(*(st.sampled_from([c.lower(), c.upper()]) for c in s)).map("".join)
    ),
)
def test_property_bakta_json_suffix_any_case_is_precomputed(
    stem: str, suffix: str, proteins: Path, genome: Path, db: Path, fake_nf: Path
) -> None:
    """For every stem and every letter case of ``.json`` / ``.json.gz``, the run is launched
    with a precomputed Bakta JSON that parses to the given annotation."""
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        ann = root / "in" / f"{stem}{suffix}"
        if suffix.lower().endswith(".gz"):
            _gzip_fixture(ann)
        else:
            ann.parent.mkdir(parents=True)
            shutil.copyfile(FIXTURE, ann)
        out = root / "OUT"
        code, text = _run(
            ["--genome", str(genome), "--annotation", str(ann), "--outdir", str(out)], proteins
        )
        assert code == ExitCode.SUCCESS, text
        _assert_precomputed_from(ann, out, _precomputed_bakta(out))


def test_local_roots_stage_a_json_gz_annotation_as_bakta_json(
    proteins: Path, genome: Path, tmp_path: Path
) -> None:
    """``--executor local`` stages a ``.json.gz`` as ``bakta_json`` (before the fix it exited 1:
    "annotation.precomputed.bakta_json is required ... not a Bakta JSON")."""
    from giap.cli import _local_roots
    from giap.config import GiapConfig

    gz = _gzip_fixture(tmp_path / "in" / "bakta.json.gz")
    out = tmp_path / "OUT"
    inputs: dict[str, Path | None] = {"proteins": proteins, "peptides": None, "abundance": None,
                                      "annotation": gz, "reference_peptides": None,
                                      "db_dir": None}  # fmt: skip
    try:
        roots = _local_roots(GiapConfig(), out, genome, inputs)
    except SystemExit as exc:  # pragma: no cover - the pre-fix failure mode
        pytest.fail(f"_local_roots exited {exc.code} for --annotation {gz.name}")
    assert "bakta_json" in roots, sorted(roots)
    _assert_precomputed_from(gz, out, str(Path(roots["bakta_json"]).absolute()))


# ---------------------------------------------------------------- run manifest


def test_run_manifest_records_sha256_of_the_gz_as_given(
    proteins: Path, genome: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    gz = _gzip_fixture(tmp_path / "in" / "bakta.json.gz")
    out = tmp_path / "OUT"
    _run(["--genome", str(genome), "--annotation", str(gz), "--outdir", str(out)], proteins)
    _assert_precomputed_from(gz, out, _precomputed_bakta(out))
    manifest = out / "run_manifest.json"
    assert manifest.is_file(), "giap run wrote no run_manifest.json"
    m = json.loads(manifest.read_text())
    by = {r["role"]: r for r in m["inputs"]}
    assert "annotation" in by, sorted(by)
    rec = by["annotation"]
    assert rec["file"] == gz.name
    assert rec["sha256"] == hashlib.sha256(gz.read_bytes()).hexdigest(), (
        "manifest annotation sha256 is not that of the compressed file as given"
    )
