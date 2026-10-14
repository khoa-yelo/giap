"""Lead decisions after the M1-12 code review (docs/specs/M1-12-preflight-cli.md, last section),
plus tests that kill the mutations which survived the review.

- Atomic outputs: a rerun removes a stale ``work/preflight/``; ``preflight/preflight.json`` is
  written last (completion marker), so after a failed run it reflects the failure and no
  normalised file of an earlier run is left behind.
- Resources are probed live by ``giap.preflight.environment.probe_resources(outdir)`` →
  ``(free_disk_bytes, ram_bytes)`` (``shutil.disk_usage`` and the Slurm/cgroup limit, never
  ``.giap/capabilities.json``); ``preflight.json`` keeps only the verdict, so it is byte-identical
  across machines, and the raw numbers go to ``OUT/logs/events.jsonl`` (§6.6, §6.7).
- Protein headers are safe IDs; ``work/preflight/protein_map.tsv`` maps them to the originals.
- ``--outdir`` naming an existing file exits 1.
- ``db_manifest.json`` entries ``{"version": ..., "sha256": ...}`` are accepted.

Every test isolates ``GIAP_REPO`` and ``GIAP_DB_DIR`` in ``tmp_path``.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import jsonschema
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner, Result

from giap.exitcodes import ExitCode
from giap.io.protein_fasta import read_protein_fasta

ROOT = Path(__file__).resolve().parents[3]
SCHEMA_PATH = ROOT / "src" / "giap" / "schemas" / "preflight.schema.json"
GIB = 1024**3
LOW = ["--set", "preflight.min_n50_warn=1", "--set", "preflight.min_largest_contig_warn=1"]
GENOME_MODE = ["--set", "annotation.promoters.mode=genome"]
LONG = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"
NORM_SUFFIXES = {".fa", ".faa", ".tsv", ".parquet"}


# ---------------------------------------------------------------- isolation and helpers


@pytest.fixture(autouse=True)
def _isolated_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Empty ``$GIAP_DB_DIR``; ``$GIAP_REPO`` holds a capabilities.json that must be ignored."""
    db = tmp_path / "_db"
    db.mkdir()
    repo = tmp_path / "_repo"
    (repo / ".giap").mkdir(parents=True)
    (repo / ".giap" / "capabilities.json").write_text(
        json.dumps({"allocation": {"mem_gb": 1, "node_mem_gb": 1}, "giap_root_free_gb": 0.001}),
        encoding="utf-8",
    )
    monkeypatch.setenv("GIAP_DB_DIR", str(db))
    monkeypatch.setenv("GIAP_REPO", str(repo))
    return db


def patch_resources(monkeypatch: pytest.MonkeyPatch, disk_bytes: int, ram_bytes: int) -> None:
    """Pretend to be a machine with ``disk_bytes`` free and ``ram_bytes`` of memory."""
    import giap.preflight.environment as env_mod
    import giap.preflight.run as run_mod

    def fake(*_a: object, **_k: object) -> tuple[int, int]:
        return disk_bytes, ram_bytes

    monkeypatch.setattr(env_mod, "probe_resources", fake)  # the decision: probe lives here
    monkeypatch.setattr(run_mod, "probe_resources", fake, raising=False)


def acgt(n: int, seed: int = 0) -> str:
    rng = np.random.default_rng(seed)
    return "".join(np.array(list("ACGT"))[rng.integers(0, 4, n)])


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def genome(tmp_path: Path, header: str = "c circular=true", n: int = 120_000) -> Path:
    return write(tmp_path / "g.fa", f">{header}\n{acgt(n)}\n")


def proteins(tmp_path: Path) -> Path:
    return write(tmp_path / "proteins.faa", f">P1\n{LONG}\n")


def peptides(tmp_path: Path, intensity: bool = True) -> Path:
    if intensity:
        return write(
            tmp_path / "peptides.tsv",
            "peptide\tprotein_ids\trun\tintensity\nPEPTIDEK\tP1\tr1\t100\n",
        )
    return write(tmp_path / "peptides.tsv", "peptide\tprotein_ids\nPEPTIDEK\tP1\n")


def invoke(*args: str | Path) -> Result:
    from giap.cli import cli

    return CliRunner().invoke(cli, ["preflight", *(str(a) for a in args)])


def assert_code(result: Result, code: int) -> None:
    assert result.exit_code == code, (result.exit_code, result.output, result.exception)


def load_report(out: Path) -> dict[str, Any]:
    path = out / "preflight" / "preflight.json"
    assert path.is_file(), f"{path} not written"
    data = json.loads(path.read_text(encoding="utf-8"))
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema).validate(data)
    assert isinstance(data, dict)
    return data


def by_rule(report: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {c["rule_id"]: c for c in report["checks"]}


def numbers(obj: object) -> Iterator[float]:
    if isinstance(obj, bool):
        return
    if isinstance(obj, int | float):
        yield float(obj)
    elif isinstance(obj, dict):
        for v in obj.values():
            yield from numbers(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from numbers(v)


def normalised_files(out: Path) -> list[Path]:
    return [p for p in out.rglob("*") if p.is_file() and p.suffix in NORM_SUFFIXES]


# ---------------------------------------------------------------- atomic outputs


def test_failing_rerun_leaves_no_work_files_and_report_reflects_failure(tmp_path: Path) -> None:
    out = tmp_path / "out"
    g, prot, pep = genome(tmp_path), proteins(tmp_path), peptides(tmp_path)
    r = invoke("--genome", g, "--proteins", prot, "--peptides", pep, *LOW, "--outdir", out)
    assert_code(r, ExitCode.SUCCESS)
    assert (out / "work" / "preflight" / "genome.norm.fa").is_file()
    assert load_report(out)["summary"]["n_error"] == 0

    write(g, ">c\nACGU\n")
    r = invoke("--genome", g, "--proteins", prot, "--peptides", pep, *LOW, "--outdir", out)
    assert_code(r, ExitCode.INVALID_INPUT)
    assert not (out / "work" / "preflight").exists()
    assert normalised_files(out) == []
    report = load_report(out)
    assert report["summary"]["n_error"] >= 1
    assert by_rule(report)["genome.parse"]["level"] == "error"
    md = (out / "preflight" / "preflight.md").read_text(encoding="utf-8")
    assert "0 errors" not in md


def test_rerun_removes_stale_work_files_of_an_earlier_run(tmp_path: Path) -> None:
    out = tmp_path / "out"
    g, prot = genome(tmp_path), proteins(tmp_path)
    ab = write(tmp_path / "abundance.tsv", "protein_id\tvalue\nP1\t1000\n")
    r = invoke("--genome", g, "--proteins", prot, "--abundance", ab, *LOW, "--outdir", out)
    assert_code(r, ExitCode.SUCCESS)
    assert (out / "work" / "preflight" / "abundance.norm.parquet").is_file()

    r = invoke("--genome", g, "--proteins", prot, *LOW, "--outdir", out)
    assert_code(r, ExitCode.SUCCESS)
    assert not (out / "work" / "preflight" / "abundance.norm.parquet").exists()
    assert (out / "work" / "preflight" / "genome.norm.fa").is_file()


def test_preflight_json_is_written_last(tmp_path: Path) -> None:
    out = tmp_path / "out"
    r = invoke(
        "--genome", genome(tmp_path), "--proteins", proteins(tmp_path),
        "--peptides", peptides(tmp_path), *LOW, "--outdir", out,
    )  # fmt: skip
    assert_code(r, ExitCode.SUCCESS)
    marker = (out / "preflight" / "preflight.json").stat().st_mtime_ns
    others = [
        p for p in out.rglob("*")
        if p.is_file() and p.name != "preflight.json" and "logs" not in p.relative_to(out).parts
    ]  # fmt: skip
    assert others
    for p in others:
        assert p.stat().st_mtime_ns <= marker, p


def test_crash_while_writing_leaves_no_marker_and_no_normalised_files(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def boom(*_a: object, **_k: object) -> None:
        raise RuntimeError("disk full")

    monkeypatch.setattr(pd.DataFrame, "to_parquet", boom)
    out = tmp_path / "out"
    r = invoke(
        "--genome", genome(tmp_path), "--proteins", proteins(tmp_path),
        "--peptides", peptides(tmp_path), *LOW, "--outdir", out,
    )  # fmt: skip
    assert r.exit_code != ExitCode.SUCCESS
    assert not (out / "preflight" / "preflight.json").exists()
    assert not (out / "work" / "preflight").exists()
    assert normalised_files(out) == []


# ---------------------------------------------------------------- live resources probe


@pytest.mark.parametrize("mode", [[], GENOME_MODE], ids=["targeted", "genome"])
def test_preflight_json_identical_on_two_machines(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: list[str]
) -> None:
    g, prot = genome(tmp_path), proteins(tmp_path)
    machines = {"a": (123 * GIB, 37 * GIB), "b": (987 * GIB, 251 * GIB)}
    for name, (disk, ram) in machines.items():
        patch_resources(monkeypatch, disk, ram)
        r = invoke("--genome", g, "--proteins", prot, *LOW, *mode, "--outdir", tmp_path / name)
        assert_code(r, ExitCode.SUCCESS)
    a = (tmp_path / "a" / "preflight" / "preflight.json").read_bytes()
    b = (tmp_path / "b" / "preflight" / "preflight.json").read_bytes()
    assert a == b
    res = by_rule(load_report(tmp_path / "a"))["environment.resources"]
    raw = {37.0, 123.0, float(37 * GIB), float(123 * GIB)}
    assert raw.isdisjoint(numbers(res["value"])), res


def test_raw_resource_numbers_go_to_events_jsonl(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    disk, ram = 123 * GIB, 37 * GIB
    patch_resources(monkeypatch, disk, ram)
    out = tmp_path / "out"
    r = invoke("--genome", genome(tmp_path), "--proteins", proteins(tmp_path), *LOW,
               "--outdir", out)  # fmt: skip
    assert_code(r, ExitCode.SUCCESS)
    path = out / "logs" / "events.jsonl"
    assert path.is_file()
    events = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line]
    assert events
    for e in events:
        assert {"ts", "step", "level", "event", "data"} <= set(e), e
    found = {n for e in events if e["step"] == "preflight" for n in numbers(e["data"])}
    assert found & {float(ram), 37.0, ram / 1e9}, events
    assert found & {float(disk), 123.0, disk / 1e9}, events


def test_resources_come_from_the_probe_not_capabilities_json(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """capabilities.json (fixture) says 1 GB RAM; the live probe says 64 GiB → no RAM warning."""
    patch_resources(monkeypatch, 500 * GIB, 64 * GIB)
    out = tmp_path / "out"
    r = invoke("--genome", genome(tmp_path), "--proteins", proteins(tmp_path), *LOW,
               *GENOME_MODE, "--outdir", out)  # fmt: skip
    assert_code(r, ExitCode.SUCCESS)
    assert by_rule(load_report(out))["environment.resources"]["level"] == "ok"

    patch_resources(monkeypatch, 500 * GIB, 8 * GIB)
    out2 = tmp_path / "out2"
    r = invoke("--genome", genome(tmp_path), "--proteins", proteins(tmp_path), *LOW,
               *GENOME_MODE, "--outdir", out2, "--strict")  # fmt: skip
    assert_code(r, ExitCode.WARNINGS_STRICT)
    assert by_rule(load_report(out2))["environment.resources"]["level"] == "warn"


def test_probe_resources_uses_disk_usage_of_outdir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import shutil

    from giap.preflight.environment import probe_resources

    seen: list[Path] = []

    def fake_usage(path: object) -> Any:
        seen.append(Path(str(path)))
        return shutil._ntuple_diskusage(1000 * GIB, 1, 77 * GIB)  # type: ignore[attr-defined]

    monkeypatch.setattr(shutil, "disk_usage", fake_usage)
    disk, ram = probe_resources(tmp_path)
    assert disk == 77 * GIB
    assert seen and all(p == tmp_path or tmp_path.is_relative_to(p) for p in seen)
    assert isinstance(ram, int) and ram > 0


# ---------------------------------------------------------------- safe protein IDs


def test_protein_headers_are_safe_ids_with_protein_map(tmp_path: Path) -> None:
    prot = write(
        tmp_path / "proteins.faa",
        f">A;B group of two\n{LONG}\n>sp|P12345|ABC_ECOLI\n{LONG[::-1]}\n>Q:9\n{LONG[1:]}\n",
    )
    out = tmp_path / "out"
    r = invoke("--genome", genome(tmp_path), "--proteins", prot, *LOW, "--outdir", out)
    assert_code(r, ExitCode.SUCCESS)
    work = out / "work" / "preflight"
    text = (work / "proteins.long.faa").read_text(encoding="utf-8")
    headers = [line[1:] for line in text.splitlines() if line.startswith(">")]
    assert headers == ["A_B", "P12345", "Q_9"]
    assert [p.id for p in read_protein_fasta(work / "proteins.long.faa")] == headers

    with (work / "protein_map.tsv").open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    assert "original_id" in rows[0]

    def row_for(safe: str) -> dict[str, str]:
        hits = [row for row in rows if safe in row.values()]
        assert len(hits) == 1, (safe, rows)
        return hits[0]

    assert row_for("A_B")["original_id"] == "A;B"
    assert {"A", "B"} <= set(";".join(row_for("A_B").values()).replace(",", ";").split(";"))
    assert row_for("P12345")["original_id"] == "P12345"
    assert row_for("Q_9")["original_id"] == "Q:9"


# ---------------------------------------------------------------- --outdir is a file


def test_outdir_that_is_a_file_exits_1(tmp_path: Path) -> None:
    out = write(tmp_path / "out", "not a directory\n")
    r = invoke("--genome", genome(tmp_path), "--proteins", proteins(tmp_path), *LOW,
               "--outdir", out)  # fmt: skip
    assert_code(r, ExitCode.INVALID_INPUT)
    assert out.read_text(encoding="utf-8") == "not a directory\n"


# ---------------------------------------------------------------- provisional db_manifest


def test_provisional_db_manifest_format_is_accepted(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, _isolated_env: Path
) -> None:
    from giap.config import resolve_config
    from giap.preflight.environment import requirements
    from giap.preflight.run import PreflightInputs, exit_code, run_preflight

    cfg = resolve_config(
        overrides=["preflight.min_n50_warn=1", "preflight.min_largest_contig_warn=1"]
    )
    manifest = {
        r.name: {"version": r.min_version or "1.0", "sha256": "ab" * 32, "path": f"{r.name}/"}
        for r in requirements(cfg.modules, cfg.annotation.promoters.mode)
    }
    manifest["bakta_db"]["version"] = "6.0"
    (_isolated_env / "db_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    patch_resources(monkeypatch, 500 * GIB, 64 * GIB)
    inp = PreflightInputs(genome=genome(tmp_path), proteins=proteins(tmp_path))
    report, _ = run_preflight(inp, cfg, None, "run")
    tools = [c for c in report.checks if c.rule_id == "environment.tools"]
    assert [c.level for c in tools] == ["ok"], tools
    assert exit_code(report, strict=False) == ExitCode.SUCCESS


# ---------------------------------------------------------------- surviving mutations


@pytest.mark.parametrize("which", ["genome", "proteins", "peptides"])
def test_parse_error_messages_hold_file_names_not_paths(tmp_path: Path, which: str) -> None:
    """``_scrub``: decoding errors name the file; ``preflight.json`` must not hold its directory."""
    inputs = tmp_path / "deep" / "input_dir"
    inputs.mkdir(parents=True)
    g, prot = genome(inputs), proteins(inputs)
    args: list[str | Path] = ["--genome", g, "--proteins", prot]
    if which == "genome":
        g.write_bytes(b">c\nAC\xff\xfeGT\n")
    elif which == "proteins":
        prot.write_bytes(b">P1\nMK\xff\xfeTAY\n")
    else:
        pep = inputs / "peptides.tsv"
        pep.write_bytes(b"")
        args += ["--peptides", pep]
    out = tmp_path / "out"
    r = invoke(*args, *LOW, "--outdir", out)
    assert_code(r, ExitCode.INVALID_INPUT)
    text = (out / "preflight" / "preflight.json").read_text(encoding="utf-8")
    assert str(inputs) not in text
    report = load_report(out)
    errors = [c for c in report["checks"] if c["level"] == "error" and c["input"] == which]
    assert errors, report["checks"]


def test_abundance_error_message_holds_file_name_not_path(tmp_path: Path) -> None:
    inputs = tmp_path / "deep" / "input_dir"
    inputs.mkdir(parents=True)
    ab = inputs / "abundance.tsv"
    ab.write_bytes(b"")  # read_tsv: "<path>: empty TSV (no header)"
    out = tmp_path / "out"
    r = invoke("--genome", genome(tmp_path), "--proteins", proteins(tmp_path),
               "--abundance", ab, *LOW, "--outdir", out)  # fmt: skip
    assert_code(r, ExitCode.INVALID_INPUT)
    shape = by_rule(load_report(out))["abundance.shape"]
    assert shape["level"] == "error"
    assert "abundance.tsv" in shape["message"]
    text = (out / "preflight" / "preflight.json").read_text(encoding="utf-8")
    assert str(inputs) not in text


@pytest.mark.parametrize(("intensity", "used"), [(True, False), (False, True)])
def test_abundance_used_only_without_peptide_intensities(
    tmp_path: Path, intensity: bool, used: bool
) -> None:
    ab = write(tmp_path / "abundance.tsv", "protein_id\tvalue\nP1\t1000\n")
    out = tmp_path / "out"
    r = invoke("--genome", genome(tmp_path), "--proteins", proteins(tmp_path),
               "--peptides", peptides(tmp_path, intensity), "--abundance", ab, *LOW,
               "--outdir", out)  # fmt: skip
    assert_code(r, ExitCode.SUCCESS)
    assert load_report(out)["routing"]["abundance_used"] is used


def test_abundance_join_uses_the_protein_list(tmp_path: Path) -> None:
    ab = write(tmp_path / "abundance.tsv", "protein_id\tvalue\nP1\t1000\nX9\t5\n")
    out = tmp_path / "out"
    r = invoke("--genome", genome(tmp_path), "--proteins", proteins(tmp_path),
               "--abundance", ab, *LOW, "--outdir", out)  # fmt: skip
    assert_code(r, ExitCode.SUCCESS)
    join = by_rule(load_report(out))["abundance.join"]
    assert join["level"] == "warn", join
    assert "X9" in json.dumps(join["value"])


@pytest.mark.parametrize(
    ("length", "kind"), [(1_000_000, "chromosome"), (999_999, "contig")], ids=["1Mb", "below"]
)
def test_linear_contig_of_1_mb_becomes_the_chromosome(
    tmp_path: Path, length: int, kind: str
) -> None:
    out = tmp_path / "out"
    r = invoke("--genome", genome(tmp_path, "big", length), "--proteins", proteins(tmp_path),
               "--outdir", out)  # fmt: skip
    assert_code(r, ExitCode.SUCCESS)
    assert (out / "work" / "preflight" / "replicons.tsv").read_text(encoding="utf-8") == (
        f"big\tbig\t{kind}\tl\t-\n"
    )


def test_chromosome_min_bp_follows_min_largest_contig_warn(tmp_path: Path) -> None:
    out = tmp_path / "out"
    r = invoke("--genome", genome(tmp_path, "big", 300_000), "--proteins", proteins(tmp_path),
               "--set", "preflight.min_largest_contig_warn=200000", "--outdir", out)  # fmt: skip
    assert_code(r, ExitCode.SUCCESS)
    assert (out / "work" / "preflight" / "replicons.tsv").read_text(encoding="utf-8") == (
        "big\tbig\tchromosome\tl\t-\n"
    )
