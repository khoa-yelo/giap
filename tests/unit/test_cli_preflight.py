"""``giap preflight`` end to end through click (issue M1-12; GIAP_SPEC.md §5.1, §5.3, §6.2,
§7.0 Outputs and Tests, §7.1.2; ADR-0005 item 3).

Inputs are tiny files in ``tmp_path`` (genomes from ``giapsim`` or random ACGT). ``GIAP_DB_DIR``
points at an empty directory, so ``db_manifest.json`` is missing: under ``giap preflight`` that is
info, never an environment error (M1-11 question 1), and the exit codes below come from the
input rules only. §1.1 makes evidence mandatory (``inputs.evidence_present``), so the spec's
genome-only examples 1-3 add a tiny valid ``--proteins`` file to reach their stated exit code.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

import jsonschema
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner, Result

from giap.exitcodes import ExitCode
from giap.io.fasta import read_genome_fasta
from giap.io.protein_fasta import read_protein_fasta
from giapsim import SimConfig, simulate

ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = ROOT / "src" / "giap" / "schemas" / "preflight.schema.json"
LOW = ["--set", "preflight.min_n50_warn=1", "--set", "preflight.min_largest_contig_warn=1"]
WORK_FILES = (
    "genome.norm.fa",
    "contig_map.tsv",
    "replicons.tsv",
    "peptides.norm.parquet",
    "peptide_runs.parquet",
)
LONG_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"


# ---------------------------------------------------------------- helpers


@pytest.fixture(autouse=True)
def _no_db_manifest(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    db = tmp_path / "db_empty"
    db.mkdir()
    monkeypatch.setenv("GIAP_DB_DIR", str(db))


def random_acgt(n: int, seed: int = 0) -> str:
    rng = np.random.default_rng(seed)
    return "".join(np.array(list("ACGT"))[rng.integers(0, 4, n)])


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def genome_120kb(tmp_path: Path, header: str = "c circular=true") -> Path:
    return write(tmp_path / "g.fa", f">{header}\n{random_acgt(120_000)}\n")


def proteins(tmp_path: Path) -> Path:
    return write(tmp_path / "proteins.faa", f">P1\n{LONG_PROTEIN}\n")


def peptides(tmp_path: Path) -> Path:
    return write(
        tmp_path / "peptides.tsv",
        "peptide\tprotein_ids\trun\tintensity\n"
        "PEPTIDEK\tP1\tr1\t100\n"
        "PEPTIDEK\tP1\tr2\t300\n"
        "AAGLIEVQAPILSR\tP1\tr1\t50\n",
    )


def invoke(*args: str | Path) -> Result:
    from giap.cli import cli

    return CliRunner().invoke(cli, ["preflight", *(str(a) for a in args)])


def schema() -> dict[str, Any]:
    data = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def load_report(out: Path) -> dict[str, Any]:
    """Read ``OUT/preflight/preflight.json`` and validate it against its schema (M1 gate)."""
    path = out / "preflight" / "preflight.json"
    assert path.is_file(), f"{path} not written"
    data = json.loads(path.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema()).validate(data)
    return data


def by_rule(report: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {c["rule_id"]: c for c in report["checks"]}


def tree(root: Path) -> list[str]:
    return sorted(str(p.relative_to(root)) for p in root.rglob("*"))


def assert_code(result: Result, code: int) -> None:
    assert result.exit_code == code, (result.exit_code, result.output, result.exception)


# ---------------------------------------------------------------- spec examples 1-5


def test_example_1_clean_genome_exits_0(tmp_path: Path) -> None:
    r = invoke("--genome", genome_120kb(tmp_path), "--proteins", proteins(tmp_path), *LOW)
    assert_code(r, ExitCode.SUCCESS)
    assert "0 errors" in r.stdout


def test_example_1_strict_with_only_info_and_ok_exits_0(tmp_path: Path) -> None:
    r = invoke(
        "--genome", genome_120kb(tmp_path), "--proteins", proteins(tmp_path), *LOW, "--strict"
    )
    assert_code(r, ExitCode.SUCCESS)


def test_example_2_contiguity_warning_exits_0(tmp_path: Path) -> None:
    out = tmp_path / "out"
    r = invoke(
        "--genome", genome_120kb(tmp_path), "--proteins", proteins(tmp_path), "--outdir", out
    )
    assert_code(r, ExitCode.SUCCESS)
    assert by_rule(load_report(out))["genome.contiguity"]["level"] == "warn"


def test_example_3_contiguity_warning_with_strict_exits_3(tmp_path: Path) -> None:
    r = invoke("--genome", genome_120kb(tmp_path), "--proteins", proteins(tmp_path), "--strict")
    assert_code(r, ExitCode.WARNINGS_STRICT)
    assert r.exit_code == 3


@pytest.mark.parametrize("strict", [False, True], ids=["plain", "strict"])
def test_examples_4_5_invalid_genome_exits_1_errors_win_over_strict(
    tmp_path: Path, strict: bool
) -> None:
    g = write(tmp_path / "g.fa", ">c\nACGU\n")
    r = invoke("--genome", g, "--proteins", proteins(tmp_path), *(["--strict"] if strict else []))
    assert_code(r, ExitCode.INVALID_INPUT)


def test_genome_only_lacks_evidence_and_exits_1(tmp_path: Path) -> None:
    out = tmp_path / "out"
    r = invoke("--genome", genome_120kb(tmp_path), *LOW, "--outdir", out)
    assert_code(r, ExitCode.INVALID_INPUT)
    assert by_rule(load_report(out))["inputs.evidence_present"]["level"] == "error"


def test_unparseable_genome_is_input_error_not_crash(tmp_path: Path) -> None:
    g = write(tmp_path / "g.fa", "ACGT\n>c\nACGT\n")  # sequence before the first header
    out = tmp_path / "out"
    r = invoke("--genome", g, "--proteins", proteins(tmp_path), "--outdir", out)
    assert_code(r, ExitCode.INVALID_INPUT)
    assert by_rule(load_report(out))["genome.parse"]["level"] == "error"


# ---------------------------------------------------------------- spec examples 6-8: files


def test_example_6_bad_peptide_table_writes_report_but_no_work_files(tmp_path: Path) -> None:
    pep = write(tmp_path / "p.tsv", "sequence\tq_value\nPEPTIDEK\t0.001\n")
    out = tmp_path / "out"
    r = invoke("--genome", genome_120kb(tmp_path), *LOW, "--peptides", pep, "--outdir", out)
    assert_code(r, ExitCode.INVALID_INPUT)
    report = load_report(out)
    assert by_rule(report)["peptides.schema"]["level"] == "error"
    assert report["summary"]["n_error"] >= 1
    assert (out / "preflight" / "preflight.md").is_file()
    assert not (out / "work" / "preflight").exists()


def test_example_7_valid_peptides_write_normalised_work_files(tmp_path: Path) -> None:
    out = tmp_path / "out"
    r = invoke(
        "--genome", genome_120kb(tmp_path), *LOW, "--peptides", peptides(tmp_path), "--outdir", out
    )
    assert_code(r, ExitCode.SUCCESS)
    work = out / "work" / "preflight"
    for name in WORK_FILES:
        assert (work / name).is_file(), name
    assert (out / "preflight" / "preflight.json").is_file()
    assert (out / "preflight" / "preflight.md").is_file()
    load_report(out)

    pep = pd.read_parquet(work / "peptides.norm.parquet")
    assert sorted(pep["peptide"]) == ["AAGLIEVQAPILSR", "PEPTIDEK"]
    row = pep.set_index("peptide").loc["PEPTIDEK"]
    assert list(row["protein_ids"]) == ["P1"]
    assert row["n_runs"] == 2
    runs = pd.read_parquet(work / "peptide_runs.parquet")
    assert len(runs) == 3
    assert set(runs["run"]) == {"r1", "r2"}


def test_example_8_unsafe_id_contig_map_and_bakta_replicon_table(tmp_path: Path) -> None:
    g = write(tmp_path / "g.fa", f">chr:1 circular=true\n{random_acgt(120_000).lower()}\n")
    out = tmp_path / "out"
    r = invoke("--genome", g, "--proteins", proteins(tmp_path), *LOW, "--outdir", out)
    assert_code(r, ExitCode.SUCCESS)
    work = out / "work" / "preflight"
    assert (work / "contig_map.tsv").read_text(encoding="utf-8") == (
        "original_id\tcontig\nchr:1\tchr_1\n"
    )
    assert (work / "replicons.tsv").read_text(
        encoding="utf-8"
    ) == "chr_1\tchr_1\tchromosome\tc\t-\n"
    norm = read_genome_fasta(work / "genome.norm.fa")
    assert [x.id for x in norm] == ["chr_1"]
    assert norm[0].sequence == random_acgt(120_000)
    text = (work / "genome.norm.fa").read_text(encoding="utf-8")
    assert "acgt" not in text.split("\n", 1)[1]
    report = load_report(out)
    assert report["contig_map"] == [{"id": "chr_1", "original_id": "chr:1"}]
    assert by_rule(report)["genome.ids_renamed"]["level"] == "warn"


def test_replicons_tsv_for_complete_genome_with_two_small_plasmids(tmp_path: Path) -> None:
    """§7.0 Tests: a complete genome with two small plasmids raises no contiguity warning."""
    g = write(
        tmp_path / "g.fa",
        f">chr [topology=circular]\n{random_acgt(1_200_000, 1)}\n"
        f">p1 circular=true\n{random_acgt(60_000, 2)}\n"
        f">p2 circular=true\n{random_acgt(5_000, 3)}\n"
        f">frag\n{random_acgt(3_000, 4)}\n",
    )
    out = tmp_path / "out"
    r = invoke("--genome", g, "--proteins", proteins(tmp_path), "--outdir", out)
    assert_code(r, ExitCode.SUCCESS)
    report = load_report(out)
    assert by_rule(report)["genome.contiguity"]["level"] == "ok"
    assert (out / "work" / "preflight" / "replicons.tsv").read_text(encoding="utf-8") == (
        "chr\tchr\tchromosome\tc\t-\n"
        "p1\tp1\tplasmid\tc\t-\n"
        "p2\tp2\tplasmid\tc\t-\n"
        "frag\tfrag\tcontig\tl\t-\n"
    )
    assert (out / "work" / "preflight" / "contig_map.tsv").read_text(encoding="utf-8") == (
        "original_id\tcontig\nchr\tchr\np1\tp1\np2\tp2\nfrag\tfrag\n"
    )
    assert [x["kind"] for x in report["replicons"]] == [
        "chromosome",
        "plasmid",
        "plasmid",
        "unknown",
    ]


# ---------------------------------------------------------------- example 9: determinism (R7)


def test_example_9_rerun_is_byte_identical(tmp_path: Path) -> None:
    g, pep, prot = genome_120kb(tmp_path), peptides(tmp_path), proteins(tmp_path)
    outs = [tmp_path / "out1", tmp_path / "out2"]
    for out in outs:
        r = invoke("--genome", g, "--peptides", pep, "--proteins", prot, *LOW, "--outdir", out)
        assert_code(r, ExitCode.SUCCESS)
    for rel in (
        "preflight/preflight.json",
        "work/preflight/contig_map.tsv",
        "work/preflight/replicons.tsv",
        "work/preflight/genome.norm.fa",
        "work/preflight/proteins.long.faa",
    ):
        assert (outs[0] / rel).read_bytes() == (outs[1] / rel).read_bytes(), rel
    text = (outs[0] / "preflight" / "preflight.json").read_text(encoding="utf-8")
    assert str(tmp_path) not in text


# ---------------------------------------------------------------- outputs and stdout


def test_without_outdir_nothing_is_written(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    g, prot = genome_120kb(inputs), proteins(inputs)
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    before = tree(tmp_path)
    r = invoke("--genome", g, "--proteins", prot, *LOW)
    assert_code(r, ExitCode.SUCCESS)
    assert tree(tmp_path) == before
    assert "0 errors" in r.stdout


def test_stdout_is_the_markdown_summary(tmp_path: Path) -> None:
    out = tmp_path / "out"
    r = invoke(
        "--genome", genome_120kb(tmp_path), "--proteins", proteins(tmp_path), "--outdir", out
    )
    assert_code(r, ExitCode.SUCCESS)
    md = (out / "preflight" / "preflight.md").read_text(encoding="utf-8")
    assert md.startswith("# GIAP preflight")
    assert md.strip() in r.stdout
    assert "0 errors, 1 warning" in md  # genome.contiguity only


def test_inputs_record_roles_names_hashes_and_sizes(tmp_path: Path) -> None:
    g, prot, pep = genome_120kb(tmp_path), proteins(tmp_path), peptides(tmp_path)
    out = tmp_path / "out"
    r = invoke("--genome", g, "--proteins", prot, "--peptides", pep, *LOW, "--outdir", out)
    assert_code(r, ExitCode.SUCCESS)
    inputs = {i["role"]: i for i in load_report(out)["inputs"]}
    for role, path in (("genome", g), ("proteins", prot), ("peptides", pep)):
        assert inputs[role]["file"] == path.name
        assert inputs[role]["size"] == path.stat().st_size
        assert inputs[role]["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()


def test_proteins_routed_long_file_and_short_ones_become_peptides(tmp_path: Path) -> None:
    prot = write(
        tmp_path / "proteins.faa",
        f">sp|P12345|ABC_ECOLI long one\n{LONG_PROTEIN}\n"
        f">P2\n{LONG_PROTEIN[::-1]}\n"
        ">S1 short\nMKVLAAGIVGR\n"  # 11 aa: peptide path (7-40)
        ">T1 tiny\nMKV\n",  # 3 aa: dropped (warn)
    )
    out = tmp_path / "out"
    r = invoke("--genome", genome_120kb(tmp_path), "--proteins", prot, *LOW, "--outdir", out)
    assert_code(r, ExitCode.SUCCESS)
    work = out / "work" / "preflight"
    long = read_protein_fasta(work / "proteins.long.faa")
    assert [(p.id, p.sequence) for p in long] == [
        ("P12345", LONG_PROTEIN),
        ("P2", LONG_PROTEIN[::-1]),
    ]
    assert all(p.header == p.id for p in long)
    pep = pd.read_parquet(work / "peptides.norm.parquet")
    assert list(pep["peptide"]) == ["MKVLAAGIVGR"]
    assert list(pep["protein_ids"].iloc[0]) == ["S1"]
    routing = load_report(out)["routing"]
    assert routing["n_proteins_long"] == 2
    assert routing["n_proteins_peptide_path"] == 1
    assert routing["n_proteins_dropped"] == 1


def test_abundance_is_normalised(tmp_path: Path) -> None:
    ab = write(tmp_path / "abundance.tsv", "protein_id\tvalue\nP1\t1000\n")
    out = tmp_path / "out"
    r = invoke(
        "--genome", genome_120kb(tmp_path), "--proteins", proteins(tmp_path), "--abundance", ab,
        *LOW, "--outdir", out,
    )  # fmt: skip
    assert_code(r, ExitCode.SUCCESS)
    assert (out / "work" / "preflight" / "abundance.norm.parquet").is_file()
    report = load_report(out)
    assert report["routing"]["abundance_used"] is True
    assert {"abundance"} <= {c["input"] for c in report["checks"]}


# ---------------------------------------------------------------- simulated genome end to end


def test_simulated_genome_with_bakta_annotation_and_evidence(tmp_path: Path) -> None:
    sim = simulate(SimConfig(genome_length=120_000), seed=7, outdir=tmp_path / "sim")
    feats = json.loads(sim.bakta_json.read_text())["features"][:5]
    prot = write(tmp_path / "proteins.faa", "".join(f">{f['locus']}\n{f['aa']}\n" for f in feats))
    rows = []
    for f in feats:
        for piece in f["aa"].replace("K", "K\n").replace("R", "R\n").split("\n"):
            if 7 <= len(piece) <= 40:
                rows.append(f"{piece}\t{f['locus']}\n")
    assert rows
    pep = write(tmp_path / "peptides.tsv", "peptide\tprotein_ids\n" + "".join(rows))
    out = tmp_path / "out"
    r = invoke(
        "--genome", sim.genome_fasta, "--proteins", prot, "--peptides", pep,
        "--annotation", sim.bakta_json, *LOW, "--outdir", out,
    )  # fmt: skip
    assert_code(r, ExitCode.SUCCESS)
    report = load_report(out)
    assert report["summary"]["n_error"] == 0
    assert report["replicons"] == [
        {"id": "chr1", "length": 120_000, "topology": "circular", "kind": "chromosome"}
    ]
    annotation = {i["role"]: i for i in report["inputs"]}["annotation"]
    assert annotation["file"] == sim.bakta_json.name
    assert annotation["sha256"] == hashlib.sha256(sim.bakta_json.read_bytes()).hexdigest()
    assert annotation["size"] == sim.bakta_json.stat().st_size
    work = out / "work" / "preflight"
    for name in (*WORK_FILES, "proteins.long.faa"):
        assert (work / name).is_file(), name
    assert (work / "replicons.tsv").read_text() == "chr1\tchr1\tchromosome\tc\t-\n"


# ---------------------------------------------------------------- --annotation kinds


@pytest.mark.parametrize("name", ["a.json", "a.gff3", "a.gff", "a.gbff", "a.gbk", "a.gff3.gz"])
def test_annotation_known_kinds_are_accepted(tmp_path: Path, name: str) -> None:
    path = tmp_path / name
    if name.endswith(".gz"):
        path.write_bytes(gzip.compress(b"##gff-version 3\n", mtime=0))
    else:
        path.write_text("{}\n" if name.endswith(".json") else "##gff-version 3\n")
    r = invoke(
        "--genome", genome_120kb(tmp_path), "--proteins", proteins(tmp_path),
        "--annotation", path, *LOW,
    )  # fmt: skip
    assert_code(r, ExitCode.SUCCESS)


def test_annotation_of_unknown_kind_is_input_error(tmp_path: Path) -> None:
    path = write(tmp_path / "annotation.txt", "hello\n")
    out = tmp_path / "out"
    r = invoke(
        "--genome", genome_120kb(tmp_path), "--proteins", proteins(tmp_path),
        "--annotation", path, *LOW, "--outdir", out,
    )  # fmt: skip
    assert_code(r, ExitCode.INVALID_INPUT)
    report = load_report(out)
    assert report["summary"]["n_error"] >= 1
    assert "annotation" in {i["role"] for i in report["inputs"]}


# ---------------------------------------------------------------- environment and internal errors


def test_unreadable_db_manifest_is_environment_error_exit_2(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    db = tmp_path / "db_bad"
    db.mkdir()
    write(db / "db_manifest.json", "{not json")
    monkeypatch.setenv("GIAP_DB_DIR", str(db))
    out = tmp_path / "out"
    r = invoke(
        "--genome", genome_120kb(tmp_path), "--proteins", proteins(tmp_path), *LOW, "--outdir", out
    )
    assert_code(r, ExitCode.ENVIRONMENT)
    env = [c for c in load_report(out)["checks"] if c["rule_id"] == "environment.tools"]
    assert len(env) == 1
    assert env[0]["level"] == "error"
    assert env[0]["exit_code"] == 2


def test_input_error_wins_over_environment_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    db = tmp_path / "db_bad"
    db.mkdir()
    write(db / "db_manifest.json", "{not json")
    monkeypatch.setenv("GIAP_DB_DIR", str(db))
    g = write(tmp_path / "g.fa", ">c\nACGU\n")
    r = invoke("--genome", g, "--proteins", proteins(tmp_path), "--strict")
    assert_code(r, ExitCode.INVALID_INPUT)


def test_missing_db_manifest_is_info_under_preflight(tmp_path: Path) -> None:
    out = tmp_path / "out"
    r = invoke(
        "--genome", genome_120kb(tmp_path), "--proteins", proteins(tmp_path), *LOW,
        "--outdir", out, "--strict",
    )  # fmt: skip
    assert_code(r, ExitCode.SUCCESS)
    env = [c for c in load_report(out)["checks"] if c["input"] == "environment"]
    assert env, "environment checks must be recorded"
    assert all(c["level"] in ("ok", "info") for c in env)
    assert any("db_manifest" in (c["message"] + json.dumps(c["value"])) for c in env)


def test_exception_inside_a_rule_exits_5_with_report_message(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    import giap.preflight.genome as genome_mod
    from giap.cli import REPORT_MESSAGE, main

    def boom(lengths: object) -> int:
        raise RuntimeError("rule bug")

    monkeypatch.setattr(genome_mod, "n50", boom)
    g, prot = genome_120kb(tmp_path), proteins(tmp_path)
    with pytest.raises(SystemExit) as info:
        main(["preflight", "--genome", str(g), "--proteins", str(prot)])
    assert info.value.code == ExitCode.INTERNAL
    assert REPORT_MESSAGE in capsys.readouterr().err


def test_preflight_options_match_spec() -> None:
    from giap.cli import cli

    names = {p.name for p in cli.commands["preflight"].params}
    assert {
        "genome",
        "proteins",
        "peptides",
        "abundance",
        "annotation",
        "outdir",
        "strict",
    } <= names
