"""Pool per-run GIAP findings over datasets and search arms (SPEC §11).

Synthetic work directories are written under ``tmp_path``; the GIAP tables hold a few real columns
of docs/reference/tables.md and use GIAP's own conventions (``NA`` for missing, ``true``/``false``).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from giap_campaign.aggregate import aggregate

PREFIX = ["dataset", "tax_id", "organism", "genome_accession", "arm"]
RUNS_COLUMNS = [
    *PREFIX,
    "state",
    "acquisition",
    "giap_version",
    "n_genes_detected",
    "n_novel_called",
    "n_alt_high",
    "notes",
]
HYPOTHETICAL_COLUMNS = [
    *PREFIX,
    "locus_tag",
    "product",
    "contig",
    "start",
    "end",
    "strand",
    "length_aa",
    "detection_tier",
    "gene_q_value",
    "n_peptides_unique",
    "abundance_log2",
]
OUTPUTS = [
    "campaign_novel_orfs.tsv",
    "campaign_novel_orfs.faa",
    "campaign_variants.tsv",
    "campaign_operon_flags.tsv",
    "campaign_islands.tsv",
    "campaign_detected_hypothetical.tsv",
    "campaign_runs.tsv",
]
MISSING_SPELLINGS = {"NA", "nan", "NaN", "None", "<NA>", "null"}

ECOLI = "Escherichia coli str. K-12 substr. MG1655"
SAUR = "Staphylococcus aureus subsp. aureus NCTC 8325"

# --- GIAP table fragments -----------------------------------------------------------------------

NOVEL_COLUMNS = [
    "orf_id",
    "class",
    "contig",
    "strand",
    "start_longest",
    "start_best_ostir",
    "end",
    "n_peptides",
    "best_peptide_q",
    "class_fdr",
    "called",
    "flags",
]
VARIANT_COLUMNS = [
    "variant_id",
    "contig",
    "codon_start",
    "strand",
    "gene",
    "ref_aa",
    "alt_aa",
    "alt_type",
    "call",
    "confidence",
    "n_ref",
    "n_alt",
    "class_error_rate",
    "source",
]
GENE_COLUMNS = [
    "locus_tag",
    "gene",
    "product",
    "feature_type",
    "hypothetical",
    "contig",
    "start",
    "end",
    "strand",
    "length_aa",
    "detected",
    "detection_tier",
    "detection_status",
    "gene_q_value",
    "n_peptides_unique",
    "abundance_log2",
    "cog_id",
]
OPERON_COLUMNS = ["pair_id", "gene_a", "gene_b", "uniop_class", "flag", "q_value", "evidence_note"]
ISLAND_COLUMNS = ["island_id", "direction", "contig", "start", "end", "n_tus", "p_value", "q_value"]


def _novel(orf: str, called: bool, q: str = "0.0012") -> list[str]:
    pos = orf.split("_")[-2]
    strand = "+" if orf.endswith("_p") else "-"
    return [
        orf,
        "intergenic",
        "c1",
        strand,
        pos,
        "NA",
        str(int(pos) + 299),
        "2",
        q,
        "unestimated",
        "true" if called else "false",
        "NA",
    ]


def _variant(vid: int, pos: int, strand: str, alt: str, call: str, conf: str) -> list[str]:
    return [
        f"VAR_{vid}",
        "c1",
        str(pos),
        strand,
        "b0001",
        "E",
        alt,
        "frameshift" if alt.startswith("fs") else "substitution",
        call,
        conf,
        "1",
        "3",
        "1.5e-05",
        "database",
    ]


def _gene(
    tag: str,
    hypothetical: str,
    detected: str,
    *,
    q: str = "0.004",
    abundance: str = "21.5",
) -> list[str]:
    return [
        tag,
        "NA",
        "hypothetical protein" if hypothetical == "true" else "DNA polymerase III subunit",
        "cds",
        hypothetical,
        "c1",
        "101",
        "400",
        "+",
        "99",
        detected,
        "decoy_picked",
        "detected" if detected == "true" else "not_detected",
        q,
        "3",
        abundance,
        "COG0592",
    ]


def _operon(n: int) -> list[str]:
    return [f"PAIR_{n}", f"g{n}", f"g{n + 1}", "OP", "unexplained_discordance", "0.01", "note"]


def _island(n: int) -> list[str]:
    return [f"ISL_{n}", "high", "c1", "1000", "9000", "4", "0.003", "NA"]


# --- work directory builders --------------------------------------------------------------------


def _write_tsv(path: Path, columns: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(columns), *("\t".join(r) for r in rows)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_arm(
    workdir: Path,
    arm: str,
    *,
    version: str | None = "1.0.0",
    tables: dict[str, Any] | None = None,
) -> Path:
    """``tables`` maps a GIAP file name to rows, raw ``bytes`` or FASTA ``str``."""
    out = workdir / f"giap_{arm}"
    (out / "tables").mkdir(parents=True, exist_ok=True)
    if version is not None:
        manifest = {"giap_version": version, "status": "success", "exit_code": 0}
        (out / "run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    columns = {
        "novel_orfs.tsv": NOVEL_COLUMNS,
        "variants.tsv": VARIANT_COLUMNS,
        "genes.tsv": GENE_COLUMNS,
        "operon_flags.tsv": OPERON_COLUMNS,
        "islands.tsv": ISLAND_COLUMNS,
    }
    for name, content in (tables or {}).items():
        path = out / "tables" / name
        if isinstance(content, bytes):
            path.write_bytes(content)
        elif isinstance(content, str):
            path.write_text(content, encoding="utf-8")
        else:
            _write_tsv(path, columns[name], content)
    return out


def _write_status(workdir: Path, steps: dict[str, tuple[str, str]], acquisition: Any) -> None:
    workdir.mkdir(parents=True, exist_ok=True)
    doc: dict[str, Any] = {
        "steps": {name: {"state": state, "message": msg} for name, (state, msg) in steps.items()}
    }
    if acquisition is not None:
        doc["acquisition"] = acquisition
    (workdir / "status.json").write_text(json.dumps(doc), encoding="utf-8")


def _job(
    n: int, root: Path, dataset: str, tax_id: int, organism: str, genome: str
) -> dict[str, Any]:
    return {
        "job": n,
        "dataset": dataset,
        "tax_id": tax_id,
        "organism": organism,
        "genome_accession": genome,
        "raw_format": "thermo_raw",
        "acquisition_metadata": "unknown",
        "run_files": "a.raw",
        "run_urls": "https://ftp.pride.ebi.ac.uk/a.raw",
        "run_bytes": "100",
        "run_checksums": "",
        "run_gb": 0.0,
        "denovo_runs": "a.raw",
        "workdir": str(root / f"{dataset}_{tax_id}"),
    }


GARBAGE = bytes([0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A, 0x00]) + bytes(range(128, 256)) * 4

FAA_1_DB = (
    ">NOVEL_c1_900_p|longest\nMKVLAAGIVGLLLAAPAQA\n"
    ">NOVEL_c1_100_p|longest,best_ostir candidate protein\n"
    "MSTNPKPQRKTKRNTNRRPQDVKFPGG\nGGQIVGGVYLLPRRGPRLGVRATR\n"
)
FAA_1_DENOVO = ">NOVEL_c1_500_m|best_ostir\nMAAAKKK\n"
FAA_3_DB = ">NOVEL_c1_100_p|longest\nMQQQ\n"


def _build_campaign(root: Path) -> pd.DataFrame:
    """Five jobs (given out of order):

    - PXD000001 / 511145 / GCF_A: both arms done (status.json, acquisition DDA).
    - PXD000002 / 511145 / GCF_A: db done; giap_denovo failed (its directory has tables but no
      manifest, so it is ignored); acquisition DIA.
    - PXD000003 / 1280 / GCF_B: both arms done, no status.json; denovo has no novel_orfs.tsv or
      genes.tsv and an unreadable variants.tsv.
    - PXD000004 / 10 / GCF_C: fetch_raw failed, no giap directories.
    - PXD000004 / 9 / GCF_D: nothing at all (not_run).
    """
    root.mkdir(parents=True, exist_ok=True)
    jobs = [
        _job(0, root, "PXD000003", 1280, SAUR, "GCF_B"),
        _job(1, root, "PXD000004", 10, "Cellvibrio sp.", "GCF_C"),
        _job(2, root, "PXD000002", 511145, ECOLI, "GCF_A"),
        _job(3, root, "PXD000004", 9, "Buchnera aphidicola", "GCF_D"),
        _job(4, root, "PXD000001", 511145, ECOLI, "GCF_A"),
    ]
    w1 = root / "PXD000001_511145"
    _write_arm(
        w1,
        "db",
        tables={
            "novel_orfs.tsv": [
                _novel("NOVEL_c1_900_p", False),
                _novel("NOVEL_c1_100_p", True, q="1.5e-05"),
                _novel("NOVEL_c1_500_m", False),
            ],
            "novel_orfs.faa": FAA_1_DB,
            "variants.tsv": [
                _variant(1, 300, "+", "K", "ALT", "high"),
                _variant(2, 600, "-", "fs+1", "ALT", "candidate"),
            ],
            "genes.tsv": [
                _gene("b0001", "true", "true"),
                _gene("b0002", "true", "false"),
                _gene("b0003", "false", "true"),
                _gene("b0004", "true", "true", q="NA", abundance="NA"),
            ],
            "operon_flags.tsv": [_operon(1), _operon(2)],
            "islands.tsv": [_island(1)],
        },
    )
    _write_arm(
        w1,
        "denovo",
        version="1.1.0.dev3",
        tables={
            "novel_orfs.tsv": [_novel("NOVEL_c1_100_p", True), _novel("NOVEL_c1_500_m", True)],
            "novel_orfs.faa": FAA_1_DENOVO,
            "variants.tsv": [
                _variant(1, 300, "+", "K", "ALT", "high"),
                _variant(2, 300, "+", "R", "mixed", "high"),
            ],
            "genes.tsv": [_gene("b0001", "true", "true"), _gene("b0002", "true", "true")],
            "operon_flags.tsv": [_operon(7)],
            "islands.tsv": [_island(7)],
        },
    )
    done = {
        s: ("done", "")
        for s in (
            "fetch_genome",
            "annotate",
            "makedb",
            "fetch_raw",
            "acquisition",
            "to_mzml",
            "search_db",
            "denovo",
            "convert",
            "giap_db",
            "giap_denovo",
            "compare",
            "cleanup",
        )
    }
    _write_status(w1, done, "DDA")

    w2 = root / "PXD000002_511145"
    _write_arm(
        w2,
        "db",
        tables={
            "novel_orfs.tsv": [_novel("NOVEL_c1_100_p", True), _novel("NOVEL_c1_900_p", False)],
            "variants.tsv": [_variant(1, 300, "+", "K", "ALT", "candidate")],
            "genes.tsv": [_gene("b0001", "true", "true")],
            "islands.tsv": [_island(1)],
        },
    )
    _write_arm(
        w2,
        "denovo",
        version=None,
        tables={
            "novel_orfs.tsv": [_novel("NOVEL_c1_900_p", True)],
            "variants.tsv": [_variant(1, 300, "+", "K", "ALT", "high")],
        },
    )
    steps2 = {k: v for k, v in done.items() if k not in {"giap_denovo", "compare", "cleanup"}}
    steps2["giap_denovo"] = ("failed", "giap run exited with 1")
    _write_status(w2, steps2, "DIA")

    w3 = root / "PXD000003_1280"
    _write_arm(
        w3,
        "db",
        tables={
            "novel_orfs.tsv": [_novel("NOVEL_c1_100_p", True)],
            "novel_orfs.faa": FAA_3_DB,
            "variants.tsv": [_variant(1, 300, "+", "K", "ALT", "high")],
            "genes.tsv": [_gene("SAOUHSC_00001", "true", "true")],
            "operon_flags.tsv": [_operon(3)],
        },
    )
    _write_arm(w3, "denovo", tables={"variants.tsv": GARBAGE})

    w4 = root / "PXD000004_10"
    _write_status(
        w4,
        {"fetch_genome": ("done", ""), "fetch_raw": ("failed", "checksum mismatch: a.raw")},
        None,
    )
    return pd.DataFrame(jobs)


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _keys(df: pd.DataFrame) -> list[tuple[str, str, str]]:
    return list(zip(df["dataset"], df["tax_id"], df["arm"], strict=True))


def _fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            records.append((line[1:], ""))
        elif line:
            header, seq = records[-1]
            records[-1] = (header, seq + line)
    return records


@pytest.fixture(scope="module")
def campaign(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, dict[str, int]]:
    base = tmp_path_factory.mktemp("campaign")
    jobs = _build_campaign(base / "work")
    out = base / "out"
    counts = aggregate(jobs, out)
    return out, counts


# --- return value and files ---------------------------------------------------------------------


def test_returns_rows_written_per_table(campaign: tuple[Path, dict[str, int]]) -> None:
    out, counts = campaign
    assert counts == {
        "campaign_novel_orfs.tsv": 8,
        "campaign_novel_orfs.faa": 4,
        "campaign_variants.tsv": 6,
        "campaign_operon_flags.tsv": 3,
        "campaign_islands.tsv": 2,
        "campaign_detected_hypothetical.tsv": 4,
        "campaign_runs.tsv": 10,
    }
    for name, n in counts.items():
        if name.endswith(".tsv"):
            assert len(_read(out / name)) == n, name
        else:
            assert len(_fasta(out / name)) == n, name


def test_no_missing_value_spellings_in_any_output(campaign: tuple[Path, dict[str, int]]) -> None:
    out, _ = campaign
    for name in OUTPUTS:
        if not name.endswith(".tsv"):
            continue
        text = (out / name).read_text(encoding="utf-8")
        for line in text.splitlines():
            for cell in line.split("\t"):
                assert cell not in MISSING_SPELLINGS, (name, line)


def test_tsv_conventions(campaign: tuple[Path, dict[str, int]]) -> None:
    out, _ = campaign
    for name in OUTPUTS:
        if not name.endswith(".tsv"):
            continue
        raw = (out / name).read_bytes()
        raw.decode("utf-8")
        assert b"\r" not in raw, name
        assert raw.endswith(b"\n"), name
        header = raw.split(b"\n", 1)[0].decode().split("\t")
        assert header[: len(PREFIX)] == PREFIX, name


# --- novel ORFs ---------------------------------------------------------------------------------


def test_novel_orfs_rows_prefix_and_order(campaign: tuple[Path, dict[str, int]]) -> None:
    out, _ = campaign
    df = _read(out / "campaign_novel_orfs.tsv")
    assert list(df.columns) == [*PREFIX, *NOVEL_COLUMNS, "arms_called", "n_datasets_called"]
    assert _keys(df) == [
        ("PXD000001", "511145", "db"),
        ("PXD000001", "511145", "db"),
        ("PXD000001", "511145", "db"),
        ("PXD000001", "511145", "denovo"),
        ("PXD000001", "511145", "denovo"),
        ("PXD000002", "511145", "db"),
        ("PXD000002", "511145", "db"),
        ("PXD000003", "1280", "db"),
    ]
    # the table's own order within a run is kept
    assert list(df["orf_id"][:3]) == ["NOVEL_c1_900_p", "NOVEL_c1_100_p", "NOVEL_c1_500_m"]
    assert set(df.loc[df["dataset"] == "PXD000003", "organism"]) == {SAUR}
    assert set(df.loc[df["dataset"] == "PXD000003", "genome_accession"]) == {"GCF_B"}
    assert set(df.loc[df["dataset"] == "PXD000001", "organism"]) == {ECOLI}


def test_novel_orfs_arms_called_and_n_datasets_called(
    campaign: tuple[Path, dict[str, int]],
) -> None:
    out, _ = campaign
    df = _read(out / "campaign_novel_orfs.tsv")
    got = list(
        zip(
            df["dataset"],
            df["arm"],
            df["orf_id"],
            df["arms_called"],
            df["n_datasets_called"],
            strict=True,
        )
    )
    assert got == [
        # called in both arms of PXD000001; PXD000002 calls it too (same genome GCF_A);
        # PXD000003 calls the same orf_id on another genome and does not count
        ("PXD000001", "db", "NOVEL_c1_900_p", "", "0"),
        ("PXD000001", "db", "NOVEL_c1_100_p", "db;denovo", "2"),
        # not called in the db arm, called in denovo: arms_called holds denovo on both rows
        ("PXD000001", "db", "NOVEL_c1_500_m", "denovo", "1"),
        ("PXD000001", "denovo", "NOVEL_c1_100_p", "db;denovo", "2"),
        ("PXD000001", "denovo", "NOVEL_c1_500_m", "denovo", "1"),
        ("PXD000002", "db", "NOVEL_c1_100_p", "db", "2"),
        # PXD000002's denovo arm called it but has no manifest: ignored
        ("PXD000002", "db", "NOVEL_c1_900_p", "", "0"),
        ("PXD000003", "db", "NOVEL_c1_100_p", "db", "1"),
    ]


def test_novel_orfs_cells_kept_unchanged(campaign: tuple[Path, dict[str, int]]) -> None:
    """Per-run values are copied as text (R4): booleans stay ``true``/``false``, q-values and
    integers keep their spelling, ``unestimated`` stays, GIAP's ``NA`` becomes an empty cell."""
    out, _ = campaign
    df = _read(out / "campaign_novel_orfs.tsv")
    row = df.iloc[1]
    assert row["orf_id"] == "NOVEL_c1_100_p"
    assert row["called"] == "true"
    assert row["best_peptide_q"] == "1.5e-05"
    assert row["class_fdr"] == "unestimated"
    assert row["start_longest"] == "100"
    assert row["end"] == "399"
    assert row["start_best_ostir"] == ""
    assert row["flags"] == ""
    assert row["tax_id"] == "511145"
    assert set(df["called"]) == {"true", "false"}


def test_novel_orfs_fasta_headers_prefixed(campaign: tuple[Path, dict[str, int]]) -> None:
    out, _ = campaign
    records = _fasta(out / "campaign_novel_orfs.faa")
    assert records == [
        ("PXD000001|db|GCF_A|NOVEL_c1_900_p|longest", "MKVLAAGIVGLLLAAPAQA"),
        (
            "PXD000001|db|GCF_A|NOVEL_c1_100_p|longest,best_ostir candidate protein",
            "MSTNPKPQRKTKRNTNRRPQDVKFPGGGGQIVGGVYLLPRRGPRLGVRATR",
        ),
        ("PXD000001|denovo|GCF_A|NOVEL_c1_500_m|best_ostir", "MAAAKKK"),
        ("PXD000003|db|GCF_B|NOVEL_c1_100_p|longest", "MQQQ"),
    ]


# --- variants -----------------------------------------------------------------------------------


def test_variants_rows_and_n_datasets(campaign: tuple[Path, dict[str, int]]) -> None:
    out, _ = campaign
    df = _read(out / "campaign_variants.tsv")
    assert list(df.columns) == [*PREFIX, *VARIANT_COLUMNS, "n_datasets"]
    got = list(
        zip(
            df["dataset"],
            df["arm"],
            df["codon_start"],
            df["strand"],
            df["alt_aa"],
            df["n_datasets"],
            strict=True,
        )
    )
    assert got == [
        # same site and ALT in PXD000001 (both arms: one dataset) and PXD000002, genome GCF_A
        ("PXD000001", "db", "300", "+", "K", "2"),
        ("PXD000001", "db", "600", "-", "fs+1", "1"),
        ("PXD000001", "denovo", "300", "+", "K", "2"),
        # another ALT at the same codon is another key
        ("PXD000001", "denovo", "300", "+", "R", "1"),
        ("PXD000002", "db", "300", "+", "K", "2"),
        # same contig/codon/strand/ALT on another genome: not pooled
        ("PXD000003", "db", "300", "+", "K", "1"),
    ]
    assert set(df["class_error_rate"]) == {"1.5e-05"}


# --- db-arm-only tables -------------------------------------------------------------------------


def test_operon_flags_db_arm_only(campaign: tuple[Path, dict[str, int]]) -> None:
    out, _ = campaign
    df = _read(out / "campaign_operon_flags.tsv")
    assert list(df.columns) == [*PREFIX, *OPERON_COLUMNS]
    assert set(df["arm"]) == {"db"}
    assert list(zip(df["dataset"], df["pair_id"], strict=True)) == [
        ("PXD000001", "PAIR_1"),
        ("PXD000001", "PAIR_2"),
        ("PXD000003", "PAIR_3"),
    ]


def test_islands_db_arm_only(campaign: tuple[Path, dict[str, int]]) -> None:
    out, _ = campaign
    df = _read(out / "campaign_islands.tsv")
    assert list(df.columns) == [*PREFIX, *ISLAND_COLUMNS]
    assert list(zip(df["dataset"], df["arm"], df["island_id"], strict=True)) == [
        ("PXD000001", "db", "ISL_1"),
        ("PXD000002", "db", "ISL_1"),
    ]
    assert set(df["q_value"]) == {""}
    assert set(df["p_value"]) == {"0.003"}


def test_detected_hypothetical_filter_and_columns(campaign: tuple[Path, dict[str, int]]) -> None:
    out, _ = campaign
    df = _read(out / "campaign_detected_hypothetical.tsv")
    assert list(df.columns) == HYPOTHETICAL_COLUMNS
    assert set(df["arm"]) == {"db"}
    assert list(zip(df["dataset"], df["locus_tag"], strict=True)) == [
        ("PXD000001", "b0001"),
        ("PXD000001", "b0004"),
        ("PXD000002", "b0001"),
        ("PXD000003", "SAOUHSC_00001"),
    ]
    first = df.iloc[0].to_dict()
    assert first == {
        "dataset": "PXD000001",
        "tax_id": "511145",
        "organism": ECOLI,
        "genome_accession": "GCF_A",
        "arm": "db",
        "locus_tag": "b0001",
        "product": "hypothetical protein",
        "contig": "c1",
        "start": "101",
        "end": "400",
        "strand": "+",
        "length_aa": "99",
        "detection_tier": "decoy_picked",
        "gene_q_value": "0.004",
        "n_peptides_unique": "3",
        "abundance_log2": "21.5",
    }
    assert df.iloc[1]["gene_q_value"] == ""
    assert df.iloc[1]["abundance_log2"] == ""


# --- campaign_runs ------------------------------------------------------------------------------


def test_runs_one_row_per_job_and_arm_sorted(campaign: tuple[Path, dict[str, int]]) -> None:
    out, _ = campaign
    df = _read(out / "campaign_runs.tsv")
    assert list(df.columns) == RUNS_COLUMNS
    assert _keys(df) == [
        ("PXD000001", "511145", "db"),
        ("PXD000001", "511145", "denovo"),
        ("PXD000002", "511145", "db"),
        ("PXD000002", "511145", "denovo"),
        ("PXD000003", "1280", "db"),
        ("PXD000003", "1280", "denovo"),
        # numeric tax_id order: 9 before 10
        ("PXD000004", "9", "db"),
        ("PXD000004", "9", "denovo"),
        ("PXD000004", "10", "db"),
        ("PXD000004", "10", "denovo"),
    ]


def test_runs_states(campaign: tuple[Path, dict[str, int]]) -> None:
    out, _ = campaign
    df = _read(out / "campaign_runs.tsv")
    assert list(df["state"]) == [
        "done",
        "done",
        "done",
        "failed:giap_denovo",
        "done",
        "done",
        "not_run",
        "not_run",
        "failed:fetch_raw",
        "failed:fetch_raw",
    ]


def test_runs_versions_and_acquisition(campaign: tuple[Path, dict[str, int]]) -> None:
    out, _ = campaign
    df = _read(out / "campaign_runs.tsv").set_index(["dataset", "tax_id", "arm"])
    assert df.loc[("PXD000001", "511145", "db"), "giap_version"] == "1.0.0"
    assert df.loc[("PXD000001", "511145", "denovo"), "giap_version"] == "1.1.0.dev3"
    assert df.loc[("PXD000002", "511145", "denovo"), "giap_version"] == ""
    assert df.loc[("PXD000004", "9", "db"), "giap_version"] == ""
    assert df.loc[("PXD000001", "511145", "db"), "acquisition"] == "DDA"
    assert df.loc[("PXD000001", "511145", "denovo"), "acquisition"] == "DDA"
    assert df.loc[("PXD000002", "511145", "db"), "acquisition"] == "DIA"
    assert df.loc[("PXD000004", "10", "db"), "acquisition"] == ""


def test_runs_counts(campaign: tuple[Path, dict[str, int]]) -> None:
    """Counting rules: ``n_genes_detected`` = genes.tsv rows with ``detected`` true (any product);
    ``n_novel_called`` = novel_orfs.tsv rows with ``called`` true; ``n_alt_high`` = variants.tsv
    rows with ``confidence`` high and ``call`` ALT (mixed calls excluded). A count whose source
    table is missing or unreadable, or a run that has no manifest, is an empty cell."""
    out, _ = campaign
    df = _read(out / "campaign_runs.tsv")
    got = list(zip(df["n_genes_detected"], df["n_novel_called"], df["n_alt_high"], strict=True))
    assert got == [
        ("3", "1", "1"),  # PXD000001 db: b0001, b0003, b0004 detected
        ("2", "2", "1"),  # PXD000001 denovo: the mixed R call is not counted
        ("1", "1", "0"),  # PXD000002 db: the K call is a candidate
        ("", "", ""),  # PXD000002 denovo failed (its tables are not read)
        ("1", "1", "1"),
        ("", "", ""),  # PXD000003 denovo: genes and novel ORFs missing, variants unreadable
        ("", "", ""),
        ("", "", ""),
        ("", "", ""),
        ("", "", ""),
    ]


def test_runs_notes_unreadable_table(campaign: tuple[Path, dict[str, int]]) -> None:
    out, _ = campaign
    df = _read(out / "campaign_runs.tsv").set_index(["dataset", "tax_id", "arm"])
    assert "variants" in df.loc[("PXD000003", "1280", "denovo"), "notes"]
    # missing tables are skipped silently; the readable runs have no unreadable-table note
    assert "variants" not in df.loc[("PXD000003", "1280", "db"), "notes"]
    assert "novel_orfs" not in df.loc[("PXD000003", "1280", "denovo"), "notes"]
    assert df.loc[("PXD000003", "1280", "denovo"), "state"] == "done"


# --- edge cases ---------------------------------------------------------------------------------


def test_job_order_does_not_change_outputs(tmp_path: Path) -> None:
    jobs = _build_campaign(tmp_path / "work")
    a = aggregate(jobs, tmp_path / "a")
    b = aggregate(jobs.iloc[::-1].reset_index(drop=True), tmp_path / "b")
    assert a == b
    for name in OUTPUTS:
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes(), name


def test_only_unrun_jobs(tmp_path: Path) -> None:
    root = tmp_path / "work"
    jobs = pd.DataFrame(
        [
            _job(0, root, "PXD000002", 562, "Escherichia coli", "GCF_X"),
            _job(1, root, "PXD000001", 562, "Escherichia coli", "GCF_X"),
        ]
    )
    counts = aggregate(jobs, tmp_path / "out")
    assert counts["campaign_runs.tsv"] == 4
    for name in OUTPUTS:
        if name != "campaign_runs.tsv":
            assert counts.get(name, 0) == 0, name
    df = _read(tmp_path / "out" / "campaign_runs.tsv")
    assert list(df["state"]) == ["not_run"] * 4
    assert _keys(df) == [
        ("PXD000001", "562", "db"),
        ("PXD000001", "562", "denovo"),
        ("PXD000002", "562", "db"),
        ("PXD000002", "562", "denovo"),
    ]


def test_manifest_without_tables_is_done_with_no_rows(tmp_path: Path) -> None:
    root = tmp_path / "work"
    jobs = pd.DataFrame([_job(0, root, "PXD000001", 562, "Escherichia coli", "GCF_X")])
    _write_arm(root / "PXD000001_562", "db", version="1.0.0")
    counts = aggregate(jobs, tmp_path / "out")
    assert counts["campaign_runs.tsv"] == 2
    assert counts.get("campaign_novel_orfs.tsv", 0) == 0
    df = _read(tmp_path / "out" / "campaign_runs.tsv")
    assert list(df["state"]) == ["done", "not_run"]
    assert df.iloc[0]["giap_version"] == "1.0.0"
    assert df.iloc[0]["notes"] == ""


def test_unreadable_novel_table_leaves_other_tables(tmp_path: Path) -> None:
    root = tmp_path / "work"
    jobs = pd.DataFrame([_job(0, root, "PXD000001", 562, "Escherichia coli", "GCF_X")])
    _write_arm(
        root / "PXD000001_562",
        "db",
        tables={
            "novel_orfs.tsv": GARBAGE,
            "variants.tsv": [_variant(1, 3, "+", "K", "ALT", "high")],
        },
    )
    counts = aggregate(jobs, tmp_path / "out")
    assert counts.get("campaign_novel_orfs.tsv", 0) == 0
    assert counts["campaign_variants.tsv"] == 1
    df = _read(tmp_path / "out" / "campaign_runs.tsv")
    db = df.iloc[0]
    assert db["arm"] == "db"
    assert db["state"] == "done"
    assert "novel_orfs" in db["notes"]
    assert db["n_alt_high"] == "1"
    assert db["n_novel_called"] == ""
