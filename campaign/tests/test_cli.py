"""``python -m giap_campaign`` (campaign/SPEC.md §12).

Exit codes: 0 success; 1 invalid input or a failed job; 2 network or environment error.

Readings: ``cli.main(argv) -> int`` never calls ``sys.exit`` for these cases; ``run-job`` reads row
``--index`` of ``jobs.tsv``, builds ``env`` from ``campaign/pins.tsv`` (``$GIAP_DB_DIR`` and
``$GIAP_ROOT`` expanded from the process environment; ``env["pins"]`` = name -> sha256) and
``env["diann_cfg"]`` = ``campaign/diann.cfg``, and returns ``driver.run_job``'s code. The network
layer is ``giap_campaign.http.get_json``.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from giap_campaign import aggregate, cli, driver, http

REPO = Path(__file__).resolve().parents[2]
CAMPAIGN = REPO / "campaign"
PRIDE = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01"
ACC = "GCF_000005845.2"
PLAN_COLUMNS = [
    "job",
    "dataset",
    "tax_id",
    "organism",
    "genome_accession",
    "raw_format",
    "acquisition_metadata",
    "run_files",
    "run_urls",
    "run_bytes",
    "run_checksums",
    "run_gb",
    "denovo_runs",
    "workdir",
]


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


# ---------------------------------------------------------------- --help (subprocess)


def test_module_help_lists_the_subcommands() -> None:
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(CAMPAIGN / "src"), str(REPO / "src"), env.get("PYTHONPATH", "")]
    )
    proc = subprocess.run(
        [sys.executable, "-m", "giap_campaign", "--help"],
        capture_output=True,
        text=True,
        env=env,
        cwd=REPO,
        timeout=120,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    for sub in ("catalog", "plan", "run-job", "convert", "aggregate"):
        assert sub in proc.stdout


# ---------------------------------------------------------------- convert

SAGE_HEADER = (
    "psm_id\tpeptide\tproteins\tnum_proteins\tfilename\tscannr\trank\tlabel\tcharge\t"
    "hyperscore\tsage_discriminant_score\tposterior_error\tspectrum_q\tpeptide_q\n"
)
SAGE_ROWS = [
    "1\tPEPTIDEK\tsp|P0A7V0|RS2_ECOLI\t1\trun1.mzML\t"
    "controllerType=0 controllerNumber=1 scan=10\t1\t1\t2\t30.1\t1.2\t-3.0\t0.001\t0.001\n",
    "2\tM[+15.9949]PEPTIDER\trev_sp|P0A7V3|RS3_ECOLI\t1\t/data/run1.mzML\t"
    "controllerType=0 controllerNumber=1 scan=11\t1\t-1\t3\t10.0\t-0.5\t-0.1\t0.5\t0.5\n",
    "3\tAAAAK\tsp|P0A7V0|RS2_ECOLI\t1\trun1.mzML\t"
    "controllerType=0 controllerNumber=1 scan=10\t2\t1\t2\t9.0\t-1.0\t-0.05\t0.9\t0.9\n",
]


def test_convert_sage_writes_a_peptide_table(tmp_path: Path) -> None:
    src = tmp_path / "results.sage.tsv"
    src.write_text(SAGE_HEADER + "".join(SAGE_ROWS))
    out = tmp_path / "peptides.raw.tsv"
    assert cli.main(["convert", "--format", "sage", str(src), "--out", str(out)]) == 0
    table = _read_tsv(out)
    assert "peptide" in table.columns
    assert set(table["peptide"]) == {"PEPTIDEK", "MPEPTIDER"}  # rank 2 dropped, mods stripped


def test_convert_missing_input_is_exit_1(tmp_path: Path) -> None:
    out = tmp_path / "peptides.raw.tsv"
    rc = cli.main(["convert", "--format", "sage", str(tmp_path / "nope.tsv"), "--out", str(out)])
    assert rc == 1
    assert not out.exists()


def test_convert_missing_required_column_is_exit_1(tmp_path: Path) -> None:
    src = tmp_path / "results.sage.tsv"
    src.write_text("peptide\tcharge\nPEPTIDEK\t2\n")
    out = tmp_path / "peptides.raw.tsv"
    assert cli.main(["convert", "--format", "sage", str(src), "--out", str(out)]) == 1


# ---------------------------------------------------------------- plan


def _sheet_row(dataset: str, n: int, gb_each: int) -> dict[str, Any]:
    names = [f"{dataset}_r{i:02d}.raw" for i in range(n)]
    return {
        "dataset": dataset,
        "tax_id": 511145,
        "organism": "Escherichia coli str. K-12 substr. MG1655",
        "below_species": "true",
        "title": "t",
        "submission_type": "COMPLETE",
        "acquisition": "DDA",
        "instruments": "Q Exactive",
        "pubmed_ids": "",
        "dois": "",
        "genome_accession": ACC,
        "genome_choice": "unique",
        "raw_format": "thermo_raw",
        "n_raw": n,
        "raw_gb": n * gb_each,
        "raw_files": ";".join(names),
        "raw_urls": ";".join(f"{PRIDE}/{dataset}/{x}" for x in names),
        "raw_bytes": ";".join(str(gb_each * 10**9) for _ in names),
        "raw_checksums": ";".join("" for _ in names),
        "formats_seen": "thermo_raw",
        "flags": "",
        "status": "ready",
    }


def test_plan_writes_jobs_and_prints_the_sbatch_line(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    sheet = tmp_path / "sheet.tsv"
    rows = [
        _sheet_row("PXD000002", 4, 10),  # 40 GB
        _sheet_row("PXD000003", 2, 10),  # 20 GB
        {**_sheet_row("PXD000004", 2, 10), "status": "needs_genome"},
    ]
    pd.DataFrame(rows).to_csv(sheet, sep="\t", index=False)
    jobs_tsv = tmp_path / "jobs.tsv"
    rc = cli.main(
        [
            "plan",
            "--sheet",
            str(sheet),
            "--root",
            str(tmp_path / "root"),
            "--out",
            str(jobs_tsv),
            "--disk-budget-gb",
            "300",
        ]
    )
    assert rc == 0
    jobs = _read_tsv(jobs_tsv)
    assert list(jobs.columns) == PLAN_COLUMNS
    assert list(jobs["dataset"]) == ["PXD000002", "PXD000003"]
    out = capsys.readouterr().out
    assert "sbatch --array=0-1%2" in out  # floor(300 / (3 * 40)) = 2


def test_plan_missing_sheet_is_exit_1(tmp_path: Path) -> None:
    rc = cli.main(
        [
            "plan",
            "--sheet",
            str(tmp_path / "missing.tsv"),
            "--root",
            str(tmp_path),
            "--out",
            str(tmp_path / "jobs.tsv"),
        ]
    )
    assert rc == 1


def test_plan_sheet_without_columns_is_exit_1(tmp_path: Path) -> None:
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("dataset\tstatus\nPXD000001\tready\n")
    rc = cli.main(
        ["plan", "--sheet", str(sheet), "--root", str(tmp_path), "--out", str(tmp_path / "j.tsv")]
    )
    assert rc == 1


# ---------------------------------------------------------------- run-job


def _jobs_tsv(tmp_path: Path) -> Path:
    rows = []
    for i, ds in enumerate(("PXD000002", "PXD000003")):
        rows.append(
            {
                "job": i,
                "dataset": ds,
                "tax_id": 511145,
                "organism": "Escherichia coli",
                "genome_accession": ACC,
                "raw_format": "thermo_raw",
                "acquisition_metadata": "DDA",
                "run_files": f"{ds}_r00.raw",
                "run_urls": f"{PRIDE}/{ds}/{ds}_r00.raw",
                "run_bytes": "1000",
                "run_checksums": "",
                "run_gb": 0.0,
                "denovo_runs": f"{ds}_r00.raw",
                "workdir": str(tmp_path / "root" / f"{ds}_511145"),
            }
        )
    path = tmp_path / "jobs.tsv"
    pd.DataFrame(rows, columns=PLAN_COLUMNS).to_csv(path, sep="\t", index=False)
    return path


class RunJobSpy:
    def __init__(self, code: int) -> None:
        self.code = code
        self.calls: list[tuple[Any, dict[str, Any]]] = []

    def __call__(self, job: Any, **kwargs: Any) -> int:
        self.calls.append((job, kwargs))
        return self.code


def _patch_run_job(monkeypatch: pytest.MonkeyPatch, spy: RunJobSpy) -> None:
    monkeypatch.setattr(driver, "run_job", spy)
    monkeypatch.setattr(cli, "run_job", spy, raising=False)


@pytest.mark.parametrize("code", [0, 1, 2])
def test_run_job_returns_the_job_code(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, code: int
) -> None:
    spy = RunJobSpy(code)
    _patch_run_job(monkeypatch, spy)
    db = tmp_path / "db"
    monkeypatch.setenv("GIAP_DB_DIR", str(db))
    monkeypatch.setenv("GIAP_ROOT", str(tmp_path / "giaproot"))
    rc = cli.main(["run-job", "--jobs", str(_jobs_tsv(tmp_path)), "--index", "1", "--db-dir", str(db)])
    assert rc == code
    ((job, kwargs),) = spy.calls
    assert job["dataset"] == "PXD000003"
    assert Path(str(job["workdir"])) == tmp_path / "root" / "PXD000003_511145"
    assert Path(kwargs["db_dir"]) == db
    assert kwargs["profile"] == "farmshare_local,apptainer"
    assert int(kwargs["threads"]) == 16
    assert kwargs["keep_inputs"] is False
    assert callable(kwargs["run"])
    assert callable(kwargs["download"])
    env = kwargs["env"]
    for key in ("db_dir", "diann_sif", "cascadia_model", "casanovo_model", "diann_cfg"):
        assert key in env, key
    assert Path(env["diann_cfg"]) == CAMPAIGN / "diann.cfg"
    assert env["diann_sif"] == f"{tmp_path / 'giaproot'}/images/diann_2.3.1.sif"
    assert str(env["cascadia_model"]).startswith(str(db))
    pins = env["pins"]
    pilot = (REPO / "validation/pilot/diann_image.sha256").read_text().split()[0]
    assert pins["diann_sif"] == pilot
    assert pins["casanovo_model"] == ""


def test_run_job_options(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    spy = RunJobSpy(0)
    _patch_run_job(monkeypatch, spy)
    rc = cli.main(
        [
            "run-job",
            "--jobs",
            str(_jobs_tsv(tmp_path)),
            "--index",
            "0",
            "--db-dir",
            str(tmp_path / "db"),
            "--profile",
            "test,docker",
            "--threads",
            "8",
            "--keep-inputs",
        ]
    )
    assert rc == 0
    ((job, kwargs),) = spy.calls
    assert job["dataset"] == "PXD000002"
    assert kwargs["profile"] == "test,docker"
    assert int(kwargs["threads"]) == 8
    assert kwargs["keep_inputs"] is True


def test_run_job_index_out_of_range_is_exit_1(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    spy = RunJobSpy(0)
    _patch_run_job(monkeypatch, spy)
    jobs = str(_jobs_tsv(tmp_path))
    rc = cli.main(["run-job", "--jobs", jobs, "--index", "2", "--db-dir", str(tmp_path)])
    assert rc == 1
    assert spy.calls == []


def test_run_job_missing_jobs_file_is_exit_1(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    spy = RunJobSpy(0)
    _patch_run_job(monkeypatch, spy)
    missing = str(tmp_path / "missing.tsv")
    assert cli.main(["run-job", "--jobs", missing, "--index", "0", "--db-dir", str(tmp_path)]) == 1
    assert spy.calls == []


# ---------------------------------------------------------------- aggregate


def test_aggregate_calls_aggregate(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[tuple[Any, Any]] = []

    def fake(jobs: Any, out_dir: Any) -> dict[str, int]:
        calls.append((jobs, out_dir))
        return {"campaign_runs.tsv": 4, "campaign_novel_orfs.tsv": 0}

    monkeypatch.setattr(aggregate, "aggregate", fake)
    if callable(getattr(cli, "aggregate", None)):  # `from .aggregate import aggregate` in cli
        monkeypatch.setattr(cli, "aggregate", fake)
    jobs_tsv = _jobs_tsv(tmp_path)
    out = tmp_path / "pooled"
    assert cli.main(["aggregate", "--jobs", str(jobs_tsv), "--out", str(out)]) == 0
    ((jobs, out_dir),) = calls
    assert Path(out_dir) == out
    if isinstance(jobs, pd.DataFrame):
        assert list(jobs["dataset"]) == ["PXD000002", "PXD000003"]
    else:
        assert Path(jobs) == jobs_tsv


def test_aggregate_missing_jobs_file_is_exit_1(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(aggregate, "aggregate", lambda jobs, out_dir: {})
    rc = cli.main(["aggregate", "--jobs", str(tmp_path / "missing.tsv"), "--out", str(tmp_path)])
    assert rc == 1


# ---------------------------------------------------------------- catalog

PROJECT = {
    "accession": "PXD000001",
    "title": "Escherichia coli K-12 proteome, data-dependent acquisition",
    "submissionType": "COMPLETE",
    "organisms": [{"accession": "NEWT:511145", "name": "Escherichia coli K-12"}],
    "instruments": [{"name": "Q Exactive HF"}],
    "experimentTypes": ["Shotgun proteomics"],
    "softwares": [{"name": "MaxQuant"}],
    "keywords": ["DDA"],
    "references": [{"pubmedID": 12345678, "doi": "10.1000/xyz"}],
}
FILES = [
    {
        "fileName": name,
        "fileCategory": {"value": "RAW"},
        "fileSizeBytes": 2_000_000_000,
        "publicFileLocations": [
            {"name": "FTP Protocol", "value": f"ftp://ftp.pride.ebi.ac.uk/pride/data/{name}"},
            {"name": "Aspera Protocol", "value": f"prd_ascp@fasp.ebi.ac.uk:pride/data/{name}"},
        ],
        "checksum": "d41d8cd98f00b204e9800998ecf8427e",
    }
    for name in ("WT_1.raw", "WT_2.raw")
]
TAXON = {
    "taxonomy_nodes": [
        {
            "taxonomy": {
                "tax_id": 511145,
                "organism_name": "Escherichia coli str. K-12 substr. MG1655",
                "rank": "strain",
                "lineage": [1, 131567, 2, 1224, 1236, 91347, 543, 561, 562, 83333],
                "classification": {
                    "domain": {"name": "Bacteria", "id": 2},
                    "species": {"name": "Escherichia coli", "id": 562},
                },
            }
        }
    ]
}
ASSEMBLIES = {
    "reports": [
        {
            "accession": ACC,
            "organism": {
                "tax_id": 511145,
                "organism_name": "Escherichia coli str. K-12 substr. MG1655",
                "infraspecific_names": {"strain": "K-12 substr. MG1655"},
            },
            "assembly_info": {
                "assembly_level": "Complete Genome",
                "refseq_category": "reference genome",
                "release_date": "2013-09-26",
            },
        }
    ]
}


def fake_get_json(url: str, *_: Any, **__: Any) -> Any:
    assert url.startswith("https://")
    if "/taxonomy/taxon/511145" in url:
        return TAXON
    if "/genome/taxon/511145/dataset_report" in url:
        return ASSEMBLIES
    if url.rstrip("/").endswith("/projects/PXD000001/files/all"):
        return FILES
    if url.rstrip("/").endswith("/projects/PXD000001"):
        return PROJECT
    if "/projects?" in url:
        return [PROJECT] if "page=0" in url else []
    raise LookupError(url)


def _patch_get_json(monkeypatch: pytest.MonkeyPatch, fake: Any) -> None:
    monkeypatch.setattr(http, "get_json", fake)
    monkeypatch.setattr(cli, "get_json", fake, raising=False)


def test_catalog_writes_the_sheet(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _patch_get_json(monkeypatch, fake_get_json)
    out = tmp_path / "sheet.tsv"
    rc = cli.main(
        [
            "catalog",
            "--out",
            str(out),
            "--accessions",
            "PXD000001",
            "--cache-dir",
            str(tmp_path / "cache"),
        ]
    )
    assert rc == 0
    sheet = _read_tsv(out)
    assert list(sheet["dataset"]) == ["PXD000001"]
    row = sheet.iloc[0]
    assert row["tax_id"] == "511145"
    assert row["genome_accession"] == ACC
    assert row["raw_format"] == "thermo_raw"
    assert row["status"] == "ready"


@pytest.mark.parametrize(
    "extra", [["--accessions", "PXD000001"], ["--max-projects", "5"]], ids=["accessions", "paging"]
)
def test_catalog_network_error_is_exit_2(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, extra: list[str]
) -> None:
    def down(url: str, *_: Any, **__: Any) -> Any:
        raise OSError(f"network unreachable: {url}")

    _patch_get_json(monkeypatch, down)
    assert cli.main(["catalog", "--out", str(tmp_path / "sheet.tsv"), *extra]) == 2
