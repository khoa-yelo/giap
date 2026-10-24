"""The R11 audit over ``git ls-files`` (docs/specs/M10-14-pilot-reproducibility.md, "The R11
audit", acceptance criterion 5; GIAP_SPEC.md §0.2 R11; Khoa's data rule as the M10-14 lead
decisions restate it: no row-level DIA-NN content is committed).

It consolidates the per-issue guards (``test_m10_01_committed.py``,
``test_m10_02_recorded_diann.py``, …) without replacing them. Each of the six rules is a pure
function over ``(path, size, head)`` triples; it passes on the tracked tree and fails on a
fixture list built to break it.

Allow-lists pinned by the test author from the tree of 2026-10-24 (main cb33f03f):

* ``.parquet`` only under ``workflow/tests/data/`` (synthetic nf-test data);
* ``changelog.d/`` is a towncrier directory, not a Bruker ``.d`` run;
* ``validation/pilot/metrics/`` may hold only ``a3_feasibility.json`` (the lead extends
  ``METRICS_ALLOWED`` if Khoa's Q-K1 answer allows more);
* ``regulondb`` in a path only for ``*.py`` and ``*.md``; no ``kegg`` file under
  ``tests/data/`` or ``validation/``.

No statistical claim: no calibration test (spec §8.4).
"""

from __future__ import annotations

import fnmatch
import subprocess
from collections.abc import Iterable
from pathlib import PurePosixPath

import pytest
from helpers.slurm_sandbox import REPO

MAX_FILE = 5 * 1024 * 1024
MAX_TESTS_DATA = 50 * 1024 * 1024
ELF = b"\x7fELF"

BANNED_SUFFIXES = (
    ".sif",
    ".speclib",
    ".pt",
    ".ckpt",
    ".mzml",
    ".raw",
    ".tdf",
    ".tdf_bin",
    ".zip",
    ".ssl",
    ".mztab",
)
D_DIR_ALLOWED = frozenset({"changelog.d"})
PARQUET_ALLOWED = ("workflow/tests/data/",)
PARQUET_FORBIDDEN = ("validation/", "tests/data/recorded/")
METRICS_DIR = "validation/pilot/metrics/"
METRICS_ALLOWED = frozenset({"a3_feasibility.json"})

DIANN_NAMES = (
    "a[1-3]*.tsv*",
    "*_abundance.tsv*",
    "*.convert.json",
    "*.convert.log",
    "report.sample.*",
    "te_real_exceptions.tsv",
    "te7_files.tsv",
    "entrapment_discoveries.tsv",
    "ATTEMPTS.tsv",
    "report.html",
    "compare.html",
)
RECORDED_DIANN = "tests/data/recorded/diann/"
RECORDED_DIANN_NAMES = ("summary.json", "schema.tsv")

Entry = tuple[str, int, bytes]


def tracked() -> list[Entry]:
    out = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-z"], capture_output=True, check=True
    ).stdout.decode()
    entries = []
    for path in filter(None, out.split("\0")):
        p = REPO / path
        if p.is_symlink() or not p.is_file():
            continue  # deleted in the working tree, or a link (its target is audited)
        with p.open("rb") as fh:
            head = fh.read(4)
        entries.append((path, p.stat().st_size, head))
    return entries


@pytest.fixture(scope="module")
def tree() -> list[Entry]:
    return tracked()


# ---------------------------------------------------------------- the six rules


def rule1_sizes(entries: Iterable[Entry]) -> list[str]:
    entries = list(entries)
    bad = [f"too_big:{p}:{n}" for p, n, _ in entries if n > MAX_FILE]
    total = sum(n for p, n, _ in entries if p.startswith("tests/data/"))
    if total >= MAX_TESTS_DATA:
        bad.append(f"tests_data_total:{total}")
    return bad


def rule2_banned(entries: Iterable[Entry]) -> list[str]:
    bad = []
    for path, _, head in entries:
        parts = PurePosixPath(path).parts
        if parts[0] == "local":
            bad.append(f"local:{path}")
        if path.lower().endswith(BANNED_SUFFIXES):
            bad.append(f"suffix:{path}")
        if any(d.endswith(".d") and d not in D_DIR_ALLOWED for d in parts[:-1]):
            bad.append(f"d_dir:{path}")
        if head.startswith(ELF):
            bad.append(f"elf:{path}")
    return bad


def rule3_parquet(entries: Iterable[Entry]) -> list[str]:
    bad = []
    for path, _, _ in entries:
        if not path.lower().endswith(".parquet"):
            continue
        if path.startswith(PARQUET_FORBIDDEN) or not path.startswith(PARQUET_ALLOWED):
            bad.append(f"parquet:{path}")
    return bad


def rule4_diann_names(entries: Iterable[Entry]) -> list[str]:
    bad = []
    for path, _, _ in entries:
        name = PurePosixPath(path).name
        if any(fnmatch.fnmatchcase(name, pat) for pat in DIANN_NAMES):
            bad.append(f"diann_name:{path}")
        if path.startswith(RECORDED_DIANN) and name in RECORDED_DIANN_NAMES:
            bad.append(f"diann_recorded:{path}")
        if fnmatch.fnmatchcase(path, "*evidence/*.json"):
            bad.append(f"evidence:{path}")
    return bad


def rule5_metrics(entries: Iterable[Entry]) -> list[str]:
    return [
        f"metrics:{p}"
        for p, _, _ in entries
        if p.startswith(METRICS_DIR) and p[len(METRICS_DIR) :] not in METRICS_ALLOWED
    ]


def rule6_regulondb_kegg(entries: Iterable[Entry]) -> list[str]:
    bad = []
    for path, _, _ in entries:
        name = PurePosixPath(path).name.lower()
        if "regulondb" in path.lower() and not name.endswith((".py", ".md")):
            bad.append(f"regulondb:{path}")
        if "kegg" in path.lower() and path.startswith(("tests/data/", "validation/")):
            bad.append(f"kegg:{path}")
    return bad


RULES = {
    "1_sizes": rule1_sizes,
    "2_banned": rule2_banned,
    "3_parquet": rule3_parquet,
    "4_diann_names": rule4_diann_names,
    "5_metrics": rule5_metrics,
    "6_regulondb_kegg": rule6_regulondb_kegg,
}


@pytest.mark.parametrize("rule", sorted(RULES))
def test_tree_passes(tree: list[Entry], rule: str) -> None:
    assert tree, "git ls-files returned nothing"
    assert RULES[rule](tree) == []


# ---------------------------------------------------------------- each rule fails on its fixture

OK: list[Entry] = [
    ("src/giap/cli.py", 1000, b'"""g'),
    ("tests/data/genomes/nissle1917.fna.gz", 1_506_323, b"\x1f\x8b\x08\x00"),
    ("changelog.d/m0-scaffold.feature.md", 50, b"Scaf"),
    ("workflow/tests/data/mapping/peptides.parquet", 2000, b"PAR1"),
    ("validation/pilot/metrics/a3_feasibility.json", 300, b"{\n  "),
    ("validation/fetch_regulondb.py", 4000, b'"""F'),
    ("docs/adr/ADR-0023-kegg-overview-maps-excluded.md", 900, b"# AD"),
    ("src/giap/report/templates/report.html.j2", 900, b"<!DO"),
]

BREAKERS: dict[str, list[Entry]] = {
    "1_sizes": [
        ("tests/data/recorded/big.json.gz", MAX_FILE + 1, b"\x1f\x8b\x08\x00"),
        *[(f"tests/data/chunk{i}.bin", 4 * 1024 * 1024, b"\0\0\0\0") for i in range(13)],
    ],
    "2_banned": [
        ("local/data/x.tsv", 10, b"a\tb\n"),
        ("validation/pilot/diann_2.3.1.sif", 10, b"hsqs"),
        ("tests/data/A1.speclib", 10, b"\0\0\0\0"),
        ("tests/data/model.ckpt", 10, b"PK\x03\x04"),
        ("tests/data/cascadia.pt", 10, b"PK\x03\x04"),
        ("tests/data/WT_1.mzML", 10, b"<?xm"),
        ("tests/data/run.raw", 10, b"\x01\xa1F\x00"),
        ("tests/data/WT_1.d/analysis.tdf", 10, b"SQLi"),
        ("tests/data/WT_1.d/analysis.tdf_bin", 10, b"\0\0\0\0"),
        ("tests/data/WT_1.zip", 10, b"PK\x03\x04"),
        ("tests/data/WT_1.ssl", 10, b"file"),
        ("tests/data/WT_1.mzTab", 10, b"MTD\t"),
        ("tests/data/runs/WT_2.d/chromatography-data.sqlite", 10, b"SQLi"),
        ("scripts/diann-linux", 10, ELF + b"\x02"),
    ],
    "3_parquet": [
        ("validation/pilot/a1_sample.parquet", 10, b"PAR1"),
        ("tests/data/recorded/diann/2.3.1/report.parquet", 10, b"PAR1"),
        ("tests/data/other/x.parquet", 10, b"PAR1"),
    ],
    "4_diann_names": [
        ("validation/pilot/inputs/a1.tsv.gz", 10, b"\x1f\x8b\x08\x00"),
        ("validation/pilot/a2_peptides.tsv", 10, b"pep\t"),
        ("validation/pilot/A2_abundance.tsv.gz", 10, b"\x1f\x8b\x08\x00"),
        ("validation/pilot/inputs/a1.convert.json", 10, b"{\n  "),
        ("validation/pilot/inputs/a1.convert.log", 10, b"conv"),
        ("tests/data/recorded/diann/2.3.1/A1/report.sample.parquet", 10, b"PAR1"),
        ("tests/data/recorded/diann/2.3.1/A1/summary.json", 10, b"{\n  "),
        ("tests/data/recorded/diann/2.3.1/A1/schema.tsv", 10, b"colu"),
        ("validation/pilot/te_real_exceptions.tsv", 10, b"file"),
        ("validation/pilot/te7_files.tsv", 10, b"file"),
        ("validation/pilot/entrapment_discoveries.tsv", 10, b"db_i"),
        ("validation/pilot/ATTEMPTS.tsv", 10, b"arm\t"),
        ("validation/pilot/evidence/A2.json", 10, b"{\n  "),
        ("validation/pilot/report.html", 10, b"<!DO"),
        ("validation/pilot/compare.html", 10, b"<!DO"),
    ],
    "5_metrics": [
        ("validation/pilot/metrics/metrics.json", 10, b"{\n  "),
        ("validation/pilot/metrics/freeze-0123456/metrics.json", 10, b"{\n  "),
    ],
    "6_regulondb_kegg": [
        ("tests/data/regulondb/TUSet.txt", 10, b"# Re"),
        ("validation/regulondb_operons.tsv", 10, b"oper"),
        ("tests/data/kegg/ko00010.json", 10, b"{\n  "),
        ("validation/pilot/kegg_pathways.tsv", 10, b"path"),
    ],
}


@pytest.mark.parametrize("rule", sorted(RULES))
def test_rule_accepts_the_allowed_fixture(rule: str) -> None:
    assert RULES[rule](OK) == [], RULES[rule](OK)


@pytest.mark.parametrize("rule", sorted(BREAKERS))
def test_each_breaker_fails_its_rule(rule: str) -> None:
    fixture = BREAKERS[rule]
    if rule == "1_sizes":
        assert any(s.startswith("too_big:") for s in rule1_sizes(fixture[:1]))
        assert any(s.startswith("tests_data_total:") for s in rule1_sizes(fixture[1:]))
        return
    for entry in fixture:
        assert RULES[rule]([entry]), f"rule {rule} does not catch {entry[0]}"
