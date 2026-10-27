"""The frozen golden inputs generator (M8-12 "What is added" item 1, examples 3; §11.5, R11).

Interface assumed (the spec leaves it open):

``python scripts/golden/make_golden_inputs.py [--out DIR] [--seed INT]``

- ``--out`` defaults to ``tests/data/golden/r/inputs``; ``--seed`` defaults to ``GOLDEN_SEED``
  read from ``tests/calibration/tolerances.yaml`` (50000).
- It writes four gzip TSVs, each with one header row, tab-separated, ``\\n`` line ends, no
  index column, written so the bytes depend only on the seed (gzip mtime 0, no file name):

  ============================  ==========================================================
  ``gsea_stats.tsv.gz``         ``tu_id  score``          2,000 rows, TU-ID order
  ``gsea_sets.tsv.gz``          ``set_id  member``        long table, 120 sets
  ``ora_genes.tsv.gz``          ``gene  p_detect  detected``  3,000 rows (detected is 0/1)
  ``ora_categories.tsv.gz``     ``category  gene``        long table, 150 categories
  ============================  ==========================================================

- numpy only (plus the standard library and PyYAML to read the seed); never ``gippa``.
"""

from __future__ import annotations

import ast
import csv
import gzip
import io
import math
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import hypergeom

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "golden" / "make_golden_inputs.py"
FILES = {
    "gsea_stats.tsv.gz": ["tu_id", "score"],
    "gsea_sets.tsv.gz": ["set_id", "member"],
    "ora_genes.tsv.gz": ["gene", "p_detect", "detected"],
    "ora_categories.tsv.gz": ["category", "gene"],
}
FIVE_MB = 5 * 1024 * 1024
# KEGG (K00001, map00010, ko00010, path:, ec:), GO (GO:0008150) and COG (COG0001) identifiers.
CURATED_ID = re.compile(
    r"(GO:\d{7}|\bK\d{5}\b|\b(map|ko|eco)\d{5}\b|path:|\bec:|COG\d{4}|\bEC \d)", re.I
)


def generate(out: Path, *seed: str) -> subprocess.CompletedProcess[str]:
    assert SCRIPT.is_file(), f"missing {SCRIPT.relative_to(REPO)}"
    args = [sys.executable, str(SCRIPT), "--out", str(out)]
    if seed:
        args += ["--seed", *seed]
    proc = subprocess.run(args, capture_output=True, text=True, cwd=REPO)
    assert proc.returncode == 0, proc.stderr
    return proc


def read_tsv(path: Path) -> tuple[list[str], list[list[str]]]:
    with gzip.open(path, "rt", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return rows[0], rows[1:]


@pytest.fixture(scope="module")
def inputs(tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("golden_inputs")
    generate(out)
    return out


@pytest.fixture(scope="module")
def gsea(inputs: Path) -> tuple[list[str], np.ndarray, dict[str, list[str]]]:
    _, rows = read_tsv(inputs / "gsea_stats.tsv.gz")
    ids = [r[0] for r in rows]
    scores = np.array([float(r[1]) for r in rows])
    sets: dict[str, list[str]] = defaultdict(list)
    for set_id, member in read_tsv(inputs / "gsea_sets.tsv.gz")[1]:
        sets[set_id].append(member)
    return ids, scores, dict(sets)


@pytest.fixture(scope="module")
def ora(inputs: Path) -> tuple[list[str], np.ndarray, np.ndarray, dict[str, list[str]]]:
    _, rows = read_tsv(inputs / "ora_genes.tsv.gz")
    genes = [r[0] for r in rows]
    p = np.array([float(r[1]) for r in rows])
    det = np.array([int(r[2]) for r in rows])
    cats: dict[str, list[str]] = defaultdict(list)
    for cat, gene in read_tsv(inputs / "ora_categories.tsv.gz")[1]:
        cats[cat].append(gene)
    return genes, p, det, dict(cats)


# --- determinism (example 3) ---------------------------------------------------------------


def test_two_runs_give_byte_identical_inputs(inputs: Path, tmp_path: Path) -> None:
    generate(tmp_path)
    for name in FILES:
        assert (tmp_path / name).read_bytes() == (inputs / name).read_bytes(), name


def test_default_seed_is_golden_seed(inputs: Path, tmp_path: Path, golden_seed: int) -> None:
    generate(tmp_path, str(golden_seed))
    for name in FILES:
        assert (tmp_path / name).read_bytes() == (inputs / name).read_bytes(), name


def test_another_seed_gives_other_bytes(inputs: Path, tmp_path: Path, golden_seed: int) -> None:
    generate(tmp_path, str(golden_seed + 1))
    assert (tmp_path / "gsea_stats.tsv.gz").read_bytes() != (
        inputs / "gsea_stats.tsv.gz"
    ).read_bytes()


def test_writes_exactly_the_four_files(inputs: Path) -> None:
    assert sorted(p.name for p in inputs.iterdir()) == sorted(FILES)


# --- format and size (R11) -----------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(FILES))
def test_file_is_gzip_tsv_with_header(inputs: Path, name: str) -> None:
    path = inputs / name
    raw = path.read_bytes()
    assert raw[:2] == b"\x1f\x8b", "gzip magic"
    assert len(raw) <= FIVE_MB
    text = gzip.decompress(raw).decode("utf-8")
    assert "\r" not in text
    assert text.endswith("\n")
    header, rows = read_tsv(path)
    assert header == FILES[name]
    assert rows
    assert all(len(r) == len(header) for r in rows)


def test_numpy_only_and_no_gippa_import() -> None:
    assert SCRIPT.is_file(), f"missing {SCRIPT.relative_to(REPO)}"
    tree = ast.parse(SCRIPT.read_text())
    allowed = set(sys.stdlib_module_names) | {"numpy", "yaml", "__future__"}
    for node in ast.walk(tree):
        names: list[str] = []
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module:
            names = [node.module]
        for name in names:
            assert name.split(".")[0] in allowed, f"import of {name}"


# --- GSEA inputs ---------------------------------------------------------------------------


def test_gsea_scores_are_2000_tus_in_id_order(gsea: tuple) -> None:  # type: ignore[type-arg]
    ids, scores, _ = gsea
    assert len(ids) == 2000
    assert ids == [f"TU_{i:04d}" for i in range(1, 2001)]
    assert np.isfinite(scores).all()


def test_gsea_scores_have_ties_zeros_and_are_centred(gsea: tuple) -> None:  # type: ignore[type-arg]
    _, scores, _ = gsea
    values, counts = np.unique(scores, return_counts=True)
    n_zero = int((scores == 0.0).sum())
    assert 2 <= n_zero <= 20, "a few zeros"
    tied = [(v, c) for v, c in zip(values, counts, strict=True) if c > 1 and v != 0.0]
    assert len(tied) == 2, "two repeated (non-zero) values"
    frac = sum(c for _, c in tied) / scores.size
    assert 0.04 <= frac <= 0.06, f"5 % exact ties, got {frac:.3f}"
    assert abs(float(np.median(scores))) <= 0.05, "already (median-)centred"
    assert (scores > 0).any() and (scores < 0).any()


def test_gsea_has_120_sets_of_legal_sizes(gsea: tuple) -> None:  # type: ignore[type-arg]
    _, _, sets = gsea
    assert len(sets) == 120
    for set_id, members in sets.items():
        assert len(members) == len(set(members)), f"duplicate member in {set_id}"
        assert 10 <= len(members) <= 500, set_id
    assert sum(len(m) == 10 for m in sets.values()) >= 5, "5 tiny-but-legal sets"
    sizes = np.array([len(m) for m in sets.values()])
    assert sizes.max() >= 200, "log-uniform sizes reach the upper range"


def test_exactly_two_sets_contain_absent_names(gsea: tuple) -> None:  # type: ignore[type-arg]
    ids, _, sets = gsea
    known = set(ids)
    with_absent = [s for s, m in sets.items() if any(x not in known for x in m)]
    assert len(with_absent) == 2, with_absent


def test_planted_sets_span_strong_signals_in_both_directions(gsea: tuple) -> None:  # type: ignore[type-arg]
    ids, scores, sets = gsea
    by_id = dict(zip(ids, scores, strict=True))
    sd = float(scores.std())
    z = []
    for members in sets.values():
        vals = np.array([by_id[m] for m in members if m in by_id])
        z.append(vals.mean() / sd * math.sqrt(vals.size))
    z_arr = np.array(z)
    strong = np.abs(z_arr) >= 2.0  # the weakest planted set has true p ~1e-2 (z ~2.3)
    assert strong.sum() >= 20, "20 planted sets"
    assert (z_arr >= 3).any() and (z_arr <= -3).any(), "planted toward top and bottom"
    assert np.abs(z_arr).max() >= 9, "true p reaches ~1e-30"


# --- ORA inputs ----------------------------------------------------------------------------


def test_ora_universe_and_detection(ora: tuple) -> None:  # type: ignore[type-arg]
    genes, p, det, _ = ora
    assert len(genes) == 3000 == len(set(genes))
    assert ((p >= 0) & (p <= 1)).all()
    assert (p == 0).sum() >= 1 and (p == 1).sum() >= 1, "a few p_detect at exactly 0 and 1"
    assert set(np.unique(det)) <= {0, 1}
    assert (det[p == 0] == 0).all() and (det[p == 1] == 1).all(), "Bernoulli(p_detect)"
    assert abs(float(p.mean()) - 0.5) <= 0.05, "Beta(2, 2)"
    assert abs(float(det.mean()) - float(p.mean())) <= 0.05


def test_ora_categories_sizes_and_special_cases(ora: tuple) -> None:  # type: ignore[type-arg]
    genes, _, det, cats = ora
    universe = set(genes)
    idx = {g: i for i, g in enumerate(genes)}
    assert len(cats) == 150
    categorised = set().union(*map(set, cats.values()))
    # PR #192 B1 (lead): "all categorised genes minus one", so the uncategorised 30 % stay out.
    near_universe = [c for c, m in cats.items() if len(set(m)) == len(categorised) - 1]
    assert len(near_universe) == 1, "one category equal to the categorised genes minus one"
    k_eq_n = k_zero = 0
    for cat, members in cats.items():
        assert len(members) == len(set(members)), f"duplicate gene in {cat}"
        assert set(members) <= universe, f"{cat} has genes outside the universe"
        if cat not in near_universe:
            assert 5 <= len(members) <= 500, cat
        k = int(det[[idx[g] for g in members]].sum())
        k_eq_n += k == len(members)
        k_zero += k == 0
    assert k_eq_n >= 1, "a category with k = n"
    assert k_zero >= 1, "a category with k = 0"


def test_ora_has_genes_without_category(ora: tuple) -> None:  # type: ignore[type-arg]
    # PR #192 B1: 30 % of genes are in no category at all, so use_genes_without_cat matters.
    genes, _, _, cats = ora
    covered = set().union(*map(set, cats.values()))
    frac = 1 - len(covered) / len(genes)
    assert 0.25 <= frac <= 0.35, frac


def test_ora_weights_reach_both_extremes(ora: tuple) -> None:  # type: ignore[type-arg]
    genes, p, _, cats = ora
    idx = {g: i for i, g in enumerate(genes)}
    omegas = []
    for members in cats.values():
        mask = np.zeros(len(genes), dtype=bool)
        mask[[idx[g] for g in members]] = True
        out = p[~mask].mean()
        if out > 0:
            omegas.append(p[mask].mean() / out)  # M8-01: mean in / mean out
    assert max(omegas) >= 1.5, "ω ≫ 1"
    assert min(omegas) <= 0.5, "ω ≪ 1"


def test_ora_planted_enrichment_is_strong(ora: tuple) -> None:  # type: ignore[type-arg]
    genes, _, det, cats = ora
    idx = {g: i for i, g in enumerate(genes)}
    big_m, big_n = len(genes), int(det.sum())
    best = min(
        hypergeom(big_m, len(m), big_n).logsf(int(det[[idx[g] for g in m]].sum()) - 1)
        for m in cats.values()
    )
    assert best <= math.log(1e-30), "planted enrichments down to ~1e-40"


# --- synthetic content only (index §KEGG) --------------------------------------------------


def test_inputs_carry_no_kegg_go_or_cog_identifiers(inputs: Path) -> None:
    for name in FILES:
        text = gzip.decompress((inputs / name).read_bytes()).decode()
        for row in csv.reader(io.StringIO(text), delimiter="\t"):
            for field in row:
                assert not CURATED_ID.search(field), f"{name}: {field!r}"
