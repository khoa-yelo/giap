"""Driver and checks shared by the M5-11 tests (docs/specs/M5-11-integration-golden.md).

The M5 chain is the M4 chain of ``m4_chain`` (``run_gene_table`` -> ``run_detectability`` ->
``run_context_stats`` -> ``run_sensitivity``): since M5-09 ``context_stats`` also writes
``tables/autocorrelation.tsv``, ``tables/islands.tsv``, the genes.tsv context columns
(``island_id``, ``replichore``, ``replichore_position``) and
``context/{gene_context,island_diagnostics}.parquet``. This module adds the two M5 fixtures and
the golden comparison helpers.

Fixtures (spec "Interface": ``build_m4_fixture(..., oric_frac: float | None = 0.8,
islands: Sequence[IslandSpec] = ())``):

* ``sim50k`` -- unchanged M4-20 call ``build_m4_fixture(work, seed=, size="sim50k")``; its
  islands golden is header-only (a 50-kb genome has too few TUs; lead decision 39 keeps it to test
  the empty path).
* the island golden -- ``ISLAND`` below: one chromosome with an oriC at ``ORIC_FRAC`` and one
  planted high island of ``n_tus`` TUs shifted by ``effect_sd`` stationary SDs, starting at
  ``gippasim.scenarios.island_first_tu(n_genes)``.

**Which island case (lead decision 40).** The spec and
decision 39 name ``sim1000_island``: 1,000 genes, +3 SD, 20 TUs, gippasim default detection
(``SPEC_ISLAND``). Checked against l' before freezing (decision 39) with the merged M5 code on a
prototype of the fixture (default config, B = 1,000), that case reports **no island at all**:

* gippasim default detection (``covariates_only``, about 45 % of genes censored independently of
  abundance) -- no island at +3, +5 or +10 SD on fixture seeds 1-4 or ``GOLDEN_SEED`` (l-hat 1,
  l' 3): the censored genes score lowest whatever their abundance, so the island's TU means stay
  near the genome's.
* ``detection = "abundance"`` on 1,000 genes -- at ``GOLDEN_SEED`` (50000) l-hat 4, l' 9 and the
  island's FWER p is 0.15 at +3 SD and floors at 0.035-0.041 from +5 SD up (rank scores
  saturate), i.e. found at alpha 0.05 only marginally and not "at any seed".
* ``detection = "abundance"`` on 4,500 genes, +3 SD, 20 TUs -- FWER p 0.004-0.008 at config seeds
  1-3 and ``GOLDEN_SEED`` and 0.003-0.005 at fixture seeds 1-3 (l' 9-13 < 20 TUs); exactly one
  island, ``high``, covering all 20 planted TUs.

So the tests use ``ISLAND = PROPOSED_ISLAND`` (``sim4500_island``) and keep ``SPEC_ISLAND``
here for the lead: switching back is a one-line change of ``ISLAND``, and the golden directory
follows ``ISLAND.name``.

``expression=`` is the one builder keyword beyond the spec's Interface: ``ExpressionConfig`` fields
passed through, as ``layout=`` passes ``LayoutConfig`` fields (M4-20 lead decision 4); the island
case needs ``detection`` and the matching ``beta``.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m4_chain as mc

REPO = Path(__file__).resolve().parents[2]
GOLDEN_ROOT = REPO / "tests" / "data" / "golden"
#: the oriC of the island fixture (the builder's default in the spec's Interface)
ORIC_FRAC = 0.8
#: the genes.tsv columns M5 fills (M5-09); the golden keeps them with ``locus_tag``
CONTEXT_COLUMNS: tuple[str, ...] = ("island_id", "replichore", "replichore_position")
#: §9.1 tables M5 adds, file name -> schema stem
M5_TABLES: Mapping[str, str] = {
    "autocorrelation.tsv": "autocorrelation",
    "islands.tsv": "islands",
}
#: every ``tables/`` file of the chain, file name -> schema stem
ALL_TABLES: Mapping[str, str] = {**mc.OUTPUT_TABLES, **M5_TABLES}
#: ``context_stats`` intermediates M5-09 writes (TE7 covers them too)
CONTEXT_PARQUETS: tuple[str, ...] = ("gene_context.parquet", "island_diagnostics.parquet")
#: ADR-0021 decision 1: l' = min(L_MAX, 2 l-hat + 1, T)
L_MAX = 20
#: the 50-kb golden runs at alpha.identification 0.05 like its M4 golden (M4-20 decision 1)
SIM50K_ALPHA = 0.05
#: R11: committed golden data of one case stays under 1 MB
MAX_GOLDEN_BYTES = 1_000_000


@dataclass(frozen=True)
class IslandCase:
    """One planted-island golden fixture."""

    name: str
    size: str  # build_m4_fixture size
    n_genes: int
    n_tus: int
    effect_sd: float
    detection: str  # gippasim ExpressionConfig.detection

    @property
    def first_tu(self) -> int:
        from gippasim.scenarios import island_first_tu

        return int(island_first_tu(self.n_genes))

    def island(self) -> Any:
        from gippasim.expression import IslandSpec

        return IslandSpec("chr1", self.first_tu, self.n_tus, self.effect_sd)

    def expression(self) -> dict[str, Any] | None:
        """``expression=`` for the builder (None: gippasim defaults)."""
        if self.detection == "covariates_only":
            return None
        from gippasim.scenarios import ABUNDANCE_BETA

        return {"detection": self.detection, "beta": dict(ABUNDANCE_BETA)}


#: the spec's case (decision 39) -- finds no island, see the module docstring
SPEC_ISLAND = IslandCase("sim1000_island", "sim1000", 1000, 20, 3.0, "covariates_only")
#: the test author's proposal, measured to be found at every seed tried
PROPOSED_ISLAND = IslandCase("sim4500_island", "sim4500", 4500, 20, 3.0, "abundance")
#: the case the tests and the golden use (lead decision 40)
ISLAND = PROPOSED_ISLAND


# --------------------------------------------------------------------------- building, running


def build_sim50k(work: Path, seed: int) -> Path:
    """The M4-20 golden fixture, unchanged call."""
    return mc.build(work, seed=seed, size="sim50k")


def build_island(work: Path, seed: int, case: IslandCase = ISLAND) -> Path:
    """``build_m4_fixture`` with the oriC and the planted island of ``case``."""
    from gippasim.fixtures import build_m4_fixture

    kw: dict[str, Any] = {"oric_frac": ORIC_FRAC, "islands": (case.island(),)}
    expr = case.expression()
    if expr is not None:
        kw["expression"] = expr
    work.mkdir(parents=True, exist_ok=True)
    build_m4_fixture(work, seed=seed, size=case.size, **kw)
    return work


def sim50k_config(seed: int) -> Any:
    """The M4 golden config: shipped defaults, ``alpha.identification`` = 0.05."""
    return mc.config(seed, stats={"alpha": {"identification": SIM50K_ALPHA}})


def island_config(seed: int) -> Any:
    """The shipped defaults (B = 1,000, ``block_len: auto``) with ``seed``."""
    return mc.config(seed)


def run_m5(work: Path, cfg: Any, db: Path) -> dict[str, Path]:
    """The chain in spec order; returns every ``tables/`` path by file name."""
    mc.run_m4(work, cfg, db)
    return {name: work / "tables" / name for name in ALL_TABLES}


def digests(work: Path) -> dict[str, str]:
    """SHA-256 of every file under ``tables/`` and of the two ``context/`` parquets (TE7)."""
    files = sorted((work / "tables").glob("*")) + [work / "context" / n for n in CONTEXT_PARQUETS]
    return {
        str(p.relative_to(work)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in files
        if p.is_file()
    }


# --------------------------------------------------------------------------- golden files


def golden_dir(case: str) -> Path:
    """``tests/data/golden/<case>/m5/`` (the M4 golden of sim50k lives in ``.../sim50k/m4/``)."""
    return GOLDEN_ROOT / case / "m5"


def golden_names(case: str) -> dict[str, str]:
    """logical name -> committed file name; the island case's gene columns are gzipped (R11)."""
    gc = "genes_context.tsv" if case == "sim50k" else "genes_context.tsv.gz"
    return {"autocorrelation.tsv": "autocorrelation.tsv", "islands.tsv": "islands.tsv",
            "genes_context.tsv": gc}  # fmt: skip


def select_columns(tsv: Path, columns: Sequence[str]) -> bytes:
    """The given columns of a written TSV, cell text untouched (byte-faithful projection)."""
    lines = tsv.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    missing = [c for c in columns if c not in header]
    assert not missing, f"{tsv.name} has no column(s) {missing}"
    idx = [header.index(c) for c in columns]
    out = ["\t".join(row[i] for i in idx) for row in (ln.split("\t") for ln in lines)]
    return ("\n".join(out) + "\n").encode("utf-8")


def drop_columns(tsv: Path, columns: Sequence[str]) -> bytes:
    """A written TSV without ``columns`` (cell text untouched)."""
    header = tsv.read_text(encoding="utf-8").splitlines()[0].split("\t")
    return select_columns(tsv, [c for c in header if c not in columns])


def genes_context_bytes(genes_tsv: Path) -> bytes:
    """``locus_tag`` and the three M5 genes.tsv columns."""
    return select_columns(genes_tsv, ("locus_tag", *CONTEXT_COLUMNS))


def produced(tables: Mapping[str, Path]) -> dict[str, bytes]:
    """The run's bytes for each golden logical name."""
    return {
        "autocorrelation.tsv": tables["autocorrelation.tsv"].read_bytes(),
        "islands.tsv": tables["islands.tsv"].read_bytes(),
        "genes_context.tsv": genes_context_bytes(tables["genes.tsv"]),
    }


def read_golden(path: Path) -> bytes:
    """A committed golden file, decompressed when ``.gz`` (gzip bytes depend on zlib)."""
    raw = path.read_bytes()
    return gzip.decompress(raw) if path.suffix == ".gz" else raw


def write_golden(path: Path, data: bytes) -> None:
    """Write a golden file; ``.gz`` with mtime 0 so regeneration is reproducible."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        data = gzip.compress(data, compresslevel=9, mtime=0)
    path.write_bytes(data)


def updating() -> bool:
    return os.environ.get("GIPPA_UPDATE_GOLDEN") == "1"


def first_diff(want: bytes, got: bytes) -> list[tuple[int, str, str]]:
    a, b = want.decode().splitlines(), got.decode().splitlines()
    out = [(i + 1, x, y) for i, (x, y) in enumerate(zip(a, b, strict=False)) if x != y][:5]
    if not out and len(a) != len(b):
        out = [(min(len(a), len(b)) + 1, f"{len(a)} lines", f"{len(b)} lines")]
    return out


# --------------------------------------------------------------------------- planted island


def planted_tus(work: Path, case: IslandCase = ISLAND) -> list[str]:
    """TU ids of the planted island: chr1 TUs in replicon order, ``first_tu`` .. + ``n_tus``."""
    tus = pd.read_parquet(work / "tus.parquet")
    chr1 = tus[tus["contig"].astype(str) == "chr1"].sort_values("start", kind="stable")
    ids = chr1["tu_id"].astype(str).tolist()
    assert len(ids) >= case.first_tu + case.n_tus, (len(ids), case)
    return ids[case.first_tu : case.first_tu + case.n_tus]


def planted_genes(work: Path, case: IslandCase = ISLAND) -> set[str]:
    tus = pd.read_parquet(work / "tus.parquet").set_index("tu_id")
    return {str(g) for t in planted_tus(work, case) for g in tus.loc[t, "genes"]}


def island_genes(row: Mapping[str, Any]) -> set[str]:
    cell = str(row["genes"])
    return set() if cell in ("", "NA") else set(cell.split(";"))


def hits(islands: pd.DataFrame, genes: set[str]) -> pd.DataFrame:
    """Island rows sharing at least one gene with ``genes``."""
    keep = [bool(island_genes(r) & genes) for _, r in islands.iterrows()]
    return islands[keep]


def tu_span(islands_row: Mapping[str, Any], work: Path) -> list[str]:
    """TU ids from ``tu_first`` to ``tu_last`` in chr1 order."""
    tus = pd.read_parquet(work / "tus.parquet")
    chr1 = tus[tus["contig"].astype(str) == str(islands_row["contig"])]
    ids = chr1.sort_values("start", kind="stable")["tu_id"].astype(str).tolist()
    a, b = ids.index(str(islands_row["tu_first"])), ids.index(str(islands_row["tu_last"]))
    return ids[a : b + 1] if a <= b else ids[a:] + ids[: b + 1]


def diagnostics(work: Path) -> pd.DataFrame:
    return pd.read_parquet(work / "context" / "island_diagnostics.parquet")


def l_prime_problems(work: Path, islands: pd.DataFrame, case: IslandCase = ISLAND) -> list[str]:
    """Decision 39: what is wrong with the island golden's block length (empty = fine).

    * chr1 was scanned with ``block_len: auto`` (the only rule the ADR-0021 calibration covers);
    * l' = min(20, 2 l-hat + 1) (the chromosome has far more than 20 scored TUs);
    * l' is shorter than the planted island, so a bootstrap block cannot copy the whole planted
      window into the null (ADR-0021: blocks >= 10 TUs copy the extreme window into ~63 % of
      replicates and kill power);
    * every chr1 island row reports that l' in ``block_len``.
    """
    d = diagnostics(work)
    row = d[d["contig"].astype(str) == "chr1"]
    if len(row) != 1:
        return [f"expected one chr1 diagnostics row, got {len(row)}"]
    r = row.iloc[0]
    out: list[str] = []
    if not bool(r["scanned"]) or str(r["block_rule"]) != "auto":
        out.append(f"chr1 scanned={r['scanned']} block_rule={r['block_rule']} (want auto)")
        return out
    l_hat, lp = int(r["l_hat"]), int(r["block_len"])
    if lp != min(L_MAX, 2 * l_hat + 1):
        out.append(f"l'={lp} != min({L_MAX}, 2*{l_hat}+1)")
    if lp >= case.n_tus:
        out.append(f"l'={lp} >= the planted island's {case.n_tus} TUs")
    for _, isl in islands[islands["contig"].astype(str) == "chr1"].iterrows():
        if str(isl["block_len"]) != str(lp):
            out.append(f"{isl['island_id']}: block_len {isl['block_len']} != l' {lp}")
    return out


def island_problems(work: Path, islands: pd.DataFrame, case: IslandCase = ISLAND) -> list[str]:
    """Spec example 2: exactly one island, ``high``, overlapping the planted TUs."""
    if len(islands) != 1:
        return [f"{len(islands)} islands reported, want exactly 1: {islands.to_dict('records')}"]
    isl = islands.iloc[0]
    out: list[str] = []
    if str(isl["direction"]) != "high":
        out.append(f"direction {isl['direction']} (want high)")
    if not set(tu_span(isl, work)) & set(planted_tus(work, case)):
        out.append(f"island {isl['tu_first']}..{isl['tu_last']} misses the planted TUs")
    if not island_genes(isl) & planted_genes(work, case):
        out.append("island gene list misses every planted gene")
    return out
