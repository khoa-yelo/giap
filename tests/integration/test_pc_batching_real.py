"""M2-14b example 4: batched Promoter Calculator output is byte-identical to the unbatched
pipeline (docs/specs/M2-14b-promoter-calculator-batching.md item 3; GIPPA_SPEC.md §7.1.1).

Real tool, ``gippa-promoter_calculator`` prebuilt env (workflow/conf/prebuilt_envs.config), the
module run through a driver script written into ``tmp_path`` (the workflow/tests/real pattern:
``-c workflow/tests/real/real.config -profile farmshare_local`` publishes the declared outputs into
``<launchDir>/out``). The driver groups the regions as annotate.nf does (sorted by TU id,
``collate(100)``; annotate.nf itself is pinned by tests/unit/test_nf_pc_batching_static.py).

* sim100k: the TU upstream regions of the merged ``parse_bakta`` / ``uniop_input`` /
  ``operon_pairs`` steps on gippasim SIM100K (recorded ``uniop.pred``, as tests/helpers/
  layers_chain.py). sim100k has 20 TUs (the spec's "288 TUs -> 3 tasks" does not hold; test
  author), so this is one batch. "Unbatched" is the M2-14 per-TU task script run once per region
  (``promoter-calculator -i <TU>.fa -o <TU>.csv -j 1``, then the CRLF header when the tool wrote
  no CSV). Every ``<TU>.csv`` and the ``tx_layer`` step's ``tx_rates.parquet`` must be equal byte
  for byte.
* three batches on recorded inputs: the first 250 K-12 MG1655 regions (TU_00001 ... TU_00250) of
  ``$GIPPA_ROOT/recordings/promoter_calculator/k12_mg1655`` (recorded one region per call with the
  same command, PROVENANCE.md) -> tasks of 100, 100 and 50, every CSV equal to the recorded one;
  the two regions without a recorded CSV (TU_00058, TU_00059) get the CRLF header.

``slow`` (minutes of real tool time per 100 regions): run on Slurm (test_branch.sbatch / nightly).
"""

from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"
ROOT = Path(os.environ.get("GIPPA_ROOT", REPO / "local"))
ENV = ROOT / "micromamba" / "envs" / "gippa-promoter_calculator"
K12 = ROOT / "recordings" / "promoter_calculator" / "k12_mg1655"
HEADER_LINE = REPO / "tests" / "data" / "recorded" / "headers" / "promoter_calculator.csv.txt"
HEADER = HEADER_LINE.read_text(encoding="utf-8").splitlines()[0]
PC_BATCH = 100

pytestmark = [
    pytest.mark.slow,
    pytest.mark.needs_nextflow,
    pytest.mark.needs_conda,
    pytest.mark.skipif(not ENV.is_dir(), reason=f"missing capability: tool env {ENV}"),
]

DRIVER = """\
nextflow.enable.dsl = 2
include {{ PROMOTER_CALCULATOR }} from '{module}'

workflow {{
    PROMOTER_CALCULATOR(channel.fromPath('{regions}/*.fa')
        .toSortedList {{ a, b -> a.baseName <=> b.baseName }}
        .flatMap {{ l -> l.collate({batch}) }}
        .map {{ b -> ["pc_${{b[0].baseName}}_${{b[-1].baseName}}", b] }})
}}
"""


def _batched(tmp: Path, regions: Path) -> tuple[dict[str, str], Path]:
    """Run the batched module on ``regions``; returns ({task tag: status}, published out/)."""
    run = tmp / "nf"
    run.mkdir()
    script = run / "pc_batched.nf"
    module = WF / "modules" / "local" / "promoter_calculator" / "main.nf"
    script.write_text(DRIVER.format(module=module, regions=regions, batch=PC_BATCH))
    trace = run / "trace.tsv"
    cmd = ["nextflow", "-log", str(run / "nextflow.log"), "run", str(script)]
    cmd += ["-c", str(WF / "tests" / "real" / "real.config"), "-profile", "farmshare_local"]
    cmd += ["-work-dir", str(run / "work"), "-with-trace", str(trace)]
    env = {**os.environ, "NXF_ANSI_LOG": "false"}
    res = subprocess.run(cmd, cwd=run, capture_output=True, text=True, env=env, timeout=4 * 3600)
    assert res.returncode == 0, f"{' '.join(cmd)}\n{res.stdout[-4000:]}\n{res.stderr[-4000:]}"
    with trace.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    tags = {}
    for r in rows:
        m = re.fullmatch(r"(?P<proc>[^(]+?)\s*\((?P<tag>.*)\)\s*", r["name"])
        assert m and m["proc"].endswith("PROMOTER_CALCULATOR"), r["name"]
        tags[m["tag"]] = r["status"]
    return tags, run / "out"


def _unbatched_one(region: Path, out: Path) -> None:
    """The M2-14 per-TU task script for one region (the unbatched pipeline)."""
    csv_path = out / f"{region.stem}.csv"
    cmd = ["micromamba", "run", "-p", str(ENV), "promoter-calculator"]
    cmd += ["-i", str(region), "-o", str(csv_path), "-j", "1"]
    res = subprocess.run(cmd, cwd=out, capture_output=True, text=True, timeout=3600)
    assert res.returncode == 0, f"{region.name}: {res.stdout[-2000:]}\n{res.stderr[-2000:]}"
    if not csv_path.exists():
        csv_path.write_bytes(f"{HEADER}\r\n".encode())


def _unbatched(regions: list[Path], out: Path) -> None:
    out.mkdir()
    workers = max(1, min(8, int(os.environ.get("SLURM_CPUS_PER_TASK", "4")) - 1))
    with ThreadPoolExecutor(workers) as pool:
        list(pool.map(lambda r: _unbatched_one(r, out), regions))


def _csvs(d: Path) -> dict[str, bytes]:
    return {p.name: p.read_bytes() for p in sorted(d.glob("*.csv"))}


def _sim100k_regions(root: Path) -> tuple[Path, Path, Path]:
    """(tu_upstream dir, tus.parquet, replicons.parquet) of the merged steps on SIM100K."""
    from helpers import layers_chain as lc

    import gippasim
    from gippa.config import GippaConfig

    sim = gippasim.simulate(gippasim.SIM100K, gippasim.SIM100K_SEED, root / "sim")
    w = root / "work"
    w.mkdir()
    cfg = GippaConfig(seed=1)
    with lc.registry() as reg:
        lc._run(reg, "parse_bakta", cfg, {"bakta_json": sim.bakta_json, "genome": sim.genome_fasta},
                {k: w / v for k, v in lc.PARSE_OUT.items()})  # fmt: skip
        core = {k: w / lc.PARSE_OUT[k] for k in ("genes", "replicons")}
        lc._run(reg, "uniop_input", cfg, core,
                {"uniop_cds": w / "uniop_cds.gff", "uniop_map": w / "uniop_map.tsv"})  # fmt: skip
        shutil.copyfile(lc.UNIOP_REC / "uniop.pred", w / "uniop.pred")
        lc._run(reg, "operon_pairs", cfg,
                {**core, "genome": sim.genome_fasta, "uniop_cds": w / "uniop_cds.gff",
                 "uniop_map": w / "uniop_map.tsv", "uniop_pred": w / "uniop.pred"},
                {"pairs_raw": w / "pairs_raw.parquet", "tus": w / "tus.parquet",
                 "tu_upstream": w / "tu_upstream",
                 "operon_meta": w / "operon_meta.json"})  # fmt: skip
    return w / "tu_upstream", w / "tus.parquet", core["replicons"]


def _tx_rates(csvs: list[Path], tus: Path, replicons: Path, out: Path) -> bytes:
    from helpers import layers_chain as lc

    from gippa.config import GippaConfig

    with lc.registry() as reg:
        lc._run(reg, "tx_layer", GippaConfig(seed=1),
                {"pc": csvs, "tus": tus, "replicons": replicons},
                {"tx_rates": out})  # fmt: skip
    return out.read_bytes()


def test_example_4_sim100k_batched_equals_unbatched_csvs_and_tx_rates(tmp_path: Path) -> None:
    regions, tus, replicons = _sim100k_regions(tmp_path / "chain")
    fas = sorted(regions.glob("*.fa"))
    assert fas, "operon_pairs wrote no TU upstream region"
    n_batches = -(-len(fas) // PC_BATCH)

    tags, out = _batched(tmp_path, regions)
    assert len(tags) == n_batches, f"{len(fas)} regions -> {n_batches} task(s): {sorted(tags)}"
    batched = _csvs(out)
    assert sorted(batched) == sorted(f"{f.stem}.csv" for f in fas), sorted(batched)

    ref_dir = tmp_path / "unbatched"
    _unbatched(fas, ref_dir)
    ref = _csvs(ref_dir)
    diff = sorted(n for n in ref if batched.get(n) != ref[n])
    assert not diff, f"per-TU CSVs differ from the unbatched run: {diff}"

    a = _tx_rates(sorted(out.glob("*.csv")), tus, replicons, tmp_path / "tx_batched.parquet")
    b = _tx_rates(sorted(ref_dir.glob("*.csv")), tus, replicons, tmp_path / "tx_unbatched.parquet")
    assert a == b, "tx_rates.parquet differs between batched and unbatched"


@pytest.mark.skipif(not (K12 / "regions").is_dir(), reason=f"no K-12 recording at {K12}")
def test_example_4_three_batches_equal_the_recorded_one_region_per_call_csvs(
    tmp_path: Path,
) -> None:
    regions = tmp_path / "regions"
    regions.mkdir()
    ids = [f"TU_{i:05d}" for i in range(1, 251)]
    for tu in ids:
        (regions / f"{tu}.fa").symlink_to(K12 / "regions" / f"{tu}.fa")

    tags, out = _batched(tmp_path, regions)
    assert tags == {
        "pc_TU_00001_TU_00100": "COMPLETED",
        "pc_TU_00101_TU_00200": "COMPLETED",
        "pc_TU_00201_TU_00250": "COMPLETED",
    }, tags
    batched = _csvs(out)
    assert sorted(batched) == [f"{tu}.csv" for tu in ids], f"{len(batched)} CSVs"
    no_csv = sorted(tu for tu in ids if not (K12 / "csv" / f"{tu}.csv").is_file())
    assert no_csv == ["TU_00058", "TU_00059"], no_csv
    header_only = f"{HEADER}\r\n".encode()
    diff = []
    for tu in ids:
        rec = K12 / "csv" / f"{tu}.csv"
        want = rec.read_bytes() if rec.is_file() else header_only
        if batched[f"{tu}.csv"] != want:
            diff.append(tu)
    assert not diff, f"batched CSVs differ from the recorded one-region-per-call run: {diff}"
