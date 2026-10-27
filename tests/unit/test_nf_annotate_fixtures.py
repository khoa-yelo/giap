"""Inputs of ``workflow/tests/annotate.nf.test`` stay valid (M2-07; M2-02 precomputed formats).

``annotation.json`` must equal the default ``annotation`` section that ``gippa run`` writes into
``params.json`` (runner.write_params), so the nf-test exercises the real parameter shape; the
precomputed TSVs must load with ``load_precomputed``.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from gippa.annotate.precomputed import load_precomputed
from gippa.config import GippaConfig
from gippa.model.coords import Replicon

DATA = Path(__file__).resolve().parents[2] / "workflow" / "tests" / "data" / "annotate"


def test_annotation_json_equals_the_default_annotation_section() -> None:
    want = GippaConfig().model_dump(mode="json")["annotation"]
    got = json.loads((DATA / "annotation.json").read_text(encoding="utf-8"))
    assert got == want, "regenerate workflow/tests/data/annotate/annotation.json from the defaults"


def test_genome_fixture_is_one_circular_contig() -> None:
    lines = (DATA / "genome.fa").read_text().splitlines()
    assert lines[0] == ">chr1"
    assert sum(len(x) for x in lines[1:]) == 3000
    rows = [r.split("\t") for r in (DATA / "replicons.tsv").read_text().splitlines()]
    # M2-07 review decision 2: the preflight format (preflight/run.py), topology `c` / `l`
    assert rows == [["chr1", "chr1", "chromosome", "c", "-"]]


@pytest.mark.parametrize("key", ["terminators", "regions", "operon_pairs"])
def test_precomputed_fixtures_load(key: str) -> None:
    reps = {"chr1": Replicon("chr1", 3000, "circular", "chromosome")}
    frame = load_precomputed(key, DATA / f"{key}.tsv", reps)
    assert len(frame) == 1
