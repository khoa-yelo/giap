"""antiSMASH review tests (M2-18 lead decisions 4 and 5 after the PR #84 code review).

Each test kills a mutant that survived the M2-18 suite:

* ``databases`` is a required argument of ``antismash_command`` (decision 4; the builder used to
  drop ``--databases`` when it was omitted, unlike ``main.nf``);
* a biosynthetic CDS crossing the area's right edge is not a core gene ("wholly inside");
* a ``join{}`` CDS with one part outside the area is excluded (every part must lie inside);
* an area with ``end > length`` on a linear record is corrupt input → ``ValueError``;
* ``start > end`` on a linear record is rejected by antiSMASH's own check (message names the tool),
  not merely by the generic interval check after an unwrap.
"""

from __future__ import annotations

import inspect
import json
from pathlib import Path
from typing import Any

import pytest

from gippa.model.coords import Interval, Replicon

REPS = {
    "chr1": Replicon("chr1", 100_000, "circular", "chromosome"),
    "lin": Replicon("lin", 50_000, "linear", "unknown"),
}


def _cds(loc: str, tag: str) -> dict[str, Any]:
    return {
        "location": loc,
        "type": "CDS",
        "qualifiers": {"locus_tag": [tag], "gene_kind": ["biosynthetic"]},
    }


def _write(tmp: Path, rid: str, start: int, end: int, feats: list[dict[str, Any]]) -> Path:
    topo = REPS[rid].topology
    rec = {
        "id": rid,
        "name": rid,
        "annotations": {"topology": topo},
        "features": feats,
        "areas": [{"start": start, "end": end, "products": ["NRPS"]}],
    }
    p = tmp / "genome.json"
    p.write_text(json.dumps({"version": "8.0.4", "records": [rec], "schema": 4}))
    return p


# ---------------------------------------------------------------- decision 4: --databases


def test_databases_is_required() -> None:
    from gippa.annotate.antismash_cmd import antismash_command

    param = inspect.signature(antismash_command).parameters["databases"]
    assert param.default is inspect.Parameter.empty
    with pytest.raises(TypeError):
        antismash_command(Path("bakta.gbff"), Path("out"), 8)  # type: ignore[call-arg]


def test_databases_passed_by_keyword() -> None:
    from gippa.annotate.antismash_cmd import antismash_command

    got = antismash_command(Path("g.gbff"), Path("o"), 2, databases=Path("db_dir/antismash"))
    i = got.index("--databases")
    assert got[i + 1] == "db_dir/antismash"
    assert got.index("--minimal") < i < got.index("--cpus")


# ---------------------------------------------------------------- decision 5: surviving mutants


def test_cds_crossing_right_edge_is_not_core(tmp_path: Path) -> None:
    from gippa.annotate.antismash import parse_antismash_json

    feats = [_cds("[2000:3000](+)", "IN_1"), _cds("[9500:10500](+)", "EDGE_R")]
    (r,) = parse_antismash_json(_write(tmp_path, "lin", 1000, 10_000, feats), REPS)
    assert r.core_genes == ("IN_1",)


def test_cds_crossing_left_edge_is_not_core(tmp_path: Path) -> None:
    from gippa.annotate.antismash import parse_antismash_json

    feats = [_cds("[500:1500](-)", "EDGE_L"), _cds("[2000:3000](+)", "IN_1")]
    (r,) = parse_antismash_json(_write(tmp_path, "lin", 1000, 10_000, feats), REPS)
    assert r.core_genes == ("IN_1",)


def test_join_cds_with_one_part_outside_is_excluded(tmp_path: Path) -> None:
    from gippa.annotate.antismash import parse_antismash_json

    feats = [
        _cds("join{[2000:2600](+), [2700:3000](+)}", "JOIN_IN"),
        _cds("join{[8000:9000](+), [12000:12500](+)}", "JOIN_OUT"),
    ]
    (r,) = parse_antismash_json(_write(tmp_path, "lin", 1000, 10_000, feats), REPS)
    assert r.core_genes == ("JOIN_IN",)


def test_join_across_origin_with_part_outside_wrapped_area_is_excluded(tmp_path: Path) -> None:
    from gippa.annotate.antismash import parse_antismash_json

    feats = [
        _cds("join{[99800:100000](+), [0:400](+)}", "WRAP_IN"),
        _cds("join{[99800:100000](+), [0:6000](+)}", "WRAP_OUT"),
    ]
    (r,) = parse_antismash_json(_write(tmp_path, "chr1", 95_000, 4_000, feats), REPS)
    assert r.interval == Interval("chr1", 95_000, 104_000, ".")
    assert r.core_genes == ("WRAP_IN",)


@pytest.mark.parametrize(("start", "end"), [(45_000, 50_001), (49_000, 60_000)])
def test_area_end_beyond_linear_contig_raises(tmp_path: Path, start: int, end: int) -> None:
    from gippa.annotate.antismash import parse_antismash_json

    with pytest.raises(ValueError):
        parse_antismash_json(_write(tmp_path, "lin", start, end, []), REPS)


def test_wrapping_area_on_linear_record_raises_its_own_error(tmp_path: Path) -> None:
    from gippa.annotate.antismash import parse_antismash_json

    with pytest.raises(ValueError, match=r"antiSMASH.*lin"):
        parse_antismash_json(_write(tmp_path, "lin", 40_000, 1_000, []), REPS)


def test_area_end_at_linear_contig_length_is_valid(tmp_path: Path) -> None:
    from gippa.annotate.antismash import parse_antismash_json

    (r,) = parse_antismash_json(_write(tmp_path, "lin", 45_000, 50_000, []), REPS)
    assert r.interval == Interval("lin", 45_000, 50_000, ".")
