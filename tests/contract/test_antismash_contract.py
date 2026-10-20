"""Contract test: antiSMASH 8.0.4 JSON parser on recorded runs (M2-18; GIAP_SPEC.md §7.1.1).

Committed recordings (``tests/data/recorded/antismash/PROVENANCE.md``): trimmed ``bakta.json``
(areas complete, CDS overlapping areas plus flanking CDS) of sim100k, K-12 MG1655 and Nissle 1917,
plus ``replicons.tsv``. Input was Bakta's GBFF, so core genes carry Bakta locus tags.

* sim100k: no areas → ``[]``; K-12: 3 regions; Nissle: 6 regions (acceptance: ≥ 1 BGC).
* core genes are the ``gene_kind == biosynthetic`` CDS only (``biosynthetic-additional`` excluded),
  including the fuzzy ``[3082708:>3082864](-)`` CDS area in Nissle region 3.
* full recordings under ``$GIAP_ROOT/recordings/antismash/`` (skipped when missing): the full JSON
  parses to exactly the result of the trimmed copy, and every core gene's locus tag is in the
  Bakta GBFF that antiSMASH read.
"""

from __future__ import annotations

import csv
import os
import re
from pathlib import Path

import pytest

from giap.model.coords import Interval, Replicon, check_interval

RECORDED = Path(__file__).resolve().parents[1] / "data" / "recorded" / "antismash" / "8.0.4"
SAMPLES = ["sim100k", "k12_mg1655", "nissle1917"]
FULL = Path(os.environ.get("GIAP_ROOT", "/nonexistent")) / "recordings" / "antismash"

K12 = "U00096.3"
K12_REGIONS = [
    (430201, 451101, ("terpene-precursor",), ("KOIAIB_00460",)),
    (
        594156,
        649297,
        ("NRP-metallophore", "NRPS"),
        ("KOIAIB_00640", "KOIAIB_00647", "KOIAIB_00650"),
    ),
    (940339, 966632, ("azole-containing-RiPP",), ("KOIAIB_00973", "KOIAIB_00976")),
]
NISSLE = "CP171242.1"
NISSLE_REGIONS = [
    (880587, 913002, ("NI-siderophore",), ("JMJJOG_00885", "JMJJOG_00887")),
    (
        1944994,
        2078575,
        ("NRPS", "T1PKS", "NRP-metallophore"),
        (
            "JMJJOG_01918",
            "JMJJOG_01924",
            "JMJJOG_01925",
            "JMJJOG_01926",
            "JMJJOG_01927",
            "JMJJOG_01930",
            "JMJJOG_01953",
            "JMJJOG_01954",
            "JMJJOG_01959",
        ),
    ),
    (3081475, 3125063, ("arylpolyene",), ("JMJJOG_03054", "JMJJOG_03057")),
    (3231195, 3257488, ("azole-containing-RiPP",), ("JMJJOG_03184", "JMJJOG_03187")),
    (
        3523284,
        3578270,
        ("NRP-metallophore", "NRPS"),
        ("JMJJOG_03486", "JMJJOG_03489", "JMJJOG_03496"),
    ),
    (3688911, 3709811, ("terpene-precursor",), ("JMJJOG_03629",)),
]
BIOSYNTHETIC_ADDITIONAL = {
    "KOIAIB_00459",
    "KOIAIB_00628",
    "KOIAIB_00648",
    "KOIAIB_00649",
    "KOIAIB_00654",
    "KOIAIB_00655",
    "KOIAIB_00661",
    "KOIAIB_00663",
    "KOIAIB_00978",
}


def _replicons(sample: str) -> dict[str, Replicon]:
    with (RECORDED / sample / "replicons.tsv").open(newline="") as fh:
        return {
            r["contig"]: Replicon(r["contig"], int(r["length"]), r["topology"])  # type: ignore[arg-type]
            for r in csv.DictReader(fh, delimiter="\t")
        }


def _parse(sample: str, path: Path | None = None):  # type: ignore[no-untyped-def]
    from giap.annotate.antismash import parse_antismash_json

    return parse_antismash_json(path or RECORDED / sample / "bakta.json", _replicons(sample))


def test_sim100k_has_no_regions() -> None:
    assert _parse("sim100k") == []


@pytest.mark.parametrize(
    ("sample", "contig", "want"),
    [("k12_mg1655", K12, K12_REGIONS), ("nissle1917", NISSLE, NISSLE_REGIONS)],
)
def test_regions(
    sample: str, contig: str, want: list[tuple[int, int, tuple[str, ...], tuple[str, ...]]]
) -> None:
    from giap.annotate.antismash import BgcRegion

    got = _parse(sample)
    assert len(got) >= 1
    reps = _replicons(sample)
    for r in got:
        check_interval(r.interval, reps[r.contig])
    assert sorted(got, key=lambda r: r.interval.start) == [
        BgcRegion(contig, Interval(contig, s, e, "."), prods, core) for s, e, prods, core in want
    ]


def test_k12_biosynthetic_additional_genes_are_not_core() -> None:
    core = {g for r in _parse("k12_mg1655") for g in r.core_genes}
    assert core.isdisjoint(BIOSYNTHETIC_ADDITIONAL)


def test_nissle_fuzzy_cds_record_parses() -> None:
    """Nissle area 3 contains a CDS at ``[3082708:>3082864](-)``; parsing must not choke on it."""
    text = (RECORDED / "nissle1917" / "bakta.json").read_text()
    assert "[3082708:>3082864](-)" in text
    (r,) = [r for r in _parse("nissle1917") if r.interval.start == 3081475]
    assert r.core_genes == ("JMJJOG_03054", "JMJJOG_03057")


@pytest.mark.parametrize("sample", SAMPLES)
def test_parse_is_deterministic(sample: str) -> None:
    assert _parse(sample) == _parse(sample)


# ---------------------------------------------------------------- full recordings (not committed)


@pytest.mark.parametrize("sample", SAMPLES)
def test_full_recording_parses_like_the_excerpt(sample: str) -> None:
    full = FULL / sample / "antismash" / "bakta.json"
    if not full.is_file():
        pytest.skip(f"full antiSMASH recording missing: {full}")
    got = _parse(sample, full)
    assert got == _parse(sample)
    gbff = FULL / sample / "bakta.gbff"
    if gbff.exists():
        tags = set(re.findall(r'/locus_tag="([^"]+)"', gbff.read_text()))
        assert {g for r in got for g in r.core_genes} <= tags
