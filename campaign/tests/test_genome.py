"""NCBI genome package -> topology-tagged FASTA (SPEC §6)."""

from __future__ import annotations

import gzip
import tempfile
from pathlib import Path

import pytest
from giap_campaign.genome import find_genome_files, gbff_topologies, tag_fasta
from hypothesis import given, settings
from hypothesis import strategies as st

GBFF = """\
LOCUS       NC_000913            4641652 bp    DNA     circular CON 09-MAR-2022
DEFINITION  Escherichia coli str. K-12 substr. MG1655, complete genome.
ACCESSION   NC_000913
VERSION     NC_000913.3
KEYWORDS    RefSeq.
SOURCE      Escherichia coli str. K-12 substr. MG1655
FEATURES             Location/Qualifiers
     source          1..4641652
                     /organism="Escherichia coli str. K-12 substr. MG1655"
ORIGIN
        1 agcttttcat tctgactgca acgggcaata tgtctctgtg tggattaaaa aaagagtgtc
//
LOCUS       NZ_CP009001              2048 bp    DNA     linear   BCT 01-JAN-2020
DEFINITION  Example linear replicon; the word circular in DEFINITION is not used.
ACCESSION   NZ_CP009001
VERSION     NZ_CP009001.1
ORIGIN
        1 acgtacgtac
//
LOCUS       pEX1                     1000 bp    DNA              BCT 01-JAN-2020
DEFINITION  Example plasmid without a topology word.
ACCESSION   pEX1
VERSION     pEX1.2
ORIGIN
        1 acgt
//
"""

GBFF_MAP = {
    "NC_000913.3": "circular",
    "NC_000913": "circular",
    "NZ_CP009001.1": "linear",
    "NZ_CP009001": "linear",
    "pEX1.2": "linear",
    "pEX1": "linear",
}

FASTA = """\
>NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTC
TGATAGCAGC
>NZ_CP009001.1 Example linear replicon
ACGTACGTAC
>contig_9 not in the GenBank file
ACGT
>NC_000913
AC
"""

TAGGED = """\
>NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome [topology=circular]
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTC
TGATAGCAGC
>NZ_CP009001.1 Example linear replicon
ACGTACGTAC
>contig_9 not in the GenBank file
ACGT
>NC_000913 [topology=circular]
AC
"""


# --- gbff_topologies ----------------------------------------------------------------------------


def test_gbff_topologies_version_and_locus_keys(tmp_path: Path) -> None:
    gbff = tmp_path / "genomic.gbff"
    gbff.write_text(GBFF, encoding="utf-8")
    assert gbff_topologies(gbff) == GBFF_MAP


def test_gbff_topologies_gzip(tmp_path: Path) -> None:
    gbff = tmp_path / "genomic.gbff.gz"
    with gzip.open(gbff, "wt", encoding="utf-8") as fh:
        fh.write(GBFF)
    assert gbff_topologies(gbff) == GBFF_MAP


def test_gbff_topologies_empty_file(tmp_path: Path) -> None:
    gbff = tmp_path / "genomic.gbff"
    gbff.write_text("", encoding="utf-8")
    assert gbff_topologies(gbff) == {}


def test_gbff_topologies_minimal_records(tmp_path: Path) -> None:
    text = (
        "LOCUS       NC_000001               100 bp    DNA     circular BCT 01-JAN-2020\n"
        "VERSION     NC_000001.1\n//\n"
        "LOCUS       NC_000002               100 bp    DNA     linear   BCT 01-JAN-2020\n"
        "VERSION     NC_000002.1\n//\n"
    )
    gbff = tmp_path / "x.gbff"
    gbff.write_text(text, encoding="utf-8")
    got = gbff_topologies(gbff)
    assert got["NC_000001.1"] == "circular"
    assert got["NC_000002.1"] == "linear"


# --- tag_fasta ----------------------------------------------------------------------------------


def test_tag_fasta_appends_tag_to_circular_only(tmp_path: Path) -> None:
    src = tmp_path / "genome_genomic.fna"
    src.write_text(FASTA, encoding="utf-8")
    dst = tmp_path / "genome.fna"
    counts = tag_fasta(src, dst, GBFF_MAP)
    assert counts == {"circular": 2, "linear": 1, "unknown": 1}
    assert dst.read_text(encoding="utf-8") == TAGGED


def test_tag_fasta_gzip_input_plain_output(tmp_path: Path) -> None:
    src = tmp_path / "genome_genomic.fna.gz"
    with gzip.open(src, "wt", encoding="utf-8") as fh:
        fh.write(FASTA)
    dst = tmp_path / "genome.fna"
    counts = tag_fasta(src, dst, GBFF_MAP)
    assert counts == {"circular": 2, "linear": 1, "unknown": 1}
    raw = dst.read_bytes()
    assert raw[:2] != b"\x1f\x8b"
    assert raw.decode("utf-8") == TAGGED


def test_tag_fasta_empty_map_everything_unknown(tmp_path: Path) -> None:
    src = tmp_path / "in.fna"
    src.write_text(FASTA, encoding="utf-8")
    dst = tmp_path / "out.fna"
    counts = tag_fasta(src, dst, {})
    assert counts == {"circular": 0, "linear": 0, "unknown": 4}
    assert dst.read_text(encoding="utf-8") == FASTA


def test_tag_fasta_end_to_end_with_gbff(tmp_path: Path) -> None:
    gbff = tmp_path / "genomic.gbff"
    gbff.write_text(GBFF, encoding="utf-8")
    src = tmp_path / "in.fna"
    src.write_text(FASTA, encoding="utf-8")
    dst = tmp_path / "out.fna"
    tag_fasta(src, dst, gbff_topologies(gbff))
    headers = [ln for ln in dst.read_text(encoding="utf-8").splitlines() if ln.startswith(">")]
    assert [h.endswith(" [topology=circular]") for h in headers] == [True, False, False, True]


_ID = st.text(alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.", min_size=1, max_size=12)


@settings(max_examples=60, deadline=None)
@given(
    records=st.lists(
        st.tuples(
            _ID,
            st.sampled_from(["", " some description"]),
            st.sampled_from(["circular", "linear", None]),
            st.text(alphabet="ACGTN", min_size=1, max_size=30),
        ),
        min_size=0,
        max_size=6,
        unique_by=lambda r: r[0],
    )
)
def test_tag_fasta_property_strip_tag_restores_input(
    records: list[tuple[str, str, str | None, str]],
) -> None:
    text = "".join(f">{rid}{desc}\n{seq}\n" for rid, desc, _, seq in records)
    topo = {rid: t for rid, _, t, _ in records if t is not None}
    with tempfile.TemporaryDirectory() as d:
        src, dst = Path(d) / "in.fna", Path(d) / "out.fna"
        src.write_text(text, encoding="utf-8")
        counts = tag_fasta(src, dst, topo)
        out = dst.read_text(encoding="utf-8")
    assert out.replace(" [topology=circular]\n", "\n") == text
    assert out.count(" [topology=circular]") == counts["circular"]
    assert counts == {
        "circular": sum(t == "circular" for _, _, t, _ in records),
        "linear": sum(t == "linear" for _, _, t, _ in records),
        "unknown": sum(t is None for _, _, t, _ in records),
    }


# --- find_genome_files --------------------------------------------------------------------------

ACC = "GCF_000005845.2"


def _touch(root: Path, names: list[str]) -> list[Path]:
    paths = []
    for name in names:
        p = root / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("x", encoding="utf-8")
        paths.append(p)
    return paths


def test_find_genome_files_ncbi_layout(tmp_path: Path) -> None:
    paths = _touch(
        tmp_path,
        [
            "README.md",
            "md5sum.txt",
            "ncbi_dataset/data/assembly_data_report.jsonl",
            "ncbi_dataset/data/dataset_catalog.json",
            f"ncbi_dataset/data/{ACC}/genomic.gbff",
            f"ncbi_dataset/data/{ACC}/{ACC}_ASM584v2_genomic.fna",
        ],
    )
    fna, gbff = find_genome_files(paths)
    assert fna == tmp_path / f"ncbi_dataset/data/{ACC}/{ACC}_ASM584v2_genomic.fna"
    assert gbff == tmp_path / f"ncbi_dataset/data/{ACC}/genomic.gbff"


def test_find_genome_files_ignores_directories(tmp_path: Path) -> None:
    paths = _touch(tmp_path, ["d/a_genomic.fna", "d/genomic.gbff"])
    fna, gbff = find_genome_files([tmp_path / "d", *paths])
    assert fna.name == "a_genomic.fna"
    assert gbff.name == "genomic.gbff"


@pytest.mark.parametrize(
    "names",
    [
        ["d/genomic.gbff"],
        ["d/a_genomic.fna"],
        ["d/a_genomic.fna", "d/b_genomic.fna", "d/genomic.gbff"],
        ["d/a_genomic.fna", "d/genomic.gbff", "e/genomic.gbff"],
        [],
    ],
    ids=["no-fna", "no-gbff", "two-fna", "two-gbff", "empty"],
)
def test_find_genome_files_errors(tmp_path: Path, names: list[str]) -> None:
    with pytest.raises(ValueError):
        find_genome_files(_touch(tmp_path, names))
