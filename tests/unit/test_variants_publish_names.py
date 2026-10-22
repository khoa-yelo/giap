"""The local executor publishes the three Module 5 files under the Nextflow names
(docs/specs/M7-04-variants-step-and-run.md rule 3, lead decision 3; GIAP_SPEC.md §6.4, §9.1).

``giap.runner.local.PUBLISHED_FILES`` maps a work key to its ``OUT`` path
(tests/unit/test_local_publish_names.py on main pins it against every ``'key=path'`` of
workflow/). Here the three variants keys must be in it. Their names must equal those in
variants.nf, and ``publish`` must copy ``work/<key>`` byte for byte. ``work/pileup`` is already
gzip (the step writes gzip whatever its path is called), so ``tables/pileup.tsv.gz`` is a valid
gzip file.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NAMES = {
    "pileup": "tables/pileup.tsv.gz",
    "variants": "tables/variants.tsv",
    "variants_vcf": "tables/variants.vcf",
}


def test_published_files_name_the_three_variants_tables() -> None:
    from giap.runner.local import PUBLISHED_FILES

    assert {k: PUBLISHED_FILES.get(k) for k in NAMES} == NAMES


def test_variants_nf_uses_the_same_names() -> None:
    nf = ROOT / "workflow" / "subworkflows" / "local" / "variants.nf"
    assert nf.is_file(), "workflow/subworkflows/local/variants.nf missing"
    text = re.sub(r"(?m)//.*$", "", nf.read_text(encoding="utf-8"))
    found = dict(re.findall(r"'([a-z_0-9]+)=([a-z_]+/[^'\s]+)", text))
    assert {k: found.get(k) for k in NAMES} == NAMES


def test_publish_copies_the_three_files(tmp_path: Path) -> None:
    from giap.runner.local import publish

    work, out = tmp_path / "work", tmp_path / "OUT"
    work.mkdir()
    pile = gzip.compress(b"contig\tcodon_start\n", mtime=0)
    (work / "pileup").write_bytes(pile)
    (work / "variants").write_text("variant_id\n", encoding="utf-8")
    (work / "variants_vcf").write_text("##fileformat=VCFv4.3\n", encoding="utf-8")
    done = publish(work, out)
    assert (out / "tables" / "pileup.tsv.gz").read_bytes() == pile
    assert (out / "tables" / "variants.tsv").read_text(encoding="utf-8") == "variant_id\n"
    assert (out / "tables" / "variants.vcf").read_text(encoding="utf-8") == "##fileformat=VCFv4.3\n"
    for rel in NAMES.values():
        assert out / rel in done, rel


def test_publish_skips_absent_variants_files(tmp_path: Path) -> None:
    """No peptides: the step never ran, nothing is published under those names."""
    from giap.runner.local import publish

    work, out = tmp_path / "work", tmp_path / "OUT"
    work.mkdir()
    publish(work, out)
    for rel in NAMES.values():
        assert not (out / rel).exists(), rel
