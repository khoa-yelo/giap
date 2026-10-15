"""The GENOMAD / ANTISMASH Nextflow stubs produce what the M2-18 parsers read (M2-07, M2-18).

The stub ``printf`` headers of ``workflow/modules/local/genomad/main.nf`` are written to
``genomad/genome_summary/genome_*.tsv`` (prefix ``genome``: the process stages ``genome.fa``) and
must parse to an empty ``GenomadResult``; the ANTISMASH stub JSON ``{"records": []}`` → ``[]``.
The headers must also equal the recorded geNOMAD 1.12.0 headers. No Nextflow needed.
"""

from __future__ import annotations

import re
from pathlib import Path

from giap.model.coords import Replicon

REPO = Path(__file__).resolve().parents[2]
MODULES = REPO / "workflow" / "modules" / "local"
RECORDED = REPO / "tests" / "data" / "recorded" / "genomad" / "1.12.0" / "k12_mg1655"
REPS = {"chr1": Replicon("chr1", 1000, "circular")}


def _stub(tool: str) -> str:
    text = (MODULES / tool / "main.nf").read_text()
    return text[text.index("stub:") :]


def _genomad_stub_headers() -> dict[str, str]:
    out = {}
    for fmt, name in re.findall(r"printf '([^']*)' > (genome_\w+\.tsv)", _stub("genomad")):
        out[name] = fmt.replace("\\\\t", "\t").replace("\\\\n", "\n")
    return out


def test_genomad_stub_writes_the_three_parsed_tsvs_with_recorded_headers() -> None:
    headers = _genomad_stub_headers()
    names = {f"genome_{n}.tsv" for n in ("virus_summary", "plasmid_summary", "virus_genes")}
    assert set(headers) == names
    assert "genomad/genome_summary" in _stub("genomad")
    for name, text in headers.items():
        recorded = (RECORDED / "genome_summary" / name).read_text().splitlines()[0]
        assert text.rstrip("\n") == recorded, name


def test_genomad_stub_output_parses_empty(tmp_path: Path) -> None:
    from giap.annotate.genomad import GenomadResult, parse_genomad

    d = tmp_path / "genomad" / "genome_summary"
    d.mkdir(parents=True)
    for name, text in _genomad_stub_headers().items():
        (d / name).write_text(text)
    assert parse_genomad(tmp_path / "genomad", "genome", REPS) == GenomadResult((), (), ())


def test_antismash_stub_json_parses_empty(tmp_path: Path) -> None:
    from giap.annotate.antismash import parse_antismash_json

    m = re.search(r"echo '(\{.*\})' > antismash/\$\{meta_id\}\.json", _stub("antismash"))
    assert m, "ANTISMASH stub JSON not found"
    p = tmp_path / "sample.json"
    p.write_text(m.group(1))
    assert parse_antismash_json(p, REPS) == []
