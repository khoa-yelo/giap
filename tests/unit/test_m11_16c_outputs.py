"""M11-16c review items M1, L1 and the origin-spanning TGA (#474; docs/specs/
M11-16-translation-tables-4-25.md item 8: every gippasim module that reads the genetic code
takes it from the config).

M1: one end-to-end test per table (4, 25) reads every simulator output that carries residues or
stops and checks it against the genome translated under the table, so a site still hard-wired
to table 11 (proteomics, features, scenario_truth, strain, blocks, readthrough) shows up:

- the sample and reference Bakta JSON: ``nt`` is the genome slice, ends in TAA/TAG, has no
  internal TAA/TAG, ``aa`` is its translation; no ``*`` in any ``aa`` (pseudogenes included);
- ``proteins.faa`` equals the annotated proteins and has no ``*``;
- every correct peptide of ``truth/peptides.tsv`` occurs in the six-frame translation of the
  genome under the table (``n_loci_genome >= 1``); no peptide contains ``*``;
- ``truth/novel_orfs.tsv``: ``aa`` is the translation of the genome slice, which ends in TAA/TAG;
- ``truth/alt_sites.tsv`` and ``truth/stop_loss.tsv`` residues, and the SNV ``aa_change`` of
  ``truth/strain_differences.tsv``, are read under the table;
- read-through extensions (``gippasim.readthrough.ReadThrough.extension``, captured from
  ``simulate``) are the translation of the codons after the stop, with no ``*``;
- tRNA ``amino_acid`` is the table's residue of the anticodon's codon (TCA -> Trp at 4, Gly at 25);
- the v2 blocks (oriC ``mnmG``, prophage, BGC genes) follow the table, and prophage genes, drawn
  from the tilted sense usage, contain TGA;
- six-frame wrong peptides (``wrong_source = "sixframe_noncoding"``) occur in the genome.

Each check also asserts that the run exercised TGA (a site at a TGA codon, a TGA anticodon, an
extension with TGA), so it cannot pass vacuously. A high ``tga_sense_share`` gives the power.

L1: the intra-TU overlap of 4 bp needs the upstream stop to start with ``TG`` (the stop is
``xx`` + ``A`` ``TG`` + the start's ``TG``), which only TGA can do; at 4/25 it is infeasible. The
spec says nothing about overlap rates, so the tests pin only that 4 bp is not a feasible
overlap at 4/25: ``SimConfig(genetic_code=4|25, overlap_bp=((4, 1.0),))`` raises ``ValueError``
(as an infeasible overlap does at 11) and no simulated TU has a 4-bp overlap.

Origin: a TGA whose bases are n-2, n-1 and 0 of a circular contig is read as W (4) / G (25) by
the six-frame counter and by the expressed-protein translation.
"""

from __future__ import annotations

import csv
import itertools
import json
import re
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

HIGH_SHARE = {4: 0.9, 25: 0.5}
FAMILY = {4: "W", 25: "G"}
STOPS_4_25 = frozenset({"TAA", "TAG"})
STARTS = frozenset({"ATG", "GTG", "TTG"})
_AA3 = dict(
    zip(
        "ACDEFGHIKLMNPQRSTVWY",
        re.findall("...", "AlaCysAspGluPheGlyHisIleLysLeuMetAsnProGlnArgSerThrValTrpTyr"),
        strict=True,
    )
)
_COMP = str.maketrans("ACGT", "TGCA")


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _codons(nt: str) -> list[str]:
    return [nt[i : i + 3] for i in range(0, len(nt) - len(nt) % 3, 3)]


def _code(table: int) -> dict[str, str]:
    from gippasim.codons import code

    return code(table)


def _translate(nt: str, table: int) -> str:
    """CDS protein: initiator read as M, the final stop dropped."""
    t = _code(table)
    return "M" + "".join(t[c] for c in _codons(nt)[1:-1])


def _fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = ""
        elif line:
            out[name] += line.strip()
    return out


def _tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _slice(seq: str, s0: int, e0: int, strand: str) -> str:
    n = len(seq)
    nt = "".join(seq[i % n] for i in range(s0, e0))
    return _rc(nt) if strand == "-" else nt


def _codon_at(seq: str, c0: int, strand: str) -> str:
    return _slice(seq, c0, c0 + 3, strand)


def _contig(f: dict[str, Any]) -> str:
    return str(f.get("contig") or f.get("sequence"))


def _span(f: dict[str, Any], n: int) -> tuple[int, int]:
    s0, e0 = int(f["start"]) - 1, int(f["stop"])
    return (s0, e0 + n) if e0 <= s0 else (s0, e0)


def _coding(data: dict[str, Any], pseudo: bool = False) -> Iterator[dict[str, Any]]:
    for f in data["features"]:
        if f["type"] in ("cds", "sorf") and (pseudo or not f.get("pseudo")):
            yield f


def _check_json(path: Path, seqs: dict[str, str], table: int) -> int:
    """Every intact coding feature against the genome; returns the TGA codons seen inside."""
    data = json.loads(path.read_text())
    assert data["genome"]["translation_table"] == table
    n_tga = 0
    for f in _coding(data, pseudo=True):
        assert "*" not in f["aa"], f["locus"]
        if f.get("pseudo"):
            continue
        seq = seqs[_contig(f)]
        s0, e0 = _span(f, len(seq))
        assert f["nt"] == _slice(seq, s0, e0, f["strand"]), f["locus"]
        cod = _codons(f["nt"])
        assert cod[-1] in STOPS_4_25, (f["locus"], cod[-1])
        assert not STOPS_4_25.intersection(cod[1:-1]), f["locus"]
        assert f["aa"] == _translate(f["nt"], table), f["locus"]
        n_tga += cod[1:-1].count("TGA")
    return n_tga


def _check_correct_peptides(
    truth: list[dict[str, str]], seqs: dict[str, str], table: int, readthrough: set[str]
) -> int:
    """Correct peptides: see ``_check_peptides``; returns those read over a TGA codon in an
    annotated gene."""
    correct = [r for r in truth if r["is_wrong"] == "false"]
    assert correct
    return _check_peptides(correct, seqs, table, readthrough, genes_only=True)


def _check_peptides(
    rows: list[dict[str, str]],
    seqs: dict[str, str],
    table: int,
    readthrough: set[str] = frozenset(),  # type: ignore[assignment]
    protein: bool = True,
    genes_only: bool = False,
) -> int:
    """Each peptide is the translation under ``table`` of its truth interval, except (protein
    products only) an initiator M over GTG/TTG and a read-through Q over TAA/TAG (both are not
    six-frame residues); a peptide with neither has ``n_loci_genome >= 1``. Returns the peptides
    read over a TGA codon (in annotated ``GSIM_`` genes only with ``genes_only``)."""
    t = _code(table)
    n_tga = 0
    for r in rows:
        pep = r["peptide"]
        assert "*" not in pep, r
        nt = _slice(seqs[r["contig"]], int(r["start"]) - 1, int(r["end"]), r["strand"])
        cod = _codons(nt)
        assert len(cod) == len(pep), r
        special = False
        for i, (c, a) in enumerate(zip(cod, pep, strict=True)):
            if t[c] == a or {t[c], a} == {"I", "L"}:
                continue
            if (protein and i == 0 and a == "M" and c in STARTS) or (
                a == "Q" and c in STOPS_4_25 and r["locus_tag"] in readthrough
            ):
                special = True
            else:
                raise AssertionError((r, i, c, a))
        if not special:
            assert int(r["n_loci_genome"]) >= 1, r
        if not genes_only or r["locus_tag"].startswith("GSIM_"):
            n_tga += "TGA" in cod
    return n_tga


_ORF_ID = re.compile(
    r"GIPPAORF_(intergenic|cds_outframe|antisense|nterm_extension)_(.+)_(\d+)_(\d+)_([pm])"
)


def _check_orf_entries(peptides: Path, seqs: dict[str, str], topology: dict[str, str]) -> int:
    """Every ``GIPPAORF_`` stop-to-stop entry named in ``peptides.tsv`` (not the nterm ones)
    has no TAA/TAG inside and is bounded by TAA/TAG (or a linear contig end) under 4/25: an
    entry cut at a TGA would end before a sense codon. Returns the entries holding a TGA."""
    seen: set[str] = set()
    for r in _tsv(peptides):
        for pid in r["protein_ids"].split(";"):
            m = _ORF_ID.fullmatch(pid.removeprefix("rev_"))
            if m and m.group(1) != "nterm_extension":
                seen.add(m.group(0))
    assert seen
    n_tga = 0
    for e in seen:
        m = _ORF_ID.fullmatch(e)
        assert m is not None
        cid, s1, e1, strand = m.group(2), int(m.group(3)), int(m.group(4)), m.group(5)
        seq = seqs[cid]
        n, circ = len(seq), topology[cid] == "circular"
        t = seq if strand == "p" else _rc(seq)
        lo, hi = (s1 - 1, e1) if strand == "p" else (n - e1, n - s1 + 1)
        cod = _codons(_slice(t, lo, hi, "+"))
        assert not STOPS_4_25.intersection(cod), e
        if circ or hi + 3 <= n:
            assert _slice(t, hi, hi + 3, "+") in STOPS_4_25, e
        if circ or lo >= 3:
            assert _slice(t, lo - 3, lo, "+") in STOPS_4_25, e
        n_tga += "TGA" in cod
    return n_tga


def _protein(nt: str, table: int) -> str:
    """Translation up to the first stop of ``table``; a valid start codon 1 reads M."""
    t, out = _code(table), []
    for i, c in enumerate(_codons(nt)):
        if c in STOPS_4_25:
            break
        out.append("M" if i == 0 and c in STARTS else t[c])
    return "".join(out)


def _residue(codon: str, first: bool, table: int) -> str:
    return "M" if first and codon in STARTS else _code(table)[codon]


def _run(tmp: Path, table: int, monkeypatch: pytest.MonkeyPatch, **kw: Any) -> tuple[Any, list]:
    """``simulate`` with ``kw``; also returns the read-through genes it picked."""
    import importlib

    S = importlib.import_module("gippasim.simulate")  # the module (the package exports a function)
    picked: list[Any] = []
    real = S.pick_readthrough

    def spy(*a: Any, **k: Any) -> Any:
        out = real(*a, **k)
        picked.extend(out)
        return out

    monkeypatch.setattr(S, "pick_readthrough", spy)
    cfg = S.SimConfig(genetic_code=table, tga_sense_share=HIGH_SHARE[table], **kw)
    return S.simulate(cfg, 11, tmp), picked


@pytest.mark.parametrize("table", [4, 25])
def test_m1_outputs_follow_table(
    table: int, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippasim.proteomics import ProteomicsConfig
    from gippasim.strain import StrainConfig

    t = _code(table)
    res, picked = _run(
        tmp_path / "a", table, monkeypatch,
        genome_length=120_000, origin_crossing_gene=True, n_plasmids=1,
        plasmid_length=(8_000, 8_000), linear_contig_length=6_000, n_trna=150, n_sorfs=5,
        n_pseudogenes=4, n_nterm_extensions=5,
        n_novel_orfs=(("intergenic", 6), ("cds_outframe", 6), ("antisense", 6)),
        m7_truth=True, n_readthrough=20, strain=StrainConfig(snv_rate_per_kb=5.0, n_stop_loss=20),
        proteomics=ProteomicsConfig(),
    )  # fmt: skip
    out = res.outdir
    sample, ref = _fasta(res.genome_fasta), _fasta(out / "reference.fasta")
    sample_json = json.loads(res.bakta_json.read_text())
    ref_json = json.loads((out / "reference_bakta_like.json").read_text())

    # Bakta JSON, sample and reference
    assert _check_json(res.bakta_json, sample, table) > 0
    assert _check_json(out / "reference_bakta_like.json", ref, table) > 0

    # proteins.faa
    faa = _fasta(out / "proteins.faa")
    assert faa and not any("*" in s for s in faa.values())
    ext = {r["feature_id"] for r in _tsv(out / "truth" / "features.tsv") if r["true_start"]}
    for f in _coding(sample_json):
        if f["locus"] in faa and f["locus"] not in ext:
            assert faa[f["locus"]] == f["aa"], f["locus"]

    # peptides
    truth = _tsv(out / "truth" / "peptides.tsv")
    rt_loci = {p.locus for p in picked}
    tga_pep = _check_correct_peptides(truth, sample, table, rt_loci)
    assert tga_pep > 0, "no correct peptide over a TGA codon"
    assert not any("*" in r["peptide"] for r in _tsv(res.peptides))

    # novel ORFs
    novel = _tsv(out / "truth" / "novel_orfs.tsv")
    assert novel
    for r in novel:
        nt = _slice(sample[r["contig"]], int(r["start"]) - 1, int(r["end"]), r["strand"])
        cod = _codons(nt)
        assert cod[-1] in STOPS_4_25, r
        assert not STOPS_4_25.intersection(cod[1:-1]), r
        assert r["aa"] == _translate(nt, table), r

    # alt sites (SNVs only: sample and reference coordinates agree): exactly one row per SNV
    # codon of a coding gene whose residue changes, read under the table, no stop involved
    genes = {f["locus"]: f for f in _coding(sample_json, pseudo=True)}
    diffs = _tsv(out / "truth" / "strain_differences.tsv")
    want: dict[tuple[str, int, str, str], tuple[str, str]] = {}
    for r in diffs:
        if r["type"] != "snv" or r["gene"] not in genes:
            continue
        g, cid = genes[r["gene"]], r["contig"]
        n = len(sample[cid])
        s0, e0 = _span(g, n)
        p0 = int(r["pos_sample"]) - 1
        p0 = p0 + n if p0 < s0 else p0
        c0 = s0 + 3 * ((p0 - s0) // 3) if g["strand"] == "+" else e0 - 3 * ((e0 - 1 - p0) // 3 + 1)
        first = c0 == (s0 if g["strand"] == "+" else e0 - 3)
        ra = _residue(_codon_at(ref[cid], c0, g["strand"]), first, table)
        aa = _residue(_codon_at(sample[cid], c0, g["strand"]), first, table)
        if ra != aa and "*" not in (ra, aa):
            want[(cid, c0 % n + 1, g["strand"], r["gene"])] = (ra, aa)
    alt = _tsv(out / "truth" / "alt_sites.tsv")
    got = {
        (r["contig"], int(r["codon_start"]), r["strand"], r["gene"]): (r["ref_aa"], r["alt_aa"])
        for r in alt
        if r["gene"] in genes
    }
    assert got == want
    tga_alt = sum(
        "TGA" in (_codon_at(ref[c], s - 1, st), _codon_at(sample[c], s - 1, st))
        for c, s, st, _ in got
    )
    assert tga_alt > 0, "no alt site at a TGA codon"

    # stop losses
    sl = _tsv(out / "truth" / "stop_loss.tsv")
    assert sl
    for r in sl:
        assert r["sample_codon"] in STOPS_4_25 and r["reference_codon"] not in STOPS_4_25, r
        assert r["reference_residue"] == t[r["reference_codon"]], r
    assert any(r["reference_codon"] == "TGA" and r["planted"] == "true" for r in sl), (
        "no planted stop loss to TGA"
    )
    ref_by = {f["locus"]: f for f in _coding(ref_json, pseudo=True)}
    n_ct = 0
    for r in sl:
        if r["planted"] == "true" and r["cterm_peptide"] and r["locus_tag"] in ref_by:
            prot = _protein(ref_by[r["locus_tag"]]["nt"], table)
            assert prot.endswith(r["cterm_peptide"]), (r, prot[-40:])
            n_ct += 1
    assert n_ct > 0

    # SNV aa_change: <reference residue><codon><sample residue>
    ref_genes = {f["locus"]: f for f in _coding(ref_json, pseudo=True)}
    tga_snv = 0
    for r in _tsv(out / "truth" / "strain_differences.tsv"):
        m = re.fullmatch(r"([A-Z*])(\d+)([A-Z*])", r["aa_change"])
        if r["type"] != "snv" or not m or r["gene"] not in ref_genes:
            continue
        i = int(m.group(2)) - 1
        a = _codons(genes[r["gene"]]["nt"])[i]
        b = _codons(ref_genes[r["gene"]]["nt"])[i]
        assert m.group(1) == _residue(b, i == 0, table), (r, b)
        assert m.group(3) == _residue(a, i == 0, table), (r, a)
        tga_snv += "TGA" in (a, b)
    assert tga_snv > 0, "no SNV aa_change at a TGA codon"

    # read-through extensions
    assert picked
    tga_rt = 0
    for p in picked:
        assert p.codon in STOPS_4_25 and p.residue == "Q"
        seq, step = sample[p.contig], 3 if p.strand == "+" else -3
        cods = [_codon_at(seq, p.codon_start + step * j, p.strand) for j in range(1, p.k)]
        assert p.extension == "".join(t[c] for c in cods), p
        assert "*" not in p.extension
        assert _codon_at(seq, p.ext_stop, p.strand) in STOPS_4_25
        tga_rt += "TGA" in cods
    assert tga_rt > 0, "no read-through extension with a TGA codon"

    # tRNAs
    trnas = [f for f in sample_json["features"] if f["type"] == "tRNA"]
    assert trnas
    for f in trnas:
        assert f["amino_acid"] == _AA3[t[_rc(f["anti_codon"])]], f
    assert any(f["anti_codon"] == "TCA" for f in trnas), "no tRNA for TGA"

    # rule 20: pure read-through genes (reference protein = sample protein + the residue read
    # at the stop) are planted first, a TGA (W/G) sense codon included
    res, _ = _run(
        tmp_path / "s", table, monkeypatch, genome_length=120_000,
        strain=StrainConfig(snv_rate_per_kb=0.3, n_stop_loss=15), proteomics=ProteomicsConfig(),
    )  # fmt: skip
    sg = {f["locus"]: f for f in _coding(json.loads(res.bakta_json.read_text()))}
    rg = {f["locus"]: f for f in _coding(json.loads((res.outdir / "reference_bakta_like.json")
                                                     .read_text()), pseudo=True)}  # fmt: skip
    planted = [r for r in _tsv(res.outdir / "truth" / "stop_loss.tsv") if r["planted"] == "true"]
    assert len(planted) == 15
    for r in planted:
        loc = r["locus_tag"]
        assert _protein(rg[loc]["nt"], table) == sg[loc]["aa"] + t[r["reference_codon"]], r
    assert any(r["reference_codon"] == "TGA" for r in planted), "no pure stop loss to TGA"

    # indels: aa_change <reference residues><position><sample residues> from the two proteins
    res, _ = _run(
        tmp_path / "i", table, monkeypatch, genome_length=60_000,
        strain=StrainConfig(snv_rate_per_kb=1.0, indel_rate_per_kb=1.0),
    )  # fmt: skip
    out = res.outdir
    _check_json(out / "reference_bakta_like.json", _fasta(out / "reference.fasta"), table)
    sg = {f["locus"]: f for f in _coding(json.loads(res.bakta_json.read_text()), pseudo=True)}
    rg = {f["locus"]: f for f in _coding(json.loads((out / "reference_bakta_like.json")
                                                     .read_text()), pseudo=True)}  # fmt: skip
    n_indel = 0
    for r in _tsv(out / "truth" / "strain_differences.tsv"):
        if r["type"] == "snv" or not r["aa_change"] or r["gene"] not in sg or r["gene"] not in rg:
            continue
        a, b = _protein(rg[r["gene"]]["nt"], table), _protein(sg[r["gene"]]["nt"], table)
        i = 0
        while i < min(len(a), len(b)) and a[i] == b[i]:
            i += 1
        j = 0
        while j < min(len(a), len(b)) - i and a[-1 - j] == b[-1 - j]:
            j += 1
        assert r["aa_change"] == f"{a[i : len(a) - j] or '-'}{i + 1}{b[i : len(b) - j] or '-'}", r
        n_indel += 1
    assert n_indel > 0


@pytest.mark.parametrize("table", [4, 25])
def test_m1_v2_blocks_and_sixframe_wrong_follow_table(
    table: int, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippasim.proteomics import ProteomicsConfig

    res, _ = _run(
        tmp_path / "b", table, monkeypatch,
        genome_length=120_000, oric=True, prophage_bp=12_000, bgc_bp=10_000, crispr=True,
        rbs_strength=0.8, promoter_strength=0.8, terminator_strength=0.8, write_precomputed=True,
        proteomics=ProteomicsConfig(pi_wrong=0.05, wrong_source="sixframe_noncoding"),
    )  # fmt: skip
    out = res.outdir
    data = json.loads(res.bakta_json.read_text())
    seqs = {s["id"]: s.get("nt") or s.get("sequence") for s in data["sequences"]}
    _check_json(res.bakta_json, seqs, table)
    assert any(f.get("gene") == "mnmG" for f in data["features"])

    regions = _tsv(out / "truth" / "regions.tsv")
    proph = [
        (int(r["start"]) - 1, int(r["end"])) for r in regions if r["region_type"] == "prophage"
    ]
    assert proph
    in_proph = [
        f for f in _coding(data)
        if any(a <= int(f["start"]) - 1 and int(f["stop"]) <= b for a, b in proph)
    ]  # fmt: skip
    assert in_proph
    assert sum(_codons(f["nt"])[1:-1].count("TGA") for f in in_proph) > 0

    topology = {s["id"]: s["topology"] for s in data["sequences"]}
    _check_wrong(out, seqs, topology, table)

    # decoys = "scored_null": six-frame nulls and protein-entry nulls (low detection leaves
    # undetected products; the strain makes some recorded loci pseudo, which are nulls too)
    from gippasim.strain import StrainConfig

    res, _ = _run(
        tmp_path / "c", table, monkeypatch, genome_length=120_000, n_nterm_extensions=5,
        strain=StrainConfig(snv_rate_per_kb=2.0, indel_rate_per_kb=1.0),
        proteomics=ProteomicsConfig(pi_wrong=0.1, decoys="scored_null",
                                    detect_logit=(-4.0, 0.08, 0.9)),
    )  # fmt: skip
    data = json.loads(res.bakta_json.read_text())
    seqs = _fasta(res.genome_fasta)
    topology = {s["id"]: s["topology"] for s in data["sequences"]}
    assert _check_wrong(res.outdir, seqs, topology, table) > 0


def _check_wrong(out: Path, seqs: dict[str, str], topology: dict[str, str], table: int) -> int:
    truth = _tsv(out / "truth" / "peptides.tsv")
    assert _check_correct_peptides(truth, seqs, table, set()) > 0
    assert not any("*" in r["peptide"] for r in truth)
    six = [r for r in truth if r["is_wrong"] == "true" and r["origin"] == "sixframe_noncoding"]
    assert six
    assert _check_peptides(six, seqs, table, protein=False) > 0, "no six-frame null over a TGA"
    nulls = [r for r in truth if r["origin"] == "undetected_product"]
    _check_peptides(nulls, seqs, table)
    assert _check_orf_entries(out / "peptides.tsv", seqs, topology) > 0
    return len(nulls)


# ------------------------------------------------------------------ L1: intra-TU overlaps


@pytest.mark.parametrize("table", [4, 25])
def test_l1_four_bp_overlap_infeasible_at_4_and_25(table: int) -> None:
    from gippasim.simulate import SimConfig

    SimConfig(genetic_code=11, overlap_bp=((4, 1.0),))  # feasible at 11
    with pytest.raises(ValueError):
        SimConfig(genetic_code=table, overlap_bp=((4, 1.0),))
    SimConfig(genetic_code=table, overlap_bp=((1, 0.5), (8, 0.5)))  # others stay feasible


def _overlaps(genes_tsv: Path) -> list[int]:
    rows = _tsv(genes_tsv)
    rows.sort(key=lambda r: (r["tu_id"], int(r["start"])))
    out = []
    for a, b in itertools.pairwise(rows):
        if a["tu_id"] and a["tu_id"] == b["tu_id"]:
            k = int(a["end"]) - (int(b["start"]) - 1)
            if k > 0:
                out.append(k)
    return out


@pytest.mark.parametrize("table", [4, 11, 25])
def test_l1_no_four_bp_overlap_in_simulated_tus(table: int, tmp_path: Path) -> None:
    from gippasim.simulate import SimConfig, simulate

    ks: list[int] = []
    for seed in range(4):
        res = simulate(
            SimConfig(genetic_code=table, genome_length=100_000), seed, tmp_path / str(seed)
        )
        ks += _overlaps(res.genes_truth)
    assert ks
    if table == 11:
        assert 4 in ks  # sanity: the count works
    else:
        assert 4 not in ks


# ------------------------------------------------------------------ TGA across the origin

#: Codons of a gene whose 6th codon is TGA; rotated so that TGA occupies bases n-2, n-1, 0.
_GENE = "ATG" + "AAACTGTGGAAAGAA" + "TGA" + "AAAAAGCTGGAAAAACTG" + "TAA"
_TGA_AT = _GENE.index("TGA", 3)  # 18
_FILLER = "CCCCCCCCCCCC"


def _origin_contig() -> tuple[str, int]:
    """Contig and the 0-based start of the gene (crossing the origin)."""
    seq = _GENE[_TGA_AT + 2 :] + _FILLER + _GENE[: _TGA_AT + 2]
    n = len(seq)
    assert seq[n - 2 :] + seq[0] == "TGA"
    return seq, n - (_TGA_AT + 2)


@pytest.mark.parametrize("table", [4, 25])
def test_origin_spanning_tga_sixframe(table: int) -> None:
    from gippasim.proteomics import sixframe_counts

    seq, _ = _origin_contig()
    pep = _translate(_GENE, table)[1:]  # KLWKE + W/G + KKLEKL
    assert FAMILY[table] in pep[5]
    for contig in (seq, _rc(seq)):  # plus strand, then the same gene on the minus strand
        got = sixframe_counts([pep], {"c": contig}, {"c": "circular"}, table)
        assert got[pep] == 1
        assert sixframe_counts([pep], {"c": contig}, {"c": "linear"}, table)[pep] == 0
    other = 4 if table == 25 else 25
    assert sixframe_counts([_translate(_GENE, other)[1:]], {"c": seq}, {"c": "circular"},
                           table)[_translate(_GENE, other)[1:]] == 0  # fmt: skip


@pytest.mark.parametrize("table", [4, 25])
def test_origin_spanning_tga_expressed_protein(table: int) -> None:
    from gippasim.features import Feat
    from gippasim.proteomics import expressed_proteins

    seq, s0 = _origin_contig()
    f = Feat("cds", contig="c", start=s0, end=s0 + len(_GENE), strand="+", locus="G1")
    (p,) = expressed_proteins([f], [], {"c": seq}, None, table=table)
    assert p.aa == _translate(_GENE, table)
    assert p.aa[6] == FAMILY[table]
