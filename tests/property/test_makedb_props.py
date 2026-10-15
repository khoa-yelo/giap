"""Hypothesis properties of `giap makedb` (issue M3-23; GIAP_SPEC.md §7.9.2, §6.6).

Random small genomes (1-2 contigs, circular or linear, 150-900 nt, AT-rich so ORFs stay short)
with up to three planted genes (CDS or sORF; ATG/GTG/TTG starts; either strand; origin-crossing
on circular contigs) and a random ``min_orf_aa``:

1. six-frame ORFs and ``nterm_extension`` entries equal the brute-force oracle
   (``tests/helpers/makedb_oracle.py``), in the spec's order;
2. annotated proteins are written once each, in gene-model order, and ``db_id_map.tsv`` describes
   every record (one row per record, header coordinates = id_map coordinates);
3. acceptance criterion: every in-silico tryptic peptide (7-40 aa) of every record, except
   contaminants and the junction part of ``nterm_extension`` entries (and the initiator ``M`` of
   an alternative start), maps exactly with M3-04 to the locus given in ``db_id_map.tsv``;
4. reruns are byte-identical.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

from giap.annotate.bakta import parse_bakta_json
from giap.genetable.covariates import digest
from giap.mapping.decoys import build_decoys
from giap.mapping.exact import match_exact
from giap.mapping.sixframe import build_sixframe
from giap.model.coords import Replicon
from giap.rng import for_step

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))

from makedb_checks import (
    assert_orfs_match_oracle,
    by_id,
    fasta_records,
    id_map,
)
from makedb_genomes import Contig, feature, write_bakta_json
from makedb_oracle import (
    STOPS,
    DegenerateCircle,
    all_orfs,
    nterm_entries,
    revcomp,
    translate,
)

SETTINGS = settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
_SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS]


def _build(annotation: Path, outdir: Path, **kw: Any) -> Path:
    from giap.makedb.build import MakedbOptions, build_strain_db

    build_strain_db(MakedbOptions(annotation=annotation, genome=None, outdir=outdir, **kw))
    return outdir


@st.composite
def genomes(draw: st.DrawFn) -> tuple[list[Contig], list[dict[str, Any]]]:
    contigs: list[Contig] = []
    feats: list[dict[str, Any]] = []
    tag = 0
    for k in range(draw(st.integers(1, 2))):
        n = draw(st.integers(150, 900))
        circular = draw(st.booleans())
        rng = np.random.default_rng(draw(st.integers(0, 2**32 - 1)))
        seq = list("".join(rng.choice(np.array(list("AAATTTGC")), n).tolist()))
        c = Contig(
            f"c{k + 1}",
            n,
            "circular" if circular else "linear",
            "chromosome" if k == 0 else "plasmid",
            seq=seq,
        )
        taken: list[tuple[int, int]] = []
        for _ in range(draw(st.integers(0, 3))):
            n_codons = draw(st.integers(6, 45))
            length = 3 * n_codons
            if length >= n // 2:
                continue
            start = draw(st.integers(0, n - 1 if circular else n - length))
            end = start + length
            if not circular and end > n:
                continue
            span = [(start, end)] if end <= n else [(start, n), (0, end - n)]
            if any(a < te and ta < b for a, b in span for ta, te in taken):
                continue
            taken.extend(span)
            strand = draw(st.sampled_from("+-"))
            start_codon = draw(st.sampled_from(["ATG", "ATG", "GTG", "TTG"]))
            body = [_SENSE[int(i)] for i in rng.integers(0, len(_SENSE), n_codons - 2)]
            nt = start_codon + "".join(body) + draw(st.sampled_from(sorted(STOPS)))
            c.put(start, nt if strand == "+" else revcomp(nt))
            aa = "M" + translate("".join(body))
            kind = "sorf" if n_codons <= 30 and draw(st.booleans()) else "cds"
            tag += 5
            feats.append(
                feature(
                    kind,
                    c.id,
                    start,
                    end,
                    strand,
                    locus=f"PROP_{tag:05d}",
                    aa=aa,
                    product=draw(st.sampled_from(["", "putative protein", "Ftsz a,b"])),
                    start_type=start_codon,
                    length=n,
                )
            )
        contigs.append(c)
    return contigs, feats


def _setup(tmp: Path, data: tuple[list[Contig], list[dict[str, Any]]]) -> tuple[Path, Any]:
    contigs, feats = data
    for c in contigs:
        try:
            all_orfs(c.id, c.text, Replicon(c.id, c.length, c.topology))  # type: ignore[arg-type]
        except DegenerateCircle:
            assume(False)
    annotation = write_bakta_json(tmp / "bakta.json", contigs, feats)
    return annotation, parse_bakta_json(annotation)


@SETTINGS
@given(data=genomes(), min_aa=st.integers(5, 25))
def test_orfs_equal_the_oracle(data: Any, min_aa: int) -> None:
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        annotation, ann = _setup(tmp, data)
        out = _build(annotation, tmp / "DB", six_frame_orfs=True, min_orf_aa=min_aa)
        assert_orfs_match_oracle(out, annotation, dict(ann.sequences), min_aa)


@SETTINGS
@given(data=genomes(), min_aa=st.integers(5, 25), orfs=st.booleans())
def test_annotated_records_and_id_map(data: Any, min_aa: int, orfs: bool) -> None:
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        annotation, ann = _setup(tmp, data)
        out = _build(annotation, tmp / "DB", six_frame_orfs=orfs, min_orf_aa=min_aa)
        recs = fasta_records(out / "strain_db.faa")
        rows = id_map(out)
        assert [h.split(" ", 1)[0] for h, _ in recs] == [r["db_id"] for r in rows]
        coding = [g for g in ann.model.genes if g.feature_type in ("cds", "sorf")]
        annotated = [
            (h, s) for (h, s), r in zip(recs, rows, strict=True) if r["kind"] in ("cds", "sorf")
        ]
        assert annotated == [
            (f"{g.locus_tag} {g.product or 'hypothetical protein'}", g.aa) for g in coding
        ]
        assert recs[: len(coding)] == annotated
        if not orfs:
            assert len(recs) == len(coding)
        for (_, seq), r in zip(recs, rows, strict=True):
            assert seq and "*" not in seq
            if r["kind"] == "orf":
                cls = r["class"]
                sfx = "p" if r["strand"] == "+" else "m"
                assert r["db_id"] == f"GIAPORF_{cls}_{r['contig']}_{r['start']}_{r['end']}_{sfx}"
                assert int(r["start"] or 0) <= int(r["end"] or 0)
                if cls != "nterm_extension":
                    assert len(seq) >= min_aa and "X" not in seq
                    assert r["locus_tag"] is None


def _contains(rep: Replicon, s0: int, e0: int, hs: int, he: int) -> bool:
    shifts = (0, rep.length, -rep.length) if rep.topology == "circular" else (0,)
    return any(s0 <= hs + d and he + d <= e0 for d in shifts)


@SETTINGS
@given(data=genomes(), min_aa=st.integers(7, 25))
def test_tryptic_peptides_map_to_the_id_map_locus(data: Any, min_aa: int) -> None:
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        annotation, ann = _setup(tmp, data)
        out = _build(annotation, tmp / "DB", six_frame_orfs=True, min_orf_aa=min_aa)
        gm, genome = ann.model, dict(ann.sequences)
        idx = build_sixframe(genome, gm.replicons)
        decoys = {"reverse": build_decoys(idx, "reverse", 50, for_step("sixframe_index", 0))}
        nterms, undefined = nterm_entries(gm, genome)
        nterm_up = {n.locus_tag: len(n.aa) - n.n_protein for n in nterms}
        jobs: list[tuple[str, dict[str, str | None]]] = []
        seqs = by_id(out)
        for r in id_map(out):
            seq = seqs[r["db_id"] or ""]
            if r["kind"] == "contaminant":
                continue
            lo = 0
            if r["kind"] in ("cds", "sorf"):
                g = gm.by_locus(r["locus_tag"] or "")
                lo = 0 if g.start_codon == "ATG" else 1  # initiator M of an alternative start
            hi = len(seq)
            if r["class"] == "nterm_extension":
                if r["locus_tag"] in undefined:  # lead decision 4: not asserted
                    continue
                hi = nterm_up[r["locus_tag"] or ""]  # upstream part only (junction excluded)
            pos = 0
            for frag in digest(seq, "trypsin"):
                a, b = pos, pos + len(frag)
                pos = b
                if a >= lo and b <= hi and 7 <= len(frag) <= 40:
                    jobs.append((frag, r))
        assume(jobs)
        peps = sorted({p for p, _ in jobs})
        hits = match_exact(peps, idx, decoys, {p: "reverse" for p in peps})
        hits = hits[hits["search"] == "target"]
        by_pep: dict[str, list[tuple[str, int, int, str]]] = {}
        for row in hits.itertuples():
            by_pep.setdefault(row.peptide, []).append((row.contig, row.start, row.end, row.strand))
        for pep, r in jobs:
            rep = gm.replicons[r["contig"] or ""]
            s0, e0 = int(r["start"] or 0) - 1, int(r["end"] or 0)
            assert any(
                c == r["contig"] and st_ == r["strand"] and _contains(rep, s0, e0, hs, he)
                for c, hs, he, st_ in by_pep.get(pep, [])
            ), (pep, r)


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(data=genomes(), min_aa=st.integers(5, 25))
def test_reruns_are_byte_identical(data: Any, min_aa: int) -> None:
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        annotation, _ = _setup(tmp, data)
        a = _build(annotation, tmp / "A", six_frame_orfs=True, min_orf_aa=min_aa)
        b = _build(annotation, tmp / "B", six_frame_orfs=True, min_orf_aa=min_aa)
        for name in ("strain_db.faa", "db_id_map.tsv"):
            assert (a / name).read_bytes() == (b / name).read_bytes()
