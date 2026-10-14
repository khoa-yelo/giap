"""giapsim review decisions for M3-18 (docs/specs/M3-18-sim-genome-features.md, "Lead decisions
after code review (PR #25)", items 1-7).

Every property is checked on ``N_SEEDS`` simulations of a 50-kb chromosome with every option on
(and dense strain differences), so rare events (a start codon hit by a SNV, an in-frame indel in a
gene, an indel next to a feature end) occur often enough to kill the mutations the review listed:
``aa_change`` swapped or empty, ``clean = True`` in the lift, the indel feature-boundary guard off.

Conventions pinned here where the decisions leave room:

* "clean" (decision 1): ``len(nt) % 3 == 0``, codon 1 in ATG/GTG/TTG, last codon a stop, no
  internal stop;
* reference ``aa`` (decision 2): table-11 translation of ``nt`` codons up to the first stop; codon 1
  reads as M only when it is a valid start (otherwise its table-11 residue: no protein the ``nt``
  does not encode);
* ``ref_pseudo`` (decision 2) is a column of ``truth/features.tsv`` (or, failing that,
  ``truth/genes.tsv``) keyed by locus; true is ``true``/``True``/``1``;
* ``aa_change`` (decision 4) for a SNV reads ``<reference residue><codon number><sample residue>``
  (VCF orientation, as in M3-18's strain table); the codon number may count in either genome.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pytest
from test_sim_features import (
    CODE,
    STARTS,
    STOPS,
    Run,
    all_options,
    by_type,
    fcontig,
    fid,
    ftype,
    is_pseudo,
    read_fasta,
    seq_of,
    simulate,
    strain_config,
)

N_SEEDS = 40
TRUE = {"true", "True", "TRUE", "1"}
FALSE = {"false", "False", "FALSE", "0", "", "NA", "."}
SNV_CHANGE = re.compile(r"^([A-Z*])(\d+)([A-Z*])$")


# ---------------------------------------------------------------- helpers


def options() -> dict[str, Any]:
    return all_options(50_000) | {
        "strain": strain_config(snv_rate_per_kb=4.0, indel_rate_per_kb=1.5),
    }


def clean(nt: str) -> bool:
    cod = [nt[i : i + 3] for i in range(0, len(nt) - len(nt) % 3, 3)]
    return (
        len(nt) % 3 == 0
        and len(cod) >= 2
        and cod[0] in STARTS
        and cod[-1] in STOPS
        and not STOPS.intersection(cod[1:-1])
    )


def residue(codon: str, index: int) -> str:
    """Table-11 residue of codon ``index`` (0-based); codon 1 reads as M when it is a start."""
    return "M" if index == 0 and codon in STARTS else CODE[codon]


def to_first_stop(nt: str) -> str:
    out = []
    for i in range(0, len(nt) - len(nt) % 3, 3):
        c = nt[i : i + 3]
        if c in STOPS:
            break
        out.append(residue(c, i // 3))
    return "".join(out)


class StrainRun:
    def __init__(self, run: Run) -> None:
        self.run = run
        self.ref = {k: v[1] for k, v in read_fasta(run.out / "reference.fasta").items()}
        self.rdata: dict[str, Any] = json.loads((run.out / "reference_bakta_like.json").read_text())
        self.rfeat = {fid(f): f for f in self.rdata["features"]}
        self.sfeat = {fid(f): f for f in run.data["features"]}
        _, self.rows = run.truth("strain_differences.tsv")

    def ref_pseudo(self) -> dict[str, bool]:
        for name in ("features.tsv", "genes.tsv"):
            header, rows = self.run.truth(name)
            if "ref_pseudo" in header:
                key = "feature_id" if "feature_id" in header else "locus_tag"
                out = {}
                for r in rows:
                    v = r["ref_pseudo"]
                    assert v in TRUE | FALSE, (name, r[key], v)
                    out[r[key]] = v in TRUE
                return out
        raise AssertionError("no ref_pseudo column in truth/features.tsv or truth/genes.tsv")


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory) -> list[StrainRun]:
    root = tmp_path_factory.mktemp("review")
    return [StrainRun(simulate(root / str(s), 100 + s, **options())) for s in range(N_SEEDS)]


def _offset(f: dict[str, Any], length: int, pos0: int) -> int | None:
    """Offset of genome position ``pos0`` in feature ``f`` in transcription direction."""
    s, e = int(f["start"]) - 1, int(f["stop"])
    if s >= e:
        e += length
    for p in (pos0, pos0 + length):
        if s <= p < e:
            return p - s if f["strand"] == "+" else e - 1 - p
    return None


def _apply(seq: str, pos: int, ref: str, alt: str) -> str:
    assert seq[pos - 1 : pos - 1 + len(ref)] == ref, (pos, ref)
    return seq[: pos - 1] + alt + seq[pos - 1 + len(ref) :]


# ---------------------------------------------------------------- decisions 1 and 2


def test_reference_non_pseudo_cds_is_clean_with_valid_start(runs: list[StrainRun]) -> None:
    n = 0
    for sr in runs:
        for f in sr.rdata["features"]:
            if ftype(f) in ("cds", "sorf") and not is_pseudo(f):
                nt = seq_of(f, sr.ref)
                assert f["nt"] == nt, fid(f)
                assert nt[:3] in STARTS, (sr.run.out, fid(f), nt[:3])
                assert clean(nt), (sr.run.out, fid(f))
                n += 1
    assert n > 1000


def test_reference_aa_is_translation_of_reference_nt(runs: list[StrainRun]) -> None:
    for sr in runs:
        for f in sr.rdata["features"]:
            if ftype(f) in ("cds", "sorf"):
                assert f["aa"].rstrip("*") == to_first_stop(f["nt"]), (sr.run.out, fid(f))


def test_ref_pseudo_set_exactly_when_lift_breaks_cds(runs: list[StrainRun]) -> None:
    broken = start_lost = 0
    for sr in runs:
        flags = sr.ref_pseudo()
        for locus, f in sr.sfeat.items():
            if ftype(f) not in ("cds", "sorf"):
                assert not flags.get(locus, False), locus
                continue
            r = sr.rfeat[locus]
            if is_pseudo(f):
                assert is_pseudo(r), locus  # a sample pseudogene stays one
                continue
            ok = clean(r["nt"])
            assert flags[locus] is (not ok), (sr.run.out, locus)
            assert is_pseudo(r) is (not ok), (sr.run.out, locus)
            broken += not ok
            start_lost += r["nt"][:3] not in STARTS
    assert broken >= 20
    assert start_lost >= 1  # decision 1 is exercised: a lost start alone breaks the CDS


# ---------------------------------------------------------------- decision 3: indels


def test_alleles_anchored_and_reapply_to_sample(runs: list[StrainRun]) -> None:
    n_indel = 0
    for sr in runs:
        rebuilt = dict(sr.ref)
        for r in sorted(sr.rows, key=lambda r: -int(r["pos_reference"])):
            ref, alt = r["ref"], r["alt"]
            assert ref and alt and set(ref + alt) <= set("ACGT"), r
            if r["type"] == "snv":
                assert len(ref) == len(alt) == 1 and ref != alt, r
            else:
                n_indel += 1
                assert ref[0] == alt[0] and min(len(ref), len(alt)) == 1, r
                assert (r["type"] == "ins") == (len(alt) > len(ref)), r
            p_s = int(r["pos_sample"])
            assert sr.run.genomes[r["contig"]][p_s - 1 : p_s - 1 + len(alt)] == alt, r
            c = r["contig"]
            rebuilt[c] = _apply(rebuilt[c], int(r["pos_reference"]), ref, alt)
        assert rebuilt == sr.run.genomes, sr.run.out
    assert n_indel >= 200


def test_indels_left_normalized(runs: list[StrainRun]) -> None:
    for sr in runs:
        for r in sr.rows:
            if r["type"] == "snv":
                continue
            seq, p, ref, alt = sr.ref[r["contig"]], int(r["pos_reference"]), r["ref"], r["alt"]
            assert p >= 2, r
            here = _apply(seq, p, ref, alt)
            left = _apply(seq, p - 1, seq[p - 2 : p - 2 + len(ref)], seq[p - 2] + alt[1:])
            assert here != left, ("indel can shift left", r)


def test_indels_stay_inside_or_outside_features(runs: list[StrainRun]) -> None:
    """The boundary guard: no indel touches a feature end, so features lift cleanly."""
    for sr in runs:
        lengths = {c: len(s) for c, s in sr.run.genomes.items()}
        net: dict[str, int] = {}
        for r in sr.rows:
            if r["type"] == "snv":
                continue
            c, a0 = r["contig"], int(r["pos_sample"]) - 1  # anchor, 0-based
            lo, hi = a0, a0 + len(r["alt"])  # anchor .. first base after the event
            for f in sr.run.data["features"]:
                if fcontig(f) != c:
                    continue
                s, e = int(f["start"]) - 1, int(f["stop"])
                for b in (s, e % lengths[c]):
                    assert not lo <= b <= hi, ("indel at a feature end", fid(f), r)
                if _offset(f, lengths[c], a0) is not None:
                    net[fid(f)] = net.get(fid(f), 0) + len(r["ref"]) - len(r["alt"])
        for locus, f in sr.sfeat.items():
            lift = sr.rfeat[locus]
            span = len(seq_of(lift, sr.ref)) - len(seq_of(f, sr.run.genomes))
            assert span == net.get(locus, 0), (sr.run.out, locus)


# ---------------------------------------------------------------- decision 4: aa_change


def test_aa_change_matches_protein_difference(runs: list[StrainRun]) -> None:
    nonsyn = inframe = 0
    for sr in runs:
        lengths = {c: len(s) for c, s in sr.run.genomes.items()}
        for r in sr.rows:
            f = sr.sfeat.get(r["gene"])
            if f is None or is_pseudo(f):
                continue
            g = sr.rfeat[r["gene"]]
            o_s = _offset(f, lengths[r["contig"]], int(r["pos_sample"]) - 1)
            o_r = _offset(g, len(sr.ref[r["contig"]]), int(r["pos_reference"]) - 1)
            assert o_s is not None and o_r is not None, r
            if r["type"] != "snv":
                if (len(r["ref"]) - len(r["alt"])) % 3 == 0:
                    assert r["aa_change"] not in FALSE, ("in-frame indel without aa_change", r)
                    inframe += 1
                continue
            if (o_r - o_s) % 3:
                continue  # behind a frameshift: no codon-level comparison
            ci_s, ci_r = o_s // 3, o_r // 3
            cs, cr = f["nt"][3 * ci_s : 3 * ci_s + 3], g["nt"][3 * ci_r : 3 * ci_r + 3]
            if len(cr) < 3:
                continue
            a_s, a_r = residue(cs, ci_s), residue(cr, ci_r)
            if a_s == a_r:
                assert r["aa_change"] in FALSE | {f"{a_r}{ci_s + 1}{a_s}", f"{a_r}{ci_r + 1}{a_s}"}
                continue
            if ci_s == 0 and "M" in (a_s, a_r) and cr not in STARTS:
                continue  # start lost: residue convention for a broken start left open
            m = SNV_CHANGE.match(r["aa_change"])
            assert m, ("non-synonymous SNV without aa_change", r)
            assert (m.group(1), m.group(3)) == (a_r, a_s), (r, a_r, a_s)
            assert int(m.group(2)) in (ci_s + 1, ci_r + 1), r
            nonsyn += 1
    assert nonsyn >= 500 and inframe >= 20


def test_start_codon_snv_between_valid_starts_is_synonymous(runs: list[StrainRun]) -> None:
    """Codon 1 reads as M: ATG<->GTG<->TTG at codon 1 changes no residue."""
    for sr in runs:
        for r in sr.rows:
            f = sr.sfeat.get(r["gene"])
            if r["type"] != "snv" or f is None or is_pseudo(f):
                continue
            g = sr.rfeat[r["gene"]]
            o = _offset(f, len(sr.run.genomes[r["contig"]]), int(r["pos_sample"]) - 1)
            if o is not None and o < 3 and g["nt"][:3] in STARTS:
                assert r["aa_change"] in FALSE | {"M1M"}, r


# ---------------------------------------------------------------- decision 5: truth tables


def test_genes_truth_covers_every_cds(runs: list[StrainRun]) -> None:
    for sr in runs[:10]:
        run = sr.run
        _, rows = run.truth("genes.tsv")
        have = {r["locus_tag"]: r for r in rows}
        cds = {fid(f) for f in by_type(run.data, "cds")}
        coding = cds | {fid(f) for f in by_type(run.data, "sorf")}
        assert cds <= set(have), sorted(cds - set(have))[:5]
        assert set(have) <= coding
        contigs = {fcontig(run.feature(k)) for k in have}
        assert len(contigs) == len(run.genomes)  # chromosome, both plasmids, linear contig
        for k, r in have.items():
            f = run.feature(k)
            iv = run.interval(f)
            assert (int(r["start"]), r["strand"]) == (int(f["start"]), f["strand"]), k
            assert int(r["end"]) in (int(f["stop"]), iv.end), k
            if not is_pseudo(f):
                assert r["start_codon"] == f["nt"][:3], k
                assert int(r["length_aa"]) == len(f["aa"].rstrip("*")), k
        wrapped = [fid(f) for f in run.data["features"] if int(f["start"]) > int(f["stop"])]
        assert len(wrapped) == 1 and wrapped[0] in have
        planted = [f for f in by_type(run.data, "cds") if is_pseudo(f)]
        assert planted and all(fid(f) in have for f in planted)


def test_features_truth_end_unwrapped(runs: list[StrainRun]) -> None:
    for sr in runs[:10]:
        run = sr.run
        _, rows = run.truth("features.tsv")
        crossing = 0
        for r in rows:
            f = run.feature(r["feature_id"])
            iv = run.interval(f)
            assert (int(r["start"]), int(r["end"])) == (iv.start + 1, iv.end), r["feature_id"]
            crossing += iv.end > run.reps[r["contig"]].length
        assert crossing == 1


# ---------------------------------------------------------------- decision 6: edge, GFF3 sORFs


def test_origin_crossing_gene_has_edge(runs: list[StrainRun]) -> None:
    for sr in runs[:10]:
        for f in sr.run.data["features"]:
            crosses = int(f["start"]) > int(f["stop"])
            assert (f.get("edge") is True) is crosses, fid(f)


def _gff_rows(path: Path) -> list[tuple[str, dict[str, str]]]:
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("##FASTA"):
            break
        if not line or line.startswith("#"):
            continue
        col = line.split("\t")
        attrs = dict(kv.split("=", 1) for kv in col[8].split(";") if "=" in kv)
        out.append((col[2], attrs))
    return out


def test_gff3_sorfs_are_cds_with_sorf_attribute(runs: list[StrainRun]) -> None:
    for sr in runs[:10]:
        run = sr.run
        sorfs = {fid(f) for f in by_type(run.data, "sorf")}
        assert sorfs
        rows = _gff_rows(run.out / "annotation.gff3")
        assert not [t for t, _ in rows if t.lower() == "sorf"]
        seen = set()
        for t, a in rows:
            ident = a.get("locus_tag") or a.get("ID", "")
            if ident in sorfs:
                assert t == "CDS" and a.get("sorf") == "true", (ident, t, a)
                seen.add(ident)
            else:
                assert a.get("sorf") != "true", ident
        assert seen == sorfs


# ---------------------------------------------------------------- decision 7: hashable config


def test_simconfig_hashable() -> None:
    import giapsim

    assert hash(giapsim.SimConfig()) == hash(giapsim.SimConfig())
    opts = all_options(50_000)
    a, b = giapsim.SimConfig(**opts), giapsim.SimConfig(**opts)
    assert a == b and hash(a) == hash(b)
    assert {a: 1}[b] == 1
    assert dict(a.n_novel_orfs) == {"intergenic": 2, "cds_outframe": 2, "antisense": 2}
