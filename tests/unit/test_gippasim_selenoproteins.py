"""gippasim ``n_selenoproteins`` (M7-16c tooling; spec rule 23, example 32, parked default Q23 "Sec
only"; test plan C11).

``SimConfig.n_selenoproteins: int = 0``. With 0 nothing is drawn or edited and
``truth/selenoproteins.tsv`` is not written (every existing build byte-identical). With ``n > 0``
the simulator draws, from its own generator ``default_rng([seed, _SEC_STREAM])``, ``n`` distinct
eligible genes (annotated non-pseudo ``cds``, not crossing the origin, ``len(aa) >= 21``) and one
residue ``i`` per gene with ``10 <= i <= len(aa) - 11`` whose codon overlaps no other annotated
feature or simulator edit. That codon becomes TGA on the gene strand (only that codon changes),
residue ``i`` of the annotated protein becomes ``U`` (gene table, Bakta JSON, ``genes.parquet``)
and the expressed protein carries ``U`` there too. Truth: ``truth/selenoproteins.tsv``
(``locus_tag``, ``contig``, ``strand``, ``codon_start`` 0-based lowest coordinate, ``residue``)
sorted by (contig, ``codon_start``). Fewer eligible genes than ``n`` raise ``ValueError("n_
selenoproteins=<n> but only <k> eligible genes")``; ``n > 0`` with ``strain`` set raises
``ValueError``. No ``O``/TAG is planted (Q23).
"""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import replace
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from gippasim.codons import CODE_11
from gippasim.proteomics import ProteomicsConfig
from gippasim.simulate import SimConfig, simulate
from gippasim.strain import StrainConfig

SEED = 1
BASE = SimConfig(genome_length=100_000, proteomics=ProteomicsConfig(search_db="same"))
COLUMNS = ["locus_tag", "contig", "strand", "codon_start", "residue"]
#: digest of the float-free files of BASE at seed 1, recorded on 09eef19e (before the knob)
KEEP = (
    "genome.fasta",
    "bakta_like.json",
    "annotation.gff3",
    "proteins.faa",
    "truth/genes.tsv",
    "truth/features.tsv",
    "truth/peptides.tsv",
    "truth/novel_orfs.tsv",
)
PINNED_BASE = "f1c2a0cdba97f4a5febcf2f215f19472683580143a519e8a17793ebfa22d0aac"
_COMP = str.maketrans("ACGT", "TGCA")


def rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _digest(out: Path) -> str:
    h = hashlib.sha256()
    for rel in sorted(KEEP):
        p = out / rel
        if p.exists():
            h.update(rel.encode() + b"\0" + hashlib.sha256(p.read_bytes()).digest())
    return h.hexdigest()


def _files(out: Path) -> dict[str, bytes]:
    return {str(p.relative_to(out)): p.read_bytes() for p in sorted(out.rglob("*")) if p.is_file()}


def fasta(path: Path) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    cur = ""
    for ln in path.read_text().splitlines():
        if ln.startswith(">"):
            cur = ln[1:].split()[0]
            out[cur] = []
        elif ln:
            out[cur].append(ln.strip())
    return {k: "".join(v) for k, v in out.items()}


def features(out: Path) -> list[dict[str, Any]]:
    return json.loads((out / "bakta_like.json").read_text())["features"]  # type: ignore[no-any-return]


def interval(f: dict[str, Any], n: int) -> tuple[int, int]:
    """0-based unwrapped ``[s, e)`` of a Bakta-like feature (1-based closed, wrapped)."""
    s, e = int(f["start"]) - 1, int(f["stop"])
    return (s, e + n) if e <= s else (s, e)


def gene_codon(seq: str, s: int, e: int, strand: str, i: int) -> tuple[int, str]:
    n = len(seq)
    cs = s + 3 * i if strand == "+" else e - 3 * i - 3
    bases = "".join(seq[(cs + k) % n] for k in range(3))
    return cs % n, bases if strand == "+" else rc(bases)


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    tmp = tmp_path_factory.mktemp("sec")
    out = {}
    for name, cfg in (
        ("base", BASE),
        ("zero", replace(BASE, n_selenoproteins=0)),
        ("five", replace(BASE, n_selenoproteins=5)),
        ("twenty", replace(BASE, n_selenoproteins=20)),
    ):
        simulate(cfg, SEED, tmp / name)
        out[name] = tmp / name
    return out


def truth(out: Path) -> pd.DataFrame:
    df = pd.read_csv(
        out / "truth" / "selenoproteins.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    df["codon_start"] = df["codon_start"].astype(int)
    df["residue"] = df["residue"].astype(int)
    return df


# ---------------------------------------------------------------- the default stays byte-identical


def test_config_default_is_zero() -> None:
    assert SimConfig().n_selenoproteins == 0


def test_zero_is_byte_identical_and_writes_no_truth(runs: dict[str, Path]) -> None:
    assert _files(runs["zero"]) == _files(runs["base"])
    assert not (runs["zero"] / "truth" / "selenoproteins.tsv").exists()
    assert _digest(runs["base"]) == PINNED_BASE


# ---------------------------------------------------------------- example 32: truth and genome


def test_example_32_truth_file(runs: dict[str, Path]) -> None:
    t = truth(runs["five"])
    assert list(t.columns) == COLUMNS
    assert len(t) == 5 and t["locus_tag"].is_unique
    assert list(zip(t["contig"], t["codon_start"], strict=True)) == sorted(
        zip(t["contig"], t["codon_start"], strict=True)
    )


def test_each_truth_codon_reads_tga_and_the_protein_has_u(runs: dict[str, Path]) -> None:
    out = runs["five"]
    genome = fasta(out / "genome.fasta")
    by = {f["locus"]: f for f in features(out)}
    for x in truth(out).itertuples(index=False):
        f = by[x.locus_tag]
        seq = genome[x.contig]
        s, e = interval(f, len(seq))
        assert f["strand"] == x.strand
        cs, codon = gene_codon(seq, s, e, x.strand, x.residue)
        assert (cs, codon) == (x.codon_start, "TGA"), x
        assert f["aa"][x.residue] == "U", x
        assert len(f["aa"]) == (e - s) // 3 - 1
        if f.get("nt") and len(f["nt"]) == e - s:
            assert f["nt"][3 * x.residue : 3 * x.residue + 3] == "TGA"


def test_eligibility_and_codon_position(runs: dict[str, Path]) -> None:
    out = runs["twenty"]
    by = {f["locus"]: f for f in features(out)}
    for x in truth(out).itertuples(index=False):
        f = by[x.locus_tag]
        assert f["type"] == "cds" and not f.get("pseudo") and not f.get("edge"), f["locus"]
        n_aa = len(f["aa"])
        assert n_aa >= 21
        assert 10 <= x.residue <= n_aa - 11, (x.locus_tag, x.residue, n_aa)


def test_the_codon_overlaps_no_other_feature(runs: dict[str, Path]) -> None:
    out = runs["twenty"]
    genome = fasta(out / "genome.fasta")
    feats = features(out)
    for x in truth(out).itertuples(index=False):
        n = len(genome[x.contig])
        a, b = x.codon_start, x.codon_start + 3
        for f in feats:
            if f["locus"] == x.locus_tag or f["contig"] != x.contig:
                continue
            s, e = interval(f, n)
            for d in (-n, 0, n):
                assert not (s + d < b and a < e + d), (x.locus_tag, f["locus"])


def test_only_the_planted_codons_change(runs: dict[str, Path]) -> None:
    """No other base and no other stream moves: the genome differs from the default build only
    inside the planted codons, the other proteins and the gene truth are unchanged, and no other
    in-frame stop appears in a selenoprotein."""
    five, base = runs["five"], runs["base"]
    g5, g0 = fasta(five / "genome.fasta"), fasta(base / "genome.fasta")
    t = truth(five)
    planted = {(x.contig, x.codon_start + k) for x in t.itertuples() for k in range(3)}
    for cid, seq in g5.items():
        assert len(seq) == len(g0[cid])
        diff = {(cid, i) for i, (p, q) in enumerate(zip(seq, g0[cid], strict=True)) if p != q}
        assert diff <= planted
    f5 = {f["locus"]: f for f in features(five)}
    f0 = {f["locus"]: f for f in features(base)}
    assert set(f5) == set(f0)
    sel = dict(zip(t["locus_tag"], t["residue"], strict=True))
    for tag, f in f5.items():
        if tag not in sel:
            assert f == f0[tag], tag
            continue
        i = sel[tag]
        assert f["aa"] == f0[tag]["aa"][:i] + "U" + f0[tag]["aa"][i + 1 :]
        seq = g5[f["contig"]]
        s, e = interval(f, len(seq))
        internal = [gene_codon(seq, s, e, f["strand"], k)[1] for k in range(1, (e - s) // 3 - 1)]
        stops = [k + 1 for k, c in enumerate(internal) if CODE_11.get(c, "*") == "*"]
        assert stops == [i], (tag, stops)
    for rel in ("truth/genes.tsv", "truth/features.tsv", "truth/novel_orfs.tsv"):
        assert (five / rel).read_bytes() == (base / rel).read_bytes(), rel


def test_no_pyl_is_planted(runs: dict[str, Path]) -> None:
    for f in features(runs["twenty"]):
        assert "O" not in f.get("aa", "")


# ---------------------------------------------------------------- the expressed protein


def test_true_peptides_spanning_the_codon_carry_u(runs: dict[str, Path]) -> None:
    out = runs["twenty"]
    t = truth(out)
    pep = pd.read_csv(out / "truth" / "peptides.tsv", sep="\t", dtype=str, keep_default_na=False)
    pep = pep[(pep["is_wrong"] == "false") & (pep["origin"] == "gene")]
    spanning = 0
    for x in t.itertuples(index=False):
        rows = pep[(pep["locus_tag"] == x.locus_tag) & (pep["contig"] == x.contig)]
        for r in rows.itertuples(index=False):
            s0, e = int(r.start) - 1, int(r.end)
            if not (s0 <= x.codon_start and x.codon_start + 3 <= e):
                continue
            k = (x.codon_start - s0) // 3 if x.strand == "+" else (e - x.codon_start - 3) // 3
            assert r.peptide[k] == "U", (x.locus_tag, r.peptide, k)
            spanning += 1
        # the expressed protein is not truncated at the Sec codon: a peptide that ends just before
        # it must end at a trypsin cut (K/R), never at the codon as if it were a stop
        for r in rows.itertuples(index=False):
            s0, e = int(r.start) - 1, int(r.end)
            ends_at = (x.strand == "+" and e == x.codon_start) or (
                x.strand == "-" and s0 == x.codon_start + 3
            )
            assert not ends_at or r.peptide[-1] in "KR", (x.locus_tag, r.peptide)
    assert spanning >= 1, "no true peptide spans any of the 20 planted Sec codons"


def test_the_search_database_carries_u(runs: dict[str, Path]) -> None:
    out = runs["five"]
    prot = fasta(out / "proteins.faa")
    for x in truth(out).itertuples(index=False):
        seqs = [v for k, v in prot.items() if k.split()[0] == x.locus_tag or x.locus_tag in k]
        assert seqs and all(s[x.residue] == "U" for s in seqs), x.locus_tag


# ---------------------------------------------------------------- shortfall and strain


def test_shortfall_raises(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match=r"n_selenoproteins=10000 but only \d+ eligible genes"):
        simulate(replace(BASE, n_selenoproteins=10_000), SEED, tmp_path / "x")


def test_strain_with_selenoproteins_raises(tmp_path: Path) -> None:
    cfg = SimConfig(
        genome_length=50_000, strain=StrainConfig(snv_rate_per_kb=2), n_selenoproteins=2
    )
    with pytest.raises(ValueError):
        simulate(cfg, SEED, tmp_path / "x")


def test_deterministic_for_a_seed(tmp_path: Path, runs: dict[str, Path]) -> None:
    simulate(replace(BASE, n_selenoproteins=5), SEED, tmp_path / "again")
    assert _files(tmp_path / "again") == _files(runs["five"])


# ---------------------------------------------------------------- through gippa (example 32)


def test_parse_bakta_genes_and_the_stop_sites_event(
    runs: dict[str, Path], tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    from helpers.readthrough_index import peptides_table

    from gippa.config import resolve_config
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    out = runs["five"]
    cfg = resolve_config(overrides=[f"seed={SEED}"])
    files: dict[str, Path] = {"genome": out / "genome.fasta", "bakta_json": out / "bakta_like.json"}
    spec = reg.get("parse_bakta")
    outs = {k: tmp_path / f"{k}.out" for k in spec.outputs}
    reg.run_step("parse_bakta", cfg, {k: files[k] for k in spec.inputs}, outs)
    files.update(outs)
    genes = pd.read_parquet(files["genes"]).set_index("locus_tag")
    t = truth(out)
    for x in t.itertuples(index=False):
        assert genes.loc[x.locus_tag, "aa"][x.residue] == "U"
    files["peptides"] = peptides_table(tmp_path / "peptides.parquet")
    caplog.set_level(logging.INFO)
    spec = reg.get("sixframe_index")
    reg.run_step(
        "sixframe_index",
        cfg,
        {k: files[k] for k in spec.inputs},
        {"sixframe": tmp_path / "six.npz"},
    )
    ev = [r for r in caplog.records if getattr(r, "event", None) == "mapping.stop_sites"]
    assert len(ev) == 1
    data = ev[0].data  # type: ignore[attr-defined]
    assert data["recoded"] == 5
    assert data["recoded_skipped"] == {"length_mismatch": 0, "codon_mismatch": 0}
    from gippa.mapping.sixframe import load_stop_cells

    rec = {
        (s.contig, s.strand, s.codon_start)
        for s in load_stop_cells(tmp_path / "six.npz").sites
        if s.kind == "recoded"
    }
    assert rec == {(x.contig, x.strand, x.codon_start) for x in t.itertuples(index=False)}
